"""GENCODE sourcing and gene / transcript lookup.

A GTF (downloaded by release, or uploaded by the user) is indexed once into a
SQLite database:

    genes(name_lc, name, gene_id, gene_type, chrom, start, end, strand)
    feats(gene_id, transcript_id, feature, start, end, strand)   # exon / CDS

Lookups are then cheap. This replaces the tabix approach from the plan for v1
(no system tabix/bgzip available); the interface is the same and tabix can be
slotted in later behind `get_transcripts`.

The index filename carries `_INDEX_VERSION`; bump it whenever the schema
changes so stale indexes (including any built by SJV, which shares this cache
directory) are rebuilt rather than read with the wrong columns.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import httpx

from .classify import Exon, Transcript
from .junctions import Gene

BACKEND_DIR = Path(__file__).resolve().parents[2]
# The GENCODE cache (downloaded GTFs + SQLite indexes) is deliberately shared
# with SJV by default so a release indexed for one app is instantly available to
# the other. Override with SJVC_CACHE_DIR.
CACHE_DIR = Path(
    os.environ.get("SJVC_CACHE_DIR")
    or os.environ.get("SJV_CACHE_DIR")
    or Path.home() / ".cache" / "sjv"
)
GENCODE_DIR = CACHE_DIR / "gencode"

# Bump when the SQLite schema in build_index() changes. Baked into the index
# filename so an out-of-date index is simply not found and gets rebuilt.
_INDEX_VERSION = 2

# JSON file mapping species -> release -> local GTF path, so a pre-downloaded
# annotation is used instead of hitting the GENCODE server. Point elsewhere with
# $SJVC_GENCODE_SOURCES.
SOURCES_FILE = Path(
    os.environ.get("SJVC_GENCODE_SOURCES")
    or os.environ.get("SJV_GENCODE_SOURCES")
    or BACKEND_DIR / "gencode_sources.json"
)

_CODING_FEATURES = {"exon", "CDS"}

# Curated low-mappability / segmental-duplication / hyper-polymorphic
# gene-family name patterns (matched case-insensitively against a gene's
# lower-cased symbol) — HLA, immunoglobulin and TCR loci, the mitochondrial
# genes, olfactory receptors, and several named segmental-duplication
# clusters. These are gene *families* known to multi-map or otherwise
# confound short-read splice-junction calling; the pattern list (and its use
# as an optional exclusion filter alongside "protein-coding only") comes from
# repeated cross-cohort TCGA splice-junction analyses where these families
# dominated an unfiltered ranking with noise rather than tumour biology.
PARALOG_FAMILY_PATTERNS = [
    r"^hla-", r"^igh[vdjc]", r"^igk[vjc]", r"^igl[vjc]", r"^tr[abgd][vdjc]",
    r"^mt-", r"^or\d", r"^spdye", r"^rasa4", r"^polr2j", r"^pkd1p",
    r"^golga[68]", r"^nbpf", r"^ankrd20a", r"^usp17l", r"^pote[a-z]", r"^dux4",
]
_PARALOG_FAMILY_RE = re.compile("|".join(PARALOG_FAMILY_PATTERNS))


def is_paralog_family(name_lc: str) -> bool:
    """Whether a lower-cased gene symbol matches one of the curated
    low-mappability paralog-family patterns above."""
    return bool(_PARALOG_FAMILY_RE.match(name_lc))

_HUMAN = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{n}/gencode.v{n}.annotation.gtf.gz"
_HUMAN37 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{n}/GRCh37_mapping/gencode.v{n}lift37.annotation.gtf.gz"
_MOUSE = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M{n}/gencode.vM{n}.annotation.gtf.gz"

# Curated list shown in the UI dropdown. The backend also accepts any other
# `v<int>` (human) / `vM<int>` (mouse) release — the URL is templated on demand.
RELEASES: Dict[str, Dict[str, str]] = {
    "human": {
        "v47": _HUMAN.format(n="47"),
        "v46": _HUMAN.format(n="46"),
        "v45": _HUMAN.format(n="45"),
        "v44": _HUMAN.format(n="44"),
        "v36": _HUMAN.format(n="36"),
        "v29": _HUMAN.format(n="29"),
        "v26": _HUMAN.format(n="26"),
        "v19 (GRCh37)": _HUMAN.format(n="19"),
        "v45lift37 (GRCh37)": _HUMAN37.format(n="45"),
    },
    "mouse": {
        "vM36": _MOUSE.format(n="36"),
        "vM35": _MOUSE.format(n="35"),
        "vM34": _MOUSE.format(n="34"),
        "vM25": _MOUSE.format(n="25"),
    },
}

_HUMAN_REL_RE = re.compile(r"^v(\d+)$")
_MOUSE_REL_RE = re.compile(r"^vM(\d+)$")


class GencodeError(ValueError):
    pass


def list_releases() -> Dict[str, List[str]]:
    return {sp: list(rels) for sp, rels in RELEASES.items()}


def _load_sources() -> Dict[str, Dict[str, str]]:
    if not SOURCES_FILE.exists():
        return {}
    try:
        return json.loads(SOURCES_FILE.read_text())
    except (json.JSONDecodeError, OSError):
        return {}


def _local_source(species: str, release: str) -> Optional[Path]:
    p = _load_sources().get(species, {}).get(release)
    if not p:
        return None
    path = Path(p).expanduser()
    return path if path.exists() else None


def register_local_source(species: str, release: str, gtf_path: Path) -> None:
    data = _load_sources()
    data.setdefault(species, {})[release] = str(gtf_path)
    SOURCES_FILE.write_text(json.dumps(data, indent=2) + "\n")


def _resolve_url(species: str, release: str) -> str:
    if release in RELEASES.get(species, {}):
        return RELEASES[species][release]
    if species == "human":
        m = _HUMAN_REL_RE.match(release)
        if m:
            return _HUMAN.format(n=m.group(1))
    if species == "mouse":
        m = _MOUSE_REL_RE.match(release)
        if m:
            return _MOUSE.format(n=m.group(1))
    raise GencodeError(
        f"unknown GENCODE release {species}/{release!r} "
        "(expected e.g. v29 for human, vM25 for mouse)"
    )


# --------------------------------------------------------------------------- #
# Index building
# --------------------------------------------------------------------------- #
def _val(field: str, key: str) -> str:
    """Extract `key "value"` from a GTF attribute field without a full parse."""
    needle = key + ' "'
    i = field.find(needle)
    if i < 0:
        return ""
    i += len(needle)
    j = field.find('"', i)
    return field[i:j] if j >= 0 else ""


def build_index(gtf_path: Path, db_path: Path) -> None:
    """Single streaming pass over a (optionally gzipped) GTF -> SQLite."""
    tmp = db_path.with_suffix(".building")
    tmp.unlink(missing_ok=True)
    con = sqlite3.connect(tmp)
    con.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        CREATE TABLE genes (
            name_lc TEXT, name TEXT, gene_id TEXT, gene_type TEXT,
            chrom TEXT, start INTEGER, end INTEGER, strand TEXT
        );
        CREATE TABLE feats (
            gene_id TEXT, transcript_id TEXT, feature TEXT,
            start INTEGER, end INTEGER, strand TEXT
        );
        """
    )
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    gene_rows: List[Tuple] = []
    feat_rows: List[Tuple] = []
    with opener(gtf_path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            feature = parts[2]
            if feature == "gene":
                attr_field = parts[8]
                name = _val(attr_field, "gene_name") or _val(attr_field, "gene_id")
                gene_type = _val(attr_field, "gene_type") or _val(attr_field, "gene_biotype")
                gene_rows.append(
                    (name.lower(), name, _val(attr_field, "gene_id"), gene_type,
                     parts[0], int(parts[3]), int(parts[4]), parts[6])
                )
            elif feature in _CODING_FEATURES:
                attr_field = parts[8]
                feat_rows.append(
                    (_val(attr_field, "gene_id"), _val(attr_field, "transcript_id"),
                     feature, int(parts[3]), int(parts[4]), parts[6])
                )
                if len(feat_rows) >= 500_000:
                    con.executemany("INSERT INTO feats VALUES (?,?,?,?,?,?)", feat_rows)
                    feat_rows.clear()

    con.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?)", gene_rows)
    con.executemany("INSERT INTO feats VALUES (?,?,?,?,?,?)", feat_rows)
    con.executescript(
        """
        CREATE INDEX ix_genes_name ON genes(name_lc);
        CREATE INDEX ix_feats_gene ON feats(gene_id);
        """
    )
    con.commit()
    con.close()
    tmp.rename(db_path)


# --------------------------------------------------------------------------- #
# Annotation handle
# --------------------------------------------------------------------------- #
class Annotation:
    def __init__(self, db_path: Path, label: str):
        self.db_path = db_path
        self.label = label

    def _con(self) -> sqlite3.Connection:
        con = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True)
        con.row_factory = sqlite3.Row
        return con

    def get_gene(self, name: str) -> Optional[Gene]:
        con = self._con()
        try:
            row = con.execute(
                "SELECT * FROM genes WHERE name_lc = ? ORDER BY (end-start) DESC LIMIT 1",
                (name.strip().lower(),),
            ).fetchone()
        finally:
            con.close()
        if not row:
            return None
        return Gene(
            name=row["name"], gene_id=row["gene_id"], chrom=row["chrom"],
            start=row["start"], end=row["end"], strand=row["strand"],
        )

    def has_gene_types(self) -> bool:
        """Whether the source GTF carried gene_type / gene_biotype attributes
        at all (a minimal hand-written GTF may not)."""
        con = self._con()
        try:
            row = con.execute(
                "SELECT 1 FROM genes WHERE gene_type IS NOT NULL AND gene_type != '' LIMIT 1"
            ).fetchone()
        finally:
            con.close()
        return row is not None

    def gene_names_of_type(self, gene_type: str) -> set:
        """Lower-cased names of all genes with the given biotype, e.g.
        'protein_coding'. Lower-cased so it can be matched against a feature
        matrix's row names directly."""
        con = self._con()
        try:
            rows = con.execute(
                "SELECT DISTINCT name_lc FROM genes WHERE gene_type = ?", (gene_type,)
            ).fetchall()
        finally:
            con.close()
        return {r["name_lc"] for r in rows}

    def paralog_family_names(self) -> set:
        """Lower-cased names of every gene in the reference matching one of
        the curated low-mappability paralog-family patterns (HLA,
        immunoglobulin/TCR loci, MT-*, olfactory receptors, and the other
        named segmental-duplication clusters — see `PARALOG_FAMILY_PATTERNS`
        above), regardless of biotype. Meant to be subtracted from (or, for a
        gene-level matrix, matched directly against) a candidate gene set —
        see `sjvc.services.mad.top_features_by_mad`'s `exclude` parameter."""
        con = self._con()
        try:
            rows = con.execute("SELECT DISTINCT name_lc FROM genes").fetchall()
        finally:
            con.close()
        return {r["name_lc"] for r in rows if is_paralog_family(r["name_lc"])}

    def genes_by_prefix(self, prefix: str, limit: int = 1000) -> List[Gene]:
        """All genes whose name starts with `prefix` (case-insensitive).
        One record per gene_id (largest span wins on collision)."""
        q = prefix.strip().lower()
        if not q:
            return []
        esc = q.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
        con = self._con()
        try:
            rows = con.execute(
                "SELECT name, gene_id, chrom, start, end, strand FROM genes "
                "WHERE name_lc LIKE ? ESCAPE '\\' ORDER BY name LIMIT ?",
                (esc + "%", limit),
            ).fetchall()
        finally:
            con.close()
        out: Dict[str, Gene] = {}
        for r in rows:
            g = Gene(
                name=r["name"], gene_id=r["gene_id"], chrom=r["chrom"],
                start=r["start"], end=r["end"], strand=r["strand"],
            )
            prev = out.get(g.gene_id)
            if prev is None or (g.end - g.start) > (prev.end - prev.start):
                out[g.gene_id] = g
        return list(out.values())

    def suggest(self, prefix: str, limit: int = 20) -> List[str]:
        """Gene names starting with `prefix` (case-insensitive), for the
        query-box autocomplete."""
        q = prefix.strip().lower()
        if not q:
            return []
        con = self._con()
        try:
            rows = con.execute(
                "SELECT DISTINCT name FROM genes WHERE name_lc LIKE ? ESCAPE '\\' "
                "ORDER BY length(name), name LIMIT ?",
                (q.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "%", limit),
            ).fetchall()
        finally:
            con.close()
        return [r["name"] for r in rows]

    def near_matches(self, name: str, limit: int = 8) -> List[str]:
        q = name.strip().lower()
        con = self._con()
        try:
            rows = con.execute(
                "SELECT DISTINCT name FROM genes WHERE name_lc LIKE ? ORDER BY name LIMIT ?",
                (q + "%", limit),
            ).fetchall()
            if not rows:
                rows = con.execute(
                    "SELECT DISTINCT name FROM genes WHERE name_lc LIKE ? ORDER BY name LIMIT ?",
                    ("%" + q + "%", limit),
                ).fetchall()
        finally:
            con.close()
        return [r["name"] for r in rows]

    def get_transcripts(self, gene_id: str) -> List[Transcript]:
        con = self._con()
        try:
            rows = con.execute(
                "SELECT transcript_id, feature, start, end, strand FROM feats WHERE gene_id = ?",
                (gene_id,),
            ).fetchall()
        finally:
            con.close()

        by_tx: Dict[str, Dict[str, list]] = {}
        for r in rows:
            d = by_tx.setdefault(r["transcript_id"], {"strand": r["strand"], "exon": [], "CDS": []})
            d[r["feature"]].append((r["start"], r["end"]))

        out: List[Transcript] = []
        for txid, d in by_tx.items():
            exons_raw = sorted(d["exon"])
            cds = _merge_intervals(sorted(d["CDS"]))
            segs: List[Exon] = []
            for (es, ee) in exons_raw:
                segs.extend(_split_exon(es, ee, cds))
            if segs:
                out.append(Transcript(transcript_id=txid, strand=d["strand"], exons=tuple(segs)))
        out.sort(key=lambda t: t.transcript_id)
        return out


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[Tuple[int, int]] = []
    for s, e in intervals:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def _split_exon(es: int, ee: int, cds: List[Tuple[int, int]]) -> List[Exon]:
    """Split one exon [es, ee] into UTR / CDS / UTR segments by CDS overlap."""
    covering = [(max(es, cs), min(ee, ce)) for cs, ce in cds if cs <= ee and ce >= es]
    covering = [(s, e) for s, e in covering if s <= e]
    if not covering:
        return [Exon(es, ee, "UTR")]
    cstart = min(s for s, _ in covering)
    cend = max(e for _, e in covering)
    out: List[Exon] = []
    if es <= cstart - 1:
        out.append(Exon(es, cstart - 1, "UTR"))
    out.append(Exon(cstart, cend, "CDS"))
    if cend + 1 <= ee:
        out.append(Exon(cend + 1, ee, "UTR"))
    return out


# --------------------------------------------------------------------------- #
# Public entry points
# --------------------------------------------------------------------------- #
def _sha1(path: Path) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def annotation_from_gtf(gtf_path: Path, label: str = "custom GTF") -> Annotation:
    """Index an already-local GTF (custom upload, or a test fixture)."""
    GENCODE_DIR.mkdir(parents=True, exist_ok=True)
    db_path = GENCODE_DIR / f"custom_v{_INDEX_VERSION}_{_sha1(gtf_path)}.sqlite"
    if not db_path.exists():
        build_index(gtf_path, db_path)
    return Annotation(db_path, label)


def _rel_dir(species: str, release: str) -> Path:
    slug = f"{species}_{release}".replace(" ", "_").replace("(", "").replace(")", "")
    d = GENCODE_DIR / slug
    d.mkdir(parents=True, exist_ok=True)
    return d


def ensure_release(species: str, release: str, force: bool = False) -> Annotation:
    """Resolve a GENCODE release to an indexed handle.

    Source precedence: a local path registered in gencode_sources.json, else a
    download from the GENCODE server. The index is built once and cached.
    """
    rel_dir = _rel_dir(species, release)
    db_path = rel_dir / f"index.v{_INDEX_VERSION}.sqlite"

    if force or not db_path.exists():
        local = _local_source(species, release)
        if local is not None:
            src = local
        else:
            src = rel_dir / "annotation.gtf.gz"
            if not src.exists():
                _download(_resolve_url(species, release), src)
        build_index(src, db_path)
    return Annotation(db_path, f"GENCODE {species} {release}")


def _download(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(".part")
    try:
        with httpx.stream("GET", url, follow_redirects=True, timeout=None) as r:
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_bytes(1 << 20):
                    fh.write(chunk)
        tmp.rename(dest)
    except httpx.HTTPError as e:
        tmp.unlink(missing_ok=True)
        raise GencodeError(f"could not download {url}: {e}")
