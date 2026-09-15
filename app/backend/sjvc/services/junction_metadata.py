"""Fast gene lookup for a junction-level matrix (junction counts, RRS
scores), built from the per-junction annotation table prepTCGAdata writes
alongside them — ``TCGA_<cohort>_junction_metadata.rds``.

That table already carries, for every junction, which gene(s) it overlaps
(``gencode_gene_id`` / ``gencode_gene_name``, comma-joined when a junction
sits inside more than one gene) — precomputed once, offline, when the cohort
was built. This is an alternative to ``junctions.map_junctions_to_genes()``,
which computes the same overlap *live* from a loaded GENCODE reference: for
a typed-genes / pathway / uploaded-list gene set, looking it up in this table
is a plain lookup, needs no GENCODE release loaded at all, and is the
preferred path whenever this file is available (see ``sjvc/api/routes.py``'s
``set_geneset``).
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from .junctions import JunctionGeneMap
from .rds import MatrixError, read_table

_RSCRIPT = shutil.which("Rscript")
_JM_TO_TSV = Path(__file__).resolve().parents[2] / "scripts" / "junction_metadata_to_tsv.R"

_SEQ_ALIASES = ("seqnames", "chr", "chrom", "chromosome")
_START_ALIASES = ("start",)
_END_ALIASES = ("end",)
_STRAND_ALIASES = ("strand",)
_GENE_ID_ALIASES = ("gencode_gene_id", "gene_id")
_GENE_NAME_ALIASES = ("gencode_gene_name", "gene_name", "gene_symbol")
_MISSING_TOKENS = {"", "nan", "none", "na", "n/a", "null"}


class JunctionMetadataError(ValueError):
    pass


def _find(cols: Dict[str, str], aliases: Sequence[str]) -> Optional[str]:
    return next((cols[a] for a in aliases if a in cols), None)


@dataclass
class JunctionGeneIndex:
    """Built once when the junction metadata table loads; every lookup below
    is then O(genes asked for), never O(rows in the cohort)."""
    gene_id_to_name: Dict[str, str]
    gene_ids_by_name: Dict[str, List[str]]      # UPPER(gene name) -> [gene_id, ...]
    rownames_by_gene_id: Dict[str, List[str]]   # gene_id -> every junction rowname overlapping it
    n_rows: int

    def resolve(
        self, symbols: Sequence[str], *, prefix: bool = False
    ) -> Tuple[JunctionGeneMap, Dict[str, str], List[str]]:
        """symbols -> (JunctionGeneMap of just their junctions, {gene_id:
        name} for those genes, symbols that matched no gene)."""
        gene_ids: List[str] = []
        seen: set = set()
        unmatched: List[str] = []
        for sym in symbols:
            q = sym.strip().upper()
            if not q:
                continue
            if prefix:
                hits = [
                    gid for name, gids in self.gene_ids_by_name.items()
                    if name.startswith(q) for gid in gids
                ]
            else:
                hits = self.gene_ids_by_name.get(q, [])
            if not hits:
                unmatched.append(sym)
                continue
            for gid in hits:
                if gid not in seen:
                    seen.add(gid)
                    gene_ids.append(gid)

        by_rowname: Dict[str, List[str]] = {}
        for gid in gene_ids:
            for rn in self.rownames_by_gene_id.get(gid, []):
                by_rowname.setdefault(rn, []).append(gid)

        labels = {gid: self.gene_id_to_name.get(gid, gid) for gid in gene_ids}
        jgmap = JunctionGeneMap(by_rowname=by_rowname, n_bad_rownames=0, bad_sample=[])
        return jgmap, labels, unmatched


def _read_junction_metadata_df(path: Path, tmp_dir: Path, *, script: Path = _JM_TO_TSV) -> pd.DataFrame:
    """Prefer a dedicated R-side projection down to just the columns this
    index needs: this table's other columns (dozens, some list-like) make
    `pyreadr` dramatically slower here than R's own readRDS + a plain TSV
    write — for an 8-9M-row cohort, 100s+ vs. well under half that. Falls
    back to the generic (slower, but column-alias-tolerant) reader for a
    non-.rds file, a missing Rscript, or a genuinely differently-shaped
    table the dedicated script's fixed canonical column names don't find.

    `script` defaults to this module's own narrow (gene-index-only) column
    projection; pass a different one (e.g. sjlookup's wider projection) to
    read more columns through the same fast R path without re-implementing
    the subprocess/fallback plumbing."""
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".rds") and _RSCRIPT:
        out = tmp_dir / "junction_metadata.tsv"
        proc = subprocess.run(
            [_RSCRIPT, "--vanilla", str(script), str(path), str(out)],
            capture_output=True, text=True, timeout=600,
        )
        if proc.returncode == 0:
            try:
                return pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
            finally:
                out.unlink(missing_ok=True)
    try:
        return read_table(path, tmp_dir)
    except MatrixError as e:
        raise JunctionMetadataError(str(e))


def load_junction_gene_index(path: Path, tmp_dir: Path) -> JunctionGeneIndex:
    df = _read_junction_metadata_df(path, tmp_dir)
    cols = {str(c).strip().lower(): str(c) for c in df.columns}

    seq_c, start_c, end_c = _find(cols, _SEQ_ALIASES), _find(cols, _START_ALIASES), _find(cols, _END_ALIASES)
    strand_c = _find(cols, _STRAND_ALIASES)
    gid_c, gname_c = _find(cols, _GENE_ID_ALIASES), _find(cols, _GENE_NAME_ALIASES)
    missing = [
        label for label, c in (
            ("seqnames", seq_c), ("start", start_c), ("end", end_c), ("strand", strand_c),
            ("gencode_gene_id", gid_c), ("gencode_gene_name", gname_c),
        ) if c is None
    ]
    if missing:
        raise JunctionMetadataError(
            "the junction metadata table is missing column(s): " + ", ".join(missing)
            + ". This should be the TCGA_<cohort>_junction_metadata.rds file written by "
              "prepTCGAdata::annotate_sj()."
        )

    n_rows = len(df)

    # A single plain-Python pass over numpy arrays, not pandas .str methods —
    # pandas' vectorised-looking string accessors (.split/.strip/.upper) are
    # themselves implemented as a Python loop per call, so chaining several
    # of them (plus .explode()) over many millions of rows is *slower* than
    # one hand-written pass doing the same work in a single sweep. This
    # matters here: a cohort's junction metadata table can be 8-9M rows, and
    # the whole point of this index is to make gene lookups fast.
    seqs = df[seq_c].to_numpy()
    starts = df[start_c].to_numpy()
    ends = df[end_c].to_numpy()
    strands = df[strand_c].to_numpy()
    raw_gids = df[gid_c].to_numpy()
    raw_gnames = df[gname_c].to_numpy()

    gene_id_to_name: Dict[str, str] = {}
    rownames_by_gene_id: Dict[str, List[str]] = {}
    gene_ids_by_name: Dict[str, List[str]] = {}
    for seq, start, end, strand, gid_raw, gname_raw in zip(
        seqs, starts, ends, strands, raw_gids, raw_gnames
    ):
        # R's write.table() writes a missing value as the literal "NA"
        # (unlike pandas' NaN-as-empty-string convention) — normalise case
        # so either source of this table is recognised the same way.
        if gid_raw.strip().lower() in _MISSING_TOKENS or gname_raw.strip().lower() in _MISSING_TOKENS:
            continue
        ids = gid_raw.split(",")
        names = gname_raw.split(",")
        if len(ids) != len(names):
            continue  # a row's gene_id / gene_name lists should pair up 1:1 — a data issue upstream
        rowname = f"{seq}:{start}-{end}:{strand}"
        for gid, gname in zip(ids, names):
            gid, gname = gid.strip(), gname.strip()
            if not gid or not gname:
                continue
            gene_id_to_name.setdefault(gid, gname)
            rlist = rownames_by_gene_id.get(gid)
            if rlist is None:
                rownames_by_gene_id[gid] = [rowname]
            else:
                rlist.append(rowname)
            gname_u = gname.upper()
            bucket = gene_ids_by_name.get(gname_u)
            if bucket is None:
                gene_ids_by_name[gname_u] = [gid]
            elif gid not in bucket:
                bucket.append(gid)

    if not gene_id_to_name:
        raise JunctionMetadataError("no usable gene annotation rows in the junction metadata table")
    return JunctionGeneIndex(
        gene_id_to_name=gene_id_to_name, gene_ids_by_name=gene_ids_by_name,
        rownames_by_gene_id=rownames_by_gene_id, n_rows=n_rows,
    )
