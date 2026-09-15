"""Direct per-junction lookup against a cohort's junction metadata table
(``TCGA_<cohort>_junction_metadata.rds``) — the reverse direction of
``sjvc.services.junction_metadata``'s gene index (which answers "which
junctions overlap gene X"; this answers "what does the table say about
exactly this junction").

Reuses that module's fast R-side TSV projection (`_read_junction_metadata_df`)
with a wider column list (`scripts/junction_lookup_to_tsv.R`), and the same
alias-tolerant column resolution / "NA"-string normalisation, rather than
re-implementing either.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from sjvc.services.junction_metadata import (
    _GENE_ID_ALIASES,
    _GENE_NAME_ALIASES,
    _MISSING_TOKENS,
    _SEQ_ALIASES,
    _START_ALIASES,
    _END_ALIASES,
    _STRAND_ALIASES,
    _find,
    _read_junction_metadata_df,
    JunctionMetadataError,
)

_LOOKUP_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "junction_lookup_to_tsv.R"

# columns this feature surfaces beyond seq/start/end/strand/gene id/name —
# optional: a cohort's table missing one just leaves that field blank rather
# than failing the whole load (only the coordinate/strand columns are load-
# critical, since they're what the lookup key itself is built from)
_OPTIONAL_COLS = ("width", "annotated", "left_motif", "right_motif", "left_annotated", "right_annotated")
# kept low-cardinality on purpose -> cast to `category` to keep a ~1.6M-row
# table's memory footprint modest, matching what this app already tolerates
# for a single sparse junction matrix
_CATEGORY_COLS = ("gene_name", "left_motif", "right_motif")

_SPLIT_RE = re.compile(r"[\s,;]+")
_MAX_JUNCTIONS = 200


def parse_junction_list(text: str) -> List[str]:
    """Free text -> a deduped, order-preserving list of junction strings.
    Unlike `geneset.parse_symbols`, no case-folding — a junction's identity
    is an exact string, not a symbol to normalise."""
    seen: set = set()
    out: List[str] = []
    for tok in _SPLIT_RE.split((text or "").strip()):
        tok = tok.strip()
        if not tok or tok in seen:
            continue
        seen.add(tok)
        out.append(tok)
        if len(out) >= _MAX_JUNCTIONS:
            break
    return out


def _clean_series(s: pd.Series) -> pd.Series:
    """R's write.table() writes a missing value as the literal "NA" string
    (unlike pandas' NaN-as-empty-string convention) — normalise both to a
    real missing value, same rule `load_junction_gene_index` uses."""
    s = s.astype(str).str.strip()
    return s.mask(s.str.lower().isin(_MISSING_TOKENS))


@dataclass
class JunctionLookupIndex:
    """Built once when the junction metadata table loads; every lookup below
    is then a single indexed row access, never a scan of the cohort."""
    df: pd.DataFrame          # indexed by synthesised "chr:start-end:strand"
    n_rows: int                # rows in the source table (before any dedup)
    n_duplicate_rownames: int  # collapsed duplicate junction keys (kept first)

    def lookup(self, rowname: str) -> Optional[Dict[str, object]]:
        try:
            row = self.df.loc[rowname]
        except KeyError:
            return None
        return row.to_dict()


def load_junction_lookup_index(path: Path, tmp_dir: Path) -> JunctionLookupIndex:
    df = _read_junction_metadata_df(path, tmp_dir, script=_LOOKUP_SCRIPT)
    cols = {str(c).strip().lower(): str(c) for c in df.columns}

    seq_c, start_c, end_c = _find(cols, _SEQ_ALIASES), _find(cols, _START_ALIASES), _find(cols, _END_ALIASES)
    strand_c = _find(cols, _STRAND_ALIASES)
    missing = [
        label for label, c in (
            ("seqnames", seq_c), ("start", start_c), ("end", end_c), ("strand", strand_c),
        ) if c is None
    ]
    if missing:
        raise JunctionMetadataError(
            "the junction metadata table is missing column(s): " + ", ".join(missing)
            + ". This should be the TCGA_<cohort>_junction_metadata.rds file written by "
              "prepTCGAdata::annotate_sj()."
        )
    gid_c, gname_c = _find(cols, _GENE_ID_ALIASES), _find(cols, _GENE_NAME_ALIASES)
    opt_c = {name: cols.get(name) for name in _OPTIONAL_COLS}

    n_rows = len(df)
    rowname = (
        df[seq_c].astype(str) + ":" + df[start_c].astype(str)
        + "-" + df[end_c].astype(str) + ":" + df[strand_c].astype(str)
    )

    # every column below is assigned as a bare numpy array, not a Series —
    # `df`'s original integer index and the junction-string index `out` is
    # about to get are unrelated, and a Series-to-Series column assignment
    # aligns by index label (turning every value to NaN here, since almost
    # nothing would match); a plain array assigns positionally instead,
    # which is what's actually wanted (row i's gene_id for row i's junction)
    out = pd.DataFrame(index=pd.Index(rowname.to_numpy(), name="junction"))
    out["gene_id"] = _clean_series(df[gid_c]).to_numpy() if gid_c else None
    out["gene_name"] = _clean_series(df[gname_c]).to_numpy() if gname_c else None
    out["width"] = pd.to_numeric(df[opt_c["width"]], errors="coerce").to_numpy() if opt_c["width"] else None
    if opt_c["annotated"]:
        annotated = _clean_series(df[opt_c["annotated"]])
        out["annotated"] = annotated.map({"1": True, "0": False}).to_numpy()
    else:
        out["annotated"] = None
    for name in ("left_motif", "right_motif", "left_annotated", "right_annotated"):
        c = opt_c[name]
        out[name] = _clean_series(df[c]).to_numpy() if c else None

    for col in _CATEGORY_COLS:
        out[col] = out[col].astype("category")

    n_duplicate_rownames = int(out.index.duplicated().sum())
    if n_duplicate_rownames:
        out = out[~out.index.duplicated(keep="first")]

    return JunctionLookupIndex(df=out, n_rows=n_rows, n_duplicate_rownames=n_duplicate_rownames)
