"""Optional sample-metadata table for group visualisation.

Accepts any `.csv` / `.tsv` / `.rds` (matrix or data.frame) that has a
**`sample_id` column** whose entries line up with the junction-matrix column
names. Every other column is a candidate stratification variable: the user picks
one, then picks a value of it, and the plot shows — per splice junction — the
**median** read count across that group's samples, **excluding samples with NA or
a 0 count** for the junction.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Set

import pandas as pd

from .rds import RdsError, _dataframe_from_pyreadr, _dataframe_from_rscript

_SAMPLE_ID_ALIASES = {"sample_id", "sampleid", "sample id", "sample.id"}
_MISSING = {"", "nan", "none", "na", "n/a"}


@dataclass
class SampleMetadata:
    # sample_id -> {column: value}, only for samples present in the matrix
    rows: Dict[str, Dict[str, str]]
    columns: List[str]                                # stratification candidates
    unmatched_samples: List[str] = field(default_factory=list)

    @property
    def n_matched(self) -> int:
        return len(self.rows)

    def column_summaries(self) -> List[Dict[str, object]]:
        return [{"name": c, "n_values": len(self.strat_values(c))} for c in self.columns]

    def strat_values(self, column: str) -> Dict[str, List[str]]:
        """value -> [sample_id, ...] for one stratification column."""
        if column not in self.columns:
            return {}
        out: Dict[str, List[str]] = {}
        for sid, row in self.rows.items():
            v = row.get(column, "")
            if v is None or str(v).strip().lower() in _MISSING:
                continue
            out.setdefault(str(v), []).append(sid)
        return out


def _norm(name: str) -> str:
    return str(name).strip().lower()


def _read_named_table(path: Path, tmp_dir: Path) -> pd.DataFrame:
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".rds"):
        df = _dataframe_from_pyreadr(path)
        if df is None:
            df = _dataframe_from_rscript(path, tmp_dir)
        # bring any index (row names / the rscript "rowname" col) back as data
        if not isinstance(df.index, pd.RangeIndex):
            df = df.reset_index(drop=df.index.name in (None, "index"))
        return df.astype(str)
    sep = "," if suffix.endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def parse_sample_metadata(
    path: Path,
    tmp_dir: Path,
    known_samples: Sequence[str],
) -> SampleMetadata:
    df = _read_named_table(path, tmp_dir)
    if df.shape[1] < 2:
        raise RdsError("metadata table needs a `sample_id` column plus at least one more column")

    # locate the sample_id column
    id_col = next((c for c in df.columns if _norm(c) in _SAMPLE_ID_ALIASES), None)
    if id_col is None:
        raise RdsError(
            "no `sample_id` column found; the table has: " + ", ".join(map(str, df.columns))
        )

    known: Set[str] = {str(s) for s in known_samples}
    other_cols = [str(c) for c in df.columns if c != id_col]

    rows: Dict[str, Dict[str, str]] = {}
    unmatched: List[str] = []
    for _, r in df.iterrows():
        sid = str(r[id_col]).strip()
        if not sid or sid.lower() in _MISSING:
            continue
        if sid in known:
            rows[sid] = {c: str(r[c]).strip() for c in other_cols}
        else:
            unmatched.append(sid)

    if not rows:
        raise RdsError(
            "no `sample_id` value in the table matched a column of the uploaded matrix"
        )

    return SampleMetadata(
        rows=rows, columns=other_cols, unmatched_samples=sorted(set(unmatched))[:20]
    )


# Backwards-compatible alias (older imports/tests)
def parse_group_table(path: Path, tmp_dir: Path, known_samples: Sequence[str]) -> SampleMetadata:
    return parse_sample_metadata(path, tmp_dir, known_samples)
