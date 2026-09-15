"""Clinical table: rows = samples, must have a `sample_id` column matching the
junction-matrix column names. Every other column is a candidate clinical
feature.

Column typing
-------------
- **numeric**: every non-missing value parses as a number. Always eligible.
- **categorical**: otherwise. Eligible only when `n_unique < n_rows` (drops IDs /
  free-text / all-distinct columns).

Auto-typing can be overridden per feature via `apply_overrides`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .rds import MatrixError, read_table

_SAMPLE_ID_ALIASES = {"sample_id", "sampleid", "sample id", "sample.id"}
_MISSING = {"", "nan", "none", "na", "n/a", "null"}
# a numeric column with few distinct values reads better as categorical by default
_LOW_CARDINALITY = 10
# a categorical column qualifies only if its distinct-value count is below this
# fraction of the row count (drops IDs / near-unique free-text)
_CATEGORICAL_MAX_FRACTION = 0.9


def _is_missing(v: object) -> bool:
    return v is None or str(v).strip().lower() in _MISSING


def _looks_numeric(series: pd.Series) -> bool:
    vals = [str(v).strip() for v in series if not _is_missing(v)]
    if not vals:
        return False
    try:
        [float(v) for v in vals]
        return True
    except ValueError:
        return False


@dataclass
class Column:
    name: str
    type: str                 # "numeric" | "categorical"
    n_unique: int
    n_missing: int
    eligible: bool


@dataclass
class Clinical:
    # sample_id -> {column: raw string value}
    rows: Dict[str, Dict[str, str]]
    columns: List[Column]
    unmatched_samples: List[str] = field(default_factory=list)
    overrides: Dict[str, str] = field(default_factory=dict)

    @property
    def n_rows(self) -> int:
        return len(self.rows)

    def _col(self, name: str) -> Column:
        for c in self.columns:
            if c.name == name:
                return c
        raise MatrixError(f"unknown clinical column {name!r}")

    def col_type(self, name: str) -> str:
        return self.overrides.get(name) or self._col(name).type

    def eligible_columns(self) -> List[Column]:
        return [c for c in self.columns if c.eligible]

    def apply_overrides(self, overrides: Dict[str, str]) -> None:
        for name, t in (overrides or {}).items():
            if t not in ("numeric", "categorical"):
                continue
            col = self._col(name)
            if t == "numeric":
                raw = [self.rows[s].get(name, "") for s in self.rows]
                if not _looks_numeric(pd.Series(raw)):
                    raise MatrixError(f"column {name!r} cannot be read as numeric")
            self.overrides[name] = t

    def values_for(self, name: str, sample_order: Sequence[str]):
        """Return (values, kind) aligned to `sample_order`.

        numeric   -> float ndarray with np.nan for missing
        categorical -> list[str | None]
        """
        kind = self.col_type(name)
        raw = [self.rows.get(s, {}).get(name, "") for s in sample_order]
        if kind == "numeric":
            out = np.array(
                [np.nan if _is_missing(v) else float(str(v).strip()) for v in raw],
                dtype=float,
            )
            return out, "numeric"
        return [None if _is_missing(v) else str(v).strip() for v in raw], "categorical"


def parse_clinical(
    path: Path,
    tmp_dir: Path,
    known_samples: Sequence[str],
) -> Clinical:
    df = read_table(path, tmp_dir)
    if df.shape[1] < 2:
        raise MatrixError("clinical table needs a `sample_id` column plus at least one more")

    id_col = next((c for c in df.columns if str(c).strip().lower() in _SAMPLE_ID_ALIASES), None)
    if id_col is None:
        raise MatrixError(
            "no `sample_id` column found; table columns: " + ", ".join(map(str, df.columns))
        )

    known = {str(s) for s in known_samples}
    other_cols = [str(c) for c in df.columns if c != id_col]

    rows: Dict[str, Dict[str, str]] = {}
    unmatched: List[str] = []
    for _, r in df.iterrows():
        sid = str(r[id_col]).strip()
        if _is_missing(sid):
            continue
        if sid in known:
            rows[sid] = {c: str(r[c]).strip() for c in other_cols}
        else:
            unmatched.append(sid)

    if not rows:
        raise MatrixError("no `sample_id` value matched a column of the junction matrix")

    n_rows = len(rows)
    columns: List[Column] = []
    for c in other_cols:
        series = pd.Series([rows[s].get(c, "") for s in rows])
        present = [v for v in series if not _is_missing(v)]
        n_missing = n_rows - len(present)
        n_unique = len(set(present))
        numeric = _looks_numeric(series)
        if numeric and n_unique > _LOW_CARDINALITY:
            ctype = "numeric"
        elif numeric:
            # numeric but low-cardinality: default categorical, user can override
            ctype = "categorical"
        else:
            ctype = "categorical"
        eligible = (
            True
            if numeric
            else (1 <= n_unique < _CATEGORICAL_MAX_FRACTION * n_rows)
        )
        columns.append(
            Column(name=c, type=ctype, n_unique=n_unique, n_missing=n_missing, eligible=eligible)
        )

    return Clinical(
        rows=rows, columns=columns, unmatched_samples=sorted(set(unmatched))[:20]
    )
