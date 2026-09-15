"""Read an uploaded `.rds` junction x sample matrix into memory.

Strategy: try ``pyreadr`` first; if that cannot produce a 2-D structure with row
and column names (pyreadr does not support base matrices in every version), fall
back to an ``Rscript`` subprocess that densifies the matrix to a TSV. The raw
upload is the caller's responsibility to delete once this returns.
"""
from __future__ import annotations

import io
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .junctions import Junction, RownameError, parse_rowname

_RSCRIPT = shutil.which("Rscript")
_RDS_TO_TSV = Path(__file__).resolve().parents[2] / "scripts" / "rds_to_tsv.R"

_MAX_BAD_SAMPLE = 10


class RdsError(ValueError):
    """Raised for user-fixable problems with the uploaded RDS."""


@dataclass
class RdsMatrix:
    samples: List[str]
    junctions: List[Junction]          # aligned with rows of `values`
    values: np.ndarray                 # shape (n_junctions, n_samples), float
    n_bad_rownames: int = 0
    bad_rownames_sample: List[str] = field(default_factory=list)
    _col: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._col = {s: i for i, s in enumerate(self.samples)}

    @property
    def n_junctions(self) -> int:
        return len(self.junctions)

    def column(self, sample: str) -> np.ndarray:
        if sample not in self._col:
            raise RdsError(f"sample {sample!r} is not in the uploaded matrix")
        return self.values[:, self._col[sample]]

    def group_median(self, sample_ids: List[str]) -> np.ndarray:
        """Per-junction median across the given samples, excluding entries that
        are NA or <= 0. Rows with no qualifying sample come back as NaN (and are
        dropped downstream like any other non-positive count)."""
        import warnings

        idx = [self._col[s] for s in sample_ids if s in self._col]
        if not idx:
            raise RdsError("none of the group's samples are columns of the matrix")
        sub = self.values[:, idx].astype(float)
        masked = np.where(np.isnan(sub) | (sub <= 0), np.nan, sub)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmedian(masked, axis=1)


def _dataframe_from_pyreadr(path: Path) -> Optional[pd.DataFrame]:
    try:
        import pyreadr
    except Exception:
        return None
    try:
        result = pyreadr.read_r(str(path))
    except Exception:
        return None
    if not result:
        return None
    # RDS -> single object under key None
    df = next(iter(result.values()))
    if not isinstance(df, pd.DataFrame) or df.shape[1] == 0:
        return None
    # pyreadr puts matrix/data.frame row names in the index only when they exist
    if isinstance(df.index, pd.RangeIndex):
        # first column might be the row names (common when R rownames were a col)
        first = df.columns[0]
        if df[first].dtype == object:
            df = df.set_index(first)
        else:
            return None
    return df


def _dataframe_from_rscript(path: Path, tmp_dir: Path) -> pd.DataFrame:
    if not _RSCRIPT:
        raise RdsError(
            "could not read the RDS with pyreadr and Rscript is not installed; "
            "re-save the object as a base matrix or install R"
        )
    out = tmp_dir / "rds_dense.tsv"
    proc = subprocess.run(
        [_RSCRIPT, "--vanilla", str(_RDS_TO_TSV), str(path), str(out)],
        capture_output=True,
        text=True,
        timeout=600,
    )
    if proc.returncode != 0:
        raise RdsError(f"R could not read the RDS: {proc.stderr.strip() or proc.stdout.strip()}")
    df = pd.read_csv(out, sep="\t", index_col=0)
    out.unlink(missing_ok=True)
    return df


def load_rds(path: Path, tmp_dir: Path) -> RdsMatrix:
    df = _dataframe_from_pyreadr(path)
    if df is None:
        df = _dataframe_from_rscript(path, tmp_dir)

    if df.shape[0] == 0 or df.shape[1] == 0:
        raise RdsError("the matrix is empty")

    samples = [str(c) for c in df.columns]
    rownames = [str(r) for r in df.index]

    try:
        values_all = df.to_numpy(dtype=float)
    except (ValueError, TypeError):
        raise RdsError("matrix values are not numeric")

    junctions: List[Junction] = []
    keep_rows: List[int] = []
    bad: List[str] = []
    for i, rn in enumerate(rownames):
        try:
            junctions.append(parse_rowname(rn))
            keep_rows.append(i)
        except RownameError:
            bad.append(rn)

    if not junctions:
        raise RdsError(
            "no row names matched chr:start-end:strand "
            f"(e.g. got {rownames[0]!r})"
        )

    values = values_all[keep_rows, :]
    return RdsMatrix(
        samples=samples,
        junctions=junctions,
        values=values,
        n_bad_rownames=len(bad),
        bad_rownames_sample=bad[:_MAX_BAD_SAMPLE],
    )
