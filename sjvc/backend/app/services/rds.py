"""Read uploaded matrices and tables.

- ``load_matrix`` reads the junction count matrix (``.rds``; rows = junctions,
  columns = samples) as a plain numeric matrix. Row names are kept verbatim; the
  ``chr:start-end:strand`` form is only parsed later, for the selected subset.
- ``read_table`` reads the clinical table (``.csv`` / ``.tsv`` / ``.rds``) with a
  header row, as a string DataFrame.

RDS reading tries ``pyreadr`` first and falls back to an ``Rscript`` subprocess.
The raw upload is the caller's responsibility to delete once these return.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

_RSCRIPT = shutil.which("Rscript")
_RDS_TO_TSV = Path(__file__).resolve().parents[2] / "scripts" / "rds_to_tsv.R"


class MatrixError(ValueError):
    """User-fixable problem with an uploaded matrix or table."""


# --------------------------------------------------------------------------- #
# RDS -> DataFrame
# --------------------------------------------------------------------------- #
def _is_numeric_series(s: pd.Series) -> bool:
    nn = s.dropna()
    if nn.empty:
        return True
    return bool(pd.to_numeric(nn, errors="coerce").notna().all())


def _promote_row_labels(df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """For a **matrix** with no row names but leading text columns — e.g. a
    gene-count tibble laid out as ``gencode_gene_id, gencode_gene_name,
    <sample>, …`` — take the row labels from the last leading text column
    (preferring one named like *name*), drop the other leading text columns,
    and drop rows whose label repeats. Returns None when there is no usable
    text column (a plain headerless numeric matrix — let the Rscript path
    handle it).

    Only ``load_matrix`` uses this; a clinical table is read verbatim so its
    columns and every one of its rows are kept."""
    lead: List = []
    for c in df.columns:
        if _is_numeric_series(df[c]):
            break
        lead.append(c)
    if not lead or len(lead) == len(df.columns):
        return None
    label = next((c for c in reversed(lead) if "name" in str(c).lower()), lead[-1])
    df = df.drop(columns=[c for c in lead if c != label]).set_index(label)
    if df.index.has_duplicates:
        df = df[~df.index.duplicated(keep="first")]
    return df


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
    df = next(iter(result.values()))
    if not isinstance(df, pd.DataFrame) or df.shape[1] == 0:
        return None
    return df


def _dataframe_from_rscript(path: Path, tmp_dir: Path) -> pd.DataFrame:
    if not _RSCRIPT:
        raise MatrixError(
            "could not read the RDS with pyreadr and Rscript is not installed; "
            "re-save it as a base matrix / data.frame or install R"
        )
    out = tmp_dir / "rds_dense.tsv"
    proc = subprocess.run(
        [_RSCRIPT, "--vanilla", str(_RDS_TO_TSV), str(path), str(out)],
        capture_output=True,
        text=True,
        timeout=600,
    )
    if proc.returncode != 0:
        raise MatrixError(f"R could not read the RDS: {proc.stderr.strip() or proc.stdout.strip()}")
    df = pd.read_csv(out, sep="\t", index_col=0)
    out.unlink(missing_ok=True)
    return df


def _read_rds_dataframe(path: Path, tmp_dir: Path, *, promote_labels: bool = False) -> pd.DataFrame:
    df = _dataframe_from_pyreadr(path)
    if df is not None and promote_labels and isinstance(df.index, pd.RangeIndex):
        df = _promote_row_labels(df)  # None -> fall through to the Rscript path
    if df is None:
        df = _dataframe_from_rscript(path, tmp_dir)
    return df


# --------------------------------------------------------------------------- #
# Junction count matrix
# --------------------------------------------------------------------------- #
@dataclass
class Matrix:
    samples: List[str]
    features: List[str]              # row names, verbatim
    values: np.ndarray              # (n_features, n_samples), float (may hold NaN)
    _col: Dict[str, int] = field(default_factory=dict)
    _row: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._col = {s: i for i, s in enumerate(self.samples)}
        self._row = {f: i for i, f in enumerate(self.features)}

    @property
    def n_features(self) -> int:
        return len(self.features)

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    def submatrix(self, feature_names: List[str]) -> np.ndarray:
        idx = [self._row[f] for f in feature_names if f in self._row]
        return self.values[idx, :]


def load_matrix(path: Path, tmp_dir: Path) -> Matrix:
    df = _read_rds_dataframe(path, tmp_dir, promote_labels=True)
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise MatrixError("the matrix is empty")
    try:
        values = df.to_numpy(dtype=float)
    except (ValueError, TypeError):
        raise MatrixError("matrix values are not all numeric")
    return Matrix(
        samples=[str(c) for c in df.columns],
        features=[str(r) for r in df.index],
        values=values,
    )


# --------------------------------------------------------------------------- #
# Generic table (clinical)
# --------------------------------------------------------------------------- #
def read_table(path: Path, tmp_dir: Path) -> pd.DataFrame:
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".rds"):
        df = _read_rds_dataframe(path, tmp_dir)
        if not isinstance(df.index, pd.RangeIndex):
            df = df.reset_index(drop=df.index.name in (None, "index"))
        return df.astype(str)
    sep = "," if suffix.endswith(".csv") else "\t"
    kind = "CSV" if sep == "," else "TSV"
    try:
        return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    except Exception as e:
        raise MatrixError(
            f"could not read {path.name!r} as {kind} ({e.__class__.__name__}: {e}); "
            f"check it's really {'comma' if sep == ',' else 'tab'}-separated, or rename it "
            f"with a .{'tsv' if sep == ',' else 'csv'} extension if it's the other delimiter"
        )
