"""Read an uploaded `.rds` junction x sample matrix into memory.

Two shapes are handled:

* **dense** — a base matrix / data.frame, read with ``pyreadr`` (or an ``Rscript``
  fallback that densifies to a TSV). Small; the common non-TCGA case.
* **sparse** — a ``Matrix`` package sparse matrix (``dgTMatrix`` / ``dgCMatrix``),
  e.g. a whole-transcriptome TCGA junction matrix (millions of rows). ``pyreadr``
  can't read these, so an ``Rscript`` dumps the CSC vectors as raw binary and we
  rebuild a ``scipy.sparse.csc_matrix`` — never densifying the whole thing. Plot
  requests slice to a gene locus (hundreds of rows) and densify only that.

Either way the row names (``chr:start-end:strand``) are parsed **vectorised** into
coordinate arrays; rows that don't match (e.g. unstranded ``chr1:100-200:*``) are
flagged invalid and excluded, exactly as before.

The parsed result of an expensive sparse read is cached under
``$SJV_CACHE_DIR/matrices/<sha1>`` so re-loading the same file (a new session, a
Data-tab reload) is cheap.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .junctions import Junction, norm_chrom

_RSCRIPT = shutil.which("Rscript")
_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
_RDS_TO_TSV = _SCRIPTS / "rds_to_tsv.R"
_SPARSE_RDS_TO_CSC = _SCRIPTS / "sparse_rds_to_csc.R"
# sparse_rds_to_csc.R exits with this code when the object is a plain (dense)
# matrix / data.frame — the caller then uses the dense path.
_NOT_SPARSE_RC = 3

_MAX_BAD_SAMPLE = 10

_ROWNAME_RE = r"^(?P<chrom>[^:\s]+):(?P<a>\d+)-(?P<b>\d+):(?P<strand>[+-])$"

# bump if the cache file layout below changes
_CACHE_VERSION = 1


def _cache_dir() -> Path:
    root = Path(os.environ.get("SJV_CACHE_DIR", Path.home() / ".cache" / "sjv"))
    return root / "matrices"


class RdsError(ValueError):
    """Raised for user-fixable problems with the uploaded RDS."""


# --------------------------------------------------------------------------- #
# row-name parsing (vectorised)
# --------------------------------------------------------------------------- #
@dataclass
class _Coords:
    valid: np.ndarray            # bool, len n_rows
    start: np.ndarray            # int64, -1 where invalid
    end: np.ndarray              # int64
    strand: np.ndarray           # '<U1', '' where invalid
    chrom_codes: np.ndarray      # int32 index into chrom_cats, -1 where invalid
    chrom_cats: List[str]        # original chrom labels
    chrom_norm_cats: List[str]   # norm_chrom() of each, for locus matching
    rownames: np.ndarray         # object, the raw strings


def _parse_rownames(rownames: Sequence[str]) -> _Coords:
    s = pd.Series(list(rownames), dtype="object")
    m = s.str.extract(_ROWNAME_RE)
    valid = m["chrom"].notna().to_numpy()

    a = pd.to_numeric(m["a"], errors="coerce")
    b = pd.to_numeric(m["b"], errors="coerce")
    start = np.minimum(a, b).fillna(-1).to_numpy().astype(np.int64)
    end = np.maximum(a, b).fillna(-1).to_numpy().astype(np.int64)
    strand = m["strand"].fillna("").to_numpy().astype("<U1")

    cat = pd.Categorical(m["chrom"])
    chrom_cats = [str(c) for c in cat.categories]
    return _Coords(
        valid=valid,
        start=start,
        end=end,
        strand=strand,
        chrom_codes=cat.codes.astype(np.int32),
        chrom_cats=chrom_cats,
        chrom_norm_cats=[norm_chrom(c) for c in chrom_cats],
        rownames=np.asarray(list(rownames), dtype=object),
    )


# --------------------------------------------------------------------------- #
# in-memory matrix
# --------------------------------------------------------------------------- #
@dataclass
class RdsMatrix:
    samples: List[str]
    values: Any                       # np.ndarray (n_rows, n_samples) OR scipy.sparse.csc_matrix
    _c: _Coords
    sparse: bool = False
    _col: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._col = {s: i for i, s in enumerate(self.samples)}

    # ---- counts ---------------------------------------------------------- #
    @property
    def n_rows(self) -> int:
        return len(self._c.valid)

    @property
    def n_junctions(self) -> int:
        return int(self._c.valid.sum())

    @property
    def n_bad_rownames(self) -> int:
        return int((~self._c.valid).sum())

    @property
    def bad_rownames_sample(self) -> List[str]:
        return [str(r) for r in self._c.rownames[~self._c.valid][:_MAX_BAD_SAMPLE]]

    def chrom_uses_chr_prefix(self) -> bool:
        return any(c.lower().startswith("chr") for c in self._c.chrom_cats)

    # ---- locus queries (used by the plot route) ------------------------- #
    def locus_row_indices(self, gene) -> np.ndarray:
        """Row indices (into the full matrix) of valid junctions on the gene's
        chromosome + strand whose interval overlaps the gene span."""
        gchr = norm_chrom(gene.chrom)
        try:
            code = self._c.chrom_norm_cats.index(gchr)
        except ValueError:
            return np.empty(0, dtype=np.int64)
        mask = (
            self._c.valid
            & (self._c.chrom_codes == code)
            & (self._c.strand == gene.strand)
            & (self._c.end >= gene.start)
            & (self._c.start <= gene.end)
        )
        return np.flatnonzero(mask)

    def junctions_at(self, rows: Sequence[int]) -> List[Junction]:
        c = self._c
        return [
            Junction(
                rowname=str(c.rownames[i]),
                chrom=c.chrom_cats[c.chrom_codes[i]],
                start=int(c.start[i]),
                end=int(c.end[i]),
                strand=str(c.strand[i]),
            )
            for i in rows
        ]

    def _block(self, rows: np.ndarray, cols: Sequence[int]) -> np.ndarray:
        """Dense (len(rows), len(cols)) float block. Both axes are sliced before
        anything is densified, so a multi-million-row sparse matrix costs the
        same as a small one."""
        if self.sparse:
            sub = self.values[rows, :][:, list(cols)]
            return np.asarray(sub.todense(), dtype=float)
        return self.values[np.ix_(rows, cols)].astype(float)

    def counts_at(self, rows: Sequence[int], sample: str) -> np.ndarray:
        if sample not in self._col:
            raise RdsError(f"sample {sample!r} is not in the uploaded matrix")
        rows = np.asarray(rows, dtype=np.int64)
        return self._block(rows, [self._col[sample]]).ravel()

    def group_median_at(self, rows: Sequence[int], sample_ids: Sequence[str]) -> np.ndarray:
        cols = [self._col[s] for s in sample_ids if s in self._col]
        if not cols:
            raise RdsError("none of the group's samples are columns of the matrix")
        rows = np.asarray(rows, dtype=np.int64)
        block = self._block(rows, cols)
        masked = np.where(np.isnan(block) | (block <= 0), np.nan, block)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmedian(masked, axis=1)

    # ---- legacy list/column API (dense fixtures, tests) ---------------- #
    @property
    def _valid_rows(self) -> np.ndarray:
        return np.flatnonzero(self._c.valid)

    @property
    def junctions(self) -> List[Junction]:
        return self.junctions_at(self._valid_rows)

    def column(self, sample: str) -> np.ndarray:
        return self.counts_at(self._valid_rows, sample)

    def group_median(self, sample_ids: List[str]) -> np.ndarray:
        return self.group_median_at(self._valid_rows, sample_ids)


# --------------------------------------------------------------------------- #
# dense readers (unchanged behaviour)
# --------------------------------------------------------------------------- #
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
    if isinstance(df.index, pd.RangeIndex):
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


def _matrix_from_dense_df(df: pd.DataFrame) -> RdsMatrix:
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise RdsError("the matrix is empty")
    samples = [str(c) for c in df.columns]
    rownames = [str(r) for r in df.index]
    try:
        values = df.to_numpy(dtype=float)
    except (ValueError, TypeError):
        raise RdsError("matrix values are not numeric")
    c = _parse_rownames(rownames)
    if not c.valid.any():
        raise RdsError(f"no row names matched chr:start-end:strand (e.g. got {rownames[0]!r})")
    return RdsMatrix(samples=samples, values=values, _c=c, sparse=False)


# --------------------------------------------------------------------------- #
# sparse reader
# --------------------------------------------------------------------------- #
def _read_sparse_csc(path: Path, tmp_dir: Path):
    """Returns (csc_matrix, rownames, samples) or None if the RDS is not a sparse
    matrix (the caller falls back to the dense path)."""
    if not _RSCRIPT:
        return None
    from scipy import sparse

    pre = tmp_dir / "csc"
    proc = subprocess.run(
        [_RSCRIPT, "--vanilla", str(_SPARSE_RDS_TO_CSC), str(path), str(pre)],
        capture_output=True,
        text=True,
        timeout=1800,
    )
    if proc.returncode == _NOT_SPARSE_RC:
        return None
    if proc.returncode != 0:
        raise RdsError(
            f"R could not read the sparse RDS: {proc.stderr.strip() or proc.stdout.strip()}"
        )
    try:
        indices = np.fromfile(f"{pre}.indices", dtype=np.int32)
        indptr = np.fromfile(f"{pre}.indptr", dtype=np.int32)
        data = np.fromfile(f"{pre}.data", dtype=np.float64)
        rownames = Path(f"{pre}.rows").read_text().splitlines()
        samples = Path(f"{pre}.cols").read_text().splitlines()
        nrow, ncol = (int(x) for x in Path(f"{pre}.shape").read_text().split())
    finally:
        for ext in ("indices", "indptr", "data", "rows", "cols", "shape"):
            Path(f"{pre}.{ext}").unlink(missing_ok=True)
    m = sparse.csc_matrix((data, indices, indptr), shape=(nrow, ncol))
    return m, rownames, samples


# --------------------------------------------------------------------------- #
# parsed-matrix cache (sparse only — the dense path is already fast)
# --------------------------------------------------------------------------- #
def _sha1(path: Path) -> str:
    h = hashlib.sha1()
    h.update(f"v{_CACHE_VERSION}|".encode())
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _cache_load(key: str) -> Optional[RdsMatrix]:
    from scipy import sparse

    d = _cache_dir() / key
    meta_path = d / "meta.json"
    if not meta_path.exists():
        return None
    try:
        meta = json.loads(meta_path.read_text())
        indices = np.fromfile(d / "indices.bin", dtype=np.int32)
        indptr = np.fromfile(d / "indptr.bin", dtype=np.int32)
        data = np.fromfile(d / "data.bin", dtype=np.float64)
        npz = np.load(d / "coords.npz", allow_pickle=False)
        rownames = (d / "rownames.txt").read_text().splitlines()
    except Exception:
        return None
    c = _Coords(
        valid=npz["valid"],
        start=npz["start"],
        end=npz["end"],
        strand=npz["strand"].astype("<U1"),
        chrom_codes=npz["chrom_codes"],
        chrom_cats=list(meta["chrom_cats"]),
        chrom_norm_cats=list(meta["chrom_norm_cats"]),
        rownames=np.asarray(rownames, dtype=object),
    )
    m = sparse.csc_matrix(
        (data, indices, indptr), shape=(meta["nrow"], meta["ncol"])
    )
    return RdsMatrix(samples=list(meta["samples"]), values=m, _c=c, sparse=True)


def _cache_store(key: str, m: RdsMatrix) -> None:
    try:
        d = _cache_dir() / key
        d.mkdir(parents=True, exist_ok=True)
        csc = m.values
        csc.indices.astype(np.int32, copy=False).tofile(d / "indices.bin")
        csc.indptr.astype(np.int32, copy=False).tofile(d / "indptr.bin")
        csc.data.astype(np.float64, copy=False).tofile(d / "data.bin")
        c = m._c
        np.savez(
            d / "coords.npz",
            valid=c.valid,
            start=c.start,
            end=c.end,
            strand=c.strand.astype("S1"),
            chrom_codes=c.chrom_codes,
        )
        (d / "rownames.txt").write_text("\n".join(str(r) for r in c.rownames))
        (d / "meta.json").write_text(
            json.dumps(
                {
                    "samples": m.samples,
                    "nrow": csc.shape[0],
                    "ncol": csc.shape[1],
                    "chrom_cats": c.chrom_cats,
                    "chrom_norm_cats": c.chrom_norm_cats,
                }
            )
        )
    except Exception:
        shutil.rmtree(_cache_dir() / key, ignore_errors=True)


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #
def load_rds(path: Path, tmp_dir: Path) -> RdsMatrix:
    # dense first — pyreadr handles base matrices / data.frames quickly
    df = _dataframe_from_pyreadr(path)
    if df is not None:
        return _matrix_from_dense_df(df)

    # sparse Matrix object?  (cache the parsed result — the read is ~1 min)
    key = _sha1(path)
    cached = _cache_load(key)
    if cached is not None:
        return cached

    sp = _read_sparse_csc(path, tmp_dir)
    if sp is not None:
        m, rownames, samples = sp
        c = _parse_rownames(rownames)
        if not c.valid.any():
            raise RdsError(
                f"no row names matched chr:start-end:strand (e.g. got {rownames[0]!r})"
            )
        rds = RdsMatrix(samples=[str(s) for s in samples], values=m, _c=c, sparse=True)
        _cache_store(key, rds)
        return rds

    # last resort: densify via Rscript (small base matrices, tiny Matrix objects)
    return _matrix_from_dense_df(_dataframe_from_rscript(path, tmp_dir))
