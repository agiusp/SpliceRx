"""Load an ``sjdat`` feature matrix (features x samples) — shared by NSJCG
(the 2D-view tab) and SJSurv, the two apps that let a user pick which of a
cohort's several matrices to analyse.

Four flavours, all ``prepTCGAdata`` outputs:

* ``junction_counts``  — ``TCGA_<c>_junction_counts.rds``: the raw junction x
  sample count matrix. Millions of rows, usually a sparse ``Matrix`` object;
  read (and cached) via the SJV sparse reader and kept sparse.
* ``rrs_scores``       — ``TCGA_<c>_novel_junction_RRS_scores.rds``: a
  relative-read-support score matrix (feature x sample). Also junction-level
  and often just as large — kept sparse too.
* ``gene_matrix``      — ``TCGA_<c>_novel_junction_counts_per_gene.rds``: one
  row per gene. Dense, ~15-25 MB.
* ``pathway_matrix``   — ``TCGA_<c>_novel_junction_counts_per_pathway.rds``:
  one row per pathway / gene-set signature (e.g. an xCell cell-type score).
  Dense, and small — a few hundred rows.

Which reader a file needs is decided by its actual row names, not by the
requested ``kind`` label (a cohort's RRS file is exactly as junction-level as
its junction-count file, and both can be too big to safely load as a dense
matrix): the sparse-capable SJV reader is tried first, and only a matrix
whose rows aren't ``chr:start-end:strand`` at all (a real gene- or pathway-
level matrix) falls through to the SJVC dense reader, which promotes a
leading ``gencode_gene_id / gencode_gene_name`` pair to row labels when
present (a gene matrix) or otherwise reads whatever row names the RDS object
already carries (a pathway matrix).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set

import numpy as np

from sjv.services.rds import RdsError as _SjvRdsError, load_rds as _load_junctions
from sjvc.services.rds import MatrixError as _SjvcMatrixError, load_matrix as _load_dense

SJDAT_KINDS = ("junction_counts", "rrs_scores", "gene_matrix", "pathway_matrix")

# user-facing copy for the sjdat picker — these are prepTCGAdata outputs
SJDAT_META = {
    "junction_counts": (
        "Junction counts",
        "Raw split-read counts for every splice junction detected in the cohort "
        "(prepTCGAdata: TCGA_<cohort>_junction_counts.rds). Millions of junctions "
        "× samples — the finest-grained signal, and the slowest to load.",
    ),
    "rrs_scores": (
        "RRS scores",
        "Relative read-support scores — each junction's reads normalised against "
        "the local splicing context, so values are comparable across genes and "
        "samples (prepTCGAdata: TCGA_<cohort>_novel_junction_RRS_scores.rds).",
    ),
    "gene_matrix": (
        "Novel junction counts per gene",
        "Per gene, the number of novel (unannotated) splice junctions with support "
        "in each sample (prepTCGAdata: TCGA_<cohort>_novel_junction_counts_per_gene.rds). "
        "One row per gene — compact and quick.",
    ),
    "pathway_matrix": (
        "Novel junction counts per pathway",
        "Per pathway / gene-set signature, a novel-splicing score aggregated across its "
        "member genes (prepTCGAdata: TCGA_<cohort>_novel_junction_counts_per_pathway.rds). "
        "One row per pathway — a few hundred, the smallest and quickest of the four.",
    ),
}


class SjdatError(ValueError):
    """User-fixable problem with an sjdat matrix."""


@dataclass
class Sjdat:
    kind: str                       # "junction_counts" | "rrs_scores" | "gene_matrix"
    features: List[str]             # row labels, verbatim
    samples: List[str]
    values: object                  # np.ndarray (n_features, n_samples) OR scipy.sparse.csc_matrix
    sparse: bool = False
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

    def column_indices(self, sample_ids: List[str]) -> List[int]:
        return [self._col[s] for s in sample_ids if s in self._col]

    def dense_block(self, row_idx, col_idx) -> np.ndarray:
        """Dense ``(len(row_idx), len(col_idx))`` float block. Both axes are
        sliced before any densification, so a multi-million-row sparse matrix
        costs no more than a small one."""
        row_idx = list(row_idx)
        col_idx = list(col_idx)
        if self.sparse:
            sub = self.values[row_idx, :][:, col_idx]
            return np.asarray(sub.todense(), dtype=float)
        return np.asarray(self.values[np.ix_(row_idx, col_idx)], dtype=float)

    def submatrix(self, feature_names: Sequence[str]) -> np.ndarray:
        """Dense ``(n_found, n_samples)`` block for the named rows, in this
        matrix's own row order — a `Matrix`-compatible convenience so code
        written against the plain dense `sjvc.services.rds.Matrix` (gene-set
        selection, MAD's top-N slice, …) works unchanged against a Sjdat,
        sparse or not. Names not present are silently skipped, same as
        `Matrix.submatrix`."""
        idx = [self._row[f] for f in feature_names if f in self._row]
        return self.dense_block(idx, list(range(self.n_samples)))

    def nonzero_and_ge_counts(self, col_idx: List[int], x_min: float) -> np.ndarray:
        """Per row, over the given columns, the number of entries that are
        non-zero **and** ``>= x_min`` — the 'qualifying' entries feature
        selection counts. Computed without densifying the whole matrix."""
        col_idx = list(col_idx)
        if self.sparse:
            sub = self.values[:, col_idx].tocsr()
            out = np.zeros(sub.shape[0], dtype=np.int64)
            # iterate the stored (non-zero) entries per row
            data = sub.data
            indptr = sub.indptr
            for r in range(sub.shape[0]):
                seg = data[indptr[r]:indptr[r + 1]]
                if seg.size:
                    out[r] = int(np.count_nonzero((seg != 0) & (seg >= x_min)))
            return out
        block = np.nan_to_num(np.asarray(self.values[:, col_idx], dtype=float), nan=0.0)
        return ((block != 0) & (block >= x_min)).sum(axis=1).astype(np.int64)

    def count_nonzero_rows(self) -> int:
        """Rows with at least one non-zero, non-NaN entry across every sample
        — cheap for a small dense matrix (the gene- or pathway-level ones);
        avoid calling this on a huge sparse one (junction counts / RRS
        scores) without a coverage-style prefilter first."""
        if self.sparse:
            # nnz per row via CSR indptr, no densification
            csr = self.values.tocsr()
            return int(np.count_nonzero(np.diff(csr.indptr)))
        vals = np.nan_to_num(np.asarray(self.values, dtype=float), nan=0.0)
        return int(np.count_nonzero((vals != 0).any(axis=1)))

    def restrict_samples(self, keep: Sequence[str]) -> "Sjdat":
        """A new Sjdat with only the given sample columns, in this matrix's
        existing column order (``keep`` may be a superset, missing some of
        this matrix's own samples, or in any order). Works for both the dense
        and sparse (CSC — column slicing is cheap) cases."""
        keep_set = set(keep)
        idx = [i for i, s in enumerate(self.samples) if s in keep_set]
        return Sjdat(
            kind=self.kind,
            features=self.features,
            samples=[self.samples[i] for i in idx],
            values=self.values[:, idx],
            sparse=self.sparse,
        )


def rows_with_min_supporting_reads(
    features: Iterable[str], junction_counts: "Sjdat", sample_ids: Sequence[str], min_reads: float,
) -> Set[str]:
    """Which of `features` (row labels of some other junction-level sjdat —
    in practice, RRS scores) have a same-named row in `junction_counts`
    (matched by the exact "chr:start-end:strand" label every junction-level
    sjdat shares) whose value exceeds `min_reads` in at least one of
    `sample_ids`. Sample columns are matched by id, not position, since the
    two matrices need not share a column order or even the same full sample
    set. A feature with no corresponding row in `junction_counts` does not
    pass — there is no read-support evidence for it either way."""
    want = [f for f in features if f in junction_counts._row]
    if not want:
        return set()
    row_idx = [junction_counts._row[f] for f in want]
    col_idx = junction_counts.column_indices(list(sample_ids))
    if not col_idx:
        return set()
    block = np.nan_to_num(junction_counts.dense_block(row_idx, col_idx), nan=0.0)
    row_max = block.max(axis=1)
    return {f for f, m in zip(want, row_max) if m > min_reads}


# --------------------------------------------------------------------------- #
def load_sjdat(kind: str, path: Path, tmp_dir: Path) -> Sjdat:
    if kind not in SJDAT_KINDS:
        raise SjdatError(f"unknown sjdat kind {kind!r}")
    # Try the sparse-capable junction reader first regardless of `kind` — a
    # cohort's RRS file is exactly as junction-level (and can be exactly as
    # large) as its junction-count file. Only a matrix whose row names truly
    # aren't chr:start-end:strand at all (a real gene-level matrix) fails
    # this and falls through to the dense reader below.
    try:
        rm = _load_junctions(path, tmp_dir)
        return Sjdat(
            kind=kind,
            features=[str(r) for r in rm._c.rownames],
            samples=list(rm.samples),
            values=rm.values,
            sparse=bool(rm.sparse),
        )
    except _SjvRdsError:
        pass

    # gene-level: dense, gene-ish row labels
    try:
        m = _load_dense(path, tmp_dir)
    except _SjvcMatrixError as e:
        raise SjdatError(f"could not read {path.name!r}: {e}")
    return Sjdat(
        kind=kind,
        features=list(m.features),
        samples=list(m.samples),
        values=np.asarray(m.values, dtype=float),
        sparse=False,
    )
