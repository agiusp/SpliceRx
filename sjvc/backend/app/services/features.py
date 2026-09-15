"""Build the feature matrix fed to the heatmap / projection.

Stages:
1. select the junction sub-matrix (junctions overlapping the selected genes)
2. optionally condense to one row per gene:
       C[gene, sample] = # of that gene's junctions with count > 0   (NA -> 0)
   (the user's R:  apply(J > 0, 2, sum))
3. preprocess: NA -> 0, drop zero-variance rows, log1p, and (for projection)
   z-score each feature across samples.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .junctions import JunctionGeneMap
from .rds import Matrix

CONDENSE_SUGGEST_ABOVE = 2000


@dataclass
class FeatureMatrix:
    kind: str                    # "junction" | "gene"
    feature_ids: List[str]
    samples: List[str]
    values: np.ndarray          # (n_features, n_samples), raw counts (NaN allowed)
    n_junctions: int             # size of the selected junction subset (pre-condense)
    dropped_zero_var: int = 0
    _gene_labels: Dict[str, str] = field(default_factory=dict)  # gene_id -> symbol

    @property
    def n_features(self) -> int:
        return len(self.feature_ids)

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    def label(self, feature_id: str) -> str:
        return self._gene_labels.get(feature_id, feature_id)


def select_and_maybe_condense(
    matrix: Matrix,
    jgmap: JunctionGeneMap,
    *,
    condense: bool,
    gene_labels: Optional[Dict[str, str]] = None,
) -> FeatureMatrix:
    rownames = jgmap.rownames
    n_junctions = len(rownames)
    if n_junctions == 0:
        raise ValueError("no junctions in the uploaded matrix fall within the selected genes")

    sub = matrix.submatrix(rownames)  # (n_junctions, n_samples)

    if not condense:
        return FeatureMatrix(
            kind="junction",
            feature_ids=list(rownames),
            samples=list(matrix.samples),
            values=sub,
            n_junctions=n_junctions,
            _gene_labels=gene_labels or {},
        )

    # condense: per gene, count junctions with > 0 reads per sample
    gene_ids = jgmap.genes_present()
    expressed = np.nan_to_num(sub, nan=0.0) > 0        # (n_junctions, n_samples) bool
    rows = []
    for gid in gene_ids:
        mask = np.array([gid in jgmap.by_rowname[rn] for rn in rownames])
        rows.append(expressed[mask, :].sum(axis=0))
    values = np.vstack(rows).astype(float) if rows else np.zeros((0, matrix.n_samples))
    return FeatureMatrix(
        kind="gene",
        feature_ids=gene_ids,
        samples=list(matrix.samples),
        values=values,
        n_junctions=n_junctions,
        _gene_labels=gene_labels or {},
    )


def select_by_labels(matrix: Matrix, labels: List[str]) -> FeatureMatrix:
    """Gene-level matrix: keep exactly the named rows, in matrix order."""
    idx = [matrix._row[lab] for lab in labels if lab in matrix._row]
    if not idx:
        raise ValueError("none of the selected genes are rows in the matrix")
    return FeatureMatrix(
        kind="gene",
        feature_ids=[matrix.features[i] for i in idx],
        samples=list(matrix.samples),
        values=matrix.values[idx, :],
        n_junctions=0,
    )


@dataclass
class Prepared:
    X: np.ndarray               # samples x features  (for projection: standardized)
    feature_ids: List[str]
    samples: List[str]
    dropped: int


def preprocess(fm: FeatureMatrix, *, standardize: bool) -> Prepared:
    M = np.nan_to_num(fm.values.astype(float), nan=0.0)   # features x samples
    var = M.var(axis=1)
    keep = var > 1e-12
    dropped = int((~keep).sum())
    M = M[keep, :]
    ids = [f for f, k in zip(fm.feature_ids, keep) if k]
    if M.shape[0] == 0:
        raise ValueError("every feature is constant after filtering — widen the gene set")

    M = np.log1p(M)
    X = M.T  # samples x features

    if standardize:
        mu = X.mean(axis=0)
        sd = X.std(axis=0)
        sd[sd == 0] = 1.0
        X = (X - mu) / sd

    return Prepared(X=X, feature_ids=ids, samples=list(fm.samples), dropped=dropped)
