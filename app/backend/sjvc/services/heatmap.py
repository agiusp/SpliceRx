"""Heatmap matrix: preprocessing, hierarchical clustering, column ordering."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np

from .features import FeatureMatrix

MAX_ROWS = 500


@dataclass
class Dendrogram:
    # line segments in leaf/height space; frontend maps to pixels
    segments: List[List[List[float]]]
    order: List[int]


@dataclass
class HeatmapResult:
    feature_ids: List[str]        # in display (row) order
    samples: List[str]            # in display (column) order
    values: np.ndarray           # (n_display_rows, n_samples) in display order
    row_dendro: Optional[Dendrogram] = None
    col_dendro: Optional[Dendrogram] = None
    warnings: List[str] = field(default_factory=list)


def _linkage_order(M: np.ndarray):
    """Return (leaf_order, dendrogram) for rows of M; falls back to identity."""
    from scipy.cluster.hierarchy import dendrogram, linkage

    if M.shape[0] < 3:
        return list(range(M.shape[0])), None
    metric = "correlation" if M.shape[1] > 2 else "euclidean"
    Z = linkage(M, method="average", metric=metric)
    dd = dendrogram(Z, no_plot=True)
    order = [int(i) for i in dd["leaves"]]
    segments = [
        [[float(x) for x in xs], [float(y) for y in ys]]
        for xs, ys in zip(dd["icoord"], dd["dcoord"])
    ]
    return order, Dendrogram(segments=segments, order=order)


def build(
    fm: FeatureMatrix,
    *,
    row_zscore: bool,
    order: str = "cluster",             # "cluster" | "group"
    group_values: Optional[Sequence[Optional[str]]] = None,
) -> HeatmapResult:
    M = np.nan_to_num(fm.values.astype(float), nan=0.0)
    ids = list(fm.feature_ids)
    samples = list(fm.samples)
    warnings: List[str] = []

    M = np.log1p(M)

    # drop zero-variance rows
    keep = M.var(axis=1) > 1e-12
    if (~keep).any():
        M = M[keep, :]
        ids = [i for i, k in zip(ids, keep) if k]

    # row legibility cap: keep highest-variance rows
    if M.shape[0] > MAX_ROWS:
        top = np.argsort(M.var(axis=1))[::-1][:MAX_ROWS]
        top.sort()
        M = M[top, :]
        ids = [ids[i] for i in top]
        warnings.append(
            f"showing the {MAX_ROWS} highest-variance features of {fm.n_features}; "
            "tighten the gene set or condense per gene for the full picture"
        )

    if M.shape[0] == 0:
        raise ValueError("no variable features to show")

    if row_zscore:
        mu = M.mean(axis=1, keepdims=True)
        sd = M.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        M = (M - mu) / sd

    # row clustering (always)
    row_order, row_dendro = _linkage_order(M)

    # column ordering
    col_dendro = None
    if order == "group" and group_values is not None:
        keyed = list(enumerate(group_values))
        keyed.sort(key=lambda t: (t[1] is None, str(t[1])))
        col_order = [i for i, _ in keyed]
    else:
        col_order, col_dendro = _linkage_order(M.T)

    M = M[np.ix_(row_order, col_order)]
    ids = [ids[i] for i in row_order]
    samples = [samples[i] for i in col_order]

    return HeatmapResult(
        feature_ids=ids,
        samples=samples,
        values=M,
        row_dendro=row_dendro,
        col_dendro=col_dendro,
        warnings=warnings,
    )
