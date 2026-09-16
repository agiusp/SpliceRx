"""Feature selection for SJSurv.

Given ``sjdat`` and the set of samples in the chosen group, a feature (row) F is
kept when:

1. **coverage** — the number of samples whose entry in F is non-zero *and*
   ``>= X`` is at least ``N``. ``N`` is read as a count when ``N >= 1`` and as a
   fraction of the group's samples when ``0 < N < 1``.
2. **magnitude** — folded into (1): only entries ``>= X`` count toward coverage.
3. **variability** — of the rows passing (1)+(2), rank by descending median
   absolute deviation (MAD, on ``log1p`` values) and keep the top ``n``.

The selected sub-matrix (features x selected-samples, dense) is what the
classifier is trained on.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Set

import numpy as np

from sjvc.services.gencode import Annotation
from sjvc.services.junctions import gene_name_of_label, looks_gene_level
from sjvc.services.mad import junction_rownames_overlapping_genes

from .sjdat import Sjdat


class SelectError(ValueError):
    pass


@dataclass
class Selection:
    kind: str                       # sjdat kind: "junction_counts" | "rrs_scores" | "gene_matrix"
    feature_ids: List[str]          # kept rows, in MAD-rank order (most variable first)
    row_indices: List[int]          # their indices into the sjdat matrix
    sample_ids: List[str]           # columns used (group ∩ matrix ∩ labelled)
    values: np.ndarray              # (n_selected_features, n_samples), dense, NaN->0
    n_after_coverage: int           # passed step (1)+(2)
    n_candidates: int               # rows considered (whole matrix)

    @property
    def feature_noun(self) -> str:
        return "junction" if self.kind == "junction_counts" else "gene" if self.kind == "gene_matrix" else "feature"

    @property
    def n_features(self) -> int:
        return len(self.feature_ids)

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)


def _mad(x: np.ndarray) -> np.ndarray:
    med = np.median(x, axis=1, keepdims=True)
    return np.median(np.abs(x - med), axis=1)


def select_features(
    sjdat: Sjdat,
    sample_ids: List[str],
    *,
    n_min: Optional[float],
    x_min: Optional[float],
    top_n: int,
    restrict_to: Optional[Set[str]] = None,
    annotation: Optional[Annotation] = None,
) -> Selection:
    """Rank `sjdat`'s rows by descending MAD (on `log1p` values, over the
    group's samples) and keep the top `top_n`. `n_min`/`x_min` are an
    optional coverage prefilter — a row is a candidate only when at least
    `n_min` of the group's samples have an entry that is non-zero *and*
    `>= x_min` — meaningful for a raw per-junction matrix; passing `n_min=None`
    skips the prefilter entirely and ranks every row.

    `restrict_to` (a set of lower-cased gene names, e.g. the protein-coding
    genes of a reference) drops rows not in the set before ranking — a direct
    row-name match for a gene-level sjdat, or (via `annotation`, required in
    that case) which junctions overlap one of those genes for a junction-
    level one. Mirrors `sjvc.services.mad.top_features_by_mad`'s own
    `restrict_to`/`annotation` handling exactly."""
    if top_n < 1:
        raise SelectError("n (top features) must be >= 1")
    if n_min is not None and n_min <= 0:
        raise SelectError("N (minimum samples) must be > 0")

    col_idx = sjdat.column_indices(sample_ids)
    used_samples = [s for s in sample_ids if s in sjdat._col]
    if len(used_samples) < 4:
        raise SelectError(
            f"only {len(used_samples)} of the group's samples are columns of this sjdat "
            f"matrix — need at least 4"
        )

    if n_min is None:
        passing = np.arange(sjdat.n_features)
    else:
        need = n_min if n_min >= 1 else max(1.0, n_min * len(used_samples))
        qual = sjdat.nonzero_and_ge_counts(col_idx, x_min or 0.0)   # per row, whole matrix
        passing = np.flatnonzero(qual >= need)
    if passing.size == 0:
        raise SelectError(
            "no feature meets the coverage threshold — lower N or X, or pick a bigger group"
        )

    if restrict_to is not None:
        if looks_gene_level(sjdat.features):
            keep = np.array(
                [gene_name_of_label(sjdat.features[i]).strip().lower() in restrict_to for i in passing]
            )
        else:
            if annotation is None:
                raise SelectError(
                    "a GENCODE reference is needed to resolve junctions to genes for the "
                    "protein-coding filter"
                )
            keep_rownames = junction_rownames_overlapping_genes(
                (sjdat.features[i] for i in passing), annotation, restrict_to,
            )
            keep = np.array([sjdat.features[i] in keep_rownames for i in passing])
        passing = passing[keep]
        if passing.size == 0:
            raise SelectError(
                "none of the candidate features are in the requested protein-coding set — "
                "lower the coverage threshold or turn the filter off"
            )

    block = sjdat.dense_block(passing.tolist(), col_idx)         # (n_passing, n_samples)
    block = np.nan_to_num(block, nan=0.0)
    scores = _mad(np.log1p(np.abs(block)))
    order = np.argsort(scores)[::-1][: min(top_n, passing.size)]

    rows = passing[order].tolist()
    return Selection(
        kind=sjdat.kind,
        feature_ids=[sjdat.features[i] for i in rows],
        row_indices=rows,
        sample_ids=used_samples,
        values=block[order, :],
        n_after_coverage=int(passing.size),
        n_candidates=sjdat.n_features,
    )


def select_from_features(fm, group_sample_ids: List[str], sjdat_kind: str) -> Selection:
    """The gene-set counterpart of `select_features()`: turns an already-built
    `sjvc.services.features.FeatureMatrix` (resolved from typed genes / an
    uploaded list / a pathway, over *every* sample of the active matrix) into
    a `Selection` for one Group — no MAD ranking, no coverage filter, the
    resolved gene set *is* the selection. `fm.values` is (features, samples),
    same orientation as `Selection.values`."""
    col_pos = {s: i for i, s in enumerate(fm.samples)}
    used = [s for s in group_sample_ids if s in col_pos]
    if len(used) < 4:
        raise SelectError(
            f"only {len(used)} of the group's samples are columns of this feature set — "
            f"need at least 4"
        )
    idx = [col_pos[s] for s in used]
    values = np.nan_to_num(np.asarray(fm.values)[:, idx], nan=0.0)
    return Selection(
        kind=sjdat_kind,
        feature_ids=list(fm.feature_ids),
        row_indices=list(range(fm.n_features)),
        sample_ids=used,
        values=values,
        n_after_coverage=fm.n_features,
        n_candidates=fm.n_features,
    )
