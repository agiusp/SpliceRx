"""Feature selection for SJSurv.

Given ``sjdat`` and the set of samples in the chosen group, a feature (row) F is
kept when:

1. **coverage** — the number of samples whose entry in F is non-zero *and*
   ``>= X`` is at least ``N``. ``N`` is read as a count when ``N >= 1`` and as a
   fraction of the group's samples when ``0 < N < 1``.
2. **magnitude** — folded into (1): only entries ``>= X`` count toward coverage.
3. **variability** — of the rows passing (1)+(2), rank (on ``log1p`` values)
   by descending median absolute deviation (MAD) and keep the top ``n`` — or
   by descending plain variance instead, for an RRS-scores ``sjdat``: RRS
   scores are bounded to [0, 1] and mostly zero, so a row's median is almost
   always exactly 0 and MAD collapses to (near-)0 for most rows, tying too
   many candidates to rank meaningfully.

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


@dataclass
class Ranking:
    """The top_n-independent result of ranking an sjdat's rows for one group:
    every candidate row (after the coverage/restrict/exclude/only_rownames
    filters) reduced to its index into `sjdat.features`, ordered descending
    by score (most variable first — the order `Selection.feature_ids` is
    reported in). Slicing this to any `top_n` (see `select_from_ranking`) is
    then just a list slice + a small `dense_block` call, cheap regardless of
    how big the candidate set was — the point being that a caller re-ranking
    the same group under the same filters but a different `top_n` (e.g. a UI
    slider) should compute the expensive part, this `Ranking`, only once."""
    kind: str
    order: List[int]           # sjdat row indices, every passing candidate, descending by score
    col_idx: List[int]         # sample_ids resolved to sjdat column indices
    used_samples: List[str]    # sample_ids that are actually columns of sjdat
    n_candidates: int
    n_after_coverage: int


def rank_features(
    sjdat: Sjdat,
    sample_ids: List[str],
    *,
    n_min: Optional[float],
    x_min: Optional[float],
    restrict_to: Optional[Set[str]] = None,
    exclude: Optional[Set[str]] = None,
    only_rownames: Optional[Set[str]] = None,
    annotation: Optional[Annotation] = None,
) -> Ranking:
    """Rank `sjdat`'s rows by descending MAD (on `log1p` values, over the
    group's samples) and return every candidate as a `Ranking` — the
    `top_n`-independent part of `select_features`, split out so a caller can
    cache it and answer several different `top_n` requests against the same
    group/filters without repeating the expensive work (gene-overlap
    resolution and the per-row score computation itself). `n_min`/`x_min` are
    an optional coverage prefilter — a row is a candidate only when at least
    `n_min` of the group's samples have an entry that is non-zero *and*
    `>= x_min` — meaningful for a raw per-junction matrix; passing `n_min=None`
    skips the prefilter entirely and ranks every row.

    `restrict_to` (a set of lower-cased gene names, e.g. the protein-coding
    genes of a reference) drops rows not in the set before ranking — a direct
    row-name match for a gene-level sjdat, or (via `annotation`, required in
    that case) which junctions overlap one of those genes for a junction-
    level one. `exclude` is the subtractive counterpart (e.g. a curated
    paralog-family list), applied independently of `restrict_to`. Both mirror
    `sjvc.services.mad.top_features_by_mad`'s own handling exactly.

    `only_rownames`, when given, is a plain allowlist of exact row labels
    (e.g. junctions clearing a minimum-supporting-read-count filter against a
    sibling matrix), applied last with no gene resolution needed."""
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

    if exclude:
        if looks_gene_level(sjdat.features):
            keep = np.array(
                [gene_name_of_label(sjdat.features[i]).strip().lower() not in exclude for i in passing]
            )
        else:
            if annotation is None:
                raise SelectError(
                    "a GENCODE reference is needed to resolve junctions to genes for the "
                    "paralog-family exclusion"
                )
            drop_rownames = junction_rownames_overlapping_genes(
                (sjdat.features[i] for i in passing), annotation, exclude,
            )
            keep = np.array([sjdat.features[i] not in drop_rownames for i in passing])
        passing = passing[keep]
        if passing.size == 0:
            raise SelectError("every candidate feature is in the excluded paralog-family gene set")

    if only_rownames is not None:
        keep = np.array([sjdat.features[i] in only_rownames for i in passing])
        passing = passing[keep]
        if passing.size == 0:
            raise SelectError("no feature meets the minimum supporting-read-count filter")

    block = sjdat.dense_block(passing.tolist(), col_idx)         # (n_passing, n_samples)
    block = np.nan_to_num(block, nan=0.0)
    y = np.log1p(np.abs(block))
    # RRS scores are bounded to [0, 1] and mostly zero, so a row's median is
    # almost always exactly 0 and MAD (robust to a minority of outliers)
    # collapses to (near-)0 for the great majority of rows, tying too many
    # candidates to rank meaningfully — plain variance ranks them better.
    # MAD still suits count-like data (gene/pathway mis-splice counts) with a
    # few strong outliers, so it stays the default everywhere else. Mirrors
    # sjvc.services.mad.top_features_by_mad's identical rationale exactly.
    scores = y.var(axis=1) if sjdat.kind == "rrs_scores" else _mad(y)
    order = np.argsort(scores)[::-1]

    return Ranking(
        kind=sjdat.kind,
        order=passing[order].tolist(),
        col_idx=col_idx,
        used_samples=used_samples,
        n_candidates=sjdat.n_features,
        n_after_coverage=int(passing.size),
    )


def select_from_ranking(sjdat: Sjdat, ranking: Ranking, top_n: int) -> Selection:
    """Slice a `Ranking` (see `rank_features`) down to its `top_n` best rows —
    most-variable-first, matching `select_features`'s long-standing output
    order exactly."""
    if top_n < 1:
        raise SelectError("n (top features) must be >= 1")
    rows = ranking.order[: min(top_n, len(ranking.order))]
    values = np.nan_to_num(sjdat.dense_block(rows, ranking.col_idx), nan=0.0)
    return Selection(
        kind=ranking.kind,
        feature_ids=[sjdat.features[i] for i in rows],
        row_indices=rows,
        sample_ids=ranking.used_samples,
        values=values,
        n_after_coverage=ranking.n_after_coverage,
        n_candidates=ranking.n_candidates,
    )


def select_features(
    sjdat: Sjdat,
    sample_ids: List[str],
    *,
    n_min: Optional[float],
    x_min: Optional[float],
    top_n: int,
    restrict_to: Optional[Set[str]] = None,
    exclude: Optional[Set[str]] = None,
    only_rownames: Optional[Set[str]] = None,
    annotation: Optional[Annotation] = None,
) -> Selection:
    """Rank `sjdat`'s rows by descending MAD and keep the top `top_n` — a
    one-shot convenience that just calls `rank_features` then
    `select_from_ranking`. A caller that will ask for several different
    `top_n` values against the same group/filters (e.g. a UI slider) should
    call those two directly and cache the `Ranking` instead, to avoid
    repeating the expensive ranking work on every change; see their
    docstrings."""
    ranking = rank_features(
        sjdat, sample_ids, n_min=n_min, x_min=x_min, restrict_to=restrict_to,
        exclude=exclude, only_rownames=only_rownames, annotation=annotation,
    )
    return select_from_ranking(sjdat, ranking, top_n)


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
