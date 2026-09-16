"""Pick the top-N most variable features by median absolute deviation (MAD),
as an alternative to choosing a gene set explicitly.

- For a raw junction matrix (rows named chr:start-end:strand), the ranked
  features are the splice junctions themselves — every parseable row.
- For an already-condensed matrix (gene symbols as row names, e.g. the
  output of functions/condense_junctions_by_gene.py), every row is ranked
  as-is, and the result is reported as gene-level features.

`top_genes_by_mad` is the older in-app path that condensed genome-wide on the
fly; it is kept for the API but no longer reachable from the UI (that work
moved to the offline CLI).

MAD is computed on `log1p(NA -> 0)` counts, matching the rest of the pipeline's
preprocessing, so a handful of extreme high-count features don't dominate the
ranking. A truly constant row (zero *variance* — the same value in every
sample) carries no signal at all, and `top_features_by_mad` excludes those
before taking the top N rather than letting them pad out a request for more
features than actually vary — otherwise "top 500 by MAD" on a matrix with
only, say, 325 non-constant rows would silently return 500 anyway (175 of
them constant filler), only to have `features.preprocess()` drop those same
175 again downstream, unexpectedly short of the 500 asked for.

This is deliberately a *variance* check, not "MAD == 0": MAD is robust to a
minority of outliers, so a row where one value repeats for just over half the
samples (e.g. five 6's and one 5) can legitimately have MAD 0 while still
varying — that row is real, just not a top-ranked one, and stays eligible.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np

from .features import FeatureMatrix
from .gencode import Annotation
from .junctions import Gene, RownameError, gene_name_of_label, norm_chrom, parse_rowname
from .rds import Matrix

Groups = Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray, List[str]]]

# Matches features.preprocess()'s own "constant" threshold exactly, so a row
# top_features_by_mad() keeps is never one preprocess() would immediately
# turn around and drop.
_CONSTANT_VAR_EPS = 1e-12


def _mad(x: np.ndarray) -> np.ndarray:
    """Per-row median absolute deviation. x: (n_rows, n_samples)."""
    med = np.median(x, axis=1, keepdims=True)
    return np.median(np.abs(x - med), axis=1)


def _coverage_counts_dense(vals: np.ndarray, x_min: float) -> np.ndarray:
    """Per row, the number of entries that are non-zero *and* ``>= x_min``.
    ``vals`` is already a dense ``(n_rows, n_samples)`` block."""
    vals = np.nan_to_num(vals, nan=0.0)
    return ((vals != 0) & (vals >= x_min)).sum(axis=1).astype(np.int64)


def _coverage_counts_sparse_rows(values, row_idx: List[int], x_min: float) -> np.ndarray:
    """Sparse-matrix equivalent of `_coverage_counts_dense`, one row at a time
    (mirrors `Sjdat.nonzero_and_ge_counts` / `_mad_and_var_sparse_rows`'s style
    — never densifies more than one row at once)."""
    sub = values[row_idx, :].tocsr()
    data, indptr = sub.data, sub.indptr
    out = np.zeros(len(row_idx), dtype=np.int64)
    for r in range(len(row_idx)):
        seg = data[indptr[r]:indptr[r + 1]]
        if seg.size:
            out[r] = int(np.count_nonzero((seg != 0) & (seg >= x_min)))
    return out


def _mad_and_var_sparse_rows(values, row_idx: List[int], n_samples: int):
    """Per-row (MAD, variance) on log1p values, one row at a time, for a
    CSC/CSR sparse matrix — needed when `row_idx` is millions of rows
    (ranking every raw junction of a whole-cohort sparse matrix): densifying
    even one row is cheap (n_samples is at most a few hundred), densifying
    all of them at once is not. Numerically identical to computing
    `_mad(np.log1p(dense_block))` and `np.log1p(dense_block).var(axis=1)`."""
    sub = values[row_idx, :].tocsr()
    data, indices, indptr = sub.data, sub.indices, sub.indptr
    mad_out = np.empty(len(row_idx), dtype=float)
    var_out = np.empty(len(row_idx), dtype=float)
    row = np.empty(n_samples, dtype=float)
    for r in range(len(row_idx)):
        row.fill(0.0)
        lo, hi = indptr[r], indptr[r + 1]
        row[indices[lo:hi]] = data[lo:hi]
        y = np.log1p(row)
        med = np.median(y)
        mad_out[r] = np.median(np.abs(y - med))
        var_out[r] = y.var()
    return mad_out, var_out


def _top_indices(scores: np.ndarray, top_n: int) -> List[int]:
    n = max(1, min(top_n, len(scores)))
    order = np.argsort(scores)[::-1][:n]
    return sorted(order.tolist())  # keep the matrix's original row order


def refine_by_mad(fm: FeatureMatrix, top_n: int) -> FeatureMatrix:
    """Narrow an already-resolved `FeatureMatrix` (from typed/uploaded/pathway
    gene-set resolution — see `sjvc.api.routes.build_features()`) down to its
    `top_n` most variable rows by descending MAD — the same ranking
    `top_features_by_mad` uses, just over a small, already-selected subset
    rather than the whole matrix. A no-op when the resolved set is already
    `<= top_n`.

    Truly-constant (zero-*variance*) rows are excluded before ranking, same
    as `top_features_by_mad` and for the same reason: MAD's robustness means
    a row can score MAD == 0 while still varying (a majority-repeated value),
    and that row is real signal that must never lose out to a genuinely
    constant one just because both nominally tie at MAD == 0."""
    if top_n >= fm.n_features:
        return fm
    vals = np.nan_to_num(np.asarray(fm.values, dtype=float), nan=0.0)
    y = np.log1p(vals)
    scores = _mad(y)
    variances = y.var(axis=1)
    idx = np.flatnonzero(variances > _CONSTANT_VAR_EPS)
    if idx.size == 0:
        idx = np.arange(fm.n_features)  # every row is constant — nothing to prefer
    order = sorted(idx[np.argsort(scores[idx])[::-1][:top_n]].tolist())
    return FeatureMatrix(
        kind=fm.kind, feature_ids=[fm.feature_ids[i] for i in order],
        samples=fm.samples, values=vals[order, :], n_junctions=fm.n_junctions,
        _gene_labels=fm._gene_labels,
    )


def junction_rownames_overlapping_genes(
    rownames: Iterable[str], annotation: Annotation, gene_names: Set[str],
) -> Set[str]:
    """Which of `rownames` (chr:start-end:strand) overlap at least one gene
    in `annotation` whose name (lower-cased) is in `gene_names` — the
    junction-level equivalent of matching a gene-level matrix's row names
    directly against a gene-name set (used by `top_features_by_mad` below,
    and by SJSurv's own `services.select.select_features`, to restrict MAD
    ranking to junctions of protein-coding genes)."""
    genes = _load_all_genes(annotation)
    groups = _group_genes(genes)
    gene_name_by_id = {g.gene_id: g.name for g in genes}
    jgmap = _map_all_junctions(list(rownames), groups)
    return {
        rn for rn, gids in jgmap.items()
        if any(gene_name_by_id.get(gid, "").strip().lower() in gene_names for gid in gids)
    }


def top_features_by_mad(
    matrix: Matrix, top_n: int, *, restrict_to: Optional[Set[str]] = None,
    n_min: Optional[float] = None, x_min: float = 0.0,
    annotation: Optional[Annotation] = None,
) -> FeatureMatrix:
    """Rank the matrix's rows by descending MAD and keep the top N.

    A raw junction matrix is ranked over its parseable junction rows only
    (kind="junction"); an already-condensed matrix, whose row names are gene
    symbols rather than chr:start-end:strand, is ranked over every row as-is
    (kind="gene").

    `restrict_to` (a set of lower-cased gene names, e.g. the protein-coding
    genes of a reference) drops rows not in the set before ranking. In the
    gene-level case this is a direct row-name match; in the junction-level
    case it needs `annotation` to resolve which gene(s) each junction
    overlaps (`junction_rownames_overlapping_genes` above) — the cohort's own
    junction-metadata gene lookup, when loaded, has no gene-biotype
    information to filter on, so a live GENCODE reference is required either
    way.

    `n_min`/`x_min` are an optional coverage prefilter, applied before
    ranking: a row is a candidate only when at least `n_min` samples (a count
    if `>= 1`, else a fraction of all samples) have an entry that is
    non-zero *and* `>= x_min`. Meaningful mainly for a raw per-junction
    matrix; `n_min=None` (the default) skips the prefilter entirely, ranking
    every row exactly as before this parameter existed."""
    junction_idx: List[int] = []
    for i, rn in enumerate(matrix.features):
        try:
            parse_rowname(rn)
            junction_idx.append(i)
        except RownameError:
            continue

    if junction_idx:
        idx, kind = junction_idx, "junction"
        if restrict_to is not None:
            if annotation is None:
                raise ValueError(
                    "a GENCODE reference is needed to resolve junctions to genes for the "
                    "protein-coding filter"
                )
            keep_rownames = junction_rownames_overlapping_genes(
                (matrix.features[i] for i in idx), annotation, restrict_to,
            )
            idx = [i for i in idx if matrix.features[i] in keep_rownames]
            if not idx:
                raise ValueError(
                    "none of the matrix's junctions overlap a protein-coding gene in the reference"
                )
    else:
        idx, kind = list(range(len(matrix.features))), "gene"
        if restrict_to is not None:
            idx = [
                i for i in idx
                if gene_name_of_label(matrix.features[i]).strip().lower() in restrict_to
            ]
            if not idx:
                raise ValueError(
                    "none of the matrix's genes are in the requested set "
                    "(e.g. no protein-coding gene names matched the reference)"
                )
    if not idx:
        raise ValueError("the uploaded matrix has no rows to rank")

    if n_min is not None:
        need = n_min if n_min >= 1 else max(1.0, n_min * matrix.n_samples)
        if getattr(matrix, "sparse", False):
            qual = _coverage_counts_sparse_rows(matrix.values, idx, x_min)
        else:
            qual = _coverage_counts_dense(matrix.values[idx, :].astype(float), x_min)
        idx = [i for i, keep in zip(idx, qual >= need) if keep]
        if not idx:
            raise ValueError(
                "no feature meets the coverage threshold — lower N or X, or rank without it"
            )

    if getattr(matrix, "sparse", False):
        # never densify a (possibly millions-of-rows) block at once — rank
        # from one row at a time, then only densify the small top-N kept
        scores, variances = _mad_and_var_sparse_rows(matrix.values, idx, matrix.n_samples)
    else:
        vals = np.nan_to_num(matrix.values[idx, :].astype(float), nan=0.0)
        y = np.log1p(vals)
        scores = _mad(y)
        variances = y.var(axis=1)

    # exclude truly constant rows (zero variance) before ranking — see the
    # module docstring for why letting them pad out the top-N is wrong, and
    # why this checks variance rather than the MAD score itself
    varies = variances > _CONSTANT_VAR_EPS
    idx = [i for i, keep in zip(idx, varies) if keep]
    scores = scores[varies]
    if not idx:
        raise ValueError("every row is constant (zero variance) — nothing to rank")

    top = _top_indices(scores, top_n)
    rows = [idx[i] for i in top]
    values = (
        matrix.dense_block(rows, list(range(matrix.n_samples)))
        if getattr(matrix, "sparse", False)
        else matrix.values[rows, :]
    )

    return FeatureMatrix(
        kind=kind,
        feature_ids=[matrix.features[i] for i in rows],
        samples=list(matrix.samples),
        values=values,
        n_junctions=len(rows) if kind == "junction" else 0,
    )


def _load_all_genes(annotation: Annotation) -> List[Gene]:
    con = annotation._con()
    try:
        rows = con.execute("SELECT name, gene_id, chrom, start, end, strand FROM genes").fetchall()
    finally:
        con.close()
    best: Dict[str, Gene] = {}
    for r in rows:
        g = Gene(
            name=r["name"], gene_id=r["gene_id"], chrom=r["chrom"],
            start=r["start"], end=r["end"], strand=r["strand"],
        )
        prev = best.get(g.gene_id)
        if prev is None or (g.end - g.start) > (prev.end - prev.start):
            best[g.gene_id] = g
    return list(best.values())


def _group_genes(genes: List[Gene]) -> Groups:
    tmp: Dict[Tuple[str, str], List[Gene]] = defaultdict(list)
    for g in genes:
        tmp[(norm_chrom(g.chrom), g.strand)].append(g)
    groups: Groups = {}
    for key, gs in tmp.items():
        groups[key] = (
            np.array([g.start for g in gs]),
            np.array([g.end for g in gs]),
            [g.gene_id for g in gs],
        )
    return groups


def _map_all_junctions(rownames: List[str], groups: Groups) -> Dict[str, List[str]]:
    """Vectorised per-chromosome overlap test against every gene in the
    reference (not a preselected few) — needed to rank genes genome-wide."""
    by_rowname: Dict[str, List[str]] = {}
    for rn in rownames:
        try:
            j = parse_rowname(rn)
        except RownameError:
            continue
        grp = groups.get((norm_chrom(j.chrom), j.strand))
        if not grp:
            continue
        starts, ends, ids = grp
        hits = [ids[i] for i in np.nonzero((starts <= j.end) & (ends >= j.start))[0]]
        if hits:
            by_rowname[rn] = hits
    return by_rowname


def top_genes_by_mad(matrix: Matrix, annotation: Annotation, top_n: int) -> FeatureMatrix:
    genes = _load_all_genes(annotation)
    groups = _group_genes(genes)
    jgmap = _map_all_junctions(matrix.features, groups)
    if not jgmap:
        raise ValueError("no matrix junction overlaps any gene in the reference")

    rownames = list(jgmap.keys())
    gene_ids = sorted({gid for gids in jgmap.values() for gid in gids})
    expressed = np.nan_to_num(matrix.submatrix(rownames), nan=0.0) > 0  # (n_junctions, n_samples)

    condensed = np.zeros((len(gene_ids), matrix.n_samples))
    for gi, gid in enumerate(gene_ids):
        mask = np.array([gid in jgmap[rn] for rn in rownames])
        condensed[gi, :] = expressed[mask, :].sum(axis=0)

    scores = _mad(np.log1p(condensed))
    top = _top_indices(scores, top_n)

    labels = {g.gene_id: g.name for g in genes}
    return FeatureMatrix(
        kind="gene",
        feature_ids=[gene_ids[i] for i in top],
        samples=list(matrix.samples),
        values=condensed[top, :],
        n_junctions=len(rownames),
        _gene_labels=labels,
    )
