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
ranking.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .features import FeatureMatrix
from .gencode import Annotation
from .junctions import Gene, RownameError, norm_chrom, parse_rowname
from .rds import Matrix

Groups = Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray, List[str]]]


def _mad(x: np.ndarray) -> np.ndarray:
    """Per-row median absolute deviation. x: (n_rows, n_samples)."""
    med = np.median(x, axis=1, keepdims=True)
    return np.median(np.abs(x - med), axis=1)


def _top_indices(scores: np.ndarray, top_n: int) -> List[int]:
    n = max(1, min(top_n, len(scores)))
    order = np.argsort(scores)[::-1][:n]
    return sorted(order.tolist())  # keep the matrix's original row order


def top_features_by_mad(
    matrix: Matrix, top_n: int, *, restrict_to: Optional[Set[str]] = None
) -> FeatureMatrix:
    """Rank the matrix's rows by descending MAD and keep the top N.

    A raw junction matrix is ranked over its parseable junction rows only
    (kind="junction"); an already-condensed matrix, whose row names are gene
    symbols rather than chr:start-end:strand, is ranked over every row as-is
    (kind="gene").

    `restrict_to` (a set of lower-cased gene names, e.g. the protein-coding
    genes of a reference) only applies in the gene-level case: rows whose
    name isn't in the set are dropped before ranking."""
    junction_idx: List[int] = []
    for i, rn in enumerate(matrix.features):
        try:
            parse_rowname(rn)
            junction_idx.append(i)
        except RownameError:
            continue

    if junction_idx:
        idx, kind = junction_idx, "junction"
    else:
        idx, kind = list(range(len(matrix.features))), "gene"
        if restrict_to is not None:
            idx = [i for i in idx if matrix.features[i].strip().lower() in restrict_to]
            if not idx:
                raise ValueError(
                    "none of the matrix's genes are in the requested set "
                    "(e.g. no protein-coding gene names matched the reference)"
                )
    if not idx:
        raise ValueError("the uploaded matrix has no rows to rank")

    vals = np.nan_to_num(matrix.values[idx, :].astype(float), nan=0.0)
    scores = _mad(np.log1p(vals))
    top = _top_indices(scores, top_n)

    rows = [idx[i] for i in top]
    return FeatureMatrix(
        kind=kind,
        feature_ids=[matrix.features[i] for i in rows],
        samples=list(matrix.samples),
        values=matrix.values[rows, :],
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
