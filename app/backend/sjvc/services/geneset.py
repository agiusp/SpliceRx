"""Resolve a user gene set (typed / uploaded list / pathway).

A junction-level matrix (junction counts, RRS scores) resolves symbols one of
two ways: if the cohort's junction metadata table is loaded, straight from
its precomputed junction/gene overlaps (fast, no GENCODE needed — see
``junction_metadata.JunctionGeneIndex.resolve``); otherwise to GENCODE genes,
whose loci are then used to find overlapping junctions live. A gene-level
matrix instead matches symbols straight against the matrix's own row labels.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .gencode import Annotation
from .junctions import Gene, JunctionGeneMap, gene_name_of_label

_SPLIT = re.compile(r"[\s,;]+")
_MAX_GENES = 2000


@dataclass
class GeneSet:
    matched: List[Gene] = field(default_factory=list)          # junction-matrix path (GENCODE genes)
    matched_names: List[str] = field(default_factory=list)     # gene-level-matrix path (row labels)
    # junction-matrix path via the junction-metadata table (no GENCODE needed):
    # its junction/gene overlaps, already narrowed to just the requested genes
    matched_junction_map: Optional[JunctionGeneMap] = None
    matched_gene_labels: Dict[str, str] = field(default_factory=dict)  # gene_id -> symbol, for the above
    unmatched: List[str] = field(default_factory=list)
    source: str = ""
    warnings: List[str] = field(default_factory=list)

    @property
    def n_genes(self) -> int:
        return len(self.matched) or len(self.matched_names) or len(self.matched_gene_labels)


def _label_parts(label: str) -> List[str]:
    """Every gene symbol a matrix row label stands for. Two shapes:

    - an offline-condensed multi-gene row, several symbols joined by commas
      (e.g. 'CFH,CFHR1,CFHR3') — split into its parts;
    - anything else — a single symbol: either already bare, or a
      count_novel_sjs()-style 'NAME:gene_id' row, unwrapped by
      gene_name_of_label()."""
    if "," in label:
        return [p.strip() for p in str(label).split(",") if p.strip()]
    return [gene_name_of_label(str(label).strip())]


def match_matrix_labels(
    symbols: Sequence[str], matrix_labels: Sequence[str], *, prefix: bool = False
) -> Tuple[List[str], List[str], List[str]]:
    """Match requested gene symbols against a gene-level matrix's row labels
    (case-insensitive). A symbol matches its own single-gene row if there is
    one; otherwise it falls back to every comma-joined row that contains it
    (some matrices merge neighbouring genes into one row). Prefix mode always
    matches on label parts. Returns (matched row labels in matrix order,
    unmatched symbols, warnings)."""
    # symbol (lower) -> every row whose *entire* gene-set is just that symbol —
    # distinct gene_ids that happen to share a name (e.g. the Rfam small RNAs)
    # are separate single-gene rows and all belong here, no "merged" warning.
    exact: Dict[str, List[str]] = {}
    by_part: Dict[str, List[str]] = {}
    ordered: List[str] = []
    for lab in matrix_labels:
        ordered.append(lab)
        parts = _label_parts(lab)
        for part in parts:
            by_part.setdefault(part.lower(), []).append(lab)
        if len(parts) == 1:
            exact.setdefault(parts[0].lower(), []).append(lab)

    matched: set = set()
    unmatched: List[str] = []
    merged_only: List[str] = []
    for sym in symbols:
        q = sym.strip().lower()
        if prefix:
            hits = {lab for part, labs in by_part.items() if part.startswith(q) for lab in labs}
        elif q in exact:
            hits = set(exact[q])
        elif q in by_part:
            hits = set(by_part[q])
            merged_only.append(sym)
        else:
            hits = set()
        if hits:
            matched.update(hits)
        else:
            unmatched.append(sym)

    warnings: List[str] = []
    if merged_only:
        warnings.append(
            f"{len(merged_only)} gene(s) have no single-gene row in this matrix "
            f"({', '.join(merged_only[:6])}{'…' if len(merged_only) > 6 else ''}) — "
            f"kept the combined rows that contain them"
        )
    return [lab for lab in ordered if lab in matched], unmatched, warnings


def parse_symbols(text: str) -> List[str]:
    seen: List[str] = []
    s: set = set()
    for tok in _SPLIT.split(text or ""):
        tok = tok.strip()
        if tok and tok.upper() not in s:
            s.add(tok.upper())
            seen.append(tok)
    return seen[:_MAX_GENES]


_PREFIX_CAP = 500


def resolve(
    symbols: List[str],
    annotation: Annotation,
    source: str,
    *,
    prefix: bool = False,
) -> GeneSet:
    matched: List[Gene] = []
    unmatched: List[str] = []
    warnings: List[str] = []
    seen_ids: set = set()

    def add(g: Gene) -> None:
        if g.gene_id not in seen_ids:
            seen_ids.add(g.gene_id)
            matched.append(g)

    for sym in symbols:
        if prefix:
            hits = annotation.genes_by_prefix(sym)
            if not hits:
                unmatched.append(sym)
            for g in hits:
                add(g)
        else:
            g = annotation.get_gene(sym)
            if g is None:
                unmatched.append(sym)
            else:
                add(g)

    if len(matched) > _PREFIX_CAP:
        warnings.append(
            f"{len(matched)} genes matched — kept the first {_PREFIX_CAP}; use a longer prefix"
        )
        matched = matched[:_PREFIX_CAP]
    return GeneSet(matched=matched, unmatched=unmatched, source=source, warnings=warnings)
