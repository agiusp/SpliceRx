"""Splice-junction row-name parsing and mapping junctions to selected genes."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence

# chr:start-end:strand   e.g.  chr1:100-200:+   or   1:100-200:-
_ROWNAME_RE = re.compile(r"^(?P<chr>[^:\s]+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[+-])$")


class RownameError(ValueError):
    pass


@dataclass(frozen=True)
class Junction:
    rowname: str
    chrom: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class Gene:
    name: str
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str


def parse_rowname(rowname: str) -> Junction:
    m = _ROWNAME_RE.match(rowname.strip())
    if not m:
        raise RownameError(f"row name {rowname!r} is not of the form chr:start-end:strand")
    a, b = int(m.group("start")), int(m.group("end"))
    lo, hi = (a, b) if a <= b else (b, a)
    return Junction(rowname, m.group("chr"), lo, hi, m.group("strand"))


def norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


# A gene-level matrix row name of the form "<gene_name>:<gene_id>" — the shape
# prepTCGAdata::count_novel_sjs() writes, so distinct gene_ids that share a
# symbol (e.g. the Rfam small RNAs, see annotate_sj()) still get one row each.
_ENSEMBL_GENE_ID_RE = re.compile(r"^ENS[A-Z]*G\d+(\.\d+)?$", re.IGNORECASE)


def gene_name_of_label(label: str) -> str:
    """The gene symbol a gene-level matrix row label stands for: strips a
    trailing ``:<gene_id>`` (count_novel_sjs() output) when present; otherwise
    the label is already a bare symbol (older / hand-built matrices) and is
    returned unchanged. Comma-joined multi-gene rows (an offline-condensed
    matrix) are also returned as-is — callers that care about those parts
    individually split on ``,`` themselves (see geneset.match_matrix_labels)."""
    if "," in label or ":" not in label:
        return label
    name, _, gid = label.rpartition(":")
    return name if _ENSEMBL_GENE_ID_RE.match(gid.strip()) else label


def looks_gene_level(rownames: Sequence[str], *, threshold: float = 0.5) -> bool:
    """True when too few row names parse as chr:start-end:strand for the
    matrix to be raw junctions — i.e. it is almost certainly an
    already-condensed, feature/gene-level matrix (e.g. the output of
    functions/condense_junctions_by_gene.py)."""
    if not rownames:
        return False
    ok = 0
    for rn in rownames:
        try:
            parse_rowname(rn)
            ok += 1
        except RownameError:
            pass
    return (ok / len(rownames)) < threshold


@dataclass
class JunctionGeneMap:
    by_rowname: Dict[str, List[str]]   # rowname -> gene_ids it overlaps (selected set)
    n_bad_rownames: int
    bad_sample: List[str]

    @property
    def rownames(self) -> List[str]:
        return list(self.by_rowname)

    def n_junctions(self) -> int:
        return len(self.by_rowname)

    def genes_present(self) -> List[str]:
        seen: List[str] = []
        s: set = set()
        for gids in self.by_rowname.values():
            for g in gids:
                if g not in s:
                    s.add(g)
                    seen.append(g)
        return seen


def map_junctions_to_genes(
    rownames: Sequence[str],
    genes: Iterable[Gene],
) -> JunctionGeneMap:
    """A junction maps to a gene when they share a chromosome (naming
    normalised) and strand and the junction interval overlaps the gene span.
    A junction can map to more than one selected gene."""
    genes = list(genes)
    by_rowname: Dict[str, List[str]] = {}
    n_bad = 0
    bad: List[str] = []
    for rn in rownames:
        try:
            j = parse_rowname(rn)
        except RownameError:
            n_bad += 1
            if len(bad) < 10:
                bad.append(str(rn))
            continue
        jchr = norm_chrom(j.chrom)
        hits = [
            g.gene_id
            for g in genes
            if g.strand == j.strand
            and norm_chrom(g.chrom) == jchr
            and not (j.end < g.start or j.start > g.end)
        ]
        if hits:
            by_rowname[rn] = hits
    return JunctionGeneMap(by_rowname=by_rowname, n_bad_rownames=n_bad, bad_sample=bad)
