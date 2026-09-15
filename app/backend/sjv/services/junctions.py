"""Splice-junction row-name parsing, locus filtering and count scaling."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import List, Optional, Sequence

# chr:start-end:strand  e.g.  chr1:100-200:+   or   1:100-200:-
_ROWNAME_RE = re.compile(r"^(?P<chr>[^:\s]+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[+-])$")


@dataclass(frozen=True)
class Junction:
    rowname: str
    chrom: str
    start: int  # always <= end (genomic intron start, inclusive)
    end: int    # always >= start (genomic intron end, inclusive)
    strand: str

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


@dataclass(frozen=True)
class Gene:
    name: str
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str


class RownameError(ValueError):
    pass


def parse_rowname(rowname: str) -> Junction:
    """Parse a `chr:start-end:strand` row name into a Junction.

    Raises RownameError if the string does not match the expected format.
    start/end are normalised so that start <= end.
    """
    m = _ROWNAME_RE.match(rowname.strip())
    if not m:
        raise RownameError(f"row name {rowname!r} is not of the form chr:start-end:strand")
    a, b = int(m.group("start")), int(m.group("end"))
    lo, hi = (a, b) if a <= b else (b, a)
    return Junction(
        rowname=rowname,
        chrom=m.group("chr"),
        start=lo,
        end=hi,
        strand=m.group("strand"),
    )


def norm_chrom(c: str) -> str:
    """Strip a leading ``chr`` so ``chr1`` and ``1`` compare equal."""
    return c[3:] if c.lower().startswith("chr") else c


_norm_chrom = norm_chrom  # backwards-compatible alias


def filter_for_gene(junctions: Sequence[Junction], gene: Gene) -> List[Junction]:
    """Keep junctions on the same chromosome and strand as the gene whose
    interval overlaps the gene span. Chromosome naming (`chr1` vs `1`) is
    normalised before comparison."""
    gchr = _norm_chrom(gene.chrom)
    out = []
    for j in junctions:
        if _norm_chrom(j.chrom) != gchr:
            continue
        if j.strand != gene.strand:
            continue
        if j.end < gene.start or j.start > gene.end:
            continue
        out.append(j)
    return out


@dataclass(frozen=True)
class ScaledArc:
    junction: Junction
    count: float
    height: float


def scale_counts(
    junctions: Sequence[Junction],
    counts: Sequence[float],
    min_reads: float = 0.0,
) -> List[ScaledArc]:
    """Keep junctions whose count (or group median) is > 0 and >= `min_reads`,
    then compute the arc height that encodes read support (SCANVIS's log2 rule).
    All arcs are drawn at the same line width, so height is the only read-support
    cue.

    height = log2(c)/max(log2(c)) + 0.1, clamped to a constant when the input
    is degenerate (<=1 junction, or all counts equal).
    """
    floor = max(min_reads, 0.0)
    pairs = [
        (j, float(c))
        for j, c in zip(junctions, counts)
        if c is not None
        and not (isinstance(c, float) and math.isnan(c))
        and float(c) > 0
        and float(c) >= floor
    ]
    if not pairs:
        return []

    counts_equal = len({round(c, 6) for _, c in pairs}) == 1
    log_vals = {c: math.log2(c) for _, c in pairs}
    max_log = max(log_vals.values())

    arcs: List[ScaledArc] = []
    for j, c in pairs:
        if len(pairs) <= 1 or counts_equal or max_log <= 0:
            height = 0.6
        else:
            height = log_vals[c] / max_log + 0.1
        arcs.append(ScaledArc(junction=j, count=c, height=round(height, 4)))
    return arcs
