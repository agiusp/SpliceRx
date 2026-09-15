"""Classify each drawn junction against the GENCODE transcript models for the
query gene into one of six categories (plus `annotated`).

A junction's genomic interval [start, end] is treated as the *intron* it
represents. Splice sites are strand-aware:

    + strand:  donor (5' ss) = intron low coord,  acceptor (3' ss) = high coord
    - strand:  donor          = intron high coord, acceptor          = low coord
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Set, Tuple

from .junctions import Junction


@dataclass(frozen=True)
class Exon:
    start: int
    end: int
    kind: str  # "CDS" | "UTR"


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    strand: str
    exons: Tuple[Exon, ...]  # any order; sorted internally on use


CATEGORIES = (
    "annotated",
    "exon_skipping",
    "alt_5p",
    "alt_3p",
    "isoform_switch",
    "novel_exon",
    "novel",
)


def _introns_of(tx: Transcript) -> List[Tuple[int, int]]:
    ex = sorted(tx.exons, key=lambda e: (e.start, e.end))
    out = []
    for a, b in zip(ex, ex[1:]):
        s, e = a.end + 1, b.start - 1
        if e >= s:
            out.append((s, e))
    return out


class GeneClassifier:
    def __init__(self, transcripts: Iterable[Transcript], tol: int = 0):
        self.tol = tol
        self.transcripts: List[Transcript] = list(transcripts)

        # (intronStart, intronEnd) -> set of transcript ids
        self.annotated_introns: Dict[Tuple[int, int], Set[str]] = {}
        # splice-site coord -> set of transcript ids (strand-aware)
        self.donor_sites: Dict[int, Set[str]] = {}
        self.acceptor_sites: Dict[int, Set[str]] = {}
        # transcript id -> sorted exon list
        self._tx_exons: Dict[str, List[Exon]] = {}

        for tx in self.transcripts:
            self._tx_exons[tx.transcript_id] = sorted(tx.exons, key=lambda e: (e.start, e.end))
            for (s, e) in _introns_of(tx):
                self.annotated_introns.setdefault((s, e), set()).add(tx.transcript_id)
                if tx.strand == "-":
                    donor, acceptor = e, s
                else:
                    donor, acceptor = s, e
                self.donor_sites.setdefault(donor, set()).add(tx.transcript_id)
                self.acceptor_sites.setdefault(acceptor, set()).add(tx.transcript_id)

        self._all_introns = sorted(self.annotated_introns)

    # -- tolerance-aware lookups -------------------------------------------------
    def _match_intron(self, s: int, e: int) -> Set[str]:
        best: Set[str] = set()
        for (a, b), txs in self.annotated_introns.items():
            if abs(a - s) <= self.tol and abs(b - e) <= self.tol:
                best |= txs
        return best

    def _match_site(self, sites: Dict[int, Set[str]], coord: int) -> Set[str]:
        out: Set[str] = set()
        for c, txs in sites.items():
            if abs(c - coord) <= self.tol:
                out |= txs
        return out

    # -- classification --------------------------------------------------------
    def classify(self, j: Junction) -> str:
        s, e = j.start, j.end

        if self._match_intron(s, e):
            return "annotated"

        if j.strand == "-":
            donor_coord, acceptor_coord = e, s
        else:
            donor_coord, acceptor_coord = s, e

        donor_txs = self._match_site(self.donor_sites, donor_coord)
        acceptor_txs = self._match_site(self.acceptor_sites, acceptor_coord)
        donor_annot = bool(donor_txs)
        acceptor_annot = bool(acceptor_txs)

        if donor_annot and acceptor_annot:
            # both splice sites are known; is there a single transcript in which
            # this junction skips >= 1 whole exon?
            for txid in donor_txs & acceptor_txs:
                exons = self._tx_exons[txid]
                inside = [x for x in exons if x.start > s and x.end < e]
                if inside:
                    return "exon_skipping"
            return "isoform_switch"

        if donor_annot != acceptor_annot:
            # exactly one end is annotated -> the other end is the alternative site
            return "alt_3p" if donor_annot else "alt_5p"

        # neither end annotated
        for (a, b) in self._all_introns:
            if a - self.tol <= s and e <= b + self.tol and (a, b) != (s, e):
                return "novel_exon"
        return "novel"
