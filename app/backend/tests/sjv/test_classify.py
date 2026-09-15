"""Classification on a synthetic 3-transcript gene, tested on + and - strand."""
import pytest

from sjv.services.classify import Exon, GeneClassifier, Transcript
from sjv.services.junctions import parse_rowname

# + strand gene: exons E1..E5 at 1000s..5000s
#   tx1: E1,E2,E3     tx2: E2,E3,E4     tx3: E3,E4,E5
E = {
    1: (1000, 1200), 2: (2000, 2200), 3: (3000, 3200),
    4: (4000, 4200), 5: (5000, 5200),
}


def _tx(txid, nums, strand):
    return Transcript(txid, strand, tuple(Exon(*E[n], "CDS") for n in nums))


def plus_gene():
    return [
        _tx("t1", [1, 2, 3], "+"),
        _tx("t2", [2, 3, 4], "+"),
        _tx("t3", [3, 4, 5], "+"),
    ]


def minus_gene():
    return [
        _tx("t1", [1, 2, 3], "-"),
        _tx("t2", [2, 3, 4], "-"),
        _tx("t3", [3, 4, 5], "-"),
    ]


@pytest.mark.parametrize(
    "rowname,expected",
    [
        ("chr1:1201-1999:+", "annotated"),      # tx1 intron 1
        ("chr1:2201-2999:+", "annotated"),      # tx1/tx2 shared intron
        ("chr1:2201-3999:+", "exon_skipping"),  # skips E3 within tx2
        ("chr1:1150-1999:+", "alt_5p"),         # novel donor, annotated acceptor
        ("chr1:2201-2500:+", "alt_3p"),         # annotated donor, novel acceptor
        ("chr1:1201-4999:+", "isoform_switch"), # donor only t1, acceptor only t3
        ("chr1:2400-2600:+", "novel_exon"),     # inside intron 2201-2999
        ("chr1:1500-3500:+", "novel"),          # spans exons, no annotated site
    ],
)
def test_plus_strand_categories(rowname, expected):
    gc = GeneClassifier(plus_gene())
    assert gc.classify(parse_rowname(rowname)) == expected


def test_minus_strand_swaps_alt_sites():
    """On the - strand the donor is the high coord, so a novel *high* end is
    alt_5p and a novel *low* end is alt_3p — the opposite of + strand."""
    gc = GeneClassifier(minus_gene())
    # annotated donor(high)=1999? no. On - strand intron (1201,1999): donor=1999,
    # acceptor=1201. Junction with annotated donor 1999 + novel low end:
    assert gc.classify(parse_rowname("chr1:1050-1999:-")) == "alt_3p"
    # annotated acceptor 1201 + novel high end:
    assert gc.classify(parse_rowname("chr1:1201-2600:-")) == "alt_5p"
    # still detects the annotated intron regardless of strand:
    assert gc.classify(parse_rowname("chr1:1201-1999:-")) == "annotated"


def test_tolerance_absorbs_off_by_one():
    gc = GeneClassifier(plus_gene(), tol=2)
    assert gc.classify(parse_rowname("chr1:1202-1998:+")) == "annotated"
