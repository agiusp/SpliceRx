import pytest

from sjv.services.junctions import (
    Gene,
    RownameError,
    filter_for_gene,
    parse_rowname,
    scale_counts,
)


def test_parse_rowname_basic():
    j = parse_rowname("chr1:100-200:+")
    assert (j.chrom, j.start, j.end, j.strand) == ("chr1", 100, 200, "+")


def test_parse_rowname_normalises_order():
    j = parse_rowname("chrX:500-100:-")
    assert (j.start, j.end, j.strand) == (100, 500, "-")


@pytest.mark.parametrize("bad", ["not-a-junction", "chr1:100:+", "chr1:100-200", "chr1:100-200:*"])
def test_parse_rowname_rejects(bad):
    with pytest.raises(RownameError):
        parse_rowname(bad)


def test_filter_for_gene_overlap_strand_chrom():
    gene = Gene("G", "G.1", "chr1", 1000, 5000, "+")
    js = [
        parse_rowname("chr1:1200-1999:+"),   # inside
        parse_rowname("chr1:900-1100:+"),     # overlaps left edge
        parse_rowname("chr1:1200-1999:-"),    # wrong strand
        parse_rowname("chr2:1200-1999:+"),    # wrong chrom
        parse_rowname("chr1:6000-7000:+"),    # downstream, no overlap
        parse_rowname("1:1200-1300:+"),        # chr-naming normalised -> kept
    ]
    kept = {j.rowname for j in filter_for_gene(js, gene)}
    assert kept == {"chr1:1200-1999:+", "chr1:900-1100:+", "1:1200-1300:+"}


def test_scale_counts_drops_zero_and_nan():
    js = [parse_rowname(f"chr1:{100+i}-{500+i}:+") for i in range(3)]
    arcs = scale_counts(js, [10.0, 0.0, float("nan")])
    assert [a.junction.rowname for a in arcs] == ["chr1:100-500:+"]


def test_scale_counts_min_reads_cutoff():
    js = [parse_rowname(f"chr1:{100+i}-{500+i}:+") for i in range(4)]
    counts = [3.0, 5.0, 12.0, 40.0]
    assert len(scale_counts(js, counts, min_reads=0)) == 4
    kept = {a.junction.rowname for a in scale_counts(js, counts, min_reads=5)}
    assert kept == {"chr1:101-501:+", "chr1:102-502:+", "chr1:103-503:+"}  # 5, 12, 40
    assert scale_counts(js, counts, min_reads=100) == []


def test_scale_counts_monotonic_height_band():
    js = [parse_rowname(f"chr1:{100+i}-{500+i}:+") for i in range(4)]
    arcs = scale_counts(js, [2.0, 8.0, 32.0, 128.0])
    heights = [a.height for a in arcs]
    assert heights == sorted(heights)
    for a in arcs:
        assert 0.1 <= a.height <= 1.1 + 1e-9


def test_scale_counts_single_junction_constant_height():
    arcs = scale_counts([parse_rowname("chr1:1-9:+")], [42.0])
    assert arcs[0].height == 0.6


def test_scale_counts_all_equal():
    js = [parse_rowname(f"chr1:{100+i}-{500+i}:+") for i in range(3)]
    arcs = scale_counts(js, [7.0, 7.0, 7.0])
    assert all(a.height == 0.6 for a in arcs)
