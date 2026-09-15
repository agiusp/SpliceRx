from app.services.junctions import Gene, map_junctions_to_genes

G1 = Gene("TESTG1", "TESTG1", "chr1", 1000, 5200, "+")
G2 = Gene("TESTG2", "TESTG2", "chr1", 10000, 14000, "+")


def test_overlap_strand_chrom_and_multigene():
    rownames = [
        "chr1:1201-1999:+",     # in G1
        "1:2201-2999:+",         # in G1, chr-naming normalised
        "chr1:1201-1999:-",      # wrong strand
        "chr2:1201-1999:+",      # wrong chrom
        "chr1:6000-7000:+",      # between genes, no overlap
        "chr1:10301-10999:+",    # in G2
        "chr1:900-10500:+",      # spans both -> maps to G1 and G2
        "junk",                  # bad row name
    ]
    m = map_junctions_to_genes(rownames, [G1, G2])
    assert m.by_rowname["chr1:1201-1999:+"] == ["TESTG1"]
    assert m.by_rowname["1:2201-2999:+"] == ["TESTG1"]
    assert set(m.by_rowname["chr1:900-10500:+"]) == {"TESTG1", "TESTG2"}
    assert "chr1:1201-1999:-" not in m.by_rowname
    assert "chr1:6000-7000:+" not in m.by_rowname
    assert m.n_bad_rownames == 1
    assert set(m.genes_present()) == {"TESTG1", "TESTG2"}
