from tests.sjv.conftest import FIXTURES

from sjv.services.gencode import annotation_from_gtf


def test_gene_lookup_hit_and_near_match():
    ann = annotation_from_gtf(FIXTURES / "mini.gtf", label="mini")
    g = ann.get_gene("testg1")  # case-insensitive
    assert g is not None
    assert (g.chrom, g.start, g.end, g.strand) == ("chr1", 1000, 5200, "+")

    assert ann.get_gene("NOPE") is None
    assert "TESTG1" in ann.near_matches("test")


def test_gene_suggest_prefix():
    ann = annotation_from_gtf(FIXTURES / "mini.gtf", label="mini")
    assert ann.suggest("test") == ["TESTG1"]
    assert ann.suggest("OT") == ["OTHERG"]
    assert ann.suggest("") == []
    assert ann.suggest("zzz") == []


def test_transcript_models_and_cds_utr_split():
    ann = annotation_from_gtf(FIXTURES / "mini.gtf", label="mini")
    g = ann.get_gene("TESTG1")
    txs = {t.transcript_id: t for t in ann.get_transcripts(g.gene_id)}
    assert set(txs) == {"TESTG1.1", "TESTG1.2", "TESTG1.3"}

    t1 = txs["TESTG1.1"]
    # exon 2 (2000-2200) is split: 2000-2049 UTR, 2050-2200 CDS
    segs = sorted((e.start, e.end, e.kind) for e in t1.exons)
    assert (2000, 2049, "UTR") in segs
    assert (2050, 2200, "CDS") in segs
    # exon 1 is entirely UTR
    assert (1000, 1200, "UTR") in segs
    # exon 3 entirely CDS
    assert (3000, 3200, "CDS") in segs


def test_classifier_from_real_fixture_gtf():
    from sjv.services.classify import GeneClassifier
    from sjv.services.junctions import parse_rowname

    ann = annotation_from_gtf(FIXTURES / "mini.gtf", label="mini")
    g = ann.get_gene("TESTG1")
    gc = GeneClassifier(ann.get_transcripts(g.gene_id))

    cases = {
        "chr1:1201-1999:+": "annotated",
        "chr1:2201-3999:+": "exon_skipping",
        "chr1:1150-1999:+": "alt_5p",
        "chr1:2201-2500:+": "alt_3p",
        "chr1:1201-4999:+": "isoform_switch",
        "chr1:2400-2600:+": "novel_exon",
        "chr1:1500-3500:+": "novel",
    }
    for rn, expected in cases.items():
        assert gc.classify(parse_rowname(rn)) == expected, rn
