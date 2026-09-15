from sjvc.services import pathways
from sjvc.services.geneset import parse_symbols, resolve
from sjvc.services.junctions import Gene

GENES = {
    "TESTG1": Gene("TESTG1", "ENSG_TESTG1", "chr1", 1000, 5200, "+"),
    "TESTG2": Gene("TESTG2", "ENSG_TESTG2", "chr1", 10000, 14000, "+"),
}


class FakeAnnotation:
    def get_gene(self, name):
        return GENES.get(name.strip().upper())

    def genes_by_prefix(self, prefix, limit=1000):
        p = prefix.strip().upper()
        return [g for k, g in GENES.items() if k.startswith(p)]


def test_parse_symbols_splits_and_dedupes():
    assert parse_symbols("TP53, kras\nMYC; TP53") == ["TP53", "kras", "MYC"]


def test_resolve_matched_and_unmatched():
    gs = resolve(["TESTG1", "testg2", "NOPE"], FakeAnnotation(), source="typed")
    assert [g.name for g in gs.matched] == ["TESTG1", "TESTG2"]
    assert gs.unmatched == ["NOPE"]


def test_resolve_prefix_expands():
    gs = resolve(["TESTG"], FakeAnnotation(), source="typed", prefix=True)
    assert {g.name for g in gs.matched} == {"TESTG1", "TESTG2"}
    gs2 = resolve(["ZZZ"], FakeAnnotation(), source="typed", prefix=True)
    assert gs2.matched == [] and gs2.unmatched == ["ZZZ"]


def test_gmt_parse():
    text = "TERM_A\tdesc\tG1\tG2\tG3\nTERM_B\t\tG2,1.0\tG4,0.5\n"
    terms = pathways._parse_gmt(text)
    assert terms["TERM_A"] == ["G1", "G2", "G3"]
    assert terms["TERM_B"] == ["G2", "G4"]      # weights stripped


def test_bundled_hallmark_loads_and_searches():
    hits = pathways.search("MSigDB Hallmark", "apoptosis")
    assert any("Apoptosis" in h["term"] for h in hits)
    genes = pathways.genes("MSigDB Hallmark", hits[0]["term"])
    assert len(genes) > 10
