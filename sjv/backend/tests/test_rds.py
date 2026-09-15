import shutil

from tests.conftest import FIXTURES

from app.services.rds import load_rds


def test_load_fixture_matrix(tmp_path):
    src = tmp_path / "upload.rds"
    shutil.copy(FIXTURES / "mini.rds", src)

    m = load_rds(src, tmp_path)
    assert m.samples == ["sample_A", "sample_B", "sample_C"]
    # 11 rows in, 1 bad row name skipped -> 10 parsed junctions
    assert m.n_junctions == 10
    assert m.n_bad_rownames == 1

    col = m.column("sample_A")
    by_name = {j.rowname: col[i] for i, j in enumerate(m.junctions)}
    assert by_name["chr1:1201-1999:+"] == 50
    assert by_name["chr1:2400-2600:+"] == 12


def test_end_to_end_plot_shape(tmp_path):
    """rds + gtf -> the pieces routes.plot assembles."""
    from app.services.classify import GeneClassifier
    from app.services.gencode import annotation_from_gtf
    from app.services.junctions import filter_for_gene, scale_counts

    src = tmp_path / "upload.rds"
    shutil.copy(FIXTURES / "mini.rds", src)
    m = load_rds(src, tmp_path)
    ann = annotation_from_gtf(FIXTURES / "mini.gtf", label="mini")
    gene = ann.get_gene("TESTG1")

    in_locus = filter_for_gene(m.junctions, gene)
    idx = {j.rowname: i for i, j in enumerate(m.junctions)}
    counts = m.column("sample_A")
    arcs = scale_counts(in_locus, [counts[idx[j.rowname]] for j in in_locus])

    # off-locus (chr2) and zero-count (chr1:3201-3999) rows are gone
    drawn = {a.junction.rowname for a in arcs}
    assert "chr2:1000-2000:+" not in drawn
    assert "chr1:3201-3999:+" not in drawn
    assert len(drawn) == 8

    gc = GeneClassifier(ann.get_transcripts(gene.gene_id))
    cats = {a.junction.rowname: gc.classify(a.junction) for a in arcs}
    assert cats["chr1:1201-1999:+"] == "annotated"
    assert cats["chr1:1201-4999:+"] == "isoform_switch"
