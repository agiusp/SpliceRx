import shutil

import pytest

from tests.sjv.conftest import FIXTURES

from sjv.services.rds import load_rds


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("SJV_CACHE_DIR", str(tmp_path / "cache"))


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


def test_load_sparse_matrix_matches_dense(tmp_path):
    """A dgCMatrix RDS reads via the sparse path and gives the same numbers."""
    import numpy as np
    from sjv.services.junctions import Gene

    dense = load_rds(_copy(tmp_path, "mini.rds", "d.rds"), tmp_path)
    sparse = load_rds(_copy(tmp_path, "mini_sparse.rds", "s.rds"), tmp_path)

    assert sparse.sparse is True and dense.sparse is False
    assert sparse.samples == dense.samples
    assert sparse.n_junctions == dense.n_junctions == 10
    assert sparse.n_bad_rownames == 1

    gene = Gene(name="TESTG1", gene_id="g1", chrom="chr1", start=1000, end=5200, strand="+")
    d_rows, s_rows = dense.locus_row_indices(gene), sparse.locus_row_indices(gene)
    assert [j.rowname for j in dense.junctions_at(d_rows)] == [
        j.rowname for j in sparse.junctions_at(s_rows)
    ]
    np.testing.assert_array_equal(
        dense.counts_at(d_rows, "sample_A"), sparse.counts_at(s_rows, "sample_A")
    )
    np.testing.assert_array_equal(
        dense.group_median_at(d_rows, dense.samples),
        sparse.group_median_at(s_rows, sparse.samples),
    )


def test_sparse_load_is_cached(tmp_path, monkeypatch):
    from sjv.services import rds as rds_mod

    src = _copy(tmp_path, "mini_sparse.rds", "s.rds")
    load_rds(src, tmp_path)
    assert any(rds_mod._cache_dir().glob("*/meta.json"))

    # a second load of the same bytes must not shell out to R again
    def _boom(*a, **k):
        raise AssertionError("sparse reader ran despite a warm cache")

    monkeypatch.setattr(rds_mod, "_read_sparse_csc", _boom)
    m = load_rds(_copy(tmp_path, "mini_sparse.rds", "s2.rds"), tmp_path)
    assert m.sparse and m.n_junctions == 10


def _copy(tmp_path, name, as_):
    dst = tmp_path / as_
    shutil.copy(FIXTURES / name, dst)
    return dst


def test_end_to_end_plot_shape(tmp_path):
    """rds + gtf -> the pieces routes.plot assembles."""
    from sjv.services.classify import GeneClassifier
    from sjv.services.gencode import annotation_from_gtf
    from sjv.services.junctions import filter_for_gene, scale_counts

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
