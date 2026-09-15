import pytest
from fastapi import HTTPException

from dataload.paths import safe_path
from dataload.scan import classify, classify_one


def _roles(files):
    return {f.name: (f.role, f.status, f.target) for f in files}


def test_classify_cohort(cohort_dir):
    r = _roles(classify(cohort_dir))
    assert r["TCGA_TEST_novel_junction_counts_per_gene.rds"] == ("gene_matrix", "ok", "sjvc")
    assert r["TCGA_TEST_sample_metadata.rds"] == ("clinical", "ok", "sjvc")
    assert r["TCGA_TEST_MSI.rds"] == ("clinical_extra", "caution", "sjvc")
    assert r["TCGA_TEST_junction_metadata.rds"] == ("junction_metadata", "ok", "sjvc")
    assert r["TCGA_TEST_junction_counts.rds"] == ("junction_counts", "ok", "sjv")
    assert r["notes.txt"][0] == "unknown"


def test_size_guard_routes_big_rds_to_sjv(tmp_path):
    big = tmp_path / "TCGA_X_weird_name.rds"
    big.write_bytes(b"\x00" * (260 * 1024 * 1024))
    f = classify_one(big)
    assert f.role == "junction_counts" and f.target == "sjv" and f.status == "ok"


def test_gene_matrix_wins_over_junction_counts_substring(tmp_path):
    # the real gene file name contains "junction_counts" AND "per_gene"
    p = tmp_path / "TCGA_X_novel_junction_counts_per_gene.rds"
    p.write_bytes(b"x")
    assert classify_one(p).role == "gene_matrix"


def test_pathway_matrix_wins_over_junction_counts_substring(tmp_path):
    # the real pathway file name contains "junction_counts" AND "per_pathway"
    p = tmp_path / "TCGA_X_novel_junction_counts_per_pathway.rds"
    p.write_bytes(b"x")
    f = classify_one(p)
    assert f.role == "pathway_matrix" and f.target == "sjvc" and f.status == "ok"


def test_gtf_recognised(tmp_path):
    for n in ("gencode.v29.annotation.gtf", "ref.gtf.gz"):
        p = tmp_path / n
        p.write_bytes(b"x")
        f = classify_one(p)
        assert f.role == "gtf" and f.status == "ok"


def test_safe_path_rejects_escape(tmp_path, monkeypatch):
    root = tmp_path / "root"
    (root / "sub").mkdir(parents=True)
    monkeypatch.setenv("SJ_DATA_ROOT", str(root))

    assert safe_path(str(root / "sub")) == (root / "sub").resolve()
    with pytest.raises(HTTPException) as e:
        safe_path(str(tmp_path / "outside"))
    assert e.value.status_code == 403
    with pytest.raises(HTTPException):
        safe_path(str(root / "does-not-exist"))
    with pytest.raises(HTTPException):
        safe_path("")
