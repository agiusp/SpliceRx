import warnings

import pytest

from sjvc.services.sjdat import SJDAT_META

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from main import app
    return TestClient(app)


def test_scan_endpoint(client, cohort_dir):
    r = client.post("/api/dataload/scan", json={"path": str(cohort_dir)})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["cohort"] == "TCGA_TEST"
    roles = {f["name"]: f["role"] for f in body["files"]}
    assert roles["TCGA_TEST_novel_junction_counts_per_gene.rds"] == "gene_matrix"
    assert roles["TCGA_TEST_junction_counts.rds"] == "junction_counts"


def test_scan_rejects_outside_root(client, tmp_path, monkeypatch):
    monkeypatch.setenv("SJ_DATA_ROOT", str(tmp_path / "root"))
    (tmp_path / "root").mkdir()
    r = client.post("/api/dataload/scan", json={"path": str(tmp_path / "elsewhere")})
    assert r.status_code in (403, 404)


def test_load_gene_matrix_and_clinical_into_sjvc(client, cohort_dir):
    sid = client.post("/api/sjvc/session").json()["session_id"]

    gm = cohort_dir / "TCGA_TEST_novel_junction_counts_per_gene.rds"
    r = client.post(f"/api/dataload/sjvc/{sid}/gene-matrix", json={"path": str(gm)})
    assert r.status_code == 200, r.text
    assert r.json()["feature_kind"] == "gene"
    assert len(r.json()["samples"]) == 12

    cl = cohort_dir / "TCGA_TEST_sample_metadata.rds"
    r = client.post(f"/api/dataload/sjvc/{sid}/clinical", json={"path": str(cl)})
    assert r.status_code == 200, r.text
    assert any(c["name"] == "subtype" for c in r.json()["columns"])

    st = client.get(f"/api/sjvc/session/{sid}/state").json()
    assert st["has_junctions"] and st["feature_kind"] == "gene"
    assert st["n_junctions"] == 3
    assert st["clinical"] is not None
    assert {c["name"] for c in st["clinical"]["columns"]} >= {"subtype", "stage"}


def test_clinical_before_matrix_is_rejected(client, cohort_dir):
    sid = client.post("/api/sjvc/session").json()["session_id"]
    cl = cohort_dir / "TCGA_TEST_sample_metadata.rds"
    r = client.post(f"/api/dataload/sjvc/{sid}/clinical", json={"path": str(cl)})
    assert r.status_code == 409


def test_state_of_empty_session(client):
    sid = client.post("/api/sjvc/session").json()["session_id"]
    st = client.get(f"/api/sjvc/session/{sid}/state").json()
    assert st == {
        "has_junctions": False,
        "feature_kind": None,
        "samples": [],
        "n_junctions": 0,
        "clinical": None,
        "gencode_label": None,
        "pathways_enabled": False,
        "sjdat_options": [
            {
                "kind": kind, "label": label, "description": desc,
                "loaded": False, "n_features": 0, "n_samples": 0, "sparse": False,
                "n_nonzero_rows": None,
            }
            for kind, (label, desc) in SJDAT_META.items()
        ],
        "active_sjdat": None,
        "has_junction_metadata": False,
    }


def test_load_junction_matrix_and_metadata_into_sjv(client, cohort_dir):
    sid = client.post("/api/sjv/session").json()["session_id"]

    jc = cohort_dir / "TCGA_TEST_junction_counts.rds"
    r = client.post(f"/api/dataload/sjv/{sid}/junctions", json={"path": str(jc)})
    assert r.status_code == 200, r.text
    assert r.json()["sparse"] is True
    assert r.json()["n_junctions"] == 10
    assert r.json()["samples"] == ["sample_A", "sample_B", "sample_C"]

    md = cohort_dir / "TCGA_TEST_sj_meta.csv"
    r = client.post(f"/api/dataload/sjv/{sid}/sample-metadata", json={"path": str(md)})
    assert r.status_code == 200, r.text
    assert any(c["name"] == "condition" for c in r.json()["columns"])

    st = client.get(f"/api/sjv/session/{sid}/state").json()
    assert st["has_rds"] and st["sparse"] and st["n_junctions"] == 10
    assert "condition" in st["sample_metadata_columns"]


def test_sjv_metadata_before_matrix_is_rejected(client, cohort_dir):
    sid = client.post("/api/sjv/session").json()["session_id"]
    md = cohort_dir / "TCGA_TEST_sj_meta.csv"
    r = client.post(f"/api/dataload/sjv/{sid}/sample-metadata", json={"path": str(md)})
    assert r.status_code == 409


def test_gencode_gtf_applies_to_both_apps(client, cohort_dir):
    sjv = client.post("/api/sjv/session").json()["session_id"]
    sjvc = client.post("/api/sjvc/session").json()["session_id"]
    gtf = cohort_dir / "TCGA_TEST.gtf"

    r = client.post(
        "/api/dataload/gtf",
        json={"sjv_sid": sjv, "sjvc_sid": sjvc, "path": str(gtf)},
    )
    assert r.status_code == 200, r.text
    assert set(r.json()["applied_to"]) == {"sjv", "sjvc"}

    assert client.get(f"/api/sjv/session/{sjv}/state").json()["gencode_label"]
    assert client.get(f"/api/sjvc/session/{sjvc}/state").json()["gencode_label"]


def test_gencode_needs_a_live_session(client, cohort_dir):
    r = client.post(
        "/api/dataload/gtf",
        json={"sjv_sid": "nope", "sjvc_sid": None, "path": str(cohort_dir / "TCGA_TEST.gtf")},
    )
    assert r.status_code == 404


def test_gencode_releases_endpoint(client):
    r = client.get("/api/dataload/gencode/releases")
    assert r.status_code == 200
    assert "human" in r.json()["releases"]


def test_load_junction_counts_and_gene_matrix_into_sjvc(client, cohort_dir):
    """NSJCG can now pick any of the 3 sjdat matrices, like SJSurv."""
    sid = client.post("/api/sjvc/session").json()["session_id"]

    jc = cohort_dir / "TCGA_TEST_junction_counts.rds"
    r = client.post(f"/api/dataload/sjvc/{sid}/sjdat/junction_counts", json={"path": str(jc)})
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "junction_counts"

    gm = cohort_dir / "TCGA_TEST_novel_junction_counts_per_gene.rds"
    r = client.post(f"/api/dataload/sjvc/{sid}/sjdat/gene_matrix", json={"path": str(gm)})
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "gene_matrix"

    st = client.get(f"/api/sjvc/session/{sid}/state").json()
    opts = {o["kind"]: o for o in st["sjdat_options"]}
    assert opts["junction_counts"]["loaded"] and opts["gene_matrix"]["loaded"]
    assert st["active_sjdat"] == "junction_counts"
