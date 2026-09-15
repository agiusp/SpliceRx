import warnings

import pytest

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from main import app
    return TestClient(app)


def test_scan_tags_rrs_and_sjsurv_roles(client, sjsurv_cohort):
    r = client.post("/api/dataload/scan", json={"path": str(sjsurv_cohort)})
    assert r.status_code == 200, r.text
    roles = {f["name"]: f["role"] for f in r.json()["files"]}
    assert roles["TCGA_SURV_RRS_scores.rds"] == "rrs_scores"
    assert roles["TCGA_SURV_novel_junction_counts_per_gene.rds"] == "gene_matrix"


def test_load_cohort_into_sjsurv_and_run(client, sjsurv_cohort):
    sid = client.post("/api/sjsurv/session").json()["session_id"]

    gm = sjsurv_cohort / "TCGA_SURV_novel_junction_counts_per_gene.rds"
    r = client.post(f"/api/dataload/sjsurv/{sid}/sjdat/gene_matrix", json={"path": str(gm)})
    assert r.status_code == 200, r.text
    assert r.json()["n_features"] == 60

    rrs = sjsurv_cohort / "TCGA_SURV_RRS_scores.rds"
    r = client.post(f"/api/dataload/sjsurv/{sid}/sjdat/rrs_scores", json={"path": str(rrs)})
    assert r.status_code == 200, r.text

    md = sjsurv_cohort / "TCGA_SURV_sample_metadata.csv"
    r = client.post(f"/api/dataload/sjsurv/{sid}/metadata", json={"path": str(md)})
    assert r.status_code == 200, r.text
    assert r.json()["n_matched"] == 80
    # auto-stratified with the classic age bands, exactly as the R package used to

    st = client.get(f"/api/sjsurv/session/{sid}/state").json()
    loaded = {o["kind"]: o["loaded"] for o in st["sjdat_options"]}
    assert loaded == {
        "junction_counts": False, "rrs_scores": True, "gene_matrix": True, "pathway_matrix": False,
    }

    # pick RRS, then the gene matrix
    client.post(f"/api/sjsurv/session/{sid}/sjdat", json={"kind": "gene_matrix"})

    # this small cohort's classic age bands fragment every Group below the
    # default min_group_n — coarsen to one age band so samples get a label
    suggest = client.post(f"/api/sjsurv/session/{sid}/age-bands/suggest", json={"n_bands": 1}).json()
    r = client.post(f"/api/sjsurv/session/{sid}/stratify", json={
        "edges": suggest["edges"], "include_lowest": suggest["include_lowest"], "min_group_n": 10,
    })
    assert r.status_code == 200, r.text

    r = client.post(f"/api/sjsurv/session/{sid}/select",
                    json={"group": "__all__", "n_min": 0.25, "x_min": 0, "top_n": 20, "n_cv": 4})
    assert r.status_code == 200, r.text
    assert r.json()["n_selected"] == 20

    r = client.post(f"/api/sjsurv/session/{sid}/cross-validate", json={"n_cv": 4})
    assert r.status_code == 200, r.text
    assert r.json()["n_splits"] == 4

    r = client.post(f"/api/sjsurv/session/{sid}/model")
    assert r.status_code == 200, r.text
    assert len(r.json()["features"]) == 20
