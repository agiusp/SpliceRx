import warnings

import pytest

from tests.conftest import FIXTURES

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from app.main import app
    return TestClient(app)


@pytest.fixture
def ready_session(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    client.post(f"/api/session/{sid}/clinical",
                files={"file": ("c.csv", open(FIXTURES / "mini_clinical.csv", "rb"))})
    client.post(f"/api/session/{sid}/gtf",
                files={"file": ("mini.gtf", open(FIXTURES / "mini.gtf", "rb"))})
    client.post(f"/api/session/{sid}/geneset", json={"mode": "typed", "text": "TESTG1, TESTG2, BOGUS"})
    client.post(f"/api/session/{sid}/features", json={"condense": False})
    return sid


def test_clinical_eligibility_surface(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    cols = client.post(f"/api/session/{sid}/clinical",
                       files={"file": ("c.csv", open(FIXTURES / "mini_clinical.csv", "rb"))}).json()["columns"]
    by = {c["name"]: c for c in cols}
    assert by["age"]["eligible"] and by["age"]["type"] == "numeric"
    assert not by["notes"]["eligible"]


def test_geneset_reports_unmatched(client, ready_session):
    # done in fixture; re-run to inspect
    r = client.post(f"/api/session/{ready_session}/geneset",
                    json={"mode": "typed", "text": "TESTG1, ZZZ"}).json()
    assert r["matched"] == ["TESTG1"] and r["unmatched"] == ["ZZZ"]


def test_features_condense_toggle(client, ready_session):
    raw = client.post(f"/api/session/{ready_session}/features", json={"condense": False}).json()
    con = client.post(f"/api/session/{ready_session}/features", json={"condense": True}).json()
    assert raw["feature_kind"] == "junction" and raw["n_features"] == 10
    assert con["feature_kind"] == "gene" and con["n_features"] == 2


def test_mad_features_skip_geneset_entirely(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    client.post(f"/api/session/{sid}/gtf",
                files={"file": ("mini.gtf", open(FIXTURES / "mini.gtf", "rb"))})

    # no /geneset call at all
    r = client.post(f"/api/session/{sid}/features/mad", json={"condense": False, "top_n": 4})
    assert r.status_code == 200
    body = r.json()
    assert body["feature_kind"] == "junction" and body["n_features"] == 4

    # feeds straight into projection/heatmap like any other feature build
    p = client.post(f"/api/session/{sid}/projection", json={"method": "pca", "clinical": []})
    assert p.status_code == 200 and len(p.json()["points"]) == 12

    # condensed mode ranks genes genome-wide, including genes never typed anywhere
    rg = client.post(f"/api/session/{sid}/features/mad", json={"condense": True, "top_n": 3})
    assert rg.status_code == 200
    assert set(rg.json()["feature_preview"]) == {"TESTG1", "TESTG2", "OTHERG"}


def test_mad_features_requires_matrix_first(client):
    sid = client.post("/api/session").json()["session_id"]
    r = client.post(f"/api/session/{sid}/features/mad", json={"condense": False, "top_n": 4})
    assert r.status_code == 409


def test_already_condensed_gene_matrix(client):
    sid = client.post("/api/session").json()["session_id"]
    up = client.post(f"/api/session/{sid}/junctions",
                     files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))}).json()
    assert up["feature_kind"] == "gene"
    assert up["warnings"]  # a "treating this as gene-level" note

    # MAD ranks the gene rows directly, no GENCODE needed
    r = client.post(f"/api/session/{sid}/features/mad", json={"condense": False, "top_n": 2})
    assert r.status_code == 200
    body = r.json()
    assert body["feature_kind"] == "gene" and body["n_features"] == 2
    assert set(body["feature_preview"]) <= {"TESTG1", "TESTG2", "OTHERG"}

    # feeds a projection like any other feature build
    p = client.post(f"/api/session/{sid}/projection", json={"method": "pca", "clinical": []})
    assert p.status_code == 200 and len(p.json()["points"]) == 12


def test_gene_level_matrix_subset_by_name(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})

    # typed gene names are matched straight against the row labels — no GENCODE
    r = client.post(f"/api/session/{sid}/geneset",
                    json={"mode": "typed", "text": "TESTG1, OTHERG, NOPE"}).json()
    assert set(r["matched"]) == {"TESTG1", "OTHERG"}
    assert r["unmatched"] == ["NOPE"]

    f = client.post(f"/api/session/{sid}/features", json={"condense": False})
    assert f.status_code == 200
    body = f.json()
    assert body["feature_kind"] == "gene" and body["n_features"] == 2
    assert set(body["feature_preview"]) == {"TESTG1", "OTHERG"}

    # and it plots
    p = client.post(f"/api/session/{sid}/projection", json={"method": "pca", "clinical": []})
    assert p.status_code == 200 and len(p.json()["points"]) == 12

    # nothing matched -> clear error
    bad = client.post(f"/api/session/{sid}/geneset", json={"mode": "typed", "text": "ZZZ, QQQ"})
    assert bad.status_code == 422


def test_features_csv_export(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})

    # 409 before a feature matrix exists
    assert client.get(f"/api/session/{sid}/features.csv").status_code == 409

    client.post(f"/api/session/{sid}/features/mad", json={"top_n": 3})
    r = client.get(f"/api/session/{sid}/features.csv")
    assert r.status_code == 200
    assert "text/csv" in r.headers["content-type"]
    assert "attachment" in r.headers["content-disposition"]

    lines = r.text.strip().splitlines()
    assert lines[0].split(",")[0] == "feature"
    assert len(lines[0].split(",")) == 1 + 12  # header + 12 samples
    assert len(lines) == 1 + 3  # header + 3 gene rows
    assert {ln.split(",")[0] for ln in lines[1:]} == {"TESTG1", "TESTG2", "OTHERG"}


def test_mad_protein_coding_filter(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})

    # without a reference, the filter can't be applied
    r = client.post(f"/api/session/{sid}/features/mad",
                    json={"top_n": 10, "protein_coding_only": True})
    assert r.status_code == 409

    client.post(f"/api/session/{sid}/gtf",
                files={"file": ("mini.gtf", open(FIXTURES / "mini.gtf", "rb"))})

    # all genes: TESTG1, TESTG2, OTHERG
    allg = client.post(f"/api/session/{sid}/features/mad",
                       json={"top_n": 10, "protein_coding_only": False}).json()
    assert set(allg["feature_preview"]) == {"TESTG1", "TESTG2", "OTHERG"}

    # protein-coding only drops OTHERG (gene_type "lncRNA" in the fixture GTF)
    pc = client.post(f"/api/session/{sid}/features/mad",
                     json={"top_n": 10, "protein_coding_only": True}).json()
    assert set(pc["feature_preview"]) == {"TESTG1", "TESTG2"}
    assert any("protein-coding" in w for w in pc["warnings"])


def test_projection_encodings(client, ready_session):
    p1 = client.post(f"/api/session/{ready_session}/projection",
                     json={"method": "pca", "clinical": ["subtype"]}).json()
    assert p1["encoding"]["color"]["kind"] == "categorical"
    assert {pt["clinical"]["subtype"] for pt in p1["points"]} == {"LUAD", "LUSC"}

    p2 = client.post(f"/api/session/{ready_session}/projection",
                     json={"method": "pca", "clinical": ["age", "tmb"]}).json()
    assert p2["encoding"]["color"]["kind"] == "bivariate"

    p3 = client.post(f"/api/session/{ready_session}/projection",
                     json={"method": "pca", "clinical": ["tmb", "subtype"]}).json()
    assert p3["encoding"]["color"]["kind"] == "sequential"
    assert p3["encoding"]["shape"]["kind"] == "categorical"


def test_heatmap_cluster_and_group(client, ready_session):
    h = client.post(f"/api/session/{ready_session}/heatmap",
                    json={"clinical": ["subtype", "age"], "order": "cluster"}).json()
    assert h["col_dendro"] is not None
    assert [a["feature"] for a in h["annotations"]] == ["subtype", "age"]

    hg = client.post(f"/api/session/{ready_session}/heatmap",
                     json={"clinical": ["subtype"], "order": "group", "group_by": "subtype"}).json()
    assert hg["col_dendro"] is None
    assert hg["samples"][:6] == [f"s{i:02d}" for i in range(1, 7)]


def test_heatmap_annotation_legends_alphabetical_and_distinct(client, ready_session):
    h = client.post(f"/api/session/{ready_session}/heatmap",
                    json={"clinical": ["response", "subtype"], "order": "cluster"}).json()
    a = {t["feature"]: t for t in h["annotations"]}

    # legend entries in alphabetical (natural) order
    assert [e["value"] for e in a["response"]["legend"]] == ["CR", "PD", "PR", "SD"]
    assert [e["value"] for e in a["subtype"]["legend"]] == ["LUAD", "LUSC"]

    # the two categorical tracks use non-overlapping colours
    c1 = {e["color"] for e in a["response"]["legend"]}
    c2 = {e["color"] for e in a["subtype"]["legend"]}
    assert c1.isdisjoint(c2)

    # two numeric tracks get different ramps
    hn = client.post(f"/api/session/{ready_session}/heatmap",
                     json={"clinical": ["age", "tmb"], "order": "cluster"}).json()
    an = {t["feature"]: t for t in hn["annotations"]}
    assert an["age"]["stops"] != an["tmb"]["stops"]
