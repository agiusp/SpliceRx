import warnings

import pytest

from tests.sjvc.conftest import FIXTURES

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from sjvc.main import app
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


def test_matrix_restricted_to_clinical_intersection(client, tmp_path):
    """Only samples present in *both* the matrix and the clinical table survive
    into the session's working matrix — and the two ingest responses agree on
    exactly how many were excluded and why."""
    sid = client.post("/api/session").json()["session_id"]
    mres = client.post(f"/api/session/{sid}/junctions",
                        files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))}).json()
    assert mres["samples"] == [f"s{i:02d}" for i in range(1, 13)]  # all 12, no clinical yet
    assert mres["warnings"] == []

    # a clinical table missing s11/s12 and carrying one extra unknown sample_id
    partial = tmp_path / "partial_clinical.csv"
    lines = ["sample_id,subtype"]
    lines += [f"s{i:02d},LUAD" for i in range(1, 11)]
    lines.append("ghost,LUAD")
    partial.write_text("\n".join(lines) + "\n")

    cres = client.post(f"/api/session/{sid}/clinical", files={"file": ("c.csv", open(partial, "rb"))}).json()
    assert cres["n_matched"] == 10
    assert any("ghost" in w for w in cres["warnings"])          # clinical row not in the matrix
    assert any("2" in w and "sample(s) excluded" in w for w in cres["warnings"])  # s11, s12 dropped

    st = client.get(f"/api/session/{sid}/state").json()
    assert st["samples"] == [f"s{i:02d}" for i in range(1, 11)]
    assert st["n_junctions"] > 0

    # re-uploading the FULL clinical table recovers all 12 — the intersection is
    # always recomputed from the pristine matrix, never narrowed further
    full = client.post(f"/api/session/{sid}/clinical",
                        files={"file": ("c2.csv", open(FIXTURES / "mini_clinical.csv", "rb"))}).json()
    assert full["n_matched"] == 12
    assert full["warnings"] == []
    st2 = client.get(f"/api/session/{sid}/state").json()
    assert st2["samples"] == [f"s{i:02d}" for i in range(1, 13)]


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


def test_mad_features_coverage_filter(client):
    """n_min/x_min are an optional prefilter on top_features_by_mad(), shown
    by the frontend only for a junction-level sjdat — mini_junctions.rds has
    11 parseable junction rows (a 12th, "not_a_junction", never counts). With
    x_min=30, rows 5 ("chr1:2400-3600:+") and 7 ("chr1:11301-11999:+") have
    fewer than 6 samples at >=30 and are dropped; every other row has >=6."""
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})

    # n_min omitted (the default, and what the frontend sends for a
    # gene-level sjdat) — unchanged from before this parameter existed
    unfiltered = client.post(f"/api/session/{sid}/features/mad", json={"top_n": 20}).json()
    assert unfiltered["n_features"] == 11

    filtered = client.post(
        f"/api/session/{sid}/features/mad", json={"top_n": 20, "n_min": 6, "x_min": 30}
    ).json()
    assert filtered["n_features"] == 9

    # a threshold nothing can meet -> a clear error, same wording as SJSurv's
    impossible = client.post(
        f"/api/session/{sid}/features/mad", json={"top_n": 20, "n_min": 100, "x_min": 0}
    )
    assert impossible.status_code == 422
    assert "coverage threshold" in impossible.json()["detail"]


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


def test_features_mad_top_n_refines_resolved_gene_set(client):
    """mad_top_n narrows an already-resolved gene set down to its most
    variable members — OTHERG is constant (all 1's) so it's the one dropped
    when asked for the top 2 of the 3 typed genes."""
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})
    client.post(f"/api/session/{sid}/geneset",
                json={"mode": "typed", "text": "TESTG1, TESTG2, OTHERG"})

    f = client.post(f"/api/session/{sid}/features", json={"condense": False, "mad_top_n": 2})
    assert f.status_code == 200, f.text
    body = f.json()
    assert body["n_features"] == 2
    assert set(body["feature_preview"]) == {"TESTG1", "TESTG2"}
    assert any("ranked" in w and "kept the top 2" in w for w in body["warnings"])

    # mad_top_n omitted, or >= the resolved count, is unchanged — every
    # resolved feature is kept and no warning is added
    unfiltered = client.post(f"/api/session/{sid}/features", json={"condense": False}).json()
    assert unfiltered["n_features"] == 3 and unfiltered["warnings"] == []
    same = client.post(f"/api/session/{sid}/features", json={"condense": False, "mad_top_n": 10}).json()
    assert same["n_features"] == 3 and same["warnings"] == []


def test_features_csv_export(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})

    # 409 before a feature matrix exists
    assert client.get(f"/api/session/{sid}/features.csv").status_code == 409

    # top_n=3 asked, but OTHERG is constant (all 1's) across every sample —
    # only TESTG1 / TESTG2 actually vary and are eligible to be ranked
    mad = client.post(f"/api/session/{sid}/features/mad", json={"top_n": 3}).json()
    assert mad["n_features"] == 2
    assert any("variability" in w for w in mad["warnings"])

    r = client.get(f"/api/session/{sid}/features.csv")
    assert r.status_code == 200
    assert "text/csv" in r.headers["content-type"]
    assert "attachment" in r.headers["content-disposition"]

    lines = r.text.strip().splitlines()
    assert lines[0].split(",")[0] == "feature"
    assert len(lines[0].split(",")) == 1 + 12  # header + 12 samples
    assert len(lines) == 1 + 2  # header + the 2 non-constant gene rows
    assert {ln.split(",")[0] for ln in lines[1:]} == {"TESTG1", "TESTG2"}


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

    # of the 3 genes, OTHERG is constant (all 1's) across every sample, so it
    # never has any variability to rank regardless of the coding filter
    allg = client.post(f"/api/session/{sid}/features/mad",
                       json={"top_n": 10, "protein_coding_only": False}).json()
    assert set(allg["feature_preview"]) == {"TESTG1", "TESTG2"}
    assert any("variability" in w for w in allg["warnings"])

    # protein-coding only also drops OTHERG (gene_type "lncRNA" in the fixture
    # GTF) — same result here, but via the coding filter this time
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


def _session_with_missing_clinical(client, tmp_path):
    """Same cohort as `ready_session`, but s12's subtype is blank — a stand-in
    for a real "this sample has no value for the selected clinical feature"
    case, to exercise the missing-clinical toggle on projection/heatmap."""
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/junctions",
                files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    clinical_csv = tmp_path / "clinical_with_gap.csv"
    lines = ["sample_id,subtype,age"]
    lines += [f"s{i:02d},{'LUAD' if i <= 6 else 'LUSC'},{50 + i}" for i in range(1, 12)]
    lines.append("s12,,62")  # missing subtype
    clinical_csv.write_text("\n".join(lines) + "\n")
    client.post(f"/api/session/{sid}/clinical", files={"file": ("c.csv", open(clinical_csv, "rb"))})
    client.post(f"/api/session/{sid}/gtf", files={"file": ("mini.gtf", open(FIXTURES / "mini.gtf", "rb"))})
    client.post(f"/api/session/{sid}/geneset", json={"mode": "typed", "text": "TESTG1, TESTG2, BOGUS"})
    client.post(f"/api/session/{sid}/features", json={"condense": False})
    return sid


def test_projection_flags_missing_clinical_point(client, tmp_path):
    sid = _session_with_missing_clinical(client, tmp_path)
    p = client.post(f"/api/session/{sid}/projection",
                     json={"method": "pca", "clinical": ["subtype"]}).json()
    assert {pt["sample"] for pt in p["points"] if pt["missing"]} == {"s12"}
    s12 = next(pt for pt in p["points"] if pt["sample"] == "s12")
    assert s12["color"] == "#cfcfcf"
    assert all(pt["color"] != "#cfcfcf" for pt in p["points"] if pt["sample"] != "s12")

    # no clinical feature selected at all -> nothing reads as "missing"
    p0 = client.post(f"/api/session/{sid}/projection", json={"method": "pca", "clinical": []}).json()
    assert not any(pt["missing"] for pt in p0["points"])


def test_heatmap_drop_missing_clinical(client, tmp_path):
    sid = _session_with_missing_clinical(client, tmp_path)

    kept = client.post(f"/api/session/{sid}/heatmap",
                       json={"clinical": ["subtype"], "order": "cluster"}).json()
    assert "s12" in kept["samples"]

    dropped = client.post(f"/api/session/{sid}/heatmap",
                          json={"clinical": ["subtype"], "order": "cluster",
                                "drop_missing_clinical": True}).json()
    assert "s12" not in dropped["samples"]
    assert len(dropped["samples"]) == len(kept["samples"]) - 1
    assert any("dropped 1" in w and "subtype" in w for w in dropped["warnings"])


def test_heatmap_drop_missing_clinical_is_a_noop_with_no_feature_selected(client, ready_session):
    h = client.post(f"/api/session/{ready_session}/heatmap",
                    json={"clinical": [], "order": "cluster", "drop_missing_clinical": True}).json()
    assert len(h["samples"]) == 12
