"""SJSurv's Group/SurviverGroup -> Sashimi-plot / NSJCG (see dataload/push_groups.py)."""
import warnings

import pytest

from dataload.push_groups import push_into_clinical, push_into_sample_metadata
from sjv.services.groups import SampleMetadata
from sjvc.services.clinical import Clinical, Column

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from main import app
    return TestClient(app)


# --------------------------------------------------------------------------- #
# service level: the two per-target merge functions
# --------------------------------------------------------------------------- #
def test_push_into_sample_metadata_overwrites_and_adds_column():
    meta = SampleMetadata(
        rows={"s01": {"Group": "STALE", "condition": "tumor"}, "s02": {"condition": "normal"}},
        columns=["condition", "Group"],
    )
    group_map = {"s01": ("A | Early | (30-50]", "Good"), "s02": (None, None)}
    n = push_into_sample_metadata(meta, group_map)
    assert n == 1                                      # only s01 got a real value
    assert meta.rows["s01"]["Group"] == "A | Early | (30-50]"   # overwritten, not "STALE"
    assert meta.rows["s01"]["SurviverGroup"] == "Good"
    assert meta.rows["s02"]["Group"] == ""              # no group -> cleared, not left stale
    assert meta.columns.count("Group") == 1             # not duplicated
    assert "SurviverGroup" in meta.columns


def test_push_into_clinical_recomputes_column_typing():
    clinical = Clinical(
        rows={
            "s01": {"Group": "STALE", "age": "50"},
            "s02": {"age": "60"},
            "s03": {"age": "70"},
        },
        columns=[
            Column(name="age", type="numeric", n_unique=3, n_missing=0, eligible=True),
            Column(name="Group", type="categorical", n_unique=1, n_missing=2, eligible=True),
        ],
        overrides={"Group": "numeric"},   # a stale override that must not survive the push
    )
    group_map = {
        "s01": ("A | Early | (30-50]", "Good"),
        "s02": ("A | Early | (30-50]", "Poor"),
        "s03": (None, None),
    }
    n = push_into_clinical(clinical, group_map)
    assert n == 2
    assert clinical.rows["s01"]["Group"] == "A | Early | (30-50]"
    assert clinical.rows["s03"]["SurviverGroup"] == ""

    cols = {c.name: c for c in clinical.columns}
    assert len(clinical.columns) == 3                   # age, Group, SurviverGroup — no dupes
    assert cols["Group"].n_unique == 1                   # one real value, shared by s01/s02
    assert cols["Group"].n_missing == 1                   # s03
    assert cols["SurviverGroup"].n_unique == 2            # Good, Poor
    assert "Group" not in clinical.overrides              # stale override cleared


# --------------------------------------------------------------------------- #
# route level: the /api/dataload/push-groups contract
# --------------------------------------------------------------------------- #
def _stratified_sjsurv_session(client, sjsurv_cohort):
    sid = client.post("/api/sjsurv/session").json()["session_id"]
    gm = sjsurv_cohort / "TCGA_SURV_novel_junction_counts_per_gene.rds"
    client.post(f"/api/dataload/sjsurv/{sid}/sjdat/gene_matrix", json={"path": str(gm)})
    md = sjsurv_cohort / "TCGA_SURV_sample_metadata.csv"
    client.post(f"/api/dataload/sjsurv/{sid}/metadata", json={"path": str(md)})
    # coarsen so this small cohort actually gets Good/Poor labels (see test_pipeline.py)
    suggest = client.post(f"/api/sjsurv/session/{sid}/age-bands/suggest", json={"n_bands": 1}).json()
    client.post(f"/api/sjsurv/session/{sid}/stratify", json={
        "edges": suggest["edges"], "include_lowest": suggest["include_lowest"], "min_group_n": 10,
    })
    return sid


def test_push_groups_needs_stratified_sjsurv_metadata(client):
    sid = client.post("/api/sjsurv/session").json()["session_id"]
    r = client.post("/api/dataload/push-groups", json={"sjsurv_sid": sid})
    assert r.status_code == 409, r.text


def test_push_groups_warns_when_target_has_no_metadata_yet(client, sjsurv_cohort):
    sjsurv_sid = _stratified_sjsurv_session(client, sjsurv_cohort)
    sjv_sid = client.post("/api/sjv/session").json()["session_id"]   # nothing loaded into it

    r = client.post("/api/dataload/push-groups",
                     json={"sjsurv_sid": sjsurv_sid, "sjv_sid": sjv_sid})
    assert r.status_code == 404, r.text          # sjv has no metadata, no other target given


def test_push_groups_updates_nsjcg_session_in_place(client, sjsurv_cohort):
    """Inject a clinical table directly into an NSJCG session store (same sample
    ids SJSurv's fixture uses), then push and check it landed."""
    from sjvc.services.sessions import store as sjvc_store

    sjsurv_sid = _stratified_sjsurv_session(client, sjsurv_cohort)
    sjvc_sid = client.post("/api/sjvc/session").json()["session_id"]
    s = sjvc_store.get(sjvc_sid)
    s.clinical = Clinical(
        rows={f"S{i:02d}": {"batch": "1"} for i in range(1, 81)},
        columns=[Column(name="batch", type="categorical", n_unique=1, n_missing=0, eligible=False)],
    )

    r = client.post("/api/dataload/push-groups",
                     json={"sjsurv_sid": sjsurv_sid, "sjvc_sid": sjvc_sid})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["applied_to"] == ["sjvc"]
    assert body["n_labelled"] == body["matched"]["sjvc"] > 0

    cols = {c["name"] for c in [c.__dict__ for c in s.clinical.columns]}
    assert {"Group", "SurviverGroup"} <= cols
    assert any(v.get("SurviverGroup") in ("Good", "Poor") for v in s.clinical.rows.values())
