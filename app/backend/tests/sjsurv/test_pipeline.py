import math

import numpy as np
import pytest
from fastapi.testclient import TestClient

from sjsurv.main import app
from sjsurv.services.metadata import (
    ALL_GROUPS,
    CLASSIC_AGE_BANDS,
    DEFAULT_MIN_GROUP_N,
    NO_GROUP,
    AgeBands,
    MetadataError,
    StratifyError,
    bin_ages,
    parse_raw_metadata,
    quantile_age_bands,
    stage_label,
    stratify,
)
from sjsurv.services.model import ModelError, cross_validate, train_full
from sjsurv.services.select import Selection, SelectError, select_features
from sjsurv.services.sjdat import Sjdat

from .conftest import FIXTURES

SAMPLES = [f"S{i:02d}" for i in range(1, 81)]


@pytest.fixture
def client():
    return TestClient(app)


@pytest.fixture
def sid(client):
    return client.post("/api/session").json()["session_id"]


# --------------------------------------------------------------------------- #
# service-level: selection + model on synthetic data
# --------------------------------------------------------------------------- #
def _synthetic_sjdat(n_feat=80, n_samp=60, seed=0):
    rng = np.random.default_rng(seed)
    vals = rng.poisson(3, size=(n_feat, n_samp)).astype(float)
    labels = {f"s{i:02d}": ("Good" if i % 2 else "Poor") for i in range(n_samp)}
    poor = np.array([v == "Poor" for v in labels.values()])
    vals[:15, poor] += rng.poisson(8, size=(15, poor.sum()))   # planted signal
    vals[-5:, :] = 0                                            # dead rows
    d = Sjdat(kind="gene_matrix", features=[f"G{i}" for i in range(n_feat)],
              samples=list(labels), values=vals, sparse=False)
    return d, labels


def test_select_filters_and_ranks():
    d, labels = _synthetic_sjdat()
    sel = select_features(d, list(labels), n_min=5, x_min=1, top_n=20)
    assert sel.n_features == 20
    assert sel.n_after_coverage <= d.n_features - 5      # the 5 all-zero rows dropped
    assert len(sel.sample_ids) == 60


def test_select_percent_threshold_and_errors():
    d, labels = _synthetic_sjdat()
    sel = select_features(d, list(labels), n_min=0.5, x_min=0, top_n=10)
    assert sel.n_features == 10
    with pytest.raises(SelectError):
        select_features(d, list(labels), n_min=999, x_min=0, top_n=10)


def test_cross_validate_beats_chance_on_planted_signal():
    d, labels = _synthetic_sjdat()
    sel = select_features(d, list(labels), n_min=3, x_min=0, top_n=25)
    r = cross_validate(sel, labels, n_splits=5)
    assert r.n_splits == 5 and r.n_samples == 60
    assert 0.0 <= r.auc <= 1.0
    assert r.auc > 0.6                       # planted signal is strong
    assert len(r.messages) >= 4
    assert r.confusion and len(r.confusion) == 2


def test_cross_validate_rejects_too_many_folds():
    d, labels = _synthetic_sjdat(n_samp=12)
    sel = select_features(d, list(labels), n_min=2, x_min=0, top_n=10)
    with pytest.raises(ModelError):
        cross_validate(sel, labels, n_splits=10)


def test_train_full_ranks_features():
    d, labels = _synthetic_sjdat()
    sel = select_features(d, list(labels), n_min=3, x_min=0, top_n=25)
    m = train_full(sel, labels, sjdat_kind="gene_matrix", group=ALL_GROUPS)
    assert len(m.features) == sel.n_features
    w = [f.abs_weight for f in m.features]
    assert w == sorted(w, reverse=True)                 # ranked by |weight|
    assert m.auc_resub >= 0.5


# --------------------------------------------------------------------------- #
# Stage / age-band helpers
# --------------------------------------------------------------------------- #
def test_stage_label():
    assert stage_label("Stage I") == "Early"
    assert stage_label("Stage IIB") == "Early"
    assert stage_label("Stage IIIA") == "Late"
    assert stage_label("Stage IV") == "Late"
    assert stage_label("Stage X") is None
    assert stage_label("Not Reported") is None
    assert stage_label(None) is None


def test_classic_age_bands_leave_low_ages_unbanded():
    labels = bin_ages([25, 35, 55, 75, None], CLASSIC_AGE_BANDS)
    assert labels == [None, "(30-50]", "(50-70]", "(>70)", None]


def test_quantile_age_bands_cover_the_whole_range():
    ages = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
    bands = quantile_age_bands(ages, 2)
    assert bands.n_bands == 2
    assert bands.include_lowest is True
    binned = bin_ages(ages, bands)
    assert all(b is not None for b in binned)            # nothing left unbanded
    assert bands.edges[0] == min(ages) and bands.edges[-1] == max(ages)


def test_quantile_age_bands_do_not_strand_the_max_age():
    # a max that rounds *down* at 1 decimal (88.9117 -> 88.9) must not exclude
    # the sample(s) that actually sit above the rounded edge
    ages = [35.93977, 50.1, 60.2, 70.3, 80.4, 88.9117]
    bands = quantile_age_bands(ages, 2)
    binned = bin_ages(ages, bands)
    assert all(b is not None for b in binned), binned
    assert bands.edges[0] <= min(ages) and bands.edges[-1] >= max(ages)


def test_age_bands_reject_non_increasing_edges():
    with pytest.raises(StratifyError):
        AgeBands(edges=[50, 30, 70])


def test_quantile_age_bands_needs_enough_samples():
    with pytest.raises(StratifyError):
        quantile_age_bands([30, 40], 3)


# --------------------------------------------------------------------------- #
# raw metadata + stratify — service level, against the fixture cohort
# --------------------------------------------------------------------------- #
def test_raw_columns_win_over_stale_derived_columns(tmp_path):
    """A cohort downloaded before surv_cohort() was removed from the R package
    can still carry its old Histology/Stage/Age_at_diagnosis/Group/
    SurviverGroup output alongside the raw columns. Those legacy names must
    never be picked over the real raw ones — they're already binned/derived
    (e.g. Age_at_diagnosis = "(50-70]", not a number) and reading them instead
    silently produces nonsense."""
    csv = tmp_path / "stale.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days,"
        "Histology,Stage,Age_at_diagnosis,Group,SurviverGroup\n"
        "S01,Adenocarcinoma,Stage IIB,55.3,900,STALE,Late,(50-70],STALE | Late | (50-70],Poor\n"
        "S02,Adenocarcinoma,Stage IB,40.1,1800,STALE,Early,(30-50],STALE | Early | (30-50],Good\n"
    )
    raw = parse_raw_metadata(csv, tmp_path, ["S01", "S02"])
    assert raw.rows["S01"]["histology"] == "Adenocarcinoma"     # not "STALE"
    assert raw.rows["S01"]["stage_raw"] == "Stage IIB"          # not "Late"
    assert raw.rows["S01"]["age"] == 55.3                       # not None (from "(50-70]")
    assert stage_label(raw.rows["S01"]["stage_raw"]) == "Early"  # re-derived, not read verbatim (IIB, not "Late")


def test_age_column_picks_the_most_complete_alias(tmp_path):
    """Not every TCGA cohort download has the same age field populated — pick
    whichever candidate actually has the most usable values, not just the
    first one present."""
    csv = tmp_path / "age_variants.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "survival_days,age_at_diagnosis_years,tcga.cgc_case_age_at_diagnosis\n"
        "S01,Adenocarcinoma,Stage I,900,,70\n"
        "S02,Adenocarcinoma,Stage I,1800,,58\n"
        "S03,Adenocarcinoma,Stage I,300,,65\n"
    )
    # age_at_diagnosis_years is present but entirely empty for this cohort —
    # tcga.cgc_case_age_at_diagnosis (fully populated) must win instead
    raw = parse_raw_metadata(csv, tmp_path, ["S01", "S02", "S03"])
    assert raw.rows["S01"]["age"] == 70.0
    assert raw.rows["S02"]["age"] == 58.0


def test_survival_days_built_from_days_to_death_and_last_follow_up(tmp_path):
    """No pre-computed survival_days column at all — fall back to the
    standard TCGA construction: days-to-death for the deceased, days-to-
    last-follow-up (censored) for everyone else."""
    csv = tmp_path / "raw_survival.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,tcga.cgc_case_days_to_death,tcga.cgc_case_days_to_last_follow_up\n"
        "S01,Adenocarcinoma,Stage I,55,900,\n"      # deceased -> days_to_death
        "S02,Adenocarcinoma,Stage I,60,,1800\n"     # censored -> days_to_last_follow_up
        "S03,Adenocarcinoma,Stage I,65,,\n"          # neither -> no survival value
    )
    raw = parse_raw_metadata(csv, tmp_path, ["S01", "S02", "S03"])
    assert raw.rows["S01"]["os"] == 900.0
    assert raw.rows["S02"]["os"] == 1800.0
    assert raw.rows["S03"]["os"] is None


def test_survival_fallback_ignores_unrelated_days_to_columns(tmp_path):
    """days_to_birth / days_to_collection etc must never be mistaken for
    days-to-death or days-to-last-follow-up."""
    csv = tmp_path / "unrelated_days.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,tcga.gdc_cases.diagnoses.days_to_birth,"
        "tcga.cgc_sample_days_to_collection\n"
        "S01,Adenocarcinoma,Stage I,55,-20000,10\n"
    )
    with pytest.raises(MetadataError, match="overall survival"):
        parse_raw_metadata(csv, tmp_path, ["S01"])


def test_pre_computed_survival_days_wins_over_raw_fallback(tmp_path):
    csv = tmp_path / "both.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days,tcga.cgc_case_days_to_death\n"
        "S01,Adenocarcinoma,Stage I,55,1234,999\n"
    )
    raw = parse_raw_metadata(csv, tmp_path, ["S01"])
    assert raw.rows["S01"]["os"] == 1234.0


def test_rownames_column_used_as_sample_id_when_no_sample_id_column(tmp_path):
    """recount3/TCGA sample-metadata tables commonly carry the sample UUID
    only as the data.frame's own row names — pyreadr surfaces those as a
    column literally named 'rownames' once promoted (see sjvc.services.rds.
    read_table); this must be usable as the sample_id, not rejected as
    'no sample_id column found'."""
    csv = tmp_path / "rownames_only.csv"
    csv.write_text(
        "rownames,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days\n"
        "S01,Adenocarcinoma,Stage IIB,55.3,900\n"
        "S02,Adenocarcinoma,Stage IB,40.1,1800\n"
    )
    raw = parse_raw_metadata(csv, tmp_path, ["S01", "S02"])
    assert raw.n_rows == 2
    assert raw.rows["S01"]["histology"] == "Adenocarcinoma"


def test_real_sample_id_column_wins_over_rownames(tmp_path):
    csv = tmp_path / "both.csv"
    csv.write_text(
        "rownames,sample_id,tcga.cgc_case_histological_diagnosis,"
        "tcga.cgc_case_pathologic_stage,age_at_diagnosis_years,survival_days\n"
        "bogus,S01,Adenocarcinoma,Stage IIB,55.3,900\n"
    )
    raw = parse_raw_metadata(csv, tmp_path, ["S01"])
    assert "S01" in raw.rows and "bogus" not in raw.rows


def test_parse_raw_metadata_requires_all_four_columns(tmp_path):
    with pytest.raises(MetadataError):
        parse_raw_metadata(FIXTURES / "mini_metadata_missing_col.csv", tmp_path, SAMPLES)


def test_classic_bands_can_fragment_a_small_cohort_below_min_group_n(tmp_path):
    """Reproduces the reported PAAD problem: with the classic 3 age bands x 2
    stages, every group in this 80-sample fixture falls under the default
    min_group_n=20, so nobody gets a SurviverGroup label."""
    raw = parse_raw_metadata(FIXTURES / "mini_metadata.csv", tmp_path, SAMPLES)
    md = stratify(raw, CLASSIC_AGE_BANDS, DEFAULT_MIN_GROUP_N)
    gcs = md.group_counts()
    assert gcs[0].group == ALL_GROUPS
    assert gcs[0].n_labelled == 0                 # nobody cleared min_group_n
    assert gcs[-1].group == NO_GROUP               # "(no Group)" pinned to the bottom
    assert gcs[-1].n_total == 2                    # the 2 samples with age <= 30
    # every real group (excluding __all__ and "(no Group)") is smaller than 20
    assert all(gc.n_total < 20 for gc in gcs[1:-1])


def _heterogeneous_histology_csv(tmp_path):
    """4 histology values, one sample with none at all — a heterogeneous
    cancer with several small subtypes plus a missing value."""
    csv = tmp_path / "het_hist.csv"
    rows = ["sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
            "age_at_diagnosis_years,survival_days"]
    hist = ["Adenocarcinoma Intestinal Type", "Adenocarcinoma Diffuse Type",
            "Adenocarcinoma Mixed Type", "Adenocarcinoma NOS"]
    for i in range(1, 21):
        h = hist[i % 4] if i != 1 else ""   # S01 has no histology at all
        rows.append(f"S{i:02d},{h},Stage I,55,{900 + i}")
    csv.write_text("\n".join(rows) + "\n")
    return csv


def test_use_histology_false_drops_histology_from_group(tmp_path):
    csv = _heterogeneous_histology_csv(tmp_path)
    samples = [f"S{i:02d}" for i in range(1, 21)]
    raw = parse_raw_metadata(csv, tmp_path, samples)
    bands = quantile_age_bands(raw.ages(), 1)

    with_hist = stratify(raw, bands, min_group_n=1, use_histology=True)
    # S01 has no histology -> no Group at all when histology is required
    assert with_hist.rows["S01"]["group"] is None
    # 4 distinct histology values x 1 stage x 1 age band -> (up to) 4 groups
    groups_with = {r["group"] for r in with_hist.rows.values() if r["group"] is not None}
    assert len(groups_with) == 4

    without_hist = stratify(raw, bands, min_group_n=1, use_histology=False)
    # S01 now gets a Group too — histology is no longer required
    assert without_hist.rows["S01"]["group"] is not None
    # everyone collapses into a single Stage | age-band group
    groups_without = {r["group"] for r in without_hist.rows.values() if r["group"] is not None}
    assert len(groups_without) == 1
    assert "|" in groups_without.pop() and "Adenocarcinoma" not in without_hist.rows["S02"]["group"]


def test_histology_map_merges_selected_values(tmp_path):
    csv = _heterogeneous_histology_csv(tmp_path)
    samples = [f"S{i:02d}" for i in range(1, 21)]
    raw = parse_raw_metadata(csv, tmp_path, samples)
    bands = quantile_age_bands(raw.ages(), 1)

    hmap = {
        "Adenocarcinoma Intestinal Type": "Adenocarcinoma (merged)",
        "Adenocarcinoma Diffuse Type": "Adenocarcinoma (merged)",
    }
    merged = stratify(raw, bands, min_group_n=1, use_histology=True, histology_map=hmap)
    groups = {r["group"] for r in merged.rows.values() if r["group"] is not None}
    # the two merged values collapse into one group; the two untouched ones stay separate
    assert len(groups) == 3
    assert any("Adenocarcinoma (merged)" in g for g in groups)
    # unmapped values pass through unchanged
    assert any("Adenocarcinoma Mixed Type" in g for g in groups)


def test_histology_counts_reports_distinct_values_most_common_first(tmp_path):
    csv = _heterogeneous_histology_csv(tmp_path)
    samples = [f"S{i:02d}" for i in range(1, 21)]
    raw = parse_raw_metadata(csv, tmp_path, samples)
    counts = raw.histology_counts()
    assert sum(n for _, n in counts) == 19            # S01's missing value excluded
    assert all(counts[i][1] >= counts[i + 1][1] for i in range(len(counts) - 1))


def test_coarser_age_bands_fix_the_too_few_labelled_problem(tmp_path):
    raw = parse_raw_metadata(FIXTURES / "mini_metadata.csv", tmp_path, SAMPLES)
    bands = quantile_age_bands(raw.ages(), 1)      # one age band -> just Histology x Stage
    md = stratify(raw, bands, min_group_n=10)
    gcs = md.group_counts()
    assert gcs[0].n_labelled == 80                 # everybody now gets a label
    assert all(gc.group != NO_GROUP for gc in gcs)  # nothing left unbanded this time
    assert set(md.labelled_samples(ALL_GROUPS).values()) == {"Good", "Poor"}


# --------------------------------------------------------------------------- #
# API end to end (upload endpoints)
# --------------------------------------------------------------------------- #
def test_api_full_flow(client, sid):
    with open(FIXTURES / "mini_sjdat_genes.rds", "rb") as fh:
        r = client.post(f"/api/session/{sid}/sjdat/gene_matrix", files={"file": fh})
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "gene_matrix" and r.json()["n_features"] == 60

    with open(FIXTURES / "mini_metadata.csv", "rb") as fh:
        r = client.post(f"/api/session/{sid}/metadata",
                        files={"file": ("mini_metadata.csv", fh, "text/csv")})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_matched"] == 80
    # auto-stratified with the classic bands on load
    assert body["age_bands"]["labels"] == ["(30-50]", "(50-70]", "(>70)"]
    assert body["min_group_n"] == DEFAULT_MIN_GROUP_N

    groups = client.get(f"/api/session/{sid}/groups").json()["groups"]
    assert groups[0]["group"] == "__all__" and groups[0]["n_labelled"] == 0
    assert groups[-1]["group"] == "(no Group)"     # pinned to the bottom

    # ask for (and apply) a data-driven, coarser alternative
    r = client.post(f"/api/session/{sid}/age-bands/suggest", json={"n_bands": 1})
    assert r.status_code == 200, r.text
    suggested = r.json()
    assert suggested["include_lowest"] is True

    r = client.post(f"/api/session/{sid}/stratify", json={
        "edges": suggested["edges"], "include_lowest": suggested["include_lowest"],
        "min_group_n": 10,
    })
    assert r.status_code == 200, r.text
    assert r.json()["min_group_n"] == 10

    groups = client.get(f"/api/session/{sid}/groups").json()["groups"]
    assert groups[0]["n_labelled"] == 80
    assert all(g["group"] != "(no Group)" for g in groups)  # nothing left unbanded now

    r = client.post(f"/api/session/{sid}/select",
                    json={"group": "__all__", "n_min": 3, "x_min": 0, "top_n": 25, "n_cv": 5})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_selected"] == 25
    assert body["n_good"] + body["n_poor"] == body["n_group_samples"]

    r = client.post(f"/api/session/{sid}/cross-validate", json={})
    assert r.status_code == 200, r.text
    assert r.json()["n_splits"] == 5
    assert 0.0 <= r.json()["auc"] <= 1.0

    r = client.post(f"/api/session/{sid}/model")
    assert r.status_code == 200, r.text
    assert r.json()["saved"] is True
    assert len(r.json()["features"]) == 25

    st = client.get(f"/api/session/{sid}/state").json()
    assert st["has_model"] is True and st["active_sjdat"] == "gene_matrix"
    assert st["metadata"]["min_group_n"] == 10
    opts = {o["kind"]: o for o in st["sjdat_options"]}
    assert opts["gene_matrix"]["loaded"] and not opts["rrs_scores"]["loaded"]


def test_api_stratify_can_drop_or_merge_histology(client, sid, tmp_path):
    with open(FIXTURES / "mini_sjdat_genes.rds", "rb") as fh:
        client.post(f"/api/session/{sid}/sjdat/gene_matrix", files={"file": fh})
    csv = _heterogeneous_histology_csv(tmp_path)
    with open(csv, "rb") as fh:
        r = client.post(f"/api/session/{sid}/metadata", files={"file": ("h.csv", fh, "text/csv")})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["use_histology"] is True and body["histology_map"] == {}
    assert {c["value"] for c in body["histology_counts"]} == {
        "Adenocarcinoma Intestinal Type", "Adenocarcinoma Diffuse Type",
        "Adenocarcinoma Mixed Type", "Adenocarcinoma NOS",
    }

    suggest = client.post(f"/api/session/{sid}/age-bands/suggest", json={"n_bands": 1}).json()
    r = client.post(f"/api/session/{sid}/stratify", json={
        "edges": suggest["edges"], "include_lowest": suggest["include_lowest"], "min_group_n": 1,
        "use_histology": False,
    })
    assert r.status_code == 200, r.text
    assert r.json()["use_histology"] is False
    groups = client.get(f"/api/session/{sid}/groups").json()["groups"]
    # dropping histology collapses everyone into one Stage | age-band group
    assert sum(1 for g in groups if g["group"] not in ("__all__", "(no Group)")) == 1

    hmap = {"Adenocarcinoma Intestinal Type": "merged", "Adenocarcinoma Diffuse Type": "merged"}
    r = client.post(f"/api/session/{sid}/stratify", json={
        "edges": suggest["edges"], "include_lowest": suggest["include_lowest"], "min_group_n": 1,
        "use_histology": True, "histology_map": hmap,
    })
    assert r.status_code == 200, r.text
    assert r.json()["histology_map"] == hmap
    groups = client.get(f"/api/session/{sid}/groups").json()["groups"]
    assert sum(1 for g in groups if g["group"] not in ("__all__", "(no Group)")) == 3


def test_sjdat_narrowed_to_metadata_intersection(client, sid):
    """A matrix sample with no row in the sample-metadata table is excluded
    from the working sjdat — mirrors the SJV/NSJCG guarantee that every
    loaded data file is narrowed to the samples the metadata actually
    covers. Reloading metadata (a broader file, here) recomputes the
    narrowing from the pristine matrix rather than compounding it."""
    with open(FIXTURES / "mini_sjdat_genes.rds", "rb") as fh:
        r = client.post(f"/api/session/{sid}/sjdat/gene_matrix", files={"file": fh})
    assert r.json()["n_samples"] == 80          # no metadata yet — unnarrowed

    with open(FIXTURES / "mini_metadata_subset.csv", "rb") as fh:
        r = client.post(f"/api/session/{sid}/metadata",
                        files={"file": ("m.csv", fh, "text/csv")})
    assert r.status_code == 200, r.text
    assert r.json()["n_matched"] == 70
    assert any("10 " in w and "excluded" in w for w in r.json()["warnings"]), r.json()["warnings"]

    st = client.get(f"/api/session/{sid}/state").json()
    opt = next(o for o in st["sjdat_options"] if o["kind"] == "gene_matrix")
    assert opt["n_samples"] == 70                # narrowed, matches the metadata

    # coarser bands so this small cohort actually gets Good/Poor labels (see
    # test_classic_bands_can_fragment_a_small_cohort_below_min_group_n)
    suggested = client.post(f"/api/session/{sid}/age-bands/suggest", json={"n_bands": 1}).json()
    client.post(f"/api/session/{sid}/stratify", json={
        "edges": suggested["edges"], "include_lowest": suggested["include_lowest"], "min_group_n": 10,
    })
    # a select over "every sample" now only ever sees the 70 that have metadata
    r = client.post(f"/api/session/{sid}/select",
                    json={"group": "__all__", "n_min": 3, "x_min": 0, "top_n": 10, "n_cv": 5})
    assert r.status_code == 200, r.text
    assert r.json()["n_group_samples"] <= 70

    # reloading the *broader* (80-sample) metadata recomputes from the
    # pristine matrix, not from the already-narrowed 70
    with open(FIXTURES / "mini_metadata.csv", "rb") as fh:
        r = client.post(f"/api/session/{sid}/metadata",
                        files={"file": ("m.csv", fh, "text/csv")})
    assert r.status_code == 200, r.text
    assert r.json()["n_matched"] == 80
    st = client.get(f"/api/session/{sid}/state").json()
    opt = next(o for o in st["sjdat_options"] if o["kind"] == "gene_matrix")
    assert opt["n_samples"] == 80                # recovered — not stuck at 70


def test_api_select_needs_data(client, sid):
    r = client.post(f"/api/session/{sid}/select",
                    json={"group": "__all__", "n_min": 3})
    assert r.status_code == 409


def test_api_stratify_rejects_bad_edges(client, sid):
    with open(FIXTURES / "mini_metadata.csv", "rb") as fh:
        client.post(f"/api/session/{sid}/metadata",
                    files={"file": ("mini_metadata.csv", fh, "text/csv")})
    r = client.post(f"/api/session/{sid}/stratify",
                    json={"edges": [50, 30, None], "include_lowest": False, "min_group_n": 20})
    assert r.status_code == 422


# --------------------------------------------------------------------------- #
# gene-set feature selection (Type genes / Upload list / Pathway tabs) — the
# SJSurv-side port of sjvc's geneset/features endpoints; see
# sjsurv/api/routes.py's set_geneset()/build_features()/select_geneset().
# --------------------------------------------------------------------------- #
SJVC_FIXTURES = FIXTURES.parents[1] / "sjvc" / "fixtures"


def _loaded_gene_session(client, sid):
    with open(FIXTURES / "mini_sjdat_genes.rds", "rb") as fh:
        r = client.post(f"/api/session/{sid}/sjdat/gene_matrix", files={"file": fh})
    assert r.status_code == 200, r.text
    with open(FIXTURES / "mini_metadata.csv", "rb") as fh:
        r = client.post(f"/api/session/{sid}/metadata",
                        files={"file": ("mini_metadata.csv", fh, "text/csv")})
    assert r.status_code == 200, r.text


def test_api_geneset_typed_flow_matches_mad_flow_shape(client, sid):
    """Typed genes -> features -> select-geneset -> cross-validate -> model
    all succeed and return the exact same response shapes the MAD path does
    (test_api_full_flow, above) — model.py, CVResponse, ModelResponse are
    untouched by the gene-set path."""
    _loaded_gene_session(client, sid)
    # the classic default bands leave nobody labelled for this fixture (see
    # test_api_full_flow) — apply the same coarser, data-driven split it uses
    suggested = client.post(f"/api/session/{sid}/age-bands/suggest", json={"n_bands": 1}).json()
    client.post(f"/api/session/{sid}/stratify", json={
        "edges": suggested["edges"], "include_lowest": suggested["include_lowest"],
        "min_group_n": 10,
    })

    r = client.post(f"/api/session/{sid}/geneset",
                     json={"mode": "typed", "text": "GENE001, GENE002, GENE003, BOGUS"})
    assert r.status_code == 200, r.text
    body = r.json()
    assert set(body["matched"]) == {"GENE001", "GENE002", "GENE003"}
    assert body["unmatched"] == ["BOGUS"]

    r = client.post(f"/api/session/{sid}/features", json={"condense": False})
    assert r.status_code == 200, r.text
    assert r.json()["n_features"] == 3 and r.json()["feature_kind"] == "gene"

    r = client.post(f"/api/session/{sid}/select-geneset", json={"group": "__all__"})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_selected"] == 3
    assert body["n_candidates"] == body["n_after_coverage"] == 3   # no coverage filter applies
    assert body["n_good"] + body["n_poor"] == body["n_group_samples"]

    r = client.post(f"/api/session/{sid}/cross-validate", json={})
    assert r.status_code == 200, r.text
    assert 0.0 <= r.json()["auc"] <= 1.0

    r = client.post(f"/api/session/{sid}/model")
    assert r.status_code == 200, r.text
    assert r.json()["saved"] is True
    assert len(r.json()["features"]) == 3

    st = client.get(f"/api/session/{sid}/state").json()
    assert st["has_model"] is True and st["has_geneset"] is True and st["has_features"] is True


def test_api_features_mad_top_n_refines_resolved_gene_set(client, sid):
    """mad_top_n narrows an already-resolved gene set down to its most
    variable members — GENE058 is all-zero (constant) in the fixture, so it's
    the one dropped when asked for the top 2 of the 3 typed genes."""
    _loaded_gene_session(client, sid)

    r = client.post(f"/api/session/{sid}/geneset",
                     json={"mode": "typed", "text": "GENE001, GENE002, GENE058"})
    assert r.status_code == 200, r.text
    assert set(r.json()["matched"]) == {"GENE001", "GENE002", "GENE058"}

    r = client.post(f"/api/session/{sid}/features", json={"condense": False, "mad_top_n": 2})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_features"] == 2
    assert set(body["feature_preview"]) == {"GENE001", "GENE002"}
    assert any("ranked" in w and "kept the top 2" in w for w in body["warnings"])

    # mad_top_n omitted, or >= the resolved count, is unchanged
    unfiltered = client.post(f"/api/session/{sid}/features", json={"condense": False}).json()
    assert unfiltered["n_features"] == 3 and unfiltered["warnings"] == []


def test_api_geneset_junction_level_needs_reference(client, sid):
    """A junction-level sjdat (no GENCODE / junction-metadata loaded) refuses
    /geneset with the same 409 message sjvc uses — the two typed/list/pathway
    tabs are unusable there until a reference is loaded, same as 2D View."""
    with open(SJVC_FIXTURES / "mini_junctions.rds", "rb") as fh:
        r = client.post(f"/api/session/{sid}/sjdat/junction_counts", files={"file": fh})
    assert r.status_code == 200, r.text

    r = client.post(f"/api/session/{sid}/geneset", json={"mode": "typed", "text": "TESTG1"})
    assert r.status_code == 409
    assert "GENCODE" in r.json()["detail"]


def test_api_select_geneset_needs_features_first(client, sid):
    _loaded_gene_session(client, sid)
    r = client.post(f"/api/session/{sid}/select-geneset", json={"group": "__all__"})
    assert r.status_code == 409
    assert "gene set" in r.json()["detail"]


def test_service_select_features_n_min_none_skips_coverage_filter():
    """n_min=None (what the frontend sends when N/X are hidden for a
    gene-level sjdat) ranks every row by MAD — no coverage prefilter."""
    sjdat, labels = _synthetic_sjdat(n_feat=80, n_samp=60)
    sample_ids = list(labels)
    sel = select_features(sjdat, sample_ids, n_min=None, x_min=0, top_n=10)
    assert sel.n_candidates == sel.n_after_coverage == sjdat.n_features
    assert sel.n_features == 10
