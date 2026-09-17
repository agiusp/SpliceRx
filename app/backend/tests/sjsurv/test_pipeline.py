import math

import numpy as np
import pytest
from fastapi.testclient import TestClient

from sjsurv.main import app
from sjsurv.services.metadata import (
    ALL_GROUPS,
    CLASSIC_AGE_BANDS,
    COV_AGE,
    COV_HISTOLOGY,
    COV_STAGE,
    DEFAULT_MIN_GROUP_N,
    NO_GROUP,
    AgeBands,
    MetadataError,
    StratifyError,
    bin_ages,
    build_covariates,
    parse_raw_metadata,
    quantile_age_bands,
    stage_label,
    stratify,
)
from sjsurv.services.model import ModelError, cross_validate, train_full
from sjsurv.services.select import (
    Selection, SelectError, rank_features, select_features, select_from_ranking,
)
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


def test_rank_features_reused_across_top_n_matches_one_shot():
    # the split that lets a caller cache the ranking (e.g. a UI "top n"
    # slider) and re-slice it must return exactly what the one-shot
    # select_features(top_n=N) would, for every N, from a single ranking
    d, labels = _synthetic_sjdat()
    ranking = rank_features(d, list(labels), n_min=5, x_min=1)
    for top_n in (1, 20, 10_000):
        sliced = select_from_ranking(d, ranking, top_n)
        one_shot = select_features(d, list(labels), n_min=5, x_min=1, top_n=top_n)
        assert sliced.feature_ids == one_shot.feature_ids
        assert np.array_equal(sliced.values, one_shot.values)


def test_select_rrs_scores_ranks_by_variance_not_mad():
    """RRS scores are bounded [0, 1] and mostly zero, so a row's median is
    almost always exactly 0 and MAD collapses to (near-)0 for a row whose
    only signal is a rare spike — plain variance still picks those up.
    F0/F3 each have one large spike among mostly-zero entries (MAD 0, real
    variance); F1 has a small, consistent spread with no spike (real MAD, but
    the smallest variance of the three varying rows). Top-2 by MAD must
    therefore differ from top-2 by variance, and kind="rrs_scores" must
    select the variance ranking."""
    values = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.9],          # spike -> MAD 0, high variance
        [0.10, 0.15, 0.05, 0.20, 0.10, 0.12],    # small consistent spread -> real MAD, low variance
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],          # constant
        [0.3, 0.3, 0.3, 0.3, 0.3, 0.9],          # spike on a nonzero baseline -> MAD 0, high variance
    ])
    samples = [f"s{i}" for i in range(6)]

    d_gene = Sjdat(kind="gene_matrix", features=["F0", "F1", "F2", "F3"], samples=samples, values=values)
    sel_mad = select_features(d_gene, samples, n_min=None, x_min=0, top_n=2)
    assert "F1" in sel_mad.feature_ids   # the only row with nonzero MAD must be picked

    d_rrs = Sjdat(kind="rrs_scores", features=["F0", "F1", "F2", "F3"], samples=samples, values=values)
    sel_var = select_features(d_rrs, samples, n_min=None, x_min=0, top_n=2)
    assert set(sel_var.feature_ids) == {"F0", "F3"}   # the two spikes, by variance
    assert set(sel_var.feature_ids) != set(sel_mad.feature_ids)


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


def test_vital_status_parsed_into_event_flag(tmp_path):
    """Dead -> event=True, Alive -> event=False, missing/unrecognised -> the
    conservative event=True fallback; no vital-status column at all ->
    event=True for everyone and RawMetadata.has_vital_status is False."""
    csv = tmp_path / "vital.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days,tcga.cgc_case_vital_status\n"
        "S01,Adenocarcinoma,Stage I,55,900,Dead\n"
        "S02,Adenocarcinoma,Stage I,60,1800,Alive\n"
        "S03,Adenocarcinoma,Stage I,65,300,\n"
    )
    raw = parse_raw_metadata(csv, tmp_path, ["S01", "S02", "S03"])
    assert raw.has_vital_status is True
    assert raw.rows["S01"]["event"] is True
    assert raw.rows["S02"]["event"] is False
    assert raw.rows["S03"]["event"] is True   # missing status -> conservative fallback

    csv2 = tmp_path / "no_vital.csv"
    csv2.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days\n"
        "S01,Adenocarcinoma,Stage I,55,900\n"
    )
    raw2 = parse_raw_metadata(csv2, tmp_path, ["S01"])
    assert raw2.has_vital_status is False
    assert raw2.rows["S01"]["event"] is True


def test_stratify_uses_cox_median_and_leaves_early_censored_samples_unlabelled(tmp_path):
    """The reported bug: a plain median of observed follow-up times treats a
    still-living, early-censored sample as if its short observed time were
    its true (Poor) survival. Six same-Group samples with a known Cox/KM
    median survival of 40 (computed by hand and cross-checked against
    lifelines directly): the naive np.median of the raw follow-up times
    (10,20,30,40,50,60) would be 35, which would wrongly call the 30-day
    censored sample "Poor". The Cox-based median must instead land on 40 (the
    censoring at 30 and 50 pushes it later), and that 30-day censored sample
    must come out unlabelled (None), not "Poor" — while a death at the same
    observed time (or earlier) still correctly comes out "Poor"."""
    csv = tmp_path / "censored_group.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days,tcga.cgc_case_vital_status\n"
        "S01,Adenocarcinoma,Stage I,45,10,Dead\n"
        "S02,Adenocarcinoma,Stage I,45,20,Dead\n"
        "S03,Adenocarcinoma,Stage I,45,30,Alive\n"
        "S04,Adenocarcinoma,Stage I,45,40,Dead\n"
        "S05,Adenocarcinoma,Stage I,45,50,Alive\n"
        "S06,Adenocarcinoma,Stage I,45,60,Dead\n"
    )
    samples = [f"S0{i}" for i in range(1, 7)]
    raw = parse_raw_metadata(csv, tmp_path, samples)
    md = stratify(raw, CLASSIC_AGE_BANDS, min_group_n=6)

    rows = md.rows
    group = rows["S01"]["group"]
    assert all(rows[s]["group"] == group for s in samples)  # all one Group, as designed

    assert rows["S01"]["surviver"] == "Poor"   # 10, Dead, below the median
    assert rows["S02"]["surviver"] == "Poor"   # 20, Dead, below the median
    assert rows["S03"]["surviver"] is None     # 30, Alive/censored, below the median -> unknown
    assert rows["S04"]["surviver"] == "Good"   # 40, Dead, at the Cox-based median
    assert rows["S05"]["surviver"] == "Good"   # 50, Alive/censored, past the median
    assert rows["S06"]["surviver"] == "Good"   # 60, Dead, past the median

    labelled = md.labelled_samples(ALL_GROUPS)
    assert "S03" not in labelled
    assert len(labelled) == 5


def test_stratify_event_quantile_lowers_the_survival_threshold(tmp_path):
    """Same 6-sample cohort as the Cox-median test above, but with
    event_quantile=0.25 instead of the classic 0.5 median. Fewer deaths are
    needed to define the threshold (1.5 of 6, effectively the 2nd death
    time), so it lands much earlier — 20 instead of 40 (hand-computed and
    cross-checked directly against lifelines) — which reaches a defined
    threshold sooner for a low-mortality group and, here, leaves nobody
    unlabelled even though two samples are censored."""
    csv = tmp_path / "censored_group.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days,tcga.cgc_case_vital_status\n"
        "S01,Adenocarcinoma,Stage I,45,10,Dead\n"
        "S02,Adenocarcinoma,Stage I,45,20,Dead\n"
        "S03,Adenocarcinoma,Stage I,45,30,Alive\n"
        "S04,Adenocarcinoma,Stage I,45,40,Dead\n"
        "S05,Adenocarcinoma,Stage I,45,50,Alive\n"
        "S06,Adenocarcinoma,Stage I,45,60,Dead\n"
    )
    samples = [f"S0{i}" for i in range(1, 7)]
    raw = parse_raw_metadata(csv, tmp_path, samples)
    md = stratify(raw, CLASSIC_AGE_BANDS, min_group_n=6, event_quantile=0.25)

    rows = md.rows
    assert rows["S01"]["surviver"] == "Poor"   # 10, Dead, below the (now-earlier) threshold of 20
    assert rows["S02"]["surviver"] == "Good"   # 20, Dead, at the threshold
    assert rows["S03"]["surviver"] == "Good"   # 30, Alive, past the threshold
    assert rows["S04"]["surviver"] == "Good"
    assert rows["S05"]["surviver"] == "Good"
    assert rows["S06"]["surviver"] == "Good"

    labelled = md.labelled_samples(ALL_GROUPS)
    assert len(labelled) == 6   # nobody unlabelled, unlike the 0.5-quantile case above


def test_stratify_rejects_out_of_range_event_quantile(tmp_path):
    raw = parse_raw_metadata(FIXTURES / "mini_metadata.csv", tmp_path, SAMPLES)
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(StratifyError, match="event_quantile"):
            stratify(raw, CLASSIC_AGE_BANDS, DEFAULT_MIN_GROUP_N, event_quantile=bad)


def test_stratify_falls_back_to_plain_median_without_vital_status(tmp_path):
    """No vital-status column at all -> every sample is treated as a
    confirmed death (the pre-existing, non-censoring-aware behaviour), so the
    Cox-based median collapses to the same value a plain median would give
    and every sample gets a definite label — nobody is left unlabelled purely
    because censoring information doesn't exist for this cohort."""
    csv = tmp_path / "no_vital_group.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days\n"
        + "\n".join(
            f"S0{i},Adenocarcinoma,Stage I,45,{d}" for i, d in enumerate([10, 20, 30, 40, 50, 60], start=1)
        )
    )
    samples = [f"S0{i}" for i in range(1, 7)]
    raw = parse_raw_metadata(csv, tmp_path, samples)
    assert raw.has_vital_status is False
    md = stratify(raw, CLASSIC_AGE_BANDS, min_group_n=6)
    labelled = md.labelled_samples(ALL_GROUPS)
    assert len(labelled) == 6
    assert set(labelled.values()) == {"Good", "Poor"}


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
    # 25 molecular + the default covariates (Histology: 1 constant value,
    # Stage: Early/Late, Age At Diagnosis: 1 numeric) = 25 + 1 + 2 + 1
    assert len(r.json()["features"]) == 29

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
    # 3 molecular + the same 4 default covariates as test_api_full_flow
    assert len(r.json()["features"]) == 7

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


def test_select_features_protein_coding_filter_gene_level():
    d, labels = _synthetic_sjdat()
    restrict_to = {"g0", "g3", "g7"}
    sel = select_features(d, list(labels), n_min=None, x_min=0, top_n=10, restrict_to=restrict_to)
    assert {f.lower() for f in sel.feature_ids} == restrict_to

    with pytest.raises(SelectError):
        select_features(d, list(labels), n_min=None, x_min=0, top_n=10, restrict_to={"nope"})


def test_select_features_protein_coding_filter_junction_level_needs_annotation():
    """Mirrors sjvc's own top_features_by_mad junction-level coding filter —
    same shared helper (junction_rownames_overlapping_genes), same
    requirement that a live GENCODE reference is loaded (the fast junction-
    metadata lookup has no gene-biotype info to filter on)."""
    from sjvc.services.gencode import annotation_from_gtf

    ann = annotation_from_gtf(SJVC_FIXTURES / "mini.gtf", label="test")
    restrict_to = ann.gene_names_of_type("protein_coding")

    # chr1:2201-2999:+ sits inside TESTG1 only (protein_coding); chr2:5000-
    # 6000:+ overlaps no gene at all in this reference.
    rownames = ["chr1:2201-2999:+", "chr2:5000-6000:+"]
    rng = np.random.default_rng(0)
    vals = rng.poisson(3, size=(len(rownames), 20)).astype(float)
    d = Sjdat(
        kind="junction_counts", features=rownames, samples=[f"s{i:02d}" for i in range(20)],
        values=vals, sparse=False,
    )

    with pytest.raises(SelectError):
        select_features(d, d.samples, n_min=None, x_min=0, top_n=10, restrict_to=restrict_to)

    sel = select_features(
        d, d.samples, n_min=None, x_min=0, top_n=10, restrict_to=restrict_to, annotation=ann,
    )
    assert sel.feature_ids == ["chr1:2201-2999:+"]


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


# --------------------------------------------------------------------------- #
# clinical covariates — Stage / Age At Diagnosis / Histology / MSI (when
# present) alongside the selected molecular features
# --------------------------------------------------------------------------- #
def _metadata_with_msi_csv(tmp_path):
    """Adds an MSI-status column (categorical) to the heterogeneous-histology
    fixture, plus a free-text column that should be offered too, just not
    checked by default."""
    csv = tmp_path / "with_msi.csv"
    rows = [
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days,msi_status,notes"
    ]
    hist = ["Adenocarcinoma Intestinal Type", "Adenocarcinoma Diffuse Type"]
    msi = ["MSI-H", "MSS"]
    for i in range(1, 21):
        rows.append(
            f"S{i:02d},{hist[i % 2]},Stage I,{40 + i},{900 + i},{msi[i % 2]},note{i}"
        )
    csv.write_text("\n".join(rows) + "\n")
    return csv


def test_available_covariates_defaults_and_kinds(tmp_path):
    samples = [f"S{i:02d}" for i in range(1, 21)]
    raw = parse_raw_metadata(_metadata_with_msi_csv(tmp_path), tmp_path, samples)
    cols = {c.key: c for c in raw.available_covariates()}

    assert cols[COV_HISTOLOGY].label == "Histology" and cols[COV_HISTOLOGY].default
    assert cols[COV_STAGE].label == "Stage" and cols[COV_STAGE].default
    assert cols[COV_AGE].label == "Age At Diagnosis" and cols[COV_AGE].kind == "numeric"
    assert cols[COV_AGE].default

    # the MSI column is offered *and* checked on by default, purely from its name
    assert cols["msi_status"].kind == "categorical" and cols["msi_status"].default
    assert cols["msi_status"].label == "Msi status"

    # an unrelated free-text column is offered too, but not default-checked
    assert "notes" in cols and cols["notes"].default is False

    # survival_days (what SurviverGroup is split on) must never be offered —
    # it would let the classifier "predict" the label for free
    assert "survival_days" not in cols
    assert "sample_id" not in cols


def test_build_covariates_numeric_and_categorical_shapes(tmp_path):
    samples = [f"S{i:02d}" for i in range(1, 21)]
    raw = parse_raw_metadata(_metadata_with_msi_csv(tmp_path), tmp_path, samples)

    names, values = build_covariates(raw, [COV_AGE, COV_STAGE, "msi_status"], samples)
    assert values.shape[1] == len(samples)
    assert "Age At Diagnosis" in names                    # one numeric row
    assert sum(n.startswith("Stage=") for n in names) == 1  # Stage I -> Early only, one value
    assert sum(n.startswith("Msi status=") for n in names) == 2   # MSI-H / MSS
    age_row = values[names.index("Age At Diagnosis")]
    assert np.allclose(sorted(age_row), sorted(41 + i for i in range(20)))


def test_build_covariates_imputes_missing_numeric_with_the_sample_mean(tmp_path):
    csv = tmp_path / "missing_age.csv"
    csv.write_text(
        "sample_id,tcga.cgc_case_histological_diagnosis,tcga.cgc_case_pathologic_stage,"
        "age_at_diagnosis_years,survival_days\n"
        "S01,Adenocarcinoma,Stage I,40,900\n"
        "S02,Adenocarcinoma,Stage I,,1800\n"     # missing age
        "S03,Adenocarcinoma,Stage I,60,300\n"
    )
    raw = parse_raw_metadata(csv, tmp_path, ["S01", "S02", "S03"])
    names, values = build_covariates(raw, [COV_AGE], ["S01", "S02", "S03"])
    assert values[0].tolist() == [40.0, 50.0, 60.0]      # S02 imputed to the mean of S01/S03


def test_histology_covariate_reuses_the_applied_merge_map(tmp_path):
    csv = _heterogeneous_histology_csv(tmp_path)
    samples = [f"S{i:02d}" for i in range(1, 21)]
    raw = parse_raw_metadata(csv, tmp_path, samples)
    hmap = {
        "Adenocarcinoma Intestinal Type": "Adenocarcinoma (merged)",
        "Adenocarcinoma Diffuse Type": "Adenocarcinoma (merged)",
    }
    names, _ = build_covariates(raw, [COV_HISTOLOGY], samples, histology_map=hmap)
    # 4 raw values, 2 merged together -> 3 distinct Histology= columns, not 4
    assert sum(n.startswith("Histology=") for n in names) == 3
    assert "Histology=Adenocarcinoma (merged)" in names


def test_api_covariates_endpoint_changes_model_feature_count(client, sid):
    _loaded_gene_session(client, sid)
    suggested = client.post(f"/api/session/{sid}/age-bands/suggest", json={"n_bands": 1}).json()
    client.post(f"/api/session/{sid}/stratify", json={
        "edges": suggested["edges"], "include_lowest": suggested["include_lowest"],
        "min_group_n": 10,
    })

    st = client.get(f"/api/session/{sid}/state").json()
    default_keys = {c["key"] for c in st["covariate_columns"] if c["default"]}
    assert default_keys == set(st["selected_covariates"])
    assert {COV_HISTOLOGY, COV_STAGE, COV_AGE} <= default_keys

    r = client.post(f"/api/session/{sid}/select",
                     json={"group": "__all__", "n_min": 3, "x_min": 0, "top_n": 20})
    assert r.status_code == 200, r.text
    n_molecular = r.json()["n_selected"]

    with_defaults = client.post(f"/api/session/{sid}/model").json()
    assert with_defaults["n_features"] > n_molecular

    r = client.post(f"/api/session/{sid}/covariates", json={"columns": []})
    assert r.status_code == 200, r.text
    assert r.json()["selected_covariates"] == []
    none = client.post(f"/api/session/{sid}/model").json()
    assert none["n_features"] == n_molecular == len(none["features"])

    r = client.post(f"/api/session/{sid}/covariates", json={"columns": [COV_AGE]})
    assert r.status_code == 200, r.text
    just_age = client.post(f"/api/session/{sid}/model").json()
    assert just_age["n_features"] == n_molecular + 1
    assert any(f["feature"] == "Age At Diagnosis" for f in just_age["features"])

    r = client.post(f"/api/session/{sid}/covariates", json={"columns": ["not-a-real-key"]})
    assert r.status_code == 422
