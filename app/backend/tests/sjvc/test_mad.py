import numpy as np
import pytest

from tests.sjvc.conftest import FIXTURES

from sjvc.services.gencode import annotation_from_gtf
from sjvc.services.mad import features_from_ranking, rank_by_mad, top_features_by_mad, top_genes_by_mad
from sjvc.services.rds import Matrix, load_matrix


def _matrix(tmp_path):
    return load_matrix(FIXTURES / "mini_junctions.rds", tmp_path)


def test_top_junctions_ranks_by_mad_descending(tmp_path):
    m = _matrix(tmp_path)
    fm = top_features_by_mad(m, top_n=3)
    assert fm.kind == "junction" and fm.n_features == 3

    # cross-check against a hand-computed MAD on log1p(NA->0)
    vals = np.nan_to_num(m.values.astype(float), nan=0.0)
    scored = np.log1p(vals)
    med = np.median(scored, axis=1, keepdims=True)
    mad = np.median(np.abs(scored - med), axis=1)
    ranked = [m.features[i] for i in np.argsort(mad)[::-1]]
    # bad rowname / off-locus junctions aren't excluded from this pool by
    # anything other than parseability, so compare against the parseable ones
    ranked = [r for r in ranked if r != "not_a_junction"]
    assert set(fm.feature_ids) == set(ranked[:3])


def test_rank_by_mad_reused_across_top_n_matches_one_shot(tmp_path):
    # the split that lets a caller cache the ranking (e.g. a UI "top n"
    # slider) and re-slice it must return exactly what the one-shot
    # top_features_by_mad(top_n=N) would, for every N, from a single ranking
    m = _matrix(tmp_path)
    ranking = rank_by_mad(m)
    for top_n in (1, 3, 5, 10_000):
        sliced = features_from_ranking(m, ranking, top_n)
        one_shot = top_features_by_mad(m, top_n=top_n)
        assert sliced.feature_ids == one_shot.feature_ids
        assert np.array_equal(sliced.values, one_shot.values, equal_nan=True)


def test_top_junctions_caps_at_available_rows(tmp_path):
    m = _matrix(tmp_path)
    fm = top_features_by_mad(m, top_n=10_000)
    assert fm.n_features == len([f for f in m.features if f != "not_a_junction"])


def test_top_features_on_already_condensed_matrix_ranks_every_row(tmp_path):
    # rows are gene symbols, not chr:start-end:strand — the output shape of
    # functions/condense_junctions_by_gene.py
    m = _matrix(tmp_path)
    m.features = [f"GENE{i}" for i in range(len(m.features))]
    fm = top_features_by_mad(m, top_n=3)
    assert fm.kind == "gene" and fm.n_features == 3
    assert fm.n_junctions == 0
    assert set(fm.feature_ids).issubset(set(m.features))


def test_rrs_scores_rank_by_variance_not_mad():
    """RRS scores are bounded [0, 1] and mostly zero, so a row's median is
    almost always exactly 0 and MAD (robust to a minority of outliers)
    collapses to (near-)0 for rows whose only signal is a rare spike — plain
    variance still picks those up. GENE0/GENE3 each have one large spike
    among mostly-zero entries (MAD == 0, real variance); GENE1 has a small,
    consistent spread with no outlier (real MAD, but the smallest variance of
    the three varying rows). Requesting the top 2 by MAD must therefore
    differ from the top 2 by variance, and `matrix.kind == "rrs_scores"` must
    select the variance ranking."""
    values = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.9],   # spike -> MAD 0, high variance
        [0.10, 0.15, 0.05, 0.20, 0.10, 0.12],  # small consistent spread -> real MAD, low variance
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # constant -> excluded either way
        [0.3, 0.3, 0.3, 0.3, 0.3, 0.9],   # spike on a nonzero baseline -> MAD 0, high variance
    ])
    m = Matrix(
        samples=[f"s{i}" for i in range(6)],
        features=["GENE0", "GENE1", "GENE2", "GENE3"],
        values=values,
    )

    fm_mad = top_features_by_mad(m, top_n=2)
    assert "GENE1" in fm_mad.feature_ids   # the only row with nonzero MAD must be picked

    m.kind = "rrs_scores"
    fm_var = top_features_by_mad(m, top_n=2)
    assert set(fm_var.feature_ids) == {"GENE0", "GENE3"}   # the two spikes, by variance
    assert set(fm_var.feature_ids) != set(fm_mad.feature_ids)


def test_top_features_excludes_constant_rows_rather_than_padding_top_n():
    """Asking for more features than actually vary must not pad the result
    out with constant (MAD == 0) rows — a request for "top 500 by MAD" on a
    matrix with only, say, 3 non-constant rows should return 3, not 500 (497
    of them arbitrary all-zero rows tied at MAD 0)."""
    # each has 6 distinct values (no majority-repeated value, so MAD > 0 —
    # MAD is robust to a *minority* of outliers, so a mostly-repeated value
    # like [1,1,1,1,1,9] actually has MAD == 0, unlike ordinary variance)
    varying = [[1, 2, 3, 4, 8, 1], [10, 3, 7, 1, 12, 5], [6, 1, 9, 2, 11, 4]]
    constant = [[0, 0, 0, 0, 0, 0]] * 5 + [[3, 3, 3, 3, 3, 3]] * 2  # two different constants
    values = np.array(varying + constant, dtype=float)
    m = Matrix(
        samples=[f"s{i}" for i in range(6)],
        features=[f"GENE{i}" for i in range(len(values))],
        values=values,
    )
    fm = top_features_by_mad(m, top_n=100)  # far more than the 3 that vary at all
    assert fm.n_features == 3
    assert set(fm.feature_ids) == {"GENE0", "GENE1", "GENE2"}


def test_top_features_all_constant_raises():
    m = Matrix(
        samples=["s1", "s2", "s3"],
        features=["G1", "G2"],
        values=np.array([[0.0, 0.0, 0.0], [3.0, 3.0, 3.0]]),
    )
    with pytest.raises(ValueError, match="constant"):
        top_features_by_mad(m, top_n=5)


def test_top_features_needs_some_rows(tmp_path):
    m = _matrix(tmp_path)
    m.features = []
    m.values = m.values[:0, :]
    with pytest.raises(ValueError):
        top_features_by_mad(m, top_n=5)


def test_top_genes_ranks_genome_wide_not_just_a_preselection(tmp_path):
    m = _matrix(tmp_path)
    ann = annotation_from_gtf(FIXTURES / "mini.gtf", label="mini")
    fm = top_genes_by_mad(m, ann, top_n=3)
    assert fm.kind == "gene"
    # OTHERG was never part of any typed/pathway selection in other tests but
    # must be considered here, since MAD mode ranks every gene in the reference
    assert set(fm.feature_ids) == {"TESTG1", "TESTG2", "OTHERG"}
    assert {fm.label(g) for g in fm.feature_ids} == {"TESTG1", "TESTG2", "OTHERG"}


def test_top_genes_needs_annotation_overlap(tmp_path):
    m = _matrix(tmp_path)
    m.features = ["chr9:1-2:+"] * len(m.features)
    ann = annotation_from_gtf(FIXTURES / "mini.gtf", label="mini")
    with pytest.raises(ValueError):
        top_genes_by_mad(m, ann, top_n=3)
