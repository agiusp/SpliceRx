import numpy as np
import pytest

from app.services.features import FeatureMatrix, preprocess
from app.services.heatmap import build as build_heatmap
from app.services.projection import ProjectionError, pca, umap


def _fm(n_samples=12, n_feat=8, seed=0):
    rng = np.random.default_rng(seed)
    half = n_samples // 2
    base = rng.poisson(5, (n_feat, n_samples)).astype(float)
    base[:3, :half] += 30      # planted group signal
    base[3:6, half:] += 30
    return FeatureMatrix(
        kind="junction",
        feature_ids=[f"j{i}" for i in range(n_feat)],
        samples=[f"s{i:02d}" for i in range(n_samples)],
        values=base,
        n_junctions=n_feat,
    )


def test_pca_shape_variance_and_sign_determinism():
    prep = preprocess(_fm(), standardize=True)
    a = pca(prep.X, prep.samples)
    b = pca(prep.X, prep.samples)
    assert a.coords.shape[0] == 12
    assert a.explained_variance == sorted(a.explained_variance, reverse=True)
    assert np.allclose(a.coords, b.coords)          # deterministic
    assert a.axis_labels[0].startswith("PC1")


def test_pca_needs_two_samples():
    with pytest.raises(ProjectionError):
        pca(np.zeros((1, 5)), ["s0"])


def test_umap_2d_and_seeded():
    prep = preprocess(_fm(), standardize=True)
    e1 = umap(prep.X, prep.samples, n_neighbors=5)
    e2 = umap(prep.X, prep.samples, n_neighbors=5)
    assert e1.coords.shape == (12, 2)
    assert np.allclose(e1.coords, e2.coords)


def test_heatmap_cluster_vs_group_order():
    fm = _fm()
    clustered = build_heatmap(fm, row_zscore=True, order="cluster")
    assert clustered.col_dendro is not None
    assert set(clustered.samples) == set(fm.samples)

    groups = ["A"] * 6 + ["B"] * 6
    grouped = build_heatmap(fm, row_zscore=True, order="group", group_values=groups)
    assert grouped.col_dendro is None
    assert grouped.samples[:6] == fm.samples[:6]      # stable within group


def test_heatmap_row_cap_and_warning():
    fm = _fm(n_feat=700)
    res = build_heatmap(fm, row_zscore=False, order="cluster")
    assert len(res.feature_ids) == 500
    assert any("highest-variance" in w for w in res.warnings)
