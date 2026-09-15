import numpy as np
import pytest

from app.services.features import preprocess, select_and_maybe_condense
from app.services.junctions import JunctionGeneMap
from app.services.rds import Matrix


def _matrix():
    feats = ["jA1", "jA2", "jB1", "jB2", "jB3"]
    vals = np.array([
        [10, 0, 5],      # jA1
        [0, 3, np.nan],  # jA2
        [1, 1, 1],       # jB1
        [0, 0, 2],       # jB2
        [4, 0, 0],       # jB3
    ], dtype=float)
    return Matrix(samples=["s1", "s2", "s3"], features=feats, values=vals)


def _map():
    return JunctionGeneMap(
        by_rowname={"jA1": ["A"], "jA2": ["A"], "jB1": ["B"], "jB2": ["B"], "jB3": ["B"]},
        n_bad_rownames=0, bad_sample=[],
    )


def test_condense_matches_R_apply():
    fm = select_and_maybe_condense(_matrix(), _map(), condense=True)
    assert fm.kind == "gene"
    row = {g: fm.values[i] for i, g in enumerate(fm.feature_ids)}
    # A: jA1>0 -> [1,0,1] ; jA2>0 (NA->0) -> [0,1,0] ; sum -> [1,1,1]
    assert list(row["A"]) == [1.0, 1.0, 1.0]
    # B: jB1 [1,1,1] + jB2 [0,0,1] + jB3 [1,0,0] -> [2,1,2]
    assert list(row["B"]) == [2.0, 1.0, 2.0]


def test_no_condense_keeps_junction_rows():
    fm = select_and_maybe_condense(_matrix(), _map(), condense=False)
    assert fm.kind == "junction" and fm.n_features == 5


def test_preprocess_drops_constant_and_standardizes():
    fm = select_and_maybe_condense(_matrix(), _map(), condense=False)
    prep = preprocess(fm, standardize=True)
    # jB1 is constant -> dropped
    assert "jB1" not in prep.feature_ids and prep.dropped == 1
    # standardized: each feature (column of X) ~ mean 0
    assert np.allclose(prep.X.mean(axis=0), 0, atol=1e-9)
    assert prep.X.shape == (3, 4)
