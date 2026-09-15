import numpy as np
import pytest

from tests.conftest import FIXTURES

from app.services.gencode import annotation_from_gtf
from app.services.mad import top_features_by_mad, top_genes_by_mad
from app.services.rds import load_matrix


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
