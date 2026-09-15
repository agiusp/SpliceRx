import shutil

import numpy as np
import pytest

from tests.conftest import FIXTURES

from app.services.groups import parse_sample_metadata
from app.services.rds import RdsError, load_rds


def _matrix(tmp_path):
    src = tmp_path / "upload.rds"
    shutil.copy(FIXTURES / "mini.rds", src)
    return load_rds(src, tmp_path)


KNOWN = ["sample_A", "sample_B", "sample_C"]


def test_finds_sample_id_column_and_lists_others(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text(
        "sample_id,condition,batch\n"
        "sample_A,tumor,1\n"
        "sample_B,tumor,2\n"
        "sample_C,normal,1\n"
    )
    meta = parse_sample_metadata(p, tmp_path, KNOWN)
    assert set(meta.columns) == {"condition", "batch"}
    assert meta.n_matched == 3
    assert meta.strat_values("condition") == {
        "tumor": ["sample_A", "sample_B"],
        "normal": ["sample_C"],
    }
    assert meta.strat_values("batch") == {"1": ["sample_A", "sample_C"], "2": ["sample_B"]}


def test_sample_id_column_position_and_case_insensitive(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("condition\tSample_ID\ntumor\tsample_A\nnormal\tsample_B\n")
    meta = parse_sample_metadata(p, tmp_path, KNOWN)
    assert meta.columns == ["condition"]
    assert meta.strat_values("condition") == {"tumor": ["sample_A"], "normal": ["sample_B"]}


def test_missing_sample_id_column_raises(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("name,group\nsample_A,x\n")
    with pytest.raises(RdsError, match="sample_id"):
        parse_sample_metadata(p, tmp_path, KNOWN)


def test_unmatched_sample_ids_reported(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("sample_id,grp\nsample_A,x\nghost,x\n")
    meta = parse_sample_metadata(p, tmp_path, KNOWN)
    assert meta.strat_values("grp") == {"x": ["sample_A"]}
    assert meta.unmatched_samples == ["ghost"]


def test_no_matching_sample_id_raises(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("sample_id,grp\nnope,x\n")
    with pytest.raises(RdsError):
        parse_sample_metadata(p, tmp_path, KNOWN)


def test_missing_values_excluded_from_a_group(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("sample_id,grp\nsample_A,x\nsample_B,\nsample_C,NA\n")
    meta = parse_sample_metadata(p, tmp_path, KNOWN)
    assert meta.strat_values("grp") == {"x": ["sample_A"]}


def test_group_median_excludes_na_and_zero(tmp_path):
    m = _matrix(tmp_path)
    m.values = np.array(
        [
            [10.0, 20.0, 30.0],   # median 20
            [0.0, 8.0, np.nan],   # only 8 qualifies
            [0.0, 0.0, 0.0],      # nothing qualifies -> NaN
        ]
    )
    m.samples = ["s1", "s2", "s3"]
    m._col = {"s1": 0, "s2": 1, "s3": 2}
    med = m.group_median(["s1", "s2", "s3"])
    assert med[0] == 20
    assert med[1] == 8
    assert np.isnan(med[2])


def test_group_median_two_samples_averages(tmp_path):
    m = _matrix(tmp_path)
    m.values = np.array([[10.0, 20.0]])
    m.samples = ["s1", "s2"]
    m._col = {"s1": 0, "s2": 1}
    assert m.group_median(["s1", "s2"])[0] == 15
