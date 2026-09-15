import shutil

import numpy as np
import pytest

from tests.sjv.conftest import FIXTURES

from sjv.services.groups import parse_sample_metadata
from sjv.services.rds import RdsError, RdsMatrix, _parse_rownames, load_rds


def _matrix(tmp_path):
    src = tmp_path / "upload.rds"
    shutil.copy(FIXTURES / "mini.rds", src)
    return load_rds(src, tmp_path)


def _synth(values, samples):
    """A small in-memory RdsMatrix with one synthetic junction per row."""
    values = np.asarray(values, dtype=float)
    rownames = [f"chr1:{100 + 10 * i}-{200 + 10 * i}:+" for i in range(values.shape[0])]
    return RdsMatrix(samples=list(samples), values=values, _c=_parse_rownames(rownames), sparse=False)


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


def test_rownames_column_used_as_sample_id_when_no_sample_id_column(tmp_path):
    """recount3/TCGA sample-metadata tables commonly carry the sample UUID
    only as the data.frame's own row names — pyreadr surfaces those as a
    column literally named 'rownames' once promoted (see sjv.services.rds.
    _read_named_table); this must be usable as the sample_id."""
    p = tmp_path / "m.csv"
    p.write_text("rownames,condition\nsample_A,tumor\nsample_B,normal\n")
    meta = parse_sample_metadata(p, tmp_path, KNOWN)
    assert meta.columns == ["condition"]
    assert meta.n_matched == 2


def test_real_sample_id_column_wins_over_rownames(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("rownames,sample_id,condition\nbogus,sample_A,tumor\n")
    meta = parse_sample_metadata(p, tmp_path, KNOWN)
    assert meta.n_matched == 1 and "sample_A" in meta.rows and "bogus" not in meta.rows


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
    m = _synth(
        [
            [10.0, 20.0, 30.0],   # median 20
            [0.0, 8.0, np.nan],   # only 8 qualifies
            [0.0, 0.0, 0.0],      # nothing qualifies -> NaN
        ],
        ["s1", "s2", "s3"],
    )
    med = m.group_median(["s1", "s2", "s3"])
    assert med[0] == 20
    assert med[1] == 8
    assert np.isnan(med[2])


def test_group_median_two_samples_averages(tmp_path):
    m = _synth([[10.0, 20.0]], ["s1", "s2"])
    assert m.group_median(["s1", "s2"])[0] == 15
