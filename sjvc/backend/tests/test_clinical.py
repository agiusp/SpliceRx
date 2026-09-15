import numpy as np
import pytest

from tests.conftest import FIXTURES

from app.services.clinical import parse_clinical
from app.services.rds import MatrixError

SAMPLES = [f"s{i:02d}" for i in range(1, 13)]


def _clin(tmp_path):
    return parse_clinical(FIXTURES / "mini_clinical.csv", tmp_path, SAMPLES)


def test_column_typing_and_eligibility(tmp_path):
    clin = _clin(tmp_path)
    by = {c.name: c for c in clin.columns}
    assert by["subtype"].type == "categorical" and by["subtype"].eligible
    assert by["stage"].type == "categorical"          # low-cardinality numeric -> categorical
    assert by["stage"].eligible
    assert by["age"].type == "numeric" and by["age"].eligible          # numeric always eligible
    assert by["tmb"].type == "numeric" and by["tmb"].eligible
    # 'notes' is all-distinct text -> not eligible
    assert by["notes"].type == "categorical" and not by["notes"].eligible


def test_values_for_alignment(tmp_path):
    clin = _clin(tmp_path)
    vals, kind = clin.values_for("subtype", ["s01", "s12", "s07"])
    assert (kind, vals) == ("categorical", ["LUAD", "LUSC", "LUSC"])
    age, kind = clin.values_for("age", ["s01", "s02"])
    assert kind == "numeric" and list(age) == [52.0, 61.0]


def test_override_numeric_to_categorical(tmp_path):
    clin = _clin(tmp_path)
    clin.apply_overrides({"age": "categorical"})
    _, kind = clin.values_for("age", ["s01"])
    assert kind == "categorical"


def test_override_categorical_to_numeric_fails(tmp_path):
    clin = _clin(tmp_path)
    with pytest.raises(MatrixError):
        clin.apply_overrides({"subtype": "numeric"})


def test_malformed_table_raises_matrix_error_not_a_crash(tmp_path):
    # a comment/header line with no tabs followed by real tab-separated rows
    # (e.g. a GTF mistakenly uploaded as the clinical table) used to escape
    # parse_clinical as an unhandled pandas.errors.ParserError (-> a raw 500
    # with no useful detail); it must come back as a clean MatrixError.
    p = tmp_path / "not_really.tsv"
    p.write_text("##a comment line with no tabs\nsample_id\tgroup\ts01\tx\n")
    with pytest.raises(MatrixError, match="tab"):
        parse_clinical(p, tmp_path, SAMPLES)


def test_missing_sample_id_column(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("id,group\ns01,x\n")
    with pytest.raises(MatrixError, match="sample_id"):
        parse_clinical(p, tmp_path, SAMPLES)


def test_rds_clinical_with_leading_id_columns_keeps_every_row(tmp_path):
    # mini_clinical.rds is a TCGA-style table: leading patient_id / specimen_id
    # columns, specimen_id repeats (2 samples per patient). Every one of the 12
    # sample rows must survive parsing — none dropped, none flagged unmatched.
    clin = parse_clinical(FIXTURES / "mini_clinical.rds", tmp_path, SAMPLES)
    assert clin.n_rows == 12
    assert clin.unmatched_samples == []
    assert all(s in clin.rows for s in SAMPLES)
    vals, _ = clin.values_for("subtype", ["s08"])
    assert vals == ["LUSC"]


def test_categorical_eligibility_threshold_is_point_nine_n(tmp_path):
    # 10 samples: a categorical column with 9 distinct values (0.9*10 = 9) is NOT
    # eligible (must be strictly below); 8 distinct IS.
    ten = [f"s{i:02d}" for i in range(1, 11)]
    rows_9 = "sample_id,c9,c8\n" + "".join(
        f"{s},v{i if i < 9 else 8},w{i if i < 8 else 7}\n" for i, s in enumerate(ten)
    )
    p = tmp_path / "t.csv"
    p.write_text(rows_9)
    clin = parse_clinical(p, tmp_path, ten)
    by = {c.name: c for c in clin.columns}
    assert by["c9"].n_unique == 9 and not by["c9"].eligible
    assert by["c8"].n_unique == 8 and by["c8"].eligible
