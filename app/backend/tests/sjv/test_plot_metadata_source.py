"""Tests for the sashimi plot's `annotation_source` option: classify arcs from
the cohort's own junction-metadata table (the recount3/STAR-aligner
`annotated` column) instead of live GENCODE transcript matching, falling back
to the GENCODE classification for any junction the table has no row for.
"""
import warnings

import pytest

from tests.sjv.conftest import FIXTURES

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from sjv.main import app
    return TestClient(app)


@pytest.fixture
def ready(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/rds",
                files={"file": ("m.rds", open(FIXTURES / "mini.rds", "rb"))})
    client.post(f"/api/session/{sid}/gtf",
                files={"file": ("mini.gtf", open(FIXTURES / "mini.gtf", "rb"))})
    return sid


def _metadata_csv(tmp_path):
    """Covers two of mini.rds's junctions:
    - chr1:1201-1999:+ marked `annotated=1` (agrees with the live GENCODE
      classification, which also calls it "annotated").
    - chr1:2201-3999:+ marked not annotated but with only its donor (left, +
      strand) site known -> classify_unannotated should call it "alt_3p",
      disagreeing with the live GENCODE classification ("exon_skipping").
    Every other junction in mini.rds is deliberately absent, to exercise the
    per-arc GENCODE fallback."""
    csv = tmp_path / "junction_metadata.csv"
    csv.write_text(
        "seqnames,start,end,strand,annotated,left_annotated,right_annotated,"
        "gencode_gene_id,gencode_gene_name\n"
        "chr1,1201,1999,+,1,g29,g29,ENSGX,TESTG1\n"
        "chr1,2201,3999,+,0,g29,0,ENSGX,TESTG1\n"
    )
    return csv


def test_upload_junction_metadata(client, ready, tmp_path):
    r = client.post(f"/api/session/{ready}/junction-metadata",
                     files={"file": ("jm.csv", open(_metadata_csv(tmp_path), "rb"))})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_rows"] == 2
    assert body["has_annotation_detail"] is True

    st = client.get(f"/api/session/{ready}/state").json()
    assert st["has_junction_metadata"] is True
    assert st["junction_metadata_has_detail"] is True


def test_metadata_source_requires_the_table_loaded(client, ready):
    r = client.post(f"/api/session/{ready}/plot",
                     json={"genes": ["TESTG1"], "sample": "sample_A", "annotation_source": "metadata"})
    assert r.status_code == 409
    assert "junction metadata" in r.json()["detail"]


def test_metadata_source_overrides_and_falls_back(client, ready, tmp_path):
    client.post(f"/api/session/{ready}/junction-metadata",
                files={"file": ("jm.csv", open(_metadata_csv(tmp_path), "rb"))})

    r = client.post(f"/api/session/{ready}/plot",
                     json={"genes": ["TESTG1", "OTHERG"], "sample": "sample_A",
                           "annotation_source": "metadata"})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["annotation_source"] == "metadata"
    by_id = {a["id"]: a for a in body["arcs"]}

    # in the table, and agrees with GENCODE
    assert by_id["chr1:1201-1999:+"]["category"] == "annotated"
    assert by_id["chr1:1201-1999:+"]["category_source"] == "metadata"

    # in the table, and *disagrees* with GENCODE (which calls this exon_skipping)
    assert by_id["chr1:2201-3999:+"]["category"] == "alt_3p"
    assert by_id["chr1:2201-3999:+"]["category_source"] == "metadata"

    # not in the table -> falls back to the live GENCODE classification
    assert by_id["chr1:1150-1999:+"]["category"] == "alt_5p"
    assert by_id["chr1:1150-1999:+"]["category_source"] == "gencode"

    assert any("6 of 8 junction" in w for w in body["warnings"])


def test_gencode_source_is_still_the_default(client, ready, tmp_path):
    """Loading the metadata table must not change anything unless a caller
    actually asks for it."""
    client.post(f"/api/session/{ready}/junction-metadata",
                files={"file": ("jm.csv", open(_metadata_csv(tmp_path), "rb"))})
    r = client.post(f"/api/session/{ready}/plot",
                     json={"genes": ["TESTG1"], "sample": "sample_A"})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["annotation_source"] == "gencode"
    by_id = {a["id"]: a for a in body["arcs"]}
    assert by_id["chr1:2201-3999:+"]["category"] == "exon_skipping"
    assert by_id["chr1:2201-3999:+"]["category_source"] == "gencode"
