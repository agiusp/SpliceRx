"""SJ Lookup — direct per-junction lookup against a junction metadata table."""
import warnings

import pytest

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from sjlookup.main import app
    return TestClient(app)


@pytest.fixture
def sid(client):
    return client.post("/api/session").json()["session_id"]


def _lookup_csv(tmp_path):
    """A junction metadata table carrying every column SJ Lookup surfaces,
    plus a duplicate row (to exercise dedup) and a no-gene-overlap row."""
    csv = tmp_path / "junction_metadata.csv"
    csv.write_text(
        "seqnames,start,end,strand,width,annotated,left_motif,right_motif,"
        "left_annotated,right_annotated,gencode_gene_id,gencode_gene_name\n"
        "chr1,1000,2000,+,1000,1,GT,AG,0,0,ENSG_A.1,GENEA\n"
        "chr1,1000,2000,+,1000,1,GT,AG,0,0,ENSG_A.1,GENEA\n"
        "chr1,3000,4000,-,1000,0,CT,AC,tag1,0,ENSG_B.1,GENEB\n"
        "chr1,5000,6000,+,1000,0,GT,AG,0,0,NA,NA\n"
    )
    return csv


def test_load_and_lookup_mixed_batch(client, sid, tmp_path):
    r = client.post(f"/api/session/{sid}/junction-metadata",
                     files={"file": ("jm.csv", open(_lookup_csv(tmp_path), "rb"))})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_rows"] == 4
    assert body["n_duplicate_rownames"] == 1
    assert any("duplicate" in w for w in body["warnings"])

    st = client.get(f"/api/session/{sid}/state").json()
    assert st == {"has_index": True, "n_rows": 4}

    r = client.post(f"/api/session/{sid}/lookup", json={
        "junctions": "chr1:1000-2000:+, chr1:3000-4000:-, chr1:9000-9999:+, notajunction",
    })
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_found"] == 2 and body["n_not_found"] == 1 and body["n_invalid"] == 1
    by_j = {r["junction"]: r for r in body["results"]}

    a = by_j["chr1:1000-2000:+"]
    assert a["found"] is True and a["error"] is None
    assert a["gene_id"] == "ENSG_A.1" and a["gene_name"] == "GENEA"
    assert a["width"] == 1000 and a["annotated"] is True
    assert a["left_motif"] == "GT" and a["right_motif"] == "AG"
    assert a["left_annotated"] == "0" and a["right_annotated"] == "0"

    b = by_j["chr1:3000-4000:-"]
    assert b["found"] is True
    assert b["gene_name"] == "GENEB" and b["annotated"] is False
    assert b["left_annotated"] == "tag1"

    missing = by_j["chr1:9000-9999:+"]
    assert missing["found"] is False and missing["error"] is None

    bad = by_j["notajunction"]
    assert bad["found"] is False and "chr:start-end:strand" in bad["error"]


def test_no_gene_overlap_row_has_null_gene_fields(client, sid, tmp_path):
    client.post(f"/api/session/{sid}/junction-metadata",
                files={"file": ("jm.csv", open(_lookup_csv(tmp_path), "rb"))})
    r = client.post(f"/api/session/{sid}/lookup", json={"junctions": "chr1:5000-6000:+"})
    body = r.json()["results"][0]
    assert body["found"] is True
    assert body["gene_id"] is None and body["gene_name"] is None


def test_reversed_start_end_is_normalised(client, sid, tmp_path):
    """A user typing coordinates in descending order still matches — the
    lookup key is built from parse_rowname()'s normalised (ascending) start/end,
    the same identity the index itself was built from."""
    client.post(f"/api/session/{sid}/junction-metadata",
                files={"file": ("jm.csv", open(_lookup_csv(tmp_path), "rb"))})
    r = client.post(f"/api/session/{sid}/lookup", json={"junctions": "chr1:2000-1000:+"})
    assert r.json()["results"][0]["found"] is True


def test_lookup_requires_index_loaded(client, sid):
    r = client.post(f"/api/session/{sid}/lookup", json={"junctions": "chr1:1-2:+"})
    assert r.status_code == 409


def test_lookup_requires_at_least_one_junction(client, sid, tmp_path):
    client.post(f"/api/session/{sid}/junction-metadata",
                files={"file": ("jm.csv", open(_lookup_csv(tmp_path), "rb"))})
    r = client.post(f"/api/session/{sid}/lookup", json={"junctions": "   "})
    assert r.status_code == 422


def test_missing_required_columns_errors(client, sid, tmp_path):
    csv = tmp_path / "bad.csv"
    csv.write_text("chrom,start,end\n1,2,3\n")
    r = client.post(f"/api/session/{sid}/junction-metadata", files={"file": ("bad.csv", open(csv, "rb"))})
    assert r.status_code == 422
    assert "missing column" in r.json()["detail"]


def test_dataload_routes_load_into_sjlookup(tmp_path, monkeypatch):
    """The dataload wiring (POST /api/dataload/sjlookup/{sid}/junction-metadata)
    mirrors the sjvc/sjsurv equivalents — load by path, not upload."""
    from fastapi.testclient import TestClient
    from main import app as combined_app

    monkeypatch.setenv("SJ_DATA_ROOT", str(tmp_path))
    client = TestClient(combined_app)
    sid = client.post("/api/sjlookup/session").json()["session_id"]
    r = client.post(f"/api/dataload/sjlookup/{sid}/junction-metadata",
                     json={"path": str(_lookup_csv(tmp_path))})
    assert r.status_code == 200, r.text
    assert r.json()["n_rows"] == 4

    st = client.get(f"/api/sjlookup/session/{sid}/state").json()
    assert st["has_index"] is True and st["n_rows"] == 4
