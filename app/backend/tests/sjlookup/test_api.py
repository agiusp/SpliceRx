"""SJ Lookup — direct per-junction lookup against a junction metadata table."""
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore")

_SJV_FIXTURES = Path(__file__).resolve().parents[1] / "sjv" / "fixtures"


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
    assert st == {"has_index": True, "n_rows": 4, "sjdat_loaded": []}

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
    assert a["category"] is None  # annotated junctions carry no category

    b = by_j["chr1:3000-4000:-"]
    assert b["found"] is True
    assert b["gene_name"] == "GENEB" and b["annotated"] is False
    assert b["left_annotated"] == "tag1"
    # "-" strand: donor = end (right, unannotated), acceptor = start (left,
    # annotated) -> exactly one splice site known -> the alt category names
    # the *other* (alternative) site, which here is the donor -> alt_5p
    assert b["category"] == "alt_5p"

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


def test_dataload_routes_load_sjdat_into_sjlookup(tmp_path, monkeypatch):
    """POST /api/dataload/sjlookup/{sid}/sjdat/{kind} — the per-sample-table
    matrices, loaded by path like every other dataload route."""
    import shutil
    from fastapi.testclient import TestClient
    from main import app as combined_app

    monkeypatch.setenv("SJ_DATA_ROOT", str(tmp_path))
    mini = tmp_path / "mini.rds"
    shutil.copy(_SJV_FIXTURES / "mini.rds", mini)

    client = TestClient(combined_app)
    sid = client.post("/api/sjlookup/session").json()["session_id"]
    r = client.post(f"/api/dataload/sjlookup/{sid}/sjdat/junction_counts", json={"path": str(mini)})
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "junction_counts" and r.json()["n_samples"] == 3

    st = client.get(f"/api/sjlookup/session/{sid}/state").json()
    assert st["sjdat_loaded"] == ["junction_counts"]

    r = client.post(f"/api/dataload/sjlookup/{sid}/sjdat/gene_matrix", json={"path": str(mini)})
    assert r.status_code == 422  # not a junction-level kind SJ Lookup uses


def _classify_csv(tmp_path):
    """One annotated intron (ENSG_C.1) plus four unannotated junctions, one
    of each non-annotated category classify_unannotated() can reach without
    a live GENCODE annotation — see services.lookup.classify_unannotated."""
    csv = tmp_path / "classify.csv"
    csv.write_text(
        "seqnames,start,end,strand,width,annotated,left_motif,right_motif,"
        "left_annotated,right_annotated,gencode_gene_id,gencode_gene_name\n"
        # the annotated intron novel_exon nests inside
        "chr2,1000,5000,+,4000,1,GT,AG,rG19,rG19,ENSG_C.1,GENEC\n"
        # neither site annotated, but [2000,2500] sits strictly inside the
        # annotated [1000,5000] intron of the same gene -> novel_exon
        "chr2,2000,2500,+,500,0,GT,AG,0,0,ENSG_C.1,GENEC\n"
        # both sites individually annotated elsewhere, just not as this
        # intron -> isoform_switch (this table can't tell exon_skipping apart)
        "chr2,7000,7500,+,500,0,GT,AG,siteA,siteB,ENSG_D.1,GENED\n"
        # only the donor (left, + strand) site annotated -> alt_3p
        "chr2,8000,8500,+,500,0,GT,AG,siteA,0,ENSG_D.1,GENED\n"
        # neither site annotated, no overlapping gene to check nesting against
        "chr2,9000,9500,+,500,0,GT,AG,0,0,NA,NA\n"
    )
    return csv


def test_classify_unannotated_categories(client, sid, tmp_path):
    client.post(f"/api/session/{sid}/junction-metadata",
                files={"file": ("jm.csv", open(_classify_csv(tmp_path), "rb"))})
    r = client.post(f"/api/session/{sid}/lookup", json={"junctions": (
        "chr2:1000-5000:+, chr2:2000-2500:+, chr2:7000-7500:+, "
        "chr2:8000-8500:+, chr2:9000-9500:+"
    )})
    by_j = {row["junction"]: row for row in r.json()["results"]}
    assert by_j["chr2:1000-5000:+"]["category"] is None       # annotated
    assert by_j["chr2:2000-2500:+"]["category"] == "novel_exon"
    assert by_j["chr2:7000-7500:+"]["category"] == "isoform_switch"
    assert by_j["chr2:8000-8500:+"]["category"] == "alt_3p"
    assert by_j["chr2:9000-9500:+"]["category"] == "novel"


def test_sample_values_for_single_junction(client, sid, tmp_path):
    """Querying exactly one junction also returns every sample's read count /
    RRS score for it, sorted descending by RRS score."""
    client.post(f"/api/session/{sid}/junction-metadata",
                files={"file": ("jm.csv", open(_lookup_csv(tmp_path), "rb"))})

    # no sjdat matrix loaded yet -> sample_values is absent (None), not empty
    r = client.post(f"/api/session/{sid}/lookup", json={"junctions": "chr1:1000-2000:+"})
    assert r.json()["sample_values"] is None

    mini = _SJV_FIXTURES / "mini.rds"
    for kind in ("junction_counts", "rrs_scores"):
        r = client.post(f"/api/session/{sid}/sjdat/{kind}", files={"file": (f"{kind}.rds", open(mini, "rb"))})
        assert r.status_code == 200, r.text

    st = client.get(f"/api/session/{sid}/state").json()
    assert sorted(st["sjdat_loaded"]) == ["junction_counts", "rrs_scores"]

    # chr1:1201-1999:+ in mini.rds: sample_A=50, sample_B=36, sample_C=31
    r = client.post(f"/api/session/{sid}/lookup", json={"junctions": "chr1:1201-1999:+"})
    body = r.json()
    assert body["sample_values_warnings"] == []
    got = [(sv["sample"], sv["count"], sv["rrs_score"]) for sv in body["sample_values"]]
    assert got == [
        ("sample_A", 50.0, 50.0),
        ("sample_B", 36.0, 36.0),
        ("sample_C", 31.0, 31.0),
    ]

    # querying more than one junction at once doesn't compute sample_values
    r = client.post(f"/api/session/{sid}/lookup", json={
        "junctions": "chr1:1201-1999:+, chr1:2201-2999:+",
    })
    assert r.json()["sample_values"] is None
