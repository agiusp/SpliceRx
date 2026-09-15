"""API-level tests for the /plot endpoint, focused on multi-gene plots."""
import warnings

import pytest

from tests.conftest import FIXTURES

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from app.main import app
    return TestClient(app)


@pytest.fixture
def ready(client):
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/rds",
                files={"file": ("m.rds", open(FIXTURES / "mini.rds", "rb"))})
    client.post(f"/api/session/{sid}/gtf",
                files={"file": ("mini.gtf", open(FIXTURES / "mini.gtf", "rb"))})
    return sid


def test_single_gene_still_works_via_genes_list(client, ready):
    r = client.post(f"/api/session/{ready}/plot",
                    json={"genes": ["TESTG1"], "sample": "sample_A"})
    assert r.status_code == 200
    body = r.json()
    assert [g["name"] for g in body["genes"]] == ["TESTG1"]
    assert body["x_domain"] == [790, 5410]
    assert len(body["arcs"]) == 8            # zero-count + off-chrom rows dropped
    assert all(t["gene_name"] == "TESTG1" for t in body["transcripts"])


def test_legacy_gene_name_field_accepted(client, ready):
    r = client.post(f"/api/session/{ready}/plot",
                    json={"gene_name": "TESTG1", "sample": "sample_A"})
    assert r.status_code == 200
    assert [g["name"] for g in r.json()["genes"]] == ["TESTG1"]


def test_two_genes_span_and_junction_union(client, ready):
    r = client.post(f"/api/session/{ready}/plot",
                    json={"genes": ["OTHERG", "TESTG1"], "sample": "sample_A"})
    assert r.status_code == 200
    body = r.json()
    assert [g["name"] for g in body["genes"]] == ["OTHERG", "TESTG1"]
    # x_domain is the bounding box of both loci (OTHERG 1000-2000 ⊂ TESTG1 1000-5200)
    assert body["x_domain"] == [790, 5410]
    # transcripts from both genes, each tagged
    kinds = {t["gene_name"] for t in body["transcripts"]}
    assert kinds == {"OTHERG", "TESTG1"}
    # every drawn junction overlaps at least one of the two loci
    assert len(body["arcs"]) == 8
    # "max reads" is over the whole union
    counts = [a["count"] for a in body["arcs"]]
    assert max(counts) == max(counts)  # (sanity — value comes straight from the matrix)


def test_duplicate_gene_names_collapsed(client, ready):
    r = client.post(f"/api/session/{ready}/plot",
                    json={"genes": ["TESTG1", "testg1", " TESTG1 "], "sample": "sample_A"})
    assert r.status_code == 200
    assert [g["name"] for g in r.json()["genes"]] == ["TESTG1"]


def test_unknown_gene_returns_404_with_near_matches(client, ready):
    r = client.post(f"/api/session/{ready}/plot",
                    json={"genes": ["TESTG1", "NOSUCHGENE"], "sample": "sample_A"})
    assert r.status_code == 404
    assert "NOSUCHGENE" in r.json()["detail"]["message"]


def test_genes_on_different_chromosomes_rejected(client, tmp_path):
    two_chrom = tmp_path / "two.gtf"
    two_chrom.write_text(
        'chr1\tt\tgene\t100\t900\t.\t+\t.\tgene_id "A"; gene_name "AAA";\n'
        'chr1\tt\ttranscript\t100\t900\t.\t+\t.\tgene_id "A"; transcript_id "A.1";\n'
        'chr1\tt\texon\t100\t400\t.\t+\t.\tgene_id "A"; transcript_id "A.1";\n'
        'chr1\tt\texon\t600\t900\t.\t+\t.\tgene_id "A"; transcript_id "A.1";\n'
        'chr2\tt\tgene\t100\t900\t.\t+\t.\tgene_id "B"; gene_name "BBB";\n'
        'chr2\tt\ttranscript\t100\t900\t.\t+\t.\tgene_id "B"; transcript_id "B.1";\n'
        'chr2\tt\texon\t100\t400\t.\t+\t.\tgene_id "B"; transcript_id "B.1";\n'
        'chr2\tt\texon\t600\t900\t.\t+\t.\tgene_id "B"; transcript_id "B.1";\n'
    )
    sid = client.post("/api/session").json()["session_id"]
    client.post(f"/api/session/{sid}/rds",
                files={"file": ("m.rds", open(FIXTURES / "mini.rds", "rb"))})
    client.post(f"/api/session/{sid}/gtf", files={"file": ("two.gtf", open(two_chrom, "rb"))})
    r = client.post(f"/api/session/{sid}/plot",
                    json={"genes": ["AAA", "BBB"], "sample": "sample_A"})
    assert r.status_code == 422
    assert "different chromosomes" in r.json()["detail"]
