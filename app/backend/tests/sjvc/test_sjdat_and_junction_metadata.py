"""NSJCG's ability to pick any of the 4 sjdat matrices (like SJSurv), and the
fast junction-metadata-based gene lookup for a junction-level one."""
import warnings

import numpy as np
import pytest

from tests.sjvc.conftest import FIXTURES

from sjvc.services.junction_metadata import JunctionMetadataError, load_junction_gene_index
from sjvc.services.mad import _mad, _mad_and_var_sparse_rows
from sjvc.services.sjdat import Sjdat

warnings.filterwarnings("ignore")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from sjvc.main import app
    return TestClient(app)


@pytest.fixture
def sid(client):
    return client.post("/api/session").json()["session_id"]


# --------------------------------------------------------------------------- #
# multiple sjdat matrices, like SJSurv
# --------------------------------------------------------------------------- #
def test_load_two_sjdat_kinds_and_switch_between_them(client, sid):
    r = client.post(f"/api/session/{sid}/sjdat/junction_counts",
                     files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "junction_counts"

    r = client.post(f"/api/session/{sid}/sjdat/gene_matrix",
                     files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "gene_matrix"

    st = client.get(f"/api/session/{sid}/state").json()
    opts = {o["kind"]: o for o in st["sjdat_options"]}
    assert opts["junction_counts"]["loaded"] and opts["gene_matrix"]["loaded"]
    assert not opts["rrs_scores"]["loaded"]
    assert st["active_sjdat"] == "junction_counts"          # first one loaded stays active
    assert st["feature_kind"] == "junction"

    r = client.post(f"/api/session/{sid}/sjdat", json={"kind": "gene_matrix"})
    assert r.status_code == 200, r.text
    assert r.json()["active_sjdat"] == "gene_matrix"
    assert r.json()["feature_kind"] == "gene"


def test_pathway_matrix_kind_loads_and_ranks_like_gene_matrix(client, sid):
    """pathway_matrix is dense with non-chr:start-end:strand row names —
    structurally identical to gene_matrix — so it goes through the exact
    same dense reader and MAD/feature-building paths with no special-casing."""
    r = client.post(f"/api/session/{sid}/sjdat/pathway_matrix",
                     files={"file": ("p.rds", open(FIXTURES / "mini_genes.rds", "rb"))})
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "pathway_matrix" and not r.json()["sparse"]

    st = client.get(f"/api/session/{sid}/state").json()
    opt = next(o for o in st["sjdat_options"] if o["kind"] == "pathway_matrix")
    assert opt["loaded"] and opt["n_features"] == 3
    assert st["feature_kind"] == "gene"  # non-junction row names -> treated as gene-level

    r = client.post(f"/api/session/{sid}/features/mad", json={"top_n": 2})
    assert r.status_code == 200, r.text
    assert r.json()["n_features"] == 2


def test_count_nonzero_rows_dense_and_sparse():
    """A row of all-zero (or all-NaN) entries doesn't count; a row with at
    least one real non-zero value does — this is what the "top by MAD"
    default for pathway_matrix is seeded from."""
    dense = Sjdat(
        kind="pathway_matrix",
        features=["allzero", "hasone", "allnan", "mixed"],
        samples=["s1", "s2", "s3"],
        values=np.array([
            [0.0, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [np.nan, np.nan, np.nan],
            [np.nan, 0.0, 2.0],
        ]),
        sparse=False,
    )
    assert dense.count_nonzero_rows() == 2   # "hasone" and "mixed"

    import scipy.sparse as sp
    sparse = Sjdat(
        kind="junction_counts",
        features=["allzero", "hasone"],
        samples=["s1", "s2", "s3"],
        values=sp.csc_matrix(np.array([[0.0, 0.0, 0.0], [0.0, 5.0, 0.0]])),
        sparse=True,
    )
    assert sparse.count_nonzero_rows() == 1


def test_sjdat_options_surface_n_nonzero_rows_only_for_pathway_matrix(client, sid):
    with open(FIXTURES / "mini_genes.rds", "rb") as fh:
        client.post(f"/api/session/{sid}/sjdat/gene_matrix", files={"file": fh})
    with open(FIXTURES / "mini_genes.rds", "rb") as fh:
        client.post(f"/api/session/{sid}/sjdat/pathway_matrix", files={"file": fh})

    st = client.get(f"/api/session/{sid}/state").json()
    opts = {o["kind"]: o for o in st["sjdat_options"]}
    assert opts["gene_matrix"]["n_nonzero_rows"] is None            # not computed — not pathway
    assert opts["pathway_matrix"]["n_nonzero_rows"] == 3            # all 3 rows have a nonzero entry
    assert opts["junction_counts"]["n_nonzero_rows"] is None        # not loaded


def test_activating_unloaded_kind_errors(client, sid):
    client.post(f"/api/session/{sid}/sjdat/gene_matrix",
                files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})
    r = client.post(f"/api/session/{sid}/sjdat", json={"kind": "rrs_scores"})
    assert r.status_code == 409


def test_old_junctions_endpoint_still_works_as_gene_matrix(client, sid):
    """Back-compat: the original single-matrix endpoint keeps working exactly
    as before, now implemented as ingest_sjdat(..., "gene_matrix", ...)."""
    r = client.post(f"/api/session/{sid}/junctions",
                     files={"file": ("g.rds", open(FIXTURES / "mini_genes.rds", "rb"))})
    assert r.status_code == 200, r.text
    st = client.get(f"/api/session/{sid}/state").json()
    assert st["active_sjdat"] == "gene_matrix"
    assert next(o for o in st["sjdat_options"] if o["kind"] == "gene_matrix")["loaded"]


# --------------------------------------------------------------------------- #
# junction-metadata fast gene lookup
# --------------------------------------------------------------------------- #
def _junction_metadata_csv(tmp_path):
    """Gene annotation for a few of mini_junctions.rds's real rows, plus one
    multi-gene overlap and one row with a name/id list-length mismatch (a
    data issue that must be skipped, not crash the loader)."""
    csv = tmp_path / "junction_metadata.csv"
    csv.write_text(
        "seqnames,start,end,strand,gencode_gene_id,gencode_gene_name\n"
        "chr1,1201,1999,+,ENSG_A.1,GENEA\n"
        "chr1,2201,2999,+,ENSG_A.1,GENEA\n"
        "chr1,3201,3999,+,\"ENSG_A.1,ENSG_B.1\",\"GENEA,GENEB\"\n"
        "chr1,4201,4999,+,ENSG_B.1,GENEB\n"
        "chr1,1500,2800,+,\"ENSG_C.1,ENSG_C.2\",GENEC\n"     # id/name length mismatch -> skipped
        "chr2,5000,6000,+,NA,NA\n"                            # no gene -> skipped
    )
    return csv


def test_junction_metadata_index_resolve(tmp_path):
    idx = load_junction_gene_index(_junction_metadata_csv(tmp_path), tmp_path)
    assert idx.n_rows == 6
    jgmap, labels, unmatched = idx.resolve(["GENEA", "nope"])
    assert unmatched == ["nope"]
    assert labels == {"ENSG_A.1": "GENEA"}
    assert set(jgmap.rownames) == {"chr1:1201-1999:+", "chr1:2201-2999:+", "chr1:3201-3999:+"}
    # the mismatched-length row never made it into the index at all
    assert "chr1:1500-2800:+" not in idx.rownames_by_gene_id.get("ENSG_C.1", [])


def test_junction_metadata_requires_expected_columns(tmp_path):
    csv = tmp_path / "bad.csv"
    csv.write_text("chrom,start,end\n1,2,3\n")
    with pytest.raises(JunctionMetadataError, match="missing column"):
        load_junction_gene_index(csv, tmp_path)


def test_geneset_via_junction_metadata_needs_no_gencode(client, sid, tmp_path):
    client.post(f"/api/session/{sid}/sjdat/junction_counts",
                files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    r = client.post(f"/api/session/{sid}/junction-metadata",
                     files={"file": ("jm.csv", open(_junction_metadata_csv(tmp_path), "rb"))})
    assert r.status_code == 200, r.text
    assert r.json() == {"n_rows": 6, "n_genes": 2}

    # no GENCODE / GTF loaded at all — the junction-metadata lookup is enough
    r = client.post(f"/api/session/{sid}/geneset", json={"mode": "typed", "text": "GENEA, GENEB"})
    assert r.status_code == 200, r.text
    body = r.json()
    assert set(body["matched"]) == {"GENEA", "GENEB"}
    assert body["unmatched"] == []

    r = client.post(f"/api/session/{sid}/features", json={"condense": False})
    assert r.status_code == 200, r.text
    fbody = r.json()
    assert fbody["feature_kind"] == "junction"
    assert fbody["n_features"] == 4   # the 4 junctions overlapping GENEA and/or GENEB

    r = client.post(f"/api/session/{sid}/features", json={"condense": True})
    assert r.status_code == 200, r.text
    assert r.json()["feature_kind"] == "gene"
    assert r.json()["n_features"] == 2   # GENEA, GENEB


def test_geneset_without_gencode_or_junction_metadata_errors(client, sid):
    client.post(f"/api/session/{sid}/sjdat/junction_counts",
                files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    r = client.post(f"/api/session/{sid}/geneset", json={"mode": "typed", "text": "GENEA"})
    assert r.status_code == 409
    assert "junction metadata" in r.json()["detail"]


# --------------------------------------------------------------------------- #
# sparse-safe "Top by MAD" ranking
# --------------------------------------------------------------------------- #
def test_sparse_row_mad_matches_dense():
    rng = np.random.default_rng(0)
    dense = rng.poisson(2, size=(30, 12)).astype(float)
    dense[dense < 1] = 0
    from scipy import sparse
    sp = sparse.csc_matrix(dense)

    y = np.log1p(dense)
    dense_scores, dense_vars = _mad(y), y.var(axis=1)
    sparse_scores, sparse_vars = _mad_and_var_sparse_rows(sp, list(range(30)), 12)
    assert np.allclose(dense_scores, sparse_scores)
    assert np.allclose(dense_vars, sparse_vars)


def test_mad_route_ranks_a_sparse_junction_matrix(client, sid):
    r = client.post(f"/api/session/{sid}/sjdat/junction_counts",
                     files={"file": ("j.rds", open(FIXTURES / "mini_junctions.rds", "rb"))})
    assert r.status_code == 200, r.text
    r = client.post(f"/api/session/{sid}/features/mad", json={"top_n": 3})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["feature_kind"] == "junction"
    assert body["n_features"] == 3
