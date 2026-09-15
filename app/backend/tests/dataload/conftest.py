import shutil
import sys
from pathlib import Path

import pytest

# app/backend/ — so `import main`, `import dataload`, `import sjvc` resolve.
BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

SJVC_FIXTURES = BACKEND / "tests" / "sjvc" / "fixtures"
SJV_FIXTURES = BACKEND / "tests" / "sjv" / "fixtures"
SJSURV_FIXTURES = BACKEND / "tests" / "sjsurv" / "fixtures"


@pytest.fixture
def cohort_dir(tmp_path, monkeypatch):
    """A fake TCGA cohort directory under an SJ_DATA_ROOT we control, with one
    file of every role the scanner recognises."""
    root = tmp_path / "Data"
    d = root / "TCGA_TEST"
    d.mkdir(parents=True)

    shutil.copy(SJVC_FIXTURES / "mini_genes.rds", d / "TCGA_TEST_novel_junction_counts_per_gene.rds")
    shutil.copy(SJVC_FIXTURES / "mini_clinical.rds", d / "TCGA_TEST_sample_metadata.rds")
    shutil.copy(SJVC_FIXTURES / "mini_clinical.rds", d / "TCGA_TEST_MSI.rds")
    (d / "TCGA_TEST_junction_metadata.rds").write_bytes(b"\x00" * 2048)
    # a real (tiny, sparse) junction count matrix for the SJV load path
    shutil.copy(SJV_FIXTURES / "mini_sparse.rds", d / "TCGA_TEST_junction_counts.rds")
    shutil.copy(SJV_FIXTURES / "mini.gtf", d / "TCGA_TEST.gtf")
    (d / "TCGA_TEST_sj_meta.csv").write_text(
        "sample_id,condition\nsample_A,tumor\nsample_B,tumor\nsample_C,normal\n"
    )
    (d / "notes.txt").write_text("ignore me")

    monkeypatch.setenv("SJ_DATA_ROOT", str(root))
    monkeypatch.setenv("SJV_CACHE_DIR", str(tmp_path / "sjvcache"))
    return d


@pytest.fixture
def sjsurv_cohort(tmp_path, monkeypatch):
    """A cohort directory whose sample metadata carries the Group / SurviverGroup
    columns SJSurv needs, plus an RRS-named matrix so the scanner tags one."""
    root = tmp_path / "Data"
    d = root / "TCGA_SURV"
    d.mkdir(parents=True)

    shutil.copy(SJSURV_FIXTURES / "mini_sjdat_genes.rds",
                d / "TCGA_SURV_novel_junction_counts_per_gene.rds")
    shutil.copy(SJSURV_FIXTURES / "mini_sjdat_genes.rds", d / "TCGA_SURV_RRS_scores.rds")
    shutil.copy(SJSURV_FIXTURES / "mini_metadata.csv", d / "TCGA_SURV_sample_metadata.csv")

    monkeypatch.setenv("SJ_DATA_ROOT", str(root))
    monkeypatch.setenv("SJV_CACHE_DIR", str(tmp_path / "sjvcache"))
    return d
