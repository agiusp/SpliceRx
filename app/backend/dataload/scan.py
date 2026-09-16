"""Classify the files in a TCGA cohort directory.

Filename-based, against the naming the ``prepTCGAdata`` package writes
(``TCGA_<COHORT>_<kind>.rds``), plus a size guard: any ``.rds`` bigger than
``_BIG_RDS`` that isn't clearly a gene-level matrix is treated as the raw
junction count matrix — which the app cannot load in its current form. The big
files are never opened here.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import List

from .models import ScanFile

# A raw junction count matrix for a TCGA cohort is hundreds of MB (sparse) and
# ~8M rows; densified it does not fit in memory on a typical machine. A
# gene-level matrix is ~15-25 MB. Anything between is ambiguous -> treat large as
# junction-level.
_BIG_RDS = 250 * 1024 * 1024

_GENE_MATRIX = re.compile(r"per[_-]?gene|gene[_-]?level|by[_-]?gene|condensed", re.I)
_PATHWAY_MATRIX = re.compile(r"per[_-]?pathway|pathway[_-]?level|by[_-]?pathway", re.I)
_JUNCTION_META = re.compile(r"junction[_-]?(metadata|meta)\b", re.I)
_JUNCTION_COUNTS = re.compile(r"junction[_-]?counts?", re.I)
_RRS = re.compile(r"(^|[_-])rrs([_-]|\.|$)|rrs[_-]?scores?|relative[_-]?read[_-]?support", re.I)
_MSI = re.compile(r"(^|[_-])msi([_-]|\.|$)", re.I)
_SAMPLE_META = re.compile(r"sample[_-]?(metadata|meta|sheet)|clinical|phenotype|pheno", re.I)
_GTF = re.compile(r"\.gtf(\.gz)?$", re.I)


def classify_one(path: Path) -> ScanFile:
    name = path.name
    lower = name.lower()
    size = path.stat().st_size

    def mk(role: str, target, status: str, note: str) -> ScanFile:
        return ScanFile(name=name, size=size, role=role, target=target, status=status, note=note)

    if _GTF.search(lower):
        return mk(
            "gtf", "sjvc", "ok",
            "GENCODE GTF — optional; only needed for the 'protein-coding genes only' MAD filter",
        )

    if not lower.endswith(".rds"):
        return mk("unknown", None, "skip", "not a recognised data file")

    # order matters: the gene-level, pathway-level, and RRS file names usually
    # also contain "junction_counts", so check the more specific patterns first
    if _RRS.search(lower):
        return mk(
            "rrs_scores", "sjsurv", "ok",
            "RRS score matrix (feature x sample) — an sjdat choice for 2D View and SJSurv, and "
            "(with junction counts) the per-sample table on a single-junction SJ Lookup",
        )
    if _JUNCTION_META.search(lower):
        return mk(
            "junction_metadata", "sjvc", "ok",
            "per-junction gene annotation — lets 2D View look up typed genes / pathways "
            "against a junction-level matrix without a GENCODE release loaded",
        )
    if _GENE_MATRIX.search(lower):
        return mk(
            "gene_matrix", "sjvc", "ok",
            "gene-level count matrix — an sjdat choice for 2D View and SJSurv",
        )
    if _PATHWAY_MATRIX.search(lower):
        return mk(
            "pathway_matrix", "sjvc", "ok",
            "pathway-level count matrix — an sjdat choice for 2D View and SJSurv",
        )
    if _JUNCTION_COUNTS.search(lower) or size > _BIG_RDS:
        return mk(
            "junction_counts", "sjv", "ok",
            "junction-level count matrix — loads into the Sashimi plot, an sjdat choice for "
            "2D View and SJSurv, and (with RRS scores) the per-sample table on a single-"
            "junction SJ Lookup; the first load of a big sparse matrix takes ~1 min to index, "
            "then it's cached",
        )
    if _MSI.search(lower):
        return mk(
            "clinical_extra", "sjvc", "caution",
            "MSI table — its sample IDs are TCGA barcodes and may not match the matrix "
            "columns; load it only if you know the IDs line up",
        )
    if _SAMPLE_META.search(lower):
        return mk(
            "clinical", "sjvc", "ok",
            "sample / clinical table — the 2D View clinical table and the SJSurv sample "
            "metadata (needs the Group / SurviverGroup columns for SJSurv)",
        )
    return mk("unknown", None, "skip", "unrecognised .rds — not loaded")


def classify(directory: Path) -> List[ScanFile]:
    return [
        classify_one(p)
        for p in sorted(directory.iterdir())
        if p.is_file() and not p.name.startswith(".")
    ]
