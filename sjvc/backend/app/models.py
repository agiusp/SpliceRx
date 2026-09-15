"""Request/response schemas."""
from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel


class SessionCreated(BaseModel):
    session_id: str


class MatrixUploaded(BaseModel):
    samples: List[str]
    n_junctions: int
    feature_kind: str = "junction"   # "junction" | "gene" (already-condensed matrix)
    warnings: List[str] = []


class ClinicalColumn(BaseModel):
    name: str
    type: str            # "numeric" | "categorical"
    n_unique: int
    n_missing: int
    eligible: bool


class ClinicalUploaded(BaseModel):
    columns: List[ClinicalColumn]
    n_matched: int
    warnings: List[str] = []


class GencodeSelect(BaseModel):
    species: str
    release: str


class GencodeReady(BaseModel):
    label: str
    pathways_enabled: bool
    warnings: List[str] = []


class GeneSetRequest(BaseModel):
    mode: str                       # "typed" | "list" | "pathway"
    text: Optional[str] = None      # typed / list contents
    library: Optional[str] = None   # pathway
    term: Optional[str] = None      # pathway
    prefix: bool = False            # treat typed/list entries as name prefixes


class GeneSetResponse(BaseModel):
    matched: List[str]
    unmatched: List[str]
    n_genes: int
    source: str
    warnings: List[str] = []


class FeaturesRequest(BaseModel):
    condense: bool = False


class MadFeaturesRequest(BaseModel):
    condense: bool = False
    top_n: int = 50
    protein_coding_only: bool = False   # gene-level matrices: restrict to protein_coding biotype


class FeaturesResponse(BaseModel):
    feature_kind: str               # "junction" | "gene"
    n_features: int
    n_samples: int
    n_junctions: int
    dropped_zero_variance: int
    suggest_condense: bool
    feature_preview: List[str] = []  # first ~20 feature labels, for display
    warnings: List[str] = []


class ProjectionRequest(BaseModel):
    method: str                                   # "pca" | "umap"
    clinical: List[str] = []                      # 0..2 feature names
    overrides: Dict[str, str] = {}                # feature -> "numeric"|"categorical"
    pc_x: int = 1
    pc_y: int = 2
    n_neighbors: int = 15
    min_dist: float = 0.1


class HeatmapRequest(BaseModel):
    clinical: List[str] = []
    order: str = "cluster"                        # "cluster" | "group"
    group_by: Optional[str] = None
    row_zscore: bool = True
    overrides: Dict[str, str] = {}
