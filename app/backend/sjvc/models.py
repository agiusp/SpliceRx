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


class SjdatOption(BaseModel):
    kind: str                        # "junction_counts" | "rrs_scores" | "gene_matrix"
    label: str
    description: str
    loaded: bool
    n_features: int = 0
    n_samples: int = 0
    sparse: bool = False
    # rows with at least one non-zero entry — only populated for a small
    # dense matrix (pathway_matrix today) where it's cheap; None otherwise
    n_nonzero_rows: Optional[int] = None


class SjdatLoaded(BaseModel):
    kind: str
    n_features: int
    n_samples: int
    sparse: bool = False
    warnings: List[str] = []


class ActivateSjdat(BaseModel):
    kind: str


class JunctionMetadataLoaded(BaseModel):
    n_rows: int
    n_genes: int


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


class SessionState(BaseModel):
    """Reconstructed step-1 state of a session — used by the Data tab hand-off
    and to survive a browser refresh."""
    has_junctions: bool
    feature_kind: Optional[str] = None       # "junction" | "gene"
    samples: List[str] = []
    n_junctions: int = 0
    clinical: Optional[ClinicalUploaded] = None
    gencode_label: Optional[str] = None
    pathways_enabled: bool = False
    sjdat_options: List[SjdatOption] = []
    active_sjdat: Optional[str] = None
    has_junction_metadata: bool = False


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
    # optionally narrow the resolved gene set down to its top_n most variable
    # rows by MAD — see mad.refine_by_mad(); None (the default) keeps every
    # resolved feature, today's behavior
    mad_top_n: Optional[int] = None


class MadFeaturesRequest(BaseModel):
    condense: bool = False
    top_n: int = 50
    protein_coding_only: bool = False   # restrict to protein_coding biotype, excluding MT-* genes
    # exclude any gene in the curated low-mappability paralog-family list
    # (HLA, immunoglobulin/TCR loci, MT-*, olfactory receptors, and other
    # named segmental-duplication clusters — see gencode.PARALOG_FAMILY_PATTERNS),
    # independently of protein_coding_only (either, both, or neither may be set)
    exclude_paralog_families: bool = False
    # coverage prefilter — meaningful (and only shown by the frontend) for a
    # junction-level matrix (junction_counts / rrs_scores); None skips it
    n_min: Optional[float] = None       # N — count if >= 1, else fraction of all samples
    x_min: float = 0.0                  # X — minimum magnitude of a counted entry
    # RRS-scores-only: keep a junction only when its corresponding row in the
    # (separately loaded) junction_counts sjdat has a max supporting read
    # count above this, across the matrix's own samples; None skips it
    min_supporting_reads: Optional[float] = None


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
    # drop samples missing any of `clinical`'s features before clustering,
    # rather than showing them with a grey annotation cell — unlike the
    # projection's equivalent toggle, this needs the heatmap rebuilt, since
    # removing a column changes the clustering itself.
    drop_missing_clinical: bool = False
