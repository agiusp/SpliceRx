"""Request / response schemas for the SJSurv API."""
from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel

# Gene-set feature selection reuses sjvc's request/response shapes verbatim —
# SJSurv's active sjdat matrix is a `sjvc.services.sjdat.Sjdat`, the same type
# 2D View resolves gene sets against, so there's nothing SJSurv-specific to
# say in these shapes. See sjsurv/api/routes.py for the endpoints.
from sjvc.models import (  # noqa: F401 (re-exported)
    FeaturesRequest,
    FeaturesResponse,
    GeneSetRequest,
    GeneSetResponse,
    JunctionMetadataLoaded,
)


class SessionCreated(BaseModel):
    session_id: str


class SjdatLoaded(BaseModel):
    kind: str
    n_features: int
    n_samples: int
    sparse: bool = False
    warnings: List[str] = []


class AgeBandsOut(BaseModel):
    edges: List[Optional[float]]   # ascending; null = open-ended (only meaningful as the last edge)
    include_lowest: bool
    labels: List[str]


class HistologyCountOut(BaseModel):
    value: str            # a raw histology value
    n: int                # samples carrying it


class MetadataLoaded(BaseModel):
    n_matched: int
    n_unmatched: int
    n_with_age: int = 0
    age_min: Optional[float] = None
    age_max: Optional[float] = None
    age_bands: Optional[AgeBandsOut] = None
    min_group_n: Optional[int] = None
    event_quantile: float = 0.5
    histology_counts: List[HistologyCountOut] = []   # raw values + counts, for the merge editor
    use_histology: bool = True                       # whether Group currently includes histology
    histology_map: Dict[str, str] = {}                # raw value -> merged label, currently applied
    # whether a recognisable vital-status column (Alive/Dead) was found — a
    # group's MedianSurvival is only actually censoring-adjusted when this is
    # True; without it, every sample is conservatively treated as a confirmed
    # death (see services/metadata.py's RawMetadata.has_vital_status)
    has_vital_status: bool = False
    warnings: List[str] = []


class CovariateColumnOut(BaseModel):
    key: str              # "__histology__" | "__stage__" | "__age__", or a raw column name
    label: str
    kind: str              # "numeric" | "categorical"
    n_available: int
    default: bool


class CovariatesRequest(BaseModel):
    columns: List[str]     # covariate keys to use alongside the selected molecular features


class SuggestAgeBandsRequest(BaseModel):
    n_bands: int = 3


class StratifyRequest(BaseModel):
    edges: List[Optional[float]]   # ascending; null as the last edge = open-ended (>edges[-2])
    include_lowest: bool = False
    min_group_n: int = 20
    use_histology: bool = True     # False drops Histology from Group entirely
    histology_map: Dict[str, str] = {}   # raw value -> merged label; values not listed pass through
    # fraction of a group required to have died (an "event") for its
    # Cox-estimated survival-time threshold to be defined; 0.5 is the classic
    # median. Lower it (e.g. 0.25) for a low-mortality cohort where the 50%
    # median is never reached within follow-up for most groups — the
    # threshold time is then reached sooner, at the cost of a more lopsided
    # Good/Poor split (few, especially-early deaths vs. everyone else).
    event_quantile: float = 0.5


class GroupCountOut(BaseModel):
    group: str            # the raw Group value, or "__all__"
    label: str            # display label
    n_total: int
    n_good: int
    n_poor: int
    n_labelled: int


class GroupsResponse(BaseModel):
    groups: List[GroupCountOut]
    warnings: List[str] = []


class SjdatOption(BaseModel):
    kind: str
    label: str
    description: str
    loaded: bool
    n_features: int = 0
    n_samples: int = 0
    sparse: bool = False
    # rows with at least one non-zero entry — only populated for a small
    # dense matrix (pathway_matrix today) where it's cheap; None otherwise
    n_nonzero_rows: Optional[int] = None


class SessionState(BaseModel):
    """Everything the frontend needs to rebuild the UI after a Data-tab hand-off
    or a browser refresh."""
    has_raw_metadata: bool = False
    has_metadata: bool = False
    metadata: Optional[MetadataLoaded] = None
    sjdat_options: List[SjdatOption] = []
    active_sjdat: Optional[str] = None
    # gene-set feature selection — mirrors sjvc.models.SessionState's fields
    gencode_label: Optional[str] = None
    pathways_enabled: bool = False
    has_junction_metadata: bool = False
    has_geneset: bool = False
    has_features: bool = False
    selected_group: Optional[str] = None
    has_selection: bool = False
    has_model: bool = False
    # clinical covariates (from the sample-metadata file) offered alongside
    # the selected molecular features — see services/metadata.py's
    # available_covariates()/build_covariates()
    covariate_columns: List[CovariateColumnOut] = []
    selected_covariates: List[str] = []


class ActivateSjdat(BaseModel):
    kind: str             # "junction_counts" | "rrs_scores" | "gene_matrix"


class SelectRequest(BaseModel):
    group: str                     # a Group value, or "__all__"
    # N / X are the coverage prefilter — only meaningful (and only shown by
    # the frontend) for a junction-level sjdat (junction_counts / rrs_scores);
    # omitted/None skips it entirely and ranks every row by MAD.
    n_min: Optional[float] = None  # N — count if >= 1, else fraction of the group
    x_min: Optional[float] = 0.0   # X — minimum magnitude of a counted entry
    top_n: int = 100               # n — keep this many by MAD
    protein_coding_only: bool = False   # restrict to protein_coding biotype, excluding MT-* genes
    # exclude any gene in the curated low-mappability paralog-family list
    # (HLA, immunoglobulin/TCR loci, MT-*, olfactory receptors, and other
    # named segmental-duplication clusters), independently of
    # protein_coding_only (either, both, or neither may be set)
    exclude_paralog_families: bool = False
    # RRS-scores-only: keep a junction only when its corresponding row in the
    # (separately loaded) junction_counts sjdat has a max supporting read
    # count above this, across the group's samples; None skips it
    min_supporting_reads: Optional[float] = None


class SelectGenesetRequest(BaseModel):
    """Turn an already-built gene-set FeatureMatrix (`POST .../features`) into
    a Selection for the chosen Group — the gene-set counterpart of
    `SelectRequest`/`POST .../select`, minus the MAD-specific fields since the
    resolved gene set *is* the selection, unranked."""
    group: str                     # a Group value, or "__all__"


class SelectResponse(BaseModel):
    group: str
    sjdat_kind: str
    n_group_samples: int           # labelled samples in the group that are in the matrix
    n_good: int
    n_poor: int
    n_candidates: int              # rows in sjdat
    n_after_coverage: int          # passed the N / X filter
    n_selected: int                # kept after MAD top-n
    feature_preview: List[str] = []
    warnings: List[str] = []


class CVRequest(BaseModel):
    n_cv: Optional[int] = None     # folds; defaults to 5 if not given


class CVResponse(BaseModel):
    n_samples: int
    n_good: int
    n_poor: int
    n_features: int
    n_splits: int
    auc: float
    fold_aucs: List[float]
    fold_auc_mean: float
    fold_auc_sd: float
    accuracy: float
    sensitivity: float
    specificity: float
    confusion: List[List[int]]     # [[TN, FP], [FN, TP]]
    baseline_accuracy: float
    messages: List[str]


class FeatureWeightOut(BaseModel):
    feature: str
    weight: float
    abs_weight: float
    direction: str
    mean_good: float
    mean_poor: float


class ModelResponse(BaseModel):
    sjdat_kind: str
    group: str
    n_samples: int
    n_good: int
    n_poor: int
    n_features: int
    auc_resub: float
    accuracy_resub: float
    confusion_resub: List[List[int]]
    features: List[FeatureWeightOut]
    saved: bool = True
    messages: List[str] = []
    warnings: List[str] = []
