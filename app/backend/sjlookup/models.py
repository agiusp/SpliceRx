"""Request / response schemas for the SJ Lookup API."""
from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel


class SessionCreated(BaseModel):
    session_id: str


class JunctionMetadataLoaded(BaseModel):
    n_rows: int
    n_duplicate_rownames: int
    warnings: List[str] = []


class SjdatLoaded(BaseModel):
    kind: str
    n_features: int
    n_samples: int
    sparse: bool = False
    warnings: List[str] = []


class SessionState(BaseModel):
    has_index: bool = False
    n_rows: int = 0
    sjdat_loaded: List[str] = []   # which of junction_counts / rrs_scores are loaded


class LookupRequest(BaseModel):
    junctions: str   # free text — one or more, comma / space / newline separated


class JunctionInfo(BaseModel):
    junction: str                       # echoes the queried string, as given
    found: bool
    error: Optional[str] = None         # set instead of the fields below when malformed
    gene_id: Optional[str] = None
    gene_name: Optional[str] = None
    width: Optional[int] = None
    annotated: Optional[bool] = None
    # set only when annotated is False — the same category terminology the
    # Sashimi plot uses (see sjv.services.classify.CATEGORIES), derived from
    # this row's own left_annotated/right_annotated/strand rather than a live
    # GENCODE transcript lookup (SJ Lookup carries no GENCODE annotation) —
    # see services.lookup.classify_unannotated for the exact rules and the
    # one case (exon_skipping vs. isoform_switch) it can't tell apart.
    category: Optional[str] = None
    left_motif: Optional[str] = None
    right_motif: Optional[str] = None
    left_annotated: Optional[str] = None
    right_annotated: Optional[str] = None


class SampleValue(BaseModel):
    sample: str
    count: Optional[float] = None
    rrs_score: Optional[float] = None


class LookupResponse(BaseModel):
    results: List[JunctionInfo]
    n_found: int
    n_not_found: int
    n_invalid: int
    # populated only when exactly one junction was queried and found — every
    # sample's read count / RRS score for that junction, sorted descending by
    # RRS score (nulls last). None when a matrix isn't loaded at all.
    sample_values: Optional[List[SampleValue]] = None
    sample_values_warnings: List[str] = []
