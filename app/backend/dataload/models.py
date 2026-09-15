"""Request / response schemas for the data-load endpoints."""
from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel


class PathBody(BaseModel):
    path: str


class GencodeApply(BaseModel):
    sjv_sid: Optional[str] = None
    sjvc_sid: Optional[str] = None
    sjsurv_sid: Optional[str] = None
    species: str = "human"
    release: str


class GtfApply(BaseModel):
    sjv_sid: Optional[str] = None
    sjvc_sid: Optional[str] = None
    sjsurv_sid: Optional[str] = None
    path: str


class GencodeApplied(BaseModel):
    label: str
    pathways_enabled: bool
    applied_to: List[str]
    warnings: List[str] = []


class PushGroupsBody(BaseModel):
    sjsurv_sid: str
    sjv_sid: Optional[str] = None
    sjvc_sid: Optional[str] = None


class PushGroupsResult(BaseModel):
    applied_to: List[str]           # "sjv" | "sjvc" — which sessions got updated
    n_sjsurv_samples: int           # samples in SJSurv's current Group/SurviverGroup
    n_group: int                    # of those, how many have a (non-null) Group
    n_labelled: int                 # of those, how many have a Good/Poor SurviverGroup
    matched: Dict[str, int] = {}    # per target: its own samples that got a value
    warnings: List[str] = []


class ScanFile(BaseModel):
    name: str
    size: int
    # gene_matrix | rrs_scores | clinical | clinical_extra | junction_counts
    #   | junction_metadata | gtf | unknown
    role: str
    target: Optional[str] = None       # which app it loads into, e.g. "sjvc"
    status: str                        # ok | caution | skip | unsupported
    note: str = ""


class ScanResponse(BaseModel):
    cohort: str
    dir: str
    files: List[ScanFile]
