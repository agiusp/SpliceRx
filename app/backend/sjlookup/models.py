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


class SessionState(BaseModel):
    has_index: bool = False
    n_rows: int = 0


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
    left_motif: Optional[str] = None
    right_motif: Optional[str] = None
    left_annotated: Optional[str] = None
    right_annotated: Optional[str] = None


class LookupResponse(BaseModel):
    results: List[JunctionInfo]
    n_found: int
    n_not_found: int
    n_invalid: int
