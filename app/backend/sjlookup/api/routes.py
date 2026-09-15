"""HTTP API for SJ Lookup — mounted at /api/sjlookup."""
from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from fastapi import APIRouter, File, HTTPException, UploadFile

from ..models import (
    JunctionInfo,
    JunctionMetadataLoaded,
    LookupRequest,
    LookupResponse,
    SessionCreated,
    SessionState,
)
from ..services.lookup import load_junction_lookup_index, parse_junction_list
from ..services.sessions import store
from sjvc.services.junction_metadata import JunctionMetadataError
from sjvc.services.junctions import RownameError, parse_rowname

router = APIRouter()

_CAP = 2 * 1024 * 1024 * 1024


def _session(sid: str):
    s = store.get(sid)
    if s is None:
        raise HTTPException(404, "unknown or expired session")
    return s


async def _save_upload(file: UploadFile, dest: Path) -> None:
    size = 0
    with open(dest, "wb") as fh:
        while chunk := await file.read(1 << 20):
            size += len(chunk)
            if size > _CAP:
                fh.close()
                dest.unlink(missing_ok=True)
                raise HTTPException(413, "file exceeds the 2 GiB limit")
            fh.write(chunk)


def _clean(v):
    """Normalise a value pulled out of the lookup DataFrame for the response
    model: pandas/numpy missing markers -> None, numpy scalars -> native
    Python types."""
    if v is None:
        return None
    if pd.isna(v):
        return None
    if isinstance(v, np.generic):
        return v.item()
    return v


@router.get("/health")
def health() -> dict:
    return {"ok": True}


@router.post("/session", response_model=SessionCreated)
def create_session() -> SessionCreated:
    return SessionCreated(session_id=store.create().id)


def ingest_junction_metadata(s, src: Path) -> JunctionMetadataLoaded:
    try:
        idx = load_junction_lookup_index(src, s.tmp_dir)
    except JunctionMetadataError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)
    s.index = idx
    warnings: List[str] = []
    if idx.n_duplicate_rownames:
        warnings.append(
            f"{idx.n_duplicate_rownames} duplicate junction row(s) in the table were collapsed "
            f"(kept the first)"
        )
    return JunctionMetadataLoaded(
        n_rows=idx.n_rows, n_duplicate_rownames=idx.n_duplicate_rownames, warnings=warnings,
    )


@router.post("/session/{sid}/junction-metadata", response_model=JunctionMetadataLoaded)
async def upload_junction_metadata(sid: str, file: UploadFile = File(...)) -> JunctionMetadataLoaded:
    s = _session(sid)
    suffix = "".join(Path(file.filename or "junction_metadata.rds").suffixes) or ".rds"
    dest = s.tmp_dir / f"junction_metadata{suffix}"
    await _save_upload(file, dest)
    return ingest_junction_metadata(s, dest)


@router.get("/session/{sid}/state", response_model=SessionState)
def session_state(sid: str) -> SessionState:
    s = _session(sid)
    return SessionState(has_index=s.index is not None, n_rows=s.index.n_rows if s.index else 0)


@router.post("/session/{sid}/lookup", response_model=LookupResponse)
def lookup(sid: str, body: LookupRequest) -> LookupResponse:
    s = _session(sid)
    if s.index is None:
        raise HTTPException(409, "load the junction metadata table on the Data tab first")

    junctions = parse_junction_list(body.junctions)
    if not junctions:
        raise HTTPException(422, "no junction given")

    results: List[JunctionInfo] = []
    n_found = n_not_found = n_invalid = 0
    for j in junctions:
        try:
            parsed = parse_rowname(j)
        except RownameError as e:
            results.append(JunctionInfo(junction=j, found=False, error=str(e)))
            n_invalid += 1
            continue

        # look up by the *normalised* key (parse_rowname reorders a reversed
        # start/end into ascending order) — same identity the index itself
        # was built from
        key = f"{parsed.chrom}:{parsed.start}-{parsed.end}:{parsed.strand}"
        row = s.index.lookup(key)
        if row is None:
            results.append(JunctionInfo(junction=j, found=False))
            n_not_found += 1
            continue

        results.append(JunctionInfo(
            junction=j, found=True,
            gene_id=_clean(row.get("gene_id")), gene_name=_clean(row.get("gene_name")),
            width=_clean(row.get("width")), annotated=_clean(row.get("annotated")),
            left_motif=_clean(row.get("left_motif")), right_motif=_clean(row.get("right_motif")),
            left_annotated=_clean(row.get("left_annotated")), right_annotated=_clean(row.get("right_annotated")),
        ))
        n_found += 1

    return LookupResponse(results=results, n_found=n_found, n_not_found=n_not_found, n_invalid=n_invalid)
