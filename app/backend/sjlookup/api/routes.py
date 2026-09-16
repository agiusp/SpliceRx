"""HTTP API for SJ Lookup — mounted at /api/sjlookup."""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
from fastapi import APIRouter, File, HTTPException, UploadFile

from ..models import (
    JunctionInfo,
    JunctionMetadataLoaded,
    LookupRequest,
    LookupResponse,
    SampleValue,
    SessionCreated,
    SessionState,
    SjdatLoaded,
)
from ..services.lookup import classify_unannotated, load_junction_lookup_index, parse_junction_list
from ..services.sessions import store
from sjvc.services.junction_metadata import JunctionMetadataError
from sjvc.services.junctions import RownameError, parse_rowname
from sjvc.services.sjdat import SjdatError, load_sjdat

router = APIRouter()

# the two junction-level sjdat matrices — the ones a "chr:start-end:strand"
# lookup key can actually address; gene_matrix/pathway_matrix are condensed
# past individual junctions, so SJ Lookup has no use for them
_SJDAT_KINDS = ("junction_counts", "rrs_scores")

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


def ingest_sjdat(s, kind: str, src: Path) -> SjdatLoaded:
    """Load the junction-counts or RRS-scores matrix at ``src`` into the
    session. ``src`` is consumed. Shared by the browser-upload endpoint below
    and the dataload package — mirrors sjvc/sjsurv's own ``ingest_sjdat``, but
    without their clinical-sample intersection step: SJ Lookup never narrows
    to a sample set, it just reports whatever samples the matrix has."""
    if kind not in _SJDAT_KINDS:
        src.unlink(missing_ok=True)
        raise HTTPException(422, f"unknown sjdat kind {kind!r} — SJ Lookup only uses junction_counts / rrs_scores")
    try:
        d = load_sjdat(kind, src, s.tmp_dir)
    except SjdatError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)
    s.sjdat[kind] = d
    return SjdatLoaded(kind=kind, n_features=d.n_features, n_samples=d.n_samples, sparse=d.sparse)


@router.post("/session/{sid}/sjdat/{kind}", response_model=SjdatLoaded)
async def upload_sjdat(sid: str, kind: str, file: UploadFile = File(...)) -> SjdatLoaded:
    s = _session(sid)
    dest = s.tmp_dir / f"sjdat_{kind}.rds"
    await _save_upload(file, dest)
    return ingest_sjdat(s, kind, dest)


@router.get("/session/{sid}/state", response_model=SessionState)
def session_state(sid: str) -> SessionState:
    s = _session(sid)
    return SessionState(
        has_index=s.index is not None, n_rows=s.index.n_rows if s.index else 0,
        sjdat_loaded=sorted(s.sjdat.keys()),
    )


def _sample_values(s, key: str) -> "tuple[Optional[List[SampleValue]], List[str]]":
    """Every sample's read count / RRS score for the single junction ``key``
    (``chr:start-end:strand``), sorted descending by RRS score (missing RRS
    scores last, then by sample id). ``None`` — not ``[]`` — when neither
    matrix is loaded at all, so the frontend can tell "nothing to show" apart
    from "loaded, but this junction isn't in either matrix"."""
    counts_d = s.sjdat.get("junction_counts")
    rrs_d = s.sjdat.get("rrs_scores")
    if counts_d is None and rrs_d is None:
        return None, []

    warnings: List[str] = []
    if counts_d is None:
        warnings.append("load the junction-counts matrix on the Data tab for per-sample read counts")
    if rrs_d is None:
        warnings.append("load the RRS-scores matrix on the Data tab for per-sample RRS scores")

    counts_row = counts_d._row.get(key) if counts_d is not None else None
    rrs_row = rrs_d._row.get(key) if rrs_d is not None else None
    if counts_d is not None and counts_row is None:
        warnings.append("this junction is not a row of the loaded junction-counts matrix")
    if rrs_d is not None and rrs_row is None:
        warnings.append("this junction is not a row of the loaded RRS-scores matrix")
    if counts_row is None and rrs_row is None:
        return [], warnings

    # sample universe: whichever matrix has this junction as a row; when both
    # do and their sample lists differ (a mismatched pair of files), the union
    if counts_row is not None and rrs_row is not None and counts_d.samples != rrs_d.samples:
        samples = sorted(set(counts_d.samples) | set(rrs_d.samples))
    else:
        samples = counts_d.samples if counts_row is not None else rrs_d.samples

    counts_vals = counts_d.dense_block([counts_row], list(range(counts_d.n_samples)))[0] if counts_row is not None else None
    rrs_vals = rrs_d.dense_block([rrs_row], list(range(rrs_d.n_samples)))[0] if rrs_row is not None else None
    counts_col = counts_d._col if counts_d is not None else {}
    rrs_col = rrs_d._col if rrs_d is not None else {}

    def _val(vals, col, sample):
        if vals is None or sample not in col:
            return None
        v = float(vals[col[sample]])
        return None if np.isnan(v) else v

    out = [
        SampleValue(sample=sample, count=_val(counts_vals, counts_col, sample), rrs_score=_val(rrs_vals, rrs_col, sample))
        for sample in samples
    ]
    out.sort(key=lambda sv: (sv.rrs_score is None, -(sv.rrs_score or 0.0), sv.sample))
    return out, warnings


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
    single_key = None  # set only when exactly one (valid) junction was queried
    for j in junctions:
        try:
            parsed = parse_rowname(j)
        except RownameError as e:
            results.append(JunctionInfo(junction=j, found=False, error=str(e)))
            n_invalid += 1
            continue

        # look up by the *normalised* key (parse_rowname reorders a reversed
        # start/end into ascending order) — same identity the index itself
        # was built from, and the sjdat matrices below use the same convention
        key = f"{parsed.chrom}:{parsed.start}-{parsed.end}:{parsed.strand}"
        if len(junctions) == 1:
            single_key = key
        row = s.index.lookup(key)
        if row is None:
            results.append(JunctionInfo(junction=j, found=False))
            n_not_found += 1
            continue

        results.append(JunctionInfo(
            junction=j, found=True,
            gene_id=_clean(row.get("gene_id")), gene_name=_clean(row.get("gene_name")),
            width=_clean(row.get("width")), annotated=_clean(row.get("annotated")),
            category=classify_unannotated(s.index, row),
            left_motif=_clean(row.get("left_motif")), right_motif=_clean(row.get("right_motif")),
            left_annotated=_clean(row.get("left_annotated")), right_annotated=_clean(row.get("right_annotated")),
        ))
        n_found += 1

    sample_values, sample_values_warnings = (None, [])
    if single_key is not None:
        sample_values, sample_values_warnings = _sample_values(s, single_key)

    return LookupResponse(
        results=results, n_found=n_found, n_not_found=n_not_found, n_invalid=n_invalid,
        sample_values=sample_values, sample_values_warnings=sample_values_warnings,
    )
