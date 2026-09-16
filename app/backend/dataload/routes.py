"""Data-load endpoints — mounted at /api/dataload.

- POST /scan                       classify a cohort directory
- POST /gencode                    pick a GENCODE release for every live app session at once
- POST /gtf                        use a GTF file for every live app session at once
- POST /sjvc/{sid}/gene-matrix       load a gene-level matrix into an SJVC session (legacy —
                                     prefer /sjvc/{sid}/sjdat/gene_matrix below)
- POST /sjvc/{sid}/clinical          load a clinical table into an SJVC session
- POST /sjvc/{sid}/sjdat/{kind}      load an sjdat matrix into an SJVC (2D View) session
- POST /sjvc/{sid}/junction-metadata load the per-junction gene annotation table into SJVC
- POST /sjv/{sid}/junctions          load a junction count matrix into an SJV session
- POST /sjv/{sid}/sample-metadata    load a sample-metadata table into an SJV session
- POST /sjsurv/{sid}/sjdat/{kind}    load an sjdat matrix into an SJSurv session
- POST /sjsurv/{sid}/metadata        load the sample metadata into an SJSurv session
- POST /sjsurv/{sid}/junction-metadata load the per-junction gene annotation table into SJSurv
- POST /sjlookup/{sid}/junction-metadata load the per-junction metadata table into SJ Lookup
- POST /sjlookup/{sid}/sjdat/{kind}      load junction_counts / rrs_scores into an SJ Lookup
                                     session — the per-sample table for a single-junction lookup
"""
from __future__ import annotations

import shutil
from typing import List

from fastapi import APIRouter, HTTPException

from sjsurv.api.routes import (
    ingest_junction_metadata as sjsurv_ingest_junction_metadata,
    ingest_metadata as sjsurv_ingest_metadata,
    ingest_sjdat as sjsurv_ingest_sjdat,
)
from sjsurv.models import (
    JunctionMetadataLoaded as SjsurvJunctionMetadataLoaded,
    MetadataLoaded as SjsurvMetadataLoaded,
    SjdatLoaded as SjsurvSjdatLoaded,
)
from sjsurv.services.sessions import store as sjsurv_store
from sjlookup.api.routes import (
    ingest_junction_metadata as sjlookup_ingest_junction_metadata,
    ingest_sjdat as sjlookup_ingest_sjdat,
)
from sjlookup.models import (
    JunctionMetadataLoaded as SjlookupJunctionMetadataLoaded,
    SjdatLoaded as SjlookupSjdatLoaded,
)
from sjlookup.services.sessions import store as sjlookup_store
from sjv.api.routes import ingest_rds, ingest_sample_metadata
from sjv.models import RdsUploaded, SampleMetadataUploaded
from sjv.services import gencode as sjv_gencode
from sjv.services.sessions import store as sjv_store
from sjvc.api.routes import (
    ingest_clinical,
    ingest_junction_metadata as sjvc_ingest_junction_metadata,
    ingest_junctions,
    ingest_sjdat as sjvc_ingest_sjdat,
)
from sjvc.models import ClinicalUploaded, JunctionMetadataLoaded, MatrixUploaded, SjdatLoaded as SjvcSjdatLoaded
from sjvc.services import gencode as sjvc_gencode
from sjvc.services.sessions import store as sjvc_store

from .models import (
    GencodeApplied,
    GencodeApply,
    GtfApply,
    PathBody,
    PushGroupsBody,
    PushGroupsResult,
    ScanResponse,
)
from .paths import safe_path
from .push_groups import push_into_clinical, push_into_sample_metadata
from .scan import classify

router = APIRouter()


@router.post("/scan", response_model=ScanResponse)
def scan(body: PathBody) -> ScanResponse:
    d = safe_path(body.path)
    if not d.is_dir():
        raise HTTPException(400, f"not a directory: {d}")
    return ScanResponse(cohort=d.name, dir=str(d), files=classify(d))


# --------------------------------------------------------------------------- #
# GENCODE reference — chosen once here, applied to both app sessions
# --------------------------------------------------------------------------- #
def _apply_annotation(body, make_sjv, make_sjvc, make_sjsurv) -> GencodeApplied:
    """`make_sjv` / `make_sjvc` / `make_sjsurv` each build an Annotation for
    their own app (the GENCODE index cache is shared, so a release is
    downloaded once even though each app keeps its own SQLite schema — and
    SJSurv reuses sjvc's own gencode module directly, it has none of its
    own)."""
    applied: List[str] = []
    label = ""
    warnings: List[str] = []
    species = getattr(body, "species", "human")

    try:
        if body.sjv_sid and (s := sjv_store.get(body.sjv_sid)):
            ann = make_sjv()
            s.annotation = ann
            label = ann.label
            applied.append("sjv")
        if body.sjvc_sid and (s := sjvc_store.get(body.sjvc_sid)):
            ann = make_sjvc()
            s.annotation = ann
            s.species = species
            label = ann.label
            applied.append("sjvc")
        if body.sjsurv_sid and (s := sjsurv_store.get(body.sjsurv_sid)):
            ann = make_sjsurv()
            s.annotation = ann
            s.species = species
            label = ann.label
            applied.append("sjsurv")
    except (sjv_gencode.GencodeError, sjvc_gencode.GencodeError) as e:
        raise HTTPException(502, str(e))

    if not applied:
        raise HTTPException(404, "no live SJV, SJVC, or SJSurv session to apply the reference to")
    if species != "human":
        warnings.append("pathway libraries are human-only")
    return GencodeApplied(
        label=label,
        pathways_enabled=species == "human",
        applied_to=applied,
        warnings=warnings,
    )


@router.post("/gencode", response_model=GencodeApplied)
def apply_gencode(body: GencodeApply) -> GencodeApplied:
    return _apply_annotation(
        body,
        lambda: sjv_gencode.ensure_release(body.species, body.release),
        lambda: sjvc_gencode.ensure_release(body.species, body.release),
        lambda: sjvc_gencode.ensure_release(body.species, body.release),
    )


@router.post("/gtf", response_model=GencodeApplied)
def apply_gtf(body: GtfApply) -> GencodeApplied:
    src = _file_arg(body.path)
    return _apply_annotation(
        body,
        lambda: sjv_gencode.annotation_from_gtf(src, label=src.name),
        lambda: sjvc_gencode.annotation_from_gtf(src, label=src.name),
        lambda: sjvc_gencode.annotation_from_gtf(src, label=src.name),
    )


@router.get("/gencode/releases")
def gencode_releases() -> dict:
    return {"releases": sjv_gencode.list_releases()}


def _file_arg(path: str):
    src = safe_path(path)
    if not src.is_file():
        raise HTTPException(400, f"not a file: {src}")
    return src


def _session(store, sid: str):
    s = store.get(sid)
    if s is None:
        raise HTTPException(404, "unknown or expired session")
    return s


def _copy_in(session, name: str, src) -> "object":
    dest = session.tmp_dir / name
    shutil.copy(src, dest)
    return dest


# --------------------------------------------------------------------------- #
# push SJSurv's Group/SurviverGroup to the other two apps
# --------------------------------------------------------------------------- #
@router.post("/push-groups", response_model=PushGroupsResult)
def push_groups(body: PushGroupsBody) -> PushGroupsResult:
    """Push SJSurv's currently-computed Group/SurviverGroup into whichever of
    the Sashimi-plot / 2D View sessions already have a metadata / clinical table
    loaded — see push_groups.py for why this overwrites same-named columns."""
    s_surv = _session(sjsurv_store, body.sjsurv_sid)
    if s_surv.metadata is None:
        raise HTTPException(409, "stratify the sample metadata on the SJSurv tab first")
    group_map = {sid: (r["group"], r["surviver"]) for sid, r in s_surv.metadata.rows.items()}

    applied: List[str] = []
    matched: dict = {}
    warnings: List[str] = []

    if body.sjv_sid and (s := sjv_store.get(body.sjv_sid)):
        if s.metadata is None:
            warnings.append(
                "Sashimi plot: no sample-metadata table loaded there yet — load one on the "
                "Data tab first, then push again"
            )
        else:
            matched["sjv"] = push_into_sample_metadata(s.metadata, group_map)
            applied.append("sjv")

    if body.sjvc_sid and (s := sjvc_store.get(body.sjvc_sid)):
        if s.clinical is None:
            warnings.append(
                "2D View: no clinical table loaded there yet — load one on the Data tab first, "
                "then push again"
            )
        else:
            matched["sjvc"] = push_into_clinical(s.clinical, group_map)
            applied.append("sjvc")

    if not applied:
        raise HTTPException(
            404, "no live Sashimi-plot or 2D View session with a metadata/clinical table loaded"
        )

    return PushGroupsResult(
        applied_to=applied,
        n_sjsurv_samples=len(group_map),
        n_group=sum(1 for g, _ in group_map.values() if g),
        n_labelled=sum(1 for _, sv in group_map.values() if sv in ("Good", "Poor")),
        matched=matched,
        warnings=warnings,
    )


# --- SJVC (cohort explorer) ------------------------------------------------- #
@router.post("/sjvc/{sid}/gene-matrix", response_model=MatrixUploaded)
def load_gene_matrix(sid: str, body: PathBody) -> MatrixUploaded:
    s = _session(sjvc_store, sid)
    return ingest_junctions(s, _copy_in(s, "junctions.rds", _file_arg(body.path)))


@router.post("/sjvc/{sid}/clinical", response_model=ClinicalUploaded)
def load_clinical(sid: str, body: PathBody) -> ClinicalUploaded:
    s = _session(sjvc_store, sid)
    src = _file_arg(body.path)
    return ingest_clinical(s, _copy_in(s, f"clinical{''.join(src.suffixes) or '.rds'}", src))


@router.post("/sjvc/{sid}/sjdat/{kind}", response_model=SjvcSjdatLoaded)
def load_sjvc_sjdat(sid: str, kind: str, body: PathBody) -> SjvcSjdatLoaded:
    s = _session(sjvc_store, sid)
    return sjvc_ingest_sjdat(s, kind, _copy_in(s, f"sjdat_{kind}.rds", _file_arg(body.path)))


@router.post("/sjvc/{sid}/junction-metadata", response_model=JunctionMetadataLoaded)
def load_sjvc_junction_metadata(sid: str, body: PathBody) -> JunctionMetadataLoaded:
    s = _session(sjvc_store, sid)
    src = _file_arg(body.path)
    return sjvc_ingest_junction_metadata(
        s, _copy_in(s, f"junction_metadata{''.join(src.suffixes) or '.rds'}", src)
    )


# --- SJV (sashimi plot) ---------------------------------------------------- #
@router.post("/sjv/{sid}/junctions", response_model=RdsUploaded)
def load_junctions(sid: str, body: PathBody) -> RdsUploaded:
    s = _session(sjv_store, sid)
    return ingest_rds(s, _copy_in(s, "upload.rds", _file_arg(body.path)))


@router.post("/sjv/{sid}/sample-metadata", response_model=SampleMetadataUploaded)
def load_sample_metadata(sid: str, body: PathBody) -> SampleMetadataUploaded:
    s = _session(sjv_store, sid)
    src = _file_arg(body.path)
    return ingest_sample_metadata(s, _copy_in(s, f"metadata{''.join(src.suffixes) or '.rds'}", src))


# --- SJSurv (survivor-group classification) ------------------------------- #
@router.post("/sjsurv/{sid}/sjdat/{kind}", response_model=SjsurvSjdatLoaded)
def load_sjsurv_sjdat(sid: str, kind: str, body: PathBody) -> SjsurvSjdatLoaded:
    s = _session(sjsurv_store, sid)
    return sjsurv_ingest_sjdat(s, kind, _copy_in(s, f"sjdat_{kind}.rds", _file_arg(body.path)))


@router.post("/sjsurv/{sid}/metadata", response_model=SjsurvMetadataLoaded)
def load_sjsurv_metadata(sid: str, body: PathBody) -> SjsurvMetadataLoaded:
    s = _session(sjsurv_store, sid)
    src = _file_arg(body.path)
    return sjsurv_ingest_metadata(s, _copy_in(s, f"metadata{''.join(src.suffixes) or '.rds'}", src))


@router.post("/sjsurv/{sid}/junction-metadata", response_model=SjsurvJunctionMetadataLoaded)
def load_sjsurv_junction_metadata(sid: str, body: PathBody) -> SjsurvJunctionMetadataLoaded:
    s = _session(sjsurv_store, sid)
    src = _file_arg(body.path)
    return sjsurv_ingest_junction_metadata(
        s, _copy_in(s, f"junction_metadata{''.join(src.suffixes) or '.rds'}", src)
    )


# --- SJ Lookup (per-junction lookup) --------------------------------------- #
@router.post("/sjlookup/{sid}/junction-metadata", response_model=SjlookupJunctionMetadataLoaded)
def load_sjlookup_junction_metadata(sid: str, body: PathBody) -> SjlookupJunctionMetadataLoaded:
    s = _session(sjlookup_store, sid)
    src = _file_arg(body.path)
    return sjlookup_ingest_junction_metadata(
        s, _copy_in(s, f"junction_metadata{''.join(src.suffixes) or '.rds'}", src)
    )


@router.post("/sjlookup/{sid}/sjdat/{kind}", response_model=SjlookupSjdatLoaded)
def load_sjlookup_sjdat(sid: str, kind: str, body: PathBody) -> SjlookupSjdatLoaded:
    s = _session(sjlookup_store, sid)
    return sjlookup_ingest_sjdat(s, kind, _copy_in(s, f"sjdat_{kind}.rds", _file_arg(body.path)))
