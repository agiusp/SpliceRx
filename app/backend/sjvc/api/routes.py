"""HTTP API for SJVC."""
from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Dict, List

import numpy as np
from fastapi import APIRouter, File, HTTPException, UploadFile
from fastapi.responses import Response

from .. import palette
from ..models import (
    ActivateSjdat,
    ClinicalColumn,
    ClinicalUploaded,
    FeaturesRequest,
    FeaturesResponse,
    GencodeReady,
    GencodeSelect,
    GeneSetRequest,
    GeneSetResponse,
    HeatmapRequest,
    JunctionMetadataLoaded,
    MadFeaturesRequest,
    MatrixUploaded,
    ProjectionRequest,
    SessionCreated,
    SessionState,
    SjdatLoaded,
    SjdatOption,
)
from ..services import encoding, features as feat_mod, gencode, geneset as gs_mod, heatmap as hm_mod, mad as mad_mod, pathways, projection as proj_mod
from ..services.clinical import parse_clinical
from ..services.junction_metadata import JunctionMetadataError, load_junction_gene_index
from ..services.junction_type import load_junction_type_index
from ..services.junctions import gene_name_of_label, looks_gene_level, map_junctions_to_genes
from ..services.rds import MatrixError
from ..services.sjdat import SJDAT_KINDS, SJDAT_META, SjdatError, load_sjdat
from ..services.sessions import store

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


@router.get("/health")
def health() -> dict:
    return {"ok": True}


@router.post("/session", response_model=SessionCreated)
def create_session() -> SessionCreated:
    return SessionCreated(session_id=store.create().id)


def _sync_sjdat_intersection(s) -> None:
    """Narrow every loaded sjdat matrix to the samples that also have a
    clinical row, recomputed from ``s.sjdat_raw`` (the pristine, as-loaded
    matrices) so re-uploading the clinical table never narrows the working
    sample set beyond what the *current* files actually share. Mirrors
    ``sjsurv.api.routes._sync_sjdat_intersection``."""
    known = set(s.clinical.rows) if s.clinical is not None else None
    for kind, raw in s.sjdat_raw.items():
        s.sjdat[kind] = raw.restrict_samples(known) if known is not None else raw
    _sync_active_junctions(s)


def _sync_active_junctions(s) -> None:
    """Refresh the `junctions`/`junctions_raw` "current view" fields from
    `sjdat`/`sjdat_raw[active_sjdat]` — every function written against a
    single matrix (gene-set matching, MAD, feature building, PCA/UMAP,
    heatmap) reads these two and needs no further change."""
    s.junctions = s.sjdat.get(s.active_sjdat)
    s.junctions_raw = s.sjdat_raw.get(s.active_sjdat)


def _sjdat_loaded(s, kind: str, n_excluded_label: str) -> SjdatLoaded:
    d = s.sjdat[kind]
    warnings: List[str] = []
    if looks_gene_level(d.features):
        warnings.append(
            "row names aren't chr:start-end:strand — treating this as an already-condensed, "
            "gene-level matrix. Use 'Top by MAD' to pick features from it."
        )
    n_excluded = s.sjdat_raw[kind].n_samples - d.n_samples
    if n_excluded:
        warnings.append(f"{n_excluded} matrix sample(s) excluded — {n_excluded_label}")
    return SjdatLoaded(kind=kind, n_features=d.n_features, n_samples=d.n_samples, sparse=d.sparse, warnings=warnings)


def ingest_sjdat(s, kind: str, src: Path) -> SjdatLoaded:
    """Load one of the 3 sjdat matrices (junction counts / RRS scores / gene
    matrix) at ``src`` into the session. ``src`` is consumed. Shared by the
    browser-upload endpoint and the dataload package.

    If a clinical table is already loaded, every loaded matrix (this one
    included) is immediately narrowed to the samples they all have in common
    — the working sample set for every downstream step (feature building,
    MAD, PCA/UMAP, heatmap) is always the intersection, never a raw matrix
    alone."""
    if kind not in SJDAT_KINDS:
        src.unlink(missing_ok=True)
        raise HTTPException(422, f"unknown sjdat kind {kind!r}")
    try:
        d = load_sjdat(kind, src, s.tmp_dir)
    except SjdatError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)

    s.sjdat_raw[kind] = d
    first_load = s.active_sjdat is None
    if first_load:
        s.active_sjdat = kind
    _sync_sjdat_intersection(s)
    if first_load or kind == s.active_sjdat:
        # (re)loading the *active* matrix invalidates any gene set / feature
        # matrix built against its old contents; a different, inactive kind
        # loading alongside it doesn't touch what's currently being analysed
        s.geneset = None
        s.features = None
    return _sjdat_loaded(s, kind, "not in the loaded clinical table")


# Old single-gene-matrix endpoint — kept working for any existing caller by
# always loading as "gene_matrix" (its one real-world use). New code should
# use `ingest_sjdat` / `POST /session/{sid}/sjdat/{kind}` directly.
def ingest_junctions(s, src: Path) -> MatrixUploaded:
    loaded = ingest_sjdat(s, "gene_matrix", src)
    return MatrixUploaded(
        samples=s.junctions.samples, n_junctions=loaded.n_features,
        feature_kind="gene" if looks_gene_level(s.junctions.features) else "junction",
        warnings=loaded.warnings,
    )


def ingest_clinical(s, src: Path) -> ClinicalUploaded:
    """Load the clinical table at ``src`` (a path inside ``s.tmp_dir``, extension
    preserved so the reader can dispatch on it) into the session. ``src`` is
    consumed. Shared by the browser-upload endpoint and the dataload package.

    Every loaded sjdat matrix is then narrowed to the samples common to it
    and this table — see ``_sync_sjdat_intersection``."""
    if not s.sjdat_raw:
        raise HTTPException(409, "load a matrix first")
    known_samples = sorted({sid for d in s.sjdat_raw.values() for sid in d.samples})
    try:
        clin = parse_clinical(src, s.tmp_dir, known_samples)
    except MatrixError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)
    s.clinical = clin
    s.features = None
    _sync_sjdat_intersection(s)
    warnings: List[str] = []
    if clin.unmatched_samples:
        warnings.append(
            f"{len(clin.unmatched_samples)} sample_id value(s) not in any loaded matrix "
            f"were ignored (e.g. {', '.join(clin.unmatched_samples[:3])})"
        )
    for kind, raw in s.sjdat_raw.items():
        n_excluded = raw.n_samples - s.sjdat[kind].n_samples
        if n_excluded:
            warnings.append(
                f"{n_excluded} {SJDAT_META[kind][0].lower()} sample(s) excluded — "
                f"no matching clinical row"
            )
    return ClinicalUploaded(
        columns=[ClinicalColumn(**c.__dict__) for c in clin.columns],
        n_matched=clin.n_rows,
        warnings=warnings,
    )


def ingest_junction_metadata(s, src: Path) -> JunctionMetadataLoaded:
    """Load the per-junction annotation table (``TCGA_<cohort>_junction_
    metadata.rds``) — enables fast typed-gene / pathway lookups against a
    junction-level matrix (junction counts, RRS scores) without needing a
    GENCODE release at all. See ``services/junction_metadata.py``.

    Also builds the gene-name/splice-type index the junction-level heatmap's
    alternate row-label view uses (``services/junction_type.py``) from the
    same file, before it's removed below — this is best-effort: a table this
    app can otherwise use just fine but that lacks a gene-name column simply
    doesn't get that view, rather than failing the whole load."""
    try:
        idx = load_junction_gene_index(src, s.tmp_dir)
    except JunctionMetadataError as e:
        src.unlink(missing_ok=True)
        raise HTTPException(422, str(e))
    try:
        type_idx = load_junction_type_index(src, s.tmp_dir)
    except JunctionMetadataError:
        type_idx = None
    finally:
        src.unlink(missing_ok=True)
    s.junction_gene_index = idx
    s.junction_type_index = type_idx
    s.geneset = None
    s.features = None
    return JunctionMetadataLoaded(n_rows=idx.n_rows, n_genes=len(idx.gene_id_to_name))


@router.post("/session/{sid}/sjdat/{kind}", response_model=SjdatLoaded)
async def upload_sjdat(sid: str, kind: str, file: UploadFile = File(...)) -> SjdatLoaded:
    s = _session(sid)
    dest = s.tmp_dir / f"sjdat_{kind}.rds"
    await _save_upload(file, dest)
    return ingest_sjdat(s, kind, dest)


@router.post("/session/{sid}/sjdat", response_model=SessionState)
def activate_sjdat(sid: str, body: ActivateSjdat) -> SessionState:
    s = _session(sid)
    if body.kind not in s.sjdat:
        raise HTTPException(409, f"the {body.kind!r} matrix has not been loaded")
    if s.active_sjdat != body.kind:
        s.active_sjdat = body.kind
        _sync_active_junctions(s)
        s.geneset = None
        s.features = None
    return session_state(sid)


@router.post("/session/{sid}/junction-metadata", response_model=JunctionMetadataLoaded)
async def upload_junction_metadata(sid: str, file: UploadFile = File(...)) -> JunctionMetadataLoaded:
    s = _session(sid)
    suffix = "".join(Path(file.filename or "junction_metadata.rds").suffixes) or ".rds"
    dest = s.tmp_dir / f"junction_metadata{suffix}"
    await _save_upload(file, dest)
    return ingest_junction_metadata(s, dest)


@router.post("/session/{sid}/junctions", response_model=MatrixUploaded)
async def upload_junctions(sid: str, file: UploadFile = File(...)) -> MatrixUploaded:
    s = _session(sid)
    dest = s.tmp_dir / "junctions.rds"
    await _save_upload(file, dest)
    return ingest_junctions(s, dest)


@router.post("/session/{sid}/clinical", response_model=ClinicalUploaded)
async def upload_clinical(sid: str, file: UploadFile = File(...)) -> ClinicalUploaded:
    s = _session(sid)
    if not s.sjdat:
        raise HTTPException(409, "upload a matrix first")
    suffix = "".join(Path(file.filename or "clinical.csv").suffixes) or ".csv"
    dest = s.tmp_dir / f"clinical{suffix}"
    await _save_upload(file, dest)
    return ingest_clinical(s, dest)


def _sjdat_options(s) -> List[SjdatOption]:
    opts = []
    for kind in SJDAT_KINDS:
        label, desc = SJDAT_META[kind]
        d = s.sjdat.get(kind)
        opts.append(SjdatOption(
            kind=kind, label=label, description=desc,
            loaded=d is not None,
            n_features=d.n_features if d else 0,
            n_samples=d.n_samples if d else 0,
            sparse=bool(d and d.sparse),
            # cheap only for the small dense pathway matrix — the frontend
            # defaults "top n by MAD" to this when that's the active sjdat
            n_nonzero_rows=d.count_nonzero_rows() if (d is not None and kind == "pathway_matrix") else None,
        ))
    return opts


@router.get("/session/{sid}/state", response_model=SessionState)
def session_state(sid: str) -> SessionState:
    """Everything the frontend needs to rebuild step 1 of the UI from a session
    it did not populate itself (the Data tab hand-off) or after a page refresh."""
    s = _session(sid)
    j = s.junctions
    clinical = None
    if s.clinical is not None:
        warnings: List[str] = []
        if s.junctions_raw is not None and j is not None:
            n_excluded = s.junctions_raw.n_samples - j.n_samples
            if n_excluded:
                warnings.append(
                    f"{n_excluded} matrix sample(s) excluded — no matching clinical row"
                )
        clinical = ClinicalUploaded(
            columns=[ClinicalColumn(**c.__dict__) for c in s.clinical.columns],
            n_matched=s.clinical.n_rows,
            warnings=warnings,
        )
    return SessionState(
        has_junctions=j is not None,
        feature_kind=(
            ("gene" if looks_gene_level(j.features) else "junction") if j is not None else None
        ),
        samples=j.samples if j is not None else [],
        n_junctions=j.n_features if j is not None else 0,
        clinical=clinical,
        gencode_label=s.annotation.label if s.annotation is not None else None,
        pathways_enabled=s.species == "human" and s.annotation is not None,
        sjdat_options=_sjdat_options(s),
        active_sjdat=s.active_sjdat,
        has_junction_metadata=s.junction_gene_index is not None,
    )


@router.get("/gencode/releases")
def releases() -> dict:
    return {"releases": gencode.list_releases()}


@router.post("/session/{sid}/gencode", response_model=GencodeReady)
def select_gencode(sid: str, body: GencodeSelect) -> GencodeReady:
    s = _session(sid)
    try:
        ann = gencode.ensure_release(body.species, body.release)
    except gencode.GencodeError as e:
        raise HTTPException(502, str(e))
    s.annotation = ann
    s.species = body.species
    return GencodeReady(
        label=ann.label,
        pathways_enabled=body.species == "human",
        warnings=[] if body.species == "human" else ["pathway libraries are human-only"],
    )


@router.post("/session/{sid}/gtf", response_model=GencodeReady)
async def upload_gtf(sid: str, file: UploadFile = File(...)) -> GencodeReady:
    s = _session(sid)
    suffix = ".gtf.gz" if (file.filename or "").endswith(".gz") else ".gtf"
    dest = s.tmp_dir / f"annotation{suffix}"
    await _save_upload(file, dest)
    try:
        ann = gencode.annotation_from_gtf(dest, label=file.filename or "custom GTF")
    except gencode.GencodeError as e:
        raise HTTPException(422, str(e))
    s.annotation = ann
    s.species = "human"
    return GencodeReady(label=ann.label, pathways_enabled=True)


@router.get("/session/{sid}/genes/suggest")
def suggest_genes(sid: str, q: str, limit: int = 20) -> dict:
    s = _session(sid)
    if s.annotation is None:
        return {"names": []}
    return {"names": s.annotation.suggest(q, min(limit, 50))}


@router.get("/pathways/libraries")
def pathway_libraries() -> dict:
    return {"libraries": pathways.list_libraries()}


@router.get("/pathways/search")
def pathway_search(library: str, q: str = "", limit: int = 25) -> dict:
    try:
        return {"library": library, "terms": pathways.search(library, q, min(limit, 50))}
    except pathways.PathwayError as e:
        raise HTTPException(502, str(e))


@router.post("/session/{sid}/geneset", response_model=GeneSetResponse)
def set_geneset(sid: str, body: GeneSetRequest) -> GeneSetResponse:
    s = _session(sid)
    if s.junctions is None:
        raise HTTPException(409, "load a matrix first")
    gene_level = looks_gene_level(s.junctions.features)
    fast_lookup = not gene_level and s.junction_gene_index is not None
    if not gene_level and not fast_lookup and s.annotation is None:
        raise HTTPException(
            409, "select a GENCODE release or upload a GTF first (or load the junction "
                 "metadata file on the Data tab for a faster, no-GENCODE-needed lookup)"
        )

    # resolve the mode to a plain list of symbols
    try:
        if body.mode == "pathway":
            if not (body.library and body.term):
                raise HTTPException(422, "pathway mode needs `library` and `term`")
            symbols = pathways.genes(body.library, body.term)
            source = f"{body.library}: {body.term}"
            prefix = False
        elif body.mode in ("typed", "list"):
            symbols = gs_mod.parse_symbols(body.text or "")
            source = body.mode
            prefix = body.prefix
        else:
            raise HTTPException(422, f"unknown mode {body.mode!r}")
    except pathways.PathwayError as e:
        raise HTTPException(502, str(e))

    if gene_level:
        # match symbols straight against the matrix's gene row labels
        matched_names, unmatched, warns = gs_mod.match_matrix_labels(
            symbols, s.junctions.features, prefix=prefix
        )
        if not matched_names:
            raise HTTPException(422, f"none of the {len(symbols)} gene(s) are rows in this matrix")
        gs = gs_mod.GeneSet(
            matched_names=matched_names, unmatched=unmatched, source=source, warnings=warns
        )
        resp_matched = matched_names
    elif fast_lookup:
        # junction-level matrix, junction metadata loaded: a plain lookup,
        # no GENCODE release needed at all (see services/junction_metadata.py)
        jgmap, labels, unmatched = s.junction_gene_index.resolve(symbols, prefix=prefix)
        if not labels:
            raise HTTPException(
                422, f"none of the {len(symbols)} gene(s) have a junction in the loaded "
                     f"junction metadata"
            )
        warns: List[str] = []
        n_dropped = sum(1 for rn in jgmap.rownames if rn not in s.junctions._row)
        if n_dropped:
            warns.append(
                f"{n_dropped} junction(s) from the metadata aren't columns of the active "
                f"matrix and were skipped"
            )
        gs = gs_mod.GeneSet(
            matched_junction_map=jgmap, matched_gene_labels=labels,
            unmatched=unmatched, source=source, warnings=warns,
        )
        resp_matched = list(labels.values())
    else:
        gs = gs_mod.resolve(symbols, s.annotation, source=source, prefix=prefix)
        if not gs.matched:
            raise HTTPException(422, f"none of the {len(gs.unmatched)} gene(s) matched the reference")
        resp_matched = [g.name for g in gs.matched]

    s.geneset = gs
    s.features = None
    return GeneSetResponse(
        matched=resp_matched,
        unmatched=gs.unmatched,
        n_genes=gs.n_genes,
        source=gs.source,
        warnings=gs.warnings,
    )


@router.post("/session/{sid}/features", response_model=FeaturesResponse)
def build_features(sid: str, body: FeaturesRequest) -> FeaturesResponse:
    s = _session(sid)
    if s.junctions is None or s.geneset is None:
        raise HTTPException(409, "upload a matrix and choose a gene set first")

    # gene-level matrix: the gene set is a plain list of row labels to keep
    if getattr(s.geneset, "matched_names", None):
        try:
            fm = feat_mod.select_by_labels(s.junctions, s.geneset.matched_names)
        except ValueError as e:
            raise HTTPException(422, str(e))
        mad_note = None
        if body.mad_top_n and fm.n_features > body.mad_top_n:
            mad_note = f"ranked {fm.n_features} resolved feature(s) by MAD, kept the top {body.mad_top_n}"
            fm = mad_mod.refine_by_mad(fm, body.mad_top_n)
        s.features = fm
        return FeaturesResponse(
            feature_kind="gene",
            n_features=fm.n_features,
            n_samples=fm.n_samples,
            n_junctions=0,
            dropped_zero_variance=0,
            suggest_condense=False,
            feature_preview=[fm.label(i) for i in fm.feature_ids[:20]],
            warnings=[mad_note] if mad_note else [],
        )

    # junction-level matrix, resolved via the junction-metadata fast lookup:
    # already have exactly which junctions overlap which (selected) genes,
    # no GENCODE overlap computation needed
    if getattr(s.geneset, "matched_junction_map", None) is not None:
        jgmap = s.geneset.matched_junction_map
        try:
            fm = feat_mod.select_and_maybe_condense(
                s.junctions, jgmap, condense=body.condense, gene_labels=s.geneset.matched_gene_labels,
            )
        except ValueError as e:
            raise HTTPException(422, str(e))
        mad_note = None
        if body.mad_top_n and fm.n_features > body.mad_top_n:
            mad_note = f"ranked {fm.n_features} resolved feature(s) by MAD, kept the top {body.mad_top_n}"
            fm = mad_mod.refine_by_mad(fm, body.mad_top_n)
        s.features = fm
        return FeaturesResponse(
            feature_kind=fm.kind,
            n_features=fm.n_features,
            n_samples=fm.n_samples,
            n_junctions=fm.n_junctions,
            dropped_zero_variance=0,
            suggest_condense=jgmap.n_junctions() > feat_mod.CONDENSE_SUGGEST_ABOVE,
            feature_preview=[fm.label(i) for i in fm.feature_ids[:20]],
            warnings=[mad_note] if mad_note else [],
        )

    if s.annotation is None:
        raise HTTPException(
            409, "select a GENCODE release or upload a GTF first (or load the junction "
                 "metadata file on the Data tab for a faster, no-GENCODE-needed lookup)"
        )

    jgmap = map_junctions_to_genes(s.junctions.features, s.geneset.matched)
    if jgmap.n_junctions() == 0:
        raise HTTPException(
            422, "no junction in the matrix overlaps the selected genes (check genome build / naming)"
        )

    labels = {g.gene_id: g.name for g in s.geneset.matched}
    try:
        fm = feat_mod.select_and_maybe_condense(
            s.junctions, jgmap, condense=body.condense, gene_labels=labels
        )
    except ValueError as e:
        raise HTTPException(422, str(e))

    warnings: List[str] = []
    if jgmap.n_bad_rownames:
        warnings.append(f"{jgmap.n_bad_rownames} matrix row name(s) were not chr:start-end:strand and skipped")
    if body.mad_top_n and fm.n_features > body.mad_top_n:
        warnings.append(f"ranked {fm.n_features} resolved feature(s) by MAD, kept the top {body.mad_top_n}")
        fm = mad_mod.refine_by_mad(fm, body.mad_top_n)
    s.features = fm
    return FeaturesResponse(
        feature_kind=fm.kind,
        n_features=fm.n_features,
        n_samples=fm.n_samples,
        n_junctions=fm.n_junctions,
        dropped_zero_variance=0,
        suggest_condense=jgmap.n_junctions() > feat_mod.CONDENSE_SUGGEST_ABOVE,
        feature_preview=[fm.label(i) for i in fm.feature_ids[:20]],
        warnings=warnings,
    )


@router.post("/session/{sid}/features/mad", response_model=FeaturesResponse)
def build_features_mad(sid: str, body: MadFeaturesRequest) -> FeaturesResponse:
    """Alternative to picking a gene set: rank the matrix's rows (junctions,
    or genes if it is already condensed) by descending MAD and keep the top
    N. Bypasses gene-set resolution entirely — the ranking itself is the
    selection."""
    s = _session(sid)
    if s.junctions is None:
        raise HTTPException(409, "upload the junction matrix first")
    if body.top_n < 1:
        raise HTTPException(422, "top_n must be >= 1")

    restrict_to = None
    if body.protein_coding_only:
        if s.annotation is None:
            raise HTTPException(
                409, "select a GENCODE release or upload a GTF to filter to protein-coding genes"
            )
        if not s.annotation.has_gene_types():
            raise HTTPException(
                422, "the selected reference has no gene_type / biotype attributes to filter on"
            )
        restrict_to = s.annotation.gene_names_of_type("protein_coding")

    try:
        if body.condense:
            if s.annotation is None:
                raise HTTPException(409, "select a GENCODE release or upload a GTF first")
            fm = mad_mod.top_genes_by_mad(s.junctions, s.annotation, body.top_n)
        else:
            fm = mad_mod.top_features_by_mad(
                s.junctions, body.top_n, restrict_to=restrict_to,
                n_min=body.n_min, x_min=body.x_min, annotation=s.annotation,
            )
    except ValueError as e:
        raise HTTPException(422, str(e))

    s.features = fm
    s.geneset = None  # this selection mode doesn't go through a resolved gene set

    warnings: List[str] = []
    if restrict_to is not None:
        if fm.kind == "gene":
            n_coding = sum(
                1 for f in s.junctions.features if gene_name_of_label(f).strip().lower() in restrict_to
            )
            warnings.append(
                f"{n_coding} of {len(s.junctions.features)} matrix genes are protein-coding "
                f"in the reference — ranked those"
            )
        else:
            warnings.append("restricted ranking to junctions overlapping a protein-coding gene in the reference")
    if fm.n_features < body.top_n:
        warnings.append(
            f"only {fm.n_features} feature(s) have any variability (nonzero MAD) — asked for "
            f"top {body.top_n}"
        )
    return FeaturesResponse(
        feature_kind=fm.kind,
        n_features=fm.n_features,
        n_samples=fm.n_samples,
        n_junctions=fm.n_junctions,
        dropped_zero_variance=0,
        suggest_condense=False,
        feature_preview=[fm.label(i) for i in fm.feature_ids[:20]],
        warnings=warnings,
    )


@router.get("/session/{sid}/features.csv")
def features_csv(sid: str) -> Response:
    """The current feature matrix as CSV — features (rows) x samples
    (columns), raw values before preprocessing. This is the data the
    heatmap / projection are built from."""
    s = _session(sid)
    if s.features is None:
        raise HTTPException(409, "build the feature matrix first")
    fm = s.features

    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["feature", *fm.samples])
    for i, fid in enumerate(fm.feature_ids):
        w.writerow(
            [fm.label(fid), *("" if not np.isfinite(v) else f"{float(v):g}" for v in fm.values[i])]
        )

    kind = "genes" if fm.kind == "gene" else "junctions"
    return Response(
        content=buf.getvalue(),
        media_type="text/csv",
        headers={"content-disposition": f'attachment; filename="heatmap_data_{kind}.csv"'},
    )


def _apply_overrides(s, overrides: Dict[str, str]) -> None:
    if overrides and s.clinical is not None:
        try:
            s.clinical.apply_overrides(overrides)
        except MatrixError as e:
            raise HTTPException(422, str(e))


@router.post("/session/{sid}/projection")
def projection(sid: str, body: ProjectionRequest) -> dict:
    s = _session(sid)
    if s.features is None:
        raise HTTPException(409, "build the feature matrix first")
    _apply_overrides(s, body.overrides)

    prep = feat_mod.preprocess(s.features, standardize=True)
    try:
        if body.method == "umap":
            emb = proj_mod.umap(
                prep.X, prep.samples,
                n_neighbors=body.n_neighbors, min_dist=body.min_dist,
            )
            xi, yi = 0, 1
        else:
            emb = proj_mod.pca(prep.X, prep.samples)
            xi = max(0, min(body.pc_x - 1, emb.coords.shape[1] - 1))
            yi = max(0, min(body.pc_y - 1, emb.coords.shape[1] - 1))
    except proj_mod.ProjectionError as e:
        raise HTTPException(422, str(e))

    clinical_feats = [f for f in body.clinical[:2]]
    enc = encoding.build(s.clinical, clinical_feats, prep.samples)

    clin_vals: Dict[str, list] = {}
    if s.clinical is not None:
        for f in clinical_feats:
            vals, kind = s.clinical.values_for(f, prep.samples)
            clin_vals[f] = [None if (kind == "numeric" and not np.isfinite(v)) else v for v in (vals.tolist() if kind == "numeric" else vals)]

    points = []
    for i, sample in enumerate(prep.samples):
        points.append({
            "sample": sample,
            "x": float(emb.coords[i, xi]),
            "y": float(emb.coords[i, yi]),
            "color": enc["colors"][i],
            "shape": enc["shapes"][i],
            # True when this point is drawn grey — missing the selected
            # clinical feature(s) — so the frontend can offer "hide these"
            # without another request (the projection itself is unaffected).
            "missing": bool(enc["missing"][i]),
            "clinical": {f: clin_vals.get(f, [None] * len(prep.samples))[i] for f in clinical_feats},
        })

    return {
        "method": emb.method,
        "axis_labels": [emb.axis_labels[xi], emb.axis_labels[yi]] if emb.method == "pca" else emb.axis_labels,
        "explained_variance": emb.explained_variance,
        "n_components": int(emb.coords.shape[1]),
        "pc_x": xi + 1, "pc_y": yi + 1,
        "points": points,
        "encoding": enc["encoding"],
        "legend": enc["legend"],
        "warnings": ([f"{prep.dropped} constant feature(s) dropped"] if prep.dropped else []),
    }


@router.post("/session/{sid}/heatmap")
def heatmap(sid: str, body: HeatmapRequest) -> dict:
    s = _session(sid)
    if s.features is None:
        raise HTTPException(409, "build the feature matrix first")
    _apply_overrides(s, body.overrides)

    # Drop samples missing any selected clinical annotation before clustering
    # — unlike the projection's equivalent toggle (a post-hoc point filter,
    # no recompute needed), a heatmap column is baked into the clustering
    # itself, so this has to rebuild it from a narrower feature matrix.
    fm = s.features
    extra_warnings: List[str] = []
    if body.drop_missing_clinical and body.clinical:
        if s.clinical is None:
            raise HTTPException(422, "no clinical table loaded to check for missing values against")
        keep = np.ones(fm.n_samples, dtype=bool)
        for f in body.clinical:
            vals, kind = s.clinical.values_for(f, fm.samples)
            present = np.isfinite(vals) if kind == "numeric" else np.array([v is not None for v in vals])
            keep &= present
        n_dropped = int((~keep).sum())
        if n_dropped:
            if not keep.any():
                raise HTTPException(
                    422,
                    "every sample is missing at least one of the selected clinical feature(s) "
                    "— nothing left to show",
                )
            fm = feat_mod.subset_samples(fm, keep)
            extra_warnings.append(
                f"dropped {n_dropped} of {len(keep)} sample(s) missing "
                + ", ".join(body.clinical)
            )

    group_values = None
    if body.order == "group":
        if not body.group_by or s.clinical is None:
            raise HTTPException(422, "group order needs `group_by` and a clinical table")
        gv, _ = s.clinical.values_for(body.group_by, fm.samples)
        group_values = [None if v is None else str(v) for v in gv]

    try:
        res = hm_mod.build(
            fm, row_zscore=body.row_zscore, order=body.order, group_values=group_values,
        )
    except ValueError as e:
        raise HTTPException(422, str(e))
    res.warnings = extra_warnings + res.warnings

    # clinical annotation tracks, aligned to the displayed column order.
    # Successive categorical tracks take non-overlapping hues, and successive
    # numeric tracks take different single-hue ramps, so the stacked bars and
    # their legends stay visually distinct.
    tracks = []
    hue_offset = 0
    ramp_idx = 0
    if s.clinical is not None:
        for f in body.clinical:
            vals, kind = s.clinical.values_for(f, res.samples)
            if kind == "numeric":
                finite = vals[np.isfinite(vals)]
                lo = float(finite.min()) if finite.size else 0.0
                hi = float(finite.max()) if finite.size else 1.0
                span = hi - lo or 1.0
                ramp = ramp_idx
                ramp_idx += 1
                colors = [
                    "#cfcfcf" if not np.isfinite(v) else palette.sequential_color((v - lo) / span, ramp)
                    for v in vals
                ]
                tracks.append({
                    "feature": f, "type": "numeric",
                    "values": [None if not np.isfinite(v) else float(v) for v in vals],
                    "colors": colors, "min": lo, "max": hi,
                    "stops": palette.SEQUENTIAL_RAMPS[ramp % len(palette.SEQUENTIAL_RAMPS)],
                })
            else:
                distinct = sorted({v for v in vals if v is not None}, key=palette.natural_key)
                cmap = palette.categorical_colors(distinct, offset=hue_offset)
                hue_offset += len(distinct)
                tracks.append({
                    "feature": f, "type": "categorical",
                    "values": vals,
                    "colors": ["#cfcfcf" if v is None else cmap.get(v, palette.OTHER[0]) for v in vals],
                    "legend": [{"value": v, "color": cmap[v]} for v in distinct],
                })

    # Alternate row labels for a junction-level heatmap: "<gene name>:<novel-
    # splicing-event type>" instead of the raw "chr:start-end:strand" id,
    # from the cohort's own junction metadata table (see
    # services/junction_type.py) — only offered when that table was loaded
    # and this heatmap is junction-level (a gene/pathway matrix is already
    # one row per gene, nothing to map). The frontend toggles between this
    # and `row_labels` without another request, since only the label text
    # differs — row order and values are identical either way.
    row_labels_gene_type = None
    if fm.kind == "junction" and s.junction_type_index is not None:
        row_labels_gene_type = s.junction_type_index.labels_for(res.feature_ids)
        n_unlabelled = sum(1 for label in row_labels_gene_type if label is None)
        if n_unlabelled:
            res.warnings.append(
                f"{n_unlabelled} of {len(row_labels_gene_type)} displayed junction(s) have no "
                "gene annotation in the junction metadata table — shown by coordinate in the "
                "gene/type row-label view"
            )

    return {
        "feature_kind": fm.kind,
        # trim a "<gene_name>:<gene_id>" row label (count_novel_sjs() output,
        # e.g. TP53:ENSG00000141510.19) down to just the gene name for display;
        # junction (chr:start-end:strand) and pathway row labels are left
        # alone — gene_name_of_label only strips a trailing real Ensembl id
        "row_labels": [gene_name_of_label(fm.label(i)) for i in res.feature_ids],
        "row_labels_gene_type": row_labels_gene_type,
        "row_ids": res.feature_ids,
        "samples": res.samples,
        "values": [[float(x) for x in row] for row in res.values],
        "row_zscore": body.row_zscore,
        "row_dendro": res.row_dendro.segments if res.row_dendro else None,
        "col_dendro": res.col_dendro.segments if res.col_dendro else None,
        "annotations": tracks,
        "warnings": res.warnings,
    }
