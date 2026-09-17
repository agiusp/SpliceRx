"""HTTP API for SJSurv — mounted at /api/sjsurv."""
from __future__ import annotations

from pathlib import Path
from typing import List

from fastapi import APIRouter, File, HTTPException, UploadFile

from ..models import (
    ActivateSjdat,
    AgeBandsOut,
    CovariateColumnOut,
    CovariatesRequest,
    CVRequest,
    CVResponse,
    FeaturesRequest,
    FeaturesResponse,
    FeatureWeightOut,
    GeneSetRequest,
    GeneSetResponse,
    GroupCountOut,
    GroupsResponse,
    HistologyCountOut,
    JunctionMetadataLoaded,
    MetadataLoaded,
    ModelResponse,
    SelectGenesetRequest,
    SelectRequest,
    SelectResponse,
    SessionCreated,
    SessionState,
    SjdatLoaded,
    SjdatOption,
    StratifyRequest,
    SuggestAgeBandsRequest,
)
from ..services import model as model_mod
from ..services.metadata import (
    ALL_GROUPS,
    ALL_GROUPS_LABEL,
    CLASSIC_AGE_BANDS,
    DEFAULT_MIN_GROUP_N,
    AgeBands,
    MetadataError,
    StratifyError,
    build_covariates,
    parse_raw_metadata,
    quantile_age_bands,
    stratify,
)
from ..services.select import SelectError, select_features, select_from_features
from ..services.sessions import SJDAT_KINDS, store
from ..services.sjdat import SJDAT_META, SjdatError, load_sjdat

# Gene-set resolution reuses sjvc's services directly — SJSurv's active sjdat
# matrix is a `sjvc.services.sjdat.Sjdat`, a deliberate drop-in replacement
# for `sjvc.services.rds.Matrix` (see that module's docstring), so this code
# is byte-for-byte what these functions already do for 2D View.
from sjvc.services import features as feat_mod, geneset as gs_mod, mad as mad_mod
from sjvc.services.junction_metadata import JunctionMetadataError, load_junction_gene_index
from sjvc.services.junctions import looks_gene_level, map_junctions_to_genes
from sjvc.services import pathways

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


# --------------------------------------------------------------------------- #
# ingest (shared with the dataload package)
# --------------------------------------------------------------------------- #
def _sync_sjdat_intersection(s) -> None:
    """Narrow every loaded sjdat matrix to the samples that also have a
    sample-metadata row, recomputed from ``s.sjdat_raw`` (the pristine,
    as-loaded matrices) so loading a *different* metadata file later is never
    narrowed beyond what the current metadata actually covers. Mirrors
    ``sjvc.api.routes._sync_sample_intersection``. A sample missing from the
    metadata (no row at all — as opposed to a row with an unlabelled group)
    is dropped here; feature selection already scopes further, per group, to
    only the *labelled* samples of ``s.metadata``."""
    known = set(s.raw_metadata.rows) if s.raw_metadata is not None else None
    for kind, raw in s.sjdat_raw.items():
        s.sjdat[kind] = raw.restrict_samples(known) if known is not None else raw


def ingest_sjdat(s, kind: str, src: Path) -> SjdatLoaded:
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
    if s.active_sjdat is None or s.active_sjdat not in s.sjdat_raw:
        s.active_sjdat = kind
    # a new matrix invalidates any downstream work
    s.selection = s.model = s.labels = s.selected_group = None
    s.geneset = s.features = None
    _sync_sjdat_intersection(s)
    narrowed = s.sjdat[kind]

    warnings: List[str] = []
    if s.raw_metadata is not None:
        n_excluded = d.n_samples - narrowed.n_samples
        if narrowed.n_samples == 0:
            warnings.append(
                "none of this matrix's columns match a sample_id in the loaded metadata"
            )
        elif n_excluded:
            warnings.append(
                f"{n_excluded} matrix sample(s) excluded — no matching metadata row"
            )
    return SjdatLoaded(
        kind=kind, n_features=narrowed.n_features, n_samples=narrowed.n_samples,
        sparse=narrowed.sparse, warnings=warnings,
    )


def ingest_metadata(s, src: Path) -> MetadataLoaded:
    """Load the *raw* sample metadata (histology / stage / age / survival —
    the file prepTCGAdata::get_tcga_data() writes) and immediately stratify it
    with the classic default age bands, so the pipeline works out of the box.
    The Stratify panel then lets the user try different age bands without
    reloading the file."""
    known = _known_samples(s)
    try:
        raw = parse_raw_metadata(src, s.tmp_dir, known)
    except MetadataError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)
    s.raw_metadata = raw
    s.selected_covariates = [c.key for c in raw.available_covariates() if c.default]
    _sync_sjdat_intersection(s)

    warnings: List[str] = []
    if raw.n_unmatched:
        warnings.append(
            f"{raw.n_unmatched} metadata row(s) had a sample_id not in the sjdat matrix"
        )
    if not known:
        warnings.append("no sjdat matrix loaded yet — every metadata row was kept")
    for kind, raw_d in s.sjdat_raw.items():
        n_excluded = raw_d.n_samples - s.sjdat[kind].n_samples
        if n_excluded:
            warnings.append(
                f"{n_excluded} {SJDAT_META[kind][0].lower()} sample(s) excluded — "
                f"no matching metadata row"
            )

    _apply_stratification(s, CLASSIC_AGE_BANDS, DEFAULT_MIN_GROUP_N)
    return _metadata_loaded(s, warnings)


def _apply_stratification(
    s, age_bands: AgeBands, min_group_n: int,
    *, use_histology: bool = True, histology_map: dict | None = None,
) -> None:
    histology_map = histology_map or {}
    s.metadata = stratify(
        s.raw_metadata, age_bands, min_group_n,
        use_histology=use_histology, histology_map=histology_map,
    )
    s.age_bands = age_bands
    s.min_group_n = min_group_n
    s.use_histology = use_histology
    s.histology_map = histology_map
    # the Group / SurviverGroup definitions just changed under it
    s.selection = s.model = s.labels = s.selected_group = None


def _age_bands_out(bands: AgeBands) -> AgeBandsOut:
    return AgeBandsOut(
        edges=[None if e == float("inf") else e for e in bands.edges],
        include_lowest=bands.include_lowest,
        labels=bands.labels,
    )


def _metadata_loaded(s, warnings: List[str]) -> MetadataLoaded:
    raw = s.raw_metadata
    summary = raw.age_summary()
    return MetadataLoaded(
        n_matched=raw.n_rows, n_unmatched=raw.n_unmatched,
        n_with_age=summary["n_with_age"], age_min=summary["age_min"], age_max=summary["age_max"],
        age_bands=_age_bands_out(s.age_bands) if s.age_bands else None,
        min_group_n=s.min_group_n,
        histology_counts=[HistologyCountOut(value=v, n=n) for v, n in raw.histology_counts()],
        use_histology=s.use_histology,
        histology_map=s.histology_map,
        warnings=warnings,
    )


def _known_samples(s) -> List[str]:
    """Every sample any *pristine* (not yet metadata-narrowed) sjdat matrix
    has, so re-loading a metadata file always sees every candidate sample —
    never just whatever the previous metadata happened to narrow down to."""
    seen: List[str] = []
    for d in s.sjdat_raw.values():
        for c in d.samples:
            seen.append(c)
    return seen


def ingest_junction_metadata(s, src: Path) -> JunctionMetadataLoaded:
    """Load the per-junction annotation table (``TCGA_<cohort>_junction_
    metadata.rds``) — enables fast typed-gene / pathway lookups against a
    junction-level sjdat (junction counts, RRS scores) without needing a
    GENCODE release at all. Identical to sjvc's own `ingest_junction_metadata`
    — see ``sjvc/services/junction_metadata.py``."""
    try:
        idx = load_junction_gene_index(src, s.tmp_dir)
    except JunctionMetadataError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)
    s.junction_gene_index = idx
    s.geneset = None
    s.features = None
    return JunctionMetadataLoaded(n_rows=idx.n_rows, n_genes=len(idx.gene_id_to_name))


# --------------------------------------------------------------------------- #
# upload endpoints (the frontend loads via the Data tab; these are for parity)
# --------------------------------------------------------------------------- #
@router.post("/session/{sid}/sjdat/{kind}", response_model=SjdatLoaded)
async def upload_sjdat(sid: str, kind: str, file: UploadFile = File(...)) -> SjdatLoaded:
    s = _session(sid)
    dest = s.tmp_dir / f"sjdat_{kind}.rds"
    await _save_upload(file, dest)
    return ingest_sjdat(s, kind, dest)


@router.post("/session/{sid}/metadata", response_model=MetadataLoaded)
async def upload_metadata(sid: str, file: UploadFile = File(...)) -> MetadataLoaded:
    s = _session(sid)
    suffix = "".join(Path(file.filename or "metadata.rds").suffixes) or ".rds"
    dest = s.tmp_dir / f"metadata{suffix}"
    await _save_upload(file, dest)
    return ingest_metadata(s, dest)


@router.post("/session/{sid}/junction-metadata", response_model=JunctionMetadataLoaded)
async def upload_junction_metadata(sid: str, file: UploadFile = File(...)) -> JunctionMetadataLoaded:
    s = _session(sid)
    suffix = "".join(Path(file.filename or "junction_metadata.rds").suffixes) or ".rds"
    dest = s.tmp_dir / f"junction_metadata{suffix}"
    await _save_upload(file, dest)
    return ingest_junction_metadata(s, dest)


# --------------------------------------------------------------------------- #
# state / options
# --------------------------------------------------------------------------- #
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
    s = _session(sid)
    md = _metadata_loaded(s, []) if s.raw_metadata is not None else None
    return SessionState(
        has_raw_metadata=s.raw_metadata is not None,
        has_metadata=s.metadata is not None,
        metadata=md,
        sjdat_options=_sjdat_options(s),
        active_sjdat=s.active_sjdat,
        gencode_label=s.annotation.label if s.annotation is not None else None,
        pathways_enabled=s.species == "human" and s.annotation is not None,
        has_junction_metadata=s.junction_gene_index is not None,
        has_geneset=s.geneset is not None,
        has_features=s.features is not None,
        selected_group=s.selected_group,
        has_selection=s.selection is not None,
        has_model=s.model is not None,
        covariate_columns=(
            [CovariateColumnOut(**c.__dict__) for c in s.raw_metadata.available_covariates()]
            if s.raw_metadata is not None else []
        ),
        selected_covariates=s.selected_covariates,
    )


# --------------------------------------------------------------------------- #
# age-band stratification
# --------------------------------------------------------------------------- #
@router.post("/session/{sid}/age-bands/suggest", response_model=AgeBandsOut)
def suggest_age_bands(sid: str, body: SuggestAgeBandsRequest) -> AgeBandsOut:
    """A data-driven alternative to the classic fixed cut points: `n_bands`
    equal-count bands computed from this cohort's own age distribution."""
    s = _session(sid)
    if s.raw_metadata is None:
        raise HTTPException(409, "load the sample metadata on the Data tab first")
    try:
        bands = quantile_age_bands(s.raw_metadata.ages(), body.n_bands)
    except StratifyError as e:
        raise HTTPException(422, str(e))
    return _age_bands_out(bands)


@router.post("/session/{sid}/stratify", response_model=MetadataLoaded)
def stratify_route(sid: str, body: StratifyRequest) -> MetadataLoaded:
    """Recompute Histology/Stage/Age_at_diagnosis/Group/MedianSurvival/
    SurviverGroup with the given age bands and minimum group size. Invalidates
    any feature selection / model made under the previous stratification."""
    s = _session(sid)
    if s.raw_metadata is None:
        raise HTTPException(409, "load the sample metadata on the Data tab first")
    if body.min_group_n < 1:
        raise HTTPException(422, "min_group_n must be >= 1")

    edges = [float("inf") if e is None else float(e) for e in body.edges]
    try:
        bands = AgeBands(edges=edges, include_lowest=body.include_lowest)
    except StratifyError as e:
        raise HTTPException(422, str(e))

    _apply_stratification(
        s, bands, body.min_group_n,
        use_histology=body.use_histology, histology_map=body.histology_map,
    )
    return _metadata_loaded(s, [])


@router.post("/session/{sid}/covariates", response_model=SessionState)
def set_covariates(sid: str, body: CovariatesRequest) -> SessionState:
    """Choose which sample aspects from the loaded metadata file ride along
    as covariates at cross-validate/model time (see `build_covariates()`).
    Doesn't touch the molecular feature selection itself — only the saved
    MODEL, which would otherwise misreport what it was actually trained on."""
    s = _session(sid)
    if s.raw_metadata is None:
        raise HTTPException(409, "load the sample metadata on the Data tab first")
    valid = {c.key for c in s.raw_metadata.available_covariates()}
    unknown = [k for k in body.columns if k not in valid]
    if unknown:
        raise HTTPException(422, f"unknown covariate column(s): {', '.join(unknown)}")
    s.selected_covariates = list(dict.fromkeys(body.columns))
    s.model = None
    return session_state(sid)


@router.post("/session/{sid}/sjdat", response_model=SessionState)
def activate_sjdat(sid: str, body: ActivateSjdat) -> SessionState:
    s = _session(sid)
    if body.kind not in s.sjdat:
        raise HTTPException(409, f"the {body.kind!r} matrix has not been loaded")
    if s.active_sjdat != body.kind:
        s.active_sjdat = body.kind
        s.selection = s.model = s.labels = s.selected_group = None
        s.geneset = s.features = None
    return session_state(sid)


# --------------------------------------------------------------------------- #
# groups
# --------------------------------------------------------------------------- #
@router.get("/session/{sid}/groups", response_model=GroupsResponse)
def groups(sid: str) -> GroupsResponse:
    s = _session(sid)
    if s.metadata is None:
        raise HTTPException(409, "load the sample metadata on the Data tab first")
    d = s.sjdat.get(s.active_sjdat) if s.active_sjdat else None

    out = []
    for gc in s.metadata.group_counts():
        label = ALL_GROUPS_LABEL if gc.group == ALL_GROUPS else gc.group
        out.append(GroupCountOut(
            group=gc.group, label=label, n_total=gc.n_total,
            n_good=gc.n_good, n_poor=gc.n_poor, n_labelled=gc.n_labelled,
        ))

    warnings: List[str] = []
    if d is not None:
        in_matrix = sum(1 for sid_ in s.metadata.rows if sid_ in d._col)
        if in_matrix < s.metadata.n_rows:
            warnings.append(
                f"{s.metadata.n_rows - in_matrix} metadata sample(s) are not columns of the "
                f"active sjdat matrix and will be dropped from the analysis"
            )
    return GroupsResponse(groups=out, warnings=warnings)


# --------------------------------------------------------------------------- #
# gene-set feature resolution (Type genes / Upload list / Pathway tabs) — a
# byte-for-byte port of sjvc.api.routes's set_geneset()/build_features(), see
# that module's docstrings for the "why" of each branch. The only difference
# is the matrix source: `s.sjdat.get(s.active_sjdat)` (SJSurv's own active
# matrix) in place of sjvc's `s.junctions` view.
# --------------------------------------------------------------------------- #
def _active_sjdat_or_409(s):
    if not s.active_sjdat or s.active_sjdat not in s.sjdat:
        raise HTTPException(409, "load an sjdat matrix on the Data tab and pick one first")
    return s.sjdat[s.active_sjdat]


@router.post("/session/{sid}/geneset", response_model=GeneSetResponse)
def set_geneset(sid: str, body: GeneSetRequest) -> GeneSetResponse:
    s = _session(sid)
    d = _active_sjdat_or_409(s)
    gene_level = looks_gene_level(d.features)
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
            symbols, d.features, prefix=prefix
        )
        if not matched_names:
            raise HTTPException(422, f"none of the {len(symbols)} gene(s) are rows in this matrix")
        gs = gs_mod.GeneSet(
            matched_names=matched_names, unmatched=unmatched, source=source, warnings=warns
        )
        resp_matched = matched_names
    elif fast_lookup:
        # junction-level matrix, junction metadata loaded: a plain lookup,
        # no GENCODE release needed at all (see sjvc/services/junction_metadata.py)
        jgmap, labels, unmatched = s.junction_gene_index.resolve(symbols, prefix=prefix)
        if not labels:
            raise HTTPException(
                422, f"none of the {len(symbols)} gene(s) have a junction in the loaded "
                     f"junction metadata"
            )
        warns: List[str] = []
        n_dropped = sum(1 for rn in jgmap.rownames if rn not in d._row)
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
    d = _active_sjdat_or_409(s)
    if s.geneset is None:
        raise HTTPException(409, "choose a gene set first")

    # gene-level matrix: the gene set is a plain list of row labels to keep
    if getattr(s.geneset, "matched_names", None):
        try:
            fm = feat_mod.select_by_labels(d, s.geneset.matched_names)
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
                d, jgmap, condense=body.condense, gene_labels=s.geneset.matched_gene_labels,
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

    jgmap = map_junctions_to_genes(d.features, s.geneset.matched)
    if jgmap.n_junctions() == 0:
        raise HTTPException(
            422, "no junction in the matrix overlaps the selected genes (check genome build / naming)"
        )

    labels = {g.gene_id: g.name for g in s.geneset.matched}
    try:
        fm = feat_mod.select_and_maybe_condense(
            d, jgmap, condense=body.condense, gene_labels=labels
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


# --------------------------------------------------------------------------- #
# feature selection
# --------------------------------------------------------------------------- #
@router.post("/session/{sid}/select", response_model=SelectResponse)
def select(sid: str, body: SelectRequest) -> SelectResponse:
    s = _session(sid)
    if s.metadata is None:
        raise HTTPException(409, "load the sample metadata on the Data tab first")
    if not s.active_sjdat or s.active_sjdat not in s.sjdat:
        raise HTTPException(409, "load an sjdat matrix on the Data tab and pick one first")
    d = s.sjdat[s.active_sjdat]

    labels = s.metadata.labelled_samples(body.group)
    if not labels:
        raise HTTPException(422, f"no Good/Poor-labelled sample in group {body.group!r}")

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

    sample_ids = [sid_ for sid_ in labels if sid_ in d._col]
    try:
        sel = select_features(
            d, sample_ids, n_min=body.n_min, x_min=body.x_min, top_n=body.top_n,
            restrict_to=restrict_to, annotation=s.annotation,
        )
    except SelectError as e:
        raise HTTPException(422, str(e))

    labels = {sid_: labels[sid_] for sid_ in sel.sample_ids}
    n_good = sum(1 for v in labels.values() if v == "Good")
    n_poor = len(labels) - n_good

    s.selection = sel
    s.selected_group = body.group
    s.labels = labels
    s.model = None

    warnings: List[str] = []
    if restrict_to is not None:
        warnings.append(
            "restricted ranking to protein-coding genes"
            if looks_gene_level(d.features)
            else "restricted ranking to junctions overlapping a protein-coding gene in the reference"
        )
    dropped = len(labels) - len(sel.sample_ids) if len(labels) > len(sel.sample_ids) else 0
    if dropped:
        warnings.append(f"{dropped} labelled sample(s) were not columns of the sjdat matrix")
    if sel.n_features < body.top_n:
        warnings.append(
            f"only {sel.n_features} feature(s) passed the coverage filter (asked for top {body.top_n})"
        )
    if min(n_good, n_poor) < 2:
        warnings.append(
            f"the smaller class has only {min(n_good, n_poor)} sample(s) — cross-validation "
            f"needs at least 2 of each"
        )
    return SelectResponse(
        group=body.group, sjdat_kind=s.active_sjdat,
        n_group_samples=len(sel.sample_ids), n_good=n_good, n_poor=n_poor,
        n_candidates=sel.n_candidates, n_after_coverage=sel.n_after_coverage,
        n_selected=sel.n_features,
        feature_preview=sel.feature_ids[:15], warnings=warnings,
    )


@router.post("/session/{sid}/select-geneset", response_model=SelectResponse)
def select_geneset(sid: str, body: SelectGenesetRequest) -> SelectResponse:
    """The Type-genes/Upload-list/Pathway counterpart of `/select`: turns the
    already-built `s.features` (see `/features` above) into a Selection for
    the chosen Group — no MAD ranking, no N/X coverage filter, the resolved
    gene set *is* the selection."""
    s = _session(sid)
    if s.metadata is None:
        raise HTTPException(409, "load the sample metadata on the Data tab first")
    if s.features is None:
        raise HTTPException(409, "resolve a gene set and build features first")
    if not s.active_sjdat:
        raise HTTPException(409, "load an sjdat matrix on the Data tab and pick one first")

    labels = s.metadata.labelled_samples(body.group)
    if not labels:
        raise HTTPException(422, f"no Good/Poor-labelled sample in group {body.group!r}")

    try:
        sel = select_from_features(s.features, list(labels), s.active_sjdat)
    except SelectError as e:
        raise HTTPException(422, str(e))

    labels = {sid_: labels[sid_] for sid_ in sel.sample_ids}
    n_good = sum(1 for v in labels.values() if v == "Good")
    n_poor = len(labels) - n_good

    s.selection = sel
    s.selected_group = body.group
    s.labels = labels
    s.model = None

    warnings: List[str] = []
    if min(n_good, n_poor) < 2:
        warnings.append(
            f"the smaller class has only {min(n_good, n_poor)} sample(s) — cross-validation "
            f"needs at least 2 of each"
        )
    return SelectResponse(
        group=body.group, sjdat_kind=s.active_sjdat,
        n_group_samples=len(sel.sample_ids), n_good=n_good, n_poor=n_poor,
        n_candidates=sel.n_candidates, n_after_coverage=sel.n_after_coverage,
        n_selected=sel.n_features,
        feature_preview=sel.feature_ids[:15], warnings=warnings,
    )


def _covariates_for(s, sample_ids: List[str]):
    if s.raw_metadata is None or not s.selected_covariates:
        return [], None
    return build_covariates(
        s.raw_metadata, s.selected_covariates, sample_ids, histology_map=s.histology_map,
    )


@router.post("/session/{sid}/cross-validate", response_model=CVResponse)
def cross_validate(sid: str, body: CVRequest) -> CVResponse:
    s = _session(sid)
    if s.selection is None or s.labels is None:
        raise HTTPException(409, "select features first")
    n_cv = body.n_cv or 5
    cov_names, cov_values = _covariates_for(s, s.selection.sample_ids)
    try:
        r = model_mod.cross_validate(
            s.selection, s.labels, int(n_cv), cov_names=cov_names, cov_values=cov_values,
        )
    except model_mod.ModelError as e:
        raise HTTPException(422, str(e))
    return CVResponse(**r.__dict__)


@router.post("/session/{sid}/model", response_model=ModelResponse)
def train_model(sid: str) -> ModelResponse:
    s = _session(sid)
    if s.selection is None or s.labels is None:
        raise HTTPException(409, "select features first")
    cov_names, cov_values = _covariates_for(s, s.selection.sample_ids)
    try:
        m = model_mod.train_full(
            s.selection, s.labels,
            sjdat_kind=s.active_sjdat, group=s.selected_group or "",
            cov_names=cov_names, cov_values=cov_values,
        )
    except model_mod.ModelError as e:
        raise HTTPException(422, str(e))
    s.model = m

    better = m.auc_resub > 0.5
    messages = [
        f"Trained on all {m.n_samples} samples (Good={m.n_good}, Poor={m.n_poor}); "
        f"tested on the same samples (resubstitution — optimistic).",
        f"Resubstitution AUC {m.auc_resub:.3f}, accuracy {m.accuracy_resub:.3f}.",
        "MODEL saved to this session for further examination.",
    ]
    if not better:
        messages.append(
            "Resubstitution AUC is not above 0.5 — these features don't separate the "
            "survivor groups; try another sjdat, group, or looser feature filters."
        )
    return ModelResponse(
        sjdat_kind=m.sjdat_kind, group=m.group,
        n_samples=m.n_samples, n_good=m.n_good, n_poor=m.n_poor, n_features=m.n_features,
        auc_resub=m.auc_resub, accuracy_resub=m.accuracy_resub,
        confusion_resub=m.confusion_resub,
        features=[FeatureWeightOut(**fw.__dict__) for fw in m.features],
        saved=True, messages=messages,
    )
