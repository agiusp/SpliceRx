"""HTTP API for SJV."""
from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, File, HTTPException, UploadFile

from ..models import (
    ArcModel,
    ExonModel,
    GencodeReady,
    GencodeSelect,
    GeneLookupResponse,
    GeneRecord,
    GeneSuggestResponse,
    JunctionMetadataLoaded,
    LegendEntry,
    PlotLayout,
    PlotRequest,
    PlotResponse,
    RdsUploaded,
    ReleasesResponse,
    SampleMetadataUploaded,
    SessionCreated,
    SessionState,
    StratColumn,
    StratValue,
    StratValuesResponse,
    TranscriptModel,
)
from ..palette import legend_entries
from ..services import gencode
from ..services.classify import GeneClassifier
from ..services.groups import parse_sample_metadata
from ..services.junctions import norm_chrom, scale_counts
from ..services.rds import RdsError, load_rds
from ..services.sessions import store
from sjlookup.services.lookup import (
    JunctionMetadataError as LookupMetadataError,
    classify_unannotated,
    load_junction_lookup_index,
)

router = APIRouter()

_UPLOAD_CAP = 2 * 1024 * 1024 * 1024  # 2 GiB


def _session(session_id: str):
    s = store.get(session_id)
    if s is None:
        raise HTTPException(404, "unknown or expired session")
    return s


@router.post("/session", response_model=SessionCreated)
def create_session() -> SessionCreated:
    return SessionCreated(session_id=store.create().id)


def ingest_rds(s, src: Path) -> RdsUploaded:
    """Parse the junction matrix at ``src`` (a path inside ``s.tmp_dir``) into the
    session. ``src`` is consumed. Shared by the browser-upload endpoint and the
    dataload package. A large sparse matrix takes ~1 min the first time and is
    then cached under ``$SJV_CACHE_DIR/matrices``."""
    try:
        matrix = load_rds(src, s.tmp_dir)
    except RdsError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)  # the raw matrix is never kept

    s.rds = matrix
    warnings = []
    if matrix.n_bad_rownames:
        warnings.append(
            f"{matrix.n_bad_rownames} row name(s) did not match chr:start-end:strand "
            f"and were skipped (e.g. {', '.join(matrix.bad_rownames_sample[:3])})"
        )
    return RdsUploaded(
        samples=matrix.samples,
        n_junctions=matrix.n_junctions,
        n_unstranded=matrix.n_bad_rownames,
        sparse=matrix.sparse,
        warnings=warnings,
    )


@router.post("/session/{session_id}/rds", response_model=RdsUploaded)
async def upload_rds(session_id: str, file: UploadFile = File(...)) -> RdsUploaded:
    s = _session(session_id)
    dest = s.tmp_dir / "upload.rds"
    size = 0
    with open(dest, "wb") as fh:
        while chunk := await file.read(1 << 20):
            size += len(chunk)
            if size > _UPLOAD_CAP:
                fh.close()
                dest.unlink(missing_ok=True)
                raise HTTPException(413, "RDS file exceeds the 2 GiB limit")
            fh.write(chunk)
    return ingest_rds(s, dest)


@router.get("/session/{session_id}/state", response_model=SessionState)
def session_state(session_id: str) -> SessionState:
    """Everything the frontend needs to rebuild step 1 from a session it did not
    populate itself (the Data tab hand-off) or after a page refresh."""
    s = _session(session_id)
    return SessionState(
        has_rds=s.rds is not None,
        n_junctions=s.rds.n_junctions if s.rds is not None else 0,
        samples=s.rds.samples if s.rds is not None else [],
        sparse=bool(s.rds and s.rds.sparse),
        sample_metadata_columns=list(s.metadata.columns) if s.metadata is not None else [],
        gencode_label=s.annotation.label if s.annotation is not None else None,
        has_junction_metadata=s.junction_lookup is not None,
        junction_metadata_has_detail=bool(s.junction_lookup and s.junction_lookup.has_annotation_detail),
    )


@router.get("/gencode/releases", response_model=ReleasesResponse)
def releases() -> ReleasesResponse:
    return ReleasesResponse(releases=gencode.list_releases())


@router.post("/session/{session_id}/gencode", response_model=GencodeReady)
def select_gencode(session_id: str, body: GencodeSelect) -> GencodeReady:
    s = _session(session_id)
    try:
        ann = gencode.ensure_release(body.species, body.release)
    except gencode.GencodeError as e:
        raise HTTPException(502, str(e))
    s.annotation = ann
    return GencodeReady(label=ann.label, warnings=_naming_warnings(s))


@router.post("/session/{session_id}/gtf", response_model=GencodeReady)
async def upload_gtf(session_id: str, file: UploadFile = File(...)) -> GencodeReady:
    s = _session(session_id)
    suffix = ".gtf.gz" if file.filename and file.filename.endswith(".gz") else ".gtf"
    dest = s.tmp_dir / f"annotation{suffix}"
    with open(dest, "wb") as fh:
        while chunk := await file.read(1 << 20):
            fh.write(chunk)
    try:
        ann = gencode.annotation_from_gtf(dest, label=file.filename or "custom GTF")
    except gencode.GencodeError as e:
        raise HTTPException(422, str(e))
    s.annotation = ann
    return GencodeReady(label=ann.label, warnings=_naming_warnings(s))


@router.get("/session/{session_id}/gene", response_model=GeneLookupResponse)
def lookup_gene(session_id: str, name: str) -> GeneLookupResponse:
    s = _session(session_id)
    if s.annotation is None:
        raise HTTPException(409, "select a GENCODE release or upload a GTF first")
    gene = s.annotation.get_gene(name)
    if gene is None:
        return GeneLookupResponse(gene=None, near_matches=s.annotation.near_matches(name))
    return GeneLookupResponse(gene=GeneRecord(**gene.__dict__))


@router.get("/session/{session_id}/genes/suggest", response_model=GeneSuggestResponse)
def suggest_genes(session_id: str, q: str, limit: int = 20) -> GeneSuggestResponse:
    """Autocomplete for the query-gene box: gene names starting with `q`."""
    s = _session(session_id)
    if s.annotation is None:
        return GeneSuggestResponse(names=[])
    return GeneSuggestResponse(names=s.annotation.suggest(q, min(limit, 50)))


def ingest_sample_metadata(s, src: Path) -> SampleMetadataUploaded:
    """Parse the sample-metadata table at ``src`` (a path inside ``s.tmp_dir``,
    extension preserved) into the session. ``src`` is consumed."""
    if s.rds is None:
        raise HTTPException(409, "load the junction matrix first")
    try:
        meta = parse_sample_metadata(src, s.tmp_dir, s.rds.samples)
    except RdsError as e:
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)

    s.metadata = meta
    warnings = []
    if meta.unmatched_samples:
        warnings.append(
            f"{len(meta.unmatched_samples)} sample_id value(s) in the table are not "
            f"columns of the matrix and were ignored (e.g. {', '.join(meta.unmatched_samples[:3])})"
        )
    return SampleMetadataUploaded(
        columns=[StratColumn(**c) for c in meta.column_summaries()],
        n_matched=meta.n_matched,
        warnings=warnings,
    )


@router.post("/session/{session_id}/sample-metadata", response_model=SampleMetadataUploaded)
async def upload_sample_metadata(
    session_id: str, file: UploadFile = File(...)
) -> SampleMetadataUploaded:
    s = _session(session_id)
    if s.rds is None:
        raise HTTPException(409, "upload an RDS file first")
    suffix = "".join(Path(file.filename or "metadata.csv").suffixes) or ".csv"
    dest = s.tmp_dir / f"metadata{suffix}"
    with open(dest, "wb") as fh:
        while chunk := await file.read(1 << 20):
            fh.write(chunk)
    return ingest_sample_metadata(s, dest)


def ingest_junction_metadata(s, src: Path) -> JunctionMetadataLoaded:
    """Load the cohort's per-junction annotation table (``TCGA_<cohort>_
    junction_metadata.rds``) so the sashimi plot can classify arcs from its
    ``annotated``/``left_annotated``/``right_annotated`` columns — the
    recount3/STAR-aligner annotation recorded when the cohort's junctions
    were originally called — as an alternative to live GENCODE transcript
    matching. See ``PlotRequest.annotation_source`` and
    ``sjlookup.services.lookup``, which this reuses rather than
    re-implementing."""
    try:
        idx = load_junction_lookup_index(src, s.tmp_dir)
    except LookupMetadataError as e:
        src.unlink(missing_ok=True)
        raise HTTPException(422, str(e))
    finally:
        src.unlink(missing_ok=True)
    s.junction_lookup = idx
    warnings = []
    if idx.n_duplicate_rownames:
        warnings.append(
            f"{idx.n_duplicate_rownames} duplicate junction row(s) in the table — the first of "
            f"each was kept"
        )
    if not idx.has_annotation_detail:
        warnings.append(
            "no left_annotated/right_annotated columns — a junction the table doesn't mark "
            "annotated will just read \"novel\" (no exon_skipping/alt_5p/alt_3p/novel_exon split)"
        )
    return JunctionMetadataLoaded(
        n_rows=idx.n_rows, n_duplicate_rownames=idx.n_duplicate_rownames,
        has_annotation_detail=idx.has_annotation_detail, warnings=warnings,
    )


@router.post("/session/{session_id}/junction-metadata", response_model=JunctionMetadataLoaded)
async def upload_junction_metadata(session_id: str, file: UploadFile = File(...)) -> JunctionMetadataLoaded:
    s = _session(session_id)
    suffix = "".join(Path(file.filename or "junction_metadata.rds").suffixes) or ".rds"
    dest = s.tmp_dir / f"junction_metadata{suffix}"
    with open(dest, "wb") as fh:
        while chunk := await file.read(1 << 20):
            fh.write(chunk)
    return ingest_junction_metadata(s, dest)


@router.get(
    "/session/{session_id}/sample-metadata/values", response_model=StratValuesResponse
)
def strat_values(session_id: str, column: str) -> StratValuesResponse:
    s = _session(session_id)
    if s.metadata is None:
        raise HTTPException(409, "upload a sample-metadata table first")
    vals = s.metadata.strat_values(column)
    if not vals and column not in s.metadata.columns:
        raise HTTPException(404, f"column {column!r} is not in the uploaded table")
    return StratValuesResponse(
        column=column,
        values=[
            StratValue(value=v, n_samples=len(ids))
            for v, ids in sorted(vals.items(), key=lambda kv: kv[0])
        ],
    )


@router.post("/session/{session_id}/plot", response_model=PlotResponse)
def plot(session_id: str, body: PlotRequest) -> PlotResponse:
    s = _session(session_id)
    if s.rds is None:
        raise HTTPException(409, "upload an RDS file first")
    if s.annotation is None:
        raise HTTPException(409, "select a GENCODE release or upload a GTF first")
    if body.annotation_source == "metadata" and s.junction_lookup is None:
        raise HTTPException(
            409,
            "load the cohort's junction metadata table (TCGA_<cohort>_junction_metadata.rds) on "
            "the Data tab first, or switch the annotation source back to \"computed live\"",
        )

    # resolve every requested gene name
    genes = []
    seen_ids: set = set()
    for name in body.genes:
        g = s.annotation.get_gene(name)
        if g is None:
            raise HTTPException(
                404,
                {
                    "message": f"gene {name!r} not found in {s.annotation.label}",
                    "near_matches": s.annotation.near_matches(name),
                },
            )
        if g.gene_id not in seen_ids:
            seen_ids.add(g.gene_id)
            genes.append(g)

    if len({norm_chrom(g.chrom) for g in genes}) > 1:
        raise HTTPException(
            422,
            "the genes are on different chromosomes ("
            + ", ".join(sorted({g.chrom for g in genes}))
            + ") — plot them one chromosome at a time",
        )

    # row indices (into the full matrix) of the junctions overlapping any gene —
    # each gene filters on its own strand; dedupe, keep matrix order. Only this
    # handful of rows is ever densified, so a multi-million-row sparse matrix is
    # as cheap as a small one.
    import numpy as np

    rows = np.unique(np.concatenate(
        [s.rds.locus_row_indices(g) for g in genes] or [np.empty(0, dtype=np.int64)]
    )).astype(np.int64)
    in_locus = s.rds.junctions_at(rows)

    try:
        if body.strat_column:
            if s.metadata is None:
                raise HTTPException(409, "upload a sample-metadata table first")
            ids = s.metadata.strat_values(body.strat_column).get(body.strat_value or "")
            if not ids:
                raise HTTPException(
                    404,
                    f"{body.strat_column}={body.strat_value!r} has no matching samples",
                )
            locus_counts = s.rds.group_median_at(rows, ids)
            series_label = f"{body.strat_column} = {body.strat_value} (median of {len(ids)})"
            count_kind = "median"
        else:
            locus_counts = s.rds.counts_at(rows, body.sample)
            series_label = body.sample or ""
            count_kind = "count"
    except RdsError as e:
        raise HTTPException(422, str(e))

    arcs_scaled = scale_counts(in_locus, list(locus_counts), min_reads=body.min_reads)

    # one classifier over every gene's transcripts, so a junction that connects
    # sites in different genes is still checked against all annotated sites
    all_transcripts = []
    tx_models = []
    for g in genes:
        txs = s.annotation.get_transcripts(g.gene_id)
        all_transcripts.extend(txs)
        tx_models.extend(
            TranscriptModel(
                transcript_id=t.transcript_id,
                strand=t.strand,
                exons=[ExonModel(start=e.start, end=e.end, kind=e.kind) for e in t.exons],
                gene_name=g.name,
            )
            for t in txs
        )
    classifier = GeneClassifier(all_transcripts, tol=0)

    n_metadata_fallback = 0
    arcs = []
    for a in arcs_scaled:
        category, source = classifier.classify(a.junction), "gencode"
        if body.annotation_source == "metadata":
            row = s.junction_lookup.lookup(a.junction.rowname)
            if row is None:
                n_metadata_fallback += 1  # not in the table — keep the GENCODE fallback above
            elif row.get("annotated"):
                category, source = "annotated", "metadata"
            else:
                sub = classify_unannotated(s.junction_lookup, row)
                category, source = (sub or "novel"), "metadata"
        arcs.append(
            ArcModel(
                id=a.junction.rowname,
                start=a.junction.start,
                end=a.junction.end,
                strand=a.junction.strand,
                count=a.count,
                height=a.height,
                category=category,
                category_source=source,
            )
        )

    lo = min(g.start for g in genes)
    hi = max(g.end for g in genes)
    pad = max(50, int((hi - lo) * 0.05))
    x_domain = [lo - pad, hi + pad]

    warnings = _naming_warnings(s)
    if not arcs:
        kind = "group (after excluding NA/0)" if body.strat_column else "sample"
        cut = f" with >= {body.min_reads:g} reads" if body.min_reads > 0 else ""
        where = "these genes" if len(genes) > 1 else "this gene"
        warnings.append(f"no junctions{cut} fall within {where} for this {kind}")
    if n_metadata_fallback:
        warnings.append(
            f"{n_metadata_fallback} of {len(arcs)} junction(s) have no row in the loaded "
            "junction metadata table and were classified from the live GENCODE reference instead"
        )

    return PlotResponse(
        genes=[GeneRecord(**g.__dict__) for g in genes],
        x_domain=x_domain,
        layout=PlotLayout(),
        transcripts=tx_models,
        arcs=arcs,
        legend=[LegendEntry(**e) for e in legend_entries([a.category for a in arcs])],
        series_label=series_label,
        annotation_source=body.annotation_source,
        count_kind=count_kind,
        warnings=warnings,
    )


def _naming_warnings(s) -> list:
    """Warn if RDS chromosome naming disagrees with the annotation's."""
    out: list = []
    if s.rds is None or s.annotation is None:
        return out
    rds_has_chr = s.rds.chrom_uses_chr_prefix()
    gene = None
    # cheap probe: look at one gene row
    con = s.annotation._con()
    try:
        row = con.execute("SELECT chrom FROM genes LIMIT 1").fetchone()
        gene = row["chrom"] if row else None
    finally:
        con.close()
    if gene is not None:
        ann_has_chr = gene.lower().startswith("chr")
        if rds_has_chr != ann_has_chr:
            out.append(
                "chromosome naming differs between the RDS "
                f"({'chr1' if rds_has_chr else '1'} style) and the annotation "
                f"({'chr1' if ann_has_chr else '1'} style); names are normalised "
                "when matching but check the genome build matches"
            )
    return out
