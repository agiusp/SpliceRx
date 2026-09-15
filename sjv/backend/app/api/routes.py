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
    LegendEntry,
    PlotLayout,
    PlotRequest,
    PlotResponse,
    RdsUploaded,
    ReleasesResponse,
    SampleMetadataUploaded,
    SessionCreated,
    StratColumn,
    StratValue,
    StratValuesResponse,
    TranscriptModel,
)
from ..palette import legend_entries
from ..services import gencode
from ..services.classify import GeneClassifier
from ..services.groups import parse_sample_metadata
from ..services.junctions import filter_for_gene, norm_chrom, scale_counts
from ..services.rds import RdsError, load_rds
from ..services.sessions import store

router = APIRouter(prefix="/api")

_UPLOAD_CAP = 2 * 1024 * 1024 * 1024  # 2 GiB


def _session(session_id: str):
    s = store.get(session_id)
    if s is None:
        raise HTTPException(404, "unknown or expired session")
    return s


@router.post("/session", response_model=SessionCreated)
def create_session() -> SessionCreated:
    return SessionCreated(session_id=store.create().id)


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

    try:
        matrix = load_rds(dest, s.tmp_dir)
    except RdsError as e:
        raise HTTPException(422, str(e))
    finally:
        dest.unlink(missing_ok=True)  # raw upload is never kept

    s.rds = matrix
    warnings = []
    if matrix.n_bad_rownames:
        warnings.append(
            f"{matrix.n_bad_rownames} row name(s) did not match chr:start-end:strand "
            f"and were skipped (e.g. {', '.join(matrix.bad_rownames_sample[:3])})"
        )
    return RdsUploaded(samples=matrix.samples, n_junctions=matrix.n_junctions, warnings=warnings)


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
    try:
        meta = parse_sample_metadata(dest, s.tmp_dir, s.rds.samples)
    except RdsError as e:
        raise HTTPException(422, str(e))
    finally:
        dest.unlink(missing_ok=True)

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
            counts = s.rds.group_median(ids)
            series_label = f"{body.strat_column} = {body.strat_value} (median of {len(ids)})"
            count_kind = "median"
        else:
            counts = s.rds.column(body.sample)
            series_label = body.sample or ""
            count_kind = "count"
    except RdsError as e:
        raise HTTPException(422, str(e))

    # union of the junctions overlapping any gene (each gene filters on its own
    # strand); dedupe on row name, keep matrix order
    idx = {j.rowname: i for i, j in enumerate(s.rds.junctions)}
    in_locus = []
    picked: set = set()
    for g in genes:
        for j in filter_for_gene(s.rds.junctions, g):
            if j.rowname not in picked:
                picked.add(j.rowname)
                in_locus.append(j)
    in_locus.sort(key=lambda j: idx[j.rowname])

    locus_counts = [counts[idx[j.rowname]] for j in in_locus]
    arcs_scaled = scale_counts(in_locus, locus_counts, min_reads=body.min_reads)

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

    arcs = [
        ArcModel(
            id=a.junction.rowname,
            start=a.junction.start,
            end=a.junction.end,
            strand=a.junction.strand,
            count=a.count,
            height=a.height,
            category=classifier.classify(a.junction),
        )
        for a in arcs_scaled
    ]

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

    return PlotResponse(
        genes=[GeneRecord(**g.__dict__) for g in genes],
        x_domain=x_domain,
        layout=PlotLayout(),
        transcripts=tx_models,
        arcs=arcs,
        legend=[LegendEntry(**e) for e in legend_entries([a.category for a in arcs])],
        series_label=series_label,
        count_kind=count_kind,
        warnings=warnings,
    )


def _naming_warnings(s) -> list:
    """Warn if RDS chromosome naming disagrees with the annotation's."""
    out: list = []
    if s.rds is None or s.annotation is None:
        return out
    rds_has_chr = any(j.chrom.lower().startswith("chr") for j in s.rds.junctions[:50])
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
