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
    ClinicalColumn,
    ClinicalUploaded,
    FeaturesRequest,
    FeaturesResponse,
    GencodeReady,
    GencodeSelect,
    GeneSetRequest,
    GeneSetResponse,
    HeatmapRequest,
    MadFeaturesRequest,
    MatrixUploaded,
    ProjectionRequest,
    SessionCreated,
)
from ..services import encoding, features as feat_mod, gencode, geneset as gs_mod, heatmap as hm_mod, mad as mad_mod, pathways, projection as proj_mod
from ..services.clinical import parse_clinical
from ..services.junctions import looks_gene_level, map_junctions_to_genes
from ..services.rds import MatrixError, load_matrix
from ..services.sessions import store

router = APIRouter(prefix="/api")
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


@router.post("/session/{sid}/junctions", response_model=MatrixUploaded)
async def upload_junctions(sid: str, file: UploadFile = File(...)) -> MatrixUploaded:
    s = _session(sid)
    dest = s.tmp_dir / "junctions.rds"
    await _save_upload(file, dest)
    try:
        m = load_matrix(dest, s.tmp_dir)
    except MatrixError as e:
        raise HTTPException(422, str(e))
    finally:
        dest.unlink(missing_ok=True)
    s.junctions = m
    s.features = None
    kind = "gene" if looks_gene_level(m.features) else "junction"
    warnings: List[str] = []
    if kind == "gene":
        warnings.append(
            "row names aren't chr:start-end:strand — treating this as an already-condensed, "
            "gene-level matrix. Use 'Top by MAD' to pick features from it."
        )
    return MatrixUploaded(
        samples=m.samples, n_junctions=m.n_features, feature_kind=kind, warnings=warnings
    )


@router.post("/session/{sid}/clinical", response_model=ClinicalUploaded)
async def upload_clinical(sid: str, file: UploadFile = File(...)) -> ClinicalUploaded:
    s = _session(sid)
    if s.junctions is None:
        raise HTTPException(409, "upload the junction matrix first")
    suffix = "".join(Path(file.filename or "clinical.csv").suffixes) or ".csv"
    dest = s.tmp_dir / f"clinical{suffix}"
    await _save_upload(file, dest)
    try:
        clin = parse_clinical(dest, s.tmp_dir, s.junctions.samples)
    except MatrixError as e:
        raise HTTPException(422, str(e))
    finally:
        dest.unlink(missing_ok=True)
    s.clinical = clin
    warnings: List[str] = []
    if clin.unmatched_samples:
        warnings.append(
            f"{len(clin.unmatched_samples)} sample_id value(s) not in the junction matrix "
            f"were ignored (e.g. {', '.join(clin.unmatched_samples[:3])})"
        )
    return ClinicalUploaded(
        columns=[ClinicalColumn(**c.__dict__) for c in clin.columns],
        n_matched=clin.n_rows,
        warnings=warnings,
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
        raise HTTPException(409, "upload the junction matrix first")
    gene_level = looks_gene_level(s.junctions.features)
    if not gene_level and s.annotation is None:
        raise HTTPException(409, "select a GENCODE release or upload a GTF first")

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
        s.features = fm
        return FeaturesResponse(
            feature_kind="gene",
            n_features=fm.n_features,
            n_samples=fm.n_samples,
            n_junctions=0,
            dropped_zero_variance=0,
            suggest_condense=False,
            feature_preview=[fm.label(i) for i in fm.feature_ids[:20]],
            warnings=[],
        )

    if s.annotation is None:
        raise HTTPException(409, "select a GENCODE release or upload a GTF first")

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
    s.features = fm

    warnings: List[str] = []
    if jgmap.n_bad_rownames:
        warnings.append(f"{jgmap.n_bad_rownames} matrix row name(s) were not chr:start-end:strand and skipped")
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
        if not looks_gene_level(s.junctions.features):
            raise HTTPException(
                422, "the protein-coding filter applies to gene-level matrices only"
            )
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
            fm = mad_mod.top_features_by_mad(s.junctions, body.top_n, restrict_to=restrict_to)
    except ValueError as e:
        raise HTTPException(422, str(e))

    s.features = fm
    s.geneset = None  # this selection mode doesn't go through a resolved gene set

    warnings: List[str] = []
    if restrict_to is not None:
        n_coding = sum(1 for f in s.junctions.features if f.strip().lower() in restrict_to)
        warnings.append(
            f"{n_coding} of {len(s.junctions.features)} matrix genes are protein-coding "
            f"in the reference — ranked those"
        )
    if fm.n_features < body.top_n:
        warnings.append(f"only {fm.n_features} feature(s) were available (asked for {body.top_n})")
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

    group_values = None
    if body.order == "group":
        if not body.group_by or s.clinical is None:
            raise HTTPException(422, "group order needs `group_by` and a clinical table")
        gv, _ = s.clinical.values_for(body.group_by, s.features.samples)
        group_values = [None if v is None else str(v) for v in gv]

    try:
        res = hm_mod.build(
            s.features, row_zscore=body.row_zscore, order=body.order, group_values=group_values,
        )
    except ValueError as e:
        raise HTTPException(422, str(e))

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

    return {
        "feature_kind": s.features.kind,
        "row_labels": [s.features.label(i) for i in res.feature_ids],
        "row_ids": res.feature_ids,
        "samples": res.samples,
        "values": [[float(x) for x in row] for row in res.values],
        "row_zscore": body.row_zscore,
        "row_dendro": res.row_dendro.segments if res.row_dendro else None,
        "col_dendro": res.col_dendro.segments if res.col_dendro else None,
        "annotations": tracks,
        "warnings": res.warnings,
    }
