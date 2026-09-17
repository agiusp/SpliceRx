"""Pydantic response/request schemas for the API."""
from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import BaseModel, model_validator


class SessionCreated(BaseModel):
    session_id: str


class RdsUploaded(BaseModel):
    samples: List[str]
    n_junctions: int
    n_unstranded: int = 0     # rows dropped because the name wasn't chr:start-end:strand
    sparse: bool = False      # read via the sparse (Matrix) path
    warnings: List[str] = []


class SessionState(BaseModel):
    """Reconstructed step-1 state of a session — for the Data tab hand-off and to
    survive a browser refresh."""
    has_rds: bool
    n_junctions: int = 0
    samples: List[str] = []
    sparse: bool = False
    sample_metadata_columns: List[str] = []
    gencode_label: Optional[str] = None
    has_junction_metadata: bool = False
    junction_metadata_has_detail: bool = False


class ReleasesResponse(BaseModel):
    releases: dict  # species -> [release, ...]


class GencodeSelect(BaseModel):
    species: str
    release: str


class GencodeReady(BaseModel):
    label: str
    warnings: List[str] = []


class GeneRecord(BaseModel):
    name: str
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str


class GeneLookupResponse(BaseModel):
    gene: Optional[GeneRecord] = None
    near_matches: List[str] = []


class GeneSuggestResponse(BaseModel):
    names: List[str]


class StratColumn(BaseModel):
    name: str
    n_values: int


class SampleMetadataUploaded(BaseModel):
    columns: List[StratColumn]      # candidate stratification columns
    n_matched: int                  # sample_id rows that matched the matrix
    warnings: List[str] = []


class StratValue(BaseModel):
    value: str
    n_samples: int


class StratValuesResponse(BaseModel):
    column: str
    values: List[StratValue]


class JunctionMetadataLoaded(BaseModel):
    n_rows: int
    n_duplicate_rownames: int
    has_annotation_detail: bool  # left_annotated/right_annotated columns were present —
    # without them, a junction the table doesn't call "annotated" can't be split into
    # exon_skipping/alt_5p/alt_3p/novel_exon/novel and just reads "novel"
    warnings: List[str] = []


class PlotRequest(BaseModel):
    genes: List[str] = []           # one or more gene names; the plot spans them all
    gene_name: Optional[str] = None  # legacy single-gene field, still accepted
    sample: Optional[str] = None
    strat_column: Optional[str] = None
    strat_value: Optional[str] = None
    min_reads: float = 0.0  # a junction is drawn only if its count (or group
    #                         median) is >= this; 0 draws every non-zero junction
    # "gencode" (default): classify every arc live against the loaded GENCODE
    # transcript models for the query gene(s) — see services/classify.py.
    # "metadata": prefer the cohort's own junction-metadata table (the
    # recount3/STAR-aligner `annotated` column prepTCGAdata carries alongside
    # the junction counts) when it's loaded, falling back to the GENCODE
    # classification for any arc the table has no row for.
    annotation_source: Literal["gencode", "metadata"] = "gencode"

    @model_validator(mode="after")
    def _validate(self) -> "PlotRequest":
        if not self.genes and self.gene_name:
            self.genes = [self.gene_name]
        # de-dupe, keep order, drop blanks
        seen: set = set()
        self.genes = [
            g for g in (n.strip() for n in self.genes)
            if g and not (g.lower() in seen or seen.add(g.lower()))
        ]
        if not self.genes:
            raise ValueError("provide at least one gene in `genes`")
        group_mode = bool(self.strat_column) and bool(self.strat_value)
        if bool(self.sample) == group_mode:
            raise ValueError(
                "provide either `sample`, or both `strat_column` and `strat_value`"
            )
        if self.min_reads < 0:
            raise ValueError("min_reads must be >= 0")
        return self


class ExonModel(BaseModel):
    start: int
    end: int
    kind: str  # "CDS" | "UTR"


class TranscriptModel(BaseModel):
    transcript_id: str
    strand: str
    exons: List[ExonModel]
    gene_name: str = ""


class ArcModel(BaseModel):
    id: str
    start: int
    end: int
    strand: str
    count: float
    height: float
    category: str
    # which source actually produced `category` for this arc — normally
    # matches the request's `annotation_source`, except in "metadata" mode
    # when this junction has no row in the loaded table and the live GENCODE
    # classification was used as a fallback instead
    category_source: Literal["gencode", "metadata"] = "gencode"


class LegendEntry(BaseModel):
    category: str
    label: str
    color: str


class PlotLayout(BaseModel):
    baseline_fraction: float = 0.55
    track_row_px: int = 16


class PlotResponse(BaseModel):
    genes: List[GeneRecord]      # every gene the plot spans, in request order
    x_domain: List[int]
    layout: PlotLayout
    transcripts: List[TranscriptModel]
    arcs: List[ArcModel]
    legend: List[LegendEntry]
    series_label: str = ""       # sample id, or "<group> (median of N)"
    count_kind: str = "count"    # "count" or "median"
    annotation_source: Literal["gencode", "metadata"] = "gencode"
    warnings: List[str] = []
