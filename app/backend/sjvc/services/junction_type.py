"""Per-junction "<gene name>:<novel-splicing-event type>" labels, for the
junction-level heatmap's alternate row-label view (see ``sjvc/api/routes.py``'s
``heatmap`` endpoint) — an alternative to the raw ``chr:start-end:strand``
row id, built from the same cohort's junction metadata table used elsewhere.

This reuses ``sjvc.services.junction_metadata``'s fast R-side TSV projection
and alias-tolerant column resolution, but reads the *wider* column set
(``sjlookup``'s ``junction_lookup_to_tsv.R``) so the ``annotated`` /
``left_annotated`` / ``right_annotated`` detail columns are available too.

The classification itself mirrors ``sjlookup.services.lookup.
classify_unannotated`` (same strand-aware donor/acceptor split against this
table's own annotation columns, not a live GENCODE transcript lookup) —
duplicated rather than imported, since ``sjlookup`` already depends on
``sjvc`` the other way round (see that module's own docstring); importing it
back from here would make a package cycle. One difference: SJ Lookup
classifies a single queried junction at a time and can afford to filter the
whole table per query, while a heatmap needs this for up to
``heatmap.MAX_ROWS`` junctions per request, so the "does this sit inside an
already-annotated intron of the same gene" check here is backed by a
per-gene intron index built once (lazily, on first use) instead of a fresh
table scan per row.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from .junction_metadata import (
    _END_ALIASES,
    _GENE_ID_ALIASES,
    _GENE_NAME_ALIASES,
    _MISSING_TOKENS,
    _SEQ_ALIASES,
    _START_ALIASES,
    _STRAND_ALIASES,
    _find,
    _read_junction_metadata_df,
    JunctionMetadataError,
)
from .rds import MatrixError

_LOOKUP_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "junction_lookup_to_tsv.R"
_DETAIL_COLS = ("annotated", "left_annotated", "right_annotated")


def _clean(s: pd.Series) -> pd.Series:
    """Same "NA"-string-as-missing normalisation as sjlookup's lookup index."""
    s = s.astype(str).str.strip()
    return s.mask(s.str.lower().isin(_MISSING_TOKENS))


def _site_annotated(v: object) -> bool:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return False
    s = str(v).strip()
    return s != "" and s != "0"


@dataclass
class JunctionTypeIndex:
    """Built once when a cohort's junction metadata table loads; `labels_for`
    is then a handful of indexed lookups, never a scan of the cohort."""

    df: pd.DataFrame  # indexed by "chr:start-end:strand"; see columns below
    has_annotation_detail: bool
    n_rows: int
    _introns_by_gene: Optional[Dict[str, List[Tuple[float, float]]]] = field(
        default=None, repr=False, compare=False
    )

    def _annotated_introns_by_gene(self) -> Dict[str, List[Tuple[float, float]]]:
        """gene_id -> this table's own annotated (start, end) introns for that
        gene — the same "already a known intron?" check
        ``sjvc.services.classify.GeneClassifier`` makes against live GENCODE
        transcripts (its ``novel_exon`` category), done here against this
        table instead. Computed once, lazily (only heatmaps that actually hit
        the "neither splice site individually annotated" case need it)."""
        if self._introns_by_gene is None:
            ann = self.df.loc[self.df["annotated"] == True, ["start", "end", "gene_id"]]
            ann = ann.dropna(subset=["gene_id"])
            introns: Dict[str, List[Tuple[float, float]]] = {}
            if not ann.empty:
                exploded = ann.assign(gene_id=ann["gene_id"].str.split(",")).explode("gene_id")
                exploded["gene_id"] = exploded["gene_id"].str.strip()
                exploded = exploded[exploded["gene_id"] != ""]
                for gid, g in exploded.groupby("gene_id", sort=False):
                    introns[gid] = list(zip(g["start"], g["end"]))
            self._introns_by_gene = introns
        return self._introns_by_gene

    def _classify(self, row: pd.Series) -> str:
        if row.get("annotated") is True:
            return "annotated"
        if not self.has_annotation_detail:
            return "unclassified"
        strand = row.get("strand")
        left = _site_annotated(row.get("left_annotated"))
        right = _site_annotated(row.get("right_annotated"))
        donor, acceptor = (right, left) if strand == "-" else (left, right)
        if donor and acceptor:
            return "isoform_switch"
        if donor != acceptor:
            return "alt_3p" if donor else "alt_5p"
        # neither splice site individually annotated — novel_exon if [start,
        # end] sits strictly inside one of this gene's own known introns
        s, e = row.get("start"), row.get("end")
        if pd.notna(s) and pd.notna(e):
            introns_by_gene = self._annotated_introns_by_gene()
            for gid in str(row.get("gene_id") or "").split(","):
                gid = gid.strip()
                for (a, b) in introns_by_gene.get(gid, ()):
                    if a <= s and e <= b and (a, b) != (s, e):
                        return "novel_exon"
        return "novel"

    def labels_for(self, rownames: Sequence[str]) -> List[Optional[str]]:
        """"<gene name>:<junction type>" for each rowname, in order — `None`
        where the table has no gene annotation for that junction (intergenic,
        or just absent from the table), so the caller can fall back to
        showing that row's plain junction id instead."""
        sub = self.df.reindex(rownames)
        out: List[Optional[str]] = []
        for row in sub.to_dict("records"):
            gname = row.get("gene_name")
            if gname is None or (isinstance(gname, float) and pd.isna(gname)):
                out.append(None)
                continue
            out.append(f"{gname}:{self._classify(row)}")
        return out


def load_junction_type_index(path: Path, tmp_dir: Path) -> JunctionTypeIndex:
    try:
        df = _read_junction_metadata_df(path, tmp_dir, script=_LOOKUP_SCRIPT)
    except MatrixError as e:
        raise JunctionMetadataError(str(e))
    cols = {str(c).strip().lower(): str(c) for c in df.columns}

    seq_c = _find(cols, _SEQ_ALIASES)
    start_c = _find(cols, _START_ALIASES)
    end_c = _find(cols, _END_ALIASES)
    strand_c = _find(cols, _STRAND_ALIASES)
    missing = [
        label for label, c in (
            ("seqnames", seq_c), ("start", start_c), ("end", end_c), ("strand", strand_c),
        ) if c is None
    ]
    if missing:
        raise JunctionMetadataError(
            "the junction metadata table is missing column(s): " + ", ".join(missing)
        )
    gid_c = _find(cols, _GENE_ID_ALIASES)
    gname_c = _find(cols, _GENE_NAME_ALIASES)
    if gname_c is None:
        raise JunctionMetadataError("the junction metadata table has no gene-name column")
    detail_c = {name: cols.get(name) for name in _DETAIL_COLS}

    n_rows = len(df)
    rowname = (
        df[seq_c].astype(str) + ":" + df[start_c].astype(str)
        + "-" + df[end_c].astype(str) + ":" + df[strand_c].astype(str)
    )

    out = pd.DataFrame(index=pd.Index(rowname.to_numpy(), name="junction"))
    out["start"] = pd.to_numeric(df[start_c], errors="coerce").to_numpy()
    out["end"] = pd.to_numeric(df[end_c], errors="coerce").to_numpy()
    out["strand"] = df[strand_c].astype(str).to_numpy()
    out["gene_id"] = _clean(df[gid_c]).to_numpy() if gid_c else None
    out["gene_name"] = _clean(df[gname_c]).to_numpy()
    if detail_c["annotated"]:
        out["annotated"] = _clean(df[detail_c["annotated"]]).map({"1": True, "0": False}).to_numpy()
    else:
        out["annotated"] = None
    out["left_annotated"] = _clean(df[detail_c["left_annotated"]]).to_numpy() if detail_c["left_annotated"] else None
    out["right_annotated"] = _clean(df[detail_c["right_annotated"]]).to_numpy() if detail_c["right_annotated"] else None

    out["gene_name"] = out["gene_name"].astype("category")
    out = out[~out.index.duplicated(keep="first")]

    return JunctionTypeIndex(
        df=out,
        has_annotation_detail=bool(detail_c["left_annotated"] and detail_c["right_annotated"]),
        n_rows=n_rows,
    )
