"""Cohort stratification for SJSurv.

`prepTCGAdata::get_tcga_data()` writes the *raw* sample metadata — histology,
pathologic stage, age at diagnosis, overall survival — but no longer bakes a
Good/Poor survivor label into it (that used to be `surv_cohort()`, an R
function baked into the package). Stratifying the cohort now happens here,
interactively, so the age bands can be tuned to the cohort at hand instead of
being fixed at download time:

1. ``parse_raw_metadata`` reads the raw table (rows = samples, keyed by
   `sample_id`, matched against the sjdat matrix's columns). Age and overall
   survival aren't always under the same column name across cohorts/downloads
   — the most complete candidate present is used for each (see
   ``_most_complete_col``); if no pre-computed survival-time column exists at
   all, it's built from days-to-death (deceased) / days-to-last-follow-up
   (censored) instead.
2. ``stratify`` derives, for a chosen set of age bands and a minimum group
   size, what `surv_cohort()` used to plus a proper accounting for
   right-censoring that it never did:

   * ``Histology``       — copy of the histological-diagnosis column, with
     any user-chosen merges applied (see `histology_map`), or dropped from
     `Group` entirely (`use_histology=False`) — some cohorts are
     heterogeneous enough that splitting on histology leaves cohorts too
     small to label.
   * ``Stage``            — `"Early"` (I/II) / `"Late"` (III/IV) / `None`.
   * ``Age_at_diagnosis`` — the age, binned (see `AgeBands` below).
   * ``Group``            — `"Histology | Stage | Age_at_diagnosis"` (or just
     `"Stage | Age_at_diagnosis"` with histology dropped), or `None` when any
     remaining component is `None`.
   * ``MedianSurvival``   — for every `Group` with at least `min_group_n`
     samples, the group's median survival time estimated from a *stratified
     Cox proportional-hazards model* (one call to `lifelines.CoxPHFitter`
     across every qualifying group at once, `strata=["Group"]`, no other
     covariates — with nothing to regress on this reduces to each group's own
     Breslow-estimated baseline survival curve, the Cox-model analogue of a
     Kaplan-Meier curve), not a plain median of observed follow-up times. That
     distinction matters whenever a group's longest survivors are still alive
     (censored) rather than dead: a plain median of raw follow-up days is
     biased low in that case, since a censored sample's *observed* time
     under-states its true (unknown, but at-least-this-long) survival. A
     group whose survival curve never drops to 0.5 within its observed
     follow-up (fewer than half its members have died) has no defined median
     and, like a group below `min_group_n`, gets no labels. The 0.5 split
     point itself is configurable (`stratify`'s `event_quantile`) — a
     low-mortality cohort can lower it (e.g. to 0.25, "died before 25% of the
     group had died") so a threshold is reached for more groups, at the cost
     of a more lopsided split.
   * ``SurviverGroup``    — `"Good"` when a sample's own observed follow-up
     time is at least its group's Cox-estimated `MedianSurvival` (true
     regardless of whether that follow-up ended in death or is still
     ongoing); `"Poor"` when it is shorter *and* the sample is a confirmed
     death (vital status "Dead") before that point; and `None` (unlabelled,
     not forced to "Poor") when it is shorter but the sample was last known
     alive — its eventual outcome relative to the group median is genuinely
     unknown, and mislabelling a still-living, merely-early-censored patient
     "Poor" is exactly the bias this Cox-based approach exists to remove.
     Vital status is read from the metadata file's own vital-status column
     (see `_VITAL_STATUS_ALIASES`); a sample with no recognisable vital
     status — or every sample, when the file carries no vital-status column
     at all — is conservatively treated as a confirmed death (the same
     "every observed time is a real event" assumption this module used
     unconditionally before censoring was accounted for).

`AgeBands` defaults to the classic fixed cut points `(30-50], (50-70], (>70)`
(ages of 30 or under are left unbanded — i.e. excluded from `Group`, exactly as
`surv_cohort()` did). `quantile_age_bands` offers a data-driven alternative:
equal-count bands computed from the cohort's own age distribution, which often
balances group sizes better than fixed cut points on a small or skewed cohort.
"""
from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

from sjvc.services.rds import MatrixError, read_table

# Checked in this order — a real sample_id-named column always wins. "rownames"
# is last: pyreadr names a data.frame's own (non-default) R row names that way
# when it reads one — recount3/TCGA sample-metadata tables commonly carry the
# sample UUID only as row names, with no separate sample_id column at all.
_SAMPLE_ID_ALIASES = ("sample_id", "sampleid", "sample id", "sample.id", "rownames")
# Exact prepTCGAdata::get_tcga_data() raw column names first, with only a
# couple of low-risk synonyms as fallback. Deliberately *not* aliased:
# "histology", "stage", "age_at_diagnosis" — those are exactly the *derived*
# column names the old, now-removed R-side surv_cohort() used to write into
# this same file. A cohort downloaded before that removal can still carry
# them (already binned / stringy, e.g. Age_at_diagnosis = "(50-70]"); matching
# them here would silently stratify on stale output instead of the raw data.
_HIST_ALIASES = ("tcga.cgc_case_histological_diagnosis", "histological_diagnosis")
_STAGE_ALIASES = ("tcga.cgc_case_pathologic_stage", "pathologic_stage")
# Age in *years* — which of these a cohort actually has (and how complete it
# is) varies: some downloads only get the harmonized prepTCGAdata column,
# others only the raw TCGA fields (also years, unlike GDC's days-based
# age_at_diagnosis — deliberately not aliased here to avoid silently reading
# days as years). The most complete one present is used, not just the first
# — see _most_complete_col().
_AGE_ALIASES = (
    "age_at_diagnosis_years",
    "tcga.cgc_case_age_at_diagnosis",
    "tcga.xml_age_at_initial_pathologic_diagnosis",
)
_OS_ALIASES = ("survival_days", "os_days", "survival_time")
# Fallback when a cohort has no pre-computed survival-time column at all: the
# standard TCGA construction is days-to-death for the deceased, days-to-last-
# follow-up for everyone else (censored/alive). Each half is picked by
# _most_complete_col() too, since — like age — the raw field name varies by
# cohort/download. Deliberately narrow (exact field names only, not a broad
# "days_to_*" pattern): several other TCGA "days_to_*" fields exist
# (days_to_birth, days_to_collection, days_to_initial_pathologic_diagnosis,
# …) that measure something other than time-from-diagnosis-to-event, and
# must never be swept in here.
_DAYS_TO_DEATH_ALIASES = (
    "days_to_death",
    "tcga.cgc_case_days_to_death",
    "tcga.gdc_cases.diagnoses.days_to_death",
    "tcga.cgc_follow_up_days_to_death",
    "tcga.xml_days_to_death",
)
_DAYS_TO_LAST_FOLLOWUP_ALIASES = (
    "days_to_last_follow_up",
    "tcga.cgc_case_days_to_last_follow_up",
    "tcga.gdc_cases.diagnoses.days_to_last_follow_up",
    "tcga.cgc_follow_up_days_to_last_follow_up",
    "tcga.xml_days_to_last_followup",
    "tcga.xml_days_to_last_known_alive",
)
# Whether a sample is a confirmed death (an "event", in survival-analysis
# terms) or was last known alive ("censored" — its true survival is at least
# this long, but unknown beyond that) — needed to compute a group's median
# survival properly (see `_cox_group_medians`) instead of naively treating
# every observed follow-up time as if it were a death. Several equivalent
# columns exist across TCGA download vintages; the most complete one present
# is used, same pattern as age/survival-time above.
_VITAL_STATUS_ALIASES = (
    "vital_status",
    "tcga.gdc_cases.diagnoses.vital_status",
    "tcga.cgc_case_vital_status",
    "tcga.cgc_follow_up_vital_status",
    "tcga.xml_vital_status",
)
_ALIVE_VALUES = {"alive", "living"}
_LEGACY_DERIVED_COLS = {"histology", "stage", "age_at_diagnosis", "group", "survivergroup", "mediansurvival"}
_MISSING = {"", "nan", "none", "na", "n/a", "null"}
_ABC_SUFFIX = re.compile(r"[abcd][0-9]?$")
# A column whose name matches this is offered as a covariate *and* checked on
# by default (with Stage / Age At Diagnosis / Histology) — same pattern
# dataload/scan.py uses to spot an `*_MSI.rds` file, applied here to a column
# name instead of a filename.
_MSI_COL = re.compile(r"(^|[_.\-])(msi(sensor)?|mantis)([_.\-]|$)", re.I)


def _looks_like_identifier_col(col: str) -> bool:
    """A column whose name is itself an id/uuid/barcode segment (`sample_id`,
    `tcga.gdc_file_id`, a joined table's own `msi.sample_id`, ...) — never a
    useful covariate, and if offered, a categorical one one-hots into nearly
    as many columns as there are samples."""
    segments = re.split(r"[_.\-]+", col.lower())
    return any(seg in ("id", "ids", "uuid", "barcode") for seg in segments)


class MetadataError(ValueError):
    pass


class StratifyError(ValueError):
    pass


def _is_missing(v: object) -> bool:
    return v is None or str(v).strip().lower() in _MISSING


def _to_float(v: object) -> Optional[float]:
    if _is_missing(v):
        return None
    try:
        return float(str(v).strip())
    except ValueError:
        return None


# --------------------------------------------------------------------------- #
# covariates — sample aspects (from the loaded file) a classifier can use
# alongside the selected molecular features. Histology/Stage/Age At Diagnosis
# are keyed specially since they're already parsed into their own `rows`
# fields under cryptic real column names; everything else is keyed by its
# raw column name directly.
# --------------------------------------------------------------------------- #
COV_HISTOLOGY = "__histology__"
COV_STAGE = "__stage__"
COV_AGE = "__age__"


@dataclass
class CovariateColumn:
    key: str                # COV_HISTOLOGY/COV_STAGE/COV_AGE, or a raw column name
    label: str               # display label
    kind: str                 # "numeric" | "categorical"
    n_available: int          # samples with a non-missing value
    default: bool             # checked by default (Stage/Age/Histology, or an MSI-named column)


# --------------------------------------------------------------------------- #
# raw sample metadata
# --------------------------------------------------------------------------- #
@dataclass
class RawMetadata:
    # sample_id -> {"histology": str|None, "stage_raw": str|None,
    #               "age": float|None, "os": float|None, "event": bool}
    # "event" is True for a confirmed death, False for a sample last known
    # alive (censored) — see _event_value(); always True when the file has no
    # recognisable vital-status column at all.
    rows: Dict[str, Dict[str, object]]
    n_unmatched: int = 0
    # every other column in the sample-metadata file (i.e. not the id column,
    # not histology/stage/age/survival, and never a survival-time-derived one
    # — see the exclusion list built in parse_raw_metadata) — sample_id ->
    # {raw column name -> raw value}. The candidate pool `available_covariates()`
    # offers alongside the three specials below.
    covariate_table: Dict[str, Dict[str, object]] = field(default_factory=dict)
    covariate_columns: List[str] = field(default_factory=list)  # in file order
    # whether a recognisable vital-status column was found at all — when
    # False, every sample's "event" above is the True-by-default fallback,
    # so the Cox-based MedianSurvival below is not actually censoring-
    # adjusted for this file (surfaced as a warning by the caller).
    has_vital_status: bool = False

    @property
    def n_rows(self) -> int:
        return len(self.rows)

    def ages(self) -> List[Optional[float]]:
        return [r["age"] for r in self.rows.values()]

    def age_summary(self) -> Dict[str, object]:
        finite = [a for a in self.ages() if a is not None and math.isfinite(a)]
        return {
            "n_with_age": len(finite),
            "age_min": float(min(finite)) if finite else None,
            "age_max": float(max(finite)) if finite else None,
        }

    def histology_counts(self) -> List[Tuple[str, int]]:
        """Distinct raw histology values and how many samples carry each,
        most-common first — what the merge-editor shows the user so they can
        decide which values to fold together (or drop histology entirely; see
        stratify()'s `use_histology`)."""
        c = Counter(r["histology"] for r in self.rows.values() if r["histology"] is not None)
        return sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))

    def available_covariates(self) -> List["CovariateColumn"]:
        """Every sample aspect in the loaded file a classifier could use as a
        covariate alongside the selected molecular features: the three
        columns already read for `Group` (under friendly names, since their
        real column names are the cryptic prepTCGAdata ones), plus every
        other column in the file — excluding the id column and anything
        survival-time-derived (`build_covariates()`'s caller never even sees
        those; including the very quantity `SurviverGroup` is split on would
        let the classifier "predict" it for free instead of testing whether
        the molecular data carries the signal)."""
        out = [
            CovariateColumn(
                key=COV_HISTOLOGY, label="Histology", kind="categorical",
                n_available=sum(1 for r in self.rows.values() if r["histology"] is not None),
                default=True,
            ),
            CovariateColumn(
                key=COV_STAGE, label="Stage", kind="categorical",
                n_available=sum(1 for r in self.rows.values() if stage_label(r["stage_raw"]) is not None),
                default=True,
            ),
            CovariateColumn(
                key=COV_AGE, label="Age At Diagnosis", kind="numeric",
                n_available=sum(1 for r in self.rows.values() if r["age"] is not None),
                default=True,
            ),
        ]
        for col in self.covariate_columns:
            if _looks_like_identifier_col(col):
                continue
            present = [row.get(col) for row in self.covariate_table.values() if not _is_missing(row.get(col))]
            if not present:
                continue
            n_numeric = sum(1 for v in present if _to_float(v) is not None)
            kind = "numeric" if n_numeric >= max(1, round(0.9 * len(present))) else "categorical"
            if kind == "categorical":
                n_distinct = len({str(v) for v in present})
                if n_distinct > max(20, round(0.5 * len(present))):
                    continue  # near-unique per sample — free-text/id-like, not a usable covariate
            out.append(CovariateColumn(
                key=col, label=_prettify_col(col), kind=kind, n_available=len(present),
                default=bool(_MSI_COL.search(col)),
            ))
        return out


def _prettify_col(name: str) -> str:
    """A raw column name (often a dotted prepTCGAdata field like
    `tcga.cgc_case_icd_10`) as a short display label — the part after the
    last `.`, underscores/hyphens turned to spaces, capitalized."""
    base = name.rsplit(".", 1)[-1]
    base = re.sub(r"[_\-]+", " ", base).strip()
    return (base[:1].upper() + base[1:]) if base else name


def _find_col(cols: Dict[str, str], aliases: Sequence[str]) -> Optional[str]:
    """First alias (checked in order — the canonical name goes first) that is
    both present and not one of the old surv_cohort()-derived column names."""
    return next(
        (cols[a] for a in aliases if a in cols and a not in _LEGACY_DERIVED_COLS), None
    )


def _most_complete_col(
    df: pd.DataFrame, cols: Dict[str, str], aliases: Sequence[str],
    *, usable=lambda v: _to_float(v) is not None,
) -> Optional[str]:
    """Among the given aliases that are present, the one with the most usable
    values — unlike _find_col, not just the first one present. Different TCGA
    cohort downloads expose different subsets of the equivalent fields, and a
    present column isn't necessarily a populated one. `usable` decides what
    counts as a populated cell — the default (parseable as a float) suits the
    numeric age/survival-time columns; pass `lambda v: not _is_missing(v)` for
    a categorical column such as vital status."""
    present = [cols[a] for a in aliases if a in cols and a not in _LEGACY_DERIVED_COLS]
    if not present:
        return None
    if len(present) == 1:
        return present[0]
    return max(present, key=lambda c: sum(1 for v in df[c] if usable(v)))


def _event_value(v: object) -> bool:
    """Whether a vital-status cell means a confirmed death (an "event") for
    survival-analysis purposes. Anything not recognisably "alive" — missing,
    "Dead", "Deceased", or an unexpected value — is conservatively treated as
    a death: the whole point of tracking vital status is to stop treating a
    still-living, early-censored sample as if it were a death, and that only
    works by erring toward "event" whenever status is otherwise unknown,
    exactly the assumption this module made unconditionally before censoring
    was accounted for."""
    return not (not _is_missing(v) and str(v).strip().lower() in _ALIVE_VALUES)


def parse_raw_metadata(path: Path, tmp_dir: Path, known_samples: Sequence[str]) -> RawMetadata:
    try:
        df = read_table(path, tmp_dir)
    except MatrixError as e:
        raise MetadataError(str(e))
    cols = {str(c).strip().lower(): str(c) for c in df.columns}

    id_col = _find_col(cols, _SAMPLE_ID_ALIASES)
    if id_col is None:
        raise MetadataError(
            "no `sample_id` column found in the sample metadata; columns: "
            + ", ".join(map(str, df.columns))
        )
    hist_col = _find_col(cols, _HIST_ALIASES)
    stage_col = _find_col(cols, _STAGE_ALIASES)
    age_col = _most_complete_col(df, cols, _AGE_ALIASES)
    os_col = _find_col(cols, _OS_ALIASES)
    death_col = fu_col = None
    if os_col is None:
        # no pre-computed survival-time column — build one from the raw
        # days-to-death / days-to-last-follow-up fields instead (see the
        # alias comment above for why these two are picked separately)
        death_col = _most_complete_col(df, cols, _DAYS_TO_DEATH_ALIASES)
        fu_col = _most_complete_col(df, cols, _DAYS_TO_LAST_FOLLOWUP_ALIASES)
    vital_col = _most_complete_col(
        df, cols, _VITAL_STATUS_ALIASES, usable=lambda v: not _is_missing(v),
    )

    missing = [
        label for label, c in (
            ("histological diagnosis (tcga.cgc_case_histological_diagnosis)", hist_col),
            ("pathologic stage (tcga.cgc_case_pathologic_stage)", stage_col),
            ("age at diagnosis (" + " / ".join(_AGE_ALIASES) + ")", age_col),
        ) if c is None
    ]
    if os_col is None and death_col is None and fu_col is None:
        missing.append(
            "overall survival (survival_days, or days_to_death / days_to_last_follow_up)"
        )
    if missing:
        raise MetadataError(
            "the sample metadata is missing column(s) needed to stratify the cohort: "
            + "; ".join(missing)
            + ". This should be the TCGA_<cohort>_sample_metadata.rds file written by "
              "prepTCGAdata::get_tcga_data()."
        )

    def _os_value(r: pd.Series) -> Optional[float]:
        if os_col is not None:
            return _to_float(r[os_col])
        # standard TCGA construction: time-to-death for the deceased,
        # time-to-last-follow-up (censored) for everyone else
        v = _to_float(r[death_col]) if death_col is not None else None
        if v is None and fu_col is not None:
            v = _to_float(r[fu_col])
        return v

    # every other column is a covariate candidate — never the id column, never
    # histology/stage/age (already exposed as the Histology/Stage/Age At
    # Diagnosis specials), and never a survival-time or vital-status column
    # (os_col / death_col / fu_col / vital_col): those are exactly the
    # quantities SurviverGroup is computed from, and offering them as a
    # covariate would let the classifier "predict" the label for free instead
    # of testing the molecular data.
    _excluded_cov_cols = {
        c for c in (id_col, hist_col, stage_col, age_col, os_col, death_col, fu_col, vital_col) if c
    }
    covariate_columns = [
        str(c) for c in df.columns
        if str(c) not in _excluded_cov_cols and str(c).strip().lower() not in _LEGACY_DERIVED_COLS
    ]

    known = {str(s) for s in known_samples}
    filter_to_known = bool(known)   # no sjdat loaded yet -> keep every row
    rows: Dict[str, Dict[str, object]] = {}
    covariate_table: Dict[str, Dict[str, object]] = {}
    n_unmatched = 0
    for _, r in df.iterrows():
        sid = str(r[id_col]).strip()
        if _is_missing(sid):
            continue
        if filter_to_known and sid not in known:
            n_unmatched += 1
            continue
        rows[sid] = {
            "histology": None if _is_missing(r[hist_col]) else str(r[hist_col]).strip(),
            "stage_raw": None if _is_missing(r[stage_col]) else str(r[stage_col]).strip(),
            "age": _to_float(r[age_col]),
            "os": _os_value(r),
            "event": _event_value(r[vital_col]) if vital_col is not None else True,
        }
        covariate_table[sid] = {c: r[c] for c in covariate_columns}

    if not rows:
        raise MetadataError(
            "no `sample_id` in the metadata matched a column of the sjdat matrix"
        )
    return RawMetadata(
        rows=rows, n_unmatched=n_unmatched, has_vital_status=vital_col is not None,
        covariate_table=covariate_table, covariate_columns=covariate_columns,
    )


# --------------------------------------------------------------------------- #
# Stage: I/II -> Early, III/IV -> Late (same rule as the old surv_cohort())
# --------------------------------------------------------------------------- #
def stage_label(raw: Optional[str]) -> Optional[str]:
    if raw is None:
        return None
    st = str(raw).strip().lower()
    if not st:
        return None
    st = re.sub(r"^stage\s+", "", st)          # "stage iia" -> "iia"
    st = _ABC_SUFFIX.sub("", st)                # "iia" -> "ii", "iv" unchanged
    if st in ("i", "ii"):
        return "Early"
    if st in ("iii", "iv"):
        return "Late"
    return None


# --------------------------------------------------------------------------- #
# Age_at_diagnosis bands
# --------------------------------------------------------------------------- #
def _fmt_edge(x: float) -> str:
    if math.isinf(x):
        return "∞"
    return f"{x:g}"


def _make_labels(edges: List[float], include_lowest: bool) -> List[str]:
    labels = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if math.isinf(hi):
            labels.append(f"(>{_fmt_edge(lo)})")
        elif i == 0 and include_lowest:
            labels.append(f"[{_fmt_edge(lo)}-{_fmt_edge(hi)}]")
        else:
            labels.append(f"({_fmt_edge(lo)}-{_fmt_edge(hi)}]")
    return labels


@dataclass
class AgeBands:
    edges: List[float]                    # ascending, length n_bands + 1; last may be inf
    include_lowest: bool = False          # left-inclusive first band (True for quantile bands)
    labels: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if len(self.edges) < 2:
            raise StratifyError("age bands need at least 2 edges (1 band)")
        if any(b <= a for a, b in zip(self.edges, self.edges[1:])):
            raise StratifyError("age-band edges must be strictly increasing")
        if not self.labels:
            self.labels = _make_labels(self.edges, self.include_lowest)

    @property
    def n_bands(self) -> int:
        return len(self.edges) - 1


# The default, exactly as `surv_cohort()` had it: (30-50], (50-70], (>70) —
# ages of 30 or under fall below the first edge and are left unbanded (None).
CLASSIC_AGE_BANDS = AgeBands(edges=[30.0, 50.0, 70.0, math.inf], include_lowest=False)
DEFAULT_MIN_GROUP_N = 20


def quantile_age_bands(ages: Sequence[Optional[float]], n_bands: int) -> AgeBands:
    """A data-driven alternative to the classic fixed cut points: `n_bands`
    equal-count bands spanning this cohort's own age distribution (nothing is
    left unbanded). Useful when the fixed 30/50/70 cut points fragment a small
    or unusually-distributed cohort into groups too small to label."""
    if n_bands < 1:
        raise StratifyError("the number of age bands must be >= 1")
    finite = sorted(a for a in ages if a is not None and math.isfinite(a))
    if len(finite) < max(n_bands, 4):
        raise StratifyError(
            f"only {len(finite)} sample(s) have a usable age — too few to compute "
            f"{n_bands} band(s)"
        )
    qs = np.quantile(finite, np.linspace(0, 1, n_bands + 1))
    # Round the outer edges outward (floor the min, ceil the max) rather than to
    # the nearest 0.1 — rounding the true max *down* would strand the cohort's
    # oldest sample(s) above the last edge and leave them unbanded, defeating
    # the point of a "covers everyone" alternative to the fixed cut points.
    edges = [math.floor(float(qs[0]) * 10) / 10]
    last = len(qs) - 1
    for i, q in enumerate(qs[1:], start=1):
        v = math.ceil(float(q) * 10) / 10 if i == last else round(float(q), 1)
        if v <= edges[-1]:            # ties (many equal ages) -> keep strictly increasing
            v = round(edges[-1] + 0.1, 1)
        edges.append(v)
    return AgeBands(edges=edges, include_lowest=True)


def bin_ages(ages: Sequence[Optional[float]], bands: AgeBands) -> List[Optional[str]]:
    s = pd.Series([np.nan if a is None else a for a in ages], dtype=float)
    cats = pd.cut(s, bins=bands.edges, right=True, include_lowest=bands.include_lowest,
                   labels=bands.labels)
    return [None if pd.isna(v) else str(v) for v in cats]


# --------------------------------------------------------------------------- #
# stratified output (what the rest of SJSurv consumes)
# --------------------------------------------------------------------------- #
ALL_GROUPS = "__all__"
ALL_GROUPS_LABEL = "All labelled samples (any Group)"
NO_GROUP = "(no Group)"


@dataclass
class GroupCount:
    group: str
    n_total: int
    n_good: int
    n_poor: int
    n_labelled: int


@dataclass
class Metadata:
    # sample_id -> {"group": str|None, "surviver": "Good"|"Poor"|None}
    rows: Dict[str, Dict[str, object]]
    n_unmatched: int = 0

    @property
    def n_rows(self) -> int:
        return len(self.rows)

    def group_counts(self) -> List[GroupCount]:
        buckets: Dict[str, List[int]] = {}
        for r in self.rows.values():
            g = r["group"] or NO_GROUP
            b = buckets.setdefault(g, [0, 0, 0])
            b[0] += 1
            if r["surviver"] == "Good":
                b[1] += 1
            elif r["surviver"] == "Poor":
                b[2] += 1

        no_group = buckets.pop(NO_GROUP, None)
        ordered = sorted(buckets.items(), key=lambda kv: (-sum(kv[1][1:]), kv[0]))
        out = [
            GroupCount(group=g, n_total=b[0], n_good=b[1], n_poor=b[2], n_labelled=b[1] + b[2])
            for g, b in ordered
        ]
        if no_group is not None:
            out.append(GroupCount(
                group=NO_GROUP, n_total=no_group[0], n_good=no_group[1],
                n_poor=no_group[2], n_labelled=no_group[1] + no_group[2],
            ))

        good = sum(gc.n_good for gc in out)
        poor = sum(gc.n_poor for gc in out)
        out.insert(0, GroupCount(
            group=ALL_GROUPS, n_total=sum(gc.n_total for gc in out),
            n_good=good, n_poor=poor, n_labelled=good + poor,
        ))
        return out

    def labelled_samples(self, group: str) -> Dict[str, str]:
        """``sample_id -> "Good"|"Poor"`` for the chosen group (or every
        labelled sample when ``group == ALL_GROUPS``)."""
        out: Dict[str, str] = {}
        for sid, r in self.rows.items():
            if r["surviver"] not in ("Good", "Poor"):
                continue
            if group != ALL_GROUPS and (r["group"] or NO_GROUP) != group:
                continue
            out[sid] = r["surviver"]
        return out


def stratify(
    raw: RawMetadata,
    age_bands: AgeBands,
    min_group_n: int,
    *,
    use_histology: bool = True,
    histology_map: Optional[Dict[str, str]] = None,
    event_quantile: float = 0.5,
) -> Metadata:
    """`use_histology=False` drops histology from Group entirely (Group
    becomes just Stage | age band, and a sample no longer needs a histology
    value to get one) — some cancers are heterogeneous enough that splitting
    on histology fragments the cohort into cohorts too small to label.
    `histology_map` (raw value -> merged label) folds selected histology
    values together before they go into Group, e.g. several closely-related
    subtypes that are individually too small; values not present in the map
    are used as-is.

    `event_quantile` is the fraction of a group that must have died for its
    Cox-estimated survival-time threshold (`MedianSurvival`) to be defined —
    0.5 is the classic median. A low-mortality cohort's 50% threshold can go
    unreached for most groups within their observed follow-up (a `None`
    result from `_cox_group_survival_quantile_times`, same as a group too
    small to fit at all); lowering this toward, say, 0.25 reaches a defined
    threshold sooner (fewer deaths needed to define it), at the cost of a
    more lopsided Good/Poor split — "Poor" then means "died unusually early
    relative to the rest of the group," not "died before the halfway point."
    """
    if min_group_n < 1:
        raise StratifyError("min_group_n must be >= 1")
    if not (0.0 < event_quantile < 1.0):
        raise StratifyError("event_quantile must be strictly between 0 and 1")

    sample_ids = list(raw.rows)
    age_band = bin_ages([raw.rows[s]["age"] for s in sample_ids], age_bands)
    hmap = histology_map or {}

    group: List[Optional[str]] = []
    for i, s in enumerate(sample_ids):
        st = stage_label(raw.rows[s]["stage_raw"])
        ab = age_band[i]
        if not use_histology:
            group.append(None if (st is None or ab is None) else f"{st} | {ab}")
            continue
        h = raw.rows[s]["histology"]
        h = hmap.get(h, h) if h is not None else h
        group.append(None if (h is None or st is None or ab is None) else f"{h} | {st} | {ab}")

    os_vals = [raw.rows[s]["os"] for s in sample_ids]
    event_vals = [raw.rows[s]["event"] for s in sample_ids]
    by_group: Dict[str, List[int]] = {}
    for i, g in enumerate(group):
        if g is not None:
            by_group.setdefault(g, []).append(i)

    # Cox-based survival-quantile threshold per qualifying group (>=
    # min_group_n members, regardless of how many actually have a usable os
    # value — same gate the old plain-median code used) — one stratified fit
    # across every such group at once, see _cox_group_survival_quantile_times().
    fit_rows: List[Tuple[float, bool, str]] = [
        (os_vals[i], event_vals[i], g)
        for g, idxs in by_group.items() if len(idxs) >= min_group_n
        for i in idxs if os_vals[i] is not None
    ]
    medians_by_group = _cox_group_survival_quantile_times(fit_rows, event_quantile)

    median_survival: List[Optional[float]] = [
        medians_by_group.get(g) if g is not None else None for g in group
    ]

    rows: Dict[str, Dict[str, object]] = {}
    for i, s in enumerate(sample_ids):
        surviver = None
        if median_survival[i] is not None and os_vals[i] is not None:
            if os_vals[i] >= median_survival[i]:
                surviver = "Good"
            elif event_vals[i]:
                surviver = "Poor"
            # else: shorter than the group median, but last known alive
            # (censored) — genuinely unknown whether they'd have reached it,
            # so left unlabelled rather than forced to "Poor"
        rows[s] = {"group": group[i], "surviver": surviver}

    return Metadata(rows=rows, n_unmatched=raw.n_unmatched)


def _cox_group_survival_quantile_times(
    rows: List[Tuple[float, bool, str]], event_quantile: float = 0.5,
) -> Dict[str, Optional[float]]:
    """Per-group survival-time threshold, adjusted for right-censoring, from
    one stratified Cox proportional-hazards fit (`lifelines.CoxPHFitter`,
    `strata=["group"]`, no other covariates) across every `(duration, event,
    group)` triple given. With no covariates to regress on, this reduces to
    each group's own Breslow-estimated baseline survival curve — the
    Cox-model analogue of a Kaplan-Meier curve — read off at the first time
    point where survival drops to `1 - event_quantile` or below (the classic
    median is `event_quantile=0.5`, i.e. survival <= 0.5). A group whose
    curve never reaches that level within its observed follow-up (fewer than
    `event_quantile` of its members have died) has no defined threshold and
    is simply absent from the result, same as if it had been dropped from
    `rows` entirely."""
    if not rows:
        return {}
    survival_threshold = 1.0 - event_quantile
    df = pd.DataFrame(rows, columns=["duration", "event", "group"])
    try:
        cph = CoxPHFitter()
        cph.fit(df, duration_col="duration", event_col="event", strata=["group"])
        baseline_survival = cph.baseline_survival_
    except Exception:
        # a genuinely degenerate input (e.g. every duration identical and
        # non-positive) — no thresholds rather than a hard failure; the
        # caller already only invokes this on data that passed min_group_n,
        # so this is a rare defensive fallback, not the expected path
        return {}
    out: Dict[str, Optional[float]] = {}
    for g in df["group"].unique():
        if g not in baseline_survival.columns:
            continue
        below = baseline_survival[g][baseline_survival[g] <= survival_threshold]
        out[g] = float(below.index[0]) if not below.empty else None
    return out


# --------------------------------------------------------------------------- #
# covariate matrix — turns the chosen covariate keys into classifier-ready
# rows, aligned to `sample_ids` (the order `select_features()`/`select_from_
# features()` fixed for the molecular feature matrix, so the two stack
# directly). A numeric column becomes one row (missing imputed to the mean
# over `sample_ids`); a categorical one becomes one 0/1 row per distinct
# value *seen anywhere in the cohort* (not just `sample_ids`) so a category
# absent from one cross-validation fold still gets a stable, zero-filled
# column there rather than shifting the feature layout fold to fold.
# --------------------------------------------------------------------------- #
def build_covariates(
    raw: RawMetadata,
    keys: Sequence[str],
    sample_ids: Sequence[str],
    *,
    histology_map: Optional[Dict[str, str]] = None,
) -> Tuple[List[str], np.ndarray]:
    hmap = histology_map or {}
    names: List[str] = []
    feature_rows: List[np.ndarray] = []

    def add_numeric(label: str, values: Dict[str, Optional[float]]) -> None:
        vals = [values.get(s) for s in sample_ids]
        present = [v for v in vals if v is not None]
        fill = float(np.mean(present)) if present else 0.0
        names.append(label)
        feature_rows.append(np.array([v if v is not None else fill for v in vals], dtype=float))

    def add_categorical(label: str, values: Dict[str, Optional[str]]) -> None:
        distinct = sorted({v for v in values.values() if v is not None})
        for v in distinct:
            names.append(f"{label}={v}")
            feature_rows.append(
                np.array([1.0 if values.get(s) == v else 0.0 for s in sample_ids], dtype=float)
            )

    for key in keys:
        if key == COV_HISTOLOGY:
            values = {
                s: (hmap.get(r["histology"], r["histology"]) if r["histology"] is not None else None)
                for s, r in raw.rows.items()
            }
            add_categorical("Histology", values)
        elif key == COV_STAGE:
            add_categorical("Stage", {s: stage_label(r["stage_raw"]) for s, r in raw.rows.items()})
        elif key == COV_AGE:
            add_numeric("Age At Diagnosis", {s: r["age"] for s, r in raw.rows.items()})
        else:
            label = _prettify_col(key)
            raw_values = {s: row.get(key) for s, row in raw.covariate_table.items()}
            present = [v for v in raw_values.values() if not _is_missing(v)]
            n_numeric = sum(1 for v in present if _to_float(v) is not None)
            if present and n_numeric >= max(1, round(0.9 * len(present))):
                add_numeric(label, {s: _to_float(v) for s, v in raw_values.items()})
            else:
                add_categorical(
                    label, {s: (None if _is_missing(v) else str(v).strip()) for s, v in raw_values.items()}
                )

    if not feature_rows:
        return [], np.zeros((0, len(sample_ids)))
    return names, np.vstack(feature_rows)
