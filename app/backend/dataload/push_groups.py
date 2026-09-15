"""Push SJSurv's computed ``Group`` / ``SurviverGroup`` labels into the
Sashimi-plot and 2D View sessions' own sample-metadata tables, so they become
pickable there too — a stratification variable in Sashimi, a clinical feature
(PCA/UMAP/heatmap colouring) in 2D View.

SJSurv recomputes ``Group``/``SurviverGroup`` itself from the raw metadata,
tunable via its age-band / min-group-size controls (``sjsurv.services.
metadata.stratify``) — they are *not* read from whatever ``Group``/
``SurviverGroup`` columns the metadata file happens to already carry (a cohort
downloaded before ``surv_cohort()`` was removed from the R package can still
have stale ones baked in). Pushing overwrites same-named columns in the target
sessions with SJSurv's current, freshly-computed values, on the same
reasoning: the file's own baked columns are obsolete once SJSurv has
recomputed them.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

from sjvc.services.clinical import _CATEGORICAL_MAX_FRACTION, _is_missing, Clinical, Column

GroupMap = Dict[str, Tuple[Optional[str], Optional[str]]]  # sample_id -> (Group, SurviverGroup)

_PUSHED_COLUMNS = ("Group", "SurviverGroup")


def push_into_sample_metadata(meta, group_map: GroupMap) -> int:
    """Sashimi plot's ``sjv.services.groups.SampleMetadata`` — a plain
    ``{column: value}}`` dict per row, no per-column typing to maintain.
    Returns how many of its own samples got a non-empty Group or
    SurviverGroup value."""
    n = 0
    for sid, row in meta.rows.items():
        g, sv = group_map.get(sid, (None, None))
        row["Group"] = g or ""
        row["SurviverGroup"] = sv or ""
        if g or sv:
            n += 1
    for col in _PUSHED_COLUMNS:
        if col not in meta.columns:
            meta.columns.append(col)
    return n


def push_into_clinical(clinical: Clinical, group_map: GroupMap) -> int:
    """2D View's ``sjvc.services.clinical.Clinical`` — rows plus a typed
    ``Column`` per field, so the pushed columns need their type/eligibility
    recomputed the same way ``parse_clinical()`` does for any other column.
    Returns how many of its own samples got a non-empty Group or
    SurviverGroup value."""
    n = 0
    for sid, row in clinical.rows.items():
        g, sv = group_map.get(sid, (None, None))
        row["Group"] = g or ""
        row["SurviverGroup"] = sv or ""
        if g or sv:
            n += 1

    n_rows = clinical.n_rows
    for name in _PUSHED_COLUMNS:
        values = [clinical.rows[s].get(name, "") for s in clinical.rows]
        present = [v for v in values if not _is_missing(v)]
        n_unique = len(set(present))
        n_missing = n_rows - len(present)
        eligible = n_rows > 0 and 1 <= n_unique < _CATEGORICAL_MAX_FRACTION * n_rows
        col = Column(name=name, type="categorical", n_unique=n_unique, n_missing=n_missing, eligible=eligible)
        clinical.columns = [c for c in clinical.columns if c.name != name] + [col]
        clinical.overrides.pop(name, None)   # a stale numeric/categorical override no longer applies
    return n
