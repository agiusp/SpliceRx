"""Turn 0-2 selected clinical features into per-sample colour + shape + legend
for the 2D projection, following the agreed encoding table."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .. import palette
from .clinical import Clinical

NEUTRAL = "#8aa0b4"


def _robust01(x: np.ndarray) -> np.ndarray:
    v = x.astype(float)
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        return np.zeros_like(v)
    lo, hi = np.percentile(finite, [2, 98])
    if hi <= lo:
        lo, hi = float(finite.min()), float(finite.max())
    if hi <= lo:
        return np.zeros_like(v)
    return np.clip((v - lo) / (hi - lo), 0, 1)


def build(
    clinical: Optional[Clinical],
    feature_names: Sequence[str],
    sample_order: Sequence[str],
) -> Dict:
    n = len(sample_order)
    feats = list(feature_names)[:2]

    if clinical is None or not feats:
        return {
            "colors": [NEUTRAL] * n,
            "shapes": ["circle"] * n,
            "missing": [False] * n,
            "encoding": {"color": {"kind": "none"}, "shape": {"kind": "none"}},
            "legend": {},
        }

    typed: List[Tuple[str, str, object]] = []
    for f in feats:
        vals, kind = clinical.values_for(f, sample_order)
        typed.append((f, kind, vals))

    colors = [NEUTRAL] * n
    shapes = ["circle"] * n
    # True wherever `colors[i]` below is the "#cfcfcf" grey-for-missing-data
    # marker — i.e. exactly the points a viewer would call "grey" and the
    # projection's "hide samples missing this data" toggle removes. A missing
    # *shape*-only feature (both features picked, only the categorical one
    # missing) doesn't turn a point grey today, so it isn't "missing" here
    # either — this mirrors that existing visual, not a stricter definition.
    missing = [False] * n
    encoding: Dict = {"color": {"kind": "none"}, "shape": {"kind": "none"}}
    legend: Dict = {}

    def color_categorical(name: str, vals: List[Optional[str]]):
        cmap = palette.categorical_colors([v for v in vals if v is not None])
        for i, v in enumerate(vals):
            if v is None:
                colors[i] = "#cfcfcf"
                missing[i] = True
            else:
                colors[i] = cmap.get(v, palette.OTHER[0])
        encoding["color"] = {"kind": "categorical", "feature": name}
        legend["color"] = {
            "feature": name,
            "kind": "categorical",
            "items": [{"value": v, "color": c} for v, c in cmap.items()],
        }

    def color_numeric(name: str, vals: np.ndarray):
        t = _robust01(vals)
        for i in range(n):
            if not np.isfinite(vals[i]):
                colors[i] = "#cfcfcf"
                missing[i] = True
            else:
                colors[i] = palette.sequential_color(t[i])
        finite = vals[np.isfinite(vals)]
        encoding["color"] = {"kind": "sequential", "feature": name}
        legend["color"] = {
            "feature": name,
            "kind": "sequential",
            "min": float(finite.min()) if finite.size else 0.0,
            "max": float(finite.max()) if finite.size else 0.0,
            "stops": palette.SEQUENTIAL,
        }

    def shape_categorical(name: str, vals: List[Optional[str]]):
        smap = palette.shape_for([v for v in vals if v is not None])
        for i, v in enumerate(vals):
            shapes[i] = smap.get(v, "dot-small") if v is not None else "circle"
        encoding["shape"] = {"kind": "categorical", "feature": name}
        legend["shape"] = {
            "feature": name,
            "items": [{"value": v, "shape": s} for v, s in smap.items()],
        }

    if len(typed) == 1:
        name, kind, vals = typed[0]
        color_numeric(name, vals) if kind == "numeric" else color_categorical(name, vals)
        return {"colors": colors, "shapes": shapes, "missing": missing, "encoding": encoding, "legend": legend}

    (n1, k1, v1), (n2, k2, v2) = typed
    if k1 == "categorical" and k2 == "categorical":
        color_categorical(n1, v1)
        shape_categorical(n2, v2)
    elif k1 == "numeric" and k2 == "numeric":
        t1, t2 = _robust01(v1), _robust01(v2)
        for i in range(n):
            if not (np.isfinite(v1[i]) and np.isfinite(v2[i])):
                colors[i] = "#cfcfcf"
                missing[i] = True
            else:
                colors[i] = palette.bivariate_color(float(t1[i]), float(t2[i]))
        encoding["color"] = {"kind": "bivariate", "features": [n1, n2]}
        legend["bivariate"] = {
            "features": [n1, n2],
            "corners": {
                "lo_lo": palette.BIV_LO_LO,
                "hi_lo": palette.BIV_HI_LO,
                "lo_hi": palette.BIV_LO_HI,
                "hi_hi": palette.BIV_HI_HI,
            },
            "ranges": [
                [float(np.nanmin(v1)), float(np.nanmax(v1))],
                [float(np.nanmin(v2)), float(np.nanmax(v2))],
            ],
        }
    else:  # one numeric, one categorical
        (num_n, _, num_v) = (n1, k1, v1) if k1 == "numeric" else (n2, k2, v2)
        (cat_n, _, cat_v) = (n2, k2, v2) if k1 == "numeric" else (n1, k1, v1)
        color_numeric(num_n, num_v)
        shape_categorical(cat_n, cat_v)

    return {"colors": colors, "shapes": shapes, "missing": missing, "encoding": encoding, "legend": legend}
