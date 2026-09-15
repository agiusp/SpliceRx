"""Colour + shape system for SJVC.

Hues are the dataviz skill's validated colourblind-safe categorical palette
(references/palette.md). The frontend mirror is frontend/src/plot/palette.ts —
keep them in sync.
"""
from __future__ import annotations

import re
from typing import Dict, List, Sequence, Tuple

# ordered categorical hues (light, dark)
CATEGORICAL: List[Tuple[str, str]] = [
    ("#2a78d6", "#3987e5"),  # blue
    ("#eb6834", "#d95926"),  # orange
    ("#1baf7a", "#199e70"),  # aqua
    ("#eda100", "#c98500"),  # yellow
    ("#e87ba4", "#d55181"),  # magenta
    ("#008300", "#008300"),  # green
    ("#4a3aa7", "#9085e9"),  # violet
    ("#e34948", "#e66767"),  # red
]
OTHER = ("#898781", "#898781")

# single-hue sequential ramp (blue), light->dark stops
SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

# a small set of distinct single-hue ramps, cycled so two numeric annotation
# tracks in the same heatmap don't both come out blue
SEQUENTIAL_RAMPS: List[List[str]] = [
    SEQUENTIAL,
    ["#fbe4d3", "#f6c09a", "#ef9760", "#e87036", "#c9551f", "#9c3f16", "#6f2c0f"],  # orange
    ["#d0efe1", "#a3ddc3", "#6fc7a1", "#3aab7d", "#22855f", "#166344", "#0d4530"],  # green
    ["#f6dce8", "#eab4cd", "#dd88ae", "#cf5b8e", "#a94270", "#7f3153", "#582138"],  # magenta
]

# ordered point-shape names (frontend renders the glyph)
SHAPES = ["circle", "square", "triangle-up", "diamond", "triangle-down", "plus", "cross", "star"]

# bivariate blend corners  (feature1 low->high on one axis, feature2 on the other)
BIV_LO_LO = "#e8e8e8"
BIV_HI_LO = "#d03b3b"   # feature 1 high  -> red
BIV_LO_HI = "#2a78d6"   # feature 2 high  -> blue
BIV_HI_HI = "#7a3b9c"   # both high       -> purple


def _hex_to_rgb(h: str) -> Tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _rgb_to_hex(rgb: Sequence[float]) -> str:
    return "#" + "".join(f"{max(0, min(255, round(c))):02x}" for c in rgb)


def _distinct(values: Sequence[str]) -> List[str]:
    seen: set = set()
    order: List[str] = []
    for v in values:
        if v is not None and v not in seen:
            seen.add(v)
            order.append(v)
    return order


def natural_key(s: str) -> list:
    """Sort key that orders '2' before '10' while still sorting text
    alphabetically (case-insensitively)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(s))]


def categorical_colors(values: Sequence[str], dark: bool = False, offset: int = 0) -> Dict[str, str]:
    """Assign hues in fixed order starting at `offset`; past the last hue ->
    'Other'. `offset` lets successive annotation tracks use non-overlapping
    hues so their legends don't collide."""
    idx = 1 if dark else 0
    out: Dict[str, str] = {}
    for i, v in enumerate(_distinct(values)):
        j = i + offset
        out[v] = CATEGORICAL[j][idx] if j < len(CATEGORICAL) else OTHER[idx]
    return out


def sequential_color(t: float, ramp: int = 0) -> str:
    """t in [0,1] -> interpolated colour on one of the single-hue ramps."""
    stops = SEQUENTIAL_RAMPS[ramp % len(SEQUENTIAL_RAMPS)]
    t = 0.0 if t != t else max(0.0, min(1.0, t))  # NaN -> 0
    pos = t * (len(stops) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(stops) - 1)
    f = pos - lo
    a, b = _hex_to_rgb(stops[lo]), _hex_to_rgb(stops[hi])
    return _rgb_to_hex([a[k] + (b[k] - a[k]) * f for k in range(3)])


def bivariate_color(f1: float, f2: float) -> str:
    """f1, f2 in [0,1]. Bilinear blend of the four corner colours."""
    f1 = 0.0 if f1 != f1 else max(0.0, min(1.0, f1))
    f2 = 0.0 if f2 != f2 else max(0.0, min(1.0, f2))
    ll, hl = _hex_to_rgb(BIV_LO_LO), _hex_to_rgb(BIV_HI_LO)
    lh, hh = _hex_to_rgb(BIV_LO_HI), _hex_to_rgb(BIV_HI_HI)
    top = [ll[k] + (hl[k] - ll[k]) * f1 for k in range(3)]       # f2 = 0 edge
    bot = [lh[k] + (hh[k] - lh[k]) * f1 for k in range(3)]       # f2 = 1 edge
    return _rgb_to_hex([top[k] + (bot[k] - top[k]) * f2 for k in range(3)])


def shape_for(values: Sequence[str]) -> Dict[str, str]:
    return {
        v: (SHAPES[i] if i < len(SHAPES) else "dot-small")
        for i, v in enumerate(_distinct(values))
    }
