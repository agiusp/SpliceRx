"""Single source of truth for junction-category colors.

Hues are the dataviz skill's validated colorblind-safe categorical palette
(references/palette.md). Grey is reserved for annotated junctions. The frontend
mirror lives at frontend/src/plot/palette.ts and must be kept in sync.
"""
from __future__ import annotations

# category -> (light hex, dark hex)
CATEGORY_COLORS = {
    "annotated": ("#898781", "#898781"),
    "exon_skipping": ("#2a78d6", "#3987e5"),
    "alt_5p": ("#eb6834", "#d95926"),
    "alt_3p": ("#1baf7a", "#199e70"),
    "isoform_switch": ("#eda100", "#c98500"),
    "novel_exon": ("#008300", "#008300"),
    "novel": ("#4a3aa7", "#9085e9"),
}

CATEGORY_LABELS = {
    "annotated": "Annotated",
    "exon_skipping": "Exon skipping",
    "alt_5p": "Alt 5' splice site",
    "alt_3p": "Alt 3' splice site",
    "isoform_switch": "Isoform switch",
    "novel_exon": "Novel exon",
    "novel": "Novel / unknown",
}

# Order used for the legend and for deterministic rendering.
CATEGORY_ORDER = [
    "annotated",
    "exon_skipping",
    "alt_5p",
    "alt_3p",
    "isoform_switch",
    "novel_exon",
    "novel",
]


def light_color(category: str) -> str:
    return CATEGORY_COLORS.get(category, CATEGORY_COLORS["novel"])[0]


def legend_entries(categories):
    """Return legend rows (category, label, light color) for the given
    categories, in canonical order, de-duplicated."""
    present = set(categories)
    rows = []
    for cat in CATEGORY_ORDER:
        if cat in present:
            rows.append(
                {
                    "category": cat,
                    "label": CATEGORY_LABELS[cat],
                    "color": light_color(cat),
                }
            )
    return rows
