"""SJ Lookup — direct per-junction lookup against a cohort's junction metadata.

Loads the same ``TCGA_<cohort>_junction_metadata.rds`` file the 2D View /
SJSurv "fast gene lookup" path uses (``sjvc.services.junction_metadata``),
but the other direction: given one or more splice junctions
(``chr:start-end:strand``), return whatever that table records for each —
gene overlap, but also the annotation/novelty and splice-motif columns the
gene-index path never surfaces.

No sample metadata, no matrix, no model — one small table, looked up by key.
Only separated from ``sjv`` / ``sjvc`` / ``sjsurv`` by URL prefix
(``/api/sjlookup``); it reuses their RDS-reading machinery.
"""
