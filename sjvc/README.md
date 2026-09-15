# SJVC — Splice Junction Visualizer with Clinical

A cohort-level companion to [SJV](../sjv). Upload a splice-junction count matrix
and a clinical table, pick a **gene set**, and explore the samples as a
**heatmap** or a **2D projection** (PCA / UMAP) built from the junctions in those
genes — with clinical variables layered on.

## Pipeline

```
junction matrix (.rds)  +  clinical table (.csv/.tsv/.rds, needs a sample_id column)
        +  GENCODE reference (release or custom GTF)
                      │
   gene set  ─ typed genes (exact, or "match as name prefix": MUC → MUC1, MUC2…)
             / uploaded list / a pathway (MSigDB/KEGG/Reactome/…)
             / top-N by MAD — ranks every row of the matrix by descending
               median absolute deviation; no GENCODE resolution involved,
               executes only when you hit **Run** in the View step (nothing
               fires while you're still typing N)

   If the uploaded matrix is gene-level rather than junction-level — row
   names are gene symbols (e.g. the output of
   functions/condense_junctions_by_gene.py), or the row labels sit in
   leading `gene_id` / `gene_name` columns of a gene-count table — the app
   detects that on upload. GENCODE then isn't needed: typed genes / uploaded
   list / pathway match straight against the row labels (preferring an exact
   single-gene row, falling back to comma-joined rows that contain the gene),
   and "top-N by MAD" ranks the rows directly. GENCODE stays optional — only
   the "protein-coding genes only" MAD option uses it, to restrict the
   ranking to genes with gene_type "protein_coding" in the reference.
                      ▼
   GENCODE resolves genes → loci → junctions in the matrix overlapping them
   → the feature matrix is built immediately    (skipped for the top-N-by-MAD
     path and for a gene-level matrix, which select rows directly)
                      ▼
   preprocess: NA→0, drop constant features, log1p, z-score
                      ▼
   ┌─ Heatmap ────────────────┐   ┌─ PCA / UMAP ─────────────┐
   │ hierarchical clustering  │   │ 2D scatter               │
   │ (rows always; columns    │   │ 1 feature  → colour      │
   │  cluster OR group-by a    │   │ 2 categorical → colour + │
   │  categorical feature)     │   │   point shape            │
   │ clinical colour bars on   │   │ 2 numeric → red↔blue     │
   │  top of the columns       │   │   bivariate blend        │
   │ compact toggle; copy /    │   │ mixed → colour + shape   │
   │  download PNG or SVG      │   └──────────────────────────┘
   └───────────────────────────┘
```

The cell grid is fitted into a slide-friendly box (a 400-sample cohort no
longer renders thousands of px wide), with cells going sub-pixel for a very
large matrix. The toolbar has a **compact** toggle (squeezes further, drops
labels) and **Copy image** / **Download PNG** / **Download SVG** — the PNG/SVG
bundle the cell grid, dendrograms, clinical bars, and labels into one file.
The projection has **Download SVG**.

In "top by MAD" mode a **download data (CSV)** link gives the feature matrix
(features × samples, raw values before preprocessing) that the plots are
built from — `GET /api/session/{id}/features.csv`.

## Run it

```bash
# backend  (port 8001)
cd backend
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
SJVC_CACHE_DIR=.cache .venv/bin/python -m uvicorn app.main:app --port 8001 --reload

# frontend (port 5174, proxies /api -> :8001)
cd frontend && npm install && npm run dev
```

The GENCODE cache is shared with SJV by default (`~/.cache/sjv`), so a release
indexed there — e.g. the prebuilt v29 — is available here immediately.
`umap-learn` is imported lazily, so the PCA / heatmap paths stay fast; the first
UMAP call pays a one-time JIT warm-up (~5–8 s).

To work with genes instead of junctions, pre-condense the junction matrix
**offline** with `~/Work/SJ/functions/condense_junctions_by_gene.py` (see that
directory's README) rather than in the app — condensing genome-wide on every
request was expensive, so the in-app "condense to one row per gene" step has
been removed for now. Upload that condensed `.rds` (or any gene-count matrix
— gene symbols as row names, or in leading `gene_id` / `gene_name` columns)
as the junction matrix: the app treats it as gene-level, and typed genes /
list / pathway / top-N-by-MAD all operate on the gene rows directly. GENCODE
stays optional — pick a release only for the "protein-coding genes only" MAD
toggle, which drops non-coding rows (by the reference's gene_type).

## Clinical column handling

- Shown to the user: **all numeric columns**, plus **categorical columns whose
  distinct-value count is below 0.9 × the row count** (drops IDs / near-unique
  free-text). Chosen from a dropdown; picks become removable chips.
- Typing is auto-detected (numeric dtype + many distinct → numeric; text or
  ≤ 10 distinct → categorical) and **overridable per feature** — click the
  `numeric`/`categorical` chip on a selected feature to flip it.
- Projections take **1 or 2** features; the heatmap takes any number.

## Pathway libraries

MSigDB Hallmark ships in `backend/data/gmt/`; KEGG, Reactome, WikiPathways and
GO-BP are fetched once from Enrichr and cached under `~/.cache/sjv/gmt`. Human
symbols only — pathway mode is disabled for mouse.

## Tests

```bash
cd backend && SJVC_CACHE_DIR=.cache_test .venv/bin/python -m pytest
```

Fixtures (`backend/fixtures/`) describe a 12-sample cohort split into two groups
with a planted junction signature (`TESTG1` up in one group, `TESTG2` in the
other) and a matching `mini_clinical.csv`. Regenerate with
`Rscript backend/scripts/make_fixtures.R`.

## Reused from SJV

`gencode.py` (verbatim, shared cache), the `.rds` reader, the session store, the
`sample_id`-column table parser, junction row-name parsing, the validated colour
palette, and the frontend autocomplete / upload patterns.

## Deferred

Survival endpoints (KM / Cox) from a time/event pair in the clinical table;
enrichment testing on the selected set; saved views; single-origin Docker deploy.
