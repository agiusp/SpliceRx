# SpliceRx — combined

One app, four tabs:

| Tab | What it is | Was |
|---|---|---|
| **Data** | point at a prepared TCGA cohort directory; the files the app can use are loaded straight into the other tabs. **This is the only way to load cohort data** — no app tab has a file upload any more. | new |
| **Sashimi plot** | [SJV](../sjv) — SCANVIS-style sashimi plot for one sample or a sample group | its own app on :8000 / :5173 |
| **NSJCG** (Novel Splice Junction Counts per Gene) | [SJVC](../sjvc) — cohort heatmap / PCA / UMAP with clinical variables layered on. Analytics run on the **gene matrix only**: `*_novel_junction_counts_per_gene.rds`. | its own app on :8001 / :5174 |
| **SJSurv** (Survivor groups) | `backend/sjsurv/` — cross-validated classification of the `SurviverGroup` (Good / Poor) label from a splice-junction feature matrix. | new |

`NSJCG` is the display name; the code, routes (`/api/sjvc/*`), CSS scope
(`.sjvc-scope`) and directories still use `sjvc`.

Each app keeps its own code and its own in-memory session store; cohort data
comes in only through the Data tab. They are merged only at the edges:

- **Backend** — one FastAPI process (`backend/main.py`) that mounts the two
  routers under distinct prefixes: `/api/sjv/*` and `/api/sjvc/*` (plus
  `/api/health`). Each app's package (`backend/sjv/`, `backend/sjvc/`) is the
  original `app/` tree unchanged apart from the router no longer hard-coding its
  `/api` prefix. `backend/sjv/main.py` / `backend/sjvc/main.py` still build the
  standalone single-app servers (used by the tests).
- **Frontend** — one Vite app. `src/Shell.tsx` is a tab switch that mounts
  `<SjvApp/>` and `<SjvcApp/>` (both stay mounted, so in-progress work on the
  other tab survives a switch). Each app's stylesheet is wrapped in a
  `.sjv-scope` / `.sjvc-scope` selector so their identically-named classes
  (`.app`, `.panel`, `.row`, …) don't collide; the shared `:root` design tokens
  live in `src/styles/base.css`. Each `src/<app>/api.ts` calls its own
  `/api/sjv` or `/api/sjvc` prefix.

- **Data tab** — `backend/dataload/` (mounted at `/api/dataload`). `POST /scan`
  classifies the files in a directory (by name + a 250 MB size guard, without
  opening the big ones); the `POST /sjv/{sid}/*`, `/sjvc/{sid}/*` and
  `/sjsurv/{sid}/*` endpoints copy a file into the session tmp dir and run the
  same loaders the upload endpoints use (`sjv.api.routes.ingest_rds` /
  `ingest_sample_metadata`, `sjvc.api.routes.ingest_junctions` /
  `ingest_clinical`, `sjsurv.api.routes.ingest_sjdat` / `ingest_metadata`).
  Paths are restricted to `$SJ_DATA_ROOT` (default: `$HOME`). The Shell owns all
  three app session ids (persisted in `localStorage`), so a cohort loaded here
  shows up in the app tab and survives a page refresh
  (`GET /api/<app>/session/{sid}/state`).
  **The GENCODE reference is chosen once on the Data tab** (`POST
  /api/dataload/gencode` / `.../gtf`), applied to both app sessions at once — the
  SJV / SJVC tabs no longer have their own release picker, they just show the
  chosen reference (hydrated from the session state). The GENCODE index cache is
  shared, so the release is downloaded once even though each app keeps its own
  SQLite schema.

- **Sparse junction matrices** — `backend/sjv/services/rds.py` reads a `Matrix`
  sparse RDS (`dgTMatrix` / `dgCMatrix`) without densifying it:
  `scripts/sparse_rds_to_csc.R` dumps the CSC vectors as raw binary and Python
  rebuilds a `scipy.sparse.csc_matrix`. Plot requests slice to the gene locus
  (a few hundred rows) before densifying, so an 8-million-junction TCGA matrix
  plots as fast as a small one. First load takes ~50 s (R read + row-name parse)
  and is cached under `$SJV_CACHE_DIR/matrices/<sha1>` (~5 s on a cache hit).
  Steady RSS ≈ 2–3 GB per loaded cohort. Unstranded (`:*`) junctions are dropped
  with a warning (unchanged behaviour).

The GENCODE cache is still `~/.cache/sjv` (or `$SJV_CACHE_DIR` / `$SJVC_CACHE_DIR`),
shared exactly as before.

- **SJSurv** — `backend/sjsurv/` (mounted at `/api/sjsurv`), its own session
  store. A session holds the *raw* sample metadata plus up to three `sjdat`
  feature matrices (`junction_counts` / `rrs_scores` / `gene_matrix`).
  Stratifying the cohort — histology / pathologic stage / an age band, into a
  `Group`, then a `Good`/`Poor` `SurviverGroup` label from each group's median
  survival — used to be baked into the file at download time by a
  `surv_cohort()` R function; that function is gone from the `prepTCGAdata`
  package, and the equivalent computation now lives here
  (`services/metadata.py`), run interactively:

  - `POST /metadata` (or the Data-tab loader) reads the raw columns
    (`tcga.cgc_case_histological_diagnosis`, `tcga.cgc_case_pathologic_stage`,
    `age_at_diagnosis_years`, `survival_days` — exactly what
    `get_tcga_data()` writes) and immediately stratifies with the classic
    fixed age bands `(30-50], (50-70], (>70)` and `min_group_n=20`, so the
    pipeline works with no extra step. Column resolution deliberately never
    matches the *old* `Histology`/`Stage`/`Age_at_diagnosis`/`Group`/
    `SurviverGroup` names — a cohort downloaded before the R-side removal can
    still carry that stale, already-binned output alongside the raw columns.
  - `POST /age-bands/suggest` computes a data-driven alternative: `n_bands`
    equal-count (quantile) bands spanning the cohort's own age distribution,
    so nothing is left unbanded the way ages ≤30 are under the classic bands.
  - `POST /stratify` recomputes `Group`/`MedianSurvival`/`SurviverGroup` for a
    chosen set of age-band edges + `min_group_n`, invalidating any downstream
    feature selection / model. This is the fix for a cohort (e.g. PAAD) where
    the fixed 3 age bands × several histologies/stages fragments every group
    below `min_group_n` and nobody ends up labelled — coarsening the age bands
    (fewer, wider bands) directly grows the groups.
  - `GroupCount`s are always returned with the `"(no Group)"` bucket (samples
    with no histology, stage, or age band) pinned to the bottom of the list.

  The user then picks a `Group` (or all labelled samples) and which `sjdat` is
  active. `POST /select` keeps features whose non-zero-and-`≥X` coverage over
  the group is `≥ N` (a count, or a fraction when `N < 1`), then the top `n` by
  MAD. `POST /cross-validate` runs stratified `ncv`-fold logistic regression
  and reports out-of-fold ROC AUC / accuracy / confusion; `POST /model` fits
  one model on every selected sample (resubstitution) and returns the features
  ranked by `|weight|`, keeping the fitted model on the session as `MODEL`.
  Junction matrices are read (and cached) through the SJV sparse reader; the
  coverage filter runs on the sparse matrix so only the passing rows are ever
  densified. `rrs_scores` is not a `prepTCGAdata` output yet — the option is
  there and works the moment such a file exists.

### What the Data tab can load (TCGA cohorts from `prepTCGAdata`)

| File | Loads into | Notes |
|---|---|---|
| `*_junction_counts.rds` | Sashimi plot · SJSurv | junction-level matrix; first load of a big sparse one takes ~1 min to index, then cached |
| `*_novel_junction_counts_per_gene.rds` | NSJCG · SJSurv | gene-level feature matrix — the only matrix NSJCG analyses |
| `*_novel_junction_RRS_scores.rds` | SJSurv | relative-read-support score matrix, an sjdat choice |
| `*_sample_metadata.rds` | any | its `sample_id` matches the matrix columns directly; SJSurv reads the *raw* histology/stage/age/survival columns and stratifies them itself (see above) |
| `*_MSI.rds` | NSJCG · SJSurv | ⚠️ sample IDs are TCGA barcodes — only load if they line up |
| `*_junction_metadata.rds` | — | not read by the app |

## Run it

```bash
# backend — one process, port 8000
cd backend
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
# SJ_DATA_ROOT bounds where the Data tab may read files from (default: $HOME)
SJ_DATA_ROOT="$HOME/Work/SJ.Sep2026/Data" SJV_CACHE_DIR=.cache \
  .venv/bin/python -m uvicorn main:app --port 8000 --reload

# frontend — one dev server, port 5173 (separate terminal)
cd frontend
npm install
npm run dev            # http://localhost:5173  (proxies /api -> :8000)
```

Open <http://localhost:5173>. `#data` / `#sjv` / `#sjvc` in the URL deep-links a
tab; the last-used tab is also remembered in `localStorage`.

### Pointing the dev proxy elsewhere

`frontend/vite.config.ts` proxies `/api` to `http://localhost:8000`. Edit that
line if the backend runs on another port.

## Tests

```bash
cd backend
SJV_CACHE_DIR=.cache_test  .venv/bin/python -m pytest tests/sjv
SJVC_CACHE_DIR=.cache_test .venv/bin/python -m pytest tests/sjvc
.venv/bin/python -m pytest tests/sjsurv
.venv/bin/python -m pytest tests/dataload
# or everything: SJV_CACHE_DIR=.cache_test SJVC_CACHE_DIR=.cache_test pytest tests
```

Fixtures moved to `backend/tests/<app>/fixtures/`. Regenerate with
`Rscript backend/scripts/make_fixtures_sjv.R` /
`…_sjvc.R` / `…_sjsurv.R` (needs R).

## Follow-ups / not done yet

- **SJSurv predictive power is only as good as the labels.** The classifier is
  L2 logistic regression on `log1p` features; swap in a non-linear estimator
  (RF / gradient boosting) behind the same API if the linear model underfits.
  No streaming progress — the CV messages are returned in one shot.
- **RRS scores** aren't a `prepTCGAdata` output yet; SJSurv already accepts a
  `*_RRS*` matrix if one turns up.
- **Reuse the saved SJSurv `MODEL`** — it is kept on the session but nothing
  reads it back yet (score a held-out cohort, inspect it, etc.).

- **Junction-level NSJCG on a TCGA cohort.** NSJCG currently analyses the
  gene matrix only. SJV now reads sparse junction matrices; the `sjvc`
  junction-level path (`sjvc/services/rds.py`, a separate `Matrix` class) still
  doesn't — it would reuse the same Rscript. Re-enabling it also means bringing
  back a way to load a junction matrix into the tab.
- **Full rename `sjvc` → `nsjcg`.** Only the display name changed; the Python
  package, API prefix, CSS scope and frontend directory are still `sjvc`.
- **`MSI.rds` sample-ID bridge** (TCGA barcodes → GDC file IDs).
- The in-RAM CSC could be memory-mapped (lower RSS, more I/O per query) if many
  large cohorts need to be open at once.
- `gencode.py` / `classify.py` are still duplicated per package.
- No single-process production build (FastAPI serving the built frontend).
- The original `../sjv` and `../sjvc` trees are left in place as reference.
