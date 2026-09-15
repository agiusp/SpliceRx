# SJ Web Apps

Interactive web tools for exploring splice-junction data from TCGA cohorts
prepared with the companion `prepTCGAdata` R package — sashimi plots,
cohort-level novel-splice-junction heatmaps, survival-group classification,
and direct junction lookup, all in one app.

📖 **[Project page / overview](https://GITHUB-USERNAME.github.io/sj_webapps/)**
&nbsp;·&nbsp; run it yourself with the instructions below.

## What's in here

| Tab | What it does |
|---|---|
| **Data** | Point the app at a local directory of prepared TCGA cohort files; this is the only way to load data — no per-tab file upload. |
| **Sashimi plot** (SJV) | SCANVIS-style sashimi plot for one sample or a sample group. |
| **NSJCG** (SJVC) | Cohort-level heatmap / PCA / UMAP of novel-splice-junction-per-gene counts, with clinical variables layered on. |
| **SJSurv** | Cross-validated logistic-regression classification of a `Good` / `Poor` survival-group label from a splice-junction feature matrix. |
| **SJ Lookup** | Direct per-junction lookup against junction metadata. |

The actively developed app lives in [`app/`](app) — one FastAPI backend, one
React/Vite frontend, four tabs. [`sjv/`](sjv) and [`sjvc/`](sjvc) are the
original standalone sashimi-plot and NSJCG apps that `app/` was merged from;
they're kept here for reference and still run on their own. See
[`app/README.md`](app/README.md) for the full architecture writeup.

## Running it locally

You'll need **Python 3.10+** and **Node.js 18+**. (R is only needed if you
want to regenerate test fixtures or prepare your own TCGA cohort data with
`prepTCGAdata` — not to run the app itself.)

### 1. Clone the repo

```bash
git clone https://github.com/GITHUB-USERNAME/sj_webapps.git
cd sj_webapps/app
```

### 2. Start the backend (FastAPI, port 8000)

```bash
cd backend
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# SJ_DATA_ROOT bounds which directory the Data tab is allowed to read from
# (defaults to $HOME if unset)
SJ_DATA_ROOT="$HOME/path/to/your/tcga-cohort-data" SJV_CACHE_DIR=.cache \
  .venv/bin/python -m uvicorn main:app --port 8000 --reload
```

### 3. Start the frontend (Vite, port 5173) — in a second terminal

```bash
cd app/frontend
npm install
npm run dev
```

Open **<http://localhost:5173>**. The frontend proxies `/api` requests to
`http://localhost:8000` (edit `frontend/vite.config.ts` if your backend runs
elsewhere). `#data` / `#sjv` / `#sjvc` / `#sjsurv` / `#sjlookup` in the URL
deep-link a tab; the last-used tab is also remembered in `localStorage`.

### Getting cohort data to point the app at

The Data tab expects the outputs of the `prepTCGAdata` R package — files
like `*_junction_counts.rds`, `*_novel_junction_counts_per_gene.rds`,
`*_sample_metadata.rds`, etc. Prepare a cohort with that package first, then
point `SJ_DATA_ROOT` at the directory containing it (or a parent of it) and
select the cohort directory from the Data tab.

### Running tests

```bash
cd app/backend
SJV_CACHE_DIR=.cache_test .venv/bin/python -m pytest tests/sjv
SJVC_CACHE_DIR=.cache_test .venv/bin/python -m pytest tests/sjvc
.venv/bin/python -m pytest tests/sjsurv tests/dataload
```

## Repo layout

```
sj_webapps/
├── app/    combined app — one backend, one frontend, all tabs (start here)
├── sjv/    standalone sashimi-plot app (reference / kept runnable)
└── sjvc/   standalone NSJCG app (reference / kept runnable)
```
