# SJV — Splice Junction Visualizer

Reproduces a [SCANVIS](https://pmc.ncbi.nlm.nih.gov/articles/PMC6853764/)-style
sashimi plot from user data:

1. a **splice-junction count matrix** (`.rds`; rows named `chr:start-end:strand`,
   columns = samples),
2. a **GENCODE reference** (choose a release, or upload a custom `.gtf` / `.gtf.gz`),
3. a **series** — either one sample, or a **sample group**: upload any
   `.csv` / `.tsv` / `.rds` that has a **`sample_id` column** matching the matrix
   columns, choose one of the other columns to stratify by, then choose a value
   of it. The plot then uses the per-junction **median** across that group's
   samples, excluding samples with NA or 0 for the junction,
4. one or more **query genes** (the box autocompletes as you type; each name
   becomes a chip). With more than one gene the plot spans the bounding box of
   all their loci and draws every junction overlapping *any* of them — so a
   read-through junction between two neighbouring genes (e.g. `MUC21` +
   `MUC22`) is shown. All genes must be on one chromosome. Junctions are still
   classified strand- and annotation-aware, against the **combined** transcript
   models of every gene in the plot,
5. optionally a **minimum read count** — junctions below it are hidden (default
   `0`, i.e. every junction with non-zero support is shown).

It draws one arc per junction that falls in the gene locus (or in any of the
loci) — all arcs at the same line width, **arc height ∝ log2(read count)** (or
log2 group median) as the read-support cue, color by how the junction relates to
the annotation — over an IGV-style transcript track (one merged model per gene
when transcripts are collapsed). Annotated (grey) arcs are drawn underneath; the
colored arcs are drawn on top in decreasing order of read support so the smaller
ones stay legible. Each unannotated arc is labelled with its read count, and the
tallest arc with the region maximum (over every junction in the figure).

The figure is condensed to fit the screen by default; untick **Fit width** for
a natural-scale figure that scrolls horizontally, with `−` / `+` to zoom.

## Layout

```
backend/    FastAPI service: RDS parsing, GENCODE indexing, gene lookup,
            junction filtering + classification + scaling, plot JSON
frontend/   React + Vite + TypeScript; SVG sashimi renderer
```

## Run it

**Backend**

```bash
cd backend
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
SJV_CACHE_DIR=.cache .venv/bin/python -m uvicorn app.main:app --port 8000 --reload
```

**Frontend** (separate terminal)

```bash
cd frontend
npm install
npm run dev            # http://localhost:5173  (proxies /api to :8000)
```

Then open http://localhost:5173, upload a matrix, pick a reference, choose a
sample and type a gene.

### Pre-indexing a release from a local GTF

First use of a full release downloads + indexes it (~6 s per GB once parsed). To
prime one from a file you already have:

```bash
cd backend
SJV_CACHE_DIR=.cache .venv/bin/python scripts/prebuild_index.py \
    human v29 ~/Work/SJ/Data/gencode.v29.annotation.gtf
```

This builds the SQLite index now and records the path in
`backend/gencode_sources.json` (git-ignored, machine-specific), so the app uses
the local file instead of downloading — while still accepting any other
`v<N>` / `vM<N>` release on demand.

### Test fixtures

`backend/fixtures/mini.rds` + `backend/fixtures/mini.gtf` describe a synthetic
gene `TESTG1` whose junctions cover every classification category. Regenerate the
RDS with `Rscript backend/scripts/make_fixtures.R` (needs R).

## Tests

```bash
cd backend && SJV_CACHE_DIR=.cache_test .venv/bin/python -m pytest
```

## Junction classification

Each drawn junction is one of:

| category | meaning | color |
|---|---|---|
| `annotated` | matches an annotated intron | grey |
| `exon_skipping` | both splice sites known; skips ≥1 exon within a transcript | blue |
| `alt_5p` | 3′ site known, 5′ (donor) site novel | orange |
| `alt_3p` | 5′ site known, 3′ (acceptor) site novel | aqua |
| `isoform_switch` | both sites known but never combined in one transcript | yellow |
| `novel_exon` | neither site known; junction lies inside an annotated intron | green |
| `novel` | neither site known; not inside an intron | violet |

Splice sites are strand-aware. Colors are the `dataviz` skill's validated
colorblind-safe categorical palette; the legend + hover labels carry identity so
color is never the only cue.

## Notes / v1 limitations

- GENCODE releases are downloaded once and indexed into SQLite
  (`~/.cache/sjv/`, or `$SJV_CACHE_DIR`). First use of a full release takes ~1 min
  to build the index. (The plan's tabix path is deferred — no system tabix here.)
- RDS is read with `pyreadr`, falling back to an `Rscript` subprocess for base
  matrices / `Matrix` sparse objects. The raw upload is deleted immediately after
  parsing; derived arrays live in an in-memory session cache with a 2 h TTL
  (single-process only).
- Genome build (GRCh37 vs 38) is **not** auto-verified — pick the release that
  matches your junction coordinates. Chromosome naming (`chr1` vs `1`) is
  normalised and mismatches are surfaced as a warning.

## Deferred (see `~/.claude/plans/reactive-mapping-engelbart.md`)

Dashed line-type for frame-shift status; inverted read-coverage profile;
variant/SNV dots; multi-sample comparison.
