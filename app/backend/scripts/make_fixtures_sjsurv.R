#!/usr/bin/env Rscript
# Build SJSurv test fixtures:
#   tests/sjsurv/fixtures/mini_sjdat_genes.rds   — gene x sample "sjdat" matrix
#   tests/sjsurv/fixtures/mini_metadata.csv      — RAW sample metadata: the
#                                                   columns prepTCGAdata::
#                                                   get_tcga_data() writes
#                                                   (histology / stage / age /
#                                                   survival) — no Group or
#                                                   SurviverGroup; SJSurv
#                                                   derives those itself now.
#   tests/sjsurv/fixtures/mini_metadata.rds      — same, as a data.frame
#   tests/sjsurv/fixtures/mini_metadata_missing_col.csv — drops `survival_days`
#                                                   (error path)
#
# 80 samples, one Histology (so Group cardinality is driven by Stage x age band
# only), half Stage I / half Stage III, ages spanning ~22-88 so the classic
# (30-50]/(50-70]/(>70) bands fragment the cohort into small groups (some ages
# <=30 fall out entirely) — the scenario that motivated making the age bands
# adjustable in the app instead of fixed in the R package.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grep("^--file=", args0)])
here <- if (length(file_arg)) dirname(normalizePath(file_arg)) else "backend/scripts"
fx <- file.path(dirname(here), "tests", "sjsurv", "fixtures")
dir.create(fx, recursive = TRUE, showWarnings = FALSE)

set.seed(11)
n <- 80
samples <- sprintf("S%02d", seq_len(n))

n_genes <- 60
genes <- sprintf("GENE%03d", seq_len(n_genes))
M <- matrix(rpois(n_genes * n, lambda = 3), nrow = n_genes, dimnames = list(genes, samples))
M[1:12, ] <- M[1:12, ] + matrix(rpois(12 * n, lambda = 4), nrow = 12)  # some variable rows
M[58:60, ] <- 0                                                        # dead rows (coverage filter)
storage.mode(M) <- "double"
saveRDS(M, file.path(fx, "mini_sjdat_genes.rds"))

meta <- data.frame(
  sample_id = samples,
  tcga.cgc_case_histological_diagnosis = "Adenocarcinoma",
  tcga.cgc_case_pathologic_stage = rep(c("Stage I", "Stage IIIA"), length.out = n),
  age_at_diagnosis_years = round(runif(n, 22, 88), 1),
  survival_days = round(runif(n, 100, 4000)),
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(fx, "mini_metadata.csv"), row.names = FALSE, na = "")
saveRDS(meta, file.path(fx, "mini_metadata.rds"))

missing_col <- meta[, setdiff(colnames(meta), "survival_days")]
write.csv(missing_col, file.path(fx, "mini_metadata_missing_col.csv"), row.names = FALSE)

n_over30 <- sum(meta$age_at_diagnosis_years > 30)
cat(sprintf("wrote fixtures to %s (%d/%d samples have age > 30)\n", fx, n_over30, n))
