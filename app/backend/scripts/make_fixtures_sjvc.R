#!/usr/bin/env Rscript
# Build SJVC test fixtures:
#   fixtures/mini_junctions.rds  — junction x sample count matrix (12 samples)
#   fixtures/mini_genes.rds      — the same, pre-condensed to one row per gene
#                                  (gene symbols as row names, not coordinates)
#   fixtures/mini_clinical.csv   — clinical table keyed by sample_id
#   fixtures/mini_clinical.rds   — same, as a tibble with leading patient_id /
#                                  specimen_id columns (specimen_id repeats),
#                                  mimicking a TCGA sample-metadata layout
# Junctions live in TESTG1 / TESTG2 of fixtures/mini.gtf, with a planted split:
#   samples s01..s06 have stronger TESTG1 signal, s07..s12 stronger TESTG2.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grep("^--file=", args0)])
here <- if (length(file_arg)) dirname(normalizePath(file_arg)) else "backend/scripts"
fx <- file.path(dirname(here), "fixtures")

set.seed(7)
samples <- sprintf("s%02d", 1:12)
grpA <- 1:6

g1 <- c("chr1:1201-1999:+", "chr1:2201-2999:+", "chr1:3201-3999:+", "chr1:4201-4999:+",
        "chr1:1500-2800:+", "chr1:2400-3600:+")
g2 <- c("chr1:10301-10999:+", "chr1:11301-11999:+", "chr1:12301-12999:+", "chr1:10500-11800:+")
extra <- c("chr2:5000-6000:+", "not_a_junction")
rn <- c(g1, g2, extra)

M <- matrix(0, nrow = length(rn), ncol = length(samples), dimnames = list(rn, samples))
pois <- function(lam) rpois(length(samples), lam)

for (j in g1) {
  base <- pois(ifelse(seq_along(samples) %in% grpA, 40, 8))
  M[j, ] <- base
}
for (j in g2) {
  base <- pois(ifelse(seq_along(samples) %in% grpA, 6, 35))
  M[j, ] <- base
}
M["chr2:5000-6000:+", ] <- pois(50)
M["not_a_junction", ]   <- pois(5)

# sprinkle a few dropouts / NAs
M["chr1:2400-3600:+", c(2, 5)] <- 0
M["chr1:11301-11999:+", 9] <- NA
storage.mode(M) <- "double"
saveRDS(M, file.path(fx, "mini_junctions.rds"))

# gene-level version: C[gene, sample] = # of that gene's junctions with count > 0
G <- rbind(
  TESTG1 = colSums(!is.na(M[g1, ]) & M[g1, ] > 0),
  TESTG2 = colSums(!is.na(M[g2, ]) & M[g2, ] > 0),
  OTHERG = colSums(!is.na(M["chr2:5000-6000:+", , drop = FALSE]) & M["chr2:5000-6000:+", , drop = FALSE] > 0)
)
storage.mode(G) <- "double"
saveRDS(G, file.path(fx, "mini_genes.rds"))

clin <- data.frame(
  sample_id = samples,
  subtype   = c(rep("LUAD", 6), rep("LUSC", 6)),
  stage     = rep(c(1, 2, 3), length.out = 12),
  age       = c(52, 61, 47, 70, 66, 58, 73, 49, 64, 55, 78, 45),
  tmb       = round(runif(12, 1, 20), 2),
  response  = rep(c("CR", "PR", "SD", "PD"), length.out = 12),
  notes     = sprintf("case-%s-note", samples),
  stringsAsFactors = FALSE
)
write.csv(clin, file.path(fx, "mini_clinical.csv"), row.names = FALSE, quote = FALSE)

# TCGA-style: leading patient_id / specimen_id columns, one numeric column, then
# sample_id + features. specimen_id deliberately repeats (2 samples per patient)
# so a reader that (wrongly) promoted it to a de-duplicated row index would drop
# half the rows.
clin_rds <- data.frame(
  patient_id  = sprintf("P%02d", rep(1:6, each = 2)),
  specimen_id = sprintf("P%02d", rep(1:6, each = 2)),
  age_years   = clin$age,
  sample_id   = samples,
  subtype     = clin$subtype,
  stage       = clin$stage,
  stringsAsFactors = FALSE
)
saveRDS(clin_rds, file.path(fx, "mini_clinical.rds"))

cat("wrote mini_junctions.rds, mini_genes.rds, mini_clinical.csv, mini_clinical.rds\n")
print(M)
print(G)
print(clin)
print(clin_rds)
