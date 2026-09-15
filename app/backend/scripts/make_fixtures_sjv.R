#!/usr/bin/env Rscript
# Build backend/fixtures/mini.rds — a junction x sample matrix whose row names
# exercise every classification category against fixtures/mini.gtf (gene TESTG1).
#
#   rownames: chr:start-end:strand
#   sample_A carries the "designed" counts; other samples are noise.
#
# Run:  Rscript backend/scripts/make_fixtures.R

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grep("^--file=", args0)])
here <- if (length(file_arg)) dirname(normalizePath(file_arg)) else "backend/scripts"
out <- file.path(dirname(here), "fixtures", "mini.rds")

rows <- rbind(
  c("chr1:1201-1999:+", 50),   # annotated      (TESTG1.1 intron 1)
  c("chr1:2201-2999:+", 30),   # annotated      (TESTG1.1 / .2 intron)
  c("chr1:2201-3999:+", 20),   # exon_skipping  (skips exon 3 in TESTG1.2 frame)
  c("chr1:1150-1999:+", 8),    # alt_5p         (novel donor, annotated acceptor)
  c("chr1:2201-2500:+", 6),    # alt_3p         (annotated donor, novel acceptor)
  c("chr1:1201-4999:+", 4),    # isoform_switch (donor only in .1, acceptor only in .3)
  c("chr1:2400-2600:+", 12),   # novel_exon     (inside annotated intron 2201-2999)
  c("chr1:1500-3500:+", 3),    # novel          (spans exons, no annotated site)
  c("chr1:3201-3999:+", 0),    # annotated but zero count -> dropped
  c("chr2:1000-2000:+", 99),   # off-locus      -> filtered out
  c("not-a-junction",   7)     # bad rowname    -> skipped with a warning
)

rn <- rows[, 1]
a  <- as.numeric(rows[, 2])
set.seed(1)
b <- pmax(0, round(a * runif(length(a), 0.4, 1.6)))
cc <- pmax(0, round(a * runif(length(a), 0.4, 1.6)))

m <- cbind(sample_A = a, sample_B = b, sample_C = cc)
rownames(m) <- rn
storage.mode(m) <- "double"

saveRDS(m, out)
cat("wrote", out, "\n")
print(m)

# mini_sparse.rds — same content as a Matrix dgCMatrix, to exercise the sparse
# reader (services/rds.py sparse path + scripts/sparse_rds_to_csc.R).
suppressWarnings(suppressMessages(library(Matrix)))
sp_out <- file.path(dirname(here), "fixtures", "mini_sparse.rds")
sp <- as(Matrix(m, sparse = TRUE), "CsparseMatrix")
saveRDS(sp, sp_out)
cat("wrote", sp_out, "\n")
