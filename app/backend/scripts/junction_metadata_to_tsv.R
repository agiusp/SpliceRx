#!/usr/bin/env Rscript
# Project a prepTCGAdata junction-metadata data.frame down to just the
# columns the junction-gene index needs, before it ever reaches Python.
#
# pyreadr parses every column of this data.frame (dozens of TCGA cohorts'
# worth of annotation, some list-like) and is dramatically slower than R's
# own readRDS + a plain TSV write for a table this wide — for an 8.8M-row
# cohort, ~100s+ via pyreadr vs. ~20s this way.
#
# Usage: junction_metadata_to_tsv.R <in.rds> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: junction_metadata_to_tsv.R <in.rds> <out.tsv>")

x <- readRDS(args[1])
if (!is.data.frame(x)) stop("junction metadata RDS is not a data.frame")

wanted <- c("seqnames", "start", "end", "strand", "gencode_gene_id", "gencode_gene_name")
have <- intersect(wanted, colnames(x))
if (length(have) == 0L) stop(paste("none of the expected columns are present:", paste(wanted, collapse = ", ")))

write.table(x[, have, drop = FALSE], args[2], sep = "\t", quote = FALSE, row.names = FALSE)
