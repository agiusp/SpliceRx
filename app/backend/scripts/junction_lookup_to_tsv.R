#!/usr/bin/env Rscript
# Project a prepTCGAdata junction-metadata data.frame down to the columns
# SJ Lookup needs — a wider set than junction_metadata_to_tsv.R's (that one
# only serves the gene-index path and drops the annotation/motif columns a
# per-junction lookup wants to show).
#
# Usage: junction_lookup_to_tsv.R <in.rds> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: junction_lookup_to_tsv.R <in.rds> <out.tsv>")

x <- readRDS(args[1])
if (!is.data.frame(x)) stop("junction metadata RDS is not a data.frame")

wanted <- c("seqnames", "start", "end", "strand", "width", "annotated",
            "left_motif", "right_motif", "left_annotated", "right_annotated",
            "gencode_gene_id", "gencode_gene_name")
have <- intersect(wanted, colnames(x))
if (length(have) == 0L) stop(paste("none of the expected columns are present:", paste(wanted, collapse = ", ")))

write.table(x[, have, drop = FALSE], args[2], sep = "\t", quote = FALSE, row.names = FALSE)
