#!/usr/bin/env Rscript
# Fallback RDS reader: densify a junction x sample matrix to a TSV.
# Usage: rds_to_tsv.R <in.rds> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: rds_to_tsv.R <in.rds> <out.tsv>")

x <- readRDS(args[1])

# Handle Matrix package sparse matrices and plain matrices / data.frames.
if (methods::is(x, "Matrix")) x <- as.matrix(x)
if (is.data.frame(x)) x <- as.matrix(x)
if (!is.matrix(x)) stop("RDS object is not a matrix / data.frame")

if (is.null(rownames(x))) stop("matrix has no row names (need chr:start-end:strand)")
if (is.null(colnames(x))) stop("matrix has no column names (need sample names)")

df <- data.frame(rowname = rownames(x), x, check.names = FALSE, stringsAsFactors = FALSE)
write.table(df, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
