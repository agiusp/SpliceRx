#!/usr/bin/env Rscript
# Dump a sparse Matrix RDS as CSC vectors, WITHOUT densifying it.
#
#   sparse_rds_to_csc.R <in.rds> <out_prefix>
#
# Writes:  <pre>.indices  int32   row indices   (dgCMatrix @i)
#          <pre>.indptr   int32   column pointers (dgCMatrix @p, length ncol+1)
#          <pre>.data     float64 non-zero values (dgCMatrix @x)
#          <pre>.rows     text    row names, one per line
#          <pre>.cols     text    column names, one per line
#          <pre>.shape    text    "<nrow> <ncol>"
#
# Exit code 3 = the object is a plain (dense) matrix / data.frame; the caller
# should use its dense reader instead.
suppressWarnings(suppressMessages(library(Matrix)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: sparse_rds_to_csc.R <in.rds> <out_prefix>")
pre <- args[2]

x <- readRDS(args[1])

if (!methods::is(x, "sparseMatrix")) {
  if (is.matrix(x) || is.data.frame(x)) quit(status = 3L)
  stop("RDS object is neither a sparse matrix nor a dense matrix / data.frame")
}

x <- as(x, "CsparseMatrix")          # dgCMatrix: @i, @p, @x

if (is.null(rownames(x))) stop("matrix has no row names (need chr:start-end:strand)")
if (is.null(colnames(x))) stop("matrix has no column names (need sample names)")

writeBin(as.integer(x@i), paste0(pre, ".indices"), size = 4L)
writeBin(as.integer(x@p), paste0(pre, ".indptr"), size = 4L)
writeBin(as.double(x@x), paste0(pre, ".data"))
writeLines(rownames(x), paste0(pre, ".rows"))
writeLines(colnames(x), paste0(pre, ".cols"))
writeLines(paste(nrow(x), ncol(x)), paste0(pre, ".shape"))
