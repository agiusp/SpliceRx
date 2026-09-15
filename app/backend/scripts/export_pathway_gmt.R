#!/usr/bin/env Rscript
# Dump xCell's cell-type signatures and ConsensusTME's per-cohort tumour-
# microenvironment gene sets to GMT files, bundled alongside the existing
# MSigDB Hallmark library at backend/data/gmt/ (see sjvc/services/pathways.py
# LIBRARIES). Both packages are R-only (xCell: dviraran/xCell on GitHub,
# ConsensusTME: cansysbio/ConsensusTME on GitHub) with no Python equivalent,
# so - like Hallmark - they're exported once here rather than fetched live by
# the Python backend.
#
# GMT format (one gene set per line, tab-separated):
#   <term>\t<description>\t<gene>\t<gene>\t...
#
# Re-run whenever xCell / ConsensusTME are updated:
#   Rscript backend/scripts/export_pathway_gmt.R

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grep("^--file=", args0)])
here <- if (length(file_arg)) dirname(normalizePath(file_arg)) else "backend/scripts"
gmt_dir <- file.path(dirname(here), "data", "gmt")
dir.create(gmt_dir, recursive = TRUE, showWarnings = FALSE)

write_gmt <- function(terms, path) {
  con <- file(path, "w")
  on.exit(close(con))
  for (term in names(terms)) {
    genes <- unique(terms[[term]])
    genes <- genes[nzchar(genes)]
    if (!length(genes)) next
    writeLines(paste(c(term, "", genes), collapse = "\t"), con)
  }
}

# --- xCell: one gene set per cell-type signature ---------------------------
sigs <- xCell::xCell.data$signatures
xcell_terms <- lapply(sigs, GSEABase::geneIds)
names(xcell_terms) <- vapply(sigs, GSEABase::setName, character(1))
write_gmt(xcell_terms, file.path(gmt_dir, "xCell.gmt"))
cat(sprintf("xCell.gmt: %d gene sets\n", length(xcell_terms)))

# --- ConsensusTME: every cohort's cell-type gene sets, term = "<COHORT>_<cell type>" ---
cons <- ConsensusTME::consensusGeneSets
cons_terms <- list()
for (cohort in names(cons)) {
  for (cell_type in names(cons[[cohort]])) {
    cons_terms[[paste0(cohort, "_", cell_type)]] <- cons[[cohort]][[cell_type]]
  }
}
write_gmt(cons_terms, file.path(gmt_dir, "ConsensusTME.gmt"))
cat(sprintf("ConsensusTME.gmt: %d gene sets across %d cohorts\n",
            length(cons_terms), length(cons)))
