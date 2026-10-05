# ============================================================
# GSE213001 — Gene annotation
#
# W6 annotation layer.
#
# Frozen input:
#   15,012 Ensembl gene IDs from the completed W5 primary DE.
#
# Annotation source:
#   org.Hs.eg.db
#
# Primary identifier:
#   Ensembl gene ID
#
# IMPORTANT:
#   This file currently declares the W6 annotation dependencies
#   only. No gene mapping or biological interpretation is
#   performed at this stage.
# ============================================================

source("renv/activate.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

stopifnot(
  as.character(packageVersion("AnnotationDbi")) == "1.70.0",
  as.character(packageVersion("org.Hs.eg.db")) == "3.21.0"
)

cat("\n========================================\n")
cat("GSE213001 W6.2 — ANNOTATION ENVIRONMENT\n")
cat("========================================\n\n")

cat(
  "AnnotationDbi: ",
  as.character(packageVersion("AnnotationDbi")),
  "\n",
  sep = ""
)

cat(
  "org.Hs.eg.db:  ",
  as.character(packageVersion("org.Hs.eg.db")),
  "\n",
  sep = ""
)

cat("\nNo gene mapping performed.\n")
cat("No biological interpretation performed.\n")

cat("\n========================================\n")
cat("W6.2 ANNOTATION ENVIRONMENT: PASSED\n")
cat("========================================\n")
