# ============================================================
# GSE213001 - W6 enrichment environment preflight
#
# Purpose:
#   Declare and verify the enrichment packages required by W6.
#
# IMPORTANT:
#   This script performs no enrichment analysis and inspects
#   no pathway-level biological results.
# ============================================================

source("renv/activate.R")

required_packages <- c(
  "clusterProfiler",
  "ReactomePA",
  "enrichplot",
  "DOSE",
  "fgsea",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

# Explicit declarations are intentional so renv can detect
# these packages as project dependencies.
stopifnot(
  requireNamespace("clusterProfiler", quietly = TRUE),
  requireNamespace("ReactomePA", quietly = TRUE),
  requireNamespace("enrichplot", quietly = TRUE),
  requireNamespace("DOSE", quietly = TRUE),
  requireNamespace("fgsea", quietly = TRUE),
  requireNamespace("AnnotationDbi", quietly = TRUE),
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
)

versions <- vapply(
  required_packages,
  function(pkg) as.character(packageVersion(pkg)),
  character(1)
)

cat("\n============================================================\n")
cat("GSE213001 W6.9 - ENRICHMENT ENVIRONMENT PREFLIGHT\n")
cat("============================================================\n\n")

for (pkg in required_packages) {
  cat(
    sprintf(
      "%-18s : %s\n",
      pkg,
      versions[[pkg]]
    )
  )
}

cat("\nR:           ", R.version.string, "\n", sep = "")
cat(
  "Bioconductor: ",
  as.character(BiocManager::version()),
  "\n",
  sep = ""
)

cat("\nNo enrichment analysis performed.\n")
cat("No pathway-level results inspected.\n")

cat("\n============================================================\n")
cat("W6.9 ENRICHMENT ENVIRONMENT PREFLIGHT: PASSED\n")
cat("============================================================\n")
