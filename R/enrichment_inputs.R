# ============================================================
# GSE213001 — W6 functional-enrichment input construction
#
# IMPORTANT:
#   This script only constructs and audits enrichment inputs.
#   No pathway enrichment is performed here.
#
# Frozen source:
#   W5 donor-aware primary DE + W6 annotation.
# ============================================================

source("renv/activate.R")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- paste0(
  "results/differential_expression/",
  "primary_IPF_vs_NDC_annotated.csv"
)

output_dir <- "results/enrichment/input"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read frozen annotated DE table
# ------------------------------------------------------------

de <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "gene_id",
  "logFC",
  "t",
  "adj.P.Val",
  "SYMBOL",
  "ENTREZID",
  "mapping_status"
)

stopifnot(all(required_cols %in% names(de)))

de$ENTREZID <- as.character(de$ENTREZID)

cat("\n============================================================\n")
cat("GSE213001 W6.8b — ENRICHMENT INPUT CONSTRUCTION\n")
cat("============================================================\n")

cat("\n=== FROZEN DE INPUT ===\n")
cat("Rows:               ", nrow(de), "\n", sep = "")
cat(
  "Unique gene IDs:    ",
  length(unique(de$gene_id)),
  "\n",
  sep = ""
)

stopifnot(
  nrow(de) == 15012,
  length(unique(de$gene_id)) == 15012
)

# ------------------------------------------------------------
# 3. Enrichment-eligible universe
# ------------------------------------------------------------

eligible <- (
  de$mapping_status == "single_mapping_record" &
  !is.na(de$ENTREZID) &
  nzchar(trimws(de$ENTREZID))
)

background <- de[eligible, , drop = FALSE]

cat("\n=== ENRICHMENT-ELIGIBLE BACKGROUND ===\n")
cat("Eligible genes:      ", nrow(background), "\n", sep = "")
cat(
  "Unique ENTREZIDs:    ",
  length(unique(background$ENTREZID)),
  "\n",
  sep = ""
)

stopifnot(
  nrow(background) == 14827,
  length(unique(background$ENTREZID)) == 14827
)

# ------------------------------------------------------------
# 4. Frozen directional ORA sets
# ------------------------------------------------------------

higher <- background[
  background$adj.P.Val < 0.05 &
  background$logFC >= 1,
  ,
  drop = FALSE
]

lower <- background[
  background$adj.P.Val < 0.05 &
  background$logFC <= -1,
  ,
  drop = FALSE
]

cat("\n=== FROZEN ORA FOREGROUNDS ===\n")
cat("Higher in IPF:       ", nrow(higher), "\n", sep = "")
cat("Lower in IPF:        ", nrow(lower), "\n", sep = "")

stopifnot(
  nrow(higher) == 1372,
  nrow(lower) == 465
)

stopifnot(
  length(intersect(higher$gene_id, lower$gene_id)) == 0
)

# ------------------------------------------------------------
# 5. Ranked enrichment input
# ------------------------------------------------------------

ranked <- background[
  order(background$t, decreasing = TRUE),
  c(
    "gene_id",
    "ENTREZID",
    "SYMBOL",
    "t",
    "logFC",
    "adj.P.Val"
  ),
  drop = FALSE
]

stopifnot(
  nrow(ranked) == 14827,
  !anyNA(ranked$t),
  length(unique(ranked$ENTREZID)) == 14827
)

cat("\n=== RANKED ENRICHMENT INPUT ===\n")
cat("Ranked genes:        ", nrow(ranked), "\n", sep = "")
cat(
  "Maximum moderated t: ",
  signif(max(ranked$t), 7),
  "\n",
  sep = ""
)
cat(
  "Minimum moderated t: ",
  signif(min(ranked$t), 7),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Write deterministic inputs
# ------------------------------------------------------------

background_out <- background[
  ,
  c(
    "gene_id",
    "ENTREZID",
    "SYMBOL",
    "logFC",
    "t",
    "adj.P.Val"
  ),
  drop = FALSE
]

higher_out <- higher[
  ,
  c(
    "gene_id",
    "ENTREZID",
    "SYMBOL",
    "logFC",
    "t",
    "adj.P.Val"
  ),
  drop = FALSE
]

lower_out <- lower[
  ,
  c(
    "gene_id",
    "ENTREZID",
    "SYMBOL",
    "logFC",
    "t",
    "adj.P.Val"
  ),
  drop = FALSE
]

write.csv(
  background_out,
  file.path(output_dir, "W6_enrichment_background.csv"),
  row.names = FALSE
)

write.csv(
  higher_out,
  file.path(output_dir, "W6_ORA_higher_in_IPF.csv"),
  row.names = FALSE
)

write.csv(
  lower_out,
  file.path(output_dir, "W6_ORA_lower_in_IPF.csv"),
  row.names = FALSE
)

write.csv(
  ranked,
  file.path(output_dir, "W6_ranked_enrichment_moderated_t.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. QC certificate
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "frozen_DE_genes",
    "enrichment_background",
    "unique_background_ENTREZID",
    "ORA_higher_IPF",
    "ORA_lower_IPF",
    "ranked_enrichment_genes",
    "higher_lower_overlap",
    "missing_rank_statistics"
  ),
  value = c(
    nrow(de),
    nrow(background),
    length(unique(background$ENTREZID)),
    nrow(higher),
    nrow(lower),
    nrow(ranked),
    length(intersect(higher$gene_id, lower$gene_id)),
    sum(is.na(ranked$t))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  "metadata/qc_W6_enrichment_inputs.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Final certificate
# ------------------------------------------------------------

cat("\n=== OUTPUTS ===\n")
cat("results/enrichment/input/W6_enrichment_background.csv\n")
cat("results/enrichment/input/W6_ORA_higher_in_IPF.csv\n")
cat("results/enrichment/input/W6_ORA_lower_in_IPF.csv\n")
cat("results/enrichment/input/W6_ranked_enrichment_moderated_t.csv\n")
cat("metadata/qc_W6_enrichment_inputs.csv\n")

cat("\nIMPORTANT:\n")
cat("No enrichment analysis was performed.\n")
cat("No pathway names were inspected.\n")
cat("The frozen W5 statistical results were not recalculated.\n")
cat("The frozen 15,012-gene DE universe was not modified.\n")

cat("\n============================================================\n")
cat("W6.8b ENRICHMENT INPUT CONSTRUCTION: PASSED\n")
cat("============================================================\n")
