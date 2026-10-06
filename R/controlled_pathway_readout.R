# ============================================================
# GSE213001 — W6.13c Controlled pathway readout
# ============================================================
#
# Purpose:
# Construct the deterministic pathway-level biological readout
# using the reporting order frozen in docs/DECISIONS.md.
#
# IMPORTANT:
# - No pathway names are used for selection.
# - No biological interpretation is performed here.
# - No enrichment analysis is recalculated.
# - No W5 statistics are modified.
# - No semantic redundancy filtering is applied.
#
# Frozen reporting rules:
#
# ORA:
#   1. p.adjust ascending
#   2. pvalue ascending
#   3. Count descending
#   4. ID ascending
#
# Ranked enrichment:
#   split by NES sign
#   1. p.adjust ascending
#   2. abs(NES) descending
#   3. ID ascending
#
# Primary controlled readout:
#   top 20 per resource / direction
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

files <- list(
  GO_BP_ORA_higher = paste0(
    "results/enrichment/GO_BP/",
    "W6_GO_BP_ORA_higher_in_IPF_significant.csv"
  ),
  GO_BP_ORA_lower = paste0(
    "results/enrichment/GO_BP/",
    "W6_GO_BP_ORA_lower_in_IPF_significant.csv"
  ),
  Reactome_ORA_higher = paste0(
    "results/enrichment/Reactome/",
    "W6_Reactome_ORA_higher_in_IPF_significant.csv"
  ),
  Reactome_ORA_lower = paste0(
    "results/enrichment/Reactome/",
    "W6_Reactome_ORA_lower_in_IPF_significant.csv"
  ),
  GO_BP_GSEA = paste0(
    "results/enrichment/GO_BP/",
    "W6_GO_BP_GSEA_significant.csv"
  )
)

stopifnot(all(file.exists(unlist(files))))

# ------------------------------------------------------------
# 2. Read frozen significant-result tables
# ------------------------------------------------------------

go_up <- read.csv(
  files$GO_BP_ORA_higher,
  check.names = FALSE
)

go_down <- read.csv(
  files$GO_BP_ORA_lower,
  check.names = FALSE
)

reactome_up <- read.csv(
  files$Reactome_ORA_higher,
  check.names = FALSE
)

reactome_down <- read.csv(
  files$Reactome_ORA_lower,
  check.names = FALSE
)

gsea <- read.csv(
  files$GO_BP_GSEA,
  check.names = FALSE
)

# ------------------------------------------------------------
# 3. Structural expectations
# ------------------------------------------------------------

ora_required <- c(
  "ID",
  "Description",
  "pvalue",
  "p.adjust",
  "Count"
)

gsea_required <- c(
  "ID",
  "Description",
  "NES",
  "pvalue",
  "p.adjust"
)

stopifnot(
  all(ora_required %in% names(go_up)),
  all(ora_required %in% names(go_down)),
  all(ora_required %in% names(reactome_up)),
  all(ora_required %in% names(reactome_down)),
  all(gsea_required %in% names(gsea))
)

# Frozen significant-result counts already established by QC.
stopifnot(
  nrow(go_up) == 399L,
  nrow(go_down) == 147L,
  nrow(reactome_up) == 65L,
  nrow(reactome_down) == 20L,
  nrow(gsea) == 1054L
)

# No missing ranking quantities.
stopifnot(
  sum(is.na(go_up$p.adjust)) == 0L,
  sum(is.na(go_up$pvalue)) == 0L,
  sum(is.na(go_up$Count)) == 0L,

  sum(is.na(go_down$p.adjust)) == 0L,
  sum(is.na(go_down$pvalue)) == 0L,
  sum(is.na(go_down$Count)) == 0L,

  sum(is.na(reactome_up$p.adjust)) == 0L,
  sum(is.na(reactome_up$pvalue)) == 0L,
  sum(is.na(reactome_up$Count)) == 0L,

  sum(is.na(reactome_down$p.adjust)) == 0L,
  sum(is.na(reactome_down$pvalue)) == 0L,
  sum(is.na(reactome_down$Count)) == 0L,

  sum(is.na(gsea$p.adjust)) == 0L,
  sum(is.na(gsea$NES)) == 0L
)

# Significant tables must contain only BH-adjusted P < 0.05.
stopifnot(
  all(go_up$p.adjust < 0.05),
  all(go_down$p.adjust < 0.05),
  all(reactome_up$p.adjust < 0.05),
  all(reactome_down$p.adjust < 0.05),
  all(gsea$p.adjust < 0.05)
)

# ------------------------------------------------------------
# 4. Deterministic ORA ranking
# ------------------------------------------------------------

rank_ora <- function(x, resource, direction) {

  ord <- order(
    x$p.adjust,
    x$pvalue,
    -x$Count,
    x$ID
  )

  out <- x[ord, , drop = FALSE]

  out$resource <- resource
  out$method <- "ORA"
  out$direction <- direction
  out$reporting_rank <- seq_len(nrow(out))

  out
}

go_up_ranked <- rank_ora(
  go_up,
  resource = "GO:BP",
  direction = "higher_in_IPF"
)

go_down_ranked <- rank_ora(
  go_down,
  resource = "GO:BP",
  direction = "lower_in_IPF"
)

reactome_up_ranked <- rank_ora(
  reactome_up,
  resource = "Reactome",
  direction = "higher_in_IPF"
)

reactome_down_ranked <- rank_ora(
  reactome_down,
  resource = "Reactome",
  direction = "lower_in_IPF"
)

# ------------------------------------------------------------
# 5. Deterministic ranked-enrichment ranking
# ------------------------------------------------------------

gsea_positive <- gsea[
  gsea$NES > 0,
  ,
  drop = FALSE
]

gsea_negative <- gsea[
  gsea$NES < 0,
  ,
  drop = FALSE
]

stopifnot(
  nrow(gsea_positive) >= 20L,
  nrow(gsea_negative) >= 20L
)

rank_gsea <- function(x, direction) {

  ord <- order(
    x$p.adjust,
    -abs(x$NES),
    x$ID
  )

  out <- x[ord, , drop = FALSE]

  out$resource <- "GO:BP"
  out$method <- "ranked_enrichment"
  out$direction <- direction
  out$reporting_rank <- seq_len(nrow(out))

  out
}

gsea_positive_ranked <- rank_gsea(
  gsea_positive,
  direction = "higher_in_IPF"
)

gsea_negative_ranked <- rank_gsea(
  gsea_negative,
  direction = "lower_in_IPF"
)

# ------------------------------------------------------------
# 6. Primary controlled Top-20 readouts
# ------------------------------------------------------------

top_n <- 20L

top_go_up <- head(go_up_ranked, top_n)
top_go_down <- head(go_down_ranked, top_n)

top_reactome_up <- head(
  reactome_up_ranked,
  top_n
)

top_reactome_down <- head(
  reactome_down_ranked,
  top_n
)

top_gsea_positive <- head(
  gsea_positive_ranked,
  top_n
)

top_gsea_negative <- head(
  gsea_negative_ranked,
  top_n
)

stopifnot(
  nrow(top_go_up) == 20L,
  nrow(top_go_down) == 20L,
  nrow(top_reactome_up) == 20L,
  nrow(top_reactome_down) == 20L,
  nrow(top_gsea_positive) == 20L,
  nrow(top_gsea_negative) == 20L
)

# ------------------------------------------------------------
# 7. Output directories
# ------------------------------------------------------------

output_dir <- "results/interpretation/pathways"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 8. Write individual controlled readouts
# ------------------------------------------------------------

write.csv(
  top_go_up,
  file.path(
    output_dir,
    "W6_GO_BP_ORA_higher_in_IPF_top20.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_go_down,
  file.path(
    output_dir,
    "W6_GO_BP_ORA_lower_in_IPF_top20.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_reactome_up,
  file.path(
    output_dir,
    "W6_Reactome_ORA_higher_in_IPF_top20.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_reactome_down,
  file.path(
    output_dir,
    "W6_Reactome_ORA_lower_in_IPF_top20.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_gsea_positive,
  file.path(
    output_dir,
    "W6_GO_BP_GSEA_positive_NES_top20.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_gsea_negative,
  file.path(
    output_dir,
    "W6_GO_BP_GSEA_negative_NES_top20.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Combined audit-friendly readout
# ------------------------------------------------------------

standardize <- function(x) {

  keep <- c(
    "resource",
    "method",
    "direction",
    "reporting_rank",
    "ID",
    "Description",
    "pvalue",
    "p.adjust"
  )

  extra <- intersect(
    c(
      "Count",
      "NES",
      "setSize",
      "GeneRatio",
      "BgRatio",
      "RichFactor",
      "FoldEnrichment",
      "zScore",
      "leading_edge",
      "core_enrichment"
    ),
    names(x)
  )

  x[, c(keep, extra), drop = FALSE]
}

tables <- list(
  standardize(top_go_up),
  standardize(top_go_down),
  standardize(top_reactome_up),
  standardize(top_reactome_down),
  standardize(top_gsea_positive),
  standardize(top_gsea_negative)
)

all_names <- unique(
  unlist(
    lapply(tables, names)
  )
)

tables <- lapply(
  tables,
  function(x) {

    missing_cols <- setdiff(
      all_names,
      names(x)
    )

    for (nm in missing_cols) {
      x[[nm]] <- NA
    }

    x[, all_names, drop = FALSE]
  }
)

combined <- do.call(
  rbind,
  tables
)

rownames(combined) <- NULL

stopifnot(nrow(combined) == 120L)

combined_file <- file.path(
  output_dir,
  "W6_controlled_pathway_readout_top20.csv"
)

write.csv(
  combined,
  combined_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. QC summary
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "GO_BP_ORA_higher_significant",
    "GO_BP_ORA_lower_significant",
    "Reactome_ORA_higher_significant",
    "Reactome_ORA_lower_significant",
    "GO_BP_GSEA_significant",
    "GO_BP_GSEA_positive_NES_significant",
    "GO_BP_GSEA_negative_NES_significant",
    "top20_tables",
    "rows_per_top_table",
    "combined_readout_rows"
  ),
  value = c(
    nrow(go_up),
    nrow(go_down),
    nrow(reactome_up),
    nrow(reactome_down),
    nrow(gsea),
    nrow(gsea_positive),
    nrow(gsea_negative),
    6L,
    top_n,
    nrow(combined)
  ),
  stringsAsFactors = FALSE
)

qc_file <- "metadata/qc_W6_controlled_pathway_readout.csv"

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.13c — CONTROLLED PATHWAY READOUT\n")
cat("============================================================\n")

cat("\n=== FROZEN INPUT COUNTS ===\n")
cat("GO:BP ORA higher:       ", nrow(go_up), "\n")
cat("GO:BP ORA lower:        ", nrow(go_down), "\n")
cat("Reactome ORA higher:    ", nrow(reactome_up), "\n")
cat("Reactome ORA lower:     ", nrow(reactome_down), "\n")
cat("GO:BP ranked enrichment:", nrow(gsea), "\n")

cat("\n=== RANKED ENRICHMENT DIRECTIONS ===\n")
cat("Positive NES terms:     ", nrow(gsea_positive), "\n")
cat("Negative NES terms:     ", nrow(gsea_negative), "\n")

cat("\n=== CONTROLLED READOUT ===\n")
cat("Top tables:             6\n")
cat("Rows per table:         20\n")
cat("Combined rows:          ", nrow(combined), "\n")

cat("\n=== OUTPUT ===\n")
cat(output_dir, "\n")
cat(combined_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No enrichment analysis was recalculated.\n")
cat("No pathway descriptions were used for selection.\n")
cat("No semantic redundancy filtering was applied.\n")
cat("Age robustness did not alter pathway ranking.\n")
cat("No biological interpretation was performed.\n")

cat("\n============================================================\n")
cat("W6.13c CONTROLLED PATHWAY READOUT: COMPLETED\n")
cat("============================================================\n")
