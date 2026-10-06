# ============================================================
# GSE213001 — W6.12c GO Biological Process ranked enrichment
# ============================================================
#
# Frozen ranked-enrichment settings:
#
#   method: preranked GSEA
#   implementation: clusterProfiler::gseGO()
#   engine: fgsea
#   ontology: BP
#   identifier: ENTREZID
#   ranking statistic: frozen W5 moderated t
#   exponent: 1
#   minGSSize: 10
#   maxGSSize: 500
#   pAdjustMethod: BH
#   pvalueCutoff: 1
#   eps: 0
#   seed: 20261006
#
# IMPORTANT:
#   - no W5 statistics are recalculated
#   - no genes are selected using pathway results
#   - no pathway names are used to alter settings
#   - the complete returned result table is retained
# ============================================================

source("renv/activate.R")

stopifnot(
  requireNamespace("clusterProfiler", quietly = TRUE),
  requireNamespace("org.Hs.eg.db", quietly = TRUE),
  requireNamespace("fgsea", quietly = TRUE)
)

input_file <-
  "results/enrichment/input/W6_ranked_enrichment_moderated_t.csv"

output_dir <-
  "results/enrichment/GO_BP"

all_file <- file.path(
  output_dir,
  "W6_GO_BP_GSEA_all_terms.csv"
)

sig_file <- file.path(
  output_dir,
  "W6_GO_BP_GSEA_significant.csv"
)

summary_file <-
  "metadata/qc_W6_GO_BP_GSEA.csv"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Read frozen ranked input
# ------------------------------------------------------------

stopifnot(file.exists(input_file))

x <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "gene_id",
  "ENTREZID",
  "t"
)

stopifnot(
  all(required_cols %in% names(x)),
  nrow(x) == 14827L,
  length(unique(x$gene_id)) == 14827L,
  length(unique(x$ENTREZID)) == 14827L,
  sum(is.na(x$ENTREZID) | x$ENTREZID == "") == 0L,
  sum(is.na(x$t)) == 0L,
  sum(!is.finite(x$t)) == 0L
)

# ------------------------------------------------------------
# 2. Reconstruct frozen ranked vector
# ------------------------------------------------------------

gene_list <- as.numeric(x$t)

names(gene_list) <- as.character(x$ENTREZID)

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

stopifnot(
  length(gene_list) == 14827L,
  length(unique(names(gene_list))) == 14827L,
  sum(gene_list > 0) == 7334L,
  sum(gene_list < 0) == 7493L,
  sum(gene_list == 0) == 0L,
  anyDuplicated(gene_list) == 0L,
  all(diff(gene_list) <= 0),
  isTRUE(
    all.equal(
      max(gene_list),
      15.61064,
      tolerance = 1e-5
    )
  ),
  isTRUE(
    all.equal(
      min(gene_list),
      -13.26736,
      tolerance = 1e-5
    )
  )
)

# ------------------------------------------------------------
# 3. Frozen execution certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.12c — GO:BP RANKED ENRICHMENT\n")
cat("============================================================\n")

cat("\n=== FROZEN RANKED INPUT ===\n")
cat("Genes:                    ", length(gene_list), "\n", sep = "")
cat("Unique ENTREZID:          ", length(unique(names(gene_list))), "\n", sep = "")
cat("Positive moderated t:     ", sum(gene_list > 0), "\n", sep = "")
cat("Negative moderated t:     ", sum(gene_list < 0), "\n", sep = "")
cat("Zero moderated t:         ", sum(gene_list == 0), "\n", sep = "")
cat("Duplicated t values:      ", anyDuplicated(gene_list), "\n", sep = "")
cat("Maximum moderated t:      ", max(gene_list), "\n", sep = "")
cat("Minimum moderated t:      ", min(gene_list), "\n", sep = "")

cat("\n=== FROZEN SETTINGS ===\n")
cat("Method:                   preranked GSEA\n")
cat("Implementation:           clusterProfiler::gseGO\n")
cat("Engine:                   fgsea\n")
cat("Ontology:                 BP\n")
cat("Identifier:               ENTREZID\n")
cat("Exponent:                 1\n")
cat("minGSSize:                10\n")
cat("maxGSSize:                500\n")
cat("Adjustment:               BH\n")
cat("pvalueCutoff:             1\n")
cat("eps:                      0\n")
cat("Seed:                     20261006\n")

# ------------------------------------------------------------
# 4. GO Biological Process ranked enrichment
# ------------------------------------------------------------

set.seed(20261006)

go_gsea <- clusterProfiler::gseGO(
  geneList = gene_list,
  ont = "BP",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENTREZID",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE,
  by = "fgsea"
)

go_all <- as.data.frame(go_gsea)

stopifnot(
  nrow(go_all) > 0L,
  all(
    c(
      "ID",
      "Description",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust"
    ) %in% names(go_all)
  ),
  length(unique(go_all$ID)) == nrow(go_all),
  all(
    go_all$pvalue[!is.na(go_all$pvalue)] >= 0 &
      go_all$pvalue[!is.na(go_all$pvalue)] <= 1
  ),
  all(
    go_all$p.adjust[!is.na(go_all$p.adjust)] >= 0 &
      go_all$p.adjust[!is.na(go_all$p.adjust)] <= 1
  )
)

# ------------------------------------------------------------
# 5. Direction and reporting significance
# ------------------------------------------------------------

go_all$direction <- ifelse(
  is.na(go_all$NES),
  NA_character_,
  ifelse(
    go_all$NES > 0,
    "higher_in_IPF",
    ifelse(
      go_all$NES < 0,
      "lower_in_IPF",
      "neutral"
    )
  )
)

go_all$BH_significant <-
  !is.na(go_all$p.adjust) &
  go_all$p.adjust < 0.05

go_sig <- go_all[
  go_all$BH_significant,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 6. Write complete and significant results
# ------------------------------------------------------------

write.csv(
  go_all,
  all_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  go_sig,
  sig_file,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 7. Numerical execution audit
# ------------------------------------------------------------

summary_df <- data.frame(
  metric = c(
    "ranked_input_genes",
    "positive_moderated_t",
    "negative_moderated_t",
    "zero_moderated_t",
    "GO_BP_terms_returned",
    "BH_adjusted_P_lt_0.05",
    "positive_NES_terms",
    "negative_NES_terms",
    "missing_NES",
    "missing_raw_P",
    "missing_adjusted_P",
    "minimum_raw_P",
    "minimum_adjusted_P",
    "minGSSize",
    "maxGSSize",
    "seed"
  ),
  value = c(
    length(gene_list),
    sum(gene_list > 0),
    sum(gene_list < 0),
    sum(gene_list == 0),
    nrow(go_all),
    nrow(go_sig),
    sum(go_all$NES > 0, na.rm = TRUE),
    sum(go_all$NES < 0, na.rm = TRUE),
    sum(is.na(go_all$NES)),
    sum(is.na(go_all$pvalue)),
    sum(is.na(go_all$p.adjust)),
    if (
      all(is.na(go_all$pvalue))
    ) NA_real_ else min(go_all$pvalue, na.rm = TRUE),
    if (
      all(is.na(go_all$p.adjust))
    ) NA_real_ else min(go_all$p.adjust, na.rm = TRUE),
    10,
    500,
    20261006
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  summary_file,
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 8. Final console certificate
# ------------------------------------------------------------

cat("\n=== GO:BP RANKED ENRICHMENT RESULT SUMMARY ===\n")
cat("GO:BP terms returned:      ", nrow(go_all), "\n", sep = "")
cat("BH-adjusted P < 0.05:      ", nrow(go_sig), "\n", sep = "")
cat(
  "Positive NES terms:       ",
  sum(go_all$NES > 0, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Negative NES terms:       ",
  sum(go_all$NES < 0, na.rm = TRUE),
  "\n",
  sep = ""
)
cat("Missing NES:               ", sum(is.na(go_all$NES)), "\n", sep = "")

if (any(!is.na(go_all$pvalue))) {
  cat(
    "Minimum raw P:           ",
    format(
      min(go_all$pvalue, na.rm = TRUE),
      scientific = TRUE,
      digits = 5
    ),
    "\n",
    sep = ""
  )
}

if (any(!is.na(go_all$p.adjust))) {
  cat(
    "Minimum adjusted P:      ",
    format(
      min(go_all$p.adjust, na.rm = TRUE),
      scientific = TRUE,
      digits = 5
    ),
    "\n",
    sep = ""
  )
}

cat("\n=== OUTPUTS ===\n")
cat(all_file, "\n")
cat(sig_file, "\n")
cat(summary_file, "\n")

cat("\nIMPORTANT:\n")
cat("No W5 differential-expression statistics were recalculated.\n")
cat("The frozen 15,012-gene DE universe was not modified.\n")
cat("Ranked enrichment used 14,827 one-to-one ENTREZID genes.\n")
cat("No pathway names were used to alter analysis settings.\n")
cat("No biological interpretation is performed by this script.\n")

cat("\n============================================================\n")
cat("W6.12c GO:BP RANKED ENRICHMENT: COMPLETED\n")
cat("============================================================\n")
