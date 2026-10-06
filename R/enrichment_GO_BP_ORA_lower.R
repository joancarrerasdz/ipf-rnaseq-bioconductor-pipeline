# ============================================================
# GSE213001 — W6.10c GO Biological Process ORA
#
# Direction:
#   Lower expression in IPF relative to NDC
#
# Frozen enrichment settings:
#   ontology: Biological Process
#   foreground: frozen lower-in-IPF ORA set
#   background: frozen enrichment-eligible universe
#   identifier: ENTREZID
#   minGSSize: 10
#   maxGSSize: 500
#   adjustment: Benjamini-Hochberg
#   significance: adjusted P < 0.05
#
# IMPORTANT:
#   W5 differential-expression statistics are not recalculated.
#   No pathway results are used to modify the statistical model,
#   gene universe, gene ranking, or enrichment settings.
# ============================================================

source("renv/activate.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

background_file <-
  "results/enrichment/input/W6_enrichment_background.csv"

foreground_file <-
  "results/enrichment/input/W6_ORA_lower_in_IPF.csv"

output_dir <- "results/enrichment/GO_BP"

all_file <- file.path(
  output_dir,
  "W6_GO_BP_ORA_lower_in_IPF_all_terms.csv"
)

sig_file <- file.path(
  output_dir,
  "W6_GO_BP_ORA_lower_in_IPF_significant.csv"
)

summary_file <-
  "metadata/qc_W6_GO_BP_ORA_lower_in_IPF.csv"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Read frozen inputs
# ------------------------------------------------------------

background <- read.csv(
  background_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

foreground <- read.csv(
  foreground_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

required_columns <- c(
  "gene_id",
  "ENTREZID",
  "SYMBOL",
  "logFC",
  "t",
  "adj.P.Val"
)

stopifnot(
  all(required_columns %in% names(background)),
  all(required_columns %in% names(foreground))
)

# ------------------------------------------------------------
# 3. Freeze identifier vectors
# ------------------------------------------------------------

background_entrez <- trimws(background$ENTREZID)
foreground_entrez <- trimws(foreground$ENTREZID)

background_entrez <- background_entrez[
  !is.na(background_entrez) &
    background_entrez != ""
]

foreground_entrez <- foreground_entrez[
  !is.na(foreground_entrez) &
    foreground_entrez != ""
]

stopifnot(
  length(background_entrez) == 14827,
  length(foreground_entrez) == 465,
  length(unique(background_entrez)) == 14827,
  length(unique(foreground_entrez)) == 465,
  all(foreground_entrez %in% background_entrez)
)

# Confirm compatibility with frozen OrgDb.

orgdb_entrez <- keys(
  org.Hs.eg.db,
  keytype = "ENTREZID"
)

background_not_in_orgdb <-
  setdiff(background_entrez, orgdb_entrez)

foreground_not_in_orgdb <-
  setdiff(foreground_entrez, orgdb_entrez)

stopifnot(
  length(background_not_in_orgdb) == 0,
  length(foreground_not_in_orgdb) == 0
)

# ------------------------------------------------------------
# 4. Input certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.10c — GO:BP ORA LOWER IN IPF\n")
cat("============================================================\n")

cat("\n=== FROZEN INPUT ===\n")
cat("Background genes:        ", length(background_entrez), "\n", sep = "")
cat("Lower-in-IPF genes:     ", length(foreground_entrez), "\n", sep = "")
cat("Identifier:               ENTREZID\n")
cat("Foreground in background:", all(foreground_entrez %in% background_entrez), "\n")
cat("Background OrgDb-valid:  ", length(background_not_in_orgdb) == 0, "\n")
cat("Foreground OrgDb-valid:  ", length(foreground_not_in_orgdb) == 0, "\n")

cat("\n=== FROZEN SETTINGS ===\n")
cat("Ontology:                 BP\n")
cat("minGSSize:                10\n")
cat("maxGSSize:                500\n")
cat("Adjustment:               BH\n")
cat("Significance criterion:   adjusted P < 0.05\n")
cat("Result retention:         permissive cutoffs; filter only for reporting\n")

# ------------------------------------------------------------
# 5. GO Biological Process ORA
#
# pvalueCutoff = 1 and qvalueCutoff = 1 intentionally preserve
# the complete enrichment result table returned by enrichGO()
# under the frozen gene-set-size constraints.
#
# Statistical significance is assigned downstream using the
# frozen BH-adjusted P < 0.05 criterion.
# ------------------------------------------------------------

ego <- enrichGO(
  gene = foreground_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)

go_all <- as.data.frame(ego)

stopifnot(
  nrow(go_all) > 0,
  all(c(
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "geneID",
    "Count"
  ) %in% names(go_all))
)

# ------------------------------------------------------------
# 6. Deterministic reporting fields
# ------------------------------------------------------------

go_all$BH_significant <-
  !is.na(go_all$p.adjust) &
  go_all$p.adjust < 0.05

go_all$direction <- "Higher_in_IPF"
go_all$ontology <- "GO_Biological_Process"

# Deterministic ordering:
# adjusted P, raw P, GO identifier.

go_all <- go_all[
  order(
    go_all$p.adjust,
    go_all$pvalue,
    go_all$ID,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

go_sig <- go_all[
  go_all$BH_significant,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 7. Write outputs
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
# 8. Quantitative audit
# ------------------------------------------------------------

summary_df <- data.frame(
  metric = c(
    "background_genes",
    "foreground_genes",
    "background_unique_ENTREZID",
    "foreground_unique_ENTREZID",
    "background_not_in_orgdb",
    "foreground_not_in_orgdb",
    "GO_BP_terms_returned",
    "GO_BP_BH_significant",
    "minGSSize",
    "maxGSSize"
  ),
  value = c(
    length(background_entrez),
    length(foreground_entrez),
    length(unique(background_entrez)),
    length(unique(foreground_entrez)),
    length(background_not_in_orgdb),
    length(foreground_not_in_orgdb),
    nrow(go_all),
    nrow(go_sig),
    10,
    500
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  summary_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Final console certificate
# ------------------------------------------------------------

cat("\n=== ORA RESULT SUMMARY ===\n")
cat("GO:BP terms returned:    ", nrow(go_all), "\n", sep = "")
cat("BH-adjusted P < 0.05:    ", nrow(go_sig), "\n", sep = "")

if (nrow(go_all) > 0) {
  cat(
    "Minimum raw P:         ",
    format(
      min(go_all$pvalue, na.rm = TRUE),
      scientific = TRUE,
      digits = 5
    ),
    "\n",
    sep = ""
  )

  cat(
    "Minimum adjusted P:    ",
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
cat("The enrichment universe contains 14,827 one-to-one ENTREZID genes.\n")
cat("No pathway names were used to alter analysis settings.\n")
cat("No biological interpretation is performed by this script.\n")

cat("\n============================================================\n")
cat("W6.10c GO:BP ORA LOWER IN IPF: COMPLETED\n")
cat("============================================================\n")
