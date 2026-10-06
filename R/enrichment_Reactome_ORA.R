# ============================================================
# GSE213001 — W6.11a Reactome ORA
#
# Direction:
#   Higher expression in IPF relative to NDC
#
# Frozen enrichment settings:
#   resource: Reactome
#   foreground: frozen higher-in-IPF ORA set
#   background: frozen enrichment-eligible universe
#   identifier: ENTREZID
#   minGSSize: 10
#   maxGSSize: 500
#   adjustment: Benjamini-Hochberg
#   significance: adjusted P < 0.05
#
# This script is specified before inspection of Reactome results.
# ============================================================

source("renv/activate.R")

stopifnot(
  requireNamespace("ReactomePA", quietly = TRUE)
)

background_file <-
  "results/enrichment/input/W6_enrichment_background.csv"

foreground_file <-
  "results/enrichment/input/W6_ORA_higher_in_IPF.csv"

output_dir <- "results/enrichment/Reactome"

all_file <- file.path(
  output_dir,
  "W6_Reactome_ORA_higher_in_IPF_all_terms.csv"
)

sig_file <- file.path(
  output_dir,
  "W6_Reactome_ORA_higher_in_IPF_significant.csv"
)

summary_file <-
  "metadata/qc_W6_Reactome_ORA_higher_in_IPF.csv"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Frozen input
# ------------------------------------------------------------

background <- read.csv(
  background_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

foreground <- read.csv(
  foreground_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  "ENTREZID" %in% names(background),
  "ENTREZID" %in% names(foreground)
)

background_entrez <- trimws(as.character(background$ENTREZID))
foreground_entrez <- trimws(as.character(foreground$ENTREZID))

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
  length(unique(background_entrez)) == 14827,
  length(foreground_entrez) == 1372,
  length(unique(foreground_entrez)) == 1372,
  all(foreground_entrez %in% background_entrez)
)

# ------------------------------------------------------------
# 2. Frozen settings certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.11a — REACTOME ORA HIGHER IN IPF\n")
cat("============================================================\n")

cat("\n=== FROZEN INPUT ===\n")
cat("Background genes:       ", length(background_entrez), "\n", sep = "")
cat("Higher-in-IPF genes:    ", length(foreground_entrez), "\n", sep = "")
cat("Identifier:              ENTREZID\n")
cat(
  "Foreground in background: ",
  all(foreground_entrez %in% background_entrez),
  "\n",
  sep = ""
)

cat("\n=== FROZEN SETTINGS ===\n")
cat("Resource:                Reactome\n")
cat("Organism:                human\n")
cat("minGSSize:               10\n")
cat("maxGSSize:               500\n")
cat("Adjustment:              BH\n")
cat("Significance criterion:  adjusted P < 0.05\n")
cat("Result retention:        permissive cutoffs; filter only for reporting\n")

# ------------------------------------------------------------
# 3. Reactome ORA
# ------------------------------------------------------------

reactome <- ReactomePA::enrichPathway(
  gene = foreground_entrez,
  organism = "human",
  universe = background_entrez,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)

reactome_all <- as.data.frame(reactome)

# ------------------------------------------------------------
# 4. Result integrity
# ------------------------------------------------------------

if (nrow(reactome_all) > 0) {

  stopifnot(
    "ID" %in% names(reactome_all),
    "pvalue" %in% names(reactome_all),
    "p.adjust" %in% names(reactome_all),
    !anyDuplicated(reactome_all$ID),
    all(
      reactome_all$pvalue >= 0 &
        reactome_all$pvalue <= 1,
      na.rm = TRUE
    ),
    all(
      reactome_all$p.adjust >= 0 &
        reactome_all$p.adjust <= 1,
      na.rm = TRUE
    )
  )

  reactome_all$BH_significant <-
    reactome_all$p.adjust < 0.05

  reactome_sig <- reactome_all[
    reactome_all$BH_significant %in% TRUE,
    ,
    drop = FALSE
  ]

} else {

  reactome_all$BH_significant <- logical(0)
  reactome_sig <- reactome_all
}

# ------------------------------------------------------------
# 5. Outputs
# ------------------------------------------------------------

write.csv(
  reactome_all,
  all_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  reactome_sig,
  sig_file,
  row.names = FALSE,
  na = ""
)

minimum_raw_p <-
  if (nrow(reactome_all) > 0) {
    min(reactome_all$pvalue, na.rm = TRUE)
  } else {
    NA_real_
  }

minimum_adjusted_p <-
  if (nrow(reactome_all) > 0) {
    min(reactome_all$p.adjust, na.rm = TRUE)
  } else {
    NA_real_
  }

summary_df <- data.frame(
  metric = c(
    "background_genes",
    "foreground_genes",
    "reactome_terms_returned",
    "BH_adjusted_P_lt_0.05",
    "minimum_raw_P",
    "minimum_adjusted_P",
    "minGSSize",
    "maxGSSize"
  ),
  value = c(
    length(background_entrez),
    length(foreground_entrez),
    nrow(reactome_all),
    nrow(reactome_sig),
    minimum_raw_p,
    minimum_adjusted_p,
    10,
    500
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
# 6. Console certificate
# ------------------------------------------------------------

cat("\n=== REACTOME ORA RESULT SUMMARY ===\n")
cat("Reactome terms returned: ", nrow(reactome_all), "\n", sep = "")
cat("BH-adjusted P < 0.05:    ", nrow(reactome_sig), "\n", sep = "")

if (nrow(reactome_all) > 0) {

  cat(
    "Minimum raw P:          ",
    format(
      minimum_raw_p,
      scientific = TRUE,
      digits = 5
    ),
    "\n",
    sep = ""
  )

  cat(
    "Minimum adjusted P:     ",
    format(
      minimum_adjusted_p,
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
cat("No Reactome pathway names were used to alter analysis settings.\n")
cat("No biological interpretation is performed by this script.\n")

cat("\n============================================================\n")
cat("W6.11a REACTOME ORA HIGHER IN IPF: COMPLETED\n")
cat("============================================================\n")
