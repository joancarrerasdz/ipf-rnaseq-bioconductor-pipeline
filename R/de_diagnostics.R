# ============================================================
# GSE213001 — W5.6 differential-expression diagnostics
#
# Diagnostics only.
# No model selection or gene-level biological interpretation.
# ============================================================

source("renv/activate.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

cat("\n============================================================\n")
cat("GSE213001 W5.6 — DE DIAGNOSTICS\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------

primary <- read.csv(
  "results/differential_expression/primary_IPF_vs_NDC_all_genes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

age_concordance <- read.csv(
  "results/differential_expression/sensitivity_age_gene_concordance.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(primary) == 15012,
  nrow(age_concordance) == 15012,
  !anyDuplicated(primary$gene_id),
  !anyDuplicated(age_concordance$gene_id),
  all(is.finite(primary$logFC)),
  all(is.finite(primary$P.Value)),
  all(is.finite(primary$adj.P.Val))
)

dir.create(
  "figures/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Reporting classifications
#
# FDR < 0.05 is the frozen inferential criterion.
# logFC sign is used only to describe direction.
# ------------------------------------------------------------

primary$DE_status <- "Not FDR-significant"

primary$DE_status[
  primary$adj.P.Val < 0.05 &
    primary$logFC > 0
] <- "Higher in IPF"

primary$DE_status[
  primary$adj.P.Val < 0.05 &
    primary$logFC < 0
] <- "Lower in IPF"

primary$DE_status <- factor(
  primary$DE_status,
  levels = c(
    "Not FDR-significant",
    "Higher in IPF",
    "Lower in IPF"
  )
)

# Avoid Inf values for visualization only.

primary$minus_log10_fdr <-
  -log10(
    pmax(
      primary$adj.P.Val,
      .Machine$double.xmin
    )
  )

# ------------------------------------------------------------
# 3. P-value histogram
# ------------------------------------------------------------

p_pvalue <- ggplot(
  primary,
  aes(x = P.Value)
) +
  geom_histogram(
    bins = 50,
    boundary = 0
  ) +
  labs(
    title = "Primary IPF vs NDC differential expression",
    subtitle = "Raw p-value distribution across 15,012 tested genes",
    x = "Raw p-value",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "figures/differential_expression/primary_pvalue_histogram.png",
  p_pvalue,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 4. Volcano plot
# ------------------------------------------------------------

p_volcano <- ggplot(
  primary,
  aes(
    x = logFC,
    y = minus_log10_fdr,
    shape = DE_status
  )
) +
  geom_point(
    alpha = 0.35,
    size = 1.1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Primary donor-aware differential expression",
    subtitle = "IPF vs NDC; adjusted for lung region and donor correlation",
    x = "log2 fold change (IPF relative to NDC)",
    y = expression(-log[10]("BH FDR")),
    shape = "FDR status"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "figures/differential_expression/primary_volcano.png",
  p_volcano,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 5. MA plot
# ------------------------------------------------------------

p_ma <- ggplot(
  primary,
  aes(
    x = AveExpr,
    y = logFC,
    shape = DE_status
  )
) +
  geom_point(
    alpha = 0.30,
    size = 1.0
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "MA plot — primary IPF vs NDC model",
    subtitle = "FDR < 0.05 defines statistical significance",
    x = "Average log-expression",
    y = "log2 fold change (IPF relative to NDC)",
    shape = "FDR status"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "figures/differential_expression/primary_MA_plot.png",
  p_ma,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 6. Age-adjustment concordance
# ------------------------------------------------------------

age_concordance$robustness_status <-
  "Other"

age_concordance$robustness_status[
  age_concordance$significant_reduced &
    age_concordance$significant_age_adjusted
] <- "FDR < 0.05 in both"

age_concordance$robustness_status[
  age_concordance$significant_reduced &
    !age_concordance$significant_age_adjusted
] <- "Reduced only"

age_concordance$robustness_status[
  !age_concordance$significant_reduced &
    age_concordance$significant_age_adjusted
] <- "Age-adjusted only"

age_concordance$robustness_status <- factor(
  age_concordance$robustness_status,
  levels = c(
    "Other",
    "FDR < 0.05 in both",
    "Reduced only",
    "Age-adjusted only"
  )
)

pearson_logfc <- cor(
  age_concordance$logFC_reduced,
  age_concordance$logFC_age_adjusted,
  method = "pearson"
)

spearman_logfc <- cor(
  age_concordance$logFC_reduced,
  age_concordance$logFC_age_adjusted,
  method = "spearman"
)

p_age <- ggplot(
  age_concordance,
  aes(
    x = logFC_reduced,
    y = logFC_age_adjusted,
    shape = robustness_status
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.3
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  coord_equal() +
  labs(
    title = "Effect-size robustness to age adjustment",
    subtitle = sprintf(
      "Same 99-sample cohort; Pearson r = %.3f, Spearman rho = %.3f",
      pearson_logfc,
      spearman_logfc
    ),
    x = "log2 fold change — reduced model",
    y = "log2 fold change — age-adjusted model",
    shape = "Sensitivity status"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "figures/differential_expression/age_adjustment_logFC_concordance.png",
  p_age,
  width = 9,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Diagnostic summary
# ------------------------------------------------------------

diagnostic_summary <- data.frame(
  metric = c(
    "tested_genes",
    "FDR_lt_0.05",
    "FDR_positive_logFC",
    "FDR_negative_logFC",
    "primary_logFC_median",
    "primary_logFC_IQR",
    "pearson_age_concordance",
    "spearman_age_concordance"
  ),
  value = c(
    nrow(primary),
    sum(primary$adj.P.Val < 0.05),
    sum(primary$adj.P.Val < 0.05 & primary$logFC > 0),
    sum(primary$adj.P.Val < 0.05 & primary$logFC < 0),
    median(primary$logFC),
    IQR(primary$logFC),
    pearson_logfc,
    spearman_logfc
  ),
  stringsAsFactors = FALSE
)

write.csv(
  diagnostic_summary,
  "results/differential_expression/W5_DE_diagnostic_summary.csv",
  row.names = FALSE
)

cat("\n=== PRIMARY DE ===\n")
cat("Genes tested:         ", nrow(primary), "\n", sep = "")
cat(
  "FDR < 0.05:           ",
  sum(primary$adj.P.Val < 0.05),
  "\n",
  sep = ""
)
cat(
  "Higher in IPF:        ",
  sum(primary$adj.P.Val < 0.05 & primary$logFC > 0),
  "\n",
  sep = ""
)
cat(
  "Lower in IPF:         ",
  sum(primary$adj.P.Val < 0.05 & primary$logFC < 0),
  "\n",
  sep = ""
)

cat("\n=== AGE ROBUSTNESS ===\n")
cat(
  "Pearson logFC:        ",
  sprintf("%.6f", pearson_logfc),
  "\n",
  sep = ""
)
cat(
  "Spearman logFC:       ",
  sprintf("%.6f", spearman_logfc),
  "\n",
  sep = ""
)

cat("\n=== FIGURES WRITTEN ===\n")
cat("figures/differential_expression/primary_pvalue_histogram.png\n")
cat("figures/differential_expression/primary_volcano.png\n")
cat("figures/differential_expression/primary_MA_plot.png\n")
cat("figures/differential_expression/age_adjustment_logFC_concordance.png\n")

cat("\nNo gene-level biological interpretation was performed.\n")

cat("\n============================================================\n")
cat("W5.6 DE DIAGNOSTICS: COMPLETED\n")
cat("============================================================\n")
