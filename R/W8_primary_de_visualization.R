# ============================================================
# GSE213001
# W8 primary differential-expression portfolio visualization
# ============================================================
#
# PURPOSE
# -------
# Generate a presentation-oriented volcano plot using only the
# analytically frozen W5 IPF-vs-NDC differential-expression result.
#
# IMPORTANT
# ---------
# This script:
# - does NOT fit a differential-expression model,
# - does NOT recompute P values,
# - does NOT recompute adjusted P values,
# - does NOT introduce an effect-size significance threshold,
# - does NOT overwrite any W5-W7 analytical output,
# - and does NOT overwrite the original W5 diagnostic volcano plot.
#
# Frozen inferential criterion:
#   FDR < 0.05
#
# Gene-label selection rule:
#   Among FDR-significant genes with an unambiguous mapped symbol,
#   select the five smallest adjusted P values in each logFC
#   direction. Ties are resolved by larger absolute logFC and then
#   by gene identifier.
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

suppressPackageStartupMessages(
  library(ggplot2)
)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- paste0(
  "results/differential_expression/",
  "primary_IPF_vs_NDC_all_genes.csv"
)

annotation_file <- "metadata/qc_annotation_gene_mapping.csv"

figure_dir <- "figures/final"

output_png <- file.path(
  figure_dir,
  "W8_primary_IPF_vs_NDC_volcano.png"
)

output_pdf <- file.path(
  figure_dir,
  "W8_primary_IPF_vs_NDC_volcano.pdf"
)

label_qc_file <- paste0(
  "metadata/qc_W8_primary_DE_visualization_labels.csv"
)

qc_file <- paste0(
  "metadata/qc_W8_primary_DE_visualization.csv"
)

# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------

select_column <- function(
  column_names,
  preferred,
  fallback_pattern,
  label
) {

  exact <- preferred[
    preferred %in% column_names
  ]

  if (length(exact) >= 1L) {
    return(exact[[1L]])
  }

  hits <- grep(
    fallback_pattern,
    column_names,
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(hits) != 1L) {
    stop(
      paste0(
        "Could not identify exactly one ",
        label,
        " column. Candidates: ",
        paste(hits, collapse = ", ")
      )
    )
  }

  hits[[1L]]
}


normalise_ensembl <- function(x) {

  x <- trimws(as.character(x))

  sub(
    "\\.[0-9]+$",
    "",
    x
  )
}


metric_row <- function(metric, value) {

  data.frame(
    metric = metric,
    value = as.character(value),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 3. Input checks
# ------------------------------------------------------------

stopifnot(
  file.exists(input_file),
  file.exists(annotation_file)
)

x <- read.csv(
  input_file,
  check.names = FALSE
)

cat("\n===== FROZEN W5 INPUT =====\n")
cat("File:", input_file, "\n")
cat("Rows:", nrow(x), "\n")
cat("Columns:", ncol(x), "\n")

gene_col <- select_column(
  names(x),
  preferred = c(
    "gene_id",
    "ensembl_gene_id",
    "ENSEMBL"
  ),
  fallback_pattern = "gene.*id|ensembl",
  label = "gene identifier"
)

logfc_col <- select_column(
  names(x),
  preferred = c(
    "logFC",
    "log2FC",
    "log2FoldChange"
  ),
  fallback_pattern = "^log.*fc$|fold.*change",
  label = "logFC"
)

raw_p_col <- select_column(
  names(x),
  preferred = c(
    "P.Value",
    "pvalue",
    "p_value"
  ),
  fallback_pattern = "^p[._]?value$|^pvalue$",
  label = "raw P value"
)

adj_p_col <- select_column(
  names(x),
  preferred = c(
    "adj.P.Val",
    "padj",
    "FDR"
  ),
  fallback_pattern = "adj.*p|padj|fdr",
  label = "adjusted P value"
)

cat("\nDetected columns:\n")
cat("gene ID:", gene_col, "\n")
cat("logFC:", logfc_col, "\n")
cat("raw P:", raw_p_col, "\n")
cat("adjusted P:", adj_p_col, "\n")

w5 <- data.frame(
  gene_id = as.character(x[[gene_col]]),
  logFC = as.numeric(x[[logfc_col]]),
  P.Value = as.numeric(x[[raw_p_col]]),
  adj.P.Val = as.numeric(x[[adj_p_col]]),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(w5) == 15012L,
  length(unique(w5$gene_id)) == 15012L,
  !anyDuplicated(w5$gene_id),
  all(is.finite(w5$logFC)),
  all(is.finite(w5$P.Value)),
  all(is.finite(w5$adj.P.Val)),
  all(w5$P.Value >= 0 & w5$P.Value <= 1),
  all(w5$adj.P.Val >= 0 & w5$adj.P.Val <= 1)
)

# ------------------------------------------------------------
# 4. Frozen significance classification
# ------------------------------------------------------------

fdr_threshold <- 0.05

w5$FDR_significant <- (
  w5$adj.P.Val < fdr_threshold
)

stopifnot(
  sum(
    w5$FDR_significant &
      w5$logFC == 0
  ) == 0L
)

w5$status <- ifelse(
  !w5$FDR_significant,
  "Not FDR-significant",
  ifelse(
    w5$logFC > 0,
    "Higher in IPF",
    "Lower in IPF"
  )
)

w5$status <- factor(
  w5$status,
  levels = c(
    "Not FDR-significant",
    "Higher in IPF",
    "Lower in IPF"
  )
)

n_genes <- nrow(w5)

n_fdr <- sum(
  w5$FDR_significant
)

n_higher <- sum(
  w5$status == "Higher in IPF"
)

n_lower <- sum(
  w5$status == "Lower in IPF"
)

cat("\n===== FROZEN W5 COUNTS =====\n")
cat("Genes:", n_genes, "\n")
cat("FDR < 0.05:", n_fdr, "\n")
cat("Higher in IPF:", n_higher, "\n")
cat("Lower in IPF:", n_lower, "\n")

# Known frozen-W5 invariants.

stopifnot(
  n_genes == 15012L,
  n_fdr == 7898L,
  n_higher == 3964L,
  n_lower == 3934L
)

# ------------------------------------------------------------
# 5. Plotting significance
# ------------------------------------------------------------

positive_adjusted_p <- w5$adj.P.Val[
  w5$adj.P.Val > 0
]

stopifnot(
  length(positive_adjusted_p) > 0L
)

minimum_positive_adjusted_p <- min(
  positive_adjusted_p
)

zero_adjusted_p <- sum(
  w5$adj.P.Val == 0
)

# Zero values, if present, are placed one order of magnitude beyond
# the smallest positive frozen adjusted P value solely to obtain a
# finite plotting coordinate. The frozen adjusted P values themselves
# are not modified.

w5$plot_adj_P <- ifelse(
  w5$adj.P.Val == 0,
  minimum_positive_adjusted_p / 10,
  w5$adj.P.Val
)

w5$minus_log10_FDR <- -log10(
  w5$plot_adj_P
)

stopifnot(
  all(is.finite(w5$minus_log10_FDR))
)

fdr_line <- -log10(
  fdr_threshold
)

# ------------------------------------------------------------
# 6. Frozen annotation mapping
# ------------------------------------------------------------

annotation <- read.csv(
  annotation_file,
  check.names = FALSE
)

annotation_gene_col <- select_column(
  names(annotation),
  preferred = c(
    "gene_id",
    "ensembl_gene_id",
    "ENSEMBL"
  ),
  fallback_pattern = "^gene.*id$|^ensembl.*id$|^ensembl$",
  label = "annotation gene identifier"
)

annotation_symbol_col <- select_column(
  names(annotation),
  preferred = c(
    "SYMBOL",
    "symbol",
    "gene_symbol"
  ),
  fallback_pattern = "^symbol$|gene.*symbol",
  label = "gene symbol"
)

mapping <- data.frame(
  map_gene_id = normalise_ensembl(
    annotation[[annotation_gene_col]]
  ),
  symbol = trimws(
    as.character(
      annotation[[annotation_symbol_col]]
    )
  ),
  stringsAsFactors = FALSE
)

mapping <- mapping[
  !is.na(mapping$map_gene_id) &
    mapping$map_gene_id != "" &
    !is.na(mapping$symbol) &
    mapping$symbol != "",
  ,
  drop = FALSE
]

mapping <- unique(mapping)

symbols_per_gene <- tapply(
  mapping$symbol,
  mapping$map_gene_id,
  function(z) {
    length(unique(z))
  }
)

unambiguous_gene_ids <- names(
  symbols_per_gene[
    symbols_per_gene == 1L
  ]
)

ambiguous_gene_ids <- names(
  symbols_per_gene[
    symbols_per_gene > 1L
  ]
)

mapping_unique <- mapping[
  mapping$map_gene_id %in%
    unambiguous_gene_ids,
  ,
  drop = FALSE
]

mapping_unique <- mapping_unique[
  !duplicated(mapping_unique$map_gene_id),
  ,
  drop = FALSE
]

w5$map_gene_id <- normalise_ensembl(
  w5$gene_id
)

w5$SYMBOL <- mapping_unique$symbol[
  match(
    w5$map_gene_id,
    mapping_unique$map_gene_id
  )
]

cat("\n===== ANNOTATION =====\n")
cat(
  "Genes with mapped unambiguous symbol:",
  sum(!is.na(w5$SYMBOL)),
  "\n"
)

cat(
  "Ambiguous annotation gene IDs excluded:",
  length(ambiguous_gene_ids),
  "\n"
)

# ------------------------------------------------------------
# 7. Deterministic label selection
# ------------------------------------------------------------

label_candidates <- w5[
  w5$FDR_significant &
    !is.na(w5$SYMBOL) &
    w5$SYMBOL != "",
  ,
  drop = FALSE
]

higher_candidates <- label_candidates[
  label_candidates$logFC > 0,
  ,
  drop = FALSE
]

lower_candidates <- label_candidates[
  label_candidates$logFC < 0,
  ,
  drop = FALSE
]

higher_candidates <- higher_candidates[
  order(
    higher_candidates$adj.P.Val,
    -abs(higher_candidates$logFC),
    higher_candidates$gene_id
  ),
  ,
  drop = FALSE
]

lower_candidates <- lower_candidates[
  order(
    lower_candidates$adj.P.Val,
    -abs(lower_candidates$logFC),
    lower_candidates$gene_id
  ),
  ,
  drop = FALSE
]

labels <- rbind(
  head(higher_candidates, 5L),
  head(lower_candidates, 5L)
)

labels <- labels[
  order(
    labels$logFC
  ),
  ,
  drop = FALSE
]

stopifnot(
  nrow(labels) == 10L,
  sum(labels$logFC > 0) == 5L,
  sum(labels$logFC < 0) == 5L,
  all(labels$FDR_significant),
  !anyDuplicated(labels$gene_id),
  !any(is.na(labels$SYMBOL))
)

cat("\n===== LABEL RULE =====\n")
cat(
  "Five lowest-FDR mapped genes per direction.\n"
)

print(
  labels[
    ,
    c(
      "gene_id",
      "SYMBOL",
      "logFC",
      "adj.P.Val",
      "status"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Portfolio figure
# ------------------------------------------------------------

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

status_colours <- c(
  "Not FDR-significant" = "#B8B8B8",
  "Higher in IPF" = "#D55E00",
  "Lower in IPF" = "#0072B2"
)

status_alpha <- c(
  "Not FDR-significant" = 0.22,
  "Higher in IPF" = 0.58,
  "Lower in IPF" = 0.58
)

status_sizes <- c(
  "Not FDR-significant" = 0.85,
  "Higher in IPF" = 1.10,
  "Lower in IPF" = 1.10
)

subtitle_text <- paste0(
  "Donor-aware model adjusted for lung region",
  "  |  FDR < 0.05",
  "  |  ",
  format(
    n_genes,
    big.mark = ","
  ),
  " genes"
)

caption_text <- paste0(
  "Frozen W5 result. ",
  format(
    n_fdr,
    big.mark = ","
  ),
  " genes meet FDR < 0.05: ",
  format(
    n_higher,
    big.mark = ","
  ),
  " higher and ",
  format(
    n_lower,
    big.mark = ","
  ),
  " lower in IPF. ",
  "No effect-size threshold applied."
)

p <- ggplot(
  w5,
  aes(
    x = logFC,
    y = minus_log10_FDR
  )
) +
  geom_point(
    aes(
      colour = status,
      alpha = status,
      size = status
    )
  ) +
  geom_hline(
    yintercept = fdr_line,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey30"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    colour = "grey45"
  ) +
  ggrepel::geom_text_repel(
    data = labels,
    aes(
      x = logFC,
      y = minus_log10_FDR,
      label = SYMBOL
    ),
    inherit.aes = FALSE,
    size = 3.15,
    fontface = "bold",
    box.padding = 0.50,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.colour = "grey55",
    segment.size = 0.30,
    max.overlaps = Inf,
    force = 1.4,
    seed = 213001
  ) +
  scale_colour_manual(
    values = status_colours,
    drop = FALSE
  ) +
  scale_alpha_manual(
    values = status_alpha,
    drop = FALSE,
    guide = "none"
  ) +
  scale_size_manual(
    values = status_sizes,
    drop = FALSE,
    guide = "none"
  ) +
  labs(
    title = "Primary differential expression - IPF vs NDC",
    subtitle = subtitle_text,
    x = "log2 fold change (IPF relative to NDC)",
    y = expression(
      -log[10]("BH-adjusted P value")
    ),
    colour = "W5 FDR status",
    caption = caption_text
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        alpha = 1,
        size = 3
      )
    )
  ) +
  theme_minimal(
    base_size = 12.5
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    plot.subtitle = element_text(
      size = 11.5,
      margin = margin(
        b = 12
      )
    ),
    plot.caption = element_text(
      size = 9,
      colour = "grey35",
      hjust = 0,
      margin = margin(
        t = 12
      )
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.35,
      colour = "grey90"
    ),
    plot.margin = margin(
      15,
      35,
      15,
      15
    )
  )

ggsave(
  filename = output_png,
  plot = p,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 10,
  height = 7,
  units = "in",
  device = "pdf"
)

# ------------------------------------------------------------
# 9. Label audit
# ------------------------------------------------------------

label_audit <- labels[
  ,
  c(
    "gene_id",
    "SYMBOL",
    "logFC",
    "P.Value",
    "adj.P.Val",
    "status"
  ),
  drop = FALSE
]

label_audit$selection_rule <- paste0(
  "Top 5 smallest adjusted P values within ",
  "each FDR-significant direction; ",
  "ties by |logFC| then gene_id"
)

write.csv(
  label_audit,
  label_qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Figure QC record
# ------------------------------------------------------------

qc <- do.call(
  rbind,
  list(
    metric_row(
      "source_file",
      input_file
    ),
    metric_row(
      "annotation_file",
      annotation_file
    ),
    metric_row(
      "tested_genes",
      n_genes
    ),
    metric_row(
      "unique_gene_ids",
      length(unique(w5$gene_id))
    ),
    metric_row(
      "fdr_threshold",
      fdr_threshold
    ),
    metric_row(
      "fdr_significant_genes",
      n_fdr
    ),
    metric_row(
      "fdr_higher_in_IPF",
      n_higher
    ),
    metric_row(
      "fdr_lower_in_IPF",
      n_lower
    ),
    metric_row(
      "zero_adjusted_p_values",
      zero_adjusted_p
    ),
    metric_row(
      "unambiguous_symbols_mapped",
      sum(!is.na(w5$SYMBOL))
    ),
    metric_row(
      "ambiguous_mapping_ids_excluded",
      length(ambiguous_gene_ids)
    ),
    metric_row(
      "displayed_gene_labels",
      nrow(labels)
    ),
    metric_row(
      "displayed_higher_IPF_labels",
      sum(labels$logFC > 0)
    ),
    metric_row(
      "displayed_lower_IPF_labels",
      sum(labels$logFC < 0)
    ),
    metric_row(
      "effect_size_threshold_applied",
      FALSE
    ),
    metric_row(
      "differential_expression_refitted",
      FALSE
    ),
    metric_row(
      "P_values_recalculated",
      FALSE
    ),
    metric_row(
      "adjusted_P_values_recalculated",
      FALSE
    ),
    metric_row(
      "output_png",
      output_png
    ),
    metric_row(
      "output_pdf",
      output_pdf
    ),
    metric_row(
      "label_audit_file",
      label_qc_file
    )
  )
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Output verification
# ------------------------------------------------------------

stopifnot(
  file.exists(output_png),
  file.exists(output_pdf),
  file.exists(label_qc_file),
  file.exists(qc_file),
  file.info(output_png)$size > 0,
  file.info(output_pdf)$size > 0,
  nrow(label_audit) == 10L
)

cat("\n===== OUTPUTS =====\n")
cat(output_png, "\n")
cat(output_pdf, "\n")
cat(label_qc_file, "\n")
cat(qc_file, "\n")

cat("\n===== INFERENTIAL BOUNDARY =====\n")
cat("Differential-expression model refitted: FALSE\n")
cat("P values recalculated: FALSE\n")
cat("Adjusted P values recalculated: FALSE\n")
cat("Effect-size significance threshold introduced: FALSE\n")
cat("Frozen W5 FDR threshold retained: TRUE\n")

cat("\n============================================================\n")
cat("W8.2e PRIMARY DE PORTFOLIO VISUALIZATION: COMPLETED\n")
cat("============================================================\n")
