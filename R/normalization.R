# ==============================================================================
# GSE213001 W4 — low-count filtering and normalization
# ==============================================================================
#
# Current stage:
#   W4.1 — frozen low-count filtering rule
#
# IMPORTANT:
#   No normalization or differential-expression testing is performed yet.
#   The filtering rule is evaluated on the pre-specified primary cohort:
#
#       IPF vs NDC, known lung region
#
#   with lung location included in the design.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

required_packages <- c(
  "edgeR",
  "limma"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------------------

metadata_file <- file.path(
  "metadata",
  "sample_metadata.csv"
)

counts_file <- file.path(
  "data",
  "source",
  "GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz"
)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

if (!file.exists(counts_file)) {
  stop("Raw-count file not found: ", counts_file)
}


# ------------------------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("NA", "")
)


# ------------------------------------------------------------------------------
# Read raw counts
# ------------------------------------------------------------------------------

counts <- read.csv(
  gzfile(counts_file),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_id <- counts[[1]]

count_matrix <- as.matrix(
  counts[, -1, drop = FALSE]
)

storage.mode(count_matrix) <- "numeric"


# ------------------------------------------------------------------------------
# Raw-count integrity
# ------------------------------------------------------------------------------

if (anyNA(count_matrix)) {
  stop("Missing values detected in raw-count matrix.")
}

if (any(count_matrix < 0)) {
  stop("Negative raw counts detected.")
}

if (
  any(
    abs(count_matrix - round(count_matrix)) >
      .Machine$double.eps^0.5
  )
) {
  stop("Non-integer raw counts detected.")
}

if (anyDuplicated(gene_id)) {
  stop("Duplicated gene IDs detected.")
}

if (anyDuplicated(colnames(count_matrix))) {
  stop("Duplicated sample columns detected.")
}


# ------------------------------------------------------------------------------
# Verify full metadata/count alignment
# ------------------------------------------------------------------------------

if (!identical(
  colnames(count_matrix),
  metadata$sample_title
)) {
  stop(
    "Full metadata order does not match raw-count column order."
  )
}


# ------------------------------------------------------------------------------
# Define primary analysis cohort
# ------------------------------------------------------------------------------

primary_metadata <- metadata[
  metadata$diseasegroup %in% c(
    "IPF",
    "NDC"
  ) &
    !is.na(metadata$lunglocation),
  ,
  drop = FALSE
]

if (nrow(primary_metadata) != 101L) {
  stop(
    "Expected 101 primary-cohort samples; found ",
    nrow(primary_metadata),
    "."
  )
}

if (
  length(unique(primary_metadata$donorid)) != 34L
) {
  stop(
    "Expected 34 primary-cohort donors; found ",
    length(unique(primary_metadata$donorid)),
    "."
  )
}


# ------------------------------------------------------------------------------
# Extract matching raw-count columns
# ------------------------------------------------------------------------------

sample_index <- match(
  primary_metadata$sample_title,
  colnames(count_matrix)
)

if (anyNA(sample_index)) {
  stop(
    "At least one primary-cohort sample is absent ",
    "from the raw-count matrix."
  )
}

primary_counts <- count_matrix[
  ,
  sample_index,
  drop = FALSE
]

if (!identical(
  colnames(primary_counts),
  primary_metadata$sample_title
)) {
  stop(
    "Primary metadata/count alignment failed."
  )
}


# ------------------------------------------------------------------------------
# Pre-specified design variables
# ------------------------------------------------------------------------------

primary_metadata$diseasegroup <- factor(
  primary_metadata$diseasegroup,
  levels = c(
    "NDC",
    "IPF"
  )
)

primary_metadata$lunglocation <- factor(
  primary_metadata$lunglocation,
  levels = c(
    "Base",
    "Apex"
  )
)

design <- model.matrix(
  ~ diseasegroup + lunglocation,
  data = primary_metadata
)

rownames(design) <- primary_metadata$sample_title

if (qr(design)$rank != ncol(design)) {
  stop(
    "Primary design matrix is not full rank."
  )
}


# ------------------------------------------------------------------------------
# Construct edgeR object
# ------------------------------------------------------------------------------

y <- edgeR::DGEList(
  counts = primary_counts,
  genes = data.frame(
    gene_id = gene_id,
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------------------------
# Candidate low-count filtering rule
# ------------------------------------------------------------------------------

# Explicitly state edgeR filterByExpr defaults so the rule is reproducible.
#
# This filtering rule was reviewed and frozen on 2026-09-24.
# See docs/DECISIONS.md for the decision record.

filter_min_count <- 10
filter_min_total_count <- 15
filter_large_n <- 10
filter_min_prop <- 0.70

keep <- edgeR::filterByExpr(
  y,
  design = design,
  min.count = filter_min_count,
  min.total.count = filter_min_total_count,
  large.n = filter_large_n,
  min.prop = filter_min_prop
)


# ------------------------------------------------------------------------------
# Gene-level filtering diagnostics
# ------------------------------------------------------------------------------

raw_cpm <- edgeR::cpm(
  y,
  normalized.lib.sizes = FALSE
)

gene_filtering <- data.frame(
  gene_id = gene_id,
  keep_primary = keep,
  total_count_primary = rowSums(primary_counts),
  max_cpm_primary = apply(
    raw_cpm,
    1,
    max
  ),
  samples_cpm_ge1_primary = rowSums(
    raw_cpm >= 1
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gene_filtering,
  file.path(
    "metadata",
    "qc_gene_filtering.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Filtering summary
# ------------------------------------------------------------------------------

n_input <- length(keep)
n_retained <- sum(keep)
n_removed <- sum(!keep)

filter_summary <- data.frame(
  metric = c(
    "primary_samples",
    "primary_donors",
    "input_genes",
    "retained_genes",
    "removed_genes",
    "retained_percent",
    "min_count",
    "min_total_count",
    "large_n",
    "min_prop"
  ),
  value = c(
    nrow(primary_metadata),
    length(unique(primary_metadata$donorid)),
    n_input,
    n_retained,
    n_removed,
    round(
      100 * n_retained / n_input,
      2
    ),
    filter_min_count,
    filter_min_total_count,
    filter_large_n,
    filter_min_prop
  ),
  stringsAsFactors = FALSE
)

write.csv(
  filter_summary,
  file.path(
    "metadata",
    "qc_gene_filtering_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Console report
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("GSE213001 W4.1 — FROZEN LOW-COUNT FILTERING\n")
cat("============================================================\n")


cat("\nPrimary cohort:\n")
cat(
  "Samples: ",
  nrow(primary_metadata),
  "\n",
  sep = ""
)
cat(
  "Donors:  ",
  length(unique(primary_metadata$donorid)),
  "\n",
  sep = ""
)


cat("\nDisease group x lung location:\n")

print(
  addmargins(
    table(
      primary_metadata$diseasegroup,
      primary_metadata$lunglocation
    )
  )
)


cat("\nDesign matrix columns:\n")
print(
  colnames(design)
)

cat(
  "Design rank: ",
  qr(design)$rank,
  " / ",
  ncol(design),
  "\n",
  sep = ""
)


cat("\nfilterByExpr candidate settings:\n")
cat(
  "min.count       = ",
  filter_min_count,
  "\n",
  sep = ""
)
cat(
  "min.total.count = ",
  filter_min_total_count,
  "\n",
  sep = ""
)
cat(
  "large.n         = ",
  filter_large_n,
  "\n",
  sep = ""
)
cat(
  "min.prop        = ",
  filter_min_prop,
  "\n",
  sep = ""
)


cat("\nGene filtering result:\n")
cat(
  "Input genes:    ",
  n_input,
  "\n",
  sep = ""
)
cat(
  "Retained genes: ",
  n_retained,
  "\n",
  sep = ""
)
cat(
  "Removed genes:  ",
  n_removed,
  "\n",
  sep = ""
)
cat(
  "Retention:      ",
  round(
    100 * n_retained / n_input,
    2
  ),
  "%\n",
  sep = ""
)


cat("\nRaw-expression diagnostics:\n")

cat(
  "Genes with total count = 0 in primary cohort: ",
  sum(
    gene_filtering$total_count_primary == 0
  ),
  "\n",
  sep = ""
)

cat(
  "Genes with CPM >= 1 in at least 1 sample:      ",
  sum(
    gene_filtering$samples_cpm_ge1_primary >= 1
  ),
  "\n",
  sep = ""
)

cat(
  "Genes with CPM >= 1 in at least 10 samples:    ",
  sum(
    gene_filtering$samples_cpm_ge1_primary >= 10
  ),
  "\n",
  sep = ""
)


cat("\nOutputs:\n")
cat("  metadata/qc_gene_filtering.csv\n")
cat("  metadata/qc_gene_filtering_summary.csv\n")


cat("\n")
cat("Filtering rule status: FROZEN / ACCEPTED\n")

cat("No normalization was performed.\n")
cat("No differential-expression testing was performed.\n")

cat("============================================================\n")


# ==============================================================================
# W4.2 — TMM normalization
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("GSE213001 W4.2 — TMM NORMALIZATION\n")
cat("============================================================\n")


# ------------------------------------------------------------------------------
# Retain genes passing the frozen filter
# ------------------------------------------------------------------------------

y_filtered <- y[
  keep,
  ,
  keep.lib.sizes = FALSE
]

if (nrow(y_filtered) != 15012L) {
  stop(
    "Frozen filtering result expected 15012 genes; found ",
    nrow(y_filtered),
    "."
  )
}


# ------------------------------------------------------------------------------
# TMM normalization
# ------------------------------------------------------------------------------

y_tmm <- edgeR::calcNormFactors(
  y_filtered,
  method = "TMM"
)

norm_factors <- y_tmm$samples$norm.factors

if (
  anyNA(norm_factors) ||
  any(!is.finite(norm_factors)) ||
  any(norm_factors <= 0)
) {
  stop(
    "Invalid TMM normalization factors detected."
  )
}


# ------------------------------------------------------------------------------
# Library-size diagnostics
# ------------------------------------------------------------------------------

full_library_size <- colSums(
  primary_counts
)

filtered_library_size <- y_tmm$samples$lib.size

effective_library_size <-
  filtered_library_size * norm_factors

normalization_metrics <- data.frame(
  sample_title = primary_metadata$sample_title,
  donorid = primary_metadata$donorid,
  diseasegroup = primary_metadata$diseasegroup,
  lunglocation = primary_metadata$lunglocation,
  full_library_size = as.numeric(full_library_size),
  filtered_library_size = as.numeric(filtered_library_size),
  retained_count_fraction =
    as.numeric(filtered_library_size / full_library_size),
  tmm_norm_factor = as.numeric(norm_factors),
  effective_library_size =
    as.numeric(effective_library_size),
  stringsAsFactors = FALSE
)

write.csv(
  normalization_metrics,
  file.path(
    "metadata",
    "qc_tmm_normalization.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Summary metrics
# ------------------------------------------------------------------------------

normalization_summary <- data.frame(
  metric = c(
    "samples",
    "retained_genes",
    "tmm_factor_min",
    "tmm_factor_q1",
    "tmm_factor_median",
    "tmm_factor_mean",
    "tmm_factor_q3",
    "tmm_factor_max",
    "tmm_factor_product",
    "retained_count_fraction_min",
    "retained_count_fraction_median",
    "retained_count_fraction_max"
  ),
  value = c(
    ncol(y_tmm),
    nrow(y_tmm),
    min(norm_factors),
    unname(
      quantile(
        norm_factors,
        0.25
      )
    ),
    median(norm_factors),
    mean(norm_factors),
    unname(
      quantile(
        norm_factors,
        0.75
      )
    ),
    max(norm_factors),
    prod(norm_factors),
    min(
      normalization_metrics$retained_count_fraction
    ),
    median(
      normalization_metrics$retained_count_fraction
    ),
    max(
      normalization_metrics$retained_count_fraction
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  normalization_summary,
  file.path(
    "metadata",
    "qc_tmm_normalization_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Console diagnostics
# ------------------------------------------------------------------------------

cat("\nFiltered DGEList:\n")
cat(
  "Genes:   ",
  nrow(y_tmm),
  "\n",
  sep = ""
)
cat(
  "Samples: ",
  ncol(y_tmm),
  "\n",
  sep = ""
)


cat("\nTMM normalization factors:\n")
print(
  summary(norm_factors)
)

cat(
  "\nProduct of normalization factors: ",
  format(
    prod(norm_factors),
    digits = 8
  ),
  "\n",
  sep = ""
)


cat("\nSmallest TMM factors:\n")

print(
  head(
    normalization_metrics[
      order(
        normalization_metrics$tmm_norm_factor
      ),
      c(
        "sample_title",
        "donorid",
        "diseasegroup",
        "lunglocation",
        "full_library_size",
        "tmm_norm_factor",
        "effective_library_size"
      )
    ],
    10
  ),
  row.names = FALSE
)


cat("\nLargest TMM factors:\n")

print(
  head(
    normalization_metrics[
      order(
        -normalization_metrics$tmm_norm_factor
      ),
      c(
        "sample_title",
        "donorid",
        "diseasegroup",
        "lunglocation",
        "full_library_size",
        "tmm_norm_factor",
        "effective_library_size"
      )
    ],
    10
  ),
  row.names = FALSE
)


cat("\nFiltered-count retention per sample:\n")

print(
  summary(
    normalization_metrics$retained_count_fraction
  )
)


cat("\nEffective library-size summary:\n")

print(
  summary(
    normalization_metrics$effective_library_size
  )
)


cat("\nOutputs:\n")
cat("  metadata/qc_tmm_normalization.csv\n")
cat("  metadata/qc_tmm_normalization_summary.csv\n")

cat("\n")
cat(
  "TMM factors were computed independently from raw counts; ",
  "GEO-reported normalization factors were not reused.\n"
)

cat("No PCA/MDS was performed.\n")
cat("No differential-expression testing was performed.\n")

cat("============================================================\n")


# ------------------------------------------------------------------------------
# W4.2 decision
# ------------------------------------------------------------------------------

# TMM normalization reviewed and accepted on 2026-09-24.
#
# ALF009B has an unusually low TMM factor because its library is strongly
# compositionally concentrated. This sample is retained and remains on the
# technical watchlist for downstream ordination/model diagnostics.
#
# See docs/DECISIONS.md.


# ==============================================================================
# W4.3 — normalized logCPM and unsupervised ordination
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("GSE213001 W4.3 — LOGCPM + PCA/MDS\n")
cat("============================================================\n")


# ------------------------------------------------------------------------------
# TMM-normalized log2 CPM
# ------------------------------------------------------------------------------

logcpm <- edgeR::cpm(
  y_tmm,
  log = TRUE,
  prior.count = 2,
  normalized.lib.sizes = TRUE
)

if (
  nrow(logcpm) != 15012L ||
  ncol(logcpm) != 101L
) {
  stop(
    "Unexpected normalized logCPM dimensions."
  )
}

if (anyNA(logcpm) || any(!is.finite(logcpm))) {
  stop(
    "Non-finite values detected in normalized logCPM matrix."
  )
}


# ------------------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------------------

pca <- prcomp(
  t(logcpm),
  center = TRUE,
  scale. = FALSE
)

pca_variance <- 100 *
  (pca$sdev^2 / sum(pca$sdev^2))

pca_coordinates <- data.frame(
  sample_title = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4],
  stringsAsFactors = FALSE
)

pca_coordinates <- merge(
  primary_metadata,
  pca_coordinates,
  by = "sample_title",
  sort = FALSE
)

pca_coordinates <- pca_coordinates[
  match(
    primary_metadata$sample_title,
    pca_coordinates$sample_title
  ),
  ,
  drop = FALSE
]

if (!identical(
  pca_coordinates$sample_title,
  primary_metadata$sample_title
)) {
  stop("PCA coordinate/metadata alignment failed.")
}


# ------------------------------------------------------------------------------
# edgeR MDS
# ------------------------------------------------------------------------------

mds <- limma::plotMDS(
  y_tmm,
  top = 500,
  gene.selection = "common",
  plot = FALSE
)

mds_coordinates <- data.frame(
  sample_title = colnames(y_tmm),
  MDS1 = mds$x,
  MDS2 = mds$y,
  stringsAsFactors = FALSE
)

mds_coordinates <- merge(
  primary_metadata,
  mds_coordinates,
  by = "sample_title",
  sort = FALSE
)

mds_coordinates <- mds_coordinates[
  match(
    primary_metadata$sample_title,
    mds_coordinates$sample_title
  ),
  ,
  drop = FALSE
]

if (!identical(
  mds_coordinates$sample_title,
  primary_metadata$sample_title
)) {
  stop("MDS coordinate/metadata alignment failed.")
}


# ------------------------------------------------------------------------------
# W3 technical watchlist
# ------------------------------------------------------------------------------

technical_watchlist <- c(
  "ALF009B",
  "ALF012F",
  "ALF030A",
  "ALF048D",
  "ALF017A",
  "ALF024A"
)

pca_coordinates$technical_watchlist <-
  pca_coordinates$sample_title %in% technical_watchlist

mds_coordinates$technical_watchlist <-
  mds_coordinates$sample_title %in% technical_watchlist


# ------------------------------------------------------------------------------
# Export ordination data
# ------------------------------------------------------------------------------

write.csv(
  pca_coordinates,
  file.path(
    "metadata",
    "qc_pca_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  mds_coordinates,
  file.path(
    "metadata",
    "qc_mds_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    component = paste0(
      "PC",
      seq_along(pca_variance)
    ),
    percent_variance = pca_variance,
    stringsAsFactors = FALSE
  ),
  file.path(
    "metadata",
    "qc_pca_variance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Console diagnostics
# ------------------------------------------------------------------------------

cat("\nNormalized logCPM matrix:\n")
cat(
  "Genes:   ",
  nrow(logcpm),
  "\n",
  sep = ""
)
cat(
  "Samples: ",
  ncol(logcpm),
  "\n",
  sep = ""
)


cat("\nPCA variance explained:\n")

for (i in 1:6) {
  cat(
    sprintf(
      "PC%d: %.2f%%\n",
      i,
      pca_variance[i]
    )
  )
}


cat("\nPCA coordinate ranges:\n")

print(
  summary(
    pca_coordinates[
      ,
      c(
        "PC1",
        "PC2",
        "PC3",
        "PC4"
      )
    ]
  )
)


cat("\nMDS coordinate ranges:\n")

print(
  summary(
    mds_coordinates[
      ,
      c(
        "MDS1",
        "MDS2"
      )
    ]
  )
)


cat("\nTechnical watchlist in primary cohort:\n")

print(
  pca_coordinates[
    pca_coordinates$technical_watchlist,
    c(
      "sample_title",
      "donorid",
      "diseasegroup",
      "lunglocation",
      "rin",
      "PC1",
      "PC2"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


cat("\nOutputs:\n")
cat("  metadata/qc_pca_coordinates.csv\n")
cat("  metadata/qc_mds_coordinates.csv\n")
cat("  metadata/qc_pca_variance.csv\n")

cat("\n")
cat("These ordinations are exploratory QC outputs.\n")
cat("No sample exclusion was triggered.\n")
cat("No differential-expression testing was performed.\n")
cat("============================================================\n")


# ==============================================================================
# W4.3b — ordination figures
# ==============================================================================

dir.create(
  file.path("figures", "qc"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Plot helper
# ------------------------------------------------------------------------------

save_qc_plot <- function(plot_object, filename) {
  ggplot2::ggsave(
    filename = file.path(
      "figures",
      "qc",
      filename
    ),
    plot = plot_object,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# ------------------------------------------------------------------------------
# PCA — disease group
# ------------------------------------------------------------------------------

p_pca_disease <- ggplot2::ggplot(
  pca_coordinates,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    shape = lunglocation,
    color = diseasegroup
  )
) +
  ggplot2::geom_point(
    size = 3,
    alpha = 0.80
  ) +
  ggplot2::labs(
    title = "PCA of TMM-normalized RNA-seq data",
    subtitle = "Primary cohort: IPF vs NDC",
    x = sprintf(
      "PC1 (%.2f%%)",
      pca_variance[1]
    ),
    y = sprintf(
      "PC2 (%.2f%%)",
      pca_variance[2]
    ),
    color = "Disease group",
    shape = "Lung region"
  ) +
  ggplot2::theme_bw(base_size = 12)

save_qc_plot(
  p_pca_disease,
  "pca_disease_region.png"
)


# ------------------------------------------------------------------------------
# PCA — technical watchlist
# ------------------------------------------------------------------------------

p_pca_watchlist <- ggplot2::ggplot(
  pca_coordinates,
  ggplot2::aes(
    x = PC1,
    y = PC2
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = diseasegroup
    ),
    size = 3,
    alpha = 0.65
  ) +
  ggplot2::geom_point(
    data = pca_coordinates[
      pca_coordinates$technical_watchlist,
      ,
      drop = FALSE
    ],
    shape = 21,
    size = 5,
    stroke = 1.2
  ) +
  ggplot2::geom_text(
    data = pca_coordinates[
      pca_coordinates$technical_watchlist,
      ,
      drop = FALSE
    ],
    ggplot2::aes(
      label = sample_title
    ),
    nudge_y = 4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "PCA with technical watchlist",
    subtitle = "Watchlist samples are outlined and labelled",
    x = sprintf(
      "PC1 (%.2f%%)",
      pca_variance[1]
    ),
    y = sprintf(
      "PC2 (%.2f%%)",
      pca_variance[2]
    ),
    color = "Disease group"
  ) +
  ggplot2::theme_bw(base_size = 12)

save_qc_plot(
  p_pca_watchlist,
  "pca_technical_watchlist.png"
)


# ------------------------------------------------------------------------------
# PCA — processing date
# ------------------------------------------------------------------------------

p_pca_processing <- ggplot2::ggplot(
  pca_coordinates,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    color = processingdate
  )
) +
  ggplot2::geom_point(
    size = 3,
    alpha = 0.80
  ) +
  ggplot2::labs(
    title = "PCA by processing date",
    subtitle = "Screen for technical/batch structure",
    x = sprintf(
      "PC1 (%.2f%%)",
      pca_variance[1]
    ),
    y = sprintf(
      "PC2 (%.2f%%)",
      pca_variance[2]
    ),
    color = "Processing date"
  ) +
  ggplot2::theme_bw(base_size = 12)

save_qc_plot(
  p_pca_processing,
  "pca_processing_date.png"
)


# ------------------------------------------------------------------------------
# PCA — RIN
# ------------------------------------------------------------------------------

p_pca_rin <- ggplot2::ggplot(
  pca_coordinates,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    color = rin
  )
) +
  ggplot2::geom_point(
    size = 3,
    alpha = 0.80
  ) +
  ggplot2::labs(
    title = "PCA by RNA integrity",
    subtitle = "Continuous RIN overlay",
    x = sprintf(
      "PC1 (%.2f%%)",
      pca_variance[1]
    ),
    y = sprintf(
      "PC2 (%.2f%%)",
      pca_variance[2]
    ),
    color = "RIN"
  ) +
  ggplot2::theme_bw(base_size = 12)

save_qc_plot(
  p_pca_rin,
  "pca_rin.png"
)


# ------------------------------------------------------------------------------
# edgeR MDS — disease group / lung region
# ------------------------------------------------------------------------------

p_mds_disease <- ggplot2::ggplot(
  mds_coordinates,
  ggplot2::aes(
    x = MDS1,
    y = MDS2,
    color = diseasegroup,
    shape = lunglocation
  )
) +
  ggplot2::geom_point(
    size = 3,
    alpha = 0.80
  ) +
  ggplot2::labs(
    title = "edgeR multidimensional scaling",
    subtitle = "Top 500 common genes",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2",
    color = "Disease group",
    shape = "Lung region"
  ) +
  ggplot2::theme_bw(base_size = 12)

save_qc_plot(
  p_mds_disease,
  "mds_disease_region.png"
)


# ------------------------------------------------------------------------------
# edgeR MDS — technical watchlist
# ------------------------------------------------------------------------------

p_mds_watchlist <- ggplot2::ggplot(
  mds_coordinates,
  ggplot2::aes(
    x = MDS1,
    y = MDS2
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = diseasegroup
    ),
    size = 3,
    alpha = 0.65
  ) +
  ggplot2::geom_point(
    data = mds_coordinates[
      mds_coordinates$technical_watchlist,
      ,
      drop = FALSE
    ],
    shape = 21,
    size = 5,
    stroke = 1.2
  ) +
  ggplot2::geom_text(
    data = mds_coordinates[
      mds_coordinates$technical_watchlist,
      ,
      drop = FALSE
    ],
    ggplot2::aes(
      label = sample_title
    ),
    nudge_y = 0.15,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "MDS with technical watchlist",
    subtitle = "Technical-screening samples remain included",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2",
    color = "Disease group"
  ) +
  ggplot2::theme_bw(base_size = 12)

save_qc_plot(
  p_mds_watchlist,
  "mds_technical_watchlist.png"
)


cat("\nOrdination figures written:\n")
cat("  figures/qc/pca_disease_region.png\n")
cat("  figures/qc/pca_technical_watchlist.png\n")
cat("  figures/qc/pca_processing_date.png\n")
cat("  figures/qc/pca_rin.png\n")
cat("  figures/qc/mds_disease_region.png\n")
cat("  figures/qc/mds_technical_watchlist.png\n")

cat("\nNo sample exclusions were made from ordination results.\n")


# ==============================================================================
# W4.3c — quantitative ordination audit
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("GSE213001 W4.3c — QUANTITATIVE ORDINATION AUDIT\n")
cat("============================================================\n")


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

r2_from_model <- function(model) {
  summary(model)$r.squared
}


partial_r2 <- function(reduced_model, full_model) {

  sse_reduced <- sum(residuals(reduced_model)^2)
  sse_full <- sum(residuals(full_model)^2)

  if (sse_reduced == 0) {
    return(NA_real_)
  }

  1 - (sse_full / sse_reduced)
}


# ------------------------------------------------------------------------------
# Donor-level PCA centroids
# ------------------------------------------------------------------------------

donors <- unique(pca_coordinates$donorid)

donor_pca <- do.call(
  rbind,
  lapply(
    donors,
    function(donor) {

      x <- pca_coordinates[
        pca_coordinates$donorid == donor,
        ,
        drop = FALSE
      ]

      disease_values <- unique(x$diseasegroup)

      if (length(disease_values) != 1L) {
        stop(
          "Disease-group inconsistency within donor: ",
          donor
        )
      }

      age_values <- unique(
        x$age[
          !is.na(x$age)
        ]
      )

      if (length(age_values) > 1L) {
        stop(
          "Age inconsistency within donor: ",
          donor
        )
      }

      age_value <- if (
        length(age_values) == 0L
      ) {
        NA_real_
      } else {
        age_values[1]
      }

      data.frame(
        donorid = donor,
        diseasegroup = disease_values[1],
        age = age_value,
        PC1 = mean(x$PC1),
        PC2 = mean(x$PC2),
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------------------------
# Disease structure at donor level
# ------------------------------------------------------------------------------

disease_pc1 <- lm(
  PC1 ~ diseasegroup,
  data = donor_pca
)

disease_pc2 <- lm(
  PC2 ~ diseasegroup,
  data = donor_pca
)


cat("\nDONOR-LEVEL DISEASE ASSOCIATION\n")

cat(
  sprintf(
    "PC1 variance explained by disease group: %.3f\n",
    r2_from_model(disease_pc1)
  )
)

cat(
  sprintf(
    "PC2 variance explained by disease group: %.3f\n",
    r2_from_model(disease_pc2)
  )
)


# ------------------------------------------------------------------------------
# Age audit at donor level
# ------------------------------------------------------------------------------

donor_age <- donor_pca[
  complete.cases(
    donor_pca[
      ,
      c(
        "age",
        "diseasegroup",
        "PC1",
        "PC2"
      )
    ]
  ),
  ,
  drop = FALSE
]


age_pc1 <- lm(
  PC1 ~ age,
  data = donor_age
)

age_pc2 <- lm(
  PC2 ~ age,
  data = donor_age
)


age_disease_pc1 <- lm(
  PC1 ~ age + diseasegroup,
  data = donor_age
)

age_disease_pc2 <- lm(
  PC2 ~ age + diseasegroup,
  data = donor_age
)


cat("\nDONOR-LEVEL AGE ASSOCIATION\n")

cat(
  sprintf(
    "PC1 variance explained by age alone: %.3f\n",
    r2_from_model(age_pc1)
  )
)

cat(
  sprintf(
    "PC2 variance explained by age alone: %.3f\n",
    r2_from_model(age_pc2)
  )
)


cat("\nDISEASE CONTRIBUTION AFTER AGE\n")

cat(
  sprintf(
    "PC1 partial R2 for disease after age: %.3f\n",
    partial_r2(
      age_pc1,
      age_disease_pc1
    )
  )
)

cat(
  sprintf(
    "PC2 partial R2 for disease after age: %.3f\n",
    partial_r2(
      age_pc2,
      age_disease_pc2
    )
  )
)


# ------------------------------------------------------------------------------
# Within-donor associations
# ------------------------------------------------------------------------------

within_donor_partial_r2 <- function(
  response,
  covariate,
  data
) {

  variables <- c(
    response,
    "donorid",
    covariate
  )

  x <- data[
    complete.cases(
      data[
        ,
        variables,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  reduced_formula <- reformulate(
    "donorid",
    response = response
  )

  full_formula <- reformulate(
    c(
      "donorid",
      covariate
    ),
    response = response
  )

  reduced <- lm(
    reduced_formula,
    data = x
  )

  full <- lm(
    full_formula,
    data = x
  )

  if (full$rank == reduced$rank) {
    return(NA_real_)
  }

  partial_r2(
    reduced,
    full
  )
}


variables_to_screen <- c(
  "lunglocation",
  "processingdate",
  "rin"
)

ordination_audit <- data.frame(
  variable = character(),
  PC1_partial_R2 = numeric(),
  PC2_partial_R2 = numeric(),
  stringsAsFactors = FALSE
)


for (variable in variables_to_screen) {

  ordination_audit <- rbind(
    ordination_audit,
    data.frame(
      variable = variable,
      PC1_partial_R2 =
        within_donor_partial_r2(
          "PC1",
          variable,
          pca_coordinates
        ),
      PC2_partial_R2 =
        within_donor_partial_r2(
          "PC2",
          variable,
          pca_coordinates
        ),
      stringsAsFactors = FALSE
    )
  )
}


cat("\nWITHIN-DONOR CONTRIBUTIONS\n")
print(
  ordination_audit,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Descriptive processing-date screen
# ------------------------------------------------------------------------------

processing_pc1 <- lm(
  PC1 ~ processingdate,
  data = pca_coordinates
)

processing_pc2 <- lm(
  PC2 ~ processingdate,
  data = pca_coordinates
)


cat("\nSAMPLE-LEVEL PROCESSING-DATE SCREEN\n")
cat("(descriptive only; repeated samples are not independent)\n")

cat(
  sprintf(
    "PC1 R2: %.3f\n",
    r2_from_model(processing_pc1)
  )
)

cat(
  sprintf(
    "PC2 R2: %.3f\n",
    r2_from_model(processing_pc2)
  )
)


# ------------------------------------------------------------------------------
# Save audit results
# ------------------------------------------------------------------------------

donor_results <- data.frame(
  metric = c(
    "PC1_disease_R2_donor_level",
    "PC2_disease_R2_donor_level",
    "PC1_age_R2_donor_level",
    "PC2_age_R2_donor_level",
    "PC1_disease_partial_R2_after_age",
    "PC2_disease_partial_R2_after_age"
  ),
  value = c(
    r2_from_model(disease_pc1),
    r2_from_model(disease_pc2),
    r2_from_model(age_pc1),
    r2_from_model(age_pc2),
    partial_r2(
      age_pc1,
      age_disease_pc1
    ),
    partial_r2(
      age_pc2,
      age_disease_pc2
    )
  )
)


write.csv(
  donor_results,
  file.path(
    "metadata",
    "qc_ordination_donor_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  ordination_audit,
  file.path(
    "metadata",
    "qc_ordination_within_donor_associations.csv"
  ),
  row.names = FALSE
)


cat("\nOutputs:\n")
cat(
  "  metadata/qc_ordination_donor_associations.csv\n"
)
cat(
  "  metadata/qc_ordination_within_donor_associations.csv\n"
)

cat("\n")
cat(
  "These metrics are exploratory QC summaries, not formal ",
  "differential-expression tests.\n",
  sep = ""
)

cat("============================================================\n")
