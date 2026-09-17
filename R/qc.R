# ==============================================================================
# GSE213001 metadata and library-level QC
# ==============================================================================
#
# Purpose:
#   Perform descriptive metadata and library-level quality checks before
#   filtering, normalization, PCA/MDS or differential-expression analysis.
#
# Inputs:
#   metadata/sample_metadata.csv
#   data/source/GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz
#
# Outputs:
#   metadata/qc_missingness.csv
#   metadata/qc_numeric_summary.csv
#   metadata/qc_categorical_summary.csv
#   metadata/qc_library_metrics.csv
#
# Important:
#   This script is descriptive.
#   It does NOT exclude samples.
#
# ==============================================================================


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
# Read data
# ------------------------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("NA", "")
)

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

sample_names <- colnames(count_matrix)


# ------------------------------------------------------------------------------
# Integrity checks
# ------------------------------------------------------------------------------

if (!identical(
  sample_names,
  metadata$sample_title
)) {
  stop(
    "Metadata order does not match raw-count column order."
  )
}

if (anyNA(count_matrix)) {
  stop("Missing values detected in raw-count matrix.")
}

if (any(count_matrix < 0)) {
  stop("Negative raw counts detected.")
}

if (anyDuplicated(sample_names)) {
  stop("Duplicated sample names detected.")
}

if (anyDuplicated(gene_id)) {
  stop("Duplicated gene identifiers detected.")
}


# ------------------------------------------------------------------------------
# Missingness summary
# ------------------------------------------------------------------------------

missingness <- data.frame(
  variable = names(metadata),
  n_missing = vapply(
    metadata,
    function(x) sum(is.na(x)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

missingness$percent_missing <- round(
  100 * missingness$n_missing / nrow(metadata),
  2
)

missingness <- missingness[
  order(
    -missingness$n_missing,
    missingness$variable
  ),
  ,
  drop = FALSE
]

write.csv(
  missingness,
  file.path(
    "metadata",
    "qc_missingness.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Numeric-variable summary
# ------------------------------------------------------------------------------

numeric_variables <- names(metadata)[
  vapply(metadata, is.numeric, logical(1))
]

numeric_summary_list <- lapply(
  numeric_variables,
  function(variable) {

    x <- metadata[[variable]]

    observed <- x[!is.na(x)]

    data.frame(
      variable = variable,
      n = length(observed),
      n_missing = sum(is.na(x)),
      mean = if (length(observed)) mean(observed) else NA_real_,
      sd = if (length(observed) > 1) sd(observed) else NA_real_,
      median = if (length(observed)) median(observed) else NA_real_,
      min = if (length(observed)) min(observed) else NA_real_,
      max = if (length(observed)) max(observed) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
)

numeric_summary <- do.call(
  rbind,
  numeric_summary_list
)

write.csv(
  numeric_summary,
  file.path(
    "metadata",
    "qc_numeric_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Categorical-variable summary
# ------------------------------------------------------------------------------

categorical_variables <- names(metadata)[
  !vapply(metadata, is.numeric, logical(1))
]

categorical_summary_list <- lapply(
  categorical_variables,
  function(variable) {

    values <- metadata[[variable]]

    values_display <- ifelse(
      is.na(values),
      "<NA>",
      values
    )

    tab <- table(values_display)

    data.frame(
      variable = variable,
      level = names(tab),
      n = as.integer(tab),
      percent = round(
        100 * as.integer(tab) / length(values),
        2
      ),
      stringsAsFactors = FALSE
    )
  }
)

categorical_summary <- do.call(
  rbind,
  categorical_summary_list
)

write.csv(
  categorical_summary,
  file.path(
    "metadata",
    "qc_categorical_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Library-level QC metrics
# ------------------------------------------------------------------------------

library_size <- colSums(count_matrix)

detected_genes_gt0 <- colSums(
  count_matrix > 0
)

detected_genes_ge10 <- colSums(
  count_matrix >= 10
)

zero_fraction <- colMeans(
  count_matrix == 0
)

library_metrics <- data.frame(
  sample_title = sample_names,
  geo_accession = metadata$geo_accession,
  donorid = metadata$donorid,
  diseasegroup = metadata$diseasegroup,
  lunglocation = metadata$lunglocation,
  rin = metadata$rin,
  library_size = as.numeric(library_size),
  log10_library_size = log10(
    as.numeric(library_size)
  ),
  detected_genes_gt0 = as.integer(
    detected_genes_gt0
  ),
  detected_genes_ge10 = as.integer(
    detected_genes_ge10
  ),
  zero_fraction = as.numeric(
    zero_fraction
  ),
  stringsAsFactors = FALSE
)

write.csv(
  library_metrics,
  file.path(
    "metadata",
    "qc_library_metrics.csv"
  ),
  row.names = FALSE,
  na = "NA"
)


# ------------------------------------------------------------------------------
# Primary cohort summary
# ------------------------------------------------------------------------------

primary <- metadata[
  metadata$diseasegroup %in% c(
    "IPF",
    "NDC"
  ),
  ,
  drop = FALSE
]

primary_known <- primary[
  !is.na(primary$lunglocation),
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# Console report
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("GSE213001 W3 metadata/library QC\n")
cat("============================================================\n")

cat("\nDataset:\n")
cat("Samples:          ", nrow(metadata), "\n")
cat(
  "Unique donors:    ",
  length(unique(metadata$donorid)),
  "\n"
)
cat("Genes:            ", nrow(count_matrix), "\n")

cat("\nPrimary analysis cohort:\n")
cat(
  "IPF + NDC:        ",
  nrow(primary),
  "samples\n"
)
cat(
  "Known region:     ",
  nrow(primary_known),
  "samples\n"
)
cat(
  "Unique donors:    ",
  length(unique(primary_known$donorid)),
  "\n"
)

cat("\nVariables with missing values:\n")

missing_nonzero <- missingness[
  missingness$n_missing > 0,
  ,
  drop = FALSE
]

if (nrow(missing_nonzero) == 0L) {
  cat("None\n")
} else {
  print(
    missing_nonzero,
    row.names = FALSE
  )
}


cat("\nDisease groups:\n")
print(
  table(
    metadata$diseasegroup,
    useNA = "ifany"
  )
)


cat("\nLung locations:\n")
print(
  table(
    metadata$lunglocation,
    useNA = "ifany"
  )
)


cat("\nLibrary-size summary:\n")
print(
  summary(
    library_metrics$library_size
  )
)


cat("\nDetected genes (>0 counts) per sample:\n")
print(
  summary(
    library_metrics$detected_genes_gt0
  )
)


cat("\nDetected genes (>=10 counts) per sample:\n")
print(
  summary(
    library_metrics$detected_genes_ge10
  )
)


cat("\nZero-count fraction per sample:\n")
print(
  summary(
    library_metrics$zero_fraction
  )
)


cat("\nSmallest libraries:\n")

print(
  head(
    library_metrics[
      order(
        library_metrics$library_size
      ),
      c(
        "sample_title",
        "donorid",
        "diseasegroup",
        "lunglocation",
        "rin",
        "library_size",
        "detected_genes_gt0",
        "zero_fraction"
      ),
      drop = FALSE
    ],
    10
  ),
  row.names = FALSE
)


cat("\nLargest libraries:\n")

print(
  head(
    library_metrics[
      order(
        -library_metrics$library_size
      ),
      c(
        "sample_title",
        "donorid",
        "diseasegroup",
        "lunglocation",
        "rin",
        "library_size",
        "detected_genes_gt0",
        "zero_fraction"
      ),
      drop = FALSE
    ],
    10
  ),
  row.names = FALSE
)


cat("\nOutputs:\n")
cat("  metadata/qc_missingness.csv\n")
cat("  metadata/qc_numeric_summary.csv\n")
cat("  metadata/qc_categorical_summary.csv\n")
cat("  metadata/qc_library_metrics.csv\n")

cat("\n")
cat("No samples were excluded by this script.\n")

cat("============================================================\n")


# ------------------------------------------------------------------------------
# Primary-cohort technical balance
# ------------------------------------------------------------------------------

primary_metrics <- library_metrics[
  library_metrics$diseasegroup %in% c("IPF", "NDC") &
    !is.na(library_metrics$lunglocation),
  ,
  drop = FALSE
]

summarise_metric <- function(data, group_variable, metric) {

  groups <- unique(data[[group_variable]])

  result <- lapply(
    sort(groups),
    function(group_value) {

      x <- data[
        data[[group_variable]] == group_value,
        metric
      ]

      x <- x[!is.na(x)]

      data.frame(
        grouping_variable = group_variable,
        group = group_value,
        metric = metric,
        n = length(x),
        mean = mean(x),
        sd = if (length(x) > 1) sd(x) else NA_real_,
        median = median(x),
        q1 = as.numeric(quantile(x, 0.25)),
        q3 = as.numeric(quantile(x, 0.75)),
        min = min(x),
        max = max(x),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, result)
}


balance_metrics <- c(
  "library_size",
  "rin",
  "detected_genes_gt0",
  "detected_genes_ge10",
  "zero_fraction"
)

balance_list <- list()

for (group_variable in c(
  "diseasegroup",
  "lunglocation"
)) {

  for (metric in balance_metrics) {

    balance_list[[length(balance_list) + 1L]] <-
      summarise_metric(
        primary_metrics,
        group_variable,
        metric
      )
  }
}

technical_balance <- do.call(
  rbind,
  balance_list
)

write.csv(
  technical_balance,
  file.path(
    "metadata",
    "qc_primary_technical_balance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Robust technical-outlier screening
# ------------------------------------------------------------------------------

robust_z <- function(x) {

  med <- median(
    x,
    na.rm = TRUE
  )

  mad_value <- mad(
    x,
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (is.na(mad_value) || mad_value == 0) {
    return(
      rep(NA_real_, length(x))
    )
  }

  (x - med) / mad_value
}


outlier_metrics <- library_metrics

outlier_metrics$rz_log10_library_size <- robust_z(
  outlier_metrics$log10_library_size
)

outlier_metrics$rz_detected_genes_gt0 <- robust_z(
  outlier_metrics$detected_genes_gt0
)

outlier_metrics$rz_zero_fraction <- robust_z(
  outlier_metrics$zero_fraction
)

outlier_metrics$flag_library_size <-
  abs(outlier_metrics$rz_log10_library_size) > 3.5

outlier_metrics$flag_detected_genes <-
  abs(outlier_metrics$rz_detected_genes_gt0) > 3.5

outlier_metrics$flag_zero_fraction <-
  abs(outlier_metrics$rz_zero_fraction) > 3.5

outlier_metrics$n_technical_flags <- rowSums(
  cbind(
    outlier_metrics$flag_library_size,
    outlier_metrics$flag_detected_genes,
    outlier_metrics$flag_zero_fraction
  ),
  na.rm = TRUE
)

write.csv(
  outlier_metrics,
  file.path(
    "metadata",
    "qc_technical_outlier_screen.csv"
  ),
  row.names = FALSE,
  na = "NA"
)


# ------------------------------------------------------------------------------
# Additional console report
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PRIMARY-COHORT TECHNICAL BALANCE\n")
cat("============================================================\n")

cat("\nBy disease group:\n")

print(
  technical_balance[
    technical_balance$grouping_variable == "diseasegroup",
    ,
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nBy lung location:\n")

print(
  technical_balance[
    technical_balance$grouping_variable == "lunglocation",
    ,
    drop = FALSE
  ],
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("ROBUST TECHNICAL-OUTLIER SCREEN\n")
cat("============================================================\n")

flagged <- outlier_metrics[
  outlier_metrics$n_technical_flags > 0,
  c(
    "sample_title",
    "donorid",
    "diseasegroup",
    "lunglocation",
    "rin",
    "library_size",
    "detected_genes_gt0",
    "zero_fraction",
    "rz_log10_library_size",
    "rz_detected_genes_gt0",
    "rz_zero_fraction",
    "n_technical_flags"
  ),
  drop = FALSE
]

if (nrow(flagged) == 0L) {

  cat("\nNo samples exceeded the robust screening threshold.\n")

} else {

  cat(
    "\nSamples with at least one robust technical flag:\n"
  )

  print(
    flagged[
      order(
        -flagged$n_technical_flags,
        flagged$sample_title
      ),
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
}

cat(
  "\nThese flags are screening signals only; ",
  "they do not trigger automatic exclusion.\n",
  sep = ""
)

cat("============================================================\n")


# ------------------------------------------------------------------------------
# QC figures
# ------------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required for QC figures. ",
    "Install it before running R/qc.R."
  )
}

library(ggplot2)

dir.create(
  file.path("figures", "qc"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 1. Metadata missingness
# ------------------------------------------------------------------------------

missing_plot_data <- missingness[
  missingness$n_missing > 0,
  ,
  drop = FALSE
]

missing_labels <- c(
  dlco_perc_corrected = "Corrected DLCO % predicted",
  dlco_perc = "DLCO % predicted",
  fvcprebd_perc = "FVC % predicted",
  severity = "Disease severity",
  smokingstatus = "Smoking status",
  age = "Age",
  lunglocation = "Lung location"
)

missing_plot_data$label <- ifelse(
  missing_plot_data$variable %in% names(missing_labels),
  missing_labels[missing_plot_data$variable],
  missing_plot_data$variable
)

missing_plot_data$label <- factor(
  missing_plot_data$label,
  levels = rev(missing_plot_data$label)
)

p_missing <- ggplot(
  missing_plot_data,
  aes(
    x = label,
    y = percent_missing
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(percent_missing, "%")
    ),
    hjust = -0.15,
    size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Metadata missingness",
    subtitle = "GSE213001 complete metadata",
    x = NULL,
    y = "Missing values (%)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    "figures",
    "qc",
    "metadata_missingness.png"
  ),
  p_missing,
  width = 8,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 2. Library size by disease group
# ------------------------------------------------------------------------------

p_library <- ggplot(
  library_metrics,
  aes(
    x = diseasegroup,
    y = library_size / 1e6
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.65
  ) +
  labs(
    title = "Library size by disease group",
    subtitle = "Raw-count library sizes before filtering and normalization",
    x = "Disease group",
    y = "Library size (millions of counts)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    "figures",
    "qc",
    "library_size_by_disease.png"
  ),
  p_library,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 3. RIN by disease group
# ------------------------------------------------------------------------------

p_rin <- ggplot(
  library_metrics,
  aes(
    x = diseasegroup,
    y = rin
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.65
  ) +
  labs(
    title = "RNA integrity by disease group",
    subtitle = "RIN is available for all 139 samples",
    x = "Disease group",
    y = "RIN"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    "figures",
    "qc",
    "rin_by_disease.png"
  ),
  p_rin,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 4. Donor-level age in primary cohort
# ------------------------------------------------------------------------------

primary_for_age <- metadata[
  metadata$diseasegroup %in% c(
    "IPF",
    "NDC"
  ),
  ,
  drop = FALSE
]

donor_age <- primary_for_age[
  !duplicated(primary_for_age$donorid),
  c(
    "donorid",
    "diseasegroup",
    "age"
  ),
  drop = FALSE
]

donor_age <- donor_age[
  !is.na(donor_age$age),
  ,
  drop = FALSE
]

p_age <- ggplot(
  donor_age,
  aes(
    x = diseasegroup,
    y = age
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.75
  ) +
  labs(
    title = "Donor age by primary disease group",
    subtitle = "One observation per donor; age is missing for one NDC donor",
    x = "Disease group",
    y = "Age (years)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    "figures",
    "qc",
    "donor_age_ipf_vs_ndc.png"
  ),
  p_age,
  width = 6,
  height = 5,
  dpi = 300
)


cat("\nQC figures written:\n")
cat("  figures/qc/metadata_missingness.png\n")
cat("  figures/qc/library_size_by_disease.png\n")
cat("  figures/qc/rin_by_disease.png\n")
cat("  figures/qc/donor_age_ipf_vs_ndc.png\n")
