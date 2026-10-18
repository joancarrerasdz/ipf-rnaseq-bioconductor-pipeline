# ============================================================
# GSE213001
# W7.4c — Combined sensitivity synthesis
#
# Purpose:
#   Integrate the already frozen continuous-age and gender
#   W7 sensitivity summaries into one reproducible comparison
#   table.
#
# IMPORTANT:
#   - No differential-expression model is fitted here.
#   - No voom transformation is recalculated.
#   - No duplicateCorrelation estimate is recalculated.
#   - No W5, W6, or previous W7 output is modified.
#   - The frozen W5 primary model remains the reference.
# ============================================================


# ------------------------------------------------------------
# 1. Frozen inputs
# ------------------------------------------------------------

age_summary_file <- file.path(
  "results",
  "sensitivity",
  "W7_continuous_age_vs_W5_summary.csv"
)

gender_summary_file <- file.path(
  "results",
  "sensitivity",
  "W7_gender_vs_W5_summary.csv"
)

age_qc_file <- file.path(
  "metadata",
  "qc_W7_continuous_age_sensitivity.csv"
)

gender_qc_file <- file.path(
  "metadata",
  "qc_W7_gender_sensitivity.csv"
)


# ------------------------------------------------------------
# 2. Outputs
# ------------------------------------------------------------

output_file <- file.path(
  "results",
  "sensitivity",
  "W7_combined_sensitivity_summary.csv"
)

qc_file <- file.path(
  "metadata",
  "qc_W7_combined_sensitivity_synthesis.csv"
)


# ------------------------------------------------------------
# 3. Input existence
# ------------------------------------------------------------

stopifnot(
  file.exists(age_summary_file),
  file.exists(gender_summary_file),
  file.exists(age_qc_file),
  file.exists(gender_qc_file)
)


# ------------------------------------------------------------
# 4. Read frozen inputs
# ------------------------------------------------------------

age_summary <- read.csv(
  age_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gender_summary <- read.csv(
  gender_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

age_qc <- read.csv(
  age_qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gender_qc <- read.csv(
  gender_qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 5. Structural checks
# ------------------------------------------------------------

stopifnot(
  all(c("metric", "value") %in% names(age_summary)),
  all(c("metric", "value") %in% names(gender_summary)),
  all(c("metric", "value") %in% names(age_qc)),
  all(c("metric", "value") %in% names(gender_qc))
)


metric_value <- function(x, metric_name) {

  z <- x$value[x$metric == metric_name]

  stopifnot(length(z) == 1L)

  z
}


metric_numeric <- function(x, metric_name) {

  z <- suppressWarnings(
    as.numeric(metric_value(x, metric_name))
  )

  stopifnot(
    length(z) == 1L,
    is.finite(z)
  )

  z
}


# ------------------------------------------------------------
# 6. Required frozen summary metrics
# ------------------------------------------------------------

required_summary_metrics <- c(
  "genes",
  "W5_FDR_0_05",
  "W7_FDR_0_05",
  "shared_FDR",
  "W5_only_FDR",
  "W7_only_FDR",
  "neither_FDR",
  "FDR_Jaccard",
  "Pearson_logFC",
  "Spearman_logFC",
  "direction_concordance_all",
  "shared_FDR_direction_concordance",
  "median_abs_delta_logFC",
  "q95_abs_delta_logFC",
  "max_abs_delta_logFC"
)

stopifnot(
  all(required_summary_metrics %in% age_summary$metric),
  all(required_summary_metrics %in% gender_summary$metric)
)


# ------------------------------------------------------------
# 7. Frozen-value consistency checks
# ------------------------------------------------------------

stopifnot(
  metric_numeric(age_summary, "genes") == 15012,
  metric_numeric(gender_summary, "genes") == 15012,

  metric_numeric(age_summary, "W5_FDR_0_05") == 7898,
  metric_numeric(gender_summary, "W5_FDR_0_05") == 7898,

  metric_numeric(age_summary, "W7_FDR_0_05") == 5277,
  metric_numeric(gender_summary, "W7_FDR_0_05") == 7745,

  metric_numeric(age_summary, "shared_FDR") == 5015,
  metric_numeric(gender_summary, "shared_FDR") == 7596,

  metric_numeric(age_summary, "W5_only_FDR") == 2883,
  metric_numeric(gender_summary, "W5_only_FDR") == 302,

  metric_numeric(age_summary, "W7_only_FDR") == 262,
  metric_numeric(gender_summary, "W7_only_FDR") == 149,

  metric_numeric(age_summary, "neither_FDR") == 6852,
  metric_numeric(gender_summary, "neither_FDR") == 6965
)


# ------------------------------------------------------------
# 8. Cohort information from frozen QC
# ------------------------------------------------------------

age_samples <- metric_numeric(
  age_qc,
  "age_complete_samples"
)

age_donors <- metric_numeric(
  age_qc,
  "age_complete_donors"
)

gender_samples <- metric_numeric(
  gender_qc,
  "analysis_samples"
)

gender_donors <- metric_numeric(
  gender_qc,
  "analysis_donors"
)

stopifnot(
  age_samples == 99,
  age_donors == 33,
  gender_samples == 101,
  gender_donors == 34
)


# ------------------------------------------------------------
# 9. Derived descriptive metrics
# ------------------------------------------------------------

build_row <- function(
  sensitivity_name,
  summary_data,
  analysis_samples,
  analysis_donors
) {

  W5_FDR <- metric_numeric(
    summary_data,
    "W5_FDR_0_05"
  )

  W7_FDR <- metric_numeric(
    summary_data,
    "W7_FDR_0_05"
  )

  shared_FDR <- metric_numeric(
    summary_data,
    "shared_FDR"
  )

  data.frame(
    sensitivity = sensitivity_name,

    analysis_samples = analysis_samples,
    analysis_donors = analysis_donors,

    genes = metric_numeric(
      summary_data,
      "genes"
    ),

    W5_FDR_0_05 = W5_FDR,

    sensitivity_FDR_0_05 = W7_FDR,

    shared_FDR = shared_FDR,

    W5_only_FDR = metric_numeric(
      summary_data,
      "W5_only_FDR"
    ),

    sensitivity_only_FDR = metric_numeric(
      summary_data,
      "W7_only_FDR"
    ),

    neither_FDR = metric_numeric(
      summary_data,
      "neither_FDR"
    ),

    W5_FDR_retention = shared_FDR / W5_FDR,

    sensitivity_FDR_preexisting_in_W5 =
      shared_FDR / W7_FDR,

    FDR_Jaccard = metric_numeric(
      summary_data,
      "FDR_Jaccard"
    ),

    Pearson_logFC = metric_numeric(
      summary_data,
      "Pearson_logFC"
    ),

    Spearman_logFC = metric_numeric(
      summary_data,
      "Spearman_logFC"
    ),

    direction_concordance_all =
      metric_numeric(
        summary_data,
        "direction_concordance_all"
      ),

    shared_FDR_direction_concordance =
      metric_numeric(
        summary_data,
        "shared_FDR_direction_concordance"
      ),

    median_abs_delta_logFC =
      metric_numeric(
        summary_data,
        "median_abs_delta_logFC"
      ),

    q95_abs_delta_logFC =
      metric_numeric(
        summary_data,
        "q95_abs_delta_logFC"
      ),

    max_abs_delta_logFC =
      metric_numeric(
        summary_data,
        "max_abs_delta_logFC"
      ),

    stringsAsFactors = FALSE
  )
}


combined <- rbind(
  build_row(
    "continuous_age",
    age_summary,
    age_samples,
    age_donors
  ),

  build_row(
    "gender",
    gender_summary,
    gender_samples,
    gender_donors
  )
)


# ------------------------------------------------------------
# 10. Combined-table assertions
# ------------------------------------------------------------

stopifnot(
  nrow(combined) == 2L,
  identical(
    combined$sensitivity,
    c("continuous_age", "gender")
  ),

  all(combined$genes == 15012),

  all(combined$W5_FDR_0_05 == 7898),

  combined$shared_FDR[
    combined$sensitivity == "continuous_age"
  ] == 5015,

  combined$shared_FDR[
    combined$sensitivity == "gender"
  ] == 7596,

  combined$shared_FDR_direction_concordance[
    combined$sensitivity == "continuous_age"
  ] == 1,

  combined$shared_FDR_direction_concordance[
    combined$sensitivity == "gender"
  ] == 1
)


# ------------------------------------------------------------
# 11. Write combined synthesis
# ------------------------------------------------------------

write.csv(
  combined,
  output_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. QC metadata
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "input_sensitivity_analyses",
    "output_rows",
    "frozen_gene_universe",
    "frozen_W5_FDR_0_05",
    "continuous_age_samples",
    "continuous_age_donors",
    "gender_samples",
    "gender_donors",
    "continuous_age_shared_FDR",
    "gender_shared_FDR",
    "continuous_age_shared_FDR_direction_concordance",
    "gender_shared_FDR_direction_concordance"
  ),

  value = c(
    2,
    nrow(combined),
    15012,
    7898,
    age_samples,
    age_donors,
    gender_samples,
    gender_donors,
    combined$shared_FDR[
      combined$sensitivity == "continuous_age"
    ],
    combined$shared_FDR[
      combined$sensitivity == "gender"
    ],
    combined$shared_FDR_direction_concordance[
      combined$sensitivity == "continuous_age"
    ],
    combined$shared_FDR_direction_concordance[
      combined$sensitivity == "gender"
    ]
  ),

  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W7.4c — COMBINED SENSITIVITY SYNTHESIS\n")
cat("============================================================\n")

cat("\n===== INPUTS =====\n")
cat(age_summary_file, "\n")
cat(gender_summary_file, "\n")
cat(age_qc_file, "\n")
cat(gender_qc_file, "\n")

cat("\n===== COMBINED ROBUSTNESS TABLE =====\n")
print(
  combined,
  row.names = FALSE
)

cat("\n===== KEY CROSS-SENSITIVITY METRICS =====\n")

cat(
  "Age W5 FDR retention: ",
  signif(
    combined$W5_FDR_retention[
      combined$sensitivity == "continuous_age"
    ],
    7
  ),
  "\n",
  sep = ""
)

cat(
  "Gender W5 FDR retention: ",
  signif(
    combined$W5_FDR_retention[
      combined$sensitivity == "gender"
    ],
    7
  ),
  "\n",
  sep = ""
)

cat(
  "Age Pearson logFC: ",
  signif(
    combined$Pearson_logFC[
      combined$sensitivity == "continuous_age"
    ],
    7
  ),
  "\n",
  sep = ""
)

cat(
  "Gender Pearson logFC: ",
  signif(
    combined$Pearson_logFC[
      combined$sensitivity == "gender"
    ],
    7
  ),
  "\n",
  sep = ""
)

cat(
  "Age shared-FDR direction concordance: ",
  combined$shared_FDR_direction_concordance[
    combined$sensitivity == "continuous_age"
  ],
  "\n",
  sep = ""
)

cat(
  "Gender shared-FDR direction concordance: ",
  combined$shared_FDR_direction_concordance[
    combined$sensitivity == "gender"
  ],
  "\n",
  sep = ""
)

cat("\n===== OUTPUTS =====\n")
cat(output_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No differential-expression model was fitted.\n")
cat("No voom transformation was recalculated.\n")
cat("No duplicateCorrelation estimate was recalculated.\n")
cat("No W5, W6, or previous W7 result was modified.\n")
cat("The frozen W5 primary analysis remains the reference.\n")

cat("\n============================================================\n")
cat("W7.4c COMBINED SENSITIVITY SYNTHESIS: COMPLETED\n")
cat("============================================================\n")
