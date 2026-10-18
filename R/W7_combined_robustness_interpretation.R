
# ============================================================
# GSE213001
# W7.4d - Combined robustness interpretation
#
# Purpose:
#   Integrate the frozen continuous-age and gender sensitivity
#   summaries into a deterministic W7 robustness report.
#
# Important:
#   - No differential-expression model is fitted.
#   - No voom transformation is recalculated.
#   - No duplicateCorrelation estimate is recalculated.
#   - No W5, W6, or previous W7 result is modified.
#   - No new statistical threshold is introduced.
#   - The frozen W5 primary analysis remains the reference.
# ============================================================


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <-
  "results/sensitivity/W7_combined_sensitivity_summary.csv"

input_qc_file <-
  "metadata/qc_W7_combined_sensitivity_synthesis.csv"

output_dir <-
  "results/sensitivity"

output_csv <-
  file.path(
    output_dir,
    "W7_combined_robustness_interpretation.csv"
  )

output_md <-
  file.path(
    output_dir,
    "W7_combined_robustness_interpretation.md"
  )

qc_file <-
  "metadata/qc_W7_combined_robustness_interpretation.csv"


# ------------------------------------------------------------
# 2. Input existence
# ------------------------------------------------------------

stopifnot(
  file.exists(input_file),
  file.exists(input_qc_file)
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Read frozen inputs
# ------------------------------------------------------------

x <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

source_qc <- read.csv(
  input_qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Required structure
# ------------------------------------------------------------

required_columns <- c(
  "sensitivity",
  "analysis_samples",
  "analysis_donors",
  "genes",
  "W5_FDR_0_05",
  "sensitivity_FDR_0_05",
  "shared_FDR",
  "W5_only_FDR",
  "sensitivity_only_FDR",
  "neither_FDR",
  "FDR_Jaccard",
  "W5_FDR_retention",
  "Pearson_logFC",
  "Spearman_logFC",
  "direction_concordance_all",
  "shared_FDR_direction_concordance",
  "median_abs_delta_logFC",
  "q95_abs_delta_logFC",
  "max_abs_delta_logFC"
)

stopifnot(
  all(required_columns %in% names(x)),
  nrow(x) == 2L,
  !anyDuplicated(x$sensitivity),
  setequal(
    x$sensitivity,
    c("continuous_age", "gender")
  )
)


# ------------------------------------------------------------
# 5. Stable ordering
# ------------------------------------------------------------

x <- x[
  match(
    c("continuous_age", "gender"),
    x$sensitivity
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    x$sensitivity,
    c("continuous_age", "gender")
  )
)


# ------------------------------------------------------------
# 6. Frozen-value integrity checks
# ------------------------------------------------------------

age <- x[x$sensitivity == "continuous_age", , drop = FALSE]
gender <- x[x$sensitivity == "gender", , drop = FALSE]

stopifnot(
  nrow(age) == 1L,
  nrow(gender) == 1L,

  age$analysis_samples == 99L,
  age$analysis_donors == 33L,

  gender$analysis_samples == 101L,
  gender$analysis_donors == 34L,

  age$genes == 15012L,
  gender$genes == 15012L,

  age$W5_FDR_0_05 == 7898L,
  gender$W5_FDR_0_05 == 7898L,

  age$sensitivity_FDR_0_05 == 5277L,
  gender$sensitivity_FDR_0_05 == 7745L,

  age$shared_FDR == 5015L,
  gender$shared_FDR == 7596L,

  age$W5_only_FDR == 2883L,
  gender$W5_only_FDR == 302L,

  age$sensitivity_only_FDR == 262L,
  gender$sensitivity_only_FDR == 149L,

  age$neither_FDR == 6852L,
  gender$neither_FDR == 6965L,

  age$shared_FDR_direction_concordance == 1,
  gender$shared_FDR_direction_concordance == 1,

  nrow(source_qc) == 12L
)


# ------------------------------------------------------------
# 7. Deterministic descriptive statements
# ------------------------------------------------------------

fmt4 <- function(z) {
  formatC(
    z,
    format = "f",
    digits = 4
  )
}

fmt_int <- function(z) {
  format(
    as.integer(z),
    big.mark = ",",
    scientific = FALSE
  )
}


age_discovery_statement <- paste0(
  "Continuous-age adjustment changed the number of genes meeting ",
  "FDR < 0.05 from ",
  fmt_int(age$W5_FDR_0_05),
  " in frozen W5 to ",
  fmt_int(age$sensitivity_FDR_0_05),
  ". ",
  fmt_int(age$shared_FDR),
  " W5 FDR-significant genes remained significant, corresponding to ",
  "W5 FDR retention of ",
  fmt4(age$W5_FDR_retention),
  "."
)

age_effect_statement <- paste0(
  "For continuous-age adjustment, genome-wide logFC agreement with W5 was ",
  "Pearson = ",
  fmt4(age$Pearson_logFC),
  " and Spearman = ",
  fmt4(age$Spearman_logFC),
  "; genome-wide directional concordance was ",
  fmt4(age$direction_concordance_all),
  ". Median absolute change in logFC was ",
  fmt4(age$median_abs_delta_logFC),
  ", the 95th percentile absolute change was ",
  fmt4(age$q95_abs_delta_logFC),
  ", and shared-FDR directional concordance was ",
  fmt4(age$shared_FDR_direction_concordance),
  "."
)

gender_discovery_statement <- paste0(
  "Gender adjustment changed the number of genes meeting ",
  "FDR < 0.05 from ",
  fmt_int(gender$W5_FDR_0_05),
  " in frozen W5 to ",
  fmt_int(gender$sensitivity_FDR_0_05),
  ". ",
  fmt_int(gender$shared_FDR),
  " W5 FDR-significant genes remained significant, corresponding to ",
  "W5 FDR retention of ",
  fmt4(gender$W5_FDR_retention),
  "."
)

gender_effect_statement <- paste0(
  "For gender adjustment, genome-wide logFC agreement with W5 was ",
  "Pearson = ",
  fmt4(gender$Pearson_logFC),
  " and Spearman = ",
  fmt4(gender$Spearman_logFC),
  "; genome-wide directional concordance was ",
  fmt4(gender$direction_concordance_all),
  ". Median absolute change in logFC was ",
  fmt4(gender$median_abs_delta_logFC),
  ", the 95th percentile absolute change was ",
  fmt4(gender$q95_abs_delta_logFC),
  ", and shared-FDR directional concordance was ",
  fmt4(gender$shared_FDR_direction_concordance),
  "."
)


# ------------------------------------------------------------
# 8. Prespecified cross-sensitivity comparison
# ------------------------------------------------------------

stopifnot(
  age$W5_FDR_retention <
    gender$W5_FDR_retention
)

cross_statement <- paste0(
  "Continuous-age adjustment produced the larger change in the ",
  "FDR-significant discovery set: W5 FDR retention was ",
  fmt4(age$W5_FDR_retention),
  " after age adjustment versus ",
  fmt4(gender$W5_FDR_retention),
  " after gender adjustment. ",
  "At the same time, shared-FDR directional concordance was 1.0000 ",
  "for both sensitivity analyses."
)

final_statement <- paste0(
  "Taken together, the prespecified W7 sensitivity analyses support ",
  "robustness of the principal frozen W5 IPF-vs-NDC differential-expression ",
  "effect pattern to the evaluated age and gender adjustments. ",
  "The FDR-significant discovery set was more sensitive to continuous-age ",
  "adjustment than to gender adjustment. This distinction does not redefine ",
  "the frozen W5 primary analysis and does not imply absence of age or gender ",
  "effects."
)


# ------------------------------------------------------------
# 9. Machine-readable interpretation table
# ------------------------------------------------------------

interpretation <- data.frame(
  sensitivity = c(
    "continuous_age",
    "gender"
  ),

  analysis_samples = c(
    age$analysis_samples,
    gender$analysis_samples
  ),

  analysis_donors = c(
    age$analysis_donors,
    gender$analysis_donors
  ),

  W5_FDR_0_05 = c(
    age$W5_FDR_0_05,
    gender$W5_FDR_0_05
  ),

  sensitivity_FDR_0_05 = c(
    age$sensitivity_FDR_0_05,
    gender$sensitivity_FDR_0_05
  ),

  shared_FDR = c(
    age$shared_FDR,
    gender$shared_FDR
  ),

  W5_FDR_retention = c(
    age$W5_FDR_retention,
    gender$W5_FDR_retention
  ),

  Pearson_logFC = c(
    age$Pearson_logFC,
    gender$Pearson_logFC
  ),

  Spearman_logFC = c(
    age$Spearman_logFC,
    gender$Spearman_logFC
  ),

  direction_concordance_all = c(
    age$direction_concordance_all,
    gender$direction_concordance_all
  ),

  median_abs_delta_logFC = c(
    age$median_abs_delta_logFC,
    gender$median_abs_delta_logFC
  ),

  q95_abs_delta_logFC = c(
    age$q95_abs_delta_logFC,
    gender$q95_abs_delta_logFC
  ),

  shared_FDR_direction_concordance = c(
    age$shared_FDR_direction_concordance,
    gender$shared_FDR_direction_concordance
  ),

  discovery_set_statement = c(
    age_discovery_statement,
    gender_discovery_statement
  ),

  effect_pattern_statement = c(
    age_effect_statement,
    gender_effect_statement
  ),

  stringsAsFactors = FALSE
)

stopifnot(
  nrow(interpretation) == 2L,
  !anyNA(interpretation)
)


# ------------------------------------------------------------
# 10. Write interpretation CSV
# ------------------------------------------------------------

write.csv(
  interpretation,
  output_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 11. Write final Markdown report
# ------------------------------------------------------------

report_lines <- c(
  "# W7 combined robustness interpretation",
  "",
  "## Scope",
  "",
  paste0(
    "This report integrates the frozen continuous-age and gender ",
    "sensitivity analyses of the frozen W5 IPF-vs-NDC ",
    "differential-expression result."
  ),
  "",
  "No new differential-expression model was fitted.",
  "",
  "The frozen W5 primary analysis remains the reference analysis.",
  "",
  "## Continuous-age sensitivity",
  "",
  age_discovery_statement,
  "",
  age_effect_statement,
  "",
  "## Gender sensitivity",
  "",
  gender_discovery_statement,
  "",
  gender_effect_statement,
  "",
  "## Cross-sensitivity interpretation",
  "",
  cross_statement,
  "",
  "## W7 robustness conclusion",
  "",
  final_statement,
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "W7 evaluates robustness to the prespecified covariate adjustments. ",
    "It does not establish that age or gender has no biological or ",
    "statistical influence."
  ),
  "",
  paste0(
    "No W5 or W6 result was recalculated, overwritten, re-ranked, ",
    "or modified during this interpretation."
  )
)

writeLines(
  report_lines,
  output_md,
  useBytes = TRUE
)


# ------------------------------------------------------------
# 12. QC metadata
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "input_rows",
    "source_qc_rows",
    "genes",
    "W5_FDR_0_05",
    "continuous_age_samples",
    "continuous_age_donors",
    "continuous_age_FDR_0_05",
    "continuous_age_shared_FDR",
    "continuous_age_W5_FDR_retention",
    "continuous_age_shared_FDR_direction_concordance",
    "gender_samples",
    "gender_donors",
    "gender_FDR_0_05",
    "gender_shared_FDR",
    "gender_W5_FDR_retention",
    "gender_shared_FDR_direction_concordance"
  ),

  value = c(
    nrow(x),
    nrow(source_qc),
    age$genes,
    age$W5_FDR_0_05,
    age$analysis_samples,
    age$analysis_donors,
    age$sensitivity_FDR_0_05,
    age$shared_FDR,
    age$W5_FDR_retention,
    age$shared_FDR_direction_concordance,
    gender$analysis_samples,
    gender$analysis_donors,
    gender$sensitivity_FDR_0_05,
    gender$shared_FDR,
    gender$W5_FDR_retention,
    gender$shared_FDR_direction_concordance
  ),

  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 13. Output integrity
# ------------------------------------------------------------

stopifnot(
  file.exists(output_csv),
  file.exists(output_md),
  file.exists(qc_file),

  nrow(interpretation) == 2L,
  nrow(qc) == 16L,

  interpretation$shared_FDR_direction_concordance[
    interpretation$sensitivity == "continuous_age"
  ] == 1,

  interpretation$shared_FDR_direction_concordance[
    interpretation$sensitivity == "gender"
  ] == 1
)


# ------------------------------------------------------------
# 14. Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W7.4d - COMBINED ROBUSTNESS INTERPRETATION\n")
cat("============================================================\n")

cat("\n===== SOURCE =====\n")
cat(input_file, "\n")
cat(input_qc_file, "\n")

cat("\n===== CONTINUOUS AGE =====\n")
cat(age_discovery_statement, "\n")
cat(age_effect_statement, "\n")

cat("\n===== GENDER =====\n")
cat(gender_discovery_statement, "\n")
cat(gender_effect_statement, "\n")

cat("\n===== CROSS-SENSITIVITY =====\n")
cat(cross_statement, "\n")

cat("\n===== FINAL INTERPRETATION =====\n")
cat(final_statement, "\n")

cat("\n===== OUTPUTS =====\n")
cat(output_csv, "\n")
cat(output_md, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No differential-expression model was fitted.\n")
cat("No voom transformation was recalculated.\n")
cat("No duplicateCorrelation estimate was recalculated.\n")
cat("No W5, W6, or previous W7 result was modified.\n")
cat("No new statistical threshold was introduced.\n")
cat("The frozen W5 primary analysis remains the reference.\n")

cat("\n============================================================\n")
cat("W7.4d COMBINED ROBUSTNESS INTERPRETATION: COMPLETED\n")
cat("============================================================\n")

