# ============================================================
# GSE213001 — W6.18b Final biological report
#
# Source of truth:
#   results/interpretation/synthesis/
#   W6_controlled_biological_synthesis.csv
#
# Purpose:
#   Produce the final human-readable W6 biological report
#   from the frozen controlled biological synthesis.
#
# IMPORTANT:
#   - No differential-expression analysis is recalculated.
#   - No enrichment analysis is recalculated.
#   - No pathway-level statistical value is modified.
#   - No pathway ranking is modified.
#   - No theme assignment is modified.
#   - No evidence classification is modified.
#   - No opposite directions are combined.
#   - No new biological inference is introduced.
# ============================================================

input_file <-
  "results/interpretation/synthesis/W6_controlled_biological_synthesis.csv"

output_dir <-
  "results/interpretation/final"

report_file <- file.path(
  output_dir,
  "W6_final_biological_report.md"
)

table_file <- file.path(
  output_dir,
  "W6_final_biological_report.csv"
)

qc_file <-
  "metadata/qc_W6_final_biological_report.csv"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Read frozen controlled synthesis
# ------------------------------------------------------------

x <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "reporting_index",
  "theme",
  "theme_label",
  "direction",
  "direction_label",
  "GO_BP_ORA_n",
  "Reactome_ORA_n",
  "GO_BP_ranked_enrichment_n",
  "evidence_streams_present",
  "evidence_class",
  "evidence_priority_rank",
  "interpretive_role",
  "biological_statement"
)

stopifnot(
  all(required_columns %in% names(x)),
  nrow(x) == 16L,
  length(unique(x$theme)) == 8L,
  length(unique(x$direction)) == 2L,
  identical(x$reporting_index, 1:16),
  sum(is.na(x$biological_statement)) == 0L,
  all(nzchar(x$biological_statement))
)

# ------------------------------------------------------------
# 2. Preserve frozen ordering
# ------------------------------------------------------------

x <- x[
  order(x$reporting_index),
  ,
  drop = FALSE
]

theme_order <- unique(x$theme)

stopifnot(
  length(theme_order) == 8L
)

# ------------------------------------------------------------
# 3. Export structured final table
# ------------------------------------------------------------

final_table <- x[, required_columns, drop = FALSE]

write.csv(
  final_table,
  table_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Deterministic Markdown report
# ------------------------------------------------------------

lines <- c(
  "# GSE213001 — W6 Final Biological Report",
  "",
  "## Scope",
  "",
  paste(
    "This report presents the final controlled biological synthesis",
    "of pathway-level differences between IPF and NDC."
  ),
  "",
  paste(
    "The report is generated exclusively from the frozen W6",
    "controlled biological synthesis."
  ),
  "",
  paste(
    "No differential-expression, enrichment, pathway-ranking,",
    "theme-assignment, evidence-classification, or multiple-testing",
    "procedure is recalculated at this stage."
  ),
  "",
  paste(
    "`higher_in_IPF` and `lower_in_IPF` are reported separately",
    "throughout."
  ),
  "",
  "## Evidence hierarchy",
  "",
  "Interpretive emphasis follows the frozen hierarchy:",
  "",
  "1. `primary_convergent_evidence`",
  "2. `secondary_convergent_evidence`",
  "3. `contextual_single_stream_evidence`",
  "4. `no_controlled_pathway_support`",
  "",
  paste(
    "Absence of controlled pathway support must not be interpreted",
    "as evidence of biological absence."
  ),
  "",
  "## Biological themes",
  ""
)

for (theme_value in theme_order) {

  z <- x[x$theme == theme_value, , drop = FALSE]

  stopifnot(
    nrow(z) == 2L,
    setequal(
      z$direction,
      c("higher_in_IPF", "lower_in_IPF")
    )
  )

  theme_label <- unique(z$theme_label)

  stopifnot(
    length(theme_label) == 1L
  )

  lines <- c(
    lines,
    paste0("### ", theme_label),
    ""
  )

  for (direction_value in c("higher_in_IPF", "lower_in_IPF")) {

    r <- z[z$direction == direction_value, , drop = FALSE]

    stopifnot(
      nrow(r) == 1L
    )

    direction_heading <- if (
      direction_value == "higher_in_IPF"
    ) {
      "Higher in IPF relative to NDC"
    } else {
      "Lower in IPF relative to NDC"
    }

    lines <- c(
      lines,
      paste0("#### ", direction_heading),
      "",
      paste0(
        "**Evidence class:** `",
        r$evidence_class,
        "`"
      ),
      "",
      paste0(
        "**Interpretive role:** `",
        r$interpretive_role,
        "`"
      ),
      "",
      paste0(
        "**Controlled pathway support:** ",
        "GO:BP ORA = ",
        r$GO_BP_ORA_n,
        "; Reactome ORA = ",
        r$Reactome_ORA_n,
        "; GO:BP ranked enrichment = ",
        r$GO_BP_ranked_enrichment_n,
        "."
      ),
      "",
      r$biological_statement,
      ""
    )
  }
}

lines <- c(
  lines,
  "## Interpretation boundary",
  "",
  paste(
    "These statements describe coordinated pathway-level patterns",
    "and should not be interpreted as causal mechanisms."
  ),
  "",
  paste(
    "No pathway was promoted, removed, re-ranked, or declared",
    "biologically absent during final reporting."
  ),
  ""
)

writeLines(
  lines,
  report_file,
  useBytes = TRUE
)

# ------------------------------------------------------------
# 5. QC metadata
# ------------------------------------------------------------

role_levels <- c(
  "primary_convergent_evidence",
  "secondary_convergent_evidence",
  "contextual_single_stream_evidence",
  "no_controlled_pathway_support"
)

qc <- data.frame(
  metric = c(
    "input_rows",
    "output_table_rows",
    "themes",
    "directions",
    "reporting_indices_unique",
    "missing_biological_statements",
    "primary_convergent_evidence",
    "secondary_convergent_evidence",
    "contextual_single_stream_evidence",
    "no_controlled_pathway_support"
  ),
  value = c(
    nrow(x),
    nrow(final_table),
    length(unique(x$theme)),
    length(unique(x$direction)),
    length(unique(x$reporting_index)),
    sum(is.na(x$biological_statement) |
          x$biological_statement == ""),
    sum(x$interpretive_role == role_levels[1]),
    sum(x$interpretive_role == role_levels[2]),
    sum(x$interpretive_role == role_levels[3]),
    sum(x$interpretive_role == role_levels[4])
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Frozen expectations
# ------------------------------------------------------------

stopifnot(
  nrow(final_table) == 16L,
  length(unique(final_table$theme)) == 8L,
  length(unique(final_table$direction)) == 2L,
  length(unique(final_table$reporting_index)) == 16L,

  sum(
    final_table$interpretive_role ==
      "primary_convergent_evidence"
  ) == 4L,

  sum(
    final_table$interpretive_role ==
      "secondary_convergent_evidence"
  ) == 4L,

  sum(
    final_table$interpretive_role ==
      "contextual_single_stream_evidence"
  ) == 2L,

  sum(
    final_table$interpretive_role ==
      "no_controlled_pathway_support"
  ) == 6L
)

# ------------------------------------------------------------
# 7. Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.18b — FINAL BIOLOGICAL REPORT\n")
cat("============================================================\n")

cat("\n=== INPUT ===\n")
cat("Rows:       ", nrow(x), "\n", sep = "")
cat("Themes:     ", length(unique(x$theme)), "\n", sep = "")
cat("Directions: ", length(unique(x$direction)), "\n", sep = "")

cat("\n=== INTERPRETIVE ROLES ===\n")
print(table(x$interpretive_role))

cat("\n=== OUTPUTS ===\n")
cat(report_file, "\n")
cat(table_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No enrichment analysis was recalculated.\n")
cat("No pathway-level statistical value was modified.\n")
cat("No pathway reporting rank was modified.\n")
cat("No theme assignment was modified.\n")
cat("No evidence classification was modified.\n")
cat("No opposite directions were combined.\n")
cat("No new biological inference was introduced.\n")

cat("\n============================================================\n")
cat("W6.18b FINAL BIOLOGICAL REPORT: COMPLETED\n")
cat("============================================================\n")
