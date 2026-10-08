# ============================================================
# GSE213001 — W6.17d Controlled Biological Synthesis
#
# Purpose:
#   Generate a deterministic, descriptive biological synthesis
#   from the frozen W6 biological-synthesis scaffold.
#
# IMPORTANT:
#   - no enrichment analysis is recalculated
#   - no pathway-level statistics are modified
#   - no pathway ranking is modified
#   - no theme assignment is modified
#   - no evidence classification is modified
#   - higher_in_IPF and lower_in_IPF remain separate
#   - no causal biological claims are generated
# ============================================================

input_file <-
  "results/interpretation/synthesis/W6_biological_synthesis_scaffold.csv"

output_dir <-
  "results/interpretation/synthesis"

output_file <- file.path(
  output_dir,
  "W6_controlled_biological_synthesis.csv"
)

qc_file <-
  "metadata/qc_W6_controlled_biological_synthesis.csv"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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
  "interpretive_role"
)

stopifnot(
  all(required_columns %in% names(x)),
  nrow(x) == 16L,
  identical(x$reporting_index, 1:16),
  length(unique(x$theme)) == 8L,
  length(unique(x$direction)) == 2L
)

expected_classes <- c(
  "triangulated",
  "cross_resource",
  "cross_method",
  "single_stream",
  "not_represented"
)

stopifnot(
  all(x$evidence_class %in% expected_classes),
  sum(x$evidence_class == "triangulated") == 4L,
  sum(x$evidence_class == "cross_resource") == 3L,
  sum(x$evidence_class == "cross_method") == 1L,
  sum(x$evidence_class == "single_stream") == 2L,
  sum(x$evidence_class == "not_represented") == 6L
)

direction_phrase <- function(direction) {
  if (direction == "higher_in_IPF") {
    return("higher in IPF relative to NDC")
  }

  if (direction == "lower_in_IPF") {
    return("lower in IPF relative to NDC")
  }

  stop("Unexpected direction: ", direction)
}

evidence_statement <- function(role, theme_label, direction) {

  d <- direction_phrase(direction)

  if (role == "primary_convergent_evidence") {
    return(
      paste0(
        "The ", theme_label,
        " theme shows convergent pathway-level evidence ",
        d,
        " across all three frozen evidence streams."
      )
    )
  }

  if (role == "secondary_convergent_evidence") {
    return(
      paste0(
        "The ", theme_label,
        " theme shows convergent pathway-level evidence ",
        d,
        " across two frozen evidence streams."
      )
    )
  }

  if (role == "contextual_single_stream_evidence") {
    return(
      paste0(
        "The ", theme_label,
        " theme has limited pathway-level evidence ",
        d,
        " from a single frozen evidence stream."
      )
    )
  }

  if (role == "no_controlled_pathway_support") {
    return(
      paste0(
        "The controlled pathway readout does not provide pathway-level ",
        "support for the ", theme_label,
        " theme ", d,
        ". This is not evidence of biological absence."
      )
    )
  }

  stop("Unexpected interpretive role: ", role)
}

x$biological_statement <- mapply(
  evidence_statement,
  role = x$interpretive_role,
  theme_label = x$theme_label,
  direction = x$direction,
  USE.NAMES = FALSE
)

stopifnot(
  length(x$biological_statement) == 16L,
  sum(is.na(x$biological_statement)) == 0L,
  all(nzchar(x$biological_statement))
)

out <- x[
  ,
  c(
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
]

write.csv(
  out,
  output_file,
  row.names = FALSE
)

qc <- data.frame(
  metric = c(
    "input_rows",
    "output_rows",
    "themes",
    "directions",
    "primary_convergent_evidence",
    "secondary_convergent_evidence",
    "contextual_single_stream_evidence",
    "no_controlled_pathway_support"
  ),
  value = c(
    nrow(x),
    nrow(out),
    length(unique(out$theme)),
    length(unique(out$direction)),
    sum(out$interpretive_role == "primary_convergent_evidence"),
    sum(out$interpretive_role == "secondary_convergent_evidence"),
    sum(out$interpretive_role == "contextual_single_stream_evidence"),
    sum(out$interpretive_role == "no_controlled_pathway_support")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("GSE213001 W6.17d — CONTROLLED BIOLOGICAL SYNTHESIS\n")
cat("============================================================\n")

cat("\n=== INPUT ===\n")
cat("Rows:       ", nrow(x), "\n", sep = "")
cat("Themes:     ", length(unique(x$theme)), "\n", sep = "")
cat("Directions: ", length(unique(x$direction)), "\n", sep = "")

cat("\n=== INTERPRETIVE ROLES ===\n")
print(table(out$interpretive_role))

cat("\n=== OUTPUTS ===\n")
cat(output_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No enrichment analysis was recalculated.\n")
cat("No pathway-level statistical value was modified.\n")
cat("No pathway reporting rank was modified.\n")
cat("No opposite directions were combined.\n")
cat("Statements describe pathway-level patterns, not causal mechanisms.\n")

cat("\n============================================================\n")
cat("W6.17d CONTROLLED BIOLOGICAL SYNTHESIS: COMPLETED\n")
cat("============================================================\n")
