# ============================================================
# GSE213001 — W6.16c biological evidence classification
# ============================================================
#
# Purpose:
#   Deterministically classify each frozen theme x direction
#   combination according to the prespecified W6 biological
#   interpretation rubric.
#
# Input:
#   results/interpretation/themes/W6_thematic_evidence_summary.csv
#
# Outputs after execution:
#   results/interpretation/themes/
#     W6_biological_evidence_classification.csv
#
#   metadata/
#     qc_W6_biological_evidence_classification.csv
#
# IMPORTANT:
#   - no enrichment analysis is recalculated
#   - no pathway-level statistical value is modified
#   - no pathway is removed
#   - no pathway reporting rank is modified
#   - higher_in_IPF and lower_in_IPF remain separate
#   - classification uses only frozen evidence-stream presence
#   - no biological conclusion is generated
#
# Frozen evidence-class precedence:
#
#   1. triangulated
#   2. cross_resource
#   3. cross_method
#   4. dual_stream_other
#   5. single_stream
#   6. not_represented
#
# ============================================================


# ------------------------------------------------------------
# 1. Files
# ------------------------------------------------------------

input_file <-
  "results/interpretation/themes/W6_thematic_evidence_summary.csv"

output_file <-
  "results/interpretation/themes/W6_biological_evidence_classification.csv"

qc_file <-
  "metadata/qc_W6_biological_evidence_classification.csv"


# ------------------------------------------------------------
# 2. Read frozen W6.15 input
# ------------------------------------------------------------

x <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 3. Frozen schema
# ------------------------------------------------------------

required <- c(
  "theme",
  "direction",
  "GO_BP_ORA_n",
  "Reactome_ORA_n",
  "GO_BP_ranked_enrichment_n",
  "evidence_streams_present",
  "cross_method_support",
  "cross_resource_support",
  "triangulated_support",
  "total_theme_memberships"
)

missing <- setdiff(required, names(x))

stopifnot(
  length(missing) == 0L,
  nrow(x) == 16L,
  length(unique(x$theme)) == 8L,
  length(unique(x$direction)) == 2L
)


# ------------------------------------------------------------
# 4. Frozen theme and direction definitions
# ------------------------------------------------------------

expected_themes <- c(
  "ECM_fibrotic_remodeling",
  "cilium_axoneme_motility",
  "tissue_structural_remodeling",
  "mucosal_antimicrobial_innate_immunity",
  "leukocyte_neutrophil_processes",
  "surfactant_lipid_sterol_metabolism",
  "mitochondrial_OXPHOS_ATP",
  "broad_signalling_or_mixed"
)

expected_directions <- c(
  "higher_in_IPF",
  "lower_in_IPF"
)

stopifnot(
  setequal(unique(x$theme), expected_themes),
  setequal(unique(x$direction), expected_directions),
  !anyDuplicated(x[c("theme", "direction")])
)


# ------------------------------------------------------------
# 5. Reconstruct stream presence independently
# ------------------------------------------------------------

GO_BP_ORA_present <-
  x$GO_BP_ORA_n > 0

Reactome_ORA_present <-
  x$Reactome_ORA_n > 0

GO_BP_ranked_present <-
  x$GO_BP_ranked_enrichment_n > 0

computed_stream_count <-
  as.integer(GO_BP_ORA_present) +
  as.integer(Reactome_ORA_present) +
  as.integer(GO_BP_ranked_present)

computed_cross_method <-
  GO_BP_ORA_present &
  GO_BP_ranked_present

computed_cross_resource <-
  GO_BP_ORA_present &
  Reactome_ORA_present

computed_triangulated <-
  GO_BP_ORA_present &
  Reactome_ORA_present &
  GO_BP_ranked_present


# ------------------------------------------------------------
# 6. Validate frozen W6.15 indicators
# ------------------------------------------------------------

stopifnot(
  identical(
    computed_stream_count,
    as.integer(x$evidence_streams_present)
  ),
  identical(
    computed_cross_method,
    as.logical(x$cross_method_support)
  ),
  identical(
    computed_cross_resource,
    as.logical(x$cross_resource_support)
  ),
  identical(
    computed_triangulated,
    as.logical(x$triangulated_support)
  )
)


# ------------------------------------------------------------
# 7. Deterministic evidence classification
# ------------------------------------------------------------

evidence_class <- rep(NA_character_, nrow(x))

evidence_class[
  computed_triangulated
] <- "triangulated"

evidence_class[
  is.na(evidence_class) &
  computed_cross_resource
] <- "cross_resource"

evidence_class[
  is.na(evidence_class) &
  computed_cross_method
] <- "cross_method"

evidence_class[
  is.na(evidence_class) &
  computed_stream_count == 2L
] <- "dual_stream_other"

evidence_class[
  is.na(evidence_class) &
  computed_stream_count == 1L
] <- "single_stream"

evidence_class[
  is.na(evidence_class) &
  computed_stream_count == 0L
] <- "not_represented"

stopifnot(
  !anyNA(evidence_class)
)


# ------------------------------------------------------------
# 8. Frozen class-count expectations
# ------------------------------------------------------------

class_levels <- c(
  "triangulated",
  "cross_resource",
  "cross_method",
  "dual_stream_other",
  "single_stream",
  "not_represented"
)

class_counts <- table(
  factor(
    evidence_class,
    levels = class_levels
  )
)

expected_class_counts <- c(
  triangulated = 4L,
  cross_resource = 3L,
  cross_method = 1L,
  dual_stream_other = 0L,
  single_stream = 2L,
  not_represented = 6L
)

stopifnot(
  identical(
    as.integer(class_counts),
    as.integer(expected_class_counts)
  )
)


# ------------------------------------------------------------
# 9. Construct immutable descriptive output
# ------------------------------------------------------------

classified <- x

classified$evidence_class <-
  evidence_class

classified$GO_BP_ORA_present <-
  GO_BP_ORA_present

classified$Reactome_ORA_present <-
  Reactome_ORA_present

classified$GO_BP_ranked_enrichment_present <-
  GO_BP_ranked_present


# ------------------------------------------------------------
# 10. Final integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(classified) == 16L,
  identical(classified$theme, x$theme),
  identical(classified$direction, x$direction),
  identical(
    classified$GO_BP_ORA_n,
    x$GO_BP_ORA_n
  ),
  identical(
    classified$Reactome_ORA_n,
    x$Reactome_ORA_n
  ),
  identical(
    classified$GO_BP_ranked_enrichment_n,
    x$GO_BP_ranked_enrichment_n
  ),
  identical(
    classified$total_theme_memberships,
    x$total_theme_memberships
  )
)


# ------------------------------------------------------------
# 11. Write classification
# ------------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  classified,
  output_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. QC metadata
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "input_rows",
    "themes",
    "directions",
    "triangulated",
    "cross_resource",
    "cross_method",
    "dual_stream_other",
    "single_stream",
    "not_represented",
    "cross_method_support_rows",
    "cross_resource_support_rows",
    "triangulated_support_rows"
  ),
  value = c(
    nrow(x),
    length(unique(x$theme)),
    length(unique(x$direction)),
    unname(class_counts["triangulated"]),
    unname(class_counts["cross_resource"]),
    unname(class_counts["cross_method"]),
    unname(class_counts["dual_stream_other"]),
    unname(class_counts["single_stream"]),
    unname(class_counts["not_represented"]),
    sum(computed_cross_method),
    sum(computed_cross_resource),
    sum(computed_triangulated)
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
cat("GSE213001 W6.16c — BIOLOGICAL EVIDENCE CLASSIFICATION\n")
cat("============================================================\n")

cat("\n=== FROZEN INPUT ===\n")
cat("Theme x direction rows: ", nrow(x), "\n", sep = "")
cat("Themes:                 ", length(unique(x$theme)), "\n", sep = "")
cat("Directions:             ", length(unique(x$direction)), "\n", sep = "")

cat("\n=== EVIDENCE CLASSES ===\n")
print(class_counts)

cat("\n=== SUPPORT INDICATORS ===\n")
cat("Cross-method support:   ", sum(computed_cross_method), "\n", sep = "")
cat("Cross-resource support: ", sum(computed_cross_resource), "\n", sep = "")
cat("Triangulated support:   ", sum(computed_triangulated), "\n", sep = "")

cat("\n=== CLASSIFICATION TABLE ===\n")
print(
  classified[
    ,
    c(
      "theme",
      "direction",
      "GO_BP_ORA_n",
      "Reactome_ORA_n",
      "GO_BP_ranked_enrichment_n",
      "evidence_streams_present",
      "evidence_class"
    )
  ],
  row.names = FALSE
)

cat("\n=== OUTPUTS ===\n")
cat(output_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No enrichment analysis was recalculated.\n")
cat("No pathway-level statistical value was modified.\n")
cat("No pathway was removed.\n")
cat("No pathway reporting rank was modified.\n")
cat("Opposite directions were not combined.\n")
cat("No biological conclusion was generated.\n")

cat("\n============================================================\n")
cat("W6.16c BIOLOGICAL EVIDENCE CLASSIFICATION: COMPLETED\n")
cat("============================================================\n")
