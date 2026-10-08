# ============================================================
# GSE213001 — W6.17b Biological synthesis scaffold
#
# Purpose:
#   Construct a deterministic, post-inference biological
#   synthesis scaffold from the frozen W6 evidence
#   classification.
#
# IMPORTANT:
#   - no enrichment analysis is recalculated
#   - no pathway-level statistic is modified
#   - no pathway ranking is modified
#   - no biological conclusion is generated
# ============================================================

input_file <- paste0(
  "results/interpretation/themes/",
  "W6_biological_evidence_classification.csv"
)

output_dir <- "results/interpretation/synthesis"

output_file <- file.path(
  output_dir,
  "W6_biological_synthesis_scaffold.csv"
)

qc_file <- "metadata/qc_W6_biological_synthesis_scaffold.csv"

stopifnot(file.exists(input_file))

x <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required <- c(
  "theme",
  "direction",
  "GO_BP_ORA_n",
  "Reactome_ORA_n",
  "GO_BP_ranked_enrichment_n",
  "evidence_streams_present",
  "evidence_class"
)

stopifnot(all(required %in% names(x)))

stopifnot(
  nrow(x) == 16L,
  length(unique(x$theme)) == 8L,
  length(unique(x$direction)) == 2L
)

# ------------------------------------------------------------
# Frozen biological-theme labels
# ------------------------------------------------------------

theme_labels <- c(
  ECM_fibrotic_remodeling =
    "extracellular-matrix / collagen / fibrotic remodeling",

  cilium_axoneme_motility =
    "cilium / axoneme / microtubule-associated motility",

  tissue_structural_remodeling =
    "tissue-development / structural remodeling",

  mucosal_antimicrobial_innate_immunity =
    "mucosal / antimicrobial / innate immune response",

  leukocyte_neutrophil_processes =
    "leukocyte / neutrophil-associated processes",

  surfactant_lipid_sterol_metabolism =
    "surfactant / lipid / sterol metabolism",

  mitochondrial_OXPHOS_ATP =
    "mitochondrial respiration / oxidative phosphorylation / ATP production",

  broad_signalling_or_mixed =
    "broad signalling or mixed processes requiring pathway-specific caution"
)

stopifnot(all(x$theme %in% names(theme_labels)))

x$theme_label <- unname(theme_labels[x$theme])

# ------------------------------------------------------------
# Frozen evidence hierarchy
# ------------------------------------------------------------

priority_rank <- c(
  triangulated = 1L,
  cross_resource = 2L,
  cross_method = 2L,
  dual_stream_other = 3L,
  single_stream = 4L,
  not_represented = 5L
)

interpretive_role <- c(
  triangulated = "primary_convergent_evidence",
  cross_resource = "secondary_convergent_evidence",
  cross_method = "secondary_convergent_evidence",
  dual_stream_other = "secondary_multi_stream_evidence",
  single_stream = "contextual_single_stream_evidence",
  not_represented = "no_controlled_pathway_support"
)

stopifnot(all(x$evidence_class %in% names(priority_rank)))

x$evidence_priority_rank <-
  unname(priority_rank[x$evidence_class])

x$interpretive_role <-
  unname(interpretive_role[x$evidence_class])

# ------------------------------------------------------------
# Explicit direction label
# ------------------------------------------------------------

direction_label <- c(
  higher_in_IPF = "higher in IPF relative to NDC",
  lower_in_IPF = "lower in IPF relative to NDC"
)

stopifnot(all(x$direction %in% names(direction_label)))

x$direction_label <-
  unname(direction_label[x$direction])

# ------------------------------------------------------------
# Deterministic reporting order
# ------------------------------------------------------------

theme_order <- c(
  "ECM_fibrotic_remodeling",
  "cilium_axoneme_motility",
  "tissue_structural_remodeling",
  "mucosal_antimicrobial_innate_immunity",
  "leukocyte_neutrophil_processes",
  "surfactant_lipid_sterol_metabolism",
  "mitochondrial_OXPHOS_ATP",
  "broad_signalling_or_mixed"
)

direction_order <- c(
  "higher_in_IPF",
  "lower_in_IPF"
)

x$theme_order <- match(x$theme, theme_order)
x$direction_order <- match(x$direction, direction_order)

stopifnot(
  !anyNA(x$theme_order),
  !anyNA(x$direction_order)
)

x <- x[
  order(
    x$theme_order,
    x$direction_order
  ),
  ,
  drop = FALSE
]

x$reporting_index <- seq_len(nrow(x))

# ------------------------------------------------------------
# Final scaffold
# ------------------------------------------------------------

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
    "interpretive_role"
  ),
  drop = FALSE
]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  out,
  output_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

class_counts <- table(out$evidence_class)

qc <- data.frame(
  metric = c(
    "input_rows",
    "output_rows",
    "themes",
    "directions",
    "triangulated",
    "cross_resource",
    "cross_method",
    "dual_stream_other",
    "single_stream",
    "not_represented"
  ),
  value = c(
    nrow(x),
    nrow(out),
    length(unique(out$theme)),
    length(unique(out$direction)),
    sum(out$evidence_class == "triangulated"),
    sum(out$evidence_class == "cross_resource"),
    sum(out$evidence_class == "cross_method"),
    sum(out$evidence_class == "dual_stream_other"),
    sum(out$evidence_class == "single_stream"),
    sum(out$evidence_class == "not_represented")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# Frozen expectations
# ------------------------------------------------------------

stopifnot(
  nrow(out) == 16L,
  length(unique(out$theme)) == 8L,
  length(unique(out$direction)) == 2L,

  sum(out$evidence_class == "triangulated") == 4L,
  sum(out$evidence_class == "cross_resource") == 3L,
  sum(out$evidence_class == "cross_method") == 1L,
  sum(out$evidence_class == "dual_stream_other") == 0L,
  sum(out$evidence_class == "single_stream") == 2L,
  sum(out$evidence_class == "not_represented") == 6L,

  identical(out$reporting_index, seq_len(16L))
)

# ------------------------------------------------------------
# Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.17b — BIOLOGICAL SYNTHESIS SCAFFOLD\n")
cat("============================================================\n")

cat("\n=== INPUT ===\n")
cat("Theme × direction rows: ", nrow(x), "\n", sep = "")
cat("Themes:                 ", length(unique(x$theme)), "\n", sep = "")
cat("Directions:             ", length(unique(x$direction)), "\n", sep = "")

cat("\n=== EVIDENCE CLASS COUNTS ===\n")
print(table(out$evidence_class))

cat("\n=== SYNTHESIS SCAFFOLD ===\n")

print(
  out[
    ,
    c(
      "reporting_index",
      "theme",
      "direction",
      "evidence_streams_present",
      "evidence_class",
      "interpretive_role"
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
cat("No pathway reporting rank was modified.\n")
cat("No opposite directions were combined.\n")
cat("No biological conclusion was generated.\n")

cat("\n============================================================\n")
cat("W6.17b BIOLOGICAL SYNTHESIS SCAFFOLD: COMPLETED\n")
cat("============================================================\n")
