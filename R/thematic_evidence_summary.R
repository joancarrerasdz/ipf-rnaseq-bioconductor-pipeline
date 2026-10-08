# ============================================================
# GSE213001 — W6.15c Deterministic thematic evidence summary
# ============================================================
#
# Purpose:
#   Aggregate the frozen W6.14 theme assignments into the
#   complete prespecified theme × direction evidence grid.
#
# IMPORTANT:
#   - no enrichment analysis is recalculated
#   - no pathway-level statistical value is modified
#   - no pathway is removed
#   - no reporting rank is modified
#   - no biological conclusion is generated
#
# Input:
#   results/interpretation/themes/W6_pathway_theme_assignments.csv
#
# Outputs:
#   results/interpretation/themes/W6_thematic_evidence_summary.csv
#   metadata/qc_W6_thematic_evidence_summary.csv
# ============================================================

input_file <-
  "results/interpretation/themes/W6_pathway_theme_assignments.csv"

output_dir <-
  "results/interpretation/themes"

output_file <- file.path(
  output_dir,
  "W6_thematic_evidence_summary.csv"
)

qc_file <-
  "metadata/qc_W6_thematic_evidence_summary.csv"

stopifnot(file.exists(input_file))

x <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "resource",
  "method",
  "direction",
  "ID",
  "theme"
)

stopifnot(all(required_cols %in% names(x)))

# ------------------------------------------------------------
# Frozen theme and direction order
# ------------------------------------------------------------

themes <- c(
  "ECM_fibrotic_remodeling",
  "cilium_axoneme_motility",
  "tissue_structural_remodeling",
  "mucosal_antimicrobial_innate_immunity",
  "leukocyte_neutrophil_processes",
  "surfactant_lipid_sterol_metabolism",
  "mitochondrial_OXPHOS_ATP",
  "broad_signalling_or_mixed"
)

directions <- c(
  "higher_in_IPF",
  "lower_in_IPF"
)

# ------------------------------------------------------------
# Expand multi-theme rows
# ------------------------------------------------------------

split_theme <- strsplit(
  x$theme,
  ";",
  fixed = TRUE
)

expanded <- do.call(
  rbind,
  lapply(seq_len(nrow(x)), function(i) {

    th <- trimws(split_theme[[i]])
    th <- th[nzchar(th)]

    if (length(th) == 0L) {
      th <- "broad_signalling_or_mixed"
    }

    data.frame(
      resource = rep(x$resource[i], length(th)),
      method = rep(x$method[i], length(th)),
      direction = rep(x$direction[i], length(th)),
      ID = rep(x$ID[i], length(th)),
      theme = th,
      stringsAsFactors = FALSE
    )
  })
)

stopifnot(
  all(expanded$theme %in% themes),
  all(expanded$direction %in% directions)
)

# A pathway may contribute at most once within each
# theme × direction × evidence-stream combination.

expanded <- unique(
  expanded[
    c(
      "theme",
      "direction",
      "resource",
      "method",
      "ID"
    )
  ]
)

# ------------------------------------------------------------
# Frozen evidence streams
# ------------------------------------------------------------

expanded$evidence_stream <- NA_character_

expanded$evidence_stream[
  expanded$resource == "GO:BP" &
    expanded$method == "ORA"
] <- "GO_BP_ORA"

expanded$evidence_stream[
  expanded$resource == "Reactome" &
    expanded$method == "ORA"
] <- "Reactome_ORA"

expanded$evidence_stream[
  expanded$resource == "GO:BP" &
    expanded$method == "ranked_enrichment"
] <- "GO_BP_ranked_enrichment"

stopifnot(
  sum(is.na(expanded$evidence_stream)) == 0L
)

# ------------------------------------------------------------
# Complete 8 × 2 reporting grid
# ------------------------------------------------------------

grid <- expand.grid(
  theme = themes,
  direction = directions,
  stringsAsFactors = FALSE
)

grid$theme <- factor(
  grid$theme,
  levels = themes
)

grid$direction <- factor(
  grid$direction,
  levels = directions
)

grid <- grid[
  order(
    grid$theme,
    grid$direction
  ),
]

grid$theme <- as.character(grid$theme)
grid$direction <- as.character(grid$direction)

# ------------------------------------------------------------
# Count evidence-stream support
# ------------------------------------------------------------

count_stream <- function(theme_value, direction_value, stream_value) {

  z <- expanded[
    expanded$theme == theme_value &
      expanded$direction == direction_value &
      expanded$evidence_stream == stream_value,
    ,
    drop = FALSE
  ]

  length(unique(z$ID))
}

grid$GO_BP_ORA_n <- mapply(
  count_stream,
  grid$theme,
  grid$direction,
  MoreArgs = list(
    stream_value = "GO_BP_ORA"
  )
)

grid$Reactome_ORA_n <- mapply(
  count_stream,
  grid$theme,
  grid$direction,
  MoreArgs = list(
    stream_value = "Reactome_ORA"
  )
)

grid$GO_BP_ranked_enrichment_n <- mapply(
  count_stream,
  grid$theme,
  grid$direction,
  MoreArgs = list(
    stream_value = "GO_BP_ranked_enrichment"
  )
)

# ------------------------------------------------------------
# Deterministic concordance indicators
# ------------------------------------------------------------

grid$evidence_streams_present <-
  (grid$GO_BP_ORA_n > 0) +
  (grid$Reactome_ORA_n > 0) +
  (grid$GO_BP_ranked_enrichment_n > 0)

grid$cross_method_support <-
  grid$GO_BP_ORA_n > 0 &
  grid$GO_BP_ranked_enrichment_n > 0

grid$cross_resource_support <-
  grid$GO_BP_ORA_n > 0 &
  grid$Reactome_ORA_n > 0

grid$triangulated_support <-
  grid$GO_BP_ORA_n > 0 &
  grid$Reactome_ORA_n > 0 &
  grid$GO_BP_ranked_enrichment_n > 0

grid$total_theme_memberships <-
  grid$GO_BP_ORA_n +
  grid$Reactome_ORA_n +
  grid$GO_BP_ranked_enrichment_n

# ------------------------------------------------------------
# Frozen expectations
# ------------------------------------------------------------

stopifnot(
  nrow(grid) == 16L,
  length(unique(grid$theme)) == 8L,
  length(unique(grid$direction)) == 2L,
  all(grid$evidence_streams_present >= 0L),
  all(grid$evidence_streams_present <= 3L),
  all(grid$total_theme_memberships >= 0L)
)

# ------------------------------------------------------------
# Write result
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  grid,
  output_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "input_pathway_rows",
    "expanded_theme_memberships",
    "summary_rows",
    "themes",
    "directions",
    "rows_with_GO_BP_ORA",
    "rows_with_Reactome_ORA",
    "rows_with_GO_BP_ranked_enrichment",
    "cross_method_support_rows",
    "cross_resource_support_rows",
    "triangulated_support_rows"
  ),
  value = c(
    nrow(x),
    nrow(expanded),
    nrow(grid),
    length(unique(grid$theme)),
    length(unique(grid$direction)),
    sum(grid$GO_BP_ORA_n > 0),
    sum(grid$Reactome_ORA_n > 0),
    sum(grid$GO_BP_ranked_enrichment_n > 0),
    sum(grid$cross_method_support),
    sum(grid$cross_resource_support),
    sum(grid$triangulated_support)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.15c — DETERMINISTIC THEMATIC EVIDENCE SUMMARY\n")
cat("============================================================\n")

cat("\n=== INPUT ===\n")
cat("Controlled pathway rows: ", nrow(x), "\n", sep = "")
cat("Expanded memberships:    ", nrow(expanded), "\n", sep = "")

cat("\n=== REPORTING GRID ===\n")
cat("Themes:                   ", length(themes), "\n", sep = "")
cat("Directions:               ", length(directions), "\n", sep = "")
cat("Summary rows:             ", nrow(grid), "\n", sep = "")

cat("\n=== EVIDENCE COVERAGE ===\n")
cat(
  "GO:BP ORA combinations:             ",
  sum(grid$GO_BP_ORA_n > 0),
  "\n",
  sep = ""
)
cat(
  "Reactome ORA combinations:          ",
  sum(grid$Reactome_ORA_n > 0),
  "\n",
  sep = ""
)
cat(
  "GO:BP ranked enrichment combinations: ",
  sum(grid$GO_BP_ranked_enrichment_n > 0),
  "\n",
  sep = ""
)

cat("\n=== CONCORDANCE ===\n")
cat(
  "Cross-method support:    ",
  sum(grid$cross_method_support),
  "\n",
  sep = ""
)
cat(
  "Cross-resource support:  ",
  sum(grid$cross_resource_support),
  "\n",
  sep = ""
)
cat(
  "Triangulated support:     ",
  sum(grid$triangulated_support),
  "\n",
  sep = ""
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
cat("W6.15c THEMATIC EVIDENCE SUMMARY: COMPLETED\n")
cat("============================================================\n")
