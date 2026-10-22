#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source_file <- "results/interpretation/final/W6_final_biological_report.md"
output_file <- "results/interpretation/final/W8_biological_evidence_matrix.csv"
qc_file <- "metadata/qc_W8_biological_evidence_matrix.csv"

if (!file.exists(source_file)) {
  stop("Frozen W6 biological report not found: ", source_file)
}

lines <- readLines(source_file, warn = FALSE, encoding = "UTF-8")

expected_themes <- c(
  "extracellular-matrix / collagen / fibrotic remodeling",
  "cilium / axoneme / microtubule-associated motility",
  "tissue-development / structural remodeling",
  "mucosal / antimicrobial / innate immune response",
  "leukocyte / neutrophil-associated processes",
  "surfactant / lipid / sterol metabolism",
  "mitochondrial respiration / oxidative phosphorylation / ATP production",
  "broad signalling or mixed processes requiring pathway-specific caution"
)

direction_map <- c(
  "Higher in IPF relative to NDC" = "higher_in_IPF",
  "Lower in IPF relative to NDC"  = "lower_in_IPF"
)

records <- list()
current_theme <- NA_character_
current_direction <- NA_character_
current_evidence_class <- NA_character_
current_interpretive_role <- NA_character_

for (line in lines) {

  if (grepl("^### ", line)) {
    candidate <- sub("^### ", "", line)

    if (candidate %in% expected_themes) {
      current_theme <- candidate
      current_direction <- NA_character_
      current_evidence_class <- NA_character_
      current_interpretive_role <- NA_character_
    }
  }

  if (grepl("^#### ", line) && !is.na(current_theme)) {
    candidate_direction <- sub("^#### ", "", line)

    if (candidate_direction %in% names(direction_map)) {
      current_direction <- unname(direction_map[candidate_direction])
      current_evidence_class <- NA_character_
      current_interpretive_role <- NA_character_
    }
  }

  if (
    !is.na(current_theme) &&
    !is.na(current_direction) &&
    grepl("^\\*\\*Evidence class:\\*\\*", line)
  ) {
    current_evidence_class <- sub(
      "^\\*\\*Evidence class:\\*\\* `([^`]+)`.*$",
      "\\1",
      line
    )
  }

  if (
    !is.na(current_theme) &&
    !is.na(current_direction) &&
    grepl("^\\*\\*Interpretive role:\\*\\*", line)
  ) {
    current_interpretive_role <- sub(
      "^\\*\\*Interpretive role:\\*\\* `([^`]+)`.*$",
      "\\1",
      line
    )
  }

  if (
    !is.na(current_theme) &&
    !is.na(current_direction) &&
    grepl("^\\*\\*Controlled pathway support:\\*\\*", line)
  ) {

    pattern <- paste0(
      "^\\*\\*Controlled pathway support:\\*\\* ",
      "GO:BP ORA = ([0-9]+); ",
      "Reactome ORA = ([0-9]+); ",
      "GO:BP ranked enrichment = ([0-9]+)\\.$"
    )

    m <- regexec(pattern, line)
    parts <- regmatches(line, m)[[1]]

    if (length(parts) != 4L) {
      stop(
        "Could not parse controlled pathway support for:\n",
        current_theme, " / ", current_direction, "\n",
        line
      )
    }

    if (
      is.na(current_evidence_class) ||
      is.na(current_interpretive_role)
    ) {
      stop(
        "Evidence metadata missing before pathway support for:\n",
        current_theme, " / ", current_direction
      )
    }

    counts <- as.integer(parts[2:4])

    records[[length(records) + 1L]] <- data.frame(
      theme = current_theme,
      direction = current_direction,
      evidence_class = current_evidence_class,
      interpretive_role = current_interpretive_role,
      go_bp_ora = counts[1],
      reactome_ora = counts[2],
      go_bp_ranked = counts[3],
      stringsAsFactors = FALSE
    )
  }
}

wide <- do.call(rbind, records)

if (is.null(wide) || nrow(wide) == 0L) {
  stop("No frozen W6 biological evidence records were extracted.")
}

theme_order <- setNames(seq_along(expected_themes), expected_themes)
direction_order <- c(
  higher_in_IPF = 1L,
  lower_in_IPF = 2L
)

wide$theme_order <- unname(theme_order[wide$theme])
wide$direction_order <- unname(direction_order[wide$direction])

wide <- wide[
  order(wide$theme_order, wide$direction_order),
  ,
  drop = FALSE
]

row.names(wide) <- NULL

expected_classes <- c(
  "triangulated",
  "cross_method",
  "cross_resource",
  "single_stream",
  "not_represented"
)

expected_roles <- c(
  "primary_convergent_evidence",
  "secondary_convergent_evidence",
  "contextual_single_stream_evidence",
  "no_controlled_pathway_support"
)

stopifnot(
  nrow(wide) == 16L,
  length(unique(wide$theme)) == 8L,
  setequal(unique(wide$theme), expected_themes),
  all(table(wide$theme) == 2L),
  all(table(wide$direction) == 8L),
  setequal(
    unique(wide$direction),
    c("higher_in_IPF", "lower_in_IPF")
  ),
  all(wide$evidence_class %in% expected_classes),
  all(wide$interpretive_role %in% expected_roles),
  all(wide$go_bp_ora >= 0L),
  all(wide$reactome_ora >= 0L),
  all(wide$go_bp_ranked >= 0L)
)

long <- rbind(
  data.frame(
    theme = wide$theme,
    direction = wide$direction,
    evidence_class = wide$evidence_class,
    interpretive_role = wide$interpretive_role,
    evidence_stream = "GO:BP ORA",
    support_count = wide$go_bp_ora,
    theme_order = wide$theme_order,
    direction_order = wide$direction_order,
    stream_order = 1L,
    stringsAsFactors = FALSE
  ),
  data.frame(
    theme = wide$theme,
    direction = wide$direction,
    evidence_class = wide$evidence_class,
    interpretive_role = wide$interpretive_role,
    evidence_stream = "Reactome ORA",
    support_count = wide$reactome_ora,
    theme_order = wide$theme_order,
    direction_order = wide$direction_order,
    stream_order = 2L,
    stringsAsFactors = FALSE
  ),
  data.frame(
    theme = wide$theme,
    direction = wide$direction,
    evidence_class = wide$evidence_class,
    interpretive_role = wide$interpretive_role,
    evidence_stream = "GO:BP ranked enrichment",
    support_count = wide$go_bp_ranked,
    theme_order = wide$theme_order,
    direction_order = wide$direction_order,
    stream_order = 3L,
    stringsAsFactors = FALSE
  )
)

long <- long[
  order(
    long$theme_order,
    long$direction_order,
    long$stream_order
  ),
  ,
  drop = FALSE
]

row.names(long) <- NULL

stopifnot(
  nrow(long) == 48L,
  length(unique(
    paste(
      long$theme,
      long$direction,
      long$evidence_stream,
      sep = "|||"
    )
  )) == 48L
)

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(qc_file),
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  long,
  output_file,
  row.names = FALSE,
  quote = TRUE
)

qc <- data.frame(
  metric = c(
    "source_file",
    "biological_themes",
    "directions",
    "evidence_streams",
    "theme_direction_records",
    "long_matrix_rows",
    "duplicate_theme_direction_stream_rows",
    "negative_support_counts",
    "missing_support_counts",
    "new_enrichment_analysis_performed",
    "new_differential_expression_performed",
    "frozen_W6_source_modified"
  ),
  value = c(
    source_file,
    "8",
    "2",
    "3",
    "16",
    "48",
    as.character(
      48L -
        length(unique(
          paste(
            long$theme,
            long$direction,
            long$evidence_stream,
            sep = "|||"
          )
        ))
    ),
    as.character(sum(long$support_count < 0L)),
    as.character(sum(is.na(long$support_count))),
    "FALSE",
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE,
  quote = TRUE
)

cat("\n===== W8 BIOLOGICAL EVIDENCE MATRIX =====\n")
cat("Frozen source:", source_file, "\n")
cat("Themes:", length(unique(long$theme)), "\n")
cat("Directions:", length(unique(long$direction)), "\n")
cat("Evidence streams:", length(unique(long$evidence_stream)), "\n")
cat("Theme-direction records:", nrow(wide), "\n")
cat("Long-format rows:", nrow(long), "\n")

cat("\n===== EVIDENCE CLASSES =====\n")
print(table(wide$evidence_class))

cat("\n===== INTERPRETIVE ROLES =====\n")
print(table(wide$interpretive_role))

cat("\n===== FROZEN W6 EVIDENCE =====\n")
print(
  wide[
    ,
    c(
      "theme",
      "direction",
      "evidence_class",
      "go_bp_ora",
      "reactome_ora",
      "go_bp_ranked"
    )
  ],
  row.names = FALSE
)

cat("\n===== OUTPUTS =====\n")
cat(output_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No enrichment analysis was recalculated.\n")
cat("No differential-expression analysis was recalculated.\n")
cat("No W6 evidence class or interpretive role was reclassified.\n")
cat("The matrix is a structured transcription of the frozen W6 report.\n")

cat("\n============================================================\n")
cat("W8.3e FROZEN BIOLOGICAL EVIDENCE MATRIX: COMPLETED\n")
cat("============================================================\n")
