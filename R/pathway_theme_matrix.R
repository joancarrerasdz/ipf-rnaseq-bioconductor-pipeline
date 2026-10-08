# ============================================================
# GSE213001 — W6.14c Theme × direction × method evidence matrix
# ============================================================
#
# Input:
#   deterministic 120-row controlled pathway readout
#
# Purpose:
#   assign the previously frozen broad biological themes
#   without modifying pathway-level statistics or rankings.
#
# IMPORTANT:
# - no pathway is removed
# - no statistical value is recalculated
# - p-values are not used for theme assignment
# - reporting rank is not used for theme assignment
# - direction is preserved
# - method/resource are preserved
# ============================================================

input_file <-
  "results/interpretation/pathways/W6_controlled_pathway_readout_top20.csv"

output_dir <-
  "results/interpretation/themes"

assignment_file <- file.path(
  output_dir,
  "W6_pathway_theme_assignments.csv"
)

matrix_file <- file.path(
  output_dir,
  "W6_theme_direction_method_matrix.csv"
)

qc_file <-
  "metadata/qc_W6_theme_direction_method_matrix.csv"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

x <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 1. Structural integrity
# ------------------------------------------------------------

required_columns <- c(
  "resource",
  "method",
  "direction",
  "reporting_rank",
  "ID",
  "Description"
)

stopifnot(
  all(required_columns %in% names(x)),
  nrow(x) == 120L,
  sum(is.na(x$ID) | x$ID == "") == 0L,
  sum(is.na(x$Description) | x$Description == "") == 0L
)

description <- tolower(x$Description)

# ------------------------------------------------------------
# 2. Frozen biological-theme rules
# ------------------------------------------------------------

theme_patterns <- list(

  ECM_fibrotic_remodeling =
    paste(
      "extracellular matrix",
      "collagen",
      "matrix organization",
      "matrix organisation",
      "matrix assembly",
      "matrix disassembly",
      "fibro",
      "wound",
      "connective tissue",
      "integrin",
      "cell adhesion",
      sep = "|"
    ),

  cilium_axoneme_motility =
    paste(
      "cilium",
      "cilia",
      "ciliary",
      "axoneme",
      "axonemal",
      "microtubule.*motility",
      "sperm motility",
      "flagell",
      sep = "|"
    ),

  tissue_structural_remodeling =
    paste(
      "tissue development",
      "morphogenesis",
      "epitheli",
      "developmental",
      "structural",
      "cell junction",
      "cell migration",
      sep = "|"
    ),

  mucosal_antimicrobial_innate_immunity =
    paste(
      "antimicrobial",
      "defense response",
      "defence response",
      "innate immune",
      "mucosal",
      "bacter",
      "pathogen",
      "host defense",
      "host defence",
      sep = "|"
    ),

  leukocyte_neutrophil_processes =
    paste(
      "neutrophil",
      "leukocyte",
      "granulocyte",
      "myeloid",
      "chemotaxis",
      "phagocyt",
      "immune cell",
      sep = "|"
    ),

  surfactant_lipid_sterol_metabolism =
    paste(
      "surfactant",
      "lipid",
      "sterol",
      "cholesterol",
      "fatty acid",
      "phospholipid",
      "lipoprotein",
      sep = "|"
    ),

  mitochondrial_OXPHOS_ATP =
    paste(
      "mitochond",
      "oxidative phosphorylation",
      "respiratory chain",
      "electron transport",
      "cellular respiration",
      "aerobic respiration",
      "ATP synthesis",
      "ATP biosynthetic",
      "oxidation of organic",
      sep = "|"
    )

)

theme_names <- names(theme_patterns)

hits <- sapply(
  theme_patterns,
  function(pattern) grepl(pattern, description, perl = TRUE)
)

if (is.null(dim(hits))) {
  hits <- matrix(
    hits,
    ncol = length(theme_patterns),
    dimnames = list(NULL, theme_names)
  )
}

n_theme_hits <- rowSums(hits)

theme_assignment <- apply(
  hits,
  1,
  function(z) {
    matched <- theme_names[z]

    if (length(matched) == 0L) {
      return("broad_signalling_or_mixed")
    }

    if (length(matched) == 1L) {
      return(matched)
    }

    paste(matched, collapse = ";")
  }
)

assignment_status <- ifelse(
  n_theme_hits == 0L,
  "fallback_mixed",
  ifelse(
    n_theme_hits == 1L,
    "single_theme",
    "multiple_theme_hits"
  )
)

annotated <- x

annotated$theme <- theme_assignment
annotated$theme_hit_count <- n_theme_hits
annotated$theme_assignment_status <- assignment_status

# ------------------------------------------------------------
# 3. Evidence matrix
# ------------------------------------------------------------

matrix <- as.data.frame(
  xtabs(
    ~ theme + resource + method + direction,
    data = annotated
  ),
  stringsAsFactors = FALSE
)

names(matrix)[names(matrix) == "Freq"] <- "pathway_count"

matrix <- matrix[
  matrix$pathway_count > 0,
  ,
  drop = FALSE
]

matrix <- matrix[
  order(
    matrix$theme,
    matrix$resource,
    matrix$method,
    matrix$direction
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 4. QC
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "controlled_readout_rows",
    "theme_assignment_rows",
    "single_theme_rows",
    "multiple_theme_rows",
    "fallback_mixed_rows",
    "unique_resources",
    "unique_methods",
    "unique_directions",
    "matrix_nonzero_cells"
  ),
  value = c(
    nrow(x),
    nrow(annotated),
    sum(assignment_status == "single_theme"),
    sum(assignment_status == "multiple_theme_hits"),
    sum(assignment_status == "fallback_mixed"),
    length(unique(x$resource)),
    length(unique(x$method)),
    length(unique(x$direction)),
    nrow(matrix)
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(annotated) == 120L,
  identical(x$ID, annotated$ID),
  identical(x$Description, annotated$Description),
  identical(x$reporting_rank, annotated$reporting_rank)
)

# ------------------------------------------------------------
# 5. Outputs
# ------------------------------------------------------------

write.csv(
  annotated,
  assignment_file,
  row.names = FALSE
)

write.csv(
  matrix,
  matrix_file,
  row.names = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W6.14c — THEME × DIRECTION × METHOD MATRIX\n")
cat("============================================================\n")

cat("\n=== CONTROLLED INPUT ===\n")
cat("Rows:                    ", nrow(x), "\n", sep = "")
cat("Resources:               ",
    paste(sort(unique(x$resource)), collapse = ", "), "\n", sep = "")
cat("Methods:                 ",
    paste(sort(unique(x$method)), collapse = ", "), "\n", sep = "")
cat("Directions:              ",
    paste(sort(unique(x$direction)), collapse = ", "), "\n", sep = "")

cat("\n=== THEME ASSIGNMENT ===\n")
cat("Single-theme rows:       ",
    sum(assignment_status == "single_theme"), "\n", sep = "")
cat("Multiple-theme rows:     ",
    sum(assignment_status == "multiple_theme_hits"), "\n", sep = "")
cat("Fallback/mixed rows:     ",
    sum(assignment_status == "fallback_mixed"), "\n", sep = "")

cat("\n=== THEME COUNTS ===\n")
print(
  sort(
    table(annotated$theme),
    decreasing = TRUE
  )
)

cat("\n=== EVIDENCE MATRIX ===\n")
print(
  matrix,
  row.names = FALSE
)

cat("\n=== OUTPUTS ===\n")
cat(assignment_file, "\n")
cat(matrix_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("The frozen 120-row controlled pathway readout was not modified.\n")
cat("No enrichment analysis was recalculated.\n")
cat("No pathway-level statistical value was modified.\n")
cat("No pathway was removed.\n")
cat("Theme assignment did not use P values or reporting rank.\n")
cat("No biological conclusion is made by this script.\n")

cat("\n============================================================\n")
cat("W6.14c THEME × DIRECTION × METHOD MATRIX: COMPLETED\n")
cat("============================================================\n")
