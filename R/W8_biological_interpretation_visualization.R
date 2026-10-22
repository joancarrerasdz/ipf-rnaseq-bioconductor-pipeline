#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

library(ggplot2)

source_file <- paste0(
  "results/interpretation/final/",
  "W8_biological_evidence_matrix.csv"
)

png_file <- paste0(
  "figures/final/",
  "W8_biological_interpretation_matrix.png"
)

pdf_file <- paste0(
  "figures/final/",
  "W8_biological_interpretation_matrix.pdf"
)

qc_file <- paste0(
  "metadata/",
  "qc_W8_biological_interpretation_visualization.csv"
)

if (!file.exists(source_file)) {
  stop("Frozen biological evidence matrix not found.")
}

x <- read.csv(
  source_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "theme",
  "direction",
  "evidence_class",
  "interpretive_role",
  "evidence_stream",
  "support_count",
  "theme_order",
  "direction_order",
  "stream_order"
)

stopifnot(
  identical(names(x), required_columns),
  nrow(x) == 48L,
  length(unique(x$theme)) == 8L,
  length(unique(x$direction)) == 2L,
  length(unique(x$evidence_stream)) == 3L,
  !anyNA(x$support_count),
  all(x$support_count >= 0),
  all(x$support_count == as.integer(x$support_count))
)

x$support_count <- as.integer(x$support_count)

cell_key <- paste(
  x$theme,
  x$direction,
  x$evidence_stream,
  sep = "|||"
)

stopifnot(
  length(unique(cell_key)) == 48L
)

expected_streams <- c(
  "GO:BP ORA",
  "Reactome ORA",
  "GO:BP ranked enrichment"
)

stopifnot(
  setequal(
    unique(x$evidence_stream),
    expected_streams
  )
)

direction_totals <- aggregate(
  support_count ~ direction,
  x,
  sum
)

stream_totals <- aggregate(
  support_count ~ evidence_stream,
  x,
  sum
)

value_for <- function(
  data,
  key_column,
  key,
  value_column = "support_count"
) {

  z <- data[
    data[[key_column]] == key,
    value_column
  ]

  stopifnot(length(z) == 1L)

  z[[1]]
}

stopifnot(
  sum(x$support_count) == 121L,

  value_for(
    direction_totals,
    "direction",
    "higher_in_IPF"
  ) == 61L,

  value_for(
    direction_totals,
    "direction",
    "lower_in_IPF"
  ) == 60L,

  value_for(
    stream_totals,
    "evidence_stream",
    "GO:BP ORA"
  ) == 40L,

  value_for(
    stream_totals,
    "evidence_stream",
    "Reactome ORA"
  ) == 40L,

  value_for(
    stream_totals,
    "evidence_stream",
    "GO:BP ranked enrichment"
  ) == 41L,

  sum(x$support_count == 0L) == 26L,
  sum(x$support_count > 0L) == 22L
)

theme_lookup <- data.frame(
  theme = c(
    "extracellular-matrix / collagen / fibrotic remodeling",
    "cilium / axoneme / microtubule-associated motility",
    "tissue-development / structural remodeling",
    "mucosal / antimicrobial / innate immune response",
    "leukocyte / neutrophil-associated processes",
    "surfactant / lipid / sterol metabolism",
    paste0(
      "mitochondrial respiration / oxidative phosphorylation / ",
      "ATP production"
    ),
    paste0(
      "broad signalling or mixed processes requiring ",
      "pathway-specific caution"
    )
  ),

  theme_label = c(
    "ECM / collagen / fibrotic\nremodeling",
    "Cilium / axoneme /\nmicrotubule motility",
    "Tissue development /\nstructural remodeling",
    "Mucosal / antimicrobial /\ninnate immune response",
    "Leukocyte / neutrophil-\nassociated processes",
    "Surfactant / lipid /\nsterol metabolism",
    "Mitochondrial respiration /\nOXPHOS / ATP production",
    "Broad signalling / mixed\nprocesses (caution)"
  ),

  theme_order = 1:8,
  stringsAsFactors = FALSE
)

idx <- match(
  x$theme,
  theme_lookup$theme
)

stopifnot(
  !anyNA(idx),
  identical(
    as.integer(x$theme_order),
    theme_lookup$theme_order[idx]
  )
)

x$theme_label <- theme_lookup$theme_label[idx]

x$theme_label <- factor(
  x$theme_label,
  levels = rev(theme_lookup$theme_label)
)

x$evidence_stream <- factor(
  x$evidence_stream,
  levels = expected_streams
)

direction_lookup <- c(
  higher_in_IPF = "Higher in IPF relative to NDC",
  lower_in_IPF = "Lower in IPF relative to NDC"
)

x$direction_label <- unname(
  direction_lookup[x$direction]
)

stopifnot(
  !anyNA(x$direction_label)
)

x$direction_label <- factor(
  x$direction_label,
  levels = unname(direction_lookup)
)

zero_df <- x[
  x$support_count == 0L,
  ,
  drop = FALSE
]

support_df <- x[
  x$support_count > 0L,
  ,
  drop = FALSE
]

stopifnot(
  nrow(zero_df) == 26L,
  nrow(support_df) == 22L
)

p <- ggplot(
  x,
  aes(
    x = evidence_stream,
    y = theme_label
  )
) +

  geom_point(
    data = zero_df,
    shape = 1,
    size = 2.6,
    stroke = 0.7,
    colour = "grey60"
  ) +

  geom_point(
    data = support_df,
    aes(
      size = support_count,
      colour = direction_label
    ),
    alpha = 0.82
  ) +

  facet_grid(
    . ~ direction_label
  ) +

  scale_size_area(
    name = "Frozen pathway findings",
    max_size = 11,
    breaks = c(1, 3, 6, 9, 13),
    limits = c(1, 13)
  ) +

  scale_colour_manual(
    values = c(
      "Higher in IPF relative to NDC" = "#D95F02",
      "Lower in IPF relative to NDC" = "#0072B2"
    ),
    guide = "none"
  ) +

  labs(
    title = paste0(
      "Biological interpretation - ",
      "convergent pathway evidence"
    ),
    subtitle = paste0(
      "Frozen W6 pathway support across ",
      "three evidence streams"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Bubble area denotes the number of frozen pathway-level ",
      "findings supporting each theme. ",
      "Open circles denote zero support. ",
      "Absence of pathway-level support is not evidence of ",
      "biological absence."
    )
  ) +

  theme_minimal(
    base_size = 14
  ) +

  theme(
    plot.title = element_text(
      size = 23,
      face = "bold",
      margin = margin(b = 6)
    ),

    plot.subtitle = element_text(
      size = 15,
      margin = margin(b = 18)
    ),

    plot.caption = element_text(
      size = 10.5,
      hjust = 0,
      colour = "grey35",
      margin = margin(t = 18)
    ),

    axis.text.x = element_text(
      size = 12,
      face = "bold",
      margin = margin(t = 8)
    ),

    axis.text.y = element_text(
      size = 11.5,
      colour = "grey15",
      lineheight = 1.0
    ),

    strip.text = element_text(
      size = 14,
      face = "bold"
    ),

    panel.grid.major.x = element_line(
      linewidth = 0.45,
      colour = "grey88"
    ),

    panel.grid.major.y = element_line(
      linewidth = 0.45,
      colour = "grey90"
    ),

    panel.grid.minor = element_blank(),

    panel.spacing.x = grid::unit(
      1.4,
      "lines"
    ),

    legend.position = "bottom",

    legend.title = element_text(
      face = "bold"
    ),

    plot.margin = margin(
      20,
      25,
      20,
      25
    )
  )

dir.create(
  dirname(png_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(qc_file),
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = png_file,
  plot = p,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 16,
  height = 9,
  units = "in",
  device = "pdf",
  bg = "white"
)

qc <- data.frame(
  metric = c(
    "source_file",
    "source_md5",
    "matrix_rows",
    "biological_themes",
    "directions",
    "evidence_streams",
    "zero_support_cells",
    "nonzero_support_cells",
    "total_support_count",
    "higher_in_IPF_support_count",
    "lower_in_IPF_support_count",
    "GO_BP_ORA_support_count",
    "Reactome_ORA_support_count",
    "GO_BP_ranked_support_count",
    "figure_width_inches",
    "figure_height_inches",
    "figure_dpi",
    "enrichment_recalculated",
    "differential_expression_recalculated",
    "gene_ranking_modified",
    "W6_interpretation_modified"
  ),

  value = c(
    source_file,
    unname(tools::md5sum(source_file)),
    "48",
    "8",
    "2",
    "3",
    "26",
    "22",
    "121",
    "61",
    "60",
    "40",
    "40",
    "41",
    "16",
    "9",
    "300",
    "FALSE",
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

cat("\n===== W8 BIOLOGICAL INTERPRETATION FIGURE =====\n")
cat("Input matrix:", source_file, "\n")
cat("Matrix rows:", nrow(x), "\n")
cat("Themes:", length(unique(x$theme)), "\n")
cat("Directions:", length(unique(x$direction)), "\n")
cat(
  "Evidence streams:",
  length(unique(x$evidence_stream)),
  "\n"
)
cat("Zero-support cells:", nrow(zero_df), "\n")
cat("Non-zero cells:", nrow(support_df), "\n")
cat("Frozen support total:", sum(x$support_count), "\n")

cat("\n===== OUTPUTS =====\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat(qc_file, "\n")

cat("\n===== INTERPRETIVE BOUNDARY =====\n")
cat("No enrichment analysis was recalculated.\n")
cat("No differential-expression analysis was recalculated.\n")
cat("No gene ranking was modified.\n")
cat("No W6 biological interpretation was modified.\n")
cat(
  "Bubble size represents pathway-support count, ",
  "not biological effect magnitude.\n",
  sep = ""
)
cat(
  "Zero support does not imply biological absence.\n"
)

cat("\n============================================================\n")
cat("W8.3 BIOLOGICAL INTERPRETATION VISUALIZATION: COMPLETED\n")
cat("============================================================\n")
