#!/usr/bin/env Rscript

# ============================================================
# GSE213001
# W8 robustness portfolio visualization
# ============================================================
#
# PURPOSE
# -------
# Generate a presentation-oriented summary of the two frozen W7
# sensitivity analyses:
#
# - continuous-age adjustment
# - gender adjustment
#
# The frozen W5 IPF-vs-NDC differential-expression result remains
# the principal inferential reference.
#
# IMPORTANT
# ---------
# This script:
#
# - does NOT refit the W5 primary differential-expression model,
# - does NOT refit either W7 sensitivity model,
# - does NOT recompute P values,
# - does NOT recompute adjusted P values,
# - does NOT recompute logFC estimates,
# - does NOT modify the FDR threshold,
# - does NOT introduce a new robustness threshold,
# - does NOT calculate a composite robustness score,
# - does NOT overwrite any W5-W8.3 analytical result.
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

library(ggplot2)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

source_file <- paste0(
  "results/sensitivity/",
  "W7_combined_sensitivity_summary.csv"
)

interpretation_csv <- paste0(
  "results/sensitivity/",
  "W7_combined_robustness_interpretation.csv"
)

interpretation_md <- paste0(
  "results/sensitivity/",
  "W7_combined_robustness_interpretation.md"
)

figure_dir <- "figures/final"

output_png <- file.path(
  figure_dir,
  "W8_robustness_summary.png"
)

output_pdf <- file.path(
  figure_dir,
  "W8_robustness_summary.pdf"
)

qc_file <- paste0(
  "metadata/",
  "qc_W8_robustness_visualization.csv"
)

# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------

metric_row <- function(metric, value) {

  data.frame(
    metric = metric,
    value = as.character(value),
    stringsAsFactors = FALSE
  )
}


value_for <- function(
  data,
  sensitivity_name,
  column_name
) {

  z <- data[
    data$sensitivity == sensitivity_name,
    column_name
  ]

  stopifnot(length(z) == 1L)

  z[[1L]]
}


make_metric <- function(
  data,
  column_name,
  metric_label,
  metric_group
) {

  data.frame(
    sensitivity = data$sensitivity,
    metric = metric_label,
    group = metric_group,
    value = as.numeric(
      data[[column_name]]
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 3. Frozen source checks
# ------------------------------------------------------------

stopifnot(
  file.exists(source_file),
  file.exists(interpretation_csv),
  file.exists(interpretation_md)
)

x <- read.csv(
  source_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

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
  "W5_FDR_retention",
  "sensitivity_FDR_preexisting_in_W5",
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
  identical(
    names(x),
    required_columns
  ),
  nrow(x) == 2L,
  identical(
    x$sensitivity,
    c(
      "continuous_age",
      "gender"
    )
  ),
  all(x$genes == 15012L),
  all(x$W5_FDR_0_05 == 7898L)
)

cat("\n===== FROZEN W7 INPUT =====\n")
cat("File:", source_file, "\n")
cat("Rows:", nrow(x), "\n")
cat(
  "Sensitivities:",
  paste(x$sensitivity, collapse = ", "),
  "\n"
)

# ------------------------------------------------------------
# 4. Known frozen W7 invariants
# ------------------------------------------------------------

age <- x[
  x$sensitivity == "continuous_age",
  ,
  drop = FALSE
]

gender <- x[
  x$sensitivity == "gender",
  ,
  drop = FALSE
]

stopifnot(
  nrow(age) == 1L,
  nrow(gender) == 1L,

  age$analysis_samples == 99L,
  age$analysis_donors == 33L,
  age$sensitivity_FDR_0_05 == 5277L,
  age$shared_FDR == 5015L,

  gender$analysis_samples == 101L,
  gender$analysis_donors == 34L,
  gender$sensitivity_FDR_0_05 == 7745L,
  gender$shared_FDR == 7596L,

  abs(
    age$W5_FDR_retention -
      0.634970878703469
  ) < 1e-12,

  abs(
    gender$W5_FDR_retention -
      0.961762471511775
  ) < 1e-12,

  abs(
    age$sensitivity_FDR_preexisting_in_W5 -
      0.950350577979913
  ) < 1e-12,

  abs(
    gender$sensitivity_FDR_preexisting_in_W5 -
      0.980761781794706
  ) < 1e-12,

  abs(
    age$FDR_Jaccard -
      0.614583333333333
  ) < 1e-12,

  abs(
    gender$FDR_Jaccard -
      0.943954268671555
  ) < 1e-12,

  abs(
    age$Pearson_logFC -
      0.965564642719243
  ) < 1e-12,

  abs(
    gender$Pearson_logFC -
      0.994647016447780
  ) < 1e-12,

  abs(
    age$Spearman_logFC -
      0.950892167969733
  ) < 1e-12,

  abs(
    gender$Spearman_logFC -
      0.996466528038178
  ) < 1e-12,

  abs(
    age$direction_concordance_all -
      0.912736477484679
  ) < 1e-12,

  abs(
    gender$direction_concordance_all -
      0.988542499333866
  ) < 1e-12,

  age$shared_FDR_direction_concordance == 1,
  gender$shared_FDR_direction_concordance == 1,

  abs(
    age$median_abs_delta_logFC -
      0.071034039050256
  ) < 1e-12,

  abs(
    gender$median_abs_delta_logFC -
      0.008607578494061
  ) < 1e-12,

  abs(
    age$q95_abs_delta_logFC -
      0.407189418327317
  ) < 1e-12,

  abs(
    gender$q95_abs_delta_logFC -
      0.111034080291588
  ) < 1e-12,

  abs(
    age$max_abs_delta_logFC -
      1.99744907285352
  ) < 1e-12,

  abs(
    gender$max_abs_delta_logFC -
      1.71138549393392
  ) < 1e-12
)

cat("Frozen W7 invariants: PASSED\n")

# ------------------------------------------------------------
# 5. Display data
# ------------------------------------------------------------

plot_data <- rbind(

  make_metric(
    x,
    "Pearson_logFC",
    "Pearson logFC",
    "Effect-pattern agreement"
  ),

  make_metric(
    x,
    "Spearman_logFC",
    "Spearman logFC",
    "Effect-pattern agreement"
  ),

  make_metric(
    x,
    "direction_concordance_all",
    "Genome-wide direction",
    "Effect-pattern agreement"
  ),

  make_metric(
    x,
    "shared_FDR_direction_concordance",
    "Shared-FDR direction",
    "Effect-pattern agreement"
  ),

  make_metric(
    x,
    "W5_FDR_retention",
    "W5 FDR retention",
    "Discovery-set stability"
  ),

  make_metric(
    x,
    "FDR_Jaccard",
    "FDR Jaccard",
    "Discovery-set stability"
  ),

  make_metric(
    x,
    "sensitivity_FDR_preexisting_in_W5",
    "Sensitivity FDR genes\nalready present in W5",
    "Discovery-set stability"
  ),

  make_metric(
    x,
    "median_abs_delta_logFC",
    "Median |delta logFC|",
    "Effect-size perturbation"
  ),

  make_metric(
    x,
    "q95_abs_delta_logFC",
    "95th percentile |delta logFC|",
    "Effect-size perturbation"
  ),

  make_metric(
    x,
    "max_abs_delta_logFC",
    "Maximum |delta logFC|",
    "Effect-size perturbation"
  )
)

stopifnot(
  nrow(plot_data) == 20L,
  !anyNA(plot_data$value),
  all(is.finite(plot_data$value))
)

plot_data$sensitivity_label <- ifelse(
  plot_data$sensitivity == "continuous_age",
  "Continuous age",
  "Gender"
)

plot_data$sensitivity_label <- factor(
  plot_data$sensitivity_label,
  levels = c(
    "Continuous age",
    "Gender"
  )
)

group_levels <- c(
  "Effect-pattern agreement",
  "Discovery-set stability",
  "Effect-size perturbation"
)

plot_data$group <- factor(
  plot_data$group,
  levels = group_levels
)

metric_levels <- c(
  "Pearson logFC",
  "Spearman logFC",
  "Genome-wide direction",
  "Shared-FDR direction",
  "W5 FDR retention",
  "FDR Jaccard",
  "Sensitivity FDR genes\nalready present in W5",
  "Median |delta logFC|",
  "95th percentile |delta logFC|",
  "Maximum |delta logFC|"
)

plot_data$metric <- factor(
  plot_data$metric,
  levels = rev(metric_levels)
)

plot_data$value_label <- sprintf(
  "%.3f",
  plot_data$value
)

# Opposing vertical offsets prevent near-identical values from
# obscuring each other's numeric labels.

plot_data$label_vjust <- ifelse(
  plot_data$sensitivity == "continuous_age",
  1.8,
  -0.8
)

# ------------------------------------------------------------
# 6. Dumbbell segments
# ------------------------------------------------------------

segment_split <- split(
  plot_data,
  interaction(
    plot_data$group,
    plot_data$metric,
    drop = TRUE
  )
)

segments <- do.call(
  rbind,
  lapply(
    segment_split,
    function(z) {

      data.frame(
        group = z$group[[1L]],
        metric = z$metric[[1L]],
        xmin = min(z$value),
        xmax = max(z$value),
        stringsAsFactors = FALSE
      )
    }
  )
)

segments$group <- factor(
  segments$group,
  levels = group_levels
)

segments$metric <- factor(
  segments$metric,
  levels = rev(metric_levels)
)

stopifnot(
  nrow(segments) == 10L
)

# Invisible anchors enforce scientifically interpretable panel scales:
#
# - agreement and discovery metrics: 0 to 1
# - absolute delta-logFC metrics: 0 to 2.1 logFC units
#
# These are display ranges, not statistical thresholds.

axis_anchors <- data.frame(
  group = factor(
    c(
      "Effect-pattern agreement",
      "Effect-pattern agreement",
      "Discovery-set stability",
      "Discovery-set stability",
      "Effect-size perturbation",
      "Effect-size perturbation"
    ),
    levels = group_levels
  ),

  metric = factor(
    c(
      "Pearson logFC",
      "Pearson logFC",
      "W5 FDR retention",
      "W5 FDR retention",
      "Median |delta logFC|",
      "Median |delta logFC|"
    ),
    levels = rev(metric_levels)
  ),

  value = c(
    0,
    1,
    0,
    1,
    0,
    2.1
  )
)

# ------------------------------------------------------------
# 7. Portfolio figure
# ------------------------------------------------------------

sensitivity_colours <- c(
  "Continuous age" = "#7A5195",
  "Gender" = "#2A9D8F"
)

p <- ggplot(
  plot_data,
  aes(
    x = value,
    y = metric,
    colour = sensitivity_label
  )
) +

  geom_blank(
    data = axis_anchors,
    aes(
      x = value,
      y = metric
    ),
    inherit.aes = FALSE
  ) +

  geom_segment(
    data = segments,
    aes(
      x = xmin,
      xend = xmax,
      y = metric,
      yend = metric
    ),
    inherit.aes = FALSE,
    linewidth = 1.0,
    colour = "grey82"
  ) +

  geom_point(
    size = 4.0,
    alpha = 0.92
  ) +

  geom_text(
    aes(
      label = value_label,
      vjust = label_vjust
    ),
    size = 3.25,
    fontface = "bold",
    show.legend = FALSE
  ) +

  facet_wrap(
    ~ group,
    nrow = 1,
    scales = "free"
  ) +

  scale_colour_manual(
    values = sensitivity_colours,
    drop = FALSE
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.02,
        0.06
      )
    )
  ) +

  labs(
    title = paste0(
      "Sensitivity and robustness of the frozen ",
      "IPF-vs-NDC result"
    ),

    subtitle = paste0(
      "Continuous-age and gender adjustments compared ",
      "with the frozen W5 primary analysis"
    ),

    x = "Frozen W7 metric value",
    y = NULL,

    colour = "Sensitivity adjustment",

    caption = paste0(
      "Agreement and discovery-set metrics use natural 0-1 scales ",
      "(higher = stronger agreement).\n",
      "|delta logFC| metrics are shown in logFC units ",
      "(lower = less perturbation).\n",
      "Shared-FDR direction concordance is 1.000 for both analyses.\n",
      "No composite robustness score or new statistical threshold ",
      "was introduced."
    )
  ) +

  theme_minimal(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      size = 22,
      face = "bold",
      margin = margin(
        b = 6
      )
    ),

    plot.subtitle = element_text(
      size = 14,
      margin = margin(
        b = 18
      )
    ),

    plot.caption = element_text(
      size = 9.5,
      colour = "grey35",
      hjust = 0,
      lineheight = 1.05,
      margin = margin(
        t = 16
      )
    ),

    axis.title.x = element_text(
      face = "bold",
      margin = margin(
        t = 10
      )
    ),

    axis.text.x = element_text(
      size = 10
    ),

    axis.text.y = element_text(
      size = 10.5,
      colour = "grey15"
    ),

    strip.text = element_text(
      size = 13,
      face = "bold",
      margin = margin(
        b = 8
      )
    ),

    panel.grid.major.y = element_blank(),

    panel.grid.major.x = element_line(
      linewidth = 0.4,
      colour = "grey88"
    ),

    panel.grid.minor = element_blank(),

    panel.spacing.x = grid::unit(
      1.2,
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

# ------------------------------------------------------------
# 8. Write outputs
# ------------------------------------------------------------

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(qc_file),
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = output_png,
  plot = p,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 16,
  height = 9,
  units = "in",
  device = "pdf",
  bg = "white"
)

# ------------------------------------------------------------
# 9. QC record
# ------------------------------------------------------------

qc <- do.call(
  rbind,
  list(
    metric_row(
      "source_file",
      source_file
    ),

    metric_row(
      "source_md5",
      unname(
        tools::md5sum(source_file)
      )
    ),

    metric_row(
      "source_rows",
      nrow(x)
    ),

    metric_row(
      "sensitivities",
      length(unique(x$sensitivity))
    ),

    metric_row(
      "displayed_metrics",
      length(unique(plot_data$metric))
    ),

    metric_row(
      "tested_genes",
      unique(x$genes)
    ),

    metric_row(
      "W5_FDR_0_05",
      unique(x$W5_FDR_0_05)
    ),

    metric_row(
      "continuous_age_samples",
      age$analysis_samples
    ),

    metric_row(
      "continuous_age_donors",
      age$analysis_donors
    ),

    metric_row(
      "gender_samples",
      gender$analysis_samples
    ),

    metric_row(
      "gender_donors",
      gender$analysis_donors
    ),

    metric_row(
      "continuous_age_W5_FDR_retention",
      age$W5_FDR_retention
    ),

    metric_row(
      "gender_W5_FDR_retention",
      gender$W5_FDR_retention
    ),

    metric_row(
      "continuous_age_FDR_Jaccard",
      age$FDR_Jaccard
    ),

    metric_row(
      "gender_FDR_Jaccard",
      gender$FDR_Jaccard
    ),

    metric_row(
      "continuous_age_Pearson_logFC",
      age$Pearson_logFC
    ),

    metric_row(
      "gender_Pearson_logFC",
      gender$Pearson_logFC
    ),

    metric_row(
      "continuous_age_Spearman_logFC",
      age$Spearman_logFC
    ),

    metric_row(
      "gender_Spearman_logFC",
      gender$Spearman_logFC
    ),

    metric_row(
      "continuous_age_direction_concordance",
      age$direction_concordance_all
    ),

    metric_row(
      "gender_direction_concordance",
      gender$direction_concordance_all
    ),

    metric_row(
      "continuous_age_shared_FDR_direction_concordance",
      age$shared_FDR_direction_concordance
    ),

    metric_row(
      "gender_shared_FDR_direction_concordance",
      gender$shared_FDR_direction_concordance
    ),

    metric_row(
      "continuous_age_median_abs_delta_logFC",
      age$median_abs_delta_logFC
    ),

    metric_row(
      "gender_median_abs_delta_logFC",
      gender$median_abs_delta_logFC
    ),

    metric_row(
      "continuous_age_q95_abs_delta_logFC",
      age$q95_abs_delta_logFC
    ),

    metric_row(
      "gender_q95_abs_delta_logFC",
      gender$q95_abs_delta_logFC
    ),

    metric_row(
      "continuous_age_max_abs_delta_logFC",
      age$max_abs_delta_logFC
    ),

    metric_row(
      "gender_max_abs_delta_logFC",
      gender$max_abs_delta_logFC
    ),

    metric_row(
      "figure_width_inches",
      "16"
    ),

    metric_row(
      "figure_height_inches",
      "9"
    ),

    metric_row(
      "figure_dpi",
      "300"
    ),

    metric_row(
      "composite_robustness_score_calculated",
      FALSE
    ),

    metric_row(
      "primary_DE_model_refitted",
      FALSE
    ),

    metric_row(
      "sensitivity_models_refitted",
      FALSE
    ),

    metric_row(
      "P_values_recalculated",
      FALSE
    ),

    metric_row(
      "adjusted_P_values_recalculated",
      FALSE
    ),

    metric_row(
      "FDR_threshold_modified",
      FALSE
    ),

    metric_row(
      "W7_results_modified",
      FALSE
    ),

    metric_row(
      "output_png",
      output_png
    ),

    metric_row(
      "output_pdf",
      output_pdf
    )
  )
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE,
  quote = TRUE
)

# ------------------------------------------------------------
# 10. Output verification
# ------------------------------------------------------------

stopifnot(
  file.exists(output_png),
  file.exists(output_pdf),
  file.exists(qc_file),

  file.info(output_png)$size > 0,
  file.info(output_pdf)$size > 0,

  nrow(plot_data) == 20L,
  length(unique(plot_data$metric)) == 10L,

  all(
    plot_data$value[
      plot_data$group !=
        "Effect-size perturbation"
    ] >= 0
  ),

  all(
    plot_data$value[
      plot_data$group !=
        "Effect-size perturbation"
    ] <= 1
  )
)

# ------------------------------------------------------------
# 11. Console summary
# ------------------------------------------------------------

cat("\n===== W8 ROBUSTNESS VISUALIZATION =====\n")
cat("Frozen source:", source_file, "\n")
cat("Sensitivities: 2\n")
cat("Displayed metrics: 10\n")

cat("\n===== EFFECT-PATTERN AGREEMENT =====\n")
cat(
  "Continuous-age Pearson:",
  signif(age$Pearson_logFC, 7),
  "\n"
)
cat(
  "Gender Pearson:",
  signif(gender$Pearson_logFC, 7),
  "\n"
)
cat(
  "Continuous-age direction concordance:",
  signif(age$direction_concordance_all, 7),
  "\n"
)
cat(
  "Gender direction concordance:",
  signif(gender$direction_concordance_all, 7),
  "\n"
)
cat(
  "Shared-FDR direction concordance:",
  "1 for both sensitivities\n"
)

cat("\n===== DISCOVERY-SET STABILITY =====\n")
cat(
  "Continuous-age W5 FDR retention:",
  signif(age$W5_FDR_retention, 7),
  "\n"
)
cat(
  "Gender W5 FDR retention:",
  signif(gender$W5_FDR_retention, 7),
  "\n"
)
cat(
  "Continuous-age FDR Jaccard:",
  signif(age$FDR_Jaccard, 7),
  "\n"
)
cat(
  "Gender FDR Jaccard:",
  signif(gender$FDR_Jaccard, 7),
  "\n"
)

cat("\n===== EFFECT-SIZE PERTURBATION =====\n")
cat(
  "Continuous-age median |delta logFC|:",
  signif(age$median_abs_delta_logFC, 7),
  "\n"
)
cat(
  "Gender median |delta logFC|:",
  signif(gender$median_abs_delta_logFC, 7),
  "\n"
)
cat(
  "Continuous-age q95 |delta logFC|:",
  signif(age$q95_abs_delta_logFC, 7),
  "\n"
)
cat(
  "Gender q95 |delta logFC|:",
  signif(gender$q95_abs_delta_logFC, 7),
  "\n"
)

cat("\n===== OUTPUTS =====\n")
cat(output_png, "\n")
cat(output_pdf, "\n")
cat(qc_file, "\n")

cat("\n===== ANALYTICAL BOUNDARY =====\n")
cat("Primary differential-expression model refitted: FALSE\n")
cat("Sensitivity models refitted: FALSE\n")
cat("P values recalculated: FALSE\n")
cat("Adjusted P values recalculated: FALSE\n")
cat("FDR threshold modified: FALSE\n")
cat("Composite robustness score calculated: FALSE\n")
cat("Frozen W7 results modified: FALSE\n")

cat("\n============================================================\n")
cat("W8.4 ROBUSTNESS PORTFOLIO VISUALIZATION: COMPLETED\n")
cat("============================================================\n")
