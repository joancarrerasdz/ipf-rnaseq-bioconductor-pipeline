library(targets)

tar_option_set(
  packages = c(
    "edgeR",
    "limma",
    "ggplot2"
  )
)

run_stage <- function(script, outputs) {
  message("Running pipeline stage: ", script)

  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  missing_outputs <- outputs[!file.exists(outputs)]

  if (length(missing_outputs) > 0L) {
    stop(
      "Stage did not create expected output(s): ",
      paste(missing_outputs, collapse = ", ")
    )
  }

  outputs
}

list(

  # ------------------------------------------------------------------
  # SOURCE DATA
  # ------------------------------------------------------------------

  tar_target(
    raw_counts_source,
    "data/source/GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz",
    format = "file"
  ),

  tar_target(
    geo_series_matrix_source,
    "data/source/GSE213001_series_matrix.txt.gz",
    format = "file"
  ),

  # Tracked study metadata specification created during W2.
  tar_target(
    data_dictionary,
    "metadata/data_dictionary.csv",
    format = "file"
  ),

  # ------------------------------------------------------------------
  # ANALYSIS SCRIPTS
  # ------------------------------------------------------------------

  tar_target(
    import_script,
    "R/import.R",
    format = "file"
  ),

  tar_target(
    qc_script,
    "R/qc.R",
    format = "file"
  ),

  tar_target(
    normalization_script,
    "R/normalization.R",
    format = "file"
  ),

  # ------------------------------------------------------------------
  # W1-W2: IMPORT + METADATA ALIGNMENT
  # ------------------------------------------------------------------

  tar_target(
    import_outputs,
    {
      raw_counts_source
      geo_series_matrix_source
      data_dictionary
      import_script

      run_stage(
        import_script,
        c(
          "metadata/sample_metadata.csv",
          "metadata/donor_sample_map.csv"
        )
      )
    },
    format = "file"
  ),

  # ------------------------------------------------------------------
  # W3: METADATA + LIBRARY QC
  # ------------------------------------------------------------------

  tar_target(
    qc_outputs,
    {
      raw_counts_source
      geo_series_matrix_source
      import_outputs
      qc_script

      run_stage(
        qc_script,
        c(
          "metadata/qc_missingness.csv",
          "metadata/qc_numeric_summary.csv",
          "metadata/qc_categorical_summary.csv",
          "metadata/qc_library_metrics.csv",
          "metadata/qc_primary_technical_balance.csv",
          "metadata/qc_technical_outlier_screen.csv",
          "figures/qc/metadata_missingness.png",
          "figures/qc/library_size_by_disease.png",
          "figures/qc/rin_by_disease.png",
          "figures/qc/donor_age_ipf_vs_ndc.png"
        )
      )
    },
    format = "file"
  ),

  # ------------------------------------------------------------------
  # W4: FILTERING + TMM + ORDINATION QC
  # ------------------------------------------------------------------

  tar_target(
    normalization_outputs,
    {
      raw_counts_source
      import_outputs
      qc_outputs
      normalization_script

      run_stage(
        normalization_script,
        c(
          "metadata/qc_gene_filtering.csv",
          "metadata/qc_gene_filtering_summary.csv",
          "metadata/qc_tmm_normalization.csv",
          "metadata/qc_tmm_normalization_summary.csv",
          "metadata/qc_pca_coordinates.csv",
          "metadata/qc_mds_coordinates.csv",
          "metadata/qc_pca_variance.csv",
          "metadata/qc_ordination_donor_associations.csv",
          "metadata/qc_ordination_within_donor_associations.csv",
          "figures/qc/pca_disease_region.png",
          "figures/qc/pca_technical_watchlist.png",
          "figures/qc/pca_processing_date.png",
          "figures/qc/pca_rin.png",
          "figures/qc/mds_disease_region.png",
          "figures/qc/mds_technical_watchlist.png"
        )
      )
    },
    format = "file"
  )
)
