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


run_quarto_report <- function(source, outputs) {

  message("Rendering final portfolio report: ", source)

  quarto <- Sys.which("quarto")

  if (!nzchar(quarto)) {
    stop("Quarto executable not found.")
  }

  status <- system2(
    quarto,
    c(
      "render",
      source
    )
  )

  if (!identical(status, 0L)) {
    stop(
      "Quarto rendering failed with status ",
      status
    )
  }

  missing_outputs <- outputs[!file.exists(outputs)]

  if (length(missing_outputs) > 0L) {
    stop(
      "Quarto did not create expected output(s): ",
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

  tar_target(
    differential_expression_script,
    "R/differential_expression.R",
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
  ),

  # ------------------------------------------------------------------
  # W5: FROZEN PRIMARY DONOR-AWARE DIFFERENTIAL EXPRESSION
  # ------------------------------------------------------------------

  tar_target(
    w5_primary_de_outputs,
    {
      raw_counts_source
      import_outputs
      normalization_outputs
      differential_expression_script

      run_stage(
        differential_expression_script,
        c(
          "results/differential_expression/primary_IPF_vs_NDC_all_genes.csv",
          "results/differential_expression/primary_IPF_vs_NDC_summary.csv"
        )
      )
    },
    format = "file"
  ),

  # ------------------------------------------------------------------
  # W6: FROZEN ANNOTATION, ENRICHMENT AND BIOLOGICAL INTERPRETATION
  # ------------------------------------------------------------------

  # W6 script provenance ------------------------------------------------

  tar_target(
    annotation_mapping_script,
    "R/annotation_mapping.R",
    format = "file"
  ),

  tar_target(
    annotation_finalize_script,
    "R/annotation_finalize.R",
    format = "file"
  ),

  tar_target(
    enrichment_inputs_script,
    "R/enrichment_inputs.R",
    format = "file"
  ),

  tar_target(
    enrichment_go_bp_ora_higher_script,
    "R/enrichment_GO_BP_ORA.R",
    format = "file"
  ),

  tar_target(
    enrichment_go_bp_ora_lower_script,
    "R/enrichment_GO_BP_ORA_lower.R",
    format = "file"
  ),

  tar_target(
    enrichment_reactome_ora_higher_script,
    "R/enrichment_Reactome_ORA.R",
    format = "file"
  ),

  tar_target(
    enrichment_reactome_ora_lower_script,
    "R/enrichment_Reactome_ORA_lower.R",
    format = "file"
  ),

  tar_target(
    enrichment_go_bp_gsea_script,
    "R/enrichment_GO_BP_GSEA.R",
    format = "file"
  ),

  tar_target(
    enrichment_go_bp_gsea_qc_script,
    "R/enrichment_GO_BP_GSEA_QC.R",
    format = "file"
  ),

  tar_target(
    controlled_pathway_readout_script,
    "R/controlled_pathway_readout.R",
    format = "file"
  ),

  tar_target(
    pathway_theme_matrix_script,
    "R/pathway_theme_matrix.R",
    format = "file"
  ),

  tar_target(
    thematic_evidence_summary_script,
    "R/thematic_evidence_summary.R",
    format = "file"
  ),

  tar_target(
    biological_evidence_classification_script,
    "R/biological_evidence_classification.R",
    format = "file"
  ),

  tar_target(
    biological_synthesis_scaffold_script,
    "R/biological_synthesis_scaffold.R",
    format = "file"
  ),

  tar_target(
    controlled_biological_synthesis_script,
    "R/controlled_biological_synthesis.R",
    format = "file"
  ),

  tar_target(
    final_biological_report_script,
    "R/final_biological_report.R",
    format = "file"
  ),

  # W6.3 annotation mapping ---------------------------------------------

  tar_target(
    w6_annotation_mapping_outputs,
    {
      w5_primary_de_outputs
      annotation_mapping_script

      run_stage(
        annotation_mapping_script,
        c(
          "metadata/qc_annotation_raw_mapping.csv",
          "metadata/qc_annotation_gene_mapping.csv",
          "metadata/qc_annotation_shared_symbols.csv",
          "metadata/qc_annotation_summary.csv"
        )
      )
    },
    format = "file"
  ),

  # W6.4 frozen annotated DE table --------------------------------------

  tar_target(
    w6_annotated_de_outputs,
    {
      w5_primary_de_outputs
      w6_annotation_mapping_outputs
      annotation_finalize_script

      run_stage(
        annotation_finalize_script,
        c(
          "results/differential_expression/primary_IPF_vs_NDC_annotated.csv",
          "metadata/qc_annotation_join_audit.csv"
        )
      )
    },
    format = "file"
  ),

  # W6.8 enrichment input construction ---------------------------------

  tar_target(
    w6_enrichment_input_outputs,
    {
      w6_annotated_de_outputs
      enrichment_inputs_script

      run_stage(
        enrichment_inputs_script,
        c(
          "results/enrichment/input/W6_enrichment_background.csv",
          "results/enrichment/input/W6_ORA_higher_in_IPF.csv",
          "results/enrichment/input/W6_ORA_lower_in_IPF.csv",
          "results/enrichment/input/W6_ranked_enrichment_moderated_t.csv",
          "metadata/qc_W6_enrichment_inputs.csv"
        )
      )
    },
    format = "file"
  ),

  # W6 functional enrichment -------------------------------------------

  tar_target(
    w6_go_bp_ora_higher_outputs,
    {
      w6_enrichment_input_outputs
      enrichment_go_bp_ora_higher_script

      run_stage(
        enrichment_go_bp_ora_higher_script,
        c(
          "results/enrichment/GO_BP/W6_GO_BP_ORA_higher_in_IPF_all_terms.csv",
          "results/enrichment/GO_BP/W6_GO_BP_ORA_higher_in_IPF_significant.csv",
          "metadata/qc_W6_GO_BP_ORA_higher_in_IPF.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_go_bp_ora_lower_outputs,
    {
      w6_enrichment_input_outputs
      enrichment_go_bp_ora_lower_script

      run_stage(
        enrichment_go_bp_ora_lower_script,
        c(
          "results/enrichment/GO_BP/W6_GO_BP_ORA_lower_in_IPF_all_terms.csv",
          "results/enrichment/GO_BP/W6_GO_BP_ORA_lower_in_IPF_significant.csv",
          "metadata/qc_W6_GO_BP_ORA_lower_in_IPF.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_reactome_ora_higher_outputs,
    {
      w6_enrichment_input_outputs
      enrichment_reactome_ora_higher_script

      run_stage(
        enrichment_reactome_ora_higher_script,
        c(
          "results/enrichment/Reactome/W6_Reactome_ORA_higher_in_IPF_all_terms.csv",
          "results/enrichment/Reactome/W6_Reactome_ORA_higher_in_IPF_significant.csv",
          "metadata/qc_W6_Reactome_ORA_higher_in_IPF.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_reactome_ora_lower_outputs,
    {
      w6_enrichment_input_outputs
      enrichment_reactome_ora_lower_script

      run_stage(
        enrichment_reactome_ora_lower_script,
        c(
          "results/enrichment/Reactome/W6_Reactome_ORA_lower_in_IPF_all_terms.csv",
          "results/enrichment/Reactome/W6_Reactome_ORA_lower_in_IPF_significant.csv",
          "metadata/qc_W6_Reactome_ORA_lower_in_IPF.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_go_bp_gsea_outputs,
    {
      w6_enrichment_input_outputs
      enrichment_go_bp_gsea_script

      run_stage(
        enrichment_go_bp_gsea_script,
        c(
          "results/enrichment/GO_BP/W6_GO_BP_GSEA_all_terms.csv",
          "results/enrichment/GO_BP/W6_GO_BP_GSEA_significant.csv",
          "metadata/qc_W6_GO_BP_GSEA.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_go_bp_gsea_qc_outputs,
    {
      w6_go_bp_gsea_outputs
      enrichment_go_bp_gsea_qc_script

      run_stage(
        enrichment_go_bp_gsea_qc_script,
        c(
          "metadata/qc_W6_GO_BP_GSEA_result_qc.csv"
        )
      )
    },
    format = "file"
  ),

  # W6 controlled pathway readout --------------------------------------

  tar_target(
    w6_controlled_pathway_outputs,
    {
      w6_go_bp_ora_higher_outputs
      w6_go_bp_ora_lower_outputs
      w6_reactome_ora_higher_outputs
      w6_reactome_ora_lower_outputs
      w6_go_bp_gsea_outputs
      w6_go_bp_gsea_qc_outputs
      controlled_pathway_readout_script

      run_stage(
        controlled_pathway_readout_script,
        c(
          "results/interpretation/pathways/W6_GO_BP_ORA_higher_in_IPF_top20.csv",
          "results/interpretation/pathways/W6_GO_BP_ORA_lower_in_IPF_top20.csv",
          "results/interpretation/pathways/W6_Reactome_ORA_higher_in_IPF_top20.csv",
          "results/interpretation/pathways/W6_Reactome_ORA_lower_in_IPF_top20.csv",
          "results/interpretation/pathways/W6_GO_BP_GSEA_positive_NES_top20.csv",
          "results/interpretation/pathways/W6_GO_BP_GSEA_negative_NES_top20.csv",
          "results/interpretation/pathways/W6_controlled_pathway_readout_top20.csv",
          "metadata/qc_W6_controlled_pathway_readout.csv"
        )
      )
    },
    format = "file"
  ),

  # W6 pathway-theme evidence ------------------------------------------

  tar_target(
    w6_pathway_theme_outputs,
    {
      w6_controlled_pathway_outputs
      pathway_theme_matrix_script

      run_stage(
        pathway_theme_matrix_script,
        c(
          "results/interpretation/themes/W6_pathway_theme_assignments.csv",
          "results/interpretation/themes/W6_theme_direction_method_matrix.csv",
          "metadata/qc_W6_theme_direction_method_matrix.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_thematic_evidence_outputs,
    {
      w6_pathway_theme_outputs
      thematic_evidence_summary_script

      run_stage(
        thematic_evidence_summary_script,
        c(
          "results/interpretation/themes/W6_thematic_evidence_summary.csv",
          "metadata/qc_W6_thematic_evidence_summary.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_biological_evidence_outputs,
    {
      w6_thematic_evidence_outputs
      biological_evidence_classification_script

      run_stage(
        biological_evidence_classification_script,
        c(
          "results/interpretation/themes/W6_biological_evidence_classification.csv",
          "metadata/qc_W6_biological_evidence_classification.csv"
        )
      )
    },
    format = "file"
  ),

  # W6 biological synthesis --------------------------------------------

  tar_target(
    w6_biological_scaffold_outputs,
    {
      w6_biological_evidence_outputs
      biological_synthesis_scaffold_script

      run_stage(
        biological_synthesis_scaffold_script,
        c(
          "results/interpretation/synthesis/W6_biological_synthesis_scaffold.csv",
          "metadata/qc_W6_biological_synthesis_scaffold.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w6_controlled_synthesis_outputs,
    {
      w6_biological_scaffold_outputs
      controlled_biological_synthesis_script

      run_stage(
        controlled_biological_synthesis_script,
        c(
          "results/interpretation/synthesis/W6_controlled_biological_synthesis.csv",
          "metadata/qc_W6_controlled_biological_synthesis.csv"
        )
      )
    },
    format = "file"
  ),

  # W6 final biological interpretation ---------------------------------

  tar_target(
    w6_final_biological_report_outputs,
    {
      w6_controlled_synthesis_outputs
      final_biological_report_script

      run_stage(
        final_biological_report_script,
        c(
          "results/interpretation/final/W6_final_biological_report.csv",
          "results/interpretation/final/W6_final_biological_report.md",
          "metadata/qc_W6_final_biological_report.csv"
        )
      )
    },
    format = "file"
  ),

  # ------------------------------------------------------------------
  # W7: FROZEN AGE AND GENDER SENSITIVITY / ROBUSTNESS
  #
  # W7 branches from the frozen W5 primary analysis.
  # It does not depend scientifically on W6.
  #
  # The sensitivity models reuse the frozen W1-W4 cohort/filter
  # foundation and preserve W5 as the primary inferential reference.
  # ------------------------------------------------------------------

  # W7 script provenance ---------------------------------------------

  tar_target(
    w7_gender_sensitivity_script,
    "R/W7_gender_sensitivity.R",
    format = "file"
  ),

  tar_target(
    w7_continuous_age_sensitivity_script,
    "R/W7_continuous_age_sensitivity.R",
    format = "file"
  ),

  tar_target(
    w7_gender_vs_w5_script,
    "R/W7_compare_W5_vs_gender.R",
    format = "file"
  ),

  tar_target(
    w7_continuous_age_vs_w5_script,
    "R/W7_compare_W5_vs_continuous_age.R",
    format = "file"
  ),

  tar_target(
    w7_combined_sensitivity_script,
    "R/W7_combined_sensitivity_synthesis.R",
    format = "file"
  ),

  tar_target(
    w7_combined_robustness_script,
    "R/W7_combined_robustness_interpretation.R",
    format = "file"
  ),

  # W7 gender-adjusted sensitivity model -----------------------------

  tar_target(
    w7_gender_sensitivity_outputs,
    {
      raw_counts_source
      import_outputs
      normalization_outputs
      w5_primary_de_outputs
      w7_gender_sensitivity_script

      run_stage(
        w7_gender_sensitivity_script,
        c(
          "results/sensitivity/W7_gender_IPF_vs_NDC_all_genes.csv",
          "results/sensitivity/W7_gender_IPF_vs_NDC_summary.csv",
          "metadata/qc_W7_gender_sensitivity.csv"
        )
      )
    },
    format = "file"
  ),

  # W7 continuous-age sensitivity model ------------------------------

  tar_target(
    w7_continuous_age_sensitivity_outputs,
    {
      raw_counts_source
      import_outputs
      normalization_outputs
      w5_primary_de_outputs
      w7_continuous_age_sensitivity_script

      run_stage(
        w7_continuous_age_sensitivity_script,
        c(
          "results/sensitivity/W7_continuous_age_IPF_vs_NDC_all_genes.csv",
          "results/sensitivity/W7_continuous_age_IPF_vs_NDC_summary.csv",
          "metadata/qc_W7_continuous_age_sensitivity.csv"
        )
      )
    },
    format = "file"
  ),

  # Frozen W5 vs gender sensitivity comparison -----------------------

  tar_target(
    w7_gender_vs_w5_outputs,
    {
      w5_primary_de_outputs
      w7_gender_sensitivity_outputs
      w7_gender_vs_w5_script

      run_stage(
        w7_gender_vs_w5_script,
        c(
          "results/sensitivity/W7_gender_vs_W5_gene_comparison.csv",
          "results/sensitivity/W7_gender_vs_W5_summary.csv",
          "metadata/qc_W7_gender_vs_W5_comparison.csv"
        )
      )
    },
    format = "file"
  ),

  # Frozen W5 vs continuous-age sensitivity comparison ---------------

  tar_target(
    w7_continuous_age_vs_w5_outputs,
    {
      w5_primary_de_outputs
      w7_continuous_age_sensitivity_outputs
      w7_continuous_age_vs_w5_script

      run_stage(
        w7_continuous_age_vs_w5_script,
        c(
          "results/sensitivity/W7_continuous_age_vs_W5_gene_comparison.csv",
          "results/sensitivity/W7_continuous_age_vs_W5_summary.csv",
          "metadata/qc_W7_continuous_age_vs_W5_comparison.csv"
        )
      )
    },
    format = "file"
  ),

  # Combined sensitivity synthesis -----------------------------------

  tar_target(
    w7_combined_sensitivity_outputs,
    {
      w7_gender_sensitivity_outputs
      w7_continuous_age_sensitivity_outputs
      w7_gender_vs_w5_outputs
      w7_continuous_age_vs_w5_outputs
      w7_combined_sensitivity_script

      run_stage(
        w7_combined_sensitivity_script,
        c(
          "results/sensitivity/W7_combined_sensitivity_summary.csv",
          "metadata/qc_W7_combined_sensitivity_synthesis.csv"
        )
      )
    },
    format = "file"
  ),

  # Final W7 robustness interpretation -------------------------------

  tar_target(
    w7_combined_robustness_outputs,
    {
      w7_combined_sensitivity_outputs
      w7_combined_robustness_script

      run_stage(
        w7_combined_robustness_script,
        c(
          "results/sensitivity/W7_combined_robustness_interpretation.csv",
          "results/sensitivity/W7_combined_robustness_interpretation.md",
          "metadata/qc_W7_combined_robustness_interpretation.csv"
        )
      )
    },
    format = "file"
  ),

  # ------------------------------------------------------------------
  # W8: PORTFOLIO VISUALIZATION AND FINAL REPORTING
  #
  # W8 consumes frozen W5-W7 analytical outputs.
  # It does not redefine scientific inference.
  # ------------------------------------------------------------------

  tar_target(
    w8_primary_visualization_script,
    "R/W8_primary_de_visualization.R",
    format = "file"
  ),

  tar_target(
    w8_biological_evidence_matrix_script,
    "R/W8_biological_evidence_matrix.R",
    format = "file"
  ),

  tar_target(
    w8_biological_visualization_script,
    "R/W8_biological_interpretation_visualization.R",
    format = "file"
  ),

  tar_target(
    w8_robustness_visualization_script,
    "R/W8_robustness_visualization.R",
    format = "file"
  ),

  tar_target(
    w8_final_report_source,
    "reports/analysis.qmd",
    format = "file"
  ),

  # W5 presentation layer --------------------------------------------

  tar_target(
    w8_primary_visualization_outputs,
    {
      w5_primary_de_outputs
      w6_annotation_mapping_outputs
      w8_primary_visualization_script

      run_stage(
        w8_primary_visualization_script,
        c(
          "figures/final/W8_primary_IPF_vs_NDC_volcano.png",
          "figures/final/W8_primary_IPF_vs_NDC_volcano.pdf",
          "metadata/qc_W8_primary_DE_visualization_labels.csv",
          "metadata/qc_W8_primary_DE_visualization.csv"
        )
      )
    },
    format = "file"
  ),

  # W6 structured evidence layer -------------------------------------

  tar_target(
    w8_biological_evidence_matrix_outputs,
    {
      w6_final_biological_report_outputs
      w8_biological_evidence_matrix_script

      run_stage(
        w8_biological_evidence_matrix_script,
        c(
          "results/interpretation/final/W8_biological_evidence_matrix.csv",
          "metadata/qc_W8_biological_evidence_matrix.csv"
        )
      )
    },
    format = "file"
  ),

  tar_target(
    w8_biological_visualization_outputs,
    {
      w8_biological_evidence_matrix_outputs
      w8_biological_visualization_script

      run_stage(
        w8_biological_visualization_script,
        c(
          "figures/final/W8_biological_interpretation_matrix.png",
          "figures/final/W8_biological_interpretation_matrix.pdf",
          "metadata/qc_W8_biological_interpretation_visualization.csv"
        )
      )
    },
    format = "file"
  ),

  # W7 robustness presentation layer ---------------------------------

  tar_target(
    w8_robustness_visualization_outputs,
    {
      w7_combined_sensitivity_outputs
      w7_combined_robustness_outputs
      w8_robustness_visualization_script

      run_stage(
        w8_robustness_visualization_script,
        c(
          "figures/final/W8_robustness_summary.png",
          "figures/final/W8_robustness_summary.pdf",
          "metadata/qc_W8_robustness_visualization.csv"
        )
      )
    },
    format = "file"
  ),

  # Final report ------------------------------------------------------
  #
  # reports/analysis.qmd has execute.enabled = false.
  # Quarto therefore consumes frozen scientific outputs and W8
  # presentation figures without rerunning analytical models.

  tar_target(
    w8_final_report_outputs,
    {
      w5_primary_de_outputs
      w6_final_biological_report_outputs
      w7_combined_robustness_outputs

      w8_primary_visualization_outputs
      w8_biological_evidence_matrix_outputs
      w8_biological_visualization_outputs
      w8_robustness_visualization_outputs

      w8_final_report_source

      run_quarto_report(
        w8_final_report_source,
        c(
          "reports/final_report.html",
          "reports/final_report.pdf"
        )
      )
    },
    format = "file"
  )
)
