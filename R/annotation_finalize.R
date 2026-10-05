# ============================================================
# GSE213001 W6.4 — Annotated DE table
#
# Purpose:
#   Attach the frozen W6 annotation layer to the frozen W5
#   primary differential-expression table.
#
# Critical invariant:
#   The 15,012-row W5 DE table must remain unchanged apart
#   from the addition of annotation/QC columns.
#
# No biological interpretation is performed.
# ============================================================

source("renv/activate.R")

# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------

de <- read.csv(
  "results/differential_expression/primary_IPF_vs_NDC_all_genes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annotation <- read.csv(
  "metadata/qc_annotation_gene_mapping.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("")
)

de_original <- de

cat("\n============================================================\n")
cat("GSE213001 W6.4 — ANNOTATED DE TABLE\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 2. Input invariants
# ------------------------------------------------------------

stopifnot(
  "gene_id" %in% names(de),
  "gene_id" %in% names(annotation),

  nrow(de) == 15012,
  nrow(annotation) == 15012,

  length(unique(de$gene_id)) == 15012,
  length(unique(annotation$gene_id)) == 15012,

  !anyDuplicated(de$gene_id),
  !anyDuplicated(annotation$gene_id)
)

cat("DE rows:                    ", nrow(de), "\n", sep = "")
cat("Annotation rows:            ", nrow(annotation), "\n", sep = "")
cat("Unique DE gene IDs:         ", length(unique(de$gene_id)), "\n", sep = "")
cat("Unique annotation gene IDs: ", length(unique(annotation$gene_id)), "\n", sep = "")

# ------------------------------------------------------------
# 3. Universe equivalence
# ------------------------------------------------------------

missing_from_annotation <- setdiff(
  de$gene_id,
  annotation$gene_id
)

extra_in_annotation <- setdiff(
  annotation$gene_id,
  de$gene_id
)

stopifnot(
  length(missing_from_annotation) == 0L,
  length(extra_in_annotation) == 0L
)

cat("\n=== GENE-UNIVERSE CHECK ===\n")
cat("Missing from annotation:    ", length(missing_from_annotation), "\n", sep = "")
cat("Extra in annotation:        ", length(extra_in_annotation), "\n", sep = "")

# ------------------------------------------------------------
# 4. Order-preserving annotation join
# ------------------------------------------------------------

idx <- match(
  de$gene_id,
  annotation$gene_id
)

stopifnot(
  !anyNA(idx),
  identical(
    annotation$gene_id[idx],
    de$gene_id
  )
)

annotation_fields <- c(
  "SYMBOL",
  "ENTREZID",
  "GENENAME",
  "n_symbol",
  "n_entrezid",
  "n_genename",
  "n_mapping_records",
  "mapping_status"
)

stopifnot(
  all(annotation_fields %in% names(annotation))
)

annotation_aligned <- annotation[
  idx,
  annotation_fields,
  drop = FALSE
]

annotated_de <- cbind(
  de,
  annotation_aligned
)

# ------------------------------------------------------------
# 5. Structural integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(annotated_de) == 15012,
  identical(annotated_de$gene_id, de_original$gene_id),
  identical(
    annotated_de[, names(de_original), drop = FALSE],
    de_original
  )
)

statistical_columns_unchanged <- identical(
  annotated_de[, names(de_original), drop = FALSE],
  de_original
)

gene_order_unchanged <- identical(
  annotated_de$gene_id,
  de_original$gene_id
)

# ------------------------------------------------------------
# 6. Annotation summary
# ------------------------------------------------------------

n_unmapped <- sum(
  annotated_de$mapping_status == "unmapped",
  na.rm = TRUE
)

n_single <- sum(
  annotated_de$mapping_status == "single_mapping_record",
  na.rm = TRUE
)

n_multiple <- sum(
  annotated_de$mapping_status == "multiple_mapping_records",
  na.rm = TRUE
)

n_missing_symbol <- sum(is.na(annotated_de$SYMBOL))

# ------------------------------------------------------------
# 7. Output
# ------------------------------------------------------------

output_file <- paste0(
  "results/differential_expression/",
  "primary_IPF_vs_NDC_annotated.csv"
)

write.csv(
  annotated_de,
  output_file,
  row.names = FALSE,
  na = ""
)

audit <- data.frame(
  metric = c(
    "input_DE_rows",
    "annotation_rows",
    "output_rows",
    "unique_output_gene_ids",
    "missing_from_annotation",
    "extra_in_annotation",
    "gene_order_unchanged",
    "statistical_columns_unchanged",
    "single_mapping_record_genes",
    "multiple_mapping_record_genes",
    "unmapped_genes",
    "missing_SYMBOL"
  ),
  value = c(
    nrow(de_original),
    nrow(annotation),
    nrow(annotated_de),
    length(unique(annotated_de$gene_id)),
    length(missing_from_annotation),
    length(extra_in_annotation),
    gene_order_unchanged,
    statistical_columns_unchanged,
    n_single,
    n_multiple,
    n_unmapped,
    n_missing_symbol
  ),
  stringsAsFactors = FALSE
)

audit_file <- "metadata/qc_annotation_join_audit.csv"

write.csv(
  audit,
  audit_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Final console certificate
# ------------------------------------------------------------

cat("\n=== STRUCTURAL INTEGRITY ===\n")

cat("Input DE rows:               ", nrow(de_original), "\n", sep = "")
cat("Output annotated rows:       ", nrow(annotated_de), "\n", sep = "")
cat(
  "Unique output gene IDs:      ",
  length(unique(annotated_de$gene_id)),
  "\n",
  sep = ""
)

cat(
  "Gene order unchanged:        ",
  gene_order_unchanged,
  "\n",
  sep = ""
)

cat(
  "W5 columns unchanged:        ",
  statistical_columns_unchanged,
  "\n",
  sep = ""
)

cat("\n=== ANNOTATION STATUS ===\n")

cat("Single mapping record:       ", n_single, "\n", sep = "")
cat("Multiple mapping records:    ", n_multiple, "\n", sep = "")
cat("Unmapped genes:              ", n_unmapped, "\n", sep = "")
cat("Missing SYMBOL:              ", n_missing_symbol, "\n", sep = "")

cat("\n=== OUTPUTS ===\n")
cat(output_file, "\n")
cat(audit_file, "\n")

cat("\nIMPORTANT:\n")
cat("The frozen W5 DE results were not recalculated.\n")
cat("The statistical gene universe remains 15,012 genes.\n")
cat("Annotation was attached by exact Ensembl gene_id matching.\n")
cat("No gene was excluded or aggregated.\n")
cat("No biological interpretation was performed.\n")

cat("\n============================================================\n")
cat("W6.4 ANNOTATED DE TABLE: PASSED\n")
cat("============================================================\n")
