# ============================================================
# GSE213001 W6.3 — Gene annotation mapping audit
#
# Frozen universe:
#   15,012 Ensembl gene IDs from the W5 primary DE analysis.
#
# Annotation:
#   AnnotationDbi 1.70.0
#   org.Hs.eg.db 3.21.0
#   keytype = ENSEMBL
#
# This stage performs annotation mapping and mapping QC only.
# No biological interpretation is performed.
# ============================================================

source("renv/activate.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# 1. Frozen W5 gene universe
# ------------------------------------------------------------

de <- read.csv(
  "results/differential_expression/primary_IPF_vs_NDC_all_genes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  "gene_id" %in% names(de),
  nrow(de) == 15012,
  length(unique(de$gene_id)) == 15012,
  all(grepl("^ENSG[0-9]+$", de$gene_id))
)

gene_ids <- de$gene_id

cat("\n============================================================\n")
cat("GSE213001 W6.3 — GENE ANNOTATION MAPPING AUDIT\n")
cat("============================================================\n\n")

cat("Frozen W5 genes:       ", length(gene_ids), "\n", sep = "")
cat("Unique Ensembl IDs:     ", length(unique(gene_ids)), "\n", sep = "")

# ------------------------------------------------------------
# 2. Raw annotation mapping
# ------------------------------------------------------------

raw_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENSEMBL",
  columns = c(
    "SYMBOL",
    "ENTREZID",
    "GENENAME"
  )
)

names(raw_map)[names(raw_map) == "ENSEMBL"] <- "gene_id"

# Remove exact duplicated mapping records only.
raw_map <- unique(raw_map)

stopifnot(
  all(raw_map$gene_id %in% gene_ids)
)

# ------------------------------------------------------------
# 3. Mapping helper functions
# ------------------------------------------------------------

valid_values <- function(x) {
  sort(unique(x[!is.na(x) & nzchar(x)]))
}

collapse_values <- function(x) {
  x <- valid_values(x)

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(x, collapse = ";")
}

count_values <- function(x) {
  length(valid_values(x))
}

# ------------------------------------------------------------
# 4. One-row-per-Ensembl mapping audit
# ------------------------------------------------------------

map_split <- split(raw_map, raw_map$gene_id)

gene_audit <- lapply(
  gene_ids,
  function(id) {

    x <- map_split[[id]]

    if (is.null(x)) {
      return(
        data.frame(
          gene_id = id,
          SYMBOL = NA_character_,
          ENTREZID = NA_character_,
          GENENAME = NA_character_,
          n_symbol = 0L,
          n_entrezid = 0L,
          n_genename = 0L,
          n_mapping_records = 0L,
          mapping_status = "unmapped",
          stringsAsFactors = FALSE
        )
      )
    }

    has_annotation <- (
      !is.na(x$SYMBOL) |
      !is.na(x$ENTREZID) |
      !is.na(x$GENENAME)
    )

    informative <- unique(x[has_annotation, , drop = FALSE])

    n_records <- nrow(informative)

    status <- if (n_records == 0L) {
      "unmapped"
    } else if (n_records == 1L) {
      "single_mapping_record"
    } else {
      "multiple_mapping_records"
    }

    data.frame(
      gene_id = id,
      SYMBOL = collapse_values(x$SYMBOL),
      ENTREZID = collapse_values(x$ENTREZID),
      GENENAME = collapse_values(x$GENENAME),
      n_symbol = count_values(x$SYMBOL),
      n_entrezid = count_values(x$ENTREZID),
      n_genename = count_values(x$GENENAME),
      n_mapping_records = n_records,
      mapping_status = status,
      stringsAsFactors = FALSE
    )
  }
)

gene_audit <- do.call(rbind, gene_audit)

stopifnot(
  nrow(gene_audit) == 15012,
  length(unique(gene_audit$gene_id)) == 15012,
  identical(gene_audit$gene_id, gene_ids)
)

# ------------------------------------------------------------
# 5. Shared SYMBOL audit
# ------------------------------------------------------------

symbol_pairs <- unique(
  raw_map[
    !is.na(raw_map$SYMBOL) &
      nzchar(raw_map$SYMBOL),
    c("gene_id", "SYMBOL"),
    drop = FALSE
  ]
)

symbol_gene_counts <- table(symbol_pairs$SYMBOL)

shared_symbol_names <- names(
  symbol_gene_counts[symbol_gene_counts > 1L]
)

if (length(shared_symbol_names) > 0L) {

  shared_symbols <- data.frame(
    SYMBOL = shared_symbol_names,
    n_ensembl_gene_ids = as.integer(
      symbol_gene_counts[shared_symbol_names]
    ),
    stringsAsFactors = FALSE
  )

  shared_symbols <- shared_symbols[
    order(
      -shared_symbols$n_ensembl_gene_ids,
      shared_symbols$SYMBOL
    ),
    ,
    drop = FALSE
  ]

} else {

  shared_symbols <- data.frame(
    SYMBOL = character(),
    n_ensembl_gene_ids = integer(),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 6. Quantitative mapping summary
# ------------------------------------------------------------

n_unmapped <- sum(
  gene_audit$mapping_status == "unmapped"
)

n_single <- sum(
  gene_audit$mapping_status == "single_mapping_record"
)

n_multiple <- sum(
  gene_audit$mapping_status == "multiple_mapping_records"
)

n_symbol_missing <- sum(is.na(gene_audit$SYMBOL))
n_entrez_missing <- sum(is.na(gene_audit$ENTREZID))
n_genename_missing <- sum(is.na(gene_audit$GENENAME))

n_multi_symbol <- sum(gene_audit$n_symbol > 1L)
n_multi_entrez <- sum(gene_audit$n_entrezid > 1L)
n_multi_genename <- sum(gene_audit$n_genename > 1L)

summary_df <- data.frame(
  metric = c(
    "input_genes",
    "unique_input_genes",
    "unmapped_genes",
    "single_mapping_record_genes",
    "multiple_mapping_record_genes",
    "missing_SYMBOL",
    "missing_ENTREZID",
    "missing_GENENAME",
    "genes_with_multiple_SYMBOLs",
    "genes_with_multiple_ENTREZIDs",
    "genes_with_multiple_GENENAMEs",
    "symbols_shared_by_multiple_ensembl_ids"
  ),
  value = c(
    length(gene_ids),
    length(unique(gene_ids)),
    n_unmapped,
    n_single,
    n_multiple,
    n_symbol_missing,
    n_entrez_missing,
    n_genename_missing,
    n_multi_symbol,
    n_multi_entrez,
    n_multi_genename,
    nrow(shared_symbols)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 7. Outputs
# ------------------------------------------------------------

dir.create(
  "metadata",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  raw_map,
  "metadata/qc_annotation_raw_mapping.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  gene_audit,
  "metadata/qc_annotation_gene_mapping.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  shared_symbols,
  "metadata/qc_annotation_shared_symbols.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  summary_df,
  "metadata/qc_annotation_summary.csv",
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 8. Console audit
# ------------------------------------------------------------

cat("\n=== MAPPING SUMMARY ===\n")

cat("Input genes:                  ", length(gene_ids), "\n", sep = "")
cat("Unmapped genes:               ", n_unmapped, "\n", sep = "")
cat("Single mapping record:        ", n_single, "\n", sep = "")
cat("Multiple mapping records:     ", n_multiple, "\n", sep = "")

cat("\n=== ANNOTATION COMPLETENESS ===\n")

cat("Missing SYMBOL:               ", n_symbol_missing, "\n", sep = "")
cat("Missing ENTREZID:             ", n_entrez_missing, "\n", sep = "")
cat("Missing GENENAME:             ", n_genename_missing, "\n", sep = "")

cat("\n=== AMBIGUITY AUDIT ===\n")

cat("Genes with >1 SYMBOL:         ", n_multi_symbol, "\n", sep = "")
cat("Genes with >1 ENTREZID:       ", n_multi_entrez, "\n", sep = "")
cat("Genes with >1 GENENAME:       ", n_multi_genename, "\n", sep = "")
cat(
  "SYMBOLs shared across ENSGs:  ",
  nrow(shared_symbols),
  "\n",
  sep = ""
)

cat("\n=== OUTPUTS ===\n")

cat("metadata/qc_annotation_raw_mapping.csv\n")
cat("metadata/qc_annotation_gene_mapping.csv\n")
cat("metadata/qc_annotation_shared_symbols.csv\n")
cat("metadata/qc_annotation_summary.csv\n")

cat("\nIMPORTANT:\n")
cat("The frozen 15,012-gene DE universe was not modified.\n")
cat("No genes were excluded because of annotation status.\n")
cat("No biological interpretation was performed.\n")

cat("\n============================================================\n")
cat("W6.3 GENE ANNOTATION MAPPING AUDIT: COMPLETED\n")
cat("============================================================\n")
