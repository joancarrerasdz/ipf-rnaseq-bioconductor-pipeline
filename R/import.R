# ==============================================================================
# GSE213001 metadata import
# ==============================================================================
#
# Purpose:
#   Parse sample metadata directly from the GEO Series Matrix for GSE213001.
#
# Inputs:
#   data/source/GSE213001_series_matrix.txt.gz
#
# Outputs:
#   metadata/sample_metadata.csv
#   metadata/donor_sample_map.csv
#
# Notes:
#   - Source files are never modified.
#   - Parsing uses base R only at this stage.
#   - Repeated samples from the same donor are explicitly preserved.
#
# ==============================================================================

series_file <- file.path(
  "data",
  "source",
  "GSE213001_series_matrix.txt.gz"
)

metadata_dir <- "metadata"

if (!file.exists(series_file)) {
  stop("Series Matrix file not found: ", series_file)
}

dir.create(metadata_dir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

strip_quotes <- function(x) {
  x <- trimws(x)
  x <- sub('^"', "", x)
  x <- sub('"$', "", x)
  x
}


split_geo_line <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  strip_quotes(fields)
}


get_sample_field <- function(lines, field) {

  idx <- which(startsWith(lines, paste0(field, "\t")))

  if (length(idx) != 1L) {
    stop(
      "Expected exactly one GEO field '",
      field,
      "', found ",
      length(idx)
    )
  }

  values <- split_geo_line(lines[idx])

  values[-1]
}


clean_field_name <- function(x) {

  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  x
}


replace_missing <- function(x) {

  x <- trimws(x)

  x[x %in% c(
    "",
    "NA",
    "N/A",
    "na",
    "n/a",
    "Unknown",
    "unknown"
  )] <- NA_character_

  x
}


to_numeric_checked <- function(x, column_name) {

  original <- x
  converted <- suppressWarnings(as.numeric(x))

  bad <- !is.na(original) & is.na(converted)

  if (any(bad)) {
    warning(
      "Non-numeric values found in numeric column '",
      column_name,
      "': ",
      paste(unique(original[bad]), collapse = ", ")
    )
  }

  converted
}


# ------------------------------------------------------------------------------
# Read GEO Series Matrix
# ------------------------------------------------------------------------------

con <- gzfile(series_file, open = "rt")
lines <- readLines(con, warn = FALSE)
close(con)

sample_accession <- get_sample_field(
  lines,
  "!Sample_geo_accession"
)

sample_title <- get_sample_field(
  lines,
  "!Sample_title"
)

n_samples <- length(sample_accession)

if (n_samples != 139L) {
  stop(
    "Expected 139 GEO samples, found ",
    n_samples
  )
}

if (length(sample_title) != n_samples) {
  stop("Sample title/accession length mismatch.")
}

if (anyDuplicated(sample_accession)) {
  stop("Duplicated GEO accessions detected.")
}

if (anyDuplicated(sample_title)) {
  stop("Duplicated sample titles detected.")
}


# ------------------------------------------------------------------------------
# Parse !Sample_characteristics_ch1 fields
# ------------------------------------------------------------------------------

characteristic_lines <- lines[
  startsWith(lines, "!Sample_characteristics_ch1\t")
]

if (length(characteristic_lines) == 0L) {
  stop("No !Sample_characteristics_ch1 fields found.")
}

metadata <- data.frame(
  sample_title = sample_title,
  geo_accession = sample_accession,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

seen_names <- character(0)

for (line in characteristic_lines) {

  fields <- split_geo_line(line)

  values <- fields[-1]

  if (length(values) != n_samples) {
    stop(
      "Characteristics row has ",
      length(values),
      " values; expected ",
      n_samples
    )
  }

  raw_keys <- sub(":.*$", "", values)
  raw_keys <- trimws(raw_keys)

  unique_keys <- unique(raw_keys)

  if (length(unique_keys) != 1L) {
    stop(
      "Inconsistent characteristic names within GEO row: ",
      paste(unique_keys, collapse = ", ")
    )
  }

  field_name <- clean_field_name(unique_keys)

  if (field_name %in% seen_names) {
    stop("Duplicated metadata characteristic: ", field_name)
  }

  values <- sub("^[^:]+:[[:space:]]*", "", values)
  values <- replace_missing(values)

  metadata[[field_name]] <- values
  seen_names <- c(seen_names, field_name)
}


# ------------------------------------------------------------------------------
# Convert known numeric metadata fields
# ------------------------------------------------------------------------------

numeric_fields <- c(
  "lib_size",
  "norm_factors",
  "age",
  "lungweight_mg",
  "rnaconcentration_ngul",
  "rin",
  "volum_ul",
  "fvcprebd_perc",
  "dlco_perc",
  "dlco_perc_corrected"
)

for (field in intersect(numeric_fields, names(metadata))) {
  metadata[[field]] <- to_numeric_checked(
    metadata[[field]],
    field
  )
}


# ------------------------------------------------------------------------------
# Core metadata validation
# ------------------------------------------------------------------------------

required_fields <- c(
  "sample_title",
  "geo_accession",
  "donorid",
  "diseasegroup",
  "lunglocation"
)

missing_required <- setdiff(
  required_fields,
  names(metadata)
)

if (length(missing_required) > 0L) {
  stop(
    "Required metadata fields missing: ",
    paste(missing_required, collapse = ", ")
  )
}


if (any(is.na(metadata$donorid))) {
  stop("Missing donor IDs detected.")
}


# Each donor should belong to only one disease group.
donor_group_n <- tapply(
  metadata$diseasegroup,
  metadata$donorid,
  function(x) length(unique(na.omit(x)))
)

if (any(donor_group_n > 1L)) {
  stop(
    "At least one donor is associated with multiple disease groups."
  )
}


# ------------------------------------------------------------------------------
# Align metadata to raw-count column order
# ------------------------------------------------------------------------------

counts_file <- file.path(
  "data",
  "source",
  "GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz"
)

if (!file.exists(counts_file)) {
  stop("Raw-count matrix not found: ", counts_file)
}

count_header <- read.csv(
  gzfile(counts_file),
  nrows = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

count_sample_names <- names(count_header)[-1]

if (length(count_sample_names) != n_samples) {
  stop(
    "Count matrix contains ",
    length(count_sample_names),
    " sample columns; expected ",
    n_samples
  )
}

missing_from_metadata <- setdiff(
  count_sample_names,
  metadata$sample_title
)

missing_from_counts <- setdiff(
  metadata$sample_title,
  count_sample_names
)

if (length(missing_from_metadata) > 0L) {
  stop(
    "Count-matrix samples missing from metadata: ",
    paste(missing_from_metadata, collapse = ", ")
  )
}

if (length(missing_from_counts) > 0L) {
  stop(
    "Metadata samples missing from count matrix: ",
    paste(missing_from_counts, collapse = ", ")
  )
}

metadata <- metadata[
  match(
    count_sample_names,
    metadata$sample_title
  ),
  ,
  drop = FALSE
]

row.names(metadata) <- NULL

if (!identical(
  metadata$sample_title,
  count_sample_names
)) {
  stop("Failed to align metadata to count-matrix column order.")
}

# ------------------------------------------------------------------------------
# Donor/sample map
# ------------------------------------------------------------------------------

map_fields <- intersect(
  c(
    "donorid",
    "sample_title",
    "geo_accession",
    "diseasegroup",
    "diagnosis",
    "diseasesubtype",
    "lunglocation",
    "leftright",
    "gender",
    "age",
    "rin",
    "smokingstatus",
    "severity",
    "fvcprebd_perc",
    "dlco_perc",
    "dlco_perc_corrected"
  ),
  names(metadata)
)

donor_sample_map <- metadata[, map_fields, drop = FALSE]

donor_sample_map <- donor_sample_map[
  order(
    donor_sample_map$donorid,
    donor_sample_map$sample_title
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------------------

write.csv(
  metadata,
  file = file.path(
    metadata_dir,
    "sample_metadata.csv"
  ),
  row.names = FALSE,
  na = "NA"
)

write.csv(
  donor_sample_map,
  file = file.path(
    metadata_dir,
    "donor_sample_map.csv"
  ),
  row.names = FALSE,
  na = "NA"
)


# ------------------------------------------------------------------------------
# Summary / sanity checks
# ------------------------------------------------------------------------------

donor_counts <- table(metadata$donorid)

donor_level <- unique(
  metadata[, c("donorid", "diseasegroup")]
)

cat("\n")
cat("============================================================\n")
cat("GSE213001 metadata import completed\n")
cat("============================================================\n")

cat("Samples:                  ", nrow(metadata), "\n")
cat(
  "Unique donors:            ",
  length(unique(metadata$donorid)),
  "\n"
)
cat(
  "Donors with >1 sample:    ",
  sum(donor_counts > 1),
  "\n"
)
cat(
  "Samples per donor range:  ",
  min(donor_counts),
  "-",
  max(donor_counts),
  "\n"
)

cat("\nSamples by disease group:\n")
print(
  table(
    metadata$diseasegroup,
    useNA = "ifany"
  )
)

cat("\nUnique donors by disease group:\n")
print(
  table(
    donor_level$diseasegroup,
    useNA = "ifany"
  )
)

cat("\nSamples by lung location:\n")
print(
  table(
    metadata$lunglocation,
    useNA = "ifany"
  )
)

cat("\nMetadata columns extracted:\n")
cat(
  paste(names(metadata), collapse = "\n"),
  "\n"
)

cat("\nOutputs:\n")
cat("  metadata/sample_metadata.csv\n")
cat("  metadata/donor_sample_map.csv\n")

cat("============================================================\n")
