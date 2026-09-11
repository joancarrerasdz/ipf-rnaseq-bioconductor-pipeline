# Decision Log

This file records methodological and scope decisions that affect the interpretation, validity or reproducibility of the project.

---

## 2026-09-10 — Use GSE213001 as the primary dataset

**Status:** ACCEPTED

GSE213001 will be used as the public bulk RNA-seq dataset for this project.

The complete downstream analysis will use the public raw-count matrix and GEO sample metadata.

Full FASTQ reprocessing is outside the core analysis scope. A small FASTQ subset may later be used as an end-to-end workflow demonstration.

---

## 2026-09-10 — Preserve source files unchanged

**Status:** ACCEPTED

Files downloaded from GEO are stored locally under:

`data/source/`

and excluded from Git.

Source files must not be manually edited.

Derived data must be generated programmatically.

Source URLs, download information and SHA-256 checksums are recorded in:

`data/provenance.tsv`

---

## 2026-09-10 — Raw counts are the primary expression input

**Status:** ACCEPTED

The primary expression input is:

`GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz`

Filtering and normalization will be performed reproducibly within this project.

The normalized expression file distributed by GEO is not the primary input for differential-expression analysis.

---

## 2026-09-10 — Count-matrix identifiers are Ensembl-like IDs

**Status:** ACCEPTED / DOCUMENTED

Although the GEO supplementary filename contains the string `Entrez-IDs`, inspection of the first column shows identifiers of the form:

`ENSG...`

All 15,065 gene identifiers follow an Ensembl-like pattern.

Gene annotation and enrichment steps must therefore verify identifier type explicitly rather than infer it from the supplementary filename.

---

## 2026-09-10 — Raw-count matrix passes basic integrity checks

**Status:** ACCEPTED

The raw-count matrix contains:

- 15,065 gene rows
- 139 sample columns
- no missing count values
- no negative count values
- no non-integer count values
- no duplicated sample-column names
- no missing gene identifiers
- no duplicated gene identifiers

The reported GEO library sizes were independently validated against raw-count column sums.

After sample alignment:

- all library sizes matched exactly
- maximum library-size difference was 0

---

## 2026-09-10 — Metadata and counts contain the same 139 samples

**Status:** ACCEPTED

The raw-count matrix contains 139 sample columns.

The GEO Series Matrix contains 139 sample records.

All 139 count-matrix sample names match the GEO `sample_title` values exactly as a set.

The GSM accessions are retained as external identifiers but are not the column identifiers used in the raw-count matrix.

---

## 2026-09-10 — Metadata must be explicitly reordered to count-matrix columns

**Status:** CORRECTED

Initial validation showed that the GEO Series Matrix and raw-count matrix contain the same 139 samples but in different orders.

Before correction:

- 139 of 139 sample names matched as a set
- 138 positions differed in ordering

Therefore, positional matching is prohibited.

`R/import.R` now explicitly aligns metadata to raw-count column order using `sample_title`.

Post-correction validation showed:

- exact metadata/count order: TRUE
- all library sizes exact: TRUE
- maximum library-size difference: 0

---

## 2026-09-10 — Complete metadata contain four disease groups

**Status:** DOCUMENTED

Parsing the complete GEO Series Matrix identified:

| Disease group | Samples | Donors |
|---|---:|---:|
| IPF | 62 | 20 |
| NDC | 41 | 14 |
| ILD | 26 | 9 |
| CLAD | 10 | 3 |
| **Total** | **139** | **46** |

The CLAD samples are retained in the complete metadata for provenance and reproducibility but are outside the primary IPF analysis scope.

---

## 2026-09-10 — Use `diseasegroup` as the canonical disease variable

**Status:** ACCEPTED

Several related disease annotations are available.

Cross-checking demonstrated that the major disease grouping is consistent across:

- `diseasegroup`
- `diagnosis`
- `diseasegroup1`

`diseasegroup` is used as the canonical disease variable because it directly represents the analytical disease groups required for this project.

`group` is not used as the primary disease variable because it combines disease group and anatomical lung region.

`diseasenormal` is retained as a broader disease-versus-normal descriptor.

`diseasesubtype` is retained as secondary clinical metadata.

---

## 2026-09-10 — Donor defines the repeated-measures blocking unit

**Status:** ACCEPTED

The complete dataset contains:

- 139 tissue samples
- 46 unique donors

44 of 46 donors contribute more than one tissue sample.

Individual donors contribute between 1 and 5 samples.

Treating all 139 tissue samples as independent biological replicates would therefore produce pseudo-replication.

`donorid` will be used to represent within-donor dependence in downstream statistical models.

---

## 2026-09-10 — Primary contrast is IPF vs NDC

**Status:** PRE-SPECIFIED

The primary disease comparison is:

**IPF vs NDC**

Before accounting for missing lung-region information, this cohort contains:

- 62 IPF samples / 20 donors
- 41 NDC samples / 14 donors
- 103 samples / 34 donors

ILD and CLAD samples are not part of the primary contrast.

---

## 2026-09-10 — Lung region is a pre-specified primary-model variable

**Status:** PRE-SPECIFIED

The primary disease analysis will account for anatomical lung region.

Two samples in the IPF/NDC cohort have unknown lung-region annotations:

- `ALF018E` / `GSM6568369` / IPF
- `ALF026E` / `GSM6568412` / NDC

No Apex/Base value will be imputed for these samples.

The currently planned primary-analysis cohort therefore contains:

- 101 samples
- 34 donors
- 61 IPF samples
- 40 NDC samples

Region distribution:

| Disease group | Apex | Base | Total |
|---|---:|---:|---:|
| IPF | 30 | 31 | 61 |
| NDC | 22 | 18 | 40 |
| **Total** | **52** | **49** | **101** |

---

## 2026-09-10 — Secondary contrast is Apex vs Base within IPF

**Status:** PRE-SPECIFIED

The secondary analysis will evaluate transcriptional differences associated with anatomical lung region within IPF.

Among the 20 IPF donors:

- 19 have samples from both Apex and Base
- 1 donor (`ALF001`) has only a Base sample

Some donors contribute multiple samples within the same anatomical region.

Therefore, the regional comparison is not a simple one-to-one paired design and requires a donor-aware repeated-measures approach.

---

## 2026-09-10 — Do not model donor as a conventional fixed effect in the primary disease model

**Status:** PRE-SPECIFIED

Disease group is defined at donor level.

A conventional fixed-effect term for every donor would therefore be confounded with the IPF-versus-NDC disease effect.

Within-donor dependence will instead be handled through an appropriate blocking, correlation or random-effect strategy.

The exact statistical implementation will be finalized before formal differential-expression testing.

---

## 2026-09-10 — Use one principal differential-expression framework

**Status:** PRE-SPECIFIED

The project will not run DESeq2, edgeR and limma simply to demonstrate multiple software packages.

The current preferred framework is a donor-aware limma-voom analysis using an appropriate blocking/correlation strategy or a statistically equivalent formulation.

The final implementation will depend on QC, filtering, normalization and model diagnostics.

---

## 2026-09-10 — Do not reuse source normalization factors automatically

**Status:** ACCEPTED

The GEO metadata contain source-study normalization factors.

These values are retained for provenance and validation.

They will not automatically be reused for the primary analysis.

Filtering and normalization choices will be implemented and justified reproducibly within this project.

---

## 2026-09-10 — QC exclusions cannot be chosen to improve biological separation

**Status:** ACCEPTED

Samples will not be excluded because their removal improves PCA separation, differential-expression significance or visual agreement with an expected biological pattern.

Any QC-based exclusion must be supported by independently defined technical or data-quality evidence and documented before formal differential-expression testing.

---

# Open decisions

The following decisions remain intentionally unresolved:

- final low-count filtering rule
- normalization method and parameters
- final donor-correlation / repeated-measures implementation
- treatment of technical or biological covariates beyond lung region
- any additional QC-based sample exclusions
- exact Ensembl-to-gene annotation strategy
- handling of any unmapped or ambiguous gene identifiers
- enrichment gene universe
- exact GSEA and/or over-representation workflow

These items will be resolved at the appropriate downstream stage and recorded in this decision log before they affect formal results.
