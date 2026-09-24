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

---

## 2026-09-17 — W3 metadata and library-level QC completed without sample exclusion

**Status:** ACCEPTED

Metadata and library-level QC were completed before filtering, normalization, PCA/MDS or differential-expression analysis.

The complete dataset remains:

- 139 samples
- 46 donors
- 15,065 genes

No sample was excluded during W3.

The QC workflow is descriptive and does not automatically remove observations.

---

## 2026-09-17 — Technical screening flags are a watchlist, not exclusion criteria

**Status:** ACCEPTED

Robust screening identified six samples with at least one technical flag:

- `ALF009B`
- `ALF012F`
- `ALF030A`
- `ALF048D`
- `ALF017A`
- `ALF024A`

The flags were based on extreme values in library size, detected genes and/or zero-count fraction.

These samples remain in the dataset.

Within-donor review did not provide sufficiently strong independent evidence to justify exclusion.

In particular, `ALF009B` and `ALF012F` showed elevated zero-count fractions relative to their donor medians, but only modest reductions in detected genes and no extreme within-donor library-size deficit.

Therefore, all six samples are retained as a QC watchlist for later PCA/MDS and model diagnostics.

---

## 2026-09-17 — Low RIN is not an automatic exclusion criterion

**Status:** ACCEPTED

RIN is available for all 139 samples.

Seven samples have RIN < 5.

However, sample-level Spearman correlations between RIN and technical metrics were weak:

- RIN vs detected genes: approximately -0.05
- RIN vs zero-count fraction: approximately 0.05
- RIN vs log10 library size: approximately 0.14

Only two of the seven low-RIN samples were identified by the robust technical screen.

Therefore, no fixed RIN cutoff will be used as an automatic exclusion rule.

RIN will remain available for QC interpretation and possible sensitivity analysis.

---

## 2026-09-17 — Library-size variation does not currently justify exclusion

**Status:** ACCEPTED

Raw-count library sizes range approximately from 7.6 million to 19.8 million counts.

The disease groups show substantial overlap in library-size distributions.

The lowest-depth sample (`ALF017A`) and highest-depth sample (`ALF024A`) do not show sufficient independent evidence of technical failure to justify removal.

Library-size differences will instead be addressed through the downstream normalization workflow.

---

## 2026-09-17 — Age is an important candidate confounder

**Status:** PRE-SPECIFIED FOR SENSITIVITY ANALYSIS

Age is constant within donor.

At donor level in the primary IPF-vs-NDC cohort:

- IPF: 20 donors, mean age 62.55 years, median 62.5 years, range 43–73
- NDC: 13 donors with known age, mean age 48.0 years, median 44 years, range 28–69
- one NDC donor (`ALF017`) has missing age

The observed age distributions therefore differ substantially between disease groups.

Age will not replace the pre-specified primary model.

The planned primary model remains conceptually:

    expression ~ diseasegroup + lunglocation

with donor-aware repeated-measures handling.

A sensitivity model including age will be evaluated after filtering, normalization and model diagnostics.

---

## 2026-09-17 — Gender does not currently require forced adjustment

**Status:** ACCEPTED

Gender is constant within donor and complete in the primary cohort.

At donor level:

- IPF: 7 female / 13 male
- NDC: 6 female / 8 male

No strong imbalance was identified that currently requires gender to be included automatically in the primary model.

Gender remains available for descriptive and sensitivity analyses if later diagnostics justify its use.

---

## 2026-09-17 — Clinical variables with substantial or structural missingness are not primary covariates

**Status:** ACCEPTED

The complete metadata show substantial missingness for:

- corrected DLCO % predicted
- DLCO % predicted
- FVC % predicted
- disease severity
- smoking status

Disease severity is structurally unavailable for NDC controls in the primary cohort.

These variables will therefore not be added automatically to the primary IPF-vs-NDC model.

They may be used in explicitly labeled secondary or exploratory analyses when scientifically appropriate.

---

## 2026-09-17 — Processing date is retained as a technical QC variable

**Status:** ACCEPTED

Processing date is complete and constant within donor.

All three observed processing dates contain both IPF and NDC samples and are therefore not perfectly confounded with disease group.

Processing date will not be added automatically to the primary model.

It will be evaluated against PCA/MDS and other downstream QC diagnostics before any adjustment decision is made.

---

## 2026-09-17 — W3 QC figures accepted for the repository

**Status:** ACCEPTED

The following W3 figures passed visual QA:

- `figures/qc/metadata_missingness.png`
- `figures/qc/library_size_by_disease.png`
- `figures/qc/rin_by_disease.png`
- `figures/qc/donor_age_ipf_vs_ndc.png`

These figures document metadata completeness, raw library-size distributions, RIN distributions and donor-level age imbalance before filtering and normalization.

---

## 2026-09-24 — W4 low-count filtering rule frozen

**Status:** ACCEPTED

Low-count filtering is defined using `edgeR::filterByExpr()` on the
pre-specified primary analysis cohort.

The filtering cohort contains:

- 101 samples
- 34 donors
- 61 IPF samples
- 40 NDC samples
- only samples with known lung region

The filtering design is:

    ~ diseasegroup + lunglocation

The donor identifier is not introduced as a fixed effect during filtering.
Within-donor correlation will be handled during downstream limma/voom model
fitting.

The frozen `filterByExpr()` settings are:

- `min.count = 10`
- `min.total.count = 15`
- `large.n = 10`
- `min.prop = 0.70`

Filtering result:

- input genes: 15,065
- retained genes: 15,012
- removed genes: 53
- retention: 99.65%

The relatively small number of removed genes is accepted and will not be
increased by introducing an arbitrary more aggressive threshold.

Audit of the removed genes confirmed that they are concentrated toward the
low-expression end of the dataset. Removed genes reached CPM >= 1 in only
7–19 primary-cohort samples, and their maximum-CPM distribution was
substantially lower than that of retained genes.

The filtering decision was made before normalization, PCA/MDS or formal
differential-expression testing.

---

## 2026-09-24 — W4 TMM normalization accepted

**Status:** ACCEPTED

TMM normalization is performed with `edgeR::calcNormFactors()` after the
frozen low-count filtering step.

Normalization factors are calculated independently from the retained raw
counts. GEO-reported normalization factors are preserved only as source
metadata and are not reused for the primary analysis.

Across the 101 primary-cohort samples, TMM factors were:

- minimum: approximately 0.417
- first quartile: approximately 0.976
- median: approximately 1.014
- mean: approximately 1.005
- third quartile: approximately 1.048
- maximum: approximately 1.154

The product of normalization factors was 1, as expected after scaling.

The lowest TMM factor belongs to `ALF009B` (NDC, Base), which was already
included in the W3 technical watchlist.

A targeted composition audit showed that `ALF009B` has exceptionally strong
library composition bias:

- most abundant gene: approximately 12.1% of the library
- 10 most abundant genes: approximately 60.0% of the library
- only 7 genes are required to account for 50% of all counts

Across all 101 samples, the corresponding median values were approximately:

- top-1 gene fraction: 1.68%
- top-10 gene fraction: 8.51%
- genes required to account for 50% of counts: 915

TMM factors were negatively associated with library concentration:

- Spearman TMM vs top-1 fraction: approximately -0.44
- Spearman TMM vs top-10 fraction: approximately -0.60

The extreme factor is therefore consistent with a strong composition effect
rather than, by itself, evidence of sample failure.

`ALF009B` remains in the analysis and remains flagged for downstream
PCA/MDS and model-diagnostic review.

No sample is excluded based solely on its TMM normalization factor.

---

## 2026-09-24 — W4 filtering, normalization and ordination decisions

**Status:** ACCEPTED / FROZEN before differential-expression testing

### Low-count filtering

The primary IPF-vs-NDC cohort was filtered using `edgeR::filterByExpr()` with the pre-specified disease-group and lung-location design.

The retained expression matrix contains:

- 15,012 genes
- 101 samples
- 34 donors

The filtering decision was made before differential-expression testing.

### Normalization

TMM normalization was performed independently from the raw counts using edgeR.

GEO-reported normalization factors were retained only for provenance and were not reused.

No sample was excluded on the basis of its TMM normalization factor.

### PCA / MDS quality-control assessment

PCA of TMM-normalized logCPM values showed:

- PC1: 31.11% variance explained
- PC2: 13.95% variance explained

PCA and edgeR MDS showed concordant large-scale structure.

Disease group was strongly associated with donor-level ordination:

- PC1 disease-group R2: 0.773
- PC2 disease-group R2: 0.594

No sample was excluded on the basis of PCA or MDS position.

Technical-watchlist samples remain included.

### Age assessment

Age showed measurable association with ordination:

- PC1 age-only R2: 0.213
- PC2 age-only R2: 0.298

However, disease group remained strongly associated with ordination after accounting for age:

- PC1 partial R2 for disease after age: 0.702
- PC2 partial R2 for disease after age: 0.426

The primary model will therefore retain the pre-specified disease and lung-region structure.

Age will be evaluated in a pre-specified secondary sensitivity analysis.

This avoids redefining the primary analysis after exploratory ordination and avoids excluding an otherwise valid donor solely because age is unavailable.

### Lung region

Lung region showed substantial within-donor association with ordination:

- PC1 partial R2: 0.237
- PC2 partial R2: 0.212

Lung region therefore remains a required adjustment variable in the primary model.

### RIN

Within-donor RIN association with the first two PCs was negligible:

- PC1 partial R2: 0.006
- PC2 partial R2: 0.001

RIN will not be added to the primary model solely on the basis of this QC screen.

### Processing date

Processing date showed little sample-level association with the first two PCs:

- PC1 R2: 0.009
- PC2 R2: 0.054

A within-donor processing-date effect is not estimable because processing date does not provide sufficient within-donor variation.

Processing date will not be included in the primary model based on the current evidence.

### Frozen primary analysis framework

Primary biological contrast:

`IPF vs NDC`

Primary fixed-effects structure:

`~ diseasegroup + lunglocation`

Repeated samples from the same donor will be handled explicitly through donor-aware correlation modelling.

Pre-specified sensitivity model:

`~ diseasegroup + lunglocation + age`

No differential-expression results were inspected when these decisions were made.
