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

# Open decisions — historical register

The list originally created during study-design pre-specification is
superseded by the dated decisions below.

Low-count filtering, normalization, QC-based sample retention and the
primary covariate framework were subsequently resolved and frozen during W4.

Current unresolved decisions are listed at the end of this document.

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

---

## 2026-09-27 — Current open decisions after W4

**Status:** ACTIVE REGISTER

The following methodological decisions remain intentionally unresolved
after completion of W4:

- final donor-correlation / repeated-measures implementation for formal
  differential-expression modelling
- exact Ensembl-to-gene annotation strategy and annotation version
- handling of unmapped or ambiguous gene identifiers
- enrichment gene universe
- exact GSEA and/or over-representation workflow

The following items are no longer open:

- low-count filtering rule — frozen on 2026-09-24
- normalization method and parameters — TMM accepted on 2026-09-24
- W4 QC-based sample exclusions — no samples excluded
- primary fixed-effects structure — frozen as
  `~ diseasegroup + lunglocation`
- age — pre-specified sensitivity model
- RIN and processing date — not added to the primary model based on W4 QC

No formal differential-expression results had been inspected when these
decisions were frozen.

---

## 2026-10-02 — W5 donor-aware differential-expression implementation

**Status:** ACCEPTED / FROZEN BEFORE FORMAL DE

The primary differential-expression framework was frozen before inspection of
formal gene-level DE results.

### Primary analysis

- biological contrast: **IPF vs NDC**
- primary cohort: **101 samples / 34 donors**
- retained genes after the frozen W4 low-count filter: **15,012**
- normalization: **TMM**
- fixed-effects design: `~ diseasegroup + lunglocation`
- reference disease group: **NDC**
- reference lung region: **Base**
- primary coefficient: `diseasegroupIPF`
- repeated-measures blocking unit: `donorid`

Repeated samples from the same donor will be handled using the
`limma-voom` / `duplicateCorrelation` framework.

A two-pass donor-correlation dry run produced:

- first consensus correlation: **0.315821**
- second consensus correlation: **0.315828**
- absolute change: approximately **6.46e-06**

The near-identical estimates support a stable consensus within-donor
correlation estimate.

The final primary model will therefore use the donor-aware voom object and
the second consensus correlation estimate for formal model fitting.

### Pre-specified sensitivity analysis

Age remains reserved for the pre-specified sensitivity model:

`~ diseasegroup + lunglocation + age`

The single donor with missing age is `ALF017`, corresponding to two samples.
The age-complete sensitivity cohort therefore contains:

- **99 samples**
- **33 donors**

RIN and processing date are not added to the primary model, consistent with
the W4 QC decisions.

No formal `lmFit`, empirical-Bayes moderation, multiple-testing results or
gene-level differential-expression tables had been inspected when this
implementation was frozen.

---

## 2026-10-02 — W5 formal-inference specification

**Status:** ACCEPTED / FROZEN BEFORE FIRST FORMAL DE TEST

The exact inferential settings for the primary IPF-vs-NDC analysis were
specified before execution of `lmFit()` or `eBayes()`.

### Model fitting

The final donor-aware voom object from the two-pass correlation workflow
will be fitted with:

`lmFit(v2, design, block = donorid, correlation = rho2)`

where `rho2` is the second consensus correlation estimate obtained from
`duplicateCorrelation()`.

The frozen fixed-effects design is:

`~ diseasegroup + lunglocation`

The primary coefficient is:

`diseasegroupIPF`

Therefore:

- positive logFC = higher expression in IPF relative to NDC
- negative logFC = lower expression in IPF relative to NDC

with lung region adjusted in the model.

### Empirical-Bayes moderation

Empirical-Bayes moderation will use:

`eBayes(fit, robust = TRUE, trend = FALSE)`

Robust moderation is pre-specified to reduce sensitivity to genes with
atypical residual variances. Mean-variance dependence is already handled
through voom weights, therefore `trend = FALSE`.

### Multiple testing

P-values for the primary coefficient will be corrected using the
Benjamini-Hochberg procedure.

The primary multiplicity criterion is:

`FDR < 0.05`

No minimum absolute log-fold-change threshold is used to define primary
statistical significance. Effect sizes will be reported separately and
will not be selected post hoc to alter significance calls.

`treat()` is not part of the primary analysis.

### Reporting order

The complete gene-level result table will retain all 15,012 tested genes.
Initial validation will focus on model-level and multiplicity-level
summaries before biological interpretation of individual genes.

No formal gene-level DE results had been generated when these settings
were frozen.

---

## 2026-10-02 — W5 primary donor-aware DE execution

**Status:** COMPLETED / POST-INFERENCE RECORD

The pre-specified primary donor-aware differential-expression model was
executed after the analytical framework and inferential settings had been
frozen and committed.

### Primary model

- contrast: **IPF vs NDC**
- tested genes: **15,012**
- samples: **101**
- donors: **34**
- fixed effects: `~ diseasegroup + lunglocation`
- repeated-measures block: `donorid`
- final consensus within-donor correlation: **0.3158276**
- empirical-Bayes moderation: `robust = TRUE, trend = FALSE`
- multiple-testing adjustment: Benjamini-Hochberg

### Global inferential results

- raw P < 0.05: **8,679 genes**
- FDR < 0.05: **7,898 genes**
- FDR < 0.05 with positive logFC: **3,964**
- FDR < 0.05 with negative logFC: **3,934**
- |logFC| >= 1: **1,884 genes**
- FDR < 0.05 and |logFC| >= 1: **1,861 genes**
- minimum raw P: **1.099e-28**
- minimum FDR: **1.649e-24**

The large number of FDR-significant genes indicates broad transcriptomic
separation between the primary disease groups. Biological interpretation
of individual genes is intentionally deferred until completion of the
pre-specified age sensitivity analysis.

No gene-level ranking or biological interpretation was used to modify the
primary analysis.

---

## 2026-10-02 — W5 age-sensitivity implementation

**Status:** ACCEPTED / FROZEN BEFORE AGE-ADJUSTED RESULTS

The pre-specified age sensitivity analysis will be performed on the
age-complete subset of the primary cohort.

### Sensitivity cohort

The donor with missing age is `ALF017`, corresponding to two samples.

The sensitivity cohort therefore contains:

- **99 samples**
- **33 donors**

The frozen W4 low-count filtering decision will be retained unchanged:

- genes evaluated: **15,065**
- genes retained: **15,012**
- genes removed: **53**

No gene filtering rule will be re-estimated using the sensitivity results.

### Same-cohort comparison

To separate the effect of age adjustment from the effect of removing the
two samples with missing age, two models will be fitted to exactly the same
99-sample cohort.

Reduced sensitivity model:

`~ diseasegroup + lunglocation`

Age-adjusted sensitivity model:

`~ diseasegroup + lunglocation + age`

The primary coefficient in both models remains:

`diseasegroupIPF`

### Normalization and repeated measurements

TMM normalization will be recalculated within the 99-sample sensitivity
cohort.

The within-donor consensus correlation will also be estimated separately
for each sensitivity design using the same two-pass
`voom` / `duplicateCorrelation` procedure used for the primary analysis.

Repeated measurements will continue to use:

`block = donorid`

### Inference

Both sensitivity models will use:

`eBayes(..., robust = TRUE, trend = FALSE)`

and Benjamini-Hochberg multiple-testing correction.

The significance criterion remains:

`FDR < 0.05`

No absolute logFC threshold will be introduced as a significance criterion.

### Concordance assessment

Age robustness will be assessed using the same 15,012 genes and will include:

- Pearson correlation of disease logFC estimates
- Spearman correlation of disease logFC estimates
- concordance of logFC direction
- median absolute change in disease logFC
- number of FDR-significant genes in each model
- overlap of FDR-significant genes
- genes significant only before age adjustment
- genes significant only after age adjustment
- concordance among genes with |logFC| >= 1

The comparison of the reduced and age-adjusted models on the identical
99-sample cohort will be the primary assessment of the impact of age.

No age-adjusted gene-level results had been generated when this sensitivity
framework was frozen.

---

## 2026-10-02 — W5 age-sensitivity results

**Status:** COMPLETED / POST-INFERENCE SENSITIVITY RECORD

The pre-specified age sensitivity analysis was executed on the identical
age-complete primary cohort.

### Sensitivity cohort

- samples: **99**
- donors: **33**
- missing-age donor removed: `ALF017`
- frozen tested gene universe: **15,012 genes**

TMM normalization and within-donor consensus correlation were re-estimated
within this cohort.

### Same-cohort models

Reduced model:

`~ diseasegroup + lunglocation`

Age-adjusted model:

`~ diseasegroup + lunglocation + age`

The coefficient of interest remained `diseasegroupIPF` in both models.

Final within-donor consensus correlations were:

- reduced model: **0.3082255**
- age-adjusted model: **0.3089274**

### Effect-estimate concordance

Comparison of the disease coefficient across all 15,012 genes produced:

- Pearson logFC correlation: **0.975307**
- Spearman logFC correlation: **0.959262**
- direction concordance: **92.14%**
- median absolute logFC change: **0.061244**

These results indicate strong overall stability of the estimated IPF-vs-NDC
effect after adjustment for age.

### Multiple-testing sensitivity

FDR < 0.05:

- reduced model: **7,662 genes**
- age-adjusted model: **5,357 genes**
- significant in both models: **5,131 genes**
- significant only before age adjustment: **2,531 genes**
- significant only after age adjustment: **226 genes**
- FDR-set Jaccard index: **0.650482**

Age adjustment therefore materially reduced the number of genes meeting
the FDR threshold, while most genes significant after age adjustment were
also significant in the reduced model.

### Large-effect concordance

For genes with |logFC| >= 1 in either sensitivity model:

- genes: **2,068**
- direction concordance: **99.90%**

For genes with |logFC| >= 1 in both models:

- genes: **1,651**
- direction concordance: **100.00%**

Thus, the largest estimated disease effects are highly stable in direction
after age adjustment.

### Interpretation for downstream analysis

Age is retained as an important pre-specified sensitivity covariate because
it materially affects statistical significance for a substantial subset of
genes.

However, the high genome-wide logFC correlations, small median effect-size
change and near-perfect direction concordance among large effects indicate
that the broad IPF-vs-NDC transcriptomic signal is not explained solely by
age imbalance.

The frozen primary analysis remains the primary model:

`~ diseasegroup + lunglocation`

The age-adjusted model remains a pre-specified sensitivity analysis and
will be used when assessing robustness of downstream biological findings.

No gene-level biological interpretation was used to alter either model.

---

## 2026-10-05 — W6 gene-annotation strategy

**Status:** ACCEPTED / FROZEN BEFORE GENE MAPPING

The gene-annotation strategy was specified before inspection of gene symbols,
gene names or biological interpretation of the W5 differential-expression
results.

### Frozen analytical universe

The annotation layer will be applied to the frozen W5 primary
differential-expression universe:

- **15,012 genes**
- **15,012 unique Ensembl gene IDs**
- **0 duplicated Ensembl gene IDs**
- all identifiers are unversioned `ENSG...` identifiers

Annotation will not alter the statistical universe used for differential
expression.

### Annotation source and software

Gene annotation will use:

- Bioconductor: **3.21**
- `AnnotationDbi`: **1.70.0**
- `org.Hs.eg.db`: **3.21.0**
- input key type: `ENSEMBL`

The annotation package versions are recorded in `renv.lock`.

### Annotation fields

The primary annotation fields will be:

- `SYMBOL`
- `ENTREZID`
- `GENENAME`

The original Ensembl identifier will remain the primary analytical identifier
and will be retained as `gene_id`.

### Mapping rules

Annotation is a downstream descriptive layer and must not modify the frozen
W5 differential-expression results.

The following rules are frozen:

1. Every original Ensembl gene ID remains represented in the annotated output.
2. Genes without an annotation mapping will be retained with missing
   annotation fields.
3. An unmapped gene will not be excluded from the DE universe.
4. Multiple mappings will be detected and audited explicitly.
5. Mapping will not silently duplicate analytical rows.
6. Different Ensembl gene IDs mapping to the same gene symbol will not be
   aggregated solely on the basis of the shared symbol.
7. Gene symbols will not replace Ensembl IDs as the primary identifier.
8. Annotation status will be reported quantitatively before biological
   interpretation begins.

### Separation from inference

The following W5 results remain frozen and unchanged:

- primary IPF-vs-NDC contrast
- donor-aware correlation modelling
- lung-region adjustment
- TMM normalization
- low-count filtering
- empirical-Bayes settings
- Benjamini-Hochberg multiple-testing correction
- age sensitivity analysis

No gene identity, annotation result, biological function or pathway information
was inspected when these annotation rules were frozen.


---

## 2026-10-05 — W6.3 gene-annotation mapping audit

**Status:** COMPLETED / ACCEPTED

The frozen W5 differential-expression universe of 15,012 unique Ensembl
gene identifiers was mapped using the annotation strategy frozen before
inspection of annotation results.

### Mapping coverage

Using:

- key type: `ENSEMBL`
- `AnnotationDbi` 1.70.0
- `org.Hs.eg.db` 3.21.0

the annotation audit produced:

- input genes: **15,012**
- unique input genes: **15,012**
- genes with a single mapping record: **14,827**
- genes with multiple mapping records: **112**
- unmapped genes: **73**

The 73 unmapped genes had no available `SYMBOL`, `ENTREZID`, or
`GENENAME` under the frozen annotation source.

### Mapping ambiguity

The mapping audit identified:

- genes with multiple `SYMBOL` values: **112**
- genes with multiple `ENTREZID` values: **112**
- genes with multiple `GENENAME` values: **111**
- gene symbols shared by multiple Ensembl gene IDs: **16**

These cases will not be resolved by selecting a mapping post hoc.

For the one-row-per-Ensembl analytical annotation table, multiple annotation
values will be retained explicitly rather than silently duplicating DE rows.

### Analytical invariants

The annotation process did not modify the W5 statistical universe:

- **15,012 genes remain represented**
- no gene was excluded because annotation was missing
- no gene was excluded because annotation was ambiguous
- Ensembl `gene_id` remains the primary identifier
- shared gene symbols will not trigger aggregation
- the ordering and inferential results of the frozen W5 DE table remain
  unchanged

The annotation mapping was audited quantitatively before biological
interpretation of individual genes.


---

## 2026-10-05 — W6.4 annotated differential-expression table

**Status:** COMPLETED / ACCEPTED

The frozen W5 primary differential-expression table was joined to the
frozen W6 annotation layer using exact Ensembl `gene_id` matching.

### Structural integrity

The join was explicitly required to preserve the W5 inferential result.

Audit results:

- input DE rows: **15,012**
- annotation rows: **15,012**
- output annotated rows: **15,012**
- unique output Ensembl gene IDs: **15,012**
- genes missing from the annotation table: **0**
- annotation-only genes: **0**
- gene order preserved: **TRUE**
- all original W5 columns preserved unchanged: **TRUE**

The annotated table therefore represents a descriptive annotation layer
over the frozen W5 result and not a new statistical analysis.

### Annotation status

Among the 15,012 tested genes:

- **14,827** have a single annotation mapping record
- **112** have multiple annotation mapping records
- **73** are unmapped
- **73** have no available gene symbol

No gene was removed, aggregated or reordered because of annotation status.

### Frozen inferential boundary

The following remain unchanged:

- W5 differential-expression universe
- IPF-vs-NDC contrast
- lung-region adjustment
- donor-aware correlation modelling
- TMM normalization
- empirical-Bayes settings
- Benjamini-Hochberg multiple-testing correction
- age sensitivity analysis
- all gene-level W5 statistical values

No biological interpretation or gene prioritization was used to construct
or modify the annotated table.


---

## 2026-10-05 — W6.5 biological interpretation and enrichment framework

**Status:** ACCEPTED / FROZEN BEFORE BIOLOGICAL INTERPRETATION

The biological-interpretation framework was specified after completion of
the annotation integrity audit but before inspection or ranking of individual
gene identities.

The frozen W5 primary differential-expression analysis remains the sole
primary inferential analysis.

### Primary statistical evidence

Primary statistical significance remains defined as:

`BH FDR < 0.05`

No absolute log-fold-change threshold is introduced as an additional
statistical significance criterion.

Effect size and statistical significance will therefore be reported
separately.

For descriptive biological prioritization, an absolute effect size of:

`|logFC| >= 1`

may be used to identify genes with relatively large estimated effects.

This threshold is descriptive and must not alter the primary FDR-based
significance calls.

### Direction of effect

The frozen primary coefficient is:

`diseasegroupIPF`

Therefore:

- positive `logFC` = higher expression in IPF relative to NDC
- negative `logFC` = lower expression in IPF relative to NDC

All biological interpretation must preserve this direction convention.

### Gene-level interpretation

The original Ensembl `gene_id` remains the primary analytical identifier.

`SYMBOL`, `ENTREZID`, and `GENENAME` are descriptive annotation fields.

Genes lacking annotation remain in the statistical result and must not be
removed from the W5 DE table.

Genes with ambiguous annotation will be flagged rather than resolved by
post-hoc selection of a preferred biological label.

Shared symbols must not trigger aggregation of differential-expression rows.

### Gene prioritization

Gene prioritization will only begin after this framework is frozen.

The primary descriptive priority set will consist of genes satisfying:

`BH FDR < 0.05` and `|logFC| >= 1`

Genes meeting `BH FDR < 0.05` but with smaller effect sizes remain valid
statistically significant findings and will not be reclassified as
non-significant.

Gene ranking for reporting will use frozen statistical quantities and will
not be used to redefine the statistical analysis.

### Age-sensitivity interpretation

Age remains a pre-specified sensitivity analysis and does not replace the
primary model.

Robustness to age adjustment will be evaluated using the previously generated
99-sample same-cohort comparison between:

Reduced model:

`~ diseasegroup + lunglocation`

and age-adjusted model:

`~ diseasegroup + lunglocation + age`

Gene-level age robustness will be treated descriptively.

Genes with concordant effect direction and FDR < 0.05 in both same-cohort
models will be considered strongly robust to age adjustment.

Genes whose significance materially changes after age adjustment will be
flagged as age-sensitive rather than removed from the primary W5 result.

### Mapping rules for enrichment

Pathway enrichment is downstream of the frozen differential-expression
analysis and must not modify the DE universe or statistical results.

For identifier-dependent enrichment analyses, only genes with an
unambiguous usable annotation identifier will contribute to the mapped
enrichment universe.

Genes excluded from a specific enrichment analysis because of missing or
ambiguous annotation will remain present in the primary DE table and their
number will be reported explicitly.

No arbitrary mapping will be selected from one-to-many annotation cases.

### Over-representation analysis

Over-representation analysis (ORA) will be performed separately for genes
with higher expression in IPF and genes with lower expression in IPF.

The descriptive foreground definition will be:

`BH FDR < 0.05` and `|logFC| >= 1`

The enrichment background will be derived from the full frozen set of
15,012 tested genes after application of the same identifier-mapping
eligibility rules used for the foreground.

The background must therefore represent tested genes, not the complete
genome.

Multiple-testing correction for enrichment results will use the
Benjamini-Hochberg procedure.

### Ranked gene-set analysis

Ranked enrichment will complement threshold-based ORA.

The ranking statistic will be the moderated disease-effect test statistic
from the frozen W5 model rather than raw p-value alone.

Positive ranking values represent enrichment toward higher expression in IPF;
negative ranking values represent enrichment toward lower expression in IPF.

No DE-significance threshold will be used to construct the ranked input.

Identifier ambiguity will be resolved only through deterministic,
pre-specified eligibility rules and never by selecting genes because they
produce stronger biological enrichment.

### Biological databases

The primary functional interpretation will focus on:

- Gene Ontology Biological Process
- Reactome pathways

Additional databases may be used only as explicitly labelled secondary or
exploratory analyses.

Enrichment databases and software versions must be recorded in the
reproducible environment before execution.

### Reporting principles

Biological interpretation will distinguish explicitly between:

1. statistical significance,
2. magnitude and direction of differential expression,
3. robustness to age adjustment,
4. annotation certainty,
5. pathway-level enrichment.

No individual gene or pathway will be selected to modify filtering,
normalization, model specification, covariate choice, or inferential
thresholds.

At the time this framework was frozen, individual gene identities had not
been inspected for biological prioritization and no enrichment analysis had
been performed.


---

## 2026-10-05 — W6.7 deterministic gene-reporting order

**Status:** ACCEPTED / FROZEN BEFORE FIRST GENE-IDENTITY READOUT

No individual gene identities had been inspected for biological
prioritisation when this reporting order was specified.

### Primary descriptive gene set

The primary descriptive set remains:

`BH FDR < 0.05` and `|logFC| >= 1`

This definition does not replace or modify the primary statistical
significance criterion of `BH FDR < 0.05`.

### Direction-specific reporting

Genes will be reported separately as:

- higher in IPF: `logFC >= 1`
- lower in IPF: `logFC <= -1`

### Deterministic ordering

Within each direction, genes will be ordered by:

1. `adj.P.Val` ascending,
2. absolute `logFC` descending,
3. absolute moderated `t` descending,
4. `gene_id` ascending as a deterministic final tie-breaker.

Annotation content, gene name, known biological function, pathway membership,
publication history or perceived biological relevance will not influence
this ordering.

### Annotation eligibility

All genes in the frozen descriptive set remain represented.

Genes with missing or ambiguous annotation will be flagged explicitly and
will not be silently removed from the gene-level statistical reporting.

A separate annotation-clean subset may be used for presentation when a
human-readable gene symbol is required, but the complete frozen statistical
set will remain available.

### Age robustness

Age sensitivity will be attached descriptively after the deterministic
primary ranking has been established.

Age robustness will not be used to reorder or redefine the primary ranked
list.

A gene will be labelled strongly age-robust when, in the pre-specified
99-sample sensitivity comparison:

- `FDR_reduced < 0.05`,
- `FDR_age_adjusted < 0.05`, and
- effect direction is concordant.

This label does not alter the frozen 101-sample primary DE result.


---

## 2026-10-05 — W6 functional enrichment framework

**Status:** ACCEPTED / FROZEN BEFORE PATHWAY INSPECTION

Functional enrichment will be performed only after completion of the
deterministic gene-level readout and annotation audit.

No pathway result had been inspected when this framework was frozen.

### Statistical separation

Functional enrichment is a downstream descriptive interpretation layer.

It will not modify:

- the frozen 15,012-gene W5 differential-expression universe
- the primary IPF-vs-NDC statistical model
- TMM normalization
- donor-aware correlation modelling
- lung-region adjustment
- empirical-Bayes settings
- Benjamini-Hochberg correction
- age sensitivity analysis
- gene-level ranking

### Identifier universe

Enrichment analyses requiring Entrez identifiers will use only genes with a
single, unambiguous Entrez mapping.

The frozen enrichment-eligible universe therefore contains:

- **14,827 genes**

Genes without an Entrez identifier or with ambiguous Entrez mapping remain in
the complete gene-level DE results but are not used for Entrez-based
enrichment.

### Directional over-representation analysis

Two independent foreground sets will be analysed:

- **Higher in IPF:** 1,372 genes
- **Lower in IPF:** 465 genes

Both sets are defined by:

- primary-model FDR < 0.05
- absolute log2 fold change >= 1
- single unambiguous Entrez mapping

The common ORA background is the complete 14,827-gene enrichment-eligible
universe.

Higher- and lower-in-IPF genes will not be combined for directional ORA.

### Ranked enrichment

Ranked enrichment will use all **14,827 enrichment-eligible genes**.

The ranking statistic will be the moderated t statistic from the frozen W5
primary differential-expression model.

Positive ranking values represent higher expression in IPF relative to NDC;
negative values represent lower expression in IPF relative to NDC.

No FDR or fold-change threshold will be used to construct the ranked list.

### Interpretation principles

Pathway interpretation will consider:

1. statistical significance,
2. enrichment magnitude,
3. direction,
4. gene-set size,
5. consistency between ORA and ranked enrichment,
6. robustness of contributing genes to the pre-specified age sensitivity
   analysis,
7. biological coherence across related pathways.

Individual pathways will not be used to redefine the gene universe or
statistical model.

Redundant pathways will be interpreted as biological themes rather than as
independent discoveries.

No enrichment results had been generated or inspected when these decisions
were frozen.
