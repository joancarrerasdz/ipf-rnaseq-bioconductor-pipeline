# Study Design

## Study

**GEO accession:** GSE213001
**Data type:** Human bulk RNA-seq
**Disease context:** Idiopathic pulmonary fibrosis (IPF)

This document pre-specifies the main analytical questions, cohorts and experimental-design principles before differential-expression testing.

---

## Experimental structure

The complete GEO dataset contains:

| Disease group | Samples | Unique donors |
|---|---:|---:|
| IPF | 62 | 20 |
| NDC | 41 | 14 |
| ILD | 26 | 9 |
| CLAD | 10 | 3 |
| **Total** | **139** | **46** |

Repeated tissue samples are present for most donors.

- 44 of 46 donors have more than one sample.
- Donors contribute between 1 and 5 samples.
- Samples may originate from different anatomical lung regions.
- Samples from the same donor must therefore not be treated as independent biological replicates.

The donor identifier `donorid` defines the repeated-measures / blocking unit.

---

# Primary analysis

## Scientific question

**Which transcriptional programs distinguish idiopathic pulmonary fibrosis from non-diseased lung tissue while accounting for anatomical lung region and repeated samples from the same donor?**

## Primary contrast

**IPF vs NDC**

The disease variable used to define the primary contrast is:

`diseasegroup`

This field was selected because it is internally consistent with `diagnosis` and `diseasegroup1` across all samples and provides the clearest disease-level grouping.

## Primary cohort

Before considering lung-region completeness:

- IPF: 62 samples / 20 donors
- NDC: 41 samples / 14 donors
- Total: 103 samples / 34 donors

Two samples have missing lung-region information:

- `ALF018E` / `GSM6568369` / IPF
- `ALF026E` / `GSM6568412` / NDC

Because anatomical lung region is a pre-specified adjustment variable in the primary analysis, these two samples will not be included in the primary model unless the design is explicitly revised before differential-expression testing.

The planned primary analysis cohort is therefore:

- **101 samples**
- **34 donors**
- IPF: 61 samples
- NDC: 40 samples

Region distribution:

| Disease group | Apex | Base | Total |
|---|---:|---:|---:|
| IPF | 30 | 31 | 61 |
| NDC | 22 | 18 | 40 |
| **Total** | **52** | **49** | **101** |

## Planned model structure

The primary model will estimate the IPF-versus-NDC disease effect while accounting for lung region and within-donor correlation.

Conceptually:

    expression ~ diseasegroup + lunglocation

with:

    block = donorid

or an equivalent repeated-measures / random-effect formulation.

The exact implementation will be finalized after filtering, normalization and QC diagnostics.

The current preferred framework is a donor-aware **limma-voom** analysis using an appropriate correlation/blocking strategy or a statistically equivalent formulation.

`donorid` must not be added as a conventional fixed-effect term alongside `diseasegroup`, because disease group is defined at donor level and would be confounded with donor-specific fixed effects.

---

# Secondary analysis

## Scientific question

**Are transcriptional differences associated with anatomical lung region within IPF?**

## Secondary contrast

**Apex vs Base within IPF**

The analysis will be restricted to IPF samples with known lung region and will explicitly account for repeated observations from the same donor.

Region coverage within IPF:

- 19 of 20 IPF donors have both Apex and Base samples.
- 1 IPF donor (`ALF001`) has only a Base sample.
- Some donors contribute multiple samples within the same anatomical region.

Therefore, this is not a simple one-to-one paired-sample design.

A repeated-measures / donor-aware model is required.

---

# Other disease groups

## ILD

The 26 ILD samples from 9 donors will be retained in the complete metadata but are not part of the primary contrast.

Any future ILD analysis is exploratory and will only be added after the primary and secondary analyses are complete.

## CLAD

The GEO Series Matrix contains 10 CLAD samples from 3 donors.

These samples are retained for provenance and metadata completeness but are outside the scope of the primary IPF question.

They will not be included in the primary or secondary analyses.

---

# Covariates

## Pre-specified

The following variables have a defined role before differential-expression testing:

- `diseasegroup`: primary exposure
- `lunglocation`: anatomical-region covariate / secondary exposure
- `donorid`: repeated-measures blocking unit

## Candidate covariates

The metadata also contain variables such as:

- age
- gender
- smoking status
- RIN
- disease severity
- pulmonary-function measurements

These variables will be examined during metadata/QC assessment.

They will not be added automatically to the primary model.

Additional covariates will only be included when scientifically and statistically justified, considering missingness, confounding, collinearity and available effective sample size.

---

# Differential-expression framework

The project will use one principal differential-expression framework rather than benchmarking multiple packages.

Current planned workflow:

    raw counts
    -> low-count filtering
    -> normalization
    -> voom transformation
    -> donor-aware model
    -> primary contrast
    -> multiple-testing correction

Differential expression will use false-discovery-rate control.

The precise implementation will be frozen before the first formal DE test.

---

# Pre-specified exclusions

At the current stage:

1. CLAD samples are outside the primary scientific scope.
2. ILD samples are outside the primary IPF-vs-NDC contrast.
3. Samples with unknown lung region are excluded from the planned primary model because lung region is a pre-specified adjustment variable.

No sample will be excluded based on PCA position, differential-expression results or improvement of visual group separation.

Any additional QC exclusion must be justified independently of the desired biological result and recorded in the decision log.

---

# Analysis principles

1. The experimental unit is the donor, not an individual tissue sample.
2. Repeated samples from the same donor are correlated observations.
3. Raw source files remain unchanged.
4. Sample metadata are aligned explicitly to raw-count columns by `sample_title`.
5. Filtering and normalization will be scripted.
6. QC exclusions must be documented before formal differential-expression testing.
7. One principal DE framework will be used.
8. Primary and secondary contrasts are pre-specified before inspecting DE results.
9. Exploratory analyses will be clearly labeled as exploratory.
10. Biological interpretation will distinguish association from causal or clinical claims.

---

## Status

**Design pre-specification — W2**

The cohort structure and scientific contrasts are defined.

The final statistical implementation remains conditional on completion of metadata QC, count QC, filtering and normalization.
