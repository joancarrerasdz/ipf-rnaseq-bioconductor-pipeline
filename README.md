# Reproducible RNA-seq Analysis of Idiopathic Pulmonary Fibrosis

A reproducible bulk RNA-seq analysis of idiopathic pulmonary fibrosis (IPF) using public data from GEO accession **GSE213001**.

## Project objective

The primary objective is to identify transcriptional programs that distinguish idiopathic pulmonary fibrosis from non-diseased lung tissue while accounting for repeated samples from the same donor.

A secondary analysis will investigate transcriptional variation associated with anatomical lung region within IPF.

## Dataset

- **GEO accession:** GSE213001
- **Data type:** Human bulk RNA-seq
- **Disease:** Idiopathic pulmonary fibrosis (IPF)
- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Analysis input:** Public raw count matrix and sample metadata
- **Raw sequencing data:** Available through SRA; a small subset may later be used for an end-to-end FASTQ demonstration

Raw/source data are not committed directly to this repository. Data provenance and acquisition instructions will be documented in `data/README.md`.

## Core design principle

Multiple tissue samples may originate from the same donor. These observations must therefore not be treated as independent biological replicates.

The statistical design will explicitly account for donor-level dependence and repeated measurements.

## Planned workflow

1. Data acquisition and provenance
2. Metadata validation
3. Quality control
4. Low-count filtering
5. Normalization
6. PCA / MDS and sample relationships
7. Donor-aware experimental design
8. Differential expression
9. Multiple-testing correction
10. Functional enrichment
11. Pathway analysis
12. Biological interpretation
13. Reproducible Quarto report

## Reproducibility

The project will use:

- R
- Bioconductor
- `renv`
- `targets`
- Quarto
- Git / GitHub

Package versions and the executable workflow will be added progressively as the analysis is implemented.

## Status

**Development / setup phase**

Current focus: dataset provenance, metadata structure, experimental design and repository initialization.

No differential-expression results are reported at this stage.
