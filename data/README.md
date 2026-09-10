# Data

This directory documents the data sources used in the project and how they are acquired.

## Primary dataset

- **Study:** GSE213001
- **Repository:** NCBI Gene Expression Omnibus (GEO)
- **Organism:** Homo sapiens
- **Data type:** Bulk RNA-seq
- **Disease context:** Idiopathic pulmonary fibrosis (IPF)
- **GEO accession:** GSE213001

GEO series page:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213001

## Study groups

The study contains samples from:

- idiopathic pulmonary fibrosis (IPF)
- non-diseased control lungs
- non-IPF interstitial lung disease (ILD)

Multiple tissue samples may originate from the same donor and from different anatomical lung regions.

Therefore, sample-level observations cannot automatically be treated as independent biological replicates.

## Analysis data

The complete downstream analysis will use the public raw-count matrix and associated sample metadata distributed through GEO.

The analysis will preserve the original source files unchanged.

Derived datasets, cleaned metadata and analysis-ready objects will be generated reproducibly by scripts in this repository.

## Raw sequencing data

Raw sequencing data are available through the Sequence Read Archive (SRA).

Full FASTQ reprocessing is not required for the primary downstream analysis.

A small subset of samples may later be used in `fastq-demo/` to demonstrate:

FASTQ
→ FastQC
→ MultiQC
→ optional trimming if justified
→ alignment / quantification
→ gene counts

## Repository policy

Raw/source files are stored locally under:

`data/source/`

and are excluded from Git through `.gitignore`.

Source files must not be manually edited.

Any cleaning, transformation or filtering must be performed programmatically and documented in the analysis workflow.

## Provenance record

For every downloaded source file, the following information will be recorded before analysis:

- original filename
- source repository
- accession
- source URL
- download date
- file size
- checksum where appropriate
- role in the analysis

Exact filenames and download commands will be added after the GEO supplementary files and metadata have been inspected.

## Reproducibility

The goal is for a fresh clone of the repository to be able to reconstruct the analysis inputs from documented public sources without relying on undocumented local files.
