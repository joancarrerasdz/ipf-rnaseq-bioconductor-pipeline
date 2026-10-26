# Data

This directory documents the public source data used by the project and the
boundary between source inputs, tracked metadata, and frozen analytical
outputs.

## Primary dataset

- **Study:** GSE213001
- **Repository:** NCBI Gene Expression Omnibus (GEO)
- **Organism:** Homo sapiens
- **Data type:** Bulk RNA-seq
- **Disease context:** Idiopathic pulmonary fibrosis (IPF)
- **Primary comparison:** IPF versus non-diseased control (NDC)

The study also contains non-IPF interstitial lung disease samples. The primary
portfolio analysis uses the prespecified IPF-versus-NDC comparison.

Multiple tissue samples can originate from the same donor and from different
anatomical lung regions. The primary model therefore accounts for lung region
and within-donor correlation rather than treating every sample as an
independent biological replicate.

## Source files used by the workflow

The tracked `_targets.R` workflow expects the following public GSE213001 files:

~~~text
data/source/GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz
data/source/GSE213001_series_matrix.txt.gz
~~~

The raw-count matrix provides gene-level expression counts.

The GEO Series Matrix provides study/sample metadata used during metadata
alignment.

Both source files originate from the public GSE213001 GEO record.

## Source-data policy

`data/source/` is the local landing area for public source files.

The source-data files themselves are intentionally excluded from Git. The
tracked `.gitkeep` preserves the expected directory structure.

Source files must remain unmodified after download.

Cleaning, transformation, filtering, normalization, and metadata alignment are
performed programmatically by the analytical workflow rather than by manually
editing the source files.

## Reconstructing the tracked pipeline foundation

After cloning the repository, restore the R environment from `renv.lock`:

~~~r
renv::restore()
~~~

Place the two expected GSE213001 source files in `data/source/`.

The tracked target graph can then be executed with:

~~~r
targets::tar_make()
~~~

The current `_targets.R` definition reconstructs the W1-W4 analytical
foundation:

1. source-file validation;
2. data import and metadata alignment;
3. metadata and library quality control;
4. gene filtering;
5. TMM normalization;
6. PCA/MDS and ordination quality control.

## Downstream analytical boundary

The completed portfolio continues beyond the current target graph.

W5-W8 contain the frozen primary differential-expression analysis, biological
interpretation, sensitivity and robustness analyses, final visualizations, and
final reporting.

Those downstream stages are preserved through versioned R scripts, frozen
result tables, QC metadata, figures, analytical decision records, and the final
HTML/PDF report.

The current repository does not claim that W5-W8 are rebuilt automatically by
a single `targets::tar_make()` call.

This boundary is documented explicitly so that reproducibility claims match
the actual repository implementation.

## Derived data and metadata

Tracked derived information is stored primarily under:

- `metadata/` — sample metadata, QC summaries, mappings, and audit records;
- `results/` — frozen analytical and interpretation outputs;
- `figures/` — analytical and final portfolio visualizations;
- `reports/` — final report source and rendered reports.

The public source files under `data/source/` remain separate from these tracked
derived outputs.

## Raw sequencing data

Raw sequencing reads are available through the public study archives, but full
FASTQ reprocessing is outside the scope of this portfolio workflow.

The primary project begins from the public gene-level raw-count matrix and
associated GEO sample metadata.

## Reproducibility principle

The repository separates:

- immutable public source inputs;
- scripted analytical transformations;
- tracked QC and provenance records;
- frozen inferential checkpoints;
- presentation-oriented final reporting.

This separation allows the analytical history and final portfolio conclusions
to be audited without silently redefining previously frozen results.
