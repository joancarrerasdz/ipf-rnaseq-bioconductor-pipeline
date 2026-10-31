# P2-16 reproducibility evidence

**Project:** GSE213001 IPF RNA-seq / Bioconductor
**QA gate:** P2-16 — downstream reproducibility
**Validation date:** 2026-10-31
**Status:** PASS
**Validated integration commit:** `2f74c5be92a414db16e7ea1eb8dfa424cfe06e4d`

## Purpose

This record documents the clean-room validation performed to verify that the
versioned `targets` graph and `renv` environment reproduce the complete
downstream workflow from the public GSE213001 source inputs through the final
W8 figures and report.

The validation was performed without changing the frozen scientific
interpretation or the historical W6-W8 analytical checkpoints.

## Environment

- R: 4.5.1
- Quarto: 1.10.20
- Dependency management: `renv`
- Workflow orchestration: `targets`
- Repository branch tested: `qa/p2-16-reproducibility`
- Fresh-clone validation: yes

## Public source inputs

The raw GEO source files are intentionally not versioned in Git and must be
provisioned according to `data/README.md`.

Validated input hashes:

- `GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz`
  - SHA-256: `832ba7db504d0fb34405833b5c271b364b9c663970974b70d19a0d78232e3bef`
- `GSE213001_series_matrix.txt.gz`
  - SHA-256: `b76293ab154b4adc21cd6cd840b71323f7849c586bea9ca1efff2f0053af71c8`

## Historical checkpoints preserved

- `w6-complete`: `3534b6c9b4f500ed470c399ce4433b330df65c54`
- `w7-complete`: `46173b04bec2383e8b8451f73660408794d09027`
- `w8-complete`: `c44bc0d6e2996cc15ce1e4cebc37066bfd83d701`

These checkpoints remained unchanged throughout P2-16 remediation and
clean-room validation.

## Clean-room procedure

A fresh clone of the QA branch was created outside the working repository.
Only the two public GEO inputs listed above were provisioned manually.

The R environment was restored from the committed lockfile:

```r
source("renv/activate.R")
renv::restore(prompt = FALSE)
```

The complete dependency graph was then rebuilt:

```r
library(targets)
tar_make()
```

## Clean-room results

The complete W1-W8 graph rebuilt successfully:

- total targets rebuilt: **65/65**
- skipped targets: **0**
- outdated targets after build: **0**
- W1-W8 target graph current: **PASS**

### Scientific outputs

Tracked scientific CSV and Markdown outputs were compared with the frozen
baseline after reconstruction.

Result:

- scientific CSV/MD outputs: **byte-identical**
- frozen scientific interpretation: **unchanged**

### HTML report

The final HTML report was compared using visible document content rather than
renderer-generated HTML infrastructure.

Results:

- visible scientific content: **exact**
- heading structure: **15/15 identical**
- figure semantic descriptions: **3/3 identical**
- required scientific terminology: **preserved**

Full-file HTML hashes may differ because Quarto embeds renderer-dependent
infrastructure that is not part of the scientific content.

### PDF report

The reconstructed PDF was assessed semantically and structurally.

Results:

- pages: **10/10 identical**
- page geometry: **identical**
- semantic tokens: **1976/1976 identical**
- numeric tokens: **133/133 identical**
- scientific content equivalence: **exact**

Binary or pixel-level PDF differences are classified as renderer-only drift
when the semantic, numerical and structural scientific content is preserved.

### Final figures

Final PNG dimensions were preserved:

- `W8_primary_IPF_vs_NDC_volcano.png`: **3000 × 2100**
- `W8_biological_interpretation_matrix.png`: **4800 × 2700**
- `W8_robustness_summary.png`: **4800 × 2700**

Renderer-dependent binary differences are permitted only when scientific
content, geometry, labels, data and interpretation remain unchanged.

## Dependency invalidation test

A temporary non-destructive modification was introduced into the raw-count
source file.

Without executing analytical targets, `tar_outdated()` detected downstream
invalidation across the workflow.

Observed result:

- outdated targets after temporary input change: **32**
- W5 primary DE invalidated: **PASS**
- W6 downstream interpretation invalidated: **PASS**
- W7 sensitivity and robustness invalidated: **PASS**
- W8 reporting invalidated: **PASS**

No `tar_make()` was executed against the modified source input.

The raw input was subsequently restored byte-for-byte. After exact
restoration:

- outdated targets: **0**
- target graph returned to current state: **PASS**

## P2-16 decision

The technical reproducibility requirements of P2-16 are satisfied.

The repository now demonstrates that:

1. the committed `renv` environment can be restored;
2. the committed `targets` graph reconstructs W1-W8 from a fresh clone;
3. upstream input changes invalidate dependent downstream targets;
4. the final report is reconstructed through the dependency graph;
5. structured scientific results reproduce byte-for-byte;
6. renderer-only presentation drift does not alter scientific content;
7. the original frozen scientific results and historical checkpoints remain
   unchanged.

**P2-16 technical reproducibility status: PASS.**
