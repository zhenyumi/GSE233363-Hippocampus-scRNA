# R/scrna/ — scRNA-seq Reproduction Scripts

This directory contains the scRNA-seq ferroDG stage-2 reproduction pipeline.

## Scripts

| Script | Purpose |
|--------|---------|
| `00_download_public_data.R` | Download official author-provided RDS files from Zenodo for reference |
| `01_load_and_subset.R` | Load custom scRNA-seq data, collapse age groups, and subset neurogenic cells |
| `02_module_scores.R` | Calculate ferroptosis module scores |
| `03_figures_0_to_4.R` | Generate UMAP, violin, line, and trajectory plots |
| `04_de_and_volcano.R` | Run differential expression and volcano plots |
| `05_correlation_heatmaps.R` | Generate correlation heatmaps |

## Runner

`run_pipeline.R` runs the main scRNA-seq reproduction scripts in order.

## Boundaries

- Do not add spatial transcriptomics scripts here; use `R/spatial/`.
- Do not write generated outputs here.
- scRNA-seq outputs go to `data/processed/` and `figures/stage-2/`.
