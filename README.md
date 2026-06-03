# ferroDG Reproduction

Reproduction and extension of the ferroDG / GSE233363 analysis in one renv-managed R project.

The repository currently has two related analysis branches:

- **scRNA-seq ferroDG stage-2 reproduction**: ferroptosis analysis in mouse dentate gyrus neurogenesis.
- **Spatial transcriptomics reproduction**: hippocampus/DG Visium reproduction from author-provided RDS objects, with later CA1/CA3/DG mitochondrial analysis planned only after the hippocampus workflow is validated.

## Quick Start

```bash
# Restore the project R environment
Rscript -e "install.packages('renv'); renv::restore()"

# Run the scRNA-seq reproduction pipeline
Rscript run_pipeline.R
```

The spatial branch is run stage-by-stage from `R/spatial/`; it is not part of `run_pipeline.R`. Optional Python work is reserved for a future Tangram/scanpy branch and lives under `python/spatial/` only if approved.

## Repository Structure

The canonical structure is documented in `docs/repository_structure.md`.

| Path | Purpose |
|------|---------|
| `R/` | R source root; contains branch-specific subdirectories only |
| `R/scrna/` | scRNA-seq R scripts and pipeline code |
| `R/spatial/` | spatial transcriptomics scripts, using `s##_` names |
| `python/spatial/` | optional future Python spatial workflows such as Tangram/scanpy |
| `docs/` | plans, reference registries, result summaries, and copied author-code reference |
| `gene_lists/` | small reviewed source gene lists |
| `data/` | raw and processed local data; gitignored |
| `figures/` | generated figures; gitignored |
| `results/`, `reports/`, `cache/` | generated outputs, reports, and scratch files; gitignored |

`analysis/` is retired in this repository. Historical upstream paths such as `analysis/stage-2/` are documented only as source references; maintained code lives under `R/scrna/`, `R/spatial/`, and optionally `python/spatial/`.

## scRNA-seq Branch

Input:

- `data/raw/GSE233363_custom/Seurat_combined_with_Celltype_Article.rds`

Scripts:

| Script | Purpose |
|--------|---------|
| `R/scrna/00_download_public_data.R` | Download official author-provided RDS files for reference |
| `R/scrna/01_load_and_subset.R` | Load custom data and subset neurogenic cells |
| `R/scrna/02_module_scores.R` | Calculate ferroptosis module scores |
| `R/scrna/03_figures_0_to_4.R` | UMAP, violin, line, and trajectory plots |
| `R/scrna/04_de_and_volcano.R` | Differential expression and volcano plots |
| `R/scrna/05_correlation_heatmaps.R` | Spearman correlation heatmaps |

Outputs:

- `data/processed/`
- `figures/stage-2/`

## Spatial Branch

Current input:

- `data/raw/GSE233363_official/seurat_Visium_DG_All.rds`
- `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds`

Current scope:

- Stages Spatial-01 through Spatial-05 are implemented for DG/Hippo.
- Phase Spatial-05 is labeled as RDS-based reproduction/approximation because the raw Space Ranger/Visium files required by the author's STutility route are not available locally.
- WholeBrain remains inventory-only unless a separate user-approved plan says otherwise.

Scripts:

- See `R/spatial/README.md`.

Optional Python:

- `python/spatial/` is reserved for Phase Spatial-08 Tangram/scanpy if reproduction validation requires it.
- No Python packages or environments are required for the current RDS-based DG/Hippo stages.

Outputs:

- `data/processed/spatial/`
- `figures/spatial/`

## Data and Git Hygiene

Do not commit generated data, figures, caches, large RDS/H5/MTX files, or root-level plot artifacts such as `Rplots.pdf`.

Keep source and generated artifacts separate:

- Source code: `R/scrna/`, `R/spatial/`, and optionally `python/spatial/`
- Committed documentation: `docs/`
- Generated outputs: `data/`, `figures/`, `results/`, `reports/`, `cache/`

## Environment

- R 4.5.1
- Seurat 5.5.0
- renv-managed project library
- macOS / Apple Silicon development environment

Only update `renv.lock` after intentional dependency changes.

## References

- [ferroDG GitHub](https://github.com/ShilongZhang116/ferroDG)
- [GSE233363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363)
- `docs/reference_sources.md`
- `docs/spatial_reproduction_plan.md`
