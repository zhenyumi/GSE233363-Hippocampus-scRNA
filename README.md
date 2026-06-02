# ferroDG Reproduction

Reproduction of ferroDG stage-2 analysis: ferroptosis in mouse dentate gyrus neurogenesis.

## Overview
- **Source**: https://github.com/ShilongZhang116/ferroDG
- **Data**: GSE233363 (scRNA-seq, Young vs Old mouse hippocampus)
- **Analysis**: Ferroptosis module scores, differential expression, correlation heatmaps

## Quick Start
```bash
# 1. Install R packages (first time only)
Rscript -e "install.packages('renv'); renv::restore()"

# 2. Place your Seurat object
# data/raw/GSE233363_custom/Seurat_combined_with_Celltype_Article.rds

# 3. Run pipeline
Rscript run_pipeline.R
```

## Output
- **Figures**: `figures/stage-2/` (Fig 0-7)
- **Data**: `data/processed/` (RDS + CSV files)
- **Gene lists**: `gene_lists/`

## Scripts
| Script | Purpose |
|--------|---------|
| `01_load_and_subset.R` | Load custom data, subset neurogenic cells |
| `02_module_scores.R` | Ferroptosis module scores |
| `03_figures_0_to_4.R` | UMAP, violin, line, trajectory plots |
| `04_de_and_volcano.R` | Differential expression + volcano |
| `05_correlation_heatmaps.R` | Correlation heatmaps |

## Data Organization
- `data/raw/GSE233363_custom/` — Independently aligned from FASTQ (used by pipeline)
- `data/raw/GSE233363_official/` — Author-provided Zenodo files (for reference, download via `R/00_download_public_data.R`)

## Environment
- R 4.5.1, Seurat 5.5.0
- renv for package management
- macOS (tested on Apple Silicon)

## References
- [ferroDG GitHub](https://github.com/ShilongZhang116/ferroDG)
- [GSE233363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363)
