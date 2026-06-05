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
- `docs/target_genes.csv` — 251 mitochondrial/ribosomal/translation genes (13 categories, GBK encoding). **This CSV is the canonical pipeline input.** `docs/target_genes.xlsx` is the original manual file (gitignored, not consumed by scripts).

Current scope:

- Stages Spatial-01 through Spatial-11/11b are implemented for DG/Hippo.
- Phase Spatial-05 is labeled as RDS-based reproduction/approximation because the raw Space Ranger/Visium files required by the author's STutility route are not available locally.
- Phase Spatial-06: Hippocampal regional atlas (CA1/CA2/CA3/ML/GCL/Hilus/DG).
- Phase Spatial-07: Reproduction validation and handoff (GO_WITH_CAVEATS).
- Phase Spatial-09: Mitochondrial target gene analysis plan.
- Phase Spatial-10: Target gene audit and region-aware expression summaries.
- Phase Spatial-11: CA1 vs CA3 pseudobulk DE, module scores, and coupling analysis.
- Phase Spatial-11b: Age-grouping stratified DE and cross-age pattern classification (PASS_WITH_WARNINGS).
- Phase Spatial-12: **Planned but NOT executed**. Young vs Aged regional comparison across CA1 / CA3 / DG.
- WholeBrain remains inventory-only unless a separate user-approved plan says otherwise.

### Reference-First Methodology Rules

This project uses two skill systems for planning:

- **ref-bio** (`.opencode/skills/ref-bio/`): Routes to authoritative upstream docs (DESeq2 vignettes, Seurat spatial vignettes, OSTA). Plans must open and review the actual documents — link-only YAML entries are not sufficient proof of review. Every consulted source must be recorded with `reviewed_status`, URL, specific section, and design impact.
- **bio-* skills** (`.opencode/skills/bio-*/`): Skill-specific guidance for bioinformatics tasks. There is NO single `bioskill/` directory. Plan phases select relevant bio-* skills by task type (e.g., `bio-differential-expression-deseq2-basics` for DE design). Python/Squidpy skills provide conceptual guidance only for this R/Seurat project.

See `AGENTS.md` for full rules and `docs/agent_reference_workflow.md` for the reference-first planning workflow.

Scripts:

| Script | Purpose |
|--------|---------|
| `R/spatial/s01_inspect_objects.R` | Inspect DG/Hippo spatial objects |
| `R/spatial/s02_dg_hippo_metadata.R` | Verify DG/Hippo metadata and structure |
| `R/spatial/s03_pseudobulk_deseq2_dg.R` | Pseudobulk DESeq2 DG old vs young |
| `R/spatial/s04_ifng_module_score.R` | IFNγ module score and FuncRegion handoff |
| `R/spatial/s05_rds_inflammatory_gradient.R` | RDS-based inflammatory gradient approximation |
| `R/spatial/s05b_plot_spatial_labels.R` | Plot spatial label maps |
| `R/spatial/s06_hippo_regional_atlas.R` | Hippocampal regional atlas (CA1/CA2/CA3/ML/GCL/Hilus/DG) |
| `R/spatial/s07_figure_validation.R` | Reproduction validation and handoff report |
| `R/spatial/s10_target_gene_audit_region_summary.R` | Target gene audit and region-aware expression summaries |
| `R/spatial/s11_ca1_ca3_de_module_coupling.R` | CA1 vs CA3 pseudobulk DE, module scores, coupling |
| `R/spatial/s11b_age_grouping_ca1_ca3_de.R` | Age-stratified CA1 vs CA3 DE and cross-age patterns |

See `R/spatial/README.md` for details.

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
