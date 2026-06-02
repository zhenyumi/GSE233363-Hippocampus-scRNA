# AGENTS.md - ferroDG Reproduction Project

## Project Overview
Reproduce ferroDG (https://github.com/ShilongZhang116/ferroDG) stage-2 analysis:
ferroptosis analysis in mouse dentate gyrus neurogenesis using scRNA-seq (GSE233363).

## Key Context
- **Species**: Mouse (gene symbols: First letter uppercase, rest lowercase)
- **Seurat version**: v5.5.0 (handles `avg_log2FC` vs `avg_logFC`)
- **Cell types**: qNSC, nIPC, Neuroblast, GC (neurogenic lineage)
- **Comparison**: Young vs Old (Middle timepoint excluded)
- **Environment**: renv-managed R 4.5.1, macOS

## Scripts (run in order)
1. `R/01_load_and_subset.R` - Load custom data, collapse age, subset neurogenic cells
2. `R/02_module_scores.R` - Calculate ferroptosis module scores
3. `R/03_figures_0_to_4.R` - UMAP, violin, line, trajectory plots
4. `R/04_de_and_volcano.R` - Differential expression + volcano plots
5. `R/05_correlation_heatmaps.R` - Spearman correlation heatmaps
6. `R/00_download_public_data.R` - Download official author-provided RDS from Zenodo

## When Modifying Figures
- Color schemes defined in each script (look for `scale_color_manual` or `color =`)
- Use `cairo_pdf()` with `tryCatch` fallback to `pdf()` (X11 not installed)
- Save both PNG (always works) and PDF (best effort)
- Call `gc()` after each large plot to manage memory

## When Adding New Analysis
- Create new script: `R/06_<descriptive_name>.R`
- Follow pattern: load data -> compute -> save -> plot -> gc()
- Add to `run_pipeline.R` scripts vector
- Update this file with new script description

## Gene Symbol Convention
- Original paper uses Human symbols (uppercase)
- Convert to Mouse: `tolower()` + capitalize first letter
- Example: "TFRC" -> "Tfrc", "SLC7A11" -> "Slc7a11"

## Known Missing Genes (expected)
- Inhibitor: Mt1g, Akr1c1-3, G6pd (not in mouse genome under these symbols)
- Regulator: Tp53, Ras (not in mouse genome under these symbols)
- Scripts handle these with guards (skip if missing)

## Data Paths
- Custom data: `data/raw/GSE233363_custom/` (gitignored, independently aligned from FASTQ — USED by pipeline)
- Official data: `data/raw/GSE233363_official/` (gitignored, author-provided Zenodo files — for reference)
- Processed: `data/processed/` (gitignored, generated)
- Figures: `figures/stage-2/` (gitignored, generated)
- Gene lists: `gene_lists/` (committed)

## File Structure and Temporary Outputs
- Keep source code and generated artifacts separate: scripts belong in `R/`, committed docs in `docs/`, and generated outputs in `data/processed/`, `figures/`, `results/`, `reports/`, or `cache/`.
- Do not write temporary files, downloaded data, logs, plots, or intermediate tables into the repository root, `R/`, `docs/`, or `gene_lists/` unless they are intentionally reviewed source/docs.
- For scRNA-seq outputs, use existing conventions: processed objects/tables in `data/processed/` and figures in `figures/stage-2/`.
- For spatial outputs, use spatial-specific subdirectories: `data/processed/spatial/`, `figures/spatial/`, `results/spatial/`, `reports/spatial/`, and `cache/spatial/`.
- Short-lived scratch files should go under `cache/` or `cache/spatial/` and be removed when no longer needed.
- Lightweight audit logs and validation reports may go under `reports/` or `reports/spatial/`; committed planning/reference docs stay under `docs/`.
- Large or generated files must remain gitignored. Do not commit RDS/RData, H5/H5AD, MTX, image tiles, compressed raw data, cache files, or generated result directories.
- Before adding a new output path, check existing project conventions and `.gitignore`; prefer adding a narrow subdirectory over scattering files across the tree.

## Memory and Large-Object Management
- Prefer staged scripts over monolithic analyses; each script should load only the inputs needed for that stage, write a checkpoint or lightweight summary, then release memory.
- Avoid loading multiple large Seurat or Visium objects at the same time. For spatial work, do not load WholeBrain, Hippo, and DG objects together unless there is a documented need.
- After saving a stage output, remove large temporary objects with `rm()` and call `gc()` before the next major operation.
- Avoid unnecessary object copying, especially Seurat object copies created by repeated assignment, broad subsetting, or adding temporary metadata columns.
- Preserve sparse matrices where possible; do not coerce large count matrices to dense matrices unless explicitly justified and memory-safe.
- Save lightweight CSV/TXT summaries for inspection and validation when full object copies are not required.
- Do not save duplicate large RDS files without a clear reason, provenance note, and gitignored destination.
- On macOS, assume available RAM may be the limiting factor for large Seurat/Visium objects; record memory failures and retry with smaller staged inspections rather than forcing the full workflow.

## Running
```bash
# Individual scripts
Rscript R/01_load_and_subset.R
Rscript R/02_module_scores.R
# ...

# Full pipeline
Rscript run_pipeline.R
```

## Debugging Tips
- If `FetchData` fails: check gene names exist in `rownames(seurat_obj)`
- If `pheatmap` errors on NA: add `cor_mat[is.na(cor_mat)] <- 0`
- If PDF fails: X11/cairo not installed, PNGs still work
- Memory: call `gc()` after large operations

## Spatial Transcriptomics Extension (Phase Spatial-00)
**Status**: SCAFFOLD ONLY — no analysis implemented yet.

### Key Distinctions
- **scRNA-seq cells ≠ Visium spots**: Visium spots may capture signal from multiple cells and are not equivalent to single cells.
- **Spots are not independent biological replicates**: spatial neighbors share local microenvironment
- **Metadata must be inspected from actual objects**: do not assume assay names, spatial keys, image slots, or anatomical labels
- **Consult .opencode/skills/ref-bio/ for authoritative reference routing** when planning methods

### Data Paths (expected, unverified)
- Spatial input: `data/raw/spatial/` (gitignored, user must supply Visium RDS files)
- Processed spatial: `data/processed/spatial/` (gitignored, generated)
- Spatial figures: `figures/spatial/` (gitignored, generated)
- Spatial scripts: `R/spatial/` (tracked, `s##_` prefix)

### Script Naming
- Spatial scripts use `s##_` prefix (e.g., `s01_inspect_objects.R`), not `0##_`
- Script pattern: load → inspect → record provenance/object structure → save lightweight summary → gc()
- Before writing any spatial script: load one object, record its structure, consult ref-bio for reference routing

### Memory Constraints (macOS 16GB)
- Load one spatial object at a time
- Use lightweight summaries, not duplicate RDS copies
- Call `gc()` after large operations
- Memory management specifics deferred to official Seurat docs and actual object structure

### Reference Routing
- Project-specific refs (paper, GitHub, GEO, Zenodo): see `docs/reference_sources.md`
- Authoritative reference catalog: `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml`
- Seurat spatial vignettes: https://satijalab.org/seurat/
- OSTA (Orchestrating Spatial Transcriptomics Analysis): https://bioconductor.org/books/release/OSTA/

### Future Object Inspection Protocol
When spatial RDS files become available, the first script (`s01_inspect_objects.R`) must:
1. Record file path and size
2. Record Seurat version used to create the object
3. Record `DefaultAssay()` and all assay names
4. Record spatial key (`SpatialKey()`), image names, coordinate ranges
5. Record all metadata columns and their types/levels
6. Record spot count and gene count
7. Record provenance (GEO accession, alignment method, author)
