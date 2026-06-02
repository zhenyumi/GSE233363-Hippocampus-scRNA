# R/spatial/ — Spatial Transcriptomics Scripts

## Status

**Stage Spatial-01 implemented** — object inspection only. No downstream spatial analysis scripts have been implemented.

This directory is reserved for spatial transcriptomics analysis scripts that extend
the ferroDG reproduction project to Visium spatial data.

Current implementation scope is DG and Hippo only. `seurat_Visium_WholeBrain.rds`
is inventoried but should not be loaded or analyzed unless the user explicitly
approves a separate WholeBrain plan.

## Naming Convention

Scripts use the `s##_` prefix (not `0##_` which is reserved for scRNA-seq pipeline):

| Prefix | Purpose | Status |
|--------|---------|--------|
| `s01_` | Inspect spatial objects | Implemented: DG/Hippo inspection, WholeBrain gated |

Actual script names and purposes will be determined after inspecting the spatial objects.

## Script Pattern

Every spatial script should follow the established pattern:

```r
## =========================================================
## Script: s##_descriptive_name.R
## Purpose: [what this script does]
## Usage: Rscript R/spatial/s##_descriptive_name.R
## Prerequisites: [what must run first]
## =========================================================

# 1. Load
# 2. Inspect (record provenance/object structure)
# 3. Save lightweight summary
# 4. gc()
```

## Prerequisites Before Writing Scripts

Before creating any spatial script:

1. **Obtain spatial objects** — place RDS files in `data/raw/spatial/`
   or use the current canonical author-provided path `data/raw/GSE233363_official/`
2. **Consult ref-bio** — check `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml`
   for authoritative reference routing on the method you plan to implement
3. **Record object structure** — run `R/spatial/s01_inspect_objects.R` for DG and Hippo;
   outputs are written under `data/processed/spatial/inspection/`
4. **Choose methods based on actual data** — do not assume assay structure, spatial keys,
   or metadata columns

## Integration with Main Pipeline

Spatial scripts are NOT added to `run_pipeline.R` (the scRNA-seq pipeline runner).
A separate `run_spatial_pipeline.R` may be created when spatial analysis is ready.
