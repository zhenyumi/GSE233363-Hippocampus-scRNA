# R/spatial/ — Spatial Transcriptomics Scripts

## Status

Stages Spatial-01 through Spatial-04 have been implemented for DG/Hippo spatial
reproduction. Phase Spatial-05 is planned as an RDS-based reproduction/approximation
because the raw Space Ranger/Visium files required by the author STutility code are
not currently available.

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
| `s02_` | Verify DG/Hippo metadata and structure | Implemented |
| `s03_` | DG pseudobulk DESeq2 via author Script 7 route | Implemented |
| `s04_` | IFN-gamma module score and FuncRegion handoff | Implemented with WARNING handoff |
| `s05_` | RDS-based inflammatory gradient approximation | Plan next; not implemented |

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

## Phase Spatial-05 Boundary

The author Script 8 STutility section uses external raw Visium files
(`filtered_feature_bc_matrix.h5`, `tissue_positions_list.csv`,
`tissue_hires_image.png`, and `scalefactors_json.json`) that are not currently
available locally. The available Hippo/DG RDS objects do contain Seurat `VisiumV1`
image slots and coordinates, so the next phase should be planned as
RDS-based reproduction/approximation.

Requirements for the Phase Spatial-05 plan:

1. Label outputs as RDS-based approximation, not strict author-code reproduction.
2. Use `seurat_Visium_Hippo_All.rds` only; do not load DG, WholeBrain, or Chromium unless a later approved plan requires it.
3. Use `metadata_hippo_func_region.csv` from Phase Spatial-04 only after explicit user acknowledgement of its WARNING status.
4. Derive neighborhood labels per image/sample; do not combine all image coordinate systems into one plane.
5. Save only lightweight metadata, tables, and figures; do not save modified full Seurat objects.
6. Do not install or require STutility for this RDS-based route unless raw Visium files are later found and the user explicitly approves a strict raw-file branch.

## Integration with Main Pipeline

Spatial scripts are NOT added to `run_pipeline.R` (the scRNA-seq pipeline runner).
A separate `run_spatial_pipeline.R` may be created when spatial analysis is ready.
