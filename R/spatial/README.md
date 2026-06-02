# R/spatial/ — Spatial Transcriptomics Scripts

## Status

**SCAFFOLD ONLY** — no analysis scripts have been implemented.

This directory is reserved for spatial transcriptomics analysis scripts that extend
the ferroDG reproduction project to Visium spatial data.

## Naming Convention

Scripts use the `s##_` prefix (not `0##_` which is reserved for scRNA-seq pipeline):

| Prefix | Purpose | Status |
|--------|---------|--------|
| `s01_` | Inspect spatial objects | Not created |

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
2. **Consult ref-bio** — check `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml`
   for authoritative reference routing on the method you plan to implement
3. **Record object structure** — run inspection and document in `docs/data_manifest_spatial.md`
4. **Choose methods based on actual data** — do not assume assay structure, spatial keys,
   or metadata columns

## Integration with Main Pipeline

Spatial scripts are NOT added to `run_pipeline.R` (the scRNA-seq pipeline runner).
A separate `run_spatial_pipeline.R` may be created when spatial analysis is ready.
