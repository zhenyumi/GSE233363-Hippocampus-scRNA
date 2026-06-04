# Phase Spatial-06: Hippocampus Regional Atlas Reproduction

## Status
- **Plan version**: v2 (approved 2026-05-29)
- **Script**: `R/spatial/s06_hippo_regional_atlas.R`
- **Execution mode**: RDS-based reproduction/approximation

## Objective
Validate the hippocampal regional framework (CA1, CA2, CA3, ML, GCL, Hilus, derived DG) using the inspected `seurat_Visium_Hippo_All.rds` object. Produce per-region summary tables, spatial maps, and reduction plots to confirm data integrity before downstream regional analyses.

## Input Files
| File | Purpose |
|------|---------|
| `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` | Primary Hippo Seurat object |
| `data/processed/spatial/phase02_metadata/Hippo/hippo_region_distribution.csv` | Phase 02 Region counts (validation reference) |
| `data/processed/spatial/phase02_metadata/Hippo/hippo_age_region_crosstab.csv` | Phase 02 Age×Region crosstab |
| `data/processed/spatial/phase02_metadata/Hippo/hippo_metadata_summary.csv` | Phase 02 metadata summary |
| `data/processed/spatial/phase02_metadata/cross_object_comparison.csv` | Cross-object DG validation |
| `data/processed/spatial/phase02_metadata/phase02_validation_summary.txt` | Phase 02 validation summary |
| `data/processed/spatial/inspection/Hippo/spatial_info.csv` | Spatial image info |
| `docs/original_code_from_paper/Scripts/R/6. Spatial-Data process.R` | Author reference (read-only) |

## Output Dirs
- `data/processed/spatial/phase06_hippo_regional_atlas/`
- `figures/spatial/phase06_hippo_regional_atlas/`

## Required Outputs

### CSV/TXT (9 files)
1. `hippo_region_counts.csv` — Region spot counts and proportions
2. `hippo_region_age_gemgroup_counts.csv` — Region×Age×GEMgroup counts
3. `hippo_region_image_counts.csv` — Region×Image counts
4. `hippo_region_age_crosstab.csv` — Region×Age crosstab
5. `hippo_region_coordinate_summary.csv` — Per-region coordinate ranges
6. `hippo_derived_dg_validation.csv` — Derived DG (ML+GCL+Hilus) vs cross-object DG
7. `hippo_region_metadata.csv` — Full per-spot metadata with region_group
8. `phase06_provenance.csv` — Script provenance (version, paths, timestamps)
9. `phase06_validation_summary.txt` — Human-readable validation report

### PNG (26 files)
- 16× `phase06_spatial_region_<image>.png` — Per-slice spatial region maps
- 1× `phase06_spatial_region_all_slices.png` — Overview 4×4 grid
- 9× reduction plots: `phase06_{umap,tsne,pca}_by_{region,age,gemgroup}.png`

## Validation Matrix (V01–V25)

### Object Invariants (STOP on fail)
| ID | Check | Threshold |
|----|-------|-----------|
| V01 | Hippo spot count | == 9921 |
| V02 | Feature count (SCT assay) | == 17096 |
| V03 | Feature count (Spatial assay) | == 32285 |
| V04 | Image count | == 16 |
| V05 | Region column exists | TRUE |
| V06 | Region levels | 6 (CA1, CA2, CA3, ML, GCL, Hilus) |
| V07 | SCT assay present | TRUE |
| V08 | Spatial assay present | TRUE |
| V09 | Age column exists | TRUE |

### Region Counts vs Phase 02 (STOP if >5% off)
| ID | Region | Expected |
|----|--------|----------|
| V10 | CA1 | 4671 |
| V11 | CA2 | 565 |
| V12 | CA3 | 2321 |
| V13 | ML | 1480 |
| V14 | GCL | 645 |
| V15 | Hilus | 239 |
| V16 | Total | 9921 |
| V17 | Derived DG (ML+GCL+Hilus) | 2364 |

### Runtime Checks
| ID | Check | On Fail |
|----|-------|---------|
| V18 | Reductions (pca, umap, tsne) available | WARNING |
| V19 | GetTissueCoordinates succeeds ≥1 image | STOP if all fail |
| V20–V25 | Required output files exist | FAIL (cannot PASS) |

## Stop Conditions (SC1–SC10)
| ID | Condition | Action |
|----|-----------|--------|
| SC1 | RDS load failure | STOP |
| SC2 | Dimension mismatch (V01/V02/V03) | STOP |
| SC3 | Assay missing (V07/V08) | STOP |
| SC4 | Region unexpected levels | STOP |
| SC5 | Region counts off >5% (V10–V17) | STOP |
| SC6 | All image coordinates fail | STOP; single fail → record FAIL, continue |
| SC7 | Memory allocation error | STOP |
| SC8 | Package load failure | STOP |
| SC9 | Rplots.pdf created mid-script | WARNING + log |
| SC10 | Output file missing | FAIL (blocks overall PASS) |

## Color Scheme
| Region | Hex |
|--------|-----|
| CA1 | #1f77b4 |
| CA2 | #ff7f0e |
| CA3 | #2ca02c |
| ML | #d62728 |
| GCL | #9467bd |
| Hilus | #8c564b |

## Memory Strategy
- Peak: ~2.5 GB (single object load)
- gc() before and after readRDS
- Record object.size() after load
- Single object at all times; no WholeBrain, no DG loading
- Call gc() after each major plot section

## Constraints
- No PDF generation (default disabled)
- No package installation
- No renv.lock modification
- No writes to Seurat object metadata
- region_group only in independent metadata_df
- Use existing reductions only (pca, umap, tsne)
