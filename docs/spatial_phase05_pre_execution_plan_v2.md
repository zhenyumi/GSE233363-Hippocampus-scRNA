# Phase Spatial-05 Pre-Execution Plan v2 — Strict Author-Code Route

> **Version**: 2.0 (2026-06-03)
> **Route**: Strict author-code reproduction of Script 8 lines 71–993
> **Status**: BLOCKED — see Section 9

---

## 1. Scope: Strict Author-Code Only

Phase Spatial-05 follows **author Script 8 lines 71–993** as the strict reproduction target. The implementation priority is:

1. STutility `InputFromTable` — load per-sample Visium data from raw directories
2. `AddMetaData` — attach sample-level metadata
3. `LoadImages` — load tissue images per sample
4. `RegionNeighbours` — identify spatial neighbors per sample
5. Per-sample `Edge` / `nbs_Edge` / `RS1` / `RS2` labels — derive distance-based spatial hierarchy
6. Merge labels back to `Hippo_All` — combine per-sample results
7. Recode to `IS` / `NNS` / `ENS1` / `ENS2` — final inflammatory gradient categories
8. Gradient DESeq2 and figures — pseudobulk differential expression across gradient

**Input from Phase Spatial-04**: `metadata_hippo_func_region.csv` (FuncRegion + is_edge labels). Phase Spatial-04 completed with WARNING status (Middle 2.00% vs expected 2.50%; Old 6.29% vs expected 5.94%). This WARNING requires explicit user acknowledgement before the file can be consumed by Phase Spatial-05.

The current plan does **not** authorize any alternative method. If strict author-code prerequisites are missing, Phase Spatial-05 remains blocked.

### Excluded from Strict Route

The following are **explicitly excluded** and must not appear as equal next-action options:

- Custom distance-based neighbor detection
- Seurat coordinate-neighbor replacement
- Hippo_All image-slot shortcut (using existing RDS image slots instead of raw Visium files)
- IS vs non-IS shortcut DESeq2 (bypassing the full IS/NNS/ENS1/ENS2 hierarchy)
- Monocle 2 pseudotime implementation (paper-code discrepancy; not part of Script 8)
- Tangram
- GO/KEGG enrichment
- Biological reinterpretation
- WholeBrain / Chromium / DG / scRNA objects

---

## 2. Dependency Table (Corrected)

### Present in renv.lock

| Package | Version | Status | Usage in Script 8 |
|---------|---------|--------|-------------------|
| DESeq2 | 1.50.2 | ✓ present | Gradient pseudobulk DE (lines 759-835) |
| Matrix.utils | 0.9.8 | ✓ present | aggregate.Matrix for pseudobulk |
| SingleCellExperiment | 1.32.0 | ✓ present | DESeq2 input construction |
| ggplot2 | 3.5.2 | ✓ present | PCA, volcano plots |
| ggrepel | 0.9.6 | ✓ present | Gene labels on volcano |
| pheatmap | 1.0.13 | ✓ present | Heatmap of gradient DEGs |
| RColorBrewer | 1.1-3 | ✓ present | Color palettes |
| dplyr | present | ✓ verify at runtime | Data manipulation |
| stringr | present | ✓ verify at runtime | String operations |
| purrr | present | ✓ verify at runtime | Functional programming |
| Seurat | 5.5.0 | ✓ present | Object loading, metadata access |
| msigdbr | 26.1.0 | ✓ present | IFNγ module (used in Phase 04, inherited) |

### Missing from renv.lock

| Package | Status | Impact |
|---------|--------|--------|
| **STutility** | **MISSING — strict-route blocker** | Required for `InputFromTable`, `LoadImages`, `RegionNeighbours` (Script 8 lines 71-697) |
| Monocle 2 | Not used in Script 8 | Paper-method discrepancy only; not required for strict Script 8 route |

**Note**: STutility installation requires separate user approval. If approved later, record: install source, package version, renv.lock changes, and any system dependency issues.

---

## 3. Data-Blocked Status

### Conclusion

The strict author-code route is **data-blocked** because all raw Visium STutility inputs are missing.

### Current Data Inventory

`data/raw/GSE233363_official/` contains author RDS files only:
- `seurat_Visium_DG_All.rds`
- `seurat_Visium_Hippo_All.rds`
- `seurat_Visium_WholeBrain.rds`

### Missing Raw Directories

The 8 Old sample raw directories are entirely missing:

| Sample | Directory | Status |
|--------|-----------|--------|
| O_1/A1 | `data/raw/spatial/O_1/A1/` | MISSING |
| O_1/B1 | `data/raw/spatial/O_1/B1/` | MISSING |
| O_1/C1 | `data/raw/spatial/O_1/C1/` | MISSING |
| O_1/D1 | `data/raw/spatial/O_1/D1/` | MISSING |
| O_2/A1 | `data/raw/spatial/O_2/A1/` | MISSING |
| O_2/B1 | `data/raw/spatial/O_2/B1/` | MISSING |
| O_2/C1 | `data/raw/spatial/O_2/C1/` | MISSING |
| O_2/D1 | `data/raw/spatial/O_2/D1/` | MISSING |

### Required Files Per Sample

Each sample directory needs 4 files for STutility:

1. `filtered_feature_bc_matrix.h5` — gene expression matrix
2. `spatial/tissue_positions_list.csv` — spot coordinates
3. `spatial/tissue_hires_image.png` — tissue image
4. `spatial/scalefactors_json.json` — image scale factors

**Total**: 32 files (8 samples × 4 files). Currently 0 of 32 found.

### Consequence

Without these files, STutility `InputFromTable`, `LoadImages`, and `RegionNeighbours` cannot run. The Hippo_All image-slot shortcut is explicitly excluded from this strict plan.

---

## 4. Dependency-Blocked Status

STutility is missing from renv.lock. Strict-route dependency update requires separate user approval. **Do not install STutility in this plan.**

If STutility installation is later approved, record:
- Install source (CRAN, GitHub, Bioconductor)
- Package version
- renv.lock changes (diff)
- Any system dependency issues (e.g., imagemagick, fftw)

---

## 5. Author Method Extraction

### Script 8 Code Structure (lines 71–993)

| Lines | Function | Description |
|-------|----------|-------------|
| 71-110 | `InputFromTable` | Load per-sample Visium data from raw directories |
| 111-130 | `AddMetaData` | Attach sample-level metadata (age, sample ID) |
| 131-170 | `LoadImages` | Load tissue images for visualization |
| 171-250 | `RegionNeighbours` | Identify spatial neighbors per spot per sample |
| 251-350 | Edge/nbs_Edge/RS1/RS2 | Derive per-sample spatial hierarchy labels |
| 351-697 | Per-sample loop | Process all 8 Old samples, collect labels |
| 704-752 | Merge + Recode | Combine per-sample labels, recode IS/NNS/ENS1/ENS2 |
| 759-835 | DESeq2 | Pseudobulk differential expression across gradient |
| 845-993 | Figures | PCA, heatmap, volcano (IS vs NNS: 1326 down, 2357 up) |

### Label Hierarchy

```
Edge (TRUE) → IS (Inflammatory Spot)
nbs_Edge (TRUE) → NNS (Nearest Neighbor of IS)
RS1 (TRUE) → ENS1 (Edge Neighbor Step 1)
RS2 (TRUE) → ENS2 (Edge Neighbor Step 2)
All FALSE → Background (excluded from DESeq2)
```

### DESeq2 Parameters

- Design: `~ FuncRegion` (IS vs NNS, IS vs ENS1, IS vs ENS2)
- Pseudobulk: `aggregate.Matrix` by sample × FuncRegion
- FDR threshold: 0.05
- LFC threshold: 1.0

---

## 6. Data Inventory Per Sample

### Old Sample GEMgroup Mapping (Script 8)

| GEMgroup | Sample | Raw Directory | STutility InputTable Row |
|----------|--------|---------------|--------------------------|
| 9 | O_1_A1 | `O_1/A1/` | row 1 |
| 10 | O_1_B1 | `O_1/B1/` | row 2 |
| 11 | O_1_C1 | `O_1/C1/` | row 3 |
| 12 | O_1_D1 | `O_1/D1/` | row 4 |
| 13 | O_2_A1 | `O_2/A1/` | row 5 |
| 14 | O_2_B1 | `O_2/B1/` | row 6 |
| 15 | O_2_C1 | `O_2/C1/` | row 7 |
| 16 | O_2_D1 | `O_2/D1/` | row 8 |

**All 8 raw directories**: MISSING (0 of 32 files found).

---

## 7. STutility Call Inventory

### Required STutility Functions

| Function | Script 8 Lines | Purpose | Prerequisites |
|----------|---------------|---------|---------------|
| `InputFromTable` | 71-110 | Load per-sample Visium data | Raw H5 + spatial files per sample |
| `AddMetaData` | 111-130 | Attach sample metadata | Loaded STutility object |
| `LoadImages` | 131-170 | Load tissue images | Raw image files per sample |
| `RegionNeighbours` | 171-250 | Identify spatial neighbors | Loaded images + coordinates |

**All 4 functions require raw Visium files that are currently missing.**

---

## 8. Paper Figure Mapping

| Paper Panel | Script 8 Lines | Content | Phase 05 Output |
|-------------|---------------|---------|-----------------|
| Fig. 7f | 845-880 | Gene module spatial patterns | Spatial plots of gradient modules |
| Extended Data Fig. 9a | 880-910 | IS/NNS/ENS1/ENS2 spatial layout | Gradient category spatial plot |
| Extended Data Fig. 9b-d | 910-950 | Gradient DEG heatmaps | Heatmap of gradient DEGs |
| Extended Data Fig. 9i | 950-993 | CAS-Up scores in IS spots | Module score spatial plot |

---

## 9. Blocking Status

### STRICT AUTHOR-CODE ROUTE STATUS

| Blocker | Status | Resolution Required |
|---------|--------|---------------------|
| Phase 04 handoff | **USER-ACKNOWLEDGEMENT-BLOCKED** | User must acknowledge Phase 04 WARNING |
| STutility | **DEPENDENCY-BLOCKED** | Separate dependency update plan required |
| Raw Visium files | **DATA-BLOCKED** | 32 files for 8 Old samples missing |
| Method branch | Strict author-code selected | Alternatives excluded |

### Overall Status

**BLOCKED** until all three conditions are satisfied:
1. User acknowledges Phase 04 WARNING handoff
2. Raw Visium files are available for 8 Old samples
3. STutility dependency plan is approved and executed

---

## 10. Recommended Next Actions

### Action A: Confirm Phase 04 WARNING Handoff Acknowledgement

Phase Spatial-04 completed with WARNING (Middle 2.00% vs expected 2.50%; Old 6.29% vs expected 5.94%). The output `metadata_hippo_func_region.csv` is valid but has deviation from expected values. User must explicitly acknowledge this WARNING before Phase Spatial-05 can consume the file.

### Action B: Acquire/Check Raw Visium Files for 8 Old Samples

Locate or download the 32 required files:
- 8 samples × 4 files each
- `filtered_feature_bc_matrix.h5`, `tissue_positions_list.csv`, `tissue_hires_image.png`, `scalefactors_json.json`
- Check Zenodo, GEO supplementary files, or contact authors

### Action C: Approve STutility Dependency Update Plan

After raw data availability is confirmed (or in parallel), approve a separate dependency update plan for STutility. This plan must record:
- Install source
- Package version
- renv.lock changes
- Any system dependency issues

### Action D: Draft Phase 05 Implementation Plan

Once Actions A–C are satisfied, draft the Phase Spatial-05 implementation plan with:
- Script structure (`R/spatial/s05_inflammatory_gradient.R`)
- Per-sample STutility processing loop
- Label merge and recode logic
- DESeq2 gradient analysis
- Figure generation
- Validation checks

---

## 11. Future Outputs (If Prerequisites Satisfied)

**No outputs are produced during planning.** The following paths are proposed for future implementation:

### Data Outputs

```
data/processed/spatial/phase05_inflammatory_gradient/
├── phase05_provenance.csv
├── phase05_validation_summary.txt
├── per_sample_labels/
│   ├── O_1_A1_labels.csv
│   ├── O_1_B1_labels.csv
│   ├── O_1_C1_labels.csv
│   ├── O_1_D1_labels.csv
│   ├── O_2_A1_labels.csv
│   ├── O_2_B1_labels.csv
│   ├── O_2_C1_labels.csv
│   └── O_2_D1_labels.csv
├── merged_gradient_labels.csv
├── gradient_deseq2_IS_vs_NNS.csv
├── gradient_deseq2_IS_vs_ENS1.csv
└── gradient_deseq2_IS_vs_ENS2.csv
```

### Figure Outputs

```
figures/spatial/phase05/
├── phase05_gradient_categories_spatial.pdf
├── phase05_pca_gradient.pdf
├── phase05_heatmap_gradient_DEGs.pdf
├── phase05_volcano_IS_vs_NNS.pdf
├── phase05_volcano_IS_vs_ENS1.pdf
└── phase05_volcano_IS_vs_ENS2.pdf
```

---

## 12. Validation Criteria

Record PCA, heatmap, and volcano outputs and compare file existence / layout / reported DEG counts against author Script 8 comments and target paper panels. **Do not interpret biology.**

Specific checks:
- Confirm Phase Spatial-04 WARNING handoff has explicit user acknowledgement
- Confirm the Hippo RDS still has 16 images and 9921 spots before neighborhood approximation
- Confirm coordinate extraction is performed per image/sample, not across all images as one shared coordinate plane
- Compare RDS-derived IS/NNS/ENS-like counts against author Script 8 comments where available
- Compare resulting patterns against paper Fig. 7 and Extended Data Fig. 9 as qualitative/quantitative comparison
- Record all deviations from author Script 8

---

## Not Part of Strict Route

The following were considered but excluded from the strict author-code route:

- **Paper methods branch**: Paper describes distance-based hierarchy and Monocle 2 pseudotime, but Script 8 uses STutility. Strict route follows Script 8.
- **Custom coordinate neighbor approach**: Would bypass STutility's `RegionNeighbours` function.
- **Hippo_All image slots approach**: Would use existing RDS image slots instead of raw Visium files.
- **IS vs non-IS shortcut DESeq2**: Would bypass the full IS/NNS/ENS1/ENS2 hierarchy.

These can only be reconsidered if the strict route is permanently abandoned with user approval.
