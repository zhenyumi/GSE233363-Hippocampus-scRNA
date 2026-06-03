# Data Manifest — Spatial Transcriptomics

**Status**: Author-provided spatial RDS objects have been downloaded to the official
local data directory. DG and Hippo have been inspected in Stages Spatial-01/02.
Raw Space Ranger/Visium files required by the author STutility code are not
currently available locally.

## Author Spatial RDS Files

| # | Expected File | Source | Status | Path | Size |
|---|---------------|--------|--------|------|------|
| 1 | seurat_Visium_DG_All.rds | JessbergerLab / GSE233363 / Zenodo | Available; inspected | `data/raw/GSE233363_official/seurat_Visium_DG_All.rds` | ~445 MB |
| 2 | seurat_Visium_Hippo_All.rds | JessbergerLab / GSE233363 / Zenodo | Available; inspected | `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` | ~1.3 GB |
| 3 | seurat_Visium_WholeBrain.rds | JessbergerLab / GSE233363 / Zenodo | Available; inventory only | `data/raw/GSE233363_official/seurat_Visium_WholeBrain.rds` | ~4.9 GB |

Current analysis focus is DG and Hippo. WholeBrain must not be loaded without a
separate user-approved plan because of memory risk and scope.

## Raw Visium / Space Ranger Files

The author Script 8 STutility route references per-sample files such as:

- `filtered_feature_bc_matrix.h5`
- `spatial/tissue_positions_list.csv`
- `spatial/tissue_hires_image.png`
- `spatial/scalefactors_json.json`

These files are **not currently available locally**. The current Phase Spatial-05
route is therefore RDS-based reproduction/approximation using the image slots and
coordinates already present in `seurat_Visium_Hippo_All.rds`. If raw Visium files
are later found, record them here and reopen a separate strict raw-file/STutility
plan.

## Metadata Checklist for First Load

For new objects or any future reinspection, record the following:

### Object Structure
- [ ] File path on disk
- [ ] File size (bytes/MB/GB)
- [ ] Seurat version used to create object
- [ ] `DefaultAssay()` value
- [ ] All assay names (`Assays(obj)`)
- [ ] All reduction names (`Reductions(obj)`)

### Spatial Information
- [ ] Spatial key (`SpatialKey(obj)`)
- [ ] Image names (`Images(obj)`)
- [ ] Coordinate ranges (x_min, x_max, y_min, y_max)
- [ ] Spot count (`ncol(obj)`)
- [ ] Gene count (`nrow(obj)`)

### Metadata Columns
- [ ] Sample IDs (column name, unique values)
- [ ] Age/timepoint labels (column name, unique values)
- [ ] Anatomical region labels (column name, unique values, if present)
- [ ] Any clustering or cell type annotations (column names, unique values)

### Provenance
- [ ] GEO accession (should be GSE233363)
- [ ] Alignment method (Space Ranger or custom)
- [ ] Reference genome used
- [ ] Author attribution

## Expected Directory Layout

```
data/raw/GSE233363_official/
├── seurat_Visium_DG_All.rds
├── seurat_Visium_Hippo_All.rds
├── seurat_Visium_WholeBrain.rds
└── METADATA.json

data/raw/spatial/
└── [reserved for future raw Space Ranger/Visium files if found]
```

## Notes

- Spatial object file sizes are large; measure file size and memory behavior at first load
- Do NOT assume the objects use SCT assay, specific spatial key names, or specific metadata column names
- The JessbergerLab GitHub repo (https://github.com/JessbergerLab/AgingNeurogenesis_Transcriptomics)
  may provide additional context on data processing
- Age/sample metadata has been recorded for DG and Hippo; still re-check critical columns at runtime
- RDS-contained image slots and coordinates can support RDS-based approximation, but do not claim strict STutility/raw-file reproduction without raw Visium files
