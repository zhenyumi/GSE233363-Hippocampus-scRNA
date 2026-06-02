# Data Manifest — Spatial Transcriptomics (UNVERIFIED)

**Status**: This manifest lists EXPECTED files based on the GSE233363 dataset.
All entries are UNVERIFIED until actual files are inspected.

## Expected Author Spatial RDS Files

| # | Expected File | Source | Status | Path | Size |
|---|---------------|--------|--------|------|------|
| 1 | seurat_Visium_DG_All.rds | JessbergerLab / GSE233363 | UNVERIFIED | TBD | TBD |
| 2 | seurat_Visium_Hippo_All.rds | JessbergerLab / GSE233363 | UNVERIFIED | TBD | TBD |
| 3 | seurat_Visium_WholeBrain.rds | JessbergerLab / GSE233363 | UNVERIFIED | TBD | TBD |

Age/sample metadata TBD after first object inspection.

## Metadata Checklist for First Load

When RDS files are available, record the following for each object:

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
data/raw/spatial/
├── [seurat_Visium_DG_All.rds]        # UNVERIFIED - name TBD
├── [seurat_Visium_Hippo_All.rds]     # UNVERIFIED - name TBD
└── [seurat_Visium_WholeBrain.rds]    # UNVERIFIED - name TBD
```

## Notes

- Spatial object file sizes are unverified and may be large; measure file size and memory behavior at first load
- Do NOT assume the objects use SCT assay, specific spatial key names, or specific metadata column names
- The JessbergerLab GitHub repo (https://github.com/JessbergerLab/AgingNeurogenesis_Transcriptomics)
  may provide additional context on data processing
- Age/sample metadata must be recorded after first object inspection; do not assume labels
