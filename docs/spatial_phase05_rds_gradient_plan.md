# Phase Spatial-05 Pre-Execution Plan v3 — RDS-Based Coordinate Neighbor Graph Approximation

**Route**: RDS-based coordinate neighbor graph approximation  
**Status**: SELECTED / PROPOSED — pending user approval and Phase 04 WARNING acknowledgement  
**Target**: Author Script 8 lines 71–993 (STutility → RegionNeighbours → IS/NNS/ENS1/ENS2 → DESeq2 + figures)  

---

## Scope

Phase 05 approximates the IS/NNS/ENS1/ENS2 neighbor classification from author Script 8 lines 71–993 using **only** `seurat_Visium_Hippo_All.rds` and the Phase 04 handoff.

### Included
- Load and validate `seurat_Visium_Hippo_All.rds`
- Runtime verification of 16 image ↔ GEMgroup 1:1 mapping (never hardcoded)
- Per-image coordinate neighbor graph construction
- Ring expansion: Edge (Phase 04 seed) → nbs_Edge → RS1 → RS2
- Author-style label recoding: Edge → IS, nbs_Edge → NNS, RS1 → ENS1, RS2 → ENS2
- Pseudobulk DESeq2 on Old subset (IS vs NNS)
- PCA, heatmap, volcano figures
- All outputs labelled "RDS-approximated"

### Excluded
- STutility (NOT in renv.lock, NOT installed, NOT used)
- Monocle 2 (NOT used in Script 8; paper-code discrepancy only)
- Tangram, GO/KEGG enrichment, biological reinterpretation
- WholeBrain / Chromium / DG / scRNA objects
- Custom distance-based neighbor detection as a separate method branch
- Seurat coordinate-neighbor replacement as a separate method branch
- Hippo_All image-slot shortcut as a separate method branch
- IS vs non-IS shortcut DESeq2 as a separate method branch
- Raw Visium files (0 of 32 available; RDS-based route bypasses)

---

## Input Files

| File | Status | Notes |
|---|---|---|
| `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` | Verified exists | Phase 02 confirmed 9921 spots, 17096 genes (SCT), 16 VisiumV1 images |
| `data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt` | Verified exists | Phase 04 handoff — WARNING on Middle and Old |
| `data/processed/spatial/phase04_ifng_module/metadata_hippo_func_region.csv` | Verified exists | 9921 rows, 9 columns: barcode, id_cleaned, id_unique, Age, Region, GEMgroup, Hallmark_IFNg, FuncRegion, is_edge |
| `data/processed/spatial/phase04_ifng_module/phase04_provenance.csv` | Verified exists | Phase 04 provenance |
| `data/processed/spatial/phase02_metadata/Hippo/hippo_object_structure_summary.csv` | Verified exists | Object structure from Phase 02 |
| `data/processed/spatial/phase02_metadata/Hippo/hippo_image_spot_counts.csv` | Verified exists | 16 images, per-image spot counts |
| `data/processed/spatial/phase02_metadata/Hippo/hippo_metadata_summary.csv` | Verified exists | Metadata column summary |
| `data/processed/spatial/inspection/Hippo/spatial_info.csv` | Verified exists | Image names, coordinate keys, VisiumV1 class |

**No cache paths used.** All inputs from `data/processed/spatial/`.

---

## Key RDS Data Points

- **Object**: `seurat_Visium_Hippo_All.rds` (Seurat V5, 9921 spots × 17096 genes in SCT assay)
- **DefaultAssay**: `SCT`
- **Spatial assay**: 32285 features (genes), present in `Assays(hippo)`
- **Spatial key**: `slice1_` (all 16 images share this prefix)
- **16 images**: `slice1` through `slice1_16` (all VisiumV1)
- **Per-image coordinate ranges are unique** — never merge across images
- **Coordinate columns**: `imagecol` (pixel column), `imagerow` (pixel row)
- **GEMgroup → Age**: 1–4 = Young, 5–8 = Middle, 9–16 = Old
- **Phase 04 handoff**: `FuncRegion` column (IS/NS), `is_edge` boolean (426 TRUE / 9495 FALSE)

### Image Slot Clarification
- RDS contains Seurat VisiumV1 image slots and per-image coordinate ranges.
- What is **available**: image slots, `imagecol`/`imagerow` coordinates, spot-to-image membership via `Cells()`.
- What is **unavailable**: raw Space Ranger files (h5, tissue_positions_list.csv, hires PNG, scalefactors), STutility `LoadImages` + `RegionNeighbours` image-mask-based neighbor calling.
- Therefore the RDS route uses RDS-contained coordinates but **cannot claim equivalence** to STutility image-mask-based neighbor calling. All outputs labelled "RDS-approximated".

---

## Dependency Table

| Package | renv.lock | Role | Notes |
|---|---|---|---|
| Seurat (v5.5.0) | Present | Object loading, coordinate extraction | DefaultAssay = SCT |
| DESeq2 | v1.50.2 | Pseudobulk DE on Old subset | |
| Matrix.utils | v0.9.8 | aggregate.Matrix for pseudobulk | |
| SingleCellExperiment | v1.32.0 | SCE container for DESeq2 | |
| SummarizedExperiment | (SCE dependency) | assay(), colData() accessors | |
| ggplot2 | Present | Figures | |
| ggrepel | Present | Volcano labels | |
| pheatmap | Present | Heatmap | |
| RColorBrewer | Present | Heatmap palette | |
| dplyr | Present | Data manipulation | Verify load at runtime |
| stringr | Present | String operations | Verify load at runtime |
| purrr | Present | Functional programming | Verify load at runtime |
| **STutility** | **NOT present** | **NOT installed, NOT used** | |
| **Monocle 2** | **NOT present** | **NOT used** | |

**No new packages. No renv.lock changes.** STOP if any required package fails to load.

---

## Output Files

### Data: `data/processed/spatial/phase05_rds_gradient/`

| File | Description |
|---|---|
| `phase05_provenance.csv` | Provenance record |
| `image_gemgroup_mapping.csv` | Runtime-verified 16 image ↔ GEMgroup 1:1 mapping |
| `phase05_per_image_label_counts.csv` | Label counts per image (Edge/nbs_Edge/RS1/RS2) |
| `phase05_final_label_counts.csv` | Final IS/NNS/ENS1/ENS2 counts by Age |
| `phase05_old_pseudobulk_counts.rds` | Pseudobulk count matrix for Old (RDS, sparse) |
| `phase05_deseq2_old_results.csv` | DESeq2 results (IS vs NNS, Old) |
| `phase05_deseq2_old_all_genes.csv` | Full gene table with padj and log2FC |
| `phase05_neighbor_param_log.csv` | Per-image graph parameters (mode NN distance, margin) |
| `phase05_validation_summary.txt` | Validation results |

### Figures: `figures/spatial/phase05/`

| File | Description |
|---|---|
| `phase05_pca_old_by_region.png` | PCA colored by FuncRegion |
| `phase05_heatmap_top_degs.png` | Top DEG heatmap |
| `phase05_volcano_is_vs_nns.png` | Volcano IS vs NNS |

Both PNG (always) + PDF (best effort via tryCatch).

**No .RData. No save.image().**

---

## Algorithm

### Step 1: Load and Validate Hippo Object

```r
hippo <- readRDS("data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds")

stopifnot(DefaultAssay(hippo) == "SCT")
stopifnot("Spatial" %in% Assays(hippo))
stopifnot(ncol(hippo) == 9921)
stopifnot(length(Images(hippo)) == 16)

counts_spatial <- Seurat::GetAssayData(
  hippo,
  assay = "Spatial",
  layer = "counts"
)
stopifnot(all(dim(counts_spatial) == c(32285, 9921)))
```

### Step 2: Runtime Image ↔ GEMgroup Mapping

```r
spots_per_image <- lapply(Images(hippo), function(img) {
  Cells(hippo)[hippo@images[[img]]$boundaries$centroids@cells]
})
names(spots_per_image) <- names(hippo@images)

gemgroup_per_image <- lapply(spots_per_image, function(cells) {
  unique(hippo$GEMgroup[cells])
})

GEMGROUPS_PER_IMAGE <- lengths(gemgroup_per_image)
stopifnot(all(GEMGROUPS_PER_IMAGE == 1))

image_gemgroup_map <- data.frame(
  image = names(gemgroup_per_image),
  GEMgroup = unlist(gemgroup_per_image),
  Age = unname(unlist(lapply(gemgroup_per_image, function(g) {
    unique(hippo$Age[hippo$GEMgroup == g])
  }))),
  n_spots = unname(unlist(lapply(spots_per_image, length))),
  stringsAsFactors = FALSE
)
stopifnot(sum(image_gemgroup_map$n_spots) == 9921)

write.csv(image_gemgroup_map, file.path(out_dir, "image_gemgroup_mapping.csv"), row.names = FALSE)
```

STOP if any image maps to >1 GEMgroup, or any GEMgroup spans >1 image (SC-2).

### Step 3: Per-Image Spot Counts from RDS

Per-image spot counts are derived from RDS image-cell membership (Step 2), verified against Phase 02 output. Never hardcoded.

```r
phase02_counts <- read.csv("data/processed/spatial/phase02_metadata/Hippo/hippo_image_spot_counts.csv")
stopifnot(sum(phase02_counts$n_spots) == 9921)
# Cross-check runtime counts against Phase 02 counts
```

STOP if counts cannot be derived or do not sum to 9921 (SC-3).

### Step 4: Phase 04 WARNING Gate (SC-1)

Read Phase 04 validation summary. STOP unless user has acknowledged WARNING.

```r
phase04_val <- readLines("data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt", warn = FALSE)
# Must contain "USER ACKNOWLEDGED" or equivalent marker
```

### Step 5: Load Phase 04 Handoff

```r
func_region <- read.csv("data/processed/spatial/phase04_ifng_module/metadata_hippo_func_region.csv")
stopifnot(all(c("barcode", "FuncRegion", "is_edge") %in% colnames(func_region)))
stopifnot(nrow(func_region) == 9921)

hippo$FuncRegion <- func_region$FuncRegion[match(Cells(hippo), func_region$barcode)]
hippo$is_edge    <- func_region$is_edge[match(Cells(hippo), func_region$barcode)]

edge_spots_all <- Cells(hippo)[hippo$is_edge == TRUE]
cat(sprintf("Phase 04 handoff: %d edge seeds loaded\n", length(edge_spots_all)))
```

### Step 6: Per-Image Coordinate Neighbor Graph

For each of 16 images:

**6.1 Extract coordinates**

```r
coords <- Seurat::GetTissueCoordinates(hippo, image = img_name)
# -> data.frame with columns: imagecol, imagerow
# rownames = cell barcodes belonging to this image
```

STOP if `GetTissueCoordinates()` returns 0 rows or fails.

**6.2 Compute characteristic NN distance**

```r
dmat <- as.matrix(dist(coords[, c("imagecol", "imagerow")]))
diag(dmat) <- NA
nn_dists <- apply(dmat, 1, min, na.rm = TRUE)

# Mode-based characteristic distance
dens <- density(nn_dists)
char_dist <- dens$x[which.max(dens$y)]
margin   <- 1.2  # fixed geometric margin
thresh   <- char_dist * margin
```

**6.3 Build adjacency graph**

```r
adj <- (dmat <= thresh & dmat > 0)
rownames(adj) <- rownames(coords)
colnames(adj) <- rownames(coords)
```

**6.4 Ring expansion**

```r
# Identify Edge seeds for this image
edge_spots <- Cells(hippo)[hippo$is_edge & Cells(hippo) %in% rownames(coords)]

if (length(edge_spots) == 0) {
  # Record zero for this image, continue
  label_df <- data.frame(
    barcode = rownames(coords),
    spot_label = "unlabelled",
    image = img_name,
    stringsAsFactors = FALSE
  )
  image_labels[[img_name]] <- label_df
  next  # skip ring expansion for this image
}

# Edge = Phase 04 seeds (already labelled)
edge_set <- edge_spots

# nbs_Edge = immediate graph neighbors of Edge
nbs_neighbors <- unique(unlist(lapply(edge_set, function(sp) {
  names(which(adj[sp, ]))
})))
nbs_edge_set <- setdiff(nbs_neighbors, edge_set)

if (length(nbs_edge_set) > 0) {
  # RS1 = graph neighbors of nbs_Edge, excluding Edge and nbs_Edge
  rs1_neighbors <- unique(unlist(lapply(nbs_edge_set, function(sp) {
    names(which(adj[sp, ]))
  })))
  rs1_set <- setdiff(rs1_neighbors, c(edge_set, nbs_edge_set))
} else {
  rs1_set <- character(0)
}

if (length(rs1_set) > 0) {
  # RS2 = graph neighbors of RS1, excluding Edge, nbs_Edge, RS1
  rs2_neighbors <- unique(unlist(lapply(rs1_set, function(sp) {
    names(which(adj[sp, ]))
  })))
  rs2_set <- setdiff(rs2_neighbors, c(edge_set, nbs_edge_set, rs1_set))
} else {
  rs2_set <- character(0)
}

# Assemble labels for this image
spot_label_vec <- rep("unlabelled", nrow(coords))
names(spot_label_vec) <- rownames(coords)
spot_label_vec[edge_set]     <- "Edge"
spot_label_vec[nbs_edge_set] <- "nbs_Edge"
spot_label_vec[rs1_set]      <- "RS1"
spot_label_vec[rs2_set]      <- "RS2"

label_df <- data.frame(
  barcode = names(spot_label_vec),
  spot_label = spot_label_vec,
  image = img_name,
  stringsAsFactors = FALSE
)
image_labels[[img_name]] <- label_df
```

Record per-image label counts and graph parameters (char_dist, margin, thresh) to `phase05_neighbor_param_log.csv`.

### Step 7: Merge Per-Image Labels

```r
all_labels <- do.call(rbind, image_labels)
stopifnot(nrow(all_labels) == 9921)

hippo$spot_label <- all_labels$spot_label[match(Cells(hippo), all_labels$barcode)]
```

### Step 8: Create NBS, Group, FuncRegion

```r
# Author-style labels: Edge, nbs_Edge, RS1, RS2 -> final NBS
hippo$NBS <- dplyr::recode(
  hippo$spot_label,
  "Edge"     = "IS",
  "nbs_Edge" = "NNS",
  "RS1"      = "RS1",
  "RS2"      = "RS2"
)

# Group = GEMgroup_NBS (sample_id for pseudobulk)
hippo$Group <- paste0(hippo$GEMgroup, "_", hippo$NBS)

# Final FuncRegion for DESeq2 group_id
hippo$FuncRegion <- dplyr::recode(
  hippo$spot_label,
  "Edge"     = "IS",
  "nbs_Edge" = "NNS",
  "RS1"      = "ENS1",
  "RS2"      = "ENS2"
)
```

**STOP if** Old samples have zero Edge/IS seeds (SC-5).
**STOP if** Old samples have zero nbs_Edge/NNS labels (SC-6).

### Step 9: Pseudobulk DESeq2 on Old (IS vs NNS)

**9.1 Subset to labelled spots**

```r
hippo_labeled <- subset(
  hippo,
  subset = spot_label %in% c("Edge", "nbs_Edge", "RS1", "RS2")
)
```

**9.2 Synchronize counts and metadata**

```r
labeled_cells <- Cells(hippo_labeled)
counts_labeled <- counts_spatial[, labeled_cells, drop = FALSE]
metadata_labeled <- hippo_labeled@meta.data

stopifnot(identical(colnames(counts_labeled), rownames(metadata_labeled)))  # SC-8
stopifnot(ncol(counts_labeled) == nrow(metadata_labeled))                  # SC-9
```

**9.3 Build SingleCellExperiment**

```r
metadata_labeled$cluster_id <- factor(metadata_labeled$Age)
metadata_labeled$group_id   <- factor(metadata_labeled$FuncRegion)
metadata_labeled$sample_id  <- factor(metadata_labeled$Group)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts_labeled),
  colData = metadata_labeled
)
```

**9.4 Pseudobulk aggregation (author Script 8 aligned)**

```r
groups <- SingleCellExperiment::colData(sce)[, c("cluster_id", "sample_id")]

pb <- Matrix.utils::aggregate.Matrix(
  t(SummarizedExperiment::assay(sce, "counts")),
  groupings = groups,
  fun = "sum"
)
```

**9.5 Split by cluster_id (Age)**

```r
splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 2), `[`, 1)

pb_list <- split.data.frame(pb, factor(splitf))
pb_list <- lapply(pb_list, function(u) {
  out <- t(u)
  colnames(out) <- stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")
  out
})
pb <- pb_list
```

STOP if `pb[["Old"]]` is absent (SC-11).

```r
counts_old <- pb[["Old"]]
cat(sprintf("Pseudobulk Old: %d genes x %d samples\n", nrow(counts_old), ncol(counts_old)))
```

**9.6 Build sample-level metadata**

```r
sample_meta <- metadata_labeled[, c("cluster_id", "sample_id", "group_id")]
sample_meta <- sample_meta[!duplicated(sample_meta$sample_id), , drop = FALSE]
rownames(sample_meta) <- as.character(sample_meta$sample_id)

cluster_metadata <- sample_meta[
  colnames(counts_old),
  c("cluster_id", "sample_id", "group_id"),
  drop = FALSE
]

stopifnot(identical(rownames(cluster_metadata), colnames(counts_old)))  # SC-12
stopifnot(all(cluster_metadata$cluster_id == "Old"))                    # SC-13
```

**9.7 DESeq2**

```r
dds_old <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_old,
  colData   = cluster_metadata,
  design    = ~ group_id
)

# Set reference level to NNS
dds_old$group_id <- stats::relevel(dds_old$group_id, ref = "NNS")

dds_old <- DESeq2::DESeq(dds_old)

# IS vs NNS contrast
res_is_nns <- DESeq2::results(dds_old, contrast = c("group_id", "IS", "NNS"))
res_is_nns <- res_is_nns[order(res_is_nns$pvalue), ]

cat(sprintf("DESeq2 IS vs NNS: %d genes tested\n", nrow(res_is_nns)))
cat(sprintf("  padj < 0.05: %d up, %d down\n",
  sum(res_is_nns$padj < 0.05 & res_is_nns$log2FoldChange > 0, na.rm = TRUE),
  sum(res_is_nns$padj < 0.05 & res_is_nns$log2FoldChange < 0, na.rm = TRUE)
))
```

### Step 10: Figures

- **PCA**: vst(dds_old), plotPCA colored by group_id
- **Heatmap**: top 50 DEGs by padj, scaled by row
- **Volcano**: IS vs NNS with EnhancedVolcano or ggplot2 + ggrepel

All figures save PNG + best-effort PDF. Labelled "RDS-approximated".

### Step 11: Validation

Record and compare against author Script 8 reference values:
- Per-image label counts (reference: line 135 has O_1_A1 RS2=129, RS1=140, nbs_Edge=119, Edge=45)
- IS/NNS/ENS1/ENS2 counts by Age
- DESeq2 IS vs NNS DEG counts (reference: Script 8 reports 1326 Down, 2357 Up)

**Do NOT interpret biology.** Write: "DESeq2 found X DEGs vs author reported Y at padj < 0.05. Figure layout corresponds to Script 8 lines 845–993."

### Step 12: Cleanup

- `rm()` all temporary objects
- `gc()` after each major operation
- No `.RData` — no `save.image()`
- Full spot-level counts remain sparse

---

## Stop Conditions

| ID | Condition | Action |
|---|---|---|
| SC-1 | Phase 04 WARNING not acknowledged | STOP |
| SC-2 | Image ↔ GEMgroup mapping not clean 1:1 | STOP |
| SC-3 | Per-image spot counts cannot be derived or do not sum to 9921 | STOP |
| SC-4 | Spatial counts not 32285 × 9921 | STOP |
| SC-5 | Old samples have zero Edge/IS seeds | STOP |
| SC-6 | Old samples have zero nbs_Edge/NNS labels | STOP |
| SC-7 | Required package fails to load | STOP; do NOT modify renv.lock |
| SC-8 | `identical(colnames(counts_labeled), rownames(metadata_labeled))` is FALSE | STOP |
| SC-9 | `counts_labeled` and `metadata_labeled` have different dimensions | STOP |
| SC-10 | Pseudobulk rownames cannot be split into cluster/sample IDs | STOP |
| SC-11 | `pb[["Old"]]` is absent | STOP |
| SC-12 | `cluster_metadata` not exactly aligned to `counts_old` | STOP |
| SC-13 | `cluster_metadata$cluster_id` not all "Old" | STOP |
| SC-14 | `GetTissueCoordinates()` fails or returns empty for any image | STOP |

---

## Memory Management

- Peak estimated memory: ~2–3 GB
- Sequential per-image processing; never hold all 16 distance matrices simultaneously
- Load only `seurat_Visium_Hippo_All.rds` (no WholeBrain/DG/Chromium)
- Full spot-level counts remain sparse (`dgCMatrix`)
- Per-image distance matrix (≤756 × 756) computed via `dist()`

---

## Route Decision

| Route | Status |
|---|---|
| Strict author-code | BLOCKED (STutility missing; 0/32 raw Visium files) |
| RDS-based approximation | **SELECTED / PROPOSED** |

---

## Implementation

Implementation via `R/spatial/s05_rds_inflammatory_gradient.R`.

Save this approved plan as `docs/spatial_phase05_rds_gradient_plan.md` before execution.

---

*Plan v3 — Last updated 2026-06-03*
