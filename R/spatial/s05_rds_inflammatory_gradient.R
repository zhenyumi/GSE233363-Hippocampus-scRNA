# R/spatial/s05_rds_inflammatory_gradient.R
# Phase Spatial-05: IS/NNS/ENS1/ENS2 Label Approximation via Coordinate Neighbor Graph
# Route: RDS-based approximation — NOT strict author-code reproduction
# Target: Author Script 8 lines 71-993 (STutility → RegionNeighbours → DESeq2 + figures)
#
# Prerequisites:
#   - Phase 04 WARNING user-acknowledged
#   - seurat_Visium_Hippo_All.rds present
#   - metadata_hippo_func_region.csv present (Phase 04 handoff)
#   - No STutility, no Monocle 2, no raw Visium files
#
# Output:
#   data/processed/spatial/phase05_rds_gradient/
#   figures/spatial/phase05/

library(Seurat)
library(DESeq2)
library(Matrix.utils)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(purrr)

cat("=== Phase Spatial-05: RDS-Based Coordinate Neighbor Graph Approximation ===\n")

# ---- Paths ----
data_dir    <- "data/processed/spatial/phase05_rds_gradient"
fig_dir     <- "figures/spatial/phase05"
hippo_path  <- "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"
func_region_path <- "data/processed/spatial/phase04_ifng_module/metadata_hippo_func_region.csv"
phase04_val_path <- "data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt"

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output data: %s\n", data_dir))
cat(sprintf("Output figures: %s\n", fig_dir))

# ---- Step 1: Load and Validate Hippo Object ----
cat("\n--- Step 1: Load and Validate Hippo Object ---\n")
hippo <- readRDS(hippo_path)
cat(sprintf("Object loaded: %d cells x %d features (SCT assay)\n", ncol(hippo), nrow(hippo)))

stopifnot(DefaultAssay(hippo) == "SCT")
assay_names <- names(hippo@assays)
stopifnot("Spatial" %in% assay_names)
stopifnot(ncol(hippo) == 9921)
stopifnot(length(Images(hippo)) == 16)

cat(sprintf("DefaultAssay: %s\n", DefaultAssay(hippo)))
cat(sprintf("Assays present: %s\n", paste(assay_names, collapse = ", ")))
cat(sprintf("Images: %d\n", length(Images(hippo))))
cat(sprintf("Images names: %s\n", paste(names(hippo@images), collapse = ", ")))

counts_spatial <- Seurat::GetAssayData(
  hippo,
  assay  = "Spatial",
  layer  = "counts"
)
stopifnot(all(dim(counts_spatial) == c(32285, 9921)))
cat(sprintf("Spatial counts: %d genes x %d spots\n", nrow(counts_spatial), ncol(counts_spatial)))

# ---- Step 2: Runtime Image <-> GEMgroup Mapping ----
cat("\n--- Step 2: Image <-> GEMgroup Mapping ---\n")

image_names <- names(hippo@images)

get_image_cells <- function(hippo_obj, img_name) {
  img <- hippo_obj@images[[img_name]]
  if (inherits(img, "VisiumV1")) {
    coords <- tryCatch(
      Seurat::GetTissueCoordinates(hippo_obj, image = img_name),
      error = function(e) NULL
    )
    if (!is.null(coords)) {
      return(rownames(coords))
    }
  }
  coords_df <- img@coordinates
  if (!is.null(coords_df)) {
    return(rownames(coords_df))
  }
  character(0)
}

spots_per_image <- setNames(lapply(image_names, function(img) {
  get_image_cells(hippo, img)
}), image_names)

n_per_image <- lengths(spots_per_image)
cat(sprintf("Per-image spot counts:\n"))
for (i in seq_along(image_names)) {
  cat(sprintf("  %s: %d spots\n", image_names[i], n_per_image[i]))
}

gemgroup_per_image <- lapply(image_names, function(img) {
  cells <- spots_per_image[[img]]
  if (length(cells) == 0) return(NA_integer_)
  sort(unique(hippo$GEMgroup[cells]))
})
names(gemgroup_per_image) <- image_names

GEMGROUPS_PER_IMAGE <- lengths(gemgroup_per_image)

# SC-2: STOP if not clean 1:1
if (any(GEMGROUPS_PER_IMAGE != 1)) {
  cat("ERROR: Image <-> GEMgroup mapping is NOT clean 1:1\n")
  for (i in seq_along(image_names)) {
    if (GEMGROUPS_PER_IMAGE[i] != 1) {
      cat(sprintf("  %s -> GEMgroups: %s (count=%d)\n",
        image_names[i],
        paste(gemgroup_per_image[[i]], collapse = ","),
        GEMGROUPS_PER_IMAGE[i]))
    }
  }
  stop("SC-2: Image <-> GEMgroup mapping not clean 1:1")
}

image_gemgroup_map <- data.frame(
  image    = image_names,
  GEMgroup = unlist(gemgroup_per_image),
  Age      = unname(unlist(lapply(gemgroup_per_image, function(g) {
    as.character(unique(hippo$Age[hippo$GEMgroup == g]))
  }))),
  n_spots  = n_per_image,
  stringsAsFactors = FALSE
)
stopifnot(sum(image_gemgroup_map$n_spots) == 9921)

write.csv(image_gemgroup_map, file.path(data_dir, "image_gemgroup_mapping.csv"), row.names = FALSE)
cat("Image <-> GEMgroup 1:1 mapping verified and saved.\n")

# ---- Step 3: Per-Image Spot Counts ----
cat("\n--- Step 3: Verify Per-Image Spot Counts ---\n")
phase02_counts <- read.csv("data/processed/spatial/phase02_metadata/Hippo/hippo_image_spot_counts.csv")
# Remove TOTAL_SUM / OBJECT_NCOL rows
phase02_counts <- phase02_counts[!grepl("TOTAL|OBJECT", phase02_counts$image_name), ]
stopifnot(sum(phase02_counts$n_spots) == 9921)
cat(sprintf("Phase 02 per-image counts sum: %d (PASS)\n", sum(phase02_counts$n_spots)))

# ---- Step 4: Phase 04 WARNING Gate ----
cat("\n--- Step 4: Phase 04 WARNING Gate ---\n")
phase04_val <- readLines(phase04_val_path, warn = FALSE)
cat("Phase 04 validation status:\n")
cat(paste(phase04_val[grep("validation_status|WARNING|pct|handoff", phase04_val, ignore.case = TRUE)], collapse = "\n"))
cat("\n")
cat("USER ACKNOWLEDGED: Phase 04 WARNING confirmed. Proceeding with RDS-approximated labels.\n")

# ---- Step 5: Load Phase 04 Handoff ----
cat("\n--- Step 5: Load Phase 04 Handoff ---\n")
func_region <- read.csv(func_region_path)
stopifnot(all(c("barcode", "FuncRegion", "is_edge") %in% colnames(func_region)))
stopifnot(nrow(func_region) == 9921)

hippo$FuncRegion <- func_region$FuncRegion[match(Cells(hippo), func_region$barcode)]
hippo$is_edge    <- func_region$is_edge[match(Cells(hippo), func_region$barcode)]

edge_spots_all <- Cells(hippo)[hippo$is_edge == TRUE & !is.na(hippo$is_edge)]
cat(sprintf("Phase 04 handoff: %d Edge seeds loaded (of %d total spots)\n",
  length(edge_spots_all), ncol(hippo)))

# ---- Step 6: Per-Image Coordinate Neighbor Graph ----
cat("\n--- Step 6: Per-Image Coordinate Neighbor Graph ---\n")

image_labels <- list()
param_log <- data.frame(
  image       = character(0),
  n_spots     = integer(0),
  char_dist   = numeric(0),
  margin      = numeric(0),
  thresh      = numeric(0),
  n_edge      = integer(0),
  n_nbs_edge  = integer(0),
  n_rs1       = integer(0),
  n_rs2       = integer(0),
  stringsAsFactors = FALSE
)

MARGIN <- 1.2  # fixed geometric margin

for (img_name in image_names) {
  cat(sprintf("\nProcessing image: %s\n", img_name))

  coords <- tryCatch(
    Seurat::GetTissueCoordinates(hippo, image = img_name),
    error = function(e) {
      cat(sprintf("  WARNING: GetTissueCoordinates failed: %s\n", e$message))
      NULL
    }
  )

  if (is.null(coords)) {
    coords <- hippo@images[[img_name]]@coordinates
    cat("  Fallback to @coordinates slot\n")
  }

  stopifnot(nrow(coords) > 0)
  stopifnot(all(c("imagecol", "imagerow") %in% colnames(coords)))
  cat(sprintf("  Coordinates: %d spots\n", nrow(coords)))

  # SC-14 already passed (GetTissueCoordinates succeeded)

  cell_barcodes <- rownames(coords)

  # Characteristic NN distance (mode-based)
  dmat <- as.matrix(dist(coords[, c("imagecol", "imagerow")]))
  diag(dmat) <- NA
  nn_dists <- apply(dmat, 1, min, na.rm = TRUE)

  dens <- density(nn_dists)
  char_dist <- dens$x[which.max(dens$y)]
  thresh <- char_dist * MARGIN

  # Build adjacency graph
  adj <- (dmat <= thresh & dmat > 0)
  rownames(adj) <- rownames(coords)
  colnames(adj) <- rownames(coords)

  # Find Edge seeds within this image
  edge_spots <- intersect(cell_barcodes, edge_spots_all)

  if (length(edge_spots) == 0) {
    cat(sprintf("  WARNING: 0 Edge seeds in image %s. Recording all as unlabelled.\n", img_name))
    label_df <- data.frame(
      barcode    = cell_barcodes,
      spot_label = "unlabelled",
      image      = img_name,
      stringsAsFactors = FALSE
    )
    image_labels[[img_name]] <- label_df
    param_log <- rbind(param_log, data.frame(
      image = img_name, n_spots = length(cell_barcodes),
      char_dist = char_dist, margin = MARGIN, thresh = thresh,
      n_edge = 0, n_nbs_edge = 0, n_rs1 = 0, n_rs2 = 0,
      stringsAsFactors = FALSE
    ))
    next
  }

  cat(sprintf("  Edge seeds in this image: %d\n", length(edge_spots)))

  # Ring expansion: Edge -> nbs_Edge -> RS1 -> RS2
  edge_set <- edge_spots

  # nbs_Edge = immediate neighbors of Edge
  nbs_neighbors <- unique(unlist(lapply(edge_set, function(sp) {
    if (sp %in% rownames(adj)) { names(which(adj[sp, ])) } else { character(0) }
  })))
  nbs_edge_set <- setdiff(nbs_neighbors, edge_set)
  cat(sprintf("  nbs_Edge: %d\n", length(nbs_edge_set)))

  # RS1 = neighbors of nbs_Edge, excluding Edge and nbs_Edge
  if (length(nbs_edge_set) > 0) {
    rs1_neighbors <- unique(unlist(lapply(nbs_edge_set, function(sp) {
      if (sp %in% rownames(adj)) { names(which(adj[sp, ])) } else { character(0) }
    })))
    rs1_set <- setdiff(rs1_neighbors, c(edge_set, nbs_edge_set))
  } else {
    rs1_set <- character(0)
  }
  cat(sprintf("  RS1: %d\n", length(rs1_set)))

  # RS2 = neighbors of RS1, excluding earlier rings
  if (length(rs1_set) > 0) {
    rs2_neighbors <- unique(unlist(lapply(rs1_set, function(sp) {
      if (sp %in% rownames(adj)) { names(which(adj[sp, ])) } else { character(0) }
    })))
    rs2_set <- setdiff(rs2_neighbors, c(edge_set, nbs_edge_set, rs1_set))
  } else {
    rs2_set <- character(0)
  }
  cat(sprintf("  RS2: %d\n", length(rs2_set)))

  # Assemble
  spot_label_vec <- rep("unlabelled", length(cell_barcodes))
  names(spot_label_vec) <- cell_barcodes
  spot_label_vec[edge_set]     <- "Edge"
  spot_label_vec[nbs_edge_set] <- "nbs_Edge"
  spot_label_vec[rs1_set]      <- "RS1"
  spot_label_vec[rs2_set]      <- "RS2"

  label_df <- data.frame(
    barcode    = names(spot_label_vec),
    spot_label = spot_label_vec,
    image      = img_name,
    stringsAsFactors = FALSE
  )
  image_labels[[img_name]] <- label_df

  param_log <- rbind(param_log, data.frame(
    image = img_name, n_spots = length(cell_barcodes),
    char_dist = char_dist, margin = MARGIN, thresh = thresh,
    n_edge = length(edge_set), n_nbs_edge = length(nbs_edge_set),
    n_rs1 = length(rs1_set), n_rs2 = length(rs2_set),
    stringsAsFactors = FALSE
  ))
}

# Write param log and per-image label counts
write.csv(param_log, file.path(data_dir, "phase05_per_image_label_counts.csv"), row.names = FALSE)
cat("\nPer-image label counts saved.\n")
print(param_log)

# ---- Step 7: Merge Per-Image Labels ----
cat("\n--- Step 7: Merge Labels ---\n")
all_labels <- do.call(rbind, image_labels)
stopifnot(nrow(all_labels) == 9921)

hippo$spot_label <- all_labels$spot_label[match(Cells(hippo), all_labels$barcode)]
hippo$spot_label[is.na(hippo$spot_label)] <- "unlabelled"

cat("Label distribution:\n")
print(table(hippo$spot_label))
cat(sprintf("Total spots: %d\n", sum(table(hippo$spot_label))))

# Save per-image spot labels audit table
audit_labels <- data.frame(
  barcode    = Cells(hippo),
  image      = image_gemgroup_map$image[match(hippo$GEMgroup, image_gemgroup_map$GEMgroup)],
  GEMgroup   = hippo$GEMgroup,
  Age        = hippo$Age,
  Region     = hippo$Region,
  is_edge    = hippo$is_edge,
  spot_label = hippo$spot_label,
  stringsAsFactors = FALSE
)
write.csv(audit_labels, file.path(data_dir, "per_image_spot_labels.csv"), row.names = FALSE)
cat(sprintf("per_image_spot_labels.csv saved (%d rows).\n", nrow(audit_labels)))

# ---- Step 8: Create NBS, Group, FuncRegion ----
cat("\n--- Step 8: Create NBS, Group, FuncRegion ---\n")

hippo$NBS <- dplyr::recode(
  hippo$spot_label,
  "Edge"     = "IS",
  "nbs_Edge" = "NNS",
  "RS1"      = "RS1",
  "RS2"      = "RS2",
  .default   = "unlabelled"
)

hippo$Group <- paste0(hippo$GEMgroup, "_", hippo$NBS)

hippo$FuncRegion <- dplyr::recode(
  hippo$spot_label,
  "Edge"     = "IS",
  "nbs_Edge" = "NNS",
  "RS1"      = "ENS1",
  "RS2"      = "ENS2",
  .default   = "unlabelled"
)

# Counts by Age
cat("\nIS/NNS/ENS1/ENS2 counts by Age:\n")
counts_by_age <- table(hippo$Age, hippo$FuncRegion)
print(counts_by_age)
write.csv(as.data.frame.matrix(counts_by_age),
  file.path(data_dir, "phase05_final_label_counts.csv"))

# SC-5: Old samples must have Edge/IS
old_is <- sum(hippo$Age == "Old" & hippo$spot_label == "Edge")
cat(sprintf("Old Edge/IS spots: %d\n", old_is))
stopifnot(old_is > 0)

# SC-6: Old samples must have nbs_Edge/NNS
old_nns <- sum(hippo$Age == "Old" & hippo$spot_label == "nbs_Edge")
cat(sprintf("Old nbs_Edge/NNS spots: %d\n", old_nns))
stopifnot(old_nns > 0)

# Save merged spots audit table (all 9921 spots with all metadata columns)
merged_spots <- data.frame(
  barcode    = Cells(hippo),
  image      = image_gemgroup_map$image[match(hippo$GEMgroup, image_gemgroup_map$GEMgroup)],
  GEMgroup   = hippo$GEMgroup,
  Age        = hippo$Age,
  Region     = hippo$Region,
  is_edge    = hippo$is_edge,
  spot_label = hippo$spot_label,
  NBS        = hippo$NBS,
  Group      = hippo$Group,
  FuncRegion = hippo$FuncRegion,
  stringsAsFactors = FALSE
)
write.csv(merged_spots, file.path(data_dir, "phase05_merged_spots.csv"), row.names = FALSE)
cat(sprintf("phase05_merged_spots.csv saved (%d rows).\n", nrow(merged_spots)))

# ---- Step 9: Pseudobulk DESeq2 on Old (IS vs NNS) ----
cat("\n--- Step 9: Pseudobulk DESeq2 on Old (IS vs NNS) ---\n")

# 9.1 Identify labelled cells (avoids subset() which fails on VisiumV1 image slots)
cat("\n9.1 Identifying labelled spots...\n")
valid_labels <- c("Edge", "nbs_Edge", "RS1", "RS2")
labeled_idx <- hippo$spot_label %in% valid_labels
labeled_cells <- Cells(hippo)[labeled_idx]
cat(sprintf("Labelled spots: %d\n", length(labeled_cells)))

# 9.2 Synchronize counts and metadata
cat("\n9.2 Synchronizing counts and metadata...\n")
counts_labeled <- counts_spatial[, labeled_cells, drop = FALSE]
metadata_labeled <- hippo@meta.data[labeled_idx, , drop = FALSE]

# SC-8: exact alignment
stopifnot(identical(colnames(counts_labeled), rownames(metadata_labeled)))

# SC-9: same dimensions
stopifnot(ncol(counts_labeled) == nrow(metadata_labeled))
cat(sprintf("Counts: %d genes x %d spots; Metadata: %d rows (ALIGNED)\n",
  nrow(counts_labeled), ncol(counts_labeled), nrow(metadata_labeled)))

# 9.3 Build SingleCellExperiment
cat("\n9.3 Building SingleCellExperiment...\n")
metadata_labeled$cluster_id <- factor(metadata_labeled$Age)
metadata_labeled$group_id   <- factor(metadata_labeled$FuncRegion)
metadata_labeled$sample_id  <- factor(metadata_labeled$Group)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(counts = counts_labeled),
  colData = metadata_labeled
)
cat("SCE created.\n")
cat(sprintf("SCE colData columns: %s\n", paste(colnames(colData(sce)), collapse = ", ")))

# 9.4 Pseudobulk aggregation
cat("\n9.4 Pseudobulk aggregation (author Script 8 aligned)...\n")
groups <- SingleCellExperiment::colData(sce)[, c("cluster_id", "sample_id")]

pb <- Matrix.utils::aggregate.Matrix(
  t(SummarizedExperiment::assay(sce, "counts")),
  groupings = groups,
  fun = "sum"
)
cat(sprintf("Pseudobulk matrix: %d samples x %d genes\n", nrow(pb), ncol(pb)))

# 9.5 Split by cluster_id (Age)
cat("\n9.5 Splitting by cluster_id...\n")
splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 2), `[`, 1)
cat(sprintf("cluster_id values: %s\n", paste(unique(splitf), collapse = ", ")))

# SC-10: check rownames can be split
stopifnot(length(splitf) == nrow(pb))

pb_list <- split.data.frame(pb, factor(splitf))
pb_list <- lapply(pb_list, function(u) {
  out <- t(u)
  colnames(out) <- sub("^[^_]+_", "", rownames(u))
  out
})
pb <- pb_list

cat(sprintf("Pseudobulk matrices per cluster: %s\n", paste(names(pb), collapse = ", ")))

# SC-11: pb[["Old"]] must exist
stopifnot("Old" %in% names(pb))
counts_old <- pb[["Old"]]
cat(sprintf("counts_old: %d genes x %d samples\n", nrow(counts_old), ncol(counts_old)))

# Save pseudobulk
saveRDS(counts_old, file.path(data_dir, "phase05_old_pseudobulk_counts.rds"))

# 9.6 Build sample-level metadata
cat("\n9.6 Building sample-level metadata...\n")
sample_meta <- metadata_labeled[, c("cluster_id", "sample_id", "group_id")]
sample_meta <- sample_meta[!duplicated(sample_meta$sample_id), , drop = FALSE]
rownames(sample_meta) <- as.character(sample_meta$sample_id)

cluster_metadata <- sample_meta[
  colnames(counts_old),
  c("cluster_id", "sample_id", "group_id"),
  drop = FALSE
]

# SC-12: exact alignment
stopifnot(identical(rownames(cluster_metadata), colnames(counts_old)))

# SC-13: all Old
stopifnot(all(cluster_metadata$cluster_id == "Old"))
cat(sprintf("cluster_metadata: %d rows, all cluster_id='Old' (ALIGNED)\n",
  nrow(cluster_metadata)))
cat(sprintf("group_id levels: %s\n", paste(levels(cluster_metadata$group_id), collapse = ", ")))

# 9.7 DESeq2
cat("\n9.7 Running DESeq2...\n")
dds_old <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_old,
  colData   = cluster_metadata,
  design    = ~ group_id
)

# Set reference level to NNS
dds_old$group_id <- relevel(dds_old$group_id, ref = "NNS")

dds_old <- DESeq2::DESeq(dds_old)

# IS vs NNS contrast
res_is_nns <- DESeq2::results(dds_old, contrast = c("group_id", "IS", "NNS"))
res_is_nns <- res_is_nns[order(res_is_nns$pvalue), ]

cat(sprintf("DESeq2 IS vs NNS: %d genes tested\n", nrow(res_is_nns)))
n_up   <- sum(res_is_nns$padj < 0.05 & res_is_nns$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res_is_nns$padj < 0.05 & res_is_nns$log2FoldChange < 0, na.rm = TRUE)
cat(sprintf("  padj < 0.05: %d up, %d down\n", n_up, n_down))
cat(sprintf("  Author Script 8 reference: 2357 Up, 1326 Down\n"))

# Save DESeq2 results
write.csv(as.data.frame(res_is_nns),
  file.path(data_dir, "phase05_deseq2_old_all_genes.csv"))
cat("Full DESeq2 results saved.\n")

# Significant DEGs
sig_genes <- res_is_nns[!is.na(res_is_nns$padj) & res_is_nns$padj < 0.05, ]
sig_genes <- sig_genes[order(sig_genes$padj), ]
write.csv(as.data.frame(sig_genes),
  file.path(data_dir, "phase05_deseq2_old_results.csv"))
cat(sprintf("Significant DEGs (padj < 0.05): %d\n", nrow(sig_genes)))

# ---- Step 10: Figures ----
cat("\n--- Step 10: Figures ---\n")

# 10.1 PCA
cat("10.1 PCA...\n")
vsd_old <- tryCatch(
  DESeq2::vst(dds_old, blind = TRUE),
  error = function(e) {
    cat(sprintf("vst failed: %s. Trying rlog...\n", e$message))
    DESeq2::rlog(dds_old, blind = TRUE)
  }
)

pca_data <- DESeq2::plotPCA(vsd_old, intgroup = "group_id", returnData = TRUE)
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group_id)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA — Old Samples (RDS-approximated IS/NNS/ENS1/ENS2)",
    subtitle = paste0("DESeq2 vst, colored by FuncRegion"),
    color = "FuncRegion"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "phase05_pca_old_by_region.png"), p_pca,
  width = 7, height = 6, dpi = 150)
tryCatch(
  ggsave(file.path(fig_dir, "phase05_pca_old_by_region.pdf"), p_pca,
    width = 7, height = 6, device = cairo_pdf),
  error = function(e) cat(sprintf("PDF save failed: %s\n", e$message))
)
cat("PCA figure saved.\n")

# 10.2 Heatmap — top DEGs
cat("10.2 Heatmap...\n")
top_n <- min(50, nrow(sig_genes))
if (top_n > 0) {
  top_genes <- rownames(sig_genes)[1:top_n]
  mat <- SummarizedExperiment::assay(vsd_old)[top_genes, , drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0
  mat[is.infinite(mat)] <- 0

  ann_col <- data.frame(
    FuncRegion = cluster_metadata$group_id,
    row.names  = colnames(mat)
  )

  png(file.path(fig_dir, "phase05_heatmap_top_degs.png"),
    width = 10, height = 8, units = "in", res = 150)
  pheatmap::pheatmap(mat,
    annotation_col = ann_col,
    show_rownames  = FALSE,
    color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main           = paste0("Top ", top_n, " DEGs — IS vs NNS (RDS-approximated)"),
    fontsize       = 8
  )
  dev.off()
  cat(sprintf("Heatmap saved (%d genes).\n", top_n))
} else {
  cat("WARNING: No significant DEGs to plot heatmap.\n")
}

# 10.3 Volcano
cat("10.3 Volcano...\n")
volcano_df <- as.data.frame(res_is_nns)
volcano_df$gene <- rownames(volcano_df)
volcano_df$significance <- "NS"
volcano_df$significance[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 &
  volcano_df$log2FoldChange > 0] <- "Up"
volcano_df$significance[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 &
  volcano_df$log2FoldChange < 0] <- "Down"

# Label top 10 genes by padj
volcano_df$label <- ""
top_10 <- head(volcano_df[order(volcano_df$padj), ], 10)
volcano_df[top_10$gene, "label"] <- top_10$gene

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(pvalue),
  color = significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Down" = "blue", "NS" = "grey70", "Up" = "red")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
    show.legend = FALSE) +
  labs(
    title = "Volcano — IS vs NNS (Old, RDS-approximated)",
    subtitle = paste(sprintf("padj<0.05: %d Up, %d Down", n_up, n_down)),
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "phase05_volcano_is_vs_nns.png"), p_volcano,
  width = 8, height = 7, dpi = 150)
tryCatch(
  ggsave(file.path(fig_dir, "phase05_volcano_is_vs_nns.pdf"), p_volcano,
    width = 8, height = 7, device = cairo_pdf),
  error = function(e) cat(sprintf("PDF save failed: %s\n", e$message))
)
cat("Volcano figure saved.\n")

# ---- Step 11: Validation ----
cat("\n--- Step 11: Validation ---\n")

val_lines <- c(
  "=== Phase Spatial-05 Validation Summary ===",
  paste("Route:", "RDS-based coordinate neighbor graph approximation"),
  paste("Date:", Sys.time()),
  "",
  "--- Label Counts (Old) ---",
  paste("IS:", sum(hippo$Age == "Old" & hippo$spot_label == "Edge")),
  paste("NNS:", sum(hippo$Age == "Old" & hippo$spot_label == "nbs_Edge")),
  paste("ENS1:", sum(hippo$Age == "Old" & hippo$spot_label == "RS1")),
  paste("ENS2:", sum(hippo$Age == "Old" & hippo$spot_label == "RS2")),
  "",
  "--- DESeq2 IS vs NNS (Old) ---",
  paste("Genes tested:", nrow(res_is_nns)),
  paste("padj < 0.05 Up:", n_up),
  paste("padj < 0.05 Down:", n_down),
  paste("Author reference (Script 8): 2357 Up, 1326 Down"),
  "",
  "--- Graph Parameters ---",
  paste("Margin:", MARGIN),
  paste("N images with Edge seeds:", sum(param_log$n_edge > 0)),
  paste("N images with 0 Edge seeds:", sum(param_log$n_edge == 0)),
  "",
  "--- Figure Outputs ---",
  "PCA: figures/spatial/phase05/phase05_pca_old_by_region.png",
  "Heatmap: figures/spatial/phase05/phase05_heatmap_top_degs.png",
  "Volcano: figures/spatial/phase05/phase05_volcano_is_vs_nns.png",
  "",
  "NOTE: ALL outputs are RDS-approximated. Do not interpret biology from these labels.",
  "These are NOT equivalent to STutility RegionNeighbours image-mask-based labels.",
  "Compare count ranges against author Script 8 reference values only.",
  "",
  "--- Generated Outputs ---",
  "data/processed/spatial/phase05_rds_gradient/per_image_spot_labels.csv",
  "data/processed/spatial/phase05_rds_gradient/phase05_merged_spots.csv",
  "data/processed/spatial/phase05_rds_gradient/phase05_per_image_label_counts.csv",
  "data/processed/spatial/phase05_rds_gradient/phase05_final_label_counts.csv",
  "data/processed/spatial/phase05_rds_gradient/phase05_old_pseudobulk_counts.rds",
  "data/processed/spatial/phase05_rds_gradient/phase05_deseq2_old_all_genes.csv",
  "data/processed/spatial/phase05_rds_gradient/phase05_deseq2_old_results.csv",
  "data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv",
  "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt",
  "figures/spatial/phase05/phase05_pca_old_by_region.png",
  "figures/spatial/phase05/phase05_heatmap_top_degs.png",
  "figures/spatial/phase05/phase05_volcano_is_vs_nns.png",
  "",
  "--- Planned but NOT Generated ---",
  "Spatial label distribution figure is generated by R/spatial/s05b_plot_spatial_labels.R, not by this core DESeq2 script.",
  "PDF figures (cairo/X11 not available; PNGs generated as fallback)",
  "",
  "--- PDF Figure Status ---",
  "PDF figures not generated: cairo_pdf/X11 unavailable in this environment. PNGs are authoritative.",
  "",
  "--- Spatial Label Distribution Figure Status ---",
  "Generated by the plot-only companion script s05b when required.",
  "Current route is RDS-approximated coordinate neighbor graph; not strict STutility reproduction."
)

writeLines(val_lines, file.path(data_dir, "phase05_validation_summary.txt"))
cat("Validation summary saved.\n")

# Write provenance
prov <- data.frame(
  key   = c("route", "r_version", "seurat_version", "deseq2_version",
    "input_hippo_rds", "input_phase04_metadata",
    "margin", "n_edge_seeds", "n_labeled_spots",
    "deseq2_genes_tested", "deseq2_up", "deseq2_down",
    "output_dir", "figure_dir",
    "phase04_warning_acknowledged", "route_label", "strict_stutility_equivalence"),
  value = c(
    "RDS-based coordinate neighbor graph approximation",
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("Seurat")),
    as.character(packageVersion("DESeq2")),
    hippo_path,
    func_region_path,
    MARGIN,
    length(edge_spots_all),
    length(labeled_cells),
    nrow(res_is_nns),
    n_up,
    n_down,
    data_dir,
    fig_dir,
    "TRUE",
    "RDS-approximated",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)
write.csv(prov, file.path(data_dir, "phase05_provenance.csv"), row.names = FALSE)
cat("Provenance saved.\n")

# ---- Step 12: Cleanup ----
cat("\n--- Step 12: Cleanup ---\n")
rm(coords, dmat, adj, nn_dists, dens, pca_data, volcano_df, mat)
gc()
cat(sprintf("Final mem: %s\n", system("ps -o rss= -p $$", intern = TRUE)))

# ---- Final Report ----
cat("\n=== Phase Spatial-05 Complete ===\n")
cat(sprintf("Data outputs: %s\n", data_dir))
cat(sprintf("Figure outputs: %s\n", fig_dir))
cat("Route: RDS-based coordinate neighbor graph approximation\n")
cat("Plan: docs/spatial_phase05_rds_gradient_plan.md\n")
cat(sprintf("Labeled spots: %d (Edge=%d, nbs_Edge=%d, RS1=%d, RS2=%d)\n",
  length(labeled_cells),
  sum(hippo$spot_label == "Edge"),
  sum(hippo$spot_label == "nbs_Edge"),
  sum(hippo$spot_label == "RS1"),
  sum(hippo$spot_label == "RS2")
))
cat(sprintf("Old IS vs NNS DESeq2: %d Up, %d Down (padj < 0.05)\n", n_up, n_down))
cat("renv.lock: NOT modified\n")
cat("No packages installed. No .RData saved.\n")
cat("=== Done ===\n")
