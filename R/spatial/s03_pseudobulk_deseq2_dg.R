#!/usr/bin/env Rscript
# Phase Spatial-03: Pseudobulk DESeq2 for DG (Old vs Young)
#
# Author: Claude (reproduction of Script 7: Spatial-DESeq2.R)
# Date:   2026-06-03
#
# Overview:
#   Load Hippo_All.rds, build metadata with ML/GCL/Hilus→DG recode,
#   aggregate Spatial raw counts to pseudobulk, subset to DG,
#   run DESeq2 with 3 contrasts (OY, MY, OM), output tables and figures.
#
# Input:
#   data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds
#
# Output:
#   data/processed/spatial/phase03_pseudobulk_deseq2/ (13 required files)
#   figures/spatial/phase03/ (7 figure PNGs + best-effort PDFs)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(SingleCellExperiment)
  library(purrr)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ── CONFIG ──────────────────────────────────────────────────────────────────────
INPUT_RDS      <- "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"
OUT_TABLES     <- "data/processed/spatial/phase03_pseudobulk_deseq2"
OUT_FIGURES    <- "figures/spatial/phase03"

EXPECTED_SPOTS <- 9921
EXPECTED_GENES_SCT <- 17096
EXPECTED_GENES_SPATIAL <- 32285

# Author Script 7 colors: Young=#c6dbef, Middle=#4292c6, Old=#08306b
AGE_COLORS <- c("Young" = "#c6dbef", "Middle" = "#4292c6", "Old" = "#08306b")

# ── HELPERS ─────────────────────────────────────────────────────────────────────
dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIGURES, recursive = TRUE, showWarnings = FALSE)

start_time <- Sys.time()
r_version  <- paste0(R.version$major, ".", R.version$minor)

safe_pdf <- function(filename, width = 7, height = 7, ...) {
  tryCatch(
    pdf(filename, width = width, height = height, ...),
    error = function(e) {
      message(sprintf("  PDF failed for %s: %s. Skipping PDF.", basename(filename), e$message))
      NULL
    }
  )
}

message("=== Phase Spatial-03: Pseudobulk DESeq2 for DG ===")
message(sprintf("Start time: %s", start_time))
message(sprintf("R version: %s", r_version))
message(sprintf("Seurat version: %s", packageVersion("Seurat")))
message(sprintf("DESeq2 version: %s", packageVersion("DESeq2")))
message(sprintf("Matrix version: %s", packageVersion("Matrix")))

# ── STEP 1: LOAD AND VALIDATE ──────────────────────────────────────────────────
message("\n--- Step 1: Loading Hippo object ---")

stopifnot(file.exists(INPUT_RDS))
file_size <- file.size(INPUT_RDS)
file_mtime <- file.info(INPUT_RDS)$mtime
message(sprintf("  Input: %s", INPUT_RDS))
message(sprintf("  File size: %s bytes (%.2f GB)", format(file_size, big.mark = ","), file_size / 1e9))
message(sprintf("  File mtime: %s", file_mtime))

hippo <- readRDS(INPUT_RDS)

# Validate Seurat class
stopifnot(inherits(hippo, "Seurat"))

# Validate assays
assay_names <- names(hippo@assays)
stopifnot("SCT" %in% assay_names)
stopifnot("Spatial" %in% assay_names)

# Validate dimensions
n_sct <- nrow(GetAssayData(hippo, assay = "SCT", layer = "counts"))
n_spatial <- nrow(GetAssayData(hippo, assay = "Spatial", layer = "counts"))
n_spots <- ncol(hippo)

message(sprintf("  Default assay: %s", DefaultAssay(hippo)))
message(sprintf("  SCT features: %d (expected %d)", n_sct, EXPECTED_GENES_SCT))
message(sprintf("  Spatial features: %d (expected %d)", n_spatial, EXPECTED_GENES_SPATIAL))
message(sprintf("  Spots: %d (expected %d)", n_spots, EXPECTED_SPOTS))

# STOP if Spatial nrow mismatch
if (n_spatial != EXPECTED_GENES_SPATIAL) {
  stop(sprintf("Spatial assay has %d genes, expected %d. STOP.", n_spatial, EXPECTED_GENES_SPATIAL))
}

# Validate metadata columns
meta_orig <- hippo@meta.data
stopifnot("Region" %in% colnames(meta_orig))
stopifnot("Age" %in% colnames(meta_orig))
stopifnot("GEMgroup" %in% colnames(meta_orig))

region_levels <- levels(factor(meta_orig$Region))
age_levels    <- levels(factor(meta_orig$Age))
gem_levels    <- sort(unique(as.integer(meta_orig$GEMgroup)))

message(sprintf("  Region levels: %s", paste(region_levels, collapse = ", ")))
message(sprintf("  Age levels: %s", paste(age_levels, collapse = ", ")))
message(sprintf("  GEMgroup: %s", paste(gem_levels, collapse = ", ")))

# Validate expected values
stopifnot(all(c("CA1", "CA2", "CA3", "ML", "GCL", "Hilus") %in% region_levels))
stopifnot(all(c("Young", "Middle", "Old") %in% age_levels))
stopifnot(length(gem_levels) == 16)
stopifnot(all(gem_levels == 1:16))

message("  Validation PASSED")
gc()

# ── STEP 2: BUILD METADATA ─────────────────────────────────────────────────────
message("\n--- Step 2: Building metadata with ML/GCL/Hilus -> DG recode ---")

meta <- data.frame(
  row.names = colnames(hippo),
  stringsAsFactors = FALSE
)

# cluster_id_final: ML/GCL/Hilus -> DG, CA1/CA2/CA3 stay
meta$cluster_id_final <- as.character(meta_orig$Region)
meta$cluster_id_final[meta$cluster_id_final %in% c("ML", "GCL", "Hilus")] <- "DG"
meta$cluster_id_final <- factor(meta$cluster_id_final, levels = c("CA1", "CA2", "CA3", "DG"))

# group_id: Age with explicit factor levels
meta$group_id <- factor(meta_orig$Age, levels = c("Young", "Middle", "Old"))

# sample_id: GEMgroup as character
meta$sample_id <- as.character(meta_orig$GEMgroup)

# pb_group_id: paste(cluster_id_final, sample_id, sep = "__")
meta$pb_group_id <- paste(meta$cluster_id_final, meta$sample_id, sep = "__")

# Audit
dg_spots <- sum(meta$cluster_id_final == "DG")
cat(sprintf("  cluster_id_final levels: %s\n", paste(levels(meta$cluster_id_final), collapse = ", ")))
cat(sprintf("  DG spots: %d (expected 2364)\n", dg_spots))
cat(sprintf("  group_id levels: %s\n", paste(levels(meta$group_id), collapse = ", ")))
cat(sprintf("  sample_id levels: %s\n", paste(sort(unique(meta$sample_id)), collapse = ", ")))

# Spot count audit
cat("\n  Spot counts by cluster_id_final:\n")
print(table(meta$cluster_id_final))

cat("\n  Spot counts by cluster_id_final × group_id:\n")
print(table(meta$cluster_id_final, meta$group_id))

cat("\n  Spot counts by cluster_id_final × sample_id:\n")
print(table(meta$cluster_id_final, meta$sample_id))

# Expected DG spot counts
dg_meta <- meta[meta$cluster_id_final == "DG", ]
dg_y <- sum(dg_meta$group_id == "Young")
dg_m <- sum(dg_meta$group_id == "Middle")
dg_o <- sum(dg_meta$group_id == "Old")
cat(sprintf("\n  DG spot counts: Young=%d (exp 797), Middle=%d (exp 519), Old=%d (exp 1048)\n", dg_y, dg_m, dg_o))

if (dg_y != 797 || dg_m != 519 || dg_o != 1048) {
  stop(sprintf("DG spot counts mismatch: Y=%d M=%d O=%d. STOP.", dg_y, dg_m, dg_o))
}

gc()

# ── STEP 3: EXTRACT SPATIAL COUNTS ─────────────────────────────────────────────
message("\n--- Step 3: Extracting Spatial@counts (sparse) ---")

counts <- GetAssayData(hippo, assay = "Spatial", layer = "counts")
message(sprintf("  Counts matrix: %d genes × %d spots", nrow(counts), ncol(counts)))
message(sprintf("  Class: %s", class(counts)))
message(sprintf("  Non-zero entries: %s", format(nnzero(counts), big.mark = ",")))

# Validate dimensions
stopifnot(nrow(counts) == EXPECTED_GENES_SPATIAL)
stopifnot(ncol(counts) == EXPECTED_SPOTS)

gc()

# ── STEP 4: CHECK aggregate.Matrix DEPENDENCY ──────────────────────────────────
message("\n--- Step 4: Checking aggregate.Matrix dependency ---")

has_aggregate_matrix <- requireNamespace("Matrix.utils", quietly = TRUE) &&
  exists("aggregate.Matrix", where = asNamespace("Matrix.utils"))

if (has_aggregate_matrix) {
  message("  Matrix.utils::aggregate.Matrix available - using author's method")
} else {
  stop(
    "Matrix.utils::aggregate.Matrix is not available. ",
    "Strict author reproduction requires aggregate.Matrix, so Phase Spatial-03 is blocked. ",
    "Do not install packages or use fallback without a separate approved dependency plan."
  )
}

# ── STEP 5: PSEUDOBULK AGGREGATION (HIPPO-LEVEL) ──────────────────────────────
message("\n--- Step 5: Pseudobulk aggregation (hippo-level) ---")

groups <- data.frame(
  pb_group_id = meta$pb_group_id,
  stringsAsFactors = FALSE
)

message("  Using Matrix.utils::aggregate.Matrix")
# aggregate.Matrix(t(counts), groupings, sum)
# Input: t(counts) = spots x genes
# Output: groups x genes (rows = pb_group_id)
pb_raw <- Matrix.utils::aggregate.Matrix(t(counts), groupings = groups, fun = "sum")
# Transpose: genes x groups
pb_hippo <- t(pb_raw)

# Validate pseudobulk dimensions
message(sprintf("  Pseudobulk matrix: %d genes × %d groups", nrow(pb_hippo), ncol(pb_hippo)))
stopifnot(nrow(pb_hippo) == EXPECTED_GENES_SPATIAL)

# Build manifest from original metadata (before alignment)
manifest_hippo <- meta %>%
  distinct(cluster_id_final, sample_id, group_id, pb_group_id) %>%
  arrange(sample_id, cluster_id_final)

spot_counts <- as.data.frame(table(meta$pb_group_id), stringsAsFactors = FALSE)
colnames(spot_counts) <- c("pb_group_id", "n_spots")
spot_counts$n_spots <- as.integer(spot_counts$n_spots)
manifest_hippo <- manifest_hippo %>%
  left_join(spot_counts, by = "pb_group_id")

# Validate: all pb_group_id in manifest match colnames(pb_hippo)
if (!setequal(manifest_hippo$pb_group_id, colnames(pb_hippo))) {
  stop("manifest_hippo pb_group_id does not match colnames(pb_hippo). STOP.")
}

# Align manifest to pb_hippo column order using match()
reorder_idx <- match(colnames(pb_hippo), manifest_hippo$pb_group_id)
manifest_hippo <- manifest_hippo[reorder_idx, ]
if (any(is.na(manifest_hippo$n_spots)) || any(manifest_hippo$n_spots <= 0)) {
  stop("manifest_hippo has missing or non-positive n_spots values. STOP.")
}

cat(sprintf("\n  Hippo manifest: %d pseudobulk samples\n", nrow(manifest_hippo)))
cat(sprintf("  pb_hippo dimensions: %d genes × %d groups\n", nrow(pb_hippo), ncol(pb_hippo)))
cat(sprintf("  Alignment check: %s\n",
  ifelse(identical(manifest_hippo$pb_group_id, colnames(pb_hippo)), "PASSED", "FAILED")))

# Save hippo-level intermediate
message("  Saving hippo-level pseudobulk intermediate...")
saveRDS(pb_hippo, file.path(OUT_TABLES, "pseudobulk_counts_hippo.rds"))
write.csv(manifest_hippo, file.path(OUT_TABLES, "manifest_hippo_pseudobulk.csv"), row.names = FALSE)
message(sprintf("    -> %s", file.path(OUT_TABLES, "pseudobulk_counts_hippo.rds")))
message(sprintf("    -> %s", file.path(OUT_TABLES, "manifest_hippo_pseudobulk.csv")))

gc()

# ── STEP 6: SUBSET TO DG ───────────────────────────────────────────────────────
message("\n--- Step 6: Subsetting to DG pseudobulk ---")

dg_idx <- manifest_hippo$cluster_id_final == "DG"
pb_dg <- pb_hippo[, dg_idx]
manifest_dg <- manifest_hippo[dg_idx, ]

# Validate DG dimensions
cat(sprintf("  DG pseudobulk: %d genes × %d samples\n", nrow(pb_dg), ncol(pb_dg)))
cat(sprintf("  Expected: %d genes × 16 samples\n", EXPECTED_GENES_SPATIAL))

if (ncol(pb_dg) != 16) {
  stop(sprintf("DG pseudobulk has %d samples, expected 16. STOP.", ncol(pb_dg)))
}

# Validate age distribution
dg_y_pb <- sum(manifest_dg$group_id == "Young")
dg_m_pb <- sum(manifest_dg$group_id == "Middle")
dg_o_pb <- sum(manifest_dg$group_id == "Old")
cat(sprintf("  DG age distribution: Young=%d (exp 4), Middle=%d (exp 4), Old=%d (exp 8)\n",
  dg_y_pb, dg_m_pb, dg_o_pb))

if (dg_y_pb != 4 || dg_m_pb != 4 || dg_o_pb != 8) {
  stop(sprintf("DG age distribution mismatch: Y=%d M=%d O=%d. STOP.", dg_y_pb, dg_m_pb, dg_o_pb))
}

# Validate alignment
if (!identical(manifest_dg$pb_group_id, colnames(pb_dg))) {
  stop("DG manifest not aligned with pb_dg colnames. STOP.")
}

cat("\n  DG manifest:\n")
print(manifest_dg)

# Save DG pseudobulk
saveRDS(pb_dg, file.path(OUT_TABLES, "pseudobulk_counts_dg.rds"))
write.csv(manifest_dg, file.path(OUT_TABLES, "manifest_dg.csv"), row.names = FALSE)
message(sprintf("    -> %s", file.path(OUT_TABLES, "pseudobulk_counts_dg.rds")))
message(sprintf("    -> %s", file.path(OUT_TABLES, "manifest_dg.csv")))

gc()

# ── STEP 7: DESeq2 ANALYSIS ────────────────────────────────────────────────────
message("\n--- Step 7: DESeq2 analysis (~ group_id, 3 contrasts) ---")

# Convert pseudobulk to dense matrix for DESeq2
# pb_dg is genes × samples, DESeq2 expects samples × genes
count_matrix <- as.matrix(pb_dg)
cat(sprintf("  Count matrix for DESeq2: %d genes × %d samples\n", nrow(count_matrix), ncol(count_matrix)))

# Build colData from manifest
col_data <- data.frame(
  row.names = colnames(count_matrix),
  group_id = manifest_dg$group_id,
  sample_id = manifest_dg$sample_id,
  stringsAsFactors = FALSE
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ group_id
)

# Set reference level to Young
dds$group_id <- relevel(dds$group_id, ref = "Young")

message("  Running DESeq2...")
dds <- DESeq(dds)

# Extract results for 3 contrasts
# 1. Old vs Young (OY)
res_OY <- results(dds, contrast = c("group_id", "Old", "Young"))
cat(sprintf("\n  Old vs Young: %d DEGs (padj < 0.05)\n", sum(res_OY$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  Old vs Young: %d DEGs (padj < 0.05, |log2FC| > %.2f)\n",
  sum(res_OY$padj < 0.05 & abs(res_OY$log2FoldChange) > log2(1.5), na.rm = TRUE), log2(1.5)))

# 2. Middle vs Young (MY)
res_MY <- results(dds, contrast = c("group_id", "Middle", "Young"))
cat(sprintf("  Middle vs Young: %d DEGs (padj < 0.05)\n", sum(res_MY$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  Middle vs Young: %d DEGs (padj < 0.05, |log2FC| > %.2f)\n",
  sum(res_MY$padj < 0.05 & abs(res_MY$log2FoldChange) > log2(1.5), na.rm = TRUE), log2(1.5)))

# 3. Old vs Middle (OM)
res_OM <- results(dds, contrast = c("group_id", "Old", "Middle"))
cat(sprintf("  Old vs Middle: %d DEGs (padj < 0.05)\n", sum(res_OM$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  Old vs Middle: %d DEGs (padj < 0.05, |log2FC| > %.2f)\n",
  sum(res_OM$padj < 0.05 & abs(res_OM$log2FoldChange) > log2(1.5), na.rm = TRUE), log2(1.5)))

# Expected: 32285 rows per contrast
cat(sprintf("  Expected rows per contrast: %d\n", EXPECTED_GENES_SPATIAL))
cat(sprintf("  Actual rows: OY=%d, MY=%d, OM=%d\n", nrow(res_OY), nrow(res_MY), nrow(res_OM)))

gc()

# ── STEP 8: DEG SUMMARY AND OUTPUT ─────────────────────────────────────────────
message("\n--- Step 8: DEG summary and output tables ---")

# Convert results to data frames
res_OY_df <- as.data.frame(res_OY)
res_OY_df$gene <- rownames(res_OY_df)
res_OY_df <- res_OY_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

res_MY_df <- as.data.frame(res_MY)
res_MY_df$gene <- rownames(res_MY_df)
res_MY_df <- res_MY_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

res_OM_df <- as.data.frame(res_OM)
res_OM_df$gene <- rownames(res_OM_df)
res_OM_df <- res_OM_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Save full results
write.csv(res_OY_df, file.path(OUT_TABLES, "deseq2_dg_OY_results.csv"), row.names = FALSE)
write.csv(res_MY_df, file.path(OUT_TABLES, "deseq2_dg_MY_results.csv"), row.names = FALSE)
write.csv(res_OM_df, file.path(OUT_TABLES, "deseq2_dg_OM_results.csv"), row.names = FALSE)
message(sprintf("    -> %s", file.path(OUT_TABLES, "deseq2_dg_OY_results.csv")))
message(sprintf("    -> %s", file.path(OUT_TABLES, "deseq2_dg_MY_results.csv")))
message(sprintf("    -> %s", file.path(OUT_TABLES, "deseq2_dg_OM_results.csv")))

# Two-tier DEG filtering
# Tier 1: padj < 0.05 (padj files)
deg_OY_padj <- res_OY_df[res_OY_df$padj < 0.05 & !is.na(res_OY_df$padj), ]
deg_MY_padj <- res_MY_df[res_MY_df$padj < 0.05 & !is.na(res_MY_df$padj), ]
deg_OM_padj <- res_OM_df[res_OM_df$padj < 0.05 & !is.na(res_OM_df$padj), ]

write.csv(deg_OY_padj, file.path(OUT_TABLES, "deg_OY_padj.csv"), row.names = FALSE)
write.csv(deg_MY_padj, file.path(OUT_TABLES, "deg_MY_padj.csv"), row.names = FALSE)
write.csv(deg_OM_padj, file.path(OUT_TABLES, "deg_OM_padj.csv"), row.names = FALSE)

# Tier 2: padj < 0.05 AND |log2FC| > log2(1.5) (strong files)
lfc_threshold <- log2(1.5)
deg_OY_strong <- deg_OY_padj[abs(deg_OY_padj$log2FoldChange) > lfc_threshold, ]
deg_MY_strong <- deg_MY_padj[abs(deg_MY_padj$log2FoldChange) > lfc_threshold, ]
deg_OM_strong <- deg_OM_padj[abs(deg_OM_padj$log2FoldChange) > lfc_threshold, ]

write.csv(deg_OY_strong, file.path(OUT_TABLES, "deg_OY_strong.csv"), row.names = FALSE)
write.csv(deg_MY_strong, file.path(OUT_TABLES, "deg_MY_strong.csv"), row.names = FALSE)
write.csv(deg_OM_strong, file.path(OUT_TABLES, "deg_OM_strong.csv"), row.names = FALSE)

# Summary table
deg_summary <- data.frame(
  contrast = c("OY", "MY", "OM"),
  n_padj = c(nrow(deg_OY_padj), nrow(deg_MY_padj), nrow(deg_OM_padj)),
  n_strong = c(nrow(deg_OY_strong), nrow(deg_MY_strong), nrow(deg_OM_strong)),
  lfc_threshold = rep(lfc_threshold, 3),
  stringsAsFactors = FALSE
)
write.csv(deg_summary, file.path(OUT_TABLES, "deg_summary.csv"), row.names = FALSE)

cat("\n  DEG Summary:\n")
print(deg_summary)

gc()

# ── STEP 9: FIGURES ────────────────────────────────────────────────────────────
message("\n--- Step 9: Figures (best-effort) ---")

# 9a: PCA on rlog-transformed counts
message("  9a: PCA plot...")
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "group_id", returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group_id)) +
  geom_point(size = 3) +
  scale_color_manual(values = AGE_COLORS) +
  labs(
    title = "PCA of DG Pseudobulk Samples",
    x = sprintf("PC1 (%d%% variance)", pca_var[1]),
    y = sprintf("PC2 (%d%% variance)", pca_var[2]),
    color = "Age"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_FIGURES, "01_pca_dg.png"), p_pca, width = 6, height = 5, dpi = 150)
message(sprintf("    -> %s", file.path(OUT_FIGURES, "01_pca_dg.png")))
tryCatch({
  pdf(file.path(OUT_FIGURES, "01_pca_dg.pdf"), width = 6, height = 5)
  print(p_pca)
  dev.off()
  message(sprintf("    -> %s", file.path(OUT_FIGURES, "01_pca_dg.pdf")))
}, error = function(e) message(sprintf("    PDF failed: %s", e$message)))

gc()

# 9b: Correlation heatmap
message("  9b: Correlation heatmap...")
cor_mat <- cor(assay(vsd), method = "spearman")

# Annotation
annotation_col <- data.frame(
  Age = col_data$group_id,
  row.names = colnames(cor_mat)
)

png(file.path(OUT_FIGURES, "02_correlation_heatmap.png"), width = 800, height = 700, res = 150)
pheatmap(
  cor_mat,
  annotation_col = annotation_col,
  annotation_colors = list(Age = AGE_COLORS),
  color = colorRampPalette(brewer.pal(9, "Blues"))(100),
  main = "Spearman Correlation of DG Pseudobulk Samples",
  clustering_method = "ward.D2"
)
dev.off()
message(sprintf("    -> %s", file.path(OUT_FIGURES, "02_correlation_heatmap.png")))
tryCatch({
  pdf(file.path(OUT_FIGURES, "02_correlation_heatmap.pdf"), width = 8, height = 7)
  pheatmap(
    cor_mat,
    annotation_col = annotation_col,
    annotation_colors = list(Age = AGE_COLORS),
    color = colorRampPalette(brewer.pal(9, "Blues"))(100),
    main = "Spearman Correlation of DG Pseudobulk Samples",
    clustering_method = "ward.D2"
  )
  dev.off()
  message(sprintf("    -> %s", file.path(OUT_FIGURES, "02_correlation_heatmap.pdf")))
}, error = function(e) message(sprintf("    PDF failed: %s", e$message)))

gc()

# 9c: Dispersion plot
message("  9c: Dispersion plot...")
png(file.path(OUT_FIGURES, "03_dispersion.png"), width = 700, height = 600, res = 150)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
dev.off()
message(sprintf("    -> %s", file.path(OUT_FIGURES, "03_dispersion.png")))
tryCatch({
  pdf(file.path(OUT_FIGURES, "03_dispersion.pdf"), width = 7, height = 6)
  plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
  dev.off()
  message(sprintf("    -> %s", file.path(OUT_FIGURES, "03_dispersion.pdf")))
}, error = function(e) message(sprintf("    PDF failed: %s", e$message)))

gc()

# 9d-f: Volcano plots (3 contrasts)
make_volcano <- function(res_df, title, lfc_thresh = log2(1.5), padj_thresh = 0.05) {
  res_df$significant <- ifelse(
    !is.na(res_df$padj) & res_df$padj < padj_thresh & abs(res_df$log2FoldChange) > lfc_thresh,
    "Significant", "Not Significant"
  )
  res_df$label <- ""
  top_genes <- res_df[res_df$significant == "Significant", ]
  top_genes <- top_genes[order(top_genes$padj), ]
  if (nrow(top_genes) > 20) top_genes <- top_genes[1:20, ]
  res_df$label[res_df$gene %in% top_genes$gene] <- top_genes$gene

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "red")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)") +
    theme_classic() +
    theme(legend.position = "bottom")
}

message("  9d: Volcano OY...")
p_vol_OY <- make_volcano(res_OY_df, "DG: Old vs Young")
ggsave(file.path(OUT_FIGURES, "04_volcano_OY.png"), p_vol_OY, width = 7, height = 6, dpi = 150)
tryCatch({
  pdf(file.path(OUT_FIGURES, "04_volcano_OY.pdf"), width = 7, height = 6)
  print(p_vol_OY)
  dev.off()
}, error = function(e) message(sprintf("    PDF failed: %s", e$message)))

message("  9e: Volcano MY...")
p_vol_MY <- make_volcano(res_MY_df, "DG: Middle vs Young")
ggsave(file.path(OUT_FIGURES, "05_volcano_MY.png"), p_vol_MY, width = 7, height = 6, dpi = 150)
tryCatch({
  pdf(file.path(OUT_FIGURES, "05_volcano_MY.pdf"), width = 7, height = 6)
  print(p_vol_MY)
  dev.off()
}, error = function(e) message(sprintf("    PDF failed: %s", e$message)))

message("  9f: Volcano OM...")
p_vol_OM <- make_volcano(res_OM_df, "DG: Old vs Middle")
ggsave(file.path(OUT_FIGURES, "06_volcano_OM.png"), p_vol_OM, width = 7, height = 6, dpi = 150)
tryCatch({
  pdf(file.path(OUT_FIGURES, "06_volcano_OM.pdf"), width = 7, height = 6)
  print(p_vol_OM)
  dev.off()
}, error = function(e) message(sprintf("    PDF failed: %s", e$message)))

gc()

# 9g: DEG heatmap (top significant genes from OY contrast)
message("  9g: DEG heatmap...")
top_OY <- deg_OY_padj[order(deg_OY_padj$padj), ]
if (nrow(top_OY) > 50) top_OY <- top_OY[1:50, ]

if (nrow(top_OY) > 0) {
  # Get normalized counts for top genes
  norm_counts <- counts(dds, normalized = TRUE)
  top_genes <- top_OY$gene[top_OY$gene %in% rownames(norm_counts)]
  heatmap_mat <- norm_counts[top_genes, ]

  # Scale rows
  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

  annotation_col_hm <- data.frame(
    Age = col_data$group_id,
    row.names = colnames(heatmap_mat_scaled)
  )

  png(file.path(OUT_FIGURES, "07_deg_heatmap.png"), width = 800, height = 900, res = 150)
  pheatmap(
    heatmap_mat_scaled,
    annotation_col = annotation_col_hm,
    annotation_colors = list(Age = AGE_COLORS),
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    main = sprintf("Top %d DEGs (OY, padj < 0.05)", length(top_genes)),
    clustering_method = "ward.D2",
    show_rownames = length(top_genes) <= 50
  )
  dev.off()
  message(sprintf("    -> %s", file.path(OUT_FIGURES, "07_deg_heatmap.png")))
  tryCatch({
    pdf(file.path(OUT_FIGURES, "07_deg_heatmap.pdf"), width = 8, height = 9)
    pheatmap(
      heatmap_mat_scaled,
      annotation_col = annotation_col_hm,
      annotation_colors = list(Age = AGE_COLORS),
      color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
      main = sprintf("Top %d DEGs (OY, padj < 0.05)", length(top_genes)),
      clustering_method = "ward.D2",
      show_rownames = length(top_genes) <= 50
    )
    dev.off()
    message(sprintf("    -> %s", file.path(OUT_FIGURES, "07_deg_heatmap.pdf")))
  }, error = function(e) message(sprintf("    PDF failed: %s", e$message)))
} else {
  message("    No significant DEGs for heatmap. Skipping.")
}

gc()

# ── STEP 10: PROVENANCE AND VALIDATION ─────────────────────────────────────────
message("\n--- Step 10: Provenance and validation ---")

# Save provenance
provenance <- data.frame(
  key = c(
    "input_path",
    "input_file_size_bytes",
    "input_file_mtime",
    "r_version",
    "seurat_version",
    "deseq2_version",
    "matrix_version",
    "pheatmap_version",
    "n_spots_total",
    "n_genes_spatial",
    "dg_spots_young",
    "dg_spots_middle",
    "dg_spots_old",
    "dg_spots_total",
    "dg_pseudobulk_samples",
    "pseudobulk_method",
    "design_formula",
    "contrasts",
    "lfc_threshold",
    "padj_threshold",
    "checksum_note"
  ),
  value = c(
    INPUT_RDS,
    as.character(file_size),
    as.character(file_mtime),
    r_version,
    as.character(packageVersion("Seurat")),
    as.character(packageVersion("DESeq2")),
    as.character(packageVersion("Matrix")),
    as.character(packageVersion("pheatmap")),
    as.character(EXPECTED_SPOTS),
    as.character(EXPECTED_GENES_SPATIAL),
    as.character(dg_y),
    as.character(dg_m),
    as.character(dg_o),
    as.character(dg_spots),
    as.character(ncol(pb_dg)),
    "Matrix.utils::aggregate.Matrix",
    "~ group_id",
    "OY, MY, OM",
    as.character(round(lfc_threshold, 4)),
    "0.05",
    "checksum skipped unless explicitly requested"
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(OUT_TABLES, "phase03_provenance.csv"), row.names = FALSE)
message(sprintf("    -> %s", file.path(OUT_TABLES, "phase03_provenance.csv")))

# Final validation
message("\n=== FINAL VALIDATION ===")
cat("Required output files:\n")
required_files <- c(
  "pseudobulk_counts_hippo.rds",
  "manifest_hippo_pseudobulk.csv",
  "pseudobulk_counts_dg.rds",
  "manifest_dg.csv",
  "deseq2_dg_OY_results.csv",
  "deseq2_dg_MY_results.csv",
  "deseq2_dg_OM_results.csv",
  "deg_summary.csv",
  "deg_OY_padj.csv",
  "deg_MY_padj.csv",
  "deg_OM_padj.csv",
  "deg_OY_strong.csv",
  "deg_MY_strong.csv",
  "deg_OM_strong.csv",
  "phase03_provenance.csv"
)

all_exist <- TRUE
for (f in required_files) {
  fpath <- file.path(OUT_TABLES, f)
  exists <- file.exists(fpath)
  if (!exists) all_exist <- FALSE
  cat(sprintf("  %s %s\n", ifelse(exists, "PASS", "FAIL"), f))
}

cat("\nFigure files:\n")
figure_files <- c(
  "01_pca_dg.png",
  "02_correlation_heatmap.png",
  "03_dispersion.png",
  "04_volcano_OY.png",
  "05_volcano_MY.png",
  "06_volcano_OM.png",
  "07_deg_heatmap.png"
)

for (f in figure_files) {
  fpath <- file.path(OUT_FIGURES, f)
  exists <- file.exists(fpath)
  if (!exists) all_exist <- FALSE
  cat(sprintf("  %s %s\n", ifelse(exists, "PASS", "FAIL"), f))
}

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat(sprintf("Total genes in analysis: %d\n", EXPECTED_GENES_SPATIAL))
cat(sprintf("DG pseudobulk samples: %d\n", ncol(pb_dg)))
cat(sprintf("DEG summary (padj < 0.05): OY=%d, MY=%d, OM=%d\n",
  nrow(deg_OY_padj), nrow(deg_MY_padj), nrow(deg_OM_padj)))
cat(sprintf("DEG summary (strong): OY=%d, MY=%d, OM=%d\n",
  nrow(deg_OY_strong), nrow(deg_MY_strong), nrow(deg_OM_strong)))

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
cat(sprintf("\nElapsed time: %.2f minutes\n", elapsed))
cat("=== Phase Spatial-03 COMPLETE ===\n")
