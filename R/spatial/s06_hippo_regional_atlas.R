#!/usr/bin/env Rscript
# s06_hippo_regional_atlas.R
# Phase Spatial-06: Hippocampus Regional Atlas Reproduction
# RDS-based approximation; no PDF generation; no Seurat object modification
# Plan version: v2 (approved 2026-05-29)

# ============================================================================
# E0: SETUP
# ============================================================================

cat("=== Phase Spatial-06: Hippocampus Regional Atlas ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# -- Rplots.pdf guard (start) --
rplots_start <- file.exists("Rplots.pdf")
if (rplots_start) {
  cat("[SC-9 WARNING] Rplots.pdf exists at script start; will check again at end.\n")
}

# -- Git status --
git_status <- tryCatch({
  system("git status --short", intern = TRUE)
}, error = function(e) "git not available")
cat("Git status at start:", if (length(git_status) == 0) "(clean)" else paste(git_status, collapse = "; "), "\n\n")

# -- Required packages --
pkgs <- c("Seurat", "ggplot2", "dplyr", "cowplot", "patchwork")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[SC-8 STOP] Package not available:", pkg, "\n")
    stop(paste("Package", pkg, "not installed"))
  }
}
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)

# -- Paths --
hippo_rds  <- "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"
out_data   <- "data/processed/spatial/phase06_hippo_regional_atlas"
out_figs   <- "figures/spatial/phase06_hippo_regional_atlas"
dir.create(out_data, showWarnings = FALSE, recursive = TRUE)
dir.create(out_figs, showWarnings = FALSE, recursive = TRUE)

# -- Phase 02 reference paths --
p02_dir       <- "data/processed/spatial/phase02_metadata/Hippo"
p02_region    <- file.path(p02_dir, "hippo_region_distribution.csv")
p02_cross     <- "data/processed/spatial/phase02_metadata/cross_object_comparison.csv"
p02_inspect   <- "data/processed/spatial/inspection/Hippo/spatial_info.csv"

# -- Color scheme --
region_colors <- c(
  CA1   = "#1f77b4",
  CA2   = "#ff7f0e",
  CA3   = "#2ca02c",
  ML    = "#d62728",
  GCL   = "#9467bd",
  Hilus = "#8c564b"
)

# -- Phase 02 expected counts (from phase02_validation_summary.txt) --
expected_counts <- c(CA1 = 4671L, CA2 = 565L, CA3 = 2321L, ML = 1480L, GCL = 645L, Hilus = 239L)
expected_total  <- 9921L
expected_derived_dg <- 2364L  # ML + GCL + Hilus

# -- Provenance collection --
provenance <- list(
  script        = "R/spatial/s06_hippo_regional_atlas.R",
  plan_version  = "v2",
  seurat_version = as.character(packageVersion("Seurat")),
  r_version     = R.version.string,
  run_start     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  git_status    = if (length(git_status) == 0) "(clean)" else paste(git_status, collapse = "; ")
)

# -- Validation tracking --
validation_log <- list()

log_check <- function(id, label, status, detail = "") {
  validation_log[[id]] <<- list(id = id, label = label, status = status, detail = detail)
  cat(sprintf("  [%s] %s: %s %s\n", id, status, label, detail))
}

# ============================================================================
# E1: LOAD HIPPO RDS
# ============================================================================

cat("\n--- E1: Load Hippo RDS ---\n")

if (!file.exists(hippo_rds)) {
  log_check("SC-1", "Hippo RDS exists", "STOP", hippo_rds)
  stop("[SC-1] Hippo RDS file not found: ", hippo_rds)
}

gc()
hippo <- tryCatch({
  readRDS(hippo_rds)
}, error = function(e) {
  cat("[SC-1 STOP] RDS load failed:", conditionMessage(e), "\n")
  stop("[SC-1] RDS load failed")
})
gc()

obj_size <- as.numeric(object.size(hippo))
cat("Object size:", format(object.size(hippo), units = "MB"), "\n")
cat("Dimensions:", paste(dim(hippo), collapse = " x "), "\n")
cat("Assays:", paste(Assays(hippo), collapse = ", "), "\n")

# Record image names
img_names <- tryCatch({
  names(hippo@images)
}, error = function(e) {
  cat("[WARNING] Could not read image names:", conditionMessage(e), "\n")
  character(0)
})
cat("Images:", length(img_names), "—", paste(head(img_names, 4), collapse = ", "), if (length(img_names) > 4) "..." else "", "\n")

# Record reductions
red_names <- names(hippo@reductions)
cat("Reductions:", paste(red_names, collapse = ", "), "\n\n")

# ============================================================================
# E2: RUNTIME INVARIANT VERIFICATION
# ============================================================================

cat("--- E2: Runtime Invariant Verification ---\n")

# V01: Spot count
n_spots <- ncol(hippo)
if (n_spots != expected_total) {
  log_check("V01", "Spot count == 9921", "STOP", paste("got", n_spots))
  stop("[SC-2] Spot count mismatch: expected 9921, got ", n_spots)
} else {
  log_check("V01", "Spot count == 9921", "PASS", paste("=", n_spots))
}

# V02: SCT feature count
sct_features <- tryCatch(nrow(hippo[["SCT"]]), error = function(e) NA)
if (is.na(sct_features) || sct_features != 17096) {
  log_check("V02", "SCT features == 17096", "STOP", paste("got", sct_features))
  stop("[SC-2] SCT feature count mismatch")
} else {
  log_check("V02", "SCT features == 17096", "PASS", paste("=", sct_features))
}

# V03: Spatial feature count
sp_features <- tryCatch(nrow(hippo[["Spatial"]]), error = function(e) NA)
if (is.na(sp_features) || sp_features != 32285) {
  log_check("V03", "Spatial features == 32285", "STOP", paste("got", sp_features))
  stop("[SC-2] Spatial feature count mismatch")
} else {
  log_check("V03", "Spatial features == 32285", "PASS", paste("=", sp_features))
}

# V04: Image count
n_images <- length(img_names)
if (n_images != 16) {
  log_check("V04", "Image count == 16", "STOP", paste("got", n_images))
  stop("[SC-2] Image count mismatch")
} else {
  log_check("V04", "Image count == 16", "PASS", paste("=", n_images))
}

# V05: Region column exists
if (!"Region" %in% colnames(hippo@meta.data)) {
  log_check("V05", "Region column exists", "STOP", "not found")
  stop("[SC-3] Region column not found in metadata")
} else {
  log_check("V05", "Region column exists", "PASS")
}

# V06: Region levels
region_vals <- as.character(hippo$Region)
region_levels <- sort(unique(region_vals))
expected_levels <- sort(c("CA1", "CA2", "CA3", "ML", "GCL", "Hilus"))
if (!identical(region_levels, expected_levels)) {
  log_check("V06", "Region levels == 6 expected", "STOP",
            paste("got", length(region_levels), ":", paste(region_levels, collapse = ", ")))
  stop("[SC-4] Unexpected Region levels: ", paste(region_levels, collapse = ", "))
} else {
  log_check("V06", "Region levels == 6 expected", "PASS", paste(region_levels, collapse = ", "))
}

# V07: SCT assay present
assay_names <- tryCatch(Assays(hippo), error = function(e) names(hippo@assays))
if (!"SCT" %in% assay_names) {
  log_check("V07", "SCT assay present", "STOP")
  stop("[SC-3] SCT assay not found")
} else {
  log_check("V07", "SCT assay present", "PASS")
}

# V08: Spatial assay present
if (!"Spatial" %in% assay_names) {
  log_check("V08", "Spatial assay present", "STOP")
  stop("[SC-3] Spatial assay not found")
} else {
  log_check("V08", "Spatial assay present", "PASS")
}

# V09: Age column exists
if (!"Age" %in% colnames(hippo@meta.data)) {
  log_check("V09", "Age column exists", "STOP", "not found")
  stop("[SC-3] Age column not found")
} else {
  log_check("V09", "Age column exists", "PASS")
}

# V18: Reductions available
missing_reds <- setdiff(c("pca", "umap", "tsne"), red_names)
if (length(missing_reds) > 0) {
  log_check("V18", "Reductions (pca, umap, tsne)", "WARNING",
            paste("missing:", paste(missing_reds, collapse = ", ")))
} else {
  log_check("V18", "Reductions (pca, umap, tsne)", "PASS")
}

cat("\n")

# ============================================================================
# E3: CREATE METADATA_DF AND VALIDATE REGION COUNTS
# ============================================================================

cat("--- E3: Create metadata_df and validate Region counts ---\n")

# Build independent metadata_df (region_group goes here, NOT on Seurat object)
coords <- tryCatch({
  GetTissueCoordinates(hippo, image = img_names[1])
}, error = function(e) {
  cat("[WARNING] GetTissueCoordinates failed for first image:", conditionMessage(e), "\n")
  NULL
})

# Build metadata_df with coordinates from all images
coord_list <- list()
for (img in img_names) {
  cc <- tryCatch({
    GetTissueCoordinates(hippo, image = img)
  }, error = function(e) {
    NULL
  })
  if (!is.null(cc)) {
    cc$image <- img
    coord_list[[img]] <- cc
  }
}

if (length(coord_list) == 0) {
  log_check("V19", "GetTissueCoordinates >= 1 image", "STOP", "all failed")
  stop("[SC-6] All image coordinate fetches failed")
} else {
  log_check("V19", "GetTissueCoordinates >= 1 image", "PASS",
            paste(length(coord_list), "of", length(img_names), "succeeded"))
}

# Build metadata_df from Seurat metadata + coordinates
metadata_df <- hippo@meta.data
metadata_df$cell_id <- rownames(metadata_df)

# Add region_group (CA1/CA2/CA3 preserved, ML/GCL/Hilus -> DG)
metadata_df$region_group <- ifelse(
  metadata_df$Region %in% c("ML", "GCL", "Hilus"), "DG", as.character(metadata_df$Region)
)

# Validate region counts vs Phase 02
region_counts <- table(as.character(metadata_df$Region))
cat("\nRegion counts (Phase 06):\n")
print(region_counts)

for (reg in names(expected_counts)) {
  actual <- as.integer(region_counts[reg])
  expected <- expected_counts[[reg]]
  pct_diff <- abs(actual - expected) / expected * 100
  if (pct_diff > 5) {
    log_check(paste0("V1", match(reg, names(expected_counts)) + 9),
              paste(reg, "count within 5%"), "STOP",
              paste("expected", expected, "got", actual, sprintf("(%.1f%% off)", pct_diff)))
    stop("[SC-5] Region count mismatch: ", reg, " expected ", expected, " got ", actual)
  } else {
    log_check(paste0("V1", match(reg, names(expected_counts)) + 9),
              paste(reg, "count within 5%"), "PASS",
              paste("expected", expected, "got", actual, sprintf("(%.1f%% off)", pct_diff)))
  }
}

# V16: Total
total_actual <- sum(as.integer(region_counts))
if (total_actual != expected_total) {
  log_check("V16", "Total == 9921", "STOP", paste("got", total_actual))
  stop("[SC-5] Total count mismatch")
} else {
  log_check("V16", "Total == 9921", "PASS", paste("=", total_actual))
}

# V17: Derived DG
dg_count <- sum(as.integer(region_counts[c("ML", "GCL", "Hilus")]))
if (dg_count != expected_derived_dg) {
  log_check("V17", "Derived DG (ML+GCL+Hilus) == 2364", "STOP", paste("got", dg_count))
  stop("[SC-5] Derived DG count mismatch")
} else {
  log_check("V17", "Derived DG (ML+GCL+Hilus) == 2364", "PASS", paste("=", dg_count))
}

cat("\n")

# ============================================================================
# E4: EXPORT TABULAR CSVS
# ============================================================================

cat("--- E4: Export tabular CSVs ---\n")

# 4.1 hippo_region_counts.csv
region_count_df <- data.frame(
  Region = names(region_counts),
  Count = as.integer(region_counts),
  Proportion = as.integer(region_counts) / sum(as.integer(region_counts)),
  stringsAsFactors = FALSE
)
write.csv(region_count_df, file.path(out_data, "hippo_region_counts.csv"), row.names = FALSE)
cat("  Wrote hippo_region_counts.csv\n")

# 4.2 hippo_region_age_gemgroup_counts.csv
rag_counts <- metadata_df %>%
  group_by(Region, Age, GEMgroup) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(Region, Age, GEMgroup)
write.csv(rag_counts, file.path(out_data, "hippo_region_age_gemgroup_counts.csv"), row.names = FALSE)
cat("  Wrote hippo_region_age_gemgroup_counts.csv\n")

# 4.3 hippo_region_image_counts.csv
# Map cells to images via coordinates (join by cell_id)
# First build image mapping from coord_list
image_map <- do.call(rbind, lapply(names(coord_list), function(img) {
  cc <- coord_list[[img]]
  # GetTissueCoordinates returns a data.frame; first col is usually cell barcode
  data.frame(cell_id = rownames(cc), image = img, stringsAsFactors = FALSE)
}))
if (is.null(image_map) || nrow(image_map) == 0) {
  cat("[WARNING] Could not build image map; hippo_region_image_counts.csv will be approximate\n")
  image_map <- data.frame(cell_id = metadata_df$cell_id, image = "unknown", stringsAsFactors = FALSE)
}

# Join with metadata
meta_with_img <- metadata_df %>%
  left_join(image_map, by = "cell_id")
ri_counts <- meta_with_img %>%
  group_by(Region, image) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(Region, image)
write.csv(ri_counts, file.path(out_data, "hippo_region_image_counts.csv"), row.names = FALSE)
cat("  Wrote hippo_region_image_counts.csv\n")

# 4.4 hippo_region_age_crosstab.csv
ra_crosstab <- metadata_df %>%
  group_by(Region, Age) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Age, values_from = Count, values_fill = 0)
write.csv(ra_crosstab, file.path(out_data, "hippo_region_age_crosstab.csv"), row.names = FALSE)
cat("  Wrote hippo_region_age_crosstab.csv\n")

# 4.5 hippo_region_coordinate_summary.csv
# Use coordinates from the first available image for summary
coord_summary <- do.call(rbind, lapply(unique(metadata_df$Region), function(reg) {
  sub <- metadata_df[metadata_df$Region == reg, ]
  # Try to get spatial coordinates
  sp_coords <- tryCatch({
    GetTissueCoordinates(hippo, image = img_names[1])
  }, error = function(e) NULL)
  if (!is.null(sp_coords)) {
    common_cells <- intersect(rownames(sp_coords), rownames(sub))
    if (length(common_cells) > 0) {
      cc <- sp_coords[common_cells, ]
      return(data.frame(
        Region = reg,
        n_spots = length(common_cells),
        x_min = min(cc[[1]], na.rm = TRUE),
        x_max = max(cc[[1]], na.rm = TRUE),
        y_min = min(cc[[2]], na.rm = TRUE),
        y_max = max(cc[[2]], na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
  data.frame(Region = reg, n_spots = nrow(sub), x_min = NA, x_max = NA, y_min = NA, y_max = NA,
             stringsAsFactors = FALSE)
}))
write.csv(coord_summary, file.path(out_data, "hippo_region_coordinate_summary.csv"), row.names = FALSE)
cat("  Wrote hippo_region_coordinate_summary.csv\n")

# 4.6 hippo_derived_dg_validation.csv
# Compare derived DG (ML+GCL+Hilus from Hippo) vs cross-object DG
cross_obj <- tryCatch(read.csv(p02_cross, stringsAsFactors = FALSE), error = function(e) NULL)
dg_cross <- if (!is.null(cross_obj) && "DG_total_spots" %in% cross_obj$check) {
  as.integer(cross_obj$value[cross_obj$check == "DG_total_spots"])
} else {
  NA_integer_
}
derived_dg_df <- data.frame(
  Metric = c("Hippo_ML", "Hippo_GCL", "Hippo_Hilus", "Derived_DG_sum", "Cross_object_DG", "Match"),
  Value = c(
    as.character(as.integer(region_counts["ML"])),
    as.character(as.integer(region_counts["GCL"])),
    as.character(as.integer(region_counts["Hilus"])),
    as.character(dg_count),
    if (is.na(dg_cross)) "NA" else as.character(dg_cross),
    if (is.na(dg_cross)) "NA" else as.character(dg_count == dg_cross)
  ),
  stringsAsFactors = FALSE
)
write.csv(derived_dg_df, file.path(out_data, "hippo_derived_dg_validation.csv"), row.names = FALSE)
cat("  Wrote hippo_derived_dg_validation.csv\n")

# 4.7 hippo_region_metadata.csv
# Full metadata with region_group (independent df, NOT written to Seurat object)
write.csv(metadata_df, file.path(out_data, "hippo_region_metadata.csv"), row.names = FALSE)
cat("  Wrote hippo_region_metadata.csv\n")

# 4.8 phase06_provenance.csv
provenance$run_end <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
provenance$hippo_spots <- n_spots
provenance$hippo_sct_features <- sct_features
provenance$hippo_spatial_features <- sp_features
provenance$hippo_images <- n_images
provenance$object_size_mb <- round(obj_size / 1e6, 1)
prov_df <- data.frame(
  Key = names(provenance),
  Value = as.character(provenance),
  stringsAsFactors = FALSE
)
write.csv(prov_df, file.path(out_data, "phase06_provenance.csv"), row.names = FALSE)
cat("  Wrote phase06_provenance.csv\n\n")

# ============================================================================
# E5: PER-IMAGE SPATIAL REGION MAPS (16 plots)
# ============================================================================

cat("--- E5: Per-image spatial region maps ---\n")

# Function to create spatial plot for one image
make_spatial_plot <- function(obj, image_name, region_vec, color_map) {
  coords <- tryCatch({
    GetTissueCoordinates(obj, image = image_name)
  }, error = function(e) NULL)

  if (is.null(coords)) {
    cat("  [FAIL]", image_name, "- GetTissueCoordinates failed\n")
    return(NULL)
  }

  # coords has cell barcodes as rownames; first two cols are spatial coords
  coord_df <- data.frame(
    x = coords[[1]],
    y = coords[[2]],
    Region = region_vec[rownames(coords)],
    stringsAsFactors = FALSE
  )
  coord_df <- coord_df[!is.na(coord_df$Region), ]

  p <- ggplot(coord_df, aes(x = x, y = y, color = Region)) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_color_manual(values = color_map) +
    coord_fixed() +
    labs(title = image_name, x = "Spatial_1", y = "Spatial_2") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )

  return(p)
}

region_vec <- setNames(as.character(hippo$Region), rownames(hippo@meta.data))
img_plot_status <- data.frame(image = img_names, status = "PASS", stringsAsFactors = FALSE)

for (i in seq_along(img_names)) {
  img <- img_names[i]
  p <- make_spatial_plot(hippo, img, region_vec, region_colors)

  if (is.null(p)) {
    img_plot_status$status[i] <- "FAIL"
    next
  }

  fname <- file.path(out_figs, paste0("phase06_spatial_region_", img, ".png"))
  ggsave(fname, plot = p, width = 6, height = 5, dpi = 150, bg = "white")
  cat("  Saved:", basename(fname), "\n")
}

n_img_passed <- sum(img_plot_status$status == "PASS")
cat(sprintf("\n  Image plots: %d/%d passed\n\n", n_img_passed, nrow(img_plot_status)))

# ============================================================================
# E6: OVERVIEW GRID (4x4)
# ============================================================================

cat("--- E6: Overview 4x4 grid ---\n")

# Collect all successful image plots
all_img_plots <- list()
for (i in seq_along(img_names)) {
  img <- img_names[i]
  p <- make_spatial_plot(hippo, img, region_vec, region_colors)
  if (!is.null(p)) {
    all_img_plots[[img]] <- p + theme(legend.position = "none")
  }
}

if (length(all_img_plots) > 0) {
  # Use patchwork to create grid
  n_plots <- length(all_img_plots)
  nrow_grid <- ceiling(n_plots / 4)
  overview <- wrap_plots(all_img_plots, ncol = 4) +
    plot_annotation(
      title = "Hippocampus - All Slices - Region Labels",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  fname_overview <- file.path(out_figs, "phase06_spatial_region_all_slices.png")
  ggsave(fname_overview, plot = overview, width = 20, height = 5 * nrow_grid, dpi = 150, bg = "white")
  cat("  Saved: phase06_spatial_region_all_slices.png\n\n")
} else {
  cat("  [FAIL] No image plots available for overview grid\n\n")
}

gc()

# ============================================================================
# E7: UMAP PLOTS (Region, Age, GEMgroup)
# ============================================================================

cat("--- E7: UMAP plots ---\n")

if ("umap" %in% red_names) {
  umap_coords <- as.data.frame(Embeddings(hippo, "umap"))
  umap_coords$Region   <- as.character(hippo$Region)
  umap_coords$Age      <- as.character(hippo$Age)
  umap_coords$GEMgroup <- as.character(hippo$GEMgroup)

  for (var in c("Region", "Age", "GEMgroup")) {
    color_map <- if (var == "Region") region_colors else NULL
    p <- ggplot(umap_coords, aes_string(x = "UMAP_1", y = "UMAP_2", color = var)) +
      geom_point(size = 0.8, alpha = 0.7) +
      {if (!is.null(color_map)) scale_color_manual(values = color_map) else theme()} +
      labs(title = paste("UMAP by", var)) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"))

    fname <- file.path(out_figs, paste0("phase06_umap_by_", tolower(var), ".png"))
    ggsave(fname, plot = p, width = 7, height = 6, dpi = 150, bg = "white")
    cat("  Saved:", basename(fname), "\n")
  }
} else {
  cat("  [WARNING] UMAP reduction not available; skipping UMAP plots\n")
}

gc()

# ============================================================================
# E8: tSNE PLOTS (Region, Age, GEMgroup)
# ============================================================================

cat("\n--- E8: tSNE plots ---\n")

if ("tsne" %in% red_names) {
  tsne_coords <- as.data.frame(Embeddings(hippo, "tsne"))
  tsne_coords$Region   <- as.character(hippo$Region)
  tsne_coords$Age      <- as.character(hippo$Age)
  tsne_coords$GEMgroup <- as.character(hippo$GEMgroup)

  for (var in c("Region", "Age", "GEMgroup")) {
    color_map <- if (var == "Region") region_colors else NULL
    p <- ggplot(tsne_coords, aes_string(x = "tSNE_1", y = "tSNE_2", color = var)) +
      geom_point(size = 0.8, alpha = 0.7) +
      {if (!is.null(color_map)) scale_color_manual(values = color_map) else theme()} +
      labs(title = paste("tSNE by", var)) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"))

    fname <- file.path(out_figs, paste0("phase06_tsne_by_", tolower(var), ".png"))
    ggsave(fname, plot = p, width = 7, height = 6, dpi = 150, bg = "white")
    cat("  Saved:", basename(fname), "\n")
  }
} else {
  cat("  [WARNING] tSNE reduction not available; skipping tSNE plots\n")
}

gc()

# ============================================================================
# E9: PCA PLOTS (Region, Age, GEMgroup)
# ============================================================================

cat("\n--- E9: PCA plots ---\n")

if ("pca" %in% red_names) {
  pca_coords <- as.data.frame(Embeddings(hippo, "pca")[, 1:2])
  colnames(pca_coords) <- c("PC_1", "PC_2")
  pca_coords$Region   <- as.character(hippo$Region)
  pca_coords$Age      <- as.character(hippo$Age)
  pca_coords$GEMgroup <- as.character(hippo$GEMgroup)

  for (var in c("Region", "Age", "GEMgroup")) {
    color_map <- if (var == "Region") region_colors else NULL
    p <- ggplot(pca_coords, aes_string(x = "PC_1", y = "PC_2", color = var)) +
      geom_point(size = 0.8, alpha = 0.7) +
      {if (!is.null(color_map)) scale_color_manual(values = color_map) else theme()} +
      labs(title = paste("PCA by", var)) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"))

    fname <- file.path(out_figs, paste0("phase06_pca_by_", tolower(var), ".png"))
    ggsave(fname, plot = p, width = 7, height = 6, dpi = 150, bg = "white")
    cat("  Saved:", basename(fname), "\n")
  }
} else {
  cat("  [WARNING] PCA reduction not available; skipping PCA plots\n")
}

gc()

# ============================================================================
# E10: VALIDATION SUMMARY
# ============================================================================

cat("\n--- E10: Validation summary ---\n")

# V20+: Check required output files exist (CSV/TXT, then PNG, then validation summary)
required_csvs <- c(
  "hippo_region_counts.csv",
  "hippo_region_age_gemgroup_counts.csv",
  "hippo_region_image_counts.csv",
  "hippo_region_age_crosstab.csv",
  "hippo_region_coordinate_summary.csv",
  "hippo_derived_dg_validation.csv",
  "hippo_region_metadata.csv",
  "phase06_provenance.csv"
)

required_pngs <- c(
  paste0("phase06_spatial_region_", img_names, ".png"),
  "phase06_spatial_region_all_slices.png",
  paste0("phase06_", rep(c("umap", "tsne", "pca"), each = 3),
         "_by_", rep(c("region", "age", "gemgroup"), 3), ".png")
)

v_idx <- 20
for (f in required_csvs) {
  exists_f <- file.exists(file.path(out_data, f))
  log_check(paste0("V", v_idx), paste("CSV exists:", f),
            if (exists_f) "PASS" else "FAIL")
  v_idx <- v_idx + 1
}

for (f in required_pngs) {
  exists_f <- file.exists(file.path(out_figs, f))
  log_check(paste0("V", v_idx), paste("PNG exists:", f),
            if (exists_f) "PASS" else "FAIL")
  v_idx <- v_idx + 1
}

# Overall status
all_statuses <- sapply(validation_log, function(x) x$status)
n_pass   <- sum(all_statuses == "PASS")
n_fail   <- sum(all_statuses == "FAIL")
n_warn   <- sum(all_statuses == "WARNING")
n_stop   <- sum(all_statuses == "STOP")

# Determine overall: any STOP -> STOPPED; any FAIL -> FAIL; else PASS
overall <- if (n_stop > 0) "STOPPED" else if (n_fail > 0) "FAIL" else "PASS"

cat("\n")
cat("========================================\n")
cat("  Phase Spatial-06 Validation Summary\n")
cat("========================================\n")
cat(sprintf("  PASS: %d | FAIL: %d | WARNING: %d | STOP: %d\n", n_pass, n_fail, n_warn, n_stop))
cat(sprintf("  Overall: %s\n", overall))
cat("========================================\n\n")

# Write validation summary
val_lines <- c(
  "=== Phase Spatial-06 Validation Summary ===",
  paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Overall status:", overall),
  "",
  "--- Check Results ---"
)
for (vl in validation_log) {
  val_lines <- c(val_lines, sprintf("[%s] %s: %s %s", vl$id, vl$status, vl$label, vl$detail))
}
val_lines <- c(val_lines, "", sprintf("PASS: %d | FAIL: %d | WARNING: %d | STOP: %d", n_pass, n_fail, n_warn, n_stop))

writeLines(val_lines, file.path(out_data, "phase06_validation_summary.txt"))
cat("  Wrote phase06_validation_summary.txt\n")

# Check validation summary file exists (after writing) — uses next sequential v_idx
val_file_exists <- file.exists(file.path(out_data, "phase06_validation_summary.txt"))
log_check(paste0("V", v_idx), "TXT exists: phase06_validation_summary.txt",
          if (val_file_exists) "PASS" else "FAIL")
v_idx <- v_idx + 1

# Recompute overall after all checks including V28
all_statuses <- sapply(validation_log, function(x) x$status)
n_pass   <- sum(all_statuses == "PASS")
n_fail   <- sum(all_statuses == "FAIL")
n_warn   <- sum(all_statuses == "WARNING")
n_stop   <- sum(all_statuses == "STOP")
overall <- if (n_stop > 0) "STOPPED" else if (n_fail > 0) "FAIL" else "PASS"

# Rewrite validation summary with final status
val_lines_final <- c(
  "=== Phase Spatial-06 Validation Summary ===",
  paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Overall status:", overall),
  "",
  "--- Check Results ---"
)
for (vl in validation_log) {
  val_lines_final <- c(val_lines_final, sprintf("[%s] %s: %s %s", vl$id, vl$status, vl$label, vl$detail))
}
val_lines_final <- c(val_lines_final, "", sprintf("PASS: %d | FAIL: %d | WARNING: %d | STOP: %d", n_pass, n_fail, n_warn, n_stop))
writeLines(val_lines_final, file.path(out_data, "phase06_validation_summary.txt"))

# ============================================================================
# E11: Rplots.pdf GUARD + CLEANUP
# ============================================================================

cat("\n--- E11: Cleanup ---\n")

# Rplots.pdf guard (end)
rplots_end <- file.exists("Rplots.pdf")
if (rplots_end && !rplots_start) {
  cat("[SC-9 WARNING] Rplots.pdf was created during script execution. Removing.\n")
  file.remove("Rplots.pdf")
} else if (rplots_end && rplots_start) {
  cat("[SC-9 WARNING] Rplots.pdf existed before and after execution. Not removing.\n")
}

# Final gc
gc()

cat("\n=== Phase Spatial-06 complete ===\n")
cat("Overall status:", overall, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

# Exit with appropriate code
if (overall == "STOPPED") {
  quit(status = 1)
} else if (overall == "FAIL") {
  quit(status = 2)
}
