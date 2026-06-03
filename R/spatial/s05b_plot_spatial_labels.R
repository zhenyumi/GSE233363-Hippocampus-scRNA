#!/usr/bin/env Rscript
# ---- s05b_plot_spatial_labels.R ----
# Plot spatial label distribution on Visium slices.
# Uses Phase 05 outputs (phase05_merged_spots.csv) + RDS coordinates.
# No DESeq2, no neighbor computation, no module score.
# RDS-based approximate spatial label map; not strict STutility reproduction.
#
# ---- Rplots.pdf guard ----
# Records root-level Rplots.pdf state at start; cleans up at end.
# Does NOT use cairo_pdf; only PNG devices.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(42)

# Null-coalescing operator (base R >= 4.4 does not guarantee `%||%`)
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- 0. Record Rplots.pdf pre-existing state ----
rplots_path <- file.path(getwd(), "Rplots.pdf")
rplots_pre <- list(
  exists   = file.exists(rplots_path),
  size     = NA_real_,
  mtime    = NA_character_
)
if (rplots_pre$exists) {
  fi <- file.info(rplots_path)
  rplots_pre$size  <- fi$size
  rplots_pre$mtime <- as.character(fi$mtime)
  cat(sprintf("[s05b] Root-level Rplots.pdf EXISTS at start: %d bytes, mtime=%s\n",
    rplots_pre$size, rplots_pre$mtime))
} else {
  cat("[s05b] No root-level Rplots.pdf at start.\n")
}

# ---- 1. Paths ----
data_dir   <- "data/processed/spatial/phase05_rds_gradient"
fig_dir    <- file.path("figures/spatial/phase05/spatial_label_maps")
prov_path  <- file.path(data_dir, "phase05_provenance.csv")
val_path   <- file.path(data_dir, "phase05_validation_summary.txt")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

merged_path <- file.path(data_dir, "phase05_merged_spots.csv")
rds_path    <- "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"

stopifnot(file.exists(merged_path))
stopifnot(file.exists(rds_path))

# ---- 2. Load data ----
cat("[s05b] Loading merged spots...\n")
merged <- read.csv(merged_path, stringsAsFactors = FALSE)
cat(sprintf("[s05b]   %d rows, columns: %s\n", nrow(merged), paste(colnames(merged), collapse=", ")))

cat("[s05b] Loading Hippo RDS (coordinates only)...\n")
hippo <- readRDS(rds_path)
cat(sprintf("[s05b]   %d cells, %d images\n", ncol(hippo), length(Images(hippo))))

stopifnot(nrow(merged) == ncol(hippo))
stopifnot(all(c("barcode", "image", "FuncRegion") %in% colnames(merged)))

# ---- 3. Validate image set ----
expected_images <- c("slice1", paste0("slice1_", 2:16))
actual_images   <- Images(hippo)
stopifnot(setequal(actual_images, expected_images))
cat(sprintf("[s05b] %d images confirmed.\n", length(actual_images)))

# ---- 4. Color scheme ----
region_colors <- c(
  "IS"        = "#E74C3C",   # red
  "NNS"       = "#3498DB",   # blue
  "ENS1"      = "#F39C12",   # orange
  "ENS2"      = "#27AE60",   # green
  "unlabelled" = "#BDC3C7"   # grey
)

# ---- 5. Per-image coordinate extraction and plotting ----
cat("[s05b] Extracting coordinates and plotting per-image spatial label maps...\n")

plot_list     <- list()
manifest_rows <- list()

for (img in actual_images) {
  cat(sprintf("[s05b]   %s ...", img))

  # Extract coordinates from Seurat
  coords <- tryCatch(
    Seurat::GetTissueCoordinates(hippo, image = img),
    error = function(e) NULL
  )

  if (is.null(coords) || nrow(coords) == 0) {
    cat(" FAIL (no coordinates)\n")
    manifest_rows[[img]] <- data.frame(
      image = img, n_spots = 0L, n_IS = 0L, n_NNS = 0L,
      n_ENS1 = 0L, n_ENS2 = 0L, n_Unlabelled = 0L,
      png_path = NA_character_, plot_status = "FAIL_no_coordinates",
      stringsAsFactors = FALSE
    )
    next
  }

  coords_df <- data.frame(
    barcode   = rownames(coords),
    imagecol  = coords$imagecol,
    imagerow  = coords$imagerow,
    stringsAsFactors = FALSE
  )

  # Merge with label data
  img_labels <- merged[merged$image == img, c("barcode", "FuncRegion")]
  coords_df <- merge(coords_df, img_labels, by = "barcode", all.x = TRUE)
  coords_df$FuncRegion[is.na(coords_df$FuncRegion)] <- "unlabelled"

  # Count labels
  tab <- table(coords_df$FuncRegion)

  # Plot
  p <- ggplot(coords_df, aes(x = imagecol, y = imagerow, color = FuncRegion)) +
    geom_point(size = 1, alpha = 0.8) +
    scale_color_manual(values = region_colors, drop = FALSE) +
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title = paste("Spatial labels —", img),
      subtitle = "RDS-based approximate spatial label map; not strict STutility reproduction.",
      x = "imagecol", y = "imagerow", color = "FuncRegion"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8, face = "italic", colour = "grey40")
    )

  png_name <- sprintf("phase05_spatial_labels_%s.png", img)
  png_path <- file.path(fig_dir, png_name)

  tryCatch({
    ggsave(png_path, p, width = 6, height = 5, dpi = 150)
    cat(sprintf(" OK (%d IS, %d NNS, %d ENS1, %d ENS2, %d unlab)\n",
      tab["IS"] %||% 0L, tab["NNS"] %||% 0L,
      tab["ENS1"] %||% 0L, tab["ENS2"] %||% 0L,
      tab["unlabelled"] %||% 0L))
    plot_status <- "OK"
  }, error = function(e) {
    cat(sprintf(" FAIL (%s)\n", e$message))
    png_path <<- NA_character_
    plot_status <<- "FAIL_ggsave"
  })

  plot_list[[img]] <- p
  manifest_rows[[img]] <- data.frame(
    image       = img,
    n_spots     = nrow(coords_df),
    n_IS        = as.integer(tab["IS"] %||% 0L),
    n_NNS       = as.integer(tab["NNS"] %||% 0L),
    n_ENS1      = as.integer(tab["ENS1"] %||% 0L),
    n_ENS2      = as.integer(tab["ENS2"] %||% 0L),
    n_Unlabelled = as.integer(tab["unlabelled"] %||% 0L),
    png_path    = png_path,
    plot_status = plot_status,
    stringsAsFactors = FALSE
  )

  rm(coords, coords_df, img_labels, p)
  gc()
}

# ---- 6. Overview figure (4x4 grid) ----
cat("[s05b] Creating overview figure (4x4 grid)...\n")

if (requireNamespace("cowplot", quietly = TRUE)) {
  # Remove legend from individual plots for overview
  plot_list_noleg <- lapply(plot_list, function(p) {
    p + theme(legend.position = "none")
  })

  overview <- cowplot::plot_grid(
    plotlist = plot_list_noleg, ncol = 4, nrow = 4,
    labels = actual_images, label_size = 7, label_x = 0.05, hjust = 0
  )

  # Add shared legend at bottom
  legend_b <- cowplot::get_legend(
    plot_list[[1]] + theme(legend.position = "bottom")
  )
  overview_with_legend <- cowplot::plot_grid(overview, legend_b, ncol = 1, rel_heights = c(1, 0.08))

  overview_path <- file.path(fig_dir, "phase05_spatial_labels_all_slices.png")
  tryCatch({
    ggsave(overview_path, overview_with_legend, width = 16, height = 17, dpi = 150)
    cat(sprintf("[s05b] Overview saved: %s\n", overview_path))
  }, error = function(e) {
    cat(sprintf("[s05b] Overview save FAILED: %s\n", e$message))
    overview_path <<- NA_character_
  })
} else {
  cat("[s05b] WARNING: cowplot not available. Skipping overview figure.\n")
  overview_path <- NA_character_
}

# ---- 7. Save manifest ----
cat("[s05b] Saving spatial plot manifest...\n")
manifest <- do.call(rbind, manifest_rows)
manifest$overview_path <- overview_path

manifest_out_path <- file.path(data_dir, "phase05_spatial_plot_manifest.csv")
write.csv(manifest, manifest_out_path, row.names = FALSE)
cat(sprintf("[s05b] Manifest saved: %s (%d rows)\n", manifest_out_path, nrow(manifest)))

# ---- 8. Update provenance ----
cat("[s05b] Updating provenance...\n")
prov <- read.csv(prov_path, stringsAsFactors = FALSE)

new_prov_rows <- data.frame(
  key   = c("spatial_plot_route", "spatial_label_maps_generated",
            "coordinates_source", "spatial_plot_date"),
  value = c("RDS_coordinates", "TRUE",
            "GetTissueCoordinates", as.character(Sys.time())),
  stringsAsFactors = FALSE
)
# Avoid duplicating keys already present
new_prov_rows <- new_prov_rows[!new_prov_rows$key %in% prov$key, ]

if (nrow(new_prov_rows) > 0) {
  prov <- rbind(prov, new_prov_rows)
  write.csv(prov, prov_path, row.names = FALSE)
  cat(sprintf("[s05b] Provenance updated: %d rows (was %d)\n", nrow(prov), nrow(prov) - nrow(new_prov_rows)))
} else {
  cat("[s05b] Provenance already up-to-date.\n")
}

# ---- 9. Update validation summary ----
cat("[s05b] Updating validation summary...\n")
val_lines <- readLines(val_path, warn = FALSE)

# Check if spatial label section already exists
if (!any(grepl("Spatial Label Maps", val_lines))) {
  spatial_section <- c(
    "",
    "--- Spatial Label Distribution Figures ---",
    sprintf("Generated: %s", as.character(Sys.time())),
    sprintf("Overview: %s", ifelse(is.na(overview_path), "FAILED", overview_path)),
    sprintf("Per-slice PNGs: %d generated", sum(!is.na(manifest$png_path))),
    sprintf("Per-slice figure dir: %s", fig_dir),
    sprintf("Manifest: %s", manifest_out_path),
    "Coordinates source: GetTissueCoordinates (RDS-contained)",
    "Route: RDS-approximated; NOT equivalent to STutility RegionNeighbours",
    "PDF figures: not generated (cairo/X11 unavailable); PNGs are authoritative",
    "Rplots.pdf: cleaned up by s05b if generated during this run"
  )
  val_lines <- c(val_lines, spatial_section)
  writeLines(val_lines, val_path)
  cat("[s05b] Validation summary updated with spatial label section.\n")
} else {
  cat("[s05b] Validation summary already has spatial label section.\n")
}

# ---- 10. Rplots.pdf cleanup ----
cat("[s05b] Rplots.pdf cleanup check...\n")
rplots_post <- list(
  exists = file.exists(rplots_path),
  size   = NA_real_,
  mtime  = NA_character_
)
if (rplots_post$exists) {
  fi <- file.info(rplots_path)
  rplots_post$size  <- fi$size
  rplots_post$mtime <- as.character(fi$mtime)
}

if (rplots_post$exists && !rplots_pre$exists) {
  # Rplots.pdf was created during this run
  file.remove(rplots_path)
  cat(sprintf("[s05b] REMOVED root-level Rplots.pdf (generated during this run, %d bytes)\n",
    rplots_post$size))
} else if (rplots_post$exists && rplots_pre$exists &&
           rplots_post$mtime != rplots_pre$mtime) {
  # Rplots.pdf existed before but was modified during this run
  file.remove(rplots_path)
  cat(sprintf("[s05b] REMOVED root-level Rplots.pdf (mtime changed during this run: %s -> %s)\n",
    rplots_pre$mtime, rplots_post$mtime))
} else if (rplots_post$exists) {
  cat(sprintf("[s05b] Rplots.pdf exists but was NOT modified during this run. NOT deleted. Size: %d bytes, mtime: %s\n",
    rplots_post$size, rplots_post$mtime))
} else {
  cat("[s05b] No root-level Rplots.pdf found. Clean.\n")
}

# ---- 11. Final summary ----
cat("\n=== s05b Complete ===\n")
cat(sprintf("Per-slice PNGs: %d / 16\n", sum(!is.na(manifest$png_path))))
cat(sprintf("Overview PNG: %s\n", ifelse(is.na(overview_path), "NOT generated", overview_path)))
cat(sprintf("Manifest: %s\n", manifest_out_path))
cat(sprintf("Total spots plotted: %d\n", sum(manifest$n_spots)))
cat(sprintf("  IS: %d | NNS: %d | ENS1: %d | ENS2: %d | unlabelled: %d\n",
  sum(manifest$n_IS), sum(manifest$n_NNS),
  sum(manifest$n_ENS1), sum(manifest$n_ENS2), sum(manifest$n_Unlabelled)))
cat("No DESeq2 run. No package install. No renv.lock change.\n")
cat("Done.\n")
