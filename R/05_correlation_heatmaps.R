## =========================================================
## Script: 05_correlation_heatmaps.R
## Purpose: Spearman correlation heatmaps (Fig 6 & 7)
## Input: data/processed/neuro_with_ferroptosis_scores.rds
## Output: figures/stage-2/Fig6*.png + .pdf, Fig7*.png + .pdf
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(grid)
})

message("=== Step 05: Correlation Heatmaps ===")
message("Start time: ", Sys.time())

## =========================================================
## 0. Paths
## =========================================================
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "figures", "stage-2")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_pheatmap <- function(ph, filename, width = 10, height = 8) {
  out_png <- file.path(fig_dir, paste0(filename, ".png"))
  out_pdf <- file.path(fig_dir, paste0(filename, ".pdf"))
  
  png(out_png, width = width, height = height, units = "in", res = 300)
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
  
  tryCatch({
    cairo_pdf(out_pdf, width = width, height = height)
    grid.newpage()
    grid.draw(ph$gtable)
    dev.off()
  }, error = function(e) {
    warning("cairo_pdf failed, using pdf(): ", e$message)
    pdf(out_pdf, width = width, height = height)
    grid.newpage()
    grid.draw(ph$gtable)
    dev.off()
  })
  
  message("  [Saved] ", filename)
}

## =========================================================
## 1. Load data
## =========================================================
message("[1/6] Loading Seurat object...")
neuro <- readRDS(file.path(data_dir, "neuro_with_ferroptosis_scores.rds"))
DefaultAssay(neuro) <- "RNA"
message("  Loaded: ", ncol(neuro), " cells")

## =========================================================
## 2. Define gene sets
## =========================================================
message("[2/6] Defining gene sets...")

# Neurogenesis markers
neuro_markers_raw <- c("Sox2", "Hopx", "Ascl1", "Dcx", "Neurod1", "Prox1", "Calb2", "Nestin", "Nes")
neuro_markers <- neuro_markers_raw[neuro_markers_raw %in% rownames(neuro)]
message("  Neurogenesis markers found: ", paste(neuro_markers, collapse = ", "))

# Ferroptosis genes (Human → Mouse conversion)
human_to_mouse <- function(genes) {
  gl <- tolower(genes)
  paste0(toupper(substr(gl, 1, 1)), substr(gl, 2, nchar(gl)))
}

ferro_promoter_raw <- c("TFRC","NCOA4","HMOX1","PHKG2","ACO1","IREB2",
                         "ACSL4","LPCAT3","ALOX5","ALOX12","ALOX15",
                         "DPP4","ACSF2","ZEB1","PTGS2","PEBP1",
                         "CARS","CHAC1","NOX1","ABCC1","NOX4",
                         "GLS2","ATP5G3","VDAC3","SAT1")

ferro_inhibitor_raw <- c("FTH1","NFS1","MT1G","STEAP3","HSBP1","FANCD2","TERF1",
                          "GPX4","AKR1C1","AKR1C2","AKR1C3","HMGCR","SQLE","CISD1","ACSL3",
                          "CS","CBS","FDFT1","FADS2","ACACA",
                          "SLC7A11","GCLC","GCLM","GSS","NFE2L2","KEAP1","SLC38A1",
                          "G6PD","PGD","SLC1A5","GOT1",
                          "CD44","CRYAB","EMC2","HSPB1","AIFM2","RPL8")

ferro_promoter  <- human_to_mouse(ferro_promoter_raw)
ferro_inhibitor <- human_to_mouse(ferro_inhibitor_raw)
ferro_all       <- c(ferro_promoter, ferro_inhibitor)
ferro_all       <- ferro_all[ferro_all %in% rownames(neuro)]
message("  Ferroptosis genes found: ", length(ferro_all))

## =========================================================
## 3. Compute correlation matrix (helper)
## =========================================================
compute_corr_matrix <- function(neuro_obj, genes_row, genes_col, subset_expr = NULL) {
  # Get expression data
  all_genes <- unique(c(genes_row, genes_col))
  all_genes <- all_genes[all_genes %in% rownames(neuro_obj)]
  
  expr <- as.matrix(FetchData(neuro_obj, vars = all_genes))
  
  # Subset cells if needed
  if (!is.null(subset_expr)) {
    cells_use <- eval(parse(text = subset_expr), envir = neuro_obj@meta.data)
    expr <- expr[cells_use, , drop = FALSE]
  }
  
  # Compute Spearman correlation
  genes_row <- genes_row[genes_row %in% colnames(expr)]
  genes_col <- genes_col[genes_col %in% colnames(expr)]
  
  cor_mat <- cor(expr[, genes_row, drop = FALSE],
                 expr[, genes_col, drop = FALSE],
                 method = "spearman", use = "pairwise.complete.obs")
  
  # Replace any remaining NA/NaN/Inf with 0
  cor_mat[is.na(cor_mat) | is.infinite(cor_mat)] <- 0
  
  return(cor_mat)
}

## =========================================================
## 4. Figure 6: Correlation heatmaps
## =========================================================
message("[3/6] Figure 6: Correlation heatmaps...")

# Fig 6A: All cells
message("  Fig 6A: All cells...")
cor_all <- compute_corr_matrix(neuro, ferro_all, neuro_markers)

ph6a <- pheatmap(cor_all, main = "Figure 6A. Correlation: All cells",
                  color = colorRampPalette(c("#3182BD", "white", "#D73027"))(100),
                  fontsize = 8, cluster_rows = TRUE, cluster_cols = TRUE)
save_pheatmap(ph6a, "Fig6A_Correlation_Heatmap_All")
gc()

# Fig 6B: Young only
message("  Fig 6B: Young only...")
cor_young <- compute_corr_matrix(neuro, ferro_all, neuro_markers,
                                  subset_expr = "Age_collapsed == 'Young'")

ph6b <- pheatmap(cor_young, main = "Figure 6B. Correlation: Young",
                  color = colorRampPalette(c("#3182BD", "white", "#D73027"))(100),
                  fontsize = 8, cluster_rows = TRUE, cluster_cols = TRUE)
save_pheatmap(ph6b, "Fig6B_Correlation_Heatmap_Young")
gc()

# Fig 6C: Old only
message("  Fig 6C: Old only...")
cor_old <- compute_corr_matrix(neuro, ferro_all, neuro_markers,
                                subset_expr = "Age_collapsed == 'Old'")

ph6c <- pheatmap(cor_old, main = "Figure 6C. Correlation: Old",
                  color = colorRampPalette(c("#3182BD", "white", "#D73027"))(100),
                  fontsize = 8, cluster_rows = TRUE, cluster_cols = TRUE)
save_pheatmap(ph6c, "Fig6C_Correlation_Heatmap_Old")
gc()

## =========================================================
## 5. Figure 7: Fixed-order heatmaps (unified color scale)
## =========================================================
message("[4/6] Figure 7: Fixed-order heatmaps...")

# Fixed gene order
ferro_promoter_use  <- ferro_promoter[ferro_promoter %in% rownames(neuro)]
ferro_inhibitor_use <- ferro_inhibitor[ferro_inhibitor %in% rownames(neuro)]
ferro_fixed <- c(ferro_promoter_use, ferro_inhibitor_use)
neuro_fixed <- neuro_markers

# Unified color scale
breaks_seq <- seq(-0.5, 0.5, length.out = 101)

# Fig 7A: All cells
message("  Fig 7A: Fixed-order All...")
cor_fixed_all <- compute_corr_matrix(neuro, ferro_fixed, neuro_fixed)
# Reorder
cor_fixed_all <- cor_fixed_all[ferro_fixed[ferro_fixed %in% rownames(cor_fixed_all)],
                                neuro_fixed[neuro_fixed %in% colnames(cor_fixed_all)], drop = FALSE]

ph7a <- pheatmap(cor_fixed_all, main = "Figure 7A. Fixed-order: All cells",
                  color = colorRampPalette(c("#3182BD", "white", "#D73027"))(100),
                  breaks = breaks_seq,
                  fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE)
save_pheatmap(ph7a, "Fig7A_FixedOrder_Heatmap_All")
gc()

# Fig 7B: Young
message("  Fig 7B: Fixed-order Young...")
cor_fixed_young <- compute_corr_matrix(neuro, ferro_fixed, neuro_fixed,
                                        subset_expr = "Age_collapsed == 'Young'")
cor_fixed_young <- cor_fixed_young[ferro_fixed[ferro_fixed %in% rownames(cor_fixed_young)],
                                    neuro_fixed[neuro_fixed %in% colnames(cor_fixed_young)], drop = FALSE]

ph7b <- pheatmap(cor_fixed_young, main = "Figure 7B. Fixed-order: Young",
                  color = colorRampPalette(c("#3182BD", "white", "#D73027"))(100),
                  breaks = breaks_seq,
                  fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE)
save_pheatmap(ph7b, "Fig7B_FixedOrder_Heatmap_Young")
gc()

# Fig 7C: Old
message("  Fig 7C: Fixed-order Old...")
cor_fixed_old <- compute_corr_matrix(neuro, ferro_fixed, neuro_fixed,
                                      subset_expr = "Age_collapsed == 'Old'")
cor_fixed_old <- cor_fixed_old[ferro_fixed[ferro_fixed %in% rownames(cor_fixed_old)],
                                neuro_fixed[neuro_fixed %in% colnames(cor_fixed_old)], drop = FALSE]

ph7c <- pheatmap(cor_fixed_old, main = "Figure 7C. Fixed-order: Old",
                  color = colorRampPalette(c("#3182BD", "white", "#D73027"))(100),
                  breaks = breaks_seq,
                  fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE)
save_pheatmap(ph7c, "Fig7C_FixedOrder_Heatmap_Old")
gc()

## =========================================================
## 6. Summary
## =========================================================
message("\n[5/6] Listing generated figures...")
fig_files <- list.files(fig_dir, pattern = "^Fig[67]", full.names = FALSE)
message("  Generated ", length(fig_files), " files:")
for (f in fig_files) message("    ", f)

message("\n=== Step 05 Complete ===")
message("End time: ", Sys.time())
message("Output directory: ", fig_dir)
