## =========================================================
## Script: 04_de_and_volcano.R
## Purpose: Differential expression per cell type + volcano plots (Fig 5)
## Input: data/processed/neuro_with_ferroptosis_scores.rds
## Output: figures/stage-2/Fig5_Volcano_*.png + .pdf, data/processed/DE_results_*.csv
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

message("=== Step 04: DE and Volcano Plots ===")
message("Start time: ", Sys.time())

## =========================================================
## 0. Paths and helpers
## =========================================================
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "figures", "stage-2")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot_local <- function(p, filename, width = 8, height = 6, dpi = 300) {
  out_png <- file.path(fig_dir, paste0(filename, ".png"))
  out_pdf <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggsave(filename = out_png, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  tryCatch(
    ggsave(filename = out_pdf, plot = p, width = width, height = height, device = cairo_pdf),
    error = function(e) {
      warning("cairo_pdf failed, using pdf(): ", e$message)
      ggsave(filename = out_pdf, plot = p, width = width, height = height, device = "pdf")
    }
  )
  message("  [Saved] ", filename)
}

## =========================================================
## 1. Load data
## =========================================================
message("[1/5] Loading Seurat object...")
neuro <- readRDS(file.path(data_dir, "neuro_with_ferroptosis_scores.rds"))
DefaultAssay(neuro) <- "RNA"
message("  Loaded: ", ncol(neuro), " cells")

neuro_order <- c("qNSC", "nIPC", "Neuroblast", "GC")

## =========================================================
## 2. DE function with Seurat v5 compatibility
## =========================================================
run_de_and_volcano <- function(neuro_obj, ct, fig_dir) {
  message("\n  [DE] Processing: ", ct)
  
  # Subset to this cell type
  sub <- subset(neuro_obj, subset = Celltype_Article == ct)
  
  # Check we have both age groups
  age_counts <- table(sub$Age_collapsed)
  if (!("Young" %in% names(age_counts)) || !("Old" %in% names(age_counts))) {
    warning("  Skipping ", ct, ": missing Young or Old cells")
    return(NULL)
  }
  
  message("    Cells: Young=", age_counts["Young"], ", Old=", age_counts["Old"])
  
  # Set Idents
  sub$Age_collapsed <- factor(sub$Age_collapsed, levels = c("Young", "Old"))
  Idents(sub) <- sub$Age_collapsed
  
  # Run FindMarkers
  message("    Running FindMarkers (wilcox)...")
  de <- FindMarkers(sub, ident.1 = "Old", ident.2 = "Young", test.use = "wilcox", slot = "data")
  
  # Handle Seurat v5 vs v3/v4 column name
  if ("avg_log2FC" %in% colnames(de)) {
    colnames(de)[colnames(de) == "avg_log2FC"] <- "avg_logFC"
  }
  
  # Add gene column
  de$gene <- rownames(de)
  
  # Save CSV
  csv_file <- file.path(data_dir, paste0("DE_results_", ct, ".csv"))
  write.csv(de, csv_file, row.names = TRUE)
  message("    Saved: ", basename(csv_file))
  
  # Volcano plot
  de$significance <- "NS"
  de$significance[de$p_val_adj < 0.05 & de$avg_logFC > 0] <- "Up"
  de$significance[de$p_val_adj < 0.05 & de$avg_logFC < 0] <- "Down"
  
  # Cap -log10(p_val_adj) for visualization
  de$neg_log10_padj <- -log10(de$p_val_adj)
  de$neg_log10_padj[is.infinite(de$neg_log10_padj)] <- max(de$neg_log10_padj[!is.infinite(de$neg_log10_padj)], na.rm = TRUE) + 1
  
  # Top genes for labeling (after adding neg_log10_padj)
  top_up   <- de %>% filter(significance == "Up")   %>% arrange(p_val_adj) %>% head(10)
  top_down <- de %>% filter(significance == "Down") %>% arrange(p_val_adj) %>% head(10)
  top_genes <- bind_rows(top_up, top_down)
  
  sig_colors <- c("Up" = "#D73027", "Down" = "#3182BD", "NS" = "grey60")
  
  p <- ggplot(de, aes(x = avg_logFC, y = neg_log10_padj, color = significance)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = sig_colors, name = "Significance") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20) +
    ggtitle(paste0("Figure 5. Volcano: ", ct, " (Old vs Young)")) +
    xlab("avg_logFC") + ylab("-log10(adjusted p-value)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  save_plot_local(p, paste0("Fig5_Volcano_", ct), width = 8, height = 6)
  
  # Cleanup
  rm(sub)
  gc()
  
  return(de)
}

## =========================================================
## 3. Run DE for each cell type
## =========================================================
message("[2/5] Running DE analysis per cell type...")
de_results <- list()

for (ct in neuro_order) {
  de_results[[ct]] <- run_de_and_volcano(neuro, ct, fig_dir)
}

gc()

## =========================================================
## 4. Summary
## =========================================================
message("\n[3/5] DE Summary:")
for (ct in neuro_order) {
  if (!is.null(de_results[[ct]])) {
    de <- de_results[[ct]]
    n_up   <- sum(de$p_val_adj < 0.05 & de$avg_logFC > 0, na.rm = TRUE)
    n_down <- sum(de$p_val_adj < 0.05 & de$avg_logFC < 0, na.rm = TRUE)
    message(sprintf("  %s: %d up, %d down (FDR < 0.05)", ct, n_up, n_down))
  }
}

## =========================================================
## 5. List generated files
## =========================================================
message("\n[4/5] Generated files:")
volcano_files <- list.files(fig_dir, pattern = "^Fig5_Volcano", full.names = FALSE)
csv_files     <- list.files(data_dir, pattern = "^DE_results_", full.names = FALSE)
for (f in c(volcano_files, csv_files)) message("  ", f)

message("\n=== Step 04 Complete ===")
message("End time: ", Sys.time())
