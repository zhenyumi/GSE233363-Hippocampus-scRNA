## =========================================================
## Script: 03_figures_0_to_4.R
## Purpose: Generate Fig 0 (global UMAP), Fig 1-4 (neurogenic analysis)
## Input: data/processed/neuro_with_ferroptosis_scores.rds
## Output: figures/stage-2/Fig0-4*.png + .pdf
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

message("=== Step 03: Figures 0-4 ===")
message("Start time: ", Sys.time())

## =========================================================
## 0. Paths and helpers
## =========================================================
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "figures", "stage-2")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot_local <- function(p, filename, width = 8, height = 5, dpi = 300) {
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
message("[1/8] Loading Seurat object...")
neuro <- readRDS(file.path(data_dir, "neuro_with_ferroptosis_scores.rds"))
DefaultAssay(neuro) <- "RNA"
message("  Loaded: ", ncol(neuro), " cells")

# Get score column names
promoter_col  <- colnames(neuro@meta.data)[grepl("^FerroPromoter", colnames(neuro@meta.data))][1]
inhibitor_col <- colnames(neuro@meta.data)[grepl("^FerroInhibitor", colnames(neuro@meta.data))][1]
message("  Promoter col: ", promoter_col)
message("  Inhibitor col: ", inhibitor_col)

# Get regulator genes
gene_dir      <- file.path(project_dir, "gene_lists")
regulator_genes <- readLines(file.path(gene_dir, "ferroptosis_regulator.txt"))
regulator_genes <- regulator_genes[regulator_genes %in% rownames(neuro)]
message("  Regulator genes: ", paste(regulator_genes, collapse = ", "))

neuro_order <- c("qNSC", "nIPC", "Neuroblast", "GC")

# Color scheme
celltype_colors <- c(
  "qNSC"       = "#FDAE6B",
  "nIPC"       = "#3182BD",
  "Neuroblast" = "#08519C",
  "GC"         = "#9ECAE1"
)
age_colors <- c("Young" = "#4DBBD5", "Old" = "#D73027")

## =========================================================
## 2. Figure 0: Global UMAP (all cell types)
## =========================================================
message("[2/8] Figure 0: Global UMAP...")

# Load full combined object for Fig 0
combined_file <- file.path(project_dir, "data", "processed", "combined_processed.rds")
if (file.exists(combined_file)) {
  combined <- readRDS(combined_file)
  
  if ("umap" %in% names(combined@reductions)) {
    p_fig0 <- DimPlot(combined, reduction = "umap", group.by = "Celltype_Article", label = TRUE) +
      ggtitle("Figure 0. All cell types UMAP") +
      theme(plot.title = element_text(hjust = 0.5))
    save_plot_local(p_fig0, "Fig0_UMAP_AllCelltypes", width = 10, height = 8)
  } else {
    warning("No UMAP reduction in combined object. Skipping Fig 0.")
  }
  
  rm(combined)
  gc()
} else {
  warning("combined_processed.rds not found. Skipping Fig 0. Run 01_load_and_subset.R first.")
}

## =========================================================
## 3. Figure 1: Neurogenic lineage UMAP split by age
## =========================================================
message("[3/8] Figure 1: Neurogenic lineage UMAP...")

if ("umap" %in% names(neuro@reductions)) {
  # Add cell counts to age labels
  age_n <- neuro@meta.data %>%
    dplyr::count(Age_collapsed, name = "n_cells") %>%
    dplyr::mutate(Age_strip = paste0(as.character(Age_collapsed), " (n=", n_cells, ")"))
  
  neuro$Age_collapsed_strip <- as.character(neuro$Age_collapsed)
  for (i in seq_len(nrow(age_n))) {
    neuro$Age_collapsed_strip[neuro$Age_collapsed_strip == as.character(age_n$Age_collapsed[i])] <-
      age_n$Age_strip[i]
  }
  neuro$Age_collapsed_strip <- factor(
    neuro$Age_collapsed_strip,
    levels = age_n$Age_strip[match(c("Young", "Old"), as.character(age_n$Age_collapsed))]
  )
  
  # Add cell counts to cell type labels
  ct_age_n <- neuro@meta.data %>%
    dplyr::count(Celltype_Article, Age_collapsed, name = "n_cells") %>%
    tidyr::pivot_wider(names_from = Age_collapsed, values_from = n_cells, values_fill = 0) %>%
    dplyr::mutate(legend_label = paste0(Celltype_Article, " (Y=", Young, ", O=", Old, ")"))
  
  neuro$Celltype_Article_legend <- as.character(neuro$Celltype_Article)
  for (i in seq_len(nrow(ct_age_n))) {
    ct <- as.character(ct_age_n$Celltype_Article[i])
    lb <- as.character(ct_age_n$legend_label[i])
    neuro$Celltype_Article_legend[neuro$Celltype_Article_legend == ct] <- lb
  }
  neuro$Celltype_Article_legend <- factor(
    neuro$Celltype_Article_legend,
    levels = ct_age_n$legend_label[match(neuro_order, ct_age_n$Celltype_Article)]
  )
  
  celltype_colors_legend <- celltype_colors[neuro_order]
  names(celltype_colors_legend) <- ct_age_n$legend_label[match(neuro_order, ct_age_n$Celltype_Article)]
  
  p_fig1 <- DimPlot(
    neuro, reduction = "umap", group.by = "Celltype_Article_legend",
    split.by = "Age_collapsed_strip", label = FALSE
  ) +
    scale_color_manual(values = celltype_colors_legend, name = "Cell type") +
    ggtitle("Figure 1. Neurogenic lineage UMAP split by age") +
    theme(plot.title = element_text(hjust = 0.5), strip.text = element_text(face = "bold"))
  
  save_plot_local(p_fig1, "Fig1_UMAP_Neuro_Celltype_SplitByAge", width = 10, height = 5)
} else {
  warning("No UMAP reduction. Skipping Fig 1.")
}

gc()

## =========================================================
## 4. Figure 2: Violin plots
## =========================================================
message("[4/8] Figure 2: Violin plots...")

# Fig 2A: Promoter score
p_fig2A <- VlnPlot(neuro, features = promoter_col, group.by = "Celltype_Article",
                    split.by = "Age_collapsed", pt.size = 0) +
  scale_fill_manual(values = age_colors, drop = TRUE) +
  ggtitle("Figure 2A. Promoter module score") +
  theme(plot.title = element_text(hjust = 0.5))
save_plot_local(p_fig2A, "Fig2A_Violin_PromoterScore", width = 8, height = 5)

# Fig 2B: Inhibitor score
p_fig2B <- VlnPlot(neuro, features = inhibitor_col, group.by = "Celltype_Article",
                    split.by = "Age_collapsed", pt.size = 0) +
  scale_fill_manual(values = age_colors, drop = TRUE) +
  ggtitle("Figure 2B. Inhibitor module score") +
  theme(plot.title = element_text(hjust = 0.5))
save_plot_local(p_fig2B, "Fig2B_Violin_InhibitorScore", width = 8, height = 5)

# Fig 2C: Regulator genes (with guard)
if (length(regulator_genes) > 0) {
  p_fig2C <- VlnPlot(neuro, features = regulator_genes, group.by = "Celltype_Article",
                      split.by = "Age_collapsed", pt.size = 0,
                      ncol = min(2, length(regulator_genes))) +
    scale_fill_manual(values = age_colors, drop = TRUE) +
    theme(plot.title = element_text(hjust = 0.5))
  save_plot_local(p_fig2C, "Fig2C_Violin_RegulatorGenes", width = 10, height = 5)
}

gc()

## =========================================================
## 5. Figure 3: Line plots
## =========================================================
message("[5/8] Figure 3: Line plots...")

plot_line_summary <- function(neuro_obj, value_col, title_text, ylab_text, out_name,
                              width = 6, height = 5) {
  df <- FetchData(neuro_obj, vars = c(value_col, "Age_collapsed", "Celltype_Article"))
  colnames(df)[1] <- "Value"
  
  df2 <- df %>%
    mutate(Age_collapsed = factor(Age_collapsed, levels = c("Young", "Old"))) %>%
    group_by(Celltype_Article, Age_collapsed) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), n_cells = n(), .groups = "drop")
  
  p <- ggplot(df2, aes(x = Age_collapsed, y = mean_value, group = Celltype_Article, color = Celltype_Article)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    scale_color_manual(values = celltype_colors, name = "Cell type") +
    theme_bw() + labs(x = "Age group", y = ylab_text, title = title_text) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())
  
  save_plot_local(p, out_name, width = width, height = height)
}

# Fig 3A: Promoter
plot_line_summary(neuro, promoter_col,
                  "Figure 3A. Promoter score: Young vs Old",
                  "Mean promoter module score",
                  "Fig3A_Line_PromoterScore")

# Fig 3B: Inhibitor
plot_line_summary(neuro, inhibitor_col,
                  "Figure 3B. Inhibitor score: Young vs Old",
                  "Mean inhibitor module score",
                  "Fig3B_Line_InhibitorScore")

# Fig 3C: Regulator genes (with guard)
if (length(regulator_genes) > 0) {
  for (g in regulator_genes) {
    plot_line_summary(neuro, g,
                      paste0("Figure 3C. ", g, " expression: Young vs Old"),
                      paste0("Mean expression of ", g),
                      paste0("Fig3C_Line_RegulatorGene_", g))
  }
}

gc()

## =========================================================
## 6. Figure 4: Trajectory plots (median)
## =========================================================
message("[6/8] Figure 4: Trajectory plots...")

plot_traj_median <- function(neuro_obj, value_col, title_text, ylab_text, out_name,
                             width = 6.5, height = 4.5) {
  df <- FetchData(neuro_obj, vars = c(value_col, "Age_collapsed", "Celltype_Article"))
  colnames(df)[1] <- "Value"
  df <- df %>%
    filter(Celltype_Article %in% neuro_order) %>%
    mutate(
      Celltype_Article = factor(Celltype_Article, levels = neuro_order),
      Age_collapsed    = factor(Age_collapsed, levels = c("Young", "Old"))
    ) %>%
    group_by(Celltype_Article, Age_collapsed) %>%
    summarise(median_value = median(Value, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(df, aes(x = Celltype_Article, y = median_value, color = Age_collapsed, group = Age_collapsed)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
    scale_color_manual(values = age_colors, name = "Age") +
    theme_bw() + labs(title = title_text, x = "Neurogenic lineage", y = ylab_text) +
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(face = "bold"),
          panel.grid = element_blank())
  
  save_plot_local(p, out_name, width = width, height = height)
}

# Fig 4A: Promoter trajectory
plot_traj_median(neuro, promoter_col,
                 "Figure 4A. Trajectory: Promoter activity",
                 "Median promoter module score",
                 "Fig4A_Trajectory_PromoterScore")

# Fig 4B: Inhibitor trajectory
plot_traj_median(neuro, inhibitor_col,
                 "Figure 4B. Trajectory: Inhibitor activity",
                 "Median inhibitor module score",
                 "Fig4B_Trajectory_InhibitorScore")

# Fig 4C: Regulator genes trajectory (with guard)
if (length(regulator_genes) > 0) {
  for (g in regulator_genes) {
    if (!g %in% rownames(neuro)) next
    plot_traj_median(neuro, g,
                     paste0("Figure 4C. Trajectory: ", g),
                     paste0("Median expression of ", g),
                     paste0("Fig4C_Trajectory_Regulator_", g))
  }
}

gc()

## =========================================================
## 7. Summary
## =========================================================
message("[7/8] Listing generated figures...")
fig_files <- list.files(fig_dir, pattern = "^Fig[0-4]", full.names = FALSE)
message("  Generated ", length(fig_files), " files:")
for (f in fig_files) message("    ", f)

message("\n=== Step 03 Complete ===")
message("End time: ", Sys.time())
message("Output directory: ", fig_dir)
