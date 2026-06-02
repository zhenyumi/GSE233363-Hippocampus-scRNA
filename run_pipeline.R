## =========================================================
## Script: run_pipeline.R
## Purpose: Initialize renv and run the full analysis pipeline
## Usage: Rscript run_pipeline.R
## =========================================================

message("========================================")
message("ferroDG Reproduction Pipeline")
message("========================================")
message("Start time: ", Sys.time())
message("Working directory: ", getwd())
message("R version: ", R.version.string)

## =========================================================
## 1. renv initialization
## =========================================================
message("\n[1/3] Setting up renv environment...")

if (!requireNamespace("renv", quietly = TRUE)) {
  message("  Installing renv...")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialize renv if not already done
if (!file.exists("renv.lock")) {
  message("  Initializing renv project...")
  renv::init(bare = TRUE)
}

# Install required packages
required_pkgs <- c("Seurat", "dplyr", "tidyr", "ggplot2", "ggrepel", "pheatmap")
installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs, installed)

if (length(to_install) > 0) {
  message("  Installing missing packages: ", paste(to_install, collapse = ", "))
  renv::install(to_install)
} else {
  message("  All required packages already installed")
}

# Snapshot current state
message("  Creating renv.lock snapshot...")
renv::snapshot(prompt = FALSE)
message("  renv setup complete")

## =========================================================
## 2. Run analysis pipeline
## =========================================================
message("\n[2/3] Running analysis pipeline...")

scripts <- c(
  "R/01_load_and_subset.R",
  "R/02_module_scores.R",
  "R/03_figures_0_to_4.R",
  "R/04_de_and_volcano.R",
  "R/05_correlation_heatmaps.R"
)

for (script in scripts) {
  if (!file.exists(script)) {
    stop("Script not found: ", script)
  }
  message("\n--- Running: ", script, " ---")
  source(script)
  message("--- Completed: ", script, " ---\n")
}

## =========================================================
## 3. Final summary
## =========================================================
message("\n[3/3] Pipeline Summary")
message("========================================")

# Check outputs
processed_files <- list.files("data/processed", pattern = "\\.rds$|\\.csv$", full.names = FALSE)
figure_files    <- list.files("figures/stage-2", pattern = "\\.png$", full.names = FALSE)
gene_files      <- list.files("gene_lists", pattern = "\\.txt$", full.names = FALSE)

message("Processed data files (", length(processed_files), "):")
for (f in processed_files) message("  ", f)

message("\nFigure files (", length(figure_files), "):")
for (f in figure_files) message("  ", f)

message("\nGene list files (", length(gene_files), "):")
for (f in gene_files) message("  ", f)

# Final renv snapshot
message("\nUpdating renv.lock...")
renv::snapshot(prompt = FALSE)

message("\n========================================")
message("Pipeline complete!")
message("End time: ", Sys.time())
message("========================================")
