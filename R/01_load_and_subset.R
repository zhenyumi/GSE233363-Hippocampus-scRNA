## =========================================================
## Script: 01_load_and_subset.R
## Purpose: Load Seurat object, collapse age groups, subset neurogenic lineage
## Input: data/raw/GSE233363_custom/Seurat_combined_with_Celltype_Article.rds
##        Note: data/raw/GSE233363_official/ contains author-provided Zenodo files (for reference)
## Output: data/processed/combined_processed.rds, neuro_processed.rds
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

message("=== Step 01: Load and Subset ===")
message("Start time: ", Sys.time())

## =========================================================
## 0. Paths
## =========================================================
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data", "raw", "GSE233363_custom")
output_dir  <- file.path(project_dir, "data", "processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1. Load Seurat object
## =========================================================
seurat_file <- file.path(data_dir, "Seurat_combined_with_Celltype_Article.rds")
if (!file.exists(seurat_file)) {
  stop("Seurat object not found: ", seurat_file,
       "\nPlease place the file in data/raw/GSE233363_custom/")
}

message("[1/5] Loading Seurat object...")
combined <- readRDS(seurat_file)
DefaultAssay(combined) <- "RNA"
message("  Loaded: ", ncol(combined), " cells x ", nrow(combined), " genes")

## =========================================================
## 2. Runtime checks
## =========================================================
message("[2/5] Validating metadata...")

# Check required columns
required_cols <- c("Celltype_Article", "timepoint")
for (col in required_cols) {
  if (!col %in% colnames(combined@meta.data)) {
    stop("Required metadata column missing: ", col,
         "\nAvailable columns: ", paste(colnames(combined@meta.data), collapse = ", "))
  }
}

# Check UMAP reduction
if (!"umap" %in% names(combined@reductions)) {
  warning("UMAP reduction not found. Fig 0 and Fig 1 will fail.",
          "\nAvailable reductions: ", paste(names(combined@reductions), collapse = ", "))
}

# Print timepoint values
message("  timepoint values:")
print(table(combined$timepoint, useNA = "ifany"))

# Print cell type values
message("  Celltype_Article values:")
print(table(combined$Celltype_Article, useNA = "ifany"))

# Check gene symbol format
message("  First 10 gene names: ", paste(head(rownames(combined), 10), collapse = ", "))

## =========================================================
## 3. Age grouping
## =========================================================
message("[3/5] Creating Age_collapsed...")

combined$Age_collapsed <- NA_character_
combined$Age_collapsed[combined$timepoint %in% c("Young", "Young1", "Young2")] <- "Young"
combined$Age_collapsed[combined$timepoint %in% c("Old", "Old1", "Old2")] <- "Old"
combined$Age_collapsed <- factor(combined$Age_collapsed, levels = c("Young", "Old"))

message("  Age distribution:")
print(table(combined$Age_collapsed, useNA = "ifany"))

# Save full combined object
saveRDS(combined, file.path(output_dir, "combined_processed.rds"))
message("  Saved: combined_processed.rds")

## =========================================================
## 4. Subset neurogenic lineage
## =========================================================
message("[4/5] Subsetting neurogenic lineage...")

neuro_ct <- c("qNSC", "nIPC", "Neuroblast", "GC")

# Check which cell types exist
ct_present <- neuro_ct[neuro_ct %in% unique(combined$Celltype_Article)]
ct_missing <- neuro_ct[!neuro_ct %in% unique(combined$Celltype_Article)]

if (length(ct_missing) > 0) {
  warning("Missing neurogenic cell types: ", paste(ct_missing, collapse = ", "))
}

message("  Neurogenic cell types found: ", paste(ct_present, collapse = ", "))
message("  Total cells before subset: ", ncol(combined))

neuro <- subset(combined, subset = Celltype_Article %in% neuro_ct)
neuro <- subset(neuro, subset = !is.na(Age_collapsed) & Age_collapsed %in% c("Young", "Old"))
neuro$Age_collapsed    <- factor(neuro$Age_collapsed, levels = c("Young", "Old"))
neuro$Celltype_Article <- factor(neuro$Celltype_Article, levels = neuro_ct)

message("  Neurogenic cells after subset: ", ncol(neuro))
message("  Cell type x Age breakdown:")
print(table(neuro$Celltype_Article, neuro$Age_collapsed, useNA = "ifany"))

## =========================================================
## 5. Save and cleanup
## =========================================================
message("[5/5] Saving neuro_processed.rds...")
saveRDS(neuro, file.path(output_dir, "neuro_processed.rds"))
message("  Saved: neuro_processed.rds")

# Release memory
rm(combined)
gc()

message("\n=== Step 01 Complete ===")
message("End time: ", Sys.time())
message("Output: ", file.path(output_dir, "neuro_processed.rds"))
