## =========================================================
## Script: 02_module_scores.R
## Purpose: Calculate ferroptosis module scores (Promoter/Inhibitor)
## Input: data/processed/neuro_processed.rds
## Output: data/processed/neuro_with_ferroptosis_scores.rds, gene_lists/
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

message("=== Step 02: Module Scores ===")
message("Start time: ", Sys.time())

## =========================================================
## 0. Paths
## =========================================================
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data", "processed")
output_dir  <- file.path(project_dir, "data", "processed")
gene_dir    <- file.path(project_dir, "gene_lists")
dir.create(gene_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1. Load processed Seurat object
## =========================================================
message("[1/6] Loading neuro_processed.rds...")
neuro <- readRDS(file.path(data_dir, "neuro_processed.rds"))
DefaultAssay(neuro) <- "RNA"
message("  Loaded: ", ncol(neuro), " cells")

## =========================================================
## 2. Define ferroptosis gene sets (Human symbols)
## =========================================================
message("[2/6] Defining ferroptosis gene sets...")

ferroptosis_genes <- list(
  Ferro_Promoter = c(
    "TFRC", "NCOA4", "HMOX1", "PHKG2", "ACO1", "IREB2",
    "ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15",
    "DPP4", "ACSF2", "ZEB1", "PTGS2", "PEBP1",
    "CARS", "CHAC1", "NOX1", "ABCC1", "NOX4",
    "GLS2", "ATP5G3", "VDAC3", "SAT1"
  ),
  Ferro_Inhibitor = c(
    "FTH1", "NFS1", "MT1G", "STEAP3", "HSBP1", "FANCD2", "TERF1",
    "GPX4", "AKR1C1", "AKR1C2", "AKR1C3", "HMGCR", "SQLE", "CISD1", "ACSL3",
    "CS", "CBS", "FDFT1", "FADS2", "ACACA",
    "SLC7A11", "GCLC", "GCLM", "GSS", "NFE2L2", "KEAP1", "SLC38A1",
    "G6PD", "PGD", "SLC1A5", "GOT1",
    "CD44", "CRYAB", "EMC2", "HSPB1", "AIFM2", "RPL8"
  ),
  Ferro_Regulator = c("NQO1", "VDAC2", "TP53", "RAS")
)

## =========================================================
## 3. Human → Mouse gene symbol conversion
## =========================================================
message("[3/6] Converting gene symbols...")

human_to_mouse_symbol <- function(genes) {
  genes_lower <- tolower(genes)
  paste0(toupper(substr(genes_lower, 1, 1)), substr(genes_lower, 2, nchar(genes_lower)))
}

ferroptosis_genes_mouse <- lapply(ferroptosis_genes, human_to_mouse_symbol)

# Intersect with genes present in dataset
genes_present <- rownames(neuro)
ferroptosis_genes_use <- lapply(ferroptosis_genes_mouse, function(g) intersect(g, genes_present))

# Report
message("  Gene retention after intersection:")
for (name in names(ferroptosis_genes_use)) {
  orig <- length(ferroptosis_genes_mouse[[name]])
  kept <- length(ferroptosis_genes_use[[name]])
  message(sprintf("    %s: %d/%d genes retained", name, kept, orig))
  
  missing <- setdiff(ferroptosis_genes_mouse[[name]], genes_present)
  if (length(missing) > 0) {
    message("      Missing: ", paste(missing, collapse = ", "))
  }
}

# Validation
if (length(ferroptosis_genes_use$Ferro_Promoter) < 20) {
  warning("Less than 20 Promoter genes found. Check gene symbol format.")
}
if (length(ferroptosis_genes_use$Ferro_Inhibitor) < 25) {
  warning("Less than 25 Inhibitor genes found. Check gene symbol format.")
}

## =========================================================
## 4. Calculate module scores
## =========================================================
message("[4/6] Calculating module scores...")
set.seed(1)

# Promoter score
if (!any(grepl("^FerroPromoter", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Promoter),
    name     = "FerroPromoter",
    assay    = "RNA"
  )
  message("  FerroPromoter score calculated")
} else {
  message("  FerroPromoter score already exists, skipping")
}

# Inhibitor score
if (!any(grepl("^FerroInhibitor", colnames(neuro@meta.data)))) {
  neuro <- AddModuleScore(
    object   = neuro,
    features = list(ferroptosis_genes_use$Ferro_Inhibitor),
    name     = "FerroInhibitor",
    assay    = "RNA"
  )
  message("  FerroInhibitor score calculated")
} else {
  message("  FerroInhibitor score already exists, skipping")
}

# Get column names (with numeric suffix from AddModuleScore)
promoter_col  <- colnames(neuro@meta.data)[grepl("^FerroPromoter", colnames(neuro@meta.data))][1]
inhibitor_col <- colnames(neuro@meta.data)[grepl("^FerroInhibitor", colnames(neuro@meta.data))][1]
message("  Promoter score column: ", promoter_col)
message("  Inhibitor score column: ", inhibitor_col)

## =========================================================
## 5. Summary statistics
## =========================================================
message("[5/6] Summary statistics...")

message("  Promoter score: mean=", round(mean(neuro@meta.data[[promoter_col]], na.rm=TRUE), 4),
        " sd=", round(sd(neuro@meta.data[[promoter_col]], na.rm=TRUE), 4))
message("  Inhibitor score: mean=", round(mean(neuro@meta.data[[inhibitor_col]], na.rm=TRUE), 4),
        " sd=", round(sd(neuro@meta.data[[inhibitor_col]], na.rm=TRUE), 4))

# By cell type and age
score_summary <- neuro@meta.data %>%
  group_by(Celltype_Article, Age_collapsed) %>%
  summarise(
    mean_promoter  = mean(.data[[promoter_col]], na.rm = TRUE),
    mean_inhibitor = mean(.data[[inhibitor_col]], na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
print(as.data.frame(score_summary))

## =========================================================
## 6. Save outputs
## =========================================================
message("[6/6] Saving outputs...")

# Save Seurat object
saveRDS(neuro, file.path(output_dir, "neuro_with_ferroptosis_scores.rds"))
message("  Saved: neuro_with_ferroptosis_scores.rds")

# Save gene lists
writeLines(ferroptosis_genes_use$Ferro_Promoter,
           file.path(gene_dir, "ferroptosis_promoter.txt"))
writeLines(ferroptosis_genes_use$Ferro_Inhibitor,
           file.path(gene_dir, "ferroptosis_inhibitor.txt"))
writeLines(unique(ferroptosis_genes_use$Ferro_Regulator),
           file.path(gene_dir, "ferroptosis_regulator.txt"))
message("  Saved: gene_lists/")

gc()

message("\n=== Step 02 Complete ===")
message("End time: ", Sys.time())
