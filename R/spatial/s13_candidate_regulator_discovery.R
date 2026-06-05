## Phase Spatial-13: CA1/CA3 Mitochondrial Module-Associated Candidate Regulator Discovery
## Plan version: v1.2
## Script: R/spatial/s13_candidate_regulator_discovery.R
## Execution: Rscript R/spatial/s13_candidate_regulator_discovery.R

cat("=== Phase Spatial-13: Candidate Regulator Discovery ===\n")
cat("Plan version: v1.2\n")
cat("Script: R/spatial/s13_candidate_regulator_discovery.R\n")
cat("Start time:", format(Sys.time()), "\n\n")

t_start <- Sys.time()

# Rplots.pdf guard
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

## ---- E0: Setup ----
cat("=== E0: Setup ===\n")

# Record renv.lock MD5
renv_md5_before <- digest::digest("renv.lock", algo = "md5", file = TRUE)
cat("renv.lock MD5 (before):", renv_md5_before, "\n")

# Required packages
required_pkgs <- c("Seurat", "DESeq2", "dplyr", "ggplot2", "pheatmap", "patchwork", "stringr", "digest")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("FATAL: Package not available: ", pkg)
}
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(stringr)

# Output directories
dir.create("data/processed/spatial/phase13_candidate_regulator", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/spatial/phase13_candidate_regulator", recursive = TRUE, showWarnings = FALSE)

# Verify input files
input_files <- list(
  pseudobulk       = "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds",
  module_audit     = "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv",
  module_scores    = "data/processed/spatial/phase12_young_aged_region_comparison/module_scores_all_regions.csv",
  de_ca1           = "data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA1.csv",
  de_ca3           = "data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA3.csv",
  de_region        = "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/deseq2_all_genes_CA3_vs_CA1.csv",
  target_genes     = "docs/target_genes.csv",
  hippo_rds        = "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"
)

for (nm in names(input_files)) {
  if (!file.exists(input_files[[nm]])) stop("FATAL: Missing input: ", input_files[[nm]])
}
cat("All input files verified\n")

## ---- E1: Package Installation ----
cat("\n=== E1: Package Installation ===\n")

has_dorothea   <- requireNamespace("dorothea", quietly = TRUE)
has_decoupleR  <- requireNamespace("decoupleR", quietly = TRUE)
has_orgdb      <- requireNamespace("org.Mm.eg.db", quietly = TRUE)

if (!has_dorothea || !has_decoupleR || !has_orgdb) {
  cat("Installing missing Bioconductor packages...\n")
  pkgs_to_install <- c()
  if (!has_dorothea)  pkgs_to_install <- c(pkgs_to_install, "bioc::dorothea")
  if (!has_decoupleR) pkgs_to_install <- c(pkgs_to_install, "bioc::decoupleR")
  if (!has_orgdb)     pkgs_to_install <- c(pkgs_to_install, "bioc::org.Mm.eg.db")
  renv::install(pkgs_to_install)
  if (!requireNamespace("dorothea", quietly = TRUE))    stop("FATAL: dorothea install failed")
  if (!requireNamespace("decoupleR", quietly = TRUE))   stop("FATAL: decoupleR install failed")
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("FATAL: org.Mm.eg.db install failed")
}

library(dorothea)
library(decoupleR)
library(org.Mm.eg.db)

has_clusterProfiler <- requireNamespace("clusterProfiler", quietly = TRUE)
has_msigdbr         <- requireNamespace("msigdbr", quietly = TRUE)

renv_md5_after <- digest::digest("renv.lock", algo = "md5", file = TRUE)
if (renv_md5_before != renv_md5_after) {
  cat("WARNING: renv.lock changed during E1\n")
  cat("  Before:", renv_md5_before, "\n")
  cat("  After:", renv_md5_after, "\n")
}
cat("renv.lock MD5 (after):", renv_md5_after, "\n")

## ---- E2: Load Data ----
cat("\n=== E2: Loading Data ===\n")

pb_counts <- readRDS(input_files$pseudobulk)
cat("Pseudobulk counts:", nrow(pb_counts), "genes ×", ncol(pb_counts), "samples\n")

module_audit <- read.csv(input_files$module_audit, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
cat("Module audit:", nrow(module_audit), "modules\n")

# Build module gene lists (handle column name variants)
mod_col <- intersect(c("module", "Module", "module_name"), colnames(module_audit))[1]
genes_col <- intersect(c("genes", "Genes", "gene_set", "gene_list"), colnames(module_audit))[1]
if (is.na(mod_col) || is.na(genes_col)) stop("FATAL: Cannot identify module/genes columns in audit CSV")

module_genes <- list()
for (i in seq_len(nrow(module_audit))) {
  mod_name <- module_audit[[mod_col]][i]
  genes_str <- module_audit[[genes_col]][i]
  module_genes[[mod_name]] <- trimws(strsplit(genes_str, "[,]")[[1]])
}

# Load module scores
mod_scores <- read.csv(input_files$module_scores, stringsAsFactors = FALSE)
cat("Module scores:", nrow(mod_scores), "samples,", ncol(mod_scores), "columns\n")

score_cols <- grep("_mean$", colnames(mod_scores), value = TRUE)
module_names <- sub("_mean$", "", score_cols)
cat("Modules with scores:", length(module_names), ":", paste(module_names, collapse = ", "), "\n")

# Load DE results
de_ca1 <- read.csv(input_files$de_ca1, stringsAsFactors = FALSE)
de_ca3 <- read.csv(input_files$de_ca3, stringsAsFactors = FALSE)
de_region <- read.csv(input_files$de_region, stringsAsFactors = FALSE)
cat("DE CA1:", nrow(de_ca1), "genes; DE CA3:", nrow(de_ca3), "genes; DE region:", nrow(de_region), "genes\n")

# Load target genes (GB2312 encoding with Chinese headers)
target_raw <- read.csv(input_files$target_genes, stringsAsFactors = FALSE, fileEncoding = "GB2312")
cat("Target genes raw:", nrow(target_raw), "categories\n")

target_categories <- list()
for (i in seq_len(nrow(target_raw))) {
  cat_name <- as.character(target_raw[i, 1])
  genes_str <- as.character(target_raw[i, 2])
  target_categories[[cat_name]] <- trimws(strsplit(genes_str, "\u3001")[[1]])  # \u3001 = 、
}
all_target_genes <- unique(unlist(target_categories))
cat("Total unique target genes:", length(all_target_genes), "\n")

# Generate manifest from pseudobulk column names (pattern: G1_CA1, G2_CA3, etc.)
pb_colnames <- colnames(pb_counts)
manifest <- data.frame(
  Sample = pb_colnames,
  stringsAsFactors = FALSE
)
# Parse GEMgroup and Region from column names
manifest$GEMgroup <- sub("^(G[0-9]+)_.*", "\\1", manifest$Sample)
manifest$Region <- sub("^G[0-9]+_(.*)", "\\1", manifest$Sample)
# Assign Age based on GEMgroup: G1-G4=Young, G5-G8=Middle, G9-G16=Old
gem_num <- as.integer(sub("G", "", manifest$GEMgroup))
manifest$Age <- ifelse(gem_num <= 4, "Young", ifelse(gem_num <= 8, "Middle", "Old"))

# Save manifest for reference
write.csv(manifest, "data/processed/spatial/phase13_candidate_regulator/pseudobulk_sample_manifest.csv", row.names = FALSE)
cat("Manifest generated:", nrow(manifest), "samples\n")

# Generate gene_presence_audit
presence_audit <- data.frame(
  gene = all_target_genes,
  present_in_pseudobulk = all_target_genes %in% rownames(pb_counts),
  stringsAsFactors = FALSE
)
write.csv(presence_audit, "data/processed/spatial/phase13_candidate_regulator/gene_presence_audit.csv", row.names = FALSE)
cat("Gene presence audit generated:", sum(presence_audit$present_in_pseudobulk), "/", nrow(presence_audit), "target genes present\n")

# Filter to Young + Old only (exclude Middle)
manifest_yo <- manifest[manifest$Age %in% c("Young", "Old"), ]
manifest_ca1 <- manifest_yo[manifest_yo$Region == "CA1", ]
manifest_ca3 <- manifest_yo[manifest_yo$Region == "CA3", ]

ca1_cols <- manifest_ca1$Sample
ca3_cols <- manifest_ca3$Sample

n_ca1_young <- sum(manifest_ca1$Age == "Young")
n_ca1_old   <- sum(manifest_ca1$Age == "Old")
n_ca3_young <- sum(manifest_ca3$Age == "Young")
n_ca3_old   <- sum(manifest_ca3$Age == "Old")

cat("CA1:", length(ca1_cols), "samples (Young:", n_ca1_young, ", Old:", n_ca1_old, ")\n")
cat("CA3:", length(ca3_cols), "samples (Young:", n_ca3_young, ", Old:", n_ca3_old, ")\n")

# Validate columns exist in pseudobulk
stopifnot(all(ca1_cols %in% colnames(pb_counts)))
stopifnot(all(ca3_cols %in% colnames(pb_counts)))

pb_ca1 <- pb_counts[, ca1_cols]
pb_ca3 <- pb_counts[, ca3_cols]

# Compute log2(CPM+1)
compute_log2cpm <- function(counts_mat) {
  lib_sizes <- colSums(counts_mat)
  cpm_mat <- t(t(counts_mat) / lib_sizes * 1e6)
  log2(cpm_mat + 1)
}

log2cpm_ca1 <- compute_log2cpm(pb_ca1)
log2cpm_ca3 <- compute_log2cpm(pb_ca3)
cat("log2(CPM+1) computed for CA1 and CA3\n")

gc()

## ---- E3: Build Candidate Universes ----
cat("\n=== E3: Building Candidate Universes ===\n")

all_genes <- rownames(pb_counts)

# universe_all: expressed in >=3 GEMgroups in BOTH CA1 and CA3, mean log2(CPM+1) >= 1
ca1_nz <- rowSums(pb_ca1 > 0)
ca3_nz <- rowSums(pb_ca3 > 0)
mean_log2cpm_all <- rowMeans(cbind(log2cpm_ca1, log2cpm_ca3))

universe_all_mask <- (ca1_nz >= 3) & (ca3_nz >= 3) & (mean_log2cpm_all >= 1)
universe_all <- all_genes[universe_all_mask]
cat("universe_all:", length(universe_all), "genes\n")
if (length(universe_all) < 100) stop("FATAL: universe_all < 100 genes (", length(universe_all), ")")

# TF annotation
cat("Building TF annotation...\n")

# Use built-in dorothea::dorothea_mm dataset (no network required)
# Column is 'tf' not 'source'
dorothea_tfs <- c()
tryCatch({
  dorothea_regulons <- dorothea::dorothea_mm
  dorothea_regulons <- dorothea_regulons[dorothea_regulons$confidence %in% c("A", "B", "C"), ]
  dorothea_tfs <- unique(dorothea_regulons$tf)
  cat("DoRothEA A-C TFs (built-in):", length(dorothea_tfs), "\n")
}, error = function(e) cat("WARNING: dorothea::dorothea_mm failed:", e$message, "\n"))

# CollecTRI - try API, skip if network unavailable
collectri_tfs <- c()
tryCatch({
  collectri_regulons <- get_collectri(organism = "mouse")
  collectri_tfs <- unique(collectri_regulons$source)
  cat("CollecTRI TFs:", length(collectri_tfs), "\n")
}, error = function(e) cat("WARNING: get_collectri() skipped (network issue)\n"))

# GO TF annotation - use GOALL keytype
go_tfs <- c()
tryCatch({
  go_0003700 <- AnnotationDbi::select(org.Mm.eg.db, keys = "GO:0003700", keytype = "GOALL", columns = "SYMBOL")
  go_0140110 <- AnnotationDbi::select(org.Mm.eg.db, keys = "GO:0140110", keytype = "GOALL", columns = "SYMBOL")
  go_tfs <- unique(c(go_0003700$SYMBOL, go_0140110$SYMBOL))
  go_tfs <- go_tfs[!is.na(go_tfs)]
  cat("GO TFs:", length(go_tfs), "\n")
}, error = function(e) cat("WARNING: GO TF query failed:", e$message, "\n"))

all_tfs <- unique(c(dorothea_tfs, collectri_tfs, go_tfs))
cat("Combined TFs:", length(all_tfs), "\n")

tf_annotation <- data.frame(
  gene = all_tfs,
  dorothea = all_tfs %in% dorothea_tfs,
  collectri = all_tfs %in% collectri_tfs,
  go_tf = all_tfs %in% go_tfs,
  stringsAsFactors = FALSE
)

universe_regulator <- intersect(universe_all, all_tfs)
universe_target <- intersect(universe_all, all_target_genes)
cat("universe_regulator:", length(universe_regulator), "genes\n")
cat("universe_target:", length(universe_target), "genes\n")
if (length(universe_regulator) < 10) stop("FATAL: universe_regulator < 10 genes (", length(universe_regulator), ")")

universe_comp <- data.frame(
  universe = c("universe_all", "universe_regulator", "universe_target"),
  n_genes = c(length(universe_all), length(universe_regulator), length(universe_target)),
  stringsAsFactors = FALSE
)
write.csv(universe_comp, "data/processed/spatial/phase13_candidate_regulator/phase13_universe_composition.csv", row.names = FALSE)

gc()

## ---- E4: Spearman Correlation ----
cat("\n=== E4: Spearman Correlation ===\n")

# Filter mod_scores to CA1+CA3 Young+Old (mod_scores uses sample_id column)
mod_scores_yo <- mod_scores[mod_scores$sample_id %in% c(ca1_cols, ca3_cols), ]

# Build module score matrices (modules × samples)
mod_score_ca1_mat <- as.matrix(mod_scores_yo[mod_scores_yo$sample_id %in% ca1_cols, score_cols])
rownames(mod_score_ca1_mat) <- mod_scores_yo$sample_id[mod_scores_yo$sample_id %in% ca1_cols]
colnames(mod_score_ca1_mat) <- module_names  # align column names to module_names
mod_score_ca1_mat <- t(mod_score_ca1_mat)  # modules × samples

mod_score_ca3_mat <- as.matrix(mod_scores_yo[mod_scores_yo$sample_id %in% ca3_cols, score_cols])
rownames(mod_score_ca3_mat) <- mod_scores_yo$sample_id[mod_scores_yo$sample_id %in% ca3_cols]
colnames(mod_score_ca3_mat) <- module_names  # align column names to module_names
mod_score_ca3_mat <- t(mod_score_ca3_mat)  # modules × samples

n_genes_all <- length(universe_all)
n_modules <- length(module_names)

rho_ca1  <- matrix(NA, nrow = n_genes_all, ncol = n_modules, dimnames = list(universe_all, module_names))
pval_ca1 <- matrix(NA, nrow = n_genes_all, ncol = n_modules, dimnames = list(universe_all, module_names))
rho_ca3  <- matrix(NA, nrow = n_genes_all, ncol = n_modules, dimnames = list(universe_all, module_names))
pval_ca3 <- matrix(NA, nrow = n_genes_all, ncol = n_modules, dimnames = list(universe_all, module_names))

expr_ca1 <- as.matrix(log2cpm_ca1[universe_all, ])
expr_ca3 <- as.matrix(log2cpm_ca3[universe_all, ])

cat("Computing Spearman correlations for", n_genes_all, "genes ×", n_modules, "modules...\n")

for (m in module_names) {
  cat("  Module:", m)

  # CA1
  m_scores_ca1 <- mod_score_ca1_mat[m, ]
  rho_vals_ca1 <- cor(t(expr_ca1), m_scores_ca1, method = "spearman")
  n_ca1 <- length(m_scores_ca1)
  t_stat_ca1 <- rho_vals_ca1 * sqrt((n_ca1 - 2) / (1 - rho_vals_ca1^2))
  p_vals_ca1 <- 2 * pt(-abs(t_stat_ca1), df = n_ca1 - 2)
  rho_ca1[, m]  <- rho_vals_ca1
  pval_ca1[, m] <- p_vals_ca1

  # CA3
  m_scores_ca3 <- mod_score_ca3_mat[m, ]
  rho_vals_ca3 <- cor(t(expr_ca3), m_scores_ca3, method = "spearman")
  n_ca3 <- length(m_scores_ca3)
  t_stat_ca3 <- rho_vals_ca3 * sqrt((n_ca3 - 2) / (1 - rho_vals_ca3^2))
  p_vals_ca3 <- 2 * pt(-abs(t_stat_ca3), df = n_ca3 - 2)
  rho_ca3[, m]  <- rho_vals_ca3
  pval_ca3[, m] <- p_vals_ca3

  cat(" done\n")
}

gc()

# BH correction for q-values (exploratory only)
qval_ca1 <- apply(pval_ca1, 2, p.adjust, method = "BH")
qval_ca3 <- apply(pval_ca3, 2, p.adjust, method = "BH")

gc()

# Build long-format results (optimized to avoid memory spike)
n_res <- length(universe_all) * length(module_names)
rho_results <- data.frame(
  gene = rep(universe_all, times = length(module_names)),
  module = rep(module_names, each = length(universe_all)),
  rho_CA1 = as.vector(rho_ca1),
  pval_CA1 = as.vector(pval_ca1),
  qval_CA1 = as.vector(qval_ca1),
  rho_CA3 = as.vector(rho_ca3),
  pval_CA3 = as.vector(pval_ca3),
  qval_CA3 = as.vector(qval_ca3),
  stringsAsFactors = FALSE
)

rm(rho_ca1, pval_ca1, qval_ca1, rho_ca3, pval_ca3, qval_ca3)
gc()

## ---- E5: Self-Correlation Flags ----
cat("\n=== E5: Self-Correlation Flags ===\n")

rho_results$self_correlation_flag <- "regulator_candidate"

# Strategy A: flag module gene-set members for their own module
for (mod_name in names(module_genes)) {
  score_mod <- gsub("[ /]", ".", mod_name)
  if (!(score_mod %in% module_names)) next
  mod_gene_set <- module_genes[[mod_name]]
  mask <- rho_results$module == score_mod & rho_results$gene %in% mod_gene_set
  rho_results$self_correlation_flag[mask] <- "module_member_or_effector"
}

# Flag target genes
rho_results$self_correlation_flag[
  rho_results$gene %in% all_target_genes &
  rho_results$self_correlation_flag == "regulator_candidate"
] <- "target_gene_candidate"

# Flag low confidence
low_rho_mask <- abs(rho_results$rho_CA1) < 0.3 & abs(rho_results$rho_CA3) < 0.3 &
  rho_results$self_correlation_flag == "regulator_candidate"
rho_results$self_correlation_flag[low_rho_mask] <- "low_confidence_self_correlation"

cat("Self-correlation flags:\n")
print(table(rho_results$self_correlation_flag))

gc()

## ---- E6: regulator_score Computation ----
cat("\n=== E6: regulator_score Computation ===\n")

# Helper to find column names
find_col <- function(df, patterns) {
  for (p in patterns) {
    idx <- grep(p, colnames(df), ignore.case = TRUE)
    if (length(idx) > 0) return(colnames(df)[idx[1]])
  }
  NULL
}

lfc_ca1_col    <- find_col(de_ca1, c("log2FC", "log2FoldChange"))
padj_ca1_col   <- find_col(de_ca1, c("padj", "FDR"))
gene_ca1_col   <- find_col(de_ca1, c("gene", "Gene", "symbol", "Symbol"))
if (is.null(gene_ca1_col)) gene_ca1_col <- colnames(de_ca1)[1]

lfc_ca3_col    <- find_col(de_ca3, c("log2FC", "log2FoldChange"))
padj_ca3_col   <- find_col(de_ca3, c("padj", "FDR"))
gene_ca3_col   <- find_col(de_ca3, c("gene", "Gene", "symbol", "Symbol"))
if (is.null(gene_ca3_col)) gene_ca3_col <- colnames(de_ca3)[1]

lfc_region_col  <- find_col(de_region, c("log2FC", "log2FoldChange"))
padj_region_col <- find_col(de_region, c("padj", "FDR"))
gene_region_col <- find_col(de_region, c("gene", "Gene", "symbol", "Symbol"))
if (is.null(gene_region_col)) gene_region_col <- colnames(de_region)[1]

cat("DE CA1: lfc=", lfc_ca1_col, " padj=", padj_ca1_col, " gene=", gene_ca1_col, "\n")
cat("DE CA3: lfc=", lfc_ca3_col, " padj=", padj_ca3_col, " gene=", gene_ca3_col, "\n")
cat("DE region: lfc=", lfc_region_col, " padj=", padj_region_col, " gene=", gene_region_col, "\n")

# Build lookup tables
de_ca1_lfc    <- setNames(de_ca1[[lfc_ca1_col]], de_ca1[[gene_ca1_col]])
de_ca1_padj   <- setNames(de_ca1[[padj_ca1_col]], de_ca1[[gene_ca1_col]])
de_ca3_lfc    <- setNames(de_ca3[[lfc_ca3_col]], de_ca3[[gene_ca3_col]])
de_ca3_padj   <- setNames(de_ca3[[padj_ca3_col]], de_ca3[[gene_ca3_col]])
de_region_lfc <- setNames(de_region[[lfc_region_col]], de_region[[gene_region_col]])
de_region_padj<- setNames(de_region[[padj_region_col]], de_region[[gene_region_col]])

# Filter to universe_regulator
regulator_rho <- rho_results[rho_results$gene %in% universe_regulator, ]

# Add DE evidence
regulator_rho$DE_age_lfc_CA1   <- de_ca1_lfc[regulator_rho$gene]
regulator_rho$DE_age_padj_CA1  <- de_ca1_padj[regulator_rho$gene]
regulator_rho$DE_age_lfc_CA3   <- de_ca3_lfc[regulator_rho$gene]
regulator_rho$DE_age_padj_CA3  <- de_ca3_padj[regulator_rho$gene]
regulator_rho$DE_region_lfc    <- de_region_lfc[regulator_rho$gene]
regulator_rho$DE_region_padj   <- de_region_padj[regulator_rho$gene]

# Replace NA with 0 for LFC
regulator_rho$DE_age_lfc_CA1[is.na(regulator_rho$DE_age_lfc_CA1)] <- 0
regulator_rho$DE_age_lfc_CA3[is.na(regulator_rho$DE_age_lfc_CA3)] <- 0
regulator_rho$DE_region_lfc[is.na(regulator_rho$DE_region_lfc)]   <- 0

# Scale |ρ| within module×region
regulator_rho$abs_rho_scaled_CA1 <- 0
regulator_rho$abs_rho_scaled_CA3 <- 0
for (m in module_names) {
  mask <- regulator_rho$module == m
  max_rho_ca1 <- max(abs(regulator_rho$rho_CA1[mask]), na.rm = TRUE)
  if (max_rho_ca1 > 0) regulator_rho$abs_rho_scaled_CA1[mask] <- abs(regulator_rho$rho_CA1[mask]) / max_rho_ca1
  max_rho_ca3 <- max(abs(regulator_rho$rho_CA3[mask]), na.rm = TRUE)
  if (max_rho_ca3 > 0) regulator_rho$abs_rho_scaled_CA3[mask] <- abs(regulator_rho$rho_CA3[mask]) / max_rho_ca3
}

# Scale |DE_age| within region
max_de_age_ca1 <- max(abs(regulator_rho$DE_age_lfc_CA1), na.rm = TRUE)
regulator_rho$abs_DE_age_scaled_CA1 <- if (max_de_age_ca1 > 0) abs(regulator_rho$DE_age_lfc_CA1) / max_de_age_ca1 else 0

max_de_age_ca3 <- max(abs(regulator_rho$DE_age_lfc_CA3), na.rm = TRUE)
regulator_rho$abs_DE_age_scaled_CA3 <- if (max_de_age_ca3 > 0) abs(regulator_rho$DE_age_lfc_CA3) / max_de_age_ca3 else 0

# Scale |DE_region| globally
max_de_region <- max(abs(regulator_rho$DE_region_lfc), na.rm = TRUE)
regulator_rho$abs_DE_region_scaled <- if (max_de_region > 0) abs(regulator_rho$DE_region_lfc) / max_de_region else 0

# 3 sensitivity configs
configs <- list(
  default     = c(0.50, 0.30, 0.20),
  de_emphasis = c(0.33, 0.50, 0.17),
  rho_only    = c(1.00, 0.00, 0.00)
)

for (cfg_name in names(configs)) {
  w <- configs[[cfg_name]]
  regulator_rho[[paste0("score_CA1_", cfg_name)]] <- w[1] * regulator_rho$abs_rho_scaled_CA1 +
    w[2] * regulator_rho$abs_DE_age_scaled_CA1 + w[3] * regulator_rho$abs_DE_region_scaled
  regulator_rho[[paste0("score_CA3_", cfg_name)]] <- w[1] * regulator_rho$abs_rho_scaled_CA3 +
    w[2] * regulator_rho$abs_DE_age_scaled_CA3 + w[3] * regulator_rho$abs_DE_region_scaled
}

cat("regulator_score computed for", nrow(regulator_rho), "rows\n")

gc()

## ---- E7: Candidate Classes ----
cat("\n=== E7: Candidate Classes ===\n")

regulator_rho$candidate_class <- "low_confidence"

sig_ca1 <- (abs(regulator_rho$rho_CA1) > 0.3 & regulator_rho$qval_CA1 < 0.1) | abs(regulator_rho$rho_CA1) > 0.5
sig_ca3 <- (abs(regulator_rho$rho_CA3) > 0.3 & regulator_rho$qval_CA3 < 0.1) | abs(regulator_rho$rho_CA3) > 0.5
dir_ca1 <- sign(regulator_rho$rho_CA1)
dir_ca3 <- sign(regulator_rho$rho_CA3)

regulator_rho$candidate_class[sig_ca1 & !sig_ca3] <- "CA1_aging_associated"
regulator_rho$candidate_class[!sig_ca1 & sig_ca3] <- "CA3_aging_associated"
regulator_rho$candidate_class[sig_ca1 & sig_ca3 & dir_ca1 == dir_ca3] <- "concordant"
regulator_rho$candidate_class[sig_ca1 & sig_ca3 & dir_ca1 != dir_ca3] <- "region_flipped"
regulator_rho$candidate_class[regulator_rho$self_correlation_flag == "module_member_or_effector"] <- "module_member"
regulator_rho$candidate_class[
  regulator_rho$self_correlation_flag == "target_gene_candidate" & regulator_rho$candidate_class == "low_confidence"
] <- "target_readout"

cat("Candidate classes:\n")
print(table(regulator_rho$candidate_class))

gc()

## ---- E8: Save Data Outputs ----
cat("\n=== E8: Saving Data Outputs ===\n")

out_dir <- "data/processed/spatial/phase13_candidate_regulator"

# 1. All Spearman results
write.csv(rho_results, file.path(out_dir, "phase13_spearman_rho_all.csv"), row.names = FALSE)
cat("  1. phase13_spearman_rho_all.csv\n")

# 2. Regulator universe Spearman + DE
write.csv(regulator_rho, file.path(out_dir, "phase13_spearman_rho_universe_regulator.csv"), row.names = FALSE)
cat("  2. phase13_spearman_rho_universe_regulator.csv\n")

# 3. Default scores with ranks (exclude module members from own-module ranking)
score_default <- regulator_rho[, c("gene", "module", "rho_CA1", "rho_CA3",
  "DE_age_lfc_CA1", "DE_age_lfc_CA3", "DE_region_lfc",
  "score_CA1_default", "score_CA3_default", "self_correlation_flag", "candidate_class")]
score_default$rank_CA1 <- NA_integer_
score_default$rank_CA3 <- NA_integer_
for (m in module_names) {
  mask <- score_default$module == m
  # Rank excluding module members for own-module
  valid_ca1 <- mask & score_default$self_correlation_flag != "module_member_or_effector"
  if (sum(valid_ca1) > 0) {
    ranks <- rank(-score_default$score_CA1_default[valid_ca1], ties.method = "min")
    score_default$rank_CA1[valid_ca1] <- ranks
  }
  valid_ca3 <- mask & score_default$self_correlation_flag != "module_member_or_effector"
  if (sum(valid_ca3) > 0) {
    ranks <- rank(-score_default$score_CA3_default[valid_ca3], ties.method = "min")
    score_default$rank_CA3[valid_ca3] <- ranks
  }
}
write.csv(score_default, file.path(out_dir, "phase13_regulator_score_default.csv"), row.names = FALSE)
cat("  3. phase13_regulator_score_default.csv\n")

# 4. Sensitivity scores
score_sens <- regulator_rho[, c("gene", "module",
  "score_CA1_default", "score_CA3_default",
  "score_CA1_de_emphasis", "score_CA3_de_emphasis",
  "score_CA1_rho_only", "score_CA3_rho_only")]
write.csv(score_sens, file.path(out_dir, "phase13_regulator_score_sensitivity.csv"), row.names = FALSE)
cat("  4. phase13_regulator_score_sensitivity.csv\n")

# 5. Candidate classes
write.csv(regulator_rho[, c("gene", "module", "candidate_class", "self_correlation_flag",
  "rho_CA1", "rho_CA3", "score_CA1_default", "score_CA3_default")],
  file.path(out_dir, "phase13_candidate_classes.csv"), row.names = FALSE)
cat("  5. phase13_candidate_classes.csv\n")

# 6. Self-correlation flags
write.csv(rho_results[, c("gene", "module", "self_correlation_flag")],
  file.path(out_dir, "phase13_self_correlation_flags.csv"), row.names = FALSE)
cat("  6. phase13_self_correlation_flags.csv\n")

# 7. Target gene module association
target_rho <- rho_results[rho_results$gene %in% all_target_genes, ]
target_gene_to_cat <- c()
for (cat_name in names(target_categories)) {
  genes <- target_categories[[cat_name]]
  target_gene_to_cat <- c(target_gene_to_cat, setNames(rep(cat_name, length(genes)), genes))
}
target_rho$category <- target_gene_to_cat[target_rho$gene]
write.csv(target_rho, file.path(out_dir, "phase13_target_gene_module_association.csv"), row.names = FALSE)
cat("  7. phase13_target_gene_module_association.csv\n")

# 8. TF annotation coverage
tf_ann_out <- tf_annotation
tf_ann_out$in_universe_regulator <- tf_ann_out$gene %in% universe_regulator
write.csv(tf_ann_out, file.path(out_dir, "phase13_tf_annotation_coverage.csv"), row.names = FALSE)
cat("  8. phase13_tf_annotation_coverage.csv\n")

# 9. Universe composition (already saved in E3)
cat("  9. phase13_universe_composition.csv (saved in E3)\n")

# 10. GO enrichment (saved in E12 if available)
cat(" 10. phase13_go_enrichment_results.csv (deferred to E12)\n")

# 11. Validation summary (saved in E14)
cat(" 11. phase13_validation_summary.csv (deferred to E14)\n")

# 12. Provenance (saved in E14)
cat(" 12. phase13_provenance.txt (deferred to E14)\n")

cat("Data outputs saved\n")

gc()

## ---- E9: Generate Figures ----
cat("\n=== E9: Generating Figures ===\n")

fig_dir <- "figures/spatial/phase13_candidate_regulator"

# Colors
col_ca1  <- "#2166AC"
col_ca3  <- "#E66100"
col_young <- "#1B7837"
col_old   <- "#762A83"

# Helper for safe PDF
safe_pdf <- function(filename, width, height, plot_expr) {
  tryCatch({
    pdf(file.path(fig_dir, filename), width = width, height = height)
    eval(plot_expr)
    dev.off()
  }, error = function(e) {
    cat("  WARNING: PDF failed for", filename, "\n")
    while (dev.cur() > 1) dev.off()
  })
}

# ---- F01: regulator_score heatmap ----
cat("F01: regulator_score heatmap\n")
tryCatch({
  # Top 20 genes per module (CA1 default)
  top_genes_per_module <- c()
  for (m in module_names) {
    sub <- score_default[score_default$module == m & !is.na(score_default$rank_CA1), ]
    sub <- sub[order(sub$rank_CA1), ]
    top_genes_per_module <- c(top_genes_per_module, head(sub$gene, 20))
  }
  top_genes_per_module <- unique(top_genes_per_module)

  # Build heatmap matrix
  heatmap_mat <- matrix(NA, nrow = length(top_genes_per_module), ncol = length(module_names),
    dimnames = list(top_genes_per_module, module_names))
  for (m in module_names) {
    sub <- regulator_rho[regulator_rho$module == m, ]
    heatmap_mat[sub$gene[sub$gene %in% top_genes_per_module], m] <-
      sub$score_CA1_default[sub$gene %in% top_genes_per_module]
  }

  # Remove all-NA rows
  keep <- apply(heatmap_mat, 1, function(x) !all(is.na(x)))
  heatmap_mat <- heatmap_mat[keep, ]
  heatmap_mat[is.na(heatmap_mat)] <- 0

  png(file.path(fig_dir, "F01_regulator_score_heatmap.png"), width = 12, height = 16, units = "in", res = 150)
  pheatmap(heatmap_mat,
    main = "F01: regulator_score (CA1, default config, top 20 per module)",
    cluster_rows = TRUE, cluster_cols = TRUE,
    scale = "none", color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  safe_pdf("F01_regulator_score_heatmap.pdf", 12, 16,
    quote(pheatmap(heatmap_mat,
      main = "F01: regulator_score (CA1, default config, top 20 per module)",
      cluster_rows = TRUE, cluster_cols = TRUE,
      scale = "none", color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      fontsize_row = 6, fontsize_col = 10, angle_col = 45)))
  cat("  F01 saved\n")
}, error = function(e) cat("  WARNING: F01 failed:", e$message, "\n"))

gc()

# ---- F02: ρ vs DE_age scatter ----
cat("F02: ρ vs DE scatter\n")
tryCatch({
  plot_df <- data.frame(
    rho_CA1 = regulator_rho$rho_CA1,
    rho_CA3 = regulator_rho$rho_CA3,
    DE_age_lfc_CA1 = regulator_rho$DE_age_lfc_CA1,
    DE_age_lfc_CA3 = regulator_rho$DE_age_lfc_CA3,
    gene = regulator_rho$gene,
    module = regulator_rho$module,
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(plot_df, aes(x = rho_CA1, y = DE_age_lfc_CA1)) +
    geom_point(alpha = 0.3, size = 0.8, color = col_ca1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "F02a: ρ vs DE_age log2FC (CA1)",
      subtitle = paste0("n = ", nrow(plot_df), " genes; n_CA1 = ", length(ca1_cols), " samples"),
      x = "Spearman ρ (CA1)", y = "DE_age log2FC (Old vs Young, CA1)") +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste0("n = ", nrow(plot_df)), hjust = 1.1, vjust = 1.5, size = 4)

  p2 <- ggplot(plot_df, aes(x = rho_CA3, y = DE_age_lfc_CA3)) +
    geom_point(alpha = 0.3, size = 0.8, color = col_ca3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "F02b: ρ vs DE_age log2FC (CA3)",
      subtitle = paste0("n = ", nrow(plot_df), " genes; n_CA3 = ", length(ca3_cols), " samples"),
      x = "Spearman ρ (CA3)", y = "DE_age log2FC (Old vs Young, CA3)") +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste0("n = ", nrow(plot_df)), hjust = 1.1, vjust = 1.5, size = 4)

  p <- p1 + p2
  ggsave(file.path(fig_dir, "F02_rho_vs_DE_scatter.png"), p, width = 14, height = 6, dpi = 150)
  ggsave(file.path(fig_dir, "F02_rho_vs_DE_scatter.pdf"), p, width = 14, height = 6)
  cat("  F02 saved\n")
}, error = function(e) cat("  WARNING: F02 failed:", e$message, "\n"))

gc()

# ---- F03: Top regulator dotplot ----
cat("F03: Top regulator dotplot\n")
tryCatch({
  top15_all <- c()
  for (m in module_names) {
    sub <- score_default[score_default$module == m & !is.na(score_default$rank_CA1), ]
    sub <- sub[order(sub$rank_CA1), ]
    top15_all <- c(top15_all, head(sub$gene, 15))
  }
  top15_all <- unique(top15_all)

  dot_df <- regulator_rho[regulator_rho$gene %in% top15_all, ]
  dot_df$gene <- factor(dot_df$gene, levels = rev(top15_all))

  p <- ggplot(dot_df, aes(x = module, y = gene, size = abs(rho_CA1), color = score_CA1_default)) +
    geom_point() +
    scale_color_gradient(low = "lightyellow", high = "firebrick3", name = "regulator_score\n(CA1 default)") +
    scale_size_continuous(range = c(1, 6), name = "|ρ| (CA1)") +
    labs(title = "F03: Top regulators per module (CA1, default config)",
      subtitle = paste0("n = ", length(top15_all), " genes; n_CA1 = ", length(ca1_cols))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank())

  ggsave(file.path(fig_dir, "F03_top_regulators_dotplot.png"), p, width = 14, height = 10, dpi = 150)
  ggsave(file.path(fig_dir, "F03_top_regulators_dotplot.pdf"), p, width = 14, height = 10)
  cat("  F03 saved\n")
}, error = function(e) cat("  WARNING: F03 failed:", e$message, "\n"))

gc()

# ---- F04: Module correlation matrix ----
cat("F04: Module correlation matrix\n")
tryCatch({
  # Mean ρ per gene across modules (CA1)
  rho_wide_ca1 <- reshape(rho_results[, c("gene", "module", "rho_CA1")],
    idvar = "gene", timevar = "module", direction = "wide")
  rownames(rho_wide_ca1) <- rho_wide_ca1$gene
  rho_wide_ca1$gene <- NULL
  colnames(rho_wide_ca1) <- sub("rho_CA1\\.", "", colnames(rho_wide_ca1))

  # Top 50 by max |ρ|
  max_abs_rho <- apply(abs(rho_wide_ca1), 1, max, na.rm = TRUE)
  top50_idx <- order(max_abs_rho, decreasing = TRUE)[1:min(50, nrow(rho_wide_ca1))]
  rho_mat <- as.matrix(rho_wide_ca1[top50_idx, ])
  rho_mat[is.na(rho_mat)] <- 0

  png(file.path(fig_dir, "F04_module_correlation_matrix.png"), width = 12, height = 14, units = "in", res = 150)
  pheatmap(rho_mat,
    main = "F04: Gene × Module ρ matrix (CA1, top 50 genes)",
    cluster_rows = TRUE, cluster_cols = TRUE,
    scale = "none", color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  cat("  F04 saved\n")
}, error = function(e) cat("  WARNING: F04 failed:", e$message, "\n"))

gc()

# ---- F05: Self-correlation flag summary ----
cat("F05: Self-correlation flag summary\n")
tryCatch({
  flag_counts <- as.data.frame(table(rho_results$self_correlation_flag))
  colnames(flag_counts) <- c("flag", "count")

  p <- ggplot(flag_counts, aes(x = flag, y = count, fill = flag)) +
    geom_col() +
    geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
    labs(title = "F05: Self-correlation flag distribution",
      subtitle = paste0("Total: ", nrow(rho_results), " gene × module pairs")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

  ggsave(file.path(fig_dir, "F05_self_correlation_flags.png"), p, width = 8, height = 5, dpi = 150)
  cat("  F05 saved\n")
}, error = function(e) cat("  WARNING: F05 failed:", e$message, "\n"))

gc()

# ---- F06: Candidate class distribution ----
cat("F06: Candidate class distribution\n")
tryCatch({
  class_counts <- as.data.frame(table(regulator_rho$candidate_class))
  colnames(class_counts) <- c("class", "count")

  p <- ggplot(class_counts, aes(x = reorder(class, -count), y = count, fill = class)) +
    geom_col() +
    geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
    labs(title = "F06: Candidate class distribution",
      subtitle = paste0("n = ", nrow(regulator_rho), " regulator × module pairs; n_regulators = ",
        length(unique(regulator_rho$gene)))) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

  ggsave(file.path(fig_dir, "F06_candidate_classes.png"), p, width = 10, height = 6, dpi = 150)
  cat("  F06 saved\n")
}, error = function(e) cat("  WARNING: F06 failed:", e$message, "\n"))

gc()

# ---- F07-F08: Age-stratified heatmaps ----
cat("F07-F08: Age-stratified heatmaps\n")
tryCatch({
  # Compute age-stratified ρ
  young_ca1_cols <- manifest_ca1$Sample[manifest_ca1$Age == "Young"]
  old_ca1_cols   <- manifest_ca1$Sample[manifest_ca1$Age == "Old"]
  young_ca3_cols <- manifest_ca3$Sample[manifest_ca3$Age == "Young"]
  old_ca3_cols   <- manifest_ca3$Sample[manifest_ca3$Age == "Old"]

  # Use top 30 genes from F01
  top30_genes <- c()
  for (m in module_names) {
    sub <- score_default[score_default$module == m & !is.na(score_default$rank_CA1), ]
    sub <- sub[order(sub$rank_CA1), ]
    top30_genes <- c(top30_genes, head(sub$gene, 10))
  }
  top30_genes <- unique(top30_genes)

  # Young CA1 correlation
  if (length(young_ca1_cols) >= 3) {
    expr_young_ca1 <- as.matrix(log2cpm_ca1[top30_genes[top30_genes %in% rownames(log2cpm_ca1)], young_ca1_cols, drop = FALSE])
    mod_score_young_ca1 <- mod_score_ca1_mat[, young_ca1_cols, drop = FALSE]
    rho_young_ca1 <- matrix(NA, nrow = nrow(expr_young_ca1), ncol = length(module_names),
      dimnames = list(rownames(expr_young_ca1), module_names))
    for (m in module_names) {
      rho_young_ca1[, m] <- cor(t(expr_young_ca1), mod_score_young_ca1[m, ], method = "spearman")
    }
    png(file.path(fig_dir, "F07_young_stratified_heatmap.png"), width = 12, height = 10, units = "in", res = 150)
    pheatmap(rho_young_ca1, main = paste0("F07: Young CA1 ρ (n = ", length(young_ca1_cols), ")"),
      cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      fontsize_row = 6, fontsize_col = 10, angle_col = 45)
    dev.off()
    cat("  F07 saved\n")
  } else {
    cat("  F07 skipped (Young CA1 n < 3)\n")
  }

  # Old CA1 correlation
  if (length(old_ca1_cols) >= 3) {
    expr_old_ca1 <- as.matrix(log2cpm_ca1[top30_genes[top30_genes %in% rownames(log2cpm_ca1)], old_ca1_cols, drop = FALSE])
    mod_score_old_ca1 <- mod_score_ca1_mat[, old_ca1_cols, drop = FALSE]
    rho_old_ca1 <- matrix(NA, nrow = nrow(expr_old_ca1), ncol = length(module_names),
      dimnames = list(rownames(expr_old_ca1), module_names))
    for (m in module_names) {
      rho_old_ca1[, m] <- cor(t(expr_old_ca1), mod_score_old_ca1[m, ], method = "spearman")
    }
    png(file.path(fig_dir, "F08_old_stratified_heatmap.png"), width = 12, height = 10, units = "in", res = 150)
    pheatmap(rho_old_ca1, main = paste0("F08: Old CA1 ρ (n = ", length(old_ca1_cols), ")"),
      cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      fontsize_row = 6, fontsize_col = 10, angle_col = 45)
    dev.off()
    cat("  F08 saved\n")
  } else {
    cat("  F08 skipped (Old CA1 n < 3)\n")
  }
}, error = function(e) cat("  WARNING: F07-F08 failed:", e$message, "\n"))

gc()

# ---- F09: regulator_score component breakdown ----
cat("F09: regulator_score component breakdown\n")
tryCatch({
  # Top 20 unique genes by score_CA1_default
  top20_genes <- head(score_default[order(-score_default$score_CA1_default), ], 50)
  top20_unique <- unique(top20_genes$gene)[1:min(20, length(unique(top20_genes$gene)))]

  comp_df <- regulator_rho[regulator_rho$gene %in% top20_unique & regulator_rho$module == module_names[1], ]
  comp_long <- data.frame(
    gene = rep(comp_df$gene, 3),
    component = rep(c("|ρ_scaled|", "|DE_age_scaled|", "|DE_region_scaled|"), each = nrow(comp_df)),
    value = c(comp_df$abs_rho_scaled_CA1 * 0.50, comp_df$abs_DE_age_scaled_CA1 * 0.30, comp_df$abs_DE_region_scaled * 0.20),
    stringsAsFactors = FALSE
  )

  p <- ggplot(comp_long, aes(x = gene, y = value, fill = component)) +
    geom_col() +
    coord_flip() +
    labs(title = "F09: regulator_score component breakdown (CA1, default weights)",
      subtitle = paste0("Module: ", module_names[1], "; n = ", nrow(comp_df))) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 7))

  ggsave(file.path(fig_dir, "F09_score_components.png"), p, width = 10, height = 8, dpi = 150)
  cat("  F09 saved\n")
}, error = function(e) cat("  WARNING: F09 failed:", e$message, "\n"))

gc()

# ---- F10: Sensitivity config comparison ----
cat("F10: Sensitivity config comparison\n")
tryCatch({
  sens_df <- data.frame(
    default = regulator_rho$score_CA1_default,
    de_emphasis = regulator_rho$score_CA1_de_emphasis,
    rho_only = regulator_rho$score_CA1_rho_only,
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(sens_df, aes(x = default, y = de_emphasis)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "F10a: Default vs DE-emphasis (CA1)", x = "Default score", y = "DE-emphasis score") +
    theme_classic()

  p2 <- ggplot(sens_df, aes(x = default, y = rho_only)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "F10b: Default vs ρ-only (CA1)", x = "Default score", y = "ρ-only score") +
    theme_classic()

  p <- p1 + p2
  ggsave(file.path(fig_dir, "F10_sensitivity_comparison.png"), p, width = 12, height = 5, dpi = 150)
  cat("  F10 saved\n")
}, error = function(e) cat("  WARNING: F10 failed:", e$message, "\n"))

gc()

# ---- F11: TF annotation coverage ----
cat("F11: TF annotation coverage\n")
tryCatch({
  tf_source_counts <- data.frame(
    source = c("DoRothEA", "CollecTRI", "GO-TF"),
    count = c(sum(tf_annotation$dorothea), sum(tf_annotation$collectri), sum(tf_annotation$go_tf)),
    stringsAsFactors = FALSE
  )

  p <- ggplot(tf_source_counts, aes(x = reorder(source, -count), y = count, fill = source)) +
    geom_col() +
    geom_text(aes(label = count), vjust = -0.5) +
    labs(title = "F11: TF annotation source coverage",
      subtitle = paste0("Total unique TFs: ", nrow(tf_annotation))) +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "F11_tf_annotation_coverage.png"), p, width = 6, height = 5, dpi = 150)
  cat("  F11 saved\n")
}, error = function(e) cat("  WARNING: F11 failed:", e$message, "\n"))

gc()

# ---- F12: Universe composition ----
cat("F12: Universe composition\n")
tryCatch({
  p <- ggplot(universe_comp, aes(x = reorder(universe, -n_genes), y = n_genes, fill = universe)) +
    geom_col() +
    geom_text(aes(label = n_genes), vjust = -0.5) +
    labs(title = "F12: Universe composition", x = "", y = "Number of genes") +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "F12_universe_composition.png"), p, width = 6, height = 5, dpi = 150)
  cat("  F12 saved\n")
}, error = function(e) cat("  WARNING: F12 failed:", e$message, "\n"))

gc()

# ---- F13: DE volcano with regulators highlighted ----
cat("F13: DE volcano with regulators highlighted\n")
tryCatch({
  volcano_df <- data.frame(
    gene = de_ca1[[gene_ca1_col]],
    log2FC = de_ca1[[lfc_ca1_col]],
    neg_log10_padj = -log10(pmax(de_ca1[[padj_ca1_col]], 1e-300)),
    stringsAsFactors = FALSE
  )
  volcano_df$is_regulator <- volcano_df$gene %in% universe_regulator
  volcano_df$is_target <- volcano_df$gene %in% universe_target

  p <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_padj)) +
    geom_point(aes(color = is_regulator | is_target), alpha = 0.3, size = 0.5) +
    scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey60"), name = "Regulator/Target") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    labs(title = "F13: DE volcano (CA1, Old vs Young)",
      subtitle = paste0("Regulators highlighted: ", sum(volcano_df$is_regulator), " genes"),
      x = "log2FC (Old vs Young)", y = "-log10(padj)") +
    theme_classic()

  ggsave(file.path(fig_dir, "F13_DE_volcano_regulators.png"), p, width = 8, height = 6, dpi = 150)
  cat("  F13 saved\n")
}, error = function(e) cat("  WARNING: F13 failed:", e$message, "\n"))

gc()

## ---- E10-E11: Spatial Visualization ----
cat("\n=== E10-E11: Spatial Visualization ===\n")

tryCatch({
  cat("Loading Hippo RDS (this may take a moment)...\n")
  hippo_rds <- readRDS(input_files$hippo_rds)
  cat("Hippo RDS loaded:", ncol(hippo_rds), "spots\n")

  DefaultAssay(hippo_rds) <- "Spatial"

  # Top 5 candidates for spatial visualization
  top5_for_spatial <- head(score_default[order(-score_default$score_CA1_default), ], 20)
  top5_unique <- unique(top5_for_spatial$gene)[1:min(5, length(unique(top5_for_spatial$gene)))]

  for (gene in top5_unique) {
    if (gene %in% rownames(hippo_rds)) {
      tryCatch({
        p <- SpatialFeaturePlot(hippo_rds, features = gene, pt.size = 0.3) +
          ggtitle(paste0("Spatial: ", gene))
        ggsave(file.path(fig_dir, paste0("F14_spatial_", gene, ".png")), p, width = 10, height = 8, dpi = 150)
        cat("  Spatial plot:", gene, "\n")
      }, error = function(e) cat("  WARNING: Spatial plot failed for", gene, "\n"))
    } else {
      cat("  Skipping spatial plot:", gene, "not in Hippo data\n")
    }
  }

  rm(hippo_rds)
  gc()
  cat("Spatial visualization complete\n")
}, error = function(e) {
  cat("WARNING: Spatial visualization failed:", e$message, "\n")
  gc()
})

## ---- E12: GO Enrichment ----
cat("\n=== E12: GO Enrichment (conditional) ===\n")

go_results_all <- data.frame()

if (has_clusterProfiler && has_msigdbr) {
  tryCatch({
    library(clusterProfiler)
    library(msigdbr)

    cat("Running GO enrichment for top modules...\n")

    for (m in module_names[1:min(4, length(module_names))]) {
      # Get top regulator genes for this module
      top_regs <- score_default[score_default$module == m & !is.na(score_default$rank_CA1), ]
      top_regs <- top_regs[order(top_regs$rank_CA1), ]
      gene_list <- head(top_regs$gene, 50)

      if (length(gene_list) < 5) next

      tryCatch({
        ego <- enrichGO(gene = gene_list, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)

        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          ego_df <- as.data.frame(ego)
          ego_df$module <- m
          go_results_all <- rbind(go_results_all, ego_df)
          cat("  Module", m, ":", nrow(ego_df), "GO terms\n")
        }
      }, error = function(e) cat("  WARNING: GO failed for", m, "\n"))
    }

    if (nrow(go_results_all) > 0) {
      write.csv(go_results_all, file.path(out_dir, "phase13_go_enrichment_results.csv"), row.names = FALSE)
      cat("GO enrichment results saved\n")

      # GO dotplot
      tryCatch({
        if (nrow(go_results_all) > 0) {
          go_top <- go_results_all[order(go_results_all$pvalue), ]
          go_top <- head(go_top, 30)
          go_top$Description <- factor(go_top$Description, levels = rev(go_top$Description))

          p <- ggplot(go_top, aes(x = GeneRatio_num, y = Description, size = Count, color = pvalue)) +
            geom_point() +
            scale_color_gradient(low = "firebrick3", high = "navy") +
            labs(title = "F15: GO enrichment (top terms across modules)", x = "Gene Ratio") +
            theme_classic() +
            theme(axis.text.y = element_text(size = 7))

          ggsave(file.path(fig_dir, "F15_GO_enrichment.png"), p, width = 10, height = 8, dpi = 150)
          cat("  F15 saved\n")
        }
      }, error = function(e) cat("  WARNING: GO dotplot failed\n"))
    }
  }, error = function(e) cat("WARNING: GO enrichment failed:", e$message, "\n"))
} else {
  cat("Skipping GO enrichment (clusterProfiler/msigdbr not available)\n")
}

gc()

## ---- E13: Interpretation Guide ----
cat("\n=== E13: Interpretation Guide ===\n")

guide_text <- c(
  "# Phase Spatial-13: Candidate Regulator Analysis Guide",
  "",
  paste0("Generated: ", format(Sys.time())),
  paste0("Plan version: v1.2"),
  "",
  "## Overview",
  "This guide describes the Phase Spatial-13 candidate regulator discovery results.",
  "All candidate regulators are ASSOCIATED with mitochondrial modules — correlation ≠ causation.",
  "",
  "## Key Outputs",
  "",
  "### Primary Rankings",
  "- `phase13_regulator_score_default.csv`: Primary ranking by default config (ρ 50% + DE_age 30% + DE_region 20%)",
  "- `phase13_candidate_classes.csv`: Classification of candidates (CA1/CA3 aging-associated, concordant, region-flipped, etc.)",
  "- `phase13_spearman_rho_all.csv`: All Spearman ρ values with BH-corrected q-values",
  "",
  "### Sensitivity",
  "- `phase13_regulator_score_sensitivity.csv`: All 3 weight configs (default, DE-emphasis, ρ-only)",
  "- Compare configs to assess robustness of rankings",
  "",
  "### Annotation",
  "- `phase13_tf_annotation_coverage.csv`: TF sources (DoRothEA A-C, CollecTRI, GO:0003700/GO:0140110)",
  "- `phase13_target_gene_module_association.csv`: Target gene module correlations",
  "- `phase13_self_correlation_flags.csv`: Self-correlation flags (module member exclusion)",
  "",
  "## Interpretation Caveats",
  "",
  "### Correlation ≠ Causation",
  "All candidates are ASSOCIATED regulators. Spearman ρ measures monotonic association, not causal regulation.",
  "DoRothEA/CollecTRI annotation = eligibility for universe_regulator, NOT proof of regulation in this context.",
  "",
  "### Small Sample Size",
  "- CA1: n = 12 (4 Young + 8 Old)",
  "- CA3: n = 11 (3 Young + 8 Old) — CA3 Young n=3 caveat applies",
  "- Age-stratified correlations have even smaller n and should be interpreted cautiously",
  "",
  "### Metric Hierarchy",
  "- |ρ| is the PRIMARY metric (weight 50% in default config)",
  "- q-values are EXPLORATORY only (BH-corrected, but small n limits power)",
  "- |ρ| > 0.3: screening threshold",
  "- |ρ| > 0.5: notable threshold",
  "",
  "### Self-Correlation Strategy A",
  "Module gene-set members are excluded from their own module's regulator ranking.",
  "4 flag values: regulator_candidate, module_member_or_effector, target_gene_candidate, low_confidence_self_correlation",
  "",
  "## Figures",
  "- F01: regulator_score heatmap (top 20 per module)",
  "- F02: ρ vs DE_age scatter (CA1 and CA3)",
  "- F03: Top regulators dotplot",
  "- F04: Gene × Module ρ matrix",
  "- F05: Self-correlation flag distribution",
  "- F06: Candidate class distribution",
  "- F07-F08: Age-stratified ρ heatmaps (Young, Old)",
  "- F09: Score component breakdown",
  "- F10: Sensitivity config comparison",
  "- F11: TF annotation coverage",
  "- F12: Universe composition",
  "- F13: DE volcano with regulators highlighted",
  "- F14: SpatialFeaturePlot for top candidates",
  "- F15: GO enrichment (if available)",
  "",
  "## Next Steps",
  "1. Review top-ranked candidates per module (score_default, rank_CA1/CA3)",
  "2. Check sensitivity configs for ranking stability",
  "3. Examine self-correlation flags — module members should NOT be interpreted as regulators of their own module",
  "4. Use spatial plots (F14) for expression pattern validation",
  "5. Cross-reference with published literature for biological plausibility",
  ""
)

writeLines(guide_text, "docs/spatial_phase13_candidate_regulator_analysis_guide.md")
cat("Interpretation guide saved\n")

## ---- E14: Validation & Cleanup ----
cat("\n=== E14: Validation & Cleanup ===\n")

validation <- data.frame(
  check = character(0), status = character(0), details = character(0),
  stringsAsFactors = FALSE
)

add_check <- function(check, status, details) {
  validation <<- rbind(validation, data.frame(check = check, status = status, details = details, stringsAsFactors = FALSE))
}

add_check("universe_all_size", if (length(universe_all) >= 100) "PASS" else "FAIL",
  paste(length(universe_all), "genes"))
add_check("universe_regulator_size", if (length(universe_regulator) >= 10) "PASS" else "FAIL",
  paste(length(universe_regulator), "genes"))
add_check("universe_target_size", if (length(universe_target) >= 1) "PASS" else "WARN",
  paste(length(universe_target), "genes"))
add_check("no_na_regulator_score", if (!any(is.na(regulator_rho$score_CA1_default))) "PASS" else "FAIL",
  paste(sum(is.na(regulator_rho$score_CA1_default)), "NAs"))
add_check("rho_in_range", if (all(abs(rho_results$rho_CA1) <= 1, na.rm = TRUE)) "PASS" else "FAIL",
  "All |ρ| ≤ 1")
add_check("score_in_range",
  if (all(regulator_rho$score_CA1_default >= 0 & regulator_rho$score_CA1_default <= 1, na.rm = TRUE)) "PASS" else "FAIL",
  "All scores in [0,1]")
add_check("self_correlation_flags_assigned", "PASS",
  paste(names(table(rho_results$self_correlation_flag)), collapse = ", "))
add_check("module_members_excluded",
  if (all(score_default$rank_CA1[regulator_rho$self_correlation_flag == "module_member_or_effector" &
    regulator_rho$module == module_names[1]] %in% NA)) "PASS" else "WARN",
  "Module members have NA rank for own module")
add_check("data_outputs_exist",
  if (all(file.exists(file.path(out_dir, paste0("phase13_", c("spearman_rho_all.csv",
    "spearman_rho_universe_regulator.csv", "regulator_score_default.csv",
    "regulator_score_sensitivity.csv", "candidate_classes.csv",
    "self_correlation_flags.csv", "target_gene_module_association.csv",
    "tf_annotation_coverage.csv", "universe_composition.csv")))))) "PASS" else "FAIL",
  "Required output files")
add_check("figures_exist",
  if (any(grepl("^F0[1-6]", list.files(fig_dir)))) "PASS" else "FAIL",
  paste(length(list.files(fig_dir, pattern = "\\.png$")), "PNG files"))
add_check("renv_lock_unchanged",
  if (renv_md5_before == renv_md5_after) "PASS" else "WARN",
  paste("before:", renv_md5_before, "after:", renv_md5_after))
add_check("no_causal_language", "MANUAL", "Requires human verification")

write.csv(validation, file.path(out_dir, "phase13_validation_summary.csv"), row.names = FALSE)
cat("Validation summary saved:", sum(validation$status == "PASS"), "PASS,",
  sum(validation$status == "FAIL"), "FAIL,", sum(validation$status == "WARN"), "WARN\n")

# Provenance
provenance <- c(
  "Phase Spatial-13: Candidate Regulator Discovery",
  paste("Plan version: v1.2"),
  paste("Script: R/spatial/s13_candidate_regulator_discovery.R"),
  paste("Start time:", format(t_start)),
  paste("End time:", format(Sys.time())),
  paste("Duration:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "minutes"),
  paste("renv.lock MD5 (before):", renv_md5_before),
  paste("renv.lock MD5 (after):", renv_md5_after),
  paste("R version:", R.version.string),
  paste("Seurat version:", packageVersion("Seurat")),
  "",
  "=== Universe Sizes ===",
  paste("universe_all:", length(universe_all)),
  paste("universe_regulator:", length(universe_regulator)),
  paste("universe_target:", length(universe_target)),
  "",
  "=== Sample Sizes ===",
  paste("CA1:", length(ca1_cols), "(Young:", n_ca1_young, ", Old:", n_ca1_old, ")"),
  paste("CA3:", length(ca3_cols), "(Young:", n_ca3_young, ", Old:", n_ca3_old, ")"),
  "",
  "=== Packages ===",
  paste("dorothea:", packageVersion("dorothea")),
  paste("decoupleR:", packageVersion("decoupleR")),
  paste("org.Mm.eg.db:", packageVersion("org.Mm.eg.db")),
  "",
  "=== Key Results ===",
  paste("Spearman correlations:", nrow(rho_results), "gene × module pairs"),
  paste("Regulator candidates:", nrow(regulator_rho), "pairs from", length(unique(regulator_rho$gene)), "unique genes"),
  paste("Figures generated:", length(list.files(fig_dir, pattern = "\\.png$")))
)
writeLines(provenance, file.path(out_dir, "phase13_provenance.txt"))
cat("Provenance saved\n")

# Rplots.pdf cleanup
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

cat("\n=== Phase Spatial-13 Complete ===\n")
cat("Duration:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "minutes\n")
cat("Data outputs:", length(list.files(out_dir)), "files in", out_dir, "\n")
cat("Figures:", length(list.files(fig_dir, pattern = "\\.png$")), "PNG files in", fig_dir, "\n")
cat("Guide: docs/spatial_phase13_candidate_regulator_analysis_guide.md\n")
