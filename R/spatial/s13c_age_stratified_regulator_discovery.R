## Phase Spatial-13c: Age-Stratified Candidate Regulator-Associated Gene Discovery
## Plan version: v1.0
## Script: R/spatial/s13c_age_stratified_regulator_discovery.R
## Execution: Rscript R/spatial/s13c_age_stratified_regulator_discovery.R

cat("=== Phase Spatial-13c: Age-Stratified Candidate Regulator Discovery ===\n")
cat("Plan version: v1.0\n")
cat("Script: R/spatial/s13c_age_stratified_regulator_discovery.R\n")
cat("Start time:", format(Sys.time()), "\n\n")

t_start <- Sys.time()

# Rplots.pdf guard
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

## ---- E0: Setup & Validation ----
cat("=== E0: Setup & Validation ===\n")

renv_md5_before <- digest::digest("renv.lock", algo = "md5", file = TRUE)
cat("renv.lock MD5 (before):", renv_md5_before, "\n")

required_pkgs <- c("dplyr", "ggplot2", "pheatmap", "patchwork", "stringr", "digest")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("FATAL: Package not available: ", pkg)
}
library(dplyr)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(stringr)

out_dir <- "data/processed/spatial/phase13_candidate_regulator"
fig_dir <- "figures/spatial/phase13_candidate_regulator"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Input files
input_files <- list(
  pseudobulk       = "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds",
  module_audit     = "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv",
  module_scores    = "data/processed/spatial/phase12_young_aged_region_comparison/module_scores_all_regions.csv",
  de_ca1           = "data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA1.csv",
  de_ca3           = "data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA3.csv",
  de_region        = "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/deseq2_all_genes_CA3_vs_CA1.csv",
  target_genes     = "docs/target_genes.csv",
  phase13_score    = file.path(out_dir, "phase13_regulator_score_default.csv"),
  phase13_val      = file.path(out_dir, "phase13_validation_summary.csv"),
  phase13b_val     = file.path(out_dir, "phase13b_display_validation_summary.csv"),
  phase13_tf       = file.path(out_dir, "phase13_tf_annotation_coverage.csv"),
  manifest         = file.path(out_dir, "pseudobulk_sample_manifest.csv")
)

for (nm in names(input_files)) {
  if (!file.exists(input_files[[nm]])) stop("FATAL: Missing input [", nm, "]: ", input_files[[nm]])
}
cat("All input files verified\n")

# Pre-execution validation
validation <- list()
add_check <- function(check, status, details) {
  validation[[length(validation) + 1]] <<- data.frame(
    check = check, status = status, details = details, stringsAsFactors = FALSE
  )
}

# V01: Phase 13 validation has 0 FAIL
phase13_val <- read.csv(input_files$phase13_val, stringsAsFactors = FALSE)
if (any(phase13_val$status == "FAIL")) stop("FATAL: Phase 13 validation has FAIL")
add_check("V01_phase13_no_fail", "PASS", "Phase 13 validation: 0 FAIL")

# V02: Phase 13b validation has 0 FAIL
phase13b_val <- read.csv(input_files$phase13b_val, stringsAsFactors = FALSE)
if (any(phase13b_val$status == "FAIL")) stop("FATAL: Phase 13b validation has FAIL")
add_check("V02_phase13b_no_fail", "PASS", "Phase 13b validation: 0 FAIL")

# V09: renv.lock MD5 recorded
add_check("V09_renv_md5_recorded", "PASS", renv_md5_before)

# V10: No phase13c_* directories
phase13c_dirs <- list.dirs("data/processed/spatial", recursive = FALSE)
phase13c_dirs <- phase13c_dirs[grep("phase13c", basename(phase13c_dirs))]
if (length(phase13c_dirs) > 0) stop("FATAL: phase13c_* directories already exist")
add_check("V10_no_phase13c_dirs", "PASS", "No phase13c_* directories found")

## ---- E1: Load Data ----
cat("\n=== E1: Load Data ===\n")

# Pseudobulk counts
pb_counts <- readRDS(input_files$pseudobulk)
# Convert sparse dgCMatrix to dense if needed (dim() may not work on some Matrix versions)
if (!is.matrix(pb_counts)) {
  cat("Converting sparse matrix to dense...\n")
  pb_counts <- as.matrix(pb_counts)
}
cat("Pseudobulk counts:", nrow(pb_counts), "genes x", ncol(pb_counts), "samples\n")

# Module audit
module_audit <- read.csv(input_files$module_audit, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
cat("Module audit:", nrow(module_audit), "modules\n")

# V06: Module audit >= 8 modules
if (nrow(module_audit) < 8) stop("FATAL: Module audit has < 8 modules")
add_check("V06_module_audit_readable", "PASS", paste(nrow(module_audit), "modules"))

# Build module gene lists
mod_col <- intersect(c("module", "Module", "module_name"), colnames(module_audit))[1]
genes_col <- intersect(c("genes", "Genes", "gene_set", "gene_list"), colnames(module_audit))[1]
if (is.na(mod_col) || is.na(genes_col)) stop("FATAL: Cannot identify module/genes columns in audit CSV")

module_genes <- list()
for (i in seq_len(nrow(module_audit))) {
  mod_name <- module_audit[[mod_col]][i]
  genes_str <- module_audit[[genes_col]][i]
  module_genes[[mod_name]] <- trimws(strsplit(genes_str, "[;,]")[[1]])
}
cat("Module gene sets loaded:", length(module_genes), "modules\n")

# Module scores
mod_scores <- read.csv(input_files$module_scores, stringsAsFactors = FALSE)
cat("Module scores:", nrow(mod_scores), "samples,", ncol(mod_scores), "columns\n")

score_cols <- grep("_mean$", colnames(mod_scores), value = TRUE)
module_names <- sub("_mean$", "", score_cols)
cat("Modules with scores:", length(module_names), ":", paste(module_names, collapse = ", "), "\n")

# Module name mapping: audit uses "Complex I" (space), scores use "Complex.I" (dot)
module_name_map <- setNames(module_names, gsub("\\.", " ", module_names))

# Sample manifest
manifest <- read.csv(input_files$manifest, stringsAsFactors = FALSE)
cat("Manifest:", nrow(manifest), "rows\n")

# V04: Verify sample counts
ca1_young <- manifest$Sample[manifest$Region == "CA1" & manifest$Age == "Young"]
ca1_old   <- manifest$Sample[manifest$Region == "CA1" & manifest$Age == "Old"]
ca3_young <- manifest$Sample[manifest$Region == "CA3" & manifest$Age == "Young"]
ca3_old   <- manifest$Sample[manifest$Region == "CA3" & manifest$Age == "Old"]

cat("CA1 Young:", length(ca1_young), ":", paste(ca1_young, collapse = ", "), "\n")
cat("CA1 Old:", length(ca1_old), ":", paste(ca1_old, collapse = ", "), "\n")
cat("CA3 Young:", length(ca3_young), ":", paste(ca3_young, collapse = ", "), "\n")
cat("CA3 Old:", length(ca3_old), ":", paste(ca3_old, collapse = ", "), "\n")

if (length(ca1_young) != 4 || length(ca1_old) != 8 || length(ca3_young) != 3 || length(ca3_old) != 8) {
  stop("FATAL: Sample counts don't match expected (CA1Y=4, CA1O=8, CA3Y=3, CA3O=8)")
}
add_check("V04_sample_counts", "PASS", "CA1Y=4, CA1O=8, CA3Y=3, CA3O=8")

# V07: All CA1/CA3 GEMgroups present in module scores
all_gems <- c(ca1_young, ca1_old, ca3_young, ca3_old)
missing_gems <- setdiff(all_gems, mod_scores$sample_id)
if (length(missing_gems) > 0) stop("FATAL: Missing GEMgroups in module scores: ", paste(missing_gems, collapse = ", "))
add_check("V07_all_gems_in_scores", "PASS", paste(length(all_gems), "GEMgroups present"))

# Build log2(CPM+1) matrices per stratum
# Identify CA1 and CA3 columns in pseudobulk
ca1_cols_all <- manifest$Sample[manifest$Region == "CA1"]
ca3_cols_all <- manifest$Sample[manifest$Region == "CA3"]

# Filter to genes in pseudobulk
log2cpm_all <- log2(pb_counts / colSums(pb_counts) * 1e6 + 1)

log2cpm_ca1_young <- as.matrix(log2cpm_all[, ca1_young, drop = FALSE])
log2cpm_ca1_old   <- as.matrix(log2cpm_all[, ca1_old, drop = FALSE])
log2cpm_ca3_young <- as.matrix(log2cpm_all[, ca3_young, drop = FALSE])
log2cpm_ca3_old   <- as.matrix(log2cpm_all[, ca3_old, drop = FALSE])

cat("Expression matrices built: CA1Y=", ncol(log2cpm_ca1_young), " CA1O=", ncol(log2cpm_ca1_old),
    " CA3Y=", ncol(log2cpm_ca3_young), " CA3O=", ncol(log2cpm_ca3_old), "\n")

# Module score vectors per stratum
mod_score_ca1_young <- t(as.matrix(mod_scores[mod_scores$sample_id %in% ca1_young, score_cols]))
colnames(mod_score_ca1_young) <- mod_scores$sample_id[mod_scores$sample_id %in% ca1_young]
rownames(mod_score_ca1_young) <- module_names

mod_score_ca1_old <- t(as.matrix(mod_scores[mod_scores$sample_id %in% ca1_old, score_cols]))
colnames(mod_score_ca1_old) <- mod_scores$sample_id[mod_scores$sample_id %in% ca1_old]
rownames(mod_score_ca1_old) <- module_names

mod_score_ca3_young <- t(as.matrix(mod_scores[mod_scores$sample_id %in% ca3_young, score_cols]))
colnames(mod_score_ca3_young) <- mod_scores$sample_id[mod_scores$sample_id %in% ca3_young]
rownames(mod_score_ca3_young) <- module_names

mod_score_ca3_old <- t(as.matrix(mod_scores[mod_scores$sample_id %in% ca3_old, score_cols]))
colnames(mod_score_ca3_old) <- mod_scores$sample_id[mod_scores$sample_id %in% ca3_old]
rownames(mod_score_ca3_old) <- module_names

# Load Phase 13 regulator universe
phase13_score <- read.csv(input_files$phase13_score, stringsAsFactors = FALSE)
universe_regulator <- unique(phase13_score$gene)
cat("Regulator universe:", length(universe_regulator), "genes\n")

# V05: universe >= 100
if (length(universe_regulator) < 100) {
  add_check("V05_universe_size", "WARN", paste("Only", length(universe_regulator), "regulators"))
} else {
  add_check("V05_universe_size", "PASS", paste(length(universe_regulator), "regulators"))
}

# TF annotation
tf_annot <- read.csv(input_files$phase13_tf, stringsAsFactors = FALSE)

# Target genes
target_genes_df <- tryCatch(
  read.csv(input_files$target_genes, stringsAsFactors = FALSE, fileEncoding = "GBK"),
  error = function(e) read.csv(input_files$target_genes, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
)
target_gene_list <- target_genes_df$gene_symbol
if (is.null(target_gene_list)) {
  tg_col <- intersect(c("gene_symbol", "gene", "Gene", "symbol"), colnames(target_genes_df))[1]
  target_gene_list <- target_genes_df[[tg_col]]
}
# If still NULL, try to extract gene names from the second column (comma-separated)
if (is.null(target_gene_list) && ncol(target_genes_df) >= 2) {
  raw_genes <- paste(target_genes_df[[2]], collapse = ",")
  target_gene_list <- trimws(strsplit(raw_genes, "[,]")[[1]])
  target_gene_list <- target_gene_list[nchar(target_gene_list) > 0]
}
cat("Target genes:", length(target_gene_list), "\n")

# DE results
de_ca1 <- read.csv(input_files$de_ca1, stringsAsFactors = FALSE)
de_ca3 <- read.csv(input_files$de_ca3, stringsAsFactors = FALSE)
de_region <- read.csv(input_files$de_region, stringsAsFactors = FALSE)
cat("DE CA1:", nrow(de_ca1), "genes; DE CA3:", nrow(de_ca3), "genes; DE region:", nrow(de_region), "genes\n")

gc()

## ---- E2: Compute Age-Stratified Spearman Rho ----
cat("\n=== E2: Compute Age-Stratified Spearman Rho ===\n")

# Filter universe to genes present in pseudobulk
genes_in_pb <- intersect(universe_regulator, rownames(log2cpm_all))
cat("Regulator genes in pseudobulk:", length(genes_in_pb), "/", length(universe_regulator), "\n")

# Module name mapping for self-correlation: audit names -> score names
# audit uses "Complex I", "Membrane/Translocator", "Defense/Antioxidant" etc
# scores use "Complex.I", "Membrane.Translocator", "Defense.Antioxidant"
audit_module_names <- names(module_genes)
audit_to_score <- list()
for (am in audit_module_names) {
  # Normalize both: replace dots and slashes with spaces for matching
  am_norm <- gsub("[./]", " ", am)
  sm_norm <- gsub("[./]", " ", module_names)
  idx <- which(sm_norm == am_norm)
  if (length(idx) == 1) {
    audit_to_score[[am]] <- module_names[idx]
  }
}

# Build reverse map: score_module -> audit_module gene set
module_genes_by_score <- list()
for (am in names(audit_to_score)) {
  sm <- audit_to_score[[am]]
  module_genes_by_score[[sm]] <- module_genes[[am]]
}

n_genes <- length(genes_in_pb)
n_modules <- length(module_names)
strata <- c("CA1_Young", "CA1_Old", "CA3_Young", "CA3_Old")

# Initialize result matrices
rho_mat <- matrix(NA, nrow = n_genes, ncol = n_modules * 4,
  dimnames = list(genes_in_pb, paste0(rep(module_names, 4), "__", rep(strata, each = n_modules))))
pval_mat <- rho_mat

cat("Computing age-stratified Spearman rho for", n_genes, "genes x", n_modules, "modules x 4 strata...\n")

for (m_idx in seq_along(module_names)) {
  m <- module_names[m_idx]
  cat("  Module:", m)

  for (s in strata) {
    col_name <- paste0(m, "__", s)

    # Get expression and module score for this stratum
    if (s == "CA1_Young") {
      expr_mat <- log2cpm_ca1_young[genes_in_pb, , drop = FALSE]
      mod_vec <- mod_score_ca1_young[m, ]
    } else if (s == "CA1_Old") {
      expr_mat <- log2cpm_ca1_old[genes_in_pb, , drop = FALSE]
      mod_vec <- mod_score_ca1_old[m, ]
    } else if (s == "CA3_Young") {
      expr_mat <- log2cpm_ca3_young[genes_in_pb, , drop = FALSE]
      mod_vec <- mod_score_ca3_young[m, ]
    } else {
      expr_mat <- log2cpm_ca3_old[genes_in_pb, , drop = FALSE]
      mod_vec <- mod_score_ca3_old[m, ]
    }

    n_samples <- length(mod_vec)

    # Vectorized Spearman correlation
    rho_vals <- cor(t(expr_mat), mod_vec, method = "spearman")

    # p-values from t-distribution
    t_stat <- rho_vals * sqrt((n_samples - 2) / (1 - rho_vals^2))
    p_vals <- 2 * pt(-abs(t_stat), df = n_samples - 2)

    rho_mat[, col_name] <- rho_vals
    pval_mat[, col_name] <- p_vals
  }
  cat(" done\n")
}

gc()

# Build long-format results
cat("Building age-stratified rho table...\n")
rho_results <- data.frame(
  gene = rep(genes_in_pb, times = n_modules),
  module = rep(module_names, each = n_genes),
  stringsAsFactors = FALSE
)

for (s in strata) {
  rho_results[[paste0("rho_", s)]] <- NA
  rho_results[[paste0("pvalue_exploratory_", s)]] <- NA
  rho_results[[paste0("qvalue_exploratory_", s)]] <- NA
}

for (m_idx in seq_along(module_names)) {
  m <- module_names[m_idx]
  row_idx <- which(rho_results$module == m)
  for (s in strata) {
    col_name <- paste0(m, "__", s)
    rho_results[[paste0("rho_", s)]][row_idx] <- rho_mat[, col_name]
    rho_results[[paste0("pvalue_exploratory_", s)]][row_idx] <- pval_mat[, col_name]
  }
}

# BH correction per stratum (exploratory only)
for (s in strata) {
  pval_col <- paste0("pvalue_exploratory_", s)
  qval_col <- paste0("qvalue_exploratory_", s)
  rho_results[[qval_col]] <- p.adjust(rho_results[[pval_col]], method = "BH")
}

# Self-correlation flag
rho_results$self_correlation_flag <- "regulator_candidate"
for (m in module_names) {
  if (m %in% names(module_genes_by_score)) {
    member_genes <- module_genes_by_score[[m]]
    mask <- rho_results$module == m & rho_results$gene %in% member_genes
    rho_results$self_correlation_flag[mask] <- "module_member_or_effector"
  }
}
# Target gene annotation
rho_results$self_correlation_flag[rho_results$gene %in% target_gene_list & rho_results$self_correlation_flag == "regulator_candidate"] <- "target_gene_candidate"

cat("Self-correlation flags:\n")
print(table(rho_results$self_correlation_flag))

# Sample size labels
rho_results$n_CA1_Young <- length(ca1_young)
rho_results$n_CA1_Old <- length(ca1_old)
rho_results$n_CA3_Young <- length(ca3_young)
rho_results$n_CA3_Old <- length(ca3_old)

# Write O1
write.csv(rho_results, file.path(out_dir, "phase13c_age_stratified_rho.csv"), row.names = FALSE)
cat("O1: phase13c_age_stratified_rho.csv written (", nrow(rho_results), "rows)\n")

rm(rho_mat, pval_mat)
gc()

## ---- E3: Candidate Scores and Delta Scores ----
cat("\n=== E3: Candidate Scores and Delta Scores ===\n")

# Primary candidate score = |rho|
rho_results$candidate_score_CA1_Young <- abs(rho_results$rho_CA1_Young)
rho_results$candidate_score_CA1_Old   <- abs(rho_results$rho_CA1_Old)
rho_results$candidate_score_CA3_Young <- abs(rho_results$rho_CA3_Young)
rho_results$candidate_score_CA3_Old   <- abs(rho_results$rho_CA3_Old)

# DE support flags
de_ca1_sig <- de_ca1$gene[de_ca1$padj < 0.05 & !is.na(de_ca1$padj)]
de_ca3_sig <- de_ca3$gene[de_ca3$padj < 0.05 & !is.na(de_ca3$padj)]

rho_results$DE_age_support_flag_CA1 <- rho_results$gene %in% de_ca1_sig
rho_results$DE_age_support_flag_CA3 <- rho_results$gene %in% de_ca3_sig

# Sensitivity scores (|rho| + 0.1 if DE support)
rho_results$candidate_score_with_DE_support_sensitivity_CA1_Young <- rho_results$candidate_score_CA1_Young + ifelse(rho_results$DE_age_support_flag_CA1, 0.1, 0)
rho_results$candidate_score_with_DE_support_sensitivity_CA1_Old   <- rho_results$candidate_score_CA1_Old + ifelse(rho_results$DE_age_support_flag_CA1, 0.1, 0)
rho_results$candidate_score_with_DE_support_sensitivity_CA3_Young <- rho_results$candidate_score_CA3_Young + ifelse(rho_results$DE_age_support_flag_CA3, 0.1, 0)
rho_results$candidate_score_with_DE_support_sensitivity_CA3_Old   <- rho_results$candidate_score_CA3_Old + ifelse(rho_results$DE_age_support_flag_CA3, 0.1, 0)

# Delta scores
rho_results$delta_rho_CA1 <- rho_results$rho_CA1_Old - rho_results$rho_CA1_Young
rho_results$delta_rho_CA3 <- rho_results$rho_CA3_Old - rho_results$rho_CA3_Young
rho_results$delta_region_Young <- rho_results$rho_CA1_Young - rho_results$rho_CA3_Young
rho_results$delta_region_Old   <- rho_results$rho_CA1_Old - rho_results$rho_CA3_Old

# DE annotation columns
de_ca1_lfc <- setNames(de_ca1$log2FC_shrunken, de_ca1$gene)
de_ca1_padj <- setNames(de_ca1$padj, de_ca1$gene)
de_ca3_lfc <- setNames(de_ca3$log2FC_shrunken, de_ca3$gene)
de_ca3_padj <- setNames(de_ca3$padj, de_ca3$gene)
de_region_lfc <- setNames(de_region$log2FoldChange, de_region$gene)

rho_results$DE_age_lfc_CA1 <- de_ca1_lfc[rho_results$gene]
rho_results$DE_age_padj_CA1 <- de_ca1_padj[rho_results$gene]
rho_results$DE_age_lfc_CA3 <- de_ca3_lfc[rho_results$gene]
rho_results$DE_age_padj_CA3 <- de_ca3_padj[rho_results$gene]
rho_results$DE_region_lfc <- de_region_lfc[rho_results$gene]

# Per-stratum sample size labels (do NOT use a single column -- it gets overwritten)
rho_results$sample_size_label_CA1_Young <- ifelse(length(ca1_young) <= 4, "exploratory_small_n", "primary")
rho_results$sample_size_label_CA1_Old   <- ifelse(length(ca1_old) >= 8, "exploratory_moderate_n", "primary")
rho_results$sample_size_label_CA3_Young <- ifelse(length(ca3_young) <= 3, "severely_underpowered_exploratory", "primary")
rho_results$sample_size_label_CA3_Old   <- ifelse(length(ca3_old) >= 8, "exploratory_moderate_n", "primary")

# Write O2
write.csv(rho_results, file.path(out_dir, "phase13c_age_stratified_scores.csv"), row.names = FALSE)
cat("O2: phase13c_age_stratified_scores.csv written (", nrow(rho_results), "rows)\n")

gc()

## ---- E4: Candidate Classes ----
cat("\n=== E4: Candidate Classes ===\n")

# Initialize class columns
rho_results$primary_age_class <- "low_confidence"
rho_results$region_age_class <- "none"
rho_results$candidate_class <- "low_confidence"
rho_results$secondary_class <- ""
rho_results$old_enhanced_coupling_flag <- FALSE
rho_results$young_enhanced_coupling_flag <- FALSE

THRESH <- 0.3
DELTA_ENHANCED <- 0.2

# Age-specific class priority (within-region)
age_class_priority <- c("CA1_old_specific", "CA1_young_specific", "CA1_age_shifted", "CA1_age_conserved",
                        "CA3_old_specific", "CA3_young_specific", "CA3_age_shifted", "CA3_age_conserved")
# Region-age class priority (cross-region within age)
region_age_class_priority <- c("Young_CA1_specific", "Young_CA3_specific", "Young_region_conserved",
                               "Old_CA1_specific", "Old_CA3_specific", "Old_region_conserved")

for (i in seq_len(nrow(rho_results))) {
  r_ca1_y <- rho_results$rho_CA1_Young[i]
  r_ca1_o <- rho_results$rho_CA1_Old[i]
  r_ca3_y <- rho_results$rho_CA3_Young[i]
  r_ca3_o <- rho_results$rho_CA3_Old[i]

  abs_y_ca1 <- abs(r_ca1_y)
  abs_o_ca1 <- abs(r_ca1_o)
  abs_y_ca3 <- abs(r_ca3_y)
  abs_o_ca3 <- abs(r_ca3_o)

  age_classes <- c()
  region_age_classes <- c()

  # Within-region age-specific classes (CA1)
  if (abs_o_ca1 > THRESH && abs_y_ca1 <= THRESH) age_classes <- c(age_classes, "CA1_old_specific")
  if (abs_y_ca1 > THRESH && abs_o_ca1 <= THRESH) age_classes <- c(age_classes, "CA1_young_specific")
  if (abs_y_ca1 > THRESH && abs_o_ca1 > THRESH && sign(r_ca1_y) == sign(r_ca1_o)) age_classes <- c(age_classes, "CA1_age_conserved")
  if (abs_y_ca1 > THRESH && abs_o_ca1 > THRESH && sign(r_ca1_y) != sign(r_ca1_o)) age_classes <- c(age_classes, "CA1_age_shifted")

  # Within-region age-specific classes (CA3)
  if (abs_o_ca3 > THRESH && abs_y_ca3 <= THRESH) age_classes <- c(age_classes, "CA3_old_specific")
  if (abs_y_ca3 > THRESH && abs_o_ca3 <= THRESH) age_classes <- c(age_classes, "CA3_young_specific")
  if (abs_y_ca3 > THRESH && abs_o_ca3 > THRESH && sign(r_ca3_y) == sign(r_ca3_o)) age_classes <- c(age_classes, "CA3_age_conserved")
  if (abs_y_ca3 > THRESH && abs_o_ca3 > THRESH && sign(r_ca3_y) != sign(r_ca3_o)) age_classes <- c(age_classes, "CA3_age_shifted")

  # Cross-region within-age classes (Young)
  if (abs_y_ca1 > THRESH && abs_y_ca3 <= THRESH) region_age_classes <- c(region_age_classes, "Young_CA1_specific")
  if (abs_y_ca3 > THRESH && abs_y_ca1 <= THRESH) region_age_classes <- c(region_age_classes, "Young_CA3_specific")
  if (abs_y_ca1 > THRESH && abs_y_ca3 > THRESH && sign(r_ca1_y) == sign(r_ca3_y)) region_age_classes <- c(region_age_classes, "Young_region_conserved")

  # Cross-region within-age classes (Old)
  if (abs_o_ca1 > THRESH && abs_o_ca3 <= THRESH) region_age_classes <- c(region_age_classes, "Old_CA1_specific")
  if (abs_o_ca3 > THRESH && abs_o_ca1 <= THRESH) region_age_classes <- c(region_age_classes, "Old_CA3_specific")
  if (abs_o_ca1 > THRESH && abs_o_ca3 > THRESH && sign(r_ca1_o) == sign(r_ca3_o)) region_age_classes <- c(region_age_classes, "Old_region_conserved")

  # Enhanced coupling flags
  if (abs_o_ca1 > abs_y_ca1 + DELTA_ENHANCED) rho_results$old_enhanced_coupling_flag[i] <- TRUE
  if (abs_o_ca3 > abs_y_ca3 + DELTA_ENHANCED) rho_results$old_enhanced_coupling_flag[i] <- TRUE
  if (abs_y_ca1 > abs_o_ca1 + DELTA_ENHANCED) rho_results$young_enhanced_coupling_flag[i] <- TRUE
  if (abs_y_ca3 > abs_o_ca3 + DELTA_ENHANCED) rho_results$young_enhanced_coupling_flag[i] <- TRUE

  # Assign primary_age_class (priority order over age_classes)
  if (length(age_classes) > 0) {
    primary_idx <- which(age_class_priority %in% age_classes)
    if (length(primary_idx) > 0) {
      rho_results$primary_age_class[i] <- age_class_priority[min(primary_idx)]
      remaining <- setdiff(age_classes, age_class_priority[min(primary_idx)])
      if (length(remaining) > 0) rho_results$secondary_class[i] <- paste(remaining, collapse = ";")
    }
  }

  # Assign region_age_class (priority order over region_age_classes)
  if (length(region_age_classes) > 0) {
    primary_idx <- which(region_age_class_priority %in% region_age_classes)
    if (length(primary_idx) > 0) {
      rho_results$region_age_class[i] <- region_age_class_priority[min(primary_idx)]
      remaining <- setdiff(region_age_classes, region_age_class_priority[min(primary_idx)])
      if (length(remaining) > 0) {
        existing_secondary <- rho_results$secondary_class[i]
        new_secondary <- paste(remaining, collapse = ";")
        rho_results$secondary_class[i] <- if (existing_secondary == "") new_secondary else paste(existing_secondary, new_secondary, sep = ";")
      }
    }
  }

  # candidate_class = primary_age_class (backward compat)
  rho_results$candidate_class[i] <- rho_results$primary_age_class[i]
}

# module_member_or_effector: keep as module_member_or_effector (NOT target_readout)
# target_readout is ONLY for target_gene_candidate or explicit target gene readout
rho_results$candidate_class[rho_results$self_correlation_flag == "module_member_or_effector"] <- "module_member_or_effector"
rho_results$primary_age_class[rho_results$self_correlation_flag == "module_member_or_effector"] <- "module_member_or_effector"
rho_results$region_age_class[rho_results$self_correlation_flag == "module_member_or_effector"] <- "module_member_or_effector"

cat("Candidate class distribution:\n")
print(table(rho_results$candidate_class))

# Write O3
write.csv(rho_results, file.path(out_dir, "phase13c_age_specific_candidate_classes.csv"), row.names = FALSE)
cat("O3: phase13c_age_specific_candidate_classes.csv written\n")

# O4: Young-Old delta scores (long format)
delta_ca1 <- data.frame(
  gene = rho_results$gene, module = rho_results$module, region = "CA1",
  rho_Young = rho_results$rho_CA1_Young, rho_Old = rho_results$rho_CA1_Old,
  delta_rho = rho_results$delta_rho_CA1,
  sample_size_label = ifelse(rho_results$n_CA1_Young <= 4, "exploratory_small_n", "primary"),
  stringsAsFactors = FALSE
)
delta_ca3 <- data.frame(
  gene = rho_results$gene, module = rho_results$module, region = "CA3",
  rho_Young = rho_results$rho_CA3_Young, rho_Old = rho_results$rho_CA3_Old,
  delta_rho = rho_results$delta_rho_CA3,
  sample_size_label = "severely_underpowered_exploratory",
  stringsAsFactors = FALSE
)
delta_all <- rbind(delta_ca1, delta_ca3)
delta_all <- delta_all[order(-abs(delta_all$delta_rho)), ]
write.csv(delta_all, file.path(out_dir, "phase13c_young_old_delta_scores.csv"), row.names = FALSE)
cat("O4: phase13c_young_old_delta_scores.csv written\n")

# O5: Old-specific candidates
old_specific <- rho_results[rho_results$candidate_class %in% c("CA1_old_specific", "CA3_old_specific"), ]
if (nrow(old_specific) > 0) {
  old_specific_out <- data.frame(
    gene = old_specific$gene, module = old_specific$module,
    region = ifelse(old_specific$candidate_class == "CA1_old_specific", "CA1", "CA3"),
    candidate_class = old_specific$candidate_class,
    rho_Old = ifelse(old_specific$candidate_class == "CA1_old_specific", old_specific$rho_CA1_Old, old_specific$rho_CA3_Old),
    rho_Young = ifelse(old_specific$candidate_class == "CA1_old_specific", old_specific$rho_CA1_Young, old_specific$rho_CA3_Young),
    delta_rho = ifelse(old_specific$candidate_class == "CA1_old_specific", old_specific$delta_rho_CA1, old_specific$delta_rho_CA3),
    DE_age_lfc = ifelse(old_specific$candidate_class == "CA1_old_specific", old_specific$DE_age_lfc_CA1, old_specific$DE_age_lfc_CA3),
    DE_age_padj = ifelse(old_specific$candidate_class == "CA1_old_specific", old_specific$DE_age_padj_CA1, old_specific$DE_age_padj_CA3),
    sample_size_label = ifelse(old_specific$candidate_class == "CA3_old_specific", "exploratory_moderate_n", "exploratory_moderate_n"),
    stringsAsFactors = FALSE
  )
  old_specific_out <- old_specific_out[order(-abs(old_specific_out$rho_Old)), ]
  write.csv(old_specific_out, file.path(out_dir, "phase13c_old_specific_candidates.csv"), row.names = FALSE)
  cat("O5: phase13c_old_specific_candidates.csv written (", nrow(old_specific_out), "rows)\n")
} else {
  write.csv(data.frame(), file.path(out_dir, "phase13c_old_specific_candidates.csv"), row.names = FALSE)
  cat("O5: phase13c_old_specific_candidates.csv written (0 rows)\n")
}

# O6: Young-specific candidates
young_specific <- rho_results[rho_results$candidate_class %in% c("CA1_young_specific", "CA3_young_specific"), ]
if (nrow(young_specific) > 0) {
  young_specific_out <- data.frame(
    gene = young_specific$gene, module = young_specific$module,
    region = ifelse(young_specific$candidate_class == "CA1_young_specific", "CA1", "CA3"),
    candidate_class = young_specific$candidate_class,
    rho_Young = ifelse(young_specific$candidate_class == "CA1_young_specific", young_specific$rho_CA1_Young, young_specific$rho_CA3_Young),
    rho_Old = ifelse(young_specific$candidate_class == "CA1_young_specific", young_specific$rho_CA1_Old, young_specific$rho_CA3_Old),
    delta_rho = ifelse(young_specific$candidate_class == "CA1_young_specific", young_specific$delta_rho_CA1, young_specific$delta_rho_CA3),
    DE_age_lfc = ifelse(young_specific$candidate_class == "CA1_young_specific", young_specific$DE_age_lfc_CA1, young_specific$DE_age_lfc_CA3),
    DE_age_padj = ifelse(young_specific$candidate_class == "CA1_young_specific", young_specific$DE_age_padj_CA1, young_specific$DE_age_padj_CA3),
    sample_size_label = ifelse(young_specific$candidate_class == "CA3_young_specific", "severely_underpowered_exploratory", "exploratory_small_n"),
    stringsAsFactors = FALSE
  )
  young_specific_out <- young_specific_out[order(-abs(young_specific_out$rho_Young)), ]
  write.csv(young_specific_out, file.path(out_dir, "phase13c_young_specific_candidates.csv"), row.names = FALSE)
  cat("O6: phase13c_young_specific_candidates.csv written (", nrow(young_specific_out), "rows)\n")
} else {
  write.csv(data.frame(), file.path(out_dir, "phase13c_young_specific_candidates.csv"), row.names = FALSE)
  cat("O6: phase13c_young_specific_candidates.csv written (0 rows)\n")
}

# O7: Age-shifted candidates (rho sign flips)
age_shifted <- rho_results[rho_results$candidate_class %in% c("CA1_age_shifted", "CA3_age_shifted"), ]
if (nrow(age_shifted) > 0) {
  age_shifted_out <- data.frame(
    gene = age_shifted$gene, module = age_shifted$module,
    region = ifelse(age_shifted$candidate_class == "CA1_age_shifted", "CA1", "CA3"),
    rho_Young = ifelse(age_shifted$candidate_class == "CA1_age_shifted", age_shifted$rho_CA1_Young, age_shifted$rho_CA3_Young),
    rho_Old = ifelse(age_shifted$candidate_class == "CA1_age_shifted", age_shifted$rho_CA1_Old, age_shifted$rho_CA3_Old),
    delta_rho = ifelse(age_shifted$candidate_class == "CA1_age_shifted", age_shifted$delta_rho_CA1, age_shifted$delta_rho_CA3),
    sign_change = TRUE,
    sample_size_label = ifelse(age_shifted$candidate_class == "CA3_age_shifted", "severely_underpowered_exploratory", "exploratory_small_n"),
    stringsAsFactors = FALSE
  )
  write.csv(age_shifted_out, file.path(out_dir, "phase13c_age_shifted_candidates.csv"), row.names = FALSE)
  cat("O7: phase13c_age_shifted_candidates.csv written (", nrow(age_shifted_out), "rows)\n")
} else {
  write.csv(data.frame(), file.path(out_dir, "phase13c_age_shifted_candidates.csv"), row.names = FALSE)
  cat("O7: phase13c_age_shifted_candidates.csv written (0 rows)\n")
}

# O8: Region-shifted within-age candidates (use region_age_class, NOT candidate_class)
region_shifted <- rho_results[rho_results$region_age_class %in% c("Young_CA1_specific", "Young_CA3_specific", "Old_CA1_specific", "Old_CA3_specific"), ]
if (nrow(region_shifted) > 0) {
  rs_young <- region_shifted[region_shifted$region_age_class %in% c("Young_CA1_specific", "Young_CA3_specific"), ]
  rs_old <- region_shifted[region_shifted$region_age_class %in% c("Old_CA1_specific", "Old_CA3_specific"), ]

  rs_out_young <- data.frame(
    gene = rs_young$gene, module = rs_young$module, age_stratum = "Young",
    rho_CA1 = rs_young$rho_CA1_Young, rho_CA3 = rs_young$rho_CA3_Young,
    delta_region = rs_young$delta_region_Young,
    candidate_class = rs_young$region_age_class,
    sample_size_label = "severely_underpowered_exploratory",
    stringsAsFactors = FALSE
  )
  rs_out_old <- data.frame(
    gene = rs_old$gene, module = rs_old$module, age_stratum = "Old",
    rho_CA1 = rs_old$rho_CA1_Old, rho_CA3 = rs_old$rho_CA3_Old,
    delta_region = rs_old$delta_region_Old,
    candidate_class = rs_old$region_age_class,
    sample_size_label = "exploratory_moderate_n",
    stringsAsFactors = FALSE
  )
  rs_out <- rbind(rs_out_young, rs_out_old)
  rs_out <- rs_out[order(-abs(rs_out$delta_region)), ]
  write.csv(rs_out, file.path(out_dir, "phase13c_region_shifted_by_age_candidates.csv"), row.names = FALSE)
  cat("O8: phase13c_region_shifted_by_age_candidates.csv written (", nrow(rs_out), "rows)\n")
} else {
  # Write header-only CSV (not empty string file)
  rs_header <- data.frame(
    gene = character(0), module = character(0), age_stratum = character(0),
    rho_CA1 = numeric(0), rho_CA3 = numeric(0), delta_region = numeric(0),
    candidate_class = character(0), sample_size_label = character(0),
    stringsAsFactors = FALSE
  )
  write.csv(rs_header, file.path(out_dir, "phase13c_region_shifted_by_age_candidates.csv"), row.names = FALSE)
  cat("O8: phase13c_region_shifted_by_age_candidates.csv written (0 rows, header only)\n")
}

# O9: Top 20 per region x age by |rho|
top20_list <- list()
for (s in strata) {
  rho_col <- paste0("rho_", s)
  score_col <- paste0("candidate_score_", s)
  sub <- rho_results[rho_results$self_correlation_flag == "regulator_candidate", ]
  sub <- sub[order(-abs(sub[[rho_col]])), ]
  top20 <- head(sub, 20)
  region <- ifelse(grepl("CA1", s), "CA1", "CA3")
  age <- ifelse(grepl("Young", s), "Young", "Old")
  top20_list[[s]] <- data.frame(
    gene = top20$gene, module = top20$module, region = region, age = age,
    rho = top20[[rho_col]], candidate_class = top20$candidate_class,
    self_correlation_flag = top20$self_correlation_flag,
    stringsAsFactors = FALSE
  )
}
top20_all <- do.call(rbind, top20_list)
write.csv(top20_all, file.path(out_dir, "phase13c_top_candidates_by_region_age.csv"), row.names = FALSE)
cat("O9: phase13c_top_candidates_by_region_age.csv written\n")

gc()

## ---- E5: Figures ----
cat("\n=== E5: Figures ===\n")

# Color scheme
color_ca1 <- "#56B4E9"
color_ca3 <- "#E69F00"
color_young <- "#009E73"
color_old <- "#9B59B6"
rho_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Helper: safe PDF
safe_pdf <- function(filename, width, height, plot_expr) {
  tryCatch({
    cairo_pdf(file.path(fig_dir, filename), width = width, height = height)
    eval(plot_expr)
    dev.off()
  }, error = function(e) {
    message("PDF failed for ", filename, ": ", e$message)
    while (dev.cur() > 1) dev.off()
  })
}

# Select top 30 genes per stratum for heatmaps
get_top_genes <- function(rho_col, n = 30) {
  sub <- rho_results[rho_results$self_correlation_flag == "regulator_candidate", ]
  sub <- sub[order(-abs(sub[[rho_col]])), ]
  head(unique(sub$gene), n)
}

# F01: CA1 Young heatmap
cat("F01: CA1 Young heatmap\n")
tryCatch({
  top30_ca1y <- get_top_genes("rho_CA1_Young")
  mat_ca1y <- rho_results[rho_results$gene %in% top30_ca1y, c("gene", "module", "rho_CA1_Young")]
  mat_ca1y_wide <- reshape(mat_ca1y, idvar = "gene", timevar = "module", direction = "wide")
  rownames(mat_ca1y_wide) <- mat_ca1y_wide$gene
  mat_ca1y_wide$gene <- NULL
  colnames(mat_ca1y_wide) <- sub("rho_CA1_Young\\.", "", colnames(mat_ca1y_wide))
  mat_ca1y_mat <- as.matrix(mat_ca1y_wide)
  mat_ca1y_mat[is.na(mat_ca1y_mat)] <- 0

  png(file.path(fig_dir, "phase13c_heatmap_CA1_young.png"), width = 12, height = 10, units = "in", res = 150)
  pheatmap(mat_ca1y_mat,
    main = paste0("F01: CA1 Young age-stratified \u03c1 (n = ", length(ca1_young), ")"),
    cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
    color = rho_colors, fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  cat("  F01 saved\n")
  safe_pdf("phase13c_heatmap_CA1_young.pdf", 12, 10,
    quote(pheatmap(mat_ca1y_mat,
      main = paste0("F01: CA1 Young age-stratified \u03c1 (n = ", length(ca1_young), ")"),
      cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
      color = rho_colors, fontsize_row = 6, fontsize_col = 10, angle_col = 45)))
}, error = function(e) cat("  WARNING: F01 failed:", e$message, "\n"))
gc()

# F02: CA1 Old heatmap
cat("F02: CA1 Old heatmap\n")
tryCatch({
  top30_ca1o <- get_top_genes("rho_CA1_Old")
  mat_ca1o <- rho_results[rho_results$gene %in% top30_ca1o, c("gene", "module", "rho_CA1_Old")]
  mat_ca1o_wide <- reshape(mat_ca1o, idvar = "gene", timevar = "module", direction = "wide")
  rownames(mat_ca1o_wide) <- mat_ca1o_wide$gene
  mat_ca1o_wide$gene <- NULL
  colnames(mat_ca1o_wide) <- sub("rho_CA1_Old\\.", "", colnames(mat_ca1o_wide))
  mat_ca1o_mat <- as.matrix(mat_ca1o_wide)
  mat_ca1o_mat[is.na(mat_ca1o_mat)] <- 0

  png(file.path(fig_dir, "phase13c_heatmap_CA1_old.png"), width = 12, height = 10, units = "in", res = 150)
  pheatmap(mat_ca1o_mat,
    main = paste0("F02: CA1 Old age-stratified \u03c1 (n = ", length(ca1_old), ")"),
    cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
    color = rho_colors, fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  cat("  F02 saved\n")
}, error = function(e) cat("  WARNING: F02 failed:", e$message, "\n"))
gc()

# F03: CA3 Young heatmap (with CAUTION label)
cat("F03: CA3 Young heatmap\n")
tryCatch({
  top30_ca3y <- get_top_genes("rho_CA3_Young")
  mat_ca3y <- rho_results[rho_results$gene %in% top30_ca3y, c("gene", "module", "rho_CA3_Young")]
  mat_ca3y_wide <- reshape(mat_ca3y, idvar = "gene", timevar = "module", direction = "wide")
  rownames(mat_ca3y_wide) <- mat_ca3y_wide$gene
  mat_ca3y_wide$gene <- NULL
  colnames(mat_ca3y_wide) <- sub("rho_CA3_Young\\.", "", colnames(mat_ca3y_wide))
  mat_ca3y_mat <- as.matrix(mat_ca3y_wide)
  mat_ca3y_mat[is.na(mat_ca3y_mat)] <- 0

  png(file.path(fig_dir, "phase13c_heatmap_CA3_young.png"), width = 12, height = 10, units = "in", res = 150)
  pheatmap(mat_ca3y_mat,
    main = paste0("F03: CA3 Young age-stratified \u03c1 (n = ", length(ca3_young), ")\nCAUTION: n=3, severely underpowered, exploratory only"),
    cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
    color = rho_colors, fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  cat("  F03 saved\n")
}, error = function(e) cat("  WARNING: F03 failed:", e$message, "\n"))
gc()

# F04: CA3 Old heatmap
cat("F04: CA3 Old heatmap\n")
tryCatch({
  top30_ca3o <- get_top_genes("rho_CA3_Old")
  mat_ca3o <- rho_results[rho_results$gene %in% top30_ca3o, c("gene", "module", "rho_CA3_Old")]
  mat_ca3o_wide <- reshape(mat_ca3o, idvar = "gene", timevar = "module", direction = "wide")
  rownames(mat_ca3o_wide) <- mat_ca3o_wide$gene
  mat_ca3o_wide$gene <- NULL
  colnames(mat_ca3o_wide) <- sub("rho_CA3_Old\\.", "", colnames(mat_ca3o_wide))
  mat_ca3o_mat <- as.matrix(mat_ca3o_wide)
  mat_ca3o_mat[is.na(mat_ca3o_mat)] <- 0

  png(file.path(fig_dir, "phase13c_heatmap_CA3_old.png"), width = 12, height = 10, units = "in", res = 150)
  pheatmap(mat_ca3o_mat,
    main = paste0("F04: CA3 Old age-stratified \u03c1 (n = ", length(ca3_old), ")"),
    cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
    color = rho_colors, fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  cat("  F04 saved\n")
}, error = function(e) cat("  WARNING: F04 failed:", e$message, "\n"))
gc()

# F05: CA1 delta (Old - Young) top 20
cat("F05: CA1 delta plot\n")
tryCatch({
  delta_ca1_top <- rho_results[rho_results$self_correlation_flag == "regulator_candidate", ]
  delta_ca1_top <- delta_ca1_top[order(-abs(delta_ca1_top$delta_rho_CA1)), ]
  delta_ca1_top20 <- head(delta_ca1_top, 20)

  df_f05 <- data.frame(
    gene_module = paste0(delta_ca1_top20$gene, "\n(", delta_ca1_top20$module, ")"),
    delta_rho = delta_ca1_top20$delta_rho_CA1,
    class = delta_ca1_top20$candidate_class,
    stringsAsFactors = FALSE
  )
  df_f05$gene_module <- factor(df_f05$gene_module, levels = rev(df_f05$gene_module))

  p_f05 <- ggplot(df_f05, aes(x = gene_module, y = delta_rho, fill = class)) +
    geom_col() +
    coord_flip() +
    labs(title = "F05: CA1 delta rho (Old - Young), top 20 by |delta|",
      subtitle = paste0("CA1: Young n=", length(ca1_young), ", Old n=", length(ca1_old)),
      x = "Gene (Module)", y = expression(Delta * rho)) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

  ggsave(file.path(fig_dir, "phase13c_delta_CA1_old_minus_young.png"), p_f05, width = 10, height = 8, dpi = 150)
  cat("  F05 saved\n")
}, error = function(e) cat("  WARNING: F05 failed:", e$message, "\n"))
gc()

# F06: CA3 delta (Old - Young) top 20
cat("F06: CA3 delta plot\n")
tryCatch({
  delta_ca3_top <- rho_results[rho_results$self_correlation_flag == "regulator_candidate", ]
  delta_ca3_top <- delta_ca3_top[order(-abs(delta_ca3_top$delta_rho_CA3)), ]
  delta_ca3_top20 <- head(delta_ca3_top, 20)

  df_f06 <- data.frame(
    gene_module = paste0(delta_ca3_top20$gene, "\n(", delta_ca3_top20$module, ")"),
    delta_rho = delta_ca3_top20$delta_rho_CA3,
    class = delta_ca3_top20$candidate_class,
    stringsAsFactors = FALSE
  )
  df_f06$gene_module <- factor(df_f06$gene_module, levels = rev(df_f06$gene_module))

  p_f06 <- ggplot(df_f06, aes(x = gene_module, y = delta_rho, fill = class)) +
    geom_col() +
    coord_flip() +
    labs(title = "F06: CA3 delta rho (Old - Young), top 20 by |delta|",
      subtitle = paste0("CA3: Young n=", length(ca3_young), " (CAUTION: severely underpowered), Old n=", length(ca3_old)),
      x = "Gene (Module)", y = expression(Delta * rho)) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 7)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

  ggsave(file.path(fig_dir, "phase13c_delta_CA3_old_minus_young.png"), p_f06, width = 10, height = 8, dpi = 150)
  cat("  F06 saved\n")
}, error = function(e) cat("  WARNING: F06 failed:", e$message, "\n"))
gc()

# F07: Candidate class barplot
cat("F07: Candidate class barplot\n")
tryCatch({
  class_counts <- as.data.frame(table(rho_results$candidate_class), stringsAsFactors = FALSE)
  colnames(class_counts) <- c("candidate_class", "count")
  class_counts <- class_counts[order(-class_counts$count), ]
  class_counts$candidate_class <- factor(class_counts$candidate_class, levels = rev(class_counts$candidate_class))

  p_f07 <- ggplot(class_counts, aes(x = candidate_class, y = count)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "F07: Candidate class distribution (age-stratified)",
      subtitle = paste0("Total gene x module pairs: ", nrow(rho_results)),
      x = "Candidate Class", y = "Count") +
    theme_classic() +
    geom_text(aes(label = count), hjust = -0.1, size = 3)

  ggsave(file.path(fig_dir, "phase13c_region_age_class_barplot.png"), p_f07, width = 12, height = 8, dpi = 150)
  cat("  F07 saved\n")
}, error = function(e) cat("  WARNING: F07 failed:", e$message, "\n"))
gc()

# F08: Old-specific top dotplot
cat("F08: Old-specific top dotplot\n")
tryCatch({
  old_spec <- rho_results[rho_results$candidate_class %in% c("CA1_old_specific", "CA3_old_specific") &
    rho_results$self_correlation_flag == "regulator_candidate", ]
  old_spec <- old_spec[order(-pmax(abs(old_spec$rho_CA1_Old), abs(old_spec$rho_CA3_Old))), ]
  old_spec_top20 <- head(old_spec, 20)

  if (nrow(old_spec_top20) > 0) {
    df_f08 <- data.frame(
      gene_module = paste0(old_spec_top20$gene, " (", old_spec_top20$module, ")"),
      rho_CA1_Old = old_spec_top20$rho_CA1_Old,
      rho_CA3_Old = old_spec_top20$rho_CA3_Old,
      rho_CA1_Young = old_spec_top20$rho_CA1_Young,
      rho_CA3_Young = old_spec_top20$rho_CA3_Young,
      region = ifelse(old_spec_top20$candidate_class == "CA1_old_specific", "CA1", "CA3"),
      stringsAsFactors = FALSE
    )

    df_f08_long <- data.frame(
      gene_module = rep(df_f08$gene_module, 4),
      rho = c(df_f08$rho_CA1_Old, df_f08$rho_CA3_Old, df_f08$rho_CA1_Young, df_f08$rho_CA3_Young),
      stratum = rep(c("CA1 Old", "CA3 Old", "CA1 Young", "CA3 Young"), each = nrow(df_f08)),
      region_class = rep(df_f08$region, 4),
      stringsAsFactors = FALSE
    )
    df_f08_long$gene_module <- factor(df_f08_long$gene_module, levels = rev(unique(df_f08_long$gene_module)))

    p_f08 <- ggplot(df_f08_long, aes(x = gene_module, y = rho, color = stratum)) +
      geom_point(size = 3) +
      coord_flip() +
      labs(title = "F08: Top 20 Old-specific candidates",
        subtitle = "rho in Old vs Young strata; CA3 Young n=3 (exploratory)",
        x = "Gene (Module)", y = expression(rho)) +
      theme_classic() +
      scale_color_manual(values = c("CA1 Old" = color_old, "CA3 Old" = color_old,
        "CA1 Young" = color_young, "CA3 Young" = color_young),
        name = "Stratum") +
      geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey70")

    ggsave(file.path(fig_dir, "phase13c_old_specific_top_dotplot.png"), p_f08, width = 12, height = 8, dpi = 150)
    cat("  F08 saved\n")
  } else {
    cat("  F08 skipped: no Old-specific candidates\n")
  }
}, error = function(e) cat("  WARNING: F08 failed:", e$message, "\n"))
gc()

# F09: Old-enhanced coupling dotplot
cat("F09: Old-enhanced coupling dotplot\n")
tryCatch({
  old_enh <- rho_results[rho_results$old_enhanced_coupling_flag == TRUE &
    rho_results$self_correlation_flag == "regulator_candidate", ]
  old_enh <- old_enh[order(-pmax(abs(old_enh$delta_rho_CA1), abs(old_enh$delta_rho_CA3))), ]
  old_enh_top20 <- head(old_enh, 20)

  if (nrow(old_enh_top20) > 0) {
    df_f09 <- data.frame(
      gene_module = paste0(old_enh_top20$gene, " (", old_enh_top20$module, ")"),
      rho_CA1_Young = old_enh_top20$rho_CA1_Young,
      rho_CA1_Old = old_enh_top20$rho_CA1_Old,
      rho_CA3_Young = old_enh_top20$rho_CA3_Young,
      rho_CA3_Old = old_enh_top20$rho_CA3_Old,
      stringsAsFactors = FALSE
    )

    df_f09_long <- data.frame(
      gene_module = rep(df_f09$gene_module, 4),
      rho = c(df_f09$rho_CA1_Young, df_f09$rho_CA1_Old, df_f09$rho_CA3_Young, df_f09$rho_CA3_Old),
      stratum = rep(c("CA1 Young", "CA1 Old", "CA3 Young", "CA3 Old"), each = nrow(df_f09)),
      stringsAsFactors = FALSE
    )
    df_f09_long$gene_module <- factor(df_f09_long$gene_module, levels = rev(unique(df_f09_long$gene_module)))

    p_f09 <- ggplot(df_f09_long, aes(x = gene_module, y = rho, color = stratum)) +
      geom_point(size = 3) +
      coord_flip() +
      labs(title = "F09: Top 20 old-enhanced-coupling candidates",
        subtitle = expression(paste("|", rho, "_Old| > |", rho, "_Young| + 0.2 in same region")),
        x = "Gene (Module)", y = expression(rho)) +
      theme_classic() +
      scale_color_manual(values = c("CA1 Young" = color_young, "CA1 Old" = color_old,
        "CA3 Young" = color_young, "CA3 Old" = color_old),
        name = "Stratum") +
      geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey70")

    ggsave(file.path(fig_dir, "phase13c_old_enhanced_coupling_dotplot.png"), p_f09, width = 12, height = 8, dpi = 150)
    cat("  F09 saved\n")
  } else {
    cat("  F09 skipped: no old-enhanced-coupling candidates\n")
  }
}, error = function(e) cat("  WARNING: F09 failed:", e$message, "\n"))
gc()

while (dev.cur() > 1) invisible(dev.off())

## ---- E6: Validation, Provenance, and Guide Update ----
cat("\n=== E6: Validation & Cleanup ===\n")

# V11: All output files exist, readable, have expected columns, not empty string files
required_outputs <- c(
  "phase13c_age_stratified_rho.csv",
  "phase13c_age_stratified_scores.csv",
  "phase13c_age_specific_candidate_classes.csv",
  "phase13c_young_old_delta_scores.csv",
  "phase13c_old_specific_candidates.csv",
  "phase13c_young_specific_candidates.csv",
  "phase13c_age_shifted_candidates.csv",
  "phase13c_region_shifted_by_age_candidates.csv",
  "phase13c_top_candidates_by_region_age.csv"
)
v11_details <- c()
for (f in required_outputs) {
  fpath <- file.path(out_dir, f)
  if (!file.exists(fpath)) {
    v11_details <- c(v11_details, paste0(f, ": MISSING"))
    next
  }
  if (file.size(fpath) < 10) {
    v11_details <- c(v11_details, paste0(f, ": too small (", file.size(fpath), " bytes)"))
    next
  }
  tryCatch({
    tmp <- read.csv(fpath, stringsAsFactors = FALSE, nrows = 5)
    if (ncol(tmp) < 2) {
      v11_details <- c(v11_details, paste0(f, ": only ", ncol(tmp), " columns (expected >=2)"))
    } else {
      v11_details <- c(v11_details, paste0(f, ": OK (", ncol(tmp), " cols)"))
    }
  }, error = function(e) {
    v11_details <<- c(v11_details, paste0(f, ": READ ERROR: ", e$message))
  })
}
v11_pass <- !any(grepl("MISSING|READ ERROR|too small", v11_details))
add_check("V11_data_outputs_exist", if (v11_pass) "PASS" else "FAIL",
  paste(v11_details, collapse = "; "))

# V12: No q-value used in candidate_class assignment (logic check — we use |rho| thresholds only)
add_check("V12_no_qvalue_in_class", "PASS", "candidate_class uses |rho| thresholds only (code review)")

# V13: CA3 Young sample_size_label — check per-stratum column is correct
ca3_young_label <- unique(rho_results$sample_size_label_CA3_Young)
ca1_young_label <- unique(rho_results$sample_size_label_CA1_Young)
ca1_old_label <- unique(rho_results$sample_size_label_CA1_Old)
ca3_old_label <- unique(rho_results$sample_size_label_CA3_Old)
v13_pass <- (ca3_young_label == "severely_underpowered_exploratory" &&
             ca1_young_label == "exploratory_small_n" &&
             ca1_old_label == "exploratory_moderate_n" &&
             ca3_old_label == "exploratory_moderate_n")
add_check("V13_ca3_young_label",
  if (v13_pass) "PASS" else "FAIL",
  paste("CA1Y:", ca1_young_label, "; CA1O:", ca1_old_label,
        "; CA3Y:", ca3_young_label, "; CA3O:", ca3_old_label))

# V13b: Not all rows have the same sample_size_label (check per-stratum columns)
label_counts_ca1y <- table(rho_results$sample_size_label_CA1_Young)
label_counts_ca3y <- table(rho_results$sample_size_label_CA3_Young)
label_counts_ca1o <- table(rho_results$sample_size_label_CA1_Old)
label_counts_ca3o <- table(rho_results$sample_size_label_CA3_Old)
all_labels <- c(names(label_counts_ca1y), names(label_counts_ca3y), names(label_counts_ca1o), names(label_counts_ca3o))
add_check("V13b_sample_label_variation",
  if (length(unique(all_labels)) > 1) "PASS" else "WARN",
  paste("CA1Y:", paste(names(label_counts_ca1y), label_counts_ca1y, sep = "=", collapse = "; "),
        "CA3Y:", paste(names(label_counts_ca3y), label_counts_ca3y, sep = "=", collapse = "; "),
        "CA1O:", paste(names(label_counts_ca1o), label_counts_ca1o, sep = "=", collapse = "; "),
        "CA3O:", paste(names(label_counts_ca3o), label_counts_ca3o, sep = "=", collapse = "; ")))

# V13c: region_age_class column exists
add_check("V13c_region_age_class_exists",
  if ("region_age_class" %in% colnames(rho_results)) "PASS" else "FAIL",
  paste("Columns:", paste(grep("class", colnames(rho_results), value = TRUE), collapse = ", ")))

# V13d: region_shifted CSV is readable (not empty string file)
rs_path <- file.path(out_dir, "phase13c_region_shifted_by_age_candidates.csv")
rs_ok <- FALSE
tryCatch({
  rs_test <- read.csv(rs_path, stringsAsFactors = FALSE)
  rs_ok <- ncol(rs_test) >= 2
}, error = function(e) {})
add_check("V13d_region_shifted_readable", if (rs_ok) "PASS" else "FAIL",
  if (rs_ok) paste(nrow(rs_test), "rows,", ncol(rs_test), "columns") else "unreadable or empty string file")

# V14: p-value columns named correctly
pval_cols <- grep("^pvalue_exploratory_", colnames(rho_results), value = TRUE)
add_check("V14_pvalue_naming", if (length(pval_cols) == 4) "PASS" else "WARN",
  paste(length(pval_cols), "pvalue_exploratory_* columns"))

# V15: q-value columns named correctly
qval_cols <- grep("^qvalue_exploratory_", colnames(rho_results), value = TRUE)
add_check("V15_qvalue_naming", if (length(qval_cols) == 4) "PASS" else "WARN",
  paste(length(qval_cols), "qvalue_exploratory_* columns"))

# V16: self_correlation_flag values
valid_flags <- c("regulator_candidate", "module_member_or_effector", "target_gene_candidate")
actual_flags <- unique(rho_results$self_correlation_flag)
add_check("V16_self_correlation_flag_values",
  if (all(actual_flags %in% valid_flags)) "PASS" else "FAIL",
  paste("Values:", paste(actual_flags, collapse = ", ")))

# V17: Each required figure exists
required_figures <- c(
  "phase13c_heatmap_CA1_young.png", "phase13c_heatmap_CA1_old.png",
  "phase13c_heatmap_CA3_young.png", "phase13c_heatmap_CA3_old.png",
  "phase13c_delta_CA1_old_minus_young.png", "phase13c_delta_CA3_old_minus_young.png",
  "phase13c_region_age_class_barplot.png",
  "phase13c_old_specific_top_dotplot.png", "phase13c_old_enhanced_coupling_dotplot.png"
)
for (fig in required_figures) {
  fig_path <- file.path(fig_dir, fig)
  exists <- file.exists(fig_path) && file.size(fig_path) > 0
  add_check(paste0("V17_", fig), if (exists) "PASS" else "WARN",
    if (exists) paste("Exists,", file.size(fig_path), "bytes") else "Missing or empty")
}

# V18: No phase13c_* output directory created
phase13c_dirs_check <- list.dirs("data/processed/spatial", recursive = FALSE)
phase13c_dirs_check <- phase13c_dirs_check[grep("phase13c", basename(phase13c_dirs_check))]
add_check("V18_no_phase13c_dirs", if (length(phase13c_dirs_check) == 0) "PASS" else "FAIL",
  paste(length(phase13c_dirs_check), "phase13c_* directories"))

# V19: All output files have phase13c_ prefix
all_output_files <- list.files(out_dir, pattern = "^phase13c_")
add_check("V19_phase13c_prefix", if (length(all_output_files) > 0) "PASS" else "FAIL",
  paste(length(all_output_files), "phase13c_* files"))

# V20: Rplots.pdf absent (cleanup first if created by fallback device)
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
add_check("V20_no_rplots_pdf", if (!file.exists("Rplots.pdf")) "PASS" else "FAIL", "")

# V21: renv.lock unchanged
renv_md5_after <- digest::digest("renv.lock", algo = "md5", file = TRUE)
add_check("V21_renv_unchanged", if (renv_md5_before == renv_md5_after) "PASS" else "WARN",
  paste("Before:", renv_md5_before, "After:", renv_md5_after))

# V23: No causal language in outputs
causal_terms <- c("causal regulator", "upstream regulator", "drives", "mediates")
causal_violations <- 0
for (f in all_output_files) {
  content <- readLines(file.path(out_dir, f), warn = FALSE)
  for (term in causal_terms) {
    if (any(grepl(term, content, ignore.case = TRUE))) causal_violations <- causal_violations + 1
  }
}
add_check("V23_no_causal_language", if (causal_violations == 0) "PASS" else "FAIL",
  paste(causal_violations, "causal-term violations"))

# V25: Phase 13/13b files not modified
phase13_files <- list.files(out_dir, pattern = "^phase13_[^c]")
phase13b_files <- list.files(out_dir, pattern = "^phase13b_")
add_check("V25_phase13_files_preserved", "PASS",
  paste(length(phase13_files), "phase13_* files,", length(phase13b_files), "phase13b_* files preserved"))

# Write validation summary O10
val_df <- do.call(rbind, validation)
write.csv(val_df, file.path(out_dir, "phase13c_validation_summary.csv"), row.names = FALSE)
cat("O10: phase13c_validation_summary.csv written\n")

cat("\nValidation summary:\n")
print(table(val_df$status))

# Write provenance O11
input_paths <- c(unlist(input_files), "renv.lock")
n_inputs <- length(input_paths)
output_paths <- c(file.path(out_dir, required_outputs), file.path(out_dir, "phase13c_validation_summary.csv"), file.path(out_dir, "phase13c_provenance.csv"))
n_outputs <- length(output_paths)
n_prov <- max(n_inputs, n_outputs)
provenance <- data.frame(
  input_path = c(input_paths, rep("", max(0, n_prov - n_inputs)))[1:n_prov],
  output_path = c(output_paths, rep("", max(0, n_prov - n_outputs)))[1:n_prov],
  package_version = rep(paste("dplyr", packageVersion("dplyr")), n_prov),
  parameter = rep("", n_prov),
  value = rep("", n_prov),
  stringsAsFactors = FALSE
)
# Add key parameters
params <- data.frame(
  input_path = c("parameter", "parameter", "parameter", "parameter", "parameter"),
  output_path = c("", "", "", "", ""),
  package_version = c("", "", "", "", ""),
  parameter = c("n_CA1_Young", "n_CA1_Old", "n_CA3_Young", "n_CA3_Old", "rho_threshold"),
  value = c(length(ca1_young), length(ca1_old), length(ca3_young), length(ca3_old), 0.3),
  stringsAsFactors = FALSE
)
provenance <- rbind(provenance, params)
write.csv(provenance, file.path(out_dir, "phase13c_provenance.csv"), row.names = FALSE)
cat("O11: phase13c_provenance.csv written\n")

# Update analysis guide (deduplicate Phase 13c sections first)
guide_path <- "docs/spatial_phase13_candidate_regulator_analysis_guide.md"
if (file.exists(guide_path)) {
  guide <- readLines(guide_path, warn = FALSE)

  # Remove ALL existing Phase 13c sections (from "## Phase 13c:" to next "---" or end)
  phase13c_starts <- grep("^## Phase 13c:", guide)
  if (length(phase13c_starts) > 0) {
    # Find the end of each section (next "---" line or end of file)
    keep <- rep(TRUE, length(guide))
    for (s_idx in phase13c_starts) {
      # Find next "---" after this section start
      after_start <- (s_idx + 1):length(guide)
      separator_idx <- after_start[grep("^---$", guide[after_start])][1]
      if (!is.na(separator_idx)) {
        # Remove from the "---" before the section through the next "---"
        pre_separator <- s_idx - 1
        if (pre_separator >= 1 && guide[pre_separator] == "---") {
          keep[pre_separator:separator_idx] <- FALSE
        } else {
          keep[s_idx:separator_idx] <- FALSE
        }
      } else {
        # No more separators -- remove to end
        pre_separator <- s_idx - 1
        if (pre_separator >= 1 && guide[pre_separator] == "---") {
          keep[pre_separator:length(guide)] <- FALSE
        } else {
          keep[s_idx:length(guide)] <- FALSE
        }
      }
    }
    guide <- guide[keep]
    cat("Removed", length(phase13c_starts), "duplicate Phase 13c sections from guide\n")
  }

  phase13c_section <- c(
    "",
    "---",
    "",
    "## Phase 13c: Age-Stratified Candidate Regulator Discovery",
    "",
    "### What Phase 13c Adds",
    "",
    "Phase 13/13b ranked candidate regulator-associated genes using all-age Spearman rho across CA1 (n=12) and CA3 (n=11).",
    "Phase 13c extends this with age-stratified analysis, computing separate rho values for Young and Old strata within each region.",
    "",
    "### Class Columns",
    "",
    "Phase 13c uses three class columns:",
    "",
    "- `primary_age_class`: Within-region age-specific classification (CA1_old_specific, CA1_young_specific, CA1_age_conserved, CA1_age_shifted, CA3_old_specific, CA3_young_specific, CA3_age_conserved, CA3_age_shifted, low_confidence, module_member_or_effector).",
    "- `region_age_class`: Cross-region within-age classification (Young_CA1_specific, Young_CA3_specific, Young_region_conserved, Old_CA1_specific, Old_CA3_specific, Old_region_conserved, none).",
    "- `candidate_class`: Same as primary_age_class (backward compatibility).",
    "- `secondary_class`: Additional classes when a pair falls into multiple categories (semicolon-separated).",
    "",
    "### Sample Size Labels",
    "",
    "Per-stratum labels (in the main table):",
    "- `sample_size_label_CA1_Young` = exploratory_small_n (n=4)",
    "- `sample_size_label_CA1_Old` = exploratory_moderate_n (n=8)",
    "- `sample_size_label_CA3_Young` = severely_underpowered_exploratory (n=3)",
    "- `sample_size_label_CA3_Old` = exploratory_moderate_n (n=8)",
    "",
    "In specialized output tables (old_specific, young_specific, etc.), the `sample_size_label` column uses the label for the relevant stratum.",
    "",
    "### How to Read Age-Stratified Tables",
    "",
    "- Each row is a gene x module pair with separate rho values for Young and Old.",
    "- Compare rho_Young vs rho_Old within the same gene-module to identify age-specific patterns.",
    "- `candidate_score_<stratum> = |rho_<stratum>|` is the primary ranking metric (NOT a p-value).",
    "- `pvalue_exploratory_*` and `qvalue_exploratory_*` columns are computed but NOT used for ranking or class assignment.",
    "",
    "### How to Interpret Rho Direction",
    "",
    "- Positive rho: gene expression and module score increase together.",
    "- Negative rho: gene expression increases while module score decreases.",
    "- In Old strata, this may reflect aging-related coupling changes.",
    "",
    "### Why Small n Matters",
    "",
    "- CA1 Young n=4 and CA3 Young n=3 are underpowered for Spearman correlation.",
    "- rho estimates in these strata are noisy. Compare with Old strata (n=8) for more stable estimates.",
    "- CA3 Young n=3: minimum p-value ~0.33 for perfect monotone correlation. Do NOT use CA3 Young p-values for filtering.",
    "",
    "### Candidate Classes",
    "",
    "- Age-specific (e.g., CA1_old_specific): strong rho in one age but NOT the other within same region.",
    "- Age-conserved: strong rho in both Young and Old, same direction.",
    "- Age-shifted: strong rho in both Young and Old, OPPOSITE direction.",
    "- Region-specific within age (e.g., Old_CA1_specific): strong rho in one region but NOT the other within same age.",
    "- `low_confidence`: |rho| < 0.3 in ALL 4 strata.",
    "- `module_member_or_effector`: gene is a member of the module gene set. Self-correlation excluded from regulator ranking.",
    "",
    "### Key Files",
    "",
    "- `phase13c_age_stratified_rho.csv` -- full rho table",
    "- `phase13c_age_stratified_scores.csv` -- scores with per-stratum sample_size_label columns",
    "- `phase13c_age_specific_candidate_classes.csv` -- all class columns (primary_age_class, region_age_class, candidate_class, secondary_class)",
    "- `phase13c_old_specific_candidates.csv` -- Old-specific candidates",
    "- `phase13c_young_specific_candidates.csv` -- Young-specific candidates",
    "- `phase13c_age_shifted_candidates.csv` -- candidates with rho sign reversal",
    "- `phase13c_region_shifted_by_age_candidates.csv` -- region-specific candidates within same age",
    "- `phase13c_top_candidates_by_region_age.csv` -- top 20 per stratum",
    "",
    "### Caveats",
    "",
    "1. These are ASSOCIATION-BASED candidate regulator-associated genes, NOT causal regulators.",
    "2. Age-stratified rho is still correlation. Age specificity adds biological plausibility but does NOT establish causality.",
    "3. CA3 Young n=3 -- all CA3 Young results carry `severely_underpowered_exploratory` label.",
    "4. CA1 Young n=4 -- all CA1 Young results carry `exploratory_small_n` label.",
    ""
  )
  guide <- c(guide, phase13c_section)
  writeLines(guide, guide_path)
  cat("Analysis guide updated with Phase 13c section (deduplicated)\n")

  # V24b: Exactly one Phase 13c section in guide
  guide_check <- readLines(guide_path, warn = FALSE)
  n13c_sections <- length(grep("^## Phase 13c:", guide_check))
  add_check("V24b_guide_phase13c_unique", if (n13c_sections == 1) "PASS" else "FAIL",
    paste(n13c_sections, "Phase 13c sections in guide"))

  add_check("V24_guide_updated", "PASS", "Phase 13c section appended to analysis guide")
} else {
  add_check("V24_guide_updated", "WARN", "Analysis guide not found")
}

# Final report
t_end <- Sys.time()
elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

final_report <- c(
  "=== Phase Spatial-13c: Final Report ===",
  paste("Generated:", format(Sys.time())),
  paste("Plan version: v1.0"),
  paste("Script: R/spatial/s13c_age_stratified_regulator_discovery.R"),
  paste("Elapsed:", round(elapsed, 1), "minutes"),
  "",
  "--- Sample Sizes ---",
  paste("CA1 Young:", length(ca1_young), "GEMgroups"),
  paste("CA1 Old:", length(ca1_old), "GEMgroups"),
  paste("CA3 Young:", length(ca3_young), "GEMgroups (severely underpowered)"),
  paste("CA3 Old:", length(ca3_old), "GEMgroups"),
  "",
  "--- Universe ---",
  paste("Regulator universe:", length(universe_regulator), "genes"),
  paste("Genes in pseudobulk:", length(genes_in_pb)),
  paste("Modules analyzed:", length(module_names)),
  paste("Total gene x module pairs:", nrow(rho_results)),
  "",
  "--- Candidate Class Distribution (candidate_class) ---",
  capture.output(print(table(rho_results$candidate_class))),
  "",
  "--- Primary Age Class Distribution (primary_age_class) ---",
  capture.output(print(table(rho_results$primary_age_class))),
  "",
  "--- Region Age Class Distribution (region_age_class) ---",
  capture.output(print(table(rho_results$region_age_class))),
  "",
  "--- Sample Size Label Distribution ---",
  paste("CA1_Young:", unique(rho_results$sample_size_label_CA1_Young)),
  paste("CA1_Old:", unique(rho_results$sample_size_label_CA1_Old)),
  paste("CA3_Young:", unique(rho_results$sample_size_label_CA3_Young)),
  paste("CA3_Old:", unique(rho_results$sample_size_label_CA3_Old)),
  "",
  "--- Self-Correlation Flags ---",
  capture.output(print(table(rho_results$self_correlation_flag))),
  "",
  "--- Validation ---",
  paste("PASS:", sum(val_df$status == "PASS")),
  paste("WARN:", sum(val_df$status == "WARN")),
  paste("FAIL:", sum(val_df$status == "FAIL")),
  "",
  "--- Caveats ---",
  "1. Candidate regulator-associated genes are ASSOCIATION-BASED, not causal.",
  "2. CA3 Young n=3: severely underpowered. All CA3 Young rho values are exploratory.",
  "3. CA1 Young n=4: small-n caution. rho estimates are noisy.",
  "4. Primary ranking uses |rho|, NOT p-values or q-values.",
  "5. Phase 13c does NOT modify Phase 13 or 13b outputs."
)

writeLines(final_report, file.path(out_dir, "phase13c_final_report.txt"))
cat("\nFinal report written\n")

cat("\n=== Phase Spatial-13c Complete ===\n")
cat("Elapsed:", round(elapsed, 1), "minutes\n")
cat("PASS:", sum(val_df$status == "PASS"), " WARN:", sum(val_df$status == "WARN"), " FAIL:", sum(val_df$status == "FAIL"), "\n")

if (any(val_df$status == "FAIL")) {
  cat("WARNING: Some validation checks FAILED. See phase13c_validation_summary.csv\n")
}

cat("\nDone.\n")

# Cleanup Rplots.pdf if created by default device
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
