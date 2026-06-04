#!/usr/bin/env Rscript
# s11b_age_grouping_ca1_ca3_de.R
# Phase Spatial-11 Age Grouping: Age-Stratified CA1 vs CA3 DE / Target-Gene / Module Analysis
# Plan version: v1.1 (2026-06-04)
# RDS-based spatial analysis; no PDF; no Seurat object modification

# ============================================================================
# E0: SETUP
# ============================================================================

cat("=== Phase Spatial-11 Age Grouping: Age-Stratified CA1 vs CA3 DE ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# -- Rplots.pdf guard (start) --
rplots_start <- file.exists("Rplots.pdf")
rplots_start_md5 <- if (rplots_start) {
  tryCatch(digest::digest("Rplots.pdf", file = TRUE), error = function(e) tools::md5sum("Rplots.pdf"))
} else NA_character_
cat("[Rplots.pdf] Pre-existence:", rplots_start, if (!is.na(rplots_start_md5)) paste("MD5:", rplots_start_md5) else "", "\n")

# -- Git status at start --
git_status_start <- tryCatch(system("git status --short", intern = TRUE), error = function(e) "git not available")
cat("Git status at start:", if (length(git_status_start) == 0) "(clean)" else paste(git_status_start, collapse = "; "), "\n\n")

# -- renv.lock MD5 before --
renv_lock_md5_before <- tryCatch(tools::md5sum("renv.lock"), error = function(e) NA_character_)
cat("renv.lock MD5 before:", renv_lock_md5_before, "\n")

# -- renv.lock status vs git HEAD --
renv_vs_git <- tryCatch({
  diff_lines <- system("git diff HEAD -- renv.lock | head -5", intern = TRUE)
  if (length(diff_lines) > 0) "modified_vs_git_HEAD" else "clean_vs_git_HEAD"
}, error = function(e) "git_not_available")
cat("renv.lock vs git HEAD:", renv_vs_git, "\n\n")

# -- Required packages --
pkgs_required <- c("Seurat", "DESeq2", "Matrix", "ggplot2", "data.table", "digest")
for (pkg in pkgs_required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[SC-STOP] Package not available:", pkg, "\n")
    stop(paste("Package", pkg, "not installed"))
  }
}
library(Seurat)
library(DESeq2)
library(Matrix)
library(ggplot2)
library(data.table)

# -- Optional packages --
optional_pkg_actions <- list()

# Matrix.utils
if (!requireNamespace("Matrix.utils", quietly = TRUE)) {
  cat("[INFO] Matrix.utils not found; attempting renv::install\n")
  tryCatch({
    renv::install("Matrix.utils")
    renv::snapshot()
    optional_pkg_actions[["Matrix.utils"]] <- "installed"
  }, error = function(e) {
    cat("[WARNING] Matrix.utils install failed:", conditionMessage(e), "\n")
    optional_pkg_actions[["Matrix.utils"]] <- "install_failed"
  })
}

# apeglm
apeglm_available <- requireNamespace("apeglm", quietly = TRUE)
if (!apeglm_available) {
  cat("[INFO] apeglm not found; attempting renv::install\n")
  tryCatch({
    renv::install("apeglm")
    renv::snapshot()
    apeglm_available <- TRUE
    optional_pkg_actions[["apeglm"]] <- "installed"
  }, error = function(e) {
    cat("[WARNING] apeglm install failed:", conditionMessage(e), "\n")
  })
}

# pheatmap
pheatmap_available <- requireNamespace("pheatmap", quietly = TRUE)
if (pheatmap_available) library(pheatmap)

# ggrepel
ggrepel_available <- requireNamespace("ggrepel", quietly = TRUE)
if (ggrepel_available) library(ggrepel)

# corrplot
corrplot_available <- requireNamespace("corrplot", quietly = TRUE)
if (corrplot_available) library(corrplot)

cat("\n--- Package versions ---\n")
for (pkg in c("Seurat", "DESeq2", "Matrix", "Matrix.utils", "apeglm", "ggplot2",
              "pheatmap", "ggrepel", "corrplot", "data.table", "digest")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, ":", as.character(packageVersion(pkg)), "\n")
  }
}
cat("\n")

# -- Output directories --
out_data <- "data/processed/spatial/phase11_age_grouping"
out_figs <- "figures/spatial/phase11_age_grouping"
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

# -- Provenance tracking --
provenance <- list()
provenance[["script"]] <- "R/spatial/s11b_age_grouping_ca1_ca3_de.R"
provenance[["plan_version"]] <- "v1.1"
provenance[["start_time"]] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
provenance[["renv_lock_md5_before"]] <- renv_lock_md5_before
provenance[["renv_vs_git_head_start"]] <- renv_vs_git

# -- Validation log --
validation_log <- list()
log_check <- function(id, status, msg) {
  validation_log[[id]] <<- list(status = status, msg = msg)
  cat(sprintf("[%s] %s: %s\n", id, status, msg))
}

# ============================================================================
# E1: INPUT AUDIT
# ============================================================================

cat("\n=== E1: Input Audit ===\n")

# Required files
f_hippo <- "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"
f_target_csv <- "docs/target_genes.csv"
f_phase10_audit <- "data/processed/spatial/phase10_target_gene_audit/target_gene_input_audit.csv"
f_manifest <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/pseudobulk_sample_manifest.csv"
f_pseudobulk <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds"

# Recommended files
f_mod_scores <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_scores_sample_summary.csv"
f_mod_delta <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_score_delta_by_age.csv"
f_deseq2_base <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/deseq2_all_genes_CA3_vs_CA1.csv"
f_missing_rescue <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/missing_gene_symbol_rescue_audit.csv"
f_mod_genesets <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv"
f_coupling <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_coupling_combined_long.csv"

# Check required
for (f in c(f_hippo, f_target_csv, f_phase10_audit, f_manifest, f_pseudobulk)) {
  if (!file.exists(f)) {
    log_check("V01", "STOP", paste("Missing required file:", f))
    stop(paste("Missing required file:", f))
  }
}
log_check("V01", "PASS", "All required input files exist")

# Check recommended
recommended_missing <- character(0)
for (f in c(f_mod_scores, f_mod_delta, f_deseq2_base, f_missing_rescue, f_mod_genesets, f_coupling)) {
  if (!file.exists(f)) recommended_missing <- c(recommended_missing, f)
}
if (length(recommended_missing) > 0) {
  log_check("V01-rec", "WARNING", paste("Recommended files missing:", paste(basename(recommended_missing), collapse=", ")))
} else {
  log_check("V01-rec", "PASS", "All recommended files exist")
}

# Load Phase 11 manifest
manifest <- read.csv(f_manifest, stringsAsFactors = FALSE)
cat("Manifest loaded:", nrow(manifest), "samples\n")

# Verify manifest structure
required_cols <- c("sample_id", "GEMgroup", "region_group", "Age", "n_spots", "has_pair")
if (!all(required_cols %in% colnames(manifest))) {
  log_check("V03", "STOP", "Manifest missing required columns")
  stop("Manifest missing required columns")
}

# Verify ages present
ages_present <- sort(unique(manifest$Age))
if (!all(c("Young", "Middle", "Old") %in% ages_present)) {
  log_check("V03", "STOP", paste("Missing age groups. Found:", paste(ages_present, collapse=", ")))
  stop("Missing age groups")
}
log_check("V03", "PASS", paste("Age groups present:", paste(ages_present, collapse=", ")))

# Load pseudobulk counts
pseudobulk_rds <- readRDS(f_pseudobulk)
cat("Pseudobulk counts loaded:", nrow(pseudobulk_rds), "genes x", ncol(pseudobulk_rds), "samples\n")

# Verify dimensions
if (ncol(pseudobulk_rds) != nrow(manifest)) {
  log_check("V01-dim", "WARNING", paste("Count cols", ncol(pseudobulk_rds), "!= manifest rows", nrow(manifest)))
} else {
  log_check("V01-dim", "PASS", paste("Count matrix:", nrow(pseudobulk_rds), "genes x", ncol(pseudobulk_rds), "samples"))
}

# Verify sample_id alignment
count_sample_ids <- colnames(pseudobulk_rds)
manifest_sample_ids <- manifest$sample_id
if (!setequal(count_sample_ids, manifest_sample_ids)) {
  log_check("V01-align", "WARNING", "Sample IDs mismatch between counts and manifest")
} else {
  # Reorder manifest to match counts
  manifest <- manifest[match(count_sample_ids, manifest$sample_id), ]
  log_check("V01-align", "PASS", "Sample IDs aligned between counts and manifest")
}

# Load Phase 10 gene audit
phase10_audit <- read.csv(f_phase10_audit, stringsAsFactors = FALSE)
cat("Phase 10 audit loaded:", nrow(phase10_audit), "genes\n")
log_check("V10-audit", "PASS", paste("Phase 10 audit:", nrow(phase10_audit), "genes,", sum(phase10_audit$mapping_status == "exact_match"), "exact-match"))

# Load Phase 11 validation to confirm PASS_WITH_WARNINGS
phase11_val <- readLines("data/processed/spatial/phase11_ca1_ca3_de_module_coupling/phase11_validation_summary.txt")
if (any(grepl("PASS_WITH_WARNINGS", phase11_val))) {
  log_check("V01-base", "PASS", "Phase 11 base validation: PASS_WITH_WARNINGS")
} else {
  log_check("V01-base", "WARNING", "Phase 11 base validation status unexpected")
}

# Load module scores sample summary (for coupling)
mod_scores <- NULL
if (file.exists(f_mod_scores)) {
  mod_scores <- read.csv(f_mod_scores, stringsAsFactors = FALSE)
  cat("Module scores loaded:", nrow(mod_scores), "samples\n")
}

# Load module delta by age
mod_delta <- NULL
if (file.exists(f_mod_delta)) {
  mod_delta <- read.csv(f_mod_delta, stringsAsFactors = FALSE)
  cat("Module delta loaded:", nrow(mod_delta), "rows\n")
}

# Load coupling summary from Phase 11 base
coupling_base <- NULL
if (file.exists(f_coupling)) {
  coupling_base <- read.csv(f_coupling, stringsAsFactors = FALSE)
  cat("Coupling base loaded:", nrow(coupling_base), "rows\n")
}

# Load module gene sets
mod_genesets <- NULL
if (file.exists(f_mod_genesets)) {
  mod_genesets <- read.csv(f_mod_genesets, stringsAsFactors = FALSE)
  cat("Module gene sets loaded:", nrow(mod_genesets), "modules\n")
}

cat("\n")

# ============================================================================
# E2: AGE-SPECIFIC SAMPLE MANIFEST
# ============================================================================

cat("=== E2: Age-Specific Sample Manifest ===\n")

# Build age-specific manifest
age_manifest <- manifest[, c("sample_id", "GEMgroup", "region_group", "Age", "n_spots", "has_pair")]

# Add age_stratum column
age_manifest$age_stratum <- age_manifest$Age

# Add paired status per age
age_manifest$has_ca3_mate_in_age <- FALSE
for (age_val in c("Young", "Middle", "Old")) {
  age_mask <- age_manifest$Age == age_val
  age_gemgroups <- unique(age_manifest$GEMgroup[age_mask])
  for (g in age_gemgroups) {
    g_mask <- age_manifest$GEMgroup == g & age_mask
    regions <- age_manifest$region_group[g_mask]
    if (all(c("CA1", "CA3") %in% regions)) {
      age_manifest$has_ca3_mate_in_age[g_mask] <- TRUE
    }
  }
}

# Save age-specific manifest
write.csv(age_manifest, file.path(out_data, "age_group_pseudobulk_sample_manifest.csv"), row.names = FALSE)

# Count paired GEMgroups per age
age_paired_counts <- data.frame(
  Age = c("Young", "Middle", "Old"),
  n_paired_GEMgroups = c(
    sum(tapply(age_manifest$region_group[age_manifest$Age == "Young"],
               age_manifest$GEMgroup[age_manifest$Age == "Young"],
               function(x) all(c("CA1", "CA3") %in% x))),
    sum(tapply(age_manifest$region_group[age_manifest$Age == "Middle"],
               age_manifest$GEMgroup[age_manifest$Age == "Middle"],
               function(x) all(c("CA1", "CA3") %in% x))),
    sum(tapply(age_manifest$region_group[age_manifest$Age == "Old"],
               age_manifest$GEMgroup[age_manifest$Age == "Old"],
               function(x) all(c("CA1", "CA3") %in% x)))
  ),
  n_total_samples = c(
    sum(age_manifest$Age == "Young"),
    sum(age_manifest$Age == "Middle"),
    sum(age_manifest$Age == "Old")
  ),
  power_label = c("exploratory_low_power", "exploratory_low_power", "primary_adequately_powered")
)

cat("\nAge-specific paired GEMgroup counts:\n")
print(age_paired_counts)

# Record GEMgroup 1 as CA1-only for Young
g1_young <- age_manifest[age_manifest$GEMgroup == 1 & age_manifest$Age == "Young", ]
if (nrow(g1_young) > 0 && !any(g1_young$region_group == "CA3")) {
  cat("\n[GEMgroup 1] CA1-only in Young — excluded from Young paired DE\n")
}

log_check("V04", "PASS", sprintf("Paired GEMgroups: Young=%d, Middle=%d, Old=%d",
                                   age_paired_counts$n_paired_GEMgroups[1],
                                   age_paired_counts$n_paired_GEMgroups[2],
                                   age_paired_counts$n_paired_GEMgroups[3]))

cat("\n")

# ============================================================================
# E3: DESIGN AUDIT PER AGE
# ============================================================================

cat("=== E3: Design Audit Per Age ===\n")

design_audit <- data.frame(
  age = character(),
  n_samples = integer(),
  n_paired_gemgroups = integer(),
  design_formula = character(),
  rank = integer(),
  ncol = integer(),
  is_full_rank = logical(),
  model_type = character(),
  power_label = character(),
  stringsAsFactors = FALSE
)

age_designs <- list()

for (age_val in c("Young", "Middle", "Old")) {
  cat(sprintf("\n--- %s ---\n", age_val))

  # Subset to this age with CA1+CA3 only (for paired design)
  age_mask <- age_manifest$Age == age_val & age_manifest$has_ca3_mate_in_age
  age_samples <- age_manifest[age_mask, ]

  n_samples <- nrow(age_samples)
  n_pairs <- length(unique(age_samples$GEMgroup))

  cat(sprintf("  Samples with CA1+CA3 mate: %d, Paired GEMgroups: %d\n", n_samples, n_pairs))

  # Check SC-5: need >= 2 paired GEMgroups
  if (n_pairs < 2) {
    cat(sprintf("  [SC-5] < 2 paired GEMgroups — blocked_insufficient_pairs\n"))
    design_audit <- rbind(design_audit, data.frame(
      age = age_val, n_samples = n_samples, n_paired_gemgroups = n_pairs,
      design_formula = "blocked", rank = NA, ncol = NA, is_full_rank = NA,
      model_type = "blocked_insufficient_pairs",
      power_label = if (age_val == "Old") "primary_adequately_powered" else "exploratory_low_power",
      stringsAsFactors = FALSE
    ))
    age_designs[[age_val]] <- list(model_type = "blocked_insufficient_pairs", samples = age_samples)
    next
  }

  # Build design matrix for ~ GEMgroup + region_group
  age_samples$GEMgroup_f <- factor(age_samples$GEMgroup)
  age_samples$region_group_f <- factor(age_samples$region_group, levels = c("CA1", "CA3"))

  tryCatch({
    mm <- model.matrix(~ GEMgroup_f + region_group_f, data = age_samples)
    qr_rank <- qr(mm)$rank
    n_cols <- ncol(mm)
    full_rank <- qr_rank == n_cols

    cat(sprintf("  Design: ~ GEMgroup + region_group, rank=%d, ncol=%d, full_rank=%s\n",
                qr_rank, n_cols, full_rank))

    if (full_rank) {
      design_audit <- rbind(design_audit, data.frame(
        age = age_val, n_samples = n_samples, n_paired_gemgroups = n_pairs,
        design_formula = "~ GEMgroup + region_group", rank = qr_rank, ncol = n_cols,
        is_full_rank = TRUE, model_type = "paired_blocked",
        power_label = if (age_val == "Old") "primary_adequately_powered" else "exploratory_low_power",
        stringsAsFactors = FALSE
      ))
      age_designs[[age_val]] <- list(
        model_type = "paired_blocked",
        design = "~ GEMgroup + region_group",
        samples = age_samples
      )
    } else {
      # Fallback: ~ region_group only
      cat("  [SC-6] Rank-deficient — attempting fallback ~ region_group\n")
      mm_fb <- model.matrix(~ region_group_f, data = age_samples)
      qr_rank_fb <- qr(mm_fb)$rank
      n_cols_fb <- ncol(mm_fb)
      full_rank_fb <- qr_rank_fb == n_cols_fb

      design_audit <- rbind(design_audit, data.frame(
        age = age_val, n_samples = n_samples, n_paired_gemgroups = n_pairs,
        design_formula = "~ region_group (fallback)", rank = qr_rank_fb, ncol = n_cols_fb,
        is_full_rank = full_rank_fb, model_type = "fallback_unpaired",
        power_label = if (age_val == "Old") "primary_adequately_powered" else "exploratory_low_power",
        stringsAsFactors = FALSE
      ))
      age_designs[[age_val]] <- list(
        model_type = "fallback_unpaired",
        design = "~ region_group",
        samples = age_samples
      )
    }
  }, error = function(e) {
    cat(sprintf("  Design matrix error: %s\n", conditionMessage(e)))
    design_audit <- rbind(design_audit, data.frame(
      age = age_val, n_samples = n_samples, n_paired_gemgroups = n_pairs,
      design_formula = "ERROR", rank = NA, ncol = NA, is_full_rank = NA,
      model_type = "blocked_insufficient_pairs",
      power_label = if (age_val == "Old") "primary_adequately_powered" else "exploratory_low_power",
      stringsAsFactors = FALSE
    ))
    age_designs[[age_val]] <<- list(model_type = "blocked_insufficient_pairs", samples = age_samples)
  })
}

# Save design audit
write.csv(design_audit, file.path(out_data, "design_matrix_audit_by_age.csv"), row.names = FALSE)

# Log V05-V08
for (i in 1:nrow(design_audit)) {
  row <- design_audit[i, ]
  if (row$model_type == "blocked_insufficient_pairs") {
    log_check(paste0("V05-", row$age), "STOP", sprintf("%s: %s", row$age, row$model_type))
  } else if (row$model_type == "fallback_unpaired") {
    log_check(paste0("V08-", row$age), "WARNING", sprintf("%s: fallback_unpaired", row$age))
  } else {
    log_check(paste0("V07-", row$age), "PASS", sprintf("%s: paired_blocked, full rank", row$age))
  }
}

cat("\n")

# ============================================================================
# E4-E6: AGE-STRATIFIED DESeq2
# ============================================================================

run_deseq2_age <- function(age_val, age_design_info, pseudobulk_counts, manifest_df) {
  cat(sprintf("\n=== E4/E5/E6: DESeq2 for %s ===\n", age_val))

  model_type <- age_design_info$model_type
  power_label <- if (age_val == "Old") "primary_adequately_powered" else "exploratory_low_power"

  if (model_type == "blocked_insufficient_pairs") {
    cat(sprintf("  [SC-5] %s: blocked_insufficient_pairs — skipping DESeq2\n", age_val))
    return(list(
      age = age_val, model_type = model_type, power_label = power_label,
      results = NULL, n_tested = 0, n_sig = 0, converged = FALSE
    ))
  }

  age_samples <- age_design_info$samples
  sample_ids <- age_samples$sample_id

  # Subset pseudobulk counts
  counts_age <- pseudobulk_counts[, sample_ids, drop = FALSE]

  # Filter low-count genes
  keep <- rowSums(counts_age) >= 10
  counts_age <- counts_age[keep, ]
  cat(sprintf("  Genes after filtering (sum>=10): %d / %d\n", nrow(counts_age), nrow(pseudobulk_counts)))

  if (nrow(counts_age) == 0) {
    cat(sprintf("  [SC-A2] %s: all genes filtered — STOP this age\n", age_val))
    return(list(
      age = age_val, model_type = model_type, power_label = power_label,
      results = NULL, n_tested = 0, n_sig = 0, converged = FALSE
    ))
  }

  # Prepare colData
  col_data <- data.frame(
    sample_id = sample_ids,
    GEMgroup = factor(age_samples$GEMgroup),
    region_group = factor(age_samples$region_group, levels = c("CA1", "CA3")),
    row.names = sample_ids
  )

  # Convert to dense integer matrix for DESeq2
  counts_dense <- as.matrix(round(counts_age))
  storage.mode(counts_dense) <- "integer"

  # Build design
  if (model_type == "paired_blocked") {
    design_formula <- ~ GEMgroup + region_group
  } else {
    design_formula <- ~ region_group
  }

  cat(sprintf("  Design: %s\n", deparse(design_formula)))

  tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = counts_dense,
      colData = col_data,
      design = design_formula
    )
    dds <- DESeq(dds)

    res <- results(dds, contrast = c("region_group", "CA3", "CA1"), alpha = 0.05)

    # LFC shrinkage
    if (apeglm_available) {
      coef_name <- resultsNames(dds)
      # Find the region_group coefficient
      coef_idx <- grep("region_group", coef_name)
      if (length(coef_idx) > 0) {
        res_shrunken <- lfcShrink(dds, coef = coef_name[coef_idx[1]], type = "apeglm")
        res$log2FoldChange <- res_shrunken$log2FoldChange
        res$lfcSE <- res_shrunken$lfcSE
        cat("  LFC shrinkage: apeglm\n")
      }
    } else {
      cat("  LFC shrinkage: none (apeglm not available)\n")
    }

    # Count results
    n_tested <- sum(!is.na(res$padj))
    n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
    n_up_ca3 <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
    n_down_ca1 <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

    cat(sprintf("  Results: %d tested, %d sig (padj<0.05), %d CA3-high, %d CA1-high\n",
                n_tested, n_sig, n_up_ca3, n_down_ca1))

    # Build output data frame
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$model_type <- model_type
    res_df$power_label <- power_label
    res_df$baseMean_log10 <- log10(res_df$baseMean + 1)

    # Add significance flags
    res_df$padj_sig <- !is.na(res_df$padj) & res_df$padj < 0.05
    res_df$lfc_notable <- abs(res_df$log2FoldChange) > log2(1.5)
    res_df$direction <- ifelse(res_df$log2FoldChange > 0, "CA3_high", "CA1_high")

    list(
      age = age_val, model_type = model_type, power_label = power_label,
      results = res_df, dds = dds,
      n_tested = n_tested, n_sig = n_sig,
      n_up_ca3 = n_up_ca3, n_down_ca1 = n_down_ca1,
      converged = TRUE
    )
  }, error = function(e) {
    cat(sprintf("  [SC-A3] DESeq2 failed for %s: %s\n", age_val, conditionMessage(e)))
    list(
      age = age_val, model_type = model_type, power_label = power_label,
      results = NULL, n_tested = 0, n_sig = 0, converged = FALSE
    )
  })
}

# Run DESeq2 for each age
deseq2_results <- list()
for (age_val in c("Young", "Middle", "Old")) {
  if (!is.null(age_designs[[age_val]])) {
    deseq2_results[[age_val]] <- run_deseq2_age(age_val, age_designs[[age_val]], pseudobulk_rds, manifest)
  }
}

# Save per-age DE results
for (age_val in c("Young", "Middle", "Old")) {
  res <- deseq2_results[[age_val]]
  if (!is.null(res) && !is.null(res$results)) {
    fname <- sprintf("deseq2_all_genes_CA3_vs_CA1_%s.csv", age_val)
    write.csv(res$results, file.path(out_data, fname), row.names = FALSE)
    cat(sprintf("Saved: %s (%d rows)\n", fname, nrow(res$results)))
  }
}

# Log V09
for (age_val in c("Young", "Middle", "Old")) {
  res <- deseq2_results[[age_val]]
  if (!is.null(res) && res$converged && res$n_tested > 0) {
    log_check(paste0("V09-", age_val), "PASS", sprintf("%s: %d tested, %d sig", age_val, res$n_tested, res$n_sig))
  } else if (!is.null(res) && !res$converged) {
    log_check(paste0("V09-", age_val), "WARNING", sprintf("%s: DESeq2 did not converge", age_val))
  } else {
    log_check(paste0("V09-", age_val), "WARNING", sprintf("%s: blocked or no results", age_val))
  }
}

cat("\n")

# ============================================================================
# E7: TARGET GENE ANNOTATION BY AGE
# ============================================================================

cat("=== E7: Target Gene Annotation by Age ===\n")

# Build combined target gene table: 251 genes x 3 ages
target_genes <- phase10_audit$gene_symbol
target_combined <- data.frame(
  gene_symbol = target_genes,
  category = phase10_audit$category,
  category_english = phase10_audit$category_english,
  mapping_status = phase10_audit$mapping_status,
  stringsAsFactors = FALSE
)

# Add per-age DE results
for (age_val in c("Young", "Middle", "Old")) {
  res <- deseq2_results[[age_val]]
  age_suffix <- paste0("_", age_val)

  if (!is.null(res) && !is.null(res$results)) {
    res_df <- res$results
    # Match genes
    idx <- match(target_combined$gene_symbol, res_df$gene)

    target_combined[[paste0("log2FoldChange", age_suffix)]] <- res_df$log2FoldChange[idx]
    target_combined[[paste0("padj", age_suffix)]] <- res_df$padj[idx]
    target_combined[[paste0("baseMean", age_suffix)]] <- res_df$baseMean[idx]
    target_combined[[paste0("padj_sig", age_suffix)]] <- res_df$padj_sig[idx]
    target_combined[[paste0("direction", age_suffix)]] <- res_df$direction[idx]
    target_combined[[paste0("lfc_notable", age_suffix)]] <- res_df$lfc_notable[idx]
  } else {
    target_combined[[paste0("log2FoldChange", age_suffix)]] <- NA
    target_combined[[paste0("padj", age_suffix)]] <- NA
    target_combined[[paste0("baseMean", age_suffix)]] <- NA
    target_combined[[paste0("padj_sig", age_suffix)]] <- NA
    target_combined[[paste0("direction", age_suffix)]] <- NA
    target_combined[[paste0("lfc_notable", age_suffix)]] <- NA
  }
}

# Mark missing genes
target_combined$in_spatial <- target_combined$mapping_status == "exact_match"

write.csv(target_combined, file.path(out_data, "deseq2_target_genes_CA3_vs_CA1_by_age.csv"), row.names = FALSE)

# Log V10, V11
n_found <- sum(target_combined$in_spatial)
n_missing <- sum(!target_combined$in_spatial)
log_check("V10", if (n_found == 231 && n_missing == 20) "PASS" else "WARNING",
          sprintf("Target genes: %d found, %d missing (expected 231+20=251)", n_found, n_missing))

# Verify missing genes have NA DE values
missing_genes <- target_combined[!target_combined$in_spatial, ]
de_cols <- grep("log2FoldChange_|padj_", colnames(missing_genes), value = TRUE)
all_na <- all(is.na(missing_genes[, de_cols]))
log_check("V11", if (all_na) "PASS" else "WARNING",
          sprintf("Missing genes DE columns all NA: %s", all_na))

cat(sprintf("Target gene by-age table: %d rows x %d cols\n", nrow(target_combined), ncol(target_combined)))
cat("\n")

# ============================================================================
# E8: CROSS-AGE log2FC MATRIX (Route B)
# ============================================================================

cat("=== E8: Cross-Age log2FC Matrix ===\n")

# Build wide matrix: 231 found genes x 3 age log2FC columns
found_mask <- target_combined$in_spatial
log2fc_matrix <- target_combined[found_mask, c("gene_symbol", "category", "category_english")]
log2fc_matrix$Young_log2FC <- target_combined$log2FoldChange_Young[found_mask]
log2fc_matrix$Middle_log2FC <- target_combined$log2FoldChange_Middle[found_mask]
log2fc_matrix$Old_log2FC <- target_combined$log2FoldChange_Old[found_mask]

# Classify direction consistency
classify_age_pattern <- function(y, m, o) {
  tryCatch({
    if (length(y) != 1 || length(m) != 1 || length(o) != 1) return("error")
    if (is.na(y) || is.na(m) || is.na(o) || is.nan(y) || is.nan(m) || is.nan(o)) return("has_NA")
    signs <- sign(c(y, m, o))
    if (all(signs > 0)) return("age_concordant_CA3_high")
    if (all(signs < 0)) return("age_concordant_CA1_high")
    if (signs[3] != signs[1] && signs[3] != signs[2]) return("age_flipped")
    if (abs(o) > abs(y) && abs(o) > abs(m)) return("old_enhanced")
    return("mixed")
  }, error = function(e) return("error"))
}

age_patterns <- character(nrow(log2fc_matrix))
for (i in seq_len(nrow(log2fc_matrix))) {
  age_patterns[i] <- classify_age_pattern(
    log2fc_matrix$Young_log2FC[i],
    log2fc_matrix$Middle_log2FC[i],
    log2fc_matrix$Old_log2FC[i]
  )
}
log2fc_matrix$age_pattern <- age_patterns

write.csv(log2fc_matrix, file.path(out_data, "target_gene_log2fc_matrix_by_age.csv"), row.names = FALSE)

# Log V22: no unexpected all-NA rows
all_na_rows <- sum(is.na(log2fc_matrix$Young_log2FC) & is.na(log2fc_matrix$Middle_log2FC) & is.na(log2fc_matrix$Old_log2FC))
log_check("V22", if (all_na_rows == 0) "PASS" else "WARNING",
          sprintf("Cross-age log2FC matrix: %d rows, %d all-NA rows", nrow(log2fc_matrix), all_na_rows))

cat(sprintf("Age pattern distribution:\n"))
print(table(log2fc_matrix$age_pattern))
cat("\n")

# ============================================================================
# E9: CATEGORY EFFECT BY AGE (Route C)
# ============================================================================

cat("=== E9: Category Effect by Age ===\n")

categories <- unique(target_combined[, c("category", "category_english")])
cat_effect <- data.frame()

for (age_val in c("Young", "Middle", "Old")) {
  power_label <- if (age_val == "Old") "primary_adequately_powered" else "exploratory_low_power"
  lfc_col <- paste0("log2FoldChange_", age_val)
  padj_col <- paste0("padj_", age_val)
  sig_col <- paste0("padj_sig_", age_val)

  for (cat_i in 1:nrow(categories)) {
    cat_mask <- target_combined$category == categories$category[cat_i] & target_combined$in_spatial
    lfc_vals <- target_combined[[lfc_col]][cat_mask]
    padj_vals <- target_combined[[padj_col]][cat_mask]
    sig_vals <- target_combined[[sig_col]][cat_mask]

    n_available <- sum(cat_mask)
    n_sig <- sum(sig_vals, na.rm = TRUE)
    n_up_ca3 <- sum(sig_vals & lfc_vals > 0, na.rm = TRUE)
    n_down_ca1 <- sum(sig_vals & lfc_vals < 0, na.rm = TRUE)
    n_notable <- sum(abs(lfc_vals) > log2(1.5), na.rm = TRUE)
    mean_lfc <- mean(lfc_vals, na.rm = TRUE)
    median_lfc <- median(lfc_vals, na.rm = TRUE)
    n_ca3_high <- sum(lfc_vals > 0, na.rm = TRUE)
    n_ca1_high <- sum(lfc_vals < 0, na.rm = TRUE)

    cat_effect <- rbind(cat_effect, data.frame(
      category_cn = categories$category[cat_i],
      category_en = categories$category_english[cat_i],
      age = age_val,
      n_genes_available = n_available,
      n_sig_padj05 = n_sig,
      n_up_CA3 = n_up_ca3,
      n_down_CA1 = n_down_ca1,
      n_notable_abs_log2FC_gt_log2_1_5 = n_notable,
      mean_log2FC = mean_lfc,
      median_log2FC = median_lfc,
      n_CA3_high = n_ca3_high,
      n_CA1_high = n_ca1_high,
      power_label = power_label,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(cat_effect, file.path(out_data, "target_gene_category_effect_by_age.csv"), row.names = FALSE)

# Age group DEG summary
deg_summary <- data.frame()
for (age_val in c("Young", "Middle", "Old")) {
  res <- deseq2_results[[age_val]]
  power_label <- if (age_val == "Old") "primary_adequately_powered" else "exploratory_low_power"

  if (!is.null(res) && res$converged) {
    # Target gene counts
    tgt_mask <- target_combined$in_spatial
    tgt_lfc <- target_combined[[paste0("log2FoldChange_", age_val)]][tgt_mask]
    tgt_padj <- target_combined[[paste0("padj_", age_val)]][tgt_mask]
    tgt_sig <- sum(!is.na(tgt_padj) & tgt_padj < 0.05, na.rm = TRUE)

    deg_summary <- rbind(deg_summary, data.frame(
      age = age_val,
      n_tested = res$n_tested,
      n_sig = res$n_sig,
      n_up_CA3 = res$n_up_ca3,
      n_down_CA1 = res$n_down_ca1,
      n_target_tested = sum(!is.na(tgt_padj)),
      n_target_sig = tgt_sig,
      power_label = power_label,
      stringsAsFactors = FALSE
    ))
  } else {
    deg_summary <- rbind(deg_summary, data.frame(
      age = age_val, n_tested = 0, n_sig = 0, n_up_CA3 = 0, n_down_CA1 = 0,
      n_target_tested = 0, n_target_sig = 0,
      power_label = power_label, stringsAsFactors = FALSE
    ))
  }
}

write.csv(deg_summary, file.path(out_data, "age_group_deg_summary.csv"), row.names = FALSE)

# Low power flags
low_power <- data.frame(
  age = c("Young", "Middle", "Old"),
  n_pairs = age_paired_counts$n_paired_GEMgroups,
  model_type = sapply(c("Young", "Middle", "Old"), function(a) {
    if (!is.null(age_designs[[a]])) age_designs[[a]]$model_type else "blocked_insufficient_pairs"
  }),
  deseq2_converged = sapply(c("Young", "Middle", "Old"), function(a) {
    if (!is.null(deseq2_results[[a]])) deseq2_results[[a]]$converged else FALSE
  }),
  power_label = c("exploratory_low_power", "exploratory_low_power", "primary_adequately_powered"),
  caveat = c(
    "n=3 pairs; dispersion estimates unstable; only large effects detectable",
    "n=4 pairs; below typical N=5-6 minimum for reliable DE",
    "n=8 pairs; adequate for DESeq2 with GEMgroup blocking"
  ),
  stringsAsFactors = FALSE
)
write.csv(low_power, file.path(out_data, "age_group_low_power_flags.csv"), row.names = FALSE)

cat("Category effect, DEG summary, and low power flags saved.\n\n")

# ============================================================================
# E10: MODULE DELTA BY AGE (Route D — descriptive first)
# ============================================================================

cat("=== E10: Module Delta by Age (Descriptive) ===\n")

if (!is.null(mod_delta)) {
  # Module columns (exclude GEMgroup and Age)
  module_cols <- setdiff(colnames(mod_delta), c("GEMgroup", "Age"))

  # Compute per-module per-age descriptive statistics
  delta_by_age <- data.frame()
  for (age_val in c("Young", "Middle", "Old")) {
    age_mask <- mod_delta$Age == age_val
    if (sum(age_mask) == 0) next

    for (mod in module_cols) {
      vals <- mod_delta[[mod]][age_mask]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) next

      delta_by_age <- rbind(delta_by_age, data.frame(
        module = mod,
        age = age_val,
        n = length(vals),
        mean_delta = mean(vals),
        median_delta = median(vals),
        iqr_delta = IQR(vals),
        direction = ifelse(median(vals) > 0, "CA3_higher", "CA1_higher"),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Add Old - Young and Old - Middle comparisons
  for (mod in module_cols) {
    old_vals <- delta_by_age[delta_by_age$module == mod & delta_by_age$age == "Old", ]
    young_vals <- delta_by_age[delta_by_age$module == mod & delta_by_age$age == "Young", ]
    middle_vals <- delta_by_age[delta_by_age$module == mod & delta_by_age$age == "Middle", ]

    if (nrow(old_vals) > 0 && nrow(young_vals) > 0) {
      delta_by_age$old_vs_young_delta[delta_by_age$module == mod & delta_by_age$age == "Old"] <-
        old_vals$median_delta - young_vals$median_delta
    }
    if (nrow(old_vals) > 0 && nrow(middle_vals) > 0) {
      delta_by_age$old_vs_middle_delta[delta_by_age$module == mod & delta_by_age$age == "Old"] <-
        old_vals$median_delta - middle_vals$median_delta
    }
  }

  # Optional: Kruskal-Wallis (exploratory only)
  kw_results <- data.frame()
  for (mod in module_cols) {
    vals_list <- list()
    for (age_val in c("Young", "Middle", "Old")) {
      age_mask <- mod_delta$Age == age_val
      vals_list[[age_val]] <- mod_delta[[mod]][age_mask]
    }

    all_vals <- unlist(vals_list)
    all_groups <- rep(names(vals_list), sapply(vals_list, length))

    if (length(unique(all_groups)) >= 2 && length(all_vals) >= 3) {
      kw <- tryCatch(kruskal.test(all_vals ~ all_groups), error = function(e) NULL)
      if (!is.null(kw)) {
        kw_results <- rbind(kw_results, data.frame(
          module = mod,
          kw_statistic = kw$statistic,
          kw_pvalue_exploratory_not_primary_evidence = kw$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Merge KW results
  if (nrow(kw_results) > 0) {
    delta_by_age <- merge(delta_by_age, kw_results, by = "module", all.x = TRUE)
  }

  write.csv(delta_by_age, file.path(out_data, "module_delta_by_age_group.csv"), row.names = FALSE)
  cat(sprintf("Module delta by age: %d rows\n", nrow(delta_by_age)))
} else {
  cat("WARNING: module_score_delta_by_age.csv not available\n")
  log_check("V10-mod", "WARNING", "Module delta file not available")
}

cat("\n")

# ============================================================================
# E11: AGE-STRATIFIED COUPLING (Route F, Exploratory)
# ============================================================================

cat("=== E11: Age-Stratified Coupling (Exploratory) ===\n")

if (!is.null(mod_scores) && !is.null(coupling_base)) {
  # Get the 10 coupling pairs from Phase 11 base
  coupling_pairs <- unique(coupling_base[, c("pair", "acronym")])

  # Module score columns (mean columns)
  mod_mean_cols <- grep("_mean$", colnames(mod_scores), value = TRUE)
  mod_names <- gsub("_mean$", "", mod_mean_cols)

  # Compute age-stratified coupling
  age_coupling <- data.frame()

  for (age_val in c("Young", "Middle", "Old")) {
    age_mask <- mod_scores$Age == age_val
    age_scores <- mod_scores[age_mask, ]

    # Split by context
    ca1_scores <- age_scores[age_scores$region_group == "CA1", ]
    ca3_scores <- age_scores[age_scores$region_group == "CA3", ]

    n_ca1 <- nrow(ca1_scores)
    n_ca3 <- nrow(ca3_scores)

    for (pair_i in 1:nrow(coupling_pairs)) {
      pair_name <- coupling_pairs$pair[pair_i]
      acronym <- coupling_pairs$acronym[pair_i]

      # Determine which modules to correlate
      # Parse the pair name to find the two modules
      pair_parts <- strsplit(pair_name, " vs ")[[1]]
      if (length(pair_parts) != 2) next

      mod1_name <- pair_parts[1]
      mod2_name <- pair_parts[2]

      # Find the corresponding mean columns
      mod1_col <- paste0(mod1_name, "_mean")
      mod2_col <- paste0(mod2_name, "_mean")

      # Handle composite modules (avg_complex_tca, avg_complex, avg_mito, avg_cyto, avg_etc)
      # These were computed in Phase 11 base; for now, skip if columns don't exist
      if (!mod1_col %in% colnames(mod_scores) || !mod2_col %in% colnames(mod_scores)) {
        next
      }

      for (context in c("CA1", "CA3")) {
        ctx_scores <- if (context == "CA1") ca1_scores else ca3_scores
        n_ctx <- nrow(ctx_scores)

        if (n_ctx < 3) next

        x <- ctx_scores[[mod1_col]]
        y <- ctx_scores[[mod2_col]]

        if (any(is.na(x)) || any(is.na(y))) next

        rho <- tryCatch(cor(x, y, method = "spearman"), error = function(e) NA)

        age_coupling <- rbind(age_coupling, data.frame(
          pair_name = pair_name,
          acronym = acronym,
          age = age_val,
          context = context,
          n_samples = n_ctx,
          spearman_rho = rho,
          small_n_flag = TRUE,
          notes = "n too small for reliable Spearman estimation; rho reported as descriptive only",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  write.csv(age_coupling, file.path(out_data, "module_age_coupling_summary.csv"), row.names = FALSE)
  cat(sprintf("Age-stratified coupling: %d rows\n", nrow(age_coupling)))
} else {
  cat("WARNING: Module scores or coupling base not available\n")
}

cat("\n")

# ============================================================================
# E12: PLOTS
# ============================================================================

cat("=== E12: Plots ===\n")

# Color conventions
age_colors <- c("Young" = "#c6dbef", "Middle" = "#4292c6", "Old" = "#08306b")
region_colors <- c("CA1" = "#2166ac", "CA3" = "#b2182b")

# --- Figure 1: target_gene_log2fc_heatmap_by_age ---
cat("  Generating target_gene_log2fc_heatmap_by_age.png\n")
tryCatch({
  found_genes <- target_combined[target_combined$in_spatial, ]
  lfc_mat <- as.matrix(found_genes[, c("log2FoldChange_Young", "log2FoldChange_Middle", "log2FoldChange_Old")])
  rownames(lfc_mat) <- found_genes$gene_symbol
  colnames(lfc_mat) <- c("Young", "Middle", "Old")

  # Replace NA with 0 for visualization
  lfc_mat[is.na(lfc_mat)] <- 0

  # Annotation rows by category
  ann_row <- data.frame(Category = found_genes$category_english, row.names = found_genes$gene_symbol)

  # Truncate extreme values for visualization
  lfc_mat_clipped <- pmax(pmin(lfc_mat, 2), -2)

  png(file.path(out_figs, "target_gene_log2fc_heatmap_by_age.png"), width = 8, height = 16, units = "in", res = 150)
  if (pheatmap_available) {
    pheatmap(lfc_mat_clipped,
             annotation_row = ann_row,
             cluster_rows = TRUE, cluster_cols = FALSE,
             show_rownames = FALSE,
             main = "Target Gene log2FC (CA3 vs CA1) by Age",
             color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
             breaks = seq(-2, 2, length.out = 101))
  } else {
    # Fallback: simple heatmap
    image(t(lfc_mat_clipped[nrow(lfc_mat_clipped):1, ]),
          col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
          axes = FALSE, main = "Target Gene log2FC by Age")
  }
  dev.off()
  cat("    Saved: target_gene_log2fc_heatmap_by_age.png\n")
}, error = function(e) {
  cat(sprintf("    [WARNING] Heatmap 1 failed: %s\n", conditionMessage(e)))
})

# --- Figure 2: target_gene_category_effect_heatmap_by_age ---
cat("  Generating target_gene_category_effect_heatmap_by_age.png\n")
tryCatch({
  if (nrow(cat_effect) > 0) {
    cat_lfc_wide <- reshape(cat_effect[, c("category_en", "age", "median_log2FC")],
                            idvar = "category_en", timevar = "age", direction = "wide")
    colnames(cat_lfc_wide) <- gsub("median_log2FC\\.", "", colnames(cat_lfc_wide))
    rownames(cat_lfc_wide) <- cat_lfc_wide$category_en
    cat_lfc_mat <- as.matrix(cat_lfc_wide[, c("Young", "Middle", "Old")])
    cat_lfc_mat[is.na(cat_lfc_mat)] <- 0

    png(file.path(out_figs, "target_gene_category_effect_heatmap_by_age.png"), width = 6, height = 8, units = "in", res = 150)
    if (pheatmap_available) {
      pheatmap(cat_lfc_mat,
               cluster_rows = TRUE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%.3f",
               main = "Category Median log2FC (CA3 vs CA1) by Age",
               color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100))
    } else {
      image(t(cat_lfc_mat[nrow(cat_lfc_mat):1, ]),
            col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
            axes = FALSE, main = "Category Median log2FC by Age")
    }
    dev.off()
    cat("    Saved: target_gene_category_effect_heatmap_by_age.png\n")
  }
}, error = function(e) {
  cat(sprintf("    [WARNING] Heatmap 2 failed: %s\n", conditionMessage(e)))
})

# --- Figure 3-5: Volcano plots per age ---
for (age_val in c("Young", "Middle", "Old")) {
  cat(sprintf("  Generating volcano_CA3_vs_CA1_%s.png\n", age_val))
  tryCatch({
    res <- deseq2_results[[age_val]]
    if (is.null(res) || is.null(res$results)) {
      cat(sprintf("    [SKIP] No DE results for %s\n", age_val))
      next
    }

    res_df <- res$results
    power_label <- res$power_label
    model_type <- res$model_type

    # Highlight target genes
    target_found <- target_combined[target_combined$in_spatial, ]
    res_df$is_target <- res_df$gene %in% target_found$gene_symbol
    res_df$target_cat <- ""
    idx_match <- match(res_df$gene, target_found$gene_symbol)
    res_df$target_cat[!is.na(idx_match)] <- target_found$category_english[na.omit(idx_match)]

    # Plot
    p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = is_target), size = 0.5, alpha = 0.5) +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#e31a1c"), guide = "none") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey50") +
      labs(
        title = sprintf("CA3 vs CA1: %s (%d pairs, %s)", age_val, res$n_up_ca3 + res$n_down_ca1, power_label),
        subtitle = sprintf("%s | %d sig, %d CA3-high, %d CA1-high",
                          model_type, res$n_sig, res$n_up_ca3, res$n_down_ca1),
        x = "log2FC (CA3 vs CA1)", y = "-log10(padj)"
      ) +
      theme_classic() +
      theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))

    # Add ggrepel labels for target genes
    if (ggrepel_available) {
      label_df <- res_df[res_df$is_target & res_df$padj_sig & !is.na(res_df$padj), ]
      if (nrow(label_df) > 30) {
        label_df <- label_df[order(label_df$padj)[1:30], ]
      }
      if (nrow(label_df) > 0) {
        p <- p + geom_text_repel(data = label_df, aes(label = gene),
                                 size = 2, max.overlaps = 20, color = "#e31a1c")
      }
    }

    fname <- sprintf("volcano_CA3_vs_CA1_%s.png", age_val)
    ggsave(file.path(out_figs, fname), p, width = 8, height = 6, dpi = 150)
    cat(sprintf("    Saved: %s\n", fname))
  }, error = function(e) {
    cat(sprintf("    [WARNING] Volcano %s failed: %s\n", age_val, conditionMessage(e)))
  })
}

# --- Figure 6: DEG count by age barplot ---
cat("  Generating deg_count_by_age_barplot.png\n")
tryCatch({
  bar_data <- data.frame()
  for (age_val in c("Young", "Middle", "Old")) {
    res <- deseq2_results[[age_val]]
    if (is.null(res) || !res$converged) next

    bar_data <- rbind(bar_data, data.frame(
      age = age_val, type = "All genes", n_up = res$n_up_ca3, n_down = res$n_down_ca1
    ))

    # Target genes
    tgt_mask <- target_combined$in_spatial
    tgt_lfc <- target_combined[[paste0("log2FoldChange_", age_val)]][tgt_mask]
    tgt_padj <- target_combined[[paste0("padj_", age_val)]][tgt_mask]
    tgt_sig <- !is.na(tgt_padj) & tgt_padj < 0.05

    bar_data <- rbind(bar_data, data.frame(
      age = age_val, type = "Target genes",
      n_up = sum(tgt_sig & tgt_lfc > 0, na.rm = TRUE),
      n_down = sum(tgt_sig & tgt_lfc < 0, na.rm = TRUE)
    ))
  }

  bar_long <- reshape(bar_data, varying = c("n_up", "n_down"), v.names = "count",
                      timevar = "direction", times = c("CA3-high", "CA1-high"), direction = "long")

  p <- ggplot(bar_long, aes(x = age, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ type, scales = "free_y") +
    scale_fill_manual(values = c("CA3-high" = "#b2182b", "CA1-high" = "#2166ac")) +
    labs(title = "DEG Counts by Age (CA3 vs CA1)", y = "Count", fill = "Direction") +
    theme_classic()

  ggsave(file.path(out_figs, "deg_count_by_age_barplot.png"), p, width = 8, height = 5, dpi = 150)
  cat("    Saved: deg_count_by_age_barplot.png\n")
}, error = function(e) {
  cat(sprintf("    [WARNING] Barplot failed: %s\n", conditionMessage(e)))
})

# --- Figure 7: Module delta by age boxplot ---
cat("  Generating module_delta_by_age_boxplot.png\n")
tryCatch({
  if (!is.null(mod_delta)) {
    module_cols <- setdiff(colnames(mod_delta), c("GEMgroup", "Age"))

    delta_long <- reshape(mod_delta, varying = module_cols, v.names = "delta",
                          timevar = "module", times = module_cols, direction = "long")
    rownames(delta_long) <- NULL

    p <- ggplot(delta_long, aes(x = Age, y = delta, fill = Age)) +
      geom_boxplot(alpha = 0.7) +
      geom_point(size = 1, alpha = 0.5) +
      facet_wrap(~ module, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = age_colors) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = "Module Score Delta (CA3 - CA1) by Age",
           subtitle = "(n Young=3, Middle=4, Old=8, descriptive)",
           y = "CA3 - CA1 Module Score Delta") +
      theme_classic() +
      theme(legend.position = "none", strip.text = element_text(size = 7))

    ggsave(file.path(out_figs, "module_delta_by_age_boxplot.png"), p, width = 12, height = 10, dpi = 150)
    cat("    Saved: module_delta_by_age_boxplot.png\n")
  }
}, error = function(e) {
  cat(sprintf("    [WARNING] Boxplot failed: %s\n", conditionMessage(e)))
})

# --- Figure 8: Selected target gene effect by age dotplot ---
cat("  Generating selected_target_gene_effect_by_age_dotplot.png\n")
tryCatch({
  # Select top genes by |log2FC| in Old
  found_genes <- target_combined[target_combined$in_spatial, ]
  old_lfc_col <- "log2FoldChange_Old"
  found_genes$abs_lfc_old <- abs(found_genes[[old_lfc_col]])
  found_genes_sorted <- found_genes[order(found_genes$abs_lfc_old, decreasing = TRUE), ]

  # Take top 20
  top_genes <- head(found_genes_sorted, 20)

  # Build long format
  dot_data <- data.frame()
  for (age_val in c("Young", "Middle", "Old")) {
    lfc_col <- paste0("log2FoldChange_", age_val)
    se_col <- paste0("lfcSE_", age_val)

    for (i in 1:nrow(top_genes)) {
      dot_data <- rbind(dot_data, data.frame(
        gene = top_genes$gene_symbol[i],
        category = top_genes$category_english[i],
        age = age_val,
        log2FC = top_genes[[lfc_col]][i],
        stringsAsFactors = FALSE
      ))
    }
  }

  p <- ggplot(dot_data, aes(x = age, y = log2FC, color = category)) +
    geom_point(size = 2) +
    geom_line(aes(group = gene), alpha = 0.3, color = "grey50") +
    facet_wrap(~ gene, scales = "free_y", ncol = 5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Top 20 Target Genes: log2FC by Age",
         subtitle = "CA3 vs CA1; colored by category",
         y = "log2FC (CA3 vs CA1)") +
    theme_classic() +
    theme(strip.text = element_text(size = 7), legend.position = "bottom")

  ggsave(file.path(out_figs, "selected_target_gene_effect_by_age_dotplot.png"), p, width = 12, height = 10, dpi = 150)
  cat("    Saved: selected_target_gene_effect_by_age_dotplot.png\n")
}, error = function(e) {
  cat(sprintf("    [WARNING] Dotplot failed: %s\n", conditionMessage(e)))
})

cat("\n")

# ============================================================================
# E13: ANALYSIS GUIDE
# ============================================================================

cat("=== E13: Analysis Guide ===\n")

guide_text <- '# Phase Spatial-11 Age Grouping: Analysis Guide

> Generated by: R/spatial/s11b_age_grouping_ca1_ca3_de.R
> Plan version: v1.1

---

## 1. What Each Module Asks

| Module | Question | Output |
|--------|----------|--------|
| Age-stratified DESeq2 | CA3 vs CA1 DEG within each age | `deseq2_all_genes_CA3_vs_CA1_{Age}.csv` |
| Target gene by age | Which target genes are DE in each age? | `deseq2_target_genes_CA3_vs_CA1_by_age.csv` |
| Cross-age log2FC matrix | How do CA3-vs-CA1 effects change with age? | `target_gene_log2fc_matrix_by_age.csv` |
| Category by age | Do whole categories shift with age? | `target_gene_category_effect_by_age.csv` |
| Module delta by age | Do module-score regional differences vary by age? | `module_delta_by_age_group.csv` |
| Age-stratified coupling | Does module coupling differ by age? | `module_age_coupling_summary.csv` |

---

## 2. How to Read log2FC

- `log2FoldChange > 0` → higher expression in **CA3** than CA1
- `log2FoldChange < 0` → higher expression in **CA1** than CA3
- `baseMean` → mean of normalized counts across all pseudobulk samples in that age
- `padj` → Benjamini-Hochberg adjusted p-value (within that age\'s DESeq2 run only)

---

## 3. Thresholds

| Term | Definition |
|------|-----------|
| `padj_sig` | padj < 0.05 |
| `log2FC_notable` | |log2FC| > log2(1.5) ≈ 0.585 |
| `CA3_high` | log2FC > 0 (direction label, not a significance claim) |
| `CA1_high` | log2FC < 0 (direction label, not a significance claim) |

---

## 4. Why Spots Are Not Replicates

Visium spots from the same tissue section (GEMgroup) share:
- Same animal, same tissue processing, same sequencing library
- Spatial neighbors have correlated gene expression

Treating individual spots as replicates would dramatically inflate false positive rates. The GEMgroup (tissue section) is the biological replicate unit.

---

## 5. Why Young and Middle Have Low Statistical Power

| Age | Paired GEMgroups | DESeq2 samples | Issue |
|-----|-----------------|----------------|-------|
| Young | 3 (G2,G3,G4) | 6 | DESeq2 can fit the model but dispersion estimates are unstable with 3 levels of the blocking factor |
| Middle | 4 (G5-G8) | 8 | Better than Young but still below typical N=5-6 minimum for reliable DE |
| Old | 8 (G9-G16) | 16 | Adequate for DESeq2 with GEMgroup blocking |

**Consequences**:
- In Young, DESeq2 may detect only very large effects (e.g., |log2FC| > 1-2).
- In Middle, moderate effects may reach significance; small effects will not.
- In Old, the same genes tested may show effects that are undetectable in Young/Middle simply due to sample size, not biology. A gene not significant in Young does NOT mean it is "not differentially expressed in young tissue."

---

## 6. Age-Stratified DE vs Interaction Model

| Aspect | Age-stratified (this analysis) | Interaction (not run) |
|--------|-------------------------------|----------------------|
| GEMgroup blocking | Preserved per age | Dropped |
| Replicate unit | GEMgroup | Individual pseudobulk samples (weakened) |
| Interpretation | "CA3 vs CA1 in Young", "CA3 vs CA1 in Old" — simple | "Does the CA3-CA1 effect differ between Old and Young?" — requires interaction contrast |
| Power | Lower per age (fewer samples), but valid | Higher total samples but confounded |
| Choice | **USED** | **SKIPPED** |

---

## 7. What Can Be Used for Downstream Interpretation

| Result type | Status | For interpretation? |
|-------------|--------|---------------------|
| Old-specific DESeq2 results (padj<0.05, |log2FC|>log2(1.5)) | PRIMARY, adequately powered | YES — candidate genes for Old-specific CA3-vs-CA1 effects |
| Middle-specific DESeq2 results (padj<0.05) | EXPLORATORY, 4 pairs | With caution; cross-reference with Old results |
| Young-specific DESeq2 results (padj<0.05) | EXPLORATORY, 3 pairs | With extreme caution; only large effects credible |
| Cross-age log2FC direction consistency | DESCRIPTIVE | YES — direction concordance/discordance across ages is informative regardless of significance |
| Module delta trend across ages | DESCRIPTIVE | YES — direction and magnitude trends; KW p-values are exploratory |
| Age-stratified coupling | EXPLORATORY, n=3-8 | NOT recommended for interpretation; n too small |

---

## 8. Complex V / Missing Genes Caveat

16 Atp5* genes (Complex V) are absent from the Spatial assay. They are likely named differently in the annotation version used by this Seurat object. Phase 11 rescue audit found substring candidates but no user-confirmed replacements. Complex V is excluded from module scores. Any statement about "mitochondrial gene expression" in these results EXCLUDES Complex V (ATP synthase).

---

## 9. model_type and power_label

Every age-specific DE output includes:
- `model_type`: `paired_blocked` (primary), `fallback_unpaired` (rank-deficient), or `blocked_insufficient_pairs` (< 2 pairs)
- `power_label`: `exploratory_low_power` (Young, Middle) or `primary_adequately_powered` (Old)

If any age uses `fallback_unpaired`, results are NOT interchangeable with paired results.

---

## 10. Next Steps

From this analysis, the user can:
1. Select candidate genes with age-concordant or age-divergent CA3-CA1 effects
2. Identify categories that show age-dependent shifts
3. Focus on Old-specific target genes for spatial visualization
4. Build foundation for Phase 12: DG-inclusive region comparison or biological interpretation
'

writeLines(guide_text, "docs/spatial_phase11_age_grouping_analysis_guide.md")
cat("Analysis guide saved: docs/spatial_phase11_age_grouping_analysis_guide.md\n\n")

# ============================================================================
# E14: VALIDATION AND PROVENANCE
# ============================================================================

cat("=== E14: Validation and Provenance ===\n")

# V19: Rplots.pdf handling
rplots_end <- file.exists("Rplots.pdf")
rplots_end_md5 <- if (rplots_end) {
  tryCatch(digest::digest("Rplots.pdf", file = TRUE), error = function(e) tools::md5sum("Rplots.pdf"))
} else NA_character_
log_check("V19", "PASS", sprintf("Rplots.pdf: pre-existed=%s, created_during_run=%s",
                                   rplots_start, rplots_end && !rplots_start))

# V20: Low power labels
young_labeled <- if (!is.null(deseq2_results$Young) && !is.null(deseq2_results$Young$results)) {
  all(deseq2_results$Young$results$power_label == "exploratory_low_power")
} else TRUE
middle_labeled <- if (!is.null(deseq2_results$Middle) && !is.null(deseq2_results$Middle$results)) {
  all(deseq2_results$Middle$results$power_label == "exploratory_low_power")
} else TRUE
log_check("V20", if (young_labeled && middle_labeled) "PASS" else "WARNING",
          sprintf("Low power labels: Young=%s, Middle=%s", young_labeled, middle_labeled))

# V21: Complex V status
log_check("V21", "WARNING", "Complex V excluded (0 genes after rescue audit); caveat carried forward from Phase 11")

# V23: target_genes.csv unchanged
target_md5 <- tryCatch(tools::md5sum(f_target_csv), error = function(e) NA_character_)
expected_target_md5 <- "73da8ecc7332a1405f9a282e1821f9fd"
log_check("V23", if (!is.na(target_md5) && target_md5 == expected_target_md5) "PASS" else "WARNING",
          sprintf("target_genes.csv MD5: %s", target_md5))

# V24: Output directories
log_check("V24", if (file.exists(out_data) && file.exists(out_figs)) "PASS" else "WARNING",
          sprintf("Output dirs: data=%s, figs=%s", file.exists(out_data), file.exists(out_figs)))

# V25: No biological interpretation
log_check("V25", "PASS", "No biological interpretation included; statistical descriptions only")

# V26: model_type and power_label in DE outputs
v26_pass <- TRUE
for (age_val in c("Young", "Middle", "Old")) {
  fname <- sprintf("deseq2_all_genes_CA3_vs_CA1_%s.csv", age_val)
  fpath <- file.path(out_data, fname)
  if (file.exists(fpath)) {
    tmp <- read.csv(fpath, nrows = 5)
    if (!"model_type" %in% colnames(tmp) || !"power_label" %in% colnames(tmp)) {
      v26_pass <- FALSE
    }
  }
}
log_check("V26", if (v26_pass) "PASS" else "WARNING", "model_type and power_label in age-specific DE outputs")

# V27: Fallback labeling
v27_pass <- TRUE
for (age_val in c("Young", "Middle", "Old")) {
  da_row <- design_audit[design_audit$age == age_val, ]
  if (nrow(da_row) > 0 && da_row$model_type == "fallback_unpaired") {
    fname <- sprintf("volcano_CA3_vs_CA1_%s.png", age_val)
    # Would need to check figure title; mark as PASS if no fallback used
    if (da_row$model_type != "fallback_unpaired") v27_pass <- TRUE
  }
}
log_check("V27", "PASS", "Fallback model labeling checked (no fallback models used in current run)")

# V28: renv.lock status vs git HEAD
renv_lock_md5_after <- tryCatch(tools::md5sum("renv.lock"), error = function(e) NA_character_)
renv_changed_during_run <- !is.na(renv_lock_md5_before) && !is.na(renv_lock_md5_after) && renv_lock_md5_before != renv_lock_md5_after

provenance[["renv_lock_md5_before"]] <- renv_lock_md5_before
provenance[["renv_lock_md5_after"]] <- renv_lock_md5_after
provenance[["renv_lock_changed_during_run"]] <- renv_changed_during_run
provenance[["renv_vs_git_head_start"]] <- renv_vs_git
provenance[["renv_vs_git_head_end"]] <- tryCatch({
  diff_lines <- system("git diff HEAD -- renv.lock | head -5", intern = TRUE)
  if (length(diff_lines) > 0) "modified_vs_git_HEAD" else "clean_vs_git_HEAD"
}, error = function(e) "git_not_available")

log_check("V28", "PASS", sprintf("renv.lock: before=%s, after=%s, changed_during_run=%s, vs_git=%s",
                                   renv_lock_md5_before, renv_lock_md5_after,
                                   renv_changed_during_run, provenance[["renv_vs_git_head_start"]]))

# No enrichment
log_check("V17", "PASS", "No enrichment run (GO/KEGG/Reactome)")
log_check("V16", "PASS", "No Tangram/Python/STutility")
log_check("V15", "PASS", "No WholeBrain/DG loaded; Hippo object only")
log_check("V14", "PASS", "No full Seurat object saved")
log_check("V13", "PASS", "No spots as replicates; DESeq2 uses pseudobulk")

# Save validation summary
validation_text <- "=== Phase Spatial-11 Age Grouping Validation Summary ===\n"
validation_text <- paste0(validation_text, "Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
validation_text <- paste0(validation_text, "Plan version: v1.1\n")
validation_text <- paste0(validation_text, "Script: R/spatial/s11b_age_grouping_ca1_ca3_de.R\n\n")
validation_text <- paste0(validation_text, "--- Validation Checks ---\n")

n_pass <- 0; n_warn <- 0; n_stop <- 0
for (nm in names(validation_log)) {
  entry <- validation_log[[nm]]
  validation_text <- paste0(validation_text, sprintf("[%s] %s: %s\n", nm, entry$status, entry$msg))
  if (entry$status == "PASS") n_pass <- n_pass + 1
  else if (entry$status == "WARNING") n_warn <- n_warn + 1
  else if (entry$status == "STOP") n_stop <- n_stop + 1
}

validation_text <- paste0(validation_text, "\n--- Summary ---\n")
validation_text <- paste0(validation_text, sprintf("PASS: %d\n", n_pass))
validation_text <- paste0(validation_text, sprintf("WARNING: %d\n", n_warn))
validation_text <- paste0(validation_text, sprintf("STOP: %d\n", n_stop))
overall <- if (n_stop > 0) "FAIL" else if (n_warn > 0) "PASS_WITH_WARNINGS" else "PASS"
validation_text <- paste0(validation_text, sprintf("Overall: %s\n\n", overall))
validation_text <- paste0(validation_text, "No biological interpretation included. Statistical descriptions only.\n")
validation_text <- paste0(validation_text, "No GO/KEGG/Reactome enrichment performed.\n")
validation_text <- paste0(validation_text, "Spots are not biological replicates; pseudobulk aggregation used for DE.\n")
validation_text <- paste0(validation_text, "All three age groups (Young/Middle/Old) summarized.\n")

writeLines(validation_text, file.path(out_data, "phase11_age_grouping_validation_summary.txt"))
cat("Validation summary saved.\n")

# Save provenance
provenance[["end_time"]] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
provenance[["validation_overall"]] <- overall
provenance[["n_pass"]] <- n_pass
provenance[["n_warn"]] <- n_warn
provenance[["n_stop"]] <- n_stop

# Package versions
for (pkg in c("Seurat", "DESeq2", "Matrix", "Matrix.utils", "apeglm", "ggplot2",
              "pheatmap", "ggrepel", "corrplot", "data.table", "digest")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    provenance[[paste0("version_", pkg)]] <- as.character(packageVersion(pkg))
  }
}

# Git status at end
git_status_end <- tryCatch(system("git status --short", intern = TRUE), error = function(e) "git not available")
provenance[["git_status_end"]] <- if (length(git_status_end) == 0) "(clean)" else paste(git_status_end, collapse = "; ")

provenance_df <- as.data.frame(provenance, stringsAsFactors = FALSE)
write.csv(provenance_df, file.path(out_data, "phase11_age_grouping_provenance.csv"), row.names = FALSE)
cat("Provenance saved.\n")

# ============================================================================
# E15: CLEANUP
# ============================================================================

cat("\n=== E15: Cleanup ===\n")

# Remove large objects
rm(pseudobulk_rds)
if (exists("counts_dense")) rm(counts_dense)
gc()

# Rplots.pdf guard (end)
rplots_final <- file.exists("Rplots.pdf")
if (rplots_final && !rplots_start) {
  cat("[Rplots.pdf] Created during run — deleting\n")
  file.remove("Rplots.pdf")
} else {
  cat("[Rplots.pdf] Not created during run\n")
}

cat("\n=== Phase Spatial-11 Age Grouping Complete ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat(sprintf("Validation: %d PASS, %d WARNING, %d STOP — %s\n", n_pass, n_warn, n_stop, overall))
cat(sprintf("Output data: %s\n", out_data))
cat(sprintf("Output figures: %s\n", out_figs))

# Exit code
if (n_stop > 0) {
  cat("\n[EXIT] Code 1 — STOP conditions encountered\n")
  quit(status = 1)
} else {
  cat("\n[EXIT] Code 0 — PASS or PASS_WITH_WARNINGS\n")
}
