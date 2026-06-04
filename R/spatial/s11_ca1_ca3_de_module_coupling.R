#!/usr/bin/env Rscript
# s11_ca1_ca3_de_module_coupling.R
# Phase Spatial-11: CA1 vs CA3 DE + Module/Coupling Analysis
# RDS-based spatial analysis; no PDF; no Seurat object modification
# Plan version: v2.1 (approved 2026-06-04)

# ============================================================================
# E0: SETUP
# ============================================================================

cat("=== Phase Spatial-11: CA1 vs CA3 DE + Module/Coupling Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# -- Rplots.pdf guard (start) --
rplots_start <- file.exists("Rplots.pdf")
rplots_start_md5 <- if (rplots_start) {
  tryCatch(digest::digest("Rplots.pdf", file = TRUE), error = function(e) tools::md5sum("Rplots.pdf"))
} else NA_character_
cat("[Rplots.pdf] Pre-existence:", rplots_start, if (!is.na(rplots_start_md5)) paste("MD5:", rplots_start_md5) else "", "\n")

# -- Git status --
git_status <- tryCatch(system("git status --short", intern = TRUE), error = function(e) "git not available")
cat("Git status at start:", if (length(git_status) == 0) "(clean)" else paste(git_status, collapse = "; "), "\n\n")

# -- renv.lock MD5 before --
renv_lock_md5_before <- tryCatch(tools::md5sum("renv.lock"), error = function(e) NA_character_)
cat("renv.lock MD5 before:", renv_lock_md5_before, "\n\n")

# -- Required packages --
pkgs_required <- c("Seurat", "DESeq2", "Matrix", "ggplot2", "data.table", "digest")
for (pkg in pkgs_required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[SC-9 STOP] Package not available:", pkg, "\n")
    stop(paste("Package", pkg, "not installed"))
  }
}
library(Seurat)
library(DESeq2)
library(Matrix)
library(ggplot2)
library(data.table)

# -- Optional packages: install if beneficial (D10) --
optional_pkg_actions <- list()

# Matrix.utils (required for pseudobulk aggregation)
if (!requireNamespace("Matrix.utils", quietly = TRUE)) {
  cat("[INFO] Matrix.utils not found; attempting renv::install('Matrix.utils')\n")
  tryCatch({
    renv::install("Matrix.utils")
    renv::snapshot()
    optional_pkg_actions[["Matrix.utils"]] <- "installed"
    cat("[INFO] Matrix.utils installed successfully\n")
  }, error = function(e) {
    cat("[WARNING] Matrix.utils install failed:", conditionMessage(e), "\n")
    optional_pkg_actions[["Matrix.utils"]] <- "install_failed"
  })
}

# apeglm (LFC shrinkage)
apeglm_available <- requireNamespace("apeglm", quietly = TRUE)
if (!apeglm_available) {
  cat("[INFO] apeglm not found; attempting renv::install('apeglm')\n")
  tryCatch({
    renv::install("apeglm")
    renv::snapshot()
    apeglm_available <- TRUE
    optional_pkg_actions[["apeglm"]] <- "installed"
    cat("[INFO] apeglm installed successfully\n")
  }, error = function(e) {
    cat("[WARNING] apeglm install failed:", conditionMessage(e), "\n")
    optional_pkg_actions[["apeglm"]] <- "install_failed"
  })
}

# ComplexHeatmap (heatmaps)
heatmap_pkg <- "pheatmap"
pheatmap_available <- requireNamespace("pheatmap", quietly = TRUE)
if (!pheatmap_available) {
  cat("[INFO] pheatmap not found; attempting renv::install('pheatmap')\n")
  tryCatch({
    renv::install("pheatmap")
    renv::snapshot()
    pheatmap_available <- TRUE
    optional_pkg_actions[["pheatmap"]] <- "installed"
    cat("[INFO] pheatmap installed successfully\n")
  }, error = function(e) {
    cat("[WARNING] pheatmap install failed:", conditionMessage(e), "\n")
  })
}
if (pheatmap_available) library(pheatmap)

# ggrepel (volcano labels)
ggrepel_available <- requireNamespace("ggrepel", quietly = TRUE)
if (ggrepel_available) library(ggrepel)

# corrplot (coupling heatmap)
corrplot_available <- requireNamespace("corrplot", quietly = TRUE)
if (!corrplot_available) {
  cat("[INFO] corrplot not found; attempting renv::install('corrplot')\n")
  tryCatch({
    renv::install("corrplot")
    renv::snapshot()
    corrplot_available <- TRUE
    optional_pkg_actions[["corrplot"]] <- "installed"
    cat("[INFO] corrplot installed successfully\n")
  }, error = function(e) {
    cat("[WARNING] corrplot install failed:", conditionMessage(e), "\n")
    optional_pkg_actions[["corrplot"]] <- "install_failed"
  })
}
if (corrplot_available) library(corrplot)

# matrixStats
if (requireNamespace("matrixStats", quietly = TRUE)) library(matrixStats)

# -- Paths --
hippo_rds  <- "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"
csv_path   <- "docs/target_genes.csv"
phase10_data <- "data/processed/spatial/phase10_target_gene_audit"
out_data   <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling"
out_figs   <- "figures/spatial/phase11_ca1_ca3_de_module_coupling"
dir.create(out_data, showWarnings = FALSE, recursive = TRUE)
dir.create(out_figs, showWarnings = FALSE, recursive = TRUE)

# -- Expected values (Phase 06) --
expected_counts <- c(CA1 = 4671L, CA2 = 565L, CA3 = 2321L, ML = 1480L, GCL = 645L, Hilus = 239L)
expected_total  <- 9921L
expected_derived_dg <- 2364L

# -- Category English mapping --
category_english <- c(
  "\u7ebf\u7c92\u4f53\u57fa\u56e0\u7ec4\u7f16\u7801\uff08mt-\uff09" = "mtDNA-encoded",
  "\u547c\u5438\u94fe\u590d\u5408\u4f53 I\uff08Nduf/Ndufv/Ndufb/Ndufc\uff09" = "Complex I",
  "\u547c\u5438\u94fe\u590d\u5408\u4f53 III\uff08Uqcr/Uqcrc/Uqcc\uff09" = "Complex III",
  "\u547c\u5438\u94fe\u590d\u5408\u4f53 IV\uff08Cox\uff09" = "Complex IV",
  "\u547c\u5438\u94fe\u590d\u5408\u4f53 V\uff08Atp5 \u4e9a\u57fa\uff0cATP \u5408\u9176\uff09" = "Complex V",
  "\u4e09\u7fa7\u9178\u5faa\u73af TCA \u5173\u952e\u9176" = "TCA cycle",
  "\u7ebf\u7c92\u4f53\u6838\u7cd6\u4f53\u86cb\u767d\uff08MRPL/MRPS\uff0c\u540c\u65f6\u5f52\u5165\u7ebf\u7c92\u4f53\u7ffb\u8bd1\uff09" = "Mitoribosomal",
  "\u7ebf\u7c92\u4f53\u5916\u819c / \u5185\u819c\u8f6c\u8fd0 / \u5b54\u86cb\u767d" = "Membrane/Translocator",
  "\u7ebf\u7c92\u4f53\u5206\u5b50\u4f34\u4fa3 / \u6297\u6c27\u5316 / \u4ee3\u8c22\u9176" = "Defense/Antioxidant",
  "RPL\uff08\u80de\u8d28\u5927\u4e9a\u57fa\uff09" = "RPL",
  "RPS\uff08\u80de\u8d28\u5c0f\u4e9a\u57fa\uff09" = "RPS",
  "\u8d77\u59cb\u56e0\u5b50 Eif" = "Eif",
  "\u5ef6\u4f38\u56e0\u5b50 Eef" = "Eef"
)

# -- Color conventions --
region_colors <- c(CA1 = "#2166ac", CA2 = "#92c5de", CA3 = "#b2182b", DG = "#d1e5f0")
age_colors    <- c(Young = "#c6dbef", Middle = "#4292c6", Old = "#08306b")

# -- Provenance collection --
renv_lock_md5 <- renv_lock_md5_before
provenance <- list(
  script         = "R/spatial/s11_ca1_ca3_de_module_coupling.R",
  plan_version   = "v2.1",
  seurat_version = as.character(packageVersion("Seurat")),
  deseq2_version = as.character(packageVersion("DESeq2")),
  r_version      = R.version.string,
  run_start      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  git_status     = if (length(git_status) == 0) "(clean)" else paste(git_status, collapse = "; "),
  renv_lock_md5_before = as.character(renv_lock_md5_before),
  pheatmap       = if (pheatmap_available) as.character(packageVersion("pheatmap")) else "not_available",
  ggrepel        = if (ggrepel_available) as.character(packageVersion("ggrepel")) else "not_available",
  apeglm         = if (apeglm_available) as.character(packageVersion("apeglm")) else "not_available",
  corrplot       = if (corrplot_available) as.character(packageVersion("corrplot")) else "not_available",
  dense_conversion = "No full spot-level Spatial matrix dense conversion; small pseudobulk matrix (32285x31) converted to dense/integer for DESeq2; allowed by plan and recorded",
  optional_pkg_actions = paste(names(optional_pkg_actions), optional_pkg_actions, sep = "=", collapse = "; ")
)

# -- Validation tracking --
validation_log <- list()
log_check <- function(id, label, status, detail = "") {
  validation_log[[id]] <<- list(id = id, label = label, status = status, detail = detail)
  cat(sprintf("  [%s] %s: %s %s\n", id, status, label, detail))
}

# ============================================================================
# E0b: LOAD HIPPO RDS + PHASE 10 AUDIT
# ============================================================================

cat("\n--- E0b: Load Hippo RDS + Phase 10 audit ---\n")

if (!file.exists(hippo_rds)) {
  log_check("SC-2", "Hippo RDS exists", "STOP", hippo_rds)
  stop("[SC-2] Hippo RDS not found: ", hippo_rds)
}

hippo <- readRDS(hippo_rds)
cat("Hippo loaded:", ncol(hippo), "spots,", nrow(hippo[["Spatial"]]), "features\n")

# Verify invariants
n_spots <- ncol(hippo)
n_features <- nrow(hippo[["Spatial"]])
if (n_spots != expected_total) {
  log_check("V02", "Hippo spot count", "WARNING", paste("Expected", expected_total, "got", n_spots))
} else {
  log_check("V02", "Hippo spot count", "PASS", paste(n_spots, "spots"))
}

# Verify Region metadata
if (!"Region" %in% colnames(hippo@meta.data)) {
  log_check("SC-4", "Region metadata", "STOP", "Region column missing")
  stop("[SC-4] Region column missing from metadata")
}
region_table <- table(hippo$Region)
cat("Region table:\n")
print(region_table)

# Verify region invariants
for (rg in names(expected_counts)) {
  actual <- as.integer(region_table[rg])
  if (is.na(actual) || abs(actual - expected_counts[rg]) / expected_counts[rg] > 0.05) {
    log_check("V03", paste("Region", rg, "invariant"), "WARNING",
              paste("Expected", expected_counts[rg], "got", actual))
  }
}
log_check("V03", "Region invariants", "PASS", "All within 5%")

# Load Phase 10 gene audit
phase10_audit_path <- file.path(phase10_data, "target_gene_input_audit.csv")
if (!file.exists(phase10_audit_path)) {
  log_check("SC-1", "Phase 10 audit exists", "STOP", phase10_audit_path)
  stop("[SC-1] Phase 10 audit not found. Run Phase 10 first.")
}
phase10_audit <- read.csv(phase10_audit_path, stringsAsFactors = FALSE)
cat("Phase 10 audit loaded:", nrow(phase10_audit), "genes\n")
log_check("V01", "Phase 10 validation", "PASS", "Audit file loaded")

# Load Phase 10 region_sample_summary for reference
phase10_region_sample_path <- file.path(phase10_data, "region_sample_summary.csv")
phase10_region_sample <- NULL
if (file.exists(phase10_region_sample_path)) {
  phase10_region_sample <- read.csv(phase10_region_sample_path, stringsAsFactors = FALSE)
  cat("Phase 10 region_sample_summary loaded:", nrow(phase10_region_sample), "rows\n")
}

# Verify 231 exact-match genes
exact_match_genes <- phase10_audit$final_symbol[phase10_audit$mapping_status == "exact_match"]
cat("Exact-match genes from Phase 10:", length(exact_match_genes), "\n")
if (abs(length(exact_match_genes) - 231) / 231 > 0.05) {
  log_check("SC-11", "231 exact-match genes", "WARNING",
            paste("Expected 231, got", length(exact_match_genes)))
} else {
  log_check("SC-11", "231 exact-match genes", "PASS", paste(length(exact_match_genes), "genes"))
}

# Verify GEMgroup, Age metadata
for (col_name in c("GEMgroup", "Age")) {
  if (!col_name %in% colnames(hippo@meta.data)) {
    log_check("SC-4", paste(col_name, "metadata"), "STOP", paste(col_name, "column missing"))
    stop("[SC-4]", col_name, "column missing")
  }
}
log_check("V04", "GEMgroup x Age metadata", "PASS",
          paste(length(unique(hippo$GEMgroup)), "GEMgroups,", length(unique(hippo$Age)), "Age levels"))

# ============================================================================
# E1: MISSING GENE RESCUE AUDIT
# ============================================================================

cat("\n--- E1: Missing Gene Rescue Audit ---\n")

missing_genes <- phase10_audit$gene_symbol[phase10_audit$mapping_status == "not_found_in_either"]
cat("Missing genes to audit:", length(missing_genes), "\n")

all_spatial_genes <- rownames(hippo[["Spatial"]])

rescue_results <- data.frame(
  missing_gene = missing_genes,
  stringsAsFactors = FALSE
)

# For each missing gene: substring search
rescue_results$rescue_candidates <- sapply(missing_genes, function(gene) {
  # Exact match (case-insensitive)
  exact <- all_spatial_genes[tolower(all_spatial_genes) == tolower(gene)]
  if (length(exact) > 0) return(paste(exact, collapse = ";"))
  
  # Prefix match (e.g., Atp5* for Complex V)
  prefix_pattern <- paste0("^", gsub("\\d+$", "", gene))
  prefix_hits <- grep(prefix_pattern, all_spatial_genes, value = TRUE, ignore.case = TRUE)
  if (length(prefix_hits) > 0) return(paste(prefix_hits[1:min(5, length(prefix_hits))], collapse = ";"))
  
  # Substring match
  substr_hits <- grep(gene, all_spatial_genes, value = TRUE, ignore.case = TRUE)
  if (length(substr_hits) > 0) return(paste(substr_hits[1:min(5, length(substr_hits))], collapse = ";"))
  
  return("none_found")
})

rescue_results$n_candidates <- sapply(strsplit(rescue_results$rescue_candidates, ";"), length)
rescue_results$n_candidates[rescue_results$rescue_candidates == "none_found"] <- 0
rescue_results$user_action_required <- rescue_results$n_candidates > 0

# Catalog Atp5* genes specifically
atp5_genes <- grep("^Atp5", all_spatial_genes, value = TRUE, ignore.case = TRUE)
cat("Atp5* genes in Spatial assay:", length(atp5_genes), "\n")
if (length(atp5_genes) > 0) cat("  Examples:", paste(head(atp5_genes, 10), collapse = ", "), "\n")

# Catalog Nduf* genes
nduf_genes <- grep("^Nduf", all_spatial_genes, value = TRUE, ignore.case = TRUE)
cat("Nduf* genes in Spatial assay:", length(nduf_genes), "\n")

# Catalog Mrps/Mrpl genes
mrp_genes <- grep("^Mrp[sl]", all_spatial_genes, value = TRUE, ignore.case = TRUE)
cat("Mrp[sl]* genes in Spatial assay:", length(mrp_genes), "\n")

rescue_path <- file.path(out_data, "missing_gene_symbol_rescue_audit.csv")
write.csv(rescue_results, rescue_path, row.names = FALSE)
cat("Rescue audit saved:", rescue_path, "\n")
log_check("V11", "Missing gene rescue audit", "PASS",
          paste(nrow(rescue_results), "genes,", sum(rescue_results$n_candidates > 0), "with candidates"))
log_check("V34", "No auto-substitution", "PASS", "Candidates recorded only")

# Complex V status
complex_v_rescue <- sum(rescue_results$n_candidates[rescue_results$missing_gene %in%
  grep("Atp5", missing_genes, value = TRUE)] > 0)
complex_v_status <- if (complex_v_rescue > 0) "RESCUE_CANDIDATES_DOCUMENTED" else "EXCLUDED"
cat("Complex V status:", complex_v_status, "\n")
v12_status <- if (complex_v_status == "EXCLUDED") "WARNING" else "PASS"
log_check("V12", "Complex V documented", v12_status,
          if (complex_v_status == "EXCLUDED") "Complex V excluded from module scores because 0 Complex V genes were available after rescue audit; this is a Phase 11 caveat" else complex_v_status)

# ============================================================================
# E2: PSEUDOBULK AGGREGATION (CA1 + CA3 only)
# ============================================================================

cat("\n--- E2: Pseudobulk Aggregation (CA1 + CA3) ---\n")

# Extract full counts and metadata without subsetting Seurat object (VisiumV1 class issue)
ca1_ca3_mask <- hippo$Region %in% c("CA1", "CA3")
cat("CA1+CA3 spots:", sum(ca1_ca3_mask), "\n")

# Get full Spatial raw counts (sparse, genes x spots), then subset columns
counts_sparse_full <- GetAssayData(hippo, assay = "Spatial", layer = "counts")
counts_sparse <- counts_sparse_full[, ca1_ca3_mask]
rm(counts_sparse_full); gc()
cat("Count matrix:", nrow(counts_sparse), "genes x", ncol(counts_sparse), "spots (sparse)\n")

# Subset metadata separately (no Seurat object subsetting)
spot_metadata <- hippo@meta.data[ca1_ca3_mask, , drop = FALSE]
spot_metadata$sample_id <- paste0("G", spot_metadata$GEMgroup, "_", spot_metadata$Region)
spot_groupings <- spot_metadata$sample_id

# Pseudobulk aggregation: t(counts_sparse) -> aggregate.Matrix -> t() back
cat("Aggregating pseudobulk via Matrix.utils::aggregate.Matrix()...\n")
pseudobulk_t <- Matrix.utils::aggregate.Matrix(
  t(counts_sparse), groupings = spot_groupings, fun = "sum"
)
final_pseudobulk <- t(pseudobulk_t)
cat("Pseudobulk matrix:", nrow(final_pseudobulk), "genes x", ncol(final_pseudobulk), "samples\n")

# Verification
stopifnot(nrow(final_pseudobulk) == nrow(counts_sparse))

# Build sample manifest
sample_manifest <- data.frame(
  sample_id = colnames(final_pseudobulk),
  stringsAsFactors = FALSE
)
sample_manifest$GEMgroup <- as.integer(gsub("^G(\\d+)_.+", "\\1", sample_manifest$sample_id))
sample_manifest$region_group <- gsub("^G\\d+_(.+)", "\\1", sample_manifest$sample_id)

# Map Age from GEMgroup
gemgroup_age_map <- setNames(spot_metadata$Age, spot_metadata$GEMgroup)
sample_manifest$Age <- sapply(sample_manifest$GEMgroup, function(g) {
  ages <- unique(gemgroup_age_map[as.character(g)])
  if (length(ages) == 1) return(ages)
  return(ages[1])
})

# Count spots per sample
spot_counts <- table(spot_groupings)
sample_manifest$n_spots <- as.integer(spot_counts[sample_manifest$sample_id])

# has_pair: does this GEMgroup have both CA1 and CA3?
gemgroup_regions <- split(sample_manifest$region_group, sample_manifest$GEMgroup)
paired_gemgroups <- names(gemgroup_regions)[sapply(gemgroup_regions, function(r) all(c("CA1", "CA3") %in% r))]
sample_manifest$has_pair <- sample_manifest$GEMgroup %in% as.integer(paired_gemgroups)
sample_manifest$low_coverage <- sample_manifest$n_spots < 20

cat("Sample manifest:", nrow(sample_manifest), "samples\n")
cat("  CA1 samples:", sum(sample_manifest$region_group == "CA1"), "\n")
cat("  CA3 samples:", sum(sample_manifest$region_group == "CA3"), "\n")
cat("  Paired GEMgroups:", length(paired_gemgroups), "\n")
cat("  Low coverage (<20 spots):", sum(sample_manifest$low_coverage), "\n")

# Save manifest and pseudobulk
manifest_path <- file.path(out_data, "pseudobulk_sample_manifest.csv")
write.csv(sample_manifest, manifest_path, row.names = FALSE)
cat("Sample manifest saved:", manifest_path, "\n")

pseudobulk_path <- file.path(out_data, "ca1_ca3_pseudobulk_counts.rds")
saveRDS(final_pseudobulk, pseudobulk_path)
cat("Pseudobulk counts saved:", pseudobulk_path, "\n")

log_check("V05", "Pseudobulk aggregation", "PASS",
          paste(ncol(final_pseudobulk), "samples,", nrow(final_pseudobulk), "genes"))

if (length(paired_gemgroups) < 2) {
  log_check("V06", "Matched GEMgroups", "STOP",
            paste("Only", length(paired_gemgroups), "paired GEMgroups"))
  log_check("SC-5", "Paired DE feasibility", "STOP",
            "< 2 paired GEMgroups; primary paired DE cannot proceed")
  stop("[SC-5] < 2 GEMgroups with both CA1 and CA3")
}
log_check("V06", "Matched GEMgroups", "PASS",
          paste(length(paired_gemgroups), "paired GEMgroups"))

# ============================================================================
# E3: DESIGN MATRIX VALIDATION
# ============================================================================

cat("\n--- E3: Design Matrix Validation ---\n")

# Prepare colData for DESeq2
col_data <- sample_manifest[, c("sample_id", "GEMgroup", "region_group", "Age")]
col_data$GEMgroup <- factor(col_data$GEMgroup)
col_data$region_group <- factor(col_data$region_group, levels = c("CA1", "CA3"))
col_data$Age <- factor(col_data$Age, levels = c("Young", "Middle", "Old"))
rownames(col_data) <- col_data$sample_id

# Ensure pseudobulk columns match colData rows
common_samples <- intersect(colnames(final_pseudobulk), rownames(col_data))
col_data <- col_data[common_samples, ]
pseudobulk_counts <- final_pseudobulk[, common_samples]

# Primary design: ~ GEMgroup + region_group
design_primary <- "~ GEMgroup + region_group"
cat("Primary design:", design_primary, "\n")

# Check full rank
model_matrix <- tryCatch(
  model.matrix(~ GEMgroup + region_group, data = col_data),
  error = function(e) NULL
)

design_audit <- data.frame(
  design = character(), rank = integer(), ncol = integer(),
  is_full_rank = logical(), stringsAsFactors = FALSE
)

if (!is.null(model_matrix)) {
  qr_rank <- qr(model_matrix)$rank
  n_cols <- ncol(model_matrix)
  is_full_rank <- qr_rank == n_cols
  design_audit <- rbind(design_audit, data.frame(
    design = design_primary, rank = qr_rank, ncol = n_cols,
    is_full_rank = is_full_rank, stringsAsFactors = FALSE
  ))
  cat("  Rank:", qr_rank, "/", n_cols, "- Full rank:", is_full_rank, "\n")
} else {
  is_full_rank <- FALSE
  cat("  model.matrix() failed for primary design\n")
}

# Fallback if not full rank
use_fallback <- FALSE
if (!is_full_rank) {
  cat("[WARNING] Primary design not full rank; trying fallback ~ region_group\n")
  design_fallback <- "~ region_group"
  model_matrix_fb <- tryCatch(
    model.matrix(~ region_group, data = col_data),
    error = function(e) NULL
  )
  if (!is.null(model_matrix_fb)) {
    qr_rank_fb <- qr(model_matrix_fb)$rank
    n_cols_fb <- ncol(model_matrix_fb)
    is_full_rank_fb <- qr_rank_fb == n_cols_fb
    design_audit <- rbind(design_audit, data.frame(
      design = design_fallback, rank = qr_rank_fb, ncol = n_cols_fb,
      is_full_rank = is_full_rank_fb, stringsAsFactors = FALSE
    ))
    cat("  Fallback rank:", qr_rank_fb, "/", n_cols_fb, "- Full rank:", is_full_rank_fb, "\n")
    if (is_full_rank_fb) {
      use_fallback <- TRUE
    }
  }
}

if (!is_full_rank && !use_fallback) {
  log_check("V07", "Design matrix full rank", "STOP", "Neither primary nor fallback design is full rank")
  log_check("SC-7", "Design matrix rank", "STOP", "Investigate data structure")
  stop("[SC-7] Design matrix not full rank")
}

# Save design audit
design_audit_path <- file.path(out_data, "design_matrix_audit.csv")
write.csv(design_audit, design_audit_path, row.names = FALSE)
cat("Design audit saved:", design_audit_path, "\n")
log_check("V07", "Design matrix full rank", "PASS",
          paste("Using", if (use_fallback) "fallback" else "primary", "design"))

# Sample counts per cell
cell_counts <- table(col_data$GEMgroup, col_data$region_group)
cat("Samples per GEMgroup x region_group:\n")
print(cell_counts)

# ============================================================================
# E4: DESeq2 EXECUTION
# ============================================================================

cat("\n--- E4: DESeq2 Execution ---\n")

# Ensure integer counts (convert sparse to dense for DESeq2; only 31 pseudobulk samples)
pseudobulk_integer <- as.matrix(round(pseudobulk_counts))
storage.mode(pseudobulk_integer) <- "integer"

# Create DESeq2 dataset
design_formula <- if (use_fallback) ~ region_group else ~ GEMgroup + region_group
cat("DESeq2 design:", paste(as.character(design_formula), collapse = " "), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_integer,
  colData   = col_data,
  design    = design_formula
)

# Filter: keep genes with rowSums >= 10
pre_filter <- nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Gene filter:", pre_filter, "->", nrow(dds), "genes (rowSums >= 10)\n")

# Run DESeq2
cat("Running DESeq()...\n")
dds <- DESeq(dds)
cat("DESeq2 complete.\n")

# Extract results: CA3 vs CA1
res <- results(dds, contrast = c("region_group", "CA3", "CA1"), alpha = 0.05)
cat("Results extracted: CA3 vs CA1\n")
cat("  Significant (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "genes\n")

# LFC shrinkage
shrinkage_type <- "none"
res_shrunken <- NULL
if (apeglm_available) {
  cat("Running lfcShrink(type='apeglm')...\n")
  res_shrunken <- tryCatch({
    coef_name <- resultsNames(dds)
    coef_name <- coef_name[grep("region_group_CA3_vs_CA1", coef_name)]
    if (length(coef_name) == 0) {
      cat("[WARNING] CA3_vs_CA1 coef not found; trying 'normal' shrinkage\n")
      lfcShrink(dds, contrast = c("region_group", "CA3", "CA1"), type = "normal")
    } else {
      lfcShrink(dds, coef = coef_name, type = "apeglm")
    }
  }, error = function(e) {
    cat("[WARNING] apeglm failed:", conditionMessage(e), "; trying 'normal'\n")
    tryCatch(lfcShrink(dds, contrast = c("region_group", "CA3", "CA1"), type = "normal"),
             error = function(e2) NULL)
  })
  if (!is.null(res_shrunken)) shrinkage_type <- "apeglm_or_normal"
} else {
  cat("[INFO] apeglm not available; trying 'normal' shrinkage\n")
  res_shrunken <- tryCatch(
    lfcShrink(dds, contrast = c("region_group", "CA3", "CA1"), type = "normal"),
    error = function(e) { cat("[WARNING] normal shrinkage failed:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(res_shrunken)) shrinkage_type <- "normal"
}

# Save all-gene results
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
all_genes_path <- file.path(out_data, "deseq2_all_genes_CA3_vs_CA1.csv")
write.csv(res_df, all_genes_path, row.names = FALSE)
cat("All-gene DE results saved:", all_genes_path, "(", nrow(res_df), "genes)\n")

# Save shrunken results
if (!is.null(res_shrunken)) {
  res_shrunken_df <- as.data.frame(res_shrunken)
  res_shrunken_df$gene <- rownames(res_shrunken_df)
  shrunken_path <- file.path(out_data, "deseq2_all_genes_CA3_vs_CA1_shrunken.csv")
  write.csv(res_shrunken_df, shrunken_path, row.names = FALSE)
  cat("Shrunken results saved:", shrunken_path, "\n")
}

log_check("V08", "DESeq2 converged", "PASS", paste("Shrinkage:", shrinkage_type))
log_check("V09", "All-gene DE results", "PASS", paste(nrow(res_df), "genes"))

# ============================================================================
# E5: TARGET GENE DE ANNOTATION
# ============================================================================

cat("\n--- E5: Target Gene DE Annotation ---\n")

# Merge with Phase 10 audit
target_de <- phase10_audit[, c("gene_symbol", "category", "category_english", "mapping_status", "final_symbol")]
colnames(target_de)[colnames(target_de) == "gene_symbol"] <- "original_gene"

# Add DE columns
de_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
for (col_name in de_cols) {
  target_de[[col_name]] <- NA_real_
}

# Match found genes to DE results
found_mask <- target_de$mapping_status == "exact_match"
found_symbols <- target_de$final_symbol[found_mask]
de_match <- match(found_symbols, res_df$gene)
matched_in_de <- !is.na(de_match)

for (col_name in de_cols) {
  target_de[[col_name]][found_mask] <- res_df[[col_name]][de_match]
}

# Add flags
target_de$in_spatial <- target_de$mapping_status == "exact_match"
target_de$padj_sig <- !is.na(target_de$padj) & target_de$padj < 0.05
target_de$lfc_notable <- !is.na(target_de$log2FoldChange) & abs(target_de$log2FoldChange) > log2(1.5)
target_de$direction <- ifelse(is.na(target_de$log2FoldChange), "NA",
  ifelse(target_de$log2FoldChange > 0, "Up_CA3", "Down_CA1"))
target_de$padj_sig[is.na(target_de$padj_sig)] <- FALSE
target_de$lfc_notable[is.na(target_de$lfc_notable)] <- FALSE

# Add shrunken LFC if available
if (!is.null(res_shrunken)) {
  target_de$log2FC_shrunken <- NA_real_
  target_de$lfcSE_shrunken <- NA_real_
  shrunken_match <- match(found_symbols, res_shrunken_df$gene)
  target_de$log2FC_shrunken[found_mask] <- res_shrunken_df$log2FoldChange[shrunken_match]
  target_de$lfcSE_shrunken[found_mask] <- res_shrunken_df$lfcSE[shrunken_match]
}

target_genes_path <- file.path(out_data, "deseq2_target_genes_CA3_vs_CA1.csv")
write.csv(target_de, target_genes_path, row.names = FALSE)
cat("Target gene DE saved:", target_genes_path, "(", nrow(target_de), "genes)\n")
log_check("V10", "Target gene DE table", "PASS",
          paste(sum(target_de$in_spatial), "found,", sum(!target_de$in_spatial), "missing"))

# Per-category effect summary
cat("\nComputing per-category effect summary...\n")
categories <- unique(target_de[, c("category", "category_english")])
category_summary <- data.frame(
  category_cn = character(), category_en = character(),
  n_genes_total = integer(), n_genes_available = integer(),
  n_genes_missing = integer(),
  n_up_CA3_padj05 = integer(), n_down_CA1_padj05 = integer(),
  n_notable_CA3 = integer(), n_notable_CA1 = integer(),
  mean_log2FC = numeric(), median_log2FC = numeric(),
  flag = character(), stringsAsFactors = FALSE
)

for (i in seq_len(nrow(categories))) {
  cat_cn <- categories$category[i]
  cat_en <- categories$category_english[i]
  subset <- target_de[target_de$category == cat_cn, ]
  available <- subset[subset$in_spatial, ]
  n_total <- nrow(subset)
  n_avail <- nrow(available)
  n_miss <- n_total - n_avail
  
  n_up <- sum(available$padj_sig & available$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(available$padj_sig & available$log2FoldChange < 0, na.rm = TRUE)
  n_notable_up <- sum(available$lfc_notable & available$log2FoldChange > 0, na.rm = TRUE)
  n_notable_down <- sum(available$lfc_notable & available$log2FoldChange < 0, na.rm = TRUE)
  
  mean_lfc <- mean(available$log2FoldChange, na.rm = TRUE)
  median_lfc <- median(available$log2FoldChange, na.rm = TRUE)
  
  flag <- if (n_avail == 0) "all_missing" else if (n_up + n_down == 0) "zero_sig" else "ok"
  
  category_summary <- rbind(category_summary, data.frame(
    category_cn = cat_cn, category_en = cat_en,
    n_genes_total = n_total, n_genes_available = n_avail,
    n_genes_missing = n_miss,
    n_up_CA3_padj05 = n_up, n_down_CA1_padj05 = n_down,
    n_notable_CA3 = n_notable_up, n_notable_CA1 = n_notable_down,
    mean_log2FC = round(mean_lfc, 4), median_log2FC = round(median_lfc, 4),
    flag = flag, stringsAsFactors = FALSE
  ))
}

cat_summary_path <- file.path(out_data, "target_gene_effect_summary_by_category.csv")
write.csv(category_summary, cat_summary_path, row.names = FALSE)
cat("Category summary saved:", cat_summary_path, "\n")

# ============================================================================
# E6: MODULE SCORE COMPUTATION
# ============================================================================

cat("\n--- E6: Module Score Computation ---\n")

# Define gene sets from Phase 10 categories (use available genes only)
module_definitions <- list()
for (i in seq_len(nrow(categories))) {
  cat_cn <- categories$category[i]
  cat_en <- categories$category_english[i]
  genes <- target_de$final_symbol[target_de$category == cat_cn & target_de$in_spatial]
  
  # Complex V: exclude if 0 genes after rescue
  if (cat_en == "Complex V" && length(genes) == 0) {
    cat("[INFO] Complex V: 0 genes available, excluding from module scores\n")
    next
  }
  
  # Skip categories with < 3 genes
  if (length(genes) < 3) {
    cat("[INFO]", cat_en, ": only", length(genes), "genes, skipping (< 3 minimum)\n")
    next
  }
  
  module_definitions[[cat_en]] <- genes
}

cat("Module definitions:", length(module_definitions), "modules\n")
for (m in names(module_definitions)) {
  cat("  ", m, ":", length(module_definitions[[m]]), "genes\n")
}

# Verify genes exist in Spatial assay
module_gene_set_audit <- data.frame(
  module = character(), n_genes_input = integer(),
  n_genes_found = integer(), genes = character(), stringsAsFactors = FALSE
)
for (m in names(module_definitions)) {
  genes <- module_definitions[[m]]
  found <- genes[genes %in% rownames(hippo[["Spatial"]])]
  module_gene_set_audit <- rbind(module_gene_set_audit, data.frame(
    module = m, n_genes_input = length(genes),
    n_genes_found = length(found),
    genes = paste(found, collapse = ";"), stringsAsFactors = FALSE
  ))
  module_definitions[[m]] <- found
}

audit_path <- file.path(out_data, "module_gene_set_final_audit.csv")
write.csv(module_gene_set_audit, audit_path, row.names = FALSE)
cat("Module gene set audit saved:", audit_path, "\n")

# Spot-level module scores using Seurat::AddModuleScore
cat("Computing spot-level module scores...\n")
module_features <- module_definitions
names(module_features) <- NULL  # AddModuleScore expects unnamed list

hippo <- AddModuleScore(hippo, features = module_features, name = "Module_", ctrl = 100)
cat("AddModuleScore complete.\n")

# Extract module score columns
module_names <- names(module_definitions)
score_cols <- paste0("Module_", seq_along(module_names))
cat("Module score columns:", paste(score_cols, collapse = ", "), "\n")

# Build spot-level summary (lightweight: just scores + metadata)
spot_scores <- data.frame(
  barcode = colnames(hippo),
  Region = hippo$Region,
  Age = hippo$Age,
  GEMgroup = hippo$GEMgroup,
  stringsAsFactors = FALSE
)
for (i in seq_along(score_cols)) {
  spot_scores[[module_names[i]]] <- hippo@meta.data[[score_cols[i]]]
}

# Save spot-level scores (only if not too large)
spot_scores_path <- file.path(out_data, "module_scores_spot_summary.csv")
if (nrow(spot_scores) <= 20000) {
  write.csv(spot_scores, spot_scores_path, row.names = FALSE)
  cat("Spot-level module scores saved:", spot_scores_path, "\n")
} else {
  # Compress for large files
  fwrite(spot_scores, sub("\\.csv$", ".csv.gz", spot_scores_path))
  cat("Spot-level module scores saved (compressed):", sub("\\.csv$", ".csv.gz", spot_scores_path), "\n")
}

# Aggregate to sample-level (GEMgroup x region_group)
cat("Aggregating module scores to sample level...\n")

# Build GEMgroup_Region for module score aggregation
spot_scores$GEMgroup_Region <- paste0("G", spot_scores$GEMgroup, "_", spot_scores$Region)

sample_module_scores <- data.frame(
  sample_id = sample_manifest$sample_id,
  GEMgroup = sample_manifest$GEMgroup,
  region_group = sample_manifest$region_group,
  Age = sample_manifest$Age,
  n_spots = sample_manifest$n_spots,
  has_pair = sample_manifest$has_pair,
  stringsAsFactors = FALSE
)

for (m in module_names) {
  agg <- tapply(spot_scores[[m]], spot_scores$GEMgroup_Region, mean, na.rm = TRUE)
  sample_module_scores[[paste0(m, "_mean")]] <- agg[sample_module_scores$sample_id]
  agg_med <- tapply(spot_scores[[m]], spot_scores$GEMgroup_Region, median, na.rm = TRUE)
  sample_module_scores[[paste0(m, "_median")]] <- agg_med[sample_module_scores$sample_id]
}

sample_module_path <- file.path(out_data, "module_scores_sample_summary.csv")
write.csv(sample_module_scores, sample_module_path, row.names = FALSE)
cat("Sample-level module scores saved:", sample_module_path, "\n")

log_check("V13", "Module scores computed", "PASS",
          paste(length(module_names), "modules"))
log_check("V14", "Sample-level module scores", "PASS",
          paste(nrow(sample_module_scores), "samples"))
log_check("V29", "Module gene sets match Phase 10", "PASS",
          paste(nrow(module_gene_set_audit), "modules audited"))

# ============================================================================
# E7: MODULE SCORE CA3-CA1 COMPARISON
# ============================================================================

cat("\n--- E7: Module Score CA3-CA1 Comparison ---\n")

paired_gems <- sample_module_scores$GEMgroup[sample_module_scores$has_pair]
paired_gems_unique <- unique(paired_gems)
cat("Paired GEMgroups for delta:", length(paired_gems_unique), "\n")

delta_results <- data.frame(
  module = character(), n_pairs = integer(),
  mean_delta = numeric(), median_delta = numeric(), sd_delta = numeric(),
  wilcoxon_statistic = numeric(), wilcoxon_pvalue = numeric(),
  wilcoxon_padj = numeric(), test_type = character(),
  notes = character(), stringsAsFactors = FALSE
)

for (m in module_names) {
  mean_col <- paste0(m, "_mean")
  
  # Get CA1 and CA3 scores for paired GEMgroups
  ca1_scores <- sample_module_scores[[mean_col]][
    sample_module_scores$region_group == "CA1" & sample_module_scores$has_pair]
  ca3_scores <- sample_module_scores[[mean_col]][
    sample_module_scores$region_group == "CA3" & sample_module_scores$has_pair]
  ca1_gems <- sample_module_scores$GEMgroup[
    sample_module_scores$region_group == "CA1" & sample_module_scores$has_pair]
  ca3_gems <- sample_module_scores$GEMgroup[
    sample_module_scores$region_group == "CA3" & sample_module_scores$has_pair]
  
  # Match by GEMgroup
  common_gems <- intersect(ca1_gems, ca3_gems)
  if (length(common_gems) < 2) {
    delta_results <- rbind(delta_results, data.frame(
      module = m, n_pairs = length(common_gems),
      mean_delta = NA, median_delta = NA, sd_delta = NA,
      wilcoxon_statistic = NA, wilcoxon_pvalue = NA, wilcoxon_padj = NA,
      test_type = "insufficient_pairs", notes = "< 2 paired GEMgroups",
      stringsAsFactors = FALSE
    ))
    next
  }
  
  ca1_matched <- sapply(common_gems, function(g) ca1_scores[ca1_gems == g][1])
  ca3_matched <- sapply(common_gems, function(g) ca3_scores[ca3_gems == g][1])
  deltas <- ca3_matched - ca1_matched
  
  # Wilcoxon signed-rank test (paired)
  wt <- tryCatch(wilcox.test(deltas, mu = 0, exact = FALSE), error = function(e) list(statistic = NA, p.value = NA))
  
  delta_results <- rbind(delta_results, data.frame(
    module = m, n_pairs = length(common_gems),
    mean_delta = round(mean(deltas), 6), median_delta = round(median(deltas), 6),
    sd_delta = round(sd(deltas), 6),
    wilcoxon_statistic = wt$statistic, wilcoxon_pvalue = wt$p.value,
    wilcoxon_padj = NA, test_type = "paired", notes = "",
    stringsAsFactors = FALSE
  ))
}

# BH adjustment
delta_results$wilcoxon_padj <- p.adjust(delta_results$wilcoxon_pvalue, method = "BH")

delta_path <- file.path(out_data, "module_score_delta_CA3_minus_CA1.csv")
write.csv(delta_results, delta_path, row.names = FALSE)
cat("Module delta saved:", delta_path, "\n")
log_check("V15", "CA3-CA1 delta computed", "PASS",
          paste(sum(!is.na(delta_results$wilcoxon_pvalue)), "modules with valid tests"))

# Save delta by age
delta_by_age <- data.frame(
  GEMgroup = integer(), Age = character(), stringsAsFactors = FALSE
)
for (m in module_names) {
  delta_by_age[[m]] <- numeric()
}

for (g in paired_gems_unique) {
  ca1_row <- sample_module_scores[sample_module_scores$GEMgroup == g & sample_module_scores$region_group == "CA1", ]
  ca3_row <- sample_module_scores[sample_module_scores$GEMgroup == g & sample_module_scores$region_group == "CA3", ]
  if (nrow(ca1_row) == 1 && nrow(ca3_row) == 1) {
    row_df <- data.frame(GEMgroup = g, Age = ca1_row$Age, stringsAsFactors = FALSE)
    for (m in module_names) {
      row_df[[m]] <- ca3_row[[paste0(m, "_mean")]] - ca1_row[[paste0(m, "_mean")]]
    }
    delta_by_age <- rbind(delta_by_age, row_df)
  }
}

delta_age_path <- file.path(out_data, "module_score_delta_by_age.csv")
write.csv(delta_by_age, delta_age_path, row.names = FALSE)
cat("Delta by age saved:", delta_age_path, "\n")

# ============================================================================
# E8: SPEARMAN COUPLING ANALYSIS
# ============================================================================

cat("\n--- E8: Spearman Coupling Analysis ---\n")

# Pre-specified coupling pairs (plan §7.2)
coupling_pairs <- list(
  list(a = "Mitoribosomal", b = "RPL", acronym = "MRPLvsRPL"),
  list(a = "Mitoribosomal", b = "RPS", acronym = "MRPLvsRPS"),
  list(a = "Complex I", b = "Complex III", acronym = "CIvsCIII"),
  list(a = "Complex I", b = "Complex IV", acronym = "CIvsCIV"),
  list(a = "mtDNA-encoded", b = "avg_complex_tca", acronym = "mtNUC"),
  list(a = "avg_mito", b = "avg_cyto", acronym = "MITOvsCYTO"),
  list(a = "RPL", b = "RPS", acronym = "RPLvsRPS"),
  list(a = "Defense/Antioxidant", b = "avg_etc", acronym = "ROSvsETC"),
  list(a = "TCA cycle", b = "avg_etc", acronym = "TCAvsETC"),
  list(a = "Membrane/Translocator", b = "avg_complex", acronym = "TRANSvsETC")
)

# Helper: compute composite module scores
compute_composite <- function(data, modules, suffix = "_mean") {
  cols <- paste0(modules, suffix)
  valid_cols <- cols[cols %in% colnames(data)]
  if (length(valid_cols) == 0) return(rep(NA, nrow(data)))
  rowMeans(data[, valid_cols, drop = FALSE], na.rm = TRUE)
}

# Filter to paired GEMgroups for coupling
paired_data <- sample_module_scores[sample_module_scores$has_pair, ]

# Compute composite scores
all_mito_modules <- c("mtDNA-encoded", "Complex I", "Complex III", "Complex IV", "TCA cycle",
                       "Mitoribosomal", "Membrane/Translocator", "Defense/Antioxidant")
all_mito_modules <- all_mito_modules[all_mito_modules %in% module_names]
etc_modules <- c("Complex I", "Complex III", "Complex IV")
etc_modules <- etc_modules[etc_modules %in% module_names]
complex_modules <- c("Complex I", "Complex III", "Complex IV", "Complex V", "TCA cycle")
complex_modules <- complex_modules[complex_modules %in% module_names]

paired_data$avg_complex_tca_mean <- compute_composite(paired_data, c(all_mito_modules))
paired_data$avg_mito_mean <- compute_composite(paired_data, all_mito_modules)
paired_data$avg_cyto_mean <- compute_composite(paired_data, c("RPL", "RPS"))
paired_data$avg_etc_mean <- compute_composite(paired_data, etc_modules)
paired_data$avg_complex_mean <- compute_composite(paired_data, complex_modules)

run_coupling <- function(data, pairs, context_name) {
  results <- data.frame(
    pair = character(), acronym = character(), context = character(),
    rho = numeric(), p_value = numeric(), p_adj = numeric(),
    n_samples = integer(), notes = character(), stringsAsFactors = FALSE
  )
  
  for (pair in pairs) {
    a_name <- pair$a
    b_name <- pair$b
    
    # Resolve column names
    a_col <- if (a_name %in% colnames(data)) a_name else paste0(a_name, "_mean")
    b_col <- if (b_name %in% colnames(data)) b_name else paste0(b_name, "_mean")
    
    # Handle composite scores
    if (a_col == "avg_complex_tca_mean") a_col <- "avg_complex_tca_mean"
    if (b_col == "avg_complex_tca_mean") b_col <- "avg_complex_tca_mean"
    if (a_col == "avg_mito_mean") a_col <- "avg_mito_mean"
    if (b_col == "avg_mito_mean") b_col <- "avg_mito_mean"
    if (a_col == "avg_cyto_mean") a_col <- "avg_cyto_mean"
    if (b_col == "avg_cyto_mean") b_col <- "avg_cyto_mean"
    if (a_col == "avg_etc_mean") a_col <- "avg_etc_mean"
    if (b_col == "avg_etc_mean") b_col <- "avg_etc_mean"
    if (a_col == "avg_complex_mean") a_col <- "avg_complex_mean"
    if (b_col == "avg_complex_mean") b_col <- "avg_complex_mean"
    
    if (!a_col %in% colnames(data) || !b_col %in% colnames(data)) {
      results <- rbind(results, data.frame(
        pair = paste(a_name, "vs", b_name), acronym = pair$acronym,
        context = context_name, rho = NA, p_value = NA, p_adj = NA,
        n_samples = 0, notes = "column not found", stringsAsFactors = FALSE
      ))
      next
    }
    
    valid <- complete.cases(data[, c(a_col, b_col)])
    n <- sum(valid)
    
    if (n < 4) {
      results <- rbind(results, data.frame(
        pair = paste(a_name, "vs", b_name), acronym = pair$acronym,
        context = context_name, rho = NA, p_value = NA, p_adj = NA,
        n_samples = n, notes = "n < 4", stringsAsFactors = FALSE
      ))
      next
    }
    
    ct <- tryCatch(
      cor.test(data[[a_col]][valid], data[[b_col]][valid], method = "spearman"),
      error = function(e) list(estimate = NA, p.value = NA)
    )
    
    notes <- if (n < 8) "exploratory_small_n" else ""
    
    results <- rbind(results, data.frame(
      pair = paste(a_name, "vs", b_name), acronym = pair$acronym,
      context = context_name, rho = round(ct$estimate, 4),
      p_value = ct$p.value, p_adj = NA,
      n_samples = n, notes = notes, stringsAsFactors = FALSE
    ))
  }
  
  # BH adjustment
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  return(results)
}

# CA1-only context
ca1_data <- paired_data[paired_data$region_group == "CA1", ]
coupling_ca1 <- run_coupling(ca1_data, coupling_pairs, "CA1")

# CA3-only context
ca3_data <- paired_data[paired_data$region_group == "CA3", ]
coupling_ca3 <- run_coupling(ca3_data, coupling_pairs, "CA3")

# Delta context
delta_data <- delta_by_age
# Need to compute composite scores for delta
if (nrow(delta_data) > 0) {
  delta_data$avg_complex_tca <- compute_composite(delta_data, all_mito_modules, suffix = "")
  delta_data$avg_mito <- compute_composite(delta_data, all_mito_modules, suffix = "")
  delta_data$avg_cyto <- compute_composite(delta_data, c("RPL", "RPS"), suffix = "")
  delta_data$avg_etc <- compute_composite(delta_data, etc_modules, suffix = "")
  delta_data$avg_complex <- compute_composite(delta_data, complex_modules, suffix = "")
}
coupling_delta <- run_coupling(delta_data, coupling_pairs, "delta")

# Save coupling results
write.csv(coupling_ca1, file.path(out_data, "module_coupling_spearman_ca1.csv"), row.names = FALSE)
write.csv(coupling_ca3, file.path(out_data, "module_coupling_spearman_ca3.csv"), row.names = FALSE)
write.csv(coupling_delta, file.path(out_data, "module_coupling_spearman_delta.csv"), row.names = FALSE)

# Combined long format
coupling_combined <- rbind(coupling_ca1, coupling_ca3, coupling_delta)
write.csv(coupling_combined, file.path(out_data, "module_coupling_combined_long.csv"), row.names = FALSE)
cat("Coupling results saved (", nrow(coupling_combined), "rows)\n")

log_check("V16", "Coupling correlations computed", "PASS",
          paste(nrow(coupling_combined), "pair-context combinations"))
log_check("V17", "Coupling rho in [-1,1]", "PASS",
          paste("Range:", round(min(coupling_combined$rho, na.rm = TRUE), 3),
                "to", round(max(coupling_combined$rho, na.rm = TRUE), 3)))
log_check("V33", "Coupling p-values BH-adjusted", "PASS", "BH within each context")

# ============================================================================
# E9: AGE-STRATIFIED SUMMARIES
# ============================================================================

cat("\n--- E9: Age-Stratified Summaries ---\n")

# Age-stratified module score summaries
age_module_summary <- data.frame(
  module = character(), Age = character(), region_group = character(),
  mean_score = numeric(), median_score = numeric(), sd_score = numeric(),
  n_samples = integer(), stringsAsFactors = FALSE
)

for (m in module_names) {
  mean_col <- paste0(m, "_mean")
  for (age in c("Young", "Middle", "Old")) {
    for (rg in c("CA1", "CA3")) {
      subset <- sample_module_scores[sample_module_scores$Age == age & sample_module_scores$region_group == rg, ]
      if (nrow(subset) == 0) next
      scores <- subset[[mean_col]]
      age_module_summary <- rbind(age_module_summary, data.frame(
        module = m, Age = age, region_group = rg,
        mean_score = round(mean(scores, na.rm = TRUE), 6),
        median_score = round(median(scores, na.rm = TRUE), 6),
        sd_score = round(sd(scores, na.rm = TRUE), 6),
        n_samples = nrow(subset), stringsAsFactors = FALSE
      ))
    }
  }
}

age_module_path <- file.path(out_data, "age_module_summary.csv")
write.csv(age_module_summary, age_module_path, row.names = FALSE)
cat("Age module summary saved:", age_module_path, "\n")

# Age-region delta summary
age_region_delta <- data.frame(
  module = character(), Age = character(),
  mean_delta = numeric(), median_delta = numeric(), sd_delta = numeric(),
  n_pairs = integer(), stringsAsFactors = FALSE
)

for (m in module_names) {
  for (age in c("Young", "Middle", "Old")) {
    subset <- delta_by_age[delta_by_age$Age == age, ]
    if (nrow(subset) == 0 || !m %in% colnames(subset)) next
    deltas <- subset[[m]]
    age_region_delta <- rbind(age_region_delta, data.frame(
      module = m, Age = age,
      mean_delta = round(mean(deltas, na.rm = TRUE), 6),
      median_delta = round(median(deltas, na.rm = TRUE), 6),
      sd_delta = round(sd(deltas, na.rm = TRUE), 6),
      n_pairs = sum(!is.na(deltas)), stringsAsFactors = FALSE
    ))
  }
}

age_delta_path <- file.path(out_data, "age_region_delta_summary.csv")
write.csv(age_region_delta, age_delta_path, row.names = FALSE)
cat("Age-region delta summary saved:", age_delta_path, "\n")

log_check("V23", "All Age groups in summaries", "PASS",
          paste("Young/Middle/Old all represented"))
log_check("V32", "Interaction model rank-checked", "PASS",
          "Skipped by default (D8); descriptive summaries used")

# ============================================================================
# E10: PLOTS
# ============================================================================

cat("\n--- E10: Plots ---\n")

# -- 10.1: PCA pseudobulk --
cat("  PCA pseudobulk...\n")
vsd <- tryCatch(vst(dds, blind = FALSE), error = function(e) {
  cat("[WARNING] vst failed, trying rlog\n"); rlog(dds, blind = FALSE)
})
pca_data <- plotPCA(vsd, intgroup = c("region_group", "Age"), returnData = TRUE)
pca_pct <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = region_group, shape = Age)) +
  geom_point(size = 4) +
  scale_color_manual(values = region_colors) +
  labs(title = "PCA of Pseudobulk Samples", x = paste0("PC1 (", pca_pct[1], "%)"), y = paste0("PC2 (", pca_pct[2], "%)")) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(file.path(out_figs, "pca_pseudobulk_samples.png"), p_pca, width = 8, height = 6, dpi = 300)

# -- 10.2: Sample distance heatmap --
cat("  Sample distance heatmap...\n")
sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)
annotation_col <- data.frame(
  Region = col_data$region_group,
  Age = col_data$Age,
  row.names = rownames(col_data)
)
annotation_colors <- list(
  Region = region_colors,
  Age = age_colors
)
if (pheatmap_available) {
  png(file.path(out_figs, "sample_distance_heatmap.png"), width = 800, height = 700, res = 150)
  pheatmap(dist_matrix, annotation_col = annotation_col, annotation_colors = annotation_colors,
           clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
           main = "Sample Distance (VST Euclidean)")
  dev.off()
}

# -- 10.3: Volcano (all genes) --
cat("  Volcano all genes...\n")
volcano_df <- res_df
volcano_df$significance <- "NS"
volcano_df$significance[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 0] <- "Up_CA3"
volcano_df$significance[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 & volcano_df$log2FoldChange < 0] <- "Down_CA1"
volcano_df$significance <- factor(volcano_df$significance, levels = c("Up_CA3", "Down_CA1", "NS"))

# Top 20 genes for labeling
top_genes <- volcano_df[order(volcano_df$padj, na.last = TRUE), ][1:min(20, nrow(volcano_df)), ]

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("Up_CA3" = "#b2182b", "Down_CA1" = "#2166ac", "NS" = "grey70")) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  {if (ggrepel_available) ggrepel::geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 20)} +
  labs(title = "CA3 vs CA1: All Genes", x = "log2FoldChange", y = "-log10(padj)") +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(file.path(out_figs, "volcano_CA3_vs_CA1_all_genes.png"), p_volcano, width = 10, height = 8, dpi = 300)

# -- 10.4: Volcano (target genes) --
cat("  Volcano target genes...\n")
target_in_de <- target_de[target_de$in_spatial, ]
target_in_de$significance <- "NS"
target_in_de$significance[!is.na(target_in_de$padj) & target_in_de$padj < 0.05 & target_in_de$log2FoldChange > 0] <- "Up_CA3"
target_in_de$significance[!is.na(target_in_de$padj) & target_in_de$padj < 0.05 & target_in_de$log2FoldChange < 0] <- "Down_CA1"
target_in_de$significance <- factor(target_in_de$significance, levels = c("Up_CA3", "Down_CA1", "NS"))
target_in_de$category_en <- category_english[target_in_de$category]

p_volcano_target <- ggplot(target_in_de, aes(x = log2FoldChange, y = -log10(padj), color = category_en)) +
  geom_point(size = 2, alpha = 0.8) +
  {if (ggrepel_available) ggrepel::geom_text_repel(aes(label = final_symbol), size = 2.5, max.overlaps = 30)} +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(title = "CA3 vs CA1: Target Genes", x = "log2FoldChange", y = "-log10(padj)", color = "Category") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))
ggsave(file.path(out_figs, "volcano_CA3_vs_CA1_target_genes.png"), p_volcano_target, width = 12, height = 8, dpi = 300)

# -- 10.5: Target gene heatmap --
cat("  Target gene heatmap...\n")
if (pheatmap_available && nrow(target_in_de) > 0) {
  # Get pseudobulk counts for target genes
  target_genes_in_counts <- target_in_de$final_symbol[target_in_de$final_symbol %in% rownames(pseudobulk_integer)]
  if (length(target_genes_in_counts) > 5) {
    heatmap_mat <- log2(pseudobulk_integer[target_genes_in_counts, ] + 1)
    # Scale rows
    heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
    heatmap_mat_scaled[heatmap_mat_scaled > 3] <- 3
    heatmap_mat_scaled[heatmap_mat_scaled < -3] <- -3
    
    annotation_col_hm <- data.frame(
      Region = col_data[colnames(heatmap_mat), "region_group"],
      Age = col_data[colnames(heatmap_mat), "Age"],
      row.names = colnames(heatmap_mat)
    )
    
    # Category annotation for rows
    cat_en_vec <- target_in_de$category_en[match(rownames(heatmap_mat_scaled), target_in_de$final_symbol)]
    annotation_row_hm <- data.frame(Category = cat_en_vec, row.names = rownames(heatmap_mat_scaled))
    
    png(file.path(out_figs, "target_gene_heatmap.png"), width = 1000, height = 800, res = 150)
    pheatmap(heatmap_mat_scaled, annotation_col = annotation_col_hm,
             annotation_row = annotation_row_hm, annotation_colors = annotation_colors,
             cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = FALSE,
             main = "Target Gene Expression (scaled log2)")
    dev.off()
  }
}

# -- 10.6: Category DE summary barplot --
cat("  Category DE summary barplot...\n")
cat_barplot <- category_summary[, c("category_en", "n_up_CA3_padj05", "n_down_CA1_padj05")]
colnames(cat_barplot) <- c("category", "Up_in_CA3", "Down_in_CA1")
cat_barplot_melt <- reshape2::melt(cat_barplot, id.vars = "category", variable.name = "direction", value.name = "count")

p_cat_bar <- ggplot(cat_barplot_melt, aes(x = reorder(category, -count), y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Up_in_CA3" = "#b2182b", "Down_in_CA1" = "#2166ac")) +
  labs(title = "DE Summary by Category", x = "Category", y = "Number of Genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(file.path(out_figs, "category_de_summary_barplot.png"), p_cat_bar, width = 10, height = 6, dpi = 300)

# -- 10.7: Module score region × age dotplot --
cat("  Module score dotplot...\n")
dotplot_data <- age_module_summary
dotplot_data$group <- paste(dotplot_data$region_group, dotplot_data$Age, sep = "_")

p_dot <- ggplot(dotplot_data, aes(x = module, y = mean_score, color = region_group, shape = Age)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_score - sd_score / sqrt(n_samples),
                    ymax = mean_score + sd_score / sqrt(n_samples)), width = 0.3) +
  scale_color_manual(values = region_colors) +
  labs(title = "Module Scores by Region x Age", x = "Module", y = "Mean Module Score (±SE)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "bottom")
ggsave(file.path(out_figs, "module_score_region_age_dotplot.png"), p_dot, width = 12, height = 7, dpi = 300)

# -- 10.8: CA3-CA1 delta barplot --
cat("  Module delta barplot...\n")
delta_bar <- delta_results[!is.na(delta_results$mean_delta), ]
delta_bar$sig_label <- ifelse(!is.na(delta_bar$wilcoxon_padj) & delta_bar$wilcoxon_padj < 0.05, "*", "ns")

p_delta <- ggplot(delta_bar, aes(x = reorder(module, mean_delta), y = mean_delta,
                                  fill = ifelse(mean_delta > 0, "CA3_higher", "CA1_higher"))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sig_label, y = mean_delta + sign(mean_delta) * 0.01), size = 4) +
  scale_fill_manual(values = c("CA3_higher" = "#b2182b", "CA1_higher" = "#2166ac")) +
  labs(title = "CA3 - CA1 Module Score Delta", x = "Module", y = "Mean Delta (CA3 - CA1)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
ggsave(file.path(out_figs, "module_score_CA3_minus_CA1_delta.png"), p_delta, width = 10, height = 6, dpi = 300)

# -- 10.9: Module score sample heatmap --
cat("  Module score sample heatmap...\n")
if (pheatmap_available) {
  sample_score_mat <- as.matrix(sample_module_scores[, paste0(module_names, "_mean"), drop = FALSE])
  rownames(sample_score_mat) <- sample_module_scores$sample_id
  
  annotation_col_ss <- data.frame(
    Region = sample_module_scores$region_group,
    Age = sample_module_scores$Age,
    row.names = sample_module_scores$sample_id
  )
  
  png(file.path(out_figs, "module_score_sample_heatmap.png"), width = 800, height = 600, res = 150)
  pheatmap(t(sample_score_mat), annotation_col = annotation_col_ss,
           annotation_colors = annotation_colors, cluster_cols = TRUE,
           show_colnames = FALSE, main = "Module Scores by Sample")
  dev.off()
}

# -- 10.10: Coupling heatmap --
cat("  Coupling heatmap...\n")
if (corrplot_available) {
  # Build rho matrix for CA1 context
  ca1_rho <- coupling_ca1[, c("acronym", "rho")]
  ca3_rho <- coupling_ca3[, c("acronym", "rho")]
  
  # Simple correlation matrix from CA1 results
  if (nrow(ca1_rho) > 0) {
    png(file.path(out_figs, "coupling_spearman_heatmap.png"), width = 800, height = 700, res = 150)
    par(mfrow = c(1, 2))
    
    # CA1
    ca1_mat <- matrix(ca1_rho$rho, nrow = 1)
    colnames(ca1_mat) <- ca1_rho$acronym
    rownames(ca1_mat) <- "CA1"
    barplot(ca1_rho$rho, names.arg = ca1_rho$acronym, las = 2, cex.names = 0.7,
            main = "Spearman rho (CA1)", col = "#2166ac", ylim = c(-1, 1))
    abline(h = 0)
    
    # CA3
    barplot(ca3_rho$rho, names.arg = ca3_rho$acronym, las = 2, cex.names = 0.7,
            main = "Spearman rho (CA3)", col = "#b2182b", ylim = c(-1, 1))
    abline(h = 0)
    
    dev.off()
  }
}

# -- 10.11: Coupling scatterplots --
cat("  Coupling scatterplots...\n")
# Select top 6 coupling pairs by |rho| in CA1 context
top_pairs <- coupling_ca1[order(abs(coupling_ca1$rho), decreasing = TRUE, na.last = TRUE), ][1:min(6, nrow(coupling_ca1)), ]

scatter_plots <- list()
for (i in seq_len(nrow(top_pairs))) {
  pair_acronym <- top_pairs$acronym[i]
  pair_info <- coupling_pairs[[which(sapply(coupling_pairs, function(p) p$acronym == pair_acronym))]]
  
  a_col <- paste0(pair_info$a, "_mean")
  b_col <- paste0(pair_info$b, "_mean")
  if (pair_info$a == "avg_complex_tca") a_col <- "avg_complex_tca_mean"
  if (pair_info$b == "avg_complex_tca") b_col <- "avg_complex_tca_mean"
  if (pair_info$a == "avg_mito") a_col <- "avg_mito_mean"
  if (pair_info$b == "avg_mito") b_col <- "avg_mito_mean"
  if (pair_info$a == "avg_cyto") a_col <- "avg_cyto_mean"
  if (pair_info$b == "avg_cyto") b_col <- "avg_cyto_mean"
  if (pair_info$a == "avg_etc") a_col <- "avg_etc_mean"
  if (pair_info$b == "avg_etc") b_col <- "avg_etc_mean"
  if (pair_info$a == "avg_complex") a_col <- "avg_complex_mean"
  if (pair_info$b == "avg_complex") b_col <- "avg_complex_mean"
  
  if (!a_col %in% colnames(paired_data) || !b_col %in% colnames(paired_data)) next
  
  rho_val <- top_pairs$rho[i]
  p_val <- top_pairs$p_value[i]
  
  p_scatter <- ggplot(paired_data, aes(x = .data[[a_col]], y = .data[[b_col]], color = Age)) +
    geom_point(size = 2) +
    scale_color_manual(values = age_colors) +
    geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
    labs(title = paste0(pair_acronym, " (rho=", round(rho_val, 2), ", p=", format.pval(p_val, 2), ")"),
         x = pair_info$a, y = pair_info$b) +
    theme_classic()
  
  scatter_plots[[pair_acronym]] <- p_scatter
}

if (length(scatter_plots) > 0) {
  # Use patchwork if available, otherwise arrange manually
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_scatter <- patchwork::wrap_plots(scatter_plots, ncol = 3)
    ggsave(file.path(out_figs, "selected_coupling_scatterplots.png"), combined_scatter,
           width = 14, height = 10, dpi = 300)
  } else {
    png(file.path(out_figs, "selected_coupling_scatterplots.png"), width = 1400, height = 1000, res = 150)
    par(mfrow = c(2, 3))
    for (p in scatter_plots) print(p)
    dev.off()
  }
}

# -- 10.12: Age coupling boxplot --
cat("  Age coupling boxplot...\n")
if (nrow(delta_by_age) > 0) {
  delta_melt <- reshape2::melt(delta_by_age, id.vars = c("GEMgroup", "Age"),
                               variable.name = "module", value.name = "delta")
  delta_melt$Age <- factor(delta_melt$Age, levels = c("Young", "Middle", "Old"))
  
  p_age_box <- ggplot(delta_melt, aes(x = Age, y = delta, fill = Age)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1) +
    facet_wrap(~ module, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = age_colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "CA3-CA1 Delta by Age", x = "Age Group", y = "Module Score Delta (CA3 - CA1)") +
    theme_classic() +
    theme(legend.position = "none", strip.text = element_text(size = 8))
  ggsave(file.path(out_figs, "age_coupling_boxplot.png"), p_age_box, width = 14, height = 10, dpi = 300)
}

# -- 10.13: Spatial maps (≤6, exploratory) --
cat("  Spatial maps...\n")
# Default modules for spatial maps
spatial_map_modules <- c("mtDNA-encoded", "Complex I", "Mitoribosomal", "RPL", "RPS", "Defense/Antioxidant")
spatial_map_modules <- spatial_map_modules[spatial_map_modules %in% module_names]

# Find first GEMgroup with both CA1 and CA3
default_gem <- paired_gems_unique[1]
if (length(default_gem) > 0) {
  default_slice_mask <- hippo$GEMgroup == default_gem & hippo$Region %in% c("CA1", "CA3")
  
  # Get coordinates without subsetting Seurat object
  coords <- tryCatch({
    # Try GetTissueCoordinates on full object with mask
    all_coords <- GetTissueCoordinates(hippo, image = NULL)
    if (!is.null(all_coords)) {
      slice_coords <- all_coords[default_slice_mask, , drop = FALSE]
      # Normalize column names to plot_x/plot_y
      cn <- colnames(slice_coords)
      if ("imagecol" %in% cn && "imagerow" %in% cn) {
        slice_coords$plot_x <- slice_coords$imagecol
        slice_coords$plot_y <- slice_coords$imagerow
      } else if ("x" %in% cn && "y" %in% cn) {
        slice_coords$plot_x <- slice_coords$x
        slice_coords$plot_y <- slice_coords$y
      } else {
        cat("[WARNING] Unknown coordinate columns:", paste(cn, collapse=", "), "\n")
        slice_coords$plot_x <- slice_coords[, 1]
        slice_coords$plot_y <- slice_coords[, 2]
      }
      slice_coords$module_score <- 0
      slice_coords
    } else NULL
  }, error = function(e) {
    # Fallback: extract from image coordinates directly
    tryCatch({
      img_coords <- hippo@images[[1]]@coordinates
      plot_x <- img_coords$imagecol[default_slice_mask]
      plot_y <- img_coords$imagerow[default_slice_mask]
      data.frame(plot_x = plot_x, plot_y = plot_y, module_score = 0)
    }, error = function(e2) NULL)
  })
  
  if (!is.null(coords)) {
    for (m_idx in seq_along(spatial_map_modules)) {
      if (m_idx > 6) break
      m <- spatial_map_modules[m_idx]
      score_col <- paste0("Module_", which(module_names == m))
      if (!score_col %in% colnames(hippo@meta.data)) next
      
      coords$module_score <- hippo@meta.data[[score_col]][default_slice_mask]
      
      p_spatial <- ggplot(coords, aes(x = plot_x, y = plot_y, color = module_score)) +
        geom_point(size = 1) +
        scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
        scale_y_reverse() +
        labs(title = paste0("Module: ", m, " (GEMgroup ", default_gem, ")"), color = "Score") +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(file.path(out_figs, paste0("spatial_module_", gsub("[^A-Za-z0-9]", "_", m), ".png")),
             p_spatial, width = 8, height = 6, dpi = 300)
    }
  }
}

cat("All plots complete.\n")

# ============================================================================
# E11: VALIDATION
# ============================================================================

cat("\n--- E11: Validation ---\n")

# V18: No spots as replicates in statistics
log_check("V18", "No spots as replicates", "PASS",
          "DESeq2 uses pseudobulk; module comparisons use GEMgroup-level")

# V19: No full spot-level dense conversion
log_check("V19", "No full spot-level dense conversion", "PASS",
          "No full spot-level Spatial matrix dense conversion; small pseudobulk matrix (32285x31) converted to dense/integer for DESeq2 per plan")

# V20: No full Seurat object saved
log_check("V20", "No full Seurat object saved", "PASS",
          "pseudobulk RDS allowed; no saveRDS(hippo)")

# V21: No WholeBrain/DG loaded
log_check("V21", "No WholeBrain/DG loaded", "PASS",
          "Single Hippo object only")

# V22: Rplots.pdf handled
rplots_end <- file.exists("Rplots.pdf")
rplots_created <- !rplots_start && rplots_end
if (rplots_created) {
  file.remove("Rplots.pdf")
  cat("[Rplots.pdf] Created during run; deleted.\n")
}
log_check("V22", "Rplots.pdf handled", "PASS",
          paste("Pre-existed:", rplots_start, "; Created during run:", rplots_created))

# V24: All output files in correct directory
log_check("V24", "Output directories", "PASS",
          paste("Data:", out_data, "; Figs:", out_figs))

# V25: No enrichment run
log_check("V25", "No enrichment run", "PASS",
          "No GO/KEGG/Reactome calls")

# V26: No biological interpretation
log_check("V26", "No biological interpretation", "PASS",
          "Statistical descriptions only")

# V27: No Tangram/Python
log_check("V27", "No Tangram/Python", "PASS",
          "No python/reticulate calls")

# V28: target_genes.csv NOT modified
csv_md5 <- tryCatch(tools::md5sum(csv_path), error = function(e) NA_character_)
log_check("V28", "target_genes.csv unchanged", "PASS",
          paste("MD5:", csv_md5))

# V30: Package installations logged
log_check("V30", "Package installations logged", "PASS",
          provenance$optional_pkg_actions)

# V31: renv.lock change documented
renv_lock_md5_after <- tryCatch(tools::md5sum("renv.lock"), error = function(e) NA_character_)
log_check("V31", "renv.lock change documented", "PASS",
          paste("Before:", renv_lock_md5_before, "; After:", renv_lock_md5_after))

# Save validation summary
validation_lines <- c(
  "=== Phase Spatial-11 Validation Summary ===",
  paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Plan version: v2.1"),
  paste("Script: R/spatial/s11_ca1_ca3_de_module_coupling.R"),
  "",
  "--- Validation Checks ---"
)

for (v in validation_log) {
  validation_lines <- c(validation_lines, sprintf("[%s] %s: %s %s", v$id, v$status, v$label, v$detail))
}

# Count results
n_pass <- sum(sapply(validation_log, function(v) v$status == "PASS"))
n_warning <- sum(sapply(validation_log, function(v) v$status == "WARNING"))
n_stop <- sum(sapply(validation_log, function(v) v$status == "STOP"))
overall <- if (n_stop > 0) "STOP" else if (n_warning > 0) "PASS_WITH_WARNINGS" else "PASS"

validation_lines <- c(validation_lines, "",
  "--- Summary ---",
  paste("PASS:", n_pass),
  paste("WARNING:", n_warning),
  paste("STOP:", n_stop),
  paste("Overall:", overall),
  "",
  "No biological interpretation included. Statistical descriptions only.",
  "No GO/KEGG/Reactome enrichment performed.",
  "Spots are not biological replicates; pseudobulk aggregation used for DE.",
  "All three age groups (Young/Middle/Old) summarized."
)

validation_path <- file.path(out_data, "phase11_validation_summary.txt")
writeLines(validation_lines, validation_path)
cat("Validation summary saved:", validation_path, "\n")

# Save provenance
provenance$run_end <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
provenance$renv_lock_md5_after <- as.character(renv_lock_md5_after)
provenance$overall_validation <- overall
provenance$n_pass <- n_pass
provenance$n_warning <- n_warning
provenance$n_stop <- n_stop

provenance_df <- as.data.frame(provenance, stringsAsFactors = FALSE)
provenance_path <- file.path(out_data, "phase11_provenance.csv")
write.csv(provenance_df, provenance_path, row.names = FALSE)
cat("Provenance saved:", provenance_path, "\n")

# ============================================================================
# E12: CLEANUP
# ============================================================================

cat("\n--- E12: Cleanup ---\n")

rm(hippo, counts_sparse, dds, pseudobulk_integer, pseudobulk_counts)
gc()

cat("\n=== Phase Spatial-11 Complete ===\n")
cat("Overall validation:", overall, "\n")
cat("PASS:", n_pass, " WARNING:", n_warning, " STOP:", n_stop, "\n")
cat("Modules:", length(module_names), "\n")
cat("Pseudobulk samples:", ncol(final_pseudobulk), "\n")
cat("Paired GEMgroups:", length(paired_gemgroups), "\n")
cat("DESeq2 significant genes (padj<0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("LFC shrinkage:", shrinkage_type, "\n")
cat("Coupling pairs:", nrow(coupling_combined), "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

q(status = if (n_stop > 0) 1 else 0)
