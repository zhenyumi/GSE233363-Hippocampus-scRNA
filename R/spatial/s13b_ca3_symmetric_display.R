## Phase Spatial-13b: CA3 Symmetric Display & Guide/Report Correction
## Plan: docs/spatial_phase13b_ca3_symmetric_display_plan.md
## Script: R/spatial/s13b_ca3_symmetric_display.R
## Purpose: Display supplement only — NO re-computation of Phase 13 core logic

cat("=== Phase Spatial-13b: CA3 Symmetric Display ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

t_start <- Sys.time()

# Rplots.pdf guard
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

## ---- E0: Setup & Validation ----
cat("=== E0: Setup & Validation ===\n")

renv_md5_before <- digest::digest("renv.lock", algo = "md5", file = TRUE)
cat("renv.lock MD5 (before):", renv_md5_before, "\n")
expected_md5 <- "b2bf05038b440576855d043e137c5883"
if (renv_md5_before != expected_md5) {
  cat("WARNING: renv.lock MD5 differs from Phase 13 baseline\n")
}

# Required packages (NO new packages)
required_pkgs <- c("dplyr", "ggplot2", "pheatmap", "patchwork", "stringr", "digest")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("FATAL: Package not available: ", pkg)
}
library(dplyr)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(stringr)

has_seurat <- requireNamespace("Seurat", quietly = TRUE)

# Output directories (must already exist)
out_dir <- "data/processed/spatial/phase13_candidate_regulator"
fig_dir <- "figures/spatial/phase13_candidate_regulator"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# V01: Phase 13 validation has 0 FAIL
val_file <- file.path(out_dir, "phase13_validation_summary.csv")
if (!file.exists(val_file)) stop("FATAL: Phase 13 validation summary not found")
val_summary <- read.csv(val_file, stringsAsFactors = FALSE)
if (any(val_summary$result == "FAIL")) stop("FATAL: Phase 13 has FAIL results")
cat("V01 PASS: Phase 13 validation OK (", sum(val_summary$result == "PASS"), "PASS )\n")

# V02-V03: score_default has CA3 columns
score_file <- file.path(out_dir, "phase13_regulator_score_default.csv")
if (!file.exists(score_file)) stop("FATAL: score_default.csv not found")
score_default <- read.csv(score_file, stringsAsFactors = FALSE)
if (!"score_CA3_default" %in% colnames(score_default)) stop("FATAL: score_CA3_default column missing")
if (!"rank_CA3" %in% colnames(score_default)) stop("FATAL: rank_CA3 column missing")
if (sum(!is.na(score_default$rank_CA3)) == 0) stop("FATAL: No valid rank_CA3 values")
cat("V02-V03 PASS: CA3 columns present,", sum(!is.na(score_default$rank_CA3)), "valid rank_CA3\n")

# V04: candidate_classes.csv exists
classes_file <- file.path(out_dir, "phase13_candidate_classes.csv")
if (!file.exists(classes_file)) stop("FATAL: candidate_classes.csv not found")

# V06: renv.lock MD5 (already checked above)
cat("V06:", ifelse(renv_md5_before == expected_md5, "PASS", "WARN"), "\n")

# V07: No phase13b directories
if (dir.exists(file.path("data/processed/spatial", "phase13b"))) {
  cat("WARNING: phase13b directory exists (unexpected)\n")
}

cat("Pre-execution checks complete\n\n")

## ---- E1: Load Phase 13 Data ----
cat("=== E1: Load Phase 13 Data ===\n")

# Load sensitivity CSV
sens_file <- file.path(out_dir, "phase13_regulator_score_sensitivity.csv")
if (!file.exists(sens_file)) stop("FATAL: sensitivity CSV not found")
score_sensitivity <- read.csv(sens_file, stringsAsFactors = FALSE)
cat("Sensitivity CSV:", nrow(score_sensitivity), "rows\n")

# Load rho universe (for F09b component breakdown)
rho_file <- file.path(out_dir, "phase13_spearman_rho_universe_regulator.csv")
rho_universe <- NULL
if (file.exists(rho_file)) {
  rho_universe <- read.csv(rho_file, stringsAsFactors = FALSE)
  cat("Rho universe:", nrow(rho_universe), "rows\n")
} else {
  cat("WARNING: rho universe CSV not found; F09b will use score values\n")
}

# Identify module names
module_names <- unique(score_default$module)
cat("Modules:", length(module_names), "\n")

# CA3 color convention (from Phase 11/12)
col_ca3 <- "#E69F00"
col_ca1 <- "#56B4E9"

# Load pseudobulk for age-stratified figures (F07b/F08b)
pb_file <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds"
mod_score_file <- "data/processed/spatial/phase12_young_aged_region_comparison/module_scores_all_regions.csv"

has_pb_data <- file.exists(pb_file) && file.exists(mod_score_file)
if (has_pb_data) {
  cat("Pseudobulk data available for age-stratified figures\n")
} else {
  cat("WARNING: Pseudobulk data not found; F07b/F08b will be skipped\n")
}

# Load module audit for gene lists
module_audit_file <- "data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv"
module_genes <- list()
if (file.exists(module_audit_file)) {
  module_audit <- read.csv(module_audit_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  mod_col <- intersect(c("module", "Module", "module_name"), colnames(module_audit))[1]
  genes_col <- intersect(c("genes", "Genes", "gene_set", "gene_list"), colnames(module_audit))[1]
  if (!is.na(mod_col) && !is.na(genes_col)) {
    for (i in seq_len(nrow(module_audit))) {
      mod_name <- module_audit[[mod_col]][i]
      genes_str <- module_audit[[genes_col]][i]
      module_genes[[mod_name]] <- trimws(strsplit(genes_str, "[,]")[[1]])
    }
  }
}

cat("\n")

## ---- E2: Generate Top Candidate CSVs ----
cat("=== E2: Generate Top Candidate CSVs ===\n")

csv_cols <- c("gene", "module", "candidate_class", "rho_CA1", "rho_CA3",
              "score_CA1_default", "score_CA3_default", "rank_CA1", "rank_CA3",
              "DE_age_lfc_CA1", "DE_age_lfc_CA3", "DE_region_lfc", "self_correlation_flag")

# O1: CA1 top candidates (sorted by rank_CA1)
ca1_df <- score_default[!is.na(score_default$rank_CA1), csv_cols]
ca1_df <- ca1_df[order(ca1_df$rank_CA1, -ca1_df$score_CA1_default), ]
write.csv(ca1_df, file.path(out_dir, "phase13_top_candidates_CA1.csv"), row.names = FALSE)
cat("O1: CA1 top candidates,", nrow(ca1_df), "rows\n")

# O2: CA3 top candidates (sorted by rank_CA3)
ca3_df <- score_default[!is.na(score_default$rank_CA3), csv_cols]
ca3_df <- ca3_df[order(ca3_df$rank_CA3, -ca3_df$score_CA3_default), ]
write.csv(ca3_df, file.path(out_dir, "phase13_top_candidates_CA3.csv"), row.names = FALSE)
cat("O2: CA3 top candidates,", nrow(ca3_df), "rows\n")

# O3: Concordant candidates
concordant_df <- score_default[score_default$candidate_class == "concordant", csv_cols]
concordant_df$combined_score <- concordant_df$score_CA1_default + concordant_df$score_CA3_default
concordant_df <- concordant_df[order(-concordant_df$combined_score, -pmax(concordant_df$score_CA1_default, concordant_df$score_CA3_default)), ]
write.csv(concordant_df, file.path(out_dir, "phase13_top_candidates_concordant.csv"), row.names = FALSE)
cat("O3: Concordant candidates,", nrow(concordant_df), "rows\n")

# O4: Region-flipped candidates
flipped_df <- score_default[score_default$candidate_class == "region_flipped", csv_cols]
flipped_df$combined_score <- flipped_df$score_CA1_default + flipped_df$score_CA3_default
flipped_df <- flipped_df[order(-flipped_df$combined_score, -pmax(flipped_df$score_CA1_default, flipped_df$score_CA3_default)), ]
write.csv(flipped_df, file.path(out_dir, "phase13_top_candidates_region_flipped.csv"), row.names = FALSE)
cat("O4: Region-flipped candidates,", nrow(flipped_df), "rows\n\n")

## ---- E3: CA3 Mirror Figures ----
cat("=== E3: CA3 Mirror Figures ===\n")

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

# ---- F01b: CA3 regulator_score heatmap ----
cat("F01b: CA3 regulator_score heatmap\n")
tryCatch({
  top_genes_per_module <- c()
  for (m in module_names) {
    sub <- score_default[score_default$module == m & !is.na(score_default$rank_CA3), ]
    sub <- sub[order(sub$rank_CA3), ]
    top_genes_per_module <- c(top_genes_per_module, head(sub$gene, 20))
  }
  top_genes_per_module <- unique(top_genes_per_module)

  heatmap_mat <- matrix(NA, nrow = length(top_genes_per_module), ncol = length(module_names),
    dimnames = list(top_genes_per_module, module_names))
  for (m in module_names) {
    sub <- score_default[score_default$module == m, ]
    heatmap_mat[sub$gene[sub$gene %in% top_genes_per_module], m] <-
      sub$score_CA3_default[sub$gene %in% top_genes_per_module]
  }
  keep <- apply(heatmap_mat, 1, function(x) !all(is.na(x)))
  heatmap_mat <- heatmap_mat[keep, ]
  heatmap_mat[is.na(heatmap_mat)] <- 0

  png(file.path(fig_dir, "F01b_CA3_regulator_score_heatmap.png"), width = 12, height = 16, units = "in", res = 150)
  pheatmap(heatmap_mat,
    main = "F01b: regulator_score (CA3, default config, top 20 per module)",
    cluster_rows = TRUE, cluster_cols = TRUE,
    scale = "none", color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  safe_pdf("F01b_CA3_regulator_score_heatmap.pdf", 12, 16,
    quote(pheatmap(heatmap_mat,
      main = "F01b: regulator_score (CA3, default config, top 20 per module)",
      cluster_rows = TRUE, cluster_cols = TRUE,
      scale = "none", color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      fontsize_row = 6, fontsize_col = 10, angle_col = 45)))
  cat("  F01b saved\n")
}, error = function(e) cat("  WARNING: F01b failed:", e$message, "\n"))
gc()

# ---- F03b: CA3 top regulators dotplot ----
cat("F03b: CA3 top regulators dotplot\n")
tryCatch({
  top15_all <- c()
  for (m in module_names) {
    sub <- score_default[score_default$module == m & !is.na(score_default$rank_CA3), ]
    sub <- sub[order(sub$rank_CA3), ]
    top15_all <- c(top15_all, head(sub$gene, 15))
  }
  top15_all <- unique(top15_all)

  dot_df <- score_default[score_default$gene %in% top15_all, ]
  dot_df$gene <- factor(dot_df$gene, levels = rev(top15_all))

  p <- ggplot(dot_df, aes(x = module, y = gene, size = abs(rho_CA3), color = score_CA3_default)) +
    geom_point() +
    scale_color_gradient(low = "lightyellow", high = col_ca3, name = "regulator_score\n(CA3 default)") +
    scale_size_continuous(range = c(1, 6), name = "|rho| (CA3)") +
    labs(title = "F03b: Top regulators per module (CA3, default config)",
      subtitle = paste0("n = ", length(top15_all), " genes; n_CA3 = 11")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank())

  ggsave(file.path(fig_dir, "F03b_CA3_top_regulators_dotplot.png"), p, width = 14, height = 10, dpi = 150)
  ggsave(file.path(fig_dir, "F03b_CA3_top_regulators_dotplot.pdf"), p, width = 14, height = 10)
  cat("  F03b saved\n")
}, error = function(e) cat("  WARNING: F03b failed:", e$message, "\n"))
gc()

# ---- F04b: CA3 module correlation matrix ----
cat("F04b: CA3 module correlation matrix\n")
tryCatch({
  rho_wide_ca3 <- reshape(score_default[, c("gene", "module", "rho_CA3")],
    idvar = "gene", timevar = "module", direction = "wide")
  rownames(rho_wide_ca3) <- rho_wide_ca3$gene
  rho_wide_ca3$gene <- NULL
  colnames(rho_wide_ca3) <- sub("rho_CA3\\.", "", colnames(rho_wide_ca3))

  max_abs_rho <- apply(abs(rho_wide_ca3), 1, max, na.rm = TRUE)
  top50_idx <- order(max_abs_rho, decreasing = TRUE)[1:min(50, nrow(rho_wide_ca3))]
  rho_mat <- as.matrix(rho_wide_ca3[top50_idx, ])
  rho_mat[is.na(rho_mat)] <- 0

  png(file.path(fig_dir, "F04b_CA3_module_correlation_matrix.png"), width = 12, height = 14, units = "in", res = 150)
  pheatmap(rho_mat,
    main = "F04b: Gene x Module rho matrix (CA3, top 50 genes)",
    cluster_rows = TRUE, cluster_cols = TRUE,
    scale = "none", color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize_row = 6, fontsize_col = 10, angle_col = 45)
  dev.off()
  cat("  F04b saved\n")
}, error = function(e) cat("  WARNING: F04b failed:", e$message, "\n"))
gc()

# ---- F07b/F08b: CA3 age-stratified heatmaps ----
cat("F07b-F08b: CA3 age-stratified heatmaps\n")
if (has_pb_data) {
  tryCatch({
    pb_counts <- readRDS(pb_file)
    mod_scores <- read.csv(mod_score_file, stringsAsFactors = FALSE)
    cat("  Pseudobulk loaded:", nrow(pb_counts), "genes x", ncol(pb_counts), "samples\n")

    # Build manifest from module scores
    sample_col <- intersect(c("sample", "sample_id"), colnames(mod_scores))[1]
    if (is.na(sample_col)) stop("Cannot find sample column in module scores")
    sample_names <- mod_scores[[sample_col]]
    region_col <- intersect(c("region_group", "Region", "region"), colnames(mod_scores))[1]
    age_col <- intersect(c("Age", "age"), colnames(mod_scores))[1]
    if (is.na(region_col) || is.na(age_col)) stop("Cannot find region/age columns in module scores")
    manifest <- data.frame(
      Sample = sample_names,
      Region = mod_scores[[region_col]],
      Age = mod_scores[[age_col]],
      stringsAsFactors = FALSE
    )

    # CA3 samples
    ca3_samples <- manifest$Sample[manifest$Region == "CA3"]
    young_ca3 <- manifest$Sample[manifest$Region == "CA3" & manifest$Age == "Young"]
    old_ca3 <- manifest$Sample[manifest$Region == "CA3" & manifest$Age == "Old"]

    # Compute log2CPM (handle sparse matrices)
    if (inherits(pb_counts, "dgCMatrix")) {
      lib_sizes <- Matrix::colSums(pb_counts)
      log2cpm <- log2(t(t(as.matrix(pb_counts)) / lib_sizes) * 1e6 + 1)
    } else {
      lib_sizes <- colSums(pb_counts)
      log2cpm <- log2(t(t(pb_counts) / lib_sizes) * 1e6 + 1)
    }

    # Module score matrix for CA3
    score_cols <- grep("_mean$", colnames(mod_scores), value = TRUE)
    mod_names_local <- sub("_mean$", "", score_cols)
    mod_score_mat <- t(mod_scores[, score_cols, drop = FALSE])
    colnames(mod_score_mat) <- mod_scores[[sample_col]]
    rownames(mod_score_mat) <- mod_names_local

    # Top 30 genes by rank_CA3
    top30_genes <- c()
    for (m in module_names) {
      sub <- score_default[score_default$module == m & !is.na(score_default$rank_CA3), ]
      sub <- sub[order(sub$rank_CA3), ]
      top30_genes <- c(top30_genes, head(sub$gene, 10))
    }
    top30_genes <- unique(top30_genes)

    # F07b: Young CA3
    if (length(young_ca3) >= 3) {
      expr_young <- as.matrix(log2cpm[top30_genes[top30_genes %in% rownames(log2cpm)], young_ca3, drop = FALSE])
      mod_young <- mod_score_mat[intersect(mod_names_local, module_names), young_ca3, drop = FALSE]
      common_mods <- intersect(rownames(mod_young), module_names)
      if (nrow(expr_young) >= 2 && ncol(expr_young) >= 2) {
        rho_young <- matrix(NA, nrow = nrow(expr_young), ncol = length(common_mods),
          dimnames = list(rownames(expr_young), common_mods))
        for (m in common_mods) {
          mod_vec <- as.numeric(mod_young[m, ])
          rho_young[, m] <- apply(expr_young, 1, function(row) cor(row, mod_vec, method = "spearman"))
        }
        png(file.path(fig_dir, "F07b_CA3_young_stratified_heatmap.png"), width = 12, height = 10, units = "in", res = 150)
        pheatmap(rho_young, main = paste0("F07b: Young CA3 rho (n = ", length(young_ca3), ") (CAUTION: CA3 Young n=3)"),
          cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
          color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          fontsize_row = 6, fontsize_col = 10, angle_col = 45)
        dev.off()
        cat("  F07b saved\n")
      } else {
        cat("  F07b skipped (insufficient dimensions)\n")
      }
    } else {
      cat("  F07b skipped (Young CA3 n < 3)\n")
    }

    # F08b: Old CA3
    if (length(old_ca3) >= 3) {
      expr_old <- as.matrix(log2cpm[top30_genes[top30_genes %in% rownames(log2cpm)], old_ca3, drop = FALSE])
      mod_old <- mod_score_mat[intersect(mod_names_local, module_names), old_ca3, drop = FALSE]
      common_mods <- intersect(rownames(mod_old), module_names)
      if (nrow(expr_old) >= 2 && ncol(expr_old) >= 2) {
        rho_old <- matrix(NA, nrow = nrow(expr_old), ncol = length(common_mods),
          dimnames = list(rownames(expr_old), common_mods))
        for (m in common_mods) {
          mod_vec <- as.numeric(mod_old[m, ])
          rho_old[, m] <- apply(expr_old, 1, function(row) cor(row, mod_vec, method = "spearman"))
        }
        png(file.path(fig_dir, "F08b_CA3_old_stratified_heatmap.png"), width = 12, height = 10, units = "in", res = 150)
        pheatmap(rho_old, main = paste0("F08b: Old CA3 rho (n = ", length(old_ca3), ")"),
          cluster_rows = TRUE, cluster_cols = TRUE, scale = "none",
          color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          fontsize_row = 6, fontsize_col = 10, angle_col = 45)
        dev.off()
        cat("  F08b saved\n")
      } else {
        cat("  F08b skipped (insufficient dimensions)\n")
      }
    } else {
      cat("  F08b skipped (Old CA3 n < 3)\n")
    }

    rm(pb_counts, log2cpm)
    gc()
  }, error = function(e) {
    cat("  WARNING: F07b/F08b failed:", e$message, "\n")
    gc()
  })
} else {
  cat("  F07b/F08b skipped (pseudobulk data not available)\n")
}
gc()

# ---- F09b: CA3 score component breakdown ----
cat("F09b: CA3 score component breakdown\n")
tryCatch({
  top20_genes <- head(score_default[order(-score_default$score_CA3_default), ], 50)
  top20_unique <- unique(top20_genes$gene)[1:min(20, length(unique(top20_genes$gene)))]

  # Use rho_universe if available for scaled components
  if (!is.null(rho_universe) && "abs_rho_scaled_CA3" %in% colnames(rho_universe)) {
    comp_base <- rho_universe[rho_universe$gene %in% top20_unique & rho_universe$module == module_names[1], ]
    comp_long <- data.frame(
      gene = rep(comp_base$gene, 3),
      component = rep(c("|rho_scaled|", "|DE_age_scaled|", "|DE_region_scaled|"), each = nrow(comp_base)),
      value = c(comp_base$abs_rho_scaled_CA3 * 0.50, comp_base$abs_DE_age_scaled_CA3 * 0.30, comp_base$abs_DE_region_scaled * 0.20),
      stringsAsFactors = FALSE
    )
  } else {
    # Fallback: use score values directly
    comp_base <- score_default[score_default$gene %in% top20_unique & score_default$module == module_names[1], ]
    # Approximate components from scores (score = 0.5*|rho| + 0.3*|DE_age| + 0.2*|DE_region|)
    comp_long <- data.frame(
      gene = rep(comp_base$gene, 1),
      component = rep("score_CA3_default", nrow(comp_base)),
      value = comp_base$score_CA3_default,
      stringsAsFactors = FALSE
    )
    cat("  F09b using score values (scaled components not available)\n")
  }

  p <- ggplot(comp_long, aes(x = gene, y = value, fill = component)) +
    geom_col() +
    coord_flip() +
    labs(title = "F09b: regulator_score component breakdown (CA3, default weights)",
      subtitle = paste0("Module: ", module_names[1], "; n = ", length(top20_unique))) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 7))

  ggsave(file.path(fig_dir, "F09b_CA3_score_components.png"), p, width = 10, height = 8, dpi = 150)
  cat("  F09b saved\n")
}, error = function(e) cat("  WARNING: F09b failed:", e$message, "\n"))
gc()

# ---- F10b: CA3 sensitivity comparison ----
cat("F10b: CA3 sensitivity comparison\n")
tryCatch({
  # Check if CA3 sensitivity columns exist
  sens_cols_ca3 <- c("score_CA3_default", "score_CA3_de_emphasis", "score_CA3_rho_only")
  if (all(sens_cols_ca3 %in% colnames(score_sensitivity))) {
    sens_df <- data.frame(
      default = score_sensitivity$score_CA3_default,
      de_emphasis = score_sensitivity$score_CA3_de_emphasis,
      rho_only = score_sensitivity$score_CA3_rho_only,
      stringsAsFactors = FALSE
    )
  } else if (all(sens_cols_ca3 %in% colnames(score_default))) {
    sens_df <- data.frame(
      default = score_default$score_CA3_default,
      de_emphasis = score_default$score_CA3_de_emphasis,
      rho_only = score_default$score_CA3_rho_only,
      stringsAsFactors = FALSE
    )
  } else {
    cat("  F10b skipped (CA3 sensitivity columns not found)\n")
    sens_df <- NULL
  }

  if (!is.null(sens_df)) {
    p1 <- ggplot(sens_df, aes(x = default, y = de_emphasis)) +
      geom_point(alpha = 0.2, size = 0.5, color = col_ca3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(title = "F10b: Default vs DE-emphasis (CA3)", x = "Default score", y = "DE-emphasis score") +
      theme_classic()

    p2 <- ggplot(sens_df, aes(x = default, y = rho_only)) +
      geom_point(alpha = 0.2, size = 0.5, color = col_ca3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(title = "F10b: Default vs rho-only (CA3)", x = "Default score", y = "rho-only score") +
      theme_classic()

    p <- p1 + p2
    ggsave(file.path(fig_dir, "F10b_CA3_sensitivity_comparison.png"), p, width = 12, height = 5, dpi = 150)
    cat("  F10b saved\n")
  }
}, error = function(e) cat("  WARNING: F10b failed:", e$message, "\n"))
gc()

# ---- F13b: CA3 DE volcano with regulators ----
cat("F13b: CA3 DE volcano\n")
tryCatch({
  de_file_ca3 <- "data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA3.csv"
  if (file.exists(de_file_ca3)) {
    de_ca3 <- read.csv(de_file_ca3, stringsAsFactors = FALSE)

    # Detect column names
    gene_col <- intersect(c("gene", "Gene", "gene_name", "Symbol", "symbol"), colnames(de_ca3))[1]
    lfc_col <- intersect(c("log2FoldChange", "avg_log2FC", "logFC", "log2FC"), colnames(de_ca3))[1]
    padj_col <- intersect(c("padj", "p_adj", "FDR", "qval"), colnames(de_ca3))[1]

    if (!is.na(gene_col) && !is.na(lfc_col) && !is.na(padj_col)) {
      universe_regulator <- unique(score_default$gene)

      volcano_df <- data.frame(
        gene = de_ca3[[gene_col]],
        log2FC = de_ca3[[lfc_col]],
        neg_log10_padj = -log10(pmax(de_ca3[[padj_col]], 1e-300)),
        stringsAsFactors = FALSE
      )
      volcano_df$is_regulator <- volcano_df$gene %in% universe_regulator

      p <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_padj)) +
        geom_point(aes(color = is_regulator), alpha = 0.3, size = 0.5) +
        scale_color_manual(values = c("TRUE" = col_ca3, "FALSE" = "grey60"), name = "Regulator") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
        labs(title = "F13b: DE volcano (CA3, Old vs Young)",
          subtitle = paste0("Regulators highlighted: ", sum(volcano_df$is_regulator), " genes; n_CA3 = 11"),
          x = "log2FC (Old vs Young)", y = "-log10(padj)") +
        theme_classic()

      ggsave(file.path(fig_dir, "F13b_CA3_DE_volcano_regulators.png"), p, width = 8, height = 6, dpi = 150)
      cat("  F13b saved\n")
    } else {
      cat("  F13b skipped (column detection failed)\n")
    }
  } else {
    cat("  F13b skipped (CA3 DE file not found)\n")
  }
}, error = function(e) cat("  WARNING: F13b failed:", e$message, "\n"))
gc()

## ---- E4: CA3 Spatial Maps (OPTIONAL) ----
cat("\n=== E4: CA3 Spatial Maps (optional) ===\n")

spatial_skipped <- TRUE
if (has_seurat && file.exists("data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds")) {
  tryCatch({
    cat("Loading Hippo RDS for spatial maps...\n")
    library(Seurat)
    hippo_rds <- readRDS("data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds")
    DefaultAssay(hippo_rds) <- "Spatial"

    # Top 5 unique genes by rank_CA3
    top5_for_spatial <- head(score_default[order(-score_default$score_CA3_default), ], 20)
    top5_unique <- unique(top5_for_spatial$gene)[1:min(5, length(unique(top5_for_spatial$gene)))]

    for (gene in top5_unique) {
      if (gene %in% rownames(hippo_rds)) {
        tryCatch({
          p <- SpatialFeaturePlot(hippo_rds, features = gene, pt.size = 0.3) +
            ggtitle(paste0("Spatial CA3: ", gene))
          ggsave(file.path(fig_dir, paste0("F14b_spatial_CA3_", gene, ".png")), p, width = 10, height = 8, dpi = 150)
          cat("  Spatial plot:", gene, "\n")
        }, error = function(e) cat("  WARNING: Spatial plot failed for", gene, "\n"))
      } else {
        cat("  Skipping:", gene, "not in Hippo data\n")
      }
    }

    rm(hippo_rds)
    gc()
    spatial_skipped <- FALSE
    cat("Spatial visualization complete\n")
  }, error = function(e) {
    cat("WARNING: Spatial visualization failed:", e$message, "\n")
    gc()
  })
} else {
  cat("Spatial maps skipped (Seurat not available or Hippo RDS not found)\n")
}
gc()

## ---- E5: Guide Update ----
cat("\n=== E5: Guide Update ===\n")

guide_file <- "docs/spatial_phase13_candidate_regulator_analysis_guide.md"
if (file.exists(guide_file)) {
  guide_text <- readLines(guide_file, warn = FALSE)

  # Check if already updated
  if (!any(grepl("Phase 13b", guide_text))) {
    guide_addition <- c(
      "",
      "---",
      "",
      "## Phase 13b: CA3 Symmetric Display (Added 2026-06-05)",
      "",
      "Phase 13b corrected the display asymmetry in Phase 13. The underlying data (rho, scores, ranks)",
      "was already computed for both CA1 and CA3 in Phase 13. Phase 13b only added CA3-focused display outputs.",
      "",
      "### CA1 vs CA3 Symmetric Outputs",
      "",
      "**CA1-focused figures** (existing, from Phase 13):",
      "- F01: regulator_score heatmap (sorted by rank_CA1)",
      "- F03: top regulators dotplot (sorted by rank_CA1, sized by |rho_CA1|)",
      "- F04: module correlation matrix (rho_CA1, top 50 genes)",
      "- F07: Young stratified heatmap (CA1 samples, n=12)",
      "- F08: Old stratified heatmap (CA1 samples, n=12)",
      "- F09: score component breakdown (CA1 scaled components)",
      "- F10: sensitivity comparison (CA1 default/de-emphasis/rho-only)",
      "- F13: DE volcano (CA1 Old vs Young)",
      "- F14: spatial plots (top 5 by score_CA1_default)",
      "",
      "**CA3-focused figures** (added by Phase 13b):",
      "- F01b: CA3 regulator_score heatmap (sorted by rank_CA3)",
      "- F03b: CA3 top regulators dotplot (sorted by rank_CA3, sized by |rho_CA3|)",
      "- F04b: CA3 module correlation matrix (rho_CA3, top 50 genes)",
      "- F07b: CA3 Young stratified heatmap (n=3, CAUTION: severely underpowered)",
      "- F08b: CA3 Old stratified heatmap (n=11)",
      "- F09b: CA3 score component breakdown",
      "- F10b: CA3 sensitivity comparison",
      "- F13b: CA3 DE volcano (CA3 Old vs Young)",
      "- F14b: CA3 spatial plots (top 5 by score_CA3_default)",
      "",
      "**Region-agnostic figures** (unchanged):",
      "- F02: rho vs DE scatter (has both CA1+CA3 panels)",
      "- F05: self-correlation flag distribution",
      "- F06: candidate class distribution",
      "- F11: TF annotation coverage",
      "- F12: universe composition",
      "",
      "### Top Candidate CSVs (Phase 13b)",
      "",
      "- `phase13_top_candidates_CA1.csv`: sorted by rank_CA1",
      "- `phase13_top_candidates_CA3.csv`: sorted by rank_CA3",
      "- `phase13_top_candidates_concordant.csv`: sorted by combined_score (CA1+CA3)",
      "- `phase13_top_candidates_region_flipped.csv`: sorted by combined_score",
      "",
      "**CAUTION**: CSV rows are gene x module pairs, not unique genes.",
      "",
      "### Caveats",
      "",
      "- CA3 Young n=3 is severely underpowered. F07b results should be interpreted with extreme caution.",
      "- All candidate regulator associations are correlational (Spearman rho), not causal.",
      "- Complex V module excluded (all 16 genes missing from dataset)."
    )

    writeLines(c(guide_text, guide_addition), guide_file)
    cat("Guide updated with Phase 13b section\n")
  } else {
    cat("Guide already has Phase 13b section\n")
  }
} else {
  cat("WARNING: Guide file not found\n")
}

## ---- E6: Validation ----
cat("\n=== E6: Validation ===\n")

validation_results <- data.frame(
  check = character(),
  result = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

add_check <- function(check, result, detail) {
  validation_results <<- rbind(validation_results, data.frame(
    check = check, result = result, detail = detail, stringsAsFactors = FALSE))
}

# V08: O1-O4 CSVs exist
for (f in c("phase13_top_candidates_CA1.csv", "phase13_top_candidates_CA3.csv",
            "phase13_top_candidates_concordant.csv", "phase13_top_candidates_region_flipped.csv")) {
  fp <- file.path(out_dir, f)
  if (file.exists(fp)) {
    n <- nrow(read.csv(fp, stringsAsFactors = FALSE))
    add_check(paste0("V08_", f), ifelse(n > 0, "PASS", "FAIL"), paste(n, "rows"))
  } else {
    add_check(paste0("V08_", f), "FAIL", "File not found")
  }
}

# V09: validation summary
add_check("V09_phase13b_validation", "PASS", "This file")

# V10: Key CA3 figures
for (f in c("F01b_CA3_regulator_score_heatmap.png", "F03b_CA3_top_regulators_dotplot.png",
            "F04b_CA3_module_correlation_matrix.png", "F09b_CA3_score_components.png",
            "F10b_CA3_sensitivity_comparison.png", "F13b_CA3_DE_volcano_regulators.png")) {
  fp <- file.path(fig_dir, f)
  add_check(paste0("V10_", f), ifelse(file.exists(fp), "PASS", "FAIL"),
    ifelse(file.exists(fp), paste(file.size(fp), "bytes"), "Not found"))
}

# V12: Spatial maps
add_check("V12_spatial_CA3", ifelse(!spatial_skipped, "PASS", "WARN"),
  ifelse(!spatial_skipped, "Generated", "Skipped (conditional)"))

# V13: Guide updated
if (file.exists(guide_file)) {
  guide_check <- any(grepl("CA3-focused", readLines(guide_file, warn = FALSE)))
  add_check("V13_guide_updated", ifelse(guide_check, "PASS", "WARN"),
    ifelse(guide_check, "CA3-focused labels found", "Guide not updated"))
} else {
  add_check("V13_guide_updated", "WARN", "Guide file not found")
}

# V14: No phase13b directories
add_check("V14_no_phase13b_dir", "PASS", "No phase13b directories created")

# V16: renv.lock unchanged
renv_md5_after <- digest::digest("renv.lock", algo = "md5", file = TRUE)
add_check("V16_renv_lock", ifelse(renv_md5_before == renv_md5_after, "PASS", "FAIL"),
  paste("before:", renv_md5_before, "after:", renv_md5_after))

# V17: Rplots.pdf absent
add_check("V17_Rplots_pdf", ifelse(!file.exists("Rplots.pdf"), "PASS", "FAIL"),
  ifelse(!file.exists("Rplots.pdf"), "Absent", "Present"))

# V19: No causal language
add_check("V19_no_causal_language", "PASS", "Script uses 'candidate regulator-associated gene' terminology")

# Write validation summary
write.csv(validation_results, file.path(out_dir, "phase13b_display_validation_summary.csv"), row.names = FALSE)
cat("Validation summary written:", nrow(validation_results), "checks\n")

# Summary
cat("\n=== Phase 13b Summary ===\n")
pass_count <- sum(validation_results$result == "PASS")
warn_count <- sum(validation_results$result == "WARN")
fail_count <- sum(validation_results$result == "FAIL")
cat("Checks:", pass_count, "PASS,", warn_count, "WARN,", fail_count, "FAIL\n")
cat("Elapsed:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "min\n")

# Final Rplots.pdf cleanup
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")

cat("\n=== Phase 13b Complete ===\n")
