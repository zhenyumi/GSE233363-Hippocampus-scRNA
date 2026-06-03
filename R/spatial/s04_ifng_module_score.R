#!/usr/bin/env Rscript
# Phase Spatial-04: IFNγ Module Score + IS/Edge Labeling
# Strict reproduction of Script 8 lines 1-66
# Author: ferroDG paper Script 8 (lines 1-66 only)
# Scope: IFNγ module score, IS label, FuncRegion, Phase 05 handoff CSV
# Excludes: STutility, gradient DESeq2, Monocle, Tangram, GO/KEGG

cat("=== Phase Spatial-04: IFNγ Module Score + IS/Edge Labeling ===\n")
cat("Start:", format(Sys.time()), "\n\n")

library(Seurat)
library(ggplot2)

# ── Paths ─────────────────────────────────────────────────────────────────────
input_path  <- "data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds"
out_dir     <- "data/processed/spatial/phase04_ifng_module"
fig_dir     <- "figures/spatial/phase04"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Step 1: Load Hippo object ─────────────────────────────────────────────────
cat("[Step 1] Loading Hippo object...\n")
if (!file.exists(input_path)) stop("STOP: Hippo object not found at ", input_path)
Hippo_All <- readRDS(input_path)
cat("  Loaded:", ncol(Hippo_All), "spots,", nrow(Hippo_All), "genes\n")

# ── Step 2: Validate Seurat object ────────────────────────────────────────────
cat("[Step 2] Validating Seurat object...\n")
stopifnot("SCT" %in% Assays(Hippo_All))
cat("  DefaultAssay:", DefaultAssay(Hippo_All), "\n")

required_meta <- c("Age", "Region", "GEMgroup")
missing_meta  <- setdiff(required_meta, colnames(Hippo_All@meta.data))
if (length(missing_meta) > 0) stop("STOP: missing metadata columns: ", paste(missing_meta, collapse=", "))
cat("  Metadata columns OK: Age, Region, GEMgroup\n")

# Record spot counts by Age
age_table <- table(Hippo_All$Age)
cat("  Age table:", paste(names(age_table), age_table, sep="=", collapse=", "), "\n")

# ── Step 3: Retrieve IFNγ gene set from msigdbr ───────────────────────────────
cat("[Step 3] Retrieving HALLMARK_INTERFERON_GAMMA_RESPONSE from msigdbr...\n")
library(msigdbr)
library(dplyr)

msigdbr_ver <- as.character(packageVersion("msigdbr"))
cat("  msigdbr version:", msigdbr_ver, "\n")

# msigdbr v26 uses 'collection' instead of 'category'
gene_set_df <- msigdbr(species = "Mus musculus", collection = "H")
ifng_genes <- gene_set_df %>%
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  pull(gene_symbol) %>%
  unique()
cat("  IFNγ genes retrieved:", length(ifng_genes), "\n")

if (length(ifng_genes) < 10) stop("STOP: fewer than 10 IFNγ genes retrieved")

# Save gene set
write.csv(
  data.frame(gene_symbol = ifng_genes, stringsAsFactors = FALSE),
  file.path(out_dir, "ifng_gene_set_resolved.csv"),
  row.names = FALSE
)
cat("  Saved ifng_gene_set_resolved.csv\n")

# ── Step 4: Check overlap with SCT assay ──────────────────────────────────────
cat("[Step 4] Checking overlap with SCT assay...\n")
sct_genes     <- rownames(Hippo_All)
overlap_genes <- intersect(ifng_genes, sct_genes)
missing_genes <- setdiff(ifng_genes, sct_genes)
cat("  SCT assay genes:", length(sct_genes), "\n")
cat("  IFNγ genes in SCT:", length(overlap_genes), "/", length(ifng_genes), "\n")

# Save missing genes
write.csv(
  data.frame(gene_symbol = missing_genes, stringsAsFactors = FALSE),
  file.path(out_dir, "ifng_missing_genes.csv"),
  row.names = FALSE
)
cat("  Saved ifng_missing_genes.csv (", length(missing_genes), " missing)\n")

if (length(overlap_genes) < 10) stop("STOP: fewer than 10 IFNγ genes overlap SCT assay")

# ── Step 5: AddModuleScore ─────────────────────────────────────────────────────
cat("[Step 5] Running AddModuleScore (assay=SCT, name=Hallmark_IFNg)...\n")
Hippo_All <- AddModuleScore(
  Hippo_All,
  features = list(ifng_genes),
  name     = "Hallmark_IFNg",
  assay    = "SCT"
)

# Rename Hallmark_IFNg1 → Hallmark_IFNg (AddModuleScore appends "1")
if (!"Hallmark_IFNg1" %in% colnames(Hippo_All@meta.data))
  stop("STOP: Hallmark_IFNg1 not found after AddModuleScore")
Hippo_All$Hallmark_IFNg <- Hippo_All$Hallmark_IFNg1
Hippo_All$Hallmark_IFNg1 <- NULL
cat("  Renamed Hallmark_IFNg1 → Hallmark_IFNg\n")

# ── Step 6: Define IS / Edge labeling ─────────────────────────────────────────
cat("[Step 6] Defining IS (Hallmark_IFNg > 0) and FuncRegion...\n")
is_spots  <- colnames(Hippo_All)[Hippo_All$Hallmark_IFNg > 0]
n_is      <- length(is_spots)
n_total   <- ncol(Hippo_All)
cat("  IS spots:", n_is, "/", n_total, "\n")

# FuncRegion: Region for non-IS, "Edge" for IS
Hippo_All$FuncRegion <- as.character(Hippo_All$Region)
Hippo_All$FuncRegion[is_spots] <- "Edge"
Hippo_All$FuncRegion <- factor(Hippo_All$FuncRegion)
Hippo_All$is_edge    <- Hippo_All$FuncRegion == "Edge"

# ── Step 7: Compute IS proportions by Age ─────────────────────────────────────
cat("[Step 7] Computing IS proportions by Age...\n")
is_prop <- data.frame(
  Age    = c("Young", "Middle", "Old"),
  n_IS   = c(
    sum(Hippo_All$Age == "Young"  & Hippo_All$is_edge),
    sum(Hippo_All$Age == "Middle" & Hippo_All$is_edge),
    sum(Hippo_All$Age == "Old"    & Hippo_All$is_edge)
  ),
  n_total = c(
    sum(Hippo_All$Age == "Young"),
    sum(Hippo_All$Age == "Middle"),
    sum(Hippo_All$Age == "Old")
  ),
  stringsAsFactors = FALSE
)
is_prop$proportion <- is_prop$n_IS / is_prop$n_total

cat("  IS proportions:\n")
for (i in seq_len(nrow(is_prop))) {
  cat(sprintf("    %s: %d/%d = %.2f%%\n",
    is_prop$Age[i], is_prop$n_IS[i], is_prop$n_total[i],
    is_prop$proportion[i] * 100))
}

write.csv(is_prop, file.path(out_dir, "is_proportion_by_age.csv"), row.names = FALSE)
cat("  Saved is_proportion_by_age.csv\n")

# ── Step 8: Tiered IS proportion validation ────────────────────────────────────
cat("[Step 8] Validating IS proportions against author values...\n")
expected <- data.frame(
  Age  = c("Young", "Middle", "Old"),
  Exp  = c(61/2483, 60/2402, 299/5036),
  stringsAsFactors = FALSE
)

validation_results <- data.frame(
  Age = c("Young", "Middle", "Old"),
  observed_n    = is_prop$n_IS,
  observed_total = is_prop$n_total,
  observed_pct  = is_prop$proportion * 100,
  expected_pct  = expected$Exp * 100,
  abs_diff_pp   = abs(is_prop$proportion - expected$Exp) * 100,
  tier          = character(3),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(validation_results))) {
  dd <- validation_results$abs_diff_pp[i]
  validation_results$tier[i] <- if (dd <= 0.25) "PASS" else if (dd <= 1.0) "WARNING" else "FAIL"
}
cat("  Validation:\n")
for (i in seq_len(nrow(validation_results))) {
  cat(sprintf("    %s: obs=%.2f%% exp=%.2f%% diff=%.2fpp [%s]\n",
    validation_results$Age[i],
    validation_results$observed_pct[i],
    validation_results$expected_pct[i],
    validation_results$abs_diff_pp[i],
    validation_results$tier[i]))
}

overall_status <- if (any(validation_results$tier == "FAIL")) "FAIL" else
                  if (any(validation_results$tier == "WARNING")) "WARNING" else "PASS"
cat("  Overall:", overall_status, "\n")

if (overall_status == "FAIL") {
  cat("STOP: IS proportion deviation >1pp detected.\n")
  cat("Diagnostic CSVs saved, but metadata_hippo_func_region.csv is NOT valid for Phase 05.\n")
  cat("Handoff invalid: overall_status = FAIL\n")
}

# ── Step 9: Save module scores CSV ─────────────────────────────────────────────
cat("[Step 9] Saving module scores...\n")
scores_df <- data.frame(
  barcode       = colnames(Hippo_All),
  Hallmark_IFNg = Hippo_All$Hallmark_IFNg,
  Age           = Hippo_All$Age,
  Region        = Hippo_All$Region,
  FuncRegion    = as.character(Hippo_All$FuncRegion),
  is_edge       = Hippo_All$is_edge,
  stringsAsFactors = FALSE
)
write.csv(scores_df, file.path(out_dir, "ifng_module_scores.csv"), row.names = FALSE)
cat("  Saved ifng_module_scores.csv (", nrow(scores_df), " rows)\n")

# ── Step 10: Phase 05 handoff CSV ─────────────────────────────────────────────
cat("[Step 10] Creating Phase 05 handoff CSV...\n")
# Barcode cleaning: strip sample prefix (author Script 8 lines 48-64)
id_cleaned <- colnames(Hippo_All)
for (prefix in paste0(c(10:16, 1:9), "_")) {
  id_cleaned <- gsub(prefix, "", id_cleaned, fixed = TRUE)
}
id_unique <- make.unique(id_cleaned)

handoff <- data.frame(
  barcode       = colnames(Hippo_All),
  id_cleaned    = id_cleaned,
  id_unique     = id_unique,
  Age           = Hippo_All$Age,
  Region        = Hippo_All$Region,
  GEMgroup      = Hippo_All$GEMgroup,
  Hallmark_IFNg = Hippo_All$Hallmark_IFNg,
  FuncRegion    = as.character(Hippo_All$FuncRegion),
  is_edge       = Hippo_All$is_edge,
  stringsAsFactors = FALSE
)
if (overall_status == "FAIL") {
  cat("  SKIP: metadata_hippo_func_region.csv not saved (FAIL blocks handoff)\n")
} else {
  write.csv(handoff, file.path(out_dir, "metadata_hippo_func_region.csv"), row.names = FALSE)
  cat("  Saved metadata_hippo_func_region.csv (", nrow(handoff), " rows)\n")
}

# ── Step 11: Provenance ────────────────────────────────────────────────────────
cat("[Step 11] Saving provenance...\n")
input_info <- file.info(input_path)

# renv.lock before/after
renv_before_md5 <- "321a75eba050c30bfec1614b14f2b2c5"
renv_after_md5  <- tools::md5sum("renv.lock")
renv_pkgs_before <- 169
renv_pkgs_after  <- 172

provenance <- data.frame(
  key   = c(
    "input_path", "input_file_size_bytes", "input_file_mtime",
    "r_version", "seurat_version", "msigdbr_version",
    "msigdbr_species", "msigdbr_collection", "msigdbr_gs_name",
    "ifng_gene_count", "ifng_overlap_sct", "ifng_missing_count",
    "addmodulescore_assay", "addmodulescore_name",
    "is_threshold", "is_spots", "total_spots",
    "young_is_n", "young_total", "young_pct",
    "middle_is_n", "middle_total", "middle_pct",
    "old_is_n", "old_total", "old_pct",
    "validation_status",
    "renv_lock_before_md5", "renv_lock_after_md5",
    "renv_pkgs_before", "renv_pkgs_after",
    "renv_packages_added",
    "checksum_note"
  ),
  value = c(
    input_path,
    as.character(input_info$size),
    as.character(input_info$mtime),
    R.version.string,
    as.character(packageVersion("Seurat")),
    msigdbr_ver,
    "Mus musculus",
    "H (Hallmark)",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    as.character(length(ifng_genes)),
    as.character(length(overlap_genes)),
    as.character(length(missing_genes)),
    "SCT",
    "Hallmark_IFNg",
    "Hallmark_IFNg > 0",
    as.character(n_is),
    as.character(n_total),
    as.character(is_prop$n_IS[is_prop$Age=="Young"]),
    as.character(is_prop$n_total[is_prop$Age=="Young"]),
    sprintf("%.2f", is_prop$proportion[is_prop$Age=="Young"] * 100),
    as.character(is_prop$n_IS[is_prop$Age=="Middle"]),
    as.character(is_prop$n_total[is_prop$Age=="Middle"]),
    sprintf("%.2f", is_prop$proportion[is_prop$Age=="Middle"] * 100),
    as.character(is_prop$n_IS[is_prop$Age=="Old"]),
    as.character(is_prop$n_total[is_prop$Age=="Old"]),
    sprintf("%.2f", is_prop$proportion[is_prop$Age=="Old"] * 100),
    overall_status,
    renv_before_md5,
    as.character(renv_after_md5),
    as.character(renv_pkgs_before),
    as.character(renv_pkgs_after),
    "assertthat 0.2.1, babelgene 22.9, msigdbr 26.1.0",
    "checksum skipped unless explicitly requested"
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "phase04_provenance.csv"), row.names = FALSE)
cat("  Saved phase04_provenance.csv\n")

# ── Step 12: Validation summary ────────────────────────────────────────────────
cat("[Step 12] Writing validation summary...\n")
val_lines <- c(
  "=== Phase Spatial-04 Validation Summary ===",
  paste("Date:", Sys.time()),
  "",
  "--- Input ---",
  paste("Input:", input_path),
  paste("Spots:", n_total),
  paste("Genes (SCT):", length(sct_genes)),
  "",
  "--- msigdbr ---",
  paste("msigdbr version:", msigdbr_ver),
  paste("Species: Mus musculus"),
  paste("Collection: H (Hallmark)"),
  paste("Gene set: HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  paste("IFNγ genes retrieved:", length(ifng_genes)),
  paste("IFNγ genes in SCT:", length(overlap_genes)),
  paste("IFNγ genes missing:", length(missing_genes)),
  "",
  "--- AddModuleScore ---",
  paste("Assay: SCT"),
  paste("Name: Hallmark_IFNg"),
  "",
  "--- IS Proportions ---",
  paste("Threshold: Hallmark_IFNg > 0"),
  paste("IS spots:", n_is, "/", n_total),
  "",
  sprintf("%-10s %6s %6s %8s %8s %8s %s", "Age", "IS_n", "Total", "Obs%", "Exp%", "Diff(pp)", "Tier"),
  sprintf("%-10s %6d %6d %7.2f%% %7.2f%% %7.2f   %s",
    validation_results$Age,
    validation_results$observed_n,
    validation_results$observed_total,
    validation_results$observed_pct,
    validation_results$expected_pct,
    validation_results$abs_diff_pp,
    validation_results$tier),
  "",
  paste("Overall validation:", overall_status),
  "",
  "--- renv.lock ---",
  paste("Before MD5:", renv_before_md5),
  paste("After MD5:", renv_after_md5),
  paste("Packages before:", renv_pkgs_before),
  paste("Packages after:", renv_pkgs_after),
  paste("Added: assertthat 0.2.1, babelgene 22.9, msigdbr 26.1.0"),
  "",
  "--- Outputs ---",
  paste("ifng_gene_set_resolved.csv:", length(ifng_genes), "genes"),
  paste("ifng_missing_genes.csv:", length(missing_genes), "genes"),
  paste("ifng_module_scores.csv:", nrow(scores_df), "rows"),
  paste("is_proportion_by_age.csv:", nrow(is_prop), "rows"),
  paste("metadata_hippo_func_region.csv:", if (overall_status == "FAIL") "NOT SAVED (FAIL)" else paste(nrow(handoff), "rows (Phase 05 handoff)")),
  paste("phase04_provenance.csv: saved"),
  paste("Validation status:"),
  paste("  Young:  expected 2.46% (61/2483)"),
  paste("  Middle: expected 2.50% (60/2402)"),
  paste("  Old:    expected 5.94% (299/5036)"),
  "",
  "--- Phase 05 Handoff Status ---",
  if (overall_status == "PASS") {
    "Phase 04 completed with PASS. metadata_hippo_func_region.csv is available for Phase 05."
  } else if (overall_status == "WARNING") {
    "Phase 04 completed with WARNING. Phase 05 use of metadata_hippo_func_region.csv requires explicit user acknowledgement/approval."
  } else {
    "Phase 04 completed with FAIL. metadata_hippo_func_region.csv is NOT valid for Phase 05. Diagnostic outputs only."
  },
  "",
  if (overall_status != "PASS") {
    "Possible contributors to deviation include missing IFNγ genes in the SCT assay, msigdbr version, Seurat/AddModuleScore version, or SCT feature differences."
  }
)
writeLines(val_lines, file.path(out_dir, "phase04_validation_summary.txt"))
cat("  Saved phase04_validation_summary.txt\n")

# ── Step 13: Figures ────────────────────────────────────────────────────────────
cat("[Step 13] Generating figures...\n")

# Histogram (required — Script 8 lines 21-27)
p_hist <- ggplot(Hippo_All@meta.data, aes(x = Hallmark_IFNg)) +
  geom_histogram(bins = 50) +
  geom_vline(aes(xintercept = 0), color = "red", linetype = 2) +
  theme_bw() +
  labs(title = "IFNγ Module Score Distribution",
       x = "Hallmark_IFNg Module Score",
       y = "Count")

tryCatch({
  ggsave(file.path(fig_dir, "ifng_score_histogram.png"), p_hist,
         width = 4.5, height = 4.5)
  cat("  Saved ifng_score_histogram.png\n")
}, error = function(e) cat("  PNG histogram failed:", conditionMessage(e), "\n"))

tryCatch({
  ggsave(file.path(fig_dir, "ifng_score_histogram.pdf"), p_hist,
         width = 4.5, height = 4.5)
  cat("  Saved ifng_score_histogram.pdf\n")
}, error = function(e) cat("  PDF histogram failed:", conditionMessage(e), "\n"))

# Spatial FeaturePlot (optional best-effort)
tryCatch({
  DefaultAssay(Hippo_All) <- "SCT"
  p_spatial <- SpatialFeaturePlot(Hippo_All, features = "Hallmark_IFNg",
                                   pt.size.factor = 1.5)
  png(file.path(fig_dir, "ifng_spatial_featureplot.png"), width = 900, height = 800)
  print(p_spatial)
  dev.off()
  cat("  Saved ifng_spatial_featureplot.png\n")
}, error = function(e) cat("  Spatial featureplot PNG failed:", conditionMessage(e), "\n"))

tryCatch({
  cairo_pdf(file.path(fig_dir, "ifng_spatial_featureplot.pdf"), width = 9, height = 8)
  print(p_spatial)
  dev.off()
  cat("  Saved ifng_spatial_featureplot.pdf\n")
}, error = function(e) cat("  Spatial featureplot PDF failed:", conditionMessage(e), "\n"))

# ── Step 14: Cleanup ──────────────────────────────────────────────────────────
cat("[Step 14] Cleaning up...\n")
rm(Hippo_All)
gc()
cat("\n=== Phase Spatial-04 COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
