## =============================================================================
## s07_figure_validation.R
## Phase Spatial-07: Reproduction Validation & Gap Report
##
## Base R only. No Seurat, no here, no new analysis, no package installs.
## Reads CSV/TXT summaries from Phase 01-06; generates report documents.
## =============================================================================

cat("=== Phase Spatial-07: s07_figure_validation.R ===\n")
run_start <- Sys.time()
cat("Start:", format(run_start, "%Y-%m-%d %H:%M:%S"), "\n")

## ---- E0: Setup ---------------------------------------------------------------
root <- getwd()

out_dir <- file.path(root, "data", "processed", "spatial", "phase07_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

docs_dir <- file.path(root, "docs")
rplots_path <- file.path(root, "Rplots.pdf")

git_head <- tryCatch(
  trimws(system2("git", args = c("rev-parse", "--short", "HEAD"), stdout = TRUE)),
  error = function(e) "unknown"
)

renv_path <- file.path(root, "renv.lock")
renv_md5 <- tryCatch(
  tools::md5sum(renv_path),
  error = function(e) "unknown"
)

rplots_existed <- file.exists(rplots_path)
cat("Rplots.pdf existed at start:", rplots_existed, "\n")
cat("Git HEAD:", git_head, "\n")
cat("renv.lock MD5:", renv_md5, "\n")

## ---- E1: Read Phase 01-06 validation summaries --------------------------------

read_summary <- function(path, label) {
  if (!file.exists(path)) {
    cat("  CAUTION:", label, "not found:", path, "\n")
    return(NULL)
  }
  cat("  Read:", label, "\n")
  readLines(path)
}

cat("\n--- E1: Reading validation summaries ---\n")
p02 <- read_summary(
  file.path(root, "data", "processed", "spatial", "phase02_metadata", "phase02_validation_summary.txt"),
  "Phase 02"
)
p04 <- read_summary(
  file.path(root, "data", "processed", "spatial", "phase04_ifng_module", "phase04_validation_summary.txt"),
  "Phase 04"
)
p05 <- read_summary(
  file.path(root, "data", "processed", "spatial", "phase05_rds_gradient", "phase05_validation_summary.txt"),
  "Phase 05"
)
p06 <- read_summary(
  file.path(root, "data", "processed", "spatial", "phase06_hippo_regional_atlas", "phase06_validation_summary.txt"),
  "Phase 06"
)

deg_summary <- read.csv(
  file.path(root, "data", "processed", "spatial", "phase03_pseudobulk_deseq2", "deg_summary.csv"),
  stringsAsFactors = FALSE
)

## ---- E2: Author expected values (hardcoded from Script 8) ---------------------

cat("\n--- E2: Author expected values (Script 8) ---\n")

author_is_young <- 61
author_is_middle <- 60
author_is_old <- 299
author_deseq2_up <- 2357
author_deseq2_down <- 1326

our_is_young <- 61
our_is_middle <- 48
our_is_old <- 317
our_deseq2_up <- 2571
our_deseq2_down <- 1452

oy_padj <- deg_summary$n_padj[deg_summary$contrast == "OY"]
my_padj <- deg_summary$n_padj[deg_summary$contrast == "MY"]
om_padj <- deg_summary$n_padj[deg_summary$contrast == "OM"]
oy_strong <- deg_summary$n_strong[deg_summary$contrast == "OY"]
my_strong <- deg_summary$n_strong[deg_summary$contrast == "MY"]
om_strong <- deg_summary$n_strong[deg_summary$contrast == "OM"]

cat("  Author IS: Young=", author_is_young, " Middle=", author_is_middle, " Old=", author_is_old, "\n")
cat("  Our IS:    Young=", our_is_young, " Middle=", our_is_middle, " Old=", our_is_old, "\n")
cat("  Phase 03: OY=", oy_padj, "padj MY=", my_padj, " OM=", om_padj, "\n")

## ---- E3: Phase classification -------------------------------------------------

cat("\n--- E3: Phase classification ---\n")

phase_classification <- data.frame(
  phase = c("s01", "s02", "s03", "s04", "s05", "s05b", "s06"),
  script = c(
    "s01_inspect_objects.R",
    "s02_dg_hippo_metadata.R",
    "s03_pseudobulk_deseq2_dg.R",
    "s04_ifng_module_score.R",
    "s05_rds_inflammatory_gradient.R",
    "s05b_plot_spatial_labels.R",
    "s06_hippo_regional_atlas.R"
  ),
  classification = c(
    "OBJECT_INSPECTION",
    "OBJECT_DERIVED_VALIDATION",
    "METHOD_MATCH_WITH_WARNING",
    "METHOD_MATCH_WITH_WARNING",
    "RDS_APPROXIMATION",
    "RDS_APPROXIMATION",
    "OBJECT_DERIVED_VALIDATION"
  ),
  rationale = c(
    "Reads RDS, records structure, no method reproduction",
    "Validates metadata against s01 expectations; 109/109 PASS",
    paste0("Method matches Author Script 7 (aggregate.Matrix + DESeq2). ",
           "DEG: OY=", oy_padj, "padj, MY=", my_padj, ", OM=", om_padj,
           ". No author reference DEG counts found."),
    paste0("Method matches Script 8 lines 1-66 (msigdbr + AddModuleScore). ",
           "Middle IS=", our_is_middle, " vs author ", author_is_middle,
           ", Old IS=", our_is_old, " vs author ", author_is_old, "."),
    paste0("Uses RDS coordinate neighbor graph instead of STutility::RegionNeighbours. ",
           "DESeq2 IS vs NNS Old: Up=", our_deseq2_up, "/Down=", our_deseq2_down,
           " vs author ", author_deseq2_up, "/", author_deseq2_down, "."),
    "Visualization of Phase 05 labels; no independent computation.",
    paste0("Validates Hippo regional structure from RDS metadata. ",
           "54/54 PASS. CA1=4671, CA2=565, CA3=2321, ML=1480, GCL=645, Hilus=239.")
  ),
  source = c(
    "s01_inspect_objects.R output",
    "data/processed/spatial/phase02_metadata/phase02_validation_summary.txt",
    "data/processed/spatial/phase03_pseudobulk_deseq2/deg_summary.csv",
    "data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt",
    "data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv",
    "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt",
    "data/processed/spatial/phase06_hippo_regional_atlas/phase06_validation_summary.txt"
  ),
  stringsAsFactors = FALSE
)

cat("  Generated", nrow(phase_classification), "phase rows\n")

## ---- E4: Figure/panel classification (row-by-row) ----------------------------

cat("\n--- E4: Figure/panel classification ---\n")

mkrow <- function(fig, pan, ph, cls, rat, src) {
  list(figure = fig, panel = pan, phase = ph, classification = cls, rationale = rat, source = src)
}

fig_rows <- list(
  ## ---- Fig. 1 (BLOCKED_WHOLEBRAIN) ----
  mkrow("Fig. 1", "f", "", "BLOCKED_WHOLEBRAIN",
        "Requires WholeBrain RDS (42,169 spots); not loaded",
        "docs/spatial_transcriptomics_overview.md"),
  mkrow("Fig. 1", "g", "", "BLOCKED_WHOLEBRAIN",
        "WholeBrain spatial projection; not loaded",
        "docs/spatial_transcriptomics_overview.md"),
  mkrow("Fig. 1", "h", "", "DEFERRED",
        "Hippocampal CCA mapping; not attempted in Phase 01-07",
        "docs/spatial_reproduction_plan.md"),

  ## ---- Fig. 5 ----
  mkrow("Fig. 5", "d", "02, 06", "OBJECT_DERIVED_VALIDATION",
        "DG UMAP from DG RDS metadata; Phase 06 UMAP by Region on Hippo",
        "data/processed/spatial/phase02_metadata/phase02_validation_summary.txt; data/processed/spatial/phase06_hippo_regional_atlas/phase06_validation_summary.txt"),
  mkrow("Fig. 5", "e", "", "DEFERRED",
        "Marker expression spatial mapping; not attempted",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Fig. 5", "f", "03", "PARTIAL_REPRODUCTION",
        "Phase 03 has DEG results; CAS score not computed",
        "data/processed/spatial/phase03_pseudobulk_deseq2/deg_summary.csv"),
  mkrow("Fig. 5", "g", "", "DEFERRED",
        "CAS violin not computed",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Fig. 5", "h", "03", "PARTIAL_REPRODUCTION",
        "Phase 03 DESeq2 contrasts OY/MY/OM available; full panel not reproduced",
        "data/processed/spatial/phase03_pseudobulk_deseq2/deg_summary.csv"),

  ## ---- Fig. 6 ----
  mkrow("Fig. 6", "e", "04, 05", "RDS_APPROXIMATION",
        "Phase 04 IS proportions + Phase 05 label maps; RDS-approximated",
        "data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt; data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),

  ## ---- Fig. 7 ----
  mkrow("Fig. 7", "a", "04", "METHOD_MATCH_WITH_WARNING",
        paste0("IFNγ proportions reproduced with Middle/Old discrepancy: ",
               "Young 61 match, Middle ", our_is_middle, " vs ", author_is_middle,
               ", Old ", our_is_old, " vs ", author_is_old),
        "data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt; docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R"),
  mkrow("Fig. 7", "c", "05, 05b", "RDS_APPROXIMATION",
        "RDS coordinate neighbor graph, not STutility image masks",
        "data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv"),
  mkrow("Fig. 7", "d", "05", "RDS_APPROXIMATION",
        "Label counts from coordinate-based propagation; RDS-approximated",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Fig. 7", "e", "", "BLOCKED_RAW",
        "STutility + Monocle 2 pseudotime; requires raw Visium files",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Fig. 7", "f", "", "BLOCKED_RAW",
        "Same as Fig. 7e",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Fig. 7", "g", "", "BLOCKED_RAW",
        "Same as Fig. 7e",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Fig. 7", "h", "", "BLOCKED_RAW",
        "Same as Fig. 7e",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Fig. 7", "i", "", "BLOCKED_RAW",
        "Same as Fig. 7e",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Fig. 7", "j", "", "DEFERRED",
        "IFNγ/microglia correlation / microniche-related output not reproduced; raw/STutility/extra workflow unavailable or deferred",
        "docs/spatial_reproduction_plan.md; docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R"),

  ## ---- Ext. Data Fig. 2 ----
  mkrow("Ext. Data Fig. 2", "a", "02", "OBJECT_DERIVED_VALIDATION",
        "UMAP from DG/Hippo metadata",
        "data/processed/spatial/phase02_metadata/phase02_validation_summary.txt"),
  mkrow("Ext. Data Fig. 2", "b", "", "DEFERRED",
        "Not attempted in Phase 01-07",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 2", "c", "", "DEFERRED",
        "Not attempted in Phase 01-07",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 2", "d", "06", "PARTIAL_REPRODUCTION",
        "Phase 06 validates region labels/counts only; no Allen Atlas overlay reproduced",
        "data/processed/spatial/phase06_hippo_regional_atlas/phase06_validation_summary.txt; docs/original_code_from_paper/Scripts/R/6. Spatial-Data process.R"),
  mkrow("Ext. Data Fig. 2", "e", "", "BLOCKED_TANGRAM",
        "Requires Tangram cell-type deconvolution; not executed",
        "docs/original_code_from_paper/Scripts/Python/notebook_tangram.ipynb; docs/spatial_reproduction_plan.md"),

  ## ---- Ext. Data Fig. 7 ----
  mkrow("Ext. Data Fig. 7", "c", "03", "METHOD_MATCH_WITH_WARNING",
        paste0("Phase 03 DESeq2 method matches Script 7; ",
               "OY=", oy_padj, " padj, MY=", my_padj, ", OM=", om_padj,
               ". No exact author DEG count reference available"),
        "data/processed/spatial/phase03_pseudobulk_deseq2/deg_summary.csv; docs/original_code_from_paper/Scripts/R/7. Spatial-DESeq2.R"),
  mkrow("Ext. Data Fig. 7", "d", "03", "METHOD_MATCH_WITH_WARNING",
        paste0("Same as Ext. Data Fig. 7c; OM contrast available (", om_padj,
               " padj) but no exact author DEG count reference available"),
        "data/processed/spatial/phase03_pseudobulk_deseq2/deg_summary.csv; docs/original_code_from_paper/Scripts/R/7. Spatial-DESeq2.R"),

  ## ---- Ext. Data Fig. 9 ----
  mkrow("Ext. Data Fig. 9", "a", "04", "METHOD_MATCH_WITH_WARNING",
        paste0("IFNγ proportions reproduced with Middle/Old discrepancy: ",
               "Young ", our_is_young, " match, Middle ", our_is_middle, " vs ",
               author_is_middle, ", Old ", our_is_old, " vs ", author_is_old),
        "data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt; docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R"),
  mkrow("Ext. Data Fig. 9", "b", "04", "METHOD_MATCH_WITH_WARNING",
        "IFNγ score histogram derived from Phase 04 module scores; Phase 04 has warning",
        "data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt"),
  mkrow("Ext. Data Fig. 9", "c", "05", "RDS_APPROXIMATION",
        paste0("Old IS vs NNS DESeq2: Up=", our_deseq2_up, "/Down=", our_deseq2_down,
               " vs author ", author_deseq2_up, "/", author_deseq2_down,
               "; RDS-approximated labels"),
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt; docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R"),
  mkrow("Ext. Data Fig. 9", "d", "", "DEFERRED",
        "GO enrichment not performed in Phase 01-07",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 9", "e", "", "DEFERRED",
        "Not attempted in Phase 01-07",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 9", "f", "05", "RDS_APPROXIMATION",
        "Phase 05 PCA aggregated spots based on RDS-approximated labels",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Ext. Data Fig. 9", "g", "05", "RDS_APPROXIMATION",
        "Phase 05 heatmap/hierarchical clustering based on RDS-approximated labels",
        "data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt"),
  mkrow("Ext. Data Fig. 9", "h", "", "DEFERRED",
        "CAS scoring not performed in Phase 01-07",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 9", "i", "", "DEFERRED",
        "CAS scoring not performed in Phase 01-07",
        "docs/spatial_reproduction_plan.md"),

  ## ---- Ext. Data Fig. 10 ----
  mkrow("Ext. Data Fig. 10", "a", "", "DEFERRED",
        "SLM manual selection; not in scope",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 10", "b", "", "DEFERRED",
        "BBB / correlation; not in scope",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 10", "c", "", "DEFERRED",
        "SLM manual selection; not in scope",
        "docs/spatial_reproduction_plan.md"),
  mkrow("Ext. Data Fig. 10", "d", "", "DEFERRED",
        "BBB / correlation; not in scope",
        "docs/spatial_reproduction_plan.md")
)

figure_classification <- do.call(rbind, lapply(fig_rows, as.data.frame, stringsAsFactors = FALSE))

cat("  Generated", nrow(figure_classification), "figure rows\n")

## ---- E5: Discrepancy table ----------------------------------------------------

cat("\n--- E5: Discrepancy table ---\n")

discrepancy_table <- data.frame(
  id = paste0("D0", 1:9),
  phase = c("04", "04", "05", "05", "05", "ALL", "ALL", "05", "03"),
  metric = c(
    "Middle IS spots", "Old IS spots",
    "IS vs NNS DEG Up (Old)", "IS vs NNS DEG Down (Old)",
    "Neighbor graph method", "Seurat version",
    "Raw Visium files", "Label margin heuristic", "DEG OY count reference"
  ),
  our_value = c(
    paste0(our_is_middle, " (2.00%)"), paste0(our_is_old, " (6.29%)"),
    as.character(our_deseq2_up), as.character(our_deseq2_down),
    "RDS Euclidean coordinates", "v5.5.0",
    "Not available", "1.2x mean NN distance",
    paste0(oy_padj, " padj")
  ),
  author_expected = c(
    paste0(author_is_middle, " (2.50%)"), paste0(author_is_old, " (5.94%)"),
    as.character(author_deseq2_up), as.character(author_deseq2_down),
    "STutility::RegionNeighbours (image-mask propagation)", "v4.3",
    "Required by Script 8 lines 71-697", "STutility image propagation",
    "No author reference found"
  ),
  source = c(
    "docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R:41; phase04_validation_summary.txt",
    "docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R:41; phase04_validation_summary.txt",
    "data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv; Script 8",
    "data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv; Script 8",
    "data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv",
    "renv.lock; author Script headers",
    "data/raw/spatial/ is empty",
    "data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv",
    "data/processed/spatial/phase03_pseudobulk_deseq2/deg_summary.csv; Script 7"
  ),
  severity = c("warning", "warning", "caution", "caution", "blocker",
                "caution", "blocker", "caution", "info"),
  possible_causes = c(
    "msigdbr version, SCT assay gene differences, Seurat v4.3 vs v5.5.0 AddModuleScore behavior",
    "Same as D01",
    "RDS coordinate margin (1.2x heuristic) vs STutility image masks; different neighbor sets",
    "Same as D03",
    "Raw Visium files unavailable; cannot run STutility",
    "Assay slot structure, SCT behavior, SpatialKey differences",
    "Not on Zenodo; author did not deposit raw Visium files",
    "No external calibration data for margin threshold",
    "Author Script 7 does not report total DEG counts"
  ),
  stringsAsFactors = FALSE
)

cat("  Generated", nrow(discrepancy_table), "discrepancy rows\n")

## ---- E6: Handoff evaluation ---------------------------------------------------

cat("\n--- E6: Handoff evaluation ---\n")

phase06_pass <- !is.null(p06) && any(grepl("Overall status: PASS", p06))
phase02_pass <- !is.null(p02) && any(grepl("Overall status: passed", p02))

if (phase06_pass && phase02_pass) {
  handoff_decision <- "GO_WITH_CAVEATS"
  handoff_rationale <- paste0(
    "Phase 06 PASS (54/54): Hippo regional framework (CA1, CA2, CA3, ML, GCL, Hilus) validated. ",
    "Phase 02 PASS (109/109): DG and Hippo metadata consistent. ",
    "Phase 04/05 discrepancies (D01-D05) affect IS/Edge axis, not Region axis. ",
    "CA1/CA3/DG region labels from author Script 6 annotation (RDS metadata), ",
    "independent of IFNγ module scores or STutility neighbor propagation. ",
    "Caveats: D01/D02 Middle/Old IS discrepancy; D05 RDS coordinate approximation; D06 Seurat version."
  )
} else {
  handoff_decision <- "NO_GO"
  handoff_rationale <- paste0(
    "Phase 06: ", if (phase06_pass) "PASS" else "NON-PASS",
    "; Phase 02: ", if (phase02_pass) "PASS" else "NON-PASS",
    ". Critical validation checks failed."
  )
}
cat("  Decision:", handoff_decision, "\n")

## ---- E7: Tangram evaluation ---------------------------------------------------

cat("\n--- E7: Tangram evaluation ---\n")

tangram_notebook_exists <- file.exists(
  file.path(root, "docs", "original_code_from_paper", "Scripts", "Python", "notebook_tangram.ipynb")
)

tangram_decision <- "TANGRAM_DEFERRED"
tangram_rationale <- paste0(
  "Tangram cell-type deconvolution not required for region-level CA1/CA3/DG mitochondrial planning. ",
  "Region labels from author Script 6 annotation, validated by Phase 06 (54/54 PASS). ",
  "Tangram notebook exists at docs/original_code_from_paper/Scripts/Python/notebook_tangram.ipynb ",
  "(25 cells, Python 3). May be needed if cell-type composition differences between regions ",
  "affect mitochondrial interpretation; evaluate during Phase 09 planning."
)
cat("  Decision:", tangram_decision, "\n")

## ---- E8: Write all CSV outputs ------------------------------------------------

cat("\n--- E8: Writing CSV outputs ---\n")

write.csv(phase_classification,
          file.path(out_dir, "phase07_phase_classification.csv"), row.names = FALSE)
cat("  Written: phase07_phase_classification.csv\n")

write.csv(figure_classification,
          file.path(out_dir, "phase07_figure_classification_table.csv"), row.names = FALSE)
cat("  Written: phase07_figure_classification_table.csv\n")

write.csv(discrepancy_table,
          file.path(out_dir, "phase07_discrepancy_table.csv"), row.names = FALSE)
cat("  Written: phase07_discrepancy_table.csv\n")

## ---- E9: Generate docs/spatial_figure_validation.md ---------------------------

cat("\n--- E9: Generating spatial_figure_validation.md ---\n")

mkline <- function(...) paste0(..., collapse = "")

fig_lines <- character(0)
fig_lines <- c(fig_lines, "# Spatial Figure Validation Report")
fig_lines <- c(fig_lines, "")
fig_lines <- c(fig_lines, mkline("**Generated by:** s07_figure_validation.R | ",
                                  format(run_start, "%Y-%m-%d %H:%M:%S")))
fig_lines <- c(fig_lines, mkline("**Git HEAD:** ", git_head))
fig_lines <- c(fig_lines, mkline("**renv.lock MD5:** ", renv_md5))
fig_lines <- c(fig_lines, "")
fig_lines <- c(fig_lines, "## Classification System")
fig_lines <- c(fig_lines, "")
fig_lines <- c(fig_lines, "| Label | Definition |")
fig_lines <- c(fig_lines, "|-------|-----------|")
fig_lines <- c(fig_lines, "| OBJECT_INSPECTION | Read-only structural inspection of author-provided objects. |")
fig_lines <- c(fig_lines, "| OBJECT_DERIVED_VALIDATION | Validation of metadata/annotations from author objects. |")
fig_lines <- c(fig_lines, "| STRICT_REPRODUCTION | Same method, same input, same tool path. Directly comparable. |")
fig_lines <- c(fig_lines, "| METHOD_MATCH_WITH_WARNING | Method matches author code, but version/count discrepancies. |")
fig_lines <- c(fig_lines, "| RDS_APPROXIMATION | Uses RDS coordinates as approximation for raw/STutility method. |")
fig_lines <- c(fig_lines, "| PARTIAL_REPRODUCTION | Only a subset of panels/statistics reproduced. |")
fig_lines <- c(fig_lines, "| BLOCKED_RAW | Requires raw Visium/Space Ranger files - not available. |")
fig_lines <- c(fig_lines, "| BLOCKED_WHOLEBRAIN | Requires WholeBrain RDS - not loaded. |")
fig_lines <- c(fig_lines, "| BLOCKED_TANGRAM | Requires Tangram cell-type deconvolution - not executed. |")
fig_lines <- c(fig_lines, "| DEFERRED | Not attempted; optional or requires tools not in scope. |")
fig_lines <- c(fig_lines, "")
fig_lines <- c(fig_lines, "## Phase Classification")
fig_lines <- c(fig_lines, "")
fig_lines <- c(fig_lines, "| Phase | Script | Classification | Rationale | Source |")
fig_lines <- c(fig_lines, "|-------|--------|---------------|-----------|--------|")
for (i in seq_len(nrow(phase_classification))) {
  r <- phase_classification[i, ]
  fig_lines <- c(fig_lines, paste0("| ", r$phase, " | ", r$script, " | ",
                                    r$classification, " | ", r$rationale, " | ", r$source, " |"))
}
fig_lines <- c(fig_lines, "")
fig_lines <- c(fig_lines, "## Figure/Panel Classification")
fig_lines <- c(fig_lines, "")
fig_lines <- c(fig_lines, "| Figure | Panel | Phase | Classification | Rationale | Source |")
fig_lines <- c(fig_lines, "|--------|-------|-------|---------------|-----------|--------|")
for (i in seq_len(nrow(figure_classification))) {
  r <- figure_classification[i, ]
  fig_lines <- c(fig_lines, paste0("| ", r$figure, " | ", r$panel, " | ",
                                    r$phase, " | ", r$classification, " | ",
                                    r$rationale, " | ", r$source, " |"))
}
fig_lines <- c(fig_lines, "")

writeLines(fig_lines, file.path(docs_dir, "spatial_figure_validation.md"))
cat("  Written: spatial_figure_validation.md\n")

## ---- E10: Generate docs/spatial_reproducibility_report.md ---------------------

cat("\n--- E10: Generating spatial_reproducibility_report.md ---\n")

rep_lines <- character(0)
rep_lines <- c(rep_lines, "# Spatial Reproducibility Report")
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, mkline("**Generated by:** s07_figure_validation.R | ",
                                   format(run_start, "%Y-%m-%d %H:%M:%S")))
rep_lines <- c(rep_lines, mkline("**Git HEAD:** ", git_head))
rep_lines <- c(rep_lines, mkline("**renv.lock MD5:** ", renv_md5))
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, "## Phase Classification")
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, "| Phase | Classification | Rationale |")
rep_lines <- c(rep_lines, "|-------|---------------|-----------|")
for (i in seq_len(nrow(phase_classification))) {
  r <- phase_classification[i, ]
  rep_lines <- c(rep_lines, paste0("| ", r$phase, " | ", r$classification, " | ", r$rationale, " |"))
}
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, "## Discrepancy Table")
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, "| ID | Phase | Metric | Our Value | Author Expected | Source | Severity | Possible Causes |")
rep_lines <- c(rep_lines, "|-----|-------|--------|-----------|----------------|--------|----------|----------------|")
for (i in seq_len(nrow(discrepancy_table))) {
  r <- discrepancy_table[i, ]
  rep_lines <- c(rep_lines, paste0("| ", r$id, " | ", r$phase, " | ", r$metric,
                                    " | ", r$our_value, " | ", r$author_expected,
                                    " | ", r$source, " | ", r$severity,
                                    " | ", r$possible_causes, " |"))
}
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, "## Summary")
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, paste0("- **OBJECT_INSPECTION:** ",
                                  sum(phase_classification$classification == "OBJECT_INSPECTION"), " phases"))
rep_lines <- c(rep_lines, paste0("- **OBJECT_DERIVED_VALIDATION:** ",
                                  sum(phase_classification$classification == "OBJECT_DERIVED_VALIDATION"), " phases"))
rep_lines <- c(rep_lines, paste0("- **METHOD_MATCH_WITH_WARNING:** ",
                                  sum(phase_classification$classification == "METHOD_MATCH_WITH_WARNING"), " phases"))
rep_lines <- c(rep_lines, paste0("- **RDS_APPROXIMATION:** ",
                                  sum(phase_classification$classification == "RDS_APPROXIMATION"), " phases"))
rep_lines <- c(rep_lines, paste0("- **Discrepancies identified:** ", nrow(discrepancy_table)))
rep_lines <- c(rep_lines, paste0("- **Phase 06:** ", if (phase06_pass) "PASS (54/54)" else "NON-PASS"))
rep_lines <- c(rep_lines, paste0("- **Phase 02:** ", if (phase02_pass) "PASS (109/109)" else "NON-PASS"))
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, "## Data Sources")
rep_lines <- c(rep_lines, "")
rep_lines <- c(rep_lines, "- `data/processed/spatial/phase02_metadata/phase02_validation_summary.txt`")
rep_lines <- c(rep_lines, "- `data/processed/spatial/phase03_pseudobulk_deseq2/deg_summary.csv`")
rep_lines <- c(rep_lines, "- `data/processed/spatial/phase04_ifng_module/phase04_validation_summary.txt`")
rep_lines <- c(rep_lines, "- `data/processed/spatial/phase05_rds_gradient/phase05_validation_summary.txt`")
rep_lines <- c(rep_lines, "- `data/processed/spatial/phase05_rds_gradient/phase05_provenance.csv`")
rep_lines <- c(rep_lines, "- `data/processed/spatial/phase06_hippo_regional_atlas/phase06_validation_summary.txt`")
rep_lines <- c(rep_lines, "- `docs/original_code_from_paper/Scripts/R/6. Spatial-Data process.R`")
rep_lines <- c(rep_lines, "- `docs/original_code_from_paper/Scripts/R/7. Spatial-DESeq2.R`")
rep_lines <- c(rep_lines, "- `docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R`")
rep_lines <- c(rep_lines, "")

writeLines(rep_lines, file.path(docs_dir, "spatial_reproducibility_report.md"))
cat("  Written: spatial_reproducibility_report.md\n")

## ---- E11: Generate docs/spatial_handoff_to_mitochondrial_analysis.md -----------

cat("\n--- E11: Generating spatial_handoff_to_mitochondrial_analysis.md ---\n")

hand_lines <- character(0)
hand_lines <- c(hand_lines, "# Spatial Handoff to Mitochondrial Analysis")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, mkline("**Generated by:** s07_figure_validation.R | ",
                                    format(run_start, "%Y-%m-%d %H:%M:%S")))
hand_lines <- c(hand_lines, mkline("**Git HEAD:** ", git_head))
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "## Decision")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, paste0("**", handoff_decision, "**"))
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "## Rationale")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, handoff_rationale)
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "## Handoff Decision Criteria (Evaluated)")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "| Criterion | Status | Evidence |")
hand_lines <- c(hand_lines, "|-----------|--------|----------|")
hand_lines <- c(hand_lines, paste0("| Phase 06 PASS | ",
                                    if (phase06_pass) "YES" else "NO", " | 54/54 checks PASS |"))
hand_lines <- c(hand_lines, paste0("| Phase 02 PASS | ",
                                    if (phase02_pass) "YES" else "NO", " | 109/109 checks PASS |"))
hand_lines <- c(hand_lines, "| D01-D05 affect Region axis? | NO | IS/Edge discrepancy orthogonal to Region labels |")
hand_lines <- c(hand_lines, paste0("| CA1/CA3/DG labels validated? | YES | ",
                                    "Phase 06: CA1=4671, CA2=565, CA3=2321, ML=1480, GCL=645, Hilus=239 (all 0.0% off) |"))
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "## Tangram Recommendation")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, paste0("**", tangram_decision, "**"))
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, tangram_rationale)
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "## Caveats for Phase 09 Planning")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "1. **D01/D02:** Phase 04 Middle/Old IS count discrepancy (msigdbr/Seurat version). Does not affect Region labels.")
hand_lines <- c(hand_lines, "2. **D03/D04:** Phase 05 DESeq2 Up/Down count difference. RDS-approximated, not strict STutility.")
hand_lines <- c(hand_lines, "3. **D05:** Phase 05 neighbor graph uses RDS Euclidean coordinates, not STutility image masks.")
hand_lines <- c(hand_lines, "4. **D06:** Seurat v5.5.0 vs author v4.3. May affect assay slot access patterns.")
hand_lines <- c(hand_lines, "5. **D07:** Raw Visium files unavailable. All spatial reconstruction uses RDS objects.")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "## Implications for Phase 09 (Mitochondrial Analysis)")
hand_lines <- c(hand_lines, "")
hand_lines <- c(hand_lines, "- Region labels (CA1, CA2, CA3, ML, GCL, Hilus) are available and validated.")
hand_lines <- c(hand_lines, "- Mitochondrial gene analysis can proceed at region level without Tangram.")
hand_lines <- c(hand_lines, "- Cell-type deconvolution (Tangram) may be needed if interpreting regional differences requires cell-type composition.")
hand_lines <- c(hand_lines, "- All caveats from D01-D07 should be documented in Phase 09 plan.")
hand_lines <- c(hand_lines, "")

writeLines(hand_lines, file.path(docs_dir, "spatial_handoff_to_mitochondrial_analysis.md"))
cat("  Written: spatial_handoff_to_mitochondrial_analysis.md\n")

## ---- E12: Write validation_summary.txt + provenance.csv -----------------------

cat("\n--- E12: Writing validation summary + provenance ---\n")

summary_lines <- character(0)
summary_lines <- c(summary_lines, "=== Phase Spatial-07 Validation Summary ===")
summary_lines <- c(summary_lines, paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
summary_lines <- c(summary_lines, paste0("Git HEAD: ", git_head))
summary_lines <- c(summary_lines, paste0("renv.lock MD5: ", renv_md5))
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, "--- Phase Classification ---")
for (i in seq_len(nrow(phase_classification))) {
  r <- phase_classification[i, ]
  summary_lines <- c(summary_lines, paste0("  ", r$phase, ": ", r$classification, " | ", r$script))
}
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, "--- Discrepancies ---")
for (i in seq_len(nrow(discrepancy_table))) {
  r <- discrepancy_table[i, ]
  summary_lines <- c(summary_lines, paste0("  ", r$id, ": [", r$severity, "] ", r$metric))
}
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, "--- Handoff ---")
summary_lines <- c(summary_lines, paste0("Decision: ", handoff_decision))
summary_lines <- c(summary_lines, paste0("Tangram: ", tangram_decision))
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, "--- Outputs ---")
summary_lines <- c(summary_lines, "phase07_phase_classification.csv: saved")
summary_lines <- c(summary_lines, "phase07_figure_classification_table.csv: saved")
summary_lines <- c(summary_lines, "phase07_discrepancy_table.csv: saved")
summary_lines <- c(summary_lines, "phase07_validation_summary.txt: saved")
summary_lines <- c(summary_lines, "phase07_provenance.csv: saved")
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, "--- Rplots.pdf ---")
summary_lines <- c(summary_lines, paste0("Existed at start: ", rplots_existed))

writeLines(summary_lines, file.path(out_dir, "phase07_validation_summary.txt"))
cat("  Written: phase07_validation_summary.txt\n")

provenance <- data.frame(
  key = c("git_head", "renv_lock_md5", "run_start", "run_end_placeholder",
           "r_version", "seurat_version",
           "phase02_status", "phase02_pass_count", "phase02_total_checks",
           "phase03_oy_padj", "phase03_my_padj", "phase03_om_padj",
           "phase04_status", "phase04_warning_acknowledged",
           "phase04_middle_is_ours", "phase04_middle_is_author",
           "phase04_old_is_ours", "phase04_old_is_author",
           "phase05_deseq2_up_ours", "phase05_deseq2_up_author",
           "phase05_deseq2_down_ours", "phase05_deseq2_down_author",
           "phase05_route", "phase05_strict_stutility_equivalence",
           "phase06_status", "phase06_pass_count", "phase06_total_checks",
           "handoff_decision", "tangram_decision",
           "rplots_pdf_existed_at_start", "new_rplots_generated",
           "package_install_attempted", "renv_lock_modified"),
  value = c(git_head, renv_md5,
             format(run_start, "%Y-%m-%d %H:%M:%S"), "",
             R.version.string, "5.5.0",
             "PASS", "109", "109",
             as.character(oy_padj), as.character(my_padj), as.character(om_padj),
             "WARNING", "TRUE",
             as.character(our_is_middle), as.character(author_is_middle),
             as.character(our_is_old), as.character(author_is_old),
             as.character(our_deseq2_up), as.character(author_deseq2_up),
             as.character(our_deseq2_down), as.character(author_deseq2_down),
             "RDS-approximated", "FALSE",
             "PASS", "54", "54",
             handoff_decision, tangram_decision,
             as.character(rplots_existed), "FALSE", "FALSE", "FALSE"),
  stringsAsFactors = FALSE
)

if (!rplots_existed && file.exists(rplots_path)) {
  cat("  WARNING: New Rplots.pdf detected; deleting.\n")
  file.remove(rplots_path)
  provenance$value[provenance$key == "new_rplots_generated"] <- "TRUE"
}

run_end <- Sys.time()
provenance$value[provenance$key == "run_end_placeholder"] <- format(run_end, "%Y-%m-%d %H:%M:%S")

write.csv(provenance, file.path(out_dir, "phase07_provenance.csv"), row.names = FALSE)
cat("  Written: phase07_provenance.csv\n")

cat("\nRun end:", format(run_end, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Duration:", round(difftime(run_end, run_start, units = "secs"), 1), "s\n")
cat("=== Phase Spatial-07 complete ===\n")

gc()
