#!/usr/bin/env Rscript

# =========================================================================
# Phase Spatial-14: Candidate Regulator-Associated Gene Report Generator
# Report integration only — NO new statistical analysis
# Reads Phase 13/13b/13c outputs, generates report bundle
# =========================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ---- E0: Setup & Validation ----
cat("Phase 14: Candidate Regulator Report Generator\n")
cat("================================================\n")

SRC_DIR <- "data/processed/spatial/phase13_candidate_regulator"
FIG_DIR <- "figures/spatial/phase13_candidate_regulator"
RPT_DIR <- "data/report/phase14_candidate_regulator_report"

# Record renv.lock MD5 at script start
RENV_MD5 <- trimws(system("md5 -q renv.lock 2>/dev/null", intern = TRUE))
if (length(RENV_MD5) == 0) RENV_MD5 <- "NO_RENV_LOCK"
cat(sprintf("renv.lock MD5: %s\n", RENV_MD5))

# Verify no Rplots.pdf exists at start
if (file.exists("Rplots.pdf")) {
  stop("FATAL: Rplots.pdf already exists at script start")
}
cat("Rplots.pdf absent at start: OK\n")

# Verify source directory exists
stopifnot(dir.exists(SRC_DIR))
stopifnot(dir.exists(FIG_DIR))

# Verify Phase 13/13b/13c validation all PASS
v13 <- read.csv(file.path(SRC_DIR, "phase13_validation_summary.csv"), stringsAsFactors = FALSE)
v13b <- read.csv(file.path(SRC_DIR, "phase13b_display_validation_summary.csv"), stringsAsFactors = FALSE)
v13c <- read.csv(file.path(SRC_DIR, "phase13c_validation_summary.csv"), stringsAsFactors = FALSE)

stopifnot(sum(grepl("FAIL", v13$status, ignore.case = TRUE)) == 0)
stopifnot(sum(grepl("FAIL", v13b$status, ignore.case = TRUE)) == 0)
stopifnot(sum(grepl("FAIL", v13c$status, ignore.case = TRUE)) == 0)
cat("Phase 13 validation: PASS\nPhase 13b validation: PASS\nPhase 13c validation: PASS\n")

# Create report bundle directories
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "tables"), showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "tables", "copied_source_tables"), showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "figures", "CA1_CA3_overview"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "figures", "age_stratified"), showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "figures", "spatial_maps"), showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "provenance"), showWarnings = FALSE)
cat("Report bundle directories created.\n")

# Read required source data
cat("\nReading source data...\n")
top_ca1 <- read.csv(file.path(SRC_DIR, "phase13_top_candidates_CA1.csv"), stringsAsFactors = FALSE)
top_ca3 <- read.csv(file.path(SRC_DIR, "phase13_top_candidates_CA3.csv"), stringsAsFactors = FALSE)
top_concordant <- read.csv(file.path(SRC_DIR, "phase13_top_candidates_concordant.csv"), stringsAsFactors = FALSE)
top_flipped <- read.csv(file.path(SRC_DIR, "phase13_top_candidates_region_flipped.csv"), stringsAsFactors = FALSE)
age_class <- read.csv(file.path(SRC_DIR, "phase13c_age_specific_candidate_classes.csv"), stringsAsFactors = FALSE)
old_specific <- read.csv(file.path(SRC_DIR, "phase13c_old_specific_candidates.csv"), stringsAsFactors = FALSE)
age_shifted <- read.csv(file.path(SRC_DIR, "phase13c_age_shifted_candidates.csv"), stringsAsFactors = FALSE)
region_shifted <- read.csv(file.path(SRC_DIR, "phase13c_region_shifted_by_age_candidates.csv"), stringsAsFactors = FALSE)
cat("All source data loaded.\n")

manifest_entries <- list()

add_manifest <- function(orig, copied, ftype, reason, phase, caption = "") {
  manifest_entries[[length(manifest_entries) + 1]] <<- data.frame(
    original_path = orig,
    copied_path = copied,
    file_type = ftype,
    reason_included = reason,
    source_phase = phase,
    caption = caption,
    stringsAsFactors = FALSE
  )
}

# ---- E1: Generate Report Summary Tables (T1-T11) ----

cat("\nGenerating report summary tables...\n")

report_header <- "# Report Summary — Phase Spatial-14. Derived from Phase 13/13b/13c outputs. NOT a new statistical analysis."

# T1: Top CA1 all-age (top 50, exclude module_member_or_effector)
t1 <- top_ca1 %>%
  filter(self_correlation_flag != "module_member_or_effector") %>%
  slice_head(n = 50)
write.csv(t1, file.path(RPT_DIR, "tables", "phase14_top_CA1_all_age.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13_top_candidates_CA1.csv"),
  "tables/phase14_top_CA1_all_age.csv",
  "report_table", "Top 50 CA1 candidates by rank_CA1", "Phase 13",
  "CA1 all-age candidates filtered to exclude module_member_or_effector"
)
cat(sprintf("  T1: Top 50 CA1 all-age — %d rows\n", nrow(t1)))

# T2: Top CA3 all-age (top 50, exclude module_member_or_effector)
t2 <- top_ca3 %>%
  filter(self_correlation_flag != "module_member_or_effector") %>%
  slice_head(n = 50)
write.csv(t2, file.path(RPT_DIR, "tables", "phase14_top_CA3_all_age.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13_top_candidates_CA3.csv"),
  "tables/phase14_top_CA3_all_age.csv",
  "report_table", "Top 50 CA3 candidates by rank_CA3", "Phase 13b",
  "CA3 all-age candidates filtered to exclude module_member_or_effector"
)
cat(sprintf("  T2: Top 50 CA3 all-age — %d rows\n", nrow(t2)))

# T3: Top concordant (top 50, exclude module_member_or_effector)
t3 <- top_concordant %>%
  filter(self_correlation_flag != "module_member_or_effector") %>%
  slice_head(n = 50)
write.csv(t3, file.path(RPT_DIR, "tables", "phase14_top_concordant.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13_top_candidates_concordant.csv"),
  "tables/phase14_top_concordant.csv",
  "report_table", "Top 50 concordant candidates by combined_score", "Phase 13",
  "Concordant in both CA1 and CA3"
)
cat(sprintf("  T3: Top 50 concordant — %d rows\n", nrow(t3)))

# T4: Top region-flipped (top 30)
t4 <- head(top_flipped, 30)
write.csv(t4, file.path(RPT_DIR, "tables", "phase14_top_region_flipped.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13_top_candidates_region_flipped.csv"),
  "tables/phase14_top_region_flipped.csv",
  "report_table", "Top 30 region-flipped candidates by combined_score", "Phase 13",
  "Opposite signs in CA1 vs CA3"
)
cat(sprintf("  T4: Top %d region-flipped\n", nrow(t4)))

# T5: Top Old-specific CA1 (from old_specific, region=="CA1", top 30 by |rho_Old|)
t5 <- old_specific %>%
  filter(region == "CA1") %>%
  arrange(desc(abs(rho_Old))) %>%
  head(30)
write.csv(t5, file.path(RPT_DIR, "tables", "phase14_top_old_specific_CA1.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13c_old_specific_candidates.csv"),
  "tables/phase14_top_old_specific_CA1.csv",
  "report_table", "Top 30 CA1 old-specific candidates by |rho_Old|", "Phase 13c",
  "CA1 genes with strong mitochondrial module association in Old but not Young"
)
cat(sprintf("  T5: Top %d Old-specific CA1\n", nrow(t5)))

# T6: Top Old CA3-specific (from age_class, region_age_class == "Old_CA3_specific")
t6 <- age_class %>%
  filter(region_age_class == "Old_CA3_specific") %>%
  arrange(desc(abs(rho_CA3_Old))) %>%
  head(30)
write.csv(t6, file.path(RPT_DIR, "tables", "phase14_top_old_specific_CA3.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13c_old_specific_candidates.csv"),
  "tables/phase14_top_old_specific_CA3.csv",
  "report_table", "Top 30 CA3 old-specific candidates by |rho_Old|", "Phase 13c",
  "CA3 genes with strong mitochondrial module association in Old but not Young"
)
cat(sprintf("  T6: Top %d Old-specific CA3\n", nrow(t6)))

# T7: Top Old-enhanced coupling (age_class csv, old_enhanced_coupling_flag==TRUE)
t7 <- age_class %>%
  filter(old_enhanced_coupling_flag == TRUE) %>%
  mutate(max_rho_Old = pmax(abs(rho_CA1_Old), abs(rho_CA3_Old), na.rm = TRUE)) %>%
  arrange(desc(max_rho_Old)) %>%
    head(30)
write.csv(t7, file.path(RPT_DIR, "tables", "phase14_top_old_enhanced_coupling.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13c_age_specific_candidate_classes.csv"),
  "tables/phase14_top_old_enhanced_coupling.csv",
  "report_table", "Top 30 old-enhanced coupling candidates by max(|rho_Old|)", "Phase 13c",
  "Genes where Old rho > Young rho by threshold"
)
cat(sprintf("  T7: Top %d old-enhanced coupling\n", nrow(t7)))

# T8: Top region-shifted in Old (from region_shifted, age_stratum=="Old", top 30 by |delta_region|)
t8 <- region_shifted %>%
  filter(age_stratum == "Old") %>%
  arrange(desc(abs(delta_region))) %>%
    head(30)
write.csv(t8, file.path(RPT_DIR, "tables", "phase14_top_region_shifted_old.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13c_region_shifted_by_age_candidates.csv"),
  "tables/phase14_top_region_shifted_old.csv",
  "report_table", "Top 30 region-shifted candidates in Old by |delta_region|", "Phase 13c",
  "Genes showing CA1/CA3 difference in Old age"
)
cat(sprintf("  T8: Top %d region-shifted Old\n", nrow(t8)))

# T9: Top age-shifted (exclude CA3 Young-only, top 30 by |delta_rho|)
t9 <- age_shifted %>%
  filter(!(region == "CA3" & sample_size_label == "severely_underpowered_exploratory")) %>%
  arrange(desc(abs(delta_rho))) %>%
    head(30)
write.csv(t9, file.path(RPT_DIR, "tables", "phase14_top_age_shifted.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13c_age_shifted_candidates.csv"),
  "tables/phase14_top_age_shifted.csv",
  "report_table", "Top 30 age-shifted candidates by |delta_rho|, excluding CA3 Young-only", "Phase 13c",
  "Genes with rho sign reversal between Young and Old, CA3 Young-only excluded"
)
cat(sprintf("  T9: Top %d age-shifted (excl CA3 Young-only)\n", nrow(t9)))

# ---- T10: Prioritized Short List ----
# Tier criteria from plan §4.1-K
# 1. |rho| > 0.5 in Old CA1 or Old CA3
# 2. NOT module_member_or_effector
# 3. NOT solely driven by CA3 Young
# 4. self_correlation_flag == "regulator_candidate"

cat("\nBuilding prioritized short list (T10)...\n")

# Create a joined data frame: all-age top candidates + age-stratified classes
age_class_for_join <- age_class %>%
  select(gene, module, primary_age_class, region_age_class, candidate_class, secondary_class,
         rho_CA1_Young, rho_CA1_Old, rho_CA3_Young, rho_CA3_Old,
         candidate_score_CA1_Young, candidate_score_CA1_Old,
         candidate_score_CA3_Young, candidate_score_CA3_Old,
         DE_age_support_flag_CA1, DE_age_support_flag_CA3,
         old_enhanced_coupling_flag, young_enhanced_coupling_flag,
         self_correlation_flag)

# Use age_class as primary source with |rho| filter
candidates <- age_class_for_join %>%
  filter(self_correlation_flag == "regulator_candidate") %>%
  mutate(
    max_rho_Old = pmax(abs(rho_CA1_Old), abs(rho_CA3_Old), na.rm = TRUE),
    has_strong_Old = max_rho_Old > 0.5,
    is_CA3_young_only = (abs(rho_CA3_Young) > 0.3 & abs(rho_CA3_Old) <= 0.3 & abs(rho_CA1_Old) <= 0.3)
  ) %>%
  filter(has_strong_Old == TRUE) %>%
  filter(is_CA3_young_only == FALSE)

# Determine priority tier
candidates <- candidates %>%
  mutate(
    priority_tier = case_when(
      # Tier 1: old_region_conserved AND |rho|>0.6 in both Old
      region_age_class == "Old_region_conserved" &
        abs(rho_CA1_Old) > 0.6 & abs(rho_CA3_Old) > 0.6 ~ "Tier 1",
      # Tier 1b: Old-specific AND concordant (both regions)
      grepl("old_specific", primary_age_class, ignore.case = TRUE) &
        (grepl("Old_region_conserved|Old_CA1_specific|Old_CA3_specific", region_age_class)) &
        abs(rho_CA1_Old) > 0.6 & abs(rho_CA3_Old) > 0.6 ~ "Tier 1",
      # Tier 2: Old-specific (either region) AND |rho|>0.6
      grepl("old_specific", primary_age_class, ignore.case = TRUE) &
        max_rho_Old > 0.6 ~ "Tier 2",
      # Tier 2b: concordant (all-age) AND appears in Old-specific
      grepl("age_conserved", primary_age_class, ignore.case = TRUE) &
        abs(rho_CA1_Old) > 0.6 ~ "Tier 2",
      # Tier 3: old_enhanced_coupling AND |rho_Old|>0.5
      old_enhanced_coupling_flag == TRUE &
        max_rho_Old > 0.5 ~ "Tier 3",
      # Tier 4: region_flipped AND |rho|>0.6 in both
      grepl("Old_CA1_specific|Old_CA3_specific", region_age_class) &
        abs(rho_CA1_Old) > 0.6 & abs(rho_CA3_Old) > 0.6 &
        sign(rho_CA1_Old) != sign(rho_CA3_Old) ~ "Tier 4",
      # Tier 5: old_specific or age_conserved with 0.3<|rho|<0.5
      (grepl("old_specific", primary_age_class, ignore.case = TRUE) |
         grepl("age_conserved", primary_age_class, ignore.case = TRUE)) &
        max_rho_Old > 0.3 & max_rho_Old <= 0.5 ~ "Tier 5",
      TRUE ~ "Unranked"
    ),
    notes_for_reviewer = case_when(
      priority_tier == "Tier 1" ~ "Highest priority: strong Old signal confirmed in both regions",
      priority_tier == "Tier 2" ~ "High priority: strong Old signal in at least one region",
      priority_tier == "Tier 3" ~ "Moderate priority: rho enhanced in Old vs Young",
      priority_tier == "Tier 4" ~ "Interesting biology: opposite signs across regions in Old",
      priority_tier == "Tier 5" ~ "Low priority: weaker rho that may still warrant checking",
      TRUE ~ "Below ranking threshold"
    )
  ) %>%
  arrange(factor(priority_tier, levels = c("Tier 1","Tier 2","Tier 3","Tier 4","Tier 5","Unranked")),
          desc(max_rho_Old))

# Add report-level metadata columns
t10 <- candidates %>%
  select(gene, module, primary_age_class, region_age_class, candidate_class,
         rho_CA1_Young, rho_CA1_Old, rho_CA3_Young, rho_CA3_Old,
         max_rho_Old, priority_tier, notes_for_reviewer,
         DE_age_support_flag_CA1, DE_age_support_flag_CA3,
         old_enhanced_coupling_flag, young_enhanced_coupling_flag)

# Add metadata as column comments (write a separate metadata header)
write.csv(t10, file.path(RPT_DIR, "tables", "phase14_prioritized_short_list.csv"), row.names = FALSE)
add_manifest(
  file.path(SRC_DIR, "phase13c_age_specific_candidate_classes.csv"),
  "tables/phase14_prioritized_short_list.csv",
  "report_table", "Prioritized short list for literature review (Tier 1-5)", "Phase 13+13c",
  "report_level_prioritization_only=TRUE; no_new_statistical_test=TRUE. Excludes module_member_or_effector and CA3 Young-only evidence."
)
cat(sprintf("  T10: Prioritized short list — %d rows\n", nrow(t10)))
cat(sprintf("    Tier 1: %d, Tier 2: %d, Tier 3: %d, Tier 4: %d, Tier 5: %d\n",
            sum(t10$priority_tier == "Tier 1"),
            sum(t10$priority_tier == "Tier 2"),
            sum(t10$priority_tier == "Tier 3"),
            sum(t10$priority_tier == "Tier 4"),
            sum(t10$priority_tier == "Tier 5")))

# T11: Unique genes (deduplicated)
t11 <- t10 %>%
  group_by(gene) %>%
  slice_max(order_by = max_rho_Old, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(factor(priority_tier, levels = c("Tier 1","Tier 2","Tier 3","Tier 4","Tier 5","Unranked")),
          desc(max_rho_Old))
write.csv(t11, file.path(RPT_DIR, "tables", "phase14_prioritized_short_list_unique_genes.csv"), row.names = FALSE)
add_manifest(
  file.path(RPT_DIR, "tables", "phase14_prioritized_short_list.csv"),
  "tables/phase14_prioritized_short_list_unique_genes.csv",
  "report_table", "Prioritized short list deduplicated to unique genes", "Phase 13+13c",
  "report_level_prioritization_only=TRUE; no_new_statistical_test=TRUE. One row per gene (best module per gene)."
)
cat(sprintf("  T11: Unique gene prioritized list — %d genes\n", nrow(t11)))

# Top 10 unique genes for quick reference
cat("\nTop 10 candidate genes by priority:\n")
t11_top <- t11 %>% slice_head(n = 10)
for (i in seq_len(nrow(t11_top))) {
  cat(sprintf("  %s (%s) | rho_Old=%.2f | %s | %s\n",
              t11_top$gene[i], t11_top$module[i],
              t11_top$max_rho_Old[i],
              t11_top$priority_tier[i],
              t11_top$primary_age_class[i]))
}

# ---- E2: Copy Figures ----

cat("\nCopying selected figures...\n")

fig_copies <- list(
  # Phase 13 CA1-focused
  list("F01_regulator_score_heatmap.png", "figures/CA1_CA3_overview/F01_CA1_regulator_score_heatmap.png", "Phase 13"),
  list("F03_top_regulators_dotplot.png", "figures/CA1_CA3_overview/F03_CA1_top_regulators_dotplot.png", "Phase 13"),
  list("F02_rho_vs_DE_scatter.png", "figures/CA1_CA3_overview/F02_rho_vs_DE_scatter.png", "Phase 13"),
  list("F06_candidate_classes.png", "figures/CA1_CA3_overview/F06_candidate_classes.png", "Phase 13"),
  list("F09_score_components.png", "figures/CA1_CA3_overview/F09_CA1_score_components.png", "Phase 13"),
  # Phase 13b CA3-focused
  list("F01b_CA3_regulator_score_heatmap.png", "figures/CA1_CA3_overview/F01b_CA3_regulator_score_heatmap.png", "Phase 13b"),
  list("F03b_CA3_top_regulators_dotplot.png", "figures/CA1_CA3_overview/F03b_CA3_top_regulators_dotplot.png", "Phase 13b"),
  list("F09b_CA3_score_components.png", "figures/CA1_CA3_overview/F09b_CA3_score_components.png", "Phase 13b"),
  list("F10b_CA3_sensitivity_comparison.png", "figures/CA1_CA3_overview/F10b_CA3_sensitivity_comparison.png", "Phase 13b"),
  # Spatial maps
  list("F14_spatial_Dach1.png", "figures/spatial_maps/F14_spatial_CA1_Dach1.png", "Phase 13"),
  list("F14b_spatial_CA3_Elk1.png", "figures/spatial_maps/F14b_spatial_CA3_Elk1.png", "Phase 13b"),
  # Phase 13c Age-stratified (7 figures)
  list("phase13c_heatmap_CA1_young.png", "figures/age_stratified/phase13c_heatmap_CA1_young.png", "Phase 13c"),
  list("phase13c_heatmap_CA1_old.png", "figures/age_stratified/phase13c_heatmap_CA1_old.png", "Phase 13c"),
  list("phase13c_heatmap_CA3_old.png", "figures/age_stratified/phase13c_heatmap_CA3_old.png", "Phase 13c"),
  list("phase13c_heatmap_CA3_young.png", "figures/age_stratified/phase13c_heatmap_CA3_young.png", "Phase 13c"),
  list("phase13c_delta_CA1_old_minus_young.png", "figures/age_stratified/phase13c_delta_CA1_old_minus_young.png", "Phase 13c"),
  list("phase13c_region_age_class_barplot.png", "figures/age_stratified/phase13c_region_age_class_barplot.png", "Phase 13c"),
  list("phase13c_old_specific_top_dotplot.png", "figures/age_stratified/phase13c_old_specific_top_dotplot.png", "Phase 13c")
)

copied_figs <- 0
skipped_figs <- 0

for (fc in fig_copies) {
  src <- file.path(FIG_DIR, fc[[1]])
  dst <- file.path(RPT_DIR, fc[[2]])
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    add_manifest(src, fc[[2]], "figure", paste("Copied for report figure guide"), fc[[3]], fc[[1]])
    copied_figs <- copied_figs + 1
  } else {
    cat(sprintf("  WARN: Missing figure: %s\n", src))
    skipped_figs <- skipped_figs + 1
    add_manifest(src, fc[[2]], "figure_missing", "Figure not found at source path", fc[[3]], fc[[1]])
  }
}
cat(sprintf("  Copied %d figures, %d skipped (missing)\n", copied_figs, skipped_figs))

# ---- E3: Copy Source Tables (Optional Head Copies) ----

cat("\nCopying source table head samples...\n")

# Head of key source tables for reference
source_tables_to_copy <- c(
  "phase13_regulator_score_default.csv",
  "phase13c_age_specific_candidate_classes.csv",
  "phase13c_validation_summary.csv"
)

for (st in source_tables_to_copy) {
  src <- file.path(SRC_DIR, st)
  dst <- file.path(RPT_DIR, "tables", "copied_source_tables", st)
  if (file.exists(src)) {
    data <- read.csv(src, stringsAsFactors = FALSE, nrows = 100)
    write.csv(data, dst, row.names = FALSE)
    add_manifest(src, file.path("tables/copied_source_tables", st),
                 "source_table_head", "First 100 rows of source table for reference",
                 "Phase 13/13c", "Head copy only. Full file linked in report by original path.")
  }
}
cat("  Source table heads copied.\n")

# ---- E4: Write Manifest CSV ----

cat("\nWriting manifest...\n")
manifest_df <- do.call(rbind, manifest_entries)
write.csv(manifest_df, file.path(RPT_DIR, "manifest.csv"), row.names = FALSE)
cat(sprintf("  Manifest: %d entries\n", nrow(manifest_df)))

# ---- E5: Write Provenance Summary ----

cat("\nWriting provenance summary...\n")

prov_lines <- c(
  sprintf("Phase Spatial-14 Report Generation Summary"),
  sprintf("========================================="),
  sprintf("Generated: %s", Sys.time()),
  sprintf("Script: R/spatial/s14_generate_candidate_regulator_report.R"),
  sprintf("renv.lock MD5: %s", RENV_MD5),
  sprintf(""),
  sprintf("Source Phases:"),
  sprintf("  Phase 13: %d PASS / 0 FAIL", sum(v13$status == "PASS")),
  sprintf("  Phase 13b: %d PASS / 0 FAIL", sum(v13b$status == "PASS")),
  sprintf("  Phase 13c: %d PASS / 0 FAIL", sum(v13c$status == "PASS")),
  sprintf(""),
  sprintf("Sample Sizes:"),
  sprintf("  CA1 Young: n=4 (exploratory_small_n)"),
  sprintf("  CA1 Old: n=8 (exploratory_moderate_n)"),
  sprintf("  CA3 Young: n=3 (severely_underpowered_exploratory)"),
  sprintf("  CA3 Old: n=8 (exploratory_moderate_n)"),
  sprintf(""),
  sprintf("Report Bundle Contents:"),
  sprintf("  Report tables: %d (T1-T11)", 11),
  sprintf("  Copied figures: %d (skipped: %d)", copied_figs, skipped_figs),
  sprintf("  Source table heads: %d", length(source_tables_to_copy)),
  sprintf("  Manifest: manifest.csv (%d entries)", nrow(manifest_df)),
  sprintf(""),
  sprintf("Prioritized Short List (T10):"),
  sprintf("  Total candidates: %d", nrow(t10)),
  sprintf("  Tier 1: %d", sum(t10$priority_tier == "Tier 1")),
  sprintf("  Tier 2: %d", sum(t10$priority_tier == "Tier 2")),
  sprintf("  Tier 3: %d", sum(t10$priority_tier == "Tier 3")),
  sprintf("  Tier 4: %d", sum(t10$priority_tier == "Tier 4")),
  sprintf("  Tier 5: %d", sum(t10$priority_tier == "Tier 5")),
  sprintf("  Unique genes (T11): %d", nrow(t11)),
  sprintf(""),
  sprintf("CAVEATS:"),
  sprintf("  - No new statistical analysis performed. All data derived from Phase 13/13b/13c."),
  sprintf("  - CA3 Young n=3: severely_underpowered_exploratory."),
  sprintf("  - No causal regulatory claims made. Candidates are hypothesis-generating only."),
  sprintf("  - T10/T11: report_level_prioritization_only=TRUE, no_new_statistical_test=TRUE."),
  sprintf(""),
  sprintf("Validation: See phase14_validation_checks below.")
)

writeLines(prov_lines, file.path(RPT_DIR, "provenance", "phase14_report_generation_summary.txt"))
add_manifest(
  "R/spatial/s14_generate_candidate_regulator_report.R",
  "provenance/phase14_report_generation_summary.txt",
  "provenance", "Report generation summary", "Phase 14",
  "Auto-generated by Phase 14 script"
)

# Write source file manifest
src_manifest <- data.frame(
  source_file = c(
    "phase13_top_candidates_CA1.csv",
    "phase13_top_candidates_CA3.csv",
    "phase13_top_candidates_concordant.csv",
    "phase13_top_candidates_region_flipped.csv",
    "phase13c_age_specific_candidate_classes.csv",
    "phase13c_old_specific_candidates.csv",
    "phase13c_age_shifted_candidates.csv",
    "phase13c_region_shifted_by_age_candidates.csv",
    "phase13_validation_summary.csv",
    "phase13b_display_validation_summary.csv",
    "phase13c_validation_summary.csv"
  ),
  source_phase = c(rep("Phase 13", 4), rep("Phase 13c", 4), "Phase 13", "Phase 13b", "Phase 13c"),
  rows = c(nrow(top_ca1), nrow(top_ca3), nrow(top_concordant), nrow(top_flipped),
           nrow(age_class), nrow(old_specific), nrow(age_shifted), nrow(region_shifted),
           nrow(v13), nrow(v13b), nrow(v13c)),
  status = "used for report generation"
)
write.csv(src_manifest, file.path(RPT_DIR, "provenance", "source_file_manifest.csv"), row.names = FALSE)

cat("\nPhase 14 script complete.\n")
cat(sprintf("Report bundle: %s/\n", RPT_DIR))
