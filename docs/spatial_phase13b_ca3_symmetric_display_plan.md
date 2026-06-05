# Phase Spatial-13b: CA3 Symmetric Display & Guide/Report Correction Plan

**Version**: 1.0
**Date**: 2026-06-05
**Status**: PLAN-ONLY (no execution; no R scripts; no data/figures generated)
**Prerequisite**: Phase 13 executed and validated (12 PASS / 0 FAIL / 0 WARN)
**Relationship**: Phase 13b is a **supplementary display/report correction** on top of Phase 13. It does NOT re-run Phase 13's core computation (rho, regulator_score, DE logic). All outputs **merge back into existing Phase 13 directories**.

---

## 1. Background & Problem Statement

### 1.1 What Phase 13 Did Correctly

Phase 13 (`s13_candidate_regulator_discovery.R`, 1470 lines) computed CA1 and CA3 symmetrically:

- `rho_CA1` and `rho_CA3` — Spearman correlation per region
- `score_CA1_default` and `score_CA3_default` — composite regulator_score per region
- `rank_CA1` and `rank_CA3` — per-module ranking per region
- `DE_age_lfc_CA1` and `DE_age_lfc_CA3` — Old-vs-Young log2FC per region
- `candidate_class` — includes both `CA1_aging_associated` (3,027 pairs) and `CA3_aging_associated` (3,290 pairs), plus `concordant` (1,573), `region_flipped` (674)

The **underlying data tables are symmetric**: every gene×module pair has both CA1 and CA3 columns.

### 1.2 The CA1 Display Bias

The display layer (figures, spatial plots, GO enrichment, top-candidate discussion) is **heavily CA1-biased**:

1. **F01** (regulator_score heatmap, line 614-651): Sorts by `rank_CA1`, uses `score_CA1_default`. Title says "CA1" but filename does not.
2. **F03** (top regulators dotplot, line 696-724): Uses `rank_CA1`, `abs(rho_CA1)`, `score_CA1_default`. Title says "(CA1, default config)" but filename says `F03_top_regulators_dotplot` — implying "global" top regulators.
3. **F04** (module correlation matrix, line 728-753): Uses `rho_CA1`. "CA1, top 50 genes" — CA1 only. No CA3 equivalent.
4. **F07** (Young stratified heatmap, line 797-854): Gene selection by `rank_CA1`, CA1 samples only. No CA3 equivalent.
5. **F08** (Old stratified heatmap, line 835-854): Gene selection by `rank_CA1`, CA1 samples only. No CA3 equivalent.
6. **F09** (regulator_score component breakdown, line 858-883): Uses `abs_rho_scaled_CA1`, `abs_DE_age_scaled_CA1`, `score_CA1_default` — CA1 only.
7. **F10** (sensitivity comparison, line 887-913): Uses `score_CA1_default`, `score_CA1_de_emphasis`, `score_CA1_rho_only` — CA1 only.
8. **F13** (DE volcano, line 955-979): Uses CA1 DE results (`de_ca1`). Title says "CA1, Old vs Young" but filename is `F13_DE_volcano_regulators` — ambiguous.
9. **F14** (spatial plots, line 994-1006): Top 5 by `score_CA1_default`. All 5 spatial genes (Dach1, Dbp, Elk1, Pou3f1, Xbp1) are CA1-prioritized.
10. **GO enrichment** (E12, line 1020-1079): Uses `rank_CA1` to select top genes per module. No CA3 GO enrichment.
11. **Analysis guide** (`spatial_phase13_candidate_regulator_analysis_guide.md`, line 1-102): Never explicitly distinguishes CA1-focused from CA3-focused outputs. The "top regulators" language is ambiguous. The guide mentions `rank_CA1/CA3` columns exist but doesn't note which figures use which.

### 1.3 What This Phase Fixes

Phase 13b is a **display/report supplement only**. It:

- Adds CA3 mirror figures for each CA1-only figure
- Generates explicit top-candidate CSVs for CA1, CA3, concordant, and region_flipped classes
- Updates the analysis guide to clearly label CA1-focused vs CA3-focused outputs
- Optionally appends a Phase 13b section to the existing guide

It does **NOT**:

- Re-compute rho, regulator_score, or candidate_class
- Create new output directories
- Install new packages
- Modify renv.lock
- Delete or rename existing Phase 13 outputs
- Load the Hippo RDS unnecessarily (only if spatial maps for CA3 top candidates are needed)
- Change any core ranking logic

---

## 2. Input Files

All input comes from existing Phase 13 outputs. No Phase 11 or Phase 12 file re-reading is needed (DE values are already in the Phase 13 score table).

| # | File | Purpose |
|---|------|---------|
| I1 | `data/processed/spatial/phase13_candidate_regulator/phase13_regulator_score_default.csv` | Primary ranking: 13,179 rows × 12 columns. Contains `rho_CA1`, `rho_CA3`, `score_CA1_default`, `score_CA3_default`, `rank_CA1`, `rank_CA3`, `DE_age_lfc_CA1`, `DE_age_lfc_CA3`, `DE_region_lfc`, `candidate_class`, `self_correlation_flag`. |
| I2 | `data/processed/spatial/phase13_candidate_regulator/phase13_regulator_score_sensitivity.csv` | All 3 weight configs per gene×module. Used for CA3 sensitivity comparison figure. |
| I3 | `data/processed/spatial/phase13_candidate_regulator/phase13_candidate_classes.csv` | Candidate class per gene×module pair. Used to cross-validate concordant/flipped classes. |
| I4 | `data/processed/spatial/phase13_candidate_regulator/phase13_spearman_rho_universe_regulator.csv` | Full regulator rho table with per-row DE. Used if additional columns (e.g., abs_rho_scaled_CA3) are needed for F09b, F10b. |
| I5 | `data/processed/spatial/phase13_candidate_regulator/phase13_final_report.txt` | Existing report text. Phase 13b will NOT overwrite this; will read for reference only. Optionally appends an addendum file. |
| I6 | `data/processed/spatial/phase13_candidate_regulator/phase13_validation_summary.csv` | 12 PASS / 0 WARN / 0 FAIL. Must verify before proceeding. |
| I7 | (Optional) `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` | Hippo RDS. Loaded ONLY if generating CA3 spatial maps (F14b). If loaded, unloaded immediately after. See §7 (Stop Conditions) for memory guard. |
| I8 | (Optional) `docs/target_genes.csv` | Target gene list. Used if generating target-gene-context annotations for top-candidate CSVs. |

---

## 3. Output Files (All Merge Into Existing Phase 13 Directories)

**CRITICAL**: All output paths point to existing Phase 13 directories. Do NOT create `phase13b_*` directories.

### 3.1 Data Outputs

All under `data/processed/spatial/phase13_candidate_regulator/`:

| # | File | Description |
|---|------|-------------|
| O1 | `phase13_top_candidates_CA1.csv` | CA1 top candidates: sorted by `rank_CA1`, then `score_CA1_default`. Columns: gene, module, candidate_class, rho_CA1, rho_CA3, score_CA1_default, score_CA3_default, rank_CA1, rank_CA3, DE_age_lfc_CA1, DE_age_lfc_CA3, DE_region_lfc, self_correlation_flag. |
| O2 | `phase13_top_candidates_CA3.csv` | CA3 top candidates: sorted by `rank_CA3`, then `score_CA3_default`. Same columns as O1. |
| O3 | `phase13_top_candidates_concordant.csv` | Concordant candidates: `candidate_class == "concordant"`, sorted by `pmax(score_CA1_default, score_CA3_default *)`. Same columns as O1, plus `combined_score = score_CA1_default + score_CA3_default`. |
| O4 | `phase13_top_candidates_region_flipped.csv` | Region-flipped candidates: `candidate_class == "region_flipped"`, sorted by `pmax(score_CA1_default, score_CA3_default)`. Same columns as O1, plus `combined_score`. Flag `rho_CA1` and `rho_CA3` signs prominently. |
| O5 | `phase13b_display_validation_summary.csv` | Phase 13b validation checks. Separate from Phase 13's validation (does NOT overwrite `phase13_validation_summary.csv`). |
| O6 | (Optional) `phase13_final_report_phase13b_addendum.txt` | If appending to report: a separate phase13b addendum file. Does NOT overwrite `phase13_final_report.txt`. |

**CSV schema for O1-O4**:

```
gene, module, candidate_class, rho_CA1, rho_CA3, score_CA1_default, score_CA3_default,
rank_CA1, rank_CA3, DE_age_lfc_CA1, DE_age_lfc_CA3, DE_region_lfc, self_correlation_flag
```

O3 and O4 add: `combined_score` (defined as `score_CA1_default + score_CA3_default`).

### 3.2 Figure Outputs

All under `figures/spatial/phase13_candidate_regulator/`. Existing CA1-focused figures remain untouched. New CA3-focused figures use `b` suffix:

**Group A — Existing CA1-only figures (rename plan is OPTIONAL):**

| Existing File | (Optional) Renamed | Content |
|---------------|-------------------|---------|
| `F01_regulator_score_heatmap.png` | `F01a_regulator_score_heatmap_CA1.png` | CA1 regulatory score heatmap |
| `F03_top_regulators_dotplot.png` | `F03a_top_regulators_dotplot_CA1.png` | CA1 top regulators dotplot |
| `F07_young_stratified_heatmap.png` | `F07a_young_stratified_heatmap_CA1.png` | CA1 Young stratified |
| `F08_old_stratified_heatmap.png` | `F08a_old_stratified_heatmap_CA1.png` | CA1 Old stratified |
| `F09_score_components.png` | `F09a_score_components_CA1.png` | CA1 score component breakdown |
| `F10_sensitivity_comparison.png` | `F10a_sensitivity_comparison_CA1.png` | CA1 sensitivity comparison |
| `F13_DE_volcano_regulators.png` | `F13a_DE_volcano_regulators_CA1.png` | CA1 DE volcano |
| `F14_spatial_<gene>.png` (5 files) | (keep as-is with note in guide: CA1 top-5) | CA1 spatial features |

**Default plan**: Keep existing file names as-is (they are already generated and gitignored). Do NOT rename. Add new CA3 versions with explicit `b` suffix. The analysis guide will list which figures are CA1 and which are CA3.

**Group B — NEW CA3 mirror figures (required):**

| # | File | Content |
|---|------|---------|
| F01b | `F01b_CA3_regulator_score_heatmap.png` | CA3 regulator_score heatmap (top 20 per module by `rank_CA3`, using `score_CA3_default`) |
| F03b | `F03b_CA3_top_regulators_dotplot.png` | CA3 top regulators dotplot (by `rank_CA3` / `score_CA3_default`) |
| F04b | `F04b_CA3_module_correlation_matrix.png` | CA3 gene×module ρ matrix (top 50 by max `|rho_CA3|`) |
| F07b | `F07b_CA3_young_stratified_heatmap.png` | CA3 Young stratified ρ heatmap (n=3 caveat displayed) |
| F08b | `F08b_CA3_old_stratified_heatmap.png` | CA3 Old stratified ρ heatmap |
| F09b | `F09b_CA3_score_components.png` | CA3 regulator_score component breakdown |
| F10b | `F10b_CA3_sensitivity_comparison.png` | CA3 sensitivity: default vs DE-emphasis vs ρ-only |
| F13b | `F13b_CA3_DE_volcano_regulators.png` | CA3 DE volcano (Old vs Young, CA3) with regulators highlighted |
| F14b | `F14b_spatial_CA3_<gene>.png` (5 files) | SpatialFeaturePlot for top 5 CA3-ranked candidates |

**Group C — Existing figures that already have both CA1+CA3 panels (rename optional, keep as-is):**

| File | Content |
|------|---------|
| `F02_rho_vs_DE_scatter.png` | Has F02a (CA1) + F02b (CA3) panels — no change needed |
| `F05_self_correlation_flags.png` | Region-agnostic — no change needed |
| `F06_candidate_classes.png` | Region-agnostic — no change needed |
| `F11_tf_annotation_coverage.png` | Region-agnostic — no change needed |
| `F12_universe_composition.png` | Region-agnostic — no change needed |
| `F15_GO_enrichment.png` | CA1 gene list (E12 uses `rank_CA1`). Could optionally add `F15b_CA3_GO_enrichment.png` |

### 3.3 Documentation Outputs

| # | File | Action |
|---|------|--------|
| D1 | `docs/spatial_phase13_candidate_regulator_analysis_guide.md` | **UPDATE** (not overwrite). Add a "Phase 13b: CA3 Symmetric Display" section at the end, or replace ambiguous "top regulators" language with explicit CA1/CA3 labeling. |
| D2 | `docs/spatial_phase13b_ca3_symmetric_display_plan.md` | **This file**. Authoritative Phase 13b plan. |
| D3 | (Optional) `docs/spatial_phase13_candidate_regulator_discovery_plan.md` | Optionally append a §19 Phase 13b section referencing this plan. |

---

## 4. Required Display/Report Corrections

### 4.1 Sorting and Selection Rules

| Output | Sort By (Primary) | Sort By (Secondary) | Cutoff |
|--------|-------------------|---------------------|--------|
| `phase13_top_candidates_CA1.csv` | `rank_CA1` ascending | `score_CA1_default` descending | All valid (non-NA rank) rows |
| `phase13_top_candidates_CA3.csv` | `rank_CA3` ascending | `score_CA3_default` descending | All valid (non-NA rank) rows |
| `phase13_top_candidates_concordant.csv` | `combined_score = score_CA1_default + score_CA3_default` descending | `pmax(score_CA1, score_CA3)` descending | `candidate_class == "concordant"` |
| `phase13_top_candidates_region_flipped.csv` | `combined_score` descending | `pmax(score_CA1, score_CA3)` descending | `candidate_class == "region_flipped"` |

**Pair-level vs unique-gene clarification**: Every row in these CSVs is a gene×module pair. A gene appearing in rows for 5 modules counts as 5 records. The header comment or first row note should state: "CAUTION: This table is gene×module pairs, not unique genes. See `phase13_final_report.txt` for unique gene counts."

### 4.2 New Figure Specifications

All new figures mirror the Phase 13 script's figure-generating code exactly, substituting:
- `rank_CA1` → `rank_CA3`
- `score_CA1_default` → `score_CA3_default`
- `rho_CA1` → `rho_CA3`
- `DE_age_lfc_CA1` → `DE_age_lfc_CA3`
- `abs_rho_scaled_CA1` → `abs_rho_scaled_CA3`
- `abs_DE_age_scaled_CA1` → `abs_DE_age_scaled_CA3`
- CA1 sample columns → CA3 sample columns (for F07b, F08b age-stratified)
- CA3 color (orange, per Phase 11/12 convention)
- Title: "CA1" → "CA3"
- n label: "n_CA3 = 11" instead of "n_CA1 = 12"

**CA3 Young n=3 caveat (F07b)**: The F07b title must include `(CAUTION: CA3 Young n=3)` and the subtitle must state `n = 3 (severely underpowered)`.

**Spatial maps (F14b)**: Top 5 unique genes by `rank_CA3` (across all modules, highest `score_CA3_default`). Use `SpatialFeaturePlot` with pt.size=0.3, same as Phase 13's E11.

### 4.3 Guide Corrections

The analysis guide must be updated to:

1. Add a section "CA1 vs CA3 Symmetric Outputs" that lists:
   - Which figures are CA1-focused (existing F01, F03, F04, F07, F08, F09, F10, F13, F14)
   - Which figures are CA3-focused (new F01b, F03b, F04b, F07b, F08b, F09b, F10b, F13b, F14b)
   - Which figures are region-agnostic (F05, F06, F11, F12)
   - Which figures have both regions (F02)

2. **Replace ambiguous language**: Where the guide says "top regulators" without qualification, add "(CA1-ranked)" or "(CA3-ranked)".

3. **Add Phase 13b caveats**:
   - "Phase 13b added CA3 mirror figures and top-candidate CSVs. The underlying data (rho, scores, ranks) was already computed for both regions in Phase 13. Phase 13b only corrected the display asymmetry."
   - "Top candidates in CSV files are gene×module pairs, not unique genes."

4. **Update the "Key Outputs" section to include**:
   - `phase13_top_candidates_CA1.csv`, `_CA3.csv`, `_concordant.csv`, `_region_flipped.csv`
   - CA3 mirror figures

### 4.4 Target Gene Context (Optional Enhancement)

If target gene annotation is desired, add a `target_gene_category` column to the top-candidate CSVs derived from `docs/target_genes.csv`, and add a `is_target_gene` annotation to the dotplot and heatmap figures. This is optional and does NOT block Phase 13b.

---

## 5. Validation Checks

### 5.1 Pre-execution checks (Phase 13b script start)

| # | Check | Expected | Severity |
|---|-------|----------|----------|
| V01 | `phase13_validation_summary.csv` exists and has 0 FAIL | All PASS (12/0/0) | FATAL |
| V02 | `phase13_regulator_score_default.csv` exists and has columns `score_CA3_default`, `rank_CA3` | Both columns present | FATAL |
| V03 | `phase13_regulator_score_default.csv` has > 0 valid `rank_CA3` values | > 0 non-NA values | FATAL |
| V04 | `phase13_candidate_classes.csv` exists and has `candidate_class` column | Present | FATAL |
| V05 | No `qval` used for `candidate_class` (spot-check) | qval columns present but `candidate_class` assignment uses ρ thresholds only (per Phase 13 E6: lines 481-488) | FATAL |
| V06 | renv.lock MD5 unchanged from Phase 13 | `b2bf05038b440576855d043e137c5883` | FATAL |
| V07 | No `phase13b_*` directory exists (output destination is `phase13_candidate_regulator/`) | `data/processed/spatial/phase13b_*` absent | WARN |

### 5.2 Post-execution checks (Phase 13b script end)

| # | Check | Expected | Severity |
|---|-------|----------|----------|
| V08 | O1-O4 CSVs exist and have > 0 rows | All 4 files present, non-zero rows | FATAL |
| V09 | `phase13b_display_validation_summary.csv` exists | Present | FATAL |
| V10 | F01b, F03b, F04b, F09b, F10b, F13b PNG files exist | At minimum: these 6 CA3 figures exist | FATAL |
| V11 | All new figures have `CA3` in filename or title | grep confirms | WARN |
| V12 | CA3 spatial maps (F14b) exist OR spatial skipped with WARNING | Either generated or explicitly skipped in log | WARN |
| V13 | Guide updated with CA1/CA3 labeling | grep for "CA3-focused" or "CA3 mirror" in guide | WARN |
| V14 | No `phase13b_*` output directory created | Directory check | FATAL |
| V15 | No existing Phase 13 files deleted or renamed | File existence pre/post check | FATAL |
| V16 | renv.lock unchanged after execution | MD5: `b2bf05038b440576855d043e137c5883` | FATAL |
| V17 | `Rplots.pdf` absent from repo root | `file.exists("Rplots.pdf") == FALSE` | FATAL |
| V18 | No full Seurat object saved | No new RDS in output directory | WARN |
| V19 | No causal language in new outputs | grep forbidden terms | FATAL |

---

## 6. Stop Conditions

### 6.1 Fatal Stop Conditions

| # | Condition | Action |
|---|----------|--------|
| FATAL-01 | `rank_CA3` or `score_CA3_default` column missing from `phase13_regulator_score_default.csv` | STOP — Cannot generate CA3-focused outputs. Phase 13 data may need re-run. |
| FATAL-02 | Phase 13 validation has FAIL | STOP — Fix Phase 13 issues first. |
| FATAL-03 | Any new output labels a gene as "causal regulator" or "upstream regulator" | STOP — Reword. All output must say "candidate regulator-associated gene". |
| FATAL-04 | Any attempt to re-compute rho, regulator_score, or candidate_class | STOP — Phase 13b is display-only. |
| FATAL-05 | Rplots.pdf detected at script start or generated during script | STOP — Fix PDF device handling. |
| FATAL-06 | renv.lock modified during execution | STOP — Phase 13b must NOT change packages. |

### 6.2 Non-Fatal Warnings

| # | Condition | Action |
|---|----------|--------|
| WARN-01 | CA3 Young n=3 in F07b: severely underpowered | Proceed; label `(CAUTION: CA3 Young n=3)` in figure title and subtitle. |
| WARN-02 | Hippo RDS spatial map generation exceeds memory or fails | Skip F14b spatial maps; proceed with all tabular/display corrections. Log WARNING. |
| WARN-03 | `phase13_spearman_rho_universe_regulator.csv` missing columns needed for F09b/F10b | Read columns from `phase13_regulator_score_default.csv` as fallback; F09b can use score values directly instead of scaled components. |
| WARN-04 | CA3 GO enrichment genes insufficient (< 5 per module) | Skip F15b; note in log. |

---

## 7. Execution Plan (for Future `s13b_ca3_symmetric_display.R`)

### 7.1 Script Outline

```
R/spatial/s13b_ca3_symmetric_display.R
```

**E0: Setup & Validation**
- Record renv.lock MD5 (must be `b2bf05038b440576855d043e137c5883`)
- Load dplyr, ggplot2, pheatmap, patchwork, stringr (no new packages)
- Verify all input files (I1-I6)
- Run pre-execution checks V01-V07
- If any FATAL, stop

**E1: Load Phase 13 Data**
- Read `phase13_regulator_score_default.csv`
- Read `phase13_regulator_score_sensitivity.csv`
- Read `phase13_candidate_classes.csv`
- Verify CA3 columns exist
- Verify no p-value/q-value used for candidate_class

**E2: Generate CA3 Top Candidates CSVs (O1-O4)**
- O1: Filter CA1 rows, sort by rank_CA1
- O2: Filter CA3 rows, sort by rank_CA3  
- O3: Filter `candidate_class == "concordant"`, sort by combined_score desc
- O4: Filter `candidate_class == "region_flipped"`, sort by combined_score desc
- Write 4 CSVs

**E3: CA3 Mirror Figures**
- See §4.2 specifications for each figure
- Use explicit `png()` with file path; optional `safe_pdf()` with tryCatch
- `gc()` between major figures

**E4: CA3 Spatial Maps (OPTIONAL — skip on memory failure)**
- Load Hippo RDS ONLY if generating
- Top 5 unique genes by rank_CA3 across all modules
- SpatialFeaturePlot per gene
- rm(hippo_rds); gc()

**E5: Guide Update**
- Read existing guide
- Insert Phase 13b section
- Fix ambiguous "top regulators" language
- Write updated guide

**E6: Validation & Cleanup**
- Run post-execution checks V08-V19
- Write `phase13b_display_validation_summary.csv`
- Rplots.pdf cleanup
- Report pass/warn/fail counts

### 7.2 Memory Management

- Phase 13b reads CSVs only (kilobytes). No large Seurat objects except optionally in E4.
- If E4 runs: load Hippo → plot → unload → gc(). Exactly one load/unload cycle.
- Each figure group: plot → ggsave → dev.off() → gc().

### 7.3 Figure Conventions

- **PNG authoritative**. PDF optional with tryCatch+cairo_pdf and explicit file paths.
- **Color scheme**: CA3 = orange (`#E69F00` or similar) — consistent with Phase 11/12.
- **All CA3 figures**: subtitle must include `n_CA3 = 11` (sample count).
- **F07b subtitle**: `n_CA3_Young = 3 (CAUTION: severely underpowered)`.
- **All correlation figures**: display n prominently.
- **No Rplots.pdf**: all device calls use explicit file paths; `while (dev.cur() > 1) invisible(dev.off())` after each tryCatch block.

---

## 8. Package Policy

**NO new packages**. Phase 13b uses only packages already in the Phase 13 environment:

| Package | Usage |
|---------|-------|
| dplyr | Data manipulation, filtering, sorting |
| ggplot2 | All new figures |
| pheatmap | Heatmaps (F01b, F04b) |
| patchwork | Multi-panel composition |
| stringr | String operations |
| Seurat (optional, E4 only) | `SpatialFeaturePlot` |
| digest | renv.lock MD5 check |

**renv.lock**: Current MD5 `b2bf05038b440576855d043e137c5883`. Must remain unchanged throughout Phase 13b.

---

## 9. Relationship to Phase 13

### 9.1 What is NOT Changed

- `phase13_regulator_score_default.csv` — untouched
- `phase13_regulator_score_sensitivity.csv` — untouched
- `phase13_candidate_classes.csv` — untouched
- `phase13_spearman_rho_all.csv` and `_universe_regulator.csv` — untouched
- `phase13_self_correlation_flags.csv` — untouched
- `phase13_target_gene_module_association.csv` — untouched
- `phase13_tf_annotation_coverage.csv` — untouched
- `phase13_final_report.txt` — untouched (Phase 13b may create a separate addendum file)
- `phase13_validation_summary.csv` — untouched (Phase 13b creates its own validation file)
- All existing figures — untouched (but renamed in guide for clarity)

### 9.2 What is ADDED

- 4 top-candidate CSVs (O1-O4)
- ~15 CA3 mirror figures (Group B in §3.2)
- 1 Phase 13b validation CSV (O5)
- Guide update (D1)
- Optionally: CA3 GO enrichment figure (if clusterProfiler available and per-module gene count ≥ 5)
- Optionally: report addendum file (O6)

---

## 10. Plan Compliance Report

### 10.1 File Creation

| Check | Status |
|-------|--------|
| Only created `docs/spatial_phase13b_ca3_symmetric_display_plan.md` | YES |
| No R scripts created | YES |
| No R analysis executed | YES |
| No RDS loaded | YES |
| No packages installed | YES |
| No renv.lock modified | YES (MD5: `b2bf05038b440576855d043e137c5883`) |
| No data/figures generated | YES |

### 10.2 Files Inspected

| # | File | Reviewed |
|---|------|----------|
| 1 | `README.md` | Reviewed — confirmed spatial branch conventions, output paths, memory rules |
| 2 | `AGENTS.md` | Reviewed — confirmed Phase 12 guardrails (remains reference), script naming convention (s##_), Phase 13b is spatial supplement, CA3 Young n=3 caveat, Complex V caveat |
| 3 | `docs/repository_structure.md` | Reviewed — confirmed output directory policy (data/processed/spatial/, figures/spatial/) |
| 4 | `docs/spatial_phase13_candidate_regulator_discovery_plan.md` | Reviewed (984 lines) — confirmed Phase 13 design (3 universes, composite score formula, self-correlation Strategy A, 12 FATAL stop conditions, validation schema) |
| 5 | `docs/spatial_phase13_candidate_regulator_analysis_guide.md` | Reviewed (102 lines) — identified CA1-bias in guide language, confirmed guide structure for update planning |
| 6 | `R/spatial/s13_candidate_regulator_discovery.R` | Reviewed (1470 lines) — identified 11 specific CA1-bias locations (F01 line 614, F03 line 696, F04 line 729, F07 line 806, F08 line 837, F09 line 862, F10 line 891, F13 line 958, F14 line 994, E12 line 1032, guide line 1108) |
| 7 | `data/.../phase13_regulator_score_default.csv` | Reviewed header — confirmed 12 columns including `score_CA3_default`, `rank_CA3` (13,178 data rows) |
| 8 | `data/.../phase13_candidate_classes.csv` | Reviewed — confirmed 6 classes: low_confidence (4608), CA3_aging (3290), CA1_aging (3027), concordant (1573), region_flipped (674), target_readout (6) |
| 9 | `data/.../phase13_final_report.txt` | Reviewed — confirmed 11 modules scored, 2 modules excluded (Complex V + mtDNA-encoded), 12/0/0 validation pass |
| 10 | `data/.../phase13_validation_summary.csv` | Reviewed — 12 PASS / 0 WARN / 0 FAIL, renv.lock unchanged confirmed |

### 10.3 Key Design Decisions

| Decision | Rationale |
|----------|----------|
| No re-computation of Phase 13 core logic | Phase 13b is display supplement only. All rho/score/rank values already exist for both regions. |
| Output back into `phase13_candidate_regulator/` directory | Per task instructions. Phase 13b is a supplement, not a new phase requiring separate directories. |
| CA3 Young n=3 warning on F07b | Must be explicit. CA3 Young is severely underpowered and Age/GEMgroup are colinear. |
| `combined_score = score_CA1 + score_CA3` for concordant/flipped | Simple, interpretable, and does not create new weighting. Alternative `pmax()` was considered but `sum` better captures the "strong in both regions" intent. `pmax` kept as secondary sort tiebreaker. |
| Keep existing filenames; add `b` suffix for CA3 | Safer than renaming existing gitignored files. The analysis guide maps old names to their CA1 focus. |
| Guide update appends, does not overwrite core content | Phase 13's guide content (caveats, metric hierarchy, limitations) is correct and should not be replaced. |
| Spatial maps conditional (not blocking) | Loading 1.3GB Hippo RDS for 5 CA3 spatial plots is optional. Tabular corrections are the primary deliverable. |

---

## 11. Execution Safety Assessment

Phase 13b is **safe to execute** after user approval because:

1. **Read-only on Phase 13 core data**: All inputs are CSVs; no RDS is modified or regenerated.
2. **No new packages**: Uses only existing packages.
3. **No renv.lock change**: No `renv::install()` or package additions.
4. **Additive outputs**: New files only; no existing files deleted or overwritten (except guide, which is text).
5. **Fatal stop conditions**: Script stops before any output if critical inputs are missing.
6. **Memory guard**: Hippo RDS loading is optional, guarded by tryCatch with WARN fallback.
7. **Rplots.pdf guard**: Explicit file paths for all devices, cleanup at end.

---

## 12. References

- Phase 13 Plan v1.2: `docs/spatial_phase13_candidate_regulator_discovery_plan.md`
- Phase 13 Script: `R/spatial/s13_candidate_regulator_discovery.R`
- Phase 13 Results: `data/processed/spatial/phase13_candidate_regulator/`
- Phase 11/12 guardrails: `AGENTS.md` §"Phase 12 Guardrails"
- Repository structure: `docs/repository_structure.md`
