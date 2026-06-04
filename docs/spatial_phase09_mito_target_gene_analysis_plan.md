# Phase Spatial-09: CA1/CA3/DG Target-Gene Spatial Analysis Planning

> **Plan version**: v2.1 (2026-06-04)
> **Status**: APPROVED — Phase 09 execution
> **Git HEAD at plan time**: 6ca46a2
> **Handoff**: GO_WITH_CAVEATS | Tangram: TANGRAM_DEFERRED

---

## 1. Read-Only Audit Summary

### 1.1 Input Object — Confirmed

| Field | Value |
|-------|-------|
| Object | `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` |
| Spots | 9,921 (Phase 06 validated, 0.0% off) |
| Images | 16 (slice1 through slice1_16) |
| Assays | Spatial (32,285 features), SCT (17,096 features) |
| Region levels | CA1, CA2, CA3, ML, GCL, Hilus |
| Region counts | CA1=4671, CA3=2321, ML=1480, GCL=645, Hilus=239 (DG=2364) |
| Age groups | Young (N=4 GEMgroups), Middle (N=4), Old (N=8) |
| Metadata | Region, Age, GEMgroup columns present and validated |

### 1.2 `docs/target_genes.xlsx` — Status

| Property | Value |
|----------|-------|
| Path | `docs/target_genes.xlsx` |
| Size | 12,036 bytes |
| Last modified | 2026-06-04 09:58 |
| `readxl` in renv | **NO** |
| `openxlsx` in renv | **NO** |
| Sheets / columns / rows | **Unknown** |
| Human vs mouse symbols | **Unknown** |
| Category/source columns | **Unknown** |
| Status | **Exists but not inspected** |

### 1.3 Caveat Boundary — Confirmed

Phase 04/05 discrepancies (D01-D09) affect IFNγ/Edge/STutility axis, NOT Region axis. CA1/CA3/DG labels are from author Script 6 RDS metadata, validated by Phase 06 (6 regions at 0.0% off). This boundary is clean.

### 1.4 ref-bio Routing

| Priority | Source ID | Upstream URL |
|----------|-----------|-------------|
| must_read | `seurat-spatial` | https://satijalab.org/seurat/articles/spatial_vignette.html |
| must_read | `seurat` | https://satijalab.org/seurat/ |
| optional_read | `pseudobulk-de-guidance` | (upstream DESeq2/edgeR vignettes) |

*Note: OSTA, SpatialExperiment, nnSVG, SPARK-X are matched in the index but not directly used by Phase 09/10. They remain available as optional reads if spatial statistics become relevant in later phases.*

---

## 2. Region Definitions

### 2.1 Primary Regions

| Region Group | Definition | Spot Count | Role |
|-------------|-----------|------------|------|
| **CA1** | `Region == "CA1"` | 4,671 | Primary question region |
| **CA3** | `Region == "CA3"` | 2,321 | Primary question region |
| **DG** | `Region %in% c("ML", "GCL", "Hilus")` | 2,364 | Primary question region (combined) |
| **CA2** | `Region == "CA2"` | 565 | Context/QC only; not a primary comparison target |

### 2.2 Subregion Preservation

DG subregions (ML=1480, GCL=645, Hilus=239) are preserved as separate metadata values AND included in subregion-level output tables alongside the combined DG grouping. Phase 10 outputs both grouped and subregion-level summaries.

### 2.3 CA2

CA2 is included in audit tables for completeness and as a sanity-check region (e.g., does a gene show uniform expression across all hippocampal regions?). CA2 is NOT a primary comparison target. If the user later requests CA2-specific analysis, that requires a separate plan update.

---

## 3. Input Object Strategy

| Object | Phase 09 | Phase 10 |
|--------|----------|----------|
| `seurat_Visium_Hippo_All.rds` | Document structure | **Load as sole spatial object** |
| `seurat_Visium_DG_All.rds` | Note existence only | Do NOT load (DG defined from Hippo metadata) |
| `seurat_Visium_WholeBrain.rds` | Not referenced | Do NOT load (BLOCKED_WHOLEBRAIN) |
| Phase 04/05 IFNγ Edge/NNS/ENS labels | Not referenced | Explicitly excluded from Phase 09/10 region analysis; recorded as orthogonal functional labels, not region labels |
| Phase 04 `metadata_hippo_func_region.csv` | Not referenced | Do NOT use |

**Memory rule**: Load Hippo only. Release via `rm()` + `gc()` after each analysis stage.

---

## 4. Gene List Strategy

### 4.1 Tiered Reading Strategy for `target_genes.xlsx`

```
Tier 1 (PREFERRED — no new R dependency):
  User exports docs/target_genes.xlsx -> docs/target_genes.csv
  Phase 10 reads CSV via read.csv() (base R)

Tier 2 (OPTIONAL — user-approved dependency):
  User explicitly approves readxl or openxlsx installation
  Phase 10 installs into renv, records renv.lock change
  Phase 10 reads .xlsx directly

Tier 3 (FALLBACK — structure-only inspection without new package):
  Use system `file` / `unzip -l` to inspect ZIP structure of .xlsx
  This reveals internal XML file list (e.g., sheet1.xml, sharedStrings.xml)
  Does NOT reliably reveal sheet names, column headers, row counts,
  gene symbols, or content dimensions
  Tier 3 is for confirming the file is a valid .xlsx only
```

**Phase 09 stance**: Record that `target_genes.xlsx` exists but is uninspectable. Do NOT infer sheet names, column names, row counts, or gene symbols from the binary file. Do NOT install packages.

### 4.2 Category/Source Attribution

**Primary rule**: Gene category and source must come from columns within `target_genes.xlsx` itself. If the workbook has no `category` or `source` column, Phase 10 marks `category_status = missing_metadata` for all genes.

**Secondary annotation** (Phase 10 or later, requires separate approval):
- MitoCarta 3.0 mouse gene lists
- GO:0005739 (mitochondrion) gene sets
- Existing `gene_lists/` ferroptosis gene lists in the project

These are **supplemental annotations** only. They do NOT auto-expand the gene list, do NOT auto-reclassify user-provided categories, and do NOT add or remove genes. A gene can have both `category_user` (from xlsx) and `annotation_mitocarta` (supplemental) columns independently.

### 4.3 Gene Symbol Mapping — Safe Audit

Phase 10 must produce a gene symbol mapping audit with these columns:

| Column | Source/Logic |
|--------|-------------|
| `original_symbol` | As read from target_genes.xlsx (verbatim) |
| `original_case` | UPPERCASE / lowercase / Mixed / TitleCase |
| `exact_match_spatial` | Does `original_symbol` exist verbatim in `rownames(hippo[["Spatial"]])`? TRUE/FALSE |
| `exact_match_sct` | Does `original_symbol` exist verbatim in `rownames(hippo[["SCT"]])`? TRUE/FALSE |
| `candidate_mouse_symbol` | If no exact match and `original_case` is UPPERCASE → derive via base R helper (see below). If already TitleCase or mixed → set to `original_symbol` |
| `candidate_match_spatial` | Does `candidate_mouse_symbol` exist in `rownames(hippo[["Spatial"]])`? TRUE/FALSE |
| `candidate_match_sct` | Does `candidate_mouse_symbol` exist in `rownames(hippo[["SCT"]])`? TRUE/FALSE |
| `mapping_status` | `exact_match` / `candidate_match` / `not_found_in_either` / `ambiguous_multiple_matches` |
| `manual_review` | TRUE if `not_found_in_either` or `ambiguous_multiple_matches`; FALSE otherwise |
| `final_symbol_for_analysis` | `NA` until user approves mapping; filled after review |

**Base R helper for candidate mouse symbol generation** (no `stringr` dependency):

```r
# Converts "SLC7A11" -> "Slc7a11", "TP53" -> "Tp53"
# Only applied when original_case is UPPERCASE and no exact match found
to_title_case_base <- function(x) {
  x_lower <- tolower(x)
  substr(x_lower, 1, 1) <- toupper(substr(x_lower, 1, 1))
  return(x_lower)
}
```

**Rules**:
1. Do NOT silently apply `tolower()` + capitalize as a global conversion
2. Always prefer `exact_match` over `candidate_match`
3. Never rename a gene in the output without recording the mapping
4. All `manual_review = TRUE` genes must be resolved before Phase 11 (statistical inference)
5. Human-to-mouse ortholog mapping (if needed for genes not present under either symbol) requires separate user approval and is NOT part of Phase 10

### 4.4 Prohibited Gene-List Operations

- Do NOT auto-add genes (e.g., "all mt- genes") to the user's list
- Do NOT auto-remove genes flagged as missing
- Do NOT auto-rename genes without audit trail
- Do NOT copy `target_genes.xlsx` to `data/` or generate derivative gene-list files without explicit approval
- Do NOT use ferroptosis gene lists from `gene_lists/` to expand the target list unless the user requests it

---

## 5. Assay / Layer Strategy

| Use Case | Assay / Slot | Rationale |
|----------|-------------|-----------|
| Gene presence audit | `rownames(hippo[["Spatial"]])` | Covers all 32,285 features |
| Expression extraction | `GetAssayData(hippo, assay = "Spatial", layer = "counts")` | Raw UMI counts, sparse |
| Visualization (spatial feature plots) | `Spatial` counts or SCT `data` | SCT `data` is author-normalized; note source in figure caption |
| Summary statistics (mean, median, % expressed) | `Spatial` `counts` | Raw counts |
| Pseudobulk aggregation (Phase 11+) | `Spatial` `counts` → `aggregate.Matrix()` | Required for DESeq2 |

**Prohibited**:
- Do NOT use SCT `data` for count-based inference
- Do NOT rerun `SCTransform()`
- Do NOT coerce sparse count matrices to dense

---

## 6. Replicate / Statistical Design (for Phase 10 & forward)

### 6.1 Core Principle

Visium spots are not independent biological replicates. All per-spot expression is summarized to sample-level (GEMgroup × Region × Age) for any inferential step.

### 6.2 Sample Structure

| Age | GEMgroups | N samples |
|-----|-----------|-----------|
| Young | 1, 2, 3, 4 | 4 |
| Middle | 5, 6, 7, 8 | 4 |
| Old | 9, 10, 11, 12, 13, 14, 15, 16 | 8 |

**Imbalance**: Old has 2× the samples of Young/Middle. This is a design constraint inherited from the author's experimental design and must be documented in all output.

### 6.3 Age Group Strategy

- Phase 10 outputs summaries for all three Age groups (Young, Middle, Old)
- Middle is preserved in all tables
- If the user narrows to Old vs Young for Phase 11+ statistical tests, Middle remains available as a reference
- The imbalance in sample counts (N=4 vs N=8) is noted but not corrected in Phase 10

### 6.4 Spot-Count Limitations

Some Region × GEMgroup combinations have very small spot counts:
- CA2 Young GEMgroup 2: 31 spots; CA2 Young GEMgroup 4: 50 spots
- Hilus Old GEMgroup 9: 6 spots; Hilus Middle GEMgroup 5: 8 spots

Phase 10 must tabulate per-GEMgroup spot counts and flag combinations with < 20 spots as `low_coverage`. Phase 11+ pseudobulk aggregation must address whether to exclude or pool low-coverage GEMgroups.

---

## 7. Phase 09 / Phase 10 Boundary

### 7.1 Phase 09 Scope (Current)

| Activity | Allowed? |
|----------|----------|
| Write plan documents to `docs/` | YES — `docs/spatial_phase09_mito_target_gene_analysis_plan.md` + `docs/spatial_mito_target_gene_strategy.md` |
| Create R scripts | NO |
| Load RDS objects | NO |
| Install packages | NO |
| Modify renv.lock | NO |
| Generate data/figures | NO |
| Read target_genes.xlsx | NO (structure-only audit in plan text) |

### 7.2 Phase 10 Scope (Proposed — Execution Gate)

| Activity | Allowed? |
|----------|----------|
| Load `seurat_Visium_Hippo_All.rds` | YES |
| Read gene list (CSV preferred; xlsx if user-approved) | YES |
| Audit gene presence in Spatial/SCT | YES |
| Generate gene symbol mapping audit | YES |
| Compute region-level expression summaries | YES |
| Generate spatial feature plots for selected genes | YES (selected only) |
| Run DESeq2 / statistical tests | NO |
| Perform biological interpretation | NO |
| Run enrichment (GO/KEGG) | NO |
| Generate heatmaps / volcano plots | NO |
| Save Seurat objects | NO |
| Install packages without approval | NO |

---

## 8. Phase 10 Output Design (Revised)

### 8.1 Always Produced

| File | Path | Description |
|------|------|-------------|
| Script | `R/spatial/s10_target_gene_audit_region_summary.R` | Phase 10 execution script |
| Gene presence audit | `data/processed/spatial/phase10_target_gene_audit/gene_presence_audit.csv` | Per-gene, per-assay presence |
| Symbol mapping audit | `data/processed/spatial/phase10_target_gene_audit/gene_symbol_mapping_audit.csv` | Columns per §4.3 |
| Category mapping | `data/processed/spatial/phase10_target_gene_audit/gene_category_mapping.csv` | From xlsx columns; `missing_metadata` if absent |
| Region × Sample summary | `data/processed/spatial/phase10_target_gene_audit/region_sample_summary.csv` | Per gene × region_group × Age × GEMgroup: n_spots, pct_expressed, mean_count, sum_count |
| Region × Age summary | `data/processed/spatial/phase10_target_gene_audit/region_age_summary.csv` | Per gene × region_group × Age: collapsed from sample-level |
| Subregion summary | `data/processed/spatial/phase10_target_gene_audit/region_subregion_summary.csv` | Per gene × subregion (ML/GCL/Hilus) × Age |
| Low-coverage flags | `data/processed/spatial/phase10_target_gene_audit/low_coverage_gemgroups.csv` | Region × GEMgroup combinations with < 20 spots |
| Validation | `data/processed/spatial/phase10_target_gene_audit/phase10_validation_summary.txt` | Check results |

### 8.2 Conditionally Produced

| File | Condition | Description |
|------|-----------|-------------|
| Missing genes list | Only if genes not found | `missing_genes.csv` — genes absent in both assays |
| Manual review queue | Only if ambiguous symbols | `manual_review_genes.csv` — genes requiring user decision |
| Spot-level expression (long) | Only if total target genes ≤ 30 OR user approves | `spot_expression_long.csv` — one row per spot×gene with count, region, age |
| Spatial feature plots | Selected genes only by default | `figures/spatial/phase10_target_gene_audit/` — up to 12 plots (e.g., top 4 genes × 3 categories) |
| All-gene spatial plots | Only if user explicitly requests | May generate many PNGs; explicit user confirmation required |

### 8.3 Plot Selection Logic (Default)

For spatial feature plots, select up to **4 genes per category** (if categories exist) or **6 genes total** (if no categories), prioritized by:
1. Genes with highest variance across regions
2. Mix of detected and borderline-detected genes
3. User-flagged genes if any notes indicate priority

---

## 9. Validation Checks (Phase 10)

| ID | Check | Criterion |
|----|-------|-----------|
| V01 | Gene list readable | CSV loaded with ≥ 1 gene row, or xlsx with user approval |
| V02 | Gene column identified | Exactly 1 column mapped to gene symbols (or user-specified column name) |
| V03 | Gene presence audit complete | Every gene from list checked against `rownames(hippo[["Spatial"]])` |
| V04 | Symbol mapping audit complete | Every gene has mapping_status and manual_review flag |
| V05 | Region counts preserved | CA1=4671, CA3=2321, DG=2364 |
| V06 | No dense conversion | All count extractions preserve sparse matrices |
| V07 | No Seurat object saved | No `.rds` file written |
| V08 | GEMgroup × Age mapping verified | 16 GEMgroups, Age mapping matches Phase 06 |
| V09 | Missing genes documented | Count and names recorded in `missing_genes.csv` |
| V10 | All found genes summarized | Every found gene present in `region_sample_summary.csv` |
| V11 | Output files are lightweight CSVs | No file > 50 MB |
| V12 | renv.lock unchanged | MD5 verified (unless readxl user-approved) |
| V13 | No WholeBrain loaded | Provenance documents single-object load |
| V14 | No CA2 auto-excluded | CA2 present in audit tables as context |

---

## 10. Proposed Phase 10 Outline (Script Steps)

```
Phase Spatial-10: Target Gene Audit + Region-Aware Expression Summary

E0: Setup
  - Verify gene list source (CSV preferred; xlsx fallback with user approval)
  - Load seurat_Visium_Hippo_All.rds
  - Verify invariants: 9921 spots, 6 regions, 16 images
  - Rplots.pdf guard start

E1: Read Gene List
  - Read CSV or xlsx
  - Record: source file, sheet name (if xlsx), column names, row count
  - Identify gene symbol column (by name heuristics: "gene", "symbol", "Gene", "Symbol";
    if ambiguous, stop and request user specification)
  - Extract category/source columns if present; otherwise fill as "missing_metadata"

E2: Gene Symbol Mapping Audit
  - For each gene symbol:
    * Exact match in Spatial rownames
    * Exact match in SCT rownames
    * If uppercase → derive candidate_mouse_symbol via to_title_case_base()
    * Candidate match in Spatial/SCT rownames
  - Assign mapping_status
  - Flag genes for manual_review
  - Output gene_symbol_mapping_audit.csv

E3: Gene Presence Audit
  - Using final_symbol_for_analysis where resolved, original_symbol otherwise
  - For each gene: check presence in Spatial counts, SCT data
  - Output gene_presence_audit.csv
  - Output missing_genes.csv (if any)

E4: Spot-Level Expression Extraction
  - For genes found in Spatial assay:
    * Extract raw counts via GetAssayData(assay = "Spatial", layer = "counts")
    * Preserve sparse matrix format
    * Build metadata-aligned table
    * Only for resolved genes with mapping_status in {exact_match, candidate_match}
  - Output spot_expression_long.csv (only if gene count ≤ 30 OR user approves)

E5: Region-Level Summaries
  - Per gene × region_group (CA1, CA3, DG, CA2) × Age × GEMgroup:
    * n_spots, pct_expressed, mean_count, median_count, sum_count
  - Collapse to region_group × Age
  - Also produce subregion-level (ML, GCL, Hilus) summaries
  - Flag low-coverage GEMgroups
  - Output: region_sample_summary.csv, region_age_summary.csv,
    region_subregion_summary.csv, low_coverage_gemgroups.csv

E6: Spatial Feature Plots (Selected Genes)
  - Select up to 4 genes per category or 6 total
  - Generate SpatialFeaturePlot per image (16 images per gene)
  - Save as PNG in figures/spatial/phase10_target_gene_audit/

E7: Validation
  - Run checks V01-V14
  - Write phase10_validation_summary.txt

E8: Cleanup
  - rm(hippo); gc()
  - Rplots.pdf guard end
  - Verify renv.lock MD5
```

---

## 11. Stop Conditions

| ID | Condition | Action |
|----|-----------|--------|
| SC-1 | `seurat_Visium_Hippo_All.rds` not found | STOP — object unavailable |
| SC-2 | Gene list source unreadable (no CSV, no xlsx reader, user declines readxl) | STOP — request CSV export |
| SC-3 | Gene list has 0 gene rows after parsing | STOP — empty list |
| SC-4 | No gene symbol column identifiable | STOP — request user specification |
| SC-5 | Phase 06 validation summary missing or non-PASS | STOP — region framework not verified |
| SC-6 | Phase 07 handoff is not GO_WITH_CAVEATS | STOP — recheck handoff state |
| SC-7 | `readxl` installation proposed but user declines | STOP — fall back to CSV path |
| SC-8 | User requests WholeBrain loading during Phase 10 | STOP — requires separate plan |
| SC-9 | User requests Tangram during Phase 10 | STOP — requires separate plan |
| SC-10 | User requests DESeq2 / statistical inference during Phase 10 | STOP — Phase 11+ scope |
| SC-11 | User requests enrichment during Phase 10 | STOP — Phase 11+ scope |
| SC-12 | User requests biological interpretation before Phase 10 audit complete | STOP — premature |

---

## 12. Acceptance Criteria (Phase 09 Plan)

| ID | Criterion |
|----|-----------|
| AC1 | Plan document saved as `docs/spatial_phase09_mito_target_gene_analysis_plan.md` |
| AC2 | Gene source strategy documented in `docs/spatial_mito_target_gene_strategy.md` |
| AC3 | Tiered reading strategy for target_genes.xlsx explicitly stated |
| AC4 | No auto-categorization of genes; defer to workbook columns |
| AC5 | Safe symbol mapping audit specified (exact match -> candidate -> manual review) |
| AC6 | Phase 09/10 boundary explicitly documented |
| AC7 | Phase 10 output defaults are lightweight |
| AC8 | CA2 defined as context/QC, not primary |
| AC9 | DG subregion preservation specified |
| AC10 | All stop conditions listed |
| AC11 | Current plan-mode writes no files; approved Phase 09 execution writes only the two docs listed in Section 14 |
| AC12 | All caveats from D01-D07 acknowledged |
| AC13 | ref-bio source routing report included |

---

## 13. Open Questions for User Decision

| # | Question | Context |
|---|----------|---------|
| Q1 | **CSV export**: Are you willing to export `docs/target_genes.xlsx` as `docs/target_genes.csv` to avoid needing `readxl`? | Preferred path; no new R dependency |
| Q2 | **readxl approval**: If CSV is not feasible, do you approve installing `readxl` (and possibly its system dependencies) into renv for Phase 10? | Will be recorded in renv.lock |
| Q3 | **Category columns**: Does `target_genes.xlsx` have `category` and/or `source` columns, or is it a flat gene list? | Determines whether Phase 10 uses workbook metadata or marks `missing_metadata` |
| Q4 | **MitoCarta / GO annotation**: Should Phase 10 add supplemental annotation columns from MitoCarta 3.0 and/or GO:0005739 as separate columns (NOT replacing user categories)? | Requires separate annotation step |
| Q5 | **Spot-level table**: Do you want `spot_expression_long.csv` for all target genes (~1 row per spot per gene), only if gene count is small, or skip entirely? | ~10K spots × N genes; can be large for N > 30 |
| Q6 | **Spatial plots scope**: All genes or selected genes only? | All genes could generate 16 × N PNGs |
| Q7 | **Age contrast**: Should the primary design focus on Old vs Young (with Middle preserved but not contrasted)? | Matches scRNA-seq pipeline convention |
| Q8 | **Manual review**: For genes that cannot be auto-mapped to mouse symbols, do you want to review them before Phase 10 proceeds, or should Phase 10 flag them and continue with resolved genes only? | Affects Phase 10 completeness |
| Q9 | **Phase 09 documents**: Should I create `docs/spatial_phase09_mito_target_gene_analysis_plan.md` and `docs/spatial_mito_target_gene_strategy.md` now, or wait for your review of this plan first? | This is the only Phase 09 file-creation action |

---

## 14. Proposed File Creation (Phase 09 Execution)

When the user approves this plan, I will create exactly these two files:

1. **`docs/spatial_phase09_mito_target_gene_analysis_plan.md`** — This full plan document
2. **`docs/spatial_mito_target_gene_strategy.md`** — Condensed gene-list strategy reference:
   - Tiered reading strategy
   - Symbol mapping audit specification
   - Category/source attribution rules
   - Prohibited operations list
   - MitoCarta/GO/future-annotation gating rules

No other files. No R scripts. No data. No figures. No package installs. No renv.lock changes.
