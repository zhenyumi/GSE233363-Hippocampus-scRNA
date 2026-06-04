# Spatial Target-Gene Strategy Reference

> **Version**: 2.1 (2026-06-04)
> **Phase**: Spatial-09 (planning) → Spatial-10 (execution)
> **Scope**: CA1/CA3/DG target-gene analysis using Visium spatial transcriptomics

---

## 1. Tiered Reading Strategy for `target_genes.xlsx`

### Tier 1 — CSV Export (PREFERRED, no new R dependency)

```
User action: Export docs/target_genes.xlsx -> docs/target_genes.csv
R code:      gene_df <- read.csv("docs/target_genes.csv", stringsAsFactors = FALSE)
Advantage:   Base R only; no renv.lock change; no package install
```

### Tier 2 — Direct xlsx Read (OPTIONAL, requires user approval)

```
User approval: Explicitly approve readxl or openxlsx installation
R code:        install.packages("readxl"); gene_df <- readxl::read_excel("docs/target_genes.xlsx")
Renv impact:   renv.lock will be modified; record change
Phase 10 only: Do not install during Phase 09
```

### Tier 3 — Structure-Only Inspection (FALLBACK, no content extraction)

```
System command: unzip -l docs/target_genes.xlsx
What it reveals: Internal ZIP/XML file list (e.g., xl/worksheets/sheet1.xml, xl/sharedStrings.xml)
What it does NOT reveal:
  - Sheet names (require parsing XML)
  - Column headers (require parsing sharedStrings.xml)
  - Row counts (require parsing sheet XML)
  - Gene symbols (require full XML parse)
  - Content dimensions (not reliably extractable)
Purpose: Confirm file is valid .xlsx only; NOT a content inspection method
```

**Phase 09 stance**: Record that `target_genes.xlsx` exists but is uninspectable. Do NOT infer sheet names, column names, row counts, or gene symbols. Do NOT install packages.

---

## 2. Gene Symbol Mapping Audit — Safe Audit Specification

### 2.1 Audit Columns

Phase 10 produces `gene_symbol_mapping_audit.csv` with these columns:

| Column | Type | Description |
|--------|------|-------------|
| `original_symbol` | character | Gene symbol as read from workbook (verbatim, preserved) |
| `original_case` | character | UPPERCASE / lowercase / Mixed / TitleCase |
| `exact_match_spatial` | logical | TRUE if `original_symbol` in `rownames(hippo[["Spatial"]])` |
| `exact_match_sct` | logical | TRUE if `original_symbol` in `rownames(hippo[["SCT"]])` |
| `candidate_mouse_symbol` | character | If no exact match and case is UPPERCASE: apply `to_title_case_base()`. Otherwise: `original_symbol` |
| `candidate_match_spatial` | logical | TRUE if `candidate_mouse_symbol` in `rownames(hippo[["Spatial"]])` |
| `candidate_match_sct` | logical | TRUE if `candidate_mouse_symbol` in `rownames(hippo[["SCT"]])` |
| `mapping_status` | character | `exact_match` / `candidate_match` / `not_found_in_either` / `ambiguous_multiple_matches` |
| `manual_review` | logical | TRUE if `not_found_in_either` or `ambiguous_multiple_matches` |
| `final_symbol_for_analysis` | character | NA until user approves; filled after review |

### 2.2 Base R Helper for Candidate Mouse Symbol

```r
# No stringr dependency. Base R only.
# Converts "SLC7A11" -> "Slc7a11", "TP53" -> "Tp53"
# Only applied when original_case is UPPERCASE and no exact match found
to_title_case_base <- function(x) {
  x_lower <- tolower(x)
  substr(x_lower, 1, 1) <- toupper(substr(x_lower, 1, 1))
  return(x_lower)
}
```

### 2.3 Mapping Rules

1. **Exact match first**: Always prefer verbatim match over candidate
2. **No silent global conversion**: Do NOT apply `tolower()` + capitalize to all genes
3. **Audit trail**: Every gene gets a `mapping_status`; no gene silently dropped
4. **Manual review queue**: Genes with `not_found_in_either` or `ambiguous_multiple_matches` go to `manual_review = TRUE`
5. **No auto-resolution**: `final_symbol_for_analysis` remains NA until user approves
6. **Ortholog mapping separate**: If a gene needs human-to-mouse ortholog lookup (e.g., BioMART), that requires separate user approval and is NOT part of Phase 10

---

## 3. Category / Source Attribution Rules

### 3.1 Primary Rule

Gene category and source come from columns within `target_genes.xlsx` itself.

| Scenario | Action |
|----------|--------|
| Workbook has `category` column | Use it as `category_user` |
| Workbook has `source` column | Use it as `source_user` |
| Workbook has both | Use both independently |
| Workbook has neither | Fill `category_status = "missing_metadata"` for all genes |

### 3.2 Secondary Annotation (Future — Requires Approval)

These annotations are **supplemental only** and do NOT replace or override user categories:

| Annotation Source | Column Name | Notes |
|------------------|-------------|-------|
| MitoCarta 3.0 (mouse) | `annotation_mitocarta` | Boolean or category string |
| GO:0005739 (mitochondrion) | `annotation_go_mito` | Boolean or GO term |
| Ferroptosis gene lists (`gene_lists/`) | `annotation_ferroptosis` | Boolean |

**Gating rule**: MitoCarta/GO/future annotation is NOT applied automatically. It requires:
1. User approval to add annotation step
2. Separate script or code block (not embedded in main audit logic)
3. Annotation columns added as new columns, not replacing existing ones
4. No gene expansion — only annotating genes already in the user's list

### 3.3 Prohibited Operations

- Do NOT auto-add genes from MitoCarta/GO to the user's target list
- Do NOT auto-remove genes flagged as missing from Spatial/SCT
- Do NOT auto-rename genes without audit trail entry
- Do NOT copy `target_genes.xlsx` to `data/` or generate derivative gene-list files
- Do NOT use `gene_lists/` ferroptosis lists to expand target list unless user requests
- Do NOT merge multiple gene lists without explicit user approval

---

## 4. Phase 09 / Phase 10 Boundary

### Phase 09 (Current — Plan-Only)

- Creates exactly 2 files in `docs/`
- No R scripts, no RDS loads, no package installs, no renv.lock changes
- No data/figures generation

### Phase 10 (Next — Execution Gate)

- Creates `R/spatial/s10_target_gene_audit_region_summary.R`
- Loads `seurat_Visium_Hippo_All.rds` as sole spatial object
- Reads gene list (CSV or user-approved xlsx)
- Performs gene audit + region-aware expression summary
- Outputs CSVs to `data/processed/spatial/phase10_target_gene_audit/`
- Selected spatial feature plots to `figures/spatial/phase10_target_gene_audit/`
- No DESeq2, no enrichment, no biological interpretation

---

## 5. Region Definitions Summary

| Region Group | Definition | Spots | Role |
|-------------|-----------|-------|------|
| CA1 | `Region == "CA1"` | 4,671 | Primary |
| CA3 | `Region == "CA3"` | 2,321 | Primary |
| DG | `Region %in% c("ML", "GCL", "Hilus")` | 2,364 | Primary (combined) |
| CA2 | `Region == "CA2"` | 565 | Context/QC only |
| ML (subregion) | `Region == "ML"` | 1,480 | Subregion audit |
| GCL (subregion) | `Region == "GCL"` | 645 | Subregion audit |
| Hilus (subregion) | `Region == "Hilus"` | 239 | Subregion audit |

---

## 6. Assay / Layer Reference

| Use Case | Assay | Layer | Notes |
|----------|-------|-------|-------|
| Gene presence check | Spatial | rownames | 32,285 features |
| Count extraction | Spatial | counts | Raw UMI, sparse |
| Visualization | Spatial or SCT | counts or data | Note source in caption |
| Pseudobulk (Phase 11+) | Spatial | counts | `aggregate.Matrix()` |

---

## 7. Replicate Design Summary

| Age | GEMgroups | N samples | Notes |
|-----|-----------|-----------|-------|
| Young | 1-4 | 4 | Balanced with Middle |
| Middle | 5-8 | 4 | Preserved in all tables |
| Old | 9-16 | 8 | 2× imbalance vs Young/Middle |

Spots are NOT independent replicates. Summary to sample-level (GEMgroup × Region × Age) for any inferential step.

---

## 8. Open Questions Pending User Decision

| # | Question |
|---|----------|
| Q1 | CSV export of target_genes.xlsx? (preferred; avoids readxl) |
| Q2 | If not CSV: approve readxl installation? |
| Q3 | Does target_genes.xlsx have category/source columns? |
| Q4 | Add supplemental MitoCarta/GO annotation columns? |
| Q5 | Spot-level long table: all genes, small gene sets only, or skip? |
| Q6 | Spatial plots: selected genes only or all? |
| Q7 | Primary age contrast: Old vs Young? |
| Q8 | Manual review before Phase 10 proceeds, or flag-and-continue? |
| Q9 | Phase 09 docs creation timing: now or after user review? |

---

## References

- Phase 06: `docs/spatial_phase06_spot_region_assignment_validation.md`
- Phase 07: `docs/spatial_phase07_handoff.md`
- Seurat spatial vignette: https://satijalab.org/seurat/articles/spatial_vignette.html
- MitoCarta 3.0: https://www.broadinstitute.org/mitocarta
