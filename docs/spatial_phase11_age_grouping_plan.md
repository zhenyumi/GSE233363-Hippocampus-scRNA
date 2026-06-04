# Phase Spatial-11 Age Grouping: Age-Stratified CA1 vs CA3 DEG / Target-Gene / Module Analysis Plan

> **Plan version**: v1.1 (2026-06-04)
> **Status**: READ-ONLY PLAN — no analysis execution; only this plan file is created/modified; no analysis outputs generated
> **Handoff**: Phase 11 base PASS_WITH_WARNINGS (34 PASS, 1 WARNING) → Phase 11 Age Grouping planning
> **Phase 11 base**: 15 paired GEMgroups, 17,970 genes tested, 5,921 sig, 12 modules, 30 coupling pairs
> **bioskill**: NOT installed / NOT inspected — no content referenced or invented
> **Authority**: This saved `.md` file is the sole execution authority. Execution phase must re-read this file first.

---

## 0. ref-bio Pre-Answer Report

### Task Classification

| Dimension | Detected |
|-----------|----------|
| Data modality | **spatial transcriptomics** (Visium, RDS-based) |
| Workflow stage | **differential expression** (age-stratified) + **module score** + **coupling** + **cross-age comparison** |
| Package/tool | **DESeq2**, **Seurat**, **ggplot2**, **pheatmap**, **corrplot** (installed), **matrixStats** |
| Platform | **Visium** (RDS-based reproduction) |

### Matched Source IDs

| Priority | Source ID | Upstream URL |
|----------|-----------|-------------|
| must_read | `deseq2` | https://bioconductor.org/packages/release/bioc/html/DESeq2.html |
| must_read | `seurat-spatial` | https://satijalab.org/seurat/articles/spatial_vignette.html |
| must_read | `seurat` | https://satijalab.org/seurat/ |
| must_read | `pseudobulk-de-guidance` | upstream DESeq2 pseudobulk vignettes |
| optional_read | `limma` | https://bioconductor.org/packages/release/bioc/html/limma.html |

---

## 1. Scientific Goals (Explicitly Stated)

### 1.1 What Phase 11 Base Already Did

Phase 11 base (`s11_ca1_ca3_de_module_coupling.R`, v2.1) performed CA3 vs CA1 paired pseudobulk DESeq2 **across all ages merged**, with GEMgroup blocking (`~ GEMgroup + region_group`). No age stratification was performed on the DESeq2 step. Age appeared only in descriptive summaries and module delta tables.

### 1.2 What This Phase Answers

| Question | Unit | Method |
|----------|------|--------|
| **Q1: Young-specific CA3 vs CA1 DEG?** | 3 paired GEMgroups (G2,G3,G4), 6 pseudobulk samples | Age-stratified paired DESeq2: `~ GEMgroup + region_group` within Young only |
| **Q2: Middle-specific CA3 vs CA1 DEG?** | 4 paired GEMgroups (G5-G8), 8 pseudobulk samples | As above, within Middle only |
| **Q3: Old-specific CA3 vs CA1 DEG?** | 8 paired GEMgroups (G9-G16), 16 pseudobulk samples | As above, within Old only |
| **Q4: Do CA3-vs-CA1 target-gene effects change with age?** | 231 target genes × 3 age groups | Cross-age log2FC matrix; direction consistency check |
| **Q5: Which target gene categories/modules show age-dependent regional differences?** | 12 modules × 3 age groups | Module delta by age (descriptive first: n, mean_delta, median_delta, IQR per Age; Kruskal-Wallis only as optional/sensitivity) |
| **Q6: Does module coupling differ by age?** | 10 coupling pairs × 3 age groups | Age-stratified Spearman coupling (very exploratory at n=3,4) |

### 1.3 Explicit Constraints

- **Units are GEMgroup-level pseudobulk samples** (aggregated from Visium spots), not cells and not individual spots.
- **Spots are NOT biological replicates**. The replicate unit for all statistics is GEMgroup (tissue section).
- **CA3 vs CA1 contrast code**: positive log2FC = higher in CA3; negative log2FC = higher in CA1.
- **DG excluded from DE**: DG (= ML + GCL + Hilus) is context/secondary, not part of the CA3-vs-CA1 DE contrast.
- **CA2 included as context/QC**: present in summaries but not a primary comparison region.
- **Young and Middle are underpowered**: 3 and 4 paired GEMgroups respectively. All results from Young/Middle are labeled `exploratory_low_power`. Statistical significance in these groups should be interpreted with extreme caution.
- **Old (8 pairs) is the only group with adequate power** for DESeq2 DEG detection.

---

## 2. Analysis Routes: Comparison and Selection

### Route A: Age-Stratified Pseudobulk DESeq2 (PRIMARY)

**Method**: Run independent DESeq2 per age group, each using the paired design `~ GEMgroup + region_group`.

| Age | GEMgroups with CA1+CA3 | Pseudobulk samples | Paired design | Power label |
|-----|------------------------|-------------------|---------------|-------------|
| Young | G2, G3, G4 (3 pairs) | 6 (3 CA1 + 3 CA3) | `~ GEMgroup + region_group` | **exploratory_low_power** |
| Middle | G5, G6, G7, G8 (4 pairs) | 8 (4 CA1 + 4 CA3) | `~ GEMgroup + region_group` | **exploratory_low_power** |
| Old | G9, G10, G11, G12, G13, G14, G15, G16 (8 pairs) | 16 (8 CA1 + 8 CA3) | `~ GEMgroup + region_group` | **primary_adequately_powered** |

**GEMgroup inclusion rule per age**:
- GEMgroup 1: CA1 only in Young → excluded from Young paired DE (no CA3 mate). Included in Young two-group sensitivity if needed.
- All other GEMgroups have both CA1 and CA3 spots → included.

**Design validation per age** (MANDATORY before DESeq2):
1. Build `model.matrix(~ GEMgroup + region_group)` for that age only
2. Check `qr(model_matrix)$rank == ncol(model_matrix)` (full rank)
3. Count samples per GEMgroup × region_group cell
4. If full rank fails: fall back to `~ region_group` (two-group unpaired), labeled `fallback_unpaired`
5. Save `design_matrix_audit_by_age.csv`

**DESeq2 execution per age**:
```r
# Per age group, independently:
dds_age <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts_age,  # genes × samples for that age only
  colData   = sample_metadata_age,    # sample × (GEMgroup, region_group)
  design    = ~ GEMgroup + region_group
)
dds_age <- dds_age[rowSums(counts(dds_age)) >= 10, ]
dds_age <- DESeq(dds_age)
res_age <- results(dds_age, contrast = c("region_group", "CA3", "CA1"), alpha = 0.05)
```

**LFC shrinkage**: apeglm if available (installed in Phase 11 base), else "normal" fallback.

**Stop conditions**:
- SC-A1: `< 2 paired GEMgroups per age` → do NOT run paired DE for that age; mark as blocked/insufficient_pairs
- SC-A2: `pseudobulk count matrix all-zero` for that age → STOP that age only
- SC-A3: `DESeq2 fails to converge` for an age → STOP that age; continue other ages if independent and documented
- SC-A4: `design matrix rank-deficient` for an age → document; fall back to unpaired `~ region_group` with WARNING

**Output labeling requirement (MANDATORY per age)**:
Every age-specific DESeq2 output file, every figure, and the analysis guide MUST include:
- `model_type`: one of `paired_blocked` (primary design succeeded), `fallback_unpaired` (rank-deficient → `~ region_group` used), or `blocked_insufficient_pairs` (< 2 pairs, no DE run)
- `power_label`: `primary_adequately_powered` (Old only) or `exploratory_low_power` (Young, Middle)

**Fallback results must NOT be presented as interchangeable with primary paired results**. If an age uses `fallback_unpaired`, volcano titles, figure captions, and guide text must explicitly state "unpaired fallback model; GEMgroup pairing not preserved; interpret with caution".

**Why Route A is PRIMARY**:
- Directly answers Q1-Q3 with the same design as Phase 11 base, just subsetted by age.
- Preserves GEMgroup pairing within each age.
- DESeq2 internal dispersion estimation uses age-internal GEMgroups only.
- Each analysis is independent; failure of one age group does not block the others.

### Route B: Cross-Age log2FC Matrix Comparison (PRIMARY)

**Method**: Collect log2FC from each age-stratified DESeq2 result and construct a comparison matrix.

**Target gene log2FC matrix**:
```
rows = 231 target genes
columns = Young_log2FC_CA3vsCA1, Middle_log2FC_CA3vsCA1, Old_log2FC_CA3vsCA1
```

**Category-level log2FC matrix**:
```
rows = 13 categories
columns = same as above
values = median_log2FC of genes in that category for that age
```

**Direction consistency analysis** (descriptive, not statistical):
- **Age-concordant**: Same sign across all 3 ages (e.g., all CA3-high)
- **Age-flipped**: Sign changes between ages (e.g., CA3-high in Young, CA1-high in Old)
- **Old-enhanced**: Same direction as Young/Middle but larger magnitude in Old
- **Old-specific**: Only significant in Old, not in Young/Middle
- **Young-specific**: Only significant/exhibiting effect in Young

**Output**: `target_gene_log2fc_matrix_by_age.csv`, `target_gene_category_effect_by_age.csv`

**Why Route B is PRIMARY**: Directly answers Q4-Q5 with quantitative comparisons across ages. No new statistical tests — descriptive organization of existing DESeq2 results.

### Route C: Target Gene Category Summaries by Age (SECONDARY)

**Method**: For each of the 13 gene categories, compute per-age DE summary:

| Column | Description |
|--------|-------------|
| category_cn | Chinese category |
| category_en | English category |
| age | Young / Middle / Old |
| n_genes_available | Genes from this category found in Spatial |
| n_sig_padj05 | Number with padj < 0.05 |
| n_up_CA3 | Number with padj < 0.05 AND log2FC > 0 |
| n_down_CA1 | Number with padj < 0.05 AND log2FC < 0 |
| n_notable_abs_log2FC_gt_log2_1_5 | Number with |log2FC| > log2(1.5) |
| mean_log2FC | Mean of non-NA log2FC |
| median_log2FC | Median of non-NA log2FC |
| n_CA3_high | Number with log2FC > 0 |
| n_CA1_high | Number with log2FC < 0 |
| power_label | exploratory_low_power / primary_adequately_powered |

**Why SECONDARY**: Aggregates Route A results; critical for interpretation but no new modeling.

### Route D: Module Score Delta by Age (SECONDARY — Descriptive First)

**Method**: Use Phase 11 base `module_score_delta_by_age.csv` (already computed: 15 rows × 12 modules with Age column) or recompute from existing sample-level module scores.

**Default analysis (descriptive first)**:
1. Per module, per Age: n_GEMgroups, mean_delta, median_delta, IQR, direction (sign of median_delta)
2. Effect-size style comparison: Old - Young delta, Old - Middle delta (no p-values, no significance claims)
3. Present as table + boxplot; the table is the authoritative output

**Optional/sensitivity (only if explicitly requested or as supplementary)**:
- **Kruskal-Wallis test**: Delta ~ Age (nonparametric, 3 groups)
  - n_Young=3, n_Middle=4, n_Old=8
  - If computed: label all p-values as `exploratory_small_n`, column name `kw_pvalue_exploratory_not_primary_evidence`
  - **NOT used as primary evidence**; does NOT support mechanism interpretation
- **No pairwise post-hoc** due to n=3 in Young
- **No repeated-measures or mixed models**: GEMgroups differ across ages (GEMgroup and Age are colinear)

**Output**:
- `module_delta_by_age_group.csv`: columns = module, age, n, mean_delta, median_delta, IQR, direction, old_vs_young_delta, old_vs_middle_delta
- If KW computed: add `kw_statistic`, `kw_pvalue_exploratory_not_primary_evidence`
- `module_delta_by_age_boxplot.png`: 12 panels, boxplot CA3-CA1 delta by Age; title annotation "(n Young=3, Middle=4, Old=8, descriptive)"

**Why SECONDARY**: Descriptive trends across ages are informative. Any statistical test is underpowered and exploratory only.

### Route E: Region × Age Interaction Model (EXPLORATORY, Default SKIP)

**Method**: 
```
design = ~ Age + region_group + Age:region_group
```

**Critical limitation**: This design drops GEMgroup blocking. GEMgroup and Age are perfectly colinear (GEMgroups 1-4 = Young, 5-8 = Middle, 9-16 = Old). You cannot simultaneously block by GEMgroup AND test Age:region_group interaction in a standard DESeq2 model. Attempting `~ GEMgroup + Age + region_group + Age:region_group` would be rank-deficient because Age is a linear combination of GEMgroups.

**If attempted**:
1. Use `design = ~ Age + region_group + Age:region_group` (NO GEMgroup blocking)
2. Check full rank via `qr(model.matrix(...))$rank`
3. All 31 pseudobulk samples (all ages) for dispersion estimation — higher power for dispersion, but no within-GEMgroup pairing
4. Results: `results(dds, name = "AgeOld.region_groupCA3")` tests whether the CA3-vs-CA1 effect in Old differs from the reference (Young)
5. Label as `sensitivity_no_GEMgroup_blocking` with explicit caveat: "not the primary result; GEMgroup pairing lost; potential false confidence from spot-level correlation"

**Default decision (D_E)**: **SKIP the interaction model**. The age-stratified paired DESeq2 (Route A) + cross-age log2FC comparison (Route B) are sufficient and interpretable. The interaction model sacrifices within-GEMgroup pairing and conflates Age with GEMgroup.

**Rationale**: Three independent paired analyses (Young, Middle, Old) are more interpretable than one interaction model with confounded blocking. Users can compare effect sizes and direction across ages without the statistical machinery of an interaction term.

### Route F: Age-Stratified Coupling (EXPLORATORY)

**Method**: Compute Spearman coupling of sample-level module scores within each age separately.

**Sample sizes**:
- Young CA1: 4 GEMgroups (G1-G4) → 4 samples for coupling
- Young CA3: 3 GEMgroups (G2-G4) → 3 samples for coupling
- Middle CA1: 4 GEMgroups → 4 samples
- Middle CA3: 4 GEMgroups → 4 samples
- Old CA1: 8 GEMgroups → 8 samples
- Old CA3: 8 GEMgroups → 8 samples

**With n=3 or 4, Spearman correlations are extremely noisy**. A single outlier GEMgroup can flip the sign of rho.

**Default decision (D_F)**: Compute age-stratified coupling but label ALL results as `exploratory_n_3_to_8`. Do NOT apply BH adjustment within age (too few pairs). Do NOT compare coupling strengths between ages statistically. Only describe rho magnitude and direction per pair × age × context.

**Output table schema** (`module_age_coupling_summary.csv`):
| Column | Description |
|--------|-------------|
| pair_name | e.g., RPLvsRPS, CIvsCIII |
| age | Young / Middle / Old |
| context | CA1 / CA3 / delta |
| n_samples | number of GEMgroups used for cor (e.g., 3, 4, 8) |
| spearman_rho | Spearman correlation coefficient |
| small_n_flag | TRUE for all (n ≤ 8) |
| notes | "n too small for reliable Spearman estimation; rho reported as descriptive only" |

**p-value policy**:
- If a raw Spearman p-value is included at all, the column MUST be named `raw_pvalue_not_interpretable_small_n` and the notes column MUST state "p-value not interpretable due to n ≤ 8; do not use for inference".
- **Default**: omit p-values entirely. Output rho, n, context, age, pair_name, small_n_flag only.

**No dedicated coupling plot** for this route (n too small for meaningful visualization beyond the summary table).

**Why EXPLORATORY**: n per age is too small for reliable Spearman estimation. Any rho near ±1 with n=3 is trivially achievable by chance. Rho values are reported as descriptive patterns, not as statistical evidence of coupling.

---

## 3. Package and Environment Strategy

### 3.1 Current State (from Phase 11 Base)

| Package | Status | Version (if known) |
|---------|--------|-------------------|
| DESeq2 | Installed | (from renv) |
| Seurat 5.5.0 | Installed | 5.5.0 |
| Matrix | Installed | 1.7-3 |
| Matrix.utils | Installed (Phase 11) | (from renv) |
| apeglm | Installed (Phase 11) | (from renv) |
| corrplot | Installed (Phase 11) | (from renv) |
| ggplot2 | Installed | 4.0.3 |
| ggrepel | Installed | (from renv) |
| pheatmap | Installed | 1.0.13 |
| data.table | Installed | 1.18.4 |
| digest | Installed | 0.6.39 |
| matrixStats | Installed? | check in renv |
| ComplexHeatmap | NOT installed | may install if needed |

### 3.2 Package Policy

| Rule | Description |
|------|-------------|
| Install allowed | YES — use renv::install() + renv::snapshot() |
| Record protocol | Before/after renv.lock MD5, packageVersion(), log in provenance |
| Prohibited | STutility, Tangram, Python, clusterProfiler, enrichR, cell-type deconvolution tools |
| May install | ComplexHeatmap (if pheatmap insufficient for category-split heatmaps), scales (for color), cowplot/patchwork (for multi-panel) |

### 3.3 renv.lock Baseline

Current renv.lock MD5: `b2bf05038b440576855d043e137c5883` (post-Phase 11 with apeglm + corrplot installed). This is the execution starting point.

---

## 4. Input Files

### 4.1 Required (STOP if missing)

| # | File | Use |
|---|------|-----|
| 1 | `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` | Source of Spatial raw counts for pseudobulk aggregation |
| 2 | `docs/target_genes.csv` | Target gene list (251 genes, 13 categories, GBK encoding) |
| 3 | `data/processed/spatial/phase10_target_gene_audit/target_gene_input_audit.csv` | Phase 10 gene audit (category mapping, mapping_status) |
| 4 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/pseudobulk_sample_manifest.csv` | GEMgroup-region-age mapping, paired status per GEMgroup |
| 5 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds` | Existing pseudobulk counts (genes × 31 samples) — can be subset by age instead of recomputing |

### 4.2 Recommended (WARNING if missing)

| # | File | Use |
|---|------|-----|
| 6 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_scores_sample_summary.csv` | Sample-level module scores (reuse to avoid recomputing AddModuleScore) |
| 7 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_score_delta_by_age.csv` | Module delta per GEMgroup × Age |
| 8 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/deseq2_all_genes_CA3_vs_CA1.csv` | Phase 11 base all-age DE results (reference comparison) |
| 9 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/missing_gene_symbol_rescue_audit.csv` | Phase 11 missing gene rescue audit |
| 10 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv` | Module gene set definitions |

### 4.3 Reference (not for loading)

| # | File | Use |
|---|------|-----|
| 11 | `docs/spatial_phase11_ca1_ca3_de_module_coupling_plan.md` | Phase 11 base plan (design reference) |
| 12 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/phase11_final_report.txt` | Phase 11 base results (caveats, boundaries) |
| 13 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/phase11_validation_summary.txt` | Phase 11 base validation status |
| 14 | `data/processed/spatial/phase10_target_gene_audit/phase10_validation_summary.txt` | Phase 10 validation status |

---

## 5. Output Design

### 5.1 Output Directories

```
data/processed/spatial/phase11_age_grouping/
figures/spatial/phase11_age_grouping/
```

### 5.2 Data Outputs

| # | File | Description |
|---|------|-------------|
| 1 | `age_group_pseudobulk_sample_manifest.csv` | Per-age sample manifest: GEMgroup, region_group, Age, n_spots, has_pair, age_stratum |
| 2 | `design_matrix_audit_by_age.csv` | Per-age design rank check, model_type (paired_blocked / fallback_unpaired / blocked_insufficient_pairs), power_label, sample counts. Critically: if an age uses fallback_unpaired, all output files, figures, and the analysis guide must explicitly label it as such and must NOT present it as interchangeable with the primary paired result |
| 3 | `deseq2_all_genes_CA3_vs_CA1_Young.csv` | Full DESeq2 results for Young (n=6, 3 pairs). Must contain columns: `model_type` = paired_blocked / fallback_unpaired / blocked_insufficient_pairs; `power_label` = exploratory_low_power |
| 4 | `deseq2_all_genes_CA3_vs_CA1_Middle.csv` | Full DESeq2 results for Middle (n=8, 4 pairs). Must contain columns: `model_type`, `power_label` = exploratory_low_power |
| 5 | `deseq2_all_genes_CA3_vs_CA1_Old.csv` | Full DESeq2 results for Old (n=16, 8 pairs). Must contain columns: `model_type`, `power_label` = primary_adequately_powered |
| 6 | `deseq2_target_genes_CA3_vs_CA1_by_age.csv` | 251 target genes × 3 ages: log2FC, padj, direction per age |
| 7 | `target_gene_log2fc_matrix_by_age.csv` | 231 genes × 3 age columns (wide format log2FC matrix) |
| 8 | `target_gene_category_effect_by_age.csv` | Per-category × per-age: n_sig, n_up, n_down, mean_log2FC, median_log2FC |
| 9 | `age_group_deg_summary.csv` | Per-age DEG counts: n_tested, n_sig, n_up_CA3, n_down_CA1, n_target_sig |
| 10 | `age_group_low_power_flags.csv` | Per-age power assessment: n_pairs, DESeq2 convergence, caveats |
| 11 | `module_delta_by_age_group.csv` | Per-module delta by age: n, mean_delta, median_delta, IQR, direction, old_vs_young_delta, old_vs_middle_delta. Optional KW columns: kw_statistic, kw_pvalue_exploratory_not_primary_evidence |
| 12 | `module_age_coupling_summary.csv` | Per-coupling-pair × per-age: n_samples, spearman_rho, context, age, pair_name, small_n_flag. No p-values by default; if included, column MUST be named `raw_pvalue_not_interpretable_small_n` |
| 13 | `optional_interaction_model_audit.csv` | If interaction attempted: rank check, convergence, caveats (else "SKIPPED_BY_DEFAULT") |
| 14 | `phase11_age_grouping_validation_summary.txt` | Check results (PASS/WARNING/STOP) |
| 15 | `phase11_age_grouping_provenance.csv` | Run metadata: versions, timings, renv.lock MD5, packages |

### 5.3 Figure Outputs

| # | File | Description |
|---|------|-------------|
| 1 | `target_gene_log2fc_heatmap_by_age.png` | 231 genes × 3 ages: log2FC heatmap, rows split by category, annotated for significance |
| 2 | `target_gene_category_effect_heatmap_by_age.png` | 13 categories × 3 ages: median/mean log2FC heatmap |
| 3 | `volcano_CA3_vs_CA1_Young.png` | Volcano: Young DESeq2 results, target genes highlighted by category |
| 4 | `volcano_CA3_vs_CA1_Middle.png` | Volcano: Middle DESeq2 results |
| 5 | `volcano_CA3_vs_CA1_Old.png` | Volcano: Old DESeq2 results |
| 6 | `deg_count_by_age_barplot.png` | Barplot: n_up_CA3, n_down_CA1, n_sig per age, split by all-genes / target-genes |
| 7 | `module_delta_by_age_boxplot.png` | 12 panels (one per module): CA3-CA1 delta boxplot by Age (Young=3, Middle=4, Old=8 points) |
| 8 | `selected_target_gene_effect_by_age_dotplot.png` | Selected genes (top N by |log2FC| or user-specified): dot = log2FC, x = Age, color=category |

### 5.4 Documentation Output

| # | File | Description |
|---|------|-------------|
| 1 | `docs/spatial_phase11_age_grouping_analysis_guide.md` | User-facing analysis guide (generated during execution) |

---

## 6. Documentation: `spatial_phase11_age_grouping_analysis_guide.md` (Template)

This guide must be generated during execution and must contain:

### 6.1 What Each Module Asks

| Module | Question | Input | Output |
|--------|----------|-------|--------|
| Age-stratified DESeq2 | CA3 vs CA1 DEG within each age | Pseudobulk counts subset by age | `deseq2_all_genes_CA3_vs_CA1_{Age}.csv` |
| Target gene by age | Which target genes are DE in each age? | DESeq2 results + Phase 10 audit | `deseq2_target_genes_CA3_vs_CA1_by_age.csv` |
| Cross-age log2FC matrix | How do CA3-vs-CA1 effects change with age? | 3 age-specific DESeq2 results | `target_gene_log2fc_matrix_by_age.csv` |
| Category by age | Do whole categories shift with age? | Cross-age matrix + categories | `target_gene_category_effect_by_age.csv` |
| Module delta by age | Do module-score regional differences vary by age? | Phase 11 `module_score_delta_by_age.csv` | `module_delta_by_age_group.csv` |

### 6.2 How to Read log2FC

- `log2FoldChange > 0` → higher expression in **CA3** than CA1
- `log2FoldChange < 0` → higher expression in **CA1** than CA3
- `baseMean` → mean of normalized counts across all pseudobulk samples in that age
- `padj` → Benjamini-Hochberg adjusted p-value (within that age's DESeq2 run only)

### 6.3 Thresholds

| Term | Definition |
|------|-----------|
| `padj_sig` | padj < 0.05 |
| `log2FC_notable` | |log2FC| > log2(1.5) ≈ 0.585 |
| `CA3_high` | log2FC > 0 (direction label, not a significance claim) |
| `CA1_high` | log2FC < 0 (direction label, not a significance claim) |

### 6.4 Why Spots Are Not Replicates

Visium spots from the same tissue section (GEMgroup) share:
- Same animal, same tissue processing, same sequencing library
- Spatial neighbors have correlated gene expression

Treating individual spots as replicates would dramatically inflate false positive rates. The GEMgroup (tissue section) is the biological replicate unit.

### 6.5 Why Young and Middle Have Low Statistical Power

| Age | Paired GEMgroups | DESeq2 samples | Issue |
|-----|-----------------|----------------|-------|
| Young | 3 (G2,G3,G4) | 6 | DESeq2 can fit the model but dispersion estimates are unstable with 3 levels of the blocking factor |
| Middle | 4 (G5-G8) | 8 | Better than Young but still below typical N=5-6 minimum for reliable DE |
| Old | 8 (G9-G16) | 16 | Adequate for DESeq2 with GEMgroup blocking |

**Consequences**:
- In Young, DESeq2 may detect only very large effects (e.g., |log2FC| > 1-2).
- In Middle, moderate effects may reach significance; small effects will not.
- In Old, the same genes tested may show effects that are undetectable in Young/Middle simply due to sample size, not biology. A gene not significant in Young does NOT mean it is "not differentially expressed in young tissue."

### 6.6 Age-Stratified DE vs Interaction Model

| Aspect | Age-stratified (Routes A+B) | Interaction (Route E) |
|--------|---------------------------|----------------------|
| GEMgroup blocking | Preserved per age | Dropped |
| Replicate unit | GEMgroup | Individual pseudobulk samples (weakened) |
| Interpretation | "CA3 vs CA1 in Young", "CA3 vs CA1 in Old" — simple | "Does the CA3-CA1 effect differ between Old and Young?" — requires interaction contrast |
| Power | Lower per age (fewer samples), but valid | Higher total samples but confounded |
| Default choice | **PRIMARY** | **SKIP** |

### 6.7 What Can Be Used for Downstream Biological Interpretation

| Result type | Status | For interpretation? |
|-------------|--------|---------------------|
| Old-specific DESeq2 results (padj<0.05, |log2FC|>log2(1.5)) | PRIMARY, adequately powered | YES — candidate genes for Old-specific CA3-vs-CA1 effects |
| Middle-specific DESeq2 results (padj<0.05) | EXPLORATORY, 4 pairs | With caution; cross-reference with Old results |
| Young-specific DESeq2 results (padj<0.05) | EXPLORATORY, 3 pairs | With extreme caution; only large effects credible |
| Cross-age log2FC direction consistency | DESCRIPTIVE | YES — direction concordance/discordance across ages is informative regardless of significance |
| Module delta trend across ages | DESCRIPTIVE | YES — direction and magnitude trends; KW p-values are exploratory |
| Age-stratified coupling | EXPLORATORY, n=3-8 | NOT recommended for interpretation; n too small |

### 6.8 Complex V / Missing Genes Caveat

16 Atp5* genes (Complex V) are absent from the Spatial assay. They are likely named differently in the annotation version used by this Seurat object (e.g., `Atp5a1` instead of `Atp5f1a`). Phase 11 rescue audit found substring candidates but no user-confirmed replacements. Complex V is excluded from module scores (11 non-Complex-V modules). Any statement about "mitochondrial gene expression" in these results EXCLUDES Complex V (ATP synthase).

### 6.9 Next Steps After Phase 11 Age Grouping

From this phase, the user can select:
1. Candidate genes with age-concordant or age-divergent CA3-CA1 effects for CA3/DG mitochondrial interpretation
2. Categories that show age-dependent shifts for biological hypothesis generation
3. Old-specific target genes for focused spatial visualization
4. A foundation for Phase 12: DG-inclusive region comparison or biological interpretation phase

---

## 7. Validation Checks

| ID | Check | Criterion | Severity |
|----|-------|-----------|----------|
| V01 | Phase 11 base outputs readable | pseudobulk_sample_manifest.csv, deseq2_all_genes_CA3_vs_CA1.csv, module_score_delta_by_age.csv all loadable | STOP if not |
| V02 | Hippo RDS exists and is loadable | 9,921 spots, 32,285 Spatial features | STOP if not |
| V03 | Age groups present in manifest | Young, Middle, Old all have ≥1 CA1 and ≥1 CA3 sample | STOP if not |
| V04 | Paired GEMgroups per age recorded | Young=3, Middle=4, Old=8 (from Phase 11 manifest) | WARNING if mismatch |
| V05 | Young paired DE attempted only if ≥2 pairs | G2,G3,G4 all have CA3 → 3 pairs, proceed | STOP if < 2 |
| V06 | Middle paired DE attempted only if ≥2 pairs | G5-G8 all have CA3 → 4 pairs, proceed | STOP if < 2 |
| V07 | Old paired DE attempted only if ≥2 pairs | G9-G16 all have CA3 → 8 pairs, proceed | STOP if < 2 |
| V08 | Design matrix full rank per age | qr()$rank == ncol() for each age's design | STOP if not (for that age) |
| V09 | DESeq2 output has rows | nrow(res) ≈ nrow(pseudobulk_counts) - filtered | WARNING if empty |
| V10 | Target genes merged correctly | 231 found + 20 missing = 251 rows in target gene by-age table | WARNING if mismatch |
| V11 | Missing genes retained with NA DE values | 20 genes with in_spatial=FALSE, DE columns NA | WARNING if missing |
| V12 | No full spot-level dense conversion | Per-gene sparse extraction or pseudobulk-only dense conversion (≤31 cols) | STOP if violated |
| V13 | No spots as replicates | All DESeq2 uses pseudobulk aggregation by GEMgroup | STOP if violated |
| V14 | No full Seurat object saved | No saveRDS(hippo), save.image(), .RData | PASS/FAIL |
| V15 | No WholeBrain/DG loaded | Single Hippo object only | PASS/FAIL |
| V16 | No enrichment run | No GO/KEGG/Reactome calls | PASS/FAIL |
| V17 | No enrichment/Tangram/Python | No python/reticulate/STutility/cell-type deconvolution | PASS/FAIL |
| V18 | renv.lock changes recorded | Before/after MD5, packages installed, snapshot() logged | WARNING if unlogged |
| V19 | Rplots.pdf handled | Pre-existence recorded; deleted only if created during run | PASS/FAIL |
| V20 | Young/Middle results labeled exploratory_low_power | In CSV columns and plot annotations | WARNING if not labeled |
| V21 | Complex V status documented | Caveat carried forward from Phase 11 | WARNING if missing |
| V22 | Cross-age log2FC matrix has values for all 231 genes | No row should be all-NA (unless gene absent from all ages) | WARNING if unexpected NAs |
| V23 | target_genes.csv NOT modified | MD5 unchanged from Phase 11 run | PASS/FAIL |
| V24 | All output files in correct directories | data/ and figures/ subdirectories exist and populated | WARNING if missing |
| V25 | No biological interpretation language in validation/provenance | "Statistical description only" | PASS/FAIL |
| V26 | model_type and power_label present in all age-specific DE outputs | Each DESeq2 CSV contains model_type and power_label columns; design_matrix_audit_by_age.csv records model_type per age | WARNING if missing |
| V27 | Fallback model results explicitly labeled in outputs and plots | If any age uses fallback_unpaired, volcano titles and figure captions contain "unpaired fallback" label | WARNING if not labeled |
| V28 | renv.lock status vs git HEAD recorded in provenance | provenance.csv has columns: renv_lock_status_vs_git_head_start, renv_lock_status_vs_git_head_end, renv_lock_changed_during_run | WARNING if missing |

---

## 8. Stop/Warning Conditions

| ID | Condition | Action |
|----|-----------|--------|
| SC-1 | Hippo RDS missing or unreadable | STOP — object unavailable |
| SC-2 | target_genes.csv missing | STOP — Phase 10 prerequisite |
| SC-3 | Phase 11 pseudobulk manifest missing | STOP — cannot determine age-GEMgroup mapping |
| SC-4 | Phase 11 pseudobulk counts RDS missing | STOP — rebuild from Hippo RDS, must document |
| SC-5 | Any age group with < 2 paired GEMgroups | Do NOT run paired DE for that age; mark as `blocked_insufficient_pairs`; continue other ages |
| SC-6 | Design matrix rank-deficient for an age | STOP that age's paired model; attempt unpaired fallback; label as `fallback_unpaired` |
| SC-7 | DESeq2 failure for a specific age | STOP that age only; continue other ages if independent and documented |
| SC-8 | Young paired DE produces 0 significant genes | NOT a stop condition — expected for n=3; record and continue |
| SC-9 | Middle paired DE produces 0 significant genes | NOT a stop condition — possible with n=4; record and continue |
| SC-10 | Young/Middle low n (<5 pairs) | WARNING — label results exploratory_low_power; do NOT stop |
| SC-11 | Complex V still missing (0 genes) | WARNING — documented caveat; do NOT stop |
| SC-12 | Interaction model rank-deficient | SKIP with WARNING — expected behavior per Route E |
| SC-13 | Coupling n ≤ 4 per age group | WARNING — label coupling exploratory; do NOT stop |
| SC-14 | User requests biological interpretation in Phase 11 Age Grouping output | STOP — Phase 12+ scope |
| SC-15 | User requests cell-type deconvolution (RCTD, cell2location, etc.) | STOP — requires separate plan |
| SC-16 | User requests WholeBrain loading | STOP — requires separate plan |
| SC-17 | User requests Seurat object saving | STOP — against AGENTS.md |
| SC-18 | User requests enrichment analysis (GO/KEGG/Reactome) | STOP — Phase 12+ scope |

---

## 9. Boundaries

### Phase 11 Age Grouping DOES

- Load `seurat_Visium_Hippo_All.rds` (same object as Phase 11 base) **OR** reuse `ca1_ca3_pseudobulk_counts.rds` from Phase 11 base (avoid recomputation if input is identical)
- Build age-stratified pseudobulk count matrices by subsetting existing pseudobulk counts by age
- Run independent DESeq2 per age: Young (3 pairs), Middle (4 pairs), Old (8 pairs)
- Design: `~ GEMgroup + region_group` per age (paired); fallback `~ region_group` if rank-deficient
- Annotate target genes from Phase 10 audit onto each age's DE results
- Compute cross-age log2FC matrix for 231 target genes (Routes A+B)
- Compute per-category effect summary by age (Route C)
- Analyze module delta by age from existing module_score_delta_by_age.csv (Route D)
- Compute age-stratified coupling descriptively (Route F, exploratory)
- Label ALL Young and Middle results as `exploratory_low_power`
- Generate all QC/DE/module/category figures
- Generate `docs/spatial_phase11_age_grouping_analysis_guide.md`
- Log all package changes, renv.lock status, provenance
- Respect all Phase 11 base caveats (Complex V excluded, spots ≠ replicates, etc.)

### Phase 11 Age Grouping DOES NOT

- Load new spatial objects beyond Hippo RDS
- Run enrichment analysis (GO/KEGG/Reactome)
- Run Tangram, Python, STutility, or cell-type deconvolution
- Claim biological mechanism, causation, or disease relevance
- Treat spots as biological replicates
- Run interaction model (Route E) unless user explicitly overrides D_E
- Save full Seurat objects
- Generate Rplots.pdf
- Install packages not listed in §3.2
- Perform DG-specific DE (deferred to Phase 12+)
- Auto-substitute missing genes without user review
- Perform any single-cell-level analysis (pseudobulk only)

---

## 10. Execution Outline

```
Phase Spatial-11 Age Grouping

E0: Setup
  - Re-read this plan file
  - Rplots.pdf guard (start)
  - renv.lock MD5 before this run
  - renv.lock status vs git HEAD before this run: record whether already modified relative to git HEAD, and which packages differ
  - Git status at start: record `git status --short` output (tracked modifications, untracked files)
  - Load required packages (same set as Phase 11 base)
  - Check optional packages; install if beneficial and log
  - Load Phase 11 base outputs (manifest, pseudobulk counts, module deltas, gene audit)
  - Verify Phase 11 base validation status

E1: Input Audit
  - Verify Hippo RDS exists (for rebuild if needed)
  - Verify pseudobulk_sample_manifest.csv: 31 rows, 3 ages, paired counts
  - Verify ca1_ca3_pseudobulk_counts.rds: dimensions, sample_id alignment
  - Verify Phase 10 gene audit: 251 genes, 231 exact-match
  - Log V01-V03

E2: Age-Specific Sample Manifest
  - Extract per-age sample metadata from Phase 11 manifest
  - Build age_group_pseudobulk_sample_manifest.csv
  - Record paired GEMgroup counts: Young=3, Middle=4, Old=8
  - Flag GEMgroup 1 as "CA1_only" for Young
  - Save: age_group_pseudobulk_sample_manifest.csv
  - Log V04

E3: Design Audit Per Age
  - For each age: build design matrix for ~ GEMgroup + region_group
  - Check full rank via qr()
  - If rank-deficient: attempt ~ region_group fallback
  - Save: design_matrix_audit_by_age.csv
  - Check SC-5, SC-6
  - Log V05-V08

E4: Age-Stratified DESeq2 (Young)
  - Subset pseudobulk counts to Young samples (G2_CA1, G2_CA3, G3_CA1, G3_CA3, G4_CA1, G4_CA3)
  - Note: G1_CA1 excluded (no CA3 mate); if using 6-sample matrix, ensure design handles it
  - Run DESeq2: ~ GEMgroup + region_group
  - Filter: rowSums(counts) >= 10
  - LFC shrinkage: apeglm (if available) else "normal"
  - Write output with columns: model_type = paired_blocked (or fallback_unpaired / blocked_insufficient_pairs), power_label = exploratory_low_power
  - Save: deseq2_all_genes_CA3_vs_CA1_Young.csv
  - Check SC-7, SC-8

E5: Age-Stratified DESeq2 (Middle)
  - Subset pseudobulk counts to Middle samples (G5-G8, 8 samples)
  - Run DESeq2 as above
  - Write output with columns: model_type, power_label = exploratory_low_power
  - Save: deseq2_all_genes_CA3_vs_CA1_Middle.csv
  - Check SC-7, SC-9

E6: Age-Stratified DESeq2 (Old)
  - Subset pseudobulk counts to Old samples (G9-G16, 16 samples)
  - Run DESeq2 as above
  - Write output with columns: model_type, power_label = primary_adequately_powered
  - Save: deseq2_all_genes_CA3_vs_CA1_Old.csv

E7: Target Gene Annotation by Age
  - Merge each age's DE results with Phase 10 gene audit
  - Create combined table: 251 genes × 3 ages
  - Annotate with padj_sig, direction, lfc_notable per age
  - Save: deseq2_target_genes_CA3_vs_CA1_by_age.csv
  - Log V10, V11

E8: Cross-Age log2FC Matrix (Route B)
  - Build wide matrix: rows=231 genes, cols=Young/Middle/Old log2FC
  - Classify each gene: age-concordant, age-flipped, old-enhanced, old-specific, young-specific
  - Save: target_gene_log2fc_matrix_by_age.csv
  - Log V22

E9: Category Effect by Age (Route C)
  - For each of 13 categories × 3 ages:
    - n_available, n_sig, n_up_CA3, n_down_CA1, mean_log2FC, median_log2FC
  - Compile age_group_deg_summary.csv
  - Save: target_gene_category_effect_by_age.csv, age_group_deg_summary.csv
  - Save: age_group_low_power_flags.csv

E10: Module Delta by Age (Route D — descriptive first)
  - Load Phase 11 module_score_delta_by_age.csv or recompute from sample-level scores
  - Compute per module × per Age: n, mean_delta, median_delta, IQR, direction
  - Compute Old - Young delta, Old - Middle delta (effect-size comparison, no p-values)
  - Optionally, if explicitly requested: Kruskal-Wallis per module, label p-value as kw_pvalue_exploratory_not_primary_evidence
  - Save: module_delta_by_age_group.csv

E11: Age-Stratified Coupling (Route F, Exploratory)
  - Load Phase 11 module_scores_sample_summary.csv
  - For each age, for each of 10 coupling pairs:
    - Spearman correlation in CA1 context, CA3 context, delta context
    - Output: n_samples, spearman_rho, small_n_flag (TRUE for all), context, age, pair_name
    - Do NOT output raw p-values; if included, column MUST be raw_pvalue_not_interpretable_small_n
  - Do NOT apply BH adjustment within age (too few pairs)
  - Do NOT compare coupling strengths between ages statistically
  - Label all as exploratory_n_3_to_8
  - Save: module_age_coupling_summary.csv (rho-only output, no p-values by default)

E12: Plots
  - Figure 1: target_gene_log2fc_heatmap_by_age.png — pheatmap: rows=231 genes split by category, cols=3 ages, colored by log2FC
  - Figure 2: target_gene_category_effect_heatmap_by_age.png — 13 categories × 3 ages, colored by median_log2FC
  - Figure 3: volcano_CA3_vs_CA1_Young.png — ggrepel labels for target genes, title="Young (3 paired GEMgroups, exploratory)"
  - Figure 4: volcano_CA3_vs_CA1_Middle.png — as above
  - Figure 5: volcano_CA3_vs_CA1_Old.png — as above
  - Figure 6: deg_count_by_age_barplot.png — stacked bars per age: n_up_CA3, n_down_CA1, split all-genes/target-genes
  - Figure 7: module_delta_by_age_boxplot.png — 12 panels, boxplot CA3-CA1 delta by Age
  - Figure 8: selected_target_gene_effect_by_age_dotplot.png — dot=log2FC, errorbar=SE, x=Age, color=category, faceted or single panel
  - Plot conventions: CA1=#2166ac, CA3=#b2182b, Age=Young#c6dbef/Middle#4292c6/Old#08306b

E13: Analysis Guide
  - Generate docs/spatial_phase11_age_grouping_analysis_guide.md from template (§6)
  - Include all sections: module descriptions, log2FC reading guide, thresholds, replicate explanation, power analysis, caveats

E14: Validation and Provenance
  - Run checks V01-V25
  - Save: phase11_age_grouping_validation_summary.txt
  - Save: phase11_age_grouping_provenance.csv
  - renv.lock MD5 after run; log whether renv.lock changed during this run
  - Distinguish "unchanged during this run" from "clean relative to git HEAD":
    - Record renv.lock status vs git HEAD at start and end
    - If renv.lock was already modified vs git HEAD before this run, record that fact separately
    - If renv.lock changed during this run, record which packages installed and snapshot() called
  - Record git status at end; note any untracked output files
  - Record packageVersion() for all loaded packages

E15: Cleanup
  - Remove large objects; gc()
  - Rplots.pdf guard (end)
  - Final exit code: 0 if PASS or PASS_WITH_WARNINGS; 1 if any STOP
```

---

## 11. Recommended: Primary, Secondary, Exploratory Summary

| Priority | Route | What | Power |
|----------|-------|------|-------|
| **PRIMARY** | A | Age-stratified paired DESeq2 (3 independent runs) | Old=adequate, Middle=low, Young=very low |
| **PRIMARY** | B | Cross-age log2FC matrix + direction consistency | Descriptive, no new tests |
| **SECONDARY** | C | Target gene category summaries by age | Aggregates A+B |
| **SECONDARY** | D | Module delta by age (descriptive first; KW optional/sensitivity) | Descriptive trends informative; KW exploratory_small_n, not primary evidence |
| **EXPLORATORY** | E | Region × Age interaction model | SKIP by default (confounded with GEMgroup) |
| **EXPLORATORY** | F | Age-stratified coupling (rho only, no p-values) | n=3-8 per context; Spearman rho descriptive only |

---

## 12. Execution Command (Expected)

```bash
Rscript R/spatial/s11b_age_grouping_ca1_ca3_de.R
```

Execution phase MUST:
1. Re-read this plan file: `docs/spatial_phase11_age_grouping_plan.md`
2. Verify plan version matches expectations
3. Proceed checkpoint-by-checkpoint (not all at once)

Checkpoint plan:
- CP1 (E0-E1): Setup + input audit
- CP2 (E2-E3): Age-specific manifest + design audit per age
- CP3 (E4-E6): Age-stratified DESeq2 (Young → Middle → Old, sequentially)
- CP4 (E7-E9): Target gene annotation + log2FC matrix + category summaries
- CP5 (E10-E11): Module delta by age + age-stratified coupling
- CP6 (E12-E13): Plots + analysis guide
- CP7 (E14-E15): Validation + provenance + cleanup

---

> **Plan status**: v1.1. No analysis execution; only this plan file modified. No packages installed, no R scripts executed, no objects loaded, no renv.lock changes, no data/figures generated.
> **Next step**: User review → approval → execution.
> **Execution prerequisite**: Phase 11 base must be PASS_WITH_WARNINGS (confirmed). All required inputs must exist.
> **Script to create**: `R/spatial/s11b_age_grouping_ca1_ca3_de.R` (following E0-E15 outline above)
