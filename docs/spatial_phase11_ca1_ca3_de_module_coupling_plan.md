# Phase Spatial-11: CA1 vs CA3 Pseudobulk DE + Target-Gene Module/Coupling Analysis

> **Plan version**: v2.1 (2026-06-04)
> **Status**: READ-ONLY PLAN REVISION — no execution, no file writes
> **Handoff**: Phase 10 PASS (25 PASS, 0 FAIL, 1 WARNING) → Phase 11 planning
> **Phase 10 validation summary**: 231 exact-match genes, 20 missing (16 Atp5*, 1 Ndufb1, 3 Mrps), 2,130,552 nonzero spot-gene records/rows, DG=2364
> **bioskill**: NOT installed / NOT inspected — no content referenced or invented
> **Authority**: This saved `.md` file is the sole execution authority. If chat context conflicts with this plan, this plan wins. Execution phase must re-read this file first.

---

## 0. ref-bio Pre-Answer Report

### Task Classification

| Dimension | Detected |
|-----------|----------|
| Data modality | **spatial transcriptomics** (Visium, RDS-based, not raw Space Ranger) |
| Workflow stage | **differential expression** + **module score** + **coupling** |
| Package/tool | **Seurat**, **DESeq2**, **ggplot2**, **pheatmap/ComplexHeatmap**, **Matrix** |
| Platform | **Visium** (RDS-based reproduction) |

### Matched Source IDs

| Priority | Source ID | Upstream URL |
|----------|-----------|-------------|
| must_read | `seurat` | https://satijalab.org/seurat/ |
| must_read | `seurat-spatial` | https://satijalab.org/seurat/articles/spatial_vignette.html |
| must_read | `deseq2` | https://bioconductor.org/packages/release/bioc/html/DESeq2.html |
| must_read | `pseudobulk-de-guidance` | (upstream DESeq2 vignettes — pseudobulk workflow, `~ condition + sample` paired designs, BioC 2024 Seurat-DESeq2 pseudobulk vignette) |
| optional_read | `limma` | https://bioconductor.org/packages/release/bioc/html/limma.html |
| optional_read | `edger` | https://bioconductor.org/packages/release/bioc/html/edgeR.html |

---

## 1. Scientific Goals (Explicitly Stated)

**Primary goal**: Characterize CA1 vs CA3 regional differences in the expression of 231 user-provided target genes (mitochondrial, ribosomal, translation-factor) in the mouse hippocampus, using Visium spatial transcriptomics.

**Units of analysis**:
- **Spots** (n = 9,921): For visualization, spatial mapping, and spot-level module scores (exploratory only — spots are not independent replicates).
- **GEMgroup-level pseudobulk samples**: For all inferential statistics (DESeq2, module score comparisons, coupling correlations). A "sample" = one GEMgroup × region_group pseudobulk aggregate.
- **Region groups**: CA1, CA3 (primary comparison), DG (user interest + context), CA2 (context/QC only).

**DG role**: DG (= ML + GCL + Hilus, 2,364 spots) is included in summaries, plots, and exploratory analyses as a user-requested region of interest. It is NOT part of the primary CA1-vs-CA3 contrast. DG-specific DE (CA1 vs DG, CA3 vs DG) is **secondary** — computed only if the primary analysis succeeds and execution time permits.

**What this phase is NOT**: This is a descriptive and exploratory regional comparison. It does NOT claim cellular mechanisms, disease causation, aging mechanisms, or therapeutic targets. All coupling and age results are **exploratory, unadjusted or BH-adjusted, and labeled as such**.

---

## 2. Inputs

### 2.1 Spatial Object

| Field | Value |
|-------|-------|
| Object | `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` |
| Spots | 9,921 |
| Assays | Spatial (32,285 features, raw UMI counts), SCT (17,096 features, normalized) |
| Assay for DE | **Spatial `counts`** (raw UMI, sparse integer) |
| Region labels | Phase 06 validated: CA1=4671, CA3=2321, DG=2364 (ML+GCL+Hilus), CA2=565 |
| GEMgroups | 16 (Young=1-4, Middle=5-8, Old=9-16) |

### 2.2 Phase 10 Gene Audit

| Property | Value |
|----------|-------|
| Total input genes | 251 |
| exact_match (found in Spatial) | 231 |
| Missing genes | 20 (16 Atp5*, 1 Ndufb1, 3 Mrps) |
| Categories | 13 (Chinese + English) |
| Audit files | `target_gene_input_audit.csv` (251 rows), `gene_symbol_mapping_audit.csv` (252 rows) |

### 2.3 Phase 10 Summary Tables (for cross-reference)

| File | Rows | Use in Phase 11 |
|------|------|-----------------|
| `region_sample_summary.csv` | 13,167 | Cross-check pseudobulk aggregation |
| `region_age_summary.csv` | 2,772 | Reference for age-stratified summaries |
| `spot_expression_nonzero_long.csv.gz` | 2,130,552 | Spot-level reference (not for inference) |

### 2.4 Package Policy

**Phase 11 ALLOWS package installation via renv**. Do NOT force base-R-only or lightweight-only strategies. Prioritize correctness and clarity over minimal dependencies. The following packages may be installed during execution if needed:

| Package | Status in renv.lock | Use |
|---------|---------------------|-----|
| DESeq2 | **INSTALLED** | Core DE engine |
| Seurat 5.5.0 | Installed | Object loading, AddModuleScore |
| Matrix 1.7-3 | Installed | Sparse matrix support (dgCMatrix) |
| Matrix.utils | **MAY NEED INSTALL** | `aggregate.Matrix()` for pseudobulk aggregation (Phase 03 pattern) |
| ggplot2 4.0.3 | Installed | All plots |
| pheatmap 1.0.13 | Installed | Heatmaps (may supplement with ComplexHeatmap) |
| ggrepel | Installed | Volcano labels |
| data.table 1.18.4 | Installed | CSV output |
| digest 0.6.39 | Installed | Rplots.pdf MD5 |
| cowplot / patchwork | Installed | Multi-panel figures |
| scales 1.4.0 | Installed | Color scales |

**May be installed during execution** (if beneficial):

| Package | Justification | Impact |
|---------|--------------|--------|
| **apeglm** | LFC shrinkage for volcano plots; standard DESeq2 workflow | `renv::install("apeglm")`; renv.lock modified |
| **ComplexHeatmap** | Better heatmap control than pheatmap; supports row-split by category, annotations, row labels | `renv::install("ComplexHeatmap")`; renv.lock modified |
| **corrplot** | Correlation matrix visualization (alternative to pheatmap for coupling) | `renv::install("corrplot")`; renv.lock modified |
| **circlize** | Circular coupling plots (optional, if coupling network visualization is requested) | `renv::install("circlize")`; renv.lock modified |
| **matrixStats** | Efficient row/col medians, SDs for pseudobulk matrix | `renv::install("matrixStats")`; renv.lock modified |
| **limma** | `removeBatchEffect` for visualization-only batch assessment | `renv::install("limma")`; renv.lock modified |

**Installation protocol** (every time):
1. Record `renv.lock` MD5 before install
2. `renv::install("package")` 
3. Record `packageVersion("package")`
4. `renv::snapshot()`
5. Record new `renv.lock` MD5
6. Log in `phase11_provenance.csv`

**Prohibited**: STutility, Tangram, Python packages, DESeq2 alternatives that don't use pseudobulk aggregation.

---

## 3. Module A: Missing Gene Rescue Audit

### 3.1 Rationale

Phase 10 identified 20 genes not found in `rownames(hippo[["Spatial"]])`. Complex V (16 Atp5* genes) must NOT be permanently excluded without a rescue attempt. Gene symbols evolve across genome annotations (e.g., `Atp5f1a` may be `Atp5a1` in the Seurat object's annotation version).

### 3.2 Rescue Protocol

For each of the 20 missing genes:
1. Search `rownames(hippo[["Spatial"]])` for substrings matching the core gene name (e.g., grep for `Atp5` → collect all Atp5-family genes present in the object)
2. For `Ndufb*` genes: search all `Nduf*` genes present
3. For `Mrps*` and `Mrpl*` genes: search all mitoribosomal genes present
4. **Manual review candidates**: List genes found by substring search that are NOT in the user's original list but belong to the same family
5. **Never auto-add** genes to the user's list. Only record candidates.

### 3.3 Output

`data/processed/spatial/phase11_ca1_ca3_de_module_coupling/missing_gene_symbol_rescue_audit.csv`

| Column | Description |
|--------|-------------|
| `original_missing_gene` | Gene from user list that was not found |
| `category` | Chinese + English category |
| `rescue_candidates_in_object` | Comma-separated list of genes found by substring search (e.g., "Atp5a1,Atp5b,...") |
| `rescue_status` | `has_candidates` / `no_candidates_found` / `requires_ortholog_mapping` |
| `user_action_required` | TRUE for all — user must decide whether any candidate is an acceptable replacement |
| `notes` | Explanation of search strategy and caveats |

### 3.4 Complex V Decision Tree

```
After rescue audit:
  - If ≥ 3 rescue candidates confirmed by user → Complex V module computed with rescued genes, with provenance note
  - If 0 candidates after rescue AND user declines ortholog mapping → Complex V EXCLUDED with WARNING, module count = 11 (not 12)
  - If user requests ortholog mapping → separate Phase 11b sub-task
```

**Default** (if user doesn't respond): Document rescue candidates, mark `user_action_required = TRUE`, compute module scores for 12 modules (Complex V excluded), record as WARNING V_MODULE_COMPLEX_V_EXCLUDED. Do NOT block Phase 11.

---

## 4. Module B: CA1 vs CA3 Pseudobulk Differential Expression

### 4.1 Design

**Comparison**: CA3 vs CA1 (positive log2FC = higher in CA3 than CA1).
**Replicate unit**: GEMgroup (tissue section). Each GEMgroup contributes one pseudobulk sample per region_group.
**Assay**: Spatial `counts` (raw UMI, sparse integer, via `GetAssayData(hippo, assay = "Spatial", layer = "counts")`).

### 4.2 Pseudobulk Aggregation

**Preferred method**: `Matrix.utils::aggregate.Matrix()` — the same as Phase 03 (`s03_pseudobulk_deseq2_dg.R`). `Matrix.utils` provides `aggregate.Matrix()`; `Matrix` provides sparse dgCMatrix support. This preserves sparsity and correctly handles the dgCMatrix format.

```r
# Phase 03 pattern (reference), corrected for gene×spot orientation:
counts_sparse <- GetAssayData(hippo, assay = "Spatial", layer = "counts")  # genes × spots
colnames(counts_sparse) <- colnames(hippo)  # cell barcodes match metadata
spot_groupings <- metadata$GEMgroup_Region  # length = ncol(counts_sparse), e.g., "1_CA1", "1_CA3", ...

# aggregate.Matrix aggregates rows by grouping, so transpose to spots × genes first
pseudobulk_t <- Matrix.utils::aggregate.Matrix(
  t(counts_sparse), groupings = spot_groupings, fun = "sum"
)  # samples × genes (rows = unique GEMgroup_Region)
final_pseudobulk <- t(pseudobulk_t)  # genes × samples

# Verification:
stopifnot(nrow(final_pseudobulk) == nrow(counts_sparse))
stopifnot(setequal(colnames(final_pseudobulk), pseudobulk_sample_manifest$sample_id))
stopifnot(all(colnames(final_pseudobulk) == pseudobulk_sample_manifest$sample_id[
  match(colnames(final_pseudobulk), pseudobulk_sample_manifest$sample_id)]))
```

**Matrix orientation note**: `GetAssayData()` returns genes × spots (standard Seurat v5). `aggregate.Matrix()` aggregates rows by a grouping vector, so the count matrix must be transposed to spots × genes before aggregation, then transposed back. All verification uses the final genes × samples orientation.

**Fallback** (if `Matrix.utils` is unavailable or fails):
```r
# Manual sparse-column aggregation via split + rowSums
# Per-gene, per-GEMgroup_Region: sum raw UMI
# Same transpose logic: t(counts_sparse) → aggregate by spot_groupings → t() back
# Never convert full count matrix to dense
```

**Prohibited**: `as.matrix(counts_sparse)` on the full 32285 × 9921 matrix (would be ~2.6 GB dense). Per-gene numeric vector extraction is allowed. `Matrix.utils` is the aggregation package; `Matrix` is the sparse matrix engine.

### 4.3 Matched GEMgroup Audit

Before DESeq2, build and save `pseudobulk_sample_manifest.csv`:

| Column | Description |
|--------|-------------|
| sample_id | `G{gemgroup}_{region_group}` (e.g., G1_CA1, G1_CA3) |
| GEMgroup | 1-16 |
| region_group | CA1 or CA3 |
| Age | Young/Middle/Old |
| n_spots | Number of Visium spots aggregated |
| has_pair | TRUE if the GEMgroup has BOTH CA1 and CA3 spots |
| low_coverage | TRUE if n_spots < 20 |

**Stop condition**: If < 2 GEMgroups have `has_pair == TRUE`, STOP primary paired DE.

### 4.4 GEMgroup Inclusion Rule

**Expected**: Most GEMgroups (2-16) have spots in both CA1 and CA3. GEMgroup 1 may have only CA1 (no CA3). The design `~ GEMgroup + region_group` handles this:
- GEMgroups with only one region contribute to dispersion estimation but not to the region contrast.
- This is standard DESeq2 behavior for unbalanced paired designs.

### 4.5 Design Matrix

**Primary model** (paired, preferred):
```
design = ~ GEMgroup + region_group
```
- GEMgroup: blocking factor (unordered factor, up to 16 levels)
- region_group: CA1 (reference level) vs CA3
- Contrast of interest: `results(dds, contrast = c("region_group", "CA3", "CA1"))`

**Alternative model** (if primary model rank-deficient):
```
design = ~ region_group
```
- Simple two-group, no pairing. Less powerful but always valid.

**Design validation** (MANDATORY before DESeq2):
1. Build `model.matrix(~ GEMgroup + region_group, data = colData)`
2. Check `qr(model_matrix)$rank == ncol(model_matrix)` (full rank)
3. Check for zero-columns (e.g., all samples in one level)
4. Count samples per GEMgroup × region_group cell
5. Save `design_matrix_audit.csv`

### 4.6 DESeq2 Execution

```r
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,  # genes × samples (integer)
  colData   = sample_metadata,    # sample × (GEMgroup, region_group, Age)
  design    = ~ GEMgroup + region_group
)
# Filter: keep genes with rowSums(counts(dds)) >= 10 across all pseudobulk samples
# (threshold inherited from Phase 03; for ~32 samples, 10 total counts is defensible)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
```

### 4.7 Results Extraction

```r
res <- results(dds, contrast = c("region_group", "CA3", "CA1"), alpha = 0.05)
# log2FoldChange > 0 → higher in CA3
# log2FoldChange < 0 → higher in CA1

# Shrinkage (if apeglm available):
res_shrunken <- lfcShrink(dds, coef = "region_group_CA3_vs_CA1", type = "apeglm")
# If apeglm not available: use type = "normal" as fallback, record WARNING
```

**Significance thresholds**:
- `padj < 0.05` for statistical significance
- `|log2FC| > log2(1.5)` (~0.585) for **notable** genes (NOT a statistical filter; used only for flagging/annotation)
- Both criteria recorded; genes can be `padj_sig`, `log2FC_notable`, both, or neither

### 4.8 Target Gene Annotation

Merge Phase 10 gene audit with DESeq2 results:
- `in_target_list` = TRUE for the 231 exact-match genes
- `category_cn`, `category_en` from Phase 10
- `mapping_status` from Phase 10
- Missing genes (20) get their own row in `deseq2_target_genes_CA3_vs_CA1.csv` with `in_spatial = FALSE` and all DE columns NA

---

## 5. Module C: Target Gene Focused DE

### 5.1 Output Tables

| File | Description |
|------|-------------|
| `pseudobulk_sample_manifest.csv` | All pseudobulk samples: GEMgroup, region, Age, n_spots, has_pair, low_coverage |
| `design_matrix_audit.csv` | Design matrix columns, rank check, sample-to-column mapping |
| `ca1_ca3_pseudobulk_counts.rds` | Pseudobulk count matrix (genes × samples, sparse where possible) |
| `deseq2_all_genes_CA3_vs_CA1.csv` | Full DESeq2 results: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj |
| `deseq2_all_genes_CA3_vs_CA1_shrunken.csv` | LFC-shrunken version (only if apeglm or normal shrinkage succeeds) |
| `deseq2_target_genes_CA3_vs_CA1.csv` | 231 found genes + 20 missing genes, annotated with category |
| `target_gene_effect_summary_by_category.csv` | Per-category: n_genes_available, n_up_CA3, n_down_CA1, n_sig, mean_log2FC |
| `missing_gene_symbol_rescue_audit.csv` | Rescue candidates for 20 missing genes (see §3) |

### 5.2 Per-Category Effect Summary

| Column | Description |
|--------|-------------|
| category_cn | Chinese category |
| category_en | English category |
| n_genes_total | Genes in this category from user list |
| n_genes_available | Found in Spatial assay |
| n_genes_missing | Not found (0 for most, 16 for Complex V, 1 for Complex I, 3 for Mitoribosomal) |
| n_up_CA3_padj05 | Number with padj < 0.05 AND log2FC > 0 |
| n_down_CA1_padj05 | Number with padj < 0.05 AND log2FC < 0 |
| n_notable_CA3 | Number with |log2FC| > log2(1.5) AND log2FC > 0 (regardless of padj) |
| n_notable_CA1 | Number with |log2FC| > log2(1.5) AND log2FC < 0 (regardless of padj) |
| mean_log2FC | Mean of non-NA log2FC for available genes |
| median_log2FC | Median of non-NA log2FC for available genes |
| flag | `all_missing` / `zero_sig` / `ok` |

---

## 6. Module D: Module Score Analysis

### 6.1 Module Definitions (from Phase 10 Categories)

| Module | Source genes | Available | Min required | Status |
|--------|-------------|-----------|-------------|--------|
| mtDNA-encoded | 9 | 9 | 3 | OK |
| Complex I | 34 | 33 (1 missing: Ndufb1) | 3 | OK |
| Complex III | 11 | 11 | 3 | OK |
| Complex IV | 13 | 13 | 3 | OK |
| Complex V | 16 | **0** (before rescue audit) | **3** | **PENDING RESCUE** (§3) |
| TCA cycle | 12 | 12 | 3 | OK |
| Mitoribosomal | 26 | 23 (3 missing) | 3 | OK |
| Membrane/Translocator | 15 | 15 | 3 | OK |
| Defense/Antioxidant | 29 | 29 | 3 | OK |
| RPL | 37 | 37 | 3 | OK |
| RPS | 27 | 27 | 3 | OK |
| Eif | 16 | 16 | 3 | OK |
| Eef | 6 | 6 | 3 | OK |

**Total modules**: 12 (if Complex V excluded after rescue audit) or 13 (if Complex V genes rescued and confirmed).

### 6.2 Spot-Level Module Scores (Exploratory / Visualization Only)

```r
# Seurat::AddModuleScore on the full Hippo object
# CAUTION: computationally expensive for 12+ modules on 9921 spots
# Batch the calls if memory pressure is high
hippo <- AddModuleScore(hippo, features = module_gene_sets, name = "Module_", ctrl = 100)
# Results stored in hippo@meta.data as Module_1, Module_2, ...
```

**Use for**: spatial maps of module activity, spot-level distribution plots, exploratory visualization.

**NOT for**: any statistical inference. Spots are not independent replicates.

### 6.3 Sample-Level Module Scores (Inferential)

From spot-level module scores, compute per-GEMgroup × region_group summary:

```r
# Aggregate to sample level:
# mean_module_score = mean of per-spot scores within the GEMgroup × region_group
# median_module_score = median of per-spot scores
```

**Sample-level scores are the input for**: CA3-CA1 comparisons, coupling correlations, age trends.

### 6.4 CA3-CA1 Module Score Comparison

For each module:
1. Compute **paired delta**: `CA3_score - CA1_score` per GEMgroup (only for GEMgroups with both regions)
2. **Wilcoxon signed-rank test** on the paired deltas (default, nonparametric)
3. If insufficient paired GEMgroups: **Wilcoxon rank-sum test** (two-group, unpaired)
4. **No paired t-test** unless n ≥ 30 samples and normality is plausible (unlikely with 16 GEMgroups)

Output: `module_score_delta_CA3_minus_CA1.csv`

| Column | Description |
|--------|-------------|
| module | Module name |
| n_pairs | Number of GEMgroups with both CA1 and CA3 |
| mean_delta | Mean (CA3 - CA1) across GEMgroups |
| median_delta | Median (CA3 - CA1) |
| sd_delta | SD of deltas |
| wilcoxon_statistic | V statistic |
| wilcoxon_pvalue | Raw p-value |
| wilcoxon_padj | BH-adjusted across modules |
| test_type | `paired` or `two_group` |
| notes | Caveats (e.g., "Module excluded: 0 genes") |

---

## 7. Module E: Spearman Coupling Analysis

### 7.1 Statistical Unit

**GEMgroup-level module scores** (NOT spot-level). Each row in the coupling matrix = one matched GEMgroup with both CA1 and CA3 scores.

### 7.2 Pre-Specified Coupling Pairs

| # | Pair | Acronym | Rationale |
|---|------|---------|-----------|
| 1 | Mitoribosomal vs RPL | MRPLvsRPL | Mitochondrial vs cytosolic large subunit |
| 2 | Mitoribosomal vs RPS | MRPLvsRPS | Mitochondrial vs cytosolic small subunit |
| 3 | Complex I vs Complex III | CIvsCIII | ETC coordination |
| 4 | Complex I vs Complex IV | CIvsCIV | ETC coordination across complexes |
| 5 | mtDNA-encoded vs avg(Complex I-V, TCA) | mtNUC | mtDNA-nuclear coordination |
| 6 | avg(all mito modules) vs avg(RPL, RPS) | MITOvsCYTO | Global mito-cytosolic |
| 7 | RPL vs RPS | RPLvsRPS | Cytosolic ribosome subunit coordination |
| 8 | Defense/Antioxidant vs avg(Complex I, III, IV) | ROSvsETC | Oxidative stress vs respiration |
| 9 | TCA vs avg(Complex I, III, IV) | TCAvsETC | Substrate supply vs consumption |
| 10 | Membrane/Translocator vs avg(Complex I-V) | TRANSvsETC | Transport vs energy production |

### 7.3 Coupling Computation

```r
# For each coupling pair (module_a, module_b):
# Spearman correlation of sample-level module scores
# Computed in three contexts:
#   1. CA1 only (CA1 scores across GEMgroups)
#   2. CA3 only (CA3 scores across GEMgroups)
#   3. CA3-CA1 delta (delta scores across GEMgroups)
cor_test <- cor.test(scores_a, scores_b, method = "spearman")
```

**Multiple testing**: Apply BH correction within each context (CA1, CA3, delta) separately. Report both raw p and BH-adjusted p.

**If n < 8 GEMgroups for a context**: Label correlation as `exploratory_small_n`, still compute but flag.

### 7.4 Outputs

| File | Description |
|------|-------------|
| `module_coupling_spearman_ca1.csv` | Spearman rho, p, padj for CA1-only context |
| `module_coupling_spearman_ca3.csv` | Same for CA3-only context |
| `module_coupling_spearman_delta.csv` | Same for CA3-CA1 delta context |
| `module_coupling_combined_long.csv` | All three contexts in melted long format |

---

## 8. Module F: Age Coupling / Age Trend (Secondary)

### 8.1 Design

Age is a **secondary** analysis axis. The primary DE contrast is CA3 vs CA1 (blocked by GEMgroup). Age is explored descriptively.

### 8.2 Sample Structure (Inherited from Experimental Design)

| Age | GEMgroups | N pairs (GEMgroups with both CA1 and CA3) |
|-----|-----------|---|
| Young | 1, 2, 3, 4 | 4? (depends on audit; GEMgroup 1 may lack CA3) |
| Middle | 5, 6, 7, 8 | 4? |
| Old | 9-16 | 8? |

**Imbalance**: Old has ~2× GEMgroups. Documented in all output. Not corrected in Phase 11.

### 8.3 Age-Stratified Summaries

For each module:
- Mean/median module score by Age × region_group
- CA3-CA1 delta stratified by Age

### 8.4 Optional Interaction Model (High-Risk, Skip if Rank-Deficient)

```
design = ~ Age + region_group + Age:region_group
```

- **Prerequisite**: ≥ 2 GEMgroups per Age × region_group cell
- **Rank check**: `qr(model.matrix(...))$rank` must equal ncol
- **Default decision**: Unless the rank check clearly passes, **skip the interaction model** and rely on descriptive age-stratified summaries instead
- Do NOT force an interaction model that violates design constraints

### 8.5 Outputs

| File | Description |
|------|-------------|
| `age_module_summary.csv` | Mean module score by Age × region_group |
| `age_region_delta_summary.csv` | CA3-CA1 delta by GEMgroup with Age column |
| `age_coupling_delta_plots.png` | Dot/box plots of delta by Age, one panel per module |

---

## 9. Module G: Marker / Marker-Coupling Analysis (Optional, Deferred by Default)

### 9.1 Current Status

The user's target gene list (251 genes, 13 categories) contains functional gene modules — NOT cell-type marker genes. No external marker gene list has been provided.

### 9.2 Marker Analysis Gate

Phase 11 does NOT auto-generate marker lists from databases. If marker analysis is desired:

1. User must provide a marker gene source (author code, literature, database export)
2. Source must be documented in Phase 11 plan update
3. Marker genes are audited against `rownames(hippo[["Spatial"]])` (same pattern as Phase 10)
4. Marker coupling = Spearman correlation of sample-level marker module scores vs mitochondrial module scores

**Default**: Marker analysis is **DEFERRED** to Phase 12+. Phase 11 proceeds with the 13 functional categories as defined in Phase 10.

### 9.3 If Deferred, No Marker Coupling

The coupling pairs in §7.2 are sufficient for Phase 11 scope. Tracing additional coupling (e.g., cell-type markers vs mitochondrial modules) requires a marker source and a plan update.

---

## 10. Module H: Plot Plan

### 10.1 QC / Diagnostic Plots

| Plot | File | Description |
|------|------|-------------|
| Pseudobulk sample PCA | `pca_pseudobulk_samples.png` | VST-transformed pseudobulk counts, color=region_group, shape=Age |
| Sample distance heatmap | `sample_distance_heatmap.png` | VST Euclidean distance, annotated by region_group and Age |

### 10.2 DE Plots

| Plot | File | Description |
|------|------|-------------|
| All-gene volcano | `volcano_CA3_vs_CA1_all_genes.png` | log2FC vs -log10(padj), color: Up_CA3 / Down_CA1 / NS, label top 20 genes |
| Target-gene volcano | `volcano_CA3_vs_CA1_target_genes.png` | Same as above but only 231 target genes shown, labeled with ggrepel, colored by category |
| Target-gene heatmap | `target_gene_heatmap.png` | pheatmap/ComplexHeatmap: rows=target genes, cols=pseudobulk samples, split by category, annotated by region_group |
| Category summary barplot | `category_de_summary_barplot.png` | Per-category: n_up / n_down / n_notable, side-by-side bars |

### 10.3 Module Score Plots

| Plot | File | Description |
|------|------|-------------|
| Module score dotplot | `module_score_region_age_dotplot.png` | 12 modules × (CA1, CA3) × (Young, Middle, Old), dot=mean score, errorbar=SE |
| CA3-CA1 delta barplot | `module_score_CA3_minus_CA1_delta.png` | Per-module delta (CA3-CA1), bar height=mean delta, color=padj significance |
| Module score heatmap | `module_score_sample_heatmap.png` | GEMgroups × modules, split by region_group |

### 10.4 Coupling Plots

| Plot | File | Description |
|------|------|-------------|
| Coupling heatmap | `coupling_spearman_heatmap.png` | Spearman rho matrix, lower triangle: CA1, upper triangle: CA3, or two separate panels |
| Coupling scatter plots | `selected_coupling_scatterplots.png` | Multi-panel (patchwork): one scatter per pre-specified pair (§7.2), point=GEMgroup, color=Age |

### 10.5 Age Plots

| Plot | File | Description |
|------|------|-------------|
| Age coupling boxplot | `age_coupling_boxplot.png` | Module score delta (CA3-CA1) by Age, one panel per module, boxplot per Age |

### 10.6 Spatial Maps (≤12 total, Exploratory)

| Plot | File | Description |
|------|------|-------------|
| Module score spatial maps | `spatial_module_{Name}_GEMgroup_{G}.png` | Up to 6 spatial map panels, one module per panel, default slice=first GEMgroup with both CA1+CA3 |

### 10.7 Plot Conventions

- **Format**: PNG, 300 dpi, width 8-14, height 6-12 (adjusted per plot type)
- **Colors**: CA1 = `#2166ac`, CA3 = `#b2182b` (blue/red consistent with Phase 03)
- **Age colors**: Young = `#c6dbef`, Middle = `#4292c6`, Old = `#08306b` (author Script 7 palette)
- **PDF**: Optional, tryCatch fallback. PNG is the authoritative format.
- **Rplots.pdf**: Guard at script start (record pre-existence, MD5 via `digest::digest()`). Delete ONLY if created by this run. Never delete a pre-existing Rplots.pdf that this run didn't modify.
- **Saving**: `ggsave()` for ggplot2; `png()`/`dev.off()` for pheatmap; `Cairo::CairoPNG()` not required

---

## 11. Module I: Output Manifest

### 11.1 Script

`R/spatial/s11_ca1_ca3_de_module_coupling.R` (~E0-E13 following Phase 03/s10 patterns)

### 11.2 Data Output Directory

`data/processed/spatial/phase11_ca1_ca3_de_module_coupling/`

| # | File | Description |
|---|------|-------------|
| 1 | `pseudobulk_sample_manifest.csv` | Per-sample metadata: GEMgroup, region, Age, n_spots, has_pair, low_coverage |
| 2 | `design_matrix_audit.csv` | Design validation: rank, columns, sample mapping |
| 3 | `ca1_ca3_pseudobulk_counts.rds` | Pseudobulk integer matrix (genes × samples, sparse if possible) |
| 4 | `deseq2_all_genes_CA3_vs_CA1.csv` | Full DESeq2 results (all genes passing count filter) |
| 5 | `deseq2_all_genes_CA3_vs_CA1_shrunken.csv` | LFC-shrunken results (only if apeglm or normal shrinkage succeeds) |
| 6 | `deseq2_target_genes_CA3_vs_CA1.csv` | Target genes subset (231 found + 20 missing annotated) |
| 7 | `target_gene_effect_summary_by_category.csv` | Per-category DE summary |
| 8 | `missing_gene_symbol_rescue_audit.csv` | Rescue candidates for 20 missing genes |
| 9 | `module_gene_set_final_audit.csv` | Per-module: gene list, counts, exclusions, notes |
| 10 | `module_scores_spot_summary.csv` | Per-spot module scores (12 modules × 9,921 rows if computed) |
| 11 | `module_scores_sample_summary.csv` | Per-GEMgroup × region_group aggregated module scores + metadata |
| 12 | `module_score_delta_CA3_minus_CA1.csv` | Per-module delta with paired Wilcoxon results |
| 13 | `module_score_delta_by_age.csv` | Delta × GEMgroup × Age for age-stratified analysis |
| 14 | `module_coupling_spearman_ca1.csv` | Spearman rho/p-values for coupling pairs (CA1 context) |
| 15 | `module_coupling_spearman_ca3.csv` | Spearman rho/p-values (CA3 context) |
| 16 | `module_coupling_spearman_delta.csv` | Spearman rho/p-values (delta context) |
| 17 | `age_module_summary.csv` | Mean module score by Age × region_group |
| 18 | `age_region_delta_summary.csv` | CA3-CA1 delta by Age for each module |
| 19 | `phase11_validation_summary.txt` | Check results (PASS/WARNING/STOP) |
| 20 | `phase11_provenance.csv` | Run metadata: versions, timings, renv.lock MD5, packages installed |

### 11.3 Figure Output Directory

`figures/spatial/phase11_ca1_ca3_de_module_coupling/`

| # | File | Description |
|---|------|-------------|
| 1 | `pca_pseudobulk_samples.png` | PCA of VST-transformed pseudobulk counts |
| 2 | `sample_distance_heatmap.png` | Euclidean distance heatmap |
| 3 | `volcano_CA3_vs_CA1_all_genes.png` | All-gene volcano |
| 4 | `volcano_CA3_vs_CA1_target_genes.png` | Target-gene highlighted volcano |
| 5 | `target_gene_heatmap.png` | Target gene DE heatmap (rows=genes, cols=samples, split by category) |
| 6 | `category_de_summary_barplot.png` | Per-category DE summary |
| 7 | `module_score_region_age_dotplot.png` | Module scores by region × age |
| 8 | `module_score_CA3_minus_CA1_delta.png` | Delta barplot |
| 9 | `module_score_sample_heatmap.png` | Module score sample-level heatmap |
| 10 | `coupling_spearman_heatmap.png` | Coupling Spearman heatmap |
| 11 | `selected_coupling_scatterplots.png` | Top coupling scatterplots |
| 12 | `age_coupling_boxplot.png` | Age trend plots |
| 13-18 | `spatial_module_*.png` | Up to 6 spatial module maps (exploratory) |

---

## 12. Module J: Validation Checks (V01-V35)

| ID | Check | Criterion | Severity |
|----|-------|-----------|----------|
| V01 | Phase 10 validation PASS | 25 PASS, 0 FAIL | STOP if not |
| V02 | Hippo RDS loaded | 9,921 spots, 32,285 Spatial features | STOP if not |
| V03 | Region invariants preserved | CA1=4671, CA3=2321, DG=2364 | STOP if >5% off |
| V04 | GEMgroup × Age mapping valid | 16 GEMgroups, 3 Age levels | STOP if not |
| V05 | Pseudobulk aggregation complete | ≥ 1 sample per CA1 and CA3 (overall) | STOP if not |
| V06 | Matched GEMgroups exist | ≥ 2 GEMgroups with has_pair=TRUE | STOP if not |
| V07 | Design matrix full rank | qr()$rank == ncol() | STOP if not |
| V08 | DESeq2 converged | No errors, results() succeeds | STOP if not |
| V09 | All-gene DE results have rows | nrow ≈ nrow(counts(dds)) | WARNING if empty |
| V10 | Target gene DE table has rows | 231 found + 20 missing = 251 | WARNING if mismatch |
| V11 | Missing gene rescue audit created | 20 rows, rescue_candidates populated | WARNING if missing |
| V12 | Complex V documented | exclusion or rescue status clear | WARNING if ambiguous |
| V13 | Module scores computed | 11-13 modules with spot-level scores | WARNING if < 11 |
| V14 | Sample-level module scores exist | Per GEMgroup × region_group | WARNING if missing |
| V15 | CA3-CA1 delta computed | Per module, per GEMgroup | WARNING if missing |
| V16 | Coupling correlations computed | Per pair, per context (CA1/CA3/delta) | WARNING if missing |
| V17 | Coupling rho in valid range | All rho ∈ [-1, 1] | WARNING if out of range |
| V18 | No spots as replicates in statistics | DESeq2 uses pseudobulk; module comparisons use GEMgroup-level | STOP if violated |
| V19 | No full dense conversion | Provenance confirms per-gene extraction only | STOP if violated |
| V20 | No full Seurat object saved | No `saveRDS(hippo)`, `save.image()`, or `.RData`. Small pseudobulk count RDS (`ca1_ca3_pseudobulk_counts.rds`) is allowed and expected | PASS/FAIL |
| V21 | No WholeBrain/DG loaded | Single Hippo object only | PASS/FAIL |
| V22 | Rplots.pdf handled | Pre-existence recorded; deleted only if created during run | PASS/FAIL |
| V23 | All Age groups present in summaries | Young, Middle, Old all represented | WARNING if not |
| V24 | All output files in correct directory | Data in data/processed/spatial/phase11_*, figs in figures/spatial/phase11_* | WARNING if not |
| V25 | No enrichment run | No GO/KEGG/Reactome calls | PASS/FAIL |
| V26 | No biological interpretation language in phase11_validation_summary.txt | Only statistical descriptions | PASS/FAIL |
| V27 | No Tangram/Python invoked | No python/reticulate calls | PASS/FAIL |
| V28 | target_genes.csv NOT modified | MD5 unchanged from Phase 10 | PASS/FAIL |
| V29 | Module gene sets match Phase 10 audit | Verified per category | WARNING if mismatch |
| V30 | Package installations logged | Each install has: command, packageVersion(), snapshot() | WARNING if unlogged |
| V31 | renv.lock change documented | Before/after MD5 recorded | PASS/FAIL |
| V32 | Interaction model rank-checked before fitting | Only attempted if full rank | STOP if forced |
| V33 | Coupling p-values BH-adjusted or labeled unadjusted | Consistent labeling in output | WARNING if ambiguous |
| V34 | Missing gene rescue does NOT auto-add genes | Only candidates recorded, no automatic substitution | STOP if violated |
| V35 | Plan file re-read at execution start | Script verifies plan MD5 or logs plan version | WARNING if skipped |

---

## 13. Module K: Stop Conditions (SC1-SC18)

| ID | Condition | Action |
|----|-----------|--------|
| SC-1 | Phase 10 validation not PASS | STOP — run Phase 10 first |
| SC-2 | Hippo RDS not found or unreadable | STOP — object unavailable |
| SC-3 | Spatial assay missing or corrupt | STOP — no count data |
| SC-4 | CA1, CA3, Region, GEMgroup, or Age metadata missing | STOP — region framework not available |
| SC-5 | < 2 GEMgroups with both CA1 and CA3 | STOP — primary paired CA1-vs-CA3 DE cannot proceed with < 2 matched pairs. An unpaired sensitivity model `~ region_group` may be run only if explicitly labeled as a fallback/sensitivity analysis and must NOT be described as the primary paired result |
| SC-6 | Pseudobulk aggregation produces all-zero pseudobulk samples | STOP — data integrity issue |
| SC-7 | Design matrix not full rank | STOP — investigate and report |
| SC-8 | DESeq2 fails to converge or errors | STOP — investigate data structure |
| SC-9 | DESeq2 unavailable (package removed from renv) | STOP — reinstall via renv::install |
| SC-10 | target_genes.csv not readable | STOP — Phase 10 audit prerequisite |
| SC-11 | 231 exact-match genes vs Phase 10 audit mismatch (>5% different) | STOP — re-run Phase 10 first |
| SC-12 | User requests WholeBrain loading | STOP — requires separate plan |
| SC-13 | User requests enrichment analysis (GO/KEGG/Reactome) | STOP — Phase 12+ scope |
| SC-14 | User requests biological interpretation in Phase 11 output | STOP — premature; Phase 12+ |
| SC-15 | User requests Tangram or Python | STOP — requires separate plan |
| SC-16 | User requests Seurat object saving | STOP — against AGENTS.md |
| SC-17 | User requests auto-substitution of missing genes without review | STOP — against plan §3 |
| SC-18 | User requests cell-type deconvolution (RCTD, cell2location, etc.) | STOP — requires separate plan |

---

## 14. Module L: Boundaries

### Phase 11 DOES

- Load `seurat_Visium_Hippo_All.rds` as sole spatial object
- Audit and rescue missing gene symbols (20 genes, especially Atp5*)
- Build CA1+CA3 pseudobulk count matrix from Spatial raw UMI
- Run paired DESeq2: CA3 vs CA1, blocked by GEMgroup
- Annotate target gene DE results with Phase 10 categories
- Compute spot-level module scores for 11-13 gene modules (AddModuleScore)
- Aggregate to sample-level (GEMgroup × region_group) module scores
- Compute CA3-CA1 module score deltas (paired Wilcoxon)
- Compute Spearman coupling for 10 pre-specified module pairs
- Compute age-stratified summaries and CA3-CA1 delta by Age
- Generate all QC, DE, module, coupling, and age plots (PNG)
- Generate exploratory spatial maps for selected modules (≤6)
- Validate all outputs against V01-V35
- Log all package installations, renv.lock changes, run metadata
- Include DG and CA2 as context regions in summaries and plots

### Phase 11 DOES NOT

- Claim biological mechanism, causation, or disease relevance
- Run enrichment analysis (GO/KEGG/Reactome)
- Use spots as biological replicates in any statistical test
- Use IFNγ/Edge/NNS/ENS labels as region labels
- Run Tangram or Python code
- Load WholeBrain or DG RDS objects
- Modify `docs/target_genes.csv`
- Save full Seurat objects to .rds
- Generate Rplots.pdf (guard at start and end)
- Perform cell-type deconvolution on Visium spots
- Auto-generate marker gene lists from external databases
- Force interaction models with insufficient data
- Auto-substitute missing genes without user review
- Perform trajectory inference (Monocle3, Slingshot, scVelo)

---

## 15. Execution Outline (s11 Script Steps)

```
Phase Spatial-11: CA1 vs CA3 DE + Module/Coupling Analysis

E0: Setup
  - Re-read this plan file and record plan_version = "v2.0"
  - Rplots.pdf guard: record pre-existence, MD5 via digest::digest()
  - Git status at start
  - renv.lock MD5 before (reference: 139bcfaf7ce24b633c9f034380a15ef5)
  - Load packages: Seurat, DESeq2, Matrix, ggplot2, pheatmap, ggrepel, data.table, digest
  - Check optional packages: apeglm, ComplexHeatmap, matrixStats → install if beneficial, log if installed
  - Load seurat_Visium_Hippo_All.rds
  - Verify invariants: 9921 spots, CA1=4671, CA3=2321, DG=2364
  - Load Phase 10 gene audit files

E1: Missing Gene Rescue Audit
  - For each of 20 missing genes: substring search in rownames(hippo[["Spatial"]])
  - For Atp5*: grep("^Atp5", rownames) → catalog all Atp5* genes in object
  - For Ndufb1: grep("^Nduf", rownames) → catalog Nduf* family
  - For Mrps*: grep("^Mrps|^Mrpl", rownames) → catalog mitoribosomal genes
  - Save: missing_gene_symbol_rescue_audit.csv
  - Log: V11

E2: Pseudobulk Aggregation (CA1 + CA3 only)
  - Subset hippo to CA1 + CA3 spots (n ≈ 6,992)
  - Get Spatial raw counts (sparse)
  - t(counts_sparse) + aggregate.Matrix by spot_groupings = GEMgroup_Region, transpose back to genes × samples
  - Verify: nrow == nrow(counts_sparse), sample_id alignment, no dense conversion
  - Build pseudobulk count matrix (genes × samples)
  - Build sample manifest: GEMgroup, region_group, Age, n_spots, has_pair, low_coverage
  - Audit matched GEMgroups
  - Save: pseudobulk_sample_manifest.csv, ca1_ca3_pseudobulk_counts.rds
  - Log: V05, V06

E3: Design Matrix Validation
  - Build primary design: ~ GEMgroup + region_group
  - Check full rank via qr()
  - If fail: try ~ region_group (two-group fallback)
  - Count samples per cell
  - Save: design_matrix_audit.csv
  - Log: V07
  - Check: SC-7

E4: DESeq2 Execution
  - DESeqDataSetFromMatrix with primary design
  - Filter: rowSums(counts(dds)) >= 10
  - DESeq()
  - results(contrast = c("region_group", "CA3", "CA1"))
  - lfcShrink(type = "apeglm") or fallback "normal"
  - Save: deseq2_all_genes_CA3_vs_CA1.csv, deseq2_all_genes_CA3_vs_CA1_shrunken.csv
  - Log: V08, V09

E5: Target Gene DE Annotation
  - Merge all-gene results with Phase 10 gene audit
  - Add category annotations
  - Per-category summary: n_available, n_up, n_down, n_sig, mean_log2FC
  - Save: deseq2_target_genes_CA3_vs_CA1.csv, target_gene_effect_summary_by_category.csv
  - Log: V10, V12

E6: Module Score Computation
  - Define gene sets from Phase 10 categories (12 categories, excluding Complex V unless rescued)
  - Spot-level: Seurat::AddModuleScore (batch if needed for memory)
  - Save spot-level: module_scores_spot_summary.csv
  - Aggregate to sample-level: mean/median per GEMgroup × region_group
  - Save: module_scores_sample_summary.csv
  - Save: module_gene_set_final_audit.csv
  - Log: V13, V14, V29

E7: Module Score CA3-CA1 Comparison
  - Paired delta: CA3 - CA1 per GEMgroup (only has_pair=TRUE)
  - Wilcoxon signed-rank per module
  - BH adjustment across modules
  - Save: module_score_delta_CA3_minus_CA1.csv, module_score_delta_by_age.csv
  - Log: V15

E8: Spearman Coupling Analysis
  - For each pre-specified pair (§7.2):
    - Spearman rho in CA1-only, CA3-only, CA3-CA1 delta contexts
    - BH adjustment within each context
  - Save: module_coupling_spearman_ca1.csv, ca3.csv, delta.csv, combined_long.csv
  - Log: V16, V17, V33

E9: Age-Stratified Summaries
  - Module score by Age × region_group
  - Delta by Age (descriptive, no interaction model unless full rank)
  - Optional interaction model: only if design supports (≥2 GEMgroups/cell, full rank)
  - Save: age_module_summary.csv, age_region_delta_summary.csv
  - Log: V32

E10: Plots
  - PCA pseudobulk, sample distance heatmap
  - Volcano (all genes), volcano (target genes highlighted)
  - Target gene heatmap, category barplot
  - Module score dotplot, delta barplot, heatmap
  - Coupling heatmap, scatterplots
  - Age boxplot
  - Spatial maps (≤6, exploratory)
  - All to figures/spatial/phase11_ca1_ca3_de_module_coupling/

E11: Validation
  - Run checks V01-V35
  - Save: phase11_validation_summary.txt
  - Save: phase11_provenance.csv (run metadata, packages, renv.lock MD5, timings)
  - Check V22 (Rplots.pdf), V30 (packages), V31 (renv.lock)

E12: Cleanup
  - rm(hippo, counts, dds); gc()
  - Rplots.pdf guard end: verify existence, delete only if created by this run
  - Verify renv.lock MD5 (record change if packages installed)
  - Final exit code: 0 if PASS overall, 1 if any STOP condition hit
```

---

## 16. Key Design Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| GEMgroup 1 lacks CA3 spots | Medium | Low | Drop from paired analysis; two-group fallback if < 2 pairs remain |
| < 3 matched GEMgroups | Low | High | Switch to two-group ~ region_group; lose pairing but all samples kept |
| Complex V 0 genes (after rescue) | High (16 Atp5* genes) | Medium | Exclude with WARNING; 11 modules instead of 12 |
| Age interaction model rank-deficient | High | Low | Skip interaction; use descriptive age-stratified summaries only |
| apeglm unavailable | Low | Low | Use "normal" shrinkage; record WARNING |
| ComplexHeatmap install fails | Low | Low | Use pheatmap (already installed) |
| Module AddModuleScore memory pressure | Medium | Medium | Batch calls, gc() between modules |
| Pseudobulk samples with < 20 spots | Medium | Medium | Flag low_coverage; run sensitivity with/without these samples |

---

## 17. Default Decisions (Formerly "Open Questions")

All decisions below are defaults for Phase 11 execution. Any can be overridden by user instruction before execution begins.

| # | Question | Default Decision |
|---|----------|-----------------|
| D1 | **Age contrast**: All three ages or Old vs Young only? | All three ages (Young, Middle, Old) in all summaries, plots, and age-stratified analysis. Primary DE is CA3 vs CA1 (not age-stratified). Middle is preserved. |
| D2 | **GEMgroup 1**: If missing CA3, drop or keep? | Drop from PAIRED analysis (no mate). Keep in two-group models and in all sample-level summaries. |
| D3 | **Complex V rescue**: Attempt rescue and wait for user review, or exclude immediately? | Attempt rescue audit (§3). Document candidates. Exclude from module scores (11 modules). Open `user_action_required` for future resolution. |
| D4 | **Module score aggregation**: Mean or median? | Both. mean for dotplots and heatmaps (standard); median for robustness checks in validation. |
| D5 | **Coupling p-value adjustment**: BH or unadjusted? | BH adjustment within each coupling context (CA1, CA3, delta). Label exploratory in output notes. |
| D6 | **Spatial maps for modules**: Include or defer? | Include ≤6 spatial maps (selected modules). Default: mtDNA-encoded, Complex I, Mitoribosomal, RPL, RPS, Defense/Antioxidant. |
| D7 | **DG DE**: Run CA1 vs DG and CA3 vs DG DESeq2? | **Secondary**. Run only if primary CA1 vs CA3 is successful AND execution time permits. Otherwise defer to Phase 12. |
| D8 | **Interaction model**: Attempt Age:region_group interaction? | Skip by default (high risk of rank deficiency). Replace with descriptive age-stratified summaries and delta-by-Age boxplots. |
| D9 | **Output compression**: gzip large CSVs? | Yes, for CSVs > 10 MB (consistent with Phase 10). Use `gzip` or `data.table::fwrite(..., compress = "gzip")`. |
| D10 | **Package installation**: Install apeglm, ComplexHeatmap, etc.? | Install if beneficial. Always use renv protocol. Record in provenance. Do NOT install packages not listed in §2.4 without user approval. |
| D11 | **Coupling network plot**: Circular/chord diagram? | Skip by default. Heatmaps and scatterplots are sufficient for Phase 11. If requested: install `circlize` and add as conditional output. |
| D12 | **DG role**: Visualizations and summaries? | Include DG in all non-statistical summaries and plots as context. Mark as `context_region` (not primary comparison). Do NOT include DG pseudobulk in the CA1-vs-CA3 DESeq2. |

---

## 18. Execution Command (Expected)

```bash
Rscript R/spatial/s11_ca1_ca3_de_module_coupling.R
```

Execution phase MUST:
1. Re-read this plan file: `docs/spatial_phase11_ca1_ca3_de_module_coupling_plan.md`
2. Verify plan version matches expectations
3. If chat context conflicts with this plan, this plan wins

---

> **Plan status**: REVISED v2.1. No analysis outputs generated during plan revision. No packages installed, no R scripts executed, no objects loaded.
> **Next step**: User review of v2.0 revisions → approval → execution.
> **Execution prerequisite**: Phase 10 must be PASS. All required inputs must exist.
