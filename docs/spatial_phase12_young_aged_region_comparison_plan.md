# Phase Spatial-12: Young vs Aged Regional Comparison across CA1 / CA3 / DG

> **Plan version**: v1.3 (2026-06-04, revised — DESeq2 QC diagnostics: padj=NA + p-value histogram + alpha=0.05 + significance vs effect-size)
> **Status**: READ-ONLY PLAN — no analysis execution; only this plan file is created/modified; no analysis outputs generated
> **v1.3 changes**: Added Route A-i (padj=NA diagnostics per region), Route A-ii (p-value histogram QC per region), alpha=0.05 explicit throughout, significance vs effect-size flag distinction, new output files (D16, F10-F12), new validation checks (V36-V40), QC backfill plan reference. No analysis design changed.
> **v1.2 changes**: Complete rewrite of §0 with Bio-Skill Guidance Audit + Authoritative Reference Review.
> **v1.1 changes**: Route A Pairs→sample counts; CA3 power label unified; module score provenance; optional model language strengthened.
> **Handoff**: Phase 11/11b PASS_WITH_WARNINGS (11: 34 PASS+1 WARNING; 11b: 30 PASS+1 WARNING) → Phase 12 planning
> **Phase 11**: CA1-vs-CA3 all-age DESeq2: 31 pseudobulk samples, 17,970 genes, 5,921 sig, 12 modules, 30 coupling pairs
> **Phase 11b**: Age-stratified CA1-vs-CA3: Young=3 GEMgroups/25 target sig, Middle=4 GEMgroups/14 target sig, Old=8 GEMgroups/123 target sig, 104 concordant, 65 CA1-high, 31 old-enhanced, 24 age-flipped
> **bioskill**: No single `bioskill` directory. 8 bio-* skills inspected (§0.2). 3 DESeq2 skills confirm Phase 12 design; 5 spatial skills are Python/Squidpy — not directly applicable.
> **/ref-bio**: Triggered. Authority: This saved `.md` file is the sole execution authority. If chat context conflicts with this plan, this plan wins. Execution phase must re-read this file first.

---

## 0. ref-bio / Bio-Skill Guidance / Authoritative Reference Review

### 0.1 Task Classification

| Dimension | Detected |
|-----------|----------|
| Data modality | **spatial transcriptomics** (Visium, RDS-based, not raw Space Ranger) |
| Workflow stage | **differential expression** (region-stratified Old vs Young pseudobulk) + **module score** + **coupling** + **cross-region comparison** |
| Package/tool | **DESeq2**, **Seurat**, **ggplot2**, **pheatmap**, **corrplot**, **data.table** |
| Platform | **Visium** (RDS-based reproduction/approximation) |
| Topic | **pseudobulk_de**, **spatial_transcriptomics**, **differential_expression** |

### 0.2 Bio-Skill Guidance Audit

No single directory named `bioskill` exists. The `.opencode/skills/` directory contains many `bio-*` skills. The following were inspected for Phase 12 relevance:

| # | Skill Path | Status | Primary Tool | Key Guidance Used | Impact on Phase 12 Plan |
|---|-----------|--------|-------------|-------------------|------------------------|
| 1 | `.opencode/skills/bio-differential-expression-deseq2-basics/SKILL.md` | **Inspected** | DESeq2 | Two-group `~ condition` design; `relevel()` BEFORE `DESeq()`; `results(name=...)` or `contrast=c(...)` for explicit contrasts; `lfcShrink(coef=..., type='apeglm')` for visualization; pre-filtering `rowSums >= 10`; Wald vs LRT; `padj=NA` three causes (independent filtering, Cook's, all-zero); `resultsNames(dds)` trap; n≤3 caveat (Schurch 2016: 20-40% false negatives) | **Confirmed**: `~ Age` two-group design; `contrast=c("Age","Old","Young")`; apeglm shrinkage; pre-filtering `rowSums >= 10`; CA3 Young n=3 → `primary_with_caution_young_n3` label; padj=NA diagnostic in validation |
| 2 | `.opencode/skills/bio-differential-expression-de-results/SKILL.md` | **Inspected** | DESeq2 | `padj=NA` diagnostic table (independent filtering / Cook's / all-zero); TREAT vs post-hoc LFC filtering with FDR caveat; p-value histogram diagnostics (U-shape = hidden batch); IHW for +5-20% power; `lfcThreshold=` for magnitude claims; Schurch 2016 small-n replication reality | **Confirmed**: `padj<0.05` significance threshold; `\|log2FC\|>log2(1.5)` as notable (NOT significance claim — "post-hoc filtering does NOT control FDR for magnitude"); p-value histogram as QC diagnostic; n≤3 labeling |
| 3 | `.opencode/skills/bio-differential-expression-de-visualization/SKILL.md` | **Inspected** | DESeq2/ggplot2 | Volcano with shrunken LFC (apeglm): p-values unchanged, LFCs compressed; PCA on VST/rlog (never raw counts); sample distance heatmap; p-value histogram; `blind=TRUE` vs `FALSE` for VST; `ggrepel max.overlaps=Inf`; row-scaling trap in heatmaps; top-variable gene selection via `rowMads` | **Confirmed**: Volcano plots use apeglm-shrunken LFC; VST for PCA/sample distance; p-value histogram as QC; plot conventions carried forward |
| 4 | `.opencode/skills/bio-spatial-transcriptomics-spatial-statistics/SKILL.md` | **Inspected** (Python/squidpy) | Squidpy | Spatial autocorrelation (Moran's I), co-occurrence, neighborhood enrichment — all Python-based | **Not directly applicable** to Phase 12 R/Seurat pipeline; spatial statistics deferred to future phase |
| 5 | `.opencode/skills/bio-spatial-transcriptomics-spatial-preprocessing/SKILL.md` | **Inspected** (Python/squidpy) | Squidpy | QC, filtering, normalization for spatial data | **Not directly applicable** (Python tool); Phase 12 uses Seurat RDS metadata, not raw processing |
| 6 | `.opencode/skills/bio-spatial-transcriptomics-spatial-data-io/SKILL.md` | **Inspected** (Python/squidpy) | Squidpy | Loading Visium data from Space Ranger output | **Not directly applicable**; Phase 12 loads author-provided Seurat RDS, not raw Visium |
| 7 | `.opencode/skills/bio-spatial-transcriptomics-spatial-visualization/SKILL.md` | **Inspected** (Python/squidpy) | Squidpy | Spatial plots (tissue overlay) | **Not directly applicable** (Python); Phase 12 uses Seurat `SpatialFeaturePlot` or ggplot2 |
| 8 | `.opencode/skills/bio-workflows-spatial-pipeline/SKILL.md` | **Inspected** (Python/squidpy) | Squidpy | End-to-end spatial workflow (QC → clustering → statistics → visualization) | **Conceptual reference only** — QC checkpoint philosophy (after_loading, after_qc) adapted as Phase 12 E1 audit checkpoint. Actual pipeline is R/Seurat, not Python/Squidpy. |

**Key takeaways from bio-skill audit**:
- The three DESeq2-related skills (deseq2-basics, de-results, de-visualization) directly validate and confirm the Phase 12 DE design.
- All spatial skills are Python/Squidpy-based — not directly applicable but their conceptual framework (spot-level data, spatial autocorrelation) informs the "spots ≠ replicates" principle.
- No bio-skill contradicts any Phase 12 design decision. The `~ Age` two-group design, apeglm shrinkage, explicit contrast naming, pre-filtering, and small-n labeling are all explicitly recommended by the DESeq2 bio-skills.

### 0.3 Authoritative Reference Review

The following authoritative sources were accessed to validate the Phase 12 analysis design. URLs were fetched from the ref-bio reference catalog.

| # | Source Title | URL | Reviewed Status | Specific Section/Topic Consulted | Design Decision Affected |
|---|-------------|-----|-----------------|----------------------------------|------------------------|
| R1 | DESeq2 Official Vignette (Bioconductor) | https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html | **Reviewed** | "Standard workflow": `DESeqDataSetFromMatrix` → `DESeq()` → `results(contrast=...)` → `lfcShrink(type='apeglm')`; pre-filtering `rowSums >= 10`; `relevel()` BEFORE `DESeq()`; "Note on factor levels": alphabetical reference level trap; "p-values and adjusted p-values": BH default, independent filtering; "Multi-factor designs": two-group design `~ condition` | **Confirmed**: `~ Age` two-group design is the standard DESeq2 workflow. `contrast=c("Age","Old","Young")` for explicit direction control. DESeq2 vignette explicitly supports `contrast` argument for naming the comparison. Pre-filtering `rowSums >= 10` exactly matches vignette recommendation. |
| R2 | Seurat Spatial Vignette | https://satijalab.org/seurat/articles/spatial_vignette.html | **Reviewed** | Overview: Visium data stored as "Spatial" assay with images slot; `SpatialFeaturePlot` for expression on tissue; "At ~50um, spots from Visium assay will encompass expression profiles of multiple cells"; `SCTransform` normalization; `Load10X_Spatial` vs RDS loading; `GetTissueCoordinates` for spatial coordinates; "Working with multiple slices" → merge function | **Confirmed**: Spatial assay = "Spatial" with raw UMI counts; spots NOT single cells; Seurat's `AddModuleScore` works on Spatial assay; pseudobulk aggregation by GEMgroup (tissue section) is the correct replicate unit per the vignette's multi-slice design. Note: Phase 12 loads author RDS, not raw Space Ranger — consistent with Phase 05 RDS-based route. |
| R3 | OSTA (Orchestrating Spatial Transcriptomics Analysis) | https://bioconductor.org/books/release/OSTA/ | **Partially reviewed** (TOC inspected; full chapters not fetched) | Relevant chapters identified: Ch.10 "Quality control", Ch.26 "Normalization", Ch.29 "Feature selection & testing", Ch.30 "Feature-set signatures", Ch.31 "Spatial statistics", Ch.33 "Differential spatial patterns", Ch.34 "Differential colocalization" | **Confirmed**: OSTA provides Bioconductor-native spatial analysis framework. Ch.33 (Differential spatial patterns) would be the direct reference for spatial DE but is deferred to execution-phase review. OSTA's existence validates that spatial pseudobulk DE is a recognized analysis pattern. |
| R4 | Seurat Differential Expression Testing | https://satijalab.org/seurat/articles/differential_expression.html | **Failed to access** (HTTP 404) | Attempted both `/differential_expression.html` and `/de_vignette.html` | **No design impact**. Seurat's DE vignette is secondary to DESeq2 vignette for pseudobulk DE. The plan relies on DESeq2 official vignette (R1) as the primary DE authority. |
| R5 | Pseudobulk DE Guidance (BioC 2024) | Referenced via ref-bio `pseudobulk-de-guidance` source ID | **Not directly reviewed** (link-only; upstream URL not explicitly listed in references.link-only.yaml) | Pseudobulk aggregation pattern; GEMgroup-level replicate unit; `aggregate.Matrix` usage; Crowell 2020 pseudobulk recommendation | **Confirmed by bio-skill**: DESeq2-basics skill (§"Single-cell pseudobulk") explicitly states Crowell 2020 *Nat Commun* 11:6077 finding that pseudobulk avoids FDR inflation of cell-level DE. The Phase 12 GEMgroup-level pseudobulk design is the recommended pattern. |

**Additional sources consulted (not in ref-bio catalog)**:
- **limma User's Guide**: Not fetched. The bio-skill confirms limma/voom as a valid sensitivity alternative for two-group comparisons. Phase 12 Route F correctly labels it as sensitivity only.
- **edgeR User's Guide**: Not fetched. Same status as limma.

**Key design decisions confirmed by authoritative review**:

| Design Decision | Source(s) Confirming | Confidence |
|----------------|---------------------|-----------|
| Region-stratified `~ Age` two-group unpaired DESeq2 | R1 (DESeq2 vignette: "Standard workflow"), Bio-skill #1 (DESeq2-basics) | **High** |
| Old-vs-Young contrast: `contrast=c("Age","Old","Young")`, positive LFC = Old-high | R1 ("Standard workflow" contrast syntax), Bio-skill #1 (Wald vs LRT, contrast naming) | **High** |
| `relevel(Age, ref="Young")` BEFORE `DESeq()` | R1 ("Note on factor levels"), Bio-skill #1 | **High** |
| `lfcShrink(type='apeglm')` for visualization/ranking | R1 ("Log fold change shrinkage"), Bio-skills #1/#3 | **High** |
| `rowSums(counts(dds)) >= 10` pre-filtering | R1 ("Pre-filtering"), Bio-skill #1 | **High** |
| `padj < 0.05` significance; `\|log2FC\| > log2(1.5)` as notable (NOT FDR-controlled threshold) | R1 ("p-values and adjusted p-values"), Bio-skill #2 (TREAT vs post-hoc) | **High** |
| GEMgroup is biological replicate unit (pseudobulk), NOT Visium spot | R2 (Visium spots encompass multiple cells), Bio-skill #1 (Crowell 2020 pseudobulk), R5 | **High** |
| Spots for visualization/exploratory only, no spot-level inference | R2 (spot-level data), Bio-skills #4-#8 (spatial context) | **High** |
| VST for PCA/sample distance (visualization only, never DE input) | R1 ("Data transformations and visualization"), Bio-skill #3 | **High** |
| P-value histogram as QC diagnostic | Bio-skills #1/#2/#3 (p-value histogram diagnostics) | **High** |
| CA3 Young n=3 → `primary_with_caution_young_n3` label | Bio-skill #1 (Schurch 2016: n=3 misses 20-40% true positives), Bio-skill #2 (small-n replication reality) | **High** |
| limma/edgeR as sensitivity only, NOT primary | Bio-skills #1/#2 (DESeq2 as default; cross-tool concordance for validation only) | **Medium** |
| Interaction model `region_group:Age` SKIP by default | Bio-skill #1 (design matrix rank check, LRT for interactions) | **High** |

**No design decision was contradicted by any authoritative source reviewed.**

### 0.4 Project Context Source Priority (per AGENTS.md)

| Priority | Source | Role |
|----------|--------|------|
| 1 | This plan file | Execution authority |
| 2 | AGENTS.md | Project conventions, memory constraints, stop conditions |
| 3 | Phase 10/11/11b outputs | Gene audit, pseudobulk manifests, DE results, module scores |
| 4 | Ref-bio source IDs | Method validation, design assurance |
| 5 | Seurat/Hippo RDS metadata | Ground truth for region/age/GEMgroup |

---

## 1. Scientific Goals and Scope

### 1.1 Primary Analysis (PRIMARY)

| # | Question | Unit | Method |
|---|----------|------|--------|
| Q1 | **Old vs Young DEG within CA1?** | GEMgroup-level pseudobulk | Region-stratified DESeq2: `~ Age` within CA1 |
| Q2 | **Old vs Young DEG within CA3?** | GEMgroup-level pseudobulk | Region-stratified DESeq2: `~ Age` within CA3 |
| Q3 | **Old vs Young DEG within DG?** | GEMgroup-level pseudobulk | Region-stratified DESeq2: `~ Age` within DG |
| Q4 | **Which target genes show age-related changes in CA1/CA3/DG?** | 231 target genes × 3 regions | Merge DESeq2 results with Phase 10 gene audit |
| Q5 | **Are age effects concordant across regions or region-specific?** | Cross-region log2FC matrix | Direction consistency classification (Route B) |

### 1.2 Secondary Analysis (SECONDARY)

| # | Question | Method | Caveat |
|---|----------|--------|--------|
| Q6 | **Do module scores differ Old vs Young within each region?** | Wilcoxon/linear model per module per region | GEMgroup-level; small n per region-age cell |
| Q7 | **Do mitochondrial/ribosomal/translation/stress-response category effects differ by region?** | Per-category Old-vs-Young log2FC summary | Aggregates DESeq2 results; no new statistical tests |

### 1.3 Exploratory Analysis (EXPLORATORY)

| # | Question | Method | Caveat |
|---|----------|--------|--------|
| Q8 | **Does module-module coupling change Young vs Old within each region?** | Age-stratified Spearman coupling per region | n ≤ 8 per age-region cell; rho descriptive only |
| Q9 | **Region × Age interaction model?** | `~ region_group + Age + region_group:Age` | SKIP by default; confounded with GEMgroup; sensitivity only if requested |

### 1.4 Scope Boundaries

**Middle is NOT a primary contrast**. Young vs Old only. Middle can be preserved in manifests and summaries as context but will NOT be used for primary DESeq2. The rationale:
- The user's scientific question targets Young vs Old.
- Phase 11b already characterized CA3-vs-CA1 effects across all three ages.
- Adding Middle would triple the number of DESeq2 runs and complicate cross-region comparison.
- Middle data remains available for reference and can be added in a future phase if needed.

**DG definition**: DG = ML + GCL + Hilus (2,364 spots), consistent with Phase 06/10/11.

**CA2**: Context/QC only. Not a primary comparison region.

**Units**: GEMgroup-level pseudobulk for all inferential statistics. Spots for visualization/exploratory maps only.

**No spots as replicates**: All DESeq2 uses pseudobulk aggregation. No spot-level statistical testing.

---

## 2. Recommended Analysis Routes

### Route A: Region-Stratified Pseudobulk DESeq2 (PRIMARY)

**Method**: Run independent DESeq2 Old vs Young within each region.

| Region | Young GEMgroups (n_Young) | Old GEMgroups (n_Old) | total_samples | Design | Power label |
|--------|---------------------------|------------------------|--------------|--------|-------------|
| CA1 | G1-G4 (4) | G9-G16 (8) | 12 | `~ Age` | **primary_adequately_powered** |
| CA3 | G2-G4 (3; G1 has no CA3) | G9-G16 (8) | 11 | `~ Age` | **primary_with_caution_young_n3** |
| DG | G1-G4 (expected 4) | G9-G16 (expected 8) | expected 12 (TBD after inspection) | `~ Age` | **primary_adequately_powered** (expected, pending inspection) |

**Note**: This is an unpaired two-group comparison (`~ Age`). Age and GEMgroup are perfectly colinear, so `~ GEMgroup + Age` is rank-deficient and cannot be used. There are no "pairs" — each GEMgroup contributes one pseudobulk sample per region. This is NOT a paired design.

**Critical design decision — why `~ Age` not `~ GEMgroup + Age`**:

Age and GEMgroup are **perfectly colinear** in this experimental design:
- GEMgroups 1-4 = Young
- GEMgroups 5-8 = Middle (excluded from Phase 12)
- GEMgroups 9-16 = Old

A design `~ GEMgroup + Age` is rank-deficient because Age is a linear combination of GEMgroup. DESeq2 will drop the Age term. The correct model for region-stratified Old-vs-Young is `~ Age` (simple two-group comparison within each region). This is the identical approach used in Phase 11 Age Grouping's design validation.

**Paired alternative**: A formal paired design `~ GEMgroup + Age` is not possible when Age = f(GEMgroup). Phase 11b's age-stratified CA3-vs-CA1 used `~ GEMgroup + region_group` because region_group varies *within* each GEMgroup, making paired design possible. Here, Age varies *across* GEMgroups, so within-region Old-vs-Young is inherently an unpaired comparison.

**Note on GEMgroup 1**: GEMgroup 1 has CA1 spots but no CA3 spots. It CAN be included in CA1 Old-vs-Young and DG Old-vs-Young, but NOT in CA3 (no spots). This is accounted for in the sample counts above.

**Design validation per region** (MANDATORY before DESeq2):
1. Build `model.matrix(~ Age)` for that region's pseudobulk samples
2. Check `qr(model_matrix)$rank == ncol(model_matrix)` (should be 2: intercept + AgeOld)
3. Check Age has exactly 2 levels: Young and Old
4. Count samples per Age cell (≥ 3 per cell minimum)
5. Save `design_matrix_audit_region_age.csv`

**Contrast**:
```r
# For each region:
dds_region <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts_region,  # genes × GEMgroups in that region
  colData   = sample_metadata_region,
  design    = ~ Age
)
dds_region <- dds_region[rowSums(counts(dds_region)) >= 10, ]
dds_region <- DESeq(dds_region)
res_region <- results(dds_region, contrast = c("Age", "Old", "Young"), alpha = 0.05)
```

**log2FC direction**:
- **positive log2FC** → higher in **Old/Aged**
- **negative log2FC** → higher in **Young**
- This is consistent with the user's specified direction convention.

**LFC shrinkage**: apeglm (installed in Phase 11); fallback "normal" if unavailable.

**Power caveats**:
- CA1 Old-vs-Young: 4 Young GEMgroups, 8 Old GEMgroups. Adequate for detecting moderate-to-large effects.
- CA3 Old-vs-Young: 3 Young GEMgroups, 8 Old GEMgroups. Slightly underpowered in Young arm; label as `primary_with_caution_young_n3`.
- DG Old-vs-Young: TBD after inspection; expected 4 Young, 8 Old → adequate.

**Alpha (significance level)**: ALL `results()` calls MUST use `alpha = 0.05` explicitly. ALL `summary(res)` calls MUST use `summary(res, alpha = 0.05)`. Do NOT rely on DESeq2 default `alpha` values.

**Significance vs effect-size flag**:
- **`padj < 0.05`** is the **significance threshold**. Only `padj < 0.05` supports FDR-controlled claims of differential expression (i.e., "Old ≠ Young").
- **`|log2FC| > log2(1.5)`** (~0.585) is a **post-hoc notable / effect-size flag**. It does NOT control FDR for the magnitude claim. This is NOT a TREAT / `lfcThreshold=` test. It is a descriptive filter only — used for annotation, flagging, and selecting candidates for visualization. Genes meeting both criteria may be described as "significant AND notable" but MUST NOT be described as having "FDR-controlled effect > 1.5x". The distinction must be recorded in all output files, figure captions, and the analysis guide.

### Route A-i: DESeq2 padj=NA Diagnostics (MANDATORY after each region's DESeq2)

**Goal**: Diagnose and document `padj = NA` causes per region before proceeding to target gene or cross-region analysis.

**For each region (CA1, CA3, DG)**, after `results()`:

```r
res_df <- as.data.frame(res_region)
# padj=NA with pvalue present → independent filtering (inferred from result table)
n_padj_na_with_pvalue <- sum(is.na(res_df$padj) & !is.na(res_df$pvalue))
# pvalue=NA → likely all-zero counts in group, Cook's outlier, or other DESeq2 behavior
n_pvalue_na <- sum(is.na(res_df$pvalue))
# Total padj=NA
n_padj_na <- sum(is.na(res_df$padj))
```

**Record the following per region**:

| Column | Description |
|--------|-------------|
| region | CA1 / CA3 / DG |
| n_genes_total | Number of genes passing pre-filter (`rowSums >= 10`) |
| n_pvalue_na | Genes with `pvalue = NA` (all-zero, Cook's outlier, or other DESeq2 behavior) |
| n_padj_na | Total genes with `padj = NA` |
| n_padj_na_with_pvalue_present | Genes where `pvalue` is present but `padj = NA` (inferred independent filtering — labeled `inferred_from_result_table`) |
| n_target_genes_padj_na | Number of target genes with `padj = NA` in this region |
| target_genes_padj_na | Comma-separated list of target gene symbols with `padj = NA` |
| n_target_genes_pvalue_na | Number of target genes with `pvalue = NA` |
| diagnostic_note | Free-text: e.g., "n_padj_na_with_pvalue > 0 consistent with independent filtering; n_pvalue_na > 0 — possible all-zero or Cook's; not asserting specific cause without further inspection" |

**Important caveats**:
- `pvalue` present but `padj` missing is **consistent with** independent filtering (Bourgon 2010), but when diagnosed from the result table alone, the cause is **inferred, not confirmed**. Label as `inferred_from_result_table`.
- `pvalue` missing can be from: all-zero counts, Cook's distance outlier flagging, or other DESeq2 internal behavior. Do NOT assert a specific cause without inspecting `mcols(res)$description`, `cooksCutoff`, and the count matrix for individual genes.
- If target genes appear in the `padj=NA` list, record them for the user. Common examples: low-count transcription factors or genes with near-zero expression in one age group.

**Output**: `data/processed/spatial/phase12_young_aged_region_comparison/deseq2_padj_na_diagnostics_by_region.csv`

### Route A-ii: DESeq2 p-Value Histogram QC (MANDATORY per region)

**Goal**: Detect model misspecification or hidden batch effects before trusting any gene list.

**For each region (CA1, CA3, DG)**, generate a raw p-value histogram:

```r
library(ggplot2)
p_hist <- ggplot(res_df, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = 'steelblue', color = 'white') +
    labs(x = 'Raw p-value', y = 'Frequency',
         title = paste('P-value distribution:', region, 'Old vs Young'),
         subtitle = paste('n_tested =', sum(!is.na(res_df$pvalue)))) +
    theme_bw()
ggsave(file.path(out_figs, paste0('pvalue_histogram_', region, '.png')),
       p_hist, width = 8, height = 6, dpi = 300)
```

**Output files**: `figures/spatial/phase12_young_aged_region_comparison/pvalue_histogram_CA1.png`, `pvalue_histogram_CA3.png`, `pvalue_histogram_DG.png`

**These are QC diagnostics, NOT biological results.** The guide (§7.12) must explain how to read them.

**Pattern interpretation** (from bio-skill de-visualization / de-results):

| Shape | Meaning | Action |
|-------|---------|--------|
| Uniform + spike near 0 | Correctly specified model; enrichment near 0 = true DE genes | Proceed |
| U-shape (spikes at 0 AND 1) | **Anti-conservative** — hidden batch effect or unmodeled covariate | Review PCA; consider adding batch to design (if not colinear) |
| Depleted near 0, spike near 1 | **Conservative** — over-modeled or wrong dispersion | Simplify model; check `plotDispEsts(dds)` |
| Spike only at p = 1 | Discrete artifact from very-low-count genes | Pre-filter more aggressively |
| Bimodal / unusual shape | Possible data or model issue | Flag for manual review; record in validation |

**p-value histograms are NOT used for biological interpretation.** They inform QC status and model trustworthiness. Recording them in the validation summary is sufficient.

**Method**: Collect Old-vs-Young log2FC from each region and construct a cross-region comparison matrix.

**Target gene log2FC matrix**:
```
rows = 231 target genes
cols = CA1_OldvsYoung_log2FC, CA3_OldvsYoung_log2FC, DG_OldvsYoung_log2FC
```

**Direction consistency classification** (descriptive, not statistical):

| Pattern | Definition | Interpretation |
|---------|-----------|----------------|
| **Age-concordant across regions** | Same sign in CA1, CA3, DG | Consistent aging direction across hippocampus |
| **CA1-specific aging** | Significant (padj<0.05) AND |log2FC|>log2(1.5) only in CA1 | Region-specific aging in CA1 |
| **CA3-specific aging** | Significant AND notable only in CA3 | Region-specific aging in CA3 |
| **DG-specific aging** | Significant AND notable only in DG | Region-specific aging in DG |
| **Direction-flipped between regions** | Different signs in ≥2 regions | Potential region-type switching with age |
| **Old-concordant** | Positive log2FC in ≥2 regions (Old-high) | Consistent Old-upregulation across regions |
| **Young-concordant** | Negative log2FC in ≥2 regions (Young-high) | Consistent Young-upregulation across regions |

**Category-level effect matrix**:
```
rows = 13 categories
cols = CA1_OldvsYoung_median_log2FC, CA3_OldvsYoung_median_log2FC, DG_OldvsYoung_median_log2FC
```

### Route C: Target Gene / Mitochondrial Category Summary (PRIMARY for user question)

**Method**: For each of the 13 gene categories, compute per-region Old-vs-Young summary:

| Column | Description |
|--------|-------------|
| region | CA1 / CA3 / DG |
| n_genes_total | Genes in this category from user list |
| n_genes_available | Found in Spatial assay |
| n_genes_missing | Not found (0 for most categories) |
| n_sig_padj05 | Number with padj < 0.05 in Old-vs-Young |
| n_notable_abs_log2FC_gt_log2_1_5 | Number with |log2FC| > log2(1.5) |
| n_old_high | Number with log2FC > 0 |
| n_young_high | Number with log2FC < 0 |
| mean_log2FC | Mean of non-NA log2FC for available genes |
| median_log2FC | Median of non-NA log2FC |
| direction_balance | ratio: n_old_high / (n_old_high + n_young_high) |

**Missing genes retained**: 20 missing genes get NA values in all DE columns, flagged as `in_spatial = FALSE`. The 16 Atp5* Complex V genes are marked with the Complex V caveat carried forward from Phase 11.

### Route D: Module Score Regional Aging Analysis (SECONDARY)

**Method**: Compare module scores between Old and Young within each region.

**Statistical unit**: GEMgroup-level module scores. CA1/CA3 scores may be reused from Phase 11 `module_scores_sample_summary.csv` (31 samples, CA1+CA3 only). DG scores must be newly computed from the Hippo object using `Seurat::AddModuleScore` and then aggregated to GEMgroup level. **Module score provenance must be recorded** per §5.5: gene set source, AddModuleScore parameters, score source per sample (`phase11_reused` vs `phase12_recomputed`), and a comparability caveat if CA1/CA3 and DG scores come from different batches.

**Modules**: 12 modules (Complex V excluded), consistent with Phase 11:
mtDNA-encoded, Complex I, Complex III, Complex IV, TCA cycle, Mitoribosomal, Membrane/Translocator, Defense/Antioxidant, RPL, RPS, Eif, Eef

**Per region, per module**:
1. Subset to Young and Old GEMgroups in that region
2. **Wilcoxon rank-sum test** (unpaired, nonparametric): Old vs Young module scores
   - Rationale: n ≤ 8 per age group per region; normality not assumed; independent groups (different GEMgroups)
3. If n per group < 3: label as `blocked_insufficient_samples`; compute descriptive statistics only
4. Report: n_Young, n_Old, mean_Young, mean_Old, mean_Old_minus_Young, Wilcoxon V statistic, raw p-value, BH-adjusted p-value across modules within each region

**Alternative** (if normality is defensible per ref-bio guidance):
- Simple linear model: `module_score ~ Age`, t-test for AgeOld coefficient
- Only if Shapiro-Wilk p > 0.05 in both groups (unlikely with n=3-8)

**Output**: `module_age_effect_by_region.csv`

| Column | Description |
|--------|-------------|
| module | Module name |
| region | CA1 / CA3 / DG |
| n_Young | Number of Young GEMgroups with scores |
| n_Old | Number of Old GEMgroups with scores |
| mean_Young | Mean module score (Young) |
| mean_Old | Mean module score (Old) |
| mean_delta_Old_minus_Young | Old - Young difference |
| wilcoxon_statistic | W statistic |
| wilcoxon_pvalue | Raw p-value |
| wilcoxon_padj_region | BH-adjusted within region |
| test_label | "wilcoxon_rank_sum" or "blocked_insufficient_samples" |

### Route E: Coupling Analysis (EXPLORATORY)

**Method**: Compute Spearman correlation of GEMgroup-level module scores within each region × age combination.

**Sample sizes per context**:
| Region | Age | n_samples (CA1/CA3 from Phase 11; DG TBD) |
|--------|-----|------------------------------------------|
| CA1 | Young | 4 GEMgroups (G1-G4) |
| CA1 | Old | 8 GEMgroups (G9-G16) |
| CA3 | Young | 3 GEMgroups (G2-G4) |
| CA3 | Old | 8 GEMgroups (G9-G16) |
| DG | Young | TBD (~4 GEMgroups) |
| DG | Old | TBD (~8 GEMgroups) |

**Coupling pairs**: Same 10 pre-specified pairs as Phase 11 (RPLvsRPS, CIvsCIII, CIvsCIV, MRPLvsRPL, MRPLvsRPS, mtNUC, MITOvsCYTO, ROSvsETC, TCAvsETC, TRANSvsETC).

**Policy**:
- n ≤ 4: Spearman rho reported as descriptive only; p-values omitted
- n = 8: Spearman rho with raw p-value, labeled `exploratory_n8`
- No BH adjustment across pairs within small-n contexts
- No statistical comparison of coupling between ages or regions
- Output: `module_coupling_young_old_by_region.csv`

**Columns**:
| Column | Description |
|--------|-------------|
| pair_name | e.g., RPLvsRPS |
| region | CA1 / CA3 / DG |
| age | Young / Old |
| n_samples | Number of GEMgroups for cor |
| spearman_rho | Spearman correlation coefficient |
| raw_pvalue_exploratory_small_n | Raw p-value (NA if n ≤ 4) |
| small_n_flag | TRUE for all (n ≤ 8) |

### Route F: Optional Advanced Models (SENSITIVITY / EXPLORATORY only; Default SKIP)

**IMPORTANT**: All models in this route are **sensitivity or exploratory only**. They must NOT replace the primary DESeq2 region-stratified `~ Age` results in Route A. The primary DE tables and analysis guide must use DESeq2 results. limma/edgeR/interaction outputs, if generated, must be in separate files clearly labeled as `_sensitivity_` and must not be cited as primary evidence.

**limma/voom**:
- Alternative DE framework with empirical Bayes moderation
- Can handle `~ Age` design identically to DESeq2
- If used: label all outputs as `sensitivity_limma_voom`; compare concordance with DESeq2 results in a separate audit table (`limma_vs_deseq2_concordance.csv`)
- Default: SKIP. DESeq2 is the primary DE engine. limma/voom adds complexity without clear benefit for simple two-group comparison but can serve as sensitivity if requested.
- Do NOT present limma/voom results in the primary volcano plots or analysis guide as if interchangeable with DESeq2.

**edgeR**:
- Analogous to limma/voom: alternative NB GLM framework
- If used: label as `sensitivity_edger`; compare with DESeq2; do NOT replace primary results

**variancePartition / dream**:
- Mixed-model approach for repeated measures
- NOT applicable here: no repeated measures; GEMgroup and Age are colinear
- Default: SKIP (not appropriate for this design)

**Region × Age interaction model**:
```r
design = ~ region_group + Age + region_group:Age
```
- **Critical limitation**: Drops GEMgroup as blocking factor. GEMgroup and Age are colinear. All pseudobulk samples (across regions) are independent — no within-GEMgroup pairing across regions.
- Requires: ≥ 2 samples per region_group × Age cell
- Sample structure: CA1_Young=4, CA1_Old=8, CA3_Young=3, CA3_Old=8, DG_Young=TBD, DG_Old=TBD
- Interpretation: tests whether the Old-vs-Young effect differs across regions
- **Default decision**: SKIP. The region-stratified DESeq2 (Route A) + cross-region log2FC comparison (Route B) are more interpretable and don't sacrifice the simplicity of within-region analysis.
- If attempted: label as `sensitivity_interaction_no_GEMgroup_blocking` in a separate output file; do NOT present as primary; do NOT include interaction results in the analysis guide as primary evidence.

---

## 3. Package and Environment Strategy

### 3.1 Current State

| Package | Status | Version |
|---------|--------|---------|
| DESeq2 | Installed (renv) | (from renv) |
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

### 3.2 Package Policy

| Rule | Description |
|------|-------------|
| Install allowed | YES — use `renv::install()` + `renv::snapshot()` |
| Record protocol | Before/after renv.lock MD5, `packageVersion()`, log in `phase12_provenance.csv` |
| May install | `ComplexHeatmap` (if pheatmap insufficient for category-split heatmaps), `scales` (for color), `cowplot`/`patchwork` (for multi-panel) |
| Prohibited | STutility, Tangram, Python, clusterProfiler, enrichR, cell-type deconvolution tools, variancePartition, dream |

### 3.3 renv.lock Baseline

Current renv.lock MD5: `b2bf05038b440576855d043e137c5883` (post-Phase 11/11b with apeglm + corrplot installed). This is the execution starting point. Record before/after MD5.

---

## 4. Input Files

### 4.1 Required (STOP if missing)

| # | File | Use |
|---|------|-----|
| I1 | `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` | Source of Spatial raw counts for pseudobulk aggregation (CA1, CA3, DG) |
| I2 | `docs/target_genes.csv` | Target gene list (251 genes, 13 categories, GBK encoding) |
| I3 | `data/processed/spatial/phase10_target_gene_audit/target_gene_input_audit.csv` | Phase 10 gene audit (category mapping, mapping_status: 231 exact-match, 20 missing) |
| I4 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/pseudobulk_sample_manifest.csv` | GEMgroup-region-age mapping for CA1+CA3, paired status |
| I5 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_scores_sample_summary.csv` | Sample-level module scores for CA1+CA3 (reuse to avoid recomputing AddModuleScore for CA1/CA3) |
| I6 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv` | Module gene set definitions (12 modules, Complex V excluded) |

### 4.2 Recommended (WARNING if missing)

| # | File | Use |
|---|------|-----|
| I7 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds` | Existing pseudobulk counts for CA1+CA3 (32,285 genes × 31 samples). Can be subset by region + age instead of recomputing |
| I8 | `data/processed/spatial/phase11_age_grouping/deseq2_target_genes_CA3_vs_CA1_by_age.csv` | Phase 11b age-stratified target gene results (reference for cross-phase comparison) |
| I9 | `data/processed/spatial/phase11_age_grouping/age_group_deg_summary.csv` | Phase 11b DEG counts per age (reference) |
| I10 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/missing_gene_symbol_rescue_audit.csv` | Phase 11 missing gene rescue audit (Complex V candidates) |
| I11 | `data/processed/spatial/phase10_target_gene_audit/region_sample_summary.csv` | Region × GEMgroup spot counts (for DG sample audit) |

### 4.3 Reference (not for loading)

| # | File | Use |
|---|------|-----|
| I12 | `docs/spatial_phase11_ca1_ca3_de_module_coupling_plan.md` | Phase 11 base plan (design reference) |
| I13 | `docs/spatial_phase11_age_grouping_plan.md` | Phase 11b plan (design reference for age stratification) |
| I14 | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/phase11_final_report.txt` | Phase 11 caveats and boundaries |
| I15 | `data/processed/spatial/phase11_age_grouping/phase11_age_grouping_validation_summary.txt` | Phase 11b validation status |
| I16 | `data/processed/spatial/phase10_target_gene_audit/phase10_validation_summary.txt` | Phase 10 validation status |

---

## 5. Data Reconstruction Strategy

### 5.1 CA1 and CA3 Pseudobulk

**Preferred**: Reuse existing `ca1_ca3_pseudobulk_counts.rds` from Phase 11 (32,285 genes × 31 samples). Subset by `Age %in% c("Young", "Old")` AND by `region_group`. This avoids reloading the full Seurat object if existing pseudobulk is valid.

**Fallback**: Rebuild from `seurat_Visium_Hippo_All.rds` using the same `Matrix.utils::aggregate.Matrix()` pattern established in Phase 11.

### 5.2 DG Pseudobulk (MUST build fresh)

DG is NOT in the Phase 11 pseudobulk counts (which only include CA1+CA3). The Phase 12 execution MUST:

1. Load `seurat_Visium_Hippo_All.rds`
2. Subset spots to DG = ML + GCL + Hilus AND Age %in% c("Young", "Old")
3. Extract Spatial raw counts (sparse) for DG spots
4. Aggregate by GEMgroup using `aggregate.Matrix()`: aggregate spots × genes by GEMgroup
5. Verify: genes × GEMgroups orientation, sample ID alignment
6. Save: `dg_pseudobulk_counts.rds`

**No full dense conversion**: All aggregation operates on sparse dgCMatrix. Per-gene numeric extraction is allowed. Full `as.matrix()` is prohibited.

**No full Seurat object saved**: Only pseudobulk count matrix and sample manifest are saved.

### 5.3 DG Sample Manifest

Build from Hippo metadata:
- Which GEMgroups have DG spots?
- Spot counts per GEMgroup × Age
- Young GEMgroups: G1-G4; Old GEMgroups: G9-G16
- Expected: all 12 GEMgroups (G1-G4 + G9-G16) have DG spots
- Record n_spots per GEMgroup; flag low_coverage if n_spots < 20

### 5.4 DG Module Scores

DG module scores are NOT in `module_scores_sample_summary.csv`. The Phase 12 execution MUST:
1. Subset Hippo object to DG spots
2. Run `Seurat::AddModuleScore()` with the same 12 gene sets from Phase 11 `module_gene_set_final_audit.csv`
3. Aggregate spot-level scores to GEMgroup-level (mean and median)
4. Save: `module_scores_dg_spot_summary.csv` and merge into an extended module scores dataset

**Memory**: Subset Hippo to DG (~2,364 spots) before AddModuleScore. gc() between module computation steps.

### 5.5 Module Score Provenance (MANDATORY)

CA1/CA3 module scores **may be reused** from Phase 11 `module_scores_sample_summary.csv`. DG module scores **must be newly computed** from the Hippo RDS. Because CA1/CA3 and DG scores may come from different computation batches, the execution phase MUST record:

1. **Gene sets used**: file path and MD5 of `module_gene_set_final_audit.csv`
2. **AddModuleScore parameters**: `ctrl`, `nbin`, `name` args used for each subset/run
3. **Source of CA1/CA3 scores**: either `"reused_from_phase11_module_scores_sample_summary.csv"` with MD5 of that file, or `"recomputed_from_hippo_rds"` if recomputed
4. **Source of DG scores**: `"recomputed_from_hippo_rds_DG_subset"`
5. **Comparability caveat**: If CA1/CA3 scores are from Phase 11 and DG scores are newly computed, record in `phase12_provenance.csv` and all module output files that CA1/CA3 scores come from a different computation batch (Phase 11, with potentially different Seurat session state, random seed, etc.) than DG scores (Phase 12). Label module score outputs with column `score_source` = `"phase11_reused"` or `"phase12_recomputed"` per sample.

**Alternative — recompute all together** (for full consistency): 
The execution phase may optionally recompute CA1, CA3, and DG module scores in a single `AddModuleScore` run on the full Hippo object (all 9,921 spots). This guarantees consistent computation parameters, random seed, and Seurat version. However:
- Memory cost: full Hippo object ~9921 spots may be significant.
- Time cost: 12 `AddModuleScore` calls on 9,921 spots may take 5-15 minutes.
- **Default**: reuse CA1/CA3 from Phase 11 + compute DG fresh (cost-effective). If the provenance audit reveals meaningful discrepancies or the user requests full consistency, fall back to recompute-all.

---

## 6. Output Design

### 6.1 Output Directories

```
data/processed/spatial/phase12_young_aged_region_comparison/
figures/spatial/phase12_young_aged_region_comparison/
```

### 6.2 Data Outputs

| # | File | Description |
|---|------|-------------|
| D1 | `region_age_pseudobulk_sample_manifest.csv` | Per-region × age sample metadata: GEMgroup, region_group, Age, n_spots, age_stratum (Young/Old) |
| D2 | `design_matrix_audit_region_age.csv` | Per-region design rank check, model_type, sample counts per Age cell |
| D3 | `deseq2_all_genes_Old_vs_Young_CA1.csv` | Full DESeq2 results for CA1: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj. With apeglm-shrunken LFC column if available. Columns: `model_type` = `two_group_unpaired`, `power_label` = `primary_adequately_powered`, `region` = `CA1` |
| D4 | `deseq2_all_genes_Old_vs_Young_CA3.csv` | Full DESeq2 results for CA3. Columns: `model_type`, `power_label` = `primary_with_caution_young_n3`, `region` = `CA3` |
| D5 | `deseq2_all_genes_Old_vs_Young_DG.csv` | Full DESeq2 results for DG. Columns: `model_type`, `power_label`, `region` = `DG` |
| D6 | `deseq2_target_genes_Old_vs_Young_by_region.csv` | 251 target genes (231 found + 20 missing) × 3 regions: log2FC, padj, baseMean, direction, lfc_notable, padj_sig per region. Columns: gene_symbol, category, category_english, mapping_status, log2FoldChange_{region}, padj_{region}, direction_{region}, sig_{region}, notable_{region}, in_spatial |
| D7 | `target_gene_log2fc_matrix_by_region.csv` | 231 genes × 3 region columns (wide format log2FC matrix: CA1_OldvsYoung, CA3_OldvsYoung, DG_OldvsYoung). With direction consistency classification: pattern_class, notes |
| D8 | `target_gene_category_effect_by_region.csv` | Per-category × per-region: n_available, n_sig, n_old_high, n_young_high, mean_log2FC, median_log2FC, direction_balance, power_label |
| D9 | `region_age_deg_summary.csv` | Per-region DEG counts: n_tested, n_sig, n_old_high, n_young_high, n_target_sig, n_target_old_high, n_target_young_high, power_label |
| D10 | `region_age_power_flags.csv` | Per-region power assessment: n_Young, n_Old, n_spots_Young, n_spots_Old, DESeq2 convergence, caveats |
| D11 | `module_age_effect_by_region.csv` | Per-module × per-region: n_Young, n_Old, mean_Young, mean_Old, mean_delta_Old_minus_Young, Wilcoxon W, p-value, padj within region |
| D12 | `module_coupling_young_old_by_region.csv` | Per-coupling-pair × per-region × per-age: n_samples, spearman_rho, small_n_flag. No p-values by default for n ≤ 4; raw p-value labeled `exploratory_small_n` for n=8 |
| D13 | `optional_region_age_interaction_audit.csv` | If interaction model attempted: rank check, convergence, caveats. Default: "SKIPPED_BY_DEFAULT" |
| D14 | `phase12_validation_summary.txt` | Check results (PASS/WARNING/STOP) |
| D15 | `phase12_provenance.csv` | Run metadata: versions, timings, renv.lock MD5 before/after, packages installed, git status |
| D16 | `deseq2_padj_na_diagnostics_by_region.csv` | Per-region padj=NA diagnostics: n_genes_total, n_pvalue_na, n_padj_na, n_padj_na_with_pvalue_present, target genes with NA, diagnostic_note |

### 6.3 Figure Outputs

| # | File | Description |
|---|------|-------------|
| F1 | `volcano_Old_vs_Young_CA1.png` | CA1 volcano: x=log2FC, y=-log10(padj), target genes highlighted by category, top 20 genes labeled. Title: "CA1: Old vs Young" |
| F2 | `volcano_Old_vs_Young_CA3.png` | CA3 volcano: as above. Title includes "CA3: Old vs Young (Young n=3, caution)" |
| F3 | `volcano_Old_vs_Young_DG.png` | DG volcano: as above. Title: "DG: Old vs Young" |
| F4 | `target_gene_log2fc_heatmap_by_region.png` | Heatmap: rows = 231 target genes (split by category), columns = 3 regions (CA1, CA3, DG), fill = log2FC Old-vs-Young, cells annotated * for padj<0.05 |
| F5 | `target_gene_category_effect_heatmap_by_region.png` | Heatmap: rows = 13 categories, columns = 3 regions, fill = median_log2FC |
| F6 | `deg_count_by_region_barplot.png` | Barplot: n_up_Old, n_down_Young per region, split all-genes/target-genes |
| F7 | `module_age_effect_by_region_boxplot.png` | 12 panels (one per module) × 3 regions (color: CA1=#2166ac, CA3=#b2182b, DG=#4daf4a), boxplot module score by Age (Young vs Old) |
| F8 | `selected_target_gene_age_effect_dotplot.png` | Dotplot: rows = selected genes (top 30 by max |log2FC| across regions), x = region, y = log2FC Old-vs-Young, color = direction, size = -log10(padj) |
| F9 | `coupling_spearman_young_old_heatmap.png` | (Exploratory, only if data supports) Heatmap: coupling rho for Old vs Young, faceted by region |
| F10 | `pvalue_histogram_CA1.png` | **QC diagnostic** (not biological result): raw p-value histogram for CA1 Old-vs-Young DESeq2 |
| F11 | `pvalue_histogram_CA3.png` | **QC diagnostic**: raw p-value histogram for CA3 Old-vs-Young DESeq2 |
| F12 | `pvalue_histogram_DG.png` | **QC diagnostic**: raw p-value histogram for DG Old-vs-Young DESeq2 |

### 6.4 Documentation Output

| # | File | Description |
|---|------|-------------|
| G1 | `docs/spatial_phase12_young_aged_region_comparison_guide.md` | User-facing analysis guide (generated during execution) |

---

## 7. Analysis Guide Template (`spatial_phase12_young_aged_region_comparison_guide.md`)

The execution phase MUST generate this guide. It must contain:

### 7.1 What Each Analysis Module Asks

| Module | Question | Input | Output |
|--------|----------|-------|--------|
| Region-stratified DESeq2 (CA1) | Old vs Young DEG in CA1? | CA1 pseudobulk counts (Young+Old only) | `deseq2_all_genes_Old_vs_Young_CA1.csv` |
| Region-stratified DESeq2 (CA3) | Old vs Young DEG in CA3? | CA3 pseudobulk counts (Young+Old only) | `deseq2_all_genes_Old_vs_Young_CA3.csv` |
| Region-stratified DESeq2 (DG) | Old vs Young DEG in DG? | DG pseudobulk counts (must build fresh) | `deseq2_all_genes_Old_vs_Young_DG.csv` |
| Target gene by region | Which target genes change with age in each region? | DESeq2 results + Phase 10 audit | `deseq2_target_genes_Old_vs_Young_by_region.csv` |
| Cross-region log2FC matrix | Are age effects region-specific or concordant? | 3 region-specific DESeq2 results | `target_gene_log2fc_matrix_by_region.csv` |
| Category by region | Do whole categories show age effects by region? | Cross-region matrix + categories | `target_gene_category_effect_by_region.csv` |
| Module delta by region | Do module scores change with age differently by region? | GEMgroup-level module scores (extended with DG) | `module_age_effect_by_region.csv` |
| Coupling by region and age | Does coupling differ Young vs Old within regions? | Module scores per GEMgroup × region × age | `module_coupling_young_old_by_region.csv` (exploratory) |

### 7.2 How to Read log2FC

- **`log2FoldChange > 0`** → higher expression in **Old/Aged**
- **`log2FoldChange < 0`** → higher expression in **Young**
- `baseMean` → mean of normalized counts across all pseudobulk samples in that region
- `padj` → Benjamini-Hochberg adjusted p-value (within that region's DESeq2 run only)

### 7.3 Thresholds

| Term | Definition | FDR control? |
|------|-----------|-------------|
| `padj_sig` | `padj < 0.05` | **YES** — BH-adjusted, FDR-controlled for "Old ≠ Young" |
| `log2FC_notable` | `|log2FC| > log2(1.5)` ≈ 0.585 | **NO** — post-hoc flag only. Does NOT control FDR for "\|effect\| > 1.5x". NOT equivalent to TREAT / `lfcThreshold=`. |
| `Old_high` | log2FC > 0 (direction label only; NOT necessarily significant) | N/A |
| `Young_high` | log2FC < 0 (direction label only; NOT necessarily significant) | N/A |
| `alpha` | 0.05 for `results(alpha=0.05)` and `summary(res, alpha=0.05)` | Used for independent filtering optimization; `summary()` default `alpha` must NOT be relied upon |

**IMPORTANT**: "Significant AND notable" (`padj < 0.05` AND `|log2FC| > log2(1.5)`) must be described as "genes with significant differential expression AND effect size above the 1.5-fold notable flag" — NOT as "genes with FDR-controlled effect > 1.5x". The latter requires TREAT (`lfcThreshold=`) or `glmTreat()`. If a reviewer challenges the magnitude claim, re-run with `results(dds, lfcThreshold = log2(1.5), altHypothesis = 'greaterAbs')` and report those results separately.

### 7.4 Why Spots Are Not Replicates

Visium spots from the same tissue section (GEMgroup) share the same animal and tissue processing. Within-section spots have correlated expression due to spatial proximity. Treating individual spots as replicates inflates false positive rates dramatically. The GEMgroup (tissue section) is the biological replicate unit. All DESeq2 analysis uses pseudobulk aggregation by GEMgroup.

### 7.5 Why Middle Is Not a Primary Contrast

Phase 12 targets the user's Young vs Old scientific question. Middle was fully characterized in Phase 11b (age-stratified CA3-vs-CA1) and can be revisited in a future phase if needed. Excluding Middle simplifies the comparison to 2 ages × 3 regions = 6 analysis cells instead of 3 ages × 3 regions = 9 cells, and reduces the number of DESeq2 runs from 9 to 3.

### 7.6 Region-Stratified DE vs Region×Age Interaction Model

| Aspect | Region-stratified (Route A+B) | Interaction (Route F) |
|--------|---------------------------|----------------------|
| Design per region | `~ Age` within each region | `~ region_group + Age + region_group:Age` across all samples |
| GEMgroup blocking | Not applicable (Age = f(GEMgroup)) | Dropped; GEMgroup not in model |
| Replicate unit | GEMgroup (preserved) | Pseudobulk sample (same) |
| Interpretation | "Old vs Young in CA1", "Old vs Young in CA3", "Old vs Young in DG" — simple | "Does the Old-vs-Young effect differ by region?" — requires interaction contrast |
| Sample handling | Independent DESeq2 per region; each region's dispersion estimated from its own samples | Pooled dispersion across all regions (higher total n, but samples are not independent across regions) |
| Default choice | **PRIMARY** | **SKIP by default** |

### 7.7 What Can Be Used for Downstream Biological Interpretation

| Result type | Status | For interpretation? |
|-------------|--------|---------------------|
| CA1 Old-vs-Young DESeq2 (padj<0.05, |log2FC|>log2(1.5)) | PRIMARY, 4 Young + 8 Old | YES — candidate genes for CA1 aging |
| CA3 Old-vs-Young DESeq2 (padj<0.05, \|log2FC\|>log2(1.5)) | PRIMARY, 3 Young + 8 Old (Young arm n=3) | YES but with caution for small Young n |
| DG Old-vs-Young DESeq2 (padj<0.05, \|log2FC\|>log2(1.5)) | PRIMARY, expected 4 Young + 8 Old | YES — candidate genes for DG aging |
| Cross-region log2FC direction consistency | DESCRIPTIVE | YES — direction patterns across regions are informative regardless of significance |
| Module score Old-vs-Young within region | SECONDARY, small n (3-8 per age) | DESCRIPTIVE — direction and magnitude trends |
| Age-stratified coupling | EXPLORATORY, n=3-8 | NOT recommended for interpretation; n too small |
| Interaction model (Route F) | EXPLORATORY, default SKIP | NOT recommended; confounded with GEMgroup |

### 7.8 Complex V / Missing Genes Caveat

16 Atp5* genes (Complex V) are absent from the Spatial assay. Phase 11 rescue audit found substring candidates but no user-confirmed replacements. Complex V is **excluded from module scores** (12 modules, not 13). Any statement about "mitochondrial gene expression" in these results EXCLUDES Complex V (ATP synthase). 20 missing genes in total (16 Atp5*, 1 Ndufb1, 3 Mrps) are retained with NA values in all DE columns. Their `in_spatial` flag is `FALSE`.

### 7.9 DG = ML + GCL + Hilus

DG is a derived region combining three subregions from the Hippo object:
- ML (Molecular Layer): 1,480 spots
- GCL (Granule Cell Layer): 645 spots
- Hilus: 239 spots
- Total DG: 2,364 spots

All DG analyses use this combined definition, consistent with Phase 06/10/11. Subregion-level analysis is deferred to future phases unless explicitly requested.

### 7.10 CA2 Role

CA2 (565 spots) is included as context/QC in summaries but is NOT a primary comparison region. CA2-specific Old-vs-Young or CA2-vs-other analyses are deferred unless explicitly requested.

### 7.11 How to Move from Phase 12 to CA1/CA3/DG Mitochondrial Interpretation

Phase 12 provides the region-aware Old-vs-Young framework. From this, the user can:
1. Select candidate target genes with region-specific or region-concordant aging effects
2. Identify mitochondrial/ribosomal categories that show region-specific aging
3. Compare module score aging patterns across regions
4. Formulate hypotheses about:
   - Why certain genes age differently in CA1 vs CA3 vs DG
   - Whether mitochondrial dysfunction markers are region-specific or pan-hippocampal
   - Whether ribosomal/translation module aging differs by region

**Do NOT** claim mechanism or causation from Phase 12 results. This phase identifies **patterns** for hypothesis generation. Biological interpretation requires additional evidence (literature, orthogonal validation, functional experiments).

### 7.12 Understanding padj=NA and p-Value Histogram QC

#### padj=NA Diagnostics

The file `deseq2_padj_na_diagnostics_by_region.csv` records `padj=NA` counts per region. There are three distinct causes:

| Scenario | In result table | Likely cause | Action |
|----------|----------------|--------------|--------|
| `pvalue` present, `padj = NA` | baseMean below data-driven threshold | Independent filtering (inferred) | Gene was excluded from FDR adjustment to maximize power. If a target gene of interest is here, record and note. Can re-run `results(dds, independentFiltering = FALSE)` to recover if needed. |
| `pvalue = NA`, `padj = NA`, all-zero in one group | baseMean very low | Insufficient counts to test | Expected; nothing to do. |
| `pvalue = NA`, `padj = NA`, nonzero counts | One sample flagged as outlier | Cook's distance outlier | If a target gene of interest is here, inspect with `plotCounts(dds, ...)`. Can re-run `results(dds, cooksCutoff = FALSE)`. |

These causes are **inferred from the result table alone**. Do NOT claim definitive causality without inspecting `mcols(res)$description`, `cooksCutoff`, or running `independentFilteringResults(res)`. The diagnostics table is for awareness and QC — not for mechanistic claims.

Target genes with `padj=NA` are recorded by name. If a user-important gene (e.g., a mitochondrial gene of interest) appears in this list, the gene was either filtered or flagged as an outlier — the user should be informed rather than silently dropping it.

#### p-Value Histogram QC

The files `pvalue_histogram_CA1.png`, `pvalue_histogram_CA3.png`, `pvalue_histogram_DG.png` are **QC diagnostics** — not biological results. They show the distribution of raw (unadjusted) p-values for each region.

**How to read**:
- **Uniform + near-zero spike**: Healthy. Null genes are uniform; true DE genes cluster near 0. Proceed.
- **Enrichment near 0 (large spike)**: Many DE genes. Expected with large effect sizes or high power.
- **U-shape (spikes at both 0 and 1)**: Anti-conservative. Suggests unmodeled batch/covariate. Review PCA; consider adding batch terms if not colinear with Age.
- **Depleted near 0 / spike near 1**: Conservative. Over-modeled or dispersion trend wrong. Check `plotDispEsts(dds)`.
- **Spike only at p = 1**: Low-count artifact. Pre-filter more aggressively.

**These histograms are NOT used for biological interpretation.** They inform QC status only. Recording them as PASS/WARNING in the validation summary is sufficient.

### 7.13 Relationship to Phase 03/11/11b DESeq2 QC Backfill

Phase 03 (DG DESeq2), Phase 11 (CA1-vs-CA3 all-age DESeq2), and Phase 11b (age-stratified CA1-vs-CA3 DESeq2) were run without padj=NA diagnostics or p-value histogram QC. A separate backfill plan, `docs/spatial_deseq2_qc_backfill_plan.md`, tracks these missing diagnostics for future completion. Phase 12 includes padj=NA diagnostics and p-value histogram QC **natively** from the start, so no Phase 12 backfill will be needed.

---

## 8. Validation Checks (V01-V40)

| ID | Check | Criterion | Severity |
|----|-------|-----------|----------|
| V01 | All required input files exist | I1-I6 all readable | STOP if missing |
| V02 | Hippo RDS loadable | 9,921 spots, CA1=4671, CA3=2321, DG=2364 | STOP if not |
| V03 | Region invariants preserved | Counts within 5% of Phase 06 | STOP if >5% off |
| V04 | Age groups: Young and Old present | At least 1 sample per region × age cell | STOP if missing |
| V05 | Middle excluded from primary contrast | Only Young/Old in DESeq2; Middle absent from colData for DESeq2 | PASS/FAIL |
| V06 | DG definition correct | DG = ML + GCL + Hilus, 2,364 spots | STOP if mismatch |
| V07 | CA1 pseudobulk: 4 Young + 8 Old = 12 samples | Count from manifest | WARNING if mismatch |
| V08 | CA3 pseudobulk: 3 Young + 8 Old = 11 samples (G1 excluded) | G1 has no CA3 spots | WARNING if mismatch |
| V09 | DG pseudobulk: GEMgroup count TBD | Built from Hippo metadata | WARNING if < 3 per Age |
| V10 | Design matrix full rank per region | qr()$rank == ncol() for each region's `~ Age` design | STOP if rank-deficient |
| V11 | DESeq2 converged per region | No errors in DESeq() | STOP if fails |
| V12 | DESeq2 output has rows | nrow ≈ nrow(pseudobulk_counts) - filtered | WARNING if empty |
| V13 | Target genes merged correctly | 231 found + 20 missing = 251 rows in target gene by-region table | WARNING if mismatch |
| V14 | Missing genes retained with NA DE values | 20 genes with in_spatial=FALSE, DE columns NA | WARNING if missing |
| V15 | Cross-region log2FC matrix complete | 231 rows, no all-NA rows | WARNING if unexpected NAs |
| V16 | Module scores available for CA1, CA3, DG | All 12 modules, all GEMgroups with scores | WARNING if missing |
| V17 | Module age effect computed per region | Wilcoxon or descriptive statistics per module × region | WARNING if missing |
| V18 | Coupling computed per region × age | Spearman rho for 10 pairs × 3 regions × 2 ages | WARNING if missing |
| V19 | Coupling rho in valid range | All rho ∈ [-1, 1] | WARNING if out of range |
| V20 | No spots as replicates | All DESeq2 uses GEMgroup-level pseudobulk | STOP if violated |
| V21 | No full spot-level dense conversion | Sparse aggregation only; per-gene numeric extraction allowed | STOP if violated |
| V22 | No full Seurat object saved | No `saveRDS(hippo)`, `save.image()`, or `.RData`. Pseudobulk count RDS allowed | PASS/FAIL |
| V23 | No WholeBrain loaded | Single Hippo object only | PASS/FAIL |
| V24 | No enrichment run | No GO/KEGG/Reactome calls | PASS/FAIL |
| V25 | No Tangram/Python/STutility/cell-type deconvolution | No python/reticulate/STutility calls | PASS/FAIL |
| V26 | No biological interpretation in validation summary | Statistical descriptions only | PASS/FAIL |
| V27 | target_genes.csv NOT modified | MD5 unchanged from Phase 10 run | PASS/FAIL |
| V28 | All output files in correct directories | Data in `data/processed/spatial/phase12_*`, figs in `figures/spatial/phase12_*` | WARNING if not |
| V29 | Package installations logged | Each install: command, packageVersion(), snapshot() | WARNING if unlogged |
| V30 | renv.lock change documented | Before/after MD5 recorded | PASS/FAIL |
| V31 | Rplots.pdf handled | Pre-existence recorded; deleted only if created during run | PASS/FAIL |
| V32 | Complex V caveat carried forward | Documented in all relevant outputs | WARNING if missing |
| V33 | log2FC direction annotated in all DE outputs | Positive = Old-high, Negative = Young-high | WARNING if missing |
| V34 | Guide doc generated | `docs/spatial_phase12_young_aged_region_comparison_guide.md` exists | WARNING if missing |
| V35 | Interaction model rank-checked or documented as SKIPPED | `optional_region_age_interaction_audit.csv` records decision | PASS/FAIL |
| V36 | padj=NA diagnostics file exists | `deseq2_padj_na_diagnostics_by_region.csv` readable, 3 regions recorded | WARNING if missing |
| V37 | p-value histograms exist for all 3 regions | `pvalue_histogram_CA1.png`, `CA3.png`, `DG.png` all exist | WARNING if missing |
| V38 | alpha = 0.05 recorded in provenance and output files | `results(alpha=0.05)`, `summary(res, alpha=0.05)` verified; not relying on defaults | WARNING if unverified |
| V39 | Target genes with padj=NA are counted and listed | `deseq2_padj_na_diagnostics_by_region.csv` has `n_target_genes_padj_na` and `target_genes_padj_na` columns | WARNING if missing |
| V40 | Guide explains padj=NA and p-value histogram QC | `§7.12` present in `docs/spatial_phase12_young_aged_region_comparison_guide.md` | WARNING if missing |

---

## 9. Stop/Warning Conditions (SC01-SC20)

| ID | Condition | Action |
|----|-----------|--------|
| SC-01 | Hippo RDS missing or unreadable | STOP — object unavailable |
| SC-02 | target_genes.csv missing | STOP — Phase 10 prerequisite |
| SC-03 | Phase 10 audit missing | STOP — rebuild plan required |
| SC-04 | CA1, CA3, DG missing from metadata | STOP that region; continue independent regions only if documented |
| SC-05 | Young or Old missing in any primary region | STOP that region |
| SC-06 | < 3 biological samples (GEMgroups) per age in a region | **Block DE for that region**; record as `blocked_insufficient_samples`; continue other regions |
| SC-07 | Design matrix full-rank failure | STOP that region's model |
| SC-08 | DESeq2 failure in one region | STOP that region; continue other independent regions only if documented |
| SC-09 | Pseudobulk aggregation produces all-zero samples | STOP — data integrity issue |
| SC-10 | DG pseudobulk rebuild fails | STOP DG DE; continue CA1 and CA3 if independent |
| SC-11 | CA3 Young < 3 GEMgroups (G1 excluded) | This is expected behavior (3 Young GEMgroups: G2,G3,G4). NOT a stop condition. Flag as `primary_with_caution_young_n3` but proceed. |
| SC-12 | Module score computation fails for DG | WARNING — exclude DG from module analyses; continue CA1/CA3 |
| SC-13 | Coupling n ≤ 3 per context | WARNING — label coupling as `exploratory_small_n`; continue |
| SC-14 | Complex V still missing (16 Atp5* genes) | WARNING — documented caveat; do NOT stop |
| SC-15 | User requests Middle in primary contrast | STOP — Phase 12 is Young vs Old only; Middle can be added in a future phase |
| SC-16 | User requests biological interpretation in Phase 12 output | STOP — Phase 13+ scope |
| SC-17 | User requests cell-type deconvolution | STOP — requires separate plan |
| SC-18 | User requests WholeBrain loading | STOP — requires separate plan |
| SC-19 | User requests enrichment analysis (GO/KEGG/Reactome) | STOP — Phase 13+ scope |
| SC-20 | User requests Seurat object saving | STOP — against AGENTS.md |

---

## 10. Execution Outline

```
Phase Spatial-12: Young vs Aged Regional Comparison across CA1/CA3/DG

E0: Setup
  - Re-read this plan file and record plan_version = "v1.3"
  - Rplots.pdf guard (start)
  - renv.lock MD5 before
  - Git status at start
  - Load required packages: Seurat, DESeq2, Matrix, Matrix.utils, ggplot2, data.table, digest
  - Check optional packages: apeglm, pheatmap, ggrepel, corrplot → install if beneficial, log if installed

E1: Input Audit (V01-V04)
  - Verify Hippo RDS exists and is loadable
  - Verify all required input files exist
  - Verify Phase 10 audit: 231 exact-match, 20 missing
  - Verify Phase 11 manifest: 31 samples, CA1=16, CA3=15
  - Load Hippo RDS: verify 9,921 spots, CA1=4671, CA3=2321, DG=2364

E2: Region × Age Sample Manifest
  - Build region_age_pseudobulk_sample_manifest.csv
  - For CA1: extract GEMgroups 1-4 (Young) and 9-16 (Old)
  - For CA3: extract GEMgroups 2-4 (Young; G1 excluded) and 9-16 (Old)
  - For DG: derive from ML+GCL+Hilus; extract GEMgroups 1-4 (Young) and 9-16 (Old)
  - Record n_spots per GEMgroup × region_group × Age
  - Flag low_coverage if n_spots < 20
  - Log V04-V09

E3: DG Pseudobulk Construction
  - Subset Hippo to DG spots (ML+GCL+Hilus)
  - Subset to Young + Old spots only
  - Get Spatial raw counts (sparse: 32,285 features × DG_YoungOld spots)
  - Aggregate by GEMgroup: t(counts) → aggregate.Matrix by GEMgroup → t()
  - Build dg_pseudobulk_counts.rds
  - Build DG sample manifest rows
  - GC after saving
  - Log V06, V09

E4: CA1/CA3 Pseudobulk (from existing or rebuild)
  - If ca1_ca3_pseudobulk_counts.rds exists and is valid: subset by Age (Young, Old) and by region_group
  - If missing: rebuild from Hippo RDS using same pattern as Phase 11
  - Create region-specific count matrices: CA1_only, CA3_only
  - Each: genes × GEMgroups (Young + Old only)
  - Log V07, V08

E5: Design Audit Per Region
  - For CA1, CA3, DG: build design matrix for ~ Age
  - Check full rank via qr()
  - Count samples per Age cell
  - Save: design_matrix_audit_region_age.csv
  - Check SC-07

E6: DESeq2 — CA1 Old vs Young
  - DESeqDataSetFromMatrix: CA1 pseudobulk counts + ~ Age
  - Filter: rowSums(counts) >= 10
  - DESeq()
  - results(contrast = c("Age", "Old", "Young"), alpha = 0.05)
  - summary(res, alpha = 0.05)
  - lfcShrink(type = "apeglm") or fallback "normal"
  - padj=NA diagnostics (Route A-i): record n_genes_total, n_pvalue_na, n_padj_na, n_padj_na_with_pvalue_present, target genes with NA
  - P-value histogram (Route A-ii): generate pvalue_histogram_CA1.png
  - Save: deseq2_all_genes_Old_vs_Young_CA1.csv
  - Add columns: model_type = "two_group_unpaired", power_label = "primary_adequately_powered", region = "CA1"
  - Check SC-08

E7: DESeq2 — CA3 Old vs Young
  - Same as E6 for CA3 pseudobulk counts
  - padj=NA diagnostics for CA3
  - pvalue_histogram_CA3.png
  - Save: deseq2_all_genes_Old_vs_Young_CA3.csv
  - power_label = "primary_with_caution_young_n3"
  - Check SC-08

E8: DESeq2 — DG Old vs Young
  - Same as E6 for DG pseudobulk counts
  - padj=NA diagnostics for DG
  - pvalue_histogram_DG.png
  - Save: deseq2_all_genes_Old_vs_Young_DG.csv
  - power_label = "primary_adequately_powered" (expected)
  - Check SC-08

E9: Target Gene Annotation by Region
  - Merge each region's DE results with Phase 10 gene audit
  - Create combined table: 251 genes × 3 regions
  - Annotate: padj_sig, direction, lfc_notable per region
  - Save: deseq2_target_genes_Old_vs_Young_by_region.csv
  - Log V13, V14

E10: Cross-Region log2FC Matrix (Route B)
  - Build wide matrix: rows = 231 genes, cols = CA1/CA3/DG Old-vs-Young log2FC
  - Classify each gene into direction consistency pattern
  - Save: target_gene_log2fc_matrix_by_region.csv
  - Log V15

E11: Category Effect by Region (Route C)
  - Per category × per region: n_available, n_sig, n_old_high, n_young_high, mean_log2FC, median_log2FC, direction_balance
  - Compile region_age_deg_summary.csv
  - Save: target_gene_category_effect_by_region.csv, region_age_deg_summary.csv
  - Save: region_age_power_flags.csv

E12: DG Module Scores (if not already computed)
  - Subset Hippo to DG spots
  - AddModuleScore with 12 gene sets from Phase 11 module_gene_set_final_audit.csv
  - Aggregate spot-level to GEMgroup-level (mean and median)
  - Merge into extended module_scores_sample_summary.csv
  - GC after computation

E13: Module Age Effect by Region (Route D)
  - For each module × region:
    - Subset GEMgroup-level scores to Young and Old
    - Wilcoxon rank-sum test (unpaired)
    - BH adjustment across modules within each region
  - Save: module_age_effect_by_region.csv
  - Log V16, V17

E14: Age-Stratified Coupling by Region (Route E, EXPLORATORY)
  - For each region, for each age (Young, Old):
    - Spearman correlation for 10 coupling pairs
    - n ≤ 4: rho only, no p-value
    - n = 8: rho + raw p-value, labeled exploratory_small_n
  - Save: module_coupling_young_old_by_region.csv
  - Log V18, V19

E15: Optional Interaction Model (Route F, Default SKIP)
  - Document decision in optional_region_age_interaction_audit.csv
  - If attempted: build ~ region_group + Age + region_group:Age, check rank
  - Log V35

E16: Plots
  - F1-F3: Volcano plots per region (CA1, CA3, DG)
  - F4: Target gene log2FC heatmap by region (231 genes × 3 regions)
  - F5: Category effect heatmap by region (13 categories × 3 regions)
  - F6: DEG count barplot by region
  - F7: Module age effect boxplot (12 panels × 3 regions)
  - F8: Selected target gene age effect dotplot (top 30 genes)
  - F9: Optional coupling heatmap (exploratory)
  - Plot conventions:
    - CA1 = #2166ac (blue), CA3 = #b2182b (red), DG = #4daf4a (green)
    - Age: Young = #c6dbef, Old = #08306b
    - Format: PNG, 300 dpi, width 8-14, height 6-14
    - PDF: cairo_pdf() with tryCatch fallback

E17: Analysis Guide Generation
  - Generate docs/spatial_phase12_young_aged_region_comparison_guide.md from template (§7)
  - Include all required sections

E18: Validation and Provenance
  - Run checks V01-V35
  - Save: phase12_validation_summary.txt
  - Save: phase12_provenance.csv
  - renv.lock MD5 after; log changes

E19: Cleanup
  - rm(hippo, counts, dds); gc()
  - Rplots.pdf guard (end)
  - Final exit code: 0 if PASS or PASS_WITH_WARNINGS; 1 if any STOP
```

---

## 11. Recommended: Primary, Secondary, Exploratory Summary

| Priority | Route | What | Power |
|----------|-------|------|-------|
| **PRIMARY** | A | Region-stratified DESeq2 Old vs Young within CA1, CA3, DG (3 independent runs, `~ Age` design) | CA1=adequate (4+8), CA3=cautious (3+8), DG=TBD |
| **PRIMARY** | B | Cross-region Old-vs-Young log2FC matrix + direction consistency classification | Descriptive, no new tests |
| **PRIMARY** | C | Target gene category summaries by region (aggregates A+B) | Aggregates DESeq2 results |
| **SECONDARY** | D | Module score Old-vs-Young comparison within each region (Wilcoxon rank-sum) | Small n (3-8 per age-region cell); descriptive trends |
| **EXPLORATORY** | E | Age-stratified module coupling within each region (Spearman rho only, no p-values for n≤4) | n=3-8 per context; descriptive patterns only |
| **EXPLORATORY** | F | Region × Age interaction model (SKIP by default) | Confounded with GEMgroup; sensitivity only |

---

## 12. Execution Command (Expected)

```bash
Rscript R/spatial/s12_young_aged_region_comparison.R
```

Execution phase MUST:
1. Re-read this plan file: `docs/spatial_phase12_young_aged_region_comparison_plan.md`
2. Verify plan version matches expectations
3. Proceed checkpoint-by-checkpoint:
   - CP1 (E0-E1): Setup + input audit
   - CP2 (E2-E5): Sample manifest + DG pseudobuild + CA1/CA3 pseudobulk + design audit
   - CP3 (E6-E8): DESeq2 per region (CA1 → CA3 → DG, sequentially)
   - CP4 (E9-E11): Target gene annotation + log2FC matrix + category summaries
   - CP5 (E12-E14): DG module scores + module age effect + coupling
   - CP6 (E15): Optional interaction model (default skip)
   - CP7 (E16-E17): Plots + analysis guide
   - CP8 (E18-E19): Validation + provenance + cleanup

---

> **Plan status**: v1.3 (revised — DESeq2 QC: padj=NA diagnostics + p-value histograms + alpha=0.05 + significance vs effect-size). No analysis outputs generated. No packages installed, no R scripts executed, no objects loaded, no renv.lock changes, no data/figures generated.
> **v1.3 additions**: Route A-i (padj=NA diagnostics), Route A-ii (p-value histogram QC), alpha=0.05 explicit, significance vs notable distinction, D16 + F10-F12 outputs, V36-V40 checks, §7.12 guide content, §7.13 QC backfill plan reference.
> **Next step**: User review → approval → execution.
> **Execution prerequisite**: Phase 11 base and 11b must be PASS_WITH_WARNINGS (confirmed). All required inputs must exist.
> **Script to create**: `R/spatial/s12_young_aged_region_comparison.R` (following E0-E19 outline above)
> **Ref-bio**: Accessible. DESeq2 vignette and Seurat spatial vignette reviewed in full; OSTA TOC inspected; Seurat DE returned 404.
> **bio-skills**: 3 DESeq2 skills confirm design; 5 spatial skills are Python/Squidpy (not directly applicable).
> **DESeq2 QC backfill**: Phase 03/11/11b tracked in `docs/spatial_deseq2_qc_backfill_plan.md`. Phase 12 includes diagnostics natively.
