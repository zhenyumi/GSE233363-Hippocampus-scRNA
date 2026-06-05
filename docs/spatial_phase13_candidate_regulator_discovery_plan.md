# Phase Spatial-13: CA1/CA3 Mitochondrial Module-Associated Candidate Regulator Discovery Plan

**Version**: 1.2
**Date**: 2026-06-05
**Changes**: See §18.6 v1.2 changelog
**Status**: PLAN-ONLY (no execution; no R scripts; no data/figures generated)
**Prerequisite**: Phase 03, 11, 11b, 12 completed and outputs verified

---

## 1. Phase Overview & Objectives

### 1.1 Goal

Using annotated CA1 and CA3 spatial transcriptomics data, with existing mitochondrial module scores as functional readouts, identify candidate genes whose expression associates with regional mitochondrial transcriptional programs. The output is a ranked, heuristic list of **candidate regulators**, not causal regulators.

### 1.2 Core Questions

1. Which genes' expression (at the GEMgroup pseudobulk level) correlates with mitochondrial module scores in CA1 and CA3?
2. How does prior DE evidence (CA3-vs-CA1, Old-vs-Young) inform candidate ranking?
3. Which candidates are annotated as transcription factors / regulators?
4. Which candidates overlap with the 251-gene target gene list (13 mitochondrial categories)?
5. What spatial expression patterns do top candidates show?

### 1.3 Non-goals

- NO causal regulator inference
- NO TF activity estimation as primary analysis
- NO network inference as primary analysis
- NO spot-level statistical inference
- NO WholeBrain analysis
- NO treating target gene membership as positive regulator evidence
- NO using module gene-set members as upstream regulator candidates for that module (circular evidence)

---

## 2. Source Review Table

All sources were actually accessed and reviewed during the plan phase. No source is cited without review.

### 2.1 Authoritative Documentation

| # | Source Title | URL / Path | Reviewed Status | Section/Topic Read | Design Decision Affected |
|---|-------------|-----------|-----------------|-------------------|------------------------|
| S1 | DESeq2 Official Vignette | https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html | reviewed | Standard workflow, design formulas, Wald test, LFC shrinkage (apeglm), independent filtering, multi-factor designs, contrasts | Pseudobulk DE design confirmed (`~ Age` within region); lfcShrink with apeglm for ranking; padj=NA diagnostics |
| S2 | Seurat Spatial Vignette (sequencing-based) | https://satijalab.org/seurat/articles/spatial_vignette | reviewed | Normalization (SCTransform), SpatialFeaturePlot, SpatialDimPlot, FindSpatiallyVariableFeatures (moransi), multiple slices, deconvolution vs integration | SpatialFeaturePlot usage confirmed; spot-level visualization only; region subsetting approach validated |
| S3 | Seurat AddModuleScore Reference | https://satijalab.org/seurat/reference/addmodulescore | reviewed | Function signature, parameters (features, pool, nbin=24, ctrl=100, name), control feature binning, Tirosh et al 2016 reference | Module score methodology confirmed; existing Phase 11/12 module scores are valid readouts |
| S4 | DoRothEA Homepage | https://saezlab.github.io/dorothea/ | reviewed | Mouse support via ortholog conversion (Holland 2019), confidence levels A-E, signed TF-target interactions, Bioconductor availability, CollecTRI as successor | TF annotation source selected; mouse regulons available on Bioconductor; signed interactions usable for direction annotation |
| S5 | DoRothEA Vignette | https://saezlab.github.io/dorothea/articles/dorothea.html | reviewed | `decoupleR::get_dorothea()` API, confidence filtering, mode of regulation (mor), regulon size distributions, CollecTRI access method | TF annotation API confirmed; primary high-confidence filter A-C, secondary extended filter A-D; mor column provides activation/repression annotation |
| S6 | OSCA Multisample (redirect) | https://bioconductor.org/books/release/OSCA/multi-sample-comparisons.html | partially_reviewed | Redirected to OSCA.multisample; pseudobulk chapter confirms DESeq2 as primary method | Confirms DESeq2 pseudobulk approach is standard; edgeR also listed as alternative |

### 2.2 ref-bio Routing Files

| # | File | Reviewed Status | Key Entries Used |
|---|------|----------------|-----------------|
| R1 | `.opencode/skills/ref-bio/SKILL.md` | reviewed | Trigger rules, routing-only policy, no summarization rule |
| R2 | `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml` | reviewed | DESeq2 (source_id: deseq2), Seurat (seurat), OSCA (osca), OSTA (osta), msigdbr (msigdbr), clusterProfiler (clusterprofiler), SpatialExperiment (spatialexperiment) |
| R3 | `.opencode/skills/ref-bio/reference-pack/indexes/topic-map.yaml` | reviewed | `differential_expression` → deseq2/limma/edger; `spatial_transcriptomics` → osta/seurat; `enrichment_analysis` → clusterprofiler/msigdbr; `pseudobulk_de` → pseudobulk-de-guidance/deseq2 |
| R4 | `.opencode/skills/ref-bio/reference-pack/indexes/workflow-stage-map.yaml` | reviewed | `differential_expression` → deseq2; `spatial_analysis` → seurat-spatial; `enrichment_and_pathway` → clusterprofiler |
| R5 | `.opencode/skills/ref-bio/reference-pack/indexes/package-map.yaml` | reviewed | Package→source mappings for Seurat, DESeq2, edgeR, limma |
| R6 | `.opencode/skills/ref-bio/reference-pack/AGENTS.reference.md` | reviewed | Authority statement; indexes are routing-only |
| R7 | `.opencode/skills/ref-bio/reference-pack/policies/source-priority.md` | reviewed | Priority: project docs → paper/author code → package docs → OSCA/OSTA → peer-reviewed → tutorials |
| R8 | `.opencode/skills/ref-bio/reference-pack/MANIFEST.yaml` | reviewed | 50 sources, all link-only, no upstream content acquired |

### 2.3 bio-* Skills

| # | Skill Path | Reviewed Status | Key Guidance Used | Design Decision Affected |
|---|-----------|----------------|-------------------|------------------------|
| B1 | `bio-differential-expression-deseq2-basics/SKILL.md` | reviewed | Negative-binomial GLM, Wald test, LFC shrinkage (apeglm/ashr), design formula with interaction trap, padj=NA diagnosis, pre-filtering, contrast syntax | Confirms Phase 11/12 DESeq2 approach; padj=NA diagnostics used in Phase 12; design `~ Age` within region is correct for two-group unpaired |
| B2 | `bio-differential-expression-de-results/SKILL.md` | reviewed | padj=NA causes (independent filtering, Cook's outliers, all-zero), BH FDR, IHW, Storey q-value, gene annotation, TREAT vs post-hoc, p-value histogram diagnostics | DE results from Phase 12 will be used as ranking weights; padj thresholds; p-value QC already done in Phase 12 |
| B3 | `bio-differential-expression-de-visualization/SKILL.md` | reviewed | MA plot (shrunken LFC), volcano (apeglm caveat: p-values unchanged), PCA on VST/rlog, sample distance heatmap, row-scaling trap, UpSet plots, ggrepel | Visualization conventions for candidate ranking plots; volcano and MA plots for DE-informed module |
| B4 | `bio-gene-regulatory-networks-coexpression-networks/SKILL.md` | reviewed | WGCNA signed networks, soft-threshold selection, bicor correlation, eigengenes, kME hubs, module preservation, scale-free topology caveat (Broido & Clauset 2019) | WGCNA noted as exploratory option but NOT primary — small n (CA1=12, CA3=11) below recommended minimum (~15 samples) |
| B5 | `bio-gene-regulatory-networks-grn-inference/SKILL.md` | reviewed | ARACNe-AP, GENIE3, GRNBoost2 for edge inference; VIPER/msVIPER for protein activity; DREAM5 wisdom-of-crowds; activity-not-edges paradigm | VIPER activity inference noted as exploratory only — pseudobulk (n≤12) is below recommended sample size for robust regulon activity estimation |
| B6 | `bio-gene-regulatory-networks-scenic-regulons/SKILL.md` | reviewed | pySCENIC three-step: GRNBoost2 → cisTarget motif pruning → AUCell; motif pruning as directionality; RSS specificity scores; database matching | pySCENIC is for single-cell data, not pseudobulk; deferred to exploratory if sufficient per-spot data is available (likely not needed for current goals) |
| B7 | `bio-data-visualization-heatmaps-clustering/SKILL.md` | reviewed | ComplexHeatmap, pheatmap, ward.D vs ward.D2, OLO, robust bounds, annotation tracks, row-scaling trap | Heatmap rendering for correlation matrix and candidate expression patterns; use ward.D2 + robust bounds |
| B8 | `bio-data-visualization-volcano-and-ma-plots/SKILL.md` | reviewed | EnhancedVolcano, ggplot2 volcano, shrunken LFC, label selection | Volcano for DE-informed candidate visualization |
| B9 | `bio-data-visualization-distribution-plots/SKILL.md` | reviewed | Box/violin/beeswarm with N-based decision tree, Weissgerber 2015 sample-size honesty | Distribution plots for candidate expression with N-visible encoding |
| B10 | `bio-data-visualization-ggplot2-fundamentals/SKILL.md` | reviewed | grammar of graphics, cairo_pdf, theme_classic, scales, faceting, programmatic aes | All output figures follow ggplot2 conventions with cairo_pdf export |
| B11 | `bio-pathway-analysis-go-enrichment/SKILL.md` | reviewed | enrichGO, bitr, background universe, simplify, GOseq for bias correction | GO enrichment on top candidates as exploratory module (not primary) |
| B12 | `bio-pathway-analysis-gsea/SKILL.md` | reviewed | GSEA vs ORA, ranking metrics (log2FC, signed p-value), gseGO, NES interpretation, MSigDB Hallmark | GSEA exploratory option if correlation ranking produces ordered gene list |
| B13 | `bio-spatial-transcriptomics-spatial-preprocessing/SKILL.md` | reviewed | Python/Squidpy — conceptual only for R pipeline; QC metrics, MT%, normalization, HVG | Conceptual spatial QC framework; NOT executable in R; confirms spatial preprocessing steps |
| B14 | `bio-spatial-transcriptomics-spatial-visualization/SKILL.md` | reviewed | Python/Squidpy — conceptual only; spatial_scatter, gene expression on tissue, multi-panel | Conceptual visualization framework; R equivalent is SpatialFeaturePlot in Seurat |

---

## 3. Route Comparison and Convergence

### 3.1 Route 1: Expression-Module Association Ranking (PRIMARY)

**Method**: For each gene in the candidate universe, compute Spearman correlation (ρ) between GEMgroup-level expression and each mitochondrial module score, stratified by region (CA1, CA3). Rank by |ρ|. Optionally, age-stratified with small-n flag.

**Strengths**:
- Simple, transparent, requires no additional packages
- Uses existing Phase 11 pseudobulk count matrices
- Directly answers "which genes covary with mitochondrial programs"
- Spearman robust to outliers and non-linearity
- Naturally supports region stratification

**Limitations**:
- Correlation ≠ causation
- Small n: CA1 n=12 (4Y+8O), CA3 n=11 (3Y+8O); age-stratified CA3 Young n=3
- **Circular evidence risk**: If gene g is a member of module m's gene set, then ρ(g, m) is self-correlation and CANNOT be interpreted as upstream regulator evidence. Such genes MUST be excluded from regulator ranking for that module (or use leave-one-out module scores with explicit `readout_member_flag`).
- Module scores and gene expression share same source data (circularity risk mitigated by module score being aggregate of many genes — but not eliminated)
- No directionality (regulator → module vs module → regulator)

**Adopted**: YES — Primary route. Provides foundational ranking that Route 2 weights refine.

### 3.2 Route 2: DE-Informed Candidate Ranking (PRIMARY)

**Method**: Augment Route 1 correlation score with DE evidence:
- Phase 11 CA3-vs-CA1 log2FC (region bias)
- Phase 12 Old-vs-Young log2FC in CA1 and CA3 (aging association)
- Combined via weighted scoring (see §8.1)

Target gene category membership is reported as annotation and in a separate target/readout table (§9.3); it does NOT add `regulator_score` weight.

**Strengths**:
- Utilizes existing validated DE outputs
- Aging direction informs biological relevance
- Region-specific DE distinguishes CA1-biased vs CA3-biased candidates

**Limitations**:
- Composite score is heuristic, not a formal test
- Weights are arbitrary and must be documented
- Correlation and DE are not independent (genes highly DE will have constrained correlation)

**Adopted**: YES — Primary route. DE weights refine Route 1 correlation ranking.

### 3.3 Route 3: TF/Regulator Annotation-Based Ranking (SECONDARY)

**Method**: Annotate candidate universe with TF/regulator status using:
- **DoRothEA mouse regulons** (Bioconductor, `dorothea` package) — Primary annotation source. Confidence levels A-D. Signed interactions. Mouse orthologs from human. Primary high-confidence filter: A-C; secondary extended filter: A-D.
- **CollecTRI** via `decoupleR::get_collectri()` — Upgraded alternative with better coverage. Can be used as cross-validation.
- **GO annotations** (`GO:0003700` DNA-binding transcription factor activity, `GO:0140110` transcription regulator activity via `org.Mm.eg.db`) — Supplementary annotation for genes not in DoRothEA/CollecTRI.

**Critical caveat**: DoRothEA/CollecTRI annotation only proves the gene has known TF/regulator annotation or regulon membership in general (across tissues/contexts). It does NOT prove the gene regulates mitochondrial modules in this hippocampal dataset. TF annotation is used as an eligibility/annotation flag, not as regulatory evidence.

**Why DoRothEA/CollecTRI**:
- Available on Bioconductor → renv-compatible
- Supports mouse (via ortholog conversion, validated in Holland et al 2019)
- Signed interactions (activation/repression) via `mor` column
- Confidence levels allow filtering (recommend A-D, or A-C for high confidence)
- No additional infrastructure needed beyond R/Bioconductor

**Limitations**:
- TF annotation ≠ regulatory activity in this tissue context
- DoRothEA is based on human orthologs → mouse gene symbols may not perfectly map (need name convention check: human UPPERCASE vs mouse first-letter uppercase)
- CollecTRI has broader coverage but is literature-based (not tissue-specific)

**Adopted**: YES — Secondary annotation layer. TF annotation used as eligibility flag and secondary priority annotation. NOT used as a major positive weight in `regulator_score`. NOT used for causal claims.

### 3.4 Route 4: Network/Enrichment/Activity Exploratory (EXPLORATORY)

**Method**: Secondary exploratory analyses:
- **VIPER / decoupleR TF activity estimation**: Use DoRothEA/CollecTRI regulons with decoupleR to estimate TF activity from pseudobulk expression. **DEFAULT: NOT RUN**. Only executed if explicitly approved after primary results are reviewed. Even if run, results MUST be labeled EXPLORATORY and MUST carry the caveat that pseudobulk n≤12 is below recommended sample size for robust activity estimation.
- **GO enrichment** (clusterProfiler): GO BP/MF enrichment on top-ranked candidates per module per region. Provides functional context for candidate lists. ONLY installed and run if exploratory enrichment is approved after primary analysis.
- **WGCNA**: Co-expression modules from spot-level expression data. **CAVEAT**: n=9,921 spots but pseudobulk n=11-12; WGCNA on pseudobulk is invalid. Spot-level WGCNA is possible but must flag "spots ≠ replicates". **NOT recommended as default.**
- **GSEA**: If correlation ranking produces ordered gene list, run GSEA against MSigDB Hallmark gene sets for functional interpretation. ONLY if exploratory enrichment is approved.

**Limitations**:
- WGCNA invalid at pseudobulk level (too few samples); spot-level WGCNA suffers from spatial autocorrelation
- VIPER/decoupleR under-powered at n≤12
- All results labeled EXPLORATORY, not used in primary ranking

**Adopted**: YES — Exploratory only. Clearly separated from primary ranking. Not used to select candidates.

### 3.5 Route 5: Spatial Exploratory (EXPLORATORY)

**Method**: SpatialFeaturePlot for top 10-20 candidates in each class, showing expression overlaid on tissue. Visual confirmation of CA1-vs-CA3 expression patterns.

**Adopted**: YES — Exploratory visualization only. No spot-level p-value inference.

### 3.6 Convergence Summary

| Route | Status | Role in Ranking |
|-------|--------|-----------------|
| Route 1: Correlation | PRIMARY | Effect size (|ρ|) — 50% weight |
| Route 2: DE-informed | PRIMARY | Aging/region log2FC magnitude — 30% weight |
| Route 3: TF annotation | SECONDARY | Eligibility/annotation flag — 10% secondary priority weight only if in primary high-confidence A-C; D-level = annotation only |
| Route 2b: Target gene membership | SECONDARY | Separate target/readout table; NO positive weight in regulator_score; membership is a RED FLAG for self-correlation |
| Route 4: Network/Enrichment | EXPLORATORY | Separate output; default NOT run |
| Route 5: Spatial maps | EXPLORATORY | Separate output |

---

## 4. Module A: Inputs and Prerequisites

### 4.1 Required Inputs

| Input | Source | Format | Description |
|-------|--------|--------|-------------|
| Pseudobulk count matrix (CA1/CA3) | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds` | RDS | GEMgroup-level pseudobulk counts for CA1 and CA3 regions. **Primary expression data source for Phase 13.** From Phase 11. Must verify exact file name during E0. |
| Module gene set audit | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv` | CSV | Authoritative list of module names and gene-set members. Read during E0 to determine actual module count and names. |
| Module scores (all regions) | `data/processed/spatial/phase12_young_aged_region_comparison/module_scores_all_regions.csv` | CSV | GEMgroup-level mean module scores per region×age from Phase 12. Reading this file during E0 determines the actual module name list. |
| Module scores (spot-level) | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_scores_spot_summary.csv` | CSV | Per-spot module scores (for spatial visualization only) |
| Phase 11 DE results (CA3-vs-CA1) | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/deseq2_all_genes_CA3_vs_CA1.csv` | CSV | Full gene DE results from Phase 11 |
| Phase 12 DE results (Old-vs-Young) | `data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA1.csv` and `..._CA3.csv` | CSV | Full gene DE by region from Phase 12. Must verify exact file names during E0. |
| Target gene list | `docs/target_genes.csv` | CSV | 251 genes in 13 categories |
| Target gene audit | `data/processed/spatial/phase10_target_gene_audit/gene_presence_audit.csv` | CSV | Which target genes are present in Hippo data |
| Sample manifest | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/pseudobulk_sample_manifest.csv` | CSV | GEMgroup → Age, Region mapping |
| Hippo Seurat object | `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` | RDS | Loaded ONLY for spatial visualization (Route 5); loaded last; no spot-level inference; unloaded immediately after figures |

### 4.2 Derived Inputs (Generated by Phase 13)

| Input | How Generated | Description |
|-------|--------------|-------------|
| Candidate universe | Three separate universes (see §6.1) | `universe_all`: expressed genes; `universe_regulator`: DoRothEA/CollecTRI/GO TF subset; `universe_target`: target_genes.csv subset |
| TF annotation table | DoRothEA/CollecTRI + GO | Which candidates are annotated TFs/regulators |
| Module name list | Read `module_gene_set_final_audit.csv` and `module_scores_all_regions.csv` column headers | Actual module names (do NOT hardcode module count or names) |
| Module score matrix (CA1/CA3 GEMgroup-level) | From `module_scores_all_regions.csv` | Module readout per sample |
| Gene expression matrix (CA1/CA3 GEMgroup-level) | From pseudobulk count RDS | log2(CPM+1) normalized per-sample expression |

---

## 5. Module B: Mitochondrial Readout Definition

### 5.1 Module Readout Definition

**CRITICAL: Do NOT hardcode module names or counts.** During E0, read:
- `module_gene_set_final_audit.csv` from Phase 11 — determines exact module names, gene-set members per module, and module count
- `module_scores_all_regions.csv` from Phase 12 — determines which module scores are available as columns

**Expected modules** (based on Phase 11/12 outputs; must verify during E0):

| Expected Module Name | Gene Count (approximate) | Category |
|---------------------|-------------------------|----------|
| Module_mtDNA_encoded | 9 | mtDNA-encoded genes |
| Module_Complex_I | ~34 | Nduf/Ndufv/Ndufb/Ndufc |
| Module_Complex_III | ~11 | Uqcr/Uqcrc/Uqcc |
| Module_Complex_IV | ~13 | Cox subunits |
| Module_Complex_V (EXCLUDED) | 16 | Atp5 — all 16 missing |
| Module_TCA_cycle | ~12 | TCA key enzymes |
| Module_Mitoribosomal | ~26 | MRPL/MRPS |
| Module_Mito_Transport | ~15 | VDAC/TOMM/TIMM/SLC25 |
| Module_Defense_Antioxidant | ~24 | PRDX/GPX4/SOD1/GLRX/CHCHD |
| Module_RPL | ~37 | Cytosolic large ribosomal subunit |
| Module_RPS | ~27 | Cytosolic small ribosomal subunit |
| Module_Eif* | Varies | Translation initiation factors |
| Module_Eef* | Varies | Translation elongation factors |

**WARNING**: If `Module_Translation_Factors` appears in the actual audit file (combining Eif+Eef), use that name. If Eif and Eef are separate modules, keep them separate. Do NOT merge modules that are distinct in the source files.

**Primary readout modules** (direct mitochondrial function — exact count depends on actual module list): mtDNA_encoded, Complex_I, Complex_III, Complex_IV, TCA_cycle, Mitoribosomal, Mito_Transport, Defense_Antioxidant.

**Secondary readout modules** (translation machinery): RPL, RPS, Eif*, Eef*.

Complex V is EXCLUDED from all analyses (16/16 Atp5* genes missing).

### 5.2 Module Score Source and Caveats

- Scores computed via `Seurat::AddModuleScore()` using `slot = "data"` (log-normalized counts)
- nbin=24, ctrl=100 (default parameters per Seurat documentation)
- Scores are spot-level; aggregated to GEMgroup-level by taking the mean per GEMgroup
- **Complex V caveat**: All 16 Atp5* genes missing from Hippo data. Module_Complex_V excluded from all analyses. This caveat must be carried forward in all Phase 13 outputs.
- No batch correction was applied to module scores (consistent with Phase 11/12)

### 5.3 Module Score Comparability and Circular Evidence Safeguards

- Module scores are NOT comparable across modules (different gene counts, different expression levels)
- Module scores ARE comparable across regions and ages within the same module (same gene set used)
- For correlation analysis: use per-module scores directly (within-module comparisons)
- For cross-module summary: rank-transform module scores before combining

**Circular evidence safeguard (CRITICAL)**:

If gene g is a member of module m's gene set (per `module_gene_set_final_audit.csv`), then the correlation ρ(g, m) is **self-correlation** and CANNOT be interpreted as upstream regulator evidence. Two strategies:

**Strategy A (default, recommended)**: Exclude all module gene-set members from the regulator ranking for that specific module. The gene can still be ranked for other modules where it is NOT a member.

**Strategy B (only if explicitly approved)**: Use leave-one-out module scores where gene g is removed from module m's gene set before recomputing the module score. Flag such genes with `readout_member_flag = TRUE`. This is only viable if the execution phase can legitimately recompute module scores on the spot-level data and re-aggregate.

**Regardless of strategy**: All outputs MUST include a `self_correlation_flag` column indicating:
- `regulator_candidate` — gene is NOT a member of the correlated module's gene set
- `module_member_or_effector` — gene IS a member; ρ is self-correlation
- `target_gene_candidate` — gene is in `target_genes.csv` but NOT in the correlated module's gene set (may be an effector in a different pathway)
- `low_confidence_self_correlation` — gene is a member; ρ cannot be interpreted as regulator evidence

**Target gene membership**: Genes in `target_genes.csv` must NOT receive a positive weight in the `regulator_score`. They can be reported in a separate target/readout table and used as an annotation column. Target membership provides context ("this gene is involved in mitochondrial biology") but does not constitute upstream regulatory evidence.

---

## 6. Module C: Candidate Regulator Universe

### 6.1 Universe Definitions (Three Separate Universes)

**`universe_all`** — Broad association screen (for correlation only, NOT for regulator ranking):
- Expressed in >= 3 GEMgroups in BOTH CA1 and CA3
- Expression threshold: mean log2(CPM+1) >= 1
- Purpose: Run correlation against all 11 modules to understand global association landscape
- NOT used for candidate regulator ranking directly — only for background distribution and exploratory patterns

**`universe_regulator`** — Candidate regulator ranking (PRIMARY ranking universe):
- Must pass the `universe_all` expression filter
- Must be annotated in at least one of:
  - DoRothEA mouse regulons (high-confidence A-C for primary ranking; D for secondary)
  - CollecTRI mouse regulons
  - GO `GO:0003700` (DNA-binding transcription factor activity) or `GO:0140110` (transcription regulator activity) via `org.Mm.eg.db`
- Purpose: This is the universe for primary candidate regulator ranking
- DoRothEA/CollecTRI annotation is an eligibility criterion for this universe, NOT a weight

**`universe_target`** — Target/readout analysis (SEPARATE, NOT regulator ranking):
- Must pass the `universe_all` expression filter
- Must be in `docs/target_genes.csv`
- Purpose: Separate target-gene readout table; assess which target genes associate with which modules
- Target genes do NOT receive positive `regulator_score` weight
- Target genes that are also module members are flagged for self-correlation

### 6.2 Subsets and Flags

| Subset/Flag | Definition | Purpose |
|-------------|-----------|---------|
| `universe_all` | All genes passing expression filter | Broad correlation screen; background distribution |
| `universe_regulator` | `universe_all` ∩ (DoRothEA A-C ∪ CollecTRI ∪ GO-TF) | PRIMARY regulator ranking universe |
| `universe_regulator_extended` | `universe_all` ∩ (DoRothEA A-D ∪ CollecTRI ∪ GO-TF) | Extended ranking with D-level TFs |
| `universe_target` | `universe_all` ∩ target_genes.csv | Separate target/readout table |
| `flag_tf_high_confidence` | DoRothEA A-C or CollecTRI | High-confidence TF annotation flag |
| `flag_tf_extended` | DoRothEA D or GO-TF only | Lower-confidence TF annotation flag |
| `flag_target_gene` | In target_genes.csv | Target gene annotation (not a positive weight) |
| `flag_module_member(m)` | In module m's gene set | Self-correlation flag for module m |
| `universe_de_ca1` | padj < 0.05 in Phase 12 CA1 Old-vs-Young | CA1 aging-associated subset |
| `universe_de_ca3` | padj < 0.05 in Phase 12 CA3 Old-vs-Young | CA3 aging-associated subset |
| `universe_de_region` | padj < 0.05 in Phase 11 CA3-vs-CA1 | Region-biased subset |

### 6.3 Expected Universe Size

- Total detected genes in CA1/CA3: ~15,000 (from Phase 11 data)
- `universe_all` (after expression filter): ~12,000 genes
- `universe_regulator` (expression filter + DoRothEA A-C/CollecTRI/GO-TF): ~1,500-2,000 TFs
- `universe_regulator_extended` (+ DoRothEA D): ~2,000-2,500 TFs
- `universe_target` (expression filter + target_genes.csv): ~231 genes (exact match after Phase 10 audit)
- These are three separate universes with different purposes; they are NOT merged into one "criterion 1 AND any of 2-5" list

### 6.4 Gene Symbol Convention

Mouse gene symbols follow the "First letter uppercase, rest lowercase" convention (e.g., "Tfrc" not "TFRC").

- DoRothEA mouse regulons should already use mouse symbols (confirm during execution)
- If human symbols are returned, convert via `stringr::str_to_title()` (consistent with Phase 10 gene symbol mapping)
- Cross-check against `phase10/gene_symbol_mapping_audit.csv` for known naming issues

---

## 7. Module D: Association Models

### 7.1 Primary Model: GEMgroup-Level Spearman Correlation

**Design**: For each gene g and each module m, compute Spearman ρ(g, m) in:
- CA1 (n=12 GEMgroups: 4 Young + 8 Old)
- CA3 (n=11 GEMgroups: 3 Young + 8 Old)

**Formula**: `cor.test(expr[g, CA1_samples], module_score[m, CA1_samples], method = "spearman")`

**Input data**:
- Expression: log2(CPM + 1) from pseudobulk counts, per gene, per GEMgroup
- Module score: GEMgroup-level mean module score from `module_scores_all_regions.csv`

**Output**: Per-gene, per-module, per-region: Spearman ρ, p-value, n

### 7.2 Age-Stratified Model (Secondary)

Same as primary but stratified:
- CA1 Young (n=4), CA1 Old (n=8)
- CA3 Young (n=3), CA3 Old (n=8)

**CAVEAT**: CA3 Young n=3 is severely underpowered. All CA3 Young results labeled `primary_with_caution_young_n3`.

### 7.3 Age-Interaction Model (Exploratory)

For genes where Young and Old show opposing correlation directions:
- Compute ρ_young - ρ_old for each region
- Flag genes with |Δρ| > 0.5 and opposite signs
- Label "potential age-dependent association"
- No statistical test; purely descriptive

### 7.4 Multiple Testing and Small-n Constraints

- Number of tests per region: ~12,000 genes × 11 modules ≈ 132,000 tests (for `universe_all`)
- **CRITICAL**: At n=11-12 per region, correlation p-values and q-values are **exploratory/supporting only**. They have severely limited power and cannot distinguish moderate from strong effects reliably.
- Apply Benjamini-Hochberg correction per region for transparency, but report as `q_exploratory` (not `padj`)
- **Primary ranking MUST use effect size |ρ| as the primary metric**, combined with DE evidence and direction consistency
- BH q-values are reported for completeness but MUST NOT be used as a primary ranking criterion or as evidence of statistical significance
- **All correlation figures and tables MUST display n** (sample count) prominently
- Do NOT threshold by q-value for candidate selection; use |ρ| > 0.3 as a screening threshold and |ρ| > 0.5 as a "notable" threshold

### 7.5 Key Guardrails

- **Replicate unit = GEMgroup, NOT spots**. No spot-level correlation.
- **No causal direction**: correlation is symmetric; we cannot distinguish "gene regulates module" from "module correlates with gene expression"
- **No formal model comparison**: CA1 vs CA3 correlations are reported descriptively; no test for "differential correlation"
- **Small-n explicit**: Every plot and table must display n.

---

## 8. Module E: Candidate Ranking

### 8.1 Composite Score Formula

For each gene g and module m in `universe_regulator`, compute:

```
regulator_score(g, m, region) = w1 * |ρ(g,m)|_scaled + w2 * |DE_age(g, region)|_scaled + w3 * |DE_region(g)|_scaled
```

Where:
- `|ρ(g,m)|_scaled` = absolute Spearman ρ scaled to [0,1] by dividing by max |ρ| in that module×region
- `|DE_age(g, region)|_scaled` = absolute log2FC from Phase 12 Old-vs-Young in region, scaled to [0,1]
- `|DE_region(g)|_scaled` = absolute log2FC from Phase 11 CA3-vs-CA1, scaled to [0,1]

**Default weights**: `w1=0.50, w2=0.30, w3=0.20`

**What is NOT in the score**:
- **TF annotation is NOT a weight** — it is an eligibility criterion for `universe_regulator` and an annotation flag
- **Target gene membership is NOT a weight** — it is a separate annotation column and may indicate self-correlation risk
- **p-values/q-values are NOT in the score** — they are exploratory only at n=11-12

### 8.1b Annotation Columns (NOT Weights)

For each ranked gene, report these annotation columns separately:
- `tf_confidence`: "high" (DoRothEA A-C / CollecTRI), "extended" (DoRothEA D / GO-TF only), "none"
- `target_gene_category`: category from `target_genes.csv` (e.g., "Complex I", "Mitoribosomal"), or "none"
- `self_correlation_flag`: per §5.3 classification
- `de_ca1_sig`: TRUE if padj < 0.05 in Phase 12 CA1
- `de_ca3_sig`: TRUE if padj < 0.05 in Phase 12 CA3
- `de_region_sig`: TRUE if padj < 0.05 in Phase 11

These annotations DO NOT contribute to `regulator_score`. They appear in output tables for interpretation.

### 8.1c Sensitivity Analysis

Run ranking with at least 3 alternative weight configurations:
- **Default**: w=(0.50, 0.30, 0.20) — ρ-dominant
- **DE-emphasis**: w=(0.33, 0.50, 0.17) — prioritize aging/region DE evidence
- **ρ-only**: w=(1.00, 0.00, 0.00) — pure correlation ranking (no DE integration)

Report top 20 candidate overlap across configurations. If overlap < 50%, flag "ranking sensitive to weight choice" and report both the default and the consensus (appears in top 20 in all 3 configurations) rankings.

### 8.2 Region-Specific Rankings

Separate rankings for:
1. **CA1 candidates**: genes ranked by `regulator_score(g, m, "CA1")` → CA1 aging-associated mitochondrial program candidates
2. **CA3 candidates**: genes ranked by `regulator_score(g, m, "CA3")` → CA3 aging-associated mitochondrial program candidates
3. **Combined candidates**: genes ranked by max(regulator_score CA1, regulator_score CA3) → strongest overall candidates

### 8.3 Per-Module Rankings

Within each module m:
- Top 10 candidates per region
- Report: gene, ρ, DE_age_log2FC, DE_region_log2FC, tf_confidence, target_gene_category, regulator_score, self_correlation_flag

### 8.4 Cross-Module Consensus

For genes appearing in top 10% (by regulator_score) for ≥3 modules:
- Flag as "cross-module candidate"
- Higher confidence for broad mitochondrial program association
- Label "robust across modules" vs "module-specific"

### 8.5 Ranking Caveats

- **HEURISTIC ONLY**: Regulator score has no statistical meaning beyond relative ranking within a module×region
- **NOT a p-value**: Regulator score must never be reported as a significance level
- **NOT a regulator probability**: High score does NOT imply the gene regulates the module
- **NOT causal**: No directionality is established
- **Weights are subjective**: Default weights are documented; sensitivity analysis (§8.1c) assesses stability
- **Self-correlation excluded**: Module gene-set members are excluded from regulator ranking for that module (Strategy A, §5.3)
- **Target genes annotated separately**: Target gene membership is a separate annotation column, not a ranking weight
- **TF annotation = eligibility**: Being in DoRothEA/CollecTRI qualifies a gene for `universe_regulator`; it does not add weight
- **Output file `composite_score_components.csv`** must include every component (ρ, DE_age, DE_region, each scaled value) so the provenance of every score is traceable

---

## 9. Module F: Candidate Classes

### 9.1 Class Definitions

| Class | Criteria | Interpretation |
|-------|----------|---------------|
| **CA1 aging-associated** | Top regulator_score in CA1; |ρ| > 0.3; DE CA1 Old-vs-Young padj < 0.05 preferred | Candidate associated with CA1 mitochondrial aging |
| **CA3 aging-associated** | Top regulator_score in CA3; |ρ| > 0.3; DE CA3 Old-vs-Young padj < 0.05 preferred | Candidate associated with CA3 mitochondrial aging |
| **Region-biased** | |ρ_CA1 - ρ_CA3| > 0.3; DE CA3-vs-CA1 padj < 0.05 | CA1-specific or CA3-specific mitochondrial program association |
| **Concordant** | Top regulator_score in BOTH CA1 and CA3; same direction of ρ | Pan-regional mitochondrial program association |
| **Region-flipped** | ρ_CA1 and ρ_CA3 of opposite sign; |ρ| > 0.4 in both | Potentially region-specific mitochondrial program regulation |
| **Target-readout** | Gene in `target_genes.csv`; NOT used in regulator_score | Separate target-gene readout table |
| **Module-member** | Gene is member of the correlated module's gene set | Self-correlation; excluded from regulator ranking for that module |
| **Low-confidence** | |ρ| < 0.3 in all regions or n < 8 | Below screening threshold; reported but not highlighted |

### 9.2 Expected Class Sizes (Estimated)

| Class | Expected Count | Rationale |
|-------|---------------|----------|
| CA1 aging-associated | ~50-100 | CA1 has ~316 DEGs (Phase 12) |
| CA3 aging-associated | ~20-50 | CA3 has ~60 DEGs (Phase 12) |
| Region-biased | ~50-100 | CA3-vs-CA1 has many DEGs |
| Concordant | ~20-50 | 60+3 genes concordant Old-high/Young-high (Phase 12) |
| Region-flipped | ~10-20 | 168 direction-flipped genes (Phase 12); only some will associate with modules |

### 9.3 Target Gene Candidates (Separate, NOT Regulator Ranking)

Target genes (231 expressed from `target_genes.csv`) are handled separately:

- **They are NOT in `universe_regulator`** for the purpose of composite scoring (unless they also have independent TF annotation)
- **Report in a separate table**: `target_gene_module_association.csv`
  - Columns: gene, category, ρ per module per region, self_correlation_flag, DE evidence
- **Self-correlation check**: For each target gene, check whether it is a member of the module's gene set. If yes, flag as `module_member_or_effector` and note that ρ is not independent evidence.
- **Target genes that also have TF annotation** (e.g., a TF in target_genes.csv that regulates mitochondrial genes): These can appear in `universe_regulator` for OTHER modules where they are NOT gene-set members, but the TF annotation must be independently verified (not inferred from target gene status).

Key Phase 12 reference:
- Phase 12 CA1 DEG target genes: 5 (Old-vs-Young)
- Phase 12 CA3 DEG target genes: 0 (Old-vs-Young)
- Phase 12 DG DEG target genes: 0 (Old-vs-Young)

### 9.4 TF Candidate Enrichment

- Test enrichment of DoRothEA TFs among top-ranked candidates vs background using Fisher's exact test
- Report: OR, 95% CI, p-value
- Interpret as "TFs are overrepresented among top candidates" (descriptive, not causal)

---

## 10. Module G: Visualizations

All figures saved to `figures/spatial/phase13_candidate_regulator/`. PNG is authoritative; PDF is optional/best-effort (with cairo_pdf tryCatch fallback and explicit file paths; see §10.6).

### 10.1 Primary Screening/Ranking Figures

These figures support candidate screening and ranking. They display association signals (correlation, DE direction, score components) but are NOT strong inferential proof — all rankings are heuristic, and small n (11-12 per region) limits statistical confidence. Every correlation figure MUST display n.

| Figure | Type | Content | Priority |
|--------|------|---------|----------|
| `regulator_module_correlation_heatmap_CA1.png` | Heatmap | Top 50 candidates × 11 modules Spearman ρ, CA1 | HIGH |
| `regulator_module_correlation_heatmap_CA3.png` | Heatmap | Top 50 candidates × 11 modules Spearman ρ, CA3 | HIGH |
| `top_candidate_ranking_dotplot.png` | Dotplot | Top 20 candidates per class, showing ρ, DE, TF status, target gene status | HIGH |
| `candidate_ca1_vs_ca3_scatter.png` | Scatter | ρ_CA1 vs ρ_CA3 for top candidates, colored by class | HIGH |
| `candidate_expression_vs_module_scatter.png` | Scatter | Selected candidate expression vs module score, colored by age | HIGH |
| `regulator_score_breakdown.png` | Bar | Top 20 candidates score component decomposition | MEDIUM |

### 10.2 Age-Stratified Figures

| Figure | Type | Content | Priority |
|--------|------|---------|----------|
| `age_stratified_correlation_dotplot_CA1.png` | Dotplot | ρ_Young vs ρ_Old for top candidates, CA1 | MEDIUM |
| `age_stratified_correlation_dotplot_CA3.png` | Dotplot | ρ_Young vs ρ_Old for top candidates, CA3 (with n=3 CAVEAT) | MEDIUM |

### 10.3 Supporting Figures

| Figure | Type | Content | Priority |
|--------|------|---------|----------|
| `candidate_expression_heatmap_by_class.png` | Heatmap | Expression (log2CPM) of top candidates across GEMgroups, annotated by class | MEDIUM |
| `module_score_by_region_age_boxplot.png` | Boxplot | Module scores faceted by region×age (actual module count read from audit file) | LOW |
| `candidate_universe_upset.png` | UpSet | Overlap between universe subsets (TF, target, DE-CA1, DE-CA3) | LOW |
| `rho_distribution_histogram.png` | Histogram | Distribution of all |ρ| values by region | LOW |
| `sensitivity_weights_heatmap.png` | Heatmap | Top 20 candidates under different weight configurations | LOW |

### 10.4 Spatial Exploratory Figures (Route 5)

| Figure | Type | Content | Priority |
|--------|------|---------|----------|
| `spatial_top_candidates_CA1.png` | SpatialFeaturePlot | Top 5 CA1 candidates | EXPLORATORY |
| `spatial_top_candidates_CA3.png` | SpatialFeaturePlot | Top 5 CA3 candidates | EXPLORATORY |
| `spatial_concordant_candidates.png` | SpatialFeaturePlot | Top 3 concordant candidates | EXPLORATORY |
| `spatial_target_gene_candidates.png` | SpatialFeaturePlot | Target gene candidates with high ρ | EXPLORATORY |

### 10.5 Exploratory Figures (Route 4)

| Figure | Type | Content | Priority |
|--------|------|---------|----------|
| `go_enrichment_top_candidates.png` | Dotplot | GO BP enrichment for top 50 candidates per region | EXPLORATORY |
| `tf_enrichment_fisher_barplot.png` | Bar | Fisher test TF enrichment across candidate classes | EXPLORATORY |

### 10.6 Figure Conventions and PDF Safety

- **PNG is authoritative**. PDF is optional/best-effort.
- All correlation plots: display n (sample count) and Spearman ρ
- All boxplots/violins: show individual points (Weissgerber 2015)
- Color scheme: CA1=blue, CA3=orange (consistent with Phase 11/12)
- Age colors: Young=green, Old=purple (consistent with Phase 12)
- Module colors: distinct palette per module (read actual module count from audit file)
- Call `gc()` after each large plot

**PDF safety protocol (MANDATORY)**:
```
# Use explicit output path for EVERY PDF device
cairo_pdf(file.path(fig_dir, "figure_name.pdf"), width = 10, height = 8)
# ... plotting code ...
invisible(dev.off())
# Safety cleanup after dev.off():
while (dev.cur() > 1) invisible(dev.off())
```

- Use `cairo_pdf()` with explicit `file` path; NEVER allow default device to write to Rplots.pdf
- Use `tryCatch` for cairo_pdf fallback:
```
tryCatch({
  cairo_pdf(file.path(fig_dir, "figure_name.pdf"), width = 10, height = 8)
  # ... plotting code ...
  invisible(dev.off())
}, error = function(e) {
  message("cairo_pdf failed for figure_name: ", e$message)
  # Fallback: also use explicit path
  pdf(file.path(fig_dir, "figure_name.pdf"), width = 10, height = 8)
  # ... plotting code ...
  invisible(dev.off())
})
while (dev.cur() > 1) invisible(dev.off())
```
- After each figure with tryCatch, call `while (dev.cur() > 1) invisible(dev.off())` to close any orphaned devices
- At script end: check for and delete any generated `Rplots.pdf` (record size and mtime before deletion):
```
rplots_path <- "Rplots.pdf"
if (file.exists(rplots_path)) {
  rplots_info <- file.info(rplots_path)
  message("Removing orphaned Rplots.pdf (size: ", rplots_info$size, " bytes, mtime: ", rplots_info$mtime, ")")
  unlink(rplots_path)
}
```

---

## 11. Module H: Interpretation Guide

### 11.1 Guide Document

After execution, generate: `docs/spatial_phase13_candidate_regulator_analysis_guide.md`

The guide SHALL contain:

1. **What each analysis asks**:
   - "Correlation analysis measures whether GEMgroup-level expression of a gene covaries with mitochondrial module scores in CA1 and CA3"
   - "DE-informed ranking integrates prior knowledge about which genes change with aging or differ between regions"
   - "TF annotation tells us whether a candidate has known regulatory function, but does NOT prove it regulates mitochondrial genes here"

 2. **How to read correlation / DE / regulator score**:
    - ρ ∈ [-1, 1]; |ρ| > 0.5 considered "moderate-to-strong" at n≥8; p-values are exploratory only at n≤12
    - DE log2FC: positive = Old higher in Phase 12, CA3 higher in Phase 11
    - Regulator score: relative ranking within a module×region; higher = more lines of evidence; NOT a p-value or regulator probability

3. **Why candidate regulator ≠ causal regulator**:
   - Correlation is symmetric and could reflect:
     - The candidate regulating mitochondrial genes
     - Mitochondrial status affecting candidate expression
     - A third variable (e.g., cell composition changes) driving both
     - Technical artifact (co-expression in same cell types)
   - No perturbation experiment or mechanistic validation is performed
   - The word "candidate" is intentional: these are hypotheses, not conclusions

4. **Which outputs are primary vs exploratory**:
   - PRIMARY: Correlation ranking, DE-informed regulator_score, candidate class assignment
   - SECONDARY: TF annotation filtering, target gene overlap (separate table)
   - EXPLORATORY: GO enrichment, TF activity estimation, WGCNA (if run), spatial maps

5. **Limitations**:
   - All correlations at GEMgroup (pseudobulk) level: n=11-12 per region
   - CA3 Young n=3: correlations in this stratum are unreliable (reported with caveat)
   - Visium spots capture multiple cells: spatial heterogeneity within spots is unmeasured
   - 16 Complex V (Atp5*) genes missing: mitochondrial module coverage is incomplete
   - Module scores and gene expression share the same underlying count data
   - RDS-based workflow (not raw Space Ranger/Visium files): limited to available slots
   - Age and GEMgroup are colinear: cannot separate age from individual effects
   - Mouse genome: gene name mapping from human reference databases may have gaps

---

## 12. Package Policy

### 12.1 Existing Packages (Reused from Phase 11/12)

| Package | Version (from renv.lock) | Usage |
|---------|-------------------------|-------|
| Seurat | 5.5.0 | `SpatialFeaturePlot()`, object loading |
| DESeq2 | (from renv.lock) | Existing DE results only (read CSVs); no new DESeq2 run |
| dplyr | (from renv.lock) | Data manipulation |
| ggplot2 | (from renv.lock) | All figures |
| pheatmap | (from renv.lock) | Heatmaps (primary) |
| ComplexHeatmap | (optional — verify in renv.lock during E0; use if installed/approved) | Heatmaps (optional; fallback to pheatmap if not available) |
| patchwork | (from renv.lock) | Figure composition |
| stringr | (from renv.lock) | String manipulation |

### 12.2 New Packages (Need Installation in Execution Phase)

| Package | Source | Why Needed | Mouse Support | R Workflow | Priority |
|---------|--------|-----------|---------------|------------|----------|
| dorothea | Bioconductor (data-experiment) | Mouse TF regulon database for `universe_regulator` annotation | YES | YES (Bioconductor) | SECONDARY |
| decoupleR | Bioconductor | API to access DoRothEA/CollecTRI regulons | YES (via dorothea) | YES (Bioconductor) | SECONDARY |
| org.Mm.eg.db | Bioconductor | GO annotation for TF annotation fallback | YES (mouse) | YES (Bioconductor) | SECONDARY |
| clusterProfiler | Bioconductor | GO enrichment (exploratory only) | YES (via org.Mm.eg.db) | YES (Bioconductor) | EXPLORATORY — NOT installed by default |
| msigdbr | CRAN | MSigDB gene sets for GSEA (exploratory only) | YES (mouse gene sets) | YES (CRAN) | EXPLORATORY — NOT installed by default |

**clusterProfiler and msigdbr**: NOT installed by default during E1. Only installed if exploratory enrichment (Route 4) is explicitly approved after primary results (E4-E6) are reviewed. Default execution installs only dorothea, decoupleR, and org.Mm.eg.db.

### 12.3 Installation Protocol (Execution Phase)

**Default installation (E1)** — only the packages needed for primary/secondary analysis:

```
# 1. Record renv.lock MD5 BEFORE
tools::md5sum("renv.lock")

# 2. Install using renv (NOT BiocManager::install / install.packages)
renv::install("bioc::dorothea")
renv::install("bioc::decoupleR")
renv::install("bioc::org.Mm.eg.db")

# 3. Record versions
packageVersion("dorothea")
packageVersion("decoupleR")
packageVersion("org.Mm.eg.db")

# 4. Snapshot
renv::snapshot()

# 5. Record renv.lock MD5 AFTER
tools::md5sum("renv.lock")
```

**Conditional installation (only if exploratory enrichment approved)**:

```
renv::install("bioc::clusterProfiler")
renv::install("msigdbr")
packageVersion("clusterProfiler")
packageVersion("msigdbr")
renv::snapshot()
tools::md5sum("renv.lock")
```

Current renv.lock MD5 baseline: `b2bf05038b440576855d043e137c5883`

### 12.4 Packages NOT Needed

| Package | Why NOT Needed |
|---------|---------------|
| WGCNA | Small n at pseudobulk level; spot-level invalid due to spatial autocorrelation |
| GENIE3 / ARACNe | Network inference requires >50 samples; only 11-12 pseudobulk samples |
| pySCENIC | Python dependency; single-cell focused; not applicable to pseudobulk |
| STutility | Raw Visium files unavailable; RDS-based workflow |
| edgeR / limma | Phase 13 does no new DE; uses existing DESeq2 results |

---

## 13. Stop Conditions

### 13.1 Fatal Stop Conditions (Execution MUST STOP)

| # | Condition | Action |
|---|----------|--------|
| FATAL-01 | Phase 11 `deseq2_all_genes_CA3_vs_CA1.csv` missing or unreadable | Abort; rerun Phase 11 |
| FATAL-02 | Phase 12 `deseq2_all_genes_Old_vs_Young_CA1.csv` or `_CA3.csv` missing or unreadable | Abort; rerun Phase 12 |
| FATAL-03 | `module_scores_all_regions.csv` missing or column mismatch with `module_gene_set_final_audit.csv` | Abort; verify Phase 12 module score provenance |
| FATAL-04 | CA1/CA3 pseudobulk count RDS missing or corrupted | Abort; regenerate from Phase 11 |
| FATAL-05 | Fewer than 3 GEMgroups per region for CA1 OR CA3 | Abort; cannot compute meaningful correlation |
| FATAL-06 | Any proposed method would treat spots as replicates | Redesign; this is a hard constraint |
| FATAL-07 | Module score provenance from `module_gene_set_final_audit.csv` inconsistent with `module_scores_all_regions.csv` column names | Reconcile or regenerate module scores |
| FATAL-08 | Any source in §2 marked `not_reviewed` | Review source before proceeding |
| FATAL-09 | Any primary output labeled as "causal" or "regulates" | Reword to "candidate"/"associates"/"correlates" |
| FATAL-10 | Candidate universe has < 10 genes in `universe_regulator` | Relax expression filter; record decision |
| FATAL-11 | Module gene-set members NOT excluded from regulator ranking for own module, and no leave-one-out strategy documented | Implement Strategy A or B per §5.3 |
| FATAL-12 | Target gene membership used as a positive weight in regulator_score | Remove; target genes are annotation only (§5.3, §8.1) |

### 13.2 Non-Fatal Warnings (Proceed with Caveat)

| # | Condition | Action |
|---|----------|--------|
| WARN-01 | DoRothEA/CollecTRI mouse regulons unavailable from Bioconductor | Proceed without TF annotation; `universe_regulator` uses GO-TF only; flag in all outputs |
| WARN-02 | CA3 Young n < 3 for age-stratified correlation | Proceed with `primary_with_caution_young_n3` label |
| WARN-03 | Fewer than 50 genes in `universe_regulator` after filtering | Proceed; note limited scope in output |
| WARN-04 | Sensitivity analysis (§8.1c) shows < 50% top-20 overlap between weight configurations | Flag "ranking sensitive to weight choice"; report both default and consensus rankings |
| WARN-05 | Fewer than 10 genes with |ρ| > 0.5 in any module×region | Proceed; note limited strong associations |
| WARN-06 | PDF device fails (cairo_pdf error) | PNG is authoritative; PDF optional |
| WARN-07 | `org.Mm.eg.db` unavailable | Proceed with DoRothEA/CollecTRI only for TF annotation; no GO fallback |

---

## 14. Validation Checks (Pre-Execution)

Before running any R code, verify:

| # | Check | Expected Result |
|---|-------|----------------|
| V01 | This plan references all sources in §2 with actually-reviewed status | All `reviewed` or `partially_reviewed` (S6); zero `not_reviewed` |
| V02 | Candidate universe definition references specific data files | Input paths match Phase 11/12 output manifest |
| V03 | Module readout definition matches Phase 11 `module_gene_set_final_audit.csv` and Phase 12 `module_scores_all_regions.csv` | Actual module names and counts read from audit files at E0; Complex V excluded; mismatch between audit and scores file is FATAL (FATAL-07) |
| V04 | CA1 GEMgroup count documented | n=12 (4 Young, 8 Old) |
| V05 | CA3 GEMgroup count documented | n=11 (3 Young, 8 Old) |
| V06 | DG GEMgroup count documented | n=12 (4 Young, 8 Old) — NOT used in primary but available |
| V07 | Primary (Route 1+2), secondary (Route 3), exploratory (Route 4+5) clearly separated | Sections marked accordingly |
| V08 | No causal language in plan | "candidate regulator", "association", "correlation", "heuristic" |
| V09 | No spot-level inference planned | "GEMgroup is the pseudobulk replicate unit" stated everywhere |
| V10 | Output manifest complete | All figures and tables listed in Module G and Module E |
| V11 | Guide document planned | Module H specifies guide format |
| V12 | Package policy documented | Section 12 lists existing and new packages with rationale |
| V13 | renv protocol documented | Section 12.3 |
| V14 | Stop conditions documented | Section 13 |
| V15 | Complex V caveat carried forward | Mentioned in Module B and Module H |
| V16 | Small-n caveats explicit | CA3 Young n=3, CA1/CA3 all-age n≤12 |

---

## 15. Execution Outline (E0-E14)

**NOT executed during plan phase. Outline for future execution only.**

| Step | Description | Input | Output |
|------|-------------|-------|--------|
| E0 | Setup: load libraries, verify inputs, record renv.lock MD5 | Phase 11/12 outputs | Validation log |
| E1 | Install new packages: dorothea, decoupleR, org.Mm.eg.db (renv::install protocol, §12.3) | Bioconductor | Package version audit |
| E1b | (CONDITIONAL) Install clusterProfiler, msigdbr ONLY if exploratory enrichment approved | Bioconductor/CRAN | Package version audit |
| E2 | Load Phase 11 pseudobulk counts, read module_gene_set_final_audit.csv, Phase 12 module scores, Phase 11/12 DE results | RDS/CSV | In-memory data |
| E3 | Define three universes (§6.1); build TF annotation table from DoRothEA/CollecTRI/GO | E2 data + Bioconductor packages | `candidate_universe.csv`, `tf_annotation_audit.csv` |
| E4 | Read actual module name list from audit files; compute GEMgroup-level Spearman ρ (gene × module × region) for all three universes; apply self-correlation exclusion (§5.3) | E2 data | `correlation_all_genes.csv`, `correlation_age_stratified.csv` |
| E5 | Compute regulator_score with 3 sensitivity configurations (§8.1); rank genes in `universe_regulator` | E4 + DE results | `candidate_rankings.csv`, `composite_score_components.csv`, `sensitivity_overlap.csv` |
| E6 | Assign candidate classes (§9.1); generate `target_gene_module_association.csv` | E5 rankings + target gene audit | `candidate_classes.csv`, `target_gene_module_association.csv` |
| E7 | Primary figures (correlation heatmap, ranking dotplot, scatter, score breakdown, sensitivity) | E4-E6 | 7 PNG figures (PDF optional) |
| E8 | Age-stratified figures (with small-n labels) | E4 (age-stratified) | 2 PNG figures |
| E9 | Supporting figures (expression heatmap, boxplot, UpSet, rho histogram) | E2-E6 | 4 PNG figures |
| E10 | Load Hippo RDS ONLY for spatial visualization | Hippo RDS | Seurat object in memory (unload after E11) |
| E11 | Spatial exploratory figures (Route 5) | E6 + Hippo RDS | 4 PNG figures; rm(Hippo_obj); gc() |
| E12 | (CONDITIONAL) Exploratory enrichment (Route 4): GO, GSEA ONLY if approved after E6 review | E6 rankings | 2 figures + enrichment CSV |
| E13 | Generate interpretation guide | All outputs | `spatial_phase13_candidate_regulator_analysis_guide.md` |
| E14 | Final validation, provenance, renv snapshot, Rplots.pdf cleanup (§10.6) | All outputs | `phase13_provenance.csv`, `phase13_validation_summary.txt` |

**Memory management**:
- Load Hippo RDS ONLY in E10-E11 (last steps requiring large object; spatial visualization only)
- Unload immediately after E11: `rm(hippo_obj); gc()`
- No inference from spots; spatial maps are exploratory visualization only
- Call `gc()` between major steps

---

## 16. Output Manifest

### 16.1 Data Outputs (`data/processed/spatial/phase13_candidate_regulator/`)

| File | Description |
|------|-------------|
| `candidate_universe.csv` | All genes with universe membership flags (`in_universe_all`, `in_universe_regulator`, `in_universe_target`, `tf_confidence`, `target_gene_category`) |
| `correlation_all_genes.csv` | Per-gene, per-module, per-region Spearman ρ, p-value, q_exploratory, n, self_correlation_flag |
| `correlation_age_stratified.csv` | Age-stratified ρ with n labels |
| `candidate_rankings.csv` | Ranked candidates in `universe_regulator` with regulator_score (3 weight configurations) |
| `candidate_classes.csv` | Candidate class assignments (§9.1) |
| `composite_score_components.csv` | Per-gene component scores (raw ρ, scaled ρ, raw DE log2FC, scaled DE log2FC) — every component traceable |
| `sensitivity_overlap.csv` | Top-20 overlap across 3 weight configurations |
| `target_gene_module_association.csv` | Target genes (separate from regulator ranking) with ρ per module, self_correlation_flag, DE evidence |
| `tf_annotation_audit.csv` | DoRothEA/CollecTRI/GO TF mapping audit with confidence levels |
| `phase13_provenance.csv` | Provenance record |
| `phase13_validation_summary.txt` | Validation check results |
| `phase13_final_report.txt` | Summary statistics: N genes in each universe, top candidates per class, caveat flags |

### 16.2 Figure Outputs (`figures/spatial/phase13_candidate_regulator/`)

~18-22 PNG figures (see Module G for complete list). PDF optional, best-effort with cairo_pdf+tryCatch+dev.off safety protocol (§10.6).

### 16.3 Documentation Outputs

| File | Description |
|------|-------------|
| `docs/spatial_phase13_candidate_regulator_analysis_guide.md` | User-facing interpretation guide |

---

## 17. Relationship to Future Work

### 17.1 What Phase 13 Enables

- A ranked list of candidate genes for follow-up literature review
- TF candidates for targeted perturbation hypothesis generation
- Region-specific mitochondrial program signatures
- Input to a future causal mediation or perturbation study (outside current scope)

### 17.2 What Phase 13 Does NOT Enable

- Causal claims about any candidate
- "This TF regulates Complex I in CA1 during aging"
- Clinical or therapeutic recommendations
- Cross-species validation
- Integration with scRNA-seq data (separate modality)

### 17.3 Potential Future Phases (Out of Scope)

- Phase 14: Literature validation of top candidates
- Phase 15: TF binding site enrichment (motif analysis) in mitochondrial gene promoters
- Phase 16: Multi-omics integration (scRNA-seq + Visium) for cell-type deconvolution-based refinement
- Phase 17: Targeted experimental validation design

---

## 18. Plan Compliance Report

### 18.1 File Creation

| Check | Status |
|-------|--------|
| Only created/modified `docs/spatial_phase13_candidate_regulator_discovery_plan.md` | YES (v1.0 → v1.1 → v1.2) |
| No R scripts created | YES |
| No R analysis executed | YES |
| No RDS loaded | YES |
| No packages installed | YES |
| No renv.lock modified | YES (MD5: `b2bf05038b440576855d043e137c5883`) |
| No data/figures generated | YES |

### 18.2 Source Review Compliance (v1.2)

| Category | Count Reviewed |
|----------|---------------|
| Authoritative documentation (URLs) | 6 (S1-S6) |
| ref-bio routing files | 8 (R1-R8) |
| bio-* skills | 14 (B1-B14) |
| **Total sources reviewed** | **28** |
| Sources claimed reviewed but actually not accessed | 0 |
| Sources inaccessible and recorded as `failed_to_access` | 0 |
| Sources marked `partially_reviewed` | 1 (S6 — OSCA multisample redirect; S6 is supporting-only and does not anchor any primary design decision) |

### 18.3 Methodology Compliance

| Requirement | Status |
|-------------|--------|
| 5 routes compared | YES (Section 3) |
| Routes converged into single plan | YES (Section 3.6) |
| Primary methods identified | YES (Route 1+2) |
| Secondary methods identified | YES (Route 3) |
| Exploratory methods identified | YES (Route 4+5) |
| Rejected methods with reasons | YES (Section 12.4) |
| Every method decision cited to source | YES (Source Review Table) |
| Stop conditions defined | YES (Section 13) |
| Validation checks defined | YES (Section 14) |
| Package policy documented | YES (Section 12) |
| Interpretation guide planned | YES (Section 11/Module H) |
| No causal language in primary outputs | YES |
| No spot-level inference in primary outputs | YES |

### 18.4 Key Decisions Summary

| Decision | Rationale | Source |
|----------|----------|--------|
| Spearman correlation as primary method | Robust to outliers, non-parametric, simple, transparent | B1 (DESeq2), standard statistical practice |
| DoRothEA for TF annotation | Bioconductor, mouse support, signed interactions, peer-reviewed | S4, S5 |
| Composite heuristic scoring | Integrates multiple evidence lines; weights documented and sensitivity-tested | B1, B2, B3 |
| Pseudobulk (GEMgroup) as replicate unit | Spots are NOT biological replicates; consistent with Phase 11/12 | S1, B1, Phase 11/12 precedent |
| lfcShrink apeglm for DE ranking | DESeq2 recommended shrinkage for ranking; p-values unchanged | S1, B1, B3 |
| No WGCNA / network inference | Insufficient sample size at pseudobulk level | B4, B5 |
| No causal claims | Correlation-based; no perturbation data | S4, B5 (activity ≠ causality) |
| Exploratory GO/GSEA via clusterProfiler | Standard Bioconductor workflow for enrichment; conditional only | B11, B12 |

### 18.5 v1.2 Changelog (2026-06-05)

1. **Phase 12 pseudobulk input** (§4.1): Removed assertion that Phase 12 provides CA1/CA3 pseudobulk RDS files. Phase 13 expression data comes exclusively from Phase 11 `ca1_ca3_pseudobulk_counts.rds`. Phase 12 contributes only DE CSVs, module_scores_all_regions.csv, guide, and provenance.

2. **Target gene in Route 2** (§3.2): Removed "Target gene category membership (mitochondrial relevance)" from ranking evidence list. Added explicit clarification: target gene category is reported as annotation and in a separate table; it does not add `regulator_score` weight.

3. **Naming consistency** (§8.2, §8.3): Changed all remaining `candidate_score` references to `regulator_score`. No more mix of terms.

4. **Figure evidence language** (§10.1): Renamed "Primary Figures (Inferential)" to "Primary Screening/Ranking Figures". Added explicit caveat: figures support candidate screening, not strong inferential proof; rankings are heuristic; small n limits statistical confidence. Every correlation figure must display n.

5. **PNG/PDF rule** (§10 header, §16.2): Changed from "Both PNG and cairo_pdf" (which could imply PDF is mandatory) to "PNG is authoritative; PDF is optional/best-effort". PDF safety protocol retained.

6. **Validation V03** (§14): Removed hardcoded "12 modules, 13 categories mapped". Now reads: actual module names/counts from `module_gene_set_final_audit.csv` and `module_scores_all_regions.csv`; Complex V excluded; mismatch is FATAL (FATAL-07).

7. **Package table** (§12.1): ComplexHeatmap changed from "assumed in renv.lock" to "optional — verify during E0; use if installed/approved; fallback to pheatmap". clusterProfiler/msigdbr remain conditional only.

8. **Compliance report** (§18): Updated to v1.2. Marked as v1.0 → v1.1 → v1.2. No R execution, no RDS loaded, no packages installed, no renv.lock modified, no data/figures generated.

---

## References

### Papers Referenced in Source Review (S1-S6) and Skills (B1-B14)

These are the primary sources cited in the Source Review Table (§2) and bio-* skills. They were consulted during plan development.

1. Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biology*, 15:550.
2. Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J (2019). "Benchmark and integration of resources for the estimation of human transcription factor activities." *Genome Research*.
3. Holland CH et al. (2019). "Robustness and applicability of transcription factor and pathway analysis tools on single-cell RNA-seq data." *BBA Gene Regulatory Mechanisms*.
4. Badia-i-Mompel P et al. (2022). "decoupleR: Ensemble of computational methods to infer biological activities from omics data." *Bioinformatics Advances*.
5. Tirosh I et al. (2016). "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq." *Science*.
6. Hafemeister C, Satija R (2019). "Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression." *Genome Biology*.

### Background References (Not Directly Reviewed During Plan Phase)

These papers are cited in documentation, bio-* skills, or are standard references for methods used. They were not individually reviewed during the plan phase.

7. Stuart T et al. (2019). "Comprehensive Integration of Single-Cell Data." *Cell*.
8. Ignatiadis N et al. (2016). "Data-driven hypothesis weighting increases detection power in genome-scale multiple testing." *Nature Methods*.
9. Zhu A, Ibrahim JG, Love MI (2018). "Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences." *Bioinformatics*.
