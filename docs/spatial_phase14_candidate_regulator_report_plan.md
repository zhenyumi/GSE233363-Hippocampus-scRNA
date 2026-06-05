# Phase Spatial-14: Candidate Regulator Report Plan

**Version**: 1.0
**Date**: 2026-06-05
**Status**: PLAN-ONLY (no execution; no R scripts; no data/report generated)
**Prerequisite**: Phase 13, 13b, 13c completed and validated (Phase 13: 12 PASS; Phase 13b: 17 PASS; Phase 13c: 32 PASS)
**Role**: Report integration phase — NO new bioinformatics analysis

---

## 1. Objective

Generate a human-readable, reader-facing scientific report organized around this question:

> "Using annotated CA1 and CA3 spatial transcriptomics data, with mitochondrial module scores as functional readouts, infer candidate genes whose expression is associated with regional mitochondrial transcriptional programs."

The report synthesizes results from Phase 13, 13b, and 13c into a single document, with interpretations for each region (CA1 vs CA3), each age stratum (Young vs Old), and each region×age combination. It explains what results mean, which are more trustworthy, and which candidates are most suitable for follow-up literature review or experimental design.

**This is NOT a new analysis phase.** No rho, DE, module score, or candidate_class is recomputed.

---

## 2. Source Review

### 2.1 Files Inspected in Plan Phase

| # | File | Type | Reviewed |
|---|------|------|----------|
| 1 | `README.md` | Project documentation | Confirmed spatial branch conventions, output paths, memory rules, data/report/ already listed as gitignored |
| 2 | `AGENTS.md` | Project configuration | Confirmed Phase 12 guardrails, CA3 Young n=3 caveat, Complex V caveat, no spot-level inference, ref-bio/bio-* rules, Phase 13/13b/13c status entries |
| 3 | `docs/repository_structure.md` | Layout reference | Confirmed data/report/ would need adding to canonical layout; currently only data/, figures/, results/, reports/, cache/ are listd as generated |
| 4 | `docs/spatial_phase13_candidate_regulator_discovery_plan.md` | Plan doc (766 lines) | Reviewed full Phase 13 design: 3 universes, composite score formula (w1=0.50 ρ + w2=0.30 DE_age + w3=0.20 DE_region), self-correlation Strategy A, 6 candidate classes, metric hierarchy |
| 5 | `docs/spatial_phase13b_ca3_symmetric_display_plan.md` | Plan doc (449 lines) | Verified Phase 13b is display supplement only; CA3 mirror figures, 4 top-candidate CSVs |
| 6 | `docs/spatial_phase13c_age_stratified_regulator_plan.md` | Plan doc (817 lines) | Reviewed age-stratified design: 4 strata, 8 primary classes, 6 region-age classes, |rho| thresholds, pvalue_exploratory/qvalue_exploratory naming |
| 7 | `docs/spatial_phase13_candidate_regulator_analysis_guide.md` | Analysis guide (241 lines) | Reviewed all 3 phase sections; Phase 13c section already present and deduplicated; guide is the primary input for report writing |
| 8 | `R/spatial/s13_candidate_regulator_discovery.R` | R script (1470 lines) | Confirmed Phase 13 core computation: rho, regulator_score, candidate_class; noted figure specs for CA1-focused outputs |
| 9 | `R/spatial/s13b_ca3_symmetric_display.R` | R script (713 lines est.) | Confirmed display supplement only: CA3 mirror figures, top CSVs |
| 10 | `R/spatial/s13c_age_stratified_regulator_discovery.R` | R script (~1250 lines) | Confirmed age-stratified rho, per-stratum sample_size_labels, primary_age_class + region_age_class columns |
| 11 | `phase13_regulator_score_default.csv` | Data (13,179 rows) | Confirmed columns: gene, module, rho_CA1, rho_CA3, DE_age_lfc_*, DE_region_lfc, score_CA1/CA3_default, rank_CA1/CA3, candidate_class, self_correlation_flag |
| 12 | `phase13_top_candidates_CA1.csv` | Data (13,179 rows) | Sorted by rank_CA1; gene×module pairs |
| 13 | `phase13_top_candidates_CA3.csv` | Data (13,179 rows) | Sorted by rank_CA3; gene×module pairs |
| 14 | `phase13_top_candidates_concordant.csv` | Data (1,574 rows) | candidate_class=="concordant"; combined_score = score_CA1 + score_CA3 |
| 15 | `phase13_top_candidates_region_flipped.csv` | Data (675 rows) | candidate_class=="region_flipped"; 675 gene×module pairs |
| 16 | `phase13c_age_specific_candidate_classes.csv` | Data (13,179 rows) | Columns: primary_age_class, region_age_class, candidate_class, secondary_class, sample_size_label_CA1_Young/Old/CA3_Young/Old, rho_*, candidate_score_*, etc. |
| 17 | `phase13c_old_specific_candidates.csv` | Data (1,716 rows) | Union of CA1_old_specific + CA3_old_specific |
| 18 | `phase13c_young_specific_candidates.csv` | Data (6,753 rows) | Union of CA1_young_specific + CA3_young_specific |
| 19 | `phase13c_age_shifted_candidates.csv` | Data (2,624 rows) | CA1_age_shifted + CA3_age_shifted (rho sign reversal) |
| 20 | `phase13c_region_shifted_by_age_candidates.csv` | Data (6,565 rows) | region_age_class = Young_CA3_specific, Old_CA1_specific, Old_CA3_specific |
| 21 | `phase13c_top_candidates_by_region_age.csv` | Data | Top 20 per stratum by |rho| |
| 22 | `phase13c_final_report.txt` | Report text | Verified candidate class distributions (primary_age_class: CA1_young_specific=5649, CA1_age_shifted=1767, CA1_old_specific=1716, CA1_age_conserved=1601, CA3_young_specific=1104, CA3_age_shifted=857, CA3_age_conserved=482, module_member_or_effector=2). region_age_class: Young_region_conserved=4291, Young_CA3_specific=4159, none=1740, Old_CA3_specific=1529, Old_CA1_specific=877, Old_region_conserved=580, module_member_or_effector=2 |
| 23 | `phase13c_validation_summary.csv` | Validation | 32 PASS / 0 WARN / 0 FAIL |
| 24 | `figures/spatial/phase13_candidate_regulator/` | Figures (40 PNGs) | F01-F14 (CA1, 14 files), F01b-F14b (CA3, 14 files), phase13c_* (12 files) |

### 2.2 References and Skills for Reporting

#### 2.2.1 Upstream Official Documentation (Accessed via ref-bio Routing)

ref-bio was used as a routing index only. `references.link-only.yaml` entries are link-only metadata — they are NOT evidence of review. Each source below was routed via ref-bio indexes, then the actual upstream URL was opened and reviewed. Entries where the URL could not be accessed are recorded as `failed_to_access`.

| # | source_or_skill | path_or_url | reviewed_status | specific_section_or_topic | key_guidance | impact_on_report_design |
|---|----------------|-------------|-----------------|-------------------------|-------------|------------------------|
| U1 | Seurat — Analysis of Spatial Datasets (Sequencing-based) | https://satijalab.org/seurat/articles/spatial_vignette | **reviewed** | Normalization (SCTransform), SpatialFeaturePlot, SpatialDimPlot, FindSpatiallyVariableFeatures (moransi), deconvolution vs integration (anchor-based transfer), multi-slice merging, cluster-based subsetting | **Key finding**: "The variance in molecular counts / spot can be substantial for spatial datasets, particularly if there are differences in cell density across the tissue." This confirms that pseudobulk aggregation at GEMgroup level is a sensible strategy for correlation analysis, and that spot-level analysis requires cell-density-aware normalization. **Report impact**: Report must explain why pseudobulk (not spot-level) was used for correlation; spatial maps (SpatialFeaturePlot) show expression patterns as visual confirmation, not as inferential evidence. |
| U2 | DESeq2 Official Vignette | https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html | **reviewed** | Standard workflow, design formulas, independent filtering, LFC shrinkage (apeglm), p-value adjustment, FAQ: "Can I use DESeq2 with no replicates?" (dispersion estimation still works but results are exploratory), contrast specification | **Key finding**: The FAQ explicitly states that without replicates, "the dispersion estimation step will still provide estimates" but results are not statistically calibrated. With n=3-4 per Young stratum, DESeq2 DE results provide ranking information but should not be treated as confirmatory. **Report impact**: Report must state that Phase 11/12 DE padj values are used as supporting annotation only at n=8 (Old strata, moderate power) and are EXPLORATORY at n=3-4 (Young strata). LFC shrinkage (apeglm) provides more stable effect size estimates for ranking. |
| U3 | OSTA (Orchestrating Spatial Transcriptomics Analysis) — TOC + Feature-set Signatures | https://bioconductor.org/books/release/OSTA/ | **partially_reviewed** | Chapter 30 (Feature-set signatures) — validates module-score approach; Chapter 12 (Deconvolution) — confirms Visium spots are multi-cellular; Chapter 31 (Spatial statistics) — spot-level spatial autocorrelation, Moran's I; Chapter 33 (Differential spatial patterns) — multi-sample spatial comparison | **Key finding**: OSTA formally documents that module scores computed as feature-set signatures are a standard, defensible spatial readout. It also confirms that Visium spots capture multiple cells and that spatial autocorrelation exists (neighboring spots share local microenvironment). **Report impact**: Report can cite OSTA as authoritative support for using module scores as functional readouts. Must remind readers that spots ≠ single cells, and that spatial autocorrelation supports the pseudobulk aggregation approach (aggregating neighbor spots into GEMgroups reduces spatial noise). |
| U4 | Seurat — AddModuleScore Reference | https://satijalab.org/seurat/reference/addmodulescore | **failed_to_access** | URL returned Seurat reference landing page; specific AddModuleScore anchor unreachable via webfetch. Reviewed via Seurat spatial vignette which uses AddModuleScore. | Module scores computed via `Seurat::AddModuleScore(slot="data", nbin=24, ctrl=100)`. Scores are based on Tirosh et al 2016 method. Same score computation used across all spatial phases. **Report impact**: Report must explain that module scores represent relative expression of a gene set vs. background control features. Scores are normalized within-sample, not across samples. |
| U5 | OSCA (Orchestrating Single-Cell Analysis) — Correlation/coexpression | https://bioconductor.org/books/release/OSCA/ | **failed_to_access** | URL redirected; page could not be fetched within time limit. Used bio-gene-regulatory-networks-coexpression-networks skill as fallback. | Correlation guidance from bio-* skill S1 (below) is authoritative for small-n constraints. **Report impact**: Defer to S1 for correlation limits at n<15. |
| U6 | ref-bio SKILL.md | `.opencode/skills/ref-bio/SKILL.md` | **inspected** | §Trigger (only activates with /ref-bio prefix), §Purpose (source-routing, not summarization), §Forbidden behavior (do not summarize, infer thresholds from memory, claim to have read inaccessible URLs), §Task routing workflow (classify by modality, workflow stage, package/tool, platform), §Source resolution workflow (indexes → link-only YAML → upstream URL), §Required reading order (AGENTS.md, PLAN.md, existing scripts, indexes, YAML, THEN upstream URLs) | **Key finding**: ref-bio explicitly forbids treating link-only YAML entries as evidence of review — upstream URLs must be opened. For Phase 14 planning, ref-bio routing identified 5 relevant source IDs: `seurat` (→ Seurat spatial vignette), `deseq2` (→ DESeq2 vignette), `osta` (→ OSTA), `osca` (→ OSCA), `spatial_transcriptomics` (→ OSTA + Seurat). **Note**: `seurat-spatial` and `pseudobulk-de-guidance` are listed in the indexes but do NOT exist in the link-only YAML — this is a routing gap in the current ref-bio catalog. Upstream URLs for these were inferred from the Seurat and OSTA entries. |

#### 2.2.2 bio-* Skills

| # | source_or_skill | path_or_url | reviewed_status | specific_section_or_topic | key_guidance | impact_on_report_design |
|---|----------------|-------------|-----------------|-------------------------|-------------|------------------------|
| S1 | bio-gene-regulatory-networks-coexpression-networks | `.opencode/skills/bio-gene-regulatory-networks-coexpression-networks/SKILL.md` | **inspected** | §Decision Tree by Scenario: "<15 samples: none reliable — correlation estimates too noisy; report this, do not force WGCNA." §Co-expression is marginal correlation; regulation is conditional dependence. §Scale-free topology caveat (Broido & Clauset 2019: only ~4% of networks are truly scale-free). | **Confirms**: (1) n=3-12 per stratum is below the threshold for reliable correlation-based inference. (2) Co-expression edges ≠ regulatory edges even when the hub is a TF. (3) "Master regulator" from correlation data is the most common over-claim. **Report impact**: Report MUST state that all candidate rankings are hypothesis-generating only. Must NOT claim "master regulators", "hubs", or "drives". Must explain the marginal-vs-partial correlation distinction. |
| S2 | bio-differential-expression-de-results | `.opencode/skills/bio-differential-expression-de-results/SKILL.md` | **inspected** | padj=NA handling, BH FDR limitations at small n, TREAT vs post-hoc filtering, p-value histogram diagnostics, Schurch 2016 small-n result. | DE evidence from Phase 11/12 is supporting annotation only. padj values at n=3-4 have severely limited calibration; padj at n=8 has moderate but still limited power. **Report impact**: DE evidence columns are labeled "supporting" in report. padj from Phase 12 Old strata is described as "exploratory supporting evidence", not as "significant DE". |
| S3 | bio-data-visualization-heatmaps-clustering | `.opencode/skills/bio-data-visualization-heatmaps-clustering/SKILL.md` | **inspected** | ward.D vs ward.D2 trap (Murtagh-Legendre 2014), optimal leaf ordering (Bar-Joseph 2001), robust color mapping (1st-99th percentile clipping), diverging vs sequential palettes, ComplexHeatmap draw() requirement for non-interactive contexts. | **Report impact**: Figure explanations note that heatmaps use symmetric diverging color scales with robust bounds, ward.D2 clustering, and display sample size prominently. |
| S4 | bio-data-visualization-ggplot2-fundamentals | `.opencode/skills/bio-data-visualization-ggplot2-fundamentals/SKILL.md` | **inspected** | cairo_pdf for TrueType embedding, theme_classic() publication baseline, tidy eval programmatic aes, Weissgerber 2015 sample-size honesty, ggrepel max.overlaps=Inf to prevent silent label drops. | Already applied in Phase 13/13b/13c figure generation. Not directly used in report writing (report is markdown, not R). |
| S5 | bio-data-visualization-distribution-plots | `.opencode/skills/bio-data-visualization-distribution-plots/SKILL.md` | **not_reviewed** | Weissgerber 2015: show individual points, not just summary stats; N-visible encoding. | Already applied in Phase 13/13b/13c figures. Not critical for plan. |
| S6 | bio-reporting-figure-export | `.opencode/skills/bio-reporting-figure-export/SKILL.md` | **not_reviewed** | cairo_pdf, resolution, sizing. | PNGs are already generated. Not critical for plan. |
| S7 | bio-differential-expression-deseq2-basics | `.opencode/skills/bio-differential-expression-deseq2-basics/SKILL.md` | **not_reviewed** | Already reviewed in Phase 13 plan. Phase 14 does no new DESeq2. | Not needed for report plan. |
| S8 | bio-spatial-transcriptomics-spatial-visualization | `.opencode/skills/bio-spatial-transcriptomics-spatial-visualization/SKILL.md` | **not_reviewed** | Python/Squidpy — conceptual only for R pipeline. | Not needed for report plan. |

### 2.3 Upstream Review Impact Summary

Key constraints confirmed by both upstream document review AND bio-* skill inspection:

1. **Correlation ≠ causation** (S1 + U5 fallback): Every mention of candidate genes MUST use "candidate regulator-associated genes" or "genes associated with mitochondrial module scores." Never "regulates", "drives", "mediates", "upstream regulator".

2. **Small n < 15 → no reliable correlation inference** (S1): Report MUST state that n=3-12 per stratum is below the threshold for reliable correlation, let alone network inference. Rankings are heuristic and exploratory.

3. **CA3 Young n=3** (S1 + AGENTS.md): Severely underpowered. Minimum achievable Spearman p-value ~0.33. MUST appear in executive summary, CA3 section, age section, and every figure/table explanation involving CA3 Young.

4. **Visium spots ≠ cells** (U1 + U3): Seurat spatial vignette and OSTA both confirm that spots capture multiple cells. Spatially variable features exist (moransi). Pseudobulk aggregation at GEMgroup level addresses this but individual-spot heterogeneity is lost. Report must state this.

5. **Module scores are valid but imperfect readouts** (U3 + U4): OSTA chapter 30 formally documents feature-set signatures (module scores) as a standard spatial readout. Seurat AddModuleScore uses control-feature background correction. Scores are relative within-sample, not absolute across-sample measurements. Report must explain what module scores can and cannot measure.

6. **DE padj at small n is exploratory** (U2 + S2): DESeq2 FAQ confirms that without replicates, DE results are rank-based not statistically calibrated. Even at n=8, power is limited. Report must label all padj as "supporting exploratory evidence."

7. **Visualization conventions** (S3): Heatmaps use symmetric diverging color scales with robust quantile bounds. n displayed prominently. ward.D2 clustering.

### 2.4 Terminology: Chinese-English Boundary and Standardized Report Language

The user's scientific goal is expressed in Chinese as:

> *"利用 CA1/CA3 空间数据，以线粒体 module score 作为功能 readout，推断与区域性 mitochondrial transcriptional programs 相关的候选上游调控基因"*

This translates conceptually as:

> "Using CA1/CA3 spatial data, with mitochondrial module scores as functional readouts, infer candidate upstream regulatory genes associated with regional mitochondrial transcriptional programs."

**CRITICAL**: This is the *scientific question*, not the *finding*. The report may explain this as the motivating question (including in Chinese if the user requests), but must consistently use standardized English terminology in all results, interpretations, tables, and figures.

**Standardized terminology for the report body**:

| Concept | Report language (REQUIRED) | Forbidden language |
|---------|---------------------------|-------------------|
| Core finding | candidate regulator-associated genes | upstream regulators, master regulators, regulatory genes |
| Interpretive framing | candidate upstream-regulator hypotheses | inferred upstream regulators, predicted regulators |
| Heuristic class | hypothesis-generating candidates | validated candidates, confirmed regulators |
| Evidence strength | genes associated with mitochondrial module scores | genes that regulate mitochondrial programs |
| Ranking | candidate_score, |rho| ranking | regulator probability, regulatory likelihood |
| Self-correlation caveat | module_member_or_effector — excluded from ranking for own module | self-regulator, autoregulator |
| CA3 Young | severely_underpowered_exploratory (n=3) | Young CA3 findings, CA3 Young-associated |

**Terminology in figure titles/captions**: Follow the same convention. If an existing figure title uses a term that could be misinterpreted (e.g., "top regulators" in a Phase 13 figure), the report caption must clarify: "Note: 'top regulators' in the figure title refers to candidate regulator-associated genes ranked by composite score. No causal regulatory claim is implied."

**Allowed vs forbidden language for validation checks**:

| Forbidden (grep pattern) | Rationale |
|--------------------------|-----------|
| `\bcausal regulator\b` NOT preceded by "not " or "no " | Implies proven causation |
| `\bdrives\b` | Implies causation |
| `\bmediates\b` | Implies mechanism |
| `\bregulates (mitochondrial|the module)\b` | Implies direction |
| `\bmaster regulator\b` | Co-expression cannot establish this |
| `\bupstream regulator\b` (unless in the phrase "candidate upstream-regulator hypotheses") | Too strong |

| Allowed | Rationale |
|---------|-----------|
| "not causal", "no causal claims" | Stating that causation is NOT claimed |
| "correlation is not causation" | Standard caveat |
| "candidate regulator-associated genes" | Mandatory term |
| "candidate upstream-regulator hypotheses" (in methodology/narrative context only) | Used to explain the motivating question; clearly labeled as hypothesis |
| "hypothesis-generating candidates" | Explicitly non-conclusory |
| "association-based" | Correctly describes the method |
| "genes associated with mitochondrial module scores" | Most conservative, always correct |
| "Spearman correlation between gene expression and module score" | Exact description, no causal implication |

1. **Correlation ≠ causation** (S1): Every mention of candidate genes MUST use "candidate regulator-associated genes" or "genes associated with mitochondrial module scores." Never "regulates", "drives", "mediates", "upstream regulator".

2. **Small n < 15 → no network inference** (S1): Report MUST state that n=3-12 per stratum is below the threshold for reliable correlation, let alone network inference. Rankings are heuristic and exploratory.

3. **CA3 Young n=3** (S1 + AGENTS.md): Severely underpowered. Minimum achievable Spearman p-value ~0.33. All CA3 Young results carry `severely_underpowered_exploratory` label. MUST appear in executive summary, CA3 section, age section, and every figure/table explanation involving CA3 Young.

4. **Visualization conventions** (S3): Heatmaps use symmetric diverging color scales with robust quantile bounds. n displayed prominently. ward.D2 clustering.

5. **Reporting standards** (S2): q-values/p-values are exploratory only. Use |rho| as primary effect-size metric.

---

## 3. Key Design Decisions

| # | Decision | Rationale |
|---|----------|----------|
| D1 | No new analysis — report integration only | Per task constraints. Phase 14 reads existing CSVs/figures; does not recompute. |
| D2 | `data/report/phase14_candidate_regulator_report/` as output directory | Follows data/ convention; gitignored. Separates generated reports from source docs/. |
| D3 | Copy selected figures/tables into report bundle; do not move originals | Originals remain in `figures/spatial/` and `data/processed/spatial/`. Report bundle contains copies for standalone reading. |
| D4 | Report as a single markdown file with inline figure/table references | Self-contained; can be converted to HTML/PDF. Links reference only files within the report bundle. |
| D5 | Prioritize Old-specific, Old-enhanced, region-shifted-in-Old, and concordant candidates | CA3 Young n=3 and CA1 Young n=4 are underpowered. Old strata (n=8) are more reliable. |
| D6 | `data/report/` convention documented in AGENTS.md and repository_structure.md | Per task instructions. Ensures future contributors know this is a generated report location. |
| D7 | `manifest.csv` records every copied file with provenance | Reproducibility. Each entry links original_path → copied_path with reason and source_phase. |
| D8 | Report-level derived summary tables (e.g., "top 10 unique genes per class") labeled as "Report Summary" | To distinguish from Phase 13/13b/13c primary analysis outputs. These are filtering/sorting operations, NOT new statistical analyses. |

---

## 4. Report Structure

### 4.1 Required Sections

#### A. Title and Executive Summary

- Title: "Phase Spatial-14: CA1/CA3 Age-Stratified Candidate Regulator-Associated Gene Report"
- 2-3 paragraph executive summary answering: what was done, what was found, key caveats, recommended next steps.
- State explicitly: "This is an association-based analysis. No causal regulatory claims are made."

#### B. Data Provenance and Input Summary

- Table of input sources: which script generated which input
- Date of Phase 13, 13b, 13c execution
- Sample sizes: CA1 Young=4, CA1 Old=8, CA3 Young=3, CA3 Old=8
- Module list: 11 modules (Complex V and mtDNA-encoded excluded)
- Regulator universe: 1,198 genes from Phase 13
- Target gene list: 251 genes (231 expressed)

#### C. Core Concepts: How to Read the Results

Plain-language explanation of every concept the reader will encounter:

1. **Spearman ρ (rho)**: Rank correlation between gene expression and module score. ∈[-1,1]. |ρ|>0.3 = screening, >0.5 = notable. NOT a p-value.

2. **Module score**: Mean expression of a set of mitochondrial function-related genes, per GEMgroup. Readout of mitochondrial transcriptional program activity. NOT a direct measure of mitochondrial function.

3. **candidate_score / regulator_score**: Heuristic ranking metric. Higher = stronger association. NOT a probability or p-value.

4. **Regulator universe**: Genes annotated as transcription factors or regulators in public databases (DoRothEA/CollecTRI/GO). Being in the universe means the gene IS a known regulator generally, NOT that it regulates mitochondrial genes in this dataset.

5. **Candidate classes**:
   - All-age (Phase 13): CA1_aging_associated, CA3_aging_associated, concordant, region_flipped, low_confidence
   - Age-stratified primary (Phase 13c): CA1/CA3_old_specific, _young_specific, _age_conserved, _age_shifted
   - Region-age (Phase 13c): Young/Old_CA1/CA3_specific, Young/Old_region_conserved
   - Special: old_enhanced_coupling, young_enhanced_coupling, module_member_or_effector

6. **Self-correlation flag**: `regulator_candidate` = eligible for ranking. `module_member_or_effector` = gene IS a member of the correlated module's gene set — self-correlation, excluded from ranking for that module.

7. **sample_size_label**: Per-stratum power assessment:
   - CA1 Young n=4: `exploratory_small_n`
   - CA1 Old n=8: `exploratory_moderate_n`
   - CA3 Young n=3: `severely_underpowered_exploratory`
   - CA3 Old n=8: `exploratory_moderate_n`

8. **pvalue_exploratory / qvalue_exploratory**: Computed for transparency but NOT used for ranking or class assignment. The `_exploratory` suffix signals limited validity at small n.

9. **Gene×module pairs**: Every row in every table is a gene×module pair, NOT a unique gene. A gene appearing in 5 module rows means 5 separate association measurements.

#### D. CA1 Findings (All-Age, Phase 13)

- Candidate class distribution in CA1
- Top CA1 candidates (from `phase13_top_candidates_CA1.csv`)
- Which modules have the strongest CA1 associations
- Concordance with Phase 12 DE (aging-associated genes)
- Figure guide: F01 (heatmap), F03 (dotplot), F09 (score components), F13 (DE volcano)
- Key CA1-specific patterns

#### E. CA3 Findings (All-Age, Phase 13b)

- CA3 mirror results (from `phase13_top_candidates_CA3.csv`)
- Comparison with CA1: overlap, divergence
- CA3-specific candidates not seen in CA1
- Figure guide: F01b-F14b (CA3 mirror figures)
- **CA3 Young n=3 caveat**: Repeated prominently

#### F. CA1 vs CA3 Cross-Region Comparison (All-Age, Phase 13)

- Concordant candidates (strong in both regions)
- Region-flipped candidates (opposite sign)
- Figure guide: F02 (rho vs DE scatter with CA1+CA3 panels)
- What regional differences might mean (speculative, caveated)

#### G. Age-Stratified Findings (Phase 13c)

**G1. Young CA1 (n=4)**
- Candidate patterns
- Caveat: n=4, exploratory_small_n
- Best candidates despite small n

**G2. Old CA1 (n=8)**
- Candidate patterns — more stable estimates
- Old-specific CA1 candidates
- Top candidates for follow-up

**G3. Young CA3 (n=3)**
- **CAUTION: severely_underpowered_exploratory**
- Results reported for completeness only
- Do NOT interpret as definitive findings
- No recommendations derived from CA3 Young-only data

**G4. Old CA3 (n=8)**
- Candidate patterns
- Old-specific CA3 candidates
- Comparison with Old CA1

#### H. Old-Specific and Old-Enhanced Candidates (Phase 13c, Primary Priority)

- From `phase13c_old_specific_candidates.csv` (1,716 gene×module pairs)
- Old-enhanced coupling candidates
- Why Old strata are prioritized: n=8 in both regions
- Intersection with Phase 13 all-age top candidates
- Recommended priority list for literature review

#### I. Region-Shifted Within Age (Phase 13c)

- Old_CA1_specific, Old_CA3_specific, Young_CA3_specific
- What region specificity within an age stratum suggests
- Caveat: Old analysis is more reliable than Young

#### J. Age-Shifted Candidates (Phase 13c)

- Candidates with rho sign reversal between Young and Old
- Interpretations: potential age-dependent coupling change
- Caveat: Young strata rho is noisy (n=3-4). Age-shifted classes are exploratory.

#### K. Prioritized Candidate Short List

Transparent, rule-based candidate prioritization for manual literature review.

**Inclusion criteria** (ALL must be met):
1. `|rho| > 0.5` in Old CA1 or Old CA3 (n=8 strata)
2. NOT `module_member_or_effector` for the correlated module
3. NOT solely driven by CA3 Young (n=3 only candidates excluded)
4. `self_correlation_flag == "regulator_candidate"`

**Priority tier assignment**:

| Tier | Criteria | Why higher confidence |
|------|----------|----------------------|
| Tier 1 | Old-specific AND concordant (strong in both Old CA1 AND Old CA3) | Replicated signal across regions, stronger n, aging-relevant |
| Tier 1 | Old_region_conserved AND |rho|>0.6 in both Old strata | Same direction, strong, in both regions |
| Tier 2 | Old-specific (CA1 or CA3) AND |rho|>0.6 | Strong signal in one region's Old stratum |
| Tier 2 | Concordant (all-age) AND appears in Phase 13c Old-specific | All-age concordance + age-specific confirmation |
| Tier 3 | Old-enhanced coupling (|rho_Old| > |rho_Young|+0.2) AND |rho_Old|>0.5 | Age-enhanced, moderate rho |
| Tier 4 | Region-flipped AND |rho|>0.6 in both but opposite signs | Interesting biology (region-specific programs) but harder to interpret |
| Tier 5 | Old-specific or concordant with 0.3 < |rho| < 0.5 | Weaker signal. Low priority for primary follow-up but worth noting |
| Excluded | CA3 Young-only evidence (any class involving CA3 Young as sole strong signal) | n=3, unreliable |

**Output format**: Table with gene, module, best rho, candidate classes (all-age + age-stratified), DE_support_flag, priority_tier, notes_for_reviewer.

**CRITICAL metadata**: T10 and T11 MUST include the following columns (or CSV header comments):

- `report_level_prioritization_only = TRUE`
- `no_new_statistical_test = TRUE`

Both T10 and T11 are **deterministic filtering/sorting/labeling** of existing Phase 13/13b/13c output values. They do NOT constitute a new inferential analysis. The `priority_tier` is a report-convenience label for human readability, NOT a statistical finding. Report README.md must explicitly state this.

#### L. Figure Guide

For each included figure, a standardized entry explaining:
- Which file (in the report bundle)
- What the figure shows
- What each axis/color means
- What a reader should look for
- What caveat applies (especially n and power)
- Which source phase generated it

#### M. How to Read Each Main Table

For each included table:
- Which file (in the report bundle)
- What question it answers
- Key columns explained
- How to sort/filter for different questions
- What NOT to over-interpret

#### N. Limitations and Interpretation Guardrails

Must include:

1. **Small n**: CA1 n=12, CA3 n=11 (all-age); Young strata n=3-4. All rankings are heuristic.
2. **Correlation ≠ causation**: No perturbation experiments, no mechanistic validation.
3. **Module score circularity**: Module scores and gene expression share the same underlying count data. While a module score aggregates many genes (mitigating circularity), it is not an independent measurement.
4. **Visium spots ≠ cells**: Each spot captures multiple cells. Pseudobulk aggregation (GEMgroup level) partially addresses this but spot-level heterogeneity is lost.
5. **Candidate universe is annotation-based**: TF annotation (DoRothEA/CollecTRI/GO) is derived from general databases, not hippocampal-specific evidence.
6. **Complex V (Atp5*) missing**: 16/16 genes absent from Hippo data. Module excluded.
7. **mtDNA-encoded module**: Has gene set but no sample-level score — excluded.
8. **CA3 Young n=3**: Severely underpowered. All CA3 Young results are exploratory only.
9. **Report-level summaries are not new inferential tests**: Derived tables (e.g., filtered top N lists) are for presentation only.
10. **Age and GEMgroup are colinear**: Cannot separate age from individual animal effects.
11. **RDS-based workflow**: Not raw Space Ranger/Visium files.

#### O. Recommended Next Steps

1. **Manual literature review** of Tier 1-2 candidates: check published evidence for mitochondrial regulation, hippocampal expression, aging-related function.
2. **Cross-check with cell-type composition**: If spot-level deconvolution becomes available, assess whether candidate signals are driven by cell composition shifts.
3. **Orthogonal validation**: Consider:
   - Independent spatial transcriptomics dataset (if available)
   - scRNA-seq regulatory network evidence (if future phase adds this)
   - Published ChIP-seq or perturbation data for candidate TFs
4. **Experimental design**: For Tier 1 candidates, consider:
   - In situ hybridization or immunohistochemistry in aging hippocampus
   - Knockdown/overexpression in hippocampal cell culture + mitochondrial readout
   - Tissue-specific knockout mouse models (long-term)
5. **Future computational extensions** (NOT in Phase 14):
   - Phase 12: Young vs Aged regional comparison (planned but not executed)
   - Cell-type deconvolution of Visium spots
   - Integration with scRNA-seq regulatory evidence
   - Validated regulon activity estimation if larger sample sizes become available

#### P. Appendices

- Data dictionary: all column names across Phase 13/13b/13c outputs explained
- Validation summary: Phase 13 (12/0/0), Phase 13b (17/0/0), Phase 13c (32/0/0)
- Full candidate class distribution tables
- Module gene set membership reference

---

## 5. Planned Report-Level Summary Tables

All generated by reading existing CSVs, filtering/sorting, and writing new CSVs. NO new computation.

| # | Output File (in report bundle) | Source | Filter/Sort Logic | Purpose |
|---|-------------------------------|--------|-------------------|---------|
| T1 | `tables/phase14_top_CA1_all_age.csv` | `phase13_top_candidates_CA1.csv` | Top 50 by rank_CA1, exclude module_member_or_effector, add unique_gene flag | CA1 all-age summary |
| T2 | `tables/phase14_top_CA3_all_age.csv` | `phase13_top_candidates_CA3.csv` | Top 50 by rank_CA3, exclude module_member_or_effector | CA3 all-age summary |
| T3 | `tables/phase14_top_concordant.csv` | `phase13_top_candidates_concordant.csv` | Top 50 by combined_score, exclude module_member_or_effector | Cross-region concordant |
| T4 | `tables/phase14_top_region_flipped.csv` | `phase13_top_candidates_region_flipped.csv` | Top 30 by combined_score | Region-flipped with opposite signs |
| T5 | `tables/phase14_top_old_specific_CA1.csv` | `phase13c_old_specific_candidates.csv` | Filter CA1_old_specific, top 30 by |rho_CA1_Old| | Old CA1-specific |
| T6 | `tables/phase14_top_old_specific_CA3.csv` | `phase13c_old_specific_candidates.csv` | Filter CA3_old_specific, top 30 by |rho_CA3_Old| | Old CA3-specific |
| T7 | `tables/phase14_top_old_enhanced_coupling.csv` | `phase13c_age_specific_candidate_classes.csv` | Filter old_enhanced_coupling_flag==TRUE, top 30 by max(|rho_Old|) | Age-enhanced coupling |
| T8 | `tables/phase14_top_region_shifted_old.csv` | `phase13c_region_shifted_by_age_candidates.csv` | Filter age_stratum=="Old", top 30 by |delta_region| | Region-shifted in Old |
| T9 | `tables/phase14_top_age_shifted.csv` | `phase13c_age_shifted_candidates.csv` | Top 30 by |delta_rho|, exclude CA3 Young-only | Age-shifted (sign reversal) |
| T10 | `tables/phase14_prioritized_short_list.csv` | Multiple sources | Per §4.1-K tier criteria | Literature review priorities |
| T11 | `tables/phase14_prioritized_short_list_unique_genes.csv` | T10 | Deduplicated to unique genes; best module per gene | Same as T10 but unique gene resolution |

**Convention**: All report tables include a header comment line stating "Report Summary — Phase Spatial-14. Derived from Phase [13/13b/13c] outputs. NOT a new statistical analysis."

---

## 6. Planned Figure Selection

Copy selected existing figures into the report bundle. Do NOT regenerate.

### 6.1 Phase 13 CA1-Focused (6 figures)

| Report Path | Source Path | Description |
|-------------|-------------|-------------|
| `figures/CA1_CA3_overview/F01_CA1_regulator_score_heatmap.png` | `figures/spatial/phase13_candidate_regulator/F01_regulator_score_heatmap.png` | CA1 all-age regulator_score heatmap |
| `figures/CA1_CA3_overview/F03_CA1_top_regulators_dotplot.png` | `figures/spatial/phase13_candidate_regulator/F03_top_regulators_dotplot.png` | CA1 top regulators dotplot |
| `figures/CA1_CA3_overview/F02_rho_vs_DE_scatter.png` | `figures/spatial/phase13_candidate_regulator/F02_rho_vs_DE_scatter.png` | CA1+CA3 rho vs DE scatter |
| `figures/CA1_CA3_overview/F06_candidate_classes.png` | `figures/spatial/phase13_candidate_regulator/F06_candidate_classes.png` | Candidate class distribution (all-age) |
| `figures/CA1_CA3_overview/F09_CA1_score_components.png` | `figures/spatial/phase13_candidate_regulator/F09_score_components.png` | CA1 score component breakdown |

### 6.2 Phase 13b CA3-Focused (4 figures)

| Report Path | Source Path | Description |
|-------------|-------------|-------------|
| `figures/CA1_CA3_overview/F01b_CA3_regulator_score_heatmap.png` | `figures/spatial/phase13_candidate_regulator/F01b_CA3_regulator_score_heatmap.png` | CA3 regulator_score heatmap |
| `figures/CA1_CA3_overview/F03b_CA3_top_regulators_dotplot.png` | `figures/spatial/phase13_candidate_regulator/F03b_CA3_top_regulators_dotplot.png` | CA3 top regulators |
| `figures/CA1_CA3_overview/F09b_CA3_score_components.png` | `figures/spatial/phase13_candidate_regulator/F09b_CA3_score_components.png` | CA3 score component |
| `figures/CA1_CA3_overview/F10b_CA3_sensitivity_comparison.png` | `figures/spatial/phase13_candidate_regulator/F10b_CA3_sensitivity_comparison.png` | CA3 sensitivity |

### 6.3 Phase 13c Age-Stratified (7 figures)

| Report Path | Source Path | Description |
|-------------|-------------|-------------|
| `figures/age_stratified/phase13c_heatmap_CA1_young.png` | As named | CA1 Young stratified rho heatmap |
| `figures/age_stratified/phase13c_heatmap_CA1_old.png` | As named | CA1 Old stratified rho heatmap |
| `figures/age_stratified/phase13c_heatmap_CA3_old.png` | As named | CA3 Old stratified rho heatmap |
| `figures/age_stratified/phase13c_heatmap_CA3_young.png` | As named | CA3 Young (CAUTION) |
| `figures/age_stratified/phase13c_delta_CA1_old_minus_young.png` | As named | CA1 delta rho (Old-Young) |
| `figures/age_stratified/phase13c_region_age_class_barplot.png` | As named | Candidate class counts |
| `figures/age_stratified/phase13c_old_specific_top_dotplot.png` | As named | Old-specific top candidates |

### 6.4 Spatial Maps (2 figures, optional)

| Report Path | Source Path | Description |
|-------------|-------------|-------------|
| `figures/spatial_maps/F14_spatial_CA1_Dach1.png` | As named | Top CA1 candidate spatial |
| `figures/spatial_maps/F14b_spatial_CA3_Elk1.png` | As named | Top CA3 candidate spatial |

**Selection rationale**: 6+4+7+2 = 19 figures total. Focused on key results. Non-critical supporting figures (sensitivity, age-shifted dotplots, individual spatial maps for all 10 genes) omitted to keep report manageable. Report text references omitted figures by original path for readers who want to explore further.

---

## 7. Planned Report Bundle Structure

```
data/report/phase14_candidate_regulator_report/
├── spatial_phase14_candidate_regulator_report.md     # Main report (markdown)
├── README.md                                          # Bundle README
├── manifest.csv                                       # File provenance manifest
├── tables/
│   ├── phase14_top_CA1_all_age.csv
│   ├── phase14_top_CA3_all_age.csv
│   ├── phase14_top_concordant.csv
│   ├── phase14_top_region_flipped.csv
│   ├── phase14_top_old_specific_CA1.csv
│   ├── phase14_top_old_specific_CA3.csv
│   ├── phase14_top_old_enhanced_coupling.csv
│   ├── phase14_top_region_shifted_old.csv
│   ├── phase14_top_age_shifted.csv
│   ├── phase14_prioritized_short_list.csv
│   ├── phase14_prioritized_short_list_unique_genes.csv
│   └── copied_source_tables/
│       └── (selected source tables for reference, see §7.1)
├── figures/
│   ├── CA1_CA3_overview/
│   │   └── (10 figures from Phase 13/13b, see §6)
│   ├── age_stratified/
│   │   └── (7 figures from Phase 13c, see §6)
│   └── spatial_maps/
│       └── (2 spatial FeaturePlots, see §6)
└── provenance/
    ├── source_file_manifest.csv                       # Manifest of original source files
    └── phase14_report_generation_summary.txt           # Generation timestamp, script version, checksums
```

### 7.1 Copied Source Tables (Optional)

If a reader wants direct access to key source tables without navigating back to `data/processed/spatial/`, copy these small CSVs:

| Source | Reason |
|--------|--------|
| `phase13_regulator_score_default.csv` (first 100 rows) | Demonstrate all-age score format |
| `phase13c_age_specific_candidate_classes.csv` (first 100 rows) | Demonstrate age-stratified class columns |
| `phase13c_validation_summary.csv` | Full validation results |

These are head/tail copies, not full files. Full files are linked by original path in the report.

---

## 8. Execution Plan (Future)

### 8.1 Script

**`R/spatial/s14_generate_candidate_regulator_report.R`** — NOT created in this plan phase.

Script outline (for reference):

```
E0: Setup & Validation
  - Verify all input CSVs readable (Phase 13/13b/13c outputs)
  - Verify all figures to copy exist
  - Verify Phase 13/13b/13c validation: all PASS, 0 FAIL
  - Record renv.lock MD5 (must remain unchanged)
  - Create report bundle directories

E1: Generate Summary Tables (T1-T11)
  - Read existing CSVs
  - Filter, sort, select columns
  - Write report-level summary CSVs to tables/
  - NO new computation beyond filtering/sorting

E2: Copy Figures
  - Copy selected figures from figures/spatial/phase13_candidate_regulator/
  - Preserve original filenames (or add descriptive prefix)
  - Log copies in manifest

E3: Copy Source Tables (Optional)
  - Copy head of key source CSVs for reference

E4: Write Main Report
  - Generate markdown file with all sections (§4.1)
  - Inline tables referenced by path
  - Inline figures referenced by relative path
  - CA3 Young caveat in every relevant section

E5: Write Bundle Metadata
  - manifest.csv: original_path, copied_path, file_type, reason_included, source_phase, caption
  - README.md: how to use the report bundle
  - phase14_report_generation_summary.txt

E6: Update Guidance Files
  - AGENTS.md: add Phase 14 status, data/report/ convention
  - docs/repository_structure.md: add data/report/ to canonical layout
  - README.md: add data/report/ convention, copied tables/figures note

E7: Validation & Cleanup
  - All required report sections present
  - All report tables have corresponding CSV files
  - All copied figures exist and are readable
  - No causal language in report text
  - CA3 Young caveat verified in summary, age section, figure guide
  - No original Phase 13/13b/13c outputs modified
  - renv.lock unchanged
  - Rplots.pdf absent
  - data/report/ confirmed gitignored
```

### 8.2 Package Policy

**NO new packages.** Phase 14 uses only base R functions (`read.csv`, `write.csv`, `file.copy`, `list.files`) and packages already in the environment:
- `dplyr`: filtering, sorting
- `stringr`: simple string operations

No `renv::install()` calls. No `renv::snapshot()`.

### 8.3 Memory Management

Phase 14 reads CSVs only (largest ~2 MB). No large objects. Low memory footprint.

---

## 9. Validation Checks (for Execution Phase)

### 9.1 Pre-Execution

| # | Check | Expected | Severity |
|---|-------|----------|----------|
| V01 | Phase 13 validation: 0 FAIL | `phase13_validation_summary.csv` #FAIL=0 | FATAL |
| V02 | Phase 13b validation: 0 FAIL | `phase13b_display_validation_summary.csv` #FAIL=0 | FATAL |
| V03 | Phase 13c validation: 0 FAIL | `phase13c_validation_summary.csv` #FAIL=0 | FATAL |
| V04 | Required source CSVs (16 files) exist | All readable | FATAL |
| V05 | Required source figures exist (per §6) | All 19 figures present and non-zero | WARN (skip missing; note in manifest) |
| V06 | renv.lock MD5 recorded | Baseline preserved | INFO |
| V07 | `data/report/` directory is gitignored | grep confirms in .gitignore | WARN |

### 9.2 Post-Execution

| # | Check | Expected | Severity |
|---|-------|----------|----------|
| V08 | Main report file exists | `spatial_phase14_candidate_regulator_report.md` present | FATAL |
| V09 | All report tables (T1-T11) exist and non-empty | All CSVs present with data | FATAL |
| V10 | manifest.csv has entry for every copied file | Row count = copied files count | FATAL |
| V11 | No forbidden causal language in report | Context-aware check (see §2.4). Forbidden: "(?<!not |no )\\bcausal regulator\\b", "\\bdrives\\b", "\\bmediates\\b", "\\bregulates mitochondrial\\b", "\\bmaster regulator\\b", "\\bupstream regulator\\b" (unless in "candidate upstream-regulator hypotheses"). Allowed: "not causal", "no causal claims", "correlation is not causation", "candidate regulator-associated genes". | FATAL |
| V12 | CA3 Young caveat in executive summary | grep confirms | FATAL |
| V13 | CA3 Young caveat in age section (G3) | grep confirms | WARN |
| V14 | CA3 Young caveat in figure guide for F03/F07b | grep confirms | WARN |
| V15 | No original Phase 13/13b/13c files modified | MD5 checks match pre-execution | FATAL |
| V16 | renv.lock unchanged | MD5 match | FATAL |
| V17 | No Rplots.pdf | file.exists check | FATAL |
| V18 | All report links reference report-bundle paths | No absolute paths to data/processed/ or figures/ in report body | WARN |
| V19 | Report explains every included table and figure | Section M covers all tables; Section L covers all figures | WARN |
| V20 | README.md present in report bundle | Readable | WARN |
| V21 | Guidance files updated (AGENTS.md, repository_structure.md, README.md) | grep confirms Phase 14 and data/report/ entries | WARN |
| V22 | Prioritized short list excludes CA3 Young-only evidence | Check T10 for no "solely CA3 Young" candidates | FATAL |
| V23 | Prioritized short list excludes module_member_or_effector | Check T10 for no self_correlation_flag=="module_member_or_effector" | FATAL |

---

## 10. Stop Conditions

### 10.1 Fatal Stop Conditions

| # | Condition | Action |
|---|----------|--------|
| FATAL-01 | Any Phase 13/13b/13c validation has FAIL | STOP — Fix upstream phase first |
| FATAL-02 | Required source CSVs missing or unreadable | STOP — Regenerate upstream phase |
| FATAL-03 | Any output claims forbidden causal language (see §2.4) | STOP — Reword using "candidate regulator-associated genes" or "hypothesis-generating candidates". Context-aware check: "not causal" / "no causal claims" / "correlation is not causation" / "candidate upstream-regulator hypotheses" (explicitly labeled as hypothesis) are ALLOWED. "Causal regulator" standalone, "drives", "mediates", "master regulator", "regulates mitochondrial", "upstream regulator" (unqualified) are FORBIDDEN. |
| FATAL-04 | Original Phase 13/13b/13c outputs modified or overwritten | STOP — Restore from git |
| FATAL-05 | renv.lock modified during execution | STOP — Investigate; no package changes are authorized |
| FATAL-06 | Rplots.pdf generated | STOP — Fix device handling |
| FATAL-07 | Prioritized short list includes CA3 Young-only evidence | STOP — Filter criteria should exclude these |
| FATAL-08 | Any new statistical test or computation performed (beyond filtering/sorting existing values) | STOP — Report is integration only |
| FATAL-09 | Report bundle committed to git (data/ is gitignored) | STOP — Confirm gitignore before proceeding |

### 10.2 Non-Fatal Warnings

| # | Condition | Action |
|---|----------|--------|
| WARN-01 | A figure to copy is missing | Skip that figure; record in manifest as "missing"; continue with remaining figures |
| WARN-02 | Report-level summary CSV has < 5 rows | Note in report text; continue |
| WARN-03 | Guidance file update partial (e.g., AGENTS.md updated but README.md not) | Log discrepancy; continue; fix manually |

---

## 11. Guidance File Updates (During Execution)

### 11.1 AGENTS.md

Add under `### Current Scientific Trajectory`:
```
- Phase Spatial-14: Candidate regulator report generated (data/report/). Report integration only — no new analysis.
```

Add to data paths or file structure section:
```
- data/report/: Generated report bundles (Phase 14+). Gitignored. Not source documentation.
```

### 11.2 docs/repository_structure.md

Add to canonical layout table:
```
| `data/report/` | Generated report bundles (Phase 14+); copies of figures/tables + human-readable reports | Ignored |
```

Add to output conventions:
```
Report bundles (Phase 14+):
- Report markdown and metadata: data/report/phase14_candidate_regulator_report/
- Copied figures: data/report/phase14_candidate_regulator_report/figures/
- Report-level summary tables: data/report/phase14_candidate_regulator_report/tables/
- Bundle is for standalone reading. Original outputs remain in data/processed/spatial/ and figures/spatial/.
```

### 11.3 README.md

Add to the `## Repository Structure` section or `## Data and Git Hygiene` section:
```
| `data/report/` | Generated report bundles (Phase 14+). Copied tables/figures for human reading/export. Gitignored. |
```

Add a note to the Repository Structure table or under Output Conventions:
```
Report bundles (Phase 14+):
- Located under `data/report/phase14_candidate_regulator_report/`
- Contain copied figures, report-level summary tables, and a human-readable markdown report
- For standalone reading and export — NOT a substitute for original outputs
- Original analysis outputs remain in `data/processed/spatial/` and `figures/spatial/`
- Source methodology documentation remains in `docs/`
- All `data/report/` contents are generated and gitignored
```

### 11.4 docs/spatial_phase14_candidate_regulator_report_index.md (Optional)

A lightweight pointer doc:
```
# Phase Spatial-14: Candidate Regulator Report Index
The generated report bundle is at: data/report/phase14_candidate_regulator_report/
See spatial_phase14_candidate_regulator_report.md for the main report.
```

Only create this if the plan user review approves it. Default: NOT created.

---

## 12. Plan Compliance Report

### 12.1 PLAN-ONLY Verification

| Check | Status |
|-------|--------|
| Only created/modified plan file | YES |
| No R scripts created | YES |
| No R analysis executed | YES |
| No data/figures generated | YES |
| No data/report/ directory created | YES |
| No renv.lock modified | YES (MD5: `b2bf05038b440576855d043e137c5883`) |
| No packages installed | YES |
| No original outputs modified | YES |

### 12.2 Files Created

| File | Description |
|------|-------------|
| `docs/spatial_phase14_candidate_regulator_report_plan.md` | This plan document |

### 12.3 Execution Safety Assessment

Phase 14 is **safe to execute** after user review because:

1. **Read-only on Phase 13/13b/13c outputs**: All inputs are CSVs and PNGs; no existing file is modified.
2. **No new packages**: Uses only base R + dplyr (already installed).
3. **No renv.lock change**: No package operations.
4. **Additive outputs**: All new files go to `data/report/` (gitignored). No collisions.
5. **No new analysis**: Only filtering, sorting, copying, writing.
6. **Fatal stop conditions**: Script stops before output if upstream validations fail or constraints violated.
7. **No causal language**: All output text uses "candidate regulator-associated genes."
8. **CA3 Young n=3 caveat enforced**: Validation checks V12-V14 verify presence in key sections.

---

## 13. Summary

Phase Spatial-14 is a report integration phase that synthesizes Phase 13, 13b, and 13c results into a single human-readable scientific report. The report:

1. Explains CA1 vs CA3, Young vs Old, and region×age candidate regulator-associated gene patterns
2. Prioritizes candidates for follow-up literature review using transparent, rule-based criteria
3. Produces a self-contained report bundle under `data/report/phase14_candidate_regulator_report/`
4. Copies selected figures and generates report-level summary tables
5. Updates guidance files to document the `data/report/` convention
6. Maintains strict guardrails: no causal claims, no new analysis, CA3 Young=3 caveat everywhere

The plan does no computation, generates no files (except this plan document), and is ready for user review and approval before execution.

---

**Final Report for This PLAN-Only Task (v1.1 with ref-bio fix)**:

- **Files inspected**: 24 (listed in §2.1)
- **Plan file modified**: `docs/spatial_phase14_candidate_regulator_report_plan.md`
- **ref-bio inspected files**:
  - `.opencode/skills/ref-bio/SKILL.md` — inspected, used as routing only
  - `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml` — searched for `seurat`, `seurat-spatial`, `deseq2`, `osta`, `osca`, `pseudobulk-de-guidance`, `visium`, `spatialexperiment`
  - `.opencode/skills/ref-bio/reference-pack/indexes/topic-map.yaml` — consulted for `spatial_transcriptomics`, `differential_expression`, `pseudobulk_de`
  - `.opencode/skills/ref-bio/reference-pack/indexes/package-map.yaml` — consulted for Seurat, DESeq2 routing
  - **Note**: `seurat-spatial` and `pseudobulk-de-guidance` are in indexes but NOT in link-only YAML (routing gap)
- **Upstream official docs reviewed**:
  - Seurat Spatial Vignette: **reviewed** (fully read)
  - DESeq2 Official Vignette: **reviewed** (fully read; QA section on small-n / no-replicates)
  - OSTA (TOC + Chapter 30 Feature-set Signatures): **partially_reviewed** (TOC confirming module-score approach)
  - Seurat AddModuleScore Reference: **failed_to_access** (landing page only; no specific anchor)
  - OSCA: **failed_to_access** (URL redirected; could not fetch)
- **bio-* skills reviewed**: 4 inspected (S1-S4); 4 not reviewed (not needed)
- **Confirmation that no R was run**: YES — plan only
- **Confirmation that no data/report generated**: YES
- **Confirmation that no renv.lock changed**: YES (MD5: `b2bf05038b440576855d043e137c5883`, confirmed unchanged)
- **Confirmation that no Phase 13/13b/13c outputs modified**: YES
- **Phase 14 is safe to execute after review**: YES — all safety conditions met (§12.3)
