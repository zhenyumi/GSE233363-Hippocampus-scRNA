# Phase Spatial-13c: Age-Stratified Candidate Regulator-Associated Gene Discovery Plan

**Version**: 1.0
**Date**: 2026-06-05
**Status**: PLAN-ONLY (no execution; no R scripts; no data/figures generated)
**Prerequisite**: Phase 13 executed and validated (12 PASS / 0 FAIL / 0 WARN), Phase 13b executed and validated
**Relationship**: Phase 13c is a NEW age-stratified candidate regulator discovery plan. It does NOT modify Phase 13 or Phase 13b core outputs. All new outputs use `phase13c_` prefix and merge into existing Phase 13 directories.

---

## 1. Background and Rationale

### 1.1 What Phase 13/13b Covers

Phase 13 (`s13_candidate_regulator_discovery.R`) computed CA1 and CA3 symmetrically:

- `rho_CA1` and `rho_CA3` — Spearman correlation between gene expression and mitochondrial module scores across all Young+Old GEMgroups per region (CA1 n=12, CA3 n=11)
- `score_CA1_default` and `score_CA3_default` — composite regulator scores
- `candidate_class` — includes `CA1_aging_associated` (3,027 pairs), `CA3_aging_associated` (3,290 pairs), `concordant` (1,573), `region_flipped` (674), `low_confidence` (4,608), `target_readout` (6)

Phase 13b (`s13b_ca3_symmetric_display.R`) added CA3 mirror figures and top-candidate CSVs.

### 1.2 What Phase 13/13b Does NOT Cover

Phase 13 and 13b do **NOT** create:

- Young-only or Old-only `regulator_score` tables
- Age-stratified Spearman rho within CA1 or CA3
- Young-specific or Old-specific candidate classes
- Age-shifted candidate identification (rho change direction between Young and Old)

The age-stratified heatmaps in Phase 13 (F07, F08, F07b, F08b) are **display-only exploratory figures**. They select genes by `rank_CA1` or `rank_CA3` (all-age ranking) and show age-split ρ values — they do NOT perform formal age-stratified candidate discovery, do NOT compute per-age rho from scratch, and do NOT classify candidates by age specificity.

### 1.3 What Phase 13c Adds

Phase 13c performs **formal age-stratified candidate regulator-associated gene discovery**:

1. Compute Spearman rho per gene × module × region × age (Young vs Old)
2. Generate age-specific candidate rankings and classes
3. Identify Old-specific, Young-specific, age-conserved, age-shifted, and region-shifted candidates
4. Produce age-stratified output tables with `phase13c_` prefix, merged into existing Phase 13 directories
5. Produce age-stratified figures
6. Update the analysis guide with a Phase 13c section

**CRITICAL**: Phase 13c does NOT change Phase 13 or 13b core outputs. It is a standalone discovery layer that reads Phase 11/12 inputs independently and produces its own age-stratified results.

---

## 2. Scientific Questions

Phase 13c addresses:

1. **Which candidate regulator-associated genes couple with mitochondrial module scores in Young CA1?** (n=4)
2. **Which in Old CA1?** (n=8)
3. **Which in Young CA3?** (n=3, severely underpowered, exploratory only)
4. **Which in Old CA3?** (n=8)
5. **Which candidates are Old-specific or Young-specific within a region?** (e.g., strong ρ in Old but absent in Young)
6. **Which candidates show age-conserved coupling?** (strong ρ in both Young and Old, same direction)
7. **Which candidates show descriptive age-shifted coupling?** (large effect-size shift in rho between Young and Old; no statistical test of rho differences is performed)
8. **Which candidates show CA1/CA3 region-shifted coupling within the same age stratum?** (e.g., a gene that couples in Old CA1 but not Old CA3)
9. **Do Old-specific candidates overlap with Phase 12 Old-vs-Young DEGs?** (supporting evidence only, not primary ranking)
10. **How does the age-stratified candidate landscape differ from the all-age Phase 13 ranking?**

---

## 3. Input Data Plan

### 3.1 Required Inputs

| # | Input | Source Path | Description |
|---|-------|-------------|-------------|
| I1 | Pseudobulk count matrix (CA1/CA3) | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/ca1_ca3_pseudobulk_counts.rds` | GEMgroup-level pseudobulk counts. Primary expression data source. |
| I2 | Module gene set audit | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/module_gene_set_final_audit.csv` | Module names and gene-set members. Self-correlation reference. |
| I3 | Module scores (all regions) | `data/processed/spatial/phase12_young_aged_region_comparison/module_scores_all_regions.csv` | GEMgroup-level mean module scores per region × age. |
| I4 | Phase 11 DE (CA3-vs-CA1) | `data/processed/spatial/phase11_ca1_ca3_de_module_coupling/deseq2_all_genes_CA3_vs_CA1.csv` | Supporting DE evidence for region-shifted candidates. |
| I5 | Phase 12 DE (Old-vs-Young CA1) | `data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA1.csv` | Supporting DE evidence for CA1 age-shifted candidates. |
| I6 | Phase 12 DE (Old-vs-Young CA3) | `data/processed/spatial/phase12_young_aged_region_comparison/deseq2_all_genes_Old_vs_Young_CA3.csv` | Supporting DE evidence for CA3 age-shifted candidates. |
| I7 | Target gene list | `docs/target_genes.csv` | 251 genes in 13 categories. Annotation only. |
| I8 | Phase 13 regulator score (default) | `data/processed/spatial/phase13_candidate_regulator/phase13_regulator_score_default.csv` | Reference ranking for comparison (NOT used as input for ranking). |
| I9 | Phase 13 validation summary | `data/processed/spatial/phase13_candidate_regulator/phase13_validation_summary.csv` | Must show 0 FAIL. |
| I10 | Phase 13b validation summary | `data/processed/spatial/phase13_candidate_regulator/phase13b_display_validation_summary.csv` | Must show 0 FAIL. |
| I11 | Sample manifest | `data/processed/spatial/phase13_candidate_regulator/pseudobulk_sample_manifest.csv` | Age/Region assignment per GEMgroup. |
| I12 | TF annotation coverage | `data/processed/spatial/phase13_candidate_regulator/phase13_tf_annotation_coverage.csv` | TF annotation for candidate genes (reuse Phase 13 annotation; do NOT recompute). |

### 3.2 Derived Inputs (Generated by Phase 13c at Runtime)

| Input | How Generated | Description |
|-------|--------------|-------------|
| Age-stratified sample lists | From manifest: split CA1 into Young (G1-G4) and Old (G9-G16); CA3 into Young (G2-G4) and Old (G9-G16) | 4 sample lists |
| Age-stratified expression matrices | Subset `ca1_ca3_pseudobulk_counts.rds` columns by age stratum | 4 log2(CPM+1) matrices |
| Age-stratified module score vectors | Subset `module_scores_all_regions.csv` rows by age stratum | 4 metadata data frames |
| Regulator universe | Reuse Phase 13 `universe_regulator` (1,198 genes) — do NOT re-derive | Gene list |

---

## 4. Sample Structure

### 4.1 Confirmed Sample Sizes (from pseudobulk_sample_manifest.csv)

| Stratum | GEMgroups | n | Interpretation Label |
|---------|-----------|---|---------------------|
| CA1 Young | G1, G2, G3, G4 | 4 | `exploratory_small_n` |
| CA1 Old | G9, G10, G11, G12, G13, G14, G15, G16 | 8 | `exploratory_moderate_n` |
| CA3 Young | G2, G3, G4 | 3 | `severely_underpowered_exploratory` |
| CA3 Old | G9, G10, G11, G12, G13, G14, G15, G16 | 8 | `exploratory_moderate_n` |
| CA1 All (Young+Old) | 12 | Phase 13 reference |
| CA3 All (Young+Old) | 11 | Phase 13 reference |

**Middle (G5-G8)**: Excluded from all Phase 13c analyses.

**G1_CA3**: Does NOT exist in the manifest. CA3 has no G1 sample. This is expected and does not affect analysis — CA3 Young is defined as G2, G3, G4 only.

### 4.2 Interpretation Labels in Outputs

Every output table and figure MUST carry the appropriate interpretation label from §4.1:

- `exploratory_small_n` for CA1 Young (n=4)
- `exploratory_moderate_n` for CA1 Old and CA3 Old (n=8)
- `severely_underpowered_exploratory` for CA3 Young (n=3)

These labels appear as:
1. A `sample_size_label` column in all CSV outputs
2. In figure titles/subtitles
3. In the updated analysis guide
4. In provenance and validation files

---

## 5. Methods

### 5.1 Primary Model: Age-Stratified Spearman Correlation

For each gene g in `universe_regulator` (1,198 genes from Phase 13), for each module m in the module list (from Phase 11 `module_gene_set_final_audit.csv`), compute Spearman ρ in each of 4 strata:

```
cor.test(
  expr[g, ca1_young_samples],
  module_score[m, ca1_young_samples],
  method = "spearman"
)
```

This produces `rho_CA1_Young`, `rho_CA1_Old`, `rho_CA3_Young`, `rho_CA3_Old` per gene × module pair.

**Expression normalization**: log2(CPM + 1) from pseudobulk counts, consistent with Phase 13.

**Module score source**: GEMgroup-level mean module scores from `module_scores_all_regions.csv` (Phase 12). Score columns are `<module>_mean`.

**Complex V exclusion**: Complex V module is EXCLUDED (16/16 Atp5* genes missing). Carried forward from Phase 11/12/13.

**mtDNA-encoded exclusion**: mtDNA-encoded module has a gene set but NO sample-level score in Phase 12 — excluded from correlation analysis. Carried forward from Phase 13.

### 5.2 Primary Ranking: Effect-Size-First Approach

**CRITICAL**: At n=3-4 in Young strata and n=8 in Old strata, correlation p-values and q-values have severely limited power and cannot distinguish moderate from strong effects reliably.

**Primary ranking metric**: `|rho|` (absolute Spearman correlation coefficient)

**Screening thresholds**:

| Threshold | Label | Interpretation |
|-----------|-------|---------------|
| `|rho| > 0.3` | `screening` | Gene passes screening threshold in this stratum |
| `|rho| > 0.5` | `notable` | Notable association at moderate n; notable-but-exploratory at small n |
| `|rho| > 0.7` | `strong` | Potentially strong association (rare at n=3-8) |

**p-value and q-value**:

- P-values from `cor.test()` are retained as `pvalue_exploratory` column only.
- Benjamini-Hochberg correction is applied per stratum for completeness.
- BH-adjusted values are reported as `qvalue_exploratory` (NOT `padj`, NOT `qval`).
- **Neither p-value nor q-value is used for candidate ranking or candidate class assignment.**
- Both columns carry the `_exploratory` suffix to signal their limited validity at small n.

**CA3 Young n=3 additional caveat**:

- `pvalue_exploratory` and `qvalue_exploratory` for CA3 Young strata must NOT be used for any screening, filtering, ranking, or candidate class assignment. With n=3, Spearman correlation p-values are mechanically constrained (the smallest achievable two-sided p-value is ~0.33 when all 3 rank pairs are perfectly monotonic). Any q-value computed from these p-values is equally uninformative.
- All CA3 Young candidate classes must carry the `severely_underpowered_exploratory` label in `sample_size_label`.

### 5.3 Primary Candidate Score

For each gene × module pair in each age stratum, the primary candidate score is simply the effect size:

```
candidate_score_<stratum> = |rho_<stratum>|
```

This is simple, transparent, and does not mix correlation with DE evidence.

**DE evidence**: DE log2FC and padj from Phase 11/12 are retained as supporting annotation columns only:
- `DE_age_lfc_<region>` — Old-vs-Young log2FC from Phase 12
- `DE_age_padj_<region>` — corresponding padj
- `DE_region_lfc` — CA3-vs-CA1 log2FC from Phase 11
- `DE_age_support_flag` — TRUE if Phase 12 padj < 0.05 in matching region (annotation only)

These columns appear in all relevant output tables but do NOT enter the primary candidate_score or candidate class assignment. They provide supporting biological context for interpretation.

**What candidate_score is NOT**:
- NOT a regulator_score (Phase 13's composite score uses CA1 vs CA3 DE and all-age rho)
- NOT a p-value or q-value
- NOT a probability of regulation
- NOT comparable across strata (different n, different score ranges)

### 5.4 Sensitivity Score (Optional, Separate Column)

If a DE-informed ranking is desired for exploratory comparison, a separate sensitivity column can be computed:

```
candidate_score_with_DE_support_sensitivity_<stratum> = |rho_<stratum>| + 0.1
```

Where the +0.1 bonus applies if the gene has `DE_age_support_flag == TRUE`. This column:

- Is named with `_sensitivity` suffix to clearly distinguish it from the primary score
- Is NOT used for primary ranking, candidate class assignment, or candidate selection
- Is described as "heuristic sensitivity score" in all outputs
- Can be used for optional sensitivity comparison figures but is NOT required

If sensitivity comparison is generated, add a single optional figure (`phase13c_sensitivity_rho_vs_rho_DE.png`) showing primary |rho| ranking vs sensitivity ranking overlap for top 20 candidates. Skip if not generated.

### 5.5 Delta Scores (Age Effect)

For each gene × module pair, compute within-region age delta:

```
delta_rho_CA1 = rho_CA1_Old - rho_CA1_Young
delta_rho_CA3 = rho_CA3_Old - rho_CA3_Young
```

These describe how the gene-module association changes with age.

**Interpretation**:
- Positive delta: association stronger (or more positive) in Old
- Negative delta: association stronger (or more positive) in Young
- |delta| > 0.5 with opposite signs → potential age-dependent coupling reversal
- |delta| > 0.3 → notable age shift

**Delta scores are descriptive only.** They are NOT used for formal ranking or candidate selection, because Young strata have n=3-4 (unreliable rho estimates).

### 5.6 Region-Shifted Within Same Age

For each gene × module pair, compute within-age region delta:

```
delta_region_Young = rho_CA1_Young - rho_CA3_Young
delta_region_Old   = rho_CA1_Old - rho_CA3_Old
```

These describe region specificity within the same age.

---

## 6. Self-Correlation Safeguards (Carried Forward from Phase 13)

### 6.1 Strategy A (Identical to Phase 13)

All rules from Phase 13 §5.3 apply:

1. **`module_member_or_effector`**: If gene g is a member of module m's gene set (per `module_gene_set_final_audit.csv`), the age-stratified ρ(g,m) is self-correlation and CANNOT be interpreted as upstream regulator evidence. Gene g is excluded from ranking for module m. Its `self_correlation_flag` is `module_member_or_effector`.

2. **`regulator_candidate`**: Gene is NOT a member of the correlated module's gene set. Eligible for ranking.

3. **`target_gene_candidate`**: Gene is in `target_genes.csv` but NOT a module member. Separate annotation — does NOT receive positive ranking weight. Reported in separate target gene readout.

4. **Low-confidence classification**: `candidate_class = low_confidence` when `|rho| < 0.3` in all strata. This is a `candidate_class` value, NOT a `self_correlation_flag` value.

5. **Target gene membership**: Does NOT add weight to ranking. Annotation column only.

### 6.2 Per-Module Exclusion

The self-correlation exclusion is per-module. A gene can be:
- `regulator_candidate` for module A (not a member)
- `module_member_or_effector` for module B (is a member)
- In the same output row, the `self_correlation_flag` refers to the specific gene × module pair

---

## 7. Candidate Classifications

### 7.1 Primary Classes (Within-Region, Age-Specific)

Classes are assigned per gene × module pair based on age-stratified rho patterns:

| Class | Criteria | Interpretation |
|-------|----------|---------------|
| `CA1_old_specific` | `|rho_CA1_Old| > 0.3` AND `|rho_CA1_Young| < 0.3` | Candidate associated with mitochondrial module in Old CA1 but NOT Young CA1 |
| `CA1_young_specific` | `|rho_CA1_Young| > 0.3` AND `|rho_CA1_Old| < 0.3` | Candidate associated with mitochondrial module in Young CA1 but NOT Old CA1 |
| `CA1_age_conserved` | `|rho_CA1_Young| > 0.3` AND `|rho_CA1_Old| > 0.3` AND same sign | Candidate consistently associated across CA1 age groups |
| `CA1_age_shifted` | `|rho_CA1_Young| > 0.3` AND `|rho_CA1_Old| > 0.3` AND opposite sign | Candidate coupling with opposite direction between Young and Old CA1 |
| `CA3_old_specific` | `|rho_CA3_Old| > 0.3` AND `|rho_CA3_Young| < 0.3` | Candidate associated in Old CA3 but NOT Young CA3 |
| `CA3_young_specific` | `|rho_CA3_Young| > 0.3` AND `|rho_CA3_Old| < 0.3` | Candidate associated in Young CA3 but NOT Old CA3 |
| `CA3_age_conserved` | `|rho_CA3_Young| > 0.3` AND `|rho_CA3_Old| > 0.3` AND same sign | Candidate consistently associated across CA3 age groups |
| `CA3_age_shifted` | `|rho_CA3_Young| > 0.3` AND `|rho_CA3_Old| > 0.3` AND opposite sign | Candidate coupling with opposite direction between Young and Old CA3 |

**CA3 Young (n=3) caveat**: All `CA3_young_specific`, `CA3_age_conserved`, and `CA3_age_shifted` classes involving CA3 Young rho are EXPLORATORY ONLY. The label `severely_underpowered_exploratory` MUST be attached to these classes in all output tables.

### 7.2 Cross-Region Classes (Within Age)

| Class | Criteria | Interpretation |
|-------|----------|---------------|
| `Young_CA1_specific` | `|rho_CA1_Young| > 0.3` AND `|rho_CA3_Young| < 0.3` | Young candidate specific to CA1 |
| `Young_CA3_specific` | `|rho_CA3_Young| > 0.3` AND `|rho_CA1_Young| < 0.3` | Young candidate specific to CA3 (EXPLORATORY, n=3) |
| `Young_region_conserved` | `|rho_CA1_Young| > 0.3` AND `|rho_CA3_Young| > 0.3` AND same sign | Young candidate conserved across CA1/CA3 (EXPLORATORY, n=3-4) |
| `Old_CA1_specific` | `|rho_CA1_Old| > 0.3` AND `|rho_CA3_Old| < 0.3` | Old candidate specific to CA1 |
| `Old_CA3_specific` | `|rho_CA3_Old| > 0.3` AND `|rho_CA1_Old| < 0.3` | Old candidate specific to CA3 |
| `Old_region_conserved` | `|rho_CA1_Old| > 0.3` AND `|rho_CA3_Old| > 0.3` AND same sign | Old candidate conserved across CA1/CA3 |

### 7.3 Special Classes

| Class | Criteria | Interpretation |
|-------|----------|---------------|
| `old_enhanced_coupling` | `|rho_Old| > |rho_Young| + 0.2` in same region | Candidate coupling is stronger in Old than Young |
| `young_enhanced_coupling` | `|rho_Young| > |rho_Old| + 0.2` in same region | Candidate coupling is stronger in Young than Old (CA3 EXPLORATORY) |
| `low_confidence` | `|rho| < 0.3` in ALL 4 strata | Below screening threshold in all strata; reported but not highlighted |
| `insufficient_n` | n < 3 in any required stratum | Cannot compute class (e.g., CA3 Young candidate classes when G1_CA3 missing — checked per §13.0) |

### 7.4 Priority Order for Class Assignment

When a gene × module pair qualifies for multiple classes, assign the MOST SPECIFIC class:

1. Age-specific classes (old_specific, young_specific) > age_conserved > low_confidence
2. Within-age region-specific classes > region-conserved
3. Age-shifted classes > age_conserved (age-shifted is more informative)
4. Special classes (enhanced coupling) are annotations, not primary classifications

Add a `secondary_class` column for any additional classes the pair also qualifies for.

### 7.5 Expected Class Sizes (Rough Estimates)

These are rough estimates only. Actual counts depend on data:

| Class | Expected Range | Confidence Rationale |
|-------|---------------|---------------------|
| CA1_old_specific | 50-150 | CA1 n=8 Old vs n=4 Young |
| CA1_young_specific | 20-80 | CA1 Young n=4 (small-n caution) |
| CA1_age_conserved | 100-300 | Both strata have detectable rho |
| CA3_old_specific | 50-150 | CA3 n=8 Old |
| CA3_young_specific | 10-50 | CA3 Young n=3 (severely underpowered) |
| CA3_age_conserved | 30-100 | CA3 Young n=3 limits detection |
| Old_CA1_specific | 50-150 | CA1 Old vs CA3 Old |
| Old_CA3_specific | 50-150 | CA3 Old vs CA1 Old |
| Old_region_conserved | 100-300 | Both Old regions |
| low_confidence | 3000-6000 | Below screening in all strata |

---

## 8. Reuse vs Re-computation Rules

### 8.1 Reuse from Phase 13 (Do NOT Re-compute)

| Item | How to Reuse |
|------|--------------|
| `universe_regulator` gene list (1,198 genes) | Filter by gene presence in pseudobulk. Do NOT re-derive TF annotation. |
| TF annotation table | Read `phase13_tf_annotation_coverage.csv`. Do NOT re-query DoRothEA/CollecTRI/GO. |
| Target gene annotation | Read `docs/target_genes.csv`. Do NOT re-audit gene presence. |
| Module gene-set members | Read `module_gene_set_final_audit.csv`. Do NOT re-define modules. |
| DE results (Phase 11/12) | Read existing CSVs. No new DESeq2 run. |
| Module scores | Read `module_scores_all_regions.csv`. Do NOT re-compute module scores. |

### 8.2 New Computation

| Item | Computation |
|------|-------------|
| Age-stratified Spearman ρ | NEW: 4 strata × ~1,198 genes × ~11 modules ≈ 52,712 tests |
| Age-stratified p-value / q-value | NEW: computed from `cor.test()` and BH correction |
| Candidate classes | NEW: based on age-stratified rho patterns |
| Delta scores | NEW: rho_Old - rho_Young per region |
| Region-shifted within age | NEW: rho_CA1 - rho_CA3 per age stratum |

### 8.3 What Phase 13c Does NOT Do

- Does NOT re-run Phase 11 or 12 DESeq2
- Does NOT re-compute module scores
- Does NOT re-derive the regulator universe
- Does NOT re-run TF annotation (DoRothEA/CollecTRI/GO)
- Does NOT load the Hippo Seurat RDS object (unless spatial visualization is needed; see §11.7)
- Does NOT perform spot-level inference
- Does NOT modify Phase 13 or Phase 13b output files

---

## 9. Planned Output Files

**CRITICAL**: All outputs merge into EXISTING Phase 13 directories. Do NOT create `phase13c_*` new directories. File names use `phase13c_` prefix.

### 9.1 Data Outputs

All under `data/processed/spatial/phase13_candidate_regulator/`:

| # | File | Description | Key Columns |
|---|------|-------------|-------------|
| O1 | `phase13c_age_stratified_rho.csv` | Full age-stratified rho table for `universe_regulator` | gene, module, rho_CA1_Young, pvalue_exploratory_CA1_Young, qvalue_exploratory_CA1_Young, rho_CA1_Old, pvalue_exploratory_CA1_Old, qvalue_exploratory_CA1_Old, rho_CA3_Young, pvalue_exploratory_CA3_Young, qvalue_exploratory_CA3_Young, rho_CA3_Old, pvalue_exploratory_CA3_Old, qvalue_exploratory_CA3_Old, n_CA1_Young, n_CA1_Old, n_CA3_Young, n_CA3_Old, self_correlation_flag |
| O2 | `phase13c_age_stratified_scores.csv` | Age-stratified candidate scores and delta scores | O1 columns + candidate_score_CA1_Young, candidate_score_CA1_Old, candidate_score_CA3_Young, candidate_score_CA3_Old, candidate_score_with_DE_support_sensitivity_CA1_Young (optional), candidate_score_with_DE_support_sensitivity_CA1_Old (optional), candidate_score_with_DE_support_sensitivity_CA3_Young (optional), candidate_score_with_DE_support_sensitivity_CA3_Old (optional), DE_age_support_flag_CA1, DE_age_support_flag_CA3, delta_rho_CA1, delta_rho_CA3, delta_region_Young, delta_region_Old, DE_age_lfc_CA1, DE_age_padj_CA1, DE_age_lfc_CA3, DE_age_padj_CA3, DE_region_lfc, sample_size_label |
| O3 | `phase13c_age_specific_candidate_classes.csv` | Primary age-specific candidate class per gene × module | O2 columns + candidate_class, secondary_class, old_enhanced_coupling_flag, young_enhanced_coupling_flag |
| O4 | `phase13c_young_old_delta_scores.csv` | Genes ranked by delta_rho (Old - Young) per region | gene, module, region, rho_Young, rho_Old, delta_rho, sample_size_label |
| O5 | `phase13c_old_specific_candidates.csv` | Old-specific candidates across regions (union CA1 + CA3) | gene, module, region, candidate_class (CA1_old_specific / CA3_old_specific), rho_Old, rho_Young, delta_rho, DE_age_lfc, DE_age_padj |
| O6 | `phase13c_young_specific_candidates.csv` | Young-specific candidates (CA3 flagged EXPLORATORY) | Same columns as O5; CA3 Young n=3 caveat in header comment |
| O7 | `phase13c_age_shifted_candidates.csv` | Age-shifted candidates (rho sign flips) | gene, module, region, rho_Young, rho_Old, delta_rho, sign_change, sample_size_label |
| O8 | `phase13c_region_shifted_by_age_candidates.csv` | Region-shifted within-age candidates | gene, module, age_stratum, rho_CA1, rho_CA3, delta_region, sample_size_label |
| O9 | `phase13c_top_candidates_by_region_age.csv` | Top 20 per region×age by |rho| | gene, module, region, age, rho, candidate_class, self_correlation_flag |
| O10 | `phase13c_validation_summary.csv` | Phase 13c validation checks (separate from Phase 13) | check, status, details |
| O11 | `phase13c_provenance.csv` | Machine-readable provenance record | input_path, output_path, package_version, parameter, value |

**CSV header conventions**:
- All CSVs include a comment line (if write.csv supports it) or a separate README noting: "CAUTION: Rows are gene × module pairs, not unique genes."
- rho columns retain full precision (decimal values from `cor.test()`)
- p-value columns named `pvalue_exploratory_<stratum>` with the `_exploratory` suffix
- q-value columns named `qvalue_exploratory_<stratum>` with the `_exploratory` suffix
- All CA3 Young-related columns include `severely_underpowered_exploratory` in sample_size_label

### 9.2 Figure Outputs

All under `figures/spatial/phase13_candidate_regulator/`:

**Required figures**:

| # | File | Type | Description | Priority |
|---|------|------|-------------|----------|
| F01 | `phase13c_heatmap_CA1_young.png` | Heatmap | Age-stratified rho for top 30 candidates (CA1 Young), genes × modules. n=4 displayed. | HIGH |
| F02 | `phase13c_heatmap_CA1_old.png` | Heatmap | Age-stratified rho for top 30 candidates (CA1 Old), genes × modules. n=8 displayed. | HIGH |
| F03 | `phase13c_heatmap_CA3_young.png` | Heatmap | Age-stratified rho for top 30 candidates (CA3 Young), genes × modules. n=3 displayed. **Title MUST include "CAUTION: n=3, severely underpowered, exploratory only".** | HIGH |
| F04 | `phase13c_heatmap_CA3_old.png` | Heatmap | Age-stratified rho for top 30 candidates (CA3 Old), genes × modules. n=8 displayed. | HIGH |
| F05 | `phase13c_delta_CA1_old_minus_young.png` | Scatter or bar | Delta rho (Old - Young) for CA1, top 20 by `|delta|`. Annotated with Old-specific / Young-specific class. | HIGH |
| F06 | `phase13c_delta_CA3_old_minus_young.png` | Scatter or bar | Delta rho for CA3, top 20 by `|delta|`. CAVEAT: CA3 Young n=3. Annotated. | MEDIUM |
| F07 | `phase13c_region_age_class_barplot.png` | Barplot | Candidate class counts by candidate_class, colored by region×age. n labels. | HIGH |
| F08 | `phase13c_old_specific_top_dotplot.png` | Dotplot | Top 20 Old-specific candidates (across CA1/CA3), showing rho_Old, rho_Young, DE evidence. | HIGH |
| F09 | `phase13c_old_enhanced_coupling_dotplot.png` | Dotplot | Top 20 old-enhanced-coupling candidates (|rho_Old| > |rho_Young| + 0.2). | MEDIUM |

**Optional figures** (memory-safe, skip if dependencies unavailable, may skip if data insufficient — must WARN with reason):

| # | File | Type | Description |
|---|------|------|-------------|
| F10 | `phase13c_young_specific_top_dotplot.png` | Dotplot | Top 20 Young-specific candidates with small-n caveat. Skip if < 5 candidates qualify. |
| F11 | `phase13c_CA1_Old_spatial_<gene>.png` | SpatialFeaturePlot | Top 3 Old-specific CA1 candidates (if Hippo RDS loaded). Skip on memory failure. |
| F12 | `phase13c_CA3_Old_spatial_<gene>.png` | SpatialFeaturePlot | Top 3 Old-specific CA3 candidates (if Hippo RDS loaded). Skip on memory failure. |

### 9.3 Documentation Outputs

| # | File | Action |
|---|------|--------|
| D1 | `docs/spatial_phase13c_age_stratified_regulator_plan.md` | **THIS FILE**. Authoritative Phase 13c plan. |
| D2 | `docs/spatial_phase13_candidate_regulator_analysis_guide.md` | **APPEND** a "Phase 13c: Age-Stratified Candidate Regulator Discovery" section explaining how to read age-stratified results. |

---

## 10. Figure Conventions

### 10.1 General Rules

- **PNG is authoritative**. PDF optional/best-effort with `cairo_pdf()` + `tryCatch` fallback.
- All correlation figures: display n prominently in title or subtitle.
- Color scheme per Phase 11/12/13 convention:
  - CA1 = `#56B4E9` (blue)
  - CA3 = `#E69F00` (orange)
  - Young = green
  - Old = purple
- Use `ggplot2` + `pheatmap` (no new package dependencies).
- Call `gc()` after each large figure.
- Explicit file paths for all device calls; NO default `Rplots.pdf`.

### 10.2 PDF Safety Protocol (Identical to Phase 13)

```
safe_pdf <- function(filename, width, height, plot_expr) {
  tryCatch({
    cairo_pdf(file.path(fig_dir, filename), width = width, height = height)
    eval(plot_expr)
    dev.off()
  }, error = function(e) {
    message("PDF failed for ", filename, ": ", e$message)
    while (dev.cur() > 1) dev.off()
  })
}

# ... after each figure block ...
while (dev.cur() > 1) invisible(dev.off())
```

### 10.3 CA3 Young n=3 Label in Figures

Every figure that includes CA3 Young data MUST display:
- Title: `(CAUTION: CA3 Young n=3 — severely underpowered)`
- OR subtitle: `n_CA3_Young = 3 (exploratory only)`

---

## 11. Validation Checks

### 11.1 Pre-Execution (Within Script, E0)

| # | Check | Expected | Severity |
|---|-------|----------|----------|
| V01 | Phase 13 validation summary has 0 FAIL | All PASS | FATAL |
| V02 | Phase 13b validation summary has 0 FAIL | All PASS | FATAL |
| V03 | All required inputs (I1-I12) exist | File paths readable | FATAL |
| V04 | Sample manifest matches expected counts (CA1 Young=4, CA1 Old=8, CA3 Young=3, CA3 Old=8) | Exact match | FATAL |
| V05 | `universe_regulator` gene count ≥ 100 | > 100 genes | WARN |
| V06 | Module gene-set audit file readable and has ≥ 8 modules | Readable, ≥ 8 rows | FATAL |
| V07 | Module scores CSV has rows for all CA1/CA3 GEMgroups | All 23 CA1/CA3 Young+Old GEMgroups present | FATAL |
| V08 | No causal language in plan or script | grep check | FATAL |
| V09 | renv.lock MD5 recorded before script execution | Recorded and logged | INFO |
| V10 | No `phase13c_*` new directories created | Output goes to `phase13_candidate_regulator/` | FATAL |

### 11.2 Post-Execution (Within Script, E6)

| # | Check | Expected | Severity |
|---|-------|----------|----------|
| V11 | O1-O11 exist and have > 0 rows | All required CSVs present, non-empty | FATAL |
| V12 | No q-value used for `candidate_class` column | `candidate_class` assignment logic uses |rho| thresholds only | FATAL |
| V13 | CA3 Young candidate classes contain `severely_underpowered_exploratory` in `sample_size_label` | grep confirms | FATAL |
| V14 | All p-value columns named `pvalue_exploratory_*` | grep confirms | WARN |
| V15 | All q-value columns named `qvalue_exploratory_*` | grep confirms | WARN |
| V16 | `self_correlation_flag` column has 3 possible values only | `regulator_candidate`, `module_member_or_effector`, `target_gene_candidate` | FATAL |
| V17 | Each required figure F01-F09 exists individually | F01–F09 each present and non-zero size. If any required figure cannot be generated due to insufficient data (e.g., < 5 candidates for a dotplot), log WARN with specific reason. CA3 Young n=3 is NOT a valid skip reason — F03 must be generated. | FATAL (per-figure); WARN if skipped with documented reason |
| V18 | No `phase13c_*` output directory created | Directory check | FATAL |
| V19 | All output files have `phase13c_` prefix | Filename check | FATAL |
| V20 | `Rplots.pdf` absent from repo root | file.exists check | FATAL |
| V21 | renv.lock unchanged from pre-execution MD5 (unless approved package install) | MD5 match | WARN |
| V22 | No full Seurat object saved | No new RDS files in output directories | WARN |
| V23 | No causal language in any output | `grep -i 'causal\|upstream regulator\|drives\|mediates'` | FATAL |
| V24 | Analysis guide updated with Phase 13c section | grep for "Phase 13c" or "Age-Stratified" in guide | WARN |
| V25 | No Phase 13 or 13b files overwritten or deleted | File existence pre/post check | FATAL |

---

## 12. Stop Conditions

### 12.1 Fatal Stop Conditions

| # | Condition | Action |
|---|----------|--------|
| FATAL-01 | Phase 13 core outputs missing (`phase13_regulator_score_default.csv`, `phase13_validation_summary.csv`) | STOP — Re-run Phase 13 |
| FATAL-02 | Phase 13 validation has any FAIL | STOP — Fix Phase 13 issues first |
| FATAL-03 | Module score file missing or unreadable | STOP — Verify Phase 12 output |
| FATAL-04 | Pseudobulk count RDS missing or corrupted | STOP — Regenerate from Phase 11 |
| FATAL-05 | CA3 Young has n < 3 (after manifest check) | STOP or mark CA3 Young stratum UNAVAILABLE; produce CA1-only age-stratified analysis |
| FATAL-06 | CA1 Young has n < 3 | STOP — Cannot compute meaningful Young CA1 correlation |
| FATAL-07 | Any attempt to treat spots as replicates | STOP — Redesign. This is a hard constraint. |
| FATAL-08 | Any output labeled as "causal regulator" or "upstream regulator" | STOP — Reword to "candidate regulator-associated gene" |
| FATAL-09 | Any attempt to use q-value or p-value as primary ranking criterion or candidate class assignment | STOP — Redesign using |rho| thresholds |
| FATAL-10 | `candidate_class` assignment logic uses q-value or p-value | STOP — Redesign |
| FATAL-11 | New `phase13c_*` directories created for data or figures | STOP — Output must go to existing `phase13_candidate_regulator/` paths |
| FATAL-12 | Any Phase 13 or 13b output file modified or deleted | STOP — Phase 13c is additive only |
| FATAL-13 | Module gene-set members NOT excluded from ranking for own module | STOP — Implement Strategy A per Phase 13 §5.3 |

### 12.2 Non-Fatal Warnings

| # | Condition | Action |
|---|----------|--------|
| WARN-01 | CA3 Young n=3 — all CA3 Young results carry `severely_underpowered_exploratory` label | PROCEED with caveat in all CA3 Young outputs |
| WARN-02 | CA1 Young n=4 — all CA1 Young results carry `exploratory_small_n` label | PROCEED with caveat |
| WARN-03 | Fewer than 10 genes with |rho| > 0.5 in any stratum | PROCEED; note limited strong associations |
| WARN-04 | Fewer than 20 genes classified as Old-specific across all modules | PROCEED; note limited age-specific signal |
| WARN-05 | < 50% overlap between Phase 13c Old-specific candidates and Phase 13 all-age ranking top 100 | PROCEED; note that age stratification reveals different candidates than all-age analysis — this is expected, not a failure |
| WARN-06 | Memory failure during optional spatial map generation (Hippo RDS loading) | SKIP optional spatial maps; proceed with core table outputs |
| WARN-07 | PDF device fails | PNG is authoritative; PDF is optional |

---

## 13. Execution Script Plan (Future)

### 13.1 Script Outline

Future execution script: **`R/spatial/s13c_age_stratified_regulator_discovery.R`**

Structure (for reference only — do NOT create now):

```
E0: Setup & Validation
  - Record renv.lock MD5
  - Load libraries (dplyr, ggplot2, pheatmap, patchwork, stringr, digest)
  - Verify all inputs (I1-I12)
  - Run pre-execution checks V01-V10
  - If FATAL, stop

E1: Load Data
  - Read pseudobulk counts, subset CA1/CA3 columns
  - Compute log2(CPM+1)
  - Read module scores, subset by age strata
  - Read module audit, build module_genes list
  - Read Phase 13 regulator universe, TF annotation
  - Read Phase 11/12 DE results
  - Read target genes (annotation only)

E2: Compute Age-Stratified Spearman Rho
  - For each gene × module pair in universe_regulator:
    - cor.test() for CA1 Young (4 samples)
    - cor.test() for CA1 Old (8 samples)
    - cor.test() for CA3 Young (3 samples)
    - cor.test() for CA3 Old (8 samples)
  - Apply self-correlation exclusion per module
  - Apply BH correction per stratum
  - Write phase13c_age_stratified_rho.csv

E3: Compute Candidate Scores and Delta Scores
  - Compute primary candidate_score per stratum = |rho|
  - Add DE_age_support_flag annotation column (padj < 0.05 in Phase 12 region match)
  - Optionally compute candidate_score_with_DE_support_sensitivity (= |rho| + 0.1 if DE support)
  - Compute delta_rho_CA1, delta_rho_CA3
  - Compute delta_region_Young, delta_region_Old
  - Write phase13c_age_stratified_scores.csv

E4: Assign Candidate Classes
  - Classify per gene × module per §7.1-7.3
  - Prioritize most specific class (§7.4)
  - Write phase13c_age_specific_candidate_classes.csv
  - Write specialized CSV outputs (O4-O9)

E5: Generate Figures
  - F01-F09 per §9.2 (9 required figures)
  - F10-F12 optional (skip if data insufficient or memory fails)
  - PNG authoritative; PDF optional
  - All n labels, CA3 Young caveats
  - gc() after each figure group

E6: Validation & Cleanup
  - Run post-execution checks V11-V25
  - Write phase13c_validation_summary.csv
  - Write phase13c_provenance.csv
  - Update analysis guide with Phase 13c section
  - Rplots.pdf cleanup
  - Report pass/warn/fail counts
```

### 13.2 Package Policy

**NO new packages required.** Phase 13c reuses packages already in the Phase 13 environment:

| Package | Usage |
|---------|-------|
| dplyr | Data manipulation, filtering, sorting |
| ggplot2 | All figures |
| pheatmap | Heatmaps |
| patchwork | Multi-panel composition |
| stringr | String manipulation |
| digest | renv.lock MD5 check |
| Seurat (optional) | SpatialFeaturePlot ONLY if generating optional F11/F12 |

**renv.lock**: Must remain unchanged throughout Phase 13c unless an approved package installation is needed (unexpected). Current MD5 baseline: `b2bf05038b440576855d043e137c5883`.

### 13.3 Memory Management

- Phase 13c reads CSVs primarily (kilobytes to megabytes)
- Pseudobulk count RDS (~2.4 MB) — loaded once in E1
- Module scores CSV (~20 KB)
- NO large Seurat object load unless optional spatial maps (F11/F12) are generated
- If spatial maps: load Hippo RDS → plot → unload → gc(). Exactly one load/unload cycle.
- `gc()` after each major step

---

## 14. Guide Update Plan (Phase 13c Section)

After Phase 13c execution, append to `docs/spatial_phase13_candidate_regulator_analysis_guide.md`:

### 14.1 Section Content

The Phase 13c section must explain:

1. **What Phase 13c adds over Phase 13/13b**:
   - "Phase 13/13b ranked candidate regulator-associated genes using all-age Spearman rho across CA1 and CA3. Phase 13c extends this with age-stratified analysis, computing separate rho values for Young and Old strata within each region."

2. **How to read age-stratified tables**:
   - "Each row is a gene × module pair with separate rho values for Young and Old. Compare rho_Young vs rho_Old within the same gene×module to identify age-specific patterns."

3. **How to interpret rho direction**:
   - "Positive rho: gene expression and module score increase together. Negative rho: gene expression increases while module score decreases. In Old strata, this may reflect aging-related coupling changes."

4. **Why small n matters**:
   - "CA1 Young n=4 and CA3 Young n=3 are underpowered for Spearman correlation. rho estimates in these strata are noisy. Compare with Old strata (n=8) for more stable estimates."

5. **Why CA3 Young is weak**:
   - "CA3 has only 3 Young GEMgroups (G2, G3, G4). All CA3 Young results are labeled 'severely_underpowered_exploratory'. Do NOT interpret these as definitive findings."

6. **Difference between region-specific and age-specific candidates**:
   - "Age-specific candidates (e.g., CA1_old_specific) show strong rho in one age but not the other WITHIN the same region. Region-specific candidates (e.g., Old_CA1_specific) show strong rho in one region but not the other WITHIN the same age stratum."

7. **Why these are candidate regulator-associated genes, not causal regulators**:
   - "Age-stratified rho is still correlation. Even when rho is strong and age-specific, we cannot determine whether the gene regulates mitochondrial programs, the programs regulate the gene, or a third factor drives both. Age specificity adds biological plausibility but does NOT establish causality."

8. **Key files to look at**:
   - `phase13c_age_stratified_rho.csv` — full rho table
   - `phase13c_old_specific_candidates.csv` — Old-specific candidates
   - `phase13c_age_shifted_candidates.csv` — candidates with rho sign reversal between Young and Old
   - `phase13c_region_shifted_by_age_candidates.csv` — region-specific candidates within same age

---

## 15. References and Skills Reviewed

### 15.1 Authoritative Documentation

| # | Source Title | URL / Path | Reviewed Status | Section Reviewed | Design Impact |
|---|-------------|-----------|-----------------|-----------------|---------------|
| S1 | Phase 13 Candidate Regulator Discovery Plan v1.2 | `docs/spatial_phase13_candidate_regulator_discovery_plan.md` | reviewed | Full document (984 lines): §5.3 self-correlation safeguards, §7 multiple testing constraints, §8 ranking formula, §9 candidate classes, §13 stop conditions | All self-correlation rules, ranking philosophy, and stop conditions inherited directly; age stratification added as new layer |
| S2 | Phase 13b CA3 Symmetric Display Plan v1.0 | `docs/spatial_phase13b_ca3_symmetric_display_plan.md` | reviewed | Full document (449 lines): §4.2 figure specs, §5 validation, §7 execution plan | CA3 mirror figure conventions adopted; `phase13c_` prefix naming convention to avoid collision |
| S3 | DESeq2 Official Vignette | https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html | not_reviewed (already reviewed in Phase 13 plan) | Standard workflow, design formulas, LFC shrinkage | Phase 13c does no new DESeq2 — reuses existing Phase 11/12 results. No re-review needed. |
| S4 | Seurat AddModuleScore Reference | https://satijalab.org/seurat/reference/addmodulescore | not_reviewed (already reviewed in Phase 13 plan) | Module score methodology | Module scores are reused readouts from Phase 11/12; no re-computation |

### 15.2 ref-bio Routing Files

| # | File | Reviewed Status | Key Entries |
|---|------|----------------|-------------|
| R1 | `.opencode/skills/ref-bio/SKILL.md` | not_reviewed (already reviewed in Phase 13 plan) | Routing-only policy, no summarization |
| R2 | `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml` | not_reviewed (already reviewed in Phase 13 plan) | Same entries as Phase 13 |

### 15.3 bio-* Skills

| # | Skill Path | Reviewed Status | Key Guidance | Design Impact |
|---|-----------|----------------|-------------|---------------|
| B1 | `bio-gene-regulatory-networks-coexpression-networks/SKILL.md` | inspected | §"<15 samples: none reliable — correlation estimates too noisy; report this, do not force WGCNA". WGCNA requires ≥15-20 samples. Marginal (WGCNA/Pearson) vs partial (GGM) correlation distinction. | Confirms WGCNA is NOT appropriate for n=3-8 per stratum. Simple Spearman rho is the correct approach. Phase 13c explicitly does NOT build co-expression networks. |
| B2 | `bio-differential-expression-de-results/SKILL.md` | inspected | padj=NA handling, BH FDR, IHW, TREAT vs post-hoc filtering, p-value histogram diagnostics | Phase 13c only reads DE results as supporting annotation; no new DE analysis. padj threshold from Phase 11/12 is reused as-is. |
| B3 | `bio-differential-expression-deseq2-basics/SKILL.md` | not_reviewed (already reviewed in Phase 13 plan) | Negative-binomial GLM, Wald test, design formulas | Phase 13c performs NO new DESeq2. Reuses Phase 11/12 results. |
| B4 | `bio-data-visualization-heatmaps-clustering/SKILL.md` | not_reviewed (already reviewed in Phase 13 plan) | ComplexHeatmap, pheatmap, ward.D2, robust bounds | Heatmap conventions for F01-F04. Reuse existing Phase 13 conventions. |
| B5 | `bio-data-visualization-distribution-plots/SKILL.md` | not_reviewed (already reviewed in Phase 13 plan) | Weissgerber 2015 sample-size honesty, N-visible encoding | All figures display n prominently. Consistent with Phase 13. |

### 15.4 Project-Specific Documents

| # | File | Reviewed Status | Contribution |
|---|------|----------------|-------------|
| P1 | `README.md` | reviewed | Confirmed spatial branch conventions, output paths, memory rules |
| P2 | `AGENTS.md` | reviewed | Confirmed Phase 12 guardrails (replicate unit, CA3 Young n=3, Complex V caveat, no spot-level inference), script naming convention, ref-bio/bio-* rules |
| P3 | `docs/repository_structure.md` | reviewed | Confirmed output directory policy, file naming conventions, cache/scratch path rules |

---

## 16. Execution Plan Compliance Report

### 16.1 Plan-Only Verification

| Check | Status |
|-------|--------|
| Only created/modified `docs/spatial_phase13c_age_stratified_regulator_plan.md` | YES |
| No R scripts created | YES |
| No R analysis executed | YES |
| No RDS loaded | YES |
| No packages installed | YES |
| No renv.lock modified | YES |
| No data/figures generated | YES |
| No `phase13c_*` directories created | YES |

### 16.2 Files Inspected

| # | File | Status |
|---|------|--------|
| 1 | `README.md` | Reviewed — spatial conventions, output paths |
| 2 | `AGENTS.md` | Reviewed — guardrails, Phase 12 rules, ref-bio/bio-* rules |
| 3 | `docs/repository_structure.md` | Reviewed — output directory policy |
| 4 | `docs/spatial_phase13_candidate_regulator_discovery_plan.md` | Reviewed (984 lines) — inherited self-correlation, ranking philosophy, stop conditions |
| 5 | `docs/spatial_phase13b_ca3_symmetric_display_plan.md` | Reviewed (449 lines) — CA3 mirror conventions, naming patterns |
| 6 | `docs/spatial_phase13_candidate_regulator_analysis_guide.md` | Reviewed (155 lines) — guide structure, update location identified |
| 7 | `R/spatial/s13_candidate_regulator_discovery.R` | Reviewed (1470 lines) — confirmed no age-stratified rho computation, confirmed F07/F08/F07b/F08b are display-only |
| 8 | `R/spatial/s13b_ca3_symmetric_display.R` | Reviewed (713 lines) — confirmed display-supplement only |
| 9 | `data/processed/spatial/phase13_candidate_regulator/phase13_regulator_score_default.csv` | Reviewed header — 13,179 rows, 12 columns, includes `score_CA3_default` and `rank_CA3` |
| 10 | `data/processed/spatial/phase13_candidate_regulator/phase13_spearman_rho_universe_regulator.csv` | Reviewed header — 13,179 rows, 26 columns, includes `qval_CA1`, `qval_CA3` but NOT age-stratified |
| 11 | `data/processed/spatial/phase13_candidate_regulator/phase13_candidate_classes.csv` | Reviewed — 6 classes as expected; NO age_stratified classes |
| 12 | `data/processed/spatial/phase13_candidate_regulator/phase13_validation_summary.csv` | Reviewed — 12 PASS / 0 FAIL / 0 WARN |
| 13 | `data/processed/spatial/phase13_candidate_regulator/phase13b_display_validation_summary.csv` | Reviewed — 17 PASS / 0 WARN / 0 FAIL |
| 14 | `data/processed/spatial/phase12_young_aged_region_comparison/module_scores_all_regions.csv` | Reviewed (44 rows) — confirmed GEMgroup-level mean module scores with Age/Region columns |
| 15 | `data/processed/spatial/phase13_candidate_regulator/pseudobulk_sample_manifest.csv` | Reviewed (32 rows) — confirmed CA1 Young=4 (G1-G4), CA1 Old=8 (G9-G16), CA3 Young=3 (G2-G4), CA3 Old=8 (G9-G16), Middle=8 |

### 16.3 Key Design Decisions

| Decision | Rationale |
|----------|----------|
| Effect-size-first ranking (|rho| thresholds) | n=3-8 is too small for reliable p/q-values. Consistent with Phase 13 philosophy (§7.4). |
| pvalue_exploratory / qvalue_exploratory naming | Explicitly signals limited validity. Prevents misinterpretation as formal statistical significance. |
| No q-value in candidate_class logic | Enforced by FATAL-09, FATAL-10. Candidate classes use |rho| thresholds only. |
| `phase13c_` prefix on all new outputs | Prevents collision with Phase 13/13b outputs. Files merge into existing directories. |
| No new directories for data/figures | Outputs go to `phase13_candidate_regulator/` data/figs dirs. Consistent with Phase 13b merge-back rule. |
| Reuse Phase 13 regulator universe and TF annotation | No re-derivation. Phase 13 already has curated universe_regulator (1,198 genes). |
| No new DESeq2 run | Reuse Phase 11/12 DE results as supporting annotation only. |
| CA3 Young n=3 → `severely_underpowered_exploratory` | Consistent with Phase 12 guardrails ("Must be labeled `primary_with_caution_young_n3`"). Then `severely_underpowered_exploratory` is more appropriate since no DE is run. |
| `module_member_or_effector` exclusion from ranking | Identical to Phase 13 Strategy A. Self-correlation cannot be interpreted as upstream evidence. |
| Target genes NOT in regulator_score | Consistent with Phase 13 §5.3, §8.1. Target gene membership is annotation only. |
| Spatial maps optional, not blocking | Hippo RDS is 1.3GB; core tabular outputs do not depend on it. |

### 16.4 Safety Assessment

Phase 13c is **safe to execute** after user review because:

1. **Read-only on Phase 13/13b core data**: All inputs are CSVs/RDS files; no existing file is modified or deleted.
2. **No new packages**: Uses only packages already in the Phase 13 environment.
3. **No renv.lock change expected**: No `renv::install()` or package additions.
4. **Additive outputs**: All new files use `phase13c_` prefix; no collisions with existing files.
5. **Fatal stop conditions**: Script stops before output if critical inputs are missing or constraints violated.
6. **Memory guard**: Hippo RDS loading is optional, guarded by tryCatch with WARN fallback.
7. **Rplots.pdf guard**: Explicit file paths for all devices; cleanup at end.
8. **No spot-level inference**: All analysis at GEMgroup pseudobulk level.
9. **No causal language**: All outputs use "candidate regulator-associated gene" terminology.

---

## 17. Summary

Phase 13c fills the age-stratified gap in Phase 13/13b by:

1. Computing per-age Spearman rho (4 strata: CA1 Young, CA1 Old, CA3 Young, CA3 Old)
2. Assigning age-specific candidate classes (Old-specific, Young-specific, age-conserved, age-shifted, region-shifted within age)
3. Producing 11 data CSVs and 9-12 figures, all with `phase13c_` prefix in existing directories
4. Updating the analysis guide with an age-stratified interpretation section
5. Maintaining all Phase 13 safeguards: self-correlation exclusion, effect-size-first ranking, exploratory q-values, no causal language, no spot-level inference

The plan is self-contained, does not modify Phase 13/13b, and is ready for execution after review.

---

**Final Report for This Plan-Only Task**:

- **Files inspected**: 15 (listed in §16.2)
- **File created**: `docs/spatial_phase13c_age_stratified_regulator_plan.md`
- **ref-bio inspected**: Not reviewed (already done in Phase 13 plan; Phase 13c reuses same methodology with age stratification)
- **bio-* skills inspected**: 2 reviewed for Phase 13c (`bio-gene-regulatory-networks-coexpression-networks`, `bio-differential-expression-de-results`); 3 carried forward from Phase 13 without re-review (not needed for this plan)
- **Confirmation that no R was run**: YES — plan only
- **Confirmation that no data/figures were generated**: YES — plan only
- **Confirmation that renv.lock was unchanged**: YES — no package operations
- **Phase 13c is safe to execute after review**: YES — all safety conditions met (§16.4)
