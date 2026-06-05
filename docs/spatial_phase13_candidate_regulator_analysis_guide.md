# Phase Spatial-13: Candidate Regulator Analysis Guide

Generated: 2026-06-05 13:30:37
Plan version: v1.2

## Overview

This guide describes the Phase Spatial-13 candidate regulator discovery results.
All outputs identify **candidate regulator-associated genes** — genes whose expression
correlates with mitochondrial module scores at the GEMgroup pseudobulk level.

**This is NOT a causal regulator analysis.** Correlation does not imply regulation.
No perturbation experiments, mechanistic validation, or directionality tests were performed.

## What Each Analysis Asks

1. **Spearman correlation**: Does gene g's expression covary with module m's score across GEMgroups in CA1 and/or CA3?
2. **DE-informed ranking**: Do genes that change with aging (Old vs Young) or differ between regions (CA3 vs CA1) rank higher?
3. **TF annotation**: Is the gene annotated as a transcription factor or regulator in DoRothEA, CollecTRI, or GO? (Eligibility only — NOT proof of regulation in this dataset.)

## Key Outputs

### Primary Rankings
- `phase13_regulator_score_default.csv`: Primary ranking by default config (ρ 50% + DE_age 30% + DE_region 20%). Includes rank_CA1 and rank_CA3 (module members have NA rank for own module).
- `phase13_candidate_classes.csv`: Classification of candidate regulator-associated genes.
- `phase13_spearman_rho_all.csv`: All Spearman ρ values with BH-corrected q-values (exploratory).

### Sensitivity
- `phase13_regulator_score_sensitivity.csv`: All 3 weight configs (default, DE-emphasis, ρ-only).
- Compare configs to assess robustness of rankings.

### Annotation
- `phase13_tf_annotation_coverage.csv`: TF sources (DoRothEA A-C, CollecTRI, GO:0003700/GO:0140110).
- `phase13_target_gene_module_association.csv`: Target gene module correlations (separate from regulator ranking).
- `phase13_self_correlation_flags.csv`: Self-correlation flags (module member exclusion).

### Reports
- `phase13_final_report.txt`: Universe sizes, module info, candidate class counts (pair + unique gene), forbidden-term check, validation summary.
- `phase13_provenance.csv`: Machine-readable provenance record.
- `phase13_provenance.txt`: Human-readable provenance record.

## How to Read Correlation / DE / Regulator Score

- **Spearman ρ** ∈ [-1, 1]. |ρ| > 0.5 is "moderate-to-strong" at n≥8. At n=11-12 per region, p-values and q-values are **exploratory only** — they have severely limited power.
- **DE log2FC**: Positive = Old higher (Phase 12) or CA3 higher (Phase 11).
- **regulator_score**: A relative ranking heuristic within a module×region. Higher = more lines of evidence. It is NOT a p-value, NOT a probability, NOT a measure of statistical significance.

## Why Candidate Regulator ≠ Causal Regulator

Correlation is symmetric. A gene correlated with a module score could mean:
- The gene regulates mitochondrial genes (hypothesis)
- Mitochondrial status affects the gene's expression
- A third variable (e.g., cell composition) drives both
- Technical artifact (co-expression in same cell types)

No perturbation experiment or mechanistic validation was performed. The word "candidate" is intentional: these are hypotheses, not conclusions.

## Self-Correlation Strategy A

Module gene-set members are excluded from their own module's regulator ranking.
A gene can still be ranked for modules where it is NOT a member.

**3 self_correlation_flag values**:
- `regulator_candidate`: Gene is NOT a member of the correlated module's gene set. Eligible for ranking.
- `module_member_or_effector`: Gene IS a member of the module's gene set. ρ is self-correlation and cannot be interpreted as upstream regulator evidence. Rank is NA for own module.
- `target_gene_candidate`: Gene is in target_genes.csv but NOT a module member. Reported in separate target gene table.

**candidate_class = low_confidence**: Applied when |ρ| < 0.3 in both regions. This is NOT a self_correlation_flag value — it is a separate classification for genes below the screening threshold.

## Metric Hierarchy

| Metric | Role | Threshold |
|--------|------|-----------|
| |ρ| | PRIMARY ranking metric (50% weight) | >0.3 screening; >0.5 notable |
| DE_age |log2FC| | Secondary ranking (30% weight) | Used as ranking weight only |
| DE_region |log2FC| | Secondary ranking (20% weight) | Used as ranking weight only |
| q-value (BH) | EXPLORATORY only | NOT used for candidate selection |
| TF annotation | Eligibility for universe_regulator | NOT a ranking weight |
| Target gene membership | Annotation only | NOT a ranking weight |

## Limitations

1. All correlations at GEMgroup (pseudobulk) level: CA1 n=12, CA3 n=11.
2. CA3 Young n=3: severely underpowered; results carry caveat label.
3. Visium spots capture multiple cells: spatial heterogeneity within spots is unmeasured.
4. 16 Complex V (Atp5*) genes missing from Hippo data.
5. mtDNA-encoded module has gene set but NO sample-level score — excluded from correlation analysis.
6. Module scores and gene expression share the same underlying count data.
7. RDS-based workflow (not raw Space Ranger/Visium files).
8. Age and GEMgroup are colinear: cannot separate age from individual effects.
9. Mouse genome: gene name mapping from human reference databases may have gaps.
10. CollecTRI unavailable (Ensembl SSL certificate expiry); GO-TF used as fallback.

## Next Steps

1. Review candidate regulator-associated genes per module (score_default, rank_CA1/CA3).
2. Check sensitivity configs for ranking stability.
3. Examine self_correlation_flag — module members should NOT be interpreted as regulators of their own module.
4. Use spatial plots (F14) for expression pattern validation.
5. Cross-reference with published literature for biological plausibility.
6. Remember: these are **hypotheses for follow-up**, not conclusions.


---

## Phase 13b: CA3 Symmetric Display (Added 2026-06-05)

Phase 13b corrected the display asymmetry in Phase 13. The underlying data (rho, scores, ranks)
was already computed for both CA1 and CA3 in Phase 13. Phase 13b only added CA3-focused display outputs.

### CA1 vs CA3 Symmetric Outputs

**CA1-focused figures** (existing, from Phase 13):
- F01: regulator_score heatmap (sorted by rank_CA1)
- F03: top regulators dotplot (sorted by rank_CA1, sized by |rho_CA1|)
- F04: module correlation matrix (rho_CA1, top 50 genes)
- F07: Young stratified heatmap (CA1 samples, n=12)
- F08: Old stratified heatmap (CA1 samples, n=12)
- F09: score component breakdown (CA1 scaled components)
- F10: sensitivity comparison (CA1 default/de-emphasis/rho-only)
- F13: DE volcano (CA1 Old vs Young)
- F14: spatial plots (top 5 by score_CA1_default)

**CA3-focused figures** (added by Phase 13b):
- F01b: CA3 regulator_score heatmap (sorted by rank_CA3)
- F03b: CA3 top regulators dotplot (sorted by rank_CA3, sized by |rho_CA3|)
- F04b: CA3 module correlation matrix (rho_CA3, top 50 genes)
- F07b: CA3 Young stratified heatmap (n=3, CAUTION: severely underpowered)
- F08b: CA3 Old stratified heatmap (n=11)
- F09b: CA3 score component breakdown
- F10b: CA3 sensitivity comparison
- F13b: CA3 DE volcano (CA3 Old vs Young)
- F14b: CA3 spatial plots (top 5 by score_CA3_default)

**Region-agnostic figures** (unchanged):
- F02: rho vs DE scatter (has both CA1+CA3 panels)
- F05: self-correlation flag distribution
- F06: candidate class distribution
- F11: TF annotation coverage
- F12: universe composition

### Top Candidate CSVs (Phase 13b)

- `phase13_top_candidates_CA1.csv`: sorted by rank_CA1
- `phase13_top_candidates_CA3.csv`: sorted by rank_CA3
- `phase13_top_candidates_concordant.csv`: sorted by combined_score (CA1+CA3)
- `phase13_top_candidates_region_flipped.csv`: sorted by combined_score

**CAUTION**: CSV rows are gene x module pairs, not unique genes.

### Caveats

- CA3 Young n=3 is severely underpowered. F07b results should be interpreted with extreme caution.
- All candidate regulator associations are correlational (Spearman rho), not causal.
- Complex V module excluded (all 16 genes missing from dataset).

---




---


---


---

## Phase 13c: Age-Stratified Candidate Regulator Discovery

### What Phase 13c Adds

Phase 13/13b ranked candidate regulator-associated genes using all-age Spearman rho across CA1 (n=12) and CA3 (n=11).
Phase 13c extends this with age-stratified analysis, computing separate rho values for Young and Old strata within each region.

### Class Columns

Phase 13c uses three class columns:

- `primary_age_class`: Within-region age-specific classification (CA1_old_specific, CA1_young_specific, CA1_age_conserved, CA1_age_shifted, CA3_old_specific, CA3_young_specific, CA3_age_conserved, CA3_age_shifted, low_confidence, module_member_or_effector).
- `region_age_class`: Cross-region within-age classification (Young_CA1_specific, Young_CA3_specific, Young_region_conserved, Old_CA1_specific, Old_CA3_specific, Old_region_conserved, none).
- `candidate_class`: Same as primary_age_class (backward compatibility).
- `secondary_class`: Additional classes when a pair falls into multiple categories (semicolon-separated).

### Sample Size Labels

Per-stratum labels (in the main table):
- `sample_size_label_CA1_Young` = exploratory_small_n (n=4)
- `sample_size_label_CA1_Old` = exploratory_moderate_n (n=8)
- `sample_size_label_CA3_Young` = severely_underpowered_exploratory (n=3)
- `sample_size_label_CA3_Old` = exploratory_moderate_n (n=8)

In specialized output tables (old_specific, young_specific, etc.), the `sample_size_label` column uses the label for the relevant stratum.

### How to Read Age-Stratified Tables

- Each row is a gene x module pair with separate rho values for Young and Old.
- Compare rho_Young vs rho_Old within the same gene-module to identify age-specific patterns.
- `candidate_score_<stratum> = |rho_<stratum>|` is the primary ranking metric (NOT a p-value).
- `pvalue_exploratory_*` and `qvalue_exploratory_*` columns are computed but NOT used for ranking or class assignment.

### How to Interpret Rho Direction

- Positive rho: gene expression and module score increase together.
- Negative rho: gene expression increases while module score decreases.
- In Old strata, this may reflect aging-related coupling changes.

### Why Small n Matters

- CA1 Young n=4 and CA3 Young n=3 are underpowered for Spearman correlation.
- rho estimates in these strata are noisy. Compare with Old strata (n=8) for more stable estimates.
- CA3 Young n=3: minimum p-value ~0.33 for perfect monotone correlation. Do NOT use CA3 Young p-values for filtering.

### Candidate Classes

- Age-specific (e.g., CA1_old_specific): strong rho in one age but NOT the other within same region.
- Age-conserved: strong rho in both Young and Old, same direction.
- Age-shifted: strong rho in both Young and Old, OPPOSITE direction.
- Region-specific within age (e.g., Old_CA1_specific): strong rho in one region but NOT the other within same age.
- `low_confidence`: |rho| < 0.3 in ALL 4 strata.
- `module_member_or_effector`: gene is a member of the module gene set. Self-correlation excluded from regulator ranking.

### Key Files

- `phase13c_age_stratified_rho.csv` -- full rho table
- `phase13c_age_stratified_scores.csv` -- scores with per-stratum sample_size_label columns
- `phase13c_age_specific_candidate_classes.csv` -- all class columns (primary_age_class, region_age_class, candidate_class, secondary_class)
- `phase13c_old_specific_candidates.csv` -- Old-specific candidates
- `phase13c_young_specific_candidates.csv` -- Young-specific candidates
- `phase13c_age_shifted_candidates.csv` -- candidates with rho sign reversal
- `phase13c_region_shifted_by_age_candidates.csv` -- region-specific candidates within same age
- `phase13c_top_candidates_by_region_age.csv` -- top 20 per stratum

### Caveats

1. These are ASSOCIATION-BASED candidate regulator-associated genes, NOT causal regulators.
2. Age-stratified rho is still correlation. Age specificity adds biological plausibility but does NOT establish causality.
3. CA3 Young n=3 -- all CA3 Young results carry `severely_underpowered_exploratory` label.
4. CA1 Young n=4 -- all CA1 Young results carry `exploratory_small_n` label.

