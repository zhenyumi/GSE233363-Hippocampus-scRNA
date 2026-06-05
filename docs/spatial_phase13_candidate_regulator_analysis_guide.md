# Phase Spatial-13: Candidate Regulator Analysis Guide

Generated: 2026-06-05 12:29:59
Plan version: v1.2

## Overview
This guide describes the Phase Spatial-13 candidate regulator discovery results.
All candidate regulators are ASSOCIATED with mitochondrial modules — correlation ≠ causation.

## Key Outputs

### Primary Rankings
- `phase13_regulator_score_default.csv`: Primary ranking by default config (ρ 50% + DE_age 30% + DE_region 20%)
- `phase13_candidate_classes.csv`: Classification of candidates (CA1/CA3 aging-associated, concordant, region-flipped, etc.)
- `phase13_spearman_rho_all.csv`: All Spearman ρ values with BH-corrected q-values

### Sensitivity
- `phase13_regulator_score_sensitivity.csv`: All 3 weight configs (default, DE-emphasis, ρ-only)
- Compare configs to assess robustness of rankings

### Annotation
- `phase13_tf_annotation_coverage.csv`: TF sources (DoRothEA A-C, CollecTRI, GO:0003700/GO:0140110)
- `phase13_target_gene_module_association.csv`: Target gene module correlations
- `phase13_self_correlation_flags.csv`: Self-correlation flags (module member exclusion)

## Interpretation Caveats

### Correlation ≠ Causation
All candidates are ASSOCIATED regulators. Spearman ρ measures monotonic association, not causal regulation.
DoRothEA/CollecTRI annotation = eligibility for universe_regulator, NOT proof of regulation in this context.

### Small Sample Size
- CA1: n = 12 (4 Young + 8 Old)
- CA3: n = 11 (3 Young + 8 Old) — CA3 Young n=3 caveat applies
- Age-stratified correlations have even smaller n and should be interpreted cautiously

### Metric Hierarchy
- |ρ| is the PRIMARY metric (weight 50% in default config)
- q-values are EXPLORATORY only (BH-corrected, but small n limits power)
- |ρ| > 0.3: screening threshold
- |ρ| > 0.5: notable threshold

### Self-Correlation Strategy A
Module gene-set members are excluded from their own module's regulator ranking.
4 flag values: regulator_candidate, module_member_or_effector, target_gene_candidate, low_confidence_self_correlation

## Figures
- F01: regulator_score heatmap (top 20 per module)
- F02: ρ vs DE_age scatter (CA1 and CA3)
- F03: Top regulators dotplot
- F04: Gene × Module ρ matrix
- F05: Self-correlation flag distribution
- F06: Candidate class distribution
- F07-F08: Age-stratified ρ heatmaps (Young, Old)
- F09: Score component breakdown
- F10: Sensitivity config comparison
- F11: TF annotation coverage
- F12: Universe composition
- F13: DE volcano with regulators highlighted
- F14: SpatialFeaturePlot for top candidates
- F15: GO enrichment (if available)

## Next Steps
1. Review top-ranked candidates per module (score_default, rank_CA1/CA3)
2. Check sensitivity configs for ranking stability
3. Examine self-correlation flags — module members should NOT be interpreted as regulators of their own module
4. Use spatial plots (F14) for expression pattern validation
5. Cross-reference with published literature for biological plausibility

