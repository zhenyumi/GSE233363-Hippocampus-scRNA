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

