# Spatial DESeq2 QC Backfill Plan

> **Status**: PLAN / MEMO ONLY — not executed.
> **Scope**: Future lightweight QC backfill for Phase 03, Phase 11, and Phase 11b DESeq2 outputs.
> **Do not run from this document**: This file records a recommended follow-up only. It does not change prior conclusions by itself.

## Purpose

Phase 03, Phase 11, and Phase 11b used GEMgroup-level pseudobulk DESeq2 designs and remain valid exploratory/statistical outputs. A later review against local bio-* differential-expression guidance identified one missing QC layer: explicit diagnostics for `padj = NA` and p-value distribution checks.

This backfill is intended to improve auditability, especially for user target genes that may have missing adjusted p-values.

## Phases Covered

| Phase | Script | Existing DE scope | Backfill priority |
|---|---|---|---|
| Phase 03 | `R/spatial/s03_pseudobulk_deseq2_dg.R` | DG pseudobulk DESeq2 reproduction | Medium |
| Phase 11 | `R/spatial/s11_ca1_ca3_de_module_coupling.R` | CA3 vs CA1 across all ages | High |
| Phase 11b | `R/spatial/s11b_age_grouping_ca1_ca3_de.R` | CA3 vs CA1 within Young, Middle, Old | High |

## Recommended Backfill Outputs

Use a new generated output directory, for example:

`data/processed/spatial/deseq2_qc_backfill/`

Recommended files:

| Output | Description |
|---|---|
| `phase03_padj_na_diagnostics.csv` | Counts of total genes, `pvalue=NA`, `padj=NA`, independent-filtering-like rows, and affected target genes if applicable |
| `phase11_padj_na_diagnostics.csv` | Same diagnostics for all-age CA3 vs CA1 |
| `phase11b_padj_na_diagnostics_by_age.csv` | Same diagnostics for Young, Middle, and Old CA3 vs CA1 |
| `phase03_pvalue_histogram.png` | DESeq2 p-value histogram QC |
| `phase11_pvalue_histogram.png` | DESeq2 p-value histogram QC |
| `phase11b_pvalue_histogram_Young.png` | Age-specific p-value histogram QC |
| `phase11b_pvalue_histogram_Middle.png` | Age-specific p-value histogram QC |
| `phase11b_pvalue_histogram_Old.png` | Age-specific p-value histogram QC |
| `deseq2_qc_backfill_summary.txt` | Human-readable summary and caveats |

## Diagnostic Logic

For each DESeq2 result table, record:

| Field | Meaning |
|---|---|
| `n_genes_total` | Number of result rows |
| `n_pvalue_na` | Number of rows with `is.na(pvalue)` |
| `n_padj_na` | Number of rows with `is.na(padj)` |
| `n_padj_na_with_pvalue_present` | Rows with p-value present but adjusted p-value absent; often consistent with independent filtering |
| `n_pvalue_na` | Rows where p-value is absent; may reflect all-zero rows, Cook's outlier handling, or other DESeq2 filtering behavior |
| `n_target_genes_padj_na` | Number of user target genes with missing adjusted p-values |
| `target_genes_padj_na` | Semicolon-separated target genes with `padj=NA`, if small enough to report |

Do not over-interpret the cause of missing values without the DESeq2 object or sufficient saved metadata. If the exact cause cannot be established from result tables alone, label it as `inferred_from_result_table`.

## Plot QC

For each DE result table, plot a p-value histogram using rows with non-missing p-values.

Interpretation should remain cautious:

- Roughly uniform p-values can be acceptable when few genes are truly differential.
- Enrichment near zero is expected when many genes differ.
- A strong U-shape, spike at one, or unusual distribution should be flagged for review.
- Histogram QC is diagnostic only; it is not a biological result.

## Execution Boundaries

Future execution of this backfill should:

- Not rerun heavy spatial analyses unless required.
- Prefer using existing DE result CSV files.
- Not load full Seurat objects.
- Not save full Seurat or DESeq2 objects.
- Not modify prior Phase 03, Phase 11, or Phase 11b result files.
- Write only new QC outputs under gitignored generated-output directories.
- Record provenance, including input file paths, file modification times, package versions if R is used, and whether `renv.lock` changed.

## Why This Matters

The local bio-* differential-expression guidance emphasizes that `padj=NA` rows can arise from multiple DESeq2 mechanisms, including independent filtering, Cook's outlier handling, or all-zero rows. Without explicit diagnostics, an important target gene could appear absent from interpreted results without a clear reason.

This backfill is not expected to overturn Phase 03/11/11b results. It improves traceability before biological interpretation.

## Relationship to Phase 12

Phase 12 should include these diagnostics from the start:

- explicit `alpha = 0.05` reporting,
- `padj=NA` diagnostics by region,
- p-value histogram QC by region,
- clear distinction between significance (`padj < 0.05`) and post-hoc effect-size flags such as `abs(log2FC) > log2(1.5)`.

