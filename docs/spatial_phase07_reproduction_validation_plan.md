# Phase Spatial-07: Reproduction Validation & Gap Report — Plan v2.1

**Version:** 2.1 (approved 2026-06-04)
**Status:** APPROVED — EXECUTING
**Git HEAD at plan time:** 6ca46a2
**renv.lock MD5 at plan time:** 139bcfaf7ce24b633c9f034380a15ef5

---

## 1. Purpose

Phase Spatial-07 is the final validation phase of the spatial transcriptomics reproduction pipeline (Phases 01–06). It:

- Classifies every Phase 01–06 output using a revised fine-grained system
- Classifies every figure/panel listed in the spatial reproduction plan
- Documents all identified discrepancies with source citations
- Evaluates handoff criteria for Phase 09 (mitochondrial analysis)
- Evaluates Tangram decision criteria

Phase 07 performs **no new analysis**. It reads only CSV/TXT summaries from prior phases and generates report documents.

---

## 2. Classification Labels (Revised System)

| Label | Definition |
|-------|-----------|
| OBJECT_INSPECTION | Read-only structural inspection of author-provided objects. No method reproduction attempted. Validates format, metadata, assays, images, reductions. |
| OBJECT_DERIVED_VALIDATION | Validation of metadata/annotations extracted from author-provided objects against expected values from prior inspection or author code. No new processing pipeline. |
| STRICT_REPRODUCTION | Same method, same input source, same tool path. Outputs directly comparable to author results with no critical warnings. |
| METHOD_MATCH_WITH_WARNING | Method lineage matches author code, but version differences (Seurat, msigdbr), gene availability, count discrepancies, or threshold differences introduce a documented warning. |
| RDS_APPROXIMATION | Uses RDS-contained coordinates/metadata as approximation for a method that originally required raw Space Ranger/Visium files, STutility, or image-mask processing. Results directionally consistent but not directly equivalent. |
| PARTIAL_REPRODUCTION | Only a subset of panels or statistics were reproduced; full figure panel not generated. |
| BLOCKED_RAW | Requires raw Visium/Space Ranger files — not available. |
| BLOCKED_WHOLEBRAIN | Requires WholeBrain RDS — not loaded per current plan scope. |
| BLOCKED_TANGRAM | Requires Tangram cell-type deconvolution — not executed. |
| DEFERRED | Not attempted; may be optional, low priority, or requires tools not in scope (e.g., Monocle 2, GO enrichment, CAS scoring). |
| NOT_IN_SCOPE | Panel is not spatial or not relevant to hippocampus/DG reproduction. |

---

## 3. Phase Classification

| Phase | Script | Classification | Rationale |
|-------|--------|---------------|-----------|
| 01 | s01_inspect_objects.R | OBJECT_INSPECTION | Reads RDS, records structure, no method reproduction |
| 02 | s02_dg_hippo_metadata.R | OBJECT_DERIVED_VALIDATION | Validates metadata against s01 expectations; 109/109 PASS |
| 03 | s03_pseudobulk_deseq2_dg.R | METHOD_MATCH_WITH_WARNING | Method matches Author Script 7 (aggregate.Matrix + DESeq2). DEG: OY=677, MY=253, OM=102 padj. No author reference DEG counts found. |
| 04 | s04_ifng_module_score.R | METHOD_MATCH_WITH_WARNING | Method matches Script 8 lines 1-66. Middle IS=48 vs author 60, Old IS=317 vs author 299. |
| 05 | s05_rds_inflammatory_gradient.R | RDS_APPROXIMATION | Uses RDS coordinate neighbor graph instead of STutility::RegionNeighbours. |
| 05b | s05b_plot_spatial_labels.R | RDS_APPROXIMATION | Visualization of Phase 05 labels. |
| 06 | s06_hippo_regional_atlas.R | OBJECT_DERIVED_VALIDATION | Validates Hippo regional structure from RDS metadata. 54/54 PASS. |

---

## 4. Figure/Panel Classification

| Figure | Panel | Phase | Classification | Rationale |
|--------|-------|-------|---------------|-----------|
| Fig. 1 | f | — | BLOCKED_WHOLEBRAIN | Requires WholeBrain RDS (42,169 spots) |
| Fig. 1 | g | — | BLOCKED_WHOLEBRAIN | WholeBrain spatial projection |
| Fig. 1 | h | — | DEFERRED | Hippocampal CCA mapping |
| Fig. 5 | d | 02, 06 | OBJECT_DERIVED_VALIDATION | DG UMAP from DG RDS metadata; Phase 06 UMAP by Region on Hippo |
| Fig. 5 | e | — | DEFERRED | Marker expression spatial mapping |
| Fig. 5 | f | 03 | PARTIAL_REPRODUCTION | Phase 03 has DEG; CAS score not computed |
| Fig. 5 | g | — | DEFERRED | CAS violin not computed |
| Fig. 5 | h | 03 | PARTIAL_REPRODUCTION | Phase 03 DESeq2 contrasts OY/MY/OM |
| Fig. 6 | e | 04, 05 | RDS_APPROXIMATION | Phase 04 IS proportions; Phase 05 label maps |
| Fig. 7 | a | 04 | METHOD_MATCH_WITH_WARNING | Same method (msigdbr + AddModuleScore), but Phase 04 has IS count discrepancy |
| Fig. 7 | c | 05, 05b | RDS_APPROXIMATION | RDS coordinate graph, not STutility image masks |
| Fig. 7 | d | 05 | RDS_APPROXIMATION | Label counts from coordinate-based propagation |
| Fig. 7 | e | — | BLOCKED_RAW | STutility + Monocle 2 pseudotime |
| Fig. 7 | f | — | BLOCKED_RAW | Same as 7e |
| Fig. 7 | g | — | BLOCKED_RAW | Same as 7e |
| Fig. 7 | h | — | BLOCKED_RAW | Same as 7e |
| Fig. 7 | i | — | BLOCKED_RAW | Same as 7e |
| Fig. 7 | j | — | BLOCKED_RAW | Same as 7e |
| Ext. Data Fig. 2 | a | 02 | OBJECT_DERIVED_VALIDATION | UMAP from metadata |
| Ext. Data Fig. 2 | b | — | DEFERRED | |
| Ext. Data Fig. 2 | c | — | DEFERRED | |
| Ext. Data Fig. 2 | d | 06 | PARTIAL_REPRODUCTION | Region validation only; no Allen overlay |
| Ext. Data Fig. 2 | e | — | BLOCKED_TANGRAM | |
| Ext. Data Fig. 7 | c | 03 | METHOD_MATCH_WITH_WARNING | DESeq2 method matches; no author reference DEG count |
| Ext. Data Fig. 7 | d | 03 | METHOD_MATCH_WITH_WARNING | Same as above |
| Ext. Data Fig. 9 | a | 04 | METHOD_MATCH_WITH_WARNING | IS proportions; Middle/Old discrepancy |
| Ext. Data Fig. 9 | b | 04 | METHOD_MATCH_WITH_WARNING | IS histogram; Phase 04 has Middle/Old IS count discrepancy, Seurat/msigdbr differences |
| Ext. Data Fig. 9 | c | 05 | RDS_APPROXIMATION | Old IS vs NNS DESeq2 (Up=2571/Down=1452 vs 2357/1326) |
| Ext. Data Fig. 9 | d | — | DEFERRED | GO enrichment |
| Ext. Data Fig. 9 | e | — | DEFERRED | |
| Ext. Data Fig. 9 | f | 05 | RDS_APPROXIMATION | |
| Ext. Data Fig. 9 | g | 05 | RDS_APPROXIMATION | |
| Ext. Data Fig. 9 | h | — | DEFERRED | CAS scoring |
| Ext. Data Fig. 9 | i | — | DEFERRED | CAS scoring |
| Ext. Data Fig. 10 | a-d | — | DEFERRED | SLM manual selection, BBB, correlation — not in scope |

---

## 5. Discrepancy Table

| ID | Phase | Metric | Our Value | Author/Expected Value | Source | Severity | Possible Causes |
|----|-------|--------|-----------|----------------------|--------|----------|----------------|
| D01 | 04 | Middle IS spots | 48 (2.00%) | 60 (2.50%) | Script 8 line 41; Phase 04 summary | warning | msigdbr version, SCT assay gene differences, Seurat v4.3 vs v5.5.0 |
| D02 | 04 | Old IS spots | 317 (6.29%) | 299 (5.94%) | Script 8 line 41; Phase 04 summary | warning | Same as D01 |
| D03 | 05 | IS vs NNS DEG Up (Old) | 2571 | 2357 | Phase 05 provenance; Script 8 reference | caution | RDS coordinate margin vs STutility image masks |
| D04 | 05 | IS vs NNS DEG Down (Old) | 1452 | 1326 | Phase 05 provenance; Script 8 reference | caution | Same as D03 |
| D05 | 05 | Neighbor graph method | RDS Euclidean coordinates | STutility::RegionNeighbours | Phase 05 provenance | blocker | Raw Visium files unavailable |
| D06 | ALL | Seurat version | v5.5.0 | v4.3 | renv.lock; author headers | caution | Assay slot structure, SCT behavior |
| D07 | ALL | Raw Visium files | Not available | Required by Script 8 | data/raw/spatial/ is empty | blocker | Not on Zenodo |
| D08 | 05 | Label margin heuristic | 1.2x mean NN distance | STutility image propagation | Phase 05 provenance | caution | No external calibration |
| D09 | 03 | DEG OY count reference | 677 padj | No author reference | Phase 03 deg_summary.csv; Script 7 | info | Author Script 7 does not report total DEG counts |

---

## 6. Handoff Decision Criteria

Possible conclusions:
- GO_WITH_CAVEATS: Phase 06 Hippo regional framework sufficient for Phase 09 planning
- NO_GO: Critical checks failed
- NEED_PHASE08_TANGRAM: Cell-type deconvolution required before interpretation
- NEED_RAW_VISIUM: Strict STutility route required

Evaluation logic:
- If Phase 06 PASS AND Phase 04/05 discrepancies do not affect CA1/CA3/DG region labels → GO_WITH_CAVEATS
- If Phase 06 non-PASS → NO_GO
- If region-level mitochondrial planning is only near-term goal → Tangram may be deferred
- If any critical data path unavailable → NO_GO

---

## 7. Tangram Decision Framework

Criteria to evaluate:
- Is Tangram required for CA1/CA3/DG region-level mitochondrial analysis? (Likely no — regions from metadata)
- Is Tangram needed to validate the Hippo regional framework itself? (No — Phase 06 validates)
- Would Tangram output change interpretation? (Possibly — if cell-type composition differs by region)

---

## 8. Rplots.pdf Handling

1. Record whether Rplots.pdf exists at script start
2. If absent at start but appears during execution: delete the new file, log as WARNING
3. If present before s07 starts: do not delete; report and request user decision

---

## 9. s07 Script Structure (12 Steps, Base R)

```
E0:  Setup (paths, git HEAD, renv.lock MD5, Rplots.pdf guard start)
E1:  Read Phase 01-06 validation summaries + deg_summary.csv
E2:  Read author Script 8 expected values (hardcoded IS counts)
E3:  Phase classification → phase07_phase_classification.csv
E4:  Figure/panel classification → phase07_figure_classification_table.csv
E5:  Discrepancy table D01-D09 → phase07_discrepancy_table.csv
E6:  Handoff evaluation → one of GO_WITH_CAVEATS / NO_GO / NEED_PHASE08 / NEED_RAW
E7:  Tangram evaluation → TANGRAM_DEFERRED / TANGRAM_RECOMMENDED
E8:  Write all CSV outputs
E9:  Generate docs/spatial_figure_validation.md
E10: Generate docs/spatial_reproducibility_report.md
E11: Generate docs/spatial_handoff_to_mitochondrial_analysis.md
E12: Write validation_summary.txt + provenance.csv; Rplots.pdf guard end; gc()
```

---

## 10. Outputs

| Output | Path | Type |
|--------|------|------|
| Plan doc | docs/spatial_phase07_reproduction_validation_plan.md | Plan |
| Script | R/spatial/s07_figure_validation.R | Script (base R) |
| Figure validation | docs/spatial_figure_validation.md | Report |
| Reproducibility report | docs/spatial_reproducibility_report.md | Report |
| Handoff document | docs/spatial_handoff_to_mitochondrial_analysis.md | Report |
| Phase classification | data/processed/spatial/phase07_validation/phase07_phase_classification.csv | CSV |
| Figure classification | data/processed/spatial/phase07_validation/phase07_figure_classification_table.csv | CSV |
| Discrepancy table | data/processed/spatial/phase07_validation/phase07_discrepancy_table.csv | CSV |
| Validation summary | data/processed/spatial/phase07_validation/phase07_validation_summary.txt | TXT |
| Provenance | data/processed/spatial/phase07_validation/phase07_provenance.csv | CSV |

---

## 11. Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|-------------|
| AC1 | Plan doc saved before script created | File timestamps |
| AC2 | Script parse-checks | Rscript parse check |
| AC3 | Script runs to completion | Exit code 0 |
| AC4 | spatial_figure_validation.md produced | Complete matrix |
| AC5 | spatial_reproducibility_report.md produced | Phase + discrepancy tables with citations |
| AC6 | spatial_handoff_to_mitochondrial_analysis.md produced | Exactly one decision with rationale |
| AC7 | phase07_discrepancy_table.csv has D01-D09 | 9 rows with source citations |
| AC8a | phase07_phase_classification.csv has 7 phases | 7 rows (s01-s06, s05b) |
| AC8b | phase07_figure_classification_table.csv covers all panels | No missing panels |
| AC9 | No new analysis performed | Provenance confirms read-only |
| AC10 | No package installed, renv.lock unchanged | MD5 verified |
| AC11 | No RDS objects loaded | Provenance records zero RDS loads |
| AC12 | Every conclusion cites source | Spot-check |

---

## 12. Stop Conditions

| SC | Condition | Action |
|----|-----------|--------|
| SC-1 | Phase 01-06 validation file missing | STOP |
| SC-2 | Author Scripts 6-8 unreadable | CAUTION; continue with available evidence |
| SC-3 | Phase 06 non-PASS | STOP |
| SC-4 | Phase 03 missing validation summary | CAUTION; classify as METHOD_MATCH_WITH_WARNING |
| SC-5 | renv.lock MD5 mismatch | CAUTION; record, continue |
