# AGENTS.md - ferroDG Reproduction Project

## Project Overview
Reproduce ferroDG (https://github.com/ShilongZhang116/ferroDG) stage-2 analysis:
ferroptosis analysis in mouse dentate gyrus neurogenesis using scRNA-seq (GSE233363).

## Key Context
- **Species**: Mouse (gene symbols: First letter uppercase, rest lowercase)
- **Seurat version**: v5.5.0 (handles `avg_log2FC` vs `avg_logFC`)
- **Cell types**: qNSC, nIPC, Neuroblast, GC (neurogenic lineage)
- **Comparison**: Young vs Old (Middle timepoint excluded)
- **Environment**: renv-managed R 4.5.1, macOS

## Scripts (run in order)
1. `R/scrna/01_load_and_subset.R` - Load custom data, collapse age, subset neurogenic cells
2. `R/scrna/02_module_scores.R` - Calculate ferroptosis module scores
3. `R/scrna/03_figures_0_to_4.R` - UMAP, violin, line, trajectory plots
4. `R/scrna/04_de_and_volcano.R` - Differential expression + volcano plots
5. `R/scrna/05_correlation_heatmaps.R` - Spearman correlation heatmaps
6. `R/scrna/00_download_public_data.R` - Download official author-provided RDS from Zenodo

## When Modifying Figures
- Color schemes defined in each script (look for `scale_color_manual` or `color =`)
- Use `cairo_pdf()` with `tryCatch` fallback to `pdf()` (X11 not installed)
- Save both PNG (always works) and PDF (best effort)
- Call `gc()` after each large plot to manage memory

## When Adding New Analysis
- For scRNA-seq pipeline additions, create `R/scrna/06_<descriptive_name>.R`.
- For spatial transcriptomics additions, create `R/spatial/s##_descriptive_name.R`; do not add spatial scripts to `run_pipeline.R`.
- For optional Python spatial work, create source under `python/spatial/` only after an approved Python/Tangram plan.
- Follow pattern: load data -> compute -> save -> plot -> gc()
- Add only scRNA-seq pipeline scripts to `run_pipeline.R`.
- Update this file and the relevant README/plan document with new script descriptions.

## Gene Symbol Convention
- Original paper uses Human symbols (uppercase)
- Convert to Mouse: `tolower()` + capitalize first letter
- Example: "TFRC" -> "Tfrc", "SLC7A11" -> "Slc7a11"

## Known Missing Genes (expected)
- Inhibitor: Mt1g, Akr1c1-3, G6pd (not in mouse genome under these symbols)
- Regulator: Tp53, Ras (not in mouse genome under these symbols)
- Scripts handle these with guards (skip if missing)

## Data Paths
- Custom data: `data/raw/GSE233363_custom/` (gitignored, independently aligned from FASTQ — USED by pipeline)
- Official data: `data/raw/GSE233363_official/` (gitignored, author-provided Zenodo files — for reference)
- Processed: `data/processed/` (gitignored, generated)
- Report bundles: `data/report/` (gitignored, generated markdown/PDF integration reports)
- Figures: `figures/stage-2/` (gitignored, generated)
- Gene lists: `gene_lists/` (committed)

## File Structure and Temporary Outputs
- `docs/repository_structure.md` is the canonical repository-layout reference. Update it when directory conventions change.
- Keep source code and generated artifacts separate: executable scripts belong in `R/`, committed docs in `docs/`, and generated outputs in `data/`, `figures/`, `results/`, `reports/`, or `cache/`.
- Active code roots are `R/scrna/`, `R/spatial/`, and the optional `python/spatial/` directory. Do not place executable scripts directly in `R/`.
- Spatial scripts belong in `R/spatial/` and use the `s##_` prefix.
- Optional Python scripts belong in `python/spatial/`; Python virtual environments, caches, notebooks checkpoints, H5AD/AnnData files, and generated outputs must not be committed.
- `analysis/` is retired and must not be used for new scripts, notebooks, scratch files, or generated outputs. If an empty local `analysis/` directory appears, delete it rather than adding `.gitkeep`.
- Do not write temporary files, downloaded data, logs, plots, or intermediate tables into the repository root, `R/`, `R/scrna/`, `R/spatial/`, `python/`, `docs/`, or `gene_lists/` unless they are intentionally reviewed source/docs.
- For scRNA-seq outputs, use existing conventions: processed objects/tables in `data/processed/` and figures in `figures/stage-2/`.
- For spatial outputs, use spatial-specific subdirectories: `data/processed/spatial/`, `figures/spatial/`, `results/spatial/`, `reports/spatial/`, and `cache/spatial/`.
- Short-lived scratch files should go under `cache/` or `cache/spatial/` and be removed when no longer needed.
- Lightweight audit logs and validation reports may go under `reports/` or `reports/spatial/`; committed planning/reference docs stay under `docs/`.
- Large or generated files must remain gitignored. Do not commit RDS/RData, H5/H5AD, MTX, image tiles, compressed raw data, cache files, generated result directories, or root-level `Rplots.pdf`.
- Before adding a new output path, check existing project conventions and `.gitignore`; prefer adding a narrow subdirectory over scattering files across the tree.

## Memory and Large-Object Management
- Prefer staged scripts over monolithic analyses; each script should load only the inputs needed for that stage, write a checkpoint or lightweight summary, then release memory.
- Avoid loading multiple large Seurat or Visium objects at the same time. For spatial work, do not load WholeBrain, Hippo, and DG objects together unless there is a documented need.
- After saving a stage output, remove large temporary objects with `rm()` and call `gc()` before the next major operation.
- Avoid unnecessary object copying, especially Seurat object copies created by repeated assignment, broad subsetting, or adding temporary metadata columns.
- Preserve sparse matrices where possible; do not coerce large count matrices to dense matrices unless explicitly justified and memory-safe.
- Save lightweight CSV/TXT summaries for inspection and validation when full object copies are not required.
- Do not save duplicate large RDS files without a clear reason, provenance note, and gitignored destination.
- On macOS, assume available RAM may be the limiting factor for large Seurat/Visium objects; record memory failures and retry with smaller staged inspections rather than forcing the full workflow.

## Running
```bash
# Individual scripts
Rscript R/scrna/01_load_and_subset.R
Rscript R/scrna/02_module_scores.R
# ...

# Full pipeline
Rscript run_pipeline.R
```

## Debugging Tips
- If `FetchData` fails: check gene names exist in `rownames(seurat_obj)`
- If `pheatmap` errors on NA: add `cor_mat[is.na(cor_mat)] <- 0`
- If PDF fails: X11/cairo not installed, PNGs still work
- Memory: call `gc()` after large operations

## Spatial Transcriptomics Extension (Phase Spatial-00)
**Status**: Stages Spatial-01 through Spatial-11/11b implemented for DG/Hippo. Phase Spatial-05 is an RDS-based reproduction/approximation, not strict STutility reproduction. Phase Spatial-07 handoff: GO_WITH_CAVEATS. Phase Spatial-10/11: CA1 vs CA3 target gene, module score, and age-stratified DE analysis complete. Phase Spatial-13/13b/13c: candidate regulator-associated gene discovery complete (all-age + age-stratified).

### Key Distinctions
- **scRNA-seq cells ≠ Visium spots**: Visium spots may capture signal from multiple cells and are not equivalent to single cells.
- **Spots are not independent biological replicates**: spatial neighbors share local microenvironment
- **Metadata must be inspected from actual objects**: do not assume assay names, spatial keys, image slots, or anatomical labels
- **Consult ref-bio + bio-* skills for authoritative reference routing** when planning methods. ref-bio YAML is a routing index only — plans must open and read the actual upstream docs, or record `failed_to_access`. Relevant `.opencode/skills/bio-*` skills must be inspected by task type. See "ref-bio Usage Rules" and "Bioinformatics Skills (bio-* skills) Usage Rules" below.
- **Route labels matter**: if raw Space Ranger/Visium files are unavailable and analysis uses author-provided Seurat RDS image slots/coordinates instead, label the work as RDS-based reproduction/approximation, not strict author-code reproduction.

### Data Paths (expected, unverified)
- Spatial input: `data/raw/GSE233363_official/` for current author-provided RDS objects; `data/raw/spatial/` is reserved for future raw Space Ranger/Visium files if they are found
- Processed spatial: `data/processed/spatial/` (gitignored, generated)
- Spatial figures: `figures/spatial/` (gitignored, generated)
- Spatial scripts: `R/spatial/` (tracked, `s##_` prefix)
- Optional Python spatial scripts: `python/spatial/` (tracked source only; no environment or generated outputs)

### Script Naming
- Spatial scripts use `s##_` prefix (e.g., `s01_inspect_objects.R`), not `0##_`
- Script pattern: load → inspect → record provenance/object structure → save lightweight summary → gc()
- Before writing any spatial script: load one object, record its structure, consult ref-bio for reference routing (follow URLs to actual docs, not just YAML entries), and inspect relevant bio-* skills

### Memory Constraints (macOS 16GB)
- Load one spatial object at a time
- Current spatial focus is DG and Hippo; do not load WholeBrain unless a separate user-approved plan says to do so
- Use lightweight summaries, not duplicate RDS copies
- Call `gc()` after large operations
- Memory management specifics deferred to official Seurat docs and actual object structure

### Current Spatial Reproduction Route
- Stages Spatial-01 through Spatial-05 use author-provided RDS objects from `data/raw/GSE233363_official/`.
- Raw Visium/Space Ranger files required by the author's STutility section in Script 8 are not currently available.
- Phase Spatial-05 was therefore implemented as RDS-based reproduction/approximation using the inspected `seurat_Visium_Hippo_All.rds` image slots and coordinates.
- Do not install or require STutility for the current RDS-based Phase Spatial-05 unless raw Visium files are later found and the user explicitly approves reopening a strict raw-file branch.
- Phase Spatial-04 completed with WARNING; use of `metadata_hippo_func_region.csv` for Phase Spatial-05 requires explicit user acknowledgement.

### Current Scientific Trajectory
- First complete and validate the paper's hippocampus-related spatial workflow using DG/Hippo RDS objects and author Scripts 6-8 as references.
- Treat Phase Spatial-04/05 inflammatory-gradient outputs as partial/RDS-approximated unless validation proves strict agreement.
- The next reproduction priority is a hippocampal regional atlas stage for CA1, CA2, CA3, ML, GCL, Hilus, and derived DG; do not start mitochondrial analysis before this regional framework is validated.
- After hippocampus reproduction is validated, migrate the same region-aware framework to the user's CA1/CA3/DG mitochondrial-gene question.
- Phase Spatial-06: Hippocampal regional atlas completed (CA1=4671, CA2=565, CA3=2321, DG=2364).
- Phase Spatial-07: Reproduction validation completed (GO_WITH_CAVEATS).
- Phase Spatial-09: Target gene analysis plan completed (251 genes, 13 categories).
- Phase Spatial-10: Target gene audit and region-aware expression summaries completed (231 exact-match, 20 missing).
- Phase Spatial-11: CA1 vs CA3 pseudobulk DESeq2, module scores, and coupling analysis completed.
- Phase Spatial-11b: Age-grouping stratified DE and cross-age pattern classification completed (Young/Middle/Old).
- Phase Spatial-12: **Planned but NOT executed**. Young vs Aged regional comparison across CA1 / CA3 / DG. See `docs/spatial_phase12_young_aged_region_comparison_plan.md`.
- Phase Spatial-13: All-age candidate regulator-associated gene discovery completed (1,198 regulators x 11 modules, composite regulator_score, 6 candidate classes).
- Phase Spatial-13b: CA3 mirror display figures and top-candidate CSVs completed.
- Phase Spatial-13c: Age-stratified candidate regulator discovery completed (4 strata: CA1Y n=4, CA1O n=8, CA3Y n=3, CA3O n=8; 16 candidate classes; 9 figures; 29/29 validation PASS).
- Phase Spatial-14: Candidate regulator integration report completed. Report bundle in `data/report/phase14_candidate_regulator_report/` (markdown + PDF). 11 summary tables, 18 figures, 1175 unique prioritized genes. CA3 Young n=3 caveat enforced. No causal claims. report-level only, no new bioinformatics analysis.
- Do not define mitochondrial gene lists, run mitochondrial DE, or interpret mitochondrial biology until a dedicated mitochondrial planning phase selects references, gene sources, assay/layer choices, and sample-level aggregation strategy.
- Report bundles use `data/report/` (gitignored, generated). This is separate from `reports/` (audit logs) and `docs/` (committed plan docs).

### Phase 12 Guardrails (pre-execution rules)

When Phase 12 is executed, the following rules apply:

- **Primary analysis**: Region-stratified Old vs Young pseudobulk DESeq2 within CA1 / CA3 / DG.
- **Model**: `~ Age` (two-group unpaired) within each region. Age and GEMgroup are colinear, making `~ GEMgroup + Age` rank-deficient.
- **log2FC direction**: positive = Old/Aged higher, negative = Young higher.
- **Replicate unit**: GEMgroup-level pseudobulk. Spots are NOT biological replicates.
- **DG definition**: ML + GCL + Hilus (2,364 spots).
- **CA3 Young n=3**: Must be labeled `primary_with_caution_young_n3` in all outputs.
- **Middle excluded**: Middle is excluded from Phase 12 primary contrast.
- **Complex V / Atp5 missing caveat**: Must be carried forward from Phase 11.
- **limma/edgeR/interaction models**: Sensitivity/exploratory only. Primary remains DESeq2 `~ Age`.
- **No spot-level inference**: All DESeq2 uses pseudobulk aggregation.
- **No WholeBrain**: Only `seurat_Visium_Hippo_All.rds` is loaded.

### ref-bio Usage Rules (authoritative reference routing)

**ref-bio is NOT a read-only YAML index.** The `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml` is a routing table only — it maps source IDs to URLs. A link-only entry is NOT proof of review.

When planning a method or checking documentation:
1. Search `references.link-only.yaml` or `indexes/topic-map.yaml` for the relevant source_id.
2. Follow the upstream URL to the official documentation.
3. Actually open and read the relevant sections of the official docs/vignettes/guides.
4. Record findings in a structured review table.

**Plan documents MUST record for each consulted source:**

| Field | Required | Description |
|-------|----------|-------------|
| source_title | Yes | e.g., "DESeq2 Official Vignette" |
| url | Yes | Full URL to the specific page consulted |
| reviewed_status | Yes | `reviewed` / `partially_reviewed` / `failed_to_access` / `not_reviewed` |
| specific_section | Yes | e.g., "Standard workflow", "Multi-factor designs" |
| design_decision_affected | Yes | Which decisions this source validated or changed |

**Failure handling**:
- If a URL returns 404 or is inaccessible: record `reviewed_status = failed_to_access`, record the failure reason, and do NOT claim it was reviewed.
- If a URL requires network access that is not available: record `failed_to_access_network_unavailable`.
- link-only YAML entries serve as routing/index only; they are NOT a substitute for actual document review.

**Reference routing locations**:
- Project-specific refs (paper, GitHub, GEO, Zenodo): see `docs/reference_sources.md`
- Authoritative reference catalog: `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml`
- ref-bio indexes: `.opencode/skills/ref-bio/reference-pack/indexes/`
- Seurat spatial vignettes: https://satijalab.org/seurat/
- OSTA (Orchestrating Spatial Transcriptomics Analysis): https://bioconductor.org/books/release/OSTA/

### Bioinformatics Skills (bio-* skills) Usage Rules

**There is NO single `.opencode/skills/bioskill/` directory.** Do not assume a unified "bioskill" folder exists.

The project's bioinformatics skills are a collection of directories under `.opencode/skills/bio-*`. These are skill-specific instruction files (SKILL.md) with optional `references/usage-guide.md` companion documents.

**Plan phase workflow for bio-* skills:**
1. Identify the task type (e.g., differential expression, spatial statistics, visualization).
2. List relevant bio-* skills by matching task to skill name. Key skills for this project include:
   - `bio-differential-expression-deseq2-basics` — DESeq2 design, contrasts, shrinkage
   - `bio-differential-expression-de-results` — result extraction, padj handling, filtering
   - `bio-differential-expression-de-visualization` — volcano, MA, PCA plots
   - `bio-spatial-transcriptomics-spatial-preprocessing` — Python/Squidpy (conceptual only for R pipeline)
   - `bio-spatial-transcriptomics-spatial-data-io` — Python/Squidpy (conceptual only)
   - `bio-spatial-transcriptomics-spatial-visualization` — Python/Squidpy (conceptual only)
   - `bio-spatial-transcriptomics-spatial-statistics` — Python/Squidpy (conceptual only)
   - `bio-workflows-spatial-pipeline` — Python/Squidpy (conceptual QC checkpoint philosophy)
   - `bio-data-visualization-volcano-and-ma-plots` — volcano plot conventions
   - `bio-data-visualization-heatmaps-clustering` — heatmap rendering
3. Read only the relevant SKILL.md files (NOT bulk-read all skills).
4. Read companion `references/usage-guide.md` if the skill references one.
5. Record in the plan document:
   - skill path
   - status: `inspected` / `missing` / `not_relevant`
   - key guidance used
   - impact on plan

**IMPORTANT**: Skills that are Python/Squidpy-based (spatial-*) provide only conceptual guidance for this R/Seurat project. They should NOT be treated as direct execution instructions. For example, `bio-spatial-transcriptomics-spatial-statistics` documents Python Squidpy functions — it cannot be used directly in R. Its conceptual framework (Moran's I, spatial autocorrelation) may still inform the "spots ≠ replicates" principle.

### Optional Python / Tangram Route
- Python is not required for current Spatial-01 through Spatial-07 RDS-based stages.
- Python may be used only for an optional Tangram/scanpy branch after Phase Spatial-07 determines that cross-modality mapping is required.
- Python source belongs in `python/spatial/`; do not create Python environments or install packages without an approved dependency plan.
- If Python is used, record Python version, package versions, CPU/GPU mode, input paths, output paths, and deviations from the author notebook.

### Future Object Inspection Protocol
When spatial RDS files become available, the first script (`s01_inspect_objects.R`) must:
1. Record file path and size
2. Record Seurat version used to create the object
3. Record `DefaultAssay()` and all assay names
4. Record spatial key (`SpatialKey()`), image names, coordinate ranges
5. Record all metadata columns and their types/levels
6. Record spot count and gene count
7. Record provenance (GEO accession, alignment method, author)
