# Spatial Transcriptomics Reproduction Plan

> **Project**: ferroDG — ferroptosis in mouse dentate gyrus neurogenesis  
> **Source paper**: Wu et al. (2025) *Nature Neuroscience* 28:415–430  
> **DOI**: https://doi.org/10.1038/s41593-024-01848-4  
> **GEO**: GSE233363  
> **Plan version**: 3.1 (2026-06-03)  
> **Status**: Stages Spatial-01 through Spatial-04 implemented for DG/Hippo. Phase Spatial-05 is now planned as an RDS-based reproduction/approximation because raw Space Ranger/Visium files required by the author STutility code are not currently available.

> **Figures output convention**: All spatial figures use `figures/spatial/` (not `figures/stage-2/spatial/`). Rationale: AGENTS.md defines `figures/spatial/` as the spatial figures directory; `figures/stage-2/` is reserved for the scRNA-seq pipeline. Spatial analysis is a distinct pipeline branch and gets its own top-level spatial directory.

> **Current spatial object scope**: Future implementation phases should focus on `seurat_Visium_DG_All.rds` and `seurat_Visium_Hippo_All.rds`. `seurat_Visium_WholeBrain.rds` is inventoried only and should not be loaded or analyzed for now because of memory risk and current project focus. WholeBrain work requires a separate user confirmation and plan update.

> **Phase Spatial-05 route decision**: The author's Script 8 uses external raw Visium/Space Ranger files with `STutility::InputFromTable()`, `LoadImages()`, and `RegionNeighbours()`. Those files are not available locally, while the author-provided Hippo/DG Seurat RDS objects already contain `VisiumV1` image slots and spot coordinates. Therefore Phase Spatial-05 should proceed as an **RDS-based reproduction/approximation** using the inspected Seurat objects. This route must be labeled as an approximation and must not be reported as a strict author-code reproduction of the STutility section. If raw Visium files are later found, restore a separate strict-author-code plan.

---

## Read-Only Audit Summary

### Repository Key Files / Scaffold Status

| Path | Status |
|------|--------|
| `R/spatial/README.md` | Exists; Stage Spatial-01 inspection script implemented |
| `docs/spatial_transcriptomics_overview.md` | Exists |
| `docs/data_manifest_spatial.md` | Exists |
| `docs/reference_sources.md` | Exists |
| `docs/environment_renv_notes.md` | Exists |
| `data/raw/GSE233363_official/` | Contains author-provided official RDS files; current spatial focus is DG and Hippo |
| `data/raw/spatial/` | Empty; retained as a future spatial-data convention, not the current canonical source |
| `data/processed/spatial/inspection/` | Stage Spatial-01 inspection outputs for DG and Hippo |
| `data/processed/neuro_processed.rds` | Exists (scRNA-seq, potential reference) |
| `data/processed/neuro_with_ferroptosis_scores.rds` | Exists (scRNA-seq with scores) |
| `docs/original_paper.pdf` | Exists |

### Author Spatial Scripts and Responsibilities

| Script | Path | Lines | Responsibility |
|--------|------|-------|----------------|
| Script 6 | `docs/original_code_from_paper/Scripts/R/6. Spatial-Data process.R` | 933 | Raw H5 → SCTransform → merge → cluster → annotate → subset Hippo/DG → Allen Atlas |
| Script 7 | `docs/original_code_from_paper/Scripts/R/7. Spatial-DESeq2.R` | 326 | Pseudobulk DESeq2 on Hippo_All.rds |
| Script 8 | `docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R` | 993 | IFNγ module score, Edge detection, STutility RegionNeighbours, gradient DESeq2. Lines 71-697 require external raw Visium files; current Phase Spatial-05 will approximate this section from the Hippo RDS coordinates. |
| Tangram | `docs/original_code_from_paper/Scripts/Python/notebook_tangram.ipynb` | — | Cell-type deconvolution to Visium spots |

### Paper PDF Spatial Figures / Methods Targets

**PDF successfully read** via `pdftotext` (4740 lines extracted).

#### Main Figures — Spatial Panels

| Figure | Panels | Content | Author Code Source |
|--------|--------|---------|-------------------|
| **Fig. 1** | f, g, h | Whole-brain ST UMAP (42,169 spots); spatial projection by age; hippocampal cell-type mapping via Seurat CCA | Script 6 |
| **Fig. 5** | d–h | CAS (cellular aging signature) confirmed in ST data; pairwise DEG comparisons | Script 7 |
| **Fig. 6** | e | Spatial distribution of inflammatory signatures peaked at SLM | Script 8 |
| **Fig. 7** | a, c–f, g–j | IFNγ response spatial feature plots; distance-based spot hierarchy (IS/NNS/ENS1/ENS2); pseudospatial alignment; gene modules; microniche analysis | Script 8 |

#### Extended Data Figures — Spatial Panels

| Figure | Panels | Content | Author Code Source |
|--------|--------|---------|-------------------|
| **Ext. Data Fig. 2** | a–e | DG spatial transcriptomics: UMAP by age, hierarchical clustering of 6 regions, expression patterns, Allen Atlas validation, Tangram cell-type mapping | Script 6 + Tangram |
| **Ext. Data Fig. 3** | h | Igfbpl1 spatial expression (neuroblast marker in ST) | Script 6 |
| **Ext. Data Fig. 6** | b | Deconvoluted ST data: Astro 1 vs Astro 2 regional distribution in ventro-lateral DG | Script 6 |
| **Ext. Data Fig. 7** | c, d | Pairwise DEG comparisons in ST dataset (whole DG) | Script 7 |
| **Ext. Data Fig. 9** | a–i | IFNγ signature proportion (young 2.4%, middle 2.5%, old 5.9%); IS vs NNS volcano/GO; spatial correlation of T cells/STAT1+/microglia; PCA and hierarchical clustering of aggregated spots; CAS scores | Script 8 |
| **Ext. Data Fig. 10** | a–d | SLM spots: manual selection, BBB gene expression, IFNγ/BBB correlation | Script 8 |

### Data Objects Needed

| Object | Source | Expected Path | Status |
|--------|--------|---------------|--------|
| `seurat_Visium_DG_All.rds` | Zenodo / author-provided | `data/raw/GSE233363_official/` | Available; inspected in Stage Spatial-01 |
| `seurat_Visium_Hippo_All.rds` | Zenodo / author-provided | `data/raw/GSE233363_official/` | Available; inspected in Stage Spatial-01 |
| `seurat_Visium_WholeBrain.rds` | Zenodo / author-provided | `data/raw/GSE233363_official/` | Available; inventoried only; do not load for now |
| Raw Visium H5 files (16 samples) | Not on Zenodo / not found locally | `data/raw/spatial/` | Not available; required only for strict STutility route |
| Space Ranger spatial files (`tissue_positions*`, hires image, scalefactors JSON) | Author code | `data/raw/spatial/` | Not available; required only for strict STutility route |
| Hippocampus coordinate CSVs | Author code | `data/raw/spatial/` | Not available |
| `neuro_processed.rds` | Stage-2 pipeline | `data/processed/` | Exists |
| `neuro_with_ferroptosis_scores.rds` | Stage-2 pipeline | `data/processed/` | Exists |

### Unknowns / Gaps

1. **Zenodo RDS internal structure**: DG and Hippo have now been inspected in Stages Spatial-01 and Spatial-02. Both contain `VisiumV1` image slots and spot coordinates. WholeBrain remains inventory-only.
2. **Author uses Seurat v4.3; project uses v5.5.0**: Assay naming (`Spatial` vs `Spatial.01`) and SCT slot behavior may differ
3. **STutility package and raw inputs**: Used in author Script 8 for `RegionNeighbours`, but the required external raw Visium files are not currently available. Current Phase Spatial-05 should not install/use STutility unless raw files are found.
4. **Raw Visium H5 + spatial files**: Not available locally; strict author-code STutility reproduction requires these files.
5. **Tangram**: Paper uses v1.0.4 + scanpy 1.10.1 + CUDA; requires separate Python environment (not planned here)

---

## Phase Spatial-01: Object Availability and First-Load Inspection

### Goal
Inventory all author-provided spatial RDS objects, then inspect DG and Hippo sequentially with lightweight summaries. No file copies of full Seurat objects. WholeBrain remains inventoried only and is not loaded for now.

### Input Files
- `data/raw/GSE233363_official/seurat_Visium_DG_All.rds`
- `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds`
- `data/raw/GSE233363_official/seurat_Visium_WholeBrain.rds` (inventory only; not loaded)

### Output Files / Directories
- `data/processed/spatial/inspection/object_inventory.csv` — spatial object availability and inspection status
- `data/processed/spatial/inspection/DG/` — lightweight DG inspection outputs
- `data/processed/spatial/inspection/Hippo/` — lightweight Hippo inspection outputs

### Scripts
- `R/spatial/s01_inspect_objects.R`

### Dependencies to Verify
- Seurat v5.5.0 (check `renv.lock` entry for `Seurat`)
- `readRDS()` works on file
- `DefaultAssay()`, `Assays()`, `Images()`, `SpatialKey()` available

### Validation Checks
- Object class is `Seurat`
- At least one assay named `Spatial` or `Spatial.01`
- At least one image slot with coordinates
- Metadata has columns: `Age`, `Region`, `GEMgroup` (or equivalents)

### Memory Cautions
- Load one object at a time
- Call `gc()` immediately after inspection
- Do NOT load all three RDS simultaneously
- Do NOT load WholeBrain in current phases; it remains inventory-only until a separate user confirmation and plan update

### Stop Conditions
- If DG or Hippo RDS fails to load: stop, record the error, and do not automatically try WholeBrain
- If no spatial image slot found: object may be pre-processed differently than expected
- If assay names differ from author code: record discrepancy, continue
- If WholeBrain is requested without explicit confirmation: stop

---

## Phase Spatial-02: Reproduce Author Spatial Preprocessing / Object-Derived Metadata

### Goal
Reproduce the metadata and clustering structure relevant to DG and Hippo as described in Script 6, using the author-provided RDS objects (not raw H5). Verify that age groups, regional annotations, and hippocampal/DG subsetting match the author's reported structure. WholeBrain is out of scope for now.

### Input Files
- `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds`
- `data/raw/GSE233363_official/seurat_Visium_DG_All.rds`
- `docs/original_code_from_paper/Scripts/R/6. Spatial-Data process.R` (reference only)

### Output Files / Directories
- `data/processed/spatial/hippo_metadata_summary.csv` — hippocampal subregion counts
- `data/processed/spatial/dg_metadata_summary.csv` — DG spot counts by age

### Scripts to Create Later
- `R/spatial/s02_dg_hippo_metadata.R`

### Dependencies to Verify
- Seurat v5.5.0
- `SCTransform` (check if installed via renv.lock)
- `ggplot2` for basic QC plots

### Validation Checks
- DG and Hippo spot counts per age group match Stage Spatial-01 inspection and author-reported object structure
- DG regions present: ML, GCL, Hilus
- Hippo subregions present: ML, GCL, Hilus, CA1, CA2, CA3
- QC metrics: percentMito (`^mt-`), percentRibo (`^Rps|^Rpl`)

### Memory Cautions
- Load DG and Hippo sequentially, not together
- Do not load WholeBrain in Phase Spatial-02
- After object inspection and official Seurat documentation review, consider memory-light object handling if needed

### Stop Conditions
- If metadata columns don't match expected names: record all column names, continue with available columns
- If spot counts differ significantly (>10%): record discrepancy, check if Middle timepoint included
- If a requested action requires WholeBrain: stop and request a separate WholeBrain-specific confirmation and plan update

---

## Phase Spatial-03: Reproduce Spatial Differential / Regional Analyses

### Goal
Reproduce the pseudobulk DESeq2 analysis from Script 7, comparing Old vs Young (and optionally Middle vs Young) in hippocampal subregions.

### Input Files
- `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds`
- `docs/original_code_from_paper/Scripts/R/7. Spatial-DESeq2.R` (reference only)

### Output Files / Directories
- `data/processed/spatial/pseudobulk_counts_hippo.rds` — aggregated counts by Region×GEMgroup
- `data/processed/spatial/deseq2_OY_results.csv` — Old vs Young DEGs
- `data/processed/spatial/deseq2_MY_results.csv` — Middle vs Young DEGs (if Middle included)
- `figures/spatial/` — PCA, volcano, heatmap plots

### Scripts to Create Later
- `R/spatial/s03_pseudobulk_deseq2.R`

### Dependencies to Verify
- `DESeq2` (check renv.lock)
- `aggregate.Matrix` from `Matrix` package
- `ggplot2`, `pheatmap` for visualization

### Validation Checks
- PCA shows age-group separation (compare against paper Fig. 1f panels, author Script 7 outputs)
- Number of DEGs (adj P < 0.05) comparable to author's reported values in Extended Data Fig. 7c,d
- Top DEGs compared against author Script 7 output gene lists and paper Extended Data Fig. 7c,d panels

### Memory Cautions
- Pseudobulk aggregation reduces memory footprint significantly
- DESeq2 object is small; no special memory concerns

### Stop Conditions
- If `aggregate.Matrix` fails: check metadata columns for Region and GEMgroup
- If DESeq2 convergence issues: check sample sizes per group

---

## Phase Spatial-04: Reproduce IFNγ Module Score and Edge Detection

### Goal
Reproduce the IFNγ response module score and inflammatory spot (IS) identification from Script 8, using the Hippo_All.rds object.

### Input Files
- `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds`
- `docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R` (reference only)
- `msigdbr` package (for HALLMARK_INTERFERON_GAMMA_RESPONSE gene set)

### Output Files / Directories
- `data/processed/spatial/ifng_module_scores.csv` — per-spot IFNγ module scores
- `data/processed/spatial/is_nns_ens_labels.csv` — IS/NNS/ENS1/ENS2 labels per spot
- `figures/spatial/ifng_spatial_featureplot.pdf` — spatial feature plot of IFNγ response

### Scripts to Create Later
- `R/spatial/s04_ifng_module_and_edge_detection.R`

### Dependencies to Verify
- `msigdbr` (check renv.lock)
- `Seurat::AddModuleScore()`
- `ggplot2` for spatial feature plots

### Validation Checks
- IFNγ module scores compared against paper Extended Data Fig. 9a reported proportions (young 2.4%, middle 2.5%, old 5.9%)
- IS spatial distribution compared against paper Fig. 7a panels and author Script 8 outputs
- IS proportion in old hippocampus compared against paper Extended Data Fig. 9a panel

### Memory Cautions
- AddModuleScore operates on existing assay; no new large objects
- gc() after module score computation

### Stop Conditions
- If `msigdbr` not installed: record as dependency, defer to user approval
- If the gene set name is not found in `msigdbr`: first cross-check the exact gene set name against author Script 8, paper Methods, and ref-bio (`HALLMARK_INTERFERON_GAMMA_RESPONSE`). If the name still cannot be resolved after these checks: **stop**, document the mismatch, and request user decision. Do not substitute an alternative gene set.

---

## Phase Spatial-05: RDS-Based Inflammatory Gradient Approximation

### Goal
Approximate the inflammatory gradient neighborhood analysis from Script 8 using the author-provided `seurat_Visium_Hippo_All.rds`, its inspected image slots/coordinates, and the Phase Spatial-04 `FuncRegion` handoff. This phase should reproduce as much of the author/paper result as possible from RDS-contained spatial data, while explicitly recording deviations from the strict STutility/raw-file route.

This phase is **not** a strict author-code reproduction of Script 8 lines 71-697 because the required raw Visium/Space Ranger inputs are unavailable. Any output from this phase must be labeled **RDS-based approximation**.

### Input Files
- `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds`
- `data/processed/spatial/phase04_ifng_module/metadata_hippo_func_region.csv` (from Phase Spatial-04; current Phase 04 status is WARNING and requires explicit user acknowledgement before use)
- `docs/original_code_from_paper/Scripts/R/8. Spatial-Inflammatory gradient.R` (reference only)
- `data/processed/spatial/inspection/Hippo/spatial_info.csv` and Phase Spatial-02 metadata summaries (for validated image/coordinate expectations)

### Output Files / Directories
- `data/processed/spatial/phase05_rds_gradient/` - lightweight metadata, neighbor labels, validation tables, and provenance
- `figures/spatial/phase05/` - diagnostic plots and paper-comparison figures

Exact file names must be defined in the Phase Spatial-05 pre-execution plan after reviewing the inspected RDS coordinate structure and official Seurat spatial documentation.

### Scripts to Create Later
- `R/spatial/s05_rds_inflammatory_gradient.R`

### Dependencies to Verify
- `Seurat` - load Hippo RDS, access metadata, images, and coordinates
- `DESeq2` and `Matrix.utils` - only if the RDS-based phase includes the Script 8 gradient pseudobulk DESeq2 section
- `ggplot2` / existing plotting packages - for diagnostic and paper-comparison figures
- `STutility` - **not required for the current RDS-based route**. Do not install it for Phase Spatial-05 unless raw Visium files are later found and a strict STutility branch is explicitly approved.
- `Monocle 2` - **not required for the current RDS-based route** unless a separate paper-methods branch is explicitly approved.

### Validation Checks
- Confirm Phase Spatial-04 WARNING handoff has explicit user acknowledgement before consuming `metadata_hippo_func_region.csv`
- Confirm the Hippo RDS still has 16 images and 9921 spots before neighborhood approximation
- Confirm coordinate extraction is performed per image/sample, not across all images as one shared coordinate plane
- Compare RDS-derived IS/NNS/ENS-like counts against author Script 8 comments where available
- Compare resulting patterns against paper Fig. 7 and Extended Data Fig. 9 as qualitative/quantitative approximation, not strict equivalence
- Record all deviations from author Script 8, especially lack of external raw image files and lack of STutility `RegionNeighbours()`

### Memory Cautions
- Load only `seurat_Visium_Hippo_All.rds`; do not load DG, WholeBrain, or Chromium in this phase unless a later plan explicitly requires it
- Work per image/sample where possible
- Save lightweight metadata/CSV outputs; do not save modified full Seurat objects
- Preserve sparse matrices and avoid dense conversion of spot-level count matrices

### Stop Conditions
- If Phase Spatial-04 handoff status is FAIL: stop; do not run Phase Spatial-05.
- If Phase Spatial-04 handoff status is WARNING and user acknowledgement is not documented: stop and request confirmation.
- If Hippo RDS lacks image slots or coordinates at runtime: stop and report that RDS-based approximation cannot proceed.
- If neighborhood labels cannot be derived per image without guessing coordinate systems or sample IDs: stop and request user decision.
- If the implementation starts to require raw Visium files or STutility behavior that cannot be reproduced from the RDS: stop and document the boundary rather than silently substituting.

---

## Phase Spatial-06: Optional Tangram / Cross-Modality Branch

### Goal
Reproduce the Tangram cell-type deconvolution from the Python notebook, projecting scRNA-seq cell types onto Visium spots. **This phase is optional** — only proceed if paper figure mapping confirms it is required for target figures.

### Input Files
- `data/processed/neuro_processed.rds` (scRNA-seq reference)
- `data/raw/GSE233363_official/seurat_Visium_DG_All.rds` (or `seurat_Visium_Hippo_All.rds`)
- `docs/original_code_from_paper/Scripts/Python/notebook_tangram.ipynb` (reference only)

### Output Files / Directories
- `data/processed/spatial/tangram_mapped_annotations.csv` — per-spot cell type probabilities
- `figures/spatial/tangram_celltype_overlay.pdf` — spatial overlay of mapped cell types

### Scripts to Create Later
- `R/spatial/s06_tangram_celltype_mapping.R` (or Python equivalent)

### Dependencies to Verify
- Python environment with: `tangram==1.0.4`, `scanpy==1.10.1`, `squidpy`
- CUDA GPU (paper specifies `cells` mode, 1000 epochs, `rna_count_based` density prior)
- If no CUDA: CPU fallback exists but will be slower; document runtime differences

### Validation Checks
- Mapped cell types compared against paper Extended Data Fig. 2e panels and author Tangram notebook outputs
- Per-region cell type distributions compared against paper Extended Data Fig. 2e and author Script 6 outputs
- Only lightweight final summaries saved (not intermediate AnnData files)

### Memory Cautions
- Python memory separate from R
- Do not load full AnnData in R; use CSV exchange

### Stop Conditions
- If CUDA not available: record as implementation constraint, defer to user decision
- If Python environment not set up: record as prerequisite, do not install in this plan
- If paper figures do not require Tangram output: skip this phase entirely

### Implementation Constraints (to log, not resolve now)
- CUDA vs CPU differences in runtime and convergence
- Tangram v1.0.4 vs current version compatibility
- scanpy version compatibility with tangram

---

## Phase Spatial-07: Figure/Output Validation and Documentation

### Goal
Validate all reproduced figures against paper targets, document discrepancies, and produce a final reproducibility report.

### Input Files
- All outputs from Phases Spatial-01 through Spatial-06
- `docs/original_paper.pdf` (for figure comparison)

### Output Files / Directories
- `docs/spatial_figure_validation.md` — per-figure validation status
- `docs/spatial_reproducibility_report.md` — final report with known limitations
- `figures/spatial/` — all reproduced figures

### Scripts to Create Later
- `R/spatial/s07_figure_validation.R`

### Validation Checklist
| Paper Figure | Panels | Reproduction Status | Notes |
|-------------|--------|-------------------|-------|
| Fig. 1f-h | Whole-brain UMAP, spatial projection, cell-type mapping | Pending | Depends on Phase 02 |
| Fig. 5d-h | CAS in ST | Pending | Depends on Phase 03 |
| Fig. 7a | IFNγ spatial feature plots | Pending | Depends on Phase 04 |
| Fig. 7c-f | Spatial hierarchy, pseudotime, modules | Pending | Depends on Phase 05 |
| Ext. Data Fig. 2a-e | DG spatial transcriptomics | Pending | Depends on Phase 02 |
| Ext. Data Fig. 7c,d | Pairwise DEG in ST | Pending | Depends on Phase 03 |
| Ext. Data Fig. 9a-i | IFNγ gradient analysis | Pending | Depends on Phase 04-05 |

### Documentation Requirements
- Record all discrepancies between reproduced and original figures
- Note where author code uses tools not mentioned in paper (e.g., STutility)
- Note Seurat version differences (v4.3 in paper vs v5.5.0 in project)
- Note missing raw data (H5 files, CSVs) that prevent full reprocessing

---

## Do Not Do Yet

The following actions are explicitly excluded from this plan:

1. **Do not install any packages** — no `renv::install()`, no `install.packages()`, no `BiocManager::install()`
2. **Do not modify `renv.lock`** — no changes to dependency versions
3. **Do not download data** — no Zenodo downloads, no GEO downloads, no SRA downloads
4. **Do not run R scripts** — no `Rscript` execution, no interactive R sessions
5. **Do not create analysis scripts** — no `R/spatial/s01*.R` etc. until implementation phase
6. **Do not set up Python environment** — no `pip install`, no conda, no virtualenv
7. **Do not modify existing scRNA-seq pipeline** — no changes to `R/01_*.R` through `R/05_*.R`
8. **Do not save full Seurat object duplicates** — only lightweight CSV/text summaries
9. **Do not perform biological reinterpretation** — reproduce author methods, do not add new analyses
10. **Do not create `docs/spatial/` directory** — defer to implementation phase for validation/report files

---

## Execution Order

```
Phase 01 (Inspect)
    ↓
Phase 02 (QC/Metadata)  ← depends on 01
    ↓
Phase 03 (DESeq2)       ← depends on 02
    ↓
Phase 04 (IFNγ Module)  ← depends on 02
    ↓
Phase 05 (Gradient)     ← depends on 04
    ↓
Phase 06 (Tangram)      ← optional, depends on 02 + Python env
    ↓
Phase 07 (Validation)   ← depends on all above
```

Phases 03 and 04 can run in parallel after Phase 02 completes.

---

## Ref-Bio Source Routing

> **`.opencode/skills/ref-bio/`** is the authoritative reference registry. The table below is an execution-phase short routing summary for quick lookup during implementation. For full reference details, provenance, and upstream URLs, consult `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml` and the index files in `.opencode/skills/ref-bio/reference-pack/`.

| Source ID | Package/Tool | URL | Priority |
|-----------|-------------|-----|----------|
| `seurat-spatial` | Seurat Spatial | https://satijalab.org/seurat/articles/spatial_vignette.html | must_read |
| `osta` | OSTA | https://bioconductor.org/books/release/OSTA/ | must_read |
| `seurat` | Seurat v5 | https://satijalab.org/seurat/ | must_read |
| `deseq2` | DESeq2 | https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html | must_read |
| `stutility` | STutility | https://github.com/jbergenstrahle/STutility | strict_raw_route_only |
| `visium` | 10x Visium | https://www.10xgenomics.com/products/spatial-gene-expression | optional_read |
| `spatialexperiment` | SpatialExperiment | https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html | optional_read |
| `tangram` | Tangram | https://github.com/broadinstitute/Tangram | optional_read |

**Note**: STutility reference URL is from the author's GitHub. The author Script 8 STutility section requires external raw Visium files that are not currently available. Current Phase Spatial-05 is an RDS-based approximation and should consult STutility only to document the strict-route discrepancy, not to implement the current route.

---

## Open Questions

1. **Data download timing**: When should Zenodo RDS files be downloaded? Before Phase 01 or as needed?
2. **Tangram requirement**: Is Tangram output needed for any target figure? If not, Phase 06 can be skipped.
3. **Middle timepoint**: Paper includes Middle-aged group. Author code compares OY, MY, OM. Should Middle be included in reproduction?
4. **Full reprocessing vs Zenodo RDS**: Current decision is to proceed with RDS-based approximation because raw Visium/Space Ranger files are not available. Revisit only if those files are found.
5. **STutility dependency**: Current decision is not to install/use STutility for Phase Spatial-05. A strict STutility branch can be reopened only if raw Visium files are found.

---

## Package Dependencies Summary

> **Note**: All package statuses listed below are **"To verify"** — check `renv.lock` and `renv::status()` before implementation. "Installed" is not confirmed until renv.lock entry is inspected.

| Package | Required By | Status | Phase |
|---------|------------|--------|-------|
| Seurat (v5.5.0) | All phases | To verify (check renv.lock) | All |
| SCTransform | Phase 02 | To verify | 02 |
| DESeq2 | Phases 03, 05 | To verify | 03, 05 |
| msigdbr | Phase 04 | To verify | 04 |
| STutility | Strict raw-file branch only | Not required for current RDS-based Phase 05 | 05 optional/future |
| Monocle 2 | Paper-methods branch only | Not required for current RDS-based Phase 05 | 05 optional/future |
| ggplot2 | All visualization | To verify (check renv.lock) | All |
| pheatmap | Heatmaps | To verify | 03 |
| Matrix | aggregate.Matrix | To verify (check renv.lock) | 03 |
