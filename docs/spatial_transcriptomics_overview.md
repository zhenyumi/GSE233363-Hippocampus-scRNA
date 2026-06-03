# Spatial Transcriptomics Overview — ferroDG Extension

## Status

Spatial scaffold plus hippocampus/DG reproduction stages are active. Stages
Spatial-01 through Spatial-05 use author-provided DG/Hippo Seurat RDS objects.
Phase Spatial-05 is implemented as an RDS-based reproduction/approximation because
the raw Space Ranger/Visium files required by the author STutility code are not
currently available.

The current project trajectory is:

1. Reproduce and validate the paper's hippocampus-related spatial workflow first.
2. Use that validated regional framework as the basis for a later CA1/CA3/DG
   mitochondrial-gene analysis.
3. Do not begin mitochondrial gene selection, differential testing, or biological
   interpretation until the hippocampal regional atlas and reproduction validation
   are complete.

## Data Source Status

- **Source**: JessbergerLab / Wu et al. / GSE233363
- **Current canonical input**: `data/raw/GSE233363_official/`
- **Available and inspected**: `seurat_Visium_DG_All.rds` and `seurat_Visium_Hippo_All.rds`
- **Available but not loaded for now**: `seurat_Visium_WholeBrain.rds`
- **Not currently available**: raw Visium/Space Ranger directories required by the author STutility code
- **Do NOT assume**: raw-file availability, exact coordinate semantics beyond inspected image slots, or equivalence between RDS-based approximation and strict author-code reproduction

## Critical Distinctions

| Concept | scRNA-seq | Visium Spatial |
|---------|-----------|----------------|
| Resolution | Single cell | Spot (captures signal from multiple cells) |
| Biological replicates | Cells are independent | Spots share local microenvironment |
| Metadata | Cell-level | Spot-level + image-level |
| Assay structure | RNA (counts/data) | Spatial (counts/data + coordinates + images) |

## Object Inspection Protocol

When spatial RDS files become available, the first script must record:

1. **File path and size** (measure actual size on disk)
2. **Seurat version** used to create the object
3. **DefaultAssay()** and all assay names (e.g., "Spatial", "SCT")
4. **Spatial key** (SpatialKey()), image names, coordinate ranges
5. **All metadata columns** and their types/levels
6. **Spot count and gene count**
7. **Provenance** (GEO accession, alignment method, author)

DG and Hippo have completed this inspection in Stage Spatial-01/02. Future stages
should reuse those inspection summaries and re-check critical structure at runtime.

## Current Reproduction Route

- Use `seurat_Visium_Hippo_All.rds` as the main hippocampus input object.
- Use Phase Spatial-04 `metadata_hippo_func_region.csv` only after explicit user acknowledgement of its WARNING status.
- Derive IS/NNS/ENS-like neighborhood labels from RDS-contained image slots and coordinates where possible.
- Process coordinates per image/sample; do not merge all image coordinate systems into one plane.
- Record all deviations from the author STutility route.
- Label outputs as **RDS-based approximation**, not strict reproduction.
- Next prioritize a hippocampal regional atlas reproduction for CA1, CA2, CA3,
  ML, GCL, Hilus, and derived DG before downstream mitochondrial planning.

If raw Visium files are later found, create a separate strict raw-file/STutility plan
instead of silently changing the RDS-based route.

## Script Convention

- Spatial scripts use `s##_` prefix (e.g., `s01_inspect_objects.R`)
- Located in `R/spatial/`
- Pattern: load → inspect → record provenance/object structure → save lightweight summary → gc()
- Optional Python/Tangram source, if approved later, belongs in `python/spatial/`

## Repository Organization

- The canonical repository-layout reference is `docs/repository_structure.md`.
- `R/spatial/` is the active spatial code directory.
- `python/spatial/` is reserved for optional future Tangram/scanpy source only.
- `analysis/` is retired and should not contain scripts, scratch files, or generated outputs.
- Spatial tables and checkpoints go under `data/processed/spatial/`.
- Spatial figures go under `figures/spatial/`.
- Spatial reports/logs go under `reports/spatial/` when they are generated artifacts; curated planning and result summaries stay under `docs/`.
- Root-level plot artifacts such as `Rplots.pdf` must be removed and kept gitignored.

## Optional Python Boundary

- Python is not needed for current RDS-based Stages Spatial-01 through Spatial-07.
- Python may be used later only for the optional Phase Spatial-08 Tangram/scanpy route if validation shows it is needed.
- Do not create Python environments or install Python packages without an approved plan.
- If Python is used, record the exact Python version, package versions, CPU/GPU mode, environment path, and all deviations from the author notebook.

## Reference Routing

- Project-specific refs: `docs/reference_sources.md`
- Authoritative ref catalog: `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml`
- Seurat spatial vignettes: https://satijalab.org/seurat/
- OSTA: https://bioconductor.org/books/release/OSTA/

## This Document Does NOT Include

- Spatial autocorrelation analysis plan
- Spatial domain detection plan
- Enrichment analysis plan
- Deconvolution strategy
- Label transfer strategy
- Biological interpretation
- Assumptions about data structure or metadata
- Mitochondrial gene-list definitions or CA1/CA3/DG mitochondrial analysis choices
