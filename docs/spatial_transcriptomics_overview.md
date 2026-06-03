# Spatial Transcriptomics Overview — ferroDG Extension

## Status

Spatial scaffold plus early reproduction stages are active. Stages Spatial-01 through
Spatial-04 use author-provided DG/Hippo Seurat RDS objects. Phase Spatial-05 is planned
as an RDS-based reproduction/approximation because the raw Space Ranger/Visium files
required by the author STutility code are not currently available.

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

## Current Phase Spatial-05 Route

- Use `seurat_Visium_Hippo_All.rds` as the input object.
- Use Phase Spatial-04 `metadata_hippo_func_region.csv` only after explicit user acknowledgement of its WARNING status.
- Derive IS/NNS/ENS-like neighborhood labels from RDS-contained image slots and coordinates where possible.
- Process coordinates per image/sample; do not merge all image coordinate systems into one plane.
- Record all deviations from the author STutility route.
- Label outputs as **RDS-based approximation**, not strict reproduction.

If raw Visium files are later found, create a separate strict raw-file/STutility plan
instead of silently changing the RDS-based route.

## Script Convention

- Spatial scripts use `s##_` prefix (e.g., `s01_inspect_objects.R`)
- Located in `R/spatial/`
- Pattern: load → inspect → record provenance/object structure → save lightweight summary → gc()

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
