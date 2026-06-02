# Spatial Transcriptomics Overview — ferroDG Extension

## Status

**SCAFFOLD ONLY** — no analysis scripts, gene sets, DE models, enrichment, deconvolution,
label transfer, or biological interpretation have been implemented.

## Unverified Data Source Expectations

- **Source**: JessbergerLab / Wu et al. / GSE233363
- **Expected**: Author spatial RDS objects (seurat_Visium_DG_All.rds, seurat_Visium_Hippo_All.rds, seurat_Visium_WholeBrain.rds)
- **Status**: NOT VERIFIED — file paths, sizes, assay structure, metadata columns,
  spatial keys, image slots, and anatomical labels are all unknown
- **Action required**: Run `s01_inspect_objects.R` (to be created) when files are available
- **Do NOT assume**: SCT assay, spatial key name, images slot, anatomical region labels,
  sample IDs, or age labels

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
- Differential expression analysis plan
- Enrichment analysis plan
- Gene set definitions
- Deconvolution strategy
- Label transfer strategy
- Biological interpretation
- Analysis plan of any kind
- Assumptions about data structure or metadata
