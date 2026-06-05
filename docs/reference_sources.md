# Reference Sources — ferroDG Spatial Extension

This document is a **reference routing guide**, not an exhaustive bibliography.
For authoritative reference routing, consult `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml`.

## How to Use This Registry

This file is a **registry and routing list** — it maps project-specific references and tool names to source IDs and upstream URLs. It is NOT a record of what has been reviewed.

**A source listed here is NOT "reviewed" unless it was actually opened and read during a planning or execution task.** To claim a source was reviewed:

1. Follow the URL to the official documentation.
2. Read the relevant sections.
3. Record in the plan document:
   - **source_title** (e.g., "DESeq2 Official Vignette")
   - **url** (full URL to the specific page consulted)
   - **reviewed_status**: `reviewed` / `partially_reviewed` / `failed_to_access` / `not_reviewed`
   - **specific_section** (e.g., "Standard workflow", "Multi-factor designs")
   - **design_decision_affected** (which decisions this source validated or changed)
4. If a URL is inaccessible: record `reviewed_status = failed_to_access` and the reason. Do NOT claim it was reviewed.
5. Link-only entries and this registry are routing aids — they do not substitute for actual source review.

See `AGENTS.md` "ref-bio Usage Rules" and `docs/agent_reference_workflow.md` for the full reference-first planning workflow.

## Project-Specific References

| Reference | Source | Notes |
|-----------|--------|-------|
| Wu et al. (original paper) | Published article | ferroDG methodology, gene sets, figures |
| ferroDG GitHub | https://github.com/ShilongZhang116/ferroDG | Stage-2 analysis code |
| JessbergerLab GitHub | https://github.com/JessbergerLab/AgingNeurogenesis_Transcriptomics | Author-provided analysis code |
| GSE233363 (GEO) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363 | scRNA-seq + Visium data |
| Zenodo RDS | Author-provided | Official Seurat objects |

## Reference Routing via ref-bio

The `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml` catalog provides
link-only references to upstream official sources. No source content has been acquired.

When planning methods or checking tool documentation:
1. Search `references.link-only.yaml` for the relevant source_id
2. Follow the upstream URL to the official documentation
3. Record the source_id used in analysis scripts

## Tool Routing Table

| Tool/Method | ref-bio source_id | Upstream URL |
|-------------|-------------------|--------------|
| Seurat (spatial vignettes) | seurat | https://satijalab.org/seurat/ |
| SpatialExperiment | spatialexperiment | https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html |
| SpatialFeatureExperiment | spatialfeatureexperiment | https://bioconductor.org/packages/release/bioc/html/SpatialFeatureExperiment.html |
| DESeq2 | deseq2 | https://bioconductor.org/packages/release/bioc/html/DESeq2.html |
| clusterProfiler | clusterprofiler | https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html |
| OSTA (guide) | osta | https://bioconductor.org/books/release/OSTA/ |
| Visium (10x platform) | visium | https://www.10xgenomics.com/platforms/visium |
| STutility | stutility | https://github.com/ludvigla/STUtility |
| BayesSpace | bayesspace | https://bioconductor.org/packages/release/bioc/html/BayesSpace.html |
| SPARK-X | spark-x | https://github.com/xzhoulab/SPARK |
| nnSVG | nnsvg | https://github.com/lmweber/nnSVG |
| Squidpy | squidpy | https://github.com/scverse/squidpy |
| Giotto | giotto | https://github.com/drieslab/Giotto |

## Scope Notes

- **Tool/method choices**: Must go through ref-bio reference routing and actual object inspection first; not pre-decided here
- **Visium-only scope**: No MERFISH, Slide-seq, Xenium, CosMx, or HD placeholders
- **ref-bio is reference routing**: Not a method summary; links to upstream official sources
- **No bundled content**: All references are link-only; no source content has been acquired
