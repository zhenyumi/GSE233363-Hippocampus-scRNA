# Agent Reference-First Planning Workflow

> **Purpose**: Step-by-step guide for using ref-bio and bio-* skills when planning analysis phases. Avoids overloading AGENTS.md with workflow details. Linked from AGENTS.md, README.md, and reference_sources.md.

## Step-by-Step Reference-First Workflow

### 1. Identify Task Type

Classify the task by:
- **Data modality**: spatial transcriptomics (Visium), scRNA-seq, bulk RNA-seq
- **Workflow stage**: differential expression, module score, visualization, coupling
- **Package/tool**: DESeq2, Seurat, ggplot2, etc.
- **Platform**: Visium RDS-based, raw Space Ranger, etc.

### 2. Query ref-bio for Authoritative Sources

```bash
# List available source IDs
ls .opencode/skills/ref-bio/reference-pack/indexes/

# Search for relevant source IDs by topic
rg "deseq2|spatial|pseudobulk|seurat" .opencode/skills/ref-bio/reference-pack/references.link-only.yaml

# Look up by topic category
rg "differential.expression|spatial" .opencode/skills/ref-bio/reference-pack/indexes/topic-map.yaml
```

### 3. Inspect Relevant bio-* Skills

List skills by task category:

| Task | Relevant Skills |
|------|----------------|
| Differential expression (DESeq2) | `bio-differential-expression-deseq2-basics`, `bio-differential-expression-de-results` |
| DE visualization | `bio-differential-expression-de-visualization`, `bio-data-visualization-volcano-and-ma-plots` |
| Spatial preprocessing/QC | `bio-spatial-transcriptomics-spatial-preprocessing` (Python/Squidpy — conceptual only) |
| Spatial statistics | `bio-spatial-transcriptomics-spatial-statistics` (Python/Squidpy — conceptual only) |
| Heatmaps/clustering | `bio-data-visualization-heatmaps-clustering` |
| Spatial pipeline overview | `bio-workflows-spatial-pipeline` (conceptual QC checkpoints) |

Read each relevant skill's `SKILL.md` and any companion `references/usage-guide.md`.

### 4. Open and Review Authoritative Documentation

For each source identified in step 2:

1. Open the URL in a browser (or use webfetch).
2. Find the specific section relevant to the design decision.
3. Extract key guidance, caveats, and recommendations.
4. Record in the plan document.

### 5. Record Review Status in Plan Document

Use this table format in every plan document:

| # | Source Title | URL | Reviewed Status | Specific Section Consulted | Design Decision Affected |
|---|-------------|-----|-----------------|---------------------------|------------------------|
| R1 | DESeq2 Official Vignette | https://bioconductor.org/.../DESeq2.html | reviewed | "Standard workflow", "Multi-factor designs" | Confirmed `~ Age` two-group design |
| R2 | Seurat Spatial Vignette | https://satijalab.org/.../spatial_vignette.html | reviewed | "Working with multiple slices" | Confirmed pseudobulk by GEMgroup |
| R3 | OSTA | https://bioconductor.org/books/release/OSTA/ | partially_reviewed | TOC inspected; Ch.33 identified | Validated spatial DE as recognized pattern |

And for bio-* skills:

| # | Skill Path | Status | Key Guidance Used | Impact on Plan |
|---|-----------|--------|------------------|----------------|
| 1 | `.opencode/skills/bio-differential-expression-deseq2-basics/SKILL.md` | inspected | Two-group `~ condition` design; apeglm shrinkage; pre-filtering | Confirmed Phase 12 DESeq2 design |
| 2 | `.opencode/skills/bio-spatial-transcriptomics-spatial-statistics/SKILL.md` | inspected (Python/squidpy) | Conceptual: spatial autocorrelation, "spots ≠ replicates" | Not directly applicable but informs replicate unit |
| 3 | `.opencode/skills/bio-differential-expression-de-results/SKILL.md` | inspected | padj=NA diagnostic (independent filtering / Cook's / all-zero), TREAT vs post-hoc LFC filtering caveat, p-value histogram QC, small-n replication reality (Schurch 2016), IHW for +power | Confirmed Phase 12 result handling: padj<0.05 threshold, post-hoc |log2FC|>log2(1.5) as notable (not FDR-controlled), p-value histogram as QC |

### 6. Finalize Plan

Only after steps 1-5 are complete:
- Lock the design decisions with confidence levels.
- Record any sources that failed to access.
- Record any skills that are Python/Squidpy and provide conceptual guidance only.
- Do NOT claim a source was reviewed based solely on a link-only YAML entry.

## Failure Handling

| Scenario | Action |
|----------|--------|
| URL returns 404 | Record `reviewed_status = failed_to_access` with reason; do NOT claim reviewed |
| URL requires network unavailable | Record `failed_to_access_network_unavailable` |
| ref-bio source_id not found | Note in plan; search official docs directly by package name |
| bio-* skill is Python/Squidpy but task is R/Seurat | Record as `inspected (Python/tool)` with note "conceptual guidance only" |
| bio-* skill not found for needed topic | Record as `missing`; proceed with direct authoritative source review |

## Example: Phase 12 Review Table

This is how the Phase 12 plan (`docs/spatial_phase12_young_aged_region_comparison_plan.md` §0) records its review:

**Authoritative Reference Review**:

| # | Source Title | URL | Reviewed Status | Specific Section | Design Decision Affected |
|---|-------------|-----|-----------------|------------------|------------------------|
| R1 | DESeq2 Official Vignette | https://bioconductor.org/.../DESeq2.html | reviewed | "Standard workflow", "Note on factor levels" | Confirmed `~ Age` two-group design |
| R2 | Seurat Spatial Vignette | https://satijalab.org/.../spatial_vignette.html | reviewed | "Working with multiple slices" | Confirmed pseudobulk by GEMgroup |
| R3 | OSTA | https://bioconductor.org/books/release/OSTA/ | partially_reviewed | TOC inspected | Validated spatial DE as recognized pattern |
| R4 | Seurat DE Testing | https://satijalab.org/.../differential_expression.html | failed_to_access | HTTP 404 | No design impact; DESeq2 vignette is primary |

**Bio-Skill Guidance Audit**:

| # | Skill Path | Status | Key Guidance | Impact |
|---|-----------|--------|-------------|--------|
| 1 | `bio-differential-expression-deseq2-basics` | inspected | Two-group design, apeglm, pre-filtering | Confirmed Phase 12 design |
| 2 | `bio-differential-expression-de-results` | inspected | padj=NA diagnostic, TREAT vs post-hoc | Confirmed significance thresholds |
| 3 | `bio-spatial-transcriptomics-spatial-statistics` | inspected (Python) | Spatial autocorrelation framework | Conceptual: reinforces "spots ≠ replicates" |
