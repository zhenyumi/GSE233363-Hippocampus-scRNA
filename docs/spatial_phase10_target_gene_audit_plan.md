# Phase Spatial-10: Target-Gene Audit + Region-Aware Expression Summary

> **Plan version**: v2.2 (2026-06-04)
> **Status**: APPROVED — Phase 10 execution
> **Phase 09 docs**: `docs/spatial_phase09_mito_target_gene_analysis_plan.md` (v2.1), `docs/spatial_mito_target_gene_strategy.md` (v2.1)
> **Handoff**: GO_WITH_CAVEATS | Tangram: TANGRAM_DEFERRED

---

## 1. Input Summary

### 1.1 Spatial Object

| Field | Value |
|-------|-------|
| Object | `data/raw/GSE233363_official/seurat_Visium_Hippo_All.rds` |
| Spots | 9,921 |
| Images | 16 |
| Assays | Spatial (32,285 features), SCT (17,096 features) |
| Regions | CA1=4671, CA2=565, CA3=2321, ML=1480, GCL=645, Hilus=239 |
| DG (derived) | 2,364 (ML+GCL+Hilus) |

### 1.2 Gene List — `docs/target_genes.csv`

| Property | Value |
|----------|-------|
| Path | `docs/target_genes.csv` |
| Encoding | GBK (Chinese Windows) |
| Columns | `分类亚类` (category), `基因名` (gene names) |
| Gene delimiter | `、` (U+3001, CJK ideographic comma) |
| Total unique genes | 251 |
| Duplicates | 0 |
| Categories | 13 (see §2.2) |
| Case distribution | 242 TitleCase, 9 mixed (mt-Co1 etc.), 0 UPPERCASE |

### 1.3 Category English Mapping

| Chinese Category | English Reference |
|-----------------|-------------------|
| 线粒体基因组编码（mt-） | mtDNA-encoded |
| 呼吸链复合体 I（Nduf/Ndufv/Ndufb/Ndufc） | Complex I |
| 呼吸链复合体 III（Uqcr/Uqcrc/Uqcc） | Complex III |
| 呼吸链复合体 IV（Cox） | Complex IV |
| 呼吸链复合体 V（Atp5 亚基，ATP 合酶） | Complex V |
| 三羧酸循环 TCA 关键酶 | TCA cycle |
| 线粒体核糖体蛋白（MRPL/MRPS，同时归入线粒体翻译） | Mitoribosomal |
| 线粒体外膜 / 内膜转运 / 孔蛋白 | Membrane/Translocator |
| 线粒体分子伴侣 / 抗氧化 / 代谢酶 | Defense/Antioxidant |
| RPL（胞质大亚基） | RPL |
| RPS（胞质小亚基） | RPS |
| 起始因子 Eif | Eif |
| 延伸因子 Eef | Eef |

---

## 2. CSV Parsing Strategy

### 2.1 Multi-Section Structure

The CSV has 3 sections separated by empty rows, each with its own header:

- **Section A** (lines 1-10): Mitochondrial-related categories (9 categories)
- **Section B** (lines 12-14): Ribosomal categories (RPL, RPS)
- **Section C** (lines 17-19): Translation factor categories (Eif, Eef)

**Parsing approach**:
1. Read raw lines with `readLines()` using GBK encoding
2. Skip empty lines and repeated header lines
3. Split each row on `,` to get category + gene string
4. Split gene string on `、` (U+3001) to get individual gene symbols
5. Trim whitespace from each gene symbol
6. Build unified data.frame with columns: `category`, `gene_symbol`

### 2.2 Encoding Handling

```r
raw_lines <- readLines("docs/target_genes.csv", encoding = "GBK")
# If GBK fails, try locale default
```

Record encoding used in provenance. Do NOT auto-convert file to UTF-8.

---

## 3. Gene Symbol Mapping

### 3.1 Expected Outcome

Since 242/251 genes are already TitleCase and 9 are mixed-case (mt-Co1 etc.), ~100% exact match in `rownames(hippo[["Spatial"]])` is expected. The UPPERCASE → TitleCase candidate path is unlikely to be needed.

### 3.2 Audit Columns

| Column | Source |
|--------|--------|
| `original_symbol` | Verbatim from CSV |
| `original_case` | TitleCase / Mixed / UPPERCASE / lowercase |
| `exact_match_spatial` | `original_symbol %in% rownames(hippo[["Spatial"]])` |
| `exact_match_sct` | `original_symbol %in% rownames(hippo[["SCT"]])` |
| `candidate_mouse_symbol` | If UPPERCASE + no exact match: `to_title_case_base()`. Else: `original_symbol` |
| `candidate_match_spatial` | `candidate_mouse_symbol %in% rownames(hippo[["Spatial"]])` |
| `candidate_match_sct` | `candidate_mouse_symbol %in% rownames(hippo[["SCT"]])` |
| `mapping_status` | `exact_match` / `candidate_match` / `not_found_in_either` |
| `manual_review` | TRUE if `not_found_in_either` |

### 3.3 Base R Helper

```r
to_title_case_base <- function(x) {
  x_lower <- tolower(x)
  substr(x_lower, 1, 1) <- toupper(substr(x_lower, 1, 1))
  return(x_lower)
}
```

No `stringr` dependency.

---

## 4. Region Group Construction

| Group | Definition | Spots |
|-------|-----------|-------|
| CA1 | `Region == "CA1"` | 4,671 |
| CA3 | `Region == "CA3"` | 2,321 |
| DG | `Region %in% c("ML", "GCL", "Hilus")` | 2,364 |
| CA2 | `Region == "CA2"` | 565 |

CA2 is included as context/QC in summaries, but not treated as a primary comparison region.

DG subregions (ML, GCL, Hilus) are preserved in subregion-level summaries.

---

## 5. Region-Aware Expression Summaries

### 5.1 Output Tables

| File | Granularity | Description |
|------|-------------|-------------|
| `region_sample_summary.csv` | gene × region_group × Age × GEMgroup | n_spots, pct_expressed, mean_count, sum_count |
| `region_age_summary.csv` | gene × region_group × Age | Collapsed from sample-level |
| `region_subregion_summary.csv` | gene × subregion × Age | ML/GCL/Hilus separate |
| `low_coverage_gemgroups.csv` | region_group × GEMgroup | Combinations with < 20 spots |

### 5.2 Extraction Method

```r
# Per-gene numeric vector extraction (NOT full dense matrix)
counts_mat <- GetAssayData(hippo, assay = "Spatial", layer = "counts")
# For each gene: gene_vec <- as.numeric(counts_mat[gene, ])
# Sparse indexing, no dense conversion of full matrix
```

**Prohibited**: Full count matrix dense conversion. Allow per-gene numeric vector extraction only. Record in provenance.

### 5.3 Spot-Level Long Table

`spot_expression_nonzero_long.csv.gz` — only nonzero counts are written, zero counts are omitted. Columns: cell_id, gene, count, Region, region_group, Age, GEMgroup, image.

pct_expressed and zero-aware summaries come from region summary tables.

---

## 6. Plot Design

### 6.1 QC Plots (Default: 3)

| Plot | Type | Description |
|------|------|-------------|
| Gene detection dot plot | ggplot2 | Genes (y) × region_group (x), dot size = pct_expressed, color = mean_count |
| Category bar plot | ggplot2 | Category (x) × n_genes_detected (y), filled by detection status |
| Expression heatmap | pheatmap | Top N genes by variance, clustered by gene and region_group |

### 6.2 Spatial Maps (Default: ≤12 genes)

- Select 3-6 genes by deterministic smoke-test (first 2-3 genes per top 4 categories, or CSV priority column if present)
- Default: 1 slice (first image), up to 3 genes per spatial map panel
- Coordinate extraction via `GetTissueCoordinates()`
- Plot with `ggplot2::geom_point()` + `scale_y_reverse()` (see §6.3)

### 6.3 Coordinate Fix

```r
coords <- GetTissueCoordinates(hippo, image = img_name)
# coords columns may be named imagecol/imagerow or x/y
if ("imagecol" %in% colnames(coords) && "imagerow" %in% colnames(coords)) {
  plot_x <- coords[["imagecol"]]
  plot_y <- coords[["imagerow"]]
  reverse_y <- TRUE
} else if ("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
  plot_x <- coords[["x"]]
  plot_y <- coords[["y"]]
  reverse_y <- FALSE
} else {
  stop("Unknown coordinate columns: ", paste(colnames(coords), collapse = ", "))
}
# Plot: ggplot(...) + geom_point() + if (reverse_y) scale_y_reverse()
```

### 6.4 Plot Limits

- Maximum 12 spatial map PNGs
- If pheatmap fails: fallback to `geom_tile()` heatmap
- Record which genes were plotted in validation summary
- No Rplots.pdf generation (guard: record pre-existence, remove only if generated during run)

---

## 7. Validation Checks

| ID | Check | Criterion |
|----|-------|-----------|
| V01 | CSV readable | ≥ 1 gene row parsed |
| V02 | Gene column identified | 251 unique genes extracted |
| V03 | Encoding handled | CSV loaded without encoding error |
| V04 | Gene presence audit complete | All 251 genes checked against Spatial rownames |
| V05 | Symbol mapping audit complete | All 251 genes have mapping_status |
| V06 | Region counts preserved | CA1=4671, CA3=2321, DG=2364 |
| V07 | No dense conversion | Provenance records per-gene extraction only |
| V08 | No Seurat object saved | No .rds written |
| V09 | GEMgroup × Age mapping | 16 GEMgroups verified |
| V10 | Missing genes documented | Count and names in missing_genes.csv (if any) |
| V11 | All found genes summarized | Every found gene in region_sample_summary.csv |
| V12 | Output files lightweight | No file > 50 MB |
| V13 | renv.lock unchanged | MD5 = `139bcfaf7ce24b633c9f034380a15ef5` |
| V14 | No WholeBrain loaded | Provenance documents single-object load |
| V15 | CA2 present as context | CA2 in summary tables |
| V16 | DG = ML+GCL+Hilus | Derived count = 2364 |
| V17 | spot_expression_nonzero_long.csv.gz | Only nonzero rows, compressed |
| V18 | Plotted gene limit | ≤ 12 genes in spatial maps |
| V19 | Rplots.pdf handling | Pre-existence recorded; removed only if generated during run |
| V20 | Category mapping complete | All 13 categories mapped to English reference |

---

## 8. Stop Conditions

| ID | Condition | Action |
|----|-----------|--------|
| SC-1 | Hippo RDS not found | STOP |
| SC-2 | CSV unreadable (encoding + parse failure) | STOP |
| SC-3 | 0 gene rows after parsing | STOP |
| SC-4 | Phase 06 validation not PASS | STOP |
| SC-5 | Region count mismatch > 5% | STOP |
| SC-6 | User requests WholeBrain loading | STOP — requires separate plan |
| SC-7 | User requests Tangram | STOP — requires separate plan |
| SC-8 | User requests DESeq2/enrichment | STOP — Phase 11+ scope |
| SC-9 | User requests biological interpretation | STOP — premature |
| SC-10 | User requests package not in renv | STOP — requires approval |
| SC-11 | Spatial plot failures | WARNING (not STOP) |
| SC-12 | Missing genes found | WARNING — document, continue with found genes |

---

## 9. Scope Boundary

### Allowed

- Load `seurat_Visium_Hippo_All.rds`
- Read gene list from CSV
- Gene symbol mapping audit
- Region-level expression summaries
- Selected spatial feature plots (≤12 genes)
- Output CSVs + PNGs

### Prohibited

- DESeq2, enrichment, biological interpretation
- WholeBrain/DG/Chromium loading
- Seurat object saving
- PDF generation
- New Rplots.pdf creation
- Full count matrix dense conversion
- Package installation without approval (all needed packages pre-installed in renv)

---

## 10. Packages

All pre-installed in renv (MD5 `139bcfaf7ce24b633c9f034380a15ef5`):

| Package | Version | Use |
|---------|---------|-----|
| Seurat | 5.5.0 | Object loading, expression extraction |
| Matrix | 1.7-3 | Sparse matrix handling |
| ggplot2 | 4.0.3 | All plots |
| pheatmap | 1.0.13 | Heatmap (install if missing, `geom_tile` fallback) |
| digest | 0.6.39 | Rplots.pdf MD5 check |
| data.table | 1.18.4 | fwrite for CSV output |

No package installations expected. If pheatmap missing: install via `renv::install("pheatmap")`, run `renv::snapshot()`, record in provenance.

---

## 11. Rplots.pdf Handling

1. Record pre-existence via `file.exists("Rplots.pdf")` at script start
2. Record MD5 via `digest::digest()` (preferred) or `tools::md5sum()` (fallback)
3. If Rplots.pdf does NOT exist at start AND is created during run: delete it, log deletion
4. If Rplots.pdf EXISTS at start AND is modified during run: log warning, do NOT delete
5. If Rplots.pdf EXISTS at start AND is NOT modified: log no-change
