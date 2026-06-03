# ferroDG 复现计划（v3 - 方案B完整版 + renv）

## 1. 项目概述

复现 ferroDG (https://github.com/ShilongZhang116/ferroDG) 的阶段二分析和全部图表。
研究主题：小鼠齿状回(Dentate Gyrus)神经发生过程中的铁死亡(Ferroptosis)分析，使用单细胞RNA测序(GSE233363)，比较 Young vs Old 组。

## 2. 数据来源与 Seurat Object 区分

### 2.0 当前仓库结构说明

- 当前维护的 scRNA-seq 复现脚本位于 `R/scrna/01_load_and_subset.R` 到 `R/scrna/05_correlation_heatmaps.R`。
- `run_pipeline.R` 是 scRNA-seq 分支的入口。
- `analysis/` 目录已退役；不要在其中新增脚本或输出。
- 仓库结构规范见 `docs/repository_structure.md`。

### 2.1 涉及的 Seurat Object

| 文件名 | 来源 | 用途 | 是否需要 |
|--------|------|------|----------|
| `Seurat_combined_with_Celltype_Article.rds` | 用户从原始 FASTQ 重新比对生成 | **阶段二主管道输入** | ✅ 已有 |
| `seurat_Chromium_All.rds` | 原始文章作者提供 (GitHub Article_RDS/) | 阶段一标签迁移参考 | ❌ 不需要 |
| `GSE233363_RAW.tar` | GEO 原始数据 (MTX 格式) | 构建 Seurat object 的原始数据 | ❌ 不需要 |

### 2.2 关键说明

- 用户提供的 `Seurat_combined_with_Celltype_Article.rds` 是**独立从 FASTQ 比对**的结果，与原文章作者提供的数据处理流程可能有差异
- 上游 ferroDG 原始仓库曾使用 `analysis/stage-2/` 路径；本仓库的当前实现已整理为 `R/scrna/01_*` 到 `R/scrna/05_*`，并使用 `Seurat_combined_with_Celltype_Article.rds`
- 本复现**不涉及阶段一**（标签迁移/QC/去双细胞），直接从已标注的 Seurat object 开始

## 3. 系统环境

| 组件 | 版本 | 备注 |
|------|------|------|
| OS | macOS (darwin, aarch64) | Apple Silicon |
| R | 4.5.1 | |
| Seurat | 5.5.0 | **v5！需处理兼容性** |
| dplyr | 1.1.4 | |
| ggplot2 | 4.0.3 | |
| tidyr | 1.3.1 | |
| pheatmap | 1.0.13 | |
| ggrepel | 0.9.6 | |
| renv | 需安装 | 用于环境隔离 |
| 内存 | 16GB | 接受较慢执行 |

### 3.1 Seurat v5 兼容性处理

Seurat v5 的主要变化：
- `FindMarkers()` 返回 `avg_log2FC`（v3/v4 为 `avg_logFC`）
- Layer 结构变化（`counts` → `layers`），但 `@meta.data` 访问方式不变
- `AddModuleScore()` 接口不变
- 代码中需用 `case_when` 兼容两种列名

## 4. 复现范围

复现上游 `analysis/stage-2/Step1-New_Cell_Group.R`（797 行）的阶段二逻辑；当前仓库实现位于 `R/scrna/01_*` 到 `R/scrna/05_*`，生成 Fig 0-7：

| 图表 | 描述 | 关键函数 |
|------|------|----------|
| Fig 0 | 全细胞类型 UMAP | `DimPlot` |
| Fig 1 | 神经发生谱系 UMAP（按年龄分面） | `DimPlot` + `split.by` |
| Fig 2A | Promoter 模块评分小提琴图 | `VlnPlot` |
| Fig 2B | Inhibitor 模块评分小提琴图 | `VlnPlot` |
| Fig 2C | 调控基因小提琴图 | `VlnPlot` |
| Fig 3A | Promoter 评分折线图（Young vs Old） | `ggplot` |
| Fig 3B | Inhibitor 评分折线图 | `ggplot` |
| Fig 3C | 调控基因折线图 | `ggplot` |
| Fig 4A | 神经发生轨迹 - Promoter | `ggplot` |
| Fig 4B | 神经发生轨迹 - Inhibitor | `ggplot` |
| Fig 4C | 神经发生轨迹 - 调控基因 | `ggplot` |
| Fig 5 | 火山图（每个细胞类型 DE） | `FindMarkers` + `ggplot` |
| Fig 6A/B/C | 铁死亡基因 × 神经发生标志物相关性热图 | `cor` + `pheatmap` |
| Fig 7A/B/C | 固定顺序相关性热图（统一颜色标尺 [-0.5, 0.5]） | `pheatmap` |

## 5. 执行步骤

### Step 1: renv 环境初始化

```r
# 安装 renv（如未安装）
install.packages("renv")

# 初始化项目 renv
renv::init()

# 安装依赖包
renv::install(c("Seurat", "dplyr", "tidyr", "ggplot2", "ggrepel", "pheatmap"))

# 创建 renv.lock
renv::snapshot()
```

**输出**: `renv/` 目录、`renv.lock`、`.Rprofile`

### Step 2: 数据验证（运行时检查）

在 R 中执行验证脚本，确认 Seurat object 兼容性：

```r
obj <- readRDS("data/raw/GSE233363_custom/Seurat_combined_with_Celltype_Article.rds")

# 必需的 metadata 列
stopifnot("Celltype_Article" %in% colnames(obj@meta.data))
stopifnot("timepoint" %in% colnames(obj@meta.data))

# UMAP reduction（Fig 0 和 Fig 1 需要）
stopifnot("umap" %in% names(obj@reductions))

# timepoint 值验证
cat("timepoint values:\n")
print(table(obj$timepoint, useNA = "ifany"))

# 细胞类型验证
cat("Cell types:\n")
print(table(obj$Celltype_Article, useNA = "ifany"))

# 基因符号格式检查（应为 Mouse 格式：首字母大写，其余小写）
cat("First 20 gene names:\n")
print(head(rownames(obj), 20))

# 检查关键铁死亡基因是否存在
key_genes <- c("Tfrc", "Gpx4", "Slc7a11", "Acsl4", "Fth1", "Nfe2l2")
cat("Key ferroptosis gene check:\n")
print(key_genes %in% rownames(obj))
```

### Step 3: 运行分析管道（拆分脚本）

按顺序执行拆分后的脚本：

```
R/scrna/01_load_and_subset.R      → data/processed/neuro_processed.rds
R/scrna/02_module_scores.R         → data/processed/neuro_with_ferroptosis_scores.rds
R/scrna/03_figures_0_to_4.R        → figures/stage-2/Fig0-4*.png + .pdf
R/scrna/04_de_and_volcano.R        → figures/stage-2/Fig5*.png + .pdf + CSV
R/scrna/05_correlation_heatmaps.R  → figures/stage-2/Fig6-7*.png + .pdf
```

### Step 4: 输出验证

检查以下目录：
- `data/processed/`: 应有 `neuro_processed.rds`, `neuro_with_ferroptosis_scores.rds`
- `figures/stage-2/`: 应有 Fig0-7 的 PNG 和 PDF 文件
- `gene_lists/`: 应有 promoter/inhibitor/regulator 基因列表

## 6. 脚本详细设计

### 6.1 `R/scrna/01_load_and_subset.R`

**输入**: `data/raw/GSE233363_custom/Seurat_combined_with_Celltype_Article.rds`
**输出**: `data/processed/combined_processed.rds`, `data/processed/neuro_processed.rds`

**逻辑**:
1. 加载 Seurat object，设置 `DefaultAssay(combined) <- "RNA"`
2. 创建 `Age_collapsed`: Young/Young1/Young2 → Young, Old/Old1/Old2 → Old
3. 定义神经发生谱系: qNSC → nIPC → Neuroblast → GC
4. 子集化: 只保留 `Celltype_Article %in% neuro_ct` 且 `Age_collapsed %in% c("Young", "Old")` 的细胞
5. 设置 factor levels
6. 保存中间对象

**运行时检查**:
- 验证 `timepoint` 列存在
- 验证 `Celltype_Article` 列存在
- 打印子集化前后的细胞数
- 打印每个细胞类型 × 年龄组的细胞数

**内存管理**:
- 使用 `gc()` 在子集化后释放内存
- 保存中间对象后释放 `combined`（`rm(combined); gc()`）

### 6.2 `R/scrna/02_module_scores.R`

**输入**: `data/processed/neuro_processed.rds`
**输出**: `data/processed/neuro_with_ferroptosis_scores.rds`, `gene_lists/`

**逻辑**:
1. 加载 `neuro_processed.rds`
2. 定义铁死亡基因集（Human 符号）:
   - Promoter: 25 genes
   - Inhibitor: 40 genes (实际交叉后约 37 个)
   - Regulator: 4 genes
3. Human→Mouse 转换: `tolower()` + 首字母大写
4. 与 `rownames(neuro)` 取交集
5. `AddModuleScore()` 计算 FerroPromoter 和 FerroInhibitor
6. 保存基因列表到 `gene_lists/`
7. 保存 Seurat object

**运行时检查**:
- 验证至少 20 个 Promoter 基因存在
- 验证至少 30 个 Inhibitor 基因存在
- 打印缺失的基因列表
- 验证模块评分列已添加到 metadata

**内存管理**:
- `set.seed(1)` 确保可重复性
- 使用 `gc()` 在计算后释放内存

### 6.3 `R/scrna/03_figures_0_to_4.R`

**输入**: `data/processed/neuro_with_ferroptosis_scores.rds`
**输出**: `figures/stage-2/Fig0*.png/pdf`, `Fig1*.png/pdf`, `Fig2*.png/pdf`, `Fig3*.png/pdf`, `Fig4*.png/pdf`

**逻辑**:
1. 加载带评分的 Seurat object
2. 定义颜色方案:
   - Cell types: qNSC=#FDAE6B, nIPC=#3182BD, Neuroblast=#08519C, GC=#9ECAE1
   - Age: Young=#4DBBD5, Old=#D73027
3. Fig 0: DimPlot 全细胞类型 UMAP
4. Fig 1: DimPlot 神经发生谱系 UMAP，按年龄分面，添加细胞数标注
5. Fig 2A/B: VlnPlot Promoter/Inhibitor 评分
6. Fig 2C: VlnPlot 调控基因（带 guard）
7. Fig 3A/B: 折线图（Young vs Old 均值）
8. Fig 3C: 调控基因折线图（带 guard）
9. Fig 4A/B: 神经发生轨迹（中位数）
10. Fig 4C: 调控基因轨迹（带 guard）

**Seurat v5 兼容性**:
- `FetchData()` 在 v5 中行为不变
- `VlnPlot()` 的 `split.by` 参数兼容
- `DimPlot()` 兼容

**内存管理**:
- 每个图保存后调用 `gc()`
- 使用 `save_plot_local()` helper 统一保存

### 6.4 `R/scrna/04_de_and_volcano.R`

**输入**: `data/processed/neuro_with_ferroptosis_scores.rds`
**输出**: `figures/stage-2/Fig5_*.png/pdf`, `data/processed/DE_results_*.csv`

**逻辑**:
1. 加载 Seurat object
2. 对每个细胞类型（qNSC, nIPC, Neuroblast, GC）:
   - 子集化该细胞类型的细胞
   - 设置 Idents 为 Age_collapsed
   - `FindMarkers(test.use = "wilcox", ident.1 = "Old", ident.2 = "Young")`
   - 处理 Seurat v5 的 `avg_log2FC` 列名兼容
   - 保存 CSV
   - 绘制火山图（标注 top 10 上调/下调基因）

**Seurat v5 兼容性**:
- `FindMarkers()` 返回的列名在 v5 中为 `avg_log2FC`
- 使用 `case_when` 兼容 `avg_logFC` (v3/v4) 和 `avg_log2FC` (v5)

**内存管理**:
- 每个细胞类型子集化后用完即释放
- DE 结果保存为 CSV 后释放内存
- 这是最耗时的步骤，可能需要 5-10 分钟

### 6.5 `R/scrna/05_correlation_heatmaps.R`

**输入**: `data/processed/neuro_with_ferroptosis_scores.rds`
**输出**: `figures/stage-2/Fig6*.png/pdf`, `Fig7*.png/pdf`

**逻辑**:
1. 加载 Seurat object
2. 定义神经发生标志物: Sox2, Hopx, Ascl1, Dcx, Neurod1, Prox1, Calb2, Nestin, Nes
3. 定义铁死亡基因列表（从 gene_lists/ 或硬编码）
4. 获取基因表达矩阵（FetchData）
5. 计算 Spearman 相关系数矩阵
6. Fig 6A: pheatmap - All cells
7. Fig 6B: pheatmap - Young only
8. Fig 6C: pheatmap - Old only
9. Fig 7A: 固定顺序热图 - All（颜色标尺 [-0.5, 0.5]）
10. Fig 7B: 固定顺序热图 - Young
11. Fig 7C: 固定顺序热图 - Old

**内存管理**:
- 相关性矩阵计算可能占用大量内存
- 使用 `gc()` 在每个热图保存后释放

## 7. 内存管理策略

| 阶段 | 预估内存 | 策略 |
|------|----------|------|
| 加载完整 Seurat object | 2-4 GB | 立即子集化，释放完整对象 |
| 神经发生谱系子集 | 500MB-1GB | 保存后释放 |
| 模块评分计算 | 1-2 GB | 计算后 gc() |
| DE 分析（每个细胞类型） | 1-2 GB | 逐个处理，用完释放 |
| 相关性热图 | 1-3 GB | 逐个计算，gc() |
| **峰值内存** | **~4-6 GB** | 16GB 系统足够 |

## 8. 输出文件清单

```
data/processed/
├── combined_processed.rds          # 完整对象（带 Age_collapsed）
├── neuro_processed.rds             # 神经发生谱系子集
├── neuro_with_ferroptosis_scores.rds  # 带铁死亡评分
├── DE_results_qNSC.csv             # DE 结果
├── DE_results_nIPC.csv
├── DE_results_Neuroblast.csv
└── DE_results_GC.csv

figures/stage-2/
├── Fig0_UMAP_AllCelltypes.png/pdf
├── Fig1_UMAP_Neuro_Celltype_SplitByAge.png/pdf
├── Fig2A_Violin_PromoterScore.png/pdf
├── Fig2B_Violin_InhibitorScore.png/pdf
├── Fig2C_Violin_RegulatorGenes.png/pdf
├── Fig3A_Line_PromoterScore.png/pdf
├── Fig3B_Line_InhibitorScore.png/pdf
├── Fig3C_Line_RegulatorGene_*.png/pdf
├── Fig4A_Trajectory_PromoterScore.png/pdf
├── Fig4B_Trajectory_InhibitorScore.png/pdf
├── Fig4C_Trajectory_Regulator_*.png/pdf
├── Fig5_Volcano_qNSC.png/pdf
├── Fig5_Volcano_nIPC.png/pdf
├── Fig5_Volcano_Neuroblast.png/pdf
├── Fig5_Volcano_GC.png/pdf
├── Fig6A_Correlation_Heatmap_All.png/pdf
├── Fig6B_Correlation_Heatmap_Young.png/pdf
├── Fig6C_Correlation_Heatmap_Old.png/pdf
├── Fig7A_FixedOrder_Heatmap_All.png/pdf
├── Fig7B_FixedOrder_Heatmap_Young.png/pdf
└── Fig7C_FixedOrder_Heatmap_Old.png/pdf

gene_lists/
├── ferroptosis_promoter.txt
├── ferroptosis_inhibitor.txt
└── ferroptosis_regulator.txt
```

## 9. 风险与缓解措施

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|----------|
| Seurat v5 列名变化 | 高 | 中 | `case_when` 兼容两种列名 |
| 基因符号不匹配 | 中 | 中 | 运行时验证，打印缺失基因 |
| UMAP reduction 不存在 | 低 | 高 | 运行时 stopifnot 检查 |
| 内存不足 | 低 | 中 | 逐脚本执行，gc() 管理 |
| cairo_pdf 不可用 | 低 | 低 | macOS 默认支持 |
| DE 分析耗时过长 | 中 | 低 | 接受较慢执行，添加进度消息 |
| 调控基因缺失（如 Ras） | 高 | 低 | 带 guard 的条件执行 |

## 10. 与原始代码的差异

| 方面 | 原始代码 | 本复现 | 原因 |
|------|----------|--------|------|
| 脚本结构 | 单文件 797 行 | 5 个拆分脚本 | 便于调试和内存管理 |
| 环境管理 | 无 | renv | 确保可重复性 |
| 基因列表 | 硬编码 Human 符号 | 硬编码 + 运行时验证 | 匹配原始代码 |
| 工作目录 | `file.path(dirname(getwd()))` | 项目根目录相对路径 | 更可靠 |
| Seurat 版本 | v3/v4 | v5 兼容 | 系统为 v5.5.0 |
| 错误处理 | 无 | 运行时检查 + 消息 | 便于调试 |

## 11. 空间转录组扩展 — 脚手架阶段 (Phase Spatial-00)

### 11.1 目的

为 ferroDG 项目的空间转录组扩展搭建目录结构、文档和配置框架。本阶段**不包含任何分析实现**，仅建立脚手架。

### 11.2 核心共识（已确认）

| 共识 | 说明 |
|------|------|
| scRNA-seq 细胞 ≠ Visium spot | Visium spot 可能捕获多个细胞的混合信号，不能视为单细胞分辨率 |
| Spot 不是独立生物学重复 | 相邻 spot 共享局部微环境，空间自相关性违反独立性假设 |
| 元数据必须从实际对象中检查 | 不能假设 assay 名称、spatial key、图像插槽或解剖标签 |
| 引用路由 | 参考 `.opencode/skills/ref-bio/reference-pack/references.link-only.yaml` |
| 基因符号约定 | 小鼠格式（首字母大写，其余小写） |

### 11.3 预期但未验证的数据来源

- **来源**: JessbergerLab / Wu et al. / GSE233363
- **当前输入**: `data/raw/GSE233363_official/` 中的作者提供空间 RDS 对象（seurat_Visium_DG_All.rds, seurat_Visium_Hippo_All.rds, seurat_Visium_WholeBrain.rds）
- **状态**: DG 和 Hippo 已完成 Stage Spatial-01/02 检查；WholeBrain 仅清点，不在当前阶段加载
- **raw Visium 状态**: 作者 Script 8 的 STutility 段需要 `filtered_feature_bc_matrix.h5`、`tissue_positions_list.csv`、`tissue_hires_image.png`、`scalefactors_json.json` 等外部 raw Visium/Space Ranger 文件；这些文件当前本地不可用
- **当前路线决策**: Phase Spatial-05 已基于作者 RDS 中已检查到的 image slots 和坐标进行 RDS-based reproduction/approximation。该路线必须明确标注为近似复现，不能报告为严格作者代码复现。
- **当前科学路线**: 先复现并验证论文中的 hippocampus 相关空间流程，确认区域框架和作者结果的对应关系；随后再把同一套区域框架迁移到 CA1/CA3/DG 线粒体相关基因表达问题。
- **下一优先级**: Phase Spatial-06 应聚焦 hippocampal regional atlas reproduction（CA1、CA2、CA3、ML、GCL、Hilus 和派生 DG），而不是直接进入线粒体分析。

### 11.4 脚手架步骤

| 步骤 | 产出 | 状态 |
|------|------|------|
| 1. 更新 AGENTS.md | 添加空间转录组扩展章节 | ✅ |
| 2. 更新 reproduction_plan.md | 添加本节（第 11 节） | ✅ |
| 3. 创建 docs/reference_sources.md | 引用路由文档 | ✅ |
| 4. 创建 docs/spatial_transcriptomics_overview.md | 脚手架概览 | ✅ |
| 5. 创建 docs/data_manifest_spatial.md | 未验证数据清单 | ✅ |
| 6. 创建 docs/environment_renv_notes.md | 环境和 renv 状态 | ✅ |
| 7. 创建 docs/memory_management.md | 内存管理约束 | ✅ |
| 8. 创建 R/spatial/README.md | 脚手架说明 | ✅ |
| 9. 定义目录约定 | `R/scrna/`、`R/spatial/`、可选 `python/spatial/` 为源码目录；`data/`、`figures/` 等为 gitignored 本地产物目录；不再保留空占位目录 | ✅ |
| 10. 更新 .gitignore | 最小化非冗余添加 | ✅ |
| 11. 最终报告 | 列出所有变更 | ✅ |

### 11.5 约束

- **不安装任何 R 包**
- **不修改 renv.lock**
- **不下载数据**
- **不在复现验证完成前实现 CA1/CA3/DG 线粒体分析**
- **不假设元数据（assay、spatial key、图像、解剖标签）**
- **不在 gitignore 目录中创建 .gitkeep 文件**
- 工具/方法选择须先经过 ref-bio 参考路由和实际对象结构检查，不在本脚手架阶段预设
- **Visium-only**（无 MERFISH、Slide-seq、HD 占位符）
- 如果使用 RDS 替代 raw Visium 文件推进后续空间邻域分析，必须记录与作者 STutility 路线的偏离点，并保留将来找到 raw 文件后重开严格复现分支的可能性
