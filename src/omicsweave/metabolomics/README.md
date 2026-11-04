# 代谢组学分析模块

本模块处理供应商已经完成峰识别和定量的代谢物峰面积表。流程以实验配置驱动，支持 QC 特征门控、样本归一化、探索性分析、经验证的 PLS-DA、差异代谢物、KEGG ORA、GSEA 和自动报告。

## 运行入口

请在实验项目目录中运行（或用 `omicsweave --project-dir <目录>` 指定）：

```powershell
omicsweave metabolomics run --experiment 20250425_cimr_wt_mice_dio
```

也可指定配置或只重跑部分步骤：

```powershell
omicsweave metabolomics run --config experiments/PROJECT/metabolomics_data/config.yaml
omicsweave metabolomics run --experiment PROJECT --steps 1,2,3
python -m omicsweave.metabolomics.metabolomics_pipeline --experiment PROJECT
```

未指定 `--steps` 时，流程根据输入文件、配置参数和各步骤 manifest 增量执行；显式传入 `--steps` 会强制重跑所选步骤，并在必要时自动补齐缺失依赖。

### 按统一项目模型运行

```powershell
omicsweave run --experiment <id> --omics metabolomics            # 步骤 qc、differential、enrichment
```

读取 `experiment.yaml` 中 `omics.metabolomics.input` 声明的峰面积表，分组与比较来自 `subjects.tsv` / `samples.tsv` 和 `design.contrasts`：

```yaml
input:
  level: table
  format: peak_table
  tissue_files: {liver: experiments/<id>/data/liver.xlsx}   # 每个组织一张表；或用 files 列出含多个组织的表
  feature_id_column: Index          # 可选，默认 Index，没有时用名称列
  feature_name_column: Compounds    # 可选，用于 KEGG / HMDB 注释
```

数据表中的样本列名即 `sample_id`，不同时在 samples.tsv 的 `data_column` 列中写出；`sample_type` 为 `qc` / `pool` 的列参与 QC 过滤，`blank` / `standard` 不进入分析。
预处理、PLS-DA、差异检验、KEGG ORA 与代谢物集 GSEA 沿用本模块的函数，同一份数据与下面的旧流程结果逐值一致；图件与解读报告暂时只在旧流程中提供。
结果写入结果包的 `omics/metabolomics/`：`data/<组织>/`（预处理后的数据集）、`qc/<组织>/`、`tables/differential.parquet`、`tables/enrichment.parquet`。

旧版 `--input/--output/--legacy-mode` 已移除。旧实现硬编码小鼠编号与分组，且重复了一套不完整统计逻辑；调用者应迁移到 `config.yaml + groups.yaml`。

## 六个分析步骤

| 步骤 | 功能 | 主要方法 |
|---|---|---|
| 1 | 预处理和 QC | 按组检出率、pooled-QC RSD、逐代谢物最小值插补、PQN、log2 |
| 2 | QC 与多变量分析 | PCA、one-hot PLS-DA、分层交叉验证、标签置换、contrast-specific VIP |
| 3 | 差异分析 | 向量化 Welch t 检验、FDR、多层显著性输出、数据库 ID 注释 |
| 4 | 热图 | 默认按无监督方差选 Top N，再逐代谢物 Z-score 和层次聚类 |
| 5 | KEGG ORA | 统一 FDR 显著集、单侧 Fisher、通路级 FDR、严格注释 universe |
| 6 | GSEA | 默认按 Welch t 统计量排序，KEGG 代谢通路集，1000 次置换 |

CLI 完成后还会生成 Step 7 HTML/Markdown 解读报告；传 `--skip-report` 可跳过。

## 输入约定

每个 `paths.sample_data_<tissue>` 指向一个 `.xlsx` 峰面积表。每行是代谢物，每列是样本或注释。最低需要：

- `Index`：唯一代谢物/峰 ID；
- `Compounds`：代谢物名称；
- 生物样本列：`<mouse_id>_<tissue>`，例如 `A35_sq_adipose`；旧版 `<mouse_id>-<tissue>` 仍可读取；
- pooled QC 列：`QC01`、`QC_02` 等；
- `kegg_map`：逗号分隔的 KEGG pathway ID，用于 ORA/GSEA。

样本 ID 必须存在于 `groups.yaml` 的 `mouse_groups`。未映射样本、零个生物样本或缺失比较组均不会再被静默忽略。

Blank 列可被识别并排除，但当前峰表流程不执行 blank subtraction；如需 blank/sample 比例过滤，应在上游峰处理阶段完成并记录。

## 配置参数

推荐在每个实验的 `metabolomics_data/config.yaml` 中显式保留以下参数：

```yaml
params:
  normalization_method: "pqn"       # pqn / median / tic / none
  detection_rate: 0.5                # 任一组达到阈值即保留
  imputation_fraction: 0.5           # 逐代谢物最小正值 × 0.5
  qc_filter_enabled: true
  qc_rsd_cutoff: 0.30
  qc_min_detection_rate: 0.80
  batch_correction: "none"

  pval_cutoff: 0.05                  # 仅用于名义发现披露
  fdr_cutoff: 0.05                   # 主显著性阈值
  log2fc_cutoff: 0.58
  fdr_method: "bh"
  min_observations_per_group: 2

  vip_cutoff: 1.0
  require_vip_for_significance: false
  plsda_permutation_num: 200
  random_seed: 42

  heatmap_top_n: 50
  heatmap_selection_method: "variance"  # variance / anova

  ora_min_set_size: 3
  ora_max_set_size: 500
  pathway_scope: "metabolism"        # metabolism（不含 ko015xx 耐药图）/ all
  exclude_overview_pathways: true    # 去掉 KEGG Global and overview maps（含 ko012xx）

  gsea_min_set_size: 3
  gsea_max_set_size: 500
  gsea_permutation_num: 1000
  gsea_rank_metric: "welch_t"        # welch_t / signed_logp
```

`batch_correction` 当前只允许 `none`。流程没有足够的样本元数据和设计矩阵来保证 ComBat 等方法不会删除真实处理效应；多批次应先在各批次内完成同口径分析，再做结果级整合。

## 显著性的唯一口径

差异结果同时保留多个证据层，但只有 `significant` / `Significant_<comparison>` 是下游主判定：

- `nominal_significant`：原始 P 与 fold change 通过，仅用于披露；
- `significant_univariate`：FDR 与 fold change 通过；
- `vip_pass`：该 comparison 的 PLS-DA VIP 达到阈值；
- `significant_with_vip`：FDR、fold change、VIP 均通过；
- `significant`：默认等于 `significant_univariate`；仅当 `require_vip_for_significance: true` 时等于 `significant_with_vip`。

火山图着色、DM CSV、ORA 输入和解读报告均消费同一个主判定。PLS-DA 无法交叉验证时，VIP 写为 NaN，避免小样本监督模型被误作证据。

当任一比较组少于 `min_observations_per_group` 时，流程仍输出组均值差和样本量，但 P、FDR 与 VIP 判定为不可推断，不产生“显著”结论。

## 输出结构

```text
metabolomics_results/
  00_dm_results/<tissue>/
    differential_analysis.csv
    dm_<tissue>_<comparison>.csv
    Pathway_<comparison>.csv
  01_data_preprocessing/<tissue>/
    filtered_data.csv
    preprocessing_summary.json
    preprocessing_summary.txt
    normalization_factors.csv
    sample_total_signal.csv
    qc_feature_metrics.csv
  02_qc_visualization/<tissue>/
    PCA.png
    PCA_data.csv
    QC_RSD_distribution.png
    sample_total_signal.png
    PLS-DA.png
    PLS-DA_validation.csv
    PLS-DA_<comparison>_vip.csv
    PLS-DA_<comparison>_validation.csv
    PLSDA_validation_summary.csv
  03_differential_analysis/<tissue>/
    differential_analysis.csv
    differential_design_report.csv
    dm_<tissue>_<comparison>.csv
    Volcano_<comparison>.png
  04_heatmap/<tissue>/
    Heatmap_Top<N>.png
    Heatmap_Top<N>_features.csv
  05_pathway_enrichment/<tissue>/
    Pathway_<comparison>.csv
    Pathway_Bar_<comparison>.png
    Pathway_Dot_<comparison>.png
  06_gsea/<tissue>/
    gsea_<comparison>_results.csv
    gsea_<comparison>.png
  07_interpretation_report/
    metabolomics_report_<experiment>.html
    per_tissue/*.md
```

`00_dm_results` 保留给 multi-omics loader；`03_differential_analysis` 是完整统计输出。两个位置的 per-comparison DM 文件由同一函数写出，不再有两份实现。

## 方法边界

- 输入是相对峰面积，不等同于绝对浓度；
- 代谢物鉴定可信度取决于供应商注释，流程不会把名称映射当作 Level 1 鉴定；
- 默认 PQN 假设大多数代谢物不发生同方向整体变化；若实验确有全局代谢池改变，应结合内标和绝对定量判断；
- PLS-DA 是监督探索工具，必须结合 Q²、交叉验证准确率和置换 P，不可只依据 score plot；
- ORA/GSEA 默认只分析 KEGG metabolism 范围（ko00001–ko01999，不含 ko015xx 耐药图），并排除全部 Global and overview maps（ko01100/ko01110/ko01120 与 ko012xx、ko01310/ko01320）；可通过配置放宽，但疾病/信号图的共享代谢物注释容易造成过度解释；
- n 很小时 FDR 功效有限；原始 P 发现会单独报告，但不能替代多重检验结论。

## 验证

```powershell
python -m pytest tests/metabolomics -v
python -m compileall -q src/omicsweave/metabolomics
```

聚焦测试覆盖：PQN、QC RSD、按组检出过滤、逐代谢物插补、FDR 方法与 NaN、one-hot PLS-DA、VIP、Welch 检验、ORA、GSEA 排名、热图选择、报告口径以及临时 Excel 全链路集成。
