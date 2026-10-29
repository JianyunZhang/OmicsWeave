# 蛋白质组学分析流程（v6.0）

本模块面向 DIA-NN/Spectronaut 蛋白定量报告，提供配置驱动的质量控制、差异丰度分析、可视化、ssGSEA、pre-ranked GSEA 和证据感知报告。v6.0 修复了旧版 FDR 校正静默失效、n=1 仍产出无效 p 值、火山图用原始 p 值判定显著性、报告/ssGSEA 文件名不一致等问题，并加入协变量设计与完整 QC。

## 统计契约

流程对“用于统计推断的数据”和“用于展示的数据”作严格区分：

1. DIA-NN 报告按 `PG.Q.Value < 0.01`、`Global.PG.Q.Value < 0.01`（存在时）过滤。
2. 同一基因/Feature 下多个 protein group 默认在 log2 空间求均值，即在线性空间取几何均值。
3. 线性强度做 log2 转换及逐样本中位数中心化。
4. **差异分析使用未插补矩阵**；每个蛋白按实际观测计算样本量和残差自由度。
5. 下移正态分布插补有固定随机种子，只用于 PCA、热图和样本级 ssGSEA。
6. 每组少于 2 个有效观测、设计矩阵不满秩或残差自由度不足时，只报告效应量，`evidence_level=insufficient`，p 值与 padj 留空，`significant=False`。
7. FDR 方法必须是 statsmodels 支持的明确名称（默认 `fdr_bh`）；非法值直接报错，不允许退回原始 p 值。
8. 火山图、差异计数和报告统一采用 padj 与 |Log2FC| 双阈值。
9. **火山图纵轴画原始 p，配色仍然只由 padj 决定**（`classify_significance()` 不变）。
   BH 是阶梯函数：单调性修正会把成片蛋白绑在同一个校正值上，本项目实测单个
   panel 里最多有 8474 个蛋白（占全部可检验蛋白的 100%）共享同一个 padj，直接
   把 padj 画在纵轴上会把点云压成几条横线，那是校正算法的形状，不是数据的形状。
   BH 对原始 p 单调，所以 FDR 阈值线画在“padj 仍然通过的最大原始 p”处
   （`fdr_cut_on_raw_p()`），纵轴与配色的判定边界完全一致。没有任何蛋白通过时
   不画线，改在副标题写明阈值不可达——空图必须能读出“达不到”而不是“漏画了”。
10. **不可检验的蛋白按数量披露，不静默丢弃。** 每组观测数不足的蛋白没有 p 值，
    落在图上就是 NaN 而整点消失；本项目每个 panel 有 500–1600 个这样的蛋白，
    且它们携带最大的倍数变化（最低到 log2FC = −9.8）。数量写进副标题。

## 模块结构

| 文件 | 功能 |
|------|------|
| `run_proteomics_pipeline.py` | 主编排：可选 DIA-NN 搜库、导入、预处理、QC、差异分析、作图、富集 |
| `diann_search.py` | Thermo `.raw` → DIA-NN 1.9 library-free；样本名预检与 provenance 跳过 |
| `differential_analysis.py` | 设计矩阵、缺失值感知 OLS、limma 方差压缩、FDR、证据等级 |
| `proteomics_qc.py` | 鉴定数、缺失率、CV、重复相关性、强度分布、PCA-协变量关联 |
| `experiment_config.py` | YAML 配置、路径、分组、动物 ID、样本元数据与参数校验 |
| `enrichment/gsea.py` | 样本级 ssGSEA 与比较级 pre-ranked GSEA |
| `plotting/` | PCA、热图、火山图及独立重绘 |
| `generate_interpretation_report.py` | FDR/证据等级感知的 Markdown 报告 |
| `test_proteomics_pipeline.py` | 聚焦回归测试 |

`run_proteomics_analysis.py` 仅是兼容入口，旧版统计实现已经删除；新分析应直接使用 `run_proteomics_pipeline.py`。

## 运行方式

仓库统一使用 Windows conda 环境：

```powershell
omicsweave proteomics run `
  --experiment 20250425_cimr_wt_mice_dio

# 只验证一个组织，不运行在线基因集分析
omicsweave proteomics run `
  --experiment 20250425_cimr_wt_mice_dio --tissues liver --skip-gsea

# 统计实现或配置改变后强制重算。--force 只重算统计，不会重跑 DIA-NN
omicsweave proteomics run `
  --experiment 20250425_cimr_wt_mice_dio --force

# 从 Thermo .raw 搜库后再做差异分析。重跑搜库必须同时带 --search 和 --force
omicsweave proteomics run `
  --experiment <id> --search

# 只搜库，不进入 QC / 差异分析
omicsweave proteomics run `
  --experiment <id> --search-only

# 从 v6 结果独立补图/补富集
omicsweave proteomics plot `
  --experiment 20250425_cimr_wt_mice_dio

# 生成报告
omicsweave proteomics report `
  --experiment 20250425_cimr_wt_mice_dio

# 回归测试
python -m pytest tests/proteomics
```

旧自动化仍可使用 `--experiment_dir <proteomics_data目录> --output_dir <目录>`；它们现在是兼容别名，内部仍读取同一份 `config.yaml`。新脚本应使用 `--experiment/--outdir`。

旧结果没有 `analysis_metadata.json`，可能包含错误的 Padj。主流程、独立绘图和报告模块会拒绝复用它们，并要求使用 `--force` 重跑。

多组学 `ProteomicsLoader` 同样默认拒绝旧结果，并保留 `evidence_*` / `inference_status_*`。蛋白组评分优先使用 padj；只有经过人工审计、显式设置 `integration.inputs.proteomics.allow_legacy_results: true` 时才允许加载没有 v6 元数据的历史表。

## 从 Thermo .raw 搜库

未加 `--search` 时，流程仍然只读 `paths.data_dirs` 里已有的 DIA-NN `report_*.tsv`。两个现有实验不用改配置。

`--search` 在统计步骤之前调用 Windows 上的 DIA-NN 1.9，直接读 Thermo `.raw`（依赖 Thermo MSFileReader），写出 `report_<tissue>.tsv`，再把这份报告交给原来的 `load_and_clean_data()`。`--search-only` 停在搜库。`--force` 单独使用只重算 QC 和差异分析；要重跑 DIA-NN，必须同时给出 `--search`（或 `--search-only`）和 `--force`。

搜库默认关闭。打开时在实验 `config.yaml` 中写：

```yaml
search:
  enabled: true
  engine: diann
  expected_version: "1.9"
  diann_exe: "tools/diann/1.9/diann.exe"          # 相对于实验项目目录
  fasta: "reference/UP000000589_10090_SwissProt_202406.fasta"
  threads: 8
  qvalue: 0.01
  reanalyse: true
  smart_profiling: true
paths:
  raw_dirs:
    brain: "\\\\server\\...\\MSF250153\\raw"
```

`diann_exe` 和 `fasta` 只在 `enabled: true` 时必须存在。FASTA 不自动下载。要和 2024 年 6 月那批 library-free 结果可比，用同一份小鼠 Swiss-Prot；换库会改变蛋白集合和定量值。

每个组织一个目录，只读取该目录下一层的 `*.raw`，不进入子目录。文件名末尾必须是 `_<部位>_<动物ID>.raw`，例如 `MSF250153_zzk_Astral_23min_DIA_100ng_Brain_52.raw`。DIA-NN 写进 `Run` 列的是去掉扩展名的文件名，流程把它解析成 `Brain52`，再按最长后缀对上 `groups.yaml` 里的 `52`。对不上的文件会在启动 DIA-NN 之前失败。

命令对应方法学里的 library-free 默认设置：酶切 `K*,R*`、最多 1 个漏切、肽段长度 7–30、前体电荷 1–4、Cys 固定脲甲基化（`--unimod4`）、Met 氧化与 N 端乙酰化且每个肽最多 1 个可变修饰、precursor q-value 0.01。不加 `--relaxed-prot-inf`，蛋白组 FDR 保持 1%。日志里的版本必须以 `1.9` 开头（接受 `1.9.x`）；对不上时不把报告交给差异分析。

指纹记录 raw 的路径、大小和修改时间，以及 FASTA 的 SHA-256 和上述搜库参数。指纹一致且报告已在时跳过 DIA-NN。线程数不进指纹。更换 DIA-NN 程序本身不会自动被发现，换版本后用 `--search --force` 重跑。

运行前需要：DIA-NN 1.9 的 `diann.exe`（不装进 conda）、Thermo MSFileReader，以及 FASTA。raw 放在 UNC 或本地绝对路径，不复制进仓库。

## 配置

`experiments/<id>/proteomics_data/config.yaml` 的核心参数如下：

```yaml
params:
  pval_cutoff: 0.05
  log2fc_cutoff: 0.58
  heatmap_top_n: 50
  fdr_method: "fdr_bh"
  dp_method: "limma"

  preprocessing:
    min_total_observations: 2
    protein_aggregation: "log2_mean"   # 推荐；linear_mean 仅用于旧结果回溯
    imputation_seed: 42
    imputation_shift: 1.8
    imputation_width: 0.3

  differential_design:
    min_observations_per_group: 2
    covariates: []                      # 示例：[batch, sex]

  enrichment:
    run_ssgsea: true
    run_prerank_gsea: true
    min_set_size: 10
    max_set_size: 500
    permutation_num: 1000
    ssgsea_permutation_num: 0
    seed: 42
    threads: 1
    gene_sets: [kegg, go]
```

### 批次与协变量

协变量必须显式提供，流程不会从文件夹名猜测批次。可在 `groups.yaml` 中写通用元数据：

```yaml
mouse_groups:
  "52": hfd_14w
  "63": lfd_14w

sample_metadata:
  "52":
    batch: batch_1
    sex: male
  "63":
    batch: batch_1
    sex: male
```

键既可以是动物 ID，也可以是解析后的完整样本名；完整样本名优先，可用于同一动物的技术重复处于不同采集批次的情况。

仅配置批次时也可使用简写：

```yaml
sample_batches:
  "52": batch_1
  "63": batch_1
```

随后在 `params.differential_design.covariates` 中列出 `batch`。分类变量会采用参考水平编码，连续数值变量会标准化。若 batch 与 condition 完全混杂，流程会在 `differential_design_report.csv` 中标记 `design_not_identifiable` 并拒绝 p 值推断；不会自动删除 batch 项。

动物 ID 支持纯数字和字母数字（如 `52`、`A35`）。解析时按已配置 ID 的最长后缀匹配，未映射样本会使流程直接失败。

## 差异分析实现

每个比较建立 `~0 + condition + covariates` 设计矩阵，并按蛋白的实际非缺失样本拟合 OLS。残差方差通过 limma 的 scaled-F 经验贝叶斯框架压缩：

- 残差自由度相同时，采用 `fitFDist` 的 log-F 矩估计（digamma/trigamma 修正）。
- 缺失值导致自由度不同时，采用新版 limma `fitFDistUnequalDF1` 的 profile maximum-likelihood 核心。
- 后验方差为 `(df × s² + d0 × s0²) / (df + d0)`，moderated t 使用 `df + d0`。

这是一套纯 Python 实现，输出中保留 `standard_error`、`tstat`、实际组内观测数、残差自由度、`evidence_level` 和 `inference_status`，便于审计。对于需要随机效应、肽段层级模型或复杂重复测量的设计，建议改用 R `MSstats`、`limma`/`dream` 或 `proDA`，而不是强行套用本模块。

## 富集分析

- `ssgsea_<library>/`：每个样本的通路分数，用于描述和聚类，不等同于组间显著性检验。
- `gsea_prerank/<comparison>/<library>/`：使用全部有限 moderated t 统计量的比较级 GSEA，输出 NES、名义 p 值、FDR 和 leading-edge 蛋白。
- 每次运行写入 `ssgsea_status.csv` 与 `gsea_prerank/prerank_status.csv`；某个在线库失败时会明确记录错误，不再静默吞掉。

## 输出结构

```text
proteinomics_results/                    # 沿用项目已有目录名
├── 00_diann_search/<tissue>/            # 仅 --search / --search-only 写出
│   ├── report_<tissue>.tsv
│   ├── report_<tissue>.pg_matrix.tsv
│   ├── report_<tissue>.pr_matrix.tsv
│   ├── diann.log
│   ├── diann_cmd.txt
│   └── search_provenance.json
├── 00_qc/<tissue>/
│   ├── sample_qc_metrics.csv
│   ├── protein_missingness.csv
│   ├── replicate_correlations.csv
│   ├── group_qc_metrics.csv
│   ├── pc_covariate_association.csv
│   ├── intensity_distributions.png
│   ├── sample_detection_missingness.png
│   ├── correlation_heatmap.png
│   ├── missingness_heatmap.png
│   └── qc_summary.json
├── 01_dp_results/<tissue>/
│   ├── normalized_matrix.csv             # 未插补；差异分析输入
│   ├── imputed_matrix.csv                # 仅展示/ssGSEA
│   ├── differential_proteins.csv
│   ├── differential_design_report.csv
│   ├── dp_<tissue>_<comparison>.csv
│   └── analysis_metadata.json
├── 02_result_plot/<tissue>/
│   ├── pca_plot.png
│   ├── heatmap_top50.png
│   ├── volcano_<comparison>.png
│   ├── ssgsea_<library>/
│   └── gsea_prerank/<comparison>/<library>/
└── 03_interpretation_report/
    ├── index.md
    └── <tissue>/<tissue>_report.md
```

## 当前限制

- 当前输入是蛋白/蛋白组层 MaxLFQ；尚未实现肽段层级混合模型。直接读 Thermo `.raw` 只通过外部 DIA-NN 1.9，不在 Python 里做谱图鉴定。
- MNAR 缺失不会被当作真实低丰度观测参与检验；因此“只在一组出现”的蛋白通常会被标为证据不足。此策略偏保守，但避免插补值制造显著性。
- 单批次内固定效应协变量已支持；跨独立采集批次的效应合并、批次内 DE + meta-analysis 尚未在本模块实现。
- Enrichr 在线基因集名称会随服务更新；状态表用于记录不可用库和网络错误。

## 版本迁移

v5 及更早结果必须使用 `--force` 重跑。尤其是旧配置中的 `fdr_method: "bh"` 曾触发 statsmodels 异常并被宽泛 `except` 吞掉，使 Padj 等于原始 p 值；v6 已通过配置白名单、无静默回退和分析元数据三层防护消除此问题。
