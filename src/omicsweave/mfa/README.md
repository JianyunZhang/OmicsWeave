# 代谢流量分析 (MFA) 分析流程说明书

> **位置**：`omicsweave/mfa/`
> **适用平台**：Windows 10/11 原生（无需 WSL/Linux 服务器）
> **运行环境**：Python ≥ 3.10，安装方式见仓库根目录 README
> **配置驱动**：按 experiment 复用，只改 YAML 配置文件，不改代码
---

## 概述

本模块提供了一套**可复用的 13C 代谢通量分析 (Metabolic Flux Analysis, MFA) 下游数据处理框架**，用于对稳定同位素（如 13C-葡萄糖）示踪实验产生的 MFA 拟合结果进行系统性统计分析与可视化。

模块采用**计算与可视化分离**的三层架构：
- **配置层** (`experiments/<id>/mfa_data/`) — 支持 experiment/protocol/batch 分层配置（`experiment.yaml` + `protocol.yaml` + `<batch>/config.yaml` + `<batch>/groups.yaml`）
- **计算层** (`mfa_analysis.py` + `df_analysis.py`) — 项目无关的核心分析引擎，仅负责数值计算与表格输出
- **可视化层** (`plotting/`) — 独立的绘图子包，从已保存的表格重新生成所有图表

新实验只需在 `mfa_data/` 下维护 YAML 配置与批次数据目录，通过 `--experiment <id>` 即可运行完整分析，无需编写任何 Python 代码。图表可独立重新生成，无需重跑分析。

### 适用的上游拟合工具

本模块处理的输入数据为 MFA 拟合软件的输出结果，兼容以下主流工具产生的数据格式：

| 上游拟合工具                        | 典型输出                                  | 兼容性     |
|-------------------------------|---------------------------------------|---------|
| **自研 MFA 软件** (Shiyu Liu Lab) | flux_raw_data.xlsx, mid_raw_data.xlsx | 完全兼容    |
| **INCA** (Vanderbilt)         | 类似 Excel 格式                           | 格式适配后兼容 |
| **Metran** (Northwestern)     | 类似 Excel 格式                           | 格式适配后兼容 |
| **13CFLUX2**                  | 文本/CSV 格式                             | 需预处理转换  |
| **OpenFLUX2**                 | CSV 格式                                | 需预处理转换  |

### 支持的数据类型

| 数据类型 | 输入文件 | 说明 |
|---|---|---|
| **代谢通量 (Flux)** | `flux_raw_data.xlsx` | 多组织、多条件的代谢反应通量拟合结果 |
| **质量同位素分布 (MID)** | `mid_raw_data.xlsx` | 各代谢物的质量同位素分布（M+0, M+1, M+2, ...）拟合结果 |

---

## 改进措施汇总

本节记录流程当前的关键设计决定与它们解决的问题。**这些都是不变式**，
修改相关代码前请先读这一节。

### 1. 推断单位：生物学重复，不是拟合迭代

每只动物先被汇总为**一个**净通量估计（其 IQR 过滤后拟合点的中位数），
统计检验只在动物之间进行。模型拟合的 ~150 次迭代描述的是同一只动物上
参数的不确定度，不是独立观测；当作重复会把样本量虚增两个数量级。

`differential_flux()` 与 `flux_heatmap()` 都必须经由
`_compute_biological_net_flux()` 取值。

### 2. 秩检验的样本量下界与自动回退

Mann-Whitney 的最小双侧 p 值恒为 `2/C(n1+n2, n1)`，**与数据无关**：
n=3 vs 3 时为 0.1，n=3 vs 4 时为 0.057——`p<0.05` 在数学上不可达。
BH 校正后还有第二道下界。`statistic: auto` 在两条路径都不可达时自动回退到
`effect_size_only`（|Hedges' g| + |log2FC|），并打印 `[POWER]` / `[AUTO]`。

否则「0 个显著通量」这一方法学上限会被误读成生物学结论。

### 3. 三层拟合质量门控（结论之前必须先看）

| 层 | 回答的问题 | 模块 |
|---|---|---|
| 灌注稳态 QC | 采样时系统到底稳没稳？ | `infusion_qc.py` |
| MID 拟合残差 | 模型能不能复现**实测**标记模式？ | `mid_fit_residuals()` |
| 通量可辨识性 | 这条通量到底被数据约束住了吗？ | `sample_flux_matrix()` |

三者都独立于下游统计——统计做得再精细，也无法补救未达稳态、失拟或
未被约束的通量。

### 4. 多批次整合：先判路线，再合并

批次不自动等于重复。`build_integration_design_report()` 在合并前把每个
(comparison × tissue) 分类为 `random_effects_meta` / `union_single_batch` /
`union_disjoint_fluxes` / `excluded_legacy_run`，**FDR 按路线分池**，
并输出 `is_cross_batch_supported`（最强档）与 `cross_batch_direction_agree`
（弱档）两级证据。

### 5. 技术噪声地板

同一只动物在两个批次里各拟合一次，其差值不含生物学信号，构成跨批次比较的
噪声地板。`bridging_qc.py` 自动识别桥接样本并估计该地板，`is_cross_batch_supported`
已把它纳入判定。前提是 `sample_id` 能跨批次唯一标识同一只动物——
做不到时结果会被标注为上界而非地板。

### 6. 配置真正生效

`analysis_defaults` 的三级覆盖链（experiment → protocol → batch）、
`infusion_qc_defaults`、`batch_metadata`、`allow_cross_protocol`、
`tissue_aliases` 此前均已在 YAML 中声明但**从未被代码读取**，现已全部接通。
仍未实现的字段集中列在文末，不要假设它们在生效。

### 7. 图注与文档必须与代码一致

火山图 / 热图的统计说明由 `primary_statistic` 驱动动态生成，不再写死
「150 次迭代」「Cohen's d」等与实际实现不符的表述。

---

## 项目结构

```
omicsweave/
  common/                          # 跨组学共享模块
    reaction_pathway_mapping.py      MFA 反应 → KEGG/Reactome/EC 数据库 ID
    metabolite_id_mapping.py         代谢物 → HMDB/KEGG/ChEBI
    protein_id_mapping.py            蛋白质 → UniProt/EC
    mouse_gene_id_mapping.py         基因 → NCBI/Ensembl
    nonparametric.py                 秩检验的样本量下界（p / BH-FDR）
    meta_analysis.py                 DerSimonian-Laird 随机效应 Meta
    stability.py                     组内 bootstrap + 标签置换校准

  mfa/
    mfa_analysis.py                # 核心引擎（MFAAnalyser）— 数据加载 + 流程编排
    df_analysis.py                 # 分析函数（计算 + CSV → 01_df_results/）
    mfa_experiment_config.py       # 配置加载器（MFAConfig + load_mfa_config）
    infusion_qc.py                 # 灌注期血糖稳态 QC
    provenance.py                  # 运行溯源指纹 + 跨批次可合并性断言
    bridging_qc.py                 # 桥接样本一致性 → 技术噪声地板
    mfa_integration.py             # 多批次 Meta + 跨 protocol 一致性 + 整合报告
    mfa_integration_plotting.py    # 整合结果绘图
    run_mfa_analysis.py            # CLI 入口（单批次 / --all-batches / --integrate）
    plot_mfa_results.py            # 独立绘图入口
    generate_mfa_interpretation_report.py   # 按组织解读报告（Step 9）
    plotting/                      # 绘图子包（与计算解耦）
      __init__.py                    公共 API: generate_all_plots()
      qc.py                          拟合 loss / 血糖稳态 / MID 残差 / 可辨识性
      volcano.py                     火山图 + 组织差异汇总
      heatmap.py                     差异通量聚类热图
      pca.py                         PCA 散点 + Scree
      pathway.py                     通路箱线图
      overview.py                    组织通量概览
      mid.py                         MID 标记分析图
    README.md                      # 本说明文档

experiments/{experiment_id}/mfa_data/
  experiment.yaml                  # 实验级：tracer_protocols 登记、整合策略、analysis_defaults
  {protocol_id}/
    protocol.yaml                  # protocol 级：示踪剂、灌注参数、稳态 QC 阈值
    {batch_id}/
      config.yaml                  # 批次运行入口：batch/groups 引用、comparisons、analysis 覆盖
      batch.yaml                   # 批次级：data_files、检测元数据、batch_qc
      groups.yaml                  # 条件定义 + 样本→条件映射
      flux_raw_data.xlsx           # 通量拟合结果
      mid_raw_data.xlsx            # MID 拟合结果（含实测 experiments 列）
      *血糖记录*.xlsx               # 灌注期血糖，供稳态 QC

analysis_results/{experiment_id}/mfa_results/
  runs/{protocol_id}/{batch_id}/
    run_provenance.json            # 技术链指纹
    01_df_results/
      all_tissues/                 # 跨组织：qc/ flux_matrix/ pca/ overview/ mid/ pathway_flux/
      {tissue}/                    # 组织特异：differential_flux/ pathway_flux/ heatmap/
    02_result_plot/
      all_tissues/                 # 0_qc/ 3_pca_clustering/ 6_tissue_overview/ 7_mid_analysis/
      {tissue}/                    # 1_volcano/ 2_heatmap/ 4_pathway_flux/
    03_interpretation_report/      # index.md + 每组织 .md
  integration/
    tracer_harmonized/{protocol}/  # 路线判定表 + 合并结果 + 汇总
    bridging_qc/{protocol}/        # 桥接对 + 技术噪声地板
    cross_tracer_scoring/          # 跨 protocol 方向印证
    02_result_plot/                # 森林图 / 异质性 / 一致性（可选）
    integration_report.md
```

### 架构设计

MFA 模块采用**配置层 + 计算层 + 可视化层**三层架构，与实验数据完全解耦：

```
┌─────────────────────────────────┐
│     配置层（实验特定，随数据存放）   │
│  experiments/<id>/mfa_data/     │
│                                 │
│  experiment.yaml                │
│  <protocol>/protocol.yaml       │
│  <protocol>/<batch>/config.yaml │
│  <protocol>/<batch>/groups.yaml │
│  flux_raw_data.xlsx / mid_raw_data.xlsx │
└──────────────┬──────────────────┘
               │ load_mfa_config()
               ▼
┌─────────────────────────────────┐     ┌─────────────────────────────────┐
│     计算层（可复用，项目无关）       │     │     可视化层（独立，可单独运行）    │
│  mfa_analysis.py                │     │  plotting/                      │
│                                 │     │                                 │
│  MFAAnalyser                    │     │  volcano.py  heatmap.py         │
│    load → compute → save CSV    │────►│  pca.py      pathway.py         │
│                                 │     │  overview.py  qc.py             │
│  run_mfa_analysis.py            │     │  mid.py                         │
│    CLI 入口（--experiment）      │     │                                 │
│                                 │     │  plot_mfa_results.py            │
└─────────────────────────────────┘     │    独立绘图入口（--experiment）  │
                                        └─────────────────────────────────┘
```

**关键设计原则**：
- 数据加载（`MFAAnalyser`）与分析计算（`df_analysis` 模块）解耦
- `df_analysis` 中的独立函数仅输出 CSV 到 `01_df_results/` 对应子目录，不直接调用 matplotlib/seaborn
- 绘图模块从 `01_df_results/` 读取 CSV 生成图表到 `02_result_plot/`
- 组织特异性分析（差异通量、通路比较、热图）按 `<tissue>/<category>/` 分层存储，避免跨组织文件混淆
- 跨组织分析（PCA、QC、概览、MID）统一存储在 `01_df_results/all_tissues/`
- TCA 相关信息统一由 `pathway_flux_reaction_stats.csv` 中 `pathway == "TCA_Cycle"` 派生，不再单独输出 `tca_cycle/` 目录
- MFA 输出目录按 `mfa_results/runs/<protocol>/<batch>/` 分层，避免不同 tracer/batch 结果互相覆盖
- `mfa_results/integration/` 作为跨 tracer 汇总与联合评分固定出口
- `all_tissues/` 目录按设计仅包含 `0_qc/3_pca_clustering/6_tissue_overview/7_mid_analysis`，不包含 `1_volcano/2_heatmap/4_pathway_flux`
- 多组学集成 loader 通过递归 glob (`**/`) 自动发现所有组织子目录的文件并合并
- 所有差异通量输出自动附带公共数据库通路注释（ec_number, kegg_pathway, reactome_id），由 `omicsweave/common/reaction_pathway_mapping.py` 提供统一映射
- 绘图层兼容新旧 `01_df_results` 布局（`all_tissues/` 新路径 / 历史根层级旧路径），重绘时会统一输出到标准 `02_result_plot/<tissue>/...` 路径
- 两个流程可独立运行：`run_mfa_analysis.py`（支持多运行模式），`plot_mfa_results.py` 仅重新生成图表

---

## 安装依赖

```bash
pip install pandas numpy scipy scikit-learn matplotlib seaborn openpyxl statsmodels
```

### Python 核心软件包列表

| 包名 | 版本要求 | 用途 | 在 MFA 分析中的应用 |
|---|---|---|---|
| **pandas** | ≥1.3.0 | 数据处理 | Excel 多 sheet 读取、通量矩阵构建、分组统计 |
| **numpy** | ≥1.20.0 | 数值计算 | 矩阵运算、统计量计算 |
| **scipy** | ≥1.7.0 | 统计检验 | Welch's t-test、Mann-Whitney U 检验 |
| **scikit-learn** | ≥1.0.0 | 降维分析 | PCA 主成分分析、StandardScaler 数据标准化 |
| **matplotlib** | ≥3.5.0 | 基础绘图 | 火山图、柱状图、散点图 |
| **seaborn** | ≥0.11.0 | 高级可视化 | 聚类热图 (Clustermap) |
| **openpyxl** | ≥3.0.0 | Excel 读写 | 读取 .xlsx 格式拟合结果 |
| **statsmodels** | ≥0.13.0 | 多重检验校正 | Benjamini-Hochberg FDR 校正 |

---

## 输入数据格式规范

### 1. flux_raw_data.xlsx

Excel 工作簿，包含以下结构：

**`sheet_information` 工作表**（必需）：

| 列名 | 类型 | 说明 |
|---|---|---|
| `complete_result_name` | str | 样本唯一标识符，如 `glucose__male__hfd__1` |
| `sheet_name` | str | 对应数据工作表名称，如 `glucose_male_hfd_1` |
| `labeling` | str | 示踪底物类型，如 `glucose` |
| `strain` | str | 品系/性别等分类信息 |
| `condition` | str | 实验条件分组，如 `hfd`, `lfd`, `hsd` |
| `index` | int | 组内样本编号 |

**数据工作表**（每个样本一个 sheet）：

| 列名 | 说明 |
|---|---|
| `Index` | 拟合迭代序号 |
| `loss` | 残差平方和 (SSR)，用于选取最优拟合 |
| `{tissue}--{reaction}` | 通量值，格式为 `组织名--反应名`，如 `liver--CS_m` |
| `{tissue}--{reaction}__R` | 逆反应通量（可逆反应），后缀 `__R` 表示反向 |
| `MIX:{tissue}--{metabolite}_{idx}` | 代谢物混合比例参数 |

**通量列命名规则**：
- `{tissue}--{reaction_name}` — 正向反应通量
- `{tissue}--{reaction_name}__R` — 逆向反应通量（可逆反应）
- `MIX:{tissue}--{metabolite}_{index}` — 代谢物池混合比例

### 2. mid_raw_data.xlsx

Excel 工作簿，结构与 flux 文件类似：

**`sheet_information` 工作表**：与 flux 文件格式相同。

**数据工作表**：

| 列名 | 说明 |
|---|---|
| `metabolite_name` | 代谢物标识符，如 `liver--PYR_c`、`serum--GLC_e` |
| `mid_index` | 质量同位素编号（0 = M+0, 1 = M+1, ...） |
| `experiments` | 实验测量值 |
| `1`, `2`, ... `N` | 各拟合迭代的模拟 MID 值 |

第一行 `metabolite_name = "loss"` 存储各迭代的拟合 loss 值，用于选取最优迭代。

---

## 分析流程详解

### 流程总览

```
输入（experiments/<id>/mfa_data/<protocol>/<batch>/）
  flux_raw_data.xlsx   mid_raw_data.xlsx   血糖记录.xlsx
  experiment.yaml + protocol.yaml + batch.yaml + config.yaml + groups.yaml
                              │
                              ▼
┌─ Step 1  数据加载 ──────────────────────────────────────────────────┐
│   load_flux_data()   每样本取 loss 最小的拟合行 + 保留全部迭代        │
│   load_mid_data()    保留模拟 MID 与**实测 MID**（experiments 列）    │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ Step 2  前提检验与质量控制（结论之前必须先看）───────────────────────┐
│   2.1 infusion_steady_state_qc()  灌注期血糖是否达稳态                │
│   2.2 quality_control()           拟合 loss 分布                      │
│   2.3 sample_flux_matrix()        样本级净通量 + 通量可辨识性          │
│   2.4 mid_fit_residuals()         实测 vs 模拟 MID 残差               │
│   2.5 write_provenance()          技术链指纹 run_provenance.json      │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ Step 3-7  分析计算（CSV → 01_df_results/）──────────────────────────┐
│   3. flux_pca()                 通量 PCA（按组织）                    │
│      flux_heatmap()             差异通量热图数据（生物学重复 + 下界）  │
│   4. differential_flux()        差异通量                              │
│        ├─ 生物学重复为推断单位                                        │
│        ├─ 样本量下界判定 → auto 回退 effect_size_only                 │
│        ├─ Hedges' g / Cohen's d / log2FC(+SE)                        │
│        ├─ BH-FDR 校正                                                 │
│        ├─ 可辨识性标注                                                │
│        └─ 通路注释（EC / KEGG / Reactome）                            │
│   5. pathway_flux_comparison()  10 个通路的反应级比较                 │
│   6. tissue_flux_overview()     组织通量概览                          │
│   7. mid_analysis()             MID 标记度                            │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ Step 8  可视化（读 01_df_results/ → 写 02_result_plot/）────────────┐
│   qc.py       → all_tissues/0_qc/   （含 3 张新增质控图）             │
│   volcano.py  → <tissue>/1_volcano/                                  │
│   heatmap.py  → <tissue>/2_heatmap/                                  │
│   pca.py      → all_tissues/3_pca_clustering/                        │
│   pathway.py  → <tissue>/4_pathway_flux/                             │
│   overview.py → all_tissues/6_tissue_overview/                       │
│   mid.py      → all_tissues/7_mid_analysis/                          │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ Step 9  解读报告 → 03_interpretation_report/（按组织 .md + index）──┐
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ Step 10  多批次 / 多协议整合（--integrate，输出到 integration/）────┐
│   assert_batch_compatibility()      技术链指纹比对                    │
│   bridging_sample_qc()              技术噪声地板                      │
│   build_integration_design_report() 路线判定（先读这张表）            │
│   within_protocol_meta_analysis()   批次间随机效应 Meta，FDR 按路线分池│
│   cross_protocol_consensus()        跨示踪剂方向印证（默认关闭）      │
│   generate_integration_report()     integration_report.md            │
└─────────────────────────────────────────────────────────────────────┘
```

**分层存储约定**：组织特异性分析按 `<tissue>/<category>/` 存放；跨组织分析
（QC / PCA / overview / MID / flux_matrix）统一在 `all_tissues/`。

### Step 1: 数据加载

#### 1.1 通量数据加载 (`load_flux_data`)

**目的**：从 `flux_raw_data.xlsx` 中提取每个样本的最优拟合通量向量。

**处理流程**：
1. 读取 `sheet_information` 工作表，获取样本元数据（条件分组、样本ID）
2. 遍历每个样本对应的数据工作表
3. 根据 `loss` 列选取残差最小（最优拟合）的行作为该样本的代表通量
4. 排除非通量列（`Index`, `loss`），构建 **样本 × 通量** 矩阵
5. 解析通量列名，自动识别组织列表和条件分组

**关键选择**：
- **Best-fit 选取策略**：取 `loss` 最小值对应行。若无 loss 列则取所有迭代均值
- **通量列解析**：`tissue--reaction` 格式自动分离为组织和反应名

#### 1.2 MID 数据加载 (`load_mid_data`)

**目的**：从 `mid_raw_data.xlsx` 中提取每个样本的最优拟合 MID 向量。

**处理流程**：
1. 读取 `sheet_information` 获取样本映射
2. 对每个样本工作表：提取 `loss` 行以确定最优迭代列
3. 选取最优迭代列对应的 MID 模拟值
4. 构建 **代谢物(metabolite, mid_index) × 样本** 矩阵
5. 同时保留 `experiments` 列（**实测** MID），供 Step 2.4 残差 QC 使用

### Step 2: 前提检验与质量控制

这一组检查回答的是下游统计无法回答的问题：系统稳没稳、模型拟没拟上、
通量约束住没有。**任何结论产生之前都应先看这四张表。**

#### 2.1 灌注稳态 QC (`infusion_steady_state_qc`)

**目的**：检验 MFA 最根本的前提——采样时代谢与同位素是否已达稳态。

13C-MFA 的全部推断都建立在这个假设上：只有稳态成立，「一组固定的通量」
才是对该时刻代谢状态的有效描述。若灌注期内血糖持续上升或剧烈波动，
拟合出的通量是一段漂移过程的平均，既不对应任何真实生理状态，也不应
拿来做组间比较。

**方法**：读取批次的血糖记录表（`batch.yaml: data_files.glycemia_record_file`），
在 `protocol.yaml: infusion_qc_defaults.steady_state_window_min` 窗口内逐动物
计算变异系数、极差与线性漂移速率，按 `glucose_cv_max` / `glucose_range_max_mmol_l`
判定。这些阈值早已在配置中声明，但此前从未被任何代码读取。

**输出**：
- `01_df_results/all_tissues/qc/infusion_steady_state_qc.csv` — 每动物一行
- `01_df_results/all_tissues/qc/infusion_glycemia_timecourse.csv` — 长表时间曲线
- `01_df_results/all_tissues/qc/infusion_steady_state_by_condition.csv`
- `01_df_results/all_tissues/qc/infusion_steady_state_unmapped.csv` — 在血糖表中但不在动物分组里的动物
- `02_result_plot/all_tissues/0_qc/qc_infusion_glycemia.png`

**动物分组**：血糖表按动物编号记录，条件取自 `groups.yaml: mouse_groups`
（动物编号 → 条件），缺失时退回 `sample_groups`。后者在部分实验里以 MFA
样本名为键（`20250613_cimr_wt_mice_dio_glp`），不能当动物分组用。不在映射中
的动物（如 GLP 葡萄糖批次里改输 U-13C 谷氨酰胺的 85 号）照常计算，但不计入
稳态判定与按条件汇总，单独写入 `infusion_steady_state_unmapped.csv`；若映射与
血糖表完全没有交集，则视为配置错误，所有动物留在主表并打印 `[WARN]`。

**样本对应的局限**：血糖表用动物编号（`A12`），MFA 结果用「组名+序号」
（`glucose__male__hfd_sham_8w__1`），两者无法逐一对应，因此 QC 结论只按
动物和条件两级报告，不下推到单个 MFA 样本（`on_fail: exclude` 因此降级为
`flag`）。要做样本级门控，需上游在 `complete_result_name` 写入动物编号。

#### 2.2 拟合质量评估 (`quality_control`)

**目的**：评估 MFA 模型拟合质量，检测潜在的异常样本。

**方法**：
- 统计每个样本的所有拟合迭代 loss 值（min, max, mean, n_fits）
- 按实验条件分组，评估组内 loss 一致性
- 异常 loss 值可能提示：模型拓扑与实际代谢网络不匹配、实验数据质量问题、或代谢状态剧烈波动

**输出**：
- `01_df_results/all_tissues/qc/qc_fitting_loss.csv` — 拟合质量汇总表
- `02_result_plot/all_tissues/0_qc/qc_fitting_loss.png` — 柱状图 + 直方图

**文献依据**：
- Antoniewicz MR. (2015). "Methods and advances in metabolic flux analysis: a mini-review". *Journal of Industrial Microbiology & Biotechnology*, 42(3), 317-325.

#### 2.3 样本级通量矩阵与可辨识性 (`sample_flux_matrix`)

**目的**：输出整个下游唯一保留**样本粒度**的表，供三类用途：

1. 跨批次桥接样本一致性检查（技术噪声地板估计，见「多批次整合」一节）；
2. 推断单位审计——确认差异分析与热图确实建立在生物学重复上；
3. 后续通量置信区间 / 可辨识性工作的输入。

**方法**：每只动物 × 每条反应输出一个净通量代表值（IQR 过滤后拟合点的
中位数），同时记录该动物在该反应上的拟合点离散度（`fit_p25` / `fit_p75` /
`fit_std` / `n_fit_points`）作为约束程度的粗略代理。

同时输出**每通量可辨识性**：`fit_rel_iqr = (p75 − p25) / |net_flux|`。
拟合迭代不是生物学重复，但它们确实刻画了「在这只动物的数据下，这条通量
被约束得有多紧」。相对 IQR 接近或超过 1，意味着该通量几乎不被数据决定，
其组间差异不携带信息。

| median_rel_iqr | identifiability |
|---|---|
| < 0.2 | well_constrained |
| < 1.0 | moderate |
| ≥ 1.0 | poorly_constrained |

差异通量结果会自动带上 `median_rel_iqr` 与 `identifiability` 两列；
`analysis.differential_flux.require_identifiable: true`（默认 `false`）可进一步
把 `poorly_constrained` 的通量强制置为不显著。

**输出**：
- `01_df_results/all_tissues/flux_matrix/sample_net_flux.csv`
- `01_df_results/all_tissues/flux_matrix/flux_identifiability.csv`
- `02_result_plot/all_tissues/0_qc/qc_flux_identifiability.png`

#### 2.4 MID 拟合残差 QC (`mid_fit_residuals`)

**目的**：回答 13C-MFA 最基本的问题——**模型到底拟合得好不好**。

`mid_raw_data.xlsx` 的 `experiments` 列是**实测** MID，此前被读入后直接丢弃。
没有它，只剩 `loss` 这一个汇总标量，看不出模型在哪个代谢物、哪个同位素体上
失拟。

**方法**：逐 (样本 × 代谢物) 比较实测与最优拟合 MID，计算 SSR / RMSE /
最大绝对残差，再按代谢物汇总。MID 是合计约为 1 的分数，因此阈值直接用
绝对分数表示（典型质谱 MID 测量误差 0.003–0.02）：

| RMSE | 判定 |
|---|---|
| < 0.02 | good |
| < 0.05 | acceptable |
| ≥ 0.05 | poor（模型无法复现实测标记模式） |

同时检查实测 MID 是否归一到 1（`measured_sum_off_by`），偏离过大提示上游
归一或自然丰度校正有问题。

**输出**：
- `01_df_results/all_tissues/mid/mid_fit_residuals.csv` — 每 (样本 × 代谢物)
- `01_df_results/all_tissues/mid/mid_fit_quality_by_metabolite.csv` — 按代谢物汇总
- `01_df_results/all_tissues/mid/mid_residuals_long.csv` — 逐同位素体明细
- `02_result_plot/all_tissues/0_qc/qc_mid_fit_residuals.png`

**为什么必须先看这个**：模型复现不了实测 MID 时，由它导出的通量不具解释力——
无论其统计量多漂亮。应优先检查网络拓扑、稳态假设与自然丰度校正。

#### 2.5 运行溯源 (`write_provenance`)

**目的**：记录本次运行的技术链指纹，供整合阶段判断两个批次是否真的可合并。

**记录内容**：示踪剂种类与同位素、**反应集合哈希**（代谢网络指纹）、
提取方法、检测仪器与检测单位、拟合迭代次数、组织与条件列表、样本数、
本次生效的统计阈值。

**输出**：`runs/<protocol>/<batch>/run_provenance.json`

整合阶段由 `provenance.assert_batch_compatibility()` 比对：示踪剂 / 网络哈希 /
提取方法不一致属**关键**不一致（合并前需人工确认），仪器与迭代次数差异仅提示。
行为由 `integration.provenance_mismatch_policy`（`warn` 默认 | `fail` | `ignore`）控制。

### Step 3: 全局模式分析

#### 3.1 通量 PCA (`flux_pca`)

**目的**：通过无监督降维评估样本间的全局代谢通量模式差异。

**方法**：
1. 对通量矩阵进行 Z-score 标准化（`StandardScaler`）
2. 执行 PCA（默认提取前 5 个主成分）
3. 以 PC1 vs PC2 散点图展示条件间分离度
4. Scree plot 展示各主成分方差贡献率

**输出**：
- `01_df_results/all_tissues/pca/flux_pca_scores.csv` — PCA 得分矩阵（含条件标签、组织）
- `01_df_results/all_tissues/pca/flux_pca_loadings_{tissue}.csv` — 各组织 PCA 载荷矩阵
- `02_result_plot/all_tissues/3_pca_clustering/flux_pca_{tissue}.png` — 各组织 PCA 散点图 + Scree plot

**解读指南**：
- 条件组在 PCA 空间中的聚集/分散程度反映代谢通量模式的一致性/异质性
- 载荷矩阵可用于识别驱动条件间分离的关键通量

#### 3.2 聚类热图 (`flux_heatmap`)

**目的**：展示条件间显著差异的通量的样本级别聚类模式。

**方法**（与 `differential_flux` 使用同一套推断单位与样本量下界判定）：
1. 每只动物先汇总为一个净通量估计（forward - reverse，取拟合点中位数）
2. p 值路径可达时按 Mann-Whitney p < `pvalue_threshold` 筛选；
   不可达时（如 n=3 vs 3）保留全部反应并按 |Hedges' g| 取 top-N
3. Z-score 标准化后进行行聚类（样本保持按条件排列）
4. 以条件颜色标注样本列

`heatmap_meta_{comp}.csv` 记录 `effect_size_mode` / `rank_test_underpowered` /
`min_achievable_pvalue`，绘图层据此选择坐标轴与图注文案。

**输出**：
- `01_df_results/heatmap/flux_heatmap_raw_{tissue}_{comp}.csv` — 热图原始数据矩阵
- `01_df_results/heatmap/flux_heatmap_stats_{tissue}_{comp}.csv` — 显著通量统计摘要
- `02_result_plot/<tissue>/2_heatmap/flux_heatmap_{tissue}_{comp}.png` — 聚类热图

### Step 4: 差异通量分析 (`differential_flux`)

**目的**：鉴定条件间统计学显著差异的代谢通量。

**推断单位：生物学重复（不是拟合迭代点）**

每只动物先被汇总为**一个**净通量估计（其 IQR 过滤后拟合点的中位数），
统计检验只在这些动物级估计之间进行。模型拟合的 150 次迭代描述的是
**同一只动物上参数的不确定度**，不是独立的生物学观测；把它们当作重复
会把样本量虚增两个数量级，产生大量假显著。`_compute_biological_net_flux()`
负责这一汇总，`differential_flux` 与 `flux_heatmap` 都必须经由它取值。

**方法**：
- 对每一对比较条件 (A vs B)、每个通量执行：
  1. 计算净通量（forward - reverse for 可逆反应），每只动物取中位数
  2. 计算均值差 (diff)、log2FC 及其 delta-method 标准误
  3. 计算 Cohen's d 与小样本校正的 **Hedges' g**（效应量主指标）
  4. **Mann-Whitney U 检验** — 非参数检验
  5. **Welch's t-test** — 作为参考补充
  6. **Benjamini-Hochberg FDR 校正** — 控制假发现率
- 排除 `MIX:` 前缀列（代谢物混合比例参数，非代谢通量）

**样本量下界与自动回退（重要）**

秩检验的 p 值是离散的：`n1` vs `n2` 时共有 `C(n1+n2, n1)` 种秩排列，
最小双侧 p 值恒为 `2/C`，**与数据无关**。MFA 常见样本量下：

| n1 vs n2 | 最小可达双侧 p | p<0.05 是否可能 |
|---|---|---|
| 1 vs 1 | 1.0 | 否 |
| 2 vs 2 | 0.333 | 否 |
| 3 vs 3 | **0.1** | **否** |
| 3 vs 4 | **0.057** | **否** |
| 4 vs 4 | 0.029 | 是 |
| 5 vs 5 | 0.0079 | 是 |

此外 BH 校正后的最小可达 FDR 还要再乘以 `n_tests / n_at_floor`——
在 474 条通量上做校正时，即便 n=4 vs 4 的未校正 p 可达 0.029，
最小可达 FDR 仍可能远高于 0.05。

`statistic: auto`（默认）会在**两条路径都不可达**时自动回退到
`effect_size_only`，用基于生物学重复的 |Hedges' g| + |log2FC| 选择，
并在日志中打印 `[POWER]` / `[AUTO]` 提示。这样避免把
「0 个显著通量」这一组合数学上限误读成生物学结论。
实现见 `omicsweave/common/nonparametric.py`。

**显著性阈值**（由 `analysis.differential_flux` 配置，见下文覆盖链）：
- `pvalue_threshold` / `fdr_threshold` — p 值路径
- `cohen_d_threshold`（`hedges_g_threshold` 未配置时沿用该值）+ `log2fc_threshold` — 效应量路径
- 任一组生物学重复 < 2 时 `evidence_level = insufficient`，`is_significant` 强制为 False

**输出**：
- `01_df_results/<tissue>/differential_flux/differential_flux_{A}_vs_{B}.csv` — 按组织存储的差异分析结果表
- `02_result_plot/<tissue>/1_volcano/volcano_{comp}.png` — 火山图
- `02_result_plot/<tissue>/1_volcano/tissue_diff_summary_{comp}.png` — 单组织差异通量上调/下调汇总

**差异分析结果表列说明**：

| 列名 | 说明 |
|---|---|
| `flux_id` | 通量完整标识符 (`tissue--reaction`) |
| `tissue` | 组织名 |
| `reaction` | 反应名 |
| `n_points_{cond}` | 该条件 IQR 过滤后的数据点数 |
| `mean_{cond}` / `std_{cond}` | 各条件的净通量均值/标准差 |
| `diff` | 均值差 (mean_A - mean_B) |
| `log2FC` | log2(|mean_A| / |mean_B|) |
| `pvalue_mannwhitney` | Mann-Whitney U p 值（主要指标） |
| `pvalue_ttest` | Welch's t-test p 值（参考） |
| `fdr_mannwhitney` / `fdr_ttest` | BH-FDR 校正值 |
| `inference_unit` | 恒为 `biological_replicate` |
| `cohen_d` / `hedges_g` | 效应量；`hedges_g` 为小样本校正后的主指标 |
| `log2fc_se` / `diff_se` | 效应量标准误（delta method / 组间合并） |
| `min_achievable_pvalue` | 该样本量下秩检验的最小可达双侧 p（与数据无关） |
| `min_achievable_fdr` | BH 校正后的最小可达 FDR |
| `rank_test_underpowered` | p 值路径（含 FDR）在该样本量下是否根本不可达 |
| `fdr_unreachable` | 未校正 p 可达但 BH-FDR 不可达 |
| `evidence_level` | `statistical_test` / `effect_size_only` / `insufficient` |
| `effect_size_forced_by_power_floor` | 是否因样本量下界被自动切到效应量路径 |
| `primary_statistic` | 本次配置实际采用的主统计方法 |
| `is_significant` | 按配置阈值综合判断的显著性布尔值 |
| `significance_rule` | 生成 `is_significant` 所用规则 |
| `stat_note` | biological_replicates / fitting_variance_only / low_n_configured_statistical_test |
| `ec_number` | 酶编号 (EC Number)，如 `2.7.1.1` |
| `kegg_reaction` | KEGG 反应 ID，如 `R00771` |
| `kegg_pathway` | KEGG 通路 map ID，如 `map00010` |
| `reactome_id` | Reactome 反应 ID，如 `R-MMU-70493` |
| `description` | 反应英文描述 |

**通路注释说明**：以上 5 个注释列由 `omicsweave/common/reaction_pathway_mapping.py` 中的 `annotate_flux_dataframe()` 函数自动添加，无需手动干预。映射表包含 137 个 MFA 反应 ID，覆盖糖酵解、TCA 循环、磷酸戊糖途径、脂肪酸代谢、氨基酸代谢、酮体代谢等主要代谢通路。注释列极大地方便了多组学集成分析时按 KEGG/Reactome 通路进行跨组学 join 操作。

**文献依据**：
- Crown SB, Antoniewicz MR. (2013). "Publishing 13C metabolic flux analysis studies: a review and future perspectives". *Metabolic Engineering*, 20, 42-48.
- Benjamini Y, Hochberg Y. (1995). "Controlling the false discovery rate". *Journal of the Royal Statistical Society*, 57(1), 289-300.

### Step 5: 通路级通量比较

#### 5.1 代谢通路通量比较 (`pathway_flux_comparison`)

**目的**：按代谢通路分组，系统比较各组织中核心代谢反应的通量水平。

**预定义通路分组** (10 个)：

| 通路 | 包含反应 | 生物学功能 |
|---|---|---|
| Glycolysis (糖酵解) | HEX, PGI, PFK, FBA, TPI, GAPD, PGK, PGM, ENO, PYK | 葡萄糖氧化分解 |
| TCA Cycle (TCA 循环) | CS, ACONT, ICDH, AKGD, SUCOAS, SUCD, FUMH, MDH | 乙酰辅酶A氧化 |
| Pentose Phosphate (磷酸戊糖途径) | G6PDH2R_PGL_GND, RPI, RPE, TKT1, TKT2, TALA | NADPH 生产 + 核苷酸前体 |
| Pyruvate Metabolism (丙酮酸代谢) | PDH, LDH, PC, PEPCK, ME2 | 丙酮酸分流节点 |
| Amino Acid (氨基酸代谢) | GLUD, GLND, GLNS, ASPTA, GPT | 氨基酸转氨与脱氨 |
| Lipid Synthesis (脂质合成) | ACITL, LIPID, ACSS2, ACOT, FA_usage | 从头脂肪合成 |
| Gluconeogenesis (糖异生) | G6Pase, FBPase, PEPCK | 肝脏/肾脏糖输出 |
| Ketogenesis (酮体生成) | ACAC_production, BHBDH, BHB_trans, BHB_output | 脂肪酸衍生酮体 |
| Glycerol Metabolism (甘油代谢) | GK, G3PP, GPS | 甘油-3-磷酸通路 |
| Transport (转运反应) | PYR_trans, ASPGLU, AKGMAL, CIT_trans, GLN_trans, GLU_trans | 线粒体-胞浆穿梭 |

**注意**：通路分组可通过 `pathway_groups` 配置参数自定义覆盖。分析时自动排除逆反应 (`__R` 后缀) 以避免重复计数。

**输出**：
- `01_df_results/<tissue>/pathway_flux/pathway_flux_samples.csv` — 按组织存储的通量数据
- `01_df_results/<tissue>/pathway_flux/pathway_flux_reaction_stats.csv` — 按组织 × 通路 × 反应 × 条件统计（`mean/std/n/n_points`，TCA 风格）
- `01_df_results/<tissue>/pathway_flux/pathway_flux_summary_stats.csv` — 按组织 × 通路 × 条件汇总统计（`mean/std/n/n_points/n_reactions`）
- `01_df_results/all_tissues/pathway_flux/pathway_flux_reaction_stats.csv` — 跨组织反应统计总表
- `01_df_results/all_tissues/pathway_flux/pathway_flux_summary_stats.csv` — 跨组织通路汇总总表
- `02_result_plot/<tissue>/4_pathway_flux/pathway_flux_{tissue}.png` — 每个组织的通路箱线图
- `02_result_plot/<tissue>/4_pathway_flux/pathway_flux_{tissue}_individual.png` — 组织内全部反应独立分布图（总览）
- `02_result_plot/<tissue>/4_pathway_flux/by_pathway/<pathway_name>/pathway_flux_{tissue}_{pathway_name}.png` — 单一通路总反应箱线图
- `02_result_plot/<tissue>/4_pathway_flux/by_pathway/<pathway_name>/pathway_flux_{tissue}_{pathway_name}_individual.png` — 单一通路反应拆分箱线图

### Step 6: 组织通量概览 (`tissue_flux_overview`)

**目的**：比较各组织在不同条件下的总体代谢活跃程度。

**方法**：
1. 对每个组织，计算所有正向反应通量的绝对值之和（总通量）和均值（平均通量/反应）
2. 按条件分组比较

**输出**：
- `01_df_results/all_tissues/overview/tissue_flux_overview.csv` — 组织通量概览表
- `02_result_plot/all_tissues/6_tissue_overview/tissue_flux_overview.png` — 总通量 + 平均通量柱状图

### Step 7: MID 分析 (`mid_analysis`)

**目的**：分析同位素标记在各代谢物中的分布模式，评估底物碳流入各代谢途径的程度。

**方法**：
1. 按条件分组计算每个代谢物各 MID（M+0, M+1, ...）的均值和标准差
2. 计算**总标记度** (Total Labeling)：sum(M+1, M+2, ..., M+n) = 1 - M+0
3. 识别条件间总标记度差异最大的代谢物
4. 对高标记代谢物展示完整 MID 分布

**输出**：
- `01_df_results/all_tissues/mid/mid_summary.csv` — 完整 MID 汇总表
- `01_df_results/all_tissues/mid/mid_total_labeling.csv` — 各代谢物总标记度
- `02_result_plot/all_tissues/7_mid_analysis/mid_total_labeling.png` — Top 差异代谢物标记度柱状图
- `02_result_plot/all_tissues/7_mid_analysis/mid_distributions.png` — 高标记代谢物的 MID 分布图

**生物学意义**：
- 总标记度越高，表示示踪底物（13C-葡萄糖）的碳原子流入该代谢物越多
- 条件间总标记度的差异反映代谢路径使用的改变
- M+0 的增加意味着来自非标记碳源（如脂肪酸、氨基酸）的贡献增大

### Step 8: 报告生成 (`generate_report`)

自动生成 Markdown 格式的分析摘要报告，列出数据概况和所有产出的表格、图表文件。

### Step 9: 解读报告生成 (`generate_mfa_interpretation_report`)

**目的**：按组织（tissue）生成独立的 Markdown 解读报告，包含统计摘要、关键发现和解读模板。

**处理流程**：
1. 从 `01_df_results/` 和 `02_result_plot/` 读取已有分析结果
2. 自动发现组织列表（从 `tissue_flux_overview.csv` 或 PCA scores）
3. 自动发现比较对（从 `differential_flux/` 目录的文件名）
4. 对每个组织生成一个独立的 `.md` 文件
5. 生成 `index.md` 目录页（含实验概况、文件清单、各组织报告链接）
6. 若存在旧的 `mfa_report.md`，自动删除（内容已合并到 `index.md`）

**每个组织报告页面包含**：

| 板块 | 数据来源 | 内容 |
|---|---|---|---|
| 组织通量概览 | `all_tissues/overview/` | 总通量、平均通量 × 条件对比表 |
| PCA 降维概况 | `all_tissues/pca/` | PC1 均值描述 + 散点图链接 |
| 差异通量分析 | `<tissue>/differential_flux/` | 每对比较的显著通量数、Top 上调/下调通量表 |
| TCA 循环通量 | `<tissue>/pathway_flux/pathway_flux_reaction_stats.csv` | 由 `pathway == TCA_Cycle` 派生 |
| 代谢通路对比 | `<tissue>/pathway_flux/` + `all_tissues/pathway_flux/` | 10 个通路的样本级分布 + 反应统计 + 通路汇总 |
| MID 同位素标记 | `all_tissues/mid/` | 各代谢物总标记度 × 条件对比 |
| 解读模板 | — | 填空式结论模板，便于写论文/汇报 |
| 下一步建议 | — | 验证实验建议 |

**输出**：
```
03_interpretation_report/
  index.md          # 目录页
  liver.md          # 每组织独立解读报告
  heart.md
  kidney.md
  ...
```

**单独运行**：
```bash
# 仅重新生成解读报告（不重新计算）
omicsweave mfa report --experiment <id>

# 或通过 run_mfa_analysis.py 的 --report-only 参数
omicsweave mfa run --experiment <id> --report-only
```

### Step 10: 多批次 / 多协议整合 (`--integrate`)

单批次分析完成后的独立步骤，输出到 `mfa_results/integration/`。执行顺序有
依赖关系，不可调换：

| 顺序 | 操作 | 为什么在这个位置 |
|---|---|---|
| 1 | `assert_batch_compatibility()` | 技术链不兼容时，后面所有合并都无意义 |
| 2 | `bridging_sample_qc()` | 噪声地板要先算出来，Meta 才能用它判断效应是否超出技术重现性 |
| 3 | `build_integration_design_report()` | 合并之前先判定每个 (comparison × tissue) 该走哪条路线 |
| 4 | `within_protocol_meta_analysis()` | 批次间随机效应 Meta，FDR 按路线分池 |
| 5 | `cross_protocol_consensus()` | 跨示踪剂方向印证（需 `allow_cross_protocol: true`） |
| 6 | `generate_integration_report()` | 汇总为 `integration_report.md` |
| 7 | `generate_integration_plots()` | 需 `integration.plotting.enabled: true` |

各步的判定规则、输出列含义与方法学局限见下文「方法 6: 多批次整合」。

---

## 使用方法

### 标准执行流程（速查）

```bash
cd <repo_root>
```

**单批次实验**（一个 protocol、一个 batch）：

```bash
omicsweave mfa run --experiment <id>
```

**多批次实验**（推荐顺序，两步）：

```bash
omicsweave mfa run --experiment <id> --all-batches --with-plots
```

```bash
omicsweave mfa run --experiment <id> --integrate
```

> 整合必须在各批次单批分析**完成之后**运行——它读取的是各批次已写出的
> `differential_flux_*.csv` 与 `sample_net_flux.csv`。

**结果阅读顺序**（这个顺序不是可选的）：

| 顺序 | 看什么 | 位置 |
|---|---|---|
| 1 | 灌注是否达稳态 | `01_df_results/all_tissues/qc/infusion_steady_state_qc.csv` |
| 2 | 模型是否拟合上 | `01_df_results/all_tissues/mid/mid_fit_quality_by_metabolite.csv` |
| 3 | 通量是否被约束 | `01_df_results/all_tissues/flux_matrix/flux_identifiability.csv` |
| 4 | 本次用了哪条统计路径 | 差异结果的 `evidence_level` / `primary_statistic` 列 |
| 5 | 差异通量 | `01_df_results/<tissue>/differential_flux/` |
| 6 | 多批次：路线判定与噪声地板 | `integration/tracer_harmonized/<p>/integration_design_report_<p>.csv`、`integration/bridging_qc/<p>/bridge_noise_floor.csv` |
| 7 | 整合结论 | `integration/integration_report.md` |

前三项任一不通过时，第 5–7 项的数值不具解释力——应先回到实验或拟合环节。

### 方法 1: 完整分析（计算 + 可视化 + 解读报告）

```bash
# 激活 conda 环境

# 使用实验 ID 运行完整分析流程（计算 → 存表 → 绘图 → 解读报告）
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio
```

这是标准入口。模块根据实验 ID 自动定位配置文件 `experiments/<experiment_id>/mfa_data/config.yaml`，读取数据路径、组别映射、颜色等参数。若 `config.yaml` 中未指定 `comparisons`，则自动从 `groups.yaml` 的 `conditions` 生成所有两两比较。

### 方法 2: 仅重新生成解读报告（不重新计算）

```bash
# 从已有分析结果重新生成按组织拆分的解读报告
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio --report-only
```

或直接调用报告生成器：
```bash
omicsweave mfa report --experiment 20250425_cimr_wt_mice_dio
```

### 方法 3: 仅输出 differential_flux（快速迭代）

```bash
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio --differential-only
```

适合参数调试和快速比较，仅输出 `01_df_results/<tissue>/differential_flux/*.csv`，跳过其他分析步骤、绘图与解读报告。

### 方法 4: 全量计算但不绘图（节省绘图耗时）

```bash
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio --no-plots
```

会执行数据分析并输出各类 CSV（含 QC/PCA/heatmap/pathway/overview/MID），但跳过 `02_result_plot` 绘图和 `03_interpretation_report` 生成。

### 方法 5: 仅重新生成图表（不重新计算）

```bash
# 从已保存的 01_df_results/ CSV 重新生成所有可视化图表
omicsweave mfa plot --experiment 20250425_cimr_wt_mice_dio
```

适用于以下场景：
- 修改了绘图参数（颜色、阈值等）后重新出图
- 添加了新的绘图类型后补生成
- 原始数据未变，仅需更新图表样式

### 方法 6: 多批次整合

```bash
# 仅执行整合（假设各 batch 的单批结果已完成）
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio --integrate

# 先批量跑各 batch（仅算表），再执行整合
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio --all-batches --integrate

# 每个批次都出图 + 解读报告，再执行整合（耗时较长）
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio --all-batches --with-plots --integrate
```

> 当 `experiment.yaml` 中 `integration.plotting.enabled: true` 时，
> 整合流程会在生成 `integration_report.md` 后自动读取整合 CSV 并出图到
> `mfa_results/integration/02_result_plot/`。

#### 6.0 两个整合轴：批次 vs 示踪剂

同一个 experiment 下可以同时存在多个示踪剂 protocol、每个 protocol 下多个
batch。流程对两条轴都支持，但它们**不是同一种操作**：

| 轴 | 前提假设 | 做法 | 入口 |
|---|---|---|---|
| batch × batch（同 protocol） | 各批次是同一测量的重复 | 随机效应 Meta，按 (comparison × tissue) 分路线 | `within_protocol_meta_analysis()` |
| protocol × protocol（不同示踪剂） | 不同示踪剂探测网络的**不同部分**，互补而非重复 | 方向一致性 + 证据门槛 | `cross_protocol_consensus()` |

跨 protocol 分析默认**关闭**，由 `integration.allow_cross_protocol` 控制。

**为什么方向一致不等于证据**：符号只有两种取值，纯噪声也有约一半概率一致。
因此 `consensus_pass` 还要求：

| 条件 | 配置项 | 默认 |
|---|---|---|
| 方向完全一致 | `cross_protocol.require_same_direction` | `true` |
| 至少 N 个 protocol 判定显著 | `cross_protocol.min_protocols_significant` | 同 `min_protocols` |
| 该通量在各 protocol 中均非 `poorly_constrained` | `cross_protocol.exclude_poorly_constrained` | `true` |

输出列 `n_protocols_significant` / `n_protocols_well_constrained` /
`excluded_poorly_constrained` 记录每条通量被哪一条挡下。

> **方法学局限**：并行标记实验的正确处理方式是**在拟合阶段联合多示踪剂数据**
> （Crown & Antoniewicz 2013），而非对各自独立拟合出的通量图事后投票。
> 本步骤的结果应作为「相互印证」，不能当作「相互验证」。

#### 6.1 整合路线判定（合并之前必须先读）

在做任何合并之前，`build_integration_design_report()` 会先回答一个问题：
**这些批次到底是彼此的重复，还是互补的组织板块？** 判定结果按
(comparison × tissue) 写入：

`integration/tracer_harmonized/<protocol>/integration_design_report_<protocol>.csv`

| route | 触发条件 | 含义 |
|---|---|---|
| `random_effects_meta` | ≥`min_batches` 个批次覆盖该组合且反应名可对齐 | 随机效应合并有意义，构成跨批次验证 |
| `union_single_batch` | 只有一个批次覆盖 | **只能并集拼接**，不是 pooled effect，不构成跨批次证据 |
| `union_disjoint_fluxes` | 多批次覆盖同一组织但反应集不相交 | 通常是模型网络/命名不一致，需人工核对 |
| `excluded_legacy_run` | 组织目录在磁盘上，但本次运行的 `sample_net_flux.csv` 里没有 | 更早一版输入数据 + 更早一版代码的产物。**文件完整保留、不删除**，只是不进入统计池；如需纳入，把该组织放回 `flux_raw_data.xlsx` 后重跑该批次 |

**批次技术链兼容性**：每次批次运行会写出 `run_provenance.json`，记录示踪剂、
反应集合哈希（代谢网络指纹）、提取方法、检测仪器与拟合迭代次数。整合前
`provenance.assert_batch_compatibility()` 比对各批次指纹：示踪剂 / 网络 /
提取方法不一致属**关键**不一致（合并前需人工确认），仪器与迭代次数差异仅提示。
由 `integration.provenance_mismatch_policy`（`warn` 默认 | `fail` | `ignore`）控制。

**组织命名归一**：不同批次可能用不同名字指同一生理隔室（本项目 `data_20260312`
记作 `serum`、`data_20260319_new` 记作 `plasma`，实为同一份循环血样本）。在
`experiment.yaml` 的 `integration.tissue_aliases` 里声明：

```yaml
integration:
  tissue_aliases:
    serum: plasma
```

别名在**载入后、分组前**统一重写 `tissue` 与 `flux_id`（原名保留在
`tissue_original` 列），Meta 分析与桥接 QC 都会套用。不归一则
`serum--LDH_e` 与 `plasma--LDH_e` 永远无法对齐，该隔室在两处都白白丢失。

该表与转录组的 `joint_design_report.csv` 作用相同。若各批次组织集完全
不重叠，日志会显式警告「这些批次是互补的组织板块，不是彼此的技术重复」。

**FDR 按路线分池**：真正合并过的 pooled effect 与仅单批次直通的结果是
两类强度不同的证据，混在同一个 BH 池里会让后者搭便车。`pooled_fdr`
在每个 `integration_route` 内部单独计算。

**跨批次证据分两档**：

| 列 | 判定条件 |
|---|---|
| `is_cross_batch_supported` | `route == random_effects_meta` **且** pooled 显著 **且** `se_quality` 不是 `imputed_from_floor_pvalue`/`unavailable` **且** `above_noise_floor_p95` |
| `cross_batch_direction_agree` | 仅要求 `route == random_effects_meta` 且两批次效应同向 |

`se_quality` 按 Meta 分析惯例（Cochrane Handbook §6.3）分档：从 p 值反解 SE
是被认可的插补，**前提是该 p 值是真实的检验结果**。当 p 值压在秩检验的样本量
下界上（n=1 vs 3 时恒为 0.5），它不含精度信息，反解出的 SE 只是 |效应量| 的
常数倍——这种「方差」不能支撑显著性主张，单独标为
`imputed_from_floor_pvalue`。

> **注意**：单批次直通行的 `pooled_pvalue` 由 `|log2FC| / SE` 算出，而该 SE
> 又是从原始 p 值反解的——正向算回去必然复现原值。它被标注为
> `pvalue_origin = single_batch_original`，不要当作合并证据。

#### 6.2 桥接样本与技术噪声地板

当同一只动物的同一组织在两个批次里各被拟合一次（不同上机日期、不同
仪器，或仅仅是两次独立的模型运行），两次结果之差**不含任何生物学信号**
——它就是跨批次比较的**技术噪声地板**。`bridging_qc.py` 自动从各批次的
`sample_net_flux.csv` 中找出重复出现的 (sample_id, flux_id) 并量化：

> **前提：`sample_id` 必须能跨批次唯一标识同一只动物。**
> 若上游 `sheet_information.complete_result_name` 写成
> `glucose__male__hfd_sham_8w__1` 这类「组名+序号」标签（本项目目前如此），
> 两个批次里的同名样本未必是同一只鼠。`classify_sample_identity()` 会自动
> 识别并把结果标注为 `upper_bound_confounded_with_biology`——那是
> 「技术重现性 + 生物学个体差异 + 可能的配对错误」之和的上界，
> **不是**技术噪声地板。可用 `integration.bridging.sample_id_is_animal_id`
> 显式覆盖自动判别。

| 输出 | 说明 |
|---|---|
| `integration/bridging_qc/<protocol>/bridge_sample_overview.csv` | 哪些样本出现在多个批次中 |
| `integration/bridging_qc/<protocol>/bridge_flux_pairs.csv` | 每个桥接对的两次估计、差值与对数比 |
| `integration/bridging_qc/<protocol>/bridge_noise_floor.csv` | 每个批次对的方向一致率与噪声地板 |

噪声地板取同号且量级足够的桥接对上 |log2 比值| 的 95 分位数
（ISO 5725 意义上的重复性临界差）。**小于该值的跨批次效应量应视为落在
技术噪声内。** 桥接 QC 在 Meta 分析**之前**运行，得到的地板会写入整合结果的
`noise_floor_log2fc_p95` / `noise_floor_log2fc_median` 列，并以
`above_noise_floor_p95` / `above_noise_floor_median` 标注每条通量是否越过门槛；
P95 门槛已纳入 `is_cross_batch_supported` 的判定。

若**中位**技术差异就已超过实验配置的生物学效应阈值，说明拟合重现性本身
是瓶颈——此时应优先排查模型收敛性与测量重复性，而不是继续在结果层面
做统计。报告会在这种情况下给出显式提示。
对数比只在两次估计同号、且较小绝对值达到中位绝对通量的
`bridging.magnitude_fraction`（默认 5%）时计算——净通量在零附近时
`log2(a/b)` 会爆炸，那属于数值噪声。被排除的对在
`exclusion_reason` 列单独计数，不会被悄悄丢弃。

若实验本身没有刻意留桥接样本，任何被两个批次共同覆盖的组织都会被
自动识别出来充当这一角色。**建议后续批次固定纳入 2–3 个桥接样本**
（同一动物分装或混合 QC 池），否则跨批次差异中的技术成分永远无法与
生物学成分分离。

### 方法 7: Python 编程调用

```python
# 方式 A: 使用 YAML 配置
from omicsweave.mfa.mfa_experiment_config import load_mfa_config
from omicsweave.mfa.run_mfa_analysis import run_from_mfa_config

cfg = load_mfa_config(experiment_id="20250425_cimr_wt_mice_dio")
run_from_mfa_config(cfg)

# 方式 B: 逐步调用（高级用法）
from omicsweave.mfa.mfa_analysis import MFAAnalyser

analyser = MFAAnalyser(output_dir="path/to/output")
analyser.load_flux_data("path/to/flux_raw_data.xlsx")
analyser.load_mid_data("path/to/mid_raw_data.xlsx")

# 可单独执行任一分析步骤
analyser.quality_control()
analyser.flux_pca()
analyser.differential_flux(comparisons=[("hfd", "lfd")])
analyser.flux_heatmap()
analyser.pathway_flux_comparison()
analyser.tissue_flux_overview()
analyser.mid_analysis()

# 手动触发绘图
analyser._generate_all_plots()
analyser.generate_report()
```

---

## 配置参数一览

### YAML 配置文件

YAML 配置由 `mfa_experiment_config.py` 自动识别两种布局：

- 旧版：`experiments/<experiment_id>/mfa_data/config.yaml`
- 新版：`experiments/<experiment_id>/mfa_data/experiment.yaml` + `<protocol>/<batch>/config.yaml`

#### config.yaml

```yaml
experiment_id: 20250425_cimr_wt_mice_dio   # 实验唯一标识
omics: mfa                                   # 组学类型标识

paths:
  flux_file: "flux_raw_data.xlsx"            # 通量数据文件（相对于 config.yaml 的路径）
  mid_file: "mid_raw_data.xlsx"              # MID 数据文件（相对于 config.yaml 的路径）
  output_dir: "analysis_results/${experiment_id}/mfa_results"  # 建议填 mfa_results 根目录；程序自动扩展到 runs/<protocol>/<batch>

# 组别映射文件（默认 groups.yaml，与 config.yaml 同目录）
groups_file: "groups.yaml"

# 可选：显式指定 protocol/batch（未配置时会从 config.yaml 路径推断）
isotope_protocol:
  protocol_id: "u-13c_glucose_infusion"
  batch_id: "data_20251125"

# 可选：显式指定比较组；不配置时自动生成所有两两比较
comparisons:
  - treat: hfd
    control: lfd
  - treat: hsd
    control: lfd

# 可选：统计方法、阈值和数据前处理参数
analysis:
  preprocessing:
    use_net_flux: true          # true: forward - reverse；false: 仅使用正向通量
    iqr_filter: true            # 是否对每个样本多次拟合值做 IQR 过滤
    zero_threshold: 1.0e-10     # log2FC 计算中的近零阈值
    scale_method: zscore        # 热图标准化：zscore | none
    missing_value: drop         # 热图缺失值处理：drop | median_fill

  differential_flux:
    statistic: auto             # auto | mannwhitney | welch_ttest | effect_size_only
                                # auto: p 值/FDR 下界不可达时自动回退 effect_size_only
    multiple_testing: fdr_bh    # fdr_bh | none
    pvalue_threshold: 0.05
    fdr_threshold: 0.10
    log2fc_threshold: 0.5
    cohen_d_threshold: 0.5
    hedges_g_threshold: 0.8     # 未配置时沿用 cohen_d_threshold
    require_identifiable: false # true: poorly_constrained 通量强制置为不显著
    stability:                  # 可选：组内 bootstrap + 标签置换校准
      enabled: false
      use_for_selection: false  # true 时并入 is_significant 判定
      n_bootstrap: 2000
      n_permutations: 5000
      alpha: 0.05
      min_selection_probability: 0.80
      min_direction_stability: 0.80

  annotation:
    require_standard_reaction_ids: true  # 强制输出标准注释列（默认 true）

  heatmap:
    top_n: 50
    pvalue_threshold: 0.05
    rank_by: auto               # auto | abs_log2fc | abs_cohen_d | pvalue

  volcano:
    log2fc_threshold: 0.5
    pvalue_threshold: 0.05
    cohen_d_threshold: 0.5
    label_top_n: 10
    symmetric_x: true         # x 轴强制对称于 0（默认 true）
    x_clip_quantile: 0.98     # x 轴范围取 |log2FC| 的该分位数
    x_max_abs: 10.0           # x 轴绝对上限（2^10 = 1024 倍）
    show_infinite: true       # 越界/无穷点钉在边界而非丢弃
    axis_pad: 0.08
    # xlim: [-6, 6]           # 可选：显式固定范围，便于跨图比较
    # ylim: [0, 5]
```

#### 配置发现契约（新增批次必读）

protocol / batch 的发现**完全由 `experiment.yaml: tracer_protocols` 驱动**，
而不是扫描磁盘。目录存在但未在 YAML 中声明时，既不会报错也不会被分析。
为避免这种静默遗漏，`discover_experiment_protocols()` 会对未声明的目录打印
`[WARN]`。

因此新增一个批次需要三步，缺一不可：

1. 建立 `mfa_data/<protocol>/<batch>/`，放入 `flux_raw_data.xlsx`、
   `mid_raw_data.xlsx`、`batch.yaml`、`config.yaml`、`groups.yaml`；
2. 把批次 ID 加入 `tracer_protocols.<protocol>.batches`；
3. 该 protocol 的 `enabled` 置为 `true`。

`has_data` 的判定会沿 `config.yaml → batch.yaml → data_files.flux_file` 逐级
解析并检查文件是否真实存在，所以「已声明但没数据」的批次会被报告而非静默跳过。

`defaults.active_batch` 可写单值或列表；取第一个 `config.yaml` 存在的批次作为
`--experiment` 单批次运行的默认目标。

**路径解析规则**：
- 以 `\\` 开头的 UNC 路径或 Windows 盘符路径（`C:\...`）视为绝对路径
- 其他路径优先按仓库根目录解析，若文件不存在则按 config.yaml 所在目录解析
- 支持 `${experiment_id}`、`${protocol_id}`、`${batch_id}`、`${output_root}`、`${output_dir}` 变量替换

| 字段 | 类型 | 必需 | 说明 |
|---|---|---|---|
| `experiment_id` | str | 是 | 实验唯一标识，可由目录名自动推断 |
| `paths.flux_file` | str | 否 | flux_raw_data.xlsx 路径（支持绝对/相对路径） |
| `paths.mid_file` | str | 否 | mid_raw_data.xlsx 路径（支持绝对/相对路径） |
| `paths.output_dir` | str | 是 | 输出根目录（建议 `analysis_results/<id>/mfa_results`，程序会自动分层到 `runs/<protocol>/<batch>`） |
| `groups_file` | str | 否 | 组别映射文件，默认 `groups.yaml` |
| `isotope_protocol.protocol_id` | str | 否 | 当前运行的同位素 protocol 标识（用于输出分层） |
| `isotope_protocol.batch_id` | str | 否 | 当前运行批次标识（用于输出分层） |
| `comparisons` | list | 否 | 比较对列表，支持 `{treat, control}` 或 `[treat, control]`，未配置则自动从 conditions 生成两两比较 |
| `analysis.preprocessing` | dict | 否 | 数据前处理参数：净通量、IQR 过滤、近零阈值、热图标准化/缺失值处理 |
| `analysis.differential_flux` | dict | 否 | 差异通量统计方法和显著性阈值 |

**分析参数的三级覆盖链**：

```
DEFAULT_ANALYSIS_CONFIG（代码内置）
  ← experiment.yaml : analysis_defaults
    ← <protocol>/protocol.yaml : analysis_defaults
      ← <protocol>/<batch>/config.yaml : analysis
```

逐层 deep-merge，越靠后优先级越高。批次 config 写 `analysis: {}` 时，
实际生效的是 experiment/protocol 级的 `analysis_defaults`。
| `analysis.annotation` | dict | 否 | 差异通量结果注释策略，`require_standard_reaction_ids` 默认 `true`（强制包含 `kegg_reaction/kegg_pathway/reactome_id/ec_number/description` 列） |
| `analysis.heatmap` | dict | 否 | 热图筛选数量、P 值阈值和排序指标 |
| `analysis.volcano` | dict | 否 | 火山图阈值和标注数量 |
| `condition_colors` | dict | 否 | 条件→颜色映射，未配置则从 groups.yaml 的 conditions.color 自动提取 |
| `tissue_colors` | dict | 否 | 组织→颜色映射，未配置则使用内置默认配色 |
| `pathway_groups` | dict | 否 | 自定义通路分组，未配置则使用 10 个预定义通路 |

#### experiment.yaml（整合与整合绘图）

```yaml
integration:
  enabled: true
  primary_strategy: "within_protocol"
  allow_cross_protocol: false      # 跨示踪剂整合默认关闭
  provenance_mismatch_policy: warn # warn | fail | ignore

  # 同一生理隔室在不同批次里的不同叫法 → 规范名
  tissue_aliases:
    serum: plasma

  bridging:
    sample_id_is_animal_id: true   # 同名样本确为同一只动物时才是技术噪声地板
    magnitude_fraction: 0.05       # 对数比的量级门槛（相对中位绝对通量）
    floor_quantile: 0.95           # 噪声地板取的分位数

  within_protocol:
    method: "random_effects_meta"
    min_batches: 2
    heterogeneity_metric: "i2"
    heterogeneity_threshold: 0.50

  cross_protocol:
    method: "directional_consensus"
    min_protocols: 2
    require_same_direction: true
    min_protocols_significant: 2   # 仅方向一致不算证据
    exclude_poorly_constrained: true

  plotting:
    enabled: true
    output_subdir: "02_result_plot"
    top_n: 20
    max_label_len: 46
    plots:
      meta_forest: true
      heterogeneity_hist: true
      consensus_bar: true
```

`integration.plotting` 字段说明：

| 字段 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `enabled` | bool | `false` | 是否在整合结束后自动绘图 |
| `output_subdir` | str | `02_result_plot` | 相对 `mfa_results/integration/` 的输出子目录 |
| `top_n` | int | `20` | 森林图和一致性条形图展示的 Top 条目数 |
| `max_label_len` | int | `46` | 反应名标签最大字符数（超长会截断） |
| `plots.meta_forest` | bool | `true` | 生成每个 protocol/比较组的 pooled effect 森林图 |
| `plots.heterogeneity_hist` | bool | `true` | 生成每个 protocol 的 I² 异质性分布图 |
| `plots.consensus_bar` | bool | `true` | 生成跨 protocol 一致性 Top 条形图及汇总图 |

#### protocol.yaml（示踪剂与灌注 QC）

```yaml
protocol_id: "u-13c_glucose_infusion"

tracer:                         # 进入 run_provenance.json，用于跨批次兼容性判定
  type: "glucose"
  isotope: "U-13C6"
  isotopic_enrichment_atom_percent: 99

infusion:
  route: "jugular_vein"
  duration_min: 180
  rate_value: 0.05
  rate_unit: "ul/g/min"

infusion_qc_defaults:           # 灌注稳态 QC 阈值（已生效）
  enabled: true
  steady_state_window_min: [30, 180]
  glucose_cv_max: 0.15
  glucose_range_max_mmol_l: 2.0
  # glucose_drift_max_mmol_l_per_h: 1.0   # 可选，配置后才参与判定
  on_fail: "flag"               # flag | exclude

analysis_defaults:              # 覆盖 experiment.yaml，被 batch config 覆盖
  differential_flux:
    statistic: auto
```

> `on_fail: exclude` 目前降级为 `flag`：血糖表用动物编号、MFA 结果用
> 「组名+序号」，两者无法逐一对应，无法做样本级排除。

#### batch.yaml（数据文件与检测元数据）

```yaml
batch_id: "data_20260319_new"
protocol_id: "u-13c_glucose_infusion"

data_files:
  flux_file: "flux_raw_data.xlsx"
  mid_file: "mid_raw_data.xlsx"
  glycemia_record_file: "26.1.2-mice_blood_glucose_record_glucose_infusion.xlsx"

batch_metadata:                 # 进入 run_provenance.json
  detection_device: "sciex_7500"          # 提示级差异
  metabolites_extraction_protocal: "hydrophilic_metabolites"  # 关键字段
  detection_division: "cimr_massspec_core"
  detection_date: "2026-03-19"

batch_qc:
  inherit_protocol_glycemia_qc: true
  on_qc_fail: "flag"            # 映射为 infusion_qc.on_fail
```

#### groups.yaml

```yaml
# 实验条件定义
conditions:
  hfd:
    display: "HFD"                # 显示名称（用于图表标签）
    diet: "D12492"                # 饲料编号
    diet_desc: "High-fat diet DIO 18-20 week"  # 饲料描述
    color: "#FF8C00"              # 绘图颜色
  lfd:
    display: "LFD"
    diet: "D12450J"
    diet_desc: "Low-fat diet control 18-20 week"
    color: "#2E8B57"
  hsd:
    display: "HSD"
    diet: "AIN-76A"
    diet_desc: "High-sugar purified diet 18-20 week"
    color: "#8B008B"

# MFA 样本名 → 条件代码映射
# 样本名对应 flux/mid Excel 的 sheet_information.complete_result_name
sample_groups:
  "glucose__male__hfd__1": hfd
  "glucose__male__hfd__2": hfd
  "glucose__male__hfd__3": hfd
  "glucose__male__hfd__4": hfd
  "glucose__male__lfd__1": lfd
  "glucose__male__lfd__2": lfd
  "glucose__male__lfd__3": lfd
  "glucose__male__hsd__1": hsd
  "glucose__male__hsd__2": hsd
  "glucose__male__hsd__3": hsd
```

---

## 输出文件清单

### 数据表（`runs/<protocol>/<batch>/01_df_results/`）

**跨组织分析（始终保存在 `all_tissues/`）**：

| 文件路径 | 说明 |
|---|---|
| `run_provenance.json`（批次根目录） | **技术链指纹**：示踪剂、反应集合哈希、提取方法、仪器、迭代次数 |
| `01_df_results/all_tissues/qc/qc_fitting_loss.csv` | 各样本拟合质量（loss 值统计） |
| `01_df_results/all_tissues/qc/infusion_steady_state_qc.csv` | **灌注稳态 QC**：每动物 CV / 极差 / 漂移速率 / 通过与否 |
| `01_df_results/all_tissues/qc/infusion_glycemia_timecourse.csv` | 血糖时间曲线长表 |
| `01_df_results/all_tissues/qc/infusion_steady_state_by_condition.csv` | 稳态 QC 按条件汇总 |
| `01_df_results/all_tissues/flux_matrix/sample_net_flux.csv` | **样本级净通量矩阵**（每动物每反应一行，含拟合离散度与 `fit_rel_iqr`） |
| `01_df_results/all_tissues/flux_matrix/flux_identifiability.csv` | **通量可辨识性**：每通量的相对 IQR 与分档 |
| `01_df_results/all_tissues/mid/mid_fit_residuals.csv` | **MID 拟合残差**：每 (样本 × 代谢物) 的 SSR / RMSE / 最大残差 |
| `01_df_results/all_tissues/mid/mid_fit_quality_by_metabolite.csv` | 按代谢物汇总的失拟情况 |
| `01_df_results/all_tissues/mid/mid_residuals_long.csv` | 逐同位素体的实测 vs 模拟明细 |
| `01_df_results/all_tissues/pca/flux_pca_scores.csv` | PCA 得分矩阵 + 条件标签 + 组织 |
| `01_df_results/all_tissues/pca/flux_pca_loadings_{tissue}.csv` | 各组织 PCA 载荷矩阵 |
| `01_df_results/all_tissues/overview/tissue_flux_overview.csv` | 各组织总通量和平均通量 |
| `01_df_results/all_tissues/mid/mid_summary.csv` | 完整 MID 汇总（代谢物 × MID_index × 条件） |
| `01_df_results/all_tissues/mid/mid_total_labeling.csv` | 各代谢物总标记度 |
| `01_df_results/all_tissues/pathway_flux/pathway_flux_reaction_stats.csv` | 跨组织通路反应统计总表 |
| `01_df_results/all_tissues/pathway_flux/pathway_flux_summary_stats.csv` | 跨组织通路汇总统计总表 |

**组织特异性分析（按 `<tissue>/<category>/` 分层存储）**：

| 文件路径 | 说明 |
|---|---|
| `01_df_results/<tissue>/differential_flux/differential_flux_{A}_vs_{B}.csv` | 差异通量完整结果表（含通路注释列） |
| `01_df_results/<tissue>/heatmap/flux_heatmap_raw_{tissue}_{comp}.csv` | 热图原始数据矩阵 |
| `01_df_results/<tissue>/heatmap/flux_heatmap_stats_{tissue}_{comp}.csv` | 热图显著通量统计 |
| `01_df_results/<tissue>/pathway_flux/pathway_flux_samples.csv` | 通路 × 条件 × 样本通量数据 |
| `01_df_results/<tissue>/pathway_flux/pathway_flux_reaction_stats.csv` | 通路 × 反应 × 条件统计（mean/std/n/n_points） |
| `01_df_results/<tissue>/pathway_flux/pathway_flux_summary_stats.csv` | 通路 × 条件汇总统计（mean/std/n/n_points/n_reactions） |

### 可视化图表（`runs/<protocol>/<batch>/02_result_plot/`）

**跨组织图表（始终保存在 `all_tissues/` 层级）**：

#### 02_result_plot/all_tissues/0_qc/（质量控制）

| 文件名 | 说明 |
|---|---|
| `qc_fitting_loss.png` | 拟合质量柱状图 + 直方图 |
| `qc_infusion_glycemia.png` | 灌注期血糖时间曲线 + 稳态窗口（实线通过 / 虚线未通过） |
| `qc_mid_fit_residuals.png` | 实测 vs 模拟 MID 散点 + 失拟最严重的代谢物 |
| `qc_flux_identifiability.png` | 通量可辨识性分布与分档计数 |

#### 02_result_plot/all_tissues/3_pca_clustering/（PCA 聚类）

| 文件名 | 说明 |
|---|---|
| `flux_pca_{tissue}.png` | 各组织 PCA 散点图 + Scree plot |

#### 02_result_plot/all_tissues/6_tissue_overview/（组织概览）

| 文件名 | 说明 |
|---|---|
| `tissue_flux_overview.png` | 组织总通量/平均通量对比图 |

#### 02_result_plot/all_tissues/7_mid_analysis/（MID 标记分析）

| 文件名 | 说明 |
|---|---|
| `mid_total_labeling.png` | Top 差异代谢物标记度柱状图 |
| `mid_distributions.png` | 高标记代谢物的 MID 分布图 |

**组织特异性图表（按 `<tissue>/<category>/` 分层存储）**：

#### 02_result_plot/<tissue>/1_volcano/（差异通量火山图）

| 文件名 | 说明 |
|---|---|
| `volcano_{comp}.png` | 单组织差异通量火山图（x 对称于 0、y 自 0 起、越界点钉在边界） |
| `tissue_diff_summary_{comp}.png` | 单组织差异通量上调/下调汇总 |

#### 02_result_plot/<tissue>/2_heatmap/（差异通量聚类热图）

| 文件名 | 说明 |
|---|---|
| `flux_heatmap_{tissue}_{comp}.png` | 聚类热图（seaborn clustermap） |

#### 02_result_plot/<tissue>/4_pathway_flux/（通路通量）

| 文件名 | 说明 |
|---|---|
| `pathway_flux_{tissue}.png` | 各组织代谢通路通量分面图 |
| `pathway_flux_{tissue}_individual.png` | 全部反应独立分布总览图 |
| `by_pathway/<pathway_name>/pathway_flux_{tissue}_{pathway_name}.png` | 单一通路总反应箱线图 |
| `by_pathway/<pathway_name>/pathway_flux_{tissue}_{pathway_name}_individual.png` | 单一通路反应拆分箱线图 |

### 整合输出（`mfa_results/integration/`，`--integrate` 产生）

| 文件路径 | 说明 |
|---|---|
| `tracer_harmonized/<protocol>/integration_design_report_<protocol>.csv` | **路线判定表——先读这张**：每个 (comparison × tissue) 走哪条整合路线及原因 |
| `tracer_harmonized/<protocol>/differential_flux_meta_<comp>.csv` | 合并结果（含 `integration_route` / `se_quality` / `above_noise_floor_*` / `is_cross_batch_supported`） |
| `tracer_harmonized/<protocol>/meta_analysis_summary_<protocol>.csv` | 各比较的合并统计汇总 |
| `bridging_qc/<protocol>/bridge_noise_floor.csv` | **技术噪声地板**：每批次对的方向一致率、|log2 比| 中位数与 P95 |
| `bridging_qc/<protocol>/bridge_flux_pairs.csv` | 每个桥接对的两次估计与对数比 |
| `bridging_qc/<protocol>/bridge_sample_overview.csv` | 哪些样本出现在多个批次中 |
| `cross_tracer_scoring/directional_consensus_<comp>.csv` | 跨 protocol 方向印证（含证据门槛列） |
| `cross_tracer_scoring/protocol_comparison_summary.csv` | 跨 protocol 汇总 |
| `integration_report.md` | 整合分析报告（含路线判定与噪声地板解读） |
| `02_result_plot/` | 森林图 / 异质性分布 / 一致性条形图 |

### 解读报告 (03_interpretation_report/)

| 文件名 | 说明 |
|---|---|
| `index.md` | 目录页（实验概况 + 文件清单 + 各组织报告链接） |
| `{tissue}.md` | 每组织独立解读报告（差异通量、TCA、通路、MID 摘要 + 解读模板） |

---

## 统计方法说明

### 推断单位：生物学重复

每只动物先被汇总为**一个**净通量估计——其 IQR 过滤后拟合点的**中位数**，
统计检验只在这些动物级估计之间进行。

模型拟合的多次迭代（本项目各批次为 24–150 次不等）描述的是**同一只动物上
参数的不确定度**，不是独立的生物学观测。把它们当作重复会把样本量虚增两个
数量级，产生大量假显著。

### IQR 过滤

对每个样本的全部拟合迭代，按排序位置保留中间 50%（去除上下各 25%），
排除极端拟合值的影响。过滤在**每个样本内部**进行，不跨样本合并。

### 净通量计算

对可逆反应（同时有正向 `reaction` 与逆向 `reaction__R`），净通量 =
forward − reverse，以反映反应的实际方向与强度。

> **注意**：净通量跨零时 `log2FC` 不可靠——通量从 +5 变到 −5 是**方向反转**，
> 不是倍数变化。可逆 / 净反应优先看 `hedges_g`。

### 秩检验的样本量下界

Mann-Whitney 精确检验的最小双侧 p 值为 `2/C(n1+n2, n1)`，由样本量唯一决定：

| n1 vs n2 | 最小可达双侧 p | p<0.05 可能吗 |
|---|---|---|
| 1 vs 1 | 1.0 | 否 |
| 2 vs 2 | 0.333 | 否 |
| 3 vs 3 | **0.1** | **否** |
| 3 vs 4 | **0.057** | **否** |
| 4 vs 4 | 0.029 | 是 |
| 5 vs 5 | 0.0079 | 是 |

BH 校正后还有第二道下界（`min_achievable_fdr`）：即使未校正 p 可达，
在数百条通量上做校正后仍可能整体不可达。

实现见 `omicsweave/common/nonparametric.py`；结果列
`min_achievable_pvalue` / `min_achievable_fdr` / `rank_test_underpowered` /
`fdr_unreachable`。

### 统计路径与证据等级

`statistic: auto`（默认）在两条 p 值路径都不可达时自动回退到效应量路径：

| evidence_level | 触发条件 | 判定依据 |
|---|---|---|
| `statistical_test` | p 值与 FDR 均可达 | Mann-Whitney p + BH-FDR + \|log2FC\| |
| `effect_size_only` | 下界超阈值，但两组 n≥2 | \|Hedges' g\| + \|log2FC\| |
| `insufficient` | 任一组 n<2 | 无法推断，`is_significant` 强制为 False |

### 效应量：Hedges' g

Cohen's d 在小样本下系统性高估效应，Hedges' g 乘以校正因子
`1 − 3/(4·df − 1)`。两者都输出，判定以 `hedges_g` 为主。
任一组 n<2 时组内方差无定义，两者均为 NaN。

### Welch's t-test

不假设两组方差相等的 t 检验变体，作为参考统计量并列输出。

### Benjamini-Hochberg FDR

对同一比较内的全部通量做多重检验校正。整合阶段的 `pooled_fdr` **按整合路线
分池**计算——真正合并过的 pooled effect 与仅单批次直通的结果是两类强度不同的
证据，混在同一个 BH 池里会让后者搭便车。

### 通量可辨识性

`fit_rel_iqr = (p75 − p25) / |net_flux|`。拟合迭代不是生物学重复，但确实刻画了
「在这只动物的数据下，这条通量被约束得有多紧」。相对 IQR 接近或超过 1，
意味着该通量几乎不被数据决定，其组间差异不携带信息。

### 稳定性校准（可选）

`analysis.differential_flux.stability.enabled: true` 时，对每条通量执行组内
bootstrap 与标签置换，输出 `selection_probability` / `direction_stability` /
`permutation_pvalue`。任一组少于 3 个重复时只返回描述性效应。
`stability.use_for_selection: true` 可把它并入显著性判定。

### 随机效应 Meta 分析

跨批次合并使用 DerSimonian-Laird 估计量，效应量尺度为 `log2FC`，
输出 `I²`（批次间异质性占比）与 `τ²`（批次间方差）。
SE 来源分档见 `se_quality`：从**压在样本量下界上**的 p 值反解出的 SE
不含精度信息，按 Meta 分析惯例不支撑显著性主张。

### 火山图坐标轴约定

火山图的坐标轴**不采用自适应范围**，原因有二：

**1. x 轴必须对称于 0。** 火山图的全部意义在于比较上调与下调；自适应范围会把
「一侧数据碰巧更极端」渲染成「一侧变化更多」的视觉印象。实测本项目 175 张
火山图，**73 张（42%）的 x 轴正负范围偏离对称超过 1.5 倍**，最极端达 11.7 倍。
对称之后左右两半面积可直接比较。

**2. 范围按分位数取，并设绝对上限。** 净通量的分母接近零时 `log2FC` 会冲到
±20 以上，那是数值假象而非更强的生物学变化。全库 13603 个有限 `|log2FC|` 的
分布印证了这一点：

| 分位 | \|log2FC\| | 倍数 |
|---|---|---|
| 50% | 1.22 | 2.3× |
| 90% | 5.37 | 41× |
| 98% | 9.15 | 567× |
| 99% | **21.71** | **3.4×10⁶** |
| 100% | 41.58 | 3.3×10¹² |

98 分位到 99 分位之间存在明显断裂——数值假象就从那里开始。因此 x 轴范围取
`min(quantile(|log2FC|, 0.98), 10)`，同时保证判定阈值线始终可见。上限 10
（1024 倍）使 1.49% 的点越界，与 `±inf` 点占比（1.5%）吻合。

**3. y 轴自 0 起。** `-log10(p)` 的 0 对应 p=1，`|Hedges' g|` 的 0 对应无效应，
都是天然基线；浮动基线会把微小差异放大成视觉上的显著分离。

**4. 越界点钉在边界，不丢弃。** 早期实现直接过滤掉 `log2FC = ±inf` 的通量
（全库 206 个，占 1.5%）——那正是「一组有通量、另一组几乎为零」的开/关式变化，
往往最值得看。现在它们与其他越界点一起以**空心三角**画在 x 轴边界上，
并在图例中计数（`off-scale |log2FC| (n)`）。

如需跨组织 / 跨比较直接对比，可在 `analysis.volcano` 中显式设定 `xlim` / `ylim`
固定所有图的范围；默认按图自适应（但仍对称、仍自 0 起），以免离群值把某些图
压得无法阅读。

### PCA

对 Z-score 标准化后的通量矩阵按组织分别做特征值分解，提取最大方差方向。

---

## 添加新实验 / 新批次的步骤

### 新增一个批次（最常见）

发现流程**完全由 `experiment.yaml: tracer_protocols` 驱动，不扫描磁盘**。
以下三步缺一不可，只建目录不登记等于什么都没做：

1. **建目录放数据**
   ```
   experiments/{experiment_id}/mfa_data/{protocol_id}/{batch_id}/
     flux_raw_data.xlsx
     mid_raw_data.xlsx
     <血糖记录>.xlsx          # 可选，用于灌注稳态 QC
     config.yaml  batch.yaml  groups.yaml
   ```
   三个 YAML 从同 protocol 的既有批次复制后修改即可。

2. **在 `experiment.yaml` 登记批次**
   ```yaml
   tracer_protocols:
     {protocol_id}:
       enabled: true
       batches:
         - {既有批次}
         - {batch_id}          # ← 新增
   ```

3. **确认该 protocol 的 `enabled: true`**

未登记的目录会在运行时打印 `[WARN] ... 未声明的批次目录（不会被分析）`。

### 各 YAML 需要修改的字段

| 文件 | 关键字段 |
|---|---|
| `batch.yaml` | `batch_id`、`data_files`（flux/mid/glycemia）、`batch_metadata`（仪器、提取方法、检测日期）、`batch_qc` |
| `config.yaml` | `experiment_id`、`comparisons`（treat/control 顺序决定输出名）、`analysis`（留 `{}` 即继承上级） |
| `groups.yaml` | `conditions`（display / diet / color）、`sample_groups`（样本名→条件，样本名须与 Excel `sheet_information.complete_result_name` 一致） |

### 新增一个示踪剂 protocol

在 `mfa_data/` 下建 `{protocol_id}/protocol.yaml`（声明 `tracer`、`infusion`、
`sampling`、`infusion_qc_defaults`），并在 `experiment.yaml` 的
`tracer_protocols` 中登记。跨 protocol 整合默认关闭，确认可比后再把
`integration.allow_cross_protocol` 置为 `true`。

### 新增一个实验

复制既有实验的 `mfa_data/` 目录结构，修改 `experiment.yaml` 的
`experiment_id`、`defaults`、`tracer_protocols`，再按上面的步骤配置
protocol 与 batch。

### 运行

```bash
omicsweave mfa run --experiment {experiment_id} --all-batches --with-plots
```

```bash
omicsweave mfa run --experiment {experiment_id} --integrate
```

### 可选：自定义通路分组

如果新实验的代谢网络模型与默认 10 个通路分组不同，在 `config.yaml` 中添加 `pathway_groups` 字段：

```yaml
pathway_groups:
  Glycolysis:
    - HEX_c
    - PGI_c
    - PFK_c
  Custom_Pathway:
    - REACTION1_c
    - REACTION2_m
```

---

## 通路注释映射

### 概述

MFA 分析模块通过 `omicsweave/common/reaction_pathway_mapping.py` 为所有自定义反应 ID 提供公共数据库通路注释。该模块是项目级别的单一数据源（single source of truth），被 MFA 管线、多组学集成和任何需要反应-通路映射的工具共享使用。

### 映射表内容

`REACTION_REGISTRY` 字典包含 **137 个 MFA 反应 ID**，每个条目提供以下注释字段：

| 字段 | 示例 | 说明 |
|---|---|---|
| `kegg_reaction` | `R00771` | KEGG Reaction ID |
| `kegg_pathway` | `map00010` | KEGG Pathway map ID |
| `reactome_id` | `R-MMU-70493` | Reactome stable ID（小鼠前缀 `R-MMU-`） |
| `ec_number` | `2.7.1.1` | EC 酶编号 |
| `description` | `Hexokinase: Glc + ATP -> G6P + ADP` | 英文描述 |

### 覆盖的代谢通路

| KEGG 通路 | 反应数 | 示例反应 |
|---|---|---|
| `map00010` Glycolysis / Gluconeogenesis | 20+ | HEX_c, PFK_c, PK_c, PEPCK_c |
| `map00020` TCA Cycle | 15+ | CS_m, ICDH_m, MDH_m |
| `map00030` Pentose Phosphate | 8 | G6PDH_c, TKT1_c, RPE_c |
| `map00071` Fatty Acid Degradation | 6 | FAO_m, CPT1, CARNITINE_shuttle |
| `map00061` Fatty Acid Biosynthesis | 3 | ACC_c, FAS_c, ACOT_c |
| `map00072` Ketone Body Metabolism | 7 | BHBDH_m, SCOT_m |
| `map00251` Glutamate/Glutamine | 7 | GLUD_m, GLNS_c, GLND_m |
| `map00250` Ala/Asp/Glutamate | 8 | ASPTA_m, GPT_c |
| `map00561` Glycerolipid Metabolism | 6 | GLYCK_c, G3PDH_c |
| `map00500` Starch/Sucrose | 3 | PGMT_c, GLYCOGENOLYSIS |
| `map00620` Pyruvate Metabolism | 5 | PDH_m, ME_c, PYR_trans |
| `map00330` Arg/Pro Metabolism | 5 | CPS1_m, UREA_cycle |
| Others | 20+ | HMGCR, CK_c, GSH_synthesis, etc. |

### 使用示例

```python
from omicsweave.common.reaction_pathway_mapping import (
    get_kegg_pathway, get_ec_number, annotate_flux_dataframe,
)

# 单个反应查询
pathway = get_kegg_pathway("HEX_c")  # → "map00010"
ec = get_ec_number("CS_m")           # → "2.3.3.1"

# 批量注释 DataFrame
import pandas as pd
df = pd.DataFrame({"reaction": ["HEX_c", "CS_m"]})
annotated = annotate_flux_dataframe(df)
# → df 新增 ec_number, kegg_pathway, reactome_id, description 列
```

### 设计原则

- **单一数据源**：定义在 `omicsweave/common/`，所有模块共享同一份映射
- **自动注释**：MFA 差分通量 CSV 输出时自动 join 注释列，无需手动调用
- **向后兼容**：`omicsweave/multiomics/mappings/reaction_registry.py` 保留为 re-export，旧导入路径仍可用
- **缺失容忍**：聚合反应或转运反应对应的注释字段为 `None`，下游代码应做空值检查

---

## 参考文献

1. Antoniewicz MR. (2015). "Methods and advances in metabolic flux analysis: a mini-review". *J Ind Microbiol Biotechnol*, 42(3), 317-325.
2. Crown SB, Antoniewicz MR. (2013). "Publishing 13C metabolic flux analysis studies: a review and future perspectives". *Metab Eng*, 20, 42-48.
3. Young JD. (2014). "INCA: a computational platform for isotopically non-stationary metabolic flux analysis". *Bioinformatics*, 30(9), 1333-1335.
4. Wiechert W. (2001). "13C metabolic flux analysis". *Metab Eng*, 3(3), 195-206.
5. Benjamini Y, Hochberg Y. (1995). "Controlling the false discovery rate: a practical and powerful approach to multiple testing". *J R Stat Soc Ser B*, 57(1), 289-300.
6. Virtanen P, et al. (2020). "SciPy 1.0: fundamental algorithms for scientific computing in Python". *Nature Methods*, 17(3), 261-272.
