# 普通转录组学 (Bulk RNA-Seq) 分析流程说明书

> **位置**：`omicsweave/transcriptomics/`
> **适用平台**：Windows 10/11 原生（无需 WSL/Linux 服务器）
> **运行环境**：Python ≥ 3.10，安装方式见仓库根目录 README
> **配置驱动**：按 experiment 复用，只改 YAML 配置文件，不改代码

## 按统一项目模型运行（阶段 2）

```bash
omicsweave run --experiment <id> --omics transcriptomics      # 步骤 qc、differential、enrichment
```

读取 `experiment.yaml` 中 `omics.transcriptomics.input` 声明的表格层输入，分组、协变量与比较来自 `subjects.tsv` / `samples.tsv` 和 `design`：

```yaml
input: {level: table, format: counts_matrix, files: [experiments/<id>/gene_counts.tsv]}   # 基因 × 样本计数
# 或已有的 kallisto 定量结果：quant/<样本>/abundance.tsv
input: {level: table, format: kallisto, files: [experiments/<id>/quant], tx2gene: experiments/<id>/tx2gene.tsv}
```

数据中的样本名即 `sample_id`，不同时在 samples.tsv 的 `data_column` 列写出。差异表达沿用本模块的 PyDESeq2（apeGLM 收缩）与秩和检验，
设计为 `~[协变量 +] condition`；同一份计数与旧流程的结果逐值一致。

富集（`enrichment`）：默认仅执行本地 Fisher ORA（上调、下调基因分别检验，背景为有 p 值的基因，不上传 Enrichr）；预排序 GSEA（log2FC × −log10 p，默认 6 个 Enrichr 库合并检验）由 `analysis.run_gsea` 控制，默认 `false`（关闭以消除 ~76%–94% 的置换耗时），可通过配置 `analysis.run_gsea: true` 或命令行 `--gsea` 开启并获得与旧流程逐值一致的结果。基因集库用 `analysis.gene_sets` 指定，可写本地 GMT 路径或 Enrichr 库名；按库名时首次下载后缓存在项目的 `references/gene_sets/` 下。

FASTQ 的质控、修剪与 kallisto 定量（阶段 3）、图件与多批次联合模型，暂时仍用下文的旧流程。结果写入结果包的 `omics/transcriptomics/`，可用 `omicsweave export` 导出为 Excel。

---

## 版本更新说明（2026-08）

本次更新修复了定量与统计层面的若干实质问题，并补齐了多批次合并所需的
前置措施。**已有的 DE 结果需要重跑**（见下方「迁移指引」）。

### 改进措施一览

| # | 改进 | 问题 | 落点 |
|---|---|---|---|
| 1 | 基因层聚合纳入转录本有效长度 | 直接对 `est_counts` 求和，在异构体切换时产生假阳性 | `normalization.py` |
| 2 | Wilcoxon 前做文库归一化 | 在原始 count 上做秩检验，测序深度差异造成全基因 log2FC 偏移 | `run_wilcoxon_de.py` |
| 3 | 比对率落盘 + QC 阈值真正执行 | 阈值「声明但不执行」；比对率只 print 不保存 | `qc_gate.py` |
| 4 | PCA 批次着色 + 主成分协变量关联检验 | 完全没有批次效应诊断 | `plotting/pca.py` |
| 5 | 技术参数记录 + 跨批次一致性断言 | 参考版本/链特异性不同的批次会被静默合并 | `provenance.py` |
| 6 | LFC 收缩（apeGLM） | 低表达基因 log2FC 膨胀，恰好被 `|log2FC|>1` 优先选中 | `rnaseq_pipeline.py` |
| 7 | 跨批次联合模型（merge_runs + `~batch + condition`） | 只有事后 Meta，统计效力偏低 | `joint_model.py` |

### 每项的具体说明

**1. 有效长度校正（lengthScaledTPM）**

旧实现把同一基因所有转录本的 `est_counts` 直接相加。当两组之间发生
**异构体切换**（长短转录本互换）时，基因总 count 会在真实表达量不变的
情况下变化——这正是 tximport 存在的理由。现按 tximport 的
`lengthScaledTPM` 计算：基因 TPM × 基因有效长度（丰度加权），再按列
重标定回原始文库大小。

验证（合成数据，长短转录本有效长度 500 / 2000）：假 log2FC 从 **2.0 降到 <0.1**。

> 注意：`lengthScaledTPM` 会**刻意保留原始文库大小**（供 DESeq2 估计离散度），
> 因此长度校正的效果要在 size factor 归一化之后才体现。

**2. Wilcoxon 文库归一化 + 样本量分辨率上限**

秩检验前先用 median-of-ratios 大小因子归一化。同时新增
`min_achievable_pvalue()`：Mann-Whitney 精确检验的最小双侧 p 为
`2/C(n1+n2, n1)`——**n=3 vs 3 时为 0.1**，经 BH 校正后 `padj<0.05`
在数学上不可达。流程会打印该上限，避免把「0 个显著基因」误读为生物学结论。

**3. QC 门控**

`experiment.yaml` 的 `analysis_defaults.qc` / `.alignment` 与 `batch.yaml`
的 `batch_qc` 现在真正生效。kallisto `run_info.json` 的 `p_pseudoaligned`
落盘到 `03_kallisto_quant/quant_summary.csv`，判定结果写入
`01_qc/qc_gate_report.csv`。`batch_qc.on_qc_fail` 决定处置方式：
`flag`（仅告警）/ `exclude`（移出下游 DE）/ `fail`（中止流程）。

R1/R2 取**较差的一侧**作为样本级指标——门控应由最差的一端决定。

**4. 批次效应诊断**

`05_result_plot/3_pca_clustering/` 新增：
- `pca_global.png` 增加「按批次着色」面板（配置了 `sample_batches` 时）
- `pc_covariate_association.csv` / `.png`：各主成分与 condition/tissue/batch
  的 Kruskal-Wallis 关联检验

PCA 改用 log2-CPM（此前是 log2 原始 count），避免文库大小主导 PC1。

**5. 技术参数与一致性断言**

`batch.yaml` 新增 `technical_params`（建库试剂盒、链特异性、读长、RIN 等）。
每次运行写出 `run_provenance.json`（参考文件指纹、kallisto 版本、链特异性、
聚合方法）。整合前比对各批次，不一致时按
`integration.reference_mismatch_policy`（`fail` 默认 / `warn` / `skip`）处置。

**链特异性现在会传给 kallisto**（`--rf-stranded` / `--fr-stranded`），不再是死配置。

**6. LFC 收缩**

apeGLM 收缩，p 值不变。

> ⚠ **实现要点**：PyDESeq2 0.5.4 的参考水平按字母序决定（`ref_level` 已废弃且无效），
> 因此 `lfc_shrink()` 可用的系数方向可能与 contrast **恰好相反**——实测符号
> 100% 翻转。`_resolve_shrinkage_coeff()` 会识别并取反，并有符号一致性断言兜底。
> **不要移除该断言。**

**7. 跨批次联合模型**

详见下方 [Step 10](#step-10跨批次联合建模joint_modelpy)。

### 迁移指引

| 改动 | 影响 | 需要做什么 |
|---|---|---|
| 有效长度聚合（1） | `gene_counts_matrix.csv` 数值变化 | 重跑 Step 4 起（定量已缓存，几分钟） |
| LFC 收缩（6） | `log2FoldChange` 变化，显著基因数**可能下降** | 同上；下游 GSEA / 多组学评分需一并重跑 |
| 链特异性（5） | 若需从 `unstranded` 改为 `rf-stranded` | **必须删除 `03_kallisto_quant/` 重新定量**（数小时） |

```bash
# 重跑（定量结果自动跳过，直接从合并计数开始）
omicsweave transcriptomics run \
    --experiment 20250425_cimr_wt_mice_dio
```

> **待核实项**：三个实验的 `technical_params.strandedness` 目前填的是占位值
> `unstranded`（保持与既有结果一致）。目前主流 mRNA 建库为 dUTP 反向链特异，
> 真值很可能是 `rf-stranded`，请对照送检报告确认。

### 回归测试

```bash
python -m pytest tests/transcriptomics -v
```

共 47 个测试：28 个覆盖改进 1–6，14 个覆盖联合建模，5 个覆盖既有批次设计守卫。

---

## 0. 配置结构层级

配置加载器（`experiment_config.py`）支持三种配置布局，**向后兼容**：

### Layout C：多批次层级（推荐，2026-05+）

```
experiments/<experiment_id>/transcriptomics_data/
├── experiment.yaml          # 顶层入口：defaults、tracer_protocols、integration
└── <protocol_id>/
    ├── protocol.yaml         # 协议级共享配置（可选）
    └── <batch_id>/
        ├── config.yaml       # 运行入口（batch_file/groups_file/comparisons）
        ├── batch.yaml        # paths/reference/params/batch metadata（可选）
        └── groups.yaml       # conditions + sample_groups
```

- `rnaseq_pipeline.py --experiment <id>` 会按 `experiment.yaml` 的 `defaults.active_protocol/active_batch` 自动定位当前批次。
- 单批次输出默认规范到：`transcriptomics_results/runs/<protocol_id>/<batch_id>/...`。
- 新增多批次调度与整合：
  - `--all-batches`：运行 experiment 下全部已启用 protocol/batch
  - `--integrate`：运行多批次/多协议整合，输出到 `transcriptomics_results/integration/`

### Layout B：批次入口配置

```
experiments/<experiment_id>/transcriptomics_data/
├── config.yaml              # 运行入口（包含 batch_file/groups_file/comparisons）
├── batch.yaml               # paths/reference/params（可选，被 config.yaml 引用）
└── groups.yaml              # conditions + sample_groups
```

### Layout A：传统单文件（Legacy）

```
experiments/<experiment_id>/transcriptomics_data/
├── config.yaml              # 完整配置（paths/reference/params 全部在此）
└── groups.yaml              # conditions + mouse_groups
```

> 老实验可继续兼容运行。三种布局的优先级：C > B > A。

---

## 一、流程总览

本流程实现从原始 FASTQ 文件到差异表达基因（DEG）、通路富集分析（GO/KEGG/GSEA）、基因注释增强及自动化解读报告的端到端分析，全部在 Windows 下以 Python 脚本完成。

```
原始 FASTQ（网络盘）
  │
  ▼
Step 1 ── QC（Pure Python，采样前 100k reads）──→ 01_qc/
  │       ├─ <sample>_qc_report.json、qc_summary.csv
  │       └─ qc_gate_report.csv（阈值判定；on_qc_fail 决定处置）★新增
  │
  ▼
Step 2 ── Adapter Trimming（cutadapt）──→ 02_trimmed/（仅保留报告，中间文件即时删除）
  │
  ▼
Step 3 ── 伪比对 & 定量（kallisto quant，按配置传链特异性）──→ 03_kallisto_quant/
  │       ├─ abundance.tsv（转录本级别：est_counts + eff_length + tpm）
  │       └─ quant_summary.csv（比对率 p_pseudoaligned）★新增
  │
  ▼
Step 3a ── 转录本 → 基因聚合（lengthScaledTPM，纳入有效长度）★改进
  │        └─ gene_counts_matrix.csv / transcript_tpm_matrix.csv
  │        └─ run_provenance.json（参考指纹 + 工具版本 + 技术参数）★新增
  │
  ▼
Step 4 ── 差异表达分析 ──→ 04_de_results/
  │       ├─ PyDESeq2（首选，需 n≥2/组）+ apeGLM LFC 收缩 ★改进
  │       ├─ 多批次时用 ~batch + condition；混杂设计直接拒绝
  │       ├─ Wilcoxon 秩和检验（自动降级；先做文库归一化）★改进
  │       ├─ de_<tissue>_<comp>.csv（PyDESeq2 结果）
  │       └─ de_<tissue>_<comp>_wilcoxon.csv（Wilcoxon 降级结果）
  │
  ▼
Step 4a ── 基因注释增强（gene_annotate.py）──→ 04_de_results/（同目录，字段追加）
  │        ├─ Ensembl Gene ID + NCBI Gene ID 交叉引用
  │        └─ 为多组学对接提供标准化 ID
  │
  ▼
Step 5 ── 下游可视化 ──→ 05_result_plot/
  │       ├─ 1_volcano/（火山图，颜色可配置）
  │       ├─ 2_heatmap/（热力图，colormap 可配置）
  │       └─ 3_pca_clustering/（PCA 含批次着色 + PC~协变量关联检验）★改进
  │
  ▼
Step 6 ── 通路富集分析（Enrichr）──→ 05_result_plot/4_enrichment/
  │       ├─ GO/KEGG/Reactome/WikiPathways ORA 富集
  │       └─ 每组织 × 比较对 × 上调/下调 的 CSV + 条形图
  │
  ▼
Step 7 ── GSEA 预排序分析 ──→ 05_result_plot/5_gsea/
  │       ├─ 排序指标: log2FC × -log10(pvalue)
  │       └─ 每组织 × 比较对的 GSEA 结果 + NES 条形图
  │
  ▼
Step 8 ── 自动化解读报告 ──→ 06_interpretation_report/
  │       ├─ index.md（目录页）
  │       └─ <tissue>.md（每个器官 1 页证据表格 + 解读模板）
  │
  ▼
【多批次专用，≥2 批次时才执行】
  │
  ├─ Step 10 ── 跨批次联合建模（主路径）★新增 ──→ integration/joint_model/
  │             ├─ 溯源一致性断言 → merge_runs → 设计可识别性预检
  │             └─ ~batch + condition 联合 DE
  │
  └─ Step 9 ─── 多批次 Meta 分析（并行路径）──→ integration/
                ├─ within_protocol/<protocol>/<tissue>/  （随机效应 Meta + I²）
                ├─ cross_protocol/<tissue>/              （跨协议 consensus）
                └─ integration_report.md（并列呈现两条路径）
```

> **为什么两条路径都跑？** 联合模型效力更高（批次作为固定效应被吸收、
> 离散度共享估计），但看不到批次间异质性；Meta 分析给出 I²/τ²。
> 二者互补，报告中并列呈现并标注主路径。

---

## 二、技术栈与版本

| 工具/库 | 用途 | 版本 |
| :--- | :--- | :--- |
| **Python** | 脚本运行环境 | 3.12 |
| **kallisto** | 伪比对 & 转录本定量（Windows 原生二进制） | 0.51.1 |
| **cutadapt** | Adapter 修剪 | ≥4.0 |
| **PyDESeq2** | 差异表达分析（DESeq2 的 Python 实现） | ≥0.4 |
| **scipy** | Wilcoxon / Mann-Whitney U 检验（降级方案） | ≥1.11 |
| **statsmodels** | BH 多重检验校正（Wilcoxon padj） | ≥0.14 |
| **gseapy** | GSEA prerank & Enrichr 富集 | ≥1.0 |
| **pandas / numpy** | 数据处理 | ≥2.0 / ≥1.26 |
| **matplotlib / seaborn** | 绑图（颜色全部可配置） | ≥3.8 / ≥0.13 |
| **scikit-learn** | PCA、标准化 | ≥1.3 |
| **PyYAML** | 实验配置读取 | ≥6.0 |

### 参考基因组 / 转录组

- **物种**：小鼠（*Mus musculus*）
- **参考转录组**：Ensembl GRCm39 cDNA (`Mus_musculus.GRCm39.cdna.all.fa.gz`)，或 GENCODE vM35
- **基因注释来源**：Ensembl 转录本 ID → 基因 symbol（从 FASTA header 提取）
- **多组学 ID 映射**：`omicsweave/core/idmap/mouse_gene_id_mapping.py` 提供 gene symbol → Ensembl Gene ID / NCBI Gene ID

---

## 三、目录结构

```
omicsweave/transcriptomics/
├── README.md                           ← 本文件（详细说明书）
├── experiment_config.py                ← 实验配置加载器（YAML → Python dataclass）
├── rnaseq_pipeline.py                  ← 主流程脚本（Step 1–4 全自动）
├── rerun_de_analysis.py                ← 仅重跑 Step 4 差异表达（含 Wilcoxon 降级）
├── run_wilcoxon_de.py                  ← Wilcoxon DE 独立工具（也可被 pipeline 调用）
├── plot_results.py                     ← Step 5–7 统一 CLI 入口（薄包装）
├── run_enrichment.py                   ← Step 6 独立 CLI 入口（薄包装）
├── run_gsea.py                         ← Step 7 独立 CLI 入口（薄包装）
├── gene_annotate.py                    ← 基因注释增强（Ensembl/NCBI ID 交叉引用）
├── normalization.py                    ← ★基因层聚合（lengthScaledTPM）+ median-of-ratios 大小因子
├── qc_gate.py                          ← ★比对率采集 + QC 阈值门控
├── provenance.py                       ← ★参考/工具指纹采集（读写与跨批次比对见 core/provenance.py）
├── joint_model.py                      ← ★跨批次联合建模（merge_runs + ~batch + condition）
├── transcriptomics_integration.py      ← 多批次/多协议 DE 结果整合（Meta 分析 + consensus）
├── generate_interpretation_report.py   ← 自动化解读报告（Step 8）
├── finish_gsea.py                      ← 补跑缺失的 GSEA（断点续跑）
├── summarize_gsea.py                   ← 按代谢关键词过滤 GSEA 结果汇总
├── setup_kallisto.py                   ← 下载 kallisto Windows 二进制
├── download_transcriptome.py           ← 下载小鼠参考转录组 FASTA
│
├── plotting/                           ← Step 5: 纯可视化子包（火山图/热力图/PCA）
│   ├── __init__.py                     ← 导出 generate_volcano/heatmap/pca_plots
│   ├── _config.py                      ← 共享配置（apply_config, load_de_results, PLOT_COLORS 等）
│   ├── volcano.py                      ← 火山图（单图 + 汇总面板）
│   ├── heatmap.py                      ← DEG 热力图（按组别排序，clustermap）
│   └── pca.py                          ← PCA（全局/组织级别）+ 样本相关性热力图
│
├── enrichment/                         ← Step 6: 通路富集分析子包（Enrichr ORA）
│   ├── __init__.py                     ← 导出 run_enrichment
│   └── _run.py                         ← 自包含（配置加载 + API 调用 + 绘图）
│
├── gsea/                               ← Step 7: GSEA 预排序分析子包
│   ├── __init__.py                     ← 导出 run_gsea
│   └── _run.py                         ← 自包含（配置加载 + prerank + 绘图）
│
```

参考转录组与 kallisto 可执行文件不随包分发：`setup_kallisto` 与 `download_transcriptome` 会把它们下载到项目目录下的 `tools/` 与 `reference/`，实验配置中再填写对应路径。

---

## 四、实验配置（YAML 驱动）

### 4.1 配置目录约定

每个实验在仓库 `experiments/` 目录下拥有独立的配置子目录：

**Legacy 布局（Layout A）**：

```
experiments/<experiment_id>/transcriptomics_data/
├── config.yaml    ← 路径、参数、参考资源、绘图颜色
├── groups.yaml    ← 实验分组、动物编号映射、组别颜色
└── README.md      ← 本实验转录组配置说明
```

**多批次布局（Layout C）**：

```
experiments/<experiment_id>/transcriptomics_data/
├── experiment.yaml     ← 顶层：defaults、tracer_protocols、integration
└── <protocol_id>/
    └── <batch_id>/
        ├── config.yaml     ← 运行入口（batch_file/groups_file/comparisons）
        ├── batch.yaml      ← paths/reference/params（可选）
        └── groups.yaml     ← conditions + sample_groups
```

### 4.2 config.yaml 字段说明

#### 传统完整配置（Layout A）

```yaml
experiment_id: 20250425_cimr_wt_mice_dio   # 实验唯一标识
omics: transcriptomics

paths:
  raw_data_root: "\\\\server\\path\\01.raw_data"  # FASTQ 原始数据根目录（网络盘 UNC 路径）
  sample_list_tsv: "\\\\server\\path\\sample_list.tsv"  # 样本列表（TSV）
  output_dir: "analysis_results/${experiment_id}_output/transcriptomics_results"

reference:
  kallisto_exe: "omicsweave/transcriptomics/tools/kallisto/kallisto.exe"
  transcriptome_fasta: "omicsweave/transcriptomics/reference/Mus_musculus.GRCm39.cdna.all.fa.gz"
  kallisto_index: "${output_dir}/kallisto_index/mouse_cdna.idx"

params:
  threads: 8                # cutadapt / kallisto 并行线程
  qc_sample_size: 100000    # QC 采样读取数（null = 全部读取）
  tissues:                  # 要分析的组织列表
    - brain
    - liver
    - sq_adipose

  # 组织中文标签（可选，用于报告标题和图表轴标签）
  tissue_labels:
    brain: "大脑"
    liver: "肝脏"
    sq_adipose: "皮下脂肪"

  # 绘图颜色（可选，未配置的键自动使用内置默认值）
  plot_colors:
    volcano_up: "#D32F2F"        # 火山图上调基因
    volcano_down: "#1565C0"      # 火山图下调基因
    volcano_ns: "#BDBDBD"        # 火山图不显著基因
    heatmap_cmap: "RdBu_r"       # 热力图 colormap
    heatmap_vmin: -2             # 热力图 Z-score 下限
    heatmap_vmax: 2              # 热力图 Z-score 上限
    gsea_positive: "#D32F2F"     # GSEA 正向 NES 颜色
    gsea_negative: "#1565C0"     # GSEA 负向 NES 颜色
    enrichment_cmap: "RdYlBu_r"  # 富集分析条形图 colormap
    pca_tissue_cmap: "Set2"      # PCA 组织分色 colormap

  # 显式比较对（可选，不配置则自动生成所有两两配对）
  comparisons:
    - treat: DIO
      control: WT
    - treat: HSD
      control: WT

groups_file: "groups.yaml"   # 相对于 config.yaml 所在目录
```

#### 批次运行入口配置（Layout B/C）

```yaml
# config.yaml（运行入口）
experiment_id: 20250425_cimr_wt_mice_dio
batch_file: "batch.yaml"        # 可选，引用同目录下的 batch.yaml
groups_file: "groups.yaml"

# 显式指定比较对（支持多种格式）
comparisons:
  # 格式 1：dict（推荐）
  - treat: DIO
    control: WT
  # 格式 2：等效别名
  - numerator: HSD
    denominator: WT
  # 格式 3：list/tuple
  - ["DIO", "HSD"]
```

```yaml
# batch.yaml（paths/reference/params，可选）
paths:
  raw_data_root: "\\\\server\\path\\01.raw_data"
  sample_list_tsv: "\\\\server\\path\\sample_list.tsv"
  output_dir: "analysis_results/${experiment_id}_output/transcriptomics_results"

reference:
  kallisto_exe: "omicsweave/transcriptomics/tools/kallisto/kallisto.exe"
  transcriptome_fasta: "omicsweave/transcriptomics/reference/Mus_musculus.GRCm39.cdna.all.fa.gz"
  kallisto_index: "${output_dir}/kallisto_index/mouse_cdna.idx"

params:
  threads: 8
  qc_sample_size: 100000
  tissues: [brain, liver, sq_adipose]
  tissue_labels: {brain: "大脑", liver: "肝脏", sq_adipose: "皮下脂肪"}

# 批次元数据（可选）
protocol_id: glucose_tracer
batch_id: batch_2025_04

# ── 2026-08 新增 ────────────────────────────────────────────────
# 技术参数：决定本批次的计数能否与其他批次直接比较。
# 整合阶段会连同参考文件指纹一起做一致性断言。
technical_params:
  library_kit: "Illumina Stranded mRNA Prep"
  strandedness: "unstranded"       # unstranded | fr-stranded | rf-stranded ← 会传给 kallisto
  read_length: 150
  layout: "paired"                 # paired | single
  adapter_r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  adapter_r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  rin_min: 7.5                     # 本批次最低 RIN，用于评估降解差异

# QC 门控
batch_qc:
  enabled: true
  on_qc_fail: "flag"               # flag（仅告警）| exclude（移出下游 DE）| fail（中止）
  thresholds: {}                   # 可选：覆盖 experiment.yaml 的阈值，只写要改的项
```

#### 顶层实验配置（Layout C 专用）

```yaml
# experiment.yaml
experiment_id: 20250425_cimr_wt_mice_dio

defaults:
  active_protocol: glucose_tracer     # 当前活跃协议
  active_batch: batch_2025_04         # 当前活跃批次
  output_root: "analysis_results/${experiment_id}_output/transcriptomics_results"

tracer_protocols:
  glucose_tracer:
    enabled: true
    batches:
      - batch_2025_04
      - batch_2025_05

integration:
  enabled: true
  within_protocol:
    min_batches: 2                    # meta-analysis 最少批次数
  weighting:
    batch: inverse_variance           # 权重方法：inverse_variance 或 equal
    protocol: equal

  # ── 2026-08 新增 ──────────────────────────────────────────────
  # 整合前校验各批次的参考版本 / kallisto 版本 / 链特异性 / 聚合方法
  reference_mismatch_policy: "fail"   # fail（默认）| warn | skip

  # 跨批次联合建模（≥2 批次时才执行）
  joint_model:
    enabled: true
    gene_policy: "intersection"       # intersection（默认）| union
    duplicate_policy: "fail"          # fail（默认）| suffix
    min_samples_per_group: 2
    shrink_lfc: true
    check_provenance: true

# ── 2026-08 起 qc / alignment 两节真正执行（此前声明但不生效）──
analysis_defaults:
  qc:
    min_q20_pct: 80
    min_q30_pct: 70
    max_adapter_pct: 5
  alignment:
    min_pseudoaligned_pct: 60
  quantification:
    gene_aggregation: "lengthScaledTPM"   # lengthScaledTPM（默认）| simple_sum
```

#### 2026-08 新增配置字段速查

| 字段 | 所在文件 | 默认值 | 作用 |
|---|---|---|---|
| `analysis_defaults.qc.*` | experiment.yaml | 见上 | 碱基质量门控阈值（**现已生效**） |
| `analysis_defaults.alignment.min_pseudoaligned_pct` | experiment.yaml | 60 | 比对率门控（**现已生效**） |
| `analysis_defaults.quantification.gene_aggregation` | experiment.yaml | `lengthScaledTPM` | 基因层聚合方式 |
| `integration.reference_mismatch_policy` | experiment.yaml | `fail` | 跨批次参考不一致时的处置 |
| `integration.joint_model.*` | experiment.yaml | 见上 | 跨批次联合建模 |
| `technical_params.strandedness` | batch.yaml | `unstranded` | 传给 kallisto 的链特异性 |
| `technical_params.*`（其余） | batch.yaml | — | 一致性断言依据 + 溯源记录 |
| `batch_qc.on_qc_fail` | batch.yaml | `flag` | 不合格样本的处置方式 |
| `batch_qc.thresholds` | batch.yaml | `{}` | 批次级阈值覆盖 |

#### `plot_colors` 内置默认值

| 键 | 默认值 | 用途 |
|---|---|---|
| `volcano_up` | `#D32F2F` | 火山图上调（红色） |
| `volcano_down` | `#1565C0` | 火山图下调（蓝色） |
| `volcano_ns` | `#BDBDBD` | 火山图不显著（灰色） |
| `heatmap_cmap` | `RdBu_r` | 热力图 colormap |
| `heatmap_vmin` | `-2` | 热力图 Z-score 最小值 |
| `heatmap_vmax` | `2` | 热力图 Z-score 最大值 |
| `gsea_positive` | `#D32F2F` | GSEA 正 NES（红色） |
| `gsea_negative` | `#1565C0` | GSEA 负 NES（蓝色） |
| `enrichment_cmap` | `RdYlBu_r` | 富集分析 colormap |
| `pca_tissue_cmap` | `Set2` | PCA 组织分色 colormap |

> 用户只需在 `plot_colors` 中覆盖想要修改的键，其余自动回退到上表默认值。

### 4.3 groups.yaml 字段说明

```yaml
conditions:
  WT:
    display: "LFD-14W"           # 展示/汇报用名称（用于图表标签）
    diet: "D12450J"              # 饲料编号
    diet_desc: "Low-fat diet control"
    color: "#4C72B0"             # 该组别的绘图颜色（可选）
  DIO:
    display: "HFD-14W"
    diet: "D12492"
    diet_desc: "High-fat diet"
    color: "#DD8452"

# sample_groups 和 mouse_groups 两种键名均可使用（二选一）
sample_groups:          # 推荐键名
  "52": DIO
  "53": DIO
  "63": WT

# mouse_groups:         # 旧键名（兼容）
#   "52": DIO
```

- **`display`**：用于图表标签、比较名称和解读报告。不配置时回退到组别编码本身。
- **`color`**：用于 PCA、热力图列注释等所有按组别分色的图表。不配置时使用内置色板（`_DEFAULT_PALETTE`，色盲友好，最多 8 组）。
- **比较对生成**：若 `config.yaml` 未指定 `comparisons` 字段，则从 `conditions` 中所有组别自动两两配对，格式为 `{a}_vs_{b}`（按字母序排列）。若指定了 `comparisons`，则只运行指定的比较对。

### 4.4 变量替换

config.yaml / batch.yaml 支持以下变量：
- `${experiment_id}`：自动替换为 `experiment_id` 字段的值
- `${output_dir}`：自动替换为 `paths.output_dir` 解析后的绝对路径
- `${output_root}`：自动替换为顶层输出根目录
- `${protocol_id}` / `${batch_id}`：在多批次布局中替换为当前协议/批次 ID

---

## 五、配置加载器（experiment_config.py）

### 5.1 TranscriptomicsConfig 数据类

```python
@dataclass(frozen=True)
class TranscriptomicsConfig:
    # ---- 标识与路径 ----
    experiment_id: str
    config_path: Path
    output_root: Path                      # 输出根目录（Layout C: .../transcriptomics_results）
    output_dir: Path                       # 实际输出目录（Layout C: .../runs/<prot>/<batch>）
    protocol_id: str | None                # 当前协议 ID（Layout C 专用）
    batch_id: str | None                   # 当前批次 ID（Layout C 专用）
    experiment_path: Path | None           # experiment.yaml 路径（Layout C）
    protocol_path: Path | None             # protocol.yaml 路径（Layout C）
    batch_path: Path | None                # batch.yaml 路径（Layout C）

    # ---- 数据与参考 ----
    raw_data_root: str
    sample_list_tsv: str
    kallisto_exe: Path
    transcriptome_fasta: Path
    kallisto_index: Path

    # ---- 参数 ----
    threads: int
    qc_sample_size: int | None
    tissues: list[str]

    # ---- 分组 ----
    groups_path: Path
    conditions: dict[str, dict[str, Any]]
    mouse_groups: dict[str, str]           # 也接受 sample_groups 键名

    # ---- 比较对（可选） ----
    comparisons: list[tuple[str, str, str]] | None = None
                                         # 显式指定 (numerator, denominator, comp_name)
                                         # None = 自动生成所有两两配对

    # ---- 可选配置 ----
    tissue_labels: dict[str, str] | None = None      # 组织中文标签
    plot_colors: dict[str, Any] | None = None          # 绘图颜色覆盖
    integration: dict[str, Any] | None = None          # 多批次整合配置
```

加载方式：`load_transcriptomics_config(experiment_id=..., config_path=..., outdir_override=...)`

**配置解析流程**：

```
输入: --experiment <id> 或 --config <path>
  │
  ├─ experiment.yaml 存在? ──→ 解析 defaults.active_protocol/active_batch
  │   └─ 定位 <protocol>/<batch>/config.yaml
  │
  ├─ config.yaml 含 paths/reference/params? ──→ Layout A（Legacy 完整配置）
  │
  ├─ config.yaml 含 batch_file/groups_file? ──→ Layout B（批次入口）
  │   └─ 加载 batch.yaml 获取 paths/reference/params
  │
  └─ 兜底: 直接使用 config.yaml 全部字段
```

### 5.2 辅助函数

| 函数 | 签名 | 用途 |
|---|---|---|
| `generate_comparisons()` | `(conditions) → list[tuple[str, str, str]]` | 自动生成所有两两比较 `(numerator, denominator, comp_name)` |
| `resolve_comparisons()` | `(cfg) → list[tuple[str, str, str]]` | 优先使用 `cfg.comparisons`，否则自动两两配对 |
| `get_group_colors()` | `(conditions) → dict[str, str]` | 从 `conditions[].color` 提取组别颜色，未配置的使用默认色板 |
| `get_tissue_labels()` | `(cfg) → dict[str, str]` | 从 `params.tissue_labels` 提取，未配置的自动 capitalize |
| `get_comp_labels()` | `(conditions, comparisons) → dict[str, str]` | 生成 `{comp_name: "DisplayA vs DisplayB"}` 映射 |
| `get_cond_map()` | `(conditions) → dict[str, str]` | 生成 `{condition_code: display_name}` 映射 |
| `get_plot_colors()` | `(cfg) → dict[str, Any]` | 合并用户 `plot_colors` 与内置默认值 |
| `discover_experiment_protocols()` | `(experiment_id) → list[dict]` | 发现所有 protocol/batch 配置（多批次调度用） |
| `load_integration_config()` | `(experiment_id) → dict` | 加载 `experiment.yaml` 中的整合设置 |
| `resolve_output_root()` | `(experiment_id) → Path` | 解析输出根目录（兼容所有布局） |
| `config_arg_help()` | `() → str` | 返回 `--experiment / --config` 用法帮助文本 |

---

## 六、各步骤详解

### Step 1：Quality Control（纯 Python）

- **脚本**：`rnaseq_pipeline.py` → `fastq_qc()`
- **输入**：原始 `.fastq.gz` 文件（从网络盘读取）
- **输出**：`01_qc/<sample>_qc.json`
- **方法**：采样前 N 条 Reads（默认 100,000），统计：
  - 总 Reads 数、总碱基数
  - GC 含量、N 含量
  - Phred+33 质量评分分布（Q20%、Q30%、各位置均值）
  - 读长分布
- **说明**：等效于 FastQC，但无需安装 Java，在 Windows 上直接运行。

#### QC 阈值门控（2026-08 起真正执行）

此前 `analysis_defaults.qc` / `.alignment` / `batch_qc` 里的阈值没有任何代码读取，
QC 纯粹是描述性的。现由 `qc_gate.py` 执行：

```
experiment.yaml: analysis_defaults.qc        → min_q20_pct / min_q30_pct / max_adapter_pct
experiment.yaml: analysis_defaults.alignment → min_pseudoaligned_pct
batch.yaml:      batch_qc.thresholds         → 批次级覆盖（只写要改的项）
batch.yaml:      batch_qc.on_qc_fail         → 不合格样本的处置方式
```

| `on_qc_fail` | 行为 |
|---|---|
| `flag`（默认） | 仅告警，样本仍参与分析 |
| `exclude` | 不合格样本移出下游 DE（同时从计数矩阵中剔除） |
| `fail` | 直接抛异常中止流程 |

- **判定结果**：`01_qc/qc_gate_report.csv`（含每个失败项的具体数值，如 `q30_pct=55.00<70`）
- **R1/R2 取最差一侧**作为样本级指标——否则 R2 质量崩溃会被 R1 掩盖
- 阈值未配置的项**不参与判定**（避免省略某项时误判全部样本不合格）
- `batch_qc.enabled: false` 可整体关闭门控

### Step 2：Adapter Trimming（cutadapt）

- **脚本**：`rnaseq_pipeline.py` → `run_cutadapt()`
- **输入**：原始 `.fastq.gz`
- **输出**：`02_trimmed/<sample>_trim_report.txt`
- **方法**：
  - 使用 cutadapt 去除 Illumina TruSeq 接头序列
  - 质量修剪（Q20）
  - 最短读长过滤（20 bp）
- **磁盘优化**：修剪后的 `.fq.gz` 在完成 kallisto 定量后**立即删除**，仅保留修剪报告，避免占用本地磁盘。

### Step 3：伪比对与定量（kallisto）

- **脚本**：`rnaseq_pipeline.py` → `run_kallisto_quant()`
- **输入**：修剪后的 `.fq.gz`（临时文件）
- **输出**：`03_kallisto_quant/<sample>/abundance.tsv`
- **方法**：
  - 首次运行自动构建 kallisto index（`kallisto index`）
  - 对每个样本运行 `kallisto quant`，输出：
    - `target_id`：Ensembl 转录本 ID
    - `length`：有效转录本长度
    - `eff_length`：有效长度
    - `est_counts`：估计 Count 数
    - `tpm`：Transcripts Per Million
- **链特异性**（2026-08 起可配置）：
  - 由 `batch.yaml` 的 `technical_params.strandedness` 决定，翻译成 kallisto 开关：

    | 配置值 | kallisto 参数 | 适用 |
    |---|---|---|
    | `unstranded`（默认） | 无 | 非链特异性文库 |
    | `rf-stranded` | `--rf-stranded` | dUTP 反向链特异（目前主流 mRNA 建库） |
    | `fr-stranded` | `--fr-stranded` | 正向链特异 |

  - ⚠ 若文库实为链特异而按非链特异处理，**反义链与重叠基因的计数会被错误归属**。
    修改该值后必须删除 `03_kallisto_quant/` 重新定量。

- **基因聚合**（2026-08 改进，见下节）：
  - 从参考 FASTA header 提取 `transcript → gene symbol` 映射
  - 按 tximport `lengthScaledTPM` 聚合（纳入转录本有效长度）
  - 输出合并矩阵：`gene_counts_matrix.csv`、`transcript_tpm_matrix.csv`

- **运行溯源**：写出 `run_provenance.json`（参考文件 sha256、kallisto 版本、
  链特异性、聚合方法、技术参数），供跨批次合并前做一致性断言。

#### Step 3 计数类型说明：Raw Estimated Counts（非 CPM、非 TPM）

`04_de_results/gene_counts_matrix.csv` 中的数值为 **Raw Estimated Counts（原始估计计数）**，既不是 CPM 也不是 TPM。

**完整计算过程如下：**

```
1. kallisto quant（伪比对 + EM 算法）
   输入: 修剪后的 paired-end FASTQ
   输出: abundance.tsv（转录本级别）
     - est_counts: kallisto 通过 EM 算法估计的每个转录本的 read count
     - tpm:        Transcripts Per Million（标准化表达量）

2. 转录本 → 基因聚合（normalization.build_gene_level_counts，默认 lengthScaledTPM）
   ① 合并所有样本的 abundance.tsv（保留未取整的 est_counts，长度校正需要原值）：
      counts[sample]  = df['est_counts']
      tpm[sample]     = df['tpm']
      efflen[sample]  = df['eff_length']      ← 2026-08 起读取，此前从未使用

   ② 构建 tx2gene 映射（从 FASTA header 提取）：
      FASTA header 格式: >ENSMUST00000103301.3 cdna ... gene:ENSMUSG00000076500.3 gene_symbol:Actb
      → tx2gene = { "ENSMUST00000103301.3": { "gene_id": "ENSMUSG00000076500.3", "gene_symbol": "Actb" } }
      未精确命中的转录本会去掉版本号后缀再查（O(1) 索引，不再全表扫描）

   ③ 基因层 TPM 与有效长度：
      gene_tpm    = 该基因所有转录本 TPM 之和
      gene_efflen = 以转录本 TPM 为权重的有效长度加权平均
                    （丰度全为 0 的基因退回简单平均，避免 0/0）

   ④ lengthScaledTPM：
      counts = gene_tpm × gene_efflen 在样本间的均值
      再按列重标定，使每样本总 count 与原始求和一致（保留文库大小信息）

3. 输出矩阵
   - gene_counts_matrix.csv    ← 长度校正后的基因级计数（整数）
   - transcript_tpm_matrix.csv ← TPM（转录本级别，浮点数）
```

#### 为什么必须做长度校正

旧实现（`simple_sum`）直接把同一基因所有转录本的 `est_counts` 相加。
`est_counts` 正比于「摩尔丰度 × 有效长度」，因此当两组之间发生**异构体切换**
（长短转录本互换）时，基因总 count 会在**真实表达量完全不变**的情况下变化，
产生假阳性差异表达。这正是 tximport 存在的全部理由。

合成数据验证（同一基因的长短转录本有效长度 500 / 2000，两样本摩尔表达量相同）：

| 聚合方法 | 归一化后的假 log2FC |
|---|---|
| `simple_sum`（旧） | ≈ 2.0（即 4 倍假差异） |
| `lengthScaledTPM`（现默认） | < 0.1 |

**配置**（`experiment.yaml`）：

```yaml
analysis_defaults:
  quantification:
    gene_aggregation: "lengthScaledTPM"   # lengthScaledTPM（默认） | simple_sum（仅用于回溯比对）
```

> **注意聚合结果的读法**：`lengthScaledTPM` 会刻意把每列总量还原成原始文库大小，
> 因此聚合矩阵里仍残留一个对所有基因一致的缩放——这部分由 DESeq2 的 size factor
> 归一化消除。长度校正的效果要在归一化之后才显现，不要直接对聚合矩阵取比值判断。

> **已知局限**：聚合键仍是 `gene_symbol` 而非 `gene_id`，共享 symbol 的不同
> Ensembl 基因会被合并。改动会波及富集分析、绘图与全部多组学 symbol 映射，
> 属于独立议题，暂未变更。

**关键结论：**

| 输出文件 | 数值类型 | 单位 | 用途 |
|---|---|---|---|
| `gene_counts_matrix.csv` | **Raw Estimated Counts** | 整数（无标准化） | PyDESeq2 / Wilcoxon DE 输入 |
| `transcript_tpm_matrix.csv` | **TPM** | 浮点数（per-million 标准化） | 跨样本表达量比较 |
| `transcript_counts_matrix.csv` | **Raw Estimated Counts** | 整数（无标准化） | 转录本级别分析 |

> **为什么 PyDESeq2 需要 Raw Counts 而非 CPM/TPM？**
> DESeq2 内部会自行进行 library size normalization（通过 size factors），因此输入必须是未标准化的原始计数。预先做 CPM 或 TPM 标准化会破坏 DESeq2 的统计模型假设。
>
> **kallisto est_counts 与 featureCounts 的 counts 有何区别？**
> kallisto 使用 EM 算法对 multi-mapping reads 进行概率分配，输出的 est_counts 可能为小数（已取整为整数）。对于基因级别定量，二者高度一致（Pearson r > 0.95）。

#### Raw Counts 与 CPM/TPM 在差异表达分析中的角色

本流程使用 **Raw Counts** 作为差异表达分析的输入，这是 RNA-Seq DE 分析的业界标准做法。以下说明原因及 CPM/TPM 的正确使用场景。

**DESeq2 为什么必须用 Raw Counts：**

1. **DESeq2 内部自带标准化（median of ratios）**：DESeq2 通过计算 size factors 自动校正 library size 差异，预先做 CPM/TPM 标准化会导致 size factor 计算被扭曲，产生错误的统计推断。

2. **负二项分布要求整数计数数据**：DESeq2 假设基因表达服从负二项分布（NB），该分布天然适用于离散整数 counts。CPM/TPM 是连续浮点数，不符合 NB 分布的假设。

3. **均值-方差关系估计依赖 raw counts**：RNA-Seq 的核心特征是低表达基因方差大、高表达基因方差小（mean-dispersion relationship）。DESeq2 利用 raw counts 估计 dispersion 参数，标准化后的值会破坏这一关系。

**主流 DE 工具的输入要求对比：**

| 工具 | 要求输入 | 内部标准化方法 |
|---|---|---|
| DESeq2 / PyDESeq2 | Raw Counts | median of ratios（size factors） |
| edgeR | Raw Counts | TMM（trimmed mean of M values） |
| limma-voom | Raw Counts | voom 精度权重 |
| Wilcoxon / Mann-Whitney U | Raw Counts | 无（非参数检验，不依赖分布假设） |

**CPM/TPM 的正确使用场景（不用于 DE 分析）：**

| 场景 | 推荐使用 | 原因 |
|---|---|---|
| 差异表达（DE） | **Raw Counts** | DESeq2/edgeR 统计模型要求 |
| 跨样本/跨组织表达量比较 | **TPM** | 校正了基因长度和测序深度，可直接比较 |
| 热力图/表达量可视化 | **TPM** 或 log2(Raw Counts + 1) | 消除 library size 差异，图形可读性更好 |
| 低表达基因过滤 | **CPM** | `cpm > 1 in ≥n samples` 是常用过滤标准 |
| 特定基因绝对表达水平查看 | **TPM** | 反映每百万 reads 中的转录本数 |

**本流程的分工：**

```
Raw Counts（gene_counts_matrix.csv）
  ├─→ PyDESeq2 DE 分析（Step 4）
  ├─→ Wilcoxon DE 分析（Step 4 降级方案）
  └─→ PCA 聚类（Step 5，取 Top 高方差基因，Z-score 标准化）

TPM（transcript_tpm_matrix.csv）
  ├─→ 跨样本表达量比较（按需使用）
  └─→ 特定基因表达量查看（按需使用）
```

### Step 4：差异表达分析（PyDESeq2 + Wilcoxon 自动降级）

- **脚本**：`rnaseq_pipeline.py` Step 4 / `rerun_de_analysis.py` / `run_wilcoxon_de.py`
- **输入**：`gene_counts_matrix.csv` + `sample_metadata.csv`
- **输出**：`04_de_results/de_<tissue>_<comparison>.csv`（或 `_wilcoxon.csv`）

#### 双方法策略

流程采用**首选 PyDESeq2、自动降级 Wilcoxon** 的双方法策略：

1. **PyDESeq2**（首选）：当分子组和分母组各有 ≥2 个生物学重复时，运行标准 DESeq2 分析
2. **Wilcoxon 降级**：当任一组重复数不足（n<2），或 PyDESeq2 运行失败时，自动切换为 Wilcoxon 秩和检验（Mann-Whitney U）

```
对每个 tissue × comparison:
  ├─ n_num ≥ 2 AND n_den ≥ 2?
  │   ├─ Yes → 尝试 PyDESeq2
  │   │   ├─ 成功 → de_{tissue}_{comp}.csv
  │   │   └─ 异常 → 降级到 Wilcoxon → de_{tissue}_{comp}_wilcoxon.csv
  │   └─ No → 直接走 Wilcoxon → de_{tissue}_{comp}_wilcoxon.csv
  └─ 总样本数 < 2 → 跳过
```

#### PyDESeq2 方法（标准方案）

- 按组织拆分样本
- **设计公式自动选择**：元数据含 `batch` 列且该比较涉及 ≥2 个批次时用
  `~batch + condition`，否则用 `~condition`
- **设计矩阵秩校验**（`_validate_de_design()`）：混杂/秩亏设计直接抛异常拒绝拟合，
  且此时**禁用**无校正的 Wilcoxon 回退——宁可不出结果，也不产出无效结论
- 对每个两两比较运行 PyDESeq2：
  1. 构建 `DeseqDataSet` 对象
  2. 拟合负二项分布模型（`deseq2()`）
  3. 统计检验（Wald test）
  4. **apeGLM LFC 收缩**（p 值不变）
  5. 多重检验校正（Benjamini-Hochberg，输出 `padj`）
- 输出列：`gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, significant,
  log2FoldChange_unshrunk, lfc_shrinkage, lfc_shrink_coeff, design_formula,
  batch_adjusted, n_batches, batch_ids, source_mode`

##### LFC 收缩（2026-08 新增）

未收缩的 log2FC 对低表达基因会剧烈膨胀，而 `|log2FC|>1` 的筛选阈值恰好
**优先选中这些噪声基因**。apeGLM 用重尾 Cauchy 先验收缩效应量，不改变 p 值。
收缩前的值保留在 `log2FoldChange_unshrunk` 列。

> ⚠ **符号翻转陷阱（务必了解）**
>
> PyDESeq2 0.5.4 的 `lfc_shrink()` 只能收缩**设计矩阵里实际存在的系数**，
> 而参考水平由 formulaic 按**字母序**决定（`ref_level` 参数已废弃且完全无效）。
> 当 contrast 的分母恰好是可用系数代表的水平时，收缩结果的符号与 contrast **相反**。
>
> 实测：contrast=(hfd_14w, lfd_14w) 而可用系数为 `condition[T.lfd_14w]` 时，
> 收缩后符号一致性为 **0.00**（100% 翻转），而 p 值完好无损、图表照常出图——
> 极难被发现，却会让所有上/下调结论整体颠倒。
>
> `_resolve_shrinkage_coeff()` 解析系数名判断方向并按需取反，
> 并用符号一致性断言（强效应基因 ≥80% 一致）兜底；不一致时回退到未收缩值并报错。
> **升级 PyDESeq2 后请重跑 `python -m pytest tests/transcriptomics -k ShrinkageEndToEndTest`。**

- 可通过 `run_pydeseq2_analysis(..., shrink_lfc=False)` 关闭

#### Wilcoxon 方法（低重复数方案）

- **文库归一化**（2026-08 新增）：检验前先用 median-of-ratios 大小因子归一化。
  此前直接在**原始 count** 上做秩检验——Mann-Whitney 对测序深度差异毫无免疫力，
  深度相差 2 倍会在**所有基因**上造成约 +1 的 log2FC 系统偏移。
- **统计检验**：Mann-Whitney U（双侧，`scipy.stats.mannwhitneyu`）
- **log2FC 计算**：`log2((mean_num + 1) / (mean_den + 1))`（基于归一化后的计数）
- **多重检验**：Benjamini-Hochberg FDR 校正（`statsmodels.multipletests`）
- **n=1 情况**：仅输出 log2FC，p 值和 padj 为 NaN，所有基因标记为 not significant
- 输出列：`gene, baseMean, log2FoldChange, pvalue, padj, significant,
  min_achievable_pvalue, normalization, size_factor_range`

##### 样本量决定的检验分辨率上限（务必了解）

Mann-Whitney 精确检验共有 `C(n1+n2, n1)` 种秩分配，最极端的一种对应单侧
`1/C`，故**最小双侧 p = `2/C(n1+n2, n1)`**。这个上限与数据无关，纯由样本量决定：

| n1 vs n2 | 最小双侧 p | padj<0.05 是否可达 |
|---|---|---|
| 3 vs 3 | **0.100** | ❌ 不可能 |
| 3 vs 4 | 0.057 | ❌ 不可能 |
| 4 vs 4 | 0.029 | ✅ 可达（但 BH 校正后仍极难） |
| 5 vs 5 | 0.008 | ✅ |

**主实验 `20250425_cimr_wt_mice_dio` 正是 n=3/组**：一旦回退到 Wilcoxon，
结果必然是「0 个显著基因」——这是**方法学上限，不是生物学结论**。
流程会打印该上限并写入 `min_achievable_pvalue` 列。

#### 显著性判断标准

两个条件**同时满足**时标记为 `significant = True`：
- **统计学显著性**：`padj < 0.05`（Benjamini-Hochberg 校正后的 FDR）
- **生物学显著性**：`|log2FoldChange| > 1`（即至少 2 倍差异，2^1 = 2）

#### 跳过阈值

- 组织总样本数 < 2 → 跳过该组织
- 比较对总样本数 < 2 → 跳过该比较
- 某组重复数 < 1 → 跳过该比较

### Step 4a：基因注释增强（gene_annotate.py）

- **脚本**：`gene_annotate.py`（由 `rnaseq_pipeline.py` 和 `rerun_de_analysis.py` 自动调用）
- **输入**：DE 结果 CSV（`gene` 列为 gene symbol）
- **输出**：同目录 CSV 文件，追加以下字段：
  - `ensembl_gene_id`：Ensembl Gene ID（如 `ENSMUSG00000076500`）
  - `ncbi_gene_id`：NCBI Gene ID（如 `11461`）
  - `canonical_symbol`：标准化基因符号
  - `mapping_status`：映射状态（`exact` / `case_insensitive` / `synonym` / `not_found`）
- **依赖**：`omicsweave/core/idmap/mouse_gene_id_mapping.py`
- **用途**：为多组学整合（MFA-gene mapping）提供标准化 ID 交叉引用
- **也可独立使用**：

```python
from gene_annotate import annotate_gene_dataframe, annotate_csv_file

# 增强 DataFrame
enriched_df = annotate_gene_dataframe(de_results_df)

# 直接处理 CSV 文件
annotate_csv_file("path/to/de_results.csv")
```

### Step 5：下游可视化（plot_results.py）

所有图表的**绘图颜色均从配置文件读取**（`groups.yaml` 的 `color` 字段 + `config.yaml` 的 `plot_colors` 字段），未配置的使用内置默认值。

> **注意**：Enrichment 和 GSEA 已从 `plot_results.py` 中拆分为独立脚本（Step 6/7），默认仍会由 `plot_results.py` 调用，也可单独运行或通过 `--skip-enrichment` / `--skip-gsea` 跳过。

| 子步骤 | 输出目录 | 方法 |
| :--- | :--- | :--- |
| **火山图** | `1_volcano/` | 按 log2FC 与 -log10(padj) 绘制散点图，标注上调/下调基因数。生成单图 + 汇总面板（所有 tissue × comparison）。**支持 PyDESeq2 和 Wilcoxon 结果。** |
| **热力图** | `2_heatmap/` | 收集每组织所有比较的显著 DEG，log2 + Z-score 标准化，取 Top 80 高方差基因绘制 `sns.clustermap`，列注释按组别着色，基因名称显示在右侧 |
| **PCA** | `3_pca_clustering/` | 全局 PCA（log2-CPM，Top 5000 高方差基因，含**按批次着色**面板）+ 每组织 PCA + 样本 Pearson 相关性热力图 + **PC~协变量关联检验** |

#### 批次效应诊断（2026-08 新增）

合并多批次数据前必须先看到批次效应有多大。此前本模块完全没有批次相关逻辑。

**① 按批次着色的 PCA 面板**

`pca_global.png` 在配置了 `sample_batches` 时变为三联图：按分组 / 按组织 / **按批次**。

这一面板的价值在于：仅看「按分组着色」（改进前的唯一视图）往往看不出异常，
而「按批次着色」能直接暴露某个主成分是否由批次主导。

**② PC ~ 协变量关联检验**

对前 5 个主成分，逐一用 Kruskal-Wallis 检验其得分在 condition / tissue / batch
各水平间是否存在系统差异，输出 `pc_covariate_association.csv` 与热图
（值为 `-log10(p)`，>1.30 即 p<0.05）。

合成数据验证（注入 tissue ≫ batch ≫ condition 三级效应）：

| 协变量 | 命中主成分 | 方差占比 | p |
|---|---|---|---|
| tissue | PC1 | 76.6% | 7.8e-4 |
| batch | PC2 | 15.1% | 7.8e-4 |
| condition | PC3 | 3.1% | 7.8e-4 |

每个协变量都被精确定位到与其效应量相符的主成分。

**读法**：
- `batch` 在 PC1/PC2 上显著 → 批次效应主导全局变异，合并分析**必须**在模型中加入 batch 项
- `batch` 与 `condition` 同时在同一主成分上显著 → 两者可能混杂，需回到设计层面确认可分离性
- 检出批次显著时流程会打印警告，并提示 ComBat 校正矩阵**只可用于可视化**

> **批次校正的铁律**：DE 一律走模型里的 `batch` 项（`~batch + condition`）。
> ComBat/ComBat-seq 校正后的矩阵**不能**再喂给 DESeq2——那会重复校正、
> 低估方差、破坏 FDR 控制。校正矩阵仅用于出图。

### Step 6：通路富集分析 — Enrichr（run_enrichment.py）

- **脚本**：`run_enrichment.py`（独立可运行，也可被 `plot_results.py` 调用）
- **逻辑所在**：`enrichment/_run.py`（自包含子包）
- **输出**：`05_result_plot/4_enrichment/`
- **方法**：使用 gseapy 的 Enrichr 接口对上调/下调基因分别做 ORA（Over-Representation Analysis）
- **输出文件**：
  - `enrichr_<tissue>_<comp>_<up|down>.csv`：完整富集结果
  - `enrichr_<tissue>_<comp>_<up|down>.png`：Top 15 显著条目条形图

#### 使用的基因集库（小鼠特异）

| 基因集库 | 用途 |
|---|---|
| `GO_Biological_Process_2025` | GO 生物过程 |
| `GO_Molecular_Function_2025` | GO 分子功能 |
| `GO_Cellular_Component_2025` | GO 细胞组分 |
| `KEGG_2019_Mouse` | KEGG 通路 |
| `WikiPathways_2024_Mouse` | WikiPathways |
| `Reactome_Pathways_2024` | Reactome |

### Step 7：GSEA 预排序分析（run_gsea.py）

- **脚本**：`run_gsea.py`（独立可运行，也可被 `plot_results.py` 调用）
- **逻辑所在**：`gsea/_run.py`（自包含子包）
- **输出**：`05_result_plot/5_gsea/`
- **方法**：以 `log2FC × -log10(pvalue)` 为排序指标做 prerank GSEA
- **输出文件**：
  - `gsea_<tissue>_<comp>_results.csv`：GSEA 结果汇总
  - `gsea_<tissue>_<comp>.png`：正向/负向 NES 条形图
  - `gsea_<tissue>_<comp>/`：GSEA 中间文件

#### GSEA 排序指标

排序指标采用 **信号 × 置信度** 的乘积：

```
rank_metric = log2FoldChange × (-log10(pvalue))
```

- `log2FoldChange`：效应方向和大小（信号）
- `-log10(pvalue)`：统计置信度（越显著值越大）
- 两者相乘确保既有大效应、又有高置信度的基因排在最前/最后，适合 GSEA prerank 分析

#### 使用的基因集库（小鼠特异）

| 基因集库 | 用途 |
|---|---|
| `GO_Biological_Process_2025` | GO 生物过程 |
| `KEGG_2019_Mouse` | KEGG 通路 |
| `WikiPathways_2024_Mouse` | WikiPathways |
| `Reactome_Pathways_2024` | Reactome |

### Step 8：自动化解读报告

- **脚本**：`generate_interpretation_report.py`（Step 8）
- **输出**：`06_interpretation_report/`
- **内容**（每个器官 1 页 Markdown）：
  - 组别定义（从 `groups.yaml` 动态读取 display 名称）
  - 样本量统计
  - Top 上调/下调基因表格
  - 预设 marker 基因（代谢/炎症/纤维化相关）的 DE 统计量
  - GSEA Top NES 正向/负向条目
  - 代谢/免疫/纤维化关键词自动匹配
  - Enrichr 富集 Top 条目
  - "证据优先"解读模板

### Step 9：多批次/多协议整合（transcriptomics_integration.py）

- **脚本**：`transcriptomics_integration.py`（由 `rnaseq_pipeline.py --integrate` 调用）
- **输出**：`<output_root>/integration/`

#### 整合方法

| 方法 | 函数 | 说明 |
|---|---|---|
| **批次内 Meta-analysis** | `within_protocol_meta_analysis()` | 对同一 protocol 下多个批次的 DE 结果进行逆方差加权 meta-analysis，输出 `pooled_log2FoldChange`、`pooled_pvalue`、I² 异质性指标 |
| **跨协议 Consensus** | `cross_protocol_consensus()` | 跨 protocol 取效应方向一致性 consensus（符号一致比例 + Fisher 合并 p 值） |
| **整合报告** | `generate_integration_report()` | Markdown 格式的多批次整合报告 |

#### 整合输出结构

```
<output_root>/integration/
├── within_protocol/<protocol_id>/<tissue>/
│   └── de_<comparison>_integrated.csv
├── cross_protocol/<tissue>/
│   └── de_<comparison>_consensus.csv
└── integration_report.md
```

#### 整合前的一致性断言（2026-08 新增）

合并 ≥2 个批次前，先比对各批次 `run_provenance.json` 的四个关键字段：
`transcriptome_fingerprint`、`kallisto_version`、`strandedness`、`gene_aggregation`。

不一致时按 `integration.reference_mismatch_policy` 处置：

| 取值 | 行为 |
|---|---|
| `fail`（默认） | 拒绝整合该 protocol，不产出任何合并结果 |
| `warn` | 告警但继续 |
| `skip` | 不做校验 |

**缺少 `run_provenance.json` 的批次视为不一致**——无法证明可比，不能默认放行。

#### 配置

整合参数在 `experiment.yaml` 的 `integration` 字段中配置：

```yaml
integration:
  enabled: true
  within_protocol:
    min_batches: 2                  # meta-analysis 最少批次数
  weighting:
    batch: inverse_variance         # inverse_variance（推荐）或 equal
    protocol: equal
  reference_mismatch_policy: "fail" # fail | warn | skip
  joint_model:                      # 见 Step 10
    enabled: true
```

> **已知局限**：只被单个批次覆盖的基因走 `single_batch_passthrough` 分支，
> 拿到的是单批次 p 值，却与真正 Meta 出来的 p 值放在同一个 BH 家族里校正。
> 结果表的 `meta_note` / `meta_eligible` 列忠实记录了这一点，解读时须区分。
> 当前三个实验均为单批次，`--integrate` 的输出**全部**属于这一分支。

### Step 10：跨批次联合建模（joint_model.py）

- **脚本**：`joint_model.py`（由 `--integrate` 自动调用，或 `--joint-model-only` 单独运行）
- **输出**：`<output_root>/integration/joint_model/<protocol_id>/`
- **前提**：同一 protocol 下有 ≥2 个批次（否则明确跳过，不产生空目录）

#### 为什么需要它

Step 9 的 Meta 分析是「各批次独立跑完 DE → 对 log2FC 池化」。当**批次内分组平衡**时，
联合模型的统计效力显著更高：批次作为固定效应被吸收，离散度在全部样本上共享估计，
而不是分别拟合两次再池化。

合成数据对照（2 批次 ×(3 hfd + 3 lfd)，批次内平衡；真实效应 2.2× 上调，
另有 3× 批次效应作用于同一批基因）：

| 模型 | 真阳性检出 (padj<0.05) | 假阳性 | 中位 padj |
|---|---|---|---|
| 朴素 `~condition` | 31/60（52%） | 0/340 | 4.81e-02 |
| 联合 `~batch + condition` | **60/60（100%）** | 2/340 | 2.38e-27 |

> **准确理解这个收益**：批次内平衡时，朴素模型**不会**产生假阳性——它只是把批次
> 方差当噪声，损失的是**效力**而非正确性。假阳性风险来自**混杂**设计，
> 而那种情况联合模型会直接拒绝拟合。

#### 四段式流程

```
runs/<protocol>/<batch>/04_de_results/gene_counts_matrix.csv
  ├─ ① 溯源一致性断言 ── 参考/kallisto/链特异性/聚合方法不符即拒绝
  ├─ ② merge_runs ───── batch 列取自 run 目录；基因取交集；样本 ID 冲突检测
  ├─ ③ 设计可识别性预检 ─ 逐（组织×比较）做秩检验 + 残差自由度
  └─ ④ PyDESeq2 ~batch + condition ── 仅对预检通过的组合拟合
```

#### `joint_design_report.csv` 是关键安全产物

**在信任任何联合结果之前先读它。** 每行一个（组织 × 比较），记录：

| 列 | 含义 |
|---|---|
| `batch_condition_table` | batch × condition 交叉表（混杂与否一眼可见） |
| `identifiable` | 设计矩阵是否满秩 |
| `residual_df` | 残差自由度 = n −(批次数−1) −(条件数−1) −1 |
| `route` | `joint_model` / `meta_analysis_only` / `insufficient_data` |
| `reason` | 判定理由 |

只有 `route == joint_model` 的组合才会被拟合；其余明确退回 Meta 路径。

#### 三条安全护栏

| 护栏 | 触发条件 | 行为 |
|---|---|---|
| **混杂拒绝** | 批次 1 全 HFD、批次 2 全 LFD | 标记 `meta_analysis_only`，不产出联合结果 |
| **桥接样本拒绝** | 样本 ID 跨批次重复 | 默认 `fail` 并说明原因 |
| **参考不一致阻断** | 溯源字段不符 | `fail` 时连合并矩阵都不落盘 |

桥接样本（同一份 RNA 分装到多个批次）是**技术重复**，当作独立生物学重复会
低估方差、虚高显著性。需显式设 `duplicate_policy: suffix` 才放行。

#### 输出结构

```
integration/joint_model/<protocol_id>/
├── merged_gene_counts.csv        ← 合并计数矩阵
├── merged_sample_metadata.csv    ← 合并元数据（含 batch / source_batch / orig_sample）
├── merge_info.json               ← 合并统计（基因交并集、重命名记录、告警）
├── joint_design_report.csv       ← ★设计预检报告
└── <tissue>/
    ├── de_<comparison>_joint.csv
    └── de_for_multiomics_<comparison>.csv
```

#### 配置

```yaml
integration:
  joint_model:
    enabled: true
    gene_policy: "intersection"     # intersection（默认）| union（对缺失基因补 0，慎用）
    duplicate_policy: "fail"        # fail（默认）| suffix
    min_samples_per_group: 2
    shrink_lfc: true
    check_provenance: true
```

#### 两个无法靠模型解决的失败模式

1. **完全混杂**（批次 = 分组）：任何算法都无法分离，只能在**设计阶段**预防——
   每个批次必须包含所有实验组。
2. **批次效应波及过半转录组**：DESeq2 的 median-of-ratios 大小因子估计会到达
   崩溃点，注入**样本特异**的缩放误差，`~batch` 项吸收不掉。仿真扫描：

   | 批次影响基因占比 | 批次内 size factor CV | 纯无效应基因假阳性 |
   |---|---|---|
   | 10% | 0.015 | 3/340 |
   | 25% | 0.012 | 2/300 |
   | **50%** | **0.083** | **78/200** ← 崩溃 |
   | 75% | 0.052 | 18/100 |

   现实对应：整批 RNA 严重降解、或建库方式不同的批次。
   预警信号是批次在高方差主成分上极显著（见 Step 5 的诊断）——
   此时应核查批次质量，而不是指望统计模型补救。

---

## 七、运行方式

所有脚本统一支持 `--experiment / --config / --outdir` 三种参数传递方式。

### 7.1 完整流程（推荐）

```bash
# Step 1–4：QC → Trim → Quant → DE（含 Wilcoxon 自动降级 + 基因注释增强）
omicsweave transcriptomics run \
    --experiment 20250425_cimr_wt_mice_dio

# Step 5：下游可视化（火山图/热力图/PCA）
# 默认会调用 Step 6 Enrichment 和 Step 7 GSEA
omicsweave transcriptomics plot \
    --experiment 20250425_cimr_wt_mice_dio

# Step 8：解读报告
omicsweave transcriptomics report \
    --experiment 20250425_cimr_wt_mice_dio
```

### 7.2 分步运行（推荐用于大项目 / 调试）

Enrichment 和 GSEA 是耗时步骤（需联网调用 Enrichr API），可单独运行：

```bash
# Step 5：仅绘图（跳过 enrichment 和 GSEA）
omicsweave transcriptomics plot \
    --experiment 20250425_cimr_wt_mice_dio \
    --skip-enrichment --skip-gsea

# Step 6：仅运行 Enrichment（独立脚本）
python -m omicsweave.transcriptomics.run_enrichment \
    --experiment 20250425_cimr_wt_mice_dio

# Step 7：仅运行 GSEA（独立脚本）
python -m omicsweave.transcriptomics.run_gsea \
    --experiment 20250425_cimr_wt_mice_dio
```

> **提示**：Enrichment 依赖网络（调用 Enrichr REST API），GSEA 的 `permutation_num=1000` 需要一定计算时间。将它们拆分为独立步骤可以方便断点续跑和并行化。

### 7.3 仅重跑差异表达（Step 4）

不重新运行 QC / Trim / kallisto，直接从已有的 `gene_counts_matrix.csv` 重跑 DE（含基因注释增强）：

```bash
python -m omicsweave.transcriptomics.rerun_de_analysis \
    --experiment 20250425_cimr_wt_mice_dio
```

### 7.4 仅运行 Wilcoxon DE

强制使用 Wilcoxon 秩和检验（不尝试 PyDESeq2）：

```bash
python -m omicsweave.transcriptomics.run_wilcoxon_de \
    --experiment 20260312_bjsjth_wt_mice_dio_sleeve_gastrectomy
```

### 7.5 多批次运行（Layout C 专用）

```bash
# 运行当前活跃批次（由 experiment.yaml 的 defaults 指定）
omicsweave transcriptomics run \
    --experiment 20250425_cimr_wt_mice_dio

# 运行所有已启用批次
omicsweave transcriptomics run \
    --experiment 20250425_cimr_wt_mice_dio --all-batches

# 运行所有批次 + 完整整合（联合模型 + Meta 分析）
omicsweave transcriptomics run \
    --experiment 20250425_cimr_wt_mice_dio --all-batches --integrate

# 仅运行整合（复用已有的各批次产物）
omicsweave transcriptomics run \
    --experiment 20250425_cimr_wt_mice_dio --integrate

# 仅运行跨批次联合模型（不重新定量、不跑 Meta 分析）
omicsweave transcriptomics run \
    --experiment 20250425_cimr_wt_mice_dio --joint-model-only

# 联合模型也可独立调用
python -m omicsweave.transcriptomics.joint_model \
    --experiment 20250425_cimr_wt_mice_dio
```

> 单批次实验运行 `--joint-model-only` 会明确输出 `skipped_single_batch` 并退出，
> 不产生空目录或误导性结果。

### 7.6 辅助工具

```bash
# 补跑缺失的 GSEA（断点续跑，自动跳过已完成的组织/比较）
python -m omicsweave.transcriptomics.finish_gsea \
    --experiment 20250425_cimr_wt_mice_dio

# GSEA 代谢关键词汇总
python -m omicsweave.transcriptomics.summarize_gsea \
    --experiment 20250425_cimr_wt_mice_dio

# 下载参考转录组（多镜像源：GENCODE → Ensembl）
python -m omicsweave.transcriptomics.download_transcriptome

# 下载 kallisto Windows 二进制
python -m omicsweave.transcriptomics.setup_kallisto
```

### 7.7 临时覆盖输出目录

```bash
# 方式 1：命令行参数
python rnaseq_pipeline.py --experiment ... --outdir /path/to/custom_output

# 方式 2：环境变量
set TRANSCRIPTOMICS_OUTDIR=/path/to/custom_output
python rnaseq_pipeline.py --experiment ...
```

---

## 八、输出目录结构

### 8.1 Legacy 布局（Layout A/B）

```
analysis_results/<experiment>_output/transcriptomics_results/
├── sample_metadata.csv               ← 样本元数据（sample, mouse_id, tissue, condition[, batch]）
├── run_provenance.json               ← ★参考指纹 + kallisto 版本 + 链特异性 + 技术参数
├── kallisto_index/                    ← kallisto 索引文件
│   └── mouse_cdna.idx
├── 01_qc/                            ← 每个样本的 QC JSON 报告
│   ├── <sample>_qc_report.json
│   ├── qc_summary.csv                 ← 每样本 × R1/R2 的 QC 指标汇总
│   └── qc_gate_report.csv             ← ★阈值判定（qc_pass / qc_failures）
├── 02_trimmed/                        ← cutadapt 修剪报告（无 FASTQ，已清理）
│   └── <sample>_cutadapt.json
├── 03_kallisto_quant/                 ← 每个样本的 kallisto 定量结果
│   ├── quant_summary.csv              ← ★比对率汇总（p_pseudoaligned / n_processed）
│   └── <sample>/
│       ├── abundance.tsv
│       ├── abundance.h5
│       └── run_info.json
├── 04_de_results/                     ← 差异表达结果
│   ├── gene_counts_matrix.csv         ← 基因级别 Count 矩阵（lengthScaledTPM 校正后）
│   ├── transcript_tpm_matrix.csv      ← 转录本级别 TPM 矩阵
│   ├── transcript_counts_matrix.csv   ← 转录本级别 Count 矩阵
│   ├── resolved_comparisons.json      ← 实际执行的比较对（可追溯）
│   ├── de_<tissue>_<comp>.csv         ← PyDESeq2 DE 结果（含 Ensembl/NCBI ID + 收缩前 LFC）
│   ├── de_for_multiomics_<comp>.csv   ← 多组学稳定列名接口
│   └── de_<tissue>_<comp>_wilcoxon.csv  ← Wilcoxon 降级结果（如触发）
├── 05_result_plot/                    ← 可视化结果
│   ├── 1_volcano/                     ← 火山图（每 tissue × comparison + 汇总面板）
│   │   ├── volcano_<tissue>_<comp>.png
│   │   └── volcano_summary_panel.png
│   ├── 2_heatmap/                     ← DEG 热力图（每 tissue）
│   │   └── heatmap_<tissue>.png
│   ├── 3_pca_clustering/              ← PCA + 样本相关性 + 批次诊断
│   │   ├── pca_global.png             ← 按分组/组织/★批次 三联着色
│   │   ├── pca_<tissue>.png
│   │   ├── sample_correlation_heatmap.png
│   │   ├── pc_covariate_association.csv  ← ★PC ~ condition/tissue/batch 关联检验
│   │   └── pc_covariate_association.png
│   ├── 4_enrichment/                  ← Enrichr 富集
│   │   ├── enrichr_<tissue>_<comp>_<up|down>.csv
│   │   └── enrichr_<tissue>_<comp>_<up|down>.png
│   └── 5_gsea/                        ← GSEA prerank 结果
│       ├── gsea_<tissue>_<comp>_results.csv
│       ├── gsea_<tissue>_<comp>.png
│       └── gsea_<tissue>_<comp>/      ← GSEA 中间文件
└── 06_interpretation_report/          ← 自动化解读报告
    ├── index.md
    └── <tissue>.md
```

### 8.2 多批次布局（Layout C）

```
analysis_results/<experiment>_output/transcriptomics_results/
├── runs/                              ← 各批次独立输出
│   └── <protocol_id>/
│       └── <batch_id>/
│           ├── 01_qc/
│           ├── 02_trimmed/
│           ├── 03_kallisto_quant/
│           ├── 04_de_results/
│           ├── 05_result_plot/
│           ├── 06_interpretation_report/
│           ├── sample_metadata.csv
│           └── run_provenance.json    ← ★整合时的一致性断言依据
└── integration/                       ← 多批次整合输出
    ├── joint_model/<protocol_id>/     ← ★联合模型（主路径）
    │   ├── merged_gene_counts.csv
    │   ├── merged_sample_metadata.csv
    │   ├── merge_info.json
    │   ├── joint_design_report.csv    ← ★设计预检报告（先读它）
    │   └── <tissue>/
    │       ├── de_<comparison>_joint.csv
    │       └── de_for_multiomics_<comparison>.csv
    ├── within_protocol/<protocol>/<tissue>/   ← Meta 分析（并行路径）
    │   └── de_<comparison>_integrated.csv
    ├── cross_protocol/<tissue>/
    │   └── de_<comparison>_consensus.csv
    └── integration_report.md          ← 并列呈现两条路径 + 一致性校验结论
```

#### 结果文件如何区分来源

三条路径的结果格式相近，靠列值区分，切勿混用：

| 来源 | 文件 | `source_mode` | 含义 |
|---|---|---|---|
| 单批次 DE | `04_de_results/de_*.csv` | `single_results_directory` | 单批次内的效应 |
| 联合模型 | `integration/joint_model/**/de_*_joint.csv` | `joint_count_model_batch_adjusted` | 批次校正后的联合效应（效力最高） |
| Meta 分析 | `integration/within_protocol/**/de_*_integrated.csv` | — | 池化效应，另有 `i2` / `meta_note` |

---

## 九、复用到新实验

### 第 1 步：创建实验配置

```bash
# 方式 A：Legacy 布局
cp -r experiments/20250425_cimr_wt_mice_dio/transcriptomics_data \
      experiments/<new_experiment>/transcriptomics_data

# 方式 B：多批次布局
mkdir -p experiments/<new_experiment>/transcriptomics_data/<protocol>/<batch>/
# 然后创建 experiment.yaml + config.yaml + groups.yaml
```

### 第 2 步：修改配置

编辑 `config.yaml`（或 `batch.yaml`）：
- 修改 `experiment_id`
- 修改 `paths.raw_data_root` 和 `paths.sample_list_tsv` 指向新数据
- 修改 `params.tissues` 为新实验采集的组织列表
- 可选：配置 `params.tissue_labels`（中文标签）和 `params.plot_colors`（绘图颜色）
- 可选：配置 `comparisons` 字段显式指定比较对（不配置则自动两两配对）
- `output_dir` 使用 `${experiment_id}` 变量则只需确认路径格式

编辑 `groups.yaml`：
- 修改 `conditions`（组别名称 / display / 饲料信息 / color）
- 修改 `sample_groups` 或 `mouse_groups`（动物编号 → 组别映射）

### 第 3 步：运行

```bash
omicsweave transcriptomics run \
    --experiment <new_experiment>
```

---

## 十、与传统流程的对比

| 步骤 | 传统流程 | 本流程 | 说明 |
| :--- | :--- | :--- | :--- |
| 质控 | FastQC + MultiQC | Pure Python QC | 等效，无需安装 Java |
| 修剪 | fastp / Trim Galore | cutadapt | 等效 |
| 比对 | HISAT2 → BAM | **kallisto（伪比对）** | 快 10–100 倍，不生成 BAM |
| 计数 | FeatureCounts | kallisto quant | 直接输出 Count + TPM |
| 差异分析 | R DESeq2 / edgeR | **PyDESeq2 + Wilcoxon 降级** | 低重复数自动降级，无需手动切换 |
| 基因注释 | 手动 | **自动 Ensembl/NCBI ID 映射** | 为多组学对接提供标准化 ID |
| 富集 | R clusterProfiler | **gseapy（Enrichr + GSEA）** | 调用 Enrichr REST API |
| 报告 | 手动整理 | **自动生成 Markdown** | 可直接复制到 PPT |
| 配置 | 硬编码/命令行 | **YAML 配置驱动** | 换实验只改 YAML，不改代码 |
| 基因聚合 | tximport（R） | **lengthScaledTPM（Python）** | 纳入有效长度，避免异构体切换假阳性 |
| 效应量估计 | DESeq2 lfcShrink | **apeGLM（含方向守卫）** | 收缩低表达基因的膨胀 log2FC |
| 批次校正 | 手动 ComBat | **`~batch + condition` 联合模型** | 混杂设计自动拒绝，不硬拟合 |
| 多批次 | 手动合并 | **联合模型 + meta-analysis 双路径** | 效力与异质性兼顾 |
| QC 门控 | 人工看报告 | **阈值自动执行** | `on_qc_fail` 决定 flag/exclude/fail |

### 目前未覆盖的高级分析（可按需扩展）

- **ComBat-seq 校正矩阵**：仅用于可视化（DE 一律走模型里的 `batch` 项，不可双重校正）
- **ssGSEA**：单样本通路评分（gseapy 已支持，可扩展）
- **细胞去卷积**：CIBERSORTx / xCell（需导出 TPM 上传网页端）
- **一致性聚类分型**：ConsensusClusterPlus（适用于大样本分子分型）
- **按 `gene_id` 聚合**：当前聚合键为 `gene_symbol`，共享 symbol 的基因会被合并
- **rRNA / globin / 线粒体占比**、**3' 偏好性**：kallisto 不产生基因组坐标，需换比对流程
- **样本身份核验**：尚未接入流程作为门控（查看 marker 的调试脚本 `check_markers.py` 已移到本地课题工作区）

---

## 十一、多批次实验检查清单

多批次的关键在于**测序之前**的设计，而非事后的统计补救。以下按时间顺序。

### A. 送样之前（决定成败）

| # | 措施 | 为什么 |
|---|---|---|
| A1 | **每个批次都包含所有实验组** | 批次与分组完全混杂时，没有任何算法能分离技术与生物学差异。这是唯一无替代方案的措施 |
| A2 | **设置桥接样本**（同一份 RNA 分装到各批次） | 桥接样本之间的差异按定义 100% 是技术差异，是唯一能直接测量批次偏移的手段 |
| A3 | **冻结并记录技术参数** | 建库试剂盒、链特异性、读长、RIN — 填入 `batch.yaml: technical_params` |
| A4 | **确认链特异性** | 按非链特异处理链特异文库会错误归属反义链与重叠基因的计数 |
| A5 | **随机化处理顺序** | 避免「按组顺序提取/建库」把操作漂移变成组间差异 |

> A2 的桥接样本目前**不能**直接进入联合模型（会被 `duplicate_policy: fail` 拦下）。
> 它们的正确用途是评估批次效应大小，而非充当生物学重复。

### B. 分析时（流程自动执行）

| # | 检查 | 产物 | 不通过时 |
|---|---|---|---|
| B1 | 参考/工具/链特异性一致性 | `run_provenance.json` 比对 | 按 `reference_mismatch_policy` 拒绝整合 |
| B2 | 批次效应大小 | `pc_covariate_association.csv` | 批次在高方差 PC 上显著 → 先查批次质量 |
| B3 | 设计可识别性 | `joint_design_report.csv` | 混杂组合退回 Meta 路径 |
| B4 | 样本 ID 冲突 | `merge_info.json` | 默认拒绝合并 |

### C. 声明批次

批次信息写在 `groups.yaml` 的 `sample_batches`（**单次运行内**的多批次）：

```yaml
sample_batches:
  "52": b1          # 可用小鼠编号或完整样本 ID
  "53": b1
  "63": b2
```

**跨运行**的多批次（Layout C，各批次独立跑）无需此字段——
`joint_model.py` 会自动把 run 目录的 `batch_id` 作为批次标识。

> 若某批次内部又有 `sample_batches`，两层批次会组合为 `<batch_id>:<内部批次>`。

---

## 十二、常见问题

**Q：kallisto 与 HISAT2/STAR 的结果差异大吗？**
A：对于已知基因的表达定量，二者高度一致（相关系数 >0.95）。kallisto 的优势是速度（分钟级 vs 小时级）和磁盘占用（无 BAM 文件）。

**Q：为什么不直接用 R？**
A：本项目采用 Python 统一管理多组学分析流水线（转录组、代谢组、蛋白质组、MFA），避免在 R 和 Python 之间频繁切换。PyDESeq2 是 DESeq2 作者团队的官方 Python 实现，算法一致。

**Q：磁盘空间不够怎么办？**
A：流程默认在 trim + quant 后立即删除中间 `.fq.gz` 文件。最终占用空间主要是 `03_kallisto_quant/`（每个样本约 5–20 MB）和 `04_de_results/`。

**Q：网络盘读取很慢？**
A：cutadapt 和 kallisto 直接从网络盘流式读取。若网络不稳定，可先将 FASTQ 复制到本地 SSD，在 `config.yaml` 中修改 `raw_data_root` 即可。

**Q：Wilcoxon 降级什么时候会触发？**
A：当某个比较对中分子组或分母组的生物学重复数 < 2 时，PyDESeq2 无法拟合负二项分布模型，流程自动切换为 Wilcoxon 秩和检验。这在低重复数实验（如 n=1-2/组的外科手术实验）中尤为常见。输出文件名带 `_wilcoxon` 后缀以示区分。

> **例外**：若该比较涉及多个批次（元数据含 ≥2 个 `batch` 值），
> 无校正的 Wilcoxon 回退会被**禁用**并打印 `[BLOCKED]`。
> 多批次比较必须走可识别的批次校正模型，否则宁可不出结果。

**Q：如何自定义绘图颜色？**
A：在实验的 `config.yaml` 中配置 `params.plot_colors` 字段（详见 4.2 节），或在 `groups.yaml` 的 `conditions` 中设置 `color` 字段。未配置的键自动使用内置默认值。无需修改任何代码。

**Q：如何添加新的实验分组？**
A：只需在 `groups.yaml` 的 `conditions` 中添加新组别，在 `sample_groups`（或 `mouse_groups`）中添加对应的动物编号映射。如果未显式指定 `comparisons`，所有两两比较会自动生成。

**Q：Enrichment 和 GSEA 为什么拆分为独立脚本？**
A：这两个步骤耗时较长（Enrichr 需联网调用 REST API，GSEA 需 1000 次 permutation），且经常因网络超时或服务端限制而中断。拆分为独立脚本后可以：(1) 单独重跑，不影响其他可视化；(2) 方便并行化；(3) 支持 `plot_results.py` 通过 `--skip-enrichment` / `--skip-gsea` 跳过，先完成快速的可视化步骤。

**Q：如何使用显式比较对？**
A：在 `config.yaml` 中添加 `comparisons` 字段，指定 `treat`（或 `numerator`/`case`/`group_a`）和 `control`（或 `denominator`/`reference`/`group_b`）。不配置此字段则自动生成所有两两比较。支持 dict 和 list 两种格式（详见 4.2 节）。

**Q：多批次整合如何工作？**
A：需要使用 Layout C（`experiment.yaml` + 多批次目录结构）。在 `experiment.yaml` 的 `integration` 字段中启用并配置参数后，运行 `--integrate` 标志。流程会**同时**跑两条路径：(1) 跨批次联合模型 `~batch + condition`（主路径，效力更高）；(2) 逆方差/随机效应 meta-analysis（含 I² 异质性评估）。跨 protocol 另取 consensus。整合结果输出到 `<output_root>/integration/`，报告中并列呈现。

---

### 2026-08 更新相关

**Q：为什么我的显著基因数变少了？**
A：大概率是 LFC 收缩（改进 6）生效了。未收缩的 log2FC 对低表达基因会剧烈膨胀，而 `|log2FC|>1` 的阈值恰好优先选中这些噪声基因；收缩后它们不再过线。**这是修复而非退化**。p 值完全不受收缩影响，可对比 `log2FoldChange_unshrunk` 列查看收缩前后的差异。

**Q：改进 1（有效长度校正）会让结果变多少？**
A：取决于异构体使用是否在组间变化。若各组异构体比例稳定，长度校正接近恒等变换；若存在异构体切换，旧结果中该基因的差异**可能完全是长度伪影**。无法预先判断哪些基因受影响，所以必须整体重跑。

**Q：一定要重跑吗？能不能只跑新分析？**
A：必须重跑。改进 1 改变了 `gene_counts_matrix.csv` 的数值，所有以它为输入的下游（DE、火山图、热图、PCA、富集、GSEA、多组学评分）都基于旧值。定量结果已缓存，重跑从 Step 4 开始，通常几分钟。

**Q：QC 门控把样本剔除了，如何恢复？**
A：把 `batch.yaml` 的 `batch_qc.on_qc_fail` 改回 `flag`（默认值即为 `flag`，只告警不剔除）。先读 `01_qc/qc_gate_report.csv` 的 `qc_failures` 列确认失败原因，再决定是放行还是真的应该排除。

**Q：`strandedness` 填错了会怎样？如何确认正确值？**
A：按非链特异处理链特异文库，会让反义链与重叠基因的 reads 被错误归属，影响这些基因的定量准确度。正确值应来自送检报告的建库试剂盒说明——目前主流 mRNA 建库（如 Illumina Stranded mRNA Prep、NEBNext Ultra II Directional）都是 dUTP 反向链特异，对应 `rf-stranded`。**改动后必须删除 `03_kallisto_quant/` 重新定量**（流程按输出存在与否跳过，不会自动重算）。

**Q：整合时报「参考/技术参数不一致，拒绝整合」怎么办？**
A：先看报错列出的具体字段。常见原因：(1) 两批次用了不同版本的转录组 FASTA — 应统一参考后重跑各批次；(2) kallisto 版本不同 — 通常影响很小，可临时设 `reference_mismatch_policy: warn`；(3) 某批次缺少 `run_provenance.json`（2026-08 之前跑的）— 重跑该批次即可生成。**不建议直接设 `skip`**，那等于放弃这道保护。

**Q：`joint_design_report.csv` 里全是 `meta_analysis_only`，为什么？**
A：说明批次与分组混杂，或残差自由度不足。看 `batch_condition_table` 列的交叉表：若某个批次只含一个组（如 `b1:hfd=4` / `b2:lfd=4`），就是完全混杂——这是**设计问题，不是分析问题**，无法靠算法解决，只能在下次实验设计时保证每批次包含所有组。

**Q：联合模型和 Meta 分析结果不一致，信哪个？**
A：批次内分组平衡时优先信联合模型（效力更高、直接建模）。但若 Meta 分析报告高 I²（批次间异质性大），说明两批次的效应本身就不一致，此时**两个结果都值得怀疑**，应先查明异质性来源。联合模型会把异质性混入残差而不报告出来，这正是保留 Meta 路径的原因。

**Q：桥接样本为什么不能直接用于联合模型？**
A：桥接样本是同一份 RNA 的技术重复，把它们当作独立生物学重复会低估组内方差、虚高显著性（伪重复）。正确用途是量化批次偏移的大小。如确需保留，设 `duplicate_policy: suffix` 并自行在设计中处理其非独立性。

**Q：升级 PyDESeq2 后需要注意什么？**
A：重跑 `python -m pytest tests/transcriptomics -k ShrinkageEndToEndTest`。LFC 收缩依赖 PyDESeq2 的系数命名约定（`condition[T.<level>]`）和参考水平选择规则；这些若发生变化，方向判定可能失效。代码里有符号一致性断言兜底（不一致时回退到未收缩值并报错），但测试能更早发现问题。

**Q：`gene_annotate.py` 的作用是什么？**
A：它为 DE 结果 CSV 自动添加 `ensembl_gene_id` 和 `ncbi_gene_id` 列，便于后续多组学整合时通过标准化 ID 对接 MFA-gene mapping 字典。该步骤由 pipeline 自动调用，无需手动运行。