# omicsweave 包内模块总说明

`omicsweave` 包提供统一的多组学分析工具链，覆盖以下 5 条主流程：

1. **MFA 代谢流数据分析**
2. **转录组学数据分析（Bulk RNA-Seq）**
3. **代谢组学数据分析**
4. **蛋白质组学数据分析**
5. **多组学集成打分评估（Flux-centric）**

---

## 1. 总体架构

工具采用“**已安装的工具包 + 实验项目目录中的配置文件**”模式：

- 工具代码：`omicsweave` 包（`src/omicsweave/<omics>/`），通过 `pip install` 安装
- 实验配置：`<项目目录>/experiments/<experiment_id>/<omics>_data/`（或对应目录）
- 结果输出：`<项目目录>/analysis_results/<experiment_id>/<omics>_results/`

| 流程 | 包内目录 | 命令行入口 |
|---|---|---|
| MFA | `omicsweave/mfa/` | `omicsweave mfa run` |
| 转录组 | `omicsweave/transcriptomics/` | `omicsweave transcriptomics run` |
| 代谢组 | `omicsweave/metabolomics/` | `omicsweave metabolomics run` |
| 蛋白组 | `omicsweave/proteomics/` | `omicsweave proteomics run` |
| 多组学集成 | `omicsweave/integration/` | `omicsweave integration run` |

---

## 2. 运行环境与通用约定

- Python ≥ 3.10；按需安装对应 extras，例如 `pip install -e ".[mfa]"`，详见仓库根目录 README
- 项目目录：`--project-dir` > 环境变量 `OMICSWEAVE_PROJECT_DIR` > 从当前目录向上查找含 `experiments/` 的目录
- 配置文件：默认优先按 `--experiment <id>` 自动定位；也支持 `--config <path>`

示例（通用）：

```bash
omicsweave --project-dir <项目目录> <模块> <动作> --experiment <id>
```

---

## 3. 各流程说明

## 3.1 MFA 代谢流数据分析

### 核心功能

- 加载 `flux_raw_data.xlsx` 与 `mid_raw_data.xlsx`
- 差异通量统计（含 FDR、效应量）
- PCA、火山图、热图、TCA 循环、通路通量、MID 可视化
- 自动生成按组织分解读报告

### 主要脚本

| 脚本 | 作用 |
|---|---|
| `run_mfa_analysis.py` | 主流程（分析 + 出图；支持 `--report-only` 仅报告） |
| `plot_mfa_results.py` | 基于已生成表格重新出图 |
| `generate_mfa_interpretation_report.py` | 生成按组织 Markdown 报告 |

### 运行示例

```bash
omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio
omicsweave mfa plot --experiment 20250425_cimr_wt_mice_dio
omicsweave mfa report --experiment 20250425_cimr_wt_mice_dio
```

---

## 3.2 转录组学数据分析（Bulk RNA-Seq）

### 核心功能

- FASTQ 质控、cutadapt 修剪、kallisto 定量
- 差异表达分析（PyDESeq2；低重复时可降级 Wilcoxon）
- 火山图/热图/PCA
- Enrichment（GO/KEGG/Reactome/WikiPathways）与 GSEA
- 自动化解读报告

### 主要脚本

| 脚本 | 作用 |
|---|---|
| `rnaseq_pipeline.py` | Step 1-4：QC→Trim→Quant→DE |
| `plot_results.py` | Step 5-7：可视化 + 富集 + GSEA |
| `generate_interpretation_report.py` | Step 8：解读报告 |

### 运行示例

```bash
omicsweave transcriptomics run --experiment 20250425_cimr_wt_mice_dio
omicsweave transcriptomics plot --experiment 20250425_cimr_wt_mice_dio
omicsweave transcriptomics report --experiment 20250425_cimr_wt_mice_dio
```

---

## 3.3 代谢组学数据分析

### 核心功能

- 配置驱动的按组织分析（支持多组织）
- 数据预处理、PCA/PLS-DA、差异分析（T-test + VIP + FDR）
- 热图、KEGG ORA 富集、GSEA
- 自动生成解读报告（可 `--skip-report`）

### 主要脚本

| 脚本 | 作用 |
|---|---|
| `run_metabolomics_analysis.py` | CLI 入口（推荐） |
| `metabolomics_pipeline.py` | 实际流程编排（Step 1-6） |
| `generate_interpretation_report.py` | 报告生成 |

### 运行示例

```bash
omicsweave metabolomics run --experiment 20250425_cimr_wt_mice_dio
```

可选参数：

- `--steps 1,2,3`：只跑指定步骤
- `--skip-report`：跳过报告
- `--outdir <dir>`：覆盖配置中的输出目录；旧版 `--input/--output` 已移除

---

## 3.4 蛋白质组学数据分析

### 核心功能

- 读取 DIA-NN 报告并进行矩阵构建与预处理
- 差异分析（默认 `limma` moderated t-test，可配置 `wilcoxon`）
- 火山图、热图、PCA
- ssGSEA/通路分析
- 自动化解读报告

### 主要脚本

| 脚本 | 作用 |
|---|---|
| `run_proteomics_pipeline.py` | 主流程 |
| `plot_results.py` | 基于已有结果重绘图 |
| `generate_interpretation_report.py` | 生成组织级报告与汇总 |

### 运行示例

```bash
omicsweave proteomics run --experiment 20250425_cimr_wt_mice_dio
omicsweave proteomics report --experiment 20250425_cimr_wt_mice_dio
```

常用参数：

- `--force`：强制覆盖重跑
- `--skip-gsea`：跳过 ssGSEA

---

## 3.5 多组学集成打分评估（Flux-centric）

### 核心功能

- 以 **MFA 差异反应** 为主线整合转录组/蛋白组/代谢组证据
- 反应级打分（0-4.0）+ 通路级聚合 + 可视化 + Markdown 报告
- 内置 MFA gate，减少伪显著反应进入集成评分

### 打分逻辑（当前实现）

1. **MFA gate**：优先使用 `FDR < 0.05 且 |Cohen's d| > 0.8`；缺失时回退到 `pvalue + |log2FC|` 规则。  
2. **分层评分**（通过 gate 后）：
   - MFA：0-1.0（显著性 + 效应量）
   - Transcriptomics：0-1.0
   - Proteomics：0-1.2
   - Metabolomics：0-0.8
3. **总分**：0-4.0，按 `confidence_thresholds` 分级（阈值由配置文件控制，主实验常用 High≥3.0 / Medium≥2.0 / Low≥1.0）。

### 主要脚本

| 脚本 | 作用 |
|---|---|
| `run_multiomics_analysis.py` | 多组学集成主入口 |
| `integration/flux_omics_scorer.py` | 反应级打分核心 |
| `integration/pathway_convergence.py` | 跨组学通路收敛分析 |
| `report.py` | 自动化输出 `multiomics_report.md` |

### 运行示例

```bash
omicsweave integration run --experiment 20250425_cimr_wt_mice_dio
```

---

## 4. 推荐执行顺序（项目级）

为保证多组学集成输入完整，建议顺序：

1. MFA  
2. 转录组  
3. 代谢组  
4. 蛋白组  
5. 多组学集成

---

## 5. 结果目录总览

按实验 ID 统一输出到 `analysis_results/<experiment_id>/`，典型结构：

```text
analysis_results/<experiment_id>/
  mfa_results/
  transcriptomics_results/
  metabolomics_results/
  proteinomics_results/
  multiomics_results/
```

其中 `multiomics_results/` 会汇总前 4 条流程的数据并输出：

- `tables/`：打分表与通路收敛表
- `plots/`：集成可视化图
- `multiomics_report.md`：综合报告
