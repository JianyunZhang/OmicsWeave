# OmicsWeave 架构规划与实施计划（v2）

> 文档类型：架构规划 + 实施计划
> 版本：v2.1（2026-10-29），在 [v1](omicsweave_arch_design.md) 基础上修订；v1 中“单仓库 + 单分发包 + 可选依赖 + 分阶段拆分”的仓库组织决策继续有效，本文不再重复论证
> 当前状态：阶段 0（工程化打包）已完成，见第 3 节；阶段 1 进行中，进度见第 11 节
> 配套文档：[OmicsWeave × TheraWeave 开发规划](omicsweave_theraweave_development_plan.md) 是两个系统的总纲，**结果数据格式与交换契约以总纲为准**。本文的阶段 0–4 即总纲中的 O0–O4
> v2.1 修订：与总纲 v1.1 同步——结果正式格式改为 Parquet；结果目录加入运行 ID；比较增加角色与逆转配对；包结构增加 `contract`、`evidence`、`reporting`；阶段 2 增加证据层与交换契约 v1；第 13 节改为决策记录；新增本地课题工作区、集群与联网的预留约束和阶段 5

---

## 1. 结论摘要

### 1.1 对新增需求的评估

目标是把一个课题专用的多组学分析代码，改造为通用、可共享的软件包。新增需求逐条评估如下：

| 需求 | 评估 | 采纳方式 |
|---|---|---|
| project/experiment 项目配置，统一管理同条件下的多组学数据 | **合理且必要**，是现有代码最大的缺口：目前 5 个组学各有一套配置加载器（约 3,450 行），分组写法互不兼容 | 采纳。统一为 `project.yaml` + `experiment.yaml` + 两张样本表（第 5 节） |
| 预留组别、处理方式、实验时间的标记 | **合理**。现有代码用 `conditions`、`sample_groups`、`mouse_groups` 三种写法表达分组，且没有时间维度 | 采纳，并扩展为通用“实验因素”：组别、处理、时间、组织、批次及自定义因素，另加显式的比较定义（contrast） |
| 单组学独立运行、多组学可选整合 | **合理**，与 v1 的可选依赖设计一致 | 采纳。整合分为 3 种模式，各有明确的启用前提（第 9 节） |
| 规定统一的导入格式 | **方向合理，但需要分两层**。只认原始下机文件会把软件绑定在特定厂商、平台和外部软件上 | 改为两层：**原始数据层**登记下机文件并交给外部工具处理；**分析数据层**规定标准中间格式，所有统计与整合只读取这一层（第 6、7 节） |
| 转录组统一输入 Illumina `fastq.gz` | 合理 | 作为原始层格式；同时接受 kallisto/salmon 定量结果或计数矩阵直接导入 |
| 蛋白组统一输入 Thermo `.raw` | **限制过窄**。会排除 Bruker timsTOF（`.d`）、SCIEX 等平台，以及开放格式 mzML；`.raw` 的搜库与定量只能依靠外部软件（DIA-NN 等），Python 无法自行完成 | 以 Thermo `.raw` 为首个支持的原始格式，同时接受 mzML 和 DIA-NN 等软件的输出报告 |
| 代谢组输入 `.wiff` 或 `.raw` | **需细化**。`.wiff` 需与同名 `.wiff.scan` 一起使用；`.raw` 可能是 Thermo 的单个文件，也可能是 Waters 的目录；靶向与非靶向的处理流程完全不同，峰提取和注释应交给 MS-DIAL、MZmine、Skyline 等成熟工具 | 原始层登记并调用外部工具；分析层以“峰面积矩阵 + 代谢物注释表”为标准输入 |
| MFA 输入 Excel 求解结果表 | 合理（求解在外部软件完成） | 采纳，并**规定 Excel 模板**（工作表和列名）；同时接受 CSV/TSV；为现有自研软件和 INCA 等编写转换器 |
| 统一输出格式 | 合理，原需求未给出细节 | 规定结果目录结构、统一列名的结果表、图件、报告和溯源文件（第 8 节） |

**总体判断**：规划方向正确，可以整合进现有方案。唯一需要调整的是“统一导入格式”：
- **统一的是分析数据层的数据契约**，不限定原始文件格式。
- 原始格式通过适配器逐步扩充。
- 原始数据处理排在最后实施（阶段 3），因为它最依赖平台和厂商软件、工作量最大，而且离开它，分析层照样可以独立交付和发表。

### 1.2 修订后的核心原则

1. **代码与数据分离**：工具包只通过项目目录访问数据（阶段 0 已实现）。
2. **契约优先**：先定义标准中间格式和输出格式，再迁移各组学代码。
3. **一个项目模型**：所有组学共用一套项目、实验、样本、因素、比较的定义。
4. **外部工具走适配器**：kallisto、DIA-NN、msconvert、MS-DIAL 等通过统一的适配器接口调用或导入，不重新实现。
5. **单组学可独立运行，整合按需启用**。
6. **全程可追溯**：每次运行都记录配置、输入文件校验值、软件版本和随机种子。
7. **只保留通用能力**：课题专用的脚本和图件移回实验项目仓库。

---

## 2. 定位与范围

**定位**：面向小鼠等模式生物的 bulk 多组学下游分析与整合工具包，特色是**以 MFA 通量为锚点、推断单位正确的多组学证据整合**。

**范围内**
- 项目/实验配置与校验
- 各组学标准化导入、质控、预处理、差异分析、富集分析、可视化与报告
- 多组学整合（3 种模式）
- 外部工具适配器：调用外部工具，或导入其输出

**范围外（不重新实现）**
- 测序比对与定量算法
- 质谱搜库与定量算法
- 非靶向代谢物峰提取与注释算法
- MFA 模型求解
- 单细胞与空间组学（留作以后扩展）

---

## 3. 现状（阶段 0 完成后）

| 项目 | 状态 |
|---|---|
| 打包 | 已有 `pyproject.toml`，可以 `pip install -e ".[all]"`；extras 按组学划分 |
| 导入方式 | 已统一为 `omicsweave.*`；`sys.path` 修改已全部删除 |
| 数据位置 | `common/project.py` 统一解析项目目录：`--project-dir`、环境变量 `OMICSWEAVE_PROJECT_DIR`、向上查找 `experiments/` |
| 命令行 | `omicsweave <模块> <动作>`，共 5 个模块、14 个动作 |
| 测试 | 185 个；Python 3.10 与 3.12 下均为 182 通过、3 跳过（需要真实数据） |
| CI | 已编写，尚未在 GitHub 上运行过 |
| 代码规模 | 143 个 Python 文件，约 4.4 万行（另有测试约 0.5 万行） |

阶段 0 之后仍然存在的结构问题：
- 5 套配置加载器、3 套溯源记录、5 套绘图配置、3–4 份重复的火山图、热图、PCA 实现。
- ID 映射是硬编码在代码里的小鼠专用字典。
- 课题专用脚本与通用代码混在一起。
- 各组学的输出格式互不统一，整合模块的读取器（`multiomics/loaders/`）需要分别解析每个组学的输出目录。

---

## 4. 目标架构

### 4.1 分层

```text
┌──────────────────────────────────────────────────────────────┐
│ cli：init / validate / run / integrate / status / migrate     │
├──────────────────────────────────────────────────────────────┤
│ integration：flux_anchored / pathway_convergence / sample_level│
├──────────────┬──────────────┬──────────────┬─────────────────┤
│ transcriptomics│ proteomics  │ metabolomics │ mfa             │  各组学：导入 → QC →
│              │              │              │                 │  预处理 → 差异 → 富集 → 图/报告
├──────────────┴──────────────┴──────────────┴─────────────────┤
│ adapters：kallisto / salmon / DIA-NN / msconvert / MS-DIAL …  │  原始层 → 分析层
├──────────────────────────────────────────────────────────────┤
│ visualization：主题、通用图（火山图、热图、PCA…）             │
├──────────────────────────────────────────────────────────────┤
│ core：项目模型、配置校验、数据契约读写、统计、ID 映射、溯源   │
└──────────────────────────────────────────────────────────────┘
```

**依赖规则**
- 上层可以依赖下层，下层不能依赖上层。
- 各组学模块之间不能相互导入。
- `integration` 只读取各组学按数据契约输出的结果，不导入各组学的内部实现。

**为集群部署和多客户端联网预留的约束**（当前只实现 Windows/Linux 单机运行，见第 13 节）：
- 运行函数接受显式的项目上下文（项目目录、配置对象），不依赖进程级全局状态。现有的 `set_project_dir()` 和各模块的全局变量（如 `rnaseq_pipeline` 中的 `KALLISTO_EXE`、`plot_results` 中的 `PLOT_DIR`）在阶段 1–2 逐步改为参数传递。
- 每次运行写入独立的 `<run_id>/` 目录，不同运行之间互不覆盖。
- 单个步骤可以单独调用，输入输出都是文件，便于以后交给作业调度器或工作流引擎拆分执行。
- 文件读写都经过 `core` 的路径接口，以后可以换成共享存储或对象存储。

### 4.2 目标包结构

```text
src/omicsweave/
├── core/
│   ├── project.py          # 项目目录定位 ✅
│   ├── model.py            # project.yaml / experiment.yaml 的 Pydantic 模型 ✅
│   ├── config.py           # YAML 与样本表的读取 ✅
│   ├── validate.py         # 结合样本表与文件系统的校验（omicsweave validate）✅
│   ├── scaffold.py         # 项目模板（omicsweave init）✅
│   ├── provenance.py       # 文件指纹、run_provenance.json 读写与跨批次比对 ✅
│   ├── statistics/         # nonparametric / meta_analysis / stability ✅
│   ├── idmap/              # ID 映射（暂为原有字典，阶段 2 改为数据文件 + 解析器）
│   └── safety.py           # pandas 数值安全检查 ✅
├── contract/               # 分析层数据集、结果表与结果包的 Pydantic 模型、读写和校验；导出 JSON Schema（草案 ✅）
├── adapters/
│   ├── base.py             # 适配器接口：检查工具、运行或导入、产出数据契约
│   ├── rnaseq/             # kallisto、salmon、计数矩阵导入
│   ├── proteomics/         # DIA-NN、MaxQuant/FragPipe 结果导入、raw→mzML
│   ├── metabolomics/       # msconvert、MS-DIAL/MZmine/Skyline 结果导入
│   └── mfa/                # 标准 Excel 模板、自研软件格式、INCA 等转换
├── transcriptomics/  proteomics/  metabolomics/  mfa/
│   └── 各自包含 qc.py / preprocess.py / differential.py / enrichment.py / plots.py / report.py（报告章节）
├── integration/            # 跨组学统计整合（已由 multiomics/ 改名 ✅；原 integration/ 子包改名为 methods/）
├── evidence/               # 实体标准化、同源映射、证据记录、候选通路/靶点假设
├── reporting/              # 由结果包生成 HTML/PDF 报告
├── visualization/          # theme.py：统一绘图主题 ✅；通用图形阶段 2 迁入
└── cli.py
```

交换契约的 JSON Schema 与标准样例结果包发布在仓库根目录的 `schemas/evidence-exchange/v<主版本>/`（见总纲 5.1 节）。`integration` 负责跨组学的统计分析，`evidence` 负责把统计结果转成标准化的证据与假设记录，两者不合并。

---

## 5. 项目 / 实验配置模型

### 5.1 层级

```text
Project（一个课题或队列，如 DIO 小鼠项目）
└── Experiment（同一设计、同一批动物下的一次多组学实验）
    ├── Subject（动物个体，跨组学关联的主键）
    ├── Sample（一次组学测量，属于某个 Subject，有组织、时间、批次）
    ├── Factor（组别、处理、时间等实验因素）
    ├── Contrast（要做的比较）
    └── Omics 设置（每个组学：是否启用、输入层级与格式、处理参数）
```

与现有配置的对应关系：

| 现有写法 | 新模型 |
|---|---|
| `conditions`、`sample_groups`（样本→条件）、`mouse_groups`（动物→条件） | `subjects.tsv` 的因素列 + `samples.tsv` 的 `subject_id` |
| `protocol / batch` 分层目录 | `samples.tsv` 的 `batch` 列 + 组学设置中的处理参数 |
| `estimand`（case/control、干预时间、示踪剂…） | `experiment.yaml` 的 `contrasts` 与 `metadata` |
| `multiomics/examples/sample_map.example.csv` | 已是 `samples.tsv` 的雏形（assay→omics，primary_sample_id→subject_id） |
| 各组学独立的 `<omics>_data/config.yaml` | 合并到 `experiment.yaml` 的 `omics:` 段 |

### 5.2 目录约定

```text
<项目目录>/
├── project.yaml
├── experiments/<experiment_id>/
│   ├── experiment.yaml
│   ├── subjects.tsv
│   ├── samples.tsv
│   └── raw/ 或外部路径        # 原始文件可以放在别处，samples.tsv 中写路径即可
├── references/               # 参考基因组、FASTA、通路库、MFA 网络模型
└── analysis_results/<experiment_id>/<run_id>/...   # 每次运行一个目录，结构见第 8.1 节
```

### 5.3 `project.yaml` 示例

```yaml
schema_version: 1
project:
  id: dio_mouse
  title: 饮食诱导肥胖小鼠多组学研究
  species: mus_musculus            # 决定默认的 ID 映射与通路库
  paths:                           # 均可省略，以下为默认值
    experiments: experiments/
    results: analysis_results/
    references: references/
    analysis: analysis/            # 课题专用脚本
  references:
    transcriptome: references/Mus_musculus.GRCm39.cdna.all.fa.gz
    annotation: ensembl_113
    proteome_fasta: references/UP000000589_10090.fasta
    mfa_network: references/central_carbon_v1.yaml
defaults:
  padj_method: fdr_bh
  alpha: 0.05
  random_seed: 20261029
experiments:
  - experiments/20250425_dio_hfd_lfd
  - experiments/20250613_dio_glp1
```

### 5.4 `experiment.yaml` 示例

```yaml
schema_version: 1
experiment:
  id: 20250613_dio_glp1
  description: HFD 小鼠 GLP-1 干预，多组织多组学
  metadata:
    tracer: U-13C glucose
    intervention_duration: 4 weeks

design:
  subjects: subjects.tsv
  samples: samples.tsv
  factors:                          # 声明因素：取值、参照水平、所在表格
    diet:      {level: subject, levels: [LFD, HFD], reference: LFD}
    treatment: {level: subject, levels: [Saline, GLP1], reference: Saline}
    timepoint: {level: sample, type: numeric, unit: week}
    tissue:    {level: sample}
  contrasts:
    - id: GLP1_vs_Saline_in_HFD
      factor: treatment
      numerator: GLP1
      denominator: Saline
      subset: {diet: HFD}
      role: treatment_effect        # disease_effect / treatment_effect / time_effect / other
    - id: HFD_vs_LFD_saline
      factor: diet
      numerator: HFD
      denominator: LFD
      subset: {treatment: Saline}
      role: disease_effect
  reversal_pairs:                   # 疾病效应被治疗逆转时标记 reversal（见总纲 2.3 节）
    - {disease: HFD_vs_LFD_saline, treatment: GLP1_vs_Saline_in_HFD}
  stratify_by: [tissue]             # 各组织分别分析
  covariates: [batch]               # 批次 > 1 个水平时自动加入设计公式

omics:
  transcriptomics:
    enabled: true
    input: {level: raw, format: fastq_gz, layout: paired, strandedness: reverse}
    adapter: {name: kallisto, executable: kallisto, index: references/kallisto.idx}
    analysis: {method: deseq2, min_count: 10}
  proteomics:
    enabled: true
    input: {level: raw, format: thermo_raw, acquisition: dia}
    adapter: {name: diann, executable: tools/diann/diann.exe, fasta: "${project.references.proteome_fasta}"}
    analysis: {method: limma, imputation: none}
  metabolomics:
    enabled: true
    input: {level: table, format: peak_table, source: msdial, mode: untargeted}
    analysis: {normalization: pqn, method: welch}
  mfa:
    enabled: true
    input: {level: table, format: omicsweave_mfa_xlsx, files: [mfa/flux_results.xlsx]}
    analysis: {inference_unit: subject, statistic: auto}

integration:
  enabled: true
  modes: [flux_anchored, pathway_convergence]
```

### 5.5 样本表

`subjects.tsv`：每只动物一行。

| 列 | 必填 | 说明 |
|---|---|---|
| `subject_id` | 是 | 动物/个体编号，项目内唯一，跨组学关联的主键 |
| 因素列（如 `diet`、`treatment`） | 按 `factors` 声明 | `level: subject` 的因素放在这里 |
| `sex`、`age`、`strain`、`cage` 等 | 否 | 可作为协变量 |

`samples.tsv`：每次组学测量一行（同一组织块分别做 RNA 和蛋白，记为两行）。

| 列 | 必填 | 说明 |
|---|---|---|
| `sample_id` | 是 | 项目内唯一 |
| `subject_id` | 生物学样本必填 | 对应 `subjects.tsv`；QC、空白、混样等非生物学样本可留空 |
| `omics` | 是 | `transcriptomics` / `proteomics` / `metabolomics` / `mfa` |
| `tissue` | 是 | 组织或体液 |
| `timepoint` | 否 | 数值，单位在 `factors` 中声明 |
| `batch` | 否 | 测序或上机批次 |
| `sample_type` | 否 | `sample` / `qc` / `blank` / `pool` / `standard`，默认 `sample` |
| `run_order` | 否 | 质谱进样顺序，用于 QC 漂移校正 |
| `specimen_id` | 否 | 同一组织块被多个组学测量时用于关联 |
| `files` | 原始层必填 | 相对项目目录的路径，多个文件以 `;` 分隔（如 R1;R2） |
| `md5` | 否 | 原始文件校验值，用于溯源 |

### 5.6 校验规则（`omicsweave validate`）

- 所有 ID 唯一；`samples.subject_id` 都能在 `subjects.tsv` 中找到。
- 因素取值都在 `levels` 范围内；每个比较的两组在每个分层内都至少有 2 个生物学重复，否则给出警告并说明可达到的最小 p 值。
- 启用的组学都有样本，文件存在且格式与 `input.format` 一致（例如 `.wiff` 必须有对应的 `.wiff.scan`）。
- 检查外部工具是否可执行，并记录其版本。
- 整合模式的前提条件已满足（见第 9 节）。

---

## 6. 数据输入规范（两层）

### 6.1 原始数据层（登记 + 外部工具处理）

| 组学 | 首批支持的格式 | 后续扩展 | 处理工具（外部） | 产出的分析层数据 |
|---|---|---|---|---|
| 转录组 | Illumina `fastq.gz`（单端/双端） | — | kallisto（已有）、salmon；修剪用 cutadapt | 转录本/基因计数与 TPM |
| 蛋白组 | Thermo `.raw`（DIA） | mzML、Bruker `.d`、DDA | DIA-NN（已有）；`.raw`→mzML 用 ThermoRawFileParser | 前体/蛋白强度矩阵 |
| 代谢组（靶向） | SCIEX `.wiff` + `.wiff.scan` | Thermo `.raw`、mzML | Skyline / SCIEX OS / MultiQuant | 峰面积矩阵 |
| 代谢组（非靶向） | Thermo `.raw`、mzML | Waters `.raw`（目录） | msconvert → MS-DIAL / MZmine | 特征强度矩阵 + 注释表 |
| 同位素示踪代谢组 | 同上 | — | 上述工具 + 天然丰度校正（IsoCor 等） | MID 表 |
| MFA | 求解结果 Excel（标准模板，或旧版自研软件导出） | Parquet 目录（全部拟合结果）、INCA 导出 | 外部求解软件（自研软件、INCA 等） | 通量表 + MID 拟合表 |

注意事项：
- **原始文件只登记、不复制**：大文件留在原位置，`samples.tsv` 中写路径，可选填 MD5。
- **厂商格式受平台限制**：读取 SCIEX、Waters 等厂商格式通常需要 Windows 上的厂商库（ProteoWizard msconvert 或其 Wine 容器）；Thermo `.raw` 可通过 ThermoRawFileParser 跨平台读取。
- **外部工具的许可需逐一核实**：DIA-NN、MSFragger 等是闭源软件或仅限学术使用，工具包只负责调用，不随包分发。
- **两条路径都要支持**：如果用户已经用别的软件处理完原始数据，可以直接从分析层导入，跳过原始层。

### 6.2 分析数据层（所有分析的唯一入口）

每个组学经过适配器后，都产出同一种结构的“数据集”：

```text
<omics>/data/
├── matrix.parquet        # 行=特征，列=sample_id
├── features.parquet      # 特征注释：feature_id、名称、各数据库 ID、类型
└── dataset.json          # 数值语义与来源说明
```

`dataset.json` 必填字段：

| 字段 | 示例 |
|---|---|
| `omics`、`tissue` | `proteomics`、`liver` |
| `value_type` | `counts` / `tpm` / `intensity` / `peak_area` / `mid_fraction` / `flux` |
| `scale` | `linear` / `log2` |
| `normalization` | `none` / `deseq2_size_factor` / `pqn` / `median` … |
| `feature_id_type` | `ensembl_gene` / `uniprot` / `hmdb` / `kegg_compound` / `mfa_reaction` |
| `units` | MFA 必填，如 `nmol/min/g tissue` |
| `source` | 适配器名称与版本、外部工具版本、输入文件校验值 |

MFA 另有专用约定（不套用普通矩阵；实现见 `omicsweave.contract.mfa`，阶段 1 已完成）：

| 表 | 列（加粗为必填） | 说明 |
|---|---|---|
| `info` | **`template_version`、`model_id`、`flux_units`**、`software`、`tracer`、`notes` | 模型与单位写在这里，不逐行重复 |
| `flux` | **`sample_id`、`tissue`、`reaction_id`、`direction`、`value`**、`fit_id`、`ci_lower`、`ci_upper` | `direction` 为 `forward` / `reverse` / `net`。可逆反应分别给出正、逆向，净通量 = 正向 − 逆向，**按同一次拟合配对计算**；求解软件直接给净通量时用 `net`。`fit_id` 留空表示每个样本一个估计值 |
| `mid` | **`sample_id`、`metabolite`、`isotopologue`**、`tissue`、`measured`、`simulated`、`sd`、`fit_id` | 每行至少有 `measured` 或 `simulated` 之一；`simulated` 为最优拟合的模拟值 |
| `fits` | **`sample_id`、`fit_id`、`loss`** | 每个样本取 `loss` 最小的拟合作为代表值；没有 `fits` 表时取各次拟合的均值 |

- **样本与动物、分组的对应关系**不写在 MFA 表里，统一由 `samples.tsv`（`omics = mfa`）给出；`sample_id` 必须与之一致。
- **两种载体，列定义相同**：
  - Excel（`omicsweave_mfa_xlsx`）：工作表 `info` / `flux` / `mid` / `fits`，适合每个样本一个估计值（可附置信区间）或少量拟合。
  - Parquet 目录（`mfa_tables`）：`info.json` 加 `flux` / `mid` / `fits.parquet`。保留全部拟合结果时，行数通常超过 Excel 的 104 万行上限（约 150 次拟合 × 450 个通量列 × 30 个样本 ≈ 200 万行），应使用这种载体。
- **旧版自研软件格式**（`legacy_mfa_xlsx`）可以直接作为输入，`files` 依次为 `flux_raw_data.xlsx` 与可选的 `mid_raw_data.xlsx`，并在 `mfa.input` 中补充 `model_id` 与 `flux_units`；也可以用转换器一次性转为标准格式。转换结果与旧版 `MFAAnalyser` 的读取结果逐值一致（最优拟合通量、MID 矩阵、按拟合配对的净通量均有回归测试）。
- **命令**：`omicsweave mfa template` 生成 Excel 模板（含“说明”工作表）；`omicsweave mfa check <文件>` 检查结果文件，`--experiment <ID>` 时检查实验声明的 MFA 输入并与 `samples.tsv` 核对样本；`omicsweave mfa convert` 把旧版导出转为标准格式。

---

## 7. 数据契约的实现方式

- **磁盘格式**：Parquet 为正式格式（列类型可靠，支持列表型字段），`pyarrow` 列为基础依赖；JSON 用于元数据、manifest 与假设记录。
- **导出格式**：TSV / XLSX 只供人工查看，通过 `omicsweave export --format tsv|xlsx` 生成，不属于交换契约。Excel 另作为 MFA 的输入格式；它会自动改写基因名、日期等，且列类型不可靠，所以不作为结果的正式格式。
- **内存表示**：pandas DataFrame + 元数据字典。暂不自建数据对象；需要与 scverse 生态衔接时，再导出为 AnnData/MuData（H5AD/H5MU）。
- **校验**：`omicsweave.contract` 用 Pydantic 模型定义数据集、结果表与结果包，读写时检查必填列、ID 唯一性、`value_type` 与 `scale` 是否一致等；同一套模型导出 JSON Schema 供下游使用。各组学模块只通过这一层读写数据。
- 字段取值的受控词表（证据层次、验证状态、效应方向等）与缺失状态的区分，见总纲 3.3 节。

---

## 8. 输出规范

### 8.1 目录结构

每次运行生成一个结果包（AnalysisResultBundle），完整结构与 manifest 字段见总纲第 4 节：

```text
analysis_results/<experiment_id>/<run_id>/
├── manifest.json                     # 运行 ID、命令、配置哈希、包版本、输入校验值、输出清单
├── metadata/                         # subjects / samples（Parquet）、展开后的 experiment 配置
├── omics/<omics>/
│   ├── data/                         # 分析层数据集（第 6.2 节）
│   ├── qc/                           # QC 指标表与 QC 图
│   ├── tables/
│   │   ├── differential.parquet      # 统一列名的差异结果（所有比较 × 所有组织）
│   │   └── enrichment.parquet
│   └── figures/                      # <图类型>_<组织>_<比较>.pdf/.svg/.png
├── integration/<mode>/               # 跨组学整合结果（启用整合时）
├── evidence/                         # 标准化实体、证据记录、实体关系、同源映射
├── hypotheses/                       # 候选通路与候选靶点（启用整合时）
├── reports/                          # analysis_report.html（默认）、可选 PDF
├── qc/
└── provenance/                       # 软件版本、参数、RO-Crate 描述
```

单组学运行时只生成对应的 `omics/<omics>/` 和 `manifest.json`。

### 8.2 差异结果统一列名（`differential.parquet`）

| 列 | 说明 |
|---|---|
| `contrast_id`、`contrast_role`、`tissue`、`omics` | 对应 `experiment.yaml` 中的定义 |
| `feature_id`、`feature_name` | ID 类型见 `dataset.json` |
| `log2fc` | 分子组相对分母组 |
| `effect_size`、`effect_size_type` | 如 `hedges_g`，MFA 必填 |
| `statistic`、`pvalue`、`padj`、`padj_method` | |
| `n_numerator`、`n_denominator` | 实际参与检验的**生物学重复**数 |
| `min_achievable_p` | 秩检验在当前样本量下可达到的最小 p 值（沿用 MFA 模块已有的逻辑） |
| `method`、`flags` | 检验方法或证据路径（如 `deseq2`、`ttest`、`effect_size_only`）；标记（如 `low_n`、`underpowered`、`effect_size_undefined`） |
| `mean_numerator`、`mean_denominator`、`units` | 两组均值及其单位 |
| `inference_unit`、`evidence_level` | 推断单位（如 `biological_replicate`）；证据等级 `statistical_test` / `effect_size_only` / `insufficient` |
| `selected`、`selection_rule` | 是否入选，以及入选规则的文字描述 |

**MFA 的行另有必填列**：`effect_size`、`effect_size_type`、`mean_numerator`、`mean_denominator`、`units`、`inference_unit`、`evidence_level`。
- `pvalue` / `padj` 取主检验路径的值。样本量不足以达到显著时，证据路径为 `effect_size_only`，`pvalue` 为空，`padj_method` 为 `none`。
- 效应量只在两种情况下允许为空：`evidence_level = insufficient`，或行上标记了 `effect_size_undefined`（例如所有动物都贴在约束边界，两组方差均为 0，Hedges' g 为 0/0）。
- 可辨识性、贴边比例、稳定性校准等 MFA 专有指标为可选列：`identifiability`、`frac_at_bound`、`selection_probability`、`direction_stability`、`permutation_pvalue`。

### 8.2b 富集结果统一列名（`enrichment.parquet`）

比较层面的通路富集结果，ORA 与预排序 GSEA 共用一张表。样本层面的 ssGSEA 得分属于数据集（行为通路、列为样本），不放在这里。

| 列 | 说明 |
|---|---|
| `contrast_id`、`contrast_role`、`tissue`、`omics` | 同差异结果表 |
| `method`、`query` | `ora`（`query` 为所用特征列表：`up` / `down` / `all`）或 `gsea_prerank`（`query` 为 `ranked`） |
| `gene_set_library`、`term_id`、`term_name` | 基因集库、条目 ID（如 `ko00010`、`GO:0006096`，可空）、条目名称 |
| `pvalue`、`padj`、`padj_method` | GSEA 自带的 FDR q 值记为 `gsea_fdr` |
| `effect_direction` | `up` / `down` / `mixed` / `none`：ORA 由查询列表决定（`all` → `mixed`），GSEA 由 NES 的符号决定 |
| `overlap_size`、`set_size` | ORA 必填：命中数与通路大小；GSEA 为 leading edge 数与通路大小 |
| `query_size`、`universe_size`、`odds_ratio` | ORA 可选 |
| `es`、`nes` | GSEA；`nes` 必填 |
| `hit_features`、`tool`、`flags` | 以分号分隔的命中特征；产生结果的工具 |

**从现有输出转换**（阶段 2 各组学改为直接写出统一表之前使用）：
- `omicsweave.mfa.export.differential_flux_to_contract`：MFA 差异通量宽表；
- `omicsweave.adapters.gseapy_results.enrichr_to_contract` / `prerank_to_contract`：转录组、蛋白组的 gseapy 结果；
- `omicsweave.metabolomics.export.ora_to_contract`：代谢组 ORA。

测试中这些转换器的输入来自真实函数：MFA 的 `differential_flux`、代谢组的 `calculate_pathway_enrichment`、gseapy 的 `prerank`。

候选通路与候选靶点表（`hypotheses/`）的字段见总纲 3.5 节。

### 8.3 图件、报告与溯源

- **图件**：矢量格式（PDF/SVG）+ PNG（300 dpi）；统一主题；配色由 `experiment.yaml` 中的因素水平决定。
- **报告**：HTML（默认）、PDF（可选），由结果包生成，只陈述统计事实与方法参数，**不写入课题专用的生物学结论**；候选清单单独成节，标注“假设，需实验验证”。
- **溯源（`provenance/` 与 `manifest.json`）**：包版本、Git 提交号（如可获取）、配置哈希、输入文件校验值、外部工具版本、数据库与同源映射版本、随机种子、Python 与主要依赖的版本。

---

## 9. 运行模式

### 9.1 单组学独立运行

```bash
omicsweave validate --experiment 20250613_dio_glp1
omicsweave run --experiment 20250613_dio_glp1 --omics mfa
omicsweave run --experiment 20250613_dio_glp1 --omics all        # 依次运行所有已启用的组学，彼此独立
omicsweave run --experiment <id> --omics transcriptomics --steps qc,differential
```

### 9.2 多组学整合（可选，三种模式）

| 模式 | 前提条件 | 内容 | 现有代码基础 |
|---|---|---|---|
| `flux_anchored` | MFA + 至少 1 个其他组学；有 MFA 反应→基因/蛋白/代谢物映射 | 以差异通量为锚点，汇总其他组学的方向一致性证据（多维证据画像） | `evidence_profile.py`、`mappings/`、`stability.py` |
| `pathway_convergence` | 任意 ≥2 个组学；不要求 MFA | 通路层面的跨组学、跨组织收敛分析 | `pathway_convergence.py`、`cross_tissue_pathway_scorer.py` |
| `sample_level` | ≥2 个组学测自同一批 `subject_id`（配对） | 样本级协调、残余批次诊断，以及潜变量/网络等补充分析 | `harmonization.py` |

另有项目级的**多实验随机效应 Meta 分析**（`project_integration.py`、`batch_effect_meta.py`），用于汇总多个可比实验的结果。

```bash
omicsweave integrate --experiment <id> --mode flux_anchored
omicsweave integrate --project --mode meta                       # 多个实验的 Meta 分析
```

不满足前提条件时，`validate` 和 `integrate` 会明确说明缺少什么（例如“flux_anchored 需要 MFA 数据，当前实验未启用 MFA”），而不是静默跳过。

---

## 10. 现有代码整合清单

依据模块之间的导入关系与实际用途整理（行数为当前值）。下文沿用改名前的路径：阶段 1 已把 `common/` 改为 `core/`（统计模块在 `core/statistics/`，ID 映射在 `core/idmap/`），`multiomics/` 改为 `integration/`，其中的 `multiomics/integration/` 子包改为 `integration/methods/`。

### 10.1 保留，迁入 `core` / `visualization` 并合并

| 现有代码 | 去向 | 说明 |
|---|---|---|
| `common/nonparametric.py`、`meta_analysis.py`、`stability.py` | `core/statistics/` | 推断单位、秩检验 p 值下限等方法学核心 |
| `common/pandas_safety.py`、`common/project.py` | `core/` | |
| 3 份 `provenance.py`（mfa / transcriptomics / multiomics） | `core/provenance.py` | ✅ MFA 与转录组的读写、比对已合并，各自只保留本组学的字段；multiomics 的那份实为输入审计，改名为 `integration/input_audit.py` |
| 5 套配置加载器（约 3,450 行） | `core/config.py` + 各组学设置段 | 改为第 5 节的统一模型；旧配置通过 `omicsweave migrate` 转换 |
| 5 份 `plotting/_config.py`；3–4 份火山图、热图、PCA | `visualization/` | 主题 ✅（`visualization/theme.py`）；通用图在阶段 2 合并，组学专用图保留在各模块；`_config.py` 中的全局状态在 1.10 去除 |

### 10.2 保留为各组学核心（接入新数据契约）

- **MFA**：`df_analysis`、`mfa_analysis`、`infusion_qc`、`bridging_qc`、`mfa_integration`、`plotting/`
- **转录组**：`normalization`、`run_wilcoxon_de`、`qc_gate`、`joint_model`、`transcriptomics_integration`；`rnaseq_pipeline`（1,480 行）拆分为“原始层适配器（QC、修剪、kallisto）”和“分析层（差异分析）”两部分
- **蛋白组**：`differential_analysis`（limma）、`proteomics_qc`、`enrichment/`；`diann_search` 改为适配器
- **代谢组**：`data_preprocessing`、`differential_analysis`、`enrichment/`、`gsea/`
- **整合**：`evidence_profile`、`harmonization`、`project_integration`、`batch_effect_meta`、`pathway_convergence`、`cross_tissue_pathway_scorer`（去掉硬编码的 `DEFAULT_TISSUES`）、`mappings/feature_resolver`、`audit`；`loaders/` 改为只读取数据契约

### 10.3 重写

| 现有代码 | 改为 |
|---|---|
| ID 映射字典：`metabolite_id_mapping`（2,094 行）、`mouse_gene_id_mapping`（522）、`mouse_protein_id_mapping`（273）、`reaction_pathway_mapping`（1,374） | 带版本和来源说明的数据文件（TSV）+ 解析器接口；按 `project.species` 选择；MFA 反应注册表改为网络模型定义文件 |
| 4 份解读报告生成器 | 基于模板的报告，只写统计事实，删除课题专用的叙述 |
| `multiomics/integration/flux_omics_scorer.py`（1,632 行） | 拆成两部分（第 13 节已定）：**删除**旧的加和总分（0–4.0 分、固定权重、High/Medium/Low 置信度及全部 `legacy_*` 输出）；**保留并重构**各组学的证据提取（反应→基因/蛋白/代谢物映射、方向与效应量提取），改名为 flux-anchored 证据构建器，只为多维证据画像提供输入 |

### 10.4 移出到本地课题工作区（课题专用）

以下代码迁移到 `local/projects/<课题>/analysis/`（本地，不提交；见仓库中的 `local/README.md`），并在实验项目仓库保留正式副本。迁移在课题数据的本地副本就位、且这些脚本在实验项目仓库中有备份之后进行，在此之前保留在包内：

- `multiomics/tissue_pathway_panels.py`、`plot_tissue_pathway_panels.py`、`validate_tissue_pathway_panels.py`（约 1,240 行，六组织图件与脑 HFD/LFD 参照）
- `multiomics/run_cross_tissue_integration.py`（依赖人工整理的 `cross_tissue_evidence.yaml`）

另外，`multiomics/pipeline_evaluation_report.md`（已标注为历史文档）移到 `docs/history/`。

### 10.5 删除

原则：已失效或已有替代的立即删除；仍提供实际功能的脚本，等替代命令就绪后再删。

| 文件 | 原因 | 状态 |
|---|---|---|
| `transcriptomics/rerun_enrichment.py` | 已无法运行（导入的函数已不存在） | ✅ 已删除 |
| `transcriptomics/run_rna_seq.sh` | 自动生成的 STAR 脚本，路径是占位符，代码中没有使用 | ✅ 已删除 |
| `proteomics/run_proteomics_analysis.py` | 已标注弃用的兼容入口 | ✅ 已删除 |
| 兼容转发模块：`multiomics/integration/meta_analysis.py`、`stability.py`，`mappings/reaction_registry.py`、`gene_annotations.py`、`metabolite_annotations.py`、`protein_annotations.py` | 目前没有外部用户，已改为新路径导入 | ✅ 已删除 |
| `transcriptomics/rerun_de_analysis.py`、`finish_gsea.py`、`summarize_gsea.py`、`run_enrichment.py`、`run_gsea.py` | 一次性补跑或薄包装脚本，由 `omicsweave run --steps` 取代 | 阶段 2 替代命令就绪后删除 |
| `transcriptomics/check_markers.py`、`qc_one_sample.py` | 调试脚本 | 随课题专用脚本移到本地工作区（第 10.4 节） |
| `transcriptomics/setup_kallisto.py`、`download_transcriptome.py` | 只适用于 Windows 和小鼠；SSL 校验已恢复默认开启 | 阶段 3 由外部工具适配器取代后删除 |

---

## 11. 实施路线

工期按 1 名主要开发者估算，仅供排期参考。每个阶段结束时打一个版本标签。本节阶段 0–4 对应总纲中的 O0–O4；TheraWeave 的阶段（T0–T3）及两者之间的前置条件见总纲第 8 节。

### 阶段 0：工程化打包 ✅（2026-10-29 完成）

可安装、可导入、统一命令行、测试与 CI 配置、大文件移出仓库。

### 阶段 1：核心契约与项目模型（约 6–8 周）→ v0.2（进行中）

| 任务 | 交付物 | 状态（2026-10-29） |
|---|---|---|
| 1.1 改名：`common`→`core`，`multiomics`→`integration`；删除兼容转发模块 | 新包结构 | ✅ `core` 下分出 `statistics/`、`idmap/`、`safety.py`；`integration` 内原 `integration/` 子包改名为 `methods/` |
| 1.2 项目模型与配置校验（`core/model.py`、`config.py`），含比较角色与逆转配对 | `project.yaml` / `experiment.yaml` / 样本表的规范与校验 | ✅ 含 `core/validate.py` 全部第 5.6 节规则 |
| 1.3 数据契约（`omicsweave.contract`）：`AnalysisRun`、`AnalysisResultBundle` 最小字段 | 分析层数据集、结果表与结果包的读写和校验（Parquet） | ✅ 数据集、差异结果表（含 MFA 专有必填列）、富集结果表、运行清单（草案 0.1），以及从现有输出转换的函数；见第 8.2 节 |
| 1.4 合并溯源记录、合并绘图主题 | `core/provenance.py`、`visualization/theme` | ✅ 另把 `contract` 与蛋白组 DIA-NN 搜库中的文件哈希也改用 `core.provenance.sha256_file` |
| 1.5 `omicsweave init / validate / migrate` | 新建项目模板；把现有实验的旧配置转换为新格式 | 部分：`init`、`validate` 已完成；`migrate` 需要真实课题配置做对照，等本地副本就位 |
| 1.6 MFA Excel 模板 + 自研软件格式转换器 | 模板文件、校验命令 | ✅ `omicsweave mfa template / check / convert`；标准格式见第 6.2 节 |
| 1.7 按第 10.5 节删除废弃代码；课题数据本地副本就位后，按第 10.4 节把课题专用脚本移到 `local/projects/<课题>/analysis/` | 代码量明显下降 | 部分：已失效的脚本与兼容模块已删除，其余见第 10.5 节 |
| 1.8 合成示例项目（小型、可公开） | `examples/demo_project/`，CI 中端到端运行 | ✅ CI 中运行 `validate` |
| 1.9 修复 SSL 校验；补 CITATION.cff（LICENSE 已于阶段 0 加入） | | ✅ |
| 1.10 运行函数改为接受显式项目上下文，逐步去除进程级全局状态 | 为阶段 5 预留 | 未开始 |

**验收标准**
- 示例项目能通过 `validate`。
- 用 `migrate` 转换至少一个现有实验后，MFA 结果与迁移前一致（数值逐项比对）。
- 测试全部通过。

### 阶段 2：各组学接入 + 运行模式 + 证据层与交换契约 v1（约 6–8 周）→ v0.3

| 任务 | 交付物 |
|---|---|
| 2.1 四个组学改为“读取数据契约 → 统一结果表”，分析层导入器：kallisto/salmon 定量、DIA-NN 报告、峰面积表/MS-DIAL 输出、MFA 模板 | 统一的 `differential.parquet` / `enrichment.parquet`；`omicsweave export` |
| 2.2 `omicsweave run --omics ... --steps ...`；旧命令保留为别名 | 单组学独立运行 |
| 2.3 整合模块只读取契约；实现三种整合模式及前提检查 | `omicsweave integrate --mode ...` |
| 2.4 ID 映射改为数据文件 + 解析器；物种可配置；加入小鼠→人同源映射 | `core/idmap/`、`evidence/orthologs.parquet` |
| 2.5 证据层：`BiologicalEntity` / `Association` / `Evidence` / `PathwayCandidate` / `TargetHypothesis`；候选清单含逆转标记与证据等级 | `evidence/`、`hypotheses/` |
| 2.6 报告改为由结果包生成的 HTML（含自动方法学段落）；按需 PDF | `reports/analysis_report.html` |
| 2.7 交换契约 v1：发布 JSON Schema 与标准样例结果包，CI 中运行契约测试 | `schemas/evidence-exchange/v1/` |
| 2.8 删除旧加和总分，`flux_omics_scorer` 的证据提取并入证据画像（第 10.3 节）；需要排序时，使用可配置、写明依据的排序规则，并做权重敏感性分析 | 证据画像为唯一整合结果 |

**验收标准**
- 示例项目能在单组学和三种整合模式下端到端运行。
- 用现有实验数据回归测试，主要结果与迁移前一致；如有差异，需能逐项说明原因。
- 任一候选通路/靶点可回溯到具体统计结果、输入样本和来源数据库。
- 标准样例结果包能通过 JSON Schema 校验（TheraWeave 的 T1 以此为前置条件）。

### 阶段 3：原始数据适配器（约 8–12 周）→ v0.4

| 任务 | 说明 |
|---|---|
| 3.1 转录组：`fastq.gz` → cutadapt → kallisto/salmon | 从 `rnaseq_pipeline` 拆出；跨平台 |
| 3.2 蛋白组：Thermo `.raw` → DIA-NN；`.raw` → mzML | 从 `diann_search` 演进 |
| 3.3 代谢组：`.wiff`/`.raw` → msconvert → MS-DIAL/MZmine（非靶向）或 Skyline（靶向）结果导入 | 优先“导入外部软件结果”，其次才是自动调用 |
| 3.4 同位素示踪：天然丰度校正 → MID 表 | 衔接 MFA |
| 3.5 （可选）Snakemake/Nextflow 封装，用于服务器/集群上的大批量原始数据处理 | 视实验室计算环境决定 |

**验收标准**：每种原始格式都有一个公开的小样例，能从原始文件跑到分析层数据集。

### 阶段 4：对外发布与论文（约 4–6 周）→ v1.0

- 英文文档、教程 notebook、API 文档。
- 用公开的 13C 示踪 + 多组学数据集做验证与基准比较。
- 发布到 PyPI、Zenodo 存档获取 DOI。
- 投稿 JOSS 或 Bioinformatics Application Note。

### 阶段 5（v1.0 之后）：集群部署与多客户端联网

按需求出现的先后实施，依赖第 4.1 节预留的约束：
- **多节点集群**：原始数据处理与大批量分析交给作业调度器（如 Slurm）或工作流引擎（Snakemake/Nextflow）执行，项目目录放在共享存储上。
- **多客户端联网使用**：增加服务端（任务提交、排队、状态查询的 HTTP 接口）、用户认证与项目级权限，客户端只提交配置、查看结果包。
- 前提：阶段 1–2 已去除进程级全局状态，所有步骤都能以“配置 + 输入文件 → 结果包”的方式独立运行。

---

## 12. 风险与对策

| 风险 | 对策 |
|---|---|
| 原始数据处理工作量远超预期（每种组学都有成熟的多年项目） | 只做适配器，不重新实现；排在阶段 3，不阻塞前面阶段的交付 |
| 厂商格式与平台绑定（Windows、闭源库） | 以开放格式（mzML）为中间层；分析层导入始终可用 |
| 外部工具许可限制 | 不随包分发，文档说明安装方式；在 `validate` 中检查并记录版本 |
| 迁移导致现有课题结果发生变化 | `migrate` + 数值回归测试；旧结果目录保持只读 |
| 小团队维护负担 | 严格控制范围；可选依赖分组；CI 覆盖最小安装 |
| 课题未发表数据外泄 | 示例只用合成或公开数据；公开前清理代码中的实验编号 |

---

## 13. 决策记录

2026-10-29 确认：

| 事项 | 决定 | 状态 |
|---|---|---|
| 配置校验库 | Pydantic v2：报错信息清晰，可导出 JSON Schema 作为配置文档与交换契约 | 已定 |
| 结果正式格式 | Parquet，TSV/XLSX 只作导出（与总纲 4.1 节一致） | 已定 |
| 首批支持的物种 | 小鼠（`mus_musculus`，taxon 10090）+ 人（`homo_sapiens`，taxon 9606）；ID 映射数据文件与小鼠→人同源映射都按这两个物种准备，其他物种以后通过数据文件扩展 | 已定 |
| 开源许可证 | **BSD-3-Clause**（OmicsWeave 与 TheraWeave 相同）：与主要依赖兼容，并含“未经许可不得用作者或单位名义为衍生产品背书”条款，与 scverse 生态（scanpy、anndata）一致。版权方为 “The OmicsWeave Developers”（集体署名；贡献者名单在 `CITATION.cff` 中维护） | 已定；`LICENSE` 与 `pyproject.toml` 许可证字段已加入 |
| 运行环境 | 当前实现 Windows/Linux 单机运行；以后扩展多节点集群部署与多客户端联网使用（阶段 5）。从阶段 1 起遵守第 4.1 节的预留约束 | 已定 |
| 旧打分矩阵 `flux_omics_scorer` | 尚未发表，且加和总分的策略过于粗糙：**删除加和总分及 `legacy_*` 输出，不保留 `--legacy`**。不另起一套新的总分：多维证据画像（`evidence_profile`，现已是主结果）作为唯一整合结果，按维度给出证据等级，并用 `stability` 做生物学稳定性校准；各组学的证据提取部分保留并重构（第 10.3 节） | 已定，阶段 2 实施 |
| 课题专用脚本与数据 | 暂不移动外部项目的原始数据。本仓库预留 `local/` 本地课题工作区（不提交，已建立）：课题数据副本放 `local/projects/<课题>/`，课题专用脚本放其中的 `analysis/`；通过 `--project-dir`、环境变量 `OMICSWEAVE_PROJECT_DIR` 或 `pytest --project-dir` 指定。包内现有的课题专用模块在本地副本就位后再迁移（第 10.4 节） | 已定，工作区已就绪 |
| 阶段 2 之后先做阶段 3 还是 TheraWeave MVP | 见总纲第 11 节；现有课题数据已处理成表格，倾向先做 TheraWeave MVP | 待定 |
| 文档语言 | 中文先行，阶段 4 补英文 | 已定 |

---

## 14. 近期执行清单（阶段 1 启动）

- [x] 确认第 13 节的决策事项（许可证：BSD-3-Clause）
- [x] 加入 `LICENSE`（版权方 The OmicsWeave Developers）与 `pyproject.toml` 许可证字段
- [x] 建立本地课题工作区 `local/` 与 `pytest --project-dir` 选项
- [ ] 把课题数据副本拷贝到 `local/projects/<课题>/`，用 `pytest --project-dir` 跑通需要真实数据的 3 个测试
- [x] 提交阶段 0 的改动
- [ ] 在 GitHub 上跑通 CI
- [x] 1.1 包改名、删除兼容转发模块
- [x] 1.2 `project.yaml` / `experiment.yaml` / 样本表规范与 `omicsweave validate`
- [x] 1.3 数据契约：数据集、差异结果表（含 MFA）、富集结果表、运行清单
- [x] 1.4 合并溯源记录与绘图主题
- [ ] 1.5 `omicsweave migrate`（需课题数据本地副本）
- [x] 1.6 MFA Excel 模板与自研软件格式转换器
- [ ] 1.7 课题专用脚本移到本地工作区（需本地副本）
- [x] 1.8 制作合成示例项目
- [x] 1.9 修复 SSL 校验、补 CITATION.cff
