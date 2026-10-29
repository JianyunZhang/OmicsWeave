# OmicsWeave：多组学 Python 分析工具包架构规划

> 文档类型：项目架构与实施规划  
> 项目名称：**OmicsWeave**  
> 英文定位：*A Modular Python Toolkit for Multi-Omics Data Analysis and Integration*  
> 当前建议：**Monorepo（单 Git 仓库）+ 单 Python 分发包 + 组学功能模块化 + 可选依赖**  
> 文档日期：2026-10-29

## 1. 项目定位与范围

OmicsWeave 旨在构建一套长期维护、可复用、可扩展的 Python 科研分析工具，用于多组学数据处理、统计分析、结果可视化及跨组学整合。

### 1.1 计划支持的分析方向

| 模块 | 主要工作范围 | 典型输入 | 典型输出 |
| --- | --- | --- | --- |
| `transcriptomics` | Bulk RNA-seq 表达矩阵处理、质控、差异表达、富集分析；未来可扩展单细胞/空间转录组 | Counts/TPM 矩阵、样本元数据、差异分析结果 | QC、DEG 表、PCA、火山图、通路富集结果 |
| `metabolomics` | 靶向/非靶向代谢组学数据清洗、注释、归一化、统计、同位素数据辅助处理 | 峰面积表、代谢物注释、样本信息 | 差异代谢物、统计表、富集与可视化结果 |
| `mfa` | 代谢流分析结果处理；MID 分析、拟合质量评价、通量表及网络可视化 | 同位素分布（MID）、通量估计、模型拟合结果 | 通量估计与区间、残差、MID 拟合图、代谢网络图 |
| `proteomics` | 蛋白质定量结果处理、缺失值处理、归一化、差异表达与功能分析 | 蛋白质丰度矩阵、鉴定结果、样本元数据 | 差异蛋白、统计结果、功能富集与可视化 |
| `integration` | 多组学结果标准化与关联分析、联合通路/网络分析 | 各组学的标准化结果、对应样本与特征映射 | 跨组学关联、网络、联合特征与综合报告 |
| `visualization` | 各模块可复用的图表和导出规范 | 标准结果对象、DataFrame | 静态/交互式图、出版级图表 |
| `core` | 公共数据模型、I/O、配置、样本元数据、基础统计、通用预处理 | 通用表格、元数据、配置 | 可复用底层能力与统一接口 |

> **重要边界**：`mfa` 与 `metabolomics` 应保持独立。MFA 的通量、置信区间、MID、拟合残差及代谢网络，与代谢组学的峰面积和丰度数据语义不同。两者可共享基础数据处理函数，但不应被强制塞入同一种结果对象。

## 2. 统一命名方案

### 2.1 品牌名称

推荐名称：**OmicsWeave**。

选择理由：

- **Omics** 指向广义组学分析，而非仅限转录组或代谢组。
- **Weave** 表达不同组学信息的连接、组合和交叉验证，适合多组学整合。
- 适合作为 Git 仓库名称、Python 包名称及未来系列工具的统一前缀。
- 可扩展空间代谢组、同位素示踪和质谱成像分析，不必重新命名项目。

备选名称：`OmicsFoundry`（侧重工具集合与数据处理）、`OmicsLattice`（侧重多维数据结构与系统生物学）。

> 发布前仍需分别检查 GitHub 仓库、PyPI 分发包、域名和相关商标的名称占用情况；本文件不代表名称可用性已获确认。

### 2.2 系列名称与功能映射

| 系列名称（展示层） | 推荐 Python 模块名 |
| --- | --- |
| OmicsWeave-RNA | `omicsweave.transcriptomics` |
| OmicsWeave-Metabo | `omicsweave.metabolomics` |
| OmicsWeave-Flux | `omicsweave.mfa` |
| OmicsWeave-Proteo | `omicsweave.proteomics` |
| OmicsWeave-Integrate | `omicsweave.integration` |

**注意**：上述系列名称是产品/模块命名示例，并不意味着第一阶段需要建立多个 Git 仓库或发布多个 PyPI 包。

## 3. 仓库组织决策

### 3.1 推荐结论

**第一阶段采用 Monorepo：一个 `OmicsWeave` Git 仓库、一个 `omicsweave` Python 分发包、多个边界清晰的功能模块。**

具体策略：

1. 统一进行代码审查、Issue、版本管理、测试、文档与示例维护。
2. 通过顶层模块划分不同组学分析能力，避免不同学科的代码相互缠绕。
3. 公共数据模型与工具放入 `core`，跨组学整合逻辑放入 `integration`。
4. 通过可选依赖（optional dependencies / extras）减少不必要的安装负担。
5. 在模块 API 稳定、依赖或发布周期显著分化后，优先考虑**单仓库多分发包**，最后才考虑**多仓库**。

### 3.2 Monorepo 与 Polyrepo 对比

| 维度 | Monorepo（统一 Git 仓库） | Polyrepo（每个模块独立 Git 仓库） |
| --- | --- | --- |
| 初期维护 | 仓库、CI、文档统一，管理成本较低 | 配置和维护重复较多 |
| 跨模块重构 | 一次提交可调整多个模块 | 需协调多个版本及仓库 |
| 公共代码复用 | 容易直接共享 `core` | 通常需要独立的公共包 |
| 发布节奏 | 可以统一，也可在成熟后按包独立发布 | 天然独立 |
| 依赖隔离 | 依赖设计与环境隔离需专门处理 | 每仓库可独立配置环境 |
| 多组学整合 | 便于同步修改接口 | 需要管理跨仓库兼容性 |
| 多团队协作 | 适合初期小团队；后续需定义边界 | 适合独立维护团队 |

**判断**：现阶段优先控制科研软件的工程维护开销，不为了“看起来模块化”而过早拆仓库。

## 4. 第一阶段建议目录结构

```text
OmicsWeave/
├── pyproject.toml                  # 单分发包、依赖组与构建配置
├── README.md                       # 项目概况、安装、快速开始
├── LICENSE
├── .gitignore
├── src/
│   └── omicsweave/
│       ├── __init__.py
│       ├── core/                   # 基础数据模型、I/O、配置、公共统计工具
│       │   ├── data_models/
│       │   ├── io/
│       │   ├── preprocessing/
│       │   └── statistics/
│       ├── transcriptomics/        # RNA-seq
│       │   ├── qc/
│       │   ├── differential/
│       │   └── enrichment/
│       ├── metabolomics/           # 靶向/非靶向代谢组学
│       │   ├── targeted/
│       │   ├── untargeted/
│       │   └── isotope/
│       ├── mfa/                     # MFA 结果分析
│       │   ├── mid_analysis/
│       │   ├── flux_analysis/
│       │   └── visualization/
│       ├── proteomics/
│       ├── integration/
│       └── visualization/          # 跨模块通用绘图
├── tests/
│   ├── core/
│   ├── transcriptomics/
│   ├── metabolomics/
│   ├── mfa/
│   ├── proteomics/
│   └── integration/
├── examples/                       # Notebook、演示数据与工作流程
├── docs/                           # API 文档、方法学及用户指南
└── .github/
    └── workflows/                 # lint、test、build、release
```

目录树表达的是**逻辑模块划分**，并非要求第一天就实现全部空目录；优先创建有实际功能及测试的模块。

## 5. 关键软件架构原则

### 5.1 依赖方向

推荐依赖关系：

```text
transcriptomics ─┐
metabolomics ────┤
mfa ─────────────┼──> core
proteomics ─────┘
         │
         └──> integration（消费各模块的标准化输出）
```

更准确地说：

- `core` **不依赖** `transcriptomics`、`metabolomics`、`mfa`、`proteomics` 或 `integration`。
- 组学模块可以依赖 `core`，但尽量避免相互直接引用内部实现。
- `integration` 使用**约定的公开输出接口**和特征/样本映射进行整合。
- `visualization` 的公共能力应尽量与特定组学解耦；领域专属绘图保留在相应模块中。
- `integration` 不应成为存放所有杂项分析代码的“万能模块”。

### 5.2 标准化数据接口

跨组学整合的重点是**接口兼容**，而非强行统一所有内部数据模型。

建议为每类分析输出提供：

- **样本元数据**：样本 ID、实验组、批次、组织/细胞类型、条件等；显式检查跨组学样本匹配。
- **特征标识**：基因 ID、蛋白质 ID、代谢物标识、反应 ID；保留映射来源与版本。
- **数值语义**：区分 counts、归一化丰度、log fold change、MID 分数、通量单位、统计显著性等。
- **结果溯源**：分析方法、参数、软件/数据库版本、处理日期与数据来源。
- **导出接口**：优先提供规范 DataFrame、Parquet/CSV 和可序列化元数据。

对 MFA 特别保存通量单位、模型标识、反应方向、拟合结果、置信区间与 MID 误差定义，以防与普通代谢物丰度混用。

## 6. 依赖管理与安装策略

### 6.1 第一阶段：单分发包 + Optional Dependencies

可使用 Python 项目的 `pyproject.toml` 定义不同功能的 extras，提供以下预期安装方式：

```bash
pip install "omicsweave[transcriptomics]"
pip install "omicsweave[metabolomics]"
pip install "omicsweave[mfa]"
pip install "omicsweave[proteomics]"
pip install "omicsweave[integration]"
```

> 以上是**规划中的安装接口**，不是当前已发布的安装命令。实际运行需先在 `pyproject.toml` 中声明对应 extras，并完成本地构建或发布。

建议：

- 把真正通用、轻量的包放入基本依赖。
- 将重型计算、绘图和领域专属工具放在相应 extras 中。
- CI 至少分别测试核心最小环境、各模块环境和典型联合环境。
- extras **不会**自动建立相互隔离的虚拟环境；不同依赖版本冲突时，应显式设计环境边界。
- 若涉及 R/Bioconductor 或外部命令行工具，明确其安装要求与版本，而不是假定 `pip` 能解决全部环境配置。

### 6.2 第二阶段：单仓库、多 Python 分发包

当模块成熟且具有稳定独立 API 后，可迁移到：

```text
OmicsWeave/
├── packages/
│   ├── omicsweave-core/
│   ├── omicsweave-rna/
│   ├── omicsweave-metabo/
│   ├── omicsweave-mfa/
│   ├── omicsweave-proteo/
│   └── omicsweave-integration/
├── examples/
└── docs/
```

每个子包拥有独立的 `pyproject.toml`，可以独立构建、测试与发布。可考虑 `uv workspace` 等工具管理代码库，但需要注意：**workspace 默认共享锁文件和工作环境，并不等同于依赖完全隔离**。若各分发包需互不兼容的运行时，应保留独立虚拟环境。

**迁移原则**：对外暴露稳定的导入路径与接口，尽可能减少从单分发包切换到多包后用户代码的修改。

## 7. 何时拆分为多个 Git 仓库？

建议不要按研究方向立即拆仓库；只有当下列条件中的多项成立时再考虑：

1. **开发与维护者独立**：模块有专门团队和独立的 Issue、Roadmap、评审流程。
2. **发布周期明显不同**：例如 MFA 模块需要高频发布，其余模块长期稳定。
3. **依赖/运行平台显著分化**：单仓库下的环境与 CI 维护已明显成为负担。
4. **科研成果独立化**：某一算法已形成可独立引用的软件、方法学论文或公共社区。
5. **公共接口已经稳定**：可以通过明确的包 API 和数据契约交换结果，而无需频繁跨仓库联动修改。

**注意**：单纯的依赖冲突通常优先通过包与运行环境隔离解决，不是自动拆分 Git 仓库的理由。

## 8. 分阶段实施路线

### 阶段 1：统一仓库 + 模块化单包（当前推荐）

**优先任务**：

- 初始化 `OmicsWeave` 仓库、`src/` 布局、`pyproject.toml`、测试和 CI。
- 实现 `core` 的样本元数据、特征元数据、数据读写与结果格式。
- 优先开发现有实验最需要的 `transcriptomics`、`metabolomics` 和 `mfa` 分析功能。
- 建立通用可视化风格与分析 Notebook 示例。
- 对每个模块建立可复现的小型测试数据和基础单元测试。

**阶段完成标准**：至少一个端到端分析流程可重复运行；核心数据结构有文档；各模块之间没有不必要的循环依赖。

### 阶段 2：统一仓库 + 多分发包（视需求）

**触发条件**：模块功能成熟、独立版本发布有实际价值、依赖集合开始明显分化。

**优先任务**：

- 迁移到 `packages/` 结构，提取 `omicsweave-core`。
- 为各模块独立定义构建配置、依赖和测试矩阵。
- 建立版本兼容策略及模块 API 稳定性约定。
- 保持跨组学分析示例与统一文档入口。

### 阶段 3：按需拆分独立仓库

**触发条件**：独立团队、独立研究成果或发布/治理需求已经超过 Monorepo 的收益。

**优先任务**：

- 仅拆分具备独立维护条件的成熟模块。
- 保持统一品牌、公开数据协议和兼容性测试。
- 维护 `OmicsWeave` 作为生态总览、集成示例或元包入口。

## 9. 近期执行清单

- [ ] 确认 `OmicsWeave` 品牌名称与仓库/PyPI 命名可用性。
- [ ] 创建单个 `OmicsWeave` 仓库，启用 `src/` 项目结构。
- [ ] 定义 `core` 的样本、特征、结果及溯源数据约定。
- [ ] 确定 RNA-seq、代谢组学、MFA 三个首批最小可用分析流程。
- [ ] 设计并声明各组学模块的 optional dependencies。
- [ ] 加入 `pytest`、格式检查、CI 与最小示例数据集。
- [ ] 完成各模块的第一版结果导出规范，并验证 `integration` 能消费这些结果。
- [ ] 在明确需要独立发布前，不急于拆包或拆仓库。

## 10. 最终决策摘要

**建议当前采用：**

- **品牌**：`OmicsWeave`
- **Git 仓库**：单一 `OmicsWeave` 仓库（Monorepo）
- **Python 包**：初期一个 `omicsweave` 分发包
- **功能划分**：`core`、`transcriptomics`、`metabolomics`、`mfa`、`proteomics`、`integration`、`visualization`
- **依赖管理**：公共轻量依赖 + 组学模块 extras；必要时隔离运行环境
- **长期扩展**：先考虑单仓库多分发包，再按实际维护需求拆为多个 Git 仓库

> 核心原则：**先确保模块边界、数据接口和可复现性，再决定仓库数量。**
