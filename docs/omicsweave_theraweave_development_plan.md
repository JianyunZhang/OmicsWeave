# OmicsWeave × TheraWeave：多组学分析与治疗发现系统架构规划

> **文档类型**：软件系统架构与研究开发规划（两个系统的总纲）  
> **版本**：v1.1（修订稿，2026-10-29；修订自 v1.0 讨论稿，修订内容见第 0 节）  
> **适用范围**：Python 多组学分析、代谢流分析、证据整合、候选靶点发现与下游治疗策略研究  
> **项目名称**：`OmicsWeave`（上游） / `TheraWeave`（下游，暂定；初步检索未发现重名，正式采用前需完成商标核查，见第 10 节）  
> **配套文档**：OmicsWeave 的项目配置、输入格式与运行模式等详细设计见 [omicsweave_arch_design_v2.md](omicsweave_arch_design_v2.md)。两份文档冲突时，**结果数据格式与交换契约以本文为准**，其余以 v2 为准。

## 0. 修订说明（v1.0 → v1.1）

| 修订内容 | 位置 |
|---|---|
| 新增当前进度：OmicsWeave 阶段 0 已完成；TheraWeave 名称的初步检索结果 | 1.2、10 |
| 新增“实验设计语义”：比较的角色（疾病效应 / 治疗效应）与“治疗逆转”标记，这是产生治疗假设的关键信息 | 2.3 |
| 实体模型补充“反应（Reaction）”类型、Biolink 类别、小鼠→人同源映射 | 3.1、3.6 |
| 新增受控词表：证据层次、验证状态、干预方向等字段的取值范围 | 3.3 |
| 候选清单字段补充：同源基因、比较角色、逆转、生物学重复数、代谢物鉴定等级、证据等级 | 3.5 |
| 结果包路径与 v2 的项目目录统一；明确 Parquet 为正式格式，TSV/XLSX 只作导出 | 4 |
| 交换契约补充：Schema 的存放位置、版本规则、共享的契约测试样例；TheraWeave 不需要安装 OmicsWeave 即可读取 | 5 |
| TheraWeave 补充：定位与差异化、MVP 范围、`knowledge` 模块、治疗模态初筛规则、更完整的知识来源及许可与保密注意事项、验证计划回流为 OmicsWeave 实验 | 6 |
| OmicsWeave 包结构与 v2 合并（增加 `adapters`、`visualization`、`cli`） | 7.2 |
| 路线图改为两个系统统一编号，并标出阶段之间的前置条件和工期估算 | 8 |
| 风险表补充：范围蔓延、未发表数据外泄、单人开发瓶颈、契约漂移 | 9 |
| 新增决策记录与待决策事项（2026-10-29 已确认物种、运行环境、旧打分模块、课题工作区等）；新增 O5 集群与联网阶段 | 8、11 |

---

## 1. 执行摘要与总体决策

拟构建两个**相互兼容、独立维护**的科研软件系统：

- **OmicsWeave — Biological Evidence Generation Platform**：处理转录组、靶向与非靶向代谢组、代谢流分析（MFA）、蛋白质组及多组学整合数据，输出可复现的定量结果、标准化生物学实体及关联证据、候选通路和靶点假设，以及科研报告。
- **TheraWeave — Evidence-Guided Therapeutic Discovery Platform**：读取上游标准化证据，融合外部生物医学知识，评估候选靶点的机制合理性、可干预性和潜在风险，比较不同治疗模态，并生成供研究人员审核的候选方案和验证研究规划。

**核心架构决策**：

1. **双仓库**：OmicsWeave 与 TheraWeave 分别建仓；不把药物研发功能直接整合进入 OmicsWeave。
2. **结构化数据优先**：Parquet/JSON 等机器可读成果是正式交付接口；HTML/PDF 是同一分析结果的人类阅读视图。
3. **保留组学差异**：各组学保留其专有单位、模型假设、质量指标和不确定性，不能将所有结果硬性压平成通用的差异表达表。
4. **证据和假设分离**：真实测量、统计关联、机制推断、治疗假设分别存储；候选靶点并不等于已验证治疗靶点。
5. **可追溯**：每项候选结果均能回溯至样本、输入数据、分析版本、参数、数据库版本及相关证据。
6. **按阶段落地**：首先完成可靠的分析输出与证据交换规范，再实现靶点优先级，最后扩展各类治疗设计工具。
7. **职责边界**：OmicsWeave 只回答“数据中发生了什么变化、证据有多强”；成药性、已知药物、安全性和治疗模态的判断全部属于 TheraWeave。
8. **只通过文件契约耦合**：两个系统之间不互相导入代码。契约以 JSON Schema 文件发布，TheraWeave 不需要安装 OmicsWeave 就能读取和校验结果包，也能接收文献、GWAS 等其他来源的候选清单。
9. **排期门控**：TheraWeave 的编码在交换契约 v1 定稿（OmicsWeave 阶段 O2）之后启动；在此之前只做设计和契约评审。

### 1.1 端到端总体流程

```text
实验 / 公共数据
    │
    ▼
┌──────────────────────────────────────────────────────────────┐
│ OmicsWeave                                                    │
│ RNA-seq | Targeted/Untargeted Metabolomics | MFA | Proteomics │
│ Quality Control → Quantification/Statistics → Integration     │
│ → Biological Evidence → Pathways & Target Hypotheses          │
└──────────────────────────────────────────────────────────────┘
    │
    ├─ AnalysisResultBundle：矩阵、表格、证据、假设、来源信息
    └─ Research Report：HTML（默认）/ PDF（可选）
    │
    ▼
版本化 Evidence Exchange Contract（JSON Schema；OmicsWeave 内由 Pydantic 模型生成）
    │
    ▼
┌──────────────────────────────────────────────────────────────┐
│ TheraWeave                                                    │
│ Evidence Import → Knowledge Integration → Target Ranking      │
│ → Modality Selection → Therapeutic Candidate Assessment       │
│ → Validation Planning → Decision Report                      │
└──────────────────────────────────────────────────────────────┘
    │
    ▼
人工审查 → 独立数据验证 → 功能实验 → 研究证据反馈
    │
    └─ 验证实验以 OmicsWeave 实验模板（experiment.yaml）的形式回到上游分析，形成闭环
```

> **解释边界**：该系统用于研究假设生成、候选优先级评估与研发决策支持，不应将计算预测直接表述为已经确定的临床有效性、安全性或治疗建议。

### 1.2 当前状态（2026-10-29）

| 系统 | 状态 |
|---|---|
| OmicsWeave | 阶段 O0（工程化打包）已完成：可 `pip install`，统一为 `omicsweave.*` 包内导入，有统一命令行，代码与实验数据目录已解耦；185 个测试在 Python 3.10/3.12 下 182 通过、3 跳过（需真实数据）；已建立不提交的本地课题工作区 `local/`，可用 `pytest --project-dir` 在真实数据上测试；CI 已编写但尚未在 GitHub 上运行；阶段 0 已提交（7ef7cf4）。阶段 O1 进行中：包改名、项目模型与 `validate`、`init`、`migrate`、数据契约草案、示例项目已完成，课题专用脚本已移到本地工作区 |
| TheraWeave | 未启动；名称初步检索未发现重名（见第 10 节） |

---

## 2. OmicsWeave 的职责及输出层次

### 2.1 五层结果模型（Layered Result Model）

| 层级 | 结果类别 | 核心内容 | 推荐载体 |
| --- | --- | --- | --- |
| L1 | Experimental Data | 样本元数据、实验设计、处理后的数据矩阵、MID、质谱定量信息；原始下机文件只登记路径与校验值，不复制 | Parquet、JSON；H5AD/H5MU 为可选导出 |
| L2 | Analysis Results | DEG、差异代谢物、差异蛋白、代谢通量、富集统计和质量指标 | Parquet |
| L3 | Biological Evidence | 基因—转录本—蛋白质—代谢物—反应—通路及组织/细胞类型的标准化关联和证据 | Parquet、JSON |
| L4 | Research Hypotheses | 候选通路、靶点、可能机制、证据分级、知识缺口与验证状态 | JSON、Parquet |
| L5 | Research Reports | 数据质量、主要结果、图表、候选列表、方法学、局限性与溯源链接 | HTML、可选 PDF |

**不变原则**：L1/L2 的原始结果与 L3/L4 的知识解释彼此独立，报告由上述数据生成，不能替代机器可读结果。

### 2.2 各组学模块需要保留的关键字段

| 模块 | 典型结果 | 必备上下文或质量信息 |
| --- | --- | --- |
| Transcriptomics | Gene/Transcript ID、counts、TPM（适用时）、log2FC、p value、FDR | 物种、参考基因组、对比定义、归一化方法、批次/混杂校正 |
| Targeted Metabolomics | 标准品对应代谢物、峰面积/绝对定量浓度、内标结果 | 单位、标准曲线、QC、定量上下限、加合物与峰归属 |
| Untargeted Metabolomics | feature ID、m/z、RT、峰面积、统计差异、候选注释 | feature 与 compound 分离、鉴定等级（MSI 1–4 级）、批次效应及缺失值处理 |
| MFA | Reaction ID、净/交换通量、flux estimate、95% CI 或可行范围、模型拟合 | 通量单位、代谢网络及原子映射、MID 观测/预测、残差、稳态假设、可辨识性、是否碰到拟合边界 |
| Proteomics | UniProt ID、protein group、peptide evidence、abundance、log2FC、FDR | 蛋白推断歧义、定量方式、肽段支持、批次信息 |
| Multi-Omics | 跨组学实体映射、关联网络、通路/模块结果 | 样本配对、跨组学对比条件、模型方法、多个假设检验及证据冲突 |

**所有组学的统计结果都必须带上**：实际参与检验的**生物学重复数**，以及秩检验在该样本量下可达到的最小 p 值。例如 n=3 对 3 时，Mann-Whitney 检验的 p 值不可能小于 0.05。

**特别注意**：

- MFA 的**通量**与代谢组学的**丰度/浓度**是不同物理量；联合解释时需要显式区分。
- MFA 的置信区间、拟合残差和结构/参数可辨识性不能省略。
- MFA 的推断单位是**动物**而不是拟合迭代：约 150 次拟合迭代描述的是同一只动物的参数不确定度，不能当作独立重复（现有代码已遵循这一点，迁移时必须保留）。
- 非靶向代谢组必须区分**未注释 feature**、**推定化合物**和**标准品确认化合物**，不能把不确定的峰直接映射为可靠药物靶点证据。
- 基因、转录本和蛋白同工型之间不保证一一对应。

### 2.3 实验设计语义（新增）

项目、实验、样本、因素与比较的完整定义见 v2 文档第 5 节。这里补充一项对治疗假设至关重要的设计信息：**每个比较的角色（`role`）**。

| `role` 取值 | 含义 | 示例 |
|---|---|---|
| `disease_effect` | 疾病或造模引起的变化 | HFD vs LFD |
| `treatment_effect` | 干预引起的变化 | 在 HFD 小鼠中 GLP1 vs Saline |
| `time_effect` | 时间引起的变化 | 14 周 vs 4 周 |
| `other` | 其他比较 | — |

当同一实体同时出现在一个 `disease_effect` 比较和一个 `treatment_effect` 比较中、且方向相反时，标记为**治疗逆转（reversal）**。例如：某分子在 HFD 相对 LFD 时上调，而 GLP1 处理使其下调。
- 这类“疾病引起、治疗逆转”的变化，是最有价值的治疗假设来源。
- 逆转只说明“与治疗效果相关”，不说明“是治疗起效的原因”。

`experiment.yaml` 中的写法示例：

```yaml
contrasts:
  - {id: HFD_vs_LFD, factor: diet, numerator: HFD, denominator: LFD, role: disease_effect}
  - {id: GLP1_vs_Saline_in_HFD, factor: treatment, numerator: GLP1, denominator: Saline,
     subset: {diet: HFD}, role: treatment_effect}
reversal_pairs:
  - {disease: HFD_vs_LFD, treatment: GLP1_vs_Saline_in_HFD}
```

---

## 3. 跨组学实体标准化与证据图谱

### 3.1 BiologicalEntity 标准模型

建议使用统一的内部实体记录，但为不同类型保留独立命名空间与多对多映射：

| 实体类型 | 主要外部标识符（示例） | Biolink 类别 | 主要用途 |
| --- | --- | --- | --- |
| Gene | Ensembl Gene ID、HGNC/MGI | `biolink:Gene` | 基因定位、差异表达、遗传靶点 |
| Transcript / mRNA | Ensembl Transcript ID、RefSeq | `biolink:Transcript` | 转录本层级和 RNA 相关干预 |
| Protein | UniProt accession | `biolink:Protein` | 蛋白、酶及受体靶点 |
| Metabolite | ChEBI、HMDB、InChIKey | `biolink:SmallMolecule` | 代谢物与化学结构映射 |
| Reaction（新增） | Rhea、KEGG Reaction、MFA 模型内反应 ID | `biolink:MolecularActivity` | MFA 通量与酶/基因的连接 |
| Pathway | Reactome stable ID、GO term | `biolink:Pathway` | 机制与通路富集解释 |
| Tissue / Organ | UBERON | `biolink:AnatomicalEntity` | 解剖组织定位 |
| Cell Type | Cell Ontology（CL） | `biolink:Cell` | 细胞类型定位 |
| Disease | MONDO、EFO、DOID | `biolink:Disease` | 疾病/表型关联 |

所有外部注释应记录 `namespace`、`source_version`、`organism_taxonomy_id`、映射方法及置信度（适用时）。生物学实体的关系可使用 `entity_relationships.parquet` 管理。
- 标注 Biolink 类别是为了让结果能直接接入生物医学知识图谱。类别名称以所用的 Biolink Model 版本为准，记录在 manifest 中。
- MFA 的反应 ID 来自所用的代谢网络模型，必须与模型版本一起记录；它到 Rhea、KEGG 和催化酶基因的映射，由网络模型定义文件提供，不再硬编码在代码里（对应 v2 第 10.3 节）。

### 3.2 区分组织和细胞类型证据的语义

以下字段**不能互相代替**：

- `sample_tissue`：本研究样本实际来自何种组织。
- `expression_context`：特定分子在外部数据/组织中表达的观察证据。
- `proposed_target_tissue`：依据多种信息提出的潜在干预组织，属于假设。
- `cell_type_id`：在存在细胞分辨率证据时附加；Bulk RNA-seq 通常不能直接证明细胞型特异性。若来自反卷积或外部单细胞图谱，必须同时记录 `cell_type_source`。

### 3.3 Evidence 与 TargetHypothesis 分开存储

推荐数据对象：

- `BiologicalEntity`：标准实体及其标识。
- `Association`：两个实体或实体与通路/疾病之间的关联。
- `Evidence`：单条证据（实验/算法/文献/数据库、效应量、统计数据、方向、限制）。
- `PathwayCandidate`：候选通路及命中实体、富集方法、背景集、校正后的统计结果。
- `TargetHypothesis`：可接受进一步研究的干预假设，其证据通过 `evidence_refs` 关联。
- `AnalysisRun`：分析操作的标识、运行参数、软件和数据版本信息。

证据建议至少分为三个层次：

1. **Observed Findings**：直接测量或模型估计支持的结果；记录其统计和模型条件。
2. **Inferred Mechanisms**：跨组学或数据库推断出的机制，明确标记推断方法和不确定性。
3. **Therapeutic Hypotheses**：关于干预靶点与治疗策略的候选假说，必须标记尚缺少的功能、因果及安全证据。

**受控词表（新增）**：以下字段只能取表中的值，便于下游程序可靠地解析。

| 字段 | 取值 |
|---|---|
| `evidence_level` | `observed` / `inferred` / `hypothesis` |
| `evidence_status` | `computational_hypothesis` / `literature_supported` / `experimentally_supported` / `refuted` |
| `validation_status` | `not_tested` / `planned` / `in_progress` / `supported` / `not_supported` / `inconclusive` |
| `effect_direction` | `up` / `down` / `mixed` / `none` / `not_applicable` |
| `intervention_direction` | `inhibit` / `activate` / `replace` / `edit` / `unknown`（由 TheraWeave 填写；OmicsWeave 一律输出 `unknown`） |
| `evidence_grade` | `high` / `medium` / `low` / `exploratory`（OmicsWeave 按证据画像规则给出，规则及其版本写入 manifest） |
| 缺失状态 | 字段为 `null` = 未提供；`not_analyzed` = 没有做这项分析；`negative` = 做了但没有发现。三者不可混用 |

### 3.4 TargetHypothesis 示例（字段为示意）

```json
{
  "schema_version": "1.0.0",
  "study_id": "study_example",
  "analysis_id": "run_001",
  "species": { "taxonomy_id": 10090 },
  "target_hypotheses": [
    {
      "hypothesis_id": "TH-001",
      "target": {
        "entity_type": "gene",
        "entity_id": "GENE_ID_PLACEHOLDER",
        "id_namespace": "Ensembl",
        "biolink_category": "biolink:Gene"
      },
      "human_ortholog": {
        "entity_id": "HUMAN_GENE_ID_PLACEHOLDER",
        "id_namespace": "Ensembl",
        "orthology_type": "one2one",
        "source": "SOURCE_AND_VERSION_PLACEHOLDER"
      },
      "associated_pathways": ["PATHWAY_ID_PLACEHOLDER"],
      "tissue_context": {
        "sample_tissue": "UBERON_ID_PLACEHOLDER",
        "expression_context": null,
        "proposed_target_tissue": null,
        "cell_type_id": null,
        "cell_type_source": null
      },
      "contrasts": [
        { "contrast_id": "HFD_vs_LFD", "role": "disease_effect", "effect_direction": "up" },
        { "contrast_id": "GLP1_vs_Saline_in_HFD", "role": "treatment_effect", "effect_direction": "down" }
      ],
      "reversal": true,
      "evidence_refs": ["EVIDENCE-001", "EVIDENCE-002"],
      "evidence_level": "hypothesis",
      "evidence_status": "computational_hypothesis",
      "evidence_grade": "exploratory",
      "validation_status": "not_tested",
      "intervention_direction": "unknown",
      "missing_validation": [
        "causal_validation",
        "protein_level_confirmation",
        "tissue_specificity"
      ]
    }
  ],
  "provenance": {
    "software": "OmicsWeave",
    "software_version": "0.1.0",
    "run_manifest": "manifest.json"
  }
}
```

`evidence_refs` 必须实际指向数据包中的独立证据记录。示例 ID 不对应任何真实生物学靶点。

### 3.5 候选通路与靶点表建议字段

**`pathway_candidates.parquet`**：`pathway_id`、`pathway_name`、`organism`、`method`、`background_set`、`contrast_id`、`contrast_role`、`effect_direction`（若方法支持）、`p_value`、`fdr`、`hit_entity_ids`、`supporting_omics`、`n_tissues_concordant`、`evidence_grade`、`evidence_refs`。

**`target_candidates.parquet`**：

| 字段组 | 字段 |
|---|---|
| 标识 | `target_id`、`target_type`、`entity_id`、`id_namespace`、`biolink_category`、`organism` |
| 人类映射 | `human_ortholog_id`、`orthology_type`（如 one2one / one2many）、`ortholog_source` |
| 上下文 | `pathway_ids`、`sample_tissue`、`cell_type_id`、`cell_type_source` |
| 证据 | `contrast_id`、`contrast_role`、`effect_direction`、`reversal`、`supporting_omics`、各组学的效应量与校正后 p 值、`n_bio_replicates`、`n_omics_concordant`、`n_tissues_concordant`、`identification_level`（代谢物） |
| 结论 | `evidence_grade`、`causal_evidence_status`、`validation_status`、`evidence_refs`、`limitations`、`flags`（如 `low_n`、`at_bound`、`annotation_uncertain`） |

> 不建议单凭差异表达的 FDR、组学变化方向或通路富集结果生成无条件的“治疗靶点排名”。统计相关性、机制因果性、可成药性和安全性是不同概念。OmicsWeave 只输出**证据等级与各维度指标**，不输出“治疗优先级分数”。

### 3.6 小鼠→人同源映射（新增）

药物最终作用于人类蛋白，因此由小鼠数据得到的候选必须附带人类同源基因：
- 映射来源（如 Ensembl Compara、HGNC HCOP、MGI）及其版本写入 manifest。
- 一对多、多对多的情况必须保留，不能只取第一个。
- 没有人类同源基因的候选不删除，而是标记为 `no_human_ortholog`，交由下游判断。

---

## 4. 推荐的分析结果包（AnalysisResultBundle）

每次完整分析应有**不可混淆的运行 ID**，并生成自描述结果目录。路径与 v2 的项目目录约定统一：

```text
<项目目录>/analysis_results/<experiment_id>/<run_id>/
├── manifest.json
├── metadata/
│   ├── subjects.parquet              # 由 subjects.tsv 解析校验后写出
│   ├── samples.parquet
│   └── experiment.resolved.json      # 展开后的实际配置（含比较及其角色）
├── omics/
│   ├── transcriptomics/
│   │   ├── data/                     # L1：matrix / features / dataset.json
│   │   ├── qc/
│   │   ├── tables/                   # L2：differential / enrichment
│   │   └── figures/
│   ├── metabolomics/
│   ├── mfa/
│   └── proteomics/
├── integration/                      # 跨组学分析（flux_anchored / pathway_convergence / sample_level）
├── evidence/
│   ├── biological_entities.parquet
│   ├── evidence.parquet
│   ├── entity_relationships.parquet
│   ├── orthologs.parquet
│   └── pathways.parquet
├── hypotheses/
│   ├── pathway_candidates.parquet
│   ├── target_candidates.parquet
│   └── target_hypotheses.json
├── reports/
│   ├── analysis_report.html
│   └── figures/
├── qc/
│   └── quality_metrics.json
└── provenance/
    ├── software_versions.json
    ├── parameters.json
    └── ro-crate-metadata.json
```

单组学运行时只生成对应的 `omics/<组学>/` 子目录和 `manifest.json`；`integration/`、`evidence/`、`hypotheses/` 只在启用整合时生成。

### 4.1 文件格式选择

- **Parquet（正式格式）**：用于结构化定量表、分析结果和证据清单；列类型可靠、支持列表型字段（如 `hit_entity_ids`），适合批处理与下游程序读取。`pyarrow` 列为基础依赖。
- **TSV / XLSX（导出格式）**：供人工查看，通过 `omicsweave export --format tsv|xlsx` 生成，不作为交换契约的一部分。原因是 Excel 会自动改写基因名、日期等，且列类型不可靠。
- **H5AD/H5MU**：适用于 AnnData/MuData 数据模型能准确表示的矩阵和多模态数据；作为可选导出，不要求所有 MFA/质谱数据转换为这些格式。
- **JSON**：用于 Schema、关联对象、假设、实验设计、参数与 manifest。
- **HTML**：默认报告形式，可实现可交互图表、结果筛选、证据点击跳转。
- **PDF**：论文补充材料、归档、会议汇报等场景按需生成。

> v2 文档已于 v2.1 同步本节决定（第 6.2、7、8 节）。

### 4.2 可复现性与来源追踪

`manifest.json` 建议包含：

- `bundle_schema_version`、`study_id`、`analysis_id`、`created_at`。
- 模块列表、输入文件摘要（checksum）、输出文件清单及其内容类型。
- 物种、实验条件、样本与比较组定义（含比较角色和逆转配对）。
- 代码提交版本（Git SHA）、Python/依赖版本、算法参数与随机种子（如适用）。
- 数据库/注释版本、参考基因组/代谢网络版本、同源映射来源版本、Biolink Model 版本。
- 外部知识库的检索日期（如使用在线接口）。
- QC 状态、重要告警、不可解析实体映射及数据缺失信息。

可借鉴 **RO-Crate** 描述结果包结构，借鉴 **W3C PROV-O** 记录数据、分析活动与来源的关系。第一版可以使用简单 manifest，逐步增加兼容性。

### 4.3 报告章节建议

1. 研究概览与实验设计。
2. 数据质量、批次效应与处理流程。
3. 各组学主要结果（附关键统计指标）。
4. 跨组学整合：一致性、互补证据和冲突证据。
5. 候选代谢/信号通路及机制假说。
6. 候选基因、转录本、蛋白受体/酶及组织细胞上下文（含人类同源基因、是否被治疗逆转）。
7. 证据可靠度、模型假设、局限性与待验证问题。
8. 可复现信息、参数、软件/数据库版本及结果文件索引。
9. 自动生成的方法学段落（软件版本、参数、统计方法），可直接改写进论文。

**报告与数据包必须使用同一数据源生成**，避免图表与导出表格出现版本不一致。报告只陈述统计事实与方法参数，候选清单单独成节并明确标注“假设，需实验验证”，不写入课题专用的生物学结论。

---

## 5. OmicsWeave 与 TheraWeave 的交互接口

### 5.1 Evidence Exchange Contract

建议采用**版本化文件交换契约**，起步时不要求部署服务器或数据库：

- 以 `AnalysisResultBundle` 为交换单位。
- 使用 JSON Schema 或 Pydantic 定义验证规则。
- 严格检查文件存在、字段类型、物种一致性、标识符命名空间、单位、证据引用完整性。
- 对不可解析的数据给出明确告警；禁止静默丢弃关键字段。
- 区分字段缺失（`null`）、证据为阴性、尚未进行分析三种状态（取值见 3.3 节）。
- 按 Schema Version 管理兼容范围和破坏性变更。

**契约的存放与版本规则（新增）**：

| 项目 | 约定 |
|---|---|
| 唯一来源 | OmicsWeave 中的 Pydantic 模型（`omicsweave.contract`） |
| 发布形式 | 由模型生成 JSON Schema，提交到 OmicsWeave 仓库的 `schemas/evidence-exchange/v<主版本>/` |
| 版本号 | 语义化版本：新增可选字段 → 次版本号 +1；删除或改名字段、改变取值含义 → 主版本号 +1 |
| 下游声明 | TheraWeave 在配置中声明支持的契约版本范围（如 `>=1.0,<2.0`），超出范围时拒绝导入并给出原因 |
| 契约测试 | 一组小型“标准样例结果包”（合成数据）随 Schema 一起发布；两个仓库的 CI 都对它运行读取与校验测试 |
| 单位 | 尽量使用 UCUM 单位代码，便于程序比较和换算 |

**建议的数据流**：

```text
OmicsWeave Result Bundle
   ├─ omics measurements & statistical estimates
   ├─ biological evidence
   ├─ hypotheses
   └─ provenance / QC
               │
               ▼
        Schema Validation
               │
               ▼
     TheraWeave Evidence Import
               │
               ▼
     External Knowledge Integration
               │
               ▼
        Target & Modality Review
```

### 5.2 未来 Python API 示例（设计目标，非现有可运行代码）

```python
# 上游：在 OmicsWeave 中读取自己的结果包
from omicsweave import load_result_bundle
bundle = load_result_bundle("analysis_results/20250613_dio_glp1/run_2026_001/")

# 下游：TheraWeave 只依赖 JSON Schema，不需要安装 OmicsWeave
from theraweave import DiscoveryProject
project = DiscoveryProject.from_bundle(
    "analysis_results/20250613_dio_glp1/run_2026_001/",
    validate_schema=True,
)
assessment = project.assess_targets()          # 分维度评估，不是单一分数
assessment.to_parquet("target_assessment.parquet")
```

对外同时支持纯文件接口（命令行 `theraweave import <结果包目录>`），以防下游运行环境无法安装 Python 包。

### 5.3 反向反馈而非篡改上游结论

TheraWeave 生成的人工审查、功能验证、靶点否证、药物筛选结果等应作为**新的证据记录**，通过来源引用关联原始假设，而不是覆盖已有的 OmicsWeave 结果。

**验证闭环（新增）**：TheraWeave 输出的验证计划包含一份 OmicsWeave 实验模板（`experiment.yaml` + 样本表模板，其中预填好分组、处理、时间点和比较角色）。验证实验完成后，数据按同一项目模型进入 OmicsWeave 分析，结果包再回到 TheraWeave 更新 `validation_status`。

---

## 6. TheraWeave：独立药物发现与治疗策略研究系统

### 6.1 目标与范围

**目标**：将可信的组学证据与外部知识结合，形成候选靶点优先级、干预策略比较、候选治疗方案的研究性评价，以及可审查的验证计划。

**不应作为初期目标**：独立重建所有化学、蛋白结构、CRISPR、RNA 药物设计或临床决策算法；对研究计算结果作出未经验证的有效性/安全性承诺。

**定位与差异化（新增）**：
- 每种治疗模态本身都是一个成熟领域，各有专业工具；只在这些工具外面包一层，附加价值很低。
- TheraWeave 的价值在于**“组学证据 → 干预方案 → 验证实验”这一转化环节**。
- 首个聚焦领域：**代谢性疾病**，以及由 MFA 通量锚定的代谢酶、转运体靶点。这与上游数据（DIO、GLP-1、袖状胃切除等研究）直接相关，其中 GLP-1 本身就是多肽药物。

**MVP 范围（新增）**：
- **包含**：结果包导入与校验、外部知识整合、分维度靶点评估、治疗模态初筛与比较、验证计划（含 OmicsWeave 实验模板）、决策报告。
- **不包含**：分子生成、结构设计、任何治疗模态的设计工具对接（放到阶段 T2）。

### 6.2 功能模块

| 模块 | 工作内容 | 典型输出 |
| --- | --- | --- |
| `contract` | 读取与校验上游结果包（基于 JSON Schema），也接收其他来源的候选清单 | 校验报告、内部证据对象 |
| `knowledge`（新增） | 外部知识源注册表：每个来源的版本、下载/检索日期、许可；本地快照与缓存 | 可复现的知识快照 |
| `evidence` | 合并上游证据与外部证据 | 标准化 Evidence Graph |
| `targets` | 优先级、多维证据、可干预性、因果/安全缺口 | Target Assessment |
| `modalities` | 比较不同治疗模态的机制适配与限制 | Modality Comparison |
| `repurposing` | 匹配现有药物与可引用的靶点证据 | Repurposing Candidates |
| `evaluation` | 采用适当基准集和指标验证模型 | Model Evaluation |
| `experiment_planning` | 研究假设、关键对照、终点、样本量估算、风险与验证优先级；导出 OmicsWeave 实验模板 | Reviewable Validation Plan |
| `reporting` | 导出可审查的研究决策报告 | HTML/PDF/JSON |

### 6.3 预留的治疗模态

| 模态 | 适合的研究功能 | 不可忽略的条件 |
| --- | --- | --- |
| Small molecule | 药物再利用、虚拟筛选、性质预测、候选分子评价 | 靶点结合与选择性、ADMET、实验确认 |
| DNA / Gene editing | 目标基因和编辑策略适配性、现有工具对接 | 脱靶、递送、组织特异性、持久性效应 |
| mRNA / RNA therapeutics | mRNA 补充、siRNA/ASO 等机制差异化比较 | 干预方向、RNA 稳定性、递送及免疫反应 |
| Peptide therapeutics | 多肽候选方案的结构与性质评价 | 稳定性、特异性、亲和力和递送 |

**重要区分**：RNA 表达补充与 RNA 降解抑制属于不同的干预机制；同样，基因差异表达方向不能直接决定应该激活、抑制、替代或编辑某个靶点。

**模态初筛规则（新增，仅作初筛，必须由专家审核）**：

| 靶点特征 | 优先考虑的模态 | 主要依据 |
|---|---|---|
| 膜受体或胞外蛋白 | 抗体、多肽 | 亚细胞定位（UniProt）、已有配体 |
| 具有可结合口袋的酶 | 小分子 | 已知配体与可成药性评估（ChEMBL、Open Targets） |
| 在肝脏表达、需要降低表达 | siRNA / ASO（肝靶向递送较成熟） | 组织表达（HPA、GTEx） |
| 功能缺失、需要恢复表达 | mRNA、基因治疗 | 干预方向、递送可行性 |
| 已有获批药物作用于该靶点或同一通路 | 药物再利用 | 已知药物（ChEMBL、DrugBank） |

### 6.4 靶点优先级应采用多维证据，而非不透明单一得分

建议至少评价：

1. **统计证据**：效应量、显著性、研究质量与重复性。
2. **多组学一致性**：跨组学是否互相支持，是否存在合理冲突。
3. **机制与因果性**：独立干预实验、遗传学证据等支持程度；是否被治疗逆转。
4. **组织/细胞上下文**：靶部位可及性、细胞类型和表达范围。
5. **可干预性**：分子类型、结构或功能特征、可用治疗模态。
6. **安全性**：潜在组织外效应、基本生物学功能与已有风险信息（如小鼠敲除表型）。
7. **转化相关性**：物种差异、人体疾病相关性和外部重复证据。

早期优先展示**分维度评分、来源、缺失值与证据冲突**；没有充分验证前，不应对外宣称“模型分数=成药成功概率”。

### 6.5 可集成的现有工具与知识来源（选型示例）

| 用途 | 来源 |
|---|---|
| 靶点—疾病证据、按模态的可成药性评估 | **Open Targets** |
| 已知配体、生物活性、作用机制 | **ChEMBL** |
| 已获批药物及其靶点 | **DrugBank**（使用条款按用途区分，商业用途需授权） |
| 蛋白功能、亚细胞定位、蛋白类别 | **UniProt** |
| 组织与细胞表达（安全性评估） | **Human Protein Atlas**、**GTEx** |
| 小鼠敲除表型 | **IMPC**、**MGI** |
| 小鼠→人同源映射 | **Ensembl Compara**、**HGNC HCOP** |
| 通路与生物学实体映射 | **Reactome** |
| 化学信息学基础能力 | **RDKit** |
| 分子机器学习基础组件 | **DeepChem** |
| 药物发现数据任务与基准评价 | **Therapeutics Data Commons** |
| 基因编辑相关靶位点分析 | **CHOPCHOP** |

这些均应当通过适配器与版本记录接入；**不默认提供外部平台授权、数据许可或网络访问权限**。

**许可与保密（新增）**：
- 各数据源的许可条款不同，部分限制商业用途或再分发。`knowledge` 模块须为每个来源记录许可，并在报告中注明。
- **用在线接口查询候选清单，等于把未发表的研究假设发送给第三方。** 默认应使用下载到本地的数据快照；需要在线查询时，须由用户显式开启，且只发送标识符，不发送任何组学数据。
- 使用本地快照也有利于结果复现。

### 6.6 建议的仓库结构

```text
TheraWeave/
├── pyproject.toml
├── README.md
├── src/
│   └── theraweave/
│       ├── core/
│       ├── contract/          # 结果包导入与 JSON Schema 校验
│       ├── knowledge/         # 知识源注册表、本地快照、各来源适配器
│       ├── evidence/
│       ├── targets/
│       │   ├── prioritization/
│       │   └── tractability/
│       ├── modalities/
│       │   ├── small_molecule/
│       │   ├── gene_editing/
│       │   ├── rna_therapeutics/
│       │   └── peptide/
│       ├── repurposing/
│       ├── experiment_planning/
│       ├── evaluation/
│       ├── reporting/
│       └── cli.py
├── tests/
│   └── contract/              # 使用 OmicsWeave 发布的标准样例结果包
├── examples/
└── docs/
```

依赖策略：基础安装只包含导入、知识整合和报告所需的轻量依赖；RDKit、深度学习框架等按模态放入可选依赖（extras），与 OmicsWeave 的做法一致。

---

## 7. 仓库边界、依赖策略与治理

### 7.1 为什么分成两个独立仓库

| 维度 | OmicsWeave | TheraWeave |
| --- | --- | --- |
| 目标 | 组学定量分析与生物学证据 | 治疗靶点及干预策略决策支持 |
| 核心输入 | 实验数据与组学矩阵 | 组学证据包、知识库与研究证据 |
| 关键输出 | 可复现的分析结果和研究假设 | 候选靶点、模态比较、验证计划 |
| 计算生态 | 统计分析、MFA、质谱和组学专用包 | 化学信息学、结构/分子模型、治疗模态工具 |
| 验证重点 | 算法正确性、统计与模型假设 | 基准评价、独立验证、转化风险 |
| 适用人群 | 组学分析与系统生物学研究者 | 靶点研究和治疗发现研究者 |

独立仓库可以降低复杂依赖的相互干扰、允许不同发布周期，也可以让 TheraWeave 将来接收**非 OmicsWeave** 来源的数据。

### 7.2 OmicsWeave 内部仓库策略

延续此前决策：**OmicsWeave 初期采用 Monorepo + 单 Python 分发包**，各组学使用独立模块和可选依赖（extras）；成熟后按需演化为同仓库多包，不因功能种类多就立即拆分仓库。

目标包结构（合并 v2 第 4.2 节）：

```text
OmicsWeave/
├── pyproject.toml
├── schemas/evidence-exchange/v1/   # 发布的 JSON Schema 与标准样例结果包
├── src/omicsweave/
│   ├── core/                       # 项目模型、配置校验、统计、ID 映射、溯源
│   ├── contract/                   # 结果包与交换契约的 Pydantic 模型、读写与校验
│   ├── adapters/                   # 外部工具适配器：原始数据 → 分析层
│   ├── transcriptomics/
│   ├── metabolomics/
│   ├── mfa/
│   ├── proteomics/
│   ├── integration/                # 跨组学统计整合（由现 multiomics/ 改名）
│   ├── evidence/                   # 实体标准化、证据记录、候选假设生成
│   ├── reporting/
│   ├── visualization/
│   └── cli.py
├── tests/
├── examples/
└── docs/
```

**依赖方向**：`core` 不依赖具体组学模块；各组学模块尽量仅依赖 `core`；`integration/evidence/reporting` 通过稳定接口组合模块结果。`integration` 负责跨组学的统计分析，`evidence` 负责把统计结果转成标准化的证据和假设记录，两者不合并。

### 7.3 是否建立独立公共 Schema 仓库

**初期不建议**再建第三个仓库。先为交换格式制定版本化的文件规范（存放位置见 5.1 节），并在两个项目中进行兼容性测试。当交换契约稳定、已有外部用户或第三方项目需要直接依赖 Schema 时，再考虑单独抽取 `omics-evidence-schema` 轻量包。

### 7.4 测试与质量门禁

- 单元测试：统计计算、实体映射、序列化/反序列化和单位转换。
- Schema 契约测试：版本兼容、必填字段、引用完整性及异常输入；两个仓库的 CI 都运行标准样例结果包。
- 数据集回归测试：固定公开/合成测试集上的结果稳定性；迁移现有课题时，新旧结果逐项比对。
- MFA 专项测试：通量单位、可辨识性、MID 拟合、约束与网络版本。
- 报告一致性测试：报告内数值与机器可读结果一致。
- 人工科学审查：候选解释的适用范围和因果性是否被夸大。

---

## 8. 分阶段开发路线（Roadmap）

两个系统统一编号：`O` 为 OmicsWeave，`T` 为 TheraWeave。工期按 1 名主要开发者估算，仅供排期参考。OmicsWeave 各阶段的详细任务见 v2 文档第 11 节。

### 8.1 总览与前置条件

```text
O0 打包 ✅ ─▶ O1 核心契约与项目模型 ─▶ O2 各组学接入 + 证据/假设 + 契约 v1 ─┬─▶ O3 原始数据适配器 ─▶ O4 OmicsWeave 发表
                    │                                                  │
                    └─ T0 TheraWeave 设计与契约评审（不写功能代码）      └─▶ T1 TheraWeave MVP ─▶ T2 模态扩展 ─▶ T3 评估与发表
```

| 阶段 | 前置条件 | 工期 | 版本 |
|---|---|---|---|
| O0 工程化打包 | — | 已完成 | — |
| O1 核心契约与项目模型 | O0 | 6–8 周 | OmicsWeave v0.2 |
| T0 TheraWeave 设计 | 可与 O1 并行（只写文档） | 1–2 周 | 设计文档 |
| O2 各组学接入 + 证据层 + 契约 v1 | O1 | 6–8 周 | OmicsWeave v0.3、契约 v1.0 |
| T1 TheraWeave MVP | 契约 v1.0 定稿 | 8–10 周 | TheraWeave v0.1 |
| O3 原始数据适配器 | O2 | 8–12 周 | OmicsWeave v0.4 |
| T2 模态扩展 | T1 | 每种模态 4–8 周 | TheraWeave v0.2+ |
| O4 / T3 发表 | O2（O3 可选）/ T2 | 各 4–6 周 | v1.0 |
| O5 集群部署与多客户端联网 | O4；按实际需求启动 | 视规模而定 | v1.x |

> **单人开发约束**：同一时间只推进一个“写代码”的阶段。O2 之后先做 T1 还是 O3，取决于实验室的近期目标（见第 11 节）。

### Phase O0 — OmicsWeave：工程化打包 ✅

- [x] 初始化统一仓库及模块化 Python 包（`pyproject.toml`、可选依赖、`src/` 布局）。
- [x] 统一包内导入，删除 `sys.path` 修改；代码与实验数据目录解耦。
- [x] 统一命令行入口；测试迁至 `tests/`；编写 CI 配置。
- [x] 提交改动（7ef7cf4）。
- [ ] 在 GitHub 上跑通 CI。

### Phase O1 — OmicsWeave：核心契约与项目模型

- [x] 包改名：`common`→`core`、`multiomics`→`integration`；删除兼容转发模块。
- [x] 项目模型：`project.yaml` / `experiment.yaml` / `subjects.tsv` / `samples.tsv`，含比较角色与逆转配对。
- [x] 定义 `AnalysisRun` 和 `AnalysisResultBundle` 的最小字段（`omicsweave.contract`）：数据集、差异结果表（含 MFA 专有必填列）、富集结果表与运行清单（草案 0.1）。
- [x] MFA 求解结果标准格式（Excel 模板 / Parquet 目录）与旧版自研软件格式转换器（`omicsweave mfa template / check / convert`，v2 第 6.2 节）。
- [ ] 标准化单位、对比组、物种、软件版本与 QC 记录。
- [x] 合并溯源记录（`core/provenance.py`）与绘图主题（`visualization/theme.py`）。
- [x] `omicsweave init / validate / migrate`；`migrate` 调用各组学原有的加载器把旧配置转为统一模型（v2 第 5.7 节），已用合成的四组学旧配置测试，尚未在真实课题配置上运行。
- [ ] 删除废弃脚本（v2 第 10.5 节，已失效的已删除，其余按替代进度删除）；课题专用脚本移到 `local/projects/<课题>/analysis/`（v2 第 10.4 节）：六组织通路图与两个转录组调试脚本已移出，`run_cross_tissue_integration` 待 O2 的整合命令就绪后移出。
- [ ] 运行函数改为接受显式的项目上下文，逐步去除进程级全局状态（为 O5 预留）：项目目录已改为上下文变量；旧流程的模块级状态随 O2 各组学改写去除（清单见 v2 第 4.1 节）。
- [x] 合成示例项目（`examples/demo_project`，CI 中运行 validate）；补 CITATION.cff；修复 SSL 校验。

**验收标准**：示例项目能通过 `validate`；至少一个现有课题能用 `migrate` 转换。迁移后 MFA 结果与迁移前逐项一致的比对，要等 `mfa run` 读取统一模型，移到 O2 验收。

### Phase O2 — OmicsWeave：各组学接入、证据整合与报告

- [ ] RNA-seq、代谢组、MFA、蛋白组改为读取分析层数据、输出统一结果表（Parquet）：四个组学（转录组、蛋白组为表格层输入）均已写出差异结果表与富集结果表，与旧流程逐值一致（转录组 ORA 有意改为本地计算）；`omicsweave export` 已完成。
- [x] MFA 多示踪剂整合（同一模型、多种示踪剂）：`input.files` 读取多组结果，按 samples.tsv 的 `protocol` 分层，可辨识性门控 + 精度加权固定效应合并 + 示踪剂内分层置换 + Q 检验，另出示踪剂覆盖与一致性 QC（`omicsweave/mfa/multi_tracer.py`）。同一示踪剂内多批次的整合尚未接入统一运行。
- [ ] 单组学独立运行与三种整合模式（v2 第 9 节）：`omicsweave run` 运行框架已完成并接入四个组学；`omicsweave integrate` 已完成，pathway_convergence、flux_anchored 已接入，sample_level 待做；旧打分模块的加和总分已删除，证据提取改为 flux-anchored 证据构建器。
- [ ] 引入 `BiologicalEntity` / `Association` / `Evidence` / `PathwayCandidate` / `TargetHypothesis`。
- [ ] 接入通路数据库和标准实体 ID 映射；ID 映射改为数据文件 + 解析器；加入小鼠→人同源映射。
- [ ] 实现候选通路与靶点清单，保留证据冲突、逆转标记与未验证标记。
- [ ] 生成 HTML 研究报告（含自动方法学段落）；根据需求增加 PDF。
- [ ] 定义 Evidence Exchange Contract v1，发布 JSON Schema 与标准样例结果包，编写 Schema 校验器。

**验收标准**：结果可被独立程序读取、复现其来源信息，并且没有丢失特有的科学语义；任一候选通路/靶点可回溯到具体统计结果、输入样本和来源数据库。

### Phase O3 — OmicsWeave：原始数据适配器

- [ ] 转录组 `fastq.gz` → cutadapt → kallisto/salmon。
- [ ] 蛋白组 Thermo `.raw` → DIA-NN；`.raw` → mzML。
- [ ] 代谢组 `.wiff`/`.raw` → msconvert → MS-DIAL/MZmine/Skyline 结果导入；同位素示踪 MID 的天然丰度校正。

**验收标准**：每种原始格式都有公开小样例，能从原始文件跑到分析层数据。

### Phase T0 — TheraWeave：设计与契约评审

- [ ] 编写 TheraWeave 设计文档（以本文第 6 节为基础）。
- [ ] 以下游使用者的角度评审交换契约草案，需求反馈给 O1/O2。
- [ ] 正式核查名称与商标。

### Phase T1 — TheraWeave：靶点发现 MVP

- [ ] 建立独立 `TheraWeave` 仓库，基础依赖保持轻量。
- [ ] 实现对 Evidence Bundle 的导入、验证、错误报告（只依赖 JSON Schema）。
- [ ] `knowledge` 模块：至少对接 Open Targets、ChEMBL、UniProt 三个来源的本地快照，记录版本、日期与许可。
- [ ] 实现透明的分维度靶点证据比较和缺口识别。
- [ ] 模态初筛与比较（规则 + 依据说明）。
- [ ] 验证计划：对照、终点、样本量估算；导出 OmicsWeave 实验模板。
- [ ] 输出候选靶点与研究性干预方向报告，支持人工审核。

**验收标准**：解释清楚候选结果的排序依据、冲突证据和不可推断事项；用标准样例结果包可完整跑通；导出的实验模板能通过 `omicsweave validate`。

### Phase T2 — 治疗模态扩展与验证闭环

- [ ] 优先加入已有药物检索、再利用和小分子性质评估。
- [ ] 逐步增加 RNA、基因编辑和多肽模态的独立适配器；顺序按实验室湿实验验证能力决定。
- [ ] 建立统一验证实验结果数据模型与反馈链路。
- [ ] 根据用户、协作及版本需求决定是否拆分子包或公共 Schema。

### Phase O4 / T3 — 评估与发表

- [ ] 利用独立基准数据/实验结果评估靶点筛选质量。
- [ ] 英文文档、教程、API 文档；PyPI 发布与 Zenodo 存档。
- [ ] OmicsWeave 与 TheraWeave 分别投稿（JOSS / Bioinformatics Application Note 等）。

**验收标准**：经过独立研究数据检查后，系统能稳定支持可审查的研究决策，而非仅给出不可验证的生成式结论。

### Phase O5 — 集群部署与多客户端联网（v1.0 之后，按需）

- [ ] 多节点集群：原始数据处理与大批量分析交给作业调度器（如 Slurm）或 Snakemake/Nextflow，项目目录放在共享存储上。
- [ ] 多客户端联网：服务端提供任务提交、排队与状态查询接口，加入用户认证与项目级权限；客户端只提交配置、查看结果包。
- [ ] TheraWeave 视需要以同样方式提供服务接口。

前提：O1–O2 已去除进程级全局状态，所有步骤都能以“配置 + 输入文件 → 结果包”的方式独立运行（v2 第 4.1 节）。

---

## 9. 关键风险与架构护栏

| 风险 | 可能后果 | 设计护栏 |
| --- | --- | --- |
| 把统计关联当作因果关系 | 靶点筛选方向偏差 | 证据等级分层，明确“仅相关”状态 |
| 非靶向代谢物错误注释 | 通路/靶点错误映射 | 保留鉴定等级与未注释 feature |
| 混淆 MFA 通量与浓度 | 错误解释生理机制 | 保留单位、模型前提与通量置信区间 |
| 伪重复（把拟合迭代当作生物学重复） | 显著性被严重高估 | 统计只在生物学重复间进行；结果表记录实际重复数 |
| 跨物种直接映射 | 转化研究偏差 | 记录物种和 ortholog 映射关系，保留一对多映射 |
| 混淆取样组织与治疗靶器官 | 错误定位干预部位 | 使用独立组织/细胞语义字段 |
| 大型不可解释靶点评分 | 结论无法审查 | 提供逐维度证据和规则 |
| 上下游紧耦合 | 系统升级困难 | 独立仓库、稳定的交换契约和契约测试 |
| 契约在两个仓库间漂移 | 下游静默读错数据 | 单一来源生成 Schema；语义化版本；双方 CI 运行同一组标准样例 |
| 报告脱离正式结果数据 | 图表与数据版本不一致 | 报告由同一结果包自动生成 |
| 外部知识库变更 | 不可复现与结果漂移 | 来源、版本、检索时间及许可记录；默认使用本地快照 |
| 未发表数据或假设外泄 | 竞争与合规风险 | 在线查询默认关闭，只发送标识符；示例只用合成或公开数据 |
| 药物设计预测被过度解读 | 误判安全性/有效性 | 仅定位科研决策支持，保留人工审核与验证环节 |
| 范围蔓延（同时铺开四种治疗模态） | 每种都做不深，维护负担过重 | MVP 不含设计工具；模态按湿实验能力逐个加入 |
| 单人开发瓶颈 | 两个仓库同时停滞 | 同一时间只推进一个写代码阶段；T1 在契约 v1 定稿后才启动 |

---

## 10. 项目名称与最终定位

| 系统 | 推荐名称 | 英文定位 | 仓库建议 |
| --- | --- | --- | --- |
| 多组学分析 | **OmicsWeave** | *Biological Evidence Generation Platform* | `OmicsWeave`：单仓库内模块化开发 |
| 下游治疗发现 | **TheraWeave**（暂定） | *Evidence-Guided Therapeutic Discovery Platform* | `TheraWeave`：独立 Git 仓库 |
| 公共数据契约 | （未来按需抽取） | *Evidence Exchange Schema* | 初期只保留版本化规范，不单独建仓 |

**名称检索结果（2026-10-29，初步）**：

| 名称 | PyPI | GitHub 仓库搜索 | 网页检索 |
|---|---|---|---|
| OmicsWeave | 无同名包 | 无同名仓库（不含本项目） | — |
| TheraWeave | 无同名包 | 无同名仓库 | 未发现同名产品；有名称相近的儿科康复平台 TheraWe / TheraWee，领域不同 |
| TheraFoundry（备选） | 无同名包 | 无同名仓库 | 未发现同名产品 |

**推荐 TheraWeave 而非 TheraFoundry**：
- 与 OmicsWeave 构成同一系列，“Weave”也贴合“编织证据与治疗策略”的定位。
- “Foundry”已被大量软件产品使用；“biofoundry”在合成生物学中专指自动化湿实验平台，容易引起误解。

以上只是初步检索，不等于商标核查。`TheraWeave`、`OmicsWeave` 的域名与商标可用性尚需正式核实；本文不表示这些名称已获得注册或独占使用资格。“Thera-” 前缀在生物技术公司名中较常见，商标核查时需特别留意。

### 一句话总结

> **OmicsWeave 负责把组学数据转化为可信、可追溯的生物学证据；TheraWeave 负责把这些证据与外部知识结合，提出可审查、可验证的治疗干预研究方向。**

---

## 11. 决策记录与待决策事项

2026-10-29 已确认的决定（OmicsWeave 侧的完整说明见 v2 文档第 13 节）：

| 事项 | 决定 | 影响阶段 |
|---|---|---|
| 契约模型与校验库 | Pydantic v2 作为唯一来源，导出 JSON Schema | O1 |
| 正式结果格式 | Parquet（本文 4.1 节；v2.1 已同步） | O1 |
| 首批支持的物种 | 小鼠（taxon 10090）+ 人（taxon 9606） | O2 |
| 运行环境 | 当前 Windows/Linux 单机；以后扩展多节点集群与多客户端联网（O5）。从 O1 起去除进程级全局状态，为此预留 | O1–O5 |
| 旧打分模块 `flux_omics_scorer` | 尚未发表且策略过于粗糙：删除加和总分与 `legacy_*` 输出，不保留 `--legacy`；多维证据画像作为唯一整合结果，按维度给出证据等级；证据提取部分保留并重构 | O2 |
| 课题专用脚本与数据 | 暂不移动外部项目的原始数据；本仓库已建立不提交的本地工作区 `local/projects/<课题>/`（脚本放 `analysis/`），通过 `--project-dir` / `OMICSWEAVE_PROJECT_DIR` / `pytest --project-dir` 指定 | O1 |
| 开源许可证 | BSD-3-Clause，OmicsWeave 与 TheraWeave 相同。OmicsWeave 版权方为 “The OmicsWeave Developers”，`LICENSE` 已加入；TheraWeave 建仓时按同样方式署名 “The TheraWeave Developers” | O0 / T1 |

仍待决定：

| 事项 | 建议 | 影响阶段 |
|---|---|---|
| O2 之后先做 T1 还是 O3 | 若近期目标是尽快形成治疗假设并开展验证，先做 T1；若近期要处理新的原始下机数据，先做 O3。现有课题数据已处理成表格，因此倾向于先做 T1 | O3 / T1 |
| 外部知识源及其许可 | 先用 Open Targets、ChEMBL、UniProt；DrugBank 等按用途核实条款后再接入 | T1 |
| TheraWeave 首批治疗模态 | 药物再利用 + 小分子性质评估；其余按湿实验验证能力排序 | T2 |

---

## 12. 可供实施时参考的规范与工具

以下链接用于查阅正式文档、决定后续技术选型；并非声明各项技术现已在项目中实现。

- [AnnData](https://anndata.readthedocs.io/) / [MuData](https://mudata.readthedocs.io/)：组学矩阵与多模态数据模型。
- [Apache Parquet](https://parquet.apache.org/)：列式结构化数据存储。
- [JSON Schema](https://json-schema.org/) / [Pydantic](https://docs.pydantic.dev/)：交换数据的校验规则。
- [RO-Crate](https://www.researchobject.org/ro-crate/) / [W3C PROV-O](https://www.w3.org/TR/prov-o/)：研究对象打包与来源追踪。
- [Biolink Model](https://biolink.github.io/biolink-model/)：生物医学实体与关系的类别体系。
- [UCUM](https://ucum.org/)：机器可读的计量单位代码。
- [Metabolomics Standards Initiative 报告标准（Sumner et al., 2007）](https://doi.org/10.1007/s11306-007-0082-2)：代谢物鉴定等级。
- [Reactome](https://reactome.org/) / [Open Targets Platform](https://platform.opentargets.org/)：通路与靶点关联证据。
- [ChEMBL](https://www.ebi.ac.uk/chembl/) / [DrugBank](https://go.drugbank.com/)：生物活性分子与已获批药物。
- [UniProt](https://www.uniprot.org/) / [Human Protein Atlas](https://www.proteinatlas.org/) / [GTEx](https://gtexportal.org/)：蛋白注释与组织表达。
- [IMPC](https://www.mousephenotype.org/) / [MGI](https://www.informatics.jax.org/)：小鼠表型。
- [Ensembl Compara](https://www.ensembl.org/info/genome/compara/index.html) / [HGNC HCOP](https://www.genenames.org/tools/hcop/)：同源基因映射。
- [RDKit](https://www.rdkit.org/) / [DeepChem](https://deepchem.io/) / [Therapeutics Data Commons](https://tdcommons.ai/)：化学信息学、机器学习与药物发现基准。
- [CHOPCHOP](https://chopchop.cbu.uib.no/)：基因编辑相关分析工具。
