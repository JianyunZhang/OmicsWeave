# 多组学联合分析使用指南

## 概述

本模块提供一套**解耦合、YAML 配置驱动的多组学联合分析框架**。MFA 反应作为生物学锚点，结论由连续的多维证据画像表达，这是唯一的整合结果；旧版的 0–4 加和总分与 High/Medium/Low 置信度已删除（阶段 2.8）。

旧版（v2 打分矩阵）技术文档已归档到 [`docs/history/integration_pipeline_evaluation_v2.md`](../../../docs/history/integration_pipeline_evaluation_v2.md)，仅供追溯。

### 支持的组学层

| 组学 | 输入格式 | 主要分析内容 |
|---|---|---|
| **代谢流分析 (MFA)** | 差异 Flux CSV + 通路 CSV + TCA CSV | 差异代谢通量、通路热图、TCA 循环 |
| **转录组学 (Transcriptomics)** | DE CSV + GSEA CSV | 差异基因、通路收敛 |
| **蛋白质组学 (Proteomics)** | 差异分析 CSV + 通路汇总 CSV | 差异蛋白（含比较名翻转逻辑）、通路收敛 |
| **代谢组学 (Metabolomics)** | 差异分析 CSV + 通路 CSV | 差异代谢物统计、通路收敛 |

### 核心分析功能

1. **反应级多维证据画像** — 输出 support、conflict、coverage、stability、design grade 和 evidence grade
2. **统一 EvidenceTable** — 区分 no_data/no_mapping/no_match/no_signal/opposite，并保留 effect、SE、qvalue、n_eff、映射质量和证据模式
3. **MFA 显著性判定** — 通过判定的反应作为方向锚点；判定只影响 MFA 层，不影响其他组学的证据评估
4. **跨组学通路融合** — 保留通路名称交集，并新增 Fisher 统计融合与 RRA 排名融合
5. **数据可用性热图** — 各组织 × 各组学数据覆盖
6. **差异特征汇总** — 各组织各组学上调/下调特征数
7. **组织整合看板** — 多组学层联合看板
8. **MFA 可视化与 Markdown 报告** — 火山图、通路图、TCA 图和证据摘要
9. **组学内多批次集成** — MFA/转录组读取上游 protocol 级结果，蛋白质组/代谢组按 effect+SE 做随机效应 Meta
10. **样本级协调** — 显式 sample map、组学特异归一化、缺失审计与保留 condition 的 batch 校正
11. **跨实验项目集成** — estimand 完整校验、随机效应 Meta、异质性审计与项目级唯一报告

---

## 模块结构

```
omicsweave/integration/
  config.py                              # MultiOmicsConfig dataclass + YAML 加载器
  harmonization.py                       # 样本映射、归一化、缺失与批次校正
  project_integration.py                 # estimand 校验 + 跨实验随机效应 Meta
  input_audit.py                         # 配置输入与实际批次来源审计
  mappings/
    __init__.py
    feature_resolver.py                  # 标准 ID 优先的基因/蛋白/代谢物精准匹配
    audit.py                             # 映射结构、ID 冲突与覆盖审计
    flux_gene_mapping.py                 # 101 反应 → 基因映射
    flux_protein_mapping.py              # 86 反应 → 酶/蛋白质映射
    flux_metabolite_mapping.py           # 111 反应 → 代谢物映射（底物/产物/关联）
  loaders/
    __init__.py
    mfa_loader.py                        # MFA 差异通量 + 通路 + TCA + 概览
    transcriptomics_loader.py            # DE + GSEA + Enrichr
    proteomics_loader.py                 # 蛋白质 DE + 通路（含比较名翻转）
    metabolomics_loader.py               # 代谢物 DE + 通路
  methods/
    __init__.py
    batch_effect_meta.py                 # 蛋白/代谢组结果级批次 Meta
    evidence_profile.py                  # EvidenceTable + support/conflict/coverage + 证据分级
    flux_evidence.py                     # flux-anchored 证据构建器：反应 → 基因/蛋白/代谢物的证据提取
    pathway_convergence.py               # 名称交集 + Fisher/RRA 融合 + BH FDR
    cross_tissue_pathway_scorer.py       # 跨组织策展通路评分
  modes/
    flux_anchored.py                     # omicsweave integrate：统一差异结果表 → 证据画像
    pathway_convergence.py               # omicsweave integrate：统一富集结果表 → 通路融合
  integrate.py                           # omicsweave integrate 入口
  plotting/
    __init__.py
    _config.py                           # 绑图全局状态（统一阈值、颜色、标签）
    overview_plots.py                    # 数据可用性热图 + 差异汇总 + 组织看板
    volcano_plots.py                     # MFA 火山图
    heatmap_plots.py                     # 通路热图 + TCA + 收敛 + 通量概览
    pca_plots.py                         # PCA 散点图 + scree plot
    scoring_plots.py                     # 证据等级分布 + 支持度/冲突度 + 各组织支持度
  report.py                              # Markdown 报告生成
  run_multiomics_analysis.py             # CLI 入口
  run_multiomics_harmonization.py        # 样本级协调 CLI
  run_multiomics_project_analysis.py     # 项目级跨实验 CLI
  examples/                              # 三类配置模板

# 配置文件（自包含，独立于其他组学配置）
experiments/{experiment_id}/
  multiomics_integration_setting/
    config.yaml

# 输出
analysis_results/{experiment_id}/
  multiomics_results/
    common/                              # EvidenceTable、反应画像、通路画像
      input_audit.csv                    # 配置批次与实际载入来源核对
      mapping_integrity_audit.csv        # 映射错误、警告与覆盖指标
    tables/                              # 通路证据融合与兼容交集表
    01_do_results/                       # 按组织和比较拆分的数据表
    02_result_plot/                      # 图表
    multiomics_report.md                 # 分析报告

analysis_results/{project_id}/
  project_multiomics_report.md           # 多实验项目级唯一报告
  experiment_manifest.csv
  compatibility_audit.csv
  evidence_long.csv
  meta_analysis_results.csv
  reaction_summary.csv
```

---

## 按统一项目模型整合（阶段 2.3）

```bash
omicsweave run --experiment <id>                                   # 各组学写出统一结果表
omicsweave integrate --experiment <id> --mode flux_anchored        # 读取最近一次运行的结果表
```

`omicsweave integrate` 只读取结果包中的统一结果表，结果写入新的结果包（`integration/<模式>/`），运行清单记录源结果包。
模式默认取 `experiment.yaml` 的 `integration.modes`，参数写在 `integration.parameters.<模式>`。

- `pathway_convergence`（已接入）：同一组织、同一比较下各组学通路富集结果的融合（组学内 Bonferroni、组学间 Fisher、BH 校正，
  与本模块 `PathwayConvergenceAnalyzer` 相同），附各组学的效应方向与方向一致性，以及显著通路的两两交集。
- `flux_anchored`（已接入）：以 MFA 差异通量为锚点，读取各组学的统一差异结果表（`differential.parquet`），
  按反应 → 基因/蛋白/代谢物映射提取证据，输出多维证据画像 `flux_evidence_profile`、证据长表 `flux_evidence_table`
  与通路汇总 `flux_pathway_profile`。需要 MFA 与至少 1 个其他组学。参数 `pval_cutoff`、`log2fc_cutoff`、`omics_cutoffs`、
  `fallback_mode`、`evidence_profile` 与旧配置的 `scoring` 同义（`evidence_profile.mfa_inference_unit` 默认取 MFA 结果的推断单位），
  `tissue_map` 把 MFA 的组织名对应到其他组学。与旧命令的差别：MFA 的 p 值与 FDR 取差异分析的主检验路径
  （旧加载器总是优先取 Mann-Whitney）；蛋白组、代谢组的有效样本量取自统一结果表的 `n_numerator + n_denominator`。
  证据等级随锚点状态而定（见下文“锚点状态与证据等级”）；样本 QC 标记按（组学, 组织, 比较）降权或判为不可评估，
  明细写在 `integration/flux_anchored/qc/contrast_reliability.parquet`（参数 `qc_imbalance`，见“质量信息与可靠性”）。
- `reversal`（已接入）：对 `design.reversal_pairs` 的每个配对，逐特征（含代谢组派生特征）计算疾病效应、处理效应、
  逆转指数及区间，分为 reversed / partially_reversed / over_reversed / aggravated / not_reversed / undetermined。
  “未逆转”须通过等效性检验（双单侧检验，界限为疾病效应的 `equivalence_margin` 倍，默认 0.5），低功效的不显著结果
  只能是 undetermined；疾病效应本身未建立（padj ≥ `alpha`）时不做判定。输出 `tables/reversal.parquet` 与
  `tables/reversal_summary.parquet`。通路层面的逆转与析因交互项尚未实现（需要公式模型）。
- `sample_level`：待接入，暂用下文的旧命令。

- `flux_anchored` 的同工酶加权（`isozyme.minor_share`，默认 0.2）：转录组与蛋白组的方向成员按该组织的平均丰度
  （RNA 的 `base_mean`、蛋白两组均值）加权，给出 `tx_/prot_weighted_capacity_change` 与主要成员；与通量同向的只有
  丰度占比低于阈值的次要同工酶、而主要同工酶反向时，该层不计为支持，记在 `minor_isozyme_driven_layers`。
- `flux_anchored` 的诱饵映射（`decoys: {replicates: N, seed: S}`，默认关闭）：把各反应的方向成员换成同组织中丰度相近的
  随机特征后用同一套规则判定，诱饵支持率与目标支持率之比估计假证据比例，写在 `evaluation/decoy_*.parquet`。
- 源结果包有逐只去掉动物的稳健性表时（`omicsweave stats`），与通量同向的特征全部由单只动物驱动的层记在
  `single_subject_driven_layers`。

项目级的跨实验 / 跨队列复现（`project.yaml` 的 `replication`）：

```bash
omicsweave replicate [--set <复现组>]
```

每个复现组列出估计对象相同的比较（可来自其他项目目录，组织名可用 `tissue_map` 对应），对特征、派生特征与 GSEA 通路做
方向性的部分合取检验与 BH 校正（`replicated`、`replication_fdr`），同层特征另给随机效应合并估计；`cross_layer` 允许
一个来源的蛋白对另一个来源的 RNA。`replication_baseline.parquet` 给出全部对应特征的同号比例与随机重排的基线（约 50%）。

结果包工具：

```bash
omicsweave diff <运行A> <运行B> [--experiment <id>] [--output <目录>]   # 0 相同 / 1 有差异 / 2 无法对照
omicsweave stats --experiment <id> [--source <运行>]                     # 对已有结果包重跑后处理（公式模型、派生特征、稳健性）
omicsweave refs list | build <genes|proteins|metabolites> | audit --experiment <id>   # 知识表与数据的审计
```

`diff` 按主键对齐每张表后逐格比较，报告逐值相同与否、变化的列、数值差值分布与 `evidence_grade` 等分类列的转移计数。
课题脚本读取结果包用 `omicsweave.contract.query`（`read_table(实验, "differential", omics="mfa")`、
`open_run(实验, "latest")`），读出的表在 `attrs["run_id"]` 中记录实际读取的运行，不必写死运行 ID。

---

## 快速开始

```bash
# 通过实验 ID 运行（自动查找 experiments/<id>/multiomics_integration_setting/config.yaml）
omicsweave integration run \
    --experiment 20250425_cimr_wt_mice_dio

# 或直接指定配置文件路径
omicsweave integration run \
    --config experiments/20250425_cimr_wt_mice_dio/multiomics_integration_setting/config.yaml

# 样本级协调（用于配对样本的潜变量/网络等补充分析）
omicsweave integration harmonize \
    --config src/omicsweave/integration/examples/harmonization.example.yaml

# 多个 estimand 可比实验的随机效应 Meta 与项目总报告
omicsweave integration project \
    --config src/omicsweave/integration/examples/project_integration.example.yaml
```

### 添加新实验

1. 复制现有 `experiments/<template>/multiomics_integration_setting/config.yaml`
2. 更新 `experiment_id`、`paths`、`conditions`、`mouse_groups`、`params`
3. 运行 `--experiment <new_id>`

### 配置文件格式

配置采用 `config.yaml` + `groups.yaml` 双文件结构，位于 `experiments/{experiment_id}/multiomics_integration_setting/`，与转录组管道配置模式一致。

RNA 原生多批次联合模型还可在转录组 `groups.yaml` 中显式增加：

```yaml
sample_batches:
  "52": batch_1
  "53": batch_1
  "54": batch_2
  "55": batch_2
```

只要配置了 `sample_batches`，所有载入样本都必须有映射；不能从样本号或目录名自动猜测。

**`groups.yaml`** — 组别与小鼠映射（同 `transcriptomics_data/groups.yaml`）：

```yaml
conditions:
  hfd_14w:
    display: "HFD"
    diet: "D12492"
    diet_desc: "High-fat diet DIO 14 week"
    color: "#FF8000"
  lfd_14w: ...
  hsd_14w: ...

mouse_groups:
  "52": hfd_14w
  ...
```

**`config.yaml`** — 分析参数与路径：

```yaml
experiment_id: 20250425_cimr_wt_mice_dio
omics: multiomics

paths:
  mfa_results: "analysis_results/20250425_cimr_wt_mice_dio/mfa_results"
  metabolomics_results: "analysis_results/20250425_cimr_wt_mice_dio/metabolomics_results"
  proteomics_results: "analysis_results/20250425_cimr_wt_mice_dio/proteinomics_results"
  transcriptomics_results: "analysis_results/20250425_cimr_wt_mice_dio/transcriptomics_results"
  output_dir: "analysis_results/20250425_cimr_wt_mice_dio/multiomics_results"

params:
  tissues: [liver, heart, kidney, gast_muscle, sq_adipose]
  comparisons: [HFD_vs_LFD, HSD_vs_LFD]
  mfa_tissue_map: {muscle: gast_muscle, white_adipose: sq_adipose}
  comparison_map: {HFD_vs_LFD: DIO_vs_WT, HSD_vs_LFD: DIO_Drug_vs_WT}
  scoring:
    pval_cutoff: 0.05
    log2fc_cutoff: 0.58
    effect_tiers:
      log2fc: [0.58, 1.0]     # [moderate下限, large下限]
      cohen_d: [0.8, 1.5]     # [moderate下限, large下限]
    confidence_thresholds: {high: 3.0, medium: 2.0, low: 1.0}
    evidence_profile:
      mfa_inference_unit: mcmc_fit_points  # 上游修正后改为 biological_replicate
      support_threshold: 0.15
      conflict_threshold: 0.15
      high_conflict_threshold: 0.25
      effect_scales:
        mfa: 1.0
        transcriptomics: 1.0
        proteomics: 1.0
        metabolomics: 1.0
  tissue_labels:
    heart: "心脏 Heart"
    liver: "肝脏 Liver"
    kidney: "肾脏 Kidney"
    gast_muscle: "腓肠肌 Gastrocnemius"
    sq_adipose: "皮下脂肪 SQ Adipose"
```

路径支持三种格式：绝对路径 → repo 相对路径 → config 相对路径。

多批次输入契约见 [`examples/integration_inputs.fragment.yaml`](examples/integration_inputs.fragment.yaml)。
样本级协调和项目级 Meta 分别见 [`examples/harmonization.example.yaml`](examples/harmonization.example.yaml)
与 [`examples/project_integration.example.yaml`](examples/project_integration.example.yaml)。

---

## 精准特征识别与映射审计

反应评分统一通过 `mappings/feature_resolver.py`，不允许各 scorer 自行实现名称匹配。

匹配优先级：

1. **mRNA**：Ensembl Gene ID → NCBI Gene ID → 规范 MGI 小鼠符号/历史别名；
2. **蛋白质**：UniProt 主 accession → 规范小鼠基因符号。只有目标蛋白尚无本地
   UniProt 注释时，才允许在实测行已有 accession 的情况下按符号安全回退；
3. **代谢物**：HMDB → KEGG Compound → ChEBI → 规范化后的完整别名精确匹配；
4. 数据行只要含有效标准 ID，就以 ID 为准，不再使用名称覆盖冲突；
5. 禁止 `str.contains` 任意子串匹配。`Pi`、`Glu`、`Ser`、`GAP` 等短别名只有
   在实测名称完整相等时才可命中；
6. 同一代谢物在同一反应中命中多个角色时标记为 ambiguous，从所有方向评分中排除；
7. 配置中没有对应 `log2fc_<comparison>` 的蛋白/代谢物层明确记为 `no_data`；
8. `params.comparisons` 是 MFA 差异表的白名单，目录中的历史比较不会自动进入评分。

`audit_mapping_integrity()` 在 scorer 初始化时执行。结构性错误会阻断运行；非阻断的
注册表覆盖缺口写入 `common/mapping_integrity_audit.csv`。EvidenceTable 的
`mapping_methods`、`n_ambiguous`、`ambiguous_features` 和 `data_reason` 可直接审计
每一层的匹配来源、歧义排除与缺失原因；反应表仍保留完整的 `*_details`。

方向证据采用保守规则：只有直接酶或直接运输执行者可产生 support/conflict；
转录因子、抑制因子、upstream/downstream pathway neighbor 仅记录为
`unresolved_signal`，不能因为与通量同向就获得方向分。

---

## 多维证据画像（主结果）

每个反应、比较和组织生成四行 EvidenceTable（MFA、转录组、蛋白质组、代谢组）。单层连续证据为：

```text
layer_evidence =
    方向占比 × 映射集合内信号占比 × 方向可信度
  × 效应强度 × 映射质量 × 数据可靠性
```

- `sign_certainty` 优先使用 `1-lfsr`，其次使用 effect/SE，最后才使用 qvalue 或证据模式代理。
- `effect_strength` 使用组学层独立的 `effect_scales` 饱和归一化。
- 映射集合使用 Simes 聚合概率和信号占比，避免“映射对象越多越容易取到极值”。转录层优先用蛋白策展表中的 `direct_enzymes` 加直接 transporter 定义直接基因；蛋白层只有直接酶/直接 transporter 提供方向证据；其他关联分子保持 unresolved。代谢层区分底物、产物、related，并排除跨角色歧义。
- effect-only、n<3 和 MCMC 拟合点伪重复会降低 `data_reliability` 并限制 `design_grade`。
- 缺失层不进入 support/conflict 的分母，但会降低 `coverage_score`。

主要输出：

| 字段 | 范围/等级 | 含义 |
|---|---|---|
| `support_score` | 0–1 | 可用且可评估组学中支持预期机制的连续证据 |
| `conflict_score` | 0–1 | 可用且可评估组学中与预期方向冲突的连续证据 |
| `coverage_score` | 0–1 | 四个组学层实际完成数据、映射和检测匹配的比例 |
| `stability_*` | 0–1 或缺失 | 生物学重复 bootstrap 的选择概率和方向稳定性；无原始重复值时不填充 |
| `design_grade` | A/B/C/D | 独立于生物学一致性的实验设计质量 |
| `anchor_type` / `anchor_state` | 分类 | 锚点类型（目前为 `mfa_flux`）与状态，见下文 |
| `omics_direction` | up/down/mixed/none | 各组学层提示的通量方向，与锚点符号无关 |
| `evidence_grade` | 分类 | 锚点成立时：Convergent、Single-layer、Exploratory、Conflicting；锚点不成立时：Omics-concordant、Omics-single-layer、Omics-discordant；以及 No-signal、Not-assessable |

`Replicated` 只应在至少两个 estimand 可比的实验完成跨实验复现后赋值；单实验的证据构建器不会自行产生该标签。

### 锚点状态与证据等级

其他组学只有在锚点成立时才按锚点方向判支持或冲突。`anchor_state` 依次判定：

| 状态 | 条件 | 等级词汇 |
|---|---|---|
| `qc_confounded` / `qc_insufficient` | 可靠性表中 MFA 层在该（组织, 比较）不可评估 | 组学一致性 |
| `not_identifiable` | 差异结果 `flags` 含 `not_identifiable`，或可辨识性为 poorly_constrained / bound_pinned | 组学一致性 |
| `tracer_conflict` | 多示踪剂方向不一致 | 组学一致性 |
| `insufficient` | MFA 证据等级为 insufficient，或较小组少于 2 只 | 组学一致性 |
| `not_established` | 未通过 MFA 显著性判定 | 组学一致性 |
| `effect_only` | 按效应量通过判定（设计等级 D，最多 Exploratory） | 锚定 |
| `established` | 按统计检验通过判定 | 锚定 |

组学一致性等级只看各组学层提示的通量方向是否一致：直接酶 / 转运体升 → 通量升，产物升 → 通量升，底物升 → 通量降。
≥ 2 层同向为 Omics-concordant，1 层为 Omics-single-layer，层间或层内方向相反为 Omics-discordant。这些方向与锚点的
符号无关（把锚点符号翻转，组学一致性等级不变，有测试保证）。`support_score`、`conflict_score` 仍按通量点估计的方向计算，
锚点不成立时只作参考。EvidenceTable 每层另有 `implied_direction`、`n_implied_up`、`n_implied_down` 与 `anchor_state`。

### 当前实施边界

- 已实现：EvidenceTable、多维画像（唯一整合结果）、MFA 显著性判定与其他组学的证据评估解耦、通路 Fisher/RRA 融合、comparison 隔离、跨组织策展评分降级和输入来源审计。
- 已实现：MFA 每只动物只贡献一个通量估计，输出 log2FC SE、Hedges' g；可选 biological-replicate bootstrap/标签置换。历史产物只有在重新运行 MFA 后才会具备这些字段。
- 已实现：RNA 显式 `sample_batches` 时使用可辨识的 `~batch + condition` 设计；完全混淆会拒绝推断。也可逐批分析后用 protocol 级随机效应 Meta。
- 已实现：蛋白质组/代谢组可读取多个已独立完成 QC/差异分析的批次目录，按同一效应尺度与 SE 做随机效应 Meta；缺少不确定性的批次只记为 observed，不计入 pooled 推断。
- 已实现：样本级 harmonization CLI 和跨实验 project integration CLI。样本级流程输出校正前/后的 PCA、控制 condition 后的 batch partial R²、缺失审计和 `harmonization_report.md`；项目级 `Replicated` 只由 estimand 完全一致的实验生成。
- 尚不能自动计算的字段：现有差异结果表若没有生物学重复原始值，`stability_*` 保持缺失；不会用 MCMC 拟合点填充。
- 仍需平台上游完成：代谢组 QC-RLSC/run-order 漂移校正、蛋白质组统一谱图库/iRT QC、MFA MID 与网络参数 QC。集成层不会把未经这些 QC 的结果数学“修复”为可比。
- 项目级 Meta 是显式独立入口，不会被单实验流程自动触发。物种、组织、case/control 含义、时间、示踪剂、MFA 网络、平台或主要终点不同的实验只并列报告，不合并总分。
- MOFA+/DIABLO/GEM 仍是条件性补充策略；当前 n=1–4 不满足默认潜变量建模要求。

### 样本级协调与残余批次诊断

[`examples/harmonization.example.yaml`](examples/harmonization.example.yaml) 使用显式
`sample_map` 统一不同平台的样本 ID，并可在 condition×batch 设计可辨识时移除 batch 项。
`batch_qc` 默认开启，控制 condition 后计算特征级和 PCA 级 batch partial R²：

```yaml
batch_qc:
  enabled: true
  max_components: 5
  max_features: 5000
  high_batch_partial_r2: 0.20
```

输出包括 `missingness_qc.csv`、`batch_qc.csv`、`batch_feature_qc.csv`、
`batch_pca_scores.csv`、`batch_pca_variance.csv`、`batch_pca/*.png` 和
`harmonization_report.md`。单批数据标为 `single_batch`，混淆设计标为 `confounded`
并拒绝校正，不会伪造跨批可靠性。

### 质量信息与可靠性

样本层面的 QC 结论写在样本表中，按（组学, 组织, 比较）自动变成可靠性（`omicsweave.core.reliability`）：

```text
sample_id        subject_id  omics            tissue  qc_flag  qc_reason
rna_liver_m63    m63         transcriptomics  liver   flag     线粒体读数占比 54%
```

`qc_flag` 取 `flag` 或 `exclude`（留空为正常）；各组学的 `qc/sample_quality.parquet`（样本身份与组成 QC 的输出）中
`flag` 非空的样本同样计入。规则：两组都没有标记 → `ok`；两组标记比例相差 ≥ `qc_imbalance`（默认 0.5）→ `qc_confounded`；
去掉标记样本后有一组少于 2 只 → `qc_insufficient`；否则 `downweighted`，可靠性为可用样本占比。后两种不可评估的状态使该层
不判支持或冲突（MFA 层则锚点不成立）。`omicsweave run` 结束时写出结果包的 `qc/contrast_reliability.parquet`，并在不可评估的
差异结果行的 `flags` 中加上状态；整合按当前样本表重新计算，所以改样本表后只需重跑整合。

审核人员在配置中显式给出的层级因子作为覆盖项，替换可靠性表中该层的数值（不解除不可评估）：

```yaml
integration:
  parameters:
    flux_anchored:
      evidence_profile:
        quality_reliability: {transcriptomics: 0.5}
```

该因子写入 EvidenceTable 的 `quality_reliability` 并乘入 `data_reliability`，可靠性状态写在 `qc_status`。
`batch_qc.csv` 的 `suggested_batch_reliability` 只是审核建议，不会被流程自动复制到评分配置。

## MFA 显著性判定（方向锚点）

通过判定的 MFA 差异反应在证据画像中 MFA 层记为 support，作为其他组学方向判断的锚点；未通过的反应 MFA 层记为 no_signal，
其他组学照常评估，但只给组学一致性等级。通过判定是锚点成立的必要条件而非充分条件：样本量、可辨识性、示踪剂一致性与 QC
也要满足（见“锚点状态与证据等级”）。旧版在此之上的 0–4 加和总分、固定权重与 High/Medium/Low 置信度已在阶段 2.8 删除，不再输出 `legacy_*` 字段。

**判定条件：FDR < 0.05 且 |Cohen's d| > 0.8（大效应量）**。样本量不足以给出 p 值时（如每组 1 只），|Cohen's d| ≥ 0.8 且
|log2FC| ≥ 0.58 按 `effect_size_only` 通过（`fallback_mode`）。

| 指标 | 阈值 | 说明 |
|---|---|---|
| FDR（Benjamini-Hochberg 校正后 p 值） | < 0.05 | 控制多重检验假阳性 |
| Cohen's d（效应量） | \|d\| > 0.8 | 要求大效应量（Cohen 标准），而非仅统计显著 |

**限制：** 旧版 MFA 差异分析的统计检验基于 MCMC 拟合点（每样本 ~150 次迭代，IQR 过滤后 ~75 点），而非生物学重复（n=3-4）。将同一只小鼠的拟合点作为独立样本会造成伪重复。FDR + Cohen's d 双阈值只是启发式筛选，不能修复独立观测单位错误；拟合点层面的推断在证据画像中标记为 `pseudo_replicate`，只给探索性设计等级。`omicsweave run` 的 MFA 差异分析以动物为推断单位，`flux_anchored` 据此不再降级。

Gate 结果示例（本项目 HFD/HSD/LFD 实验，1191 个反应）：

| 比较组 | 总反应数 | 通过 Gate | 比例 |
|---|---|---|---|
| HFD_vs_LFD | 397 | 52 | 13.1% |
| HSD_vs_LFD | 397 | 114 | 28.7% |
| HFD_vs_HSD | 397 | 82 | 20.7% |

---

## 列名规范

所有 loader 在加载时统一列名：

| 概念 | 规范列名 | 原始变体 |
|---|---|---|
| 通量反应 ID | `flux_id` | flux_id, Flux_ID |
| 组织 | `tissue` | tissue, Tissue, tissue_tx |
| 通量 log2FC | `flux_log2fc` | log2FC |
| 通量 p 值 | `flux_pvalue` | pvalue_mannwhitney, pvalue_ttest |
| 通量 FDR | `flux_fdr` | fdr_mannwhitney, fdr_ttest |
| 通量效应量 | `cohen_d` | cohen_d |
| 基因名 | `gene` | gene |
| 基因 log2FC | `log2fc` | log2FoldChange |
| 蛋白名 | `protein` | Protein |
| 代谢物名 | `compound` | Compounds |

---

## 输出文件

### 数据表 (`common/` 与 `tables/`)

| 文件 | 说明 |
|---|---|
| `common/flux_evidence_profile.csv` | 每个（反应, 组织, 比较）的多维证据画像 |
| `common/evidence_table.csv` | 每个反应 × 组学层的统一长表，含匹配来源、歧义和缺失原因 |
| `common/mapping_integrity_audit.csv` | 映射结构错误、ID 冲突、注册表覆盖与别名统计 |
| `common/pathway_evidence_profile.csv` | 通路级 support/conflict/coverage 汇总与各证据等级的反应数 |
| `tables/pathway_evidence_long.csv` | 通路融合前的组学来源、概率和排名长表 |
| `tables/pathway_fusion.csv` | Fisher/RRA 融合结果及 BH FDR |
| `tables/pathway_convergence.csv` | 兼容用跨组学共享通路详细列表 |
| `tables/pathway_convergence_counts.csv` | 兼容用跨组学共享通路计数 |
| `multi_omics_overview.csv` | 数据可用性矩阵 |
| `differential_summary.csv` | 差异特征计数 |
| `tissue_dashboard_{comp}.csv` | 组织看板数据 |

### 图表 (`plots/`)

| 文件 | 说明 |
|---|---|
| `multi_omics_overview.png` | 数据可用性热图 |
| `mfa_volcano_{comp}.png` | MFA 火山图 |
| `mfa_pathway_condition_heatmap.png` | 通路-条件热图 |
| `mfa_tca_cross_tissue_heatmap.png` | TCA 循环热图 |
| `evidence_grades.png` | 证据等级分布 |
| `support_conflict.png` | 每个反应的支持度与冲突度 |
| `support_by_tissue.png` | 各组织支持度箱线图 |
| `pathway_convergence_heatmap.png` | 通路收敛热图 |
| `differential_summary_{comp}.png` | 差异汇总柱形图 |
| `tissue_dashboard_{comp}.png` | 组织看板 |

### 报告

| 文件 | 说明 |
|---|---|
| `multiomics_report.md` | 完整分析报告 |

---

## 常见问题

### 1. 某组学数据缺失

模块设计为各组学可选加载。缺少某组学数据时，该组学的分析步骤自动跳过，不影响其他分析。

### 2. 比较组名称不一致

不同组学使用不同比较名（如 `HFD_vs_LFD` vs `DIO_vs_WT`），通过 `comparison_map` 建立映射。蛋白质组 loader 自动处理反向比较名（如文件中存为 `LFD_vs_HSD` 时自动翻转 log2FC 符号）。

### 3. 组织名称映射

通过 `mfa_tissue_map` 将 MFA 组织名映射到项目规范名称（如 `muscle` → `gast_muscle`）。所有组织名内部统一为小写下划线。

### 4. 图表中文显示

自动尝试 SimHei、Microsoft YaHei、DejaVu Sans 字体。如中文显示为方框，确保系统安装了中文字体。

---

## 多组学集成技术方案评估

### 项目定位

本项目以 **代谢流分析（MFA/Fluxomics）作为代谢生理学研究的主要研究方法**，采用 **flux-centric 多策略集成**：反应级证据画像保留 MFA 锚点的可解释性，通路级 Fisher/RRA 融合降低一对多映射噪声；稳定性、设计质量和覆盖度与生物学支持强度分开表达。

### 为什么选择 MFA 作为 baseline

从生物学层级看因果方向：

```
基因转录(mRNA) → 蛋白质(酶) → 酶活性(翻译后调控) → 代谢通量(MFA) → 代谢物丰度
     ↑                  ↑              ↑                     ↑                  ↑
   最间接            中等间接        无法直接测量          直接功能测量        受多层调控
```

MFA 是唯一能直接测量代谢通路**实际功能输出**的技术。以功能输出为锚点，向上游追溯调控原因（基因→蛋白→酶→通量），是因果推断方向上**从果到因**，比反过来（从 mRNA 预测通量）更合理。

各组学在 flux-centric 集成中的角色：

| 组学 | 角色 | 回答的问题 |
|---|---|---|
| **MFA (baseline)** | 锚定功能变化 | 哪些代谢通路的**实际功能**发生了改变？ |
| **蛋白质组** | 解释酶资源 | 通量改变是否由酶丰度变化驱动？ |
| **转录组** | 解释调控源头 | 酶丰度改变是否源于转录调控？ |
| **代谢组** | 验证底物/产物 | 通量改变是否反映在代谢物池的变化上？ |

### 主流多组学集成方法对比

以下对当前主流的多组学集成方法进行评估，并判断其在本项目条件下的适用性。

| 方法 | 类型 | 核心原理 | 本项目适用性 | 核心限制 |
|---|---|---|---|---|
| **连续证据画像** | Late, 知识驱动 | 以 MFA 为锚点，逐反应区分支持、冲突、覆盖和设计质量 | **默认策略** | 仍依赖映射质量；不能替代跨实验复现 |
| **约束代谢模型 (GEM/INTEGRATE)** | 知识驱动 | 将组学数据转化为代谢网络模型约束条件，通过 FBA 预测通量并区分调控层次 | **中长期最佳升级方向** | 需构建小鼠基因组规模代谢模型 |
| **通路级集成** | Late, 知识驱动 | Fisher 融合统计证据；只有排序时使用 RRA | **默认补充策略** | 分辨率仅到通路水平，依赖术语标准化 |
| **MOFA/MOFA+** | Early, 数据驱动 | 贝叶斯因子分析，将多组学矩阵分解为共享/特异潜变量因子 | **当前不作为默认** | 可处理缺失值和缺失 assay，但当前 n=1–4 难以稳定估计潜变量 |
| **SNF (Similarity Network Fusion)** | 网络融合, 无监督 | 对每个组学构建样本间相似性网络，迭代融合为统一网络 | **不适用** | 用于样本分型/亚群发现，本项目已有明确分组；n=3 太少 |
| **DIABLO (mixOmics)** | Intermediate, 监督 | 基于偏最小二乘（PLS）的监督学习，在多组学间寻找与分组变量最大协方差的特征子集 | **可作补充分析** | 可寻找跨组学特征模块，但不能回答因果方向 |
| **机器学习通量预测** | 数据驱动 | 用转录组/蛋白质组作为特征，MFA 通量作为标签训练监督学习模型 | **不适用** | 样本量严重不足（n=3-4/组） |

#### 参考文献

- **MOFA+**: Argelaguet R, et al. (2020). MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology*.
- **SNF**: Wang B, et al. (2014). Similarity network fusion for aggregating data types on a genomic scale. *Nature Methods*.
- **DIABLO**: Singh A, et al. (2019). DIABLO: an integrative approach for identifying key molecular drivers from multi-omics studies. *Bioinformatics*.
- **INTEGRATE**: Di Filippo M, et al. (2022). Model-based multi-omics data integration to characterize multi-level metabolic regulation. *PLoS Comput Biol*. PMID: 35130273
- **GEM 多组学集成综述**: Shahreen N, et al. (2026). Multi-omics integration in genome-scale metabolic models: a review of constraint-based approaches. *Mol Omics*. PMID: 41660891
- **ML 通量预测**: Gonçalves DM, et al. (2023). Predicting metabolic fluxes from omics data via machine learning. *Comput Struct Biotechnol J*. PMID: 37876626

### 方向一致性判断

对每条 MFA 差异反应 × 每个比较 × 每个组织，证据构建器在各组学层判断映射特征的变化方向是否与通量一致：
转录组与蛋白组只用直接酶 / 直接转运体，代谢组只用底物与产物，其余关联分子只计入匹配，不提供方向。

#### 方向一致性判断的生物学逻辑

对酶蛋白：酶上调 + 通量上调 = 一致性（正调控）。

对代谢物需区分角色：
- **底物**：通量上调 → 底物被消耗 → 丰度**下调**才是"一致"
- **产物**：通量上调 → 产物增加 → 丰度**上调**才是"一致"

不能简单用"同方向 = 一致"判断代谢物与通量的关系。

#### 不一致信号的生物学意义

不一致不代表错误。酶活性调控（磷酸化/别构调控）可能在不改变蛋白质丰度的情况下改变通量。多层不一致本身也是有意义的生物学信号，提示翻译后调控或底物可用性等非转录水平的调控机制。

---

## 补充样本批次校准方案

### 核心原则：预防 > 校正

批次效应一旦产生，数学方法只能缓解无法完全消除。**实验设计阶段的预防措施是投入产出比最高的。**

### 一、实验设计阶段的准备

#### 1. 样本随机化（必须做）

```
错误做法（批次 = 组别，完全混淆）：
  Batch 1 (旧): HFD_1, HFD_2, HFD_3 | LFD_1, LFD_2, LFD_3
  Batch 2 (新): HFD_4, HFD_5, HFD_6 | LFD_4, LFD_5, LFD_6

正确做法（批次 ≠ 组别，统计上可分离）：
  Batch 1: HFD_1, HFD_3, LFD_2, LFD_4  (随机分配)
  Batch 2: HFD_2, HFD_4, LFD_1, LFD_3  (随机分配)
```

新增样本必须在**每批次中同时包含各饮食组**，不要整组整批跑。

#### 2. QC 样本准备（必须做）

| 组学 | QC 样本 | 作用 |
|---|---|---|
| **RNA-Seq** | ERCC spike-in mix（等量加入每个样本） | 跨批次定量校准 |
| **代谢组** | Pooled QC（所有样本等量混合，每 5-10 针插一针） | 监测仪器漂移 |
| **蛋白质组** | iRT 标准肽 + Pooled QC | 保留时间对齐 + 信号漂移校正 |
| **MFA** | 同位素标记对照（如 U-13C 葡萄糖平行处理） | 标记效率跨批次校准 |

#### 3. Bridge Samples（强烈建议）

在新 batch 中至少包含 2-3 个旧 batch 的**相同样本**（从同一 aliquot 重新提取上机），用于量化批次差异量级。Bridge samples 是连接新旧批次最可靠的方式，比任何数学校正都更可信。

#### 4. 实验条件标准化清单

```
□ 同一操作人员或严格统一 SOP
□ 同一批试剂（裂解液、酶、抗体等）
□ 同一仪器参数（质谱条件、测序平台）
□ 样本冻融次数一致
□ 提取到上机的间隔时间一致
□ 详细记录每个样本的处理日期、操作人员、试剂批号
```

### 二、各组学批次校正方法

#### 1. RNA-Seq（转录组）

**当前 pipeline**：支持两条显式路径。

- 同一计数矩阵且设计可辨识：在转录组 `groups.yaml` 配置完整 `sample_batches`，PyDESeq2 自动使用 `~batch + condition`。
- 各批已独立分析或 batch 与 condition 混淆：逐批生成 effect/SE，再由 `transcriptomics_integration.py` 做随机效应 Meta；不会回退为忽略 batch 的 Wilcoxon 推断。

联合模型示意：

```python
design = "~batch + condition"

# rnaseq_pipeline.py 中修改
dds = DeseqDataSet(
    counts=counts_filtered,
    metadata=metadata_sub,
    design="~batch + condition",
)
```

前提：batch ≠ condition 完全混淆（即每个 batch 中有各种 condition）。

为什么用模型内校正而非 ComBat：DESeq2 在模型内处理 batch 效应比先 ComBat 再 DE 更准确，因为 ComBat 针对连续数据设计，对 count 数据需要用 ComBat-seq（NB 分布假设），不如直接在 GLM 中加协变量简洁。

方案 B（备选）— ComBat-seq：适用于极端情况（如 3 个以上批次，或批次×条件不平衡）。

配置要求：在 `groups.yaml` 增加 `sample_batches: {sample_or_mouse_id: batch_id}`；只要该映射非空，每个实际载入样本都必须有批次身份。代码会检查设计矩阵秩，完全混淆时拒绝拟合。

#### 2. 代谢组学

**当前 pipeline**：log2 变换 + 半最小值填补，无批次校正。

**新增批次后需要**（分两步）：

Step 1 — 仪器漂移校正（QC-RLSC）：
利用 Pooled QC 样本，对每个代谢物做局部加权回归（LOESS），将所有样本信号校正到 QC 水平线。这是代谢组学跨批次校正的标准操作。

Step 2 — 批次间归一化：
- 方案 A: PQN (Probabilistic Quotient Normalization) — 批次间
- 方案 B: ComBat — 在 log2 变换后应用
- 方案 C: SERRF (Systematic Error Removal using Random Forest)

需要提前准备：
- **必须在跑样时插入 Pooled QC**（每 5-10 个样本插一针）
- 记录每个样本的 run order（上机顺序）
- 保留所有原始峰面积数据（不要只保留归一化后的）

#### 3. 蛋白质组学（DIA）

**当前 pipeline**：log2 + 中位数中心化 + MNAR 填补，无批次校正。

**新增批次后需要**：

方案 A（推荐）— limma removeBatchEffect 或模型内 batch 因子：
在差异分析之前用 limma 的 `removeBatchEffect()` 去除批次效应。注意 limma 本身支持在模型中加入 batch 因子，优先用模型内方法。

方案 B — ComBat (parametric)：对 log2 变换后的蛋白质丰度矩阵应用 ComBat。DIA 蛋白质组学中 ComBat 是验证过的标准方法。

关键前提：
- 新旧批次必须用**相同的 DIA 谱图库**，或至少有足够的共同鉴定蛋白
- 如果不同批次用不同仪器，需额外的 inter-instrument normalization
- iRT 标准肽用于保留时间对齐

#### 4. MFA（代谢流分析）

**当前 pipeline**：每只动物汇总为一个通量估计后做差异推断，输出 log2FC SE、Hedges' g 和可选生物学重复稳定性；多个独立 MFA 批次可在同一 effect/SE 尺度上做随机效应 Meta。原始 MID 的跨批漂移校正与网络模型一致性检查仍必须在 MFA 上游完成。

MFA 的批次效应来源与丰度数据不同：
- 同位素标记效率差异（不同次 13C 葡萄糖配制）
- 质谱采集参数漂移（影响 MID 测量）
- INCA/Metran 拟合初值不同（影响收敛）

**推荐方案**：

Step 1 — 对原始 MID 数据做 QC 检查：比较新旧 batch 中相同条件样本的 MID 分布，如果同条件样本 MID 一致，无需校正。

Step 2 — 统一拟合参数：新旧 batch 用相同的 metabolic network model、统一 flux bounds 和 carbon mapping。

Step 3 — 如果 MID 存在系统性偏移：对 MID 做 ComBat 校正后重新拟合通量，或在通量水平做中位数/分位数归一化。

MFA 特殊注意：
- 13C 底物必须同一批配制，或至少标定标记纯度
- 如果旧 batch 样本已消耗无法重跑，新 batch 应包含重复测量（同一旧样本重新提取 + 上机），用于量化批次差异量级

### 三、多组学跨批次统一检查

各组学分别校正后，进入多组学集成前还需一步：

```
跨组学批次一致性检查：
  1. 对每个组织，画出 PCA/t-SNE 图，标注 batch 信息
     pca_plot(data, color_by="condition", shape_by="batch")
     理想结果: condition 分开，batch 内混杂
     问题信号: batch 分开，condition 内混杂

  2. 计算同条件、跨批次的配对相关性（如 HFD_1 vs HFD_4）

  3. 如果某组学的批次效应残余特别大，在打分矩阵中降权
```

### 四、补充样本操作清单

按时间顺序：

```
实验前（现在就要准备）:
  ☐ 设计随机化方案 — 新样本不在同一 batch 跑同一组全部
  ☐ 采购 QC 标准品 — ERCC spike-in, iRT 肽, U-13C 葡萄糖
  ☐ 制定 Pooled QC 方案 — 代谢组/蛋白质组混合 QC
  ☑ RNA 配置契约已支持 `sample_batches`；实验人员仍需填写真实批次映射
  ☐ 旧样本保留备份 — 至少各组织保留 1-2 份用于跨 batch 重复测量

实验中:
  ☐ 记录每个样本的完整元数据（日期、人员、试剂批号、上机顺序）
  ☐ 按计划插入 QC 样本
  ☐ 新 batch 中加入 2-3 个旧 batch 相同样本作为 bridge sample

分析前:
  ☑ RNA-Seq pipeline: 显式多批次时使用 ~batch + condition，并拒绝混淆设计
  ☐ 代谢组 pipeline: 加入 QC-RLSC 漂移校正
  ☐ 蛋白质组 pipeline: 加入 ComBat 或 limma batch 校正
  ☐ MFA pipeline: 检查 MID 跨 batch 一致性
  ☑ 样本级 harmonization 输出跨组学 PCA 与 condition-adjusted batch partial R²

分析中（多组学集成）:
  ☑ EvidenceTable 记录 n_batches、batch_ids、I²、tau² 和 source_mode
  ☑ 导出残余 batch PCA/QC 与建议可靠性；仅在人为审核后通过 `quality_reliability` 应用
```

### 五、跨批次样本合并分析策略

#### 前提条件：批次与条件不可混淆

合并分析是否可行，取决于批次（batch）与实验条件（condition）之间是否存在完全混淆：

```
场景A — 可以合并（推荐）：
  Batch 1 (旧): HFD×3, HSD×3, LFD×3   ← 每批都有各种条件
  Batch 2 (新): HFD×3, HSD×3, LFD×3   ← 每批都有各种条件
  → batch ≠ condition，统计上可分离

场景B — 不能合并：
  Batch 1 (旧): HFD×3, HSD×3, LFD×3
  Batch 2 (新): GLP_1×3, Sleeve×3     ← 新条件只在新batch
  → 对 GLP_1/Sleeve 的分析，batch 效应无法与条件效应分离
  → 但 HFD/HSD/LFD 在两个 batch 都有，这部分仍可合并
```

只有实际新增批次在每个待比较条件中都有覆盖，且设计矩阵满秩时，才能采用场景 A 的联合模型；代码不会仅凭“增加 n”自动认定可合并。

#### 联合模型与结果级 Meta 的选择

| 维度 | 合并后统一分析 | 分开分析后随机效应 Meta |
|---|---|---|
| 原始数据与平台 | 同一特征空间、模型和可联合原始数据 | 各批只保留差异结果，或平台不宜直接拼矩阵 |
| 批次处理 | 模型内显式 `~batch + condition` | 批次内估计 effect/SE，随机效应汇总 |
| 异质性 | 由设计模型吸收可辨识 batch 项 | 显式输出 I²、tau²、Q 和方向一致性 |
| 完全混淆 | 不可拟合 | 只有各批内部都含相同对照且 estimand 一致时才可 Meta |
| 当前实现 | 已实现 | 已实现 |

#### 各组学的具体合并流程

**RNA-Seq：直接合并 count 矩阵**

```
推荐流程：
  1. 旧 batch FASTQ + 新 batch FASTQ → 分别 QC/trim
  2. 用同一个 kallisto index 量化所有样本
  3. 合并 count 矩阵
  4. PyDESeq2: design = "~batch + condition"

  不需要：对不同 batch 做 ComBat 预处理
  原因：DESeq2 的 NB GLM 天然支持协变量，在模型内处理 batch 效应比任何事后校正都准确
```

Pipeline 使用点：
- `groups.yaml` 增加 `sample_batches` 映射
- `rnaseq_pipeline.py` 自动验证设计秩并选择 `~batch + condition`
- `--all-batches --integrate` 可执行逐批分析和 protocol 级结果 Meta

**代谢组学：先校正再合并分析**

```
推荐流程：
  1. 两个 batch 分别提取峰面积（原始数据）
  2. 合并峰面积表 → 同一矩阵
  3. QC-RLSC 校正仪器漂移（利用 Pooled QC）
  4. ComBat 批次校正（利用 bridge samples 或 batch label）
  5. 在校正后数据上做 log2 变换 + 缺失值填补
  6. 统一做 t-test / PLS-DA / 差异分析

  关键：必须保留原始峰面积，不能只拿已处理的归一化数据
```

**蛋白质组学（DIA）：取决于是否用同一谱图库**

```
情况1 — 同一谱图库：
  1. 合并蛋白丰度矩阵
  2. log2 变换 + 中位数中心化
  3. ComBat 或 limma removeBatchEffect
  4. 统一做 limma 差异分析

情况2 — 不同谱图库：
  1. 先用共同鉴定蛋白做交集
  2. 对交集蛋白做分位数归一化对齐
  3. 再做 ComBat 校正
  4. 统一差异分析
```

**MFA：合并 MID 后统一拟合**

```
推荐流程：
  1. 检查新旧 batch 的 MID 分布一致性
  2. 如一致 → 合并所有样本的 MID 数据
  3. 用统一的 network model + bounds 重新拟合全部样本
  4. 在生物学重复水平做差异分析（代码已实现；历史产物需重新运行）

  关键：统一拟合参数比事后校正通量更可靠
```

#### 合并分析前的必做检查

```
Step 1: 各组学分别做 PCA，标注 batch 信息
  → 检查 batch 效应是否是主要变异来源
  → 理想：condition 主导，batch 内混杂
  → 警告：batch 形成分开的 cluster

Step 2: Bridge samples 跨批次相关性
  → 同一样本在新旧 batch 的 Pearson r 应 > 0.9（RNA-Seq/蛋白）
  → 低于 0.7 说明批次效应过大，需加强校正

Step 3: 校正后重新 PCA 验证
  → 确认 batch cluster 消失，condition 分离保持
```

**警告**：千万不要在没有 batch 校正的情况下直接合并原始数据跑差异分析，否则差异结果反映的是批次差异而非生物学差异。

---

## 各组学差异分析统计方法评估

### 统计方法总览

| 组学 | 主要统计方法 | FDR 校正 | 归一化 | 评价 |
|---|---|---|---|---|
| **转录组** | PyDESeq2 (NB GLM + Wald；可选 `~batch + condition`) | BH (内置) | median-of-ratios | 联合模型与批次结果 Meta 均已实现；混淆设计拒绝推断 |
| **代谢组** | Student's t-test + PLS-DA VIP | BH (自实现) | log2 + 半最小值填补 | 基本符合，通路收敛分析已使用 BH FDR 校正 |
| **蛋白质组** | limma 经验贝叶斯 moderated t-test | BH (statsmodels) | log2 + 中位数中心化 + MNAR 下移正态填补 | 完全符合主流标准 |
| **MFA** | 每动物通量估计 + Welch/Mann-Whitney + Hedges' g | BH (statsmodels) | 样本内 MCMC 汇总 | 生物学重复推断已实现；n=1 仅描述；历史结果需重跑 |

### MFA 伪重复问题与当前修正

历史 MFA 产物曾在 MCMC 拟合点（~304 vs ~228）上运行检验，而非在生物学重复均值（n=4 vs n=3）上。同一动物的 MCMC 点不是独立生物学观测，历史 p 值可能被严重压缩。

**当前代码**：先把每只动物的拟合分布汇总为一个通量估计，再进行组间差异推断；输出 log2FC SE、Hedges' g、`inference_unit=biological_replicate` 以及可选 bootstrap/置换稳定性。n=1 只允许描述性输出。旧矩阵仍保留历史兼容字段，但不能用来重新解释尚未重跑的旧 CSV。

所有可视化模块（仪表板、火山图、差异汇总）也已统一使用相同的阈值标准，通过 `_config.py` 从 `config.yaml` 配置驱动，消除了阈值不一致问题。

**迁移要求**：重新运行 MFA 后，确认差异表具有 `log2fc_se`、`hedges_g`、`n_eff` 和 `inference_unit=biological_replicate`，再把对应多组学配置的 `evidence_profile.mfa_inference_unit` 改为 `biological_replicate`。未重跑的历史文件继续按 `mcmc_fit_points` 降级，避免伪造升级。

---

## 参考文献

1. Escofier, B., & Pagès, J. (1994). Multiple factor analysis (AFMULT package). *Computational Statistics & Data Analysis*, 18(1), 121-140.
2. Pagès, J. (2014). *Multiple factor analysis by example using R*. CRC press.
3. Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto encyclopedia of genes and genomes. *Nucleic Acids Research*, 28(1), 27-30.
4. Subramanian, A., et al. (2005). Gene set enrichment analysis. *PNAS*, 102(43), 15545-15550.
5. Antoniewicz, M. R. (2015). Methods and advances in metabolic flux analysis: a mini-review. *Journal of Industrial Microbiology & Biotechnology*, 42(3), 317-325.

---

*最后更新：2026年8月13日*
