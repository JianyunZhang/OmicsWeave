# OmicsWeave 通用多组学分析与评估能力：从 DIO 案例提炼的改进规划

> 文档类型：调研 + 设计建议
> 版本：v0.7（2026-11-24，增加 11.10：预排序 GSEA 双层开关与 `run_gsea: false` 默认关闭优化）；v0.6（2026-11-24，增加 11.9：G1 / G2 / G18 转录组样本质控双重惩罚修正、分级风险标注、统计保留与通路图去感叹号降噪）；v0.5（2026-11-23，增加 11.8：G1 通用转录组 / MFA 样本质检、边界贴边强化与评估规则硬化，更新 11.6 状态）；v0.4（2026-11-18，增加 11.7：G10 中不依赖新数据的部分；11.6 记下组学结果包早于批次 2）；v0.3（2026-11-17，增加 11.6 未完成事项，删除文末误贴的对话文字；同日随反应档案 S1–S5 更新 11.6 与第 12 节；v0.2 为 2026-11-15，增加第 11 节实施记录；v0.1 为 2026-11-14）
> 依据：两个 DIO 队列在 2026-11-03 至 11-13 的分析与映射修复（阶段 0–3，提交 `75f0bb8` 至 `d72b2c7`），以及课题工作区中的评估记录
> （`examples/dio_project/analysis/focus_targets/INTEGRATION_ASSESSMENT.md`、`midterm_report/midterm_key_pathways.md`，未纳入版本库；
> 2026-11-18 起两个课题目录在 `local/projects/` 下）
> 配套文档：[架构规划 v2](../omicsweave_arch_design_v2.md)、[OmicsWeave × TheraWeave 开发规划](../omicsweave_theraweave_development_plan.md)。
> 本文的建议按 v2 第 11 节的阶段编号归位（第 7 节），不改变已有决策
> 范围：只讨论流程能力与评估方法。文中的数字描述的是流程行为（格子数、匹配率、入选数），不含课题的生物学结论

---

## 1. 结论摘要

1. **案例中真正产出结论的步骤，多数写在课题脚本里，不在通用流程里。** 跨队列同向筛选、多臂合并模型、类别比值、逆转判定、
   重跑前后逐格对照，都是临时脚本。通用化的主线是把这些步骤做成流程的一部分（第 3 节）。
2. **通用流程目前擅长“算出统计量”，不擅长“判断统计量能不能信”。** 样本组成问题、映射错误、锚点证据不足，都是靠人工发现的。
   建议把质量信息做成一等数据，自动传到整合层（G1、G2、G5）。
3. **评估质量需要可度量的指标。** 目前只有回归测试能说明“结果没变”，没有指标能说明“证据等级可信”。
   建议每次整合运行输出一张评估卡，含 7 项指标（第 6 节）。
4. **建议实现 12 项能力，分三批**（细节见第 5、7 节）：

| 批次 | 能力 | 目的 | 估计工作量 |
|:--:|---|---|---|
| 1 | G5 证据等级与锚点状态挂钩；G2 质量信息自动进入可靠性；G4 逆转判定；G11 结果包对照与查询；E7 配置项消费者测试 | 让现有结果不再被误读 | 约 2.5–3 周 |
| 2 | G1 样本身份与组成 QC；G3 设计感知的统计；G6 复现模式；G8 派生特征规范；G7 知识层（代谢物注册表、同工酶加权）；E3、E5 | 把课题脚本里的有效分析产品化 | 约 7–8 周 |
| 3 | G9 通路层整合升级；G10 示踪层的不依赖模型读数；G12 结论台账；E1、E2 | 提高发现能力，建立基准 | 约 7–8 周 |

工作量按 1 名开发者估算，仅供排期参考。

---

## 2. 依据：案例暴露的通用问题

| # | 现象（流程行为） | 根因 | 通用教训 | 对应能力 |
|---|---|---|---|---|
| 1 | 手写的蛋白注册表与数据自带登记号一致的只有 93 / 185；重建后 99.1% | 知识表手写，从不对照数据 | 所有知识表都应由带版本的来源生成，并在每个项目上对照数据审计 | G7 |
| 2 | 修复前整合中唯一的 Convergent 来自误配；修复后有方向证据的格子由 87 降到 31 | 映射错误会制造证据，而不只是丢失证据 | 需要一个能估计“假证据率”的指标 | E3 |
| 3 | 31 个有等级的格子，锚点（MFA）层全部是“证据不足”；其中 14 个来自每组 2 对 1 的比较 | 等级的含义与锚点状态脱钩：锚点方向未确定时，其他层仍按它判支持或冲突 | 等级词汇应取决于锚点状态；样本量不足的比较不可评估 | G5 |
| 4 | 文库组成异常、离群对照、跨组织混入都由人工发现；31 个有等级的格子中 16 个因此不能采信 | QC 只检查读段质量与比对率；QC 结论不传到整合 | 组成与身份 QC 要通用化；与比较混杂的质量问题要自动降权 | G1、G2 |
| 5 | 秩检验在每组 3 只时达不到 0.05；“仅效应量”规则入选 50–58% 的通量 | 小样本下逐特征检验没有区分力 | 用方差借用、合并模型和可检出效应表，替代单一阈值 | G3 |
| 6 | 最有产出的分析是跨队列同向筛选和多臂合并模型，都在课题脚本里 | 流程只支持“一个因素、两个水平、一个子集”的比较 | 复现与多因素模型应是流程能力 | G3、G6 |
| 7 | `reversal_pairs` 已在配置模型中声明，但没有分析代码读取；`design.covariates` 只有两个组学使用 | 声明先于实现 | 每个配置项都要有消费者，并有测试保证 | G4、E7 |
| 8 | 类别总量与质量作用比加入后，新增的可复现信号主要来自“同一样本内的比值” | 派生特征只有两类，且定义写在代码里 | 派生特征应可声明、可由知识自动生成 | G8 |
| 9 | 方向成员引用了该组织中丰度很低的同工酶，而主要同工酶方向相反 | 映射不区分组织 | 成员按该组织的实测丰度加权 | G7 |
| 10 | 一个队列的通路融合产出为 0；富集在部分组织没有任何通路通过校正 | 通路按名称匹配；只有预排序 GSEA 一种检验；没有产出诊断 | 通路层需要适合小样本的检验和产出报告 | G9 |
| 11 | “组学差异表重跑前后逐值相同”与“等级如何变化”是用临时脚本证明的 | 没有结果包对照工具 | 对照应是一条命令 | G11 |
| 12 | 一次重跑使 6 处已写入报告的说法需要调整，靠人工逐条核对 | 结论与结果表之间没有可检查的对应关系 | 结论应登记为可执行的判据 | G12 |

---

## 3. 现状盘点：哪些分析还在课题脚本里

| 分析步骤 | 目前的位置 | 通用化后的归属 |
|---|---|---|
| 样本组成 QC（线粒体读数占比、留一相关、标志物） | 课题导出脚本，问题样本靠手写备注 | `qc` 步骤的通用模块（G1） |
| QC 结论 → 整合的降权 | 人工填写 `quality_reliability`（按组学层，不分组织与比较） | 自动生成的可靠性表（G2） |
| 多臂合并模型（`y ~ 因素 A + 因素 B`） | 课题脚本中的最小二乘 | `design.models`（G3） |
| 逆转判定 | 课题脚本按“方向相反且 q < 0.1”筛选 | `reversal_pairs` 的消费者（G4） |
| 跨队列同向筛选 | 课题脚本，规则手写 | `replication` 整合模式（G6） |
| 类别比值（产物 / 底物、亚类占比） | 课题脚本 | 派生特征规范（G8） |
| 映射与注册表审计 | 课题审计脚本；包内只有注册表对照函数 | `omicsweave audit`（G7、G11） |
| 重跑前后逐格对照 | 课题对照脚本 | `omicsweave diff`（G11） |
| 读取结果包 | 课题脚本中写死运行 ID 的读取函数 | 查询接口（G11） |
| 预期结果与结论口径 | Markdown 文档，人工核对 | 结论台账（G12） |

包内已经具备、可以直接复用的部件：随机效应 Meta（`core/statistics/meta_analysis.py`）、生物学重复的自助与置换校准（`stability.py`，目前只在 MFA 差异分析中使用）、
秩检验可达最小 p（`nonparametric.py`）、蛋白组的方差压缩线性模型、网络语义表与派生特征（阶段 1–3）。

---

## 4. 方法学调研

文献的题目、期刊与年份均在 PubMed 核对，编号对应第 10 节。

### 4.1 小样本推断

- **方差借用。** limma 的经验贝叶斯把各特征的方差向共同先验收缩，专为小样本设计 [1]；voom 把它扩展到计数数据 [2]；DESeq2 对离散度和倍数变化做收缩 [3]。
  本流程只有蛋白组和转录组用了这类方法，代谢组、派生特征和通量仍是逐特征的 Welch 检验或秩检验。
- **符号的确定性。** 自适应收缩由效应量和标准误估计局部假符号率（lfsr），比 p 值更直接地回答“方向可信吗” [4]；
  多条件扩展可在组织或比较之间共享信息 [5]。证据画像的设计本来就优先使用 `1 − lfsr`，但目前没有任何层提供这一列。
- **功效与夸大。** 低功效研究中，显著结果为真的比例低，效应量被系统性夸大 [6]；符号错误率与夸大倍数可以直接估计 [7]。
  多组学研究的功效应按各组学的变异分别计算 [8]。
- **推断单位。** 实验单位的定义决定 N [9]。本流程已在契约中记录 `inference_unit`，这一点应保持为硬约束。
- **加权的多重检验。** 用与检验统计量独立的协变量（丰度、可辨识性）给假设加权，可在控制 FDR 的同时提高功效 [10]。

### 4.2 逆转与“未改变”

- “药物不逆转”不能由“不显著”推出，需要等效性检验：效应落在预设界限内才算未改变 [11]。
- 用疾病签名与药物签名的反向关系寻找候选，是 Connectivity Map 的基本思想 [12]，与总纲中 TheraWeave 的用途一致。

### 4.3 复现

- 部分合取检验回答“效应是否至少在 u 个研究中存在”，并可做多重校正 [13]；r 值给出复现性声明的 FDR [14]。
  这两种方法正好对应案例中手写的“两队列同向且各自 p < 阈值”规则，但能控制错误率。

### 4.4 质量控制与组成

- 测序与质谱各有成熟的 QC 汇总工具 [15、16]；非靶向代谢组的混合 QC、系统适用性样本与基于 QC 的信号校正有明确指南 [17、18]。
- 批次效应在高通量数据中普遍存在，常与研究变量混杂 [19]；可以用对照基因或对照样本估计并去除未知的技术因素 [20]。
- mRNA 与蛋白的对应受翻译和降解调控，两者并不总是同步变化 [21]。因此 RNA 与蛋白的相关适合作为诊断指标，不适合作为硬性门槛。
- 组织块的细胞组成变化会表现为大量“差异基因”；可以用标志物或解卷积估计组成 [22]。

### 4.5 通路层整合

- ActivePathways 先合并各数据集的 p 值、再做排序富集，并指出每条通路由哪些数据集支持 [23]；其后续方法在合并时加入方向约束，
  方向不一致的基因被降权 [24]。multiGSEA 在各层分别做 GSEA 后合并 [25]。
- 竞争性检验需要考虑基因间相关，否则假阳性偏高 [26]；旋转检验适用于小样本和复杂设计 [27]。
- 脂质组有专门的本体富集工具 [28] 和基于脂质反应网络的通路分析 [29]。
- 基于代谢网络拓扑的 reporter metabolite 方法可以在没有通量数据时定位调控集中的代谢物 [30]；小鼠等模式动物已有基因组尺度代谢模型 [31]。

### 4.6 特征命名与派生特征

- 代谢物名称应归一到参考命名体系 [32]；脂质应按 LIPID MAPS 的简写规范解析 [33]，已有通用的语法解析器 [34]。
- 代谢物鉴定的可信等级有社区标准 [35]，应随特征进入结果表。

### 4.7 示踪层

- 天然丰度校正有成熟工具 [36、37]；通量的置信区间应由非线性方法得到，由局部标准差近似的区间并不可靠 [38]。
- 平行标记实验可以对同一网络联合拟合 [39]；循环营养物的周转与直接贡献可以不依赖组织网络模型求得 [40、41]。

### 4.8 评估与基准

- 目标–诱饵策略用已知为假的条目估计假发现率 [42]，可以移植到“反应 → 特征”的映射上。
- 在真实数据上加入已知效应的半合成模拟，能在保留真实噪声结构的前提下评估灵敏度 [43]。
- 多组学联合降维方法的基准比较显示，不同方法在不同任务上各有优劣 [44]；领域综述把基准比较、FAIR 与可重复性列为需要遵循的实践 [45]。
- 因子模型与监督式整合 [46、47] 以同一批样本的多组学矩阵为输入；样本不配对或样本量很小时不适用（这是本文的判断）。

### 4.9 报告规范

- 动物实验报告应包含样本量依据、实验单位、纳入与剔除标准 [48]；数据与流程应满足 FAIR 原则 [49]。

---

## 5. 建议的通用能力

每项按“问题 → 设计 → 接口与契约 → 验收”说明。验收都以合成示例项目（CI 中运行）为准，课题数据只作本地回归。

### G1 样本身份与组成 QC（跨组学）

- **问题**：现有 QC 门控只检查 Q20 / Q30、N 含量、比对率与读数（`transcriptomics/qc_gate.py`）。组成异常、跨组织混入、离群样本没有自动检查。
- **设计**：新增通用模块，输入为分析层数据集（特征 × 样本矩阵），对任何组学计算：
  ① 样本与同组织其他样本的留一相关及其稳健 z 值；② 细胞器或区室占比（如线粒体编码基因的读数占比）；
  ③ 组织身份得分与外源组织得分，依据随物种分发的组织标志物表；④ 性别标志；
  ⑤ 同一动物在同一组织有多个组学时的跨组学身份核对。阈值与处理策略沿用现有的 `flag / exclude / fail`。
- **接口与契约**：`omics/<组学>/qc/sample_quality.parquet`（`sample_id`、`metric`、`value`、`robust_z`、`flag`、`reason`）；标志物表由 G7 的注册表框架生成。
- **验收**：示例项目中人为混入另一组织信号的样本被标记；干净样本的误报率低于设定值。

### G2 质量信息自动进入可靠性

- **问题**：整合层只有人工填写的 `quality_reliability`，粒度是组学层。
- **设计**：由 G1 的样本标记和各组学 QC 表生成“组学 × 组织 × 比较”的可靠性。规则区分两种情况：标记样本在两组间均衡时按可用样本比例降权；
  标记样本集中在一组（与比较混杂）时，该比较的方向证据记为不可评估。MFA 层同样读取可辨识性与示踪剂一致性。
- **接口与契约**：`qc/contrast_reliability.parquet`（`omics`、`tissue`、`contrast_id`、`n_usable_numerator`、`n_usable_denominator`、`reliability`、`reason`）；
  差异结果表的 `flags` 增加 `qc_confounded`；`flux_anchored` 读取该表，人工系数只作为覆盖项保留。
- **验收**：把一组样本全部标记后，该比较在整合中不再产生支持或冲突；运行清单登记可靠性表为输入。

### G3 设计感知的统计

- **问题**：比较只能是“一个因素的两个水平加一个子集”。多臂设计中信息最多的估计（跨处理臂的主效应、交互项）做不了；连续型数据层没有方差借用。
- **设计**：
  1. `design.models`：在已声明的因素上写公式（如 `~ diet + treatment`），指定要报告的项。对连续型数据层用同一个线性模型引擎，计数数据仍用 PyDESeq2。
     合并估计与两两比较并列输出，不互相替代；文档中注明被处理逆转的指标合并估计会偏小。
  2. 把蛋白组已有的方差压缩线性模型提到 `core/statistics`，供代谢组、派生特征和通量使用。
  3. 可选的自适应收缩，输出 `lfsr`，供证据画像使用。
  4. 可检出效应表：每个“组学 × 组织 × 比较”给出样本量、残差标准差中位数、80% 功效下的最小可检出效应、秩检验可达最小 p。
  5. `design.covariates` 对四个组学都生效。
- **接口与契约**：差异结果表增加 `estimand_type`（`pairwise`、`pooled_main`、`interaction`）、`n_total`、`lfsr`（可空）；新增 `qc/design_power.parquet`。
- **验收**：在两因素的合成数据上，合并主效应的估计与真值一致，功效高于单臂比较；设计矩阵不满秩时拒绝拟合并说明原因。

### G4 逆转与“未改变”的判定

- **问题**：`reversal_pairs` 没有消费者；课题脚本用“方向相反且过阈值”判定，“不逆转”则由“不显著”推出。
- **设计**：对每个逆转配对、每个特征（含派生特征与通路）计算疾病效应、处理效应与逆转指数（处理效应与疾病效应之比的相反数）及其区间。分为五类：
  已逆转、部分逆转、未逆转（需通过等效性检验 [11]）、过度逆转、无法判定。疾病效应本身未建立时不做判定。设计为析因时同时给出交互项。
- **接口与契约**：`tables/reversal.parquet`（`pair_id`、`feature_id`、`disease_effect`、`treatment_effect`、`reversal_index`、`ci_lower`、`ci_upper`、`class`、`rule`）。
- **验收**：合成数据中五类情形各自被正确归类；“未逆转”不会由低功效的不显著结果产生。

### G5 证据等级与锚点状态挂钩

- **问题**：锚点层证据不足时，其他层仍按锚点的点估计方向判支持或冲突（第 2 节第 3 条）。
- **设计**：
  1. 锚点门控：锚点层的 `evidence_level` 为 `insufficient`，或任一组少于 2 只时，不产生 Convergent / Single-layer / Conflicting。
  2. 无锚点时改为评估各组学层之间的相互一致性，使用单独的等级词汇（如“组学一致”“组学单层”“组学不一致”），与有锚点的等级不混用。
  3. 锚点类型泛化：除通量外，可以是预先登记的预期方向（G12），或指定的参考层。没有 MFA 的项目因此也能使用证据画像。
  4. 把“已复现”（G6）作为独立维度写入画像，而不是并入支持度。
- **接口与契约**：证据画像增加 `anchor_type`、`anchor_state`；等级枚举扩展；旧等级的含义不变。
- **验收**：锚点证据不足的格子不再出现有方向含义的等级；把锚点的符号翻转后，无锚点等级不变（这是一条可自动运行的不变性测试）。

### G6 复现模式

- **问题**：跨队列复现是案例中最有效的筛选，但规则手写、不控制错误率，也不区分同层复现与跨层复现。
- **设计**：项目级整合模式 `replication`。在 `project.yaml` 中声明可比的实验与比较（估计对象相同），可选声明允许的跨层对应（如一个队列的蛋白对另一个队列的 RNA）。
  对每个特征、派生特征、通路计算符号一致性、随机效应合并（已有）、部分合取 p 值与校正 [13、14]。同层复现与跨层复现分别标记。
  同时报告随机特征的符号一致率作为基线。
- **接口与契约**：`integration/replication/tables/replication.parquet`；证据画像读取其中的 `replicated` 与 `replication_fdr`。
- **验收**：合成的两个实验中，共享效应被召回，单实验特有效应不被标为复现；基线一致率接近 50%。

### G7 知识层：注册表框架、网络语义、同工酶加权

- **问题**：基因与蛋白注册表已改为生成式，代谢物 ID 表仍是手写且有错配；审计只在测试与课题脚本中运行；方向成员不区分组织。
- **设计**：
  1. 注册表框架：每张知识表由“生成器 + 来源清单 + 对照数据的审计函数”组成，统一入口 `omicsweave refs build <名称>`。
     `validate` 对当前项目运行审计，报告覆盖率与冲突。待建：代谢物（按参考命名 [32]）、脂质简写解析 [33、34]、组织标志物（G1）、同源映射（v2 的 2.4）。
  2. 网络语义：阶段 1 已能从一个 MFA 网络生成反应注释。把它抽象为“模型包”接口，使其他网络和基因组尺度模型的基因–反应规则 [31] 也能提供反应语义。
  3. 同工酶加权：方向成员按该组织数据集中的平均丰度占比加权，输出容量变化的加权值；支持证据来自占比低于阈值的成员且主要成员方向相反时，标记 `minor_isozyme_driven`。
- **验收**：示例项目中登记号与数据不一致的条目被 `validate` 报出；主要与次要同工酶方向相反的合成反应不被判为支持。

### G8 派生特征规范

- **问题**：派生特征只有类别总量与固定的几个质量作用比，定义在代码中。
- **设计**：声明式规范，类型包括 `sum`（已有）、`ratio`（已有）、`step_ratio`（一步反应的产物对底物，可由网络化学计量或脂质组装步骤表自动生成）、
  `fraction`（亚类占比）、`index`（如按丰度加权的链长与不饱和度）、`module_score`（基因或蛋白集合的得分）。
  每种类型在“组织 × 比较”内单独校正（现有做法）。规范对任何组学的数据集适用，产出进入差异、逆转、复现和整合。
- **接口与契约**：`derived/features.parquet` 增加 `feature_type`、`definition`、`generated_from`；鉴定可信等级 [35] 随成员特征带入。
- **验收**：由网络自动生成的 `step_ratio` 与手写定义逐值一致；成员缺失时该特征不计算并说明原因。

### G9 通路层整合升级

- **问题**：通路按名称匹配；检验只有预排序 GSEA；各层合并时不考虑方向；没有产出诊断。
- **设计**：
  1. 匹配改用条目 ID 与集合重叠。
  2. 对连续型数据层增加适合小样本的集合检验（旋转检验与考虑相关的竞争性检验 [26、27]），与 GSEA 并列输出。
  3. 跨层合并采用带方向约束的 p 值合并 [23、24]，替代对校正后 p 值做 Fisher 合并。
  4. 脂质类别与脂质本体 [28] 作为集合库。
  5. 产出报告：每个“组学 × 组织 × 比较 × 库”的检验条目数、最小 p、通过数；产出为 0 时给出提示。
- **验收**：合成数据中方向一致的多层通路排在方向相反的通路之前；产出报告写入结果包。

### G10 示踪层的不依赖模型读数

- **问题**：示踪数据只以求解后的通量进入流程；通量可辨识性差时没有替代读数。
- **设计**：同位素分布（MID）的数据契约（v2 的 3.4）；血浆富集平台检查；按前体富集归一化的标记比作为派生特征进入同一套差异引擎；
  多示踪剂的直接贡献矩阵 [40、41]；化学计量上锁死的反应在检验前合并为通量组；通量效应量改用相对参考通量的差值，区间取自求解器的非线性置信区间 [38]。
- **验收**：合成的两示踪剂数据中，贡献矩阵还原设定的贡献比例；通量组内的反应只计一次检验。

### G11 结果包工具：对照、查询、审计

- **设计**：
  1. `omicsweave diff <运行 A> <运行 B>`：逐表判断是否逐值相同，给出特征级差值分布与证据等级的转移表。
  2. 查询接口：按实验与运行读取统一结果表、派生特征与数据集，支持“最新运行”与显式指定两种方式，返回值带运行 ID。课题脚本不再写死路径。
  3. `omicsweave audit`：注册表对照数据、映射覆盖率、没有消费者的配置项。
- **验收**：对同一输入的两次运行，`diff` 报告全部相同；修改一个参数后，`diff` 只报告受影响的表。

### G12 结论台账

- **问题**：结论写在文档里，重跑后靠人工核对口径（第 2 节第 12 条）。
- **设计**：课题目录中的 `claims.yaml`（不进入包）。每条结论有编号、表述、级别（结论或假说）和可执行的判据：数据层、组织、比较、特征或通路、预期方向、所需证据（如“两队列复现”）。
  `omicsweave claims check` 对一次运行给出每条结论的状态（支持、减弱、相反、不可检验），并与上一次运行比较，列出状态变化。包内只有判据引擎，不含任何课题内容。
- **验收**：示例项目的台账在 CI 中检查；判据引用的特征不存在时报错而不是静默通过。

---

## 6. 多组学评估质量的度量

目标是让每次整合运行都回答：这些证据等级有多可信。建议在结果包中新增 `evaluation/`，输出一张评估卡。

| # | 指标 | 定义 | 用途 | 实现要点 |
|---|---|---|---|---|
| E1 | 校准 | 按动物置换组别标签后重跑差异与整合，统计各等级的期望格子数 | 回答“这么多 Convergent 里有多少是偶然的” | 每组 3 只时不同的置换只有 10 种，需注明分辨率；置换在动物层面进行 |
| E2 | 灵敏度 | 在真实数据上加入已知大小与方向的效应 [43]，统计被召回的比例；另支持用户声明的阳性对照 | 回答“给定样本量下流程能看到多大的效应” | 加标发生在分析层数据集，流程其余部分不变 |
| E3 | 假证据率 | 给每个反应配一组丰度匹配的随机特征作为诱饵映射 [42]，统计诱饵获得支持的比例 | 度量映射层与判定规则制造证据的倾向 | 诱饵集在运行清单中登记随机种子 |
| E4 | 复现率 | 入选特征在另一队列中符号一致的比例，与随机特征的基线比较 | 外部有效性 | 依赖 G6 |
| E5 | 稳健性 | 逐只去掉动物后入选与等级是否保持；估计符号错误率与夸大倍数 [7] | 标出由单只动物驱动的结果 | 输出 `single_subject_driven` 标记 |
| E6 | 层间一致性 | 同一“组织 × 比较”中 RNA 与蛋白倍数变化的相关；各层符号一致率 | 诊断某一层是否可疑 [21] | 只作诊断，不作门槛 |
| E7 | 回归与完整性 | 金标准结果包的逐值对照（G11）；每个配置项至少被一个分析读取的测试；映射与网络的同步测试（已有） | 防止退化 | 进入 CI |

评估卡给出各指标的数值与预设的参考范围，不合成单一总分，理由与 v2 第 13 节关于加和总分的决定相同。

---

## 7. 路线与优先级

| 批次 | 项目 | 依赖 | 归位（v2 第 11 节） | 估计工作量 |
|:--:|---|---|---|---|
| 1 | G5 锚点门控与无锚点等级 | 无 | 2.3、2.8 | 2–3 天 |
| 1 | G2 可靠性表（先接入已有的 QC 表与 MFA 可辨识性） | 无 | 2.3 | 3–4 天 |
| 1 | G4 逆转判定 | 无 | 2.5 的前置 | 3 天 |
| 1 | G11 `diff` 与查询接口 | 无 | 2.2 | 3–4 天 |
| 1 | E7 配置项消费者测试 | 无 | 总纲 7.4 | 1 天 |
| 2 | G1 样本身份与组成 QC | G7 的标志物表 | 2.1 | 1–1.5 周 |
| 2 | G3 公式模型、方差借用、可检出效应表 | 无 | 2.1 | 1.5–2 周 |
| 2 | G6 复现模式 | G11 的查询接口 | 2.3 | 1 周 |
| 2 | G8 派生特征规范 | 阶段 1–2 的网络语义与派生特征 | 2.1 | 1 周 |
| 2 | G7 代谢物注册表、脂质解析、同工酶加权 | 无 | 2.4 | 1–1.5 周 |
| 2 | E3 诱饵映射、E5 留一稳健性 | G5 | 2.8 | 1 周 |
| 3 | G9 通路层整合升级 | G8 | 2.3 | 2 周 |
| 3 | G10 示踪层读数 | MID 契约 | 3.4 | 2–3 周，取决于数据到位 |
| 3 | G12 结论台账 | G6、G11 | 2.5、2.6 | 1 周 |
| 3 | E1 校准、E2 半合成加标与公开数据基准 | G3 | 4 | 2 周 |

批次 1 不需要新数据，也不改变任何已有统计量，只改变等级与标记，建议先做。批次 1 完成后，现有两个课题的整合应重跑，并用 `diff` 留存对照。

---

## 8. 不建议现在做的

| 内容 | 原因 |
|---|---|
| 把因子模型或监督式整合作为默认 | 需要样本配对且样本量足够 [46、47]；目前的课题设计不满足，`sample_level` 模式保留为条件性补充 |
| 自己实现富集算法 | v2 的原则是外部工具走适配器；G9 只增加检验的种类与合并方式 |
| 把评估卡合成一个总分 | 与删除加和总分的决定冲突；各指标含义不同 |
| 重度的细胞组成解卷积 | 需要匹配的单细胞参考 [22]；先用标志物得分做诊断即可 |
| 约束模型的通量预测 | 属于中期选项；G7 只需要基因–反应规则，不需要求解 |

---

## 9. 风险

| 风险 | 对策 |
|---|---|
| 能力增加后配置变复杂 | 新配置项都有默认值；`validate` 说明当前设计支持哪些分析、缺什么 |
| 小样本下新统计方法给出过于乐观的结果 | E1、E2 先行；每种新方法进入默认前，先在半合成数据上与现有方法对照 |
| 标志物表、诱饵集带来新的知识维护负担 | 一律走 G7 的注册表框架，由来源生成并登记版本 |
| 课题内容混入通用包 | 结论台账、预期结果、阳性对照清单都放在课题目录；包内只有引擎与合成示例 |
| 等级词汇变化影响下游（TheraWeave） | 旧等级含义不变，只新增；交换契约按主版本管理 |

---

## 10. 参考文献

以下文献的题目、期刊与年份已在 PubMed 核对。

1. Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 2015. [doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)
2. Law CW, et al. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol* 2014. [doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)
3. Love MI, et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 2014. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
4. Stephens M. False discovery rates: a new deal. *Biostatistics* 2017. [doi:10.1093/biostatistics/kxw041](https://doi.org/10.1093/biostatistics/kxw041)
5. Urbut SM, et al. Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions. *Nat Genet* 2019. [doi:10.1038/s41588-018-0268-8](https://doi.org/10.1038/s41588-018-0268-8)
6. Button KS, et al. Power failure: why small sample size undermines the reliability of neuroscience. *Nat Rev Neurosci* 2013. [doi:10.1038/nrn3475](https://doi.org/10.1038/nrn3475)
7. Gelman A, Carlin J. Beyond power calculations: assessing Type S (sign) and Type M (magnitude) errors. *Perspect Psychol Sci* 2014. [doi:10.1177/1745691614551642](https://doi.org/10.1177/1745691614551642)
8. Tarazona S, et al. Harmonization of quality metrics and power calculation in multi-omic studies. *Nat Commun* 2020. [doi:10.1038/s41467-020-16937-8](https://doi.org/10.1038/s41467-020-16937-8)
9. Lazic SE, et al. What exactly is 'N' in cell culture and animal experiments? *PLoS Biol* 2018. [doi:10.1371/journal.pbio.2005282](https://doi.org/10.1371/journal.pbio.2005282)
10. Ignatiadis N, et al. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods* 2016. [doi:10.1038/nmeth.3885](https://doi.org/10.1038/nmeth.3885)
11. Schuirmann DJ. A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. *J Pharmacokinet Biopharm* 1987. [doi:10.1007/BF01068419](https://doi.org/10.1007/BF01068419)
12. Lamb J, et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science* 2006. [doi:10.1126/science.1132939](https://doi.org/10.1126/science.1132939)
13. Benjamini Y, Heller R. Screening for partial conjunction hypotheses. *Biometrics* 2008. [doi:10.1111/j.1541-0420.2007.00984.x](https://doi.org/10.1111/j.1541-0420.2007.00984.x)
14. Heller R, et al. Deciding whether follow-up studies have replicated findings in a preliminary large-scale omics study. *Proc Natl Acad Sci U S A* 2014. [doi:10.1073/pnas.1314814111](https://doi.org/10.1073/pnas.1314814111)
15. Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 2016. [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)
16. Bielow C, et al. Proteomics Quality Control: quality control software for MaxQuant results. *J Proteome Res* 2016. [doi:10.1021/acs.jproteome.5b00780](https://doi.org/10.1021/acs.jproteome.5b00780)
17. Broadhurst D, et al. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics* 2018. [doi:10.1007/s11306-018-1367-3](https://doi.org/10.1007/s11306-018-1367-3)
18. Dunn WB, et al. Procedures for large-scale metabolic profiling of serum and plasma using gas chromatography and liquid chromatography coupled to mass spectrometry. *Nat Protoc* 2011. [doi:10.1038/nprot.2011.335](https://doi.org/10.1038/nprot.2011.335)
19. Leek JT, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nat Rev Genet* 2010. [doi:10.1038/nrg2825](https://doi.org/10.1038/nrg2825)
20. Risso D, et al. Normalization of RNA-seq data using factor analysis of control genes or samples. *Nat Biotechnol* 2014. [doi:10.1038/nbt.2931](https://doi.org/10.1038/nbt.2931)
21. Buccitelli C, Selbach M. mRNAs, proteins and the emerging principles of gene expression control. *Nat Rev Genet* 2020. [doi:10.1038/s41576-020-0258-4](https://doi.org/10.1038/s41576-020-0258-4)
22. Newman AM, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nat Biotechnol* 2019. [doi:10.1038/s41587-019-0114-2](https://doi.org/10.1038/s41587-019-0114-2)
23. Paczkowska M, et al. Integrative pathway enrichment analysis of multivariate omics data. *Nat Commun* 2020. [doi:10.1038/s41467-019-13983-9](https://doi.org/10.1038/s41467-019-13983-9)
24. Slobodyanyuk M, et al. Directional integration and pathway enrichment analysis for multi-omics data. *Nat Commun* 2024. [doi:10.1038/s41467-024-49986-4](https://doi.org/10.1038/s41467-024-49986-4)
25. Canzler S, Hackermüller J. multiGSEA: a GSEA-based pathway enrichment analysis for multi-omics data. *BMC Bioinformatics* 2020. [doi:10.1186/s12859-020-03910-x](https://doi.org/10.1186/s12859-020-03910-x)
26. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 2012. [doi:10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461)
27. Wu D, et al. ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics* 2010. [doi:10.1093/bioinformatics/btq401](https://doi.org/10.1093/bioinformatics/btq401)
28. Molenaar MR, et al. LION/web: a web-based ontology enrichment tool for lipidomic data analysis. *Gigascience* 2019. [doi:10.1093/gigascience/giz061](https://doi.org/10.1093/gigascience/giz061)
29. Gaud C, et al. BioPAN: a web-based tool to explore mammalian lipidome metabolic pathways on LIPID MAPS. *F1000Res* 2021. [doi:10.12688/f1000research.28022.2](https://doi.org/10.12688/f1000research.28022.2)
30. Patil KR, Nielsen J. Uncovering transcriptional regulation of metabolism by using metabolic network topology. *Proc Natl Acad Sci U S A* 2005. [doi:10.1073/pnas.0406811102](https://doi.org/10.1073/pnas.0406811102)
31. Wang H, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. *Proc Natl Acad Sci U S A* 2021. [doi:10.1073/pnas.2102344118](https://doi.org/10.1073/pnas.2102344118)
32. Fahy E, Subramaniam S. RefMet: a reference nomenclature for metabolomics. *Nat Methods* 2020. [doi:10.1038/s41592-020-01009-y](https://doi.org/10.1038/s41592-020-01009-y)
33. Liebisch G, et al. Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived lipid structures. *J Lipid Res* 2020. [doi:10.1194/jlr.S120001025](https://doi.org/10.1194/jlr.S120001025)
34. Kopczynski D, et al. Goslin: a grammar of succinct lipid nomenclature. *Anal Chem* 2020. [doi:10.1021/acs.analchem.0c01690](https://doi.org/10.1021/acs.analchem.0c01690)
35. Sumner LW, et al. Proposed minimum reporting standards for chemical analysis: Chemical Analysis Working Group (CAWG) Metabolomics Standards Initiative (MSI). *Metabolomics* 2007. [doi:10.1007/s11306-007-0082-2](https://doi.org/10.1007/s11306-007-0082-2)
36. Millard P, et al. IsoCor: isotope correction for high-resolution MS labeling experiments. *Bioinformatics* 2019. [doi:10.1093/bioinformatics/btz209](https://doi.org/10.1093/bioinformatics/btz209)
37. Su X, et al. Metabolite spectral accuracy on Orbitraps. *Anal Chem* 2017. [doi:10.1021/acs.analchem.7b00396](https://doi.org/10.1021/acs.analchem.7b00396)
38. Antoniewicz MR, et al. Determination of confidence intervals of metabolic fluxes estimated from stable isotope measurements. *Metab Eng* 2006. [doi:10.1016/j.ymben.2006.01.004](https://doi.org/10.1016/j.ymben.2006.01.004)
39. Leighty RW, Antoniewicz MR. COMPLETE-MFA: complementary parallel labeling experiments technique for metabolic flux analysis. *Metab Eng* 2013. [doi:10.1016/j.ymben.2013.08.006](https://doi.org/10.1016/j.ymben.2013.08.006)
40. Hui S, et al. Glucose feeds the TCA cycle via circulating lactate. *Nature* 2017. [doi:10.1038/nature24057](https://doi.org/10.1038/nature24057)
41. Hui S, et al. Quantitative fluxomics of circulating metabolites. *Cell Metab* 2020. [doi:10.1016/j.cmet.2020.07.013](https://doi.org/10.1016/j.cmet.2020.07.013)
42. Elias JE, Gygi SP. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. *Nat Methods* 2007. [doi:10.1038/nmeth1019](https://doi.org/10.1038/nmeth1019)
43. Gerard D. Data-based RNA-seq simulations by binomial thinning. *BMC Bioinformatics* 2020. [doi:10.1186/s12859-020-3450-9](https://doi.org/10.1186/s12859-020-3450-9)
44. Cantini L, et al. Benchmarking joint multi-omics dimensionality reduction approaches for the study of cancer. *Nat Commun* 2021. [doi:10.1038/s41467-020-20430-7](https://doi.org/10.1038/s41467-020-20430-7)
45. Krassowski M, et al. State of the field in multi-omics research: from computational needs to data mining and sharing. *Front Genet* 2020. [doi:10.3389/fgene.2020.610798](https://doi.org/10.3389/fgene.2020.610798)
46. Argelaguet R, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biol* 2020. [doi:10.1186/s13059-020-02015-1](https://doi.org/10.1186/s13059-020-02015-1)
47. Singh A, et al. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics* 2019. [doi:10.1093/bioinformatics/bty1054](https://doi.org/10.1093/bioinformatics/bty1054)
48. Percie du Sert N, et al. The ARRIVE guidelines 2.0: updated guidelines for reporting animal research. *PLoS Biol* 2020. [doi:10.1371/journal.pbio.3000410](https://doi.org/10.1371/journal.pbio.3000410)
49. Wilkinson MD, et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data* 2016. [doi:10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18)

---

## 11. 实施记录

### 11.1 实施条件评估（2026-11-15）

| 批次 | 项目 | 依赖是否满足 | 结论 |
|:--:|---|---|---|
| 1 | G5、G2、G4、G11（diff、查询）、E7 | 无依赖，不需要新数据 | 已实施（11.2） |
| 1 | G11 第 3 项 `omicsweave audit` | 与 G7 的注册表审计重叠 | 随 G7 实施 |
| 2 | G3 公式模型、方差借用、可检出效应表 | 无依赖 | 具备条件；会改变统计量，按第 9 节应先在半合成数据上与现有方法对照，宜以可选方式接入 |
| 2 | G7 代谢物注册表、脂质解析、同工酶加权 | 无依赖 | 具备条件；需下载 RefMet 与 LIPID MAPS 等来源 |
| 2 | G8 派生特征规范 | 阶段 1–2 已完成 | 具备条件 |
| 2 | G6 复现模式 | G11 查询接口 | 本批完成后具备 |
| 2 | E3 诱饵映射、E5 留一稳健性 | G5 | 本批完成后具备 |
| 2 | G1 样本身份与组成 QC | G7 的标志物表 | 部分具备：留一相关、线粒体占比不需要标志物；组织身份与性别需先建标志物注册表 |
| 3 | G9 通路层整合升级 | G8 | 不具备 |
| 3 | G10 示踪层读数 | MID 契约与示踪批次数据导入 | 不具备（数据未到位） |
| 3 | G12 结论台账 | G6、G11 | 不具备（G6 未完成） |
| 3 | E1 校准、E2 半合成加标 | G3 | 不具备 |
| — | E4 复现率 | G6 | 随 G6；E6 层间一致性不在第 7 节路线中，无依赖，可随时补 |

### 11.2 批次 1 的实现

| 项目 | 实现 | 未做的部分（原因） |
|---|---|---|
| G5 | `evidence_profile.anchor_state`：锚点状态 established / effect_only / not_established / insufficient / not_identifiable / tracer_conflict / qc_confounded / qc_insufficient；锚点不成立时等级为 Omics-concordant / Omics-single-layer / Omics-discordant，由各层不依赖锚点符号的方向提示（`implied_direction`）决定；画像增加 `anchor_type`、`anchor_state`、`omics_direction` | 第 3 项锚点类型泛化（依赖 G12）；第 4 项“已复现”维度（依赖 G6） |
| G2 | `core/reliability.py`：样本表 `qc_flag` / `qc_reason` 与各组学 `qc/sample_quality.parquet`（G1 的契约）→ `contrast_reliability`（ok / downweighted / qc_confounded / qc_insufficient）；`omicsweave run` 写出 `qc/contrast_reliability.parquet` 并给差异结果加 `qc_confounded` 标记；整合按当前样本表重算；人工 `quality_reliability` 只作覆盖项；MFA 的可辨识性与示踪剂一致性进入锚点状态 | 由现有 QC 指标（文库大小、缺失率、拟合损失）自动生成样本标记属于 G1，需先确定误报率 |
| G4 | 整合模式 `reversal`：逆转指数与 delta 法区间，类别 reversed / partially_reversed / over_reversed / not_reversed（双单侧等效性检验）/ undetermined，另加 aggravated（处理效应显著且与疾病同向，第 5 节未列）；含代谢组派生特征；QC 不可评估的比较归入 undetermined | 通路层面（需要通路效应的标准误，G9）；析因交互项（需要 G3 的公式模型）；两个比较共用一组动物时的协方差（区间按独立计算） |
| G11 | `omicsweave diff`（逐表、逐格、差值分布、分类列转移，退出码 0 / 1 / 2）；`contract/query.py`（`open_run`、`read_table`，返回值带 `run_id`） | `omicsweave audit`（见 11.1） |
| E7 | `tests/core/test_config_consumers.py`：每个配置字段须登记读取它的模块、或说明为描述性、或列为已知缺口（缺口有探针，被实现后测试提示移出） | 发现的缺口：`defaults.padj_method`、`defaults.random_seed`、`factors.*.reference`、`design.stratify_by`；`design.covariates` 在代谢组与 MFA 中未接入 |

### 11.3 现有课题的重跑与对照

按第 7 节的要求，批次 1 完成后两个课题的整合已重跑，并用 `omicsweave diff` 与 2026-11-12 的整合结果对照
（源组学结果包不变：A `run_20261112_102147_5d76`，B `run_20261112_113117_1269`）：

| | 队列 A | 队列 B |
|---|---|---|
| 新整合结果包 | `run_20261115_010915_d97b` | `run_20261115_011321_316d`（另含 `reversal`） |
| 逐值相同的表 | 通路收敛 3 张、元数据 2 张 | 同左 |
| 证据长表 | 只新增列，原有值不变 | 同左 |
| 锚点状态 | not_identifiable 807、not_established 624、effect_only 162、tracer_conflict 10、insufficient 2；established 0 | not_identifiable 1696、not_established 888、effect_only 744、insufficient 544；established 0 |
| 有方向含义的锚定等级（Convergent / Single-layer / Conflicting） | 15 → 0 | 37 → 0 |
| 等级转移 | Single-layer → Omics-single-layer 11；Conflicting → Omics-discordant 3；Convergent → Omics-concordant 1；No-signal → Omics-* 60 | Single-layer → Omics-single-layer 37；Exploratory → No-signal 493（锚点不可辨识或一组少于 2 只，且组学没有方向提示）；Exploratory / No-signal → Omics-* 38 |
| 可靠性表 | 59 个（组学, 组织, 比较）全部 ok | 122 个全部 ok |

说明：

1. 两个队列都没有锚点成立（established）的格子，有方向含义的结论只能来自 effect_only 锚点（设计等级 D，最多 Exploratory）或组学一致性等级。
2. 可靠性表全部为 ok，是因为课题已知的样本问题（第 2 节第 4 条）尚未登记到样本表的 `qc_flag`；登记后只需重跑整合即可生效，不必重跑组学。
3. 队列 B 的 `reversal` 共 213,212 行，其中 210,892 行因疾病效应未建立（HFD_vs_LFD 未达 padj < 0.05，MFA 为每组 2 对 1）归入 undetermined；
   通过等效性检验判为 not_reversed 的 1,541 行，低功效而不能判定的 546 行。课题配置的 `integration.modes` 还没有列出 `reversal`，
   `validate` 会给出提示。

### 11.4 批次 2 的实现（2026-11-15）

批次 2 中具备条件的 G3、G7、G8、G6、E3、E5 已实施；G1 仍待 G7 的标志物表。

| 项目 | 实现 | 未做的部分（原因） |
|---|---|---|
| G3 | `core/statistics/moderated.py`（limma 方差压缩，从蛋白组移出，蛋白组结果不变）；`core/statistics/linear_model.py`（多特征最小二乘 + 方差压缩，支持缺失值，设计不满秩时说明原因）；`design.models`（公式 `~ a + b`、`~ a * b`，协变量自动加入；连续型数据层用上述模型，计数用 PyDESeq2）；有交互项时主效应对另一因素等权平均，`estimand_type` 才是真正的 pooled_main；差异结果增加 `estimand_type`、`n_total`；`design.adaptive_shrinkage` 时增加 `lfsr`（ash 正态混合版，证据画像据此计算方向可信度）；`qc/design_power.parquet`（样本量、组内标准差中位数、80% 功效下的最小可检出效应、秩检验可达最小 p） | 代谢组与 MFA 的两两比较仍不读协变量（公式模型中已读取）；转录组的公式模型不报告交互项 |
| E5 | `stats/robustness.py`：入选特征逐只去掉样本后反号或减半记 `single_subject_driven`（写入 `flags` 与 `robustness.parquet`，并给出影响最大的样本）；以正态—正态收缩估计作为可信真实效应，给出功效、Type S 与 Type M；`flux_anchored` 增加 `single_subject_driven_layers` | 只看入选特征；计数数据用 log2 CPM 的组均值差，与 DESeq2 的估计不同 |
| G8 | `derived` 配置与 `stats/derived.py`：sum、ratio、step_ratio（`source: network` 由 MFA 网络生成）、fraction、index（按丰度加权的链长 / 双键数）、module_score；各类型在“组织 × 比较”内单独校正；features 记 `definition`、`generated_from`、成员的最差鉴定等级；成员不足时不计算并说明原因；网络生成的 step_ratio 只作为关联读数进入 `flux_anchored`；派生特征同样进入公式模型、稳健性、逆转与复现 | 通路层面的派生特征待 G9 |
| G7 | `core/refs.py` 与 `omicsweave refs list / build / audit`：基因、蛋白、代谢物三张知识表的生成器与审计；`validate` 在有结果包时报告登记号与数据不一致；`core/idmap/lipids.py`（LIPID MAPS 简写解析，含羟基、醚键、N-酰基与常见厂商写法）；代谢物注册表由 RefMet 生成（`build_metabolite_registry`）；“模型包”接口（`mappings/model_package.py`），证据构建器的反应知识改由模型包提供；同工酶加权（`isozyme.minor_share`） | 组织标志物表与同源映射（G1 的前置）；基因组尺度模型的基因—反应规则需另写模型包；代谢物名称查找仍以手写表为主，RefMet 注册表（整合用到的 642 个名称中 337 个匹配到 ID）先用于审计：与手写表有 22 处不一致。2026-11-16 逐条在 KEGG、HMDB（经 Metabolomics Workbench 与 PubChem）、ChEBI 核对：10 个条目是明确错配并已修正（如 HMG-CoA 原为 C00418 甲羟戊酸与 (S)-甲基丙二酰辅酶 A 的 ChEBI、甘油醛-3-磷酸原为 C00354 果糖-1,6-二磷酸、叶酸原为尿黑酸、苹果酸原为 D-苹果酸、糖原的 HMDB 原为 3-己酮），“Malic acid”随之从 FUMH_m / MDH_m / MDH_c 的 related 移到苹果酸所在的角色。旧结果包中运行器按旧表写入的错误 ID 由 `flux_anchored` 按 `CORRECTED_IDS` 只更正这些确定错误的值；在队列 A 上，证据等级不变，甲羟戊酸不再被当作 HMG-CoA（8 个格子少匹配 1 个），测得的苹果酸在苹果酸酶、柠檬酸合酶、PEPCK 等 6 个反应中被匹配到（84 个格子多匹配 1 个，此前被错误的 D-苹果酸 ID 挡住）；其余 12 处为同一化合物的立体异构体或异头物，登记在 `refs.ACCEPTED_METABOLITE_VARIANTS`，审计不再报告 |
| E3 | `integration/evaluation/decoys.py`：丰度匹配的诱饵映射，按层与锚点状态汇总目标与诱饵支持率；`flux_anchored` 参数 `decoys`，种子登记在运行清单 | 只做反应 → 特征映射层；评估卡的其他指标（E1、E2）待 G3 的半合成框架 |
| G6 | `omicsweave replicate` 与 `project.yaml` 的 `replication`：跨项目来源、组织名对应、方向性部分合取检验（Benjamini 与 Heller）与 BH、随机效应合并、跨层（蛋白对 RNA）对应、随机重排基线；结果写入项目级结果包 `_project/` | 证据画像读取 `replicated`（G5 第 4 项）尚未接入 |

另外，`omicsweave run` 在各组学运行后统一做后处理（派生特征、公式模型、可检出效应、稳健性、可靠性）；`omicsweave stats --source <运行>`
对已有结果包重跑后处理，不必重跑预处理与两两比较。

### 11.5 批次 2 在现有课题上的验证

验证在 scratch 目录中进行：读取两个课题现有的结果包（队列 A 组学 `run_20261112_102147_5d76`、队列 B 组学
`run_20261112_113117_1269`），在内存中加上验证用的配置（`~ diet` 公式模型、`adaptive_shrinkage`、派生特征声明、诱饵 3 次），
不改动课题的配置与结果目录。

| 能力 | 队列 A 上的结果 |
|---|---|
| 后处理（G3、E5、G8） | 67 秒。公式模型在蛋白组 6 个、代谢组 6 个、MFA 10 个组织拟合成功；全部差异结果行有 `lfsr`。可检出效应：蛋白组每组 3 只时 80% 功效的最小可检出 log2 倍数变化为 0.53–1.61，代谢组 1.01–2.42，MFA 2.4–5.8 |
| E5 稳健性 | 转录组 10,086 个入选基因中 1,017 个去掉一只动物后反号或减半；蛋白组 610 个、代谢组 61 个中没有。MFA 313 个入选反应中 163 个由单只动物驱动，Type S 中位数 0.35、夸大倍数中位数 18 |
| G8 派生特征 | 网络生成 step_ratio 108 个（组织 × 特征），另有链长 / 双键指数 12 个、醚 PE 占比 6 个 |
| G7 同工酶加权 | 7 个格子的支持因只来自次要同工酶而撤销，如肝脏 PYK_c（主要为 Pklr）、HEX_c（Gck）、ICDH_m（Idh1） |
| E3 诱饵映射 | 转录组目标支持率 15.7%，诱饵 15.4%（比值 0.98）；蛋白组 0.8% 与 0.9%；代谢组 0% 与 1.5%。即在 3 对 3 的设计下，组学层“支持”的出现率与随机映射没有差别 |
| G3 锚点 | 两两比较中仍没有成立的锚点。加上 `~ diet` 模型后，ND_vs_LFD 的合并估计在心脏有 3 个反应锚点成立（只有 MFA 有 ND 组，等级为 Single-layer）。验证中发现并修正：合并估计最初未校正示踪剂、也没有继承两两比较的可辨识性标记，曾使肾脏 GLNS_c（不可辨识的通量）成为锚点；现在 MFA 模型自动加入 `protocol` 协变量并沿用可辨识性与示踪剂一致性标记 |
| 运行时间 | `flux_anchored` 加上合并比较约 5.5 分钟；诱饵每次重复约增加 4.5 分钟（3 次共约 20 分钟），所以默认关闭 |

| 能力 | 队列 A × 队列 B（HFD_vs_LFD；B 的 renal_cortex 对应 A 的 kidney） |
|---|---|
| G6 复现 | 146,323 个对应上的特征与通路，184 个复现（特征 178、通路 6），其中心脏 RNA 167 个。随机重排的同号基线为 0.47–0.53（小集合除外）；实际同号比例：代谢物 0.61–0.69、心脏 RNA 0.69、其余组织 RNA 0.45–0.57、蛋白对 RNA 的跨层对应 0.45–0.53 |

E3 的结果说明，在现有样本量下，组学层的单层“支持”大多不能与随机映射区分，引用证据等级时应以锚点状态、复现（G6）
与稳健性（E5）为准。队列 B 心脏 RNA 有已知的离群样本（第 2 节第 4 条），其复现结果在登记 `qc_flag` 之前不宜采信。

### 11.6 未完成的事项（2026-11-17 核对）

| 项目 | 现状 | 出处 |
|---|---|---|
| 两个 DIO 课题登记 `qc_flag` | 两个课题的 `samples.tsv` 暂未手工登记 `qc_flag`，但 11.8 已在通用流程层面实现转录组与 MFA 的自动 `sample_quality.parquet` 生成与严重异常样本过滤，重跑 `omicsweave run` 即可自动驱动可靠性表（两个课题维持各自独立分析，不合并为双来源） | 11.3 说明 2、11.8 |
| G1 样本身份与组成 QC | 转录组与 MFA 已落地（11.8）：`transcriptomics/sample_qc.py` 与 `mfa/sample_qc.py` 自动产出 `qc/sample_quality.parquet` 并接入 `core/reliability.py`；蛋白组 / 代谢组样本级组成 QC 与跨组学性别/身份一致性核验待后续按需补齐 | 11.1、11.8 |
| G5 第 3、4 项 | 锚点类型泛化待 G12；证据画像尚未读取复现结果（`replicated`） | 11.2、11.4 |
| E7 登记的配置缺口 | 剩 3 项：`defaults.padj_method`、`defaults.random_seed`、`design.stratify_by`（11.2 中的 `factors.*.reference` 已不在缺口表中）；`design.covariates` 在代谢组与 MFA 的两两比较中仍未接入 | `tests/core/test_config_consumers.py` |
| 代谢物名称查找 | 仍以手写表为主，RefMet 注册表用于审计。审计报出的 22 处不一致已全部处理：10 处修正，12 处登记为可接受的立体异构体或异头物 | 11.4 |
| 批次 3（G9、G10、G12、E1、E2） | G10 中不依赖新数据的部分已实施（11.7）；G10 的 MID 契约、血浆富集平台、标记比与贡献矩阵待示踪批次数据。G9、G12、E1、E2 未开始 | 11.1、11.7 |
| 通路层面的完整展示（G13–G17、E8） | S0（G17，示例面板的 6 个反应）、S1（成员行与 `reaction_dossier` 模式）、S2（其他来源、逐只动物表、通路面板、图 A–C 与导出表）、S3（分轴综合评估与结论等级）、S4（E8 校准：诱饵、对照集、不变性检查）、S5（G16：节点分流比、模块旋转检验、交叉点分析、图 D）已实施；S4 的组别置换待 E1，S5 的标记层分流比待 G10，S6 待实施。G9 可复用 S5 的旋转检验（`core/statistics/rotation.py`） | [反应档案与分轴综合评估](omicsweave_reaction_dossier_and_assessment.md) 8.1 |
| 两个 DIO 课题的组学结果包 | 2026-11-12 生成，早于批次 2：没有 E5 的 `robustness.parquet`，任何特征都不会被标为单只动物驱动。在 scratch 中只重跑 MFA 后，队列 A 档案的 L2 由 144 个降为 99 个（48 个格子带上 `single_subject_driven`）。引用档案或 `flux_anchored` 的结论前，宜在课题目录中重跑 `omicsweave run` | 11.7 |

### 11.7 G10 中不依赖新数据的部分（2026-11-18）

两个 DIO 课题的示踪文件只有逐次拟合的通量（队列 A 每只动物约 150–175 次拟合，队列 B 1–50 次），没有同位素分布与血浆富集，
G10 的 MID 部分不具备条件。以下三项只用到拟合结果，先行实施：

| 项目 | 实现 | 与 G10 设计的差别 |
|---|---|---|
| 通量组 | `mfa/flux_groups.py`：把一个模型全部样本的逐次拟合按组织堆叠，净通量在所有拟合中相关系数绝对值 ≥ 1 − 1e−6 的反应视为化学计量上锁死，按传递闭包合并；所有拟合中不变（常贴边界）的通量不参与；多示踪剂时一对反应须在每个两者都有变化的模型中都锁死。`tables/flux_groups.parquet` 记录组、成员、相对首个成员的斜率与截距。单示踪剂与合并结果的 BH 校正中每组只计一次（组内取最小 p），差异结果多 `flux_group`、`flux_group_size`；参数 `analysis.flux_groups`（`enabled`，默认 true；`min_abs_correlation`） | 设计说“化学计量上锁死的反应”。这里由拟合结果识别，而不是由网络的零空间推出：几百上千次独立拟合中的精确线性关系只能来自模型结构，且不需要完整的网络与交换反应定义 |
| 标准化尺度的示踪剂合并 | `multi_tracer`：各示踪剂的 Hedges' g 按 1/Var(g) 合并（`smd`、`smd_se`、`smd_heterogeneity_*`、`tracer_smd`），合并的 `effect_size` 改为该值；示踪剂冲突（`multi_tracer_conflict`）要求**方向相反**且标准化尺度的异质性 p < 0.10，方向一致而幅度不一致时只标 `tracer_magnitude_heterogeneity`。主检验（分层置换 p）不变 | 设计说“通量效应量改用相对参考通量的差值”。核对数据后改用标准化效应：同一组织内葡萄糖与其他示踪剂的通量之比在反应之间差异很大（log2 比的四分位距约 2），不是一个共同的尺度因子，除以某个参考通量只能去掉中位偏移；心脏糖酵解的净通量还在两组间变号，比值没有定义。g 在变号时有定义，是不同尺度测同一量时的标准合并方式。置换检验对任何权重都精确，所以主检验不改 |
| 逐只动物的拟合分布 | `sample_flux_matrix` 另计全部拟合点的四分位数；`tables/sample_flux_fits.parquet` 写出每只动物每条反应的代表值、拟合次数、四分位数、IQR 过滤后的离散度与是否贴边（`at_bound`）。反应档案的逐只动物表读入 `fit_low` / `fit_high` / `at_bound`，图 B 的通量面板画误差线、贴边的动物画空心点 | 设计中的区间取自求解器的非线性置信区间。旧格式的导出没有求解器区间；拟合集合的四分位距是“数据下这条通量被约束得多紧”的代理，不是置信区间 |

另外修正了反应档案 F 轴的一处错误：净通量在两组间变号时，旧流程的 log2FC 为两组绝对值之比，符号不一定与净通量的变化同向
（甘油示踪下心脏 ENO_c 由 −28.5 升到 +8.0，log2FC 为 −1.84）。S3 用 log2FC 的符号判断各示踪剂的方向，把同向的葡萄糖与甘油判成了
相反。档案中的通量行现在按组均值之差定号（幅度不变），标记 `net_flux_sign_change`。

在两个课题上（scratch 中只重跑 MFA，与原组学结果包拼成新的运行）：

| 读数 | 结果 |
|---|---|
| 通量组 | 队列 A（三种示踪剂合并）133 组、336 条反应；队列 B 87 组、221 条反应。例如 PFK_c 与 FBA_c、AKGD_m 到 FUMH_m 四步、心脏的 PC_m 与 GLUD_m |
| 合并结果的 FDR | 队列 A 有 p 值的 153 行中，144 行的 padj 下降、6 行上升（同一组内重复的小 p 不再抬高其他检验的秩）；入选的 313 条通量没有变化。队列 B 每组 3 只，全部只给效应量，不受影响 |
| 示踪剂冲突 | 队列 A 由 10 个降为 1 个。解除的 9 个都是两种示踪剂同向（心脏 ENO_c、GAPD_c、PGK_c、PGM_c、TPI_c、GPS_c、GLYCEROL_exchange；肝脏 ACONT_m、ICDH_m，后两者甘油的 g 为 −0.14，接近 0）。保留的是肝脏 CS_m（葡萄糖 g +2.32，甘油 −0.40） |
| 反应档案（队列 A + B） | LX 由 3 个降为 0 个：心脏 ENO_c、PGM_c、TPI_c 在 HFD_vs_LFD 中两种示踪剂同向上升，F 轴为 suggestive，结论 L2。L2 由 144 降为 99，原因不是 G10：新的 MFA 结果包带有 E5 的稳健性表，48 个格子的合并通量由单只动物驱动（见 11.6） |
| 运行时间 | 队列 A 的 MFA 52 秒，队列 B 85 秒 |

验证：新增 9 项测试（`tests/mfa/test_tracer_readouts.py` 8 项，`test_multi_tracer.py` 1 项）：锁死关系的识别（比例与仿射、不变的通量不参与、
多模型须都锁死）、组内只计一次的 BH、尺度相差 50 倍而同向的示踪剂不判冲突、方向相反仍为冲突、结果包的两张新表与分组校正、档案读入
误差线与贴边标记、变号时的定号；导出命令的预期文件随结果包增加两张表。

### 11.8 G1 通用样本质检与评估规则硬化（2026-11-23）

在两个独立 DIO / GLP-1RA 课题（各自作为单来源独立分析，不合并双来源）的多组学核验中，暴露出四类具有跨项目普遍性的流程问题：
转录组解剖混入与组织身份缺失无自动质检过滤、虚拟/边界反应被升为 L2、集总多步/强诱导型酶被可替代同工酶丰度规则误杀、以及 MFA 大队列全局稀释掩盖单组上界贴边与跨批次整鼠通量尺度漂移。上述改进全部以零项目硬编码、可配置参数的方式写入通用流程：

| 类别 | 模块与实现 | 配置项（默认值） |
|---|---|---|
| **G1 转录组通用样本 QC 与自动过滤** | `transcriptomics/sample_qc.py` + `transcriptomics/runner.py`：由计数矩阵自动推导成对异源标记与本组织特异标记，检测 `cross_tissue_contamination`（跨组织 RNA 混入）、`own_marker_retention`（组织身份缺失，整组同时失效亦可检出）、`transcriptome_correlation`（全表达 + 高变基因成对相关离群）与 `mitochondrial_fraction`（绝对 ≥ 88% 或相对同组织中位数极端偏移，不误杀天然高线粒体组织）；写出 `omics/transcriptomics/qc/sample_quality.parquet` 与 `sample_qc_metrics.parquet` 供 `core/reliability.py` 自动计算 `contrast_reliability.parquet`；在差异分析前过滤 `flag == "exclude"` 样本（`library_sizes.parquet` 保留过滤前全量原始样本） | `omics.transcriptomics.analysis.sample_qc`（`enabled: true`、`filter_excluded: true` 等阈值） |
| **虚拟 / 不可评估反应门控** | `integration/dossier/assessment.py`：`decide(..., semantics=...)` 与 `assess()` 增设硬门控——`direction_evidence == "not_assessable"`（如 `BIOMASS_REACTION_*`、`Salvage_c`、`*_protein_input`）禁止升为 `L1` / `L2` / `L3` / `LX`，一律锁在 `L0`（`claim_direction = 0`，`phenomena = "not_assessable_reaction"`） | 内置硬规则（依据反应语义注册表） |
| **G7 `minor_isozyme` 双层修复** | 1. 生物化学成员关系：`integration/mappings/reaction_members.py` 新增 `pathway_step` 关系；`reaction_members.yaml` 将集总多步/不可替代诱导酶反应（`THERMOGENESIS_UCP`、`FA_esterification_c`、`HMGCR`、`BILE_acid`、`ETHER_LIPID`）标为 `pathway_step`，并将泛表达 `Ucp2`/`Ucp3` 移至 `THERMOGENESIS_UCP.related_isozymes`，同步更新 `flux_gene_mapping.yaml` 与 `flux_protein_mapping.yaml`；<br>2. 算法层：`reaction_dossier.py` 与 `flux_evidence.py` 新增条件活跃丰度占比 `active_abundance_share = max(pooled, share_num, share_den)`；`assessment.py` 的 `_judge_members()` 豁免 `pathway_step` 的同工酶丰度压制，普通 `isozyme` 采用 `eff_share = max(abundance_share, active_abundance_share)` 并在无主要同工酶反向时支持多个同向过 FDR 成员合计占比达标放行，同时保留对痕量反向诱饵（如 `Pkm` 对 `Pklr`）的拦截 | `integration.parameters.reaction_dossier.isozyme.minor_share`（`0.2`） |
| **G1 / G10 MFA 上界相对容差、单组贴边与尺度漂移 QC** | 1. `mfa/mfa_experiment_config.py` 与 `mfa/df_analysis.py`：新增 `bound_rel_tolerance = 0.02` 与上界点质量检测 `_has_upper_bound_pileup`（捕捉 3000 上界附近的 `2965.6`/`2998.5` 等且不误伤无盒约束拟合），输出 `at_upper_bound`；<br>2. `df_analysis.py`、`mfa/export.py` 与 `mfa/multi_tracer.py`：在每个比较的两组内分别计算 `contrast_frac_at_bound = max(num_frac, den_frac)`，任一组 ≥ `bound_pinned_min_fraction`（0.5）即升级为 `bound_pinned`（标记 `contrast_bound_pinned` / `{tracer}:contrast_bound_pinned`），消除多组全局稀释；<br>3. `mfa/sample_qc.py` 与 `mfa/runner.py`：产出样本级 `omics/mfa/qc/sample_quality.parquet`（`flux_scale_factor`、`upper_bound_saturation`、`fitting_loss`）与比较级 `qc/contrast_scale_qc.parquet`（`|Δlog2_total| ≥ 0.30` 且 `|Δmed_log2| ≥ 0.15` 同号时标记 `flux_scale_drift` 并写入差异表 `flags`）；<br>4. `assessment.py` 与 `reaction_dossier.py`：新增 `rules.flux_scale_drift_gate`（默认 `true`），对存在尺度漂移的比较将 F 轴降级为 `undetermined` 并在 `phenomena` 中明示 | `omics.mfa.analysis.identifiability.bound_rel_tolerance`（`0.02`）；`omics.mfa.analysis.sample_qc`；`integration.parameters.reaction_dossier.rules.flux_scale_drift_gate`（`true`） |

在两个现有课题结果包上的只读实测对照（不改动 `local/projects/`）：

| 核验项 | 改进前 | 改进后（通用规则实测） |
|---|---|---|
| 队列 B 转录组样本 QC（89 样本） | 无自动样本过滤与 `sample_quality` | `exclude` 8 例（含 `sq_adipose` Sema 肌混入 `m83/m84`、淋巴结混入 `m118/m125`、`renal_cortex` `m83/m120` 等）、`flag` 7 例（含 `heart` LFD `m129`、Sema 肌混入 `m84`）；`contrast_reliability` 自动阻断 `heart` 全部 LFD/Sema 比较及 `sq_adipose` Sema 比较（`qc_confounded`/`qc_insufficient`），同时保留可恢复比较（`downweighted`） |
| 队列 A 转录组样本 QC（51 样本） | 无自动样本过滤 | 自动 `exclude` 全部 3 例异常 LFD 心脏（`m63/m64/m65`，保留比 0.13–0.23），健康高线粒体 HFD 心脏零误杀，`heart HFD_vs_LFD` 自动标为 `qc_confounded` |
| 虚拟反应门控（队列 B `not_assessable`） | 23 格升为 L2 | **0 格**（全部锁为 `L0 / not_assessable_reaction`） |
| `minor_isozyme` 修复（队列 B） | `sq_adipose THERMOGENESIS_UCP` 4 比较、`liver FA_esterification_c` 3 比较均被误杀为 `L0`；全表 `minor_isozyme_driven` 316 格 | `sq_adipose THERMOGENESIS_UCP`（Mazdu vs Saline/Sema `L1(+1)`、HFD vs LFD `L1(-1)`）与 `liver FA_esterification_c`（HFD vs LFD `L1(+1)`、Sema `L1(-1)`）均恢复；全表降为 240 格，真诱饵（`AMPD_c`/`ENO_c` 等）继续拦截 |
| MFA 上界贴边与单组贴边（队列 B） | 上界命中 13 次；`liver G6Pase_c [Mazdu_3W_vs_1W]` 判 `well_constrained` → `L2(+1)` | 上界 2% 相对容差命中 **27 次**；`liver G6Pase_c [Mazdu_3W_vs_1W]` 因 3W 组 `2/3` 撞顶（`3000.0`/`2965.6`）升级 `bound_pinned` → `L0`，未贴边比较 `Mazdu_vs_Saline_1W` 保持 `well_constrained` |
| MFA 尺度漂移 QC（队列 B & A） | 队列 B `3W_vs_1W` 产出 76 个 L2 通量声明 | 队列 B 精准标出 `Mazdu_3W_vs_1W_in_HFD`（`+0.410/+0.174`）与 `Saline_3W_vs_1W_in_HFD`（`+0.342/+0.167`）为 `flux_scale_drift`（L2/L3 降为 **0**），其余 7 个同周期比较全干净；队列 A 8 项中仅标出 `lactate HFD_vs_LFD`（`+1.358/+0.217`），其余 7 项干净 |

验证：新增 `tests/transcriptomics/test_sample_qc.py` 与 `tests/mfa/test_mfa_sample_qc.py`，并扩展 `tests/mfa/test_flux_identifiability_bounds.py`、`tests/integration/test_dossier_assessment.py`、`tests/integration/test_reaction_members.py`、`tests/integration/test_isozyme_decoys.py`；14 个核心测试套件零回归通过。

### 11.9 G1 / G2 / G18 转录组样本质控双重惩罚修正、分级风险标注与通路图降噪（2026-11-24）

#### 问题背景与实测影响

对两个独立课题结果包的只读核查发现，11.8 引入自动质检与过滤后，上游差异分析与整合可靠性之间存在评价口径错位，导致大量具备强统计显著性的转录组比较被整组抹除为 `qc_blocked`，并在通路图全部反应卡片上反复绘制琥珀色感叹号（`⚠`），严重干扰整体生物学判读：

1. **剔除后双重惩罚（`core/reliability.py`）**：转录组已启用 `filter_excluded: true` 在 PyDESeq2 建模前移除 `exclude` 样本，但 `contrast_reliability()` 仍按全量 `samples.tsv` 分组计数计算 `flagged_num` / `flagged_den`，导致已洁净建模的比较（每组均 ≥ 2 只）再次因“原始设计两组标记率不均”被判 `qc_confounded` 封杀。
2. **轻度提示按严重污染同等封杀**：任何未达 `exclude` 的自动提示（`flag`）与人工剔除等价扣除 `usable` 并触发 `qc_confounded`，例如某单侧样本出现极轻微异源组织痕量表达，整个组织成百上千个完全不相关的代谢基因全部丧失评估资格。
3. **`mtRNA` 双侧误标与全池高变基因相关性混淆**：`sample_qc.py` 原先对 `mitochondrial_fraction` 做双侧 `abs(rz)` 判定，把线粒体占比偏低的正常样本（如队列 B `heart m80`，mtRNA `6.5%` 对同组织中位数 `36%`）误判 `flag`；同时 `var_pairwise_corr` 与 `joint_degradation` 在计算全组织中位相关性时未参考实验分组 `design.factors`，将强处理效应导致的组间分化误判为 `flag`（如队列 A `liver` / `kidney` HFD 对 LFD）。
4. **“阻断即抹除”覆盖真实统计，且通路图遍布 `⚠` 噪声**：`reaction_dossier.py` 对 `qc_confounded` 不区分是否存在有效 `log2fc` / `padj` 一律强制覆写为 `qc_blocked`（清空 `implied_direction` 并禁止参考）；`primitives.py`、`pathway_map_plots.py`、`pathway_map_html.py` 与 `dossier_plots.py` 对每个 `qc_blocked` 卡片均绘制高对比度琥珀色 `⚠`，淹没通量与蛋白层主信号。

受上述过严门控影响的真实显著差异转录组比较（只读统计结果包，仅列 `padj < 0.05` 且 `|log2FC| >= 0.58` 的显著 DEG 数）：

| 队列结果包（只读） | 受阻比较 | 建模样本数 | 显著 DEG 数 | 原触发原因 |
|---|---|--:|--:|---|
| 队列 A（`run_20261124_155042_tracer-views`） | `liver HFD_vs_LFD` | 3v2 | 2,740 | LFD `m65` 已被 PyDESeq2 剔除（4v3 → 3v2）仍遭双重惩罚 |
| 队列 A | `liver HSD_vs_LFD` | 3v2 | 1,111 | 同上（已剔除样本导致原始组间标记率失衡） |
| 队列 A | `kidney HFD_vs_LFD` | 3v3 | 2,319 | 强 HFD 处理效应拉低全池 `var_pairwise_corr` 触发软 `flag` |
| 队列 B（`run_20261123_124248`） | `liver Mazdu_vs_Sema_in_HFD_3W` | 2v3 | 3,915 | `m125` 已被剔除仍触发双重惩罚（部分组兼有 `mt_flag`） |
| 队列 B | `heart Sema_vs_Saline_in_LFD_3W` | 3v2 | 1,830 | `heart m80` mtRNA 偏低（6.5%）被双侧 `abs(rz)` 误标 + 已剔除混杂 |
| 队列 B | `heart Mazdu_vs_Saline_in_LFD_3W` | 2v2 | 1,257 | `heart m80` 低 mtRNA 误标与组间相关性误标 |
| 队列 B | `heart Saline_in_HFD_vs_LFD_3W` | 3v2 | 581 | 同上（已剔除后可用样本 ≥ 2 只仍被整组阻断） |
| 队列 B | `liver Sema_vs_Saline_in_HFD_3W` | 3v4 | 146 | 痕量异源标记 / 轻度 `flag` 导致整层归零 |
| 队列 B | `sq_adipose Mazdu_vs_Saline_in_LFD_3W` | 2v2 | 95 | 已剔除肌组织混入样本 `m83/m84` 后仍遭双重惩罚 |

#### 四层通用改进实施

| 层次与编号 | 模块与实现 | 配置项 / 契约变更 |
|---|---|---|
| **G2 剔除后双重惩罚修复与软硬标记分级** | `core/reliability.py`：新增 `analyzed_sample_sets(run_dir)` 读取 `omics/*/data/*/matrix.parquet` 实际参与差异分析的保留样本列；`contrast_reliability(..., *, analyzed_samples=None)` 与 `build_reliability()` 将样本标记细分为三类：<br>1. **已预过滤样本**（`pre_filtered`：在 `hard_set` 中但已不在保留矩阵中）：不参与当前比较的 `qc_confounded` 失衡判定，扣除后每组仍 ≥ 2 只时记 `status = "downweighted"`；<br>2. **未过滤严重/人工标记**（`unfiltered_hard`：`flag == "exclude"` 或 `samples.tsv` 人工指定且仍留在分析矩阵中）：继续严格执行 `qc_confounded`（比例差 ≥ `imbalance`）与 `qc_insufficient`（扣除后任一组 < 2 只）；<br>3. **轻度自动提示**（`soft`：仅来自自动 `sample_quality.parquet` 且 `flag != "exclude"`）：永不单独触发 `qc_confounded` 或 `qc_insufficient`，按每个软标记样本 `0.5` 权重温和折算 `reliability`（下限 `0.25`） | 自动推断（零项目硬编码）；`qc/contrast_reliability.parquet` 口径保持兼容 |
| **G1 转录组样本 QC 校准与异源基因级隔离** | 1. `transcriptomics/sample_qc.py`：`DEFAULT_SAMPLE_QC_CONFIG` 新增 `"mt_flag_low": False`，默认仅对同组织偏高 mtRNA（衰老/裂解/空泡）告警；`evaluate_sample_quality(..., groups=None)` 计算同组织**同实验组**生物学重复间中位相关 `r_within`，当 `r_within >= var_corr_flag_min` 时抑制跨组分化引起的 `var_pairwise_corr` 与 `joint_degradation` 软标；`SampleQCResult` 新增 `foreign_genes`（来自 `_derive_pair_marker_sets()`）；<br>2. `transcriptomics/runner.py`：由 `ctx.experiment.config.design.factors` 构建 `sample_groups` 传入 `evaluate_sample_quality()`；新增 `_flag_foreign_markers()`，当差异比较仍保留带 `cross_tissue_contamination` 轻度标记的样本时，对与混入组织均值同向且显著/名义变化的异源特异基因，在 `differential.parquet` 的 `flags` 中精准追加基因级 `"foreign_tissue_marker"`（不影响同样本内其余洁净代谢基因） | `omics.transcriptomics.analysis.sample_qc.mt_flag_low`（默认 `false`）；差异表 `flags` 新增 `foreign_tissue_marker` |
| **G18 统计状态保留与分级风险软封顶** | 1. `integration/modes/reaction_dossier.py`：`_state()` 将样本不足（`n < 2` 或 `qc_insufficient`）独立归为 `"insufficient_n"`；若行存在未过滤 `qc_confounded` 风险但已有有效差异统计（`padj` / `pvalue` / `log2fc`），保留真实统计状态（`changed` / `nominal` / `no_evidence`，禁用 `equivalent`）并在 `state_reason` 附注“含未过滤样本 QC 混杂风险”，`_finish()` 为其追加 `"qc_caution"`；<br>2. `integration/dossier/assessment.py`：定义统一风险标记 `_CAUTION_FLAGS = ("single_subject_driven", "foreign_tissue_marker", "qc_caution")` 与 `_LAYER_QC_CAUTION = ("foreign_tissue_marker", "qc_caution")`；新增 `_caution_info()` 与 `_cap_if_cautioned()`，在 `_judge_members()`、`regulation_axis()`、`pool_axis()` 中当支撑 `changed` 的成员均带风险标记时软封顶为 `suggestive` 并保留对应 `flags`；`decide()` 禁止单靠含 `_CAUTION_FLAGS` 的 C 轴证据升为 `L1`，并将 `foreign_tissue_marker`、`qc_caution` 写入 `phenomena` | 成员表 `cell_state` 新增/解耦 `insufficient_n`；`flags` / `phenomena` 新增 `qc_caution`、`foreign_tissue_marker` |
| **通路图与档案图去感叹号降噪** | `visualization/primitives.py`、`integration/plotting/pathway_map_plots.py`、`pathway_map_html.py`、`dossier_plots.py`：彻底移除针对 `qc_blocked` / `insufficient_n` 的反应卡片级、四轴小条、聚合方框与交互 HTML 琥珀色感叹号三角（`⚠`），统一改用低对比度浅灰底色（`missing_fill`）配合中性双斜杠 `//`（图例标签统一为 `Low-n (n<2)`；刻意避开 HTML 逐只动物散点图已使用的贴边空心菱形 `◇`）；`integration/pathway_map.py` 新增 `_qc_disclosure()`，在每页图注首行集中附带单处纯英文/ASCII 比较级汇总（`Sample QC: RNA XvY usable of AvB (status)`） | 视觉编码升级（图面零感叹号噪声，兼容英文出版字体与 ASCII 校验） |

验证：扩展 `tests/core/test_reliability.py`、`tests/transcriptomics/test_sample_qc.py`、`tests/integration/test_reaction_dossier.py`、`tests/integration/test_dossier_assessment.py`、`tests/visualization/test_primitives.py`，并与 `tests/integration/test_pathway_map.py`、`tests/core/test_config_consumers.py`、`tests/integration/test_candidates.py`、`tests/integration/test_reversal.py` 共 9 个单元测试套件联合运行，**93 项测试全部通过（OK）**。

### 11.10 G9 / 工程性能：预排序 GSEA 双层开关与默认关闭优化（2026-11-24）

#### 耗时特征与依赖核查

对两个课题正式结果包文件落盘时间戳的只读核查表明，基于全排序特征列表与 `1000` 次随机置换的预排序 GSEA（`gseapy.prerank`）占据了全流程绝大部分运行时间，而核心多组学评估与通路图主体绘制均不依赖它：

| 课题与运行包（只读） | 全流程 `run + integrate` 耗时 | 转录组 PyDESeq2 差异分析 | 转录组富集（ORA + GSEA） | **转录组 GSEA 占全流程比重** |
|---|--:|--:|--:|--:|
| `dio_glp1ra_project`（`run_20261123_124248`，30 个比较） | 170 分 10 秒 | 3 分 12 秒 | 159 分 13 秒 | **93.6%** |
| `dio_project`（`run_20261123_124459`，10 个有效比较） | 28 分 41 秒 | 1 分 13 秒 | 21 分 56 秒 | **76.5%** |

- **核心评估与通路总图零依赖**：`reaction_dossier`（F/C/R/M 四轴判定、`L0–L3/LX`、节点分流比、有向模块旋转检验、候选表）及通路图主画布（反应边、四格微条、代谢物节点、分专题子图、Graphical Abstract）只读 `differential.parquet` 与样本矩阵，不读取 `enrichment.parquet`。
- **探索性融合可由本地 ORA 支撑**：唯一读取 `enrichment.parquet` 的 `pathway_convergence`（及总图 `atlas` 模式右侧栏 `Fusion FDR` 提示）由本地 Fisher ORA 与预排序 GSEA 共同供数；保留毫秒级本地 ORA 并默认关闭预排序 GSEA，可将数小时运行压缩至数分钟。

#### 双层开关设计与实现

| 层级 | 模块与实现 | 默认值与契约说明 |
|---|---|---|
| **统一富集适配器** | `adapters/enrichment.py`：`EnrichmentSettings` 新增 `run_gsea: bool = False`；新增 `resolve_run_gsea(analysis, *, override=None, aliases=(), default=False)` 校验严格布尔类型并按 `命令行覆盖 > analysis.run_gsea > 兼容别名 > default` 解析 | 默认 `false`（关闭置换 GSEA、保留本地 ORA）；非布尔值报 `RunError` |
| **三大组学运行器** | 1. `transcriptomics/runner.py`：`_enrichment()` 解析 `resolve_run_gsea(analysis, override=ctx.run_gsea)`，仅当 `settings.run_gsea` 为真时调用 `prerank_rows()`；<br>2. `proteomics/runner.py`：`_DEFAULTS` 将 `"run_prerank_gsea"` 与 `"run_gsea"` 默认设为 `False`，`resolve_run_gsea(..., aliases=("run_prerank_gsea",))` 兼容旧配置键名；<br>3. `metabolomics/experiment_config.py` + `metabolomics/runner.py`：`MetabolomicsConfig` 新增 `run_gsea: bool = False`，仅开启时调用 `_gsea()` | 三大组学统一默认 `run_gsea: false`；`omics.proteomics.analysis.run_prerank_gsea` 向后兼容 |
| **运行框架、流水线与 CLI** | `core/run.py`：`RunContext`、`start_run()`、`run_experiment()` 新增 `run_gsea: bool | None = None` 并写入 `manifest.parameters["run_gsea"]`；<br>`core/pipeline.py`：`run_pipeline()` 透传 `run_gsea`并在日志命令附加 `--gsea`/`--no-gsea`，`_plan()` 在旧单组学结果包已包含 `enrichment` 步骤时允许零/单表续跑；<br>`cli.py`：`omicsweave run` 与 `omicsweave pipeline` 新增互斥开关 `--gsea` / `--no-gsea` | 命令行显式开关优先级高于 `experiment.yaml` 配置 |
| **通路融合优雅降级** | `integration/modes/pathway_convergence.py`：当源结果包中已有 $\ge 2$ 个非 MFA 组学执行了 `enrichment` 步骤、仅因关闭 GSEA 且无显著 ORA 命中使非空 `enrichment.parquet` 少于 2 份时，不抛错阻断流程，而是写出空 `pathway_fusion.parquet` / `pathway_overlaps.parquet` 并附加运行提示；未满足组学数或跳过 `enrichment` 步骤仍严格报 `RunError` | 保障 `pathway_convergence + reaction_dossier` 联合流水线在零 ORA 命中下平滑完结 |
| **契约表完整列保障** | 1. `contract/tables.py` 与 `contract/__init__.py`：新增 `empty_enrichment()` 与 `_complete_enrichment_columns()`，在 `write_enrichment()` / `read_enrichment()` 时自动补齐全部 23 个 `ENRICHMENT_COLUMNS`（纯 ORA 模式下仍稳定具备 `es`、`nes`、`flags`，纯 GSEA 具备 `query_size`、`universe_size`、`odds_ratio`）；<br>2. `adapters/enrichment.py`、`adapters/gseapy_results.py`、`metabolomics/export.py`：各转换器直接输出完整列；三大运行器执行 `enrichment` 步骤但无显著条目时均落盘 0 行标准契约表；<br>3. `integration/methods/pathway_convergence.py`：`fuse_pathway_evidence()` 空融合结果表固定输出含 `omics_probabilities`、`omics_rank_fractions` 的完整 10 列 | 彻底消除关闭 GSEA 后下游课题脚本（如 `export_key_results.py`）直读 Parquet 时因缺列/缺表触发的 `KeyError` / `FileNotFoundError` |

验证：在 `tests/transcriptomics/test_transcriptomics_runner.py`、`tests/proteomics/test_proteomics_runner.py`、`tests/metabolomics/test_metabolomics_runner.py`、`tests/integration/test_integrate.py`、`tests/core/test_pipeline.py`、`tests/core/test_tissue_aliases.py`、`tests/core/test_config_consumers.py` 7 套件通过（59 项 `OK`）基础上，补充 `tests/contract/test_result_tables.py` 全列结构契约验证并联合 `tests/contract/test_query_diff.py` 等 7 个套件回归，**58 项测试全部通过（OK）**。

---

## 12. 相关文档

| 文档 | 内容 |
|---|---|
| [架构规划 v2](../omicsweave_arch_design_v2.md) | 分层、数据契约、运行模式与阶段划分 |
| [开发规划总纲](../omicsweave_theraweave_development_plan.md) | 结果模型、证据层与交换契约 |
| [整合模块说明](../../src/omicsweave/integration/README.md) | 证据画像与整合模式的现状 |
| [v2 打分矩阵技术文档（归档）](../history/integration_pipeline_evaluation_v2.md) | 已弃用的加和总分，仅供追溯 |
| [反应档案与分轴综合评估](omicsweave_reaction_dossier_and_assessment.md) | G13–G17、E8：通路层面的完整展示与不做加和的评估方案（S0–S5 已实施） |
