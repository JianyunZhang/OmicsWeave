# OmicsWeave 反应档案与分轴综合评估：通路层面的完整展示与不做加和的评估方案

> 文档类型：设计建议（S0 已实施，S1–S6 待实施）
> 版本：v0.2（2026-11-17，增加 8.1 实施记录；v0.1 为 2026-11-16）
> 依据：对队列 A 结果包（组学 `run_20261112_102147_5d76`、整合 `run_20261115_010915_d97b`）与队列 B 结果包
> （`run_20261112_113117_1269`）的只读核查；课题工作区中的 `midterm_report/midterm_key_pathways.md` 与
> `focus_targets/INTEGRATION_ASSESSMENT.md`（未纳入版本库）
> 配套文档：[通用能力规划](omicsweave_generic_capabilities_and_evaluation.md)（本文的编号接续其中的 G1–G12、E1–E7）、
> [架构规划 v2](../omicsweave_arch_design_v2.md)
> 范围：流程能力、评估方法与图表规范。附录 A 的数字只用来说明格子能填到什么程度，不是课题结论

---

## 1. 结论摘要

1. **展示与评估要分开设计。** 展示的目标是完整：一条通路上每个反应的酶（RNA、蛋白）、调控因子、前后代谢物、各示踪剂的通量，
   测到的、没测到的、互相矛盾的都要出现在同一张图上。评估的目标是可解释：每个结论能追溯到具体的格子和规则。
2. **不用加和总分，改用“分轴判定 + 判定表”。** 通量、酶容量、调控状态、底物与产物是四个不同的量，分别给出方向与状态；
   “通量是否改变”只由示踪数据回答，酶和代谢物用来解释通量的变化，或在通量未定时给出可检验的预测（第 3 节）。
3. **需要一个总的判断时，用三样东西代替总分**：结论等级（L3 已证实、L2 有提示、L1 预测、L0 未定）、机制类型（酶量主导、
   代谢调控主导、两者拮抗等），以及一条公开的排序规则。等级的错误率用诱饵、置换和对照集校准。
4. **现有流程离这个目标不远。** 成员级的匹配结果已经算出，只是写表时丢弃；分示踪剂的通量、逐只动物的数据都在结果包里。
   缺的是成员级证据表、通路的声明、带符号的调控因子，以及图（第 2、5 节）。
5. **建议新增 5 项能力和 1 项评估指标**（G13–G17、E8），分 6 步实施，前 4 步不需要新数据（第 8 节）。

| 编号 | 能力 | 目的 |
|:--:|---|---|
| G13 | 反应档案：成员级证据表与逐只动物数据表 | 展示的数据基础 |
| G14 | 通路面板：通路的声明，证据矩阵图与档案页 | 完整展示 |
| G15 | 分轴综合评估：四轴判定、判定表、合取检验、排序规则 | 替代加和总分 |
| G16 | 节点分流比、预先定号的模块检验、交叉点分析 | 通路层的综合读数 |
| G17 | 知识层补充：带符号的调控因子、成员关系、角色统一 | 调控因子能被正确展示与解读 |
| E8 | 判定规则的校准 | 等级有可估计的错误率 |

---

## 2. 依据：现状核查

| # | 现象（流程行为） | 位置 | 对本方案的含义 |
|---|---|---|---|
| 1 | 证据构建器为每个反应生成各层的成员级明细（`tx_details`、`prot_details`、`metab_details`），`profile_table` 写表时整列删除，只留下“与通量同向的特征”字符串 | `integration/modes/flux_anchored.py` 的 `profile_table` | 成员级证据表主要是保留已有结果，不需要新的匹配逻辑 |
| 2 | 层级表只记数量。例：肝脏 PC_m 的代谢组层为“15 个映射、4 个匹配、0 个有信号”，没有记是哪几个、数值多少 | `flux_evidence_table.parquet` | 读者无法区分“没测到”与“测到但没变” |
| 3 | 锚点不成立时，“同向特征”列为空。队列 A 的 1,605 个格子中没有一个锚点为 `established`（`not_identifiable` 807、`not_established` 624、`effect_only` 162、`tracer_conflict` 10、`insufficient` 2） | `flux_evidence_profile.parquet` | 以通量为锚的展示在当前数据上几乎不显示任何成员 |
| 4 | 整合只读合并后的通量表。分示踪剂的结果（`differential_by_tracer.parquet`）、各示踪剂的可辨识性（`tracer_coverage.parquet`）已在结果包中，但没有进入整合与图 | `flux_anchored.run` 读取 `omics/*/tables/differential.parquet` | 分示踪剂的展示不需要重跑 MFA |
| 5 | 各组学都有逐只动物的数据集（`omics/<组学>/data/<组织>/matrix.parquet`）；通量矩阵的样本名含示踪剂（如 `glycerol__male__hfd__1`） | 结果包 | 逐只动物的图可以直接画 |
| 6 | 通量数据集的 `units` 为 `TBD` | `omics/mfa/data/*/dataset.json` | 图的纵轴没有单位，需在 MFA 导出时补上 |
| 7 | 调控因子没有符号。蛋白表的 `upstream_regulators` 在 91 个反应中有 24 个非空，激活与抑制不分（PDH_m 下 Pdk1–4 与 Pdp1–2 并列）；基因表把 Pdk1–4、Pdp1–2 记在 PDH_m 的 `catalytic` 下 | `mappings/flux_protein_mapping.yaml`、`flux_gene_mapping.yaml` | 运行时方向成员以蛋白表的 `direct_enzymes` 为准，没有造成方向误判；但按角色展示前需要统一角色并补符号（G17） |
| 8 | 复合体的亚基与同工酶不区分；DLD 这类多个复合体共用的亚基按普通成员处理 | 同上；同工酶加权（G7） | 复合体不应按“可互相替代的同工酶”加权 |
| 9 | 没有“通路”的定义。反应只按 KEGG 通路名分组，没有顺序，也没有分支节点 | `aggregate_by_pathway` | 需要声明式的通路面板（G14） |
| 10 | 现有整合图（通路 × 条件热图、TCA 示意图、证据画像图）都在层级汇总的粒度，没有成员级或逐只动物的图 | `integration/plotting/` | 需要新的绘图模块 |
| 11 | 另一个队列在另一个项目目录中 | `examples/dio_glp1ra_project` | 档案需要像复现模式一样读取其他来源 |
| 12 | 诱饵评估（E3）：转录组“支持”通量方向的比例为 15.7%，诱饵为 15.4%；蛋白组 0.8% 对 0.9%；代谢组 0% 对 1.5% | 通用能力规划 11.5 | 组学层与通量方向的一致性目前没有区分力，不能作为给通量加分的依据 |
| 13 | 同一化合物可能有多个特征（如肝脏的 “Aspartic Acid” 与 “Aspartic acid”）；知识表的名称与数据的特征名可能不同（Pcx 在蛋白组数据中记作 Pc，靠 UniProt 登记号匹配） | 代谢组、蛋白组数据集 | 档案要同时显示知识表的成员名与实际匹配到的特征，重复特征要标记 |

---

## 3. 评估方法：不做加和时，如何综合评价一个反应

### 3.1 先分清估计对象

“一个反应变强了”可以指三个不同的量：

| 量 | 含义 | 由什么测得 |
|---|---|---|
| 速率 J | 净通量 | 只有示踪数据 |
| 容量 E | 酶量，近似最大反应速率 | 蛋白组；转录组是间接读数 |
| 驱动与调控状态 | 底物、产物、效应物的浓度，酶的修饰状态 | 代谢组；调控因子的量；磷酸化数据（目前没有） |

三者的关系是 J = E · k(底物, 产物, 效应物, 修饰状态)。两个条件之间取对数差：

    Δln J = Δln E + Δln k

两边除以 Δln J，得到层级调控系数 ρ_h = Δln E / Δln J 与代谢调控系数 ρ_m = 1 − ρ_h [1、2]。这个恒等式有三个直接的推论：

- 酶量和代谢物的变化是通量变化的**组成部分**，不是通量变化的独立“证人”。把它们的分数加到通量上，等于把同一个量算了两遍。
- 酶量上升而通量不变，说明 Δln k 下降了同样的幅度。这是一种机制（代偿），不是应该扣分的“冲突”。
- “通量是否改变”只能由对 J 的测量回答。

实测中 ρ_h 常常远离 1：酵母糖酵解各步的通量变化有 50–80% 来自代谢调控 [3]；枯草芽孢杆菌中心代谢的多数反应，酶量不足以解释通量 [4]；
酵母中代谢物浓度对反应速率的影响是酶量的两倍以上 [5]；小鼠肝脏 PEPCK 蛋白减少 90%，糖异生通量只下降约 40% [6]。通量控制通常分散在多个步骤上 [7]。

### 3.2 四条原则

1. **只在同一个估计对象内合并。** 多种示踪剂估计同一个通量，多个队列估计同一个效应，可以按精度加权。不同的量之间不相加。
2. **不同的量之间用合取与分解。** 联合结论（“通量上升且酶量上升”）的强度由最弱的一条必要证据决定（3.6）；通量已知时，用分解说明它由什么引起（3.7）。
3. **缺测、未检出、无变化、不可辨识是四种不同的状态**，各自单独表示（4.1）。“无变化”需要等效性检验，不能由不显著推出。
4. **每条规则都要有可估计的错误率**（3.11）。几条偏倚来源不同的证据指向同一结论时，结论才更可信；偏倚相同的证据重复再多也不增加可信度 [8、9]。

### 3.3 方法一览

| # | 方法 | 回答的问题 | 需要的数据 | 现有数据能否使用 |
|---|---|---|---|---|
| M1 | 四轴判定（3.4） | 每个量各自变了没有、方向如何、证据到哪一步 | 各组学差异结果 | 能 |
| M2 | 判定表（3.5） | 关于通量能说到哪一级；属于哪种机制 | M1 的结果 | 能 |
| M3 | 同一估计对象内的合并 | 多示踪剂、多队列合起来的估计与方向概率 | 分示踪剂通量；两个队列 | 已有（多示踪剂分层合并、复现模式） |
| M4 | 合取检验（3.6） | 一个联合结论成立的 p 值 | 各分量的单侧 p | 能 |
| M5 | 调控系数 ρ_h（3.7） | 通量变化中酶量贡献了多少 | 通量已确立，且有酶量 | 通量确立后才能用 |
| M6 | 节点分流比（3.8） | 分支点上碳往哪边走 | 两个出口在同一层都有读数 | 酶层、代谢物层能；通量层待 G10 |
| M7 | 交叉点分析；预先定号的模块检验（3.9） | 通路上哪一步被限制或激活；整条通路的容量是否按预期方向重排 | 相邻中间物；通路成员的预设符号 | 能，受中间物覆盖率限制 |
| M8 | 热力学驱动力；动力学拟合；由数据学习权重（3.12） | 反应离平衡多远；哪个调控因子解释通量 | 绝对浓度；8–10 个以上条件 | 暂不能 |

### 3.4 M1 四轴判定

每个（反应, 组织, 比较）在四条轴上各给出方向、状态、效应量与区间、方向概率、应测与实测的成员数。四条轴互不依赖：
任何一条轴的结论都不读取其他轴的结果。

**F 轴：通量。** 读数按对模型的依赖从少到多排列：不依赖模型的标记读数（G10 接入后）、各示踪剂的拟合通量、多示踪剂合并估计。

| 状态 | 条件（默认值，可配置） |
|---|---|
| `established` | 至少两种示踪剂可辨识且方向一致（异质性 p ≥ 0.10），合并估计通过分层置换检验与 FDR；或一个不依赖模型的读数通过检验，且拟合通量同向、Hedges' g 的绝对值 ≥ 1。同时要求较小组 ≥ 3 只、富集已达平台（有该检查时）、不是由单只动物驱动 |
| `suggestive` | 有可辨识的示踪剂，g 的绝对值 ≥ 1，逐只去除后方向保持（`direction_stability` ≥ 0.8），但达不到检验阈值（例如单一示踪剂 3 对 3 时可达的最小 p 为 0.1） |
| `conflicting` | 至少两种可辨识的示踪剂，方向不一致（异质性 p < 0.10） |
| `undetermined` | 可辨识，但效应小或不稳定 |
| `not_identifiable` | 没有任何示踪剂能约束该反应（`poorly_constrained` 或 `bound_pinned`） |
| `insufficient_n` | 所有示踪剂都有一组少于 2 只 |
| `not_measured` | 该组织或比较没有示踪数据 |

这些条件用到的字段（`identifiability`、`tracer_evidence`、`permutation_pvalue`、`min_achievable_p`、`direction_stability`、`single_subject_driven`）
都已在 MFA 差异结果中。

对现有结果包的核查（2026-11-17）另外提示两条补充规则：

- **合并前换算到相对尺度。** 不同示踪剂的绝对通量不在同一个尺度上：葡萄糖示踪的组均值有 24% 落在其他示踪剂合并值的 1.96 倍标准误之外，中位数偏高 55%。
  按绝对差值做逆方差合并会偏向数值小的示踪剂（心脏 GPS_c 两种示踪剂同向，甘油的权重为 0.997，结果被判为冲突）。
  合并应在相对同一示踪实验内参考通量的值或通量比值上进行，与 G10 中“通量效应量改用相对参考通量的差值”一并实现。
- **声明预期的主示踪剂。** 面板声明中为每个反应写明按标记位置预期能约束它的示踪剂（第 7 节的 `expected_tracers`）。
  数据判定的可辨识性与预期不符时（例如肝脏甘油激酶在甘油示踪下约束差、在葡萄糖示踪下中等），记 `tracer_expectation_mismatch`，F 轴最高为 `undetermined`。

**C 轴：酶容量。** 成员是直接酶；转运与交换反应用载体。蛋白优先，RNA 作为第二读数单独成列，不与蛋白合并。

| 状态 | 条件 |
|---|---|
| `changed` | 蛋白：按丰度加权的变化（G7）方向明确，主要成员通过 FDR，且不是只由次要同工酶驱动。没有蛋白数据时：RNA 在两个队列复现（G6），标记 `rna_only` |
| `suggestive` | 只有名义 p，或只有单个队列的 RNA |
| `unchanged` | 等效性检验通过（界限默认为 log2 倍数变化 ±0.58） |
| `mixed` | 主要成员方向相反。记为“同工酶转换”，不求和 |
| `undetermined` | 测到，但既不能判为变化，也不能判为等效 |
| `not_detected` / `not_measured` / `qc_blocked` | 映射了但数据中没有；该层没有数据；该层被 QC 判为不可评估（G2） |

复合体（如 PDH 的 PDHA1、PDHB、DLAT）按整体判断：已定量的核心亚基中至少一半同向，且至少一个通过 FDR。多个复合体共用的亚基（如 DLD）不参与方向判断。
RNA 与蛋白不一致时，C 轴取蛋白的结论，并标记 `rna_protein_discordant`，两者都展示。

**R 轴：调控状态。** 成员是带符号的调控因子（G17）：修饰酶（激酶、磷酸酶）、别构效应物（代谢物）。成员的变化乘以它对酶活性的作用符号，得到对活性的指向。
由调控因子的量得出的结论最多到 `suggestive`；有活性状态数据（磷酸化位点、测得的效应物）时才能到 `changed`。
知识表中没有调控因子的反应记为 `no_knowledge`。转录因子通过酶量起作用，只展示，不计入 R 轴的方向，以免与 C 轴重复。

**M 轴：底物与产物。** 给出底物池、产物池各自的变化，以及同一只动物内的产物/底物比（派生特征 `step_ratio`，G8）。
辅因子（ATP、NADH 等）只展示，不参与判断。近平衡反应（`direction_evidence: capacity`）以质量作用比为主要读数。
M 轴单独不能给出通量方向；它的用途是在通量方向已知或假定时，说明这一步是被上游推动、被下游拉动，还是受到限制（3.9）。

### 3.5 M2 判定表：结论等级与机制类型

**结论等级**回答“关于这个反应的通量，现在能说到哪一步”。

| 等级 | 条件 | 可以怎么写 |
|:--:|---|---|
| L3 已证实 | F 轴为 `established` | “通量上升 / 下降” |
| L2 有提示 | F 轴为 `suggestive` | “示踪数据提示上升 / 下降，待重复” |
| L1 预测 | F 轴无法判断；C 轴为 `changed` 且稳健（通过 FDR 或跨队列复现，不是由单只动物驱动，未被 QC 阻断）；M 轴与 R 轴没有相反的指向 | “酶容量上调，预测通量上升；可由某示踪剂或读数检验” |
| L0 未定 | 其余情况，包括 C 轴与 M 轴或 R 轴指向相反 | 只描述各轴的现象 |
| LX 矛盾 | F 轴为 `conflicting` | “示踪剂之间不一致”，先查模型与稳态 |

L1 另记同向的非通量轴数（1–3）和能检验它的示踪剂（取自 `tracer_coverage`：哪种示踪剂能约束这个反应）。

**机制类型**只在 L2 及以上给出，按 C 轴相对 F 轴的方向与 ρ_h 的区间归类。

| 类型 | 条件 | 含义 |
|---|---|---|
| 酶量主导 | ρ_h 的区间包含 1 | 通量变化主要由酶量解释 |
| 代谢调控主导 | C 轴为 `unchanged`，或 ρ_h 的区间在 0 附近 | 酶量没变，通量由底物、产物或修饰状态改变 |
| 两者都有 | ρ_h 在 0 与 1 之间，区间不含两端 | |
| 拮抗（代偿） | ρ_h < 0，即酶量变化与通量变化方向相反 | 酶量的变化在抵消通量的变化 [2] |
| 未解析 | C 轴未定 | |

L1 及以下不给机制类型，只记现象标签（如 `capacity_up`、`isozyme_switch`、`regulator_inhibitory_up`、`step_restrained`、`capacity_pool_opposed`）。

与加和总分相比，这张表有两点不同：酶量与通量方向相反不再降低“可信度”，而是归入一种机制；没有通量证据的反应最高只到 L1，不会因为组学层热闹而排到前面。

### 3.6 M4 合取检验

联合结论由若干个方向明确的分量组成，例如“通量上升且酶容量上升”。它的 p 值取各分量单侧 p 的最大值（合取检验，即部分合取检验中 u = n 的情形 [10]），
在（组织, 比较）内对同一类联合结论做 BH 校正。这个 p 值在分量之间有任意相关时都有效。

同时报告“至少一个分量有该方向效应”的 p 值（Simes 合并），并注明它回答的是较弱的问题。两者并列，读者能看出结论卡在哪一个分量上。

### 3.7 M5 调控系数

条件：F 轴为 L2 及以上，且 Δln J 的区间不含 0。ρ_h = Δln E / Δln J，E 取 C 轴的加权蛋白变化。区间由分层自助法得到：
在每一层、每一组内分别对动物重抽样（通量动物与组学动物不是同一批，所以各层独立重抽）。只报告区间和归类，不报告点估计的小数位。

### 3.8 M6 节点分流比

分支点上的问题通常不是每个反应各自多强，而是底物往哪个出口走。对一个节点的两个出口 a、b，在每一层定义同一个量：

| 层 | 分流比的定义 | 现有数据 |
|---|---|---|
| 通量 | ln(J_a / J_b)，逐只动物 | 有；受可辨识性限制 |
| 标记 | 两个出口的产物标记之比（如 M+3 苹果酸对 M+2 柠檬酸） | 待 G10 |
| 酶容量 | 同一只动物内两个出口的酶的 log 强度之差 | 有 |
| 代谢物 | 同一只动物内两个出口的产物之比 | 有 |

理由有三条：示踪分析中通量比值通常比绝对通量更容易辨识 [11]；同一只动物内的比值消去了上样量和组织组成这类样本级因素；
两个蛋白的强度不能直接比较绝对量，但比值在组间的**变化**是有定义的，因为蛋白特异的响应系数在组间差中抵消。
分流比作为派生特征（G8 的 `ratio` 类型，扩展到蛋白组）进入同一套差异检验。

### 3.9 M7 通路层的两种读数

**预先定号的模块检验。** 在通路面板中预先声明每个酶属于哪个出口，赋予 +1 或 −1（例如糖异生出口为正、糖酵解出口为负），
对带符号的集合做旋转检验 [12]，得到“这条通路的容量是否按声明的方向重排”的一个 p 值。符号由生物化学知识预先确定，不由数据决定；
检验的是整个集合，避免从中挑选显著的酶。旋转检验适用于每组 3 只的设计，可与 G9 共用实现。

**交叉点分析。** 沿通路的线性片段，给定通量方向（已确立的，或写明的假定），找出相邻两个代谢物变化方向相反的位置：
上游升、下游降提示两者之间的步骤受到限制；上游降、下游升提示该步骤被激活 [13]。输出为候选的调控步骤及其依据的两个代谢物。
这条定理的前提是无分支、无守恒基团的线性通路，所以结果只作为假说来源，分支点和辅因子偶联的步骤要标注。

### 3.10 排序规则（代替总分）

需要决定先验证哪个反应时，按以下键依次排序，每个键都显示在表里：

1. 结论等级（L3 > L2 > L1 > L0）；
2. 同向的非通量轴数；
3. 是否跨队列复现；
4. 完整度（实测格子数 / 应测格子数）；
5. 最高一级证据的效应量。

这是字典序，没有权重。任何两个反应的先后都可以由表中的列直接解释。

### 3.11 E8 校准

| 方法 | 做法 | 回答的问题 |
|---|---|---|
| 诱饵映射 | 用现有的诱饵模型包（E3）重跑四轴判定，统计诱饵达到“C 轴 `changed` 且 M 轴同向”和 L1 的比例，与目标对照 | 多轴一致有多少可以由随机配对产生 |
| 组别置换 | 在动物层面置换组别标签（E1），统计各等级的期望个数 | 这么多 L1、L2 里有多少是偶然的。每组 3 只时不同的置换只有 10 种，需注明分辨率 |
| 对照集 | 课题声明阳性对照（已知方向的反应或通路）与阴性对照，统计判定结果与预期的符合率 | 规则在已知生物学上是否给出正确答案 |
| 不变性测试 | 翻转通量的符号，C、R、M 轴不变；删除一整层，其余轴不变 | 各轴是否真的互不依赖 |

校准结果写入评估卡，并显示在每张矩阵图的图注中（例如“本组织的 L1：目标 12 个，诱饵期望 9 个”）。

### 3.12 M8 条件成熟后的方法

| 方法 | 条件 | 说明 |
|---|---|---|
| 热力学驱动力 | 绝对浓度 | ΔG 同时决定浓度比与正逆通量比 [14–16]；相对峰面积只能给出质量作用比的变化，给不出离平衡的距离 |
| 动力学拟合 | 同一反应在 8–10 个以上条件下都有通量、酶量、代谢物 | 检验哪个底物或调控因子解释通量 [5、17] |
| 由数据学习权重 | 至少几十个通量已确立的反应 | 在通量已确立的反应上回归 Δ通量 对 Δ酶量、Δ质量作用比，交叉验证后用作 L1 的预测模型。届时权重来自数据，而不是预设 |
| 联合拟合 | 同位素分布数据进入流程 | 平行标记实验对同一网络联合拟合 [18]；在体实验每只动物只输一种示踪剂，联合拟合只能在组水平进行 |

---

## 4. 数据展示

### 4.1 格子状态

所有图和表共用一套状态。成员（酶、调控因子、代谢物）的格子：

| 状态 | 含义 |
|---|---|
| `changed` | 通过 FDR，方向可信。是否同时达到该层入选规则的效应量阈值，另记在 `selected` 列，图上用描边的粗细区分 |
| `nominal` | 名义 p < 0.05，未过校正 |
| `no_evidence` | 测到，既未判为变化也未判为等效 |
| `equivalent` | 等效性检验通过，可以说“未改变” |
| `not_detected` | 知识表中有该成员，该层有数据，但数据中没有这个特征 |
| `not_measured` | 该来源、该组织没有这一层数据 |
| `qc_blocked` | 该层在此组织 × 比较上被 QC 判为不可评估 |
| `no_knowledge` | 知识表中这个角色没有成员 |

通量的格子另有 `not_infused`（没有输注该示踪剂）、`insufficient_n`、`not_identifiable`。

### 4.2 图 A：通路证据矩阵

每个（通路面板, 组织, 比较）一张。行按通路顺序自上而下，代谢物节点与反应交替排列，所以每个反应的上一行是它的底物，下一行是它的产物。
列分为五块；每个来源（队列）在自己的块内各占一列。面板中声明的示踪剂即使没有数据也要占一列。

示意（PC_m 与 PDH_m 两行按附录 A 的数据填写；仅示意布局，省略了队列 B 的代谢物列）：

```
panel: liver_pyruvate   tissue: liver   contrast: HFD_vs_LFD

                 | FLUX by tracer          | ENZYME        | REGULATOR     | POOL       | ASSESSMENT
row              | Glc Gly Lac Ala pooled  | P:A  R:A  R:B | P:A  R:A  R:B | A    P/S   | F  C  R  M   level
-----------------+-------------------------+---------------+---------------+------------+----------------------
(m) pyruvate     |                         |               |               | .          |
[r] PC_m         | x   .   !   /   .       | ^^   ^^   .   | -    -    -   |      /     | .  ^^ -  .   L1 (1)
(m) oxaloacetate |                         |               |               | .          |
(m) pyruvate     |                         |               |               | .          |
[r] PDH_m        | x   x   x   /   x       | .    ^^   .   | .    ^^i  ^^i |      /     | x  .  v  .   L0
(m) acetyl-CoA   |                         |               |               | ?          |

^^ / vv  升 / 降，通过 FDR        ^ / v  升 / 降，名义 p 或仅效应量      .  测到，无变化证据      =  等效
x  不可辨识      !  样本不足      /  未测或未输注      ?  映射了但未检出      -  知识表无成员      i  抑制性调控因子
```

实际绘图时的编码见 4.5。右侧的评估列显示四条轴的结论、结论等级（括号内为同向的非通量轴数）与完整度。

### 4.3 图 B：反应档案页

每个（反应, 组织）一页，展示逐只动物的数据。小样本下只画汇总柱会掩盖数据的分布，所以所有面板都画出每只动物的点 [19]。

| 行 | 面板 | 内容 |
|---|---|---|
| 页眉 | 反应式与结论 | 网络中的底物 → 产物、是否可逆、反应类别；四条轴的结论、结论等级、标记 |
| 1 | 通量，每种示踪剂一个面板 | 横轴为实验的全部组别；每只动物一个点（拟合点的中位数，误差线为四分位距）；画出约束边界；标注可辨识性与样本数；没有数据的示踪剂保留空面板并写明原因 |
| 2 | 同位素标记（G10 接入后） | 底物与产物的同位素体分布，按血浆前体富集归一化，每种示踪剂一个面板 |
| 3 | 酶的蛋白 | 每个成员一个面板，按丰度占比排序并标出占比；复合体亚基、共用亚基分别标注 |
| 4 | 酶的 RNA | 每个来源一组面板，与蛋白面板上下对齐 |
| 5 | 调控因子 | 蛋白与 RNA；面板标题写明激活或抑制 |
| 6 | 底物与产物 | 归一化峰面积（log2），每只动物一个点；另一个面板为同一只动物的产物/底物比 |
| 页脚 | 来源 | 各来源的运行 ID、被剔除或带 QC 标记的样本、重复特征的说明 |

横轴始终是实验的全部组别（如 LFD、HFD、HSD，或加上药物臂），而不只是当前比较的两组。

### 4.4 其他图

| 图 | 内容 | 用途 |
|---|---|---|
| C 示踪剂覆盖图 | 面板内的反应 × 示踪剂，颜色为可辨识性等级，标出每组样本数 | 读矩阵之前先知道哪种示踪剂能约束哪个反应 |
| D 节点分流图 | 每个分支节点一行，三层分流比的逐只动物点图并排 | 展示 3.8 的结果 |
| E 通路示意图叠加 | 在手工布局的通路图上，节点着色表示代谢物变化，边上的小图示表示通量与酶 | 汇报用。布局坐标写在面板声明中，放在后期实现 |

### 4.5 视觉编码

- 颜色表示带符号的效应量，使用色盲可辨的发散色板，零点为中性色，色阶在图内固定并截断于预设上限。
- 点的大小表示方向的确定程度（`1 − lfsr`，没有时用效应量与标准误之比）。
- 通过 FDR 的格子加描边。
- 八种状态各有固定的图示，不可辨识、未测、未检出不能画成空白，也不能画成零效应的颜色。图例必须出现在每张图上。
- 多成员的格子显示汇总值（按 3.4 的规则），角标写“通过数 / 定量数”；成员明细在档案页。
- 一张矩阵默认不超过 15 个反应，超过时按面板声明中的分段拆成多张。
- 图注自动写入：来源与运行 ID、各列的样本数、校准结果（3.11）。

### 4.6 导出表

每张矩阵图同时导出一张宽表（TSV 与 XLSX），行与图相同，列为每个成员在每一层、每个来源的效应量、p、q、样本数与状态，可直接作为补充表。
图是 `dossier_members` 的视图：同一个（反应, 层, 成员）在矩阵图、档案页和表中的数值必须逐值相同，由测试保证。

---

## 5. 建议的通用能力

### G13 反应档案：成员级证据表与逐只动物数据表

- **问题**：成员级明细写表时丢弃；分示踪剂通量与逐只动物数据没有进入整合（第 2 节第 1–5 条）。
- **设计**：
  1. 证据构建器的 `transcriptomics_evidence`、`proteomics_evidence`、`metabolomics_evidence` 除层级摘要外，返回成员行（成员名、匹配到的特征、角色、统计量、状态）。匹配逻辑只有这一处，`flux_anchored` 与新模式共用。
  2. 新的整合模式 `reaction_dossier`：读取各组学的统一差异结果、`differential_by_tracer.parquet`、`tracer_coverage.parquet`、派生特征、可靠性表、稳健性表，以及声明的其他来源；不以通量为锚。
  3. 逐只动物数据表：只抽取成员表引用到的特征，从各数据集的 `matrix.parquet` 读取，带上组别、受试者与 QC 标记。
- **接口与契约**：`integration/reaction_dossier/tables/dossier_members.parquet`、`dossier_samples.parquet`（第 6 节）。
- **验收**：合成数据中，整层缺失、成员未检出、测到无变化分别得到三种不同的状态；成员表按层聚合后与 `flux_evidence_table` 的计数一致；对同一输入，`flux_anchored` 的输出逐值不变。

### G14 通路面板：声明、证据矩阵图与档案页

- **问题**：没有带顺序的通路定义，也没有成员级的图（第 2 节第 9、10 条）。
- **设计**：
  1. 通路面板声明（第 7 节）：有序的反应清单、分支节点、预先定号的模块、应出现的示踪剂。放在课题目录，不进入包。
  2. 代谢物节点由网络化学计量（`network_reactions.tsv`）补全：相邻反应共享的代谢物成为它们之间的节点行。
  3. `omicsweave panels suggest`：按 KEGG 通路或两个代谢物之间的网络路径生成面板草稿，由人工确认顺序。
  4. 绘图模块 `integration/plotting/dossier_plots.py`：图 A–D；图与导出表共用同一个数据视图。
- **验收**：声明中引用的反应不在网络中时报错；面板声明了某种示踪剂而数据中没有时，图上出现 `not_infused` 列；图表一致性测试通过。

### G15 分轴综合评估

- **问题**：需要综合判断，但加和总分已弃用；现有证据等级以通量为锚，锚点不成立时信息很少（第 2 节第 3、12 条）。
- **设计**：四轴判定（3.4）、判定表（3.5）、合取检验（3.6）、调控系数（3.7）、排序规则（3.10）。规则与阈值集中在一个配置对象中，默认值如第 3 节。
  现有的 `evidence_grade` 词汇与 `flux_anchored` 的输出不变；新等级写在新表中（附录 B）。
- **接口与契约**：`axis_assessment.parquet`、`reaction_assessment.parquet`。
- **验收**：3.11 的不变性测试；合成数据中，通量上升而酶量下降的反应被归为“拮抗”，不出现在 L1；没有通量证据的反应最高为 L1；
  把一条轴的全部数据删除后，其余轴逐值不变。

### G16 节点分流比、预先定号的模块检验、交叉点分析

- **设计**：见 3.8、3.9。分流比作为派生特征进入现有差异引擎；模块检验使用旋转检验；交叉点分析输出候选步骤。
- **接口与契约**：`node_split.parquet`、`pathway_tests.parquet`。
- **验收**：合成的分支节点中，设定的分流变化在酶层与代谢物层被检出，两个出口同比例变化时不被检出；带符号的模块在符号全部取反后 p 值不变、方向相反。

### G17 知识层补充

- **问题**：调控因子没有符号；基因表与蛋白表的角色不一致；复合体与同工酶不分（第 2 节第 7、8 条）。
- **设计**：
  1. 调控因子条目增加 `effect`（`activates`、`inhibits`、`unknown`）与 `mode`（`phosphorylation`、`allosteric`、`transcription` 等）；别构效应物（代谢物）从 `related` 中分出，单独成一类并带符号。
  2. 基因表与蛋白表的角色由同一张表生成，消除不一致（PDH_m 的 Pdk、Pdp 归入调控因子）。
  3. 直接酶增加成员关系：`isozyme`、`complex_subunit`、`shared_subunit`。同工酶加权只在 `isozyme` 之间进行。
  4. 审计（`omicsweave refs audit`）增加：每个调控因子都有 `effect`，或明确写 `unknown`；两张表的角色一致。
- **验收**：审计在现有映射上报出 PDH_m 的角色不一致；`unknown` 的调控因子只展示，不进入 R 轴的方向。

### E8 判定规则的校准

见 3.11。结果写入 `evaluation/assessment_calibration.parquet`，并进入评估卡。

---

## 6. 数据契约

均位于 `integration/reaction_dossier/`。

**`tables/dossier_members.parquet`**：每个（反应, 组织, 比较, 来源, 层, 角色, 成员[, 示踪剂]）一行。

| 列 | 说明 |
|---|---|
| `flux_id`、`tissue`、`contrast_id`、`source` | 键；`source` 为来源标签（队列） |
| `layer` | `mfa`、`transcriptomics`、`proteomics`、`metabolomics`、`derived` |
| `role` | `flux`、`direct_enzyme`、`carrier`、`regulator`、`transcription_factor`、`downstream_enzyme`、`substrate`、`product`、`effector`、`cofactor`、`step_ratio` |
| `member_id`、`feature_id`、`feature_name` | 知识表中的成员名；数据中匹配到的特征（未检出时为空） |
| `tracer` | 通量行：示踪剂名，合并估计为 `pooled` |
| `regulator_effect`、`member_relation`、`mapping_method`、`mapping_quality` | 知识层信息 |
| `log2fc`、`effect`、`effect_type`、`se`、`ci_lower`、`ci_upper`、`pvalue`、`padj`、`lfsr` | 统计量 |
| `n_numerator`、`n_denominator`、`mean_numerator`、`mean_denominator`、`abundance_share` | 样本数、组均值、同层成员中的丰度占比 |
| `cell_state`、`state_reason`、`selected` | 4.1 的状态与原因；是否达到该层的入选规则（FDR 与效应量阈值） |
| `identifiability`、`frac_at_bound` | 通量行 |
| `qc_status`、`flags`、`replicated`、`replication_fdr` | 可靠性、稳健性（`single_subject_driven`、`minor_isozyme`、`duplicate_feature`）、复现 |
| `implied_direction`、`direction_evidence` | 该成员对通量方向的指向（酶升为升，产物升为升，底物升为降，调控因子按符号）；反应语义 |

**`tables/dossier_samples.parquet`**：`source`、`layer`、`tissue`、`feature_id`、`tracer`、`sample_id`、`subject_id`、`group`、`value`、`value_type`、`scale`、`units`、`at_bound`、`qc_flag`。只含成员表引用到的特征，以及逐只动物的派生比值。

**`tables/axis_assessment.parquet`**：`flux_id`、`tissue`、`contrast_id`、`axis`（`flux`、`capacity`、`regulation`、`pool`）、`direction`、`state`、`effect`、`ci_lower`、`ci_upper`、`sign_probability`、`n_expected`、`n_measured`、`n_sources`、`replicated`、`basis`（依据的成员）、`rule`、`flags`。

**`tables/reaction_assessment.parquet`**：`flux_id`、`tissue`、`contrast_id`、`reaction_class`、`flux_claim`（`L3`、`L2`、`L1`、`L0`、`LX`）、`claim_direction`、`n_concordant_axes`、`discordant_axes`、`mechanism_class`、`rho_h_lower`、`rho_h_upper`、`conjunction_statement`、`conjunction_pvalue`、`conjunction_padj`、`any_pvalue`、`completeness`、`completeness_by_layer`、`testable_by`、`phenomena`、`rank`。

**`tables/panel_elements.parquet`**：`panel_id`、`order`、`element_type`（`metabolite`、`reaction`）、`element_id`、`label`、`segment`、`branch_node`、`exit_label`。

**`tables/node_split.parquet`**：`panel_id`、`node`、`exit_a`、`exit_b`、`layer`、`source`、`tissue`、`contrast_id`、`log_ratio_change`、`se`、`pvalue`、`padj`、`n_numerator`、`n_denominator`。

**`tables/pathway_tests.parquet`**：`panel_id`、`test`（`signed_module`、`crossover`）、`layer`、`source`、`tissue`、`contrast_id`、`statistic`、`pvalue`、`padj`、`direction`、`members`、`detail`。

**`evaluation/assessment_calibration.parquet`**：按（组织, 比较, 等级或轴状态）给出目标数、诱饵期望数、置换期望数、对照集符合率。

**`figures/<面板>/<组织>/`**：`matrix_<比较>`、`dossier_<反应>`、`tracer_coverage`、`node_split_<节点>`，以及与矩阵图同名的导出表。

---

## 7. 配置

实验配置中启用模式并给出参数；通路面板写在课题目录的独立文件中。

```yaml
# experiment.yaml
integration:
  modes: [flux_anchored, pathway_convergence, reaction_dossier]
  parameters:
    reaction_dossier:
      panels: panels.yaml                 # 相对实验目录
      tracers: [glucose, glycerol, lactate, alanine]   # 应出现的示踪剂；没有数据的列画为“未输注”
      sources:                            # 其他来源，只读；project、experiment、run、tissue_map 沿用 replication 的来源声明
        - label: B
          project: ../dio_glp1ra_project
          experiment: all_tissues
          run: latest
          layers: [transcriptomics, metabolomics]
          tissue_map: {renal_cortex: kidney}
          contrast_map: {HFD_vs_LFD: HFD_vs_LFD}
      rules:                              # 第 3 节的阈值，以下为默认值
        flux_min_abs_g: 1.0
        heterogeneity_alpha: 0.10
        equivalence_bound: 0.58
        alpha: 0.05
      decoys: {replicates: 0, seed: 20261116}
      figures: {formats: [pdf, png], max_reactions: 15}
```

```yaml
# panels.yaml（课题内容，不进入包）
panels:
  liver_pyruvate:
    title: 肝脏丙酮酸的去向
    tissues: [liver]
    contrasts: [HFD_vs_LFD]
    segments:                             # 自上而下的阅读顺序；代谢物节点由网络补全
      - name: 进入线粒体
        reactions: [PYR_trans]
      - name: 羧化与氧化
        reactions: [PC_m, PDH_m]
      - name: 糖异生出口
        reactions: [PEPCK_c, G6Pase_c]
      - name: 糖酵解出口
        reactions: [PYK_c]
    branch_nodes:
      - node: PYR_m
        exits: {carboxylation: [PC_m], oxidation: [PDH_m]}
    expected_tracers:                     # 按标记位置预期能约束该反应的示踪剂；与数据判定不符时标记
      PC_m: [lactate, alanine]
      PDH_m: [glucose, lactate]
      G6Pase_c: [glucose, glycerol, lactate, alanine]
    signed_modules:                       # 预先定号：+1 与 −1 由通路知识决定
      gluconeogenic_routing:
        plus: [Pcx, Pck1, G6pc1]
        minus: [Pklr, Gck]
    controls:                             # 可选：用于校准的对照
      positive: []
      negative: []
```

`validate` 检查：面板引用的反应在网络中存在；`sources` 能解析到运行；`tracers` 与 MFA 数据的示踪剂名可以对应。

---

## 8. 实施路线

| 步骤 | 内容 | 依赖 | 估计工作量 |
|:--:|---|---|---|
| S0 | G17：调控因子符号、角色统一、成员关系；先覆盖已有调控因子的 24 个反应 | 无 | 2–3 天，另加策展 |
| S1 | G13 第一部分：证据构建器返回成员行；`reaction_dossier` 模式；分示踪剂通量行；格子状态 | S0 | 4–5 天 |
| S2 | G13 第二部分与 G14：逐只动物表；面板声明与解析；图 A、B、C 与导出表 | S1 | 约 1.5 周 |
| S3 | G15：四轴判定、判定表、合取检验、排序规则 | S1 | 约 1 周 |
| S4 | E8：诱饵、置换、对照集、不变性测试 | S3 | 3 天 |
| S5 | G16：分流比、模块检验、交叉点分析；图 D | S2；旋转检验与 G9 共用 | 约 1 周 |
| S6 | 示踪层读数接入档案（标记面板、F 轴的不依赖模型读数）；调控系数 | G10 与新示踪批次 | 取决于数据 |

工作量按 1 名开发者估算，仅供排期参考。S0–S4 不需要新数据，也不改变任何已有统计量。`flux_anchored` 的输出只在 S0 有一处有意的改变：
G17 第 3 项让同工酶加权只在 `isozyme` 之间进行、共用亚基不再作为方向成员（8.1）。S1–S4 不改变 `flux_anchored` 的输出。
建议 S2 完成后先用一个面板在现有课题上生成图，确认版式后再做 S3。

与已有规划的关系：L 等级可以作为结论台账（G12）的判据来源；“跨队列复现”读取复现模式（G6）的结果，同时补上证据画像尚未读取复现结果的缺口（G5 第 4 项）；
通量数据集的单位（第 2 节第 6 条）在 MFA 导出时补上。

### 8.1 实施记录

**S0（2026-11-17）**：G17 先覆盖第 7 节示例面板 `liver_pyruvate` 的 6 个反应（PYR_trans、PC_m、PDH_m、PEPCK_c、G6Pase_c、PYK_c）。

| 项目 | 实现 | 未做的部分 |
|---|---|---|
| 成员知识表 | `integration/mappings/reaction_members.yaml` 与 `reaction_members.py`：角色键 `direct_enzymes`、`carriers`、`regulators`、`effectors`、`transcription_factors`、`upstream_enzymes`、`downstream_enzymes`、`related_isozymes`、`cofactors`、`pathway_neighbors`；直接酶与转运体的 `relation`（`isozyme` 默认、`complex_subunit`、`shared_subunit`）；调控因子与效应物必须写 `effect`（`activates`、`inhibits`、`unknown`）与 `mode`。`Member.activity_sign` 只对作用于活性的方式（磷酸化、去磷酸化、别构、蛋白互作）且方向已知的成员给出 ±1，`unknown`、`transcription`、`signaling` 为 0 | 其余 85 个反应未登记 |
| 角色统一 | 两张映射表中登记反应的角色分组由成员表生成（`omicsweave refs build reaction_roles`，`--check` 只比较），只改分组，不增删成员。本次改写：PDH_m 的 Pdk1–4、Pdp1–2 由基因表的 catalytic 移到 regulatory；PEPCK_c 的 Pck2 统一为 related_isozyme；G6Pase_c 的 Pck1、Pck2、Fbp1、Fbp2 由两张表的 catalytic / upstream_regulators 移到 upstream_enzymes（蛋白表新增分类 `upstream_enzymes`、`transcription_factors`、`pathway_neighbors`） | 别构效应物没有从代谢物表的 `related` 中移出（移出会改变 `flux_anchored` 的代谢物匹配数）；反应档案以成员表的角色为准 |
| 调控方向与文献 | 方向为激活或抑制的条目逐条附 PubMed 核对过的文献：PDK 抑制、PDP 激活 PDH（PMID 24798336）；乙酰辅酶 A、NADH 抑制、丙酮酸激活 PDH（PMID 373769）；DLD 为多个复合体共用的 E3（PMID 28579060）；乙酰辅酶 A 别构调节 PC（PMID 18613815）；果糖-1,6-二磷酸激活、丙氨酸抑制丙酮酸激酶（PMID 26277545，丙氨酸的证据针对 PKM2）；MPC1 与 MPC2 组成复合体（PMID 22628558）。PC 的天冬氨酸、G6Pase_c 的胰岛素信号成员方向未核对，记 `unknown` | 肝型丙酮酸激酶受 ATP 与 PKA 磷酸化抑制，未找到核对过的来源，暂未收录 |
| 成员关系进入 `flux_anchored` | 同工酶加权排除复合体亚基与共用亚基；共用亚基只记为关联，不作为方向成员（转录组与蛋白组相同） | 复合体的整体判断（核心亚基半数同向）属于 S3 的 C 轴 |
| 审计 | `omicsweave refs audit` 增加 `reaction_roles`（未登记反应的两表角色不一致 `role_mismatch`；登记反应的 `out_of_sync`、`unannotated_member`、`effector_unresolved`）与 `regulator_effect`（`effect_missing`）；`validate` 不随项目报告这些包内知识表的检查。当前：登记的 6 个反应无问题；其余反应有 53 处角色不一致（25 个反应，如 ACONT_m 的 Idh1–3、HEX_c 的 Gckr）与 103 个没有方向的调控因子 | 这份清单是其余反应的策展工作量 |

在队列 A 上的对照（组学结果包 `run_20261112_102147_5d76`，只运行 `flux_anchored`，scratch 目录中进行）：改动前后的证据等级分布相同；
`flux_evidence_profile` 1,605 行中 4 格变化（PYR_trans 在肾脏与皮下脂肪的 `tx_weighted_capacity_change` 与 `tx_major_member` 变为空，
因为 Mpc1 与 Mpc2 是复合体亚基，不再按同工酶加权），`flux_evidence_table` 6,420 行中 1 格变化（心脏 PDH_m 转录组的 `n_implied_up` 由 4 变为 3，
Dld 不再作为方向成员）。

---

## 9. 不建议做的

| 内容 | 原因 |
|---|---|
| 任何形式的跨轴加和或加权总分 | 3.1 的恒等式；诱饵评估显示各分量目前没有区分力 |
| 把“酶或 RNA 与通量同向”当作通量的支持证据 | 同上；它属于机制归类 |
| 对缺测的层做填补，或把缺测画成无变化 | 掩盖了数据的实际覆盖 |
| 只展示显著的成员 | 这正是本方案要纠正的问题 |
| 由 L1 自动生成结论性的文字 | L1 是待检验的预测 |
| 把通路示意图的自动布局放在前期 | 布局的工作量大，而矩阵图已能完整表达；示意图放在后期，用手工坐标 |

---

## 10. 风险

| 风险 | 对策 |
|---|---|
| 矩阵信息量过大，读不过来 | 每张不超过 15 个反应；多成员格子显示汇总值加角标；明细在档案页 |
| 展示完整后，通量列大部分是不可辨识 | 这是数据的实际状态，应当如实展示；图 C 说明原因，L1 的 `testable_by` 指出哪种示踪剂能补上 |
| 调控因子符号增加知识维护量 | 走注册表框架（G7），允许 `unknown`；先覆盖课题关注的通路 |
| 规则的阈值带来新的主观性 | 阈值集中在配置中并写入运行清单；E8 给出每套阈值下的诱饵比例 |
| 新等级与现有 `evidence_grade` 并存造成混淆 | 两套词汇写在不同的表中，附录 B 给出对应关系；报告中每次只引用一套 |
| 课题内容混入通用包 | 面板、对照集、来源声明都在课题目录；包内只有引擎与合成示例 |

---

## 11. 参考文献

以下文献的题目、期刊与年份已在 PubMed 核对。

1. ter Kuile BH, Westerhoff HV. Transcriptome meets metabolome: hierarchical and metabolic regulation of the glycolytic pathway. *FEBS Lett* 2001. [doi:10.1016/s0014-5793(01)02613-8](https://doi.org/10.1016/s0014-5793(01)02613-8)
2. Rossell S, et al. Unraveling the complexity of flux regulation: a new method demonstrated for nutrient starvation in Saccharomyces cerevisiae. *Proc Natl Acad Sci U S A* 2006. [doi:10.1073/pnas.0509831103](https://doi.org/10.1073/pnas.0509831103)
3. Daran-Lapujade P, et al. The fluxes through glycolytic enzymes in Saccharomyces cerevisiae are predominantly regulated at posttranscriptional levels. *Proc Natl Acad Sci U S A* 2007. [doi:10.1073/pnas.0707476104](https://doi.org/10.1073/pnas.0707476104)
4. Chubukov V, et al. Transcriptional regulation is insufficient to explain substrate-induced flux changes in Bacillus subtilis. *Mol Syst Biol* 2013. [doi:10.1038/msb.2013.66](https://doi.org/10.1038/msb.2013.66)
5. Hackett SR, et al. Systems-level analysis of mechanisms regulating yeast metabolic flux. *Science* 2016. [doi:10.1126/science.aaf2786](https://doi.org/10.1126/science.aaf2786)
6. Burgess SC, et al. Cytosolic phosphoenolpyruvate carboxykinase does not solely control the rate of hepatic gluconeogenesis in the intact mouse liver. *Cell Metab* 2007. [doi:10.1016/j.cmet.2007.03.004](https://doi.org/10.1016/j.cmet.2007.03.004)
7. Fell DA. Metabolic control analysis: a survey of its theoretical and experimental development. *Biochem J* 1992. [doi:10.1042/bj2860313](https://doi.org/10.1042/bj2860313)
8. Lawlor DA, et al. Triangulation in aetiological epidemiology. *Int J Epidemiol* 2016. [doi:10.1093/ije/dyw314](https://doi.org/10.1093/ije/dyw314)
9. Munafò MR, Davey Smith G. Robust research needs many lines of evidence. *Nature* 2018. [doi:10.1038/d41586-018-01023-3](https://doi.org/10.1038/d41586-018-01023-3)
10. Benjamini Y, Heller R. Screening for partial conjunction hypotheses. *Biometrics* 2008. [doi:10.1111/j.1541-0420.2007.00984.x](https://doi.org/10.1111/j.1541-0420.2007.00984.x)
11. Fischer E, Sauer U. Metabolic flux profiling of Escherichia coli mutants in central carbon metabolism using GC-MS. *Eur J Biochem* 2003. [doi:10.1046/j.1432-1033.2003.03448.x](https://doi.org/10.1046/j.1432-1033.2003.03448.x)
12. Wu D, et al. ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics* 2010. [doi:10.1093/bioinformatics/btq401](https://doi.org/10.1093/bioinformatics/btq401)
13. Chance B, et al. Localization of interaction sites in multi-component transfer systems: theorems derived from analogues. *Nature* 1958. [doi:10.1038/1821190a0](https://doi.org/10.1038/1821190a0)
14. Beard DA, Qian H. Relationship between thermodynamic driving force and one-way fluxes in reversible processes. *PLoS One* 2007. [doi:10.1371/journal.pone.0000144](https://doi.org/10.1371/journal.pone.0000144)
15. Noor E, et al. Pathway thermodynamics highlights kinetic obstacles in central metabolism. *PLoS Comput Biol* 2014. [doi:10.1371/journal.pcbi.1003483](https://doi.org/10.1371/journal.pcbi.1003483)
16. Park JO, et al. Metabolite concentrations, fluxes and free energies imply efficient enzyme usage. *Nat Chem Biol* 2016. [doi:10.1038/nchembio.2077](https://doi.org/10.1038/nchembio.2077)
17. Gerosa L, et al. Pseudo-transition analysis identifies the key regulators of dynamic metabolic adaptations from steady-state data. *Cell Syst* 2015. [doi:10.1016/j.cels.2015.09.008](https://doi.org/10.1016/j.cels.2015.09.008)
18. Crown SB, Antoniewicz MR. Parallel labeling experiments and metabolic flux analysis: past, present and future methodologies. *Metab Eng* 2013. [doi:10.1016/j.ymben.2012.11.010](https://doi.org/10.1016/j.ymben.2012.11.010)
19. Weissgerber TL, et al. Beyond bar and line graphs: time for a new data presentation paradigm. *PLoS Biol* 2015. [doi:10.1371/journal.pbio.1002128](https://doi.org/10.1371/journal.pbio.1002128)

---

## 附录 A 样例：用现有结果包手工拼出的两个反应

数据取自队列 A 组学结果包 `run_20261112_102147_5d76` 与队列 B `run_20261112_113117_1269`，肝脏，HFD vs LFD。数值为 log2 倍数变化；
“A”“B”指队列。这里只用来检验所需数据都能取到，并说明各种格子状态在实际数据中的样子。
RNA 的数值取自统一差异结果表，是 apeGLM 收缩后的估计；课题报告中的基因面板用的是未收缩的估计，两者的 p 值相同、倍数变化不同
（例如队列 B 肝脏 Mpc1 收缩后 +0.17，未收缩 +0.72）。档案的成员表应同时保留两种估计。

### A.1 PC_m（网络：PYR_m + CO2 → OAC_m，不可逆）

| 块 | 成员 | 数值 | 状态 |
|---|---|---|---|
| 通量 | 葡萄糖示踪（4 对 3） | −1.06；逐只动物 HFD 1.1、51.5、110.4、1.0，LFD 174.3、81.0、1.0 | `not_identifiable`（贴边界，一半拟合点在边界上） |
| 通量 | 甘油示踪（3 对 3） | +1.02；HFD 28.3、8.8、67.3，LFD 1.0、48.1、2.3；Cohen's d 0.62 | 可辨识（中等），`no_evidence` |
| 通量 | 乳酸示踪（1 对 2） | 0.00；三只动物都在边界值 1.0 | `insufficient_n` |
| 通量 | 丙氨酸示踪 | — | `not_infused`（批次尚未导入） |
| 通量 | 合并 | +1.02，只来自甘油示踪 | `no_evidence` |
| 酶，蛋白 A | PCX（数据中记作 Pc，按 UniProt Q05920 匹配） | +0.51，q 0.047 | `changed`（未达到蛋白层入选规则的效应量阈值 0.58，`selected` 为否） |
| 酶，RNA A | Pcx | +1.19，q 0.0008 | `changed` |
| 酶，RNA B | Pcx | +0.06，p 0.27 | `no_evidence` |
| 载体，蛋白 A | MPC1、MPC2 | +0.80（q 0.028）、+0.87（q 0.027） | `changed` |
| 载体，RNA | Mpc1、Mpc2 | A：+0.61（p 0.073）、+0.71（p 0.026）；B：+0.17（p 0.046）、−0.01 | `nominal` 或 `no_evidence` |
| 调控因子 | — | 知识表中为空（已知的别构激活物乙酰辅酶 A 未收录） | `no_knowledge` |
| 底物 | 丙酮酸 | −0.21，p 0.33 | `no_evidence` |
| 产物 | 草酰乙酸 | −0.03，p 0.88 | `no_evidence` |

按第 3 节的规则：F 轴 `undetermined`；C 轴 `changed`（上调；蛋白为单队列，RNA 在队列 B 未复现）；R 轴 `no_knowledge`；M 轴没有变化的证据。
结论等级 L1，同向的非通量轴数为 1，可检验的示踪剂为甘油（现有数据中唯一能约束该反应的），以及尚未导入的乳酸与丙氨酸批次。

逐只动物的通量值说明了为什么要画点：+1.02 来自跨度从 1 到 67 的 6 个值，其中 1 个在边界上。

### A.2 PDH_m（网络：PYR_m → ACCOA_m + CO2，不可逆）

| 块 | 成员 | 数值 | 状态 |
|---|---|---|---|
| 通量 | 葡萄糖、甘油、乳酸示踪 | 全部贴边界；23 只动物中 21 只在边界值 1.0 附近 | `not_identifiable` |
| 通量 | 丙氨酸示踪 | — | `not_infused` |
| 酶，蛋白 A | PDHA1、PDHB、DLAT | +0.04、+0.19、−0.20，均不显著 | `no_evidence` |
| 酶，蛋白 A | DLD（多个复合体共用的亚基） | +0.46，q 0.042 | `changed`，不参与方向判断 |
| 酶，RNA A | Pdha1、Pdhb、Dlat、Dld | +0.90（q 0.006）、+0.16、+0.54、+0.33 | 复合体整体：`changed`（上调） |
| 酶，RNA B | 同上 | +0.17、+0.03、−0.08、+0.11 | `no_evidence` |
| 调控因子，蛋白 A | PDK1、PDK2、PDK4、PDP2 | −0.02、+0.34（p 0.066）、+0.65（p 0.35，3 对 2）、+0.26（p 0.093） | `no_evidence` |
| 调控因子，蛋白 A | PDK3、PDP1 | — | `not_detected` |
| 调控因子，RNA A | Pdk2（抑制性） | +1.18，q 0.002 | `changed` |
| 调控因子，RNA B | Pdk1（抑制性） | +1.02，q 0.022 | `changed` |
| 底物 | 丙酮酸 | −0.21，p 0.33 | `no_evidence` |
| 产物 | 乙酰辅酶 A | — | `not_detected` |

按规则：F 轴 `not_identifiable`；C 轴 `undetermined`（蛋白无变化的证据，RNA 上调只在队列 A，标记 `rna_protein_discordant`）；
R 轴 `suggestive`（抑制方向：抑制性激酶的 RNA 在两个队列上调，但不是同一个亚型，蛋白层不显著）；M 轴不完整。结论等级 L0，
现象标签 `regulator_inhibitory_up`。现有三种示踪剂都不能约束该反应，`testable_by` 为空。

这两个例子里，现行的层级表分别只给出“转录组 6 个映射、2 个有信号”和“转录组 14 个映射、2 个有信号”。PDH_m 的 2 个信号是 Pdha1 与 Pdk2，
一个是催化亚基，一个是抑制它的激酶，在只记数量的表里无法区分。

## 附录 B 与现有证据等级的关系

`flux_anchored` 的 `evidence_grade` 与本方案的结论等级回答不同的问题，两者并存，不互相替换。

| 现有等级 | 对应到本方案 |
|---|---|
| Convergent | L3 或 L2，且机制类型为“酶量主导”或“两者都有” |
| Single-layer、Exploratory | L3 或 L2；机制类型由 ρ_h 决定 |
| Conflicting | L3 或 L2，机制类型为“拮抗”；不再表示可信度低 |
| Omics-concordant | L1（同向的非通量轴数 ≥ 2），或 L0 加现象标签 |
| Omics-single-layer | L1（轴数为 1）或 L0，取决于 C 轴是否稳健 |
| Omics-discordant | L0，现象标签 `capacity_pool_opposed` 或 `rna_protein_discordant` |
| No-signal、Not-assessable | L0；不可评估的反应不进入面板 |
