# OmicsWeave 分组织代谢通路总图（图 E）：设计方案

> 文档类型：设计方案（M0a–M3 已实现，提交 `f8ea390`、`5d64c38`，记录见第 11 节；底图版式已做第一轮调整，可继续调）
> 版本：v0.4（2026-11-24，4.5 与第 11 节补充样本 QC 图面降噪、中性双斜杠 `//` 符号统一与图注 `Sample QC:` 比较级披露）；v0.3（2026-11-19，第 11 节实施记录）；v0.2（2026-11-19，第 8 节确定：自绘总图、聚合节点上图并为未评估的酶促反应补登成员、高亮 L2 及以上与 L1 空心并标前 10；
> M0 拆为 M0a 别名核对与 M0b 补登成员，4.3 增加聚合节点被逐步反应覆盖时的处理）；v0.1（2026-11-19）
> 依据：DIO 课题 `local/projects/dio_project` 的无面板整合包 `run_20261118_184922_034f` 与交付 `b6352c7c…` 的只读核查；
> 包内网络表 `integration/mappings/data/network_reactions.tsv`（模型 `dio_multi_tissue`）
> 配套文档：[反应档案与分轴评估](omicsweave_reaction_dossier_and_assessment.md)（图 E 原列于其 4.4 与第 9 节，本文接续其编号，
> 新能力记为 G18）、[候选发现与交付计划](omicsweave_discovery_and_delivery_plan.md)
> 范围：底图、数据上图的规则、整体判断的边界、配置与输出、实施路线。文中 DIO 的数字只说明图能画到什么程度，不是课题结论

---

## 1. 结论摘要

1. **适合做，而且不需要新的统计。** 图 E 是现有评估表的又一个视图：每个（组织, 比较）一页，在一张固定的中心代谢总图上，
   标出每个反应的四轴状态（F 通量、C 酶容量、R 调控、M 底物与产物）与结论等级，以及每个代谢物节点的变化。
   它只读 `dossier_members`、`axis_assessment`、`reaction_assessment` 与校准表，不改变任何等级；同一个反应在图 A、图 B 与图 E 上的数值相同，由测试保证。
2. **底图自绘，不用 KEGG 图片。** 坐标文件与网络表放在一起，节点用 KEGG、Reactome、EC 交叉引用（3.1）。
3. **覆盖不均是这张图要如实画出的内容。** DIO 每个组织有拟合通量的网络反应为 70–75 个，集中在糖酵解与糖异生、TCA、PPP、
   丙酮酸节点、几种氨基酸和酮体。网络中的脂肪酸合成与氧化、甘油脂合成与分解、尿素循环、脯氨酸、肌酸和全部核苷酸反应，
   既没有拟合通量，也没有成员知识，目前不进入评估。图上这些区域照画，标为“未评估”；它们的组学证据目前只能由聚合节点（如 `FAO_m`、
   `FAS_c`、`UREA_cycle`）给出（2.2）。已确定在 M0b 中为这些逐步反应补登成员，使它们有逐步的组学证据（通量仍为未测）。
4. **箭头按生理净方向画，变化按生理过程表述。** 网络中的可逆反应可以在某个组织里反向运行，通量的符号相对于网络的书写方向（4.1）。
5. **关键节点沿用现有判定。** L2 及以上实心粗框，L1 空心虚框，另标出排名前 10 的反应；每页印出该组织 × 比较的 L1 数与诱饵期望。
   图上没有通路级总分；各类别的整体读数放在侧栏，用通路融合结果，标明只作探索排序（4.4、第 5 节）。
6. **分 4 步实施（M0–M3）。** M0 整理知识表，分两部分：M0a 核对与网络反应重复的别名 ID（不先处理会在图上重复画、在候选表中重复列，2.3），
   是绘图的前提；M0b 为网络中未评估的 51 个酶促反应补登成员。M0 改变评估集，完成后要重跑整合与诱饵校准。可交互版本另列为 M4。

| 编号 | 能力 | 目的 |
|:--:|---|---|
| G18 | 分组织代谢通路总图：总图布局、图 E 与配套表 | 在一张图上看一个组织的整体代谢状态，以及每一步证据的完整程度 |

---

## 2. 依据：现状核查

### 2.1 可以直接复用的部分

| 需要的东西 | 现有来源 | 情况 |
|---|---|---|
| 反应、化学计量、可逆性、区室 | `network_reactions.tsv` | 156 个反应 ID（碳网络 107 行、氮网络 53 行，`ASPTA_c`、`GLNS_c`、`GPT_c`、`PHGDH_PSAT_PSP_c` 两个网络都有）；133 个物种，85 种代谢物 |
| 数据库交叉引用 | `core/idmap/reaction_pathway_mapping.py` 的 `REACTION_REGISTRY` | KEGG 反应与通路图号、Reactome、EC。评估集的 163 个 ID 中 147 个有 KEGG 图号；网络的 156 个 ID 中 73 个没有（脂肪酸、甘油脂、氮网络几乎全部） |
| 每个（反应, 组织, 比较）的成员证据 | `tables/dossier_members.parquet` | 034f 共 83,284 行：RNA 23,842、通量 21,190、代谢物 20,518、蛋白 17,680；未测、未检出等各有状态 |
| 判定 | `tables/axis_assessment.parquet`、`tables/reaction_assessment.parquet` | 163 个反应 ID × 组织 × 比较，共 4,238 行 |
| 诱饵校准 | `evaluation/assessment_calibration.parquet` | 034f 全局 L1 为 19，诱饵期望 15.0（3 次重复） |
| 代谢物节点 → 测得特征 | `dossier/view.py` 的 `species_key`、`node_members` | 不区分区室：`PYR_c` 与 `PYR_m` 取同一个特征 |
| 锁死的通量 | 通量行的 `flux_group:<首反应>` 标记 | 21,190 个通量行中 2,244 行带此标记（化学计量上锁死，与组内反应为同一次检验） |
| 净通量在两组间变号 | 通量行的 `net_flux_sign_change` 标记 | 309 行 |
| 视觉编码 | `integration/plotting/dossier_plots.py` | 发散色板（截断于 ±2）、状态符号、轴符号、点大小，可直接复用，与图 A–D 一致 |
| 交付 | `core/delivery.py` | 交付目录已有 `figures/`、`panels/`、`summary.xlsx`、`index.json` |

### 2.2 各类别的覆盖（DIO，034f）

按第 3.3 节的区域划分统计网络中的反应。“L2 / L1”为行数（反应 × 组织 × 比较）。

| 区域 | 网络反应 | 其中酶促或转运 | 有拟合通量 | 进入评估 | L2 | L1 |
|---|--:|--:|--:|--:|--:|--:|
| 糖酵解 / 糖异生 | 20 | 14 | 20 | 20 | 17 | 4 |
| 糖原 | 3 | 3 | 3 | 3 | 0 | 0 |
| 磷酸戊糖途径 | 7 | 6 | 7 | 7 | 12 | 0 |
| 丙酮酸节点与 TCA | 21 | 20 | 21 | 21 | 12 | 5 |
| 甘油-3-磷酸与甘油脂 | 13 | 11 | 5 | 5 | 10 | 2 |
| 脂肪酸 | 21 | 15 | 4 | 5 | 3 | 2 |
| 酮体 | 7 | 5 | 7 | 7 | 5 | 0 |
| 氨基酸 | 28 | 20 | 16 | 16 | 9 | 2 |
| 尿素循环与肌酸 | 15 | 8 | 0 | 1 | 0 | 0 |
| 核苷酸 | 21 | 17 | 0 | 0 | 0 | 0 |
| 合计（网络内） | 156 | 119 | 83 | 85 | 68 | 15 |

另有 4 行 L2、4 行 L1 落在网络以外的评估 ID 上：L2 是 MFA 的伪反应（`BIOMASS_REACTION_1`、`BIOMASS_REACTION_6` 两行、`SER_protein_input`），
L1 是 `FAO_m`、`MAL_trans` 与两个别名 `G6Pase`、`PK_c`（2.3）。全部 L2 的 F 轴均为 suggestive，没有 established，各自只靠一种示踪剂；
L1 全部在肝脏 HFD_vs_LFD，与诱饵期望相当。

几点说明：

- **脂肪酸**：有通量的 4 个是 `LIPID_c`、`FA_usage`、`FA_exchange`、`FA_esterification_c`，都是聚合或交换反应；`ACC_c` 只有组学。
  逐步反应（`ACSL_PAL_c`、`CPT1_c`、`CACT_PALCAR`、`CPT2_m`、`CRAT_m`、`ELOVL6_c`、`SCD_*`、`PAL_beta_oxidation`）都不在评估中。
- **甘油脂**：有通量的是 `GPS_c`、`G3PP_c`、`GK_c`、`GLYCEROL_exchange`、`lipid_mobilization_c`；`GPAT_c` 到 `MGL_c` 的 TAG 合成与分解不在评估中。
- **氨基酸**：有通量的 16 个是丙氨酸、天冬氨酸、谷氨酸、谷氨酰胺、丝氨酸的转氨、转运与交换。网络中的脯氨酸、天冬酰胺、`SHMT_c`
  不在评估中；支链氨基酸、蛋氨酸、色氨酸、组氨酸、芳香族氨基酸在网络中没有逐步反应，只有评估集里的聚合节点。
- **尿素循环、肌酸、核苷酸**：网络中逐步反应齐全，但 DIO 拟合的模型都没有这些模块；进入评估的只有 `CK_c`（只有组学）。
- **组织间的差异来自数据**：肝脏独有通量的反应为 `FBPase_c`、`G6Pase_c`、`PGMT_c`、糖原的合成与分解、`ACAC_production_m`、`BHB_output`；
  心脏独有的是 `SCOT_m`、`ACACCOA_lyase_m`、`BHB_input`。包内没有“组织 → 模块”的对应表，按 MFA 结果中是否有该组织的拟合值判断即可。
- `white_adipose`、`muscle`、`renal_cortex` 在 DIO 中只有 MFA 层（`tissue_map` 未对应到组学），这些页的组学格子为 not_measured。

### 2.3 需要先处理的问题

1. **别名 ID 重复。** 评估集的 163 个 ID 中，78 个不在网络表里。其中 12 个没有任何测得的组学成员（`*_protein_input`、`*_infuse`、
   多数 `BIOMASS_REACTION_*`）；其余 66 个在肝脏有测得的成员，分为两类：
   - **疑似网络反应的别名**：如 `PK_c` 与 `PYK_c`、`GAPDH_c` 与 `GAPD_c`、`AKGDH_m` 与 `AKGD_m`、`MDH1_c` / `MDH2_m` 与 `MDH_c` / `MDH_m`、
     `G6PDH_c` + `PGL_c` + `GND_c` 与 `G6PDH2R_PGL_GND_c`、`G3PDH_c` 与 `GPS_c`、`GLYCK_c` 与 `GK_c`、`G6Pase` 与 `G6Pase_c`、
     `GLYCOGENOLYSIS` 与 `GLYCOGENOLYSIS_c`、`ACAC_production` 与 `ACAC_production_m`、`ASPGLu_m` 与 `ASPGLU_m`、`TRIOSE_c` 与 `TPI_c`、
     `FUM_m` / `FUM1_m` 与 `FUMH_m`、`CPT1` 与 `CPT1_c`、`CPS1_m` 与 `CPS1_c`。它们来自旧的基因映射表，在 `reaction_assessment` 中与网络反应各占一行，
     同一步骤因此可以得到两个等级。034f 中已经出现：肝脏 HFD_vs_LFD 的 19 个 L1 里，`G6Pase` 与 `G6Pase_c`、`PK_c` 与 `PYK_c` 各占两个，
     所以不同的步骤最多 17 个；候选表也把它们各列两行。诱饵也在同一个评估集上计算，两者的比较口径一致，但报告中的“L1 个数”多计了。
     成员也不一致，且旧表的质量较低：肝脏 HFD_vs_LFD 中，`PYK_c` 测到 2 个 RNA 成员（`Pklr`、`Pkm`），
     `PK_c` 测到 10 个，多出的是 PKM 剪接因子（`Ptbp1`、`Hnrnpa1`、`Srsf1` 等，在旧表中记为直接酶）与转录因子（`Hif1a`、`Myc`、`Tp53`）。
   - **网络以外的步骤或通路聚合**：如 `G3PD_m`（线粒体甘油-3-磷酸脱氢酶）、`ME2_m`（线粒体苹果酸酶），以及 `FAO_m`、`FAS_c`、`GLYCOLYSIS`、
     `UREA_cycle`、`BCAA_trans`、`MET_cycle`、`TRP_metabolism`、`HIS_metabolism`、`PHE_TYR_metabolism`、`FOLATE_cycle`、`GSH_synthesis`、
     `POLYAMINE_synthesis`、`PHOSPHOLIPID`、`SPHINGOLIPID`、`BILE_acid` 等。其中不少只有代谢物成员（如 `UREA_cycle` 7 个、`BCAA_trans` 5 个），
     `FAO_m`、`FAS_c` 有蛋白和 RNA（各约 6–8 个蛋白、20 个以上 RNA）。

   别名是否就是同一步骤，需要逐个核对；核对后的处理放在 M0a（第 7 节）。绘图时只画网络 ID，不在绘图阶段合并两行的证据，
   否则同一步骤在图 E 与图 A、候选表中的数值会不一致。
2. **通量的符号。** 通量相对网络的书写方向定号。F 轴的方向是组均值之差的符号（档案文档 8.1 的 G10 部分接入记录），不是生理过程的升降。
   例如 `BHBDH_m` 写作 BHB → AcAc，在肝脏为负（葡萄糖示踪剂：HFD −84.7、LFD −131.9），差为正，但生理上是由 AcAc 生成 BHB 的净通量减少。
   上图时必须换算（4.1）。
3. **锁死的反应不是独立证据。** 无分支线性片段（如 PFK 与 FBA、TCA 中 AKGD → FUMH）的净通量由质量守恒决定，同组反应在同一次检验中，
   图上一段红色可能只是一个检验结果。
4. **证据强度。** 034f 的 21,190 个通量格子中没有 changed，nominal 2 个；F 轴 suggestive 72 行、not_identifiable 905 行、not_measured 2,632 行、
   undetermined 626 行。若按“显著”给边着色，整张图几乎全灰；若按等级高亮，肝脏的标记大多处在随机水平。所以图注必须带诱饵比，
   高亮的门槛默认放在 L2。

---

## 3. 底图

### 3.1 为什么不直接用 KEGG 图

| 问题 | 说明 |
|---|---|
| 结构不对应 | KEGG 每张图一条通路（map00010、map00020、map00030 分开），不区分胞质与线粒体。网络中的合并反应（`G6PDH2R_PGL_GND_c`、`PHGDH_PSAT_PSP_c`、`PURINE_DE_NOVO_c`）、转运与交换反应在 KEGG 图上没有对应位置，通量画不上去 |
| 不能看全局 | 用户要的是一个组织一页看到糖、TCA、PPP、脂肪酸、氨基酸之间的连接；拼多张 KEGG 图做不到 |
| 授权 | KEGG 的通路图片与 KGML 不能随软件包分发；出版使用需按 KEGG 的规定申请许可 |
| 网络覆盖之外的基因 | 只有 RNA、蛋白层的全基因视图，KEGG 图确实更全。这个需求可以另用 pathview 式的叠加图满足，不作为图 E 的底图 |

所以图 E 用自绘的总图：节点和边与网络表一一对应，图注与导出表给出每个反应的 KEGG 图号、反应号与 EC，读者可以回到 KEGG 查看。

### 3.2 布局文件

布局与网络一起属于模型知识，不属于某个课题，所以放在网络表旁边：`integration/mappings/data/network_layout.yaml`。
换一个项目而网络相同时，底图直接复用；课题可以用自己的坐标文件覆盖（第 6 节的 `layout`）。

```yaml
model_id: dio_multi_tissue
network_sha256: <network_reactions.tsv 的 SHA-256>   # 网络表改动后校验失败，提示更新坐标
canvas: {width: 260, height: 180, unit: mm}

compartments:                 # 区室边界画成圆角框
  m: {label: 线粒体, box: [80, 40, 175, 120]}

regions:                      # 区域：底色带、标题、侧栏中匹配的通路名
  glycolysis_gng:
    label: 糖酵解 / 糖异生
    kegg: [map00010]
    terms: [glycolysis / gluconeogenesis]      # 通路融合表中规范化后的名称（4.4）
  ppp: {label: 磷酸戊糖途径, kegg: [map00030], terms: [pentose phosphate pathway]}
  # …其余区域见 3.3

hide_species: [CO2, NH4_c, NH4_e]   # 辅底物不画节点，反应式仍写在图注表中

nodes:                        # 网络物种；同一代谢物在不同区室各占一个节点
  GLC_c:   {xy: [22, 12]}
  GLC6P_c: {xy: [22, 24]}
  PYR_m:   {xy: [92, 62]}
  GLC_e:   {xy: [22, 2], kind: boundary}      # _e、_stock、_out、_store 等画成图框边缘的短箭头
  ACCOA_lipid: {xy: [210, 30], kind: sink}

edges:                        # 网络反应；默认从底物节点连到产物节点
  HEX_c:    {region: glycolysis_gng}
  FBPase_c: {region: glycolysis_gng, offset: -4}   # 与 PFK_c 平行的反向旁路，横向错开
  TKT1_c:   {region: ppp, via: [[38, 34]]}         # 多底物、多产物反应的折点

aggregates:                   # 只有组学证据、不在网络中的评估 ID，画成方框，不画箭头
  FAO_m:      {region: fatty_acid, xy: [200, 70], label: 脂肪酸 β 氧化（基因集）}
  UREA_cycle: {region: urea_creatine, xy: [60, 160], label: 尿素循环（代谢物集）}

exclude:                      # 评估集中不上图的 ID，必须写原因
  PK_c: {reason: alias, of: PYK_c}
  BIOMASS_REACTION_2: {reason: pseudo}
  GLYCEROL_infuse: {reason: tracer_input}
```

校验（M1 的测试，运行时也检查）：

- 网络表的每个反应 ID 在 `edges` 或 `exclude` 中；交换与 `pseudo` 反应可以排除，但要写原因；
- 每条边用到的物种都有坐标（`hide_species` 除外）；
- 某次整合包中有测得组学成员的评估 ID，必须出现在 `edges`、`aggregates` 或 `exclude` 中。运行时遇到不在其中的 ID，不静默丢弃：
  写入 `pathway_map_coverage.tsv` 并在图注中计数；
- `network_sha256` 与包内网络表一致，否则报错并说明需要更新坐标。

工作量的主体是手工排坐标：133 个物种中约 40 个是边界物种（`_e`、`_stock`、`_out` 等），余下约 90 个节点需要定位。

### 3.3 区域划分

总图的大致版面：左侧竖列为糖酵解主干，PPP、丝氨酸、甘油-3-磷酸作为旁支；中央是线粒体框，内含丙酮酸入口与 TCA 环；
右侧为脂肪酸合成（胞质）、肉碱穿梭与 β 氧化（跨膜）、酮体；下方为氨基酸、尿素循环与肌酸、核苷酸；交换反应画成伸向图框外（血液侧）的短箭头。

| 区域 | 网络反应（示例） | 侧栏匹配的通路名（KEGG） |
|---|---|---|
| 糖酵解 / 糖异生 | HEX→PGI→PFK→FBA→TPI→GAPD→PGK→PGM→ENO→PYK、LDH；FBPase、G6Pase、PEPCK | glycolysis / gluconeogenesis |
| 糖原 | PGMT、糖原合成与分解 | starch and sucrose metabolism（如有） |
| 磷酸戊糖途径 | G6PDH2R_PGL_GND（氧化段）、RPI、RPE、TKT1、TKT2、TALA | pentose phosphate pathway |
| 丙酮酸节点与 TCA | PYR_trans、PDH、PC、CS、ACONT、ICDH、AKGD、SUCOAS、SUCD、FUMH、MDH，ME、ACITL，各穿梭体，乙酸支路 | citrate cycle (tca cycle)、pyruvate metabolism |
| 甘油-3-磷酸与甘油脂 | GPS、G3PP、GK；GPAT→AGPAT→LIPIN→DGAT；ATGL→HSL→MGL | glycerolipid metabolism |
| 脂肪酸 | ACC、LIPID（合并的从头合成）、ELOVL6、SCD；ACSL、CPT1、CACT、CPT2、β 氧化 | fatty acid biosynthesis / elongation / degradation |
| 酮体 | ACAC_production、BHBDH、BHB 转运与交换、SCOT、ACACCOA_lyase | synthesis and degradation of ketone bodies（如有） |
| 氨基酸 | GPT、ASPTA、GLUD、GLND、GLNS、PHGDH_PSAT_PSP、SHMT、ASNS、ASPG、脯氨酸；聚合节点 BCAA、MET、TRP 等 | alanine, aspartate and glutamate metabolism；glycine, serine and threonine metabolism；arginine and proline metabolism；valine, leucine and isoleucine degradation |
| 尿素循环与肌酸 | CPS1、OTC、ASS_ASL、ARGN；GATM、GAMT、CK | arginine biosynthesis |
| 核苷酸 | 嘌呤与嘧啶从头合成、ADSS_ADSL、IMPDH_GMPS、CTPS、AMPD、ADA、NT5、PNP、CDA、RNR、降解 | purine metabolism、pyrimidine metabolism |

区域划分只用于底色、标题和侧栏，不参与任何判定。一个反应只属于一个区域（如 PC_m 画在 TCA 区，不在糖异生区重复画）。

### 3.4 组织之间的差异

每页只画一个组织。反应的画法按该组织的数据分三种：

| 画法 | 条件 |
|---|---|
| 实线 | 该组织至少有一种示踪剂的拟合通量（含不可辨识、样本不足：这些是“有模型、约束不住”，用状态符号表示） |
| 虚线 | 该组织的 MFA 结果中没有这一步，但反应有组学成员（例如心脏的 `G6Pase_c`：模型没有这一步，`G6pc` 的 RNA 仍可能测到） |
| 细灰线 | 没有进入评估：既没有通量，也没有成员知识（DIO 的尿素循环、核苷酸即如此） |

“该组织的模型没有这一步”与“有这一步但没测到”必须区分，后者用状态符号，不画成细灰线（档案文档 4.1 的规则）。

---

## 4. 数据如何上图

### 4.1 反应边

**方向。** 箭头画在生理净方向上：

| 情况 | 箭头 |
|---|---|
| 两组的组均值同号 | 按该符号的方向（为负时从网络的产物指向底物） |
| 两组变号（`net_flux_sign_change`） | 双向箭头，加“⇄”标记，图注说明 |
| 没有通量 | 按网络的书写方向；可逆反应画成细的双向箭头 |

**颜色。** 边的颜色只表示 F 轴，并换算到生理方向：生理变化 = F 轴方向 × 组均值的符号。红为生理方向上的通量增加，蓝为减少，
深浅为效应量（与图 A 同一色阶）。F 轴不是 suggestive 或 changed 时边为灰色，并在边上标状态符号（x 不可辨识、! 样本不足、/ 未测或未输注）。
以 2.3 的 `BHBDH_m` 为例：组均值为负，箭头从 AcAc 指向 BHB；差为正，F 轴若达到 suggestive，换算后为蓝色，读作 BHB 生成减少。
034f 中它在肝脏 HFD_vs_LFD 的 F 轴为 undetermined，所以实际画成带“·”的灰色边。

**四轴小条。** 每条边的中点放一个 4 格小条 `[F C R M]`，沿用图 A 评估块的符号：实心 ▲▼ 为 changed（F 轴为 established），空心为 suggestive，
= unchanged，≠ mixed，· undetermined，以及各缺测状态。F、M 两轴的方向换算到生理方向；C、R 两轴表示酶容量与活性，不随方向换号，图例写明。

**等级与排名。** 小条右侧标结论等级；L2、L3 的边加粗黑框并写出反应 ID，L1 加空心虚框；排名前 10 的反应标出名次，
与 `candidates/reactions.tsv` 的 `rank` 一致。

**锁死的反应。** 同一 `flux_group` 的边加同一个上标（如 g3），图例写“同组反应为同一次检验，不算独立证据”。

### 4.2 代谢物节点

- 填色为代谢组的 log2FC，取值规则与图 A 的池格子相同（`node_members` 选出的特征），所以两张图的数值一致；通过 FDR 的节点加描边。
- 代谢组测的是全组织，`PYR_c` 与 `PYR_m` 共用同一个值，图例写明；不把它画成两个独立的测量。
- 状态：? 映射了但未检出、/ 该组织没有代谢组、– 没有可映射的特征。缺测不画成中性灰（中性灰表示测到且无变化）。
- 边界物种（血液侧）只画短箭头，不填色；血清的测量属于 serum 页。

### 4.3 聚合节点

不在网络中的评估 ID（如 `FAO_m`、`FAS_c`、`UREA_cycle`）画成所在区域旁的方框，内含 P、R、M 三个小格，取值规则同图 A 的酶格子
（同工酶按丰度占比加权，没有占比时取中位数），角标为“通过 FDR 数 / 定量数”，另附 C 轴状态。方框不画箭头，不表示通量。
别名 ID 不画（2.3）。

M0b 为逐步反应补登成员后，聚合节点与逐步反应可能共用同一批基因（例如 `FAO_m` 的酶登记到 `CPT1_c`、`PAL_beta_oxidation` 等步骤之后）。
成员已全部被逐步反应覆盖的聚合节点改记为 `exclude`（原因 `covered_by_steps`），不再画方框，避免同一批基因在图上出现两次；
只被部分覆盖的保留方框，配套表中列出它与逐步反应共有的成员。是否覆盖由 M1 的校验按成员清单自动判断，结果写入 `pathway_map_coverage.tsv`。

### 4.4 类别侧栏

图的右侧为每个区域一行，内容全部取自已有的表，不做汇总计分：

- **覆盖**：该区域画出的反应数、该组织有拟合通量的数目、可辨识的数目、L2 与 L1 的数目；这些是完整度的描述，不是评分；
- **通路融合**：按区域声明的 `terms` 从 `pathway_fusion.parquet` 取该组织 × 比较的行，列出各层 padj 与方向、融合 FDR、方向是否一致；
  标明“按名称匹配；融合值仅作探索排序”（交付计划 P2 的结论）；
- **预先定号的模块检验**：只有在面板中声明了 `signed_modules` 且属于该组织时才显示（G16），注明面板名。

### 4.5 图注

每页自动写入：组学与整合的运行 ID；该组织 × 比较各层、各示踪剂的样本数；当该（组织, 比较）存在非 `ok` 的样本可靠性状态（如 `downweighted`、`qc_confounded`、`qc_insufficient`）时，在图注首行自动附加纯英文比较级 QC 披露（如 `Sample QC: RNA 3v2 usable of 4v2 (downweighted)`），避免在全图每一张反应卡片上重复盖感叹号造成视觉噪声；该组织 × 比较的 L1 数、诱饵期望与标准差（没有运行校准时写“未做诱饵校准”）；
`panels.yaml` 的 `selection`（pre_specified / post_hoc，没有面板时不写）；未上图的评估 ID 数（3.2）；布局文件的路径与校验值。

### 4.6 配套表

每页同时导出 `pathway_map` 表：每个画出的元素（反应边、代谢物节点、聚合节点）一行，列出组织、比较、元素 ID、区域、KEGG 图号与反应号、EC、
画法（实线、虚线、细灰线）、生理方向、四轴状态与方向、等级、名次、颜色所用的值与状态、`flux_group`。
同一个（反应, 组织, 比较）的值必须与图 A 的格子逐值相同，由测试保证（与档案文档 4.6 的要求一致）。

---

## 5. 整体判断能说到哪一步

这张图用来回答“这个组织里，哪些步骤的证据指向变化，哪些步骤约束不住，各层是否一致”。它不能单独回答“某条通路整体上调”：

| 想说的话 | 能否从图 E 读出 | 依据 |
|---|---|---|
| 某一步的通量在生理方向上增加 | 能，限于 F 轴为 suggestive 或 established 的边 | F 轴；L2 的边全部只靠一种示踪剂，应写明是哪一种 |
| 某一步的酶量与通量同向或拮抗 | 能，看小条的 C 格与 F 格 | 机制类型（档案文档 3.5），属于解释，不是通量的证据 |
| 糖酵解整体上调 | 不能只看红边多少 | 线性片段上的反应锁死在同一次检验里；整体判断需要预先声明的模块检验（G16），或把通路融合当作探索性线索 |
| 尿素循环、核苷酸代谢的通量变化 | 不能 | DIO 没有这些模块的拟合通量，只能说酶或代谢物的变化 |
| 肝脏的 L1 是可信的预测 | 要先看图注 | 034f 肝脏 L1 19 个，诱饵期望 15.0，与随机水平相当 |

建议的读法：先看图注的诱饵比与样本数，再看哪些区域是实线（有通量），然后看粗框的关键节点，最后到图 B 档案页看逐只动物的数据。

---

## 6. 配置与输出

在 `reaction_dossier` 的参数中增加一组 `pathway_map`（与现有的 `figures`、`calibration` 同级，未知键报错）：

```yaml
integration:
  parameters:
    reaction_dossier:
      pathway_map:
        enabled: true
        layout: ""                     # 空 = 包内 network_layout.yaml；可写课题目录中的坐标文件
        tissues: []                    # 空 = 有数据的全部组织；DIO 建议排除 serum（只有 4 个有通量的网络反应）
        contrasts: []                  # 空 = 全部比较
        highlight: {min_claim: L2, show_l1: hollow, top_rank: 10}
        formats: [pdf, svg, png]
```

输出：

| 位置 | 内容 |
|---|---|
| 整合包 `integration/reaction_dossier/figures/pathway_map/<组织>/map_<比较>.*` | 图 E |
| 整合包 `tables/pathway_map.parquet`、`tables/pathway_map_coverage.tsv` | 配套表；未上图的评估 ID 及原因 |
| 交付目录 `pathway_map/<组织>/` | 图与表的副本；`summary.xlsx` 增加“通路总图”页，`index.json` 登记每页 |

绘图只读整合包中的表与布局文件，所以另提供一个从已有整合包重画而不重算的入口，如 `omicsweave map --source <整合运行>`
（名称在 M3 时定，用法与现有的 `stats --source` 一致）。
改坐标或改高亮规则只需重画；改知识表（M0）会改变评估集，需要重跑整合与诱饵校准。

DIO 按上面的默认值（排除 serum）为 23 页（9 个组织，各 2–3 个比较）。

---

## 7. 实施路线

| 步骤 | 内容 | 依赖 | 估计工作量 |
|:--:|---|---|---|
| M0a | 别名核对：逐个核对 2.3 的别名 ID，确认为同一步骤的记为别名，别名独有的成员逐个核对角色后再决定去向（如 `PK_c` 多出的剪接因子与转录因子不应作为 `PYK_c` 的直接酶，需要时按 G17 记为调控因子），然后从映射中移除别名，避免候选表重复列出同一步骤；不是同一步骤的（如 `G3PD_m`、`ME2_m`）保留为聚合节点 | 无 | 1–2 天 |
| M0b | 补登成员：为网络中未评估的 51 个酶促反应登记直接酶（同工酶、复合体亚基按 G17 的 `relation` 标注）与 KEGG 反应号、通路图号、EC，涉及脂肪酸合成与氧化、甘油脂、尿素循环、脯氨酸、天冬酰胺、肌酸、核苷酸；调控因子只在有核对过的文献时登记，否则不写或记 `unknown`。完成后这些反应进入评估：F 轴为 not_measured，C、R、M 轴照常判定，最高到 L1。登记前后各列出一次聚合节点与逐步反应共有的成员（4.3） | M0a（避免把成员登记到将被移除的别名上） | 2–3 天，另加策展 |
| M1 | 布局文件与校验：`network_layout.yaml`、3.2 的校验测试与 4.3 的覆盖判断；先出一版不带数据的底图请用户确认版式 | M0a 的别名清单 | 2–3 天，另加一轮版式确认 |
| M2 | 绘图：`map_view`（只生成 4.6 的表）与 `plot_pathway_map`（只画这张表）；图例、侧栏、图注；与图 A 逐值一致的测试 | M1 | 4–5 天 |
| M3 | 接入：`pathway_map` 参数组、重画命令、交付目录的副本与摘要页、`--dry-run` 检查布局文件 | M2 | 1–2 天 |
| M4 | 可选：可交互的 HTML 版（悬停显示成员明细与逐只动物小图），或导出 Escher 格式的 JSON | M3 | 以后再定 |

工作量按 1 名开发者估算。M0b 与 M1、M2 互不依赖，可以并行。M0a、M0b 都改变评估集（别名行消失，51 个反应加入），
两步完成后重跑一次整合与诱饵校准，并在档案文档 8.1 记录前后的 L 等级与诱饵期望；DIO 的正式图 E 用重跑后的整合包生成。
M1–M3 不改变任何已有统计量。

---

## 8. 已确定事项（2026-11-19）

| 事项 | 决定 | 未采纳的选择 |
|---|---|---|
| 底图 | 自绘总图（第 3 节） | KEGG 图片叠加：不分区室，画不了合并反应与通量，且不能随包分发；自绘为主、另加 KEGG 补充图的方案也未采纳 |
| 只有组学的聚合节点（`FAO_m`、`FAS_c`、`UREA_cycle` 等） | 上图，画成方框（4.3）；同时做 M0b，为未评估的 51 个酶促反应补登成员 | 上图但暂不补成员；不上图 |
| 关键节点的高亮门槛 | L2 及以上实心粗框，L1 空心虚框，另标排名前 10，每页印诱饵比（4.1、4.5），即第 6 节 `highlight` 的默认值 | 只标 L2 及以上；L2 及以上加前 10 名、不标 L1 |

---

## 9. 不做的

| 内容 | 原因 |
|---|---|
| 通路级或区域级的总分、按颜色多少给通路定性 | 档案文档第 9 节；锁死的反应会让颜色数量失真 |
| 自动布局 | 自动布局难以保持代谢图的惯用版面（糖酵解竖列、TCA 成环），而布局只需做一次 |
| 在绘图时合并别名行的证据 | 会让图 E 与图 A、候选表的数值不一致；合并属于知识表（M0） |
| 把缺测画成无变化 | 档案文档 4.1 |
| 用 KEGG 图片做底图并随包分发 | 3.1 |

---

## 10. 风险

| 风险 | 对策 |
|---|---|
| 一页信息过多，读不过来 | 四轴小条只有 4 格；成员明细不上图，在图 B 与配套表中；侧栏与图注分开放 |
| 读者把红边多当作通路上调 | 锁死反应加组标记；第 5 节的读法写进图注的固定说明 |
| 网络表更新后坐标失效 | `network_sha256` 校验，失败时报错并列出新增或改名的反应 |
| M0b 加入 51 个只有组学的反应后，L1 数与诱饵期望都会变，和此前引用的数字不可直接比较 | 整合与校准在 M0 前后各跑一次，两组数字在实施记录中并列；报告引用时写明所用的整合包 |
| 别名核对出错，丢掉真实的成员，或把旧表中角色错误的成员并进来 | M0 逐个核对并留记录；合并前后对比每个网络 ID 的成员清单与角色，差异写入实施记录 |
| 只有 MFA 层的组织页面大面积 not_measured | 这是数据的实际状态；图注写明该组织未对应到组学（`tissue_map`） |

---

## 11. 实施记录（2026-11-19）

M0a–M3 已实现：M0b 的知识表为 `f8ea390`，其余（M0a、M1–M3 与第一轮版式调整）为 `5d64c38`。与设计不同之处在各步的
“偏差”中写明。

**M0a 别名核对。** 新增 `integration/mappings/reaction_aliases.yaml`（登记）与 `reaction_aliases.py`（读取、校验与评估集）。
评估集中 66 个不在网络表中、有测得成员的 ID 逐个核对：31 个登记为别名（同一转化、网络中合并反应的拆分部分或按方向拆开的交换反应，
如 `PK_c` → `PYK_c`、`G6PDH_c` / `PGL_c` / `GND_c` → `G6PDH2R_PGL_GND_c`、`BHB_exchange` → `BHB_output` + `BHB_input`），
其余 35 个保留（33 个聚合条目，另 2 个是 MFA 的生物量伪反应）。别名中网络 ID 没有的直接酶与载体都属于其他步骤、调控因子或无关基因（如 PKM 剪接因子、抗氧化酶），
没有需要并入的成员；核对结论写在各条的 `note` 中。
偏差：设计写的是“从映射中移除别名”，实现改为登记别名、在 `reaction_dossier` 的评估集中排除（`assessment_universe`），
映射表条目保留。原因是 `flux_anchored` 以 MFA 结果的反应 ID 为锚，旧模型命名的项目（MFA 结果用 `PK_c` 等）仍需要这些条目；
MFA 结果用到的别名 ID 因此也保留在评估集中。运行备注列出被排除的别名，参数记录写入 `excluded_aliases`。

**M0b 补登成员。** 网络中没有任何映射条目的酶促与转运反应共 53 个；登记 50 个，其余 3 个不登记：`GLND_c`、`GLUD_c`
与碳网络的 `GLND_m`、`GLUD_m` 是同一种酶（重复登记会让同一批基因在两个反应中各算一次），`CRN_formation_c` 为非酶反应。
直接酶按 EC 号与 KEGG REST（`link/mmu/ec:…`）逐个对照：KEGG 多出的基因来自同一 EC 号下的其他步骤或底物
（如 CPT1 与 CPT2 共用 2.3.1.21、各 ELOVL 共用 2.3.1.199、PLPP 与 lipin 共用 3.1.3.4），按反应拆分后不收；
`Pnp2` 不在 KEGG 的 2.4.2.1 中，未收。符号以包内小鼠基因登记表为准（`Adss1`、`Adss2`、`Ctps1`、`Nt5c3`、`Pycr3`、`Gyg1`）。
三张映射表的新条目按 PyYAML 默认格式写出（`reaction_members` 的往返测试要求），原有行不变；`REACTION_REGISTRY` 增加 50 条
KEGG 反应号、通路图号与 EC 号（合并反应没有单一反应号）。代谢物条目 37 个，取网络组织侧名称中能解析到标准 ID 的部分，
另为 `PNP_INO_c`、`PNP_GUO_c`、`URI_degradation_c` 补上产物（次黄嘌呤、鸟嘌呤、尿嘧啶）。
限制：代谢物 ID 表不认识酰基肉碱、硬脂酸、油酸、棕榈油酸、胍基乙酸、吡咯啉-5-羧酸，相应步骤没有 M 轴成员；
补上这些名称属于代谢物 ID 表的重建（尚未进行）。调控因子没有登记（设计中只在有核对过的文献时登记）。

**M1 布局。** `integration/mappings/data/network_layout.yaml`（画布 120 × 102；134 个反应画成边；133 个物种中 113 个上图：
代谢物 89、血液侧边界 19、库存 5，另 3 个辅底物隐藏，其余 17 个只出现在排除的伪反应中），
`network_layout.py` 读取与校验。偏差与补充：节点可有多个位置（克隆，`at` 列表），反应连到离反应中心最近的一个，
用于 `GLU_c`、`ASP_c`、`GLN_c` 等枢纽代谢物和线粒体膜两侧的穿梭；边增加 `strip`（四轴小条相对中心的位置）。
排除 36 个 ID：MFA 伪反应 32 个（示踪剂输入、库存、蛋白分解输入、生物量）、血清胞外反应 2 个、与碳网络重复的 2 个。
聚合条目 33 个画成方框，其中 18 个（支链氨基酸、蛋氨酸循环、肌醇等）在网络以外，集中画在顶部格子中。
校验值取 LF 换行后的 SHA-256，避免 Windows 检出（CRLF）与仓库内容不一致。底图版式尚待用户确认（见下）。

**M2 视图与绘图。** `integration/dossier/map_view.py`（`map_view`、`coverage_table`、`region_summary`、`physiological_sign`）
与 `integration/plotting/pathway_map_plots.py`。代谢物节点调用图 A 的 `node_members` 与 `summarize`，测试核对两者逐值相同。
生理方向取合并估计的组均值符号（没有时取各示踪剂多数），两组变号时双向箭头加“↔”。聚合节点是否被覆盖按已知成员
（排除 `no_knowledge` 的占位行）判断；DIO 的 034f 中被覆盖的是 `ACAC_trans`、`GLYCOLYSIS`、`MAL_trans`、`TRIOSE_c`。
锁死的通量组按 MFA 行的 `flux_group` 标记在反应名后用方括号标出（组名为组内首个反应）。

**M3 接入。** 新模块 `integration/pathway_map.py`（`render_maps`、`redraw`、参数组 `settings`）；`reaction_dossier` 的参数组
`pathway_map`（默认打开）；课题的布局文件登记为整合输入。`omicsweave map --experiment E [--source …]` 从已有整合包重画到
`<结果根>/<实验>/_redraws/<运行 ID>/pathway_map_<时间>/`，不改动原结果包。交付目录复制每页的 PDF（没有时为 PNG）与视图表到
`pathway_map/<组织>/`，`summary.xlsx` 增加「通路总图」。`pyproject.toml` 的包数据增加 `mappings/data/*.yaml`。

**DIO 验证（只在临时副本中运行，未改动课题目录）。** 用正式组学包 `run_20261118_172720_cc0e` 在课题配置的临时副本中
重跑 `pathway_convergence` 与 `reaction_dossier`（与 034f 相同：不读面板，诱饵 3 次），整合包 `run_20261119_131530_9036`，
用时 585 秒（034f 约 6 分钟，图 E 26 页占其余时间；图件共 55 MB，交付只复制 PDF 与视图表）。

| | 034f（M0 之前） | 9036（M0a、M0b 之后） |
|---|--:|--:|
| 评估集反应 ID | 163 | 181（去掉 31 个别名，加入 49 个新登记的网络反应；`GLYCOGEN_SYNTHESIS_c` 原已因 MFA 结果在集内） |
| L2 | 72 | 72 |
| L1（全部在肝脏 HFD_vs_LFD） | 19 | 27 |
| 肝脏 HFD_vs_LFD 的 L1 诱饵期望 | 15.0 | 20.7（标准差 3.8） |

两包共有的（反应, 组织, 比较）中没有一行改变等级。肝脏 L1 的变化为：去掉别名 `G6Pase`、`PK_c` 两行（它们与 `G6Pase_c`、
`PYK_c` 重复），加入 10 个新登记的反应（`ACSL_PAL_c`、`CPT1_c`、`CACT_PALCAR`、`CPT2_m`、`DGAT_c`、`MGL_c`、`PAL_TG_hydrolysis_c`、
`OTC_c`、`GAMT_c`、`PURINE_DE_NOVO_c`），都由蛋白层的 C 轴给出。L1 与诱饵期望之比仍约为 1.3，肝脏的 L1 依旧应当作待检验的预测。
总图上评估集的 ID 全部有去向：画成边 129、方框 29、被逐步反应覆盖 4（`ACAC_trans`、`GLYCOLYSIS`、`MAL_trans`、`TRIOSE_c`）、
伪反应 17、血清反应 2，没有未上图的 ID。正式结果需要用户在课题目录中重跑整合（本次没有面板图；课题的面板会照常出图）。

**底图版式的第一轮调整（2026-11-19，用户要求先调版式）。** 第一版的局部重叠集中在四处：小条压在节点与名称上、
PPP 非氧化段的交叉线、核苷酸与氨基酸区的克隆连线、效应小的 L2 边颜色接近中性。做了以下改动：

- 小条自动避让（绘图）：布局不写 `strip` 时，在反应中心四周 12 个候选位置中取与节点及名称、边界与库存、方框、区域标题、
  已放好的小条和边的重叠最少的一个；L2 及以上、L1、排名靠前的先放。小条离中心较远时画细引线，中心画一个小圆点；
  反应名带浅色底，压在线上也能读；
- 次要底物 / 产物（布局 `secondary`）：转氨中的谷氨酸与 α-酮戊二酸、合成步骤中的谷氨酰胺 / 天冬氨酸 / 甘氨酸、ELOVL6 的
  丙二酸单酰辅酶 A 不画线，在中心下方写成 “GLN, ASP → GLU”。用于 GPT、ASPTA（胞质与线粒体）、ASNS、OAT、ASS_ASL、GATM、
  ELOVL6、嘌呤与嘧啶从头合成、IMPDH_GMPS、ADSS_ADSL、CTPS 共 13 条边；去掉了因此不再需要的 6 个克隆位置
（甘氨酸 2、天冬氨酸 2、谷氨酸 1、谷氨酰胺 1）；
- PPP 非氧化段改为两行：核糖-5-磷酸、景天庚酮糖-7-磷酸、赤藓糖-4-磷酸在上，木酮糖-5-磷酸、甘油醛-3-磷酸、果糖-6-磷酸在下，
  TKT1 与 TALA 各成一个 X，TKT2 在下方单独连到赤藓糖-4-磷酸的第二个位置；网络以外的方框改为 6 × 3 的格子，右移让出位置；
- 嘧啶碱出口分成两个位置，脱氧胞苷降解不再横穿 CMP；血液侧丙酮酸移到丙酮酸左上；
- 锁死的通量组只在 F 轴有指向时标在反应名后（视图表照常保留全部）；F 轴有指向的边着色时 |log2FC| 至少取 0.8，
  避免效应小的 L2 边看起来像未评估（数值以视图表为准，图例写明）；
- 校验增加一项：节点的每个位置都必须有反应连到（防止调整后留下孤立的克隆）。

仍可再调的地方：甘油脂区较挤（MGL、HSL、DGAT 的小条只能放到区块外侧，靠引线对应）、TKT2 的木酮糖连线较长、
少数反应名压在线上。坐标都在布局文件中，调整后用 `omicsweave map` 重画即可，不需要重跑整合。

**测试。** `tests/integration/test_pathway_map.py`（11 项：别名与评估集、新成员、布局校验与克隆、参数、生理方向、
覆盖判断、整合出图与视图表、与图 A 池格子逐值相同、布局错误中止运行、交付复制与摘要表、`omicsweave map` 重画的视图表与原页
相同且不改动原结果包；版式调整后另加次要物种与孤立克隆两项，共 13 项）。全套测试 573 项通过。合成数据的整合测试默认关闭出图
（`tests/integration/test_flux_anchored.py` 的 `set_integration`），否则每个测试多 20 秒左右。

**样本 QC 降噪与低样本量符号统一（2026-11-24）。** 针对样本 QC 被判定为 `qc_blocked` 时在静态图与交互式 HTML 的每一条反应边、四轴小条与聚合方框上重复叠加高对比度琥珀色感叹号三角（`⚠`）导致全图视觉噪声过重的问题，在 `visualization/primitives.py`、`integration/plotting/pathway_map_plots.py`、`pathway_map_html.py`、`dossier_plots.py` 与 `integration/pathway_map.py` 中落实两项改动：

- 移除各反应卡片、四轴小条、聚合节点与交互式 HTML 上针对 `qc_blocked` / `insufficient_n` 的琥珀色 `⚠` 警示块，统一改为浅灰缺测底色（`missing_fill`）配合低对比度中性双斜杠 `//`（静态与 HTML 图例统一标注为 `Low-n (n<2)`；保留空心菱形 `◇` 专用于 HTML 逐只动物散点图的 `bound-pinned` 贴边标记，避免符号冲突）；
- 对于仍保有有效差异统计量但存在样本 QC 风险的比较，图面正常渲染真实效应量与软封顶后等级，同时由 `pathway_map.py` 的 `_qc_disclosure()` 在每页图注首行自动附加一处纯英文/ASCII 比较级汇总（`Sample QC: RNA XvY usable of AvB (status)`），满足英文出版字体兼容与单一位置风险披露要求。
