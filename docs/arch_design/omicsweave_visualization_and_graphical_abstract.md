# OmicsWeave 统一绘图体系与 Graphical Abstract 子系统设计方案

> 文档类型：架构与设计方案（v1.0，2026-11-20）
> 适用范围：`omicsweave.visualization` 通用绘图包、反应档案图 A–E（`integration/plotting`）、单组学图件主题收敛，以及独立解耦的 Graphical Abstract 子系统。
> 配套文档：[分组织代谢通路总图（图 E）](omicsweave_pathway_map.md)、[反应档案与分轴评估](omicsweave_reaction_dossier_and_assessment.md)、[v2 架构设计](../omicsweave_arch_design_v2.md)
> 海报与组合图版：[figure_draw 海报衍生绘图标准](visualization_artwork/figure_draw_style.md)，定量图可显式使用 `poster_editorial`。

---

## 1. 设计原则与决策边界

1. **学术出版优先（`nature_editorial` 作为全库默认主风格）**
   - 所有定量数据图（火山图、证据矩阵图 A、逐只动物档案图 B、示踪剂覆盖图 C、分支分流图 D、代谢通路图 E、跨组学汇总图）统一采用 **`nature_editorial`** 风格：纯白画布 `#FFFFFF`、无上/右多余边框、无冗余粗网格、外向细刻度、色盲安全发散与分类色板、300 dpi + 可编辑矢量字体（`pdf.fonttype=42`、`svg.fonttype="none"`）、确定性无时间戳输出。
   - **视觉通道唯一化**：高饱和红/蓝严格保留给生物学方向（上调/下调），证据等级（L3/L2/L1）改用**无彩色黑白灰阶 + 实心点阶梯**编码，消除“蓝轴 vs 蓝下调”“金排名 vs 金星标 vs 金变号”等色相冲突。
2. **概念摘要与定量数据严格解耦（独立的 Graphical Abstract 入口）**
   - 手绘科学插画风（`editorial_ink`）与像素概念风（`pixel_abstract`）**绝不直接套用于高密度定量数据图（如 134 步反应全图或火山图）**，而是收敛到一个独立的 **Graphical Abstract（GA）子系统**（`omicsweave.visualization.graphical_abstract` 与 CLI 命令 `omicsweave abstract`）。
   - GA 子系统强制执行**小规模概念摘要护栏（`MAX_HUBS = 12`）**：仅用于 4–12 个核心代谢/通路枢纽的顶层故事呈现（论文图文摘要、综述 Working Model、PPT 总结页、期刊封面、项目 README 展示）。
   - **数据诚实性不变**：GA 可直接从图 E 的视图表（`map_view` + `region_summary`）自动聚合生成，所有枢纽方向、变化幅度与证据层级均严格溯源至真实统计评估表，同时支持通过声明式 YAML/dataclass 独立绘制任意其他课题的概念摘要图。
3. **高内聚、低耦合、跨项目可直接复用**
   - 所有通用视觉组件下沉到 `src/omicsweave/visualization/`，**不依赖 `omicsweave.contract`、`pydantic`、`seaborn` 或特定课题/集成业务模块**（仅依赖 `numpy`、`pandas`、`matplotlib`、`PyYAML`）。
   - 其他项目（如 `dio_project`、`dio_glp1ra_project` 或外部仓库）可直接 `from omicsweave.visualization import ...` 复用主题预设、色板、矢量徽章原语与 Graphical Abstract 渲染器，且导入过程零全局状态污染。

---

## 2. 模块目录架构与职责划分

```text
src/omicsweave/visualization/
├── __init__.py                    # 统一对外公共 API（跨项目复用入口）
├── palettes.py                    # 语义化色板与发散色阶（色盲安全，按预设组织）
├── theme.py                       # 风格预设引擎（默认 nature_editorial）、rcParams、期刊栏宽工具
├── figures.py                     # 确定性矢量/位图保存（PDF/SVG/PNG、Type 42 字体回退、去时间戳）
├── primitives.py                  # 通用矢量图元（三点阶梯徽章、角标排名、耦合轨道底带、去噪四轴微条、贝塞尔箭头、区室卡片）
├── volcano.py                     # 通用差异火山图（消费 nature_editorial 色板）
└── graphical_abstract/            # 独立 Graphical Abstract 子系统（严格 ≤ 12 枢纽）
    ├── __init__.py                # 公共入口：AbstractSpec, build_abstract_from_pathway_view, render_graphical_abstract
    ├── spec.py                    # 声明式数据结构、MAX_HUBS=12 护栏、YAML 读取、通路视图自动聚合器
    ├── styles.py                  # 三套独立渲染皮肤：nature_editorial / editorial_ink（手绘科学插画）/ pixel_abstract（现代像素风）
    └── renderer.py                # 统一渲染调度、确定性导出与溯源注脚
```

### 依赖方向约束（单向无环）

```text
  [下游业务绘图：integration/plotting/*, omics/*/plotting/*, 课题 plot_results.py]
                                      │
                                      ▼
       [omicsweave.visualization.{theme, palettes, figures, primitives}]
                                      ▲
                                      │
           [omicsweave.visualization.graphical_abstract.{spec, styles, renderer}]
```

- `visualization` 核心层**绝不反向导入** `integration`、`mfa`、`proteomics`、`transcriptomics`、`metabolomics` 或 `contract`（`volcano.py` 的 `write_bundle_volcanoes` 保持局部/向后兼容导入）。
- `integration/plotting/__init__.py` 采用延迟导出（`__getattr__`），避免导入 `pathway_map_plots` 或 `dossier_plots` 时被动导入旧版 `seaborn` 绘图模块并触发全局样式污染。

---

## 3. `nature_editorial` 视觉规范与多预设体系

### 3.1 风格预设清单（`theme.py`）

| 预设名 | 定位场景 | 背景 / 墨色 | 核心特征 |
|---|---|---|---|
| `nature_editorial`（**默认**） | 期刊正文主图、火山图、图 A–D、图 E | `#FFFFFF` / `#111827` | 纯白纸面、开放坐标轴（去上/右框线与默认粗网格）、外向刻度、300 dpi、色盲安全 `RdBu` |
| `metro_atlas` | 完整版图 E 诊断总图 / 补充材料 S-Fig | `#FFFFFF` / `#0F172A` | 强化分区结构线与轨道对比，突出多通路拓扑层次 |
| `slide_dark` | 学术会议 PPT / 大屏汇报 | `#111827` / `#F3F4F6` | 深石板灰底面，未变化骨架退入暗背景，珊瑚红 `#FB7185` / 冰川蓝 `#38BDF8` 高亮显著变化 |

### 3.2 统一语义色板（`palettes.py`）

- **方向色（唯一保留给生物学升/降）**：
  - `up = "#C83E37"`（编辑朱红），`down = "#2166AC"`（编辑深蓝），`neutral = "#9CA3AF"`
- **连续发散色阶（`diverging_cmap`，截断 `±2.0`）**：
  - `#1B4F8A → #3B82C4 → #A8C9EA → #F7F7F7 → #F4B6A8 → #D9483B → #8F1D24`
  - 零点为极浅灰白 `#F7F7F7`，与缺测底色 `#F3F4F6` + 纹理清晰区分；两端亮度控制在可清晰承托黑色 FDR 显著描边的区间。
- **证据等级无彩色阶（彻底与红/蓝涨跌解耦）**：
  - `L3 (Confirmed Flux)`：深炭实心 `#111827` + 3 实心点 `●●●`
  - `L2 (Flux-Supported)`：石墨实心 `#374151` + 2 实心点 `●●○`
  - `L1 (Enzyme-Predicted)`：纯白底 + 虚线板岩框 `#64748B` + 1 实心点 `●○○`
- **分类色板（Okabe-Ito / ColorBrewer Dark2 调和，色盲安全）**：
  - `("#2166AC", "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#A6761D", "#666666")`
  - 替换 `scoring_plots.py` 原有的红绿对比色（`#1B7837` vs `#B2182B`）。

### 3.3 期刊栏宽标准（`figures.py`）

- `SINGLE_COL_MM = 89.0`（≈ 3.50 in），`ONE_HALF_COL_MM = 120.0`（≈ 4.72 in），`DOUBLE_COL_MM = 183.0`（≈ 7.20 in），`MAX_PAGE_H_MM = 247.0`（≈ 9.72 in）。
- 要求任何面向印刷的 `publication` 模式图件按实际印刷栏宽设计画布，确保正文与标注最终字号 ≥ 5.5–6.5 pt。

---

## 4. 图 E 符号徽章升级与双模式输出

### 4.1 符号与徽章重设计（基于 `visualization/primitives.py`）

1. **证据等级 `L3 / L2 / L1` → 三点信号阶梯徽章（`draw_tier_badge`）**：
   - 用 `●●● L3` / `●●○ L2` / `●○○ L1` 的无彩色材质与实心点数编码证据强度，无需查表即可直觉读出等级递进，且与正文表格 `L1–L3` 命名保持一致。
2. **排名 `#n` → 卡片左上角孤立角标（`draw_rank_tab`）**：
   - 仅在 `atlas` 模式下显示在关键卡片左上角（如 `#1`），物理上与右侧等级/耦合徽章彻底分离，消除 `(1) [🔗 1]` 双数字歧义；`publication` 模式默认隐藏。
3. **化学计量锁死通量组 `g` → 拓扑轨道底带 + 链标（`draw_coupled_track` + `draw_coupled_badge`）**：
   - 同一 `flux_group` 中具有方向性 F 轴信号的连续边下方，绘制柔和的半透明靛青色轨道包络底带（`#E0E7FF`），从空间上直接把 `CS–ACONT–ICDH` 等耦合步骤串成一条“同检验轨道”，同时在卡片右侧标清晰的 `G1` 链环胶囊。
4. **四轴微条 `[F|C|R|M]` 视觉降噪（`draw_axis_strip`）**：
   - 表头改用中性灰阶（不再用蓝/琥珀色与下调/变号撞色）；
   - 8 种微小字符折叠为 3 级预注意编码：
     - **有方向信号**：红/蓝实心三角（changed/established）、红/蓝空心三角（suggestive）；
     - **测得无变 / 冲突**：干净横杠 `—`（unchanged）、红蓝对角斜拼方格（mixed，直观替代 `≠`）；
     - **缺测 / 无知识**：极淡灰底 + 柔和斜纹（不再铺满 `/`、`–` 字符噪点），仅在必要格保留轻量 `×`（不可辨识）、`!`（样本不足）与 `?`（未检出）。

### 4.2 图 E 双模式（`mode = "atlas" | "publication"`）

| 维度 | `atlas` 模式（默认，完整审计总图） | `publication` 模式（期刊双栏聚焦主图） |
|---|---|---|
| 目标用途 | 补充材料 S-Fig、全量数据核查、A3/屏幕缩放 | 论文正文双栏主图（宽度锁定 183 mm 级比例） |
| 覆盖区域 | 全部 13 个区域 + 右侧覆盖/融合/G16 仪表盘 | 自动聚焦有 MFA 通量或 L1–L3 关键信号的活跃区域（也可通过 `regions` 手动指定），自动裁剪边界框；非活跃区折叠为底部紧凑单行覆盖摘要条 |
| 非关键反应 | 保留全部四轴小条与排名角标 `#n` | 背景反应仅保留干净酶名，隐藏冗余非关键小条与探索性 `#n` 排名，突出主故事线 |
| 印刷可读性 | 放大字体与碰撞避让优化（修复线粒体标题与 `Mito β-Oxidation` 卡片重叠） | 真实目标印刷字号 ≥ 5.5–6.5 pt，适配双栏整页排版 |

配置接口（`experiment.yaml` 的 `integration.parameters.reaction_dossier.pathway_map` 与 `omicsweave map` CLI）：
- 新增可选键：`mode: atlas | publication`、`style: nature_editorial | metro_atlas | slide_dark`、`regions: []`。
- CLI 新增参数：`omicsweave map --mode publication --style nature_editorial --region glycolysis_gng --region pyruvate_tca ...`。

---

## 5. 独立 Graphical Abstract 子系统设计

### 5.1 为什么独立？
手绘插画与像素风只有在**高抽象、低元素密度（4–12 个核心枢纽）**时才具备强表达力与高级美感。将其封装为独立模块，既保护了定量图件的学术严肃性，又为论文图文摘要、综述、汇报 PPT、项目主页和期刊封面提供了零手工、数据驱动的精品生成管线。

### 5.2 数据契约（`visualization/graphical_abstract/spec.py`）
- **硬护栏**：`MAX_HUBS = 12`（超出即抛出 `ValueError`，防止误用于全量代谢网络）。
- **核心对象**：`AbstractHub`、`AbstractFlow`、`AbstractCompartment`、`AbstractSpec`。
- **两种输入路径**：
  1. **数据驱动自动生成（`build_abstract_from_pathway_view`）**：直接读取图 E 的 `pathway_map` 视图表与 `region_summary`，按核心功能枢纽（如 Blood Glucose、Glycolysis、PPP、Pyruvate Hub、TCA Cycle、Ketone Bodies、Fatty Acid / FAO、Glycerolipids、Amino Acids）自动汇总方向、变化幅度与最高证据等级，自动过滤空区，产出 ≤ 10 枢纽的 `AbstractSpec`。
  2. **跨项目通用声明式 YAML（`load_abstract_spec`）**：任意非代谢课题（免疫、神经、多组织互作等）只需写一份轻量 YAML，即可复用全套三套风格皮肤渲染。

### 5.3 三套专用渲染风格（`visualization/graphical_abstract/styles.py`）
1. **`nature_editorial`**：纯白极简顶刊机制图（清晰无衬线排版、细线圆角容器、方向箭头与三点证据阶梯）。
2. **`editorial_ink`（科学手绘插画风）**：
   - 暖宣纸底 `#FAF8F2` + 深炭墨线 `#1F2937`；
   - **有机容器 + 精密数据分层**：线粒体双层膜与内嵴、过氧化物酶体、血液波浪边界用算法生成的柔和多笔触有机墨线 + 低饱和水彩矿物底色晕染；枢纽文字、定量箭头与证据徽章保持绝对清晰的矢量渲染。
3. **`pixel_abstract`（现代像素概念风）**：
   - 专为 PPT 总结页、期刊封面构思、开源项目首页与科普展示预留；
   - 基于受控克制调色板在 16×16 像素网格上程序化绘制硬边缘无抗锯齿像素图标（葡萄糖六边形、丙酮酸节点、线粒体、酮体火焰、脂滴、氨基酸、ATP、血流血管等），配合像素阶梯边框枢纽卡、正交像素流向箭头与高可读标签。

4. **`metro_atlas`（Beck 代谢地铁换乘图风，同时可用于 GA 与正文图 E）**：
   - 基于独立解耦的 `visualization/metro.py` 原语库实现八向圆角折线轨道（0°/45°/90°）、内嵌方向 chevron、圆角跑道闭环（TCA / PPP 等循环通路）、双层圆环换乘站与低存在感浅灰圆角区室底板；
   - **证据等级编码**：实心加粗轨道 = Flux-Supported（L2/L3 实测通量支撑），空心白芯双轨（hollow double-rail）= Enzyme-Predicted / Omics-only（L1 仅酶预测），极细淡灰轨 = 未评估骨架；
   - **双配色语义开关（`color_mode`）**：
     - `color_mode="semantic"`（默认）：严格沿用正文统一的红/蓝/灰生物学涨跌语义色，与同论文其他 `nature_editorial` 统计图浑然一体；
     - `color_mode="pathway"`（可选）：启用克制的低饱和莫兰迪学术通路分区色（`METRO_PATHWAY_COLORS`），适用于无实验涨跌对比的综述或机制导览图。

### 5.4 正文分专题代谢子通路图（`pathway_submaps.py`）与多组学默认输出

为解决“161 反应全景图 E 适合附件诊断、正文需要单/双栏聚焦机制图”的出版诉求，新增独立子图预设目录 `integration/plotting/pathway_submaps.py`，与图 E 共享同一份 `pathway_map` 视图表、同一套 `network_layout.yaml` 坐标及默认 `metro_atlas`（`color_mode="semantic"`）风格：

- **5 张默认核心子图**（按信号感知自动出图，无评估数据时跳过并记日志）：
  1. `central_carbon`（糖原 + 糖酵解/糖异生 + PPP + 丙酮酸/TCA）
  2. `fatty_acid_ketone`（长链脂肪酸 + 酮体 + 过氧化物酶体 β-氧化，自动紧裁线粒体空白侧）
  3. `glycerolipid`（甘油与甘油三酯/磷脂合成，单栏紧凑幅面）
  4. `amino_acid`（氨基酸与 BCAA 代谢）
  5. `urea_creatine`（尿素循环、多胺与肌酸）
- **2 张可选扩展子图**（存在显著/活跃信号时自动附带，或按需指定）：`nucleotide`（嘌呤/嘧啶）、`redox_one_carbon`（一碳/谷胱甘肽/NAD 氧化还原）。
- **多组学默认交付布局**（`reaction_dossier` 运行与 `Delivery.copy_maps()` 自动联动）：
  - `figures/pathway_map/<组织>/map_<比较>.{pdf,png,tsv}`（图 E 全景总图）
  - `figures/pathway_map/<组织>/submaps/<比较>_<子图ID>.{pdf,png}`（正文分专题子图）
  - `figures/pathway_map/<组织>/abstract_<比较>.{pdf,png}`（高层地铁概念摘要图 GA）

### 5.5 图标化微条、自解释徽章与纯矢量状态符号
- **四层多组学微条（`strip_mode="omics"` 默认，支持 `--strip-mode axes` 切回 S3 四评估轴）**：
  - `map_view.py` 为每条反应同时填充真实分层读数 `P_state`/`P_value`（蛋白）与 `R_layer_state`/`R_value`（mRNA），消除早期反应卡片中 `R = Regulation` 与聚合卡片 `R = mRNA` 的语义冲突；
  - 微条顶栏以 4 枚纯矢量微图标取代难懂缩写：`⇉`（Flux）、`◑` 活性口袋酶球（Protein）、`∿` 波形单链（mRNA）、`⬡` 苯环/六边形（Metab）。
- **状态符号全面图标化**：所有格状态均由 `draw_status_cell` 纯矢量绘制并在总图/子图图例中给出实景色块对照：实心 `▲▼`（Sig. FDR<0.05）、空心 `△▽`（Suggestive p<0.05）、红蓝斜拼（Mixed）、横杠 `—`（Unchanged）、细斜线 `/`（Not measured）、圆圈斜杠 `⊘`（Unidentifiable flux）、虚线空圈 `◌`（Not detected）、琥珀三角叹号 `⚠`（QC / Low-n）。
- **自解释示踪剂证据徽章、串联通量徽章与净通量反转徽章**：
  - 全库统一内置 16 种规范同位素示踪剂注册表（`mfa/tracers.py`）：现有 4 种明确标为 `U-13C-glucose`、`U-13C-glycerol`、`U-13C-alanine`、`U-13C-lactate`（向后兼容历史裸名），并完整预置 12 种后续示踪剂（`U-13C-palmitic acid`，以及 `alanine` / `glutamine` / `glycine` / `leucine` 的 `U-13C`、`U-15N` 与 `U-13C-15N` 双标记协议）；
  - 三点阶梯徽章直接携带实际约束该反应的示踪剂短标：`●●● <多示踪剂名>`（如 `●●● U-13C-Glc+Ala` / `●●● U-13C-Gln+U-15N-Gln`，L3 多示踪剂一致确认）、`●●○ <示踪剂短标> Flux`（如 `●●○ U-13C-Glc Flux`，L2 MFA 通量支撑）、`●○○ Enz-pred`（L1 仅酶容量预测）；
  - `map_view.py` 增补 `flux_tracers` 与 `flux_pooling_mode`（`single_tracer` / `multi_tracer_pooled` / `multi_tracer_conflict` / `unconstrained`），在静态总图/子图图例与交互 HTML 图例/Tooltip/抽屉中完整注明多示踪剂分层逆方差合并（Pooled IVW）与单示踪剂退化规则；
  - `Linked flux` 徽章在静态 PDF/PNG 与交互 HTML 中均保留，并在图例与抽屉中明确界定其统计学性质：**共享 1 个独立自由度的化学计量冗余去重提示，而非更高置信度**（`v_fwd − v_rev` 在全部独立拟合中满足 `|r| ≥ 1−1e-6`，`grouped_bh` 校正中整组仅计 1 次检验；顶栏 KPI 同步报告独立通量链数）；
  - 组间净通量方向反转以纯矢量琥珀圆上下反向双箭头徽章 `draw_sign_reversal_badge`（`⇄`）标注，与无单一净方向的可逆反应轨道双向箭头（`< >`）在总图与子图图例中并列实景区分；
  - **酶容量推断空心柔彩双轨（`omics_edge_color="tinted"` 默认，支持 `--omics-edge-color grey` 切回保守灰）**：对于缺乏方向性 MFA 通量但达到 `Enz-pred (L1)` 酶容量证据的反应，采用空心白芯 + 饱和度封顶（≤72%）的柔和红/蓝双轨表示推断方向，绝不因代谢物池浓度给反应边上色，实现与 Graphical Abstract 的视觉语法统一；
  - Graphical Abstract 底栏新增量化图例与聚合细则卡片（`AbstractSpec.rules_note`），完整展示三点证据徽章、`Linked flux` 轨道以及自动聚合量化标准（枢纽等级取最高 claim；方向优先级 `Flux (F sig./sug. 多数) > Enzyme (C 多数) > Pool (|median log2FC| ≥ 0.45)`；`hits` 取 `L2+` 计数否则 `L1`；连线继承证据等级更高一端）。
- **全库绘图统一收敛**：
  - `primitives.py` 的 `draw_status_cell` / `draw_layer_icon` / `draw_tier_badge` 同时支持 y 向下画布（图 E / 子图 / GA）与标准 y 向上非等比坐标系（图 A–D 矩阵与面板）；
  - 反应档案图 A–D（`dossier_plots.py`）、Graphical Abstract 全部 4 套皮肤（含 `pixel_abstract` 像素能量条）以及单组学 / MFA 绘图模块（`mfa`、`metabolomics`、`proteomics`、`transcriptomics`、`integration`）统一采用 `NATURE_EDITORIAL` 红蓝语义色板与 **300 dpi** 出版输出规范。

### 5.6 CLI 入口
```bash
# 正文图 E（默认输出 pdf + png + 单页交互 html，统一默认 metro_atlas 学术风格）与分专题子图
omicsweave map --experiment all_tissues --tissue liver
omicsweave map --experiment all_tissues --tissue liver --submap central_carbon --submap fatty_acid_ketone
omicsweave map --experiment all_tissues --tissue liver --style economist --format pdf --format png --format html
omicsweave map --experiment all_tissues --tissue liver --omics-edge-color grey --no-abstract

# Graphical Abstract（默认 metro_atlas；可选 economist / editorial_ink / pixel_abstract / nature_editorial）
omicsweave abstract --experiment all_tissues --tissue liver --contrast HFD_vs_LFD
omicsweave abstract --experiment all_tissues --tissue liver --contrast HFD_vs_LFD --style economist
omicsweave abstract --experiment all_tissues --tissue liver --contrast HFD_vs_LFD --style metro_atlas --color-mode pathway

# 从任意自定义 YAML 概念草图生成（跨项目复用，不依赖整合包）
omicsweave abstract --spec my_project_ga.yaml --style metro_atlas --output ./ga_out
```

---

## 7. 方案 A 单页交互式 HTML 与《经济学人》可选主题体系（M4）

1. **统一学术默认 `metro_atlas` + 一级可选主题 `economist`（全库保留自由切换）**：
   - 全库绘图（通路总图、分专题子图、Graphical Abstract、反应档案图 A–D、火山图、全部单组学绘图 `_config.py` 与交互 HTML）统一以 `metro_atlas` 为默认学术出版风格（`DEFAULT_PRESET = "metro_atlas"`）；
   - 同时在 `palettes.py` / `theme.py` / `primitives.py` / `assets/glyphs/economist/*.svg`（全套 24 枚）完整内置《经济学人》（`economist`）编辑部视觉皮肤（淡黄暖灰纸底 `#EBEBDF`、左上角编辑部红标 `#E3120B`、实验组红 `#E3120B` / 对照组深蓝 `#006BA2` / London 暖灰渐变、细水平网格与炭黑底基线），命令行通过 `--style economist`、网页端通过顶部主题下拉框即可一键切换。
2. **方案 A 单页离线可交互 HTML（`pathway_map_html.py`，默认 `formats: ["pdf", "png", "html"]`）**：
   - 共享 `build_pathway_map_scene()` 的精确正交路由、空心柔色双轨与避让卡片坐标，并按 `em_units(pt, unit_in)` 1:1 精确换算 `metro_atlas` PDF 的精细边框线宽、圆角区室、圆点标题、聚卡片 `P/R/M` 色块与徽章字号；零外部 CDN 依赖、确定性逐字节可重现；
   - 内置浏览器端 `getBBox()` 零溢出自校准守护（`fitSvgLabels`），自动适配 Windows/macOS/Linux 系统字体宽度差异，确保卡片、证据徽章、聚合卡标题与边界胶囊均无文本越界；
   - 浏览器内支持无缝实时切换主题（`Metro Atlas` / `Economist` / `Nature Editorial` / `Slide Dark`）、分专题一键缩放聚焦、按反应/基因/蛋白/代谢物/EC 瞬时搜索、证据层级过滤、`Enz-pred` 柔色/灰轨切换、悬停 `Linked flux` 整链脉冲高亮与当前视图 SVG/PNG 导出；
   - 右侧 Reaction Dossier 抽屉内置 **逐只动物/样本散点小图**（读取 `dossier_samples.parquet` 与 `dossier_sample_contrasts.parquet`，含实验组红点、对照组蓝点、其他组灰点、中位线、MFA 拟合 IQR 误差棒、盒约束空心菱形与 0 基线）及可排序全量成员明细表；交付模块 `delivery.py` 同步复制 `map_<比较>.html` 至交付目录。

---

## 6. 分步实施清单与验收标准

1. **核心通用库（`visualization/`）**：实现 `palettes.py`、`theme.py`、`figures.py`、`primitives.py`，更新 `volcano.py` 与 `__init__.py`；保持 `test_theme.py` 与 `test_volcano.py` 全部向后兼容。
2. **图 E 升级（`integration/plotting/pathway_map_plots.py`、`integration/pathway_map.py`）**：接入统一原语、修复已知标签碰撞、支持 `atlas` 与 `publication` 双模式与多预设。
3. **Graphical Abstract 子系统（`visualization/graphical_abstract/` + `cli.py`）**：实现 `spec.py`、`styles.py`、`renderer.py` 与 `omicsweave abstract` 命令。
4. **全库绘图收敛**：更新 `dossier_plots.py`、`scoring_plots.py`、`integration/plotting/__init__.py`（延迟导入解除 `seaborn` 强耦合）以及单组学 `_config.py` / `pathway_plots.py` 的色板与 300 dpi 默认值。
5. **测试与实图目视核验**：新增 `tests/visualization/` 单元测试覆盖预设、原语、双模式图 E 与三套 GA 皮肤，并用 DIO 真实肝脏数据渲染全部样图逐张目视检查。
