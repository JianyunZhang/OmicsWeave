# 多示踪剂 MFA 通量权重分配规则与稳健性分层实证评估报告

> **文档类型**：统计方案对比与实证评估报告（v1.0，2026-11-22）  
> **实证数据集**：`local/projects/dio_project` 单组学结果包 `run_20261121_174318_7322` 与整合结果包 `run_20261121_182112_fa44`（4 种成对示踪剂：glucose / glycerol / lactate / alanine + 1 种未成对示踪剂：alanine_13c15n；7 个组织 × 3 个比较；共 **902** 条可辨识反应行，其中 $k \ge 2$ 多示踪剂行 **535** 条，$k = 1$ 单示踪剂行 **367** 条）  
> **关联文档**：[01_mfa_flux_assessment_thresholds_and_improvements.md](01_mfa_flux_assessment_thresholds_and_improvements.md) · [统计方法目录总览](README.md)

---

## 1. 评估目的与核心结论

在多示踪剂 MFA 整合流程中，示踪剂权重分配直接影响四个关键环节：
1. **点估计与置换检验**（`multi_tracer._pool_reaction()` 的 `diff`、`log2fc`、`permutation_pvalue`）；
2. **标准化效应量与异质性检验**（`multi_tracer._pool_standardized()` 的 `smd` / `effect_size`、`smd_heterogeneity_pvalue`）；
3. **E5 单只动物稳健性**（`stats/robustness.py` 的 `single_subject_driven`）；
4. **F 轴方向稳定性判定**（`dossier/assessment.py` 的 `flux_axis()`）。

通过对 `dio_project` 全部 902 条可辨识 MFA 反应（535 条多示踪剂反应）进行六种候选权重方案的定量基准测试，核心结论如下：

1. **现行绝对尺度 IVW（W1：`tracer_weights`）存在显著的“低数值尺度绑架”现象**：
   - 在 **27.3%**（146/535）的多示踪剂反应中，单一示踪剂占据了 **>90%** 的绝对权重；在 **7.7%**（41/535）的反应中甚至超过 **99%**（例如心脏 `GPS_c` 甘油占 `99.3%`、心脏 `PGM_c/ENO_c` 单一示踪剂占 `98.6%`、肝脏 `PYR_trans` 4 示踪剂中单一示踪剂占 `95.8%`、心脏 `RPE_c` 葡萄糖占 `100%`）。
   - 更关键的是，在 **14.6%（78/535）** 的多示踪剂反应中，绝对 IVW 得到的 `diff` 与标准化合并得到的 `smd` **符号相反**（`log2fc` 与 `smd` 异号比例达 **16.1%**，86/535）。而当前 `flux_axis()` 的方向取自 `log2fc`（W1），效应量大小却取自 `|smd|`（W2），两套权重体系存在口径分裂。
2. **现行标准化 Hedges' $g$ 权重（W2：`smd_weights`）彻底消除了尺度绑架，但忽略了拟合方差（`fit_var`）**：
   - W2 使最大单示踪剂权重中位数从 `0.728` 降至 `0.394`，有效示踪剂利用率 $\text{median}(N_{\text{eff}}/k)$ 从 `0.628` 跃升至 **`0.977`**；
   - 然而，在 346 个可辨识 `(组织, 反应, 示踪剂)` 层中，动物内拟合方差占总方差的比例 $\bar{f}/(s^2+\bar{f})$ 中位数高达 **47.0%**（84.4% 的层超过 20%，46.5% 的层超过 50%）。W2 忽略 `fit_var` 会导致动物均值偶发相近、但拟合误差很大的反应出现极端虚高 $g$ 值（如心脏 `GLN_exchange` 动物均值 `[154.21, 154.50]` 而 `fit_std = 59.0`，W2 算出 `g = +49.13`）。
3. **权重方案选型建议（供 Phase 2 决策）**：
   - **短期（Phase 1 已落地）**：保持 `multi_tracer.py` 现有主估计向后兼容，新增输出 **`smd_weights`** 列与 `tracer_weights` 并列透明披露；在 `stats/robustness.py` 中修复多示踪剂未分层缺陷（改为仅在可辨识示踪剂内执行固定尺度标准化 SMD 分层 LOO，在 `dio_project` 上清除了 **21 个假阳性**并新捕获了 **50 个被未分层噪音掩盖的真单鼠驱动案例**）。
   - **中期推荐（Phase 2 主选方案）**：
     - **合并效应量与方向统一**：采用 **W3b（可靠度加权 SMD：`w_t = 1 / (Var(g_t) · (1 + f_bar / s_p^2))`）** 或 **W4（可辨识度先验加权 SMD：`well_constrained=1.0, moderate=0.5`）** 作为多示踪剂标准化合并权重，并让 `flux_axis()` 的方向、效应量与分层 LOO / Bootstrap `direction_stability` 统一在同一标准化权重体系上（彻底消除 `sign(log2fc) != sign(smd)` 的分裂及 `min(computed) >= 1.0` 的最弱示踪剂一票否决）；
     - **置换检验尺度升级**：将分层置换统计量由绝对尺度 $\Delta_t / V_t$ 评估升级为 **W5（层内尺度归一化相对差 $\Delta_t / \bar{S}_t$）** 或标准化 $t/g$ 统计量，避免分层置换 $p$ 值被绝对通量最小的那一种示踪剂独占。

---

## 2. 六种候选多示踪剂权重分配规则对比

### 2.1 数学定义

记第 $t$ 种可辨识示踪剂（$n_{1,t}, n_{2,t} \ge 2$，$df_t = n_{1,t} + n_{2,t} - 2$）的组均值为 $\bar{a}_t, \bar{b}_t$，组间差 $\Delta_t = \bar{a}_t - \bar{b}_t$，动物间样本方差为 $s_{1,t}^2, s_{2,t}^2$，动物内拟合方差均值为 $\bar{f}_{1,t}, \bar{f}_{2,t}$，合并动物方差为 $s_{p,t}^2 = \frac{(n_{1,t}-1)s_{1,t}^2 + (n_{2,t}-1)s_{2,t}^2}{df_t}$，Hedges 小样本校正因子 $J_t = 1 - \frac{3}{4df_t - 1}$。

| 方案 | 名称 | 合并量 | 权重公式 $w_t$ | 是否尺度无关 | 是否纳入 `fit_var` | 是否区分 `well` vs `moderate` |
|---|---|---|---|:-:|:-:|:-:|
| **W1** | 现行绝对 IVW（`tracer_weights`） | 绝对差 $\Delta_t$ | $1 / \left(\frac{s_{1,t}^2 + \bar{f}_{1,t}}{n_{1,t}} + \frac{s_{2,t}^2 + \bar{f}_{2,t}}{n_{2,t}}\right)$ | **否**（强偏向小尺度示踪剂） | 是 | 间接（通过 $\bar{f}_t$） |
| **W2** | 现行标准化 SMD（`smd_weights`） | 动物方差 Hedges' $g_t = J_t \Delta_t / s_{p,t}$ | $1 / \left(\frac{n_{1,t}+n_{2,t}}{n_{1,t}n_{2,t}} + \frac{g_t^2}{2(n_{1,t}+n_{2,t})}\right)$ | **是** | **否** | 否 |
| **W3** | 全方差标准化 SMD（Total-Var $g_{\text{tot}}$） | 全方差 Hedges' $g_{\text{tot},t} = J_t \Delta_t / \sqrt{s_{p,t}^2 + \bar{f}_{p,t}}$ | $1 / \mathrm{Var}(g_{\text{tot},t})$ | **是** | **是**（进分母） | 间接 |
| **W3b** | 可靠度加权 SMD（Reliability-Weighted） | 动物方差 Hedges' $g_t$ | $\frac{1}{\mathrm{Var}(g_t) \cdot (1 + \bar{f}_{p,t} / s_{p,t}^2)}$ | **是** | **是**（降权不压缩 $g_t$） | 间接（拟合越稳权越高） |
| **W4** | 可辨识度先验加权 SMD | 动物方差 Hedges' $g_t$ | $q_t / \mathrm{Var}(g_t)$（$q_{\text{well}}=1.0,\ q_{\text{mod}}=0.5$） | **是** | 间接（通过可辨识分层） | **是** |
| **W5** | 层内尺度归一化相对差 IVW | 相对差 $\delta_t = \Delta_t / \bar{S}_t$（$\bar{S}_t = \frac{|\bar{a}_t|+|\bar{b}_t|}{2}$） | $\bar{S}_t^2 / V_{\text{abs},t}$ | **是** | **是** | 间接 |

### 2.2 全库实证定量对比（`dio_project` · $k \ge 2$，共 535 行）

定义有效示踪剂数（Kish $N_{\text{eff}} = 1 / \sum_t \tilde{w}_t^2$，归一化后 $N_{\text{eff}} / k \in [1/k, 1]$，越接近 1 表示信息越均衡、不被单一示踪剂垄断）：

| 方案 | 最大单示踪剂权重中位数 `med(maxW)` | 单示踪剂权重 `>0.80` 占比 | 单示踪剂权重 `>0.90` 占比 | 有效示踪剂利用率 `med(Neff/k)` | 与 W2（SMD）符号一致率 | 全库 `\|g\|≥1` 数（$k=1$ / $k\ge 2$） | 极端虚高 `\|g\|≥5` 数（`max\|g\|`） |
|---|--:|--:|--:|--:|--:|--:|--:|
| **W1** 绝对 IVW | **0.728** | **38.3%** | **27.3%** | **0.628** | **85.4%** | — | — |
| **W2** 标准化 SMD | 0.394 | 0.4% | 0.2% | 0.977 | 100.0% | 64 (56 / 8) | 2 (49.13) |
| **W3** 全方差 SMD | 0.386 | 0.0% | 0.0% | 0.981 | 94.0% | **17 (14 / 3)** | **0 (2.62)** |
| **W3b** 可靠度加权 SMD | 0.525 | 11.0% | 5.0% | 0.895 | 89.3% | 65 (56 / 9) | 2 (49.13) |
| **W4** 可辨识度先验 SMD | 0.439 | 0.4% | 0.2% | 0.971 | 95.1% | 64 (56 / 8) | 2 (49.13) |
| **W5** 尺度归一化相对差 IVW | 0.663 | 24.1% | 12.5% | 0.737 | 89.9% | — | — |

### 2.3 关键面板反应在不同权重规则下的表现（`HFD_vs_LFD`）

| 组织 / 反应 | $k$ | W1 `maxW` | W2 `smd` (`maxW`) | W3 `g_tot` | W3b `g_rel` (`maxW`) | W5 `rel` (`maxW`) | 现象解读 |
|---|--:|--:|--:|--:|--:|--:|---|
| 心脏 `PGM_c` / `ENO_c` | 2 | **0.986** | **+1.31** (0.51) | **+1.11** | **+1.29** (0.54) | +2.00 (0.66) | 所有标准化/相对方案高度稳健（均 >1.1），W1 却把 98.6% 权重压给了单示踪剂 |
| 心脏 `PYK_c` | 4 | 0.750 | +0.66 (0.30) | +0.54 | +0.65 (0.37) | +0.84 (0.52) | 4 示踪剂全正但组内离散度大，所有方案下标准化效应均 < 1.0 |
| 心脏 `G3PP_c` | 3 | 0.671 | **-1.10** (0.47) | -0.37 | -0.33 (0.82) | -0.41 (0.48) | `glucose` 动物方差小（$g=-3.16$）但拟合方差较高，纳入拟合方差后大幅衰减 |
| 心脏 `GPS_c` | 2 | **0.993** | +0.38 (0.51) | +0.39 | +0.72 (0.67) | +1.37 (0.79) | 经典尺度绑架案例：W1 中甘油独占 99.3% 令 `log2fc = -1.15`，而 W2/W3/W3b/W5 均为正 |
| 心脏 `TPI_c` | 2 | **0.992** | +0.32 (0.50) | +0.36 | +0.70 (0.68) | +1.66 (0.92) | 同 `GPS_c`，W1 反向（`diff < 0`），标准化方案均为正 |
| 心脏 `GLN_exchange` | 1 | 1.000 | **+49.13** (1.00) | **+2.62** | +49.13 (1.00) | +1.85 (1.00) | `lactate` 2v2 动物均值巧合贴近令 W2 虚高至 49.13；W3 纳入 `fit_std` 后回到合理值 +2.62 |
| 肝脏 `PYR_trans` | 4 | **0.958** | -0.78 (0.34) | -0.35 | -0.55 (0.49) | -0.77 (0.51) | 4 示踪剂全负，W1 单示踪剂占 95.8%，W2 权重最均衡（0.34） |
| 肝脏 `PC_m` | 2 | 0.768 | -0.81 (0.53) | -0.66 | **-0.98** (0.71) | -1.70 (0.80) | `alanine` 拟合更紧致，W3b 提升其权重至 71% 后效应接近 -1.0 |
| 皮下脂肪 `BHB_input` | 2 | 0.548 | -0.81 (0.58) | -0.51 | -0.76 (0.54) | -1.14 (0.58) | 两种示踪剂绝对尺度接近，W1/W2/W3b 权重基本一致（~0.55） |

---

## 3. 稳健性（E5）与方向稳定性权重规则评估

### 3.1 Phase 1 已落地的修复：`stats/robustness.py` 多示踪剂分层 LOO 实测效果

在 `dio_project/run_20261121_174318_7322` 的全部 162 条入选 MFA 反应上对比修复前后：

| 指标 | 修复前（未分层全样本混算） | 修复后（仅可辨识示踪剂 + 固定尺度 SMD 分层 LOO） |
|---|--:|--:|
| 标为 `single_subject_driven` 的反应数 | 89 / 162 (54.9%) | 118 / 162 (72.8%) |
| `driver_sample` 来自不可辨识或未成对示踪剂的数量 | **48 / 89 (53.9%)**（含 `poorly_constrained` 30、`bound_pinned` 18 及未成对的 `alanine_13c15n`） | **0 / 118 (0.0%)** |
| 清除的假阳性（有效示踪剂实际稳定，原被噪音示踪剂误杀） | — | **21 条**（如心脏 `GLN_exchange` 有效 `min_ratio = 0.988`、脑 `GLNS_c` `0.915`、肝脏 `G6PDH2R_PGL_GND_c` `0.704`、肝脏 `ACSS2_m` `0.668`、肾脏 `ACE_exchange` `0.645`） |
| 新捕获的真单鼠驱动（原被噪音示踪剂稀释掩盖） | — | **50 条**（如腓肠肌 `PEPCK_c` 有效甘油 HFD=`[273.66, 1.0, 1.0]` vs LFD=`[1.0, 23.26, 1.0]`，去 1 只直接反号 `ratio = -0.089`） |

### 3.2 F 轴方向稳定性聚合规则对比（为 Phase 2 提供决策依据）

以心脏 `G3PP_c`（`HFD_vs_LFD`：`glucose` $g=-3.16, \text{LOO}=1.0$；`lactate` $g=-0.98, \text{LOO}=1.0$；`alanine` $g=-0.12, \text{LOO}=0.67$；合并 `smd = -1.10`, `p = 0.0224`）为例：

| 稳定性聚合规则 | 逻辑 | `heart/G3PP_c` 判定 | 统计特性评估 |
|---|---|---|---|
| **现行规则**：`min(per_tracer_loo) >= 1.0` | 所有可辨识示踪剂各自 LOO 均不得反号 | 不稳定 → `undetermined` | 对弱贡献示踪剂过于严苛，出现“多一种弱示踪剂反而降低稳健性”的悖论 |
| **候选 A（推荐）**：分层合并 SMD LOO（`strat_loo_same == 1.0`）+ 分层 Bootstrap `direction_stability >= 0.8` | 每次在任一有效示踪剂中剔除 1 只动物，按 SMD 权重重算合并效应，要求合并方向零反号 | **稳定**（`strat_loo_same = 1.0`, `min_ratio = 0.808`）→ `suggestive` | 与合并估计量权重完全一致；单示踪剂时自动退化为精确组均值差 LOO |
| **候选 B**：主权重示踪剂过滤（仅对 `smd_weight >= 0.2` 且 `\|g_t\| >= 0.5` 的示踪剂要求 `loo == 1.0`） | 过滤近零效应的弱示踪剂 | 稳定 → `suggestive` | 规则简单但引入额外人为截断阈值，不如候选 A 连续自然 |

---

## 4. 综合推荐结论（供 Phase 2 实施拍板）

1. **避免 W3 直接替换主 `effect_size`**：W3 将拟合方差直接并入 $g$ 的分母在小样本（`n=2–3`）且多为 `moderate` 约束的在体 MFA 上过于保守（全库 `|g| ≥ 1` 从 64 骤降至 17），更适合作为辅助质控指标（用于识别如 `g = 49.13` 这类靠偶然低动物方差放大的极端值）。
2. **推荐在 Phase 2 采用“W2/W3b/W4 统一标准化权重 + W5 尺度归一化置换 + 候选 A 分层合并稳定性”组合**：
   - 在 `multi_tracer.py` 中：
     - **选项 2a（最稳健、变动最小）**：保留 W2（`smd` / `smd_weights`）作为默认标准化效应量，增加可配置的 `weighting_mode: smd | reliability_smd | identifiability_smd`（默认 `smd` 或 `identifiability_smd`），同时输出基于同一套权重计算的分层合并 `direction_stability`（Stratified LOO + Bootstrap）与 `g_total_var` 质控列；
     - **选项 2b（彻底解决绝对尺度主导置换检验与方向分裂）**：将分层置换检验统计量与 `flux_axis()` 主方向切到与标准化效应一致的尺度（W5 相对差或 W2/W3b 标准化 $g$），避免出现 14.6%–16.1% 的 `sign(diff/log2fc) != sign(smd)` 冲突。
3. **Phase 1 & Phase 2 落地总结与实测验证**：
   - **Phase 1**：`multi_tracer.py` 输出 `smd_weights`；`stats/robustness.py` 完成多示踪剂分层与可辨识性过滤修复；图 A / 图 E / HTML 彻底区分 `Unchanged (equiv.)`（`—`）与 `Undetermined`（`--`），并在导出表中补全 `effect`、`rule`、`flags` 等字段。
   - **Phase 2（已落地）**：
     - **权重与置换模式可选化（默认完全向后兼容）**：在 `experiment.yaml` 的 `omics.mfa.analysis.multi_tracer` 下新增：
       - `smd_weighting`: `"smd"`（默认，W2）\| `"reliability_smd"`（W3b，按动物方差/总方差可靠度加权）\| `"identifiability_smd"`（W4，按 `well_constrained=1.0 / moderate=0.5 / undetermined=0.25` 质量先验加权）；
       - `permutation_scale`: `"absolute"`（默认）\| `"relative"`（W5，各示踪剂层内按置换不变的 `mean(|values|)` 归一化后再做分层置换，彻底消除跨示踪剂绝对通量量级差对置换统计量的主导）；
       - `n_bootstrap`: `2000`（分层组内自助法重抽样次数），同时落盘 `direction_stability`、`pooled_loo_stability`、`g_total_var`、`smd_weights`、`smd_weighting`、`permutation_scale`。
     - **分层合并稳定性（S2）**：`assessment.py` 新增 `stratified_pooled_loo()`（固定层内 Hedges' $g$ 尺度与权重，避免 `n=2v2` 退化），`flux_axis()` 优先采纳分层合并的 `direction_stability` + `pooled_loo_stability`（或由 `dossier_samples` 即时计算 `"__pooled__"`），取代原有的最弱示踪剂一票否决 `min(computed) >= 1.0`，并对跨尺度符号相反的多示踪剂反应标记 `scale_direction_discordant`。
     - **`dio_project` 全库实测结果（4,200 条 `(tissue, contrast, flux_id)` 反应行）**：
       - F 轴状态变化仅 **3 条**（零假阳性膨胀），全部为原本被弱次要示踪剂 LOO 一票否决的真多示踪剂信号：`heart/HFD_vs_LFD/G3PP_c`（`|g|=1.10`，`undetermined → suggestive`）、`brain/ND_vs_LFD/ME2_c`（`|g|=1.39`，`undetermined → suggestive`）、`sq_adipose/ND_vs_LFD/FA_usage`（`|g|=1.03`，`undetermined → suggestive`）；
       - 全库透明标记 **72 条** `scale_direction_discordant` 跨尺度分歧反应（含 `heart/TPI_c`、`heart/GPS_c`、`heart/GLYCEROL_exchange`）；
       - 图 A（V4）在 `ASSESSMENT -> F` 的 `undetermined` 虚线下方显示精简门控标签（如 `PYK_c` 显示 `g=0.66`，`PDH_m` 显示 `1-subj`），交互式 HTML 抽屉（V5）完整呈现 **Flux Assessment Gate Checklist**、权重对照与 MFA 成员 `g` 值。
