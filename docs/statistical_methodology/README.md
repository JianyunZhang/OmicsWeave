# OmicsWeave 统计方案选用与阈值设定文档库

> **目录定位**：集中沉淀 OmicsWeave 在单组学差异分析、多示踪剂 MFA 整合、分轴反应评估（`reaction_dossier`）与代谢通路可视化中的**统计方法选用依据、默认阈值基线、边界问题诊断及演进方案**。  
> **建立日期**：2026-11-22  
> **配套架构文档**：[反应档案与分轴综合评估](../arch_design/omicsweave_reaction_dossier_and_assessment.md) · [分组织代谢通路总图（图 E）](../arch_design/omicsweave_pathway_map.md) · [通用能力与评估框架](../arch_design/omicsweave_generic_capabilities_and_evaluation.md)

---

## 1. 文档导航

| 编号 | 文档 | 核心内容 | 状态 |
|:--:|---|---|---|
| 01 | [01_mfa_flux_assessment_thresholds_and_improvements.md](01_mfa_flux_assessment_thresholds_and_improvements.md) | 多组学整合与通路绘制中 **ASSESSMENT Flux（F 轴）** 的端到端分析流程、`Sig.` / `Suggestive` / `unchanged` / `undetermined` 判定条件与阈值基线、“多示踪剂与 `pooled` 全同号仍显示 `Unchanged`” 的 5 大根因诊断，以及**展示绘图改进（V1–V5）**、**统计方法/阈值改进（S1–S6）**与 **Phase 4 通用质检及评估规则硬化**路线。 | v1.4（Phase 1–4 全部已落地：V1–V5, S1–S5, 通用 QC/贴边/规则硬化） |
| 02 | [02_mfa_multi_tracer_weighting_evaluation.md](02_mfa_multi_tracer_weighting_evaluation.md) | **多示踪剂 MFA 权重分配规则与稳健性分层实证评估**：基于 `dio_project` 全库 902 条可辨识反应（535 条多示踪剂行）定量对比 6 种权重方案（W1–W5）、E5 稳健性修复对照，以及 Phase 2 权重模式可选化（`smd_weighting` / `permutation_scale`）与分层合并稳定性落地实测。 | v1.1（Phase 2 已落地并完成实测验证） |
| 03 | [03_mfa_phase3_calibration_and_equivalence.md](03_mfa_phase3_calibration_and_equivalence.md) | **Phase 3 统计方法深化校准报告**：基于 48 万次蒙特卡洛空/备择模拟与 14,400 次 `dio_project` 真实标签置换，验证并固化 **S3**（`fdr_scope` 按组织校正与 `padj_global` 保留）、**S4**（多示踪剂分层置换辅助 `Suggestive` Path B）与 **S5**（五重护栏 F 轴通量 TOST 等效检验与 `L1` 防冲突联动）。 | v1.0（Phase 3 已落地并完成空假设校准） |

---

## 2. 核心阈值配置位置速查

OmicsWeave 的统计与评估阈值按职责分层配置，修改基线或开展敏感性对照时应优先在 `experiment.yaml` 中声明，避免硬编码修改默认值：

| 作用域 | 配置路径（`experiment.yaml`） | 代码默认值位置 | 控制内容 |
|---|---|---|---|
| 转录组通用样本质检（G1） | `omics.transcriptomics.analysis.sample_qc` | `src/omicsweave/transcriptomics/sample_qc.py` (`DEFAULT_SAMPLE_QC_CONFIG`) | `enabled` (True)、`filter_excluded` (True)、跨组织污染（`contam_*`）、本组织标记保留度（`marker_retention_*`）、转录组相关离群（`corr_*` / `var_corr_*`）、线粒体占比（`mt_*`） |
| 单示踪剂 MFA 可辨识性与差异分析 | `omics.mfa.analysis.identifiability` / `differential_flux` | `src/omicsweave/mfa/df_analysis.py` (`selection_thresholds`) / `mfa_experiment_config.py` | `bound_tolerance` (0.05)、`bound_rel_tolerance` (0.02)、`bound_pinned_min_fraction` (0.5，全局 + 单组双核查)、`pvalue_threshold` (0.05)、`fdr_threshold` (0.10)、`fdr_scope` (`per_tissue` \| `per_contrast`)、`log2fc_threshold` (0.5)、`hedges_g_threshold` (0.5) |
| MFA 通用样本与尺度漂移质检（G1） | `omics.mfa.analysis.sample_qc` | `src/omicsweave/mfa/sample_qc.py` (`DEFAULT_MFA_QC_CONFIG`) | `enabled` (True)、样本级 `flux_scale_factor` / `upper_bound_saturation` / `fitting_loss`、比较级 `flux_scale_drift`（`contrast_scale_total_log2` 0.30、`contrast_scale_median_log2` 0.15） |
| 多示踪剂 MFA 合并 | `omics.mfa.analysis.multi_tracer` | `src/omicsweave/mfa/multi_tracer.py` (`PoolingSettings`) | `exclude_unconstrained` (True)、`bound_pinned_min_fraction` (0.5，按比较单组贴边门控)、`include_fit_variance` (True)、`heterogeneity_alpha` (0.10)、`min_tracers` (1)、`fdr_scope` (`per_tissue` \| `per_contrast`)、`smd_weighting` (`smd` \| `reliability_smd` \| `identifiability_smd`)、`permutation_scale` (`absolute` \| `relative`)、`n_bootstrap` (2000)、置换次数 |
| 化学计量耦合通量组 | `omics.mfa.analysis.flux_groups` | `src/omicsweave/mfa/flux_groups.py` | `enabled` (True)、`min_abs_correlation` ($1 - 10^{-6}$) |
| 反应档案分轴综合评估 | `integration.parameters.reaction_dossier.rules` | `src/omicsweave/integration/dossier/assessment.py` (`RULES`) | `alpha` (0.05)、`equivalence_bound` (0.58)、`flux_min_abs_g` (1.0)、`flux_concordant_min_abs_g` (0.8)、`flux_equiv_g_bound` (1.0)、`heterogeneity_alpha` (0.10)、`direction_stability_min` (0.8)、`established_min_group` (3)、`opposition_min_state` (`suggestive`)、`rna_single_source_changed` (False)、`flux_scale_drift_gate` (True) |
| 同工酶权重、活跃占比与成员关系 | `integration.parameters.reaction_dossier.isozyme` / `reaction_members.yaml` | `src/omicsweave/integration/modes/reaction_dossier.py` / `mappings/reaction_members.py` | `minor_share` (0.2)、条件活跃占比 `active_abundance_share`、成员关系 `isozyme` / `complex_subunit` / `shared_subunit` / `pathway_step` |
| 面板声明与预期示踪剂 | `panels.yaml`（由 `integration.parameters.reaction_dossier.panels` 指定） | `src/omicsweave/integration/dossier/panels.py` | `segments`、`branch_nodes`、`expected_tracers`、`signed_modules`、`controls`、`selection.mode` |
| 后处理稳健性（E5） | 自动随 `run` / `stats` 执行 | `src/omicsweave/stats/robustness.py` | `DEFAULT_MIN_RATIO = 0.5`（留一法衰减比例门槛与 `single_subject_driven` 标记） |

---

## 3. 统计方案选用与阈值调整的治理原则

在后续迭代本目录记录的任何统计方法或阈值设定时，统一遵循以下四条原则：

1. **区分“工程缺陷修复”与“统计门槛放宽”**：
   - 修复计算逻辑偏差（如多示踪剂下 `robustness.py` 未按示踪剂分层导致不可辨识样本污染 LOO）属于工程缺陷修复，可直接推进；
   - 调整显著性/效应量/稳定性阈值属于统计口径变更，必须有明确的方法学理由，并在敏感性对比中保留基线结果对照。
2. **展示清晰度优先于强行拔高结论等级**：
   - 面对“全同号但未达显著”的反应，首要解决的是**消除图例与视觉误导**（将 `Undetermined` 与 `Unchanged` 明确分离，并透明暴露效应量与未达标原因），而不是通过无约束降低门槛人为制造 `Sig.` / `Suggestive`。
3. **小样本设计必须配套误差率校准**：
   - 在每组 $n = 2\text{--}4$ 的在体 MFA 与多组学实验中，任何规则变更都应通过诱饵校准（`calibration.decoys`）、预设阳/阴性对照集（`panels.yaml` 的 `controls`）和逐只动物散点图（图 B）共同验证，防范假阳性膨胀。
4. **严格区隔 `pre_specified` 与 `post_hoc`**：
   - 任何基于结果调整后的面板或定向检验范围，必须在 `panels.yaml` 的 `selection` 字段如实记录为 `post_hoc` 及其依据，确保交付与发表可追溯。
