# 多组学集成分析流程技术文档

**版本**: v2.0 (2026-04-27)  
**范围**: `omicsweave/multiomics/` 完整管道 + MFA 映射表 (`mappings/`)  
**关联管道**: `omicsweave/mfa/`, `omicsweave/metabolomics/`, `omicsweave/proteomics/`, `omicsweave/transcriptomics/`

> **历史文档警告（2026-08-12）**：本文描述 v2 旧矩阵，仅用于追溯 `legacy_*` 输出，不代表当前主流程。尤其是“FDR + Cohen's d 能缓解或修复 MFA 伪重复”的结论已弃用：效应量方差若仍来自 MCMC 拟合点，就没有恢复正确的生物学独立观测单位。当前实现见 `README.md`、`integration/evidence_profile.py` 和 `integration/stability.py`。

---

## 一、总体流程架构

管道采用 **MFA 通量为中心的 6 步**设计：

```
加载配置 → 加载四组学数据 → 计分 → 通路收敛分析 → 绘图 → 报告
```

**设计思想**：代谢通量是表型终结点，以 MFA 差异通量作为"入口闸门"筛选后，再用转录组/蛋白组/代谢组做方向一致性验证。通过闸门的反应进入四层评分矩阵，未通过的标记为 Insufficient。

### 运行方式

```bash
# 通过实验 ID 运行
omicsweave integration run \
    --experiment 20250425_cimr_wt_mice_dio

# 或直接指定配置文件
omicsweave integration run \
    --config experiments/20250425_cimr_wt_mice_dio/multiomics_integration_setting/config.yaml
```

### 模块结构

```
omicsweave/multiomics/
  config.py                              # MultiOmicsConfig dataclass + YAML 加载器
  mappings/
    flux_gene_mapping.py                 # 87 反应 → 基因映射
    flux_protein_mapping.py              # 63 反应 → 酶/蛋白质映射
    flux_metabolite_mapping.py           # ~70 反应 → 代谢物映射（底物/产物/关联 分层）
  loaders/
    mfa_loader.py                        # MFA 差异通量 + 通路 + TCA + 概览
    transcriptomics_loader.py            # DE + GSEA + Enrichr
    proteomics_loader.py                 # 蛋白质 DE + 通路（含比较名翻转）
    metabolomics_loader.py               # 代谢物 DE + 通路
  integration/
    flux_omics_scorer.py                 # 改进打分矩阵 + MFA gate + 置信度分级
    pathway_convergence.py               # 跨组学通路收敛分析（含 BH FDR 校正）
  plotting/
    _config.py                           # 绑图全局状态（统一阈值、颜色、标签）
    overview_plots.py                    # 数据可用性热图 + 差异汇总 + 组织看板
    volcano_plots.py                     # MFA 火山图
    heatmap_plots.py                     # 通路热图 + TCA + 收敛 + 通量概览
    scoring_plots.py                     # 评分分布 + 置信度饼图
  report.py                              # Markdown 报告生成
  run_multiomics_analysis.py             # CLI 入口
```

---

## 二、评分模型总览

### 四层评分结构

| 层级 | 满分 | 权重 | 显著性 | 方向 | 效应量 |
|------|------|------|--------|------|--------|
| MFA (gate + layer) | 1.0 | 1.0 | 0.5 | — | 0–0.5 (Cohen's d) |
| 转录组 | 1.0 | 1.0 | 0.4 | 0–0.3 | 0–0.3 (log2FC) |
| 蛋白组 | 1.2 | 1.2 | 0.5 | 0–0.3 | 0–0.4 (log2FC) |
| 代谢组 | 0.8 | 0.8 | 0.3 | 0–0.3 | 0–0.2 (log2FC) |
| **总计** | **4.0** | | | | |

### 置信度分级

| 等级 | 总分范围 | 含义 |
|------|---------|------|
| **High** | ≥ 3.0 | 至少 3 个组学层提供一致性证据（满分 4.0 的 75%） |
| **Medium** | ≥ 2.0 | 至少 2 个组学层提供支持证据 |
| **Low** | ≥ 1.0 | 仅 MFA 层通过闸门但其他组学证据不足 |
| **Insufficient** | < 1.0 | 未通过 MFA 闸门，或通过但各层均无证据 |

### 权重设计原理

| 组学 | 权重 | 原理 |
|------|------|------|
| 蛋白组 | **1.2** | 最高权重。蛋白质是基因功能的直接执行者，与代谢通量物理距离最近（酶-反应对应关系） |
| 转录组 | 1.0 | 标准权重。mRNA 变化反映转录调控，但受翻译效率/蛋白稳定性等下游环节影响 |
| MFA | 1.0 | 标准权重。通量是最直接的代谢表型测量，但 MFA gate 已在入口处做了筛选 |
| 代谢组 | **0.8** | 最低权重。代谢物池浓度受多种反应汇聚/分流影响，单一反应→代谢物的对应关系特异性较低 |

---

## 三、MFA Gate（入口闸门）

所有 MFA 差异反应必须通过闸门才能进入后续评分。未通过的反应直接标记为 `Insufficient`，总分为 0。

### 闸门条件

| 路径 | 条件 | 触发场景 |
|------|------|---------|
| **主路径** | `FDR < 0.05` AND `\|Cohen's d\| > 0.8` | 默认路径，FDR 和 Cohen's d 列可用时 |
| Fallback | `pvalue < 0.05` AND `\|log2FC\| > 0.58` | 当 FDR 或 Cohen's d 列缺失时 |

### 为什么使用 Cohen's d

MFA 使用 ~150 次 bootstrap 拟合产生伪重复数据点，将实际 n=3 的生物学重复膨胀为 ~450 个统计样本。原始 p 值因此被极度压缩（72-82% 的通量 p<0.05），无法有效区分真阳性与统计假象。

Cohen's d 是标准化效应量（均值差/合并标准差），不受样本量影响。伪重复增加了数据点数但不改变效应量的真实大小，因此可有效过滤掉微小变化的假阳性反应。

| 指标 | 对伪重复的敏感度 | 作用 |
|------|---------------|------|
| P-value | **高度敏感** — n 越大 p 越小 | 判断差异是否存在 |
| FDR | 依赖 P-value | 控制多重比较假阳性率 |
| Cohen's d | **不敏感** — 与 n 无关 | 判断差异实际有多大 |

三者组合使用：FDR 确保差异确实存在（非随机噪声），Cohen's d 确保差异足够大（非统计假象）。

---

## 四、评分细则

### 4.1 MFA 评分层（0–1.0 分）

通过闸门后，MFA 本身作为独立评分层。

| 子项 | 满分 | 计算规则 |
|------|------|---------|
| **显著性 (sig_score)** | 0.5 | `FDR < 0.05` → 0.5；否则检查 `pvalue < pval_cutoff` → 0.5；均不满足 → 0.0 |
| **效应量 (effect_score)** | 0–0.5 | 基于 Cohen's d 分档：`\|d\| ≥ 1.5` → 0.5（大效应）；`0.8 ≤ \|d\| < 1.5` → 0.25（中等效应）；`< 0.8` → 0.0。Fallback 时使用 log2FC 分档 × 0.5 |

> **历史说明**: v1 版本中 MFA 层包含"生物学注释"子项（bio_score 0.4 分），基于映射表中是否存在 `pathway`/`description` 字段。该子项在 v2 中被移除，原因是 96% 的 gate-passed 反应都获得了该分数，区分度极低。0.4 分被平均分配给 sig_score（0.3→0.5）和 effect_score（0.3→0.5），使评分完全由统计指标驱动。

### 4.2 转录组评分层（0–1.0 分）

使用 `flux_gene_mapping.py` 映射反应→基因列表。

**筛选标准**: `padj < 0.05` AND `|log2FC| > 0.58`

| 子项 | 满分 | 计算规则 |
|------|------|---------|
| **显著性 (sig_score)** | 0.4 | 至少一个映射基因显著**且方向一致**时得 0.4 |
| **方向一致性 (dir_score)** | 0 / 0.3 | 基因 log2FC 方向与通量 log2FC 方向相同 → 0.3；若所有显著基因均为反方向 → **整层归零** |
| **效应量 (effect_score)** | 0–0.3 | 取方向一致基因中最大 `\|log2FC\|`，分档：`< cutoff` → 0；`cutoff ~ 2×cutoff` → 0.15；`≥ 2×cutoff` → 0.3 |

**反向归零规则**: 如果所有显著基因的 fold-change 方向与通量方向相反（如通量上调但所有基因下调），转录组整层得分强制为 0。反方向证据被视为对该通量变化的**反驳**而非支持。

**匹配方式**: 基因名精确匹配（大小写敏感），从 `flux_gene_mapping` 的 `genes` 字段获取映射列表。

### 4.3 蛋白质组评分层（0–1.2 分）

使用 `flux_protein_mapping.py` 映射反应→蛋白质列表，区分直接酶和间接蛋白。

**筛选标准**: `padj < 0.05` AND `|log2FC| > 0.58`

| 子项 | 满分 | 计算规则 |
|------|------|---------|
| **显著性 (sig_score)** | 0.5 | 至少一个映射蛋白显著时得 0.5（无论方向） |
| **方向一致性 (dir_score)** | 0 / 0.15 / 0.3 | **直接酶 (direct_enzymes)** 方向一致 → 0.3；间接蛋白 (upstream_regulators / transporters / cofactors) 方向一致 → 0.15。优先检查直接酶 |
| **效应量 (effect_score)** | 0–0.4 | 取方向一致蛋白中最大 `\|log2FC\|`，分档：`< cutoff` → 0；`cutoff ~ 2×cutoff` → 0.2；`≥ 2×cutoff` → 0.4 |

**与转录组的区别**:
- 蛋白组**不会因反向而整层归零**——翻译后修饰可能导致蛋白丰度与通量方向不一致但仍有生物学意义，因此 sig_score 保留
- 直接酶权重是间接蛋白的 2 倍（0.3 vs 0.15），因为直接催化反应的酶与通量因果关系最强

**匹配方式**: 蛋白名精确匹配，区分 `direct_enzymes` 和间接类别。

### 4.4 代谢组评分层（0–0.8 分）

使用 `flux_metabolite_mapping.py` 映射反应→代谢物列表，映射结构为 `{substrates, products, related}`。

**筛选标准**: `pval < 0.05` AND `|log2FC| > 0.58`（有 `fdr` 列时使用 `FDR < 0.05`）

| 子项 | 满分 | 计算规则 |
|------|------|---------|
| **显著性 (sig_score)** | 0.3 | 至少一个映射代谢物显著时得 0.3（无论方向） |
| **方向一致性 (dir_score)** | 0 / 0.15 / 0.3 | **产物**方向一致 → 0.3；**底物**方向一致 → 0.3；**关联**仅显著 → 0.15。优先级：products > substrates > related |
| **效应量 (effect_score)** | 0–0.2 | 取方向一致代谢物中最大 `\|log2FC\|`，分档：`< cutoff` → 0；`cutoff ~ 2×cutoff` → 0.1；`≥ 2×cutoff` → 0.2 |

**底物/产物方向规则**:

```
flux ↑ (通量上调)
  → products 丰度 ↑ (产物增加, 同向一致) → dir_score = 0.3
  → substrates 丰度 ↓ (底物被消耗, 反向一致) → dir_score = 0.3
  → related 仅显著变化 → dir_score = 0.15

flux ↓ (通量下调)
  → products 丰度 ↓ (产物减少, 同向一致) → dir_score = 0.3
  → substrates 丰度 ↑ (底物积累, 反向一致) → dir_score = 0.3
  → related 仅显著变化 → dir_score = 0.15
```

**历史匹配方式**: 本节所述名称子串匹配已经弃用。当前实现使用 HMDB/KEGG/ChEBI
标准 ID 优先、规范化完整名称精确回退，并隔离跨角色歧义；详见 `README.md` 和
`mappings/feature_resolver.py`。

### 4.5 效应量分档函数（共用）

所有组学层共用 `_effect_tier()` 函数，将 |log2FC| 映射为三档：

```
|log2FC| < cutoff           → tier = 0.0 (trivial, 微不足道)
cutoff ≤ |log2FC| < 1.0     → tier = 0.5 (moderate, 中等)
|log2FC| ≥ 1.0              → tier = 1.0 (large, 大效应, 固定为 2 倍变化)

cutoff 默认 = 0.58 (即 |log2FC| < 0.58 → 0; 0.58~1.0 → 0.5; ≥ 1.0 → 1.0)
```

各层实际 effect_score = tier × 层效应量系数:

| 组学 | 系数 | 最大 effect_score | 阈值依据 |
|------|------|-----------------|---------|
| MFA | 0.5 | 0.5 | Cohen's d: 0.8/1.5 |
| 转录组 | 0.3 | 0.3 | log2FC: 0.58/1.0 |
| 蛋白组 | 0.4 | 0.4 | log2FC: 0.58/1.0 |
| 代谢组 | 0.2 | 0.2 | log2FC: 0.58/1.0 |

---

## 五、阈值体系

### 全局阈值（配置驱动）

以下阈值从 `config.yaml` → `scoring` 节读取，通过 `_config.py` 统一分发到积分器、仪表板和绘图模块：

| 参数 | 默认值 | 配置入口 | 使用位置 |
|------|--------|---------|---------|
| `pval_cutoff` | 0.05 | `scoring.pval_cutoff` | 积分器、仪表板、差异汇总 |
| `log2fc_cutoff` | 0.58 | `scoring.log2fc_cutoff` | 积分器、仪表板、火山图、差异汇总 |
| `log2fc_tiers` | [0.58, 1.0] | `scoring.effect_tiers.log2fc` | 转录组/蛋白组/代谢组效应量分档 |
| `cohen_d_tiers` | [0.8, 1.5] | `scoring.effect_tiers.cohen_d` | MFA 层效应量分档 |
| 置信度 High | 3.0 | `scoring.confidence_thresholds.high` | 积分器置信度分级 |
| 置信度 Medium | 2.0 | `scoring.confidence_thresholds.medium` | 积分器置信度分级 |
| 置信度 Low | 1.0 | `scoring.confidence_thresholds.low` | 积分器置信度分级 |

### 效应量分档（配置驱动）

效应量分档阈值通过 `config.yaml` → `scoring.effect_tiers` 配置，格式为 `[moderate下限, large下限]`：

```yaml
effect_tiers:
  log2fc: [0.58, 1.0]     # |log2FC| < 0.58 → trivial; 0.58~1.0 → moderate; ≥1.0 → large
  cohen_d: [0.8, 1.5]     # |Cohen's d| < 0.8 → trivial; 0.8~1.5 → moderate; ≥1.5 → large
```

- `log2fc_tiers` 用于转录组/蛋白组/代谢组的 `_effect_tier()` 函数
- `cohen_d_tiers` 用于 MFA 层的 Cohen's d 效应量评分
- 修改配置文件即可自定义阈值，无需改动代码

### MFA 专用阈值（硬编码）

| 参数 | 值 | 位置 | 说明 |
|------|------|------|------|
| MFA gate FDR | 0.05 | `flux_omics_scorer._mfa_gate_passed()` | 积分器入口闸门 |
| MFA gate Cohen's d | 0.8 | `flux_omics_scorer._mfa_gate_passed()` | 效应量闸门 |
| 仪表板 MFA | FDR<0.05 + \|d\|>0.8 | `overview_plots.plot_tissue_dashboard()` | 与 gate 一致 |

### 其他阈值

| 参数 | 值 | 位置 | 说明 |
|------|------|------|------|
| GSEA FDR | 0.25 | `pathway_convergence.py` | 标准 GSEA 实践 |
| 蛋白 Q-value 预过滤 | 0.01 | `proteomics/run_proteomics_pipeline.py` | 蛋白组加载阶段 |
| 代谢组通路收敛 | BH FDR < 0.05 | `pathway_convergence.py` | Benjamini-Hochberg 校正 |

### 跨管道 log2FC 统一说明

| 管道 | log2FC cutoff | 对应 FC | 配置文件 |
|------|--------------|---------|---------|
| 多组学集成 | **0.58** | ~1.50x | `multiomics_integration_setting/config.yaml` → `scoring.log2fc_cutoff` |
| 代谢组独立分析 | **0.58** | ~1.50x | `metabolomics_data/config.yaml` → `params.log2fc_cutoff` |
| 蛋白组独立分析 | **0.58** | ~1.50x | `proteomics_data/config.yaml` → `params.log2fc_cutoff` |

**设计决策**: v2 版本统一为 0.58（~1.5 倍变化），与单组学管道保持一致。理由：

1. **统计口径一致性**：同一概念（log2FC 显著性）在所有模块使用相同阈值，方法学上更自洽
2. **MFA gate 已提供核心筛选**：通过 MFA gate（FDR < 0.05 + |Cohen's d| > 0.8）的反应已具备强统计证据，各 omics 层的 log2FC cutoff 只影响评分细节
3. **实际影响有限**：仅 4.4% 的 gate-passed 反应 MFA log2FC 在 0.5-0.58 区间，log2FC cutoff 主要影响效应量子项（满分 0.2-0.5），方向一致性和 MFA gate 不受影响

> **历史说明**: v1 版本多组学集成使用 0.5，单组学使用 0.58，基于"多层交叉验证可适度放宽单一层阈值"的考虑。v2 中统一为 0.58 以简化统计框架。

---

## 六、MFA 映射表设计

### 6.1 `flux_gene_mapping.py` — 基因映射

**结构**: `Dict[str, {genes, pathway, description}]`

- 覆盖所有主要代谢通路（TCA、糖酵解、氨基酸代谢等，共 87 个反应）
- 提供 `pathway` 和 `description` 字段用于 MFA 层的生物学注释评分
- 基因名精确匹配（大小写敏感）

### 6.2 `flux_protein_mapping.py` — 蛋白映射

**结构**: `Dict[str, {direct_enzymes, upstream_regulators, downstream_enzymes, transporters, cofactors, description}]`

- 分类为 `direct_enzymes` / `upstream_regulators` / `transporters` 等
- 积分器据此对直接酶赋予更高方向权重（0.3 vs 0.15）
- 蛋白名使用基因符号

### 6.3 `flux_metabolite_mapping.py` — 代谢物映射

**结构**: `Dict[str, {substrates: [...], products: [...], related: [...]}]`

- 每个反应的代谢物分为三类：底物（被消耗）、产物（被生成）、关联（辅因子/上下游）
- 评分器根据底物/产物角色做正确的方向一致性判断：
  - 产物：通量↑ → 产物丰度↑（同向一致）
  - 底物：通量↑ → 底物丰度↓（反向一致）
  - 关联：仅要求显著变化，无方向期望
- 历史子串匹配已弃用；当前使用标准 ID 优先、完整名称精确回退和歧义隔离
- 辅助函数：`get_metabolite_roles()` 返回结构化映射，`get_related_metabolites()` 向后兼容返回扁平列表

### 6.4 映射表已知局限

| 局限 | 影响 | 计划 |
|------|------|------|
| 部分反应缺少完整基因角色 | 已有 `catalytic/transporter` 与 `regulatory/transcription_factor` 会在 EvidenceProfile 中按直接/间接映射质量处理；无角色者按间接映射 | 持续补全并人工审核角色注释 |
| 本地 UniProt 注册表仍有覆盖缺口 | 缺少 accession 的目标只能按规范符号回退 | 持续用官方 mouse reviewed 条目补全 |
| 未区分同工酶的组织特异性 | 肌肉用 HK，肝脏用 GCK，映射表不区分 | 添加组织特异性标注 |
| 部分代谢物缺 HMDB/KEGG/ChEBI | 只能进行规范化完整名称精确匹配 | 建立供应商 feature ID 到标准 ID 的策展表 |
| 映射表缺少版本号和来源引用 | 无法追溯数据库版本 | 添加 KEGG/MetaCyc 版本注释 |

---

## 七、通路收敛分析

跨组学通路收敛分析识别在 ≥2 个组学层中同时富集的通路术语。

### 方法

1. **代谢组通路**: 从各组织的 `Pathway_*.csv` 提取显著通路，使用 Benjamini-Hochberg FDR 校正（对每个 tissue×comparison 组合内独立应用）
2. **蛋白组通路**: 从 `pathway_summary_all.csv` 提取，优先使用 FDR 列过滤
3. **转录组 GSEA**: 从 GSEA 结果提取 `FDR q-val < 0.25` 的通路（标准 GSEA 实践）

对每个组织，计算所有组学对的通路术语交集，输出共享术语列表和计数汇总。

### BH FDR 校正

`PathwayConvergenceAnalyzer._bh_fdr()` 实现标准 Benjamini-Hochberg step-up 流程：

```python
# 排序 → 计算 p×n/rank → 从后向前取 min 保证单调性 → clip [0, 1]
fdr[i] = min(fdr[i], fdr[i+1])  # step-up 单调性
```

无需外部依赖，在 numpy 上直接实现。

---

## 八、可视化模块

所有绘图模块从 `_config.py` 读取统一的阈值和配色，避免硬编码不一致。

| 模块 | 产出 | 阈值来源 |
|------|------|---------|
| `overview_plots.plot_data_availability()` | 数据可用性热图 | 使用 MFA gate (FDR+Cohen's d) 计数 |
| `overview_plots.plot_differential_summary()` | 差异特征柱状图 | `cfg.PVAL_CUTOFF` / `cfg.LOG2FC_CUTOFF` |
| `overview_plots.plot_tissue_dashboard()` | 组织整合看板 | 代谢物/蛋白: `cfg.PVAL_CUTOFF` + `cfg.LOG2FC_CUTOFF`; MFA: FDR+Cohen's d gate |
| `volcano_plots.plot_mfa_volcanoes()` | MFA 火山图 | `pval_cutoff=0.05`, `fc_cutoff=0.58` |
| `heatmap_plots` | 通路热图/TCA/通量概览 | MFA 数据驱动 |
| `scoring_plots.plot_score_distribution()` | 评分分布图 | 置信度阈值 3.0/2.0/1.0 |

---

## 九、配置文件格式

配置采用 `config.yaml` + `groups.yaml` 双文件结构，位于 `experiments/{experiment_id}/multiomics_integration_setting/`，与转录组管道的配置模式一致：

**`groups.yaml`** — 组别与小鼠映射（与 `transcriptomics_data/groups.yaml` 保持一致）：

```yaml
conditions:
  hfd_14w:
    display: "HFD"
    diet: "D12492"
    diet_desc: "High-fat diet DIO 14 week"
    color: "#FF8000"
  lfd_14w:
    display: "LFD"
    diet: "D12450J"
    diet_desc: "Low-fat diet control 14 week"
    color: "#00A550"
  hsd_14w:
    display: "HSD"
    diet: "AIN-76A"
    diet_desc: "High-sugar purified diet 14 week"
    color: "#800080"

mouse_groups:
  "52": hfd_14w
  "53": hfd_14w
  "54": hfd_14w
  "63": lfd_14w
  "64": lfd_14w
  "65": lfd_14w
  "73": hsd_14w
  "74": hsd_14w
  "75": hsd_14w
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
      log2fc: [0.58, 1.0]
      cohen_d: [0.8, 1.5]
    confidence_thresholds: {high: 3.0, medium: 2.0, low: 1.0}
  tissue_labels:
    heart: "心脏 Heart"
    liver: "肝脏 Liver"
    kidney: "肾脏 Kidney"
    gast_muscle: "腓肠肌 Gastrocnemius"
    sq_adipose: "皮下脂肪 SQ Adipose"
```

---

## 十、输出文件

```
analysis_results/{experiment_id}/multiomics_results/
  tables/
    flux_omics_scores.csv          # 每个反应 × 每个比较的评分明细
    pathway_scores.csv             # 按通路聚合的评分
    pathway_convergence.csv        # 跨组学共享通路术语
    pathway_convergence_counts.csv # 共享通路计数
    multi_omics_overview.csv       # 各组织各组学显著特征数
    differential_summary.csv       # 各比较的差异特征汇总
    tissue_dashboard_{comp}.csv    # 各比较的组织看板数据
  plots/
    multi_omics_overview.png       # 数据可用性热图
    differential_summary_{comp}.png
    tissue_dashboard_{comp}.png
    mfa_volcano_{comp}.png
    mfa_pathway_condition_heatmap.png
    tca_cycle_*.png
    tissue_flux_overview.png
    convergence_heatmap.png
    score_distribution.png
    confidence_pie.png
  multiomics_report.md             # 分析报告
```

### 输出列说明

`flux_omics_scores.csv` 主要列：

| 列名 | 含义 |
|------|------|
| `flux_id` | MFA 反应 ID |
| `tissue` | 组织 |
| `comparison` | 比较名称 |
| `flux_log2fc` | 通量 log2 fold change |
| `flux_pvalue` | 通量原始 p 值 |
| `flux_fdr` | 通量 FDR 校正 p 值 |
| `cohen_d` | Cohen's d 效应量 |
| `mfa_score` | MFA 层得分 (0–1.0) |
| `tx_score` | 转录组层得分 (0–1.0) |
| `prot_score` | 蛋白组层得分 (0–1.2) |
| `metab_score` | 代谢组层得分 (0–0.8) |
| `total_score` | 总分 (0–4.0) |
| `confidence` | 置信度 (High/Medium/Low/Insufficient) |
| `gate_passed` | 是否通过 MFA gate |
| `*_details` | 各层评分明细（JSON 格式） |

---

## 十一、当前实验结果

实验 ID: `20250425_cimr_wt_mice_dio` (HFD/HSD/LFD, 5 组织, n=3/group)

| 指标 | 数值 |
|------|------|
| 总反应数 | 1191 |
| 通过 MFA gate | 248 (20.8%) |
| High (≥3.0) | 1 |
| Medium (≥2.0) | 41 |
| Low (≥1.0) | 169 |
| Insufficient (<1.0) | 980 |
| 跨组学收敛通路 | 11 个共享术语 (3 个组织) |

---

## 十二、已知局限与未来优化

| 优先级 | 问题 | 实施建议 |
|--------|------|---------|
| 中 | MFA gate fallback 路径过宽 | 将 fallback 的 log2fc_cutoff 提升至 1.0，或确保 `fdr_mannwhitney` 和 `cohen_d` 列总是存在 |
| 中 | 代谢物名称未标准化 | 使用 HMDB ID 或 KEGG Compound ID 作为映射键 |
| 中 | 映射表缺少组织特异性 | 基因/蛋白映射中标注组织特异性同工酶 |
| 中 | 蛋白组通路数据缺失 | 在蛋白组 pipeline 中增加通路汇总输出 |
| 低 | 映射表缺少版本号 | 添加数据库来源和版本注释 |
| 低 | 映射覆盖率未验证 | 编写覆盖率测试脚本 |
| 低 | 蛋白名未映射 UniProt ID | 增加 symbol→UniProt 转换层 |
