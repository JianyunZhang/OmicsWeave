# Phase 3 统计方法深化校准报告：`fdr_scope`、多示踪剂置换辅助 `Suggestive` 路径与 F 轴通量等效性检验

> **文档状态**：v1.0（Phase 3：S3 / S4 / S5 已落地并通过空假设模拟与真实数据置换校准）  
> **建立日期**：2026-11-23  
> **关联文档**：[01_mfa_flux_assessment_thresholds_and_improvements.md](01_mfa_flux_assessment_thresholds_and_improvements.md) · [02_mfa_multi_tracer_weighting_evaluation.md](02_mfa_multi_tracer_weighting_evaluation.md)  
> **评估基准数据**：`local/projects/dio_project/analysis_results/all_tissues/run_20261121_174318_7322`（单组学 MFA）与 `run_20261121_182112_fa44`（`reaction_dossier`，全库 4,200 条 `(tissue, contrast, flux_id)` 反应行）

---

## 1. 校准方法与设计

为满足 [01 号文档](01_mfa_flux_assessment_thresholds_and_improvements.md) 第 7 节对 Phase 3 的硬性验收要求（**“必须运行模拟空假设数据与组别置换，确认假阳性率受控，并保留基线对照”**），本次在落地 S3 / S4 / S5 前后并行执行了三套独立验证：

1. **高斯蒙特卡洛空假设与备择假设校准**：针对多示踪剂数 $k \in \{2, 3, 4\}$、单示踪剂每组样本量 $n \in \{3\text{v}3, 4\text{v}4\}$、真实标准化效应 $\delta \in \{0.0, 0.8, 1.0, 1.5\}$，每个组合重复 $M = 20{,}000$ 次独立反应（共 480,000 次模拟），精确测量各规则的空假设假阳性率（FPR）、备择检验功效（Power）与错误等效率（TOST Type I Error）。
2. **真实数据组别标签分层置换空假设（Real-Data Permutation Null）**：在 `dio_project` `HFD_vs_LFD` 全部 288 条可辨识多示踪剂反应上，于各示踪剂内部随机打乱 `HFD` / `LFD` 动物标签并重复 50 次全库多组织置换（共 $N = 14{,}400$ 条反应-置换实例），在保留真实偏态分布、拟合方差与化学计量通量组相关性的前提下直接测定经验假阳性率。
3. **基线配置 vs Phase 3 全库实测对照**：在 `dio_project` 全部 4,200 条 `(tissue, contrast, flux_id)` 行上逐条对比 Phase 2 基线与 Phase 3 变更。

---

## 2. S3：多重检验校正域 `fdr_scope`（`per_tissue` vs `per_contrast`）

### 2.1 设计与默认口径

- 在 `src/omicsweave/mfa/multi_tracer.py`（`PoolingSettings.fdr_scope`）与 `src/omicsweave/mfa/df_analysis.py`（`analysis.differential_flux.fdr_scope`）统一引入 `fdr_scope: per_tissue | per_contrast`，默认值为 **`per_tissue`**（与 `design.stratify_by: [tissue]`、转录组/蛋白组/代谢组单组织 BH-FDR 及 `assessment.py` 合取检验保持统计契约一致）。
- **透明双保留**：无论选用何种 `fdr_scope`，`differential.parquet`、`dossier_members.parquet` 与导出成员表始终同时保留 **`padj`**（当前作用域校正值）、**`padj_global`**（跨全组织校正值）与 **`fdr_scope`** 列，便于无预设全组织扫描对照。

### 2.2 空假设校验与 `dio_project` 对照

- **真实数据置换空假设（50 次 7 组织全库置换）**：`padj < 0.05` 的平均假阳性数在 `per_contrast` 下为 `0.000 / 次`，在 `per_tissue` 下为 `0.080 / 次`（折合单组织 FWER ≈ **1.1%**，严格低于 5%）。
- **`dio_project` 实测对照**：
  - `HFD_vs_LFD` 检验组数从跨 7 组织 **192** 组校正回各组织实际独立通量组数（`heart=31`, `liver=35`, `kidney=35`, `sq_adipose=33`, `gast_muscle=33`, `brain=21`, `plasma=4`）；
  - `heart/ENO_c & PGM_c`（`p = 0.0129`）的 `padj` 由 `0.960` 校正为 `0.347`，`gast_muscle/G3PP_c`（`ND_vs_LFD`, `p = 0.010`）由 `0.934` 校正为 `0.220`；
  - `dio_project` 上 `padj < 0.05` 仍为 0（零假阳性膨胀），同时在 $n \ge 4\text{v}4$ 多示踪剂设计下解除 7 组织叠加对单组织 `L3 / Sig.` 的过度压制。

---

## 3. S4：多示踪剂分层置换辅助 `Suggestive` 路径（Path B）

### 3.1 规则定义与配置

在 `src/omicsweave/integration/dossier/assessment.py`（`RULES`）新增 `flux_concordant_min_abs_g: 0.8`（设为 `None` 或 `>= flux_min_abs_g` 即可关闭 Path B，退回纯基线 Path A）。当反应未达 `changed (established)` 时，`suggestive (L2)` 判定为两条显式路径的并集：

- **Path A（原基线）**：有可辨识示踪剂、`|Hedges' g| >= 1.0` 且方向稳定（分层 LOO / Bootstrap 通过且非单鼠驱动）；
- **Path B（S4 新增，标 `multi_tracer_perm_suggestive`）**：
  1. 可辨识示踪剂数 $k \ge 2$；
  2. 全部可辨识示踪剂方向一致（`len(signs) == 1`）；
  3. 无显著示踪剂间异质性（标准化 Cochran $Q$ 检验 $p \ge 0.10$）；
  4. 分层置换 $p < \alpha$（默认 `0.05`）；
  5. 合并标准化效应量 `|Hedges' g| >= 0.8`（对应经典大效应量基准）；
  6. 通过分层合并稳定性门控且非单鼠驱动。

### 3.2 蒙特卡洛（$M = 20{,}000$）与真实数据置换校准结果

| 场景 | 空假设 $\delta=0$ Path A | 空假设 $\delta=0$ Path A\|B(0.8) | 空假设增量 $\Delta$FPR | 备择 $\delta=0.8$ Power 变化 | 备择 $\delta=1.0$ Power 变化 |
|---|--:|--:|--:|--:|--:|
| 模拟 $k=2, n=3\text{v}3$ | `6.73%` | `6.73%` | `+0.00%` | `27.9% → 27.9%` | `38.3% → 38.3%` |
| 模拟 $k=2, n=4\text{v}4$ | `4.15%` | `4.15%` | `+0.00%` | `28.7% → 28.7%` | `41.4% → 41.4%` |
| 模拟 $k=3, n=3\text{v}3$ | `2.62%` | `2.67%` | `+0.05%` | `22.9% → 23.4%` | `35.3% → 35.9%` |
| 模拟 $k=3, n=4\text{v}4$ | `1.34%` | `2.72%` | `+1.39%` | `23.0% → 33.9%` (**+10.9%**) | `38.9% → 52.1%` (**+13.1%**) |
| 模拟 $k=4, n=3\text{v}3$ | `0.92%` | `1.92%` | `+1.01%` | `19.1% → 28.0%` (**+8.8%**) | `32.9% → 44.2%` (**+11.2%**) |
| 模拟 $k=4, n=4\text{v}4$ | `0.38%` | `1.37%` | `+0.98%` | `19.3% → 34.4%` (**+15.1%**) | `36.3% → 54.7%` (**+18.4%**) |
| **`dio_project` 真实标签置换（$N=14{,}400$）** | **`2.58%`** | **`2.76%`** | **`+0.19%`** | — | — |

**统计结论**：在纯 Path A 下，随着示踪剂数 $k$ 从 2 增至 4，合并 `g` 抽样方差收敛反而令真实效应 $\delta \in [0.8, 0.99]$ 的检出率下降（从 `27.9%` 跌至 `19.1%`）；Path B 在将真实数据空假设假阳性增量控制在 **`+0.19%`**（总 FPR `2.76% < 5%`）的前提下，消除了“多示踪剂一致 + 置换 $p < 0.05$ 反而检出率更低”的统计悖论，将 $\delta = 0.8\text{--}1.0$ 的功效提升 **+8.8% 至 +18.4%**。

### 3.3 `dio_project` 实测新增 `Suggestive` 明细

全库 4,200 行中仅新增 **10 条反应（5 条独立代谢通量组）**，全部满足 $p \le 0.035$、`|\log_2\text{FC}| \ge 1.11$（2.2–3.7 倍变化）、全示踪剂同向且分层 LOO = 1.0：

- `sq_adipose / HFD_vs_LFD`: `BHB_input → BHB_trans → BHBDH_m → SCOT_m → ACACCOA_lyase_m`（5 条化学计量耦合酮体利用链，$k=2, n_{\min}=5, p=0.0333, |g|=0.81, \log_2\text{FC}=-1.87$）
- `gast_muscle / HFD_vs_LFD`: `PGM_c / ENO_c`（2 条糖酵解下段耦合链，$k=4, n_{\min}=11, p=0.0352, |g|=0.90, \log_2\text{FC}=-1.11$）
- `kidney / ND_vs_LFD`: `G3PP_c`（$k=2, n_{\min}=5, p=0.0350, |g|=0.93, \log_2\text{FC}=-1.87$）
- `gast_muscle / ND_vs_LFD`: `GK_c`（$k=2, n_{\min}=5, p=0.0300, |g|=0.91, \log_2\text{FC}=+1.79$）
- `brain / ND_vs_LFD`: `GPT_c`（$k=3, n_{\min}=8, p=0.0230, |g|=0.85, \log_2\text{FC}=-1.55$）

---

## 4. S5：F 轴通量级 TOST 等效性检验（`unchanged` vs `undetermined`）

### 4.1 五重统计护栏与配置

实证分析发现，仅靠 `pooled` 单格 `cell_state == "equivalent"`（`log2FC` 90% CI 落在 `±0.58` 内）不足以防止跨尺度假等效（例如心脏 `RPI_c` 因绝对均值差近零被标为 `equivalent`，但单示踪剂 `|g|` 高达 `1.62`）。因此 `assessment.py` 新增 `flux_equiv_g_bound: 1.0`（设为 `None` 或 `<= 0` 关闭），仅当同时满足全部 **5 条护栏** 时才允许 F 轴输出真正的 **`unchanged`**（实线 `—`）：

1. **相对尺度 TOST 通过**：`head["cell_state"] == "equivalent"`（`log2FC` 的 90% $t$-CI 严格位于 `±equivalence_bound` 即 `±0.58` 内）；
2. **无异质性且非单鼠驱动**：`not heterogeneous` 且 `not single_subject_driven`；
3. **无任意单示踪剂指向信号**：所有可辨识示踪剂均不满足 `|g_t| >= flux_min_abs_g`（`1.0`），且均非 `changed` / `nominal`；
4. **标准化尺度 TOST 通过**：合并标准化效应的 90% 上界 $|g_{\text{pooled}}| + t_{0.95,\text{df}} \cdot \text{SE}(g_{\text{pooled}}) < \text{flux\_equiv\_g\_bound}$（默认 `1.0`；优先取落盘 `smd_se`，兼容历史结果包时由可辨识示踪剂样本量与 $g_t$ 精确回算）；
5. **与 `decide()` 判定表联动一致**：当 F 轴已通过等效检验判为 `unchanged` 时，即使 C 轴酶量改变，综合结论也不再外推为 `L1`（酶预测通量变化），而归为 `L0` 并在 `phenomena` 中记录 **`flux_unchanged_despite_capacity`**（通量受代谢缓冲而未随酶量改变）。

### 4.2 蒙特卡洛 TOST 误差率校验与 `dio_project` 实测

- **TOST 第一类错误率（真实效应位于等效边界 $\delta = 1.0$ 时误判为 `unchanged` 的概率）**：
  - $k=2, n=3\text{v}3$: `0.00%`（样本量不足时拒绝宣称等效）
  - $k=2, n=4\text{v}4$: `2.28%`
  - $k=3, n=3\text{v}3$: `2.65%`
  - $k=3, n=4\text{v}4$: `4.02%`
  - $k=4, n=3\text{v}3$: `2.82%`
  - $k=4, n=4\text{v}4$: `3.18%`
  - 全部严格低于名义单侧显著性水平 **`5.0%`**；当真实大效应 $\delta = 1.5$ 时误判率仅 **`0.14%–0.24%`**。
- **真实零效应（$\delta = 0.0$）下的等效识别功效**：$k=3\text{--}4$ 时为 `20.5%–45.7%`（真实数据置换空假设下为 `12.0%`）。
- **`dio_project` 实测结果**：全库 855 条 `undetermined` 行中，仅 **8 条**高把握未变反应通过五重护栏转为 `unchanged`（其余 27 条仅 `pooled.cell_state == "equivalent"` 但不满足标准化 90% 上界或单示踪剂护栏的反应如 `RPI_c`、`SER_protein_input (HFD)` 均被安全拦截留在 `undetermined`）：
  - `heart / HFD_vs_LFD`: `PFK_c`、 `FBA_c`（$k=4, \log_2\text{FC}=+0.15\pm0.24, |g|=0.14, g_{90,\text{hi}}=0.89$；其中 `FBA_c` 伴随酶容量提示下调，记 `flux_unchanged_despite_capacity`）
  - `sq_adipose / HFD_vs_LFD`: `CS_m`（$k=4, \log_2\text{FC}=-0.01\pm0.31, |g|=0.01, g_{90,\text{hi}}=0.76$，记 `flux_unchanged_despite_capacity`）
  - `gast_muscle / HFD_vs_LFD`: `GLYCOGENOLYSIS_c`（$g_{90,\text{hi}}=0.95$）、`PGMT_c`（$g_{90,\text{hi}}=0.83$）、`Salvage_c`（$g_{90,\text{hi}}=0.81$）
  - `kidney / ND_vs_LFD`: `Salvage_c`（$g_{90,\text{hi}}=0.98$）；`plasma / ND_vs_LFD`: `SER_protein_input`（$g_{90,\text{hi}}=0.94$）

---

## 5. 回退基线与项目级敏感性配置示例

若需在特定项目中复现 Phase 2 严格基线进行敏感性对照，可在 `experiment.yaml` 中显式关闭 S3/S4/S5：

```yaml
omics:
  mfa:
    analysis:
      multi_tracer:
        fdr_scope: per_contrast        # 默认 per_tissue；设为 per_contrast 回退跨组织联合校正

integration:
  parameters:
    reaction_dossier:
      rules:
        flux_concordant_min_abs_g: null  # 默认 0.8；设为 null 或 1.0 关闭 S4 路径 B
        flux_equiv_g_bound: null         # 默认 1.0；设为 null 或 0 关闭 S5 F 轴等效判定
```
