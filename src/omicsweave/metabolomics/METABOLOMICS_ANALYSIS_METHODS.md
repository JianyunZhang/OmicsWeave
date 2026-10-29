# 代谢组学分析方法说明

本文档描述 `omicsweave/metabolomics/` 当前代码实际执行的方法，可作为论文 Methods 的技术底稿。最终投稿时仍应补充仪器、色谱、离子源、采集模式、峰识别、定量和代谢物鉴定置信等级等上游信息。

## 1. 数据类型与分析单位

输入为供应商处理后的代谢物峰面积表，行表示代谢物或代谢特征，列表示生物样本、pooled QC 或注释。流程分析相对丰度，不将峰面积解释为绝对浓度。统计推断单位是生物学样本；QC 样本只用于技术质量评估和特征过滤，不进入组间检验。

样本通过 `groups.yaml` 中的 mouse ID 映射到实验条件。支持 `<mouse_id>_<tissue>` 和旧版 `<mouse_id>-<tissue>`。未配置样本会导致流程停止，而不会归入 Unknown 组。

## 2. 特征过滤和 QC

零值、负值和非数值首先视为缺失。对每个代谢物分别计算各实验组的检出率；只要该代谢物在至少一个组达到 `detection_rate`，即可保留。这一规则允许保留只在某个处理组检出的 on/off 代谢物。

存在至少两个 pooled QC 样本且 `qc_filter_enabled` 开启时，计算：

\[
RSD_j = \frac{SD(QC_j)}{Mean(QC_j)}
\]

代谢物需同时满足 QC 检出率不低于 `qc_min_detection_rate`，且 QC RSD 不高于 `qc_rsd_cutoff`。默认阈值分别为 80% 和 30%。所有特征的 QC 均值、标准差、RSD、检出率和 `QC_Pass` 都写入 `qc_feature_metrics.csv`，包括被过滤特征。

若无足够 QC，流程不会伪造 RSD 结论，而是在 summary 中记录 `qc_filter_applied: false`。

## 3. 缺失值插补

通过检出率和 QC 门控后，剩余缺失值按代谢物自身的最小正值插补：

\[
x_{ij}^{imp} = f \times \min_{k:x_{kj}>0}(x_{kj})
\]

其中默认 `imputation_fraction = 0.5`。每个代谢物使用自身尺度，不使用全矩阵最小值，因此不会把高丰度代谢物缺失值压到其他低丰度代谢物的数量级。

## 4. 样本归一化与转换

默认使用 probabilistic quotient normalization（PQN）。先以所有分析列逐代谢物中位数构建参考谱，再计算每个样本相对参考谱的逐代谢物商；样本的归一化因子是这些有限商的中位数：

\[
q_{ij}=\frac{x_{ij}}{r_i},\qquad
s_j=median_i(q_{ij}),\qquad
x^{norm}_{ij}=\frac{x_{ij}}{s_j}
\]

配置还支持 median、TIC 和 none。归一化后执行 log2 转换。每个样本的 factor、归一化前后总信号均落盘，便于审计纯上样量偏移是否被消除。

PCA 和 PLSRegression 内部的 feature scaling 仅服务于多变量模型，不能替代上述样本间归一化。

## 5. PCA

PCA 对 log2、样本归一化后的矩阵进行特征标准化，再计算前两个主成分。若 pooled QC 存在，其坐标会与生物样本一同显示。PCA 是无监督探索，用于识别离群点、技术漂移和总体谱形，不作为组间差异检验。

## 6. PLS-DA 与 VIP

多分类 PLS-DA 将类别编码为 one-hot Y 矩阵，而非把类别标签编码为 0、1、2 的有序连续变量。模型组件数同时受样本数、特征数和 `类别数−1` 限制。

每个模型报告：

- 训练集 R²Y；
- 分层 K-fold 交叉验证 Q²；
- 交叉验证分类准确率；
- 标签置换检验 P 值；
- 实际完成的有效置换次数。

默认运行 200 次标签置换。任一类别少于 2 个样本时，模型标记为 `insufficient_validation`，Q²/置换结果和供差异判定使用的 VIP 均为 NaN。

除全局 one-hot 模型外，每个预定义 comparison 单独拟合二分类模型并计算 `VIP_<comparison>`。VIP 不跨 comparison 复用。VIP 计算基于各潜变量解释的 Y 方差和 X 权重。

默认主显著性不要求 VIP，因为高维小样本监督模型可能不稳定。只有显式配置 `require_vip_for_significance: true` 时，VIP 才进入主判定；此时未通过模型验证的比较不会产生 VIP 驱动的显著结果。

## 7. 差异代谢物

差异检验在归一化、log2 转换后的数据上运行。每个 comparison 使用向量化 Welch 独立样本 t 检验（`equal_var=False`）：

\[
log_2FC_i=\overline{x}_{i,A}-\overline{x}_{i,B}
\]

其中 comparison 名称的左侧组为分子、右侧组为分母。推荐在 `contrast_groups` 中显式定义方向；自动比较保留 YAML 条件顺序，并把标记为 `reference: true` 或 `role: control` 的条件置于分母。

有限 P 值使用所配置方法校正；`bh` 规范化为 Benjamini–Hochberg (`fdr_bh`)。NaN 不计入检验数并保持原位。默认主判定为：

\[
FDR < 0.05 \quad \land \quad |log_2FC| > 0.58
\]

结果同时输出：

- `NominalSignificant_*`：原始 P + fold change；
- `SignificantUnivariate_*`：FDR + fold change；
- `SignificantWithVIP_*`：FDR + fold change + comparison-specific VIP；
- `Significant_*`：唯一主判定；
- `Evidence_*`：inferential 或 insufficient。

任一组少于 `min_observations_per_group`（默认 2）时，仅保留描述性 log2FC，P/FDR 为 NaN，主显著性为 False。

## 8. 火山图和热图

火山图纵轴默认使用 `-log10(FDR)`，颜色和标签直接消费 `Significant_<comparison>`，因此与 DM CSV、ORA 和报告一致。

热图默认按全部生物样本上的方差选择 Top N，不使用组标签，是无监督展示。每个代谢物随后做行 Z-score 并进行层次聚类。用户可明确选择 `heatmap_selection_method: anova`，但该图会标记为 supervised ANOVA，不能作为独立的组间分离证据。

## 9. KEGG ORA

ORA 直接使用 `Significant_<comparison>`。统计 universe 是具有至少一个合格 KEGG pathway 注释的已测代谢物；显著集先与同一 universe 取交集，因此列联表不会混用两个背景。

对每个符合大小阈值的通路建立：

| | 通路内 | 通路外 |
|---|---:|---:|
| 主显著代谢物 | a | b |
| 非显著背景 | c | d |

使用 `fisher_exact(..., alternative="greater")` 检验过度代表。所有合格通路，包括 overlap 为 0 的通路，都进入多重检验分母；输出 P、FDR、odds ratio、overlap、背景大小和成员列表。

默认 `pathway_scope: metabolism`，排除信号、疾病通路及 ko01100/ko01110/ko01120 总览图，并用 3–500 的集合大小门控。设置 `all` 可放宽类别，但需谨慎解释共享代谢物造成的疾病图富集。

## 10. GSEA

GSEA 默认以 Welch t 统计量对全部有限代谢物排序。若旧结果没有 t 统计量，则回退到带方向的 `-log10(P)`：

\[
rank_i=sign(log_2FC_i)\times[-log_{10}(P_i)]
\]

NaN 和无穷值在传给 gseapy 前移除；零值保留。通路成员先与排名 universe 相交，再执行集合大小过滤。默认使用同一 KEGG metabolism 范围、固定随机种子和 1000 次置换，输出 NES、nominal P 和 FDR q-value。

## 11. 报告口径

自动报告将 FDR 主发现置于首位，并单独披露原始 P 名义发现数。Top 代谢物仅从主显著集合中选择。ORA 以通路 FDR 计数；PLS-DA 同时显示 Q²、交叉验证准确率与置换 P，不根据 score plot 视觉分离直接下结论。

## 12. 批次与方法限制

当前版本不在峰面积矩阵上自动执行 ComBat/SVA。批次与处理混杂时，事后校正无法区分技术效应和真实生物效应。配置只允许 `batch_correction: none`；多批次建议在平衡设计下分别分析，再对效应量进行结果级整合。

流程也不替代：

- 原始谱图峰提取、去卷积、对齐和 blank subtraction；
- 同位素内标归一化和绝对定量；
- MS/MS 与保留时间支持的代谢物鉴定等级；
- 针对重复测量、配对样本或复杂协变量的线性/混合模型。

## 13. 主要软件

- pandas / NumPy：表格与矩阵处理；
- SciPy：Welch t 和 Fisher exact；
- statsmodels：多重检验校正；
- scikit-learn：PCA、PLSRegression、分层交叉验证；
- gseapy：pre-ranked GSEA；
- matplotlib / seaborn：图形；
- openpyxl / PyYAML：输入与配置。

## 14. 方法学参考

1. Parsons HM et al. Comparison of normalization strategies for metabolomic data. *Analytical Chemistry*. 2009;81:7668–7675.
2. Dieterle F et al. Probabilistic quotient normalization as robust method to account for dilution of complex biological mixtures. *Analytical Chemistry*. 2006;78:4281–4290.
3. Dunn WB et al. Procedures for large-scale metabolic profiling of serum and plasma using GC- and LC-MS. *Nature Protocols*. 2011;6:1060–1083.
4. Xia J, Wishart DS. Web-based inference of biological patterns, functions and pathways from metabolomic data. *Nature Protocols*. 2011;6:743–760.
5. Xia J, Wishart DS. MetPA: a web-based metabolomics tool for pathway analysis and visualization. *Bioinformatics*. 2010;26:2342–2344.
6. Trygg J et al. Chemometrics in metabonomics. *Journal of Proteome Research*. 2007;6:469–479.

## 15. 可复现性

流程保存输入文件大小/修改时间、关键配置、步骤版本、归一化因子、QC 指标、PLS-DA 验证、差异设计和各步骤完成 manifest。随机步骤使用 `random_seed`。任何输入或关键参数变化都会使相关 manifest 失效，并触发必要步骤重跑。
