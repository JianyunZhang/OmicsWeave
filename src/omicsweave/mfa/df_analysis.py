#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 差异反应通路 DataFrame 分析模块（解耦版）

从 MFAAnalyser 中拆出的独立数据分析函数。
每个函数接收原始数据（DataFrame / dict），执行统计计算，
将结果 CSV 输出到 01_df_results/ 对应子目录，并返回 DataFrame。

目录结构：
    01_df_results/
      all_tissues/
        pca/              PCA 得分与载荷
        qc/               拟合质量
        overview/         组织通量概览
        mid/              MID 分析结果
        tca_cycle/        跨组织 TCA 汇总
      <tissue>/
        differential_flux/  差异通量分析结果
        pathway_flux/       通路通量样本数据
        tca_cycle/          TCA 循环通量
        heatmap/            热图原始数据与统计摘要

Usage:
    from df_analysis import differential_flux, flux_pca, ...
"""

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
import warnings
import zlib

from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from omicsweave.common.stability import calibrate_biological_stability
from omicsweave.common.nonparametric import (
    min_achievable_bh_fdr,
    min_achievable_pvalue,
)
from omicsweave.common.reaction_pathway_mapping import annotate_flux_dataframe

warnings.filterwarnings("ignore")

_REQUIRED_REACTION_ANNOTATION_COLUMNS = [
    "kegg_reaction",
    "kegg_pathway",
    "reactome_id",
    "ec_number",
    "description",
]


# ── 工具函数 ─────────────────────────────────────────────────────────────

def _parse_flux_column(col_name: str):
    """Parse 'tissue--reaction' or 'MIX:tissue--metabolite_idx' into components."""
    if col_name.startswith("MIX:"):
        rest = col_name[4:]
        parts = rest.split("--", 1)
        if len(parts) == 2:
            return "MIX", parts[0], parts[1]
        return "MIX", "unknown", rest
    parts = col_name.split("--", 1)
    if len(parts) == 2:
        return "flux", parts[0], parts[1]
    return "flux", "unknown", col_name


def _filter_iqr(values):
    """
    Filter values to keep exactly the middle 50% by index.

    Parameters
    ----------
    values : array-like
        Input values to filter.

    Returns
    -------
    np.ndarray
        Middle 50% of values by index (sorted).
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    n = len(values)
    if n == 0:
        return values

    sorted_values = np.sort(values)
    start = n // 4
    end = n - n // 4
    return sorted_values[start:end]


def _get_all_fits_for_samples(flux_all_fits, sample_ids, col, iqr_filter=True):
    """
    收集多个样本的所有拟合值，对每个样本先进行IQR过滤，然后合并。

    Parameters
    ----------
    flux_all_fits : dict
        {sample_id: DataFrame}，每个样本的所有拟合结果。
    sample_ids : list
        样本 ID 列表。
    col : str
        列名（flux_id）。

    Returns
    -------
    np.ndarray
        合并后的 IQR 过滤值。
    """
    all_vals = []
    for sid in sample_ids:
        if flux_all_fits and sid in flux_all_fits:
            fit_df = flux_all_fits[sid]
            if col in fit_df.columns:
                vals = pd.to_numeric(fit_df[col], errors='coerce').dropna().values
                if len(vals) > 0:
                    filtered = _filter_iqr(vals) if iqr_filter else vals
                    all_vals.extend(filtered)
    return np.array(all_vals) if all_vals else np.array([])


def _paired_net_fits(fit_df, forward_col, reverse_col):
    """单个样本逐次拟合配对的净通量：同一行（同一次拟合）的 forward − reverse。

    正、逆向通量只有在同一次拟合内才构成一个净通量；任一侧为 NaN 的
    拟合行整行丢弃。``fit_df`` 缺少逆向列时退回正向值。
    """
    forward = pd.to_numeric(fit_df[forward_col], errors="coerce")
    if reverse_col in fit_df.columns:
        reverse = pd.to_numeric(fit_df[reverse_col], errors="coerce")
        values = (forward - reverse).to_numpy(dtype=float)
    else:
        values = forward.to_numpy(dtype=float)
    return values[~np.isnan(values)]


def _compute_net_flux(
    flux_all_fits,
    flux_df,
    sample_ids,
    tissue,
    reaction,
    *,
    use_net_flux=True,
    iqr_filter=True,
):
    """
    计算可逆反应的净通量：forward - reverse (__R suffix)。

    净通量按拟合迭代逐行配对：先在每个样本内对同一次拟合求
    forward − reverse，再对净通量做 IQR 过滤，最后按样本拼接（每只动物
    各自过滤）。不能先分别过滤正、逆向再相减——``_filter_iqr`` 会排序，
    相减就成了「第 k 小正向 − 第 k 小逆向」，配对的是不同拟合迭代，
    结果可与真实净通量相差数百单位甚至反号。

    Parameters
    ----------
    flux_all_fits : dict
        所有拟合数据。
    flux_df : DataFrame
        最优拟合通量矩阵（用于检查列是否存在）。
    sample_ids : list
        样本 ID 列表。
    tissue : str
        组织名。
    reaction : str
        基础反应名（不含 __R 后缀）。

    Returns
    -------
    np.ndarray
        净通量值。
    """
    forward_col = f"{tissue}--{reaction}"
    reverse_col = f"{tissue}--{reaction}__R"

    if not use_net_flux or reverse_col not in flux_df.columns:
        return _get_all_fits_for_samples(
            flux_all_fits, sample_ids, forward_col, iqr_filter=iqr_filter
        )

    all_vals = []
    for sid in sample_ids:
        if flux_all_fits and sid in flux_all_fits:
            fit_df = flux_all_fits[sid]
            if forward_col in fit_df.columns:
                vals = _paired_net_fits(fit_df, forward_col, reverse_col)
                if len(vals) > 0:
                    filtered = _filter_iqr(vals) if iqr_filter else vals
                    all_vals.extend(filtered)
    return np.array(all_vals) if all_vals else np.array([])


def _compute_biological_net_flux(
    flux_all_fits,
    flux_df,
    sample_ids,
    tissue,
    reaction,
    *,
    use_net_flux=True,
    iqr_filter=True,
):
    """每个生物学样本先汇总为一个净通量估计，避免拟合点伪重复。"""
    sample_values = []
    for sample_id in sample_ids:
        values = _compute_net_flux(
            flux_all_fits,
            flux_df,
            [sample_id],
            tissue,
            reaction,
            use_net_flux=use_net_flux,
            iqr_filter=iqr_filter,
        )
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if len(finite):
            sample_values.append(float(np.median(finite)))
            continue
        # 无拟合点时退回最优拟合行；净通量同样取该行自身的 forward − reverse。
        column = f"{tissue}--{reaction}"
        if sample_id in flux_df.index and column in flux_df.columns:
            best_row = flux_df.loc[[sample_id]]
            reverse_col = f"{column}__R"
            if use_net_flux and reverse_col in flux_df.columns:
                fallback = _paired_net_fits(best_row, column, reverse_col)
            else:
                fallback = pd.to_numeric(best_row[column], errors="coerce").dropna().to_numpy()
            if len(fallback):
                sample_values.append(float(fallback[0]))
    return np.asarray(sample_values, dtype=float)


def _section(config: dict | None, name: str) -> dict:
    return (config or {}).get(name, {}) or {}


def _effective_statistic(
    statistic: str,
    low_n: bool,
    rank_test_underpowered: bool = False,
) -> str:
    """决定本次比较实际使用的统计路径。

    ``auto`` 模式下有两种情况必须放弃 p 值路径：

    1. ``low_n`` —— 任一组生物学重复 < 2，组内方差无法估计；
    2. ``rank_test_underpowered`` —— 秩检验的最小可达双侧 p 值
       （``min_achievable_pvalue``）已经高于配置阈值。此时"无显著结果"
       是组合数学的上限而非生物学结论，继续跑 p 值只会产生误导性的
       "0 个显著通量"。

    显式配置了具体检验方法时尊重配置，由调用方负责打印警告。
    """
    statistic = str(statistic or "auto").lower()
    if statistic == "auto":
        if low_n or rank_test_underpowered:
            return "effect_size_only"
        return "mannwhitney"
    return statistic


def _primary_pvalue_column(effective_statistic: str, pvalue_column: str | None = None) -> str | None:
    if effective_statistic == "effect_size_only":
        return None
    if effective_statistic in {"ttest", "t_test", "welch", "welch_ttest"}:
        return "pvalue_ttest"
    if effective_statistic in {"mannwhitney", "mann_whitney", "mannwhitneyu"}:
        return "pvalue_mannwhitney"
    return pvalue_column or "pvalue_mannwhitney"


def _log2fc_from_means(mean_a: float, mean_b: float, diff: float, zero_threshold: float) -> float:
    a_sig = abs(mean_a) > zero_threshold
    b_sig = abs(mean_b) > zero_threshold

    if a_sig and b_sig:
        abs_ratio = abs(mean_a) / abs(mean_b)
        if (mean_a > 0 and mean_b > 0) or (mean_a < 0 and mean_b < 0):
            log2fc = np.log2(abs_ratio)
            if mean_a < 0:
                log2fc = -log2fc
        else:
            log2fc = np.log2(abs_ratio)
            if diff < 0:
                log2fc = -log2fc
    elif a_sig and not b_sig:
        log2fc = np.inf if mean_a > 0 else -np.inf
    elif not a_sig and b_sig:
        log2fc = -np.inf if mean_b > 0 else np.inf
    else:
        log2fc = np.nan
    return log2fc


def _log2fc_standard_error(
    mean_a: float,
    mean_b: float,
    std_a: float,
    std_b: float,
    n_a: int,
    n_b: int,
    zero_threshold: float,
) -> float:
    """用 delta method 估计 ``log2(|mean_a|/|mean_b|)`` 的标准误。"""
    values = (mean_a, mean_b, std_a, std_b)
    if n_a < 2 or n_b < 2 or not all(np.isfinite(value) for value in values):
        return np.nan
    if abs(mean_a) <= zero_threshold or abs(mean_b) <= zero_threshold:
        return np.nan
    ln2 = np.log(2.0)
    variance = (
        (std_a ** 2) / (n_a * mean_a ** 2 * ln2 ** 2)
        + (std_b ** 2) / (n_b * mean_b ** 2 * ln2 ** 2)
    )
    return float(np.sqrt(variance)) if variance >= 0 else np.nan


def _standardized_mean_difference(
    mean_a: float,
    mean_b: float,
    std_a: float,
    std_b: float,
    n_a: int,
    n_b: int,
) -> tuple[float, float]:
    """返回基于生物学重复的 Cohen's d 与小样本校正 Hedges' g。"""
    degrees_freedom = n_a + n_b - 2
    if n_a < 2 or n_b < 2 or degrees_freedom <= 0:
        return np.nan, np.nan
    pooled_variance = (
        (n_a - 1) * std_a ** 2 + (n_b - 1) * std_b ** 2
    ) / degrees_freedom
    if not np.isfinite(pooled_variance) or pooled_variance <= 0:
        return np.nan, np.nan
    cohen_d = float((mean_a - mean_b) / np.sqrt(pooled_variance))
    correction = 1.0 - 3.0 / (4.0 * degrees_freedom - 1.0)
    return cohen_d, float(correction * cohen_d)


def _ensure_dir(output_dir: Path, *subdirs):
    """确保子目录存在，返回子目录 Path。"""
    d = output_dir.joinpath(*subdirs)
    d.mkdir(parents=True, exist_ok=True)
    return d


# 视为"未被数据决定"的可辨识性标签。require_identifiable 与跨 protocol
# 表决（mfa_integration._UNCONSTRAINED_LABELS，需保持同步）都按此集合剔除。
UNCONSTRAINED_LABELS = ("poorly_constrained", "bound_pinned")

# 盒约束贴边判定的默认值（可由 analysis.identifiability 覆盖）。
# 本项目拟合软件的方向通量盒约束为 [1, 3000]：下界 1.0 处的点质量在
# data_20251204 中占全部拟合点的 15%。净通量 = 正向 − 逆向，两向同贴下界
# 时恰为 0，因此 0 与 1.0 是最常见的贴边值。
_DEFAULT_BOUND_VALUES = (0.0, 1.0)
_DEFAULT_BOUND_TOLERANCE = 0.05
_DEFAULT_BOUND_PINNED_MIN_FRACTION = 0.5


def _infer_fit_box_bounds(flux_all_fits) -> tuple[float, float] | None:
    """从全部拟合点的**方向**通量推断盒约束上下界。

    导出的工作簿不含通量上下界，但盒约束在拟合点分布上留有明确痕迹：
    优化器把大量迭代推到边界上，形成最小值 / 最大值处的点质量。
    这里只取全局极值；是否"贴边"由后续的容差判定决定——某只动物的
    代表净通量（拟合中位数）落在极值附近，本身就意味着至少一半拟合点
    挤在那里，偶然的尾部极值不会触发。
    """
    lows, highs = [], []
    for fit_df in (flux_all_fits or {}).values():
        cols = [c for c in fit_df.columns
                if "--" in str(c) and not str(c).startswith("MIX:")]
        if not cols:
            continue
        values = fit_df[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if len(values):
            lows.append(values.min())
            highs.append(values.max())
    if not lows:
        return None
    return float(min(lows)), float(max(highs))


def _net_bound_values(bound_values, inferred_bounds) -> list[float]:
    """汇总"净通量贴边"的候选绝对值。

    方向通量的盒约束为 [lb, ub] 时，净通量 f − r 贴边的情形有：
    不可逆反应 f = lb / ub；可逆反应两向同贴一侧得 0，分贴两侧得 ub − lb。
    """
    values = {abs(float(v)) for v in (bound_values or [])}
    if inferred_bounds is not None:
        lb, ub = inferred_bounds
        values.update({abs(lb), abs(ub), abs(ub - lb)})
    # 推断出的下界常为 0.9999998 这类浮点值，与配置的 1.0 是同一个边界。
    return sorted({round(v, 6) for v in values})


def _flag_at_bound(values, bound_values, tolerance) -> np.ndarray:
    """|value| 落在任一贴边值 ± tolerance 内即为 True。"""
    values = np.abs(np.asarray(values, dtype=float))
    bounds = np.asarray(list(bound_values), dtype=float)
    if values.size == 0 or bounds.size == 0:
        return np.zeros(values.shape, dtype=bool)
    near = np.abs(values[:, None] - bounds[None, :]) <= float(tolerance)
    return near.any(axis=1) & np.isfinite(values)


def _load_identifiability(output_dir) -> pd.DataFrame | None:
    """读取 sample_flux_matrix 写出的每通量可辨识性表。

    该表由同一次运行的 ``sample_flux_matrix()`` 生成（流程中先于差异分析），
    因此差异结果可以直接带上"这条通量到底被数据约束得有多紧"的标注。
    """
    path = (Path(output_dir) / "01_df_results" / "all_tissues"
            / "flux_matrix" / "flux_identifiability.csv")
    if not path.exists():
        return None
    try:
        table = pd.read_csv(path)
    except Exception:
        return None
    columns = [c for c in ("flux_id", "median_rel_iqr", "frac_at_bound",
                           "identifiability")
               if c in table.columns]
    if "flux_id" not in columns:
        return None
    return table[columns].drop_duplicates(subset=["flux_id"])


def _enforce_standard_reaction_annotations(
    df: pd.DataFrame,
    *,
    reaction_col: str = "reaction",
    strict: bool = True,
) -> pd.DataFrame:
    """
    为差异通量结果强制补齐标准数据库注释列。

    strict=True 时：
    - reaction 列缺失直接报错；
    - 注释列未成功写入直接报错。
    """
    if reaction_col not in df.columns:
        if strict:
            raise ValueError(f"Missing required reaction column '{reaction_col}' for annotation.")
        return df

    df = annotate_flux_dataframe(df, reaction_col=reaction_col)

    missing_cols = [c for c in _REQUIRED_REACTION_ANNOTATION_COLUMNS if c not in df.columns]
    if missing_cols and strict:
        raise RuntimeError(
            "Failed to attach mandatory reaction annotation columns: "
            + ", ".join(missing_cols)
        )

    for col in missing_cols:
        df[col] = np.nan

    return df


# ── QC ────────────────────────────────────────────────────────────────────

def quality_control(flux_all_fits, sample_conditions, output_dir):
    """
    评估拟合质量：每个样本的 Loss 分布。

    Parameters
    ----------
    flux_all_fits : dict
        {sample_id: DataFrame}。
    sample_conditions : dict
        {sample_id: condition}。
    output_dir : str or Path
        MFA 输出根目录。

    Returns
    -------
    pd.DataFrame or None
    """
    print("\n--- Quality Control ---")
    rows = []
    if flux_all_fits:
        for sid, df in flux_all_fits.items():
            cond = sample_conditions.get(sid, "unknown")
            loss_col = "loss" if "loss" in df.columns else "Loss"
            if loss_col in df.columns:
                losses = pd.to_numeric(df[loss_col], errors="coerce").dropna()
                rows.append({
                    "sample": sid,
                    "condition": cond,
                    "min_loss": losses.min(),
                    "max_loss": losses.max(),
                    "mean_loss": losses.mean(),
                    "n_fits": len(losses),
                })

    if not rows:
        print("  No loss data available for QC")
        return None

    qc_df = pd.DataFrame(rows)
    save_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "qc")
    qc_df.to_csv(save_dir / "qc_fitting_loss.csv", index=False)
    print("  Saved all_tissues/qc/qc_fitting_loss.csv")
    return qc_df


# ── PCA ───────────────────────────────────────────────────────────────────

def flux_pca(flux_df, sample_conditions, output_dir, tissues,
             n_components=None):
    """
    对通量数据进行 PCA 分析，按组织分别计算。

    Parameters
    ----------
    flux_df : pd.DataFrame
        最优拟合通量矩阵 (samples × fluxes)。
    sample_conditions : dict
        {sample_id: condition}。
    output_dir : str or Path
        MFA 输出根目录。
    tissues : list of str
        组织列表。
    n_components : int, optional
        PCA 降维维度。

    Returns
    -------
    pd.DataFrame or None
    """
    if flux_df is None:
        print("[SKIP] No flux data loaded")
        return None

    print("\n--- Flux PCA ---")
    all_scores = []

    for tissue in tissues:
        tissue_cols = [c for c in flux_df.columns
                       if _parse_flux_column(c)[1] == tissue
                       and _parse_flux_column(c)[0] == "flux"]

        if len(tissue_cols) < 2:
            print(f"  [SKIP] {tissue}: not enough flux columns ({len(tissue_cols)})")
            continue

        tissue_df = flux_df[tissue_cols].dropna(axis=0, how='all')
        if tissue_df.shape[0] < 2:
            print(f"  [SKIP] {tissue}: not enough samples ({tissue_df.shape[0]})")
            continue

        X = tissue_df.values
        n_max = min(X.shape[0], X.shape[1], 10)
        n_comp = n_components or min(5, n_max)
        if n_comp < 2:
            n_comp = 2

        scaler = StandardScaler()
        Xs = scaler.fit_transform(X)

        pca = PCA(n_components=n_comp)
        scores = pca.fit_transform(Xs)
        explained = pca.explained_variance_ratio_

        conditions = [sample_conditions.get(s, "unknown") for s in tissue_df.index]

        score_df = pd.DataFrame(scores, index=tissue_df.index,
                                columns=[f"PC{i+1}" for i in range(len(explained))])
        score_df["condition"] = conditions
        score_df["tissue"] = tissue
        all_scores.append(score_df)

        loadings = pd.DataFrame(
            pca.components_.T,
            index=tissue_cols,
            columns=[f"PC{i+1}" for i in range(len(explained))],
        )
        # 将方差解释率追加为最后一行
        loadings.loc["variance_explained"] = explained

        pca_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "pca")
        loadings.to_csv(pca_dir / f"flux_pca_loadings_{tissue}.csv")

        print(f"  {tissue}: variance explained = {', '.join(f'PC{i+1}={v:.1%}' for i, v in enumerate(explained))}")

    if all_scores:
        combined_scores = pd.concat(all_scores)
        pca_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "pca")
        combined_scores.to_csv(pca_dir / "flux_pca_scores.csv")
        return combined_scores
    return None


# ── 差异通量分析 ──────────────────────────────────────────────────────────

def differential_flux(flux_df, flux_all_fits, sample_conditions, conditions,
                      tissues, output_dir, comparisons=None, tissue_filter=None,
                      analysis_config=None):
    """
    每个生物学样本先汇总为一个通量估计，再进行差异通量分析。

    Parameters
    ----------
    flux_df : pd.DataFrame
        最优拟合通量矩阵。
    flux_all_fits : dict
        所有拟合数据。
    sample_conditions : dict
        {sample_id: condition}。
    conditions : list of str
        条件列表。
    tissues : list of str
        组织列表（用于信息提示）。
    output_dir : str or Path
        MFA 输出根目录。
    comparisons : list of tuple, optional
        比较对列表，如 [("hfd", "lfd")]。
    tissue_filter : list of str, optional
        限制特定组织。

    Returns
    -------
    dict of {str: pd.DataFrame}
    """
    if flux_df is None:
        print("[SKIP] No flux data loaded")
        return {}

    if comparisons is None:
        comparisons = list(combinations(conditions, 2))

    preprocessing_cfg = _section(analysis_config, "preprocessing")
    diff_cfg = _section(analysis_config, "differential_flux")
    use_net_flux = bool(preprocessing_cfg.get("use_net_flux", True))
    iqr_filter = bool(preprocessing_cfg.get("iqr_filter", True))
    zero_threshold = float(preprocessing_cfg.get("zero_threshold", 1.0e-10))
    statistic = str(diff_cfg.get("statistic", "auto")).lower()
    multiple_testing = str(diff_cfg.get("multiple_testing", "fdr_bh")).lower()
    pvalue_threshold = float(diff_cfg.get("pvalue_threshold", 0.05))
    fdr_threshold = float(diff_cfg.get("fdr_threshold", 0.10))
    log2fc_threshold = float(diff_cfg.get("log2fc_threshold", 0.5))
    cohen_d_threshold = float(diff_cfg.get("cohen_d_threshold", 0.5))
    # 默认只标注不过滤：可辨识性差是模型/数据问题，交由使用者判断是否剔除。
    # 开启后 UNCONSTRAINED_LABELS（poorly_constrained、bound_pinned）强制不显著。
    require_identifiable = bool(diff_cfg.get("require_identifiable", False))
    hedges_g_threshold = float(diff_cfg.get("hedges_g_threshold", cohen_d_threshold))
    stability_cfg = dict(diff_cfg.get("stability", {}) or {})
    stability_enabled = bool(stability_cfg.get("enabled", False))
    stability_use_for_selection = bool(stability_cfg.get("use_for_selection", False))
    stability_alpha = float(stability_cfg.get("alpha", 0.05))
    min_selection_probability = float(
        stability_cfg.get("min_selection_probability", 0.80)
    )
    min_direction_stability = float(
        stability_cfg.get("min_direction_stability", 0.80)
    )
    annotation_cfg = _section(analysis_config, "annotation")
    require_standard_annotations = bool(
        annotation_cfg.get("require_standard_reaction_ids", True)
    )

    filter_label = "IQR Filter" if iqr_filter else "No IQR Filter"
    flux_label = "Net Flux" if use_net_flux else "Forward Flux"
    print(f"\n--- Differential Flux Analysis ({flux_label}, {filter_label}) ---")

    flux_cols = [c for c in flux_df.columns if not c.startswith("MIX:")]
    results = {}

    for cond_a, cond_b in comparisons:
        comp_name = f"{cond_a.upper()}_vs_{cond_b.upper()}"

        sids_a = [s for s, c in sample_conditions.items()
                  if c == cond_a and s in flux_df.index]
        sids_b = [s for s, c in sample_conditions.items()
                  if c == cond_b and s in flux_df.index]

        n_bio_a, n_bio_b = len(sids_a), len(sids_b)

        if n_bio_a < 1 or n_bio_b < 1:
            print(f"  Comparison: {comp_name}")
            print(f"    [SKIP] No samples: {cond_a}={n_bio_a}, {cond_b}={n_bio_b}")
            continue

        low_n = n_bio_a < 2 or n_bio_b < 2

        # 秩检验的最小可达双侧 p 值只由样本量决定（组合数学上限）。
        # 若它已经高于配置阈值，p 值路径永远不可能显著——必须显式提示，
        # 否则"0 个显著通量"会被误读为生物学结论。
        p_floor = min_achievable_pvalue(n_bio_a, n_bio_b)
        rank_test_underpowered = bool(
            np.isfinite(p_floor) and p_floor > pvalue_threshold
        )

        print(f"  Comparison: {comp_name}")
        print(f"    n(bio): {cond_a}={n_bio_a}, {cond_b}={n_bio_b}"
              f" | Mann-Whitney min two-sided p = {p_floor:.4g}")
        if low_n:
            print(f"    [LOW-N] 任一组 < 2 只动物，组内方差无法估计")
        if rank_test_underpowered:
            print(f"    [POWER] p_floor {p_floor:.4g} > pvalue_threshold {pvalue_threshold:g}"
                  f" —— p 值路径在数学上不可能达到显著")

        rows = []
        processed_reactions = set()

        for col in flux_cols:
            kind, tissue, reaction = _parse_flux_column(col)
            if tissue_filter and tissue not in tissue_filter:
                continue

            if reaction.endswith("__R"):
                continue

            reaction_key = f"{tissue}--{reaction}"
            if reaction_key in processed_reactions:
                continue
            processed_reactions.add(reaction_key)

            vals_a = _compute_biological_net_flux(
                flux_all_fits, flux_df, sids_a, tissue, reaction,
                use_net_flux=use_net_flux, iqr_filter=iqr_filter,
            )
            vals_b = _compute_biological_net_flux(
                flux_all_fits, flux_df, sids_b, tissue, reaction,
                use_net_flux=use_net_flux, iqr_filter=iqr_filter,
            )

            if len(vals_a) == 0 or len(vals_b) == 0:
                continue

            n_obs_a, n_obs_b = len(vals_a), len(vals_b)
            mean_a, mean_b = np.nanmean(vals_a), np.nanmean(vals_b)
            std_a = np.nanstd(vals_a, ddof=1) if n_obs_a >= 2 else np.nan
            std_b = np.nanstd(vals_b, ddof=1) if n_obs_b >= 2 else np.nan
            diff = mean_a - mean_b

            log2fc = _log2fc_from_means(mean_a, mean_b, diff, zero_threshold)
            log2fc_se = _log2fc_standard_error(
                mean_a, mean_b, std_a, std_b,
                n_obs_a, n_obs_b, zero_threshold,
            )

            # 兼容保留 Cohen's d；证据判定优先使用小样本校正 Hedges' g。
            cohen_d, hedges_g = _standardized_mean_difference(
                mean_a, mean_b, std_a, std_b, n_obs_a, n_obs_b
            )

            if stability_enabled:
                base_seed = int(stability_cfg.get("random_state", 0) or 0)
                reaction_seed = zlib.crc32(
                    f"{comp_name}|{reaction_key}".encode("utf-8")
                )
                stability = calibrate_biological_stability(
                    vals_a,
                    vals_b,
                    n_bootstrap=int(stability_cfg.get("n_bootstrap", 2000)),
                    n_permutations=int(stability_cfg.get("n_permutations", 5000)),
                    alpha=stability_alpha,
                    random_state=(base_seed + reaction_seed) % (2 ** 32),
                )
                stability_values = stability.to_dict()
            else:
                stability_values = {
                    "calibrated_effect_threshold": np.nan,
                    "selection_probability": np.nan,
                    "direction_stability": np.nan,
                    "permutation_pvalue": np.nan,
                    "evidence_mode": "not_requested",
                }

            try:
                u_stat, pval_mw = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
            except Exception:
                u_stat, pval_mw = np.nan, np.nan

            try:
                t_stat, pval = stats.ttest_ind(vals_a, vals_b, equal_var=False)
            except Exception:
                t_stat, pval = np.nan, np.nan

            rows.append({
                "flux_id": col,
                "tissue": tissue,
                "reaction": reaction,
                f"n_bio_{cond_a}": n_bio_a,
                f"n_bio_{cond_b}": n_bio_b,
                f"n_points_{cond_a}": len(vals_a),
                f"n_points_{cond_b}": len(vals_b),
                "inference_unit": "biological_replicate",
                "n_eff": n_bio_a + n_bio_b,
                "diff_se": (
                    np.sqrt(std_a ** 2 / n_obs_a + std_b ** 2 / n_obs_b)
                    if n_obs_a >= 2 and n_obs_b >= 2 else np.nan
                ),
                f"mean_{cond_a}": mean_a,
                f"mean_{cond_b}": mean_b,
                f"std_{cond_a}": std_a,
                f"std_{cond_b}": std_b,
                "diff": diff,
                "log2FC": log2fc,
                "log2fc_se": log2fc_se,
                "cohen_d": cohen_d,
                "hedges_g": hedges_g,
                "effect_size_metric": "hedges_g",
                "calibrated_effect_threshold": stability_values["calibrated_effect_threshold"],
                "selection_probability": stability_values["selection_probability"],
                "direction_stability": stability_values["direction_stability"],
                "permutation_pvalue": stability_values["permutation_pvalue"],
                "stability_evidence_mode": stability_values["evidence_mode"],
                "t_statistic": t_stat,
                "pvalue_ttest": pval,
                "pvalue_mannwhitney": pval_mw,
            })

        if not rows:
            continue

        df = pd.DataFrame(rows)

        # 样本量下界（与数据无关，同一比较内为常数）
        df["min_achievable_pvalue"] = p_floor
        df["rank_test_underpowered"] = rank_test_underpowered
        if np.isfinite(p_floor):
            n_at_floor = int(df["pvalue_mannwhitney"].le(p_floor + 1e-12).sum())
            df["min_achievable_fdr"] = min_achievable_bh_fdr(
                p_floor, len(df), max(1, n_at_floor)
            )
        else:
            df["min_achievable_fdr"] = np.nan

        if multiple_testing == "fdr_bh":
            from statsmodels.stats.multitest import multipletests
            mask = df["pvalue_ttest"].notna()
            if mask.sum() > 0:
                _, fdr, _, _ = multipletests(df.loc[mask, "pvalue_ttest"], method="fdr_bh")
                df.loc[mask, "fdr_ttest"] = fdr

            mask = df["pvalue_mannwhitney"].notna()
            if mask.sum() > 0:
                _, fdr, _, _ = multipletests(df.loc[mask, "pvalue_mannwhitney"], method="fdr_bh")
                df.loc[mask, "fdr_mannwhitney"] = fdr

        # BH 校正后的可达下界只有在检验总数已知后才能算。选择规则是
        # "p < 阈值 **且** FDR <= 阈值"，所以 FDR 不可达时整条 p 值
        # 路径同样无效，auto 模式必须一并回退到效应量。
        fdr_floor = float(df["min_achievable_fdr"].iloc[0])
        fdr_unreachable = bool(
            multiple_testing == "fdr_bh"
            and np.isfinite(fdr_floor)
            and fdr_floor > fdr_threshold
        )
        if fdr_unreachable and not rank_test_underpowered:
            print(f"    [POWER] 未校正 p 可达，但 BH 校正后最小可达 FDR ="
                  f" {fdr_floor:.3g} > fdr_threshold {fdr_threshold:g}"
                  f"（{len(df)} 个检验）—— 选择规则整体不可达")
        rank_test_underpowered = bool(rank_test_underpowered or fdr_unreachable)
        df["fdr_unreachable"] = fdr_unreachable
        df["rank_test_underpowered"] = rank_test_underpowered

        effective_stat = _effective_statistic(statistic, low_n, rank_test_underpowered)
        if rank_test_underpowered and effective_stat != "effect_size_only":
            print(f"    [WARN] 配置显式指定 statistic={statistic}，在 n={n_bio_a}v{n_bio_b} 下"
                  f" 该检验无法达到 p<{pvalue_threshold:g}；结果仅供参考，"
                  f"建议改用 statistic: auto")
        if statistic == "auto" and rank_test_underpowered and not low_n:
            print(f"    [AUTO] 已自动回退到 effect_size_only（基于生物学重复的 Hedges' g）")

        # 证据等级：区分"能做统计检验"、"只能给效应量"、"连效应量都没有"。
        if low_n:
            evidence_level = "insufficient"
        elif effective_stat == "effect_size_only":
            evidence_level = "effect_size_only"
        else:
            evidence_level = "statistical_test"

        primary_p_col = _primary_pvalue_column(effective_stat, diff_cfg.get("pvalue_column"))
        primary_fdr_col = None
        if primary_p_col == "pvalue_ttest":
            primary_fdr_col = "fdr_ttest"
        elif primary_p_col == "pvalue_mannwhitney":
            primary_fdr_col = "fdr_mannwhitney"

        if effective_stat == "effect_size_only":
            df["is_significant"] = (
                df["hedges_g"].abs().ge(hedges_g_threshold)
                & df["log2FC"].abs().ge(log2fc_threshold)
            )
            df["abs_hedges_g"] = df["hedges_g"].abs()
            df = df.sort_values("abs_hedges_g", ascending=False)
            df = df.drop(columns=["abs_hedges_g"])
        else:
            if primary_p_col in df.columns:
                p_sig = df[primary_p_col].lt(pvalue_threshold)
            else:
                p_sig = pd.Series(False, index=df.index)
            if multiple_testing == "fdr_bh" and primary_fdr_col in df.columns:
                p_sig = p_sig & df[primary_fdr_col].le(fdr_threshold)
            df["is_significant"] = p_sig & df["log2FC"].abs().ge(log2fc_threshold)
            df = df.sort_values(primary_p_col)

        if stability_enabled and stability_use_for_selection:
            stability_pass = (
                df["selection_probability"].ge(min_selection_probability)
                & df["direction_stability"].ge(min_direction_stability)
                & df["permutation_pvalue"].le(stability_alpha)
            )
            df["is_significant"] = df["is_significant"] & stability_pass.fillna(False)

        df["primary_statistic"] = effective_stat
        df["primary_pvalue_column"] = primary_p_col or ""
        df["multiple_testing"] = multiple_testing
        df["pvalue_threshold"] = pvalue_threshold
        df["fdr_threshold"] = fdr_threshold
        df["log2fc_threshold"] = log2fc_threshold
        df["cohen_d_threshold"] = cohen_d_threshold
        df["hedges_g_threshold"] = hedges_g_threshold
        df["stability_enabled"] = stability_enabled
        df["stability_used_for_selection"] = stability_use_for_selection

        if low_n and effective_stat != "effect_size_only":
            stat_note = "low_n_configured_statistical_test"
        elif effective_stat == "effect_size_only":
            stat_note = "fitting_variance_only"
        else:
            stat_note = "biological_replicates"

        # 兼容旧下游：保留 stat_note，同时添加更明确的配置列。
        df["stat_note"] = stat_note
        df["evidence_level"] = evidence_level
        df["effect_size_forced_by_power_floor"] = bool(
            rank_test_underpowered
            and not low_n
            and effective_stat == "effect_size_only"
        )

        # 如果显著性列因为缺少主 P 值未能得到布尔值，回退为 False。
        if df["is_significant"].dtype != bool:
            df["is_significant"] = df["is_significant"].fillna(False).astype(bool)

        # n<2 时 Hedges' g 与 p 值均无定义，任何"显著"都是伪结论。
        if evidence_level == "insufficient":
            df["is_significant"] = False

        if effective_stat == "effect_size_only":
            summary_rule = (
                f"|Hedges' g|>={hedges_g_threshold} and |log2FC|>={log2fc_threshold}"
            )
        else:
            fdr_part = f" and {primary_fdr_col}<={fdr_threshold}" if multiple_testing == "fdr_bh" else ""
            summary_rule = (
                f"{primary_p_col}<{pvalue_threshold}{fdr_part} "
                f"and |log2FC|>={log2fc_threshold}"
            )
        df["significance_rule"] = summary_rule

        # 如果没有统计阈值显著项，仍保留完整表格，图表阶段按配置阈值展示。
        if effective_stat == "effect_size_only":
            df["abs_hedges_g"] = df["hedges_g"].abs()
            df = df.sort_values("abs_hedges_g", ascending=False)
            df = df.drop(columns=["abs_hedges_g"])

        # 可辨识性标注：拟合迭代的相对 IQR 宽度反映该通量被数据约束的程度。
        # 未被约束的通量即便"显著"也不值得解释，这里让它在结果表里可见。
        identifiability = _load_identifiability(output_dir)
        if identifiability is not None:
            df = df.merge(identifiability, on="flux_id", how="left")
            if require_identifiable:
                unconstrained = df["identifiability"].isin(UNCONSTRAINED_LABELS)
                df.loc[unconstrained, "is_significant"] = False
        else:
            df["median_rel_iqr"] = np.nan
            df["frac_at_bound"] = np.nan
            df["identifiability"] = "undetermined"

        # 强制附加标准数据库注释（用于多组学整合稳定映射）
        df = _enforce_standard_reaction_annotations(
            df,
            reaction_col="reaction",
            strict=require_standard_annotations,
        )
        mapped_mask = df[["kegg_reaction", "kegg_pathway", "reactome_id", "ec_number"]].notna().any(axis=1)
        print(f"    annotation coverage: {int(mapped_mask.sum())}/{len(df)} reactions with DB IDs")

        # 按组织分组保存到 <tissue>/differential_flux/ 子目录
        for tissue_name, group_df in df.groupby("tissue", sort=True):
            save_dir = _ensure_dir(output_dir, "01_df_results", str(tissue_name), "differential_flux")
            group_df.to_csv(save_dir / f"differential_flux_{comp_name}.csv", index=False)

        results[comp_name] = df

        sig = df[df["is_significant"]]
        if "identifiability" in df.columns and len(sig):
            poor = int(sig["identifiability"].eq("poorly_constrained").sum())
            if poor:
                print(f"    [IDENT] {poor}/{len(sig)} 条入选通量属 poorly_constrained"
                      f"（拟合 IQR 宽度已达通量量级），解释时需谨慎")
            pinned = int(sig["identifiability"].eq("bound_pinned").sum())
            if pinned:
                print(f"    [IDENT] {pinned}/{len(sig)} 条入选通量属 bound_pinned"
                      f"（多数动物的拟合贴在盒约束边界），差异可能只是边界效应")
        if effective_stat == "effect_size_only":
            sig_fc = df[df["log2FC"].abs() >= log2fc_threshold]
            print(f"    {len(df)} fluxes analyzed ({effective_stat}, sorted by |Hedges' g|)")
            print(f"    {len(sig)} fluxes significant ({summary_rule})")
            print(f"    {len(sig_fc)} fluxes with |log2FC|>={log2fc_threshold}")
        else:
            print(f"    {len(sig)} / {len(df)} fluxes significant ({summary_rule})")
            fdr_floor = float(df["min_achievable_fdr"].iloc[0])
            if (multiple_testing == "fdr_bh" and np.isfinite(fdr_floor)
                    and fdr_floor > fdr_threshold):
                print(f"    [POWER] BH 校正后最小可达 FDR = {fdr_floor:.3g}"
                      f" > fdr_threshold {fdr_threshold:g}（{len(df)} 个检验）"
                      f" —— FDR 路径同样不可能显著")

    return results


# ── 通路通量比较 ──────────────────────────────────────────────────────────

def pathway_flux_comparison(flux_df, flux_all_fits, sample_conditions,
                            conditions, tissues, pathway_groups,
                            output_dir, comparisons=None,
                            analysis_config=None):
    """
    按代谢通路分组，收集所有 IQR 过滤后的拟合数据（净通量），并输出：
    1) 样本级拟合点明细（用于箱线图）
    2) 按 reaction 的统计表（TCA 风格：mean/std/n）
    3) 按 pathway 的汇总统计表

    Parameters
    ----------
    flux_df, flux_all_fits, sample_conditions, conditions, tissues, pathway_groups :
        同上。
    output_dir : str or Path
        MFA 输出根目录。
    comparisons : list of tuple, optional
        比较对列表。

    Returns
    -------
    dict[str, pd.DataFrame] or None
    """
    if flux_df is None:
        print("[SKIP] No flux data loaded")
        return None

    if comparisons is None:
        comparisons = list(combinations(conditions, 2))

    preprocessing_cfg = _section(analysis_config, "preprocessing")
    use_net_flux = bool(preprocessing_cfg.get("use_net_flux", True))
    iqr_filter = bool(preprocessing_cfg.get("iqr_filter", True))

    flux_label = "Net Flux" if use_net_flux else "Forward Flux"
    filter_label = "IQR Filtered" if iqr_filter else "Unfiltered"
    print(f"\n--- Pathway Flux Comparison ({flux_label}, {filter_label}) ---")

    all_fitting_rows = []
    for tissue in tissues:
        tissue_cols = {col: _parse_flux_column(col)
                       for col in flux_df.columns
                       if _parse_flux_column(col)[1] == tissue
                       and _parse_flux_column(col)[0] == "flux"}

        for pw_name, pw_reactions in pathway_groups.items():
            for col, (_, _, rxn) in tissue_cols.items():
                if rxn.endswith("__R"):
                    continue

                base_rxn = rxn
                if base_rxn not in pw_reactions:
                    continue

                for sample_id in flux_df.index:
                    cond = sample_conditions.get(sample_id, "unknown")

                    net_vals = _compute_net_flux(
                        flux_all_fits, flux_df, [sample_id], tissue, base_rxn,
                        use_net_flux=use_net_flux, iqr_filter=iqr_filter,
                    )

                    for val in net_vals:
                        all_fitting_rows.append({
                            "tissue": tissue,
                            "pathway": pw_name,
                            "reaction": base_rxn,
                            "flux_id": f"{tissue}--{base_rxn}",
                            "condition": cond,
                            "sample_id": sample_id,
                            "flux": float(val),
                        })

    if not all_fitting_rows:
        print("  No pathway fluxes found")
        return None

    sample_df = pd.DataFrame(all_fitting_rows)
    sample_df = sample_df.dropna(subset=["flux"])
    if sample_df.empty:
        print("  No pathway fluxes after filtering")
        return None

    reaction_stats_df = (
        sample_df.groupby(["tissue", "pathway", "reaction", "condition"], as_index=False)
        .agg(
            mean=("flux", "mean"),
            std=("flux", "std"),
            n=("sample_id", "nunique"),
            n_points=("flux", "size"),
        )
    )
    reaction_stats_df["std"] = reaction_stats_df["std"].fillna(0.0)

    pathway_stats_df = (
        sample_df.groupby(["tissue", "pathway", "condition"], as_index=False)
        .agg(
            mean=("flux", "mean"),
            std=("flux", "std"),
            n=("sample_id", "nunique"),
            n_points=("flux", "size"),
            n_reactions=("reaction", "nunique"),
        )
    )
    pathway_stats_df["std"] = pathway_stats_df["std"].fillna(0.0)

    for tissue_name, group_df in sample_df.groupby("tissue", sort=True):
        save_dir = _ensure_dir(output_dir, "01_df_results", str(tissue_name), "pathway_flux")
        group_df.to_csv(save_dir / "pathway_flux_samples.csv", index=False)
        reaction_stats_df[reaction_stats_df["tissue"] == tissue_name].to_csv(
            save_dir / "pathway_flux_reaction_stats.csv",
            index=False,
        )
        pathway_stats_df[pathway_stats_df["tissue"] == tissue_name].to_csv(
            save_dir / "pathway_flux_summary_stats.csv",
            index=False,
        )

    # 跨组织总表：便于下游统一读取与审计
    all_tissues_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "pathway_flux")
    reaction_stats_df.to_csv(all_tissues_dir / "pathway_flux_reaction_stats.csv", index=False)
    pathway_stats_df.to_csv(all_tissues_dir / "pathway_flux_summary_stats.csv", index=False)

    print(
        f"  Saved pathway flux data for {len(tissues)} tissues "
        f"(samples + reaction stats + pathway summary, based on {len(pathway_groups)} pathway groups)"
    )
    return {
        "samples": sample_df,
        "reaction_stats": reaction_stats_df,
        "pathway_summary": pathway_stats_df,
    }


# ── 热图统计 ─────────────────────────────────────────────────────────────

def flux_heatmap(flux_df, flux_all_fits, sample_conditions, conditions,
                 tissues, output_dir, comparisons=None, top_n=50,
                 p_thresh=0.05, analysis_config=None):
    """
    差异通量热图数据计算（原始值 + 统计摘要）。

    Parameters
    ----------
    (同前)

    Returns
    -------
    dict of {str: pd.DataFrame}
    """
    if flux_df is None:
        print("[SKIP] No flux data loaded")
        return {}

    if comparisons is None:
        comparisons = []
        if "lfd" in conditions:
            for cond in conditions:
                if cond != "lfd":
                    comparisons.append((cond, "lfd"))
        else:
            comparisons = list(combinations(conditions, 2))

    if not comparisons:
        print("[SKIP] No valid comparisons for heatmap")
        return {}

    preprocessing_cfg = _section(analysis_config, "preprocessing")
    heatmap_cfg = _section(analysis_config, "heatmap")
    use_net_flux = bool(preprocessing_cfg.get("use_net_flux", True))
    iqr_filter = bool(preprocessing_cfg.get("iqr_filter", True))
    zero_threshold = float(preprocessing_cfg.get("zero_threshold", 1.0e-10))
    scale_method = str(preprocessing_cfg.get("scale_method", "zscore")).lower()
    missing_value = str(preprocessing_cfg.get("missing_value", "drop")).lower()
    rank_by = str(heatmap_cfg.get("rank_by", "auto")).lower()
    p_thresh = float(heatmap_cfg.get("pvalue_threshold", p_thresh))
    top_n = int(heatmap_cfg.get("top_n", top_n))

    print("\n--- Flux Heatmap (Differential, Individual Samples) ---")

    all_results = {}

    for cond_a, cond_b in comparisons:
        comp_name = f"{cond_a.upper()}_vs_{cond_b.upper()}"
        print(f"  Generating heatmap data for: {comp_name}")

        sids_a = [s for s, c in sample_conditions.items()
                  if c == cond_a and s in flux_df.index]
        sids_b = [s for s, c in sample_conditions.items()
                  if c == cond_b and s in flux_df.index]

        n_bio_a, n_bio_b = len(sids_a), len(sids_b)

        if n_bio_a < 1 or n_bio_b < 1:
            print(f"    [SKIP] Not enough samples")
            continue

        low_n = n_bio_a < 2 or n_bio_b < 2

        # 与 differential_flux 保持同一套推断单位与样本量下界判定：
        # 反应筛选统计必须基于生物学重复，绝不能用拟合迭代点冒充重复。
        p_floor = min_achievable_pvalue(n_bio_a, n_bio_b)
        rank_test_underpowered = bool(
            np.isfinite(p_floor) and p_floor > p_thresh
        )
        effect_size_mode = bool(low_n or rank_test_underpowered)
        if effect_size_mode:
            print(f"    [POWER] n={n_bio_a}v{n_bio_b}, min p = {p_floor:.4g} > {p_thresh:g}"
                  f" —— 热图按 |Hedges g| 排序，不施加 p 值阈值")

        all_sids = sids_a + sids_b

        for tissue in tissues:
            tissue_cols = [c for c in flux_df.columns
                           if _parse_flux_column(c)[1] == tissue
                           and _parse_flux_column(c)[0] == "flux"
                           and not c.startswith("MIX:")]

            if len(tissue_cols) < 2:
                continue

            sig_fluxes = []
            processed = set()

            for col in tissue_cols:
                kind, _, reaction = _parse_flux_column(col)

                if reaction.endswith("__R"):
                    continue

                reaction_key = f"{tissue}--{reaction}"
                if reaction_key in processed:
                    continue
                processed.add(reaction_key)

                vals_a = _compute_biological_net_flux(
                    flux_all_fits, flux_df, sids_a, tissue, reaction,
                    use_net_flux=use_net_flux, iqr_filter=iqr_filter,
                )
                vals_b = _compute_biological_net_flux(
                    flux_all_fits, flux_df, sids_b, tissue, reaction,
                    use_net_flux=use_net_flux, iqr_filter=iqr_filter,
                )

                if len(vals_a) == 0 or len(vals_b) == 0:
                    continue

                mean_a = np.nanmean(vals_a)
                mean_b = np.nanmean(vals_b)

                diff = mean_a - mean_b
                log2fc = _log2fc_from_means(mean_a, mean_b, diff, zero_threshold)

                n_a, n_b = len(vals_a), len(vals_b)
                std_a = np.nanstd(vals_a, ddof=1) if n_a >= 2 else np.nan
                std_b = np.nanstd(vals_b, ddof=1) if n_b >= 2 else np.nan
                cohen_d, hedges_g = _standardized_mean_difference(
                    mean_a, mean_b, std_a, std_b, n_a, n_b
                )

                try:
                    _, pval = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
                except Exception:
                    pval = np.nan

                # p 值路径可达时用 P-value 筛选；不可达时保留全部反应，
                # 后续按 |Hedges g| 排序截取 top-N。
                if not effect_size_mode:
                    is_sig = pval < p_thresh
                else:
                    is_sig = True

                if is_sig:
                    if pd.notna(log2fc) and not np.isinf(log2fc):
                        sig_fluxes.append({
                            "flux_id": col,
                            "reaction": reaction,
                            "log2FC": log2fc,
                            "pvalue": pval,
                            "cohen_d": cohen_d,
                            "hedges_g": hedges_g,
                            f"n_bio_{cond_a}": n_a,
                            f"n_bio_{cond_b}": n_b,
                        })

            if not sig_fluxes:
                criterion = (f"top-{top_n} by |Hedges g|" if effect_size_mode
                             else f"p<{p_thresh}")
                print(f"    {tissue}: no significant fluxes ({criterion})")
                continue

            sig_df = pd.DataFrame(sig_fluxes)

            # 排序策略可由 config 控制；auto 保持旧逻辑。
            if rank_by == "abs_cohen_d" or (rank_by == "auto" and effect_size_mode):
                # Hedges g 是小样本校正后的效应量，优先于未校正的 Cohen d
                sig_df["_sort_key"] = sig_df["hedges_g"].abs().fillna(
                    sig_df["cohen_d"].abs()
                )
            elif rank_by == "pvalue":
                sig_df["_sort_key"] = sig_df["pvalue"]
                sig_df = sig_df.sort_values("_sort_key", ascending=True).head(top_n)
                sig_df = sig_df.drop(columns=["_sort_key"])
            else:
                sig_df["_sort_key"] = sig_df["log2FC"].abs()
            if "_sort_key" in sig_df.columns:
                sig_df = sig_df.sort_values("_sort_key", ascending=False).head(top_n)
                sig_df = sig_df.drop(columns=["_sort_key"])

            if len(sig_df) < 2:
                print(f"    {tissue}: not enough significant fluxes ({len(sig_df)})")
                continue

            selected_reactions = sig_df["reaction"].tolist()

            heatmap_rows = []
            for reaction in selected_reactions:
                row = {"reaction": reaction}
                for sid in all_sids:
                    vals = _compute_biological_net_flux(
                        flux_all_fits, flux_df, [sid], tissue, reaction,
                        use_net_flux=use_net_flux, iqr_filter=iqr_filter,
                    )
                    row[sid] = float(vals[0]) if len(vals) > 0 else np.nan
                heatmap_rows.append(row)

            heatmap_df = pd.DataFrame(heatmap_rows).set_index("reaction")
            heatmap_df = heatmap_df[all_sids]

            if missing_value == "median_fill":
                heatmap_df = heatmap_df.apply(
                    lambda row: row.fillna(row.median()),
                    axis=1,
                )
            else:
                heatmap_df = heatmap_df.dropna(axis=0, how="any")

            if heatmap_df.shape[0] < 2:
                print(f"    {tissue}: not enough heatmap rows after missing-value handling ({heatmap_df.shape[0]})")
                continue

            if scale_method == "zscore":
                heatmap_df = heatmap_df.apply(
                    lambda row: (row - row.mean()) / row.std() if row.std() > 0 else row * 0,
                    axis=1
                )

            mode_tag = (" (effect-size mode, ranked by |Hedges g|)"
                        if effect_size_mode else "")
            print(f"    {tissue}: {len(sig_df)} significant fluxes, {len(all_sids)} individual samples{mode_tag}")

            save_dir = _ensure_dir(output_dir, "01_df_results", tissue, "heatmap")
            heatmap_df.to_csv(
                save_dir / f"flux_heatmap_raw_{tissue}_{comp_name}.csv"
            )
            sig_df.to_csv(
                save_dir / f"flux_heatmap_stats_{tissue}_{comp_name}.csv",
                index=False
            )

            # 保存比较元数据（用于下游绘图判断模式）
            meta_path = save_dir / f"heatmap_meta_{comp_name}.csv"
            meta_df = pd.DataFrame([{
                "comparison": comp_name,
                f"n_bio_{cond_a}": n_bio_a,
                f"n_bio_{cond_b}": n_bio_b,
                "low_n": low_n,
                "effect_size_mode": effect_size_mode,
                "rank_test_underpowered": rank_test_underpowered,
                "min_achievable_pvalue": p_floor,
                "inference_unit": "biological_replicate",
                "stat_note": ("fitting_variance_only" if effect_size_mode
                              else "biological_replicates"),
                "top_n": top_n,
                "pvalue_threshold": p_thresh,
                "rank_by": rank_by,
                "scale_method": scale_method,
                "missing_value": missing_value,
                "use_net_flux": use_net_flux,
                "iqr_filter": iqr_filter,
            }])
            meta_df.to_csv(meta_path, index=False)

            all_results[f"{tissue}_{comp_name}"] = {
                "raw": heatmap_df,
                "stats": sig_df,
            }

    return all_results


# ── TCA 循环 ─────────────────────────────────────────────────────────────

def tca_cycle_analysis(flux_df, sample_conditions, conditions, tissues,
                       pathway_groups, output_dir):
    """
    TCA 循环各反应通量的组织间比较。

    Parameters
    ----------
    (同前)

    Returns
    -------
    pd.DataFrame or None
    """
    if flux_df is None:
        print("[SKIP] No flux data loaded")
        return None

    print("\n--- TCA Cycle Analysis ---")

    tca_rxns = pathway_groups.get("TCA_Cycle", [])
    rows = []

    for tissue in tissues:
        for rxn in tca_rxns:
            col = f"{tissue}--{rxn}"
            if col not in flux_df.columns:
                continue
            for cond in conditions:
                sids = [s for s, c in sample_conditions.items()
                        if c == cond and s in flux_df.index]
                vals = flux_df.loc[sids, col].values.astype(float)
                rows.append({
                    "tissue": tissue,
                    "reaction": rxn,
                    "condition": cond,
                    "mean": np.nanmean(vals),
                    "std": np.nanstd(vals, ddof=1) if len(vals) > 1 else 0,
                    "n": len(vals),
                })

    if not rows:
        print("  No TCA cycle fluxes found")
        return None

    tca_df = pd.DataFrame(rows)

    # 为避免与 pathway_flux_reaction_stats.csv 信息冗余，
    # 不再单独输出 01_df_results/*/tca_cycle/ 目录。
    print("  TCA cycle summary computed in-memory (standalone tca_cycle CSV output disabled)")
    return tca_df


# ── MID 分析 ─────────────────────────────────────────────────────────────

def mid_analysis(mid_df, sample_conditions, conditions, output_dir,
                 top_metabolites=20):
    """
    质量同位素分布（MID）比较分析。

    Parameters
    ----------
    mid_df : pd.DataFrame
        MID 数据矩阵。
    sample_conditions : dict
        {sample_id: condition}。
    conditions : list of str
        条件列表。
    output_dir : str or Path
        MFA 输出根目录。
    top_metabolites : int
        Top 代谢物数量。

    Returns
    -------
    pd.DataFrame or None
    """
    if mid_df is None:
        print("[SKIP] No MID data loaded")
        return None

    print("\n--- MID Analysis ---")

    all_rows = []
    for (metab, midx), row in mid_df.iterrows():
        for cond in conditions:
            sids = [s for s, c in sample_conditions.items()
                    if c == cond and s in mid_df.columns]
            if not sids:
                continue
            vals = row[sids].values.astype(float)
            all_rows.append({
                "metabolite": metab,
                "mid_index": midx,
                "condition": cond,
                "mean": np.nanmean(vals),
                "std": np.nanstd(vals, ddof=1) if len(vals) > 1 else 0,
            })

    mid_summary = pd.DataFrame(all_rows)
    save_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "mid")
    mid_summary.to_csv(save_dir / "mid_summary.csv", index=False)

    total_label = mid_summary[mid_summary["mid_index"] > 0].groupby(
        ["metabolite", "condition"]
    )["mean"].sum().reset_index()
    total_label.columns = ["metabolite", "condition", "total_labeling"]
    total_label.to_csv(save_dir / "mid_total_labeling.csv", index=False)

    var_per_metab = total_label.groupby("metabolite")["total_labeling"].var()
    top_metabs = var_per_metab.nlargest(top_metabolites).index.tolist()

    print(f"  Saved MID summary and labeling data ({len(top_metabs)} top metabolites)")
    return mid_summary


# ── 组织通量概览 ─────────────────────────────────────────────────────────

def tissue_flux_overview(flux_df, sample_conditions, conditions, tissues,
                         output_dir):
    """
    各组织总通量幅度汇总。

    Parameters
    ----------
    (同前)

    Returns
    -------
    pd.DataFrame or None
    """
    if flux_df is None:
        print("[SKIP] No flux data loaded")
        return None

    print("\n--- Tissue Flux Overview ---")

    rows = []
    for tissue in tissues:
        tissue_cols = [c for c in flux_df.columns
                       if _parse_flux_column(c)[1] == tissue
                       and _parse_flux_column(c)[0] == "flux"
                       and not _parse_flux_column(c)[2].endswith("__R")]

        for cond in conditions:
            sids = [s for s, c in sample_conditions.items()
                    if c == cond and s in flux_df.index]
            if not sids:
                continue
            vals = flux_df.loc[sids, tissue_cols].values.astype(float)
            total = np.nanmean(np.abs(vals).sum(axis=1))
            n_rxn = len(tissue_cols)
            rows.append({
                "tissue": tissue,
                "condition": cond,
                "n_reactions": n_rxn,
                "total_abs_flux": total,
                "mean_abs_flux": total / n_rxn if n_rxn > 0 else 0,
            })

    if not rows:
        return None

    ov_df = pd.DataFrame(rows)
    save_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "overview")
    ov_df.to_csv(save_dir / "tissue_flux_overview.csv", index=False)
    print("  Saved tissue flux overview")
    return ov_df


# -- 样本级通量矩阵 --------------------------------------------------------

def sample_flux_matrix(flux_df, flux_all_fits, sample_conditions, tissues,
                       output_dir, analysis_config=None):
    """输出**每只动物一行一反应**的净通量代表值。

    这是整个 MFA 下游唯一保留样本粒度的表，用途有三：

    1. 跨批次桥接样本一致性检查（同一只动物在两个批次/两台仪器上被拟合
       两次时，二者之差就是可测量的技术噪声地板）；
    2. 差异分析与热图的推断单位审计（确认统计确实建立在生物学重复上）；
    3. 后续通量置信区间 / 可辨识性工作的输入。

    同时记录每个样本在该反应上的拟合点离散度（IQR 与标准差），作为
    该通量在这只动物身上是否被数据良好约束的粗略代理；以及代表值是否
    贴在盒约束边界上（``at_bound``）——贴边时离散度为 0 却并非被数据
    决定，汇总到 ``flux_identifiability.csv`` 时以 ``bound_pinned`` 优先标注。

    Returns
    -------
    pd.DataFrame or None
        长表：sample_id x flux_id。
    """
    if flux_df is None:
        print("[SKIP] No flux data loaded")
        return None

    preprocessing_cfg = _section(analysis_config, "preprocessing")
    use_net_flux = bool(preprocessing_cfg.get("use_net_flux", True))
    iqr_filter = bool(preprocessing_cfg.get("iqr_filter", True))
    ident_cfg = _section(analysis_config, "identifiability")
    bound_values = ident_cfg.get("bound_values", list(_DEFAULT_BOUND_VALUES))
    bound_tolerance = float(ident_cfg.get("bound_tolerance", _DEFAULT_BOUND_TOLERANCE))
    infer_bounds = bool(ident_cfg.get("infer_bounds_from_fits", True))
    pinned_min_fraction = float(
        ident_cfg.get("bound_pinned_min_fraction", _DEFAULT_BOUND_PINNED_MIN_FRACTION)
    )

    print("\n--- Sample-level Flux Matrix ---")

    flux_cols = [c for c in flux_df.columns if not c.startswith("MIX:")]
    seen = set()
    reaction_keys = []
    for col in flux_cols:
        kind, tissue, reaction = _parse_flux_column(col)
        if kind != "flux" or reaction.endswith("__R"):
            continue
        key = (tissue, reaction)
        if key in seen:
            continue
        seen.add(key)
        reaction_keys.append(key)

    rows = []
    for sample_id in flux_df.index:
        condition = sample_conditions.get(sample_id, "unknown")
        for tissue, reaction in reaction_keys:
            net_value = _compute_biological_net_flux(
                flux_all_fits, flux_df, [sample_id], tissue, reaction,
                use_net_flux=use_net_flux, iqr_filter=iqr_filter,
            )
            if len(net_value) == 0:
                continue

            fit_points = _compute_net_flux(
                flux_all_fits, flux_df, [sample_id], tissue, reaction,
                use_net_flux=use_net_flux, iqr_filter=iqr_filter,
            )
            fit_points = np.asarray(fit_points, dtype=float)
            fit_points = fit_points[np.isfinite(fit_points)]

            value = float(net_value[0])
            p25 = float(np.percentile(fit_points, 25)) if len(fit_points) else np.nan
            p75 = float(np.percentile(fit_points, 75)) if len(fit_points) else np.nan
            # 可辨识性代理：拟合点的 IQR 宽度相对于通量本身的大小。
            # 拟合迭代不是生物学重复，但它们**确实**刻画了「在这只动物的
            # 数据下，这条通量被约束得有多紧」。相对宽度很大意味着该通量
            # 几乎不被数据决定，其组间差异不值得解释。
            scale = abs(value)
            rel_iqr = ((p75 - p25) / scale
                       if np.isfinite(p25) and np.isfinite(p75) and scale > 0
                       else np.nan)
            rows.append({
                "sample_id": sample_id,
                "condition": condition,
                "tissue": tissue,
                "reaction": reaction,
                "flux_id": f"{tissue}--{reaction}",
                "net_flux": value,
                "n_fit_points": int(len(fit_points)),
                "fit_p25": p25,
                "fit_p75": p75,
                "fit_std": float(np.std(fit_points, ddof=1)) if len(fit_points) >= 2 else np.nan,
                "fit_rel_iqr": rel_iqr,
                "use_net_flux": use_net_flux,
                "iqr_filter": iqr_filter,
            })

    if not rows:
        print("  No sample-level fluxes to write")
        return None

    matrix_df = pd.DataFrame(rows)

    # ── 盒约束贴边 ────────────────────────────────────────────────────
    # 拟合点全部贴在边界上时 IQR = 0，相对 IQR 也是 0，若只看离散度会被
    # 判为"约束良好"；实际上是数据推不动这条通量，优化器把它压在了边界上。
    inferred_bounds = _infer_fit_box_bounds(flux_all_fits) if infer_bounds else None
    net_bounds = _net_bound_values(bound_values, inferred_bounds)
    matrix_df["at_bound"] = _flag_at_bound(
        matrix_df["net_flux"].to_numpy(), net_bounds, bound_tolerance
    )
    bound_source = (f"inferred box [{inferred_bounds[0]:.6g}, {inferred_bounds[1]:.6g}] + config"
                    if inferred_bounds is not None else "config only")
    print(f"  Bound check: |net_flux| within {bound_tolerance:g} of "
          f"{[float(f'{b:.6g}') for b in net_bounds]} ({bound_source})")

    save_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "flux_matrix")
    matrix_df.to_csv(save_dir / "sample_net_flux.csv", index=False)
    print(f"  Saved sample_net_flux.csv: "
          f"{matrix_df['sample_id'].nunique()} samples x "
          f"{matrix_df['flux_id'].nunique()} fluxes "
          f"({len(matrix_df)} rows)")

    # ── 按通量汇总可辨识性 ────────────────────────────────────────────
    identifiability = (
        matrix_df.groupby(["tissue", "reaction", "flux_id"], as_index=False)
        .agg(n_samples=("sample_id", "nunique"),
             median_rel_iqr=("fit_rel_iqr", "median"),
             max_rel_iqr=("fit_rel_iqr", "max"),
             median_abs_flux=("net_flux", lambda s: float(np.median(np.abs(s)))),
             n_fit_points=("n_fit_points", "median"),
             frac_at_bound=("at_bound", "mean"))
    )
    # 阈值取相对 IQR：<0.2 视为约束良好，<1.0 尚可，>=1.0 表示 IQR 宽度
    # 已达通量本身量级，等同于未被数据决定。
    iqr_class = pd.cut(
        identifiability["median_rel_iqr"],
        bins=[-np.inf, 0.2, 1.0, np.inf],
        labels=["well_constrained", "moderate", "poorly_constrained"],
    ).astype(object)
    iqr_class[identifiability["median_rel_iqr"].isna()] = "undetermined"
    # 仅凭 IQR 的分类单独留档，便于核对贴边判定改掉了哪些标签。
    identifiability["iqr_identifiability"] = iqr_class
    # 贴边优先于 IQR 分档：多数动物的代表值压在盒约束上，零离散度只说明
    # 拟合点全挤在边界，不说明数据决定了它。
    pinned = identifiability["frac_at_bound"] >= pinned_min_fraction
    identifiability["identifiability"] = iqr_class.where(~pinned, "bound_pinned")
    identifiability = identifiability.sort_values(
        "median_rel_iqr", ascending=False, na_position="last"
    )
    identifiability.to_csv(save_dir / "flux_identifiability.csv", index=False)

    counts = identifiability["identifiability"].value_counts().to_dict()
    print(f"  Flux identifiability (median relative IQR of model fits, "
          f"bound_pinned first): {counts}")
    n_total = len(identifiability)
    n_poor = int(counts.get("poorly_constrained", 0))
    if n_poor:
        share = n_poor / n_total
        print(f"  [!] {n_poor}/{n_total} ({share:.0%}) 条通量的拟合 IQR"
              f"宽度已达通量本身量级——这些通量基本未被数据约束，")
        print(f"      其组间差异不应作为结论。见 flux_identifiability.csv")
    n_pinned = int(counts.get("bound_pinned", 0))
    if n_pinned:
        share = n_pinned / n_total
        was_well = int(
            identifiability.loc[pinned, "iqr_identifiability"]
            .eq("well_constrained").sum()
        )
        print(f"  [!] {n_pinned}/{n_total} ({share:.0%}) 条通量在 ≥{pinned_min_fraction:.0%}"
              f" 的动物中贴在盒约束边界上——拟合点被压在边界而非由数据决定，")
        print(f"      其中 {was_well} 条若只看拟合 IQR 会被误判为 well_constrained。"
              f"见 flux_identifiability.csv 的 frac_at_bound")
    return matrix_df


# -- MID 拟合残差 QC -------------------------------------------------------

# MID 是分数（各同位素体占比，合计约 1），因此残差阈值直接用绝对分数表示。
# 典型 GC/LC-MS 的 MID 测量误差在 0.003-0.02 之间，故 0.05（5 mole%）以上的
# 均方根残差已属明显失拟。
_MID_RMSE_GOOD = 0.02
_MID_RMSE_ACCEPTABLE = 0.05


def mid_fit_residuals(mid_df, mid_measured_df, sample_conditions, output_dir):
    """实测 MID 与模型模拟 MID 的残差 QC。

    ``mid_raw_data.xlsx`` 的 ``experiments`` 列是**实测**质量同位素分布，
    此前被读入后直接丢弃。没有它，"模型到底拟合得好不好"这个 13C-MFA 最
    基本的问题就无从回答——而 loss 只是一个汇总标量，看不出是哪个代谢物、
    哪个同位素体拟合失败。

    本函数逐 (样本 x 代谢物) 计算残差，并按代谢物汇总，定位失拟来源。

    Returns
    -------
    pd.DataFrame or None
        每个 (样本, 代谢物) 一行的残差汇总表。
    """
    if mid_df is None or mid_measured_df is None:
        print("[SKIP] No MID data for residual QC")
        return None

    print("\n--- MID Fit Residuals (measured vs simulated) ---")

    shared_index = mid_df.index.intersection(mid_measured_df.index)
    shared_columns = [c for c in mid_df.columns if c in mid_measured_df.columns]
    if len(shared_index) == 0 or not shared_columns:
        print("  [SKIP] 实测与模拟 MID 没有共同的代谢物 / 样本")
        return None

    simulated = mid_df.loc[shared_index, shared_columns]
    measured = mid_measured_df.loc[shared_index, shared_columns]

    long_rows = []
    for (metabolite, mid_index), sim_row in simulated.iterrows():
        meas_row = measured.loc[(metabolite, mid_index)]
        for sample_id in shared_columns:
            sim_value = pd.to_numeric(pd.Series([sim_row[sample_id]]), errors="coerce").iloc[0]
            meas_value = pd.to_numeric(pd.Series([meas_row[sample_id]]), errors="coerce").iloc[0]
            if pd.isna(sim_value) and pd.isna(meas_value):
                continue
            long_rows.append({
                "sample_id": sample_id,
                "condition": sample_conditions.get(sample_id, "unknown"),
                "metabolite": metabolite,
                "mid_index": mid_index,
                "measured": meas_value,
                "simulated": sim_value,
                "residual": (sim_value - meas_value)
                if pd.notna(sim_value) and pd.notna(meas_value) else np.nan,
            })

    if not long_rows:
        print("  [SKIP] 无可比较的 MID 数据点")
        return None

    long_df = pd.DataFrame(long_rows)
    save_dir = _ensure_dir(output_dir, "01_df_results", "all_tissues", "mid")
    long_df.to_csv(save_dir / "mid_residuals_long.csv", index=False)

    # ── 逐 (样本 x 代谢物) 汇总 ───────────────────────────────────────────
    def _summarise(group: pd.DataFrame) -> pd.Series:
        residuals = group["residual"].dropna()
        return pd.Series({
            "n_mid": int(len(residuals)),
            "ssr": float(np.sum(residuals ** 2)) if len(residuals) else np.nan,
            "rmse": float(np.sqrt(np.mean(residuals ** 2))) if len(residuals) else np.nan,
            "max_abs_residual": float(np.max(np.abs(residuals))) if len(residuals) else np.nan,
            "measured_sum": float(group["measured"].sum(skipna=True)),
            "simulated_sum": float(group["simulated"].sum(skipna=True)),
        })

    per_sample = (
        long_df.groupby(["sample_id", "condition", "metabolite"])
        .apply(_summarise, include_groups=False)
        .reset_index()
    )
    per_sample["fit_quality"] = np.select(
        [per_sample["rmse"] < _MID_RMSE_GOOD,
         per_sample["rmse"] < _MID_RMSE_ACCEPTABLE],
        ["good", "acceptable"],
        default="poor",
    )
    # 实测 MID 应当归一到 1；明显偏离说明上游归一或自然丰度校正有问题。
    per_sample["measured_sum_off_by"] = (per_sample["measured_sum"] - 1.0).abs()
    per_sample.to_csv(save_dir / "mid_fit_residuals.csv", index=False)

    # ── 按代谢物汇总 ─────────────────────────────────────────────────────
    by_metabolite = (
        per_sample.groupby("metabolite", as_index=False)
        .agg(n_samples=("sample_id", "nunique"),
             median_rmse=("rmse", "median"),
             max_rmse=("rmse", "max"),
             median_max_abs_residual=("max_abs_residual", "median"),
             n_poor=("fit_quality", lambda s: int((s == "poor").sum())))
    )
    by_metabolite["frac_poor"] = (
        by_metabolite["n_poor"] / by_metabolite["n_samples"].replace(0, np.nan)
    )
    by_metabolite = by_metabolite.sort_values("median_rmse", ascending=False)
    by_metabolite.to_csv(save_dir / "mid_fit_quality_by_metabolite.csv", index=False)

    quality_counts = per_sample["fit_quality"].value_counts().to_dict()
    n_total = len(per_sample)
    print(f"  {n_total} 个 (样本 x 代谢物) 组合: {quality_counts}")
    print(f"  阈值: RMSE<{_MID_RMSE_GOOD:g} good, <{_MID_RMSE_ACCEPTABLE:g} acceptable")

    n_poor = int(quality_counts.get("poor", 0))
    if n_poor:
        share = n_poor / n_total if n_total else 0.0
        print(f"  [!] {n_poor}/{n_total} ({share:.0%}) 组合失拟 (RMSE>={_MID_RMSE_ACCEPTABLE:g})")
        worst = by_metabolite.head(8)
        print("      失拟最严重的代谢物:")
        for _, row in worst.iterrows():
            print(f"        {row['metabolite']}: median RMSE={row['median_rmse']:.4f}, "
                  f"{int(row['n_poor'])}/{int(row['n_samples'])} 样本失拟")
        print("      模型无法复现实测 MID 时，由它导出的通量不具解释力——")
        print("      应优先检查网络拓扑、稳态假设与自然丰度校正，而非直接比较通量。")

    off_normalized = per_sample[per_sample["measured_sum_off_by"] > 0.05]
    if len(off_normalized):
        print(f"  [!] {len(off_normalized)}/{n_total} 组合的实测 MID 合计偏离 1 超过 0.05，"
              f"提示上游归一或自然丰度校正存在问题")

    print("  Saved all_tissues/mid/mid_fit_residuals.csv, mid_fit_quality_by_metabolite.csv")
    return per_sample
