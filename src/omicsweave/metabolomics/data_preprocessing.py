"""
代谢组学数据预处理模块：加载、分组、过滤、插补、标准化。
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PreprocessingOptions:
    """代谢物预处理参数。比例参数均使用 0-1 小数。"""

    normalization_method: str = "pqn"
    detection_rate: float = 0.5
    imputation_fraction: float = 0.5
    qc_filter_enabled: bool = True
    qc_rsd_cutoff: float = 0.30
    qc_min_detection_rate: float = 0.80


@dataclass
class PreprocessingResult:
    """预处理公开结果；供流程、报告和测试共同消费。"""

    processed_df: pd.DataFrame
    sample_cols: list[str]
    qc_cols: list[str]
    qc_metrics: pd.DataFrame
    normalization_factors: pd.Series
    sample_totals: pd.DataFrame
    summary: dict[str, Any]


def load_tissue_data(path: str | Path) -> pd.DataFrame:
    """加载单个组织的 Excel 数据文件。"""
    print(f"  加载数据: {path}")
    df = pd.read_excel(path)
    print(f"    Shape: {df.shape}")
    return df


def _validate_options(options: PreprocessingOptions) -> None:
    if not 0 < options.detection_rate <= 1:
        raise ValueError("detection_rate 必须在 (0, 1] 范围内")
    if not 0 < options.imputation_fraction <= 1:
        raise ValueError("imputation_fraction 必须在 (0, 1] 范围内")
    if options.qc_rsd_cutoff <= 0:
        raise ValueError("qc_rsd_cutoff 必须大于 0")
    if not 0 < options.qc_min_detection_rate <= 1:
        raise ValueError("qc_min_detection_rate 必须在 (0, 1] 范围内")


def _normalize_abundances(
    data: pd.DataFrame,
    method: str,
) -> tuple[pd.DataFrame, pd.Series]:
    """按样本归一化，返回归一化矩阵和除数型 scaling factor。"""
    normalized_method = str(method).strip().lower().replace("-", "_")
    aliases = {
        "probabilistic_quotient": "pqn",
        "median_centering": "median",
        "total_ion_current": "tic",
        "off": "none",
    }
    normalized_method = aliases.get(normalized_method, normalized_method)

    if normalized_method == "pqn":
        reference = data.median(axis=1, skipna=True)
        reference = reference.where(reference > 0)
        quotients = data.div(reference, axis=0)
        factors = quotients.replace([np.inf, -np.inf], np.nan).median(axis=0)
    elif normalized_method == "median":
        sample_medians = data.median(axis=0, skipna=True)
        target = float(sample_medians.median())
        factors = sample_medians / target
    elif normalized_method == "tic":
        sample_totals = data.sum(axis=0, skipna=True)
        target = float(sample_totals.median())
        factors = sample_totals / target
    elif normalized_method == "none":
        factors = pd.Series(1.0, index=data.columns, dtype=float)
    else:
        raise ValueError(
            "不支持的 normalization_method: "
            f"{method!r}；可选 pqn/median/tic/none"
        )

    factors = pd.to_numeric(factors, errors="coerce").astype(float)
    invalid = ~np.isfinite(factors) | (factors <= 0)
    if invalid.any():
        bad = ", ".join(map(str, factors.index[invalid]))
        raise ValueError(f"无法为以下样本计算有效归一化因子: {bad}")
    return data.div(factors, axis=1), factors.rename("normalization_factor")


def _compute_qc_metrics(data: pd.DataFrame, qc_cols: list[str]) -> pd.DataFrame:
    if not qc_cols:
        return pd.DataFrame(
            index=data.index,
            columns=["QC_Mean", "QC_SD", "QC_RSD", "QC_DetectionRate", "QC_N"],
            dtype=float,
        )
    qc = data[qc_cols]
    means = qc.mean(axis=1, skipna=True)
    counts = qc.notna().sum(axis=1)
    std = qc.std(axis=1, skipna=True, ddof=1)
    rsd = std / means.replace(0, np.nan)
    return pd.DataFrame(
        {
            "QC_Mean": means,
            "QC_SD": std,
            "QC_RSD": rsd,
            "QC_DetectionRate": counts / len(qc_cols),
            "QC_N": counts,
        },
        index=data.index,
    )


def preprocess_metabolomics_data(
    df: pd.DataFrame,
    sample_cols: list[str],
    groups: dict[str, str] | None = None,
    *,
    qc_cols: list[str] | None = None,
    options: PreprocessingOptions | None = None,
) -> PreprocessingResult:
    """
    过滤、QC、逐代谢物插补、样本归一化并 log2 转换峰面积表。

    生物样本的检出过滤按组计算：代谢物只需在任一组达到阈值即可保留，
    从而保留组特异的 on/off 信号。若未提供组别，则回退到全样本检出率。
    """
    opts = options or PreprocessingOptions()
    _validate_options(opts)
    qc_cols = list(qc_cols or [])
    sample_cols = list(sample_cols)
    intensity_cols = list(dict.fromkeys(sample_cols + qc_cols))

    missing_cols = [col for col in intensity_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"数据中缺少强度列: {missing_cols}")
    if not sample_cols:
        raise ValueError("未识别到生物学样本列")
    if groups is not None:
        missing_groups = [col for col in sample_cols if col not in groups]
        if missing_groups:
            raise ValueError(f"以下样本缺少组别映射: {missing_groups}")

    data = df[intensity_cols].apply(pd.to_numeric, errors="coerce")
    data = data.mask(data <= 0)
    original_feature_count = len(data)

    if groups:
        keep_by_group = pd.Series(False, index=data.index)
        for group_name in dict.fromkeys(groups[col] for col in sample_cols):
            group_cols = [col for col in sample_cols if groups[col] == group_name]
            keep_by_group |= data[group_cols].notna().mean(axis=1) >= opts.detection_rate
        detection_keep = keep_by_group
        detection_rule = "at_least_one_group"
    else:
        detection_keep = data[sample_cols].notna().mean(axis=1) >= opts.detection_rate
        detection_rule = "pooled_samples"

    data = data.loc[detection_keep].copy()
    removed_detection = int((~detection_keep).sum())

    qc_metrics = _compute_qc_metrics(data, qc_cols)
    removed_qc = 0
    qc_filter_applied = bool(opts.qc_filter_enabled and len(qc_cols) >= 2)
    if qc_filter_applied:
        qc_keep = (
            (qc_metrics["QC_DetectionRate"] >= opts.qc_min_detection_rate)
            & np.isfinite(qc_metrics["QC_RSD"])
            & (qc_metrics["QC_RSD"] <= opts.qc_rsd_cutoff)
        )
        qc_metrics["QC_Pass"] = qc_keep.astype(bool)
        removed_qc = int((~qc_keep).sum())
        data = data.loc[qc_keep].copy()
    else:
        qc_metrics["QC_Pass"] = pd.Series(pd.NA, index=qc_metrics.index, dtype="boolean")

    missing_before = int(data.isna().sum().sum())
    row_min = data.min(axis=1, skipna=True)
    impute_values = row_min * opts.imputation_fraction
    data = data.T.fillna(impute_values).T
    if data.isna().any().any():
        raise ValueError("存在无法插补的全缺失代谢物；请检查检出过滤参数")

    totals_before = data.sum(axis=0)
    normalized, factors = _normalize_abundances(data, opts.normalization_method)
    totals_after = normalized.sum(axis=0)
    transformed = np.log2(normalized)

    processed_df = df.loc[transformed.index].copy()
    processed_df[intensity_cols] = transformed
    sample_totals = pd.DataFrame(
        {
            "sample": intensity_cols,
            "sample_type": ["biological" if c in sample_cols else "QC" for c in intensity_cols],
            "total_before": totals_before.reindex(intensity_cols).to_numpy(),
            "normalization_factor": factors.reindex(intensity_cols).to_numpy(),
            "total_after": totals_after.reindex(intensity_cols).to_numpy(),
        }
    ).set_index("sample")

    method = str(opts.normalization_method).strip().lower().replace("-", "_")
    method = {
        "probabilistic_quotient": "pqn",
        "median_centering": "median",
        "total_ion_current": "tic",
        "off": "none",
    }.get(method, method)
    summary: dict[str, Any] = {
        "original_feature_count": original_feature_count,
        "filtered_feature_count": len(processed_df),
        "removed_by_detection": removed_detection,
        "removed_by_qc": removed_qc,
        "imputed_cell_count": missing_before,
        "detection_rule": detection_rule,
        "detection_rate": opts.detection_rate,
        "normalization_method": method,
        "qc_sample_count": len(qc_cols),
        "qc_filter_applied": qc_filter_applied,
        "qc_rsd_cutoff": opts.qc_rsd_cutoff,
    }
    return PreprocessingResult(
        processed_df=processed_df,
        sample_cols=sample_cols,
        qc_cols=qc_cols,
        qc_metrics=qc_metrics,
        normalization_factors=factors,
        sample_totals=sample_totals,
        summary=summary,
    )


@dataclass
class DataColumnClassification:
    """峰表列分类及从样本名解析出的实验元数据。"""

    biological: list[str]
    unmapped_biological: list[str]
    qc: list[str]
    blank: list[str]
    metadata: list[str]
    mouse_ids: dict[str, str]
    tissues: dict[str, str]


# 样本列名: <alphanumeric_mouse_id>_<tissue> 或旧版连字符格式。
# 要求 mouse id 至少含一个数字，以免把 Class-I 等元数据误识别为样本。
_SAMPLE_PATTERN = re.compile(
    r"^((?=[A-Za-z0-9]*\d)[A-Za-z0-9]+)[_-](.+)$",
)
_QC_PATTERN = re.compile(r"^QC(?:[_-]?\d+|[_-].+)?$", re.IGNORECASE)
_BLANK_PATTERN = re.compile(r"^(?:blank|blk)(?:[_-]?\d+|[_-].+)?$", re.IGNORECASE)


def classify_data_columns(
    columns: pd.Index | list[object],
    known_mouse_ids: set[str] | list[str] | tuple[str, ...] | None = None,
) -> DataColumnClassification:
    """可靠地区分生物样本、QC、空白和元数据列。"""
    known = {str(value) for value in known_mouse_ids} if known_mouse_ids is not None else None
    biological: list[str] = []
    unmapped_biological: list[str] = []
    qc: list[str] = []
    blank: list[str] = []
    metadata: list[str] = []
    mouse_ids: dict[str, str] = {}
    tissues: dict[str, str] = {}

    for raw_col in columns:
        col = str(raw_col)
        if _QC_PATTERN.fullmatch(col):
            qc.append(col)
            continue
        if _BLANK_PATTERN.fullmatch(col):
            blank.append(col)
            continue
        match = _SAMPLE_PATTERN.fullmatch(col)
        if match:
            if known is None or match.group(1) in known:
                biological.append(col)
                mouse_ids[col] = match.group(1)
                tissues[col] = match.group(2)
            else:
                unmapped_biological.append(col)
        else:
            metadata.append(col)

    return DataColumnClassification(
        biological=biological,
        unmapped_biological=unmapped_biological,
        qc=qc,
        blank=blank,
        metadata=metadata,
        mouse_ids=mouse_ids,
        tissues=tissues,
    )


def identify_sample_columns(
    columns: pd.Index,
    known_mouse_ids: set[str] | list[str] | tuple[str, ...] | None = None,
) -> tuple[list[str], list[str]]:
    """兼容接口；QC/blank 仍归入非生物样本列。"""
    classified = classify_data_columns(columns, known_mouse_ids)
    non_biological = (
        classified.metadata + classified.unmapped_biological
        + classified.qc + classified.blank
    )
    return classified.biological, non_biological


def get_sample_groups(
    sample_cols: list[str],
    mouse_groups: dict[str, str],
    conditions: dict[str, dict[str, Any]],
    *,
    strict: bool = True,
) -> dict[str, str]:
    """
    根据配置文件的 mouse_groups 映射，将样本列映射到组别显示名。

    Returns: {sample_col_name: display_name}
    """
    cond_map = {}
    for code, val in conditions.items():
        cond_map[code] = val.get("display", code) if isinstance(val, dict) else code

    groups = {}
    unmapped: list[str] = []
    for col in sample_cols:
        m = _SAMPLE_PATTERN.fullmatch(col)
        if m:
            mouse_id = m.group(1)
            code = mouse_groups.get(mouse_id)
            if code:
                groups[col] = cond_map.get(code, code)
            else:
                unmapped.append(col)
                if not strict:
                    groups[col] = "Unknown"
    if strict and unmapped:
        raise ValueError(f"以下样本未在 mouse_groups 中配置: {unmapped}")
    return groups


def preprocess_data(
    df: pd.DataFrame,
    sample_cols: list[str],
    groups: dict[str, str] | None = None,
    *,
    qc_cols: list[str] | None = None,
    options: PreprocessingOptions | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """
    兼容接口；实际行为委托给 ``preprocess_metabolomics_data``。

    Returns: (processed_df, sample_cols)
    """
    result = preprocess_metabolomics_data(
        df,
        sample_cols,
        groups,
        qc_cols=qc_cols,
        options=options,
    )
    return result.processed_df, result.sample_cols


_FDR_ALIASES = {
        "bh": "fdr_bh",
        "benjamini-hochberg": "fdr_bh",
        "benjamini_hochberg": "fdr_bh",
        "by": "fdr_by",
        "benjamini-yekutieli": "fdr_by",
        "benjamini_yekutieli": "fdr_by",
        "holm-sidak": "holm-sidak",
}
_FDR_METHODS = {
    "bonferroni", "sidak", "holm-sidak", "holm", "simes-hochberg",
    "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky",
}


def normalize_fdr_method(method: str) -> str:
    """规范化并验证多重检验方法名称。"""
    normalized_method = str(method).strip().lower()
    normalized_method = _FDR_ALIASES.get(normalized_method, normalized_method)
    if normalized_method not in _FDR_METHODS:
        raise ValueError(
            f"不支持的 FDR/多重检验校正方法: {method!r}；"
            f"可选 {sorted(_FDR_METHODS)}"
        )
    return normalized_method


def apply_fdr_correction(
    p_values: np.ndarray | list[float],
    method: str = "bh",
) -> np.ndarray:
    """
    对有限 p 值进行多重检验校正，NaN 保持原位且不计入检验数。

    Parameters
    ----------
    p_values : array-like
        原始 p 值
    method : str
        statsmodels 支持的方法或常用别名，默认 "bh"。

    Returns
    -------
    np.ndarray
        校正后的 FDR 值
    """
    normalized_method = normalize_fdr_method(method)

    p_arr = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_arr.shape, np.nan, dtype=float)
    finite = np.isfinite(p_arr)
    if not finite.any():
        return adjusted
    if ((p_arr[finite] < 0) | (p_arr[finite] > 1)).any():
        raise ValueError("p 值必须位于 [0, 1] 范围内")

    from statsmodels.stats.multitest import multipletests

    adjusted[finite] = multipletests(
        p_arr[finite], method=normalized_method,
    )[1]
    return adjusted
