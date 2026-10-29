"""代谢组学差异分析：统一检验、FDR、VIP 与显著性输出口径。"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .data_preprocessing import apply_fdr_correction


@dataclass(frozen=True)
class DifferentialAnalysisOptions:
    pvalue_cutoff: float = 0.05
    fdr_cutoff: float = 0.05
    log2fc_cutoff: float = 0.58
    fdr_method: str = "bh"
    vip_cutoff: float = 1.0
    require_vip_for_significance: bool = False
    min_observations_per_group: int = 2


@dataclass
class DifferentialAnalysisResult:
    combined: pd.DataFrame
    contrasts: dict[str, pd.DataFrame]
    design_report: pd.DataFrame


_DEFAULT_METADATA_COLUMNS = (
    "Index", "Compounds", "Formula", "kegg_map", "HMDB", "KEGG_Compound",
    "ChEBI", "Lipidmaps ID", "CAS", "PubChem CID", "cpd_ID", "Mode",
    "Ionization model", "Molecular weight (Da)", "Q1 (Da)",
)


def _coerce_vip(
    value: pd.Series | np.ndarray | list[float] | None,
    index: pd.Index,
) -> pd.Series:
    if value is None:
        return pd.Series(np.nan, index=index, dtype=float)
    if isinstance(value, pd.Series):
        if value.index.equals(index):
            return pd.to_numeric(value, errors="coerce").astype(float)
        if len(value) != len(index):
            raise ValueError("VIP 分数长度与代谢物数不一致")
        return pd.Series(pd.to_numeric(value, errors="coerce").to_numpy(), index=index)
    array = np.asarray(value, dtype=float)
    if len(array) != len(index):
        raise ValueError("VIP 分数长度与代谢物数不一致")
    return pd.Series(array, index=index, dtype=float)


def analyze_differential_metabolites(
    processed_df: pd.DataFrame,
    groups: dict[str, str],
    comparisons: list[tuple[str, str, str]],
    *,
    options: DifferentialAnalysisOptions | None = None,
    contrast_vip: dict[str, pd.Series | np.ndarray | list[float]] | None = None,
    global_vip: pd.Series | np.ndarray | list[float] | None = None,
    metadata_columns: list[str] | None = None,
) -> DifferentialAnalysisResult:
    """
    在 log2 归一化矩阵上运行向量化 Welch t 检验。

    ``significant`` 是唯一供火山图、ORA、报告和下游 CSV 使用的主判定列：
    默认要求 FDR 和 |log2FC| 同时通过；仅在配置明确开启时再要求 VIP。
    """
    opts = options or DifferentialAnalysisOptions()
    if not 0 < opts.pvalue_cutoff <= 1 or not 0 < opts.fdr_cutoff <= 1:
        raise ValueError("pvalue_cutoff 和 fdr_cutoff 必须在 (0, 1] 范围内")
    if opts.log2fc_cutoff < 0 or opts.vip_cutoff < 0:
        raise ValueError("log2fc_cutoff 和 vip_cutoff 不得为负数")
    if opts.min_observations_per_group < 2:
        raise ValueError("min_observations_per_group 至少为 2")

    sample_cols = list(groups)
    missing = [col for col in sample_cols if col not in processed_df.columns]
    if missing:
        raise ValueError(f"差异分析缺少样本列: {missing}")

    if metadata_columns is None:
        metadata_columns = [col for col in _DEFAULT_METADATA_COLUMNS if col in processed_df.columns]
    if "Index" not in metadata_columns and "Index" in processed_df.columns:
        metadata_columns = ["Index", *metadata_columns]
    combined = processed_df[metadata_columns].copy()
    if global_vip is not None:
        combined["VIP_Global"] = _coerce_vip(global_vip, processed_df.index).to_numpy()

    contrast_tables: dict[str, pd.DataFrame] = {}
    design_rows: list[dict[str, Any]] = []
    vip_map = contrast_vip or {}

    for group_a, group_b, comp_name in comparisons:
        cols_a = [col for col, group in groups.items() if group == group_a]
        cols_b = [col for col, group in groups.items() if group == group_b]
        if not cols_a or not cols_b:
            design_rows.append(
                {
                    "comparison": comp_name,
                    "group_a": group_a,
                    "group_b": group_b,
                    "n_group_a": len(cols_a),
                    "n_group_b": len(cols_b),
                    "test": "not_run",
                    "route": "missing_group",
                    "fdr_method": opts.fdr_method,
                    "require_vip": opts.require_vip_for_significance,
                }
            )
            continue

        values_a = processed_df[cols_a].to_numpy(dtype=float)
        values_b = processed_df[cols_b].to_numpy(dtype=float)
        log2fc = np.nanmean(values_a, axis=1) - np.nanmean(values_b, axis=1)
        enough = (
            len(cols_a) >= opts.min_observations_per_group
            and len(cols_b) >= opts.min_observations_per_group
        )
        if enough:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                statistic, pvalues = stats.ttest_ind(
                    values_a,
                    values_b,
                    axis=1,
                    equal_var=False,
                    nan_policy="omit",
                )
            pvalues = np.asarray(pvalues, dtype=float)
            statistic = np.asarray(statistic, dtype=float)
            padj = apply_fdr_correction(pvalues, opts.fdr_method)
            route = "inferential"
            evidence = np.repeat("inferential", len(processed_df))
        else:
            statistic = np.full(len(processed_df), np.nan)
            pvalues = np.full(len(processed_df), np.nan)
            padj = np.full(len(processed_df), np.nan)
            route = "insufficient_data"
            evidence = np.repeat("insufficient", len(processed_df))

        vip = _coerce_vip(vip_map.get(comp_name), processed_df.index)
        fc_pass = np.isfinite(log2fc) & (np.abs(log2fc) > opts.log2fc_cutoff)
        nominal = np.isfinite(pvalues) & (pvalues < opts.pvalue_cutoff) & fc_pass
        univariate = np.isfinite(padj) & (padj < opts.fdr_cutoff) & fc_pass
        vip_pass = np.isfinite(vip.to_numpy()) & (vip.to_numpy() >= opts.vip_cutoff)
        with_vip = univariate & vip_pass
        primary = with_vip if opts.require_vip_for_significance else univariate

        combined[f"Tstat_{comp_name}"] = statistic
        combined[f"Pval_{comp_name}"] = pvalues
        combined[f"FDR_{comp_name}"] = padj
        combined[f"Log2FC_{comp_name}"] = log2fc
        combined[f"VIP_{comp_name}"] = vip.to_numpy()
        combined[f"NominalSignificant_{comp_name}"] = nominal
        combined[f"SignificantUnivariate_{comp_name}"] = univariate
        combined[f"SignificantWithVIP_{comp_name}"] = with_vip
        combined[f"Significant_{comp_name}"] = primary
        combined[f"Evidence_{comp_name}"] = evidence

        contrast = combined[metadata_columns].copy()
        if "VIP_Global" in combined.columns:
            contrast["VIP_Global"] = combined["VIP_Global"]
        contrast[f"VIP_{comp_name}"] = vip.to_numpy()
        contrast["log2FoldChange"] = log2fc
        contrast["statistic"] = statistic
        contrast["pvalue"] = pvalues
        contrast["padj"] = padj
        contrast["nominal_significant"] = nominal
        contrast["significant_univariate"] = univariate
        contrast["vip_pass"] = vip_pass
        contrast["significant_with_vip"] = with_vip
        contrast["significant"] = primary
        contrast["evidence_level"] = evidence
        contrast_tables[comp_name] = contrast

        design_rows.append(
            {
                "comparison": comp_name,
                "group_a": group_a,
                "group_b": group_b,
                "n_group_a": len(cols_a),
                "n_group_b": len(cols_b),
                "test": "welch_t_test" if enough else "not_run",
                "route": route,
                "fdr_method": opts.fdr_method,
                "pvalue_cutoff": opts.pvalue_cutoff,
                "fdr_cutoff": opts.fdr_cutoff,
                "log2fc_cutoff": opts.log2fc_cutoff,
                "vip_cutoff": opts.vip_cutoff,
                "require_vip": opts.require_vip_for_significance,
            }
        )

    return DifferentialAnalysisResult(
        combined=combined,
        contrasts=contrast_tables,
        design_report=pd.DataFrame(design_rows),
    )
