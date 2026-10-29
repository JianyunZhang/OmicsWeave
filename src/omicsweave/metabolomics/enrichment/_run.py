"""KEGG 通路过度代表分析（单侧 Fisher / ORA）。"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import stats

from ..data_preprocessing import apply_fdr_correction
from ..kegg_pathways import build_metabolite_sets, get_pathway_name, is_pathway_in_scope
from ..plotting.pathway_plots import plot_pathway_bar, plot_pathway_dot


def calculate_pathway_enrichment(
    diff_df: pd.DataFrame,
    *,
    significance_col: str,
    pathway_col: str = "kegg_map",
    feature_col: str = "Index",
    min_set_size: int = 3,
    max_set_size: int = 500,
    pathway_scope: str = "metabolism",
    exclude_overview: bool = True,
    fdr_method: str = "bh",
) -> pd.DataFrame:
    """
    计算 ORA；显著集与背景集严格限制在同一可注释 KEGG universe 内。

    所有合格通路（包括 SigCount=0）都纳入多重检验分母。返回表按 FDR、
    P-value 排序，便于调用者决定显示阈值而不丢失完整统计结果。
    """
    required = {significance_col, pathway_col, feature_col}
    missing = sorted(required - set(diff_df.columns))
    if missing:
        raise ValueError(f"ORA 输入缺少列: {missing}")
    if min_set_size < 1 or max_set_size < min_set_size:
        raise ValueError("ORA 集合大小参数无效")

    all_sets = build_metabolite_sets(
        diff_df[pathway_col], diff_df[feature_col].astype(str),
    )
    scoped_sets = {
        pathway: members
        for pathway, members in all_sets.items()
        if is_pathway_in_scope(
            pathway, pathway_scope, exclude_overview=exclude_overview,
        )
    }
    universe = set().union(*scoped_sets.values()) if scoped_sets else set()
    tested_sets = {
        pathway: members & universe
        for pathway, members in scoped_sets.items()
        if min_set_size <= len(members & universe) <= max_set_size
    }
    if not universe or not tested_sets:
        return pd.DataFrame(
            columns=[
                "ID", "Name", "P-value", "FDR", "OddsRatio", "SigCount",
                "BgCount", "SigTotal", "UniverseCount", "RichFactor", "Members",
            ]
        )

    significant_mask = diff_df[significance_col].fillna(False).astype(bool)
    significant_ids = set(diff_df.loc[significant_mask, feature_col].astype(str)) & universe
    total_bg = len(universe)
    total_sig = len(significant_ids)
    rows: list[dict[str, object]] = []

    for pathway, members in tested_sets.items():
        members = members & universe
        overlap = significant_ids & members
        a = len(overlap)
        b = total_sig - a
        c = len(members) - a
        d = total_bg - a - b - c
        if min(a, b, c, d) < 0:
            raise ValueError(f"ORA 列联表无效: {pathway} -> {[a, b, c, d]}")
        odds_ratio, pvalue = stats.fisher_exact(
            [[a, b], [c, d]], alternative="greater",
        )
        rows.append(
            {
                "ID": pathway,
                "Name": get_pathway_name(pathway),
                "P-value": float(pvalue),
                "OddsRatio": float(odds_ratio),
                "SigCount": a,
                "BgCount": len(members),
                "SigTotal": total_sig,
                "UniverseCount": total_bg,
                "RichFactor": a / len(members) if members else np.nan,
                "Members": ";".join(sorted(overlap)),
            }
        )

    result = pd.DataFrame(rows)
    result["FDR"] = apply_fdr_correction(result["P-value"].to_numpy(), fdr_method)
    return result.sort_values(["FDR", "P-value", "ID"], kind="stable").reset_index(drop=True)


def run_pathway_enrichment(
    diff_df: pd.DataFrame,
    comp_name: str,
    tissue: str,
    output_dir: str,
    pval_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.58,
    *,
    fdr_method: str = "bh",
    min_set_size: int = 3,
    max_set_size: int = 500,
    pathway_scope: str = "metabolism",
    exclude_overview: bool = True,
) -> pd.DataFrame | None:
    """运行并保存 ORA；优先消费差异分析生成的统一 Significant 列。"""
    print(f"    通路富集分析: {comp_name}")
    sig_col = f"Significant_{comp_name}"
    if sig_col not in diff_df.columns:
        fdr_col = f"FDR_{comp_name}"
        fc_col = f"Log2FC_{comp_name}"
        if fdr_col not in diff_df.columns or fc_col not in diff_df.columns:
            print(f"      跳过: 缺少 {sig_col} 或 FDR/Log2FC 列")
            return None
        diff_df = diff_df.copy()
        diff_df[sig_col] = (
            (diff_df[fdr_col] < pval_cutoff)
            & (diff_df[fc_col].abs() > log2fc_cutoff)
        )

    result = calculate_pathway_enrichment(
        diff_df,
        significance_col=sig_col,
        min_set_size=min_set_size,
        max_set_size=max_set_size,
        pathway_scope=pathway_scope,
        exclude_overview=exclude_overview,
        fdr_method=fdr_method,
    )
    if result.empty:
        print("      无可检验的 KEGG 代谢通路")
        return None

    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)
    csv_path = os.path.join(out_path, f"Pathway_{comp_name}.csv")
    result.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"      通路结果: {csv_path}")
    plot_pathway_bar(result, comp_name, tissue, output_dir)
    plot_pathway_dot(result, comp_name, tissue, output_dir)
    return result
