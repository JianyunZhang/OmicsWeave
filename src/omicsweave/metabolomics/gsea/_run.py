"""
代谢组学 GSEA 分析 — 基于 gseapy.prerank 的 rank-based 富集分析。

将代谢物按 log2FC * -log10(pvalue) 排序，以 KEGG 通路为代谢物集，
进行 GSEA 分析。
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from ..kegg_pathways import build_metabolite_sets, get_pathway_name, is_pathway_in_scope


def build_ranked_metabolites(
    diff_df: pd.DataFrame,
    comp_name: str,
    *,
    metric: str = "welch_t",
) -> pd.Series:
    """构建只含有限值的预排序代谢物序列，默认使用 Welch t 统计量。"""
    if "Index" not in diff_df.columns:
        raise ValueError("GSEA 输入缺少 Index 列")
    normalized_metric = str(metric).strip().lower()
    t_col = f"Tstat_{comp_name}"
    p_col = f"Pval_{comp_name}"
    fc_col = f"Log2FC_{comp_name}"

    if normalized_metric in {"welch_t", "t", "statistic"}:
        if t_col not in diff_df.columns:
            normalized_metric = "signed_logp"
        else:
            values = pd.to_numeric(diff_df[t_col], errors="coerce")
    if normalized_metric in {"signed_logp", "log2fc_logp"}:
        missing = [col for col in (p_col, fc_col) if col not in diff_df.columns]
        if missing:
            raise ValueError(f"GSEA 排名输入缺少列: {missing}")
        pvalues = pd.to_numeric(diff_df[p_col], errors="coerce").clip(
            lower=np.finfo(float).tiny,
            upper=1.0,
        )
        fold_changes = pd.to_numeric(diff_df[fc_col], errors="coerce")
        values = np.sign(fold_changes) * (-np.log10(pvalues))
    elif normalized_metric not in {"welch_t", "t", "statistic"}:
        raise ValueError("gsea_rank_metric 仅支持 welch_t 或 signed_logp")

    ranked_table = pd.DataFrame(
        {
            "feature": diff_df["Index"].astype(str).to_numpy(),
            "rank": values.to_numpy(dtype=float),
        }
    )
    ranked_table = ranked_table[np.isfinite(ranked_table["rank"])]
    ranked_table["absolute_rank"] = ranked_table["rank"].abs()
    ranked_table = ranked_table.sort_values(
        "absolute_rank", ascending=False, kind="stable",
    ).drop_duplicates("feature", keep="first")
    ranked_table = ranked_table.sort_values("rank", ascending=False, kind="stable")
    return ranked_table.set_index("feature")["rank"].rename("rank")


def run_gsea(
    diff_df: pd.DataFrame,
    comp_name: str,
    tissue: str,
    output_dir: str,
    min_set_size: int = 3,
    max_set_size: int = 500,
    permutation_num: int = 1000,
    rank_metric: str = "welch_t",
    pathway_scope: str = "metabolism",
    exclude_overview: bool = True,
) -> pd.DataFrame | None:
    """
    对一组差异分析结果运行 GSEA。

    Parameters
    ----------
    diff_df : 差异分析结果 DataFrame
    comp_name : 比较名称，如 "HFD_vs_LFD"
    tissue : 组织名
    output_dir : 输出目录

    Returns
    -------
    pd.DataFrame | None — GSEA 结果表
    """
    print(f"    GSEA 分析: {comp_name}")

    p_col = f"Pval_{comp_name}"
    fc_col = f"Log2FC_{comp_name}"

    if p_col not in diff_df.columns or fc_col not in diff_df.columns:
        print(f"      跳过: 缺少 {p_col} 或 {fc_col}")
        return None

    try:
        import gseapy as gp
    except ImportError:
        print("      gseapy 未安装，跳过 GSEA。安装: pip install gseapy")
        return None

    df = diff_df.copy()
    rank_series = build_ranked_metabolites(df, comp_name, metric=rank_metric)

    if len(rank_series) < min_set_size:
        print(f"      排名代谢物不足 ({len(rank_series)})，跳过 GSEA")
        return None

    # 构建 KEGG 通路代谢物集
    if "kegg_map" not in df.columns:
        print("      缺少 kegg_map，跳过 GSEA")
        return None
    metabolite_sets = build_metabolite_sets(df["kegg_map"], df["Index"].astype(str))
    ranked_ids = set(rank_series.index)
    metabolite_sets = {
        pathway: members & ranked_ids
        for pathway, members in metabolite_sets.items()
        if is_pathway_in_scope(
            pathway, pathway_scope, exclude_overview=exclude_overview,
        )
        and min_set_size <= len(members & ranked_ids) <= max_set_size
    }

    if not metabolite_sets:
        print("      无有效 KEGG 通路集，跳过 GSEA")
        return None

    # ID 前缀防止同名通路覆盖。
    gene_sets = {}
    for ko, members in metabolite_sets.items():
        name = get_pathway_name(ko)
        gene_sets[f"{ko} | {name}"] = members

    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)

    try:
        res = gp.prerank(
            rnk=rank_series,
            gene_sets=gene_sets,
            min_size=min_set_size,
            max_size=max_set_size,
            permutation_num=permutation_num,
            seed=42,
            no_plot=True,
            outdir=None,
            verbose=False,
        )
    except Exception as e:
        print(f"      GSEA 运行失败: {e}")
        return None

    # 结果表
    res_df = res.res2d.copy() if res.res2d is not None else None
    if res_df is None or res_df.empty:
        print("      GSEA 无结果")
        return None

    # 保存
    csv_path = os.path.join(out_path, f"gsea_{comp_name}_results.csv")
    res_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    print(f"      GSEA 结果: {csv_path}")

    # 绘制 NES 柱状图
    _plot_gsea_bar(res_df, comp_name, tissue, out_path)

    return res_df


def _plot_gsea_bar(res_df: pd.DataFrame, comp_name: str, tissue: str, out_path: str):
    """绘制 GSEA NES 柱状图。"""
    import matplotlib.pyplot as plt

    # 取 Top 10 positive NES + Top 10 negative NES
    if "NES" not in res_df.columns:
        return
    res_df["NES"] = pd.to_numeric(res_df["NES"], errors="coerce")
    top_pos = res_df.nlargest(10, "NES")
    top_neg = res_df.nsmallest(10, "NES")
    plot_df = pd.concat([top_pos, top_neg]).drop_duplicates()

    if plot_df.empty:
        return

    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df) * 0.4)))
    colors = ["#D32F2F" if nes > 0 else "#1565C0" for nes in plot_df["NES"]]
    ax.barh(range(len(plot_df)), plot_df["NES"], color=colors, alpha=0.8)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["Term"].values, fontsize=8)
    ax.set_xlabel("NES (Normalized Enrichment Score)")
    ax.invert_yaxis()

    display_comp = comp_name.replace("_vs_", " vs ")
    ax.set_title(f"GSEA: {display_comp} ({tissue})")
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, f"gsea_{comp_name}.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"      GSEA 图保存: {out_path}/gsea_{comp_name}.png")
