"""
热力图绘制。
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from . import _config  # noqa: F401


def select_heatmap_features(
    df: pd.DataFrame,
    sample_cols: list[str],
    *,
    top_n: int = 50,
    method: str = "variance",
    groups: dict[str, str] | None = None,
) -> pd.DataFrame:
    """选择热图特征；默认方差法不使用实验组标签。"""
    if top_n < 1:
        raise ValueError("heatmap_top_n 必须大于 0")
    matrix = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    normalized_method = str(method).strip().lower()
    if normalized_method == "variance":
        scores = matrix.var(axis=1, ddof=1)
        score_name = "Variance"
        order = scores.sort_values(ascending=False, kind="stable").index
    elif normalized_method == "anova":
        if not groups:
            raise ValueError("ANOVA 热图筛选需要 groups")
        unique_groups = list(dict.fromkeys(groups[col] for col in sample_cols))
        group_indices = {
            group: [col for col in sample_cols if groups[col] == group]
            for group in unique_groups
        }
        pvalues: list[float] = []
        for _, row in matrix.iterrows():
            values = [row[cols].dropna().to_numpy(dtype=float) for cols in group_indices.values()]
            if any(len(group_values) < 2 for group_values in values):
                pvalues.append(np.nan)
                continue
            try:
                with np.errstate(all="ignore"):
                    pvalues.append(float(stats.f_oneway(*values).pvalue))
            except (TypeError, ValueError):
                pvalues.append(np.nan)
        scores = pd.Series(pvalues, index=df.index, dtype=float)
        score_name = "ANOVA_P"
        order = scores.sort_values(ascending=True, kind="stable", na_position="last").index
    else:
        raise ValueError("heatmap_selection_method 仅支持 variance 或 anova")

    selected_index = order[: min(top_n, len(order))]
    selected = df.loc[selected_index].copy()
    selected[score_name] = scores.loc[selected_index]
    return selected


def plot_heatmap(
    df: pd.DataFrame,
    sample_cols: list[str],
    groups: dict[str, str],
    tissue: str,
    output_dir: str,
    top_n: int = 50,
    group_colors: dict[str, str] | None = None,
    selection_method: str = "variance",
):
    """
    默认按无监督方差筛选 Top 代谢物并绘制聚类热力图。
    """
    print(f"    绘制热力图: {tissue} (top {top_n}, {selection_method})")

    unique_groups = sorted(set(groups.values()))
    top_df = select_heatmap_features(
        df,
        sample_cols,
        top_n=top_n,
        method=selection_method,
        groups=groups,
    )

    heatmap_data = top_df[sample_cols].astype(float)
    if "Compounds" in top_df.columns:
        labels = top_df["Compounds"].where(top_df["Compounds"].notna(), top_df["Index"])
    else:
        labels = top_df["Index"]
    labels = labels.astype(str)
    if labels.duplicated().any():
        labels = labels + " [" + top_df["Index"].astype(str) + "]"
    heatmap_data.index = labels
    # Z-score 行标准化
    row_mean = heatmap_data.mean(axis=1)
    row_std = heatmap_data.std(axis=1).replace(0, np.nan)
    heatmap_data = heatmap_data.sub(row_mean, axis=0).div(row_std, axis=0)
    heatmap_data = heatmap_data.dropna(how="all")
    if heatmap_data.empty:
        print("    热力图跳过: 所选代谢物均无样本间变异")
        return

    # 列颜色条
    if group_colors:
        lut = {g: group_colors.get(g, "#888888") for g in unique_groups}
    else:
        lut = dict(zip(unique_groups, sns.color_palette("viridis", len(unique_groups))))
    col_colors = pd.Series(groups).reindex(sample_cols).map(lut)

    try:
        g = sns.clustermap(
            heatmap_data.fillna(0),
            col_colors=col_colors,
            cmap="vlag",
            center=0,
            xticklabels=False,
            yticklabels=True,
            figsize=(12, 12),
            dendrogram_ratio=(0.1, 0.2),
            cbar_pos=(0.02, 0.32, 0.03, 0.2),
        )
        for label, color in lut.items():
            g.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
        g.ax_col_dendrogram.legend(loc="center", ncol=len(unique_groups), bbox_to_anchor=(0.5, 1.2))

        out_path = os.path.join(output_dir, tissue)
        os.makedirs(out_path, exist_ok=True)
        top_df.to_csv(
            os.path.join(out_path, f"Heatmap_Top{top_n}_features.csv"),
            index=False,
            encoding="utf-8-sig",
        )
        selection_label = "unsupervised variance" if selection_method == "variance" else "supervised ANOVA"
        g.fig.suptitle(f"{tissue}: Top {top_n} features ({selection_label})", y=1.02)
        plt.savefig(os.path.join(out_path, f"Heatmap_Top{top_n}.png"), dpi=150, bbox_inches="tight")
        plt.close()
        print(f"    热力图保存: {out_path}/Heatmap_Top{top_n}.png")
    except Exception as e:
        print(f"    热力图绘制失败: {e}")
        plt.close()
