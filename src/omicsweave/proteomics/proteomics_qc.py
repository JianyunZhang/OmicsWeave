"""DIA 蛋白质组学的质量控制指标与图表。"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests


def compute_sample_qc_metrics(intensity_matrix: pd.DataFrame) -> pd.DataFrame:
    """Calculate identification, missingness, intensity and CV per sample."""

    matrix = intensity_matrix.astype(float)
    metrics = pd.DataFrame(index=matrix.columns.astype(str))
    metrics.index.name = "sample"
    metrics["n_proteins"] = matrix.notna().sum(axis=0).to_numpy(dtype=int)
    metrics["n_missing"] = matrix.isna().sum(axis=0).to_numpy(dtype=int)
    metrics["missing_pct"] = 100.0 * metrics["n_missing"] / max(len(matrix), 1)
    metrics["median_intensity"] = matrix.median(axis=0, skipna=True).to_numpy()
    metrics["mean_intensity"] = matrix.mean(axis=0, skipna=True).to_numpy()
    metrics["sd_intensity"] = matrix.std(axis=0, skipna=True, ddof=1).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        metrics["cv_pct"] = (
            100.0 * metrics["sd_intensity"] / metrics["mean_intensity"]
        )
    metrics.loc[~np.isfinite(metrics["cv_pct"]), "cv_pct"] = np.nan
    return metrics


def compute_replicate_correlations(
    normalized_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate all pairwise-complete correlations and identify within-group pairs."""

    correlation = normalized_matrix.corr(method="pearson", min_periods=3)
    rows: list[dict[str, Any]] = []
    for sample_a, sample_b in combinations(correlation.columns.astype(str), 2):
        condition_a = str(metadata.loc[sample_a, "condition"])
        condition_b = str(metadata.loc[sample_b, "condition"])
        rows.append({
            "sample_a": sample_a,
            "sample_b": sample_b,
            "condition_a": condition_a,
            "condition_b": condition_b,
            "within_condition": condition_a == condition_b,
            "pearson_r": correlation.loc[sample_a, sample_b],
            "n_pairwise_complete": int(
                normalized_matrix[[sample_a, sample_b]].notna().all(axis=1).sum()
            ),
        })
    return pd.DataFrame(rows)


def _median_impute_for_qc(matrix: pd.DataFrame) -> pd.DataFrame:
    row_medians = matrix.median(axis=1, skipna=True)
    imputed = matrix.T.fillna(row_medians).T
    finite_values = imputed.to_numpy(dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    fallback = float(np.median(finite_values)) if finite_values.size else 0.0
    return imputed.fillna(fallback)


def _pc_covariate_association(
    normalized_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    imputed = _median_impute_for_qc(normalized_matrix)
    variable = imputed.var(axis=1, ddof=0) > 0
    imputed = imputed.loc[variable]
    n_components = min(5, imputed.shape[0], max(imputed.shape[1] - 1, 0))
    if n_components < 1:
        return pd.DataFrame(), pd.DataFrame(index=normalized_matrix.columns)

    pca = PCA(n_components=n_components, random_state=seed)
    scores = pca.fit_transform(imputed.T)
    score_columns = [f"PC{i + 1}" for i in range(n_components)]
    score_df = pd.DataFrame(scores, index=imputed.columns, columns=score_columns)
    for column in metadata.columns:
        score_df[column] = metadata.loc[score_df.index, column]

    ignored = {"sample", "mouse_id", "tissue"}
    rows: list[dict[str, Any]] = []
    for pc_index, pc in enumerate(score_columns):
        for covariate in metadata.columns:
            if covariate in ignored:
                continue
            values = metadata.loc[score_df.index, covariate]
            if values.isna().any() or values.nunique(dropna=False) < 2:
                continue
            if pd.api.types.is_numeric_dtype(values) and values.nunique() > 2:
                statistic, pvalue = stats.spearmanr(
                    score_df[pc].to_numpy(dtype=float), values.to_numpy(dtype=float),
                )
                test = "spearman"
            else:
                groups = [
                    score_df.loc[values.astype(str) == level, pc].to_numpy(dtype=float)
                    for level in sorted(values.astype(str).unique())
                ]
                groups = [group for group in groups if len(group) > 0]
                if len(groups) < 2:
                    continue
                statistic, pvalue = stats.f_oneway(*groups)
                test = "one_way_anova"
            rows.append({
                "pc": pc,
                "explained_variance_ratio": float(pca.explained_variance_ratio_[pc_index]),
                "covariate": covariate,
                "test": test,
                "statistic": float(statistic),
                "pvalue": float(pvalue),
            })
    association = pd.DataFrame(rows)
    if not association.empty:
        finite = np.isfinite(association["pvalue"])
        association["padj_bh"] = np.nan
        association.loc[finite, "padj_bh"] = multipletests(
            association.loc[finite, "pvalue"], method="fdr_bh",
        )[1]
    return association, score_df


def _save_intensity_distributions(matrix: pd.DataFrame, output_path: Path) -> None:
    log_matrix = np.log2(matrix.where(matrix > 0))
    long = log_matrix.rename_axis("protein").reset_index().melt(
        id_vars="protein", var_name="sample", value_name="log2_intensity",
    ).dropna()
    plt.figure(figsize=(max(8, 0.75 * matrix.shape[1]), 6))
    sns.violinplot(
        data=long,
        x="sample",
        y="log2_intensity",
        inner="quartile",
        cut=0,
        color="#4C72B0",
    )
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("样本")
    plt.ylabel("Log2 强度")
    plt.title("蛋白强度分布")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()


def _save_detection_missingness(metrics: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(max(8, 0.75 * len(metrics)), 8), sharex=True)
    sns.barplot(x=metrics.index, y=metrics["n_proteins"], ax=axes[0], color="#4C72B0")
    axes[0].set_ylabel("鉴定蛋白数")
    axes[0].set_xlabel("")
    sns.barplot(x=metrics.index, y=metrics["missing_pct"], ax=axes[1], color="#DD8452")
    axes[1].set_ylabel("缺失比例 (%)")
    axes[1].set_xlabel("样本")
    axes[1].tick_params(axis="x", rotation=45)
    fig.suptitle("样本鉴定数与缺失率")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _save_correlation_heatmap(matrix: pd.DataFrame, output_path: Path) -> None:
    correlation = matrix.corr(method="pearson", min_periods=3)
    plt.figure(figsize=(max(7, 0.7 * len(correlation)), max(6, 0.65 * len(correlation))))
    sns.heatmap(
        correlation,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        annot=len(correlation) <= 10,
        fmt=".2f",
    )
    plt.title("样本相关性（pairwise complete Pearson r）")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()


def _save_missingness_heatmap(matrix: pd.DataFrame, output_path: Path) -> None:
    missing_fraction = matrix.isna().mean(axis=1)
    selected = missing_fraction[(missing_fraction > 0) & (missing_fraction < 1)].sort_values(
        ascending=False,
    ).head(200).index
    if len(selected) == 0:
        selected = matrix.index[: min(50, len(matrix))]
    binary = matrix.loc[selected].isna().astype(int)
    plt.figure(figsize=(max(7, 0.7 * matrix.shape[1]), 8))
    sns.heatmap(binary, cmap=["#F2F2F2", "#D32F2F"], cbar=False, yticklabels=False)
    plt.xlabel("样本")
    plt.ylabel("蛋白（最多 200 个）")
    plt.title("缺失值模式（红色=缺失）")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()


def generate_qc_outputs(
    raw_intensity_matrix: pd.DataFrame,
    normalized_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: str | Path,
    *,
    seed: int = 42,
) -> dict[str, Any]:
    """Generate reusable QC tables, plots and a machine-readable summary."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    metadata = metadata.copy()
    metadata.index = metadata.index.astype(str)
    samples = normalized_matrix.columns.astype(str).tolist()
    metadata = metadata.loc[samples]

    sample_metrics = compute_sample_qc_metrics(raw_intensity_matrix)
    sample_metrics = sample_metrics.join(metadata, how="left")
    sample_metrics.to_csv(output / "sample_qc_metrics.csv", encoding="utf-8-sig")
    metadata.to_csv(output / "sample_metadata.csv", encoding="utf-8-sig")

    protein_missingness = pd.DataFrame({
        "n_observed": normalized_matrix.notna().sum(axis=1),
        "n_missing": normalized_matrix.isna().sum(axis=1),
        "missing_pct": 100.0 * normalized_matrix.isna().mean(axis=1),
    })
    protein_missingness.index.name = "Protein"
    protein_missingness.to_csv(output / "protein_missingness.csv", encoding="utf-8-sig")

    replicate_correlations = compute_replicate_correlations(normalized_matrix, metadata)
    replicate_correlations.to_csv(
        output / "replicate_correlations.csv", index=False, encoding="utf-8-sig",
    )

    group_rows: list[dict[str, Any]] = []
    for condition, group_metadata in metadata.groupby("condition", dropna=False):
        group_samples = group_metadata.index.astype(str).tolist()
        group_matrix = raw_intensity_matrix[group_samples]
        linear_cv = 100.0 * group_matrix.std(axis=1, ddof=1) / group_matrix.mean(axis=1)
        within = replicate_correlations[
            replicate_correlations["within_condition"]
            & (replicate_correlations["condition_a"].astype(str) == str(condition))
        ]
        group_rows.append({
            "condition": condition,
            "n_samples": len(group_samples),
            "median_proteins_detected": float(group_matrix.notna().sum(axis=0).median()),
            "median_protein_cv_pct": float(linear_cv.replace([np.inf, -np.inf], np.nan).median()),
            "median_within_group_pearson_r": float(within["pearson_r"].median()) if len(within) else np.nan,
        })
    pd.DataFrame(group_rows).to_csv(
        output / "group_qc_metrics.csv", index=False, encoding="utf-8-sig",
    )

    association, scores = _pc_covariate_association(normalized_matrix, metadata, seed)
    association.to_csv(
        output / "pc_covariate_association.csv", index=False, encoding="utf-8-sig",
    )
    scores.to_csv(output / "pca_scores.csv", encoding="utf-8-sig")

    _save_intensity_distributions(raw_intensity_matrix, output / "intensity_distributions.png")
    _save_detection_missingness(sample_metrics, output / "sample_detection_missingness.png")
    _save_correlation_heatmap(normalized_matrix, output / "correlation_heatmap.png")
    _save_missingness_heatmap(normalized_matrix, output / "missingness_heatmap.png")

    all_correlations = replicate_correlations["pearson_r"].dropna()
    within_correlations = replicate_correlations.loc[
        replicate_correlations["within_condition"], "pearson_r",
    ].dropna()
    overall_missing = float(100.0 * normalized_matrix.isna().to_numpy().mean())
    summary: dict[str, Any] = {
        "n_samples": int(normalized_matrix.shape[1]),
        "n_proteins_before_filter": int(raw_intensity_matrix.shape[0]),
        "n_proteins_after_filter": int(normalized_matrix.shape[0]),
        "overall_missing_pct": overall_missing,
        "median_proteins_per_sample": float(sample_metrics["n_proteins"].median()),
        "median_all_pairwise_correlation": float(all_correlations.median()) if len(all_correlations) else None,
        "median_within_group_correlation": float(within_correlations.median()) if len(within_correlations) else None,
        "flags": [],
    }
    if overall_missing > 30:
        summary["flags"].append("high_missingness_gt_30pct")
    if len(within_correlations) and float(within_correlations.median()) < 0.8:
        summary["flags"].append("low_within_group_correlation_lt_0.8")
    if not association.empty:
        batch_hits = association[
            (association["covariate"].astype(str).str.lower() == "batch")
            & (association["padj_bh"] < 0.05)
        ]
        if len(batch_hits):
            summary["flags"].append("batch_associated_with_principal_components")

    (output / "qc_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    return summary


__all__ = [
    "compute_replicate_correlations",
    "compute_sample_qc_metrics",
    "generate_qc_outputs",
]
