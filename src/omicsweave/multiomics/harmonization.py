#!/usr/bin/env python3
"""样本级多组学协调：归一化、样本映射、缺失审计与批次校正。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import yaml


REQUIRED_SAMPLE_MAP_COLUMNS = {
    "assay", "assay_sample_id", "primary_sample_id", "batch_id", "condition",
}


@dataclass(frozen=True)
class HarmonizedAssays:
    assays: dict[str, pd.DataFrame]
    sample_metadata: pd.DataFrame
    sample_map: pd.DataFrame
    missingness: pd.DataFrame
    primary_samples: tuple[str, ...]


@dataclass(frozen=True)
class ResidualBatchDiagnostic:
    """单个 assay、单个处理阶段的残余批次审计结果。"""

    summary: dict[str, Any]
    feature_qc: pd.DataFrame
    pca_scores: pd.DataFrame
    pca_variance: pd.DataFrame


_BATCH_FEATURE_COLUMNS = [
    "assay", "stage", "feature", "n_valid", "batch_partial_r2",
    "sse_condition", "sse_condition_batch", "status",
]


def _categorical_design(metadata: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    parts = [
        pd.DataFrame(
            {"intercept": np.ones(len(metadata), dtype=float)},
            index=metadata.index,
        )
    ]
    for column in columns:
        encoded = pd.get_dummies(
            metadata[column].astype(str),
            prefix=column,
            drop_first=True,
            dtype=float,
        )
        parts.append(encoded)
    return pd.concat(parts, axis=1)


def _batch_partial_r2(
    values: np.ndarray,
    metadata: pd.DataFrame,
    *,
    batch_column: str,
    condition_column: str,
) -> tuple[float, float, float, str]:
    """计算控制 condition 后 batch 对一个连续变量的增量解释率。"""
    valid = np.isfinite(values)
    one_meta = metadata.iloc[np.flatnonzero(valid)]
    y_values = values[valid]
    if len(y_values) == 0:
        return float("nan"), float("nan"), float("nan"), "no_finite_values"
    if one_meta[batch_column].nunique(dropna=True) < 2:
        return float("nan"), float("nan"), float("nan"), "single_batch_subset"

    reduced = _categorical_design(one_meta, [condition_column])
    full = _categorical_design(one_meta, [condition_column, batch_column])
    reduced_values = reduced.to_numpy(dtype=float)
    full_values = full.to_numpy(dtype=float)
    reduced_rank = int(np.linalg.matrix_rank(reduced_values))
    full_rank = int(np.linalg.matrix_rank(full_values))
    if reduced_rank < reduced_values.shape[1] or full_rank < full_values.shape[1]:
        return float("nan"), float("nan"), float("nan"), "rank_deficient"
    if len(y_values) <= full_rank:
        return float("nan"), float("nan"), float("nan"), "insufficient_df"

    reduced_coef, *_ = np.linalg.lstsq(reduced_values, y_values, rcond=None)
    full_coef, *_ = np.linalg.lstsq(full_values, y_values, rcond=None)
    reduced_residual = y_values - reduced_values @ reduced_coef
    full_residual = y_values - full_values @ full_coef
    sse_reduced = float(np.sum(reduced_residual ** 2))
    sse_full = float(np.sum(full_residual ** 2))
    tolerance = np.finfo(float).eps * max(float(np.sum(y_values ** 2)), 1.0) * 100.0
    if sse_reduced <= tolerance:
        partial_r2 = 0.0 if sse_full <= tolerance else float("nan")
    else:
        partial_r2 = float(np.clip((sse_reduced - sse_full) / sse_reduced, 0.0, 1.0))
    return partial_r2, sse_reduced, sse_full, "assessable"


def _pca_batch_diagnostic(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    assay: str,
    stage: str,
    batch_column: str,
    condition_column: str,
    max_components: int,
    max_features: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    numeric = matrix.apply(pd.to_numeric, errors="coerce").astype(float)
    variances = numeric.var(axis=1, skipna=True, ddof=1)
    usable = variances[np.isfinite(variances) & variances.gt(0)].sort_values(
        ascending=False
    )
    if max_features > 0:
        usable = usable.head(max_features)
    if usable.empty or numeric.shape[1] < 2:
        return pd.DataFrame(), pd.DataFrame()

    sample_feature = numeric.loc[usable.index].T.copy()
    sample_feature = sample_feature.apply(
        lambda column: column.fillna(column.median()), axis=0
    )
    sample_feature = sample_feature.dropna(axis=1, how="any")
    if sample_feature.empty:
        return pd.DataFrame(), pd.DataFrame()
    centered = sample_feature - sample_feature.mean(axis=0)
    scales = sample_feature.std(axis=0, ddof=1).replace(0, np.nan)
    standardized = centered.divide(scales, axis=1).dropna(axis=1, how="any")
    if standardized.empty:
        return pd.DataFrame(), pd.DataFrame()

    u_values, singular_values, _ = np.linalg.svd(
        standardized.to_numpy(dtype=float), full_matrices=False
    )
    n_components = min(
        max(int(max_components), 1),
        len(singular_values),
        max(len(standardized) - 1, 1),
    )
    scores = u_values[:, :n_components] * singular_values[:n_components]
    score_table = pd.DataFrame(
        scores,
        index=standardized.index.astype(str),
        columns=[f"PC{index + 1}" for index in range(n_components)],
    )
    score_table.insert(0, "sample_id", score_table.index)
    score_table.insert(0, "stage", stage)
    score_table.insert(0, "assay", assay)
    score_table[condition_column] = metadata.loc[score_table.index, condition_column].to_numpy()
    score_table[batch_column] = metadata.loc[score_table.index, batch_column].to_numpy()
    score_table = score_table.reset_index(drop=True)

    total = float(np.sum(singular_values ** 2))
    variance_rows: list[dict[str, Any]] = []
    for component_index in range(n_components):
        component = f"PC{component_index + 1}"
        partial_r2, _, _, status = _batch_partial_r2(
            scores[:, component_index],
            metadata.loc[standardized.index],
            batch_column=batch_column,
            condition_column=condition_column,
        )
        variance_rows.append({
            "assay": assay,
            "stage": stage,
            "component": component,
            "explained_variance_ratio": (
                float(singular_values[component_index] ** 2 / total)
                if total > 0 else float("nan")
            ),
            "batch_partial_r2": partial_r2,
            "status": status,
        })
    return score_table, pd.DataFrame(variance_rows)


def diagnose_residual_batch(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    assay: str,
    stage: str,
    batch_column: str = "batch_id",
    condition_column: str = "condition",
    max_components: int = 5,
    max_features: int = 5000,
    high_batch_partial_r2: float = 0.20,
) -> ResidualBatchDiagnostic:
    """审计控制 condition 后的残余 batch 方差，并导出 PCA 坐标。

    ``suggested_batch_reliability`` 仅是供人工审核的 0--1 建议值；本函数
    不会自动改变 EvidenceProfile 分数。
    """
    required = {batch_column, condition_column}
    missing = required - set(sample_metadata.columns)
    if missing:
        raise ValueError("sample_metadata missing columns: " + ", ".join(sorted(missing)))
    numeric = matrix.copy()
    numeric.columns = numeric.columns.astype(str)
    metadata = sample_metadata.copy()
    metadata.index = metadata.index.astype(str)
    samples = list(numeric.columns)
    if not set(samples).issubset(set(metadata.index)):
        raise ValueError("sample_metadata does not cover all assay columns")
    metadata = metadata.loc[samples]

    base_summary: dict[str, Any] = {
        "assay": str(assay),
        "stage": str(stage),
        "n_samples": len(samples),
        "n_features": int(numeric.shape[0]),
        "n_conditions": int(metadata[condition_column].nunique(dropna=True)),
        "n_batches": int(metadata[batch_column].nunique(dropna=True)),
        "n_features_assessed": 0,
        "feature_batch_partial_r2_median": float("nan"),
        "feature_batch_partial_r2_p90": float("nan"),
        "fraction_features_high_batch": float("nan"),
        "pca_batch_partial_r2_max": float("nan"),
        "suggested_batch_reliability": float("nan"),
        "status": "unassessable",
    }
    if len(samples) < 3:
        base_summary["status"] = "insufficient_samples"
        return ResidualBatchDiagnostic(
            base_summary,
            pd.DataFrame(columns=_BATCH_FEATURE_COLUMNS),
            pd.DataFrame(),
            pd.DataFrame(),
        )
    if base_summary["n_batches"] < 2:
        base_summary["status"] = "single_batch"
        return ResidualBatchDiagnostic(
            base_summary,
            pd.DataFrame(columns=_BATCH_FEATURE_COLUMNS),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    global_design = _categorical_design(
        metadata, [condition_column, batch_column]
    ).to_numpy(dtype=float)
    if np.linalg.matrix_rank(global_design) < global_design.shape[1]:
        base_summary.update({
            "status": "confounded",
            "suggested_batch_reliability": 0.25,
        })
        return ResidualBatchDiagnostic(
            base_summary,
            pd.DataFrame(columns=_BATCH_FEATURE_COLUMNS),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    feature_rows: list[dict[str, Any]] = []
    numeric_values = numeric.apply(pd.to_numeric, errors="coerce")
    for feature, row in numeric_values.iterrows():
        values = row.to_numpy(dtype=float)
        partial_r2, sse_reduced, sse_full, status = _batch_partial_r2(
            values,
            metadata,
            batch_column=batch_column,
            condition_column=condition_column,
        )
        feature_rows.append({
            "assay": str(assay),
            "stage": str(stage),
            "feature": str(feature),
            "n_valid": int(np.isfinite(values).sum()),
            "batch_partial_r2": partial_r2,
            "sse_condition": sse_reduced,
            "sse_condition_batch": sse_full,
            "status": status,
        })
    feature_qc = pd.DataFrame(feature_rows, columns=_BATCH_FEATURE_COLUMNS)
    assessed = pd.to_numeric(
        feature_qc.loc[feature_qc["status"] == "assessable", "batch_partial_r2"],
        errors="coerce",
    ).dropna()

    pca_scores, pca_variance = _pca_batch_diagnostic(
        numeric_values,
        metadata,
        assay=str(assay),
        stage=str(stage),
        batch_column=batch_column,
        condition_column=condition_column,
        max_components=max_components,
        max_features=max_features,
    )
    pc_r2 = (
        pd.to_numeric(pca_variance.get("batch_partial_r2"), errors="coerce").dropna()
        if not pca_variance.empty else pd.Series(dtype=float)
    )
    if assessed.empty:
        base_summary["status"] = "insufficient_feature_df"
    else:
        p90 = float(assessed.quantile(0.90))
        base_summary.update({
            "status": "assessable",
            "n_features_assessed": int(len(assessed)),
            "feature_batch_partial_r2_median": float(assessed.median()),
            "feature_batch_partial_r2_p90": p90,
            "fraction_features_high_batch": float(
                assessed.ge(float(high_batch_partial_r2)).mean()
            ),
            "pca_batch_partial_r2_max": (
                float(pc_r2.max()) if not pc_r2.empty else float("nan")
            ),
            "suggested_batch_reliability": float(np.clip(1.0 - p90, 0.25, 1.0)),
        })
    return ResidualBatchDiagnostic(
        summary=base_summary,
        feature_qc=feature_qc,
        pca_scores=pca_scores,
        pca_variance=pca_variance,
    )


def _write_batch_pca_plots(
    scores: pd.DataFrame,
    variance: pd.DataFrame,
    output_dir: Path,
    *,
    batch_column: str,
    condition_column: str,
) -> list[Path]:
    """按 assay×stage 输出 batch 颜色/condition 形状的 PCA 审计图。"""
    if scores.empty or "PC1" not in scores.columns:
        return []
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib 不可用，跳过 batch PCA 图片；CSV 仍已输出")
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for (assay, stage), group in scores.groupby(["assay", "stage"], sort=True):
        y_column = "PC2" if "PC2" in group.columns else None
        batches = sorted(group[batch_column].astype(str).unique())
        conditions = sorted(group[condition_column].astype(str).unique())
        cmap = plt.get_cmap("tab10")
        batch_colors = {value: cmap(index % 10) for index, value in enumerate(batches)}
        markers = ["o", "s", "^", "D", "P", "X", "v", "<", ">"]
        condition_markers = {
            value: markers[index % len(markers)]
            for index, value in enumerate(conditions)
        }
        fig, axis = plt.subplots(figsize=(7.2, 5.5))
        for _, row in group.iterrows():
            axis.scatter(
                float(row["PC1"]),
                float(row[y_column]) if y_column else 0.0,
                color=batch_colors[str(row[batch_column])],
                marker=condition_markers[str(row[condition_column])],
                edgecolor="black",
                linewidth=0.5,
                s=70,
            )
            axis.annotate(
                str(row["sample_id"]),
                (float(row["PC1"]), float(row[y_column]) if y_column else 0.0),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )
        one_variance = variance[
            (variance["assay"].astype(str) == str(assay))
            & (variance["stage"].astype(str) == str(stage))
        ]
        variance_map = {
            str(row["component"]): float(row["explained_variance_ratio"])
            for _, row in one_variance.iterrows()
            if pd.notna(row.get("explained_variance_ratio"))
        }
        axis.set_xlabel(f"PC1 ({variance_map.get('PC1', np.nan):.1%})")
        axis.set_ylabel(
            f"PC2 ({variance_map.get('PC2', np.nan):.1%})" if y_column else "PC2 unavailable"
        )
        axis.set_title(f"{assay} — {stage} residual batch PCA")
        batch_handles = [
            plt.Line2D(
                [], [], marker="o", linestyle="", color=batch_colors[value],
                label=f"batch={value}", markeredgecolor="black",
            )
            for value in batches
        ]
        condition_handles = [
            plt.Line2D(
                [], [], marker=condition_markers[value], linestyle="", color="gray",
                label=f"condition={value}", markeredgecolor="black",
            )
            for value in conditions
        ]
        axis.legend(
            handles=batch_handles + condition_handles,
            fontsize=8,
            frameon=False,
            loc="best",
        )
        axis.axhline(0.0, color="#cccccc", linewidth=0.6)
        axis.axvline(0.0, color="#cccccc", linewidth=0.6)
        fig.tight_layout()
        path = output_dir / f"{assay}_{stage}_batch_pca.png"
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths


def _write_harmonization_report(
    output_dir: Path,
    *,
    config_file: Path,
    manifest: pd.DataFrame,
    missingness: pd.DataFrame,
    batch_qc: pd.DataFrame,
    pca_paths: list[Path],
) -> Path:
    lines = [
        "# Multi-Omics Harmonization Report",
        "",
        f"- Config: `{config_file}`",
        f"- Assays: {len(manifest)}",
        f"- Output: `{output_dir}`",
        "",
        "## Sample Alignment and Missingness",
        "",
    ]
    if missingness.empty:
        lines.append("No assay remained for missingness assessment.")
    else:
        lines.extend([
            "| Assay | Features | Samples | Missing fraction |",
            "|---|---:|---:|---:|",
        ])
        for _, row in missingness.iterrows():
            lines.append(
                f"| {row['assay']} | {int(row['n_features'])} | {int(row['n_samples'])} "
                f"| {float(row['missing_fraction']):.3f} |"
            )
    lines.extend([
        "",
        "## Residual Batch Diagnostics",
        "",
        "`batch_partial_r2` is the incremental batch variance after condition is retained. "
        "The suggested reliability is an audit aid and is **not automatically applied** "
        "to EvidenceProfile scores.",
        "",
    ])
    if batch_qc.empty:
        lines.append("Batch diagnostics were disabled or no mapped samples were available.")
    else:
        lines.extend([
            "| Assay | Stage | n | Batches | Feature partial R² p90 | PCA partial R² max | Suggested reliability | Status |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        ])
        for _, row in batch_qc.iterrows():
            def _format(value: Any) -> str:
                numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
                return f"{float(numeric):.3f}" if pd.notna(numeric) else "NA"

            lines.append(
                f"| {row['assay']} | {row['stage']} | {int(row['n_samples'])} "
                f"| {int(row['n_batches'])} | {_format(row['feature_batch_partial_r2_p90'])} "
                f"| {_format(row['pca_batch_partial_r2_max'])} "
                f"| {_format(row['suggested_batch_reliability'])} | {row['status']} |"
            )
    if pca_paths:
        lines.extend(["", "## PCA Audit Plots", ""])
        for path in pca_paths:
            relative = path.relative_to(output_dir).as_posix()
            lines.append(f"- [{path.stem}]({relative})")
    lines.extend([
        "",
        "## Interpretation Boundary",
        "",
        "- A corrected matrix is for sample-level exploratory integration; differential "
        "inference should remain in each assay's native statistical model.",
        "- A confounded condition×batch design must be analyzed within batch and combined "
        "using estimand-compatible Meta-analysis.",
        "- Single-batch data cannot demonstrate cross-batch reproducibility.",
        "",
    ])
    report_path = output_dir / "harmonization_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def normalize_assay(matrix: pd.DataFrame, assay: str) -> pd.DataFrame:
    """执行适合下游协调的组学内归一化，不负责推断模型。"""
    values = matrix.apply(pd.to_numeric, errors="coerce").astype(float)
    assay_name = str(assay).lower()
    if assay_name in {"proteomics", "protein", "metabolomics", "metabolite"}:
        positive = values.where(values > 0)
        logged = np.log2(positive)
        medians = logged.median(axis=0, skipna=True)
        target = float(medians.median(skipna=True)) if medians.notna().any() else 0.0
        return logged.subtract(medians - target, axis=1)
    if assay_name in {"transcriptomics", "rna", "rnaseq"}:
        library_sizes = values.sum(axis=0, skipna=True)
        library_sizes = library_sizes.replace(0, np.nan)
        counts_per_million = values.divide(library_sizes, axis=1) * 1_000_000.0
        return np.log2(counts_per_million + 1.0)
    if assay_name in {"mfa", "fluxomics", "flux"}:
        means = values.mean(axis=1, skipna=True)
        stds = values.std(axis=1, skipna=True, ddof=1).replace(0, np.nan)
        return values.subtract(means, axis=0).divide(stds, axis=0)
    raise ValueError(f"Unsupported assay normalization: {assay}")


def align_assay_samples(
    assays: Mapping[str, pd.DataFrame],
    sample_map: pd.DataFrame,
    *,
    min_assays: int = 1,
) -> HarmonizedAssays:
    """用显式 sample map 将各 assay 列映射到 primary sample，并记录缺失。"""
    missing_columns = REQUIRED_SAMPLE_MAP_COLUMNS - set(sample_map.columns)
    if missing_columns:
        raise ValueError("sample_map missing columns: " + ", ".join(sorted(missing_columns)))
    if min_assays < 1:
        raise ValueError("min_assays must be >= 1")

    mapping = sample_map.copy()
    missing_values = mapping[list(REQUIRED_SAMPLE_MAP_COLUMNS)].isna()
    empty_values = mapping[list(REQUIRED_SAMPLE_MAP_COLUMNS)].apply(
        lambda column: column.astype(str).str.strip().eq("")
    )
    invalid_rows = missing_values.any(axis=1) | empty_values.any(axis=1)
    if invalid_rows.any():
        invalid_columns = sorted(
            column
            for column in REQUIRED_SAMPLE_MAP_COLUMNS
            if bool(missing_values[column].any() or empty_values[column].any())
        )
        raise ValueError(
            "sample_map has missing/empty values in required columns: "
            + ", ".join(invalid_columns)
        )
    for column in REQUIRED_SAMPLE_MAP_COLUMNS:
        mapping[column] = mapping[column].astype(str).str.strip()
    duplicates = mapping.duplicated(["assay", "assay_sample_id"], keep=False)
    if duplicates.any():
        raise ValueError("Each assay_sample_id must map exactly once within an assay")
    metadata_conflicts = (
        mapping.groupby("primary_sample_id")["condition"]
        .nunique(dropna=False)
        .gt(1)
    )
    if metadata_conflicts.any():
        conflicting = ", ".join(metadata_conflicts[metadata_conflicts].index.astype(str))
        raise ValueError(
            "Primary samples have inconsistent condition across assays: "
            + conflicting
        )

    aligned: dict[str, pd.DataFrame] = {}
    assay_membership: dict[str, set[str]] = {}
    for assay, matrix in assays.items():
        one_map = mapping[mapping["assay"] == str(assay)]
        if one_map.empty:
            continue
        present = one_map[one_map["assay_sample_id"].isin(matrix.columns)].copy()
        if present.empty:
            continue
        subset = matrix[present["assay_sample_id"].tolist()].copy()
        rename = dict(zip(present["assay_sample_id"], present["primary_sample_id"], strict=True))
        subset = subset.rename(columns=rename)
        if subset.columns.duplicated().any():
            raise ValueError(f"Assay '{assay}' maps multiple columns to one primary sample")
        aligned[str(assay)] = subset
        assay_membership[str(assay)] = set(subset.columns.astype(str))

    sample_counts: dict[str, int] = {}
    for members in assay_membership.values():
        for sample in members:
            sample_counts[sample] = sample_counts.get(sample, 0) + 1
    selected = tuple(sorted(sample for sample, count in sample_counts.items() if count >= min_assays))
    aligned = {
        assay: matrix.reindex(columns=list(selected)).copy()
        for assay, matrix in aligned.items()
    }

    missing_rows: list[dict[str, object]] = []
    for assay, subset in aligned.items():
        missing_rows.append({
            "assay": assay,
            "n_features": int(subset.shape[0]),
            "n_samples": int(subset.shape[1]),
            "n_missing": int(subset.isna().sum().sum()),
            "missing_fraction": float(subset.isna().to_numpy().mean()) if subset.size else 0.0,
            "n_complete_samples": int((subset.isna().sum(axis=0) == 0).sum()),
            "n_complete_features": int((subset.isna().sum(axis=1) == 0).sum()),
        })

    selected_mapping = mapping[mapping["primary_sample_id"].isin(selected)].copy()
    condition = (
        selected_mapping.drop_duplicates("primary_sample_id")
        .set_index("primary_sample_id")["condition"]
    )
    batch_wide = selected_mapping.pivot_table(
        index="primary_sample_id",
        columns="assay",
        values="batch_id",
        aggfunc="first",
    )
    batch_wide.columns = [f"batch_id__{column}" for column in batch_wide.columns]
    metadata = condition.to_frame().join(batch_wide, how="left").reindex(selected)
    batch_columns = [column for column in metadata.columns if column.startswith("batch_id__")]
    metadata["batch_id"] = metadata[batch_columns].apply(
        lambda row: (
            next(iter(set(row.dropna().astype(str))))
            if len(set(row.dropna().astype(str))) == 1
            else np.nan
        ),
        axis=1,
    )
    return HarmonizedAssays(
        assays=aligned,
        sample_metadata=metadata,
        sample_map=mapping,
        missingness=pd.DataFrame(missing_rows),
        primary_samples=selected,
    )


def correct_batch_preserving_condition(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    batch_column: str = "batch_id",
    condition_column: str = "condition",
) -> pd.DataFrame:
    """用线性模型移除 batch 项，同时保留 condition 拟合项。

    若 batch 与 condition 完全混杂，模型不可辨识，必须分批建模后做 Meta。
    """
    required = {batch_column, condition_column}
    missing = required - set(sample_metadata.columns)
    if missing:
        raise ValueError("sample_metadata missing columns: " + ", ".join(sorted(missing)))
    samples = [str(column) for column in matrix.columns]
    if not set(samples).issubset(set(sample_metadata.index.astype(str))):
        raise ValueError("sample_metadata does not cover all assay columns")
    metadata = sample_metadata.copy()
    metadata.index = metadata.index.astype(str)
    metadata = metadata.loc[samples]

    cross = pd.crosstab(metadata[condition_column], metadata[batch_column])
    if (cross.gt(0).sum(axis=1) == 1).all() and (cross.gt(0).sum(axis=0) == 1).all():
        raise ValueError(
            "condition and batch are completely confounded; use within-batch effects and Meta-analysis"
        )

    condition_design = pd.get_dummies(
        metadata[condition_column].astype(str), prefix="condition", drop_first=True, dtype=float
    )
    batch_design = pd.get_dummies(
        metadata[batch_column].astype(str), prefix="batch", drop_first=True, dtype=float
    )
    intercept = pd.DataFrame({"intercept": np.ones(len(metadata))}, index=metadata.index)
    design = pd.concat([intercept, condition_design, batch_design], axis=1)
    if np.linalg.matrix_rank(design.to_numpy(dtype=float)) < design.shape[1]:
        raise ValueError(
            "condition and batch design is rank deficient/confounded; use within-batch Meta-analysis"
        )

    y_values = matrix.loc[:, samples].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).T
    x_values = design.to_numpy(dtype=float)
    corrected = y_values.copy()
    batch_columns = [index for index, name in enumerate(design.columns) if name.startswith("batch_")]
    for feature_index in range(y_values.shape[1]):
        y_feature = y_values[:, feature_index]
        valid = np.isfinite(y_feature)
        if valid.sum() < design.shape[1]:
            continue
        coefficients, *_ = np.linalg.lstsq(x_values[valid], y_feature[valid], rcond=None)
        if batch_columns:
            corrected[:, feature_index] = y_feature - x_values[:, batch_columns] @ coefficients[batch_columns]
    return pd.DataFrame(corrected.T, index=matrix.index, columns=samples)


def _resolve_config_path(value: str | Path, config_dir: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_dir / path).resolve()


def _read_assay_matrix(path: Path, feature_column: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Assay matrix not found: {path}")
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        table = pd.read_excel(path)
    elif suffix in {".tsv", ".txt"}:
        table = pd.read_csv(path, sep="\t")
    else:
        table = pd.read_csv(path)
    if feature_column not in table.columns:
        raise ValueError(f"Assay matrix {path} missing feature column '{feature_column}'")
    if table[feature_column].isna().any() or table[feature_column].duplicated().any():
        raise ValueError(f"Assay matrix {path} requires unique non-missing feature IDs")
    return table.set_index(feature_column)


def run_harmonization(config_path: str | Path) -> HarmonizedAssays:
    """从 YAML 配置执行样本级协调并导出可审计矩阵。

    该入口用于潜变量/相似网络等样本级补充分析；差异推断仍应在各组学
    原生模型中完成。若 batch 与 condition 混杂，将拒绝校正并要求使用
    批次内效应量与随机效应 Meta。
    """
    config_file = Path(config_path).resolve()
    if not config_file.exists():
        raise FileNotFoundError(f"Harmonization config not found: {config_file}")
    config = yaml.safe_load(config_file.read_text(encoding="utf-8")) or {}
    if not isinstance(config, dict):
        raise ValueError("Harmonization config root must be a mapping")
    config_dir = config_file.parent

    output_value = config.get("output_dir")
    sample_map_value = config.get("sample_map")
    assay_specs = config.get("assays", {}) or {}
    if not output_value or not sample_map_value or not assay_specs:
        raise ValueError("Harmonization config requires output_dir, sample_map, and assays")
    output_dir = _resolve_config_path(output_value, config_dir)
    sample_map_path = _resolve_config_path(sample_map_value, config_dir)
    if not sample_map_path.exists():
        raise FileNotFoundError(f"Sample map not found: {sample_map_path}")
    sample_map = pd.read_csv(sample_map_path)

    normalized: dict[str, pd.DataFrame] = {}
    source_paths: dict[str, Path] = {}
    for assay, raw_spec in assay_specs.items():
        if not isinstance(raw_spec, dict) or not raw_spec.get("path"):
            raise ValueError(f"Assay '{assay}' requires path and optional feature_column")
        source_path = _resolve_config_path(raw_spec["path"], config_dir)
        feature_column = str(raw_spec.get("feature_column", "feature"))
        matrix = _read_assay_matrix(source_path, feature_column)
        normalization = str(raw_spec.get("normalization", assay))
        normalized[str(assay)] = normalize_assay(matrix, normalization)
        source_paths[str(assay)] = source_path

    aligned = align_assay_samples(
        normalized,
        sample_map,
        min_assays=int(config.get("min_assays", 1)),
    )
    if not aligned.assays:
        raise ValueError("No assay samples remained after explicit sample-map alignment")

    correction = dict(config.get("batch_correction", {}) or {})
    correction_enabled = bool(correction.get("enabled", False))
    batch_qc_config = dict(config.get("batch_qc", {}) or {})
    batch_qc_enabled = bool(batch_qc_config.get("enabled", True))
    batch_column = str(correction.get("batch_column", "batch_id"))
    condition_column = str(correction.get("condition_column", "condition"))
    output_assays: dict[str, pd.DataFrame] = {}
    manifest_rows: list[dict[str, object]] = []
    batch_qc_rows: list[dict[str, Any]] = []
    batch_feature_tables: list[pd.DataFrame] = []
    batch_pca_tables: list[pd.DataFrame] = []
    batch_pca_variance_tables: list[pd.DataFrame] = []
    for assay, matrix in aligned.assays.items():
        corrected = matrix
        applied = False
        assay_mapping = aligned.sample_map[
            (aligned.sample_map["assay"] == assay)
            & aligned.sample_map["primary_sample_id"].isin(aligned.primary_samples)
        ].drop_duplicates("primary_sample_id")
        assay_metadata = assay_mapping.set_index("primary_sample_id")[[
            "condition", "batch_id",
        ]]
        assay_metadata.index = assay_metadata.index.astype(str)
        present_samples = [
            sample for sample in aligned.primary_samples
            if sample in assay_metadata.index and sample in matrix.columns
        ]

        before_qc: ResidualBatchDiagnostic | None = None
        after_qc: ResidualBatchDiagnostic | None = None
        if batch_qc_enabled and present_samples:
            before_qc = diagnose_residual_batch(
                matrix[present_samples],
                assay_metadata,
                assay=str(assay),
                stage="normalized",
                batch_column=batch_column,
                condition_column=condition_column,
                max_components=int(batch_qc_config.get("max_components", 5)),
                max_features=int(batch_qc_config.get("max_features", 5000)),
                high_batch_partial_r2=float(
                    batch_qc_config.get("high_batch_partial_r2", 0.20)
                ),
            )
            batch_qc_rows.append(before_qc.summary)
            batch_feature_tables.append(before_qc.feature_qc)
            batch_pca_tables.append(before_qc.pca_scores)
            batch_pca_variance_tables.append(before_qc.pca_variance)

        has_multiple_batches = (
            bool(present_samples)
            and assay_metadata.loc[present_samples, batch_column].nunique(dropna=True) >= 2
        )
        if correction_enabled and has_multiple_batches:
            corrected_present = correct_batch_preserving_condition(
                matrix[present_samples],
                assay_metadata,
                batch_column=batch_column,
                condition_column=condition_column,
            )
            corrected = corrected_present.reindex(columns=list(aligned.primary_samples))
            applied = True
            if batch_qc_enabled:
                after_qc = diagnose_residual_batch(
                    corrected_present,
                    assay_metadata,
                    assay=str(assay),
                    stage="corrected",
                    batch_column=batch_column,
                    condition_column=condition_column,
                    max_components=int(batch_qc_config.get("max_components", 5)),
                    max_features=int(batch_qc_config.get("max_features", 5000)),
                    high_batch_partial_r2=float(
                        batch_qc_config.get("high_batch_partial_r2", 0.20)
                    ),
                )
                batch_qc_rows.append(after_qc.summary)
                batch_feature_tables.append(after_qc.feature_qc)
                batch_pca_tables.append(after_qc.pca_scores)
                batch_pca_variance_tables.append(after_qc.pca_variance)
        output_assays[assay] = corrected
        manifest_rows.append({
            "assay": assay,
            "source_path": str(source_paths[assay]),
            "n_features": int(corrected.shape[0]),
            "n_samples": int(corrected.shape[1]),
            "normalization": str(assay_specs[assay].get("normalization", assay)),
            "batch_corrected": applied,
            "batch_qc_status_before": (
                before_qc.summary.get("status", "not_run") if before_qc else "not_run"
            ),
            "batch_qc_status_after": (
                after_qc.summary.get("status", "not_run") if after_qc else "not_run"
            ),
            "suggested_batch_reliability_before": (
                before_qc.summary.get("suggested_batch_reliability", np.nan)
                if before_qc else np.nan
            ),
            "suggested_batch_reliability_after": (
                after_qc.summary.get("suggested_batch_reliability", np.nan)
                if after_qc else np.nan
            ),
            "sample_map_path": str(sample_map_path),
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    for assay, matrix in output_assays.items():
        matrix.to_csv(output_dir / f"{assay}_harmonized.csv", encoding="utf-8-sig")
    aligned.sample_metadata.to_csv(
        output_dir / "sample_metadata.csv", encoding="utf-8-sig"
    )
    aligned.sample_map.to_csv(
        output_dir / "resolved_sample_map.csv", index=False, encoding="utf-8-sig"
    )
    aligned.missingness.to_csv(
        output_dir / "missingness_qc.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(manifest_rows).to_csv(
        output_dir / "harmonization_manifest.csv", index=False, encoding="utf-8-sig"
    )
    batch_qc_columns = [
        "assay", "stage", "n_samples", "n_features", "n_conditions",
        "n_batches", "n_features_assessed", "feature_batch_partial_r2_median",
        "feature_batch_partial_r2_p90", "fraction_features_high_batch",
        "pca_batch_partial_r2_max", "suggested_batch_reliability", "status",
    ]
    batch_qc_table = pd.DataFrame(batch_qc_rows, columns=batch_qc_columns)
    batch_qc_table.to_csv(
        output_dir / "batch_qc.csv", index=False, encoding="utf-8-sig"
    )
    nonempty_feature_tables = [
        table for table in batch_feature_tables if not table.empty
    ]
    batch_feature_qc = (
        pd.concat(nonempty_feature_tables, ignore_index=True)
        if nonempty_feature_tables else pd.DataFrame(columns=_BATCH_FEATURE_COLUMNS)
    )
    batch_feature_qc.to_csv(
        output_dir / "batch_feature_qc.csv", index=False, encoding="utf-8-sig"
    )
    nonempty_pca_tables = [table for table in batch_pca_tables if not table.empty]
    batch_pca_scores = (
        pd.concat(nonempty_pca_tables, ignore_index=True)
        if nonempty_pca_tables else pd.DataFrame()
    )
    batch_pca_scores.to_csv(
        output_dir / "batch_pca_scores.csv", index=False, encoding="utf-8-sig"
    )
    nonempty_variance_tables = [
        table for table in batch_pca_variance_tables if not table.empty
    ]
    batch_pca_variance = (
        pd.concat(nonempty_variance_tables, ignore_index=True)
        if nonempty_variance_tables else pd.DataFrame()
    )
    batch_pca_variance.to_csv(
        output_dir / "batch_pca_variance.csv", index=False, encoding="utf-8-sig"
    )
    pca_paths = _write_batch_pca_plots(
        batch_pca_scores,
        batch_pca_variance,
        output_dir / "batch_pca",
        batch_column=batch_column,
        condition_column=condition_column,
    )
    _write_harmonization_report(
        output_dir,
        config_file=config_file,
        manifest=pd.DataFrame(manifest_rows),
        missingness=aligned.missingness,
        batch_qc=batch_qc_table,
        pca_paths=pca_paths,
    )

    return HarmonizedAssays(
        assays=output_assays,
        sample_metadata=aligned.sample_metadata,
        sample_map=aligned.sample_map,
        missingness=aligned.missingness,
        primary_samples=aligned.primary_samples,
    )
