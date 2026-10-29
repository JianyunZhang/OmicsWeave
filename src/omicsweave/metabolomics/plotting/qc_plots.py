"""
PCA 和 PLS-DA 聚类可视化。
"""
from __future__ import annotations

import os
import warnings
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# 确保中文字体设置生效
from . import _config  # noqa: F401


@dataclass
class PLSDAResult:
    scores: pd.DataFrame
    vip: pd.Series
    r2y: float
    q2: float
    cv_accuracy: float
    permutation_pvalue: float
    status: str
    classes: tuple[str, ...]
    n_components: int
    n_permutations: int


def run_pca(
    df,
    sample_cols: list[str],
    tissue: str,
    output_dir: str,
    groups: dict[str, str],
    group_colors: dict[str, str] | None = None,
    *,
    qc_cols: list[str] | None = None,
):
    """PCA 散点图；若存在 pooled QC，则一并显示以评估技术聚类。"""
    qc_cols = list(qc_cols or [])
    plotted_cols = sample_cols + [col for col in qc_cols if col not in sample_cols]
    if len(plotted_cols) < 2:
        print("    PCA 跳过: 样本数不足 2")
        return None
    X = df[plotted_cols].T.astype(float)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    component_count = min(2, X_scaled.shape[0], X_scaled.shape[1])
    pca = PCA(n_components=component_count)
    components = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame(components, columns=[f"PC{i + 1}" for i in range(component_count)])
    if "PC2" not in pca_df:
        pca_df["PC2"] = 0.0
    pca_df["Sample"] = X.index
    plot_groups = {**groups, **{col: "QC" for col in qc_cols}}
    pca_df["Group"] = pca_df["Sample"].map(plot_groups)

    # 保存 PCA 坐标
    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)
    pca_df.to_csv(os.path.join(out_path, "PCA_data.csv"), index=False, encoding='utf-8-sig')

    # 配色
    unique_groups = sorted(set(plot_groups.values()))
    if group_colors:
        palette = {g: group_colors.get(g, "#888888") for g in unique_groups}
        palette["QC"] = "#111111"
    else:
        palette = dict(zip(unique_groups, sns.color_palette("viridis", len(unique_groups))))

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Group", s=100, palette=palette, ax=ax)

    for i, txt in enumerate(pca_df["Sample"]):
        label = txt.replace("-", "_").split("_")[0]
        ax.annotate(label, (pca_df.PC1[i], pca_df.PC2[i]), xytext=(5, 5), textcoords="offset points")

    ax.set_title(f"PCA - 代谢组学 - {tissue}")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%})")
    pc2_variance = pca.explained_variance_ratio_[1] if component_count > 1 else 0.0
    ax.set_ylabel(f"PC2 ({pc2_variance:.2%})")

    plt.savefig(os.path.join(out_path, "PCA.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    PCA 保存: {out_path}/PCA.png")
    return pca_df


def plot_preprocessing_qc(
    qc_metrics: pd.DataFrame,
    sample_totals: pd.DataFrame,
    tissue: str,
    output_dir: str,
    *,
    qc_rsd_cutoff: float = 0.30,
) -> None:
    """保存 QC RSD 分布和归一化前后样本总信号图。"""
    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)

    if "QC_RSD" in qc_metrics.columns:
        rsd = pd.to_numeric(qc_metrics["QC_RSD"], errors="coerce").dropna() * 100.0
        if not rsd.empty:
            fig, ax = plt.subplots(figsize=(8, 5))
            sns.histplot(rsd, bins=40, ax=ax, color="#4C72B0")
            ax.axvline(qc_rsd_cutoff * 100.0, color="#D32F2F", linestyle="--", label="RSD cutoff")
            ax.set_xlabel("QC RSD (%)")
            ax.set_ylabel("Metabolites")
            ax.set_title(f"QC feature RSD - {tissue}")
            ax.legend()
            fig.savefig(os.path.join(out_path, "QC_RSD_distribution.png"), dpi=150, bbox_inches="tight")
            plt.close(fig)

    required = {"total_before", "total_after"}
    if required.issubset(sample_totals.columns) and not sample_totals.empty:
        totals = sample_totals.copy()
        totals.index = totals.index.astype(str)
        fig, axes = plt.subplots(2, 1, figsize=(max(8, len(totals) * 0.6), 8), sharex=True)
        axes[0].bar(totals.index, totals["total_before"], color="#999999")
        axes[0].set_ylabel("Raw total signal")
        axes[0].set_title(f"Sample loading before/after normalization - {tissue}")
        axes[1].bar(totals.index, totals["total_after"], color="#4C72B0")
        axes[1].set_ylabel("Normalized total signal")
        axes[1].tick_params(axis="x", rotation=60)
        fig.tight_layout()
        fig.savefig(os.path.join(out_path, "sample_total_signal.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)


def calculate_vip(model) -> np.ndarray:
    """计算 PLS-DA VIP 分数。"""
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_
    n_features, n_components = w.shape
    explained = np.sum(t ** 2, axis=0) * np.sum(q ** 2, axis=0)
    total_explained = float(np.sum(explained))
    if not np.isfinite(total_explained) or total_explained <= 0:
        return np.full(n_features, np.nan)
    weight_norm = np.sum(w ** 2, axis=0)
    weight_norm = np.where(weight_norm > 0, weight_norm, np.nan)
    contribution = (w ** 2 / weight_norm) @ explained
    return np.sqrt(n_features * contribution / total_explained)


def _one_hot(labels: np.ndarray, classes: tuple[str, ...]) -> np.ndarray:
    class_index = {label: index for index, label in enumerate(classes)}
    encoded = np.zeros((len(labels), len(classes)), dtype=float)
    for row, label in enumerate(labels):
        encoded[row, class_index[str(label)]] = 1.0
    return encoded


def _fit_pls_model(X: np.ndarray, y: np.ndarray, n_components: int) -> PLSRegression:
    model = PLSRegression(n_components=n_components, scale=True, max_iter=1000)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.fit(X, y)
    return model


def _q2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    press = float(np.sum((y_true - y_pred) ** 2))
    total = float(np.sum((y_true - y_true.mean(axis=0, keepdims=True)) ** 2))
    return float(1.0 - press / total) if total > 0 else np.nan


def _cross_validated_predictions(
    X: np.ndarray,
    labels: np.ndarray,
    classes: tuple[str, ...],
    n_components: int,
    random_state: int,
) -> tuple[np.ndarray, float, float]:
    y = _one_hot(labels, classes)
    min_class_size = int(pd.Series(labels).value_counts().min())
    if min_class_size < 2:
        return np.full_like(y, np.nan), np.nan, np.nan
    n_splits = min(5, min_class_size)
    splitter = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=random_state,
    )
    predicted = np.full_like(y, np.nan)
    for train, test in splitter.split(X, labels):
        fold_components = min(
            n_components,
            X.shape[1],
            max(1, len(train) - 1),
            max(1, len(classes) - 1),
        )
        try:
            model = _fit_pls_model(X[train], y[train], fold_components)
            predicted[test] = model.predict(X[test])
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            return np.full_like(y, np.nan), np.nan, np.nan
    q2_value = _q2(y, predicted)
    accuracy = float(np.mean(np.argmax(predicted, axis=1) == np.argmax(y, axis=1)))
    return predicted, q2_value, accuracy


def fit_plsda(
    X: pd.DataFrame,
    labels: list[str] | np.ndarray | pd.Series,
    *,
    n_components: int = 2,
    n_permutations: int = 200,
    random_state: int = 42,
) -> PLSDAResult:
    """拟合 one-hot PLS-DA，并用分层交叉验证与标签置换评估模型。"""
    if not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(X)
    full_matrix = X.to_numpy(dtype=float)
    if full_matrix.ndim != 2 or full_matrix.shape[0] < 2 or full_matrix.shape[1] < 1:
        raise ValueError("PLS-DA 至少需要 2 个样本和 1 个代谢物")
    if not np.isfinite(full_matrix).all():
        raise ValueError("PLS-DA 输入包含 NaN 或无穷值")
    label_array = np.asarray(labels, dtype=str)
    if len(label_array) != len(X):
        raise ValueError("PLS-DA 标签数与样本数不一致")
    classes = tuple(sorted(np.unique(label_array).tolist()))
    if len(classes) < 2:
        raise ValueError("PLS-DA 至少需要两个组别")
    if n_components < 1 or n_permutations < 0:
        raise ValueError("PLS-DA 组件数须 >=1，置换次数须 >=0")

    variable_mask = np.nanstd(full_matrix, axis=0, ddof=0) > np.finfo(float).eps
    if not variable_mask.any():
        scores = pd.DataFrame(
            {"Component1": np.zeros(len(X)), "Group": label_array},
            index=X.index,
        )
        return PLSDAResult(
            scores=scores,
            vip=pd.Series(np.nan, index=X.columns, name="VIP", dtype=float),
            r2y=np.nan,
            q2=np.nan,
            cv_accuracy=np.nan,
            permutation_pvalue=np.nan,
            status="degenerate_features",
            classes=classes,
            n_components=0,
            n_permutations=0,
        )
    matrix = full_matrix[:, variable_mask]

    used_components = min(
        n_components,
        matrix.shape[1],
        matrix.shape[0] - 1,
        len(classes) - 1,
    )
    y = _one_hot(label_array, classes)
    model = _fit_pls_model(matrix, y, used_components)
    predicted = model.predict(matrix)
    scores_array = model.x_scores_[:, :used_components]
    score_columns = [f"Component{i + 1}" for i in range(used_components)]
    scores = pd.DataFrame(scores_array, index=X.index, columns=score_columns)
    scores["Group"] = label_array
    vip = pd.Series(np.nan, index=X.columns, name="VIP", dtype=float)
    vip.iloc[np.flatnonzero(variable_mask)] = calculate_vip(model)
    r2y = _q2(y, predicted)

    _, q2_value, cv_accuracy = _cross_validated_predictions(
        matrix,
        label_array,
        classes,
        used_components,
        random_state,
    )
    min_class_size = int(pd.Series(label_array).value_counts().min())
    if min_class_size < 2:
        permutation_pvalue = np.nan
        status = "insufficient_validation"
        completed_permutations = 0
    elif not np.isfinite(q2_value) or not np.isfinite(cv_accuracy):
        permutation_pvalue = np.nan
        status = "validation_failed"
        completed_permutations = 0
    else:
        permutation_q2: list[float] = []
        rng = np.random.default_rng(random_state)
        for _ in range(n_permutations):
            permuted = rng.permutation(label_array)
            _, permuted_q2, _ = _cross_validated_predictions(
                matrix,
                permuted,
                classes,
                used_components,
                random_state,
            )
            if np.isfinite(permuted_q2):
                permutation_q2.append(float(permuted_q2))
        completed_permutations = len(permutation_q2)
        if completed_permutations:
            permutation_pvalue = (
                1.0 + sum(value >= q2_value for value in permutation_q2)
            ) / (completed_permutations + 1.0)
            status = "validated"
        else:
            permutation_pvalue = np.nan
            status = "cross_validated" if n_permutations == 0 else "validation_failed"

    return PLSDAResult(
        scores=scores,
        vip=vip,
        r2y=float(r2y),
        q2=float(q2_value),
        cv_accuracy=float(cv_accuracy),
        permutation_pvalue=float(permutation_pvalue),
        status=status,
        classes=classes,
        n_components=used_components,
        n_permutations=completed_permutations,
    )


def run_plsda(
    df,
    sample_cols: list[str],
    tissue: str,
    output_dir: str,
    groups: dict[str, str],
    group_colors: dict[str, str] | None = None,
    *,
    comp_name: str | None = None,
    n_permutations: int = 200,
    random_state: int = 42,
) -> PLSDAResult:
    """拟合、验证并保存 PLS-DA 得分、VIP 和模型质量指标。"""
    X = df[sample_cols].T.astype(float)
    y_labels = [groups[col] for col in sample_cols]
    result = fit_plsda(
        X,
        y_labels,
        n_components=2,
        n_permutations=n_permutations,
        random_state=random_state,
    )

    unique_groups = sorted(set(groups.values()))
    if group_colors:
        palette = {g: group_colors.get(g, "#888888") for g in unique_groups}
    else:
        palette = dict(zip(unique_groups, sns.color_palette("viridis", len(unique_groups))))

    fig, ax = plt.subplots(figsize=(8, 6))
    score_x = result.scores["Component1"].to_numpy()
    if "Component2" in result.scores:
        score_y = result.scores["Component2"].to_numpy()
    else:
        score_y = np.zeros(len(result.scores))
    sns.scatterplot(
        x=score_x, y=score_y,
        hue=y_labels, style=y_labels, s=100,
        palette=palette, ax=ax,
    )

    for i in range(len(sample_cols)):
        label = sample_cols[i].replace("-", "_").split("_")[0]
        ax.annotate(label, (score_x[i], score_y[i]), xytext=(5, 5), textcoords="offset points")

    comparison_label = f" - {comp_name.replace('_vs_', ' vs ')}" if comp_name else ""
    validation = (
        f"R2Y={result.r2y:.3f}; Q2={result.q2:.3f}; "
        f"CV accuracy={result.cv_accuracy:.3f}; permutation P={result.permutation_pvalue:.3f}"
        if result.status == "validated"
        else "validation unavailable: at least 2 samples per class required"
    )
    ax.set_title(f"PLS-DA - {tissue}{comparison_label}\n{validation}")
    ax.set_xlabel("Component 1")
    ax.set_ylabel("Component 2")

    out_path = os.path.join(output_dir, tissue)
    os.makedirs(out_path, exist_ok=True)
    suffix = f"_{comp_name}" if comp_name else ""
    image_name = f"PLS-DA{suffix}.png"
    plt.savefig(os.path.join(out_path, image_name), dpi=150, bbox_inches="tight")
    plt.close()
    scores = result.scores.copy()
    scores.insert(0, "Sample", scores.index)
    scores.to_csv(
        os.path.join(out_path, f"PLS-DA{suffix}_scores.csv"),
        index=False,
        encoding="utf-8-sig",
    )
    result.vip.rename("VIP").to_csv(
        os.path.join(out_path, f"PLS-DA{suffix}_vip.csv"),
        index_label="Feature",
        encoding="utf-8-sig",
    )
    pd.DataFrame(
        [
            {
                "comparison": comp_name or "global_multiclass",
                "status": result.status,
                "classes": ";".join(result.classes),
                "n_components": result.n_components,
                "R2Y": result.r2y,
                "Q2": result.q2,
                "CV_Accuracy": result.cv_accuracy,
                "Permutation_P": result.permutation_pvalue,
                "Permutations": result.n_permutations,
                "random_state": random_state,
            }
        ]
    ).to_csv(
        os.path.join(out_path, f"PLS-DA{suffix}_validation.csv"),
        index=False,
        encoding="utf-8-sig",
    )
    print(f"    PLS-DA 保存: {out_path}/{image_name} ({result.status})")
    return result
