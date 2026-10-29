"""蛋白质组差异丰度分析的统计核心。

该模块只接收已经 log2 转换并归一化、但**没有插补**的蛋白 × 样本矩阵。
缺失值按蛋白、按比较保留；每个比较建立独立线性模型，并在可识别的前提下
用 limma 风格经验贝叶斯方法压缩残差方差。
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.special import digamma, gammaln, polygamma
from statsmodels.stats.multitest import multipletests

from .experiment_config import normalize_fdr_method  # type: ignore


COMBINED_CSV = "differential_proteins.csv"
DESIGN_REPORT_CSV = "differential_design_report.csv"


@dataclass(frozen=True)
class ComparisonDesign:
    """一个比较对应的全样本设计矩阵与可识别性结论。"""

    numerator: str
    denominator: str
    name: str
    samples: list[str]
    matrix: np.ndarray
    columns: list[str]
    contrast: np.ndarray
    route: str
    reason: str
    rank: int
    residual_df: int
    n_numerator: int
    n_denominator: int
    dropped_constant_covariates: tuple[str, ...] = ()


def _log_mdigamma(x: np.ndarray | float) -> np.ndarray | float:
    """limma/statmod 的 ``logmdigamma(x) = log(x) - digamma(x)``。"""

    return np.log(x) - digamma(x)


def _trigamma_inverse(value: float) -> float:
    """求解 ``trigamma(y) = value``，迭代式与 limma 相同。"""

    if not np.isfinite(value):
        return 0.0 if value > 0 else np.nan
    if value < 0:
        return np.nan
    if value > 1e7:
        return 1.0 / np.sqrt(value)
    if value < 1e-6:
        return 1.0 / max(value, np.finfo(float).tiny)

    y = 0.5 + 1.0 / value
    for _ in range(51):
        tri = float(polygamma(1, y))
        tetra = float(polygamma(2, y))
        difference = tri * (1.0 - tri / value) / tetra
        y += difference
        if -difference / y < 1e-8:
            break
    return float(y)


def _fit_f_dist_equal_df(variances: np.ndarray, dfs: np.ndarray) -> tuple[float, float]:
    """移植 limma ``fitFDist`` 的等自由度 log-F 矩估计。"""

    x = np.asarray(variances, dtype=float)
    df1 = np.asarray(dfs, dtype=float)
    valid = np.isfinite(x) & np.isfinite(df1) & (x >= 0) & (df1 > 1e-15)
    x = x[valid]
    df1 = df1[valid]
    if x.size < 3:
        positive = x[x > 0]
        return (float(np.median(positive)) if positive.size else 1.0), 0.0

    median_variance = float(np.median(x))
    if median_variance <= 0:
        median_variance = float(np.mean(x[x > 0])) if np.any(x > 0) else 1.0
    x = np.maximum(x, 1e-5 * median_variance)

    half_df = df1 / 2.0
    adjusted_log_variance = np.log(x) + _log_mdigamma(half_df)
    mean_adjusted = float(np.mean(adjusted_log_variance))
    variance_adjusted = float(np.var(adjusted_log_variance, ddof=1))
    variance_adjusted -= float(np.mean(polygamma(1, half_df)))

    if variance_adjusted > 0:
        prior_df = 2.0 * _trigamma_inverse(variance_adjusted)
        prior_variance = float(
            np.exp(mean_adjusted - _log_mdigamma(prior_df / 2.0))
        )
    else:
        prior_df = np.inf
        prior_variance = float(np.mean(x))
    return prior_variance, float(prior_df)


def _fit_f_dist_unequal_df(
    variances: np.ndarray,
    dfs: np.ndarray,
) -> tuple[float, float]:
    """移植 limma >=3.62 ``fitFDistUnequalDF1`` 的 profile MLE 核心。"""

    x = np.asarray(variances, dtype=float)
    df1 = np.asarray(dfs, dtype=float)
    valid = np.isfinite(x) & np.isfinite(df1) & (x >= 0) & (df1 >= 0.01)
    x = x[valid]
    df1 = df1[valid]
    if x.size < 3:
        positive = x[x > 0]
        return (float(np.median(positive)) if positive.size else 1.0), 0.0

    positive = x[x > 0]
    median_variance = float(np.median(positive)) if positive.size else 1.0
    x_positive = np.maximum(x, 1e-12 * median_variance)
    d1 = df1 / 2.0
    adjusted = np.log(x_positive) + _log_mdigamma(d1)
    weights = 1.0 / polygamma(1, d1)
    adjusted_mean = float(np.sum(weights * adjusted) / np.sum(weights))
    d1_times_x = d1 * x_positive

    def minus_twice_log_likelihood(parameter: float) -> float:
        d2 = parameter / (1.0 - parameter)
        d2_s20 = d2 * np.exp(adjusted_mean - _log_mdigamma(d2))
        terms = (
            -(d1 + d2) * np.log1p(d1_times_x / d2_s20)
            - d1 * np.log(d2_s20)
            + gammaln(d1 + d2)
            - gammaln(d2)
        )
        return float(-2.0 * np.sum(terms))

    fit = minimize_scalar(
        minus_twice_log_likelihood,
        bounds=(0.5, 0.9998),
        method="bounded",
        options={"xatol": 1e-10},
    )
    if not fit.success or not np.isfinite(fit.fun):
        return _fit_f_dist_equal_df(x, np.repeat(float(np.median(df1)), len(x)))

    d2 = float(fit.x / (1.0 - fit.x))
    prior_variance = float(np.exp(adjusted_mean - _log_mdigamma(d2)))
    return prior_variance, 2.0 * d2


def estimate_limma_priors(
    variances: np.ndarray,
    dfs: np.ndarray,
) -> tuple[float, float]:
    """估计 limma 的先验方差 ``s0²`` 与先验自由度 ``d0``。

    等自由度时使用官方 ``fitFDist`` 的 log-F 矩估计；缺失值导致自由度
    不等时使用 limma 新版的 profile maximum-likelihood 路径。
    """

    variances = np.asarray(variances, dtype=float)
    dfs = np.asarray(dfs, dtype=float)
    valid = np.isfinite(variances) & np.isfinite(dfs) & (variances >= 0) & (dfs > 0)
    variances = variances[valid]
    dfs = dfs[valid]
    if variances.size < 3:
        positive = variances[variances > 0]
        return (float(np.median(positive)) if positive.size else 1.0), 0.0
    if np.allclose(dfs, dfs[0], rtol=0, atol=1e-12):
        return _fit_f_dist_equal_df(variances, dfs)
    return _fit_f_dist_unequal_df(variances, dfs)


def adjust_pvalues(pvalues: Sequence[float], method: str) -> np.ndarray:
    """只校正有限 p 值，并为不可推断行保留 NaN。"""

    normalized_method = normalize_fdr_method(method)
    raw = np.asarray(pvalues, dtype=float)
    adjusted = np.full(raw.shape, np.nan, dtype=float)
    finite = np.isfinite(raw)
    if finite.any():
        adjusted[finite] = multipletests(raw[finite], method=normalized_method)[1]
    return adjusted


def _prepare_metadata(
    expression: pd.DataFrame,
    group_map: dict[str, str],
    sample_metadata: pd.DataFrame | None,
) -> pd.DataFrame:
    samples = [str(s) for s in expression.columns]
    group_map = {str(sample): str(group) for sample, group in group_map.items()}
    missing_groups = [s for s in samples if group_map.get(s, "Unknown") == "Unknown"]
    if missing_groups:
        raise ValueError(
            "Unmapped proteomics samples: " + ", ".join(missing_groups)
        )

    if sample_metadata is None:
        metadata = pd.DataFrame(
            {"condition": [group_map[s] for s in samples]},
            index=samples,
        )
    else:
        metadata = sample_metadata.copy()
        metadata.index = metadata.index.astype(str)
        absent = [s for s in samples if s not in metadata.index]
        if absent:
            raise ValueError("sample_metadata missing samples: " + ", ".join(absent))
        metadata = metadata.loc[samples].copy()
        if "condition" not in metadata.columns:
            metadata["condition"] = [group_map[s] for s in samples]
        configured = metadata["condition"].astype(str)
        conflicts = [s for s in samples if configured.loc[s] != group_map[s]]
        if conflicts:
            raise ValueError(
                "sample_metadata condition conflicts with group_map for: "
                + ", ".join(conflicts)
            )
    metadata["condition"] = metadata["condition"].astype(str)
    return metadata


def _normalize_expression_axes(expression: pd.DataFrame) -> pd.DataFrame:
    matrix = expression.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    if not matrix.index.is_unique:
        duplicates = matrix.index[matrix.index.duplicated()].unique().tolist()
        raise ValueError("Duplicate protein identifiers: " + ", ".join(duplicates[:10]))
    if not matrix.columns.is_unique:
        duplicates = matrix.columns[matrix.columns.duplicated()].unique().tolist()
        raise ValueError("Duplicate sample names: " + ", ".join(duplicates[:10]))
    return matrix.apply(pd.to_numeric, errors="coerce")


def _encode_covariates(
    metadata: pd.DataFrame,
    covariates: Sequence[str],
) -> tuple[pd.DataFrame, list[str]]:
    encoded: list[pd.DataFrame] = []
    dropped: list[str] = []
    for covariate in covariates:
        if covariate not in metadata.columns:
            raise ValueError(f"Configured covariate '{covariate}' missing from sample metadata")
        series = metadata[covariate]
        if series.isna().any():
            missing = metadata.index[series.isna()].astype(str).tolist()
            raise ValueError(
                f"Covariate '{covariate}' has missing values for: {', '.join(missing)}"
            )
        if series.nunique(dropna=False) < 2:
            dropped.append(covariate)
            continue

        if pd.api.types.is_numeric_dtype(series):
            values = series.astype(float)
            std = float(values.std(ddof=0))
            if not np.isfinite(std) or std <= 0:
                dropped.append(covariate)
                continue
            encoded.append(pd.DataFrame(
                {(f"{covariate}[numeric]"): (values - values.mean()) / std},
                index=metadata.index,
            ))
        else:
            categorical = series.astype(str)
            dummies = pd.get_dummies(
                categorical,
                prefix=covariate,
                prefix_sep="[",
                drop_first=True,
                dtype=float,
            )
            dummies.columns = [f"{column}]" for column in dummies.columns]
            encoded.append(dummies)

    if not encoded:
        return pd.DataFrame(index=metadata.index), dropped
    return pd.concat(encoded, axis=1), dropped


def _build_comparison_design(
    metadata: pd.DataFrame,
    numerator: str,
    denominator: str,
    name: str,
    covariates: Sequence[str],
    min_observations_per_group: int,
) -> ComparisonDesign:
    subset = metadata[metadata["condition"].isin([numerator, denominator])].copy()
    samples = subset.index.astype(str).tolist()
    n_num = int((subset["condition"] == numerator).sum())
    n_den = int((subset["condition"] == denominator).sum())

    condition_design = pd.DataFrame(
        {
            f"condition[{numerator}]": (subset["condition"] == numerator).astype(float),
            f"condition[{denominator}]": (subset["condition"] == denominator).astype(float),
        },
        index=subset.index,
    )
    covariate_design, dropped = _encode_covariates(subset, covariates)
    design_df = pd.concat([condition_design, covariate_design], axis=1)
    matrix = design_df.to_numpy(dtype=float)
    rank = int(np.linalg.matrix_rank(matrix)) if matrix.size else 0
    residual_df = int(matrix.shape[0] - rank)
    contrast = np.zeros(matrix.shape[1], dtype=float)
    if contrast.size >= 2:
        contrast[0] = 1.0
        contrast[1] = -1.0

    if n_num < min_observations_per_group or n_den < min_observations_per_group:
        route = "insufficient_data"
        reason = (
            f"each group requires >= {min_observations_per_group} samples "
            f"(observed {n_num}/{n_den})"
        )
    elif rank < matrix.shape[1]:
        route = "design_not_identifiable"
        reason = (
            f"design matrix is rank deficient ({rank}/{matrix.shape[1]}); "
            "condition is confounded with one or more covariates"
        )
    elif residual_df < 1:
        route = "insufficient_residual_df"
        reason = "design has no residual degrees of freedom"
    else:
        route = "inferential"
        reason = "design identifiable"

    return ComparisonDesign(
        numerator=numerator,
        denominator=denominator,
        name=name,
        samples=samples,
        matrix=matrix,
        columns=design_df.columns.astype(str).tolist(),
        contrast=contrast,
        route=route,
        reason=reason,
        rank=rank,
        residual_df=residual_df,
        n_numerator=n_num,
        n_denominator=n_den,
        dropped_constant_covariates=tuple(dropped),
    )


def _fit_one_comparison(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    design: ComparisonDesign,
    min_observations_per_group: int,
) -> pd.DataFrame:
    """逐蛋白拟合一个比较；暂不计算 moderated t 与 p 值。"""

    sample_positions = [expression.columns.get_loc(s) for s in design.samples]
    values = expression.iloc[:, sample_positions].to_numpy(dtype=float)
    sub_meta = metadata.loc[design.samples]
    numerator_mask = (sub_meta["condition"] == design.numerator).to_numpy()
    denominator_mask = (sub_meta["condition"] == design.denominator).to_numpy()

    records: list[dict[str, Any]] = []
    for protein_index, protein in enumerate(expression.index):
        y = values[protein_index]
        finite = np.isfinite(y)
        n_num = int(np.sum(finite & numerator_mask))
        n_den = int(np.sum(finite & denominator_mask))
        num_values = y[finite & numerator_mask]
        den_values = y[finite & denominator_mask]
        raw_fc = (
            float(np.mean(num_values) - np.mean(den_values))
            if num_values.size and den_values.size
            else np.nan
        )
        record: dict[str, Any] = {
            "Protein": protein,
            "log2FoldChange": raw_fc,
            "n_numerator": n_num,
            "n_denominator": n_den,
            "residual_df": np.nan,
            "standard_error": np.nan,
            "tstat": np.nan,
            "pvalue": np.nan,
            "padj": np.nan,
            "evidence_level": "insufficient",
            "inference_status": design.route,
            "_variance": np.nan,
            "_unscaled_se_sq": np.nan,
        }

        if design.route != "inferential":
            if design.route == "insufficient_data":
                record["inference_status"] = "insufficient_replicates"
            records.append(record)
            continue
        if n_num < min_observations_per_group or n_den < min_observations_per_group:
            record["inference_status"] = "insufficient_observations_per_group"
            records.append(record)
            continue

        observed_design = design.matrix[finite]
        observed_y = y[finite]
        rank = int(np.linalg.matrix_rank(observed_design))
        if rank < observed_design.shape[1]:
            record["inference_status"] = "feature_design_not_identifiable"
            records.append(record)
            continue
        residual_df = int(observed_design.shape[0] - rank)
        if residual_df < 1:
            record["inference_status"] = "insufficient_residual_df"
            records.append(record)
            continue

        coefficients, *_ = np.linalg.lstsq(observed_design, observed_y, rcond=None)
        fitted = observed_design @ coefficients
        residuals = observed_y - fitted
        residual_variance = float(np.sum(residuals ** 2) / residual_df)
        xtx_inverse = np.linalg.pinv(observed_design.T @ observed_design)
        unscaled_se_sq = float(design.contrast @ xtx_inverse @ design.contrast)
        if not np.isfinite(unscaled_se_sq) or unscaled_se_sq <= 0:
            record["inference_status"] = "contrast_not_estimable"
            records.append(record)
            continue

        record.update({
            "log2FoldChange": float(design.contrast @ coefficients),
            "residual_df": float(residual_df),
            "inference_status": "model_fit",
            "_variance": residual_variance,
            "_unscaled_se_sq": unscaled_se_sq,
        })
        records.append(record)

    return pd.DataFrame.from_records(records).set_index("Protein")


def _moderate_comparison(
    fitted: pd.DataFrame,
) -> tuple[pd.DataFrame, float, float]:
    result = fitted.copy()
    valid = (
        np.isfinite(result["_variance"])
        & np.isfinite(result["residual_df"])
        & np.isfinite(result["_unscaled_se_sq"])
        & (result["residual_df"] > 0)
        & (result["_unscaled_se_sq"] > 0)
    )
    n_valid = int(valid.sum())
    if n_valid == 0:
        return result, np.nan, 0.0

    if n_valid < 3:
        # Empirical Bayes has no theoretical advantage with fewer than three
        # variances. Preserve ordinary OLS t inference instead of returning
        # unexplained NaNs.
        observed_variance = result.loc[valid, "_variance"].to_numpy(dtype=float)
        observed_df = result.loc[valid, "residual_df"].to_numpy(dtype=float)
        variance_floor = max(
            float(np.median(observed_variance[observed_variance > 0])) * 1e-12
            if np.any(observed_variance > 0) else np.finfo(float).eps,
            np.finfo(float).tiny,
        )
        observed_variance = np.maximum(observed_variance, variance_floor)
        standard_error = np.sqrt(
            observed_variance
            * result.loc[valid, "_unscaled_se_sq"].to_numpy(dtype=float)
        )
        t_statistic = (
            result.loc[valid, "log2FoldChange"].to_numpy(dtype=float)
            / standard_error
        )
        pvalues = 2.0 * stats.t.sf(np.abs(t_statistic), df=observed_df)
        result.loc[valid, "standard_error"] = standard_error
        result.loc[valid, "tstat"] = t_statistic
        result.loc[valid, "pvalue"] = pvalues
        result.loc[valid, "evidence_level"] = "inferential"
        result.loc[valid, "inference_status"] = "ok_unmoderated_t"
        return result, np.nan, 0.0

    prior_variance, prior_df = estimate_limma_priors(
        result.loc[valid, "_variance"].to_numpy(dtype=float),
        result.loc[valid, "residual_df"].to_numpy(dtype=float),
    )
    observed_variance = result.loc[valid, "_variance"].to_numpy(dtype=float)
    observed_df = result.loc[valid, "residual_df"].to_numpy(dtype=float)
    if np.isinf(prior_df):
        posterior_variance = np.repeat(prior_variance, int(valid.sum()))
        posterior_df = np.repeat(np.inf, int(valid.sum()))
    else:
        posterior_variance = (
            observed_df * observed_variance + prior_df * prior_variance
        ) / (observed_df + prior_df)
        posterior_df = observed_df + prior_df

    variance_floor = max(prior_variance * 1e-12, np.finfo(float).tiny)
    posterior_variance = np.maximum(posterior_variance, variance_floor)
    standard_error = np.sqrt(
        posterior_variance
        * result.loc[valid, "_unscaled_se_sq"].to_numpy(dtype=float)
    )
    t_statistic = result.loc[valid, "log2FoldChange"].to_numpy(dtype=float) / standard_error
    pvalues = 2.0 * stats.t.sf(np.abs(t_statistic), df=posterior_df)

    result.loc[valid, "standard_error"] = standard_error
    result.loc[valid, "tstat"] = t_statistic
    result.loc[valid, "pvalue"] = pvalues
    result.loc[valid, "evidence_level"] = "inferential"
    result.loc[valid, "inference_status"] = "ok"
    return result, float(prior_variance), float(prior_df)


def _comparison_to_output(
    combined: pd.DataFrame,
    comparison: str,
    pval_cutoff: float,
    log2fc_cutoff: float,
    method_label: str,
    fdr_method: str,
    design_formula: str,
) -> pd.DataFrame:
    output = pd.DataFrame({
        "Protein": combined["Protein"],
        "uniprot_id": combined["uniprot_id"],
        "log2FoldChange": combined[f"Log2FC_{comparison}"],
        "standard_error": combined[f"SE_{comparison}"],
        "tstat": combined[f"Tstat_{comparison}"],
        "residual_df": combined[f"ResidualDF_{comparison}"],
        "n_numerator": combined[f"NObsNum_{comparison}"],
        "n_denominator": combined[f"NObsDen_{comparison}"],
        "pvalue": combined[f"Pval_{comparison}"],
        "padj": combined[f"Padj_{comparison}"],
        "evidence_level": combined[f"Evidence_{comparison}"],
        "inference_status": combined[f"InferenceStatus_{comparison}"],
    })
    output["method"] = method_label
    output.loc[
        output["inference_status"] == "ok_unmoderated_t", "method"
    ] = "ols_t_test_unmoderated"
    output["design_formula"] = design_formula
    output["fdr_method"] = fdr_method
    output["significant"] = (
        output["padj"].notna()
        & (output["padj"] < pval_cutoff)
        & (output["log2FoldChange"].abs() > log2fc_cutoff)
        & (output["evidence_level"] == "inferential")
    )
    return output


def run_limma_analysis(
    expression: pd.DataFrame,
    group_map: dict[str, str],
    comparisons: list[tuple[str, str, str]],
    output_dir: str,
    tissue: str,
    pval_cutoff: float,
    log2fc_cutoff: float,
    fdr_method: str,
    *,
    feature_to_uniprot: dict[str, str] | None = None,
    sample_metadata: pd.DataFrame | None = None,
    design_covariates: Sequence[str] = (),
    min_observations_per_group: int = 2,
    flush: Callable[[str], None] = print,
) -> pd.DataFrame:
    """线性模型 + limma 方差压缩，支持分类/连续协变量和蛋白级缺失值。"""

    os.makedirs(output_dir, exist_ok=True)
    expression = _normalize_expression_axes(expression)
    fdr_method = normalize_fdr_method(fdr_method)
    metadata = _prepare_metadata(expression, group_map, sample_metadata)

    proteins = expression.index.astype(str)
    combined = pd.DataFrame({
        "Protein": proteins,
        "uniprot_id": [(feature_to_uniprot or {}).get(p, "") for p in proteins],
    })
    combined.index = proteins

    groups_present = sorted(metadata["condition"].unique())
    for group in groups_present:
        samples = metadata.index[metadata["condition"] == group].tolist()
        combined[f"Mean_{group}"] = expression[samples].mean(axis=1, skipna=True)

    report_rows: list[dict[str, Any]] = []
    designs: dict[str, ComparisonDesign] = {}
    for numerator, denominator, name in comparisons:
        design = _build_comparison_design(
            metadata,
            numerator,
            denominator,
            name,
            design_covariates,
            min_observations_per_group,
        )
        designs[name] = design
        fitted = _fit_one_comparison(
            expression,
            metadata,
            design,
            min_observations_per_group,
        )
        fitted, prior_variance, prior_df = _moderate_comparison(fitted)
        fitted["padj"] = adjust_pvalues(fitted["pvalue"], fdr_method)

        combined[f"Log2FC_{name}"] = fitted["log2FoldChange"]
        combined[f"SE_{name}"] = fitted["standard_error"]
        combined[f"Tstat_{name}"] = fitted["tstat"]
        combined[f"ResidualDF_{name}"] = fitted["residual_df"]
        combined[f"NObsNum_{name}"] = fitted["n_numerator"].astype(int)
        combined[f"NObsDen_{name}"] = fitted["n_denominator"].astype(int)
        combined[f"Pval_{name}"] = fitted["pvalue"]
        combined[f"Padj_{name}"] = fitted["padj"]
        combined[f"Evidence_{name}"] = fitted["evidence_level"]
        combined[f"InferenceStatus_{name}"] = fitted["inference_status"]

        tested = int(np.isfinite(fitted["pvalue"]).sum())
        report_rows.append({
            "tissue": tissue,
            "comparison": name,
            "numerator": numerator,
            "denominator": denominator,
            "n_numerator": design.n_numerator,
            "n_denominator": design.n_denominator,
            "n_samples": len(design.samples),
            "design_formula": "~0 + condition" + "".join(
                f" + {covariate}" for covariate in design_covariates
            ),
            "design_columns": ";".join(design.columns),
            "design_rank": design.rank,
            "n_design_columns": len(design.columns),
            "residual_df_complete_cases": design.residual_df,
            "route": design.route,
            "reason": design.reason,
            "dropped_constant_covariates": ";".join(design.dropped_constant_covariates),
            "n_features_total": len(expression),
            "n_features_tested": tested,
            "n_features_insufficient": int(len(expression) - tested),
            "prior_variance": prior_variance,
            "prior_df": prior_df,
            "fdr_method": fdr_method,
            "differential_input": "log2_normalized_unimputed",
        })
        if design.route == "inferential":
            if prior_df == 0 or not np.isfinite(prior_variance):
                flush(
                    f"  [{name}] OLS t fallback: tested={tested}/{len(expression)} "
                    "(fewer than 3 estimable protein variances)"
                )
            else:
                flush(
                    f"  [{name}] limma: tested={tested}/{len(expression)}, "
                    f"s0_sq={prior_variance:.4g}, d0={prior_df:.3g}"
                )
        else:
            flush(f"  [{name}] [SKIP INFERENCE] {design.reason}")

    combined.reset_index(drop=True, inplace=True)
    combined.to_csv(os.path.join(output_dir, COMBINED_CSV), index=False)
    pd.DataFrame(report_rows).to_csv(
        os.path.join(output_dir, DESIGN_REPORT_CSV), index=False,
    )

    for _, _, name in comparisons:
        design = designs[name]
        per_comparison = _comparison_to_output(
            combined,
            name,
            pval_cutoff,
            log2fc_cutoff,
            "limma_moderated_linear_model",
            fdr_method,
            "~0 + condition" + "".join(
                f" + {covariate}" for covariate in design_covariates
            ),
        )
        per_comparison.to_csv(
            os.path.join(output_dir, f"dp_{tissue}_{name}.csv"), index=False,
        )
    return combined


def run_wilcoxon_analysis(
    expression: pd.DataFrame,
    group_map: dict[str, str],
    comparisons: list[tuple[str, str, str]],
    output_dir: str,
    tissue: str,
    pval_cutoff: float,
    log2fc_cutoff: float,
    fdr_method: str,
    *,
    feature_to_uniprot: dict[str, str] | None = None,
    sample_metadata: pd.DataFrame | None = None,
    design_covariates: Sequence[str] = (),
    min_observations_per_group: int = 2,
    flush: Callable[[str], None] = print,
) -> pd.DataFrame:
    """遗留 Mann-Whitney U 路径；不支持协变量，并显式标注不可推断行。"""

    if design_covariates:
        raise ValueError("Wilcoxon mode cannot adjust covariates; use dp_method='limma'")
    os.makedirs(output_dir, exist_ok=True)
    expression = _normalize_expression_axes(expression)
    fdr_method = normalize_fdr_method(fdr_method)
    metadata = _prepare_metadata(expression, group_map, sample_metadata)
    proteins = expression.index.astype(str)
    combined = pd.DataFrame({
        "Protein": proteins,
        "uniprot_id": [(feature_to_uniprot or {}).get(p, "") for p in proteins],
    })
    combined.index = proteins
    report_rows: list[dict[str, Any]] = []

    for numerator, denominator, name in comparisons:
        numerator_samples = metadata.index[metadata["condition"] == numerator].tolist()
        denominator_samples = metadata.index[metadata["condition"] == denominator].tolist()
        records: list[dict[str, Any]] = []
        for protein in proteins:
            x = expression.loc[protein, numerator_samples].dropna().to_numpy(dtype=float)
            y = expression.loc[protein, denominator_samples].dropna().to_numpy(dtype=float)
            fc = float(np.mean(x) - np.mean(y)) if x.size and y.size else np.nan
            if x.size >= min_observations_per_group and y.size >= min_observations_per_group:
                _, pvalue = stats.mannwhitneyu(x, y, alternative="two-sided")
                evidence = "inferential"
                status = "ok"
            else:
                pvalue = np.nan
                evidence = "insufficient"
                status = "insufficient_observations_per_group"
            records.append({
                "log2FoldChange": fc,
                "n_numerator": int(x.size),
                "n_denominator": int(y.size),
                "pvalue": pvalue,
                "evidence_level": evidence,
                "inference_status": status,
            })
        fitted = pd.DataFrame(records, index=proteins)
        fitted["padj"] = adjust_pvalues(fitted["pvalue"], fdr_method)
        combined[f"Log2FC_{name}"] = fitted["log2FoldChange"]
        combined[f"SE_{name}"] = np.nan
        combined[f"Tstat_{name}"] = np.nan
        combined[f"ResidualDF_{name}"] = np.nan
        combined[f"NObsNum_{name}"] = fitted["n_numerator"]
        combined[f"NObsDen_{name}"] = fitted["n_denominator"]
        combined[f"Pval_{name}"] = fitted["pvalue"]
        combined[f"Padj_{name}"] = fitted["padj"]
        combined[f"Evidence_{name}"] = fitted["evidence_level"]
        combined[f"InferenceStatus_{name}"] = fitted["inference_status"]
        tested = int(np.isfinite(fitted["pvalue"]).sum())
        route = "inferential" if tested else "insufficient_data"
        report_rows.append({
            "tissue": tissue,
            "comparison": name,
            "numerator": numerator,
            "denominator": denominator,
            "n_numerator": len(numerator_samples),
            "n_denominator": len(denominator_samples),
            "n_samples": len(numerator_samples) + len(denominator_samples),
            "design_formula": "Mann-Whitney U",
            "design_columns": "",
            "design_rank": np.nan,
            "n_design_columns": np.nan,
            "residual_df_complete_cases": np.nan,
            "route": route,
            "reason": "legacy non-parametric analysis",
            "dropped_constant_covariates": "",
            "n_features_total": len(expression),
            "n_features_tested": tested,
            "n_features_insufficient": len(expression) - tested,
            "prior_variance": np.nan,
            "prior_df": np.nan,
            "fdr_method": fdr_method,
            "differential_input": "log2_normalized_unimputed",
        })
        flush(f"  [{name}] Wilcoxon tested={tested}/{len(expression)}")

    combined.reset_index(drop=True, inplace=True)
    combined.to_csv(os.path.join(output_dir, COMBINED_CSV), index=False)
    pd.DataFrame(report_rows).to_csv(
        os.path.join(output_dir, DESIGN_REPORT_CSV), index=False,
    )
    for _, _, name in comparisons:
        per_comparison = _comparison_to_output(
            combined,
            name,
            pval_cutoff,
            log2fc_cutoff,
            "mann_whitney_u",
            fdr_method,
            "Mann-Whitney U",
        )
        per_comparison.to_csv(
            os.path.join(output_dir, f"dp_{tissue}_{name}.csv"), index=False,
        )
    return combined


def run_differential_analysis(
    expression: pd.DataFrame,
    group_map: dict[str, str],
    comparisons: list[tuple[str, str, str]],
    output_dir: str,
    tissue: str,
    pval_cutoff: float,
    log2fc_cutoff: float,
    fdr_method: str,
    dp_method: str = "limma",
    feature_to_uniprot: dict[str, str] | None = None,
    *,
    sample_metadata: pd.DataFrame | None = None,
    design_covariates: Sequence[str] = (),
    min_observations_per_group: int = 2,
    flush: Callable[[str], None] = print,
) -> pd.DataFrame:
    """按配置分派差异分析方法。"""

    common = dict(
        expression=expression,
        group_map=group_map,
        comparisons=comparisons,
        output_dir=output_dir,
        tissue=tissue,
        pval_cutoff=pval_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        fdr_method=fdr_method,
        feature_to_uniprot=feature_to_uniprot,
        sample_metadata=sample_metadata,
        design_covariates=design_covariates,
        min_observations_per_group=min_observations_per_group,
        flush=flush,
    )
    if dp_method == "limma":
        return run_limma_analysis(**common)
    if dp_method == "wilcoxon":
        return run_wilcoxon_analysis(**common)
    raise ValueError(f"Unsupported dp_method '{dp_method}'")


__all__ = [
    "COMBINED_CSV",
    "DESIGN_REPORT_CSV",
    "adjust_pvalues",
    "estimate_limma_priors",
    "run_differential_analysis",
    "run_limma_analysis",
    "run_wilcoxon_analysis",
]
