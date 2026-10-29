#!/usr/bin/env python3
"""跨实验多组学证据整合与项目级报告。

本模块只合并 estimand 完全可比、且提供 effect/SE 的实验级 EvidenceTable。
不同 estimand 会保留在同一项目报告中，但不会进入同一个 Meta 分析。
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml

from omicsweave.common.meta_analysis import random_effects_meta
from omicsweave.common.project import project_dir


ESTIMAND_FIELDS = (
    "species",
    "case",
    "control",
    "intervention_time",
    "tracer",
    "mfa_network",
    "primary_endpoint",
)

REQUIRED_EVIDENCE_COLUMNS = {
    "flux_id",
    "tissue",
    "comparison",
    "layer",
    "se",
}


def _bh_fdr(values: Iterable[float]) -> np.ndarray:
    pvalues = np.asarray(list(values), dtype=float)
    adjusted = np.full(pvalues.shape, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    if not valid.any():
        return adjusted
    finite = pvalues[valid]
    order = np.argsort(finite)
    ranked = finite[order]
    qvalues = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    qvalues = np.minimum.accumulate(qvalues[::-1])[::-1]
    restored = np.empty_like(qvalues)
    restored[order] = np.clip(qvalues, 0.0, 1.0)
    adjusted[valid] = restored
    return adjusted


def _as_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _resolve_path(value: str | Path, config_dir: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    repo_candidate = project_dir() / path
    if str(path).replace("\\", "/").startswith(("analysis_results/", "experiments/")):
        return repo_candidate.resolve()
    config_candidate = (config_dir / path).resolve()
    return config_candidate if config_candidate.exists() else repo_candidate.resolve()


@dataclass(frozen=True)
class ProjectExperiment:
    experiment_id: str
    results_dir: Path
    estimand: dict[str, Any]
    comparison_map: dict[str, str]
    enabled: bool = True

    @property
    def estimand_id(self) -> str:
        return _as_text(self.estimand.get("id"))

    def signature(self, layer: str) -> str:
        payload = {field: self.estimand.get(field) for field in ESTIMAND_FIELDS}
        platforms = self.estimand.get("platforms", {}) or {}
        payload["platform"] = platforms.get(layer, platforms.get("default", ""))
        return _canonical_json(payload)


@dataclass(frozen=True)
class ProjectIntegrationConfig:
    project_id: str
    config_path: Path
    output_dir: Path
    experiments: tuple[ProjectExperiment, ...]
    min_experiments: int
    fdr_threshold: float
    max_i_squared: float
    require_direction_consistency: bool
    missing_uncertainty_policy: str


@dataclass(frozen=True)
class ProjectIntegrationResult:
    config: ProjectIntegrationConfig
    manifest: pd.DataFrame
    compatibility_audit: pd.DataFrame
    evidence_long: pd.DataFrame
    meta_results: pd.DataFrame
    reaction_summary: pd.DataFrame
    report_path: Path


def _validate_estimand(experiment_id: str, estimand: dict[str, Any]) -> None:
    missing = [field for field in ("id", *ESTIMAND_FIELDS) if not _as_text(estimand.get(field))]
    if missing:
        raise ValueError(
            f"Experiment '{experiment_id}' estimand is missing required fields: "
            + ", ".join(missing)
        )
    platforms = estimand.get("platforms")
    if not isinstance(platforms, dict) or not platforms:
        raise ValueError(
            f"Experiment '{experiment_id}' estimand.platforms must be a non-empty mapping"
        )


def load_project_config(config_path: str | Path) -> ProjectIntegrationConfig:
    path = Path(config_path).resolve()
    if not path.exists():
        raise FileNotFoundError(f"Project integration config not found: {path}")
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Project config root must be a mapping")
    project_id = _as_text(raw.get("project_id"))
    if not project_id:
        raise ValueError("project_id is required")
    output_raw = raw.get("output_dir")
    if not output_raw:
        raise ValueError("output_dir is required")

    experiments: list[ProjectExperiment] = []
    seen: set[str] = set()
    for item in raw.get("experiments", []) or []:
        if not isinstance(item, dict):
            raise ValueError("Each experiments entry must be a mapping")
        experiment_id = _as_text(item.get("experiment_id"))
        if not experiment_id:
            raise ValueError("Each experiment requires experiment_id")
        if experiment_id in seen:
            raise ValueError(f"Duplicate experiment_id in project config: {experiment_id}")
        seen.add(experiment_id)
        estimand = dict(item.get("estimand", {}) or {})
        _validate_estimand(experiment_id, estimand)
        results_raw = item.get("results_dir")
        if not results_raw:
            raise ValueError(f"Experiment '{experiment_id}' requires results_dir")
        experiments.append(ProjectExperiment(
            experiment_id=experiment_id,
            results_dir=_resolve_path(results_raw, path.parent),
            estimand=estimand,
            comparison_map={
                str(key): str(value)
                for key, value in (item.get("comparison_map", {}) or {}).items()
                if value is not None
            },
            enabled=bool(item.get("enabled", True)),
        ))
    if not experiments:
        raise ValueError("At least one experiment is required")

    meta = raw.get("meta_analysis", {}) or {}
    missing_policy = _as_text(meta.get("missing_uncertainty_policy", "exclude")) or "exclude"
    if missing_policy not in {"exclude", "fail"}:
        raise ValueError("meta_analysis.missing_uncertainty_policy must be exclude or fail")
    return ProjectIntegrationConfig(
        project_id=project_id,
        config_path=path,
        output_dir=_resolve_path(output_raw, path.parent),
        experiments=tuple(experiments),
        min_experiments=max(2, int(meta.get("min_experiments", 2))),
        fdr_threshold=float(meta.get("fdr_threshold", 0.05)),
        max_i_squared=float(meta.get("max_i_squared", 75.0)),
        require_direction_consistency=bool(meta.get("require_direction_consistency", True)),
        missing_uncertainty_policy=missing_policy,
    )


def _load_experiment_evidence(experiment: ProjectExperiment) -> tuple[pd.DataFrame, dict[str, Any]]:
    evidence_path = experiment.results_dir / "common" / "evidence_table.csv"
    manifest = {
        "experiment_id": experiment.experiment_id,
        "estimand_id": experiment.estimand_id,
        "results_dir": str(experiment.results_dir),
        "evidence_table": str(evidence_path),
        "status": "disabled" if not experiment.enabled else "pending",
        "n_evidence_rows": 0,
        "n_meta_eligible_rows": 0,
        "message": "",
    }
    if not experiment.enabled:
        return pd.DataFrame(), manifest
    if not evidence_path.exists():
        manifest.update(status="not_assessable", message="evidence_table.csv not found")
        return pd.DataFrame(), manifest

    table = pd.read_csv(evidence_path)
    missing = REQUIRED_EVIDENCE_COLUMNS - set(table.columns)
    if missing:
        manifest.update(
            status="not_assessable",
            message="missing columns: " + ", ".join(sorted(missing)),
        )
        return pd.DataFrame(), manifest
    if "mechanism_aligned_effect" not in table.columns and "effect" not in table.columns:
        manifest.update(status="not_assessable", message="missing effect column")
        return pd.DataFrame(), manifest

    table = table.copy()
    table["experiment_id"] = experiment.experiment_id
    table["estimand_id"] = experiment.estimand_id
    table["canonical_comparison"] = table["comparison"].astype(str).map(
        lambda value: experiment.comparison_map.get(value, value)
    )
    table["effect_for_meta"] = pd.to_numeric(
        table.get("mechanism_aligned_effect", table.get("effect")), errors="coerce"
    )
    if "mechanism_aligned_effect" in table.columns and "effect" in table.columns:
        table["effect_for_meta"] = table["effect_for_meta"].fillna(
            pd.to_numeric(table["effect"], errors="coerce")
        )
    table["se"] = pd.to_numeric(table["se"], errors="coerce")
    table["estimand_signature"] = table["layer"].astype(str).map(experiment.signature)
    table["estimand_metadata"] = _canonical_json(experiment.estimand)
    table["meta_eligible"] = (
        table["effect_for_meta"].notna()
        & table["se"].notna()
        & np.isfinite(table["se"])
        & (table["se"] > 0)
    )
    manifest.update(
        status="loaded",
        n_evidence_rows=len(table),
        n_meta_eligible_rows=int(table["meta_eligible"].sum()),
    )
    return table, manifest


def _build_compatibility_audit(
    experiments: tuple[ProjectExperiment, ...],
    evidence: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    enabled = [experiment for experiment in experiments if experiment.enabled]
    observed_layers = sorted(evidence["layer"].dropna().astype(str).unique()) if not evidence.empty else []
    layers = observed_layers or ["mfa", "transcriptomics", "proteomics", "metabolomics"]
    for index, left in enumerate(enabled):
        for right in enabled[index + 1:]:
            for layer in layers:
                differences: list[str] = []
                for field in ESTIMAND_FIELDS:
                    if _as_text(left.estimand.get(field)) != _as_text(right.estimand.get(field)):
                        differences.append(field)
                left_platform = (left.estimand.get("platforms", {}) or {}).get(layer, "")
                right_platform = (right.estimand.get("platforms", {}) or {}).get(layer, "")
                if _as_text(left_platform) != _as_text(right_platform):
                    differences.append(f"platforms.{layer}")
                same_id = left.estimand_id == right.estimand_id
                rows.append({
                    "experiment_a": left.experiment_id,
                    "experiment_b": right.experiment_id,
                    "layer": layer,
                    "estimand_id_a": left.estimand_id,
                    "estimand_id_b": right.estimand_id,
                    "comparable": same_id and not differences,
                    "differences": ";".join(differences),
                    "note": (
                        "same estimand"
                        if same_id and not differences
                        else "different estimand_id"
                        if not same_id
                        else "estimand metadata mismatch"
                    ),
                })
    return pd.DataFrame(rows, columns=[
        "experiment_a", "experiment_b", "layer", "estimand_id_a", "estimand_id_b",
        "comparable", "differences", "note",
    ])


def _validate_estimand_ids(experiments: tuple[ProjectExperiment, ...]) -> None:
    by_id: dict[str, list[ProjectExperiment]] = {}
    for experiment in experiments:
        if experiment.enabled:
            by_id.setdefault(experiment.estimand_id, []).append(experiment)
    for estimand_id, group in by_id.items():
        for field in ESTIMAND_FIELDS:
            values = {_canonical_json(experiment.estimand.get(field)) for experiment in group}
            if len(values) > 1:
                raise ValueError(
                    f"Estimand '{estimand_id}' has inconsistent field '{field}' across experiments"
                )
        layers = {
            key
            for experiment in group
            for key in (experiment.estimand.get("platforms", {}) or {}).keys()
        }
        for layer in layers:
            values = {
                _as_text((experiment.estimand.get("platforms", {}) or {}).get(layer))
                for experiment in group
            }
            if len(values) > 1:
                raise ValueError(
                    f"Estimand '{estimand_id}' has inconsistent platform for layer '{layer}'"
                )


def _run_meta_analysis(
    evidence: pd.DataFrame,
    config: ProjectIntegrationConfig,
) -> pd.DataFrame:
    columns = [
        "estimand_id", "estimand_signature", "canonical_comparison", "tissue",
        "flux_id", "layer", "pooled_effect", "pooled_se", "ci_lower", "ci_upper",
        "pvalue", "fdr", "q_statistic", "i_squared", "tau_squared",
        "n_experiments", "experiment_ids", "direction_consistent",
        "replication_status", "meta_note",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=columns)

    eligible = evidence[evidence["meta_eligible"]].copy()
    if config.missing_uncertainty_policy == "fail" and (~evidence["meta_eligible"]).any():
        bad = evidence.loc[~evidence["meta_eligible"], [
            "experiment_id", "flux_id", "tissue", "comparison", "layer",
        ]]
        raise ValueError(
            "Evidence rows without finite effect/SE are not allowed:\n"
            + bad.head(20).to_string(index=False)
        )

    keys = [
        "estimand_id", "estimand_signature", "canonical_comparison",
        "tissue", "flux_id", "layer",
    ]
    rows: list[dict[str, Any]] = []
    for key, group in eligible.groupby(keys, dropna=False, sort=True):
        per_experiment = (
            group.sort_values("se")
            .drop_duplicates(subset=["experiment_id"], keep="first")
        )
        if per_experiment["experiment_id"].nunique() < config.min_experiments:
            continue
        meta = random_effects_meta(
            per_experiment["effect_for_meta"].to_numpy(dtype=float),
            per_experiment["se"].to_numpy(dtype=float),
            per_experiment["experiment_id"].astype(str).tolist(),
        )
        rows.append({
            **dict(zip(keys, key, strict=True)),
            "pooled_effect": meta.pooled_effect,
            "pooled_se": meta.pooled_se,
            "ci_lower": meta.ci_lower,
            "ci_upper": meta.ci_upper,
            "pvalue": meta.pvalue,
            "fdr": np.nan,
            "q_statistic": meta.q_statistic,
            "i_squared": meta.i_squared,
            "tau_squared": meta.tau_squared,
            "n_experiments": meta.k,
            "experiment_ids": ";".join(meta.source_ids),
            "direction_consistent": meta.direction_consistent,
            "replication_status": "pending",
            "meta_note": "random_effects_meta",
        })
    result = pd.DataFrame(rows, columns=columns)
    if result.empty:
        return result

    fdr_groups = ["estimand_id", "canonical_comparison", "tissue", "layer"]
    for _, indices in result.groupby(fdr_groups, dropna=False).groups.items():
        result.loc[indices, "fdr"] = _bh_fdr(result.loc[indices, "pvalue"])
    replicated = (
        (result["n_experiments"] >= config.min_experiments)
        & (result["fdr"] <= config.fdr_threshold)
        & (result["i_squared"].fillna(0.0) <= config.max_i_squared)
    )
    if config.require_direction_consistency:
        replicated &= result["direction_consistent"].astype(bool)
    result.loc[replicated, "replication_status"] = "replicated"
    result.loc[~replicated, "replication_status"] = "not_replicated"
    result.loc[result["i_squared"] > config.max_i_squared, "meta_note"] = "high_heterogeneity"
    result.loc[~result["direction_consistent"], "meta_note"] = "direction_conflict"
    return result.sort_values(
        ["replication_status", "fdr", "i_squared"],
        ascending=[True, True, True],
    ).reset_index(drop=True)


def _reaction_summary(meta: pd.DataFrame, fdr_threshold: float) -> pd.DataFrame:
    columns = [
        "estimand_id", "canonical_comparison", "tissue", "flux_id",
        "n_meta_layers", "n_replicated_layers", "support_score", "conflict_score",
        "coverage_score", "evidence_grade", "replicated_layers", "experiment_ids",
    ]
    if meta.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    keys = ["estimand_id", "canonical_comparison", "tissue", "flux_id"]
    for key, group in meta.groupby(keys, dropna=False, sort=True):
        effects = pd.to_numeric(group["pooled_effect"], errors="coerce")
        significant = pd.to_numeric(group["fdr"], errors="coerce") <= fdr_threshold
        significant_effects = effects[significant]
        positive = effects.clip(lower=0.0)
        negative = (-effects.clip(upper=0.0))
        replicated = group["replication_status"] == "replicated"
        n_layers = int(group["layer"].nunique())
        rows.append({
            **dict(zip(keys, key, strict=True)),
            "n_meta_layers": n_layers,
            "n_replicated_layers": int(group.loc[replicated, "layer"].nunique()),
            "support_score": float(np.tanh(positive[significant].mean())) if (positive[significant] > 0).any() else 0.0,
            "conflict_score": float(np.tanh(negative[significant].mean())) if (negative[significant] > 0).any() else 0.0,
            "coverage_score": min(1.0, n_layers / 4.0),
            "evidence_grade": (
                "Replicated"
                if replicated.any()
                else "Conflicting"
                if (significant_effects > 0).any() and (significant_effects < 0).any()
                else "Convergent"
                if (
                    int((significant_effects > 0).sum()) >= 2
                    or int((significant_effects < 0).sum()) >= 2
                )
                else "Single-layer"
                if significant.any()
                else "No-signal"
            ),
            "replicated_layers": ";".join(sorted(group.loc[replicated, "layer"].astype(str).unique())),
            "experiment_ids": ";".join(sorted({
                item
                for value in group["experiment_ids"].astype(str)
                for item in value.split(";")
                if item
            })),
        })
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["n_replicated_layers", "support_score", "coverage_score"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def _generate_project_report(
    config: ProjectIntegrationConfig,
    manifest: pd.DataFrame,
    audit: pd.DataFrame,
    meta: pd.DataFrame,
    summary: pd.DataFrame,
) -> Path:
    lines = [
        f"# 多组学项目整合报告：{config.project_id}\n\n",
        f"**生成时间：** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}  \n",
        f"**登记实验：** {len(manifest)}  \n",
        f"**可用于随机效应 Meta 的实验：** {(manifest['n_meta_eligible_rows'] > 0).sum()}  \n\n",
        "## 1. 实验与输入状态\n\n",
        "| 实验 | Estimand | 状态 | Evidence 行数 | Meta 可用行数 | 说明 |\n",
        "|---|---|---|---:|---:|---|\n",
    ]
    for _, row in manifest.iterrows():
        lines.append(
            f"| {row['experiment_id']} | {row['estimand_id']} | {row['status']} "
            f"| {row['n_evidence_rows']} | {row['n_meta_eligible_rows']} | {row['message']} |\n"
        )

    lines.extend(["\n## 2. 可比性审计\n\n"])
    if audit.empty:
        lines.append("只有一个启用实验，无法进行成对可比性审计。\n\n")
    else:
        n_comparable = int(audit["comparable"].sum())
        lines.append(
            f"共检查 **{len(audit)}** 个实验对×组学层，其中 **{n_comparable}** 个可比。"
            "不可比组合保留在 `compatibility_audit.csv`，未进入同一 Meta 分析。\n\n"
        )

    lines.extend(["## 3. 跨实验随机效应 Meta\n\n"])
    if meta.empty:
        lines.append("没有满足最少实验数且同时具有 effect/SE 的可比证据。\n\n")
    else:
        n_replicated = int((meta["replication_status"] == "replicated").sum())
        lines.append(
            f"完成 **{len(meta)}** 个反应层级 Meta，其中 **{n_replicated}** 个达到 Replicated 标准。\n\n"
        )
        lines.append(
            "Replicated 要求：至少两个可比实验、随机效应 Meta FDR 达标、方向一致，"
            f"且 I² ≤ {config.max_i_squared:g}%。\n\n"
        )
        lines.append("| 反应 | 组织 | 比较 | 层 | 效应 | SE | FDR | I² | 实验数 | 结论 |\n")
        lines.append("|---|---|---|---|---:|---:|---:|---:|---:|---|\n")
        for _, row in meta.head(30).iterrows():
            lines.append(
                f"| {row['flux_id']} | {row['tissue']} | {row['canonical_comparison']} "
                f"| {row['layer']} | {row['pooled_effect']:.3f} | {row['pooled_se']:.3f} "
                f"| {row['fdr']:.3g} | {row['i_squared']:.1f} | {row['n_experiments']} "
                f"| {row['replication_status']} |\n"
            )

    lines.extend(["\n## 4. 项目级反应证据画像\n\n"])
    if summary.empty:
        lines.append("无可汇总的项目级反应结果。\n\n")
    else:
        grade_counts = summary["evidence_grade"].value_counts()
        for grade in ["Replicated", "Convergent", "Single-layer", "Conflicting"]:
            lines.append(f"- {grade}: {int(grade_counts.get(grade, 0))}\n")
        lines.append("\n完整结果见 `reaction_meta_summary.csv`。\n\n")

    lines.extend([
        "## 5. 统计边界\n\n",
        "- 不同物种、case/control 定义、时间、示踪剂、MFA 网络、平台或主要终点不会被数学归一化为可比。\n",
        "- 缺失 effect/SE 的证据标记为不可用于 Meta，而不是填零。\n",
        "- 本报告整合的是独立实验的标准化效应证据；样本级 batch correction 必须在各组学差异分析前完成。\n",
        "- `Replicated` 是跨实验结论，单实验的 Convergent/Exploratory 不会被自动升级。\n\n",
        "## 6. 输出文件\n\n",
        "- `experiment_manifest.csv`\n",
        "- `compatibility_audit.csv`\n",
        "- `project_evidence_long.csv`\n",
        "- `cross_experiment_meta.csv`\n",
        "- `reaction_meta_summary.csv`\n",
    ])
    path = config.output_dir / "project_multiomics_report.md"
    path.write_text("".join(lines), encoding="utf-8")
    return path


def run_project_integration(config_path: str | Path) -> ProjectIntegrationResult:
    """运行项目级整合并输出一份集中报告。"""
    config = load_project_config(config_path)
    _validate_estimand_ids(config.experiments)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    evidence_tables: list[pd.DataFrame] = []
    manifest_rows: list[dict[str, Any]] = []
    for experiment in config.experiments:
        table, manifest = _load_experiment_evidence(experiment)
        manifest_rows.append(manifest)
        if not table.empty:
            evidence_tables.append(table)
    manifest = pd.DataFrame(manifest_rows)
    evidence = pd.concat(evidence_tables, ignore_index=True) if evidence_tables else pd.DataFrame()
    audit = _build_compatibility_audit(config.experiments, evidence)
    meta = _run_meta_analysis(evidence, config)
    summary = _reaction_summary(meta, config.fdr_threshold)

    manifest.to_csv(config.output_dir / "experiment_manifest.csv", index=False, encoding="utf-8-sig")
    audit.to_csv(config.output_dir / "compatibility_audit.csv", index=False, encoding="utf-8-sig")
    evidence.to_csv(config.output_dir / "project_evidence_long.csv", index=False, encoding="utf-8-sig")
    meta.to_csv(config.output_dir / "cross_experiment_meta.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(config.output_dir / "reaction_meta_summary.csv", index=False, encoding="utf-8-sig")
    report_path = _generate_project_report(config, manifest, audit, meta, summary)

    return ProjectIntegrationResult(
        config=config,
        manifest=manifest,
        compatibility_audit=audit,
        evidence_long=evidence,
        meta_results=meta,
        reaction_summary=summary,
        report_path=report_path,
    )
