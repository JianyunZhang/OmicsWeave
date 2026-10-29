"""
Transcriptomics multi-batch / multi-protocol integration.

Outputs are written under:
  <output_root>/integration/
    within_protocol/<protocol>/<tissue>/de_<comparison>_integrated.csv
    cross_protocol/<tissue>/de_<comparison>_consensus.csv
    integration_report.md
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

from omicsweave.core.statistics.meta_analysis import random_effects_meta

from .experiment_config import (  # type: ignore
    discover_experiment_protocols,
    load_integration_config,
    resolve_output_root,
)
from . import provenance as provenance_mod  # type: ignore
from . import joint_model as joint_model_mod  # type: ignore


def _flush_print(*args, **kwargs) -> None:
    print(*args, **kwargs, flush=True)


def _discover_batch_de_results(
    output_root: Path,
    protocol_id: str,
    batch_id: str,
) -> dict[tuple[str, str], Path]:
    """
    Discover DE result CSV files for one batch.

    Returns
    -------
    dict[(tissue, comparison)] = path_to_csv
    """
    de_root = output_root / "runs" / protocol_id / batch_id / "04_de_results"
    if not de_root.exists():
        return {}

    found: dict[tuple[str, str], Path] = {}

    # Preferred structure: 04_de_results/<tissue>/de_<comparison>[ _wilcoxon].csv
    for csv_file in sorted(de_root.glob("*/de_*.csv")):
        tissue = csv_file.parent.name
        stem = csv_file.stem
        comp = stem.replace("de_", "")
        if comp.endswith("_wilcoxon"):
            comp = comp[: -len("_wilcoxon")]
        found[(tissue, comp)] = csv_file

    return found


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    out = np.full(n, np.nan, dtype=float)
    valid = np.isfinite(pvals)
    if not valid.any():
        return out

    pv = pvals[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    q = np.empty_like(ranked)
    prev = 1.0
    for i in range(len(ranked) - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * len(ranked) / rank
        prev = min(prev, val)
        q[i] = prev
    q = np.clip(q, 0.0, 1.0)

    restored = np.empty_like(q)
    restored[order] = q
    out[valid] = restored
    return out


def _read_de_table(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if df.empty:
        return pd.DataFrame(columns=["gene", "log2FoldChange", "pvalue", "padj", "lfcSE"])

    if "gene" not in df.columns:
        first_col = str(df.columns[0])
        df = df.rename(columns={first_col: "gene"})

    keep_cols = ["gene", "log2FoldChange", "pvalue", "padj", "lfcSE"]
    for col in keep_cols:
        if col not in df.columns:
            df[col] = np.nan

    out = df[keep_cols].copy()
    out["gene"] = out["gene"].astype(str)
    for col in ["log2FoldChange", "pvalue", "padj", "lfcSE"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=["gene", "log2FoldChange"])
    return out


def _estimate_se_with_source(
    log2fc: float,
    pvalue: float,
    lfc_se: float,
) -> tuple[float, str]:
    if np.isfinite(lfc_se) and lfc_se > 0:
        return float(lfc_se), "reported"
    if not np.isfinite(log2fc) or not np.isfinite(pvalue) or pvalue <= 0 or pvalue >= 1:
        return np.nan, "unavailable"
    z = norm.isf(pvalue / 2.0)
    if not np.isfinite(z) or z <= 0:
        return np.nan, "unavailable"
    return abs(log2fc) / z, "pvalue_normal_approximation"


def _estimate_se(log2fc: float, pvalue: float, lfc_se: float) -> float:
    """兼容旧调用；新集成流程同时记录 SE 的来源。"""
    return _estimate_se_with_source(log2fc, pvalue, lfc_se)[0]


def _heterogeneity_threshold_percent(value: Any) -> float:
    """接受旧版 0--1 比例或新版 0--100 百分比阈值。"""
    threshold = float(value)
    return threshold * 100.0 if 0.0 <= threshold <= 1.0 else threshold


def _combine_effects(
    effects: np.ndarray,
    ses: np.ndarray,
    method: str = "inverse_variance",
) -> dict[str, float]:
    k = len(effects)
    if k == 0:
        return {}

    effects = effects.astype(float)
    ses = ses.astype(float)
    valid = np.isfinite(effects)
    effects = effects[valid]
    ses = ses[valid]
    k = len(effects)
    if k == 0:
        return {}

    if k == 1:
        se = ses[0] if np.isfinite(ses[0]) and ses[0] > 0 else np.nan
        pvalue = 2 * norm.sf(abs(effects[0] / se)) if np.isfinite(se) and se > 0 else np.nan
        return {
            "pooled_log2FoldChange": float(effects[0]),
            "pooled_se": float(se) if np.isfinite(se) else np.nan,
            "pooled_pvalue": float(pvalue) if np.isfinite(pvalue) else np.nan,
            "i2": np.nan,
            "k": 1.0,
        }

    if method == "random_effects_meta":
        se_valid = np.isfinite(ses) & (ses > 0)
        if not se_valid.any():
            return {}
        meta = random_effects_meta(effects[se_valid], ses[se_valid])
        return {
            "pooled_log2FoldChange": meta.pooled_effect,
            "pooled_se": meta.pooled_se,
            "pooled_pvalue": meta.pvalue,
            # 公共契约统一使用 0--100 百分比，避免与 MFA/项目级 Meta 混用单位。
            "i2": meta.i_squared,
            "tau_squared": meta.tau_squared,
            "q_statistic": meta.q_statistic,
            "k": float(meta.k),
        }

    use_inverse = method == "inverse_variance"
    se_valid = np.isfinite(ses) & (ses > 0)
    if use_inverse and se_valid.any():
        w = 1.0 / (ses[se_valid] ** 2)
        e = effects[se_valid]
        pooled = float(np.sum(w * e) / np.sum(w))
        pooled_se = float(np.sqrt(1.0 / np.sum(w)))
        z = pooled / pooled_se if pooled_se > 0 else np.nan
        pooled_p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan

        theta_fe = pooled
        q = float(np.sum(w * (e - theta_fe) ** 2))
        df = max(1, len(e) - 1)
        i2 = 100.0 * max(0.0, (q - df) / q) if q > 0 else 0.0
        return {
            "pooled_log2FoldChange": pooled,
            "pooled_se": pooled_se,
            "pooled_pvalue": pooled_p,
            "i2": float(i2),
            "tau_squared": 0.0,
            "q_statistic": q,
            "k": float(len(e)),
        }

    pooled = float(np.mean(effects))
    pooled_se = float(np.std(effects, ddof=1) / np.sqrt(k)) if k > 1 else np.nan
    z = pooled / pooled_se if np.isfinite(pooled_se) and pooled_se > 0 else np.nan
    pooled_p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return {
        "pooled_log2FoldChange": pooled,
        "pooled_se": pooled_se if np.isfinite(pooled_se) else np.nan,
        "pooled_pvalue": pooled_p if np.isfinite(pooled_p) else np.nan,
        "i2": np.nan,
        "tau_squared": np.nan,
        "q_statistic": np.nan,
        "k": float(k),
    }


def check_batch_consistency(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
    policy: str = "fail",
) -> tuple[bool, dict[str, Any], str]:
    """
    整合前校验各批次的参考与技术参数是否一致。

    跨批次合并的前提是技术链路一致：不同的转录组参考版本、kallisto 版本或
    链特异性设置会让基因层计数在定义上不可比，此时 Meta 分析只会把技术
    差异当成生物学差异。

    Parameters
    ----------
    policy : {"fail", "warn", "skip"}
        ``fail`` 不一致时拒绝整合；``warn`` 仅告警继续；``skip`` 不做校验。

    Returns
    -------
    (should_proceed, comparison_result, report_text)
    """
    resolved_policy = str(policy or "fail").strip().lower()
    if resolved_policy == "skip":
        return True, {"consistent": True, "skipped": True}, "一致性校验已按配置跳过。"

    provenances = {
        batch_id: provenance_mod.load_run_provenance(
            output_root / "runs" / protocol_id / batch_id
        )
        for batch_id in batch_ids
    }
    result = provenance_mod.compare_provenance(provenances)
    report = provenance_mod.format_provenance_report(result)

    if result["consistent"]:
        return True, result, report

    if resolved_policy == "warn":
        return True, result, report
    return False, result, report


def within_protocol_meta_analysis(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
    integration_config: dict[str, Any],
) -> dict[tuple[str, str], pd.DataFrame]:
    wp_cfg = integration_config.get("within_protocol", {}) or {}
    min_batches = int(wp_cfg.get("min_batches", 2))
    configured_method = str(wp_cfg.get("method", "random_effects_meta")).lower()
    if configured_method == "random_effects_meta":
        weight_method = "random_effects_meta"
    else:
        weight_method = str(
            (integration_config.get("weighting", {}) or {}).get("batch", "inverse_variance")
        ).lower()
        if weight_method not in {"equal", "inverse_variance"}:
            weight_method = "inverse_variance"

    log2fc_cutoff = float((wp_cfg.get("log2fc_threshold", 1.0)))
    fdr_cutoff = float((wp_cfg.get("fdr_threshold", 0.05)))
    heterogeneity_threshold = _heterogeneity_threshold_percent(
        wp_cfg.get("heterogeneity_threshold", 50.0)
    )

    batch_tables: dict[str, dict[tuple[str, str], pd.DataFrame]] = {}
    all_keys: set[tuple[str, str]] = set()
    for batch_id in batch_ids:
        paths = _discover_batch_de_results(output_root, protocol_id, batch_id)
        one_batch: dict[tuple[str, str], pd.DataFrame] = {}
        for key, path in paths.items():
            table = _read_de_table(path)
            if table.empty:
                continue
            table["_batch_id"] = batch_id
            estimated = table.apply(
                lambda r: _estimate_se_with_source(
                    r["log2FoldChange"], r["pvalue"], r["lfcSE"]
                ),
                axis=1,
            )
            table["_se"] = estimated.map(lambda value: value[0])
            table["_se_source"] = estimated.map(lambda value: value[1])
            one_batch[key] = table
            all_keys.add(key)
        if one_batch:
            batch_tables[batch_id] = one_batch

    out_root = output_root / "integration" / "within_protocol" / protocol_id
    out_root.mkdir(parents=True, exist_ok=True)
    integrated: dict[tuple[str, str], pd.DataFrame] = {}

    for key in sorted(all_keys):
        tissue, comp = key
        rows = []
        for batch_id, tables in batch_tables.items():
            if key not in tables:
                continue
            t = tables[key][["gene", "log2FoldChange", "_se", "_se_source"]].copy()
            t = t.rename(columns={"log2FoldChange": "effect", "_se": "se"})
            t["batch_id"] = batch_id
            rows.append(t)
        if not rows:
            continue

        joined = pd.concat(rows, ignore_index=True)
        out_rows = []
        for gene, grp in joined.groupby("gene"):
            observed = (
                grp.sort_values("se", na_position="last")
                .drop_duplicates("batch_id")
                .copy()
            )
            observed_ids = sorted(observed["batch_id"].astype(str).tolist())
            eligible = observed[
                np.isfinite(observed["se"]) & (observed["se"] > 0)
            ].copy()
            eligible_ids = sorted(eligible["batch_id"].astype(str).tolist())
            effects = observed["effect"].to_numpy(dtype=float)
            n_observed = len(observed_ids)
            n_eligible = len(eligible_ids)
            meta_eligible = n_eligible >= min_batches

            if n_eligible == 0:
                pooled_effect = float(np.nanmean(effects))
                meta = {
                    "pooled_log2FoldChange": pooled_effect,
                    "pooled_se": np.nan,
                    "pooled_pvalue": np.nan,
                    "i2": np.nan,
                    "tau_squared": np.nan,
                    "q_statistic": np.nan,
                    "k": 0.0,
                }
                meta_note = "not_meta_eligible_missing_uncertainty"
                se_source = "unavailable"
            elif n_eligible < min_batches:
                if n_eligible == 1:
                    meta = _combine_effects(
                        eligible["effect"].to_numpy(dtype=float),
                        eligible["se"].to_numpy(dtype=float),
                        method=weight_method,
                    )
                    meta_note = "single_batch_passthrough"
                else:
                    meta = {
                        "pooled_log2FoldChange": float(eligible["effect"].mean()),
                        "pooled_se": np.nan,
                        "pooled_pvalue": np.nan,
                        "i2": np.nan,
                        "tau_squared": np.nan,
                        "q_statistic": np.nan,
                        "k": float(n_eligible),
                    }
                    meta_note = "insufficient_batches_for_meta"
                se_source = ";".join(sorted(eligible["_se_source"].astype(str).unique()))
            else:
                meta = _combine_effects(
                    eligible["effect"].to_numpy(dtype=float),
                    eligible["se"].to_numpy(dtype=float),
                    method=weight_method,
                )
                if not meta:
                    continue
                meta_note = (
                    "high_heterogeneity"
                    if np.isfinite(meta["i2"])
                    and meta["i2"] > heterogeneity_threshold
                    else weight_method
                )
                se_source = ";".join(sorted(eligible["_se_source"].astype(str).unique()))

            nonzero_effects = effects[np.abs(effects) > 0]
            direction_consistent = (
                bool(np.all(nonzero_effects > 0) or np.all(nonzero_effects < 0))
                if len(nonzero_effects) > 0
                else True
            )
            out_rows.append({
                "gene": gene,
                "pooled_log2FoldChange": meta["pooled_log2FoldChange"],
                "pooled_se": meta["pooled_se"],
                "pooled_pvalue": meta["pooled_pvalue"],
                "pooled_fdr": np.nan,
                "i2": meta["i2"],
                "tau_squared": meta.get("tau_squared", 0.0),
                "q_statistic": meta.get("q_statistic", np.nan),
                "n_batches": int(meta.get("k", n_eligible)),
                "n_batches_observed": n_observed,
                "n_batches_eligible": n_eligible,
                "direction_consistent": direction_consistent,
                "batch_ids": ";".join(observed_ids),
                "eligible_batch_ids": ";".join(eligible_ids),
                "meta_eligible": meta_eligible,
                "meta_note": meta_note,
                "se_source": se_source,
                "effect_scale": "log2FoldChange",
            })

        if not out_rows:
            continue

        out_df = pd.DataFrame(out_rows)
        out_df["pooled_fdr"] = _bh_fdr(out_df["pooled_pvalue"].to_numpy(dtype=float))
        out_df["significant"] = (
            (out_df["pooled_fdr"] <= fdr_cutoff)
            & (out_df["pooled_log2FoldChange"].abs() >= log2fc_cutoff)
        )
        out_df = out_df.sort_values(
            by=["pooled_fdr", "pooled_pvalue", "pooled_log2FoldChange"],
            ascending=[True, True, False],
        )

        tissue_dir = out_root / tissue
        tissue_dir.mkdir(parents=True, exist_ok=True)
        out_file = tissue_dir / f"de_{comp}_integrated.csv"
        out_df.to_csv(out_file, index=False)
        integrated[key] = out_df
        _flush_print(f"  [within] {protocol_id}/{tissue}/{comp}: {len(out_df)} genes -> {out_file}")

    return integrated


def cross_protocol_consensus(
    output_root: Path,
    protocol_ids: list[str],
    integration_config: dict[str, Any],
) -> dict[tuple[str, str], pd.DataFrame]:
    cp_cfg = integration_config.get("cross_protocol", {}) or {}
    min_protocols = int(cp_cfg.get("min_protocols", 2))
    require_same_direction = bool(cp_cfg.get("require_same_direction", True))
    fdr_cutoff = float(cp_cfg.get("fdr_threshold", 0.05))

    per_protocol: dict[str, dict[tuple[str, str], pd.DataFrame]] = {}
    all_keys: set[tuple[str, str]] = set()

    for protocol_id in protocol_ids:
        root = output_root / "integration" / "within_protocol" / protocol_id
        protocol_tables: dict[tuple[str, str], pd.DataFrame] = {}
        if root.exists():
            for csv_file in sorted(root.glob("*/de_*_integrated.csv")):
                tissue = csv_file.parent.name
                comp = csv_file.stem.replace("de_", "").replace("_integrated", "")
                df = pd.read_csv(csv_file)
                if df.empty or "gene" not in df.columns:
                    continue
                protocol_tables[(tissue, comp)] = df
                all_keys.add((tissue, comp))
        if protocol_tables:
            per_protocol[protocol_id] = protocol_tables

    out_root = output_root / "integration" / "cross_protocol"
    out_root.mkdir(parents=True, exist_ok=True)
    consensus: dict[tuple[str, str], pd.DataFrame] = {}

    for key in sorted(all_keys):
        tissue, comp = key
        rows = []
        for protocol_id, tables in per_protocol.items():
            if key not in tables:
                continue
            t = tables[key][["gene", "pooled_log2FoldChange", "pooled_fdr"]].copy()
            t = t.rename(columns={"pooled_log2FoldChange": "effect", "pooled_fdr": "fdr"})
            t["protocol_id"] = protocol_id
            rows.append(t)
        if not rows:
            continue

        joined = pd.concat(rows, ignore_index=True)
        out_rows = []
        for gene, grp in joined.groupby("gene"):
            effects = grp["effect"].to_numpy(dtype=float)
            protocol_count = len(grp["protocol_id"].unique())
            if protocol_count < min_protocols:
                continue

            nonzero = effects[np.abs(effects) > 0]
            direction_consistent = (
                bool(np.all(nonzero > 0) or np.all(nonzero < 0))
                if len(nonzero) > 0
                else True
            )
            if require_same_direction and not direction_consistent:
                continue

            out_rows.append({
                "gene": gene,
                "consensus_log2FoldChange": float(np.mean(effects)),
                "n_protocols": protocol_count,
                "direction_consistent": direction_consistent,
                "significant_protocols": int((grp["fdr"] <= fdr_cutoff).sum()),
                "protocol_ids": ",".join(sorted(grp["protocol_id"].unique().tolist())),
                "protocol_effects": ";".join(
                    f"{pid}:{eff:+.4f}"
                    for pid, eff in zip(grp["protocol_id"].tolist(), grp["effect"].tolist())
                ),
            })

        if not out_rows:
            continue

        out_df = pd.DataFrame(out_rows)
        out_df = out_df.sort_values(
            by=["significant_protocols", "n_protocols", "consensus_log2FoldChange"],
            ascending=[False, False, False],
        )
        tissue_dir = out_root / tissue
        tissue_dir.mkdir(parents=True, exist_ok=True)
        out_file = tissue_dir / f"de_{comp}_consensus.csv"
        out_df.to_csv(out_file, index=False)
        consensus[key] = out_df
        _flush_print(f"  [cross] {tissue}/{comp}: {len(out_df)} genes -> {out_file}")

    return consensus


def generate_integration_report(
    output_root: Path,
    within_result: dict[str, dict[tuple[str, str], pd.DataFrame]],
    cross_result: dict[tuple[str, str], pd.DataFrame],
    experiment_id: str,
    consistency_reports: dict[str, str] | None = None,
    joint_summary: dict[str, Any] | None = None,
) -> Path:
    out_file = output_root / "integration" / "integration_report.md"
    out_file.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    lines.append(f"# Transcriptomics Integration Report: {experiment_id}")
    lines.append("")

    if joint_summary and joint_summary.get("protocols"):
        lines.append("## Joint cross-batch model (primary route)")
        lines.append("")
        lines.append(
            "批次内分组平衡时，联合 `~batch + condition` 的统计效力高于事后 Meta 分析："
            "批次作为固定效应被吸收，离散度在全部样本上共享估计。"
            "设计不可识别（批次与分组混杂）的组合会被预检拒绝并退回 Meta 路径。"
        )
        lines.append("")
        for protocol_id, summary in sorted(joint_summary["protocols"].items()):
            status = summary.get("status", "unknown")
            lines.append(f"- **{protocol_id}**: status=`{status}`")
            if summary.get("merge_info"):
                mi = summary["merge_info"]
                lines.append(
                    f"  - merged {mi['n_samples']} samples × {mi['n_genes']} genes "
                    f"from {mi['n_batches']} batches ({', '.join(mi['batch_ids'])})"
                )
                if mi.get("gene_overlap_warning"):
                    lines.append(f"  - ⚠ {mi['gene_overlap_warning']}")
            if "n_joint_eligible" in summary:
                lines.append(
                    f"  - design pre-check: {summary['n_joint_eligible']} joint-eligible, "
                    f"{summary.get('n_meta_only', 0)} meta-only, "
                    f"{summary.get('n_insufficient', 0)} insufficient"
                )
            if summary.get("design_report_path"):
                lines.append(f"  - design report: `{summary['design_report_path']}`")
            if summary.get("error"):
                lines.append(f"  - error: {summary['error']}")
            for key, n_genes in sorted((summary.get("results") or {}).items()):
                lines.append(f"  - {key}: {n_genes} genes")
        lines.append("")

    if consistency_reports:
        lines.append("## Batch consistency check")
        lines.append("")
        lines.append("跨批次合并要求参考版本与技术参数一致，否则技术差异会被当成生物学差异。")
        lines.append("")
        for protocol_id, report in sorted(consistency_reports.items()):
            lines.append(f"- **{protocol_id}**: {report}")
        lines.append("")

    lines.append("## Within-protocol integration")
    lines.append("")
    if not within_result:
        lines.append("- No within-protocol integrated outputs were generated.")
    else:
        for protocol_id, tables in sorted(within_result.items()):
            lines.append(f"- **{protocol_id}**: {len(tables)} tissue/comparison tables")
            for (tissue, comp), df in sorted(tables.items()):
                n_sig = int(df["significant"].sum()) if "significant" in df.columns else 0
                lines.append(f"  - {tissue}/{comp}: {len(df)} genes, {n_sig} significant")

    lines.append("")
    lines.append("## Cross-protocol consensus")
    lines.append("")
    if not cross_result:
        lines.append("- No cross-protocol consensus outputs were generated.")
    else:
        for (tissue, comp), df in sorted(cross_result.items()):
            lines.append(f"- {tissue}/{comp}: {len(df)} consensus genes")

    out_file.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return out_file


def run_integration(experiment_id: str) -> None:
    integration_cfg = load_integration_config(experiment_id)
    if not integration_cfg.get("enabled", False):
        _flush_print("[INFO] integration.enabled=false, skip transcriptomics integration.")
        return

    output_root = resolve_output_root(experiment_id)
    protocols = discover_experiment_protocols(experiment_id)

    protocol_batch_map: dict[str, list[str]] = {}
    for proto in protocols:
        if not proto.get("enabled", False):
            continue
        protocol_id = proto.get("protocol_id")
        if not protocol_id:
            continue
        batch_ids = [
            b["batch_id"]
            for b in proto.get("batches", [])
            if b.get("has_data", False) and b.get("batch_id")
        ]
        if batch_ids:
            protocol_batch_map[protocol_id] = batch_ids

    if not protocol_batch_map:
        _flush_print("[WARN] no enabled protocol/batch with data for integration.")
        return

    _flush_print(f"[INFO] transcriptomics integration: experiment={experiment_id}")
    _flush_print(f"[INFO] output_root={output_root}")
    for pid, bids in sorted(protocol_batch_map.items()):
        _flush_print(f"[INFO]   {pid}: {bids}")

    consistency_policy = str(
        integration_cfg.get("reference_mismatch_policy", "fail")
    ).strip().lower()

    # ---- 联合模型（主路径）----
    # 批次内分组平衡时，联合 ~batch + condition 的统计效力高于事后 Meta：
    # 批次作为固定效应被吸收，离散度在全部样本上共享估计。
    # 设计不可识别的组合会被预检拒绝，自动退回下面的 Meta 路径。
    joint_cfg = (integration_cfg.get("joint_model", {}) or {})
    joint_summary: dict[str, Any] = {}
    if joint_cfg.get("enabled", True) and any(
        len(bids) >= 2 for bids in protocol_batch_map.values()
    ):
        _flush_print("\n[integration] joint cross-batch model (~batch + condition)")
        joint_summary = joint_model_mod.run_joint_model(experiment_id)
    elif joint_cfg.get("enabled", True):
        _flush_print("\n[integration] joint model skipped: no protocol has >= 2 batches.")

    within_result: dict[str, dict[tuple[str, str], pd.DataFrame]] = {}
    consistency_reports: dict[str, str] = {}
    for protocol_id, batch_ids in sorted(protocol_batch_map.items()):
        _flush_print(f"\n[integration] within protocol: {protocol_id}")

        # 只有真正要合并多个批次时才需要校验参考一致性
        if len(batch_ids) >= 2:
            ok, result, report = check_batch_consistency(
                output_root, protocol_id, batch_ids, policy=consistency_policy
            )
            consistency_reports[protocol_id] = report
            if result.get("consistent"):
                _flush_print(f"  [consistency] {report}")
            elif ok:
                _flush_print(f"  [consistency][WARN] 参考/技术参数不一致，"
                             f"但 reference_mismatch_policy='{consistency_policy}' 允许继续:\n{report}")
            else:
                _flush_print(f"  [consistency][BLOCKED] 拒绝整合 {protocol_id}:\n{report}")
                _flush_print("  修复参考版本/技术参数后重跑各批次，或将 "
                             "integration.reference_mismatch_policy 设为 'warn' 显式接受风险。")
                within_result[protocol_id] = {}
                continue

        within_result[protocol_id] = within_protocol_meta_analysis(
            output_root=output_root,
            protocol_id=protocol_id,
            batch_ids=batch_ids,
            integration_config=integration_cfg,
        )

    cross_result: dict[tuple[str, str], pd.DataFrame] = {}
    if (
        integration_cfg.get("allow_cross_protocol", False)
        and len(within_result.keys()) >= 2
    ):
        _flush_print("\n[integration] cross protocol consensus")
        cross_result = cross_protocol_consensus(
            output_root=output_root,
            protocol_ids=sorted(within_result.keys()),
            integration_config=integration_cfg,
        )

    report_path = generate_integration_report(
        output_root=output_root,
        within_result=within_result,
        cross_result=cross_result,
        experiment_id=experiment_id,
        consistency_reports=consistency_reports,
        joint_summary=joint_summary,
    )
    _flush_print(f"\n[DONE] integration complete: {report_path}")
