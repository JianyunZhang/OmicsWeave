#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA Analysis Runner

配置驱动的 MFA 分析入口。通过实验 ID 自动加载
experiments/<id>/mfa_data/config.yaml + groups.yaml。

Usage:
    python run_mfa_analysis.py --experiment 20250425_cimr_wt_mice_dio
    python run_mfa_analysis.py --experiment 20250425_cimr_wt_mice_dio --differential-only
    python run_mfa_analysis.py --experiment 20250425_cimr_wt_mice_dio --no-plots
"""

import argparse
from pathlib import Path

from .mfa_analysis import MFAAnalyser

_RUN_MODES = {"default", "no_plots", "differential_only"}


def _normalize_comparisons(comparisons):
    """将比较对统一规范为小写 tuple 列表。"""
    if not comparisons:
        return None
    return [(str(a).lower(), str(b).lower()) for a, b in comparisons]


def _run_pipeline(analyser: MFAAnalyser, comparisons=None, mode: str = "default"):
    """按运行模式执行 MFA 流程。"""
    mode = mode.lower()
    if mode not in _RUN_MODES:
        raise ValueError(f"Unknown run mode: {mode}")

    comparisons = _normalize_comparisons(comparisons)

    if mode == "default":
        analyser.run_all(comparisons=comparisons)
        return

    if mode == "differential_only":
        print("=" * 70)
        print("MFA Differential-Only Mode")
        print("=" * 70)
        analyser.sample_flux_matrix()
        analyser.differential_flux(comparisons=comparisons)
        print("\n[DONE] differential_flux CSV 已输出，跳过其他分析与绘图。")
        return

    print("=" * 70)
    print("MFA No-Plots Mode")
    print("=" * 70)
    analyser.quality_control()
    analyser.infusion_steady_state_qc()
    analyser.write_provenance()
    analyser.sample_flux_matrix()
    analyser.flux_pca()
    analyser.differential_flux(comparisons=comparisons)
    analyser.flux_heatmap(comparisons=comparisons)
    analyser.pathway_flux_comparison(comparisons=comparisons)
    analyser.tissue_flux_overview()
    analyser.mid_analysis()
    analyser.mid_fit_residuals()
    print("\n[DONE] 数据分析与 CSV 输出完成，已跳过绘图与解读报告。")


def run_from_config(config: dict, mode: str = "default"):
    """Execute MFA analysis from a config dictionary (for programmatic use)."""
    output_dir = config.get("output_dir")
    if not output_dir:
        raise ValueError("config must include 'output_dir'")

    analyser = MFAAnalyser(
        output_dir=output_dir,
        condition_colors=config.get("condition_colors"),
        tissue_colors=config.get("tissue_colors"),
        pathway_groups=config.get("pathway_groups"),
        analysis_config=config.get("analysis"),
    )

    if config.get("flux_file"):
        analyser.load_flux_data(config["flux_file"])

    if config.get("mid_file"):
        analyser.load_mid_data(config["mid_file"])

    comparisons = config.get("comparisons")
    if comparisons:
        comparisons = [tuple(c) for c in comparisons]

    _run_pipeline(analyser, comparisons=comparisons, mode=mode)
    return analyser


def run_from_mfa_config(cfg, mode: str = "default"):
    """Execute MFA analysis from an MFAConfig object."""
    from .mfa_experiment_config import generate_comparisons

    analyser = MFAAnalyser(
        output_dir=str(cfg.output_dir),
        condition_colors=cfg.condition_colors,
        tissue_colors=cfg.tissue_colors,
        pathway_groups=cfg.pathway_groups,
        analysis_config=cfg.analysis,
        run_context={
            "experiment_id": cfg.experiment_id,
            "protocol_id": cfg.protocol_id,
            "batch_id": cfg.batch_id,
            "glycemia_file": cfg.glycemia_file,
            "infusion_qc": cfg.infusion_qc,
            "batch_metadata": cfg.batch_metadata,
            "tracer": cfg.tracer,
            "sample_groups": cfg.sample_groups,
            "mouse_groups": cfg.mouse_groups,
        },
    )

    if cfg.flux_file and cfg.flux_file.exists():
        analyser.load_flux_data(str(cfg.flux_file))
    else:
        if cfg.flux_file:
            print(f"[WARN] Flux file not found: {cfg.flux_file}")

    if cfg.mid_file and cfg.mid_file.exists():
        analyser.load_mid_data(str(cfg.mid_file))
    else:
        if cfg.mid_file:
            print(f"[WARN] MID file not found: {cfg.mid_file}, will skip MID analysis")

    comparisons = cfg.comparisons
    if not comparisons:
        comparisons = generate_comparisons(cfg.conditions)

    _run_pipeline(analyser, comparisons=comparisons, mode=mode)
    return analyser


# ══════════════════════════════════════════════════════════════════════════
# 多批次 / 多协议整合运行
# ══════════════════════════════════════════════════════════════════════════

def run_all_batches(experiment_id: str, mode: str = "default") -> list:
    """扫描 experiment.yaml，依次运行所有已启用协议的批次。

    Returns
    -------
    list of (protocol_id, batch_id, status) tuples.
    """
    from .mfa_experiment_config import (
        discover_experiment_protocols,
        load_mfa_config,
    )

    protocols = discover_experiment_protocols(experiment_id)
    results = []

    n_run = 0
    n_skip = 0

    for proto in protocols:
        protocol_id = proto.get("protocol_id")
        if not proto.get("enabled", False):
            for batch in proto.get("batches", []):
                results.append((protocol_id, batch.get("batch_id"), "disabled_protocol"))
                n_skip += 1
            continue

        for batch in proto.get("batches", []):
            batch_id = batch.get("batch_id")
            config_path = batch.get("config_path")
            has_data = batch.get("has_data", False)

            if not has_data or not config_path:
                results.append((protocol_id, batch_id, "no_data"))
                n_skip += 1
                print(f"  [skip] {protocol_id}/{batch_id}: 无数据或配置文件缺失")
                continue

            print(f"\n{'=' * 70}")
            print(f"  Running: {protocol_id}/{batch_id}")
            print(f"{'=' * 70}")

            try:
                cfg = load_mfa_config(config_path=str(config_path))
                run_from_mfa_config(cfg, mode=mode)
                results.append((protocol_id, batch_id, "success"))
                n_run += 1
            except Exception as e:
                print(f"  [ERROR] {protocol_id}/{batch_id}: {e}")
                results.append((protocol_id, batch_id, f"error: {e}"))

    print(f"\n{'=' * 70}")
    print(f"  Batch Runner Summary: {n_run} ran, {n_skip} skipped")
    for pid, bid, status in results:
        print(f"    {pid}/{bid}: {status}")
    print(f"{'=' * 70}")

    return results


def run_integration(experiment_id: str) -> None:
    """执行多批次/多协议整合分析（前提：各批次已完成单批分析）。

    整合流程：
    1. 检查 integration.enabled
    2. 发现所有协议/批次
    3. 同 protocol 内 meta 分析
    4. 跨 protocol 方向一致性（可选）
    5. 生成整合报告
    6. 生成整合绘图（可选，受 integration.plotting.enabled 控制）
    """
    from .mfa_experiment_config import (
        discover_experiment_protocols,
        load_integration_config,
        resolve_output_root,
    )
    from .mfa_integration import (
        cross_protocol_consensus,
        generate_integration_report,
        within_protocol_meta_analysis,
    )
    from .bridging_qc import bridging_sample_qc
    from .provenance import assert_batch_compatibility
    from .mfa_integration_plotting import generate_integration_plots

    integration_config = load_integration_config(experiment_id)
    if not integration_config.get("enabled", False):
        print("[INFO] integration.enabled = false，跳过整合分析。")
        print("       如需启用，请在 experiment.yaml 中设置 integration.enabled: true")
        return

    output_root = resolve_output_root(experiment_id)
    print(f"[INFO] 整合分析: experiment={experiment_id}")
    print(f"[INFO]   output_root: {output_root}")

    protocols = discover_experiment_protocols(experiment_id)

    # 按 protocol 分组有数据的 batch
    protocol_batch_map: dict[str, list[str]] = {}
    for proto in protocols:
        if not proto.get("enabled", False):
            continue
        pid = proto.get("protocol_id")
        if not pid:
            continue
        data_batches = [
            b["batch_id"] for b in proto.get("batches", [])
            if b.get("has_data", False) and b.get("batch_id")
        ]
        if data_batches:
            protocol_batch_map[pid] = data_batches

    if not protocol_batch_map:
        print("[WARN] 无任何 protocol/batch 有数据可整合。")
        return

    print(f"[INFO]   protocols with data: {list(protocol_batch_map.keys())}")
    for pid, bids in protocol_batch_map.items():
        print(f"[INFO]     {pid}: {bids}")

    # 同 protocol 内 Meta 分析
    all_protocol_ids = []
    for pid, batch_ids in sorted(protocol_batch_map.items()):
        # 桥接必须先跑：同一只动物被两个批次各拟合一次时，二者之差就是
        # 不含生物学信号的技术噪声地板，Meta 分析要用它判断合并效应
        # 是否真的超出技术重现性。
        # 合并之前先确认各批次的技术链（示踪剂、代谢网络、提取方法）一致。
        print(f"\n[integration] 批次技术链兼容性: {pid}")
        assert_batch_compatibility(
            output_root, pid, batch_ids,
            policy=integration_config.get("provenance_mismatch_policy", "warn"),
        )

        print(f"\n[integration] 桥接样本一致性检查: {pid}")
        bridging_sample_qc(output_root, pid, batch_ids, integration_config)

        print(f"\n[integration] 同 protocol 整合: {pid}")
        within_protocol_meta_analysis(output_root, pid, batch_ids, integration_config)
        all_protocol_ids.append(pid)

    # 跨 protocol 方向一致性（受 allow_cross_protocol 控制）
    allow_cross = bool(integration_config.get("allow_cross_protocol", False))
    if len(all_protocol_ids) >= 2 and allow_cross:
        print(f"\n[integration] 跨 protocol 方向一致性分析")
        cross_protocol_consensus(output_root, all_protocol_ids, integration_config)
    elif len(all_protocol_ids) >= 2:
        print(f"\n[integration] 有 {len(all_protocol_ids)} 个 protocol 有数据，"
              f"但 integration.allow_cross_protocol=false，跳过跨 protocol 分析")
        print("              不同示踪剂探测的是代谢网络的不同部分，通常互补而非重复；")
        print("              确认可比后再置为 true。")

    # 生成报告
    print(f"\n[integration] 生成整合报告")
    generate_integration_report(
        output_root, all_protocol_ids, integration_config, experiment_id,
    )

    # 可选：整合绘图
    plotting_cfg = integration_config.get("plotting", {}) or {}
    if plotting_cfg.get("enabled", False):
        print(f"\n[integration] 生成整合绘图")
        generate_integration_plots(output_root, integration_config)
    else:
        print("[integration] integration.plotting.enabled=false，跳过整合绘图")

    print("\n[DONE] 整合分析完成。")


def main():
    parser = argparse.ArgumentParser(description="MFA Analysis Pipeline")
    parser.add_argument("--experiment", type=str, required=True,
                        help="Experiment ID, e.g. 20250425_cimr_wt_mice_dio")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--report-only", action="store_true",
                            help="仅生成解读报告（不重新计算，从已有 CSV 重新生成）")
    mode_group.add_argument("--differential-only", action="store_true",
                            help="仅计算并输出 differential_flux CSV，跳过其他分析与绘图")
    mode_group.add_argument("--no-plots", action="store_true",
                            help="执行完整数据分析并输出 CSV，但跳过绘图和解读报告")
    parser.add_argument("--all-batches", action="store_true",
                        help="运行所有已启用协议的批次（需 experiment.yaml）")
    parser.add_argument("--with-plots", action="store_true",
                        help="配合 --all-batches：每个批次都生成图表与解读报告"
                             "（默认仅算表，绘图耗时较长）")
    parser.add_argument("--integrate", action="store_true",
                        help="执行多批次/多协议整合分析（前提：各批次已完成单批分析）")
    args = parser.parse_args()

    # ── 多批次模式 ──
    if args.all_batches:
        run_mode = "no_plots"
        if args.differential_only:
            run_mode = "differential_only"
        elif args.with_plots:
            run_mode = "default"
        run_all_batches(args.experiment, mode=run_mode)
        if args.integrate:
            run_integration(args.experiment)
        return

    # ── 仅整合模式 ──
    if args.integrate:
        run_integration(args.experiment)
        return

    # ── 单批次模式（默认行为，向后兼容） ──
    from .mfa_experiment_config import load_mfa_config
    cfg = load_mfa_config(experiment_id=args.experiment)
    print(f"[INFO] Loaded config for experiment: {cfg.experiment_id}")
    if cfg.protocol_id and cfg.batch_id:
        print(f"[INFO]   run:       {cfg.protocol_id}/{cfg.batch_id}")
    print(f"[INFO]   flux_file: {cfg.flux_file}")
    print(f"[INFO]   mid_file:  {cfg.mid_file}")
    print(f"[INFO]   output_root: {cfg.output_root}")
    print(f"[INFO]   output:    {cfg.output_dir}")

    if args.report_only:
        # 仅生成解读报告（Step 9），不重新运行分析计算
        from .generate_mfa_interpretation_report import generate_interpretation_report
        generate_interpretation_report(
            experiment_id=args.experiment,
            config_path=str(cfg.config_path),
        )
    else:
        run_mode = "default"
        if args.differential_only:
            run_mode = "differential_only"
        elif args.no_plots:
            run_mode = "no_plots"
        run_from_mfa_config(cfg, mode=run_mode)


if __name__ == "__main__":
    main()
