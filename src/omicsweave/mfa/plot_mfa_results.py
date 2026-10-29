#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 结果可视化

从已保存的分析表格重新生成所有图表，支持在不重新运行分析的情况下单独执行。

Usage:
    python -m omicsweave.mfa.plot_mfa_results --experiment 20250425_cimr_wt_mice_dio
"""

import argparse
from pathlib import Path

from .mfa_experiment_config import load_mfa_config
from .plotting import generate_all_plots


def main():
    parser = argparse.ArgumentParser(description="MFA Plot Generation")
    parser.add_argument("--experiment", type=str, required=True,
                        help="Experiment ID, e.g. 20250425_cimr_wt_mice_dio")
    args = parser.parse_args()

    cfg = load_mfa_config(experiment_id=args.experiment)
    print(f"[INFO] Generating plots for: {cfg.experiment_id}")
    print(f"[INFO]   output: {cfg.output_dir}")

    print(f"[INFO]   conditions: {list(cfg.conditions.keys())}")

    condition_colors = cfg.condition_colors or {}
    tissue_colors = cfg.tissue_colors or {}

    generate_all_plots(
        output_dir=cfg.output_dir,
        condition_colors=condition_colors,
        tissue_colors=tissue_colors,
        comparisons=cfg.comparisons,
        conditions=[k.lower() for k in cfg.conditions.keys()],
        sample_conditions=cfg.sample_groups,
        analysis_config=cfg.analysis,
    )


if __name__ == "__main__":
    main()
