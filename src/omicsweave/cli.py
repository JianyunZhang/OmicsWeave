#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omicsweave 命令行入口。

用法::

    omicsweave [--project-dir DIR] <模块> <动作> [该动作自身的参数 ...]

示例::

    omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio
    omicsweave --project-dir D:/lab_project integration run --experiment <id>
    omicsweave mfa run --help          # 查看该动作自身的参数

``--project-dir`` 指定实验项目目录（包含 ``experiments/`` 与 ``analysis_results/``），
省略时见 :func:`omicsweave.common.project.project_dir` 的查找顺序。
动作之后的参数原样转交给对应模块的 ``main()``。
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

from omicsweave import __version__
from omicsweave.common.project import set_project_dir

# 模块 -> 动作 -> (入口 "包.模块:函数", 简介)
COMMANDS: dict[str, dict[str, tuple[str, str]]] = {
    "mfa": {
        "run": ("omicsweave.mfa.run_mfa_analysis:main", "MFA 差异通量分析与绘图"),
        "plot": ("omicsweave.mfa.plot_mfa_results:main", "基于已生成表格重新绘图"),
        "report": ("omicsweave.mfa.generate_mfa_interpretation_report:main", "生成按组织解读报告"),
    },
    "transcriptomics": {
        "run": ("omicsweave.transcriptomics.rnaseq_pipeline:main", "QC → 修剪 → kallisto 定量 → 差异表达"),
        "plot": ("omicsweave.transcriptomics.plot_results:main", "可视化、富集与 GSEA"),
        "report": ("omicsweave.transcriptomics.generate_interpretation_report:main", "生成解读报告"),
    },
    "metabolomics": {
        "run": ("omicsweave.metabolomics.run_metabolomics_analysis:main", "代谢组学完整流程"),
    },
    "proteomics": {
        "run": ("omicsweave.proteomics.run_proteomics_pipeline:main", "蛋白质组学主流程"),
        "plot": ("omicsweave.proteomics.plot_results:main", "基于已有结果重新绘图"),
        "report": ("omicsweave.proteomics.generate_interpretation_report:main", "生成组织级报告"),
    },
    "integration": {
        "run": ("omicsweave.multiomics.run_multiomics_analysis:main", "以 MFA 为核心的多组学整合"),
        "harmonize": ("omicsweave.multiomics.run_multiomics_harmonization:main", "样本级多组学协调"),
        "project": ("omicsweave.multiomics.run_multiomics_project_analysis:main", "多实验随机效应 Meta 与项目报告"),
        "cross-tissue": ("omicsweave.multiomics.run_cross_tissue_integration:main", "跨组织通路整合打分"),
    },
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="omicsweave",
        description="OmicsWeave 多组学分析工具包",
        epilog="在动作后加 --help 查看该动作自身的参数，例如：omicsweave mfa run --help",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--project-dir", type=Path, default=None,
        help="实验项目目录（包含 experiments/ 与 analysis_results/）",
    )
    modules = parser.add_subparsers(dest="module", metavar="<模块>", required=True)
    for module, actions in COMMANDS.items():
        module_parser = modules.add_parser(module, help=f"{module} 相关动作")
        action_parsers = module_parser.add_subparsers(dest="action", metavar="<动作>", required=True)
        for action, (_, summary) in actions.items():
            # 动作自身的参数（含 --help）全部转交给目标模块解析
            action_parsers.add_parser(action, help=summary, add_help=False)
    return parser


def main(argv: list[str] | None = None) -> None:
    args, forwarded = build_parser().parse_known_args(argv)
    if args.project_dir is not None:
        set_project_dir(args.project_dir)

    target, _ = COMMANDS[args.module][args.action]
    module_name, func_name = target.split(":")
    entry = getattr(importlib.import_module(module_name), func_name)

    sys.argv = [f"omicsweave {args.module} {args.action}", *forwarded]
    entry()


if __name__ == "__main__":
    main()
