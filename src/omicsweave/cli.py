#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omicsweave 命令行入口。

用法::

    omicsweave [--project-dir DIR] init [--id 项目ID] [--species mus_musculus] [--experiment 实验ID]
    omicsweave [--project-dir DIR] validate [--experiment 实验ID ...]
    omicsweave [--project-dir DIR] <模块> <动作> [该动作自身的参数 ...]

示例::

    omicsweave --project-dir D:/lab_project init --id dio_mouse
    omicsweave --project-dir D:/lab_project validate
    omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio
    omicsweave mfa run --help          # 查看该动作自身的参数

``--project-dir`` 指定实验项目目录（包含 ``project.yaml``、``experiments/`` 与 ``analysis_results/``），
省略时见 :func:`omicsweave.core.project.project_dir` 的查找顺序。
模块动作之后的参数原样转交给对应模块的 ``main()``。
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

from omicsweave import __version__
from omicsweave.core.model import Species
from omicsweave.core.project import project_dir, set_project_dir

# 模块 -> 动作 -> (入口 "包.模块:函数", 简介)
COMMANDS: dict[str, dict[str, tuple[str, str]]] = {
    "mfa": {
        "run": ("omicsweave.mfa.run_mfa_analysis:main", "MFA 差异通量分析与绘图"),
        "plot": ("omicsweave.mfa.plot_mfa_results:main", "基于已生成表格重新绘图"),
        "report": ("omicsweave.mfa.generate_mfa_interpretation_report:main", "生成按组织解读报告"),
        "template": ("omicsweave.adapters.mfa.cli:template_main", "生成 MFA 求解结果的标准 Excel 模板"),
        "check": ("omicsweave.adapters.mfa.cli:check_main", "检查 MFA 结果文件，或检查实验的 MFA 输入"),
        "convert": ("omicsweave.adapters.mfa.cli:convert_main", "旧版自研软件导出格式 → 标准格式"),
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
        "run": ("omicsweave.integration.run_multiomics_analysis:main", "以 MFA 为核心的多组学整合"),
        "harmonize": ("omicsweave.integration.run_multiomics_harmonization:main", "样本级多组学协调"),
        "project": ("omicsweave.integration.run_multiomics_project_analysis:main", "多实验随机效应 Meta 与项目报告"),
        "cross-tissue": ("omicsweave.integration.run_cross_tissue_integration:main", "跨组织通路整合打分"),
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
        help="实验项目目录（包含 project.yaml、experiments/ 与 analysis_results/）",
    )
    modules = parser.add_subparsers(dest="module", metavar="<命令>", required=True)

    init = modules.add_parser("init", help="新建项目目录模板")
    init.add_argument("--id", dest="project_id", default=None, help="项目 ID，默认使用目录名")
    init.add_argument("--species", choices=[s.value for s in Species], default=Species.MOUSE.value)
    init.add_argument("--experiment", default="exp01", help="模板中第一个实验的 ID（默认 exp01）")
    init.add_argument("--force", action="store_true", help="覆盖已有的模板文件")

    validate = modules.add_parser("validate", help="校验项目配置、样本表与输入文件")
    validate.add_argument("--experiment", action="append", default=None,
                          help="只校验指定实验，可重复；默认校验全部实验")

    for module, actions in COMMANDS.items():
        module_parser = modules.add_parser(module, help=f"{module} 相关动作")
        action_parsers = module_parser.add_subparsers(dest="action", metavar="<动作>", required=True)
        for action, (_, summary) in actions.items():
            # 动作自身的参数（含 --help）全部转交给目标模块解析
            action_parsers.add_parser(action, help=summary, add_help=False)
    return parser


def _run_init(args: argparse.Namespace) -> int:
    from omicsweave.core.scaffold import create_project

    root = args.project_dir.resolve() if args.project_dir is not None else Path.cwd()
    project_id = args.project_id or root.name
    try:
        written = create_project(root, project_id, Species(args.species), args.experiment, force=args.force)
    except FileExistsError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    print(f"已在 {root} 生成项目模板：")
    for path in written:
        print(f"  {path.relative_to(root)}")
    print("下一步：填写 subjects.tsv / samples.tsv 并启用所需组学，然后运行 omicsweave validate")
    return 0


def _run_validate(args: argparse.Namespace) -> int:
    from omicsweave.core.validate import validate_project

    report = validate_project(project_dir(), args.experiment)
    print(report.format())
    return 0 if report.ok else 1


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args, forwarded = parser.parse_known_args(argv)
    if args.project_dir is not None:
        set_project_dir(args.project_dir)

    if args.module in ("init", "validate"):
        if forwarded:
            parser.error(f"无法识别的参数：{' '.join(forwarded)}")
        handler = _run_init if args.module == "init" else _run_validate
        sys.exit(handler(args))

    target, _ = COMMANDS[args.module][args.action]
    module_name, func_name = target.split(":")
    entry = getattr(importlib.import_module(module_name), func_name)

    sys.argv = [f"omicsweave {args.module} {args.action}", *forwarded]
    entry()


if __name__ == "__main__":
    main()
