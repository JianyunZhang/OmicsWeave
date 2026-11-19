#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omicsweave 命令行入口。

用法::

    omicsweave [--project-dir DIR] init [--id 项目ID] [--species mus_musculus] [--experiment 实验ID]
    omicsweave [--project-dir DIR] validate [--experiment 实验ID ...]
    omicsweave [--project-dir DIR] migrate --experiment 实验ID [--dry-run] [--force]
    omicsweave [--project-dir DIR] run --experiment 实验ID [--omics mfa,...|all] [--steps qc,differential]
    omicsweave [--project-dir DIR] integrate --experiment 实验ID [--mode flux_anchored,pathway_convergence] [--source 运行ID]
    omicsweave [--project-dir DIR] export --experiment 实验ID [--run 运行ID] [--format xlsx|tsv] [--output 路径]
    omicsweave [--project-dir DIR] stats --experiment 实验ID [--source 运行ID]
    omicsweave [--project-dir DIR] refs list | build <genes|proteins|metabolites|reaction_roles> [...] | audit --experiment 实验ID
    omicsweave [--project-dir DIR] replicate [--set 复现组ID]
    omicsweave [--project-dir DIR] diff 运行A 运行B [--experiment 实验ID] [--output 目录]
    omicsweave [--project-dir DIR] pipeline [--experiment 实验ID ...] [--from run|integrate|scripts] [--source 运行ID] [--dry-run]
    omicsweave [--project-dir DIR] <模块> <动作> [该动作自身的参数 ...]

示例::

    omicsweave --project-dir D:/lab_project init --id dio_mouse
    omicsweave --project-dir D:/lab_project validate
    omicsweave --project-dir D:/lab_project migrate --experiment 20250425_cimr_wt_mice_dio
    omicsweave --project-dir D:/lab_project run --experiment 20250425_cimr_wt_mice_dio --omics mfa
    omicsweave --project-dir D:/lab_project pipeline        # validate → run → integrate → 课题脚本，一条命令
    omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio      # 旧版流程，读取旧版配置
    omicsweave mfa run --help          # 查看该动作自身的参数

``--project-dir`` 指定实验项目目录（包含 ``project.yaml``、``experiments/`` 与 ``analysis_results/``），
省略时见 :func:`omicsweave.core.project.project_dir` 的查找顺序。
模块动作之后的参数原样转交给对应模块的 ``main()``。
"""
from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

from omicsweave import __version__
from omicsweave.core.model import Species
from omicsweave.core.project import project_context, project_dir

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

    migrate = modules.add_parser("migrate", help="把旧版各组学配置转换为 experiment.yaml 与样本表")
    migrate.add_argument("--experiment", action="append", required=True,
                         help="experiments/ 下的实验 ID，可重复")
    migrate.add_argument("--dry-run", action="store_true", help="只打印迁移结果，不写文件")
    migrate.add_argument("--force", action="store_true", help="覆盖实验目录中已有的 experiment.yaml 与样本表")

    run = modules.add_parser("run", help="按 experiment.yaml 运行组学分析，生成结果包")
    run.add_argument("--experiment", required=True, help="experiments/ 下的实验 ID")
    run.add_argument("--omics", action="append", default=None,
                     help="要运行的组学，可重复或用逗号分隔；默认 all（全部已启用且已接入统一运行的组学）")
    run.add_argument("--steps", action="append", default=None,
                     help="只运行这些步骤，可重复或用逗号分隔，如 qc,differential；默认全部")

    integrate = modules.add_parser("integrate", help="读取统一结果表做多组学整合，结果写入新的结果包")
    integrate.add_argument("--experiment", required=True, help="experiments/ 下的实验 ID")
    integrate.add_argument("--mode", action="append", default=None,
                           help="整合模式，可重复或用逗号分隔；默认 experiment.yaml 中 integration.modes 列出的模式")
    integrate.add_argument("--source", default=None, help="源结果包的运行 ID，默认最近一次含统一结果表的运行")

    export = modules.add_parser("export", help="把结果包中的统一结果表导出为 Excel 或 TSV")
    export.add_argument("--experiment", required=True, help="experiments/ 下的实验 ID")
    export.add_argument("--run", default=None, help="运行 ID，默认最近一次")
    export.add_argument("--format", choices=["xlsx", "tsv"], default="xlsx", help="导出格式，默认 xlsx")
    export.add_argument("--output", default=None, help="输出文件（xlsx）或目录（tsv），默认写在结果包的 exports/ 下")

    refs = modules.add_parser("refs", help="知识表：列出（list）、生成（build <名称> [生成器参数]）、对照数据审计（audit）")
    refs.add_argument("action", choices=["list", "build", "audit"])
    refs.add_argument("name", nargs="?", default=None, help="build 时的知识表名称：genes / proteins / metabolites / reaction_roles")
    refs.add_argument("--experiment", default=None, help="audit：实验 ID")
    refs.add_argument("--run", default=None, help="audit：结果包运行 ID，默认最近一次含数据集的运行")
    refs.add_argument("--output", default=None, help="audit：把明细写到该 TSV 文件")

    replicate = modules.add_parser("replicate", help="按 project.yaml 的 replication 判断跨实验 / 跨队列的复现，"
                                                      "结果写入项目级结果包")
    replicate.add_argument("--set", action="append", default=None, help="只运行这些复现组，可重复；默认全部")

    stats = modules.add_parser("stats", help="对已有结果包重新运行后处理（公式模型、派生特征、可检出效应、稳健性），"
                                              "结果写入新的结果包")
    stats.add_argument("--experiment", required=True, help="experiments/ 下的实验 ID")
    stats.add_argument("--source", default=None, help="源结果包的运行 ID，默认最近一次含统一结果表的运行")

    maps = modules.add_parser("map", help="从已有的整合结果包重画分组织代谢通路总图（图 E），写到新目录，不改动原结果包")
    maps.add_argument("--experiment", required=True, help="experiments/ 下的实验 ID")
    maps.add_argument("--source", default=None, help="整合结果包的运行 ID，默认最近一次含反应档案评估表的运行")
    maps.add_argument("--layout", default=None,
                      help="布局文件（相对路径以实验目录为基准）；默认取 experiment.yaml 的 pathway_map.layout，再默认为包内布局")
    maps.add_argument("--tissue", action="append", default=None, help="只画这些组织，可重复")
    maps.add_argument("--contrast", action="append", default=None, help="只画这些比较，可重复")
    maps.add_argument("--format", action="append", default=None, help="图件格式，可重复，默认同 experiment.yaml（pdf、png）")
    maps.add_argument("--output", type=Path, default=None,
                      help="输出目录，默认 <结果根>/<实验>/_redraws/<运行 ID>/pathway_map_<时间>/")

    pipeline = modules.add_parser("pipeline", help="标准流程：validate → run → integrate → 课题脚本"
                                                    "（project.yaml 的 pipeline 一节），各步之间传递确切的运行 ID")
    pipeline.add_argument("--experiment", action="append", default=None,
                          help="只运行这些实验，可重复；默认 pipeline.experiments，未写时为全部实验")
    pipeline.add_argument("--omics", action="append", default=None, help="同 run --omics，默认 all")
    pipeline.add_argument("--from", dest="start", choices=["run", "integrate", "scripts"], default=None,
                          help="从哪一步开始：integrate / scripts 沿用已有的单组学结果包（--source，默认最近一次）；"
                               "默认 run，给出 --source 时默认 integrate")
    pipeline.add_argument("--source", default=None, help="沿用的单组学结果包运行 ID（只能选一个实验）")
    pipeline.add_argument("--results-dir", type=Path, help="覆盖本次结果根目录，相对路径以项目为基准")
    pipeline.add_argument("--integration-run", help="从 scripts 续跑时固定与 --source 对应的整合包")
    pipeline.add_argument("--dry-run", action="store_true", help="只做检查并列出各步命令，不运行")

    diff = modules.add_parser("diff", help="逐表对照两个结果包，报告是否逐值相同、差值分布与等级转移")
    diff.add_argument("run_a", help="结果包 A：运行 ID 或结果包目录")
    diff.add_argument("run_b", help="结果包 B：运行 ID 或结果包目录")
    diff.add_argument("--experiment", default=None, help="按运行 ID 查找时限定实验，默认在全部实验中查找")
    diff.add_argument("--output", default=None, help="把对照明细写到该目录（tables / columns / transitions.tsv）")
    diff.add_argument("--atol", type=float, default=1e-12, help="数值相同的绝对容差，默认 1e-12")
    diff.add_argument("--rtol", type=float, default=1e-9, help="数值相同的相对容差，默认 1e-9")

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


def _run_migrate(args: argparse.Namespace) -> int:
    from omicsweave.core.migrate import run

    return run(args.experiment, project_dir(), dry_run=args.dry_run, force=args.force)


def _split(values: list[str] | None) -> list[str]:
    return [part.strip() for value in values or [] for part in value.split(",") if part.strip()]


def _run_run(args: argparse.Namespace) -> int:
    from omicsweave.core.run import ALL, RunError, run_experiment

    try:
        run_dir, notes = run_experiment(args.experiment, _split(args.omics) or [ALL], _split(args.steps) or None,
                                        root=project_dir(), command=["omicsweave", *args.argv])
    except RunError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    for note in notes:
        print(f"[提示] {note}")
    print(f"结果包：{run_dir}")
    return 0


def _run_integrate(args: argparse.Namespace) -> int:
    from omicsweave.core.run import RunError
    from omicsweave.integration.integrate import integrate_experiment

    try:
        run_dir, notes = integrate_experiment(args.experiment, _split(args.mode), source=args.source,
                                              root=project_dir(), command=["omicsweave", *args.argv])
    except RunError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    for note in notes:
        print(f"[提示] {note}")
    print(f"结果包：{run_dir}")
    return 0


def _run_export(args: argparse.Namespace) -> int:
    from omicsweave.contract.dataset import ContractError
    from omicsweave.contract.export import export_bundle, find_run
    from omicsweave.core.config import ConfigError, load_project

    try:
        project = load_project(project_dir())
        run_dir = find_run(project.path("results"), args.experiment, args.run)
        path = export_bundle(run_dir, args.format, args.output)
    except (ConfigError, ContractError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    print(f"已导出：{path}")
    return 0


def _resolve_bundle(value: str, experiment: str | None) -> Path:
    """结果包目录，或在 ``<paths.results>/<实验>/`` 下按运行 ID 查找。"""
    from omicsweave.contract.dataset import ContractError
    from omicsweave.contract.manifest import MANIFEST_FILE
    from omicsweave.core.config import load_project

    path = Path(value)
    # 已有的目录直接对照（不要求运行清单，便于对照临时输出）
    if path.is_dir():
        return path
    results = load_project(project_dir()).path("results")
    pattern = f"{experiment}/{value}" if experiment else f"*/{value}"
    found = [p for p in results.glob(pattern) if (p / MANIFEST_FILE).is_file()]
    if len(found) != 1:
        problem = "找不到" if not found else f"在多个实验中都有（{'、'.join(p.parent.name for p in found)}），请用 --experiment 指定"
        raise ContractError(f"结果包 {value} {problem}")
    return found[0]


def _run_diff(args: argparse.Namespace) -> int:
    """0：逐值相同；1：有差异；2：无法对照。"""
    from omicsweave.contract.dataset import ContractError
    from omicsweave.contract.diff import diff_bundles, format_report, write_report
    from omicsweave.core.config import ConfigError

    try:
        a, b = (_resolve_bundle(v, args.experiment) for v in (args.run_a, args.run_b))
        result = diff_bundles(a, b, atol=args.atol, rtol=args.rtol)
    except (ConfigError, ContractError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2
    print(format_report(result))
    if args.output:
        print(f"明细：{write_report(result, args.output)}")
    return 0 if result.identical else 1


def _run_refs(args: argparse.Namespace) -> int:
    from omicsweave.core import refs

    if args.action == "list":
        print(refs.registry_status().drop(columns=["builder"]).to_string(index=False))
        return 0
    if args.action == "build":
        if not args.name:
            print("错误：refs build 需要知识表名称，如 omicsweave refs build metabolites", file=sys.stderr)
            return 1
        return refs.build(args.name, list(args.forwarded))
    from omicsweave.core.config import ConfigError, load_experiment, load_project

    try:
        experiment = load_experiment(load_project(project_dir()), args.experiment) if args.experiment else None
    except ConfigError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    run_dir = (experiment.project.path("results") / experiment.id / args.run) if experiment and args.run else None
    table = refs.audit(experiment, run_dir)
    for message in refs.audit_messages(table):
        print(f"[不一致] {message}")
    coverage = table[table["issue"] == "coverage"]
    for row in coverage.itertuples():
        print(f"[覆盖] {row.omics} {row.tissue}：反应映射的直接酶测到 {row.data_value}")
    if args.output:
        from omicsweave.contract.export import write_excel_tsv
        write_excel_tsv(table, Path(args.output))
        print(f"明细：{args.output}")
    return 0


def _run_replicate(args: argparse.Namespace) -> int:
    from omicsweave.core.config import ConfigError
    from omicsweave.integration.replication import ReplicationError, replicate_project

    try:
        run_dir, notes = replicate_project(project_dir(), _split(args.set) or None, command=["omicsweave", *args.argv])
    except (ConfigError, ReplicationError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    for note in notes:
        print(f"[提示] {note}")
    print(f"结果包：{run_dir}")
    return 0


def _run_stats(args: argparse.Namespace) -> int:
    from omicsweave.core.run import RunError
    from omicsweave.stats.post import stats_experiment

    try:
        run_dir, notes = stats_experiment(args.experiment, args.source, root=project_dir(),
                                          command=["omicsweave", *args.argv])
    except RunError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    for note in notes:
        print(f"[提示] {note}")
    print(f"结果包：{run_dir}")
    return 0


def _run_map(args: argparse.Namespace) -> int:
    from omicsweave.contract.dataset import ContractError
    from omicsweave.contract.export import find_run
    from omicsweave.core.config import ConfigError, load_experiment, load_project
    from omicsweave.integration.mappings.network_layout import LayoutError
    from omicsweave.integration.pathway_map import redraw

    try:
        project = load_project(project_dir())
        experiment = load_experiment(project, args.experiment)
        base = project.path("results") / args.experiment
        if args.source:
            bundle = find_run(project.path("results"), args.experiment, args.source)
        else:
            found = sorted(p.parents[3] for p in base.glob(
                "run_*/integration/reaction_dossier/tables/reaction_assessment.parquet"))
            if not found:
                raise ContractError(f"{base} 下没有含反应档案评估表的整合结果包")
            bundle = found[-1]
        configured = dict((experiment.config.integration.parameters.get("reaction_dossier") or {})
                          .get("pathway_map") or {})
        overrides = {"layout": args.layout, "tissues": args.tissue, "contrasts": args.contrast, "formats": args.format}
        configured.update({k: v for k, v in overrides.items() if v})
        configured["enabled"] = True
        out_dir, notes = redraw(bundle, args.output, configured, experiment.root)
    except (ConfigError, ContractError, FileNotFoundError, LayoutError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    for note in notes:
        print(f"[提示] {note}")
    print(f"输出目录：{out_dir}")
    return 0


def _run_pipeline(args: argparse.Namespace) -> int:
    from omicsweave.core.pipeline import run_pipeline
    from omicsweave.core.run import RunError

    start = args.start or ("integrate" if args.source else "run")
    try:
        results = run_pipeline(project_dir(), args.experiment, _split(args.omics) or None, start=start,
                               source=args.source, dry_run=args.dry_run, results_dir=args.results_dir,
                               integration_run=args.integration_run, command=["omicsweave", *args.argv])
    except RunError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    if args.dry_run:
        print("（只列出了各步，没有运行）")
        return 0
    for result in results:
        print(f"[{result.experiment}] 单组学结果包：{result.run_dir}")
        if result.integration_dir:
            print(f"[{result.experiment}] 整合结果包：{result.integration_dir}")
        if result.scripts:
            print(f"[{result.experiment}] 课题脚本：{'、'.join(result.scripts)}")
        print(f"[{result.experiment}] 结果摘要：{result.delivery_dir / 'summary.xlsx'}")
        print(f"[{result.experiment}] 交付索引：{result.delivery_dir / 'index.json'}")
    return 0


_HANDLERS = {"init": _run_init, "validate": _run_validate, "migrate": _run_migrate, "run": _run_run,
             "integrate": _run_integrate, "export": _run_export, "diff": _run_diff, "stats": _run_stats,
             "refs": _run_refs, "replicate": _run_replicate, "pipeline": _run_pipeline, "map": _run_map}


def project_pipeline_main(project: Path) -> None:
    """本地课题脚本入口，保持 IDE 与重定向日志中的 Unicode 可读。"""
    if os.name == "nt" and "PYTHONIOENCODING" not in os.environ:
        for stream in (sys.stdout, sys.stderr):
            if hasattr(stream, "reconfigure"):
                stream.reconfigure(encoding="utf-8")
        os.environ["PYTHONIOENCODING"] = "utf-8"
    main(["--project-dir", str(project), "pipeline", *sys.argv[1:]])


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args, forwarded = parser.parse_known_args(argv)
    args.argv = list(sys.argv[1:] if argv is None else argv)

    if args.module in _HANDLERS:
        # refs build 的其余参数交给生成器
        if args.module == "refs" and getattr(args, "action", None) == "build":
            args.forwarded = forwarded
        elif forwarded:
            parser.error(f"无法识别的参数：{' '.join(forwarded)}")
        handler = _HANDLERS[args.module]
        with project_context(args.project_dir):
            code = handler(args)
        sys.exit(code)

    target, _ = COMMANDS[args.module][args.action]
    module_name, func_name = target.split(":")
    # 在项目上下文内导入：个别旧模块在导入时就读取项目目录
    with project_context(args.project_dir):
        entry = getattr(importlib.import_module(module_name), func_name)
        sys.argv = [f"omicsweave {args.module} {args.action}", *forwarded]
        entry()


if __name__ == "__main__":
    main()
