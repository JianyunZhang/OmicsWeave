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
    omicsweave [--project-dir DIR] pipeline [--experiment 实验ID ...] [--from run|integrate|scripts] [--source 运行ID] [--tag 标签] [--pin] [--dry-run]
    omicsweave [--project-dir DIR] latest [--experiment 实验ID ...] [--target dir|summary|index] [--path] [--json] [--open] [--list]
    omicsweave [--project-dir DIR] clean [--experiment 实验ID ...] [--keep-last N] [--prune-bundles] [--dry-run | --apply]
    omicsweave [--project-dir DIR] <模块> <动作> [该动作自身的参数 ...]

示例::

    omicsweave --project-dir D:/lab_project init --id dio_mouse
    omicsweave --project-dir D:/lab_project validate
    omicsweave --project-dir D:/lab_project migrate --experiment 20250425_cimr_wt_mice_dio
    omicsweave --project-dir D:/lab_project run --experiment 20250425_cimr_wt_mice_dio --omics mfa
    omicsweave --project-dir D:/lab_project pipeline        # validate → run → integrate → 课题脚本，一条命令
    omicsweave mfa run --experiment 20250425_cimr_wt_mice_dio      # 旧版流程，读取旧版配置
    omicsweave mfa run --help          # 查看该动作自身的参数

``--project-dir`` 指定实验项目目录（包含 ``project.yaml``、``experiments/`` 与 ``results/``），
省略时见 :func:`omicsweave.core.project.project_dir` 的查找顺序。
模块动作之后的参数原样转交给对应模块的 ``main()``。
"""
from __future__ import annotations

import argparse
import contextvars
import importlib
import json
import os
import sys
from contextlib import contextmanager, redirect_stdout
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING

from omicsweave import __version__
from omicsweave.core.model import Species
from omicsweave.core.project import project_context, project_dir
from omicsweave.service.schemas import ErrorCode, ServiceError

if TYPE_CHECKING:
    from omicsweave.core.retention import RetentionPolicy

# 入口脚本（如 run_analysis.py）透传的保留策略上下文，默认均为 None（未启用）
_ENTRY_RETENTION: contextvars.ContextVar["RetentionPolicy | None"] = contextvars.ContextVar(
    "omicsweave_entry_retention", default=None,
)
_ENTRY_CLEAN_DEFAULTS: contextvars.ContextVar["RetentionPolicy | None"] = contextvars.ContextVar(
    "omicsweave_entry_clean_defaults", default=None,
)

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


class _ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        if _JSON_PARSE.get():
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, message)
        super().error(message)


_JSON_PARSE = contextvars.ContextVar("omicsweave_json_parse", default=False)


def _json_options(parser, *, use_format=True, alias=True):
    if use_format:
        parser.add_argument("--format", dest="response_format", choices=["text", "json"],
                            default=argparse.SUPPRESS, help="CLI 响应格式")
    if alias:
        parser.add_argument("--json", dest="response_format", action="store_const", const="json",
                            default=argparse.SUPPRESS, help="输出结构化 JSON 响应")


def build_parser() -> argparse.ArgumentParser:
    parser = _ArgumentParser(
        prog="omicsweave",
        description="OmicsWeave 多组学分析工具包",
        epilog="在动作后加 --help 查看该动作自身的参数，例如：omicsweave mfa run --help",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    _json_options(parser)
    parser.add_argument(
        "--project-dir", type=Path, default=None,
        help="实验项目目录（包含 project.yaml、experiments/ 与 results/）",
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
    run_gsea_group = run.add_mutually_exclusive_group()
    run_gsea_group.add_argument(
        "--gsea", dest="run_gsea", action="store_true", default=None,
        help="启用预排序 GSEA（覆盖配置中的默认关闭设置）",
    )
    run_gsea_group.add_argument(
        "--no-gsea", dest="run_gsea", action="store_false",
        help="跳过耗时的预排序 GSEA，仅保留快速本地 ORA（默认已关闭 GSEA，此项可强制覆盖配置中的 run_gsea: true）",
    )

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
    maps.add_argument("--mode", choices=["atlas", "publication"], default=None,
                      help="出图模式：atlas（完整诊断总图，默认）或 publication（双栏聚焦出版图）")
    maps.add_argument("--style", choices=["metro_atlas", "economist", "nature_editorial", "slide_dark"], default=None,
                      help="视觉主题预设：metro_atlas（默认，Beck 学术出版风）/ economist（经济学人）/ nature_editorial / slide_dark")
    maps.add_argument("--journal", choices=["cell_press", "nature"], default=None,
                      help="期刊出版规格：cell_press（默认，85/114/174 mm + Arial ≥6pt + 净版子图 + 1200×1200 方形 GA）/ nature（89/120/183 mm）")
    maps.add_argument("--color-mode", choices=["semantic", "pathway"], default=None,
                      help="轨道着色语义：semantic（默认，正文统一红蓝涨跌色）/ pathway（低饱和莫兰迪通路分区色）")
    maps.add_argument("--strip-mode", choices=["omics", "axes"], default=None,
                      help="四格微条内容：omics（默认，Flux / Protein / mRNA / Metab 四层组学图标）/ axes（S3 四机制轴）")
    maps.add_argument("--omics-edge-color", choices=["tinted", "grey"], default=None,
                      help="无方向性 MFA 通量的 Enz-pred (L1) 反应边配色：tinted（默认，饱和度封顶的空心柔彩红/蓝双轨）/ grey（保守中性灰）")
    maps.add_argument("--region", action="append", default=None,
                      help="publication 模式下显式指定保留的通路区域 ID，可重复；省略时自动聚焦活跃区域")
    maps.add_argument("--submap", action="append", default=None,
                      help="指定生成的正文代谢子通路专题图 ID（如 central_carbon / fatty_acid_ketone / glycerolipid / "
                           "amino_acid / urea_creatine / nucleotide / redox_one_carbon / all），可重复")
    maps.add_argument("--no-submaps", action="store_true", default=False,
                      help="本次重画不生成分专题子通路图")
    maps.add_argument("--no-abstract", action="store_true", default=False,
                      help="本次重画不同步生成 Graphical Abstract")
    maps.add_argument("--ga-layout", choices=["cell_press_square", "widescreen"], default=None,
                      help="Graphical Abstract 版式：cell_press_square（默认，1200×1200 px 方形）/ widescreen（宽屏）")
    maps.add_argument("--format", action="append", default=None, help="图件格式（pdf / png / html / tiff / eps / svg），可重复，默认 pdf、png、html")
    maps.add_argument("--output", type=Path, default=None,
                      help="输出目录，默认 <结果根>/<实验>/_redraws/<运行 ID>/pathway_map_<时间>/")

    abstract = modules.add_parser(
        "abstract",
        help="生成小规模 Graphical Abstract 概念摘要图（默认 Cell Press 1200×1200 方形，支持地铁 / 经济学人 / 手绘 / 像素 / nature_editorial）",
    )
    abstract.add_argument("--spec", type=Path, default=None,
                          help="跨项目通用的自定义 YAML 摘要规格文件（提供时不依赖实验整合结果包）")
    abstract.add_argument("--experiment", default=None, help="从已有整合结果包自动聚合时的实验 ID")
    abstract.add_argument("--source", default=None, help="整合结果包的运行 ID，默认最近一次含反应档案评估表的运行")
    abstract.add_argument("--layout", default="", help="可选布局文件路径（仅配合 --experiment）")
    abstract.add_argument("--tissue", action="append", default=None, help="只生成这些组织，可重复")
    abstract.add_argument("--contrast", action="append", default=None, help="只生成这些比较，可重复")
    abstract.add_argument("--style", choices=["metro_atlas", "economist", "editorial_ink", "pixel_abstract", "nature_editorial"],
                          default="metro_atlas",
                          help="概念摘要渲染风格：metro_atlas（代谢地铁，默认）/ economist（经济学人）/ editorial_ink（手绘）/ pixel_abstract（像素）/ nature_editorial")
    abstract.add_argument("--journal", choices=["cell_press", "nature"], default="cell_press",
                          help="期刊出版规格：cell_press（默认，1200×1200 方形 + Arial 12–16 pt）/ nature（宽屏）")
    abstract.add_argument("--ga-layout", choices=["cell_press_square", "widescreen"], default=None,
                          help="显式指定版式：cell_press_square（1200×1200 方形）/ widescreen（16:10 宽屏）；默认跟随 --journal")
    abstract.add_argument("--color-mode", choices=["semantic", "pathway"], default="semantic",
                          help="metro_atlas 轨道着色语义：semantic（默认，正文统一红蓝方向色）/ pathway（低饱和莫兰迪通路分区色）")
    abstract.add_argument("--format", action="append", default=None, help="输出格式（pdf / png / tiff / jpg / eps / svg），可重复，默认 pdf,png")
    abstract.add_argument("--output", type=Path, default=None, help="输出目录或不含扩展名的输出路径前缀")

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
    pipeline.add_argument("--tag", default=None, help="给本次交付目录附加人类可读标签（如 v1、manuscript）")
    pipeline.add_argument("--pin", action="store_true", help="锁定本次交付（写入 .keep 与 pinned 标记），避免被保留策略清理")
    pipeline_gsea_group = pipeline.add_mutually_exclusive_group()
    pipeline_gsea_group.add_argument(
        "--gsea", dest="run_gsea", action="store_true", default=None,
        help="启用预排序 GSEA（覆盖配置中的默认关闭设置）",
    )
    pipeline_gsea_group.add_argument(
        "--no-gsea", dest="run_gsea", action="store_false",
        help="跳过耗时的预排序 GSEA，仅保留快速本地 ORA（默认已关闭 GSEA，此项可强制覆盖配置中的 run_gsea: true）",
    )
    pipeline.add_argument("--dry-run", action="store_true", help="只做检查并列出各步命令，不运行")

    latest = modules.add_parser("latest", help="通过 results/latest.json 指针查询、列出或打开最新一次成功交付")
    latest.add_argument("--experiment", action="append", default=None,
                        help="指定实验 ID，可重复；默认查询结果目录下已有交付的全部实验")
    latest.add_argument("--results-dir", type=Path, default=None, help="覆盖结果根目录，相对路径以项目为基准")
    latest.add_argument("--target", choices=["dir", "summary", "index"], default=None,
                        help="定位目标：dir（交付目录）、summary（summary.xlsx）或 index（index.json）")
    latest.add_argument("--path", action="store_true", help="仅输出解析后的绝对路径（便于脚本或命令行串联）")
    latest.add_argument("--json", dest="as_json", action="store_true", help="以 JSON 格式输出最新交付指针与路径")
    latest.add_argument("--open", dest="open_target", action="store_true",
                        help="用系统默认程序打开目标文件或目录（未指定 --target 时默认打开 summary.xlsx）")
    latest.add_argument("--list", dest="list_all", action="store_true",
                        help="按时间倒序列出全部交付目录及状态、大小、标签与锁定标记")

    clean = modules.add_parser("clean", help="按保留策略预览或清理历史交付目录（及可选的孤立结果包），默认仅预览")
    clean.add_argument("--experiment", action="append", default=None,
                       help="只清理这些实验，可重复；默认处理结果目录下已有交付的全部实验")
    clean.add_argument("--results-dir", type=Path, default=None, help="覆盖结果根目录，相对路径以项目为基准")
    keep_group = clean.add_mutually_exclusive_group()
    keep_group.add_argument("--keep-last", type=int, default=None,
                            help="每个实验保留最近 N 次成功交付（>= 1）；默认读取入口/配置（缺省 5）")
    keep_group.add_argument("--keep-all", action="store_true",
                            help="不按数量淘汰成功交付，仅清理失败交付或孤立结果包")
    failed_group = clean.add_mutually_exclusive_group()
    failed_group.add_argument("--drop-failed", action="store_true", default=False,
                              help="清理失败及超时未完成的交付目录")
    failed_group.add_argument("--keep-failed", action="store_true", default=False,
                              help="保留失败与中断的交付目录")
    bundle_group = clean.add_mutually_exclusive_group()
    bundle_group.add_argument("--prune-bundles", action="store_true", default=False,
                              help="联动清理不再被任何保留交付引用的孤立 run_* 结果包")
    bundle_group.add_argument("--no-prune-bundles", action="store_true", default=False,
                              help="不删除 run_* 结果包（覆盖入口或配置默认值）")
    clean.add_argument("--stale-running-hours", type=float, default=None,
                       help="超过该小时数仍处于 running 状态的交付视为已中断（默认 24）")
    mode_group = clean.add_mutually_exclusive_group()
    mode_group.add_argument("--dry-run", action="store_true", default=False,
                            help="仅预览待清理目录与预计释放空间，不删除任何文件（默认）")
    mode_group.add_argument("--apply", action="store_true", default=False,
                            help="确认执行实际删除操作")

    diff = modules.add_parser("diff", help="逐表对照两个结果包，报告是否逐值相同、差值分布与等级转移")
    diff.add_argument("run_a", help="结果包 A：运行 ID 或结果包目录")
    diff.add_argument("run_b", help="结果包 B：运行 ID 或结果包目录")
    diff.add_argument("--experiment", default=None, help="按运行 ID 查找时限定实验，默认在全部实验中查找")
    diff.add_argument("--output", default=None, help="把对照明细写到该目录（tables / columns / transitions.tsv）")
    diff.add_argument("--atol", type=float, default=1e-12, help="数值相同的绝对容差，默认 1e-12")
    diff.add_argument("--rtol", type=float, default=1e-9, help="数值相同的相对容差，默认 1e-9")
    diff.add_argument("--results-dir", type=Path, help="正式结果根目录覆盖项")
    diff.add_argument("--trial-a", help="A 侧试验 ID")
    diff.add_argument("--trial-b", help="B 侧试验 ID")

    for subparser in (run, integrate, pipeline):
        subparser.add_argument("--background", action="store_true", help="派发独立后台作业")
        subparser.add_argument("--override", action="append", help="隔离试验覆盖 key.path=yaml_value，可重复")
        subparser.add_argument("--override-file", type=Path, help="隔离试验的 YAML/JSON 参数表")
        subparser.add_argument("--intent", choices=["sensitivity_check", "post_hoc_exploration", "qc_refinement",
                                                    "visualization_tune", "baseline", "manual"])
        subparser.add_argument("--rationale", help="隔离试验的依据")
        subparser.add_argument("--parent-run", "--parent-run-id", dest="parent_run", help="对照的正式运行 ID")
    for subparser in (run, integrate):
        subparser.add_argument("--dry-run", action="store_true", help="预览运行选择及试验参数")
    for subparser in (validate, run, integrate, pipeline, clean, diff):
        _json_options(subparser)
    _json_options(latest, alias=False)
    for subparser in (maps, abstract):
        _json_options(subparser, use_format=False)

    status = modules.add_parser("status", help="读取后台作业状态")
    selector = status.add_mutually_exclusive_group()
    selector.add_argument("--job")
    selector.add_argument("--experiment")
    selector.add_argument("--active", action="store_true")
    status.add_argument("--limit", type=int, default=20)

    logs = modules.add_parser("logs", help="读取后台作业日志与结构化事件")
    logs.add_argument("job")
    logs.add_argument("--tail", type=int, default=50)
    logs.add_argument("--since-event", type=int, default=0)
    logs.add_argument("--events-only", action="store_true")

    cancel = modules.add_parser("cancel", help="在下一安全检查点停止后台作业")
    cancel.add_argument("job")

    inspect = modules.add_parser("inspect", help="读取结果包或交付摘要")
    query = modules.add_parser("query", help="读取结果表的有界 Parquet 切片")
    trace = modules.add_parser("trace", help="按确切标识符读取关联结果")
    diagnose = modules.add_parser("diagnose", help="按固定结果包汇总科研护栏、质控和校准")
    for subparser in (inspect, query, trace, diagnose):
        subparser.add_argument("--experiment", required=True)
        subparser.add_argument("--trial", help="显式选择隔离试验")
        selection = subparser.add_mutually_exclusive_group()
        selection.add_argument("--run", default="latest")
        selection.add_argument("--delivery")
    query.add_argument("--table", required=True)
    layer = query.add_mutually_exclusive_group()
    layer.add_argument("--omics")
    layer.add_argument("--mode")
    query.add_argument("--columns", action="append", help="列名，逗号分隔或重复指定")
    query.add_argument("--filters", help="JSON 条件数组，每项包含 column、op、value")
    query.add_argument("--sort", "--sort-by", dest="sort_by")
    query.add_argument("--descending", action="store_true")
    query.add_argument("--limit", type=int, default=50)
    query.add_argument("--offset", type=int, default=0)
    trace.add_argument("--entity", required=True)
    trace.add_argument("--tissue")
    trace.add_argument("--contrast", dest="contrast_id")
    trace.add_argument("--limit", type=int, default=30)
    diagnose.add_argument("--limit", type=int, default=100, help="返回诊断条数，最多 200")
    for subparser in (status, logs, cancel, inspect, query, trace, diagnose):
        subparser.add_argument("--results-dir", type=Path, help="结果根目录覆盖项，相对路径以项目为基准")
        _json_options(subparser)

    schema = modules.add_parser("schema", help="查看已注册的参数和结果 JSON Schema")
    schema.add_argument("target", nargs="?", default="targets", help="如 omics/mfa 或 integration/reversal")
    inputs = modules.add_parser("inspect-input", help="查看实验输入文件元信息")
    inputs.add_argument("paths", nargs="*")
    inputs.add_argument("--experiment")
    config = modules.add_parser("config", help="预览或提交实验配置补丁")
    config.add_argument("action", choices=["patch"])
    config.add_argument("--target", choices=["experiment_yaml", "panels_yaml"], default="experiment_yaml")
    patch = config.add_mutually_exclusive_group(required=True)
    patch.add_argument("--patch", help="JSON 合并补丁")
    patch.add_argument("--patch-file", type=Path, help="包含 JSON 合并补丁的 UTF-8 文件")
    qc = modules.add_parser("qc", help="预览或提交样本 QC 标记")
    qc.add_argument("action", choices=["flag"])
    qc.add_argument("--updates", required=True, help="JSON 更新数组，每项包含 sample_id、qc_flag 和 qc_reason")
    for edit in (config, qc):
        edit.add_argument("--experiment", required=True)
        mode = edit.add_mutually_exclusive_group()
        mode.add_argument("--dry-run", action="store_true", help="预览变更（默认）")
        mode.add_argument("--apply", action="store_true", help="提交预览变更")
        edit.add_argument("--reason", default="")
        edit.add_argument("--expected-sha256", help="预览响应中的 sha256")
    for subparser in (schema, inputs, config, qc):
        _json_options(subparser)
    mcp = modules.add_parser("mcp", help="以 stdio 提供本地 MCP 服务")
    mcp.add_argument("action", choices=["serve"])
    mcp.add_argument("--allowed-data-root", type=Path, action="append", default=[])
    mcp.add_argument("--allow-writes", action="store_true", help="授权配置、样本表和项目模板工具写入")
    mcp.add_argument("--allow-custom-scripts", action="store_true", help="信任并运行项目的课题脚本")

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
    from omicsweave.service.project import validation, validation_summary

    report = validation(project_dir(), args.experiment)
    if args.machine:
        issues = [asdict(issue) for issue in report.issues]
        args.warnings.extend({k: row[k] for k in ("code", "where", "message")}
                             for row in issues if row["level"] == "warning")
        if not report.ok:
            code = (ErrorCode.CONFIG_SCHEMA_ERROR if any(i.code == "CONFIG_SCHEMA_ERROR" for i in report.errors)
                    else ErrorCode.VALIDATION_FAILED)
            raise ServiceError(code, "项目校验未通过", details={"issues": issues})
        return {key: value for key, value in validation_summary(report).items() if key != "ok"}
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
                                        root=project_dir(), command=["omicsweave", *args.argv],
                                        run_gsea=args.run_gsea)
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


def _resolve_bundle(value: str, experiment: str | None, results_dir: Path | None = None,
                    trial_id: str | None = None) -> Path:
    """结果包目录，或在 ``<paths.results>/<实验>/`` 下按运行 ID 查找。"""
    from omicsweave.contract.dataset import ContractError
    from omicsweave.contract.manifest import MANIFEST_FILE
    path = Path(value)
    # 已有的目录直接对照（不要求运行清单，便于对照临时输出）
    if path.is_dir() and not trial_id:
        return path
    if trial_id:
        from omicsweave.core.overrides import trial_root
        from omicsweave.core.path_scope import identifier

        identifier(value)
    _, results = _resolve_results_root(results_dir)
    if trial_id:
        results = trial_root(results, trial_id)
        if experiment:
            identifier(experiment)
    pattern = f"{experiment}/{value}" if experiment else f"*/{value}"
    found = [p for p in results.glob(pattern) if (p / MANIFEST_FILE).is_file()]
    if len(found) != 1:
        problem = "找不到" if not found else f"在多个实验中都有（{'、'.join(p.parent.name for p in found)}），请用 --experiment 指定"
        if not found:
            raise FileNotFoundError(f"结果包 {value} {problem}")
        raise ContractError(f"结果包 {value} {problem}")
    return found[0]


def _run_diff(args: argparse.Namespace) -> int:
    """0：逐值相同；1：有差异；2：无法对照。"""
    from omicsweave.contract.dataset import ContractError
    from omicsweave.contract.diff import diff_bundles, format_report, write_report
    from omicsweave.core.config import ConfigError
    from omicsweave.core.overrides import OverrideError
    from omicsweave.core.path_scope import ScopeError

    try:
        a, b = (_resolve_bundle(v, args.experiment, args.results_dir, trial)
                for v, trial in ((args.run_a, args.trial_a), (args.run_b, args.trial_b)))
        result = diff_bundles(a, b, atol=args.atol, rtol=args.rtol)
    except (ConfigError, ContractError, FileNotFoundError, OverrideError, ScopeError) as exc:
        if args.machine:
            raise
        print(f"错误：{exc}", file=sys.stderr)
        return 2
    if args.machine:
        from omicsweave.service.query import diff_payload

        return diff_payload(result, args.output)
    print(format_report(result))
    from omicsweave.service.query import configuration_diff

    if report := configuration_diff(a, b)["config_diff"]:
        print(f"生效配置差异：{len(report)} 项")
        for item in report[:20]:
            print(f"  {item['op']} {item['path']}")
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
    from omicsweave.core.config import ConfigError
    from omicsweave.integration.mappings.network_layout import LayoutError
    from omicsweave.service.visualization import render_map

    overrides = {
        "layout": args.layout, "tissues": args.tissue, "contrasts": args.contrast, "formats": args.format,
        "mode": args.mode, "style": args.style, "journal": args.journal, "color_mode": args.color_mode,
        "strip_mode": args.strip_mode, "omics_edge_color": args.omics_edge_color, "regions": args.region,
    }
    settings = {key: value for key, value in overrides.items() if value}
    if args.no_submaps:
        settings["submaps"] = {"enabled": False}
    elif args.submap:
        settings["submaps"] = {"enabled": True, "panels": _split(args.submap), "min_signal": "all"}
    if args.no_abstract:
        settings["abstract"] = {"enabled": False}
    elif args.ga_layout:
        settings["abstract"] = {"enabled": True, "layout": args.ga_layout}
    try:
        data = render_map(project_dir(), args.experiment, args.source, args.output, settings)
    except (ConfigError, ContractError, FileNotFoundError, LayoutError, ValueError) as exc:
        if args.machine:
            raise
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    if args.machine:
        return data
    for note in data["notes"]:
        print(f"[提示] {note}")
    print(f"输出目录：{data['output_dir']}")
    return 0


def _run_abstract(args: argparse.Namespace) -> int:
    from omicsweave.contract.dataset import ContractError
    from omicsweave.core.config import ConfigError
    from omicsweave.integration.mappings.network_layout import LayoutError
    from omicsweave.service.visualization import render_abstract

    try:
        data = render_abstract(project_dir(), experiment=args.experiment, source=args.source, spec_path=args.spec,
                               output=args.output, style=args.style, color_mode=args.color_mode, journal=args.journal,
                               layout=args.ga_layout, formats=args.format, tissues=args.tissue,
                               contrasts=args.contrast, layout_path=args.layout)
    except (ConfigError, ContractError, FileNotFoundError, LayoutError, ValueError, ServiceError) as exc:
        if args.machine:
            raise
        print(f"错误：{exc}", file=sys.stderr)
        return 1
    if args.machine:
        return data
    if args.spec is not None:
        for path in data["files"]:
            print(f"已写出：{path}")
    else:
        print(f"Graphical Abstract（{args.style}）共写出 {len(data['files'])} 个文件，目录：{data['output_dir']}")
    return 0


def _run_pipeline(args: argparse.Namespace) -> int:
    from omicsweave.core.pipeline import run_pipeline
    from omicsweave.core.run import RunError

    start = args.start or ("integrate" if args.source else "run")
    try:
        results = run_pipeline(project_dir(), args.experiment, _split(args.omics) or None, start=start,
                               source=args.source, dry_run=args.dry_run, results_dir=args.results_dir,
                               integration_run=args.integration_run, tag=args.tag, pin=args.pin,
                               retention=_ENTRY_RETENTION.get(), run_gsea=args.run_gsea,
                               command=["omicsweave", *args.argv])
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
        if result.latest_pointer:
            print(f"[{result.experiment}] 最新指针：{result.latest_pointer}")
    return 0


def _resolve_results_root(results_dir: Path | None):
    from omicsweave.core.config import load_project

    project = load_project(project_dir())
    root = project.resolve(results_dir).resolve() if results_dir is not None else project.path("results")
    return project, root


def _run_latest(args: argparse.Namespace) -> int:
    from omicsweave.core.config import ConfigError
    from omicsweave.core.delivery import (
        build_latest_payload,
        delivery_experiments,
        latest_pointer_path,
        list_deliveries,
        open_path,
        read_latest_pointer,
        resolve_latest_delivery,
    )
    from omicsweave.core.retention import directory_size, format_bytes

    try:
        project, results_root = _resolve_results_root(args.results_dir)
    except ConfigError as exc:
        if args.machine:
            raise
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    experiments = (
        list(dict.fromkeys(args.experiment))
        if args.experiment
        else (delivery_experiments(results_root) or list(project.config.pipeline.experiments))
    )
    if not experiments:
        if args.machine:
            raise ServiceError(ErrorCode.INPUT_NOT_FOUND, f"{results_root} 下没有找到任何交付目录")
        print(f"错误：{results_root} 下没有找到任何交付目录", file=sys.stderr)
        return 1

    if args.list_all:
        rows = []
        for exp_id in experiments:
            records = list(reversed(list_deliveries(results_root, exp_id)))
            latest_rec, _ = resolve_latest_delivery(results_root, exp_id)
            latest_id = latest_rec.delivery_id if latest_rec else None
            if args.as_json or args.machine:
                for rec in records:
                    rows.append({
                        **build_latest_payload(rec),
                        "status": rec.status,
                        "is_latest": rec.delivery_id == latest_id,
                        "directory": str(rec.directory),
                        "size_bytes": directory_size(rec.directory),
                    })
                continue
            print(f"[{exp_id}] 共 {len(records)} 个交付目录（{results_root / exp_id}）：")
            if not records:
                print("  （无）")
                continue
            for rec in records:
                flags = []
                if rec.delivery_id == latest_id:
                    flags.append("latest")
                if rec.pinned:
                    flags.append("pinned")
                if rec.tag:
                    flags.append(f"tag={rec.tag}")
                badge = f" [{' '.join(flags)}]" if flags else ""
                ts = rec.completed_at or rec.created_at or "-"
                size = format_bytes(directory_size(rec.directory))
                print(f"  - {rec.delivery_id} | {rec.status} | {ts} | {size}{badge}")
        if args.machine:
            return {"deliveries": rows}
        if args.as_json:
            print(json.dumps(rows, ensure_ascii=False, indent=2))
        return 0

    found = []
    missing = []
    for exp_id in experiments:
        record, source = resolve_latest_delivery(results_root, exp_id)
        if record is None:
            missing.append(exp_id)
        else:
            found.append((record, source))

    if not found:
        names = "、".join(missing)
        if args.machine:
            raise ServiceError(ErrorCode.INPUT_NOT_FOUND, f"实验 {names} 没有已完成的交付")
        print(f"错误：实验 {names} 在 {results_root} 下没有已完成的交付", file=sys.stderr)
        return 1

    if (args.path or args.open_target) and len(found) > 1:
        names = "、".join(rec.experiment for rec, _ in found)
        if args.machine:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"请从 {names} 中指定单个实验")
        print(f"错误：存在多个实验的交付（{names}），使用 --path 或 --open 时请用 --experiment 指定单个实验",
              file=sys.stderr)
        return 1

    pointer_file = latest_pointer_path(results_root)
    if args.machine and (args.path or args.open_target):
        record, _ = found[0]
        target = args.target or ("summary" if args.open_target and record.summary_path else "dir")
        resolved = record.target_path(target)
        if args.open_target:
            open_path(resolved)
        return {"path": str(resolved), "opened": args.open_target}
    if args.as_json or args.machine:
        items = []
        for rec, source in found:
            entry_exists = read_latest_pointer(results_root, rec.experiment) is not None
            items.append({
                **build_latest_payload(rec),
                "source": source,
                "pointer_path": str(pointer_file),
                "pointer_exists": entry_exists,
                "directory": str(rec.directory),
                "summary_path": str(rec.summary_path) if rec.summary_path else None,
                "index_path": str(rec.index_path),
            })
        if args.machine:
            return {"deliveries": items, "missing_experiments": missing}
        payload = items[0] if len(items) == 1 else items
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0

    if args.path:
        rec, _ = found[0]
        target = args.target or "dir"
        try:
            print(rec.target_path(target))
        except (FileNotFoundError, ValueError) as exc:
            print(f"错误：{exc}", file=sys.stderr)
            return 1
        return 0

    if args.open_target:
        rec, _ = found[0]
        target = args.target or ("summary" if rec.summary_path and rec.summary_path.is_file() else "dir")
        try:
            resolved = rec.target_path(target)
            open_path(resolved)
        except Exception as exc:
            print(f"错误：无法打开目标：{exc}", file=sys.stderr)
            return 1
        print(f"已打开：{resolved}")
        return 0

    for rec, source in found:
        source_label = "latest.json 指针" if source == "pointer" else "目录扫描回退"
        has_pointer = read_latest_pointer(results_root, rec.experiment) is not None
        print(f"[{rec.experiment}] 最新交付：{rec.delivery_id}（来源：{source_label}）")
        print(f"  目录：{rec.directory}")
        if rec.summary_path:
            print(f"  摘要：{rec.summary_path}")
        print(f"  索引：{rec.index_path}")
        print(f"  指针：{pointer_file}{'（存在）' if has_pointer else '（未生成）'}")
        meta = [f"时间：{rec.completed_at or rec.created_at or '-'}"]
        if rec.tag:
            meta.append(f"标签：{rec.tag}")
        if rec.pinned:
            meta.append("锁定：是")
        print("  " + "  ".join(meta))
        if rec.bundles:
            pairs = "，".join(f"{stage}={run_id}" for stage, run_id in rec.bundles.items())
            print(f"  关联结果包：{pairs}")
    return 0


def _run_clean(args: argparse.Namespace) -> int:
    from omicsweave.core.config import ConfigError
    from omicsweave.core.retention import RetentionPolicy, apply_cleanup, format_cleanup_report

    try:
        project, results_root = _resolve_results_root(args.results_dir)
    except ConfigError as exc:
        if args.machine:
            raise
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    base = _ENTRY_CLEAN_DEFAULTS.get() or RetentionPolicy.from_settings(project.config.pipeline.retention)
    keep_last = None if args.keep_all else (args.keep_last if args.keep_last is not None else base.keep_last)
    drop_failed = True if args.drop_failed else (False if args.keep_failed else base.drop_failed)
    prune_bundles = True if args.prune_bundles else (False if args.no_prune_bundles else base.prune_orphan_bundles)
    stale_hours = args.stale_running_hours if args.stale_running_hours is not None else base.stale_running_hours

    try:
        policy = RetentionPolicy(
            enabled=True,
            keep_last=keep_last,
            drop_failed=drop_failed,
            prune_orphan_bundles=prune_bundles,
            stale_running_hours=stale_hours,
        )
    except ValueError as exc:
        if args.machine:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, str(exc)) from exc
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    dry_run = not args.apply
    experiments = list(dict.fromkeys(args.experiment)) if args.experiment else None
    report = apply_cleanup(results_root, policy, experiments, dry_run=dry_run)
    if args.machine:
        from pydantic import TypeAdapter

        data = TypeAdapter(dict).dump_python(asdict(report), mode="json")
        if report.failed:
            raise ServiceError(ErrorCode.EXECUTION_FAILED, "部分目录清理失败", details={"report": data})
        return data
    print(format_cleanup_report(report))
    if dry_run and report.to_remove:
        print("（预览模式，未删除任何文件；确认执行清理请加 --apply）")
    return 1 if report.failed else 0


_HANDLERS = {"init": _run_init, "validate": _run_validate, "migrate": _run_migrate, "run": _run_run,
             "integrate": _run_integrate, "export": _run_export, "diff": _run_diff, "stats": _run_stats,
             "refs": _run_refs, "replicate": _run_replicate, "pipeline": _run_pipeline, "map": _run_map,
             "abstract": _run_abstract, "latest": _run_latest, "clean": _run_clean}

_EXECUTION = {"run", "integrate", "pipeline"}
_PROJECT_OPERATIONS = {"schema", "inspect-input", "config", "qc"}
_OPERATIONS = {"status", "logs", "cancel", "inspect", "query", "trace", "diagnose"} | _PROJECT_OPERATIONS
_MACHINE = _EXECUTION | _OPERATIONS | {"validate", "diff", "latest", "clean", "map", "abstract"}


def _wants_trial(args) -> bool:
    return any(getattr(args, name, None) is not None for name in ("override", "override_file", "intent", "rationale", "parent_run"))


def _execution(args):
    from omicsweave.service.execution import ExecutionRequest, execute, submit_execution

    values = vars(args)
    experiment = args.experiment
    trial = None
    if _wants_trial(args):
        from omicsweave.core.overrides import TrialOptions, parse_overrides

        trial = TrialOptions(overrides=parse_overrides(args.override, args.override_file),
                             intent=args.intent or "sensitivity_check", rationale=args.rationale or "",
                             parent_run_id=args.parent_run)
    request = ExecutionRequest(
        kind=args.module,
        experiments=[experiment] if isinstance(experiment, str) else experiment or [],
        omics=_split(values.get("omics")) or None, steps=_split(values.get("steps")) or None,
        modes=_split(values.get("mode")) or None, source=values.get("source"),
        start=values.get("start") or ("integrate" if values.get("source") else "run"),
        results_dir=values.get("results_dir"), integration_run=values.get("integration_run"),
        tag=values.get("tag"), pin=values.get("pin", False), dry_run=values.get("dry_run", False),
        run_gsea=values.get("run_gsea"), command=["omicsweave", *args.argv],
        retention=asdict(_ENTRY_RETENTION.get()) if _ENTRY_RETENTION.get() is not None else None,
        trial=trial,
    )
    if args.background:
        return submit_execution(project_dir(), request)
    return execute(project_dir(), request)


def _operation(args):
    from pydantic import ValidationError

    if args.module in _PROJECT_OPERATIONS:
        return _project_operation(args)

    project, root = _resolve_results_root(args.results_dir)
    if args.module in {"status", "logs", "cancel"}:
        from omicsweave.service.execution import job_operation

        values = vars(args)
        return job_operation(project.root, root, args.module, job=values.get("job"),
                             experiment=values.get("experiment"), active=values.get("active", False),
                             limit=values.get("limit", 20), tail=values.get("tail", 50),
                             since_event=values.get("since_event", 0), events_only=values.get("events_only", False))

    import pyarrow as pa

    from omicsweave.service.query import ResultView, TableQuery, query_table, run_digest, trace_entity

    try:
        view = ResultView(root, args.experiment, args.run, args.delivery, trial_id=args.trial)
        if args.module == "diagnose":
            from omicsweave.service.query import diagnose_view

            return diagnose_view(view, limit=args.limit, project_alpha=project.config.defaults.alpha)
        if args.module == "inspect":
            return run_digest(view)
        if args.module == "trace":
            return trace_entity(view, args.entity, tissue=args.tissue, contrast=args.contrast_id, limit=args.limit)
        try:
            request = TableQuery(
                table=args.table, omics=args.omics, mode=args.mode, columns=_split(args.columns) or None,
                filters=json.loads(args.filters) if args.filters is not None else [],
                limit=args.limit, offset=args.offset, sort_by=args.sort_by, ascending=not args.descending,
            )
        except (ValidationError, json.JSONDecodeError) as exc:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, str(exc)) from exc
        return query_table(view, request)
    except pa.ArrowException as exc:
        raise ServiceError(ErrorCode.CONTRACT_ERROR, f"无法读取结果表 {exc}") from exc


def _project_operation(args):
    from pydantic import TypeAdapter, ValidationError

    from omicsweave.core.param_schemas import schema
    from omicsweave.core.path_scope import PathScope
    from omicsweave.service.project import EditOptions, QCUpdate, inspect_inputs, patch_config, update_sample_qc

    if args.module == "schema":
        try:
            return schema(args.target)
        except KeyError as exc:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"未知 Schema {args.target!r}，可用 schema targets 查看") from exc
    scope = PathScope.create(project_dir())
    if args.module == "inspect-input":
        return inspect_inputs(scope, args.experiment, args.paths)
    try:
        options = EditOptions(dry_run=not args.apply, reason=args.reason, expected_sha256=args.expected_sha256)
        if args.module == "config":
            raw = args.patch if args.patch is not None else args.patch_file.read_text(encoding="utf-8-sig")
            patch = TypeAdapter(dict[str, object]).validate_json(raw)
        else:
            updates = TypeAdapter(list[QCUpdate]).validate_json(args.updates)
    except ValidationError as exc:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, str(exc)) from exc
    if args.module == "config":
        return patch_config(scope, args.experiment, patch, options, args.target)
    return update_sample_qc(scope, args.experiment, updates, options)


def _guess_json(raw: list[str]) -> tuple[str | None, bool]:
    """解析失败之前识别外层协议，保留旧模块及 latest 的专属参数。"""
    options = iter(enumerate(raw))
    command, command_index = None, len(raw)
    for index, token in options:
        if token in {"--project-dir", "--format"}:
            next(options, None)
        elif not token.startswith("-"):
            command, command_index = token, index
            break
    prefix, suffix = raw[:command_index], raw[command_index + 1:]

    def requested(tokens, alias):
        return any(token == "--format=json" or (token == "--format" and tokens[i + 1:i + 2] == ["json"])
                   or (alias and token == "--json") for i, token in enumerate(tokens))

    return command, requested(prefix, True) or (
        command not in COMMANDS and requested(suffix, command != "latest")
    )


@contextmanager
def _diagnostic_output():
    """保留 stdout 给信封，隔离 Python 与原生扩展写入的诊断信息。"""
    sys.stdout.flush()
    saved = None
    try:
        if sys.stdout.fileno() == 1 and sys.stderr.fileno() == 2:
            saved = os.dup(1)
            os.dup2(2, 1)
    except (AttributeError, OSError, ValueError):
        if saved is not None:
            os.close(saved)
            saved = None
    try:
        with redirect_stdout(sys.stderr):
            yield
    finally:
        if saved is not None:
            try:
                sys.stderr.flush()
            finally:
                os.dup2(saved, 1)
                os.close(saved)


def _response(command, root, data=None, error=None, warnings=()):
    from omicsweave.service.schemas import ServiceResponse, json_safe

    response = ServiceResponse(
        ok=error is None, command=command, project_dir=str(root) if root is not None else None,
        data=json_safe(data), error=error.detail if error else None, warnings=list(warnings),
    )
    print(json.dumps(response.model_dump(mode="json"), ensure_ascii=True, allow_nan=False, separators=(",", ":")))


def _dispatch_machine(args):
    from omicsweave.service.schemas import classify_error

    root = project_dir()
    try:
        with _diagnostic_output():
            if args.module in _EXECUTION:
                data = _execution(args)
            elif args.module in _OPERATIONS:
                data = _operation(args)
            else:
                data = _HANDLERS[args.module](args)
        _response(args.module, root, data=data, warnings=args.warnings)
        return 1 if args.module == "diff" and not data["identical"] else 0
    except Exception as exc:
        error = classify_error(exc)
        _response(args.module, root, error=error, warnings=args.warnings)
        return error.exit_code


def project_pipeline_main(
    project: Path,
    *,
    retention: "RetentionPolicy | None" = None,
    allow_clean: bool = False,
    argv: list[str] | None = None,
) -> None:
    """本地课题脚本入口，保持 IDE 与重定向日志中的 Unicode 可读。

    默认运行 ``pipeline``，也支持直接传 ``latest`` / ``clean`` 子命令；保留策略与 ``clean`` 开关均默认关闭，
    需在各课题入口文件（``run_analysis.py``）中显式开启。
    """
    if os.name == "nt" and "PYTHONIOENCODING" not in os.environ:
        for stream in (sys.stdout, sys.stderr):
            if hasattr(stream, "reconfigure"):
                stream.reconfigure(encoding="utf-8")
        os.environ["PYTHONIOENCODING"] = "utf-8"

    raw = list(sys.argv[1:] if argv is None else argv)
    if raw and raw[0] == "latest":
        main(["--project-dir", str(project), *raw])
        return

    if raw and raw[0] == "clean":
        if any(flag in ("-h", "--help") for flag in raw[1:]):
            main(["--project-dir", str(project), *raw])
            return
        if not allow_clean:
            print(
                "错误：本入口的清理命令默认关闭，未执行任何操作。"
                "如需启用，请在入口文件 run_analysis.py 中将 CLEAN_COMMAND_ENABLED 设为 True。",
                file=sys.stderr,
            )
            sys.exit(1)
        token = _ENTRY_CLEAN_DEFAULTS.set(retention)
        try:
            main(["--project-dir", str(project), *raw])
        finally:
            _ENTRY_CLEAN_DEFAULTS.reset(token)
        return

    cmd = raw if (raw and raw[0] == "pipeline") else ["pipeline", *raw]
    token = _ENTRY_RETENTION.set(retention)
    try:
        main(["--project-dir", str(project), *cmd])
    finally:
        _ENTRY_RETENTION.reset(token)


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    raw = list(sys.argv[1:] if argv is None else argv)
    command, json_parse = _guess_json(raw)
    token = _JSON_PARSE.set(json_parse)
    try:
        args, forwarded = parser.parse_known_args(raw)
        args.argv = raw
        if args.module in {"map", "abstract"} and "json" in (args.format or []):
            args.format = [fmt for fmt in args.format if fmt != "json"]
            args.response_format = "json"
        args.machine = getattr(args, "response_format", None) == "json"
        if args.machine and args.module not in _MACHINE:
            raise ServiceError(ErrorCode.CLI_USAGE_ERROR, f"{args.module} 暂未支持结构化 JSON")
        args.warnings = []
        if args.module in _HANDLERS or args.module in _OPERATIONS or args.module == "mcp":
            if args.module == "refs" and args.action == "build":
                args.forwarded = forwarded
            elif forwarded:
                parser.error(f"无法识别的参数：{' '.join(forwarded)}")
    except ServiceError as exc:
        _response(command, None, error=exc)
        sys.exit(exc.exit_code)
    finally:
        _JSON_PARSE.reset(token)

    if args.module == "mcp":
        from omicsweave.mcp.server import serve

        with project_context(args.project_dir):
            sys.exit(serve(project_dir(), args.allowed_data_root, allow_writes=args.allow_writes,
                           allow_custom_scripts=args.allow_custom_scripts))

    if args.module in _HANDLERS or args.module in _OPERATIONS:
        with project_context(args.project_dir):
            if args.machine:
                code = _dispatch_machine(args)
            elif (args.module in _OPERATIONS or getattr(args, "background", False)
                  or (args.module in _EXECUTION and (_wants_trial(args)
                      or (args.module != "pipeline" and args.dry_run)))):
                from omicsweave.service.schemas import classify_error, json_safe

                try:
                    result = _operation(args) if args.module in _OPERATIONS else _execution(args)
                    if getattr(args, "background", False):
                        print(f"后台作业：{result['job_id']}\n状态文件：{result['state_path']}")
                    else:
                        print(json.dumps(json_safe(result), ensure_ascii=False, indent=2))
                    code = 0
                except Exception as exc:
                    error = classify_error(exc)
                    print(f"错误：{error}", file=sys.stderr)
                    code = error.exit_code
            else:
                code = _HANDLERS[args.module](args)
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
