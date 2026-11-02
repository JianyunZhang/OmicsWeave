#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
``omicsweave mfa template / check / convert`` 三个动作的实现。

- ``template``：生成标准 Excel 模板；
- ``check``：检查一个 MFA 结果文件，或检查某个实验声明的 MFA 输入（可为多组结果，如每种示踪剂一组）
  并与 samples.tsv 核对样本；
- ``convert``：把旧版自研软件导出的 Excel 转换为标准格式。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from omicsweave.contract.dataset import ContractError
from omicsweave.contract.mfa import MFAResults, read_mfa, validate_mfa, write_mfa, write_template
from omicsweave.core.config import ConfigError, load_experiment, load_project
from omicsweave.core.model import Omics

from .legacy import convert_legacy_mfa, detect_legacy_kind, read_input_sets, read_legacy_mid

UNSPECIFIED = "(未指定)"


def _report(errors: list[str], warnings: list[str], notes: list[str] = ()) -> int:
    for note in notes:
        print(f"[说明] {note}")
    for warning in warnings:
        print(f"[警告] {warning}")
    for error in errors:
        print(f"[错误] {error}")
    print(f"共 {len(errors)} 个错误，{len(warnings)} 个警告")
    return 1 if errors else 0


def _summary(results: MFAResults) -> str:
    flux = results.flux
    n_fits = flux.groupby("sample_id")["fit_id"].nunique().median()
    parts = [
        f"{flux['sample_id'].nunique()} 个样本",
        f"{flux['tissue'].nunique()} 个组织",
        f"{flux[['tissue', 'reaction_id']].drop_duplicates().shape[0]} 个反应",
        f"每个样本约 {int(n_fits)} 次拟合" if (flux["fit_id"] != "").any() else "每个样本一个估计值",
    ]
    if results.mid is not None:
        parts.append(f"{len(results.mid)} 行 MID")
    return "，".join(parts)


def template_main() -> None:
    parser = argparse.ArgumentParser(description="生成 MFA 求解结果的标准 Excel 模板")
    parser.add_argument("output", nargs="?", default="mfa_results_template.xlsx", help="输出路径（默认当前目录）")
    args = parser.parse_args()
    path = write_template(args.output)
    print(f"已生成模板：{path}")
    print("填写 info / flux / mid / fits 工作表（各列说明见“说明”工作表），然后运行 omicsweave mfa check <文件>")


def _check_path(path: Path) -> int:
    if path.suffix.lower() == ".xlsx" and path.is_file():
        kind = detect_legacy_kind(path)
        if kind == "mid":
            notes: list[str] = []
            mid, _ = read_legacy_mid(path, notes)
            print(f"旧版自研软件 MID 文件：{mid['sample_id'].nunique()} 个样本，{len(mid)} 行 MID。"
                  "请与通量文件一起用 omicsweave mfa convert --flux ... --mid ... 转换")
            return _report([], [], notes)
        if kind == "flux":
            conversion = convert_legacy_mfa(path, model_id=UNSPECIFIED, flux_units=UNSPECIFIED)
            print(f"旧版自研软件通量文件：{_summary(conversion.results)}。可以读取；转换为标准格式请用 omicsweave mfa convert")
            return _report(*validate_mfa(conversion.results), conversion.notes)
    results = read_mfa(path)
    print(f"标准格式（模型 {results.info.model_id}，单位 {results.info.flux_units}）：{_summary(results)}")
    return _report(*validate_mfa(results))


def _check_experiment(experiment_id: str) -> int:
    exp = load_experiment(load_project(), experiment_id)
    settings = exp.config.omics.mfa
    if not settings.enabled:
        return _report([f"实验 {exp.id} 没有启用 MFA"], [])
    spec = settings.input
    files = [exp.project.resolve(f) for f in spec.files]
    if not files:
        return _report(["experiment.yaml 中 mfa.input.files 为空"], [])
    results, notes = read_input_sets(files, spec.format, model_id=getattr(spec, "model_id", None),
                                     flux_units=getattr(spec, "flux_units", None))
    errors, warnings = validate_mfa(results)

    rows = exp.samples[exp.samples["omics"] == Omics.MFA.value]
    declared = set(rows["sample_id"])
    found = set(results.samples)
    if found - declared:
        errors.append(f"MFA 结果中有 {len(found - declared)} 个样本不在 samples.tsv（omics = mfa）中："
                      f"{', '.join(sorted(found - declared)[:5])}")
    if declared - found:
        warnings.append(f"samples.tsv 中有 {len(declared - found)} 个 MFA 样本在结果文件里找不到："
                        f"{', '.join(sorted(declared - found)[:5])}")
    if "protocol" in rows.columns:
        tracers = rows.loc[rows["sample_id"].isin(found), "protocol"].astype(str).str.strip()
        counts = tracers[tracers != ""].value_counts()
        if len(counts) > 1:
            notes.append("示踪剂（samples.tsv 的 protocol 列）：" + "，".join(f"{t} {n} 个样本"
                                                                     for t, n in sorted(counts.items()))
                         + "；omicsweave run 将按示踪剂分层合并")
    print(f"实验 {exp.id} 的 MFA 输入（{spec.format}）：{_summary(results)}")
    return _report(errors, warnings, notes)


def check_main() -> None:
    parser = argparse.ArgumentParser(description="检查 MFA 求解结果（标准格式或旧版自研软件格式）")
    parser.add_argument("path", nargs="?", help="标准格式 .xlsx、Parquet 目录，或旧版 flux_raw_data.xlsx / mid_raw_data.xlsx")
    parser.add_argument("--experiment", help="改为检查该实验 experiment.yaml 声明的 MFA 输入，并与 samples.tsv 核对样本")
    args = parser.parse_args()
    if bool(args.path) == bool(args.experiment):
        parser.error("请给出文件路径，或用 --experiment 指定实验（二选一）")
    try:
        code = _check_path(Path(args.path)) if args.path else _check_experiment(args.experiment)
    except (ContractError, ConfigError) as exc:
        code = _report([str(exc)], [])
    sys.exit(code)


def convert_main() -> None:
    parser = argparse.ArgumentParser(description="把旧版自研 MFA 软件导出的 Excel 转换为标准格式")
    parser.add_argument("--flux", required=True, help="旧版 flux_raw_data.xlsx")
    parser.add_argument("--mid", help="旧版 mid_raw_data.xlsx（可选）")
    parser.add_argument("--model-id", required=True, help="代谢网络模型标识，如 central_carbon_v1")
    parser.add_argument("--flux-units", required=True, help="通量单位，如 'nmol/min/g tissue'")
    parser.add_argument("--software", help="求解软件及版本")
    parser.add_argument("--tracer", help="示踪剂，如 'U-13C glucose'")
    parser.add_argument("-o", "--output",
                        help="输出路径：以 .xlsx 结尾写 Excel，否则写 Parquet 目录（默认：通量文件旁的 <文件名>_omicsweave/）")
    args = parser.parse_args()

    flux_path = Path(args.flux)
    output = Path(args.output) if args.output else flux_path.with_name(f"{flux_path.stem}_omicsweave")
    try:
        conversion = convert_legacy_mfa(flux_path, args.mid, model_id=args.model_id, flux_units=args.flux_units,
                                        software=args.software, tracer=args.tracer)
        errors, warnings = validate_mfa(conversion.results)
        if not errors:
            write_mfa(conversion.results, output)
    except ContractError as exc:
        sys.exit(_report([str(exc)], []))

    print(f"转换结果：{_summary(conversion.results)}")
    if not errors:
        print(f"已写出：{output}")
    if conversion.conditions:
        print("旧文件中的样本分组（请据此填写 subjects.tsv / samples.tsv，分组不写入标准格式）：")
        for sample_id, condition in conversion.conditions.items():
            print(f"  {sample_id}\t{condition}")
    sys.exit(_report(errors, warnings, conversion.notes))
