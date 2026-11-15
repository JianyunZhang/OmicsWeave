#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
项目配置校验（``omicsweave validate``）。

在 :mod:`omicsweave.core.model` 的字段校验之外，结合样本表与文件系统检查：

- ID 唯一、样本都能对应到动物；
- 因素列存在且取值合法；
- 每个比较在每个分层、每个组学中的生物学重复数，以及秩检验可达到的最小 p 值；
- 原始文件与表格输入是否存在、格式是否匹配；
- 外部程序是否可找到；
- 整合模式的前提条件。

规则对应 v2 设计文档第 5.6 节。错误（error）意味着无法运行；警告（warning）不阻止运行，但结果可能不可靠。
"""
from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pandas as pd

from .config import ConfigError, Experiment, Project, load_experiment, load_project
from .model import RAW_EXTENSIONS, IntegrationMode, Omics
from .statistics.nonparametric import min_achievable_pvalue

SUBJECT_COLUMNS = ("subject_id",)
SAMPLE_COLUMNS = ("sample_id", "subject_id", "omics", "tissue")
SAMPLE_TYPES = {"sample", "qc", "blank", "pool", "standard"}
# 样本表的人工 QC 标记列，由 omicsweave.core.reliability 读取
QC_FLAG_COLUMN = "qc_flag"
QC_FLAG_VALUES = ("flag", "exclude")
MIN_REPLICATES = 2


@dataclass
class Issue:
    level: Literal["error", "warning"]
    where: str
    message: str


@dataclass
class ValidationReport:
    issues: list[Issue] = field(default_factory=list)

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "error"]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "warning"]

    @property
    def ok(self) -> bool:
        return not self.errors

    def format(self) -> str:
        lines = [f"[{'错误' if i.level == 'error' else '警告'}] {i.where}: {i.message}" for i in self.issues]
        lines.append(f"共 {len(self.errors)} 个错误，{len(self.warnings)} 个警告")
        return "\n".join(lines)


class _Collector:
    def __init__(self, where: str):
        self.where = where
        self.issues: list[Issue] = []

    def error(self, message: str) -> None:
        self.issues.append(Issue("error", self.where, message))

    def warning(self, message: str) -> None:
        self.issues.append(Issue("warning", self.where, message))


def _preview(values, limit: int = 5) -> str:
    values = sorted(set(map(str, values)))
    shown = "、".join(values[:limit])
    return shown + (f" 等 {len(values)} 个" if len(values) > limit else "")


def _check_tables(exp: Experiment, out: _Collector) -> bool:
    """检查两张样本表的结构；结构有误时返回 False，后续依赖样本表的检查跳过。"""
    ok = True
    for name, table, required in (("subjects.tsv", exp.subjects, SUBJECT_COLUMNS),
                                  ("samples.tsv", exp.samples, SAMPLE_COLUMNS)):
        missing = [c for c in required if c not in table.columns]
        if missing:
            out.error(f"{name} 缺少必填列：{'、'.join(missing)}")
            ok = False
            continue
        id_col = required[0]
        blank = table[id_col] == ""
        if blank.any():
            out.error(f"{name} 有 {int(blank.sum())} 行的 {id_col} 为空")
        duplicated = table.loc[table[id_col].duplicated() & ~blank, id_col]
        if not duplicated.empty:
            out.error(f"{name} 中 {id_col} 重复：{_preview(duplicated)}")
    if not ok:
        return False

    # QC、空白、混样等非生物学样本可以不对应动物
    biological = _is_biological(exp.samples)
    no_subject = biological & (exp.samples["subject_id"] == "")
    if no_subject.any():
        out.error(f"samples.tsv 有 {int(no_subject.sum())} 个生物学样本的 subject_id 为空")
    linked = exp.samples.loc[exp.samples["subject_id"] != "", "subject_id"]
    unknown_subjects = set(linked) - set(exp.subjects["subject_id"])
    if unknown_subjects:
        out.error(f"samples.tsv 中的 subject_id 在 subjects.tsv 中找不到：{_preview(unknown_subjects)}")
    unknown_omics = set(exp.samples["omics"]) - {o.value for o in Omics}
    if unknown_omics:
        out.error(f"samples.tsv 中的 omics 取值无法识别：{_preview(unknown_omics)}")
    if "sample_type" in exp.samples.columns:
        bad_types = set(exp.samples["sample_type"]) - SAMPLE_TYPES - {""}
        if bad_types:
            out.error(f"samples.tsv 中的 sample_type 取值无法识别：{_preview(bad_types)}（可选：{'、'.join(sorted(SAMPLE_TYPES))}）")
    if QC_FLAG_COLUMN in exp.samples.columns:
        bad_flags = set(exp.samples[QC_FLAG_COLUMN]) - set(QC_FLAG_VALUES) - {""}
        if bad_flags:
            out.error(f"samples.tsv 中的 {QC_FLAG_COLUMN} 取值无法识别：{_preview(bad_flags)}"
                      f"（可选：{'、'.join(QC_FLAG_VALUES)}，或留空）")
    blank_tissue = exp.samples["tissue"] == ""
    if blank_tissue.any():
        out.error(f"samples.tsv 有 {int(blank_tissue.sum())} 行的 tissue 为空")
    return True


def _is_biological(samples: pd.DataFrame) -> pd.Series:
    if "sample_type" not in samples.columns:
        return pd.Series(True, index=samples.index)
    return samples["sample_type"].isin(["", "sample"])


def biological_samples(exp: Experiment) -> pd.DataFrame:
    """生物学样本（去掉 QC、空白等），并合并动物层面的因素列。"""
    samples = exp.samples[_is_biological(exp.samples)]
    subject_cols = [c for c in exp.subjects.columns if c not in samples.columns or c == "subject_id"]
    return samples.merge(exp.subjects[subject_cols], on="subject_id", how="left")


def _check_factors(exp: Experiment, merged: pd.DataFrame, out: _Collector) -> None:
    design = exp.config.design
    for name, factor in design.factors.items():
        table_name = "subjects.tsv" if factor.level == "subject" else "samples.tsv"
        table = exp.subjects if factor.level == "subject" else exp.samples
        if name not in table.columns:
            out.error(f"因素 {name} 声明在 {factor.level} 层，但 {table_name} 中没有该列")
            continue
        # 动物层面的因素直接查 subjects.tsv；样本层面的只查生物学样本（QC、空白可以不填）
        values = table[name] if factor.level == "subject" else merged[name]
        present = values[values != ""]
        if factor.type == "numeric":
            parsed = pd.to_numeric(present, errors="coerce")
            if parsed.isna().any():
                out.error(f"数值因素 {name} 中有非数值取值：{_preview(present[parsed.isna()])}")
        elif factor.levels is not None:
            unknown = set(present) - set(factor.levels)
            if unknown:
                out.error(f"因素 {name} 中有未声明的水平：{_preview(unknown)}（levels：{'、'.join(factor.levels)}）")
        if factor.level == "subject" and (exp.subjects[name] == "").any():
            out.warning(f"因素 {name} 在 {int((exp.subjects[name] == '').sum())} 只动物上为空，这些动物不参与相关比较")

    for column in [*design.stratify_by, *design.covariates]:
        if column not in merged.columns:
            out.error(f"stratify_by / covariates 中的列 {column} 在样本表中不存在")


def matches_level(series: pd.Series, value: str, numeric: bool) -> pd.Series:
    """样本表中取值等于某个因素水平的行；数值因素按数值比较（"4" 与 "4.0" 相同）。"""
    if numeric:
        return pd.to_numeric(series, errors="coerce") == float(value)
    return series == value


def _check_contrasts(exp: Experiment, merged: pd.DataFrame, out: _Collector) -> None:
    design = exp.config.design
    alpha = exp.project.config.defaults.alpha
    strata_cols = [c for c in design.stratify_by if c in merged.columns]
    for contrast in design.contrasts:
        factor = design.factors.get(contrast.factor)
        if factor is None or contrast.factor not in merged.columns:
            continue
        numeric = factor.type == "numeric"
        selected = merged
        for key, value in contrast.subset.items():
            if key not in merged.columns:
                out.error(f"比较 {contrast.id} 的 subset 列 {key} 在样本表中不存在")
                selected = selected.iloc[0:0]
                continue
            selected = selected[matches_level(selected[key], value, design.factors[key].type == "numeric")]
        for omics in exp.config.omics.enabled():
            part = selected[selected["omics"] == omics.value]
            if part.empty:
                out.warning(f"比较 {contrast.id} 在 {omics.value} 中没有任何样本")
                continue
            groups = part.groupby(strata_cols, dropna=False) if strata_cols else [((), part)]
            too_few, underpowered = [], []
            for stratum, rows in groups:
                label = "/".join(map(str, stratum if isinstance(stratum, tuple) else (stratum,))) or "全部"
                n_num = rows.loc[matches_level(rows[contrast.factor], contrast.numerator, numeric), "subject_id"].nunique()
                n_den = rows.loc[matches_level(rows[contrast.factor], contrast.denominator, numeric), "subject_id"].nunique()
                if min(n_num, n_den) < MIN_REPLICATES:
                    too_few.append(f"{label}（{n_num} vs {n_den}）")
                elif min_achievable_pvalue(n_num, n_den) >= alpha:
                    p_floor = min_achievable_pvalue(n_num, n_den)
                    underpowered.append(f"{label}（{n_num} vs {n_den}，最小 p = {p_floor:.3g}）")
            if too_few:
                out.warning(f"比较 {contrast.id} 在 {omics.value} 中生物学重复不足 {MIN_REPLICATES}："
                            f"{'；'.join(too_few)}，这些分层无法做组间检验")
            if underpowered:
                out.warning(f"比较 {contrast.id} 在 {omics.value} 中，秩检验在以下分层达不到 p < {alpha:g}："
                            f"{'；'.join(underpowered)}")


def _split_files(value: str) -> list[str]:
    return [part.strip() for part in value.split(";") if part.strip()]


def _check_inputs(exp: Experiment, out: _Collector) -> None:
    project = exp.project
    enabled = exp.config.omics.enabled()
    present_omics = set(exp.samples["omics"])
    for omics in Omics:
        if omics not in enabled and omics.value in present_omics:
            out.warning(f"samples.tsv 中有 {omics.value} 样本，但该组学未在 experiment.yaml 中启用")

    for omics in enabled:
        settings = exp.config.omics.settings(omics)
        rows = exp.samples[exp.samples["omics"] == omics.value]
        if rows.empty:
            out.error(f"{omics.value} 已启用，但 samples.tsv 中没有该组学的样本")
        spec = settings.input
        if spec.level == "raw":
            if "files" not in exp.samples.columns:
                out.error(f"{omics.value} 使用原始层输入，samples.tsv 需要 files 列")
            else:
                _check_raw_files(project, omics, spec.format, getattr(spec, "layout", None), rows, out)
        else:
            if not spec.files:
                out.error(f"{omics.value} 使用表格层输入，需要在 input.files 中列出文件")
            for f in spec.files:
                if not project.resolve(f).exists():
                    out.error(f"{omics.value} 的输入文件不存在：{f}")
            if omics is Omics.MFA:
                _check_mfa_input(project, spec, out)
                _check_mfa_tracers(rows, out)

        executable = settings.adapter.get("executable")
        if executable and shutil.which(str(executable)) is None and not project.resolve(executable).is_file():
            out.warning(f"{omics.value} 的外部程序找不到：{executable}")


def _check_mfa_input(project: Project, spec, out: _Collector) -> None:
    """MFA 输入的结构检查；文件内容与样本核对见 ``omicsweave mfa check --experiment``。"""
    files = spec.files
    if spec.format in ("omicsweave_mfa_xlsx", "legacy_mfa_xlsx"):
        not_xlsx = [f for f in files if not f.lower().endswith(".xlsx")]
        if not_xlsx:
            out.error(f"mfa 输入格式 {spec.format} 需要 .xlsx 文件：{_preview(not_xlsx)}")
    if spec.format == "mfa_tables":
        not_dir = [f for f in files if project.resolve(f).exists() and not project.resolve(f).is_dir()]
        if not_dir:
            out.error(f"mfa 输入格式 mfa_tables 需要 Parquet 目录：{_preview(not_dir)}")
    if spec.format == "legacy_mfa_xlsx":
        missing = [key for key in ("model_id", "flux_units") if not getattr(spec, key, None)]
        if missing:
            out.error(f"legacy_mfa_xlsx 需要在 mfa.input 中给出 {'、'.join(missing)}（旧文件中没有这些信息）")


def _check_mfa_tracers(rows: pd.DataFrame, out: _Collector) -> None:
    """多示踪剂时 samples.tsv 的 protocol 列区分示踪剂，必须每个 MFA 样本都填写。"""
    if "protocol" not in rows.columns:
        return
    protocol = rows["protocol"].astype(str).str.strip()
    tracers = set(protocol) - {""}
    blank = rows.loc[protocol == "", "sample_id"]
    if len(tracers) > 1 and not blank.empty:
        out.error(f"MFA 样本的 protocol（示踪剂）有 {len(tracers)} 种，但以下样本没有填写：{_preview(list(blank))}")


def _check_raw_files(project: Project, omics: Omics, fmt: str, layout: str | None,
                     rows: pd.DataFrame, out: _Collector) -> None:
    extensions = RAW_EXTENSIONS.get(fmt, ())
    missing, wrong_ext, no_scan, wrong_count, no_files = [], [], [], [], []
    for _, row in rows.iterrows():
        files = _split_files(row["files"])
        if not files:
            no_files.append(row["sample_id"])
            continue
        if layout == "paired" and len(files) != 2:
            wrong_count.append(row["sample_id"])
        for f in files:
            path = project.resolve(f)
            if not path.exists():
                missing.append(f)
            if extensions and not f.lower().endswith(extensions):
                wrong_ext.append(f)
            if fmt == "sciex_wiff" and not path.with_name(path.name + ".scan").exists():
                no_scan.append(f)
    if no_files:
        out.error(f"{omics.value} 样本的 files 为空：{_preview(no_files)}")
    if missing:
        out.error(f"{omics.value} 的原始文件不存在：{_preview(missing)}")
    if wrong_ext:
        out.warning(f"{omics.value} 的文件扩展名与格式 {fmt} 不符：{_preview(wrong_ext)}")
    if no_scan:
        out.error(f"SCIEX .wiff 缺少同名的 .wiff.scan：{_preview(no_scan)}")
    if wrong_count:
        out.warning(f"{omics.value} 为双端测序，但以下样本不是 2 个文件：{_preview(wrong_count)}")


def _check_integration(exp: Experiment, out: _Collector) -> None:
    integration = exp.config.integration
    if not integration.enabled:
        return
    enabled = exp.config.omics.enabled()
    if not integration.modes:
        out.warning("integration 已启用，但没有指定 modes")
    if exp.config.design.reversal_pairs and IntegrationMode.REVERSAL not in integration.modes:
        out.warning("声明了 design.reversal_pairs，但 integration.modes 中没有 reversal，逆转判定不会运行"
                    "（也可用 omicsweave integrate --mode reversal 单独运行）")
    for mode in integration.modes:
        if mode is IntegrationMode.FLUX_ANCHORED:
            if Omics.MFA not in enabled:
                out.error("flux_anchored 需要 MFA 数据，当前实验未启用 MFA")
            elif len(enabled) < 2:
                out.error("flux_anchored 需要 MFA 之外至少再启用 1 个组学")
        elif mode is IntegrationMode.PATHWAY_CONVERGENCE:
            if len(enabled) < 2:
                out.error("pathway_convergence 需要至少启用 2 个组学")
        elif mode is IntegrationMode.REVERSAL:
            if not exp.config.design.reversal_pairs:
                out.error("reversal 需要在 design.reversal_pairs 中声明（疾病效应, 处理效应）配对")
        elif mode is IntegrationMode.SAMPLE_LEVEL:
            subjects_by_omics = {
                o: set(exp.samples.loc[exp.samples["omics"] == o.value, "subject_id"]) for o in enabled
            }
            paired = [
                (a, b) for i, a in enumerate(enabled) for b in enabled[i + 1:]
                if len(subjects_by_omics[a] & subjects_by_omics[b]) >= MIN_REPLICATES
            ]
            if not paired:
                out.error(f"sample_level 需要至少 2 个组学测自同一批动物（共同动物 ≥ {MIN_REPLICATES}）")


def validate_experiment(exp: Experiment) -> list[Issue]:
    out = _Collector(f"实验 {exp.id}")
    if exp.root.name != exp.id:
        out.warning(f"experiment.id（{exp.id}）与目录名（{exp.root.name}）不一致")
    if _check_tables(exp, out):
        merged = biological_samples(exp)
        _check_factors(exp, merged, out)
        _check_contrasts(exp, merged, out)
        _check_inputs(exp, out)
    _check_integration(exp, out)
    _check_registries(exp, out)
    return out.issues


def _check_registries(exp: Experiment, out: _Collector) -> None:
    """实验有结果包时，用其中的数据集对照知识表（omicsweave.core.refs），不一致处给出警告。"""
    from .refs import audit, audit_messages, latest_bundle

    try:
        bundle = latest_bundle(exp)
        if bundle is None:
            return
        table = audit(exp, bundle)
        # 手写代谢物 ID 表与 RefMet 的冲突属于包内知识表，不随项目报告（用 omicsweave refs audit 查看）
        for message in audit_messages(table[table["check"] != "metabolite_table"]):
            out.warning(f"{message}（结果包 {bundle.name}；明细用 omicsweave refs audit 导出）")
    except Exception as exc:  # noqa: BLE001 — 审计失败不影响校验结论
        out.warning(f"知识表审计未能完成：{exc}")


def validate_project(root: str | Path | None = None, experiments: list[str] | None = None) -> ValidationReport:
    """校验项目及其实验；``experiments`` 为空时校验全部实验。"""
    report = ValidationReport()
    try:
        project = load_project(root)
    except ConfigError as exc:
        report.issues += [Issue("error", "project.yaml", m) for m in exc.messages]
        return report

    if experiments:
        targets: list[str | Path] = list(experiments)
    else:
        targets = list(project.experiment_dirs())
        if not targets:
            report.issues.append(Issue("warning", "项目", f"{project.path('experiments')} 下没有找到任何实验"))
    for target in targets:
        try:
            exp = load_experiment(project, target)
        except ConfigError as exc:
            name = Path(target).name
            report.issues += [Issue("error", f"实验 {name}", m) for m in exc.messages]
            continue
        report.issues += validate_experiment(exp)
    return report


__all__ = ["Issue", "ValidationReport", "biological_samples", "matches_level", "validate_experiment",
           "validate_project"]
