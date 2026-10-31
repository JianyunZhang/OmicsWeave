#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 的统一运行器（``omicsweave run --omics mfa``）：读取 ``experiment.yaml`` 声明的 MFA 输入，
按其中的比较做差异通量分析，写出统一结果表。

计算沿用 :mod:`omicsweave.mfa.df_analysis`。输入先读为标准格式，再经 :mod:`omicsweave.adapters.mfa.legacy`
还原为旧版 ``MFAAnalyser`` 的内部数据，所以同一份数据与旧流程（``omicsweave mfa run``）的结果逐值一致。
样本分组来自 samples.tsv / subjects.tsv，比较来自 ``design.contrasts``（含 ``subset``），不再读取文件中的分组。

输出（``omics/mfa/``）::

    data/<组织>/                    每只动物的代表净通量（分析层数据集；行 = 反应，列 = 样本）
    qc/fitting_loss.parquet         每个样本的拟合残差分布
    qc/mid_fit_residuals.parquet    实测与模拟 MID 的残差（有 MID 时）
    tables/differential.parquet     统一差异结果表（所有比较 × 所有组织）
    details/                        旧版格式的明细 CSV；各比较的差异通量宽表在 details/contrasts/<比较>/ 下
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from omicsweave.adapters.mfa.legacy import (
    convert_legacy_mfa,
    legacy_all_fits,
    legacy_flux_matrix,
    legacy_mid_matrices,
)
from omicsweave.contract.dataset import ContractError, Dataset, DatasetMeta, ValueType, empty_features, write_dataset
from omicsweave.contract.mfa import MFAResults, read_mfa, validate_mfa
from omicsweave.contract.tables import write_differential
from omicsweave.core.config import Project
from omicsweave.core.model import Omics
from omicsweave.core.run import OmicsRun, RunContext, RunError, contrast_groups, sample_labels

from .df_analysis import differential_flux, mid_fit_residuals, quality_control, sample_flux_matrix
from .export import differential_flux_to_contract

STEPS = ("qc", "differential")
# df_analysis.differential_flux 从输出目录读取 sample_flux_matrix 写出的可辨识性表
_IDENTIFIABILITY = Path("01_df_results", "all_tissues", "flux_matrix", "flux_identifiability.csv")


def _preview(values, limit: int = 5) -> str:
    values = sorted(map(str, values))
    return "、".join(values[:limit]) + (f" 等 {len(values)} 个" if len(values) > limit else "")


def load_results(project: Project, spec) -> tuple[MFAResults, list[Path], list[str]]:
    """按 ``omics.mfa.input`` 读取 MFA 结果，返回 (结果, 读取的文件, 转换提示)。"""
    files = [project.resolve(f) for f in spec.files]
    if not files:
        raise RunError("omics.mfa.input.files 为空")
    try:
        if spec.format == "legacy_mfa_xlsx":
            conversion = convert_legacy_mfa(files[0], files[1] if len(files) > 1 else None,
                                            model_id=spec.model_id, flux_units=spec.flux_units,
                                            software=getattr(spec, "software", None),
                                            tracer=getattr(spec, "tracer", None))
            return conversion.results, files[:2], conversion.notes
        return read_mfa(files[0]), files[:1], []
    except ContractError as exc:
        raise RunError(str(exc)) from exc


def _write_flux_datasets(matrix: pd.DataFrame, results: MFAResults, out: Path) -> None:
    """每个组织一个分析层数据集：每只动物的代表净通量（拟合点经 IQR 过滤后的中位数）。"""
    for tissue, part in matrix.groupby("tissue", sort=True):
        wide = part.pivot(index="reaction", columns="sample_id", values="net_flux")
        meta = DatasetMeta(
            omics=Omics.MFA, tissue=str(tissue), value_type=ValueType.FLUX, scale="linear",
            feature_id_type="mfa_reaction", units=results.info.flux_units,
            source={"model_id": results.info.model_id, "software": results.info.software,
                    "estimate": "per-animal median of IQR-filtered net flux fits"},
        )
        write_dataset(out / "data" / str(tissue), Dataset(wide, empty_features(wide.index), meta))


def run(ctx: RunContext) -> OmicsRun:
    experiment = ctx.experiment
    settings = experiment.config.omics.mfa
    analysis = settings.analysis
    out = ctx.omics_dir(Omics.MFA)
    details = out / "details"
    record = OmicsRun(Omics.MFA)

    results, record.inputs, notes = load_results(ctx.project, settings.input)
    record.notes += [f"MFA 输入：{n}" for n in notes]
    errors, warnings = validate_mfa(results)
    if errors:
        raise RunError("MFA 输入不符合标准格式：" + "；".join(errors))
    record.notes += [f"MFA 输入：{w}" for w in warnings]

    declared = set(experiment.samples.loc[experiment.samples["omics"] == Omics.MFA.value, "sample_id"])
    found = set(results.samples)
    if found - declared:
        record.notes.append(f"MFA 结果中有 {len(found - declared)} 个样本不在 samples.tsv 中，未参与分析："
                            f"{_preview(found - declared)}")
    if declared - found:
        record.notes.append(f"samples.tsv 中有 {len(declared - found)} 个 MFA 样本在结果中找不到："
                            f"{_preview(declared - found)}")
    used = sorted(found & declared)
    if not used:
        raise RunError("MFA 结果与 samples.tsv（omics = mfa）没有共同的样本")

    flux_df = legacy_flux_matrix(results).loc[used]
    all_fits = {s: fits for s, fits in legacy_all_fits(results).items() if s in used}
    tissues = sorted(results.flux["tissue"].unique())
    # 只用于明细表与 QC 表中的标注；差异分析按各比较的分组进行
    labels = sample_labels(experiment, used)

    if ctx.wants("qc"):
        (out / "qc").mkdir(parents=True, exist_ok=True)
        loss = quality_control(all_fits, labels, details)
        if loss is not None:
            loss.to_parquet(out / "qc" / "fitting_loss.parquet", index=False)
        if results.mid is not None:
            simulated, measured = legacy_mid_matrices(results)
            keep = [s for s in used if s in simulated.columns]
            residuals = mid_fit_residuals(simulated[keep], measured[keep], labels, details)
            if residuals is not None:
                residuals.to_parquet(out / "qc" / "mid_fit_residuals.parquet", index=False)

    contrasts_run = []
    if ctx.wants("differential"):
        matrix = sample_flux_matrix(flux_df, all_fits, labels, tissues, details, analysis)
        if matrix is not None:
            _write_flux_datasets(matrix, results, out)
        tables = []
        for contrast in experiment.config.design.contrasts:
            groups = {s: g for s, g in contrast_groups(experiment, contrast, Omics.MFA).items() if s in flux_df.index}
            if not groups:
                record.notes.append(f"比较 {contrast.id} 在 MFA 中没有样本，已跳过")
                continue
            # 每个比较单独一个明细目录，避免水平相同、subset 不同的比较互相覆盖
            contrast_dir = details / "contrasts" / contrast.id
            if (details / _IDENTIFIABILITY).is_file():
                (contrast_dir / _IDENTIFIABILITY).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(details / _IDENTIFIABILITY, contrast_dir / _IDENTIFIABILITY)
            pair = (contrast.numerator, contrast.denominator)
            found_tables = differential_flux(flux_df, all_fits, groups, list(pair), tissues, contrast_dir,
                                             comparisons=[pair], analysis_config=analysis)
            if not found_tables:
                record.notes.append(f"比较 {contrast.id} 在 MFA 中有一组没有样本，已跳过")
                continue
            tables.append(differential_flux_to_contract(
                next(iter(found_tables.values())), numerator=contrast.numerator, denominator=contrast.denominator,
                contrast_id=contrast.id, contrast_role=contrast.role.value, flux_units=results.info.flux_units))
            contrasts_run.append(contrast.id)
        if tables:
            write_differential(pd.concat(tables, ignore_index=True), out / "tables" / "differential.parquet")

    record.parameters = {
        "input_format": settings.input.format,
        "model_id": results.info.model_id,
        "flux_units": results.info.flux_units,
        "samples": len(used),
        "contrasts": contrasts_run,
        "analysis": analysis,
    }
    return record
