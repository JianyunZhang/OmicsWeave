#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 的统一运行器（``omicsweave run --omics mfa``）：读取 ``experiment.yaml`` 声明的 MFA 输入，
按其中的比较做差异通量分析，写出统一结果表。

计算沿用 :mod:`omicsweave.mfa.df_analysis`。输入先读为标准格式，再经 :mod:`omicsweave.adapters.mfa.legacy`
还原为旧版 ``MFAAnalyser`` 的内部数据，所以同一份数据与旧流程（``omicsweave mfa run``）的结果逐值一致。
样本分组来自 samples.tsv / subjects.tsv，比较来自 ``design.contrasts``（含 ``subset``），不再读取文件中的分组。

``input.files`` 可以列出多组结果（如每种示踪剂一组，要求同一模型、同一单位）。samples.tsv 的 ``protocol``
列给出每个 MFA 样本的示踪剂；出现两种及以上示踪剂时，先在每个示踪剂内单独分析，再按
:mod:`omicsweave.mfa.multi_tracer` 分层合并，``tables/differential.parquet`` 为合并结果。

输出（``omics/mfa/``）::

    data/<组织>/                    每只动物的代表净通量（分析层数据集；行 = 反应，列 = 样本）
    tables/sample_flux_fits.parquet 每只动物每条反应的拟合分布：代表值、拟合次数、全部拟合点的四分位数
                                    （fit_q25 / fit_q75）、IQR 过滤后的四分位数与标准差、是否贴在盒约束边界上（at_bound）
    tables/flux_groups.parquet      化学计量上锁死的通量组（:mod:`omicsweave.mfa.flux_groups`）；差异分析的 BH 校正中
                                    每组只计一次，结果表多 flux_group、flux_group_size 两列。``analysis.flux_groups``：
                                    ``enabled``（默认 true）、``min_abs_correlation``（默认 1 − 1e−6）
    qc/fitting_loss.parquet         每个样本的拟合残差分布
    qc/mid_fit_residuals.parquet    实测与模拟 MID 的残差（有 MID 时）
    tables/differential.parquet     统一差异结果表（所有比较 × 所有组织）
    figures/                        火山图 volcano_<组织>_<比较>.pdf / .png（``figures`` 步骤，见
                                    :mod:`omicsweave.visualization.volcano`）：纵轴为主检验的 p 值，只给效应量的
                                    比较为 |Hedges' g|；横轴半宽不超过 10（1024 倍），更大的倍数来自接近 0 的
                                    净通量分母，钉在边界上
    details/                        旧版格式的明细 CSV；各比较的差异通量宽表在 details/contrasts/<比较>/ 下

多示踪剂时另有::

    tables/differential_by_tracer.parquet   各示踪剂单独分析的差异结果（多一列 tracer）
    tables/tracer_coverage.parquet          各示踪剂对每条反应的可辨识性与是否参与合并
    tables/tracer_flux_levels.parquet       按条件合并的绝对通量及示踪剂间一致性检验
    qc/tracer_concordance.parquet           各示踪剂相对其余示踪剂的系统偏离
    details/tracers/<示踪剂>/                各示踪剂的明细 CSV（结构同 details/）
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from omicsweave.adapters.mfa.legacy import (
    legacy_all_fits,
    legacy_flux_matrix,
    legacy_mid_matrices,
    read_input_sets,
)
from omicsweave.contract.dataset import ContractError, Dataset, DatasetMeta, ValueType, empty_features, write_dataset
from omicsweave.contract.mfa import MFAResults, apply_mfa_tissue_aliases, validate_mfa
from omicsweave.contract.tables import validate_differential, write_differential
from omicsweave.core.config import Experiment, Project
from omicsweave.core.model import Omics
from omicsweave.core.run import OmicsRun, RunContext, RunError, contrast_groups, sample_labels
from omicsweave.core.tissues import format_alias_note
from omicsweave.visualization.volcano import Thresholds, write_bundle_volcanoes

from . import flux_groups as fg
from . import multi_tracer
from .df_analysis import (
    differential_flux,
    mid_fit_residuals,
    quality_control,
    sample_flux_matrix,
    selection_thresholds,
)
from .export import differential_flux_to_contract

STEPS = ("qc", "differential", "figures")
# 火山图横轴半宽的上限（log2FC = 10，即 1024 倍）：净通量的分母接近 0 时 |log2FC| 可达 20 以上，属于数值假象
_VOLCANO_X_CAP = 10.0
# df_analysis.differential_flux 从输出目录读取 sample_flux_matrix 写出的可辨识性表
_IDENTIFIABILITY = Path("01_df_results", "all_tissues", "flux_matrix", "flux_identifiability.csv")


def _preview(values, limit: int = 5) -> str:
    values = sorted(map(str, values))
    return "、".join(values[:limit]) + (f" 等 {len(values)} 个" if len(values) > limit else "")


def load_results(project: Project, spec) -> tuple[MFAResults, list[Path], list[str]]:
    """按 ``omics.mfa.input`` 读取 MFA 结果，返回 (结果, 读取的文件, 提示)；多组结果合并为一份。"""
    files = [project.resolve(f) for f in spec.files]
    if not files:
        raise RunError("omics.mfa.input.files 为空")
    try:
        results, notes = read_input_sets(files, spec.format, model_id=getattr(spec, "model_id", None),
                                         flux_units=getattr(spec, "flux_units", None),
                                         software=getattr(spec, "software", None),
                                         tracer=getattr(spec, "tracer", None))
    except ContractError as exc:
        raise RunError(str(exc)) from exc
    return results, files, notes


def tracer_strata(experiment: Experiment, samples: list[str], results: MFAResults) -> dict[str, str] | None:
    """MFA 样本 → 示踪剂（samples.tsv 的 ``protocol`` 列）；只有一种示踪剂时返回 None，按单示踪剂分析。

    以 samples.tsv 为准：``protocol`` 全部填写且只有一种取值时，即使各结果文件 ``info.tracer`` 的写法不同，
    也按单示踪剂分析。``protocol`` 没有（完整）填写时，用结果文件的 ``info.tracer`` 防止不同示踪剂被混成一组。
    """
    table = experiment.samples
    rows = table[(table["omics"] == Omics.MFA.value) & table["sample_id"].isin(samples)]
    protocol = rows["protocol"] if "protocol" in rows.columns else pd.Series("", index=rows.index)
    tracer_of = dict(zip(rows["sample_id"], protocol.fillna("").astype(str).str.strip()))
    blank = [s for s, t in tracer_of.items() if not t]
    if len({t for t in tracer_of.values() if t}) <= 1:
        declared = [t.strip() for t in (results.info.tracer or "").split(";") if t.strip()]
        if blank and len(declared) > 1:
            raise RunError(f"MFA 结果文件记录了多种示踪剂（{'、'.join(declared)}），但 samples.tsv 的 protocol 列没有为每个"
                           "样本填写；请填写每个 MFA 样本的示踪剂，否则不同示踪剂的动物会被当作同一组比较。"
                           "若实为同一种示踪剂（只是写法不同），给所有 MFA 样本填同一个 protocol 即按单示踪剂分析")
        return None
    if blank:
        raise RunError(f"samples.tsv 中有 {len(blank)} 个 MFA 样本没有填写 protocol（示踪剂）：{_preview(blank)}")
    return tracer_of


SAMPLE_FIT_COLUMNS = ["sample_id", "tissue", "reaction", "tracer", "net_flux", "n_fits", "fit_q25", "fit_q75",
                      "fit_p25", "fit_p75", "fit_std", "at_bound"]
FLUX_GROUP_SETTINGS = {"enabled": True, "min_abs_correlation": fg.DEFAULT_MIN_ABS_CORRELATION}


def _write_sample_fits(matrix: pd.DataFrame, out: Path) -> None:
    """每只动物每条反应的拟合分布（数据集只存代表值）：误差线、贴边标记等逐只动物的读数从这里取。"""
    table = matrix.copy()
    if "tracer" not in table.columns:
        table["tracer"] = ""
    for column in SAMPLE_FIT_COLUMNS:
        if column not in table.columns:
            table[column] = pd.NA
    (out / "tables").mkdir(parents=True, exist_ok=True)
    table[SAMPLE_FIT_COLUMNS].to_parquet(out / "tables" / "sample_flux_fits.parquet", index=False)


def _flux_group_settings(analysis: dict) -> dict:
    own = dict((analysis or {}).get("flux_groups") or {})
    unknown = sorted(set(own) - set(FLUX_GROUP_SETTINGS))
    if unknown:
        raise RunError(f"omics.mfa.analysis.flux_groups 中有无法识别的参数：{'、'.join(unknown)}"
                       f"（可选：{'、'.join(FLUX_GROUP_SETTINGS)}）")
    return {**FLUX_GROUP_SETTINGS, **own}


def _flux_groups(fit_sets: dict[str, dict], settings: dict, out: Path, record: OmicsRun) -> dict[str, dict] | None:
    """各模型（示踪剂）与合并分析用的通量组：返回 范围 → {(组织, 反应): 组}；关闭时返回 None。"""
    if not settings["enabled"]:
        return None
    threshold = float(settings["min_abs_correlation"])
    tables = {name: fg.locked_groups({name: fits}, min_abs_correlation=threshold, scope=name)
              for name, fits in fit_sets.items()}
    if len(fit_sets) > 1:
        tables["all"] = fg.locked_groups(fit_sets, min_abs_correlation=threshold, scope="all")
    found = [t for t in tables.values() if not t.empty]
    combined = pd.concat(found, ignore_index=True) if found else pd.DataFrame(columns=fg.COLUMNS)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    combined.to_parquet(out / "tables" / "flux_groups.parquet", index=False)
    main = tables.get("all", next(iter(tables.values())) if tables else pd.DataFrame())
    if not main.empty:
        record.notes.append(f"通量组：{main['flux_group'].nunique()} 组、{len(main)} 条反应在所有拟合中成线性关系"
                            "（化学计量上锁死），BH 校正中每组只计一次")
    return {name: fg.group_map(table) for name, table in tables.items()}


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
    _write_sample_fits(matrix, out)


def _differential_tables(experiment: Experiment, flux_df: pd.DataFrame, all_fits: dict, tissues: list[str],
                         details: Path, analysis: dict, flux_units: str, record: OmicsRun,
                         scope: str = "", groups: dict | None = None) -> dict[str, pd.DataFrame]:
    """按 ``design.contrasts`` 对 ``flux_df`` 中的样本做差异通量分析，返回 比较 → 统一差异结果表。

    ``details`` 下须已有 ``sample_flux_matrix`` 写出的可辨识性表；``scope`` 用于提示中标明示踪剂。
    """
    tables: dict[str, pd.DataFrame] = {}
    where = f"（示踪剂 {scope}）" if scope else ""
    for contrast in experiment.config.design.contrasts:
        members = {s: g for s, g in contrast_groups(experiment, contrast, Omics.MFA).items() if s in flux_df.index}
        if not members:
            record.notes.append(f"比较 {contrast.id} 在 MFA{where}中没有样本，已跳过")
            continue
        # 每个比较单独一个明细目录，避免水平相同、subset 不同的比较互相覆盖
        contrast_dir = details / "contrasts" / contrast.id
        if (details / _IDENTIFIABILITY).is_file():
            (contrast_dir / _IDENTIFIABILITY).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(details / _IDENTIFIABILITY, contrast_dir / _IDENTIFIABILITY)
        pair = (contrast.numerator, contrast.denominator)
        found = differential_flux(flux_df, all_fits, members, list(pair), tissues, contrast_dir,
                                  comparisons=[pair], analysis_config=analysis, flux_groups=groups)
        if not found:
            record.notes.append(f"比较 {contrast.id} 在 MFA{where}中有一组没有样本，已跳过")
            continue
        tables[contrast.id] = differential_flux_to_contract(
            next(iter(found.values())), numerator=contrast.numerator, denominator=contrast.denominator,
            contrast_id=contrast.id, contrast_role=contrast.role.value, flux_units=flux_units)
    return tables


def _run_single(experiment: Experiment, results: MFAResults, flux_df: pd.DataFrame, all_fits: dict,
                tissues: list[str], labels: dict[str, str], analysis: dict, out: Path, record: OmicsRun) -> list[str]:
    details = out / "details"
    matrix = sample_flux_matrix(flux_df, all_fits, labels, tissues, details, analysis)
    if matrix is not None:
        _write_flux_datasets(matrix, results, out)
    groups = _flux_groups({"all": all_fits}, _flux_group_settings(analysis), out, record)
    tables = _differential_tables(experiment, flux_df, all_fits, tissues, details, analysis,
                                  results.info.flux_units, record, groups=(groups or {}).get("all"))
    if tables:
        write_differential(pd.concat(tables.values(), ignore_index=True), out / "tables" / "differential.parquet")
    return list(tables)


def _run_multi_tracer(experiment: Experiment, results: MFAResults, flux_df: pd.DataFrame, all_fits: dict,
                      tissues: list[str], labels: dict[str, str], tracer_of: dict[str, str], analysis: dict,
                      out: Path, record: OmicsRun) -> list[str]:
    """先在每个示踪剂内单独分析（可辨识性、差异通量），再分层合并。"""
    try:
        settings = multi_tracer.PoolingSettings.from_analysis(analysis)
    except (TypeError, ValueError) as exc:
        raise RunError(str(exc)) from exc
    units = results.info.flux_units
    animals, identifiability, by_tracer = [], [], []
    tracers = sorted(set(tracer_of.values()))
    fit_sets = {t: {s: all_fits[s] for s in flux_df.index if tracer_of.get(s) == t and s in all_fits} for t in tracers}
    group_maps = _flux_groups({t: f for t, f in fit_sets.items() if f}, _flux_group_settings(analysis), out, record) or {}
    for tracer in tracers:
        samples = [s for s in flux_df.index if tracer_of.get(s) == tracer]
        details = out / "details" / "tracers" / tracer
        fits = fit_sets[tracer]
        matrix = sample_flux_matrix(flux_df.loc[samples], fits, labels, tissues, details, analysis)
        if matrix is None:
            record.notes.append(f"示踪剂 {tracer} 没有可用的通量，未参与合并")
            continue
        animals.append(matrix.assign(tracer=tracer))
        identifiability.append(pd.read_csv(details / _IDENTIFIABILITY).assign(tracer=tracer))
        tables = _differential_tables(experiment, flux_df.loc[samples], fits, tissues, details, analysis,
                                      units, record, scope=tracer, groups=group_maps.get(tracer))
        by_tracer += [t.assign(tracer=tracer) for t in tables.values()]
    if not animals:
        return []

    animals_df = pd.concat(animals, ignore_index=True)
    _write_flux_datasets(animals_df, results, out)
    coverage = multi_tracer.tracer_coverage(pd.concat(identifiability, ignore_index=True), settings)
    tables_dir = out / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_parquet(tables_dir / "tracer_coverage.parquet", index=False)
    if by_tracer:
        per_tracer = pd.concat(by_tracer, ignore_index=True)
        for tracer, part in per_tracer.groupby("tracer"):
            problems = validate_differential(part.drop(columns="tracer"))
            if problems:
                raise RunError(f"示踪剂 {tracer} 的差异结果不符合规范：" + "；".join(problems))
        per_tracer.to_parquet(tables_dir / "differential_by_tracer.parquet", index=False)

    pooled_tables, contrasts_run = [], []
    for contrast in experiment.config.design.contrasts:
        groups = contrast_groups(experiment, contrast, Omics.MFA)
        pooled = multi_tracer.pool_contrast(animals_df, groups, contrast.numerator, contrast.denominator,
                                            coverage, settings, seed_key=contrast.id)
        table = multi_tracer.to_differential(pooled, contrast_id=contrast.id, contrast_role=contrast.role.value,
                                             flux_units=units, settings=settings,
                                             groups=group_maps.get("all", group_maps.get(tracers[0])) if group_maps else None)
        if table.empty:
            continue
        pooled_tables.append(table)
        contrasts_run.append(contrast.id)
        counts = table["tracer_evidence"].value_counts().to_dict()
        record.notes.append(f"比较 {contrast.id} 多示踪剂合并：{int(table['selected'].sum())} 条通量入选；"
                            + "、".join(f"{k} {v}" for k, v in sorted(counts.items())))
    if pooled_tables:
        write_differential(pd.concat(pooled_tables, ignore_index=True), tables_dir / "differential.parquet")

    levels, deviations = multi_tracer.flux_levels(animals_df, labels, coverage, settings)
    if not levels.empty:
        levels.to_parquet(tables_dir / "tracer_flux_levels.parquet", index=False)
        concordance = multi_tracer.tracer_concordance(deviations)
        if not concordance.empty:
            (out / "qc").mkdir(parents=True, exist_ok=True)
            concordance.to_parquet(out / "qc" / "tracer_concordance.parquet", index=False)
    return contrasts_run


def _multi_tracer_notes(experiment: Experiment, tracer_of: dict[str, str]) -> list[str]:
    """合并的前提：各示踪剂是独立的动物队列；同一示踪剂内的多个批次不单独建模。"""
    rows = experiment.samples.drop_duplicates("sample_id").set_index("sample_id").loc[list(tracer_of)]
    tracer = pd.Series(tracer_of)
    notes = []
    subjects = rows["subject_id"][rows["subject_id"] != ""]
    per_subject = tracer[subjects.index].groupby(subjects).nunique()
    shared = sorted(per_subject[per_subject > 1].index)
    if shared:
        notes.append(f"{len(shared)} 只动物接受了不止一种示踪剂（{_preview(shared)}）；合并假定各示踪剂的动物相互独立，"
                     "p 值可能偏乐观")
    if "batch" in rows.columns:
        per_tracer = rows["batch"].groupby(tracer[rows.index]).nunique()
        mixed = sorted(per_tracer[per_tracer > 1].index)
        if mixed:
            notes.append(f"示踪剂 {'、'.join(mixed)} 内有多个批次，按同一组动物分析，未做批次整合")
    return notes


def run(ctx: RunContext) -> OmicsRun:
    experiment = ctx.experiment
    settings = experiment.config.omics.mfa
    analysis = settings.analysis
    out = ctx.omics_dir(Omics.MFA)
    details = out / "details"
    record = OmicsRun(Omics.MFA)

    results, record.inputs, notes = load_results(ctx.project, settings.input)
    record.notes += [f"MFA 输入：{n}" for n in notes]
    aliases = experiment.config.tissue_aliases_for(Omics.MFA)
    try:
        results, applied_aliases = apply_mfa_tissue_aliases(results, aliases)
    except ContractError as exc:
        raise RunError(str(exc)) from exc
    if applied_aliases:
        record.notes.append(format_alias_note(Omics.MFA.value, applied_aliases))
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
    tracer_of = tracer_strata(experiment, used, results)

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
        if tracer_of is None:
            contrasts_run = _run_single(experiment, results, flux_df, all_fits, tissues, labels, analysis, out, record)
        else:
            record.notes += _multi_tracer_notes(experiment, tracer_of)
            contrasts_run = _run_multi_tracer(experiment, results, flux_df, all_fits, tissues, labels, tracer_of,
                                              analysis, out, record)

    if ctx.wants("figures"):
        cut = selection_thresholds(analysis)
        record.notes += write_bundle_volcanoes(
            out, experiment.config.design.contrasts, metric="pvalue", x_cap=_VOLCANO_X_CAP,
            thresholds=Thresholds(p=cut["pvalue"], log2fc=cut["log2fc"], effect_size=cut["hedges_g"]))

    record.parameters = {
        "input_format": settings.input.format,
        "model_id": results.info.model_id,
        "flux_units": results.info.flux_units,
        "samples": len(used),
        "contrasts": contrasts_run,
        "analysis": analysis,
    }
    if tracer_of is not None:
        counts = pd.Series(tracer_of).value_counts()
        record.parameters["tracers"] = {str(t): int(counts[t]) for t in sorted(counts.index)}
    if applied_aliases:
        record.parameters["tissue_aliases_applied"] = applied_aliases
    return record
