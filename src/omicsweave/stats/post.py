#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
结果包的后处理：``omicsweave run`` 在各组学运行完后调用，``omicsweave stats`` 对已有结果包重新运行。

依次：

1. 声明式派生特征（G8，``derived`` 配置）；
2. ``design.models`` 的公式模型（G3），结果作为新的“比较”追加到差异结果，拟合记录写入 ``qc/model_fits.parquet``；
3. 差异结果补 ``estimand_type``（原有行为 pairwise）与 ``n_total``；
4. ``design.adaptive_shrinkage`` 开启时，按（组学, 组织, 比较）做自适应收缩，增加 ``lfsr``；
5. 入选特征的逐只去掉动物稳健性与 Type S / M（E5），写入 ``omics/<组学>/tables/robustness.parquet``，
   由单只动物驱动的行在 ``flags`` 中加 ``single_subject_driven``；
6. 可检出效应表 ``qc/design_power.parquet``（G3）；
7. 样本 QC 标记的可靠性表 ``qc/contrast_reliability.parquet``（G2）。
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from omicsweave.contract.tables import effect_standard_errors, read_differential, write_differential
from omicsweave.core.config import Experiment
from omicsweave.core.reliability import QC_DIR, write_run_reliability
from omicsweave.core.statistics.shrinkage import adaptive_shrinkage

from .models import fit_models
from .power import design_power
from .robustness import flag_rows, robustness

DIFFERENTIAL = "differential.parquet"
DERIVED = "derived_differential.parquet"
ROBUSTNESS = "robustness.parquet"


def _tables(run_dir: Path, name: str) -> dict[str, pd.DataFrame]:
    return {p.parent.parent.name: read_differential(p) for p in sorted(run_dir.glob(f"omics/*/tables/{name}"))}


def _complete(table: pd.DataFrame, shrink: bool) -> pd.DataFrame:
    table = table.copy()
    estimand = table["estimand_type"] if "estimand_type" in table.columns else pd.Series(index=table.index, dtype=object)
    table["estimand_type"] = estimand.fillna("pairwise").replace("", "pairwise")
    total = pd.to_numeric(table["n_numerator"], errors="coerce") + pd.to_numeric(table["n_denominator"], errors="coerce")
    existing = pd.to_numeric(table["n_total"], errors="coerce") if "n_total" in table.columns else total
    table["n_total"] = existing.fillna(total)
    if shrink:
        se, _ = effect_standard_errors(table)
        lfsr = pd.Series(float("nan"), index=table.index)
        for _, rows in table.groupby(["tissue", "contrast_id"], sort=False):
            lfsr.loc[rows.index] = adaptive_shrinkage(pd.to_numeric(rows["log2fc"], errors="coerce"),
                                                      se.loc[rows.index])[1]
        table["lfsr"] = lfsr
    return table


COPIED = ("data", "derived", "tables", "qc")


def stats_experiment(experiment_id: str, source: str | None = None, *, root: str | Path | None = None,
                     command: list[str] | None = None, run_id: str | None = None) -> tuple[Path, list[str]]:
    """``omicsweave stats``：复制源结果包中各组学的数据集与结果表到新的结果包，再按当前配置运行后处理。

    预处理与两两比较不重算，所以修改 ``design.models``、``derived``、``adaptive_shrinkage`` 或样本 QC 标记后不必重跑组学。
    """
    import shutil

    from omicsweave.core.run import finish_run, load_validated, start_run
    from omicsweave.integration.integrate import source_run

    project, experiment = load_validated(experiment_id, root)
    src = source_run(project.path("results"), experiment.id, source)
    ctx = start_run(project, experiment, run_id)
    copied = []
    for name in COPIED:
        for directory in sorted(src.glob(f"omics/*/{name}")):
            shutil.copytree(directory, ctx.run_dir / directory.relative_to(src))
            copied.append(directory)
    for stale in (ctx.run_dir / QC_DIR).glob("*.parquet"):
        stale.unlink()
    notes = postprocess(experiment, ctx.run_dir)
    omics = sorted({p.parent.name for p in copied})
    finish_run(ctx, copied, omics=[o for o in experiment.config.omics.enabled() if o.value in omics],
               command=command, parameters={"source_run": src.name, "steps": "stats"})
    return ctx.run_dir, notes


def postprocess(experiment: Experiment, run_dir: str | Path) -> list[str]:
    run_dir = Path(run_dir)
    design = experiment.config.design
    alpha = float(experiment.project.config.defaults.alpha)
    notes: list[str] = []
    qc = run_dir / QC_DIR
    qc.mkdir(parents=True, exist_ok=True)

    if experiment.config.derived:
        from .derived import write_derived
        notes += write_derived(experiment, run_dir)

    tables, derived = _tables(run_dir, DIFFERENTIAL), _tables(run_dir, DERIVED)
    if not tables and not derived:
        return notes
    if design.models:
        measured, derived_rows, fits, model_notes = fit_models(experiment, run_dir, tables, alpha)
        notes += model_notes
        fits.to_parquet(qc / "model_fits.parquet", index=False)
        for target, rows in ((tables, measured), (derived, derived_rows)):
            for omics, parts in rows.items():
                target[omics] = pd.concat([t for t in [target.get(omics), *parts] if t is not None and not t.empty],
                                          ignore_index=True)
        failed = fits[fits["status"] != "fit"]
        if not failed.empty:
            notes.append(f"公式模型有 {len(failed)} 处未拟合（组学 × 组织），原因见 {QC_DIR}/model_fits.parquet")

    tables = {o: _complete(t, design.adaptive_shrinkage) for o, t in tables.items()}
    derived = {o: _complete(t, design.adaptive_shrinkage) for o, t in derived.items()}

    robust = robustness(experiment, run_dir, tables, derived, alpha=alpha)
    for omics, table in robust.items():
        table.to_parquet(run_dir / "omics" / omics / "tables" / ROBUSTNESS, index=False)
    flagged = sum(int(t["single_subject_driven"].sum()) for t in robust.values())
    if flagged:
        notes.append(f"{flagged} 个入选特征去掉一只动物后反号或减半（single_subject_driven），"
                     f"见 omics/<组学>/tables/{ROBUSTNESS}")

    design_power(experiment, run_dir, alpha).to_parquet(qc / "design_power.parquet", index=False)

    for name, group in ((DIFFERENTIAL, tables), (DERIVED, derived)):
        for omics, table in group.items():
            write_differential(flag_rows(table, robust.get(omics)), run_dir / "omics" / omics / "tables" / name)
    return notes + write_run_reliability(experiment, run_dir)
