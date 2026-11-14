#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
组学 × 组织 × 比较的 QC 可靠性（``qc/contrast_reliability.parquet``）。

样本层面的 QC 结论来自两处：

- 样本表 ``samples.tsv`` 的 ``qc_flag`` 列（``flag`` 或 ``exclude``，可选 ``qc_reason`` 说明原因），用于登记人工复核
  发现的问题（文库组成异常、离群样本、跨组织混入等）。改样本表后重新运行整合即可生效，不必重跑组学；
- 各组学的样本质量表 ``omics/<组学>/qc/sample_quality.parquet``（列 ``sample_id``、``metric``、``value``、``robust_z``、
  ``flag``、``reason``），由样本身份与组成 QC 生成，``flag`` 非空的样本视为已标记。

对每个（组学, 组织, 比较），统计两组中被标记的样本数，按以下规则给出可靠性：

- ``ok``：两组都没有标记样本，可靠性 1；
- ``qc_confounded``：两组的标记比例相差不小于 ``imbalance``（默认 0.5），质量问题与比较混杂，方向证据不可评估，可靠性 0；
- ``qc_insufficient``：去掉标记样本后有一组少于 2 只，可靠性 0；
- ``downweighted``：标记样本在两组间大致均衡，可靠性 = 可用样本数 / 总样本数。

``omicsweave run`` 结束时写出该表，并在不可评估的比较的差异结果行的 ``flags`` 中加上 ``qc_confounded`` 或
``qc_insufficient``；``flux_anchored`` 与 ``reversal`` 整合按当前样本表重新计算，结果写入整合结果包。
用户在 ``evidence_profile.quality_reliability`` 中显式配置的层级因子作为覆盖项，只替换数值，不解除不可评估。
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from omicsweave.contract.tables import read_differential, write_differential

from .config import Experiment
from .model import Omics
from .run import contrast_groups
from .validate import QC_FLAG_COLUMN

QC_REASON_COLUMN = "qc_reason"
SAMPLE_QUALITY_FILE = "sample_quality.parquet"
QC_DIR = "qc"
RELIABILITY_FILE = "contrast_reliability.parquet"
DEFAULT_IMBALANCE = 0.5
BLOCKED = ("qc_confounded", "qc_insufficient")
COLUMNS = ("omics", "tissue", "contrast_id", "numerator", "denominator", "n_numerator", "n_denominator",
           "n_flagged_numerator", "n_flagged_denominator", "n_usable_numerator", "n_usable_denominator",
           "reliability", "status", "reason", "flagged_samples")


def sample_quality_paths(run_dir: str | Path) -> list[Path]:
    return sorted(Path(run_dir).glob(f"omics/*/{QC_DIR}/{SAMPLE_QUALITY_FILE}"))


def sample_flags(experiment: Experiment, run_dir: str | Path | None = None) -> pd.DataFrame:
    """被标记的样本：``sample_id``、``flag``、``reason``、``source``，每个样本一行（多个来源的原因合并）。"""
    rows = []
    samples = experiment.samples
    if QC_FLAG_COLUMN in samples.columns:
        reasons = samples[QC_REASON_COLUMN] if QC_REASON_COLUMN in samples.columns else pd.Series("", index=samples.index)
        for sample, flag, reason in zip(samples["sample_id"], samples[QC_FLAG_COLUMN], reasons):
            if str(flag).strip():
                rows.append((str(sample), str(flag).strip(), str(reason or "").strip(), "samples.tsv"))
    for path in sample_quality_paths(run_dir) if run_dir is not None else []:
        table = pd.read_parquet(path)
        flagged = table[table["flag"].fillna("").astype(str).str.strip() != ""]
        source = path.relative_to(Path(run_dir)).as_posix()
        for sample, part in flagged.groupby("sample_id", sort=True):
            flag = "exclude" if part["flag"].astype(str).eq("exclude").any() else str(part["flag"].iloc[0])
            reason = "；".join(dict.fromkeys(part["reason"].fillna(part["metric"]).astype(str)))
            rows.append((str(sample), flag, reason, source))
    if not rows:
        return pd.DataFrame(columns=["sample_id", "flag", "reason", "source"])
    table = pd.DataFrame(rows, columns=["sample_id", "flag", "reason", "source"])
    return table.groupby("sample_id", as_index=False, sort=True).agg(
        flag=("flag", lambda s: "exclude" if (s == "exclude").any() else s.iloc[0]),
        reason=("reason", lambda s: "；".join(r for r in dict.fromkeys(s) if r)),
        source=("source", lambda s: ";".join(dict.fromkeys(s))))


def contrast_reliability(experiment: Experiment, keys: pd.DataFrame, flags: pd.DataFrame,
                         imbalance: float = DEFAULT_IMBALANCE) -> pd.DataFrame:
    """``keys``：差异结果中出现的（``omics``, ``tissue``, ``contrast_id``）。不在当前配置中的比较跳过。"""
    contrasts = {c.id: c for c in experiment.config.design.contrasts}
    tissue_of = dict(zip(experiment.samples["sample_id"], experiment.samples["tissue"]))
    flagged = dict(zip(flags["sample_id"], flags["reason"]))
    groups_cache: dict[tuple[str, str], dict[str, str]] = {}
    rows = []
    for omics, tissue, contrast_id in keys[["omics", "tissue", "contrast_id"]].drop_duplicates().itertuples(index=False):
        contrast = contrasts.get(str(contrast_id))
        if contrast is None:
            continue
        cache_key = (str(omics), contrast.id)
        if cache_key not in groups_cache:
            groups_cache[cache_key] = contrast_groups(experiment, contrast, Omics(str(omics)))
        groups = groups_cache[cache_key]
        # 多组织共用一份样本（如 MFA 每只动物一个样本）时，样本表中的组织不是结果中的组织，按全部样本统计
        members = {s: g for s, g in groups.items() if tissue_of.get(s) == tissue} or groups
        num = [s for s, g in members.items() if g == contrast.numerator]
        den = [s for s, g in members.items() if g == contrast.denominator]
        bad_num = [s for s in num if s in flagged]
        bad_den = [s for s in den if s in flagged]
        usable_num, usable_den = len(num) - len(bad_num), len(den) - len(bad_den)
        if not bad_num and not bad_den:
            status, reliability, reason = "ok", 1.0, ""
        else:
            frac_num = len(bad_num) / len(num) if num else 0.0
            frac_den = len(bad_den) / len(den) if den else 0.0
            if abs(frac_num - frac_den) >= imbalance:
                status, reliability = "qc_confounded", 0.0
            elif (bad_num and usable_num < 2) or (bad_den and usable_den < 2):
                status, reliability = "qc_insufficient", 0.0
            else:
                status = "downweighted"
                reliability = (usable_num + usable_den) / max(len(num) + len(den), 1)
            causes = "；".join(dict.fromkeys(flagged[s] for s in bad_num + bad_den if flagged[s]))
            reason = (f"{contrast.numerator} 组 {len(bad_num)}/{len(num)}、{contrast.denominator} 组 "
                      f"{len(bad_den)}/{len(den)} 个样本有 QC 标记" + (f"（{causes}）" if causes else ""))
        rows.append((str(omics), str(tissue), contrast.id, contrast.numerator, contrast.denominator, len(num),
                     len(den), len(bad_num), len(bad_den), usable_num, usable_den, round(reliability, 4), status,
                     reason, ";".join(sorted(bad_num + bad_den))))
    return pd.DataFrame(rows, columns=list(COLUMNS))


def differential_keys(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """各组学差异结果中的（``omics``, ``tissue``, ``contrast_id``）。"""
    parts = [t[["tissue", "contrast_id"]].drop_duplicates().assign(omics=omics) for omics, t in tables.items()
             if not t.empty]
    if not parts:
        return pd.DataFrame(columns=["omics", "tissue", "contrast_id"])
    return pd.concat(parts, ignore_index=True).astype(str)


def build_reliability(experiment: Experiment, run_dir: str | Path, tables: dict[str, pd.DataFrame],
                      imbalance: float = DEFAULT_IMBALANCE) -> tuple[pd.DataFrame, list[Path]]:
    """由结果包的差异结果与当前样本表计算可靠性表；返回 (表, 读取的样本质量表)。"""
    flags = sample_flags(experiment, run_dir)
    return contrast_reliability(experiment, differential_keys(tables), flags, imbalance), sample_quality_paths(run_dir)


def summary_note(table: pd.DataFrame) -> str | None:
    """可靠性表中非 ok 的情况，一句话；全部为 ok 时返回 None。"""
    counts = table["status"].value_counts()
    parts = [f"{status} {int(counts[status])} 个" for status in ("downweighted", *BLOCKED) if status in counts]
    if not parts:
        return None
    return f"样本 QC 标记影响（组学, 组织, 比较）：{'，'.join(parts)}；明细见 {QC_DIR}/{RELIABILITY_FILE}"


def append_flags(table: pd.DataFrame, reliability: pd.DataFrame, omics: str) -> pd.DataFrame:
    """差异结果中不可评估的（组织, 比较）行，在 ``flags`` 中加上可靠性状态（qc_confounded / qc_insufficient）。"""
    blocked = reliability[(reliability["omics"] == omics) & reliability["status"].isin(BLOCKED)]
    if blocked.empty:
        return table
    status = table[["tissue", "contrast_id"]].astype(str).merge(
        blocked[["tissue", "contrast_id", "status"]], on=["tissue", "contrast_id"], how="left")["status"]
    hit = status.notna().to_numpy()
    flags = table["flags"].fillna("").astype(str) if "flags" in table.columns else pd.Series("", index=table.index)
    joined = [";".join(p for p in (f, s) if p) if h and s not in f.split(";") else f
              for f, s, h in zip(flags, status.fillna(""), hit)]
    return table.assign(flags=joined)


def write_run_reliability(experiment: Experiment, run_dir: str | Path) -> list[str]:
    """``omicsweave run`` 结束时：写出 ``qc/contrast_reliability.parquet``，并给不可评估的差异结果行加标记。"""
    run_dir = Path(run_dir)
    paths = {p.parent.parent.name: p for p in sorted(run_dir.glob("omics/*/tables/differential.parquet"))}
    if not paths:
        return []
    tables = {omics: read_differential(path) for omics, path in paths.items()}
    reliability, _ = build_reliability(experiment, run_dir, tables)
    out = run_dir / QC_DIR / RELIABILITY_FILE
    out.parent.mkdir(parents=True, exist_ok=True)
    reliability.to_parquet(out, index=False)
    for omics, table in tables.items():
        flagged = append_flags(table, reliability, omics)
        if flagged is not table:
            write_differential(flagged, paths[omics])
    note = summary_note(reliability)
    return [note] if note else []
