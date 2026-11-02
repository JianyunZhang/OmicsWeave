#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把旧版各组学配置转换为统一的项目模型（``omicsweave migrate``）。

旧版每个组学各有一套配置（``experiments/<id>/`` 下的 ``mfa_data/``、``transcriptomics_data/``、
``proteomics_data/``、``metabolomics_data/``，各自带 ``groups.yaml``），多层配置的合并规则写在
各组学的加载器里。这里直接调用这些加载器得到实际生效的配置，再映射为 ``experiment.yaml``、
``subjects.tsv`` 与 ``samples.tsv``，保证迁移结果与旧流程读到的内容一致。

只读取配置文件，以及数据文件中的样本信息（MFA 通量文件的 ``sheet_information`` 工作表、
代谢组数据表的表头、蛋白组报告的 ``Run`` 列）；不移动、不复制、不修改任何数据文件，旧配置原样保留。
无法从旧配置推断的内容（比较的 role、逆转配对、把 ``group`` 拆成饮食 / 处理等因素、
MFA 的模型 ID 与通量单位等）写进迁移报告 ``MIGRATION.md``，由用户补充。
"""
from __future__ import annotations

import copy
import os
import re
import sys
from dataclasses import asdict, dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from pydantic import ValidationError

from omicsweave import __version__

from .config import EXPERIMENT_FILE, PROJECT_FILE
from .model import ID_PATTERN, ExperimentConfig, Omics, Species
from .project import project_context
from .scaffold import PROJECT_TEMPLATE

GROUP_FACTOR = "group"
# MFA 的一个样本是一次多组织联合拟合；各组织见结果表的 tissue 列
MFA_TISSUE = "multi_tissue"
REPORT_FILE = "MIGRATION.md"
OUTPUT_FILES = (EXPERIMENT_FILE, "subjects.tsv", "samples.tsv", REPORT_FILE)

# samples.tsv 的列顺序：必填列在前，其余按出现顺序排在后面
_LEADING_COLUMNS = ["sample_id", "subject_id", "omics", "tissue", "protocol", "batch", "sample_type", "files",
                    "data_column"]


class MigrationError(Exception):
    """旧配置无法迁移。"""


@dataclass
class _Part:
    """一个组学（含全部批次）从旧配置中读出的内容。"""
    omics: Omics
    settings: dict[str, Any]
    conditions: dict[str, Any]
    # 每行至少有 sample_id、subject_id、tissue；生物学样本另有 group（迁移时移到 subjects.tsv）
    samples: list[dict[str, str]]
    comparisons: list[tuple[str, str]]
    sources: list[Path]
    covariates: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


@dataclass
class Migration:
    experiment_id: str
    experiment: dict[str, Any]
    subjects: pd.DataFrame
    samples: pd.DataFrame
    report: str
    # 项目目录还没有 project.yaml 时生成的内容；已有时为 None（不改动）
    project_yaml: str | None = None


# ── 通用工具 ────────────────────────────────────────────────────────────────

def _rel(path: str | Path | None, root: Path) -> str | None:
    """项目目录内的路径写成相对路径（正斜杠），目录外的保留绝对路径。"""
    if path is None or str(path) == "":
        return None
    p = Path(path)
    if not p.is_absolute():
        p = root / p
    p = Path(os.path.normpath(p))
    try:
        return p.relative_to(root).as_posix()
    except ValueError:
        return str(p)


def _plain(value: Any) -> Any:
    """转为 YAML 可写的基本类型：Path → str，元组 → 列表。"""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def _compact(mapping: dict[str, Any]) -> dict[str, Any]:
    """去掉值为 None 或空容器的键。"""
    return {k: v for k, v in mapping.items() if v is not None and v != {} and v != []}


def _contrast_id(numerator: str, denominator: str) -> str:
    cid = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{numerator}_vs_{denominator}")
    return cid if re.match(ID_PATTERN, cid) else f"c_{cid}"


def _batches(discovered: list[dict[str, Any]], omics: Omics, notes: list[str]) -> list[tuple[str | None, str | None, Path]]:
    """旧加载器发现的 (protocol, batch, 配置文件)；跳过停用或缺少配置的批次并记入 notes。"""
    found = []
    for protocol in discovered:
        for batch in protocol["batches"]:
            label = "/".join(str(x) for x in (protocol["protocol_id"], batch["batch_id"]) if x) or "（单批次）"
            if batch["config_path"] is None:
                notes.append(f"{omics.value} 批次 {label} 在旧配置中声明了，但找不到它的 config.yaml，已跳过")
            elif not protocol["enabled"]:
                notes.append(f"{omics.value} 协议 {protocol['protocol_id']} 在旧配置中停用（enabled: false），批次 {label} 未迁移")
            else:
                found.append((protocol["protocol_id"], batch["batch_id"], Path(batch["config_path"])))
    return found


# ── MFA ────────────────────────────────────────────────────────────────────

def _mfa_samples(cfg, notes: list[str]) -> dict[str, str]:
    """MFA 样本 → 条件。与旧版 ``MFAAnalyser.load_flux_data`` 相同，取自通量文件的 sheet_information。"""
    flux = cfg.flux_file
    if flux is not None and Path(flux).is_file():
        with pd.ExcelFile(flux) as xl:
            info = pd.read_excel(xl, sheet_name="sheet_information") if "sheet_information" in xl.sheet_names else None
        if info is not None and "sheet_name" in info.columns:
            samples = {}
            for _, row in info.iterrows():
                name = row.get("complete_result_name")
                sid = str(name if pd.notna(name) and str(name).strip() else row["sheet_name"]).strip()
                samples[sid] = str(row.get("condition", "unknown")).strip().lower()
            return samples
        notes.append(f"MFA 通量文件 {Path(flux).name} 中没有 sheet_information 工作表，样本取自 groups.yaml 的 sample_groups")
    else:
        notes.append("找不到 MFA 通量文件，样本取自 groups.yaml 的 sample_groups；"
                     "数据就位后可用 omicsweave mfa check --experiment 核对样本")
    return dict(cfg.sample_groups)


def _mfa_subjects(samples: dict[str, str], animal_groups: dict[str, str], notes: list[str]) -> dict[str, str]:
    """MFA 样本 → 动物编号。

    旧配置没有记录二者的对应关系。只有每个样本名都以某个动物编号结尾、一一对应且分组一致时
    才采用推断结果，否则把每个 MFA 样本当作一只动物。
    """
    if not animal_groups or set(samples) <= set(animal_groups):
        return {sid: sid for sid in samples}
    from omicsweave.proteomics.experiment_config import get_sample_to_mouse_id_map

    matched = get_sample_to_mouse_id_map(list(samples), animal_groups)
    animals = list(matched.values())
    consistent = all(
        animal is not None and animal_groups[animal].lower() == samples[sid]
        for sid, animal in matched.items()
    )
    if consistent and len(set(animals)) == len(animals):
        notes.append("MFA 样本对应的动物编号是按样本名末尾的编号推断的（分组一致），请在 samples.tsv 中核对 subject_id")
        return {sid: str(animal) for sid, animal in matched.items()}
    notes.append("MFA 样本名与 groups.yaml 的动物编号（mouse_groups，用于血糖 QC）对应不上，"
                 "暂时把每个 MFA 样本当作一只动物；如与其他组学是同一批动物，请把 samples.tsv 中的 subject_id 改为动物编号")
    return {sid: sid for sid in samples}


def _collect_mfa(experiment_id: str, root: Path) -> _Part | None:
    from omicsweave.mfa import mfa_experiment_config as legacy

    if not (root / "experiments" / experiment_id / "mfa_data").is_dir():
        return None
    notes: list[str] = []
    batches = _batches(legacy.discover_experiment_protocols(experiment_id), Omics.MFA, notes)
    if not batches:
        raise MigrationError("mfa_data/ 下没有可用的 MFA 配置")
    default_path = legacy.default_mfa_config_path(experiment_id).resolve()

    conditions: dict[str, Any] = {}
    rows: list[dict[str, str]] = []
    entries: list[dict[str, Any]] = []
    comparisons: list[tuple[str, str]] = []
    default_cfg = None
    for protocol_id, batch_id, path in batches:
        cfg = legacy.load_mfa_config(experiment_id=experiment_id, config_path=str(path))
        if default_cfg is None or path.resolve() == default_path:
            default_cfg = cfg
        conditions.update({k: v for k, v in cfg.conditions.items() if k not in conditions})
        samples = _mfa_samples(cfg, notes)
        subjects = _mfa_subjects(samples, cfg.mouse_groups, notes)
        for sid, condition in samples.items():
            rows.append({"sample_id": sid, "subject_id": subjects[sid], "tissue": MFA_TISSUE,
                         "protocol": protocol_id or "", "batch": batch_id or "", "group": condition})
        if cfg.comparisons:
            pairs = cfg.comparisons
        else:
            # 旧流程未配置比较时，对数据中出现的条件做全部两两组合（df_analysis.differential_flux）
            pairs = list(combinations(sorted(set(samples.values())), 2))
            if pairs:
                notes.append("MFA 旧配置没有写 comparisons，旧流程会对所有条件两两比较；已全部列为比较，请删去不需要的")
        comparisons += [p for p in pairs if p not in comparisons]
        # 血糖 QC 按动物编号分组；subjects.tsv 覆盖不到的动物保留在这里
        animal_groups = cfg.mouse_groups if set(cfg.mouse_groups) - set(subjects.values()) else {}
        entries.append(_compact({
            "protocol": protocol_id,
            "batch": batch_id,
            "config": _rel(cfg.config_path, root),
            "flux_file": _rel(cfg.flux_file, root),
            "mid_file": _rel(cfg.mid_file, root),
            "glycemia_file": _rel(cfg.glycemia_file, root),
            "tracer": _plain(cfg.tracer),
            "batch_metadata": _plain(cfg.batch_metadata),
            "infusion_qc": _plain(cfg.infusion_qc),
            "animal_groups": _plain(animal_groups),
            "analysis": _plain(cfg.analysis),
        }))

    default_analysis = _plain(default_cfg.analysis)
    for entry in entries:
        if entry.get("analysis") == default_analysis:
            entry.pop("analysis")
        else:
            notes.append(f"MFA 批次 {entry.get('protocol')}/{entry.get('batch')} 的分析参数与默认批次不同，已保留在 adapter.batches 中")
    if len(entries) > 1:
        notes.append(f"MFA 共 {len(entries)} 个批次；input.files 已写默认批次"
                     f"（{default_cfg.protocol_id}/{default_cfg.batch_id}），全部批次见 omics.mfa.adapter.batches")
        if len({e.get("protocol") for e in entries}) > 1:
            notes.append("MFA 有多个示踪剂 protocol：如需多示踪剂整合，在 input.files 中依次列出各示踪剂批次的通量文件与"
                         "MID 文件（要求同一代谢网络模型），samples.tsv 的 protocol 列已标明每个样本的示踪剂；"
                         "同一示踪剂内多批次的整合尚未接入统一运行")
    files = [f for f in (_rel(default_cfg.flux_file, root), _rel(default_cfg.mid_file, root)) if f]
    if not files:
        notes.append("MFA 旧配置中没有通量文件路径，请在 omics.mfa.input.files 中补充")
    notes.append("MFA 输入按旧版自研软件格式（legacy_mfa_xlsx）迁移，旧文件中没有模型 ID 与通量单位，"
                 "请在 omics.mfa.input 中补充 model_id 与 flux_units")
    settings = {
        "enabled": True,
        "input": {"level": "table", "format": "legacy_mfa_xlsx", "files": files},
        "adapter": {"batches": entries},
        "analysis": default_analysis,
    }
    return _Part(Omics.MFA, settings, conditions, rows, comparisons,
                 [b[2] for b in batches], notes=notes)


# ── 转录组 ──────────────────────────────────────────────────────────────────

_SAMPLE_COLUMNS = ["sample", "sample_id", "complete_result_name"]
_R1_COLUMNS = ["file1", "fastq_1", "r1", "read1"]
_R2_COLUMNS = ["file2", "fastq_2", "r2", "read2"]


def _first_column(table: pd.DataFrame, candidates: list[str]) -> str | None:
    return next((c for c in candidates if c in table.columns), None)


def _resolve_legacy(path: str, config_path: Path, root: Path) -> Path:
    """旧配置中的相对路径：先按项目目录，再按配置文件所在目录。"""
    p = Path(path)
    if p.is_absolute():
        return p
    return root / p if (root / p).exists() else config_path.parent / p


def _collect_transcriptomics(experiment_id: str, root: Path) -> _Part | None:
    from omicsweave.transcriptomics import experiment_config as legacy

    exp_dir = root / "experiments" / experiment_id
    if not ((exp_dir / "transcriptomics_data").is_dir() or (exp_dir / "transcriptomics").is_dir()):
        return None
    notes: list[str] = []
    batches = _batches(legacy.discover_experiment_protocols(experiment_id), Omics.TRANSCRIPTOMICS, notes)
    if not batches:
        raise MigrationError("transcriptomics_data/ 下没有可用的转录组配置")

    conditions: dict[str, Any] = {}
    rows: list[dict[str, str]] = []
    entries: list[dict[str, Any]] = []
    comparisons: list[tuple[str, str]] = []
    covariates: list[str] = []
    first = None
    for protocol_id, batch_id, path in batches:
        cfg = legacy.load_transcriptomics_config(experiment_id=experiment_id, config_path=str(path))
        first = first or cfg
        conditions.update({k: v for k, v in cfg.conditions.items() if k not in conditions})
        sample_list = _resolve_legacy(cfg.sample_list_tsv, cfg.config_path, root)
        if not sample_list.is_file():
            notes.append(f"找不到转录组样本列表 {_rel(sample_list, root)}，该批次没有生成样本行")
            table = pd.DataFrame()
        else:
            table = pd.read_csv(sample_list, sep="\t", dtype=str, keep_default_na=False)
        sample_col, r1_col, r2_col = (_first_column(table, c) for c in (_SAMPLE_COLUMNS, _R1_COLUMNS, _R2_COLUMNS))
        if not table.empty and not (sample_col and r1_col):
            notes.append(f"转录组样本列表 {_rel(sample_list, root)} 缺少样本名或 R1 文件列，该批次没有生成样本行")
            table = pd.DataFrame()
        outside = set()
        for _, row in table.iterrows():
            # 与旧版 rnaseq_pipeline.get_sample_metadata 相同：样本名为 "<动物编号>_<组织>"
            sid = row[sample_col].strip()
            mouse_id, _, tissue = sid.partition("_")
            tissue = tissue or "unknown"
            if cfg.tissues and tissue not in cfg.tissues:
                outside.add(tissue)
            group = cfg.mouse_groups.get(sid) or cfg.mouse_groups.get(mouse_id) or ""
            batch = cfg.sample_batches.get(sid) or cfg.sample_batches.get(mouse_id) or batch_id or ""
            reads = [row[c] for c in (r1_col, r2_col) if c and row[c]]
            files = [_rel(os.path.join(cfg.raw_data_root, f), root) for f in reads]
            rows.append({"sample_id": sid, "subject_id": mouse_id, "tissue": tissue, "protocol": protocol_id or "",
                         "batch": batch, "group": group, "files": ";".join(files)})
        if outside:
            notes.append(f"转录组样本列表中有不在 params.tissues 中的组织：{'、'.join(sorted(outside))}")
        if cfg.sample_batches and "batch" not in covariates:
            covariates.append("batch")
        pairs = [(num, den) for num, den, _ in legacy.resolve_comparisons(cfg)]
        comparisons += [p for p in pairs if p not in comparisons]
        entries.append(_compact({
            "protocol": protocol_id,
            "batch": batch_id,
            "config": _rel(cfg.config_path, root),
            "sample_list": _rel(sample_list, root),
            "raw_data_root": _rel(cfg.raw_data_root, root),
            "technical_params": _plain(cfg.technical_params),
        }))

    analysis = _compact({
        "gene_aggregation": first.gene_aggregation,
        "qc_thresholds": _plain(first.qc_thresholds),
        "on_qc_fail": first.on_qc_fail,
        "qc_sample_size": first.qc_sample_size,
    })
    adapter = _compact({
        "name": "kallisto",
        "executable": _rel(first.kallisto_exe, root),
        "index": _rel(first.kallisto_index, root),
        "transcriptome_fasta": _rel(first.transcriptome_fasta, root),
        "threads": first.threads,
        "batches": entries,
    })
    settings = {
        "enabled": True,
        "input": {"level": "raw", "format": "fastq_gz", "layout": "paired", "strandedness": first.strandedness},
        "adapter": adapter,
        "analysis": analysis,
    }
    if any(";" not in r["files"] for r in rows):
        notes.append("部分转录组样本只有 R1 文件，已按双端（layout: paired）迁移；如为单端请改为 single")
    return _Part(Omics.TRANSCRIPTOMICS, settings, conditions, rows, comparisons,
                 [b[2] for b in batches], covariates, notes)


# ── 蛋白组 ──────────────────────────────────────────────────────────────────

_PROTEOMICS_SKIP = {"experiment_id", "config_path", "output_dir", "data_dirs", "tissues", "groups_path",
                    "conditions", "mouse_groups", "sample_metadata", "tissue_labels", "plot_colors",
                    "contrast_groups", "raw_dirs", "diann_exe", "search_fasta"}


def _collect_proteomics(experiment_id: str, root: Path) -> _Part | None:
    from omicsweave.proteomics import experiment_config as legacy

    if not legacy.default_config_path(experiment_id).is_file():
        return None
    cfg = legacy.load_proteomics_config(experiment_id=experiment_id)
    notes: list[str] = []
    rows: list[dict[str, str]] = []
    for tissue, report in cfg.data_dirs.items():
        if not Path(report).is_file():
            notes.append(f"找不到蛋白组 {tissue} 的报告 {_rel(report, root)}，该组织没有生成样本行")
            continue
        # 与旧流程相同：样本名由 DIA-NN 报告的 Run 列提取，再按末尾编号对应动物
        runs = pd.read_csv(report, sep="\t", usecols=["Run"], dtype=str)["Run"].dropna().unique()
        names = sorted({legacy.extract_sample_name(r) for r in runs})
        meta = legacy.build_sample_metadata(names, cfg.mouse_groups, cfg.sample_metadata, tissue=tissue)
        unmatched = [str(s) for s, m in meta["mouse_id"].items() if pd.isna(m)]
        if unmatched:
            notes.append(f"蛋白组 {tissue} 有 {len(unmatched)} 个样本对应不到动物编号：{'、'.join(unmatched[:5])}")
        for sample, info in meta.iterrows():
            mouse_id = "" if pd.isna(info["mouse_id"]) else str(info["mouse_id"])
            row = {"sample_id": str(sample), "subject_id": mouse_id, "tissue": tissue,
                   "group": str(info["condition"]) if mouse_id else ""}
            extra = {k: "" if pd.isna(v) else str(v) for k, v in info.items()
                     if k not in ("mouse_id", "condition", "tissue")}
            rows.append({**row, **extra})

    pairs = [(num, den) for num, den, _ in legacy.generate_comparisons(cfg.conditions, cfg.contrast_groups)]
    analysis = _plain({k: v for k, v in asdict(cfg).items()
                       if k not in _PROTEOMICS_SKIP and not k.startswith(("search_", "diann_"))})
    tissue_files = {t: _rel(p, root) for t, p in cfg.data_dirs.items()}
    adapter: dict[str, Any] = {"name": "diann"}
    if cfg.search_enabled:
        adapter.update(_compact({
            "executable": _rel(cfg.diann_exe, root), "fasta": _rel(cfg.search_fasta, root),
            "expected_version": cfg.diann_expected_version, "threads": cfg.search_threads,
            "qvalue": cfg.search_qvalue, "reanalyse": cfg.search_reanalyse,
            "smart_profiling": cfg.search_smart_profiling,
            "raw_dirs": {t: _rel(p, root) for t, p in cfg.raw_dirs.items()},
        }))
        notes.append("旧配置启用了 DIA-NN 搜库；已按现有报告（表格层）迁移，搜库参数保留在 omics.proteomics.adapter 中")
    settings = {
        "enabled": True,
        "input": {"level": "table", "format": "diann_report", "files": list(tissue_files.values()),
                  "tissue_files": tissue_files},
        "adapter": adapter,
        "analysis": analysis,
    }
    covariates = [c for c in cfg.design_covariates if c]
    return _Part(Omics.PROTEOMICS, settings, dict(cfg.conditions), rows, pairs, [cfg.config_path],
                 covariates, notes)


# ── 代谢组 ──────────────────────────────────────────────────────────────────

_METABOLOMICS_SKIP = {"experiment_id", "config_path", "output_dir", "tissue_data_paths", "groups_path",
                      "conditions", "mouse_groups", "tissue_labels", "plot_colors"}


def _collect_metabolomics(experiment_id: str, root: Path) -> _Part | None:
    from omicsweave.metabolomics import experiment_config as legacy
    from omicsweave.metabolomics.data_preprocessing import classify_data_columns

    if not legacy.default_config_path(experiment_id).is_file():
        return None
    cfg = legacy.load_metabolomics_config(experiment_id=experiment_id)
    notes: list[str] = []
    rows: list[dict[str, str]] = []
    for tissue, data in cfg.tissue_data_paths.items():
        if not Path(data).is_file():
            notes.append(f"找不到代谢组 {tissue} 的数据表 {_rel(data, root)}，该组织没有生成样本行")
            continue
        # 只读表头；与旧流程相同，按列名区分生物学样本、QC 与空白
        columns = pd.read_excel(data, nrows=0).columns
        classified = classify_data_columns(columns, set(cfg.mouse_groups))
        for column in classified.biological:
            mouse_id = classified.mouse_ids[column]
            rows.append({"sample_id": str(column), "subject_id": mouse_id, "tissue": tissue,
                         "sample_type": "sample", "group": cfg.mouse_groups.get(mouse_id, "")})
        for kind, columns_of_kind in (("qc", classified.qc), ("blank", classified.blank)):
            rows += [{"sample_id": str(c), "subject_id": "", "tissue": tissue, "sample_type": kind}
                     for c in columns_of_kind]
        if classified.unmapped_biological:
            notes.append(f"代谢组 {tissue} 有样本列对应不到 mouse_groups 中的动物，未迁移："
                         f"{'、'.join(map(str, classified.unmapped_biological[:5]))}")

    groups_doc = yaml.safe_load(cfg.groups_path.read_text(encoding="utf-8")) or {}
    pairs = [(num, den) for num, den, _ in legacy.generate_comparisons(cfg.conditions, groups_doc.get("contrast_groups"))]
    tissue_files = {t: _rel(p, root) for t, p in cfg.tissue_data_paths.items()}
    settings = {
        "enabled": True,
        "input": {"level": "table", "format": "peak_table", "files": list(tissue_files.values()),
                  "tissue_files": tissue_files},
        "analysis": _plain({k: v for k, v in asdict(cfg).items() if k not in _METABOLOMICS_SKIP}),
    }
    return _Part(Omics.METABOLOMICS, settings, dict(cfg.conditions), rows, pairs,
                 [cfg.config_path], notes=notes)


# ── 汇总 ────────────────────────────────────────────────────────────────────

def _canonical(conditions: dict[str, Any]) -> dict[str, str]:
    """条件名（小写）→ 旧配置中的写法；不同组学大小写不同时以先出现的为准。"""
    canon: dict[str, str] = {}
    for key in conditions:
        canon.setdefault(str(key).lower(), str(key))
    return canon


def _merge_samples(parts: list[_Part], canon: dict[str, str], notes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    subject_groups: dict[str, str] = {}
    records: list[dict[str, str]] = []
    seen: set[str] = set()
    renamed: list[str] = []
    for part in parts:
        for row in part.samples:
            row = dict(row)
            group = row.pop("group", "")
            group = canon.get(group.lower(), group) if group else ""
            subject = row.get("subject_id", "")
            if subject:
                if subject_groups.get(subject, group) != group:
                    notes.append(f"动物 {subject} 在不同组学中的分组不一致（{subject_groups[subject]} / {group}），"
                                 f"subjects.tsv 中保留了前者")
                subject_groups.setdefault(subject, group)
            if row["sample_id"] in seen:
                # 数据表中的列名不变，写在 data_column 中供运行器对应
                renamed.append(row["sample_id"])
                row["data_column"] = row["sample_id"]
                row["sample_id"] = f"{part.omics.value}_{row['sample_id']}"
            seen.add(row["sample_id"])
            records.append({"omics": part.omics.value, **row})
    if renamed:
        notes.append(f"{len(renamed)} 个样本名在不同组学中重复，后出现的已加上组学前缀（数据表中的列名写在 data_column 列）："
                     f"{'、'.join(renamed[:5])}")

    samples = pd.DataFrame(records, columns=None if records else _LEADING_COLUMNS[:4]).fillna("")
    for column in ("protocol", "batch", "sample_type", "files", "data_column"):
        if column in samples.columns and (samples[column] == "").all():
            samples = samples.drop(columns=column)
    leading = [c for c in _LEADING_COLUMNS if c in samples.columns]
    samples = samples[leading + [c for c in samples.columns if c not in leading]]
    subjects = pd.DataFrame({"subject_id": list(subject_groups), GROUP_FACTOR: list(subject_groups.values())})
    return subjects, samples


def _contrasts(parts: list[_Part], canon: dict[str, str]) -> list[dict[str, str]]:
    contrasts: list[dict[str, str]] = []
    seen: set[str] = set()
    for part in parts:
        for num, den in part.comparisons:
            num, den = canon.get(str(num).lower(), str(num)), canon.get(str(den).lower(), str(den))
            cid = _contrast_id(num, den)
            if cid not in seen:
                seen.add(cid)
                contrasts.append({"id": cid, "factor": GROUP_FACTOR, "numerator": num,
                                  "denominator": den, "role": "other"})
    return contrasts


def _report(experiment_id: str, parts: list[_Part], sources: list[str], subjects: pd.DataFrame,
            samples: pd.DataFrame, notes: list[str], has_integration: bool) -> str:
    counts = samples["omics"].value_counts()
    lines = [
        f"# 迁移报告：{experiment_id}",
        "",
        f"由 omicsweave {__version__} 的 `omicsweave migrate` 从旧版配置生成。旧配置与数据文件均未改动。",
        "",
        "## 旧配置入口",
        "",
        "各组学的 groups.yaml 以及上层的 experiment / protocol / batch 配置由原加载器按需一并读取。",
        "",
        *[f"- `{s}`" for s in sources],
        "",
        "## 生成的内容",
        "",
        f"- `experiment.yaml`：启用 {'、'.join(p.omics.value for p in parts)}",
        f"- `subjects.tsv`：{len(subjects)} 只动物",
        f"- `samples.tsv`：{len(samples)} 个样本（"
        + "，".join(f"{o} {int(n)}" for o, n in counts.items()) + "）",
        "",
        "## 需要补充或核对",
        "",
        f"- [ ] 分组因素 `{GROUP_FACTOR}` 直接沿用旧配置的条件名。条件由饮食、处理、时间等组合而成时，"
        "建议拆成多个因素（在 subjects.tsv 中加列），并相应改写比较的 factor / subset",
        "- [ ] 比较的 `role` 全部为 `other`：疾病效应改为 `disease_effect`、治疗效应改为 `treatment_effect`，"
        "再在 `reversal_pairs` 中配对",
        "- [ ] 参照水平（`factors.group.reference`）未设置",
    ]
    if has_integration:
        lines.append("- [ ] 旧配置中有多组学整合设置，已启用 `integration`；请在 `integration.modes` 中选择整合模式")
    lines += [f"- [ ] {note}" for note in notes]
    lines += [
        "",
        "## 下一步",
        "",
        "```bash",
        f"omicsweave validate --experiment {experiment_id}",
        "```",
        "",
        "旧版各组学流程目前仍读取旧配置（阶段 2 改为读取本文件）。",
        "",
    ]
    return "\n".join(lines)


def migrate_experiment(experiment_id: str, root: str | Path) -> Migration:
    """读取 ``root/experiments/<experiment_id>/`` 下的旧版配置，返回迁移结果（不写文件）。"""
    if not re.match(ID_PATTERN, experiment_id):
        raise MigrationError(f"实验 ID {experiment_id!r} 只能包含字母、数字、下划线、点和连字符，且以字母或数字开头")
    root = Path(root).expanduser().resolve()
    exp_dir = root / "experiments" / experiment_id
    if not exp_dir.is_dir():
        raise MigrationError(f"找不到实验目录 {exp_dir}")

    # 旧加载器按项目目录解析路径
    with project_context(root):
        parts = [part for collect in (_collect_transcriptomics, _collect_proteomics, _collect_metabolomics, _collect_mfa)
                 if (part := collect(experiment_id, root)) is not None]
        from omicsweave.integration.config import default_config_path as integration_config
        has_integration = integration_config(experiment_id).is_file()
    if not parts:
        raise MigrationError(f"{exp_dir} 下没有找到旧版组学配置"
                             "（mfa_data/、transcriptomics_data/、proteomics_data/、metabolomics_data/）")

    # 多批次时同一条提示会出现多次
    notes = list(dict.fromkeys(note for part in parts for note in part.notes))
    conditions: dict[str, Any] = {}
    for part in parts:
        for key, spec in part.conditions.items():
            if str(key).lower() not in _canonical(conditions):
                conditions[key] = spec
    canon = _canonical(conditions)
    subjects, samples = _merge_samples(parts, canon, notes)

    levels = list(canon.values())
    extra = sorted(set(subjects[GROUP_FACTOR]) - set(levels) - {""})
    if extra:
        notes.append(f"以下分组没有在 groups.yaml 的 conditions 中定义：{'、'.join(extra)}，已加入 levels")
        levels += extra
    covariates = [c for part in parts for c in part.covariates if c in samples.columns]
    sources = sorted({_rel(s, root) for part in parts for s in part.sources})
    experiment = {
        "schema_version": 1,
        "experiment": {
            "id": experiment_id,
            "description": "",
            "metadata": {
                "migrated_from": sources,
                "migrated_with": f"omicsweave {__version__}",
                "conditions": _plain(conditions),
            },
        },
        "design": {
            "subjects": "subjects.tsv",
            "samples": "samples.tsv",
            "factors": {GROUP_FACTOR: {"level": "subject", "levels": levels}},
            "contrasts": _contrasts(parts, canon),
            "reversal_pairs": [],
            "stratify_by": ["tissue"],
            "covariates": list(dict.fromkeys(covariates)),
        },
        "omics": {part.omics.value: copy.deepcopy(part.settings) for part in parts},
        "integration": {"enabled": has_integration, "modes": []},
    }
    try:
        ExperimentConfig.model_validate(experiment)
    except ValidationError as exc:
        raise MigrationError(f"生成的 experiment.yaml 不符合规范：{exc}") from exc

    project_yaml = None
    if not (root / PROJECT_FILE).is_file():
        project_id = re.sub(r"[^A-Za-z0-9_.-]", "_", root.name).lstrip("_.-") or "project"
        project_yaml = PROJECT_TEMPLATE.format(project_id=project_id, species=Species.MOUSE.value)
        notes.append(f"项目目录没有 project.yaml，已按小鼠（{Species.MOUSE.value}）生成，请核对项目 ID 与物种")

    report = _report(experiment_id, parts, sources, subjects, samples, notes, has_integration)
    return Migration(experiment_id, experiment, subjects, samples, report, project_yaml)


def experiment_yaml_text(migration: Migration) -> str:
    header = (f"# 由 omicsweave migrate 从旧版配置生成，待补充项见 {REPORT_FILE}。\n"
              "# 文件路径（input.files、samples.tsv 的 files 列）均相对于项目目录。\n")
    return header + yaml.safe_dump(migration.experiment, sort_keys=False, allow_unicode=True, width=120)


def write_migration(migration: Migration, root: str | Path, force: bool = False) -> list[Path]:
    """写出迁移结果，返回写出的文件。已有新格式文件时除非 ``force`` 否则报错；不会覆盖 project.yaml。"""
    root = Path(root).expanduser().resolve()
    exp_dir = root / "experiments" / migration.experiment_id
    targets = [exp_dir / name for name in OUTPUT_FILES]
    existing = [p for p in targets if p.exists()]
    if existing and not force:
        names = "、".join(p.name for p in existing)
        raise FileExistsError(f"{exp_dir} 中已有 {names}；如需覆盖请加 --force")

    written = []
    if migration.project_yaml is not None:
        (root / PROJECT_FILE).write_text(migration.project_yaml, encoding="utf-8")
        written.append(root / PROJECT_FILE)
    contents = {
        EXPERIMENT_FILE: experiment_yaml_text(migration),
        "subjects.tsv": migration.subjects.to_csv(sep="\t", index=False, lineterminator="\n"),
        "samples.tsv": migration.samples.to_csv(sep="\t", index=False, lineterminator="\n"),
        REPORT_FILE: migration.report,
    }
    for name, text in contents.items():
        (exp_dir / name).write_text(text, encoding="utf-8")
        written.append(exp_dir / name)
    return written


def run(experiment_ids: list[str], root: Path, dry_run: bool = False, force: bool = False) -> int:
    """命令行入口的实现：逐个迁移并校验，返回退出码。"""
    from .validate import validate_project

    code = 0
    for experiment_id in experiment_ids:
        try:
            migration = migrate_experiment(experiment_id, root)
            if dry_run:
                print(f"── {experiment_id}/{EXPERIMENT_FILE}（预览，未写出）──")
                print(experiment_yaml_text(migration))
                print(migration.report)
                continue
            written = write_migration(migration, root, force=force)
        except (MigrationError, FileExistsError) as exc:
            print(f"错误：{exc}", file=sys.stderr)
            code = 1
            continue
        print(f"已迁移实验 {experiment_id}：")
        for path in written:
            print(f"  {path.relative_to(root)}")
        report = validate_project(root, [experiment_id])
        print(f"\n迁移结果的校验（需补充的内容见 {REPORT_FILE}）：")
        print(report.format())
    return code
