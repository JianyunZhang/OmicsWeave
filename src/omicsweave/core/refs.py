#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
知识表注册框架（G7，``omicsweave refs``）：每张知识表由“生成器 + 来源说明 + 对照数据的审计”组成。

- ``omicsweave refs list``：已登记的知识表、文件、生成日期与行数；
- ``omicsweave refs build <名称> [生成器参数]``：运行生成器（会联网下载来源，结果写入包内 ``core/idmap/data/``）；
- ``omicsweave refs audit --experiment <实验> [--run <运行>]``：用结果包中的数据集对照知识表；
  ``omicsweave validate`` 在实验有结果包时运行同一审计，以警告报告不一致。

审计项：

- ``protein_accession``：蛋白组数据自带的 UniProt 登记号与蛋白注册表中该基因符号的登记号不一致；
- ``metabolite_ids``：代谢组数据自带的 HMDB / KEGG 与代谢物注册表（RefMet）中同名代谢物的 ID 不一致；
- ``metabolite_table``：手写的代谢物 ID 表与 RefMet 注册表不一致（不需要数据）；
- ``mapping_coverage``：反应映射中的直接酶在各组织蛋白组数据中的覆盖率；
- ``reaction_roles`` 与 ``regulator_effect``：反应成员知识表（G17）的角色一致性与调控因子的作用方向（不需要数据）。

``metabolite_table``、``reaction_roles``、``regulator_effect`` 只检查包内知识表（:data:`KNOWLEDGE_CHECKS`），
``validate`` 不随项目报告。
"""
from __future__ import annotations

import importlib
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from omicsweave.contract.dataset import read_dataset

from .config import Experiment

AUDIT_COLUMNS = ["check", "omics", "tissue", "feature_id", "feature_name", "data_value", "registry_value", "issue"]


# 只检查包内知识表、与项目数据无关的审计项
KNOWLEDGE_CHECKS = ("metabolite_table", "reaction_roles", "regulator_effect")


@dataclass(frozen=True)
class RegistryInfo:
    name: str
    description: str
    files: tuple[str, ...]
    builder: str            # "模块:函数"，接受 argv
    sources_key: str
    package: str = "omicsweave.core.idmap"   # 文件所在的包；idmap 的表在其 data/ 子目录
    subdir: str = "data"


REGISTRIES: dict[str, RegistryInfo] = {
    "genes": RegistryInfo("genes", "小鼠基因：Ensembl GRCm39 cDNA + UniProt，按 MGI 对齐",
                          ("mouse_gene_registry.tsv.gz", "mouse_gene_aliases.tsv.gz"),
                          "omicsweave.core.idmap.build_registry:main", "ensembl"),
    "proteins": RegistryInfo("proteins", "小鼠蛋白：UniProt Swiss-Prot 登记号与 EC 号（与 genes 一起生成）",
                             ("mouse_protein_registry.tsv.gz",), "omicsweave.core.idmap.build_registry:main", "uniprot"),
    "metabolites": RegistryInfo("metabolites", "代谢物：RefMet 标准名称、类别与 HMDB / KEGG / ChEBI",
                                ("metabolite_registry.tsv.gz",),
                                "omicsweave.core.idmap.build_metabolite_registry:main", "metabolites"),
    "reaction_roles": RegistryInfo("reaction_roles", "反应成员角色（G17）：由策展表生成基因表与蛋白表中登记反应的角色分组；行数为登记的反应数",
                                   ("reaction_members.yaml",),
                                   "omicsweave.integration.mappings.reaction_members:main", "",
                                   package="omicsweave.integration.mappings", subdir=""),
}


def _data_dir() -> Path:
    return Path(importlib.import_module("omicsweave.core.idmap").__file__).parent / "data"


def _registry_dir(info: RegistryInfo) -> Path:
    return Path(importlib.import_module(info.package).__file__).parent / info.subdir


def _count_rows(info: RegistryInfo) -> int:
    from .idmap._registry import read_table

    if info.files[0].endswith((".yaml", ".yml")):
        import yaml
        return len(yaml.safe_load((_registry_dir(info) / info.files[0]).read_text(encoding="utf-8")) or {})
    return len(read_table(info.files[0]))


def registry_status() -> pd.DataFrame:
    from .idmap._registry import read_sources

    sources = read_sources() if (_data_dir() / "SOURCES.json").is_file() else {}
    rows = []
    for info in REGISTRIES.values():
        present = all((_registry_dir(info) / f).is_file() for f in info.files)
        n_rows = _count_rows(info) if present else 0
        section = sources.get(info.sources_key, {}) if isinstance(sources.get(info.sources_key), dict) else {}
        # 策展表（没有下载来源）不记生成日期
        generated = section.get("generated_on") or sources.get("generated_on", "") if present and info.sources_key \
            else ""
        rows.append((info.name, info.description, ";".join(info.files), present, generated, n_rows, info.builder))
    return pd.DataFrame(rows, columns=["name", "description", "files", "present", "generated_on", "n_rows", "builder"])


def build(name: str, argv: list[str]) -> int:
    if name not in REGISTRIES:
        raise KeyError(f"未知知识表 {name!r}，可选：{'、'.join(REGISTRIES)}")
    module, function = REGISTRIES[name].builder.split(":")
    return int(getattr(importlib.import_module(module), function)(argv) or 0)


def latest_bundle(experiment: Experiment) -> Path | None:
    from omicsweave.contract.query import list_runs

    runs = [r for r in list_runs(experiment.project.path("results"), experiment.id) if any(r.glob("omics/*/data/*"))]
    return runs[-1] if runs else None


def _normalize_hmdb(value: object) -> set[str]:
    out = set()
    for piece in re.split(r"[;,|\s]+", str(value or "")):
        match = re.fullmatch(r"HMDB0*(\d+)", piece.strip().upper())
        if match:
            out.add(f"HMDB{int(match.group(1)):07d}")
    return out


def _normalize_kegg(value: object) -> set[str]:
    return {p.strip().upper() for p in re.split(r"[;,|\s]+", str(value or "")) if re.fullmatch(r"C\d{5}", p.strip().upper())}


def _metabolite_registry() -> dict[str, dict[str, str]]:
    from .idmap.metabolite_registry import registry_by_name
    return registry_by_name()


def audit_protein_accessions(run_dir: Path) -> list[tuple]:
    from .idmap.mouse_protein_id_mapping import get_uniprot_ids

    rows = []
    for meta in sorted(run_dir.glob("omics/proteomics/data/*/dataset.json")):
        dataset = read_dataset(meta.parent)
        if "uniprot_id" not in dataset.features.columns:
            continue
        for feature, accession in dataset.features["uniprot_id"].items():
            data_ids = {a.strip().split("-")[0] for a in re.split(r"[;,]", str(accession or "")) if a.strip()}
            registry = set(get_uniprot_ids(str(feature)))
            if data_ids and registry and not data_ids & registry:
                rows.append(("protein_accession", "proteomics", meta.parent.name, str(feature), str(feature),
                             ";".join(sorted(data_ids)), ";".join(sorted(registry)), "accession_mismatch"))
    return rows


def audit_metabolite_ids(run_dir: Path) -> list[tuple]:
    registry = _metabolite_registry()
    if not registry:
        return []
    rows = []
    for meta in sorted(run_dir.glob("omics/metabolomics/data/*/dataset.json")):
        dataset = read_dataset(meta.parent)
        name_column = next((c for c in ("Compounds", "feature_name", "compound") if c in dataset.features.columns), None)
        if name_column is None:
            continue
        for feature, row in dataset.features.iterrows():
            entry = registry.get(str(row[name_column]).strip().lower())
            if not entry or entry.get("status") != "matched":
                continue
            for column, kind, normalize in (("HMDB", "hmdb_id", _normalize_hmdb), ("cpd_ID", "kegg_compound",
                                                                                   _normalize_kegg)):
                data_ids = normalize(row.get(column)) if column in row.index else set()
                expected = normalize(entry.get(kind))
                if data_ids and expected and not data_ids & expected:
                    rows.append(("metabolite_ids", "metabolomics", meta.parent.name, str(feature),
                                 str(row[name_column]), ";".join(sorted(data_ids)), ";".join(sorted(expected)),
                                 f"{kind}_mismatch"))
    return rows


# 已核对的不一致：手写表与 RefMet 指向同一化合物的不同条目（立体异构体、异头物或 HMDB 的重复条目），不再报告。
# 键为（手写表中的名称, ID 类型），值为核对结论（2026-11-16，KEGG REST / Metabolomics Workbench / PubChem）。
ACCEPTED_METABOLITE_VARIANTS: dict[tuple[str, str], str] = {
    ("Glucose", "hmdb_id"): "HMDB0000122 与 HMDB0304632 都是 D-葡萄糖",
    ("Fructose-6-phosphate", "hmdb_id"): "HMDB0003971 为 β-D-果糖-6-磷酸，同一化合物的异头物",
    ("Fructose-1,6-bisphosphate", "kegg_compound"): "C05378 为 β-D-果糖-1,6-二磷酸，C00354 的异头物",
    ("3-Phosphoglycerate", "hmdb_id"): "HMDB0000807 为外消旋体，HMDB0060180 为 D 构型",
    ("Isocitrate", "hmdb_id"): "HMDB0001874 为 D-苏式异柠檬酸，HMDB0000193 为未指定构型",
    ("Isocitric acid", "hmdb_id"): "同 Isocitrate",
    ("NADP", "hmdb_id"): "HMDB0000217 与 HMDB0304435 都是 NADP",
    ("beta-Hydroxybutyrate", "hmdb_id"): "HMDB0000011 与 HMDB0000357 都是 (R)-3-羟丁酸",
    ("Carnitine", "kegg_compound"): "C00318 为 L-肉碱（生理构型），C00487 为未指定构型",
    ("Sedoheptulose-7-phosphate", "hmdb_id"): "HMDB0001068 为 D-景天庚酮糖-7-磷酸，HMDB0258206 未指定构型",
    ("Sedoheptulose 7-Phosphate", "hmdb_id"): "同 Sedoheptulose-7-phosphate",
    ("Allantoin", "kegg_compound"): "C02350 为 (S)-尿囊素，C01551 未指定构型",
}


def audit_metabolite_table(include_accepted: bool = False) -> list[tuple]:
    """手写的代谢物 ID 表与 RefMet 注册表不一致的名称；已核对为同一化合物不同条目的（
    :data:`ACCEPTED_METABOLITE_VARIANTS`）默认不报告。"""
    from .idmap.metabolite_id_mapping import METABOLITE_DATABASE_IDS

    registry = _metabolite_registry()
    rows = []
    for name, ids in METABOLITE_DATABASE_IDS.items():
        entry = registry.get(name.strip().lower())
        if not entry or entry.get("status") != "matched":
            continue
        for kind, normalize in (("hmdb_id", _normalize_hmdb), ("kegg_compound", _normalize_kegg)):
            if not include_accepted and (name, kind) in ACCEPTED_METABOLITE_VARIANTS:
                continue
            written, expected = normalize((ids or {}).get(kind)), normalize(entry.get(kind))
            if written and expected and not written & expected:
                rows.append(("metabolite_table", "metabolomics", "", name, name, ";".join(sorted(written)),
                             ";".join(sorted(expected)), f"{kind}_conflict"))
    return rows


def audit_mapping_coverage(run_dir: Path) -> list[tuple]:
    from omicsweave.integration.mappings import REACTION_ENZYME_MAPPING

    enzymes = sorted({e for m in REACTION_ENZYME_MAPPING.values() for e in m.get("direct_enzymes", [])})
    rows = []
    for meta in sorted(run_dir.glob("omics/proteomics/data/*/dataset.json")):
        measured = {str(f).lower() for f in read_dataset(meta.parent).matrix.index}
        found = [e for e in enzymes if e.lower() in measured]
        rows.append(("mapping_coverage", "proteomics", meta.parent.name, "", "direct_enzymes",
                     f"{len(found)}/{len(enzymes)}", "", "coverage"))
    return rows


def audit(experiment: Experiment | None = None, run_dir: str | Path | None = None) -> pd.DataFrame:
    """对照知识表与数据；``run_dir`` 省略时取实验最近一次含数据集的结果包。"""
    from omicsweave.integration.mappings.reaction_members import audit_reaction_roles

    if run_dir is None and experiment is not None:
        run_dir = latest_bundle(experiment)
    rows = audit_metabolite_table() + audit_reaction_roles()
    if run_dir is not None:
        run_dir = Path(run_dir)
        rows = audit_protein_accessions(run_dir) + audit_metabolite_ids(run_dir) + rows + audit_mapping_coverage(run_dir)
    return pd.DataFrame(rows, columns=AUDIT_COLUMNS)


def audit_messages(table: pd.DataFrame, limit: int = 5) -> list[str]:
    """不一致项的汇总（每个检查与组织一条），供 validate 作为警告输出。"""
    problems = table[table["issue"] != "coverage"]
    messages = []
    for (check, tissue), rows in problems.groupby(["check", "tissue"], sort=True):
        if check in ("reaction_roles", "regulator_effect"):
            issues = "、".join(f"{issue} {n}" for issue, n in rows["issue"].value_counts().sort_index().items())
            examples = "、".join(f"{r.feature_id}.{r.feature_name}（{r.issue}：{r.data_value}"
                                + (f" / {r.registry_value}" if r.registry_value else "") + "）"
                                for r in rows.head(limit).itertuples())
            messages.append(f"{check}：{len(rows)} 处待核对（{issues}），如 {examples}")
            continue
        examples = "、".join(f"{r.feature_name}（数据 {r.data_value}，注册表 {r.registry_value}）"
                            for r in rows.head(limit).itertuples())
        where = f" {tissue}" if tissue else ""
        messages.append(f"{check}{where}：{len(rows)} 处与知识表不一致，如 {examples}")
    return messages
