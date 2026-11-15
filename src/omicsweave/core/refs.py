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
- ``mapping_coverage``：反应映射中的直接酶在各组织蛋白组数据中的覆盖率。
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


@dataclass(frozen=True)
class RegistryInfo:
    name: str
    description: str
    files: tuple[str, ...]
    builder: str            # "模块:函数"，接受 argv
    sources_key: str


REGISTRIES: dict[str, RegistryInfo] = {
    "genes": RegistryInfo("genes", "小鼠基因：Ensembl GRCm39 cDNA + UniProt，按 MGI 对齐",
                          ("mouse_gene_registry.tsv.gz", "mouse_gene_aliases.tsv.gz"),
                          "omicsweave.core.idmap.build_registry:main", "ensembl"),
    "proteins": RegistryInfo("proteins", "小鼠蛋白：UniProt Swiss-Prot 登记号与 EC 号（与 genes 一起生成）",
                             ("mouse_protein_registry.tsv.gz",), "omicsweave.core.idmap.build_registry:main", "uniprot"),
    "metabolites": RegistryInfo("metabolites", "代谢物：RefMet 标准名称、类别与 HMDB / KEGG / ChEBI",
                                ("metabolite_registry.tsv.gz",),
                                "omicsweave.core.idmap.build_metabolite_registry:main", "metabolites"),
}


def _data_dir() -> Path:
    return Path(importlib.import_module("omicsweave.core.idmap").__file__).parent / "data"


def registry_status() -> pd.DataFrame:
    from .idmap._registry import read_sources, read_table

    sources = read_sources() if (_data_dir() / "SOURCES.json").is_file() else {}
    rows = []
    for info in REGISTRIES.values():
        present = all((_data_dir() / f).is_file() for f in info.files)
        n_rows = len(read_table(info.files[0])) if present else 0
        section = sources.get(info.sources_key, {}) if isinstance(sources.get(info.sources_key), dict) else {}
        generated = section.get("generated_on") or sources.get("generated_on", "") if present else ""
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


def audit_metabolite_table() -> list[tuple]:
    """手写的代谢物 ID 表与 RefMet 注册表不一致的名称。"""
    from .idmap.metabolite_id_mapping import METABOLITE_DATABASE_IDS

    registry = _metabolite_registry()
    rows = []
    for name, ids in METABOLITE_DATABASE_IDS.items():
        entry = registry.get(name.strip().lower())
        if not entry or entry.get("status") != "matched":
            continue
        for kind, normalize in (("hmdb_id", _normalize_hmdb), ("kegg_compound", _normalize_kegg)):
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
    if run_dir is None and experiment is not None:
        run_dir = latest_bundle(experiment)
    rows = audit_metabolite_table()
    if run_dir is not None:
        run_dir = Path(run_dir)
        rows = audit_protein_accessions(run_dir) + audit_metabolite_ids(run_dir) + rows + audit_mapping_coverage(run_dir)
    return pd.DataFrame(rows, columns=AUDIT_COLUMNS)


def audit_messages(table: pd.DataFrame, limit: int = 5) -> list[str]:
    """不一致项的汇总（每个检查与组织一条），供 validate 作为警告输出。"""
    problems = table[table["issue"] != "coverage"]
    messages = []
    for (check, tissue), rows in problems.groupby(["check", "tissue"], sort=True):
        examples = "、".join(f"{r.feature_name}（数据 {r.data_value}，注册表 {r.registry_value}）"
                            for r in rows.head(limit).itertuples())
        where = f" {tissue}" if tissue else ""
        messages.append(f"{check}{where}：{len(rows)} 处与知识表不一致，如 {examples}")
    return messages
