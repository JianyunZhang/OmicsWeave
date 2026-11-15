#!/usr/bin/env python3
"""每个配置项都要有读取它的代码（E7）。

枚举 project.yaml / experiment.yaml 模型的全部字段，每个字段必须落在下面三类之一：

- ``CONSUMERS``：读取它的模块与代码片段，测试检查片段确实出现在模块源码中；
- ``DESCRIPTIVE``：只作记录或说明（写入结果包的 metadata），不驱动分析；
- ``KNOWN_GAPS``：已声明但还没有分析读取的字段，附原因与探针（正则）。测试同时检查探针在包源码中确实不存在，
  有人实现后测试会提示把它移到 ``CONSUMERS``。

新增配置字段时若没有登记，测试失败，避免“声明先于实现”。
"""
from __future__ import annotations

import fnmatch
import importlib
import inspect
import re
import types
import typing
import unittest
from pathlib import Path

from pydantic import BaseModel

import omicsweave
from omicsweave.core.model import ExperimentConfig, ProjectConfig

# 字段（fnmatch 模式）→ 读取它的“模块:代码片段”
CONSUMERS: dict[str, tuple[str, ...]] = {
    "project.project.id": ("omicsweave.core.run:ctx.project.id",),
    "project.project.species": ("omicsweave.core.genesets:project.config.project.species",),
    "project.project.paths.experiments": ('omicsweave.core.config:path("experiments")',),
    "project.project.paths.results": ('omicsweave.core.run:path("results")',),
    "project.project.paths.references": ('omicsweave.core.genesets:path("references")',),
    "project.project.references": ("omicsweave.core.config:project.config.project.references",),
    "project.defaults.alpha": ("omicsweave.core.validate:defaults.alpha",
                               "omicsweave.integration.modes.reversal:defaults.alpha"),
    "project.experiments": ("omicsweave.core.config:self.config.experiments",),
    "experiment.experiment.id": ("omicsweave.core.run:experiment.id",),
    "experiment.design.subjects": ("omicsweave.core.config:config.design.subjects",),
    "experiment.design.samples": ("omicsweave.core.config:config.design.samples",),
    "experiment.design.factors.*.level": ("omicsweave.core.run:f.level",),
    "experiment.design.factors.*.type": ('omicsweave.core.run:.type == "numeric"',),
    "experiment.design.factors.*.levels": ("omicsweave.core.validate:.levels",),
    "experiment.design.contrasts[].id": ("omicsweave.mfa.runner:contrast.id",),
    "experiment.design.contrasts[].factor": ("omicsweave.core.run:contrast.factor",),
    "experiment.design.contrasts[].numerator": ("omicsweave.core.run:contrast.numerator",),
    "experiment.design.contrasts[].denominator": ("omicsweave.core.run:contrast.denominator",),
    "experiment.design.contrasts[].subset": ("omicsweave.core.run:contrast.subset",),
    "experiment.design.contrasts[].role": ("omicsweave.mfa.runner:.role",),
    "experiment.design.reversal_pairs[].disease": ("omicsweave.integration.modes.reversal:p.disease",),
    "experiment.design.reversal_pairs[].treatment": ("omicsweave.integration.modes.reversal:p.treatment",),
    "experiment.design.covariates": ("omicsweave.transcriptomics.runner:design.covariates",
                                     "omicsweave.proteomics.runner:design.covariates",
                                     "omicsweave.stats.models:design_cfg.covariates"),
    "experiment.design.factors.*.reference": ("omicsweave.stats.data:factor.reference",),
    "experiment.design.models[].id": ("omicsweave.stats.models:model.id",),
    "experiment.design.models[].formula": ("omicsweave.stats.models:model.formula",),
    "experiment.design.models[].report": ("omicsweave.stats.models:model.reported()",),
    "experiment.design.models[].subset": ("omicsweave.stats.models:model.subset",),
    "experiment.design.models[].role": ("omicsweave.stats.models:model.role",),
    "experiment.design.models[].omics": ("omicsweave.stats.models:model.omics",),
    "experiment.design.adaptive_shrinkage": ("omicsweave.stats.post:design.adaptive_shrinkage",),
    "experiment.derived[].id": ("omicsweave.stats.derived:spec.id",),
    "experiment.derived[].type": ("omicsweave.stats.derived:spec.type",),
    "experiment.derived[].omics": ("omicsweave.stats.derived:spec.omics",),
    "experiment.derived[].members.*": ("omicsweave.stats.derived:spec.members",),
    "experiment.derived[].numerator.*": ("omicsweave.stats.derived:spec.numerator",),
    "experiment.derived[].denominator.*": ("omicsweave.stats.derived:spec.denominator",),
    "experiment.derived[].*.features": ("omicsweave.stats.derived:selector.features",),
    "experiment.derived[].*.lipid_class": ("omicsweave.stats.derived:selector.lipid_class",),
    "experiment.derived[].*.database_ids": ("omicsweave.stats.derived:selector.database_ids",),
    "experiment.derived[].index": ("omicsweave.stats.derived:spec.index",),
    "experiment.derived[].source": ("omicsweave.stats.derived:spec.source",),
    "experiment.derived[].min_members": ("omicsweave.stats.derived:spec.min_members",),
    "project.replication[].id": ("omicsweave.integration.replication:spec.id",),
    "project.replication[].sources[].project": ("omicsweave.integration.replication:source.project",),
    "project.replication[].sources[].experiment": ("omicsweave.integration.replication:source.experiment",),
    "project.replication[].sources[].contrast": ("omicsweave.integration.replication:source.contrast",),
    "project.replication[].sources[].run": ("omicsweave.integration.replication:source.run",),
    "project.replication[].sources[].tissue_map": ("omicsweave.integration.replication:source.tissue_map",),
    "project.replication[].cross_layer": ("omicsweave.integration.replication:spec.cross_layer",),
    "project.replication[].min_sources": ("omicsweave.integration.replication:spec.min_sources",),
    "project.replication[].alpha": ("omicsweave.integration.replication:spec.alpha",),
    "experiment.omics.*.enabled": ("omicsweave.core.run:config.omics.enabled()",),
    "experiment.omics.*.input.level": ('omicsweave.core.validate:spec.level == "raw"',),
    "experiment.omics.*.input.format": ("omicsweave.mfa.runner:input.format",
                                        "omicsweave.transcriptomics.runner:input.format"),
    "experiment.omics.*.input.files": ("omicsweave.mfa.runner:spec.files",
                                       "omicsweave.transcriptomics.runner:spec.files"),
    "experiment.omics.*.adapter": ("omicsweave.transcriptomics.runner:settings.adapter",),
    "experiment.omics.*.analysis": ("omicsweave.mfa.runner:settings.analysis",
                                    "omicsweave.metabolomics.runner:metabolomics.analysis"),
    "experiment.integration.enabled": ("omicsweave.core.validate:integration.enabled",),
    "experiment.integration.modes": ("omicsweave.integration.integrate:settings.modes",),
    "experiment.integration.parameters": ("omicsweave.integration.integrate:settings.parameters",),
}

# 只作记录或说明的字段 → 说明
DESCRIPTIVE: dict[str, str] = {
    "*.schema_version": "配置格式版本，由模型的 Literal 校验",
    "project.project.title": "项目标题，记录用",
    "project.project.paths.analysis": "课题脚本与报告的目录约定，init 创建，流程不写入",
    "experiment.experiment.description": "实验说明，写入结果包 metadata/experiment.json",
    "experiment.experiment.metadata": "实验背景信息（饲料、示踪剂等），写入结果包",
    "experiment.design.factors.*.unit": "数值型因素的单位，记录用",
}

# 已声明、尚无分析读取的字段 → (原因, 探针正则：在包源码中不应出现)
KNOWN_GAPS: dict[str, tuple[str, str]] = {
    "project.defaults.padj_method": ("各组学用自己的多重检验校正设置，项目默认值未接入", r"defaults\.padj_method"),
    "project.defaults.random_seed": ("各组学用自己的随机种子，项目默认值未接入", r"defaults\.random_seed"),
    "experiment.design.stratify_by": ("各组学按 tissue 分层是写死的，只有 validate 读取该字段", r"stratify_by"),
}
# 探针的扫描范围之外：模型定义、模板、迁移与校验
_PROBE_EXCLUDE = ("core/model.py", "core/scaffold.py", "core/migrate.py", "core/validate.py")
# 字段级已有读取、但部分组学未接入 → (原因, 未接入的组学运行器模块, 探针)
PARTIAL: dict[str, tuple[str, tuple[str, ...], str]] = {
    "experiment.design.covariates": ("代谢组与 MFA 的两两比较不读取协变量；公式模型（design.models）对所有组学都加入协变量",
                                     ("omicsweave.metabolomics.runner", "omicsweave.mfa.runner"), "covariates"),
}


def _unwrap(annotation):
    """去掉 Optional，返回 (容器标记, 内层类型)：'' / '[]' / '.*'。"""
    origin = typing.get_origin(annotation)
    if origin in (typing.Union, types.UnionType):
        args = [a for a in typing.get_args(annotation) if a is not type(None)]
        return _unwrap(args[0]) if len(args) == 1 else ("", annotation)
    if origin is list:
        return "[]", typing.get_args(annotation)[0]
    if origin is dict:
        return ".*", typing.get_args(annotation)[1]
    return "", annotation


def field_paths(model: type[BaseModel], prefix: str) -> list[str]:
    paths = []
    for name, info in model.model_fields.items():
        marker, inner = _unwrap(info.annotation)
        path = f"{prefix}.{name}"
        if inspect.isclass(inner) and issubclass(inner, BaseModel):
            paths += field_paths(inner, path + marker)
        else:
            paths.append(path)
    return paths


def all_paths() -> list[str]:
    paths = field_paths(ProjectConfig, "project") + field_paths(ExperimentConfig, "experiment")
    # 四个组学的设置结构相同，按 omics.* 登记
    return sorted({re.sub(r"^experiment\.omics\.(transcriptomics|proteomics|metabolomics|mfa)\.",
                          "experiment.omics.*.", p) for p in paths})


def _matches(path: str, pattern: str) -> bool:
    """fnmatch 中的 [] 是字符集；路径里的列表标记 [] 按字面匹配。"""
    return fnmatch.fnmatchcase(path, pattern.replace("[]", "[[][]]"))


def _source(module: str) -> str:
    return Path(importlib.import_module(module).__file__).read_text(encoding="utf-8")


def _package_sources() -> dict[str, str]:
    root = Path(omicsweave.__file__).parent
    return {p.relative_to(root).as_posix(): p.read_text(encoding="utf-8") for p in root.rglob("*.py")
            if not p.relative_to(root).as_posix().endswith(_PROBE_EXCLUDE)}


class ConfigConsumerTest(unittest.TestCase):
    def test_every_field_is_registered(self) -> None:
        registries = [*CONSUMERS, *DESCRIPTIVE, *KNOWN_GAPS]
        unregistered = [p for p in all_paths() if not any(_matches(p, pattern) for pattern in registries)]
        self.assertEqual(unregistered, [], "新配置字段需要在 CONSUMERS / DESCRIPTIVE / KNOWN_GAPS 中登记")

    def test_no_stale_entries(self) -> None:
        paths = all_paths()
        stale = [pattern for pattern in [*CONSUMERS, *DESCRIPTIVE, *KNOWN_GAPS]
                 if not any(_matches(p, pattern) for p in paths)]
        self.assertEqual(stale, [])

    def test_consumers_read_the_field(self) -> None:
        for pattern, consumers in CONSUMERS.items():
            for consumer in consumers:
                module, token = consumer.split(":", 1)
                with self.subTest(field=pattern, consumer=module):
                    self.assertTrue(token in _source(module), f"{module} 中找不到 {token!r}")

    def test_known_gaps_are_still_gaps(self) -> None:
        sources = _package_sources()
        for pattern, (reason, probe) in KNOWN_GAPS.items():
            hits = sorted(path for path, text in sources.items() if re.search(probe, text))
            with self.subTest(field=pattern):
                self.assertEqual(hits, [], f"{pattern} 已有读取（{reason}），请移到 CONSUMERS")

    def test_partial_consumers(self) -> None:
        for pattern, (reason, modules, probe) in PARTIAL.items():
            self.assertIn(pattern, CONSUMERS)
            for module in modules:
                with self.subTest(field=pattern, module=module):
                    self.assertFalse(probe in _source(module), f"{module} 已读取 {pattern}，请更新 PARTIAL（{reason}）")


if __name__ == "__main__":
    unittest.main()
