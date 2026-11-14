#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
项目与实验配置的数据模型（``project.yaml`` / ``experiment.yaml``）。

字段定义即配置规范：:func:`export_json_schema` 导出的 JSON Schema 可直接作为配置文档。
所有模型都拒绝未声明的键，拼写错误会在加载时报出，而不是被静默忽略。

只依赖配置文件本身就能判断的规则在这里校验（例如比较引用了未声明的因素）；
需要结合样本表或文件系统的规则见 :mod:`omicsweave.core.validate`。
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

SCHEMA_VERSION = 1
ID_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_.-]*$"


class Species(str, Enum):
    MOUSE = "mus_musculus"
    HUMAN = "homo_sapiens"


TAXON_IDS = {Species.MOUSE: 10090, Species.HUMAN: 9606}


class Omics(str, Enum):
    TRANSCRIPTOMICS = "transcriptomics"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    MFA = "mfa"


# 各组学允许的输入层级与格式（v2 设计文档第 6 节）
INPUT_FORMATS: dict[Omics, dict[str, set[str]]] = {
    Omics.TRANSCRIPTOMICS: {"raw": {"fastq_gz"}, "table": {"counts_matrix", "kallisto", "salmon"}},
    Omics.PROTEOMICS: {"raw": {"thermo_raw", "mzml"}, "table": {"diann_report", "protein_matrix"}},
    Omics.METABOLOMICS: {"raw": {"sciex_wiff", "thermo_raw", "mzml"}, "table": {"peak_table"}},
    Omics.MFA: {"raw": set(), "table": {"omicsweave_mfa_xlsx", "legacy_mfa_xlsx", "mfa_tables"}},
}

# 原始层文件的扩展名，用于检查 samples.tsv 的 files 列
RAW_EXTENSIONS: dict[str, tuple[str, ...]] = {
    "fastq_gz": (".fastq.gz", ".fq.gz"),
    "thermo_raw": (".raw",),
    "mzml": (".mzml",),
    "sciex_wiff": (".wiff",),
}


class ContrastRole(str, Enum):
    DISEASE = "disease_effect"
    TREATMENT = "treatment_effect"
    TIME = "time_effect"
    OTHER = "other"


class IntegrationMode(str, Enum):
    FLUX_ANCHORED = "flux_anchored"
    PATHWAY_CONVERGENCE = "pathway_convergence"
    SAMPLE_LEVEL = "sample_level"
    REVERSAL = "reversal"            # design.reversal_pairs 的逐特征逆转判定


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


def _as_str(value: Any) -> Any:
    """YAML 会把 4、14 读成数字；因素水平统一按字符串保存。"""
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, (int, float)):
        return format(value, "g")
    return value


# ── project.yaml ──────────────────────────────────────────────────────────

class ProjectPaths(_Strict):
    experiments: str = "experiments"
    results: str = "analysis_results"
    references: str = "references"
    analysis: str = "analysis"


class ProjectInfo(_Strict):
    id: str = Field(pattern=ID_PATTERN)
    title: str = ""
    species: Species
    paths: ProjectPaths = Field(default_factory=ProjectPaths)
    references: dict[str, str] = Field(default_factory=dict)


class ProjectDefaults(_Strict):
    padj_method: Literal["fdr_bh", "fdr_by", "bonferroni", "holm"] = "fdr_bh"
    alpha: float = Field(0.05, gt=0, lt=1)
    random_seed: int | None = None


class ProjectConfig(_Strict):
    schema_version: Literal[1] = SCHEMA_VERSION
    project: ProjectInfo
    defaults: ProjectDefaults = Field(default_factory=ProjectDefaults)
    # 相对于项目目录的实验目录；留空时自动查找 paths.experiments 下含 experiment.yaml 的目录
    experiments: list[str] = Field(default_factory=list)


# ── experiment.yaml ───────────────────────────────────────────────────────

class Factor(_Strict):
    level: Literal["subject", "sample"]
    type: Literal["categorical", "numeric"] = "categorical"
    levels: list[str] | None = None
    reference: str | None = None
    unit: str | None = None

    @field_validator("levels", mode="before")
    @classmethod
    def _levels_as_str(cls, value: Any) -> Any:
        return [_as_str(v) for v in value] if isinstance(value, list) else value

    @field_validator("reference", mode="before")
    @classmethod
    def _reference_as_str(cls, value: Any) -> Any:
        return _as_str(value)

    @model_validator(mode="after")
    def _check_reference(self) -> Factor:
        if self.levels is not None and len(set(self.levels)) != len(self.levels):
            raise ValueError("levels 中有重复值")
        if self.reference is not None and self.levels is not None and self.reference not in self.levels:
            raise ValueError(f"reference {self.reference!r} 不在 levels 中")
        return self


class Contrast(_Strict):
    id: str = Field(pattern=ID_PATTERN)
    factor: str
    numerator: str
    denominator: str
    subset: dict[str, str] = Field(default_factory=dict)
    role: ContrastRole = ContrastRole.OTHER

    @field_validator("numerator", "denominator", mode="before")
    @classmethod
    def _level_as_str(cls, value: Any) -> Any:
        return _as_str(value)

    @field_validator("subset", mode="before")
    @classmethod
    def _subset_as_str(cls, value: Any) -> Any:
        return {k: _as_str(v) for k, v in value.items()} if isinstance(value, dict) else value

    @model_validator(mode="after")
    def _check_levels(self) -> Contrast:
        if self.numerator == self.denominator:
            raise ValueError("numerator 与 denominator 不能相同")
        if self.factor in self.subset:
            raise ValueError(f"subset 不能再限定比较因素 {self.factor!r}")
        return self


class ReversalPair(_Strict):
    disease: str
    treatment: str


class Design(_Strict):
    subjects: str = "subjects.tsv"
    samples: str = "samples.tsv"
    factors: dict[str, Factor] = Field(default_factory=dict)
    contrasts: list[Contrast] = Field(default_factory=list)
    reversal_pairs: list[ReversalPair] = Field(default_factory=list)
    stratify_by: list[str] = Field(default_factory=list)
    covariates: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _check_references(self) -> Design:
        problems: list[str] = []
        seen: set[str] = set()
        for c in self.contrasts:
            if c.id in seen:
                problems.append(f"比较 id 重复：{c.id}")
            seen.add(c.id)
            factor = self.factors.get(c.factor)
            if factor is None:
                problems.append(f"比较 {c.id} 引用了未声明的因素 {c.factor!r}")
            elif factor.type == "categorical" and factor.levels is not None:
                for side in (c.numerator, c.denominator):
                    if side not in factor.levels:
                        problems.append(f"比较 {c.id} 的水平 {side!r} 不在因素 {c.factor} 的 levels 中")
            for key, value in c.subset.items():
                sub = self.factors.get(key)
                if sub is None:
                    problems.append(f"比较 {c.id} 的 subset 引用了未声明的因素 {key!r}")
                elif sub.type == "categorical" and sub.levels is not None and value not in sub.levels:
                    problems.append(f"比较 {c.id} 的 subset 水平 {value!r} 不在因素 {key} 的 levels 中")
        roles = {c.id: c.role for c in self.contrasts}
        for pair in self.reversal_pairs:
            for contrast_id, expected in ((pair.disease, ContrastRole.DISEASE), (pair.treatment, ContrastRole.TREATMENT)):
                if contrast_id not in roles:
                    problems.append(f"reversal_pairs 引用了不存在的比较 {contrast_id!r}")
                elif roles[contrast_id] is not expected:
                    problems.append(f"reversal_pairs 中的 {contrast_id} 的 role 应为 {expected.value}")
        if problems:
            raise ValueError("；".join(problems))
        return self


class OmicsInput(BaseModel):
    # layout、strandedness、acquisition、source、mode 等组学专有字段原样保留
    model_config = ConfigDict(extra="allow")

    level: Literal["raw", "table"]
    format: str
    files: list[str] = Field(default_factory=list)


class OmicsSettings(_Strict):
    enabled: bool = False
    input: OmicsInput | None = None
    adapter: dict[str, Any] = Field(default_factory=dict)
    analysis: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _check_input(self) -> OmicsSettings:
        if self.enabled and self.input is None:
            raise ValueError("已启用的组学必须提供 input")
        return self


class OmicsBlock(_Strict):
    transcriptomics: OmicsSettings = Field(default_factory=OmicsSettings)
    proteomics: OmicsSettings = Field(default_factory=OmicsSettings)
    metabolomics: OmicsSettings = Field(default_factory=OmicsSettings)
    mfa: OmicsSettings = Field(default_factory=OmicsSettings)

    @model_validator(mode="after")
    def _check_formats(self) -> OmicsBlock:
        problems = []
        for omics in Omics:
            settings: OmicsSettings = getattr(self, omics.value)
            if not settings.enabled or settings.input is None:
                continue
            allowed = INPUT_FORMATS[omics][settings.input.level]
            if settings.input.format not in allowed:
                choices = "、".join(sorted(allowed)) or "无（该组学不支持此层级）"
                problems.append(
                    f"{omics.value} 的 {settings.input.level} 层不支持格式 {settings.input.format!r}，可选：{choices}"
                )
        if problems:
            raise ValueError("；".join(problems))
        return self

    def settings(self, omics: Omics) -> OmicsSettings:
        return getattr(self, omics.value)

    def enabled(self) -> list[Omics]:
        return [o for o in Omics if self.settings(o).enabled]


class IntegrationSettings(_Strict):
    enabled: bool = False
    modes: list[IntegrationMode] = Field(default_factory=list)
    # 各整合模式的参数，如 pathway_convergence: {ora_padj_cutoff: 0.05, gsea_padj_cutoff: 0.25}
    parameters: dict[IntegrationMode, dict[str, Any]] = Field(default_factory=dict)


class ExperimentInfo(_Strict):
    id: str = Field(pattern=ID_PATTERN)
    description: str = ""
    metadata: dict[str, Any] = Field(default_factory=dict)


class ExperimentConfig(_Strict):
    schema_version: Literal[1] = SCHEMA_VERSION
    experiment: ExperimentInfo
    design: Design = Field(default_factory=Design)
    omics: OmicsBlock = Field(default_factory=OmicsBlock)
    integration: IntegrationSettings = Field(default_factory=IntegrationSettings)


def export_json_schema() -> dict[str, dict[str, Any]]:
    """返回 project.yaml 与 experiment.yaml 的 JSON Schema。"""
    return {
        "project": ProjectConfig.model_json_schema(),
        "experiment": ExperimentConfig.model_json_schema(),
    }
