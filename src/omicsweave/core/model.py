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

import re
import shlex
from enum import Enum
from pathlib import PureWindowsPath
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .tissues import merge_tissue_aliases, resolve_tissue, validate_tissue_aliases

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
    REACTION_DOSSIER = "reaction_dossier"   # 反应档案：成员级证据表与格子状态（G13）


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


class ReplicationSource(_Strict):
    """复现组中的一个来源：某个项目（默认为本项目，可为另一个项目目录的相对或绝对路径）中某个实验的某个比较。"""
    project: str | None = None
    experiment: str
    contrast: str
    run: str = "latest"
    tissue_map: dict[str, str] = Field(default_factory=dict)   # 该来源的组织名 → 共同的组织名


class ReplicationSet(_Strict):
    """估计对象相同、可以互相复现的一组比较（G6，``omicsweave replicate``）。"""
    id: str = Field(pattern=ID_PATTERN)
    sources: list[ReplicationSource] = Field(min_length=2)
    # 允许的跨层对应，如 [[proteomics, transcriptomics]]：一个来源的蛋白对另一个来源的 RNA（按基因符号）
    cross_layer: list[list[Omics]] = Field(default_factory=list)
    # 至少在几个来源中同向（部分合取检验的 u），默认为全部来源
    min_sources: int | None = Field(default=None, ge=2)
    alpha: float = Field(0.05, gt=0, lt=1)

    @model_validator(mode="after")
    def _check(self) -> ReplicationSet:
        if self.min_sources is not None and self.min_sources > len(self.sources):
            raise ValueError(f"min_sources（{self.min_sources}）大于来源数（{len(self.sources)}）")
        for pair in self.cross_layer:
            if len(pair) != 2 or pair[0] == pair[1]:
                raise ValueError("cross_layer 的每一项应为两个不同的组学，如 [proteomics, transcriptomics]")
        return self


# omicsweave pipeline 课题脚本命令中可用的占位符
PIPELINE_PLACEHOLDERS = ("project_dir", "experiment", "run", "run_dir", "integration", "integration_dir", "delivery_dir")
PIPELINE_PLACEHOLDER = re.compile(r"\{([^{}]*)\}")


class PipelineScript(_Strict):
    """``omicsweave pipeline`` 在每个实验的 run / integrate 之后运行的课题脚本（如差异图表、关键结果导出）。

    ``command`` 按 shell 规则拆分（路径用 ``/``，含空格的参数加引号）；首项为 ``.py`` 文件时用运行 omicsweave 的
    同一个 Python 解释器执行，路径相对于项目目录。占位符替换为本次流程的值：``{experiment}`` 实验 ID，``{run}`` /
    ``{run_dir}`` 单组学结果包的运行 ID / 目录，``{integration}`` / ``{integration_dir}`` 整合结果包的运行 ID / 目录，
    ``{project_dir}`` 项目目录，``{delivery_dir}`` 本次交付目录。
    """
    command: str = Field(min_length=1)
    name: str | None = None                                      # 默认为脚本文件名（不含扩展名）
    experiments: list[str] = Field(default_factory=list)         # 只对这些实验运行；默认全部
    outputs: list[str] = Field(default_factory=list)             # 交付目录内必须生成的相对文件

    @field_validator("outputs")
    @classmethod
    def _check_outputs(cls, values: list[str]) -> list[str]:
        for value in values:
            path = PureWindowsPath(value)
            if set(PIPELINE_PLACEHOLDER.findall(value)) - {"experiment"}:
                raise ValueError("脚本 outputs 只支持 {experiment} 占位符")
            if not path.parts or path.anchor or ".." in path.parts:
                raise ValueError(f"脚本预期输出应是交付目录内的相对文件：{value!r}")
        return values

    @field_validator("command")
    @classmethod
    def _check_command(cls, value: str) -> str:
        unknown = sorted(set(PIPELINE_PLACEHOLDER.findall(value)) - set(PIPELINE_PLACEHOLDERS))
        if unknown:
            raise ValueError(f"未知占位符 {'、'.join('{' + u + '}' for u in unknown)}，"
                             f"可用：{'、'.join('{' + p + '}' for p in PIPELINE_PLACEHOLDERS)}")
        try:
            tokens = shlex.split(value)
        except ValueError as exc:
            raise ValueError(f"无法拆分命令 {value!r}：{exc}") from None
        if not tokens:
            raise ValueError("命令为空")
        return value


class RetentionSettings(_Strict):
    """交付目录保留与清理策略（默认关闭，不会自动删除任何历史交付或结果包）。"""
    enabled: bool = False
    keep_last: int | None = Field(default=5, ge=1)     # 每个实验保留最近 N 次成功交付；None 表示不按数量淘汰
    drop_failed: bool = True                           # 清理时是否移除失败与超时中断的交付目录
    prune_orphan_bundles: bool = False                 # 清理时是否联动移除未被保留交付引用的孤立 run_* 结果包
    stale_running_hours: float = Field(default=24.0, ge=0)


class PipelineSettings(_Strict):
    """``omicsweave pipeline`` 的标准流程：对每个实验依次 validate → run → integrate → 课题脚本。"""
    experiments: list[str] = Field(default_factory=list)         # 实验 ID；默认项目中的全部实验
    scripts: list[PipelineScript] = Field(default_factory=list)
    retention: RetentionSettings = Field(default_factory=RetentionSettings)


class ProjectConfig(_Strict):
    schema_version: Literal[1] = SCHEMA_VERSION
    project: ProjectInfo
    defaults: ProjectDefaults = Field(default_factory=ProjectDefaults)
    # 相对于项目目录的实验目录；留空时自动查找 paths.experiments 下含 experiment.yaml 的目录
    experiments: list[str] = Field(default_factory=list)
    replication: list[ReplicationSet] = Field(default_factory=list)
    pipeline: PipelineSettings = Field(default_factory=PipelineSettings)


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


_TERM = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?:\s*[:*]\s*[A-Za-z_][A-Za-z0-9_]*)*$")


def formula_terms(formula: str) -> list[tuple[str, ...]]:
    """``~ a + b + a:b`` / ``~ a * b`` → [(a,), (b,), (a, b)]：主效应与交互项，按出现顺序去重。"""
    body = formula.strip()
    if not body.startswith("~"):
        raise ValueError(f"公式应以 ~ 开头：{formula!r}")
    terms: list[tuple[str, ...]] = []
    for part in body[1:].split("+"):
        part = part.strip()
        if not part or not _TERM.match(part):
            raise ValueError(f"公式 {formula!r} 中无法识别的项：{part!r}（支持 a、a:b、a*b）")
        if "*" in part:
            names = [n.strip() for n in part.split("*")]
            expanded = [(n,) for n in names] + [tuple(names)]
        else:
            expanded = [tuple(n.strip() for n in part.split(":"))]
        for term in expanded:
            if term not in terms:
                terms.append(term)
    return terms


class ModelSpec(_Strict):
    """设计感知的统计模型（G3）：在已声明的因素上写公式，报告其中的主效应与交互项。

    ``design.covariates`` 自动加入每个模型。主效应报告为每个非参照水平相对参照水平的合并估计
    （``estimand_type = pooled_main``），交互项为 ``interaction``；与两两比较并列输出，不互相替代。
    """
    id: str = Field(pattern=ID_PATTERN)
    formula: str
    report: list[str] = Field(default_factory=list)    # 报告的项（如 diet、diet:treatment），默认为公式中的全部项
    subset: dict[str, str] = Field(default_factory=dict)
    role: ContrastRole = ContrastRole.OTHER
    omics: list[Omics] = Field(default_factory=list)    # 默认为结果包中的全部组学

    @field_validator("subset", mode="before")
    @classmethod
    def _subset_as_str(cls, value: Any) -> Any:
        return {k: _as_str(v) for k, v in value.items()} if isinstance(value, dict) else value

    @model_validator(mode="after")
    def _check_formula(self) -> ModelSpec:
        terms = [":".join(t) for t in formula_terms(self.formula)]
        unknown = [r for r in self.report if r not in terms]
        if unknown:
            raise ValueError(f"report 中的项不在公式中：{'、'.join(unknown)}（公式中的项：{'、'.join(terms)}）")
        return self

    def terms(self) -> list[tuple[str, ...]]:
        return formula_terms(self.formula)

    def reported(self) -> list[tuple[str, ...]]:
        return [t for t in self.terms() if not self.report or ":".join(t) in self.report]


class Design(_Strict):
    subjects: str = "subjects.tsv"
    samples: str = "samples.tsv"
    factors: dict[str, Factor] = Field(default_factory=dict)
    contrasts: list[Contrast] = Field(default_factory=list)
    reversal_pairs: list[ReversalPair] = Field(default_factory=list)
    stratify_by: list[str] = Field(default_factory=list)
    covariates: list[str] = Field(default_factory=list)
    models: list[ModelSpec] = Field(default_factory=list)
    # 对各差异结果做自适应收缩，增加 lfsr 列（证据画像据此计算方向可信度）
    adaptive_shrinkage: bool = False
    # 实验级组织别名（别名 → 规范组织名），所有单组学运行与整合模式都会套用
    tissue_aliases: dict[str, str] = Field(default_factory=dict)

    @field_validator("tissue_aliases")
    @classmethod
    def _check_tissue_aliases(cls, value: dict[str, str]) -> dict[str, str]:
        return validate_tissue_aliases(value, "design.tissue_aliases")

    @model_validator(mode="after")
    def _check_references(self) -> Design:
        problems: list[str] = []
        known = set(self.factors) | set(self.covariates)
        model_ids = set()
        for model in self.models:
            if model.id in model_ids or model.id in {c.id for c in self.contrasts}:
                problems.append(f"模型 id {model.id} 与其他模型或比较重复")
            model_ids.add(model.id)
            for term in model.terms():
                missing = [name for name in term if name not in known]
                if missing:
                    problems.append(f"模型 {model.id} 的公式引用了未声明的因素：{'、'.join(missing)}")
            for key in model.subset:
                if key not in self.factors:
                    problems.append(f"模型 {model.id} 的 subset 引用了未声明的因素 {key!r}")
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
    # 组学级组织别名（别名 → 规范组织名），合并并覆写 design.tissue_aliases 的同名条目
    tissue_aliases: dict[str, str] = Field(default_factory=dict)

    @field_validator("tissue_aliases")
    @classmethod
    def _check_tissue_aliases(cls, value: dict[str, str]) -> dict[str, str]:
        return validate_tissue_aliases(value, "omics.tissue_aliases")

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


class FeatureSelector(_Strict):
    """派生特征的成员：按特征 ID 或名称、脂质类别（LIPID MAPS 简写，醚脂为 PE-O / PE-P）或 HMDB / KEGG ID 选取。"""
    features: list[str] = Field(default_factory=list)
    lipid_class: list[str] = Field(default_factory=list)
    database_ids: list[str] = Field(default_factory=list)

    def empty(self) -> bool:
        return not (self.features or self.lipid_class or self.database_ids)


DERIVED_TYPES = ("sum", "ratio", "step_ratio", "fraction", "index", "module_score")


class DerivedSpec(_Strict):
    """声明式派生特征（G8），在分析层数据集（log2 尺度）上按样本计算，与测得特征一起做差异分析、逆转与整合。

    - ``sum``：``members`` 在线性尺度上求和后取 log2；
    - ``ratio``：``numerator`` 与 ``denominator`` 之比（各自先求和）；
    - ``step_ratio``：一步反应的产物 / 底物比，``source: network`` 时由 MFA 网络的化学计量自动生成；
    - ``fraction``：``members`` 占 ``denominator``（总量）的比例；
    - ``index``：``members`` 中脂质按丰度加权的平均链长（``index: chain_length``）或双键数（``double_bonds``）；
    - ``module_score``：``members`` 的样本得分（各成员标准化后取平均）。
    """
    id: str | None = Field(default=None, pattern=ID_PATTERN)
    type: Literal["sum", "ratio", "step_ratio", "fraction", "index", "module_score"]
    omics: Omics
    members: FeatureSelector | None = None
    numerator: FeatureSelector | None = None
    denominator: FeatureSelector | None = None
    index: Literal["chain_length", "double_bonds"] | None = None
    source: Literal["network"] | None = None
    min_members: int = Field(1, ge=1)

    @model_validator(mode="after")
    def _check_type(self) -> DerivedSpec:
        need = {"sum": ("members",), "ratio": ("numerator", "denominator"), "fraction": ("members", "denominator"),
                "index": ("members", "index"), "module_score": ("members",)}.get(self.type, ())
        missing = [n for n in need if getattr(self, n) is None or (isinstance(getattr(self, n), FeatureSelector)
                                                                     and getattr(self, n).empty())]
        if missing:
            raise ValueError(f"派生特征类型 {self.type} 需要 {'、'.join(missing)}")
        if self.type == "step_ratio" and self.source is None and (self.numerator is None or self.denominator is None):
            raise ValueError("step_ratio 需要 source: network，或给出 numerator（产物）与 denominator（底物）")
        if self.id is None and not (self.type == "step_ratio" and self.source == "network"):
            raise ValueError(f"派生特征类型 {self.type} 需要 id")
        return self


class ExperimentInfo(_Strict):
    id: str = Field(pattern=ID_PATTERN)
    description: str = ""
    metadata: dict[str, Any] = Field(default_factory=dict)


class ExperimentConfig(_Strict):
    schema_version: Literal[1] = SCHEMA_VERSION
    experiment: ExperimentInfo
    design: Design = Field(default_factory=Design)
    omics: OmicsBlock = Field(default_factory=OmicsBlock)
    derived: list[DerivedSpec] = Field(default_factory=list)
    integration: IntegrationSettings = Field(default_factory=IntegrationSettings)

    @model_validator(mode="after")
    def _check_derived_ids(self) -> ExperimentConfig:
        ids = [d.id for d in self.derived if d.id]
        duplicated = sorted({i for i in ids if ids.count(i) > 1})
        if duplicated:
            raise ValueError(f"派生特征 id 重复：{'、'.join(duplicated)}")
        for omics in Omics:
            own = self.omics.settings(omics).tissue_aliases
            if own or self.design.tissue_aliases:
                merge_tissue_aliases(self.design.tissue_aliases, own, label=f"omics.{omics.value}.tissue_aliases")
        return self

    def tissue_aliases_for(self, omics: Omics | str | None = None) -> dict[str, str]:
        """某个组学生效的组织别名表（实验级 ``design.tissue_aliases`` 与组学级 ``omics.<组学>.tissue_aliases`` 合并）。"""
        if omics is None:
            return dict(self.design.tissue_aliases)
        try:
            key = omics if isinstance(omics, Omics) else Omics(str(omics))
        except ValueError:
            return dict(self.design.tissue_aliases)
        own = self.omics.settings(key).tissue_aliases
        if not own:
            return dict(self.design.tissue_aliases)
        return merge_tissue_aliases(self.design.tissue_aliases, own, label=f"omics.{key.value}.tissue_aliases")

    def canonical_tissue(self, tissue: str, omics: Omics | str | None = None) -> str:
        """按实验级与组学级组织别名表返回规范组织名。"""
        return resolve_tissue(tissue, self.tissue_aliases_for(omics))


def export_json_schema() -> dict[str, dict[str, Any]]:
    """返回 project.yaml 与 experiment.yaml 的 JSON Schema。"""
    return {
        "project": ProjectConfig.model_json_schema(),
        "experiment": ExperimentConfig.model_json_schema(),
    }
