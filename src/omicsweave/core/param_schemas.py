"""参数验证及自省。原始配置保持稀疏，运行器继续按其原有顺序解析继承。"""
from __future__ import annotations

from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, create_model

from . import parameter_defaults as defaults
from .model import ExperimentConfig, IntegrationMode, Omics, ProjectConfig


class Parameters(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True, allow_inf_nan=False)


class DossierSource(Parameters):
    label: str
    experiment: str
    project: str | None = None
    run: str = "latest"
    layers: list[Literal["mfa", "transcriptomics", "proteomics", "metabolomics"]] | None = None
    contrast_map: dict[str, str] | None = None
    tissue_map: dict[str, str] | None = None


Probability = Annotated[float, Field(ge=0, le=1)]
NonNegative = Annotated[int, Field(ge=0)]
FdrMethod = Literal["bh", "fdr_bh", "fdr_by", "bonferroni", "holm", "none"]
PlotFormat = Literal["pdf", "png", "html", "svg", "eps", "tiff", "jpg", "jpeg"]

_ENUMS = {
    "fdr_scope": Literal["per_tissue", "per_contrast"],
    "smd_weighting": Literal["smd", "reliability_smd", "identifiability_smd"],
    "permutation_scale": Literal["absolute", "relative"],
    "gene_aggregation": Literal["lengthScaledTPM", "simple_sum"],
    "protein_aggregation": Literal["log2_mean", "linear_mean"],
    "layout$abstract": Literal["", "widescreen", "cell_press_square"],
    "mode$pathway_map": Literal["atlas", "publication"],
    "journal": Literal["cell_press", "nature"],
    "color_mode": Literal["", "semantic", "pathway"],
    "strip_mode": Literal["omics", "axes"],
    "omics_edge_color": Literal["tinted", "grey"],
    "min_claim": Literal["L1", "L2", "L3"],
    "min_signal": Literal["assessed", "active", "all"],
    "multiple_testing": FdrMethod,
}
_PROBABILITIES = {
    "alpha", "heterogeneity_alpha", "pval_cutoff", "fdr_cutoff", "pvalue_threshold", "fdr_threshold",
    "ora_padj_cutoff", "gsea_padj_cutoff", "min_abs_correlation", "qc_imbalance", "minor_share",
    "direction_stability_min", "min_selection_probability", "min_direction_stability", "ci_level",
}
_INPUT_FIELDS = {"gene_sets", "layout", "sources", "panels"}


def _model(name: str, values: dict, *, path: str = "") -> type[Parameters]:
    fields = {}
    for key, default in values.items():
        nested = f"{path}.{key}".strip(".")
        kind = ("presentation" if path.startswith("pathway_map") or key in {"figures", "heatmap", "volcano"}
                else "input" if key in _INPUT_FIELDS else "method")
        if key == "sources":
            annotation = list[DossierSource]
        elif key == "formats":
            annotation = list[PlotFormat]
        elif key in {"run_gsea", "run_prerank_gsea"}:
            annotation = bool | None  # 显式 null 表示沿用旧别名或运行默认值
        elif key == "gsea_permutations":
            annotation = NonNegative | None  # 历史配置曾记录此字段，现行运行使用 gsea_permutation_num
        elif key in _PROBABILITIES:
            annotation = Probability if default is not None else Probability | None
        elif key == "tissue_map":
            annotation = dict[str, str] | None
        elif key == "omics_cutoffs":
            annotation = dict[Literal["mfa", "transcriptomics", "proteomics", "metabolomics"], dict[str, float]] | None
        elif default is None:
            annotation = float | None if key == "hedges_g_threshold" else str | None
        elif isinstance(default, dict) and default:
            annotation = _model(f"{name}_{key}", default, path=nested)
            if key in {"submaps", "abstract"}:
                annotation = annotation | bool | None
            elif key in {"sample_qc", "flux_groups", "multi_tracer", "stability", "pathway_map"}:
                annotation = annotation | None
        elif isinstance(default, dict):
            annotation = dict[str, Any] | None
        elif isinstance(default, list):
            annotation = list[float] if default and isinstance(default[0], float) else list[str]
        elif type(default) is int:
            annotation = NonNegative if key.endswith(("seed", "replicates", "permutations")) or key.startswith("min_") else int
        else:
            suffix = path.rsplit(".", 1)[-1]
            annotation = _ENUMS.get(f"{key}${suffix}", _ENUMS.get(key, type(default)))
        fields[key] = (annotation, Field(default=default, json_schema_extra={"x-tunable-class": kind}))
    return create_model(name, __base__=Parameters, **fields)


OMICS_SCHEMAS = {
    Omics.MFA: _model("MfaAnalysis", defaults.MFA),
    Omics.TRANSCRIPTOMICS: _model("TranscriptomicsAnalysis", {
        **defaults.TRANSCRIPTOMICS, "qc_thresholds": {}, "on_qc_fail": "", "qc_sample_size": 0,
        "gsea_permutations": None,
    }),
    Omics.PROTEOMICS: _model("ProteomicsAnalysis", {
        **defaults.PROTEOMICS,
        # 旧配置迁移仍保存绘图、填补和 ssGSEA 字段，供原版入口重现。
        "heatmap_top_n": 50, "imputation_shift": 1.8, "imputation_width": 0.3, "gsea_min_set_size": 15,
        "gsea_max_set_size": 500, "gsea_permutation_num": 1000, "ssgsea_permutation_num": 0,
        "gsea_seed": 42, "gsea_threads": 4, "run_ssgsea": False, "gsea_permutations": None,
    }),
    Omics.METABOLOMICS: _model("MetabolomicsAnalysis", {**defaults.METABOLOMICS, "tissue_labels": {}, "plot_colors": {}}),
}
INTEGRATION_SCHEMAS = {
    mode: _model(f"{mode.value.title().replace('_', '')}Parameters", values)
    for mode, values in {
        IntegrationMode.FLUX_ANCHORED: defaults.FLUX_ANCHORED,
        IntegrationMode.PATHWAY_CONVERGENCE: defaults.PATHWAY_CONVERGENCE,
        IntegrationMode.REVERSAL: defaults.REVERSAL,
        IntegrationMode.REACTION_DOSSIER: defaults.REACTION_DOSSIER,
    }.items()
}
PATHWAY_MAP_SCHEMA = _model("PathwayMapParameters", defaults.PATHWAY_MAP, path="pathway_map")


def parameter_errors(config: ExperimentConfig):
    """按 YAML 原始位置报告参数错误，兼顾暂时关闭的已配置参数块。"""
    blocks = [(f"omics.{layer.value}.analysis", model, config.omics.settings(layer).analysis)
              for layer, model in OMICS_SCHEMAS.items()]
    blocks.extend((f"integration.parameters.{mode.value}", INTEGRATION_SCHEMAS[mode], value)
                  for mode, value in config.integration.parameters.items() if mode in INTEGRATION_SCHEMAS)
    for prefix, model, payload in blocks:
        try:
            model.model_validate(payload)
        except ValidationError as exc:
            for row in exc.errors(include_input=False, include_url=False, include_context=False):
                locations = [str(item) for item in row["loc"]]
                where = ".".join([prefix, *locations])
                field = locations[-1] if locations else prefix.rsplit(".", 1)[-1]
                if row["type"] == "extra_forbidden":
                    if len(locations) > 1:
                        message = f"{locations[-2]} 中有无法识别的项：{field}"
                    else:
                        message = f"{prefix} 中有无法识别的参数：{field}"
                elif row["type"] == "bool_type":
                    message = f"{where} 应为 true 或 false（布尔值）"
                else:
                    message = f"{where}：{row['msg']}"
                yield where, message


def schema(target: str) -> dict:
    models = {
        "project": ProjectConfig, "experiment": ExperimentConfig, "pathway_map": PATHWAY_MAP_SCHEMA,
        **{f"omics/{layer.value}": model for layer, model in OMICS_SCHEMAS.items()},
        **{f"integration/{mode.value}": model for mode, model in INTEGRATION_SCHEMAS.items()},
    }
    if target == "targets":
        return {"targets": [*models, "manifest", "tables/differential", "tables/enrichment"]}
    if target == "manifest":
        from omicsweave.contract.manifest import RunManifest

        return RunManifest.model_json_schema()
    if target.startswith("tables/"):
        from omicsweave.contract.tables import DIFFERENTIAL_COLUMNS, ENRICHMENT_COLUMNS

        columns = {"tables/differential": DIFFERENTIAL_COLUMNS, "tables/enrichment": ENRICHMENT_COLUMNS}[target]
        field_types = {"str": "string", "float": "number", "int": "integer", "bool": "boolean"}
        return {"type": "object", "title": target,
                "required": [col.name for col in columns if col.required],
                "properties": {col.name: {"type": [field_types[col.kind], "null"],
                                         "description": col.description} for col in columns}}
    result = models[target].model_json_schema()
    if target == "experiment":
        result["x-omicsweave-parameter-schemas"] = {
            **{f"omics.{layer.value}.analysis": model.model_json_schema()
               for layer, model in OMICS_SCHEMAS.items()},
            **{f"integration.parameters.{mode.value}": model.model_json_schema()
               for mode, model in INTEGRATION_SCHEMAS.items()},
        }
    return result
