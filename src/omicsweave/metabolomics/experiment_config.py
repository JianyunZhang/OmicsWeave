"""
代谢组学配置加载器 — 从 config.yaml / groups.yaml 构建 MetabolomicsConfig。

参考 omicsweave/transcriptomics/experiment_config.py 的设计模式。
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from omicsweave.core.project import project_dir


# ---------------------------------------------------------------------------
# 数据类
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MetabolomicsConfig:
    experiment_id: str
    config_path: Path
    output_dir: Path
    # per-tissue sample data: {tissue_name: resolved Path}
    tissue_data_paths: dict[str, Path]
    groups_path: Path
    conditions: dict[str, dict[str, Any]]
    mouse_groups: dict[str, str]
    # optional display
    tissue_labels: dict[str, str] | None = None
    plot_colors: dict[str, Any] | None = None
    # analysis parameters
    pval_cutoff: float = 0.05
    fdr_cutoff: float = 0.05
    log2fc_cutoff: float = 0.58
    vip_cutoff: float = 1.0
    require_vip_for_significance: bool = False
    min_observations_per_group: int = 2
    normalization_method: str = "pqn"
    detection_rate: float = 0.5
    imputation_fraction: float = 0.5
    qc_filter_enabled: bool = True
    qc_rsd_cutoff: float = 0.30
    qc_min_detection_rate: float = 0.80
    batch_correction: str = "none"
    heatmap_top_n: int = 50
    heatmap_selection_method: str = "variance"
    fdr_method: str = "bh"
    plsda_permutation_num: int = 200
    random_seed: int = 42
    ora_min_set_size: int = 3
    ora_max_set_size: int = 500
    pathway_scope: str = "metabolism"
    exclude_overview_pathways: bool = True
    gsea_min_set_size: int = 3
    gsea_max_set_size: int = 500
    gsea_permutation_num: int = 1000
    gsea_rank_metric: str = "welch_t"


# ---------------------------------------------------------------------------
# 内部工具函数
# ---------------------------------------------------------------------------

def repo_root() -> Path:
    """项目目录（实验配置与结果所在仓库），见 :func:`omicsweave.core.project.project_dir`。"""
    return project_dir()


def _is_abs_path_like(p: str) -> bool:
    if not p:
        return False
    return p.startswith("\\\\") or bool(re.match(r"^[A-Za-z]:[\\/]", p))


def _resolve_path(p: str, base_dir: Path, must_exist: bool = True) -> Path:
    if _is_abs_path_like(p):
        return Path(p)
    pr = repo_root() / p
    if pr.exists() or not must_exist:
        return pr
    return (base_dir / p).resolve()


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"YAML 根节点不是 mapping: {path}")
    return data


def _substitute_vars(s: str, vars_map: dict[str, str]) -> str:
    out = s
    for k, v in vars_map.items():
        out = out.replace("${" + k + "}", v)
    return out


# ---------------------------------------------------------------------------
# 默认配置路径
# ---------------------------------------------------------------------------

def default_config_path(experiment_id: str) -> Path:
    root = repo_root() / "experiments" / experiment_id
    p1 = root / "metabolomics_data" / "config.yaml"
    if p1.exists():
        return p1
    return root / "metabolomics" / "config.yaml"


# ---------------------------------------------------------------------------
# 加载器
# ---------------------------------------------------------------------------

def load_metabolomics_config(
    *,
    experiment_id: str | None = None,
    config_path: str | None = None,
    outdir_override: str | None = None,
) -> MetabolomicsConfig:
    exp = experiment_id or os.environ.get("EXPERIMENT_ID")
    if config_path:
        cfg_path = Path(config_path).resolve()
    else:
        if not exp:
            raise ValueError(
                "缺少 experiment id/config，请传 --experiment <id> 或 --config <path>"
            )
        cfg_path = default_config_path(exp)
    cfg = _load_yaml(cfg_path)
    base_dir = cfg_path.parent

    # experiment_id
    if not exp:
        exp = str(cfg.get("experiment_id", "")).strip() or None
    exp = exp or str(cfg.get("experiment_id", "")).strip() or cfg_path.parents[1].name

    # paths
    paths = cfg.get("paths", {}) or {}

    # 解析 sample_data_* 键 → tissue_data_paths
    tissue_data_paths: dict[str, Path] = {}
    for key, val in paths.items():
        if key.startswith("sample_data_"):
            tissue = key[len("sample_data_"):]
            p = str(val).strip()
            p = _substitute_vars(p, {"experiment_id": exp})
            tissue_data_paths[tissue] = _resolve_path(p, base_dir)
    if not tissue_data_paths:
        raise ValueError(f"未在 paths 中找到 sample_data_* 配置: {cfg_path}")

    # output_dir
    output_dir_s = str(paths.get("output_dir", "")).strip()
    if outdir_override:
        output_dir_s = outdir_override
    if not output_dir_s:
        raise ValueError(f"paths.output_dir 缺失: {cfg_path}")
    output_dir_s = _substitute_vars(output_dir_s, {"experiment_id": exp})
    output_dir = _resolve_path(output_dir_s, base_dir, must_exist=False)

    # params
    params = cfg.get("params", {}) or {}
    pval_cutoff = float(params.get("pval_cutoff", 0.05))
    fdr_cutoff = float(params.get("fdr_cutoff", pval_cutoff))
    log2fc_cutoff = float(params.get("log2fc_cutoff", 0.58))
    vip_cutoff = float(params.get("vip_cutoff", 1.0))
    require_vip_for_significance = bool(params.get("require_vip_for_significance", False))
    min_observations_per_group = int(params.get("min_observations_per_group", 2))
    normalization_method = str(params.get("normalization_method", "pqn")).strip().lower()
    detection_rate = float(params.get("detection_rate", 0.5))
    imputation_fraction = float(params.get("imputation_fraction", 0.5))
    qc_filter_enabled = bool(params.get("qc_filter_enabled", True))
    qc_rsd_cutoff = float(params.get("qc_rsd_cutoff", 0.30))
    qc_min_detection_rate = float(params.get("qc_min_detection_rate", 0.80))
    batch_correction = str(params.get("batch_correction", "none")).strip().lower()
    heatmap_top_n = int(params.get("heatmap_top_n", 50))
    heatmap_selection_method = str(params.get("heatmap_selection_method", "variance")).strip().lower()
    fdr_method = str(params.get("fdr_method", "bh"))
    from .data_preprocessing import normalize_fdr_method
    fdr_method = normalize_fdr_method(fdr_method)
    plsda_permutation_num = int(params.get("plsda_permutation_num", 200))
    random_seed = int(params.get("random_seed", 42))
    ora_min_set_size = int(params.get("ora_min_set_size", 3))
    ora_max_set_size = int(params.get("ora_max_set_size", 500))
    pathway_scope = str(params.get("pathway_scope", "metabolism")).strip().lower()
    exclude_overview_pathways = bool(params.get("exclude_overview_pathways", True))
    gsea_min_set_size = int(params.get("gsea_min_set_size", 3))
    gsea_max_set_size = int(params.get("gsea_max_set_size", 500))
    gsea_permutation_num = int(params.get("gsea_permutation_num", 1000))
    gsea_rank_metric = str(params.get("gsea_rank_metric", "welch_t")).strip().lower()

    if normalization_method not in {"pqn", "median", "tic", "none"}:
        raise ValueError("normalization_method 仅支持 pqn/median/tic/none")
    if not 0 < pval_cutoff <= 1 or not 0 < fdr_cutoff <= 1:
        raise ValueError("pval_cutoff 和 fdr_cutoff 必须位于 (0, 1]")
    if not 0 < detection_rate <= 1 or not 0 < qc_min_detection_rate <= 1:
        raise ValueError("检出率阈值必须位于 (0, 1]")
    if not 0 < imputation_fraction <= 1 or qc_rsd_cutoff <= 0:
        raise ValueError("插补比例必须位于 (0, 1]，QC RSD 阈值必须大于 0")
    if min_observations_per_group < 2:
        raise ValueError("min_observations_per_group 至少为 2")
    if batch_correction != "none":
        raise ValueError(
            "代谢组流程尚不支持安全的批次校正；请设 batch_correction: none，"
            "并在上游提供平衡设计或分批分析后做结果级整合"
        )
    if heatmap_selection_method not in {"variance", "anova"}:
        raise ValueError("heatmap_selection_method 仅支持 variance/anova")
    if pathway_scope not in {"metabolism", "all"}:
        raise ValueError("pathway_scope 仅支持 metabolism/all")
    if gsea_rank_metric not in {"welch_t", "signed_logp"}:
        raise ValueError("gsea_rank_metric 仅支持 welch_t/signed_logp")
    if ora_min_set_size < 1 or ora_max_set_size < ora_min_set_size:
        raise ValueError("ORA 集合大小参数无效")
    if gsea_min_set_size < 1 or gsea_max_set_size < gsea_min_set_size:
        raise ValueError("GSEA 集合大小参数无效")
    if plsda_permutation_num < 0 or gsea_permutation_num < 1:
        raise ValueError("置换次数参数无效")

    tissue_labels = params.get("tissue_labels", None)
    if tissue_labels:
        tissue_labels = {str(k): str(v) for k, v in tissue_labels.items()}

    plot_colors = params.get("plot_colors", None)

    # groups
    groups_file = str(cfg.get("groups_file", "groups.yaml")).strip() or "groups.yaml"
    groups_path = _resolve_path(groups_file, base_dir)
    groups = _load_yaml(groups_path)

    conditions = groups.get("conditions", {}) or {}
    mouse_groups = groups.get("mouse_groups", {}) or {}
    mouse_groups = {str(k): str(v) for k, v in mouse_groups.items()}

    if not conditions:
        raise ValueError(f"conditions 为空: {groups_path}")
    if not mouse_groups:
        raise ValueError(f"mouse_groups 为空: {groups_path}")

    return MetabolomicsConfig(
        experiment_id=exp,
        config_path=cfg_path,
        output_dir=output_dir,
        tissue_data_paths=tissue_data_paths,
        groups_path=groups_path,
        conditions=conditions,
        mouse_groups=mouse_groups,
        tissue_labels=tissue_labels,
        plot_colors=plot_colors,
        pval_cutoff=pval_cutoff,
        fdr_cutoff=fdr_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        vip_cutoff=vip_cutoff,
        require_vip_for_significance=require_vip_for_significance,
        min_observations_per_group=min_observations_per_group,
        normalization_method=normalization_method,
        detection_rate=detection_rate,
        imputation_fraction=imputation_fraction,
        qc_filter_enabled=qc_filter_enabled,
        qc_rsd_cutoff=qc_rsd_cutoff,
        qc_min_detection_rate=qc_min_detection_rate,
        batch_correction=batch_correction,
        heatmap_top_n=heatmap_top_n,
        heatmap_selection_method=heatmap_selection_method,
        fdr_method=fdr_method,
        plsda_permutation_num=plsda_permutation_num,
        random_seed=random_seed,
        ora_min_set_size=ora_min_set_size,
        ora_max_set_size=ora_max_set_size,
        pathway_scope=pathway_scope,
        exclude_overview_pathways=exclude_overview_pathways,
        gsea_min_set_size=gsea_min_set_size,
        gsea_max_set_size=gsea_max_set_size,
        gsea_permutation_num=gsea_permutation_num,
        gsea_rank_metric=gsea_rank_metric,
    )


# ---------------------------------------------------------------------------
# 辅助函数
# ---------------------------------------------------------------------------

_DEFAULT_PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52",
    "#8172B3", "#937860", "#DA8BC3", "#8C8C8C",
]


def generate_comparisons(
    conditions: dict[str, Any],
    groups_cfg: dict[str, Any] | None = None,
) -> list[tuple[str, str, str]]:
    """
    生成两两比较列表。
    如果 groups.yaml 中定义了 contrast_groups，则使用用户指定的比较；
    否则自动生成所有两两比较。

    Returns: [(numerator, denominator, comp_name), ...]
    comp_name 格式: '{display_a}_vs_{display_b}'
    """
    if groups_cfg:
        # 使用用户定义的对比组
        comps = []
        for _label, pair in groups_cfg.items():
            if len(pair) == 2:
                a, b = str(pair[0]), str(pair[1])
                da = conditions.get(a, {}).get("display", a)
                db = conditions.get(b, {}).get("display", b)
                comp_name = f"{da}_vs_{db}"
                comps.append((a, b, comp_name))
        if comps:
            return comps

    # 自动生成所有两两比较
    # YAML 插入顺序即默认方向；若标记 reference/control，则始终置于分母。
    keys = list(conditions.keys())
    reference_keys = {
        key for key, value in conditions.items()
        if isinstance(value, dict)
        and (
            bool(value.get("reference", False))
            or str(value.get("role", "")).lower() in {"reference", "control"}
        )
    }
    comps = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            numerator, denominator = a, b
            if a in reference_keys and b not in reference_keys:
                numerator, denominator = b, a
            da = conditions.get(numerator, {}).get("display", numerator)
            db = conditions.get(denominator, {}).get("display", denominator)
            comp_name = f"{da}_vs_{db}"
            comps.append((numerator, denominator, comp_name))
    return comps


def get_group_colors(conditions: dict[str, Any]) -> dict[str, str]:
    colors = {}
    for i, (key, val) in enumerate(sorted(conditions.items())):
        if isinstance(val, dict) and val.get("color"):
            colors[key] = val["color"]
        else:
            colors[key] = _DEFAULT_PALETTE[i % len(_DEFAULT_PALETTE)]
    return colors


def get_tissue_labels(cfg: MetabolomicsConfig) -> dict[str, str]:
    if cfg.tissue_labels:
        return dict(cfg.tissue_labels)
    return {t: t.replace("_", " ").title() for t in cfg.tissue_data_paths}


def get_comp_labels(
    conditions: dict[str, Any],
    comparisons: list[tuple[str, str, str]],
) -> dict[str, str]:
    def _disp(code: str) -> str:
        val = conditions.get(code)
        if isinstance(val, dict) and val.get("display"):
            return val["display"]
        return code
    result = {}
    for a, b, comp_name in comparisons:
        result[comp_name] = f"{_disp(a)} vs {_disp(b)}"
    return result


def get_cond_map(conditions: dict[str, Any]) -> dict[str, str]:
    result = {}
    for key, val in conditions.items():
        if isinstance(val, dict) and val.get("display"):
            result[key] = val["display"]
        else:
            result[key] = key
    return result


_PLOT_COLOR_DEFAULTS: dict[str, Any] = {
    "heatmap_cmap": "RdBu_r",
    "heatmap_vmin": -2,
    "heatmap_vmax": 2,
    "volcano_up": "#D32F2F",
    "volcano_down": "#1565C0",
    "volcano_ns": "#BDBDBD",
    "gsea_positive": "#D32F2F",
    "gsea_negative": "#1565C0",
    "enrichment_cmap": "RdYlBu_r",
}


def get_plot_colors(cfg: MetabolomicsConfig) -> dict[str, Any]:
    raw = cfg.plot_colors or {}
    return {**_PLOT_COLOR_DEFAULTS, **raw}


def config_arg_help() -> str:
    return (
        "使用 experiments/<experiment>/metabolomics_data/config.yaml 配置。"
        "传 --experiment <id> 或 --config <path>。"
    )
