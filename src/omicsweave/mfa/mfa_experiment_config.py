#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 实验配置加载器

仿照 transcriptomics 的 experiment_config.py 模式，支持两类 MFA 配置布局：
1) 旧版：experiments/<id>/mfa_data/config.yaml + groups.yaml
2) 新版：experiments/<id>/mfa_data/experiment.yaml + <protocol>/<batch>/config.yaml

Usage:
    from mfa_experiment_config import load_mfa_config, MFAConfig
    cfg = load_mfa_config(experiment_id="20250425_cimr_wt_mice_dio")
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from omicsweave.core.project import project_dir


@dataclass(frozen=True)
class MFAConfig:
    experiment_id: str
    config_path: Path
    output_root: Path
    output_dir: Path
    protocol_id: str | None
    batch_id: str | None
    flux_file: Path | None
    mid_file: Path | None
    glycemia_file: Path | None
    comparisons: list[tuple[str, str]] | None
    condition_colors: dict[str, str] | None
    tissue_colors: dict[str, str] | None
    pathway_groups: dict[str, list[str]] | None
    analysis: dict[str, Any]
    conditions: dict[str, dict[str, Any]]
    sample_groups: dict[str, str]
    # groups.yaml: mouse_groups（动物编号 → 条件）。有的实验 sample_groups 以
    # MFA 样本名为键（20250613_cimr_wt_mice_dio_glp），动物编号只在这里。
    mouse_groups: dict[str, str]
    infusion_qc: dict[str, Any]
    batch_metadata: dict[str, Any]
    tracer: dict[str, Any]


# ── 公共工具函数（复用 transcriptomics 的逻辑） ──────────────────────────

def repo_root() -> Path:
    """项目目录（实验配置与结果所在仓库），见 :func:`omicsweave.core.project.project_dir`。"""
    return project_dir()


def _is_abs_path_like(p: str) -> bool:
    if not p:
        return False
    return p.startswith("\\\\") or bool(re.match(r"^[A-Za-z]:[\\/]", p))


def _resolve_path(p: str, base_dir: Path) -> Path:
    if _is_abs_path_like(p):
        return Path(p)
    pr = repo_root() / p
    if pr.exists():
        return pr
    return (base_dir / p).resolve()


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Invalid YAML root (expected mapping): {path}")
    return data


def _substitute_vars(s: str, vars_map: dict[str, str]) -> str:
    out = s
    for k, v in vars_map.items():
        out = out.replace("${" + k + "}", v)
    return out


def _deep_merge(base: dict[str, Any], override: dict[str, Any] | None) -> dict[str, Any]:
    """递归合并配置字典，保留 base 中未被覆盖的默认值。"""
    merged = dict(base)
    if not override:
        return merged
    for key, val in override.items():
        if isinstance(val, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], val)
        else:
            merged[key] = val
    return merged


DEFAULT_ANALYSIS_CONFIG: dict[str, Any] = {
    "preprocessing": {
        "use_net_flux": True,
        "iqr_filter": True,
        "zero_threshold": 1.0e-10,
        "scale_method": "zscore",
        "missing_value": "drop",
    },
    "differential_flux": {
        "statistic": "auto",
        "pvalue_column": "pvalue_mannwhitney",
        "multiple_testing": "fdr_bh",
        "pvalue_threshold": 0.05,
        "fdr_threshold": 0.10,
        "log2fc_threshold": 0.5,
        "cohen_d_threshold": 0.5,
    },
    # 盒约束贴边判定（df_analysis.sample_flux_matrix）：每只动物的代表净通量
    # |net_flux| 落在 bound_values（及从拟合点推断的盒约束上下界）± bound_tolerance
    # 内记为 at_bound；≥ bound_pinned_min_fraction 的动物贴边的通量标为
    # bound_pinned，优先于拟合 IQR 分档。
    "identifiability": {
        "bound_values": [0.0, 1.0],
        "bound_tolerance": 0.05,
        "infer_bounds_from_fits": True,
        "bound_pinned_min_fraction": 0.5,
    },
    "heatmap": {
        "top_n": 50,
        "pvalue_threshold": 0.05,
        "rank_by": "auto",
    },
    "volcano": {
        "log2fc_threshold": 0.5,
        "pvalue_threshold": 0.05,
        "cohen_d_threshold": 0.5,
        "label_top_n": 10,
    },
}


# ── 配置查找与加载 ──────────────────────────────────────────────────────

def default_mfa_config_path(experiment_id: str) -> Path:
    """
    返回默认 MFA 配置路径。

    兼容两种布局：
    1) 旧版：experiments/<id>/mfa_data/config.yaml
    2) 新版：experiments/<id>/mfa_data/experiment.yaml + <protocol>/<batch>/config.yaml
    """
    mfa_dir = repo_root() / "experiments" / experiment_id / "mfa_data"

    legacy_cfg = mfa_dir / "config.yaml"
    if legacy_cfg.exists():
        return legacy_cfg

    experiment_cfg = mfa_dir / "experiment.yaml"
    if not experiment_cfg.exists():
        return legacy_cfg

    exp_doc = _load_yaml(experiment_cfg)

    defaults = exp_doc.get("defaults", {}) or {}
    active_protocol = str(defaults.get("active_protocol", "")).strip()
    # active_batch 允许写成单值或列表（多批次实验常写列表）。
    # 早期实现直接 str() 整个列表，得到 "['a', 'b']" 这种不存在的目录名，
    # 于是静默落到后面的兜底扫描——用户声明的默认批次实际上被忽略了。
    raw_active_batch = defaults.get("active_batch", "")
    if isinstance(raw_active_batch, (list, tuple)):
        active_batches = [str(b).strip() for b in raw_active_batch if str(b).strip()]
    else:
        active_batches = [str(raw_active_batch).strip()] if str(raw_active_batch).strip() else []
    if active_protocol:
        for batch_id in active_batches:
            candidate = mfa_dir / active_protocol / batch_id / "config.yaml"
            if candidate.exists():
                return candidate

    tracer_protocols = exp_doc.get("tracer_protocols", {}) or {}
    if isinstance(tracer_protocols, dict):
        for protocol_id, spec in tracer_protocols.items():
            if not isinstance(spec, dict):
                continue
            if spec.get("enabled", True) is False:
                continue
            for batch_id in spec.get("batches", []) or []:
                candidate = mfa_dir / str(protocol_id) / str(batch_id) / "config.yaml"
                if candidate.exists():
                    return candidate

    # 最后兜底：扫描第一个存在的 <protocol>/<batch>/config.yaml
    for candidate in sorted(mfa_dir.glob("*/*/config.yaml")):
        return candidate

    return legacy_cfg


def _infer_protocol_batch(cfg: dict[str, Any], cfg_path: Path) -> tuple[str | None, str | None]:
    """从 config 内容或路径推断 protocol_id / batch_id。"""
    protocol_id: str | None = None
    batch_id: str | None = None

    isotope_protocol = cfg.get("isotope_protocol", {}) or {}
    if isinstance(isotope_protocol, dict):
        protocol_id = str(isotope_protocol.get("protocol_id", "")).strip() or None
        batch_id = str(isotope_protocol.get("batch_id", "")).strip() or None

    if protocol_id and batch_id:
        return protocol_id, batch_id

    parts = list(cfg_path.parts)
    idx = None
    for i, p in enumerate(parts):
        if p.lower() == "mfa_data":
            idx = i
            break
    if idx is not None and len(parts) >= idx + 4:
        # .../mfa_data/<protocol>/<batch>/config.yaml
        if parts[-1].lower() == "config.yaml":
            protocol_id = protocol_id or parts[idx + 1]
            batch_id = batch_id or parts[idx + 2]

    return protocol_id, batch_id


def _infer_experiment_id_from_path(cfg_path: Path) -> str | None:
    """从 .../experiments/<experiment_id>/mfa_data/.../config.yaml 推断实验 ID。"""
    parts = list(cfg_path.parts)
    for i, p in enumerate(parts):
        if p.lower() == "experiments" and i + 1 < len(parts):
            return parts[i + 1]
    for i, p in enumerate(parts):
        if p.lower() == "mfa_data" and i - 1 >= 0:
            return parts[i - 1]
    return None


def _split_output_root(path: Path) -> tuple[Path, bool]:
    """
    拆分输出根目录。

    返回 (output_root, is_run_path)：
    - is_run_path=True 表示 path 已包含 runs/ 层级（不再追加 runs/<protocol>/<batch>）
    - is_run_path=False 表示 path 是 mfa_results 根目录
    """
    parts = list(path.parts)
    for i, p in enumerate(parts):
        if p.lower() == "runs":
            if i == 0:
                return path, True
            return Path(*parts[:i]), True
    return path, False


def load_mfa_config(
    *,
    experiment_id: str | None = None,
    config_path: str | None = None,
) -> MFAConfig:
    """
    加载 MFA 实验配置。

    优先使用 config_path，否则用 experiment_id 查找默认位置。
    """
    exp = experiment_id or os.environ.get("EXPERIMENT_ID")

    if config_path:
        cfg_path = Path(config_path).resolve()
    else:
        if not exp:
            raise ValueError(
                "缺少实验 ID 或配置路径。请传入 --experiment <id> 或 --config <path>。"
            )
        cfg_path = default_mfa_config_path(exp)

    if not cfg_path.exists():
        raise FileNotFoundError(f"配置文件不存在: {cfg_path}")

    cfg = _load_yaml(cfg_path)
    base_dir = cfg_path.parent

    # experiment_id: 参数 > config 中的值 > 路径推断
    if not exp:
        exp = str(cfg.get("experiment_id", "")).strip() or None
    exp = exp or str(cfg.get("experiment_id", "")).strip() or _infer_experiment_id_from_path(cfg_path)
    if not exp:
        raise ValueError(f"无法从配置推断 experiment_id: {cfg_path}")

    # ── paths ──
    paths = cfg.get("paths", {}) or {}

    flux_file_s = str(paths.get("flux_file", "")).strip()
    mid_file_s = str(paths.get("mid_file", "")).strip()
    output_dir_s = str(paths.get("output_dir", "")).strip()
    protocol_id, batch_id = _infer_protocol_batch(cfg, cfg_path)

    # 新版布局：paths 可在 batch.yaml 或从 experiment.yaml 推导
    if not output_dir_s:
        output_dir_s = _resolve_output_dir_from_experiment(exp, protocol_id, batch_id)

    if not output_dir_s:
        raise ValueError(f"paths.output_dir 缺失: {cfg_path}")

    # 新版布局：data_files 可在 batch.yaml 中
    batch_file_ref = cfg.get("batch_file") or paths.get("batch_file")
    batch_data = _load_batch_data_files(cfg_path, batch_file_ref)
    if not flux_file_s:
        flux_file_s = batch_data.get("flux_file", "")
    if not mid_file_s:
        mid_file_s = batch_data.get("mid_file", "")

    output_dir_s = _substitute_vars(output_dir_s, {
        "experiment_id": exp,
        "protocol_id": protocol_id or "",
        "batch_id": batch_id or "",
    })
    output_dir_base = _resolve_path(output_dir_s, base_dir)
    output_root, is_run_path = _split_output_root(output_dir_base)
    output_dir = output_dir_base
    if protocol_id and batch_id and not is_run_path:
        output_dir = output_dir_base / "runs" / protocol_id / batch_id

    vars_map = {
        "experiment_id": exp,
        "protocol_id": protocol_id or "",
        "batch_id": batch_id or "",
        "output_root": str(output_root),
        "output_dir": str(output_dir),
    }

    flux_file = _resolve_path(_substitute_vars(flux_file_s, {
        **vars_map,
    }), base_dir) if flux_file_s else None

    mid_file = _resolve_path(_substitute_vars(mid_file_s, {
        **vars_map,
    }), base_dir) if mid_file_s else None

    glycemia_file_s = str(
        paths.get("glycemia_record_file", "")
        or batch_data.get("glycemia_record_file", "")
    ).strip()
    glycemia_file = _resolve_path(
        _substitute_vars(glycemia_file_s, {**vars_map}), base_dir
    ) if glycemia_file_s else None

    # ── groups.yaml ──
    groups_file = str(cfg.get("groups_file", "groups.yaml")).strip() or "groups.yaml"
    groups_file = _substitute_vars(groups_file, vars_map)
    groups_path = _resolve_path(groups_file, base_dir)
    groups = _load_yaml(groups_path)

    conditions = groups.get("conditions", {}) or {}
    sample_groups = groups.get("sample_groups", {}) or {}
    sample_groups = {str(k): str(v).lower() for k, v in sample_groups.items()}
    # 动物编号 → 条件；YAML 里写成数字的编号（86）统一转成字符串，与血糖表一致
    mouse_groups = groups.get("mouse_groups", {}) or {}
    mouse_groups = {str(k).strip(): str(v).lower() for k, v in mouse_groups.items()}

    if not conditions:
        raise ValueError(f"conditions 缺失或为空: {groups_path}")

    # ── comparisons（可选，不配则自动生成） ──
    comparisons_raw = cfg.get("comparisons", None) or groups.get("comparisons", None)
    comparisons = _parse_comparisons(comparisons_raw, conditions) if comparisons_raw else None

    # ── 分析参数：experiment.yaml → protocol.yaml → batch config.yaml ──
    # README 与三级 YAML 都声明了这条覆盖链，但此前只有 batch 级的
    # analysis 段被读取，experiment.yaml / protocol.yaml 的
    # analysis_defaults 从未生效——批次 config 写 `analysis: {}` 时
    # 实际用的是代码内置默认值，而不是实验里配置的阈值。
    analysis = dict(DEFAULT_ANALYSIS_CONFIG)
    for layer_name, layer_cfg in _analysis_default_layers(cfg_path, protocol_id):
        analysis = _deep_merge(analysis, layer_cfg)
    analysis = _deep_merge(analysis, cfg.get("analysis", None))
    _validate_analysis_config(analysis)

    # ── 灌注稳态 QC 阈值：protocol.yaml 提供默认，batch 可覆盖 ──
    infusion_qc = _load_infusion_qc(cfg_path, protocol_id, batch_file_ref)
    batch_metadata = _load_batch_metadata(cfg_path, batch_file_ref)
    tracer = _load_tracer(cfg_path, protocol_id)

    # ── 颜色与通路分组（可选） ──
    condition_colors = cfg.get("condition_colors", None)
    tissue_colors = cfg.get("tissue_colors", None)
    pathway_groups = cfg.get("pathway_groups", None)

    # 如果 groups.yaml 中 conditions 有 color 但未单独配 condition_colors，
    # 自动从 conditions 中提取
    if not condition_colors:
        condition_colors = _extract_colors_from_conditions(conditions)

    return MFAConfig(
        experiment_id=exp,
        config_path=cfg_path,
        output_root=output_root,
        output_dir=output_dir,
        protocol_id=protocol_id,
        batch_id=batch_id,
        flux_file=flux_file,
        mid_file=mid_file,
        glycemia_file=glycemia_file,
        comparisons=comparisons,
        condition_colors=condition_colors,
        tissue_colors=tissue_colors,
        pathway_groups=pathway_groups,
        analysis=analysis,
        conditions=conditions,
        sample_groups=sample_groups,
        mouse_groups=mouse_groups,
        infusion_qc=infusion_qc,
        batch_metadata=batch_metadata,
        tracer=tracer,
    )


# ── 辅助函数 ─────────────────────────────────────────────────────────────

def _find_mfa_data_dir(cfg_path: Path) -> Path | None:
    """从批次 config.yaml 向上找到 mfa_data 目录。"""
    for parent in cfg_path.parents:
        if parent.name.lower() == "mfa_data":
            return parent
    return None


def _load_batch_doc(cfg_path: Path, batch_file_ref: str | None) -> dict[str, Any]:
    """加载批次 batch.yaml 全文（找不到时返回空 dict）。"""
    if not batch_file_ref:
        return {}
    batch_path = (cfg_path.parent / batch_file_ref).resolve()
    if not batch_path.exists():
        return {}
    try:
        return _load_yaml(batch_path)
    except Exception:
        return {}


def _load_batch_metadata(cfg_path: Path, batch_file_ref: str | None) -> dict[str, Any]:
    """批次检测元数据（仪器、提取方法、检测日期等）。

    这些字段此前只是文档，从不被代码读取。现在它们进入运行溯源文件，
    用于判断两个批次是否真的可以合并。
    """
    doc = _load_batch_doc(cfg_path, batch_file_ref)
    metadata = dict(doc.get("batch_metadata", {}) or {})
    scope = doc.get("sample_scope", {}) or {}
    if scope:
        metadata["sample_scope"] = scope
    return metadata


def _load_tracer(cfg_path: Path, protocol_id: str | None) -> dict[str, Any]:
    """protocol.yaml 的 tracer 段（示踪剂种类、标记位置、富集度）。"""
    mfa_dir = _find_mfa_data_dir(cfg_path)
    if mfa_dir is None or not protocol_id:
        return {}
    protocol_yaml = mfa_dir / protocol_id / "protocol.yaml"
    if not protocol_yaml.exists():
        return {}
    try:
        return dict(_load_yaml(protocol_yaml).get("tracer", {}) or {})
    except Exception:
        return {}


def _load_infusion_qc(
    cfg_path: Path,
    protocol_id: str | None,
    batch_file_ref: str | None,
) -> dict[str, Any]:
    """灌注血糖稳态 QC 阈值：protocol.yaml 默认 ← batch.yaml 覆盖。

    ``protocol.yaml: infusion_qc_defaults`` 与 ``batch.yaml: batch_qc`` 早已
    在配置里声明，但此前没有任何代码读取它们——稳态假设从未被检验过。
    """
    config: dict[str, Any] = {}
    mfa_dir = _find_mfa_data_dir(cfg_path)
    if mfa_dir is not None and protocol_id:
        protocol_yaml = mfa_dir / protocol_id / "protocol.yaml"
        if protocol_yaml.exists():
            doc = _load_yaml(protocol_yaml)
            config = dict(doc.get("infusion_qc_defaults", {}) or {})

    batch_doc = _load_batch_doc(cfg_path, batch_file_ref)
    batch_qc = batch_doc.get("batch_qc", {}) or {}
    if not batch_qc.get("inherit_protocol_glycemia_qc", True):
        config = {}
    config = _deep_merge(config, batch_doc.get("infusion_qc", None))
    if "on_qc_fail" in batch_qc:
        config["on_fail"] = batch_qc["on_qc_fail"]
    return config


def _analysis_default_layers(
    cfg_path: Path,
    protocol_id: str | None,
) -> list[tuple[str, dict[str, Any]]]:
    """按 experiment → protocol 顺序返回 ``analysis_defaults`` 覆盖层。

    越靠后的层优先级越高；批次级 ``analysis`` 由调用方最后合并。
    """
    layers: list[tuple[str, dict[str, Any]]] = []
    mfa_dir = _find_mfa_data_dir(cfg_path)
    if mfa_dir is None:
        return layers

    experiment_yaml = mfa_dir / "experiment.yaml"
    if experiment_yaml.exists():
        defaults = _load_yaml(experiment_yaml).get("analysis_defaults") or {}
        if defaults:
            layers.append(("experiment", defaults))

    if protocol_id:
        protocol_yaml = mfa_dir / protocol_id / "protocol.yaml"
        if protocol_yaml.exists():
            defaults = _load_yaml(protocol_yaml).get("analysis_defaults") or {}
            if defaults:
                layers.append(("protocol", defaults))

    return layers


def _resolve_output_dir_from_experiment(
    experiment_id: str,
    protocol_id: str | None,
    batch_id: str | None,
) -> str:
    """从 experiment.yaml 的 defaults.output_root 推导 output_dir。

    新版布局中 batch config.yaml 不含 paths.output_dir，
    需要从 experiment.yaml 的 defaults.output_root + runs/<protocol>/<batch> 推导。

    Returns
    -------
    str: 推导出的 output_dir 路径字符串，无法推导时返回空字符串。
    """
    mfa_dir = repo_root() / "experiments" / experiment_id / "mfa_data"
    experiment_cfg = mfa_dir / "experiment.yaml"

    if not experiment_cfg.exists():
        return ""

    exp_doc = _load_yaml(experiment_cfg)
    defaults = exp_doc.get("defaults", {}) or {}
    raw = defaults.get("output_root", "")
    if not raw:
        return ""

    resolved = _substitute_vars(str(raw), {"experiment_id": experiment_id})
    output_root = _resolve_path(resolved, mfa_dir)

    if protocol_id and batch_id:
        return str(output_root / "runs" / protocol_id / batch_id)
    return str(output_root)


def _load_batch_data_files(
    cfg_path: Path,
    batch_file_ref: str | None,
) -> dict[str, str]:
    """从 batch.yaml 中加载 data_files 段。

    新版布局中数据文件路径存放在 batch.yaml 而非 config.yaml 的 paths 中。

    Returns
    -------
    dict: data_files 映射（如 flux_file, mid_file），文件不存在或无法加载时返回空 dict。
    """
    if not batch_file_ref:
        return {}

    batch_path = (cfg_path.parent / batch_file_ref).resolve()
    if not batch_path.exists():
        return {}

    try:
        batch_doc = _load_yaml(batch_path)
        data_files = batch_doc.get("data_files", {}) or {}
        return {str(k): str(v) for k, v in data_files.items() if v}
    except Exception:
        return {}


def _extract_colors_from_conditions(conditions: dict[str, dict]) -> dict[str, str] | None:
    """从 conditions 的 color 字段提取颜色映射。key 统一小写以匹配 CSV 数据。"""
    colors = {}
    for key, val in conditions.items():
        if isinstance(val, dict) and val.get("color"):
            colors[key.lower()] = val["color"]
    return colors if colors else None


def _parse_comparisons(
    comparisons_raw: list[Any],
    conditions: dict[str, Any],
) -> list[tuple[str, str]]:
    """
    解析比较列表。

    支持两种写法：
      - ["treat", "control"]
      - {treat: "...", control: "..."} / {treat_group: "...", control_group: "..."}
    """
    valid_conditions = {str(k).lower() for k in conditions.keys()}
    comparisons: list[tuple[str, str]] = []

    for item in comparisons_raw:
        if isinstance(item, dict):
            treat = (
                item.get("treat")
                or item.get("treat_group")
                or item.get("treatment")
                or item.get("case")
            )
            control = (
                item.get("control")
                or item.get("control_group")
                or item.get("reference")
            )
        elif isinstance(item, (list, tuple)) and len(item) >= 2:
            treat, control = item[0], item[1]
        else:
            raise ValueError(f"Invalid comparison item: {item!r}")

        treat = str(treat).lower()
        control = str(control).lower()
        if treat not in valid_conditions:
            raise ValueError(f"Comparison treat group not in conditions: {treat}")
        if control not in valid_conditions:
            raise ValueError(f"Comparison control group not in conditions: {control}")
        comparisons.append((treat, control))

    return comparisons


def _validate_analysis_config(analysis: dict[str, Any]) -> None:
    """校验受控枚举，避免配置拼写错误静默改变分析含义。"""
    preprocessing = analysis.get("preprocessing", {})
    diff = analysis.get("differential_flux", {})
    heatmap = analysis.get("heatmap", {})

    statistic = str(diff.get("statistic", "auto")).lower()
    allowed_stats = {"auto", "mannwhitney", "welch_ttest", "effect_size_only"}
    if statistic not in allowed_stats:
        raise ValueError(f"analysis.differential_flux.statistic must be one of {sorted(allowed_stats)}")

    mt = str(diff.get("multiple_testing", "fdr_bh")).lower()
    if mt not in {"fdr_bh", "none"}:
        raise ValueError("analysis.differential_flux.multiple_testing must be 'fdr_bh' or 'none'")

    scale_method = str(preprocessing.get("scale_method", "zscore")).lower()
    if scale_method not in {"none", "zscore"}:
        raise ValueError("analysis.preprocessing.scale_method must be 'none' or 'zscore'")

    missing_value = str(preprocessing.get("missing_value", "drop")).lower()
    if missing_value not in {"drop", "median_fill"}:
        raise ValueError("analysis.preprocessing.missing_value must be 'drop' or 'median_fill'")

    rank_by = str(heatmap.get("rank_by", "auto")).lower()
    if rank_by not in {"auto", "abs_log2fc", "abs_cohen_d", "pvalue"}:
        raise ValueError("analysis.heatmap.rank_by must be auto|abs_log2fc|abs_cohen_d|pvalue")

    ident = analysis.get("identifiability", {}) or {}
    bound_values = ident.get("bound_values", [0.0, 1.0])
    if not isinstance(bound_values, (list, tuple)):
        raise ValueError("analysis.identifiability.bound_values must be a list of numbers")
    try:
        [float(v) for v in bound_values]
    except (TypeError, ValueError):
        raise ValueError("analysis.identifiability.bound_values must be a list of numbers")
    if float(ident.get("bound_tolerance", 0.05)) < 0:
        raise ValueError("analysis.identifiability.bound_tolerance must be >= 0")
    min_fraction = float(ident.get("bound_pinned_min_fraction", 0.5))
    if not 0 < min_fraction <= 1:
        raise ValueError("analysis.identifiability.bound_pinned_min_fraction must be in (0, 1]")


def generate_comparisons(conditions: dict[str, Any]) -> list[tuple[str, str]]:
    """
    自动从 conditions 生成所有两两比较组合。
    返回 [(a, b), ...] 列表，按 key 排序。
    """
    from itertools import combinations
    keys = sorted(conditions.keys())
    return list(combinations(keys, 2))


def get_group_colors(conditions: dict[str, Any]) -> dict[str, str]:
    """从 conditions 提取颜色，缺失时使用默认调色板。key 统一小写。"""
    _DEFAULT_PALETTE = [
        "#4C72B0", "#DD8452", "#55A868", "#C44E52",
        "#8172B3", "#937860", "#DA8BC3", "#8C8C8C",
    ]
    colors = {}
    for i, (key, val) in enumerate(sorted(conditions.items())):
        k = key.lower()
        if isinstance(val, dict) and val.get("color"):
            colors[k] = val["color"]
        else:
            colors[k] = _DEFAULT_PALETTE[i % len(_DEFAULT_PALETTE)]
    return colors


# ══════════════════════════════════════════════════════════════════════════════
# 多批次 / 多协议整合：协议发现与配置加载
# ══════════════════════════════════════════════════════════════════════════════

def discover_experiment_protocols(
    experiment_id: str,
) -> list[dict[str, Any]]:
    """扫描 experiment.yaml，返回所有已启用协议及其批次的发现结果。

    Returns
    -------
    list of dict，每个 dict 包含：
        protocol_id : str
        enabled : bool
        batches : list[dict]，每个包含 batch_id, config_path (Path | None), has_data (bool)
    """
    mfa_dir = repo_root() / "experiments" / experiment_id / "mfa_data"

    # Legacy 模式：无 experiment.yaml，从扁平 config.yaml 推断
    experiment_cfg = mfa_dir / "experiment.yaml"
    if not experiment_cfg.exists():
        legacy_cfg = mfa_dir / "config.yaml"
        return [{
            "protocol_id": None,
            "enabled": True,
            "batches": [{
                "batch_id": None,
                "config_path": legacy_cfg if legacy_cfg.exists() else None,
                "has_data": legacy_cfg.exists(),
            }],
        }]

    exp_doc = _load_yaml(experiment_cfg)
    tracer_protocols = exp_doc.get("tracer_protocols", {}) or {}
    result = []

    for protocol_id, spec in tracer_protocols.items():
        is_enabled = bool(spec.get("enabled", False))
        batch_names = spec.get("batches", []) or []
        batches = []

        for batch_id in batch_names:
            config_path = mfa_dir / protocol_id / batch_id / "config.yaml"
            has_data = config_path.exists()
            # 进一步检查 flux_file 是否实际存在
            if has_data:
                try:
                    batch_cfg = _load_yaml(config_path)
                    batch_file_ref = (batch_cfg.get("batch_file")
                                      or batch_cfg.get("paths", {}).get("batch_file"))
                    if batch_file_ref:
                        batch_meta_path = (config_path.parent / batch_file_ref).resolve()
                        if batch_meta_path.exists():
                            batch_meta = _load_yaml(batch_meta_path)
                            data_files = batch_meta.get("data_files", {}) or {}
                            flux_ref = data_files.get("flux_file", "")
                            if flux_ref:
                                flux_path = (config_path.parent / flux_ref).resolve()
                                has_data = flux_path.exists()
                            else:
                                has_data = False
                        else:
                            has_data = False
                except Exception:
                    has_data = False

            batches.append({
                "batch_id": batch_id,
                "config_path": config_path if config_path.exists() else None,
                "has_data": has_data,
            })

        result.append({
            "protocol_id": protocol_id,
            "enabled": is_enabled,
            "batches": batches,
        })

    _warn_undeclared_directories(mfa_dir, tracer_protocols)
    return result


def _warn_undeclared_directories(
    mfa_dir: Path,
    tracer_protocols: dict[str, Any],
) -> None:
    """提示磁盘上存在、但 experiment.yaml 未声明的 protocol / batch 目录。

    发现流程完全由 ``tracer_protocols`` 驱动，未声明的目录不会报错、也不会
    运行——新增一个批次却忘了登记，结果就是"什么都没发生"。这类静默遗漏
    比报错更难察觉，因此显式提示。
    """
    declared_protocols = {str(k) for k in tracer_protocols}
    on_disk_protocols = {
        d.name for d in mfa_dir.iterdir()
        if d.is_dir() and not d.name.startswith((".", "_"))
    }
    undeclared_protocols = sorted(on_disk_protocols - declared_protocols)
    if undeclared_protocols:
        print(f"[WARN] experiment.yaml 未声明的 protocol 目录（不会被分析）: "
              f"{undeclared_protocols}")

    for protocol_id, spec in tracer_protocols.items():
        protocol_dir = mfa_dir / str(protocol_id)
        if not protocol_dir.is_dir():
            continue
        declared_batches = {str(b) for b in (spec.get("batches", []) or [])}
        on_disk_batches = {
            d.name for d in protocol_dir.iterdir()
            if d.is_dir() and not d.name.startswith((".", "_"))
        }
        undeclared = sorted(on_disk_batches - declared_batches)
        if undeclared:
            print(f"[WARN] {protocol_id}: experiment.yaml 未声明的批次目录"
                  f"（不会被分析）: {undeclared}")
            print(f"       如需纳入，请加入 tracer_protocols.{protocol_id}.batches")


def load_integration_config(experiment_id: str) -> dict[str, Any]:
    """加载 experiment.yaml 中的 integration 配置段。

    Legacy 实验（无 experiment.yaml）返回 ``{"enabled": False}``。
    """
    mfa_dir = repo_root() / "experiments" / experiment_id / "mfa_data"
    experiment_cfg = mfa_dir / "experiment.yaml"

    if not experiment_cfg.exists():
        return {"enabled": False}

    exp_doc = _load_yaml(experiment_cfg)
    integration = exp_doc.get("integration", {}) or {}
    # 确保顶层 enabled 字段存在
    integration.setdefault("enabled", False)
    return integration


def resolve_output_root(experiment_id: str) -> Path:
    """解析实验的 mfa_results 输出根目录。

    从 experiment.yaml 的 defaults.output_root 读取，支持 ${experiment_id} 变量替换。
    回退路径：analysis_results/<experiment_id>/mfa_results
    """
    mfa_dir = repo_root() / "experiments" / experiment_id / "mfa_data"
    experiment_cfg = mfa_dir / "experiment.yaml"

    if experiment_cfg.exists():
        exp_doc = _load_yaml(experiment_cfg)
        defaults = exp_doc.get("defaults", {}) or {}
        raw = defaults.get("output_root", "")
        if raw:
            resolved = _substitute_vars(str(raw), {"experiment_id": experiment_id})
            return _resolve_path(resolved, mfa_dir)

    return repo_root() / "analysis_results" / experiment_id / "mfa_results"
