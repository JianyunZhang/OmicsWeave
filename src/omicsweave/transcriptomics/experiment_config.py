from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from omicsweave.core.project import project_dir


@dataclass(frozen=True)
class TranscriptomicsConfig:
    experiment_id: str
    config_path: Path
    output_root: Path
    output_dir: Path
    protocol_id: str | None
    batch_id: str | None
    experiment_path: Path | None
    protocol_path: Path | None
    batch_path: Path | None
    raw_data_root: str
    sample_list_tsv: str
    kallisto_exe: Path
    transcriptome_fasta: Path
    kallisto_index: Path
    threads: int
    qc_sample_size: int | None
    tissues: list[str]
    groups_path: Path
    conditions: dict[str, dict[str, Any]]
    mouse_groups: dict[str, str]
    sample_batches: dict[str, str]
    comparisons: list[tuple[str, str, str]] | None = None
    tissue_labels: dict[str, str] | None = None
    plot_colors: dict[str, Any] | None = None
    integration: dict[str, Any] | None = None
    # QC 阈值（experiment.yaml analysis_defaults 与 batch.yaml batch_qc 合并结果）
    qc_thresholds: dict[str, Any] | None = None
    on_qc_fail: str = "flag"
    # 技术参数（batch.yaml technical_params）——跨批次一致性校验依据
    technical_params: dict[str, Any] | None = None
    strandedness: str = "unstranded"
    gene_aggregation: str = "lengthScaledTPM"


def repo_root() -> Path:
    """项目目录（实验配置与结果所在仓库），见 :func:`omicsweave.core.project.project_dir`。"""
    return project_dir()


def resolve_qc_thresholds(
    experiment_doc: dict[str, Any] | None,
    batch_doc: dict[str, Any] | None,
) -> tuple[dict[str, Any], str]:
    """
    合并 QC 阈值：experiment.yaml 提供默认值，batch.yaml 可覆盖。

    读取 ``analysis_defaults.qc``（碱基质量）、``analysis_defaults.alignment``
    （比对率）以及 ``batch_qc``（批次级覆盖与失败策略）。

    这些字段此前全是「声明但不执行」——本函数是把它们接入流程的入口。

    Returns
    -------
    (thresholds, on_qc_fail)
    """
    thresholds: dict[str, Any] = {}

    defaults = ((experiment_doc or {}).get("analysis_defaults", {}) or {})
    for section in ("qc", "alignment"):
        for key, value in (defaults.get(section, {}) or {}).items():
            thresholds[str(key)] = value

    batch_qc = ((batch_doc or {}).get("batch_qc", {}) or {})
    on_qc_fail = str(batch_qc.get("on_qc_fail", "flag")).strip().lower() or "flag"

    if not bool(batch_qc.get("enabled", True)):
        # 批次显式关闭 QC 门控：返回空阈值，等效于不做判定
        return {}, "flag"

    for key, value in (batch_qc.get("thresholds", {}) or {}).items():
        thresholds[str(key)] = value

    return thresholds, on_qc_fail


def resolve_technical_params(
    protocol_doc: dict[str, Any] | None,
    batch_doc: dict[str, Any] | None,
) -> dict[str, Any]:
    """
    合并技术参数：protocol.yaml 给出协议级默认，batch.yaml 覆盖批次实测值。

    这些参数（建库试剂盒、链特异性、读长、RIN 等）决定了两个批次的计数
    是否可比，是跨批次整合前一致性断言的判定依据。
    """
    merged: dict[str, Any] = {}
    for doc in (protocol_doc, batch_doc):
        for key, value in ((doc or {}).get("technical_params", {}) or {}).items():
            merged[str(key)] = value
    return merged


def _is_abs_path_like(p: str) -> bool:
    if not p:
        return False
    # UNC path or Windows drive
    return p.startswith("\\\\") or bool(re.match(r"^[A-Za-z]:[\\/]", p))


_REPO_LEVEL_PREFIXES = ("results/", "analysis_results/", "analysis_tools/", "experiments/")


def _resolve_path(p: str, base_dir: Path) -> Path:
    if _is_abs_path_like(p):
        return Path(p)
    # Normalize to forward slashes for prefix check
    p_norm = p.replace("\\", "/")
    # Paths starting with known repo-level prefixes are always repo-relative,
    # even if the target does not exist yet (e.g. first-run output directories).
    for prefix in _REPO_LEVEL_PREFIXES:
        if p_norm.startswith(prefix):
            return (repo_root() / p).resolve()
    # For other relative paths, try repo-relative first (if exists), then config-relative
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


def _infer_experiment_id_from_path(cfg_path: Path) -> str | None:
    """Infer experiment id from .../experiments/<id>/transcriptomics_data/..."""
    parts = list(cfg_path.parts)
    for i, part in enumerate(parts):
        if part.lower() == "experiments" and i + 1 < len(parts):
            return parts[i + 1]
    for i, part in enumerate(parts):
        if part.lower() in {"transcriptomics_data", "transcriptomics"} and i - 1 >= 0:
            return parts[i - 1]
    return None


def _find_transcriptomics_root(cfg_path: Path, exp_id: str | None = None) -> Path:
    """Return transcriptomics_data root directory for a config path."""
    parts = list(cfg_path.parts)
    for i, part in enumerate(parts):
        if part.lower() in {"transcriptomics_data", "transcriptomics"}:
            return Path(*parts[: i + 1])
    if exp_id:
        root = repo_root() / "experiments" / exp_id / "transcriptomics_data"
        if root.exists():
            return root
        return repo_root() / "experiments" / exp_id / "transcriptomics"
    return cfg_path.parent


def _split_output_root(path: Path) -> tuple[Path, bool]:
    """
    Return (output_root, is_run_path).

    is_run_path=True means `path` already contains `/runs/` segment.
    """
    parts = list(path.parts)
    for i, part in enumerate(parts):
        if part.lower() == "runs":
            if i == 0:
                return path, True
            return Path(*parts[:i]), True
    return path, False


def _discover_active_batch_from_experiment(
    transcriptomics_root: Path,
    exp_doc: dict[str, Any],
) -> Path | None:
    defaults = exp_doc.get("defaults", {}) or {}
    active_protocol = str(defaults.get("active_protocol", "")).strip()
    active_batch = str(defaults.get("active_batch", "")).strip()
    if active_protocol and active_batch:
        candidate = transcriptomics_root / active_protocol / active_batch / "config.yaml"
        if candidate.exists():
            return candidate
        candidate = transcriptomics_root / active_protocol / active_batch / "batch.yaml"
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
                candidate = transcriptomics_root / str(protocol_id) / str(batch_id) / "config.yaml"
                if candidate.exists():
                    return candidate
                candidate = transcriptomics_root / str(protocol_id) / str(batch_id) / "batch.yaml"
                if candidate.exists():
                    return candidate

    for candidate in sorted(transcriptomics_root.glob("*/*/config.yaml")):
        return candidate
    for candidate in sorted(transcriptomics_root.glob("*/*/batch.yaml")):
        return candidate
    return None


def default_config_path(experiment_id: str) -> Path:
    """
    Resolve default transcriptomics config path.

    Preference:
    1. experiment.yaml -> active batch config.yaml/batch.yaml
    2. Legacy config.yaml in transcriptomics_data/ or transcriptomics/
    3. fallback legacy path in transcriptomics_data/config.yaml
    """
    exp_root = repo_root() / "experiments" / experiment_id
    roots = [exp_root / "transcriptomics_data", exp_root / "transcriptomics"]

    for root in roots:
        experiment_yaml = root / "experiment.yaml"
        if not experiment_yaml.exists():
            continue
        exp_doc = _load_yaml(experiment_yaml)
        active = _discover_active_batch_from_experiment(root, exp_doc)
        if active:
            return active

    for root in roots:
        legacy = root / "config.yaml"
        if legacy.exists():
            return legacy

    return exp_root / "transcriptomics_data" / "config.yaml"


def _infer_protocol_batch(
    cfg_path: Path,
    batch_doc: dict[str, Any] | None = None,
) -> tuple[str | None, str | None]:
    protocol_id = None
    batch_id = None

    if batch_doc:
        protocol_id = str(batch_doc.get("protocol_id", "")).strip() or None
        batch_id = str(batch_doc.get("batch_id", "")).strip() or None
        if protocol_id and batch_id:
            return protocol_id, batch_id

    parts = list(cfg_path.parts)
    idx = None
    for i, part in enumerate(parts):
        if part.lower() in {"transcriptomics_data", "transcriptomics"}:
            idx = i
            break

    if idx is not None and len(parts) >= idx + 4:
        protocol_id = protocol_id or parts[idx + 1]
        batch_id = batch_id or parts[idx + 2]

    return protocol_id, batch_id


def _resolve_output_dir_from_experiment(
    transcriptomics_root: Path,
    experiment_id: str,
    protocol_id: str | None,
    batch_id: str | None,
) -> str:
    experiment_yaml = transcriptomics_root / "experiment.yaml"
    if not experiment_yaml.exists():
        return ""
    exp_doc = _load_yaml(experiment_yaml)
    defaults = exp_doc.get("defaults", {}) or {}
    raw = str(defaults.get("output_root", "")).strip()
    if not raw:
        return ""

    output_root = _resolve_path(
        _substitute_vars(raw, {"experiment_id": experiment_id}),
        transcriptomics_root,
    )
    if protocol_id and batch_id:
        return str(output_root / "runs" / protocol_id / batch_id)
    return str(output_root)


def _parse_comparisons(
    comparisons_raw: Any,
    conditions: dict[str, Any],
    strict: bool = False,
) -> tuple[list[tuple[str, str, str]] | None, list[str]]:
    """
    Parse comparison entries from config into (numerator, denominator, comp_name) tuples.

    Parameters
    ----------
    comparisons_raw : Any
        Raw comparisons list from config YAML.
    conditions : dict
        Conditions dict from groups.yaml, used to validate group names.
    strict : bool
        If True, raise ValueError on unrecognized condition names.
        If False (default), log warnings and skip invalid entries.

    Returns
    -------
    tuple[list[tuple[str, str, str]] | None, list[str]]
        (parsed comparisons, list of warning messages)
    """
    warnings: list[str] = []

    if not comparisons_raw:
        return None, warnings
    if not isinstance(comparisons_raw, list):
        warnings.append(f"comparisons is not a list (got {type(comparisons_raw).__name__}), ignoring")
        return None, warnings

    lookup = {str(k).lower(): str(k) for k in conditions.keys()}
    parsed: list[tuple[str, str, str]] = []
    seen: set[str] = set()

    for idx, item in enumerate(comparisons_raw):
        if isinstance(item, dict):
            num = (
                item.get("treat")
                or item.get("treatment")
                or item.get("numerator")
                or item.get("case")
                or item.get("group_a")
            )
            den = (
                item.get("control")
                or item.get("reference")
                or item.get("denominator")
                or item.get("group_b")
            )
        elif isinstance(item, (list, tuple)) and len(item) >= 2:
            num, den = item[0], item[1]
        else:
            warnings.append(f"comparisons[{idx}]: unrecognized format ({type(item).__name__}), skipped")
            continue

        if num is None or den is None:
            warnings.append(f"comparisons[{idx}]: missing numerator or denominator, skipped")
            continue

        num_str = str(num).strip()
        den_str = str(den).strip()
        num_key = lookup.get(num_str.lower())
        den_key = lookup.get(den_str.lower())

        if not num_key:
            msg = f"comparisons[{idx}]: numerator '{num_str}' not found in conditions {list(conditions.keys())}"
            if strict:
                raise ValueError(msg)
            warnings.append(msg)
            continue
        if not den_key:
            msg = f"comparisons[{idx}]: denominator '{den_str}' not found in conditions {list(conditions.keys())}"
            if strict:
                raise ValueError(msg)
            warnings.append(msg)
            continue
        if num_key == den_key:
            warnings.append(f"comparisons[{idx}]: numerator == denominator ('{num_key}'), skipped")
            continue

        comp_name = f"{num_key}_vs_{den_key}"
        if comp_name in seen:
            warnings.append(f"comparisons[{idx}]: duplicate comparison '{comp_name}', skipped")
            continue
        seen.add(comp_name)
        parsed.append((num_key, den_key, comp_name))

    return (parsed or None), warnings


def _load_groups(
    groups_path: Path,
) -> tuple[dict[str, dict[str, Any]], dict[str, str], dict[str, str]]:
    groups = _load_yaml(groups_path)
    conditions = groups.get("conditions", {}) or {}
    sample_groups = groups.get("sample_groups", {}) or groups.get("mouse_groups", {}) or {}
    sample_groups = {str(k): str(v) for k, v in sample_groups.items()}
    sample_batches = groups.get("sample_batches", {}) or groups.get("mouse_batches", {}) or {}
    sample_batches = {str(k): str(v) for k, v in sample_batches.items()}
    if not conditions:
        raise ValueError(f"conditions missing/empty in {groups_path}")
    if not sample_groups:
        raise ValueError(f"sample_groups/mouse_groups missing/empty in {groups_path}")
    return conditions, sample_groups, sample_batches


def load_transcriptomics_config(
    *,
    experiment_id: str | None = None,
    config_path: str | None = None,
    outdir_override: str | None = None,
) -> TranscriptomicsConfig:
    exp = experiment_id or os.environ.get("EXPERIMENT_ID")
    if config_path:
        cfg_path = Path(config_path).resolve()
    else:
        if not exp:
            raise ValueError("Missing experiment id/config. Please pass --experiment <id> or --config <path>.")
        cfg_path = default_config_path(exp)

    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")

    cfg = _load_yaml(cfg_path)
    base_dir = cfg_path.parent

    if not exp:
        exp = str(cfg.get("experiment_id", "")).strip() or None
    exp = exp or str(cfg.get("experiment_id", "")).strip() or _infer_experiment_id_from_path(cfg_path)
    if not exp:
        raise ValueError(f"Cannot infer experiment_id from config path: {cfg_path}")

    transcriptomics_root = _find_transcriptomics_root(cfg_path, exp)

    # ------------------------------------------------------------------
    # Layout A: Legacy full config (paths/reference/params in config.yaml)
    # Layout B: Batch run-entry config (config.yaml with batch_file/groups_file)
    # Layout C: experiment.yaml / protocol.yaml / batch.yaml entrypoint
    # ------------------------------------------------------------------
    experiment_path: Path | None = None
    protocol_path: Path | None = None
    batch_path: Path | None = None
    experiment_doc: dict[str, Any] = {}
    protocol_doc: dict[str, Any] = {}
    batch_doc: dict[str, Any] = {}
    entry_doc: dict[str, Any] = {}

    cfg_name = cfg_path.name.lower()
    has_legacy_sections = any(k in cfg for k in ("paths", "reference", "params"))
    is_entry_config = "batch_file" in cfg or "groups_file" in cfg or "comparisons" in cfg

    if cfg_name == "experiment.yaml":
        experiment_path = cfg_path
        experiment_doc = cfg
        active = _discover_active_batch_from_experiment(transcriptomics_root, experiment_doc)
        if not active:
            raise ValueError(f"No active batch config found from experiment.yaml: {cfg_path}")
        cfg_path = active
        cfg_name = cfg_path.name.lower()
        cfg = _load_yaml(cfg_path)
        base_dir = cfg_path.parent
        has_legacy_sections = any(k in cfg for k in ("paths", "reference", "params"))
        is_entry_config = "batch_file" in cfg or "groups_file" in cfg or "comparisons" in cfg

    if cfg_name == "protocol.yaml":
        protocol_path = cfg_path
        protocol_doc = cfg
        candidate = _discover_active_batch_from_experiment(
            transcriptomics_root,
            _load_yaml(transcriptomics_root / "experiment.yaml")
            if (transcriptomics_root / "experiment.yaml").exists()
            else {},
        )
        if candidate and candidate.parent.parent == cfg_path.parent:
            cfg_path = candidate
        else:
            batch_candidates = sorted(cfg_path.parent.glob("*/config.yaml"))
            if batch_candidates:
                cfg_path = batch_candidates[0]
            else:
                batch_yaml_candidates = sorted(cfg_path.parent.glob("*/batch.yaml"))
                if not batch_yaml_candidates:
                    raise ValueError(f"No batch config found under protocol: {cfg_path}")
                cfg_path = batch_yaml_candidates[0]
        cfg_name = cfg_path.name.lower()
        cfg = _load_yaml(cfg_path)
        base_dir = cfg_path.parent
        has_legacy_sections = any(k in cfg for k in ("paths", "reference", "params"))
        is_entry_config = "batch_file" in cfg or "groups_file" in cfg or "comparisons" in cfg

    if has_legacy_sections and not is_entry_config:
        # Legacy full config style
        paths = cfg.get("paths", {}) or {}
        ref = cfg.get("reference", {}) or {}
        params = cfg.get("params", {}) or {}

        raw_data_root = str(paths.get("raw_data_root", "")).strip()
        sample_list_tsv = str(paths.get("sample_list_tsv", "")).strip()
        output_dir_s = str(paths.get("output_dir", "")).strip()
        if outdir_override:
            output_dir_s = outdir_override

        if not raw_data_root:
            raise ValueError(f"paths.raw_data_root missing in {cfg_path}")
        if not sample_list_tsv:
            raise ValueError(f"paths.sample_list_tsv missing in {cfg_path}")
        if not output_dir_s:
            raise ValueError(f"paths.output_dir missing in {cfg_path}")

        output_dir_s = _substitute_vars(output_dir_s, {"experiment_id": exp})
        output_dir_base = _resolve_path(output_dir_s, base_dir)
        output_root, is_run_path = _split_output_root(output_dir_base)
        output_dir = output_dir_base

        vars_map = {
            "experiment_id": exp,
            "output_root": str(output_root),
            "output_dir": str(output_dir),
            "protocol_id": "",
            "batch_id": "",
        }

        kallisto_exe_s = _substitute_vars(str(ref.get("kallisto_exe", "")).strip(), vars_map)
        transcriptome_fasta_s = _substitute_vars(str(ref.get("transcriptome_fasta", "")).strip(), vars_map)
        kallisto_index_s = _substitute_vars(str(ref.get("kallisto_index", "")).strip(), vars_map)
        if not kallisto_exe_s:
            raise ValueError(f"reference.kallisto_exe missing in {cfg_path}")
        if not transcriptome_fasta_s:
            raise ValueError(f"reference.transcriptome_fasta missing in {cfg_path}")
        if not kallisto_index_s:
            raise ValueError(f"reference.kallisto_index missing in {cfg_path}")

        kallisto_exe = _resolve_path(kallisto_exe_s, base_dir)
        transcriptome_fasta = _resolve_path(transcriptome_fasta_s, base_dir)
        kallisto_index = _resolve_path(kallisto_index_s, base_dir)

        threads = int(params.get("threads", 8))
        qc_sample_size = params.get("qc_sample_size", None)
        qc_sample_size = int(qc_sample_size) if qc_sample_size is not None else None
        tissues = list(params.get("tissues", []))
        if not tissues:
            raise ValueError(f"params.tissues missing/empty in {cfg_path}")

        groups_file = str(cfg.get("groups_file", "groups.yaml")).strip() or "groups.yaml"
        groups_path = _resolve_path(groups_file, base_dir)
        conditions, sample_groups, sample_batches = _load_groups(groups_path)
        comparisons, comp_warnings = _parse_comparisons(cfg.get("comparisons", None), conditions)
        for w in comp_warnings:
            print(f"  [WARN] {w}")

        tissue_labels = params.get("tissue_labels", None)
        if tissue_labels:
            tissue_labels = {str(k): str(v) for k, v in tissue_labels.items()}
        plot_colors = params.get("plot_colors", None)

        return TranscriptomicsConfig(
            experiment_id=exp,
            config_path=cfg_path,
            output_root=output_root if not is_run_path else output_root,
            output_dir=output_dir,
            protocol_id=None,
            batch_id=None,
            experiment_path=None,
            protocol_path=None,
            batch_path=None,
            raw_data_root=_substitute_vars(raw_data_root, vars_map),
            sample_list_tsv=_substitute_vars(sample_list_tsv, vars_map),
            kallisto_exe=kallisto_exe,
            transcriptome_fasta=transcriptome_fasta,
            kallisto_index=kallisto_index,
            threads=threads,
            qc_sample_size=qc_sample_size,
            tissues=tissues,
            groups_path=groups_path,
            conditions=conditions,
            mouse_groups=sample_groups,
            sample_batches=sample_batches,
            comparisons=comparisons,
            tissue_labels=tissue_labels,
            plot_colors=plot_colors,
            integration=None,
            qc_thresholds=resolve_qc_thresholds(cfg, cfg)[0],
            on_qc_fail=resolve_qc_thresholds(cfg, cfg)[1],
            technical_params=resolve_technical_params(None, cfg),
            strandedness=str(params.get("strandedness", "unstranded")),
            gene_aggregation=str(params.get("gene_aggregation", "lengthScaledTPM")),
        )

    # New layout: entry config + batch/protocol/experiment hierarchy
    if cfg_name == "batch.yaml":
        batch_path = cfg_path
        batch_doc = cfg
        entry_doc = {}
    else:
        entry_doc = cfg
        batch_file = str(entry_doc.get("batch_file", "batch.yaml")).strip() or "batch.yaml"
        batch_path = _resolve_path(batch_file, base_dir)
        if not batch_path.exists():
            raise FileNotFoundError(f"batch_file not found: {batch_path}")
        batch_doc = _load_yaml(batch_path)

    if not protocol_path:
        candidate = batch_path.parent.parent / "protocol.yaml"
        if candidate.exists():
            protocol_path = candidate
            protocol_doc = _load_yaml(candidate)

    if not experiment_path:
        candidate = transcriptomics_root / "experiment.yaml"
        if candidate.exists():
            experiment_path = candidate
            experiment_doc = _load_yaml(candidate)

    protocol_id, batch_id = _infer_protocol_batch(batch_path, batch_doc)

    paths = batch_doc.get("paths", {}) or {}
    ref = batch_doc.get("reference", {}) or {}
    params = batch_doc.get("params", {}) or {}

    raw_data_root = str(paths.get("raw_data_root", "")).strip()
    sample_list_tsv = str(paths.get("sample_list_tsv", "")).strip()
    output_dir_s = str(paths.get("output_dir", "")).strip()

    if outdir_override:
        output_dir_s = outdir_override
    elif not output_dir_s and experiment_path:
        output_dir_s = _resolve_output_dir_from_experiment(transcriptomics_root, exp, protocol_id, batch_id)

    if not raw_data_root:
        raise ValueError(f"paths.raw_data_root missing in {batch_path}")
    if not sample_list_tsv:
        raise ValueError(f"paths.sample_list_tsv missing in {batch_path}")
    if not output_dir_s:
        raise ValueError(f"paths.output_dir missing in {batch_path}")

    output_dir_s = _substitute_vars(
        output_dir_s,
        {
            "experiment_id": exp,
            "protocol_id": protocol_id or "",
            "batch_id": batch_id or "",
        },
    )
    output_dir_base = _resolve_path(output_dir_s, batch_path.parent)
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

    kallisto_exe_s = _substitute_vars(str(ref.get("kallisto_exe", "")).strip(), vars_map)
    transcriptome_fasta_s = _substitute_vars(str(ref.get("transcriptome_fasta", "")).strip(), vars_map)
    kallisto_index_s = _substitute_vars(str(ref.get("kallisto_index", "")).strip(), vars_map)
    if not kallisto_exe_s:
        raise ValueError(f"reference.kallisto_exe missing in {batch_path}")
    if not transcriptome_fasta_s:
        raise ValueError(f"reference.transcriptome_fasta missing in {batch_path}")
    if not kallisto_index_s:
        raise ValueError(f"reference.kallisto_index missing in {batch_path}")

    kallisto_exe = _resolve_path(kallisto_exe_s, batch_path.parent)
    transcriptome_fasta = _resolve_path(transcriptome_fasta_s, batch_path.parent)
    kallisto_index = _resolve_path(kallisto_index_s, batch_path.parent)

    threads = int(params.get("threads", 8))
    qc_sample_size = params.get("qc_sample_size", None)
    qc_sample_size = int(qc_sample_size) if qc_sample_size is not None else None
    tissues = list(params.get("tissues", []))
    if not tissues:
        raise ValueError(f"params.tissues missing/empty in {batch_path}")

    groups_file = str(entry_doc.get("groups_file", "groups.yaml")).strip() or "groups.yaml"
    groups_path = _resolve_path(_substitute_vars(groups_file, vars_map), batch_path.parent)
    conditions, sample_groups, sample_batches = _load_groups(groups_path)
    comparisons, comp_warnings = _parse_comparisons(entry_doc.get("comparisons", None), conditions)
    for w in comp_warnings:
        print(f"  [WARN] {w}")

    tissue_labels = params.get("tissue_labels", None)
    if tissue_labels:
        tissue_labels = {str(k): str(v) for k, v in tissue_labels.items()}
    plot_colors = params.get("plot_colors", None)

    integration = {
        "experiment": experiment_doc.get("integration", {}) if experiment_doc else {},
        "protocol": protocol_doc.get("integration_context", {}) if protocol_doc else {},
        "batch": batch_doc.get("integration", {}) if batch_doc else {},
    }

    qc_thresholds, on_qc_fail = resolve_qc_thresholds(experiment_doc, batch_doc)
    technical_params = resolve_technical_params(protocol_doc, batch_doc)
    # 链特异性/聚合方法优先级：batch params > technical_params > 全局默认
    strandedness = str(
        params.get("strandedness")
        or technical_params.get("strandedness")
        or "unstranded"
    )
    gene_aggregation = str(
        params.get("gene_aggregation")
        or ((experiment_doc or {}).get("analysis_defaults", {}) or {})
        .get("quantification", {})
        .get("gene_aggregation")
        or "lengthScaledTPM"
    )

    return TranscriptomicsConfig(
        experiment_id=exp,
        config_path=cfg_path,
        output_root=output_root,
        output_dir=output_dir,
        protocol_id=protocol_id,
        batch_id=batch_id,
        experiment_path=experiment_path,
        protocol_path=protocol_path,
        batch_path=batch_path,
        raw_data_root=_substitute_vars(raw_data_root, vars_map),
        sample_list_tsv=_substitute_vars(sample_list_tsv, vars_map),
        kallisto_exe=kallisto_exe,
        transcriptome_fasta=transcriptome_fasta,
        kallisto_index=kallisto_index,
        threads=threads,
        qc_sample_size=qc_sample_size,
        tissues=tissues,
        groups_path=groups_path,
        conditions=conditions,
        mouse_groups=sample_groups,
        sample_batches=sample_batches,
        comparisons=comparisons,
        tissue_labels=tissue_labels,
        plot_colors=plot_colors,
        integration=integration,
        qc_thresholds=qc_thresholds,
        on_qc_fail=on_qc_fail,
        technical_params=technical_params,
        strandedness=strandedness,
        gene_aggregation=gene_aggregation,
    )


# ============================================================================
# Helper functions for config-driven visualization / analysis
# ============================================================================

# Default color palette (colorblind-friendly, up to 8 groups)
_DEFAULT_PALETTE = [
    "#4C72B0",  # blue
    "#DD8452",  # orange
    "#55A868",  # green
    "#C44E52",  # red
    "#8172B3",  # purple
    "#937860",  # brown
    "#DA8BC3",  # pink
    "#8C8C8C",  # gray
]


def generate_comparisons(conditions: dict[str, Any]) -> list[tuple[str, str, str]]:
    """
    Auto-generate all pairwise comparisons from conditions.
    Returns list of (numerator, denominator, comp_name) tuples.
    comp_name format: '{a}_vs_{b}'.
    """
    keys = sorted(conditions.keys())
    comps = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            comp_name = f"{a}_vs_{b}"
            comps.append((a, b, comp_name))
    return comps


def resolve_comparisons(cfg: TranscriptomicsConfig) -> list[tuple[str, str, str]]:
    """
    Return comparison list from config if provided, otherwise auto-generate pairwise.
    """
    return list(cfg.comparisons) if cfg.comparisons else generate_comparisons(cfg.conditions)


def get_group_colors(conditions: dict[str, Any]) -> dict[str, str]:
    """
    Build group -> color mapping from conditions' `color` field.
    Falls back to _DEFAULT_PALETTE for groups without a configured color.
    """
    colors = {}
    for i, (key, val) in enumerate(sorted(conditions.items())):
        if isinstance(val, dict) and val.get("color"):
            colors[key] = val["color"]
        else:
            colors[key] = _DEFAULT_PALETTE[i % len(_DEFAULT_PALETTE)]
    return colors


def get_tissue_labels(cfg: TranscriptomicsConfig) -> dict[str, str]:
    """
    Return tissue display labels.
    Priority: config params.tissue_labels > capitalize tissue name.
    """
    if cfg.tissue_labels:
        return dict(cfg.tissue_labels)
    return {t: t.replace("_", " ").title() for t in cfg.tissues}


def get_comp_labels(
    conditions: dict[str, Any],
    comparisons: list[tuple[str, str, str]],
) -> dict[str, str]:
    """
    Build comparison name -> display label mapping.
    Uses conditions' `display` field if available, else the condition code.
    Label format: '{display_a} vs {display_b}'.
    """

    def _disp(code: str) -> str:
        val = conditions.get(code)
        if isinstance(val, dict) and val.get("display"):
            return val["display"]
        return code

    return {comp_name: f"{_disp(a)} vs {_disp(b)}" for a, b, comp_name in comparisons}


def get_cond_map(conditions: dict[str, Any]) -> dict[str, str]:
    """
    Build condition code -> display name mapping.
    """
    result = {}
    for key, val in conditions.items():
        if isinstance(val, dict) and val.get("display"):
            result[key] = val["display"]
        else:
            result[key] = key
    return result


# Default color scheme for all plots — used when config has no plot_colors section
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
    "pca_tissue_cmap": "Set2",
}


def get_plot_colors(cfg: "TranscriptomicsConfig") -> dict[str, Any]:
    """
    Return plot color configuration with defaults.
    Merges user-provided plot_colors from config over built-in defaults.
    """
    raw = cfg.plot_colors or {}
    return {**_PLOT_COLOR_DEFAULTS, **raw}


def discover_experiment_protocols(experiment_id: str) -> list[dict[str, Any]]:
    """
    Discover all protocol/batch run-entry configs from experiment.yaml.

    Returns a list of:
      {
        "protocol_id": str | None,
        "enabled": bool,
        "batches": [
          {"batch_id": str | None, "config_path": Path | None, "has_data": bool}
        ]
      }
    """
    tx_dir = repo_root() / "experiments" / experiment_id / "transcriptomics_data"
    if not tx_dir.exists():
        tx_dir = repo_root() / "experiments" / experiment_id / "transcriptomics"

    experiment_yaml = tx_dir / "experiment.yaml"
    if not experiment_yaml.exists():
        legacy_cfg = tx_dir / "config.yaml"
        return [{
            "protocol_id": None,
            "enabled": True,
            "batches": [{
                "batch_id": None,
                "config_path": legacy_cfg if legacy_cfg.exists() else None,
                "has_data": legacy_cfg.exists(),
            }],
        }]

    exp_doc = _load_yaml(experiment_yaml)
    tracer_protocols = exp_doc.get("tracer_protocols", {}) or {}
    results: list[dict[str, Any]] = []

    for protocol_id, spec in tracer_protocols.items():
        if not isinstance(spec, dict):
            continue
        is_enabled = bool(spec.get("enabled", False))
        batch_ids = spec.get("batches", []) or []
        batches: list[dict[str, Any]] = []

        for batch_id in batch_ids:
            run_cfg = tx_dir / str(protocol_id) / str(batch_id) / "config.yaml"
            if not run_cfg.exists():
                run_cfg = tx_dir / str(protocol_id) / str(batch_id) / "batch.yaml"
            has_data = run_cfg.exists()
            batches.append({
                "batch_id": str(batch_id),
                "config_path": run_cfg if run_cfg.exists() else None,
                "has_data": has_data,
            })

        results.append({
            "protocol_id": str(protocol_id),
            "enabled": is_enabled,
            "batches": batches,
        })

    return results


def load_integration_config(experiment_id: str) -> dict[str, Any]:
    tx_dir = repo_root() / "experiments" / experiment_id / "transcriptomics_data"
    if not tx_dir.exists():
        tx_dir = repo_root() / "experiments" / experiment_id / "transcriptomics"

    experiment_yaml = tx_dir / "experiment.yaml"
    if not experiment_yaml.exists():
        return {"enabled": False}

    exp_doc = _load_yaml(experiment_yaml)
    integration = exp_doc.get("integration", {}) or {}
    integration.setdefault("enabled", False)
    return integration


def resolve_output_root(experiment_id: str) -> Path:
    tx_dir = repo_root() / "experiments" / experiment_id / "transcriptomics_data"
    if not tx_dir.exists():
        tx_dir = repo_root() / "experiments" / experiment_id / "transcriptomics"

    experiment_yaml = tx_dir / "experiment.yaml"
    if experiment_yaml.exists():
        exp_doc = _load_yaml(experiment_yaml)
        defaults = exp_doc.get("defaults", {}) or {}
        raw = str(defaults.get("output_root", "")).strip()
        if raw:
            resolved = _substitute_vars(raw, {"experiment_id": experiment_id})
            return _resolve_path(resolved, tx_dir)

    # fallback from legacy config
    legacy_cfg = tx_dir / "config.yaml"
    if legacy_cfg.exists():
        cfg_doc = _load_yaml(legacy_cfg)
        paths = cfg_doc.get("paths", {}) or {}
        out = str(paths.get("output_dir", "")).strip()
        if out:
            out_path = _resolve_path(
                _substitute_vars(out, {"experiment_id": experiment_id}),
                tx_dir,
            )
            root, _ = _split_output_root(out_path)
            return root

    return repo_root() / "analysis_results" / experiment_id / "transcriptomics_results"


def config_arg_help() -> str:
    return (
        "Use --experiment <id> to auto-pick active batch from "
        "experiments/<id>/transcriptomics_data/experiment.yaml, or pass --config "
        "to a run-entry config.yaml (or batch.yaml). Legacy config.yaml is still supported."
    )
