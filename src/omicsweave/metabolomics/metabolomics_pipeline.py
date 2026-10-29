"""配置驱动、可增量执行的代谢组学分析流程。"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any


import numpy as np
import pandas as pd
import yaml

from omicsweave.common.metabolite_id_mapping import annotate_metabolites
from omicsweave.metabolomics.data_preprocessing import (
    PreprocessingOptions,
    classify_data_columns,
    get_sample_groups,
    load_tissue_data,
    preprocess_metabolomics_data,
)
from omicsweave.metabolomics.differential_analysis import (
    DifferentialAnalysisOptions,
    analyze_differential_metabolites,
)
from omicsweave.metabolomics.enrichment._run import run_pathway_enrichment
from omicsweave.metabolomics.experiment_config import (
    MetabolomicsConfig,
    generate_comparisons,
    get_cond_map,
    get_group_colors,
    get_plot_colors,
    get_tissue_labels,
)
from omicsweave.metabolomics.gsea._run import run_gsea
from omicsweave.metabolomics.plotting.heatmap import plot_heatmap
from omicsweave.metabolomics.plotting.qc_plots import (
    plot_preprocessing_qc,
    run_pca,
    run_plsda,
)
from omicsweave.metabolomics.plotting.volcano import plot_volcano


PIPELINE_VERSION = 2
ALL_STEPS = {1, 2, 3, 4, 5, 6}
STEP_DEPENDENCIES = {
    1: set(),
    2: {1},
    3: {1, 2},
    4: {1},
    5: {1, 2, 3},
    6: {1, 2, 3},
}


@dataclass
class TissueContext:
    tissue: str
    tissue_label: str
    processed_df: pd.DataFrame
    sample_cols: list[str]
    qc_cols: list[str]
    groups: dict[str, str]
    qc_metrics: pd.DataFrame
    sample_totals: pd.DataFrame
    filtered_csv: Path


def _apply_config(
    experiment_id: str | None,
    config_path: str | None,
    outdir_override: str | None,
) -> MetabolomicsConfig:
    from omicsweave.metabolomics.experiment_config import load_metabolomics_config

    return load_metabolomics_config(
        experiment_id=experiment_id,
        config_path=config_path,
        outdir_override=outdir_override,
    )


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )


def _manifest_matches(path: Path, expected: dict[str, Any]) -> bool:
    if not path.exists():
        return False
    try:
        actual = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return False
    return all(actual.get(key) == value for key, value in expected.items())


def _comparison_names(comparisons: list[tuple[str, str, str]]) -> list[str]:
    return [comparison[2] for comparison in comparisons]


def _resolve_steps(steps: list[int] | None) -> tuple[set[int], set[int], set[int]]:
    selected = set(steps) if steps is not None else set(ALL_STEPS)
    unknown = selected - ALL_STEPS
    if unknown:
        raise ValueError(f"未知步骤编号: {sorted(unknown)}；可选 1-6")
    if not selected:
        raise ValueError("steps 不得为空")
    needed = set(selected)
    for step in list(selected):
        needed.update(STEP_DEPENDENCIES[step])
    forced = set(steps or [])
    return selected, needed, forced


def _load_comparisons(cfg: MetabolomicsConfig) -> list[tuple[str, str, str]]:
    groups_yaml = yaml.safe_load(cfg.groups_path.read_text(encoding="utf-8")) or {}
    if not isinstance(groups_yaml, dict):
        raise ValueError(f"groups.yaml 根节点必须为 mapping: {cfg.groups_path}")
    return generate_comparisons(cfg.conditions, groups_yaml.get("contrast_groups"))


def _preprocessing_signature(cfg: MetabolomicsConfig, data_path: Path) -> dict[str, Any]:
    stat = data_path.stat()
    return {
        "pipeline_version": PIPELINE_VERSION,
        "input_size": stat.st_size,
        "input_mtime_ns": stat.st_mtime_ns,
        "normalization_method": cfg.normalization_method,
        "detection_rate": cfg.detection_rate,
        "imputation_fraction": cfg.imputation_fraction,
        "qc_filter_enabled": cfg.qc_filter_enabled,
        "qc_rsd_cutoff": cfg.qc_rsd_cutoff,
        "qc_min_detection_rate": cfg.qc_min_detection_rate,
        "batch_correction": cfg.batch_correction,
    }


def _run_preprocessing(
    cfg: MetabolomicsConfig,
    tissue: str,
    tissue_label: str,
    data_path: Path,
    step_dir: Path,
    signature: dict[str, Any],
) -> TissueContext:
    print("  [Step 1] 数据预处理")
    raw_df = load_tissue_data(data_path)
    classified = classify_data_columns(raw_df.columns, set(cfg.mouse_groups))
    if classified.unmapped_biological:
        raise ValueError(
            f"{tissue} 存在形似生物样本但未配置 mouse ID 的列: "
            f"{classified.unmapped_biological}"
        )
    if not classified.biological:
        raise ValueError(
            f"{tissue} 未识别到生物样本列；支持 <mouse_id>_<tissue> 或旧版连字符，"
            "且 mouse_id 必须存在于 groups.yaml"
        )
    if classified.blank:
        print(f"    警告: 检测到 {len(classified.blank)} 个 blank 列；当前峰表流程不执行 blank subtraction，已从分析矩阵排除")
        raw_df = raw_df.drop(columns=classified.blank)

    groups = get_sample_groups(
        classified.biological,
        cfg.mouse_groups,
        cfg.conditions,
        strict=True,
    )
    preprocessing = preprocess_metabolomics_data(
        raw_df,
        classified.biological,
        groups,
        qc_cols=classified.qc,
        options=PreprocessingOptions(
            normalization_method=cfg.normalization_method,
            detection_rate=cfg.detection_rate,
            imputation_fraction=cfg.imputation_fraction,
            qc_filter_enabled=cfg.qc_filter_enabled,
            qc_rsd_cutoff=cfg.qc_rsd_cutoff,
            qc_min_detection_rate=cfg.qc_min_detection_rate,
        ),
    )

    step_dir.mkdir(parents=True, exist_ok=True)
    filtered_csv = step_dir / "filtered_data.csv"
    preprocessing.processed_df.to_csv(filtered_csv, index=False, encoding="utf-8-sig")

    qc_table = preprocessing.qc_metrics.copy()
    metadata_cols = [col for col in ("Index", "Compounds", "Formula") if col in raw_df.columns]
    if metadata_cols and not qc_table.empty:
        qc_table = raw_df.loc[qc_table.index, metadata_cols].join(qc_table)
    qc_table.to_csv(step_dir / "qc_feature_metrics.csv", index=False, encoding="utf-8-sig")
    preprocessing.sample_totals.reset_index().to_csv(
        step_dir / "sample_total_signal.csv", index=False, encoding="utf-8-sig",
    )
    preprocessing.sample_totals.reset_index()[
        ["sample", "sample_type", "normalization_factor"]
    ].to_csv(step_dir / "normalization_factors.csv", index=False, encoding="utf-8-sig")

    group_counts = {
        group: sum(value == group for value in groups.values())
        for group in sorted(set(groups.values()))
    }
    finite_rsd = pd.to_numeric(
        preprocessing.qc_metrics.get("QC_RSD", pd.Series(dtype=float)),
        errors="coerce",
    ).dropna()
    summary = {
        **signature,
        **preprocessing.summary,
        "tissue": tissue,
        "tissue_label": tissue_label,
        "biological_sample_count": len(preprocessing.sample_cols),
        "qc_sample_count": len(preprocessing.qc_cols),
        "group_counts": group_counts,
        "qc_median_rsd": float(finite_rsd.median()) if not finite_rsd.empty else None,
        "qc_features_below_rsd_cutoff": int((finite_rsd <= cfg.qc_rsd_cutoff).sum()),
    }
    _write_json(step_dir / "preprocessing_summary.json", summary)
    text_lines = [
        f"组织: {tissue_label}",
        f"原始代谢物数: {summary['original_feature_count']}",
        f"过滤后代谢物数: {summary['filtered_feature_count']}",
        f"检出率过滤移除: {summary['removed_by_detection']}",
        f"QC 过滤移除: {summary['removed_by_qc']}",
        f"插补单元格数: {summary['imputed_cell_count']}",
        f"归一化: {cfg.normalization_method}",
        f"批次校正: {cfg.batch_correction}",
        f"生物样本数: {len(preprocessing.sample_cols)}",
        f"QC 样本数: {len(preprocessing.qc_cols)}",
        f"组别: {group_counts}",
        f"QC 中位 RSD: {summary['qc_median_rsd']}",
    ]
    (step_dir / "preprocessing_summary.txt").write_text(
        "\n".join(text_lines) + "\n", encoding="utf-8",
    )
    print(
        f"    预处理完成: {len(preprocessing.processed_df)} 个代谢物；"
        f"归一化={cfg.normalization_method}；QC 移除={summary['removed_by_qc']}"
    )
    return TissueContext(
        tissue=tissue,
        tissue_label=tissue_label,
        processed_df=preprocessing.processed_df,
        sample_cols=preprocessing.sample_cols,
        qc_cols=preprocessing.qc_cols,
        groups=groups,
        qc_metrics=qc_table,
        sample_totals=preprocessing.sample_totals,
        filtered_csv=filtered_csv,
    )


def _load_preprocessed_context(
    cfg: MetabolomicsConfig,
    tissue: str,
    tissue_label: str,
    step_dir: Path,
) -> TissueContext:
    filtered_csv = step_dir / "filtered_data.csv"
    processed_df = pd.read_csv(filtered_csv)
    classified = classify_data_columns(processed_df.columns, set(cfg.mouse_groups))
    if classified.unmapped_biological:
        raise ValueError(
            f"缓存预处理矩阵含未配置样本列: {classified.unmapped_biological}"
        )
    groups = get_sample_groups(
        classified.biological,
        cfg.mouse_groups,
        cfg.conditions,
        strict=True,
    )
    qc_path = step_dir / "qc_feature_metrics.csv"
    totals_path = step_dir / "sample_total_signal.csv"
    qc_metrics = pd.read_csv(qc_path) if qc_path.exists() else pd.DataFrame()
    if totals_path.exists():
        sample_totals = pd.read_csv(totals_path).set_index("sample")
    else:
        sample_totals = pd.DataFrame()
    return TissueContext(
        tissue=tissue,
        tissue_label=tissue_label,
        processed_df=processed_df,
        sample_cols=classified.biological,
        qc_cols=classified.qc,
        groups=groups,
        qc_metrics=qc_metrics,
        sample_totals=sample_totals,
        filtered_csv=filtered_csv,
    )


def _plsda_signature(
    cfg: MetabolomicsConfig,
    context: TissueContext,
    comparisons: list[tuple[str, str, str]],
) -> dict[str, Any]:
    return {
        "pipeline_version": PIPELINE_VERSION,
        "filtered_mtime_ns": context.filtered_csv.stat().st_mtime_ns,
        "feature_count": len(context.processed_df),
        "sample_count": len(context.sample_cols),
        "comparisons": _comparison_names(comparisons),
        "plsda_permutation_num": cfg.plsda_permutation_num,
        "random_seed": cfg.random_seed,
    }


def _run_qc_and_plsda(
    cfg: MetabolomicsConfig,
    context: TissueContext,
    comparisons: list[tuple[str, str, str]],
    cond_map: dict[str, str],
    plot_group_colors: dict[str, str],
    step_root: Path,
    signature: dict[str, Any],
) -> tuple[np.ndarray | None, dict[str, np.ndarray]]:
    print("  [Step 2] QC 与 PLS-DA 验证")
    run_pca(
        context.processed_df,
        context.sample_cols,
        context.tissue,
        str(step_root),
        context.groups,
        plot_group_colors,
        qc_cols=context.qc_cols,
    )
    plot_preprocessing_qc(
        context.qc_metrics,
        context.sample_totals,
        context.tissue,
        str(step_root),
        qc_rsd_cutoff=cfg.qc_rsd_cutoff,
    )

    tissue_dir = step_root / context.tissue
    validation_rows: list[dict[str, Any]] = []
    global_vip: np.ndarray | None = None
    contrast_vip: dict[str, np.ndarray] = {}

    if len(set(context.groups.values())) >= 2:
        global_result = run_plsda(
            context.processed_df,
            context.sample_cols,
            context.tissue,
            str(step_root),
            context.groups,
            plot_group_colors,
            n_permutations=cfg.plsda_permutation_num,
            random_state=cfg.random_seed,
        )
        global_vip = global_result.vip.to_numpy(dtype=float)
        if global_result.status != "validated":
            global_vip[:] = np.nan
        validation_rows.append(
            {
                "comparison": "global_multiclass",
                "status": global_result.status,
                "R2Y": global_result.r2y,
                "Q2": global_result.q2,
                "CV_Accuracy": global_result.cv_accuracy,
                "Permutation_P": global_result.permutation_pvalue,
            }
        )

    for cond_a, cond_b, comp_name in comparisons:
        group_a = cond_map.get(cond_a, cond_a)
        group_b = cond_map.get(cond_b, cond_b)
        columns = [
            col for col in context.sample_cols
            if context.groups[col] in {group_a, group_b}
        ]
        if not any(context.groups[col] == group_a for col in columns) or not any(
            context.groups[col] == group_b for col in columns
        ):
            contrast_vip[comp_name] = np.full(len(context.processed_df), np.nan)
            validation_rows.append({"comparison": comp_name, "status": "missing_group"})
            continue
        result = run_plsda(
            context.processed_df,
            columns,
            context.tissue,
            str(step_root),
            context.groups,
            plot_group_colors,
            comp_name=comp_name,
            n_permutations=cfg.plsda_permutation_num,
            random_state=cfg.random_seed,
        )
        values = result.vip.to_numpy(dtype=float)
        if result.status != "validated":
            values[:] = np.nan
        contrast_vip[comp_name] = values
        validation_rows.append(
            {
                "comparison": comp_name,
                "status": result.status,
                "R2Y": result.r2y,
                "Q2": result.q2,
                "CV_Accuracy": result.cv_accuracy,
                "Permutation_P": result.permutation_pvalue,
            }
        )

    pd.DataFrame(validation_rows).to_csv(
        tissue_dir / "PLSDA_validation_summary.csv", index=False, encoding="utf-8-sig",
    )
    _write_json(tissue_dir / "step2_complete.json", signature)
    return global_vip, contrast_vip


def _read_validated_vip(path: Path, validation_path: Path, feature_count: int) -> np.ndarray:
    if not path.exists() or not validation_path.exists():
        return np.full(feature_count, np.nan)
    validation = pd.read_csv(validation_path)
    if validation.empty or validation.loc[0, "status"] != "validated":
        return np.full(feature_count, np.nan)
    values = pd.to_numeric(pd.read_csv(path)["VIP"], errors="coerce").to_numpy(dtype=float)
    if len(values) != feature_count:
        raise ValueError(f"VIP 特征数与预处理矩阵不一致: {path}")
    return values


def _load_plsda_outputs(
    context: TissueContext,
    comparisons: list[tuple[str, str, str]],
    step_root: Path,
) -> tuple[np.ndarray | None, dict[str, np.ndarray]]:
    tissue_dir = step_root / context.tissue
    feature_count = len(context.processed_df)
    global_vip = _read_validated_vip(
        tissue_dir / "PLS-DA_vip.csv",
        tissue_dir / "PLS-DA_validation.csv",
        feature_count,
    )
    contrast_vip = {
        comp_name: _read_validated_vip(
            tissue_dir / f"PLS-DA_{comp_name}_vip.csv",
            tissue_dir / f"PLS-DA_{comp_name}_validation.csv",
            feature_count,
        )
        for _, _, comp_name in comparisons
    }
    return global_vip, contrast_vip


def _differential_signature(
    cfg: MetabolomicsConfig,
    context: TissueContext,
    comparisons: list[tuple[str, str, str]],
    step2_marker: Path,
) -> dict[str, Any]:
    return {
        "pipeline_version": PIPELINE_VERSION,
        "filtered_mtime_ns": context.filtered_csv.stat().st_mtime_ns,
        "step2_mtime_ns": step2_marker.stat().st_mtime_ns,
        "feature_count": len(context.processed_df),
        "comparisons": _comparison_names(comparisons),
        "pval_cutoff": cfg.pval_cutoff,
        "fdr_cutoff": cfg.fdr_cutoff,
        "log2fc_cutoff": cfg.log2fc_cutoff,
        "vip_cutoff": cfg.vip_cutoff,
        "require_vip_for_significance": cfg.require_vip_for_significance,
        "fdr_method": cfg.fdr_method,
        "min_observations_per_group": cfg.min_observations_per_group,
    }


def _write_contrast_outputs(
    contrast_tables: dict[str, pd.DataFrame],
    targets: list[Path],
    tissue: str,
) -> None:
    for target in targets:
        target.mkdir(parents=True, exist_ok=True)
        for comp_name, table in contrast_tables.items():
            table.to_csv(
                target / f"dm_{tissue}_{comp_name}.csv",
                index=False,
                encoding="utf-8-sig",
            )


def _run_differential(
    cfg: MetabolomicsConfig,
    context: TissueContext,
    comparisons: list[tuple[str, str, str]],
    cond_map: dict[str, str],
    global_vip: np.ndarray | None,
    contrast_vip: dict[str, np.ndarray],
    step_root: Path,
    dm_root: Path,
    plot_colors: dict[str, Any],
    signature: dict[str, Any],
) -> pd.DataFrame:
    print("  [Step 3] 差异分析（Welch t + FDR）")
    annotated = annotate_metabolites(context.processed_df, name_col="Compounds")
    metadata_columns = [
        col for col in (
            "Index", "Compounds", "Formula", "kegg_map",
            "HMDB", "KEGG_Compound", "ChEBI",
        )
        if col in annotated.columns
    ]
    display_comparisons = [
        (cond_map.get(a, a), cond_map.get(b, b), comp_name)
        for a, b, comp_name in comparisons
    ]
    analysis = analyze_differential_metabolites(
        annotated,
        context.groups,
        display_comparisons,
        options=DifferentialAnalysisOptions(
            pvalue_cutoff=cfg.pval_cutoff,
            fdr_cutoff=cfg.fdr_cutoff,
            log2fc_cutoff=cfg.log2fc_cutoff,
            fdr_method=cfg.fdr_method,
            vip_cutoff=cfg.vip_cutoff,
            require_vip_for_significance=cfg.require_vip_for_significance,
            min_observations_per_group=cfg.min_observations_per_group,
        ),
        contrast_vip=contrast_vip,
        global_vip=global_vip,
        metadata_columns=metadata_columns,
    )

    tissue_dir = step_root / context.tissue
    tissue_dir.mkdir(parents=True, exist_ok=True)
    differential_csv = tissue_dir / "differential_analysis.csv"
    analysis.combined.to_csv(differential_csv, index=False, encoding="utf-8-sig")
    analysis.design_report.to_csv(
        tissue_dir / "differential_design_report.csv", index=False, encoding="utf-8-sig",
    )

    dm_dir = dm_root / context.tissue
    dm_dir.mkdir(parents=True, exist_ok=True)
    analysis.combined.to_csv(dm_dir / "differential_analysis.csv", index=False, encoding="utf-8-sig")
    _write_contrast_outputs(
        analysis.contrasts,
        [tissue_dir, dm_dir],
        context.tissue,
    )
    for _, _, comp_name in comparisons:
        if comp_name in analysis.contrasts:
            plot_volcano(
                analysis.combined,
                comp_name,
                context.tissue,
                str(step_root),
                pval_cutoff=cfg.pval_cutoff,
                fdr_cutoff=cfg.fdr_cutoff,
                log2fc_cutoff=cfg.log2fc_cutoff,
                colors=plot_colors,
            )
    _write_json(tissue_dir / "step3_complete.json", signature)
    print(f"    差异分析完成: {differential_csv}")
    return analysis.combined


def run_pipeline(cfg: MetabolomicsConfig, steps: list[int] | None = None) -> dict[str, Any]:
    """执行代谢组流程；默认增量运行，显式 ``--steps`` 会强制重跑所选步骤。"""
    selected, needed, forced = _resolve_steps(steps)
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    comparisons = _load_comparisons(cfg)
    if not comparisons:
        raise ValueError("未生成任何比较，请检查 conditions/contrast_groups")
    cond_map = get_cond_map(cfg.conditions)
    group_colors = get_group_colors(cfg.conditions)
    plot_group_colors = {cond_map[code]: color for code, color in group_colors.items()}
    plot_colors = get_plot_colors(cfg)
    tissue_labels = get_tissue_labels(cfg)

    print(f"\n{'=' * 68}")
    print("代谢组学分析流程 v2")
    print(f"实验: {cfg.experiment_id}")
    print(f"组织: {list(cfg.tissue_data_paths)}")
    print(f"比较: {_comparison_names(comparisons)}")
    print(f"归一化: {cfg.normalization_method}; 主显著性: FDR<{cfg.fdr_cutoff}")
    print(f"输出: {output_dir}")
    print(f"{'=' * 68}\n")

    completed: dict[str, list[int]] = {}
    for tissue, raw_path in cfg.tissue_data_paths.items():
        data_path = Path(raw_path)
        tissue_label = tissue_labels.get(tissue, tissue)
        completed[tissue] = []
        print(f"\n--- 组织: {tissue_label} ---")
        if not data_path.exists():
            print(f"  文件不存在: {data_path}，跳过该组织")
            continue

        step1_dir = output_dir / "01_data_preprocessing" / tissue
        step1_manifest = step1_dir / "preprocessing_summary.json"
        preprocessing_signature = _preprocessing_signature(cfg, data_path)
        step1_current = _manifest_matches(step1_manifest, preprocessing_signature)
        if 1 in needed and (1 in forced or not step1_current):
            context = _run_preprocessing(
                cfg, tissue, tissue_label, data_path, step1_dir, preprocessing_signature,
            )
            completed[tissue].append(1)
        else:
            context = _load_preprocessed_context(cfg, tissue, tissue_label, step1_dir)
            if 1 in selected:
                print("  [Step 1] 结果与输入/配置一致，跳过")

        step2_root = output_dir / "02_qc_visualization"
        step2_marker = step2_root / tissue / "step2_complete.json"
        plsda_signature = _plsda_signature(cfg, context, comparisons)
        step2_current = _manifest_matches(step2_marker, plsda_signature)
        if 2 in needed and (2 in forced or not step2_current):
            global_vip, contrast_vip = _run_qc_and_plsda(
                cfg,
                context,
                comparisons,
                cond_map,
                plot_group_colors,
                step2_root,
                plsda_signature,
            )
            completed[tissue].append(2)
        elif 2 in needed:
            global_vip, contrast_vip = _load_plsda_outputs(context, comparisons, step2_root)
            if 2 in selected:
                print("  [Step 2] 结果与预处理矩阵/配置一致，跳过")
        else:
            global_vip, contrast_vip = None, {}
            print("  [Step 2] 未选择，跳过")

        step3_root = output_dir / "03_differential_analysis"
        step3_dir = step3_root / tissue
        step3_marker = step3_dir / "step3_complete.json"
        differential_signature = _differential_signature(
            cfg, context, comparisons, step2_marker,
        )
        step3_current = _manifest_matches(step3_marker, differential_signature)
        if 3 in needed and (3 in forced or not step3_current):
            results = _run_differential(
                cfg,
                context,
                comparisons,
                cond_map,
                global_vip,
                contrast_vip,
                step3_root,
                output_dir / "00_dm_results",
                plot_colors,
                differential_signature,
            )
            completed[tissue].append(3)
        elif 3 in needed:
            results = pd.read_csv(step3_dir / "differential_analysis.csv")
            if 3 in selected:
                print("  [Step 3] 结果与输入/配置一致，跳过")
        else:
            results = pd.DataFrame()
            print("  [Step 3] 未选择，跳过")

        step4_root = output_dir / "04_heatmap"
        heatmap_path = step4_root / tissue / f"Heatmap_Top{cfg.heatmap_top_n}.png"
        if 4 in selected and (4 in forced or not heatmap_path.exists()):
            print("  [Step 4] 热力图")
            plot_heatmap(
                context.processed_df,
                context.sample_cols,
                context.groups,
                tissue,
                str(step4_root),
                top_n=cfg.heatmap_top_n,
                group_colors=plot_group_colors,
                selection_method=cfg.heatmap_selection_method,
            )
            completed[tissue].append(4)
        elif 4 in selected:
            print("  [Step 4] 已存在，跳过")
        else:
            print("  [Step 4] 未选择，跳过")

        step5_root = output_dir / "05_pathway_enrichment"
        step5_marker = step5_root / tissue / "step5_complete.json"
        step5_signature = {
            "pipeline_version": PIPELINE_VERSION,
            "differential_mtime_ns": (step3_dir / "differential_analysis.csv").stat().st_mtime_ns,
            "comparisons": _comparison_names(comparisons),
            "fdr_cutoff": cfg.fdr_cutoff,
            "fdr_method": cfg.fdr_method,
            "ora_min_set_size": cfg.ora_min_set_size,
            "ora_max_set_size": cfg.ora_max_set_size,
            "pathway_scope": cfg.pathway_scope,
            "exclude_overview_pathways": cfg.exclude_overview_pathways,
        } if 5 in needed else {}
        if 5 in selected and (5 in forced or not _manifest_matches(step5_marker, step5_signature)):
            print("  [Step 5] KEGG ORA（单侧 Fisher + FDR）")
            for _, _, comp_name in comparisons:
                enrichment = run_pathway_enrichment(
                    results,
                    comp_name,
                    tissue,
                    str(step5_root),
                    pval_cutoff=cfg.fdr_cutoff,
                    log2fc_cutoff=cfg.log2fc_cutoff,
                    fdr_method=cfg.fdr_method,
                    min_set_size=cfg.ora_min_set_size,
                    max_set_size=cfg.ora_max_set_size,
                    pathway_scope=cfg.pathway_scope,
                    exclude_overview=cfg.exclude_overview_pathways,
                )
                if enrichment is not None:
                    dm_dir = output_dir / "00_dm_results" / tissue
                    dm_dir.mkdir(parents=True, exist_ok=True)
                    enrichment.to_csv(
                        dm_dir / f"Pathway_{comp_name}.csv",
                        index=False,
                        encoding="utf-8-sig",
                    )
            _write_json(step5_marker, step5_signature)
            completed[tissue].append(5)
        elif 5 in selected:
            print("  [Step 5] 结果与差异分析/配置一致，跳过")
        else:
            print("  [Step 5] 未选择，跳过")

        step6_root = output_dir / "06_gsea"
        step6_marker = step6_root / tissue / "step6_complete.json"
        step6_signature = {
            "pipeline_version": PIPELINE_VERSION,
            "differential_mtime_ns": (step3_dir / "differential_analysis.csv").stat().st_mtime_ns,
            "comparisons": _comparison_names(comparisons),
            "min_set_size": cfg.gsea_min_set_size,
            "max_set_size": cfg.gsea_max_set_size,
            "permutation_num": cfg.gsea_permutation_num,
            "rank_metric": cfg.gsea_rank_metric,
            "pathway_scope": cfg.pathway_scope,
        } if 6 in needed else {}
        if 6 in selected and (6 in forced or not _manifest_matches(step6_marker, step6_signature)):
            print("  [Step 6] 代谢物集 GSEA")
            for _, _, comp_name in comparisons:
                run_gsea(
                    results,
                    comp_name,
                    tissue,
                    str(step6_root),
                    min_set_size=cfg.gsea_min_set_size,
                    max_set_size=cfg.gsea_max_set_size,
                    permutation_num=cfg.gsea_permutation_num,
                    rank_metric=cfg.gsea_rank_metric,
                    pathway_scope=cfg.pathway_scope,
                    exclude_overview=cfg.exclude_overview_pathways,
                )
            _write_json(step6_marker, step6_signature)
            completed[tissue].append(6)
        elif 6 in selected:
            print("  [Step 6] 结果与差异分析/配置一致，跳过")
        else:
            print("  [Step 6] 未选择，跳过")

    print(f"\n{'=' * 68}")
    print("代谢组学分析流程完成")
    print(f"结果目录: {output_dir}")
    print(f"{'=' * 68}\n")
    return {"output_dir": str(output_dir), "completed_steps": completed}


def main() -> None:
    """允许直接运行本文件；CLI 实现集中在唯一入口模块。"""
    from omicsweave.metabolomics.run_metabolomics_analysis import main as cli_main

    cli_main()


if __name__ == "__main__":
    main()
