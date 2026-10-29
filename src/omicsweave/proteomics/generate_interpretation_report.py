"""Generate config-driven, evidence-aware proteomics Markdown reports."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .enrichment.gsea import find_ssgsea_report  # noqa: E402
from .experiment_config import (  # noqa: E402
    config_arg_help,
    generate_comparisons,
    get_comp_labels,
    get_tissue_labels,
    load_proteomics_config,
)


OUTPUT_BASE = Path(".")
DP_BASE = Path(".")
PLOT_BASE = Path(".")
QC_BASE = Path(".")
REPORT_BASE = Path(".")
TISSUES: list[str] = []
COMPARISONS: list[tuple[str, str, str]] = []
TISSUE_LABELS: dict[str, str] = {}
COMP_LABELS: dict[str, str] = {}
PVAL_CUTOFF = 0.05
LOG2FC_CUTOFF = 0.58
DP_METHOD = "limma"
FDR_METHOD = "fdr_bh"


def flush_print(message: str) -> None:
    print(message, flush=True)


def apply_config(
    experiment: str | None,
    config_path: str | None,
    outdir: str | None,
) -> dict[str, str]:
    global OUTPUT_BASE, DP_BASE, PLOT_BASE, QC_BASE, REPORT_BASE
    global TISSUES, COMPARISONS, TISSUE_LABELS, COMP_LABELS
    global PVAL_CUTOFF, LOG2FC_CUTOFF, DP_METHOD, FDR_METHOD

    cfg = load_proteomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )
    OUTPUT_BASE = Path(cfg.output_dir)
    DP_BASE = OUTPUT_BASE / "01_dp_results"
    PLOT_BASE = OUTPUT_BASE / "02_result_plot"
    QC_BASE = OUTPUT_BASE / "00_qc"
    REPORT_BASE = OUTPUT_BASE / "03_interpretation_report"
    TISSUES = list(cfg.tissues)
    COMPARISONS = generate_comparisons(cfg.conditions, cfg.contrast_groups)
    TISSUE_LABELS = get_tissue_labels(cfg)
    COMP_LABELS = get_comp_labels(cfg.conditions, COMPARISONS)
    PVAL_CUTOFF = cfg.pval_cutoff
    LOG2FC_CUTOFF = cfg.log2fc_cutoff
    DP_METHOD = cfg.dp_method
    FDR_METHOD = cfg.fdr_method
    return {"experiment_id": cfg.experiment_id, "output_dir": str(OUTPUT_BASE)}


def load_dp(tissue: str, comparison: str) -> pd.DataFrame | None:
    path = DP_BASE / tissue / f"dp_{tissue}_{comparison}.csv"
    return pd.read_csv(path) if path.exists() else None


def load_combined_dp(tissue: str) -> pd.DataFrame | None:
    path = DP_BASE / tissue / "differential_proteins.csv"
    if path.exists():
        return pd.read_csv(path)
    legacy = DP_BASE / tissue / "differential_expression_wilcoxon.csv"
    return pd.read_csv(legacy) if legacy.exists() else None


def load_design_report(tissue: str) -> pd.DataFrame | None:
    path = DP_BASE / tissue / "differential_design_report.csv"
    return pd.read_csv(path) if path.exists() else None


def load_qc_summary(tissue: str) -> dict[str, Any] | None:
    path = QC_BASE / tissue / "qc_summary.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def load_ssgsea_results(tissue: str, library_key: str = "kegg") -> pd.DataFrame | None:
    """Load the version-dependent gseapy ssGSEA report filename."""

    report = find_ssgsea_report(PLOT_BASE / tissue / f"ssgsea_{library_key}")
    return pd.read_csv(report) if report is not None else None


def load_all_ssgsea_results(tissue: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    tissue_plot_dir = PLOT_BASE / tissue
    if tissue_plot_dir.exists():
        for directory in sorted(tissue_plot_dir.glob("ssgsea_*")):
            report = find_ssgsea_report(directory)
            if report is None:
                continue
            frame = pd.read_csv(report)
            frame.insert(0, "library", directory.name.removeprefix("ssgsea_"))
            frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def load_prerank_results(tissue: str, comparison: str) -> pd.DataFrame:
    root = PLOT_BASE / tissue / "gsea_prerank" / comparison
    frames: list[pd.DataFrame] = []
    if root.exists():
        for report in sorted(root.glob("*/*prerank*.report.csv")):
            frame = pd.read_csv(report)
            frame.insert(0, "library", report.parent.name)
            frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def validate_corrected_results(tissue: str) -> None:
    """Refuse to report from pre-v6 outputs that may contain unadjusted padj."""

    combined = DP_BASE / tissue / "differential_proteins.csv"
    if not combined.exists():
        return
    metadata_path = DP_BASE / tissue / "analysis_metadata.json"
    if not metadata_path.exists():
        raise RuntimeError(
            f"[{tissue}] differential results have no v6 analysis metadata. "
            "Re-run run_proteomics_pipeline.py with --force before reporting."
        )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if str(metadata.get("pipeline_version")) != "6.0":
        raise RuntimeError(
            f"[{tissue}] results were produced by pipeline "
            f"{metadata.get('pipeline_version')!r}; re-run with --force."
        )


def summarize_dp(df: pd.DataFrame) -> dict[str, Any]:
    """Summarize only FDR-controlled inferential results."""

    work = df.copy()
    padj = pd.to_numeric(work.get("padj"), errors="coerce")
    fold_change = pd.to_numeric(work.get("log2FoldChange"), errors="coerce")
    if "evidence_level" in work.columns:
        inferential = work["evidence_level"].astype(str) == "inferential"
    else:
        inferential = pd.Series(True, index=work.index)
    significant = inferential & (padj < PVAL_CUTOFF) & (fold_change.abs() > LOG2FC_CUTOFF)
    up = significant & (fold_change > 0)
    down = significant & (fold_change < 0)

    ranked = work.assign(_padj=padj, _fc=fold_change, _inferential=inferential)
    ranked = ranked[ranked["_inferential"] & ranked["_padj"].notna()]
    top_up = ranked[ranked["_fc"] > 0].nsmallest(5, "_padj")["Protein"].tolist()
    top_down = ranked[ranked["_fc"] < 0].nsmallest(5, "_padj")["Protein"].tolist()
    return {
        "total_proteins": int(len(work)),
        "n_inferential": int(inferential.sum()),
        "n_insufficient": int((~inferential).sum()),
        "sig_up": int(up.sum()),
        "sig_down": int(down.sum()),
        "sig_total": int(significant.sum()),
        "top_up": top_up,
        "top_down": top_down,
    }


def _format_number(value: Any, precision: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "N/A"
    if not np.isfinite(number):
        return "N/A"
    return f"{number:.{precision}g}"


def _prerank_summary_table(frame: pd.DataFrame, top_n: int = 10) -> list[str]:
    if frame.empty:
        return ["*未找到可用的 pre-ranked GSEA 结果。*"]
    term_col = next((c for c in ("Term", "Name", "Term_name") if c in frame.columns), None)
    nes_col = next((c for c in ("NES", "nes") if c in frame.columns), None)
    fdr_col = next(
        (c for c in ("FDR q-val", "FDR", "padj", "Adjusted P-value") if c in frame.columns),
        None,
    )
    pvalue_col = next((c for c in ("NOM p-val", "pvalue", "P-value") if c in frame.columns), None)
    if term_col is None:
        return ["*GSEA 报告缺少通路名称列。*"]
    work = frame.copy()
    sort_col = fdr_col or pvalue_col
    if sort_col:
        work[sort_col] = pd.to_numeric(work[sort_col], errors="coerce")
        work = work.sort_values(sort_col)
    lines = [
        "| 基因集库 | 通路 | NES | FDR | 名义 P 值 |",
        "|----------|------|-----|-----|------------|",
    ]
    for row in work.head(top_n).to_dict(orient="records"):
        lines.append(
            f"| {row.get('library', 'N/A')} | {str(row.get(term_col, '')).replace('|', '/')} | "
            f"{_format_number(row.get(nes_col)) if nes_col else 'N/A'} | "
            f"{_format_number(row.get(fdr_col)) if fdr_col else 'N/A'} | "
            f"{_format_number(row.get(pvalue_col)) if pvalue_col else 'N/A'} |"
        )
    return lines


def _ssgsea_variability_table(frame: pd.DataFrame, top_n: int = 10) -> list[str]:
    if frame.empty:
        return ["*未找到样本级 ssGSEA 结果。*"]
    term_col = next((c for c in ("Term", "Term_name") if c in frame.columns), None)
    score_col = next((c for c in ("NES", "ES", "enrichment_score") if c in frame.columns), None)
    if term_col is None or score_col is None:
        return ["*ssGSEA 报告缺少通路或分数列。*"]
    work = frame.copy()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    summary = (
        work.dropna(subset=[score_col])
        .groupby(["library", term_col], as_index=False)[score_col]
        .agg(["min", "max", "mean"])
        .reset_index()
    )
    summary["range"] = summary["max"] - summary["min"]
    summary = summary.sort_values("range", ascending=False).head(top_n)
    lines = [
        "| 基因集库 | 通路 | 样本间范围 | 平均分数 |",
        "|----------|------|------------|----------|",
    ]
    for row in summary.to_dict(orient="records"):
        lines.append(
            f"| {row['library']} | {str(row[term_col]).replace('|', '/')} | "
            f"{_format_number(row['range'])} | {_format_number(row['mean'])} |"
        )
    return lines


def tissue_page(tissue: str) -> str:
    label = TISSUE_LABELS.get(tissue, tissue)
    combined = load_combined_dp(tissue)
    qc = load_qc_summary(tissue)
    design = load_design_report(tissue)
    lines = [f"# {label} 蛋白质组学分析报告", "", "## 1. 数据与 QC", ""]
    if combined is None:
        lines.append("**数据未找到**：请先运行蛋白质组学主流程。")
        return "\n".join(lines)

    lines.append(f"- 进入差异分析的蛋白数：**{len(combined)}**")
    if qc:
        lines.extend([
            f"- 样本数：**{qc.get('n_samples', 'N/A')}**",
            f"- 总体缺失率：**{_format_number(qc.get('overall_missing_pct'))}%**",
            f"- 组内相关性中位数：**{_format_number(qc.get('median_within_group_correlation'))}**",
        ])
        flags = qc.get("flags", [])
        lines.append(f"- QC 标记：**{', '.join(flags) if flags else '无自动标记'}**")
    lines.append("")

    lines.extend([
        f"## 2. 差异蛋白分析（{DP_METHOD}）",
        "",
        f"判定规则：{FDR_METHOD} 校正后 P < {PVAL_CUTOFF}，且 |Log2FC| > {LOG2FC_CUTOFF}。",
        "不可识别设计、每组少于 2 个有效观测或残差自由度不足的蛋白仅报告效应量，不进入显著性计数。",
        "",
        "| 比较组 | 可推断蛋白 | 证据不足 | 显著上调 | 显著下调 | 显著总数 |",
        "|--------|------------|----------|----------|----------|----------|",
    ])
    summaries: dict[str, dict[str, Any]] = {}
    for _numerator, _denominator, comparison in COMPARISONS:
        frame = load_dp(tissue, comparison)
        comparison_label = COMP_LABELS.get(comparison, comparison)
        if frame is None:
            lines.append(f"| {comparison_label} | - | - | - | - | - |")
            continue
        summary = summarize_dp(frame)
        summaries[comparison] = summary
        lines.append(
            f"| {comparison_label} | {summary['n_inferential']} | "
            f"{summary['n_insufficient']} | {summary['sig_up']} | "
            f"{summary['sig_down']} | {summary['sig_total']} |"
        )
    lines.append("")

    if design is not None and not design.empty:
        lines.extend([
            "### 设计可识别性",
            "",
            "| 比较组 | 样本数（分子/分母） | 设计式 | 路径 | 可检验蛋白 | 原因 |",
            "|--------|----------------------|--------|------|------------|------|",
        ])
        for row in design.to_dict(orient="records"):
            lines.append(
                f"| {COMP_LABELS.get(str(row.get('comparison')), row.get('comparison'))} | "
                f"{row.get('n_numerator')}/{row.get('n_denominator')} | "
                f"{row.get('design_formula')} | {row.get('route')} | "
                f"{row.get('n_features_tested')} | {str(row.get('reason', '')).replace('|', '/')} |"
            )
        lines.append("")

    lines.extend(["## 3. Top 差异蛋白", ""])
    for _numerator, _denominator, comparison in COMPARISONS:
        frame = load_dp(tissue, comparison)
        if frame is None:
            continue
        summary = summaries.get(comparison) or summarize_dp(frame)
        lines.extend([
            f"### {COMP_LABELS.get(comparison, comparison)}",
            "",
            f"- 上调 Top 5：{', '.join(summary['top_up']) if summary['top_up'] else '无'}",
            f"- 下调 Top 5：{', '.join(summary['top_down']) if summary['top_down'] else '无'}",
            "",
        ])
        inferential = frame.get("evidence_level", "inferential")
        if isinstance(inferential, pd.Series):
            inferential_mask = inferential.astype(str) == "inferential"
        else:
            inferential_mask = pd.Series(True, index=frame.index)
        significant = frame[
            inferential_mask
            & (pd.to_numeric(frame["padj"], errors="coerce") < PVAL_CUTOFF)
            & (pd.to_numeric(frame["log2FoldChange"], errors="coerce").abs() > LOG2FC_CUTOFF)
        ].sort_values("padj")
        if significant.empty:
            lines.extend(["*无满足 FDR 与效应量双阈值的蛋白。*", ""])
            continue
        lines.extend([
            "| 蛋白 | Log2FC | P 值 | Padj | 方向 |",
            "|------|--------|------|------|------|",
        ])
        for row in significant.head(20).to_dict(orient="records"):
            direction = "Up" if float(row["log2FoldChange"]) > 0 else "Down"
            lines.append(
                f"| {row['Protein']} | {_format_number(row['log2FoldChange'])} | "
                f"{_format_number(row['pvalue'])} | {_format_number(row['padj'])} | {direction} |"
            )
        lines.append("")

    lines.extend([
        "## 4. Pre-ranked GSEA（比较层推断）",
        "",
        "GSEA 使用全部具有有限 moderated t 统计量的蛋白排序，不以差异蛋白阈值预筛选。",
        "",
    ])
    for _numerator, _denominator, comparison in COMPARISONS:
        lines.extend([
            f"### {COMP_LABELS.get(comparison, comparison)}",
            "",
            *_prerank_summary_table(load_prerank_results(tissue, comparison)),
            "",
        ])

    lines.extend([
        "## 5. 样本级 ssGSEA（描述性）",
        "",
        "下表按样本间分数范围排序，仅用于发现异质性较大的通路，不代表组间显著性。",
        "",
        *_ssgsea_variability_table(load_all_ssgsea_results(tissue)),
        "",
        "## 6. 统计与数据处理说明",
        "",
        f"- 差异分析：**{DP_METHOD}**；多重检验：**{FDR_METHOD}**。",
        "- 差异检验输入：log2 转换并按样本中位数中心化的**未插补**矩阵。",
        "- 下移正态分布插补仅用于 PCA、热图和样本级 ssGSEA。",
        "- ssGSEA 是样本级描述性分数；组间通路结论以 pre-ranked GSEA 为准。",
        "- 每个比较的样本量、设计秩、残差自由度和可检验蛋白数见设计可识别性表。",
        "",
    ])
    return "\n".join(lines)


def generate_index(
    experiment_id: str,
    tissue_summaries: dict[str, dict[str, Any]],
) -> str:
    comparison_labels = [COMP_LABELS.get(item[2], item[2]) for item in COMPARISONS]
    lines = [
        "# 蛋白质组学分析总报告",
        "",
        f"实验：**{experiment_id}**",
        "",
        "## 各组织差异蛋白概览",
        "",
        "| 组织 | 检测蛋白数 | " + " | ".join(comparison_labels) + " |",
        "|------|------------|" + "|".join(["----------"] * len(COMPARISONS)) + "|",
    ]
    for tissue in TISSUES:
        info = tissue_summaries.get(tissue, {})
        cells = [
            str(info.get("comparisons", {}).get(comparison, {}).get("sig_total", "N/A"))
            for _num, _den, comparison in COMPARISONS
        ]
        lines.append(
            f"| {TISSUE_LABELS.get(tissue, tissue)} | {info.get('total_proteins', 'N/A')} | "
            + " | ".join(cells)
            + " |"
        )
    lines.extend(["", "## 详细报告", ""])
    for tissue in TISSUES:
        lines.append(
            f"- [{TISSUE_LABELS.get(tissue, tissue)}](./{tissue}/{tissue}_report.md)"
        )
    lines.extend(["", "---", "", "*本报告由证据感知的蛋白质组学报告模块自动生成。*", ""])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate corrected proteomics interpretation reports.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment ID")
    parser.add_argument("--config", default=None, help="Path to config.yaml")
    parser.add_argument("--outdir", default=None, help="Override output directory")
    args = parser.parse_args()

    config = apply_config(args.experiment, args.config, args.outdir)
    REPORT_BASE.mkdir(parents=True, exist_ok=True)
    tissue_summaries: dict[str, dict[str, Any]] = {}
    for tissue in TISSUES:
        validate_corrected_results(tissue)
        tissue_dir = REPORT_BASE / tissue
        tissue_dir.mkdir(parents=True, exist_ok=True)
        report_path = tissue_dir / f"{tissue}_report.md"
        report_path.write_text(tissue_page(tissue), encoding="utf-8")

        combined = load_combined_dp(tissue)
        comparisons: dict[str, dict[str, Any]] = {}
        for _num, _den, comparison in COMPARISONS:
            frame = load_dp(tissue, comparison)
            if frame is not None:
                comparisons[comparison] = summarize_dp(frame)
        tissue_summaries[tissue] = {
            "total_proteins": len(combined) if combined is not None else 0,
            "comparisons": comparisons,
        }
        flush_print(f"Saved: {report_path}")

    index = REPORT_BASE / "index.md"
    index.write_text(
        generate_index(config["experiment_id"], tissue_summaries), encoding="utf-8",
    )
    flush_print(f"Saved: {index}")


if __name__ == "__main__":
    main()
