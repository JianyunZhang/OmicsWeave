#!/usr/bin/env python3
"""
Cross-tissue pathway integration runner.

Loads curated pathway evidence (cross_tissue_evidence.yaml) derived from the
five tissue-level multiomics reports, scores each systemic pathway by
cross-tissue + cross-omics convergence, and emits:

  multiomics_results/
    common/cross_tissue_evidence_matrix.csv   (pathway x tissue evidence)
    common/pathway_integration_scores.csv     (evidence profile + legacy score)
    cross_tissue_pathway_report.md            (cross-tissue curated summary)

Usage:
    python run_cross_tissue_integration.py --experiment 20250613_cimr_wt_mice_dio_glp

This is a pathway-centric integration for transcriptomics + metabolomics data
(no MFA / proteomics). The MFA-gated flux_omics_scorer is NOT used.

Author: Claude Code
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


from omicsweave.multiomics.config import repo_root
from omicsweave.multiomics.integration.cross_tissue_pathway_scorer import CrossTissuePathwayScorer


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _evidence_yaml_path(experiment_id: str) -> Path:
    return (
        repo_root()
        / "experiments"
        / experiment_id
        / "multiomics_integration_setting"
        / "cross_tissue_evidence.yaml"
    )


def _output_dir(experiment_id: str, outdir: str | Path | None = None) -> Path:
    odir = (
        Path(outdir).resolve()
        if outdir is not None
        else repo_root() / "analysis_results" / experiment_id / "multiomics_results"
    )
    (odir / "common").mkdir(parents=True, exist_ok=True)
    return odir


def _load_evidence(yaml_path: Path) -> tuple[list[str], list[dict], list[str]]:
    data = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Cross-tissue evidence must be a YAML mapping: {yaml_path}")
    tissues = data.get("tissues", [])
    pathways = data.get("pathways", [])
    if not isinstance(tissues, list) or not all(isinstance(value, str) for value in tissues):
        raise ValueError("cross_tissue_evidence.yaml: tissues must be a list of names")
    if not isinstance(pathways, list) or not pathways:
        raise ValueError("cross_tissue_evidence.yaml: pathways must be a non-empty list")
    comparisons = sorted({
        str(tissue_evidence.get("comparison"))
        for pathway in pathways
        for tissue_evidence in (pathway.get("evidence", {}) or {}).values()
        if tissue_evidence.get("comparison")
    })
    return tissues, pathways, comparisons


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def _fmt_conf_distribution(scores: pd.DataFrame) -> str:
    counts = scores["confidence"].value_counts()
    order = ["High", "Medium", "Low", "Insufficient"]
    parts = []
    for c in order:
        n = int(counts.get(c, 0))
        parts.append(f"{c}={n}")
    return ", ".join(parts)


def _tissue_evidence_table(matrix: pd.DataFrame) -> str:
    """Render the pathway x tissue matrix as a markdown table."""
    cols = ["id", "pathway"] + [c for c in matrix.columns if c not in ("id", "pathway")]
    md = matrix[cols].to_markdown(index=False)
    return md


def generate_report(
    experiment_id: str,
    scores: pd.DataFrame,
    matrix: pd.DataFrame,
    tissues: list[str],
    output_dir: Path,
    comparisons: list[str] | None = None,
) -> Path:
    report_path = output_dir / "cross_tissue_pathway_report.md"
    n_pathways = len(scores)
    n_tissues = len(tissues)
    conf_dist = _fmt_conf_distribution(scores)
    grade_dist = ", ".join(
        f"{grade}={int(count)}"
        for grade, count in scores["evidence_grade"].value_counts().items()
    )

    high = scores[scores["confidence"] == "High"]
    medium = scores[scores["confidence"] == "Medium"]
    source_reports = sorted({
        value.strip()
        for raw in scores.get("source_reports", pd.Series(dtype=str)).dropna().astype(str)
        for value in raw.split(",")
        if value.strip()
    })

    lines: list[str] = []
    lines.append("# 跨组织 multiomics 通路整合报告\n")
    lines.append(f"**实验**: {experiment_id}  \n")
    lines.append(f"**组织** ({n_tissues}): {', '.join(tissues)}  \n")
    lines.append(f"**比较范围**: {' / '.join(comparisons or ['未在证据文件中声明'])}  \n")
    lines.append(f"**整合策略**: 通路中心探索性证据画像（转录组+代谢组，无 MFA/蛋白组）  \n")
    lines.append(f"**生成时间**: {pd.Timestamp.now().strftime('%Y-%m-%d')}\n")
    lines.append("---\n")

    lines.append("## 一、概述\n")
    lines.append(
        f"本报告基于证据文件声明的 **{len(source_reports)} 个来源标识**，覆盖 "
        f"**{n_tissues} 个组织**与 **{len(comparisons or [])} 个比较**，对 "
        f"**{n_pathways} 条系统性跨组织通路**生成多维证据画像。"
        "由于输入是人工策展代表性条目，缺少 SE、qvalue、n_eff 和重采样信息，所有可评估结果"
        "均标记为 design grade D / Exploratory。旧 0–4 分仅保留为 legacy 兼容列。\n"
    )
    lines.append(f"**证据等级分布**: {grade_dist}\n")
    lines.append(f"**旧分级分布（仅兼容）**: {conf_dist}\n")
    lines.append("---\n")

    lines.append("## 二、证据画像方法\n")
    lines.append("主要字段为 support/conflict/coverage、组织复现比例、中位绝对效应、design grade 和 evidence grade。support 使用组织内及跨组织中位效应，避免单个 `max_log2fc` 异常值主导。\n")
    lines.append("以下 0–4 规则只解释 `legacy_*` 列；组织广度与新颖性重复奖励组织数量，因此不得称为独立置信度。\n")
    lines.append("| 维度 | 分值 | 计算 |")
    lines.append("|:---|:---:|:---|")
    lines.append("| 组织广度 | 0–1.5 | 每个有证据组织 +0.3，上限 1.5 |")
    lines.append("| 跨组学一致性 | 0–1.5 | 转录+代谢同向 ≥1 组织=1.0, ≥2 组织=1.5；仅单组学=0.5 |")
    lines.append("| 效应量级 | 0–0.5 | max\\|log2FC\\|: ≥2.0→0.5, ≥1.0→0.4, ≥0.58→0.3, 其他→0.1 |")
    lines.append("| 跨组织新颖性 | 0–0.5 | 跨 ≥2 组织=0.5；单组织=0.2 |")
    lines.append(
        "\n> 注：对 GSEA 通路级证据，用 NES 作为量级代理（已在 evidence.yaml 标注 feature=GSEA_NES）。\n"
    )
    lines.append("---\n")

    lines.append("## 三、系统性通路评分总表\n")
    show_cols = [
        "id", "pathway", "category", "comparison_sources", "support_score",
        "conflict_score", "coverage_score", "tissue_recurrence",
        "median_abs_log2fc", "design_grade", "evidence_grade",
        "legacy_total_score", "legacy_confidence", "n_tissues", "tissues_present",
    ]
    show_cols = [c for c in show_cols if c in scores.columns]
    lines.append(scores[show_cols].to_markdown(index=False))
    lines.append("")
    lines.append("---\n")

    lines.append("## 四、跨组织证据矩阵（通路 × 组织）\n")
    lines.append("方向 + 代表性特征（log2FC）。『—』表示该组织无显著证据。\n")
    lines.append(_tissue_evidence_table(matrix))
    lines.append("")
    lines.append("---\n")

    lines.append("## 五、Legacy High 通路详述（不代表确认性置信度）\n")
    if high.empty:
        lines.append("（无 High 置信通路）\n")
    else:
        for _, r in high.iterrows():
            lines.append(f"### {r['id']} {r['pathway']} — legacy {r['legacy_total_score']} ({r['legacy_confidence']})\n")
            lines.append(f"- **新画像**: support={r['support_score']}, conflict={r['conflict_score']}, coverage={r['coverage_score']}, {r['evidence_grade']} (design {r['design_grade']})")
            lines.append(f"- **类别**: {r['category']} | **主导药物**: {r['drug']}")
            lines.append(f"- **跨组织广度**: {r['n_tissues']} 组织 ({r['tissues_present']})")
            lines.append(f"- **评分明细**: 广度 {r['tissue_breadth']} + 一致性 {r['omics_concordance']} "
                         f"(跨组学同向 {r['n_concordant_tissues']} 组织) + 量级 {r['effect_magnitude']} "
                         f"(max log2FC {r['max_log2fc']}) + 新颖性 {r['novelty']}")
            lines.append(f"- **描述**: {r['description']}")
            lines.append(f"- **来源报告**: {r['source_reports']}\n")
    lines.append("---\n")

    lines.append("## 六、Legacy Medium 通路（仅兼容）\n")
    if medium.empty:
        lines.append("（无 Medium 置信通路）\n")
    else:
        for _, r in medium.iterrows():
            lines.append(f"- **{r['id']} {r['pathway']}** — legacy {r['legacy_total_score']} "
                         f"({r['n_tissues']} 组织: {r['tissues_present']}) — {r['description']}")
    lines.append("\n---\n")

    lines.append("## 七、按证据画像排序的探索性发现\n")
    for index, (_, row) in enumerate(scores.head(6).iterrows(), start=1):
        lines.append(
            f"{index}. **{row['id']} {row['pathway']}**：{row['description']} "
            f"（组织={row['n_tissues']}；support={row['support_score']}；"
            f"conflict={row['conflict_score']}；coverage={row['coverage_score']}；"
            f"比较={row.get('comparison_sources', '') or '未声明'}）。\n"
        )
    lines.append("---\n")

    lines.append("## 八、输入证据的引用链\n")
    lines.append("| 来源标识 | 声明位置 |")
    lines.append("|:---|:---|")
    for source in source_reports:
        lines.append(f"| {source} | `cross_tissue_evidence.yaml: source_reports` |")
    lines.append(f"\n证据源: `experiments/{experiment_id}/multiomics_integration_setting/cross_tissue_evidence.yaml`\n")
    lines.append("---\n")

    lines.append("## 九、局限性\n")
    lines.append(f"- 证据为策展定性（{len(source_reports)} 个来源标识的代表性条目），非全特征 raw 数据重整合\n")
    lines.append(f"- 本报告只覆盖证据文件声明的 {len(comparisons or [])} 个比较；未声明的比较不得推断为阴性\n")
    lines.append("- 当前入口不消费 MFA/蛋白组；涉及通量或酶资源的机制需独立数据验证\n")
    lines.append("- GSEA 通路级证据用 NES 作量级代理（非严格 log2FC）\n")
    lines.append("- 输入未提供生物学重复数、SE、qvalue 或重采样结果，不能生成 Replicated 等级\n")
    lines.append("---\n")

    lines.append("## 十、后续建议\n")
    lines.append("1. 为每个策展条目回填 effect、SE、qvalue、n_eff、mapping_quality 与样本/批次来源\n")
    lines.append("2. 在生物学重复层面执行 bootstrap 与标签置换，校准选择概率和方向稳定性\n")
    lines.append("3. 对统一 estimand 的独立实验使用项目级随机效应 Meta，而不是比较 legacy 总分\n")
    lines.append("4. 对排名靠前的探索性通路设计靶向测量或同位素示踪进行机制验证\n")
    lines.append("---\n")

    lines.append("*本报告由 `run_cross_tissue_integration.py` 自动生成（评分器: "
                 "`omicsweave/multiomics/integration/cross_tissue_pathway_scorer.py`）*\n")

    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(experiment_id: str, outdir: str | Path | None = None) -> Path:
    yaml_path = _evidence_yaml_path(experiment_id)
    if not yaml_path.exists():
        raise FileNotFoundError(f"Evidence YAML not found: {yaml_path}")

    output_dir = _output_dir(experiment_id, outdir)
    tissues, pathways, comparisons = _load_evidence(yaml_path)

    print("=" * 70)
    print("Cross-Tissue Pathway Integration")
    print("=" * 70)
    print(f"Experiment: {experiment_id}")
    print(f"Evidence:   {yaml_path}")
    print(f"Output:     {output_dir}")
    print(f"Tissues:    {', '.join(tissues)}")
    print(f"Pathways:   {len(pathways)}")
    print()

    scorer = CrossTissuePathwayScorer(tissues=tissues)
    scores_df = scorer.score_all(pathways)
    matrix_df = scorer.build_matrix(pathways)

    # Write CSVs
    scores_csv = output_dir / "common" / "pathway_integration_scores.csv"
    matrix_csv = output_dir / "common" / "cross_tissue_evidence_matrix.csv"
    scores_df.to_csv(scores_csv, index=False, encoding="utf-8-sig")
    matrix_df.to_csv(matrix_csv, index=False, encoding="utf-8-sig")
    print(f"[OK] Scores:   {scores_csv}")
    print(f"[OK] Matrix:   {matrix_csv}")

    # Report
    report_path = generate_report(
        experiment_id, scores_df, matrix_df, tissues, output_dir, comparisons
    )
    print(f"[OK] Report:   {report_path}")

    # Console summary
    print()
    print("Evidence grade distribution:")
    print(scores_df["evidence_grade"].value_counts().to_string())
    print()
    print("Top pathways:")
    print(scores_df[["id", "pathway", "support_score", "conflict_score", "coverage_score", "evidence_grade"]].head(6).to_string(index=False))
    return report_path


def main() -> None:
    p = argparse.ArgumentParser(description="Cross-tissue pathway integration scorer")
    p.add_argument("--experiment", required=True, help="Experiment ID")
    p.add_argument("--outdir", default=None, help="Optional isolated output directory")
    args = p.parse_args()
    run(args.experiment, args.outdir)


if __name__ == "__main__":
    main()
