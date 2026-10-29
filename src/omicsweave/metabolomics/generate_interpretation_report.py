"""
代谢组学分析解读报告生成器。

为每个组织生成 Markdown 摘要，并汇总为 HTML 报告。
"""
from __future__ import annotations

import base64
import json
import os
import sys
from pathlib import Path

import pandas as pd

from .experiment_config import (
    MetabolomicsConfig,
    generate_comparisons,
    get_comp_labels,
    get_cond_map,
    get_tissue_labels,
)


def _embed_image(img_path: str) -> str:
    """将图片转为 base64 嵌入 HTML。"""
    if not os.path.exists(img_path):
        return f"<p><em>图片未找到: {img_path}</em></p>"
    with open(img_path, "rb") as f:
        data = base64.b64encode(f.read()).decode("utf-8")
    ext = Path(img_path).suffix.lstrip(".")
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg"}.get(ext, "image/png")
    return f'<img src="data:{mime};base64,{data}" style="max-width:800px;width:100%;" />'


def _significance_masks(
    diff_df: pd.DataFrame,
    comp_name: str,
    pval_cutoff: float,
    fdr_cutoff: float,
    log2fc_cutoff: float,
) -> tuple[pd.Series, pd.Series]:
    """返回唯一主显著性与名义 P 显著性掩码。"""
    p_col = f"Pval_{comp_name}"
    fdr_col = f"FDR_{comp_name}"
    fc_col = f"Log2FC_{comp_name}"
    sig_col = f"Significant_{comp_name}"
    empty = pd.Series(False, index=diff_df.index)
    if p_col not in diff_df.columns or fc_col not in diff_df.columns:
        return empty, empty
    fc_pass = diff_df[fc_col].abs() > log2fc_cutoff
    nominal = (diff_df[p_col] < pval_cutoff) & fc_pass
    if sig_col in diff_df.columns:
        primary = diff_df[sig_col].fillna(False).astype(bool)
    elif fdr_col in diff_df.columns:
        primary = (diff_df[fdr_col] < fdr_cutoff) & fc_pass
    else:
        primary = nominal
    return primary, nominal


def _count_sig(
    diff_df: pd.DataFrame,
    comp_name: str,
    pval_cutoff: float,
    fdr_cutoff: float,
    log2fc_cutoff: float,
) -> tuple[int, int, int, int]:
    """统计 FDR 主发现与原始 P 名义发现。"""
    fc_col = f"Log2FC_{comp_name}"
    if fc_col not in diff_df.columns:
        return 0, 0, 0, 0
    sig, nominal = _significance_masks(
        diff_df, comp_name, pval_cutoff, fdr_cutoff, log2fc_cutoff,
    )
    up = (sig) & (diff_df[fc_col] > 0)
    down = (sig) & (diff_df[fc_col] < 0)
    return int(up.sum()), int(down.sum()), int(sig.sum()), int(nominal.sum())


def _top_metabolites_table(
    diff_df: pd.DataFrame,
    comp_name: str,
    direction: str,
    pval_cutoff: float,
    fdr_cutoff: float,
    log2fc_cutoff: float,
    n: int = 5,
) -> str:
    """生成 Top N 代谢物 HTML 表格。"""
    p_col = f"Pval_{comp_name}"
    fdr_col = f"FDR_{comp_name}"
    fc_col = f"Log2FC_{comp_name}"
    if p_col not in diff_df.columns:
        return "<p>无数据</p>"

    primary, _ = _significance_masks(
        diff_df, comp_name, pval_cutoff, fdr_cutoff, log2fc_cutoff,
    )
    direction_mask = diff_df[fc_col] > 0 if direction == "up" else diff_df[fc_col] < 0
    sort_col = fdr_col if fdr_col in diff_df.columns else p_col
    subset = diff_df[primary & direction_mask].nsmallest(n, sort_col)

    if subset.empty:
        return "<p>无显著代谢物</p>"

    rows = ""
    for _, r in subset.iterrows():
        fdr_value = r.get(fdr_col, float("nan"))
        rows += (
            f"<tr><td>{r.get('Compounds', r.get('Index', ''))}</td>"
            f"<td>{r[fc_col]:.3f}</td><td>{r[p_col]:.2e}</td>"
            f"<td>{fdr_value:.2e}</td></tr>"
        )
    return f"<table><tr><th>代谢物</th><th>Log2FC</th><th>P-value</th><th>FDR</th></tr>{rows}</table>"


def generate_tissue_summary(
    tissue: str,
    tissue_label: str,
    cfg: MetabolomicsConfig,
    output_dir: str,
    comparisons: list,
) -> str:
    """为单个组织生成 Markdown 摘要。"""
    sections = []
    sections.append(f"## {tissue_label} 代谢组学分析总结\n")

    # ── 数据概况 ──
    filtered_csv = os.path.join(output_dir, "01_data_preprocessing", tissue, "filtered_data.csv")
    summary_txt = os.path.join(output_dir, "01_data_preprocessing", tissue, "preprocessing_summary.txt")

    if os.path.exists(summary_txt):
        with open(summary_txt, "r", encoding="utf-8") as f:
            sections.append("### 1. 数据概况\n")
            sections.append("```")
            sections.append(f.read())
            sections.append("```\n")
    elif os.path.exists(filtered_csv):
        df = pd.read_csv(filtered_csv)
        sections.append("### 1. 数据概况\n")
        sections.append(f"- 过滤后代谢物数: {len(df)}\n")

    # ── QC 聚类 ──
    pca_png = os.path.join(output_dir, "02_qc_visualization", tissue, "PCA.png")
    if os.path.exists(pca_png):
        sections.append("### 2. QC/聚类质量\n")
        summary_json = os.path.join(output_dir, "01_data_preprocessing", tissue, "preprocessing_summary.json")
        if os.path.exists(summary_json):
            with open(summary_json, "r", encoding="utf-8") as handle:
                qc_summary = json.load(handle)
            median_rsd = qc_summary.get("qc_median_rsd")
            sections.append(f"- 归一化方法: {qc_summary.get('normalization_method', cfg.normalization_method)}\n")
            sections.append(f"- QC 样本数: {qc_summary.get('qc_sample_count', 0)}\n")
            sections.append(f"- QC 过滤移除代谢物: {qc_summary.get('removed_by_qc', 0)}\n")
            if median_rsd is not None:
                sections.append(f"- QC 中位 RSD: {float(median_rsd):.1%}\n")
        validation_csv = os.path.join(output_dir, "02_qc_visualization", tissue, "PLSDA_validation_summary.csv")
        if os.path.exists(validation_csv):
            validation = pd.read_csv(validation_csv)
            for _, row in validation.iterrows():
                if row.get("status") == "validated":
                    sections.append(
                        f"- PLS-DA {row['comparison']}: Q²={row.get('Q2', float('nan')):.3f}, "
                        f"CV accuracy={row.get('CV_Accuracy', float('nan')):.3f}, "
                        f"permutation P={row.get('Permutation_P', float('nan')):.3g}\n"
                    )
                else:
                    sections.append(f"- PLS-DA {row['comparison']}: {row.get('status')}\n")
        sections.append("PCA 为无监督探索；PLS-DA 仅在交叉验证和置换结果同时报告时解释。\n")

    # ── 差异代谢物 ──
    diff_csv = os.path.join(output_dir, "03_differential_analysis", tissue, "differential_analysis.csv")
    sections.append("### 3. 差异代谢物\n")

    if os.path.exists(diff_csv):
        diff_df = pd.read_csv(diff_csv)
        comp_labels = get_comp_labels(cfg.conditions, comparisons)

        for _, _, comp_name in comparisons:
            label = comp_labels.get(comp_name, comp_name)
            n_up, n_down, n_total, n_nominal = _count_sig(
                diff_df,
                comp_name,
                cfg.pval_cutoff,
                cfg.fdr_cutoff,
                cfg.log2fc_cutoff,
            )
            sections.append(f"#### {label}\n")
            sections.append(f"- FDR 显著代谢物: **{n_total}** (上调 {n_up}, 下调 {n_down})\n")
            sections.append(f"- 名义发现 (原始 P): **{n_nominal}**\n")
            vip_rule = f", VIP ≥ {cfg.vip_cutoff}" if cfg.require_vip_for_significance else ""
            sections.append(
                f"- 主阈值: FDR < {cfg.fdr_cutoff}, |Log2FC| > {cfg.log2fc_cutoff}{vip_rule}\n"
            )

            if n_total > 0:
                sections.append("\n**Top 上调代谢物:**\n")
                sections.append(_top_metabolites_table(
                    diff_df, comp_name, "up", cfg.pval_cutoff,
                    cfg.fdr_cutoff, cfg.log2fc_cutoff,
                ))
                sections.append("\n**Top 下调代谢物:**\n")
                sections.append(_top_metabolites_table(
                    diff_df, comp_name, "down", cfg.pval_cutoff,
                    cfg.fdr_cutoff, cfg.log2fc_cutoff,
                ))
            sections.append("")

        # 火山图
        for _, _, comp_name in comparisons:
            vol_png = os.path.join(output_dir, "03_differential_analysis", tissue, f"Volcano_{comp_name}.png")
            if os.path.exists(vol_png):
                sections.append(f"![火山图: {comp_name}]({vol_png})\n")
    else:
        sections.append("差异分析结果未找到。\n")

    # ── 通路富集 ──
    sections.append("### 4. KEGG 通路富集 (ORA)\n")
    for _, _, comp_name in comparisons:
        pathway_csv = os.path.join(output_dir, "05_pathway_enrichment", tissue, f"Pathway_{comp_name}.csv")
        if os.path.exists(pathway_csv):
            pw_df = pd.read_csv(pathway_csv)
            significance_col = "FDR" if "FDR" in pw_df.columns else "P-value"
            sig_pathways = pw_df[pw_df[significance_col] < cfg.fdr_cutoff]
            sections.append(f"#### {comp_name.replace('_', ' ')}\n")
            sections.append(f"- 显著通路 ({significance_col}<{cfg.fdr_cutoff}): {len(sig_pathways)}\n")
            if not sig_pathways.empty:
                for _, r in sig_pathways.head(10).iterrows():
                    sections.append(
                        f"  - {r['Name']} (P={r['P-value']:.2e}, "
                        f"FDR={r.get('FDR', float('nan')):.2e}, SigCount={r['SigCount']})\n"
                    )
            sections.append("")

    # ── GSEA ──
    sections.append("### 5. GSEA 富集分析\n")
    for _, _, comp_name in comparisons:
        gsea_csv = os.path.join(output_dir, "06_gsea", tissue, f"gsea_{comp_name}_results.csv")
        if os.path.exists(gsea_csv):
            gsea_df = pd.read_csv(gsea_csv)
            fdr_col = "FDR q-val" if "FDR q-val" in gsea_df.columns else None
            if fdr_col:
                sig = gsea_df[gsea_df[fdr_col] < 0.25]
                sections.append(f"#### {comp_name.replace('_', ' ')}\n")
                sections.append(f"- 显著通路 (FDR<0.25): {len(sig)}\n")
                if not sig.empty:
                    top_pos = sig.nlargest(min(5, len(sig)), "NES")
                    top_neg = sig.nsmallest(min(5, len(sig)), "NES")
                    sections.append("\n**Top 正向富集通路:**\n")
                    for _, r in top_pos.iterrows():
                        sections.append(f"  - {r['Term']} (NES={r['NES']:.2f}, FDR={r[fdr_col]:.2e})\n")
                    if not top_neg.empty:
                        sections.append("\n**Top 负向富集通路:**\n")
                        for _, r in top_neg.iterrows():
                            sections.append(f"  - {r['Term']} (NES={r['NES']:.2f}, FDR={r[fdr_col]:.2e})\n")
            sections.append("")

    return "\n".join(sections)


def generate_html_report(
    cfg: MetabolomicsConfig,
    output_dir: str,
) -> str:
    """
    生成汇总 HTML 报告（含嵌入图片）。

    Returns
    -------
    str — HTML 报告文件路径
    """
    print("\n  生成解读报告...")

    import yaml
    groups_yaml = yaml.safe_load(cfg.groups_path.read_text(encoding="utf-8"))
    contrast_groups = groups_yaml.get("contrast_groups", None)
    comparisons = generate_comparisons(cfg.conditions, contrast_groups)
    tissue_labels = get_tissue_labels(cfg)
    comp_labels = get_comp_labels(cfg.conditions, comparisons)

    html_sections = []

    for tissue in cfg.tissue_data_paths:
        label = tissue_labels.get(tissue, tissue)
        md = generate_tissue_summary(tissue, label, cfg, str(output_dir), comparisons)
        # 简单 Markdown → HTML 转换
        html = _md_to_html(md, output_dir, tissue)
        html_sections.append(html)

    # 组装完整 HTML
    comp_names = [cl for _, _, cn in comparisons for cl in [comp_labels.get(cn, cn)]]
    title = f"代谢组学分析报告 — {cfg.experiment_id}"

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
body {{ font-family: 'Microsoft YaHei', 'SimHei', Arial, sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; background: #fafafa; }}
h1 {{ color: #333; border-bottom: 2px solid #4C72B0; padding-bottom: 10px; }}
h2 {{ color: #4C72B0; margin-top: 40px; border-left: 4px solid #4C72B0; padding-left: 10px; }}
h3 {{ color: #555; }}
h4 {{ color: #666; }}
table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
th {{ background: #4C72B0; color: white; }}
tr:nth-child(even) {{ background: #f9f9f9; }}
code {{ background: #f0f0f0; padding: 2px 6px; border-radius: 3px; }}
img {{ display: block; margin: 10px auto; }}
.tissue-section {{ background: white; padding: 20px; margin: 20px 0; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
</style>
</head>
<body>
<h1>{title}</h1>
<p><strong>比较组:</strong> {', '.join(comp_names)}</p>
<p><strong>组织:</strong> {', '.join(tissue_labels.values())}</p>
<hr>
{''.join(f'<div class="tissue-section">{s}</div>' for s in html_sections)}
<hr>
<p style="color:#999;text-align:center;">由 metabolomics pipeline 自动生成</p>
</body>
</html>"""

    report_dir = os.path.join(output_dir, "07_interpretation_report")
    os.makedirs(report_dir, exist_ok=True)
    html_path = os.path.join(report_dir, f"metabolomics_report_{cfg.experiment_id}.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)

    # 同时生成每个组织的 Markdown 摘要
    per_tissue_dir = os.path.join(report_dir, "per_tissue")
    os.makedirs(per_tissue_dir, exist_ok=True)
    for tissue in cfg.tissue_data_paths:
        label = tissue_labels.get(tissue, tissue)
        md = generate_tissue_summary(tissue, label, cfg, str(output_dir), comparisons)
        md_path = os.path.join(per_tissue_dir, f"{tissue}_summary.md")
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(md)

    print(f"  HTML 报告: {html_path}")
    print(f"  Markdown 摘要: {per_tissue_dir}/")
    return html_path


def _md_to_html(md: str, output_dir: str | Path, tissue: str) -> str:
    """简易 Markdown → HTML（支持标题、列表、表格、图片）。"""
    lines = md.split("\n")
    html_lines = []
    in_table = False
    in_code = False

    for line in lines:
        stripped = line.strip()

        # 代码块
        if stripped == "```":
            if in_code:
                html_lines.append("</pre>")
                in_code = False
            else:
                html_lines.append('<pre style="background:#f5f5f5;padding:10px;border-radius:4px;overflow-x:auto;">')
                in_code = True
            continue
        if in_code:
            html_lines.append(line)
            continue

        # 标题
        if stripped.startswith("#### "):
            if in_table:
                html_lines.append("</table>")
                in_table = False
            html_lines.append(f"<h4>{stripped[5:]}</h4>")
        elif stripped.startswith("### "):
            if in_table:
                html_lines.append("</table>")
                in_table = False
            html_lines.append(f"<h3>{stripped[4:]}</h3>")
        elif stripped.startswith("## "):
            if in_table:
                html_lines.append("</table>")
                in_table = False
            html_lines.append(f"<h2>{stripped[3:]}</h2>")
        # 图片引用 ![alt](path)
        elif stripped.startswith("!["):
            import re
            m = re.match(r"!\[(.+?)\]\((.+?)\)", stripped)
            if m:
                img_path = m.group(2)
                html_lines.append(_embed_image(img_path))
        # 列表项
        elif stripped.startswith("- "):
            if in_table:
                html_lines.append("</table>")
                in_table = False
            html_lines.append(f"<li>{stripped[2:]}</li>")
        # 空行
        elif not stripped:
            if in_table:
                html_lines.append("</table>")
                in_table = False
            html_lines.append("<br>")
        else:
            html_lines.append(f"<p>{stripped}</p>")

    if in_table:
        html_lines.append("</table>")
    if in_code:
        html_lines.append("</pre>")

    return "\n".join(html_lines)
