#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFA 解读报告生成器

从已有的 MFA 分析结果（01_df_results/、02_result_plot/）中提取数据，
按组织（tissue）生成独立的 Markdown 解读报告。

仿照转录组学的 generate_interpretation_report.py 模式：
  - 每个组织一个 .md 文件
  - index.md 作为目录页（含分析概况与文件清单）
  - 报告包含统计摘要、关键发现、解读模板

Usage:
    # 与 run_mfa_analysis.py 集成（Step 9，自动调用）
    python run_mfa_analysis.py --experiment <id>

    # 单独运行本脚本
    python generate_mfa_interpretation_report.py --experiment <id>

Run:
    python -m omicsweave.mfa.generate_mfa_interpretation_report --experiment 20250425_cimr_wt_mice_dio
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from .mfa_experiment_config import load_mfa_config, generate_comparisons

# ── 中英文组织名映射 ───────────────────────────────────────────────────────

TISSUE_CN: dict[str, str] = {
    "liver": "肝脏（Liver）",
    "heart": "心脏（Heart）",
    "kidney": "肾脏（Kidney）",
    "brain": "脑（Brain）",
    "muscle": "骨骼肌（Skeletal Muscle）",
    "gastroc_muscle": "腓肠肌（Gastrocnemius Muscle）",
    "white_adipose": "白色脂肪（White Adipose）",
    "sub_q_fat": "皮下脂肪（Subcutaneous Fat）",
    "serum": "血清（Serum）",
    "brown_adipose": "棕色脂肪（Brown Adipose）",
    "bat": "棕色脂肪（BAT）",
    "wat": "白色脂肪（WAT）",
    "sq_adipose": "皮下脂肪（Subcutaneous Adipose）",
}

# ── 每个组织的关注代谢通路（用于解读模板中快速定位） ─────────────────────────

TISSUE_PATHWAY_FOCUS: dict[str, dict[str, list[str]]] = {
    "liver": {
        "标签": ["糖异生", "脂质合成", "酮体生成", "TCA 循环"],
        "通路": ["Gluconeogenesis", "Lipid_Synthesis", "Ketogenesis", "TCA_Cycle"],
    },
    "heart": {
        "标签": ["脂肪酸氧化", "TCA 循环", "线粒体转运"],
        "通路": ["TCA_Cycle", "Amino_Acid", "Transport"],
    },
    "kidney": {
        "标签": ["糖异生", "TCA 循环", "氨基酸代谢"],
        "通路": ["Gluconeogenesis", "TCA_Cycle", "Amino_Acid"],
    },
    "brain": {
        "标签": ["糖酵解", "TCA 循环", "氨基酸代谢"],
        "通路": ["Glycolysis", "TCA_Cycle", "Amino_Acid"],
    },
    "muscle": {
        "标签": ["糖酵解", "TCA 循环", "丙酮酸代谢"],
        "通路": ["Glycolysis", "TCA_Cycle", "Pyruvate_Metabolism"],
    },
    "gastroc_muscle": {
        "标签": ["糖酵解", "TCA 循环", "丙酮酸代谢"],
        "通路": ["Glycolysis", "TCA_Cycle", "Pyruvate_Metabolism"],
    },
    "serum": {
        "标签": ["酮体生成", "脂质合成", "甘油代谢"],
        "通路": ["Ketogenesis", "Lipid_Synthesis", "Glycerol_Metabolism"],
    },
    "white_adipose": {
        "标签": ["脂质合成", "甘油代谢", "糖酵解"],
        "通路": ["Lipid_Synthesis", "Glycerol_Metabolism", "Glycolysis"],
    },
    "sub_q_fat": {
        "标签": ["脂质合成", "甘油代谢", "糖酵解"],
        "通路": ["Lipid_Synthesis", "Glycerol_Metabolism", "Glycolysis"],
    },
}


# ── 数据类 ─────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Comp:
    """一对比较（A vs B）"""
    key: str   # 文件名中的 key，如 HFD_vs_LFD
    a: str     # 条件 A（原始 key）
    b: str     # 条件 B（原始 key）
    a_display: str  # 条件 A 显示名
    b_display: str  # 条件 B 显示名

    @property
    def label(self) -> str:
        return f"{self.a_display} vs {self.b_display}"


# ── 工具函数 ───────────────────────────────────────────────────────────────

def _safe_read_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, low_memory=False)
    except Exception:
        return None


def _read_cross_tissue_csv(results_dir: Path, category: str, filename: str) -> pd.DataFrame | None:
    """
    读取跨组织 CSV：优先新路径 01_df_results/all_tissues/<category>/，
    回退旧路径 01_df_results/<category>/。
    """
    new_path = results_dir / "01_df_results" / "all_tissues" / category / filename
    df = _safe_read_csv(new_path)
    if df is not None:
        return df
    old_path = results_dir / "01_df_results" / category / filename
    return _safe_read_csv(old_path)


def _discover_comps_from_files(results_dir: Path, cond_display: dict[str, str]) -> list[Comp]:
    """
    从各组织子目录的 differential_flux/ CSV 文件名中自动发现比较对。
    文件名格式：differential_flux_{A}_vs_{B}.csv
    """
    diff_root = results_dir / "01_df_results"
    if not diff_root.exists():
        return []

    comps: list[Comp] = []
    seen = set()
    for csv_file in sorted(diff_root.glob("**/differential_flux_*_vs_*.csv")):
        name = csv_file.stem
        comp_key = name.replace("differential_flux_", "", 1)
        if comp_key in seen:
            continue
        seen.add(comp_key)

        parts = comp_key.split("_vs_", 1)
        if len(parts) != 2:
            continue

        a_upper, b_upper = parts[0], parts[1]
        a_key = a_upper.lower()
        b_key = b_upper.lower()

        comps.append(Comp(
            key=comp_key,
            a=a_key, b=b_key,
            a_display=cond_display.get(a_key, a_upper),
            b_display=cond_display.get(b_key, b_upper),
        ))

    return comps


def _mk_table(df: pd.DataFrame, max_rows: int = 15) -> str:
    """将 DataFrame 转为 Markdown 管道表格。"""
    if df is None or df.empty:
        return "_（无）_"
    d = df.copy()
    if len(d) > max_rows:
        d = d.head(max_rows)
    d = d.fillna("")
    cols = list(d.columns)
    rows = [[str(v) for v in r] for r in d.to_numpy().tolist()]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    body = "\n".join("| " + " | ".join(r) + " |" for r in rows) if rows else ""
    return "\n".join([header, sep, body]).rstrip()


def _fmt_float(v: Any, digits: int = 3) -> str:
    """格式化浮点数。"""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "—"
    try:
        return f"{float(v):.{digits}f}"
    except (ValueError, TypeError):
        return str(v)


def _fmt_pvalue(v: Any) -> str:
    """格式化 p 值（科学计数法）。"""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "—"
    try:
        fv = float(v)
        if fv < 0.001:
            return f"{fv:.2e}"
        return f"{fv:.4f}"
    except (ValueError, TypeError):
        return str(v)


# ── 数据加载函数 ───────────────────────────────────────────────────────────

def _discover_tissues(results_dir: Path) -> list[str]:
    """从跨组织 overview 或 PCA scores 中发现组织列表。"""
    overview = _read_cross_tissue_csv(results_dir, "overview", "tissue_flux_overview.csv")
    if overview is not None and "tissue" in overview.columns:
        return sorted(overview["tissue"].unique().tolist())

    pca = _read_cross_tissue_csv(results_dir, "pca", "flux_pca_scores.csv")
    if pca is not None and "tissue" in pca.columns:
        return sorted(pca["tissue"].unique().tolist())

    return []


def _load_differential_flux(results_dir: Path, comp: Comp) -> pd.DataFrame | None:
    """加载差异通量结果（合并各组织子目录的 CSV）。"""
    df_root = results_dir / "01_df_results"
    all_frames = []
    for path in sorted(df_root.glob(f"**/differential_flux/differential_flux_{comp.key}.csv")):
        frame = _safe_read_csv(path)
        if frame is not None:
            all_frames.append(frame)
    if not all_frames:
        return None
    return pd.concat(all_frames, ignore_index=True)


def _load_pathway_flux(results_dir: Path) -> pd.DataFrame | None:
    """加载通路通量数据（合并各组织子目录的 CSV）。"""
    df_root = results_dir / "01_df_results"
    all_frames = []
    for path in sorted(df_root.glob("**/pathway_flux/pathway_flux_samples.csv")):
        frame = _safe_read_csv(path)
        if frame is not None:
            all_frames.append(frame)
    if not all_frames:
        return None
    return pd.concat(all_frames, ignore_index=True)


def _load_tca_cycle(results_dir: Path) -> pd.DataFrame | None:
    """加载 TCA 循环通量（优先从 pathway_flux_reaction_stats 推导，兼容旧 tca_cycle 文件）。"""
    df_root = results_dir / "01_df_results"
    all_frames: list[pd.DataFrame] = []
    seen_tissues: set[str] = set()

    # 新路径：由 pathway_flux_reaction_stats.csv 推导 TCA_Cycle
    for path in sorted(df_root.glob("*/pathway_flux/pathway_flux_reaction_stats.csv")):
        tissue = path.parent.parent.name
        if tissue == "all_tissues":
            continue
        frame = _safe_read_csv(path)
        if frame is None or frame.empty:
            continue
        if "pathway" not in frame.columns:
            continue
        frame = frame[frame["pathway"].astype(str) == "TCA_Cycle"].copy()
        if frame.empty:
            continue
        seen_tissues.add(str(tissue))
        all_frames.append(frame)

    cross_df = _read_cross_tissue_csv(results_dir, "pathway_flux", "pathway_flux_reaction_stats.csv")
    if cross_df is not None and not cross_df.empty:
        if "pathway" in cross_df.columns:
            cross_df = cross_df[cross_df["pathway"].astype(str) == "TCA_Cycle"].copy()
        if not cross_df.empty and ("tissue" not in cross_df.columns or not seen_tissues):
            all_frames.append(cross_df)
        elif not cross_df.empty:
            for tissue in sorted(cross_df["tissue"].dropna().astype(str).unique()):
                if tissue in seen_tissues:
                    continue
                tissue_df = cross_df[cross_df["tissue"].astype(str) == tissue]
                if not tissue_df.empty:
                    all_frames.append(tissue_df.copy())

    # 兼容旧路径：历史结果可能仍只有 tca_cycle_fluxes.csv
    if not all_frames:
        for path in sorted(df_root.glob("*/tca_cycle/tca_cycle_fluxes.csv")):
            frame = _safe_read_csv(path)
            if frame is None or frame.empty:
                continue
            tissue = path.parent.parent.name
            seen_tissues.add(str(tissue))
            all_frames.append(frame)

        legacy_cross_df = _read_cross_tissue_csv(results_dir, "tca_cycle", "tca_cycle_fluxes.csv")
        if legacy_cross_df is not None and not legacy_cross_df.empty:
            if "tissue" not in legacy_cross_df.columns or not seen_tissues:
                all_frames.append(legacy_cross_df)
            else:
                for tissue in sorted(legacy_cross_df["tissue"].dropna().astype(str).unique()):
                    if tissue in seen_tissues:
                        continue
                    tissue_df = legacy_cross_df[legacy_cross_df["tissue"].astype(str) == tissue]
                    if not tissue_df.empty:
                        all_frames.append(tissue_df.copy())

    if not all_frames:
        return None
    merged = pd.concat(all_frames, ignore_index=True)
    keep_cols = [c for c in ["tissue", "reaction", "condition", "mean", "std", "n"] if c in merged.columns]
    return merged[keep_cols] if keep_cols else merged


def _load_overview(results_dir: Path) -> pd.DataFrame | None:
    """加载组织通量概览。"""
    return _read_cross_tissue_csv(results_dir, "overview", "tissue_flux_overview.csv")


def _load_mid_labeling(results_dir: Path) -> pd.DataFrame | None:
    """加载 MID 总标记度。"""
    return _read_cross_tissue_csv(results_dir, "mid", "mid_total_labeling.csv")


def _load_pca_scores(results_dir: Path) -> pd.DataFrame | None:
    """加载 PCA 得分。"""
    return _read_cross_tissue_csv(results_dir, "pca", "flux_pca_scores.csv")


# ── 摘要计算函数 ───────────────────────────────────────────────────────────

def _summarize_diff_flux(diff_df: pd.DataFrame, tissue: str) -> dict[str, Any]:
    """计算某个组织在差异通量分析中的摘要统计。"""
    sub = diff_df[diff_df["tissue"] == tissue].copy() if "tissue" in diff_df.columns else diff_df.copy()

    if sub.empty:
        return {"total": 0, "sig_p": 0, "sig_fdr": 0, "top_up": pd.DataFrame(), "top_down": pd.DataFrame()}

    # 确保数值列
    for c in ["log2FC", "pvalue_mannwhitney", "fdr_mannwhitney", "diff", "cohen_d"]:
        if c in sub.columns:
            sub[c] = pd.to_numeric(sub[c], errors="coerce")

    total = len(sub)
    sig_p = sub[sub["pvalue_mannwhitney"] < 0.05] if "pvalue_mannwhitney" in sub.columns else pd.DataFrame()
    sig_fdr = sub[sub["fdr_mannwhitney"] < 0.1] if "fdr_mannwhitney" in sub.columns else pd.DataFrame()

    # 按 diff（算术差）排序：diff > 0 表示 A 通量比 B 高
    def _top_fluxes(df: pd.DataFrame, ascending_diff: bool, n: int = 15) -> pd.DataFrame:
        if df.empty or "diff" not in df.columns:
            return pd.DataFrame()
        s = df.dropna(subset=["diff"]).sort_values("diff", ascending=ascending_diff).head(n)
        cols = [c for c in ["reaction", "diff", "log2FC", "pvalue_mannwhitney", "fdr_mannwhitney", "cohen_d"] if c in s.columns]
        return s[cols].copy()

    return {
        "total": total,
        "sig_p": len(sig_p),
        "sig_fdr": len(sig_fdr),
        "top_up": _top_fluxes(sig_p, ascending_diff=False),    # diff 降序：最大正值在前（A > B）
        "top_down": _top_fluxes(sig_p, ascending_diff=True),   # diff 升序：最大负值在前（A < B）
    }


def _summarize_tca(tca_df: pd.DataFrame, tissue: str) -> pd.DataFrame | None:
    """提取某组织的 TCA 循环通量摘要。"""
    if "tissue" not in tca_df.columns:
        return None
    sub = tca_df[tca_df["tissue"] == tissue].copy()
    if sub.empty:
        return None
    cols = [c for c in ["tissue", "reaction", "condition", "mean", "std", "n"] if c in sub.columns]
    return sub[cols]


def _summarize_pathway(pathway_df: pd.DataFrame, tissue: str, conditions: list[str]) -> pd.DataFrame | None:
    """提取某组织各通路的条件间均值对比。"""
    if "tissue" not in pathway_df.columns:
        return None
    sub = pathway_df[pathway_df["tissue"] == tissue].copy()
    if sub.empty:
        return None

    # 按通路 × 条件 → 均值
    if all(c in sub.columns for c in ["pathway", "condition", "flux"]):
        summary = sub.groupby(["pathway", "condition"])["flux"].mean().reset_index()
        summary = summary.rename(columns={"flux": "mean_flux"})
        summary["mean_flux"] = summary["mean_flux"].round(2)
        # 透视表：行=通路，列=条件
        pivot = summary.pivot(index="pathway", columns="condition", values="mean_flux").reset_index()
        return pivot
    return None


def _summarize_mid(mid_df: pd.DataFrame, tissue: str) -> pd.DataFrame | None:
    """提取某组织的 MID 标记度摘要。"""
    if mid_df is None or mid_df.empty:
        return None
    if "metabolite" not in mid_df.columns:
        return None

    # 过滤到该组织的代谢物
    sub = mid_df[mid_df["metabolite"].str.startswith(tissue + "--")].copy()
    if sub.empty:
        return None

    # 去掉组织前缀以简化显示
    sub["metabolite_short"] = sub["metabolite"].str.replace(f"{tissue}--", "", regex=False)

    # 透视表：行=代谢物，列=条件
    if "condition" in sub.columns and "total_labeling" in sub.columns:
        pivot = sub.pivot(index="metabolite_short", columns="condition", values="total_labeling").reset_index()
        # 四舍五入
        for c in pivot.columns:
            if c != "metabolite_short":
                pivot[c] = pd.to_numeric(pivot[c], errors="coerce").round(4)
        # 排序：按第一个条件列的标记度降序
        cond_col = [c for c in pivot.columns if c != "metabolite_short"]
        if cond_col:
            pivot = pivot.sort_values(cond_col[0], ascending=False)
        return pivot
    return None


def _summarize_overview(overview_df: pd.DataFrame, tissue: str) -> dict[str, Any] | None:
    """提取某组织的通量概览。"""
    if overview_df is None or overview_df.empty:
        return None
    sub = overview_df[overview_df["tissue"] == tissue].copy()
    if sub.empty:
        return None
    return sub.to_dict("records")


def _summarize_pca(pca_df: pd.DataFrame, tissue: str) -> str:
    """简要描述某组织的 PCA 结果。"""
    if pca_df is None or pca_df.empty:
        return "—"
    sub = pca_df[pca_df["tissue"] == tissue].copy()
    if sub.empty:
        return "—"

    pc_cols = [c for c in sub.columns if c.startswith("PC")]
    if not pc_cols:
        return "—"

    # 检查 condition 列
    if "condition" not in sub.columns:
        return "—"

    # 各条件在 PC1 上的均值
    pc1 = pc_cols[0]
    means = sub.groupby("condition")[pc1].mean()
    parts = []
    for cond, val in means.items():
        parts.append(f"{cond}: {_fmt_float(val)}")
    return f"PC1 均值 — " + "; ".join(parts)


# ── 报告页面生成 ───────────────────────────────────────────────────────────

def tissue_page(
    tissue: str,
    comps: list[Comp],
    results_dir: Path,
    cond_display: dict[str, str],
    pathway_groups: dict[str, list[str]] | None,
) -> str:
    """生成单个组织的解读报告页面。"""
    title = TISSUE_CN.get(tissue, tissue)
    lines: list[str] = []

    lines.append(f"## {title} — MFA 解读报告")
    lines.append("")
    lines.append(f"- **组织**：`{tissue}`")
    lines.append(f"- **生成时间**：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    # 组别定义
    lines.append("### 组别定义")
    lines.append("")
    for code, display in sorted(cond_display.items()):
        lines.append(f"- **{display}**（内部代码: `{code}`）")
    lines.append("")

    # 阅读说明
    lines.append("### 结果阅读说明")
    lines.append("")
    lines.append("- **差异通量**：`p < 0.05`（Mann-Whitney U）认为初步显著；`FDR < 0.1`（BH 校正）为严格显著。")
    lines.append("- **log2FC**：反映条件间通量的倍数变化方向与幅度。")
    lines.append("- **Cohen's d**：效应量，|d| > 0.8 通常认为大效应。")
    lines.append("- **通路通量**：各代谢通路内反应的汇总对比。")
    lines.append("- **MID 标记度**：13C-葡萄糖碳原子流入各代谢物的比例。值越高表示标记底物贡献越大。")
    lines.append('- **注意**：n < 3 的组别，统计功效较弱，相关结论应标注为"探索性"。')
    lines.append("")

    # ── 组织通量概览 ──
    overview_df = _load_overview(results_dir)
    overview_info = _summarize_overview(overview_df, tissue)
    if overview_info:
        lines.append("### 组织通量概览")
        lines.append("")
        ov_df = pd.DataFrame(overview_info)
        display_cols = [c for c in ["tissue", "condition", "n_reactions", "total_abs_flux", "mean_abs_flux"] if c in ov_df.columns]
        if "condition" in ov_df.columns:
            ov_df["condition"] = ov_df["condition"].map(cond_display).fillna(ov_df["condition"])
        lines.append(_mk_table(ov_df[display_cols]))
        lines.append("")

    # ── PCA 概况 ──
    pca_df = _load_pca_scores(results_dir)
    pca_summary = _summarize_pca(pca_df, tissue)
    lines.append("### PCA 降维概况")
    lines.append("")
    lines.append(pca_summary)
    lines.append("")
    pca_plot = results_dir / "02_result_plot" / "all_tissues" / "3_pca_clustering" / f"flux_pca_{tissue}.png"
    if pca_plot.exists():
        lines.append(f"![PCA 散点图](../02_result_plot/all_tissues/3_pca_clustering/flux_pca_{tissue}.png)")
        lines.append("")

    # ── 逐比较差异分析 ──
    tca_df = _load_tca_cycle(results_dir)
    pathway_df = _load_pathway_flux(results_dir)
    mid_df = _load_mid_labeling(results_dir)
    conditions = list(cond_display.keys())

    for comp in comps:
        lines.append(f"### 比较：{comp.label}")
        lines.append("")

        # ── 差异通量 ──
        diff_df = _load_differential_flux(results_dir, comp)
        if diff_df is None:
            lines.append("- **差异通量分析**：未找到对应结果文件。")
            lines.append("")
        else:
            s = _summarize_diff_flux(diff_df, tissue)
            lines.append(f"- **差异通量总结**：通量总数 {s['total']}；显著（p<0.05）{s['sig_p']}；严格显著（FDR<0.1）{s['sig_fdr']}。")
            lines.append("")

            if not s["top_up"].empty:
                lines.append("#### 通量升高（Top，diff > 0，p<0.05）")
                lines.append("")
                lines.append(_mk_table(s["top_up"]))
                lines.append("")

            if not s["top_down"].empty:
                lines.append("#### 通量降低（Top，diff < 0，p<0.05）")
                lines.append("")
                lines.append(_mk_table(s["top_down"]))
                lines.append("")

            # 火山图
            volcano_plot = results_dir / "02_result_plot" / tissue / "1_volcano" / f"volcano_{comp.key}.png"
            if volcano_plot.exists():
                lines.append(f"![火山图](../02_result_plot/{tissue}/1_volcano/volcano_{comp.key}.png)")
                lines.append("")

        # ── TCA 循环 ──
        tca_sub = _summarize_tca(tca_df, tissue) if tca_df is not None else None
        if tca_sub is not None and not tca_sub.empty:
            if "condition" in tca_sub.columns:
                tca_sub["condition"] = tca_sub["condition"].map(cond_display).fillna(tca_sub["condition"])
            lines.append("#### TCA 循环通量")
            lines.append("")
            lines.append(_mk_table(tca_sub))
            lines.append("")

        # ── 通路通量 ──
        pathway_sub = _summarize_pathway(pathway_df, tissue, conditions) if pathway_df is not None else None
        if pathway_sub is not None and not pathway_sub.empty:
            # 将 condition 列名替换为显示名
            rename_map = {c: cond_display.get(c, c) for c in pathway_sub.columns if c != "pathway"}
            pathway_sub = pathway_sub.rename(columns=rename_map)
            lines.append("#### 代谢通路通量对比（各条件均值）")
            lines.append("")
            lines.append(_mk_table(pathway_sub))
            lines.append("")
            # 通路箱线图
            pathway_plot = results_dir / "02_result_plot" / tissue / "4_pathway_flux" / f"pathway_flux_{tissue}.png"
            if pathway_plot.exists():
                lines.append(f"![通路箱线图](../02_result_plot/{tissue}/4_pathway_flux/pathway_flux_{tissue}.png)")
                lines.append("")

        # ── 解读模板 ──
        lines.append("#### 证据优先解读模板（可复制到汇报 PPT）")
        lines.append("")
        focus = TISSUE_PATHWAY_FOCUS.get(tissue, {})
        focus_tags = focus.get("标签", [])
        focus_str = "、".join(focus_tags) if focus_tags else "主要代谢通路"
        lines.append(f"> **一句话结论（建议写法）**：在该组织({title})中，`{comp.a_display} vs {comp.b_display}` 引起的主要变化可概括为「{focus_str}相关通量的增强/抑制」，证据来自（Top 差异通量 + 通路对比）。")
        lines.append("")
        lines.append("- **最强证据 3 条（建议写法）**：")
        lines.append("  - 差异通量：列出 1–2 个 log2FC 绝对值最大且 p 值最小的通量（写方向与倍数）。")
        lines.append("  - 通路：列出 1–2 个条件间变化最大的通路名称及其均值。")
        lines.append("  - 若有：TCA 循环中变化最大的反应步骤。")
        lines.append("")

    # ── MID 标记分析 ──
    mid_sub = _summarize_mid(mid_df, tissue)
    if mid_sub is not None and not mid_sub.empty:
        lines.append("### MID 同位素标记分析")
        lines.append("")
        rename_map = {c: cond_display.get(c, c) for c in mid_sub.columns if c != "metabolite_short"}
        mid_sub = mid_sub.rename(columns=rename_map)
        lines.append("_总标记度 = 1 − M+0，值越高表示 13C-葡萄糖碳流入越多。_")
        lines.append("")
        lines.append(_mk_table(mid_sub))
        lines.append("")

    # ── 下一步建议 ──
    lines.append("### 下一步建议")
    lines.append("")
    lines.append('- **统计功效**：优先检查每组重复数；若某组 n < 3，该结论应标注为"探索性"。')
    lines.append("- **验证路径**：")
    lines.append("  - 与转录组学数据交叉验证（同组织差异基因 vs 差异通量）。")
    lines.append("  - 靶向代谢物定量（LC-MS）验证 MID 标记度变化。")
    lines.append("  - 13C 示踪实验重复或扩展（更多时间点 / 组织）。")
    lines.append("- **建议补充图**：按组织做「多比较并排」的通路通量柱状图。")
    lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def generate_index(
    tissues: list[str],
    comps: list[Comp],
    results_dir: Path,
    cond_display: dict[str, str],
    experiment_id: str,
) -> str:
    """生成 index.md 目录页。包含实验概况 + 文件清单 + 各组织报告链接。"""
    lines: list[str] = []

    lines.append("# MFA 解读报告（每组织 1 页）")
    lines.append("")
    lines.append(f"- **实验 ID**：`{experiment_id}`")
    lines.append(f"- **生成时间**：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    cond_names = " / ".join(cond_display[c] for c in sorted(cond_display.keys()))
    lines.append(f"- **组别命名**：{cond_names}")
    lines.append(f"- **组织数**：{len(tissues)}")
    if comps:
        lines.append(f"- **比较对**：{len(comps)}")
    lines.append("")

    # 组织目录
    lines.append("### 组织目录")
    lines.append("")
    for tissue in tissues:
        cn = TISSUE_CN.get(tissue, tissue)
        lines.append(f"- [{cn}](./{tissue}.md)")
    lines.append("")

    # 数据表清单
    lines.append("### 分析数据表（01_df_results/）")
    lines.append("")
    df_dir = results_dir / "01_df_results"
    if df_dir.exists():
        for csv_file in sorted(df_dir.glob("**/*.csv")):
            rel = csv_file.relative_to(df_dir)
            lines.append(f"- `{rel}`")
        lines.append("")

    # 图表清单
    lines.append("### 分析图表（02_result_plot/）")
    lines.append("")
    plots_dir = results_dir / "02_result_plot"
    if plots_dir.exists():
        for png_file in sorted(plots_dir.glob("**/*.png")):
            rel = png_file.relative_to(plots_dir)
            lines.append(f"- `{rel}`")
        lines.append("")

    lines.append("### 重要备注")
    lines.append("")
    lines.append('- 本报告为「**证据汇总 + 解读模板**」，优先展示可复核的统计表格。')
    lines.append('- 重复数偏低的组别/组织，相关结论请标注为「探索性」。')
    lines.append("- MFA 差异分析基于全部拟合迭代（非仅 best-fit），使用 Mann-Whitney U 检验 + BH-FDR 校正。")
    lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ── 主流程 ─────────────────────────────────────────────────────────────────

def generate_interpretation_report(
    *,
    experiment_id: str | None = None,
    config_path: str | None = None,
) -> Path:
    """
    生成 MFA 解读报告。

    返回输出目录路径。
    """
    # 加载配置
    cfg = load_mfa_config(experiment_id=experiment_id, config_path=config_path)
    results_dir = cfg.output_dir

    if not results_dir.exists():
        raise FileNotFoundError(f"分析结果目录不存在: {results_dir}")

    # 条件显示名映射
    cond_display: dict[str, str] = {}
    for code, info in cfg.conditions.items():
        if isinstance(info, dict):
            cond_display[code.lower()] = info.get("display", code)
        else:
            cond_display[code.lower()] = str(code)

    # 发现组织
    tissues = _discover_tissues(results_dir)
    if not tissues:
        print("[WARN] 未发现任何组织数据，跳过报告生成。")
        return results_dir

    # 发现比较对（从已有的差异通量文件自动扫描）
    comps = _discover_comps_from_files(results_dir, cond_display)
    if not comps:
        # 回退：从 config 的 comparisons 生成
        comp_tuples = cfg.comparisons or generate_comparisons(cfg.conditions)
        for a, b in comp_tuples:
            a_key = a.lower()
            b_key = b.lower()
            file_key = f"{a.upper()}_vs_{b.upper()}"
            comps.append(Comp(
                key=file_key,
                a=a_key, b=b_key,
                a_display=cond_display.get(a_key, a_key),
                b_display=cond_display.get(b_key, b_key),
            ))

    # 创建输出目录
    out_dir = results_dir / "03_interpretation_report"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 生成每个组织的报告
    created: list[str] = []
    for tissue in tissues:
        page = tissue_page(
            tissue=tissue,
            comps=comps,
            results_dir=results_dir,
            cond_display=cond_display,
            pathway_groups=cfg.pathway_groups,
        )
        out_path = out_dir / f"{tissue}.md"
        out_path.write_text(page, encoding="utf-8")
        created.append(out_path.name)

    # 生成 index.md（替代旧的 mfa_report.md）
    index = generate_index(
        tissues=tissues,
        comps=comps,
        results_dir=results_dir,
        cond_display=cond_display,
        experiment_id=cfg.experiment_id,
    )
    (out_dir / "index.md").write_text(index, encoding="utf-8")

    # 如果旧的 mfa_report.md 存在，删除它（内容已合并到 index.md）
    old_report = results_dir / "mfa_report.md"
    if old_report.exists():
        old_report.unlink()
        print(f"[INFO] 已移除旧报告: {old_report.name}（内容已合并至 index.md）")

    print(f"\n[OK] MFA 解读报告生成完成")
    print(f"  输出目录: {out_dir}")
    print(f"  文件列表:")
    for fn in ["index.md"] + created:
        print(f"    - {fn}")

    return out_dir


def generate_report_from_output_dir(
    results_dir: Path,
    tissues: list[str],
    comp_tuples: list[tuple[str, str]] | None = None,
    cond_display: dict[str, str] | None = None,
    experiment_id: str = "",
) -> Path:
    """
    从已有的分析结果目录直接生成解读报告（不需要 YAML 配置）。

    用于 MFAAnalyser.generate_report() 内部调用。

    Parameters
    ----------
    results_dir : Path
        分析结果根目录（包含 01_df_results/ 和 02_result_plot/）
    tissues : list[str]
        组织列表
    comp_tuples : list[tuple[str, str]] | None
        比较对列表 [(cond_a, cond_b), ...]，key 已小写。
        若为 None，则从差异通量文件名自动发现。
    cond_display : dict[str, str] | None
        条件代码 → 显示名映射
    experiment_id : str
        实验 ID（仅用于 index 页显示）

    Returns
    -------
    Path
        输出目录路径
    """
    if not results_dir.exists():
        raise FileNotFoundError(f"分析结果目录不存在: {results_dir}")

    cond_display = cond_display or {}

    if not tissues:
        tissues = _discover_tissues(results_dir)
    if not tissues:
        print("[WARN] 未发现任何组织数据，跳过报告生成。")
        return results_dir

    # 优先从文件发现比较对，回退到参数传入的比较对
    comps = _discover_comps_from_files(results_dir, cond_display)
    if not comps and comp_tuples:
        for a, b in comp_tuples:
            a_key = a.lower()
            b_key = b.lower()
            file_key = f"{a.upper()}_vs_{b.upper()}"
            comps.append(Comp(
                key=file_key,
                a=a_key, b=b_key,
                a_display=cond_display.get(a_key, a_key),
                b_display=cond_display.get(b_key, b_key),
            ))

    # 创建输出目录
    out_dir = results_dir / "03_interpretation_report"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 生成每个组织的报告
    created: list[str] = []
    for tissue in tissues:
        page = tissue_page(
            tissue=tissue,
            comps=comps,
            results_dir=results_dir,
            cond_display=cond_display,
            pathway_groups=None,
        )
        out_path = out_dir / f"{tissue}.md"
        out_path.write_text(page, encoding="utf-8")
        created.append(out_path.name)

    # 生成 index.md
    index = generate_index(
        tissues=tissues,
        comps=comps,
        results_dir=results_dir,
        cond_display=cond_display,
        experiment_id=experiment_id,
    )
    (out_dir / "index.md").write_text(index, encoding="utf-8")

    # 删除旧的 mfa_report.md
    old_report = results_dir / "mfa_report.md"
    if old_report.exists():
        old_report.unlink()
        print(f"[INFO] 已移除旧报告: {old_report.name}（内容已合并至 index.md）")

    print(f"[OK] MFA 解读报告生成完成: {out_dir}")
    print(f"  文件: index.md, {', '.join(created)}")

    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MFA 解读报告生成器（每组织 1 页）",
    )
    parser.add_argument("--experiment", default=None,
                        help="实验 ID，如 20250425_cimr_wt_mice_dio")
    parser.add_argument("--config", default=None,
                        help="配置文件路径（config.yaml）")
    args = parser.parse_args()

    if not args.experiment and not args.config:
        parser.error("请传入 --experiment <id> 或 --config <path>")

    generate_interpretation_report(
        experiment_id=args.experiment,
        config_path=args.config,
    )


if __name__ == "__main__":
    main()
