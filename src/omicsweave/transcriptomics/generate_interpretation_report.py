"""
Generate a per-organ interpretation report (1 page / organ) from existing outputs:
  - analysis_output/04_de_results (DESeq2 gene-level results)
  - analysis_output/05_result_plot/4_enrichment (Enrichr CSVs)
  - analysis_output/05_result_plot/5_gsea (GSEA prerank results CSVs)

This script focuses on:
  - consistent group naming (LFD-14W / HFD-14W / HSD-14W)
  - sample size summary per tissue
  - concise evidence tables (top DE genes, top GSEA terms, top Enrichr terms)
  - Chinese narrative guidance (non-overclaiming; evidence-first)

Run:
  python -m omicsweave.transcriptomics.generate_interpretation_report
"""

from __future__ import annotations

import os
import re
import argparse
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd

from .experiment_config import (
    load_transcriptomics_config, config_arg_help,
    resolve_comparisons, get_cond_map, get_comp_labels, get_tissue_labels,
)


# Placeholders — populated by _apply_config()
ANALYSIS_DIR: Path = Path()  # overwritten at runtime
DE_DIR: Path = Path()
ENRICH_DIR: Path = Path()
GSEA_DIR: Path = Path()
META_PATH: Path = Path()
OUT_DIR: Path = Path()

COND_MAP: dict[str, str] = {}
COMP_MAP: dict[str, tuple[str, str]] = {}
TISSUE_LABELS: dict[str, str] = {}


def _apply_config(experiment: str | None, config_path: str | None, outdir: str | None) -> dict:
    """
    Apply experiment YAML config to globals (dirs + group display names).
    Returns a small dict for logging.
    """
    cfg = load_transcriptomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )

    global ANALYSIS_DIR, DE_DIR, ENRICH_DIR, GSEA_DIR, META_PATH, OUT_DIR
    global COND_MAP, COMP_MAP, TISSUE_LABELS

    ANALYSIS_DIR = Path(str(cfg.output_dir)).resolve()
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = str(ANALYSIS_DIR)

    DE_DIR = ANALYSIS_DIR / "04_de_results"
    ENRICH_DIR = ANALYSIS_DIR / "05_result_plot" / "4_enrichment"
    GSEA_DIR = ANALYSIS_DIR / "05_result_plot" / "5_gsea"
    META_PATH = ANALYSIS_DIR / "sample_metadata.csv"
    OUT_DIR = ANALYSIS_DIR / "06_interpretation_report"

    # condition display names — from config or fall back to codes
    COND_MAP = get_cond_map(cfg.conditions or {})
    if not COND_MAP:
        COND_MAP = {k: k for k in (cfg.conditions or {})}

    # comparisons — auto-generated pairwise
    comp_tuples = resolve_comparisons(cfg)
    COMP_MAP = {}
    for a, b, name in comp_tuples:
        COMP_MAP[name] = (COND_MAP.get(a, a), COND_MAP.get(b, b))

    # tissue labels
    TISSUE_LABELS = get_tissue_labels(cfg)

    return {"experiment_id": cfg.experiment_id, "config_path": str(cfg.config_path), "output_dir": str(ANALYSIS_DIR)}

TISSUE_CN = {
    "liver": "肝脏（Liver）",
    "heart": "心脏（Heart）",
    "kidney": "肾脏（Kidney）",
    "gast_muscle": "腓肠肌（Gastrocnemius muscle）",
    "brain": "脑（Brain）",
    "sq_adipose": "皮下脂肪（Subcutaneous adipose）",
}


MARKERS: dict[str, list[str]] = {
    "liver": [
        # lipid / cholesterol / FAO
        "Srebf1",
        "Srebf2",
        "Fasn",
        "Acaca",
        "Scd1",
        "Elovl6",
        "Hmgcr",
        "Sqle",
        "Cpt1a",
        "Cpt2",
        "Ppara",
        "Pparg",
        "Ppargc1a",
        "Fgf21",
        "Cd36",
        "Fabp1",
        # glucose
        "Pck1",
        "G6pc",
        # fibrosis / inflammation
        "Tgfb1",
        "Col1a1",
        "Col3a1",
        "Acta2",
        "Il1b",
        "Tnf",
        "Ccl2",
        "Nos2",
    ],
    "heart": [
        # mito/FAO
        "Ppargc1a",
        "Cpt1b",
        "Cpt2",
        "Acadm",
        "Acadl",
        "Ucp3",
        # ETC
        "Ndufs1",
        "Cox4i1",
        "Atp5f1a",
        # contractile/remodeling
        "Myh6",
        "Myh7",
        "Tnni3",
        "Tnnt2",
        "Nppa",
        "Nppb",
        # inflammation/IFN
        "Il1b",
        "Tnf",
        "Cxcl9",
        "Stat1",
        "Isg15",
    ],
    "kidney": [
        "Havcr1",
        "Lcn2",
        "Krt8",
        "Krt18",
        "Slc22a6",
        "Slc34a1",
        "Ppara",
        "Ppargc1a",
        "Cpt1a",
        "Cpt2",
        "Il1b",
        "Tnf",
        "Ccl2",
        "Isg15",
    ],
    "gast_muscle": [
        "Ppargc1a",
        "Pdk4",
        "Cpt1b",
        "Cpt2",
        "Ucp3",
        "Slc2a4",
        "Myog",
        "Atrogin1",
        "Trim63",
        "Il1b",
        "Tnf",
        "Isg15",
    ],
    "brain": [
        "Aif1",
        "Trem2",
        "Apoe",
        "C1qa",
        "C1qb",
        "C1qc",
        "Gfap",
        "Stat1",
        "Isg15",
        "Il1b",
        "Tnf",
    ],
    "sq_adipose": [
        "Pparg",
        "Adipoq",
        "Lep",
        "Plin1",
        "Cpt1a",
        "Cpt2",
        "Ppargc1a",
        "Ucp1",
        "Tnf",
        "Il1b",
        "Ccl2",
    ],
}


METABOLIC_KWS = [
    "mitochond",
    "oxidative",
    "oxphos",
    "respir",
    "tca",
    "glycol",
    "glucone",
    "lipid",
    "fatty",
    "beta-oxid",
    "peroxisom",
    "cholesterol",
    "bile",
    "ppar",
    "insulin",
    "amino",
    "ketone",
    "ribosom",
    "proteasom",
    "er",
    "unfolded",
]

IMMUNE_KWS = [
    "interferon",
    "cytok",
    "chemok",
    "toll",
    "nf-kb",
    "inflamm",
    "complement",
    "antigen",
    "mhc",
    "leukocyte",
    "macroph",
    "t cell",
    "b cell",
]

FIBROSIS_KWS = [
    "extracellular matrix",
    "collagen",
    "fibro",
    "tgf",
    "wound",
    "epithelial",
    "mesench",
]


@dataclass(frozen=True)
class Comp:
    key: str  # e.g. DIO_vs_WT
    a: str  # group A (numerator)
    b: str  # group B (denominator)

    @property
    def label(self) -> str:
        return f"{self.a} vs {self.b}"


def _safe_read_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def _mk_table(df: pd.DataFrame, max_rows: int = 15) -> str:
    if df is None or df.empty:
        return "_（无）_"
    d = df.copy()
    if len(d) > max_rows:
        d = d.head(max_rows)
    # Avoid pandas.to_markdown() optional dependency on "tabulate"
    d = d.fillna("")
    cols = list(d.columns)
    # Convert all values to strings (keep concise)
    rows = [[str(v) for v in r] for r in d.to_numpy().tolist()]

    # Basic pipe table with header separator
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    body = "\n".join("| " + " | ".join(r) + " |" for r in rows) if rows else ""
    return "\n".join([header, sep, body]).rstrip()


def _norm_term(s: str) -> str:
    # Remove library prefix like "Reactome_Pathways_2024__"
    return re.sub(r"^[^_]+_[^_]+_\d{4}__", "", s).strip()


def _keyword_hits(terms: Iterable[str], keywords: list[str]) -> list[str]:
    out: list[str] = []
    for t in terms:
        tl = t.lower()
        if any(kw in tl for kw in keywords):
            out.append(t)
    return out


def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(META_PATH)
    meta["condition_named"] = meta["condition"].map(COND_MAP).fillna(meta["condition"])
    return meta


def sample_n_table(meta: pd.DataFrame, tissue: str) -> pd.DataFrame:
    sub = meta[meta["tissue"] == tissue].copy()
    counts = (
        sub.groupby(["condition_named"])["sample"]
        .count()
        .reset_index()
        .rename(columns={"condition_named": "组别", "sample": "n"})
        .sort_values("组别")
    )
    return counts


def load_de(tissue: str, comp_key: str) -> pd.DataFrame | None:
    path = DE_DIR / tissue / f"de_{comp_key}.csv"
    if not path.exists():
        path = DE_DIR / tissue / f"de_{comp_key}_wilcoxon.csv"
    df = _safe_read_csv(path)
    if df is None or df.empty:
        return None
    # ensure columns
    if "gene" not in df.columns:
        return None
    return df


def summarize_de(df: pd.DataFrame) -> dict[str, object]:
    # Robust numeric conversion
    for c in ["log2FoldChange", "padj", "pvalue", "baseMean"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    sig = df[df["padj"].notna() & (df["padj"] < 0.05)]
    sig_strict = sig[sig["log2FoldChange"].abs() >= 1.0]
    up = (
        sig[sig["log2FoldChange"].notna() & (sig["log2FoldChange"] > 0)]
        .sort_values(["log2FoldChange", "padj"], ascending=[False, True])
        .head(20)
    )
    down = (
        sig[sig["log2FoldChange"].notna() & (sig["log2FoldChange"] < 0)]
        .sort_values(["log2FoldChange", "padj"], ascending=[True, True])
        .head(20)
    )

    def _sel(d: pd.DataFrame) -> pd.DataFrame:
        cols = [c for c in ["gene", "log2FoldChange", "padj", "baseMean"] if c in d.columns]
        out = d[cols].copy()
        if "log2FoldChange" in out.columns:
            out["log2FoldChange"] = out["log2FoldChange"].map(lambda x: float(x) if pd.notna(x) else x)
        if "padj" in out.columns:
            out["padj"] = out["padj"].map(lambda x: float(x) if pd.notna(x) else x)
        if "baseMean" in out.columns:
            out["baseMean"] = out["baseMean"].map(lambda x: float(x) if pd.notna(x) else x)
        return out

    return {
        "n_genes": int(df.shape[0]),
        "n_sig_padj": int(sig.shape[0]),
        "n_sig_padj_lfc1": int(sig_strict.shape[0]),
        "top_up": _sel(up),
        "top_down": _sel(down),
    }


def load_gsea(tissue: str, comp_key: str) -> pd.DataFrame | None:
    path = GSEA_DIR / f"gsea_{tissue}_{comp_key}_results.csv"
    df = _safe_read_csv(path)
    if df is None or df.empty:
        return None
    # Normalize columns
    for c in ["NES", "FDR q-val", "NOM p-val"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "Term" in df.columns:
        df["Term_clean"] = df["Term"].map(_norm_term)
    else:
        df["Term_clean"] = ""
    return df


def summarize_gsea(df: pd.DataFrame) -> dict[str, object]:
    if df is None or df.empty:
        return {
            "n_terms": 0,
            "n_fdr25": 0,
            "top_pos": pd.DataFrame(),
            "top_neg": pd.DataFrame(),
            "metabolic_hits": [],
            "immune_hits": [],
            "fibrosis_hits": [],
        }
    fdr25 = df[df["FDR q-val"].notna() & (df["FDR q-val"] <= 0.25)]

    # Prefer displaying terms that pass FDR<=0.25; fall back to overall if too few.
    pos_src = fdr25[fdr25["NES"].notna() & (fdr25["NES"] > 0)]
    neg_src = fdr25[fdr25["NES"].notna() & (fdr25["NES"] < 0)]
    if len(pos_src) < 5:
        pos_src = df[df["NES"].notna() & (df["NES"] > 0)]
    if len(neg_src) < 5:
        neg_src = df[df["NES"].notna() & (df["NES"] < 0)]

    pos = pos_src.sort_values("NES", ascending=False).head(15)
    neg = neg_src.sort_values("NES", ascending=True).head(15)

    def _sel(d: pd.DataFrame) -> pd.DataFrame:
        cols = [c for c in ["Term_clean", "NES", "FDR q-val", "NOM p-val", "Lead_genes"] if c in d.columns]
        out = d[cols].copy()
        out = out.rename(
            columns={
                "Term_clean": "Term",
                "FDR q-val": "FDR",
                "NOM p-val": "NOM_p",
                "Lead_genes": "Lead_genes",
            }
        )
        return out

    terms = [t for t in df["Term_clean"].dropna().astype(str).tolist() if t.strip()]
    return {
        "n_terms": int(df.shape[0]),
        "n_fdr25": int(fdr25.shape[0]),
        "top_pos": _sel(pos),
        "top_neg": _sel(neg),
        "metabolic_hits": _keyword_hits(terms, METABOLIC_KWS)[:20],
        "immune_hits": _keyword_hits(terms, IMMUNE_KWS)[:20],
        "fibrosis_hits": _keyword_hits(terms, FIBROSIS_KWS)[:20],
    }


def load_enrichr(tissue: str, comp_key: str, direction: str) -> pd.DataFrame | None:
    # direction: "up" or "down"
    path = ENRICH_DIR / f"enrichr_{tissue}_{comp_key}_{direction}.csv"
    df = _safe_read_csv(path)
    if df is None or df.empty:
        return None
    for c in ["Adjusted P-value", "P-value", "Combined Score"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def summarize_enrichr(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    cols = [c for c in ["Gene_set", "Term", "Adjusted P-value", "Combined Score", "Overlap"] if c in df.columns]
    out = df[cols].copy()
    out = out.rename(columns={"Adjusted P-value": "adjP", "Combined Score": "score"})
    out = out.sort_values(["adjP", "score"], ascending=[True, False]).head(15)
    return out


def marker_table(de_df: pd.DataFrame | None, markers: list[str]) -> pd.DataFrame:
    if de_df is None or de_df.empty:
        return pd.DataFrame()
    df = de_df.copy()
    for c in ["log2FoldChange", "padj", "baseMean"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.set_index("gene", drop=False)
    rows = []
    for g in markers:
        if g in df.index:
            r = df.loc[g]
            rows.append(
                {
                    "gene": g,
                    "log2FC": float(r.get("log2FoldChange")) if pd.notna(r.get("log2FoldChange")) else None,
                    "padj": float(r.get("padj")) if pd.notna(r.get("padj")) else None,
                    "baseMean": float(r.get("baseMean")) if pd.notna(r.get("baseMean")) else None,
                }
            )
        else:
            rows.append({"gene": g, "log2FC": None, "padj": None, "baseMean": None})
    out = pd.DataFrame(rows)
    return out


def tissue_page(
    meta: pd.DataFrame,
    tissue: str,
    comps: list[Comp],
) -> str:
    title = TISSUE_LABELS.get(tissue, TISSUE_CN.get(tissue, tissue))
    lines: list[str] = []

    lines.append(f"## {title} — 每器官 1 页解读")
    lines.append("")
    lines.append(f"- **组织**：`{tissue}`")
    lines.append(f"- **生成时间**：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("### 组别定义")
    lines.append("")
    for code, display in COND_MAP.items():
        lines.append(f"- **{display}**（内部代码: {code}）")
    lines.append("")

    # Sample n
    n_tbl = sample_n_table(meta, tissue)
    lines.append("### 样本量（按该组织实际存在的样本计数）")
    lines.append("")
    lines.append(_mk_table(n_tbl))
    lines.append("")

    lines.append("### 结果阅读说明（建议先看这段）")
    lines.append("")
    lines.append("- **差异表达（DE）**：`padj < 0.05` 认为显著；同时报告更严格阈值 `padj < 0.05 且 |log2FC|≥1` 的基因数。")
    lines.append("- **GSEA**：重点看 **NES 方向** 与 **FDR q-val**（常用阈值 0.25；更严格可用 0.05）。")
    lines.append('- **注意**：某些组别在某组织的重复数可能较少，统计功效较弱，相关结论应标注为"探索性"。')
    lines.append("")

    for comp in comps:
        lines.append(f"### 比较：{comp.label}")
        lines.append("")

        comp_key = comp.key

        # DE
        de = load_de(tissue, comp_key)
        if de is None:
            lines.append("- **DE**：未找到对应结果文件。")
            lines.append("")
        else:
            sde = summarize_de(de)
            lines.append(
                f"- **DE 总结**：基因数 {sde['n_genes']}；显著（padj<0.05）{sde['n_sig_padj']}；显著且 |log2FC|≥1：{sde['n_sig_padj_lfc1']}。"
            )
            lines.append("")
            lines.append("#### 关键上调基因（Top，padj<0.05）")
            lines.append("")
            lines.append(_mk_table(sde["top_up"]))
            lines.append("")
            lines.append("#### 关键下调基因（Top，padj<0.05）")
            lines.append("")
            lines.append(_mk_table(sde["top_down"]))
            lines.append("")

            # Markers
            markers = MARKERS.get(tissue, [])
            if markers:
                mt = marker_table(de, markers)
                lines.append("#### 预设 marker 基因（便于快速定位代谢/炎症/重构信号）")
                lines.append("")
                lines.append(_mk_table(mt))
                lines.append("")

        # GSEA
        gsea = load_gsea(tissue, comp_key)
        if gsea is None:
            lines.append("- **GSEA**：未找到对应结果文件。")
            lines.append("")
        else:
            sg = summarize_gsea(gsea)
            lines.append(f"- **GSEA 总结**：条目数 {sg['n_terms']}；FDR≤0.25 的条目数 {sg['n_fdr25']}。")
            lines.append("")
            lines.append("#### GSEA：NES 正向（Top）")
            lines.append("")
            lines.append(_mk_table(sg["top_pos"]))
            lines.append("")
            lines.append("#### GSEA：NES 负向（Top）")
            lines.append("")
            lines.append(_mk_table(sg["top_neg"]))
            lines.append("")

            lines.append("#### 关键词聚焦（便于快速抓代谢/免疫/纤维化方向；来自全部 term 的匹配）")
            lines.append("")
            lines.append(f"- **代谢相关命中**（最多 20 条）：{'; '.join(sg['metabolic_hits']) if sg['metabolic_hits'] else '（无）'}")
            lines.append(f"- **免疫/炎症命中**（最多 20 条）：{'; '.join(sg['immune_hits']) if sg['immune_hits'] else '（无）'}")
            lines.append(f"- **纤维化/ECM 命中**（最多 20 条）：{'; '.join(sg['fibrosis_hits']) if sg['fibrosis_hits'] else '（无）'}")
            lines.append("")

        # Enrichr
        up = summarize_enrichr(load_enrichr(tissue, comp_key, "up"))
        down = summarize_enrichr(load_enrichr(tissue, comp_key, "down"))
        if up.empty and down.empty:
            lines.append("- **GO/KEGG/Reactome/WikiPathways 富集（Enrichr）**：未找到对应结果文件。")
            lines.append("")
        else:
            lines.append("#### 富集（Enrichr）：上调基因集合（Top）")
            lines.append("")
            lines.append(_mk_table(up))
            lines.append("")
            lines.append("#### 富集（Enrichr）：下调基因集合（Top）")
            lines.append("")
            lines.append(_mk_table(down))
            lines.append("")

        lines.append("#### 本比较的证据优先解读模板（你可以直接复制到汇报PPT）")
        lines.append("")
        lines.append("'- **一句话结论（建议写法）**：在该组织中，`{A} vs {B}` 引起的主要变化可概括为'(代谢/免疫/结构重塑)相关通路的(增强/抑制)'，证据来自（Top GSEA term + Top DEG/marker）。'".format(A=comp.a, B=comp.b))
        lines.append("- **最强证据 3 条（建议写法）**：")
        lines.append("  - GSEA：列出 1–2 个 NES 绝对值最大且 FDR 较小的 term（写方向）。")
        lines.append("  - DE：列出 1–2 个与你生物学假设最贴近的基因（含 log2FC 与 padj）。")
        lines.append("  - 若有：Enrichr 的 top term 作为补充。")
        lines.append("")

    lines.append("### 下一步建议（实验与再分析）")
    lines.append("")
    lines.append('- **统计功效**：优先检查每组重复数；若某组织某组 n 很少，该组织的结论应标注为"探索性"。')
    lines.append("- **验证路径**：qPCR（marker）、WB/ELISA（关键通路蛋白/细胞因子）、组织学（炎症/纤维化/脂肪变）、功能测定（线粒体呼吸、脂代谢通量）。")
    lines.append('- **建议补充图**：按组织做"**三比较并排**"，一页就能讲清方向。')
    lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate per-organ interpretation report driven by experiment config.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None, help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None, help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()

    cfg_info = _apply_config(args.experiment, args.config, args.outdir)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    meta = load_metadata()
    tissues = sorted(meta["tissue"].unique().tolist())

    comps = [Comp(k, a, b) for k, (a, b) in COMP_MAP.items()]

    # Per tissue pages
    created: list[str] = []
    for tissue in tissues:
        page = tissue_page(meta, tissue, comps)
        out_path = OUT_DIR / f"{tissue}.md"
        out_path.write_text(page, encoding="utf-8")
        created.append(out_path.name)

    # Index
    index_lines: list[str] = []
    index_lines.append("## Bulk RNA-Seq 解读报告（每器官 1 页）")
    index_lines.append("")
    index_lines.append(f"- **生成时间**：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if len(COND_MAP) >= 2:
        names = " / ".join(COND_MAP.values())
        index_lines.append(f"- **组别命名（按 groups.yaml）**：{names}")
    else:
        index_lines.append("- **组别命名**：（见 config）")
    index_lines.append("")
    index_lines.append("### 器官目录")
    index_lines.append("")
    for tissue in tissues:
        cn = TISSUE_LABELS.get(tissue, TISSUE_CN.get(tissue, tissue))
        index_lines.append(f"- [{cn}](./{tissue}.md)")
    index_lines.append("")
    index_lines.append("### 重要备注")
    index_lines.append("")
    index_lines.append('- 本报告为"**证据汇总 + 解读模板**"，优先展示可复核的表格（DE/GSEA/富集）。')
    index_lines.append('- 重复数偏低的组别/组织，相关结论请标注为"探索性"。')
    index_lines.append("")

    (OUT_DIR / "index.md").write_text("\n".join(index_lines).rstrip() + "\n", encoding="utf-8")

    print(f"[OK] Wrote {len(created)} tissue pages + index to: {OUT_DIR}")
    print("Files:")
    for fn in ["index.md"] + created:
        print(f"  - {fn}")


if __name__ == "__main__":
    main()

