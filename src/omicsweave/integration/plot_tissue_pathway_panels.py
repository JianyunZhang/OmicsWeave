"""绘制经证据表核查的组织通路、动物层箱线图及MFA拟合范围图。"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.cbook import boxplot_stats
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.text import Text
import numpy as np
import pandas as pd


INK = "#193547"
GRAY = "#667985"
UP = "#0072B2"
DOWN = "#D55E00"
LIGHT = "#F4F8FA"
BORDER = "#CEDCE3"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Microsoft YaHei", "DejaVu Sans", "Arial"],
    "axes.unicode_minus": False, "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
})


def save(fig, folder, stem):
    for extension in ["svg", "pdf", "png"]:
        fig.savefig(folder / f"{stem}.{extension}", facecolor="white", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fmt(value):
    return "NA" if pd.isna(value) else f"{value:.3g}"


def grayscale_plot_artists(fig):
    """在矢量绘制阶段生成真正的灰度QA版本，保留透明度和文字。"""
    def gray(color):
        r, g, b, a = to_rgba(color)
        light = .2126 * r + .7152 * g + .0722 * b
        return (light, light, light, a)
    for artist in fig.findobj():
        if isinstance(artist, (Text, Line2D)):
            artist.set_color(gray(artist.get_color()))
        if isinstance(artist, Patch):
            artist.set_facecolor(gray(artist.get_facecolor()))
            artist.set_edgecolor(gray(artist.get_edgecolor()))


class Pathway:
    """灰色反应箭头仅表示已知化学；彩色方向只来自对应层的FDR证据。"""

    def __init__(self, ax, evidence, monochrome=False, settings=None):
        self.ax = ax
        self.evidence = evidence.set_index(["layer", "node_id"])
        self.monochrome = monochrome
        self.settings = settings or {}
        self.annotation_log = []
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 140)
        ax.axis("off")

    def text(self, x, y, label, size=11, color=INK, ha="left", bold=False, **kwargs):
        return self.ax.text(x, y, label, fontsize=size, color=color, ha=ha,
                            fontweight="bold" if bold else "normal", **kwargs)

    def node(self, x, y, label, width=15, height=12):
        self.ax.add_patch(FancyBboxPatch((x - width / 2, y - height / 2), width, height,
                          boxstyle="round,pad=0.4,rounding_size=1.1", facecolor=LIGHT, edgecolor=BORDER, lw=1.4))
        self.text(x, y, label, 11.5, ha="center", va="center")

    def edge(self, x1, y1, x2, y2, label=None, curve=0, label_y=None):
        self.ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops={
            "arrowstyle": "-|>", "lw": 1.7, "color": "#727272" if self.monochrome else GRAY,
            "connectionstyle": f"arc3,rad={curve}"})
        if label:
            self.text((x1 + x2) / 2, label_y if label_y is not None else (y1 + y2) / 2 + 10,
                      label, 10, GRAY, ha="center")

    def obs(self, x, y, layer, node_id, label=None, size=12, ha="center", q_offset=5.2):
        row = self.evidence.loc[(layer, node_id)]
        label = label or node_id
        if pd.isna(row.padj):
            direction, color = "NA", GRAY
        elif row.padj >= .05:
            direction, color = "NS", GRAY
        else:
            direction = "↑" if row.effect > 0 else "↓"
            color = UP if row.effect > 0 else DOWN
        if self.monochrome:
            color = INK if direction in ["↑", "↓"] else GRAY
        self.text(x, y, f"{label} {direction}", size, color, ha=ha, bold=direction in ["↑", "↓"],
                  fontstyle="italic" if layer == "RNA" else "normal")
        q_label = "BH q" if layer == "lipid_class" else "padj"
        self.text(x, y - q_offset, f"{q_label} = {fmt(row.padj)}", 9, GRAY, ha=ha)
        self.annotation_log.append({"layer": layer, "node_id": node_id, "label": label,
                                    "display_direction": direction, "padj": row.padj, "effect": row.effect})

    def divider(self, y):
        self.ax.plot([2, 98], [y, y], color=BORDER, lw=1.2)


def lipid_chain(p, y, *, genes=False, start=8, end=88):
    xs = np.linspace(start, end, 5)
    for x, label in zip(xs, ["G3P", "LPA", "PA", "DAG", "TAG"]):
        p.node(x, y, label, width=13)
    for index, enzyme in enumerate(["GPAT", "AGPAT", "LPIN / PAP", "DGAT1/2"]):
        p.edge(xs[index] + 7, y, xs[index + 1] - 7, y)
        mid = (xs[index] + xs[index + 1]) / 2
        if genes and index in [2, 3]:
            p.obs(mid, y + 14, "RNA", "Lpin1" if index == 2 else "Dgat2")
        else:
            p.text(mid, y + 10, enzyme, 9.5, GRAY, ha="center")
    return xs


def draw_heart(p):
    p.text(2, 119, "Long-chain fatty-acid handling: measured carnitine signals", 12.5, bold=True)
    for x, label in [(10, "LC acyl-CoA\n+ carnitine"), (36, "Acylcarnitine"),
                     (63, "Acylcarnitine\n(matrix)"), (90, "LC acyl-CoA\n+ carnitine")]:
        p.node(x, 98, label, width=18)
    p.edge(19.5, 98, 26.5, 98, "CPT1", label_y=112)
    p.edge(45.5, 98, 53.5, 98, "SLC25A20", label_y=112)
    p.text(49.5, 107, "acyl-CAR / CAR exchange", 8.5, GRAY, ha="center")
    p.edge(72.5, 98, 80.5, 98, "CPT2", label_y=112)
    p.text(77, 90, "+ CoA", 9, GRAY, ha="center")
    p.obs(35, 81, "metabolite", "LIPID-P-0021", "C16:0 acylcarnitine feature", size=12)
    p.obs(3, 65, "lipid_class", "Acylcarnitine", "Carnitine-related class", size=12, ha="left")
    p.text(3, 54, "Includes free carnitine; class score is a relative signal.", 9.5, GRAY)
    p.text(76, 66, "Downstream beta-oxidation\nflux remains unresolved.", 10, GRAY, ha="center")
    p.divider(47)
    p.text(2, 41, "G3P / TAG context and the selected RNA panel", 12, bold=True)
    lipid_chain(p, 20)
    p.obs(7, 4, "lipid_class", "TG", "TAG class", size=10, ha="left", q_offset=4)
    genes = p.settings["selected_gene_panel"]
    sub = p.evidence.loc[[("RNA", gene) for gene in genes]]
    p.text(42, 3, f"Selected RNA panel: {sub.padj.lt(.05).sum()}/{len(genes)} FDR-significant; min padj = {sub.padj.min():.3g}", 10, GRAY)


def draw_liver(p):
    p.text(2, 119, "G3P backbone and TAG synthesis-related RNA", 12.5, bold=True)
    xs = lipid_chain(p, 98, genes=True)
    p.obs(86, 82, "lipid_class", "TG", "TAG class", size=10, q_offset=4)
    p.divider(73)
    p.text(2, 68, "Kennedy-PC branch", 12.5, bold=True)
    for x, label in [(10, "Choline"), (35, "Phospho-\ncholine"), (61, "CDP-choline"), (87, "PC")]:
        p.node(x, 47, label, width=17)
    p.edge(19, 47, 26, 47)
    p.obs(22.5, 61, "RNA", "Chka")
    p.edge(44, 47, 52, 47, "PCYT1A/B", label_y=58)
    p.edge(70, 47, 78, 47)
    p.obs(74, 61, "RNA", "Chpt1", size=11)
    p.text(87, 36, "+ DAG → PC + CMP (CHPT1 / CEPT1)", 8.5, GRAY, ha="center")
    p.divider(29)
    p.text(2, 24, "Adjacent glutamine node", 12.5, bold=True)
    p.node(20, 10, "Glutamate + NH4+", width=28, height=10)
    p.node(77, 10, "Glutamine", width=24, height=10)
    p.edge(34.5, 10, 64.5, 10)
    p.obs(49, 23, "RNA", "Glul", size=11, q_offset=4)
    p.text(49, 0, "ATP → ADP + Pi", 9, GRAY, ha="center")


def draw_adipose(p):
    p.text(2, 119, "Glycerolipid backbone: selected species-level changes", 12.5, bold=True)
    lipid_chain(p, 99)
    p.obs(28, 82, "metabolite", "MW0011785HN", "LPA(20:1) feature", size=11.5)
    p.obs(69, 82, "metabolite", "LIPID-P-0211", "DG(16:1_16:1) feature", size=11.5)
    p.text(3, 67, "These two features do not represent total LPA or total DAG pools.", 9.5, GRAY)
    p.divider(61)
    p.text(2, 55, "GPD1–GPD2 redox-shuttle context", 12.5, bold=True)
    p.node(13, 32, "DHAP", width=17)
    p.node(47, 32, "G3P", width=17)
    p.edge(22, 36, 38, 36, curve=-.3)
    p.text(30, 47, "GPD1 · NADH / NAD+", 10, GRAY, ha="center")
    p.edge(38, 28, 22, 28, curve=-.3)
    p.obs(30, 19, "RNA", "Gpd2", size=12)
    p.text(30, 7, "Electron transfer: bound FAD/FADH2\n→ Q/QH2 → ETC", 9.5, GRAY, ha="center")
    p.text(72, 55, "Thermogenic / fuel-handling RNA", 11.5, bold=True, ha="center")
    for y, gene in [(45, "Ucp1"), (32, "Cpt1b"), (19, "Cox7a1"), (6, "Cidea")]:
        p.obs(78, y, "RNA", gene, size=12, q_offset=4.5)


def draw_muscle(p):
    p.text(2, 119, "Glucose entry and a glycolytic intermediate", 12.5, bold=True)
    for x, label in [(9, "Glucose"), (31, "G6P"), (55, "Glycolytic\nintermediates"), (81, "2-PG")]:
        p.node(x, 98, label, width=16)
    p.edge(17.5, 98, 22.5, 98)
    p.obs(20, 112, "RNA", "Hk2", size=11)
    p.edge(39.5, 98, 46.5, 98, "multiple steps", label_y=108)
    p.edge(63.5, 98, 72.5, 98, "multiple steps", label_y=108)
    p.obs(81, 81, "metabolite", "MEDN2351", "2-PG relative signal", size=11)
    p.divider(71)
    p.text(2, 65, "G3P backbone and TAG synthesis-related RNA", 12.5, bold=True)
    for x, label in [(8, "DHAP"), (30, "G3P"), (57, "DAG"), (86, "TAG")]:
        p.node(x, 44, label, width=14)
    p.edge(15.5, 44, 22.5, 44, "GPD1", label_y=54)
    p.edge(37.5, 44, 49.5, 44, "via LPA / PA", label_y=54)
    p.edge(64.5, 44, 78.5, 44)
    p.obs(71.5, 58, "RNA", "Dgat2", size=11)
    p.divider(31)
    p.text(2, 25, "Adjacent glutamine node: RNA and metabolite signal", 12, bold=True)
    p.node(18, 10, "Glutamate + NH4+", width=27, height=10)
    p.node(75, 10, "Glutamine", width=25, height=10)
    p.edge(32, 10, 62, 10)
    p.obs(47, 23, "RNA", "Glul", size=11, q_offset=4)
    p.obs(77, 0, "metabolite", "MADN0217", "L-Gln signal", size=10, q_offset=4)
    p.text(45, 0, "ATP → ADP + Pi", 9, GRAY, ha="center")


def draw_renal(p):
    p.text(2, 119, "Kennedy-PC branch and phospholipid remodeling", 12.5, bold=True)
    for x, label in [(10, "Choline"), (35, "Phospho-\ncholine"), (61, "CDP-choline"), (87, "PC")]:
        p.node(x, 97, label, width=17)
    p.edge(19, 97, 26, 97)
    p.obs(22.5, 111, "RNA", "Chkb")
    p.edge(44, 97, 52, 97)
    p.obs(48, 111, "RNA", "Pcyt1a")
    p.edge(70, 97, 78, 97, "CHPT1 / CEPT1", label_y=110)
    p.text(75, 104, "+ DAG", 9, GRAY, ha="center")
    p.node(60, 73, "LPC", width=17)
    p.edge(69, 73, 84, 90)
    p.text(73, 86, "+ acyl-CoA", 9, GRAY, ha="center")
    p.obs(84, 72, "RNA", "Lpcat4")
    p.text(2, 73, "Opposing RNA directions do not\nestablish net PC synthesis.", 10, GRAY)
    p.divider(59)
    p.text(2, 53, "Adjacent sphingolipid module", 12.5, bold=True)
    p.node(17, 32, "Dihydroceramide", width=27)
    p.node(50, 32, "Ceramide", width=19)
    p.edge(31, 32, 40, 32)
    p.obs(35, 46, "RNA", "Degs1", size=11)
    p.obs(76, 42, "lipid_class", "Ceramide", "Ceramide-family signal", size=11)
    p.text(76, 29, "154 features; includes CerP", 9, GRAY, ha="center")
    p.obs(69, 13, "metabolite", "LIPID-P-0332", "SHexCer(d18:1/16:0) feature", size=10.5)
    p.text(2, 12, "Sphingolipids are an adjacent module,\nnot a direct G3P-backbone fate.", 9.5, GRAY)


def draw_brain_mfa(p):
    p.text(2, 119, "Known G3P / glycerolipid connections: biochemical background", 12, bold=True)
    for x, label, width in [(9, "DHAP", 15), (35, "G3P", 15), (64, "DAG", 16), (89, "TAG", 16)]:
        p.node(x, 99, label, width=width)
    p.edge(17, 99, 27, 99, "GPD1", label_y=110)
    p.edge(43, 99, 55, 99, "via LPA / PA", label_y=110)
    p.edge(73, 99, 80, 99, "DGAT", label_y=110)
    p.node(34, 72, "Glycerol", width=19)
    p.edge(29, 92, 29, 79, curve=.18)
    p.text(18, 84, "PGP / G3PP", 9, GRAY, ha="center")
    p.edge(39, 79, 39, 92, curve=.18)
    p.text(47, 84, "GK", 9, GRAY, ha="center")
    summary = p.settings["model_summary"]
    p.text(66, 80, "Treatment data: MFA only", 11, bold=True)
    p.text(66, 72, f"{int(summary['mfa_mw_fdr_significant'])}/{int(summary['mfa_total_reactions'])} reactions pass MW FDR < 0.05", 10, GRAY)
    p.text(66, 64, "RNA / protein / relative metabolite\ndata absent for this treatment contrast", 9.5, GRAY)
    p.divider(55)
    p.text(2, 49, "Model readouts (aggregated nodes; no enzyme-specific flux inference)", 11.5, bold=True)
    for index, node_id in enumerate(p.settings["mfa_nodes"]):
        x = [17, 50, 83][index % 3]
        y = [38, 22, 6][index // 3]
        row = p.evidence.loc[("MFA", node_id)]
        p.text(x, y, node_id, 10.5, GRAY, ha="center", bold=True)
        p.text(x, y - 5.5, f"MW FDR = {fmt(row.padj)}", 9, GRAY, ha="center")
        p.annotation_log.append({"layer": "MFA", "node_id": node_id, "label": node_id,
                                 "display_direction": "NS" if row.padj >= .05 else "FDR_model",
                                 "padj": row.padj, "effect": row.effect})


def draw_brain_reference(p):
    p.text(2, 119, "Independent HFD/LFD reference: glycerolipid-related RNA", 12, bold=True)
    for x, label in [(10, "G3P"), (40, "DAG"), (75, "TAG")]:
        p.node(x, 98, label, width=19)
    p.edge(20, 98, 30, 98, "via LPA / PA", label_y=109)
    p.edge(50, 98, 65, 98)
    p.obs(57.5, 112, "RNA", "Dgat1")
    p.divider(82)
    p.text(2, 76, "Kennedy-PC entry", 12.5, bold=True)
    for x, label in [(10, "Choline"), (42, "Phospho-\ncholine"), (72, "PC")]:
        p.node(x, 56, label, width=18)
    p.edge(21, 56, 31, 56)
    p.obs(20, 69, "RNA", "Chka", size=10.5, ha="center")
    p.obs(33, 69, "RNA", "Chkb", size=10.5, ha="center")
    p.edge(53, 56, 61, 56, "via CDP-choline + DAG", label_y=68)
    p.divider(41)
    p.text(2, 35, "Adjacent glutamine node", 12, bold=True)
    p.node(20, 19, "Glutamate + NH4+", width=28, height=10)
    p.node(75, 19, "Glutamine", width=24, height=10)
    p.edge(35, 19, 62, 19)
    p.obs(48, 33, "RNA", "Glul", size=11)
    p.text(48, 8, "ATP → ADP + Pi", 9, GRAY, ha="center")
    protein = p.evidence.loc[[('protein', g) for g in ['Chka', 'Chkb', 'Dgat1', 'Glul']]]
    p.text(2, -2, f"Corresponding proteins: {protein.padj.lt(.05).sum()}/4 FDR-significant. This is not a Mazdutide contrast.", 10, GRAY)


DRAWERS = {"heart": draw_heart, "liver": draw_liver, "adipose": draw_adipose, "muscle": draw_muscle,
           "renal": draw_renal, "brain_mfa": draw_brain_mfa}


def pathway_figure(directory, tissue, settings, config, monochrome=False):
    evidence = pd.read_csv(directory / tissue / f"{tissue}_pathway_evidence.csv")
    fig, ax = plt.subplots(figsize=(14, 9.4))
    settings = dict(settings)
    if not settings.get("molecular_available", True):
        settings["model_summary"] = pd.read_csv(directory / "evidence_summary.csv").set_index("tissue").loc[tissue].to_dict()
    p = Pathway(ax, evidence, monochrome, settings)
    p.text(2, 137, f"Panel B | {settings['display']}", 21, bold=True, va="top")
    molecular_available = settings.get("molecular_available", True)
    counts = config["expected_counts"]["molecular" if molecular_available else "mfa"]
    layer_label = "RNA / metabolite" if molecular_available else "MFA"
    p.text(2, 127, f"{config['case_label']} vs {config['control_label']}  ·  {config['scope_label']}  ·  {layer_label} animals: n = {counts['case']} vs {counts['control']}", 11, GRAY)
    DRAWERS[settings["layout"]](p)
    footer = ("Grey reaction arrows: known biochemical connections (co-substrates simplified). Coloured ↑/↓: FDR-supported RNA or relative-signal changes.\n"
              "NS: not FDR-significant; this does not establish no response. Pathway flux, enzyme activity and organ function are not inferred from these arrows.") if molecular_available else (
              "Grey connections show biochemical context. Model labels denote aggregated, exploratory MFA readouts; no treatment-direction arrows are assigned.\n"
              "Independent HFD/LFD brain omics are kept in a separate reference folder and do not supply molecular evidence for this treatment contrast.")
    fig.text(.14, .018, footer, fontsize=9, color=GRAY)
    stem = f"panelB_{tissue}_pathway" + ("_grayscale_QA" if monochrome else "")
    if monochrome:
        grayscale_plot_artists(fig)
    save(fig, directory / tissue, stem)
    if not monochrome:
        pd.DataFrame(p.annotation_log).to_csv(directory / tissue / f"{tissue}_pathway_annotations.csv", index=False, encoding="utf-8-sig")


def boxes(ax, frame, groups, summary_rows, tissue, node_id, layer):
    for index, group in enumerate(groups):
        sub = frame[frame.group.eq(group)].sort_values("sample_id")
        values = sub.value.to_numpy(float)
        assert np.isfinite(values).all()
        box = boxplot_stats(values, whis=1.5)[0]
        color = GRAY if index == 0 else UP
        ax.bxp([box], positions=[index], widths=.44, showfliers=False, patch_artist=True,
               boxprops={"facecolor": "white", "edgecolor": color, "linewidth": 1.3},
               medianprops={"color": INK, "linewidth": 1.5}, whiskerprops={"color": color}, capprops={"color": color})
        xx = index + np.linspace(-.12, .12, len(values))
        ax.scatter(xx, values, color=color, s=31, marker="o" if index == 0 else "^", zorder=4)
        row = {"tissue": tissue, "layer": layer, "node_id": node_id, "group": group, "n_animals": len(values),
               **{key: float(box[key]) for key in ["mean", "med", "q1", "q3", "whislo", "whishi"]},
               "n_tukey_outliers": len(box["fliers"]), "whisker_rule": "1.5_IQR", "all_animals_displayed": True}
        for column in ["comparison", "data_role", "is_treatment_comparison", "value_scale"]:
            if column in frame:
                row[column] = frame[column].iloc[0]
        summary_rows.append(row)
    ax.set_xticks([0, 1], [f"{g}\nn = {len(frame[frame.group.eq(g)])}" for g in groups])
    ax.set_xlim(-.5, 1.5)
    ax.grid(axis="y", alpha=.13)
    ax.tick_params(labelsize=8.5)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)


def boxplot_figure(directory, tissue, settings, config, mode, folder=None):
    folder = folder or directory / tissue
    data = pd.read_csv(folder / f"{tissue}_{mode}_plot_data.csv")
    nodes = data[["layer", "node_id", "node_order"]].drop_duplicates().sort_values("node_order")
    ncols = 3
    nrows = int(np.ceil(len(nodes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.8, 3.2 * nrows + 1.15), squeeze=False)
    fig.subplots_adjust(hspace=.68, wspace=.36, top=.88, bottom=.13)
    summaries = []
    groups = [config["control_label"], config["case_label"]]
    quality_labels = {"poorly_constrained": "约束不足", "moderate": "中等约束", "well_constrained": "较强约束", "undetermined": "未确定"}
    for index, (ax, (_, node)) in enumerate(zip(axes.ravel(), nodes.iterrows())):
        node_id = node.node_id
        sub = data[data.node_id.eq(node_id) & data.layer.eq(node.layer)]
        row = sub.iloc[0]
        boxes(ax, sub, groups, summaries, tissue, node_id, row.layer)
        label = textwrap.fill(str(row.node_label), width=29)
        ax.set_title(f"{chr(65 + index)}  {label}", loc="left", fontsize=11, fontweight="bold", pad=28)
        if mode == "MFA":
            q = row.fdr_mannwhitney
            qualifier = quality_labels.get(row.identifiability, str(row.identifiability))
            near_one = int(sub.near_one_flag.sum())
            note = f"  |  近1：{near_one}/{len(sub)}只" if near_one else ""
            axis_rule = "automatic"
            if near_one == len(sub):
                # 避免将1附近的浮点数值误差画成明显的生物学差异。
                ax.set_ylim(.9, 1.1)
                qualifier, note = "全部近1，需核查边界", ""
                axis_rule = "all_near_one_fixed_0.9_1.1"
            elif sub.near_zero_flag.all():
                ax.set_ylim(-.05, .05)
                qualifier, note = "全部近0", ""
                axis_rule = "all_near_zero_fixed_minus0.05_0.05"
            ax.text(.01, 1.02, f"MW FDR = {fmt(q)}  |  {qualifier}{note}", transform=ax.transAxes, fontsize=8.5, color=GRAY)
            ax.set_ylabel("动物层拟合净通量（原模型尺度）", fontsize=9)
            if row.identifiability == "poorly_constrained":
                ax.set_facecolor("#FFFCF8")
        else:
            axis_rule = "automatic"
            q = row.source_padj
            q_name = "七类别 BH q" if row.layer == "lipid_class" else "padj"
            ax.text(.01, 1.02, f"{row.layer}  |  {q_name} = {fmt(q)}", transform=ax.transAxes, fontsize=9,
                    color=UP if pd.notna(q) and q < .05 else GRAY)
            scales = {"log2_CPM_plus_1": "RNA log2(CPM + 1)", "mean_log2_PQN_signal": "类别内平均 log2 归一化信号", "log2_PQN_signal": "log2 归一化相对信号", "log2_normalized_unimputed_protein_signal": "蛋白log2归一化信号（未插补）"}
            ax.set_ylabel(scales[row.value_scale], fontsize=9)
        for record in summaries:
            if record["node_id"] == node_id and record["layer"] == node.layer:
                record.update({"axis_rule": axis_rule, "axis_ymin": ax.get_ylim()[0], "axis_ymax": ax.get_ylim()[1]})
    for ax in axes.ravel()[len(nodes):]:
        ax.axis("off")
    if mode == "MFA":
        title = f"{settings['chinese']} | 关键MFA节点的动物层箱线图"
        footer = ("箱体：动物层估计值的Q1–Q3；中线：中位数；须：1.5 IQR；每点一只动物（拟合迭代不作为动物重复）。\n"
                  "每动物值为现有管线过滤后拟合净通量的中位数。FDR沿用同一比较全部组织×反应的校正；源效应量入选标记不等于FDR显著。\n"
                  "图用于动物层的探索性展示；源可辨识性标签与接近1标志用于提示模型约束问题，拟合P25/P75见配套CSV和QA图。")
    else:
        title = f"{settings['chinese']} | 关键分子节点的动物层箱线图"
        footer = ("箱体：动物值Q1–Q3；中线：中位数；须：1.5 IQR；所有动物均叠加显示。RNA图使用log2(CPM+1)，p值来自原DESeq2模型。\n"
                  "RNA padj：全转录组；单代谢特征padj：该组织全部受检特征；脂质类别q：每组织每比较七类别BH。不同数据层的数值尺度不可直接相加。\n"
                  "代谢物为现有注释支持的相对信号，单个脂质种不代表整个类别；多个特征或多个组织不能当作增加的动物重复。")
        if config.get("reference_mode"):
            footer = ("独立HFD/LFD造模参照，不是Mazdutide治疗比较。箱体：动物值Q1–Q3；中线：中位数；须：1.5 IQR；全部动物点叠加显示。\n"
                      "RNA使用log2(CPM+1)，蛋白使用未插补的log2归一化信号；padj分别沿用原全转录组、可检验蛋白和全代谢特征的校正范围。\n"
                      "原始组织标签为brain，细分脑区未明确；不自动等同治疗MFA的cerebral_cortex。")
    fig.suptitle(title, fontsize=16, fontweight="bold", color=INK, y=.978)
    fig.text(.5, .94, f"{config['case_label']} vs {config['control_label']} · {config['scope_label']}", ha="center", fontsize=11, color=GRAY)
    fig.text(.09, .035, footer, fontsize=8.7, color=GRAY)
    stem = f"{tissue}_{mode}_key_nodes_boxplots"
    save(fig, folder, stem)
    pd.DataFrame(summaries).to_csv(folder / f"{tissue}_{mode}_box_statistics.csv", index=False, encoding="utf-8-sig")


def mfa_fit_qa(directory, tissue, settings, config):
    data = pd.read_csv(directory / tissue / f"{tissue}_MFA_plot_data.csv")
    nodes = data.sort_values("node_order").node_id.unique()
    fig, axes = plt.subplots(3, 3, figsize=(13.8, 10.6))
    fig.subplots_adjust(hspace=.63, wspace=.34, top=.88, bottom=.13)
    groups = [config["control_label"], config["case_label"]]
    for index, (ax, node) in enumerate(zip(axes.ravel(), nodes)):
        subset = data[data.node_id.eq(node)]
        for i, group in enumerate(groups):
            selected = subset[subset.group.eq(group)].sort_values("sample_id")
            for x, (_, row) in zip(i + np.linspace(-.13, .13, len(selected)), selected.iterrows()):
                color = GRAY if i == 0 else UP
                ax.plot([x, x], [row.fit_p25, row.fit_p75], color=color, lw=1.5, alpha=.8)
                ax.scatter(x, row.value, color=color, s=30, marker="o" if i == 0 else "^", zorder=3)
        ax.set_title(f"{chr(65 + index)}  {node}\n{subset.identifiability.iloc[0]}", loc="left", fontsize=10, fontweight="bold")
        ax.set_xticks([0, 1], [f"{group}\nn = {len(subset[subset.group.eq(group)])}" for group in groups])
        ax.set_xlim(-.5, 1.5)
        ax.set_ylabel("拟合净通量（原模型尺度）", fontsize=9)
        ax.grid(axis="y", alpha=.15)
        ax.tick_params(labelsize=8.5)
        if np.isclose(subset[["fit_p25", "fit_p75"]].to_numpy(float), 1, atol=1e-3, rtol=0).all():
            ax.set_ylim(.9, 1.1)
        elif np.isclose(subset[["fit_p25", "fit_p75"]].to_numpy(float), 0, atol=1e-6, rtol=0).all():
            ax.set_ylim(-.05, .05)
    fig.suptitle(f"{settings['chinese']} | 同一关键节点的MFA拟合范围核查", fontsize=16, fontweight="bold", color=INK, y=.977)
    fig.text(.5, .94, f"{config['case_label']} vs {config['control_label']} · {config['scope_label']}", ha="center", fontsize=11, color=GRAY)
    fig.text(.09, .043, "每点一只动物；竖线为该动物过滤后拟合点的P25–P75，不是组间置信区间，也不增加动物数。\n"
             "此图与MFA箱线图使用同一份 *_MFA_plot_data.csv；箱线图展示动物间分布，本图展示每只动物内部的模型拟合范围。", fontsize=9, color=GRAY)
    save(fig, directory / tissue, f"{tissue}_MFA_fit_QA")


def plot_brain_reference(directory, config):
    reference = config["brain_reference"]
    folder = directory / reference["output_subdir"]
    evidence = pd.read_csv(folder / "brain_reference_pathway_evidence.csv")
    settings = {"display": reference["display"], "chinese": reference["chinese"]}
    ref_config = dict(config)
    ref_config.update({key: reference[key] for key in ["case_label", "control_label", "scope_label"]})
    ref_config["reference_mode"] = True
    for mono in [False, True]:
        fig, ax = plt.subplots(figsize=(14, 9.4))
        p = Pathway(ax, evidence, mono)
        p.text(2, 137, "Brain | HFD vs LFD reference", 21, bold=True, va="top")
        p.text(2, 127, f"Independent cohort {reference['experiment_id']} · 14 weeks · n = 3 vs 3 · not a treatment contrast", 10.5, GRAY)
        draw_brain_reference(p)
        fig.text(.14, .018, "Grey arrows: known biochemical connections. Coloured arrows: RNA changes in the independent HFD/LFD cohort only.\n"
                 "Original tissue label: brain (subregion unspecified). These observations cannot substitute for Mazdutide treatment molecular data.", fontsize=9, color=GRAY)
        if mono:
            grayscale_plot_artists(fig)
        save(fig, folder, "brain_reference_pathway" + ("_grayscale_QA" if mono else ""))
        if not mono:
            pd.DataFrame(p.annotation_log).to_csv(folder / "brain_reference_pathway_annotations.csv", index=False, encoding="utf-8-sig")
    boxplot_figure(directory, "brain_reference", settings, ref_config, "molecular", folder=folder)
    (folder / "README_HFD_LFD参照.md").write_text(
        "# 脑HFD/LFD独立参照\n\n本目录全部来自20250425队列HFD-14W vs LFD-14W，非Mazdutide治疗比较。\n\n"
        "- [通路图](brain_reference_pathway.svg)\n- [9节点箱线图](brain_reference_molecular_key_nodes_boxplots.svg)\n"
        "- [逐动物绘图CSV](brain_reference_molecular_plot_data.csv)\n- [箱体统计CSV](brain_reference_molecular_box_statistics.csv)\n"
        "- [源差异证据](brain_reference_pathway_evidence.csv)\n\n"
        "9个节点包括Chka/Chkb/Dgat1/Glul的RNA、对应4个蛋白及HexCer(d18:1/25:0)注释特征。"
        "RNA和蛋白的node_id可相同，应联合layer区分，不能合并为同一个测量。所有节点设计n=3 vs 3，蛋白图使用未插补归一化矩阵。\n\n"
        "这些参考记录未加入正式治疗合并表；脑的正式治疗结果位于上一级目录，仅包含MFA。\n",
        encoding="utf-8")


def write_report(directory, config):
    """从最终绘图数据生成可追溯的说明，避免文案与CSV分别维护而漂移。"""
    summary = pd.read_csv(directory / "evidence_summary.csv").set_index("tissue")
    evidence = pd.read_csv(directory / "all_tissues_pathway_evidence.csv")
    manifest = json.loads((directory / "source_manifest.json").read_text(encoding="utf-8"))
    selected_stats = []
    descriptions = {
        "heart": "心脏主线是肉碱相关相对信号响应：C16:0酰基肉碱注释特征和肉碱相关类别上升。所选G3P/穿梭九基因无FDR显著差异。不能将肉碱池升高解释为CPT活性、脂肪酸氧化效率或心功能改善。",
        "liver": "肝脏可并列展示Lpin1/Dgat2转录下降、Kennedy入口Chka上调和相邻Glul下调。Chpt1、TAG类别和PC/LPC类别未达FDR，因此不能把入口酶上调归纳为整个PC通量增加，或把两项TAG相关酶下降归纳为TAG清除。另有Cidea转录下降（见完整RNA候选表），可作为脂滴相关的扩展候选。",
        "sq_adipose": "皮下脂肪可展示Ucp1/Cpt1b/Cox7a1/Cidea转录上调，以及单个LPA(20:1)和DG(16:1_16:1)注释特征下降；Gpd2未达FDR。这区分了产热相关RNA响应、G3P穿梭相关RNA和具体脂质种变化，不代表已测定产热、全部LPA/DAG池或穿梭通量。Slc25a13另有FDR显著转录上调，归入扩展的MAS转运候选表，不据此推断整个MAS活化。",
        "gast_muscle": "腓肠肌新增的组合证据包括Glul转录与L-glutamine相对信号同向下降、Dgat2下调，以及Hk2和2-PG相对信号下降。2-PG是2-磷酸甘油酸，不能与G3P（甘油-3-磷酸）混为同一节点。谷氨酰胺反应作为相邻氮代谢模块单列，不从G3P直接画碳流箭头通向GLUL。Pdk4/Ucp3等扩展节点也在箱线图及候选表中。当前GLNS_c模型估计反而表观上升，但未达FDR且约束不足，不能用该MFA结果为RNA/代谢物方向作通量验证。",
        "renal_cortex": "肾皮质补齐Kennedy-PC步骤间分离响应（Chkb上调、Pcyt1a下调、Lpcat4上调），并展示相邻鞘脂模块的Degs1下调、神经酰胺相关类别和SHexCer注释特征下降。Degs2/Asah2及G3P也纳入箱线图。排除CerP的类别定义敏感性另有CSV，不能将154个特征当作154只动物，也不能从这些分子读出认定肾功能获益。",
        "brain": "脑皮层的正式队列B结果仅有同期Mazdutide治疗MFA。本次补齐灰色通路背景及9个聚合模型节点的箱线图、拟合范围QA和逐动物CSV。治疗RNA/蛋白/常规代谢物箱线图标为数据不适用，不用另一队列HFD/LFD的数据填充。独立脑参考组学已在HFD_LFD_reference子目录补齐通路图、9个分子节点箱线图及逐动物CSV。",
    }
    lines = ["# 队列B六组织正式结果与独立脑参照", "",
             f"分析日期：{config['analysis_date']}。比较：{manifest['comparison']}；范围：{config['scope_label']}。", "",
             "主目录名称沿用panel_b_four_tissues_20260910，正式结果已扩展至heart、liver、sq_adipose、gast_muscle、renal_cortex、brain六个组织。脑的独立HFD/LFD三组学未混入同期治疗证据。", "",
             "六组织均提供通路SVG/PDF/PNG、9个MFA节点箱线图、同节点MFA拟合范围QA图，以及逐动物CSV和箱体统计CSV。五个具备同期分子数据的组织另有9个分子节点箱线图；脑的分子箱线图仅在独立HFD/LFD参照子目录提供。彩色方向标注只用于该数据层FDR支持的变化，灰色反应箭头表示已知生化连接。", "",
             "[完整六组织交付包：图件、CSV和复现脚本](panelB_six_tissues_figures_and_csv.zip)。旧的四组织ZIP保留为此前交付快照；本页及六组织ZIP为正式入口。", "",
             "## 图件与绘图数据", "",
             "| 组织 | 通路图 | 分子节点箱线图 | MFA节点箱线图 | MFA拟合范围QA |",
             "|---|---|---|---|---|"]
    for tissue, settings in config["tissues"].items():
        molecule_link = f"[SVG]({tissue}/{tissue}_molecular_key_nodes_boxplots.svg)" if settings.get("molecular_available", True) else "不适用：缺少同期治疗分子数据"
        lines.append(f"| {settings['chinese']} | [SVG]({tissue}/panelB_{tissue}_pathway.svg) | {molecule_link} | [SVG]({tissue}/{tissue}_MFA_key_nodes_boxplots.svg) | [SVG]({tissue}/{tissue}_MFA_fit_QA.svg) |")
    lines += ["", "| 组织 | 分子逐动物数据 | MFA逐动物数据 | 分子箱体统计 | MFA箱体统计 |", "|---|---|---|---|---|"]
    for tissue, settings in config["tissues"].items():
        molecule_link = f"[CSV]({tissue}/{tissue}_molecular_plot_data.csv)" if settings.get("molecular_available", True) else f"[数据覆盖说明]({tissue}/{tissue}_molecular_data_availability.csv)"
        box_link = f"[CSV]({tissue}/{tissue}_molecular_box_statistics.csv)" if settings.get("molecular_available", True) else "不适用"
        lines.append(f"| {settings['chinese']} | {molecule_link} | [CSV]({tissue}/{tissue}_MFA_plot_data.csv) | {box_link} | [CSV]({tissue}/{tissue}_MFA_box_statistics.csv) |")
    lines += ["", "合并表：[全部分子绘图数据](all_tissues_molecular_plot_data.csv)、[全部MFA绘图数据](all_tissues_MFA_plot_data.csv)、[通路证据](all_tissues_pathway_evidence.csv)、[RNA候选表](all_tissues_RNA_candidate_evidence.csv)。", "",
              "## 如何读箱线图和CSV", "",
              "- 分子层：Mazdutide 3只、Saline 4只；MFA：Mazdutide 3只、Saline 2只。各组织内每点一只动物，不把多个组织或多个脂质特征当作新的动物重复。",
              "- 箱体为动物值Q1–Q3，中线为中位数，须采用1.5 IQR；所有动物点均叠加显示，包括Tukey规则下的离群点。n=2/3时箱体仅是小样本描述。",
              "- RNA纵轴为从基因汇总计数计算的log2(CPM+1)，用于展示动物分布；图上的padj仍来自原DESeq2差异模型，并非对CPM重新做t检验。通路标注效应采用未收缩log2FC。",
              "- 单代谢特征纵轴沿用PQN处理后的log2相对信号；类别值先在每只动物内对成员的log2信号求均值。它们都不是绝对浓度。",
              "- MFA每动物值为现有管线过滤后拟合净通量的中位数。箱线图只使用这些动物值；QA竖线为同一动物过滤后拟合点P25–P75。两者使用同一CSV，但表达不同层次的变异。",
              "- `near_one_flag`仅提示估计贴近数值1，需结合模型边界核实；它不是自动剔除规则。负净通量保留原始符号，不把带符号源log2FC当作普通浓度倍数。",
              "- 为防止浮点误差被自动缩放夸大：某MFA节点全部动物接近1时，箱线图纵轴设为0.9–1.1；全部接近0时设为−0.05–0.05。数据不被取整或替换；箱体统计CSV附axis_rule与实际y轴范围。",
              "- `node_id`与`sample_id`标识节点和动物；`value`是实际绘图值；`value_scale`说明尺度；`source_padj`或`fdr_mannwhitney/fdr_ttest`注明统计；`source_file`可追溯源表。MFA另保留`n_fit_points`、`fit_p25`、`fit_p75`、`fit_rel_iqr`。", "",
              "## 逐组织核查结果", ""]
    for tissue, settings in config["tissues"].items():
        lines += [f"### {settings['chinese']} / {tissue}", "", descriptions[tissue], ""]
        if settings.get("molecular_available", True):
            lines += ["| 绘图节点 | 数据层 | 源效应值 | 校正p/q | 表述 |", "|---|---|---:|---:|---|"]
        else:
            lines += ["治疗分子层不适用，详见本组织的数据覆盖说明CSV；这不是0个差异分子的结论。", ""]
        subset = evidence[evidence.tissue.eq(tissue)].set_index(["layer", "node_id"])
        for node in settings["molecular_nodes"]:
            row = subset.loc[(node["layer"], node["id"])]
            direction = {"up": "上升", "down": "下降", "NS": "未达FDR", "NA": "不可报告padj"}[row.display_direction]
            lines.append(f"| {node['label']} | {node['layer']} | {row.effect:+.3f} | {fmt(row.padj)} | {direction} |")
            selected_stats.append({"tissue": tissue, "layer": node["layer"], "node_id": node["id"], "node_label": node["label"],
                                   "effect": row.effect, "effect_scale": row.effect_scale, "padj": row.padj,
                                   "fdr_family": row.fdr_family, "source_file": row.source_file})
        if settings.get("molecular_available", True):
            lines += ["", "效应尺度：RNA为未收缩log2FC；代谢特征为相对信号log2FC；脂质类别为两组每动物平均log2信号之差。", ""]
        s = summary.loc[tissue]
        lines.append(f"本组织MFA共{s.mfa_total_reactions}条反应；源效应量门槛标记{s.mfa_source_effect_size_selected}条；MW和Welch的FDR<0.05均为0条。对应最小q分别为{fmt(s.mfa_min_q_mw)}与{fmt(s.mfa_min_q_welch)}。源标记不能写成显著通量变化。")
        lines.append("")
    lines += ["## 校正范围、核查与解释边界", "",
              f"MFA的原始BH family为同一比较的全部{manifest['mfa_bh_family_size']}个组织×反应，先校正再分组织写表。本次重建整个family核对MW/Welch的BH值，未在每张图的9个节点内重新校正。各组织的单代谢特征BH和七类脂质BH也已交叉核对。详见[source_validation.csv](source_validation.csv)。", "",
              "现有MFA的自动统计路径在小样本下回退到effect_size_only，所以文件中的is_significant可以只是效应量入选。当前六组织没有FDR支持的通量变化；这一点不足以证明无药物效应，只说明这些模型估计尚不能提供强治疗推断。", "",
              "RNA和常规代谢物取自同一常规组学动物集合（Mazdutide 120/123/124；Saline 95/117/119/125）。MFA使用另一动物集合，不能把RNA、代谢物和MFA当作同一只动物的配对测量。", "",
              "心脏肉碱相关类别的七类BH q约0.0216，逐一留出动物后的方向保留；精确MW敏感性p约0.0571，表明小样本显著性仍依赖Welch模型假设。所有器官的方向敏感性表保存在各自目录；RNA留一检查是绘图尺度上的描述，不是DESeq2重新拟合。", "",
              "旧类别名Acylcarnitine实际按Carnitine字符串选择，包含游离肉碱；Ceramide按^Cer选择，包含部分CerP、二氢/植物神经酰胺。类别成员CSV保留完整定义；PA统一写为PA/LPA。", "",
              "生化映射只采用核查过的反应关系。整合用flux_gene_mapping还包含上下游和转运关联，不把所有关联基因都视作该反应的直接催化酶。例如GLUL是谷氨酰胺合成反应催化酶；GPS_c、LIPID_c、FA_usage等是模型聚合节点，不能直接替代GPD2、PC合成或CPT单酶活性。", "",
              "## 脑HFD/LFD参照：独立存放", "",
              "[参照说明](brain/HFD_LFD_reference/README_HFD_LFD参照.md)、[参照通路图](brain/HFD_LFD_reference/brain_reference_pathway.svg)、[参照9节点箱线图](brain/HFD_LFD_reference/brain_reference_molecular_key_nodes_boxplots.svg)、[逐动物CSV](brain/HFD_LFD_reference/brain_reference_molecular_plot_data.csv)、[箱体统计](brain/HFD_LFD_reference/brain_reference_molecular_box_statistics.csv)。", "",
              "该参照来自20250425队列HFD-14W vs LFD-14W；node_id应联合layer区分RNA与同名蛋白。所有记录显式标记is_treatment_comparison=False，并未加入本目录的治疗合并表。", "",
              "## 已核查的生化资料", ""]
    for reference in config["biochemical_references"]:
        lines.append(f"- [{reference['name']}]({reference['url']})")
    lines += ["", "上述资料支撑灰色生化连接；彩色实验方向来自本目录的队列B数据。", "",
              "## 重现与文件版本", "",
              "数据提取与绘图使用配置驱动脚本。输入文件SHA256记录在[source_manifest.json](source_manifest.json)，实际配置记录在[resolved_config.json](resolved_config.json)。本次未重新拟合RNA、重做原始质谱处理或运行MFA求解器。", "", "```powershell",
              f"python -m omicsweave.integration.tissue_pathway_panels --config {manifest['config']}",
              f"python -m omicsweave.integration.plot_tissue_pathway_panels --evidence-dir {config['paths']['output']}",
              "```", ""]
    (directory / "README_六组织通路与箱线图.md").write_text("\n".join(lines), encoding="utf-8")
    (directory / "README_四组织通路与箱线图.md").write_text("# 正式结果已扩展至六组织\n\n请查看[六组织正式结果索引](README_六组织通路与箱线图.md)。原四组织的关键节点及比较口径保留，新增肾皮质、脑皮层MFA和独立脑HFD/LFD参照。\n", encoding="utf-8")
    (directory / "FORMAL_RESULTS.md").write_text("# 正式结果入口\n\n[六组织结果与文件索引](README_六组织通路与箱线图.md)\n\n[完整图件与CSV包](panelB_six_tissues_figures_and_csv.zip)\n\n脑的HFD/LFD分子参照位于brain/HFD_LFD_reference，与队列B的MFA治疗证据分开。旧panel_b_renal_brain_20260909目录仅供追溯参考。\n", encoding="utf-8")
    pd.DataFrame(selected_stats).to_csv(directory / "all_tissues_molecular_key_node_statistics.csv", index=False, encoding="utf-8-sig")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--tissue", action="append", help="只重绘指定组织，可重复")
    args = parser.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    directory = args.evidence_dir.resolve()
    config = json.loads((directory / "resolved_config.json").read_text(encoding="utf-8"))
    for tissue, settings in config["tissues"].items():
        if args.tissue and tissue not in args.tissue:
            continue
        pathway_figure(directory, tissue, settings, config)
        pathway_figure(directory, tissue, settings, config, monochrome=True)
        if settings.get("molecular_available", True):
            boxplot_figure(directory, tissue, settings, config, "molecular")
        boxplot_figure(directory, tissue, settings, config, "MFA")
        mfa_fit_qa(directory, tissue, settings, config)
        print(f"{tissue}: 通路图、可用数据层箱线图、拟合QA及箱线统计CSV已生成。", flush=True)
    if config.get("brain_reference") and (not args.tissue or "brain" in args.tissue):
        plot_brain_reference(directory, config)
    write_report(directory, config)


if __name__ == "__main__":
    main()
