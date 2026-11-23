#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通路总图（图 E，G18 M3）的出图流程：读反应档案的表，按组织 × 比较生成视图、图与导出表。

由两处调用：``reaction_dossier`` 模式在评估与校准之后调用 :func:`render_maps`；``omicsweave map`` 从已有的整合结果包
重画（:func:`redraw`），只读表，不改动原结果包。

输出（相对 ``out_dir``）：

- ``figures/pathway_map/<组织>/map_<比较>.{pdf,png}``：图 E；同名 ``.tsv`` 为该页的视图（图只画这张表）；
- ``tables/pathway_map.parquet``：全部页的视图；``tables/pathway_map_coverage.parquet``：评估集中每个 ID 在总图上的去向。
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from omicsweave.contract.export import write_excel_tsv
from omicsweave.mfa.tracers import format_tracer_display
from omicsweave.visualization.graphical_abstract import GA_STYLES
from omicsweave.visualization.palettes import VALID_COLOR_MODES
from omicsweave.visualization.theme import list_presets

from .dossier import calibration as cal
from .dossier.map_view import (
    coverage_table,
    extra_pathway_fusion,
    map_view,
    region_summary,
    signed_module_summary,
)
from .mappings.network_layout import Layout, LayoutError, default_layout, load_layout, validate_layout
from .plotting.pathway_submaps import ALL_SUBMAP_IDS, VALID_MIN_SIGNALS, select_submaps_for_view

DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "layout": "",
    "tissues": [],
    "contrasts": [],
    "highlight": {"min_claim": "L2", "top_rank": 10},
    "formats": ["pdf", "png", "html"],
    "mode": "atlas",
    "style": "metro_atlas",
    "color_mode": "semantic",
    "strip_mode": "omics",
    "omics_edge_color": "tinted",
    "regions": [],
    "submaps": {
        "enabled": True,
        "panels": [],
        "min_signal": "assessed",
        "include_optional_active": True,
    },
    "abstract": {
        "enabled": True,
        "style": "",
        "color_mode": "",
    },
}
MAP_TABLE = "pathway_map.parquet"
COVERAGE_TABLE = "pathway_map_coverage.parquet"
FIGURE_DIR = ("figures", "pathway_map")
SUBMAP_SUBDIR = "submaps"
CLAIM_LEVELS = ("L1", "L2", "L3")
MAP_MODES = ("atlas", "publication")


def settings(raw: dict[str, Any] | None) -> dict[str, Any]:
    """参数组 ``pathway_map``：未知键报错，``highlight`` / ``submaps`` / ``abstract`` 按项合并。"""
    raw = dict(raw or {})
    unknown = sorted(set(raw) - set(DEFAULTS))
    if unknown:
        raise ValueError(f"pathway_map 中有无法识别的项：{'、'.join(unknown)}（可用：{'、'.join(DEFAULTS)}）")
    highlight = {**DEFAULTS["highlight"], **(raw.get("highlight") or {})}
    if set(highlight) - set(DEFAULTS["highlight"]):
        raise ValueError(f"pathway_map.highlight 只接受 {'、'.join(DEFAULTS['highlight'])}")
    if highlight["min_claim"] not in CLAIM_LEVELS:
        raise ValueError(f"pathway_map.highlight.min_claim 应为 {'、'.join(CLAIM_LEVELS)} 之一")

    raw_sub = raw.get("submaps")
    if isinstance(raw_sub, bool):
        submaps = {**DEFAULTS["submaps"], "enabled": raw_sub}
    else:
        submaps = {**DEFAULTS["submaps"], **(raw_sub or {})}
    if set(submaps) - set(DEFAULTS["submaps"]):
        raise ValueError(f"pathway_map.submaps 只接受 {'、'.join(DEFAULTS['submaps'])}")
    if submaps["min_signal"] not in VALID_MIN_SIGNALS:
        raise ValueError(f"pathway_map.submaps.min_signal 应为 {'、'.join(VALID_MIN_SIGNALS)} 之一")
    for pid in submaps["panels"] or ():
        if str(pid).strip().lower() != "all" and str(pid).strip() not in ALL_SUBMAP_IDS:
            raise ValueError(f"pathway_map.submaps.panels 包含未知子图 {pid!r}（可用：{'、'.join(ALL_SUBMAP_IDS)} 或 all）")

    raw_abs = raw.get("abstract")
    if isinstance(raw_abs, bool):
        abstract = {**DEFAULTS["abstract"], "enabled": raw_abs}
    else:
        abstract = {**DEFAULTS["abstract"], **(raw_abs or {})}
    if set(abstract) - set(DEFAULTS["abstract"]):
        raise ValueError(f"pathway_map.abstract 只接受 {'、'.join(DEFAULTS['abstract'])}")
    if abstract["style"] and abstract["style"] not in GA_STYLES:
        raise ValueError(f"pathway_map.abstract.style 应为 {'、'.join(GA_STYLES)} 之一，实际为 {abstract['style']!r}")
    if abstract["color_mode"] and abstract["color_mode"] not in VALID_COLOR_MODES:
        raise ValueError(f"pathway_map.abstract.color_mode 应为 {'、'.join(VALID_COLOR_MODES)} 之一")

    merged = {**DEFAULTS, **raw, "highlight": highlight, "submaps": submaps, "abstract": abstract}
    if merged["mode"] not in MAP_MODES:
        raise ValueError(f"pathway_map.mode 应为 {'、'.join(MAP_MODES)} 之一，实际为 {merged['mode']!r}")
    if merged["style"] not in list_presets():
        raise ValueError(f"pathway_map.style 应为 {'、'.join(list_presets())} 之一，实际为 {merged['style']!r}")
    if merged["color_mode"] not in VALID_COLOR_MODES:
        raise ValueError(f"pathway_map.color_mode 应为 {'、'.join(VALID_COLOR_MODES)} 之一，实际为 {merged['color_mode']!r}")
    if merged["strip_mode"] not in ("omics", "axes"):
        raise ValueError(f"pathway_map.strip_mode 应为 omics、axes 之一，实际为 {merged['strip_mode']!r}")
    if merged["omics_edge_color"] not in ("tinted", "grey"):
        raise ValueError(
            f"pathway_map.omics_edge_color 应为 tinted、grey 之一，实际为 {merged['omics_edge_color']!r}"
        )
    return merged


def resolve_layout(value: str, root: Path | None) -> Layout:
    if not value:
        return default_layout()
    path = Path(value)
    if not path.is_absolute() and root is not None:
        path = root / path
    return load_layout(path)


def _sample_sizes(rows: pd.DataFrame) -> str:
    parts = []
    names = {"transcriptomics": "RNA", "proteomics": "Protein", "metabolomics": "Metabolite"}
    for layer, label in names.items():
        part = rows[(rows["layer"] == layer) & rows["n_numerator"].notna()]
        if not part.empty:
            parts.append(f"{label} {int(part['n_numerator'].max())} vs {int(part['n_denominator'].max())}")
    flux = rows[(rows["layer"] == "mfa") & rows["n_numerator"].notna() & (rows["tracer"] != "pooled")]
    for tracer, part in flux.groupby("tracer"):
        t_disp = format_tracer_display(str(tracer))
        parts.append(f"MFA ({t_disp}) {int(part['n_numerator'].max())} vs {int(part['n_denominator'].max())}")
    return "; ".join(parts) or "N/A"


def _english_run_label(name: str) -> str:
    text = str(name).strip()
    if text == "整合":
        return "Integration"
    if text.startswith("来源 "):
        return "Source " + text[len("来源 "):]
    if text.startswith("来源"):
        return "Source " + text[len("来源"):].strip()
    return text


def _fusion(dossier_dir: Path) -> pd.DataFrame | None:
    path = dossier_dir.parent / "pathway_convergence" / "tables" / "pathway_fusion.parquet"
    return pd.read_parquet(path) if path.is_file() else None


def _pathway_tests(dossier_dir: Path) -> pd.DataFrame | None:
    path = dossier_dir / "tables" / "pathway_tests.parquet"
    return pd.read_parquet(path) if path.is_file() else None


def render_maps(dossier_dir: Path, out_dir: Path, raw_settings: dict[str, Any] | None, source: str,
                runs: dict[str, str], members: pd.DataFrame | None = None, axes: pd.DataFrame | None = None,
                assessed: pd.DataFrame | None = None, calibration: pd.DataFrame | None = None,
                root: Path | None = None) -> list[str]:
    """画全部（组织, 比较）的图 E，返回运行说明。``dossier_dir``：``integration/reaction_dossier``（读表与通路融合表）；
    ``runs``：写入图注的 {名称: 运行 ID}；``root``：解析相对路径的布局文件（实验目录）。"""
    from .plotting.pathway_map_plots import plot_pathway_map

    conf = settings(raw_settings)
    if not conf["enabled"]:
        return ["通路总图（图 E）：pathway_map.enabled 为 false，未出图"]
    tables = dossier_dir / "tables"
    members = pd.read_parquet(tables / "dossier_members.parquet") if members is None else members
    axes = pd.read_parquet(tables / "axis_assessment.parquet") if axes is None else axes
    assessed = pd.read_parquet(tables / "reaction_assessment.parquet") if assessed is None else assessed
    samples_path = tables / "dossier_samples.parquet"
    samples = pd.read_parquet(samples_path) if samples_path.is_file() else None
    sc_path = tables / "dossier_sample_contrasts.parquet"
    sample_contrasts = pd.read_parquet(sc_path) if sc_path.is_file() else None
    if calibration is None:
        path = dossier_dir / "evaluation" / "assessment_calibration.parquet"
        calibration = pd.read_parquet(path) if path.is_file() else pd.DataFrame()
    layout = resolve_layout(conf["layout"], root)
    problems = validate_layout(layout)
    if problems:
        raise LayoutError(f"布局文件 {layout.path} 未通过校验：" + "；".join(problems[:8])
                          + (f" 等 {len(problems)} 处" if len(problems) > 8 else ""))
    fusion = _fusion(dossier_dir)
    pathway_tests = _pathway_tests(dossier_dir)
    coverage = coverage_table(layout, members, assessed["flux_id"].unique())
    pairs = assessed[["tissue", "contrast_id"]].drop_duplicates().sort_values(["tissue", "contrast_id"])
    if conf["tissues"]:
        pairs = pairs[pairs["tissue"].isin([str(t) for t in conf["tissues"]])]
    if conf["contrasts"]:
        pairs = pairs[pairs["contrast_id"].isin([str(c) for c in conf["contrasts"]])]
    formats = [str(f) for f in conf["formats"]]
    static_formats = [f for f in formats if f.lstrip(".").lower() != "html"]
    unplaced = coverage[coverage["placement"] == "unplaced"]
    covered = coverage[coverage["placement"] == "covered_by_steps"]
    views = []
    n_submaps_rendered = 0
    skipped_submaps: list[str] = []
    n_abstracts_rendered = 0
    sub_conf = conf["submaps"]
    abs_conf = conf["abstract"]
    run_text = "; ".join(f"{_english_run_label(name)} {run}" for name, run in runs.items())
    for row in pairs.itertuples():
        view = map_view(layout, members, axes, assessed, row.tissue, row.contrast_id, source, coverage,
                        top_rank=int(conf["highlight"]["top_rank"]), min_claim=str(conf["highlight"]["min_claim"]))
        views.append(view)
        here = members[(members["tissue"] == row.tissue) & (members["contrast_id"] == row.contrast_id)
                       & (members["source"] == source)]
        from .plotting.pathway_map_plots import _independent_flux_chains
        claims = view.loc[view["flux_claim"].isin(CLAIM_LEVELS), "flux_claim"].value_counts()
        n_l23 = int(claims.get("L2", 0)) + int(claims.get("L3", 0))
        n_indep = _independent_flux_chains(view)
        indep_note = f" ({n_indep} indep. chains)" if (n_l23 > 0 and n_indep < n_l23) else ""
        caption = [
            f"Run: {run_text}  |  Sample sizes (num vs den): {_sample_sizes(here)}",
            (cal.caption(calibration, row.tissue, row.contrast_id) or "Calibration: not available")
            + f"  |  On-map claims: Multi-tracer={claims.get('L3', 0)}, Flux={claims.get('L2', 0)}{indep_note}, Enz-pred={claims.get('L1', 0)}"
            + "  (Enz-pred = enzyme-only; verify against decoy baseline before citing).",
            f"Layout: {layout.path.name} (sha256 {layout.sha256[:12]})  |  Unplaced assessed IDs: {len(unplaced)}"
            + (f" ({', '.join(unplaced['flux_id'][:6])}{', ...' if len(unplaced) > 6 else ''})" if len(unplaced) else "")
            + f"  |  Step-covered aggregates folded onto edges: {len(covered)}.",
        ]
        stem = out_dir.joinpath(*FIGURE_DIR, row.tissue, f"map_{row.contrast_id}")
        stem.parent.mkdir(parents=True, exist_ok=True)
        write_excel_tsv(view, stem.with_suffix(".tsv"))
        regions = region_summary(view, layout, fusion, row.tissue, row.contrast_id)
        extra_fusion = extra_pathway_fusion(fusion, layout, row.tissue, row.contrast_id)
        signed_mods = signed_module_summary(pathway_tests, row.tissue, row.contrast_id)
        tissue_title = str(row.tissue).replace("_", " ").upper() if len(str(row.tissue)) <= 3 else str(row.tissue).replace("_", " ").title()
        plot_pathway_map(view, layout, regions,
                         f"{tissue_title} · {row.contrast_id}  —  Tissue Multi-Omics Pathway Atlas",
                         caption, stem, formats, extra_fusion=extra_fusion, signed_modules=signed_mods,
                         top_rank=int(conf["highlight"]["top_rank"]),
                         mode=str(conf["mode"]),
                         style=str(conf["style"]),
                         color_mode=str(conf["color_mode"]),
                         strip_mode=str(conf["strip_mode"]),
                         omics_edge_color=str(conf["omics_edge_color"]),
                         focus_regions=list(conf["regions"]) or None,
                         members=members,
                         samples=samples,
                         sample_contrasts=sample_contrasts,
                         tissue=str(row.tissue),
                         contrast_id=str(row.contrast_id),
                         source=str(source))

        # 默认生成正文分专题代谢子通路图（沿用与主通路图相同的风格与色彩语义，无数据时自动跳过；静态格式）
        if sub_conf["enabled"] and static_formats:
            selected_subs, skipped = select_submaps_for_view(
                view, layout,
                panels=list(sub_conf["panels"] or ()),
                min_signal=str(sub_conf["min_signal"]),
                include_optional_active=bool(sub_conf["include_optional_active"]),
            )
            sub_dir = stem.parent / SUBMAP_SUBDIR
            if selected_subs:
                sub_dir.mkdir(parents=True, exist_ok=True)
            for sub in selected_subs:
                sub_stem = sub_dir / f"{row.contrast_id}_{sub.id}"
                sub_title = f"{tissue_title} · {row.contrast_id}  —  {sub.title}"
                sub_note = f"Sub-pathway panel ({sub.subtitle}) · Cross-pathway links attach at boundary stations; full network in Pathway Atlas."
                plot_pathway_map(
                    view, layout, regions,
                    sub_title, caption[:2], sub_stem, static_formats,
                    top_rank=int(conf["highlight"]["top_rank"]),
                    mode="publication",
                    style=str(conf["style"]),
                    color_mode=str(conf["color_mode"]),
                    strip_mode=str(conf["strip_mode"]),
                    omics_edge_color=str(conf["omics_edge_color"]),
                    focus_regions=list(sub.regions),
                    include_reactions=list(sub.include_reactions) or None,
                    exclude_reactions=list(sub.exclude_reactions) or None,
                    panel_note=sub_note,
                )
                n_submaps_rendered += 1
            for sid, reason in skipped:
                if reason != "optional_inactive":
                    skipped_submaps.append(f"{row.tissue}/{row.contrast_id}:{sid}({reason})")

        # 默认同步生成小规模 Graphical Abstract（继承主通路图风格；slide_dark 回退到 metro_atlas）
        if abs_conf["enabled"] and static_formats:
            from omicsweave.visualization.graphical_abstract import (
                build_abstract_from_pathway_view,
                render_graphical_abstract,
            )

            ga_style = str(abs_conf["style"]) or (str(conf["style"]) if str(conf["style"]) in GA_STYLES else "metro_atlas")
            ga_color_mode = str(abs_conf["color_mode"]) or str(conf["color_mode"])
            ga_spec = build_abstract_from_pathway_view(
                view, regions,
                tissue=str(row.tissue),
                contrast=str(row.contrast_id),
                provenance=f"Run: {run_text}  |  Style: {ga_style} ({ga_color_mode})  |  Auto-aggregated from Pathway Atlas view (<=10 hubs)",
            )
            ga_stem = stem.parent / f"abstract_{row.contrast_id}"
            render_graphical_abstract(ga_spec, ga_stem, style=ga_style, color_mode=ga_color_mode, formats=static_formats)
            n_abstracts_rendered += 1

    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    table = pd.concat(views, ignore_index=True) if views else pd.DataFrame()
    table.to_parquet(out_dir / "tables" / MAP_TABLE, index=False)
    coverage.to_parquet(out_dir / "tables" / COVERAGE_TABLE, index=False)
    notes = [f"通路总图（图 E）：{len(views)} 页，在 {'/'.join(FIGURE_DIR)}/；风格 {conf['style']}（{conf['color_mode']}）；布局 {layout.path.name}"
             f"（sha256 {layout.sha256[:12]}）"]
    if n_submaps_rendered:
        notes.append(f"代谢子通路图：{n_submaps_rendered} 张，在 {'/'.join(FIGURE_DIR)}/<组织>/{SUBMAP_SUBDIR}/")
    if skipped_submaps:
        preview = "、".join(skipped_submaps[:6]) + (f" 等 {len(skipped_submaps)} 项" if len(skipped_submaps) > 6 else "")
        notes.append(f"代谢子通路图（因无评估数据跳过）：{preview}")
    if n_abstracts_rendered:
        notes.append(f"概念摘要图（Graphical Abstract）：{n_abstracts_rendered} 张，在 {'/'.join(FIGURE_DIR)}/<组织>/abstract_*.pdf")
    if len(unplaced):
        notes.append(f"通路总图：评估集中 {len(unplaced)} 个 ID 在布局中没有位置：{'、'.join(unplaced['flux_id'])}")
    return notes


def redraw(bundle: Path, out_dir: Path | None, raw_settings: dict[str, Any] | None, root: Path | None) -> tuple[Path, list[str]]:
    """从已有整合结果包重画图 E，写到新目录（默认在结果包旁的 ``_redraws/<运行 ID>/pathway_map_<时间>/``），不改动原结果包。"""
    dossier_dir = bundle / "integration" / "reaction_dossier"
    if not (dossier_dir / "tables" / "reaction_assessment.parquet").is_file():
        raise FileNotFoundError(f"{bundle} 没有反应档案的评估表（integration/reaction_dossier/tables/reaction_assessment.parquet）")
    members = pd.read_parquet(dossier_dir / "tables" / "dossier_members.parquet")
    sources = list(dict.fromkeys(members["source"].dropna()))
    with_flux = members[(members["layer"] == "mfa") & ~members["cell_state"].isin(["not_measured", "not_infused"])]
    source = with_flux["source"].mode().iloc[0] if not with_flux.empty else sources[0]
    if out_dir is None:
        out_dir = bundle.parent / "_redraws" / bundle.name / f"pathway_map_{datetime.now():%Y%m%d_%H%M%S}"
    notes = render_maps(dossier_dir, out_dir, raw_settings, source, {"Integration": bundle.name}, members=members, root=root)
    return out_dir, notes


def render_abstracts(
    bundle: Path,
    out_dir: Path | None,
    *,
    tissues: list[str] | None = None,
    contrasts: list[str] | None = None,
    style: str = "metro_atlas",
    color_mode: str = "semantic",
    formats: list[str] | None = None,
    layout_path: str = "",
    root: Path | None = None,
) -> tuple[Path, list[Path]]:
    """从已有整合结果包自动聚合生成小规模 Graphical Abstract（≤ 10 核心枢纽），不改动原结果包。"""
    from omicsweave.visualization.graphical_abstract import (
        build_abstract_from_pathway_view,
        render_graphical_abstract,
    )

    dossier_dir = bundle / "integration" / "reaction_dossier"
    tables = dossier_dir / "tables"
    if not (tables / "reaction_assessment.parquet").is_file():
        raise FileNotFoundError(f"{bundle} 没有反应档案的评估表（integration/reaction_dossier/tables/reaction_assessment.parquet）")
    members = pd.read_parquet(tables / "dossier_members.parquet")
    axes = pd.read_parquet(tables / "axis_assessment.parquet")
    assessed = pd.read_parquet(tables / "reaction_assessment.parquet")
    fusion = _fusion(dossier_dir)
    layout = resolve_layout(layout_path, root)
    coverage = coverage_table(layout, members, assessed["flux_id"].unique())

    sources = list(dict.fromkeys(members["source"].dropna()))
    with_flux = members[(members["layer"] == "mfa") & ~members["cell_state"].isin(["not_measured", "not_infused"])]
    source = with_flux["source"].mode().iloc[0] if not with_flux.empty else sources[0]

    pairs = assessed[["tissue", "contrast_id"]].drop_duplicates().sort_values(["tissue", "contrast_id"])
    if tissues:
        pairs = pairs[pairs["tissue"].isin([str(t) for t in tissues])]
    if contrasts:
        pairs = pairs[pairs["contrast_id"].isin([str(c) for c in contrasts])]

    target_dir = out_dir or (bundle.parent / "_redraws" / bundle.name / f"graphical_abstract_{style}_{datetime.now():%Y%m%d_%H%M%S}")
    fmts = list(formats or ["pdf", "png"])
    written: list[Path] = []
    for row in pairs.itertuples():
        view = map_view(layout, members, axes, assessed, row.tissue, row.contrast_id, source, coverage)
        reg = region_summary(view, layout, fusion, row.tissue, row.contrast_id)
        spec = build_abstract_from_pathway_view(
            view, reg,
            tissue=str(row.tissue),
            contrast=str(row.contrast_id),
            provenance=f"Source: {bundle.name}  |  Style: {style}  |  Auto-aggregated from Pathway Atlas view (<=10 hubs)",
        )
        stem = target_dir / str(row.tissue) / f"graphical_abstract_{row.contrast_id}_{style}"
        written.extend(render_graphical_abstract(spec, stem, style=style, color_mode=color_mode, formats=fmts))
    return target_dir, written
