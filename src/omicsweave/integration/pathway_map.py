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

from .dossier import calibration as cal
from .dossier.map_view import coverage_table, map_view, region_summary
from .mappings.network_layout import Layout, LayoutError, default_layout, load_layout, validate_layout

DEFAULTS: dict[str, Any] = {"enabled": True, "layout": "", "tissues": [], "contrasts": [],
                            "highlight": {"min_claim": "L2", "top_rank": 10}, "formats": ["pdf", "png"]}
MAP_TABLE = "pathway_map.parquet"
COVERAGE_TABLE = "pathway_map_coverage.parquet"
FIGURE_DIR = ("figures", "pathway_map")
CLAIM_LEVELS = ("L1", "L2", "L3")


def settings(raw: dict[str, Any] | None) -> dict[str, Any]:
    """参数组 ``pathway_map``：未知键报错，``highlight`` 按项合并。"""
    raw = dict(raw or {})
    unknown = sorted(set(raw) - set(DEFAULTS))
    if unknown:
        raise ValueError(f"pathway_map 中有无法识别的项：{'、'.join(unknown)}（可用：{'、'.join(DEFAULTS)}）")
    highlight = {**DEFAULTS["highlight"], **(raw.get("highlight") or {})}
    if set(highlight) - set(DEFAULTS["highlight"]):
        raise ValueError(f"pathway_map.highlight 只接受 {'、'.join(DEFAULTS['highlight'])}")
    if highlight["min_claim"] not in CLAIM_LEVELS:
        raise ValueError(f"pathway_map.highlight.min_claim 应为 {'、'.join(CLAIM_LEVELS)} 之一")
    return {**DEFAULTS, **raw, "highlight": highlight}


def resolve_layout(value: str, root: Path | None) -> Layout:
    if not value:
        return default_layout()
    path = Path(value)
    if not path.is_absolute() and root is not None:
        path = root / path
    return load_layout(path)


def _sample_sizes(rows: pd.DataFrame) -> str:
    parts = []
    names = {"transcriptomics": "RNA", "proteomics": "蛋白", "metabolomics": "代谢物"}
    for layer, label in names.items():
        part = rows[(rows["layer"] == layer) & rows["n_numerator"].notna()]
        if not part.empty:
            parts.append(f"{label} {int(part['n_numerator'].max())} vs {int(part['n_denominator'].max())}")
    flux = rows[(rows["layer"] == "mfa") & rows["n_numerator"].notna() & (rows["tracer"] != "pooled")]
    for tracer, part in flux.groupby("tracer"):
        parts.append(f"MFA {tracer} {int(part['n_numerator'].max())} vs {int(part['n_denominator'].max())}")
    return "；".join(parts) or "没有样本数信息"


def _fusion(dossier_dir: Path) -> pd.DataFrame | None:
    path = dossier_dir.parent / "pathway_convergence" / "tables" / "pathway_fusion.parquet"
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
    if calibration is None:
        path = dossier_dir / "evaluation" / "assessment_calibration.parquet"
        calibration = pd.read_parquet(path) if path.is_file() else pd.DataFrame()
    layout = resolve_layout(conf["layout"], root)
    problems = validate_layout(layout)
    if problems:
        raise LayoutError(f"布局文件 {layout.path} 未通过校验：" + "；".join(problems[:8])
                          + (f" 等 {len(problems)} 处" if len(problems) > 8 else ""))
    fusion = _fusion(dossier_dir)
    coverage = coverage_table(layout, members, assessed["flux_id"].unique())
    pairs = assessed[["tissue", "contrast_id"]].drop_duplicates().sort_values(["tissue", "contrast_id"])
    if conf["tissues"]:
        pairs = pairs[pairs["tissue"].isin([str(t) for t in conf["tissues"]])]
    if conf["contrasts"]:
        pairs = pairs[pairs["contrast_id"].isin([str(c) for c in conf["contrasts"]])]
    formats = [str(f) for f in conf["formats"]]
    unplaced = coverage[coverage["placement"] == "unplaced"]
    covered = coverage[coverage["placement"] == "covered_by_steps"]
    views = []
    run_text = "；".join(f"{name} {run}" for name, run in runs.items())
    for row in pairs.itertuples():
        view = map_view(layout, members, axes, assessed, row.tissue, row.contrast_id, source, coverage,
                        top_rank=int(conf["highlight"]["top_rank"]), min_claim=str(conf["highlight"]["min_claim"]))
        views.append(view)
        here = members[(members["tissue"] == row.tissue) & (members["contrast_id"] == row.contrast_id)
                       & (members["source"] == source)]
        claims = view.loc[view["flux_claim"].isin(CLAIM_LEVELS), "flux_claim"].value_counts()
        caption = [
            f"来源：{run_text}；样本数（分子组 vs 分母组）：{_sample_sizes(here)}",
            (cal.caption(calibration, row.tissue, row.contrast_id) or "Calibration: not available")
            + "。图上 L1 / L2 / L3：" + "、".join(f"{k} {claims.get(k, 0)}" for k in CLAIM_LEVELS)
            + "。L1 是待检验的预测，引用前先对照诱饵期望。",
            f"布局 {layout.path.name}（sha256 {layout.sha256[:12]}）；评估集中未上图的 ID {len(unplaced)} 个"
            + (f"（{'、'.join(unplaced['flux_id'][:6])}{' 等' if len(unplaced) > 6 else ''}）" if len(unplaced) else "")
            + f"；成员全部被画出的反应覆盖、不另画方框的聚合节点 {len(covered)} 个。",
        ]
        stem = out_dir.joinpath(*FIGURE_DIR, row.tissue, f"map_{row.contrast_id}")
        stem.parent.mkdir(parents=True, exist_ok=True)
        write_excel_tsv(view, stem.with_suffix(".tsv"))
        regions = region_summary(view, layout, fusion, row.tissue, row.contrast_id)
        plot_pathway_map(view, layout, regions, f"图 E　{row.tissue} · {row.contrast_id}　分组织代谢通路总图",
                         caption, stem, formats)
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    table = pd.concat(views, ignore_index=True) if views else pd.DataFrame()
    table.to_parquet(out_dir / "tables" / MAP_TABLE, index=False)
    coverage.to_parquet(out_dir / "tables" / COVERAGE_TABLE, index=False)
    notes = [f"通路总图（图 E）：{len(views)} 页，在 {'/'.join(FIGURE_DIR)}/；布局 {layout.path.name}"
             f"（sha256 {layout.sha256[:12]}）"]
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
    notes = render_maps(dossier_dir, out_dir, raw_settings, source, {"整合": bundle.name}, members=members, root=root)
    return out_dir, notes
