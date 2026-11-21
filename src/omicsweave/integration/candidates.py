#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
无面板候选表：只读整合结果包中已有的全局表，生成各组织 × 比较的候选清单，不依赖 ``panels.yaml`` 选了什么。

反应（``reaction_dossier`` 的 ``reaction_assessment``，全部网络反应 × 组织 × 比较）：

- 全量表逐行保留；简表取通量结论 L1 及以上的全部，再加每个组织 × 比较排名前 :data:`REACTION_TOP_N` 的反应
  （``selected_by`` 说明入选原因）；
- ``panels`` 列标注该行属于哪些研究者面板（面板声明的组织与比较内），只作标注，不影响入选；
- ``cell_L1`` / ``decoy_L1_expected``：该组织 × 比较的 L1 数与诱饵校准下的期望（``evaluation/
  assessment_calibration.parquet``；未运行校准时为空）。两者接近时，该格的 L1 与随机水平相当；
- ``qc_status`` / ``group_sizes``：各层的样本 QC 状态（``qc/contrast_reliability.parquet``）与可用样本数。

通路（``pathway_convergence`` 的证据表与融合表），每个（组织, 比较, 规范化通路名）一行：

- 各层取该通路最显著的一条记录，列出 padj、方向、库、方法与原名；``direction_agreement`` 的规则与融合表相同；
- ``support``：``fused``（融合 FDR < :data:`FUSED_ALPHA`）、``multi_layer``（未达融合阈值但至少两层各自显著）、
  ``single_layer``（只有一层显著）。各层的显著阈值取整合运行清单中的 ``ora_padj_cutoff`` / ``gsea_padj_cutoff``；
- 全量表取融合显著或至少一层显著的通路；简表取 ``fused`` 与 ``multi_layer`` 的全部，再加每个（组织, 比较, 组学）
  单层显著中 padj 最小的 :data:`PATHWAY_TOP_N` 条；
- 跨层只按规范化名称匹配（``matched_by``），融合值用于探索排序，不是经验证的发现错误率。

同一结果包换一套面板，候选行与各列逐值相同，只有 ``panels`` 列不同。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np
import pandas as pd

from omicsweave.contract.export import write_excel_tsv
from omicsweave.core.tissues import merge_tissue_aliases, resolve_tissue
from omicsweave.contract.manifest import read_manifest
from omicsweave.core.provenance import sha256_file
from omicsweave.core.reliability import QC_DIR, RELIABILITY_FILE

from .dossier.panels import Panel, PanelError, load_panel_selection, load_panels
from .modes.pathway_convergence import DEFAULTS as CONVERGENCE_DEFAULTS
from .modes.pathway_convergence import direction_agreement

DOSSIER = "reaction_dossier"
CONVERGENCE = "pathway_convergence"
REACTION_TOP_N = 10
PATHWAY_TOP_N = 5
FUSED_ALPHA = 0.05
WORKBOOK = "candidates.xlsx"
SELECTION_LABELS = {"pre_specified": "预设（看结果之前定）", "post_hoc": "看过结果后选定", "unspecified": "未声明"}
REACTION_FRONT = ["tissue", "contrast_id", "rank", "flux_id", "flux_claim", "direction", "reaction_class",
                  "n_concordant_axes", "discordant_axes", "mechanism_class", "conjunction_statement", "conjunction_padj",
                  "any_pvalue", "completeness", "completeness_by_layer", "testable_by", "phenomena", "panels",
                  "selected_by", "cell_L1", "decoy_L1_expected", "decoy_L1_sd", "decoy_replicates", "qc_status",
                  "group_sizes"]
PATHWAY_COLUMNS = ["tissue", "contrast_id", "contrast_role", "pathway_term", "support", "fused_fdr", "fused_pvalue",
                   "n_layers_tested", "layers_tested", "n_layers_significant", "significant_layers", "best_padj",
                   "layer_padj", "layer_directions", "direction_agreement", "libraries", "methods", "source_terms",
                   "matched_by"]
SUPPORT_ORDER = {"fused": 0, "multi_layer": 1, "single_layer": 2}


@dataclass
class PanelInfo:
    """整合时使用的面板声明：文件（相对项目目录）、选定方式与依据；文件已改动时不再标注面板。"""
    file: str | None = None
    sha256: str | None = None
    selection: str = "unspecified"
    basis: str = ""
    panels: dict[str, Panel] = field(default_factory=dict)
    note: str = ""

    def record(self) -> dict[str, Any]:
        return {"file": self.file, "sha256": self.sha256, "selection": self.selection, "basis": self.basis,
                "panels": sorted(self.panels), "note": self.note}


@dataclass
class Candidates:
    """写出的候选文件与按组织 × 比较的概览。"""
    directory: Path
    overview: pd.DataFrame
    panels: PanelInfo | None
    files: list[Path] = field(default_factory=list)


def _bundle_tissue_aliases(integration_dir: Path, manifest) -> dict[str, str]:
    design_aliases: dict[str, str] = {}
    exp_json = integration_dir / "metadata" / "experiment.json"
    if exp_json.is_file():
        try:
            raw = json.loads(exp_json.read_text(encoding="utf-8"))
            design_aliases = dict(((raw.get("design") or {}).get("tissue_aliases")) or {})
        except (OSError, ValueError, TypeError):
            design_aliases = {}
    dossier_map = dict(((manifest.parameters.get(DOSSIER) or {}).get("tissue_map")) or {})
    try:
        return merge_tissue_aliases(design_aliases, dossier_map)
    except ValueError:
        return design_aliases


def panel_info(integration_dir: Path, project_root: Path) -> PanelInfo | None:
    """整合运行清单登记的 panels.yaml；内容与登记的校验值不同时不读取面板，只记一条说明。"""
    manifest = read_manifest(integration_dir)
    spec = (manifest.parameters.get(DOSSIER) or {}).get("panels")
    if not spec:
        return None
    suffix = PurePosixPath(str(spec).replace("\\", "/"))
    records = [f for f in manifest.files if f.role == "input" and PurePosixPath(f.path).parts[-len(suffix.parts):]
               == suffix.parts]
    if not records:
        return PanelInfo(note=f"运行清单中没有登记面板声明 {spec}")
    record = records[0]
    path = project_root / record.path
    info = PanelInfo(file=record.path, sha256=record.sha256)
    if not path.is_file() or sha256_file(path) != record.sha256:
        info.note = f"{record.path} 在整合之后改过或已不存在，panels 列未标注"
        return info
    try:
        info.panels = load_panels(path)
        aliases = _bundle_tissue_aliases(integration_dir, manifest)
        if aliases:
            for panel in info.panels.values():
                if panel.tissues:
                    panel.tissues = list(dict.fromkeys(resolve_tissue(t, aliases) for t in panel.tissues))
        selection = load_panel_selection(path)
    except PanelError as exc:
        info.note = f"面板声明无法读取：{exc}"
        return info
    info.selection, info.basis = selection["mode"], selection["basis"]
    return info


def _panel_column(table: pd.DataFrame, panels: dict[str, Panel]) -> pd.Series:
    labels = pd.Series([[] for _ in range(len(table))], index=table.index)
    for panel_id, panel in sorted(panels.items()):
        mask = table["flux_id"].isin(panel.reactions)
        if panel.tissues:
            mask &= table["tissue"].isin(panel.tissues)
        if panel.contrasts:
            mask &= table["contrast_id"].isin(panel.contrasts)
        for index in table.index[mask]:
            labels[index].append(panel_id)
    return labels.map(", ".join)


def _qc(base: Path) -> pd.DataFrame:
    path = base / QC_DIR / RELIABILITY_FILE
    columns = ["tissue", "contrast_id", "qc_status", "group_sizes"]
    if not path.is_file():
        return pd.DataFrame(columns=columns)
    table = pd.read_parquet(path).sort_values("omics")
    rows = []
    for (tissue, contrast), group in table.groupby(["tissue", "contrast_id"], sort=True):
        flagged = [f"{r.omics}: {r.status}" + (f"（{r.reason}）" if r.reason else "")
                   for r in group.itertuples() if r.status != "ok"]
        sizes = [f"{r.omics} {r.n_usable_numerator}v{r.n_usable_denominator}" for r in group.itertuples()]
        rows.append({"tissue": tissue, "contrast_id": contrast, "qc_status": "; ".join(flagged) or "ok",
                     "group_sizes": "; ".join(sizes)})
    return pd.DataFrame(rows, columns=columns)


def _decoys(base: Path) -> pd.DataFrame:
    columns = ["tissue", "contrast_id", "cell_L1", "decoy_L1_expected", "decoy_L1_sd", "decoy_replicates"]
    path = base / "evaluation" / "assessment_calibration.parquet"
    if not path.is_file():
        return pd.DataFrame(columns=columns)
    table = pd.read_parquet(path)
    table = table[(table["method"] == "decoy") & (table["metric"] == "L1") & (table["tissue"].fillna("") != "")]
    return table.rename(columns={"target": "cell_L1", "decoy_expected": "decoy_L1_expected",
                                 "decoy_sd": "decoy_L1_sd", "n_replicates": "decoy_replicates"})[columns]


def reaction_candidates(integration_dir: Path, panels: dict[str, Panel] | None = None,
                        top_n: int = REACTION_TOP_N) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """(全量表, 简表)；结果包没有反应档案时为 None。"""
    base = integration_dir / "integration" / DOSSIER
    path = base / "tables" / "reaction_assessment.parquet"
    if not path.is_file():
        return None
    table = pd.read_parquet(path)
    table["direction"] = table["claim_direction"].map({1: "up", -1: "down"}).fillna("")
    table["panels"] = _panel_column(table, panels or {})
    claimed = table["flux_claim"].fillna("L0") != "L0"
    top = table["rank"] <= top_n
    table["selected_by"] = np.select([claimed & top, claimed, top], ["claim, top_n", "claim", "top_n"], "")
    table = table.merge(_decoys(base), on=["tissue", "contrast_id"], how="left")
    table = table.merge(_qc(base), on=["tissue", "contrast_id"], how="left")
    table = table.sort_values(["tissue", "contrast_id", "rank"], kind="stable").reset_index(drop=True)
    front = [c for c in REACTION_FRONT if c in table.columns]   # 旧结果包可能缺少个别列
    table = table[[*front, *[c for c in table.columns if c not in front]]]
    return table, table[table["selected_by"] != ""].reset_index(drop=True)


def _joined(group: pd.DataFrame, column: str, fmt=str) -> str:
    return "; ".join(f"{o}:{fmt(v)}" for o, v in zip(group["omics"], group[column]))


def pathway_candidates(integration_dir: Path, top_n: int = PATHWAY_TOP_N) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """(全量表, 简表)；结果包没有通路收敛结果时为 None。"""
    tables = integration_dir / "integration" / CONVERGENCE / "tables"
    if not (tables / "pathway_evidence.parquet").is_file():
        return None
    evidence = pd.read_parquet(tables / "pathway_evidence.parquet")
    fusion = pd.read_parquet(tables / "pathway_fusion.parquet")
    params = {**CONVERGENCE_DEFAULTS, **(read_manifest(integration_dir).parameters.get(CONVERGENCE) or {})}
    keys = ["tissue", "comparison", "pathway_term"]
    best = (evidence.sort_values("probability", na_position="last", kind="stable")
            .drop_duplicates([*keys, "omics"]).sort_values([*keys, "omics"]))
    cutoff = np.where(best["method"] == "gsea_prerank", float(params["gsea_padj_cutoff"]),
                      float(params["ora_padj_cutoff"]))
    best["significant"] = best["probability"] < cutoff
    fused = fusion[fusion["fused_fdr"] < FUSED_ALPHA]
    wanted = (set(map(tuple, best.loc[best["significant"], keys].to_numpy()))
              | set(map(tuple, fused[keys].to_numpy())))
    best = best[[k in wanted for k in map(tuple, best[keys].to_numpy())]]
    if best.empty:
        empty = pd.DataFrame(columns=PATHWAY_COLUMNS)
        return empty, empty

    rows = []
    for key, group in best.groupby(keys, sort=True):
        significant = group[group["significant"]]
        rows.append({
            **dict(zip(keys, key)), "contrast_role": group["contrast_role"].iloc[0],
            "n_layers_tested": len(group), "layers_tested": ", ".join(group["omics"]),
            "n_layers_significant": len(significant), "significant_layers": ", ".join(significant["omics"]),
            "best_padj": significant["probability"].min() if len(significant) else group["probability"].min(),
            "layer_padj": _joined(group, "probability", lambda v: f"{v:.3g}"),
            "layer_directions": _joined(group, "effect_direction"),
            "direction_agreement": direction_agreement(list(group["effect_direction"])),
            "libraries": _joined(group, "gene_set_library"), "methods": _joined(group, "method"),
            "source_terms": _joined(group, "term_name"), "matched_by": "normalized_name",
        })
    table = pd.DataFrame(rows).merge(fusion[[*keys, "fused_pvalue", "fused_fdr"]], on=keys, how="left")
    table["support"] = np.select(
        [table["fused_fdr"] < FUSED_ALPHA, table["n_layers_significant"] >= 2], ["fused", "multi_layer"],
        "single_layer")
    table = table.rename(columns={"comparison": "contrast_id"})
    table = (table.assign(_order=table["support"].map(SUPPORT_ORDER))
             .sort_values(["tissue", "contrast_id", "_order", "fused_fdr", "best_padj"], na_position="last",
                          kind="stable")
             .drop(columns="_order").reset_index(drop=True)[PATHWAY_COLUMNS])
    single = table[table["support"] == "single_layer"]
    leading = single.groupby(["tissue", "contrast_id", "significant_layers"], sort=False).head(top_n).index
    short = table[(table["support"] != "single_layer") | table.index.isin(leading)].reset_index(drop=True)
    return table, short


def _overview(reactions: tuple | None, pathways: tuple | None) -> pd.DataFrame:
    keys = ["tissue", "contrast_id"]
    parts = []
    if reactions is not None:
        full, short = reactions
        claims = full.assign(L1=full["flux_claim"] == "L1", L2_up=full["flux_claim"].isin(["L2", "L3"]))
        summary = claims.groupby(keys).agg(reactions_L1=("L1", "sum"), reactions_L2_or_higher=("L2_up", "sum"),
                                           decoy_L1_expected=("decoy_L1_expected", "first"))
        summary["reactions_in_short_list"] = short.groupby(keys).size()
        summary["short_list_in_panels"] = short[short["panels"] != ""].groupby(keys).size()
        parts.append(summary)
    if pathways is not None:
        full, _ = pathways
        counts = pd.crosstab([full["tissue"], full["contrast_id"]], full["support"])
        counts = counts.reindex(columns=list(SUPPORT_ORDER), fill_value=0)
        counts.columns = [f"pathways_{c}" for c in counts.columns]
        conflicting = full[(full["support"] == "fused") & (full["direction_agreement"] == "conflicting")]
        counts["pathways_fused_conflicting"] = conflicting.groupby(keys).size()
        parts.append(counts)
    if not parts:
        return pd.DataFrame(columns=keys)
    table = pd.concat(parts, axis=1).reset_index()
    counts = [c for c in table.columns if c.startswith(("pathways_", "short_list", "reactions_"))]
    table[counts] = table[counts].fillna(0).astype(int)
    return table


def _readme(integration_dir: Path, info: PanelInfo | None, reactions: tuple | None, pathways: tuple | None,
            params: dict[str, Any]) -> pd.DataFrame:
    manifest = read_manifest(integration_dir)
    rows = [("整合结果包", manifest.run_id), ("单组学结果包", manifest.parameters.get("source_run", ""))]
    if info is None:
        rows.append(("面板", "本次整合没有声明面板；候选表与面板无关，照常列出"))
    else:
        rows += [("面板声明", info.file or ""), ("面板选定方式", SELECTION_LABELS.get(info.selection, info.selection)),
                 ("选定依据", info.basis or "（未写）"), ("面板", "、".join(sorted(info.panels)) or "（未标注）")]
        if info.note:
            rows.append(("面板说明", info.note))
    rows.append(("面板与入选", "panels 列只作标注，不影响入选；同一结果包换面板时，候选行与其他各列不变"))
    if reactions is not None:
        decoys = reactions[0]["decoy_L1_expected"].notna().any()
        rows += [("反应简表", f"通量结论 L1 及以上的全部，加每个组织 × 比较排名前 {REACTION_TOP_N} 的反应"
                                f"（selected_by）；全量表见 reactions_full.tsv"),
                 ("诱饵校准", "cell_L1 为该组织 × 比较的 L1 数，decoy_L1_expected 为打乱成员后的期望；两者接近时该格的 "
                              "L1 与随机水平相当，只作探索" if decoys
                  else "本次整合未运行诱饵校准（reaction_dossier.calibration.decoys 为 0），L1 无随机对照")]
    if pathways is not None:
        rows += [("通路匹配", "跨层只按规范化名称匹配（去库名前缀、小写、合并空白），没有按通路 ID 对齐"),
                 ("通路显著", f"单层：ORA padj < {params['ora_padj_cutoff']}，GSEA FDR < {params['gsea_padj_cutoff']}；"
                              f"融合：fused_fdr < {FUSED_ALPHA}（Fisher 合并后分层 BH）"),
                 ("融合值的用途", "用于探索排序，不是经验证的发现错误率；报告时同时给出各层效应和校正值"),
                 ("通路简表", f"fused 与 multi_layer 的全部，加每个组织 × 比较 × 组学单层显著中 padj 最小的 "
                              f"{PATHWAY_TOP_N} 条；全量表见 pathways_full.tsv")]
    return pd.DataFrame(rows, columns=["项目", "说明"])


def write_candidates(integration_dir: Path, out_dir: Path, project_root: Path) -> Candidates | None:
    """在 ``out_dir`` 写出 candidates.xlsx（说明、反应候选、通路候选）与各表的 TSV；没有可用的全局表时返回 None。"""
    info = panel_info(integration_dir, project_root)
    reactions = reaction_candidates(integration_dir, info.panels if info else None)
    pathways = pathway_candidates(integration_dir)
    if reactions is None and pathways is None:
        return None
    params = {**CONVERGENCE_DEFAULTS, **(read_manifest(integration_dir).parameters.get(CONVERGENCE) or {})}
    out_dir.mkdir(parents=True, exist_ok=True)
    files = [out_dir / WORKBOOK]
    with pd.ExcelWriter(files[0], engine="openpyxl") as writer:
        _readme(integration_dir, info, reactions, pathways, params).to_excel(writer, sheet_name="说明", index=False)
        for name, tables in (("反应候选", reactions), ("通路候选", pathways)):
            if tables is not None:
                tables[1].to_excel(writer, sheet_name=name, index=False)
                writer.sheets[name].freeze_panes = "A2"
                writer.sheets[name].auto_filter.ref = writer.sheets[name].dimensions
        writer.sheets["说明"].column_dimensions["A"].width = 16
        writer.sheets["说明"].column_dimensions["B"].width = 110
    for stem, tables in (("reactions", reactions), ("pathways", pathways)):
        if tables is not None:
            for path, table in ((out_dir / f"{stem}_full.tsv", tables[0]), (out_dir / f"{stem}.tsv", tables[1])):
                write_excel_tsv(table, path)
                files.append(path)
    return Candidates(out_dir, _overview(reactions, pathways), info, files)
