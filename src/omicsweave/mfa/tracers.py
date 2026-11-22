#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
标准化同位素示踪剂（Isotope Tracer）注册表与图面展示名解析。

涵盖：
- 现有 4 种 ``U-13C`` 标记底物（向后兼容历史表中的裸名 ``glucose`` / ``glycerol`` / ``alanine`` / ``lactate``）：
  ``U-13C-glucose``、``U-13C-glycerol``、``U-13C-alanine``、``U-13C-lactate``；
- 后续扩展 12 种 ``U-13C`` / ``U-15N`` / ``U-13C-15N`` 标记底物：
  ``U-13C-palmitic acid``，以及 ``alanine`` / ``glutamine`` / ``glycine`` / ``leucine`` 的全套
  ``U-13C``、``U-15N`` 与 ``U-13C-15N`` 双标记示踪剂（共 16 种规范示踪剂）。

建模原则：同一代谢物的不同同位素标记（如 ``U-13C-alanine``、``U-15N-alanine``、``U-13C-15N-alanine``）
在 ``samples.tsv`` 的 ``protocol`` 列中必须各自作为独立分层处理，不可将不同同位素标记的动物混为同一组拟合。
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class TracerSpec:
    """单种同位素示踪剂的规范元数据。"""

    canonical_id: str
    display_name: str
    short_label: str
    substrate: str
    isotope_label: str              # "U-13C" | "U-15N" | "U-13C-15N"
    tracked_atoms: tuple[str, ...]  # ("C",) | ("N",) | ("C", "N")
    entry_node: str
    target_regions: tuple[str, ...]
    aliases: tuple[str, ...] = ()
    notes: str = ""


TRACER_SPECS: tuple[TracerSpec, ...] = (
    # ---- 现有 4 种 U-13C 标记底物 ----
    TracerSpec(
        canonical_id="u13c_glucose",
        display_name="U-13C-glucose",
        short_label="U-13C-Glc",
        substrate="glucose",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="GLC_e",
        target_regions=("glycolysis_gng", "ppp", "pyruvate_tca", "glycogen"),
        aliases=(
            "glucose", "glc", "13c_glucose", "13c-glucose", "13c glucose",
            "u13c_glucose", "u13c-glucose", "u-13c-glucose", "u-13c glucose",
            "d-glucose", "u_13c_glucose",
        ),
        notes="主要约束糖酵解、PPP 与 PDH/TCA 入口碳通量",
    ),
    TracerSpec(
        canonical_id="u13c_glycerol",
        display_name="U-13C-glycerol",
        short_label="U-13C-Glyc",
        substrate="glycerol",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="GLYCEROL_e",
        target_regions=("glycerolipid", "glycolysis_gng"),
        aliases=(
            "glycerol", "glyc", "13c_glycerol", "13c-glycerol", "13c glycerol",
            "u13c_glycerol", "u13c-glycerol", "u-13c-glycerol", "u-13c glycerol",
            "u_13c_glycerol",
        ),
        notes="主要约束甘油激酶 (GK)、G3P 穿梭与糖异生入口",
    ),
    TracerSpec(
        canonical_id="u13c_alanine",
        display_name="U-13C-alanine",
        short_label="U-13C-Ala",
        substrate="alanine",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="ALA_e",
        target_regions=("amino_acid", "pyruvate_tca", "glycolysis_gng"),
        aliases=(
            "alanine", "ala", "13c_alanine", "13c-alanine", "13c alanine",
            "u13c_alanine", "u13c-alanine", "u-13c-alanine", "u-13c alanine",
            "u_13c_alanine", "alanine_13c", "alanine_u13c",
        ),
        notes="主要约束丙氨酸摄取、GPT 碳骨架、丙酮酸羧化 (PC) 与糖异生",
    ),
    TracerSpec(
        canonical_id="u13c_lactate",
        display_name="U-13C-lactate",
        short_label="U-13C-Lac",
        substrate="lactate",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="LAC_e",
        target_regions=("glycolysis_gng", "pyruvate_tca"),
        aliases=(
            "lactate", "lac", "13c_lactate", "13c-lactate", "13c lactate",
            "u13c_lactate", "u13c-lactate", "u-13c-lactate", "u-13c lactate",
            "u_13c_lactate",
        ),
        notes="主要约束 LDH、Cori 循环与丙酮酸进入 TCA 氧化/糖异生",
    ),
    # ---- 新增：脂肪酸 U-13C-palmitic acid ----
    TracerSpec(
        canonical_id="u13c_palmitic_acid",
        display_name="U-13C-palmitic acid",
        short_label="U-13C-Pal",
        substrate="palmitic acid",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="PAL_e",
        target_regions=("fatty_acid", "ketone", "peroxisome", "pyruvate_tca"),
        aliases=(
            "palmitic_acid", "palmitic acid", "palmitate", "pal",
            "13c_palmitic_acid", "13c-palmitic-acid", "13c palmitic acid",
            "13c_palmitate", "13c-palmitate", "13c palmitate",
            "u13c_palmitic_acid", "u13c_palmitate", "u-13c-palmitic acid",
            "u-13c-palmitic_acid", "u-13c palmitic acid", "u-13c-palmitate",
            "u-13c palmitate", "u_13c_palmitic_acid", "u_13c_palmitate",
        ),
        notes="主要约束长链脂肪酸摄取、酯化/脂解、CPT1/β-氧化、生酮与 TCA 乙酰辅酶 A 供碳",
    ),
    # ---- 新增：丙氨酸 U-15N 与 U-13C-15N ----
    TracerSpec(
        canonical_id="u15n_alanine",
        display_name="U-15N-alanine",
        short_label="U-15N-Ala",
        substrate="alanine",
        isotope_label="U-15N",
        tracked_atoms=("N",),
        entry_node="ALA_e",
        target_regions=("amino_acid", "urea_creatine"),
        aliases=(
            "15n_alanine", "15n-alanine", "15n alanine",
            "alanine_15n", "alanine-15n", "alanine_u15n",
            "u15n_alanine", "u15n-alanine", "u-15n-alanine", "u-15n alanine",
            "u_15n_alanine",
        ),
        notes="主要约束丙氨酸氨基氮转氨 (GPT/ASPTA) 与尿素循环氮去向",
    ),
    TracerSpec(
        canonical_id="u13c_15n_alanine",
        display_name="U-13C-15N-alanine",
        short_label="U-13C15N-Ala",
        substrate="alanine",
        isotope_label="U-13C-15N",
        tracked_atoms=("C", "N"),
        entry_node="ALA_e",
        target_regions=("amino_acid", "pyruvate_tca", "glycolysis_gng", "urea_creatine"),
        aliases=(
            "alanine_13c15n", "alanine_13c_15n", "alanine_u13c15n", "alanine_u13c_15n",
            "13c15n_alanine", "13c_15n_alanine", "13c-15n-alanine", "13c/15n alanine",
            "13c/15n_alanine", "u13c15n_alanine", "u13c_15n_alanine",
            "u-13c-15n-alanine", "u-13c-15n alanine", "u-13c15n-alanine",
            "u_13c_15n_alanine",
        ),
        notes="碳氮双标记丙氨酸；同步约束碳骨架糖异生/TCA 与氨基氮转氨（独立协议分层）",
    ),
    # ---- 新增：谷氨酰胺 U-13C / U-15N / U-13C-15N ----
    TracerSpec(
        canonical_id="u13c_glutamine",
        display_name="U-13C-glutamine",
        short_label="U-13C-Gln",
        substrate="glutamine",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="GLN_e",
        target_regions=("amino_acid", "pyruvate_tca"),
        aliases=(
            "glutamine", "gln", "13c_glutamine", "13c-glutamine", "13c glutamine",
            "u13c_glutamine", "u13c-glutamine", "u-13c-glutamine", "u-13c glutamine",
            "u_13c_glutamine", "glutamine_13c", "glutamine_u13c",
        ),
        notes="主要约束谷氨酰胺摄取、谷氨酰胺分解 (GLND/GLUD) 与 TCA 回补/还原羧化",
    ),
    TracerSpec(
        canonical_id="u15n_glutamine",
        display_name="U-15N-glutamine",
        short_label="U-15N-Gln",
        substrate="glutamine",
        isotope_label="U-15N",
        tracked_atoms=("N",),
        entry_node="GLN_e",
        target_regions=("amino_acid", "urea_creatine", "nucleotide", "glycolysis_gng"),
        aliases=(
            "15n_glutamine", "15n-glutamine", "15n glutamine",
            "glutamine_15n", "glutamine-15n", "glutamine_u15n",
            "u15n_glutamine", "u15n-glutamine", "u-15n-glutamine", "u-15n glutamine",
            "u_15n_glutamine",
        ),
        notes="主要约束谷氨酰胺酰胺/氨基供氮、尿素循环、嘌呤/嘧啶从头合成与己糖胺途径",
    ),
    TracerSpec(
        canonical_id="u13c_15n_glutamine",
        display_name="U-13C-15N-glutamine",
        short_label="U-13C15N-Gln",
        substrate="glutamine",
        isotope_label="U-13C-15N",
        tracked_atoms=("C", "N"),
        entry_node="GLN_e",
        target_regions=("amino_acid", "pyruvate_tca", "urea_creatine", "nucleotide"),
        aliases=(
            "glutamine_13c15n", "glutamine_13c_15n", "glutamine_u13c15n",
            "13c15n_glutamine", "13c_15n_glutamine", "13c-15n-glutamine", "13c/15n glutamine",
            "u13c15n_glutamine", "u13c_15n_glutamine",
            "u-13c-15n-glutamine", "u-13c-15n glutamine", "u-13c15n-glutamine",
            "u_13c_15n_glutamine",
        ),
        notes="碳氮双标记谷氨酰胺；同步追踪谷氨酰胺碳回补与氮供体通量",
    ),
    # ---- 新增：甘氨酸 U-13C / U-15N / U-13C-15N ----
    TracerSpec(
        canonical_id="u13c_glycine",
        display_name="U-13C-glycine",
        short_label="U-13C-Gly",
        substrate="glycine",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="GLY_c",
        target_regions=("other", "amino_acid", "urea_creatine", "nucleotide"),
        aliases=(
            "glycine", "gly", "13c_glycine", "13c-glycine", "13c glycine",
            "u13c_glycine", "u13c-glycine", "u-13c-glycine", "u-13c glycine",
            "u_13c_glycine", "glycine_13c", "glycine_u13c",
        ),
        notes="主要约束丝氨酸/甘氨酸互转 (SHMT)、一碳叶酸循环、谷胱甘肽、肌酸与嘌呤碳通量",
    ),
    TracerSpec(
        canonical_id="u15n_glycine",
        display_name="U-15N-glycine",
        short_label="U-15N-Gly",
        substrate="glycine",
        isotope_label="U-15N",
        tracked_atoms=("N",),
        entry_node="GLY_c",
        target_regions=("other", "amino_acid", "urea_creatine", "nucleotide"),
        aliases=(
            "15n_glycine", "15n-glycine", "15n glycine",
            "glycine_15n", "glycine-15n", "glycine_u15n",
            "u15n_glycine", "u15n-glycine", "u-15n-glycine", "u-15n glycine",
            "u_15n_glycine",
        ),
        notes="主要约束甘氨酸氨基氮进入丝氨酸、谷胱甘肽、肌酸、血红素与嘌呤环",
    ),
    TracerSpec(
        canonical_id="u13c_15n_glycine",
        display_name="U-13C-15N-glycine",
        short_label="U-13C15N-Gly",
        substrate="glycine",
        isotope_label="U-13C-15N",
        tracked_atoms=("C", "N"),
        entry_node="GLY_c",
        target_regions=("other", "amino_acid", "urea_creatine", "nucleotide"),
        aliases=(
            "glycine_13c15n", "glycine_13c_15n", "glycine_u13c15n",
            "13c15n_glycine", "13c_15n_glycine", "13c-15n-glycine", "13c/15n glycine",
            "u13c15n_glycine", "u13c_15n_glycine",
            "u-13c-15n-glycine", "u-13c-15n glycine", "u-13c15n-glycine",
            "u_13c_15n_glycine",
        ),
        notes="碳氮双标记甘氨酸；同步追踪甘氨酸裂解、SHMT 与一碳/GSH/肌酸合成",
    ),
    # ---- 新增：亮氨酸 U-13C / U-15N / U-13C-15N ----
    TracerSpec(
        canonical_id="u13c_leucine",
        display_name="U-13C-leucine",
        short_label="U-13C-Leu",
        substrate="leucine",
        isotope_label="U-13C",
        tracked_atoms=("C",),
        entry_node="LEU_c",
        target_regions=("amino_acid", "ketone", "pyruvate_tca"),
        aliases=(
            "leucine", "leu", "13c_leucine", "13c-leucine", "13c leucine",
            "u13c_leucine", "u13c-leucine", "u-13c-leucine", "u-13c leucine",
            "u_13c_leucine", "leucine_13c", "leucine_u13c",
        ),
        notes="主要约束亮氨酸/BCAA 碳骨架分解、生酮乙酰乙酰辅酶 A 与 TCA 供碳",
    ),
    TracerSpec(
        canonical_id="u15n_leucine",
        display_name="U-15N-leucine",
        short_label="U-15N-Leu",
        substrate="leucine",
        isotope_label="U-15N",
        tracked_atoms=("N",),
        entry_node="LEU_c",
        target_regions=("amino_acid", "urea_creatine"),
        aliases=(
            "15n_leucine", "15n-leucine", "15n leucine",
            "leucine_15n", "leucine-15n", "leucine_u15n",
            "u15n_leucine", "u15n-leucine", "u-15n-leucine", "u-15n leucine",
            "u_15n_leucine",
        ),
        notes="主要约束支链氨基酸转氨酶 (BCAT) 向谷氨酸/丙氨酸/谷氨酰胺的氮传递",
    ),
    TracerSpec(
        canonical_id="u13c_15n_leucine",
        display_name="U-13C-15N-leucine",
        short_label="U-13C15N-Leu",
        substrate="leucine",
        isotope_label="U-13C-15N",
        tracked_atoms=("C", "N"),
        entry_node="LEU_c",
        target_regions=("amino_acid", "ketone", "pyruvate_tca", "urea_creatine"),
        aliases=(
            "leucine_13c15n", "leucine_13c_15n", "leucine_u13c15n",
            "13c15n_leucine", "13c_15n_leucine", "13c-15n-leucine", "13c/15n leucine",
            "u13c15n_leucine", "u13c_15n_leucine",
            "u-13c-15n-leucine", "u-13c-15n leucine", "u-13c15n-leucine",
            "u_13c_15n_leucine",
        ),
        notes="碳氮双标记亮氨酸；同步解析 BCAA 转氨氮流与酮体/TCA 碳流",
    ),
)


def _norm_key(raw: str) -> str:
    text = str(raw or "").strip().lower()
    text = re.sub(r"[\s/]+", "_", text)
    text = text.replace("-", "_")
    text = re.sub(r"_+", "_", text).strip("_")
    return text


TRACER_REGISTRY: dict[str, TracerSpec] = {spec.canonical_id: spec for spec in TRACER_SPECS}
_ALIAS_INDEX: dict[str, TracerSpec] = {}
for _spec in TRACER_SPECS:
    for _tok in (_spec.canonical_id, _spec.display_name, _spec.short_label, *_spec.aliases):
        _ALIAS_INDEX[_norm_key(_tok)] = _spec


def list_tracers() -> tuple[TracerSpec, ...]:
    """返回全部已注册规范同位素示踪剂（按规范顺序，共 16 种）。"""
    return TRACER_SPECS


def resolve_tracer(raw: str | None) -> TracerSpec | None:
    """将任意协议名/别名（如 ``glucose``、``alanine_13c15n``、``U-13C-palmitic acid``）解析为 :class:`TracerSpec`。"""
    if raw is None:
        return None
    key = _norm_key(str(raw))
    if not key or key == "pooled":
        return None
    return _ALIAS_INDEX.get(key)


def is_known_tracer(raw: str | None) -> bool:
    return resolve_tracer(raw) is not None


def canonical_tracer_id(raw: str | None) -> str:
    """返回规范机读 ID（如 ``u13c_glucose``）；``pooled`` 原样返回；未注册名保留去空白原值。"""
    s = str(raw or "").strip()
    if not s or s.lower() == "pooled":
        return s.lower()
    spec = resolve_tracer(s)
    return spec.canonical_id if spec is not None else s


def format_tracer_display(raw: str | None) -> str:
    """返回面向图例、标题与报告的规范同位素名（如 ``U-13C-glucose``、``U-13C-palmitic acid``、``Pooled (IVW)``）。"""
    s = str(raw or "").strip()
    if not s:
        return ""
    if s.lower() == "pooled":
        return "Pooled (IVW)"
    spec = resolve_tracer(s)
    return spec.display_name if spec is not None else s


def format_tracer_short(raw: str | None) -> str:
    """返回面向反应卡片微徽章的紧凑同位素短标（如 ``U-13C-Glc``、``U-13C15N-Ala``、``Pooled``）。"""
    s = str(raw or "").strip()
    if not s:
        return ""
    if s.lower() == "pooled":
        return "Pooled"
    spec = resolve_tracer(s)
    return spec.short_label if spec is not None else s


def format_tracer_badge_label(tracers: Sequence[str] | str | None, *, level: str) -> str:
    """根据实际贡献该反应的示踪剂列表与证据等级，生成自解释反应徽章文字：

    - ``L1`` → 恒为 ``Enz-pred``；
    - ``L2``（单/多示踪剂 MFA 支撑）→ ``U-13C-Glc Flux`` 或 ``U-13C-Glc+Ala Flux``（无示踪剂信息时回退 ``Flux``）；
    - ``L3``（多示踪剂一致确认）→ 直接展示多重示踪剂短名（如 ``U-13C-Glc+Ala``、``U-13C-Ala+U-15N-Ala``；
      无示踪剂信息时回退 ``Confirmed``）。
    """
    if level == "L1":
        return "Enz-pred"
    if isinstance(tracers, str):
        items = [t.strip() for t in re.split(r"[;,]", tracers) if t.strip() and t.strip().lower() != "pooled"]
    elif tracers:
        items = [str(t).strip() for t in tracers if str(t).strip() and str(t).strip().lower() != "pooled"]
    else:
        items = []

    # 按规范名去重并保持稳定顺序
    seen: set[str] = set()
    resolved: list[tuple[TracerSpec | None, str]] = []
    for item in items:
        spec = resolve_tracer(item)
        key = spec.canonical_id if spec is not None else item.lower()
        if key not in seen:
            seen.add(key)
            resolved.append((spec, item))

    if not resolved:
        return "Confirmed" if level == "L3" else ("Flux" if level == "L2" else "")

    if len(resolved) == 1:
        spec, raw = resolved[0]
        short = spec.short_label if spec is not None else raw
        return f"{short} Flux" if level == "L2" else short

    # 多重示踪剂：若同属单一标记类型（如全为 U-13C），合并前缀以保持卡片紧凑且无歧义；
    # 若混用不同标记类型（如 U-13C + U-15N / U-13C-15N），保留各自完整同位素前缀绝不混淆。
    specs_only = [sp for sp, _ in resolved]
    if all(sp is not None for sp in specs_only):
        iso_set = {sp.isotope_label for sp in specs_only if sp is not None}
        if len(iso_set) == 1:
            iso = next(iter(iso_set))
            iso_prefix = iso.replace("U-13C-15N", "U-13C15N") + "-"
            suffixes = [
                sp.short_label[len(iso_prefix):] if sp.short_label.startswith(iso_prefix) else sp.short_label
                for sp in specs_only if sp is not None
            ]
            if len(suffixes) <= 3:
                multi = f"{iso_prefix}" + "+".join(suffixes)
            else:
                multi = f"{iso_prefix}{suffixes[0]}+{suffixes[1]}+{len(suffixes) - 2}"
        else:
            shorts = [sp.short_label for sp in specs_only if sp is not None]
            multi = "+".join(shorts[:2]) + (f"+{len(shorts) - 2}" if len(shorts) > 2 else "")
    else:
        shorts = [sp.short_label if sp is not None else raw for sp, raw in resolved]
        multi = "+".join(shorts[:2]) + (f"+{len(shorts) - 2}" if len(shorts) > 2 else "")

    return f"{multi} Flux" if level == "L2" else multi


__all__ = [
    "TRACER_REGISTRY",
    "TRACER_SPECS",
    "TracerSpec",
    "canonical_tracer_id",
    "format_tracer_badge_label",
    "format_tracer_display",
    "format_tracer_short",
    "is_known_tracer",
    "list_tracers",
    "resolve_tracer",
]
