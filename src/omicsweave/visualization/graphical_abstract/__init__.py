"""独立的 Graphical Abstract（概念摘要图）子系统（严格限定 ≤ 12 个核心枢纽）。

提供四套专用渲染皮肤：
- ``nature_editorial``：顶刊极简机制图；
- ``editorial_ink``：手绘科学插画风（有机细胞器墨线淡彩 + 精密矢量数据层）；
- ``pixel_abstract``：现代 16-bit 像素概念风（汇报总结 / 封面构思 / 项目首页 / 科普展示）；
- ``metro_atlas``：Beck 代谢地铁换乘图风（八向轨道 + 内嵌方向 chevron + 跑道闭环 + 圆环换乘站，支持 ``semantic`` / ``pathway`` 双配色模式）。
"""
from .renderer import GA_STYLES, VALID_COLOR_MODES, render_graphical_abstract
from .spec import (
    MAX_HUBS,
    VALID_COMPARTMENTS,
    VALID_ICONS,
    AbstractCompartment,
    AbstractFlow,
    AbstractHub,
    AbstractSpec,
    build_abstract_from_pathway_view,
    load_abstract_spec,
    spec_from_dict,
    validate_spec,
)

__all__ = [
    "AbstractCompartment",
    "AbstractFlow",
    "AbstractHub",
    "AbstractSpec",
    "GA_STYLES",
    "MAX_HUBS",
    "VALID_COLOR_MODES",
    "VALID_COMPARTMENTS",
    "VALID_ICONS",
    "build_abstract_from_pathway_view",
    "load_abstract_spec",
    "render_graphical_abstract",
    "spec_from_dict",
    "validate_spec",
]
