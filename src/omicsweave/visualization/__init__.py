"""跨组学共用的可视化：统一主题（theme）、图件保存（figures）与由统一结果表绘制的通用图。

已有的通用图：火山图（:mod:`omicsweave.visualization.volcano`）。热图、PCA 等仍在各组学的 ``plotting/`` 中，
阶段 2 陆续迁入。
"""

from .theme import apply_theme, font_family, has_cjk_font, theme_rc, themed

__all__ = ["apply_theme", "font_family", "has_cjk_font", "theme_rc", "themed"]
