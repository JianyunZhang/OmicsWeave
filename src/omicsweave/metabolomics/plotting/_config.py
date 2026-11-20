"""
绘图模块共享配置：导入时应用统一 ``nature_editorial`` 出版主题（300 dpi）。
"""
from __future__ import annotations

from omicsweave.visualization.theme import apply_theme

apply_theme(preset="nature_editorial", save_dpi=300)
