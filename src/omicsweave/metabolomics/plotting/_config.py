"""
绘图模块共享配置：导入时应用统一 ``metro_atlas`` 出版主题（300 dpi）。
"""
from __future__ import annotations

from omicsweave.visualization.theme import apply_theme

apply_theme(preset="metro_atlas", save_dpi=300)
