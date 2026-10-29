"""
绘图模块共享配置：导入时应用统一主题（代谢组图件沿用白色网格样式）。
"""
from __future__ import annotations

from omicsweave.visualization.theme import apply_theme

apply_theme(style="seaborn-v0_8-whitegrid")
