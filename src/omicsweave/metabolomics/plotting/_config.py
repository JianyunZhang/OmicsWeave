"""
绘图模块共享配置：中文字体、样式设置。
"""
from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def setup_plotting():
    """设置绘图样式和中文字体。"""
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
    plt.rcParams["axes.unicode_minus"] = False


setup_plotting()