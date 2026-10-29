#!/usr/bin/env python3
"""向后兼容重导出；统计实现位于 :mod:`omicsweave.common`。"""

from omicsweave.common.meta_analysis import MetaAnalysisResult, random_effects_meta

__all__ = ["MetaAnalysisResult", "random_effects_meta"]
