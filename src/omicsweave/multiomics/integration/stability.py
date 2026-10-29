#!/usr/bin/env python3
"""兼容导入路径；实现位于 :mod:`omicsweave.common.stability`。"""

from omicsweave.common.stability import (
    StabilityResult,
    calibrate_biological_stability,
)

__all__ = ["StabilityResult", "calibrate_biological_stability"]
