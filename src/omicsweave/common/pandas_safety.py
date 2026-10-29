"""
pandas 大矩阵广播运算的安全闸门
================================

pandas 在元素数超过 ``expressions._MIN_ELEMENTS``（默认 100 万）时，会把
``DataFrame.<op>(Series, axis=...)`` 交给 numexpr 加速。走这条路时 pandas
先用 ``np.broadcast_to`` 把 Series 摊成 **0-stride 的二维视图**再送进
numexpr。2026-09-01 在本机 numexpr 2.10.1 + pandas 2.3.3 + numpy 2.3.5 上
实测：该组合会让约 **93% 的输出元素变成 0**，而且**每次运行的错误位置都
不一样**（非确定性）。

事故是静默的。总和会被下游的重标定拉回正确量级，于是文库大小自检、总
计数自检全部通过 —— 只有「非零元素个数」会露馅。本项目因此产出过一整套
错误的转录组差异表达结果：35,983 基因 × 95 样本 = 341 万元素，刚好越过
阈值；**同一份代码在 19 样本 / 68 万元素下结果完全正确**，所以小规模的
单元测试天然测不到它。

本模块提供两件工具：

* :func:`harden_pandas` —— 关闭 pandas 的 numexpr 路径。代价是大矩阵逐
  元素运算回落到 numpy；本项目耗时由 kallisto 与 DESeq2 主导，可忽略。
* :func:`verify_broadcast_ops` / :func:`assert_broadcast_ops_sane` ——
  用一个跨过阈值的探针实测真值，把「环境又变了」变成一次明确的报错，
  而不是一批看不出问题的结果。

.. warning::
   ``pd.get_option("compute.use_numexpr")`` 为 ``True`` **不代表**正在使用
   numexpr。真正的开关是 ``expressions.USE_NUMEXPR``（= 选项 **且** 已安装）。
   卸载 numexpr 后选项仍是 True，只看选项会得出完全相反的结论。

.. note::
   这道闸门只防「环境层面的静默算错」。真正的兜底是对**输出**做不变量
   校验 —— 见 ``omicsweave/transcriptomics/normalization.py`` 里的
   ``_assert_aggregation_sane()``，它不依赖对成因的判断。
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pandas.core.computation import expressions as _expr

__all__ = [
    "numexpr_status",
    "harden_pandas",
    "verify_broadcast_ops",
    "assert_broadcast_ops_sane",
]

# 探针形状：元素数必须**大于** pandas 的 numexpr 阈值，否则根本测不到
# 这条代码路径。1,400,021 个 float64 约 11 MB，一次探测在毫秒量级。
_PROBE_SHAPE: tuple[int, int] = (200_003, 7)

_MIN_ELEMENTS_DEFAULT = 1_000_000


def _min_elements() -> int:
    return int(getattr(_expr, "_MIN_ELEMENTS", _MIN_ELEMENTS_DEFAULT))


def numexpr_status() -> dict[str, object]:
    """
    返回 pandas 当前**是否真的**在走 numexpr 路径。

    ``active`` 才是有效开关；``option_enabled`` 单独看没有意义。
    """
    try:
        import numexpr  # noqa: PLC0415 — 探测用，缺失是正常情况
        version: str | None = str(numexpr.__version__)
    except ImportError:
        version = None

    return {
        "numexpr_version": version,
        "installed": bool(getattr(_expr, "NUMEXPR_INSTALLED", False)),
        "option_enabled": bool(pd.get_option("compute.use_numexpr")),
        "active": bool(getattr(_expr, "USE_NUMEXPR", False)),
        "min_elements": _min_elements(),
    }


def harden_pandas() -> dict[str, object]:
    """
    关闭 pandas 的 numexpr 加速路径，返回**关闭前**的状态。

    幂等；numexpr 未安装时是空操作。放在模块导入期调用是刻意的：
    该路径一旦被用到就已经算错了，没有「先跑一段再关」的余地。
    """
    before = numexpr_status()
    if before["active"]:
        pd.set_option("compute.use_numexpr", False)
    return before


def verify_broadcast_ops(shape: tuple[int, int] = _PROBE_SHAPE) -> list[str]:
    """
    实测 ``DataFrame.<op>(Series, axis=0/1)`` 跨过 numexpr 阈值时是否算对。

    探针值全部 ``>= 1``，因此「输出被填成 0」这种失效必定被抓到。

    Parameters
    ----------
    shape
        探针矩阵形状，元素数必须大于 pandas 的 numexpr 阈值。

    Returns
    -------
    list[str]
        出错的运算名；全部正确时返回空列表。
    """
    n_rows, n_cols = shape
    if n_rows * n_cols <= _min_elements():
        raise ValueError(
            f"探针元素数 {n_rows * n_cols} 未超过 pandas 阈值 {_min_elements()}，"
            "测不到 numexpr 路径"
        )

    rng = np.random.default_rng(0)
    values = rng.random(shape) + 1.0
    by_row = rng.random(n_rows) + 1.0
    by_col = rng.random(n_cols) + 1.0
    frame = pd.DataFrame(values)
    row_series = pd.Series(by_row, index=frame.index)
    col_series = pd.Series(by_col, index=frame.columns)

    cases: dict[str, tuple[pd.DataFrame, np.ndarray]] = {
        "mul(axis=0)": (frame.mul(row_series, axis=0), values * by_row[:, None]),
        "div(axis=0)": (frame.div(row_series, axis=0), values / by_row[:, None]),
        "sub(axis=0)": (frame.sub(row_series, axis=0), values - by_row[:, None]),
        "mul(axis=1)": (frame.mul(col_series, axis=1), values * by_col[None, :]),
        "div(axis=1)": (frame.div(col_series, axis=1), values / by_col[None, :]),
    }

    return [
        name
        for name, (got, want) in cases.items()
        if not np.allclose(got.to_numpy(), want, rtol=1e-12, atol=0.0)
    ]


def assert_broadcast_ops_sane(shape: tuple[int, int] = _PROBE_SHAPE) -> None:
    """探针不通过就抛 ``RuntimeError``，附带足以定位的环境信息。"""
    failures = verify_broadcast_ops(shape)
    if failures:
        status = numexpr_status()
        raise RuntimeError(
            "pandas 大矩阵广播运算返回了错误结果，拒绝继续 —— 任何在此环境下"
            "产出的组学矩阵都不可信。\n"
            f"  失效运算 : {', '.join(failures)}\n"
            f"  探针形状 : {shape[0]} x {shape[1]} = {shape[0] * shape[1]} 元素"
            f"（阈值 {status['min_elements']}）\n"
            f"  numexpr  : version={status['numexpr_version']} "
            f"installed={status['installed']} active={status['active']}\n"
            "  处理方式 : 卸载 numexpr（conda remove numexpr），或确认"
            "harden_pandas() 已在流程启动时被调用。"
        )
