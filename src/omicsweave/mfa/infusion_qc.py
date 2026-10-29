#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
灌注期血糖稳态 QC（MFA 的前提假设检验）

13C-MFA 的整套推断都建立在一个假设上：**代谢与同位素在采样时已达稳态**。
只有稳态成立，"一组固定的通量"才是对该时刻代谢状态的有效描述；否则拟合
出的通量是一段漂移过程的平均，既不对应任何真实生理状态，也不该被拿来做
组间比较。

对静脉输注 U-13C 葡萄糖的在体实验，血糖时间曲线是这个假设最直接、也最
容易获得的检验：若灌注期内血糖持续上升或剧烈波动，代谢稳态就没有建立。

本模块读取批次的血糖记录表，在配置的稳态窗口内计算每只动物的
变异系数、极差与线性漂移速率，按 ``protocol.yaml: infusion_qc_defaults``
声明的阈值判定通过与否。

这些阈值和记录文件早在配置里声明，但此前没有任何代码读取——稳态假设
从未被真正检验过。

输出
----
``01_df_results/all_tissues/qc/infusion_steady_state_qc.csv``
    每只动物一行：窗口内统计量、判定结果与失败原因。
``01_df_results/all_tissues/qc/infusion_glycemia_timecourse.csv``
    长表形式的血糖时间曲线，供绘图与人工核查。
``01_df_results/all_tissues/qc/infusion_steady_state_unmapped.csv``
    在血糖表中、但不在本批次动物分组里的动物（见下文「动物分组」），
    列与主表相同；没有这类动物时只写表头，避免旧结果残留。

动物分组
--------
血糖表按动物编号（``mice_id``）记录，条件因此要从「动物编号 → 条件」的
映射里取。``groups.yaml`` 里两个字段的含义因实验而异：

- ``20260312_bjsjth_wt_mice_dio_sleeve_gastrectomy``：``sample_groups`` 以
  动物编号（``A1``、``A12``）为键，可直接用；
- ``20250613_cimr_wt_mice_dio_glp``：``sample_groups`` 以 MFA 样本名
  （``hfd_saline_q3d_1w_1``）为键，动物编号只在 ``mouse_groups`` 里。

所以管线经 :func:`select_animal_groups` 优先取 ``mouse_groups``，缺失时才退回
``sample_groups``。

**不在映射中的动物不计入稳态判定，单独报告。** 血糖表可能混入不属于本
批次的动物——例如 GLP 葡萄糖批次的血糖表里有 85 号，它输注的是 U-13C
谷氨酰胺（「灌流量和流速」列注明），不是本批次的样本。把它算进「通过
x/N」或按条件汇总里一个虚构的 ``unmapped`` 组，会歪曲本批次的稳态结论；
直接丢弃又会掩盖分组表里的笔误。因此这类动物照常计算统计量，但只写入
``infusion_steady_state_unmapped.csv`` 并打印 ``[WARN]``。

例外：没有提供映射，或映射与血糖表**完全没有交集**时不拆分，所有动物
留在主表、条件记为 ``unmapped``。整张表都不属于本批次几乎不可能，更可能是
传错了映射（例如把以样本名为键的 ``sample_groups`` 当成动物分组），这时
打印醒目的 ``[WARN]``，而不是让主表变空。

输入格式
--------
血糖时间点列接受两种表头写法：

- 双语长列名 ``120 min (blood_glucose_after_infusion_120_min) (mmol/L)``；
- 纯时间表头 ``0min`` / ``15min`` / ``120 min``（如
  ``20250613_cimr_wt_mice_dio_glp`` 的血糖表）。只匹配整个表头恰好是
  「数字 + min」的列，``灌注2h血酮`` 之类的列不会被误认。

此前只认第一种，遇到第二种时整个 QC 以 ``[SKIP]`` 静默跳过。

单元格里的录入笔误按以下规则处理，且每一处都打印 ``[WARN]``：

- **重复小数点**（``21..5``）：合并成一个小数点（→ ``21.5``）后计入。
  这是唯一合理的读法，丢弃反而会让该动物少一个点、窗口统计量偏移；
  规范化记录写入时间曲线表的 ``value_note`` 与 QC 表的 ``value_notes``。
- **其他无法解析的值**：按缺失处理（不计入统计），同样记录在案，
  不再被 ``pd.to_numeric(errors="coerce")`` 静默吞掉。

样本对应关系的局限
------------------
血糖表用动物编号（``A12``），而 MFA 拟合结果用「组名+序号」
（``glucose__male__hfd_sham_8w__1``），两者无法逐一对应。因此本模块按
**动物**和**条件**两个层级报告，不把 QC 结论下推到单个 MFA 样本。要做到
样本级门控，需要上游在 ``complete_result_name`` 里写入动物编号。
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# 血糖时间点列名形如 "120 min (blood_glucose_after_infusion_120_min) (mmol/L)"
_TIMEPOINT_PATTERN = re.compile(
    r"blood_glucose_after_infusion_(\d+)_min", re.IGNORECASE
)
# 纯时间表头 "0min" / "120 min"。表头不含括号时 _english_column_name 原样
# 返回（仅 strip + 小写），所以匹配规范化后的列名即等同于匹配原始表头。
_BARE_TIMEPOINT_PATTERN = re.compile(r"^\s*(\d+)\s*min\s*$", re.IGNORECASE)
# 双语列名统一用括号里的英文标识符
_ENGLISH_NAME_PATTERN = re.compile(r"\(([a-z0-9_/]+)\)", re.IGNORECASE)
# 录入笔误：重复小数点，如 "21..5"
_REPEATED_DECIMAL_PATTERN = re.compile(r"^\s*([+-]?\d+)\.{2,}(\d+)\s*$")

_DEFAULT_WINDOW = (30, 180)
_DEFAULT_CV_MAX = 0.15
_DEFAULT_RANGE_MAX = 2.0


def _english_column_name(column: str) -> str:
    """从双语列名中取出括号里的英文标识符。"""
    match = _ENGLISH_NAME_PATTERN.search(str(column))
    return match.group(1).lower() if match else str(column).strip().lower()


def _normalize_mouse_id(value: Any) -> str:
    """动物编号转成字符串。

    纯数字编号列只要混进一个空单元格，Excel 读出来就是 float（``86.0``），
    与 ``groups.yaml`` 里的 ``"86"`` 永远对不上，所以整数值去掉 ``.0``。
    """
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def load_glycemia_record(path: Path) -> pd.DataFrame:
    """读取血糖记录表并把双语列名规范成英文标识符。"""
    raw = pd.read_excel(path)
    raw = raw.dropna(how="all")
    raw.columns = [_english_column_name(c) for c in raw.columns]
    if "mice_id" in raw.columns:
        raw = raw[raw["mice_id"].notna()].copy()
        raw["mice_id"] = raw["mice_id"].map(_normalize_mouse_id)
    return raw


def select_animal_groups(
    mouse_groups: dict[str, str] | None,
    sample_groups: dict[str, str] | None,
) -> dict[str, str]:
    """选出「动物编号 → 条件」映射：有 ``mouse_groups`` 用它，否则退回 ``sample_groups``。

    ``sample_groups`` 在部分实验里以 MFA 样本名为键（见模块文档「动物分组」），
    只有 ``mouse_groups`` 缺失时才拿它兜底——20260312 的批次就是这种情况。
    """
    chosen = mouse_groups or sample_groups or {}
    return {str(k).strip(): str(v) for k, v in chosen.items()}


def _split_unmapped(
    qc_df: pd.DataFrame, groups: dict[str, str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """把不在动物分组里的动物从稳态判定中分出来。

    返回 ``(纳入判定的动物, 另行报告的动物)``。没有映射、或映射与血糖表
    完全没有交集时不拆分（理由见模块文档「动物分组」）。
    """
    empty = qc_df.iloc[0:0]
    if not groups:
        print("  [WARN] 未提供动物分组（groups.yaml: mouse_groups / sample_groups），"
              "所有动物的条件记为 unmapped")
        return qc_df, empty

    in_mapping = qc_df["mouse_id"].isin(list(groups))
    if not in_mapping.any():
        print("  [WARN] 动物分组与血糖表的编号完全没有交集，所有动物的条件记为 unmapped。")
        print(f"         分组键示例: {', '.join(list(groups)[:3])}；"
              f"血糖表编号示例: {', '.join(qc_df['mouse_id'].head(3))}")
        print("         若分组键是 MFA 样本名，请在 groups.yaml 中补充 mouse_groups"
              "（动物编号 → 条件）。")
        return qc_df, empty

    unmapped = qc_df[~in_mapping].reset_index(drop=True)
    if len(unmapped):
        print(f"  [WARN] {len(unmapped)} 只动物在血糖表中但不在动物分组中，不计入稳态判定，"
              f"另存 infusion_steady_state_unmapped.csv: {', '.join(unmapped['mouse_id'])}")
    missing = [m for m in groups if m not in set(qc_df["mouse_id"])]
    if missing:
        print(f"  [INFO] 动物分组中有 {len(missing)} 只不在血糖表中: {', '.join(missing)}")
    return qc_df[in_mapping].reset_index(drop=True), unmapped


def _parse_timepoint(column: str) -> int | None:
    """从列名解析灌注后分钟数；不是血糖时间点列时返回 None。"""
    text = str(column)
    for pattern in (_TIMEPOINT_PATTERN, _BARE_TIMEPOINT_PATTERN):
        match = pattern.search(text)
        if match:
            return int(match.group(1))
    return None


def _timepoint_columns(df: pd.DataFrame) -> list[tuple[int, str]]:
    """返回 [(分钟, 列名)]，按时间升序。

    同一分钟出现多列时只保留第一列——两种表头写法并存时，重复计入同一
    时间点会让 CV 与漂移斜率都被错误加权。
    """
    found: dict[int, str] = {}
    for column in df.columns:
        minute = _parse_timepoint(column)
        if minute is None:
            continue
        if minute in found:
            print(f"  [WARN] {minute} min 对应多个血糖列（{found[minute]!r}、"
                  f"{column!r}），只使用第一列")
            continue
        found[minute] = column
    return sorted(found.items())


def _coerce_glucose(value: Any) -> tuple[float, str]:
    """把一个血糖单元格转成 float，返回 ``(数值, 备注)``。

    备注为空表示原值即合法数字或空缺；``normalized_from:<原值>`` 表示按
    重复小数点笔误规范化；``unparseable:<原值>`` 表示无法解析、按缺失处理。
    """
    if pd.isna(value):
        return np.nan, ""
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.notna(numeric):
        return float(numeric), ""
    text = str(value).strip()
    if not text:
        return np.nan, ""
    match = _REPEATED_DECIMAL_PATTERN.match(text)
    if match:
        return float(f"{match.group(1)}.{match.group(2)}"), f"normalized_from:{text}"
    return np.nan, f"unparseable:{text}"


def _drift_slope(minutes: np.ndarray, values: np.ndarray) -> float:
    """血糖漂移速率（mmol/L 每小时）。少于 3 点时不估计。"""
    if len(values) < 3:
        return float("nan")
    slope_per_min = np.polyfit(minutes, values, 1)[0]
    return float(slope_per_min * 60.0)


def infusion_steady_state_qc(
    glycemia_file: Path | str | None,
    output_dir: Path | str,
    infusion_qc_config: dict[str, Any] | None = None,
    sample_groups: dict[str, str] | None = None,
) -> pd.DataFrame | None:
    """检验灌注期血糖是否达到稳态。

    Parameters
    ----------
    glycemia_file : Path
        批次的血糖记录 xlsx。
    output_dir : Path
        批次输出目录（``runs/<protocol>/<batch>``）。
    infusion_qc_config : dict
        ``protocol.yaml: infusion_qc_defaults`` 合并 batch 覆盖后的结果。
    sample_groups : dict
        **动物编号** → 条件代码，键须与血糖表的 ``mice_id`` 一致。参数名沿用
        历史写法；管线通过 :func:`select_animal_groups` 传入
        ``groups.yaml: mouse_groups``，缺失时才退回 ``sample_groups``——后者在
        部分实验里以 MFA 样本名为键，与血糖表对不上。不在映射中的动物不计入
        稳态判定，另存 ``infusion_steady_state_unmapped.csv``。

    Returns
    -------
    pd.DataFrame or None
        纳入判定的动物每只一行的 QC 结果；未启用或无数据时返回 None。
    """
    config = dict(infusion_qc_config or {})
    if not config:
        print("\n--- Infusion Steady-State QC ---")
        print("  [SKIP] 未配置 infusion_qc_defaults")
        return None
    if not config.get("enabled", True):
        print("\n--- Infusion Steady-State QC ---")
        print("  [SKIP] infusion_qc.enabled = false")
        return None

    print("\n--- Infusion Steady-State QC ---")

    if not glycemia_file:
        print("  [SKIP] 批次未声明 glycemia_record_file")
        return None
    glycemia_path = Path(glycemia_file)
    if not glycemia_path.exists():
        print(f"  [SKIP] 血糖记录文件不存在: {glycemia_path.name}")
        return None

    window = config.get("steady_state_window_min", list(_DEFAULT_WINDOW))
    window_start, window_end = int(window[0]), int(window[1])
    cv_max = float(config.get("glucose_cv_max", _DEFAULT_CV_MAX))
    range_max = float(config.get("glucose_range_max_mmol_l", _DEFAULT_RANGE_MAX))
    drift_max = config.get("glucose_drift_max_mmol_l_per_h", None)
    drift_max = float(drift_max) if drift_max is not None else None
    on_fail = str(config.get("on_fail", "flag")).lower()

    record = load_glycemia_record(glycemia_path)
    timepoints = _timepoint_columns(record)
    if not timepoints:
        print(f"  [SKIP] 未在 {glycemia_path.name} 中识别到血糖时间点列")
        return None

    in_window = [(m, c) for m, c in timepoints if window_start <= m <= window_end]
    if not in_window:
        print(f"  [SKIP] 稳态窗口 [{window_start}, {window_end}] 内没有时间点")
        return None

    groups = {str(k).strip(): str(v) for k, v in (sample_groups or {}).items()}

    rows = []
    timecourse_rows = []
    for _, mouse in record.iterrows():
        mouse_id = str(mouse.get("mice_id", "")).strip()
        if not mouse_id or mouse_id.lower() == "nan":
            continue
        condition = groups.get(mouse_id, "unmapped")

        parsed: dict[int, float] = {}
        value_notes = []
        for minute, column in timepoints:
            raw_value = mouse.get(column)
            value, note = _coerce_glucose(raw_value)
            parsed[minute] = value
            if note.startswith("normalized_from:"):
                print(f"  [WARN] 小鼠 {mouse_id} @ {minute} min: 血糖值 {raw_value!r} "
                      f"按重复小数点笔误规范化为 {value:g}")
            elif note:
                print(f"  [WARN] 小鼠 {mouse_id} @ {minute} min: 血糖值 {raw_value!r} "
                      f"无法解析，按缺失处理")
            if note:
                value_notes.append(f"{minute}min:{note}")
            timecourse_rows.append({
                "mouse_id": mouse_id,
                "condition": condition,
                "minute": minute,
                "glucose_mmol_l": value,
                "in_steady_state_window": window_start <= minute <= window_end,
                "value_note": note,
            })
        value_notes_str = ";".join(value_notes)

        values = np.array([parsed[m] for m, _ in in_window], dtype=float)
        minutes = np.array([m for m, _ in in_window], dtype=float)
        valid = ~np.isnan(values)
        values_valid = values[valid]
        minutes_valid = minutes[valid]

        n_points = int(len(values_valid))
        if n_points == 0:
            rows.append({
                "mouse_id": mouse_id,
                "condition": condition,
                "n_timepoints": 0,
                "steady_state_window_min": f"[{window_start}, {window_end}]",
                "glucose_mean": np.nan,
                "glucose_sd": np.nan,
                "glucose_cv": np.nan,
                "glucose_range_mmol_l": np.nan,
                "glucose_min": np.nan,
                "glucose_max": np.nan,
                "drift_mmol_l_per_h": np.nan,
                "passed": False,
                "fail_reasons": "no_measurements_in_window",
                "value_notes": value_notes_str,
            })
            continue

        mean = float(np.mean(values_valid))
        sd = float(np.std(values_valid, ddof=1)) if n_points >= 2 else np.nan
        cv = float(sd / mean) if (n_points >= 2 and mean != 0) else np.nan
        value_range = float(np.max(values_valid) - np.min(values_valid))
        drift = _drift_slope(minutes_valid, values_valid)

        reasons = []
        if n_points < 3:
            reasons.append(f"only_{n_points}_timepoints")
        if pd.notna(cv) and cv > cv_max:
            reasons.append(f"cv_{cv:.3f}>{cv_max:g}")
        if value_range > range_max:
            reasons.append(f"range_{value_range:.2f}>{range_max:g}mmol/L")
        if drift_max is not None and pd.notna(drift) and abs(drift) > drift_max:
            reasons.append(f"drift_{drift:+.2f}>{drift_max:g}mmol/L/h")

        rows.append({
            "mouse_id": mouse_id,
            "condition": condition,
            "n_timepoints": n_points,
            "steady_state_window_min": f"[{window_start}, {window_end}]",
            "glucose_mean": round(mean, 3),
            "glucose_sd": round(sd, 3) if pd.notna(sd) else np.nan,
            "glucose_cv": round(cv, 4) if pd.notna(cv) else np.nan,
            "glucose_range_mmol_l": round(value_range, 3),
            "glucose_min": float(np.min(values_valid)),
            "glucose_max": float(np.max(values_valid)),
            "drift_mmol_l_per_h": round(drift, 3) if pd.notna(drift) else np.nan,
            "passed": not reasons,
            "fail_reasons": ";".join(reasons),
            "value_notes": value_notes_str,
        })

    if not rows:
        print("  [SKIP] 血糖记录中没有可用的动物行")
        return None

    qc_df = pd.DataFrame(rows)
    qc_df["cv_threshold"] = cv_max
    qc_df["range_threshold_mmol_l"] = range_max
    qc_df["on_fail_policy"] = on_fail

    # 不属于本批次的动物不进入主表、时间曲线与按条件汇总
    qc_df, unmapped_df = _split_unmapped(qc_df, groups)
    timecourse_df = pd.DataFrame(timecourse_rows)
    timecourse_df = timecourse_df[timecourse_df["mouse_id"].isin(qc_df["mouse_id"])]

    save_dir = Path(output_dir) / "01_df_results" / "all_tissues" / "qc"
    save_dir.mkdir(parents=True, exist_ok=True)
    qc_df.to_csv(save_dir / "infusion_steady_state_qc.csv", index=False)
    timecourse_df.to_csv(save_dir / "infusion_glycemia_timecourse.csv", index=False)
    unmapped_df.to_csv(save_dir / "infusion_steady_state_unmapped.csv", index=False)

    n_total = len(qc_df)
    n_pass = int(qc_df["passed"].sum())
    print(f"  窗口 [{window_start}, {window_end}] min | 阈值 CV<={cv_max:g}, "
          f"极差<={range_max:g} mmol/L")
    print(f"  稳态通过: {n_pass}/{n_total} 只动物")

    if n_pass < n_total:
        failed = qc_df[~qc_df["passed"]]
        print(f"  [WARN] {n_total - n_pass} 只动物未达血糖稳态：")
        for _, row in failed.head(12).iterrows():
            print(f"    {row['mouse_id']:>6} ({row['condition']}): "
                  f"mean={row['glucose_mean']}, CV={row['glucose_cv']}, "
                  f"range={row['glucose_range_mmol_l']}, "
                  f"drift={row['drift_mmol_l_per_h']:+} mmol/L/h "
                  f"| {row['fail_reasons']}")
        if len(failed) > 12:
            print(f"    ... 另有 {len(failed) - 12} 只，详见 infusion_steady_state_qc.csv")
        print("  [!] MFA 假设采样时代谢与同位素已达稳态。未达稳态的动物，其拟合通量")
        print("      是一段漂移过程的平均，不对应任何真实生理状态，组间比较需谨慎。")
        if on_fail == "exclude":
            print("  [!] on_fail=exclude 已声明，但样本级排除需要动物编号与 MFA 样本名")
            print("      的对应关系（当前 complete_result_name 不含动物编号），暂按 flag 处理。")

    by_condition = (
        qc_df.groupby("condition")
        .agg(n_mice=("mouse_id", "nunique"),
             n_passed=("passed", "sum"),
             mean_cv=("glucose_cv", "mean"),
             mean_drift=("drift_mmol_l_per_h", "mean"))
        .reset_index()
    )
    by_condition.to_csv(save_dir / "infusion_steady_state_by_condition.csv", index=False)
    print(f"  Saved all_tissues/qc/infusion_steady_state_qc.csv")
    return qc_df
