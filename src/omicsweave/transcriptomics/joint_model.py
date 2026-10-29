"""
跨批次联合建模（merge_runs + ~batch + condition）
==================================================

现有整合走的是「各批次独立跑完 DE → 对 log2FC 做随机效应 Meta」这条路
（``transcriptomics_integration.py``）。当批次内分组平衡时，**联合模型的
统计效力显著高于事后 Meta**：DESeq2 把批次作为固定效应吸收掉，并在全部
样本上共享离散度估计，而不是分别拟合两次再池化。

本模块提供另一条路径：

    runs/<protocol>/<batch>/04_de_results/gene_counts_matrix.csv
      ──merge_runs──►  合并计数矩阵 + 合并元数据（batch 列来自 run 目录）
      ──设计可识别性预检──►  拒绝混杂设计
      ──PyDESeq2 ~batch + condition──►  integration/joint_model/...

两条路径并不互斥：联合模型给出效力更高的主结论，Meta 分析给出批次间
异质性（I²）——后者是联合模型看不到的信息。整合流程会同时产出两者，
并在报告中标注哪条是主路径。

**关键前提 1：批次与分组不能完全混杂。**
若批次 1 全是 HFD、批次 2 全是 LFD，没有任何算法能把生物学差异和技术差异
分开——此时联合模型会被 ``assess_joint_design()`` 拒绝，只能退回 Meta 分析
并明确标注局限。

**关键前提 2：批次效应不能波及过半转录组。**
仿真扫描（批次影响基因占比 10%/25%/50%/75%）显示，当受影响基因接近 50% 时，
DESeq2 的 median-of-ratios 大小因子估计正处于崩溃点——一半基因的比值为 r、
另一半为 1，中位数在两者间跳变，导致批次内大小因子的变异系数从 ~0.013
恶化到 ~0.083，并注入**样本特异**的缩放误差。这类误差不是纯批次结构，
``~batch`` 项吸收不掉，会表现为假阳性。

  批次影响占比   批次内 size factor CV   纯无效应基因假阳性
       10%              0.015                3/340
       25%              0.012                2/300
       50%              0.083               78/200   ← 崩溃
       75%              0.052               18/100

  现实中这种情形出现在整批 RNA 严重降解、或建库方式不同的批次上。
  ``05_result_plot/3_pca_clustering/pc_covariate_association.csv`` 里
  批次在高方差主成分上极显著，就是预警信号；此时应先核查批次质量，
  而不是指望统计模型救回来。

用法::

    python joint_model.py --experiment <id>
    # 或经由整合入口（推荐）：
    python rnaseq_pipeline.py --experiment <id> --all-batches --integrate
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


from .experiment_config import (  # type: ignore
    load_transcriptomics_config,
    config_arg_help,
    discover_experiment_protocols,
    generate_comparisons,
    load_integration_config,
    resolve_comparisons,
    resolve_output_root,
)
from . import provenance as provenance_mod  # type: ignore


# 合并时对基因集的处理策略
GENE_POLICY_INTERSECTION = "intersection"
GENE_POLICY_UNION = "union"

# 样本 ID 跨批次重复时的处理策略
DUPLICATE_FAIL = "fail"
DUPLICATE_SUFFIX = "suffix"

# 设计可识别性的判定结论
ROUTE_JOINT = "joint_model"
ROUTE_META = "meta_analysis_only"
ROUTE_NONE = "insufficient_data"


def _flush_print(*args, **kwargs) -> None:
    print(*args, **kwargs, flush=True)


def _de_core():
    """
    延迟导入 rnaseq_pipeline 里的 DE 实现。

    ``rnaseq_pipeline`` → ``transcriptomics_integration`` → ``joint_model``
    已构成一条导入链；在模块顶层反向导入 rnaseq_pipeline 会形成循环。
    放在函数内部导入可以打断它，且调用时机总在模块加载完成之后。
    """
    from . import rnaseq_pipeline  # type: ignore
    return rnaseq_pipeline


# ---------------------------------------------------------------------------
# 1. 发现并读取各批次产物
# ---------------------------------------------------------------------------

@dataclass
class RunOutput:
    """一次批次运行的产物位置。"""
    protocol_id: str
    batch_id: str
    run_dir: Path
    counts_path: Path
    metadata_path: Path
    provenance: dict[str, Any] | None = None

    @property
    def available(self) -> bool:
        return self.counts_path.exists() and self.metadata_path.exists()


def discover_run_outputs(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
) -> dict[str, RunOutput]:
    """定位各批次的基因计数矩阵与样本元数据。"""
    found: dict[str, RunOutput] = {}
    for batch_id in batch_ids:
        run_dir = Path(output_root) / "runs" / protocol_id / batch_id
        found[batch_id] = RunOutput(
            protocol_id=protocol_id,
            batch_id=batch_id,
            run_dir=run_dir,
            counts_path=run_dir / "04_de_results" / "gene_counts_matrix.csv",
            metadata_path=run_dir / "sample_metadata.csv",
            provenance=provenance_mod.load_run_provenance(run_dir),
        )
    return found


def load_run_matrix(run: RunOutput) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    读取单批次的计数矩阵与元数据。

    计数矩阵可能被 ``annotate_csv_file`` 追加过注释列（ensembl_gene_id 等），
    因此只保留出现在元数据索引里的列——以元数据为准而非猜测列名。
    """
    metadata = pd.read_csv(run.metadata_path, index_col=0)
    metadata.index = metadata.index.astype(str)

    counts = pd.read_csv(run.counts_path, index_col=0, low_memory=False)
    counts.index = counts.index.astype(str)

    sample_cols = [c for c in counts.columns if str(c) in set(metadata.index)]
    if not sample_cols:
        raise ValueError(
            f"{run.counts_path} 与 {run.metadata_path} 没有共同样本；"
            "两者可能来自不同的运行。"
        )
    counts = counts[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    metadata = metadata.loc[[s for s in metadata.index if s in set(sample_cols)]]

    # 批次标识以 run 目录为准；若该批次内部还配置了 sample_batches，
    # 则组合成 "<batch_id>:<内部批次>"，避免两层批次语义互相覆盖。
    if "batch" in metadata.columns and metadata["batch"].notna().any():
        metadata["batch"] = [
            f"{run.batch_id}:{inner}" if pd.notna(inner) else run.batch_id
            for inner in metadata["batch"]
        ]
    else:
        metadata["batch"] = run.batch_id
    metadata["source_batch"] = run.batch_id
    metadata["source_protocol"] = run.protocol_id
    metadata["orig_sample"] = metadata.index.astype(str)

    return counts, metadata


def merge_run_counts(
    runs: dict[str, RunOutput],
    gene_policy: str = GENE_POLICY_INTERSECTION,
    duplicate_policy: str = DUPLICATE_FAIL,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """
    把多个批次的基因计数矩阵合并为一个联合矩阵。

    Parameters
    ----------
    gene_policy : {"intersection", "union"}
        ``intersection``（默认）只保留所有批次共有的基因。参考一致时两者
        应当等价；若交集显著小于任一批次，说明参考版本存在差异——这正是
        溯源断言要拦下的情况。``union`` 会对缺失基因补 0，等于制造假的
        「零表达」，仅在明确知情时使用。
    duplicate_policy : {"fail", "suffix"}
        样本 ID 跨批次重复时的处理。默认 ``fail``：重复通常意味着桥接样本
        或配置错误，两者都不应被静默当作独立生物学重复（会低估方差、
        虚高显著性）。``suffix`` 保留两份并加批次后缀，同时告警。

    Returns
    -------
    (counts, metadata, info)
    """
    usable = {bid: run for bid, run in runs.items() if run.available}
    missing = sorted(set(runs) - set(usable))
    if len(usable) < 2:
        raise ValueError(
            "联合建模至少需要 2 个已完成的批次；当前可用: "
            f"{sorted(usable)}，缺少产物: {missing}"
        )

    counts_parts: list[pd.DataFrame] = []
    meta_parts: list[pd.DataFrame] = []
    per_batch_genes: dict[str, int] = {}

    for batch_id in sorted(usable):
        counts, metadata = load_run_matrix(usable[batch_id])
        per_batch_genes[batch_id] = int(counts.shape[0])
        counts_parts.append(counts)
        meta_parts.append(metadata)

    # ---- 样本 ID 冲突检测 ----
    all_samples: list[str] = []
    for part in counts_parts:
        all_samples.extend(str(c) for c in part.columns)
    duplicates = sorted({s for s in all_samples if all_samples.count(s) > 1})

    renamed: dict[str, str] = {}
    if duplicates:
        if duplicate_policy == DUPLICATE_FAIL:
            raise ValueError(
                f"以下样本 ID 在多个批次中重复出现: {duplicates[:10]}"
                f"{' …' if len(duplicates) > 10 else ''}\n"
                "若这些是桥接样本（同一份 RNA 分装到多个批次），它们是技术重复，"
                "不能当作独立生物学重复直接进入联合模型——那会低估方差、虚高显著性。\n"
                "处理方式：(a) 用 duplicate_policy='suffix' 保留并自行在设计中处理；"
                "(b) 只保留其中一个批次的副本；"
                "(c) 把桥接样本用于批次效应评估而不纳入 DE。"
            )
        # suffix 策略：加批次后缀区分
        for i, batch_id in enumerate(sorted(usable)):
            part = counts_parts[i]
            meta = meta_parts[i]
            mapping = {
                str(c): f"{c}__{batch_id}"
                for c in part.columns if str(c) in set(duplicates)
            }
            if mapping:
                counts_parts[i] = part.rename(columns=mapping)
                meta_parts[i] = meta.rename(index=mapping)
                renamed.update({v: k for k, v in mapping.items()})

    # ---- 基因集对齐 ----
    gene_sets = [set(part.index) for part in counts_parts]
    intersection = set.intersection(*gene_sets)
    union = set.union(*gene_sets)

    if gene_policy == GENE_POLICY_UNION:
        genes = sorted(union)
        counts_parts = [part.reindex(index=genes).fillna(0) for part in counts_parts]
    else:
        genes = sorted(intersection)
        if not genes:
            raise ValueError("各批次基因集无交集；参考版本很可能不同。")
        counts_parts = [part.loc[genes] for part in counts_parts]

    merged_counts = pd.concat(counts_parts, axis=1)
    merged_metadata = pd.concat(meta_parts, axis=0)
    merged_metadata.index = merged_metadata.index.astype(str)
    merged_counts = merged_counts.round(0).astype(int)

    # 列顺序与元数据索引严格对齐，杜绝错位
    merged_metadata = merged_metadata.loc[
        [s for s in merged_metadata.index if s in set(merged_counts.columns)]
    ]
    merged_counts = merged_counts[list(merged_metadata.index)]

    info = {
        "n_batches": len(usable),
        "batch_ids": sorted(usable),
        "missing_batches": missing,
        "n_samples": int(merged_counts.shape[1]),
        "n_genes": int(merged_counts.shape[0]),
        "gene_policy": gene_policy,
        "genes_per_batch": per_batch_genes,
        "n_genes_intersection": len(intersection),
        "n_genes_union": len(union),
        "duplicate_samples": duplicates,
        "duplicate_policy": duplicate_policy,
        "renamed_samples": renamed,
    }

    # 交集远小于并集 = 参考不一致的强信号
    if union and len(intersection) / len(union) < 0.9:
        info["gene_overlap_warning"] = (
            f"基因交集仅占并集的 {len(intersection) / len(union):.1%}，"
            "各批次参考版本很可能不同。"
        )

    return merged_counts, merged_metadata, info


# ---------------------------------------------------------------------------
# 2. 设计可识别性预检
# ---------------------------------------------------------------------------

def is_design_identifiable(
    metadata_sub: pd.DataFrame,
    design_formula: str,
) -> tuple[bool, str]:
    """
    在不拟合模型的前提下判断设计是否可识别。

    复用主流程的 ``_validate_de_design()``（设计矩阵秩检验），把异常
    转成布尔判定，用于生成预检报告而不是直接中断。
    """
    try:
        _de_core()._validate_de_design(metadata_sub, design_formula)
    except ValueError as exc:
        return False, str(exc)
    return True, ""


def assess_joint_design(
    metadata: pd.DataFrame,
    comparisons: list[tuple[str, str, str]],
    min_samples_per_group: int = 2,
) -> pd.DataFrame:
    """
    逐（组织 × 比较）评估联合模型是否可行，生成预检报告。

    这一步必须在拟合之前完成：批次与分组完全混杂时，任何算法都无法分离
    技术差异与生物学差异，此时应当明确拒绝而不是给出看似合理的结果。

    Returns
    -------
    pd.DataFrame
        每行一个（tissue, comparison），含各批次×条件的样本数、是否可识别、
        推荐路径与原因。
    """
    rows: list[dict[str, Any]] = []
    tissues = sorted(metadata["tissue"].dropna().astype(str).unique())

    for tissue in tissues:
        tissue_meta = metadata[metadata["tissue"].astype(str) == tissue]
        for numerator, denominator, comp_name in comparisons:
            sub = tissue_meta[
                tissue_meta["condition"].isin([numerator, denominator])
            ]
            n_num = int((sub["condition"] == numerator).sum())
            n_den = int((sub["condition"] == denominator).sum())
            batches = sorted(sub["batch"].dropna().astype(str).unique())

            # batch × condition 交叉表 —— 混杂与否一眼可见
            if len(sub) > 0:
                crosstab = pd.crosstab(
                    sub["batch"].astype(str), sub["condition"].astype(str)
                )
                crosstab_str = ";".join(
                    f"{b}:{'/'.join(f'{c}={crosstab.loc[b, c]}' for c in crosstab.columns)}"
                    for b in crosstab.index
                )
            else:
                crosstab_str = ""

            row: dict[str, Any] = {
                "tissue": tissue,
                "comparison": comp_name,
                "numerator": numerator,
                "denominator": denominator,
                "n_numerator": n_num,
                "n_denominator": n_den,
                "n_samples": int(len(sub)),
                "n_batches": len(batches),
                "batch_ids": ";".join(batches),
                "batch_condition_table": crosstab_str,
            }

            if n_num < min_samples_per_group or n_den < min_samples_per_group:
                row.update({
                    "design_formula": "",
                    "identifiable": False,
                    "residual_df": 0,
                    "route": ROUTE_NONE,
                    "reason": (
                        f"每组至少需要 {min_samples_per_group} 个样本 "
                        f"(实际 {n_num}/{n_den})"
                    ),
                })
                rows.append(row)
                continue

            if len(batches) < 2:
                # 只有一个批次：不需要批次项，也谈不上「联合」
                identifiable, err = is_design_identifiable(sub, "condition")
                row.update({
                    "design_formula": "~condition",
                    "identifiable": identifiable,
                    "residual_df": max(0, len(sub) - 2),
                    "route": ROUTE_JOINT if identifiable else ROUTE_NONE,
                    "reason": err or "单批次，无需批次校正",
                })
                rows.append(row)
                continue

            identifiable, err = is_design_identifiable(sub, "batch + condition")
            # 残差自由度 = n - (批次水平-1) - (条件水平-1) - 截距
            residual_df = len(sub) - (len(batches) - 1) - 1 - 1
            row.update({
                "design_formula": "~batch + condition",
                "identifiable": identifiable,
                "residual_df": int(residual_df),
            })

            if not identifiable:
                row.update({
                    "route": ROUTE_META,
                    "reason": (
                        "批次与分组混杂，联合模型不可识别；"
                        "退回批次内效应 + Meta 分析。原始信息: " + err
                    ),
                })
            elif residual_df < 1:
                row.update({
                    "route": ROUTE_META,
                    "reason": f"残差自由度不足 ({residual_df})，无法稳定估计离散度",
                })
            else:
                row.update({
                    "route": ROUTE_JOINT,
                    "reason": "批次内分组平衡，设计可识别",
                })
            rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3. 联合 DE
# ---------------------------------------------------------------------------

def run_joint_de(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design_report: pd.DataFrame,
    out_dir: Path,
    protocol_id: str,
    shrink_lfc: bool = True,
) -> dict[tuple[str, str], pd.DataFrame]:
    """
    对预检通过的（组织 × 比较）拟合 ``~batch + condition`` 联合模型。

    预检判定为 ``meta_analysis_only`` 的组合会被跳过——那些组合在联合模型
    下不可识别，硬拟合只会产出无意义的系数。
    """
    core = _de_core()
    results: dict[tuple[str, str], pd.DataFrame] = {}

    eligible = design_report[design_report["route"] == ROUTE_JOINT]
    if eligible.empty:
        _flush_print("  [joint] 没有任何组合通过设计预检，跳过联合建模。")
        return results

    for row in eligible.itertuples():
        tissue = row.tissue
        comp_name = row.comparison
        numerator = row.numerator
        denominator = row.denominator

        sub_meta = metadata[
            (metadata["tissue"].astype(str) == tissue)
            & (metadata["condition"].isin([numerator, denominator]))
        ]
        sub_counts = counts[list(sub_meta.index)]
        design_formula = str(row.design_formula).lstrip("~")

        _flush_print(
            f"\n  [joint] {tissue} / {comp_name} "
            f"(n={len(sub_meta)}, batches={row.n_batches}, ~{design_formula})"
        )

        try:
            de_results = core.run_pydeseq2_analysis(
                sub_counts,
                sub_meta,
                design_formula=design_formula,
                contrast=("condition", numerator, denominator),
                shrink_lfc=shrink_lfc,
            )
        except Exception as exc:  # noqa: BLE001
            _flush_print(f"  [joint][FAIL] {tissue}/{comp_name}: "
                         f"{type(exc).__name__}: {exc}")
            continue

        if de_results is None or de_results.empty:
            continue

        # 标注来源，避免与单批次 / Meta 结果混淆
        de_results = de_results.copy()
        de_results["source_mode"] = "joint_count_model_batch_adjusted"
        de_results["integration_route"] = ROUTE_JOINT
        de_results["protocol_id"] = protocol_id

        tissue_dir = out_dir / tissue
        tissue_dir.mkdir(parents=True, exist_ok=True)
        out_file = tissue_dir / f"de_{comp_name}_joint.csv"
        de_results.to_csv(out_file)
        _flush_print(f"  [joint] saved: {out_file}")

        # 多组学稳定接口（列名契约与单批次一致）
        try:
            core.export_multiomics_csv(de_results, str(tissue_dir), comp_name)
        except Exception as exc:  # noqa: BLE001
            _flush_print(f"  [joint][WARN] multi-omics export skipped: {exc}")

        results[(tissue, comp_name)] = de_results

    return results


# ---------------------------------------------------------------------------
# 4. 编排
# ---------------------------------------------------------------------------

def run_joint_model_for_protocol(
    output_root: Path,
    protocol_id: str,
    batch_ids: list[str],
    comparisons: list[tuple[str, str, str]],
    joint_cfg: dict[str, Any],
    check_provenance: bool = True,
    provenance_policy: str = "fail",
) -> dict[str, Any]:
    """对单个 protocol 下的多个批次执行 merge → 预检 → 联合 DE。"""
    out_dir = Path(output_root) / "integration" / "joint_model" / protocol_id
    out_dir.mkdir(parents=True, exist_ok=True)

    summary: dict[str, Any] = {
        "protocol_id": protocol_id,
        "batch_ids": list(batch_ids),
        "status": "unknown",
    }

    runs = discover_run_outputs(output_root, protocol_id, batch_ids)

    # ---- 参考/技术参数一致性（与 Meta 路径共用同一套断言）----
    if check_provenance:
        comparison = provenance_mod.compare_provenance(
            {bid: run.provenance for bid, run in runs.items()}
        )
        report = provenance_mod.format_provenance_report(comparison)
        summary["provenance"] = comparison
        summary["provenance_report"] = report
        if not comparison["consistent"]:
            if provenance_policy == "fail":
                _flush_print(f"  [joint][BLOCKED] 参考/技术参数不一致，拒绝合并:\n{report}")
                summary["status"] = "blocked_provenance_mismatch"
                return summary
            _flush_print(f"  [joint][WARN] 参考/技术参数不一致，"
                         f"但策略为 '{provenance_policy}'，继续:\n{report}")
        else:
            _flush_print(f"  [joint] {report}")

    # ---- 合并 ----
    try:
        counts, metadata, merge_info = merge_run_counts(
            runs,
            gene_policy=str(joint_cfg.get("gene_policy", GENE_POLICY_INTERSECTION)),
            duplicate_policy=str(joint_cfg.get("duplicate_policy", DUPLICATE_FAIL)),
        )
    except ValueError as exc:
        _flush_print(f"  [joint][BLOCKED] 合并失败: {exc}")
        summary["status"] = "merge_failed"
        summary["error"] = str(exc)
        return summary

    summary["merge_info"] = merge_info
    _flush_print(
        f"  [joint] merged {merge_info['n_samples']} samples × "
        f"{merge_info['n_genes']} genes from {merge_info['n_batches']} batches"
    )
    if merge_info.get("gene_overlap_warning"):
        _flush_print(f"  [joint][WARN] {merge_info['gene_overlap_warning']}")

    counts.to_csv(out_dir / "merged_gene_counts.csv")
    metadata.to_csv(out_dir / "merged_sample_metadata.csv")
    with open(out_dir / "merge_info.json", "w", encoding="utf-8") as handle:
        json.dump(merge_info, handle, indent=2, ensure_ascii=False)

    # ---- 设计预检 ----
    available_conditions = set(metadata["condition"].dropna().astype(str))
    usable_comparisons = [
        (n, d, c) for n, d, c in comparisons
        if n in available_conditions and d in available_conditions
    ]
    if not usable_comparisons:
        _flush_print("  [joint] 合并数据中没有可执行的比较（条件不匹配）。")
        summary["status"] = "no_usable_comparisons"
        return summary

    design_report = assess_joint_design(
        metadata,
        usable_comparisons,
        min_samples_per_group=int(joint_cfg.get("min_samples_per_group", 2)),
    )
    design_path = out_dir / "joint_design_report.csv"
    design_report.to_csv(design_path, index=False)
    summary["design_report_path"] = str(design_path)

    n_joint = int((design_report["route"] == ROUTE_JOINT).sum())
    n_meta = int((design_report["route"] == ROUTE_META).sum())
    n_none = int((design_report["route"] == ROUTE_NONE).sum())
    _flush_print(
        f"  [joint] 设计预检: {n_joint} 可联合建模, "
        f"{n_meta} 仅可 Meta, {n_none} 数据不足 → {design_path}"
    )
    for row in design_report[design_report["route"] != ROUTE_JOINT].itertuples():
        _flush_print(f"    [{row.route}] {row.tissue}/{row.comparison}: {row.reason}")

    summary["n_joint_eligible"] = n_joint
    summary["n_meta_only"] = n_meta
    summary["n_insufficient"] = n_none

    # ---- 联合 DE ----
    results = run_joint_de(
        counts, metadata, design_report, out_dir, protocol_id,
        shrink_lfc=bool(joint_cfg.get("shrink_lfc", True)),
    )
    summary["n_results"] = len(results)
    summary["results"] = {f"{t}/{c}": len(df) for (t, c), df in results.items()}
    summary["status"] = "completed" if results else "no_results"
    return summary


def run_joint_model(experiment_id: str, outdir: str | None = None) -> dict[str, Any]:
    """
    实验级入口：对所有启用的 protocol 执行跨批次联合建模。

    只有批次数 ≥2 的 protocol 才会进入——单批次没有「跨批次」可言，
    其结果已由常规流程产出。
    """
    integration_cfg = load_integration_config(experiment_id)
    joint_cfg = (integration_cfg.get("joint_model", {}) or {})

    if not joint_cfg.get("enabled", True):
        _flush_print("[INFO] integration.joint_model.enabled=false, skip joint modeling.")
        return {"status": "disabled"}

    output_root = Path(outdir) if outdir else resolve_output_root(experiment_id)
    protocols = discover_experiment_protocols(experiment_id)

    # 比较列表来自实验的活动配置
    try:
        cfg = load_transcriptomics_config(experiment_id=experiment_id)
        comparisons = resolve_comparisons(cfg)
    except Exception as exc:  # noqa: BLE001
        _flush_print(f"[WARN] 无法加载活动配置 ({exc})，将从合并元数据自动生成比较。")
        comparisons = []

    results: dict[str, Any] = {"experiment_id": experiment_id, "protocols": {}}

    for proto in protocols:
        if not proto.get("enabled", False):
            continue
        protocol_id = proto.get("protocol_id")
        if not protocol_id:
            continue
        batch_ids = [
            b["batch_id"] for b in proto.get("batches", [])
            if b.get("has_data", False) and b.get("batch_id")
        ]
        if len(batch_ids) < 2:
            _flush_print(
                f"[joint] {protocol_id}: 仅 {len(batch_ids)} 个批次，"
                "跨批次联合建模需要 ≥2 个批次，跳过。"
            )
            results["protocols"][protocol_id] = {
                "status": "skipped_single_batch",
                "batch_ids": batch_ids,
            }
            continue

        _flush_print(f"\n[joint] protocol: {protocol_id} ({len(batch_ids)} batches)")

        comps = comparisons
        if not comps:
            runs = discover_run_outputs(output_root, protocol_id, batch_ids)
            conditions: set[str] = set()
            for run in runs.values():
                if run.metadata_path.exists():
                    meta = pd.read_csv(run.metadata_path, index_col=0)
                    conditions.update(meta["condition"].dropna().astype(str))
            comps = generate_comparisons({c: {} for c in conditions})

        results["protocols"][protocol_id] = run_joint_model_for_protocol(
            output_root=output_root,
            protocol_id=protocol_id,
            batch_ids=batch_ids,
            comparisons=comps,
            joint_cfg=joint_cfg,
            check_provenance=bool(joint_cfg.get("check_provenance", True)),
            provenance_policy=str(
                integration_cfg.get("reference_mismatch_policy", "fail")
            ),
        )

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="跨批次联合建模（merge_runs + ~batch + condition）。",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", required=True,
                        help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--outdir", default=None,
                        help="Override output root (transcriptomics_results)")
    args = parser.parse_args()

    summary = run_joint_model(args.experiment, outdir=args.outdir)
    _flush_print("\n" + json.dumps(summary, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
