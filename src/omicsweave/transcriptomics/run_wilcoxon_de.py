"""
Wilcoxon (Mann-Whitney U) DE Analysis — config-driven
======================================================

Non-parametric alternative to PyDESeq2, suitable for experiments with
low biological replicates (n=1-2 per group).  Automatically used as
fallback when PyDESeq2 fails due to insufficient samples.

Can also be run standalone:
  python run_wilcoxon_de.py --experiment <id>
  python run_wilcoxon_de.py --config <path/to/config.yaml>
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .experiment_config import (
    load_transcriptomics_config, config_arg_help,
    resolve_comparisons, get_cond_map,
)
from .normalization import normalize_counts

# Placeholders — populated by _apply_config()
OUTPUT_BASE: str = ""
DE_DIR: Path = Path("")
MOUSE_GROUPS: dict[str, str] = {}
TISSUES: list[str] = []
COMPARISONS: list[tuple[str, str, str]] = []
COND_MAP: dict[str, str] = {}


def _apply_config(experiment: str | None, config_path: str | None, outdir: str | None) -> dict:
    cfg = load_transcriptomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )
    global OUTPUT_BASE, DE_DIR, MOUSE_GROUPS, TISSUES, COMPARISONS, COND_MAP

    OUTPUT_BASE = os.path.abspath(str(cfg.output_dir))
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = OUTPUT_BASE
    DE_DIR = Path(OUTPUT_BASE) / "04_de_results"
    MOUSE_GROUPS = dict(cfg.mouse_groups)
    TISSUES = list(cfg.tissues)
    COMPARISONS = resolve_comparisons(cfg)
    COND_MAP = get_cond_map(cfg.conditions)

    return {
        "experiment_id": cfg.experiment_id,
        "config_path": str(cfg.config_path),
        "output_dir": OUTPUT_BASE,
    }


def _validate_no_unadjusted_multibatch(metadata_df: pd.DataFrame) -> None:
    """Wilcoxon 入口不允许忽略显式多批次设计。"""
    if "batch" in metadata_df.columns and metadata_df["batch"].nunique() > 1:
        raise ValueError(
            "Unadjusted Wilcoxon is not valid for an explicit multi-batch design; "
            "fit ~batch + condition or use within-batch effects and Meta-analysis"
        )


def min_achievable_pvalue(n_case: int, n_control: int) -> float:
    """
    Mann-Whitney 精确检验在给定样本量下可能取到的**最小双侧 p 值**。

    共有 C(n1+n2, n1) 种秩分配，最极端的一种对应单侧 p = 1/C，
    双侧即 2/C。这个上限与数据无关——纯粹由样本量决定。

    实践含义：n=3 vs n=3 时 C(6,3)=20，最小双侧 p = 0.1。
    此时无论生物学效应多强，经 BH 校正后 padj 恒 ≥ 0.1，
    ``padj < 0.05`` 永远不可能成立。若不显式提示，流程会给出
    「0 个显著基因」这一看似生物学的结论，实则是方法学上限。
    """
    if n_case < 1 or n_control < 1:
        return float("nan")
    total_arrangements = math.comb(n_case + n_control, n_case)
    return min(1.0, 2.0 / total_arrangements)


def run_wilcoxon_analysis(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    numerator: str,
    denominator: str,
) -> pd.DataFrame | None:
    """
    Run Wilcoxon rank-sum test for each gene.

    Args:
        counts_df: Gene count matrix (genes x samples), genes as index.
        metadata_df: Sample metadata with 'condition' column.
        numerator: Condition name for numerator group.
        denominator: Condition name for denominator group.

    Returns:
        DataFrame with genes as index, columns: baseMean, log2FoldChange,
        pvalue, padj, significant, test_method, n_case, n_control,
        evidence_level, min_achievable_pvalue, normalization.

    Notes:
        Counts are library-size normalized (DESeq2 median-of-ratios) before
        testing.  Mann-Whitney has no immunity to sequencing-depth differences:
        on raw counts, two samples differing 2x in depth show a systematic
        log2FC offset across *every* gene.
    """
    from statsmodels.stats.multitest import multipletests

    _validate_no_unadjusted_multibatch(metadata_df)

    num_samples = metadata_df[metadata_df['condition'] == numerator].index.tolist()
    den_samples = metadata_df[metadata_df['condition'] == denominator].index.tolist()

    n_num = len(num_samples)
    n_den = len(den_samples)

    # 文库大小归一化 —— 必须在检验和 log2FC 计算之前
    comp_samples = num_samples + den_samples
    normalized_df, size_factors = normalize_counts(counts_df[comp_samples])

    # Determine evidence level
    if n_num < 2 or n_den < 2:
        test_method = "wilcoxon_effect_only"
        evidence_level = "effect_size_only"
    else:
        test_method = "wilcoxon"
        evidence_level = "statistical_test"

    # 样本量决定的检验分辨率上限（与数据无关）
    p_floor = min_achievable_pvalue(n_num, n_den)
    if evidence_level == "statistical_test" and np.isfinite(p_floor) and p_floor > 0.05:
        print(
            f"    [WARN] n={n_num} vs n={n_den}: Mann-Whitney 最小可达双侧 p = {p_floor:.3f}，"
            f"高于 0.05。本比较在方法学上不可能产生 padj<0.05 的基因——"
            f"「0 个显著基因」是样本量上限，不是生物学结论。",
            flush=True,
        )

    results = []
    for gene in normalized_df.index:
        vals_num = normalized_df.loc[gene, num_samples].values.astype(float)
        vals_den = normalized_df.loc[gene, den_samples].values.astype(float)

        base_mean = np.mean(np.concatenate([vals_num, vals_den]))
        log2fc = np.log2((np.mean(vals_num) + 1) / (np.mean(vals_den) + 1))

        if n_num >= 2 and n_den >= 2:
            try:
                _, pval = stats.mannwhitneyu(vals_num, vals_den, alternative="two-sided")
            except Exception:
                pval = np.nan
        else:
            # n=1 per group: only log2FC, no p-value
            pval = np.nan

        results.append({
            "gene": gene,
            "baseMean": base_mean,
            "log2FoldChange": log2fc,
            "pvalue": pval,
            "padj": np.nan,
        })

    res_df = pd.DataFrame(results).set_index("gene")

    # BH correction for valid p-values
    valid_mask = res_df["pvalue"].notna()
    if valid_mask.sum() > 0:
        try:
            _, padj_valid, _, _ = multipletests(
                res_df.loc[valid_mask, "pvalue"].values, method="fdr_bh"
            )
            res_df.loc[valid_mask, "padj"] = padj_valid
        except Exception:
            pass

    # Mark significance — effect_size_only results are NEVER significant
    if evidence_level == "effect_size_only":
        res_df["significant"] = False
    elif valid_mask.any():
        res_df["significant"] = (res_df["padj"] < 0.05) & (res_df["log2FoldChange"].abs() > 1)
    else:
        res_df["significant"] = False

    # Add metadata columns
    res_df["test_method"] = test_method
    res_df["n_case"] = n_num
    res_df["n_control"] = n_den
    res_df["evidence_level"] = evidence_level
    res_df["min_achievable_pvalue"] = p_floor
    res_df["normalization"] = "median_of_ratios"
    res_df["size_factor_range"] = (
        f"{size_factors.min():.3f}-{size_factors.max():.3f}" if len(size_factors) else ""
    )

    return res_df


# ---------------------------------------------------------------------------
# Standalone CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Wilcoxon (Mann-Whitney U) DE analysis for low-replicate experiments.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None,
                        help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None,
                        help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None,
                        help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()

    _apply_config(experiment=args.experiment, config_path=args.config, outdir=args.outdir)

    # Load gene counts
    gene_counts_path = os.path.join(OUTPUT_BASE, "04_de_results", "gene_counts_matrix.csv")
    if not os.path.exists(gene_counts_path):
        print(f"ERROR: gene_counts_matrix.csv not found at {gene_counts_path}")
        sys.exit(1)

    gene_counts = pd.read_csv(gene_counts_path, index_col=0)
    metadata_path = os.path.join(OUTPUT_BASE, "sample_metadata.csv")
    metadata_df = pd.read_csv(metadata_path, index_col=0)

    print(f"Loaded {len(gene_counts)} genes x {len(gene_counts.columns)} samples")

    for tissue in TISSUES:
        tissue_samples = [s for s in gene_counts.columns if s.endswith(f'_{tissue}')]
        if len(tissue_samples) < 2:
            continue

        tissue_counts = gene_counts[tissue_samples]
        tissue_meta = metadata_df.loc[tissue_samples]

        # 创建组织子目录
        tissue_de_dir = DE_DIR / tissue
        tissue_de_dir.mkdir(parents=True, exist_ok=True)

        for numerator, denominator, comp_name in COMPARISONS:
            comp_samples = tissue_meta[
                tissue_meta['condition'].isin([numerator, denominator])
            ].index.tolist()

            if len(comp_samples) < 2:
                continue

            print(f"\n  Wilcoxon: {tissue} / {comp_name}")

            res_df = run_wilcoxon_analysis(
                tissue_counts[comp_samples],
                tissue_meta.loc[comp_samples],
                numerator=numerator,
                denominator=denominator,
            )

            if res_df is None or len(res_df) == 0:
                continue

            sig = res_df[res_df["significant"]]
            up = sig[sig["log2FoldChange"] > 0]
            down = sig[sig["log2FoldChange"] < 0]
            print(f"  Genes tested: {len(res_df)}, "
                  f"Significant: {len(sig)} (Up: {len(up)}, Down: {len(down)})")

            out_file = tissue_de_dir / f"de_{comp_name}_wilcoxon.csv"
            res_df.to_csv(out_file)
            print(f"  Saved: {out_file}")

    print("\nWilcoxon DE analysis complete.")


if __name__ == "__main__":
    main()
