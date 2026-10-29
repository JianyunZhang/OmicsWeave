#!/usr/bin/env python3
"""
Re-run Differential Expression Analysis (config-driven)
========================================================

Re-runs PyDESeq2 differential expression analysis for all pairwise
comparisons defined in the experiment YAML config.

All experiment-specific values (mouse groups, tissues, comparisons)
are loaded from config.yaml / groups.yaml. No hardcoded experiment data.

Usage:
  python rerun_de_analysis.py --experiment <id>
  python rerun_de_analysis.py --config <path/to/config.yaml>
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd

from .experiment_config import (
    load_transcriptomics_config, config_arg_help, resolve_comparisons, get_cond_map,
)
from .run_wilcoxon_de import run_wilcoxon_analysis
from .gene_annotate import annotate_csv_file
from .rnaseq_pipeline import run_pydeseq2_analysis

# ============================================================================
# CONFIGURATION — all values loaded from experiment YAML config.
# ============================================================================
# Placeholders — will be overwritten by _apply_config()
OUTPUT_BASE = ""
DE_DIR = Path("")
MOUSE_GROUPS: dict[str, str] = {}
SAMPLE_BATCHES: dict[str, str] = {}
TISSUES: list[str] = []
COMPARISONS: list[tuple[str, str, str]] = []
COND_MAP: dict[str, str] = {}


def _apply_config(experiment: str | None, config_path: str | None, outdir: str | None) -> dict:
    cfg = load_transcriptomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )
    global OUTPUT_BASE, DE_DIR, MOUSE_GROUPS, SAMPLE_BATCHES, TISSUES, COMPARISONS, COND_MAP

    OUTPUT_BASE = os.path.abspath(str(cfg.output_dir))
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = OUTPUT_BASE
    DE_DIR = Path(OUTPUT_BASE) / "04_de_results"
    MOUSE_GROUPS = dict(cfg.mouse_groups)
    SAMPLE_BATCHES = dict(cfg.sample_batches)
    TISSUES = list(cfg.tissues)

    # Comparison list from config, fallback to all pairwise comparisons
    COMPARISONS = resolve_comparisons(cfg)
    COND_MAP = get_cond_map(cfg.conditions)

    return {
        "experiment_id": cfg.experiment_id,
        "config_path": str(cfg.config_path),
        "output_dir": OUTPUT_BASE,
    }


def flush_print(*args, **kwargs):
    """Print with immediate flush."""
    print(*args, **kwargs, flush=True)


def get_sample_metadata(sample_ids):
    """
    Build metadata DataFrame from sample IDs.
    Sample IDs look like: '52_liver', '63_brain', etc.

    批次列必须在此构造：主流程靠 metadata 里的 ``batch`` 列触发
    ``~batch + condition`` 设计与多批次守卫。若本入口不带批次信息，
    重跑 DE 就会绕过守卫、对多批次数据做无校正分析。
    """
    meta_rows = []
    for sid in sample_ids:
        parts = sid.split('_', 1)
        mouse_id = parts[0]
        tissue = parts[1] if len(parts) > 1 else 'unknown'
        group = MOUSE_GROUPS.get(sid) or MOUSE_GROUPS.get(mouse_id) or 'Unknown'
        row = {
            'sample': sid,
            'mouse_id': mouse_id,
            'tissue': tissue,
            'condition': group,
        }
        if SAMPLE_BATCHES:
            row['batch'] = SAMPLE_BATCHES.get(sid) or SAMPLE_BATCHES.get(mouse_id)
        meta_rows.append(row)

    meta_df = pd.DataFrame(meta_rows).set_index('sample')
    if SAMPLE_BATCHES:
        missing = meta_df['batch'].isna()
        if missing.any():
            raise ValueError(
                "sample_batches is configured but missing samples: "
                + ", ".join(meta_df.index[missing].astype(str))
            )
    return meta_df


# run_pydeseq2_analysis 统一由 rnaseq_pipeline 提供（见文件顶部 import）。
# 此处原有一份副本，缺少设计矩阵秩校验、批次元数据列、LFC 收缩与多组学导出，
# 会让主流程的修复只生效一半——已移除以保证两个入口行为一致。


def main():
    parser = argparse.ArgumentParser(
        description="Re-run DE analysis driven by experiment config.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None, help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None, help="Override output_dir (transcriptomics_results root)")
    args = parser.parse_args()

    cfg_info = _apply_config(args.experiment, args.config, args.outdir)
    start_time = datetime.now()

    flush_print("\n" + "=" * 70)
    flush_print("  Re-run Differential Expression Analysis (config-driven)")
    flush_print("=" * 70)
    flush_print(f"Start time: {start_time}")
    flush_print(f"Experiment: {cfg_info.get('experiment_id')}")
    flush_print(f"Config:     {cfg_info.get('config_path')}")
    flush_print(f"\nGroup Mapping (from groups.yaml):")
    for code, display in sorted(COND_MAP.items()):
        mice = sorted(mid for mid, grp in MOUSE_GROUPS.items() if grp == code)
        flush_print(f"  {code} ({display}): mice {mice}")
    flush_print(f"\nComparisons (auto-generated):")
    for num, den, comp_name in COMPARISONS:
        flush_print(f"  {comp_name}: {COND_MAP.get(num, num)} vs {COND_MAP.get(den, den)}")

    # Load gene counts matrix
    gene_counts_path = DE_DIR / "gene_counts_matrix.csv"
    if not gene_counts_path.exists():
        flush_print(f"\nERROR: Gene counts matrix not found: {gene_counts_path}")
        return

    flush_print(f"\nLoading gene counts matrix: {gene_counts_path}")
    gene_counts = pd.read_csv(gene_counts_path, index_col=0)
    flush_print(f"  Loaded: {gene_counts.shape[1]} samples, {gene_counts.shape[0]} genes")

    # Build correct metadata
    sample_ids = gene_counts.columns.tolist()
    metadata_df = get_sample_metadata(sample_ids)

    # Save corrected metadata
    metadata_path = DE_DIR / "sample_metadata_corrected.csv"
    metadata_df.to_csv(metadata_path)
    flush_print(f"\nCorrected metadata saved: {metadata_path}")

    # Print group summary
    flush_print(f"\nGroup summary:")
    for grp in sorted(metadata_df['condition'].unique()):
        samples = metadata_df[metadata_df['condition'] == grp].index.tolist()
        mice = sorted(metadata_df[metadata_df['condition'] == grp]['mouse_id'].unique())
        flush_print(f"  {grp}: {len(samples)} samples (mice: {', '.join(mice)})")

    # Run DE analysis per tissue
    flush_print("\n" + "=" * 70)
    flush_print("Differential Expression Analysis per Tissue")
    flush_print("=" * 70)

    for tissue in TISSUES:
        flush_print(f"\n{'='*50}")
        flush_print(f"Tissue: {tissue}")
        flush_print(f"{'='*50}")

        tissue_samples = [s for s in gene_counts.columns
                          if s.endswith(f'_{tissue}')]

        if len(tissue_samples) < 2:
            flush_print(f"  Skipping {tissue}: too few samples ({len(tissue_samples)})")
            continue

        tissue_de_dir = DE_DIR / tissue
        tissue_de_dir.mkdir(parents=True, exist_ok=True)

        tissue_counts = gene_counts[tissue_samples]
        tissue_meta = metadata_df.loc[tissue_samples]

        flush_print(f"  Samples: {len(tissue_samples)}")
        for grp in tissue_meta['condition'].unique():
            n = (tissue_meta['condition'] == grp).sum()
            flush_print(f"    {grp}: {n}")

        # Run each comparison
        for numerator, denominator, comp_name in COMPARISONS:
            comp_samples = tissue_meta[
                tissue_meta['condition'].isin([numerator, denominator])
            ].index.tolist()

            if len(comp_samples) < 2:
                flush_print(f"  [{comp_name}] Skipping: too few samples")
                continue

            n_num = (tissue_meta.loc[comp_samples, 'condition'] == numerator).sum()
            n_den = (tissue_meta.loc[comp_samples, 'condition'] == denominator).sum()

            if n_num < 1 or n_den < 1:
                flush_print(f"  [{comp_name}] Skipping: need >= 1 per group")
                continue

            flush_print(f"\n  --- {comp_name} ({tissue}) ---")
            flush_print(f"    {numerator}: {n_num}, {denominator}: {n_den}")

            de_done = False
            has_batch_design = (
                'batch' in tissue_meta.columns
                and tissue_meta.loc[comp_samples, 'batch'].nunique() > 1
            )
            design_formula = 'batch + condition' if has_batch_design else 'condition'

            # Try PyDESeq2 first (needs >= 2 per group)
            if n_num >= 2 and n_den >= 2:
                try:
                    de_results = run_pydeseq2_analysis(
                        tissue_counts[comp_samples],
                        tissue_meta.loc[comp_samples],
                        design_formula=design_formula,
                        contrast=('condition', numerator, denominator),
                    )

                    if de_results is not None:
                        de_results = de_results.reset_index()
                        de_results = de_results.rename(columns={'index': 'gene'})

                        out_file = tissue_de_dir / f"de_{comp_name}.csv"
                        de_results.to_csv(out_file, index=False)
                        flush_print(f"  Saved: {out_file}")
                        # 添加数据库注释
                        try:
                            annotate_csv_file(str(out_file), gene_col='gene')
                        except Exception as e:
                            flush_print(f"    [WARN] DE annotation skipped: {e}")
                        de_done = True
                except Exception as e:
                    flush_print(f"  PyDESeq2 failed for {comp_name}: {e}")
                    import traceback
                    traceback.print_exc()
                    if has_batch_design:
                        flush_print(
                            "  [BLOCKED] Explicit multi-batch design cannot be fit; "
                            "unadjusted Wilcoxon fallback is disabled."
                        )
                        continue

            # Fallback: Wilcoxon rank-sum (works with n >= 1 per group)
            if not de_done:
                if has_batch_design:
                    flush_print(
                        "  [BLOCKED] Multi-batch comparison requires an identifiable "
                        "batch-adjusted model; unadjusted fallback is disabled."
                    )
                    continue
                flush_print(f"  [Fallback] Using Wilcoxon rank-sum for {comp_name} ...")
                try:
                    de_results = run_wilcoxon_analysis(
                        tissue_counts[comp_samples],
                        tissue_meta.loc[comp_samples],
                        numerator=numerator,
                        denominator=denominator,
                    )
                    if de_results is not None:
                        out_file = tissue_de_dir / f"de_{comp_name}_wilcoxon.csv"
                        de_results.to_csv(out_file)
                        flush_print(f"  Saved (Wilcoxon): {out_file}")
                        # 添加数据库注释
                        try:
                            annotate_csv_file(str(out_file))
                        except Exception as e:
                            flush_print(f"    [WARN] DE annotation skipped: {e}")
                except Exception as e:
                    flush_print(f"  Wilcoxon also failed for {comp_name}: {e}")

    # Summary
    end_time = datetime.now()
    duration = end_time - start_time

    flush_print("\n" + "=" * 70)
    flush_print("ANALYSIS COMPLETE")
    flush_print("=" * 70)
    flush_print(f"Duration: {duration}")
    flush_print(f"\nOutput directory: {DE_DIR}")
    flush_print(f"\nComparisons performed:")
    for num, den, comp_name in COMPARISONS:
        flush_print(f"  - {comp_name}: {COND_MAP.get(num, num)} vs {COND_MAP.get(den, den)}")


if __name__ == "__main__":
    main()
