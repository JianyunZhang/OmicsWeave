#!/usr/bin/env python3
"""
Overview plots: data availability heatmap and differential feature summary.

Adapted from multiomics_integration_analysis.py lines 883-910 (overview) and
lines 390-440 (differential_summary).

Author: Claude Code
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from . import _config as cfg


def _truthy(series: pd.Series) -> pd.Series:
    """Parse a significance flag column that may be bool, int or string."""
    return series.map(lambda v: str(v).strip().lower() in {"true", "1", "yes"})


def plot_data_availability(
    tissues: list[str],
    mfa_tissues: set[str],
    metab_tissues: set[str],
    prot_tissues: set[str],
    tx_tissues: set[str],
    output_dir: Path | None = None,
    mfa_diff: dict | None = None,
    metab_diff: dict | None = None,
    prot_diff: dict | None = None,
    tx_de: dict | None = None,
    pval_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.58,
) -> pd.DataFrame:
    """Heatmap showing significant feature counts per omics per tissue.

    Parameters
    ----------
    tissues : list[str]
        All configured tissue names.
    mfa_tissues, metab_tissues, prot_tissues, tx_tissues : set[str]
        Tissues with each omics layer loaded.
    output_dir : Path, optional
        Override output directory.
    mfa_diff, metab_diff, prot_diff, tx_de : dict, optional
        Differential data dicts for counting significant features.
    pval_cutoff, log2fc_cutoff : float
        Significance thresholds for counting.

    Returns
    -------
    pd.DataFrame
        Data availability table with counts.
    """
    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    # Build count matrix: significant features per omics per tissue
    def _count_sig_mfa(tissue: str) -> int:
        """Count significant MFA fluxes for this tissue.

        Significance comes from the MFA pipeline's own ``is_significant``
        column when present. That column encodes the power-aware rule the MFA
        run actually applied: ``statistic: auto`` falls back to
        ``effect_size_only`` (|Hedges' g| + |log2FC|) whenever the rank-test
        p-value floor ``2/C(n1+n2, n1)`` or the BH floor makes p<0.05
        unreachable -- see omicsweave/common/nonparametric.py.

        Re-deriving significance from raw p-values here contradicts that
        verdict. At n=2-4 per group the Mann-Whitney floor alone exceeds 0.05,
        so a p-value criterion reports 0 for every tissue regardless of
        biology -- which is what this function used to do.

        The FDR / p-value branches remain for legacy MFA results predating
        the ``is_significant`` column.

        Returns -1 when the tissue has no MFA data at all (rendered as a blank
        cell), which is distinct from 0 = "data present, nothing significant".
        """
        if mfa_diff is None:
            return -1  # data not provided
        seen_tissue = False
        total = 0
        for comp_df in mfa_diff.values():
            if "tissue" not in comp_df.columns:
                continue
            sub = comp_df[comp_df["tissue"] == tissue]
            if sub.empty:
                continue
            seen_tissue = True
            if "is_significant" in sub.columns:
                mask = _truthy(sub["is_significant"])
            else:
                # Legacy gate criteria: FDR + Cohen's d
                fdr_col = "fdr_mannwhitney" if "fdr_mannwhitney" in sub.columns else None
                cd_col = "cohen_d" if "cohen_d" in sub.columns else None
                if fdr_col and cd_col:
                    mask = (sub[fdr_col].fillna(1.0) < 0.05) & (sub[cd_col].abs().fillna(0) > 0.8)
                else:
                    # Legacy fallback: pvalue + log2FC
                    pcol = "flux_pvalue" if "flux_pvalue" in sub.columns else "pvalue_mannwhitney"
                    fcol = "flux_log2fc" if "flux_log2fc" in sub.columns else "log2FC"
                    if pcol in sub.columns and fcol in sub.columns:
                        mask = (sub[pcol].fillna(1.0) < pval_cutoff) & (sub[fcol].abs().fillna(0) > log2fc_cutoff)
                    else:
                        mask = pd.Series(False, index=sub.index)
            total += int(mask.sum())
        return total if seen_tissue else -1

    def _count_sig_tx(tissue: str) -> int:
        """Count significant DE genes for this tissue across all comparisons.

        When p-values are NaN (effect-size-only mode, e.g. n=1/low-n),
        uses |log2FC| > 1.0 (2-fold) as fallback threshold.
        """
        if tx_de is None:
            return -1
        genes = set()
        seen_tissue = False
        for comp_dict in tx_de.values():
            if tissue not in comp_dict:
                continue
            seen_tissue = True
            df = comp_dict[tissue]
            padj_col = "padj" if "padj" in df.columns else "pvalue"
            fc_col = "log2fc" if "log2fc" in df.columns else "log2FoldChange"
            if "gene" not in df.columns or fc_col not in df.columns:
                continue
            has_pval = padj_col in df.columns and df[padj_col].notna().any()
            if has_pval:
                mask = (df[padj_col].fillna(1.0) < pval_cutoff) & (df[fc_col].abs().fillna(0) > log2fc_cutoff)
            else:
                # Effect-size-only: require |log2FC| > 1.0 (2-fold)
                mask = df[fc_col].abs() > 1.0
            genes.update(df.loc[mask, "gene"].dropna().tolist())
        # -1 = this tissue has no transcriptomics at all, matching the sentinel
        # used by the proteomics/metabolomics counters. Returning 0 here made
        # "no data" indistinguishable from "data present, nothing significant".
        return len(genes) if seen_tissue else -1

    def _count_sig_prot(tissue: str) -> int:
        """Count significant proteins for this tissue.

        Per-comparison logic:
        - If pval column has valid data: p < cutoff
        - If pval is all-NaN (effect-size-only): |log2FC| > 1.0
        A protein is counted if it passes in ANY comparison.
        """
        if prot_diff is None or tissue not in prot_diff:
            return -1
        df = prot_diff[tissue]
        pval_cols = [c for c in df.columns if c.startswith("pval_")]
        fc_cols = [c for c in df.columns if c.startswith("log2fc_")]
        if not pval_cols and not fc_cols:
            return -1
        sig_mask = pd.Series(False, index=df.index)
        for pc in pval_cols:
            comp_suffix = pc[len("pval_"):]
            fc_col = f"log2fc_{comp_suffix}"
            has_pval = df[pc].notna().any()
            if has_pval:
                sig_mask |= df[pc].fillna(1.0) < pval_cutoff
            elif fc_col in df.columns:
                # Effect-size-only: |log2FC| > 1.0
                sig_mask |= df[fc_col].abs() > 1.0
        return int(sig_mask.sum())

    def _count_sig_metab(tissue: str) -> int:
        """Count significant metabolites for this tissue."""
        if metab_diff is None or tissue not in metab_diff:
            return -1
        df = metab_diff[tissue]
        pval_cols = [c for c in df.columns if c.startswith("pval_") or c.startswith("fdr_")]
        if not pval_cols:
            return -1
        sig_mask = pd.Series(False, index=df.index)
        for pc in pval_cols:
            sig_mask |= df[pc].fillna(1.0) < pval_cutoff
        return int(sig_mask.sum())

    rows = []
    for tissue in sorted(set(tissues)):
        row = {"tissue": tissue}
        row["Metabolomics"] = _count_sig_metab(tissue)
        row["Proteomics"] = _count_sig_prot(tissue)
        row["Transcriptomics"] = _count_sig_tx(tissue)
        row["MFA Flux"] = _count_sig_mfa(tissue)
        rows.append(row)

    ov_df = pd.DataFrame(rows).set_index("tissue")

    # -1 is the "no data for this omics x tissue" sentinel. Render it as a
    # masked (blank) cell. Mapping it to 1 -- as this used to do -- made an
    # entirely absent omics layer read as "1 significant feature".
    display_df = ov_df.astype(float)
    display_df[display_df < 0] = np.nan

    fig, ax = plt.subplots(figsize=(max(6, len(display_df) * 1.2), 3.4))
    ax.set_facecolor("#E8E8E8")
    sns.heatmap(
        display_df.T,
        annot=True, fmt=".0f", cmap="YlOrRd",
        linewidths=1, ax=ax, cbar=False,
        mask=display_df.T.isna(),
    )
    ax.set_title(
        "Significant Features per Omics per Tissue\n"
        "metab / prot / tx: p<0.05 & |log2FC|>0.58 (effect-size-only fallback: |log2FC|>1.0)"
        "  |  MFA: pipeline is_significant\n"
        "grey = omics unavailable for that tissue; 0 = data present, none significant",
        fontweight="bold", fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()
    fig.savefig(out / "multi_omics_overview.png", dpi=300)
    plt.close(fig)

    # Save table to common/
    common_dir = cfg.COMMON_DIR
    common_dir.mkdir(parents=True, exist_ok=True)
    ov_df.to_csv(common_dir / "multi_omics_overview.csv", encoding="utf-8-sig")

    cfg.flush_print(f"  [overview] {len(ov_df)} tissues, saved overview table + heatmap")
    return ov_df


def plot_differential_summary(
    mfa_diff: dict[str, pd.DataFrame],
    metab_diff: dict[str, pd.DataFrame],
    prot_diff: dict[str, pd.DataFrame],
    tx_de: dict[str, dict[str, pd.DataFrame]],
    tissues: list[str],
    comparisons: list[str],
    output_dir: Path | None = None,
) -> pd.DataFrame:
    """Bar chart summarising differential features per omics per tissue.

    Parameters
    ----------
    mfa_diff : dict
        ``{comparison: DataFrame}`` with columns flux_log2fc, flux_pvalue.
    metab_diff : dict
        ``{tissue: DataFrame}`` with columns pval_{comp}, log2fc_{comp}.
    prot_diff : dict
        ``{tissue: DataFrame}`` with columns pval_{comp}, log2fc_{comp}.
    tx_de : dict
        ``{comparison: {tissue: DataFrame}}`` with columns log2fc, padj.
    tissues : list[str]
    comparisons : list[str]
    output_dir : Path, optional

    Returns
    -------
    pd.DataFrame
        Summary counts.
    """
    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    pval_cut = cfg.PVAL_CUTOFF
    fc_cut = cfg.LOG2FC_CUTOFF

    rows = []
    for tissue in sorted(set(tissues)):
        for comp in comparisons:
            row = {"tissue": tissue, "comparison": comp}

            # MFA
            if comp in mfa_diff:
                df = mfa_diff[comp]
                tdf = df[df["tissue"] == tissue] if "tissue" in df.columns else df
                if not tdf.empty:
                    fcol = "flux_log2fc" if "flux_log2fc" in tdf.columns else "log2FC"
                    if fcol in tdf.columns:
                        # Prefer the MFA pipeline's own verdict: it encodes the
                        # power-aware rule actually applied (effect-size-only
                        # fallback when the rank-test p floor puts p<0.05 out
                        # of reach). See _count_sig_mfa above.
                        if "is_significant" in tdf.columns:
                            mask = _truthy(tdf["is_significant"])
                        else:
                            pcol = "flux_pvalue" if "flux_pvalue" in tdf.columns else "pvalue_ttest"
                            if pcol in tdf.columns:
                                mask = (tdf[pcol].fillna(1.0) < pval_cut) & (tdf[fcol].abs() > fc_cut)
                            else:
                                mask = pd.Series(False, index=tdf.index)
                        sig = tdf[mask]
                        row["flux_up"] = int((sig[fcol] > 0).sum())
                        row["flux_down"] = int((sig[fcol] < 0).sum())

            # Metabolomics
            if tissue in metab_diff:
                df = metab_diff[tissue]
                pcol = f"pval_{comp}"
                fcol = f"log2fc_{comp}"
                if pcol in df.columns and fcol in df.columns:
                    sig = df[(df[pcol] < pval_cut) & (df[fcol].abs() > fc_cut)]
                    row["metab_up"] = int((sig[fcol] > 0).sum())
                    row["metab_down"] = int((sig[fcol] < 0).sum())

            # Proteomics
            if tissue in prot_diff:
                df = prot_diff[tissue]
                pcol = f"pval_{comp}"
                fcol = f"log2fc_{comp}"
                if pcol in df.columns and fcol in df.columns:
                    has_pval = df[pcol].notna().any()
                    if has_pval:
                        sig = df[(df[pcol].fillna(1.0) < pval_cut) & (df[fcol].abs() > fc_cut)]
                    else:
                        sig = df[df[fcol].abs() > 1.0]
                    row["prot_up"] = int((sig[fcol] > 0).sum())
                    row["prot_down"] = int((sig[fcol] < 0).sum())

            # Transcriptomics
            if comp in tx_de and tissue in tx_de[comp]:
                df = tx_de[comp][tissue]
                padj_col = "padj" if "padj" in df.columns else "pvalue"
                fc_col = "log2fc" if "log2fc" in df.columns else "log2FoldChange"
                if fc_col in df.columns:
                    has_pval = padj_col in df.columns and df[padj_col].notna().any()
                    if has_pval:
                        sig = df[(df[padj_col].fillna(1.0) < pval_cut) & (df[fc_col].abs() > fc_cut)]
                    else:
                        sig = df[df[fc_col].abs() > 1.0]
                    row["tx_up"] = int((sig[fc_col] > 0).sum())
                    row["tx_down"] = int((sig[fc_col] < 0).sum())

            rows.append(row)

    summary_df = pd.DataFrame(rows)

    if summary_df.empty:
        return summary_df

    # Save table to common/
    common_dir = cfg.COMMON_DIR
    common_dir.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(common_dir / "differential_summary.csv", index=False, encoding="utf-8-sig")

    # Plot grouped bars per comparison
    omics_colors = cfg.PLOT_COLORS.get("omics", {})
    for comp in comparisons:
        cdf = summary_df[summary_df["comparison"] == comp]
        if cdf.empty:
            continue

        pairs = [
            ("flux_up", "flux_down", "MFA Flux", omics_colors.get("mfa", "#8172B2")),
            ("metab_up", "metab_down", "Metabolomics", omics_colors.get("metabolomics", "#55A868")),
            ("prot_up", "prot_down", "Proteomics", omics_colors.get("proteomics", "#C44E52")),
            ("tx_up", "tx_down", "Transcriptomics", omics_colors.get("transcriptomics", "#4C72B0")),
        ]
        active = [p for p in pairs if p[0] in cdf.columns]
        if not active:
            continue

        n_panels = len(active)
        fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels + 2, 5), sharey=False)
        if n_panels == 1:
            axes = [axes]

        tissues_list = cdf["tissue"].values
        for ax, (up_col, dn_col, title, color) in zip(axes, active):
            ups = cdf[up_col].fillna(0).values if up_col in cdf.columns else np.zeros(len(tissues_list))
            dns = cdf[dn_col].fillna(0).values if dn_col in cdf.columns else np.zeros(len(tissues_list))
            x = np.arange(len(tissues_list))
            ax.bar(x, ups, 0.6, label="Up", color=color, alpha=0.8)
            ax.bar(x, -dns, 0.6, label="Down", color=color, alpha=0.4)
            ax.axhline(0, color="k", linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels(tissues_list, rotation=40, ha="right", fontsize=9)
            ax.set_title(title, fontweight="bold")
            ax.set_ylabel("# Up / # Down")
            ax.legend(fontsize=8)
            ax.grid(axis="y", alpha=0.3)

        label = cfg.COMP_LABELS.get(comp, comp)
        plt.suptitle(f"Differential Features: {label}", fontweight="bold", fontsize=13)
        plt.tight_layout()
        fig.savefig(out / f"differential_summary_{comp}.png", dpi=300)
        plt.close(fig)

    cfg.flush_print(f"  [diff_summary] Saved {len(comparisons)} summary plots")
    return summary_df


def plot_tissue_dashboard(
    mfa_diff: dict[str, pd.DataFrame],
    metab_diff: dict[str, pd.DataFrame],
    prot_diff: dict[str, pd.DataFrame],
    tx_gsea: dict[str, dict[str, pd.DataFrame]],
    tissues: list[str],
    comparison: str,
    comparison_map: dict[str, str] | None = None,
    output_dir: Path | None = None,
) -> pd.DataFrame:
    """Tissue integration dashboard for a single comparison.

    Adapted from multiomics_integration_analysis.py lines 530-625.
    """
    out = output_dir or cfg.PLOTS_DIR
    out.mkdir(parents=True, exist_ok=True)

    comp_map = comparison_map or {}
    omics_colors = cfg.PLOT_COLORS.get("omics", {})

    dashboard_rows = []
    for tissue in sorted(set(tissues)):
        row = {"tissue": tissue}

        # Metabolomics
        if tissue in metab_diff:
            df = metab_diff[tissue]
            fc_col = f"log2fc_{comparison}"
            p_col = f"pval_{comparison}"
            if fc_col in df.columns and p_col in df.columns:
                sig = df[(df[p_col] < cfg.PVAL_CUTOFF) & (df[fc_col].abs() > cfg.LOG2FC_CUTOFF)]
                row["metabo_up"] = int((sig[fc_col] > 0).sum())
                row["metabo_down"] = int((sig[fc_col] < 0).sum())

        # Proteomics
        if tissue in prot_diff:
            df = prot_diff[tissue]
            fc_col = f"log2fc_{comparison}"
            p_col = f"pval_{comparison}"
            if fc_col in df.columns and p_col in df.columns:
                has_pval = df[p_col].notna().any()
                if has_pval:
                    sig = df[(df[p_col].fillna(1.0) < cfg.PVAL_CUTOFF) & (df[fc_col].abs() > cfg.LOG2FC_CUTOFF)]
                else:
                    sig = df[df[fc_col].abs() > 1.0]
                row["protein_up"] = int((sig[fc_col] > 0).sum())
                row["protein_down"] = int((sig[fc_col] < 0).sum())

        # Transcriptomics GSEA
        tx_comp = comp_map.get(comparison, comparison)
        if tissue in tx_gsea and tx_comp in tx_gsea[tissue]:
            gdf = tx_gsea[tissue][tx_comp]
            if "FDR q-val" in gdf.columns and "NES" in gdf.columns:
                sig = gdf[gdf["FDR q-val"] < 0.25]
                row["gsea_up"] = int((sig["NES"] > 0).sum())
                row["gsea_down"] = int((sig["NES"] < 0).sum())

        # MFA differential fluxes — use FDR + Cohen's d gate criteria
        if comparison in mfa_diff:
            fdf = mfa_diff[comparison]
            tdf = fdf[fdf["tissue"] == tissue] if "tissue" in fdf.columns else fdf
            if not tdf.empty:
                fdr_col = "fdr_mannwhitney" if "fdr_mannwhitney" in tdf.columns else None
                cd_col = "cohen_d" if "cohen_d" in tdf.columns else None
                fcol = "flux_log2fc" if "flux_log2fc" in tdf.columns else "log2FC"
                if fdr_col and cd_col and fcol in tdf.columns:
                    sig = tdf[
                        (tdf[fdr_col].fillna(1.0) < 0.05) &
                        (tdf[cd_col].abs().fillna(0) > 0.8)
                    ]
                else:
                    pcol = "flux_pvalue" if "flux_pvalue" in tdf.columns else "pvalue_ttest"
                    if pcol in tdf.columns and fcol in tdf.columns:
                        sig = tdf[
                            (tdf[pcol].fillna(1.0) < cfg.PVAL_CUTOFF) &
                            (tdf[fcol].abs().fillna(0) > cfg.LOG2FC_CUTOFF)
                        ]
                    else:
                        sig = pd.DataFrame()
                row["flux_up"] = int((sig[fcol] > 0).sum()) if not sig.empty else 0
                row["flux_down"] = int((sig[fcol] < 0).sum()) if not sig.empty else 0

        dashboard_rows.append(row)

    dash_df = pd.DataFrame(dashboard_rows)

    # Save table to common/
    common_dir = cfg.COMMON_DIR
    common_dir.mkdir(parents=True, exist_ok=True)
    dash_df.to_csv(common_dir / f"tissue_dashboard_{comparison}.csv", index=False, encoding="utf-8-sig")

    # Stacked bar chart
    tissues_list = dash_df["tissue"].values
    pairs = [
        ("metabo_up", "metabo_down", "Metabolomics", omics_colors.get("metabolomics", "#55A868")),
        ("protein_up", "protein_down", "Proteomics", omics_colors.get("proteomics", "#C44E52")),
        ("gsea_up", "gsea_down", "Transcriptomics GSEA", omics_colors.get("transcriptomics", "#4C72B0")),
        ("flux_up", "flux_down", "MFA Flux", omics_colors.get("mfa", "#8172B2")),
    ]
    active = [p for p in pairs if p[0] in dash_df.columns or p[1] in dash_df.columns]
    if not active:
        active = pairs[:3]

    n_panels = len(active)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels + 2, 5), sharey=False)
    if n_panels == 1:
        axes = [axes]

    for ax, (up_col, dn_col, title, color) in zip(axes, active):
        ups = dash_df[up_col].fillna(0).values if up_col in dash_df.columns else np.zeros(len(tissues_list))
        dns = dash_df[dn_col].fillna(0).values if dn_col in dash_df.columns else np.zeros(len(tissues_list))
        x = np.arange(len(tissues_list))
        ax.bar(x, ups, 0.6, label="Up", color=color, alpha=0.8)
        ax.bar(x, -dns, 0.6, label="Down", color=color, alpha=0.4)
        ax.axhline(0, color="k", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(tissues_list, rotation=40, ha="right", fontsize=9)
        ax.set_title(title, fontweight="bold")
        ax.set_ylabel("# Up / # Down")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    label = cfg.COMP_LABELS.get(comparison, comparison)
    plt.suptitle(f"Multi-Omics Tissue Dashboard: {label}", fontweight="bold", fontsize=13)
    plt.tight_layout()
    fig.savefig(out / f"tissue_dashboard_{comparison}.png", dpi=300)
    plt.close(fig)

    cfg.flush_print(f"  [dashboard] Saved tissue_dashboard_{comparison}.csv + .png")
    return dash_df