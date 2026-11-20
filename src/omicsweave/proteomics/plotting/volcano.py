"""
Volcano plot module for proteomics.

Generates volcano plots from per-comparison DE results.

Plotting conventions (see README):

- **The y-axis is the raw p-value, the colour is the BH-adjusted p-value.**
  Plotting ``-log10(padj)`` directly collapses the cloud into horizontal
  stripes, because BH is a step function: its monotonicity enforcement ties
  large blocks of proteins to one adjusted value (here up to 8474 proteins —
  100% of everything tested — at a single padj). Those stripes are an artefact
  of the adjustment, not structure in the data. The FDR cut is preserved
  exactly by drawing the threshold line at the largest raw p whose padj still
  passes, which is well defined because BH is monotone in raw p.
- **Untestable proteins are counted, never silently dropped.** Features with
  fewer than ``min_observations_per_group`` samples in a group get no p-value
  and would vanish from the figure without trace — 500-1600 per panel in this
  project, and they carry the largest fold changes (down to log2FC = -9.8).
  The excluded count goes in the subtitle.
- **Draw order is explicit**: NS first, then Down/Up on top, so significant
  points are never buried under the grey cloud.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from omicsweave.visualization.palettes import NATURE_EDITORIAL

from . import _config as cfg

_PALETTE = {"Up": NATURE_EDITORIAL.up, "Down": NATURE_EDITORIAL.down, "NS": NATURE_EDITORIAL.faint}
_DRAW_ORDER = ("NS", "Down", "Up")
_MAX_LABELS = 10


def classify_significance(
    log2_fold_change: pd.Series | np.ndarray,
    adjusted_pvalue: pd.Series | np.ndarray,
    log2fc_cutoff: float,
    adjusted_pvalue_cutoff: float,
) -> np.ndarray:
    """Classify volcano points using the same adjusted-p rule as result CSVs."""

    fold_change = np.asarray(log2_fold_change, dtype=float)
    padj = np.asarray(adjusted_pvalue, dtype=float)
    return np.select(
        [
            (fold_change > log2fc_cutoff) & (padj < adjusted_pvalue_cutoff),
            (fold_change < -log2fc_cutoff) & (padj < adjusted_pvalue_cutoff),
        ],
        ["Up", "Down"],
        default="NS",
    )


def fdr_cut_on_raw_p(
    pvalue: pd.Series,
    adjusted_pvalue: pd.Series,
    adjusted_pvalue_cutoff: float,
) -> float | None:
    """Raw-p position of the FDR threshold, or None if nothing passes.

    BH is monotone in the raw p-value, so the largest raw p among features
    that clear ``padj < cutoff`` is exactly the boundary: everything above it
    on a -log10(raw p) axis passes FDR, everything below does not.
    """
    passing = adjusted_pvalue.notna() & (adjusted_pvalue < adjusted_pvalue_cutoff)
    passing &= pvalue.notna()
    if not passing.any():
        return None
    return float(pvalue[passing].max())


def _render(
    frame: pd.DataFrame,
    title: str,
    output_path: str,
    *,
    pval_cutoff: float,
    log2fc_cutoff: float,
) -> None:
    """Draw one volcano from a frame with Protein/lfc/pvalue/padj columns."""
    frame = frame.copy()
    frame["Significance"] = classify_significance(
        frame["lfc"], frame["padj"], log2fc_cutoff, pval_cutoff,
    )

    # The y-axis is the raw p-value; fall back to padj only for legacy CSVs
    # that predate the pvalue column.
    on_raw_p = frame["pvalue"].notna().any()
    y_source = "pvalue" if on_raw_p else "padj"
    frame["y"] = -np.log10(frame[y_source].clip(lower=1e-300))

    n_total = len(frame)
    testable = frame["y"].notna() & frame["lfc"].notna()
    n_excluded = int((~testable).sum())
    frame = frame[testable]

    fig, ax = plt.subplots(figsize=(8, 6))

    for level in _DRAW_ORDER:
        subset = frame[frame["Significance"] == level]
        if subset.empty:
            continue
        is_ns = level == "NS"
        ax.scatter(
            subset["lfc"], subset["y"],
            s=12 if is_ns else 18,
            c=_PALETTE[level],
            alpha=0.35 if is_ns else 0.85,
            linewidths=0,
            label=f"{level} ({len(subset)})",
            rasterized=is_ns,
            zorder=1 if is_ns else 3,
        )

    ax.axvline(x=log2fc_cutoff, color="k", linestyle=":", alpha=0.5, zorder=2)
    ax.axvline(x=-log2fc_cutoff, color="k", linestyle=":", alpha=0.5, zorder=2)

    if on_raw_p:
        threshold = fdr_cut_on_raw_p(frame["pvalue"], frame["padj"], pval_cutoff)
        fdr_note = ""
        if threshold is not None:
            ax.axhline(-np.log10(threshold), color="k", linestyle=":",
                       alpha=0.6, zorder=2)
            ax.text(
                ax.get_xlim()[1], -np.log10(threshold),
                f"BH FDR {pval_cutoff:g}", fontsize=7, ha="right", va="bottom",
                zorder=4,
                bbox={"facecolor": "white", "edgecolor": "none",
                      "pad": 1.5, "alpha": 0.85},
            )
        else:
            # No line rather than a line nothing can reach: say so explicitly,
            # so an empty panel reads as "unreachable", not "forgotten".
            fdr_note = f"no protein reaches BH FDR < {pval_cutoff:g}"
        ax.set_ylabel("-Log10 raw P-value")
    else:
        ax.axhline(y=-np.log10(pval_cutoff), color="k", linestyle=":",
                   alpha=0.6, zorder=2)
        fdr_note = "y-axis is adjusted P (no raw P in source)"
        ax.set_ylabel("-Log10 adjusted P-value")

    _label_top_hits(ax, frame)

    notes = [f"colour = BH FDR < {pval_cutoff:g} and |log2FC| > {log2fc_cutoff:g}"]
    if n_excluded:
        notes.append(f"{n_excluded}/{n_total} not testable (excluded)")
    if fdr_note:
        notes.append(fdr_note)
    ax.set_title(f"{title}\n" + "\n".join(notes), fontsize=9)
    ax.set_xlabel("Log2 Fold Change")
    ax.legend(
        title="Significance", loc="upper left", bbox_to_anchor=(1.01, 1.0),
        frameon=False, fontsize=8, title_fontsize=9, markerscale=1.8,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _label_top_hits(ax: plt.Axes, frame: pd.DataFrame) -> None:
    """Annotate the strongest hits with a deterministic tie-break.

    Sorting on padj alone is not deterministic here: BH ties whole blocks of
    proteins to one value, so a head(n) cut can land inside a tie group and
    which names survive is then decided by sort order, not by the data.
    Breaking ties on |log2FC| descending makes the label set reproducible.
    """
    hits = frame[frame["Significance"] != "NS"]
    if hits.empty:
        return
    hits = hits.assign(_abs_lfc=hits["lfc"].abs()).sort_values(
        ["padj", "_abs_lfc"], ascending=[True, False], kind="mergesort",
    ).head(_MAX_LABELS)

    for _, row in hits.iterrows():
        ax.annotate(
            row["Protein"],
            xy=(row["lfc"], row["y"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
            zorder=4,
        )


def _frame_from_per_comparison(df: pd.DataFrame) -> pd.DataFrame:
    if "padj" not in df.columns:
        raise ValueError("Volcano input requires a 'padj' column")
    return pd.DataFrame({
        "Protein": df["Protein"],
        "lfc": df["log2FoldChange"],
        "pvalue": df["pvalue"] if "pvalue" in df.columns else np.nan,
        "padj": df["padj"],
    })


def _plot_single(
    df: pd.DataFrame,
    tissue: str,
    comp_name: str,
    output_path: str,
    *,
    pval_cutoff: float | None = None,
    log2fc_cutoff: float | None = None,
    method: str | None = None,
) -> None:
    """Generate a single volcano plot from per-comparison DE CSV."""
    pv_cut = cfg.PVAL_CUTOFF if pval_cutoff is None else pval_cutoff
    fc_cut = cfg.LOG2FC_CUTOFF if log2fc_cutoff is None else log2fc_cutoff
    method_label = method or cfg.DP_METHOD

    label = cfg.COMP_LABELS.get(comp_name, comp_name)
    tissue_label = cfg.TISSUE_LABELS.get(tissue, tissue)
    _render(
        _frame_from_per_comparison(df),
        f"{tissue_label} - {label} ({method_label})",
        output_path,
        pval_cutoff=pv_cut,
        log2fc_cutoff=fc_cut,
    )


def _plot_single_combined(
    results_df: pd.DataFrame,
    comp_name: str,
    output_path: str,
    *,
    pval_cutoff: float | None = None,
    log2fc_cutoff: float | None = None,
    method: str | None = None,
) -> None:
    """Generate a single volcano plot from combined DE CSV columns."""
    pv_cut = cfg.PVAL_CUTOFF if pval_cutoff is None else pval_cutoff
    fc_cut = cfg.LOG2FC_CUTOFF if log2fc_cutoff is None else log2fc_cutoff
    method_label = method or cfg.DP_METHOD

    fc_col = f"Log2FC_{comp_name}"
    padj_col = f"Padj_{comp_name}"
    pval_col = f"Pval_{comp_name}"

    if fc_col not in results_df.columns or padj_col not in results_df.columns:
        return

    frame = pd.DataFrame({
        "Protein": results_df["Protein"],
        "lfc": results_df[fc_col],
        "pvalue": results_df[pval_col] if pval_col in results_df.columns else np.nan,
        "padj": results_df[padj_col],
    })

    parts = comp_name.split("_vs_")
    title = (
        f"{parts[0]} vs {parts[1]} ({method_label})"
        if len(parts) == 2
        else f"{comp_name} ({method_label})"
    )
    _render(
        frame, title, output_path,
        pval_cutoff=pv_cut, log2fc_cutoff=fc_cut,
    )


def generate_all() -> None:
    """Generate volcano plots for all tissues and comparisons."""
    cfg.flush_print("Generating volcano plots...")
    for tissue in cfg.TISSUES:
        tissue_dir = os.path.join(cfg.PLOT_BASE, tissue)
        os.makedirs(tissue_dir, exist_ok=True)

        for _, _, comp_name in cfg.COMPARISONS:
            df = cfg.load_dp_results(tissue, comp_name)
            if df is None:
                continue
            out_path = os.path.join(tissue_dir, f"volcano_{comp_name}.png")
            _plot_single(df, tissue, comp_name, out_path, method=cfg.DP_METHOD)

    cfg.flush_print("  Volcano plots done.")


def plot_single_from_combined(
    results_df: pd.DataFrame,
    comp_name: str,
    output_path: str,
    *,
    pval_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.58,
    method: str | None = None,
) -> None:
    """Public API: plot from combined DE DataFrame (used by pipeline)."""
    _plot_single_combined(
        results_df, comp_name, output_path,
        pval_cutoff=pval_cutoff, log2fc_cutoff=log2fc_cutoff, method=method,
    )
