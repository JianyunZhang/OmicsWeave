#!/usr/bin/env python3
"""
MFA differential flux volcano plots with reaction labels.

Per-comparison volcanoes are saved at the top level of output_dir.
Per-tissue volcanoes are saved under output_dir/<tissue>/.

Two y-axis modes
----------------
``pvalue``       classic volcano, y = -log10(p).
``effect_size``  y = |Hedges' g|, used when the p-value carries no evidence.

The second mode exists because MFA runs at n=2-4 per group. There the
Mann-Whitney two-sided p-value has a floor of ``2/C(n1+n2, n1)`` and takes
only a handful of discrete values (n=3 vs 3 -> min 0.1, so -log10(p) <= 1.0
while the p<0.05 line sits at 1.301). Plotting that as a volcano produces a
few horizontal bands of points underneath a threshold line nothing can reach.
The MFA pipeline itself already switches to an effect-size rule in this
regime (``evidence_level = effect_size_only``, ``significance_rule =
|Hedges' g|>=0.8 and |log2FC|>=0.58``); these plots now follow it.

Axes are clipped to robust percentile limits. Out-of-range points are drawn
at the boundary with a triangle marker rather than dropped -- net fluxes that
cross zero produce |log2FC| outliers (a flux going +5 -> -5 is a direction
reversal, not a fold change), and variance collapse in the fit can push
|Hedges' g| to absurd magnitudes. Both are visible at the axis edge.

Author: Claude Code
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator

try:
    from adjustText import adjust_text
    _HAS_ADJUST_TEXT = True
except ImportError:
    _HAS_ADJUST_TEXT = False


# Effect-size threshold used by the MFA significance rule and by
# FluxOmicsScorer's fallback gate (scoring.fallback_mode.mfa_gate).
_DEFAULT_EFFECT_CUTOFF = 0.8


def _truthy(series: pd.Series) -> pd.Series:
    """Parse a significance flag column that may be bool, int or string."""
    return series.map(lambda v: str(v).strip().lower() in {"true", "1", "yes"})


def _pick_effect_col(df: pd.DataFrame) -> str | None:
    """Prefer Hedges' g (small-sample corrected) over Cohen's d."""
    for col in ("hedges_g", "cohen_d"):
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col
    return None


def _resolve_y_mode(df: pd.DataFrame, pval_col: str, pval_cutoff: float) -> str:
    """Decide whether the p-value axis carries evidence for this table.

    Returns "pvalue", "effect_size" or "insufficient".

    Two independent triggers move off the p-value axis, either is sufficient:
    1. The MFA pipeline declared the regime via ``evidence_level`` -- every row
       is ``effect_size_only`` and/or ``insufficient``.
    2. Empirically no p-value clears the threshold, i.e. the significance line
       would sit above every point. That is the rank-test floor showing up.

    ``insufficient`` is returned when neither axis carries evidence: n<2 in one
    group leaves Hedges' g undefined (the MFA pipeline emits all-NaN effect
    sizes and forces ``is_significant`` False), while the p-value is still at
    the rank floor. Nothing can be inferred, so the plot degrades to a
    descriptive log2FC strip rather than a volcano nothing could ever populate.
    """
    has_effect = _pick_effect_col(df) is not None

    declared_low_power = False
    declared_insufficient = False
    if "evidence_level" in df.columns:
        levels = set(df["evidence_level"].dropna().astype(str).str.strip())
        if levels:
            declared_low_power = levels <= {"effect_size_only", "insufficient"}
            declared_insufficient = levels <= {"insufficient"}

    pval_at_floor = False
    if pval_col in df.columns:
        pv = pd.to_numeric(df[pval_col], errors="coerce").dropna()
        pval_at_floor = pv.empty or float(pv.min()) >= pval_cutoff

    if has_effect and (declared_low_power or pval_at_floor):
        return "effect_size"
    if not has_effect and (declared_insufficient or pval_at_floor):
        return "insufficient"
    return "pvalue"


def _robust_limits(values: pd.Series, lo_q: float = 0.01, hi_q: float = 0.99,
                   pad: float = 0.08) -> tuple[float, float]:
    """Percentile-based axis limits, padded. Falls back to min/max."""
    v = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if v.empty:
        return (-1.0, 1.0)
    lo, hi = float(v.quantile(lo_q)), float(v.quantile(hi_q))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(v.min()), float(v.max())
    if hi <= lo:
        lo, hi = lo - 1.0, hi + 1.0
    span = hi - lo
    return (lo - pad * span, hi + pad * span)


def _symmetric_x_limits(x: pd.Series, fc_cutoff: float,
                        keep: pd.Series | None = None,
                        q: float = 0.95, pad: float = 0.10) -> tuple[float, float]:
    """Symmetric log2FC limits that never clip a significant point.

    A plain percentile is too weak here: a per-tissue table holds ~70-140
    fluxes, so quantile(0.01) lands on the extreme value itself and the axis
    stretches to the outlier. Net fluxes crossing zero routinely produce
    |log2FC| of 30+ while the bulk sits under 3, which compresses all real
    structure into a sliver around 0.

    The limit is the widest of: |log2FC| at quantile ``q``, the largest
    |log2FC| among points marked significant, and twice the fold-change
    cutoff -- so the significance bands stay visible and no called point is
    pushed to the boundary.
    """
    a = x.abs().replace([np.inf, -np.inf], np.nan).dropna()
    if a.empty:
        return (-1.0, 1.0)
    lim = float(a.quantile(q))
    if keep is not None and keep.any():
        kept = x[keep].abs().replace([np.inf, -np.inf], np.nan).dropna()
        if not kept.empty:
            lim = max(lim, float(kept.max()))
    lim = max(lim, 2.0 * fc_cutoff)
    if not np.isfinite(lim) or lim <= 0:
        lim = 1.0
    lim *= (1.0 + pad)
    return (-lim, lim)


def _plain_number(value: float, _pos=None) -> str:
    """Tick label without mathtext.

    The log-scale default renders 10^-1 through mathtext, whose minus sign is
    missing from the Chinese fonts this project sets as font.sans-serif
    (SimHei / Microsoft YaHei) -- it came out as a currency glyph. Plain
    decimals sidestep the issue and read better for effect sizes.
    """
    if value <= 0:
        return ""
    return f"{value:,.0f}" if value >= 1 else f"{value:g}"


def plot_mfa_volcanoes(
    mfa_diff: dict[str, pd.DataFrame],
    output_dir: Path | None = None,
    pval_cutoff: float = 0.05,
    fc_cutoff: float = 0.58,
    tissues: list[str] | None = None,
    max_labels: int = 25,
    effect_cutoff: float = _DEFAULT_EFFECT_CUTOFF,
) -> list[Path]:
    """Generate volcano plots for each MFA differential flux comparison.

    Generates two sets:
    - Per-comparison (all tissues merged) at ``output_dir/mfa_volcano_<comp>.png``
    - Per-tissue per-comparison at ``output_dir/<tissue>/mfa_volcano_<comp>.png``

    Parameters
    ----------
    mfa_diff : dict
        ``{comparison: DataFrame}``. Uses ``flux_log2fc``/``log2FC`` for x and,
        depending on the resolved mode, ``flux_pvalue`` or ``hedges_g``.
        Optional: ``is_significant``, ``evidence_level``, ``flux_id``,
        ``reaction``, ``tissue``.
    output_dir : Path, optional
    pval_cutoff, fc_cutoff : float
        Thresholds for the classic p-value mode.
    tissues : list[str], optional
        When provided (and a ``tissue`` column exists), also emit per-tissue plots.
    max_labels : int
        Maximum number of reaction labels to annotate.
    effect_cutoff : float
        |Hedges' g| threshold used in effect-size mode (default 0.8, matching
        the MFA significance rule).

    Returns
    -------
    list[Path]
        Paths to saved PNG files.
    """
    out = output_dir or Path(".") / "02_result_plot"
    out.mkdir(parents=True, exist_ok=True)

    saved: list[Path] = []

    for comp, df in mfa_diff.items():
        fc_col = "flux_log2fc" if "flux_log2fc" in df.columns else "log2FC"
        pval_col = "flux_pvalue" if "flux_pvalue" in df.columns else "pvalue_ttest"

        if fc_col not in df.columns:
            print(f"[SKIP] Volcano for {comp}: missing {fc_col}")
            continue

        mode = _resolve_y_mode(df, pval_col, pval_cutoff)
        if mode == "effect_size":
            y_col = _pick_effect_col(df)
        elif mode == "insufficient":
            y_col = fc_col  # y is unused; kept non-null so dropna() works
        else:
            y_col = pval_col
        if y_col is None or y_col not in df.columns:
            print(f"[SKIP] Volcano for {comp}: no usable y-axis column")
            continue
        if mode == "effect_size":
            print(f"  [volcano] {comp}: effect-size mode "
                  f"(y=|{y_col}|) -- p-values are at the rank-test floor")
        elif mode == "insufficient":
            print(f"  [volcano] {comp}: insufficient mode -- n<2 in one group, "
                  "effect size undefined; descriptive log2FC strip only")

        label_col = _pick_label_col(df)

        # ---- Per-comparison volcano (all tissues) ----
        df_plot = df.dropna(subset=[fc_col, y_col]).copy()
        if not df_plot.empty:
            _volcano_single(
                df_plot, fc_col, y_col, mode,
                pval_cutoff, fc_cutoff, effect_cutoff,
                label_col=label_col, max_labels=max_labels,
                title=f"MFA Differential Flux: {comp}",
                save_path=out / f"mfa_volcano_{comp}.png",
            )
            saved.append(out / f"mfa_volcano_{comp}.png")

        # ---- Per-tissue volcanoes ----
        if tissues and "tissue" in df.columns:
            for tissue in tissues:
                tdf = df[df["tissue"] == tissue].dropna(subset=[fc_col, y_col]).copy()
                if tdf.empty:
                    continue
                tissue_dir = out / tissue
                tissue_dir.mkdir(parents=True, exist_ok=True)
                _volcano_single(
                    tdf, fc_col, y_col, mode,
                    pval_cutoff, fc_cutoff, effect_cutoff,
                    label_col=label_col, max_labels=max_labels,
                    title=f"MFA Flux: {tissue} - {comp}",
                    save_path=tissue_dir / f"mfa_volcano_{comp}.png",
                )
                saved.append(tissue_dir / f"mfa_volcano_{comp}.png")

    n_per_comp = len(mfa_diff)
    n_per_tissue = len(saved) - n_per_comp if len(saved) > n_per_comp else 0
    msg = f"  [volcano] Saved {n_per_comp} per-comparison volcano plots"
    if n_per_tissue > 0:
        msg += f" + {n_per_tissue} per-tissue volcano plots"
    print(msg)
    return saved


def _pick_label_col(df: pd.DataFrame) -> str | None:
    """Choose the best column for point labels.

    Priority: reaction (resolved to short display name) > flux_id.
    """
    if "reaction" in df.columns:
        return "reaction"
    if "flux_id" in df.columns:
        return "flux_id"
    return None


def _volcano_single(
    df: pd.DataFrame,
    fc_col: str,
    y_col: str,
    mode: str,
    pval_cutoff: float,
    fc_cutoff: float,
    effect_cutoff: float,
    *,
    label_col: str | None = None,
    max_labels: int = 25,
    title: str,
    save_path: Path,
) -> None:
    """Render and save a single volcano scatter plot with reaction labels."""
    df = df.copy()

    # ---- y axis ----------------------------------------------------------
    if mode == "insufficient":
        # No inference is possible; y carries nothing. Spread points vertically
        # by rank only, so overlapping log2FC values stay individually visible.
        x_raw = pd.to_numeric(df[fc_col], errors="coerce")
        df["_y"] = x_raw.rank(method="first", na_option="keep") / max(len(df), 1)
        y_label = "(no inference possible - points spread by rank)"
        y_thresh = None
        y_log = False
    elif mode == "effect_size":
        df["_y"] = pd.to_numeric(df[y_col], errors="coerce").abs()
        y_label = f"|{'Hedges g' if y_col == 'hedges_g' else 'Cohen d'}| (effect size)"
        y_thresh = effect_cutoff
        y_log = True
    else:
        pv = pd.to_numeric(df[y_col], errors="coerce").clip(lower=1e-50)
        df["_y"] = -np.log10(pv)
        y_label = "-log10(p-value)"
        y_thresh = -np.log10(pval_cutoff)
        y_log = False

    df["_x"] = pd.to_numeric(df[fc_col], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["_x", "_y"])
    if df.empty:
        return

    # ---- significance ----------------------------------------------------
    # Prefer the MFA pipeline's own verdict; it encodes the power-aware rule
    # actually applied. Re-deriving it here is only a fallback for legacy
    # results without the column.
    if mode == "insufficient":
        sig_mask = pd.Series(False, index=df.index)
        sig_note = "none callable (n<2 in one group)"
    elif "is_significant" in df.columns:
        sig_mask = _truthy(df["is_significant"])
        sig_note = "pipeline is_significant"
    elif mode == "effect_size":
        sig_mask = (df["_y"] >= effect_cutoff) & (df["_x"].abs() >= fc_cutoff)
        sig_note = f"|effect|>={effect_cutoff} & |log2FC|>={fc_cutoff}"
    else:
        pv = pd.to_numeric(df[y_col], errors="coerce")
        sig_mask = (pv.fillna(1.0) < pval_cutoff) & (df["_x"].abs() > fc_cutoff)
        sig_note = f"p<{pval_cutoff} & |log2FC|>{fc_cutoff}"
    n_sig = int(sig_mask.sum())

    # ---- robust axis limits; out-of-range points pinned to the edge -------
    x_lo, x_hi = _symmetric_x_limits(df["_x"], fc_cutoff, keep=sig_mask)
    if y_log:
        pos = df.loc[df["_y"] > 0, "_y"]
        y_lo = max(float(pos.quantile(0.01)), 1e-3) if not pos.empty else 1e-3
        y_lo = min(y_lo, effect_cutoff / 2)
        y_hi = float(df["_y"].quantile(0.99)) if len(df) > 10 else float(df["_y"].max())
        y_hi = max(y_hi, effect_cutoff * 2)
        y_hi *= 1.5
    elif mode == "insufficient":
        y_lo, y_hi = -0.05, 1.05
    else:
        y_lo, y_hi = _robust_limits(df["_y"], lo_q=0.0, hi_q=1.0)
        y_lo = min(y_lo, 0.0)
        y_hi = max(y_hi, y_thresh * 1.05)

    df["_xp"] = df["_x"].clip(x_lo, x_hi)
    df["_yp"] = df["_y"].clip(y_lo, y_hi)
    clipped = (df["_xp"] != df["_x"]) | (df["_yp"] != df["_y"])
    n_clipped = int(clipped.sum())

    fig, ax = plt.subplots(figsize=(14, 9))
    if y_log:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=15))
        ax.yaxis.set_minor_locator(
            LogLocator(base=10.0, subs=tuple(np.arange(2, 10) * 0.1), numticks=15))
        ax.yaxis.set_major_formatter(FuncFormatter(_plain_number))
        ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _p: ""))

    groups = [
        (~sig_mask, "#999999", 20, 0.5, "NS"),
        (sig_mask & (df["_x"] > 0), "#C44E52", 30, 0.75, "Up"),
        (sig_mask & (df["_x"] < 0), "#4C72B0", 30, 0.75, "Down"),
    ]
    for mask, color, size, alpha, label in groups:
        inr = mask & ~clipped
        ax.scatter(df.loc[inr, "_xp"], df.loc[inr, "_yp"],
                   c=color, s=size, alpha=alpha, label=label)
        oor = mask & clipped
        if oor.any():
            ax.scatter(df.loc[oor, "_xp"], df.loc[oor, "_yp"],
                       c=color, s=size + 25, alpha=alpha, marker="^",
                       edgecolors="black", linewidths=0.4)

    if y_thresh is not None:
        ax.axhline(y_thresh, color="grey", linestyle="--", alpha=0.6)
    if mode == "insufficient":
        ax.set_yticks([])
    ax.axvline(-fc_cutoff, color="grey", linestyle="--", alpha=0.5)
    ax.axvline(fc_cutoff, color="grey", linestyle="--", alpha=0.5)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    # ---- annotate top reactions -----------------------------------------
    texts = []
    if label_col and n_sig > 0:
        sig_df = df[sig_mask].copy()
        if mode == "effect_size":
            # Largest effect first.
            sig_df = sig_df.nlargest(max_labels, "_y")
        else:
            # Smallest p-value relative to largest |FC|.
            pv = pd.to_numeric(sig_df[y_col], errors="coerce").fillna(1.0)
            sig_df["_rank"] = pv / sig_df["_x"].abs().clip(lower=0.01)
            sig_df = sig_df.nsmallest(max_labels, "_rank")

        for _, row in sig_df.iterrows():
            raw = str(row[label_col]) if pd.notna(row[label_col]) else ""
            if not raw:
                continue
            txt = ax.text(
                row["_xp"], row["_yp"], raw,
                fontsize=6.5, fontweight="bold", color="black",
                alpha=0.85, ha="left", va="bottom",
            )
            texts.append(txt)

        if texts and _HAS_ADJUST_TEXT:
            adjust_text(
                texts, ax=ax,
                expand_points=(1.5, 1.8),
                force_text=(0.3, 0.5),
                force_points=(0.5, 0.8),
                lim=100,
            )

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel(y_label)

    subtitle = f"n_sig={n_sig}  |  significance: {sig_note}"
    if mode == "insufficient":
        subtitle += (
            "\nNO INFERENCE: one group has n<2, so Hedges' g is undefined and "
            "the pipeline forces is_significant=False."
            "\nlog2FC is shown as a descriptive distribution only - it cannot "
            "support a difference claim."
        )
    elif mode == "effect_size":
        subtitle += (
            "\np-value axis unusable at this n (Mann-Whitney floor 2/C(n1+n2,n1) "
            f">= {pval_cutoff}); y = effect size, log scale"
        )
    if n_clipped:
        subtitle += (
            f"\n{n_clipped} point(s) outside robust axis limits, "
            "drawn as triangles at the boundary"
        )
    ax.set_title(f"{title}\n{subtitle}", fontweight="bold", fontsize=10)
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
