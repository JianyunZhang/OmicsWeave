"""Proteomics ssGSEA and comparison-level pre-ranked GSEA."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd


DEFAULT_GENE_SETS: dict[str, str] = {
    "kegg": "KEGG_2019_Mouse",
    "go": "GO_Biological_Process_2025",
    "reactome": "Reactome_Pathways_2024",
    "wikipathways": "WikiPathways_2024_Mouse",
}

SSGSEA_REPORT_NAMES = (
    "gseapy.gene_set.ssgsea.report.csv",
    "gsa_report.csv",  # historical compatibility
)
PRERANK_REPORT_NAMES = (
    "gseapy.gene_set.prerank.report.csv",
    "gseapy.prerank.gene_sets.report.csv",
)


def _lib_to_folder_key(lib_name: str) -> str:
    for key, full_name in DEFAULT_GENE_SETS.items():
        if lib_name == full_name:
            return key
    path = Path(lib_name)
    candidate = path.stem if path.suffix.lower() == ".gmt" else lib_name
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", candidate).strip("._-").lower()
    return safe or "custom_gene_sets"


def _normalize_libraries(
    gene_sets: str | Sequence[str] | None,
) -> list[tuple[str, str]]:
    if gene_sets is None:
        requested = ["kegg"]
    elif isinstance(gene_sets, str):
        requested = [gene_sets]
    else:
        requested = [str(value) for value in gene_sets]
    libraries = [DEFAULT_GENE_SETS.get(value, value) for value in requested]
    return [(_lib_to_folder_key(library), library) for library in libraries]


def _find_report(directory: str | Path, names: Sequence[str], pattern: str) -> Path | None:
    folder = Path(directory)
    for name in names:
        candidate = folder / name
        if candidate.exists():
            return candidate
    matches = sorted(folder.glob(pattern)) if folder.exists() else []
    return matches[0] if matches else None


def find_ssgsea_report(directory: str | Path) -> Path | None:
    """Locate the actual report filename emitted by supported gseapy versions."""

    return _find_report(directory, SSGSEA_REPORT_NAMES, "*ssgsea*.report.csv")


def find_prerank_report(directory: str | Path) -> Path | None:
    return _find_report(directory, PRERANK_REPORT_NAMES, "*prerank*.report.csv")


def run_ssgsea(
    df: pd.DataFrame,
    output_dir: str,
    *,
    gene_sets: str | Sequence[str] | None = None,
    min_size: int = 10,
    max_size: int = 500,
    permutation_num: int = 0,
    seed: int = 42,
    threads: int = 1,
    force: bool = False,
    flush: Callable[[str], None] = print,
) -> pd.DataFrame:
    """Run sample-level ssGSEA and persist an explicit status table.

    ssGSEA scores are descriptive sample-level pathway scores. Comparison-level
    pathway inference is handled separately by :func:`run_prerank_gsea`.
    """

    import gseapy as gp

    rows: list[dict[str, Any]] = []
    for lib_key, library in _normalize_libraries(gene_sets):
        ssgsea_dir = Path(output_dir) / f"ssgsea_{lib_key}"
        existing = find_ssgsea_report(ssgsea_dir)
        if existing is not None and not force:
            flush(f"  ssGSEA ({library}) exists, skipping: {existing.name}")
            rows.append({
                "analysis": "ssgsea",
                "library": library,
                "status": "skipped_existing",
                "report_path": str(existing),
                "error": "",
            })
            continue

        flush(f"  Running ssGSEA ({library})...")
        try:
            gp.ssgsea(
                data=df,
                gene_sets=library,
                outdir=str(ssgsea_dir),
                sample_norm_method="rank",
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutation_num or None,
                no_plot=False,
                seed=seed,
                threads=threads,
                verbose=False,
            )
            report = find_ssgsea_report(ssgsea_dir)
            if report is None:
                raise RuntimeError("gseapy completed but no ssGSEA report was created")
            status = "completed"
            error = ""
            flush(f"  ssGSEA ({library}) completed.")
        except Exception as exc:  # keep other libraries running, but never hide failure
            report = None
            status = "failed"
            error = f"{type(exc).__name__}: {exc}"
            flush(f"  [FAIL] ssGSEA ({library}): {error}")
        rows.append({
            "analysis": "ssgsea",
            "library": library,
            "status": status,
            "report_path": str(report) if report else "",
            "error": error,
        })

    status_df = pd.DataFrame(rows)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    status_df.to_csv(Path(output_dir) / "ssgsea_status.csv", index=False)
    return status_df


def build_ranked_list(
    results_df: pd.DataFrame,
    comparison: str,
    *,
    metric: str = "moderated_t",
) -> pd.DataFrame:
    """Build a decreasing, finite all-protein ranking for GSEA."""

    protein_col = "Protein"
    if protein_col not in results_df.columns:
        raise ValueError("Differential result table has no 'Protein' column")
    fc_col = f"Log2FC_{comparison}"
    pvalue_col = f"Pval_{comparison}"
    t_col = f"Tstat_{comparison}"

    if metric == "moderated_t" and t_col in results_df.columns:
        score = pd.to_numeric(results_df[t_col], errors="coerce")
    elif metric in {"moderated_t", "signed_log10_p"} and {
        fc_col, pvalue_col,
    }.issubset(results_df.columns):
        pvalues = pd.to_numeric(results_df[pvalue_col], errors="coerce").clip(
            lower=np.finfo(float).tiny,
        )
        fold_changes = pd.to_numeric(results_df[fc_col], errors="coerce")
        score = -np.log10(pvalues) * np.sign(fold_changes)
    elif fc_col in results_df.columns:
        score = pd.to_numeric(results_df[fc_col], errors="coerce")
    else:
        raise ValueError(f"No ranking statistic available for comparison '{comparison}'")

    ranked = pd.DataFrame({
        "Protein": results_df[protein_col].astype(str),
        "ranking_score": score,
    })
    ranked = ranked.replace([np.inf, -np.inf], np.nan).dropna()
    ranked = ranked[ranked["Protein"].str.len() > 0]
    ranked["_abs_score"] = ranked["ranking_score"].abs()
    ranked = ranked.sort_values("_abs_score", ascending=False).drop_duplicates("Protein")
    ranked = ranked.drop(columns="_abs_score").sort_values(
        "ranking_score", ascending=False,
    ).reset_index(drop=True)
    return ranked


def run_prerank_gsea(
    results_df: pd.DataFrame,
    output_dir: str,
    comparisons: list[tuple[str, str, str]],
    *,
    gene_sets: str | Sequence[str] | None = None,
    min_size: int = 10,
    max_size: int = 500,
    permutation_num: int = 1000,
    seed: int = 42,
    threads: int = 1,
    ranking_metric: str = "moderated_t",
    force: bool = False,
    flush: Callable[[str], None] = print,
) -> pd.DataFrame:
    """Run comparison-level pre-ranked GSEA on all finite protein statistics."""

    import gseapy as gp

    root = Path(output_dir) / "gsea_prerank"
    rows: list[dict[str, Any]] = []
    for _numerator, _denominator, comparison in comparisons:
        ranked = build_ranked_list(
            results_df, comparison, metric=ranking_metric,
        )
        comparison_dir = root / comparison
        comparison_dir.mkdir(parents=True, exist_ok=True)
        rank_path = comparison_dir / "ranked_proteins.csv"
        ranked.to_csv(rank_path, index=False)
        if len(ranked) < min_size:
            message = f"only {len(ranked)} finite ranked proteins"
            flush(f"  [SKIP] pre-ranked GSEA {comparison}: {message}")
            rows.append({
                "analysis": "prerank_gsea",
                "comparison": comparison,
                "library": "",
                "n_ranked_proteins": len(ranked),
                "status": "insufficient_ranked_features",
                "report_path": "",
                "error": message,
            })
            continue

        for lib_key, library in _normalize_libraries(gene_sets):
            gsea_dir = comparison_dir / lib_key
            existing = find_prerank_report(gsea_dir)
            if existing is not None and not force:
                flush(f"  pre-ranked GSEA ({comparison}, {library}) exists, skipping.")
                rows.append({
                    "analysis": "prerank_gsea",
                    "comparison": comparison,
                    "library": library,
                    "n_ranked_proteins": len(ranked),
                    "status": "skipped_existing",
                    "report_path": str(existing),
                    "error": "",
                })
                continue

            flush(f"  Running pre-ranked GSEA ({comparison}, {library})...")
            try:
                gp.prerank(
                    rnk=ranked[["Protein", "ranking_score"]],
                    gene_sets=library,
                    outdir=str(gsea_dir),
                    min_size=min_size,
                    max_size=max_size,
                    permutation_num=permutation_num,
                    seed=seed,
                    threads=threads,
                    ascending=False,
                    no_plot=False,
                    verbose=False,
                )
                report = find_prerank_report(gsea_dir)
                if report is None:
                    raise RuntimeError("gseapy completed but no pre-ranked report was created")
                status = "completed"
                error = ""
                flush(f"  pre-ranked GSEA ({comparison}, {library}) completed.")
            except Exception as exc:
                report = None
                status = "failed"
                error = f"{type(exc).__name__}: {exc}"
                flush(f"  [FAIL] pre-ranked GSEA ({comparison}, {library}): {error}")
            rows.append({
                "analysis": "prerank_gsea",
                "comparison": comparison,
                "library": library,
                "n_ranked_proteins": len(ranked),
                "status": status,
                "report_path": str(report) if report else "",
                "error": error,
            })

    status_df = pd.DataFrame(rows)
    root.mkdir(parents=True, exist_ok=True)
    status_df.to_csv(root / "prerank_status.csv", index=False)
    return status_df


def run_ssgsea_from_config(
    df: pd.DataFrame,
    output_dir: str,
    cfg: "ProteomicsConfig",
    *,
    gene_sets: str | Sequence[str] | None = None,
    flush: Callable[[str], None] = print,
) -> pd.DataFrame:
    return run_ssgsea(
        df,
        output_dir,
        gene_sets=gene_sets or cfg.gene_sets,
        min_size=cfg.gsea_min_set_size,
        max_size=cfg.gsea_max_set_size,
        permutation_num=cfg.ssgsea_permutation_num,
        seed=cfg.gsea_seed,
        threads=cfg.gsea_threads,
        flush=flush,
    )


__all__ = [
    "DEFAULT_GENE_SETS",
    "build_ranked_list",
    "find_prerank_report",
    "find_ssgsea_report",
    "run_prerank_gsea",
    "run_ssgsea",
    "run_ssgsea_from_config",
]
