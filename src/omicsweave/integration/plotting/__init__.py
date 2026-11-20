"""Multi-omics visualization modules (lazy exports to avoid eager seaborn/global rc side-effects)."""
from __future__ import annotations

from importlib import import_module

_EXPORTS: dict[str, tuple[str, str]] = {
    "apply_config": ("._config", "apply_config"),
    "plot_data_availability": (".overview_plots", "plot_data_availability"),
    "plot_differential_summary": (".overview_plots", "plot_differential_summary"),
    "plot_tissue_dashboard": (".overview_plots", "plot_tissue_dashboard"),
    "plot_mfa_volcanoes": (".volcano_plots", "plot_mfa_volcanoes"),
    "plot_pathway_condition_heatmap": (".heatmap_plots", "plot_pathway_condition_heatmap"),
    "plot_tca_cycle": (".heatmap_plots", "plot_tca_cycle"),
    "plot_convergence_heatmap": (".heatmap_plots", "plot_convergence_heatmap"),
    "plot_tissue_flux_overview": (".heatmap_plots", "plot_tissue_flux_overview"),
    "plot_evidence_profile": (".scoring_plots", "plot_evidence_profile"),
    "plot_pathway_map": (".pathway_map_plots", "plot_pathway_map"),
    "SUBMAP_CATALOG": (".pathway_submaps", "SUBMAP_CATALOG"),
    "get_submap": (".pathway_submaps", "get_submap"),
    "select_submaps_for_view": (".pathway_submaps", "select_submaps_for_view"),
    "plot_matrix": (".dossier_plots", "plot_matrix"),
    "plot_dossier_page": (".dossier_plots", "plot_dossier_page"),
    "plot_node_split": (".dossier_plots", "plot_node_split"),
    "plot_tracer_coverage": (".dossier_plots", "plot_tracer_coverage"),
}

__all__ = list(_EXPORTS)


def __getattr__(name: str):
    if name in _EXPORTS:
        mod_name, attr = _EXPORTS[name]
        module = import_module(mod_name, __name__)
        return getattr(module, attr)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
