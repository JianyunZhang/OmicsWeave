"""
Plotting subpackage for proteomics.

Provides volcano, heatmap, PCA, and pathway visualization modules.
"""
from .volcano import generate_all as generate_volcano_plots
from .heatmap import generate_all as generate_heatmap_plots
from .pca import generate_all as generate_pca_plots
from .pathway import generate_all as generate_pathway_plots
from ._config import apply_config, load_dp_results, load_combined_dp, flush_print

__all__ = [
    "generate_volcano_plots",
    "generate_heatmap_plots",
    "generate_pca_plots",
    "generate_pathway_plots",
    "apply_config",
    "load_dp_results",
    "load_combined_dp",
    "flush_print",
]