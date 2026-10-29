"""
RNA-Seq Plotting Module
========================
Visualization subpackage for RNA-Seq downstream plots (Step 5 only).

Provides:
  - Volcano plots (per tissue x comparison + summary panel)
  - DEG heatmaps (per tissue, clustermap)
  - PCA / sample clustering (global, per-tissue, correlation)

For enrichment (Step 6) and GSEA (Step 7), see the separate subpackages:
  - omicsweave/transcriptomics/enrichment/
  - omicsweave/transcriptomics/gsea/
"""

from .volcano import generate_all as generate_volcano_plots
from .heatmap import generate_all as generate_heatmap_plots
from .pca import generate_all as generate_pca_plots
from ._config import apply_config, load_de_results, flush_print
