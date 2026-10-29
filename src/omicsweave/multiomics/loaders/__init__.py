"""Multi-omics data loaders with normalized column names."""

from .mfa_loader import MFALoader
from .transcriptomics_loader import TranscriptomicsLoader
from .proteomics_loader import ProteomicsLoader
from .metabolomics_loader import MetabolomicsLoader

__all__ = [
    "MFALoader",
    "TranscriptomicsLoader",
    "ProteomicsLoader",
    "MetabolomicsLoader",
]