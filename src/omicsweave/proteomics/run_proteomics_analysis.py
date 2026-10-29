"""Deprecated compatibility entry point for the v6 proteomics pipeline.

The historical hard-coded implementation was removed because it used an
unseeded imputation-first test path. Existing command shortcuts may keep using
this filename, but must now provide the same CLI arguments as
``run_proteomics_pipeline.py``.
"""

from __future__ import annotations

import warnings

from .run_proteomics_pipeline import main  # type: ignore


if __name__ == "__main__":
    warnings.warn(
        "run_proteomics_analysis.py is deprecated; use "
        "run_proteomics_pipeline.py instead.",
        FutureWarning,
        stacklevel=1,
    )
    main()
