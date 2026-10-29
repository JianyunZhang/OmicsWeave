#!/usr/bin/env python3
"""样本级多组学协调命令行入口。"""

from __future__ import annotations

import argparse
from pathlib import Path


from omicsweave.integration.harmonization import run_harmonization


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize, sample-map, audit missingness, and optionally correct batch"
    )
    parser.add_argument("--config", required=True, help="Path to harmonization YAML")
    args = parser.parse_args()
    result = run_harmonization(args.config)
    print(
        f"[OK] Harmonized {len(result.assays)} assay(s), "
        f"{len(result.primary_samples)} primary sample(s)"
    )


if __name__ == "__main__":
    main()
