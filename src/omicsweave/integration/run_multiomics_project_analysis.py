#!/usr/bin/env python3
"""项目级多实验多组学整合命令行入口。"""

from __future__ import annotations

import argparse
from pathlib import Path


from omicsweave.integration.project_integration import run_project_integration


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Comparable-experiment multi-omics random-effects integration"
    )
    parser.add_argument("--config", required=True, help="Path to project_config.yaml")
    args = parser.parse_args()
    result = run_project_integration(args.config)
    print(f"[OK] Project report: {result.report_path}")
    print(f"[OK] Meta results: {len(result.meta_results)}")


if __name__ == "__main__":
    main()
