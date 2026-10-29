"""代谢组学配置驱动流程的唯一 CLI 入口。"""
from __future__ import annotations

import argparse
import os
from pathlib import Path



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="代谢组学分析流程（PQN/QC/PLS-DA/Welch-FDR/ORA/GSEA）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python -m omicsweave.metabolomics.run_metabolomics_analysis --experiment 20250425_cimr_wt_mice_dio
  python -m omicsweave.metabolomics.run_metabolomics_analysis --config experiments/.../metabolomics_data/config.yaml
  python -m omicsweave.metabolomics.run_metabolomics_analysis --experiment <id> --steps 1,2,3

步骤:
  1=预处理与QC过滤  2=PCA/PLS-DA验证  3=差异分析
  4=无监督热图      5=KEGG ORA        6=GSEA
        """,
    )
    parser.add_argument(
        "--experiment",
        default=None,
        help="实验 ID，对应 experiments/<id>/metabolomics_data/config.yaml",
    )
    parser.add_argument("--config", default=None, help="config.yaml 的显式路径")
    parser.add_argument("--outdir", default=None, help="覆盖 paths.output_dir")
    parser.add_argument(
        "--steps",
        default=None,
        help="强制重跑的步骤编号，逗号分隔；未指定时按 manifest 增量执行全部步骤",
    )
    parser.add_argument(
        "--skip-report",
        action="store_true",
        help="不生成 Step 7 HTML/Markdown 解读报告",
    )

    # 仅用于给旧调用者明确迁移错误；不再执行任何硬编码分组。
    parser.add_argument("--input", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--output", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--legacy-mode", action="store_true", help=argparse.SUPPRESS)
    return parser


def parse_steps(raw: str | None) -> list[int] | None:
    if raw is None:
        return None
    try:
        values = [int(value.strip()) for value in raw.split(",") if value.strip()]
    except ValueError as exc:
        raise ValueError("--steps 必须是逗号分隔的整数 1-6") from exc
    if not values or any(value not in range(1, 7) for value in values):
        raise ValueError("--steps 只能包含 1-6，且不得为空")
    return list(dict.fromkeys(values))


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.legacy_mode or args.input or args.output:
        parser.error(
            "--input/--output/--legacy-mode 已移除，因为旧模式硬编码样本分组且统计口径不完整。"
            "请创建 experiments/<id>/metabolomics_data/config.yaml，并使用 --experiment 或 --config；"
            "可用 --outdir 覆盖输出目录。"
        )
    if not args.experiment and not args.config:
        parser.error("请指定 --experiment <id> 或 --config <path>")
    try:
        steps = parse_steps(args.steps)
    except ValueError as exc:
        parser.error(str(exc))

    from omicsweave.metabolomics.experiment_config import load_metabolomics_config
    from omicsweave.metabolomics.metabolomics_pipeline import run_pipeline

    cfg = load_metabolomics_config(
        experiment_id=args.experiment,
        config_path=args.config,
        outdir_override=args.outdir,
    )
    run_pipeline(cfg, steps=steps)

    if not args.skip_report:
        from omicsweave.metabolomics.generate_interpretation_report import generate_html_report

        generate_html_report(cfg, str(cfg.output_dir))


if __name__ == "__main__":
    main()
