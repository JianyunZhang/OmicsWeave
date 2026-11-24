#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
新建项目目录模板（``omicsweave init``）。

生成的模板可以直接通过 ``omicsweave validate``：所有组学默认关闭，样本表只有表头。
"""
from __future__ import annotations

from pathlib import Path

from .config import EXPERIMENT_FILE, PROJECT_FILE
from .model import Species

PROJECT_TEMPLATE = """\
# OmicsWeave 项目配置。字段说明见 docs/omicsweave_arch_design_v2.md 第 5 节。
schema_version: 1
project:
  id: {project_id}
  title: ""
  species: {species}                # mus_musculus / homo_sapiens
  paths:                            # 均可省略，以下为默认值（相对于本目录）
    experiments: experiments
    results: results
    references: references
    analysis: analysis              # 课题专用脚本
  references: {{}}                    # 例如 transcriptome: references/Mus_musculus.GRCm39.cdna.all.fa.gz
defaults:
  padj_method: fdr_bh
  alpha: 0.05
  random_seed: 20261029
experiments: []                     # 留空时自动查找 experiments/ 下含 experiment.yaml 的目录
# omicsweave pipeline 的标准流程（validate → run → integrate → 课题脚本），可省略：
# pipeline:
#   experiments: [{experiment_id}]
#   scripts:
#     - command: analysis/plot_results.py --experiment {{experiment}} --run {{run}}
#   retention:                        # 交付保留策略（默认关闭，也可在 run_analysis.py 中配置）
#     enabled: false
#     keep_last: 5
#     drop_failed: true
#     prune_orphan_bundles: false
#     stale_running_hours: 24
"""

RUN_ANALYSIS_TEMPLATE = """\
\"\"\"运行本课题 project.yaml 配置的标准流程，并提供 latest / clean 快捷入口。

常用示例：
  python run_analysis.py                     # 运行完整 pipeline（默认关闭自动清理）
  python run_analysis.py --dry-run           # 预检，不运行、不写文件
  python run_analysis.py --tag v1 --pin      # 给本次交付加标签并锁定，免于被清理
  python run_analysis.py latest              # 查询最新成功交付（支持 --open / --path / --json / --list）
  python run_analysis.py clean --dry-run     # 预览可清理项（需先将下方 CLEAN_COMMAND_ENABLED 设为 True）
  python run_analysis.py clean --apply       # 执行清理（需先将下方 CLEAN_COMMAND_ENABLED 设为 True）
\"\"\"
from pathlib import Path
import sys

# ── 交付保留策略与清理开关（默认全部关闭）──────────────────────────────────
# 仅当设为 True 时，pipeline 成功后才按下方规则自动淘汰旧交付
RETENTION_ENABLED = False
# 仅当设为 True 时，才允许通过本入口执行 `clean` 清理命令
CLEAN_COMMAND_ENABLED = False

# 保留策略参数（仅在启用上述任一开关后生效）
KEEP_LAST_DELIVERIES = 5          # 每个实验保留最近 N 次成功交付（>= 1）；锁定交付始终保留
DROP_FAILED_DELIVERIES = True     # 清理时是否一并移除失败及超时中断的交付目录
PRUNE_ORPHAN_BUNDLES = False      # 清理时是否联动移除不再被保留交付引用的孤立 run_* 结果包
STALE_RUNNING_HOURS = 24.0        # 超过该小时数仍处于 running 状态的交付视为已中断


if __name__ == "__main__":
    project = Path(__file__).resolve().parent
    for parent in project.parents:
        if (parent / "src" / "omicsweave" / "cli.py").is_file():
            sys.path.insert(0, str(parent / "src"))
            break
    from omicsweave.cli import project_pipeline_main
    from omicsweave.core.retention import RetentionPolicy

    retention = RetentionPolicy(
        enabled=RETENTION_ENABLED,
        keep_last=KEEP_LAST_DELIVERIES,
        drop_failed=DROP_FAILED_DELIVERIES,
        prune_orphan_bundles=PRUNE_ORPHAN_BUNDLES,
        stale_running_hours=STALE_RUNNING_HOURS,
    )
    project_pipeline_main(
        project,
        retention=retention,
        allow_clean=CLEAN_COMMAND_ENABLED,
    )
"""

EXPERIMENT_TEMPLATE = """\
# OmicsWeave 实验配置。文件路径（input.files、samples.tsv 的 files 列）均相对于项目目录。
schema_version: 1
experiment:
  id: {experiment_id}               # 必须与所在目录名一致
  description: ""
  metadata: {{}}                      # 例如 tracer: U-13C glucose

design:
  subjects: subjects.tsv            # 每只动物一行；UTF-8 制表符分隔，也可写 .csv（逗号分隔）
  samples: samples.tsv              # 每次组学测量一行
  factors:                          # 声明实验因素：所在表格（subject / sample）、取值与参照水平
    diet:      {{level: subject, levels: [LFD, HFD], reference: LFD}}
    treatment: {{level: subject, levels: [Saline, Drug], reference: Saline}}
    timepoint: {{level: sample, type: numeric, unit: week}}
  contrasts:                        # role：disease_effect / treatment_effect / time_effect / other
    - {{id: HFD_vs_LFD, factor: diet, numerator: HFD, denominator: LFD,
       subset: {{treatment: Saline}}, role: disease_effect}}
    - {{id: Drug_vs_Saline_in_HFD, factor: treatment, numerator: Drug, denominator: Saline,
       subset: {{diet: HFD}}, role: treatment_effect}}
  reversal_pairs:                   # 整合模式 reversal 逐特征判定处理是否逆转了疾病效应
    - {{disease: HFD_vs_LFD, treatment: Drug_vs_Saline_in_HFD}}
  stratify_by: [tissue]             # 各组织分别分析
  # tissue_aliases: {{}}              # 实验级组织别名 → 规范名，如 {{white_adipose: sq_adipose, cerebral_cortex: brain, renal_cortex: kidney}}
  covariates: []                    # 例如 [batch]
  models: []                        # 公式模型，如 [{{id: diet_pooled, formula: "~ diet + treatment", report: [diet]}}]
  adaptive_shrinkage: false         # true 时给差异结果加 lfsr（自适应收缩的局部假符号率）

derived: []                         # 声明式派生特征，如 [{{type: step_ratio, omics: metabolomics, source: network}}]

omics:                              # 按需启用；input.level 为 raw（下机文件）或 table（已处理的表格）
  transcriptomics:
    enabled: false
    # input: {{level: raw, format: fastq_gz, layout: paired, strandedness: reverse}}
    # adapter: {{name: kallisto, executable: kallisto}}
  proteomics:
    enabled: false
    # input: {{level: raw, format: thermo_raw, acquisition: dia}}
  metabolomics:
    enabled: false
    # input: {{level: table, format: peak_table, source: msdial, files: [experiments/{experiment_id}/data/peaks.tsv]}}
  mfa:
    enabled: false
    # 旧版自研软件导出：files 依次为通量文件与可选的 MID 文件；样本名需与 samples.tsv 中 omics = mfa 的 sample_id 一致
    # input: {{level: table, format: legacy_mfa_xlsx, model_id: central_carbon_v1, flux_units: "nmol/min/g tissue",
    #          files: [experiments/{experiment_id}/mfa/flux_raw_data.xlsx, experiments/{experiment_id}/mfa/mid_raw_data.xlsx]}}
    # 标准格式（omicsweave mfa template 生成模板）：format: omicsweave_mfa_xlsx；全部拟合结果用 Parquet 目录：format: mfa_tables
    # 单示踪剂无需额外配置；多示踪剂（同一模型）在 files 中每种示踪剂列一组结果，并在 samples.tsv 加 protocol 列写示踪剂

integration:
  enabled: false
  modes: []                         # flux_anchored / pathway_convergence / reversal / reaction_dossier
"""

SUBJECTS_HEADER = "subject_id\tdiet\ttreatment\tsex\n"
SAMPLES_HEADER = "sample_id\tsubject_id\tomics\ttissue\ttimepoint\tbatch\tsample_type\tfiles\n"


def create_project(root: str | Path, project_id: str, species: Species = Species.MOUSE,
                   experiment_id: str = "exp01", force: bool = False) -> list[Path]:
    """在 ``root`` 下生成项目模板，返回写出的文件。已有 project.yaml 时除非 ``force`` 否则报错。"""
    root = Path(root)
    project_file = root / PROJECT_FILE
    if project_file.exists() and not force:
        raise FileExistsError(f"{project_file} 已存在；如需覆盖请加 --force")

    exp_dir = root / "experiments" / experiment_id
    for directory in (exp_dir, root / "results", root / "references", root / "analysis"):
        directory.mkdir(parents=True, exist_ok=True)

    files = {
        project_file: PROJECT_TEMPLATE.format(project_id=project_id, species=Species(species).value,
                                              experiment_id=experiment_id),
        root / "run_analysis.py": RUN_ANALYSIS_TEMPLATE,
        exp_dir / EXPERIMENT_FILE: EXPERIMENT_TEMPLATE.format(experiment_id=experiment_id),
        exp_dir / "subjects.tsv": SUBJECTS_HEADER,
        exp_dir / "samples.tsv": SAMPLES_HEADER,
    }
    written = []
    for path, text in files.items():
        if path.exists() and not force:
            continue
        path.write_text(text, encoding="utf-8")
        written.append(path)
    return written
