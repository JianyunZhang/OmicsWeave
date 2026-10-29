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
    results: analysis_results
    references: references
    analysis: analysis              # 课题专用脚本
  references: {{}}                    # 例如 transcriptome: references/Mus_musculus.GRCm39.cdna.all.fa.gz
defaults:
  padj_method: fdr_bh
  alpha: 0.05
  random_seed: 20261029
experiments: []                     # 留空时自动查找 experiments/ 下含 experiment.yaml 的目录
"""

EXPERIMENT_TEMPLATE = """\
# OmicsWeave 实验配置。文件路径（input.files、samples.tsv 的 files 列）均相对于项目目录。
schema_version: 1
experiment:
  id: {experiment_id}               # 必须与所在目录名一致
  description: ""
  metadata: {{}}                      # 例如 tracer: U-13C glucose

design:
  subjects: subjects.tsv            # 每只动物一行
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
  reversal_pairs:                   # 疾病效应被治疗逆转时，候选结果会标记 reversal
    - {{disease: HFD_vs_LFD, treatment: Drug_vs_Saline_in_HFD}}
  stratify_by: [tissue]             # 各组织分别分析
  covariates: []                    # 例如 [batch]

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
    # input: {{level: table, format: legacy_mfa_xlsx, files: [experiments/{experiment_id}/mfa/flux_raw_data.xlsx]}}

integration:
  enabled: false
  modes: []                         # flux_anchored / pathway_convergence / sample_level
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
    for directory in (exp_dir, root / "analysis_results", root / "references", root / "analysis"):
        directory.mkdir(parents=True, exist_ok=True)

    files = {
        project_file: PROJECT_TEMPLATE.format(project_id=project_id, species=Species(species).value),
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
