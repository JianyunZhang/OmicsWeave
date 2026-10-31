# OmicsWeave
OmicsWeave — A Modular Python Toolkit for Multi-Omics Data Analysis and Integration

> 状态：**0.1.0.dev0，开发中**。接口尚未稳定，暂未发布到 PyPI。

OmicsWeave 提供 MFA（13C 代谢通量分析）、Bulk RNA-seq、代谢组学、蛋白质组学的下游统计与可视化，
以及以 MFA 通量为核心的多组学证据整合。各模块说明见 [src/omicsweave/README.md](src/omicsweave/README.md)，
架构规划见 [docs/omicsweave_arch_design.md](docs/omicsweave_arch_design.md)。

| 模块 | 包路径 | extras |
|---|---|---|
| MFA | `omicsweave.mfa` | `mfa` |
| 转录组学 | `omicsweave.transcriptomics` | `transcriptomics` |
| 代谢组学 | `omicsweave.metabolomics` | `metabolomics` |
| 蛋白质组学 | `omicsweave.proteomics` | `proteomics` |
| 多组学整合 | `omicsweave.integration` | `integration` |
| 项目配置与校验、统计方法、ID 映射 | `omicsweave.core` | — |
| 数据契约（分析层数据集、统一结果表、运行清单） | `omicsweave.contract` | — |

## 安装

需要 Python ≥ 3.10。在本仓库根目录执行（`-e` 为可编辑安装，修改代码后无需重装）：

```bash
pip install -e ".[all]"
```

只需要部分模块时，按需选择 extras，例如：

```bash
pip install -e ".[mfa,integration]"
```

FASTQ 修剪（cutadapt）在 `upstream` extras 中；kallisto、DIA-NN 等外部程序不随包分发，
需自行安装，并在实验配置中填写其路径。

## 项目目录

工具包与实验数据分开存放。实验配置、原始数据和分析结果放在**实验项目目录**中：

```text
<项目目录>/
├── project.yaml                           # 项目：物种、参考文件、默认参数
├── experiments/<experiment_id>/
│   ├── experiment.yaml                    # 实验因素、比较、启用的组学与输入、整合模式
│   ├── subjects.tsv                       # 每只动物一行（组别、处理等）
│   └── samples.tsv                        # 每次组学测量一行（组织、时间点、批次、文件）
└── analysis_results/<experiment_id>/...   # 分析输出
```

新建项目用 `omicsweave init` 生成模板，填好后用 `omicsweave validate` 检查；完整示例见
[examples/demo_project](examples/demo_project)（合成数据），字段说明见 v2 设计文档第 5 节。
已有旧版各组学配置（`mfa_data/`、`transcriptomics_data/` 等）的实验，可用
`omicsweave migrate --experiment <ID>` 转换为上述格式，待补充的内容写在生成的 `MIGRATION.md` 中
（只读取旧配置和数据文件中的样本信息，不改动、不复制任何文件；见 v2 设计文档第 5.7 节）。

配置里的相对路径都相对于项目目录解析。项目目录按以下顺序确定：

1. 命令行 `omicsweave --project-dir <目录> ...`
2. 环境变量 `OMICSWEAVE_PROJECT_DIR`
3. 从当前工作目录向上查找第一个包含 `experiments/` 的目录
4. 当前工作目录

## 使用

```bash
omicsweave --help
omicsweave --project-dir D:/lab_project init --id dio_mouse
omicsweave --project-dir D:/lab_project validate
omicsweave --project-dir D:/lab_project run --experiment <experiment_id> --omics mfa
omicsweave mfa run --experiment <experiment_id>      # 旧版流程，读取旧版配置
omicsweave mfa run --help
```

项目命令：

| 命令 | 作用 |
|---|---|
| `init [--id 项目ID] [--species mus_musculus\|homo_sapiens] [--experiment 实验ID]` | 生成项目模板（不会覆盖已有文件，除非加 `--force`） |
| `validate [--experiment 实验ID ...]` | 检查配置、样本表、因素取值、每个比较的生物学重复数、输入文件与整合前提；有错误时退出码为 1 |
| `migrate --experiment 实验ID [--dry-run] [--force]` | 把旧版各组学配置转换为上述格式，待补充项写在 `MIGRATION.md` |
| `run --experiment 实验ID [--omics 组学\|all] [--steps 步骤]` | 按 `experiment.yaml` 运行，结果包写到 `analysis_results/<实验ID>/<run_id>/`；目前已接入 MFA、代谢组与转录组（表格层输入），蛋白组暂用 `omicsweave proteomics run` |

各组学流程的命令格式为 `omicsweave [--project-dir DIR] <模块> <动作> [参数]`，动作之后的参数原样交给对应流程。
这些流程目前仍读取旧的各组学配置（`experiments/<id>/<组学>_data/`），阶段 2 改为读取上面的统一项目配置。
可用的模块与动作：

| 模块 | 动作 |
|---|---|
| `mfa` | `run`、`plot`、`report`；`template`、`check`、`convert`（MFA 求解结果的标准模板、检查与旧格式转换） |
| `transcriptomics` | `run`、`plot`、`report` |
| `metabolomics` | `run` |
| `proteomics` | `run`、`plot`、`report` |
| `integration` | `run`、`harmonize`、`project`、`cross-tissue` |

未列入命令行的脚本可用 `python -m` 运行，例如 `python -m omicsweave.transcriptomics.run_gsea --experiment <id>`。

在 Python 中直接调用：

```python
from omicsweave.mfa.mfa_experiment_config import load_mfa_config
from omicsweave.mfa.run_mfa_analysis import run_from_mfa_config

run_from_mfa_config(load_mfa_config(experiment_id="<experiment_id>"))
```

## 开发

```bash
pip install -e ".[all,dev]"
python -m pytest
```

少数测试需要真实实验数据，未找到时自动跳过。把课题数据副本放到本地工作区 `local/projects/<课题名>/`（不提交，见 [local/README.md](local/README.md)），
再运行 `python -m pytest --project-dir local/projects/<课题名>` 即可运行它们；也可以设置环境变量 `OMICSWEAVE_PROJECT_DIR`。
CI（`.github/workflows/ci.yml`）在 Linux 与 Windows、Python 3.10 与 3.12 上运行全部测试，并检查仅安装基础依赖时的导入。

```text
OmicsWeave/
├── pyproject.toml          # 包配置、依赖与 extras
├── src/omicsweave/         # 源代码
├── tests/<模块>/           # 测试
├── examples/demo_project/  # 合成示例项目
├── docs/                   # 架构与方法文档
├── local/                  # 本地课题工作区：课题数据副本与课题专用脚本（不提交）
└── .github/workflows/      # CI
```

## 许可证

[BSD 3-Clause License](LICENSE)，Copyright (c) 2026, The OmicsWeave Developers。
