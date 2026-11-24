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
└── results/<experiment_id>/...            # 分析输出（latest.json 在 results/ 根目录）
```

配置（`*.yaml`）与样本表（`*.tsv`）分开：配置是层级结构、手写，样本表每行一条记录、常在 Excel 中编辑。
样本表用 UTF-8 的制表符分隔文本；在 Excel 中编辑后另存为“Unicode 文本(*.txt)”再改扩展名为 `.tsv`，
或另存为“CSV UTF-8(逗号分隔)(*.csv)”并在 `experiment.yaml` 的 `design.samples` / `design.subjects` 中写 `.csv` 文件名
（读取规则见 v2 设计文档第 5.5 节）。

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
| `run --experiment 实验ID [--omics 组学\|all] [--steps 步骤]` | 按 `experiment.yaml` 运行，结果包写到 `results/<实验ID>/<run_id>/`；四个组学均已接入（转录组、蛋白组为表格层输入），写出差异与富集结果表，并为每个组织 × 比较画火山图（`figures` 步骤，`omics/<组学>/figures/`）；其他图件与解读报告暂用 `omicsweave <组学> run` |
| `integrate --experiment 实验ID [--mode 模式] [--source 运行ID]` | 读取一次运行的统一结果表做多组学整合，结果写入新的结果包；目前已接入 `flux_anchored`、`pathway_convergence` |
| `export --experiment 实验ID [--run 运行ID] [--format xlsx\|tsv]` | 把结果包中的结果表导出为 Excel（附运行信息）或 TSV，默认最近一次运行 |
| `map --experiment 实验ID [--source 整合运行ID] [--tissue 组织] [--contrast 比较] [--layout 布局文件]` | 从已有整合包重画分组织代谢通路总图（图 E），写到 `<结果根>/<实验ID>/_redraws/<运行ID>/`，不改动原结果包 |
| `pipeline [--experiment 实验ID ...] [--from run\|integrate\|scripts] [--source 运行ID] [--tag 标签] [--pin] [--dry-run]` | 标准流程：对 `project.yaml` 中 `pipeline.experiments` 的每个实验依次 validate → run → integrate → 课题脚本（见下） |
| `latest [--experiment 实验ID ...] [--target dir\|summary\|index] [--path] [--json] [--open] [--list]` | 通过 `results/latest.json` 指针查询、列出或打开最新一次成功交付（跨平台） |
| `clean [--experiment 实验ID ...] [--keep-last N\|--keep-all] [--prune-bundles] [--dry-run\|--apply]` | 按保留策略预览或清理历史交付及可选的孤立 `run_*` 结果包（默认仅预览，加 `--apply` 才实际删除） |

### 标准流程（`pipeline`）

当 MFA 样本在 `samples.tsv` 的 `protocol` 列中登记了多种示踪剂时，一次运行会按示踪剂分别做组间比较，同时输出跨示踪剂的分层合并结果。
交付的 `summary.xlsx` 在「MFA统计与图」中链接各自的差异表和火山图，在「通路总图」中链接综合图、各示踪剂通路图和对应的 TSV。
`reaction_dossier/pathway_maps.html` 可按组织、比较并排打开各套通路图。分示踪剂页只使用该队列的 MFA 拟合和动物重复，沿用同一网络布局；综合页保留 pooled MFA 与其他组学的证据。
某个比较缺少对应示踪剂时，该格显示「无比较数据」。如只需原有综合图，可将 `integration.parameters.reaction_dossier.pathway_map.by_tracer` 设为 `false`；这项设置只控制通路出图。

在 `project.yaml` 中写好 `pipeline` 一节后，一条命令完成整个项目的标准分析：

```yaml
pipeline:
  experiments: [all_tissues]          # 默认项目中的全部实验
  scripts:                            # 课题脚本，每个实验 run / integrate 之后依次运行
    - name: figures
      command: analysis/plot_results.py --experiment {experiment} --run-dir {run_dir} --output-dir {delivery_dir}/figures
    - name: key_results
      command: analysis/export_key_results.py --run-dir {run_dir} --integration-dir {integration_dir} --output-dir {delivery_dir}/key_results
      experiments: [all_tissues]      # 只对这些实验运行
      outputs: [key_results/tables.md] # 相对于本次交付目录；文件缺失时该步骤失败
```

```bash
omicsweave --project-dir D:/lab_project pipeline --dry-run           # 先检查并列出各步
omicsweave --project-dir D:/lab_project pipeline                     # validate → run → integrate → 脚本
omicsweave --project-dir D:/lab_project pipeline --source <运行ID>   # 沿用单组学结果包，重跑整合与脚本
omicsweave --project-dir D:/lab_project pipeline --from scripts      # 只重跑脚本（最近一次结果包）
omicsweave --project-dir D:/lab_project pipeline --results-dir "D:/研究结果/课题 A"
```

- 各步之间传递确切的运行 ID：integrate 读取本次 run 的结果包，脚本的 `{run}` / `{integration}` 是本次的运行 ID，
  其他占位符包括 `{experiment}`、`{run_dir}`、`{integration_dir}`、`{project_dir}` 和 `{delivery_dir}`。
- 实验名、脚本文件与占位符在开始前全部检查；校验有错误时一步也不运行。实验没有启用整合
  （`integration.enabled` 与 `modes`）时跳过 integrate。
- 脚本在项目目录下运行，`OMICSWEAVE_PROJECT_DIR` 指向项目目录；`.py` 脚本用运行 omicsweave 的同一个 Python 执行，
  其他命令（如 `Rscript ...`）按 PATH 查找。命令按 shell 规则拆分，路径用 `/`。
- 某一步失败时停止，并给出续跑命令（`--from integrate|scripts --source <运行ID>`）。
- 沿用已有结果包前核对它登记的输入：单组学结果包对照数据、样本表、参考与 `experiment.yaml` 中 `integration` 以外的各节
  （只改 `integration`，如打开 `panels`，仍可沿用）；整合结果包对照 `experiment.yaml`、样本表与 `panels.yaml`。
  有改动时拒绝沿用，并说明从哪一步重跑，避免交付旧的面板矩阵。

结果根目录默认取 `project.paths.results`（缺省为 `results`）。单次调用可用 `--results-dir` 覆盖，相对路径以课题根目录为基准。
每次正式执行在 `<结果根>/<实验 ID>/run_YYYYMMDD_HHMMSS[_<标签>]/`（同秒冲突时自动附 `_02` 等序号）新建一次自包含的运行与交付目录，
全部完成后原子更新 `<结果根>/latest.json` 跨平台相对路径指针（按实验索引，中途失败不覆盖 `latest.json`）；单组学、整合与交付表均在同一目录下（`--from integrate` / `scripts` 续跑时在同实验下新建目录并沿用指定源包）。
日常可用 `omicsweave latest`（加 `--open` 打开最新 `summary.xlsx`、`--path` 输出路径、`--json` 输出指针、`--list` 列出全部交付）快速访问最新交付；
重要里程碑可加 `--tag <标签> --pin`（或在目录内放 `.keep`）锁定以防被清理。
保留策略与清理功能默认关闭（需在课题入口 `run_analysis.py` 中开启 `RETENTION_ENABLED` / `CLEAN_COMMAND_ENABLED`，或在 `project.yaml` 的 `pipeline.retention` 中启用，或使用 `omicsweave clean --apply`）。
打开 `summary.xlsx` 可查看来源、各组织与比较的差异结果行数、可检验数和入选数，以及结果表和质控清单。
整合包有 `reaction_dossier` 面板产物时，矩阵（XLSX）、分流与模块检验表及矩阵、分流、示踪剂覆盖图件复制到交付目录的
`panels/<面板>/<组织>/`；摘要中的「面板结果」可按面板和组织筛选并链接到这些副本（逐反应档案页仍在整合包中）。
分组织代谢通路总图（图 E，不依赖面板）每页的 PDF 与视图表复制到 `pathway_map/<组织>/`，摘要中为「通路总图」。
有整合包时另写无面板候选表 `candidates/`（`candidates.xlsx` 与 `reactions*.tsv`、`pathways*.tsv`，摘要中为「候选概览」）：
反应取 L1 及以上的全部加每个组织 × 比较排名前 10，通路取融合显著、多层显著的全部加各层单层显著的前 5 条；
入选与面板无关，`panels` 列只作标注，并附诱饵期望、各层 QC 与样本量、各层 padj 与方向（见
`omicsweave.integration.candidates`）。`panels.yaml` 可写顶层 `selection: {mode: pre_specified | post_hoc, basis: 依据}`，
登记面板是预设的还是看过结果后选定的。
`index.json` 登记阶段状态、配置和脚本指纹、源结果包及其产物路径。脚本失败时可按提示传入
`--from scripts --source <单组学 ID> --integration-run <整合 ID>`，在新的交付目录续跑。
标准导出和面板 TSV 使用带 BOM 的 UTF-8。Excel 导入时选择制表符，带前导零的 ID 列按文本读取。

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
