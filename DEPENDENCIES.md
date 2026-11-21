# OmicsWeave 依赖包清单与环境配置说明

本文档列出运行 OmicsWeave 完整多组学分析、整合与可视化流程所需的 Python 版本、第三方依赖包及外部程序清单（以根目录 [`pyproject.toml`](pyproject.toml) 为权威来源）。

---

## 1. Python 版本与环境建议

- **Python 版本**：`>= 3.10`（推荐 **Python 3.12**）。
- **环境隔离**：建议使用独立的 Conda 环境或虚拟环境（例如 `conda create -n omicsweave python=3.12`），**避免直接混用系统或 `Anaconda3` 的 `base` 环境**，以防包版本冲突或缺少可选分析依赖（如 `pydeseq2`）。

---

## 2. 独立环境创建与安装（推荐）

### 方式 A：使用仓库根目录的 `environment.yml`（推荐，适用于 Conda / Miniconda / Anaconda）

直接在仓库根目录执行以下命令创建干净的 `omicsweave` 独立环境（仅安装 Python 3.12 + 统一 ABI 的完整依赖，不带入旧版 Anaconda base 中的 `dask` / `numexpr` 等冲突包）：

```bash
conda env create -f environment.yml
conda activate omicsweave
pip install -e . --no-deps
```

后续更新环境时执行：

```bash
conda env update -f environment.yml --prune
conda activate omicsweave
pip install -e . --no-deps
```

> **PyCharm 配置**：创建完成后，在 `Settings → Project → Python Interpreter → Add Interpreter → Conda Environment → Use existing environment` 中选择 `omicsweave`，之后运行 `run_analysis.py` 即可自动使用该独立环境。

### 方式 B：手动创建 Conda 环境并按 `pyproject.toml` 安装

```bash
conda create -n omicsweave python=3.12 -y
conda activate omicsweave
pip install -e ".[all]"
```

### 方式 C：不使用 Conda（Python 内置 `venv`）

```bash
# Windows PowerShell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -e ".[all]"
```

```bash
# Linux / macOS
python3.12 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e ".[all]"
```

### 按需选择可选依赖组（extras）

如需同时运行单元测试（`pytest`）或上游 FASTQ 修剪（`cutadapt`）：

```bash
pip install -e ".[all,dev]"            # 含测试依赖
pip install -e ".[all,upstream,dev]"   # 含上游 FASTQ 修剪与测试依赖
```

若仅需部分模块，可按需选择 extras（多个用逗号分隔）：

```bash
pip install -e ".[mfa,integration]"
pip install -e ".[transcriptomics,proteomics,metabolomics]"
```

若在已有干净环境中直接补装全部第三方依赖包，也可直接运行：

```bash
pip install "numpy>=1.24" "pandas>=2.0" "scipy>=1.10" "matplotlib>=3.7" "PyYAML>=6.0" "openpyxl>=3.1" "pydantic>=2.6" "pyarrow>=14" "seaborn>=0.12" "scikit-learn>=1.2" "statsmodels>=0.14" "pydeseq2>=0.4" "gseapy>=1.1" "adjustText>=0.8"
```

---

## 3. 依赖包清单

### 3.1 核心必选依赖（`pip install -e .` 默认安装）

所有模块、配置校验、数据契约、结果交付（`summary.xlsx` / Parquet）与出版级可视化均依赖以下包：

| 包名 | 版本要求 | 主要用途 |
|---|---|---|
| `numpy` | `>=1.24` | 数值计算与矩阵运算 |
| `pandas` | `>=2.0` | 数据表处理与多组学对齐 |
| `scipy` | `>=1.10` | 统计检验（t 检验、Wilcoxon、Fisher 融合、相关性等） |
| `matplotlib` | `>=3.7` | 火山图、反应档案图 A–D、通路总图 E、子通路图与概念摘要图渲染 |
| `PyYAML` | `>=6.0` | 读取 `project.yaml`、`experiment.yaml`、`panels.yaml`、布局与映射表 |
| `openpyxl` | `>=3.1` | 读写 `.xlsx` 输入矩阵、MFA 工作簿及交付摘要 `summary.xlsx` |
| `pydantic` | `>=2.6` | 项目与实验配置 Schema 严格校验 |
| `pyarrow` | `>=14` | Parquet 统一结果表与分析层数据集读写 |

### 3.2 多组学分析与整合依赖（`[all]` 包含全部）

运行完整 `pipeline`（包含转录组、蛋白组、代谢组、MFA 与多组学整合）时**全部需要安装**：

| 包名 | 版本要求 | 所属 extras | 主要用途与缺失影响 |
|---|---|---|---|
| `pydeseq2` | `>=0.4` | `transcriptomics` | **Bulk RNA-seq 差异表达与 apeGLM LFC 收缩（每组重复 ≥ 2 时必需）**；缺失会导致转录组运行直接报错终止 |
| `gseapy` | `>=1.1` | `transcriptomics`、`proteomics`、`metabolomics` | 各组学预排序 GSEA（`gsea_prerank`）与在线基因集库获取 |
| `statsmodels` | `>=0.14` | `mfa`、`transcriptomics`、`proteomics`、`metabolomics`、`integration` | 多重检验校正（BH FDR）、回归/方差分析与多组学统计 |
| `scikit-learn` | `>=1.2` | `mfa`、`transcriptomics`、`proteomics`、`metabolomics` | PCA、PLS-DA/VIP 及多元质控分析 |
| `seaborn` | `>=0.12` | `mfa`、`transcriptomics`、`proteomics`、`metabolomics`、`integration` | 热图及部分旧版统计绘图模块 |
| `adjustText` | `>=0.8` | `transcriptomics`、`integration` | 火山图显著特征标签自动防重叠（未安装时退化为无避让标注） |

### 3.3 上游原始数据处理与开发测试依赖（按需安装）

| 包名 | 版本要求 | 所属 extras | 主要用途 |
|---|---|---|---|
| `cutadapt` | `>=4.0` | `upstream` | Bulk RNA-seq 原始 FASTQ 接头与低质量序列修剪 |
| `pytest` | `>=7.4` | `dev` | 运行 `tests/` 自动化测试套件 |

### 3.4 外部可执行程序（非 Python 包，仅原始数据定量时需要）

当实验输入为表格层（`counts_matrix`、`protein_matrix`、`peak_table`、`legacy_mfa_xlsx` 等，例如 `dio_project`）时**无需安装**以下程序；仅从原始测序/质谱文件定量时需要自行安装并在 `experiment.yaml` 中指定路径：

- **`kallisto`**：Bulk RNA-seq 原始 FASTQ 伪比对与转录本定量
- **`DIA-NN`**：DIA 蛋白质组原始质谱数据搜库与定量

---

## 4. 测试机环境自检命令

部署或重跑 `run_analysis.py` 之前，可在目标 Python 环境中运行以下命令，同时核验**包已安装**且**真实导入无 NumPy ABI / 子依赖冲突**：

```bash
python -c "import importlib.metadata as m, numpy, pandas, scipy, matplotlib, pyarrow, openpyxl, pydantic, yaml, seaborn, sklearn, statsmodels, gseapy; from pydeseq2.dds import DeseqDataSet; pkgs=['numpy','pandas','scipy','matplotlib','PyYAML','openpyxl','pydantic','pyarrow','seaborn','scikit-learn','statsmodels','pydeseq2','gseapy','adjustText']; print('\n'.join(f'{p}: {m.version(p)}' for p in pkgs)); print('ALL IMPORTS OK')"
```

若全部打印出版本号并以 `ALL IMPORTS OK` 结尾、无报错或 `_ARRAY_API` 告警，说明环境满足完整流程运行要求。也可在项目目录下先运行预检确认配置与依赖状态：

```bash
python local/projects/dio_project/run_analysis.py --dry-run
```

---

## 5. 环境注意事项

1. **`numexpr` 数值安全**：在部分 `pandas` + `numexpr` 组合下，超百万元素矩阵广播存在静默置零风险（见 `src/omicsweave/core/safety.py`）；OmicsWeave 内部已自动关闭 `compute.use_numexpr`，如无必要建议不要在环境中额外强依赖旧版 `numexpr`。
2. **Windows 路径长度**：建议将项目目录放在较浅的路径下，或在 Windows 系统策略中开启长路径支持（`LongPathsEnabled = 1`），避免深层交付目录与子通路图路径触发 260 字符限制。
3. **补装 `pydeseq2` 后的重跑方式**：若此前因缺少 `pydeseq2` 产生过失败或不完整运行，安装依赖后请从头运行 `python run_analysis.py`（不要使用 `--from scripts` 沿用旧结果包），以确保转录组差异表、下游整合与通路总图均基于 PyDESeq2 重新生成。
