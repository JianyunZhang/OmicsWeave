# 本地课题工作区（不提交）

本目录用于存放**课题数据的本地副本**和**课题专用脚本**，用来在真实数据上测试 OmicsWeave。
除本说明文件外，目录下的所有内容都已被 `.gitignore` 排除，不会进入 Git 仓库，也不会打进安装包。

## 目录约定

每个课题一个子目录，子目录本身就是一个“项目目录”：

```text
local/
└── projects/
    └── <课题名>/                      # 例如 dio_mouse
        ├── experiments/<experiment_id>/...   # 从实验项目仓库拷贝的配置与输入数据
        ├── analysis_results/                 # 在本地运行产生的结果
        ├── references/                       # 参考基因组、FASTA、通路库等（可选）
        └── analysis/                         # 课题专用脚本（如六组织通路图），以库的方式调用 omicsweave
```

原始下机文件体积较大时不必拷贝，在配置或样本表中写原位置的绝对路径即可。

## 指定使用哪个课题

| 场景 | 写法 |
|---|---|
| 命令行 | `omicsweave --project-dir local/projects/<课题名> mfa run --experiment <id>` |
| 环境变量（适合写进 PyCharm 运行配置） | `OMICSWEAVE_PROJECT_DIR=local/projects/<课题名>` |
| 在课题目录内运行 | 进入 `local/projects/<课题名>/` 后直接运行，会自动向上找到 `experiments/` |
| 测试 | `python -m pytest --project-dir local/projects/<课题名>` |

不指定时，需要真实数据的测试会自动跳过，其余测试照常运行。

## 课题专用脚本

`analysis/` 中的脚本通过 `import omicsweave` 使用工具包，自身不放进 `src/omicsweave/`。
在脚本中取项目目录：

```python
from omicsweave.core.project import project_dir

experiment_dir = project_dir() / "experiments" / "<experiment_id>"
```

## 注意

- 这里的数据和脚本可能涉及未发表结果，请勿复制到仓库的其他位置，也不要用 `git add -f` 强制提交。
- 本目录只是测试用副本，正式数据仍以实验项目仓库为准。
