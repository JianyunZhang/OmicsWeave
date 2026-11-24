# 合成示例项目

一个最小的 OmicsWeave 项目目录，用于演示配置格式，并在 CI 中检查 `omicsweave validate`。
**所有数据均为随机生成，不对应任何真实实验。**

```text
demo_project/
├── project.yaml
└── experiments/demo_glp1/
    ├── experiment.yaml          # 因素、比较（含角色与逆转配对）、启用的组学、整合模式
    ├── subjects.tsv             # 12 只小鼠：LFD-Saline、HFD-Saline、HFD-GLP1 各 4 只
    ├── samples.tsv              # 肝脏转录组 12 个 + 代谢组 12 个 + 代谢组 QC 混样 2 个
    └── data/
        ├── rna_counts.tsv       # 30 个基因的计数矩阵
        └── metabolite_peaks.tsv # 20 个代谢物的峰面积
```

在仓库根目录运行：

```bash
omicsweave --project-dir examples/demo_project validate
omicsweave --project-dir examples/demo_project pipeline
```

此配置可离线执行两层差异分析与逆转评估。合成代谢物表缺少通路注释，此示例的整合模式为 `reversal`。
完成后打开终端列出的 `summary.xlsx`。`results/` 等运行输出不会被提交（见根目录 `.gitignore`）。
