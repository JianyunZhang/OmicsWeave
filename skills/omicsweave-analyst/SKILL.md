---
name: omicsweave-analyst
description: 在 OmicsWeave 中发现实验、预检参数与样本表、跟踪后台作业、读取结果证据或重新渲染图件。
---

# OmicsWeave 分析员

- 先确定实验 ID、对比方向、组织和结果根。用 `list_experiments`、`inspect_input_files`、`validate_experiment` 或对应 CLI 查看输入，参数选择以 `omicsweave schema targets` 及具体 Schema 为准。
- 从 `submit_pipeline(dry_run=true)` 开始。提交计算后记录 `job_id`，按 `get_job_status` 的事件游标观察，读取终态中登记的 `run_id`。修复失败时检查交付索引及输入指纹。 `stalled` 表示 Worker 未正常落盘终态。
- 先用 `get_run_digest` 确认表与组学，再用 `query_table` 投影必要列，按 `next_offset` 分页。报告实体证据时记录 `trace_entity` 返回的表路径、运行 ID、组织和比较。
- 样本 QC 与配置修改应提供证据、预览和当前文件 SHA-256。样本表修改通常从 `run` 重算。已有结果后的面板修改登记 `post_hoc` 和依据。
- 报告重复数、数据缺失、可靠性、方向冲突和诱饵校准。效应量路径与通过统计检验的结果分别解释。以固定比较和事先描述的敏感性方案检查稳健性。
- 试验先用 `diagnose_run` 核对基线限制，登记 `parent_run_id`、`rationale` 和 `overrides`。试验完成后用 `trial_id` 读取结果，以 `diff_runs(trial_b=...)` 对照正式运行。报告 `optimization_blocked` 和 `decoy_calibration_status`，以区分检测到的风险与未测的校准。
- 重画图件时固定源包 ID。留存图件路径与使用参数，核对标签、示踪剂分层和证据等级。

完整 CLI 与 MCP 调用约定见 `docs/arch_design/ai_native_and_mcp/`。
