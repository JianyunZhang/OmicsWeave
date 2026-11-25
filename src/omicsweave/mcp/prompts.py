"""用于本机多组学项目的 MCP 工作流提示。"""
from __future__ import annotations

from mcp.server.fastmcp import FastMCP


def register_prompts(server: FastMCP) -> None:
    @server.prompt()
    def onboard_experiment(experiment_id: str, data_paths: str = "") -> str:
        """核对实验输入、配置与首轮预检。"""
        return (f"核对实验 {experiment_id} 的输入 {data_paths}。先调用 inspect_input_files 和 validate_experiment，"
                "读取对应的参数 Schema。记录组别、统计单位及缺失数据。配置变更先预览，获授权后提交。"
                "正式计算前调用 submit_pipeline(dry_run=true)。")

    @server.prompt()
    def qc_review_and_rerun(experiment_id: str, run_id: str = "latest") -> str:
        """审查异常样本证据及重新运行的范围。"""
        return (f"审查 {experiment_id} 的运行 {run_id}。读取摘要，按组织和比较查询 sample_quality 与 contrast_reliability。"
                "说明标记的依据和样本量变化，人工确认前只预览 update_sample_qc。应用 QC 更新后从 run 重算并用 diff_runs 对照。")

    @server.prompt()
    def method_sensitivity_study(experiment_id: str, source_run_id: str, override_variants: str = "") -> str:
        """登记并比较预先声明的敏感性方案。"""
        return (f"以 {experiment_id} 的 {source_run_id} 为基线，记录候选方案 {override_variants} 和判断标准。"
                "先调用 diagnose_run 核对基线限制。使用 submit_pipeline 的 overrides、parent_run_id、rationale "
                "和 dry_run=true 预检。计算后固定返回的 trial_id，并以 diff_runs(trial_b=...) 对照方向、效应和覆盖。")

    @server.prompt()
    def iterate_publication_figures(experiment_id: str, tissue: str, contrast_id: str) -> str:
        """从已固定的结果包重新渲染出版图件。"""
        return (f"为 {experiment_id} 的 {tissue} / {contrast_id} 查找确切整合包。"
                "调用 render_figure，按需指定区域或子图，再核对反应方向、标签与显著性标记。")
