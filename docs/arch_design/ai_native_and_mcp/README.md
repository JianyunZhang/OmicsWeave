# OmicsWeave AI Native 演进与 MCP 接口架构设计

> **文档类型**：架构设计规划套件（提案 v0.1，2026-11-21）\
> **适用范围**：OmicsWeave 面向 AI Agent 协作的 CLI 机器接口、后台作业可观测性、迭代优化闭环与 Model Context Protocol（MCP）接口设计\
> **现状依据**：基于工作树提交 `6abb9cb` 的代码核查（`cli.py`、`core/`、`contract/`、`integration/`、`visualization/`）\
> **配套文档**：[OmicsWeave 架构规划 v2](../../omicsweave_arch_design_v2.md) · [OmicsWeave × TheraWeave 总纲](../../omicsweave_theraweave_development_plan.md) · [候选发现与交付计划](../omicsweave_discovery_and_delivery_plan.md)

**2026-11-25 实施更新**。工作树复核基于 `825956e` 及已有未提交修改。Phase 1 的接口、评审决定和验证边界见[实施路线](03_implementation_roadmap.md#phase-1a1-可观测后台-cli-与结构化读写底座)。下文的现状描述对应原始提案基线。

Phase 2 已加入参数 Schema、预览式配置和样本 QC 更新，以及可选的本机 MCP 服务。启动方式和接口边界见 [MCP 接口设计](02_mcp_interface_spec.md)。

Phase 3.1 与 3.2 已加入隔离试验、配置与证据对照、谱系和只读科研诊断。使用示例及验证记录见 [Phase 3 实施路线](03_implementation_roadmap.md)。

**2026-11-26 实施与规划更新**。Phase 3.3 的本机只读交换和回环 HTTP 已在合成项目验收。下一轮优先做 TheraWeave 导入与双仓库契约测试，其后推进预检锁定、提交去重、多方案试验、面板建议审阅和作业恢复。[后续实施优先级](03_implementation_roadmap.md#phase-3-后续实施优先级)

---

## 1. 核心结论

1. **当前定性**：OmicsWeave **目前还不是完整的 AI Native 应用，但已经具备高质量的 “AI-Ready / Agent-Friendly” 工程底座**。其声明式配置（Pydantic）、代码与数据解耦（`ContextVar` 项目上下文）、不可变自描述结果包（`manifest.json` + Parquet）、交付索引（`index.json` + `latest.json`）、增量续跑（`pipeline --from`）以及逐表对照（`diff`）都天然契合 Agent 编排。
2. **核心缺口**：主要集中在“Agent 闭环交互层”——缺少统一的 CLI 结构化 JSON 协议与错误码、缺少强类型可自省的分析/整合参数 Schema、缺少细粒度的异步作业运行时与运行态探活（`status/logs/cancel`）、缺少面向 LLM 上下文窗口的分级结果读取接口（`inspect/query/trace`）、缺少带防过拟合护栏的诊断与迭代谱系，以及缺少 MCP 协议入口。
3. **演进策略**：**不重写统计与整合内核**。在既有 `core` / `contract` 之上抽取无状态的应用服务层 `omicsweave.service`，统一支撑机器友好的 `--json` CLI 与薄层 MCP Server（`omicsweave.mcp`，可选依赖），分三阶段平滑落地。

---

## 2. 文档导航

| 文档 | 定位与核心内容 | 适用读者 |
|---|---|---|
| [01_ai_native_architecture.md](01_ai_native_architecture.md) | **总体架构与核心能力改造**：现状核查矩阵、分层架构、CLI 机器协议与错误码、强类型参数 Schema 收敛、后台作业状态机与心跳探活（`core/jobs.py`）、受控输入与 `--override` 零污染试验、金字塔式结果读取、科学防过拟合护栏与运行谱系。 | 架构设计评审、Phase 1–3 核心开发 |
| [02_mcp_interface_spec.md](02_mcp_interface_spec.md) | **MCP 接口详细规范**：`omicsweave.service` 共享契约、传输模式、路径权限、Resources、候选 Tools、长任务协议和科研 Prompts。 | MCP 服务端开发、Agent 集成方 |
| [03_implementation_roadmap.md](03_implementation_roadmap.md) | **分阶段实施路线与治理**：目标包目录、Phase 1–3 任务拆解与可验证验收标准、与总纲 O1–O5 / TheraWeave 的衔接、测试矩阵、风险治理与待决策事项。 | 项目排期、验收与测试编写 |

---

## 3. 关键架构决策摘要（ADR）

| 编号 | 决策项 | 选定方案 | 核心理由 |
|---|---|---|---|
| ADR-1 | CLI 与 MCP 的代码复用 | 共享 `omicsweave.service` 应用服务层（Pydantic DTO），MCP 仅作薄协议门面 | 避免 CLI 与 MCP 逻辑分叉；避免 MCP 拼接子命令字符串解析文本。 |
| ADR-2 | 后台任务执行模型 | 独立子进程 Worker + 纯文件系统作业目录（`_jobs/<job_id>/state.json` + `events.jsonl`） | 跨 Windows/Linux 零运维依赖；天然隔离旧流程残留模块级状态与高内存峰值。 |
| ADR-3 | 长任务 MCP 交互 | Hybrid Sync/Async（`wait_seconds` 窗口内直返，超时转 `job_id` + 进度通知） | 兼顾短任务低延迟与长任务（DESeq2/GSEA/诱饵校准）不触发 MCP 超时。 |
| ADR-4 | 迭代调参防污染 | `config patch`（带 `dry_run` 与失效面预览）+ `--override`（单次试验不改主配置） | 区分“正式基线变更”与“敏感性对照试验”，所有生效参数仍完整落盘至 `manifest.json`。 |
| ADR-5 | 结果读取防上下文溢出 | L0 Digest → L1 分页条件切片（`query`）→ L2 单实体跨层溯源（`trace`）→ L3 PNG 审图 | 防止几十万行 Parquet 撑爆 LLM 上下文；PyArrow 谓词下推低内存高性能。 |
| ADR-6 | 科研诚信与防 p-hacking | 诊断引擎区分“工程修复 / 敏感性对照 / 统计红线”，事后改面板强制 `post_hoc`，记录运行谱系 | 防止 Agent 通过无约束反复换方法/调阈值刷显著性或伪造预设面板。 |
| ADR-7 | 数据保密与传输默认 | MCP 默认 `stdio` 本地进程通信，锁定项目根目录沙箱；`streamable-http` 仅作可选远程模式 | 严格落实总纲“未发表数据不外泄”约束，默认不暴露任何网络端口。 |
