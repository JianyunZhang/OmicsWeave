# OmicsWeave AI Native 总体架构与核心能力改造设计

> **文档类型**：系统架构设计（v0.1，2026-11-21）\
> **所属套件**：[AI Native 与 MCP 设计套件](README.md)\
> **配套规范**：[02_mcp_interface_spec.md](02_mcp_interface_spec.md) · [03_implementation_roadmap.md](03_implementation_roadmap.md)\
> **代码基线**：`6abb9cb`（`src/omicsweave/cli.py`、`core/`、`contract/`、`integration/`、`visualization/`）

---

## 1. 背景、定位与设计原则

### 1.1 什么是科研多组学场景下的 AI Native

本方案将 OmicsWeave 的 **AI Native** 目标定义为：**人类研究者与 AI Agent 均作为软件的一等公民（First-Class User）**，支持 Agent 端到端完成以下五步科学计算闭环：

```text
数据接入与校验 (Ingest & Validate)
        │
        ▼
后台可控执行 (Async Run / Pipeline + Progress)
        │
        ▼
运行态可观测 (Status / Heartbeat / Logs / Cancel)
        │
        ▼
分级结果读取 (Digest → Query Slice → Entity Trace → Visual Inspect)
        │
        ▼
有护栏的迭代优化 (Diagnose → Safe Patch/Override → Incremental Rerun → Diff)
```

### 1.2 核心设计原则与非目标

**必须遵守的原则**：
1. **内核不重写**：统计方法、四组学 Runner、四整合模式、数据契约和可视化渲染器全部复用，改造集中在运行时、协议层与服务抽象。
2. **机器可读优先（Machine-First）**：所有正式状态、错误、诊断、对比结果均由强类型 Pydantic 模型定义，CLI JSON 与 MCP 共用同一模型。
3. **科学严谨优先于“指标好看”**：迭代优化的目标是**修正数据/工程问题、校验方法稳健性、改进图件可读性**，绝不允许 Agent 以“刷出更多显著基因或 L1 反应”为目标进行无约束调参（p-hacking）。
4. **本地默认与数据保密**：默认零外部服务、零外网上传、零端口监听，未发表实验数据始终停留在本地项目目录。
5. **向后兼容**：既有 CLI 命令、人类可读中文输出、`manifest.json`（v0.1）、`index.json`、`latest.json` 及全部现有单元测试保持完全兼容。

**本阶段非目标**：
- 不引入重型外部调度框架（如 Celery/Redis/Airflow/K8s）；单机多任务与队列以纯文件系统（JSON/JSONL + 子进程）实现，集群调度留待总纲 O5 阶段适配。
- 不在 OmicsWeave 内内置大模型推理客户端（LLM Client）；OmicsWeave 作为**工具提供方与确定性计算引擎**被外部 Agent（Claude Code、Cursor、Antigravity、自定义科研 Agent 等）调用。
- 不改变 OmicsWeave 与 TheraWeave 的职责边界（成药性与治疗决策仍归属 TheraWeave）。

---

## 2. 现状核查与能力差距矩阵

对当前代码库（提交 `6abb9cb`）逐模块核查的结果如下：

| 模块与代码位置 | 现状优势（可直接复用） | AI Native 视角下的具体缺口 |
|---|---|---|
| **CLI 入口**<br>`src/omicsweave/cli.py` | 15 个现代项目命令（`init / validate / migrate / run / integrate / export / diff / stats / refs / replicate / pipeline / map / abstract / latest / clean`）；支持 `--project-dir` 显式作用域。 | 1. 仅 `latest` 支持 `--json`，其余 14 个命令仅输出中文文本；<br>2. 几乎所有失败都返回退出码 `1`，无结构化错误码与修复提示；<br>3. 缺少 `status / logs / cancel / inspect / query / trace / diagnose / schema / config` 子命令；<br>4. 全同步阻塞执行，无 `--background`。 |
| **项目上下文**<br>`core/project.py` | 使用 `contextvars.ContextVar` 保存项目路径，支持线程/异步协程并发隔离。 | 旧版 `<组学> <动作>` 流程仍遗留少量模块级状态（v2 文档 4.1 节），长任务需通过独立子进程 Worker 隔离。 |
| **配置与模型**<br>`core/model.py`<br>`core/config.py` | `ProjectConfig` / `ExperimentConfig` 均设 `extra="forbid"`，已有 `export_json_schema()`；样本表严格检测编码（GBK 告警）、分隔符与注释行。 | 1. `OmicsSettings.analysis` 与 `IntegrationSettings.parameters` 为无类型 `dict[str, Any]`；<br>2. 各组学与整合模式把合法参数 `DEFAULTS` 分散写在各自 `runner.py` / `modes/*.py` 内部，**仅在真正运行时才校验未知键**，`validate` 与 JSON Schema 均看不到；<br>3. 无安全配置 patch 工具，也无单次试验 `--override` 机制。 |
| **运行与交付**<br>`core/run.py`<br>`core/pipeline.py`<br>`core/delivery.py`<br>`core/retention.py` | `pipeline` 拥有完整的预检（`--dry-run`）、SHA-256 过期拦截、增量续跑（`--from`）、不可变交付目录、`index.json` 阶段状态与 `latest.json` 原子指针。 | 1. 单独调用 `run` / `integrate` 时，运行中途没有状态文件，崩溃后只留无 `manifest.json` 的半成品目录；<br>2. `index.json` 仅到阶段级，无法获知正在跑哪个组学/组织/步骤；<br>3. 无 PID、无心跳、无进程存活探测，`retention.py` 只能靠 `stale_running_hours=24` 盲判中断；<br>4. 无运行日志归档与取消信号处理。 |
| **数据契约与查询**<br>`contract/manifest.py`<br>`contract/query.py`<br>`contract/tables.py`<br>`contract/diff.py` | 强类型 Parquet 表、完整输入/输出 SHA-256 清单、Python 查询对象 `Bundle`、逐表逐列 + 分类转移矩阵 `diff_bundles()`。 | 1. `query.py` 与 `diff.py` 无 CLI JSON 接口；<br>2. 缺少按条件分页裁剪的轻量化切片查询；<br>3. `RunManifest` 无迭代谱系字段（`parent_run_id`、`intent`、`rationale`、`overrides`）。 |
| **质控与候选发现**<br>`core/reliability.py`<br>`integration/candidates.py`<br>`core/validate.py` | 已有 `contrast_reliability.parquet`、`sample_quality.parquet`、无面板候选表、诱饵期望（`decoy_L1_expected`）、`panels.yaml` 的 `selection` 登记与最小可达 p 值检查。 | 各类质量与统计警告分散在 `validate` 输出、Parquet 表和 `index.json` 中，尚未汇总成统一的机器可读诊断报告和策略护栏。 |

---

## 3. 目标分层架构

在现有分层（[v2 设计文档 4.1 节](../../omicsweave_arch_design_v2.md)）中新增 **应用服务层（`service`）**、**作业运行时（`core/jobs.py` + `core/progress.py`）** 与 **MCP 协议层（`mcp`）**：

```text
┌─────────────────────────────┐    ┌─────────────────────────────────┐
│       CLI (omicsweave)      │    │   MCP Server (omicsweave.mcp)   │
│  人类可读文本  |  --json     │    │  stdio (默认) | Streamable HTTP │
└──────────────┬──────────────┘    └────────────────┬────────────────┘
               │                                    │
               └──────────────────┬─────────────────┘
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│                omicsweave.service (无状态应用服务层)                │
│  project | config | execution | query | diagnostics | visualization│
│  - 强类型 Pydantic Request/Response DTO                            │
│  - 统一语义错误码 (ServiceError) 与安全路径校验                    │
└───────────┬─────────────────────┬─────────────────────┬────────────┘
            ▼                     ▼                     ▼
┌────────────────────┐ ┌────────────────────┐ ┌──────────────────────┐
│ omicsweave.core    │ │ core/jobs.py       │ │ omicsweave.contract  │
│ config / model     │ │ 后台子进程 Worker  │ │ manifest / tables    │
│ validate / run     │ │ state.json (原子)  │ │ query / diff / export│
│ pipeline / delivery│ │ events.jsonl (流)  │ │                      │
│ diagnostics        │ │ PID探活 / 心跳/取消│ │                      │
└───────────┬────────┘ └────────────────────┘ └──────────────────────┘
            ▼
┌────────────────────────────────────────────────────────────────────┐
│  omicsweave.{mfa, transcriptomics, metabolomics, proteomics,       │
│              integration, stats, adapters, visualization}          │
│  (计算内核保持不变，仅通过可选 ProgressReporter 汇报细粒度进度)     │
└────────────────────────────────────────────────────────────────────┘
```

**依赖规则**：
1. `mcp` 只能依赖 `service` 和 `core.project` / 路径沙箱工具，**禁止直接导入具体组学或整合内部实现**。
2. `cli` 对新增及重构的机器接口命令一律委托给 `service`，保证 CLI `--json` 与 MCP Tool 返回的数据模型完全一致。
3. `core` 不反向依赖 `service`、`cli` 或 `mcp`。
4. `core/jobs.py` 仅使用 Python 标准库（`subprocess`、`os`、`signal`、`json`、`threading`），不引入第三方队列依赖。

---

## 4. 统一机器协议层（Machine-First CLI & Parameter Schema）

### 4.1 统一 JSON 响应信封与流隔离

新增全局参数 `--format {text,json}`。已接入的命令通常也接受子命令后的 `--json`。已有的 `latest --json` 保留原始对象或数组，`latest --format json` 返回版本化信封。`map/abstract` 继续接受图件 `--format png/pdf/...`，响应信封可用 `--json` 或全局 `--format json`。导出的 `export --format xlsx/tsv` 表示文件格式。

**流隔离约定**：
- 当启用 `--format json` 时：
  - `stdout` **严格且仅输出单个合法 JSON 对象**；
  - 所有运行期过程打印、第三方库警告（如 `pydeseq2`、`gseapy`、`matplotlib` 日志）全部重定向至 `stderr` 或作业日志文件，严禁污染 `stdout`。
- 统一响应信封（`ServiceResponse`）：

```json
{
  "ok": true,
  "schema_version": 1,
  "command": "validate",
  "project_dir": "/path/to/project",
  "timestamp": "2026-11-21T10:30:00Z",
  "data": {},
  "warnings": [
    {
      "code": "UNDERPOWERED_RANK_TEST",
      "where": "实验 all_tissues",
      "message": "比较 HFD_vs_LFD 在 mfa 中，秩检验在以下分层达不到 p < 0.05：liver（3 vs 3，最小 p = 0.1）"
    }
  ],
  "error": null
}
```

失败时 `ok: false`，`data: null`，`error` 必填：

```json
{
  "ok": false,
  "schema_version": 1,
  "command": "pipeline",
  "project_dir": "/path/to/project",
  "timestamp": "2026-11-21T10:31:12Z",
  "data": null,
  "warnings": [],
  "error": {
    "code": "STALE_SOURCE_BUNDLE",
    "message": "all_tissues 的整合结果包 run_20261118_195428_674e 已过期：experiments/all_tissues/panels.yaml 在整合之后改过",
    "hint": "请从 integrate 开始续跑：omicsweave pipeline --experiment all_tissues --from integrate --source run_20261118_172720_cc0e",
    "retryable": true,
    "details": {
      "experiment": "all_tissues",
      "stale_bundle": "run_20261118_195428_674e",
      "changed_files": ["experiments/all_tissues/panels.yaml"],
      "recommended_start": "integrate",
      "source_run_id": "run_20261118_172720_cc0e"
    }
  }
}
```

### 4.2 标准化语义错误码与退出码表

| 退出码 | 错误码（`error.code`） | 典型触发场景 | Agent 推荐处置动作 |
|---:|---|---|---|
| `0` | — | 执行成功（注：`diff` 仍保持 `0` 相同 / `1` 有差异 / `2` 无法对照，JSON 中通过 `data.identical` 明确区分） | 继续下游步骤 |
| `2` | `CLI_USAGE_ERROR` | 参数解析错误、非法枚举值 | 检查 `schema` 或帮助后修正调用 |
| `10` | `CONFIG_SCHEMA_ERROR` | YAML 语法错误、Pydantic 字段校验失败、未知配置键 | 读取 `details.field_errors`，按 Schema 修正配置 |
| `11` | `VALIDATION_FAILED` | `omicsweave validate` 检出错误（缺样本、ID 不匹配、样本表编码非 UTF-8） | 按 `details.issues` 修复样本表或配置 |
| `12` | `INPUT_NOT_FOUND` | 输入数据文件、参考文件或源结果包不存在 | 调用 `inspect-input` 核对路径 |
| `13` | `STALE_SOURCE_BUNDLE` | 续跑沿用的单组学或整合包因输入/配置改动而过期 | 按 `details.recommended_start` 调整 `--from` 续跑 |
| `14` | `OPTIONAL_DEP_MISSING` | 缺少 `pydeseq2`、`gseapy`、`mcp` 等可选依赖或外部工具（kallisto/DIA-NN） | 告知用户安装对应 extras，不可盲目重试 |
| `15` | `SECURITY_POLICY_VIOLATION` | 路径越界出项目根目录、未授权触发自定义脚本或破坏性清理 | 修正路径或请求人工授权 |
| `20` | `EXECUTION_FAILED` | 组学统计、整合计算或课题脚本执行异常 | 查 `logs` 与 `index.json`，定位具体阶段失败原因 |
| `21` | `CONTRACT_ERROR` | 结果包损坏、缺 `manifest.json`、Parquet 列不符合规范 | 检查源包完整性或重跑 |
| `22` | `JOB_NOT_FOUND` / `JOB_CONFLICT` | 查询不存在的 `job_id`，或取消已结束作业 | 刷新 `status` 列表 |
| `30` | `JOB_CANCELLED` | 作业被用户或 Agent 主动取消 | 无需自动重试 |

### 4.3 强类型参数 Schema 收敛与自省（`omicsweave schema`）

**当前痛点**：各模块的可调参数散落在私有字典中（如 `proteomics/runner.py` 的 `_DEFAULTS`、`metabolomics/runner.py` 的 `MetabolomicsConfig`、`integration/modes/flux_anchored.py` / `pathway_convergence.py` / `reaction_dossier.py` / `reversal.py` 的 `DEFAULTS`、`integration/pathway_map.py` 的 `DEFAULTS`），导致 `validate` 无法静态发现拼写错误，且 Agent 无法自省参数空间。

**设计方案**：
1. 在 `core/param_schemas.py`（或各子包公开模型后统一汇总到注册表）定义带文档与范围约束的 Pydantic 模型：
   - `TranscriptomicsAnalysisParams`（`gene_aggregation`、`gene_sets`、`run_prerank_gsea`、`gsea_ranking_metric`、`min_count` 等）
   - `ProteomicsAnalysisParams`（`dp_method`、`pval_cutoff`、`log2fc_cutoff`、`fdr_method`、`protein_aggregation`、`min_observations_per_group` 等）
   - `MetabolomicsAnalysisParams`（镜像 `MetabolomicsConfig` 字段）
   - `MFAAnalysisParams`（`flux_groups`、`inference_unit`、`statistic` 等）
   - `FluxAnchoredParams`、`PathwayConvergenceParams`、`ReactionDossierParams`、`ReversalParams`、`PathwayMapParams`
2. 为每个参数增加元数据标记（通过 `json_schema_extra`）：
   - `x-tunable-class`: `"engineering"`（可安全调优，如布局/别名/格式）| `"sensitivity"`（需对照验证，如收缩/协变量/聚合）| `"statistical_threshold"`（受防 p-hacking 护栏约束，如 `pval_cutoff`、`alpha`）。
3. 校验阶段过渡策略：
   - `OmicsSettings.analysis` 与 `IntegrationSettings.parameters` 在 `validate_experiment()` 阶段对照注册表校验键名、类型与值域，将原先“运行时才报未知参数”提前到 `validate` 阶段拦截。
4. 新增 `omicsweave schema [target] --json`：
   - 支持导出 `project`、`experiment`、`omics.<name>`、`integration.<mode>`、`tables.<differential|enrichment|reaction_assessment|...>` 及 `all` 目录清单。

---

## 5. 后台异步作业引擎与运行态可观测性

### 5.1 统一作业存储布局

在实际结果根目录新增 `<results_root>/_jobs/`。交付继续使用当前 `index.json` 与 `latest.json` 布局，作业记录独立于结果包枚举。

```text
<results_root>/
├── latest.json                           # 每个实验当前成功交付指针
├── <experiment_id>/
│   └── run_YYYYMMDD_HHMMSS.../            # 结果包，若用于交付可含 index.json、summary.xlsx
└── _jobs/
    └── job_<uuid>/                       # 每次后台执行的运行时记录
        ├── spec.json                      # 不可变的执行请求
        ├── launch.json                    # 派发进程记录的 PID 与进程身份
        ├── state.json                     # 原子替换：实时状态、PID、心跳、当前进度、产出指针
        ├── events.jsonl                   # 追加写：细粒度结构化生命周期事件流
        ├── output.log                     # 合并的 stdout/stderr 文本日志
        └── cancel.flag                    # 可选：软取消请求哨兵文件
```

### 5.2 作业状态机与健康探测

```text
           submit
             │
             ▼
          [queued]
             │ worker 启动并写入 pid
             ▼
         [running] ──(心跳更新 + events.jsonl 追加)──┐
         │   │   │                                  │
    正常结束 │   │ 收到 cancel / SIGTERM             │ 进程退出但状态仍为 running
         │   │   ▼                                  │ (PID 已退出或身份不符)
         │   │ [cancelled]                          ▼
         │   │                                  [stalled]
         │   ▼ 抛出异常/退出码非 0
         │ [failed]
         ▼
     [completed]
```

**`state.json` 数据结构**：

```json
{
  "schema_version": 1,
  "job_id": "job_07f09c6e1b8c4b10925aa8919a4baa86",
  "kind": "pipeline",
  "project_dir": "/path/to/project",
  "experiments": ["all_tissues"],
  "status": "running",
  "pid": 24816,
  "host": "WORKSTATION-01",
  "created_at": "2026-11-21T10:30:12Z",
  "started_at": "2026-11-21T10:30:13Z",
  "heartbeat_at": "2026-11-21T10:32:45Z",
  "finished_at": null,
  "exit_code": null,
  "progress": {
    "stage": "run",
    "experiment": "all_tissues",
    "omics": "transcriptomics",
    "tissue": "liver",
    "step": "enrichment"
  },
  "artifacts": {
    "all_tissues": {
      "experiment": "all_tissues",
      "run_id": "run_20261121_103015_c1d2",
      "delivery_id": "run_20261121_103015_c1d2"
    }
  },
  "last_event_seq": 7,
  "result": null,
  "error": null
}
```

**跨平台 PID 探活与心跳机制**：
- **心跳线程**：Worker 进程启动一个轻量守护线程（每 `5s` 原子刷新 `state.json` 的 `heartbeat_at`），避免主线程陷入重型数值计算（如 PyDESeq2 拟合、旋转检验）时误判卡死。
- **按需探活（Read-Time Liveness Check）**：每次查询 `status` 时，若 `state.json` 显示 `running`：
  1. 检查 `host` 是否为本机；
  2. Windows 用进程句柄与创建时间，Linux 用 `/proc` 状态、启动时间与系统启动 ID。其他 POSIX 平台用信号 0 探活。
  3. 确认进程退出或 PID 身份变更，且磁盘上仍是活动状态时，在查询视图报告 `stalled`。仅心跳超时且未确认退出时报告 `unresponsive`，清理仍视为活跃。
- **协作取消**。写入 `cancel.flag` 后，由 Worker 在计算、脚本或交付的安全边界确认并写入 `cancelled`。正在执行的统计库调用和课题脚本可持续至自身返回。强制取消须待进程树归属及产物一致性机制明确后再加入。

### 5.3 细粒度进度上报（`core/progress.py`）

在 [RunContext](../../../src/omicsweave/core/run.py#L53) 中增加可选字段 `reporter: ProgressReporter = NullReporter()`：

```python
class ProgressReporter(Protocol):
    job_id: str | None

    def emit(self, event: str, *, stage: str, detail: dict | None = None) -> None: ...

    def check_cancelled(self) -> None: ...
```

**埋点位置（最小侵入）**：
- `core/pipeline.py`：每个 `Delivery.stage()` 进入/退出及课题脚本启动前；
- `core/run.py`：每个组学 Runner 启动/完成、后处理 `stats.post.postprocess()`；
- 4 个组学 `runner.py`：在 `STEPS`（`qc / differential / enrichment / figures`）与按组织循环边界调用 `ctx.reporter.emit(...)` 和 `check_cancelled()`；
- `integration/integrate.py` 及 4 个 `modes/*.py`：每个整合模式启动、诱饵重复迭代、图件渲染循环。

### 5.4 新增运行态 CLI 子命令与既有模块联动

- `omicsweave run|integrate|pipeline ... --background [--json]` 派发独立 Worker，返回 `job_id` 与初始 `state.json` 路径。完整预检由 Worker 执行。
- `omicsweave status [--job <job_id> | --experiment <exp_id> | --active] [--limit 20] [--json]`：展示作业实时状态、心跳年龄、当前步骤进度、已关联的 `run_id` / `delivery_id`。
- `omicsweave logs <job_id> [--tail 50] [--since-event <seq>] [--events-only] [--json]`：按事件序号增量拉取 `events.jsonl` 或日志尾部（避免 Agent 重复读取整份大日志）。
- `omicsweave cancel <job_id> [--json]` 请求协作取消，终态经 `status` 确认。
- **与 `retention.py` 联动**。活跃作业保护对应实验及可能复用的源包。没有活跃作业的交付继续采用现有保留策略。作业记录继续留存，供追查异常。

---

## 6. 安全数据接入、事务配置变更与零污染试验

### 6.1 数据输入探查（`omicsweave inspect-input`）

Agent 接手新数据或新实验时，最易出错的是编码、列名与样本 ID 映射。新增只读探查命令：

```bash
omicsweave inspect-input <path...> [--experiment <exp_id>] [--omics <omics>] [--json]
```

输出结构化摘要：
- 文件存在性、大小、SHA-256（可选）、编码检测（UTF-8 / UTF-8-BOM / UTF-16 / GBK 告警）、分隔符识别；
- 格式推断（`counts_matrix` / `kallisto` / `diann_report` / `protein_matrix` / `peak_table` / `omicsweave_mfa_xlsx` / `legacy_mfa_xlsx`）；
- 列头与候选样本 ID 列表、组织名列表、行数与缺失率概览；
- 若指定了 `--experiment`，自动与现有 `samples.tsv` 对账：列出“已匹配样本”、“数据中有但样本表缺失”、“样本表有但数据缺失”、“建议写入的 `data_column` / `tissue_aliases`”。

### 6.2 事务性配置补丁与样本 QC 更新

1. **`omicsweave config patch --experiment <exp_id> --patch-json '<json>' [--dry-run] [--json]`**：
   - 采用 RFC 7396 JSON Merge Patch（或显式字段路径操作列表）修改 `experiment.yaml` / `project.yaml`。
   - 执行流程：读取原配置 → 内存中应用补丁 → Pydantic 强校验 → 运行 `validate_experiment()` → **调用 `_changed_inputs` / `_analysis_config_changes` 逻辑计算影响面**（明确告知此次改动会使单组学包过期还是仅使整合包过期）→ 若为 `--dry-run` 返回统一 diff 与校验报告；非 `--dry-run` 且无 error 时原子替换文件（保留 `.bak` 或记录审计事件）。
2. **`omicsweave qc flag --experiment <exp_id> --sample <id> --flag <flag|exclude|clear> --reason "<text>" [--dry-run] [--json]`**：
   - 安全更新 `samples.tsv` 的 `qc_flag` / `qc_reason` 列（保持原列顺序、UTF-8 编码、不破坏前导零 ID）；
   - 预计算变更后的 `contrast_reliability` 影响（例如提示“排除该样本后 `liver / HFD_vs_LFD` 将从 `downweighted` 变为 `qc_insufficient`”），并提醒可直接 `--from integrate` 续跑。

### 6.3 零污染参数试验（`--override`）

允许 Agent 在做敏感性分析或方法探索时**不修改主 `experiment.yaml`**：

```bash
omicsweave pipeline --experiment all_tissues \
  --from integrate --source run_20261118_172720_cc0e \
  --override integration.parameters.reaction_dossier.rules.rna_single_source_changed=true \
  --intent sensitivity_check \
  --rationale "Check sensitivity of L1/L2 count to single-source RNA rule" \
  --tag exp_rna_rule --json
```

- 解析规则：`--override key.path=yaml_value`（支持多次传入或 `--override-file <yaml/json>`）。
- 覆盖后的配置同样经 Pydantic 与参数注册表严格校验。
- 生效配置写入结果包 `metadata/experiment.json`。覆盖项、原配置指纹和对照运行保存在 `metadata/lineage.json`，其文件哈希由 v0.1 清单登记。
- 复用源包时校验分析设置、输入与重要产物。 `diff` 对照归档配置及结果变化，保留两侧来源。

---

## 7. 面向 Agent 上下文的金字塔式结果读取层

LLM 上下文窗口有限且按 token 计费，直接倾倒数万行 Parquet 或大型 Excel 既昂贵又容易引发幻觉。设计四级渐进式读取接口（全部封装在 `omicsweave.service.query`，同时提供 CLI 与 MCP）：

| 层级 | 命令 / 接口 | 数据量级 | 核心用途 |
|---|---|---|---|
| **L0 摘要** | `omicsweave inspect --experiment <id> [--run <id>\|--delivery <id>] --json` | ~2–5 KB JSON | 一眼掌握全局状态：各组学/组织/比较的检验数与入选数、`contrast_reliability` 异常、L1 数 vs `decoy_L1_expected`、通路方向冲突数、验证警告、所有 Parquet 表名与行数、图件清单。 |
| **L1 切片** | `omicsweave query --experiment <id> --table <name> [filters] [--columns ...] [--sort ...] [--limit 50] [--offset 0] --json` | 严格分页（默认 50 行，上限 200 行） | 按组织、比较、组学、`padj`、`flux_claim`、`support` 等白名单过滤条件切片读取任意标准表或候选表。 |
| **L2 溯源** | `omicsweave trace --experiment <id> --entity <symbol\|id\|flux_id> [--run <id>] --json` | ~3–10 KB JSON | 跨 MFA、蛋白、转录、代谢、派生特征、反应档案四轴和富集通路，聚合单个实体的全部证据、样本量、QC 标记与冲突说明。 |
| **L3 视觉与对照** | `omicsweave diff <run_a> <run_b> --json` & 图件 PNG 清单 | 结构化差异 + PNG 路径 | 量化对比两次迭代变化（含 `config_diff`、列差值分位数、等级转移矩阵）；提供 PNG 路径供多模态 Agent 直接审图并迭代调用 `map` / `abstract`。 |

**安全查询实现要点（`query`）**：
- **禁止 `eval` / 任意 pandas `df.query()` 表达式注入**。
- 过滤条件采用结构化条件列表：`[{"column": "tissue", "op": "eq", "value": "liver"}, {"column": "padj", "op": "lte", "value": 0.05}]`，支持白名单运算符：`eq / ne / gt / gte / lt / lte / in / contains / is_null / not_null`。
- 底层使用 `pyarrow.dataset` 谓词下推与列裁剪读取，仅反序列化命中的页；`NaN` / `Inf` 严格序列化为 JSON `null` 并在元信息中标注数值健康度。

---

## 8. 带科学护栏的迭代优化闭环与运行谱系

### 8.1 结构化诊断引擎（`core/diagnostics.py`）

新增 `omicsweave diagnose --experiment <id> [--run <id>] --json`，整合 5 处已有质量信号：
1. `validate_experiment()` 的静态与效能警告；
2. `omics/*/qc/sample_quality.parquet` 的样本级离群检测；
3. `qc/contrast_reliability.parquet` 的组学 × 组织 × 比较混杂/样本不足判定；
4. `core/refs.audit()` 的 ID 映射缺失与覆盖率；
5. `integration/reaction_dossier/evaluation/assessment_calibration.parquet` 的诱饵期望（`decoy_L1_expected`）与不变性检查。

每条诊断结论必须归入以下三类策略标签之一，明确告知 Agent **能做什么、不能做什么**：

| 诊断类别（`category`） | 含义与典型项 | 允许的 Agent 动作（`allowed_actions`） |
|---|---|---|
| `engineering_fixable` | 数据接入或展示问题：未映射的组织别名、表编码告警、图件区域拥挤、缺少子图、未登记的 `panels.yaml selection` | **可建议并自动执行（经确认或 `dry_run` 预览）**：补 `tissue_aliases`、重画 `map`/`abstract`、补录元数据 |
| `sensitivity_explorable` | 方法学稳健性检查项：存在高 Z 值样本、批次协变量未纳入、开启/关闭 `adaptive_shrinkage` 或 `flux_groups`、多示踪剂异质性 | **允许通过 `--override` 或 `qc_flag` 做对照试验**，并强制配合 `diff` 汇报结果变化与稳健性，不得单挑有利版本隐瞒基线 |
| `statistical_guardrail` | 统计学硬边界：`min_achievable_pvalue >= alpha`、`qc_confounded`、`qc_insufficient`、`cell_L1` 与 `decoy_L1_expected` 无显著区分 | **禁止通过放宽阈值或更换检验方法来制造显著性**；必须在结论中如实声明为探索性/证据不足 |

### 8.2 科研防过拟合硬护栏（Scientific Guardrails）

1. **面板事后选择强制标记**。已有正式或试验结果时，`config patch` 的预览和提交检查 `selection.mode=post_hoc` 和引用同实验运行 ID 的 `basis`。
2. **探索性试验与正式交付隔离**。覆盖试验及明确登记的敏感性试验写入 `_trials/<trial_id>/`。显式试验查询使用该目录内的交付指针。核定方法后，通过受控配置变更与一次正式计算产生正式交付。
3. **结论口径强提示**：`inspect` / `diagnose` 输出显式携带每个（组织, 比较）的 `cell_L1`、`decoy_L1_expected`、`qc_status` 与 `n_numerator/n_denominator`，并在元数据中携带 `interpretation_caveats` 列表，防止 Agent 解读报告时遗漏关键前提。

### 8.3 独立谱系文件与运行清单兼容

在 `metadata/lineage.json` 记录现行谱系。下例列出主要字段。实际记录包含试验 ID、实验 ID 和基线指纹。 [运行清单](../../../src/omicsweave/contract/manifest.py) · [试验谱系](../../../src/omicsweave/core/overrides.py)

```python
class RunLineage(BaseModel):
    model_config = ConfigDict(extra="forbid")

    parent_run_id: str | None = None
    job_id: str | None = None
    intent: Literal[
        "baseline",
        "qc_refinement",
        "sensitivity_check",
        "visualization_tune",
        "post_hoc_exploration",
        "manual",
    ] = "manual"
    rationale: str = ""
    overrides: dict[str, Any] = Field(default_factory=dict)
```

配合 `omicsweave diff <run_a> <run_b> --json`，自动附加 `config_diff`（对比两包 `metadata/experiment.json`）与 `lineage` 信息，使 Agent 的每一轮优化形成完整、可审计的科学实验记录。

### 8.4 随仓分发的领域操作手册（Agent Skill）

在仓库根目录或文档中提供面向 Agent 的标准技能文件（如 `skills/omicsweave-analyst/SKILL.md`），内置：
- 标准操作 SOP：`inspect-input` → `validate` → `pipeline --dry-run` → `submit` → `status` → `inspect/diagnose` → `query/trace` → `diff`；
- 统计与生物学红线（MFA 通量 ≠ 代谢物丰度、推断单位是动物而非拟合迭代、不可用不显著推出“无变化”、诱饵校准接近时不得夸大靶点结论）；
- 上下文节约准则（永远先读 L0 `inspect`，再用 `query` 分页，禁止全表转文本）。
