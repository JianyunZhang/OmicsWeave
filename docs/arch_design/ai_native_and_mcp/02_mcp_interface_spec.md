# OmicsWeave MCP（Model Context Protocol）接口详细设计规范

> **文档类型**：协议与接口规范（v0.1，2026-11-21）\
> **所属套件**：[AI Native 与 MCP 设计套件](README.md)\
> **前置架构**：[01_ai_native_architecture.md](01_ai_native_architecture.md) · **实施排期**：[03_implementation_roadmap.md](03_implementation_roadmap.md)\
> **协议基线**：Model Context Protocol（Python 官方 SDK `mcp` / `FastMCP`，支持 `stdio` 与 `Streamable HTTP`）

**本机接口**。使用 `omicsweave --project-dir <项目> mcp serve` 或 `omicsweave-mcp --project-dir <项目>` 启动固定版本 SDK 的 stdio 服务。服务器批准的 `--allowed-data-root` 用于读取。配置与样本表提交需要 `--allow-writes` 和预览 SHA-256。项目脚本默认跳过，显式放行后仅接受项目内 Python 脚本。清理工具提供预览，旧版迁移与实际删除走本机 CLI。

**Phase 3.1／3.2**。提交工具接受 `overrides`、`intent`、`rationale` 和 `parent_run_id`。计算完成后，用返回的 `trial_id` 查询或诊断，以 `diff_runs(trial_b=...)` 与正式基线进行对照。正式发布走已审阅配置的正式运行。网络传输与系统间互操作属于 3.3。

---

## 1. 架构定位、依赖隔离与 Service 共享

### 1.1 Thin MCP Facade 原则

`omicsweave.mcp` 定位为**薄协议适配层**，本身不包含生物信息统计逻辑，也不通过拼接 CLI 字符串和正则解析 stdout 来驱动 OmicsWeave：

- 所有业务逻辑、DTO 定义、路径校验与错误分类下沉至 `omicsweave.service`；
- CLI（`--json`）与 MCP Server 均调用同一个 `omicsweave.service` 函数；
- 长耗时任务由 `omicsweave.service.execution` 委托 `omicsweave.core.jobs` 以**独立子进程 Worker** 执行，确保 MCP Server 进程本身常驻轻量、低内存占用、无旧模块全局状态污染，且不会因某次分析崩溃导致 MCP 连接断开。

### 1.2 安装与启动入口

在 [pyproject.toml](../../../pyproject.toml) 中新增独立可选依赖组与脚本入口，不增加基础用户的安装负担：

```toml
[project.optional-dependencies]
mcp = [
    "mcp==1.26.0",
]

[project.scripts]
omicsweave = "omicsweave.cli:main"
omicsweave-mcp = "omicsweave.mcp.server:main"
```

CLI 同步提供子命令入口：

```bash
# 本地默认 stdio 模式（推荐，供桌面/IDE Agent 直接拉起）
omicsweave --project-dir /path/to/project mcp serve

# 独立入口使用同一 stdio 服务
omicsweave-mcp --project-dir /path/to/project
```

---

## 2. 传输模式、安全沙箱与权限门控

### 2.1 传输模式策略

| 模式 | 默认状态 | 适用场景 | 安全基线 |
|---|---|---|---|
| `stdio` | **默认启用** | Claude Desktop、Claude Code、Cursor、Antigravity 等本机 Agent | 服务端通过父进程标准输入输出通信。需另行审查客户端模型调用、遥测和日志转发，以及课题配置中的远程基因集访问。 |
| `streamable-http` | 显式开启 | O5 阶段实验室工作站/HPC 集群远程调用、多客户端协同 | 默认仅绑定 `127.0.0.1`；跨主机绑定时强制要求 Bearer Token（通过环境变量 `OMICSWEAVE_MCP_TOKEN` 注入），建议配合 SSH 隧道或内网 TLS 反代。 |

### 2.2 路径沙箱与越权拦截（`omicsweave.mcp.security`）

MCP Server 启动时通过 `--project-dir` 锁定单一项目根目录 `PROJECT_ROOT`（或显式允许的 `--allowed-root` 只读外部数据目录列表，用于读取位于项目外的大体积原始数据文件）：

1. **写入路径硬约束**：所有写操作（配置补丁、`qc_flag`、导出、图件重画、结果输出、作业状态）的目标路径必须经过 `Path.resolve().is_relative_to(PROJECT_ROOT)` 校验，拒绝 `..` 穿越与指向项目外的符号链接。
2. **读取路径白名单**：`inspect_input_files` 只允许访问 `PROJECT_ROOT` 及启动参数显式声明的 `--allowed-data-root` 列表。
3. **自定义脚本执行门控**：`project.yaml` 的 `pipeline.scripts` 可执行任意本地脚本或外部命令。MCP 模式下默认 `--allow-custom-scripts=false`（执行 `submit_pipeline` 时跳过自定义脚本并提示，或仅执行 `validate → run → integrate → summary`）；仅在用户显式传入 `--allow-custom-scripts` 时才允许 MCP 触发课题脚本。
4. **写入门控**。`patch_experiment_config` 与 `update_sample_qc` 默认预览。提交需要服务器写入授权、目标文件的预览 SHA-256 和变更理由。`cleanup_deliveries` 提供清理计划，实际删除由本机 `clean` 命令执行。

Phase 2 的权限测试须覆盖 `project.yaml`、`experiment.yaml`、样本表及参考资产间接解析出的路径，以及 Worker 生效的结果目录。工具调用参数只是入口。若开放给互不信任的本机用户，还需部署 OS 权限隔离，并控制可被替换的目录及文件。

### 2.3 MCP Tool Annotations 映射

每个 Tool 均声明标准 MCP `ToolAnnotations`，辅助客户端自动审批或弹出确认：

| 操作级别 | `readOnlyHint` | `destructiveHint` | `idempotentHint` | `openWorldHint` | 对应工具 |
|---|---:|---:|---:|---:|---|
| 纯只读探查/查询/校验/Diff | `true` | `false` | `true` | `false` | `list_experiments`、`validate_experiment`、`inspect_input_files`、`get_job_status`、`get_run_digest`、`query_table`、`trace_entity`、`diagnose_run`、`diff_runs` |
| 不可变计算提交/出图/导出 | `false` | `false` | `false` | `false` | `submit_pipeline`、`submit_omics_run`、`submit_integration`、`render_figure`、`export_results` |
| 作业控制 | `false` | `false` | `true` | `false` | `cancel_job` |
| 工作区配置/样本表写入 | `false` | `false` | `false` | `false` | `init_or_migrate_experiment`、`patch_experiment_config`、`update_sample_qc` |
| 不可逆磁盘清理 | `false` | `true` | `false` | `false` | `cleanup_deliveries` |

Tool 注解表达实际副作用，访问控制由服务端执行。允许调用外部数据源的工具应如实标注。离线模式应拦截未获准的远程基因集与自定义程序。

---

## 3. MCP Resources 规范（只读上下文寻址）

Resources 用于提供稳定、可缓存、只读的结构化上下文，避免 Agent 反复通过 Tool 拉取大段静态配置或 Schema。所有资源均返回 `application/json`。

| URI 模板 | 说明 | 核心返回字段 |
|---|---|---|
| `omicsweave://project/overview` | 项目全局概览 | `project_id`、`title`、`species`、`paths`、`experiments[]`、`pipeline_scripts[]`、`registries_status`、`active_jobs_count` |
| `omicsweave://schema/{target}` | 参数与契约 Schema 自省 | `target` 支持 `project`、`experiment`、`pathway_map`、`manifest`、`omics/{name}`、`integration/{mode}` 与 `tables/{name}`。参数 Schema 带 `x-tunable-class`。 |
| `omicsweave://experiments/{exp_id}/config` | 展开后的实验设计视图 | `experiment_config`、`subject_count`、`sample_count` 和按组学及组织汇总的 `sample_matrix`。 |
| `omicsweave://experiments/{exp_id}/latest` | 最新正式交付指针与索引 | 等价 `latest --json` + `index.json` 阶段状态与主要产物相对路径 |
| `omicsweave://runs/{exp_id}/{run_id}/manifest` | 不可变结果包清单 | `RunManifest`（含 `lineage`、`parameters`、`files[]`、输入/输出 SHA-256） |
| `omicsweave://runs/{exp_id}/{run_id}/digest` | L0 运行诊断摘要 | 同 `get_run_digest` 返回内容，便于客户端直接订阅/挂载到对话上下文 |
| `omicsweave://jobs/{job_id}/state` | 作业实时状态 | 经 PID 实时探活修正后的 `state.json` |

---

## 4. MCP Tools 详细契约

共计 **17 个高内聚工具**，分为 5 组，覆盖“发现校验 → 数据与配置 → 后台执行与状态 → 结果读取与出图 → 诊断与迭代对比”全闭环。

> **通用约定**：所有工具均返回结构化 JSON（对应 `ServiceResponse` 语义），发生异常时使用 MCP `isError=true` 并附带统一的 `{code, message, hint, retryable, details}` 错误结构；且不泄露宿主机项目根目录之外的敏感路径。

### 4.1 A 组：项目发现与校验（3 个）

#### 1. `list_experiments`
- **用途**：列出当前项目内全部实验、启用组学、最近单组学/整合 `run_id`、最新交付指针及运行中作业。
- **输入**：`{ include_active_jobs: bool = true }`
- **输出**：`{ project_id, species, experiments: [{ experiment_id, enabled_omics, integration_modes, latest_delivery_id, latest_run_id, latest_integration_run_id, active_job_ids }] }`

#### 2. `validate_experiment`
- **用途**：执行配置、样本表、输入文件、最小可达 p 值与整合前提校验（对应 `validate_experiment` + 参数注册表静态校验）。
- **输入**：
  ```json
  {
    "experiment_id": "all_tissues",
    "include_registry_audit": false
  }
  ```
- **输出**：`{ ok: bool, error_count: int, warning_count: int, issues: [{ level, code, where, message, hint }], contrast_power_summary: [...] }`

#### 3. `inspect_input_files`
- **用途**：在写入配置或修改样本表前探查原始/表格数据文件，避免编码与列名错配。
- **输入**：
  ```json
  {
    "paths": ["experiments/exp01/data/rna_counts.tsv"],
    "experiment_id": "exp01",
    "omics_hint": "transcriptomics"
  }
  ```
- **输出**：`{ files: [{ path, exists, size_bytes, encoding, encoding_ok, detected_format, columns_preview, detected_samples, detected_tissues, row_count_estimate }], alignment_with_samples_tsv: { matched, only_in_data, only_in_samples_tsv, suggested_tissue_aliases } }`

---

### 4.2 B 组：项目初始化与受控配置变更（3 个）

#### 4. `init_or_migrate_experiment`
- **用途**：新建实验模板（`init`）或将旧版各组学配置无损转换为统一模型（`migrate`）。
- **输入**：`{ action: "init" | "migrate", experiment_id: str, species?: "mus_musculus" | "homo_sapiens", dry_run: bool = true, force: bool = false }`
- **输出**：`{ action, dry_run, files_planned_or_written: [str], migration_todos?: [str], validation_ok: bool }`

#### 5. `patch_experiment_config`
- **用途**：事务性修改 `experiment.yaml`（或指定实验的 `panels.yaml`），强制经过 Pydantic + 校验 + 过期影响评估。
- **输入**：
  ```json
  {
    "experiment_id": "all_tissues",
    "target": "experiment_yaml",
    "patch": {
      "design": {
        "tissue_aliases": { "white_adipose": "sq_adipose" }
      }
    },
    "dry_run": true,
    "reason": "Align adipose naming between MFA and transcriptomics"
  }
  ```
- **输出**：
  ```json
  {
    "applied": false,
    "dry_run": true,
    "yaml_diff": "--- ...",
    "validation": { "ok": true, "errors": [], "warnings": [] },
    "invalidation_impact": {
      "invalidates_omics_bundles": true,
      "invalidates_integration_bundles": true,
      "recommended_rerun_from": "run"
    },
    "guardrail_notices": []
  }
  ```
- **护栏逻辑**：当 `target == "panels_yaml"` 且该实验已存在结果包时，若补丁未声明 `selection.mode` 或将其设为 `pre_specified`，自动拦截并返回 `guardrail_notices` 要求显式标为 `post_hoc` 并登记 `basis`。

#### 6. `update_sample_qc`
- **用途**：批量设置或清除 `samples.tsv` 的 `qc_flag`（`flag | exclude | ""`）和 `qc_reason`。
- **输入**：
  ```json
  {
    "experiment_id": "all_tissues",
    "updates": [
      { "sample_id": "M04_liver_rna", "qc_flag": "flag", "qc_reason": "Library composition outlier (robust_z=4.8)" }
    ],
    "dry_run": true
  }
  ```
- **输出**：`{ applied, dry_run, changed_rows, projected_reliability_changes: [{ omics, tissue, contrast_id, before_status, after_status }], recommended_rerun_from: "integrate" }`

---

### 4.3 C 组：执行与异步作业控制（5 个）

#### 7. `submit_pipeline`
- **用途**：运行 `validate → run → integrate → [scripts] → summary` 标准流程，支持增量续跑、零污染 `overrides`、标签与锁定。
- **输入**：
  ```json
  {
    "experiment_id": "all_tissues",
    "start": "integrate",
    "source_run_id": "run_20261118_172720_cc0e",
    "integration_run_id": null,
    "omics": ["all"],
    "overrides": {"integration.parameters.reversal.equivalence_margin": 0.7},
    "tag": null,
    "pin": false,
    "dry_run": true,
    "intent": "sensitivity_check",
    "rationale": "比较预先登记的等效性边界",
    "parent_run_id": null,
    "wait_seconds": 15
  }
  ```
- **本机实现**。试验独立存放，并在响应中返回 `trial_id` 和谱系。 `dry_run=true` 同步预检，正式提交后通过作业 ID 观察。试验复用正式结果时须声明确切来源，正式交付通过独立正式运行发布。

#### 8. `submit_omics_run`
- **用途**：仅运行指定组学与步骤（对应 `omicsweave run`）。
- **输入**：`{ experiment_id, omics: ["all"], steps?: ["qc", "differential", "enrichment", "figures"], overrides?: {}, intent?: str, rationale?: str, parent_run_id?: str, wait_seconds: int = 15 }`

#### 9. `submit_integration`
- **用途**：基于已有单组学结果包运行指定整合模式（对应 `omicsweave integrate`）。
- **输入**：`{ experiment_id, source_run_id?: str, modes?: [str], overrides?: {}, intent?: str, rationale?: str, parent_run_id?: str, wait_seconds: int = 15 }`

#### 10. `get_job_status`
- **用途**：查询作业实时状态、经 PID 探活后的真实存活状态、细粒度进度、最新结构化事件与日志尾部。
- **输入**：`{ job_id: str, since_event_seq: int = 0, tail_log_lines: int = 30 }`
- **输出**：
  ```json
  {
    "job_id": "job_20261121_183012_9f2a",
    "status": "running",
    "alive": true,
    "heartbeat_age_seconds": 3.2,
    "progress": { "stage": "run", "step": "omics=proteomics | tissue=liver | step=differential", "fraction": 0.55 },
    "artifacts": { "run_id": "run_...", "integration_run_id": null, "delivery_id": "delivery_..." },
    "new_events": [ ... ],
    "next_event_seq": 18,
    "log_tail": [ ... ],
    "poll_after_seconds": 10
  }
  ```

#### 11. `cancel_job`
- **用途**：向运行中作业发送取消请求并终止子进程 Worker。
- **输入**：`{ job_id: str, force: bool = false }`
- **输出**：`{ job_id, previous_status, current_status, terminated_pid }`

---

### 4.4 D 组：结果读取、切片查询、可视化与运维（5 项）

#### 12. `get_run_digest`
- **用途**：L0 全局摘要——以极低 token 开销让 Agent 掌握一次运行或交付的全部关键质量与发现信号。
- **输入**：`{ experiment_id: str, run_id: str = "latest", delivery_id?: str }`
- **输出**：
  ```json
  {
    "experiment_id": "all_tissues",
    "run_id": "run_20261118_172720_cc0e",
    "integration_run_id": "run_20261118_195428_674e",
    "delivery_id": "delivery_20261118_195501_a1b2",
    "lineage": { "intent": "baseline", "parent_run_id": null, "overrides": {} },
    "omics_summary": [
      { "omics": "mfa", "tissue": "liver", "contrast_id": "HFD_vs_LFD", "n_tested": 73, "n_selected": 12, "reliability_status": "ok" }
    ],
    "reliability_alerts": [],
    "candidate_overview": {
      "reactions_L3": 0, "reactions_L2": 72, "reactions_L1": 19,
      "decoy_calibration": [{ "tissue": "liver", "contrast_id": "HFD_vs_LFD", "cell_L1": 19, "decoy_L1_expected": 15.0 }],
      "pathways_fused": 58, "pathways_direction_conflict": 15
    },
    "tables": [{ "name": "omics/mfa/tables/differential.parquet", "rows": 450, "columns": ["..."] }],
    "figures": [{ "kind": "pathway_map", "tissue": "liver", "contrast_id": "HFD_vs_LFD", "png_path": "...", "pdf_path": "..." }],
    "interpretation_caveats": [
      "肝脏 HFD_vs_LFD 的 L1 数（19）与诱饵期望（15.0）接近，宜视为探索性排序而非已确认信号。"
    ]
  }
  ```

#### 13. `query_table`
- **用途**：L1 安全分页切片查询结果包或候选目录中的任意 Parquet/TSV 表。
- **输入**：
  ```json
  {
    "experiment_id": "all_tissues",
    "run_id": "latest",
    "table": "reaction_assessment",
    "omics": null,
    "mode": "reaction_dossier",
    "filters": [
      { "column": "tissue", "op": "eq", "value": "liver" },
      { "column": "flux_claim", "op": "in", "value": ["L1", "L2", "L3"] }
    ],
    "columns": ["flux_id", "tissue", "contrast_id", "flux_claim", "direction", "mechanism_class", "rank"],
    "sort_by": "rank",
    "ascending": true,
    "limit": 30,
    "offset": 0
  }
  ```
- **防上下文溢出约束**：`limit` 硬上限为 `200`（默认 `50`），单次响应序列化字节硬上限建议 `128 KB`，超限时截断并提示进一步加过滤条件或裁剪列。
- **输出**：`{ run_id, table_path, total_rows_before_filter, total_rows_matched, offset, limit, truncated, columns, rows: [...] }`

#### 14. `trace_entity`
- **用途**：L2 跨层单实体全景溯源，将 MFA、转录组、蛋白组、代谢组、反应档案成员/四轴状态及命中通路聚合返回。
- **输入**：`{ experiment_id: str, query: str, entity_type?: "gene" | "protein" | "metabolite" | "reaction" | "auto", run_id: str = "latest", tissue?: str, contrast_id?: str }`
- **输出**：`{ resolved_entities, differential_records, derived_records, reaction_dossier_axes, reaction_members, enrichment_hits, reliability_context }`

#### 15. `render_figure`
- **用途**：基于已有整合结果包迭代重画分组织通路总图（`map`）、子通路图（`submaps`）或 Graphical Abstract（`abstract`），不改动原结果包。
- **输入**：
  ```json
  {
    "experiment_id": "all_tissues",
    "kind": "pathway_map",
    "source_run_id": null,
    "tissues": ["liver"],
    "contrasts": ["HFD_vs_LFD"],
    "mode": "publication",
    "style": "metro_atlas",
    "color_mode": "semantic",
    "strip_mode": "omics",
    "regions": null,
    "submaps": ["central_carbon"],
    "formats": ["png", "pdf"],
    "return_preview_image": false
  }
  ```
- **多模态审图支持**：当 `return_preview_image=true` 且生成的 PNG 不超过 1 张时，除返回结构化路径清单外，可按 MCP 规范附带缩放后的 `ImageContent`（如最大边 1280px PNG），使支持视觉的 Agent 直接检查标签重叠、布局密度与配色效果并继续迭代参数。

#### 16. `export_or_clean`（或拆为 `export_results` / `cleanup_deliveries`）
- 推荐保持语义清晰，独立为 `export_results`（导出 XLSX/TSV）和 `cleanup_deliveries`（`dry_run=true` 默认预览，按 `RetentionPolicy` 清理过期交付），合计 **17–18 个工具**。

---

### 4.5 E 组：诊断与迭代对比（2 个）

#### 17. `diagnose_run`
- **用途**：汇总校验、样本质量、比较可靠性、知识表审计与诱饵校准，输出按科学护栏分级的诊断报告。
- **输入**：`{ experiment_id: str, run_id: str = "latest", include_registry_audit: bool = true }`
- **输出**：
  ```json
  {
    "experiment_id": "all_tissues",
    "run_id": "run_...",
    "overall_health": "caution",
    "findings": [
      {
        "id": "DIAG-001",
        "category": "engineering_fixable",
        "severity": "medium",
        "code": "UNMAPPED_TISSUE_ALIAS",
        "summary": "MFA 中的 white_adipose 未与转录组 sq_adipose 对齐，导致该组织组学层全为 not_measured",
        "evidence": { "mfa_tissue": "white_adipose", "other_tissues": ["sq_adipose"] },
        "suggested_tool": "patch_experiment_config",
        "suggested_args_preview": { "patch": { "design": { "tissue_aliases": { "white_adipose": "sq_adipose" } } } }
      },
      {
        "id": "DIAG-002",
        "category": "statistical_guardrail",
        "severity": "high",
        "code": "DECOY_L1_INDISTINGUISHABLE",
        "summary": "肝脏 HFD_vs_LFD 的 L1 数（19）与诱饵期望（15.0 ± 3.2）无明显分离",
        "evidence": { "cell_L1": 19, "decoy_L1_expected": 15.0 },
        "optimization_blocked": true,
        "guidance": "禁止通过放宽阈值或事后筛选子集提高 L1 计数；报告中需标明探索性。"
      }
    ]
  }
  ```

#### 18. `diff_runs`
- **用途**：结构化逐表对比两个结果包（复用 [contract/diff.py](../../../src/omicsweave/contract/diff.py)），同时对比两者 `metadata/experiment.json` 的配置差异。
- **输入**：`{ run_a: str, run_b: str, experiment_id?: str, atol: float = 1e-12, rtol: float = 1e-9, top_changed_columns: int = 20 }`
- **输出**：`{ identical: bool, config_diff: [...], table_summary: [{ path, status, rows_a, rows_b, cells_differ, max_abs_diff }], column_diffs: [...], grade_transitions: [{ path, column, from, to, count }], other_files: {...} }`

---

## 5. 长任务执行协议（Hybrid Sync/Async + Progress）

MCP 调用方对单个 Tool 调用通常有 30–120 秒的超时限制，而真实多组学分析可能运行数十分钟。采用**混合同步/异步模式**彻底消除超时风险：

```text
Agent Client                         OmicsWeave MCP Server                  Worker Subprocess (_jobs/<id>/)
     │                                         │                                           │
     │  1. submit_pipeline(wait_seconds=15)    │                                           │
     │────────────────────────────────────────▶│  2. 核对权限与请求参数                     │
     │                                         │  3. 创建 _jobs/<job_id>/ & 拉起子进程     │
     │                                         │──────────────────────────────────────────▶│
     │                                         │                                           │ 写入 pid, state=running
     │  4. [若支持] notifications/progress     │◀───────── 轮询 events.jsonl (每 1s) ──────│ 持续 emit step_progress
     │◀════════════════════════════════════════│                                           │ & 刷新 heartbeat_at
     │                                         │                                           │
     │  5a. 若 ≤15s 完成: 直接返回 completed    │                                           │
     │  5b. 若 >15s 未完: 立即返回 running +    │                                           │
     │      job_id + poll_after_seconds        │                                           │
     │◀────────────────────────────────────────│                                           │
     │                                         │                                           │ (后台独立继续运行，
     │  6. get_job_status(job_id, since_seq)   │                                           │  即使 MCP 断开也不中断)
     │────────────────────────────────────────▶│  PID 探活 + 读 state.json + 增量 events   │
     │◀────────────────────────────────────────│                                           │
     │                                         │                                           │ 结束写 completed & run_id
     │  7. get_run_digest / diagnose_run       │                                           │
     │────────────────────────────────────────▶│                                           │
```

**关键健壮性保证**：
- 普通 Tool 请求只在收到有效 `progressToken` 且请求仍在进行时发送 `notifications/progress`。返回项目自己的 `job_id` 后，通过 `get_job_status` 与事件游标继续观察。采用 MCP Tasks 时应另行协商能力并测试生命周期。[MCP 进度规范](https://modelcontextprotocol.io/specification/2025-11-25/basic/utilities/progress)
- Phase 2 固定并验证 Python SDK、协议版本与目标客户端的组合。stdio 的 stdout 仅承载 MCP 消息。未来开放 Streamable HTTP 时，应验收 Origin 检查、授权及作用域。[MCP 传输规范](https://modelcontextprotocol.io/specification/2025-11-25/basic/transports) · [MCP Python SDK](https://py.sdk.modelcontextprotocol.io/run/)
- **断线恢复**。Worker 使用独立会话或进程组和文件日志。在客户端或宿主未回收子进程的情况下，重启 Server 后可从磁盘恢复作业查询。需在发版客户端逐项验证进程回收策略。
- **并发写保护**。配置写入需要实验级锁。分析任务在独立的 `<run_id>/` 与 `_jobs/<job_id>/` 计算，共享 `latest.json` 发布和保留策略使用结果根目录下的 OS 锁。

---

## 6. MCP Prompts（内置标准化科研 SOP 模板）

通过 MCP Prompts 把最佳实践编排为可复用提示词模板，供用户一键触发：

| Prompt 名称 | 参数 | 编排步骤与内置约束 |
|---|---|---|
| `onboard_experiment` | `experiment_id`, `data_paths` | `inspect_input_files` → 读取 `omicsweave://schema/experiment` → `patch_experiment_config(dry_run=true)` → 确认后落盘 → `validate_experiment` → `submit_pipeline(dry_run=true)` → 正式基线运行 |
| `qc_review_and_rerun` | `experiment_id`, `run_id="latest"` | `get_run_digest` → `diagnose_run` → `query_table(table="sample_quality")` → 列出异常样本证据并征询用户 → `update_sample_qc(dry_run=true/false)` → `submit_pipeline(start="integrate", intent="qc_refinement")` → `diff_runs` 对比稳健性 |
| `method_sensitivity_study` | `experiment_id`, `source_run_id`, `override_variants` | 保持基线 `experiment.yaml` 不变，逐个以 `overrides` + `intent="sensitivity_check"` 提交整合或后处理 → 全部完成后调用 `diff_runs` 生成跨方案等级转移与稳健性汇总表 |
| `iterate_publication_figures` | `experiment_id`, `tissue`, `contrast_id` | 读取 `get_run_digest` 中活跃通路区域 → 调用 `render_figure(mode="publication", return_preview_image=true)` → 视觉检查拥挤度与标注 → 调整 `regions / submaps / style` 直至出版质量 |

---

## 7. 客户端集成配置与协议测试策略

### 7.1 典型客户端 `stdio` 配置示例

适用于 Claude Desktop、Cursor、Claude Code、Antigravity 等支持标准 MCP 配置的客户端：

```json
{
  "mcpServers": {
    "omicsweave": {
      "command": "C:/Users/zhang/.conda/envs/tutorial-mfa-2025/python.exe",
      "args": [
        "-m", "omicsweave",
        "--project-dir", "C:/Users/zhang/PycharmProjects/OmicsWeave/examples/demo_project",
        "mcp", "serve"
      ],
      "env": {
        "PYTHONIOENCODING": "utf-8"
      }
    }
  }
}
```

### 7.2 MCP 契约测试规范

在 `tests/mcp/` 下编写无外部网络依赖的自动化测试套件（依托 `examples/demo_project` 合成数据集与 MCP SDK 内置的内存传输 `create_connected_server_and_client_session`）：

1. **工具目录与 Schema 完整性测试**：验证所有 Tools、Resources、Prompts 的 Schema 可合法序列化，Tool Annotations 符合风险等级表。
2. **安全沙箱对抗测试**：验证 `../` 路径穿越、绝对路径越界、未开 `--allow-custom-scripts` 时触发脚本被拦截、默认 `dry_run=true` 不改写磁盘文件。
3. **端到端闭环测试**：在临时复制的 `demo_project` 上经 MCP 顺序执行 `list_experiments` → `validate_experiment` → `submit_pipeline(wait_seconds=0)` → 轮询 `get_job_status` 至 `completed` → `get_run_digest` → `query_table` → `diagnose_run` → 带 `overrides` 二次运行 → `diff_runs`。
4. **崩溃与取消恢复测试**：模拟 Worker 异常退出，验证 `get_job_status` 能通过 PID 探活识别 `stalled`；测试 `cancel_job` 能正确终止运行中任务。
5. **流隔离测试**：在 `stdio` 子进程模式下运行真实命令，确认 `stdout` 中绝不混入任何非 JSON-RPC 的中文提示或第三方库日志。
