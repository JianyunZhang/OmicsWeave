# OmicsWeave AI Native 与 MCP 分阶段实施路线、验证与治理

> **文档类型**：实施路线与验收规范（v0.1，2026-11-21）\
> **所属套件**：[AI Native 与 MCP 设计套件](README.md)\
> **参考架构**：[01_ai_native_architecture.md](01_ai_native_architecture.md) · [02_mcp_interface_spec.md](02_mcp_interface_spec.md)

---

## 1. 目标包目录结构与既有文件改动清单

### 1.1 新增与改动文件树

```text
OmicsWeave/
├── pyproject.toml                       # [改] 新增 optional-dependencies.mcp 与 omicsweave-mcp 入口
├── skills/
│   └── omicsweave-analyst/
│       └── SKILL.md                     # [新] 面向外部 Agent 的领域操作手册与科研红线规范
├── src/omicsweave/
│   ├── cli.py                           # [改] 补齐 --format/--json、status/logs/cancel、inspect/query/trace、
│   │                                    #      diagnose、schema、inspect-input、config patch、qc flag、mcp serve
│   ├── _worker.py                       # [新] 独立子进程 Worker 入口（python -m omicsweave._worker <job_dir>）
│   ├── contract/
│   │   ├── manifest.py                  # [改] RunManifest 增加向后兼容的 lineage 可选字段
│   │   └── diff.py                      # [改] 增加结构化 dict/JSON 导出与 metadata/experiment.json 配置对照
│   ├── core/
│   │   ├── progress.py                  # [新] ProgressReporter 协议、NullReporter、JsonlFileReporter
│   │   ├── jobs.py                      # [新] JobSpec/JobState 状态机、子进程派发、心跳线程、跨平台 PID 探活、取消
│   │   ├── param_schemas.py             # [新] 各组学 analysis 与整合 parameters 的强类型 Pydantic Schema 注册表
│   │   ├── diagnostics.py               # [新] 聚合 validate/reliability/sample_quality/refs/calibration 诊断引擎
│   │   ├── overrides.py                 # [新] --override 路径解析、校验、内存叠加与 config_diff 计算
│   │   ├── run.py                       # [改] RunContext 增加 reporter 与 lineage 字段；在组学循环发进度事件
│   │   ├── pipeline.py                  # [改] 接入 ProgressReporter、overrides、promote_latest 与 lineage
│   │   ├── validate.py                  # [改] Issue 增加结构化 code 字段；静态校验分析/整合参数注册表
│   │   └── retention.py                 # [改] 结合 _jobs/ 探活判定活跃交付，可选清理过期历史作业目录
│   ├── service/                         # [新] 无状态应用服务层（CLI 与 MCP 共享，Pydantic Request/Response DTO）
│   │   ├── __init__.py
│   │   ├── schemas.py                   # ServiceResponse、ServiceError、统一 ErrorCode 枚举
│   │   ├── project.py                   # list_experiments、validate、inspect_input、init/migrate、config patch、qc flag
│   │   ├── execution.py                 # submit_pipeline/run/integrate、get_job_status、get_logs、cancel_job
│   │   ├── query.py                     # get_run_digest、query_table、trace_entity、diff_runs、diagnose_run
│   │   └── visualization.py             # render_pathway_map、render_abstract（封装现有 redraw/render_abstracts）
│   └── mcp/                             # [新] MCP Server（仅在安装 omicsweave[mcp] 时加载）
│       ├── __init__.py
│       ├── security.py                  # 项目根路径沙箱校验、文件锁、脚本执行门控
│       ├── resources.py                 # omicsweave:// URI 资源路由
│       ├── tools.py                     # 17-18 个核心工具注册、风险注解与 Hybrid Sync/Async 进度转发
│       ├── prompts.py                   # 4 个科研 SOP Prompt 模板
│       └── server.py                    # FastMCP 装配与 stdio / streamable-http 启动入口
└── tests/
    ├── core/
    │   ├── test_jobs.py                 # [新] 作业状态机、后台子进程、PID 探活、心跳、取消与崩溃恢复测试
    │   ├── test_param_schemas.py        # [新] 注册表与各 runner/mode DEFAULTS 一致性守卫测试
    │   ├── test_overrides.py            # [新] --override 解析、校验、零污染与 config_diff 测试
    │   └── test_diagnostics.py          # [新] 诊断分类与科研防过拟合护栏测试
    ├── service/
    │   ├── test_query_service.py        # [新] L0 Digest、安全切片 query、trace_entity、防注入与分页限流测试
    │   └── test_cli_json.py             # [新] 全部核心命令 --json 信封、流隔离与退出码契约测试
    └── mcp/
        ├── test_mcp_protocol.py         # [新] In-Memory MCP 全链路闭环测试（基于 demo_project）
        └── test_mcp_security.py         # [新] 路径穿越、默认 dry_run、脚本门控与 stdio 洁净度测试
```

### 1.2 既有模块最小侵入守则

- **不改动四个组学统计与整合算法的核心计算函数**；仅在 `mfa/metabolomics/proteomics/transcriptomics/runner.py` 与 `integration/modes/*.py` 的步骤与组织外层循环插入 `ctx.reporter.emit(...)` 和 `ctx.reporter.check_cancelled()`。
- 通过单元测试 `test_param_schemas.py` 强制断言：`param_schemas.py` 中声明的键集合和默认值与各 Runner/Mode 内的 `DEFAULTS` 严格一致，避免后续新增参数后 Schema 漂移。

---

## 2. 分阶段实施路线（Phase 1 → Phase 3）

本路线作为 OmicsWeave 总纲（[omicsweave_theraweave_development_plan.md](../../omicsweave_theraweave_development_plan.md)）中 O2 与 O5 之间的正交能力增强线（记为 **A1–A3**），**不阻塞现有课题分析与 O2/O3 推进**，每阶段均可独立交付产生直接收益。

```text
Phase 1 (A1): 可观测后台 CLI + 机器协议     ──▶ 任何 Coding Agent 通过 Shell 即可稳定闭环操作
       │
       ▼
Phase 2 (A2): 强 Schema + Service + MCP v1 ──▶ 任意 MCP 客户端即插即用，达成单机 AI Native
       │
       ▼
Phase 3 (A3): 零污染试验 + 科学护栏 + 远程化 ──▶ 严谨科研迭代闭环，无缝对接 O5 集群与 TheraWeave MCP
```

### Phase 1（A1）— 可观测后台 CLI 与结构化读写底座

**目标**：先解决 Agent 通过后台 CLI 使用软件的 80% 痛点——异步执行、状态可观测、JSON 输出与分页读结果。

- [x] **1.1 进度与作业运行时**。独立 Worker、原子状态、心跳、进程身份、协作取消，以及流水线、四组学 Runner 和整合入口的进度检查点。
- [x] **1.2 统一 CLI 机器协议**。九个既有命令支持版本化信封与结构化失败，日志进入 stderr 或作业日志。沿用 `diff` 的 0/1 退出语义表示相同/有差异，机器模式的失败用语义退出码。
- [x] **1.3 运行态与查询子命令**。支持 `status`、`logs`、`cancel`、`inspect`、`query`、`trace`，以及 `run / integrate / pipeline --background`。首版实体溯源按表内确切标识符和富集成员匹配，保留表路径及证据字段。
- [x] **1.4 交付与保留策略联动**。活跃作业保护相关实验及潜在源结果包。派发与清理使用同一结果根目录的 OS 锁。共享 `latest.json` 使用独立 OS 锁保护读取、合并及原子替换。

**工作树复核后的实施决定**

| 复核发现 | 实施决定 |
|---|---|
| `latest --json` 已有对象/数组调用契约，图件和导出命令已有 `--format`。 | 保留旧调用，使用 `latest --format json` 取得新信封，图件命令可用 `--json --format png`。 |
| 现行交付与结果包可以同目录保存，重跑阶段也可引用其他包。 | 读取实际交付索引和 `latest.json`，沿清单的输出记录与 `source_run` 查表。 |
| PID 可复用，心跳延迟也可能来自计算压力。脚本能派生子进程。 | 用 PID 和进程身份探活。心跳迟到显示 `unresponsive`。取消采用协作检查点，强制终止留到进程树管理验证后。 |
| 完整预检会读取样本表并执行审计。 | 派发时解析项目及结果根、记录作业规格，完整预检在 Worker 执行。单独调用 `pipeline --dry-run` 可预览。 |
| 大表排序、同名表和非有限值影响机器读取。 | `query` 通过 Arrow 过滤和投影，默认 50 行、最多 200 行，支持 `--omics` 或 `--mode` 消歧。排序窗口最多 20,000 行，响应预算约 128 KiB，非有限值回传 `null`。 |
| MCP 传输方式与研究数据外传风险需要各自核验。 | Phase 2 先验证受信任客户端的本机 `stdio`，再审计配置间接路径、脚本、外联和写入授权。 |

**本地验收命令**

```powershell
omicsweave --project-dir examples/demo_project pipeline --results-dir results/a1 --background --json
omicsweave --project-dir examples/demo_project status --results-dir results/a1 --job <job-id> --json
omicsweave --project-dir examples/demo_project logs --results-dir results/a1 <job-id> --since-event 0 --json
omicsweave --project-dir examples/demo_project inspect --results-dir results/a1 --experiment demo_glp1 --json
omicsweave --project-dir examples/demo_project query --results-dir results/a1 --experiment demo_glp1 --table differential --omics transcriptomics --limit 10 --json
```

**Phase 1 可自动化验收标准**：
1. 在 `examples/demo_project` 上提交后台流水线。目标环境空闲时的冷启动时延目标为 1 秒以内，stdout 提前返回含 `job_id` 的单一合法 JSON。记录目标机器、负载和实际时延。
2. 轮询 `omicsweave status --job <job_id> --json` 可观测到阶段推进（`run → integrate → summary → completed`）及心跳更新；强杀 Worker 进程后，`status` 能在下一次查询立即报告 `stalled` 而非停留在 `running`。
3. `inspect` 汇报有来源标识的摘要，`query` 返回真实结果表的分页切片，`trace` 定位单实体。Windows 与 Linux CI 均运行 A1 测试，并与既有基线对照。

**2026-11-25 本地验收记录**

- Windows Python 3.12.9 合成项目完成 `validate → run → integrate → summary → completed`。交付 ID 为 `run_20261125_002804`，运行清单中 37 个输出文件的 SHA-256 与磁盘一致。摘要发现 20 张 Parquet 表，转录组差异表共 60 行，首个 10 行分页返回 `next_offset=10`，实体 `ENSMUSG00000000001` 返回 2 条带表路径的证据。
- 去掉派发路径中统计库的预加载后，本机一次提交约 0.925 秒。验收的延迟口径包含启动 Python 与机器响应返回，部署环境需复测。
- 5 项 A1 测试通过，覆盖真实 Worker 协作取消与外部终止、活跃作业保留保护、原生子进程输出隔离、旧 `latest --json` 兼容、结果追溯和分页字节上限。独立 `run`、`integrate`、`map` 与 `abstract` 的 JSON 调用实测成功，图件命令均实际产出 PNG。两个进程并发更新共享指针 40 次后，两个实验的最终条目均保留。
- 实施前 Windows 全量回归为 687 通过、3 跳过、1 失败，实施后为 693 通过、3 跳过、1 失败，另有 473 个子测试通过。两轮的唯一失败均为 RefMet 果糖标识符冲突 `HMDB0000660 / HMDB0062538`。最后一处空结果目录清理兼容性修正后，重新通过 16 项相关测试及 7 个子测试，同时单独验证了空目录情形。
- CI 已包含 Windows/Linux 矩阵。本轮本机验证覆盖 Windows，Linux CI 结果待验证。MCP Server 属于 Phase 2。

---

### Phase 2（A2）— 强参数 Schema、安全写接口与 MCP Server v1

**目标**：实现完整的单机 AI Native 体验，让用户可直接通过任意 MCP 客户端以自然语言驱动 OmicsWeave。

- [x] **2.1 强类型参数注册表**。四组学、四整合模式及通路总图提供独立 Schema，校验报告记录未知参数的字段路径。保留已迁移配置中的已知历史字段。
- [x] **2.2 共享应用服务**。执行、作业状态、结果读取、配置补丁、样本 QC 及渲染由 CLI 与 MCP 复用。补丁与 QC 先预览，应用时核对原文件 SHA-256、变更理由及活跃作业。
- [x] **2.3 本机 MCP v1**。可选 SDK 固定为 `mcp==1.26.0`，支持 `stdio`、路径作用域、受控写入、混合同步/异步返回及可选 PNG 预览。MCP 清理提供预览，实际删除和旧版迁移在受信任的本机 CLI 执行。
- [x] **2.4 分析员 Skill**。在 `skills/omicsweave-analyst/SKILL.md` 记录发现、执行、结果查证及科学解释步骤。

**进入 Phase 2 的依据与实施边界**

- Phase 1 后台、JSON 响应和结果查询的 Windows 定向用例复跑通过。既有 RefMet 冲突和 Linux CI 结果作为验收记录跟踪。
- MCP 使用固定授权项目。显式添加的外部数据根用于读取，输出和配置写入位于项目内。应用层检查覆盖已知间接路径及结果包文件树。部署时仍需信任本机宿主、项目文件与明确放行的脚本。
- 计算、配置补丁和 QC 默认预览。MCP 提交写入要求服务器启动授权。实际写入按原文件 SHA-256 对照，合并补丁中 `null` 删除字段，预览显示 YAML 序列化后的差异。人工 QC 修改可能影响统计，建议从 `run` 重算。
- 多客户端远程化、覆盖参数的试验谱系和诊断工具属于 Phase 3。

```powershell
python -m pip install -e ".[mcp]"
omicsweave schema omics/mfa --json
omicsweave --project-dir examples/demo_project inspect-input --experiment demo_glp1 --json
omicsweave --project-dir examples/demo_project mcp serve
```

**Phase 2 可自动化验收标准**：
1. `test_param_schemas.py` 验证注册表与所有 Runner/Mode `DEFAULTS` 100% 同步；在 `experiment.yaml` 写入拼错的 `analysis` 参数时，`omicsweave validate` 在运行前即报错（退出码 `10` 或 `11`）。
2. `tests/mcp/test_mcp_protocol.py` 通过 In-Memory MCP Client 跑通完整端到端流程（发现 → 校验 → 提交 → 探活 → 读 Digest → 分页查询 → 重画图 → 对比），`test_mcp_security.py` 全部安全越界用例被拦截。
3. 仅安装基础依赖（不装 `mcp` extras）时，现有 CI 的基础导入测试仍然通过（`mcp` 延迟导入，零侵入基础包）。

**2026-11-25 本地验收记录**

- Windows Python 3.12.9 固定使用 `mcp==1.26.0`，4 项协议与权限测试通过，覆盖内存客户端完整流水线、真实 stdio 握手、PNG 预览、直接和间接路径约束，以及配置提交授权。Windows asyncio 验收使用允许本机回环的环境，MCP 消息通过 stdio 传输。
- 真实 stdio 合成交付为 `run_20261125_125552`。协商协议版本 `2025-11-25`，客户端发现 18 个工具，收到 4 条递增进度通知。转录组差异首批 10 行返回 `next_offset=10`，同结果包对照一致。证据保存在本地 `results/phase2-validation-20261125/stdio-8c6355e5/`。
- 全量回归的 10 项 Phase 2 兼容问题来自历史参数和既有报错文本。修复后 10 项全部复测通过，最终定向组共 31 项通过。全量运行同时复现了基线 RefMet 果糖标识符冲突。完整复测与 Linux CI 结果仍待记录。
- 基础环境验证了 CLI 的可选依赖隔离、新 Schema、输入探查、配置预览与应用及 QC 预览。DIO 和 DIO-GLP1RA 共三份现有实验配置通过静态参数校验。实际计算验证使用合成项目。

---

### Phase 3（A3）— 科研护栏迭代闭环与远程/多系统协同（预计 1–2 周，按需推进）

**目标**：支撑高质量、可审计的方法迭代优化，并为总纲 O5（远程集群）与 TheraWeave 双 Agent 协同铺路。

- [x] **3.1 零污染试验与运行谱系**。CLI 与 MCP 支持带理由的单次参数覆盖。试验结果写入 `<正式结果根>/_trials/<trial_id>/`，可用显式试验 ID 阅读。谱系保存为 `metadata/lineage.json`，作为既有 v0.1 清单中的输出文件登记。对照报告包含生效配置差异、谱系和证据等级转移。
- [x] **3.2 结构化诊断与防过拟合护栏**。CLI `diagnose` 和 MCP `diagnose_run` 共享只读诊断。输出包含样本 QC、比较可靠性、最小可达 p 值、知识表审计、校准状态和 `optimization_blocked`。正式或试验结果存在后，面板改动需要登记 `post_hoc` 及同实验的运行 ID。
- [ ] **3.3 远程 MCP 与生态衔接**：按需启用 `streamable-http` 传输 + Bearer Token 鉴权；评估以只读 MCP Resource 暴露 `AnalysisResultBundle` 契约，供未来 **TheraWeave Agent** 通过 MCP 或文件契约联合消费。

**Phase 3 可自动化验收标准**：
1. 使用 `--override` 跑敏感性试验后，核验工作区 `experiment.yaml`、正式 `latest.json`、源包及清单的 SHA-256。核验正式查询的扫描回退、保留策略、`metadata/experiment.json`，以及由清单登记的 `metadata/lineage.json`。
2. `diff_runs` 同时准确列出两包的配置改动与表/列/等级转移变化；`diagnose_run` 在诱饵期望接近的合成数据上准确产出 `statistical_guardrail` 且标记 `optimization_blocked: true`。

**3.1 与 3.2 的实施决定**

- 试验需要确切的正式基线和理由。试验可从 `run` 开始，或在原始单组学配置和数据有据可复用时从 `integrate` 开始。输入路径、面板文件和外部来源通过现行受控配置入口管理。
- 试验 Worker 使用本机缓存及离线基因集，跳过课题脚本和自动保留清理。试验保留引用的正式基线及复用包。研究者核定方案后可采用配置补丁和一次正式运行生成正式交付。
- v0.1 清单顶层字段保持既有契约，实际应用的配置、原始覆盖项、基线 ID 及指纹可从试验包查证。CLI 使用 `--trial` 读取试验，`diff` 使用 `--trial-a` 或 `--trial-b` 对照。后台作业继续用正式结果根下的 `status` 查询。
- 诱饵诊断报告观察 L1 数、诱饵均值、标准差和重复数。观察值位于诱饵均值加一倍标准差以内时标记描述性风险。样本或重复不足时保留未评估状态。已经登记的统计限制触发 `optimization_blocked`，对应操作为报告限制、复核设计。

```powershell
omicsweave --project-dir examples/demo_project pipeline --experiment demo_glp1 --source "<正式运行ID>" --override integration.parameters.reversal.equivalence_margin=0.72 --rationale "检查预先登记的等效性边界" --dry-run --json
omicsweave --project-dir examples/demo_project diagnose --experiment demo_glp1 --trial "<试验ID>" --json
omicsweave --project-dir examples/demo_project diff "<正式运行ID>" "<试验运行ID>" --experiment demo_glp1 --trial-b "<试验ID>" --json
```

**2026-11-25 本地验收**

- 合成项目的正式运行 `run_20261125_133209` 派生出后台试验 `trial_57ce75c7aa054ff8b363fb315c0c1938`。验收复核正式配置、指针、清单及所用转录组差异表的 SHA-256，七个试验输出均与清单哈希匹配。试验清单可通过现有 TheraWeave v0.1 manifest Schema。证据位于本地 `results/phase3-validation-20261125/check-98504518/final_acceptance.json`。
- Windows 定向组 38 项测试及 11 个子测试通过。面板来源校验加固后，诊断与配置定向组另有 5 项通过。MCP 组 4 项通过，覆盖试验预检和后台运行、结果对照与只读诊断。独立 `run` 和 `integrate` 的试验入口完成合成验收。
- 本轮验证覆盖 Windows 合成和定向回归。Linux CI 与整仓回归结果仍待记录。

---

## 3. 测试与质量门禁矩阵

| 测试类别 | 覆盖目标 | 关键断言 |
|---|---|---|
| **跨平台进程与并发测试** | `core/jobs.py`、`_worker.py` | Windows 与 Linux 下的子进程启动、UTF-8 日志无乱码、心跳刷新、PID 存活判断、强杀后识别 `stalled`、优雅取消 `cancelled`、并行双任务目录无冲突。 |
| **CLI 机器契约快照测试** | `cli.py --json` | 每个命令的 `stdout` 必为单个合法 JSON 且符合 `ServiceResponse` Schema；错误退出码严格匹配错误码表；非 JSON 默认模式保持人类输出不变。 |
| **Schema 一致性守卫测试** | `core/param_schemas.py` | 对比各 Runner / Integration Mode 的 `DEFAULTS` 键集与默认值，杜绝“代码改了参数但 Schema 忘了更新”。 |
| **查询安全与防爆炸测试** | `service/query.py` | 非法列名/运算符被拒；`NaN/Inf` 正确转 `null`；超大表分页与 128 KB 响应体截断有效；无 `eval` 注入面。 |
| **科研护栏回归测试** | `core/diagnostics.py`、`service/project.py` | 事后改 `panels.yaml` 强标 `post_hoc`；`--override` 不污染主配置和 `latest.json`；样本量不足/诱饵无区分时产出 `statistical_guardrail`。 |
| **MCP 协议与安全测试** | `omicsweave.mcp` | `stdio` 零杂音污染；路径穿越/符号链接越界被拒；默认 `dry_run=true`；未开启 `--allow-custom-scripts` 时跳过或拦截课题脚本；`wait_seconds` 超时正确降级为异步返回。 |

---

## 4. 风险矩阵与缓解策略

| 风险 | 影响 | 缓解策略 |
|---|---|---|
| **Agent 引发 p-hacking / 过度解读** | 科研结论失真、违背 OmicsWeave 科学严谨定位 | 诊断层强制区分 `engineering_fixable` / `sensitivity_explorable` / `statistical_guardrail`；强制保留诱饵期望与样本量口径；`panels.yaml` 事后修改自动记为 `post_hoc`；探索试验默认不覆盖 `latest.json` 并要求配对 `diff`。 |
| **未发表数据通过 MCP 或网络外泄** | 实验室数据泄露风险 | MCP 默认仅 `stdio` 本机通信；远程 HTTP 默认仅绑 `127.0.0.1` 且强鉴权；不内置任何外发遥测或云端上传。 |
| **Windows 进程管理与编码差异** | PID 探活失败、GBK 乱码污染 JSON-RPC | 子进程强制 `PYTHONIOENCODING=utf-8`；Windows 使用标准库 `ctypes` 调 Win32 API 探活而非 POSIX `os.kill`；CI 在 Windows 与 Linux 双平台跑全量 `test_jobs` 与 `test_mcp`。 |
| **长任务僵尸进程或孤立半成品目录** | 磁盘膨胀、状态误报 | 独立心跳线程 + 读取时 PID 探活自动修正 `stalled`；`retention.py` 识别并按策略清理僵尸作业与未完成目录。 |
| **大结果表撑爆 LLM 上下文** | Agent 调用失败、成本飙升、幻觉增多 | 强制金字塔式阅读（L0 Digest → L1 分页裁剪 Query → L2 Trace）；`query_table` 施加行数上限（200）和字节上限（128 KB）。 |
| **CLI 与 MCP 行为漂移** | 同操作结果不一致、维护成本翻倍 | 所有逻辑统一归口 `omicsweave.service`，CLI `--json` 与 MCP 仅为薄外壳，由共享 DTO 测试保障一致。 |

---

## 5. 待项目负责人确认事项

Phase 1 已采用第 1 项的兼容双入口，并以合成项目进行首轮实测。以下保留为 Phase 2 的集成议题，涉及真实课题时需确认数据和客户端使用范围。

1. **CLI 全局 JSON 参数命名**：采用 `--format {text,json}`（同时保留 `--json` 为简写别名，向后兼容 `latest --json`），还是所有子命令统一只提供 `--json`？
2. **异步作业默认同步等待窗口**：MCP `submit_*` 工具默认 `wait_seconds` 设为 `15` 秒（短任务直返、长任务转异步），是否符合日常实验节奏？
3. **课题脚本（`pipeline.scripts`）的 MCP 默认策略**：默认 `--allow-custom-scripts=false`（安全优先，仅在用户启动参数显式放行时允许执行 `analysis/*.py`），还是允许只读白名单脚本？
4. **首期实测靶标项目**：Phase 1/2 联调除合成项目 `examples/demo_project` 外，是否以本地真实课题（如 `local/projects/dio_project` 只读副本）作为 L0 Digest、诊断引擎与性能压测基准？
