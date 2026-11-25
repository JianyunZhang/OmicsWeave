"""本地 MCP v1 启动入口，仅在启动服务时加载可选 SDK。"""
from __future__ import annotations

import argparse
import os
import sys
from contextlib import redirect_stdout
from pathlib import Path

from omicsweave.core.path_scope import PathScope
from omicsweave.service.schemas import classify_error


def build_server(project: Path, data_roots=(), *, allow_writes=False, allow_custom_scripts=False):
    from mcp.server.fastmcp import FastMCP

    from .prompts import register_prompts
    from .resources import register_resources
    from .security import ServerPolicy
    from .tools import register_tools

    scope = PathScope.create(project, data_roots)
    policy = ServerPolicy(scope, allow_writes, allow_custom_scripts)
    server = FastMCP("OmicsWeave", instructions="在授权的本机项目中发现、预检、执行和解读多组学分析。")
    register_tools(server, policy)
    register_resources(server, policy)
    register_prompts(server)
    return server


def serve(project: Path, data_roots=(), *, allow_writes=False, allow_custom_scripts=False) -> int:
    # 保留独立的 stdout 句柄给 JSON-RPC，诊断进入 stderr。
    sys.stdout.flush()
    protocol_fd = os.dup(1)
    os.dup2(2, 1)
    with os.fdopen(protocol_fd, "w", encoding="utf-8", newline="\n") as protocol, redirect_stdout(sys.stderr):
        try:
            import anyio
            from mcp.server.stdio import stdio_server

            server = build_server(project, data_roots, allow_writes=allow_writes,
                                  allow_custom_scripts=allow_custom_scripts)

            async def run():
                async with stdio_server(stdout=anyio.wrap_file(protocol)) as (receive, send):
                    core = server._mcp_server  # SDK v1 的 FastMCP 未公开注入 stdout 的入口
                    await core.run(receive, send, core.create_initialization_options())

            anyio.run(run)
            return 0
        except ModuleNotFoundError as exc:
            if exc.name in {"mcp", "anyio"}:
                print('MCP 依赖缺失，请安装 "omicsweave[mcp]"', file=sys.stderr)
                return 14
            error = classify_error(exc)
            print(error, file=sys.stderr)
            return error.exit_code
        except Exception as exc:
            error = classify_error(exc)
            print(error, file=sys.stderr)
            return error.exit_code
        finally:
            os.dup2(protocol.fileno(), 1)


def main() -> None:
    parser = argparse.ArgumentParser(prog="omicsweave-mcp")
    parser.add_argument("--project-dir", type=Path, required=True)
    parser.add_argument("--allowed-data-root", type=Path, action="append", default=[])
    parser.add_argument("--allow-writes", action="store_true")
    parser.add_argument("--allow-custom-scripts", action="store_true")
    args = parser.parse_args()
    sys.exit(serve(args.project_dir, args.allowed_data_root, allow_writes=args.allow_writes,
                   allow_custom_scripts=args.allow_custom_scripts))


if __name__ == "__main__":
    main()
