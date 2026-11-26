"""Single-package Streamable HTTP pilot bound to loopback."""
from __future__ import annotations

import hmac
import os
import re
from pathlib import Path

from mcp.types import CallToolResult, TextContent, ToolAnnotations

from omicsweave.contract.exchange import AnalysisResultExchange, ExchangeSnapshot, open_exchange
from omicsweave.service.schemas import ErrorCode, ErrorDetail, ServiceError, ServiceResponse

TOKEN_ENV = "OMICSWEAVE_MCP_TOKEN"
HOST = "127.0.0.1"
DEFAULT_PORT = 8765


def _token() -> bytes:
    value = os.environ.get(TOKEN_ENV, "")
    if not re.fullmatch(r"[A-Za-z0-9._~+/=-]{32,}", value):
        raise ServiceError(ErrorCode.SECURITY_POLICY_VIOLATION,
                           f"请通过环境变量 {TOKEN_ENV} 提供至少 32 字符的随机 Bearer Token")
    return value.encode("ascii")


def _public_columns(path: str) -> set[str]:
    from pathlib import PurePosixPath

    from omicsweave.contract.tables import DIFFERENTIAL_COLUMNS, ENRICHMENT_COLUMNS

    relative = PurePosixPath(path)
    if relative.parts[0] == "omics":
        columns = ENRICHMENT_COLUMNS if relative.stem == "enrichment" else DIFFERENTIAL_COLUMNS
        return {item.name for item in columns} | {"feature_type", "related_reaction", "tracer"}
    return {
        "omics", "tissue", "contrast_id", "comparison", "pair_id", "feature_id", "feature_type",
        "term_id", "term_name", "pathway", "pathway_term", "pvalue", "padj", "direction",
        "direction_agreement", "flux_id", "flux_claim", "claim_direction", "mechanism_class",
        "axis", "state", "effect", "ci_lower", "ci_upper", "n_sources", "n_concordant_axes",
    }


def read_exchange_page(snapshot: ExchangeSnapshot, path: str, columns: list[str] | None, limit: int, offset: int) -> dict:
    import pyarrow.parquet as pq
    from pydantic import ValidationError

    from omicsweave.service.query import TableQuery, read_table_page

    try:
        request = TableQuery(table=path, columns=columns, limit=limit, offset=offset)
    except ValidationError as exc:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "分页参数无效") from exc
    artifact = next((item for item in snapshot.manifest.artifacts if item.kind == "table" and item.path == path), None)
    if artifact is None:
        raise ServiceError(ErrorCode.INPUT_NOT_FOUND, "该路径没有可通过 HTTP 查询的已登记结果表")
    table_path = snapshot.artifact_path(artifact)
    columns_allowed = _public_columns(artifact.path)
    available = [name for name in pq.read_schema(table_path).names if name in columns_allowed]
    selected = columns if columns is not None else available
    if not selected or any(name not in available for name in selected):
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "请选择当前结果表中可交换的列")
    return read_table_page(snapshot.manifest.run.run_id, artifact.path, table_path,
                           request.model_copy(update={"columns": selected}))


def build_http_app(directory: Path, port: int):
    """Only registers published metadata and paged aggregate results."""
    from mcp.server.fastmcp import FastMCP
    from mcp.server.transport_security import TransportSecuritySettings
    from starlette.responses import JSONResponse

    if not 1 <= port <= 65535:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "--port 须在 1–65535")
    token = _token()
    snapshot = open_exchange(directory)
    hosts_allowed = [f"{HOST}:{port}", HOST] if port == 80 else [f"{HOST}:{port}"]
    origins_allowed = [f"http://{host}" for host in hosts_allowed]
    server = FastMCP(
        "OmicsWeave Exchange", instructions="只读已发布的实验结果快照。",
        host=HOST, port=port, stateless_http=True, json_response=True,
        transport_security=TransportSecuritySettings(allowed_hosts=hosts_allowed, allowed_origins=origins_allowed),
    )

    @server.resource("omicsweave://exchange/bundle", mime_type="application/json")
    def bundle() -> dict:
        return snapshot.manifest.model_dump(mode="json")

    @server.resource("omicsweave://exchange/schema", mime_type="application/json")
    def schema() -> dict:
        return AnalysisResultExchange.model_json_schema()

    @server.tool(annotations=ToolAnnotations(readOnlyHint=True, destructiveHint=False,
                                            idempotentHint=True, openWorldHint=False))
    async def query_exchange_table(path: str, columns: list[str] | None = None, limit: int = 50,
                                   offset: int = 0) -> CallToolResult:
        """读取一个已发布的汇总表，按结果行分页。"""
        from anyio import to_thread

        error = None
        data = None
        try:
            data = await to_thread.run_sync(read_exchange_page, snapshot, path, columns, limit, offset)
        except ServiceError as exc:
            error = exc.detail
        except Exception:
            # Transport errors never disclose the host filesystem or Python exception text.
            error = ErrorDetail(code=ErrorCode.CONTRACT_ERROR, message="交换表读取或完整性校验失败")
        response = ServiceResponse(command="query_exchange_table", project_dir=None,
                                   ok=error is None, data=data, error=error)
        return CallToolResult(isError=error is not None, structuredContent=response.model_dump(mode="json"),
                              content=[TextContent(type="text", text=response.model_dump_json())])

    app = server.streamable_http_app()
    host_bytes = {host.encode("ascii") for host in hosts_allowed}
    origin_bytes = {origin.encode("ascii") for origin in origins_allowed}

    async def guarded(scope, receive, send):
        if scope["type"] == "http":
            headers = scope["headers"]
            hosts = [value for name, value in headers if name == b"host"]
            origins = [value for name, value in headers if name == b"origin"]
            auth = [value for name, value in headers if name == b"authorization"]
            scheme, separator, credential = auth[0].partition(b" ") if len(auth) == 1 else (b"", b"", b"")
            if (scope.get("client") or ("",))[0] != HOST:
                response = JSONResponse({"error": "loopback_required"}, status_code=403)
            elif len(hosts) != 1 or hosts[0] not in host_bytes:
                response = JSONResponse({"error": "invalid_host"}, status_code=421)
            elif origins and (len(origins) != 1 or origins[0] not in origin_bytes):
                response = JSONResponse({"error": "invalid_origin"}, status_code=403)
            elif scheme.lower() != b"bearer" or not separator or not hmac.compare_digest(credential, token):
                response = JSONResponse({"error": "unauthorized"}, status_code=401,
                                        headers={"WWW-Authenticate": "Bearer"})
            else:
                await app(scope, receive, send)
                return
            await response(scope, receive, send)
        else:
            await app(scope, receive, send)

    return guarded


def serve_http(directory: Path, port: int) -> None:
    import uvicorn

    uvicorn.run(build_http_app(directory, port), host=HOST, port=port, proxy_headers=False, access_log=False)
