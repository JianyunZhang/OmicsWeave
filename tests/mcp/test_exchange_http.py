"""Exercise a real loopback Streamable HTTP exchange server."""
from __future__ import annotations

import json
import socket
from datetime import timedelta

import pytest

pytest.importorskip("mcp")

import anyio
import httpx
import uvicorn
from mcp import ClientSession
from mcp.client.streamable_http import streamable_http_client
from pydantic import AnyUrl

from omicsweave.contract.exchange import export_exchange
from omicsweave.mcp.exchange_http import HOST, TOKEN_ENV, build_http_app
from omicsweave.service.schemas import ErrorCode, ServiceError
from tests.contract.test_exchange import derived_run


def test_http_exchange_requires_token_and_exposes_only_paged_tables(derived_run, tmp_path, monkeypatch):
    _, _, run = derived_run
    output = tmp_path / "bundle"
    snapshot = export_exchange(run, output)
    table = next(item.path for item in snapshot.artifacts if item.path.endswith("/derived_differential.parquet"))
    matrix = next(item.path for item in snapshot.artifacts if item.kind == "matrix")
    token = "local-synthetic-test-token-1234567890"
    monkeypatch.setenv(TOKEN_ENV, token)

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        listener.bind((HOST, 0))
        port = listener.getsockname()[1]
        listener.listen()
        server = uvicorn.Server(uvicorn.Config(build_http_app(output, port), log_level="warning", proxy_headers=False))
        url = f"http://{HOST}:{port}/mcp"

        async def check():
            async with anyio.create_task_group() as tasks:
                tasks.start_soon(server.serve, [listener])
                try:
                    with anyio.fail_after(40):
                        while not server.started:
                            await anyio.sleep(0.025)
                        async with httpx.AsyncClient(trust_env=False) as denied:
                            unauthenticated = await denied.get(url)
                            assert unauthenticated.status_code == 401
                            assert unauthenticated.headers.get("WWW-Authenticate") == "Bearer"
                            wrong_origin = await denied.get(url, headers={"Origin": "https://example.invalid",
                                                                          "Authorization": f"Bearer {token}"})
                            assert wrong_origin.status_code == 403
                            wrong_host = await denied.get(url, headers={"Host": "example.invalid",
                                                                         "Authorization": f"Bearer {token}"})
                            assert wrong_host.status_code == 421
                        async with httpx.AsyncClient(trust_env=False, headers={"Authorization": f"Bearer {token}"}) as http:
                            async with streamable_http_client(url, http_client=http) as (read, write, _):
                                async with ClientSession(read, write, read_timeout_seconds=timedelta(seconds=20)) as client:
                                    await client.initialize()
                                    tools = await client.list_tools()
                                    assert {tool.name for tool in tools.tools} == {"query_exchange_table"}
                                    resource = await client.read_resource(AnyUrl("omicsweave://exchange/bundle"))
                                    descriptor = json.loads(resource.contents[0].text)
                                    assert descriptor["run"]["run_id"] == run.name
                                    assert str(run) not in resource.contents[0].text
                                    first = await client.call_tool("query_exchange_table", {
                                        "path": table, "columns": ["feature_id", "log2fc"], "limit": 1,
                                    })
                                    assert not first.isError, first.content
                                    page = first.structuredContent["data"]
                                    assert len(page["rows"]) == 1 and page["next_offset"] == 1
                                    second = await client.call_tool("query_exchange_table", {
                                        "path": table, "columns": ["feature_id"], "limit": 1, "offset": page["next_offset"],
                                    })
                                    assert second.structuredContent["data"]["rows"][0]["feature_id"] != page["rows"][0]["feature_id"]
                                    denied_matrix = await client.call_tool("query_exchange_table", {"path": matrix})
                                    assert denied_matrix.isError
                                    assert denied_matrix.structuredContent["error"]["code"] == "INPUT_NOT_FOUND"
                finally:
                    server.should_exit = True

        anyio.run(check)


def test_http_exchange_rejects_missing_token(tmp_path, monkeypatch):
    monkeypatch.delenv(TOKEN_ENV, raising=False)
    with pytest.raises(ServiceError) as error:
        build_http_app(tmp_path, 8765)
    assert error.value.detail.code == ErrorCode.SECURITY_POLICY_VIOLATION
