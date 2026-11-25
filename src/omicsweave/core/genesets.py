#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因集（通路）库：读取本地 GMT 文件，或按 Enrichr 库名下载后缓存为 GMT。

按库名使用时只下载库本身（``gseapy.get_library``），不上传任何基因列表；首次下载后缓存在
``<paths.references>/gene_sets/<库名>.gmt``，之后的运行读取缓存文件，运行清单会登记它的校验值，
保证同一结果包可以用同一版本的库复现。
"""
from __future__ import annotations

import re
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

from .config import Project
from .model import Species

GENE_SET_DIR = "gene_sets"
# Enrichr 的物种名
ENRICHR_ORGANISMS = {Species.MOUSE: "Mouse", Species.HUMAN: "Human"}
_OFFLINE = ContextVar("omicsweave_offline_genesets", default=False)


@contextmanager
def offline_libraries():
    token = _OFFLINE.set(True)
    try:
        yield
    finally:
        _OFFLINE.reset(token)


class GeneSetError(Exception):
    """基因集库无法读取或下载。"""


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """GMT：每行 ``条目<TAB>描述<TAB>成员1<TAB>成员2…``；保持文件中的条目顺序，去掉重复与空成员。"""
    gene_sets: dict[str, list[str]] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 3 or not fields[0].strip():
            continue
        gene_sets[fields[0].strip()] = list(dict.fromkeys(m.strip() for m in fields[2:] if m.strip()))
    if not gene_sets:
        raise GeneSetError(f"{path} 中没有读到任何基因集（应为 GMT 格式）")
    return gene_sets


def write_gmt(gene_sets: dict[str, list[str]], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{term}\t\t" + "\t".join(members) + "\n" for term, members in gene_sets.items()),
                    encoding="utf-8")
    return path


def library_name(source: str) -> str:
    """库名：GMT 文件取文件名（不含扩展名），否则就是 Enrichr 库名。"""
    return Path(source).stem if source.lower().endswith(".gmt") else source


def cache_path(project: Project, name: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-") or "gene_sets"
    return project.path("references") / GENE_SET_DIR / f"{safe}.gmt"


def load_libraries(project: Project, sources: list[str]) -> tuple[dict[str, dict[str, list[str]]], list[Path]]:
    """读取基因集库，返回 ({库名: {条目: 成员}}, 用到的 GMT 文件)。

    ``sources`` 中以 ``.gmt`` 结尾的按项目目录解析为本地文件，其余视为 Enrichr 库名。
    """
    libraries: dict[str, dict[str, list[str]]] = {}
    files: list[Path] = []
    for source in sources:
        source = str(source)
        if source.lower().endswith(".gmt"):
            path = project.resolve(source)
            if not path.is_file():
                raise GeneSetError(f"找不到基因集文件 {path}")
        else:
            path = cache_path(project, source)
            if not path.is_file():
                if _OFFLINE.get():
                    raise GeneSetError(f"离线执行需要已缓存的基因集 {source}")
                _download(source, ENRICHR_ORGANISMS[project.config.project.species], path)
        libraries[library_name(source)] = read_gmt(path)
        files.append(path)
    return libraries, files


def _download(name: str, organism: str, path: Path) -> None:
    try:
        import gseapy
    except ImportError as exc:
        raise GeneSetError(f"按库名使用基因集 {name} 需要 gseapy；也可以改用本地 GMT 文件") from exc
    try:
        gene_sets = gseapy.get_library(name=name, organism=organism)
    except Exception as exc:
        raise GeneSetError(f"下载基因集库 {name}（{organism}）失败：{exc}。可以改用本地 GMT 文件，"
                           f"或把下载好的库放到 {path}") from exc
    write_gmt(gene_sets, path)
