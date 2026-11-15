#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
由 RefMet（Fahy 与 Subramaniam 2020，Metabolomics Workbench REST）生成代谢物注册表 ``data/metabolite_registry.tsv.gz``。

每个查询名称：``refmet/match/<名称>`` 得到 RefMet 标准名称与类别，``refmet/name/<标准名称>/all`` 得到 RefMet 编号与
regno，``compound/regno/<regno>/all`` 得到 HMDB / KEGG / ChEBI / PubChem。状态：``matched``（有数据库 ID）、
``no_ids``（RefMet 认得但没有 ID）、``unmatched``（RefMet 不认得）。

默认查询整合实际用到的名称（反应 → 代谢物映射中的全部别名），``--all`` 时另加手写 ID 表中的全部名称。
网络响应缓存在 ``--cache``（JSON），中断后重跑会跳过已查询的名称。来源与生成日期写入 ``data/SOURCES.json``
的 ``metabolites`` 一节。

用法::

    python -m omicsweave.core.idmap.build_metabolite_registry [--all] [--cache refmet_cache.json]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path
from typing import Callable

DATA_DIR = Path(__file__).resolve().parent / "data"
TABLE_FILE = "metabolite_registry.tsv.gz"
BASE_URL = "https://www.metabolomicsworkbench.org/rest"
COLUMNS = ["query_name", "refmet_name", "refmet_id", "regno", "formula", "exact_mass", "super_class", "main_class",
           "sub_class", "hmdb_id", "kegg_compound", "chebi_id", "pubchem_cid", "inchi_key", "status"]
Fetch = Callable[[str], object]


def http_fetch(path: str, timeout: float = 30.0) -> object:
    url = f"{BASE_URL}/{urllib.parse.quote(path, safe='/')}"
    request = urllib.request.Request(url, headers={"User-Agent": "omicsweave-registry-builder"})
    with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 — 固定的公开 REST 地址
        text = response.read().decode("utf-8").strip()
    return json.loads(text) if text and text not in ("[]", "{}") else {}


def _first(value: object) -> dict:
    if isinstance(value, list):
        return value[0] if value else {}
    if isinstance(value, dict) and value and all(isinstance(v, dict) for v in value.values()):
        return next(iter(value.values()))
    return value if isinstance(value, dict) else {}


def resolve(name: str, fetch: Fetch) -> dict[str, str]:
    row = dict.fromkeys(COLUMNS, "")
    row["query_name"] = name
    match = _first(fetch(f"refmet/match/{name}"))
    standard = str(match.get("refmet_name", "") or "")
    if not standard or standard == "-":
        row["status"] = "unmatched"
        return row
    row.update(refmet_name=standard, formula=str(match.get("formula", "")), exact_mass=str(match.get("exactmass", "")),
               super_class=str(match.get("super_class", "")), main_class=str(match.get("main_class", "")),
               sub_class=str(match.get("sub_class", "")))
    record = _first(fetch(f"refmet/name/{standard}/all"))
    row.update(refmet_id=str(record.get("refmet_id", "")), regno=str(record.get("regno", "")),
               pubchem_cid=str(record.get("pubchem_cid", "")), inchi_key=str(record.get("inchi_key", "")))
    compound = _first(fetch(f"compound/regno/{row['regno']}/all")) if row["regno"] else {}
    hmdb = str(compound.get("hmdb_id", "") or "")
    kegg = str(compound.get("kegg_id", "") or "")
    chebi = str(compound.get("chebi_id", "") or "")
    row.update(hmdb_id=hmdb, kegg_compound=kegg, chebi_id=f"CHEBI:{chebi}" if chebi and not chebi.startswith("CHEBI")
               else chebi)
    row["status"] = "matched" if hmdb or kegg or chebi else "no_ids"
    return row


def query_names(include_all: bool = False) -> list[str]:
    """整合用到的名称：反应 → 代谢物映射中各角色的全部别名；``include_all`` 时加上手写 ID 表的全部名称。"""
    from omicsweave.integration.mappings.flux_metabolite_mapping import FLUX_METABOLITE_MAPPING

    names: list[str] = []
    for roles in FLUX_METABOLITE_MAPPING.values():
        if isinstance(roles, dict):
            for role in ("substrates", "products", "related"):
                names += [str(n) for n in roles.get(role, [])]
    if include_all:
        from omicsweave.core.idmap.metabolite_id_mapping import METABOLITE_DATABASE_IDS
        names += list(METABOLITE_DATABASE_IDS)
    return sorted(dict.fromkeys(n.strip() for n in names if n and n.strip()))


def build(names: list[str], fetch: Fetch = http_fetch, *, cache: Path | None = None, workers: int = 4
          ) -> list[dict[str, str]]:
    cached: dict[str, object] = json.loads(cache.read_text(encoding="utf-8")) if cache and cache.is_file() else {}

    def remembered(path: str) -> object:
        if path not in cached:
            try:
                cached[path] = fetch(path)
            except Exception as exc:  # noqa: BLE001 — 单个名称失败不中断，状态记为 error
                return {"_error": str(exc)}
        return cached[path]

    def one(name: str) -> dict[str, str]:
        row = resolve(name, remembered)
        return row

    with ThreadPoolExecutor(max_workers=max(workers, 1)) as pool:
        rows = list(pool.map(one, names))
    if cache:
        cache.write_text(json.dumps(cached, ensure_ascii=False), encoding="utf-8")
    return rows


def write(rows: list[dict[str, str]], out_dir: Path = DATA_DIR) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / TABLE_FILE
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: r["query_name"].lower()))
    sources_path = out_dir / "SOURCES.json"
    sources = json.loads(sources_path.read_text(encoding="utf-8")) if sources_path.is_file() else {}
    status = {s: sum(r["status"] == s for r in rows) for s in ("matched", "no_ids", "unmatched")}
    sources["metabolites"] = {
        "generated_by": "python -m omicsweave.core.idmap.build_metabolite_registry",
        "generated_on": date.today().isoformat(), "source": "RefMet via Metabolomics Workbench REST",
        "endpoints": [f"{BASE_URL}/refmet/match/<name>", f"{BASE_URL}/refmet/name/<refmet_name>/all",
                      f"{BASE_URL}/compound/regno/<regno>/all"],
        "n_names": len(rows), "status": status,
    }
    sources_path.write_text(json.dumps(sources, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="由 RefMet 生成代谢物注册表")
    parser.add_argument("--all", action="store_true", help="另查询手写 ID 表中的全部名称")
    parser.add_argument("--cache", type=Path, default=None, help="RefMet 响应缓存（JSON），中断后可续跑")
    parser.add_argument("--out", type=Path, default=DATA_DIR, help="输出目录（默认为本包的 data/）")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    names = query_names(args.all)
    rows = build(names, cache=args.cache, workers=args.workers)
    path = write(rows, args.out)
    counts = {s: sum(r["status"] == s for r in rows) for s in ("matched", "no_ids", "unmatched")}
    print(f"{len(rows)} 个名称：{counts}；已写出 {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
