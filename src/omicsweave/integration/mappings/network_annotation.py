"""MFA 网络的反应注释：反应类别、区室、底物 / 产物，以及各组学层能否给出方向证据。

注释表 ``data/network_reactions.tsv`` 由 ``python -m omicsweave.integration.mappings.network_annotation``
从 MFA 项目的网络源码（每个模块里形如 ``{'id': ..., 'sub': [...], 'pro': [...], 'reverse': ...}`` 的反应字典）
用 ``ast`` 静态解析生成，不执行源码；来源与 git 提交记在 ``data/network_SOURCES.json``。

反应类别（``reaction_class``）：

- ``enzymatic``：底物与产物不是同一组代谢物（含集总反应，如 LIPID_c、FA_esterification_c）；
- ``transport``：两侧是同一组代谢物、都在组织内（胞质 ⇌ 线粒体，含反向转运体）；
- ``exchange``：两侧是同一组代谢物，一侧在血浆（``_e``）；
- ``extracellular``：只涉及血浆代谢物（如 LDH_e），组织组学无法评估；
- ``pseudo``：输注、蛋白分解输入、原有储存输入、BIOMASS 等建模用的虚拟反应。

方向证据（``direction_evidence``）：``directional`` 时直接酶 / 载体与底物 / 产物可给出方向；``capacity`` 时它们的变化
只表示容量或关联，不判定支持或冲突。默认规则与逐反应的例外写在 ``flux_reaction_semantics.yaml``。

网络中每个代谢物 ID 形如 ``G3P_c``（基名 + 区室）；组织侧（``_c``、``_m``）的底物 / 产物换成标准名称写入注释表，
供代谢物映射表的底物 / 产物角色与网络化学计量保持一致（``sync_metabolite_roles``）。
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
import re
import subprocess
from dataclasses import dataclass
from datetime import date
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml

from omicsweave.core.idmap.metabolite_id_mapping import lookup_metabolite_ids

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
TABLE_FILE = "network_reactions.tsv"
SOURCES_FILE = "network_SOURCES.json"
SEMANTICS_FILE = HERE / "flux_reaction_semantics.yaml"

TISSUE_COMPARTMENTS = ("c", "m")
PSEUDO_ID = re.compile(r"(?:_infuse$|_protein_input$|_stock_input$|^BIOMASS|_intestinal_synthesis$|^NH4_input$|"
                       r"_excretion$)")
_COMPARTMENT = re.compile(r"^(?P<base>.+)_(?P<comp>[cme])$")
COLUMNS = ["flux_id", "network", "module", "substrates", "products", "reversible", "reaction_class", "compartments",
           "tissue_substrate_names", "tissue_product_names", "n_variants"]


# ── 解析与分类 ────────────────────────────────────────────────────────────────────────────────────────────────

def split_metabolite(metabolite: str) -> tuple[str, str]:
    """``'G3P_c'`` → ``('G3P', 'c')``；没有区室后缀的池（LIPID、TG_store）区室为空。"""
    match = _COMPARTMENT.match(metabolite)
    return (match["base"], match["comp"]) if match else (metabolite, "")


def _constant(node: ast.AST) -> Any:
    return node.value if isinstance(node, ast.Constant) else None


def _species(node: ast.AST) -> list[str] | None:
    if not isinstance(node, (ast.List, ast.Tuple)):
        return None
    out = []
    for element in node.elts:
        if not isinstance(element, ast.Tuple) or not element.elts:
            return None
        name = _constant(element.elts[0])
        if not isinstance(name, str):
            return None
        out.append(name)
    return out


def parse_module(source: str) -> list[dict[str, Any]]:
    """从一个网络模块的源码中取出字面量反应字典（``id`` 为字符串常量者）；f-string 生成的反应被跳过。"""
    reactions = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Dict):
            continue
        fields = {_constant(k): v for k, v in zip(node.keys, node.values) if k is not None}
        flux_id = _constant(fields["id"]) if "id" in fields else None
        if not isinstance(flux_id, str) or "sub" not in fields or "pro" not in fields:
            continue
        substrates, products = _species(fields["sub"]), _species(fields["pro"])
        if substrates is None or products is None:
            continue
        reverse = _constant(fields["reverse"]) if "reverse" in fields else False
        reactions.append({"flux_id": flux_id, "substrates": tuple(substrates), "products": tuple(products),
                          "reversible": bool(reverse)})
    return reactions


def classify(flux_id: str, substrates: Iterable[str], products: Iterable[str]) -> str:
    substrates, products = list(substrates), list(products)
    if PSEUDO_ID.search(flux_id):
        return "pseudo"
    if substrates and all(re.search(r"_(?:stock|inf)$", m) for m in substrates):
        return "pseudo"                             # 未标记的原有储存或示踪剂输入
    parts = [split_metabolite(m) for m in substrates + products]
    if not any(comp in TISSUE_COMPARTMENTS for _, comp in parts):
        return "extracellular" if any(comp == "e" for _, comp in parts) else "pseudo"
    sub_bases = sorted(split_metabolite(m)[0] for m in substrates)
    pro_bases = sorted(split_metabolite(m)[0] for m in products)
    if sub_bases and sub_bases == pro_bases:
        comps = {comp for _, comp in parts}
        return "exchange" if "e" in comps else "transport"
    return "enzymatic"


def tissue_side(substrates: Iterable[str], products: Iterable[str]) -> tuple[list[str], list[str]]:
    """组织侧的底物与产物基名；两侧都出现的基名（转运、反向转运）不计入任一侧。"""
    sub = [split_metabolite(m) for m in substrates]
    pro = [split_metabolite(m) for m in products]
    sub_bases = {b for b, c in sub if c in TISSUE_COMPARTMENTS}
    pro_bases = {b for b, c in pro if c in TISSUE_COMPARTMENTS}
    shared = sub_bases & pro_bases
    ordered = lambda items, keep: list(dict.fromkeys(b for b, c in items if c in TISSUE_COMPARTMENTS and b in keep))
    return ordered(sub, sub_bases - shared), ordered(pro, pro_bases - shared)


def read_standard_names(path: Path) -> dict[str, str]:
    """MFA 项目的模型代谢物基名 → 标准名称（``model_metabolite_to_standard_name_dict.py``）。"""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Dict):
            for key, value in zip(node.keys, node.values):
                if isinstance(_constant(key), str) and isinstance(_constant(value), str):
                    names.setdefault(_constant(key), _constant(value))
    return names


def resolvable_name(name: str) -> str:
    """标准名称的一个能被代谢物 ID 表解析的写法（``fructose 6-phosphate`` → ``fructose-6-phosphate``）；都不行时原样返回。"""
    if not name or lookup_metabolite_ids(name):
        return name
    candidates = [name.replace(" ", "-"), re.sub(r"^a-", "alpha-", name)]
    candidates.append(re.sub(r"^a-", "alpha-", name.replace(" ", "-")))
    for candidate in candidates:
        if lookup_metabolite_ids(candidate):
            return candidate
    return name


def build_table(model_dirs: Mapping[str, Path], standard_names: Mapping[str, str]) -> tuple[list[dict], dict]:
    rows: dict[tuple[str, str], dict[str, Any]] = {}
    variants: dict[tuple[str, str], set[tuple]] = {}
    for network, folder in model_dirs.items():
        for path in sorted(Path(folder).glob("*.py")):
            for reaction in parse_module(path.read_text(encoding="utf-8")):
                key = (network, reaction["flux_id"])
                variants.setdefault(key, set()).add((reaction["substrates"], reaction["products"]))
                rows.setdefault(key, {**reaction, "network": network, "module": path.stem})
    out, unnamed = [], set()
    for key, row in sorted(rows.items(), key=lambda kv: (kv[0][0] != "carbon", kv[0][1])):
        subs, pros = tissue_side(row["substrates"], row["products"])
        cls = classify(row["flux_id"], row["substrates"], row["products"])
        if cls in ("transport", "pseudo", "extracellular"):
            subs, pros = [], []                     # 转运两侧是同一个池；虚拟与血浆反应没有组织侧角色
        names_for = lambda bases: [resolvable_name(standard_names.get(b, "")) for b in bases]
        for base in subs + pros:
            if base not in standard_names:
                unnamed.add(base)
        out.append({
            "flux_id": row["flux_id"], "network": row["network"], "module": row["module"],
            "substrates": ";".join(row["substrates"]), "products": ";".join(row["products"]),
            "reversible": str(row["reversible"]).lower(), "reaction_class": cls,
            "compartments": ";".join(sorted({split_metabolite(m)[1] or "pool"
                                             for m in row["substrates"] + row["products"]})),
            "tissue_substrate_names": ";".join(n for n in names_for(subs) if n),
            "tissue_product_names": ";".join(n for n in names_for(pros) if n),
            "n_variants": len(variants[key]),
        })
    report = {
        "n_reactions": {net: sum(r["network"] == net for r in out) for net in model_dirs},
        "reaction_classes": {cls: sum(r["reaction_class"] == cls for r in out)
                             for cls in ("enzymatic", "transport", "exchange", "extracellular", "pseudo")},
        "reactions_with_several_definitions": sorted(f"{r['network']}:{r['flux_id']}" for r in out
                                                     if int(r["n_variants"]) > 1),
        "tissue_metabolites_without_standard_name": sorted(unnamed),
        "standard_names_without_database_id": sorted({n for r in out if r["network"] == "carbon"
                                                      for c in ("tissue_substrate_names", "tissue_product_names")
                                                      for n in _split(r[c]) if not lookup_metabolite_ids(n)}),
    }
    return out, report


# ── 读取 ────────────────────────────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class NetworkReaction:
    flux_id: str
    network: str
    module: str
    substrates: tuple[str, ...]
    products: tuple[str, ...]
    reversible: bool
    reaction_class: str
    compartments: tuple[str, ...]
    tissue_substrate_names: tuple[str, ...]
    tissue_product_names: tuple[str, ...]


@dataclass(frozen=True)
class ReactionSemantics:
    flux_id: str
    reaction_class: str          # enzymatic / transport / exchange / extracellular / pseudo / unknown
    direction_evidence: str      # directional / capacity / not_assessable
    reason: str
    in_network: bool

    @property
    def assessable(self) -> bool:
        return self.direction_evidence != "not_assessable"

    @property
    def directional(self) -> bool:
        return self.direction_evidence == "directional"


def _split(value: str) -> tuple[str, ...]:
    return tuple(piece for piece in (value or "").split(";") if piece)


@lru_cache(maxsize=1)
def network_reactions() -> dict[str, NetworkReaction]:
    """``flux_id`` → 网络反应；碳网络优先，只出现在氮网络中的反应也收录。"""
    path = DATA_DIR / TABLE_FILE
    if not path.exists():
        return {}
    out: dict[str, NetworkReaction] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["flux_id"] in out:
                continue
            out[row["flux_id"]] = NetworkReaction(
                row["flux_id"], row["network"], row["module"], _split(row["substrates"]), _split(row["products"]),
                row["reversible"] == "true", row["reaction_class"], _split(row["compartments"]),
                _split(row["tissue_substrate_names"]), _split(row["tissue_product_names"]))
    return out


@lru_cache(maxsize=1)
def network_sources() -> dict[str, Any]:
    path = DATA_DIR / SOURCES_FILE
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


@lru_cache(maxsize=1)
def _semantics_config() -> dict[str, Any]:
    return yaml.safe_load(SEMANTICS_FILE.read_text(encoding="utf-8")) or {}


@lru_cache(maxsize=None)
def reaction_semantics(flux_id: str) -> ReactionSemantics:
    """一个反应的类别与方向证据规则（逐反应例外 > 按类别的默认值）。不在网络注释中的反应按 ID 推断。"""
    config = _semantics_config()
    reaction = network_reactions().get(flux_id)
    if reaction is not None:
        cls = reaction.reaction_class
    elif PSEUDO_ID.search(flux_id):
        cls = "pseudo"
    else:
        cls = "unknown"
    override = (config.get("overrides") or {}).get(flux_id) or {}
    defaults = config.get("defaults") or {}
    evidence = override.get("direction_evidence") or defaults.get(cls, "directional")
    reason = override.get("reason") or (config.get("default_reasons") or {}).get(cls, "")
    return ReactionSemantics(flux_id, cls, str(evidence), str(reason), reaction is not None)


# ── 代谢物角色与网络化学计量同步 ────────────────────────────────────────────────────────────────────────────────

def _metabolite_ids(name: str) -> frozenset[str]:
    """名称解析到的全部标准 ID（HMDB / KEGG / ChEBI）；比较时按交集判断是否同一代谢物。"""
    ids = lookup_metabolite_ids(str(name)) or {}
    return frozenset(f"{field}:{ids[field]}" for field in ("hmdb_id", "kegg_compound", "chebi_id") if ids.get(field))


def sync_metabolite_roles(mapping: Mapping[str, Any],
                          reactions: Mapping[str, NetworkReaction] | None = None,
                          registered: Iterable[str] | None = None) -> tuple[dict[str, Any], dict[str, str]]:
    """让代谢物映射表的底物 / 产物与网络化学计量一致（组织侧），其余名称移入 related。

    - 网络中的底物 / 产物换成可解析到标准 ID 的标准名称，并保留映射表中解析到同一 ID 的别名；
    - 原先写在底物 / 产物中、但不属于该反应组织侧底物 / 产物的名称（含无法解析 ID 的统称）移入 related；
    - 转运反应（胞质 ⇌ 线粒体）两侧是同一个池，底物 / 产物清空；
    - 网络中有、映射表中没有、且已在 ``registered``（反应注册表）中的组织反应，新建条目。
    返回（新映射, {反应: 变化说明}）。只处理碳网络反应。
    """
    reactions = network_reactions() if reactions is None else reactions
    registered = set(registered) if registered is not None else None
    out = {key: dict(value) if isinstance(value, dict) else value for key, value in mapping.items()}
    changes: dict[str, str] = {}
    for flux_id, reaction in reactions.items():
        if reaction.network != "carbon" or reaction.reaction_class in ("pseudo", "extracellular"):
            continue
        entry = out.get(flux_id)
        if entry is None:
            if registered is not None and flux_id not in registered:
                continue
            if not (reaction.tissue_substrate_names or reaction.tissue_product_names):
                continue
            entry = {"substrates": [], "products": [], "related": []}
        if not isinstance(entry, dict):
            continue
        targets = {"substrates": reaction.tissue_substrate_names, "products": reaction.tissue_product_names}
        keys = {role: set().union(*(_metabolite_ids(n) for n in names)) for role, names in targets.items()}
        old_names = [str(n) for role in ("substrates", "products", "related") for n in entry.get(role, []) or []]
        new = {"substrates": [], "products": [], "related": []}
        for role in ("substrates", "products"):
            for name in targets[role]:
                if _metabolite_ids(name):
                    new[role].append(name)
        for name in old_names:
            ids = _metabolite_ids(name)
            role = next((r for r in ("substrates", "products") if ids & keys[r]), "related")
            if name not in new[role]:
                new[role].append(name)
        new = {role: list(dict.fromkeys(values)) for role, values in new.items()}
        old = {role: list(entry.get(role, []) or []) for role in ("substrates", "products", "related")}
        if new != old:
            changes[flux_id] = ("new entry" if flux_id not in mapping else
                                f"substrates {old['substrates']} -> {new['substrates']}; "
                                f"products {old['products']} -> {new['products']}")
            out[flux_id] = {**{k: v for k, v in entry.items() if k not in new}, **new} if flux_id in mapping else new
    return out, changes


def role_agreement(mapping: Mapping[str, Any], reactions: Mapping[str, NetworkReaction] | None = None) -> dict[str, str]:
    """``{反应: 状态}``：映射表底物 / 产物与网络组织侧底物 / 产物是否一致（consistent / reversed / partial / ...）。"""
    reactions = network_reactions() if reactions is None else reactions
    status: dict[str, str] = {}
    for flux_id, entry in mapping.items():
        reaction = reactions.get(flux_id)
        if reaction is None or reaction.network != "carbon" or not isinstance(entry, dict):
            continue
        net = {r: [_metabolite_ids(n) for n in names]
               for r, names in (("substrates", reaction.tissue_substrate_names),
                                ("products", reaction.tissue_product_names))}
        yml = {r: [ids for n in entry.get(r, []) or [] if (ids := _metabolite_ids(n))] for r in ("substrates", "products")}
        hits = lambda mine, theirs: [ids for ids in mine if any(ids & other for other in theirs)]
        if hits(yml["substrates"], net["products"]) or hits(yml["products"], net["substrates"]):
            status[flux_id] = "reversed"
        elif (len(hits(yml["substrates"], net["substrates"])) < len(yml["substrates"])
              or len(hits(yml["products"], net["products"])) < len(yml["products"])):
            status[flux_id] = "extra_roles"
        else:
            status[flux_id] = "consistent"
    return status


# ── 生成 ────────────────────────────────────────────────────────────────────────────────────────────────────

def _sha256_dir(folder: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(Path(folder).glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _git_commit(folder: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(folder), "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def write_table(rows: list[dict], out_dir: Path) -> None:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / TABLE_FILE).write_bytes(buffer.getvalue().encode("utf-8"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="从 MFA 网络源码生成反应注释表")
    parser.add_argument("--carbon", type=Path, required=True, help="碳网络模块目录")
    parser.add_argument("--nitrogen", type=Path, help="氮网络模块目录（可选）")
    parser.add_argument("--standard-names", type=Path, required=True,
                        help="model_metabolite_to_standard_name_dict.py")
    parser.add_argument("--model-id", required=True, help="写入来源说明的模型 ID，如 dio_multi_tissue")
    parser.add_argument("--out", type=Path, default=DATA_DIR)
    args = parser.parse_args(argv)

    dirs = {"carbon": args.carbon, **({"nitrogen": args.nitrogen} if args.nitrogen else {})}
    rows, report = build_table(dirs, read_standard_names(args.standard_names))
    write_table(rows, args.out)
    sources = {
        "generated_by": "python -m omicsweave.integration.mappings.network_annotation",
        "generated_on": date.today().isoformat(),
        "model_id": args.model_id,
        "networks": {name: {"folder": Path(folder).name, "git_commit": _git_commit(Path(folder)),
                            "sha256_of_modules": _sha256_dir(Path(folder))} for name, folder in dirs.items()},
        "standard_names": {"file": args.standard_names.name, "sha256": hashlib.sha256(
            args.standard_names.read_bytes()).hexdigest()},
        "note": "网络源码只做 ast 静态解析，不执行；f-string 生成的反应（BIOMASS 等）不在表中，按 ID 规则判为虚拟反应",
        "summary": report,
    }
    (args.out / SOURCES_FILE).write_text(json.dumps(sources, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
