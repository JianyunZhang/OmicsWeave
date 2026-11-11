"""生成小鼠基因 / 蛋白 ID 注册表（``data/`` 下的三张表与 ``SOURCES.json``）。

注册表不再手写：基因符号、Ensembl ID、MGI ID 来自 Ensembl GRCm39 cDNA FASTA 的序列头，
UniProt 登记号、条目名、EC 号、NCBI Gene ID 与基因别名来自 UniProtKB/Swiss-Prot 小鼠条目，
两者按 MGI 登记号对齐，规范符号一律取 MGI 符号（例如 UniProt 的 ``Pc`` 记为 ``Pcx`` 的别名）。

用法::

    python -m omicsweave.core.idmap.build_registry \\
        --uniprot uniprot_mouse_sprot.tsv.gz --uniprot-release 2026_04 --uniprot-release-date 2026-11-04 \\
        --ensembl-cdna external/transcriptomics/reference/Mus_musculus.GRCm39.cdna.all.fa.gz

UniProt 表的获取方式（字段顺序不限，列名须为 UniProt 的默认表头）::

    https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv
        &fields=accession,id,gene_primary,gene_synonym,protein_name,ec,xref_mgi,xref_geneid,xref_ensembl
        &query=(organism_id:10090) AND (reviewed:true)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable

DATA_DIR = Path(__file__).resolve().parent / "data"
GENE_FILE = "mouse_gene_registry.tsv.gz"
PROTEIN_FILE = "mouse_protein_registry.tsv.gz"
ALIAS_FILE = "mouse_gene_aliases.tsv.gz"
SOURCES_FILE = "SOURCES.json"
UNIPROT_QUERY = ("https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv"
                 "&fields=accession,id,gene_primary,gene_synonym,protein_name,ec,xref_mgi,xref_geneid,xref_ensembl"
                 "&query=(organism_id:10090) AND (reviewed:true)")

# 人工维护的别名：蛋白俗名、映射表里出现过的旧写法，以及不是基因的条目（目标为空）。
# 构建时逐条校验：键若本身是 MGI 正式符号则丢弃，目标若不是正式符号则丢弃，结果写入报告。
CURATED_ALIASES: dict[str, str | None] = {
    "Pdh1": "Pdha1", "Pdh1a": "Pdha1", "Ogdh2": "Ogdh", "Suc1": "Suclg1", "Suc2": "Sucla2",
    "Gpk": "Gk", "Gyk": "Gk", "Acats1": "Acat1", "Tkt1": "Tkt", "Rpi1": "Rpia", "Ugpdh": "Ugp2",
    "Gpdt1": "Gpam", "Bpgm1": "Bpgm",
    "Atgl": "Pnpla2", "Hsl": "Lipe", "Mct1": "Slc16a1", "Mct2": "Slc16a7", "Mct4": "Slc16a3",
    "Fatp1": "Slc27a1", "Fatp2": "Slc27a2", "Fatp3": "Slc27a3", "Fatp4": "Slc27a4", "Cd147": "Bsg",
    "Nrf2": "Nfe2l2", "Chrebp": "Mlxipl", "Mpc": None,
    "Ppar1a": "Ppara", "Ppar1b": "Ppard", "Ppar1g": "Pparg", "Ppar2": "Ppard", "Ppar2g": "Pparg",
    "Pparg1a": "Ppargc1a", "Pparg1b": "Ppargc1b", "Pparg2": "Pparg", "Srebp1": "Srebf1", "Srebp1c": "Srebf1",
    "Lxr1a": "Nr1h3", "Fx1r": "Nr1h4", "Ampk1a": "Prkaa1", "Ampk1b": "Prkaa2",
    "Gcn2": "Eif2ak4", "Phd1": "Egln2", "Phd2": "Egln1", "Nampt1": "Nampt",
    "Oaa": None, "Nadh": None, "Nadu": None, "Biotin": None, "Carnitine": None, "Pp1": None, "Ppk1": None,
    "Purat": None,
}

_ENSEMBL_HEADER = re.compile(
    r"gene:(?P<gene>ENSMUSG\d+)(?:\.\d+)?\s+gene_biotype:(?P<biotype>\S+).*?"
    r"(?:gene_symbol:(?P<symbol>\S+))?(?:\s+description:(?P<desc>.*))?$")
_CHROMOSOME = re.compile(r"chromosome:[^:]+:([^:]+):")
_MGI = re.compile(r"\[Source:MGI Symbol;Acc:(MGI:\d+)\]")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _open_text(path: Path) -> io.TextIOBase:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_ensembl_genes(path: Path) -> dict[str, dict[str, str]]:
    """Ensembl cDNA FASTA → ``{Ensembl 基因 ID: {symbol, mgi_id, gene_type, chromosome, gene_name}}``。"""
    genes: dict[str, dict[str, str]] = {}
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = _ENSEMBL_HEADER.search(line.rstrip("\n"))
            if not match or match["gene"] in genes:
                continue
            desc = match["desc"] or ""
            mgi = _MGI.search(desc)
            chromosome = _CHROMOSOME.search(line)
            genes[match["gene"]] = {
                "symbol": match["symbol"] or "",
                "mgi_id": mgi.group(1) if mgi else "",
                "gene_type": match["biotype"],
                "chromosome": chromosome.group(1) if chromosome else "",
                "gene_name": _MGI.sub("", desc).strip(),
            }
    return genes


def _split(value: str, sep: str = ";") -> list[str]:
    return [piece.strip() for piece in (value or "").split(sep) if piece.strip()]


def read_uniprot(path: Path) -> list[dict[str, object]]:
    rows = []
    with _open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({
                "uniprot_id": row["Entry"].strip(),
                "entry_name": row["Entry Name"].strip(),
                "gene_primary": (row.get("Gene Names (primary)") or "").strip(),
                "gene_synonyms": (row.get("Gene Names (synonym)") or "").split(),
                "protein_name": re.split(r" \(", (row.get("Protein names") or "").strip(), maxsplit=1)[0],
                "ec_number": "; ".join(_split(row.get("EC number") or "")),
                "mgi_ids": _split(row.get("MGI") or ""),
                "ncbi_gene_ids": _split(row.get("GeneID") or ""),
            })
    return rows


def build(uniprot_path: Path, ensembl_path: Path) -> tuple[list[dict], list[dict], list[dict], dict]:
    ensembl = read_ensembl_genes(ensembl_path)
    uniprot = read_uniprot(uniprot_path)

    # ── 基因：规范符号 = MGI 符号；一个符号可对应多个 Ensembl 基因（极少数）
    genes: dict[str, dict[str, object]] = {}
    mgi_to_symbol: dict[str, str] = {}
    for gene_id, info in sorted(ensembl.items()):
        symbol = str(info["symbol"])
        if not symbol:
            continue
        entry = genes.setdefault(symbol, {
            "symbol": symbol, "mgi_id": info["mgi_id"], "ensembl_gene_ids": [], "ncbi_gene_ids": set(),
            "gene_type": info["gene_type"], "chromosome": info["chromosome"], "gene_name": info["gene_name"],
            "source": "ensembl",
        })
        entry["ensembl_gene_ids"].append(gene_id)
        if info["mgi_id"]:
            mgi_to_symbol.setdefault(str(info["mgi_id"]), symbol)

    # ── 蛋白：按 MGI 登记号对齐到规范符号；对不上时退回 UniProt 主基因名
    proteins: list[dict[str, object]] = []
    unmatched_uniprot = 0
    for row in uniprot:
        symbol = next((mgi_to_symbol[m] for m in row["mgi_ids"] if m in mgi_to_symbol), "")
        if not symbol and row["gene_primary"] in genes:
            symbol = str(row["gene_primary"])
        if not symbol and row["gene_primary"]:
            symbol = str(row["gene_primary"])
            unmatched_uniprot += 1
            genes.setdefault(symbol, {
                "symbol": symbol, "mgi_id": (row["mgi_ids"] or [""])[0], "ensembl_gene_ids": [],
                "ncbi_gene_ids": set(), "gene_type": "protein_coding", "chromosome": "",
                "gene_name": row["protein_name"], "source": "uniprot",
            })
        if not symbol:
            continue
        genes[symbol]["ncbi_gene_ids"].update(row["ncbi_gene_ids"])
        proteins.append({
            "uniprot_id": row["uniprot_id"], "entry_name": row["entry_name"], "symbol": symbol,
            "uniprot_gene_name": row["gene_primary"], "ec_number": row["ec_number"],
            "protein_name": row["protein_name"],
        })

    official = set(genes)
    official_fold: dict[str, set[str]] = defaultdict(set)
    for symbol in official:
        official_fold[symbol.casefold()].add(symbol)

    # ── 别名：UniProt 主基因名（与 MGI 不同者）和同义名 + 人工别名；歧义者丢弃
    candidates: dict[str, set[str]] = defaultdict(set)
    alias_source: dict[str, str] = {}
    for row, protein in _pair_rows(uniprot, proteins):
        names = [row["gene_primary"], *row["gene_synonyms"]]
        for name in names:
            name = str(name)
            if name and name != protein["symbol"]:
                candidates[name].add(str(protein["symbol"]))
                alias_source.setdefault(name, "uniprot")
    dropped: dict[str, list[str]] = defaultdict(list)
    aliases: list[dict[str, str]] = []
    for alias, targets in sorted(candidates.items()):
        if alias in official or official_fold.get(alias.casefold(), set()) - targets:
            dropped["alias_is_official_symbol"].append(alias)
        elif len(targets) > 1:
            dropped["ambiguous_alias"].append(f"{alias}->{'|'.join(sorted(targets))}")
        else:
            aliases.append({"alias": alias, "symbol": next(iter(targets)), "source": alias_source[alias]})
    uniprot_aliases = {a["alias"] for a in aliases}
    for alias, target in sorted(CURATED_ALIASES.items()):
        if alias in official:
            dropped["curated_key_is_official_symbol"].append(alias)
        elif target is not None and target not in official:
            dropped["curated_target_not_official"].append(f"{alias}->{target}")
        elif alias in uniprot_aliases:
            continue
        else:
            aliases.append({"alias": alias, "symbol": target or "", "source": "curated"})
    aliases.sort(key=lambda a: a["alias"].casefold())

    gene_rows = [{
        "symbol": g["symbol"], "mgi_id": g["mgi_id"], "ensembl_gene_ids": ";".join(g["ensembl_gene_ids"]),
        "ncbi_gene_ids": ";".join(sorted(g["ncbi_gene_ids"], key=int)), "gene_type": g["gene_type"],
        "chromosome": g["chromosome"], "gene_name": g["gene_name"], "source": g["source"],
    } for g in sorted(genes.values(), key=lambda g: str(g["symbol"]))]
    proteins.sort(key=lambda p: (str(p["symbol"]), str(p["uniprot_id"])))

    collisions = sorted("|".join(sorted(v)) for v in official_fold.values() if len(v) > 1)
    report = {
        "n_genes": len(gene_rows), "n_genes_from_ensembl": sum(r["source"] == "ensembl" for r in gene_rows),
        "n_proteins": len(proteins), "n_uniprot_entries": len(uniprot),
        "n_uniprot_without_gene_name": len(uniprot) - len(proteins),
        "n_uniprot_symbol_not_in_ensembl": unmatched_uniprot,
        "n_symbols_with_several_proteins": sum(1 for n in _count(p["symbol"] for p in proteins).values() if n > 1),
        "n_aliases": len(aliases), "dropped_aliases": {k: len(v) for k, v in dropped.items()},
        "dropped_curated_aliases": {k: v for k, v in dropped.items() if k.startswith("curated")},
        "casefold_collisions": collisions,
    }
    return gene_rows, proteins, aliases, report


def _count(values: Iterable[object]) -> dict[object, int]:
    counts: dict[object, int] = defaultdict(int)
    for value in values:
        counts[value] += 1
    return counts


def _pair_rows(uniprot: list[dict], proteins: list[dict]) -> Iterable[tuple[dict, dict]]:
    by_accession = {str(p["uniprot_id"]): p for p in proteins}
    for row in uniprot:
        protein = by_accession.get(str(row["uniprot_id"]))
        if protein is not None:
            yield row, protein


def _write_tsv_gz(path: Path, rows: list[dict], columns: list[str]) -> None:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    with path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
        handle.write(buffer.getvalue().encode("utf-8"))


def write(out_dir: Path, gene_rows: list[dict], proteins: list[dict], aliases: list[dict], sources: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv_gz(out_dir / GENE_FILE, gene_rows,
                  ["symbol", "mgi_id", "ensembl_gene_ids", "ncbi_gene_ids", "gene_type", "chromosome", "gene_name",
                   "source"])
    _write_tsv_gz(out_dir / PROTEIN_FILE, proteins,
                  ["uniprot_id", "entry_name", "symbol", "uniprot_gene_name", "ec_number", "protein_name"])
    _write_tsv_gz(out_dir / ALIAS_FILE, aliases, ["alias", "symbol", "source"])
    (out_dir / SOURCES_FILE).write_text(json.dumps(sources, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--uniprot", type=Path, required=True, help="UniProt 小鼠 Swiss-Prot TSV（可为 .gz）")
    parser.add_argument("--uniprot-release", required=True, help="UniProt 版本号，如 2026_04")
    parser.add_argument("--uniprot-release-date", required=True, help="UniProt 发布日期，如 2026-11-04")
    parser.add_argument("--ensembl-cdna", type=Path, required=True, help="Ensembl GRCm39 cDNA FASTA（可为 .gz）")
    parser.add_argument("--out", type=Path, default=DATA_DIR, help="输出目录（默认为本包的 data/）")
    args = parser.parse_args(argv)

    gene_rows, proteins, aliases, report = build(args.uniprot, args.ensembl_cdna)
    sources = {
        "generated_by": "python -m omicsweave.core.idmap.build_registry",
        "generated_on": date.today().isoformat(),
        "uniprot": {
            "release": args.uniprot_release, "release_date": args.uniprot_release_date, "query": UNIPROT_QUERY,
            "file": args.uniprot.name, "sha256": _sha256(args.uniprot),
        },
        "ensembl": {
            "file": args.ensembl_cdna.name, "assembly": "GRCm39", "sha256": _sha256(args.ensembl_cdna),
            "note": "Ensembl 发布号未写入 FASTA；以文件名与 sha256 追溯",
        },
        "join": "UniProt 条目按 MGI 登记号对齐到 Ensembl 序列头中的 MGI 符号；对不上时用 UniProt 主基因名",
        "summary": report,
    }
    write(args.out, gene_rows, proteins, aliases, sources)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
