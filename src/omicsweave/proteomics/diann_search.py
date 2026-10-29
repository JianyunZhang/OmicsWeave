"""Run DIA-NN 1.9 library-free search on Thermo .raw files.

The statistical pipeline still starts from a DIA-NN report. This module only
fills that report in. It does not read spectra in-process: on Windows, DIA-NN
opens ``.raw`` through Thermo MSFileReader.

``--force`` on the parent pipeline re-runs this step only when search itself
was requested. A matching provenance file skips the executable.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Callable, Sequence

from omicsweave.core.provenance import sha256_file

from .experiment_config import ProteomicsConfig, get_sample_to_mouse_id_map


SEARCH_DIR_NAME = "00_diann_search"
PROVENANCE_NAME = "search_provenance.json"
# Stay under the Windows CreateProcess limit (~32767) after quoting.
_COMMAND_LINE_LIMIT = 30000

MET_OXIDATION = "UniMod:35,15.994915,M"
NTERM_ACETYLATION = "UniMod:1,42.010565,*n"

_VERSION_RE = re.compile(
    r"DIA-NN(?:\s+\(version)?\s+v?(\d+\.\d+(?:\.\d+)?)",
    re.IGNORECASE,
)


class DiannSearchError(RuntimeError):
    """DIA-NN did not produce a report that can enter differential analysis."""


def search_output_dir(output_dir: Path, tissue: str) -> Path:
    _require_tissue_name(tissue)
    return Path(output_dir) / SEARCH_DIR_NAME / tissue


def search_report_path(output_dir: Path, tissue: str) -> Path:
    return search_output_dir(output_dir, tissue) / f"report_{tissue}.tsv"


def list_raw_files(raw_dir: Path) -> list[Path]:
    """List ``*.raw`` files in one directory, without descending into backups."""

    raw_dir = Path(raw_dir)
    if not raw_dir.is_dir():
        raise FileNotFoundError(f"Raw directory not found: {raw_dir}")
    files = [
        path for path in raw_dir.iterdir()
        if path.is_file() and path.suffix.lower() == ".raw"
    ]
    if not files:
        raise FileNotFoundError(f"No .raw files in {raw_dir}")
    return sorted(files, key=lambda path: path.name.casefold())


def build_diann_command(
    *,
    diann_exe: Path,
    raw_files: Sequence[Path],
    fasta: Path,
    report_path: Path,
    threads: int,
    qvalue: float,
    reanalyse: bool,
    smart_profiling: bool,
) -> list[str]:
    """Argv for DIA-NN 1.9 library-free search.

    Fixed cysteine carbamidomethylation is ``--unimod4``. Variable modifications
    are methionine oxidation and protein N-terminal acetylation, at most one
    per peptide. ``--relaxed-prot-inf`` is intentionally absent so protein-group
    FDR stays at the same 1% threshold as the precursor q-value.
    """

    command = [str(diann_exe)]
    for raw_file in raw_files:
        command.extend(["--f", str(raw_file)])
    command.extend([
        "--fasta", str(fasta),
        "--fasta-search",
        "--predictor",
        "--gen-spec-lib",
        "--met-excision",
        "--cut", "K*,R*",
        "--missed-cleavages", "1",
        "--min-pep-len", "7",
        "--max-pep-len", "30",
        "--min-pr-charge", "1",
        "--max-pr-charge", "4",
        "--unimod4",
        "--var-mods", "1",
        "--var-mod", MET_OXIDATION,
        "--var-mod", NTERM_ACETYLATION,
        "--qvalue", _format_qvalue(qvalue),
        "--matrices",
    ])
    if reanalyse:
        command.append("--reanalyse")
    if smart_profiling:
        command.append("--smart-profiling")
    command.extend([
        "--threads", str(int(threads)),
        "--verbose", "1",
        "--out", str(report_path),
    ])
    return command


def validate_raw_sample_names(
    raw_files: Sequence[Path],
    mouse_groups: dict[str, str],
) -> list[dict[str, str]]:
    """Map each raw file to an animal ID before a search is launched.

    DIA-NN writes the filename stem into the report ``Run`` column. The same
    stem is parsed here, then matched to ``groups.yaml`` by longest suffix.
    """

    from .run_proteomics_pipeline import _extract_sample_name

    rows: list[dict[str, str]] = []
    seen: dict[str, str] = {}
    problems: list[str] = []
    for raw_file in raw_files:
        sample = _extract_sample_name(raw_file.stem)
        mouse_id = get_sample_to_mouse_id_map([sample], mouse_groups).get(sample)
        if not mouse_id:
            problems.append(
                f"{raw_file.name} parsed as {sample!r}, which does not end with "
                "an animal ID in groups.yaml"
            )
            continue
        previous = seen.get(sample)
        if previous is not None:
            problems.append(
                f"{raw_file.name} and {previous} both parse to sample {sample!r}"
            )
            continue
        seen[sample] = raw_file.name
        rows.append({
            "raw_file": raw_file.name,
            "sample": sample,
            "mouse_id": mouse_id,
        })
    if problems:
        raise ValueError(
            "Refusing to start DIA-NN; raw filenames do not match groups.yaml:\n"
            + "\n".join(f"  - {problem}" for problem in problems)
        )
    return rows


def extract_diann_version(log_text: str) -> str | None:
    """Read the DIA-NN version from the first lines of its log."""

    head = "\n".join(log_text.splitlines()[:40])
    match = _VERSION_RE.search(head)
    if match is None:
        return None
    return match.group(1)


def version_matches(found: str, expected: str) -> bool:
    """Accept ``1.9`` and ``1.9.x``, and reject ``1.8``, ``1.90``, and ``2.x``."""

    return found == expected or found.startswith(expected + ".")


def run_tissue_search(
    tissue: str,
    cfg: ProteomicsConfig,
    *,
    force: bool = False,
    runner: Callable[..., Any] | None = None,
) -> Path:
    """Search one tissue, or return the existing report when provenance matches.

    ``runner`` is ``subprocess.run`` in production. Tests pass a fake executable
    that must accept the same keyword arguments and return an object with
    ``returncode``, ``stdout``, and ``stderr``.
    """

    if not cfg.search_enabled:
        raise ValueError(
            "DIA-NN search is disabled. Set search.enabled: true and provide "
            "search.diann_exe, search.fasta, and paths.raw_dirs."
        )
    if cfg.diann_exe is None or cfg.search_fasta is None:
        raise ValueError("search.diann_exe and search.fasta are required")
    raw_dir = cfg.raw_dirs.get(tissue)
    if raw_dir is None:
        known = ", ".join(sorted(cfg.raw_dirs)) or "(none)"
        raise ValueError(
            f"[{tissue}] has no paths.raw_dirs entry. Configured tissues: {known}."
        )

    raw_files = list_raw_files(raw_dir)
    samples = validate_raw_sample_names(raw_files, cfg.mouse_groups)
    out_dir = search_output_dir(cfg.output_dir, tissue)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = search_report_path(cfg.output_dir, tissue)
    provenance_path = out_dir / PROVENANCE_NAME
    fingerprint = _fingerprint(cfg, raw_files)

    if not force and _can_skip(report_path, provenance_path, fingerprint, cfg):
        print(
            f"  [{tissue}] DIA-NN provenance matches {len(raw_files)} raw files; "
            "skipping search.",
            flush=True,
        )
        return report_path

    command = build_diann_command(
        diann_exe=cfg.diann_exe,
        raw_files=raw_files,
        fasta=cfg.search_fasta,
        report_path=report_path,
        threads=cfg.search_threads,
        qvalue=cfg.search_qvalue,
        reanalyse=cfg.search_reanalyse,
        smart_profiling=cfg.search_smart_profiling,
    )
    _reject_long_command(tissue, command)
    cmd_path = out_dir / "diann_cmd.txt"
    cmd_path.write_text(
        json.dumps(command, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    # Drop a successful provenance before a re-run so a failed search cannot
    # be skipped on the next invocation while a partial report is left behind.
    if provenance_path.exists():
        provenance_path.unlink()

    print(
        f"  [{tissue}] running DIA-NN on {len(raw_files)} raw files "
        f"-> {report_path}",
        flush=True,
    )
    for row in samples:
        print(f"    {row['raw_file']} -> {row['sample']} ({row['mouse_id']})", flush=True)

    run = runner or subprocess.run
    completed = run(
        command,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    log_text = f"{completed.stdout or ''}{completed.stderr or ''}"
    log_path = out_dir / "diann.log"
    log_path.write_text(log_text, encoding="utf-8")
    if int(completed.returncode) != 0:
        raise DiannSearchError(
            f"[{tissue}] DIA-NN exited with code {completed.returncode}. "
            f"See {log_path}."
        )

    found_version = extract_diann_version(log_text)
    if found_version is None or not version_matches(found_version, cfg.diann_expected_version):
        raise DiannSearchError(
            f"[{tissue}] DIA-NN version {found_version!r} does not match "
            f"expected {cfg.diann_expected_version!r}. "
            "The report was not accepted for differential analysis. "
            f"See {log_path}."
        )
    if not report_path.is_file() or report_path.stat().st_size <= 0:
        raise DiannSearchError(
            f"[{tissue}] DIA-NN finished without a report at {report_path}."
        )

    provenance = {
        "engine": cfg.search_engine,
        "diann_version": found_version,
        "fingerprint": fingerprint,
        "samples": samples,
        "command": command,
        "report": str(report_path),
    }
    provenance_path.write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"  [{tissue}] DIA-NN {found_version} wrote {report_path.name}", flush=True)
    return report_path


def _require_tissue_name(tissue: str) -> None:
    if not re.fullmatch(r"[A-Za-z0-9_]+", tissue):
        raise ValueError(
            f"Tissue name {tissue!r} cannot be used as a DIA-NN output directory"
        )


def _format_qvalue(qvalue: float) -> str:
    return format(float(qvalue), "g")


def _file_stamp(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _fingerprint(cfg: ProteomicsConfig, raw_files: Sequence[Path]) -> dict[str, Any]:
    """Identity of the inputs that change a DIA-NN result.

    Thread count is omitted: it does not change the report, and including it
    would repeat a multi-hour search after a hardware change. The FASTA is
    hashed; raw files are stamped by size and modification time.
    """

    assert cfg.search_fasta is not None
    fasta_stamp = _file_stamp(cfg.search_fasta)
    fasta_stamp["sha256"] = sha256_file(cfg.search_fasta)
    return {
        "expected_version": cfg.diann_expected_version,
        "qvalue": _format_qvalue(cfg.search_qvalue),
        "reanalyse": bool(cfg.search_reanalyse),
        "smart_profiling": bool(cfg.search_smart_profiling),
        "cut": "K*,R*",
        "missed_cleavages": 1,
        "min_pep_len": 7,
        "max_pep_len": 30,
        "min_pr_charge": 1,
        "max_pr_charge": 4,
        "fixed_modification": "unimod4",
        "variable_modifications": [MET_OXIDATION, NTERM_ACETYLATION],
        "max_variable_modifications": 1,
        "fasta": fasta_stamp,
        "raw_files": [_file_stamp(path) for path in raw_files],
    }


def _can_skip(
    report_path: Path,
    provenance_path: Path,
    fingerprint: dict[str, Any],
    cfg: ProteomicsConfig,
) -> bool:
    if not report_path.is_file() or report_path.stat().st_size <= 0:
        return False
    if not provenance_path.is_file():
        return False
    try:
        stored = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    if not isinstance(stored, dict):
        return False
    found = str(stored.get("diann_version", ""))
    if not version_matches(found, cfg.diann_expected_version):
        return False
    return stored.get("fingerprint") == fingerprint


def _reject_long_command(tissue: str, command: Sequence[str]) -> None:
    if os.name != "nt":
        return
    length = sum(len(part) + 3 for part in command)
    if length > _COMMAND_LINE_LIMIT:
        raise ValueError(
            f"[{tissue}] DIA-NN command line is {length} characters, above the "
            f"Windows limit ({_COMMAND_LINE_LIMIT}). This tissue has too many "
            ".raw files for one DIA-NN invocation."
        )
