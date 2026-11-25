"""从固定结果包汇总科学质量证据，诊断本身只读。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from math import isfinite
from pathlib import Path
from typing import Literal

import pyarrow.parquet as pq

from .reliability import BLOCKED

Category = Literal["engineering_fixable", "sensitivity_explorable", "statistical_guardrail"]
_ACTIONS = {
    "engineering_fixable": ("inspect_source", "preview_config_patch"),
    "sensitivity_explorable": ("preview_qc_update", "isolated_trial", "diff_runs"),
    "statistical_guardrail": ("report_evidence_limit", "review_design"),
}
_CONTEXT = ("omics", "tissue", "contrast_id", "n_numerator", "n_denominator")


@dataclass(frozen=True)
class Finding:
    category: Category
    code: str
    source: str
    message: str
    context: dict = field(default_factory=dict)

    def payload(self) -> dict:
        return {"category": self.category, "code": self.code, "source": self.source,
                "message": self.message, "context": self.context,
                "allowed_actions": list(_ACTIONS[self.category])}


def _rows(path: Path, columns: tuple[str, ...], required: set[str]):
    with pq.ParquetFile(path) as parquet:
        names = set(parquet.schema_arrow.names)
        if not required <= names:
            return
        for batch in parquet.iter_batches(batch_size=2048, columns=[key for key in columns if key in names]):
            yield from batch.to_pylist()


def _finite(value) -> bool:
    return isinstance(value, (int, float)) and isfinite(value)


def inspect_evidence(tables: dict[str, Path], warnings: list[str], alpha: float) -> tuple[list[Finding], dict]:
    findings = []
    availability = {"sample_quality": 0, "contrast_reliability": 0, "rank_floor": 0,
                    "decoy_strata": 0, "decoy_unassessed_strata": 0, "validation_warnings": len(warnings)}
    for warning in warnings:
        limited = "秩检验在以下分层达不到" in warning or "生物学重复不足" in warning
        findings.append(Finding("statistical_guardrail" if limited else "engineering_fixable",
                                "VALIDATION_POWER_LIMIT" if limited else "VALIDATION_WARNING",
                                "index.json#validation_warnings", warning[:700]))

    for key, path in sorted(tables.items()):
        stem = path.stem
        if stem == "contrast_reliability":
            fields = (*_CONTEXT, "status", "reason", "source", "n_usable_numerator", "n_usable_denominator")
            for row in _rows(path, fields, {"status"}):
                availability["contrast_reliability"] += 1
                if row["status"] in BLOCKED:
                    reason = str(row.pop("reason", "") or "比较受样本质控限制")[:700]
                    findings.append(Finding("statistical_guardrail", row["status"].upper(), key, reason, row))
        elif stem == "sample_quality":
            flagged = Counter()
            for row in _rows(path, ("tissue", "flag", "metric"), {"flag"}):
                availability["sample_quality"] += 1
                if row["flag"] in {"exclude", "flag"}:
                    flagged[(row.get("tissue"), row.get("metric"), row["flag"])] += 1
            for (tissue, metric, flag), count in flagged.items():
                findings.append(Finding("sensitivity_explorable", "SAMPLE_QC_FLAG", key,
                                        "复核已记录的样本质控标记，再进行有配对基线的敏感性试验",
                                        {"tissue": tissue, "metric": metric, "flag": flag, "flagged_rows": count}))
        elif stem in {"differential", "differential_by_tracer"}:
            strata = Counter()
            context = (*_CONTEXT, "tracer", "method")
            for row in _rows(path, (*context, "min_achievable_p"), {"min_achievable_p"}):
                p_floor = row["min_achievable_p"]
                if not _finite(p_floor):
                    continue
                availability["rank_floor"] += 1
                if p_floor >= alpha:
                    strata[tuple(row.get(name) for name in context) + (p_floor,)] += 1
            for values, count in strata.items():
                details = dict(zip(context, values[:-1]))
                details.update(alpha=alpha, min_achievable_p=values[-1], affected_rows=count)
                findings.append(Finding("statistical_guardrail", "RANK_TEST_POWER_LIMIT", key,
                                        "已登记的最小可达 p 值达到或高于诊断阈值", details))
        elif stem == "assessment_calibration":
            fields = ("method", "tissue", "contrast_id", "metric", "target", "decoy_expected", "decoy_sd",
                      "n_replicates", "passed", "detail")
            for row in _rows(path, fields, {"method"}):
                if row["method"] == "invariance" and row.get("passed") is False:
                    findings.append(Finding("statistical_guardrail", "CALIBRATION_INVARIANCE_FAILED", key,
                                            "校准记录显示不变性检查未通过",
                                            {name: row.get(name) for name in ("tissue", "contrast_id", "metric", "target")}))
                if row["method"] != "decoy" or row.get("metric") != "L1":
                    continue
                observed, expected, spread = (row.get(name) for name in ("target", "decoy_expected", "decoy_sd"))
                repeats = row.get("n_replicates") or 0
                if repeats < 2 or not all(map(_finite, (observed, expected, spread))):
                    availability["decoy_unassessed_strata"] += 1
                    continue
                availability["decoy_strata"] += 1
                if (observed or expected) and observed <= expected + spread:
                    findings.append(Finding(
                        "statistical_guardrail", "L1_CLOSE_TO_DECOYS", key,
                        "L1 数位于诱饵均值加一个标准差以内。此规则用于描述性筛查",
                        {"tissue": row.get("tissue"), "contrast_id": row.get("contrast_id"),
                         "cell_L1": observed, "decoy_L1_expected": expected,
                         "decoy_L1_sd": spread, "decoy_replicates": repeats},
                    ))
    return findings, availability


def inspect_registry(directories: dict[str, Path]) -> list[Finding]:
    from .refs import KNOWLEDGE_CHECKS, audit

    findings = []
    reported = set()
    for run, directory in directories.items():
        table = audit(run_dir=directory)
        for (check, layer, tissue, issue), rows in table.groupby(["check", "omics", "tissue", "issue"], sort=True):
            if issue == "coverage":
                continue
            key = (check, layer, tissue, issue)
            if check in KNOWLEDGE_CHECKS and key in reported:
                continue
            reported.add(key)
            findings.append(Finding("engineering_fixable", "REFERENCE_AUDIT_FINDING",
                                    "refs audit" if check in KNOWLEDGE_CHECKS else f"refs audit:{run}",
                                    "参考知识表审计发现需要核对的条目",
                                    {"omics": layer, "tissue": tissue, "check": check,
                                     "issue": issue, "affected_rows": len(rows)}))
    return findings
