#!/usr/bin/env python3
"""
跨批次联合建模的回归测试
==========================

本项目目前尚无符合条件的多批次转录组数据，因此这里用**合成数据**完整
模拟 ``runs/<protocol>/<batch>/`` 目录结构，覆盖：

* 平衡设计 → 联合模型运行，且批次效应被正确吸收（不产生假阳性）
* 混杂设计 → 预检拒绝，退回 Meta 路径
* 参考不一致 → 合并被阻断
* 样本 ID 跨批次重复 → 默认拒绝（桥接样本不能当独立重复）
* 基因集不一致 → 取交集并告警

运行：
  conda run -n tutorial-mfa-2025 python -m unittest test_joint_model -v
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from omicsweave.transcriptomics import joint_model


CONSISTENT_PROV = {
    "transcriptome_fingerprint": "ref-abc",
    "kallisto_version": "kallisto 0.48.0",
    "strandedness": "rf-stranded",
    "gene_aggregation": "lengthScaledTPM",
}


def _write_run(root: Path, protocol: str, batch: str,
               counts: pd.DataFrame, metadata: pd.DataFrame,
               prov: dict | None = CONSISTENT_PROV) -> None:
    """在临时目录里造出一个批次运行的产物。"""
    run_dir = root / "runs" / protocol / batch
    (run_dir / "04_de_results").mkdir(parents=True, exist_ok=True)
    counts.to_csv(run_dir / "04_de_results" / "gene_counts_matrix.csv")
    metadata.to_csv(run_dir / "sample_metadata.csv")
    if prov is not None:
        (run_dir / "run_provenance.json").write_text(
            json.dumps(prov), encoding="utf-8"
        )


def _simulate_batch(
    mice: list[str],
    conditions: list[str],
    tissue: str,
    n_genes: int,
    rng: np.random.Generator,
    batch_shift_genes: slice | None = None,
    batch_shift: float = 1.0,
    condition_effect_genes: slice | None = None,
    condition_effect: float = 1.0,
    effect_condition: str | None = None,
    base_counts: np.ndarray | None = None,
):
    """生成一个批次的计数矩阵与元数据。"""
    samples = [f"{m}_{tissue}" for m in mice]
    if base_counts is None:
        base_counts = rng.lognormal(mean=5.5, sigma=0.8, size=(n_genes, 1))
    mat = np.repeat(base_counts, len(samples), axis=1)
    mat = mat * rng.lognormal(0, 0.15, size=(n_genes, len(samples)))

    if batch_shift_genes is not None:
        mat[batch_shift_genes, :] *= batch_shift

    if condition_effect_genes is not None and effect_condition is not None:
        cols = [i for i, c in enumerate(conditions) if c == effect_condition]
        if cols:
            mat[condition_effect_genes, :][:, cols]  # noqa: B018
            mat[condition_effect_genes.start:condition_effect_genes.stop, cols] *= condition_effect

    counts = pd.DataFrame(
        np.random.default_rng(0).poisson(np.clip(mat, 1, None)).astype(int),
        index=[f"Gene{i}" for i in range(n_genes)],
        columns=samples,
    )
    counts.index.name = "gene"
    metadata = pd.DataFrame(
        {
            "mouse_id": mice,
            "tissue": [tissue] * len(mice),
            "condition": conditions,
        },
        index=pd.Index(samples, name="sample"),
    )
    return counts, metadata


class MergeRunsTest(unittest.TestCase):
    def _two_balanced_batches(self, root: Path, n_genes: int = 200):
        rng = np.random.default_rng(5)
        base = rng.lognormal(mean=5.5, sigma=0.8, size=(n_genes, 1))
        c1, m1 = _simulate_batch(
            ["1", "2", "3", "4"], ["hfd", "hfd", "lfd", "lfd"], "liver",
            n_genes, rng, base_counts=base,
        )
        c2, m2 = _simulate_batch(
            ["5", "6", "7", "8"], ["hfd", "hfd", "lfd", "lfd"], "liver",
            n_genes, rng, base_counts=base,
        )
        _write_run(root, "proto", "b1", c1, m1)
        _write_run(root, "proto", "b2", c2, m2)

    def test_merges_and_assigns_batch_from_run_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._two_balanced_batches(root)

            runs = joint_model.discover_run_outputs(root, "proto", ["b1", "b2"])
            counts, metadata, info = joint_model.merge_run_counts(runs)

            self.assertEqual(info["n_batches"], 2)
            self.assertEqual(info["n_samples"], 8)
            self.assertEqual(sorted(metadata["batch"].unique()), ["b1", "b2"])
            # 计数列顺序必须与元数据索引严格对齐
            self.assertEqual(list(counts.columns), list(metadata.index))

    def test_requires_at_least_two_batches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(1)
            c1, m1 = _simulate_batch(["1", "2"], ["hfd", "lfd"], "liver", 50, rng)
            _write_run(root, "proto", "b1", c1, m1)

            runs = joint_model.discover_run_outputs(root, "proto", ["b1", "b2"])
            with self.assertRaisesRegex(ValueError, "至少需要 2 个"):
                joint_model.merge_run_counts(runs)

    def test_duplicate_sample_ids_rejected_by_default(self) -> None:
        """桥接样本是技术重复，不能静默当作独立生物学重复。"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(2)
            # 两个批次用了相同的小鼠编号
            c1, m1 = _simulate_batch(["1", "2", "3", "4"],
                                     ["hfd", "hfd", "lfd", "lfd"], "liver", 50, rng)
            c2, m2 = _simulate_batch(["1", "2", "3", "4"],
                                     ["hfd", "hfd", "lfd", "lfd"], "liver", 50, rng)
            _write_run(root, "proto", "b1", c1, m1)
            _write_run(root, "proto", "b2", c2, m2)

            runs = joint_model.discover_run_outputs(root, "proto", ["b1", "b2"])
            with self.assertRaisesRegex(ValueError, "重复出现"):
                joint_model.merge_run_counts(runs)

    def test_duplicate_suffix_policy_keeps_both(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(3)
            c1, m1 = _simulate_batch(["1", "2", "3", "4"],
                                     ["hfd", "hfd", "lfd", "lfd"], "liver", 50, rng)
            c2, m2 = _simulate_batch(["1", "2", "3", "4"],
                                     ["hfd", "hfd", "lfd", "lfd"], "liver", 50, rng)
            _write_run(root, "proto", "b1", c1, m1)
            _write_run(root, "proto", "b2", c2, m2)

            runs = joint_model.discover_run_outputs(root, "proto", ["b1", "b2"])
            counts, metadata, info = joint_model.merge_run_counts(
                runs, duplicate_policy=joint_model.DUPLICATE_SUFFIX
            )

            self.assertEqual(info["n_samples"], 8)
            self.assertTrue(any("__b1" in s for s in metadata.index))
            self.assertTrue(any("__b2" in s for s in metadata.index))

    def test_gene_intersection_and_overlap_warning(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(4)
            c1, m1 = _simulate_batch(["1", "2", "3", "4"],
                                     ["hfd", "hfd", "lfd", "lfd"], "liver", 200, rng)
            c2, m2 = _simulate_batch(["5", "6", "7", "8"],
                                     ["hfd", "hfd", "lfd", "lfd"], "liver", 200, rng)
            # 批次 2 缺失一半基因 —— 模拟参考版本不同
            c2 = c2.iloc[:100]
            _write_run(root, "proto", "b1", c1, m1)
            _write_run(root, "proto", "b2", c2, m2)

            runs = joint_model.discover_run_outputs(root, "proto", ["b1", "b2"])
            counts, _, info = joint_model.merge_run_counts(runs)

            self.assertEqual(info["n_genes"], 100)
            self.assertIn("gene_overlap_warning", info)


class DesignAssessmentTest(unittest.TestCase):
    def _metadata(self, rows):
        return pd.DataFrame(
            rows,
            index=pd.Index([f"s{i}" for i in range(len(rows))], name="sample"),
        )

    def test_balanced_design_routed_to_joint_model(self) -> None:
        meta = self._metadata([
            {"tissue": "liver", "condition": c, "batch": b}
            for b in ("b1", "b2") for c in ("hfd", "hfd", "lfd", "lfd")
        ])

        report = joint_model.assess_joint_design(meta, [("hfd", "lfd", "hfd_vs_lfd")])

        row = report.iloc[0]
        self.assertTrue(row["identifiable"])
        self.assertEqual(row["route"], joint_model.ROUTE_JOINT)
        self.assertEqual(row["design_formula"], "~batch + condition")
        self.assertEqual(row["n_batches"], 2)

    def test_confounded_design_falls_back_to_meta(self) -> None:
        """批次 1 全是 hfd、批次 2 全是 lfd —— 技术与生物学差异不可分离。"""
        meta = self._metadata(
            [{"tissue": "liver", "condition": "hfd", "batch": "b1"}] * 4
            + [{"tissue": "liver", "condition": "lfd", "batch": "b2"}] * 4
        )

        report = joint_model.assess_joint_design(meta, [("hfd", "lfd", "hfd_vs_lfd")])

        row = report.iloc[0]
        self.assertFalse(row["identifiable"])
        self.assertEqual(row["route"], joint_model.ROUTE_META)
        self.assertIn("混杂", row["reason"])

    def test_too_few_samples_marked_insufficient(self) -> None:
        meta = self._metadata([
            {"tissue": "liver", "condition": "hfd", "batch": "b1"},
            {"tissue": "liver", "condition": "lfd", "batch": "b2"},
        ])

        report = joint_model.assess_joint_design(meta, [("hfd", "lfd", "hfd_vs_lfd")])

        self.assertEqual(report.iloc[0]["route"], joint_model.ROUTE_NONE)

    def test_crosstab_is_recorded_for_inspection(self) -> None:
        meta = self._metadata([
            {"tissue": "liver", "condition": c, "batch": b}
            for b in ("b1", "b2") for c in ("hfd", "hfd", "lfd", "lfd")
        ])

        report = joint_model.assess_joint_design(meta, [("hfd", "lfd", "hfd_vs_lfd")])

        table = report.iloc[0]["batch_condition_table"]
        self.assertIn("b1", table)
        self.assertIn("b2", table)
        self.assertIn("hfd=2", table)

    def test_per_tissue_assessment_is_independent(self) -> None:
        """一个组织混杂不应影响另一个组织的判定。"""
        meta = self._metadata(
            # liver: 平衡
            [{"tissue": "liver", "condition": c, "batch": b}
             for b in ("b1", "b2") for c in ("hfd", "hfd", "lfd", "lfd")]
            # brain: 混杂
            + [{"tissue": "brain", "condition": "hfd", "batch": "b1"}] * 4
            + [{"tissue": "brain", "condition": "lfd", "batch": "b2"}] * 4
        )

        report = joint_model.assess_joint_design(meta, [("hfd", "lfd", "hfd_vs_lfd")])
        by_tissue = report.set_index("tissue")

        self.assertEqual(by_tissue.loc["liver", "route"], joint_model.ROUTE_JOINT)
        self.assertEqual(by_tissue.loc["brain", "route"], joint_model.ROUTE_META)


class JointModelEndToEndTest(unittest.TestCase):
    """
    完整端到端：注入已知的批次效应与条件效应，验证联合模型能把两者分开。

    这是本功能最关键的验证——批次项若没有真正吸收批次效应，
    仅受批次影响的基因就会被误判为差异表达。
    """

    def _build(self, root: Path, confounded: bool = False, n_genes: int = 300):
        rng = np.random.default_rng(17)
        base = rng.lognormal(mean=6.0, sigma=0.7, size=(n_genes, 1))

        def make(mice, conds, batch_factor):
            samples = [f"{m}_liver" for m in mice]
            mat = np.repeat(base, len(samples), axis=1)
            mat = mat * rng.lognormal(0, 0.10, size=(n_genes, len(samples)))
            # 基因 0-39：仅受批次影响（不应被判为差异表达）
            mat[0:40, :] *= batch_factor
            # 基因 40-79：仅受条件影响（hfd 上调，应被检出）
            hfd_cols = [i for i, c in enumerate(conds) if c == "hfd"]
            if hfd_cols:
                mat[np.ix_(range(40, 80), hfd_cols)] *= 4.0
            counts = pd.DataFrame(
                rng.poisson(np.clip(mat, 1, None)).astype(int),
                index=[f"Gene{i}" for i in range(n_genes)],
                columns=samples,
            )
            counts.index.name = "gene"
            meta = pd.DataFrame(
                {"mouse_id": mice, "tissue": ["liver"] * len(mice), "condition": conds},
                index=pd.Index(samples, name="sample"),
            )
            return counts, meta

        if confounded:
            c1, m1 = make(["1", "2", "3", "4"], ["hfd"] * 4, 3.0)
            c2, m2 = make(["5", "6", "7", "8"], ["lfd"] * 4, 1.0)
        else:
            c1, m1 = make(["1", "2", "3", "4"], ["hfd", "hfd", "lfd", "lfd"], 3.0)
            c2, m2 = make(["5", "6", "7", "8"], ["hfd", "hfd", "lfd", "lfd"], 1.0)

        _write_run(root, "proto", "b1", c1, m1)
        _write_run(root, "proto", "b2", c2, m2)

    def test_balanced_design_recovers_condition_effect_and_absorbs_batch(self) -> None:
        try:
            import pydeseq2  # noqa: F401
        except ImportError:
            self.skipTest("pydeseq2 not installed")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._build(root)

            summary = joint_model.run_joint_model_for_protocol(
                output_root=root,
                protocol_id="proto",
                batch_ids=["b1", "b2"],
                comparisons=[("hfd", "lfd", "hfd_vs_lfd")],
                joint_cfg={},
            )

            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["n_joint_eligible"], 1)

            out = root / "integration" / "joint_model" / "proto"
            self.assertTrue((out / "merged_gene_counts.csv").exists())
            self.assertTrue((out / "merged_sample_metadata.csv").exists())
            self.assertTrue((out / "joint_design_report.csv").exists())

            de = pd.read_csv(out / "liver" / "de_hfd_vs_lfd_joint.csv", index_col=0)
            self.assertEqual(de["source_mode"].iloc[0],
                             "joint_count_model_batch_adjusted")
            self.assertTrue(bool(de["batch_adjusted"].iloc[0]))
            self.assertEqual(int(de["n_batches"].iloc[0]), 2)

            condition_genes = [f"Gene{i}" for i in range(40, 80) if f"Gene{i}" in de.index]
            batch_genes = [f"Gene{i}" for i in range(0, 40) if f"Gene{i}" in de.index]
            null_genes = [f"Gene{i}" for i in range(80, 300) if f"Gene{i}" in de.index]

            # 条件效应基因：应被检出且方向为上调
            cond_sig = de.loc[condition_genes, "significant"].mean()
            self.assertGreater(cond_sig, 0.8,
                               f"condition effect not recovered ({cond_sig:.2f} detected)")
            self.assertGreater(de.loc[condition_genes, "log2FoldChange"].mean(), 1.0)

            # 纯批次效应基因：批次项应把它们吸收掉。
            # 这里用 padj 而非 significant —— significant 叠加了 |log2FC|>1 的
            # 筛选阈值，会掩盖「p 值已经显著、只是效应量没过线」的泄漏。
            batch_leak = float((de.loc[batch_genes, "padj"] < 0.05).mean())
            self.assertLess(batch_leak, 0.15,
                            f"batch effect leaked into DE p-values ({batch_leak:.2%})")

            # 无任何效应的基因：假阳性率应受 FDR 控制
            null_fp = float((de.loc[null_genes, "padj"] < 0.05).mean())
            self.assertLess(null_fp, 0.10,
                            f"FDR control violated on null genes ({null_fp:.2%})")

    def test_confounded_design_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._build(root, confounded=True)

            summary = joint_model.run_joint_model_for_protocol(
                output_root=root,
                protocol_id="proto",
                batch_ids=["b1", "b2"],
                comparisons=[("hfd", "lfd", "hfd_vs_lfd")],
                joint_cfg={},
            )

            self.assertEqual(summary["n_joint_eligible"], 0)
            self.assertEqual(summary["n_meta_only"], 1)
            self.assertEqual(summary["status"], "no_results")

            report = pd.read_csv(
                root / "integration" / "joint_model" / "proto" / "joint_design_report.csv"
            )
            self.assertEqual(report.iloc[0]["route"], joint_model.ROUTE_META)

    def test_provenance_mismatch_blocks_merge(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._build(root)
            # 篡改 b2 的参考指纹
            prov_path = root / "runs" / "proto" / "b2" / "run_provenance.json"
            prov = json.loads(prov_path.read_text(encoding="utf-8"))
            prov["transcriptome_fingerprint"] = "ref-DIFFERENT"
            prov_path.write_text(json.dumps(prov), encoding="utf-8")

            summary = joint_model.run_joint_model_for_protocol(
                output_root=root,
                protocol_id="proto",
                batch_ids=["b1", "b2"],
                comparisons=[("hfd", "lfd", "hfd_vs_lfd")],
                joint_cfg={},
                provenance_policy="fail",
            )

            self.assertEqual(summary["status"], "blocked_provenance_mismatch")
            self.assertFalse(
                (root / "integration" / "joint_model" / "proto"
                 / "merged_gene_counts.csv").exists()
            )

    def test_provenance_warn_policy_allows_merge(self) -> None:
        try:
            import pydeseq2  # noqa: F401
        except ImportError:
            self.skipTest("pydeseq2 not installed")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._build(root)
            prov_path = root / "runs" / "proto" / "b2" / "run_provenance.json"
            prov = json.loads(prov_path.read_text(encoding="utf-8"))
            prov["kallisto_version"] = "kallisto 0.50.1"
            prov_path.write_text(json.dumps(prov), encoding="utf-8")

            summary = joint_model.run_joint_model_for_protocol(
                output_root=root,
                protocol_id="proto",
                batch_ids=["b1", "b2"],
                comparisons=[("hfd", "lfd", "hfd_vs_lfd")],
                joint_cfg={},
                provenance_policy="warn",
            )

            self.assertEqual(summary["status"], "completed")


if __name__ == "__main__":
    unittest.main(verbosity=2)
