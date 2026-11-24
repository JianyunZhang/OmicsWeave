#!/usr/bin/env python3
"""转录组通用样本质检（sample_qc.py）与运行器过滤/可靠性联动测试。"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.contract.dataset import read_dataset
from omicsweave.contract.tables import read_differential
from omicsweave.core.run import run_experiment
from omicsweave.transcriptomics.sample_qc import (
    SAMPLE_QUALITY_COLUMNS,
    evaluate_sample_quality,
)


def _make_multi_tissue_counts(
    n_per_tissue: int = 8,
    n_bg: int = 180,
    n_markers: int = 15,
    seed: int = 42,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """构造 3 组织（liver / muscle / adipose）干净计数矩阵，每组织各有 15 个高特异标记与 5 个线粒体基因。"""
    rng = np.random.default_rng(seed)
    tissues = ("liver", "muscle", "adipose")
    genes = (
        [f"Bg{i:04d}" for i in range(n_bg)]
        + [f"LivM{i:02d}" for i in range(n_markers)]
        + [f"MusM{i:02d}" for i in range(n_markers)]
        + [f"AdiM{i:02d}" for i in range(n_markers)]
        + [f"mt-Gene{i}" for i in range(5)]
    )
    bg_base = rng.lognormal(8.5, 0.9, n_bg)
    data: dict[str, np.ndarray] = {}
    tissue_map: dict[str, str] = {}
    for tissue in tissues:
        for i in range(n_per_tissue):
            sid = f"rna_{tissue}_{i}"
            tissue_map[sid] = tissue
            bg = bg_base * rng.lognormal(0.0, 0.12, n_bg)
            liv = np.full(n_markers, 12000.0 if tissue == "liver" else 0.5) * rng.lognormal(0.0, 0.10, n_markers)
            mus = np.full(n_markers, 12000.0 if tissue == "muscle" else 0.5) * rng.lognormal(0.0, 0.10, n_markers)
            adi = np.full(n_markers, 12000.0 if tissue == "adipose" else 0.5) * rng.lognormal(0.0, 0.10, n_markers)
            mt = np.full(5, 30000.0) * rng.lognormal(0.0, 0.10, 5)
            means = np.concatenate([bg, liv, mus, adi, mt])
            data[sid] = rng.poisson(np.clip(means, 0.0, 1e7)).astype(float)
    return pd.DataFrame(data, index=pd.Index(genes, name="gene")), tissue_map


class SampleQCContractAndCleanTest(unittest.TestCase):
    def test_schema_matches_core_reliability(self) -> None:
        self.assertEqual(SAMPLE_QUALITY_COLUMNS, ["sample_id", "metric", "value", "robust_z", "flag", "reason"])

    def test_clean_multi_tissue_matrix_has_zero_flags(self) -> None:
        counts, tissue_map = _make_multi_tissue_counts()
        res = evaluate_sample_quality(counts, tissue_map)
        self.assertEqual(res.excluded, set())
        self.assertEqual(res.flagged, set())
        self.assertTrue(res.quality["flag"].eq("").all())
        self.assertEqual(len(res.metrics), len(tissue_map))
        self.assertTrue(res.metrics["qc_status"].eq("ok").all())


class SampleQCAnomalyDetectionTest(unittest.TestCase):
    def test_cross_tissue_contamination_exclude_and_flag(self) -> None:
        counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=8, seed=1)
        # 2.5% 肌肉 RNA 混入脂肪 -> exclude
        counts["rna_adipose_0"] = counts["rna_adipose_0"] * 0.975 + counts["rna_muscle_0"] * 0.025
        # 0.25% 肌肉 RNA 混入脂肪 -> flag（不 exclude）
        counts["rna_adipose_1"] = counts["rna_adipose_1"] * 0.9975 + counts["rna_muscle_1"] * 0.0025

        res = evaluate_sample_quality(counts, tissue_map)
        self.assertIn("rna_adipose_0", res.excluded)
        self.assertIn("rna_adipose_1", res.flagged)
        self.assertNotIn("rna_adipose_1", res.excluded)
        # 其他干净样本不得误报
        clean_other = {s for s in tissue_map if s not in ("rna_adipose_0", "rna_adipose_1")}
        self.assertEqual((res.excluded | res.flagged) & clean_other, set())

        q0 = res.quality[(res.quality["sample_id"] == "rna_adipose_0")
                         & (res.quality["metric"] == "cross_tissue_contamination")].iloc[0]
        self.assertEqual(q0["flag"], "exclude")
        self.assertEqual(q0["foreign_tissue"], "muscle")

    def test_own_marker_depletion_catches_group_failure(self) -> None:
        counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=8, seed=2)
        # 3 只脂肪样本（模拟解剖混入淋巴结）本组织特异标记降至 15%
        adi_markers = [g for g in counts.index if g.startswith("AdiM")]
        for sid in ("rna_adipose_0", "rna_adipose_1", "rna_adipose_2"):
            counts.loc[adi_markers, sid] = counts.loc[adi_markers, sid] * 0.15

        res = evaluate_sample_quality(counts, tissue_map)
        self.assertTrue({"rna_adipose_0", "rna_adipose_1", "rna_adipose_2"} <= res.excluded)
        # 其余 5 只正常脂肪样本不受影响
        self.assertEqual(res.excluded & {f"rna_adipose_{i}" for i in range(3, 8)}, set())

    def test_single_tissue_correlation_outlier(self) -> None:
        counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=8, seed=3)
        liver_cols = [s for s, t in tissue_map.items() if t == "liver"]
        sub = counts[liver_cols].copy()
        # 严重打乱 rna_liver_0 的全转录组表达谱 -> exclude
        rng = np.random.default_rng(99)
        sub["rna_liver_0"] = rng.permutation(sub["rna_liver_0"].to_numpy())

        res = evaluate_sample_quality(sub, {s: "liver" for s in liver_cols})
        self.assertIn("rna_liver_0", res.excluded)
        self.assertEqual(res.excluded - {"rna_liver_0"}, set())

    def test_mitochondrial_outlier_without_penalizing_high_mt_tissue(self) -> None:
        counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=8, seed=4)
        mt_genes = [g for g in counts.index if g.startswith("mt-")]
        # 模拟正常高线粒体组织（如心脏：全体样本约 75% 线粒体读数）
        for s, t in tissue_map.items():
            if t == "muscle":
                non_mt_sum = counts.loc[~counts.index.isin(mt_genes), s].sum()
                counts.loc[mt_genes, s] = (non_mt_sum * 3.0) / len(mt_genes)
        # liver 中仅 1 只样本线粒体飙升至 72%（相对同组织中位数 ~5% 偏移巨大）-> exclude
        non_mt_l0 = counts.loc[~counts.index.isin(mt_genes), "rna_liver_0"].sum()
        counts.loc[mt_genes, "rna_liver_0"] = (non_mt_l0 * 2.6) / len(mt_genes)

        res = evaluate_sample_quality(counts, tissue_map)
        self.assertIn("rna_liver_0", res.excluded)
        # 全体健康高 MT 肌肉样本不应被误杀
        muscle_ids = {s for s, t in tissue_map.items() if t == "muscle"}
        self.assertEqual(res.excluded & muscle_ids, set())
        self.assertEqual(res.flagged & muscle_ids, set())

    def test_low_mitochondrial_fraction_not_flagged_by_default(self) -> None:
        counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=8, n_bg=300, seed=10)
        mt_genes = [g for g in counts.index if g.startswith("mt-")]
        for s, t in tissue_map.items():
            if t == "muscle":
                non_mt_sum = counts.loc[~counts.index.isin(mt_genes), s].sum()
                counts.loc[mt_genes, s] = (non_mt_sum * 0.56) / len(mt_genes)  # ~36% MT
        # 单只样本线粒体占比偏低（~6.5% vs 同组织中位数 ~36%）：默认不罚为异常
        non_mt_m0 = counts.loc[~counts.index.isin(mt_genes), "rna_muscle_0"].sum()
        counts.loc[mt_genes, "rna_muscle_0"] = (non_mt_m0 * 0.07) / len(mt_genes)
        res_default = evaluate_sample_quality(counts, tissue_map)
        self.assertNotIn("rna_muscle_0", res_default.flagged | res_default.excluded)
        q_def = res_default.quality[
            (res_default.quality["sample_id"] == "rna_muscle_0")
            & (res_default.quality["metric"] == "mitochondrial_fraction")
        ].iloc[0]
        self.assertEqual(q_def["flag"], "")
        res_strict = evaluate_sample_quality(counts, tissue_map, config={"mt_flag_low": True})
        self.assertIn("rna_muscle_0", res_strict.flagged)
        q_str = res_strict.quality[
            (res_strict.quality["sample_id"] == "rna_muscle_0")
            & (res_strict.quality["metric"] == "mitochondrial_fraction")
        ].iloc[0]
        self.assertEqual(q_str["flag"], "flag")

    def test_group_consistent_remodeling_not_flagged_as_var_outlier(self) -> None:
        counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=8, seed=11)
        liver_cols = [s for s, t in tissue_map.items() if t == "liver"]
        # 两组（4 ctrl vs 4 treated）在部分背景基因上存在强烈一致性生物学重塑，使全群体 r_var_pw 降低
        sub = counts[liver_cols].copy()
        rng = np.random.default_rng(12)
        shift_genes = [f"Bg{i:04d}" for i in range(90)]
        for s in liver_cols[4:]:
            sub.loc[shift_genes[:45], s] = sub.loc[shift_genes[:45], s] * 12.0
            sub.loc[shift_genes[45:], s] = sub.loc[shift_genes[45:], s] * 0.08
        # 加入组内极小噪声，使组内相关极高（r_within >= 0.95）而跨组 r_var_pw 降低
        groups = {s: ("ctrl" if i < 4 else "treated") for i, s in enumerate(liver_cols)}
        sub.loc["Bg0000", liver_cols[4]] *= rng.uniform(0.98, 1.02)
        res_grouped = evaluate_sample_quality(
            sub, {s: "liver" for s in liver_cols},
            config={"var_corr_flag_min": 0.99, "var_corr_flag_robust_z": 0.5},
            groups=groups,
        )
        self.assertEqual(res_grouped.flagged & set(liver_cols[4:]), set())


class TranscriptomicsRunnerSampleQCTest(unittest.TestCase):
    def _build_project(self, root: Path, *, filter_excluded: bool = True) -> Path:
        counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=6, seed=7)
        # 构造真实处理效应：前 15 个背景基因在 HFD 中上调 4 倍
        subjects = []
        sample_rows = []
        for i in range(6):
            diet = "LFD" if i < 3 else "HFD"
            subjects.append({"subject_id": f"m{i}", "diet": diet})
            for tissue in ("liver", "muscle", "adipose"):
                sid = f"rna_{tissue}_{i}"
                sample_rows.append({"sample_id": sid, "subject_id": f"m{i}",
                                    "omics": "transcriptomics", "tissue": tissue})
                if diet == "HFD":
                    counts.loc[[f"Bg{k:04d}" for k in range(15)], sid] *= 4.0

        # 注入 2 类异常：
        # 1) adipose HFD 中 2/3 样本混入 4% 肌肉 RNA -> exclude（且 HFD 组 2/3 被标 -> qc_confounded / qc_insufficient）
        counts["rna_adipose_3"] = counts["rna_adipose_3"] * 0.96 + counts["rna_muscle_3"] * 0.04
        counts["rna_adipose_4"] = counts["rna_adipose_4"] * 0.96 + counts["rna_muscle_4"] * 0.04
        # 2) liver HFD 中仅 1/3 样本混入 4% 肌肉 RNA -> exclude，但剩余 2 vs 3 仍可分析（reliability = downweighted）
        counts["rna_liver_3"] = counts["rna_liver_3"] * 0.96 + counts["rna_muscle_3"] * 0.04

        exp = "qc_exp"
        exp_dir = root / "experiments" / exp
        exp_dir.mkdir(parents=True)
        counts_path = exp_dir / "counts.csv"
        counts.to_csv(counts_path)
        pd.DataFrame(subjects).to_csv(exp_dir / "subjects.tsv", sep="\t", index=False)
        pd.DataFrame(sample_rows).to_csv(exp_dir / "samples.tsv", sep="\t", index=False)
        (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                           encoding="utf-8")
        (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
            "schema_version": 1,
            "experiment": {"id": exp},
            "design": {
                "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"], "reference": "LFD"}},
                "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD", "denominator": "LFD"}],
            },
            "omics": {"transcriptomics": {
                "enabled": True,
                "input": {"level": "table", "format": "counts_matrix", "files": [f"experiments/{exp}/counts.csv"]},
                "analysis": {"sample_qc": {"filter_excluded": filter_excluded}},
            }},
        }, sort_keys=False), encoding="utf-8")
        return exp_dir

    def test_runner_filters_excluded_and_feeds_contrast_reliability(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            self._build_project(root, filter_excluded=True)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                run_dir, notes = run_experiment("qc_exp", ["transcriptomics"], steps=["qc", "differential"], root=root)

            tx = run_dir / "omics" / "transcriptomics"
            sq = pd.read_parquet(tx / "qc" / "sample_quality.parquet")
            self.assertEqual(set(sq.loc[sq["flag"] == "exclude", "sample_id"]),
                             {"rna_adipose_3", "rna_adipose_4", "rna_liver_3"})
            # library_sizes 保留过滤前全部 18 个原始样本
            lib = pd.read_parquet(tx / "qc" / "library_sizes.parquet")
            self.assertEqual(len(lib), 18)

            # data/<tissue> 已剔除 exclude 样本
            liver_ds = read_dataset(tx / "data" / "liver")
            self.assertNotIn("rna_liver_3", liver_ds.matrix.columns)
            self.assertEqual(len(liver_ds.matrix.columns), 5)

            adi_ds = read_dataset(tx / "data" / "adipose")
            self.assertEqual(len(adi_ds.matrix.columns), 4)

            # contrast_reliability 自动读取 sample_quality.parquet
            rel = pd.read_parquet(run_dir / "qc" / "contrast_reliability.parquet").set_index("tissue")
            self.assertEqual(rel.at["muscle", "status"], "ok")
            self.assertEqual(rel.at["liver", "status"], "downweighted")
            self.assertIn(rel.at["adipose", "status"], {"qc_confounded", "qc_insufficient"})

            # 差异表保留干净样本量，且 adipose 自动带阻塞标记
            diff = read_differential(tx / "tables" / "differential.parquet")
            liver_diff = diff[diff["tissue"] == "liver"]
            self.assertTrue((liver_diff["n_numerator"] == 2).all() and (liver_diff["n_denominator"] == 3).all())
            self.assertTrue(any("通用质检" in n for n in notes))

    def test_runner_respects_filter_excluded_false(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            self._build_project(root, filter_excluded=False)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                run_dir, _ = run_experiment("qc_exp", ["transcriptomics"], steps=["qc", "differential"], root=root)

            tx = run_dir / "omics" / "transcriptomics"
            liver_ds = read_dataset(tx / "data" / "liver")
            # 关闭过滤时样本留在矩阵中，但 sample_quality 与 contrast_reliability 仍完整记录
            self.assertIn("rna_liver_3", liver_ds.matrix.columns)
            sq = pd.read_parquet(tx / "qc" / "sample_quality.parquet")
            self.assertIn("rna_liver_3", set(sq.loc[sq["flag"] == "exclude", "sample_id"]))

    def test_trace_contamination_downweights_contrast_and_flags_foreign_markers(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            counts, tissue_map = _make_multi_tissue_counts(n_per_tissue=6, seed=13)
            subjects, sample_rows = [], []
            for i in range(6):
                diet = "LFD" if i < 3 else "HFD"
                subjects.append({"subject_id": f"m{i}", "diet": diet})
                for tissue in ("liver", "muscle", "adipose"):
                    sid = f"rna_{tissue}_{i}"
                    sample_rows.append({"sample_id": sid, "subject_id": f"m{i}",
                                        "omics": "transcriptomics", "tissue": tissue})
                    if diet == "HFD":
                        counts.loc[[f"Bg{k:04d}" for k in range(15)], sid] *= 4.0
            # 在 liver HFD 的 2 只样本中注入 0.25% 肌肉 RNA -> 仅标 flag（不 exclude）
            counts["rna_liver_3"] = counts["rna_liver_3"] * 0.9975 + counts["rna_muscle_3"] * 0.0025
            counts["rna_liver_4"] = counts["rna_liver_4"] * 0.9975 + counts["rna_muscle_4"] * 0.0025

            exp = "trace_exp"
            exp_dir = root / "experiments" / exp
            exp_dir.mkdir(parents=True)
            counts.to_csv(exp_dir / "counts.csv")
            pd.DataFrame(subjects).to_csv(exp_dir / "subjects.tsv", sep="\t", index=False)
            pd.DataFrame(sample_rows).to_csv(exp_dir / "samples.tsv", sep="\t", index=False)
            (root / "project.yaml").write_text("schema_version: 1\nproject: {id: p, species: mus_musculus}\n",
                                               encoding="utf-8")
            (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
                "schema_version": 1,
                "experiment": {"id": exp},
                "design": {
                    "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"], "reference": "LFD"}},
                    "contrasts": [{"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD", "denominator": "LFD"}],
                },
                "omics": {"transcriptomics": {
                    "enabled": True,
                    "input": {"level": "table", "format": "counts_matrix", "files": [f"experiments/{exp}/counts.csv"]},
                    "analysis": {"sample_qc": {"filter_excluded": True}},
                }},
            }, sort_keys=False), encoding="utf-8")

            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                run_dir, _ = run_experiment(exp, ["transcriptomics"], steps=["qc", "differential"], root=root)

            rel = pd.read_parquet(run_dir / "qc" / "contrast_reliability.parquet").set_index("tissue")
            self.assertEqual(rel.at["liver", "status"], "downweighted")

            diff = read_differential(run_dir / "omics" / "transcriptomics" / "tables" / "differential.parquet")
            liver_diff = diff[diff["tissue"] == "liver"].set_index("feature_id")
            self.assertIn("foreign_tissue_marker", str(liver_diff.at["MusM00", "flags"]).split(";"))
            self.assertNotIn("foreign_tissue_marker", str(liver_diff.at["Bg0000", "flags"]).split(";"))


if __name__ == "__main__":
    unittest.main()
