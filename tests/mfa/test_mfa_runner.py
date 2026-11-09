#!/usr/bin/env python3
"""MFA 统一运行器：与旧流程逐值一致，比较按 subset 取样本，同一进程中并行运行互不影响。

核心回归：同一个旧版 MFA 实验，旧流程（load_mfa_config → run_from_mfa_config）与
“omicsweave migrate → omicsweave run --omics mfa”得到的差异结果逐值一致
（v2 设计文档阶段 2 验收标准）。
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.adapters.mfa.legacy import convert_legacy_mfa, legacy_all_fits
from omicsweave.contract.dataset import read_dataset
from omicsweave.contract.manifest import read_manifest
from omicsweave.contract.mfa import write_mfa
from omicsweave.contract.tables import read_differential
from omicsweave.core.migrate import migrate_experiment, write_migration
from omicsweave.core.project import project_context
from omicsweave.core.run import run_experiment
from omicsweave.mfa.export import differential_flux_to_contract
from omicsweave.mfa.mfa_analysis import MFAAnalyser
from omicsweave.mfa.mfa_experiment_config import load_mfa_config
from omicsweave.mfa.run_mfa_analysis import run_from_mfa_config
from tests.visualization.test_volcano import assert_volcanoes

EXP = "mfa_exp"
UNITS = "nmol/min/g"
COMPARED = ["tissue", "feature_id", "log2fc", "effect_size", "effect_size_type", "statistic", "pvalue", "padj",
            "padj_method", "n_numerator", "n_denominator", "min_achievable_p", "method", "mean_numerator",
            "mean_denominator", "evidence_level", "selected", "flags", "identifiability"]


def write_legacy_mfa(directory: Path, conditions: dict[str, str], seed: int = 0) -> tuple[Path, Path]:
    """旧版自研软件格式：每个样本 6 次拟合；处理组的 PDH、MDH 逆向与 CS 有差异。"""
    rng = np.random.default_rng(seed)
    directory.mkdir(parents=True, exist_ok=True)
    flux_path, mid_path = directory / "flux_raw_data.xlsx", directory / "mid_raw_data.xlsx"
    names = list(conditions)
    treated = {s: float("treat" in c or "glp1" in c) for s, c in conditions.items()}
    flux_sheets, mid_sheets = {}, {}
    for i, sample in enumerate(names):
        n, t = 6, treated[sample]
        flux_sheets[f"f{i}"] = pd.DataFrame({
            "Index": range(n),
            "loss": rng.uniform(1, 5, n),
            "liver--PDH_m": 10 + 3 * t + rng.normal(0, 0.5, n),
            "liver--MDH_m": 6 + rng.normal(0, 0.3, n),
            "liver--MDH_m__R": 2 + 1.0 * t + rng.normal(0, 0.3, n),
            "heart--PDH_m": 8 + rng.normal(0, 0.5, n),
            "heart--CS_m": 4 - 0.8 * t + rng.normal(0, 0.2, n),
            "MIX:liver--citrate_0": rng.uniform(0, 1, n),
        })
        measured = [0.6, 0.3, 0.1, 0.7, 0.3]
        sheet = pd.DataFrame({"metabolite_name": ["liver--citrate"] * 3 + ["liver--lactate"] * 2 + ["loss"],
                              "mid_index": [0, 1, 2, 0, 1, np.nan], "experiments": measured + [np.nan]})
        for fit in range(3):
            sheet[fit] = list(np.clip(np.array(measured) + rng.normal(0, 0.02, 5), 0, 1)) + [rng.uniform(0.1, 2)]
        mid_sheets[f"m{i}"] = sheet
    for path, sheets in ((flux_path, flux_sheets), (mid_path, mid_sheets)):
        info = pd.DataFrame({"complete_result_name": names, "sheet_name": list(sheets),
                             "condition": list(conditions.values())})
        with pd.ExcelWriter(path) as writer:
            info.to_excel(writer, sheet_name="sheet_information", index=False)
            for name, frame in sheets.items():
                frame.to_excel(writer, sheet_name=name, index=False)
    return flux_path, mid_path


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


def _sorted(table: pd.DataFrame) -> pd.DataFrame:
    return table[COMPARED].sort_values(["tissue", "feature_id"]).reset_index(drop=True)


def _set_mfa_input(root: Path, **fields) -> None:
    path = root / "experiments" / EXP / "experiment.yaml"
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    doc["omics"]["mfa"]["input"].update(fields)
    path.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False), encoding="utf-8")


class MigratedRunMatchesLegacyTest(unittest.TestCase):
    """旧版实验 → 旧流程 vs migrate + omicsweave run。"""

    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = root = Path(cls._tmp.name).resolve()
        conditions = {f"c{i}": "control" for i in range(1, 4)} | {f"t{i}": "treat" for i in range(1, 4)}
        mfa_dir = root / "experiments" / EXP / "mfa_data"
        cls.flux_path, cls.mid_path = write_legacy_mfa(mfa_dir, conditions)
        (mfa_dir / "config.yaml").write_text(yaml.safe_dump({
            "experiment_id": EXP,
            "paths": {"output_dir": "analysis_results/${experiment_id}/mfa_results",
                      "flux_file": f"experiments/{EXP}/mfa_data/flux_raw_data.xlsx",
                      "mid_file": f"experiments/{EXP}/mfa_data/mid_raw_data.xlsx"},
            "comparisons": [["treat", "control"]],
        }), encoding="utf-8")
        (mfa_dir / "groups.yaml").write_text(yaml.safe_dump({
            "conditions": {"control": {"color": "#4C72B0"}, "treat": {"color": "#C44E52"}},
            "sample_groups": conditions,
        }), encoding="utf-8")

        # 旧流程
        with project_context(root):
            cfg = load_mfa_config(experiment_id=EXP)
            analyser = _quiet(run_from_mfa_config, cfg, mode="differential_only")
            legacy = _quiet(analyser.differential_flux, comparisons=[("treat", "control")])
        cls.analyser = analyser
        cls.legacy = differential_flux_to_contract(
            legacy["TREAT_vs_CONTROL"], numerator="treat", denominator="control",
            contrast_id="treat_vs_control", contrast_role="other", flux_units=UNITS)

        # 迁移后按统一模型运行
        write_migration(_quiet(migrate_experiment, EXP, root), root)
        _set_mfa_input(root, model_id="toy_model", flux_units=UNITS)
        cls.run_dir, cls.notes = _quiet(run_experiment, EXP, ["mfa"], root=root, command=["omicsweave", "run"])
        cls.table = read_differential(cls.run_dir / "omics" / "mfa" / "tables" / "differential.parquet")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_differential_matches_legacy_pipeline(self) -> None:
        self.assertEqual(set(self.table["contrast_id"]), {"treat_vs_control"})
        self.assertEqual(len(self.table), 4)
        pd.testing.assert_frame_equal(_sorted(self.table), _sorted(self.legacy), check_dtype=False,
                                      check_exact=False, rtol=1e-12)
        # 有真实效应的通量确实被检出（结果有意义，而不只是两边一样）
        pdh = self.table.set_index(["tissue", "feature_id"]).loc[("liver", "PDH_m")]
        self.assertGreater(pdh["log2fc"], 0)
        self.assertEqual((pdh["n_numerator"], pdh["n_denominator"]), (3, 3))

    def test_all_fits_reconstruction_matches_legacy_loader(self) -> None:
        results = convert_legacy_mfa(self.flux_path, model_id="m", flux_units=UNITS).results
        rebuilt = legacy_all_fits(results)
        self.assertEqual(set(rebuilt), set(self.analyser.flux_all_fits))
        for sample, expected in self.analyser.flux_all_fits.items():
            columns = [c for c in expected.columns if "--" in str(c) and not str(c).startswith("MIX:")] + ["loss"]
            pd.testing.assert_frame_equal(rebuilt[sample][columns], expected[columns].reset_index(drop=True),
                                          check_dtype=False, check_names=False, check_column_type=False)

    def test_bundle_contents(self) -> None:
        mfa = self.run_dir / "omics" / "mfa"
        liver = read_dataset(mfa / "data" / "liver")
        self.assertEqual(sorted(liver.matrix.columns), ["c1", "c2", "c3", "t1", "t2", "t3"])
        self.assertEqual(liver.meta.units, UNITS)
        self.assertTrue((mfa / "qc" / "fitting_loss.parquet").is_file())
        self.assertTrue((mfa / "qc" / "mid_fit_residuals.parquet").is_file())
        self.assertTrue((self.run_dir / "metadata" / "samples.parquet").is_file())

        manifest = read_manifest(self.run_dir)
        self.assertEqual([o.value for o in manifest.omics], ["mfa"])
        inputs = {f.path for f in manifest.files if f.role == "input"}
        self.assertIn(f"experiments/{EXP}/mfa_data/flux_raw_data.xlsx", inputs)
        self.assertIn(f"experiments/{EXP}/experiment.yaml", inputs)
        outputs = {f.path for f in manifest.files if f.role == "output"}
        self.assertIn("omics/mfa/tables/differential.parquet", outputs)
        self.assertEqual(manifest.parameters["mfa"]["contrasts"], ["treat_vs_control"])
        # 火山图随默认步骤画出，并登记在运行清单中
        assert_volcanoes(self, mfa)
        self.assertIn("omics/mfa/figures/volcano_liver_treat_vs_control.pdf", outputs)

    def test_standard_format_input_gives_same_result(self) -> None:
        results = convert_legacy_mfa(self.flux_path, self.mid_path, model_id="toy_model", flux_units=UNITS).results
        standard = write_mfa(results, self.root / "standard_tables")
        _set_mfa_input(self.root, format="mfa_tables", files=[standard.relative_to(self.root).as_posix()])
        try:
            run_dir, _ = _quiet(run_experiment, EXP, ["mfa"], root=self.root)
        finally:
            _set_mfa_input(self.root, format="legacy_mfa_xlsx",
                           files=[f"experiments/{EXP}/mfa_data/flux_raw_data.xlsx",
                                  f"experiments/{EXP}/mfa_data/mid_raw_data.xlsx"])
        table = read_differential(run_dir / "omics" / "mfa" / "tables" / "differential.parquet")
        pd.testing.assert_frame_equal(_sorted(table), _sorted(self.table), check_dtype=False)


def build_two_factor_project(root: Path, seed: int = 0) -> None:
    """2 × 2 设计：diet（LFD/HFD）× treatment（Saline/GLP1），每格 3 只；比较按 subset 取样本。"""
    cells = [(d, t) for d in ("LFD", "HFD") for t in ("Saline", "GLP1")]
    subjects = [(f"m{i * 3 + k}", d, t) for i, (d, t) in enumerate(cells) for k in range(3)]
    exp_dir = root / "experiments" / EXP
    write_legacy_mfa(exp_dir / "mfa", {s: f"{d}_{t}".lower() for s, d, t in subjects}, seed=seed)
    (root / "project.yaml").write_text(yaml.safe_dump({
        "schema_version": 1, "project": {"id": "toy", "species": "mus_musculus"}}), encoding="utf-8")
    (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
        "schema_version": 1,
        "experiment": {"id": EXP},
        "design": {
            "factors": {"diet": {"level": "subject", "levels": ["LFD", "HFD"], "reference": "LFD"},
                        "treatment": {"level": "subject", "levels": ["Saline", "GLP1"], "reference": "Saline"}},
            "contrasts": [
                {"id": "HFD_vs_LFD", "factor": "diet", "numerator": "HFD", "denominator": "LFD",
                 "subset": {"treatment": "Saline"}, "role": "disease_effect"},
                {"id": "GLP1_vs_Saline_in_HFD", "factor": "treatment", "numerator": "GLP1", "denominator": "Saline",
                 "subset": {"diet": "HFD"}, "role": "treatment_effect"},
            ],
            "reversal_pairs": [{"disease": "HFD_vs_LFD", "treatment": "GLP1_vs_Saline_in_HFD"}],
        },
        "omics": {"mfa": {"enabled": True, "input": {
            "level": "table", "format": "legacy_mfa_xlsx", "model_id": "toy", "flux_units": UNITS,
            "files": [f"experiments/{EXP}/mfa/flux_raw_data.xlsx", f"experiments/{EXP}/mfa/mid_raw_data.xlsx"]}}},
    }, sort_keys=False), encoding="utf-8")
    pd.DataFrame(subjects, columns=["subject_id", "diet", "treatment"]).to_csv(
        exp_dir / "subjects.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": [s for s, _, _ in subjects], "subject_id": [s for s, _, _ in subjects],
                  "omics": "mfa", "tissue": "multi_tissue"}).to_csv(exp_dir / "samples.tsv", sep="\t", index=False)


class SubsetContrastTest(unittest.TestCase):
    def test_contrasts_use_subset_samples_and_roles(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            build_two_factor_project(root)
            run_dir, _ = _quiet(run_experiment, EXP, ["mfa"], root=root)
            table = read_differential(run_dir / "omics" / "mfa" / "tables" / "differential.parquet")
        roles = dict(zip(table["contrast_id"], table["contrast_role"]))
        self.assertEqual(roles, {"HFD_vs_LFD": "disease_effect", "GLP1_vs_Saline_in_HFD": "treatment_effect"})
        self.assertTrue((table["n_numerator"] == 3).all() and (table["n_denominator"] == 3).all())
        # GLP1 组 PDH 升高（处理效应），饮食比较中两组都是 Saline，PDH 无系统差异
        pdh = table[(table["tissue"] == "liver") & (table["feature_id"] == "PDH_m")].set_index("contrast_id")
        self.assertGreater(pdh.at["GLP1_vs_Saline_in_HFD", "log2fc"], 0.2)
        self.assertLess(abs(pdh.at["HFD_vs_LFD", "log2fc"]), 0.2)


class ParallelRunTest(unittest.TestCase):
    def test_two_experiments_in_parallel_match_sequential_runs(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            roots = [Path(tmp, name).resolve() for name in ("a", "b")]
            for seed, root in enumerate(roots):
                build_two_factor_project(root, seed=seed)

            def run(root: Path) -> pd.DataFrame:
                run_dir, _ = run_experiment(EXP, ["mfa"], root=root)
                return read_differential(run_dir / "omics" / "mfa" / "tables" / "differential.parquet")

            with contextlib.redirect_stdout(io.StringIO()):
                sequential = [run(r) for r in roots]
                with ThreadPoolExecutor(max_workers=2) as pool:
                    parallel = list(pool.map(run, roots))
        for expected, actual in zip(sequential, parallel):
            pd.testing.assert_frame_equal(actual, expected)
        self.assertFalse(sequential[0]["log2fc"].equals(sequential[1]["log2fc"]))


if __name__ == "__main__":
    unittest.main()
