#!/usr/bin/env python3
"""多示踪剂整合：分层置换的正确性、可辨识性门控、示踪剂冲突，以及经 omicsweave run 的端到端结果。

端到端数据：同一模型、四种示踪剂（葡萄糖、乳酸、甘油、丙氨酸），每种示踪剂一个独立的 3 vs 3 动物队列：

- liver PDH_m：四种示踪剂都约束得好、效应一致 → multi_tracer_consistent，合并后可达 p < 0.05；
- liver GK_c：只有甘油约束它（其余示踪剂拟合点极分散）→ single_tracer，只给效应量；
- liver LDH_c：乳酸给出与其余示踪剂相反的效应（如输注乳酸本身扰动了通量）→ multi_tracer_conflict，不入选；
- heart CS_m：没有效应。
"""

from __future__ import annotations

import contextlib
import io
import math
import tempfile
import unittest
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from omicsweave.adapters.mfa.legacy import read_input_sets
from omicsweave.contract.dataset import ContractError
from omicsweave.contract.manifest import read_manifest
from omicsweave.contract.mfa import MFAInfo, MFAResults, merge_mfa, write_mfa
from omicsweave.contract.tables import read_differential
from omicsweave.core.run import RunError, run_experiment
from omicsweave.mfa import multi_tracer as mt
from tests.mfa.test_mfa_runner import write_legacy_mfa

EXP = "tracer_exp"
UNITS = "nmol/min/g"
TRACERS = ("alanine", "glucose", "glycerol", "lactate")
N_FITS = 20
# (组织, 反应, 对照组均值, 各示踪剂的处理效应, 各示踪剂的拟合点标准差)
WELL, POOR = 0.2, 300.0
REACTIONS = [
    ("liver", "PDH_m", 10.0, dict.fromkeys(TRACERS, 6.0), dict.fromkeys(TRACERS, WELL)),
    ("liver", "GK_c", 5.0, dict.fromkeys(TRACERS, 2.5), {t: WELL if t == "glycerol" else POOR for t in TRACERS}),
    ("liver", "LDH_c", 20.0, {t: 8.0 if t == "lactate" else -4.0 for t in TRACERS}, dict.fromkeys(TRACERS, WELL)),
    ("heart", "CS_m", 4.0, dict.fromkeys(TRACERS, 0.0), dict.fromkeys(TRACERS, WELL)),
]


def _quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


def tracer_results(tracer: str, samples: dict[str, bool], rng: np.random.Generator,
                   model_id: str = "central_carbon_v1", label: str | None = None) -> MFAResults:
    """一种示踪剂的标准格式结果：每只动物 N_FITS 次拟合。``samples`` 为样本 → 是否处理组，
    ``label`` 为写入 ``info.tracer`` 的示踪剂名（默认 ``13C <tracer>``）。"""
    flux, fits = [], []
    for sample, treated in samples.items():
        for tissue, reaction, base, effect, fit_sd in REACTIONS:
            animal = base + effect[tracer] * treated + rng.normal(0, 0.5)
            values = animal + rng.normal(0, fit_sd[tracer], N_FITS)
            flux += [(sample, tissue, reaction, "forward", str(i), v) for i, v in enumerate(values)]
        fits += [(sample, str(i), rng.uniform(1, 5)) for i in range(N_FITS)]
    return MFAResults(
        info=MFAInfo(model_id=model_id, flux_units=UNITS, tracer=label or f"13C {tracer}"),
        flux=pd.DataFrame(flux, columns=["sample_id", "tissue", "reaction_id", "direction", "fit_id", "value"]),
        fits=pd.DataFrame(fits, columns=["sample_id", "fit_id", "loss"]),
    )


def build_project(root: Path, *, cohorts: list[tuple[str, str, str, str]] | None = None, protocol: bool = True,
                  seed: int = 1) -> list[dict[str, str]]:
    """每组结果一个 3 vs 3 动物队列、一个 Parquet 结果目录。返回 samples.tsv 的行。

    ``cohorts`` 的每项为（目录名, 生成数据所用的示踪剂, ``info.tracer``, samples.tsv 的 ``protocol``）；
    默认四种示踪剂各一组。``protocol=False`` 时 samples.tsv 的 protocol 列留空。
    """
    cohorts = cohorts or [(t, t, f"13C {t}", t) for t in TRACERS]
    rng = np.random.default_rng(seed)
    exp_dir = root / "experiments" / EXP
    rows, files = [], []
    for name, tracer, label, protocol_value in cohorts:
        samples = {f"{name}_{g}{i}": g == "t" for g in ("c", "t") for i in range(1, 4)}
        path = exp_dir / "mfa" / name
        write_mfa(tracer_results(tracer, samples, rng, label=label), path)
        files.append(path.relative_to(root).as_posix())
        rows += [{"sample_id": s, "subject_id": s, "omics": "mfa", "tissue": "multi_tissue",
                  "protocol": protocol_value if protocol else "", "group": "treat" if t else "control"}
                 for s, t in samples.items()]
    (root / "project.yaml").write_text(yaml.safe_dump({
        "schema_version": 1, "project": {"id": "toy", "species": "mus_musculus"}}), encoding="utf-8")
    (exp_dir / "experiment.yaml").write_text(yaml.safe_dump({
        "schema_version": 1,
        "experiment": {"id": EXP},
        "design": {
            "factors": {"group": {"level": "subject", "levels": ["control", "treat"], "reference": "control"}},
            "contrasts": [{"id": "treat_vs_control", "factor": "group", "numerator": "treat",
                           "denominator": "control", "role": "treatment_effect"}],
        },
        "omics": {"mfa": {"enabled": True, "input": {"level": "table", "format": "mfa_tables", "files": files}}},
    }, sort_keys=False), encoding="utf-8")
    table = pd.DataFrame(rows)
    table[["subject_id", "group"]].to_csv(exp_dir / "subjects.tsv", sep="\t", index=False)
    table.drop(columns="group").to_csv(exp_dir / "samples.tsv", sep="\t", index=False)
    return rows


def _stratum(tracer: str, num: list[float], den: list[float], fit_var: float = 0.0) -> mt._Stratum:
    values = np.array(num + den, dtype=float)
    return mt._Stratum(tracer, values, np.full(len(values), fit_var), np.array([True] * len(num) + [False] * len(den)))


def _brute_force_pvalue(strata: list[mt._Stratum]) -> float:
    """逐一枚举各示踪剂内的全部组划分，按定义计算合并 z 的双侧置换 p。"""
    def z_of(masks):
        num = den = 0.0
        for s, m in zip(strata, masks):
            a, b = s.values[m], s.values[~m]
            v = (a.var(ddof=1) + s.fit_var[m].mean()) / len(a) + (b.var(ddof=1) + s.fit_var[~m].mean()) / len(b)
            num += (a.mean() - b.mean()) / v
            den += 1.0 / v
        return num / math.sqrt(den)

    def masks_of(s):
        for chosen in combinations(range(s.n), s.n_num):
            mask = np.zeros(s.n, dtype=bool)
            mask[list(chosen)] = True
            yield mask

    observed = abs(z_of([s.is_num for s in strata]))
    values = [abs(z_of(choice)) for choice in product(*(list(masks_of(s)) for s in strata))]
    return float(np.mean(np.array(values) >= observed * (1 - 1e-9)))


class PermutationTest(unittest.TestCase):
    def test_matches_brute_force_enumeration(self) -> None:
        strata = [_stratum("glucose", [10.2, 11.0, 12.5], [9.0, 9.8, 10.1], fit_var=0.3),
                  _stratum("lactate", [8.0, 12.0, 11.5, 9.9], [9.5, 7.0, 8.8], fit_var=1.0)]
        result = mt._permutation_test(strata, mt.PoolingSettings(), seed=0)
        self.assertTrue(result["exact"])
        self.assertEqual(result["n_permutations"], 20 * 35)
        self.assertAlmostEqual(result["pvalue"], _brute_force_pvalue(strata), places=12)
        # 第二个示踪剂两组不等大，不存在镜像划分
        self.assertAlmostEqual(result["p_floor"], 1 / (20 * 35))

    def test_monte_carlo_agrees_with_exact(self) -> None:
        strata = [_stratum(t, list(rng.normal(0.8, 1, 3)), list(rng.normal(0, 1, 3)))
                  for t, rng in ((t, np.random.default_rng(i)) for i, t in enumerate(TRACERS[:3]))]
        exact = mt._permutation_test(strata, mt.PoolingSettings(), seed=0)
        sampled = mt._permutation_test(strata, mt.PoolingSettings(max_exact_permutations=10, n_permutations=40_000),
                                       seed=3)
        self.assertTrue(exact["exact"])
        self.assertFalse(sampled["exact"])
        self.assertLess(abs(exact["pvalue"] - sampled["pvalue"]), 0.01)

    def test_more_tracers_lower_the_power_floor(self) -> None:
        one = mt._permutation_test([_stratum("a", [5, 6, 7], [1, 2, 3])], mt.PoolingSettings(), seed=0)
        self.assertAlmostEqual(one["p_floor"], 0.1)       # 与 Mann-Whitney 3 vs 3 的下界相同
        self.assertAlmostEqual(one["pvalue"], 0.1)
        four = mt._permutation_test([_stratum(t, [5, 6, 7], [1, 2, 3]) for t in TRACERS], mt.PoolingSettings(), 0)
        self.assertAlmostEqual(four["p_floor"], 2 / 20 ** 4)
        self.assertAlmostEqual(four["pvalue"], 2 / 20 ** 4)

    def test_zero_variance_split_stays_finite(self) -> None:
        result = mt._permutation_test([_stratum("a", [5, 5, 5], [1, 1, 1])], mt.PoolingSettings(), seed=0)
        self.assertTrue(np.isfinite(result["z"]))
        self.assertAlmostEqual(result["pvalue"], 0.1)


def _animals(spec: dict[str, tuple[list[float], list[float]]], fit_std: float = 0.1,
             reaction: str = "R1") -> tuple[pd.DataFrame, dict[str, str]]:
    rows, groups = [], {}
    for tracer, (num, den) in spec.items():
        for group, values in (("treat", num), ("control", den)):
            for i, value in enumerate(values):
                sample = f"{tracer}_{group}{i}"
                groups[sample] = group
                rows.append({"sample_id": sample, "tracer": tracer, "tissue": "liver", "reaction": reaction,
                             "net_flux": value, "fit_std": fit_std})
    return pd.DataFrame(rows), groups


def _coverage(labels: dict[str, str], reaction: str = "R1") -> pd.DataFrame:
    table = pd.DataFrame({"tracer": list(labels), "tissue": "liver", "reaction": reaction,
                          "identifiability": list(labels.values())})
    return mt.tracer_coverage(table, mt.PoolingSettings())


class PoolingTest(unittest.TestCase):
    def test_fixed_effect_estimate_and_heterogeneity(self) -> None:
        spec = {"glucose": ([12.0, 13.0, 12.5], [10.0, 10.4, 9.9]), "lactate": ([14.0, 11.0, 13.5], [10.5, 9.0, 11.0])}
        animals, groups = _animals(spec, fit_std=0.5)
        pooled = mt.pool_contrast(animals, groups, "treat", "control",
                                  _coverage(dict.fromkeys(spec, "well_constrained")), mt.PoolingSettings())
        row = pooled.iloc[0]
        deltas, weights = [], []
        for num, den in spec.values():
            a, b = np.array(num), np.array(den)
            deltas.append(a.mean() - b.mean())
            weights.append(1 / ((a.var(ddof=1) + 0.25) / 3 + (b.var(ddof=1) + 0.25) / 3))
        deltas, weights = np.array(deltas), np.array(weights)
        expected = (weights * deltas).sum() / weights.sum()
        self.assertAlmostEqual(row["diff"], expected)
        self.assertAlmostEqual(row["diff_se"], math.sqrt(1 / weights.sum()))
        self.assertAlmostEqual(row["heterogeneity_q"], (weights * (deltas - expected) ** 2).sum())
        # 组均值用同一组权重，两者之差就是合并的组间差
        self.assertAlmostEqual(row["mean_numerator"] - row["mean_denominator"], row["diff"])
        self.assertEqual(row["tracer_evidence"], mt.CONSISTENT)
        self.assertEqual((row["n_numerator"], row["n_denominator"]), (6, 6))

    def test_unconstrained_tracers_are_left_out(self) -> None:
        spec = {t: ([7.5, 7.6, 7.4], [5.0, 5.1, 4.9]) for t in TRACERS}
        animals, groups = _animals(spec)
        labels = {"glycerol": "well_constrained", "glucose": "poorly_constrained", "lactate": "bound_pinned",
                  "alanine": "poorly_constrained"}
        row = mt.pool_contrast(animals, groups, "treat", "control", _coverage(labels), mt.PoolingSettings()).iloc[0]
        self.assertEqual(row["tracers_informative"], "glycerol")
        self.assertEqual(row["tracer_evidence"], mt.SINGLE)
        self.assertIn("lactate:bound_pinned", row["tracers_excluded"])
        self.assertEqual((row["n_numerator"], row["n_denominator"]), (3, 3))

        table = mt.to_differential(pd.DataFrame([row]), contrast_id="c", contrast_role="other", flux_units=UNITS,
                                   settings=mt.PoolingSettings())
        # 单示踪剂 3 vs 3：置换检验与秩检验同样达不到 p < 0.05，只给效应量
        self.assertEqual(table.at[0, "evidence_level"], "effect_size_only")
        self.assertTrue(np.isnan(table.at[0, "pvalue"]))
        self.assertTrue(table.at[0, "selected"])
        self.assertIn("single_tracer", table.at[0, "flags"])

    def test_no_informative_tracer(self) -> None:
        spec = {t: ([7.5, 7.6, 7.4], [5.0, 5.1, 4.9]) for t in TRACERS[:2]}
        animals, groups = _animals(spec)
        pooled = mt.pool_contrast(animals, groups, "treat", "control",
                                  _coverage(dict.fromkeys(spec, "poorly_constrained")), mt.PoolingSettings())
        table = mt.to_differential(pooled, contrast_id="c", contrast_role="other", flux_units=UNITS,
                                   settings=mt.PoolingSettings())
        self.assertEqual(pooled.at[0, "tracer_evidence"], mt.NOT_IDENTIFIABLE)
        self.assertEqual(table.at[0, "evidence_level"], "insufficient")
        self.assertFalse(table.at[0, "selected"])
        self.assertIn("identifiability:poorly_constrained", table.at[0, "flags"])

    def test_smd_weighting_modes_and_validation(self) -> None:
        # 两个同向示踪剂：glucose 拟合方差小且 well_constrained，lactate 拟合方差大且 moderate
        spec = {
            "glucose": ([12.0, 13.0, 12.5], [10.0, 10.4, 9.9]),
            "lactate": ([13.5, 11.2, 13.0], [10.5, 9.4, 10.8]),
        }
        animals, groups = _animals(spec, fit_std=0.2)
        animals.loc[animals["tracer"] == "lactate", "fit_std"] = 3.5
        cov = _coverage({"glucose": "well_constrained", "lactate": "moderate"})

        def parse_weights(text: str) -> dict[str, float]:
            return {k: float(v) for k, v in (item.split(":") for item in str(text).split(";") if ":" in item)}

        r_smd = mt.pool_contrast(animals, groups, "treat", "control", cov,
                                 mt.PoolingSettings(smd_weighting="smd")).iloc[0]
        r_rel = mt.pool_contrast(animals, groups, "treat", "control", cov,
                                 mt.PoolingSettings(smd_weighting="reliability_smd")).iloc[0]
        r_id = mt.pool_contrast(animals, groups, "treat", "control", cov,
                                mt.PoolingSettings(smd_weighting="identifiability_smd")).iloc[0]

        w_smd = parse_weights(r_smd["smd_weights"])
        w_rel = parse_weights(r_rel["smd_weights"])
        w_id = parse_weights(r_id["smd_weights"])
        for parsed in (w_smd, w_rel, w_id):
            self.assertAlmostEqual(sum(parsed.values()), 1.0, delta=0.01)
        # reliability_smd 应降低高 fit_std 的 lactate 权重，提高 glucose 权重
        self.assertGreater(w_rel["glucose"], w_smd["glucose"])
        # identifiability_smd 应降低 moderate 的 lactate 权重（q=0.5 vs 1.0）
        self.assertGreater(w_id["glucose"], w_smd["glucose"])
        self.assertEqual(r_rel["smd_weighting"], "reliability_smd")
        self.assertEqual(r_id["smd_weighting"], "identifiability_smd")
        self.assertTrue(np.isfinite(r_smd["g_total_var"]) and r_smd["g_total_var"] > 0)

        with self.assertRaises(ValueError):
            mt.PoolingSettings(smd_weighting="invalid_mode")
        with self.assertRaises(ValueError):
            mt.PoolingSettings(permutation_scale="invalid_scale")
        with self.assertRaises(ValueError):
            mt.PoolingSettings(fdr_scope="invalid_scope")

    def test_fdr_scope_per_tissue_vs_per_contrast_and_padj_global(self) -> None:
        # heart 2 条反应（p=0.01, 0.80），liver 10 条反应（p=0.80）
        def fake_row(tissue: str, rxn: str, pval: float) -> dict:
            return {
                "tissue": tissue, "reaction": rxn, "k": 2, "n_tracers": 2, "n_tracers_informative": 2,
                "tracers_informative": "glucose;lactate", "tracers_excluded": "",
                "tracer_evidence": mt.CONSISTENT, "identifiability": "well_constrained",
                "n_numerator": 6, "n_denominator": 6, "mean_numerator": 20.0, "mean_denominator": 10.0,
                "diff": 10.0, "diff_se": 2.0, "log2fc": 1.0, "log2fc_se": 0.2,
                "effect_size": 1.2, "statistic": 2.5, "permutation_pvalue": pval, "min_achievable_p": 0.005,
                "heterogeneity_q": 0.1, "heterogeneity_pvalue": 0.75, "i_squared": 0.0,
                "smd": 1.2, "smd_se": 0.4, "smd_heterogeneity_pvalue": 0.75,
                "n_tracers_up": 2, "n_tracers_down": 0, "direction_stability": 1.0, "pooled_loo_stability": 1.0,
            }

        pooled = pd.DataFrame(
            [fake_row("heart", "R_sig", 0.01), fake_row("heart", "R_null", 0.80)]
            + [fake_row("liver", f"L_{i}", 0.80) for i in range(10)]
        )
        t_tissue = mt.to_differential(pooled, contrast_id="c", contrast_role="other", flux_units=UNITS,
                                      settings=mt.PoolingSettings(fdr_scope="per_tissue"))
        t_global = mt.to_differential(pooled, contrast_id="c", contrast_role="other", flux_units=UNITS,
                                      settings=mt.PoolingSettings(fdr_scope="per_contrast"))
        r_sig_t = t_tissue[t_tissue["feature_id"] == "R_sig"].iloc[0]
        r_sig_g = t_global[t_global["feature_id"] == "R_sig"].iloc[0]
        # 按组织（heart 仅 2 检验）：0.01 * 2 = 0.02；跨全组织（12 检验）：0.01 * 12 = 0.12
        self.assertAlmostEqual(r_sig_t["padj"], 0.02, places=6)
        self.assertAlmostEqual(r_sig_t["padj_global"], 0.12, places=6)
        self.assertEqual(r_sig_t["fdr_scope"], "per_tissue")
        self.assertAlmostEqual(r_sig_g["padj"], 0.12, places=6)
        self.assertAlmostEqual(r_sig_g["padj_global"], 0.12, places=6)
        self.assertEqual(r_sig_g["fdr_scope"], "per_contrast")
        # differential_flux.fdr_scope 可被 multi_tracer 继承
        inherited = mt.PoolingSettings.from_analysis({"differential_flux": {"fdr_scope": "per_contrast"}})
        self.assertEqual(inherited.fdr_scope, "per_contrast")

    def test_relative_permutation_scale_is_invariant_to_stratum_scaling(self) -> None:
        base_strata = [
            _stratum("glucose", [12.0, 13.0, 12.5], [10.0, 10.4, 9.9], fit_var=0.2),
            _stratum("lactate", [11.0, 11.4, 10.8], [9.8, 10.1, 9.7], fit_var=0.1),
        ]
        # 把 glucose 的绝对尺度整体放大 1000 倍（fit_var 放大 10^6 倍）
        scaled_strata = [
            _stratum("glucose", [12000.0, 13000.0, 12500.0], [10000.0, 10400.0, 9900.0], fit_var=200_000.0),
            _stratum("lactate", [11.0, 11.4, 10.8], [9.8, 10.1, 9.7], fit_var=0.1),
        ]
        rel = mt.PoolingSettings(permutation_scale="relative")
        p_base = mt._permutation_test(base_strata, rel, seed=0)
        p_scaled = mt._permutation_test(scaled_strata, rel, seed=0)
        self.assertAlmostEqual(p_base["pvalue"], p_scaled["pvalue"], places=12)
        self.assertAlmostEqual(p_base["z"], p_scaled["z"], places=10)

    def test_stratified_direction_and_loo_stability(self) -> None:
        # 场景 1：两个主示踪剂强同向 + 一个弱示踪剂轻度波动（类似 heart/G3PP_c）→ 分层合并稳定
        spec_stable = {
            "glucose": ([4.0, 4.4, 4.2], [7.0, 7.5, 7.2]),
            "glycerol": ([8.0, 8.6, 8.2], [11.0, 11.5, 11.2]),
            "alanine": ([9.5, 10.5, 9.2], [10.0, 9.8, 10.4]),
        }
        animals, groups = _animals(spec_stable, fit_std=0.1)
        row_ok = mt.pool_contrast(
            animals, groups, "treat", "control",
            _coverage(dict.fromkeys(spec_stable, "well_constrained")),
            mt.PoolingSettings(n_bootstrap=800, random_state=7),
        ).iloc[0]
        self.assertEqual(row_ok["pooled_loo_stability"], 1.0)
        self.assertGreaterEqual(row_ok["direction_stability"], 0.95)
        diff_ok = mt.to_differential(pd.DataFrame([row_ok]), contrast_id="c", contrast_role="other",
                                     flux_units=UNITS, settings=mt.PoolingSettings())
        for col in ("direction_stability", "pooled_loo_stability", "g_total_var", "smd_weights", "smd_weighting",
                    "permutation_scale"):
            self.assertIn(col, diff_ok.columns)

        # 场景 2：单只极端异常值驱动全部合并信号 → 分层 LOO 与自助法稳定性均跌落
        spec_outlier = {
            "glucose": ([1.0, 1.0, 25.0], [2.0, 2.1, 1.9]),
            "lactate": ([2.0, 2.1, 1.9], [2.0, 2.05, 1.95]),
        }
        animals2, groups2 = _animals(spec_outlier, fit_std=0.05)
        row_bad = mt.pool_contrast(
            animals2, groups2, "treat", "control",
            _coverage(dict.fromkeys(spec_outlier, "moderate")),
            mt.PoolingSettings(n_bootstrap=800, random_state=7),
        ).iloc[0]
        self.assertLess(row_bad["pooled_loo_stability"], 1.0)
        self.assertLess(row_bad["direction_stability"], 0.85)

    def test_conflicting_tracers_are_not_selected(self) -> None:
        spec = {"glucose": ([16.0, 16.2, 15.9], [20.0, 20.1, 19.8]), "lactate": ([28.0, 28.3, 27.9], [20.0, 20.2, 19.9])}
        animals, groups = _animals(spec)
        pooled = mt.pool_contrast(animals, groups, "treat", "control",
                                  _coverage(dict.fromkeys(spec, "well_constrained")), mt.PoolingSettings())
        table = mt.to_differential(pooled, contrast_id="c", contrast_role="other", flux_units=UNITS,
                                   settings=mt.PoolingSettings())
        self.assertEqual(pooled.at[0, "tracer_evidence"], mt.CONFLICT)
        self.assertLess(pooled.at[0, "heterogeneity_pvalue"], 0.01)
        self.assertFalse(table.at[0, "selected"])
        self.assertIn("tracer_conflict", table.at[0, "flags"])

    def test_unknown_setting_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            mt.PoolingSettings.from_analysis({"multi_tracer": {"heterogenity_alpha": 0.1}})
        settings = mt.PoolingSettings.from_analysis({"multi_tracer": {"min_tracers": 2},
                                                     "differential_flux": {"cohen_d_threshold": 0.8}})
        self.assertEqual((settings.min_tracers, settings.hedges_g_threshold), (2, 0.8))


class MergeTest(unittest.TestCase):
    def test_merge_requires_same_model_and_unique_samples(self) -> None:
        rng = np.random.default_rng(0)
        a = tracer_results("glucose", {"g1": False, "g2": True}, rng)
        b = tracer_results("lactate", {"l1": False, "l2": True}, rng)
        merged = merge_mfa([a, b])
        self.assertEqual(merged.samples, ["g1", "g2", "l1", "l2"])
        self.assertEqual(merged.info.tracer, "13C glucose; 13C lactate")
        with self.assertRaisesRegex(ContractError, "model_id"):
            merge_mfa([a, tracer_results("lactate", {"l3": False}, rng, model_id="other_model")])
        with self.assertRaisesRegex(ContractError, "g1"):
            merge_mfa([a, a])


class LegacyInputSetsTest(unittest.TestCase):
    def test_flux_and_mid_files_are_paired_by_content(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            glc = write_legacy_mfa(root / "glucose", {"g1": "control", "g2": "treat"})
            lac = write_legacy_mfa(root / "lactate", {"l1": "control", "l2": "treat"}, seed=1)
            # 第二组只给通量文件
            results, notes = read_input_sets([glc[0], glc[1], lac[0]], "legacy_mfa_xlsx", model_id="m",
                                             flux_units=UNITS, tracer="13C tracer")
            self.assertEqual(results.samples, ["g1", "g2", "l1", "l2"])
            self.assertEqual(set(results.mid["sample_id"]), {"g1", "g2"})
            self.assertIn("2 组结果", notes[-1])
            with self.assertRaisesRegex(ContractError, "MID 文件"):
                read_input_sets([glc[1], glc[0]], "legacy_mfa_xlsx", model_id="m", flux_units=UNITS)


class MultiTracerRunTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        cls.root = root = Path(cls._tmp.name).resolve()
        build_project(root)
        cls.run_dir, cls.notes = _quiet(run_experiment, EXP, ["mfa"], root=root)
        cls.mfa = cls.run_dir / "omics" / "mfa"
        cls.table = read_differential(cls.mfa / "tables" / "differential.parquet").set_index(["tissue", "feature_id"])

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_one_pooled_row_per_reaction(self) -> None:
        self.assertEqual(len(self.table), len(REACTIONS))
        self.assertTrue((self.table["inference_unit"] == "biological_replicate").all())

    def test_pooled_table_carries_flux_groups_and_standardized_scale(self) -> None:
        # 合并结果也按通量组做 BH（G10）；标准化尺度的异质性列与各示踪剂的 g 一并写出
        for column in ("flux_group", "flux_group_size", "smd", "smd_heterogeneity_pvalue", "tracer_smd", "smd_weights"):
            self.assertIn(column, self.table.columns)
        pdh_w = self.table.loc[("liver", "PDH_m"), "smd_weights"]
        parsed = {k: float(v) for k, v in (item.split(":") for item in pdh_w.split(";"))}
        self.assertEqual(set(parsed), set(TRACERS))
        self.assertAlmostEqual(sum(parsed.values()), 1.0, delta=0.01)
        self.assertTrue((self.table["flux_group_size"] >= 1).all())
        self.assertTrue((self.mfa / "tables" / "flux_groups.parquet").is_file())
        fits = pd.read_parquet(self.mfa / "tables" / "sample_flux_fits.parquet")
        self.assertEqual(set(fits["tracer"]), set(TRACERS))

    def test_consistent_reaction_reaches_significance_after_pooling(self) -> None:
        pdh = self.table.loc[("liver", "PDH_m")]
        self.assertEqual(pdh["tracer_evidence"], mt.CONSISTENT)
        self.assertEqual(pdh["n_tracers_informative"], 4)
        self.assertEqual(pdh["evidence_level"], "statistical_test")
        self.assertLess(pdh["pvalue"], 0.001)
        self.assertTrue(pdh["selected"])
        self.assertAlmostEqual(pdh["diff"], 6.0, delta=0.6)
        self.assertEqual((pdh["n_numerator"], pdh["n_denominator"]), (12, 12))
        # 单个示踪剂 3 vs 3 只能给效应量
        by_tracer = pd.read_parquet(self.mfa / "tables" / "differential_by_tracer.parquet")
        pdh_each = by_tracer[(by_tracer["tissue"] == "liver") & (by_tracer["feature_id"] == "PDH_m")]
        self.assertEqual(sorted(pdh_each["tracer"]), sorted(TRACERS))
        self.assertTrue((pdh_each["evidence_level"] == "effect_size_only").all())

    def test_reaction_resolved_by_one_tracer(self) -> None:
        gk = self.table.loc[("liver", "GK_c")]
        self.assertEqual(gk["tracers_informative"], "glycerol")
        self.assertEqual(gk["tracer_evidence"], mt.SINGLE)
        self.assertEqual(gk["evidence_level"], "effect_size_only")
        self.assertTrue(gk["selected"])
        coverage = pd.read_parquet(self.mfa / "tables" / "tracer_coverage.parquet")
        gk_cov = coverage[coverage["reaction"] == "GK_c"].set_index("tracer")
        self.assertEqual(gk_cov.loc["glycerol", "identifiability"], "well_constrained")
        self.assertTrue((gk_cov.drop(index="glycerol")["identifiability"] == "poorly_constrained").all())
        self.assertTrue((gk_cov["n_informative_tracers"] == 1).all())

    def test_conflicting_tracer_blocks_the_conclusion(self) -> None:
        ldh = self.table.loc[("liver", "LDH_c")]
        self.assertEqual(ldh["tracer_evidence"], mt.CONFLICT)
        self.assertFalse(ldh["selected"])
        self.assertEqual((ldh["n_tracers_up"], ldh["n_tracers_down"]), (1, 3))
        # 乳酸队列处理组的 LDH 绝对水平也偏离其余示踪剂：一致性 QC 把它标出来
        concordance = pd.read_parquet(self.mfa / "qc" / "tracer_concordance.parquet")
        liver = concordance[concordance["tissue"] == "liver"].set_index("tracer")
        self.assertEqual(liver["max_abs_z"].idxmax(), "lactate")

    def test_flux_levels_and_bundle(self) -> None:
        levels = pd.read_parquet(self.mfa / "tables" / "tracer_flux_levels.parquet")
        pdh = levels[(levels["reaction"] == "PDH_m") & (levels["condition"] == "control")].iloc[0]
        self.assertEqual(pdh["n_tracers"], 4)
        self.assertAlmostEqual(pdh["pooled_flux"], 10.0, delta=0.6)
        ldh_treat = levels[(levels["reaction"] == "LDH_c") & (levels["condition"] == "treat")].iloc[0]
        self.assertFalse(ldh_treat["tracers_consistent"])

        self.assertTrue((self.mfa / "data" / "liver").is_dir())
        manifest = read_manifest(self.run_dir)
        self.assertEqual(manifest.parameters["mfa"]["tracers"], dict.fromkeys(sorted(TRACERS), 6))
        inputs = {f.path for f in manifest.files if f.role == "input"}
        self.assertTrue(any(p.startswith(f"experiments/{EXP}/mfa/lactate/") for p in inputs))


_MULTI_TRACER_OUTPUTS = ("tables/differential_by_tracer.parquet", "tables/tracer_coverage.parquet",
                         "tables/tracer_flux_levels.parquet", "qc/tracer_concordance.parquet")
GLUCOSE = [("glucose", "glucose", "U-13C glucose", "glucose")]
# 同一种示踪剂的两个批次，info.tracer 的写法不同
GLUCOSE_TWO_BATCHES = [("glucose_b1", "glucose", "U-13C glucose", "glucose"),
                       ("glucose_b2", "glucose", "13C-glucose", "glucose")]


class SingleTracerCompatibilityTest(unittest.TestCase):
    """只有一种示踪剂时走原来的单示踪剂分析：protocol 列可有可无，不产生多示踪剂的表。"""

    def _run(self, **kwargs) -> tuple[Path, dict]:
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name).resolve()
        build_project(root, **kwargs)
        run_dir, _ = _quiet(run_experiment, EXP, ["mfa"], root=root)
        return run_dir / "omics" / "mfa", read_manifest(run_dir).parameters["mfa"]

    def test_protocol_column_is_optional(self) -> None:
        with_protocol, parameters = self._run(cohorts=GLUCOSE)
        without_protocol, _ = self._run(cohorts=GLUCOSE, protocol=False)
        tables = [read_differential(d / "tables" / "differential.parquet") for d in (with_protocol, without_protocol)]
        pd.testing.assert_frame_equal(tables[0], tables[1])
        self.assertNotIn("tracer_evidence", tables[0].columns)
        self.assertTrue((tables[0]["n_numerator"] == 3).all())
        # 单示踪剂 3 vs 3 与原流程一样只给效应量
        self.assertTrue((tables[0]["evidence_level"] == "effect_size_only").all())
        self.assertNotIn("tracers", parameters)
        for directory in (with_protocol, without_protocol):
            for name in _MULTI_TRACER_OUTPUTS:
                self.assertFalse((directory / name).exists(), name)

    def test_batches_of_one_tracer_with_differently_spelled_labels(self) -> None:
        mfa, parameters = self._run(cohorts=GLUCOSE_TWO_BATCHES)
        table = read_differential(mfa / "tables" / "differential.parquet")
        self.assertTrue((table["n_numerator"] == 6).all())
        self.assertNotIn("tracers", parameters)
        self.assertFalse((mfa / "tables" / "tracer_coverage.parquet").exists())


class MultiTracerInputErrorsTest(unittest.TestCase):
    def test_differently_labelled_results_need_protocol(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            build_project(root, cohorts=GLUCOSE_TWO_BATCHES, protocol=False)
            with self.assertRaisesRegex(RunError, "同一个 protocol"):
                _quiet(run_experiment, EXP, ["mfa"], root=root)

    def test_tracers_must_be_labelled(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            build_project(root, protocol=False)
            with self.assertRaisesRegex(RunError, "protocol"):
                _quiet(run_experiment, EXP, ["mfa"], root=root)

    def test_partly_labelled_tracers_fail_validation(self) -> None:
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp).resolve()
            rows = build_project(root)
            rows[0]["protocol"] = ""
            pd.DataFrame(rows).drop(columns="group").to_csv(root / "experiments" / EXP / "samples.tsv",
                                                            sep="\t", index=False)
            with self.assertRaisesRegex(RunError, "没有填写"):
                _quiet(run_experiment, EXP, ["mfa"], root=root)


if __name__ == "__main__":
    unittest.main()
