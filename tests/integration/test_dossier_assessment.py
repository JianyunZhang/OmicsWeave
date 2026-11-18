#!/usr/bin/env python3
"""分轴综合评估（S3，G15）。

验收（反应档案与分轴综合评估，第 5 节 G15 与 3.11 的不变性测试）：
- 翻转通量的符号，C、R、M 轴不变；删除一整层，其余轴不变；
- 通量上升而酶量下降的反应被归为“拮抗”，不出现在 L1；
- 没有通量证据的反应最高为 L1。
"""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from omicsweave.integration.dossier.assessment import (
    Axis, assess, bh, capacity_axis, decide, flux_axis, loo_stability, mechanism, one_sided, pool_axis,
    regulation_axis, regulation_coefficient, simes,
)

RULES = {"alpha": 0.05, "equivalence_bound": 0.58, "flux_min_abs_g": 1.0, "heterogeneity_alpha": 0.10,
         "direction_stability_min": 0.8, "established_min_group": 3, "minor_share": 0.2, "bootstrap": 400,
         "seed": 1, "opposition_min_state": "suggestive"}


def row(**fields) -> dict:
    base = {"flux_id": "R1", "tissue": "liver", "contrast_id": "C", "source": "A", "layer": "proteomics",
            "role": "direct_enzyme", "member_id": "E1", "feature_id": "E1", "feature_name": "E1", "tracer": "",
            "member_relation": "isozyme", "regulator_effect": "", "regulator_mode": "", "log2fc": np.nan,
            "se": np.nan, "pvalue": np.nan, "padj": np.nan, "lfsr": np.nan, "n_numerator": 3, "n_denominator": 3,
            "abundance_share": np.nan, "cell_state": "no_evidence", "selected": False, "identifiability": "",
            "flags": "", "direction_evidence": "directional", "effect": np.nan, "heterogeneity_pvalue": np.nan,
            "direction_stability": np.nan}
    return {**base, **fields}


def frame(*rows: dict) -> pd.DataFrame:
    return pd.DataFrame(list(rows))


def flux(tracer: str, fc: float, state: str = "no_evidence", ident: str = "moderate", **kw) -> dict:
    return row(layer="mfa", role="flux", member_id="R1", feature_id="R1", tracer=tracer, log2fc=fc,
               identifiability=ident, cell_state=state, **kw)


def changed(fc: float, **kw) -> dict:
    return row(log2fc=fc, cell_state="changed", padj=0.01, pvalue=0.001, se=0.2, **kw)


# 一个典型反应：两种示踪剂同向且合并估计通过检验；蛋白上调；抑制性激酶 RNA 下降；产物上升
TYPICAL = [
    flux("glucose", 1.0), flux("glycerol", 0.8),
    flux("pooled", 0.9, state="changed", padj=0.01, pvalue=0.001, se=0.3, effect=1.6, direction_stability=0.95,
         n_numerator=4, n_denominator=4),
    changed(0.9, member_id="E1", abundance_share=0.8), row(member_id="E2", log2fc=0.1, abundance_share=0.2,
                                                           cell_state="no_evidence", pvalue=0.6),
    changed(-1.0, layer="transcriptomics", role="regulator", member_id="K1", regulator_effect="inhibits",
            regulator_mode="phosphorylation"),
    changed(0.7, layer="metabolomics", role="product", member_id="P"),
    row(layer="metabolomics", role="substrate", member_id="S", log2fc=-0.1, pvalue=0.7),
]


def axes(rows: pd.DataFrame, expected=None) -> dict[str, Axis]:
    return {"flux": flux_axis(rows[rows["layer"] == "mfa"], RULES, expected),
            "capacity": capacity_axis(rows, RULES, ["A", "B"]), "regulation": regulation_axis(rows, RULES),
            "pool": pool_axis(rows, RULES, "directional")}


def summary(axis: Axis) -> tuple:
    return axis.state, axis.direction, round(axis.effect, 6) if pd.notna(axis.effect) else None


class AxisTest(unittest.TestCase):
    def test_typical_reaction(self) -> None:
        result = axes(frame(*TYPICAL))
        self.assertEqual(summary(result["flux"])[:2], ("changed", 1))         # established
        self.assertEqual(summary(result["capacity"])[:2], ("changed", 1))
        self.assertEqual(summary(result["regulation"])[:2], ("suggestive", 1))  # 抑制因子下降 → 活性升，最多 suggestive
        self.assertEqual(summary(result["pool"])[:2], ("changed", 1))
        decision = decide(*result.values(), RULES)
        self.assertEqual((decision["flux_claim"], decision["claim_direction"], decision["n_concordant_axes"]),
                         ("L3", 1, 3))

    def test_flipping_flux_leaves_other_axes(self) -> None:
        rows = frame(*TYPICAL)
        flipped = rows.copy()
        is_flux = flipped["layer"] == "mfa"
        flipped.loc[is_flux, "log2fc"] = -flipped.loc[is_flux, "log2fc"]
        before, after = axes(rows), axes(flipped)
        self.assertEqual(after["flux"].direction, -before["flux"].direction)
        for name in ("capacity", "regulation", "pool"):
            self.assertEqual(summary(after[name]), summary(before[name]), name)

    def test_removing_a_layer_leaves_other_axes(self) -> None:
        rows = frame(*TYPICAL)
        before = axes(rows)
        for layer, untouched in (("mfa", ("capacity", "regulation", "pool")),
                                 ("metabolomics", ("flux", "capacity", "regulation")),
                                 ("proteomics", ("flux", "regulation", "pool"))):
            after = axes(rows[rows["layer"] != layer])
            for name in untouched:
                with self.subTest(removed=layer, axis=name):
                    self.assertEqual(summary(after[name]), summary(before[name]))

    def test_flux_states(self) -> None:
        conflicting = frame(flux("glucose", 1.0, effect=1.5), flux("glycerol", -0.9, effect=-1.2), flux("pooled", 0.1))
        self.assertEqual(flux_axis(conflicting, RULES).state, "mixed")
        # 接近 0 的效应不算方向：不判为冲突（真实数据中 RPI_c 的 −0.00 对 +0.22）
        near_zero = frame(flux("glucose", -0.001, effect=-0.01), flux("glycerol", 0.22, effect=0.3), flux("pooled", 0.1))
        self.assertNotEqual(flux_axis(near_zero, RULES).state, "mixed")
        # 同号、只有绝对尺度上的异质性显著（真实数据中心脏 GAPD_c 的 +7.53 对 +0.94）：不判冲突，也不能判为 established
        heterogeneous = frame(flux("glucose", 7.5, effect=3.0), flux("glycerol", 0.9, effect=1.1),
                              flux("pooled", 1.2, padj=0.01, pvalue=0.001, effect=2.0, heterogeneity_pvalue=0.04,
                                   direction_stability=1.0, n_numerator=4, n_denominator=4))
        axis = flux_axis(heterogeneous, RULES)
        self.assertEqual((axis.state, "tracer_heterogeneity" in axis.flags), ("suggestive", True))
        low_n = frame(flux("glucose", 1.0, state="insufficient_n"), flux("pooled", 1.0, state="insufficient_n"))
        self.assertEqual(flux_axis(low_n, RULES).state, "insufficient_n")
        pinned = frame(flux("glucose", 1.0, state="not_identifiable", ident="bound_pinned"),
                       flux("lactate", 0.0, state="insufficient_n"), flux("pooled", 1.0, state="not_identifiable"))
        self.assertEqual(flux_axis(pinned, RULES).state, "not_identifiable")
        suggestive = frame(flux("glycerol", 1.0), flux("pooled", 1.0, effect=1.3, direction_stability=0.9))
        self.assertEqual(flux_axis(suggestive, RULES).state, "suggestive")
        unstable = frame(flux("glycerol", 1.0), flux("pooled", 1.0, effect=1.3, direction_stability=0.5))
        self.assertEqual(flux_axis(unstable, RULES).state, "undetermined")
        # 差异结果没有方向稳定性时，用逐只动物的值算；算不出时不能判为 suggestive
        unknown = frame(flux("glycerol", 1.0), flux("pooled", 1.0, effect=1.3))
        self.assertEqual(flux_axis(unknown, RULES).state, "undetermined")
        self.assertIn("stability_unknown", flux_axis(unknown, RULES).flags)
        self.assertEqual(flux_axis(unknown, RULES, stability={"glycerol": 1.0}).state, "suggestive")
        self.assertEqual(flux_axis(unknown, RULES, stability={"glycerol": 0.6}).state, "undetermined")
        self.assertEqual(flux_axis(frame(flux("alanine", np.nan, state="not_infused")), RULES).state, "not_measured")

    def test_tracer_expectation(self) -> None:
        rows = frame(flux("glycerol", 1.0), flux("lactate", 0.0, state="insufficient_n"),
                     flux("pooled", 1.0, effect=1.3, direction_stability=0.9))
        # 预期 lactate（样本不足，不算不符）与 glycerol：相符
        self.assertNotIn("tracer_expectation_mismatch", flux_axis(rows, RULES, ["glycerol", "lactate"]).flags)
        # 只预期 lactate、alanine，而约束它的是 glycerol：不符，最高为 undetermined
        mismatch = flux_axis(rows, RULES, ["lactate", "alanine"])
        self.assertIn("tracer_expectation_mismatch", mismatch.flags)
        self.assertEqual(mismatch.state, "undetermined")

    def test_capacity_rules(self) -> None:
        # 次要同工酶上调、主要同工酶下调：主要成员决定方向
        pk = frame(changed(0.64, member_id="Pkm", abundance_share=0.07), changed(-0.88, member_id="Pklr",
                                                                                abundance_share=0.93))
        self.assertEqual(summary(capacity_axis(pk, RULES, ["A"]))[:2], ("changed", -1))
        # 两个主要同工酶方向相反：同工酶转换
        switch = frame(changed(1.0, member_id="X", abundance_share=0.5), changed(-1.0, member_id="Y",
                                                                               abundance_share=0.5))
        self.assertEqual(capacity_axis(switch, RULES, ["A"]).state, "mixed")
        # 只有次要同工酶通过 FDR：不判变化
        minor = frame(changed(1.0, member_id="X", abundance_share=0.1),
                      row(member_id="Y", log2fc=0.05, abundance_share=0.9, pvalue=0.8))
        self.assertIn("minor_isozyme_driven", capacity_axis(minor, RULES, ["A"]).flags)
        # 复合体：核心亚基半数以上同向且有一个通过 FDR；共用亚基不参与
        complex_ = frame(changed(0.6, member_id="A1", member_relation="complex_subunit"),
                         row(member_id="A2", log2fc=0.3, member_relation="complex_subunit", pvalue=0.3),
                         row(member_id="A3", log2fc=-0.1, member_relation="complex_subunit", pvalue=0.8),
                         changed(-2.0, member_id="Dld", member_relation="shared_subunit"))
        self.assertEqual(summary(capacity_axis(complex_, RULES, ["A"]))[:2], ("changed", 1))
        # 没有蛋白：RNA 在两个来源同向才到 changed
        rna = frame(changed(1.0, layer="transcriptomics", source="A"), changed(0.8, layer="transcriptomics", source="B"))
        axis = capacity_axis(rna, RULES, ["A", "B"])
        self.assertEqual((axis.state, axis.replicated), ("changed", True))
        self.assertIn("rna_only", axis.flags)
        one = capacity_axis(frame(changed(1.0, layer="transcriptomics")), RULES, ["A"])
        self.assertEqual(one.state, "suggestive")
        # 蛋白等效、RNA 上调：容量取蛋白，标记 RNA 与蛋白不一致
        mixed_layers = frame(row(log2fc=0.02, se=0.1, cell_state="equivalent"), changed(1.2, layer="transcriptomics"))
        axis = capacity_axis(mixed_layers, RULES, ["A"])
        self.assertEqual(axis.state, "unchanged")
        self.assertIn("rna_protein_discordant", axis.flags)

    def test_single_source_rna_switch(self) -> None:
        on = {**RULES, "rna_single_source_changed": True}
        single = frame(changed(1.0, layer="transcriptomics"))
        # 默认关闭：单一来源最多 suggestive
        self.assertEqual(capacity_axis(single, {**RULES, "rna_single_source_changed": False}, ["A"]).state, "suggestive")
        # 打开：单一来源通过 FDR 即 changed，标 rna_single_source，不算复现
        axis = capacity_axis(single, on, ["A"])
        self.assertEqual((axis.state, axis.direction, axis.replicated), ("changed", 1, False))
        self.assertTrue({"rna_only", "rna_single_source"} <= set(axis.flags))
        # 第二个来源测到 RNA 而没有复现：开关不适用
        unreplicated = frame(changed(1.0, layer="transcriptomics", source="A"),
                             row(layer="transcriptomics", source="B", log2fc=0.1, pvalue=0.7))
        self.assertEqual(capacity_axis(unreplicated, on, ["A", "B"]).state, "suggestive")
        # 第二个来源有该层、但成员没有检出：只有一个来源测到，开关适用
        undetected = frame(changed(1.0, layer="transcriptomics", source="A"),
                           row(layer="transcriptomics", source="B", cell_state="not_detected"))
        self.assertEqual(capacity_axis(undetected, on, ["A", "B"]).state, "changed")
        # 由单只动物驱动的不放宽
        driven = frame(changed(1.0, layer="transcriptomics", flags="single_subject_driven"))
        self.assertEqual(capacity_axis(driven, on, ["A"]).state, "suggestive")
        # 只有名义 p：仍为 suggestive
        nominal = frame(row(layer="transcriptomics", log2fc=1.0, pvalue=0.01, padj=0.2, cell_state="nominal"))
        self.assertEqual(capacity_axis(nominal, on, ["A"]).state, "suggestive")
        # 有蛋白时以蛋白为准，开关不起作用
        protein = frame(row(log2fc=0.02, se=0.1, cell_state="equivalent"), changed(1.2, layer="transcriptomics"))
        self.assertEqual(capacity_axis(protein, on, ["A"]).state, "unchanged")
        # 判定表：没有通量证据时成为 L1，现象标签写明依据
        decision = decide(Axis("flux", "not_identifiable"), axis, Axis("regulation", "no_knowledge"),
                          Axis("pool", "not_measured"), on)
        self.assertEqual((decision["flux_claim"], decision["claim_direction"]), ("L1", 1))
        self.assertIn("rna_single_source", decision["phenomena"])

    def test_regulation_rules(self) -> None:
        unsigned = frame(changed(1.0, role="regulator", regulator_effect="unknown", regulator_mode="signaling"))
        self.assertEqual(regulation_axis(unsigned, RULES).state, "no_knowledge")
        effector = frame(changed(-1.0, layer="metabolomics", role="effector", regulator_effect="activates",
                                 regulator_mode="allosteric"))
        self.assertEqual(summary(regulation_axis(effector, RULES))[:2], ("changed", -1))   # 测得的效应物可到 changed

    def test_claims(self) -> None:
        typical = axes(frame(*TYPICAL))
        no_flux = axes(frame(*[r for r in TYPICAL if r["layer"] != "mfa"]))
        decision = decide(*no_flux.values(), RULES)
        self.assertEqual(decision["flux_claim"], "L1")                     # 没有通量证据：最高 L1
        # 通量上升、酶量下降：L3，机制为拮抗（不在 L1）
        down = frame(*[r for r in TYPICAL if r["layer"] != "proteomics"],
                     changed(-0.9, member_id="E1", abundance_share=1.0))
        result = axes(down)
        decision = decide(*result.values(), RULES)
        self.assertEqual(decision["flux_claim"], "L3")
        self.assertIn("capacity", decision["discordant_axes"])
        self.assertEqual(mechanism("L3", result["flux"], result["capacity"], (-1.4, -0.6)), "antagonistic")
        self.assertEqual(mechanism("L3", result["flux"], typical["capacity"], (0.4, 1.3)), "enzyme_dominant")
        self.assertEqual(mechanism("L1", result["flux"], typical["capacity"], (0.4, 1.3)), "")
        # 池与容量反向：默认 suggestive 即算反向，改为 changed 后不算
        opposed = frame(*[r for r in TYPICAL if r["layer"] not in ("mfa", "metabolomics")],
                        row(layer="metabolomics", role="product", member_id="P", log2fc=-0.6, pvalue=0.02,
                            padj=0.3, cell_state="nominal"))
        result = axes(opposed)
        self.assertEqual(decide(*result.values(), RULES)["flux_claim"], "L0")
        self.assertEqual(decide(*result.values(), {**RULES, "opposition_min_state": "changed"})["flux_claim"], "L1")

    def test_leave_one_out_stability(self) -> None:
        self.assertEqual(loo_stability([3.0, 3.2, 2.9], [1.0, 1.1, 0.9]), 1.0)
        # 一只动物决定方向：去掉它后反号（6 个估计中 1 个反号）；由逐只动物算得时要求没有反号
        driven = loo_stability([1.0, 1.0, 10.0], [2.0, 2.1, 1.9])
        self.assertAlmostEqual(driven, 5 / 6)
        rows = frame(flux("glycerol", 1.0), flux("pooled", 1.0, effect=1.3))
        self.assertEqual(flux_axis(rows, RULES, stability={"glycerol": driven}).state, "undetermined")
        self.assertTrue(np.isnan(loo_stability([1.0], [2.0, 2.1])))

    def test_statistics(self) -> None:
        self.assertAlmostEqual(one_sided(0.04, 1.0, 1), 0.02)
        self.assertAlmostEqual(one_sided(0.04, -1.0, 1), 0.98)
        self.assertAlmostEqual(simes([0.01, 0.04, 0.5]), 0.03)
        np.testing.assert_allclose(bh(pd.Series([0.01, 0.04, 0.03])), [0.03, 0.04, 0.04])

    def test_regulation_coefficient(self) -> None:
        """通量加倍（Δln J = ln 2）、酶量减半（Δln E = −ln 2）：ρ_h ≈ −1。"""
        rng = np.random.default_rng(0)
        flux_samples = pd.DataFrame({"sample_id": [f"f{i}" for i in range(8)],
                                     "value": np.r_[20 + rng.normal(0, 0.5, 4), 10 + rng.normal(0, 0.5, 4)]})
        enzyme = pd.DataFrame({"sample_id": [f"p{i}" for i in range(8)], "feature_id": "E1", "layer": "proteomics",
                               "value": np.r_[9 + rng.normal(0, 0.05, 4), 10 + rng.normal(0, 0.05, 4)]})
        sides = {**{("mfa", f"f{i}"): "numerator" if i < 4 else "denominator" for i in range(8)},
                 **{("proteomics", f"p{i}"): "numerator" if i < 4 else "denominator" for i in range(8)}}
        lower, upper = regulation_coefficient(flux_samples, enzyme, sides, RULES, {"E1": 1.0})
        self.assertTrue(lower < -0.8 and upper < 0 and lower < -1 < upper + 0.3)


class AssessTest(unittest.TestCase):
    def test_tables_and_ranking(self) -> None:
        typical = frame(*TYPICAL)
        no_flux = frame(*[r for r in TYPICAL if r["layer"] != "mfa"]).assign(flux_id="R2")
        nothing = frame(row(cell_state="not_measured", log2fc=np.nan)).assign(flux_id="R3")
        members = pd.concat([typical, no_flux, nothing], ignore_index=True)
        axes_table, reactions = assess(members, "A", RULES)
        self.assertEqual(len(axes_table), 3 * 4)
        claims = reactions.set_index("flux_id")
        self.assertEqual(list(claims.loc[["R1", "R2", "R3"], "flux_claim"]), ["L3", "L1", "L0"])
        self.assertEqual(list(claims.loc[["R1", "R2", "R3"], "rank"]), [1, 2, 3])
        self.assertEqual(claims.at["R1", "conjunction_statement"], "flux up & capacity up & regulation up & pool up")
        components = axes_table[axes_table["flux_id"] == "R1"]["pvalue_one_sided"]
        self.assertAlmostEqual(claims.at["R1", "conjunction_pvalue"], components.max())
        self.assertEqual(claims.at["R2", "conjunction_statement"], "capacity up & regulation up & pool up")


if __name__ == "__main__":
    unittest.main()
