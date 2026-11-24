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
    regulation_axis, regulation_coefficient, simes, stratified_pooled_loo,
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

    def test_stratified_pooled_stability_rescues_weak_minor_tracer(self) -> None:
        # 1. stratified_pooled_loo 直接测试：强主示踪剂 + 弱次要示踪剂 → 合并 LOO 为 1.0；n=2 vs 2 固定尺度不退化
        stable = stratified_pooled_loo([
            (np.array([4.0, 4.4, 4.2]), np.array([7.0, 7.5, 7.2])),
            (np.array([9.5, 10.5, 9.2]), np.array([10.0, 9.8, 10.4])),
        ])
        self.assertEqual(stable, 1.0)
        small_n = stratified_pooled_loo([
            (np.array([10.0, 10.2]), np.array([5.0, 5.2])),
            (np.array([8.0, 8.2]), np.array([4.0, 4.1])),
        ])
        self.assertEqual(small_n, 1.0)
        # 单只极端动物驱动 → 合并 LOO < 1.0
        outlier = stratified_pooled_loo([
            (np.array([1.0, 1.0, 25.0]), np.array([2.0, 2.1, 1.9])),
            (np.array([2.0, 2.1, 1.9]), np.array([2.0, 2.05, 1.95])),
        ])
        self.assertLess(outlier, 1.0)

        # 2. flux_axis 优先使用 __pooled__，不再因次要示踪剂 alanine 单层 LOO=0.67 而否决
        multi = frame(
            flux("glucose", -0.8, effect=-1.4),
            flux("glycerol", -0.7, effect=-1.2),
            flux("alanine", -0.05, effect=-0.12),
            flux("pooled", -0.75, effect=-1.10, pvalue=0.02, padj=0.25),
        )
        ax_rescued = flux_axis(multi, RULES, stability={"glucose": 1.0, "glycerol": 1.0, "alanine": 0.67, "__pooled__": 1.0})
        self.assertEqual((ax_rescued.state, ax_rescued.direction), ("suggestive", -1))

        # 3. pooled_loo_stability < 1.0 即使 direction_stability >= 0.8 也不通过稳定性门控
        with_loo_fail = frame(
            flux("glucose", 1.0, effect=1.4),
            flux("lactate", 0.3, effect=0.4),
            flux("pooled", 0.9, effect=1.2, direction_stability=0.85, pooled_loo_stability=0.91),
        )
        ax_loo_fail = flux_axis(with_loo_fail, RULES)
        self.assertEqual(ax_loo_fail.state, "undetermined")
        self.assertIn("LOO", ax_loo_fail.rule)

        # 4. 绝对尺度 log2fc 与标准化尺度 Hedges' g 异号时记录 scale_direction_discordant
        discordant = frame(
            flux("glucose", 1.2, effect=0.9),
            flux("lactate", -0.4, effect=-1.1),
            flux("pooled", 0.6, effect=-0.3, direction_stability=0.9),
        )
        ax_disc = flux_axis(discordant, RULES)
        self.assertIn("scale_direction_discordant", ax_disc.flags)

    def test_assess_uses_stratified_pooled_loo_end_to_end(self) -> None:
        # 端到端 assess()：仅提供 dossier_samples（无预存 direction_stability），次要示踪剂单层反号但分层合并稳定
        members = frame(
            flux("glucose", -0.8, effect=-1.4),
            flux("glycerol", -0.7, effect=-1.2),
            flux("alanine", -0.05, effect=-0.12),
            flux("pooled", -0.75, effect=-1.10, pvalue=0.02, padj=0.25),
        )
        samples_rows, side_rows = [], []
        spec = {
            "glucose": ([4.0, 4.4, 4.2], [7.0, 7.5, 7.2]),
            "glycerol": ([8.0, 8.6, 8.2], [11.0, 11.5, 11.2]),
            "alanine": ([9.5, 10.5, 9.2], [10.0, 9.8, 10.4]),
        }
        for tr, (num, den) in spec.items():
            for side, vals in (("numerator", num), ("denominator", den)):
                for i, v in enumerate(vals):
                    sid = f"{tr}_{side}_{i}"
                    samples_rows.append({"source": "A", "layer": "mfa", "tissue": "liver",
                                         "feature_id": "R1", "tracer": tr, "sample_id": sid, "value": v})
                    side_rows.append({"source": "A", "layer": "mfa", "contrast_id": "C",
                                      "sample_id": sid, "side": side})
        axes_tbl, rxn_tbl = assess(members, "A", RULES, samples=pd.DataFrame(samples_rows),
                                   sides=pd.DataFrame(side_rows))
        f_row = axes_tbl[axes_tbl["axis"] == "flux"].iloc[0]
        self.assertEqual((f_row["state"], int(f_row["direction"])), ("suggestive", -1))
        self.assertEqual(rxn_tbl.iloc[0]["flux_claim"], "L2")

    def test_phase3_multi_tracer_perm_suggestive_path_b(self) -> None:
        # S4 路径 B：≥2 可辨识同向示踪剂、无异质性、置换 p < 0.05、|g| = 0.85 >= 0.8、方向稳定 → suggestive
        concordant = frame(
            flux("glucose", -1.4, effect=-0.9),
            flux("glycerol", -1.8, effect=-0.8),
            flux("pooled", -1.6, effect=-0.85, pvalue=0.033, padj=0.35, direction_stability=0.96,
                 pooled_loo_stability=1.0),
        )
        ax = flux_axis(concordant, RULES)
        self.assertEqual((ax.state, ax.direction), ("suggestive", -1))
        self.assertIn("multi_tracer_perm_suggestive", ax.flags)

        # 关闭路径 B（flux_concordant_min_abs_g=None）回退为 undetermined
        ax_off = flux_axis(concordant, {**RULES, "flux_concordant_min_abs_g": None})
        self.assertEqual(ax_off.state, "undetermined")

        # 仅 1 种示踪剂 / 有异质性 / |g| < 0.8 均不能走路径 B
        single = frame(
            flux("glucose", -1.4, effect=-0.85),
            flux("pooled", -1.4, effect=-0.85, pvalue=0.033, direction_stability=0.96, pooled_loo_stability=1.0),
        )
        self.assertEqual(flux_axis(single, RULES).state, "undetermined")
        het = frame(
            flux("glucose", -1.4, effect=-0.9),
            flux("glycerol", -1.8, effect=-0.8),
            flux("pooled", -1.6, effect=-0.85, pvalue=0.033, heterogeneity_pvalue=0.04,
                 direction_stability=0.96, pooled_loo_stability=1.0),
        )
        self.assertEqual(flux_axis(het, RULES).state, "undetermined")
        weak_g = frame(
            flux("glucose", -1.4, effect=-0.7),
            flux("glycerol", -1.8, effect=-0.7),
            flux("pooled", -1.6, effect=-0.72, pvalue=0.033, direction_stability=0.96, pooled_loo_stability=1.0),
        )
        self.assertEqual(flux_axis(weak_g, RULES).state, "undetermined")

    def test_phase3_flux_tost_equivalence_and_l1_guard(self) -> None:
        # S5：4 示踪剂均无强效应、pooled 为 equivalent、无异质性且 |g| 90% 上界 < 1.0 → unchanged
        equiv_rows = frame(
            flux("glucose", 0.20, effect=0.35, n_numerator=3, n_denominator=3),
            flux("glycerol", -0.15, effect=-0.25, n_numerator=3, n_denominator=3),
            flux("lactate", -0.10, effect=-0.20, n_numerator=3, n_denominator=3),
            flux("alanine", 0.25, effect=0.40, n_numerator=3, n_denominator=3),
            flux("pooled", 0.08, state="equivalent", se=0.18, effect=0.10, smd_se=0.35,
                 n_numerator=12, n_denominator=12, direction_stability=0.60, pooled_loo_stability=0.80),
        )
        ax_eq = flux_axis(equiv_rows, RULES)
        self.assertEqual((ax_eq.state, ax_eq.direction), ("unchanged", 0))
        self.assertIn("通量等效未变", ax_eq.rule)

        # 护栏 1：若任一可辨识示踪剂自身 |g| >= 1.0（类似 RPI_c 陷阱），不得判为 unchanged
        trap_rows = equiv_rows.copy()
        trap_rows.loc[trap_rows["tracer"] == "alanine", "effect"] = 1.45
        self.assertEqual(flux_axis(trap_rows, RULES).state, "undetermined")

        # 护栏 2：样本量不足（单示踪剂 3v3，|g| 90% 上界 > 1.0）不得判为 unchanged
        low_power = frame(
            flux("glucose", 0.05, effect=0.10, n_numerator=3, n_denominator=3),
            flux("pooled", 0.05, state="equivalent", se=0.10, effect=0.10, n_numerator=3, n_denominator=3),
        )
        self.assertEqual(flux_axis(low_power, RULES).state, "undetermined")

        # 护栏 3：关闭 S5（flux_equiv_g_bound=None）回退为 undetermined
        self.assertEqual(flux_axis(equiv_rows, {**RULES, "flux_equiv_g_bound": None}).state, "undetermined")

        # 护栏 4：decide() 当 F 轴为 unchanged 且 C 轴 changed 时，不再给 L1，而是 L0 + flux_unchanged_despite_capacity
        decision = decide(ax_eq, Axis("capacity", "changed", direction=1),
                          Axis("regulation", "no_knowledge"), Axis("pool", "not_measured"), RULES)
        self.assertEqual((decision["flux_claim"], decision["claim_direction"]), ("L0", 0))
        self.assertIn("flux_unchanged_despite_capacity", decision["phenomena"])

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

    def test_not_assessable_reaction_cannot_be_promoted_above_l0(self) -> None:
        # 虚拟/边界反应（如 BIOMASS_REACTION_1、Salvage_c）：即使通量轴达到 suggestive / changed，也禁止升为 L1/L2/L3/LX
        biomass = frame(
            flux("glucose", 1.5, effect=1.8, direction_stability=0.95, flux_id="BIOMASS_REACTION_1",
                 direction_evidence="not_assessable"),
            flux("pooled", 1.5, effect=1.8, padj=0.01, pvalue=0.001, direction_stability=0.95,
                 flux_id="BIOMASS_REACTION_1", direction_evidence="not_assessable"),
        )
        _, rxn = assess(biomass, "A", RULES)
        row_out = rxn.iloc[0]
        self.assertEqual(row_out["flux_claim"], "L0")
        self.assertEqual(int(row_out["claim_direction"]), 0)
        self.assertEqual(int(row_out["n_concordant_axes"]), 0)
        self.assertEqual(row_out["mechanism_class"], "")
        self.assertEqual(row_out["phenomena"], "not_assessable_reaction")

        # decide() 直接验证：即使传入 changed F 轴 + changed C 轴，not_assessable 也强制 L0
        dec = decide(Axis("flux", "changed", direction=1), Axis("capacity", "changed", direction=1),
                     Axis("regulation", "no_knowledge"), Axis("pool", "not_measured"), RULES,
                     semantics="not_assessable")
        self.assertEqual((dec["flux_claim"], dec["claim_direction"], dec["phenomena"]),
                         ("L0", 0, "not_assessable_reaction"))

    def test_pathway_step_active_share_and_combined_concordant_rescue(self) -> None:
        # 1. pathway_step（如 THERMOGENESIS_UCP 的 Ucp1、FA_esterification_c 的 Mogat1）：不受可替代同工酶丰度占比压制
        mogat = frame(
            changed(3.2, layer="transcriptomics", member_id="Mogat1", member_relation="pathway_step",
                    abundance_share=0.004),
            row(layer="transcriptomics", member_id="Lpin2", member_relation="pathway_step",
                log2fc=0.05, abundance_share=0.40, pvalue=0.7),
            row(layer="transcriptomics", member_id="Dgat2", member_relation="pathway_step",
                log2fc=0.10, abundance_share=0.35, pvalue=0.6),
        )
        ax_step = capacity_axis(mogat, {**RULES, "rna_single_source_changed": True}, ["A"])
        self.assertEqual((ax_step.state, ax_step.direction), ("changed", 1))
        self.assertNotIn("minor_isozyme_driven", ax_step.flags)

        # pathway_step 中若两个通过 FDR 的步骤方向相反 -> mixed
        mixed_step = frame(
            changed(2.0, member_id="StepA", member_relation="pathway_step", abundance_share=0.1),
            changed(-1.8, member_id="StepB", member_relation="pathway_step", abundance_share=0.9),
        )
        self.assertEqual(capacity_axis(mixed_step, RULES, ["A"]).state, "mixed")

        # 2. 条件丰度占比 active_abundance_share：在活跃条件达到 >=20% 的诱导型同工酶得以保留
        inducible = frame(
            changed(4.0, member_id="IsoB", abundance_share=0.06, active_abundance_share=0.28),
            row(member_id="IsoA", log2fc=0.0, abundance_share=0.94, active_abundance_share=0.72, pvalue=0.9),
        )
        ax_ind = capacity_axis(inducible, RULES, ["A"])
        self.assertEqual((ax_ind.state, ax_ind.direction), ("changed", 1))
        self.assertNotIn("minor_isozyme_driven", ax_ind.flags)

        # 真正痕量且与主要同工酶反向的诱饵仍被拦截
        decoy = frame(
            changed(1.2, member_id="Pkm", abundance_share=0.05, active_abundance_share=0.08),
            row(member_id="Pklr", log2fc=-0.6, abundance_share=0.95, active_abundance_share=0.92, pvalue=0.15),
        )
        ax_decoy = capacity_axis(decoy, RULES, ["A"])
        self.assertEqual(ax_decoy.state, "undetermined")
        self.assertIn("minor_isozyme_driven", ax_decoy.flags)

        # 3. 多成员同向合计占比达标（无主要同工酶反向）：不被 minor_isozyme_driven 误杀
        multi_ok = frame(
            changed(1.1, member_id="Iso1", abundance_share=0.12),
            changed(0.9, member_id="Iso2", abundance_share=0.11),
            row(member_id="IsoMain", log2fc=0.1, abundance_share=0.77, pvalue=0.5),
        )
        ax_multi = capacity_axis(multi_ok, RULES, ["A"])
        self.assertEqual((ax_multi.state, ax_multi.direction), ("changed", 1))
        self.assertNotIn("minor_isozyme_driven", ax_multi.flags)

        # 若主要同工酶明显反向（|log2fc| >= 0.3），即使多个次要同工酶合计 >= 20% 仍禁止放行
        multi_opp = frame(
            changed(1.1, member_id="Iso1", abundance_share=0.12),
            changed(0.9, member_id="Iso2", abundance_share=0.11),
            row(member_id="IsoMain", log2fc=-0.5, abundance_share=0.77, pvalue=0.12),
        )
        ax_opp = capacity_axis(multi_opp, RULES, ["A"])
        self.assertEqual(ax_opp.state, "undetermined")
        self.assertIn("minor_isozyme_driven", ax_opp.flags)

    def test_flux_scale_drift_gates_flux_axis_and_propagates_to_phenomena(self) -> None:
        drifted = frame(
            flux("glucose", 1.2, effect=1.8, direction_stability=1.0, flags="flux_scale_drift"),
            flux("pooled", 1.2, effect=1.8, direction_stability=1.0, flags="flux_scale_drift"),
        )
        # 默认 flux_scale_drift_gate=True：降级为 undetermined，阻止升 L2
        ax_gated = flux_axis(drifted, RULES)
        self.assertEqual(ax_gated.state, "undetermined")
        self.assertIn("flux_scale_drift", ax_gated.flags)
        dec = decide(ax_gated, Axis("capacity", "undetermined"),
                     Axis("regulation", "no_knowledge"), Axis("pool", "not_measured"), RULES)
        self.assertEqual(dec["flux_claim"], "L0")
        self.assertIn("flux_scale_drift", dec["phenomena"])

        # 关闭门控时恢复 suggestive，但仍保留 flux_scale_drift 标记提醒
        ax_open = flux_axis(drifted, {**RULES, "flux_scale_drift_gate": False})
        self.assertEqual(ax_open.state, "suggestive")
        self.assertIn("flux_scale_drift", ax_open.flags)

    def test_foreign_tissue_marker_and_qc_caution_soft_cap_capacity_and_block_l1(self) -> None:
        single_rna_rules = {**RULES, "rna_single_source_changed": True}
        for flag in ("foreign_tissue_marker", "qc_caution"):
            with self.subTest(flag=flag):
                rna = frame(changed(1.8, layer="transcriptomics", flags=flag))
                ax_c = capacity_axis(rna, single_rna_rules, ["A"])
                self.assertEqual((ax_c.state, ax_c.direction), ("suggestive", 1))
                self.assertIn(flag, ax_c.flags)
                dec = decide(Axis("flux", "not_measured"), ax_c,
                             Axis("regulation", "no_knowledge"), Axis("pool", "not_measured"), single_rna_rules)
                self.assertEqual(dec["flux_claim"], "L0")
                self.assertIn(flag, dec["phenomena"])

        # 若同时存在未标记的干净主要同工酶，则不被单个外源标记误降级
        mixed_clean = frame(
            changed(1.5, layer="transcriptomics", member_id="CleanIso", abundance_share=0.6, flags=""),
            changed(2.2, layer="transcriptomics", member_id="ForeignIso", abundance_share=0.4, flags="foreign_tissue_marker"),
        )
        ax_clean = capacity_axis(mixed_clean, single_rna_rules, ["A"])
        self.assertEqual((ax_clean.state, ax_clean.direction), ("changed", 1))
        self.assertNotIn("foreign_tissue_marker", ax_clean.flags)


if __name__ == "__main__":
    unittest.main()
