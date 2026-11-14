#!/usr/bin/env python3
"""证据等级与锚点状态挂钩（G5）：锚点不成立时不出现有方向含义的锚定等级，改用组学之间一致性的等级。"""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from omicsweave.integration.methods.evidence_profile import (
    ANCHORED_GRADES, ANCHORED_STATES, OMICS_ONLY_GRADES, EvidenceProfileScorer,
)
from omicsweave.integration.modes.flux_anchored import (
    metabolomics_frames, mfa_frames, proteomics_frames, transcriptomics_frames,
)
from tests.integration.test_flux_anchored import DIET, DRUG, build, contract_tables


def profile_from(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    profile, _ = build(mfa_frames(tables["mfa"]), transcriptomics_frames(tables["transcriptomics"]),
                       proteomics_frames(tables["proteomics"]), metabolomics_frames(tables["metabolomics"]))
    return profile.set_index(["comparison", "tissue", "flux_id"]).sort_index()


def without_anchor(tables: dict[str, pd.DataFrame], flip: bool = False) -> dict[str, pd.DataFrame]:
    """MFA 全部不显著（锚点不成立）；``flip`` 时把通量方向全部反过来。"""
    mfa = tables["mfa"].copy()
    mfa["padj"] = mfa["padj"].where(mfa["padj"].isna(), 0.6)
    if flip:
        mfa["log2fc"] = -mfa["log2fc"]
    return {**tables, "mfa": mfa}


class AnchorGatingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tables = contract_tables()
        cls.profile = profile_from(cls.tables)

    def test_anchor_states(self) -> None:
        states = self.profile["anchor_state"]
        self.assertEqual(states[(DIET, "liver", "HEX_c")], "established")
        self.assertEqual(states[(DIET, "liver", "G6Pase_c")], "not_established")
        # 药物比较每组 1 只：即使效应量判定通过，锚点也证据不足
        self.assertTrue(self.profile.loc[(DRUG, "liver", "HEX_c"), "gate_passed"])
        self.assertEqual(states[(DRUG, "liver", "HEX_c")], "insufficient")
        self.assertEqual(set(self.profile["anchor_type"]), {"mfa_flux"})

    def test_no_anchored_grade_without_anchor(self) -> None:
        for key, row in self.profile.iterrows():
            with self.subTest(key=key):
                if row["anchor_state"] in ANCHORED_STATES:
                    self.assertNotIn(row["evidence_grade"], OMICS_ONLY_GRADES)
                else:
                    self.assertNotIn(row["evidence_grade"], ANCHORED_GRADES)
        # 疾病比较中葡萄糖-6-磷酸酶：RNA 升（提示通量升），底物 G6P 升、产物葡萄糖降（提示通量降）
        self.assertEqual(self.profile.loc[(DIET, "liver", "G6Pase_c"), "evidence_grade"], "Omics-discordant")
        self.assertEqual(self.profile.loc[(DRUG, "liver", "HEX_c"), "evidence_grade"], "Omics-single-layer")

    def test_unanchored_grades_do_not_depend_on_anchor_sign(self) -> None:
        """不变性：把锚点的符号翻转后，无锚点等级与组学方向不变。"""
        plain = profile_from(without_anchor(self.tables))
        flipped = profile_from(without_anchor(self.tables, flip=True))
        self.assertFalse(plain["anchor_state"].isin(ANCHORED_STATES).any())
        pd.testing.assert_series_equal(plain["evidence_grade"], flipped["evidence_grade"])
        pd.testing.assert_series_equal(plain["omics_direction"], flipped["omics_direction"])
        self.assertIn("Omics-concordant", set(plain["evidence_grade"]))   # HEX_c：RNA、蛋白、代谢物都提示升

    def test_unidentifiable_or_conflicting_tracers_are_not_anchors(self) -> None:
        mfa = self.tables["mfa"].copy()
        hexokinase = (mfa["contrast_id"] == DIET) & (mfa["feature_id"] == "HEX_c")
        flags = mfa["flags"].fillna("").astype(str)
        mfa.loc[hexokinase, "flags"] = flags[hexokinase] + ";not_identifiable"
        profile = profile_from({**self.tables, "mfa": mfa})
        row = profile.loc[(DIET, "liver", "HEX_c")]
        self.assertTrue(row["gate_passed"])
        self.assertEqual(row["anchor_state"], "not_identifiable")
        self.assertEqual(row["evidence_grade"], "Omics-concordant")

        scorer = EvidenceProfileScorer({"mfa_inference_unit": "biological_replicate"})
        record = {"gate_passed": True, "mfa_evidence_mode": "statistical_test", "tracer_evidence": "multi_tracer_conflict"}
        self.assertEqual(scorer.anchor_state(record), "tracer_conflict")
        self.assertEqual(scorer.anchor_state({**record, "tracer_evidence": "", "n_min_group": 1}), "insufficient")
        self.assertEqual(scorer.anchor_state({**record, "tracer_evidence": "", "n_eff": 3}), "insufficient")
        self.assertEqual(scorer.anchor_state({**record, "tracer_evidence": "", "n_eff": 8}), "established")

    def test_mfa_only_rows_keep_old_grade_meaning(self) -> None:
        """只有 MFA 数据、锚点不成立时仍为 No-signal（旧等级含义不变），不是 Not-assessable。"""
        scorer = EvidenceProfileScorer({"mfa_inference_unit": "biological_replicate"})
        empty = {"evidence_mode": "no_data", "n_sig": 0}
        record = {"flux_id": "HEX_c", "tissue": "muscle", "comparison": DIET, "flux_log2fc": 0.4, "flux_fdr": 0.6,
                  "flux_se": 0.3, "n_eff": 8, "gate_passed": False, "mfa_evidence_mode": "insufficient",
                  "mfa_details": {"evidence_mode": "insufficient", "se": 0.3, "n_eff": 8},
                  "tx_details": empty, "prot_details": empty, "metab_details": empty}
        summary = scorer.summarize(scorer.build_records(record))
        self.assertEqual((summary["anchor_state"], summary["evidence_grade"]), ("not_established", "No-signal"))

    def test_evidence_table_carries_anchor_and_implied_direction(self) -> None:
        _, evidence = build(mfa_frames(self.tables["mfa"]), transcriptomics_frames(self.tables["transcriptomics"]),
                            proteomics_frames(self.tables["proteomics"]), metabolomics_frames(self.tables["metabolomics"]))
        rows = evidence.set_index(["comparison", "tissue", "flux_id", "layer"])
        self.assertEqual(rows.loc[(DIET, "liver", "HEX_c", "transcriptomics"), "implied_direction"], "up")
        self.assertEqual(rows.loc[(DIET, "liver", "G6Pase_c", "metabolomics"), "implied_direction"], "down")
        self.assertTrue((evidence["anchor_state"].notna()).all())
        self.assertTrue(np.isin(evidence["implied_direction"], ["up", "down", "mixed", "none"]).all())


if __name__ == "__main__":
    unittest.main()
