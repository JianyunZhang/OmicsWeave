#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回归测试：样本间 MIX 列不一致时，真实通量仍应正常加载。"""

import tempfile
import unittest
from pathlib import Path

import pandas as pd



from omicsweave.mfa.mfa_analysis import MFAAnalyser


class FluxLoaderVariableMixColumnsTest(unittest.TestCase):
    def test_variable_mix_columns_do_not_break_flux_loading(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            temp_path = Path(temp_dir)
            workbook = temp_path / "flux_raw_data.xlsx"

            sheet_info = pd.DataFrame(
                {
                    "complete_result_name": ["sample_1", "sample_2"],
                    "sheet_name": ["sample_1_fit", "sample_2_fit"],
                    "condition": ["control", "treatment"],
                }
            )
            sample_1 = pd.DataFrame(
                {
                    "loss": [2.0, 1.0],
                    "liver--R1": [10.0, 20.0],
                    "MIX:liver--citrate_0": [0.9, 0.8],
                }
            )
            sample_2 = pd.DataFrame(
                {
                    "loss": [3.0, 1.0],
                    "liver--R1": [30.0, 40.0],
                }
            )

            with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
                sheet_info.to_excel(writer, sheet_name="sheet_information", index=False)
                sample_1.to_excel(writer, sheet_name="sample_1_fit", index=False)
                sample_2.to_excel(writer, sheet_name="sample_2_fit", index=False)

            analyser = MFAAnalyser(str(temp_path / "output"))
            analyser.load_flux_data(str(workbook))

            self.assertEqual(["liver--R1"], analyser.flux_df.columns.tolist())
            self.assertEqual(["sample_1", "sample_2"], analyser.flux_df.index.tolist())
            self.assertEqual(20.0, analyser.flux_df.loc["sample_1", "liver--R1"])
            self.assertEqual(40.0, analyser.flux_df.loc["sample_2", "liver--R1"])


if __name__ == "__main__":
    unittest.main()
