#!/usr/bin/env python3
"""RNA-Seq 批次设计与元数据契约的回归测试。"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from omicsweave.transcriptomics import rnaseq_pipeline
from omicsweave.transcriptomics.experiment_config import _load_groups
from omicsweave.transcriptomics.run_wilcoxon_de import _validate_no_unadjusted_multibatch


class TranscriptomicsBatchDesignTest(unittest.TestCase):
    def test_groups_load_explicit_sample_batches(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            groups_path = Path(temp_dir) / "groups.yaml"
            groups_path.write_text(yaml.safe_dump({
                "conditions": {"case": {}, "control": {}},
                "sample_groups": {"m1": "case", "m2": "control"},
                "sample_batches": {"m1": "b1", "m2": "b2"},
            }), encoding="utf-8")

            _, groups, batches = _load_groups(groups_path)

            self.assertEqual(groups, {"m1": "case", "m2": "control"})
            self.assertEqual(batches, {"m1": "b1", "m2": "b2"})

    def test_sample_metadata_requires_complete_batch_mapping(self) -> None:
        old_groups = rnaseq_pipeline.MOUSE_GROUPS
        old_batches = rnaseq_pipeline.SAMPLE_BATCHES
        try:
            rnaseq_pipeline.MOUSE_GROUPS = {"m1": "case", "m2": "control"}
            rnaseq_pipeline.SAMPLE_BATCHES = {"m1": "b1"}
            samples = pd.DataFrame({
                "sample": ["m1_liver", "m2_liver"],
                "file1": ["a", "b"],
                "file2": ["c", "d"],
            })

            with self.assertRaisesRegex(ValueError, "missing samples"):
                rnaseq_pipeline.get_sample_metadata(samples)
        finally:
            rnaseq_pipeline.MOUSE_GROUPS = old_groups
            rnaseq_pipeline.SAMPLE_BATCHES = old_batches

    def test_balanced_batch_condition_design_is_identifiable(self) -> None:
        metadata = pd.DataFrame({
            "condition": ["case", "control", "case", "control"],
            "batch": ["b1", "b1", "b2", "b2"],
        }, index=["s1", "s2", "s3", "s4"])

        terms = rnaseq_pipeline._validate_de_design(metadata, "batch + condition")

        self.assertEqual(terms, ["batch", "condition"])

    def test_confounded_batch_condition_design_is_rejected(self) -> None:
        metadata = pd.DataFrame({
            "condition": ["case", "case", "control", "control"],
            "batch": ["b1", "b1", "b2", "b2"],
        }, index=["s1", "s2", "s3", "s4"])

        with self.assertRaisesRegex(ValueError, "rank deficient/confounded"):
            rnaseq_pipeline._validate_de_design(metadata, "batch + condition")
        with self.assertRaisesRegex(ValueError, "Unadjusted Wilcoxon"):
            _validate_no_unadjusted_multibatch(metadata)

    def test_single_batch_design_keeps_unadjusted_fallback_available(self) -> None:
        metadata = pd.DataFrame({
            "condition": ["case", "case", "control", "control"],
            "batch": ["b1", "b1", "b1", "b1"],
        })

        _validate_no_unadjusted_multibatch(metadata)


if __name__ == "__main__":
    unittest.main()
