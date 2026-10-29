#!/usr/bin/env python3
"""数据契约：分析层数据集、差异结果表与运行清单的读写和校验。"""

from __future__ import annotations

import json
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import ValidationError

from omicsweave.contract import (
    ContractError,
    Dataset,
    DatasetMeta,
    RunManifest,
    export_json_schema,
    file_record,
    new_run_id,
    read_dataset,
    read_differential,
    read_manifest,
    write_dataset,
    write_differential,
    write_manifest,
)


def counts_dataset() -> Dataset:
    index = pd.Index(["g1", "g2", "g3"], name="feature_id")
    matrix = pd.DataFrame(np.arange(6).reshape(3, 2), index=index, columns=["s1", "s2"])
    features = pd.DataFrame({"symbol": ["A", "B", "C"]}, index=index)
    meta = DatasetMeta(omics="transcriptomics", tissue="liver", value_type="counts", scale="linear",
                       feature_id_type="ensembl_gene")
    return Dataset(matrix=matrix, features=features, meta=meta)


def differential_table() -> pd.DataFrame:
    return pd.DataFrame({
        "contrast_id": ["HFD_vs_LFD"] * 2,
        "contrast_role": ["disease_effect"] * 2,
        "tissue": ["liver"] * 2,
        "omics": ["transcriptomics"] * 2,
        "feature_id": ["g1", "g2"],
        "log2fc": [1.2, -0.4],
        "pvalue": [0.001, 0.3],
        "padj": [0.002, 0.4],
        "padj_method": ["fdr_bh"] * 2,
        "n_numerator": [4, 4],
        "n_denominator": [4, 4],
        "method": ["deseq2"] * 2,
    })


class DatasetTest(unittest.TestCase):
    def test_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_dataset(tmp, counts_dataset())
            loaded = read_dataset(tmp)
        pd.testing.assert_frame_equal(loaded.matrix, counts_dataset().matrix, check_names=False)
        self.assertEqual(loaded.meta.value_type.value, "counts")

    def test_feature_mismatch_is_rejected(self) -> None:
        dataset = counts_dataset()
        dataset.features = dataset.features.iloc[:2]
        with self.assertRaises(ContractError):
            dataset.check()

    def test_negative_counts_are_rejected(self) -> None:
        dataset = counts_dataset()
        dataset.matrix.iloc[0, 0] = -1
        with self.assertRaisesRegex(ContractError, "负值"):
            dataset.check()

    def test_flux_requires_units(self) -> None:
        with self.assertRaises(ValidationError):
            DatasetMeta(omics="mfa", tissue="liver", value_type="flux", scale="linear",
                        feature_id_type="mfa_reaction")


class DifferentialTest(unittest.TestCase):
    def test_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = write_differential(differential_table(), Path(tmp) / "differential.parquet")
            loaded = read_differential(path)
        self.assertEqual(list(loaded["feature_id"]), ["g1", "g2"])

    def test_missing_required_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp, self.assertRaisesRegex(ContractError, "n_numerator"):
            write_differential(differential_table().drop(columns="n_numerator"), Path(tmp) / "d.parquet")

    def test_out_of_range_and_duplicates(self) -> None:
        table = differential_table()
        table.loc[0, "padj"] = 1.5
        table.loc[1, "feature_id"] = "g1"
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(ContractError) as ctx:
            write_differential(table, Path(tmp) / "d.parquet")
        self.assertIn("padj", str(ctx.exception))
        self.assertIn("重复", str(ctx.exception))


class ManifestTest(unittest.TestCase):
    def test_run_id_format(self) -> None:
        self.assertRegex(new_run_id(), r"^run_\d{8}_\d{6}_[0-9a-f]{4}$")

    def test_roundtrip_with_file_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root / "input.tsv"
            data.write_text("a\tb\n1\t2\n", encoding="utf-8")
            record = file_record(data, root, "input", kind="counts_matrix")
            manifest = RunManifest(run_id=new_run_id(), project_id="p", experiment_id="e",
                                   omics=["transcriptomics"], files=[record])
            write_manifest(manifest, root / "run")
            loaded = read_manifest(root / "run")
        self.assertEqual(loaded.files[0].path, "input.tsv")
        self.assertTrue(re.fullmatch(r"[0-9a-f]{64}", loaded.files[0].sha256))

    def test_json_schema_is_serializable(self) -> None:
        schema = export_json_schema()
        self.assertEqual(set(schema), {"dataset", "manifest"})
        json.dumps(schema)


if __name__ == "__main__":
    unittest.main()
