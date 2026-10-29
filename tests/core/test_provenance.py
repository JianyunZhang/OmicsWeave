#!/usr/bin/env python3
"""运行溯源：文件指纹、run_provenance.json 读写、跨批次比对，以及 MFA 的批次可合并性断言。"""

from __future__ import annotations

import contextlib
import hashlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from omicsweave.core import provenance
from omicsweave.mfa import provenance as mfa_provenance


class FingerprintTest(unittest.TestCase):
    def test_full_hash_matches_hashlib(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ref.fa"
            path.write_bytes(b">tx1\nACGT\n" * 1000)
            record = provenance.file_fingerprint(path)
            self.assertEqual(record["hash_mode"], "sha256_full")
            self.assertEqual(record["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual(record["sha256"], provenance.sha256_file(path))

    def test_large_file_uses_head_tail_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(provenance, "FULL_HASH_MAX_BYTES", 10):
            path = Path(tmp) / "index.idx"
            path.write_bytes(b"x" * 100)
            record = provenance.file_fingerprint(path)
        self.assertEqual(record["hash_mode"], "sha256_head_tail")
        self.assertEqual(record["size_bytes"], 100)


class RecordTest(unittest.TestCase):
    def test_write_load_round_trip_and_bad_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "runs" / "b1"
            path = provenance.write_run_provenance(run, {"strandedness": "rf-stranded", "组织": "肝脏"})
            self.assertEqual(path.name, provenance.PROVENANCE_FILE)
            self.assertEqual(provenance.load_run_provenance(run), {"strandedness": "rf-stranded", "组织": "肝脏"})
            self.assertIsNone(provenance.load_run_provenance(Path(tmp) / "missing"))
            path.write_text("{not json", encoding="utf-8")
            self.assertIsNone(provenance.load_run_provenance(run))


class CompareTest(unittest.TestCase):
    def test_advisory_mismatch_keeps_batches_consistent(self) -> None:
        result = provenance.compare_provenance(
            {"b1": {"tracer": "glc", "device": "QE"}, "b2": {"tracer": "glc", "device": "Exploris"}},
            critical=["tracer"], advisory=["device"],
        )
        self.assertTrue(result["consistent"])
        self.assertEqual(result["severity"], {"device": "advisory"})
        self.assertIn("（提示）", provenance.format_provenance_report(result))

    def test_critical_mismatch_and_missing_batch(self) -> None:
        result = provenance.compare_provenance(
            {"b1": {"tracer": "glc"}, "b2": {"tracer": "gln"}, "b3": None}, critical=["tracer"],
        )
        self.assertFalse(result["consistent"])
        self.assertEqual(result["missing"], ["b3"])
        self.assertEqual(result["n_compared"], 2)
        report = provenance.format_provenance_report(result)
        self.assertIn("b3", report)
        self.assertIn("字段 tracer 不一致:", report)


class MfaBatchCompatibilityTest(unittest.TestCase):
    """MFA 在 runs/<protocol>/<batch>/ 下写溯源，整合前按 fail / warn / ignore 策略比对。"""

    def _write(self, root: Path, batch: str, tracer_type: str, device: str) -> None:
        with contextlib.redirect_stdout(io.StringIO()):
            mfa_provenance.write_run_provenance(
                root / "runs" / "proto" / batch, experiment_id="e", protocol_id="proto", batch_id=batch,
                batch_metadata={"detection_device": device, "metabolites_extraction_protocal": "MeOH"},
                tracer={"type": tracer_type, "isotope": "13C"}, reactions=["R1", "R2"], tissues=["liver"],
                conditions=["a", "b"], n_samples=6, n_fit_iterations=150, analysis_config={},
            )

    def _check(self, root: Path, policy: str):
        with contextlib.redirect_stdout(io.StringIO()):
            return mfa_provenance.assert_batch_compatibility(root, "proto", ["b1", "b2"], policy=policy)

    def test_policies(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root, "b1", "glucose", "QE")
            self._write(root, "b2", "glucose", "Exploris")
            mismatches = self._check(root, "fail")  # 只有提示字段不同：不报错
            self.assertEqual([(m["field"], m["severity"]) for m in mismatches], [("detection_device", "advisory")])

            self._write(root, "b2", "glutamine", "QE")
            self.assertEqual(self._check(root, "warn")[0]["field"], "tracer_type")
            self.assertEqual(self._check(root, "ignore"), [])
            with self.assertRaisesRegex(RuntimeError, "tracer_type"):
                self._check(root, "fail")


if __name__ == "__main__":
    unittest.main()
