#!/usr/bin/env python3
import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import model_quality


class ModelQualityTests(unittest.TestCase):
    def test_detect_json_shape(self):
        proc = subprocess.run(
            [sys.executable, "scripts/model_quality.py", "detect", "--json"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        report = json.loads(proc.stdout)

        self.assertEqual(report["status"], "ok")
        self.assertIn("data", report)
        self.assertEqual(set(report["tools"]), {"usalign", "dockq", "foldseek"})
        self.assertEqual(set(report["capabilities"]), {"usalign", "dockq", "foldseek"})
        for name, tool in report["tools"].items():
            self.assertIn(name, report["capabilities"])
            self.assertIsInstance(tool["available"], bool)
            self.assertEqual(report["capabilities"][name], tool["available"])
            self.assertIn("path", tool)
            self.assertIn("executable", tool)
            self.assertIsInstance(tool["candidates"], list)
            self.assertGreaterEqual(len(tool["candidates"]), 1)

    def test_parse_usalign_stdout_extracts_metrics(self):
        stdout = """
Name of Chain_1: reference.pdb
Name of Chain_2: mobile.pdb
Length of Chain_1: 100 residues
Length of Chain_2: 98 residues

Aligned length=   98, RMSD=   0.85, Seq_ID=n_identical/n_aligned= 0.990
TM-score= 0.99883 (normalized by length of Chain_1: L=100, d0=4.50)
TM-score= 0.99710 (if normalized by length of Chain_2, i.e., LN=98, d0=4.46)
"""
        metrics = model_quality.parse_usalign_stdout(stdout)

        self.assertEqual(metrics["aligned_length"], 98)
        self.assertEqual(metrics["rmsd"], 0.85)
        self.assertEqual(metrics["seq_id"], 0.99)
        self.assertEqual(metrics["structures"]["structure_1"], {"name": "reference.pdb", "length": 100})
        self.assertEqual(metrics["structures"]["structure_2"], {"name": "mobile.pdb", "length": 98})
        self.assertEqual(len(metrics["tm_scores"]), 2)
        self.assertEqual(metrics["tm_scores"][0]["score"], 0.99883)
        self.assertEqual(metrics["tm_scores"][0]["normalized_by"], "Chain_1")
        self.assertEqual(metrics["tm_scores"][0]["length"], 100)
        self.assertEqual(metrics["tm_scores"][0]["d0"], 4.5)
        self.assertEqual(metrics["tm_scores"][1]["score"], 0.9971)
        self.assertEqual(metrics["tm_scores"][1]["normalized_by"], "Chain_2")
        self.assertEqual(metrics["tm_scores"][1]["length"], 98)
        self.assertEqual(metrics["tm_scores"][1]["d0"], 4.46)


if __name__ == "__main__":
    unittest.main()
