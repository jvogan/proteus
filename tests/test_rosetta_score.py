#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import rosetta_score


def run_script(*args, check=True):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


class RosettaScoreTests(unittest.TestCase):
    def test_detect_json_shape_without_running_rosetta(self):
        proc = run_script("scripts/rosetta_score.py", "detect", "--json")
        report = json.loads(proc.stdout)

        self.assertEqual(report["status"], "ok")
        self.assertIn("data", report)
        self.assertEqual(set(report["tools"]), {"rosetta_scripts", "score_jd2"})
        self.assertEqual(set(report["capabilities"]), {"rosetta_scripts", "score_jd2", "pyrosetta"})
        for name, tool in report["tools"].items():
            self.assertIn(name, report["capabilities"])
            self.assertIsInstance(tool["available"], bool)
            self.assertEqual(tool["available"], report["capabilities"][name])
            self.assertIn("path", tool)
            self.assertIn("executable", tool)
            self.assertIsInstance(tool["candidates"], list)
            self.assertGreaterEqual(len(tool["candidates"]), 1)
        self.assertIn("checked_by", report["pyrosetta"])
        self.assertEqual(report["pyrosetta"]["checked_by"], "importlib.util.find_spec")

    def test_parse_scorefile_ranks_total_score_ascending(self):
        score_text = """SEQUENCE: AAA
SCORE: total_score fa_atr fa_rep description
SCORE: -10.5 -1.0 0.2 model_0002
SCORE: -8.0 -0.8 0.1 model_0003
SCORE: -12.25 -1.2 0.3 model_0001
"""
        with tempfile.NamedTemporaryFile("w", suffix=".sc", delete=False) as handle:
            handle.write(score_text)
            score_path = Path(handle.name)
        try:
            report = rosetta_score.parse_scorefile(score_path)
        finally:
            score_path.unlink()

        self.assertEqual(report["status"], "ok")
        self.assertEqual(report["row_count"], 3)
        self.assertEqual(report["ranked_count"], 3)
        self.assertEqual([item["description"] for item in report["ranks"]], ["model_0001", "model_0002", "model_0003"])
        self.assertEqual(report["ranks"][0]["rank"], 1)
        self.assertEqual(report["ranks"][0]["total_score"], -12.25)

    def test_parse_scorefile_cli_emits_json_and_scrubs_absolute_path(self):
        score_text = """SCORE: total_score other description
SCORE: 3.0 1 decoy high
SCORE: 1.0 2 decoy low
"""
        with tempfile.TemporaryDirectory() as tmp:
            score_path = Path(tmp) / "score.sc"
            score_path.write_text(score_text)
            proc = run_script(
                "scripts/rosetta_score.py",
                "parse-scorefile",
                str(score_path),
                "--limit",
                "1",
                "--json",
            )

            report = json.loads(proc.stdout)
            rendered = json.dumps(report)

            self.assertEqual(report["status"], "ok")
            self.assertEqual(report["ranks"][0]["description"], "decoy low")
            self.assertIn("score.sc", report["file"])
            self.assertNotIn(str(score_path), rendered)
            self.assertNotIn(str(Path(tmp)), rendered)

    def test_plan_is_dry_run_and_scrubs_external_absolute_structure_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            structure_path = Path(tmp) / "model.pdb"
            structure_path.write_text("ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n")
            report = rosetta_score.plan_scoring(
                [str(structure_path)],
                backend="score_jd2",
                scorefile=str(Path(tmp) / "score.sc"),
            )

            rendered = json.dumps(report)

            self.assertEqual(report["status"], "ok")
            self.assertTrue(report["dry_run"])
            self.assertFalse(report["execute"])
            self.assertEqual(report["selected_backend"], "score_jd2")
            self.assertEqual(report["commands"][0]["backend"], "score_jd2")
            self.assertIn("-in:file:s", report["commands"][0]["command"])
            self.assertIn("model.pdb (absolute path omitted)", rendered)
            self.assertIn("score.sc (absolute path omitted)", rendered)
            self.assertNotIn(str(structure_path), rendered)
            self.assertNotIn(str(Path(tmp)), rendered)
            self.assertIn("No Rosetta or PyRosetta scoring command was executed.", report["warnings"])

    def test_plan_cli_with_protocol_builds_rosetta_scripts_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            structure_path = tmp_path / "model.cif"
            protocol_path = tmp_path / "score.xml"
            structure_path.write_text("data_model\n")
            protocol_path.write_text("<ROSETTASCRIPTS />\n")

            proc = run_script(
                "scripts/rosetta_score.py",
                "plan",
                str(structure_path),
                "--backend",
                "rosetta_scripts",
                "--protocol",
                str(protocol_path),
                "--json",
            )
            report = json.loads(proc.stdout)
            rendered = json.dumps(report)

            self.assertEqual(report["status"], "ok")
            command = report["commands"][0]["command"]
            self.assertEqual(report["selected_backend"], "rosetta_scripts")
            self.assertIn("-parser:protocol", command)
            self.assertIn("score.xml (absolute path omitted)", rendered)
            self.assertIn("model.cif (absolute path omitted)", rendered)
            self.assertNotIn(str(tmp_path), rendered)


if __name__ == "__main__":
    unittest.main()
