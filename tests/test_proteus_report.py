#!/usr/bin/env python3
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_script(*args, check=True):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ProteusReportTests(unittest.TestCase):
    def test_cli_combines_enveloped_and_plain_json_with_scrubbed_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            private_structure = tmpdir / "private" / "model.pdb"
            private_structure.parent.mkdir()
            wrapped = tmpdir / "wrapped.json"
            plain = tmpdir / "plain.json"
            outdir = tmpdir / "out"

            wrapped.write_text(json.dumps({
                "status": "ok",
                "data": {
                    "source": str(private_structure),
                    "pdb_id": "1ABC",
                    "chains": ["A", "B"],
                    "warnings": [f"Review local file {private_structure}"],
                },
            }))
            plain.write_text(json.dumps({
                "tool": "plain-helper",
                "value": 42,
                "path": str(private_structure),
            }))

            proc = run_script(
                "scripts/proteus_report.py",
                "--input",
                str(wrapped),
                "--input",
                str(plain),
                "--outdir",
                str(outdir),
                "--title",
                "Local Evidence",
                "--json",
            )

            stdout_report = json.loads(proc.stdout)
            disk_report = json.loads((outdir / "report.json").read_text())
            self.assertEqual(stdout_report, disk_report)
            self.assertEqual(disk_report["status"], "ok")
            self.assertEqual(disk_report["title"], "Local Evidence")
            self.assertEqual(disk_report["counts"]["inputs_total"], 2)
            self.assertEqual(disk_report["counts"]["inputs_ok"], 2)
            self.assertEqual(disk_report["reports"][0]["shape"], "proteus_envelope")
            self.assertEqual(disk_report["reports"][1]["shape"], "plain_json")
            self.assertEqual(disk_report["reports"][0]["provenance"]["sha256"], sha256(wrapped))
            self.assertEqual(disk_report["reports"][1]["provenance"]["sha256"], sha256(plain))
            self.assertTrue((outdir / "evidence" / "01-wrapped.json").exists())
            self.assertTrue((outdir / "evidence" / "02-plain.json").exists())

            markdown = (outdir / "REPORT.md").read_text()
            report_text = json.dumps(disk_report, sort_keys=True)
            self.assertNotIn(str(tmpdir), markdown)
            self.assertNotIn(str(tmpdir), report_text)
            self.assertIn("wrapped.json (absolute path omitted)", markdown)
            self.assertIn("model.pdb (absolute path omitted)", report_text)

    def test_no_copy_omits_evidence_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            report_path = tmpdir / "helper.json"
            outdir = tmpdir / "out"
            report_path.write_text(json.dumps({"status": "ok", "data": {"format": "pdb"}}))

            proc = run_script(
                "scripts/proteus_report.py",
                "--input",
                str(report_path),
                "--outdir",
                str(outdir),
                "--no-copy",
                "--json",
            )

            combined = json.loads(proc.stdout)
            self.assertFalse((outdir / "evidence").exists())
            self.assertFalse(combined["copy_inputs"])
            self.assertNotIn("copied_to", combined["reports"][0]["provenance"])
            self.assertIn("Input copying disabled", (outdir / "REPORT.md").read_text())

    def test_coordinate_like_contents_are_pruned_from_combined_outputs(self):
        atom_line = (
            "ATOM      1  CA  GLY A   1      11.000  12.000  13.000  "
            "1.00 20.00           C"
        )
        raw_pdb = "\n".join(atom_line for _ in range(80))
        coordinate_rows = [{"atom_name": "CA", "residue_name": "GLY", "chain": "A", "x": i, "y": i, "z": i}
                           for i in range(150)]

        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_path = tmpdir / "raw.json"
            outdir = tmpdir / "out"
            input_path.write_text(json.dumps({
                "status": "ok",
                "data": {
                    "raw_pdb": raw_pdb,
                    "atoms": coordinate_rows,
                    "path": str(tmpdir / "secret" / "structure.pdb"),
                },
            }))

            proc = run_script(
                "scripts/proteus_report.py",
                "--input",
                str(input_path),
                "--outdir",
                str(outdir),
                "--json",
            )

            combined = json.loads(proc.stdout)
            report_data = combined["reports"][0]["data"]
            self.assertTrue(report_data["raw_pdb"]["omitted"])
            self.assertTrue(report_data["atoms"]["omitted"])
            self.assertGreaterEqual(combined["counts"]["warnings"], 3)
            self.assertEqual(combined["counts"]["inputs_copied"], 0)
            self.assertEqual(combined["counts"]["inputs_copy_skipped"], 1)
            self.assertEqual(combined["reports"][0]["provenance"]["copy_skipped"], "pruned_large_content")
            self.assertFalse((outdir / "evidence" / "01-raw.json").exists())

            report_text = (outdir / "report.json").read_text()
            markdown = (outdir / "REPORT.md").read_text()
            self.assertNotIn("ATOM      1", report_text)
            self.assertNotIn("ATOM      1", markdown)
            self.assertNotIn(str(tmpdir), report_text)
            self.assertNotIn(str(tmpdir), markdown)


if __name__ == "__main__":
    unittest.main()
