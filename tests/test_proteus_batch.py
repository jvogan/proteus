#!/usr/bin/env python3
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


class ProteusBatchTests(unittest.TestCase):
    def test_dry_run_writes_planned_results_for_allowlisted_tasks(self):
        manifest = {
            "items": [{
                "id": "tiny",
                "input": "tests/fixtures/tiny.pdb",
                "tasks": [
                    "structure_info",
                    "pocket_report",
                    "interface_report",
                    "validation_report",
                    "resolve_structure",
                    "ligand_extract",
                    "interaction_report",
                    {"name": "docking_box", "args": ["--ligand", "LIG"]},
                    {"name": "mutation_triage", "args": ["G1A"]},
                ],
            }]
        }
        with tempfile.TemporaryDirectory() as tmp:
            manifest_path = Path(tmp) / "manifest.json"
            outdir = Path(tmp) / "out"
            manifest_path.write_text(json.dumps(manifest))

            proc = run_script(
                "scripts/proteus_batch.py",
                str(manifest_path),
                "--outdir",
                str(outdir),
                "--dry-run",
                "--json",
            )

            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "ok")
            self.assertTrue(summary["dry_run"])
            self.assertEqual(summary["counts"]["tasks_planned"], 9)
            self.assertEqual(summary["items"][0]["status"], "planned")

            lines = (outdir / "results.jsonl").read_text().splitlines()
            self.assertEqual(len(lines), 9)
            results = [json.loads(line) for line in lines]
            self.assertTrue(all(result["status"] == "ok" for result in results))
            self.assertTrue(all(result["data"]["task_status"] == "planned" for result in results))
            self.assertEqual(results[0]["command"][:2], ["python3", "scripts/structure_info.py"])
            self.assertEqual(results[-2]["task"], "docking_box")
            self.assertEqual(results[-2]["command"][-3:], ["--ligand", "LIG", "--json"])
            self.assertEqual(results[-1]["task"], "mutation_triage")
            self.assertIn("G1A", results[-1]["command"])
            self.assertEqual(
                results[-1]["command"][-4:],
                ["G1A", "--structure", "tests/fixtures/tiny.pdb", "--json"],
            )
            self.assertNotIn(str(ROOT), (outdir / "summary.md").read_text())

    def test_offline_run_continues_after_item_error(self):
        manifest = {
            "items": [
                {
                    "id": "missing",
                    "input": "tests/fixtures/does_not_exist.pdb",
                    "tasks": ["structure_info"],
                },
                {
                    "id": "tiny",
                    "input": "tests/fixtures/tiny.pdb",
                    "tasks": ["structure_info", "resolve_structure"],
                },
            ]
        }
        with tempfile.TemporaryDirectory() as tmp:
            manifest_path = Path(tmp) / "manifest.json"
            outdir = Path(tmp) / "out"
            manifest_path.write_text(json.dumps(manifest))

            proc = run_script(
                "scripts/proteus_batch.py",
                str(manifest_path),
                "--outdir",
                str(outdir),
                "--json",
                check=False,
            )

            self.assertEqual(proc.returncode, 1)
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "error")
            self.assertEqual(summary["counts"]["tasks_error"], 1)
            self.assertEqual(summary["counts"]["tasks_ok"], 2)
            self.assertEqual(summary["items"][0]["status"], "error")
            self.assertEqual(summary["items"][1]["status"], "ok")

            results = [
                json.loads(line)
                for line in (outdir / "results.jsonl").read_text().splitlines()
            ]
            self.assertEqual([result["status"] for result in results], ["error", "ok", "ok"])
            self.assertIn("File not found", results[0]["error"])
            self.assertEqual(results[2]["result"]["resolved_kind"], "local_file")

    def test_rejects_invalid_task_object_args(self):
        manifest = {
            "items": [{
                "id": "bad",
                "input": "tests/fixtures/tiny.pdb",
                "tasks": [{"name": "structure_info", "args": "--not-a-list"}],
            }]
        }
        with tempfile.TemporaryDirectory() as tmp:
            manifest_path = Path(tmp) / "manifest.json"
            manifest_path.write_text(json.dumps(manifest))

            proc = run_script(
                "scripts/proteus_batch.py",
                str(manifest_path),
                "--json",
                check=False,
            )

            self.assertEqual(proc.returncode, 2)
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "error")
            self.assertIn("args must be a list", summary["error"])


if __name__ == "__main__":
    unittest.main()
