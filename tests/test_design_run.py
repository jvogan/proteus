#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import design_run


def run_script(*args, check=True):
    proc = subprocess.run(
        [sys.executable, "scripts/design_run.py", *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if check and proc.returncode != 0:
        raise AssertionError(f"command failed: {args}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    return proc


class DesignRunTests(unittest.TestCase):
    def test_dry_run_writes_plan_commands_candidates_and_scrubs_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            target = tmp_path / "target.pdb"
            ligand = tmp_path / "ligand.sdf"
            target.write_text("ATOM      1  CA  GLY A   1       0.000   0.000   0.000\n", encoding="utf-8")
            ligand.write_text("ligand\n", encoding="utf-8")
            manifest = {
                "name": "local-design",
                "inputs": {
                    "backbone": {"path": str(target), "kind": "structure"},
                    "ligand": {"path": str(ligand), "kind": "ligand"},
                },
                "constraints": {
                    "candidate_count": 3,
                    "fixed_residues": ["A:1"],
                },
                "stages": [
                    {
                        "id": "design",
                        "tool": "proteinmpnn",
                        "inputs": ["backbone"],
                        "params": {"num_seq_per_target": 3, "temperature": 0.2},
                    },
                    {
                        "id": "fold",
                        "tool": "colabfold",
                        "depends_on": ["design"],
                        "params": {"num_models": 1},
                    },
                ],
            }
            manifest_path = tmp_path / "manifest.json"
            outdir = tmp_path / "run"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            proc = run_script(str(manifest_path), "--outdir", str(outdir), "--dry-run", "--json")
            summary = json.loads(proc.stdout)

            self.assertEqual(summary["status"], "ok")
            self.assertTrue(summary["dry_run"])
            self.assertFalse(summary["heavy_tools_executed"])
            self.assertEqual(summary["counts"]["stages"], 2)
            self.assertEqual(summary["counts"]["candidates"], 3)
            self.assertEqual(summary["outputs"]["plan_json"], "<run>/plan.json")
            self.assertTrue((outdir / "plan.json").exists())
            self.assertTrue((outdir / "commands.md").exists())
            self.assertTrue((outdir / "candidates.jsonl").exists())

            plan = json.loads((outdir / "plan.json").read_text(encoding="utf-8"))
            self.assertEqual(plan["mode"], "plan_only")
            self.assertEqual(plan["inputs"][0]["path"], "<abs>/target.pdb")
            self.assertEqual(plan["inputs"][1]["path"], "<abs>/ligand.sdf")
            self.assertEqual(plan["stages"][0]["command"][:4], [
                "protein_mpnn_run.py",
                "--pdb_path",
                "<abs>/target.pdb",
                "--out_folder",
            ])
            self.assertIn("colabfold", plan["tools"])

            candidates = [
                json.loads(line)
                for line in (outdir / "candidates.jsonl").read_text(encoding="utf-8").splitlines()
            ]
            self.assertEqual([item["candidate_id"] for item in candidates], [
                "candidate_0001",
                "candidate_0002",
                "candidate_0003",
            ])
            self.assertTrue(all(item["status"] == "planned" for item in candidates))

            combined_output = "\n".join([
                proc.stdout,
                (outdir / "plan.json").read_text(encoding="utf-8"),
                (outdir / "commands.md").read_text(encoding="utf-8"),
                (outdir / "candidates.jsonl").read_text(encoding="utf-8"),
            ])
            self.assertNotIn(str(tmp_path), combined_output)
            self.assertNotIn(str(ROOT), combined_output)

    def test_detect_tools_by_executable_or_import_without_importing_modules(self):
        def fake_which(executable):
            if executable == "protein_mpnn_run.py":
                return "/opt/proteinmpnn/protein_mpnn_run.py"
            return None

        def fake_find_spec(import_name):
            if import_name == "chai_lab":
                return object()
            return None

        with patch("design_run.shutil.which", side_effect=fake_which), patch(
            "design_run.importlib.util.find_spec", side_effect=fake_find_spec
        ):
            tools = design_run.detect_tools()

        self.assertTrue(tools["proteinmpnn"]["available"])
        self.assertEqual(tools["proteinmpnn"]["detection"], "executable")
        self.assertEqual(tools["proteinmpnn"]["executable"], "protein_mpnn_run.py")
        self.assertEqual(tools["proteinmpnn"]["path"], "<abs>/protein_mpnn_run.py")
        self.assertTrue(tools["chai-lab"]["available"])
        self.assertEqual(tools["chai-lab"]["detection"], "import")
        self.assertEqual(tools["chai-lab"]["import_name"], "chai_lab")
        self.assertFalse(tools["colabfold"]["available"])

    def test_rejects_invalid_manifest_as_json_error(self):
        manifest = {
            "inputs": {"backbone": "target.pdb"},
            "stages": [{"id": "design", "tool": "not-a-tool", "inputs": ["backbone"]}],
        }
        with tempfile.TemporaryDirectory() as tmp:
            manifest_path = Path(tmp) / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            proc = run_script(str(manifest_path), "--json", check=False)

            self.assertEqual(proc.returncode, 2)
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "error")
            self.assertIn("tool must be one of", summary["error"])
            self.assertNotIn(str(manifest_path), proc.stdout)

    def test_help_documents_schema_examples(self):
        proc = run_script("--help")

        self.assertIn("Example manifest:", proc.stdout)
        self.assertIn('"inputs"', proc.stdout)
        self.assertIn('"stages"', proc.stdout)
        self.assertIn('"constraints"', proc.stdout)
        self.assertIn("ProteinMPNN", proc.stdout)


if __name__ == "__main__":
    unittest.main()
