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


class TargetDossierTests(unittest.TestCase):
    def test_cli_help_lists_core_options(self):
        proc = run_script("scripts/target_dossier.py", "--help")
        self.assertIn("--gene", proc.stdout)
        self.assertIn("--uniprot", proc.stdout)
        self.assertIn("--pdb", proc.stdout)
        self.assertIn("--no-network", proc.stdout)
        self.assertIn("--analyze-local", proc.stdout)
        self.assertIn("--json", proc.stdout)

    def test_no_network_mixed_inputs_write_markdown_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp) / "out"
            proc = run_script(
                "scripts/target_dossier.py",
                "--gene",
                "KRAS",
                "--uniprot",
                "P01116",
                "--pdb",
                "4HHB",
                "--pdb",
                "tests/fixtures/tiny.pdb",
                "--out",
                str(outdir),
                "--no-network",
                "--json",
            )

            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "ok")
            self.assertFalse(summary["network"]["enabled"])
            self.assertEqual(summary["network"]["mode"], "offline")
            self.assertEqual(summary["outputs"]["markdown"], "TARGET_DOSSIER.md")
            self.assertEqual(summary["outputs"]["provenance"], "provenance.json")

            markdown_path = outdir / "TARGET_DOSSIER.md"
            provenance_path = outdir / "provenance.json"
            self.assertTrue(markdown_path.exists())
            self.assertTrue(provenance_path.exists())

            markdown = markdown_path.read_text()
            provenance_text = provenance_path.read_text()
            provenance = json.loads(provenance_text)

            self.assertIn("Network access disabled", markdown)
            self.assertIn("KRAS", markdown)
            self.assertIn("P01116", markdown)
            self.assertIn("4HHB", markdown)
            self.assertIn("tests/fixtures/tiny.pdb", markdown)
            self.assertNotIn(str(ROOT), markdown)
            self.assertNotIn(str(ROOT), provenance_text)

            commands = provenance["commands"]
            self.assertEqual(len(commands), 1)
            self.assertEqual(commands[0]["script"], "scripts/structure_info.py")
            self.assertIn("tests/fixtures/tiny.pdb", commands[0]["command"])

            statuses = {(item["kind"], item["value"]): item["status"] for item in provenance["inputs"]}
            self.assertEqual(statuses[("gene", "KRAS")], "skipped")
            self.assertEqual(statuses[("uniprot", "P01116")], "skipped")
            self.assertEqual(statuses[("pdb", "4HHB")], "skipped")
            self.assertEqual(statuses[("pdb", "tests/fixtures/tiny.pdb")], "ok")

            local = [item for item in provenance["structures"] if item["kind"] == "local_structure"][0]
            self.assertEqual(local["status"], "ok")
            self.assertEqual(local["format"], "pdb")
            self.assertEqual(local["atom_records"], 3)
            self.assertEqual(local["inspection"]["file"], "./tests/fixtures/tiny.pdb")
            self.assertNotIn("local_analyses", local)

    def test_analyze_local_adds_structural_analysis_summaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp) / "out"
            proc = run_script(
                "scripts/target_dossier.py",
                "--pdb",
                "tests/fixtures/tiny.pdb",
                "--out",
                str(outdir),
                "--no-network",
                "--analyze-local",
                "--json",
            )

            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "ok")
            structure = summary["structures"][0]
            analyses = structure["local_analyses"]
            self.assertEqual(analyses["ligands"]["status"], "ok")
            self.assertEqual(analyses["ligands"]["ligand_group_count"], 1)
            self.assertEqual(analyses["interfaces"]["status"], "ok")
            self.assertEqual(analyses["interactions"]["status"], "ok")
            self.assertGreaterEqual(analyses["interactions"]["contact_count"], 1)
            self.assertEqual(analyses["docking_box"]["status"], "ok")
            self.assertEqual(analyses["pocket"]["status"], "ok")

            markdown = (outdir / "TARGET_DOSSIER.md").read_text()
            provenance_text = (outdir / "provenance.json").read_text()
            self.assertIn("## Local Analyses", markdown)
            self.assertIn("ligands", markdown)
            self.assertNotIn(str(ROOT), markdown)
            self.assertNotIn(str(ROOT), provenance_text)

    def test_absolute_local_paths_are_scrubbed_from_outputs(self):
        pdb_text = (
            "TITLE     temporary test structure\n"
            "ATOM      1  CA  GLY A   1      12.560  13.205   2.100  1.00 20.00           C\n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            pdb_path = tmp_path / "absolute_sample.pdb"
            outdir = tmp_path / "out"
            pdb_path.write_text(pdb_text)

            proc = run_script(
                "scripts/target_dossier.py",
                "--pdb",
                str(pdb_path),
                "--out",
                str(outdir),
                "--no-network",
                "--json",
            )
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "ok")

            markdown = (outdir / "TARGET_DOSSIER.md").read_text()
            provenance_text = (outdir / "provenance.json").read_text()

            self.assertIn("absolute_sample.pdb", markdown)
            self.assertIn("absolute_sample.pdb", provenance_text)
            self.assertNotIn(str(pdb_path), markdown)
            self.assertNotIn(str(pdb_path), provenance_text)
            self.assertNotIn(str(outdir), provenance_text)

            provenance = json.loads(provenance_text)
            local = provenance["structures"][0]
            self.assertEqual(local["id"], "absolute_sample.pdb (absolute path omitted)")
            self.assertEqual(local["inspection"]["file"], "absolute_sample.pdb (absolute path omitted)")

    def test_invalid_pdb_input_returns_error_but_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp) / "out"
            proc = run_script(
                "scripts/target_dossier.py",
                "--pdb",
                "not_a_structure.pdb",
                "--out",
                str(outdir),
                "--no-network",
                "--json",
                check=False,
            )

            self.assertEqual(proc.returncode, 1)
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["status"], "error")
            self.assertIn("invalid or missing input", summary["errors"][0])
            self.assertTrue((outdir / "TARGET_DOSSIER.md").exists())
            self.assertTrue((outdir / "provenance.json").exists())


if __name__ == "__main__":
    unittest.main()
