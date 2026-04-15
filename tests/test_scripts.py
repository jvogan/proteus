#!/usr/bin/env python3
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_script(*args):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )


class ScriptTests(unittest.TestCase):
    def test_structure_info_pdb_json(self):
        proc = run_script("scripts/structure_info.py", "tests/fixtures/tiny.pdb", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["format"], "pdb")
        self.assertEqual(data["atom_records"], 3)
        self.assertEqual(data["hetatm_records"], 1)
        self.assertEqual(data["chains"], ["A", "B"])

    def test_structure_info_mmcif_json(self):
        proc = run_script("scripts/structure_info.py", "tests/fixtures/tiny.cif", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["format"], "mmcif")
        self.assertEqual(data["atom_records"], 2)
        self.assertEqual(data["hetatm_records"], 1)
        self.assertEqual(data["title"], "Tiny mmCIF test structure")

    def test_cli_help(self):
        for script in [
            "scripts/fetch_pdb.py",
            "scripts/uniprot_lookup.py",
            "scripts/structure_info.py",
            "scripts/fetch_alphafold.py",
            "scripts/pdb_info.py",
            "scripts/pymol_agent.py",
            "scripts/chimerax_agent.py",
        ]:
            proc = run_script(script, "--help")
            self.assertIn("usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
