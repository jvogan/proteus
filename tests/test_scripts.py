#!/usr/bin/env python3
import ast
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import chimerax_agent
import fetch_alphafold
import fetch_pdb
import pymol_agent
import uniprot_lookup


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
        self.assertEqual(data["status"], "ok")
        self.assertIn("data", data)
        self.assertEqual(data["format"], "pdb")
        self.assertEqual(data["atom_records"], 3)
        self.assertEqual(data["hetatm_records"], 1)
        self.assertEqual(data["chains"], ["A", "B"])

    def test_structure_info_mmcif_json(self):
        proc = run_script("scripts/structure_info.py", "tests/fixtures/tiny.cif", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertIn("data", data)
        self.assertEqual(data["format"], "mmcif")
        self.assertEqual(data["atom_records"], 2)
        self.assertEqual(data["hetatm_records"], 1)
        self.assertEqual(data["title"], "Tiny mmCIF test structure")

    def test_pdb_info_json_contract(self):
        proc = run_script("scripts/pdb_info.py", "tests/fixtures/tiny.pdb", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertIn("data", data)
        self.assertEqual(data["atom_records"], 3)

    def test_pae_report_json(self):
        proc = run_script("scripts/pae_report.py", "tests/fixtures/tiny_pae.json", "--min-segment", "1", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["size"], {"rows": 4, "columns": 4})
        self.assertIn("per_residue_mean_pae", data["data"])

    def test_resolve_local_file_json(self):
        proc = run_script("scripts/resolve_structure.py", "tests/fixtures/tiny.pdb", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["resolved_kind"], "local_file")
        self.assertEqual(data["inspection"]["format"], "pdb")

    def test_pocket_report_local_json(self):
        proc = run_script("scripts/pocket_report.py", "tests/fixtures/tiny.pdb", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["ligand_count"], 1)
        self.assertGreaterEqual(data["ligands"][0]["contact_residue_count"], 1)

    def test_chimerax_parser_keeps_continuation_lines(self):
        stdout = """INFO:\nExecuting: version\nINFO:\nUCSF ChimeraX version: 1.11.1\nINFO:\nmodel id #1 type AtomicStructure name tiny.pdb\n"""
        info, errors = chimerax_agent._parse_output(stdout)
        self.assertFalse(errors)
        self.assertNotIn("Executing: version", info)
        self.assertIn("UCSF ChimeraX version: 1.11.1", info)
        self.assertIn("model id #1 type AtomicStructure name tiny.pdb", info)

    def test_chimerax_quotes_paths_and_rejects_unsafe_chains(self):
        quoted = chimerax_agent._quote_chimerax_value('path with spaces/"quote";still-path')
        self.assertTrue(quoted.startswith('"') and quoted.endswith('"'))
        self.assertIn('\\"quote\\"', quoted)
        self.assertIn(";still-path", quoted)
        with self.assertRaises(ValueError):
            chimerax_agent._chain_spec("#1", "A; close session")
        result = chimerax_agent.find_hbonds("tests/fixtures/tiny.pdb", chain1="A; close session")
        self.assertEqual(result["status"], "error")

    def test_pymol_literals_and_color_validation(self):
        value = 'path with spaces/"quote";still-path'
        self.assertEqual(ast.literal_eval(pymol_agent._py_literal(value)), value)
        self.assertEqual(pymol_agent._validate_pymol_color("carbon"), "carbon")
        with self.assertRaises(ValueError):
            pymol_agent._validate_pymol_color("red; import os")
        result = pymol_agent.render_structure("tests/fixtures/tiny.pdb", "out.png", color="red; import os")
        self.assertEqual(result["status"], "error")

    def test_snapshot_json_files_parse(self):
        for path in (ROOT / "docs" / "snapshots").glob("*.json"):
            data = json.loads(path.read_text())
            self.assertIn(data["status"], {"ok", "error"})

    def test_error_payloads_are_json_parseable(self):
        for payload in [
            fetch_alphafold._error_payload("bad alphafold input"),
            fetch_pdb._error_payload("bad pdb input"),
            uniprot_lookup._error_payload("bad uniprot input"),
        ]:
            encoded = json.dumps(payload)
            decoded = json.loads(encoded)
            self.assertEqual(decoded["status"], "error")
            self.assertIn("error", decoded)

    def test_cli_help(self):
        for script in [
            "scripts/fetch_pdb.py",
            "scripts/uniprot_lookup.py",
            "scripts/structure_info.py",
            "scripts/fetch_alphafold.py",
            "scripts/pdb_info.py",
            "scripts/pymol_agent.py",
            "scripts/chimerax_agent.py",
            "scripts/proteus_doctor.py",
            "scripts/pae_report.py",
            "scripts/resolve_structure.py",
            "scripts/validation_report.py",
            "scripts/pocket_report.py",
            "scripts/compare_structures.py",
        ]:
            proc = run_script(script, "--help")
            self.assertIn("usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
