#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import interaction_report


def pdb_line(record, serial, atom, resname, chain, resi, x, y, z, element):
    return (
        f"{record:<6}{serial:5d} {atom:<4} {resname:>3} {chain:1}{resi:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2}\n"
    )


PDB_TEXT = "".join([
    "TITLE     INTERACTION REPORT TEST STRUCTURE\n",
    pdb_line("ATOM", 1, "C", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
    pdb_line("ATOM", 2, "O", "SER", "A", 2, 5.0, 0.0, 0.0, "O"),
    pdb_line("ATOM", 3, "N", "HIS", "A", 3, 9.0, 0.0, 0.0, "N"),
    pdb_line("HETATM", 4, "O", "HOH", "A", 99, 0.5, 0.0, 0.0, "O"),
    pdb_line("HETATM", 5, "C1", "LIG", "B", 10, 0.0, 0.0, 1.8, "C"),
    pdb_line("HETATM", 6, "O1", "LIG", "B", 10, 5.0, 3.0, 0.0, "O"),
    pdb_line("HETATM", 7, "ZN", "ZN", "C", 20, 7.1, 0.0, 0.0, "ZN"),
    "END\n",
])


NO_LIGAND_PDB = "".join([
    "TITLE     NO LIGAND TEST STRUCTURE\n",
    pdb_line("ATOM", 1, "C", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
    pdb_line("HETATM", 2, "O", "HOH", "A", 99, 0.5, 0.0, 0.0, "O"),
    "END\n",
])


def run_script(*args, check=True):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


class InteractionReportTests(unittest.TestCase):
    def write_pdb(self, text=PDB_TEXT):
        handle = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
        with handle:
            handle.write(text)
        self.addCleanup(lambda: Path(handle.name).exists() and Path(handle.name).unlink())
        return Path(handle.name)

    def test_analyze_local_pdb_classifies_contacts_by_ligand_group(self):
        path = self.write_pdb()

        report = interaction_report.analyze_interactions(str(path), max_contacts=20)

        self.assertEqual(report["status"], "ok")
        self.assertEqual(report["input"]["kind"], "local_file")
        self.assertEqual(report["source_label"], path.name)
        self.assertEqual(report["protein_atom_count"], 3)
        self.assertEqual(report["hetatm_records"], 4)
        self.assertEqual(report["water_hetatm_records"], 1)
        self.assertEqual(report["ligand_group_count"], 2)
        self.assertEqual(report["ligand_atom_count"], 3)

        groups = {group["ligand"]: group for group in report["ligand_groups"]}
        lig = groups["LIG"]
        self.assertGreaterEqual(lig["contact_counts"]["close_contact_clash"], 1)
        self.assertGreaterEqual(lig["contact_counts"]["polar_candidate"], 1)
        self.assertGreaterEqual(lig["contact_counts"]["hydrophobic_candidate"], 1)
        self.assertEqual(lig["contacting_residue_count"], 2)

        zinc = groups["ZN"]
        self.assertGreaterEqual(zinc["contact_counts"]["metal_candidate"], 1)
        self.assertIn("metal_candidate", zinc["closest_contacts"][0]["classifications"])
        self.assertGreaterEqual(report["classification_counts"]["metal_candidate"], 1)

    def test_ligand_filter_limits_selected_groups(self):
        path = self.write_pdb()

        report = interaction_report.analyze_interactions(str(path), ligand_filters=["ZN"])

        self.assertEqual(report["ligand_filter"], ["ZN"])
        self.assertEqual(report["ligand_group_count"], 1)
        self.assertEqual(report["ligand_groups"][0]["ligand"], "ZN")
        self.assertEqual(report["ligand_atom_count"], 1)

    def test_local_mmcif_classifies_contacts(self):
        report = interaction_report.analyze_interactions("tests/fixtures/tiny.cif", max_contacts=20)

        self.assertEqual(report["status"], "ok")
        self.assertEqual(report["format"], "mmcif")
        self.assertEqual(report["protein_atom_count"], 2)
        self.assertEqual(report["ligand_group_count"], 1)
        self.assertEqual(report["ligand_groups"][0]["ligand"], "LIG")
        self.assertEqual(report["ligand_groups"][0]["residue_label"], "B:2")
        self.assertGreaterEqual(report["ligand_groups"][0]["contact_count"], 1)

    def test_cli_json_uses_privacy_safe_path_labels(self):
        path = self.write_pdb()

        proc = run_script("scripts/interaction_report.py", str(path), "--json")
        data = json.loads(proc.stdout)

        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["source_label"], path.name)
        self.assertNotIn(str(path.parent), proc.stdout)
        self.assertNotIn(str(path), proc.stdout)

    def test_cli_markdown_is_readable_and_privacy_safe(self):
        path = self.write_pdb()

        proc = run_script("scripts/interaction_report.py", str(path))

        self.assertIn("# Protein-Ligand Interaction Report", proc.stdout)
        self.assertIn("## Ligand Groups", proc.stdout)
        self.assertIn("`polar_candidate`", proc.stdout)
        self.assertIn("LIG B:10", proc.stdout)
        self.assertNotIn(str(path.parent), proc.stdout)
        self.assertNotIn(str(path), proc.stdout)

    def test_optional_tool_detection_shape(self):
        tools = interaction_report.detect_optional_tools()

        self.assertEqual(set(tools), {"plip", "prolif"})
        for name in ("plip", "prolif"):
            self.assertIn("available", tools[name])
            self.assertIn("executable", tools[name])
            self.assertIn("module", tools[name])
            self.assertIsInstance(tools[name]["available"], bool)
            self.assertIsInstance(tools[name]["executable"]["ok"], bool)
            self.assertIn("path_label", tools[name]["executable"])
            self.assertIsInstance(tools[name]["module"]["ok"], bool)
            self.assertIn("origin_label", tools[name]["module"])

    def test_no_ligands_reports_json_error_without_network(self):
        path = self.write_pdb(NO_LIGAND_PDB)

        proc = run_script("scripts/interaction_report.py", str(path), "--json", check=False)
        data = json.loads(proc.stdout)

        self.assertEqual(proc.returncode, 1)
        self.assertEqual(data["status"], "error")
        self.assertIn("No non-water HETATM ligand atoms", data["error"])

    def test_missing_non_pdbid_input_errors_without_download(self):
        proc = run_script("scripts/interaction_report.py", "missing-local-file.pdb", "--json", check=False)
        data = json.loads(proc.stdout)

        self.assertEqual(proc.returncode, 1)
        self.assertEqual(data["status"], "error")
        self.assertIn("local PDB/mmCIF file", data["error"])


if __name__ == "__main__":
    unittest.main()
