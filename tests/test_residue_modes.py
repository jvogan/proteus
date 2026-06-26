#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import interface_report
import pocket_report


def run_script(*args):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )


def pdb_line(record, serial, atom, resname, chain, resi, x, y, z, element="C"):
    return (
        f"{record:<6}{serial:5d} {atom:>4} {resname:>3} {chain:1}{resi:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{50.00:6.2f}          {element:>2}\n"
    )


def write_interface_fixture(path: Path):
    path.write_text(
        "".join([
            "TITLE     RESIDUE MODE INTERFACE FIXTURE\n",
            pdb_line("ATOM", 1, "N", "ARG", "A", 10, 0.0, 0.0, 0.0, "N"),
            pdb_line("ATOM", 2, "CA", "ARG", "A", 10, 1.0, 0.0, 0.0, "C"),
            pdb_line("ATOM", 3, "CB", "ARG", "A", 10, 1.0, 1.0, 0.0, "C"),
            pdb_line("ATOM", 4, "CA", "SER", "B", 7, 1.0, 2.4, 0.0, "C"),
            pdb_line("ATOM", 5, "CA", "GLY", "B", 8, 10.0, 10.0, 0.0, "C"),
            "END\n",
        ])
    )


class ResidueModeTests(unittest.TestCase):
    def test_pocket_residue_focus_cli_reports_ligand_cutoff(self):
        proc = run_script(
            "scripts/pocket_report.py",
            "tests/fixtures/tiny.pdb",
            "--residue",
            "1",
            "--variant",
            "G1A",
            "--radius",
            "2",
            "--json",
        )
        data = json.loads(proc.stdout)

        self.assertEqual(data["status"], "ok")
        self.assertIn("SIFTS", " ".join(data["warnings"]))
        focus = data["residue_focus"]
        self.assertTrue(focus["within_pocket_cutoff"])
        self.assertEqual(focus["selection"]["chain"], None)
        self.assertEqual(focus["variant"]["short"], "G1A")

        candidate = focus["candidates"][0]
        self.assertEqual(candidate["chain"], "A")
        self.assertTrue(candidate["reference_matches_variant"])
        contact = candidate["nearest_ligand_contact"]
        self.assertEqual(contact["ligand"]["resname"], "LIG")
        self.assertTrue(contact["within_cutoff"])

    def test_pocket_residue_focus_api_can_use_variant_index(self):
        output = pocket_report.analyze_pocket(
            "tests/fixtures/tiny.pdb",
            radius=2.0,
            variant="G1A",
        )

        self.assertEqual(output["status"], "ok")
        focus = output["residue_focus"]
        self.assertEqual(focus["selection"]["source"], "variant")
        self.assertTrue(focus["within_pocket_cutoff"])

    def test_interface_residue_focus_cli_reports_opposite_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb_path = Path(tmp) / "interface.pdb"
            write_interface_fixture(pdb_path)
            proc = run_script(
                "scripts/interface_report.py",
                str(pdb_path),
                "--chains",
                "A,B",
                "--residue",
                "A:10",
                "--variant",
                "R10H",
                "--cutoff",
                "2",
                "--json",
            )

        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertIn("interface.pdb", data["file"])
        self.assertNotIn(str(pdb_path), proc.stdout)
        self.assertEqual(data["interface_count"], 1)
        self.assertIn("SIFTS", " ".join(data["warnings"]))

        focus = data["residue_focus"]
        self.assertTrue(focus["participates_in_interface"])
        candidate = focus["candidates"][0]
        self.assertTrue(candidate["reference_matches_variant"])
        self.assertEqual(candidate["opposite_contacts_within_cutoff_count"], 1)
        contact = candidate["nearest_opposite_chain_contact"]
        self.assertEqual(contact["chain"], "B")
        self.assertEqual(contact["resi"], "7")
        self.assertTrue(contact["within_cutoff"])

    def test_interface_residue_focus_api_accepts_direct_residue_number(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb_path = Path(tmp) / "interface.pdb"
            write_interface_fixture(pdb_path)
            output = interface_report.analyze_interfaces(
                str(pdb_path),
                cutoff=2.0,
                residue="10",
                variant="R10H",
            )

        self.assertEqual(output["status"], "ok")
        focus = output["residue_focus"]
        self.assertEqual(focus["selection"]["chain"], None)
        self.assertTrue(focus["participates_in_interface"])
        self.assertEqual(focus["candidate_count"], 1)


if __name__ == "__main__":
    unittest.main()
