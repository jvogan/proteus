#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import mutation_triage


def run_script(*args, check=True):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


def pdb_line(record, serial, atom, resname, chain, resi, x, y, z, element="C"):
    return (
        f"{record:<6}{serial:5d} {atom:>4} {resname:>3} {chain:1}{resi:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{50.00:6.2f}          {element:>2}\n"
    )


def write_contact_fixture(path: Path):
    path.write_text(
        "".join([
            "TITLE     MUTATION TRIAGE CONTACT FIXTURE\n",
            pdb_line("ATOM", 1, "N", "ARG", "A", 10, 0.0, 0.0, 0.0, "N"),
            pdb_line("ATOM", 2, "CA", "ARG", "A", 10, 1.0, 0.0, 0.0, "C"),
            pdb_line("ATOM", 3, "CB", "ARG", "A", 10, 1.0, 1.0, 0.0, "C"),
            pdb_line("ATOM", 4, "N", "GLY", "A", 11, 2.2, 1.0, 0.0, "N"),
            pdb_line("ATOM", 5, "CA", "GLY", "A", 11, 2.6, 1.0, 0.0, "C"),
            pdb_line("ATOM", 6, "CA", "SER", "B", 7, 1.0, 2.4, 0.0, "C"),
            pdb_line("HETATM", 7, "P", "ATP", "C", 101, 1.0, 1.5, 0.0, "P"),
            pdb_line("HETATM", 8, "O", "HOH", "D", 1, 1.0, 1.2, 0.0, "O"),
            "END\n",
        ])
    )


def write_pae_fixture(path: Path, size=12):
    matrix = [
        [abs(row - column) + 0.5 for column in range(size)]
        for row in range(size)
    ]
    path.write_text(json.dumps({"predicted_aligned_error": matrix}))


def write_sifts_fixture(path: Path):
    path.write_text(json.dumps({
        "1abc": {
            "UniProt": {
                "P04637": {
                    "identifier": "TP53_HUMAN",
                    "mappings": [
                        {
                            "entity_id": 1,
                            "chain_id": "A",
                            "struct_asym_id": "A",
                            "unp_start": 175,
                            "unp_end": 175,
                            "start": {
                                "author_residue_number": 1,
                                "author_insertion_code": "",
                                "residue_number": 1,
                            },
                            "end": {
                                "author_residue_number": 1,
                                "author_insertion_code": "",
                                "residue_number": 1,
                            },
                        }
                    ],
                }
            }
        }
    }))


class MutationTriageTests(unittest.TestCase):
    def test_help_exposes_json_and_structure_options(self):
        proc = run_script("scripts/mutation_triage.py", "--help")

        self.assertEqual(proc.returncode, 0)
        self.assertIn("--json", proc.stdout)
        self.assertIn("--structure", proc.stdout)
        self.assertIn("direct structure numbering", proc.stdout)

    def test_json_reports_contacts_ligands_absent_variant_and_pae(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb_path = Path(tmp) / "contacts.pdb"
            pae_path = Path(tmp) / "pae.json"
            write_contact_fixture(pdb_path)
            write_pae_fixture(pae_path)

            proc = run_script(
                "scripts/mutation_triage.py",
                "R10H",
                "G99A",
                "--uniprot",
                "P04637",
                "--structure",
                str(pdb_path),
                "--cutoff",
                "2.0",
                "--pae",
                str(pae_path),
                "--json",
            )

        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["structure"]["ligand_group_count"], 1)
        self.assertIn("direct structure residue numbering", " ".join(data["warnings"]))
        self.assertEqual(len(data["variants"]), 2)

        hit = data["variants"][0]
        self.assertEqual(hit["variant"]["short"], "R10H")
        self.assertEqual(hit["candidate_count"], 1)
        self.assertIn("same_chain_close_contacts", hit["flags"])
        self.assertIn("near_ligand", hit["flags"])
        self.assertIn("near_other_chain", hit["flags"])

        candidate = hit["candidates"][0]
        self.assertEqual(candidate["residue_label"], "A:10")
        self.assertTrue(candidate["reference_matches_variant"])
        self.assertEqual(candidate["ca_coordinate"], {"x": 1.0, "y": 0.0, "z": 0.0})
        self.assertEqual(candidate["nearby_ligands"][0]["ligand"], "ATP")
        self.assertEqual(candidate["same_chain_contacts"][0]["residue_label"], "A:11")
        self.assertEqual(candidate["same_chain_contacts"][0]["pae_pair"]["mean_pair_pae"], 1.5)
        self.assertEqual(candidate["other_chain_summary"][0]["chain"], "B")
        self.assertEqual(candidate["pae"]["row_mean_pae"], 4.5)

        miss = data["variants"][1]
        self.assertEqual(miss["variant"]["short"], "G99A")
        self.assertEqual(miss["candidate_count"], 0)
        self.assertIn("No structure residue matched", " ".join(miss["warnings"]))

    def test_mmcif_triage_uses_local_atom_parser(self):
        output = mutation_triage.triage_variants(
            str(ROOT / "tests/fixtures/tiny.cif"),
            ["G1A"],
            cutoff=4.0,
        )

        self.assertEqual(output["status"], "ok")
        self.assertEqual(output["structure"]["format"], "mmcif")
        candidate = output["variants"][0]["candidates"][0]
        self.assertEqual(candidate["residue_label"], "A:1")
        self.assertEqual(candidate["ca_coordinate"]["x"], 12.56)
        self.assertEqual(candidate["nearby_ligands"][0]["ligand"], "LIG")

    def test_sifts_json_maps_uniprot_index_to_structure_auth_numbering(self):
        with tempfile.TemporaryDirectory() as tmp:
            sifts_path = Path(tmp) / "sifts.json"
            write_sifts_fixture(sifts_path)
            proc = run_script(
                "scripts/mutation_triage.py",
                "G175A",
                "--uniprot",
                "P04637",
                "--structure",
                "tests/fixtures/tiny.pdb",
                "--sifts-json",
                str(sifts_path),
                "--json",
            )

        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["numbering"]["method"], "sifts_uniprot_to_structure_mapping")
        self.assertTrue(data["numbering"]["sifts_mapping"]["applied"])
        self.assertNotIn("direct structure residue numbering", " ".join(data["warnings"]))

        variant = data["variants"][0]
        self.assertTrue(variant["sifts_mapping"]["applied"])
        self.assertEqual(variant["candidate_count"], 1)
        candidate = variant["candidates"][0]
        self.assertEqual(candidate["residue_label"], "A:1")
        self.assertTrue(candidate["reference_matches_variant"])
        self.assertEqual(candidate["sifts_mapping"]["structure_residue_id"], "1")
        self.assertNotIn(str(sifts_path), proc.stdout)

    def test_default_output_is_readable_markdown(self):
        proc = run_script(
            "scripts/mutation_triage.py",
            "G1A",
            "--structure",
            "tests/fixtures/tiny.pdb",
        )

        self.assertTrue(proc.stdout.startswith("# Mutation Triage"))
        self.assertIn("## G1A", proc.stdout)
        self.assertIn("Ligands within cutoff", proc.stdout)
        self.assertIn("SIFTS mapping was not applied", proc.stdout)

    def test_absolute_paths_are_scrubbed_from_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb_path = Path(tmp) / "contacts.pdb"
            pae_path = Path(tmp) / "pae.json"
            write_contact_fixture(pdb_path)
            write_pae_fixture(pae_path)

            json_proc = run_script(
                "scripts/mutation_triage.py",
                "R10H",
                "--structure",
                str(pdb_path),
                "--pae",
                str(pae_path),
                "--json",
            )
            markdown_proc = run_script(
                "scripts/mutation_triage.py",
                "R10H",
                "--structure",
                str(pdb_path),
                "--pae",
                str(pae_path),
            )

            self.assertIn("contacts.pdb", json_proc.stdout)
            self.assertIn("pae.json", json_proc.stdout)
            self.assertIn("contacts.pdb", markdown_proc.stdout)
            self.assertIn("pae.json", markdown_proc.stdout)
            self.assertNotIn(str(pdb_path), json_proc.stdout)
            self.assertNotIn(str(pae_path), json_proc.stdout)
            self.assertNotIn(str(pdb_path), markdown_proc.stdout)
            self.assertNotIn(str(pae_path), markdown_proc.stdout)


if __name__ == "__main__":
    unittest.main()
