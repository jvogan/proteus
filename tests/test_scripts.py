#!/usr/bin/env python3
import ast
import json
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import add_helix_records
import assembly_report
import chimerax_agent
import fetch_alphafold
import fetch_pdb
import interface_report
import map_info
import proteus_doctor
import pymol_agent
import resolve_structure
import sifts_map
import uniprot_lookup
import variant_map


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
        self.assertEqual(data["file"], "./tests/fixtures/tiny.pdb")
        self.assertEqual(data["format"], "pdb")
        self.assertEqual(data["atom_records"], 3)
        self.assertEqual(data["hetatm_records"], 1)
        self.assertEqual(data["chains"], ["A", "B"])

    def test_structure_info_mmcif_json(self):
        proc = run_script("scripts/structure_info.py", "tests/fixtures/tiny.cif", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertIn("data", data)
        self.assertEqual(data["file"], "./tests/fixtures/tiny.cif")
        self.assertEqual(data["format"], "mmcif")
        self.assertEqual(data["atom_records"], 2)
        self.assertEqual(data["hetatm_records"], 1)
        self.assertEqual(data["title"], "Tiny mmCIF test structure")

    def test_pdb_info_json_contract(self):
        proc = run_script("scripts/pdb_info.py", "tests/fixtures/tiny.pdb", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertIn("data", data)
        self.assertEqual(data["file"], "./tests/fixtures/tiny.pdb")
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
        self.assertEqual(data["structure_path"], "./tests/fixtures/tiny.pdb")
        self.assertEqual(data["inspection"]["file"], "./tests/fixtures/tiny.pdb")
        self.assertEqual(data["inspection"]["format"], "pdb")

    def test_resolve_pdb_reports_biological_assembly_warning(self):
        metadata = {
            "pdb_id": "4HHB",
            "title": "Example structure",
            "assembly_count": 2,
            "polymer_entity_count": 3,
            "deposited_atom_count": 1234,
            "deposited_model_count": 1,
        }

        def fake_download(_url, path):
            path.write_text("data_4HHB\n", encoding="utf-8")

        def fake_inspect(path, force_alphafold=False):
            return {
                "status": "ok",
                "data": {
                    "file": str(Path(path).resolve()),
                    "format": "mmcif",
                    "force_alphafold": force_alphafold,
                },
                "file": str(Path(path).resolve()),
                "format": "mmcif",
            }

        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(resolve_structure.fetch_pdb, "fetch_entry_metadata", return_value=metadata), \
                    patch.object(
                        resolve_structure.fetch_pdb,
                        "build_download_url",
                        return_value=("https://files.rcsb.org/download/4HHB.cif", "4HHB.cif"),
                    ), \
                    patch.object(resolve_structure.fetch_pdb, "download", side_effect=fake_download), \
                    patch.object(resolve_structure.structure_info, "inspect_structure", side_effect=fake_inspect):
                out = resolve_structure.resolve("4hhb", source="pdb", outdir=tmp, download=True)

        self.assertEqual(out["status"], "ok")
        self.assertEqual(out["biological_assembly"]["assembly_count"], 2)
        self.assertEqual(out["biological_assembly"]["recommended_assembly"]["filename"], "4HHB-assembly1.cif")
        self.assertEqual(out["download"]["coordinate_scope"], "asymmetric_unit")
        self.assertEqual(out["structure_path"], "4HHB.cif (absolute path omitted)")
        self.assertEqual(out["inspection"]["file"], "4HHB.cif (absolute path omitted)")
        self.assertEqual(out["inspection"]["data"]["file"], "4HHB.cif (absolute path omitted)")
        self.assertIn("asymmetric unit", " ".join(out["warnings"]))
        self.assertIn("4HHB-assembly1.cif", " ".join(out["warnings"]))
        self.assertNotIn(tmp, json.dumps(out))

    def test_resolve_pdb_no_download_reports_assembly_without_warning(self):
        metadata = {
            "pdb_id": "1ABC",
            "title": "Metadata only",
            "assembly_count": 1,
        }
        with patch.object(resolve_structure.fetch_pdb, "fetch_entry_metadata", return_value=metadata), \
                patch.object(resolve_structure.fetch_pdb, "download") as download:
            out = resolve_structure.resolve("1abc", source="pdb", download=False)

        download.assert_not_called()
        self.assertEqual(out["status"], "ok")
        self.assertEqual(out["biological_assembly"]["assembly_downloads"][0]["filename"], "1ABC-assembly1.cif")
        self.assertEqual(out["structure_path"], None)
        self.assertEqual(out["warnings"], [])

    def test_variant_map_parses_supported_notations(self):
        parsed = variant_map.parse_variant("P04637:p.Arg175His")
        self.assertEqual(parsed["uniprot_id"], "P04637")
        self.assertEqual(parsed["short"], "R175H")
        self.assertEqual(parsed["hgvs_protein"], "p.Arg175His")
        self.assertEqual(parsed["residue_index"], 175)

        parsed = variant_map.parse_variant("R175H", uniprot_id="P04637")
        self.assertEqual(parsed["uniprot_id"], "P04637")
        self.assertEqual(parsed["from"]["three_letter"], "ARG")
        self.assertEqual(parsed["to"]["three_letter"], "HIS")

    def test_variant_map_local_pdb_json(self):
        proc = run_script(
            "scripts/variant_map.py",
            "G1A",
            "--uniprot",
            "P04637",
            "--structure",
            "tests/fixtures/tiny.pdb",
            "--json",
        )
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["variant"]["short"], "G1A")
        self.assertEqual(data["residue_index"], 1)
        self.assertEqual(data["structure"]["structure_path"], "./tests/fixtures/tiny.pdb")
        self.assertEqual(data["structure"]["inspection"]["file"], "./tests/fixtures/tiny.pdb")
        self.assertEqual(data["local_residue"]["format"], "pdb")
        self.assertTrue(data["local_residue"]["residue_present"])
        self.assertTrue(data["local_residue"]["ca_present"])
        selected = data["local_residue"]["selected"]
        self.assertEqual(selected["resname"], "GLY")
        self.assertEqual(selected["ca_coordinate"], {"x": 12.56, "y": 13.205, "z": 2.1})
        self.assertTrue(selected["reference_matches_variant"])
        self.assertIn("SIFTS", " ".join(data["warnings"]))

    def test_variant_map_local_mmcif_json(self):
        proc = run_script(
            "scripts/variant_map.py",
            "p.Gly1Ala",
            "--structure",
            "tests/fixtures/tiny.cif",
            "--json",
        )
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["variant"]["short"], "G1A")
        self.assertEqual(data["local_residue"]["format"], "mmcif")
        self.assertEqual(data["local_residue"]["selected"]["chain"], "A")
        self.assertEqual(data["local_residue"]["selected"]["ca_coordinate"]["x"], 12.56)

    def test_pocket_report_local_json(self):
        proc = run_script("scripts/pocket_report.py", "tests/fixtures/tiny.pdb", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["file"], "./tests/fixtures/tiny.pdb")
        self.assertEqual(data["ligand_count"], 1)
        self.assertGreaterEqual(data["ligands"][0]["contact_residue_count"], 1)

    def test_pocket_report_local_mmcif_json(self):
        proc = run_script("scripts/pocket_report.py", "tests/fixtures/tiny.cif", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["file"], "./tests/fixtures/tiny.cif")
        self.assertEqual(data["format"], "mmcif")
        self.assertEqual(data["ligand_count"], 1)
        self.assertEqual(data["ligands"][0]["ligand"]["resname"], "LIG")
        self.assertGreaterEqual(data["ligands"][0]["contact_residue_count"], 1)

    def test_pocket_report_residue_focus(self):
        proc = run_script(
            "scripts/pocket_report.py",
            "tests/fixtures/tiny.pdb",
            "--residue",
            "A:1",
            "--json",
        )
        data = json.loads(proc.stdout)
        focus = data["residue_focus"]
        self.assertTrue(focus["within_pocket_cutoff"])
        self.assertEqual(focus["candidates"][0]["nearest_ligand_contact"]["ligand"]["resname"], "LIG")
        self.assertIn("SIFTS", " ".join(data["warnings"]))

    def test_ligand_extract_local_json_filters_water(self):
        proc = run_script("scripts/ligand_extract.py", "tests/fixtures/ligands.pdb", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertIn("data", data)
        self.assertEqual(data["input"]["query"], "./tests/fixtures/ligands.pdb")
        self.assertEqual(data["file"], "./tests/fixtures/ligands.pdb")
        self.assertEqual(data["hetatm_records"], 5)
        self.assertEqual(data["water_hetatm_records"], 1)
        self.assertEqual(data["ligand_group_count"], 2)
        self.assertEqual(data["ligand_atom_count"], 4)
        self.assertEqual([group["ligand"] for group in data["ligand_groups"]], ["LIG", "SO4"])
        self.assertEqual(data["ligand_groups"][0]["residue_label"], "B:2")
        self.assertEqual(data["ligand_groups"][1]["residue_label"], "C:42A")

    def test_ligand_extract_ligand_filter(self):
        proc = run_script(
            "scripts/ligand_extract.py",
            "tests/fixtures/ligands.pdb",
            "--ligand",
            "lig",
            "--json",
        )
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["ligand_filter"], ["LIG"])
        self.assertEqual(data["ligand_group_count"], 1)
        self.assertEqual(data["ligand_components"][0]["ligand"], "LIG")
        self.assertEqual(data["ligand_components"][0]["atom_count"], 2)

    def test_ligand_extract_local_mmcif_json(self):
        proc = run_script("scripts/ligand_extract.py", "tests/fixtures/tiny.cif", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["format"], "mmcif")
        self.assertEqual(data["ligand_group_count"], 1)
        self.assertEqual(data["ligand_groups"][0]["ligand"], "LIG")
        self.assertEqual(data["ligand_groups"][0]["residue_label"], "B:2")

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
        self.assertEqual(pymol_agent._validate_pymol_color("plddt"), "plddt")
        with self.assertRaises(ValueError):
            pymol_agent._validate_pymol_color("red; import os")
        result = pymol_agent.render_structure("tests/fixtures/tiny.pdb", "out.png", color="red; import os")
        self.assertEqual(result["status"], "error")

    def test_pymol_plddt_color_script_is_layered(self):
        script = pymol_agent._color_script("plddt")
        # Broadest-first layering (no <= in PyMOL selection algebra), with
        # confidence-scale detection for normalized 0-1 or AlphaFold 0-100 B-factors.
        self.assertIn('cmd.color("orange"', script)
        self.assertIn("max(_proteus_b_values) <= 1.5", script)
        self.assertIn("(0.50, 0.70, 0.90)", script)
        self.assertIn("(50.0, 70.0, 90.0)", script)
        self.assertLess(script.index('cmd.color("orange"'), script.index('cmd.color("yellow"'))
        self.assertLess(script.index('cmd.color("yellow"'), script.index('cmd.color("blue"'))

    def test_chimerax_rest_color_validation(self):
        import chimerax_rest
        self.assertEqual(chimerax_rest._validate_color("plddt"), "plddt")
        with self.assertRaises(ValueError):
            chimerax_rest._validate_color("red; close session")
        result = chimerax_rest.rest_render("tests/fixtures/tiny.pdb", "out.png", color="red; close")
        self.assertEqual(result["status"], "error")

    def test_add_helix_detects_ideal_helix(self):
        import math
        coords = []
        for i in range(1, 15):
            ang = math.radians(100 * i)
            coords.append((i, "A", 2.3 * math.cos(ang), 2.3 * math.sin(ang), 1.5 * i))
        is_helix = add_helix_records.detect_helix_residues(coords)
        self.assertTrue(all(is_helix))
        segments = add_helix_records.helix_segments(coords, is_helix, min_len=6)
        self.assertEqual(len(segments), 1)
        self.assertEqual(segments[0], ("A", 1, 14))

    def test_interface_report_pair_detection(self):
        atoms_a = [{"atom": "CA", "resname": "ALA", "resi": "1", "x": 0, "y": 0, "z": 0},
                   {"atom": "CA", "resname": "ARG", "resi": "2", "x": 3, "y": 0, "z": 0}]
        atoms_b = [{"atom": "CA", "resname": "ASP", "resi": "1", "x": 3, "y": 3, "z": 0},
                   {"atom": "CA", "resname": "GLU", "resi": "2", "x": 10, "y": 10, "z": 0}]
        res_a, res_b, pair_min = interface_report.analyze_pair(atoms_a, atoms_b, 5.0)
        self.assertEqual(set(res_a.keys()), {("1", "ALA"), ("2", "ARG")})
        self.assertEqual(set(res_b.keys()), {("1", "ASP")})
        self.assertAlmostEqual(min(pair_min.values()), 3.0, places=2)

    def test_interface_report_single_chain_note(self):
        out = interface_report.analyze_interfaces("tests/fixtures/tiny.pdb")
        self.assertEqual(out["status"], "ok")
        self.assertEqual(out["file"], "./tests/fixtures/tiny.pdb")
        self.assertEqual(out["data"]["chains"], ["A"])
        self.assertEqual(out["data"]["interface_count"], 0)
        self.assertIn("note", out["data"])

    def test_interface_report_local_mmcif_json(self):
        proc = run_script("scripts/interface_report.py", "tests/fixtures/tiny.cif", "--json")
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["file"], "./tests/fixtures/tiny.cif")
        self.assertEqual(data["format"], "mmcif")
        self.assertEqual(data["chains"], ["A"])
        self.assertEqual(data["interface_count"], 0)
        self.assertIn("note", data)

    def test_interface_report_residue_focus_single_chain(self):
        proc = run_script(
            "scripts/interface_report.py",
            "tests/fixtures/tiny.pdb",
            "--variant",
            "G1A",
            "--json",
        )
        data = json.loads(proc.stdout)
        focus = data["residue_focus"]
        self.assertFalse(focus["participates_in_interface"])
        self.assertEqual(focus["candidate_count"], 1)
        self.assertTrue(focus["candidates"][0]["reference_matches_variant"])
        self.assertIn("SIFTS", " ".join(data["warnings"]))

    def test_assembly_report_builds_downloads_and_recommendation(self):
        entry = {
            "struct": {"title": "Example structure"},
            "rcsb_entry_info": {
                "assembly_count": 2,
                "polymer_entity_count": 3,
                "deposited_atom_count": 1234,
                "deposited_model_count": 1,
            },
        }
        out = assembly_report.build_report("4HHB", entry)
        self.assertEqual(out["status"], "ok")
        self.assertIn("data", out)
        self.assertEqual(out["assembly_count"], 2)
        self.assertEqual(out["polymer_entity_count"], 3)
        self.assertEqual(out["deposited_atom_count"], 1234)
        self.assertEqual(out["deposited_model_count"], 1)
        self.assertEqual(
            [item["filename"] for item in out["assembly_downloads"]],
            ["4HHB-assembly1.cif", "4HHB-assembly2.cif"],
        )
        self.assertEqual(out["recommended_assembly"]["assembly_id"], "1")

    def test_assembly_report_can_select_requested_assembly(self):
        entry = {"rcsb_entry_info": {"assembly_count": 2}}
        out = assembly_report.build_report("4HHB", entry, requested_assembly=2)
        self.assertEqual(out["recommended_assembly"]["assembly_id"], "1")
        self.assertEqual(out["selected_assembly"]["assembly_id"], "2")
        with self.assertRaises(assembly_report.AssemblyReportError):
            assembly_report.build_report("4HHB", entry, requested_assembly=3)

    def test_assembly_report_download_path_is_scrubbed(self):
        entry = {"rcsb_entry_info": {"assembly_count": 1}}
        report = assembly_report.build_report("4HHB", entry)
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(assembly_report, "_download", side_effect=lambda _url, path: path.write_text("x")):
                out = assembly_report.download_selected_assembly(report, tmp)
        self.assertEqual(out["download"]["path"], "4HHB-assembly1.cif (absolute path omitted)")

    def test_pymol_verify_png(self):
        err = {"status": "error", "error": "boom"}
        self.assertIs(pymol_agent._verify_png(err, "/no/such.png"), err)
        downgraded = pymol_agent._verify_png({"status": "ok", "data": {}}, "/no/such.png")
        self.assertEqual(downgraded["status"], "error")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as fh:
            fh.write(b"\x89PNG")
            path = fh.name
        try:
            kept = pymol_agent._verify_png({"status": "ok", "data": {}}, path)
            self.assertEqual(kept["status"], "ok")
        finally:
            Path(path).unlink()

    def test_pymol_pocket_missing_file(self):
        result = pymol_agent.render_pocket("tests/fixtures/does_not_exist.pdb", "out.png")
        self.assertEqual(result["status"], "error")

    def test_pymol_density_requires_map_or_simulate(self):
        # Missing model file
        self.assertEqual(
            pymol_agent.render_density("tests/fixtures/missing.pdb", "out.png", simulate=True)["status"],
            "error")
        # Neither --map nor --simulate
        result = pymol_agent.render_density("tests/fixtures/tiny.pdb", "out.png")
        self.assertEqual(result["status"], "error")
        self.assertIn("simulate", result["error"])
        # Map path that does not exist
        self.assertEqual(
            pymol_agent.render_density("tests/fixtures/tiny.pdb", "out.png", map_path="nope.mrc")["status"],
            "error")

    def test_map_info_sigma_from_synthetic_mrc(self):
        nx = ny = nz = 4
        vals = [float(i % 7) for i in range(nx * ny * nz)]
        header = bytearray(1024)
        struct.pack_into("<4i", header, 0, nx, ny, nz, 2)  # dims + float32 mode
        struct.pack_into("<i", header, 92, 0)              # nsymbt = 0
        body = b"".join(struct.pack("<f", v) for v in vals)
        with tempfile.NamedTemporaryFile(suffix=".mrc", delete=False) as fh:
            fh.write(bytes(header) + body)
            path = fh.name
        try:
            dims, mode, mean, sigma = map_info.read_map_stats(path)
        finally:
            Path(path).unlink()
        expected_mean = sum(vals) / len(vals)
        self.assertEqual(dims, (4, 4, 4))
        self.assertEqual(mode, 2)
        self.assertAlmostEqual(mean, expected_mean, places=4)
        self.assertGreater(sigma, 0)

    def test_snapshot_json_files_parse(self):
        for path in (ROOT / "docs" / "snapshots").glob("*.json"):
            data = json.loads(path.read_text())
            self.assertIn(data["status"], {"ok", "error"})

    def test_error_payloads_are_json_parseable(self):
        for payload in [
            fetch_alphafold._error_payload("bad alphafold input"),
            fetch_pdb._error_payload("bad pdb input"),
            sifts_map._error_payload("bad sifts input"),
            uniprot_lookup._error_payload("bad uniprot input"),
        ]:
            encoded = json.dumps(payload)
            decoded = json.loads(encoded)
            self.assertEqual(decoded["status"], "error")
            self.assertIn("error", decoded)

    def test_sifts_map_normalizes_pdb_payload(self):
        raw = {
            "1hsg": {
                "UniProt": {
                    "P03367": {
                        "identifier": "POL_HV1BR",
                        "mappings": [
                            {
                                "entity_id": 1,
                                "chain_id": "A",
                                "struct_asym_id": "A",
                                "unp_start": 501,
                                "unp_end": 599,
                                "start": {
                                    "author_residue_number": 1,
                                    "author_insertion_code": "",
                                    "residue_number": 1,
                                },
                                "end": {
                                    "author_residue_number": 99,
                                    "author_insertion_code": "",
                                    "residue_number": 99,
                                },
                                "identity": 1.0,
                                "coverage": 0.068,
                            }
                        ],
                    }
                }
            }
        }
        records = sifts_map.normalize_pdb_mappings(raw, "1hsg")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["pdb_id"], "1hsg")
        self.assertEqual(records[0]["chain_id"], "A")
        self.assertEqual(records[0]["uniprot_accession"], "P03367")
        self.assertEqual(records[0]["uniprot_id"], "POL_HV1BR")
        self.assertEqual(records[0]["uniprot_start"], 501)
        self.assertEqual(records[0]["pdb_end"], 99)

    def test_sifts_map_normalizes_uniprot_payload(self):
        raw = {
            "P04637": {
                "PDB": {
                    "1tup": [
                        {
                            "entity_id": 1,
                            "chain_id": "A",
                            "struct_asym_id": "A",
                            "unp_start": 94,
                            "unp_end": 289,
                            "start": {
                                "author_residue_number": 94,
                                "author_insertion_code": "",
                                "residue_number": 1,
                            },
                            "end": {
                                "author_residue_number": 289,
                                "author_insertion_code": "",
                                "residue_number": 196,
                            },
                        }
                    ]
                }
            }
        }
        records = sifts_map.normalize_uniprot_mappings(raw, "P04637")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["pdb_id"], "1tup")
        self.assertEqual(records[0]["uniprot_accession"], "P04637")
        self.assertEqual(records[0]["uniprot_start"], 94)
        self.assertEqual(records[0]["auth_end"], 289)

    def test_sifts_json_loader_maps_uniprot_residue_to_auth_residue(self):
        raw = {
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
        }
        with tempfile.TemporaryDirectory() as tmp:
            sifts_path = Path(tmp) / "sifts.json"
            sifts_path.write_text(json.dumps(raw))
            records = sifts_map.load_sifts_json(sifts_path)

        candidates = sifts_map.map_uniprot_residue_candidates(records, "P04637", 175)
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]["structure_chain_id"], "A")
        self.assertEqual(candidates[0]["structure_residue_id"], "1")
        self.assertEqual(candidates[0]["auth_residue_id"], "1")

    def test_variant_map_sifts_json_maps_uniprot_index_to_local_auth_residue(self):
        raw = {
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
        }
        with tempfile.TemporaryDirectory() as tmp:
            sifts_path = Path(tmp) / "sifts.json"
            sifts_path.write_text(json.dumps(raw))
            proc = run_script(
                "scripts/variant_map.py",
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
        self.assertTrue(data["sifts_mapping"]["applied"])
        self.assertEqual(data["local_residue"]["lookup_method"], "sifts_uniprot_to_structure_mapping")
        selected = data["local_residue"]["selected"]
        self.assertEqual(selected["chain"], "A")
        self.assertEqual(selected["residue_id"], "1")
        self.assertTrue(selected["reference_matches_variant"])
        self.assertNotIn(str(sifts_path), proc.stdout)
        self.assertNotIn("direct structure residue numbering", " ".join(data["warnings"]))

    def test_sifts_map_invalid_pdb_json_error(self):
        proc = subprocess.run(
            [sys.executable, "scripts/sifts_map.py", "pdb", "bad", "--json"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(proc.returncode, 0)
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "error")
        self.assertIn("PDB ID", data["error"])

    def test_doctor_network_capabilities_follow_checks(self):
        with patch.object(proteus_doctor, "_script_smoke", return_value={}):
            no_network = proteus_doctor.build_report(False)["data"]
        self.assertEqual(no_network["root"], ".")
        self.assertNotIn(str(ROOT), json.dumps(no_network))
        self.assertIsNone(no_network["network"])
        self.assertFalse(no_network["capabilities"]["rcsb_fetch"])
        self.assertFalse(no_network["capabilities"]["uniprot_lookup"])
        self.assertFalse(no_network["capabilities"]["alphafold_fetch"])

        def fake_network_check(url):
            if "data.rcsb.org" in url:
                return {"ok": True}
            if "uniprot.org" in url:
                return {"ok": False, "error": "blocked"}
            if "alphafold.ebi.ac.uk" in url:
                return {"ok": True}
            return {"ok": False, "error": "unexpected"}

        with patch.object(proteus_doctor, "_script_smoke", return_value={}), \
                patch.object(proteus_doctor, "_network_check", side_effect=fake_network_check):
            with_network = proteus_doctor.build_report(True)["data"]
        self.assertTrue(with_network["network"]["rcsb"]["ok"])
        self.assertFalse(with_network["network"]["uniprot"]["ok"])
        self.assertTrue(with_network["network"]["alphafold"]["ok"])
        self.assertTrue(with_network["capabilities"]["rcsb_fetch"])
        self.assertFalse(with_network["capabilities"]["uniprot_lookup"])
        self.assertTrue(with_network["capabilities"]["alphafold_fetch"])

    def test_cli_help(self):
        for script in [
            "scripts/fetch_pdb.py",
            "scripts/pdb_search.py",
            "scripts/pdb_select.py",
            "scripts/uniprot_lookup.py",
            "scripts/structure_info.py",
            "scripts/fetch_alphafold.py",
            "scripts/pdb_info.py",
            "scripts/pymol_agent.py",
            "scripts/chimerax_agent.py",
            "scripts/chimerax_rest.py",
            "scripts/proteus_doctor.py",
            "scripts/pae_report.py",
            "scripts/resolve_structure.py",
            "scripts/validation_report.py",
            "scripts/pocket_report.py",
            "scripts/compare_structures.py",
            "scripts/add_helix_records.py",
            "scripts/interface_report.py",
            "scripts/map_info.py",
            "scripts/model_quality.py",
            "scripts/design_run.py",
            "scripts/rosetta_score.py",
            "scripts/proteus_batch.py",
            "scripts/proteus_cache.py",
            "scripts/proteus_report.py",
            "scripts/target_dossier.py",
            "scripts/sifts_map.py",
            "scripts/assembly_report.py",
            "scripts/kras_dossier.py",
            "scripts/ligand_extract.py",
            "scripts/dock_prep.py",
            "scripts/docking_box.py",
            "scripts/dock_vina.py",
            "scripts/interaction_report.py",
            "scripts/variant_map.py",
            "scripts/mutation_triage.py",
        ]:
            proc = run_script(script, "--help")
            self.assertIn("usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
