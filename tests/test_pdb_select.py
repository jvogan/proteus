#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import pdb_select


def run_script(*args):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )


class PDBSelectTests(unittest.TestCase):
    def candidates(self):
        return [
            {
                "id": "2aaa",
                "title": "X-ray apo candidate",
                "resolution": [1.8],
                "method": "X-RAY DIFFRACTION",
                "ligands": ["SO4"],
                "assembly_count": 1,
                "validation": {"geometry": {"clashscore": 8.0}},
            },
            {
                "id": "3bbb",
                "title": "Electron microscopy ligand candidate",
                "resolution": [2.7],
                "method": "ELECTRON MICROSCOPY",
                "ligands": ["ATP", "MG"],
                "assembly_count": 2,
                "validation": {"geometry": {"clashscore": 4.0}},
            },
            {
                "id": "1ccc",
                "title": "Best ligand candidate",
                "resolution": [1.5],
                "method": "X-RAY DIFFRACTION",
                "ligands": ["ATP"],
                "assembly_count": 1,
                "validation": {
                    "geometry": {
                        "clashscore": 5.0,
                        "percent_ramachandran_outliers": 0.2,
                    }
                },
            },
            {
                "id": "4ddd",
                "title": "NMR ligand candidate",
                "method": "SOLUTION NMR",
                "ligands": ["ATP"],
                "assembly_count": 1,
            },
        ]

    def test_rank_candidates_prefers_lower_resolution(self):
        report = pdb_select.rank_candidates(self.candidates())
        self.assertEqual(report["best"]["pdb_id"], "1CCC")
        self.assertEqual([item["pdb_id"] for item in report["ranked"]], ["1CCC", "2AAA", "3BBB", "4DDD"])
        self.assertIn("resolution 1.5 A", " ".join(report["best"]["reasons"]))
        self.assertEqual(report["best"]["selection_key"]["method_preference_rank"], 1)

    def test_ligand_filter_excludes_missing_ligands(self):
        filters = pdb_select.normalize_ligand_filters(["atp,mg"])
        report = pdb_select.rank_candidates(self.candidates(), ligand_filters=filters)
        self.assertEqual(report["ligand_filter"], ["ATP", "MG"])
        self.assertEqual(report["best"]["pdb_id"], "1CCC")
        self.assertEqual(report["eligible_count"], 3)
        self.assertEqual(report["excluded_count"], 1)
        self.assertEqual(report["excluded"][0]["pdb_id"], "2AAA")
        self.assertIn("matched ligand filter: ATP", " ".join(report["best"]["reasons"]))

    def test_normalizes_rcsb_entry_shape_and_validation(self):
        entry = {
            "rcsb_id": "5eee",
            "struct": {"title": "RCSB shaped entry"},
            "rcsb_entry_info": {
                "resolution_combined": [2.1],
                "experimental_method": "ELECTRON MICROSCOPY",
                "nonpolymer_bound_components": ["HEM"],
                "assembly_count": 3,
            },
            "pdbx_vrpt_summary_geometry": [
                {
                    "clashscore": 6.5,
                    "percent_ramachandran_outliers": 0.0,
                }
            ],
        }
        candidate = pdb_select.normalize_candidate(entry)
        self.assertEqual(candidate["pdb_id"], "5EEE")
        self.assertEqual(candidate["title"], "RCSB shaped entry")
        self.assertEqual(candidate["resolution"], 2.1)
        self.assertEqual(candidate["method_category"], "electron microscopy")
        self.assertEqual(candidate["ligands"], ["HEM"])
        self.assertEqual(candidate["assembly_count"], 3)
        self.assertEqual(candidate["validation"]["clashscore"], 6.5)
        self.assertEqual(candidate["validation"]["ramachandran_outliers_percent"], 0.0)

    def test_method_preference_orders_nmr_before_other_and_unknown(self):
        self.assertEqual(pdb_select._method_preference(["SOLUTION NMR"]), ("NMR", 3))
        self.assertEqual(pdb_select._method_preference(["NEUTRON DIFFRACTION"]), ("other", 4))
        self.assertEqual(pdb_select._method_preference([]), ("unknown", 5))

    def test_cli_ranks_offline_repo_style_json(self):
        payload = {"status": "ok", "data": {"results": self.candidates()}}
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
            json.dump(payload, handle)
            path = Path(handle.name)
        try:
            proc = run_script("scripts/pdb_select.py", "--input", str(path), "--ligand", "ATP", "--json")
        finally:
            path.unlink()
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "ok")
        self.assertIn("data", data)
        self.assertEqual(data["best"]["pdb_id"], "1CCC")
        self.assertEqual(data["source"]["kind"], "input_file")
        self.assertEqual(data["source"]["file"], path.name)

    def test_cli_help_includes_examples(self):
        proc = run_script("scripts/pdb_select.py", "--help")
        self.assertIn("usage:", proc.stdout)
        self.assertIn("Examples:", proc.stdout)
        self.assertIn("--live", proc.stdout)


if __name__ == "__main__":
    unittest.main()
