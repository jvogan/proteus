#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import docking_box


PDB_TEXT = """\
TITLE     DOCKING BOX TEST STRUCTURE
ATOM      1  N   GLY A   1      11.104  13.207   2.100  1.00 90.00           N
HETATM    2  O   HOH A  10      12.000  14.000   2.000  1.00 20.00           O
HETATM    3  C1  LIG B   2      14.000  15.000   3.000  1.00 40.00           C
HETATM    4  O1  LIG B   2      14.500  15.250   3.200  1.00 40.00           O
HETATM    5  S   SO4 C  42A     20.000  20.000   2.000  1.00 20.00           S
HETATM    6  O1  SO4 C  42A     20.500  20.200   2.300  1.00 20.00           O
END
"""


NO_LIGAND_PDB = """\
TITLE     NO LIGAND TEST STRUCTURE
ATOM      1  N   GLY A   1      11.104  13.207   2.100  1.00 90.00           N
HETATM    2  O   HOH A  10      12.000  14.000   2.000  1.00 20.00           O
END
"""


def run_script(*args, check=True):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


class DockingBoxTests(unittest.TestCase):
    def write_pdb(self, text=PDB_TEXT):
        handle = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
        with handle:
            handle.write(text)
        self.addCleanup(lambda: Path(handle.name).exists() and Path(handle.name).unlink())
        return Path(handle.name)

    def test_local_pdb_json_computes_box_around_all_ligand_atoms(self):
        path = self.write_pdb()

        proc = run_script("scripts/docking_box.py", str(path), "--padding", "1.5", "--json")
        data = json.loads(proc.stdout)

        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["input"]["kind"], "local_file")
        self.assertEqual(data["input"]["query"], f"{path.name} (absolute path omitted)")
        self.assertNotIn(str(path.parent), data["input"]["query"])
        self.assertEqual(data["hetatm_records"], 5)
        self.assertEqual(data["water_hetatm_records"], 1)
        self.assertEqual(data["ligand_atom_count"], 4)
        self.assertEqual([group["ligand"] for group in data["ligand_groups"]], ["LIG", "SO4"])
        self.assertAlmostEqual(data["box"]["center"]["x"], 17.25)
        self.assertAlmostEqual(data["box"]["center"]["y"], 17.6)
        self.assertAlmostEqual(data["box"]["center"]["z"], 2.6)
        self.assertAlmostEqual(data["box"]["size"]["x"], 9.5)
        self.assertAlmostEqual(data["box"]["size"]["y"], 8.2)
        self.assertAlmostEqual(data["box"]["size"]["z"], 4.2)

    def test_ligand_filter_limits_selected_atoms_and_box(self):
        path = self.write_pdb()

        proc = run_script(
            "scripts/docking_box.py",
            str(path),
            "--ligand",
            "lig",
            "--padding",
            "2",
            "--json",
        )
        data = json.loads(proc.stdout)

        self.assertEqual(data["ligand_filter"], ["LIG"])
        self.assertEqual(data["ligand_group_count"], 1)
        self.assertEqual(data["ligand_components"][0]["ligand"], "LIG")
        self.assertEqual(data["ligand_atom_count"], 2)
        self.assertAlmostEqual(data["box"]["center"]["x"], 14.25)
        self.assertAlmostEqual(data["box"]["center"]["y"], 15.125)
        self.assertAlmostEqual(data["box"]["center"]["z"], 3.1)
        self.assertAlmostEqual(data["box"]["size"]["x"], 4.5)
        self.assertAlmostEqual(data["box"]["size"]["y"], 4.25)
        self.assertAlmostEqual(data["box"]["size"]["z"], 4.2)

    def test_local_mmcif_json_computes_box(self):
        proc = run_script(
            "scripts/docking_box.py",
            "tests/fixtures/tiny.cif",
            "--padding",
            "2",
            "--json",
        )
        data = json.loads(proc.stdout)

        self.assertEqual(data["status"], "ok")
        self.assertEqual(data["format"], "mmcif")
        self.assertEqual(data["ligand_atom_count"], 1)
        self.assertEqual(data["ligand_groups"][0]["ligand"], "LIG")
        self.assertEqual(data["ligand_groups"][0]["residue_label"], "B:2")
        self.assertAlmostEqual(data["box"]["center"]["x"], 14.0)
        self.assertAlmostEqual(data["box"]["center"]["y"], 15.0)
        self.assertAlmostEqual(data["box"]["center"]["z"], 3.0)
        self.assertAlmostEqual(data["box"]["size"]["x"], 4.0)
        self.assertAlmostEqual(data["box"]["size"]["y"], 4.0)
        self.assertAlmostEqual(data["box"]["size"]["z"], 4.0)

    def test_writes_vina_text_config(self):
        path = self.write_pdb()
        config_handle = tempfile.NamedTemporaryFile(suffix=".txt", delete=False)
        config_handle.close()
        config = Path(config_handle.name)
        self.addCleanup(lambda: config.exists() and config.unlink())

        proc = run_script(
            "scripts/docking_box.py",
            str(path),
            "--ligand",
            "LIG",
            "--padding",
            "2",
            "--config-out",
            str(config),
            "--json",
        )
        data = json.loads(proc.stdout)

        self.assertEqual(data["vina_config"]["path"], f"{config.name} (absolute path omitted)")
        self.assertNotIn(str(config.parent), data["vina_config"]["path"])
        self.assertEqual(
            config.read_text(),
            "\n".join([
                "center_x = 14.250",
                "center_y = 15.125",
                "center_z = 3.100",
                "size_x = 4.500",
                "size_y = 4.250",
                "size_z = 4.200",
                "",
            ]),
        )

    def test_optional_tool_detection_shape(self):
        tools = docking_box.detect_optional_tools()

        self.assertEqual(set(tools), {"vina", "obabel", "smina", "gnina", "rdkit"})
        for name in ("vina", "obabel", "smina", "gnina"):
            self.assertIn("ok", tools[name])
            self.assertIn("path", tools[name])
        self.assertIn("ok", tools["rdkit"])
        self.assertEqual(tools["rdkit"]["module"], "rdkit")

    def test_cli_help_includes_expected_options(self):
        proc = run_script("scripts/docking_box.py", "--help")

        self.assertIn("usage:", proc.stdout)
        self.assertIn("Examples:", proc.stdout)
        self.assertIn("--json", proc.stdout)
        self.assertIn("--config-out", proc.stdout)
        self.assertIn("--padding", proc.stdout)

    def test_no_ligands_reports_json_error_without_network(self):
        path = self.write_pdb(NO_LIGAND_PDB)

        proc = run_script("scripts/docking_box.py", str(path), "--json", check=False)
        data = json.loads(proc.stdout)

        self.assertEqual(proc.returncode, 1)
        self.assertEqual(data["status"], "error")
        self.assertIn("No non-water HETATM ligand atoms", data["error"])


if __name__ == "__main__":
    unittest.main()
