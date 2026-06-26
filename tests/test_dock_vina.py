#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import dock_vina


RECEPTOR_PDBQT = """\
ATOM      1  N   GLY A   1       0.000   0.000   0.000  0.00  0.00    -0.300 N
ATOM      2  CA  GLY A   1       1.400   0.000   0.000  0.00  0.00     0.100 C
END
"""


LIGAND_PDBQT = """\
REMARK  ligand fixture
ROOT
HETATM    1  C1  LIG X   1       0.000   1.000   0.000  0.00  0.00     0.000 C
ENDROOT
TORSDOF 0
END
"""


CONFIG_TEXT = """\
center_x = 1
center_y = 2
center_z = 3
size_x = 20
size_y = 21
size_z = 22
exhaustiveness = 8
"""


VINA_LOG_TEMPLATE = """\
AutoDock Vina v1.2.5
Scoring function : vina
Rigid receptor: {receptor}
Ligand: {ligand}
Exhaustiveness: 8

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -7.500      0.000      0.000
   2       -8.100      1.200      2.400
   3       -6.900      2.000      3.500
"""


def run_script(*args, check=True):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


def fake_detection():
    return {
        "executables": {
            "vina": {"available": False, "executable": "vina", "path": None},
            "smina": {"available": False, "executable": "smina", "path": None},
            "gnina": {"available": False, "executable": "gnina", "path": None},
        },
        "capabilities": {"vina": False, "vina_compatible": False},
    }


class DockVinaTests(unittest.TestCase):
    def write_file(self, directory: Path, name: str, text: str) -> Path:
        path = directory / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_detect_json_shape_without_running_vina(self):
        proc = run_script("scripts/dock_vina.py", "detect", "--json")
        report = json.loads(proc.stdout)

        self.assertEqual(report["status"], "ok")
        self.assertEqual(set(report["executables"]), {"vina", "smina", "gnina"})
        self.assertIn("vina_compatible", report["capabilities"])
        for name, tool in report["executables"].items():
            self.assertEqual(tool["executable"], name)
            self.assertIsInstance(tool["available"], bool)
            self.assertIn("path", tool)
        self.assertIn("No Vina-compatible docking command was executed.", report["notes"])

    def test_plan_is_dry_run_and_scrubs_absolute_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            receptor = self.write_file(tmpdir, "receptor.pdbqt", RECEPTOR_PDBQT)
            ligand = self.write_file(tmpdir, "ligand.pdbqt", LIGAND_PDBQT)
            config = self.write_file(tmpdir, "vina_box.txt", CONFIG_TEXT)
            report = dock_vina.plan_docking(
                receptor,
                ligand,
                config=config,
                out=tmpdir / "poses.pdbqt",
                log=tmpdir / "vina.log",
                exhaustiveness=16,
                num_modes=5,
                detection=fake_detection(),
            )
            rendered = json.dumps(report)

            self.assertEqual(report["status"], "ok")
            self.assertTrue(report["dry_run"])
            self.assertFalse(report["execute"])
            self.assertEqual(report["selected_engine"], "vina")
            self.assertEqual(report["inputs"]["config"]["values"]["center_x"], 1)
            self.assertIn("--config", report["commands"][0]["command"])
            self.assertIn("receptor.pdbqt (absolute path omitted)", rendered)
            self.assertIn("ligand.pdbqt (absolute path omitted)", rendered)
            self.assertIn("vina_box.txt (absolute path omitted)", rendered)
            self.assertIn("No Vina-compatible docking command was executed.", report["warnings"])
            self.assertNotIn(str(tmpdir), rendered)

    def test_plan_cli_accepts_inline_center_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            receptor = self.write_file(tmpdir, "receptor.pdbqt", RECEPTOR_PDBQT)
            ligand = self.write_file(tmpdir, "ligand.pdbqt", LIGAND_PDBQT)

            proc = run_script(
                "scripts/dock_vina.py",
                "plan",
                "--receptor",
                str(receptor),
                "--ligand",
                str(ligand),
                "--center",
                "1",
                "2",
                "3",
                "--size",
                "20",
                "21",
                "22",
                "--json",
            )
            report = json.loads(proc.stdout)
            rendered = json.dumps(report)

            self.assertEqual(report["status"], "ok")
            self.assertEqual(report["box"]["source"], "inline")
            self.assertEqual(report["box"]["center"], {"x": 1.0, "y": 2.0, "z": 3.0})
            self.assertIn("--center_x", report["commands"][0]["command"])
            self.assertNotIn(str(tmpdir), rendered)

    def test_parse_log_extracts_and_ranks_modes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            log_text = VINA_LOG_TEMPLATE.format(
                receptor=tmpdir / "receptor.pdbqt",
                ligand=tmpdir / "ligand.pdbqt",
            )
            log_path = self.write_file(tmpdir, "vina.log", log_text)
            report = dock_vina.parse_vina_log(log_path)
            rendered = json.dumps(report)

            self.assertEqual(report["status"], "ok")
            self.assertEqual(report["mode_count"], 3)
            self.assertEqual(report["best_mode"]["mode"], 2)
            self.assertEqual(report["best_mode"]["affinity_kcal_mol"], -8.1)
            self.assertEqual([item["mode"] for item in report["modes"]], [1, 2, 3])
            self.assertEqual([item["mode"] for item in report["ranks"]], [2, 1, 3])
            self.assertIn("receptor.pdbqt (absolute path omitted)", rendered)
            self.assertNotIn(str(tmpdir), rendered)

    def test_parse_log_cli_errors_are_path_safe(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            proc = run_script(
                "scripts/dock_vina.py",
                "parse-log",
                str(tmpdir / "missing.log"),
                "--json",
                check=False,
            )
            report = json.loads(proc.stdout)

            self.assertEqual(proc.returncode, 1)
            self.assertEqual(report["status"], "error")
            self.assertIn("missing.log", report["error"])
            self.assertNotIn(str(tmpdir), proc.stdout)

    def test_cli_help_includes_expected_subcommands(self):
        proc = run_script("scripts/dock_vina.py", "--help")

        self.assertIn("detect", proc.stdout)
        self.assertIn("plan", proc.stdout)
        self.assertIn("parse-log", proc.stdout)
        self.assertIn("Examples:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
