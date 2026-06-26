#!/usr/bin/env python3
import importlib.machinery
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import dock_prep


RECEPTOR_PDB = """\
TITLE     DOCK PREP RECEPTOR
ATOM      1  N   GLY A   1      11.104  13.207   2.100  1.00 90.00           N
ATOM      2  CA  GLY A   1      12.560  13.205   2.100  1.00 90.00           C
HETATM    3  O   HOH A 101      12.000  14.000   2.000  1.00 20.00           O
HETATM    4 ZN    ZN A 201      13.000  15.000   2.500  1.00 20.00          ZN
HETATM    5  C1  HEM A 301      14.000  15.500   3.000  1.00 20.00           C
END
"""


LIGAND_SDF = """\
Dock prep ligand
  Proteus

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
M  END
$$$$
"""


def run_script(*args, check=True):
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


class DockPrepTests(unittest.TestCase):
    def write_file(self, directory: Path, name: str, text: str) -> Path:
        path = directory / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_build_plan_warns_about_common_docking_prep_risks(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            receptor = self.write_file(tmpdir, "receptor.pdb", RECEPTOR_PDB)
            ligand = self.write_file(tmpdir, "ligand.sdf", LIGAND_SDF)

            plan = dock_prep.build_prep_plan(str(receptor), str(ligand))
            warnings = {warning["code"] for warning in plan["warnings"]}

            self.assertEqual(plan["status"], "ok")
            self.assertFalse(plan["executes_tools"])
            self.assertIn("water_removed_by_default", warnings)
            self.assertIn("metal_review_required", warnings)
            self.assertIn("cofactor_review_required", warnings)
            self.assertIn("receptor_missing_hydrogens", warnings)
            self.assertIn("ligand_missing_hydrogens", warnings)
            self.assertIn("ligand_protonation_ambiguous", warnings)
            self.assertIn("ZN", plan["inputs"]["receptor"]["elements"])
            self.assertEqual(plan["inputs"]["receptor"]["path"], "receptor.pdb (absolute path omitted)")
            self.assertEqual(plan["inputs"]["ligand"]["path"], "ligand.sdf (absolute path omitted)")
            self.assertEqual(len(plan["commands"]), 5)
            self.assertTrue(all(command["executes_by_default"] is False for command in plan["commands"]))

    def test_keep_flags_change_cleanup_policy_and_warning_codes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            receptor = self.write_file(tmpdir, "receptor.pdb", RECEPTOR_PDB)
            ligand = self.write_file(tmpdir, "ligand.sdf", LIGAND_SDF)

            plan = dock_prep.build_prep_plan(
                receptor,
                ligand,
                keep_water=True,
                keep_cofactors=True,
            )
            warnings = {warning["code"] for warning in plan["warnings"]}

            self.assertEqual(plan["cleanup"]["water_policy"], "keep")
            self.assertEqual(plan["cleanup"]["cofactor_policy"], "keep")
            self.assertIn("water_kept", warnings)
            self.assertIn("metals_kept", warnings)
            self.assertIn("cofactors_kept", warnings)
            self.assertNotIn("water_removed_by_default", warnings)
            self.assertNotIn("metal_review_required", warnings)
            self.assertNotIn("cofactor_review_required", warnings)

    def test_optional_detection_is_mockable_and_scrubs_paths(self):
        def fake_which(name):
            if name in {"obabel", "prepare_receptor", "prepare_ligand", "vina"}:
                return f"/private/tools/{name}"
            return None

        def fake_find_spec(name):
            if name in {"pdbfixer", "meeko"}:
                return importlib.machinery.ModuleSpec(
                    name,
                    loader=None,
                    origin=f"/private/site-packages/{name}/__init__.py",
                )
            return None

        with patch.object(dock_prep.shutil, "which", side_effect=fake_which), \
                patch.object(dock_prep.importlib.util, "find_spec", side_effect=fake_find_spec):
            tools = dock_prep.detect_optional_tools()

        self.assertTrue(tools["imports"]["pdbfixer"]["ok"])
        self.assertTrue(tools["imports"]["meeko"]["ok"])
        self.assertFalse(tools["imports"]["openbabel"]["ok"])
        self.assertTrue(tools["executables"]["obabel"]["ok"])
        self.assertEqual(tools["executables"]["obabel"]["path"], "obabel (absolute path omitted)")
        self.assertIn("absolute path omitted", tools["imports"]["pdbfixer"]["origin"])
        self.assertNotIn("/private", json.dumps(tools))

    def test_cli_json_outdir_writes_scrubbed_json_and_markdown(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            receptor = self.write_file(tmpdir, "private_receptor.pdb", RECEPTOR_PDB)
            ligand = self.write_file(tmpdir, "private_ligand.sdf", LIGAND_SDF)
            outdir = tmpdir / "out"

            proc = run_script(
                "scripts/dock_prep.py",
                "--receptor",
                str(receptor),
                "--ligand",
                str(ligand),
                "--outdir",
                str(outdir),
                "--json",
            )
            data = json.loads(proc.stdout)
            disk_json = (outdir / "dock_prep_plan.json").read_text(encoding="utf-8")
            disk_markdown = (outdir / "dock_prep_plan.md").read_text(encoding="utf-8")

            self.assertEqual(data["status"], "ok")
            self.assertTrue(data["outputs"]["written"])
            self.assertTrue((outdir / "dock_prep_plan.json").exists())
            self.assertTrue((outdir / "dock_prep_plan.md").exists())
            self.assertIn("Recommended Commands", disk_markdown)
            self.assertIn("prepare_receptor", disk_markdown)
            self.assertNotIn(str(tmpdir), proc.stdout)
            self.assertNotIn(str(tmpdir), disk_json)
            self.assertNotIn(str(tmpdir), disk_markdown)

    def test_cli_markdown_default_and_missing_input_json_error_are_scrubbed(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            receptor = self.write_file(tmpdir, "receptor.pdb", RECEPTOR_PDB)
            ligand = self.write_file(tmpdir, "ligand.sdf", LIGAND_SDF)

            markdown = run_script(
                "scripts/dock_prep.py",
                "--receptor",
                str(receptor),
                "--ligand",
                str(ligand),
            )
            self.assertIn("# Docking Prep Plan", markdown.stdout)
            self.assertIn("no prep or docking tools were executed", markdown.stdout)
            self.assertNotIn(str(tmpdir), markdown.stdout)

            missing = run_script(
                "scripts/dock_prep.py",
                "--receptor",
                str(tmpdir / "missing_receptor.pdb"),
                "--ligand",
                str(ligand),
                "--json",
                check=False,
            )
            error = json.loads(missing.stdout)

            self.assertEqual(missing.returncode, 1)
            self.assertEqual(error["status"], "error")
            self.assertIn("Receptor file not found", error["error"])
            self.assertNotIn(str(tmpdir), missing.stdout)


if __name__ == "__main__":
    unittest.main()
