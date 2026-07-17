#!/usr/bin/env python3
import json
import struct
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import annotation_figure
import chemical_site
import cryoem_workflow
import electrostatics_workflow
import ensemble_report
import model_quality
import pae_report
import proteus_common
import restraint_report
import scene_figure
import structure_info
import structure_qc
import sync_skill_package
import visual_common


TINY_PDB = ROOT / "tests" / "fixtures" / "tiny.pdb"


class RuntimeAndQCWorkflowTests(unittest.TestCase):
    def test_runtime_envelope_paths_and_secret_urls(self):
        report = proteus_common.ok_payload({"value": 1}, warnings=["check"], provenance={"source": "fixture"})
        self.assertEqual(report["status"], "ok")
        self.assertEqual(report["schema_version"], 1)
        self.assertNotIn(str(Path.home()), proteus_common.display_path(TINY_PDB))
        self.assertNotIn(str(Path.home()), proteus_common.scrub_text(str(Path.home()).encode()))
        with self.assertRaises(proteus_common.ProteusRuntimeError):
            proteus_common.request_bytes("https://example.org/data?token=secret", retries=0)

    def test_structure_qc_reports_selection_components_and_missing_backbone(self):
        report = structure_qc.build_report(TINY_PDB)
        self.assertEqual(report["status"], "ok")
        data = report["data"]
        self.assertEqual(data["selection"]["selected_models"], [1])
        self.assertEqual(data["component_role_counts"], {"ligand_or_unknown": 1})
        self.assertEqual(data["missing_backbone"][0]["missing"], ["O"])
        self.assertEqual(report["provenance"]["input"]["sha256"], proteus_common.sha256_file(TINY_PDB))

    def test_structure_info_uses_ca_per_residue_for_plddt(self):
        report = structure_info.inspect_structure(str(TINY_PDB), force_alphafold=True)
        self.assertEqual(report["plddt"]["basis"], "ca_per_residue")
        self.assertEqual(report["plddt"]["residues"][0]["plddt"], 91.0)

    def test_chain_aware_pae_blocks(self):
        report = pae_report.summarize_pae(
            str(ROOT / "tests" / "fixtures" / "tiny_pae.json"), min_segment=1,
            chain_lengths="A:2,B:2",
        )
        blocks = report["data"]["chain_analysis"]["blocks"]
        self.assertEqual(len(blocks), 4)
        self.assertLess(blocks[0]["mean_pae"], blocks[1]["mean_pae"])


class ReproducibleWorkflowTests(unittest.TestCase):
    def test_scene_manifest_writes_scripts_and_persisted_artifact_index(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "scene.json"
            manifest.write_text(json.dumps({
                "name": "tiny-scene",
                "tool": "pymol",
                "structures": [{"id": "tiny", "path": str(TINY_PDB)}],
                "representations": [{"style": "cartoon", "selection": "polymer", "color": "chain"}],
                "views": [{"id": "overview", "output": "overview.png"}],
            }), encoding="utf-8")
            report = scene_figure.compile_scene(str(manifest), str(root / "out"))
            disk = json.loads((root / "out" / "tiny-scene.json").read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "ok")
            self.assertTrue((root / "out" / "tiny-scene.pml").is_file())
            self.assertIn("artifacts", disk["data"])

    def test_scene_manifest_rejects_output_traversal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "scene.json"
            manifest.write_text(json.dumps({
                "structures": [{"path": str(TINY_PDB)}],
                "views": [{"output": "../outside.png"}],
            }), encoding="utf-8")
            with self.assertRaises(visual_common.VisualWorkflowError):
                scene_figure.compile_scene(str(manifest), str(root / "out"))

    def test_chimerax_scene_uses_supported_label_and_transparency_syntax(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "scene.json"
            manifest.write_text(json.dumps({
                "tool": "chimerax",
                "structures": [{"path": str(TINY_PDB)}],
                "representations": [{
                    "style": "cartoon", "selection": "#1", "color": "chain", "transparency": 0.25,
                }],
                "labels": [{"selection": "#1/A:1@CA", "text": "site"}],
            }), encoding="utf-8")
            scene_figure.compile_scene(str(manifest), str(root / "out"))
            cxc = (root / "out" / "scene.cxc").read_text(encoding="utf-8")
            self.assertIn("transparency #1 25.0 target c", cxc)
            self.assertIn("2dlabels create legend1 text \"site\"", cxc)
            self.assertNotIn("create name", cxc)

    def test_annotation_and_restraint_reports(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            annotations = root / "annotations.csv"
            annotations.write_text("chain,residue,score,label\nA,1,0.8,site\n", encoding="utf-8")
            annotation = annotation_figure.build_annotation_figure(
                str(TINY_PDB), str(annotations), str(root / "annotation"),
            )
            self.assertEqual(annotation["data"]["annotation_count"], 1)

            restraints = root / "restraints.json"
            restraints.write_text(json.dumps([{
                "chain1": "A", "residue1": "1", "atom1": "CA",
                "chain2": "B", "residue2": "2", "atom2": "C1",
                "min": 0, "max": 5,
            }]), encoding="utf-8")
            result = restraint_report.evaluate_restraints(
                str(TINY_PDB), str(restraints), str(root / "restraints"),
            )
            self.assertEqual(result["data"]["summary"]["satisfied"], 1)
            self.assertGreater(result["data"]["restraints"][0]["distance"], 0)
            cxc = (root / "restraints" / "restraint_report.cxc").read_text(encoding="utf-8")
            self.assertIn("show #1/A:1@CA atoms\nshow #1/B:2@C1 atoms", cxc)
            self.assertNotIn("@CA,#1", cxc)

    def test_chemical_site_reports_neighborhood(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = chemical_site.inspect_chemical_site(
                str(TINY_PDB), temporary, component="LIG", radius=5,
            )
            self.assertEqual(report["status"], "ok")
            self.assertEqual(report["data"]["site"]["role"], "ligand_or_unknown")
            self.assertTrue(report["data"]["neighbors"])

    def test_cryoem_local_map_generates_contour_sweep(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            density = root / "tiny.mrc"
            header = bytearray(1024)
            struct.pack_into("<4i", header, 0, 2, 2, 1, 2)
            struct.pack_into("<i", header, 92, 0)
            density.write_bytes(bytes(header) + struct.pack("<4f", 0.0, 1.0, 2.0, 3.0))
            report = cryoem_workflow.build_cryoem_workflow(
                str(TINY_PDB), str(density), str(root / "out"), resolution=3.0,
            )
            self.assertEqual(report["status"], "ok")
            self.assertEqual(len(report["data"]["contour_levels"]), 3)
            self.assertTrue((root / "out" / "cryoem_workflow.cxc").is_file())

    def test_electrostatics_uses_supported_coulombic_options(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = electrostatics_workflow.build_electrostatics_workflow(
                str(TINY_PDB), str(root / "out"), offset=1.8,
            )
            self.assertEqual(report["status"], "ok")
            cxc = (root / "out" / "electrostatics_workflow.cxc").read_text(encoding="utf-8")
            self.assertIn("coulombic #1 & protein offset 1.8", cxc)
            self.assertIn("transparency #1 & protein 8 target s", cxc)
            self.assertNotIn(" distance ", cxc)

    def test_execution_summary_omits_commands_and_scrubs_local_paths(self):
        summary = visual_common.execution_summary({
            "status": "ok",
            "returncode": 0,
            "stdout": f"opened {Path.home()}/private/input.pdb",
            "payload": {"data": {"history": [
                {"command": f"open {Path.home()}/private/input.pdb", "elapsed_seconds": 0.25},
            ]}},
        })
        self.assertEqual(summary, {"status": "ok", "returncode": 0, "commands_completed": 1, "elapsed_seconds": 0.25})
        self.assertNotIn(str(Path.home()), json.dumps(summary))

    def test_ensemble_rmsf_uses_common_ca_atoms(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            structure = root / "ensemble.pdb"
            lines = []
            serial = 1
            for model, offset in ((1, 0.0), (2, 0.4)):
                lines.append(f"MODEL     {model:4d}")
                for residue, x in ((1, 0.0), (2, 3.8), (3, 7.6)):
                    lines.append(
                        f"ATOM  {serial:5d}  CA  ALA A{residue:4d}    "
                        f"{x:8.3f}{offset * residue:8.3f}{0.0:8.3f}{1.0:6.2f}{80.0:6.2f}           C"
                    )
                    serial += 1
                lines.append("ENDMDL")
            structure.write_text("\n".join(lines) + "\n", encoding="utf-8")
            report = ensemble_report.analyze_ensemble(str(structure), str(root / "out"))
            self.assertEqual(report["data"]["model_count"], 2)
            self.assertEqual(report["data"]["common_ca_count"], 3)
            self.assertEqual(len(report["data"]["residues"]), 3)


class OptionalToolAndPackagingTests(unittest.TestCase):
    def test_foldseek_table_parser(self):
        rows = model_quality._parse_foldseek_table("q\tt\t0.8\t0.7\t0.9\t0.75\t1e-6\t120\n")
        self.assertEqual(rows[0]["query"], "q")
        self.assertEqual(rows[0]["qtmscore"], 0.7)
        self.assertEqual(rows[0]["evalue"], 1e-6)

    def test_installable_package_is_in_sync(self):
        self.assertEqual(sync_skill_package.compare()["status"], "ok")


if __name__ == "__main__":
    unittest.main()
