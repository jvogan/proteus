#!/usr/bin/env python3
"""Analyze a multi-model coordinate ensemble and build replayable views."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common
import structure_qc
import visual_common


SPEC_RE = re.compile(r"^[A-Za-z0-9_.+?-]+$")


def _coordinates_by_model(atoms: list[dict[str, Any]]) -> tuple[list[int], list[tuple[Any, ...]], list[list[list[float]]]]:
    models = sorted({int(atom["model"]) for atom in atoms})
    ca_by_model: dict[int, dict[tuple[Any, ...], list[float]]] = defaultdict(dict)
    for atom in atoms:
        if atom["atom_name"] != "CA":
            continue
        key = (atom["chain"], atom["residue_id"], atom["insertion_code"], atom["resname"])
        ca_by_model[int(atom["model"])][key] = [float(atom[axis]) for axis in ("x", "y", "z")]
    if not models:
        return [], [], []
    common = sorted(set.intersection(*(set(ca_by_model[model]) for model in models)))
    coordinates = [[ca_by_model[model][key] for key in common] for model in models]
    return models, common, coordinates


def _center(coordinates: list[list[float]]) -> tuple[list[list[float]], list[float]]:
    centroid = [sum(point[axis] for point in coordinates) / len(coordinates) for axis in range(3)]
    return [[point[axis] - centroid[axis] for axis in range(3)] for point in coordinates], centroid


def _align(coordinates: list[list[list[float]]]) -> tuple[list[list[list[float]]], str]:
    if not coordinates or not coordinates[0]:
        return coordinates, "none"
    try:
        import numpy as np

        reference, _ = _center(coordinates[0])
        reference_array = np.asarray(reference, dtype=float)
        aligned = []
        for model_coordinates in coordinates:
            centered, _ = _center(model_coordinates)
            mobile = np.asarray(centered, dtype=float)
            covariance = mobile.T @ reference_array
            left, _singular, right_t = np.linalg.svd(covariance)
            rotation = left @ right_t
            if np.linalg.det(rotation) < 0:
                left[:, -1] *= -1
                rotation = left @ right_t
            aligned.append((mobile @ rotation).tolist())
        return aligned, "kabsch_ca"
    except ImportError:
        return [_center(model)[0] for model in coordinates], "centroid_translation"


def _rmsf(coordinates: list[list[list[float]]]) -> list[float]:
    count = len(coordinates)
    output = []
    for atom_index in range(len(coordinates[0])):
        mean = [sum(model[atom_index][axis] for model in coordinates) / count for axis in range(3)]
        squared = sum(
            sum((model[atom_index][axis] - mean[axis]) ** 2 for axis in range(3))
            for model in coordinates
        ) / count
        output.append(math.sqrt(squared))
    return output


def analyze_ensemble(structure: str, outdir: str, *, execute: bool = False,
                     width: int = 1300, height: int = 900) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    text = structure_qc._read_text(source)
    fmt = structure_qc._format_from_path(source)
    raw, parser_meta = structure_qc.parse_mmcif(text) if fmt == "mmcif" else structure_qc.parse_pdb(text)
    atoms, selection = structure_qc.select_atoms(raw, model="all", altloc="highest")
    models, keys, coordinates = _coordinates_by_model(atoms)
    if len(models) < 2:
        raise visual_common.VisualWorkflowError("Ensemble analysis requires at least two coordinate models.")
    if len(keys) < 3:
        raise visual_common.VisualWorkflowError("Fewer than three common C-alpha atoms were found across models.")
    aligned, method = _align(coordinates)
    fluctuations = _rmsf(aligned)
    residues = []
    for key, value in zip(keys, fluctuations):
        chain, residue_id, insertion_code, resname = key
        residues.append({
            "chain": chain,
            "residue": f"{residue_id}{insertion_code}",
            "resname": resname,
            "rmsf_angstrom": round(value, 3),
        })
    ranked = sorted(residues, key=lambda item: item["rmsf_angstrom"], reverse=True)
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)

    pml = visual_common.pymol_base(width=width, height=height)
    pml.extend([
        f"load {visual_common.quote_pymol(source.resolve())}, ensemble",
        "hide everything, all",
        "show cartoon, ensemble",
        "intra_fit ensemble and name CA, 1",
    ])
    cxc = [
        *visual_common.chimerax_base(),
        f"open {visual_common.quote_chimerax(source.resolve())} coordsets true",
        "cartoon #1",
        "coordset #1 1",
    ]
    assigned = 0
    for item in residues:
        chain = str(item["chain"])
        residue = str(item["residue"])
        if not SPEC_RE.fullmatch(chain) or not SPEC_RE.fullmatch(residue):
            continue
        value = item["rmsf_angstrom"]
        pml.append(f"alter ensemble and chain {chain} and resi {residue}, b={value:g}")
        cxc.append(f"setattr #1/{chain}:{residue} r proteus_rmsf {value:g} create true type float")
        assigned += 1
    pml.extend([
        "spectrum b, blue_white_red, ensemble",
        "set all_states, on",
        "set cartoon_transparency, 0.35",
        *visual_common.finalize_pymol(
            destination / "ensemble.png", width=width, height=height,
            session=destination / "ensemble.pse",
        ),
    ])
    cxc.extend([
        "color byattribute r:proteus_rmsf #1 palette blue:white:red",
        "coordset #1",
        f"save {visual_common.quote_chimerax((destination / 'ensemble.cxs').resolve())}",
    ])
    data: dict[str, Any] = {
        "workflow": "ensemble",
        "structure": proteus_common.display_path(source),
        "models": models,
        "model_count": len(models),
        "common_ca_count": len(keys),
        "alignment_method": method,
        "rmsf_summary_angstrom": {
            "mean": round(sum(fluctuations) / len(fluctuations), 3),
            "max": round(max(fluctuations), 3),
        },
        "most_variable_residues": ranked[:20],
        "residues": residues,
        "visual_attributes_assigned": assigned,
        "parser": parser_meta,
        "selection": selection,
        "executed": execute,
    }
    warnings = [
        "RMSF reports coordinate dispersion across the supplied models; it is not an experimental uncertainty estimate.",
        "Only C-alpha atoms present in every model are compared.",
    ]
    if method == "centroid_translation":
        warnings.append("NumPy was unavailable, so models were centered but not rotationally superposed for RMSF.")
    report = proteus_common.ok_payload(
        data, warnings=warnings, provenance={"structure": proteus_common.file_provenance(source)},
    )
    data["artifacts"] = visual_common.write_workflow(
        destination, "ensemble_report", report=report, pymol_lines=pml, chimerax_lines=cxc,
    )
    if execute:
        data["execution"] = {
            "pymol": visual_common.run_pymol(destination / "ensemble_report.pml"),
            "chimerax": visual_common.run_chimerax(cxc),
        }
        failures = [name for name, result in data["execution"].items() if result.get("status") not in {"ok", "unavailable"}]
        if failures:
            report["status"] = "error"
            report["error"] = f"Ensemble rendering failed for: {', '.join(failures)}"
        proteus_common.write_json(destination / "ensemble_report.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Analyze multi-model structural variability and build PyMOL/ChimeraX sessions.")
    parser.add_argument("structure", help="Multi-model PDB/mmCIF file")
    parser.add_argument("--outdir", default="proteus_ensemble")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = analyze_ensemble(args.structure, args.outdir, execute=args.execute)
    except (OSError, ValueError, structure_qc.StructureQCError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
