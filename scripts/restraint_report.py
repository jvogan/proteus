#!/usr/bin/env python3
"""Evaluate residue-pair distance restraints and create visual overlays."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common
import structure_qc
import visual_common


SPEC_RE = re.compile(r"^[A-Za-z0-9_.+?-]+$")


def _records(path: str | Path) -> list[dict[str, Any]]:
    values = visual_common.load_records(path)
    output = []
    for index, value in enumerate(values, start=1):
        try:
            residue1 = str(value["residue1"])
            residue2 = str(value["residue2"])
        except KeyError as exc:
            raise visual_common.VisualWorkflowError(f"Restraint {index} requires residue1 and residue2.") from exc
        minimum = float(value.get("min", value.get("minimum", 0.0)) or 0.0)
        maximum = float(value.get("max", value.get("maximum", 30.0)) or 30.0)
        if minimum < 0 or maximum <= 0 or minimum > maximum:
            raise visual_common.VisualWorkflowError(f"Restraint {index} has invalid min/max bounds.")
        chain1 = str(value.get("chain1") or "")
        chain2 = str(value.get("chain2") or "")
        atom1 = str(value.get("atom1") or "CA")
        atom2 = str(value.get("atom2") or "CA")
        for field, candidate in (
            ("chain1", chain1), ("residue1", residue1), ("atom1", atom1),
            ("chain2", chain2), ("residue2", residue2), ("atom2", atom2),
        ):
            if candidate and not SPEC_RE.fullmatch(candidate):
                raise visual_common.VisualWorkflowError(f"Restraint {index} has an invalid {field} identifier.")
        output.append({
            "id": str(value.get("id") or f"restraint_{index}"),
            "chain1": chain1,
            "residue1": residue1,
            "atom1": atom1,
            "chain2": chain2,
            "residue2": residue2,
            "atom2": atom2,
            "minimum": minimum,
            "maximum": maximum,
            "label": str(value.get("label") or value.get("id") or f"R{index}"),
        })
    return output


def _find_atom(atoms: list[dict[str, Any]], chain: str, residue: str, atom_name: str) -> tuple[dict[str, Any] | None, int]:
    matches = [
        atom for atom in atoms
        if atom["residue_id"] == residue
        and atom["atom_name"].upper() == atom_name.upper()
        and (not chain or atom["chain"] == chain)
    ]
    return (matches[0] if matches else None), len(matches)


def _distance(left: dict[str, Any], right: dict[str, Any]) -> float:
    return math.sqrt(sum((left[axis] - right[axis]) ** 2 for axis in ("x", "y", "z")))


def _pymol_spec(chain: str, residue: str, atom: str) -> str:
    base = f"resi {residue} and name {atom}"
    return f"chain {chain} and {base}" if chain else base


def _chimerax_spec(chain: str, residue: str, atom: str) -> str:
    return f"#1/{chain}:{residue}@{atom}" if chain else f"#1:{residue}@{atom}"


def evaluate_restraints(structure: str, restraints_file: str, outdir: str, *,
                        model: str = "first", execute: bool = False,
                        width: int = 1300, height: int = 900) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    text = structure_qc._read_text(source)
    fmt = structure_qc._format_from_path(source)
    raw, _meta = structure_qc.parse_mmcif(text) if fmt == "mmcif" else structure_qc.parse_pdb(text)
    atoms, selection = structure_qc.select_atoms(raw, model=model, altloc="highest")
    restraints = _records(restraints_file)
    results = []
    warnings: list[str] = []
    for item in restraints:
        left, left_count = _find_atom(atoms, item["chain1"], item["residue1"], item["atom1"])
        right, right_count = _find_atom(atoms, item["chain2"], item["residue2"], item["atom2"])
        result = dict(item)
        if left_count > 1 or right_count > 1:
            warnings.append(f"{item['id']} matched multiple atoms; specify chain IDs to disambiguate.")
        if left is None or right is None:
            result.update({"status": "unresolved", "distance": None})
        else:
            distance = _distance(left, right)
            status = "satisfied" if item["minimum"] <= distance <= item["maximum"] else "too_short" if distance < item["minimum"] else "too_long"
            result.update({"status": status, "distance": round(distance, 3)})
        results.append(result)
    counts = Counter(item["status"] for item in results)
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)

    pml = visual_common.pymol_base(width=width, height=height)
    pml.extend([
        f"load {visual_common.quote_pymol(source.resolve())}, structure",
        "hide everything, all",
        "show cartoon, polymer",
        "color gray80, polymer",
    ])
    cxc = [
        *visual_common.chimerax_base(),
        f"open {visual_common.quote_chimerax(source.resolve())}",
        "cartoon",
        "color protein lightgray",
    ]
    for index, item in enumerate(results, start=1):
        if item["status"] == "unresolved":
            continue
        first_pml = _pymol_spec(item["chain1"], item["residue1"], item["atom1"])
        second_pml = _pymol_spec(item["chain2"], item["residue2"], item["atom2"])
        object_name = f"restraint_{index}"
        color_name = "green" if item["status"] == "satisfied" else "red"
        pml.extend([
            f"show sticks, ({first_pml}) or ({second_pml})",
            f"distance {object_name}, ({first_pml}), ({second_pml})",
            f"set dash_color, {color_name}, {object_name}",
        ])
        first_cx = _chimerax_spec(item["chain1"], item["residue1"], item["atom1"])
        second_cx = _chimerax_spec(item["chain2"], item["residue2"], item["atom2"])
        cxc.extend([
            f"show {first_cx} atoms",
            f"show {second_cx} atoms",
            f"style {first_cx} stick",
            f"style {second_cx} stick",
            f"distance {first_cx} {second_cx} color {color_name}",
        ])
    pml.extend(visual_common.finalize_pymol(
        destination / "restraints.png", width=width, height=height,
        session=destination / "restraints.pse",
    ))
    cxc.extend(visual_common.finalize_chimerax(
        destination / "restraints_chimerax.png", width=width, height=height,
        session=destination / "restraints.cxs",
    ))
    data: dict[str, Any] = {
        "workflow": "restraint_report",
        "structure": proteus_common.display_path(source),
        "model_selection": selection,
        "summary": {"total": len(results), **dict(sorted(counts.items()))},
        "restraints": results,
        "images": {
            "pymol": proteus_common.display_path(destination / "restraints.png"),
            "chimerax": proteus_common.display_path(destination / "restraints_chimerax.png"),
        },
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data,
        warnings=warnings,
        provenance={
            "structure": proteus_common.file_provenance(source),
            "restraints": proteus_common.file_provenance(restraints_file),
        },
    )
    data["artifacts"] = visual_common.write_workflow(
        destination, "restraint_report", report=report, pymol_lines=pml, chimerax_lines=cxc,
    )
    if execute:
        data["execution"] = {
            "pymol": visual_common.run_pymol(destination / "restraint_report.pml"),
            "chimerax": visual_common.run_chimerax(cxc),
        }
        failures = [name for name, result in data["execution"].items() if result.get("status") not in {"ok", "unavailable"}]
        if failures:
            report["status"] = "error"
            report["error"] = f"Restraint rendering failed for: {', '.join(failures)}"
        proteus_common.write_json(destination / "restraint_report.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate residue-pair distance restraints and render overlays.")
    parser.add_argument("structure")
    parser.add_argument("restraints", help="CSV/JSON with chain1,residue1,atom1,chain2,residue2,atom2,min,max")
    parser.add_argument("--model", default="first")
    parser.add_argument("--outdir", default="proteus_restraints")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = evaluate_restraints(args.structure, args.restraints, args.outdir,
                                     model=args.model, execute=args.execute)
    except (OSError, ValueError, json.JSONDecodeError, structure_qc.StructureQCError,
            visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
