#!/usr/bin/env python3
"""Create ChimeraX electrostatic surface workflows with explicit caveats."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common
import visual_common


def build_electrostatics_workflow(
    structure: str,
    outdir: str,
    *,
    selection: str = "protein",
    range_value: float = 10.0,
    offset: float = 1.4,
    execute: bool = False,
    width: int = 1300,
    height: int = 900,
) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    visual_common.selection(selection)
    if range_value <= 0 or offset < 0:
        raise visual_common.VisualWorkflowError("Potential range must be positive and offset cannot be negative.")
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    image = (destination / "electrostatic_surface.png").resolve()
    session = (destination / "electrostatics.cxs").resolve()
    atom_spec = f"#1 & {selection}"
    commands = [
        *visual_common.chimerax_base(),
        f"open {visual_common.quote_chimerax(source.resolve())}",
        f"surface {atom_spec}",
        f"coulombic {atom_spec} offset {offset:g} palette red-white-blue range {-range_value:g},{range_value:g} key true",
        f"transparency {atom_spec} 8 target s",
        f"view {atom_spec}",
        f"save {visual_common.quote_chimerax(image)} width {width} height {height} supersample 3",
        "wait 1",
        f"save {visual_common.quote_chimerax(session)}",
    ]
    external = {
        "pdb2pqr": bool(proteus_common.find_executable("pdb2pqr", "pdb2pqr30")),
        "apbs": bool(proteus_common.find_executable("apbs")),
    }
    data: dict[str, Any] = {
        "workflow": "electrostatics",
        "structure": proteus_common.display_path(source),
        "selection": selection,
        "potential_range_kcal_per_mol_e": [-range_value, range_value],
        "surface_offset_angstrom": offset,
        "method": "chimerax_coulombic",
        "external_poisson_boltzmann_tools": external,
        "image": proteus_common.display_path(image),
        "session": proteus_common.display_path(session),
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data,
        warnings=[
            "Coulombic coloring is a qualitative potential estimate, not a Poisson-Boltzmann calculation.",
            "Protonation, tautomer, charge, dielectric, ionic-strength, and missing-atom choices can change the interpretation.",
            "Inspect ligands, metals, modified residues, termini, and alternate conformers before drawing mechanistic conclusions.",
        ],
        provenance={"structure": proteus_common.file_provenance(source)},
    )
    data["artifacts"] = visual_common.write_workflow(
        destination, "electrostatics_workflow", report=report, chimerax_lines=commands,
    )
    if execute:
        data["execution"] = {"chimerax": visual_common.run_chimerax(commands, timeout=900)}
        if data["execution"]["chimerax"].get("status") not in {"ok", "unavailable"}:
            report["status"] = "error"
            report["error"] = "ChimeraX electrostatic rendering failed."
        proteus_common.write_json(destination / "electrostatics_workflow.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create a ChimeraX Coulombic electrostatic surface and reusable session.")
    parser.add_argument("structure")
    parser.add_argument("--selection", default="protein", help="ChimeraX atom selection (default: protein)")
    parser.add_argument("--range", dest="range_value", type=float, default=10.0, help="Symmetric potential color range")
    parser.add_argument("--offset", type=float, default=1.4, help="Distance outward from the surface for evaluating potential")
    parser.add_argument("--outdir", default="proteus_electrostatics")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = build_electrostatics_workflow(
            args.structure, args.outdir, selection=args.selection, range_value=args.range_value,
            offset=args.offset, execute=args.execute,
        )
    except (OSError, ValueError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
