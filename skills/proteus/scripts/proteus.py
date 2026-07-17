#!/usr/bin/env python3
"""Unified entry point for Proteus user workflows."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common


COMMANDS = {
    "annotate": ("annotation_figure", "Map residue annotations and scores onto figures"),
    "assembly": ("assembly_explorer", "Compare asymmetric unit, biological assembly, and crystal neighbors"),
    "chemical-site": ("chemical_site", "Inspect a ligand, cofactor, metal, or water network"),
    "compare": ("state_compare", "Create aligned state-comparison figures and contact changes"),
    "cryoem": ("cryoem_workflow", "Inspect and render a map/model pair"),
    "electrostatics": ("electrostatics_workflow", "Create a Coulombic surface workflow"),
    "ensemble": ("ensemble_report", "Analyze multi-model structural variability"),
    "figure": ("scene_figure", "Compile a declarative PyMOL or ChimeraX scene"),
    "interface": ("interface_story", "Analyze and render a protein-protein interface"),
    "pockets": ("pocket_tunnel", "Run optional pocket detection and render candidates"),
    "qc": ("structure_qc", "Preflight coordinate models, conformers, occupancy, and components"),
    "residue": ("residue_story", "Build a residue- or variant-centered structural story"),
    "restraints": ("restraint_report", "Evaluate and render residue-pair distance restraints"),
}


def _help() -> str:
    rows = [
        "Proteus structural-biology workflows",
        "",
        "Usage: python3 scripts/proteus.py <command> [arguments]",
        "",
        "Commands:",
    ]
    width = max(len(command) for command in COMMANDS)
    rows.extend(f"  {command:<{width}}  {description}" for command, (_module, description) in sorted(COMMANDS.items()))
    rows.extend(["", "Run '<command> --help' for workflow-specific options."])
    return "\n".join(rows)


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] in {"-h", "--help", "help"}:
        print(_help())
        return 0
    if arguments[0] in {"-V", "--version"}:
        print(proteus_common.PROTEUS_VERSION)
        return 0
    command = arguments.pop(0)
    target = COMMANDS.get(command)
    if target is None:
        print(f"Unknown workflow: {command}\n\n{_help()}", file=sys.stderr)
        return 2
    module = importlib.import_module(target[0])
    return int(module.main(arguments) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
