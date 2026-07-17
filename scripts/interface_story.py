#!/usr/bin/env python3
"""Create an interface analysis report with ChimeraX and PyMOL artifacts."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import chimerax_agent
import interface_report
import proteus_common
import visual_common


def _analysis_summary(result: dict[str, Any]) -> dict[str, Any]:
    cleaned = proteus_common.scrub_private(result)
    output: dict[str, Any] = {"status": cleaned.get("status", "error")}
    if cleaned.get("error"):
        output["error"] = cleaned["error"]
    info = cleaned.get("info") or cleaned.get("data", {}).get("info") or []
    text = "\n".join(str(item) for item in info)
    buried = re.findall(r"Buried area .*?=\s*([0-9.]+)", text)
    hydrogen = re.findall(r"(\d+) hydrogen bonds found", text)
    contacts = re.findall(r"(?m)^(\d+) contacts$", text)
    clashes = re.findall(r"(?m)^(\d+) clashes$", text)
    output["metrics"] = {
        "buried_area_angstrom2": float(buried[-1]) if buried else None,
        "hydrogen_bonds": int(hydrogen[-1]) if hydrogen else None,
        "contacts": max((int(value) for value in contacts), default=None),
        "clashes": max((int(value) for value in clashes), default=None),
    }
    return output


def _chain(value: str) -> str:
    if not value or not value.replace("_", "").isalnum():
        raise visual_common.VisualWorkflowError(f"Invalid chain identifier: {value!r}")
    return value


def _pymol_lines(structure: Path, chain1: str, chain2: str, outdir: Path,
                  cutoff: float, width: int, height: int) -> list[str]:
    overview = (outdir / "interface_overview.png").resolve()
    closeup = (outdir / "interface_closeup.png").resolve()
    session = (outdir / "interface_story.pse").resolve()
    lines = visual_common.pymol_base(width=width, height=height)
    lines.extend([
        f"load {visual_common.quote_pymol(structure.resolve())}, complex",
        "hide everything, all",
        f"show cartoon, chain {chain1} or chain {chain2}",
        f"color marine, chain {chain1}",
        f"color salmon, chain {chain2}",
        f"select interface1, byres (chain {chain1} within {cutoff:.2f} of chain {chain2})",
        f"select interface2, byres (chain {chain2} within {cutoff:.2f} of chain {chain1})",
        "show sticks, interface1 or interface2",
        "color cyan, interface1",
        "color orange, interface2",
        f"orient chain {chain1} or chain {chain2}",
        f"ray {width}, {height}",
        *visual_common.pymol_png(overview, width=width, height=height),
        "distance interface_polar, interface1, interface2, 3.6, 2",
        "set dash_color, gray30, interface_polar",
        "label (interface1 or interface2) and name CA, resn + resi",
        "orient interface1 or interface2",
        "zoom interface1 or interface2, 4",
        f"ray {width}, {height}",
        *visual_common.pymol_png(closeup, width=width, height=height),
        *visual_common.pymol_save(session),
        "quit",
    ])
    return lines


def _chimerax_lines(structure: Path, chain1: str, chain2: str, outdir: Path,
                    width: int, height: int) -> list[str]:
    image = (outdir / "interface_chimerax.png").resolve()
    session = (outdir / "interface_story.cxs").resolve()
    left = f"(#1/{chain1} & protein)"
    right = f"(#1/{chain2} & protein)"
    return [
        *visual_common.chimerax_base(),
        f"open {visual_common.quote_chimerax(structure.resolve())}",
        "cartoon",
        f"color #1/{chain1} marine",
        f"color #1/{chain2} salmon",
        f"measure buriedarea {left} withAtoms2 {right} listResidues true select true color magenta",
        f"hbonds #1/{chain1} restrict #1/{chain2} reveal true log true",
        f"contacts #1/{chain1} restrict #1/{chain2} reveal true showDist true log true",
        f"clashes #1/{chain1} restrict #1/{chain2} reveal true showDist true log true",
        "show sel atoms",
        "style sel stick",
        "view sel",
        f"save {visual_common.quote_chimerax(image)} width {width} height {height} supersample 3",
        "wait 1",
        f"save {visual_common.quote_chimerax(session)}",
    ]


def build_interface_story(structure: str, chains: str, outdir: str, *, cutoff: float = 5.0,
                          execute: bool = False, width: int = 1300, height: int = 900) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    parts = [_chain(item.strip()) for item in chains.split(",") if item.strip()]
    if len(parts) != 2:
        raise visual_common.VisualWorkflowError("--chains must contain exactly two identifiers, such as A,B.")
    if cutoff <= 0:
        raise visual_common.VisualWorkflowError("Cutoff must be positive.")
    local = interface_report.analyze_interfaces(str(source), cutoff, chains_filter=parts)
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    data: dict[str, Any] = {
        "workflow": "interface_story",
        "structure": proteus_common.display_path(source),
        "chains": parts,
        "cutoff_angstrom": cutoff,
        "local_geometry": local.get("data", {}),
        "images": {
            "overview": proteus_common.display_path(destination / "interface_overview.png"),
            "closeup": proteus_common.display_path(destination / "interface_closeup.png"),
            "chimerax": proteus_common.display_path(destination / "interface_chimerax.png"),
        },
        "sessions": {
            "pymol": proteus_common.display_path(destination / "interface_story.pse"),
            "chimerax": proteus_common.display_path(destination / "interface_story.cxs"),
        },
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data,
        warnings=[
            "Distance-based contacts are geometric candidates; hydrogen bonding and energetic hotspots require chemical context.",
            "Buried area excludes solvent, ions, and ligands by explicitly measuring protein atoms only.",
        ],
        provenance={"structure": proteus_common.file_provenance(source)},
    )
    pml = _pymol_lines(source, parts[0], parts[1], destination, cutoff, width, height)
    cxc = _chimerax_lines(source, parts[0], parts[1], destination, width, height)
    data["artifacts"] = visual_common.write_workflow(
        destination, "interface_story", report=report, pymol_lines=pml, chimerax_lines=cxc,
    )
    if execute:
        analysis_commands = [
            f"open {visual_common.quote_chimerax(source.resolve())}",
            f"measure buriedarea (#1/{parts[0]} & protein) withAtoms2 (#1/{parts[1]} & protein) listResidues true",
            f"hbonds #1/{parts[0]} restrict #1/{parts[1]} log true",
            f"contacts #1/{parts[0]} restrict #1/{parts[1]} log true",
            f"clashes #1/{parts[0]} restrict #1/{parts[1]} log true",
        ]
        data["chimerax_analysis"] = _analysis_summary(
            chimerax_agent.run_chimerax_command_list(analysis_commands, timeout=300)
        )
        pymol_result = visual_common.run_pymol(destination / "interface_story.pml")
        chimerax_result = visual_common.run_chimerax(cxc)
        data["execution"] = {"pymol": pymol_result, "chimerax": chimerax_result}
        failures = [name for name, result in data["execution"].items() if result.get("status") not in {"ok", "unavailable"}]
        if failures:
            report["status"] = "error"
            report["error"] = f"Interface rendering failed for: {', '.join(failures)}"
        proteus_common.write_json(destination / "interface_story.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create a protein-interface analysis and publication figure set.")
    parser.add_argument("structure", help="Local PDB/mmCIF structure")
    parser.add_argument("--chains", required=True, help="Two chain IDs, e.g. A,D")
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument("--outdir", default="proteus_interface")
    parser.add_argument("--execute", action="store_true", help="Run ChimeraX analysis and both renderers")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = build_interface_story(args.structure, args.chains, args.outdir, cutoff=args.cutoff,
                                       execute=args.execute)
    except (OSError, ValueError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
