#!/usr/bin/env python3
"""Run optional pocket detection and prepare PyMOL/ChimeraX review artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common
import visual_common


def _copy_matches(root: Path, pattern: str, destination: Path) -> list[Path]:
    destination.mkdir(parents=True, exist_ok=True)
    output = []
    for index, source in enumerate(sorted(root.rglob(pattern)), start=1):
        target = destination / f"pocket_{index}{source.suffix.lower()}"
        shutil.copy2(source, target)
        output.append(target)
    return output


def _run_fpocket(executable: str, structure: Path, destination: Path) -> tuple[dict[str, Any], list[Path]]:
    if structure.suffix.lower() not in {".pdb", ".ent"}:
        return {"status": "error", "error": "fpocket execution requires a PDB-format input."}, []
    with tempfile.TemporaryDirectory(prefix="proteus-fpocket-") as temporary:
        workdir = Path(temporary)
        local = workdir / "input.pdb"
        shutil.copy2(structure, local)
        result = proteus_common.run_command([executable, "-f", local.name], cwd=workdir, timeout=900)
        pockets = _copy_matches(workdir, "pocket*_atm.pdb", destination / "pockets") if result["status"] == "ok" else []
    result.pop("stdout", None)
    result.pop("stderr", None)
    result["pockets_found"] = len(pockets)
    return result, pockets


def _run_p2rank(executable: str, structure: Path, destination: Path) -> tuple[dict[str, Any], list[Path], list[dict[str, str]]]:
    workdir = destination / "p2rank"
    workdir.mkdir(parents=True, exist_ok=True)
    result = proteus_common.run_command(
        [executable, "predict", "-f", str(structure.resolve()), "-o", str(workdir.resolve())],
        timeout=1200,
    )
    result.pop("stdout", None)
    result.pop("stderr", None)
    pockets = _copy_matches(workdir, "*.pdb", destination / "pockets") if result["status"] == "ok" else []
    predictions: list[dict[str, str]] = []
    for table in sorted(workdir.rglob("*predictions*.csv")):
        try:
            with table.open(newline="", encoding="utf-8-sig") as handle:
                predictions.extend(dict(row) for row in csv.DictReader(handle))
        except (OSError, csv.Error):
            continue
    result["pocket_files_found"] = len(pockets)
    result["prediction_rows"] = len(predictions)
    return result, pockets, predictions


def _scripts(structure: Path, pockets: list[Path], destination: Path,
             width: int, height: int) -> tuple[list[str], list[str]]:
    pml = visual_common.pymol_base(width=width, height=height)
    pml.extend([
        f"load {visual_common.quote_pymol(structure.resolve())}, structure",
        "hide everything, all",
        "show cartoon, structure and polymer",
        "color gray80, structure and polymer",
    ])
    cxc = [
        *visual_common.chimerax_base(),
        f"open {visual_common.quote_chimerax(structure.resolve())}",
        "cartoon #1",
        "color #1 lightgray",
    ]
    for index, pocket in enumerate(pockets, start=1):
        pml.extend([
            f"load {visual_common.quote_pymol(pocket.resolve())}, pocket_{index}",
            f"show spheres, pocket_{index}",
            f"set sphere_scale, 0.35, pocket_{index}",
            f"color tv_orange, pocket_{index}",
            f"select pocket_{index}_lining, byres (structure and polymer within 4.5 of pocket_{index})",
            f"show sticks, pocket_{index}_lining",
        ])
        cxc.extend([
            f"open {visual_common.quote_chimerax(pocket.resolve())}",
            f"style #{index + 1} sphere",
            f"color #{index + 1} orange",
            f"select zone #{index + 1} 4.5 residues true",
            "show sel atoms",
            "style sel stick",
        ])
    pml.extend(visual_common.finalize_pymol(
        destination / "pockets.png", width=width, height=height,
        session=destination / "pockets.pse",
    ))
    cxc.extend(visual_common.finalize_chimerax(
        destination / "pockets_chimerax.png", width=width, height=height,
        session=destination / "pockets.cxs",
    ))
    return pml, cxc


def pocket_tunnel_workflow(structure: str, outdir: str, *, detector: str = "auto",
                           execute: bool = False, render: bool = False,
                           width: int = 1300, height: int = 900) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    fpocket = proteus_common.find_executable("fpocket")
    p2rank = proteus_common.find_executable("prank", "p2rank")
    tunnel_tools = {
        "caver": bool(proteus_common.find_executable("caver")),
        "hole": bool(proteus_common.find_executable("hole")),
    }
    if detector not in {"auto", "fpocket", "p2rank", "none"}:
        raise visual_common.VisualWorkflowError("Detector must be auto, fpocket, p2rank, or none.")
    selected = detector
    if detector == "auto":
        selected = "fpocket" if fpocket else "p2rank" if p2rank else "none"
    executable = fpocket if selected == "fpocket" else p2rank if selected == "p2rank" else None
    if execute and selected != "none" and not executable:
        raise visual_common.VisualWorkflowError(f"Requested detector is not installed: {selected}")

    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    detector_result: dict[str, Any] | None = None
    pockets: list[Path] = []
    predictions: list[dict[str, str]] = []
    if execute and selected == "fpocket" and executable:
        detector_result, pockets = _run_fpocket(executable, source, destination)
    elif execute and selected == "p2rank" and executable:
        detector_result, pockets, predictions = _run_p2rank(executable, source, destination)

    pml, cxc = _scripts(source, pockets, destination, width, height)
    warnings = [
        "Predicted pockets and tunnels are geometric hypotheses; prioritize sites using experimental, evolutionary, and chemical evidence.",
        "Tool scores are not binding-affinity estimates and are not directly comparable across detectors.",
        "Tunnel results depend on probe radius, starting point, conformational state, and treatment of waters or ligands.",
    ]
    if selected == "none":
        warnings.append("No supported pocket detector was found; install fpocket or P2Rank, or choose a detector explicitly.")
    data: dict[str, Any] = {
        "workflow": "pocket_tunnel",
        "structure": proteus_common.display_path(source),
        "detector": selected,
        "availability": {
            "fpocket": bool(fpocket),
            "p2rank": bool(p2rank),
            "tunnel_tools": tunnel_tools,
        },
        "detector_execution": detector_result,
        "pockets": [proteus_common.display_path(item) for item in pockets],
        "predictions": predictions,
        "render_requested": render,
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data, warnings=warnings, provenance={"structure": proteus_common.file_provenance(source)},
    )
    data["artifacts"] = visual_common.write_workflow(
        destination, "pocket_tunnel", report=report, pymol_lines=pml, chimerax_lines=cxc,
    )
    if execute and detector_result and detector_result.get("status") != "ok":
        report["status"] = "error"
        report["error"] = f"{selected} pocket detection failed."
    if render:
        data["rendering"] = {
            "pymol": visual_common.run_pymol(destination / "pocket_tunnel.pml"),
            "chimerax": visual_common.run_chimerax(cxc),
        }
        failures = [name for name, result in data["rendering"].items() if result.get("status") not in {"ok", "unavailable"}]
        if failures:
            report["status"] = "error"
            report["error"] = f"Pocket rendering failed for: {', '.join(failures)}"
    proteus_common.write_json(destination / "pocket_tunnel.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run fpocket/P2Rank when available and prepare pocket-review sessions.")
    parser.add_argument("structure")
    parser.add_argument("--detector", default="auto", choices=["auto", "fpocket", "p2rank", "none"])
    parser.add_argument("--execute", action="store_true", help="Run the selected pocket detector")
    parser.add_argument("--render", action="store_true", help="Render generated PyMOL and ChimeraX review artifacts")
    parser.add_argument("--outdir", default="proteus_pockets")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = pocket_tunnel_workflow(
            args.structure, args.outdir, detector=args.detector, execute=args.execute, render=args.render,
        )
    except (OSError, ValueError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
