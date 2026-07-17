#!/usr/bin/env python3
"""Compare structural states and create reproducible PyMOL comparison figures."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import compare_structures
import interaction_report
import proteus_common
import visual_common


def _contact_residues(report: dict[str, Any]) -> set[str]:
    data = report.get("data", report)
    output: set[str] = set()
    for group in data.get("ligand_groups", []):
        for residue in group.get("contact_residues", []):
            label = residue.get("residue_label")
            if label:
                output.add(str(label))
    return output


def _contact_change(reference: Path, mobile: Path, ligand: str | None) -> tuple[dict[str, Any], list[str]]:
    if not ligand:
        return {}, []
    warnings: list[str] = []
    filters = [ligand.upper()]
    reports = []
    for path in (reference, mobile):
        try:
            reports.append(interaction_report.analyze_interactions(str(path), ligand_filters=filters))
        except (OSError, ValueError, interaction_report.InteractionReportError) as exc:
            warnings.append(f"Could not compare {ligand} contacts for {proteus_common.display_path(path)}: {exc}")
            reports.append({"status": "error"})
    if any(item.get("status") != "ok" for item in reports):
        return {}, warnings
    ref = _contact_residues(reports[0])
    mob = _contact_residues(reports[1])
    return {
        "ligand": ligand.upper(),
        "reference_contacts": sorted(ref),
        "mobile_contacts": sorted(mob),
        "gained_in_mobile": sorted(mob - ref),
        "lost_in_mobile": sorted(ref - mob),
        "preserved": sorted(ref & mob),
    }, warnings


def _pymol_lines(reference: Path, mobile: Path, outdir: Path, ligand: str | None,
                  width: int, height: int) -> list[str]:
    overlay = (outdir / "comparison_overlay.png").resolve()
    exploded = (outdir / "comparison_side_by_side.png").resolve()
    session = (outdir / "comparison.pse").resolve()
    lines = visual_common.pymol_base(width=width, height=height)
    lines.extend([
        f"load {visual_common.quote_pymol(reference.resolve())}, reference",
        f"load {visual_common.quote_pymol(mobile.resolve())}, mobile",
        "hide everything, all",
        "show cartoon, reference or mobile",
        "color gray70, reference",
        "color cyan, mobile",
        "cealign reference, mobile",
    ])
    if ligand:
        ligand_selection = f"resn {ligand.upper()}"
        lines.extend([
            f"show sticks, ({ligand_selection})",
            f"color orange, reference and ({ligand_selection})",
            f"color magenta, mobile and ({ligand_selection})",
            f"show sticks, byres ((reference or mobile) within 4 of ({ligand_selection}))",
        ])
    lines.extend([
        "orient reference or mobile",
        f"ray {width}, {height}",
        *visual_common.pymol_png(overlay, width=width, height=height),
        "translate [22, 0, 0], mobile, camera=1",
        "rotate y, 22, mobile, camera=1",
        "orient reference or mobile",
        f"ray {width}, {height}",
        *visual_common.pymol_png(exploded, width=width, height=height),
        *visual_common.pymol_save(session),
        "quit",
    ])
    return lines


def compare_states(reference: str, mobile: str, outdir: str, *, ligand: str | None = None,
                   execute: bool = False, width: int = 1400, height: int = 900) -> dict[str, Any]:
    ref = Path(reference).expanduser()
    mob = Path(mobile).expanduser()
    for path in (ref, mob):
        if not path.is_file():
            raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(path)}")
    if ligand and not ligand.replace("_", "").isalnum():
        raise visual_common.VisualWorkflowError("Ligand code must be alphanumeric.")
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    contact_change, warnings = _contact_change(ref, mob, ligand)
    data: dict[str, Any] = {
        "workflow": "compare",
        "reference": proteus_common.display_path(ref),
        "mobile": proteus_common.display_path(mob),
        "ligand_contact_changes": contact_change or None,
        "images": {
            "overlay": proteus_common.display_path(destination / "comparison_overlay.png"),
            "side_by_side": proteus_common.display_path(destination / "comparison_side_by_side.png"),
        },
        "session": proteus_common.display_path(destination / "comparison.pse"),
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data,
        warnings=warnings,
        provenance={"reference": proteus_common.file_provenance(ref), "mobile": proteus_common.file_provenance(mob)},
    )
    lines = _pymol_lines(ref, mob, destination, ligand, width, height)
    artifacts = visual_common.write_workflow(destination, "state_compare", report=report, pymol_lines=lines)
    data["artifacts"] = artifacts
    if execute:
        metrics = compare_structures.compare(str(ref), str(mob), per_residue=True)
        if metrics.get("status") == "ok":
            metric_data = metrics.get("data", {})
            data["alignment"] = metric_data.get("alignment")
            data["per_residue"] = metric_data.get("per_residue")
        else:
            report.setdefault("warnings", []).append(f"Alignment metrics failed: {metrics.get('error', 'unknown error')}")
        execution = visual_common.run_pymol(destination / "state_compare.pml")
        data["execution"] = execution
        if execution.get("status") != "ok":
            report["status"] = "error"
            report["error"] = execution.get("error", "PyMOL comparison failed.")
        else:
            expected = [destination / "comparison_overlay.png", destination / "comparison_side_by_side.png"]
            if not all(visual_common.verify_nonempty(path) for path in expected):
                report["status"] = "error"
                report["error"] = "PyMOL did not create all expected comparison images."
        proteus_common.write_json(destination / "state_compare.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare two structural states and render overlay/side-by-side figures.")
    parser.add_argument("reference", help="Reference structure")
    parser.add_argument("mobile", help="Mobile structure")
    parser.add_argument("--ligand", help="Optional ligand code for contact-change analysis")
    parser.add_argument("--outdir", default="proteus_compare", help="Output directory")
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument("--execute", action="store_true", help="Run PyMOL and produce figures/metrics")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = compare_states(args.reference, args.mobile, args.outdir, ligand=args.ligand,
                                execute=args.execute, width=args.width, height=args.height)
    except (OSError, ValueError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
