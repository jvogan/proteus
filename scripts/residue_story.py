#!/usr/bin/env python3
"""Turn a residue or substitution into a local structural evidence story."""

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

import interface_report
import mutation_triage
import pocket_report
import proteus_common
import visual_common


VARIANT_RE = re.compile(r"^(?:[A-Za-z0-9_.-]+\s+)?(?:p\.)?[A-Za-z]{1,3}\d+[A-Za-z]{1,3}$", re.IGNORECASE)
RESIDUE_RE = re.compile(r"^(?:(?P<chain>[A-Za-z0-9]+):)?(?P<resi>-?\d+[A-Za-z]?)$")


def _focus(focus: str, chain: str | None) -> dict[str, Any]:
    text = focus.strip()
    if chain is not None and not re.fullmatch(r"[A-Za-z0-9]+", chain):
        raise visual_common.VisualWorkflowError("Chain identifiers must be alphanumeric.")
    if VARIANT_RE.fullmatch(text):
        match = re.search(r"(\d+)", text)
        if match is None:
            raise visual_common.VisualWorkflowError("Could not parse variant residue number.")
        return {"kind": "variant", "value": text, "resi": match.group(1), "chain": chain}
    match = RESIDUE_RE.fullmatch(text)
    if not match:
        raise visual_common.VisualWorkflowError("Focus must be a substitution such as R175H or residue such as A:175.")
    return {"kind": "residue", "value": text, "resi": match.group("resi"), "chain": chain or match.group("chain")}


def _selection(parsed: dict[str, Any]) -> str:
    base = f"resi {parsed['resi']}"
    if parsed.get("chain"):
        base = f"chain {parsed['chain']} and {base}"
    return base


def _pymol_lines(structure: Path, parsed: dict[str, Any], outdir: Path,
                  cutoff: float, width: int, height: int) -> list[str]:
    focus = _selection(parsed)
    global_png = (outdir / "residue_context.png").resolve()
    local_png = (outdir / "residue_local.png").resolve()
    session = (outdir / "residue_story.pse").resolve()
    lines = visual_common.pymol_base(width=width, height=height)
    lines.extend([
        f"load {visual_common.quote_pymol(structure.resolve())}, structure",
        "hide everything, all",
        "show cartoon, polymer",
        "color gray80, polymer",
        f"select focus, structure and ({focus})",
        f"select neighborhood, byres (structure within {cutoff:.2f} of focus)",
        f"select nearby_ligands, (organic or inorganic) within {cutoff:.2f} of focus",
        "show sticks, focus or neighborhood or nearby_ligands",
        "color magenta, focus",
        "color cyan, neighborhood and not focus",
        "color orange, nearby_ligands",
        "label focus and name CA, resn + resi",
        "orient structure",
        f"ray {width}, {height}",
        *visual_common.pymol_png(global_png, width=width, height=height),
        "orient focus or neighborhood or nearby_ligands",
        "zoom focus or neighborhood or nearby_ligands, 4",
        "distance local_contacts, focus, neighborhood or nearby_ligands, 3.6, 2",
        "set dash_color, gray40, local_contacts",
        f"ray {width}, {height}",
        *visual_common.pymol_png(local_png, width=width, height=height),
        *visual_common.pymol_save(session),
        "quit",
    ])
    return lines


def build_residue_story(structure: str, focus: str, outdir: str, *, chain: str | None = None,
                        cutoff: float = 5.0, pae: str | None = None, execute: bool = False,
                        width: int = 1200, height: int = 900) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    if cutoff <= 0:
        raise visual_common.VisualWorkflowError("Cutoff must be positive.")
    parsed = _focus(focus, chain)
    residue_arg = parsed["value"] if parsed["kind"] == "residue" else None
    variant_arg = parsed["value"] if parsed["kind"] == "variant" else None
    pocket = pocket_report.analyze_pocket(str(source), cutoff, residue=residue_arg, variant=variant_arg)
    interface = interface_report.analyze_interfaces(str(source), cutoff, residue=residue_arg, variant=variant_arg)
    triage = None
    warnings = [
        "Residue numbering is interpreted directly from the coordinate file unless an external mapping is supplied."
    ]
    if parsed["kind"] == "variant":
        try:
            triage = mutation_triage.triage_variants(
                str(source), [parsed["value"]], chain=parsed.get("chain"), cutoff=cutoff, pae_json=pae,
            )
        except (OSError, ValueError, mutation_triage.MutationTriageError) as exc:
            warnings.append(f"Mutation triage was unavailable: {exc}")
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    data: dict[str, Any] = {
        "workflow": "residue_story",
        "structure": proteus_common.display_path(source),
        "focus": parsed,
        "cutoff_angstrom": cutoff,
        "pocket_context": pocket.get("data", {}).get("residue_focus"),
        "interface_context": interface.get("data", {}).get("residue_focus"),
        "mutation_context": triage.get("data") if triage and triage.get("status") == "ok" else None,
        "images": {
            "global": proteus_common.display_path(destination / "residue_context.png"),
            "local": proteus_common.display_path(destination / "residue_local.png"),
        },
        "session": proteus_common.display_path(destination / "residue_story.pse"),
        "executed": execute,
    }
    provenance: dict[str, Any] = {"structure": proteus_common.file_provenance(source)}
    if pae:
        provenance["pae"] = proteus_common.file_provenance(pae)
    report = proteus_common.ok_payload(data, warnings=warnings, provenance=provenance)
    lines = _pymol_lines(source, parsed, destination, cutoff, width, height)
    data["artifacts"] = visual_common.write_workflow(destination, "residue_story", report=report, pymol_lines=lines)
    if execute:
        result = visual_common.run_pymol(destination / "residue_story.pml")
        data["execution"] = result
        if result.get("status") != "ok":
            report["status"] = "error"
            report["error"] = result.get("error", "PyMOL residue workflow failed.")
        elif not all(visual_common.verify_nonempty(destination / name) for name in ("residue_context.png", "residue_local.png")):
            report["status"] = "error"
            report["error"] = "PyMOL did not create both residue figures."
        proteus_common.write_json(destination / "residue_story.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create a residue-centric structural evidence report and PyMOL figures.")
    parser.add_argument("structure", help="Local PDB/mmCIF structure")
    parser.add_argument("focus", help="Substitution (R175H) or residue (A:175)")
    parser.add_argument("--chain", help="Optional chain override")
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument("--pae", help="Optional local PAE JSON for variants")
    parser.add_argument("--outdir", default="proteus_residue")
    parser.add_argument("--execute", action="store_true", help="Run PyMOL and create figures")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = build_residue_story(args.structure, args.focus, args.outdir, chain=args.chain,
                                     cutoff=args.cutoff, pae=args.pae, execute=args.execute)
    except (OSError, ValueError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
