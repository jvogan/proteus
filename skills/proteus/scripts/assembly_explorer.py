#!/usr/bin/env python3
"""Compare asymmetric-unit, biological-assembly, and crystal-neighbor views."""

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

import assembly_report
import proteus_common
import visual_common


PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def explore_assembly(input_value: str, outdir: str, *, assembly: int = 1,
                     execute: bool = False, width: int = 1300, height: int = 900) -> dict[str, Any]:
    if assembly < 1:
        raise visual_common.VisualWorkflowError("Assembly number must be positive.")
    source = Path(input_value).expanduser()
    pdb_id = input_value.upper() if PDB_ID_RE.fullmatch(input_value.strip()) and not source.exists() else None
    metadata = None
    provenance: dict[str, Any]
    warnings = [
        "Biological assemblies are author/depositor interpretations and should be checked against experimental context.",
        "Crystal neighbors show packing contacts and are not automatically biological interfaces.",
    ]
    if pdb_id:
        metadata = assembly_report.fetch_assembly_report(pdb_id, requested_assembly=assembly)
        provenance = {"source": {"kind": "pdb_id", "id": pdb_id, "api": assembly_report.RCSB_ENTRY.format(pdb_id=pdb_id)}}
        open_command = f"open {pdb_id} from pdb"
        pymol_load = f"fetch {pdb_id}, async=0"
    else:
        if not source.is_file():
            raise visual_common.VisualWorkflowError("Input must be a local structure or four-character PDB ID.")
        provenance = {"source": proteus_common.file_provenance(source)}
        open_command = f"open {visual_common.quote_chimerax(source.resolve())}"
        pymol_load = f"load {visual_common.quote_pymol(source.resolve())}, structure"
        warnings.append("Assembly availability could not be verified from a local file alone.")

    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    asu_image = (destination / "asymmetric_unit.png").resolve()
    assembly_image = (destination / f"assembly_{assembly}.png").resolve()
    crystal_image = (destination / "crystal_neighbors.png").resolve()
    cxc = [
        *visual_common.chimerax_base(),
        open_command,
        "cartoon",
        "color #1 bychain",
        "view",
        f"save {visual_common.quote_chimerax(asu_image)} width {width} height {height} supersample 3",
        "wait 1",
        f"sym #1 assembly {assembly} copies true",
        "color all bychain",
        "view",
        f"save {visual_common.quote_chimerax(assembly_image)} width {width} height {height} supersample 3",
        "wait 1",
        f"save {visual_common.quote_chimerax((destination / 'assembly_explorer.cxs').resolve())}",
    ]
    pml = visual_common.pymol_base(width=width, height=height)
    pml.extend([
        pymol_load,
        "hide everything, all",
        "show cartoon, polymer",
        "util.cbc()",
        "symexp crystal_, all, all, 5",
        "color gray80, crystal_*",
        "set cartoon_transparency, 0.35, crystal_*",
        "orient all",
        f"ray {width}, {height}",
        *visual_common.pymol_png(crystal_image, width=width, height=height),
        *visual_common.pymol_save(destination / "crystal_neighbors.pse"),
        "quit",
    ])
    data: dict[str, Any] = {
        "workflow": "assembly_explorer",
        "input": pdb_id or proteus_common.display_path(source),
        "assembly": assembly,
        "assembly_metadata": metadata.get("data", metadata) if metadata else None,
        "images": {
            "asymmetric_unit": proteus_common.display_path(asu_image),
            "biological_assembly": proteus_common.display_path(assembly_image),
            "crystal_neighbors": proteus_common.display_path(crystal_image),
        },
        "executed": execute,
    }
    report = proteus_common.ok_payload(data, warnings=warnings, provenance=provenance)
    data["artifacts"] = visual_common.write_workflow(
        destination, "assembly_explorer", report=report, pymol_lines=pml, chimerax_lines=cxc,
    )
    if execute:
        data["execution"] = {
            "chimerax": visual_common.run_chimerax(cxc),
            "pymol": visual_common.run_pymol(destination / "assembly_explorer.pml"),
        }
        failures = [name for name, result in data["execution"].items() if result.get("status") not in {"ok", "unavailable"}]
        if failures:
            report["status"] = "error"
            report["error"] = f"Assembly rendering failed for: {', '.join(failures)}"
        proteus_common.write_json(destination / "assembly_explorer.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare asymmetric unit, biological assembly, and crystal neighbors.")
    parser.add_argument("input", help="Local PDB/mmCIF or PDB ID")
    parser.add_argument("--assembly", type=int, default=1)
    parser.add_argument("--outdir", default="proteus_assembly")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = explore_assembly(args.input, args.outdir, assembly=args.assembly, execute=args.execute)
    except (OSError, ValueError, assembly_report.AssemblyReportError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
