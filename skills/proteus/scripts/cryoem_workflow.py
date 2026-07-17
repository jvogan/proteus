#!/usr/bin/env python3
"""Build a replayable ChimeraX cryo-EM map/model inspection workflow."""

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
import map_info
import proteus_common
import visual_common


PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
EMDB_ID_RE = re.compile(r"^(?:EMD[-_ ]?)?([0-9]{3,6})$", re.IGNORECASE)


def _analysis_summary(result: dict[str, Any]) -> dict[str, Any]:
    cleaned = proteus_common.scrub_private(result)
    output: dict[str, Any] = {"status": cleaned.get("status", "error")}
    if cleaned.get("error"):
        output["error"] = cleaned["error"]
    info = cleaned.get("info") or cleaned.get("data", {}).get("info") or []
    useful = []
    for line in info:
        text = str(line)
        if any(term in text.lower() for term in ("correlation", "mean", "rms", "minimum", "maximum", "standard deviation", "map value")):
            useful.append(text)
        if len(useful) >= 30:
            break
    if useful:
        output["measurements"] = useful
    return output


def _model_source(value: str) -> tuple[str, dict[str, Any]]:
    path = Path(value).expanduser()
    if path.is_file():
        return f"open {visual_common.quote_chimerax(path.resolve())}", proteus_common.file_provenance(path)
    if PDB_ID_RE.fullmatch(value.strip()):
        pdb_id = value.upper()
        return f"open {pdb_id} from pdb", {"kind": "pdb_id", "id": pdb_id}
    raise visual_common.VisualWorkflowError("Model must be a local PDB/mmCIF file or four-character PDB ID.")


def _map_source(value: str) -> tuple[str, dict[str, Any], dict[str, Any] | None]:
    path = Path(value).expanduser()
    if path.is_file():
        stats = map_info.map_info(str(path)).get("data", {})
        return (
            f"open {visual_common.quote_chimerax(path.resolve())}",
            proteus_common.file_provenance(path),
            {key: stats.get(key) for key in ("dimensions", "dtype", "mean", "sigma", "suggested_levels")},
        )
    match = EMDB_ID_RE.fullmatch(value.strip())
    if match:
        emdb_id = f"EMD-{match.group(1)}"
        return f"open {emdb_id} from emdb", {"kind": "emdb_id", "id": emdb_id}, None
    raise visual_common.VisualWorkflowError("Map must be a local MRC/CCP4 file or EMDB ID.")


def _parse_levels(value: str | None, stats: dict[str, Any] | None) -> list[float]:
    if value:
        try:
            levels = [float(item.strip()) for item in value.split(",") if item.strip()]
        except ValueError as exc:
            raise visual_common.VisualWorkflowError("--levels must be comma-separated numeric map values.") from exc
        if not levels:
            raise visual_common.VisualWorkflowError("--levels did not contain any map values.")
        return levels
    suggested = (stats or {}).get("suggested_levels") or {}
    defaults = [suggested.get(key) for key in ("1.5_sigma", "2.0_sigma", "3.0_sigma")]
    return [float(item) for item in defaults if item is not None]


def build_cryoem_workflow(
    model: str,
    density_map: str,
    outdir: str,
    *,
    resolution: float | None = None,
    levels: str | None = None,
    fit: bool = False,
    difference: bool = False,
    execute: bool = False,
    width: int = 1300,
    height: int = 900,
) -> dict[str, Any]:
    if resolution is not None and resolution <= 0:
        raise visual_common.VisualWorkflowError("Resolution must be positive.")
    if difference and resolution is None:
        raise visual_common.VisualWorkflowError("--difference requires --resolution for the model-derived map.")
    model_open, model_provenance = _model_source(model)
    map_open, map_provenance, stats = _map_source(density_map)
    contour_levels = _parse_levels(levels, stats)
    if not contour_levels:
        raise visual_common.VisualWorkflowError("Remote EMDB maps require explicit absolute --levels values.")

    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    commands = [
        *visual_common.chimerax_base(),
        model_open,
        map_open,
        "cartoon #1",
        "color #1 bychain",
        "volume #2 style surface",
        "volume #2 color #7FB3D5 transparency 55",
    ]
    if fit:
        fit_command = "fitmap #1 inMap #2"
        if resolution is not None:
            fit_command += f" resolution {resolution:g}"
        fit_command += " metric correlation"
        commands.append(fit_command)
    commands.extend([
        "measure mapstats #2",
    ])
    images: list[str] = []
    for index, level in enumerate(contour_levels, start=1):
        image = (destination / f"contour_{index}_{level:g}.png").resolve()
        images.append(proteus_common.display_path(image))
        commands.extend([
            f"volume #2 level {level:g}",
            "view #1,2",
            f"save {visual_common.quote_chimerax(image)} width {width} height {height} supersample 3",
            "wait 1",
        ])
    if difference:
        difference_image = (destination / "difference_map.png").resolve()
        span = max(abs(item) for item in contour_levels)
        commands.extend([
            f"molmap #1 {resolution:g}",
            "volume subtract #2 #3 minRms true onGrid #2",
            "hide #2 models",
            f"volume #4 level {-span:g} color red level {span:g} color blue transparency 35",
            "view #1,4",
            f"save {visual_common.quote_chimerax(difference_image)} width {width} height {height} supersample 3",
            "wait 1",
        ])
        images.append(proteus_common.display_path(difference_image))
    commands.append(f"save {visual_common.quote_chimerax((destination / 'cryoem_workflow.cxs').resolve())}")

    warnings = [
        "Contour levels are absolute map values; compare several levels and retain the map's experimental context.",
        "Rigid-body fitting can overstate agreement and should be validated independently, especially when the initial pose is poor.",
    ]
    if difference:
        warnings.append("The model-derived difference map is qualitative and depends on resolution, sharpening, masking, and map scaling.")
    data: dict[str, Any] = {
        "workflow": "cryoem",
        "model": proteus_common.display_path(model) if Path(model).expanduser().is_file() else model.upper(),
        "map": proteus_common.display_path(density_map) if Path(density_map).expanduser().is_file() else density_map.upper(),
        "map_statistics": stats,
        "contour_levels": contour_levels,
        "resolution_angstrom": resolution,
        "fit_requested": fit,
        "difference_map_requested": difference,
        "images": images,
        "session": proteus_common.display_path(destination / "cryoem_workflow.cxs"),
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data,
        warnings=warnings,
        provenance={"model": model_provenance, "map": map_provenance},
    )
    data["artifacts"] = visual_common.write_workflow(
        destination, "cryoem_workflow", report=report, chimerax_lines=commands,
    )
    if execute:
        analysis = [model_open, map_open]
        if fit:
            analysis.append(next(item for item in commands if item.startswith("fitmap ")))
        analysis.append("measure mapstats #2")
        data["analysis"] = _analysis_summary(
            chimerax_agent.run_chimerax_command_list(analysis, timeout=600)
        )
        data["execution"] = {"chimerax": visual_common.run_chimerax(commands, timeout=1200)}
        if data["execution"]["chimerax"].get("status") not in {"ok", "unavailable"}:
            report["status"] = "error"
            report["error"] = "ChimeraX map rendering failed."
        proteus_common.write_json(destination / "cryoem_workflow.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Inspect a cryo-EM map/model pair and produce contour, fit, and session artifacts.")
    parser.add_argument("model", help="Local structure or PDB ID")
    parser.add_argument("map", help="Local MRC/CCP4 map or EMDB ID")
    parser.add_argument("--resolution", type=float, help="Nominal map resolution in angstroms")
    parser.add_argument("--levels", help="Comma-separated absolute map contour values")
    parser.add_argument("--fit", action="store_true", help="Rigid-body fit the model into the map")
    parser.add_argument("--difference", action="store_true", help="Create a qualitative map-minus-model difference view")
    parser.add_argument("--outdir", default="proteus_cryoem")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = build_cryoem_workflow(
            args.model, args.map, args.outdir, resolution=args.resolution, levels=args.levels,
            fit=args.fit, difference=args.difference, execute=args.execute,
        )
    except (OSError, ValueError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
