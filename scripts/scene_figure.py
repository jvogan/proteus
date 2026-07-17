#!/usr/bin/env python3
"""Compile a declarative scene manifest into PyMOL or ChimeraX figures."""

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


def _load_manifest(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Manifest not found: {proteus_common.display_path(source)}")
    payload = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise visual_common.VisualWorkflowError("Scene manifest must be a JSON object.")
    return payload


def _structures(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    values = manifest.get("structures")
    if not isinstance(values, list) or not values:
        raise visual_common.VisualWorkflowError("Manifest requires a non-empty structures list.")
    output = []
    for index, value in enumerate(values, start=1):
        if not isinstance(value, dict) or "path" not in value:
            raise visual_common.VisualWorkflowError(f"Structure {index} requires a path.")
        path = Path(str(value["path"])).expanduser()
        if not path.is_file():
            raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(path)}")
        name = visual_common.identifier(str(value.get("id") or f"structure_{index}"), field="structure id")
        output.append({"id": name, "path": path})
    return output


def _views(manifest: dict[str, Any], name: str) -> list[dict[str, Any]]:
    values = manifest.get("views")
    if values is None:
        values = [{"id": "overview", "selection": "all"}]
    if not isinstance(values, list) or not values:
        raise visual_common.VisualWorkflowError("views must be a non-empty list.")
    output = []
    for index, value in enumerate(values, start=1):
        if not isinstance(value, dict):
            raise visual_common.VisualWorkflowError("Each view must be an object.")
        view_id = visual_common.identifier(str(value.get("id") or f"view_{index}"), field="view id")
        output_name = str(value.get("output") or f"{name}_{view_id}.png")
        output_path = Path(output_name)
        if output_path.name != output_name or output_path.suffix.lower() != ".png":
            raise visual_common.VisualWorkflowError("View output must be a PNG filename without directory traversal.")
        turn = value.get("turn") or {}
        if not isinstance(turn, dict):
            raise visual_common.VisualWorkflowError("View turn must be an object with optional x, y, and z angles.")
        output.append({
            "id": view_id,
            "selection": visual_common.selection(str(value.get("selection") or "all")),
            "turn": turn,
            "output": output_name,
        })
    return output


def _representations(manifest: dict[str, Any]) -> list[dict[str, str]]:
    values = manifest.get("representations") or [{"style": "cartoon", "selection": "all", "color": "chain"}]
    if not isinstance(values, list):
        raise visual_common.VisualWorkflowError("representations must be a list.")
    output = []
    for value in values:
        if not isinstance(value, dict):
            raise visual_common.VisualWorkflowError("Each representation must be an object.")
        item = {
            "style": visual_common.style(str(value.get("style") or "cartoon")),
            "selection": visual_common.selection(str(value.get("selection") or "all")),
            "color": str(value.get("color") or "chain"),
        }
        if item["color"] not in {"chain", "element", "plddt", "bfactor", "rainbow"}:
            visual_common.color(item["color"])
        if "transparency" in value:
            number = float(value["transparency"])
            if number < 0 or number > 1:
                raise visual_common.VisualWorkflowError("Transparency must be between 0 and 1.")
            item["transparency"] = str(number)
        output.append(item)
    return output


def _labels(manifest: dict[str, Any]) -> list[dict[str, str]]:
    values = manifest.get("labels") or []
    if not isinstance(values, list):
        raise visual_common.VisualWorkflowError("labels must be a list.")
    output = []
    for value in values:
        if not isinstance(value, dict):
            raise visual_common.VisualWorkflowError("Each label must be an object.")
        text = str(value.get("text") or "%s%s")
        if any(character in text for character in "\r\n;"):
            raise visual_common.VisualWorkflowError("Label text cannot contain newlines or semicolons.")
        output.append({
            "selection": visual_common.selection(str(value.get("selection") or "name CA")),
            "text": text,
            "color": visual_common.color(str(value.get("color") or "black")),
        })
    return output


def _pymol_script(
    structures: list[dict[str, Any]], representations: list[dict[str, str]],
    labels: list[dict[str, str]], views: list[dict[str, Any]], outdir: Path,
    *, background: str, width: int, height: int, session: Path,
) -> list[str]:
    lines = visual_common.pymol_base(background=background, width=width, height=height)
    for item in structures:
        lines.append(f"load {visual_common.quote_pymol(item['path'].resolve())}, {item['id']}")
    lines.append("hide everything, all")
    for item in representations:
        lines.append(f"show {item['style']}, {item['selection']}")
        mode = item["color"]
        if mode == "chain":
            lines.append(f"util.cbc {item['selection']}")
        elif mode == "element":
            lines.append(f"util.cnc {item['selection']}")
        elif mode == "plddt":
            lines.extend([
                f"color orange, {item['selection']}",
                f"color yellow, ({item['selection']}) and b > 50",
                f"color cyan, ({item['selection']}) and b > 70",
                f"color blue, ({item['selection']}) and b > 90",
            ])
        elif mode == "bfactor":
            lines.append(f"spectrum b, blue_white_red, {item['selection']}")
        elif mode == "rainbow":
            lines.append(f"spectrum count, rainbow, {item['selection']}")
        else:
            lines.append(f"color {mode}, {item['selection']}")
        if "transparency" in item:
            setting = "surface_transparency" if item["style"] == "surface" else f"{item['style'].rstrip('s')}_transparency"
            lines.append(f"set {setting}, {item['transparency']}, {item['selection']}")
    for item in labels:
        text = item["text"].replace('"', '\\"')
        lines.extend([
            f"label {item['selection']}, \"{text}\"",
            f"set label_color, {item['color']}, {item['selection']}",
        ])
    for view in views:
        lines.append(f"orient {view['selection']}")
        for axis in ("x", "y", "z"):
            if axis in view["turn"]:
                lines.append(f"turn {axis}, {float(view['turn'][axis])}")
        output = (outdir / view["output"]).resolve()
        lines.append(f"ray {width}, {height}")
        lines.extend(visual_common.pymol_png(output, width=width, height=height))
    lines.extend(visual_common.pymol_save(session))
    lines.append("quit")
    return lines


def _chimerax_script(
    structures: list[dict[str, Any]], representations: list[dict[str, str]],
    labels: list[dict[str, str]], views: list[dict[str, Any]], outdir: Path,
    *, background: str, width: int, height: int, session: Path,
) -> list[str]:
    lines = visual_common.chimerax_base(background=background)
    for item in structures:
        lines.append(f"open {visual_common.quote_chimerax(item['path'].resolve())} name {item['id']}")
    lines.append("hide")
    for item in representations:
        target = item["selection"]
        style = item["style"]
        if style in {"cartoon", "ribbon"}:
            lines.append(f"cartoon {target}")
        elif style == "surface":
            lines.append(f"surface {target}")
        else:
            atom_style = {"sticks": "stick", "spheres": "sphere", "lines": "wire"}[style]
            lines.extend([f"show {target} atoms", f"style {target} {atom_style}"])
        mode = item["color"]
        if mode == "chain":
            lines.append(f"color {target} bychain")
        elif mode == "element":
            lines.append(f"color {target} byelement")
        elif mode == "plddt":
            lines.append(f"color bfactor {target} palette alphafold")
        elif mode == "bfactor":
            lines.append(f"color bfactor {target}")
        elif mode == "rainbow":
            lines.append(f"rainbow {target}")
        else:
            lines.append(f"color {target} {mode}")
        if "transparency" in item:
            target_type = "s" if style == "surface" else "c" if style in {"cartoon", "ribbon"} else "a"
            lines.append(f"transparency {target} {float(item['transparency']) * 100:.1f} target {target_type}")
    for index, item in enumerate(labels, start=1):
        text = item["text"].replace('"', '\\"')
        lines.append(f"label {item['selection']} text \"{text}\" color {item['color']}")
        lines.append(f"2dlabels create legend{index} text \"{text}\" color {item['color']} xpos 0.03 ypos {0.96 - index * 0.05:.2f}")
    for view in views:
        lines.append(f"view {view['selection']}")
        for axis in ("x", "y", "z"):
            if axis in view["turn"]:
                lines.append(f"turn {axis} {float(view['turn'][axis])}")
        output = (outdir / view["output"]).resolve()
        lines.extend([f"save {visual_common.quote_chimerax(output)} width {width} height {height} supersample 3", "wait 1"])
    lines.append(f"save {visual_common.quote_chimerax(session.resolve())}")
    return lines


def compile_scene(manifest_path: str, outdir: str, *, execute: bool = False) -> dict[str, Any]:
    manifest = _load_manifest(manifest_path)
    name = proteus_common.slug(str(manifest.get("name") or "scene"))
    tool = str(manifest.get("tool") or "pymol").lower()
    if tool not in {"pymol", "chimerax"}:
        raise visual_common.VisualWorkflowError("tool must be pymol or chimerax.")
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    structures = _structures(manifest)
    representations = _representations(manifest)
    labels = _labels(manifest)
    views = _views(manifest, name)
    width = int(manifest.get("width") or 1200)
    height = int(manifest.get("height") or 900)
    background = visual_common.color(str(manifest.get("background") or "white"))
    if width < 100 or height < 100 or width > 10000 or height > 10000:
        raise visual_common.VisualWorkflowError("Image dimensions must be between 100 and 10000 pixels.")
    session = destination / f"{name}.{'pse' if tool == 'pymol' else 'cxs'}"
    if tool == "pymol":
        pml = _pymol_script(structures, representations, labels, views, destination,
                            background=background, width=width, height=height, session=session)
        cxc = None
    else:
        pml = None
        cxc = _chimerax_script(structures, representations, labels, views, destination,
                               background=background, width=width, height=height, session=session)
    data: dict[str, Any] = {
        "workflow": "figure",
        "name": name,
        "tool": tool,
        "structures": [{"id": item["id"], "path": proteus_common.display_path(item["path"])} for item in structures],
        "views": [{"id": item["id"], "output": proteus_common.display_path(destination / item["output"])} for item in views],
        "session": proteus_common.display_path(session),
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data,
        provenance={
            "manifest": proteus_common.file_provenance(manifest_path),
            "structures": [proteus_common.file_provenance(item["path"]) for item in structures],
        },
    )
    outputs = visual_common.write_workflow(destination, name, report=report, pymol_lines=pml, chimerax_lines=cxc)
    data["artifacts"] = outputs
    if execute:
        if tool == "pymol":
            result = visual_common.run_pymol(destination / f"{name}.pml")
        else:
            result = visual_common.run_chimerax(cxc or [])
        data["execution"] = result
        if result.get("status") != "ok":
            report["status"] = "error"
            report["error"] = result.get("error", f"{tool} execution failed.")
        else:
            missing = [item["output"] for item in views if not visual_common.verify_nonempty(destination / item["output"])]
            if missing:
                report["status"] = "error"
                report["error"] = f"Expected image outputs were not created: {', '.join(missing)}"
        proteus_common.write_json(destination / f"{name}.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compile a reproducible PyMOL/ChimeraX scene manifest.")
    parser.add_argument("manifest", help="JSON scene manifest")
    parser.add_argument("--outdir", default="proteus_figure", help="Output directory")
    parser.add_argument("--execute", action="store_true", help="Run the generated PyMOL/ChimeraX workflow")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = compile_scene(args.manifest, args.outdir, execute=args.execute)
    except (OSError, ValueError, json.JSONDecodeError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
