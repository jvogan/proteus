#!/usr/bin/env python3
"""Map user-supplied residue annotations onto PyMOL and ChimeraX scenes."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common
import visual_common


SPEC_RE = re.compile(r"^[A-Za-z0-9_.+?-]+$")


def _annotations(path: str | Path) -> list[dict[str, Any]]:
    records = visual_common.load_records(path)
    output = []
    for index, record in enumerate(records, start=1):
        chain = str(record.get("chain") or "").strip()
        residue = str(record.get("residue") or record.get("resi") or "").strip()
        if not residue:
            raise visual_common.VisualWorkflowError(f"Annotation {index} is missing residue/resi.")
        if any(character in residue for character in "\r\n; "):
            raise visual_common.VisualWorkflowError(f"Invalid residue identifier in annotation {index}.")
        if not SPEC_RE.fullmatch(residue) or (chain and not SPEC_RE.fullmatch(chain)):
            raise visual_common.VisualWorkflowError(f"Invalid chain or residue identifier in annotation {index}.")
        score_raw = record.get("score", record.get("value"))
        score = None
        if score_raw not in {None, ""}:
            score = float(score_raw)
            if not math.isfinite(score):
                raise visual_common.VisualWorkflowError(f"Annotation {index} has a non-finite score.")
        label = str(record.get("label") or "").strip()
        if any(character in label for character in "\r\n;"):
            raise visual_common.VisualWorkflowError(f"Annotation {index} label cannot contain newlines or semicolons.")
        item = {
            "chain": chain,
            "residue": residue,
            "score": score,
            "label": label,
            "color": str(record.get("color") or "").strip(),
        }
        if item["color"]:
            visual_common.color(item["color"])
        output.append(item)
    return output


def _pymol_selection(item: dict[str, Any]) -> str:
    base = f"resi {item['residue']}"
    return f"chain {item['chain']} and {base}" if item["chain"] else base


def _chimerax_selection(item: dict[str, Any]) -> str:
    return f"#1/{item['chain']}:{item['residue']}" if item["chain"] else f"#1:{item['residue']}"


def build_annotation_figure(structure: str, annotation_file: str, outdir: str, *,
                            execute: bool = False, surface: bool = False,
                            width: int = 1200, height: int = 900) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    annotations = _annotations(annotation_file)
    if not annotations:
        raise visual_common.VisualWorkflowError("No annotations were supplied.")
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    scores = [item["score"] for item in annotations if item["score"] is not None]
    score_range = [min(scores), max(scores)] if scores else None

    pml = visual_common.pymol_base(width=width, height=height)
    pml.extend([
        f"load {visual_common.quote_pymol(source.resolve())}, structure",
        "hide everything, all",
        "show cartoon, polymer",
        "color gray80, polymer",
    ])
    if surface:
        pml.extend(["show surface, polymer", "set surface_transparency, 0.2, polymer"])
    for index, item in enumerate(annotations, start=1):
        selector = _pymol_selection(item)
        name = f"annotation_{index}"
        pml.append(f"select {name}, structure and ({selector})")
        pml.append(f"show sticks, {name}")
        if item["score"] is not None:
            pml.append(f"alter {name}, b={item['score']:.8g}")
        if item["color"]:
            pml.append(f"color {item['color']}, {name}")
        if item["label"]:
            label = item["label"].replace('"', '\\"')
            pml.append(f"label {name} and name CA, \"{label}\"")
    if scores and not any(item["color"] for item in annotations):
        minimum, maximum = score_range or [0.0, 1.0]
        if minimum == maximum:
            maximum = minimum + 1.0
        pml.append(f"spectrum b, blue_white_red, annotation_*, minimum={minimum:.8g}, maximum={maximum:.8g}")
    pml.extend(visual_common.finalize_pymol(
        destination / "annotations.png", width=width, height=height,
        session=destination / "annotations.pse",
    ))

    cxc = [
        *visual_common.chimerax_base(),
        f"open {visual_common.quote_chimerax(source.resolve())}",
        "cartoon",
        "color protein lightgray",
    ]
    if surface:
        cxc.extend(["surface protein", "transparency protein 20 target s"])
    for item in annotations:
        selector = _chimerax_selection(item)
        if item["score"] is not None:
            cxc.append(f"setattr {selector} r proteus_score {item['score']:.8g} create true type float")
        if item["color"]:
            cxc.append(f"color {selector} {item['color']}")
        cxc.extend([f"show {selector} atoms", f"style {selector} stick"])
        if item["label"]:
            label = item["label"].replace('"', '\\"')
            cxc.append(f"label {selector} text \"{label}\"")
    if scores and not any(item["color"] for item in annotations):
        minimum, maximum = score_range or [0.0, 1.0]
        if minimum == maximum:
            maximum = minimum + 1.0
        cxc.append(f"color byattribute r:proteus_score #1 palette blue-white-red range {minimum:.8g},{maximum:.8g}")
    cxc.extend(visual_common.finalize_chimerax(
        destination / "annotations_chimerax.png", width=width, height=height,
        session=destination / "annotations.cxs",
    ))

    data: dict[str, Any] = {
        "workflow": "annotation_figure",
        "structure": proteus_common.display_path(source),
        "annotation_count": len(annotations),
        "score_range": score_range,
        "surface": surface,
        "images": {
            "pymol": proteus_common.display_path(destination / "annotations.png"),
            "chimerax": proteus_common.display_path(destination / "annotations_chimerax.png"),
        },
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data,
        warnings=["Annotation residue identifiers use coordinate-file numbering."],
        provenance={
            "structure": proteus_common.file_provenance(source),
            "annotations": proteus_common.file_provenance(annotation_file),
        },
    )
    data["artifacts"] = visual_common.write_workflow(
        destination, "annotation_figure", report=report, pymol_lines=pml, chimerax_lines=cxc,
    )
    if execute:
        data["execution"] = {
            "pymol": visual_common.run_pymol(destination / "annotation_figure.pml"),
            "chimerax": visual_common.run_chimerax(cxc),
        }
        failures = [name for name, result in data["execution"].items() if result.get("status") not in {"ok", "unavailable"}]
        if failures:
            report["status"] = "error"
            report["error"] = f"Annotation rendering failed for: {', '.join(failures)}"
        proteus_common.write_json(destination / "annotation_figure.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Color and label structures from CSV/JSON residue annotations.")
    parser.add_argument("structure")
    parser.add_argument("annotations", help="CSV/JSON with chain,residue,score,label,color fields")
    parser.add_argument("--outdir", default="proteus_annotations")
    parser.add_argument("--surface", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = build_annotation_figure(args.structure, args.annotations, args.outdir,
                                         execute=args.execute, surface=args.surface)
    except (OSError, ValueError, json.JSONDecodeError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
