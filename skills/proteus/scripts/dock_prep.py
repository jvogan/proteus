#!/usr/bin/env python3
"""Plan receptor and ligand preparation for local docking workflows.

The planner is intentionally stdlib-only. It validates local inputs, performs
lightweight structural checks, detects common optional docking-prep tools, and
emits commands for the user to review. It never executes chemistry tools.

Usage:
    python3 dock_prep.py --receptor receptor.pdb --ligand ligand.sdf
    python3 dock_prep.py --receptor receptor.pdb --ligand ligand.sdf --json
    python3 dock_prep.py --receptor receptor.pdb --ligand ligand.sdf --outdir prep-plan
"""

import argparse
import importlib.util
import json
import re
import shlex
import shutil
import sys
from pathlib import Path
from typing import Any


WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}
METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS",
    "BE", "MG", "CA", "SR", "BA",
    "AL", "GA", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG",
}
TWO_LETTER_ELEMENTS = METAL_ELEMENTS | {
    "AC", "AG", "AM", "AR", "AS", "AT", "AU", "BI", "BK", "BR", "CE", "CF",
    "CL", "CM", "CR", "DY", "ER", "ES", "EU", "FM", "FR", "GD", "HE", "HF",
    "HO", "IN", "IR", "KR", "LA", "LR", "LU", "MD", "MO", "ND", "NE", "NO",
    "NP", "OS", "PA", "PB", "PD", "PM", "PO", "PR", "PT", "PU", "RA", "RE",
    "RH", "RN", "RU", "SB", "SC", "SE", "SI", "SM", "SN", "TA", "TB", "TC",
    "TE", "TH", "TI", "TL", "TM", "XE", "YB", "ZR",
}
KNOWN_COFACTORS = {
    "ADP", "AMP", "ATP", "B12", "COA", "FAD", "FMN", "GDP", "GMP", "GTP",
    "HEA", "HEC", "HEM", "NAD", "NAP", "NDP", "PLP", "SAM", "SAH",
}
TITRATABLE_ELEMENTS = {"N", "O", "S", "P"}

IMPORT_MODULES = ("pdbfixer", "openbabel", "meeko")
EXECUTABLES = ("obabel", "openbabel", "prepare_receptor", "prepare_ligand", "vina")
RECEPTOR_EXTENSIONS = {".pdb", ".ent", ".cif", ".mmcif", ".pdbqt"}
LIGAND_EXTENSIONS = {".sdf", ".mol", ".mol2", ".pdb", ".ent", ".pdbqt"}


class DockPrepError(RuntimeError):
    pass


def _ok_payload(data: dict[str, Any]) -> dict[str, Any]:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _error_payload(message: str) -> dict[str, str]:
    return {"status": "error", "error": _scrub_text(message)}


def _path_display(path: str | Path) -> str:
    value = Path(path)
    if value.is_absolute():
        return f"{value.name} (absolute path omitted)"
    return str(value)


def _scrub_text(text: str) -> str:
    """Remove Unix-style absolute paths from public output."""

    def replace(match: re.Match[str]) -> str:
        path = Path(match.group(0))
        name = path.name or "path"
        return f"{name} (absolute path omitted)"

    return re.sub(r"(?<![\w.-])/(?:[^\s'\",)]+/?)+", replace, text)


def _scrub_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _scrub_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_scrub_value(item) for item in value]
    if isinstance(value, tuple):
        return [_scrub_value(item) for item in value]
    if isinstance(value, str):
        return _scrub_text(value)
    return value


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return slug or "input"


def _quote(value: str) -> str:
    return shlex.quote(value)


def _warning(code: str, severity: str, source: str, message: str,
             recommendation: str) -> dict[str, str]:
    return {
        "code": code,
        "severity": severity,
        "source": source,
        "message": message,
        "recommendation": recommendation,
    }


def _normalize_element(value: str, prefer_two_letter: bool = False) -> str:
    letters = re.sub(r"[^A-Za-z]", "", value).upper()
    if not letters:
        return ""
    if prefer_two_letter and letters[:2] in TWO_LETTER_ELEMENTS:
        return letters[:2]
    return letters[:1]


def _infer_element(atom_name: str, element_field: str = "",
                   prefer_two_letter: bool = False) -> str:
    if element_field:
        return _normalize_element(element_field, prefer_two_letter=True)
    return _normalize_element(atom_name, prefer_two_letter=prefer_two_letter)


def _residue_label(chain: str, residue_id: str, insertion_code: str) -> str:
    suffix = insertion_code if insertion_code else ""
    return f"{chain}:{residue_id}{suffix}"


def _sort_residue_label(label: str) -> tuple[str, int | str]:
    chain, _, residue = label.partition(":")
    number = residue.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
    try:
        return chain, int(number)
    except ValueError:
        return chain, residue


def _empty_inspection(fmt: str, inspection_available: bool = False) -> dict[str, Any]:
    return {
        "format": fmt,
        "inspection_available": inspection_available,
        "atom_records": 0,
        "hetatm_records": 0,
        "atom_count": 0,
        "hydrogen_atoms": 0,
        "water_records": 0,
        "water_residues": [],
        "metal_atoms": 0,
        "metal_residues": [],
        "cofactor_records": 0,
        "cofactor_residues": [],
        "hetero_residues": [],
        "elements": [],
        "titratable_elements": [],
    }


def _parse_pdb_like(path: Path, fmt: str) -> dict[str, Any]:
    summary = _empty_inspection(fmt, inspection_available=True)
    water_residues: set[str] = set()
    metal_residues: set[str] = set()
    cofactor_residues: set[str] = set()
    hetero_residues: set[str] = set()
    elements: set[str] = set()

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            summary["atom_count"] += 1
            if record == "ATOM":
                summary["atom_records"] += 1
            else:
                summary["hetatm_records"] += 1

            atom_name = line[12:16].strip() if len(line) >= 16 else ""
            resname = line[17:20].strip().upper() if len(line) >= 20 else ""
            chain = line[21].strip() if len(line) >= 22 else "?"
            residue_id = line[22:26].strip() if len(line) >= 26 else "?"
            insertion_code = line[26].strip() if len(line) >= 27 else ""
            element_field = line[76:78].strip() if len(line) >= 78 else ""
            element = _infer_element(atom_name, element_field)
            label = _residue_label(chain or "?", residue_id or "?", insertion_code)

            if element:
                elements.add(element)
            if element == "H":
                summary["hydrogen_atoms"] += 1
            if element in TITRATABLE_ELEMENTS:
                summary["titratable_elements"].append(element)

            if resname in WATER_NAMES:
                summary["water_records"] += 1
                water_residues.add(label)
                continue

            is_metal = element in METAL_ELEMENTS or resname in METAL_ELEMENTS
            if is_metal:
                summary["metal_atoms"] += 1
                metal_residues.add(f"{resname or element} {label}")

            if record == "HETATM" and not is_metal:
                hetero_label = f"{resname or 'UNK'} {label}"
                hetero_residues.add(hetero_label)
                if resname in KNOWN_COFACTORS or resname:
                    summary["cofactor_records"] += 1
                    cofactor_residues.add(hetero_label)

    summary["water_residues"] = sorted(water_residues, key=_sort_residue_label)
    summary["metal_residues"] = sorted(metal_residues)
    summary["cofactor_residues"] = sorted(cofactor_residues)
    summary["hetero_residues"] = sorted(hetero_residues)
    summary["elements"] = sorted(elements)
    summary["titratable_elements"] = sorted(set(summary["titratable_elements"]))
    return summary


def _parse_sdf(path: Path) -> dict[str, Any]:
    summary = _empty_inspection("sdf", inspection_available=True)
    elements: set[str] = set()
    atom_lines: list[str] = []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) >= 4:
        fields = lines[3].split()
        if fields and fields[0].isdigit():
            atom_count = int(fields[0])
            atom_lines = lines[4:4 + atom_count]

    for line in atom_lines:
        fields = line.split()
        if len(fields) < 4:
            continue
        element = _infer_element(fields[3], prefer_two_letter=True)
        if not element:
            continue
        summary["atom_count"] += 1
        elements.add(element)
        if element == "H":
            summary["hydrogen_atoms"] += 1
        if element in TITRATABLE_ELEMENTS:
            summary["titratable_elements"].append(element)

    summary["elements"] = sorted(elements)
    summary["titratable_elements"] = sorted(set(summary["titratable_elements"]))
    return summary


def _parse_mol2(path: Path) -> dict[str, Any]:
    summary = _empty_inspection("mol2", inspection_available=True)
    elements: set[str] = set()
    in_atom_section = False

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.upper().startswith("@<TRIPOS>"):
                in_atom_section = line.upper() == "@<TRIPOS>ATOM"
                continue
            if not in_atom_section:
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            atom_name = fields[1]
            atom_type = fields[5].split(".", 1)[0]
            element = _infer_element(atom_name, atom_type)
            if not element:
                continue
            summary["atom_count"] += 1
            elements.add(element)
            if element == "H":
                summary["hydrogen_atoms"] += 1
            if element in TITRATABLE_ELEMENTS:
                summary["titratable_elements"].append(element)

    summary["elements"] = sorted(elements)
    summary["titratable_elements"] = sorted(set(summary["titratable_elements"]))
    return summary


def _format_for_path(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".ent":
        return "pdb"
    if suffix == ".mmcif":
        return "cif"
    if suffix:
        return suffix.lstrip(".")
    return "unknown"


def inspect_structure(path: str | Path) -> dict[str, Any]:
    file_path = Path(path)
    fmt = _format_for_path(file_path)
    suffix = file_path.suffix.lower()
    if suffix in {".pdb", ".ent", ".pdbqt"}:
        return _parse_pdb_like(file_path, fmt)
    if suffix in {".sdf", ".mol"}:
        return _parse_sdf(file_path)
    if suffix == ".mol2":
        return _parse_mol2(file_path)
    return _empty_inspection(fmt, inspection_available=False)


def validate_input(path: str | Path, role: str) -> tuple[dict[str, Any], list[dict[str, str]]]:
    file_path = Path(path)
    display = _path_display(file_path)
    role_title = role.title()
    if not file_path.exists():
        raise DockPrepError(f"{role_title} file not found: {display}")
    if not file_path.is_file():
        raise DockPrepError(f"{role_title} input is not a file: {display}")
    try:
        size_bytes = file_path.stat().st_size
    except OSError as exc:
        raise DockPrepError(f"Could not stat {role} input: {_scrub_text(str(exc))}") from exc
    if size_bytes == 0:
        raise DockPrepError(f"{role_title} file is empty: {display}")

    suffix = file_path.suffix.lower()
    expected = RECEPTOR_EXTENSIONS if role == "receptor" else LIGAND_EXTENSIONS
    warnings: list[dict[str, str]] = []
    if suffix not in expected:
        warnings.append(_warning(
            "unsupported_extension",
            "warning",
            role,
            f"{role_title} extension '{suffix or '(none)'}' is not one of the common docking-prep inputs.",
            "Proceed only if a local prep tool can read this format or convert it first.",
        ))

    inspection = inspect_structure(file_path)
    if not inspection["inspection_available"]:
        warnings.append(_warning(
            "limited_format_inspection",
            "info",
            role,
            f"{role_title} format '{inspection['format']}' was validated as a file but not structurally inspected.",
            "Use a format-aware checker before docking prep if waters, metals, or hydrogens are uncertain.",
        ))

    data = {
        "role": role,
        "path": display,
        "name": file_path.name,
        "stem": _slug(file_path.stem),
        "size_bytes": size_bytes,
        "extension": suffix,
    }
    data.update(inspection)
    return data, warnings


def _detect_import(module: str) -> dict[str, Any]:
    try:
        spec = importlib.util.find_spec(module)
    except (ImportError, ValueError):
        spec = None
    origin = getattr(spec, "origin", None) if spec else None
    return {
        "ok": spec is not None,
        "module": module,
        "origin": _scrub_text(origin) if origin else None,
    }


def _detect_executable(name: str) -> dict[str, Any]:
    path = shutil.which(name)
    return {
        "ok": bool(path),
        "executable": name,
        "path": _path_display(path) if path else None,
    }


def detect_optional_tools() -> dict[str, Any]:
    imports = {module: _detect_import(module) for module in IMPORT_MODULES}
    executables = {name: _detect_executable(name) for name in EXECUTABLES}
    capabilities = {
        "pdbfixer_receptor_repair": imports["pdbfixer"]["ok"],
        "openbabel_conversion": (
            imports["openbabel"]["ok"]
            or executables["obabel"]["ok"]
            or executables["openbabel"]["ok"]
        ),
        "meeko_ligand_prep": imports["meeko"]["ok"],
        "prepare_receptor": executables["prepare_receptor"]["ok"],
        "prepare_ligand": executables["prepare_ligand"]["ok"],
        "vina_docking": executables["vina"]["ok"],
    }
    return {
        "imports": imports,
        "executables": executables,
        "capabilities": capabilities,
    }


def _tool_ok(tools: dict[str, Any], kind: str, name: str) -> bool:
    return bool(tools[kind][name]["ok"])


def _build_warnings(receptor: dict[str, Any], ligand: dict[str, Any],
                    keep_water: bool, keep_cofactors: bool) -> list[dict[str, str]]:
    warnings: list[dict[str, str]] = []

    for item in (receptor, ligand):
        role = item["role"]
        role_title = role.title()
        water_records = item.get("water_records", 0)
        if water_records:
            if keep_water:
                warnings.append(_warning(
                    "water_kept",
                    "info",
                    role,
                    f"{role_title} contains {water_records} water records and --keep-water was set.",
                    "Keep only waters with clear structural or mechanistic support.",
                ))
            else:
                warnings.append(_warning(
                    "water_removed_by_default",
                    "warning",
                    role,
                    f"{role_title} contains {water_records} water records; the default plan removes waters.",
                    "Review conserved waters before deletion and rerun with --keep-water if needed.",
                ))

    metal_atoms = receptor.get("metal_atoms", 0)
    if metal_atoms:
        warnings.append(_warning(
            "metals_kept" if keep_cofactors else "metal_review_required",
            "info" if keep_cofactors else "warning",
            "receptor",
            f"Receptor contains {metal_atoms} metal atoms or metal-like residues.",
            "Confirm charges, coordination, and whether the docking engine should retain these atoms.",
        ))

    cofactor_records = receptor.get("cofactor_records", 0)
    if cofactor_records:
        warnings.append(_warning(
            "cofactors_kept" if keep_cofactors else "cofactor_review_required",
            "info" if keep_cofactors else "warning",
            "receptor",
            f"Receptor contains {cofactor_records} non-water HETATM records that may be cofactors or bound ligands.",
            "Retain biologically required cofactors and remove unrelated crystallographic ligands before docking.",
        ))

    if receptor.get("atom_count", 0) and receptor.get("hydrogen_atoms", 0) == 0:
        warnings.append(_warning(
            "receptor_missing_hydrogens",
            "warning",
            "receptor",
            "Receptor has no detected hydrogen atoms.",
            "Add hydrogens and assign protonation states before generating receptor PDBQT.",
        ))

    if ligand.get("atom_count", 0) and ligand.get("hydrogen_atoms", 0) == 0:
        warnings.append(_warning(
            "ligand_missing_hydrogens",
            "warning",
            "ligand",
            "Ligand has no detected hydrogen atoms.",
            "Generate the intended protonation/tautomer state before ligand PDBQT preparation.",
        ))

    titratable = ligand.get("titratable_elements", [])
    if titratable:
        warnings.append(_warning(
            "ligand_protonation_ambiguous",
            "warning",
            "ligand",
            f"Ligand contains protonation-sensitive elements: {', '.join(titratable)}.",
            "Choose a pH-appropriate charge, protonation, and tautomer state before docking.",
        ))

    return warnings


def _command(command_id: str, title: str, command: str, requires: list[str],
             available: bool, notes: list[str] | None = None) -> dict[str, Any]:
    return {
        "id": command_id,
        "title": title,
        "command": command,
        "requires": requires,
        "available": available,
        "executes_by_default": False,
        "notes": notes or [],
    }


def build_commands(receptor: dict[str, Any], ligand: dict[str, Any],
                   tools: dict[str, Any], keep_water: bool,
                   keep_cofactors: bool) -> list[dict[str, Any]]:
    receptor_input = receptor["name"]
    ligand_input = ligand["name"]
    receptor_fixed = f"{receptor['stem']}_fixed.pdb"
    receptor_pdbqt = f"{receptor['stem']}_prepared.pdbqt"
    ligand_normalized = f"{ligand['stem']}_prepared.sdf"
    ligand_pdbqt = f"{ligand['stem']}_prepared.pdbqt"
    vina_config = "vina_box.txt"
    vina_out = f"{ligand['stem']}_docked.pdbqt"

    cleanup_notes = []
    cleanup_notes.append("Retain waters intentionally." if keep_water else "Remove waters unless they are known structural waters.")
    cleanup_notes.append(
        "Retain metals/cofactors intentionally." if keep_cofactors
        else "Review metals/cofactors before removing or retaining them."
    )

    commands = [
        _command(
            "receptor_repair",
            "Repair receptor model",
            (
                "python -m pdbfixer "
                f"{_quote(receptor_input)} --add-atoms=heavy --replace-nonstandard "
                f"--output={_quote(receptor_fixed)}"
            ),
            ["pdbfixer import"],
            _tool_ok(tools, "imports", "pdbfixer"),
            cleanup_notes,
        ),
        _command(
            "receptor_pdbqt",
            "Prepare receptor PDBQT",
            (
                "prepare_receptor "
                f"-r {_quote(receptor_fixed)} -o {_quote(receptor_pdbqt)} -A hydrogens"
            ),
            ["prepare_receptor executable"],
            _tool_ok(tools, "executables", "prepare_receptor"),
            ["Inspect retained hetero atoms and charges before using the PDBQT."],
        ),
        _command(
            "ligand_normalize",
            "Normalize ligand geometry and protonation",
            (
                "obabel "
                f"{_quote(ligand_input)} -O {_quote(ligand_normalized)} --gen3d -p 7.4"
            ),
            ["obabel or openbabel executable"],
            (
                _tool_ok(tools, "executables", "obabel")
                or _tool_ok(tools, "executables", "openbabel")
            ),
            ["Replace pH 7.4 with the assay-relevant pH when appropriate."],
        ),
        _command(
            "ligand_pdbqt",
            "Prepare ligand PDBQT",
            (
                "prepare_ligand "
                f"-l {_quote(ligand_normalized)} -o {_quote(ligand_pdbqt)} -A hydrogens"
            ),
            ["prepare_ligand executable or meeko import"],
            (
                _tool_ok(tools, "executables", "prepare_ligand")
                or _tool_ok(tools, "imports", "meeko")
            ),
            ["If using Meeko directly, generate an equivalent ligand PDBQT from the normalized ligand."],
        ),
        _command(
            "vina_dock",
            "Run docking after prep review",
            (
                "vina "
                f"--receptor {_quote(receptor_pdbqt)} --ligand {_quote(ligand_pdbqt)} "
                f"--config {_quote(vina_config)} --out {_quote(vina_out)} --log vina.log"
            ),
            ["vina executable", "reviewed docking box config"],
            _tool_ok(tools, "executables", "vina"),
            ["Generate or review vina_box.txt before running docking."],
        ),
    ]
    return commands


def build_prep_plan(receptor: str | Path, ligand: str | Path,
                    keep_water: bool = False, keep_cofactors: bool = False,
                    tools: dict[str, Any] | None = None) -> dict[str, Any]:
    receptor_data, receptor_warnings = validate_input(receptor, "receptor")
    ligand_data, ligand_warnings = validate_input(ligand, "ligand")
    detected_tools = tools if tools is not None else detect_optional_tools()
    warnings = []
    warnings.extend(receptor_warnings)
    warnings.extend(ligand_warnings)
    warnings.extend(_build_warnings(receptor_data, ligand_data, keep_water, keep_cofactors))

    commands = build_commands(receptor_data, ligand_data, detected_tools, keep_water, keep_cofactors)
    data = {
        "planner": "dock_prep",
        "executes_tools": False,
        "inputs": {
            "receptor": receptor_data,
            "ligand": ligand_data,
        },
        "cleanup": {
            "keep_water": keep_water,
            "keep_cofactors": keep_cofactors,
            "water_policy": "keep" if keep_water else "remove_by_default",
            "cofactor_policy": "keep" if keep_cofactors else "review_before_removal",
        },
        "tools": detected_tools,
        "warnings": warnings,
        "commands": commands,
        "outputs": {
            "written": False,
            "directory": None,
            "json": None,
            "markdown": None,
        },
    }
    return _ok_payload(_scrub_value(data))


def render_markdown(plan: dict[str, Any]) -> str:
    data = plan["data"] if "data" in plan else plan
    receptor = data["inputs"]["receptor"]
    ligand = data["inputs"]["ligand"]
    lines = [
        "# Docking Prep Plan",
        "",
        "This plan validates inputs and recommends commands only; no prep or docking tools were executed.",
        "",
        "## Inputs",
        "",
        f"- Receptor: {receptor['path']} ({receptor['format']}, {receptor['atom_count']} inspected atoms)",
        f"- Ligand: {ligand['path']} ({ligand['format']}, {ligand['atom_count']} inspected atoms)",
        "",
        "## Cleanup Policy",
        "",
        f"- Waters: {data['cleanup']['water_policy']}",
        f"- Metals/cofactors: {data['cleanup']['cofactor_policy']}",
        "",
        "## Warnings",
        "",
    ]

    if data["warnings"]:
        for warning in data["warnings"]:
            lines.append(
                f"- [{warning['severity'].upper()}] {warning['code']} ({warning['source']}): "
                f"{warning['message']} {warning['recommendation']}"
            )
    else:
        lines.append("- None detected by the stdlib planner.")

    lines.extend([
        "",
        "## Optional Tools",
        "",
        "| Tool | Kind | Status |",
        "| --- | --- | --- |",
    ])
    for name, result in data["tools"]["imports"].items():
        lines.append(f"| {name} | import | {'detected' if result['ok'] else 'not detected'} |")
    for name, result in data["tools"]["executables"].items():
        lines.append(f"| {name} | executable | {'detected' if result['ok'] else 'not detected'} |")

    lines.extend([
        "",
        "## Recommended Commands",
        "",
    ])
    for index, command in enumerate(data["commands"], start=1):
        status = "available" if command["available"] else "tool not detected"
        lines.extend([
            f"{index}. {command['title']} ({status})",
            "",
            "```bash",
            command["command"],
            "```",
            "",
        ])
        for note in command["notes"]:
            lines.append(f"   - {note}")
        if command["notes"]:
            lines.append("")

    if data["outputs"]["written"]:
        lines.extend([
            "## Written Outputs",
            "",
            f"- JSON: {data['outputs']['json']}",
            f"- Markdown: {data['outputs']['markdown']}",
            "",
        ])

    return "\n".join(lines).rstrip() + "\n"


def write_plan_outputs(plan: dict[str, Any], outdir: str | Path) -> dict[str, Any]:
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_name = "dock_prep_plan.json"
    markdown_name = "dock_prep_plan.md"

    plan["outputs"] = {
        "written": True,
        "directory": _path_display(output_dir),
        "json": json_name,
        "markdown": markdown_name,
    }
    plan["data"]["outputs"] = plan["outputs"]
    public_plan = _scrub_value(plan)
    markdown = render_markdown(public_plan)

    (output_dir / json_name).write_text(json.dumps(public_plan, indent=2) + "\n", encoding="utf-8")
    (output_dir / markdown_name).write_text(markdown, encoding="utf-8")
    return public_plan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Create a local, non-executing docking prep plan for receptor and ligand inputs.",
        epilog=(
            "Examples:\n"
            "  %(prog)s --receptor receptor.pdb --ligand ligand.sdf\n"
            "  %(prog)s --receptor receptor.pdb --ligand ligand.sdf --json\n"
            "  %(prog)s --receptor receptor.pdb --ligand ligand.sdf --outdir prep-plan --keep-cofactors"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--receptor", required=True, help="Local receptor file path")
    parser.add_argument("--ligand", required=True, help="Local ligand file path")
    parser.add_argument("--outdir", help="Write dock_prep_plan.json and dock_prep_plan.md to this directory")
    parser.add_argument("--keep-water", action="store_true", help="Plan to retain water records")
    parser.add_argument("--keep-cofactors", action="store_true", help="Plan to retain metals and cofactors")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args(argv)

    try:
        plan = build_prep_plan(
            args.receptor,
            args.ligand,
            keep_water=args.keep_water,
            keep_cofactors=args.keep_cofactors,
        )
        if args.outdir:
            plan = write_plan_outputs(plan, args.outdir)
    except (DockPrepError, OSError) as exc:
        if args.json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {_scrub_text(str(exc))}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(plan, indent=2))
    else:
        print(render_markdown(plan), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
