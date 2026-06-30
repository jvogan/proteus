#!/usr/bin/env python3
"""Summarize simple protein-ligand contacts from PDB/mmCIF files or explicit PDB IDs.

This helper is intentionally standard-library only. It provides a local
preflight contact report and detects optional PLIP/ProLIF integrations without
requiring or installing them.

Usage:
    python3 interaction_report.py complex.cif --json
    python3 interaction_report.py complex.pdb --ligand ATP --cutoff 4.5
    python3 interaction_report.py 1HSG --ligand MK1 --json
"""

import argparse
import importlib.util
import json
import math
import re
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import structure_atoms


RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}

CLASSIFICATION_NAMES = (
    "close_contact_clash",
    "polar_candidate",
    "hydrophobic_candidate",
    "metal_candidate",
    "contact",
)

POLAR_ELEMENTS = {"N", "O", "S", "P"}
HYDROPHOBIC_ELEMENTS = {"C", "S", "F", "CL", "BR", "I"}
METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS", "BE", "MG", "CA", "SR", "BA", "SC", "TI",
    "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "Y", "ZR", "MO", "RU",
    "RH", "PD", "AG", "CD", "W", "PT", "AU", "HG", "AL", "GA", "IN", "SN",
    "PB",
}
TWO_LETTER_ELEMENTS = {
    "HE", "LI", "BE", "NE", "NA", "MG", "AL", "SI", "CL", "AR", "CA", "SC",
    "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE",
    "BR", "KR", "RB", "SR", "ZR", "MO", "TC", "RU", "RH", "PD", "AG", "CD",
    "IN", "SN", "SB", "TE", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "SM",
    "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "RE",
    "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI",
}
ONE_LETTER_ELEMENTS = {"H", "B", "C", "N", "O", "F", "P", "S", "K", "V", "Y", "I", "W", "U"}
PROTEIN_ATOM_FIRST_LETTER_ELEMENTS = {"H", "C", "N", "O", "P", "S"}
OPTIONAL_TOOLS = {
    "plip": {
        "executables": ("plip", "plipcmd"),
        "module": "plip",
    },
    "prolif": {
        "executables": ("prolif",),
        "module": "prolif",
    },
}


class InteractionReportError(RuntimeError):
    pass


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _looks_like_pdb_id(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()))


def _normalize_ligand_filters(values: list[str] | None) -> list[str]:
    if not values:
        return []
    ligands = []
    seen = set()
    for value in values:
        for item in value.split(","):
            code = item.strip().upper()
            if not code:
                continue
            if not re.fullmatch(r"[A-Z0-9]{1,10}", code):
                raise ValueError(f"Invalid ligand code '{item.strip()}'.")
            if code not in seen:
                seen.add(code)
                ligands.append(code)
    return ligands


def _safe_path_label(path: Path, original: str | None = None) -> str:
    source = Path(original).name if original else path.name
    return source or path.name or "input.pdb"


def _download_pdb_id(pdb_id: str, outdir: str) -> tuple[Path, dict]:
    output_dir = Path(outdir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"{pdb_id}.pdb"
    cached = destination.exists()
    url = RCSB_DOWNLOAD.format(pdb_id=pdb_id)
    if not cached:
        request = urllib.request.Request(url, headers={"User-Agent": "proteus-skill/1.0"})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                destination.write_bytes(response.read())
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise InteractionReportError(f"PDB ID '{pdb_id}' was not found at RCSB.") from exc
            raise InteractionReportError(f"RCSB download returned HTTP {exc.code}: {exc.reason}") from exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise InteractionReportError(f"Failed to download PDB file for {pdb_id}: {exc}") from exc
    return destination, {
        "kind": "pdb_id",
        "query": pdb_id,
        "label": pdb_id,
        "downloaded_label": destination.name,
        "url": url,
        "cached": cached,
    }


def _prepare_input(query: str, outdir: str) -> tuple[Path, dict]:
    path = Path(query).expanduser()
    if path.exists():
        if not path.is_file():
            raise ValueError("Input path exists but is not a file.")
        return path, {
            "kind": "local_file",
            "label": _safe_path_label(path, query),
        }

    if not _looks_like_pdb_id(query):
        raise ValueError("Input must be a local PDB/mmCIF file or an explicit four-character PDB ID.")

    return _download_pdb_id(query.strip().upper(), outdir)


def _clean_element(value: str) -> str:
    letters = "".join(ch for ch in value.strip().upper() if ch.isalpha())
    return letters[:2]


def _infer_element(record: str, atom_name: str, explicit_element: str, resname: str) -> str:
    element = _clean_element(explicit_element)
    if element:
        return element

    letters = _clean_element(atom_name)
    if not letters:
        return ""

    resname_letters = _clean_element(resname)
    if record == "HETATM" and letters[:2] in METAL_ELEMENTS and letters[:2] == resname_letters[:2]:
        return letters[:2]

    if record == "ATOM" and letters[0] in PROTEIN_ATOM_FIRST_LETTER_ELEMENTS:
        return letters[0]

    if len(letters) >= 2 and letters[:2] in TWO_LETTER_ELEMENTS:
        return letters[:2]
    if letters[0] in ONE_LETTER_ELEMENTS:
        return letters[0]
    return letters[0]


def _residue_label(chain: str, residue_id: str, insertion_code: str) -> str:
    return f"{chain}:{residue_id}{insertion_code or ''}"


def _sort_residue_id(value: str) -> tuple[int, int | str]:
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)


def _group_sort_key(group: dict) -> tuple:
    return (
        group["ligand"],
        group["chain"],
        _sort_residue_id(group["residue_id"]),
        group["insertion_code"],
    )


def _protein_residue_key(atom: dict) -> tuple[str, str, str, str]:
    return (atom["chain"], atom["residue_id"], atom["insertion_code"], atom["resname"])


def _atom_ref(atom: dict, include_ligand: bool = False) -> dict:
    output = {
        "atom_name": atom["atom_name"],
        "element": atom["element"],
        "chain": atom["chain"],
        "residue_id": atom["residue_id"],
        "insertion_code": atom["insertion_code"],
        "residue_label": atom["residue_label"],
        "resname": atom["resname"],
    }
    if include_ligand:
        output["ligand"] = atom["resname"]
    return output


def _distance(a: dict, b: dict) -> float:
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _empty_classification_counts() -> dict:
    return {name: 0 for name in CLASSIFICATION_NAMES}


def _is_polar_pair(protein_atom: dict, ligand_atom: dict) -> bool:
    return protein_atom["element"] in POLAR_ELEMENTS and ligand_atom["element"] in POLAR_ELEMENTS


def _is_hydrophobic_pair(protein_atom: dict, ligand_atom: dict) -> bool:
    return (
        protein_atom["element"] in HYDROPHOBIC_ELEMENTS
        and ligand_atom["element"] in HYDROPHOBIC_ELEMENTS
    )


def _is_metal_pair(protein_atom: dict, ligand_atom: dict) -> bool:
    return protein_atom["element"] in METAL_ELEMENTS or ligand_atom["element"] in METAL_ELEMENTS


def classify_contact(protein_atom: dict, ligand_atom: dict, distance: float,
                     thresholds: dict) -> list[str]:
    classes = []
    if distance <= thresholds["close_contact_clash"]:
        classes.append("close_contact_clash")
    if distance <= thresholds["polar_candidate"] and _is_polar_pair(protein_atom, ligand_atom):
        classes.append("polar_candidate")
    if distance <= thresholds["hydrophobic_candidate"] and _is_hydrophobic_pair(protein_atom, ligand_atom):
        classes.append("hydrophobic_candidate")
    if distance <= thresholds["metal_candidate"] and _is_metal_pair(protein_atom, ligand_atom):
        classes.append("metal_candidate")
    if not classes:
        classes.append("contact")
    return classes


def parse_pdb_atoms(path: str | Path, ligand_filters: list[str] | None = None) -> dict:
    return structure_atoms.parse_pdb_atoms(
        path,
        ligand_filters=ligand_filters,
        water_names=WATER_NAMES,
    )


def parse_structure_atoms(path: str | Path, ligand_filters: list[str] | None = None) -> dict:
    return structure_atoms.parse_structure_atoms(
        path,
        ligand_filters=ligand_filters,
        water_names=WATER_NAMES,
    )


def _base_ligand_groups(ligand_atoms: list[dict]) -> list[dict]:
    groups: dict[tuple[str, str, str, str], dict] = {}
    for atom in ligand_atoms:
        key = (atom["resname"], atom["chain"], atom["residue_id"], atom["insertion_code"])
        group = groups.setdefault(key, {
            "ligand": atom["resname"],
            "chain": atom["chain"],
            "residue_id": atom["residue_id"],
            "insertion_code": atom["insertion_code"],
            "residue_label": atom["residue_label"],
            "atom_count": 0,
            "atom_names": [],
            "elements": [],
            "_atoms": [],
        })
        group["atom_count"] += 1
        group["atom_names"].append(atom["atom_name"])
        if atom["element"]:
            group["elements"].append(atom["element"])
        group["_atoms"].append(atom)

    ligand_groups = sorted(groups.values(), key=_group_sort_key)
    for group in ligand_groups:
        group["atom_names"] = sorted(set(group["atom_names"]))
        group["elements"] = sorted(set(group["elements"]))
    return ligand_groups


def _residue_contact(protein_atom: dict, distance: float, classes: list[str],
                     ligand_atom: dict) -> dict:
    return {
        "chain": protein_atom["chain"],
        "residue_id": protein_atom["residue_id"],
        "insertion_code": protein_atom["insertion_code"],
        "residue_label": protein_atom["residue_label"],
        "resname": protein_atom["resname"],
        "min_distance": distance,
        "classifications": set(classes),
        "ligand_atom": ligand_atom["atom_name"],
        "protein_atom": protein_atom["atom_name"],
    }


def summarize_ligand_contacts(protein_atoms: list[dict], ligand_atoms: list[dict],
                              thresholds: dict, max_contacts: int = 10) -> tuple[list[dict], dict]:
    ligand_groups = _base_ligand_groups(ligand_atoms)
    aggregate_counts = _empty_classification_counts()
    total_contacts = 0
    cutoff_sq = thresholds["contact_cutoff"] * thresholds["contact_cutoff"]

    for group in ligand_groups:
        contacts = []
        residue_contacts: dict[tuple[str, str, str, str], dict] = {}
        group_counts = _empty_classification_counts()

        for ligand_atom in group["_atoms"]:
            for protein_atom in protein_atoms:
                dx = ligand_atom["x"] - protein_atom["x"]
                dy = ligand_atom["y"] - protein_atom["y"]
                dz = ligand_atom["z"] - protein_atom["z"]
                distance_sq = dx * dx + dy * dy + dz * dz
                if distance_sq > cutoff_sq:
                    continue

                distance = math.sqrt(distance_sq)
                classes = classify_contact(protein_atom, ligand_atom, distance, thresholds)
                for name in classes:
                    group_counts[name] += 1
                    aggregate_counts[name] += 1
                total_contacts += 1

                contact = {
                    "distance": round(distance, 2),
                    "classifications": classes,
                    "ligand_atom": _atom_ref(ligand_atom, include_ligand=True),
                    "protein_atom": _atom_ref(protein_atom),
                }
                contacts.append(contact)

                residue_key = _protein_residue_key(protein_atom)
                residue_contact = residue_contacts.get(residue_key)
                if residue_contact is None or distance < residue_contact["min_distance"]:
                    residue_contacts[residue_key] = _residue_contact(
                        protein_atom, distance, classes, ligand_atom
                    )
                else:
                    residue_contact["classifications"].update(classes)

        contacts.sort(key=lambda item: (item["distance"], item["protein_atom"]["residue_label"],
                                        item["ligand_atom"]["atom_name"]))
        residues = []
        for item in residue_contacts.values():
            item["min_distance"] = round(item["min_distance"], 2)
            item["classifications"] = sorted(item["classifications"])
            residues.append(item)
        residues.sort(key=lambda item: (item["min_distance"], item["chain"],
                                        _sort_residue_id(item["residue_id"]), item["resname"]))

        group.pop("_atoms")
        group["contact_count"] = len(contacts)
        group["contact_counts"] = group_counts
        group["contacting_residue_count"] = len(residues)
        group["contacting_residues"] = residues
        group["closest_contacts"] = contacts[:max_contacts]

    summary = {
        "contact_count": total_contacts,
        "classification_counts": aggregate_counts,
    }
    return ligand_groups, summary


def detect_optional_tools() -> dict:
    tools = {}
    for name, config in OPTIONAL_TOOLS.items():
        executable_hits = []
        for candidate in config["executables"]:
            path = shutil.which(candidate)
            executable_hits.append({
                "name": candidate,
                "path_label": Path(path).name if path else None,
                "ok": bool(path),
            })

        module_name = config["module"]
        spec = importlib.util.find_spec(module_name)
        executable_ok = any(item["ok"] for item in executable_hits)
        module_ok = spec is not None
        tools[name] = {
            "available": executable_ok or module_ok,
            "executable": {
                "ok": executable_ok,
                "candidates": executable_hits,
                "path_label": next((item["path_label"] for item in executable_hits if item["path_label"]), None),
            },
            "module": {
                "ok": module_ok,
                "name": module_name,
                "origin_label": Path(spec.origin).name if spec and spec.origin else None,
            },
        }
    return tools


def _validate_thresholds(contact_cutoff: float, close_cutoff: float, polar_cutoff: float,
                         hydrophobic_cutoff: float, metal_cutoff: float) -> dict:
    thresholds = {
        "contact_cutoff": contact_cutoff,
        "close_contact_clash": close_cutoff,
        "polar_candidate": polar_cutoff,
        "hydrophobic_candidate": hydrophobic_cutoff,
        "metal_candidate": metal_cutoff,
    }
    for name, value in thresholds.items():
        if value < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be greater than or equal to 0.")
    if close_cutoff > contact_cutoff:
        raise ValueError("--close-cutoff must be less than or equal to --cutoff.")
    return thresholds


def analyze_interactions(query: str, outdir: str = ".", ligand_filters: list[str] | None = None,
                         cutoff: float = 4.0, close_cutoff: float = 2.0,
                         polar_cutoff: float = 3.5, hydrophobic_cutoff: float = 4.0,
                         metal_cutoff: float = 3.0, max_contacts: int = 10) -> dict:
    if max_contacts < 0:
        raise ValueError("--max-contacts must be greater than or equal to 0.")

    ligand_filters = ligand_filters or []
    thresholds = _validate_thresholds(
        cutoff,
        close_cutoff,
        polar_cutoff,
        hydrophobic_cutoff,
        metal_cutoff,
    )
    pdb_path, provenance = _prepare_input(query, outdir)
    parsed = parse_structure_atoms(pdb_path, ligand_filters)
    protein_atoms = parsed.pop("protein_atoms")
    ligand_atoms = parsed.pop("ligand_atoms")
    if not ligand_atoms:
        if ligand_filters:
            raise InteractionReportError(
                f"No non-water HETATM atoms matched ligand filter: {', '.join(ligand_filters)}."
            )
        raise InteractionReportError("No non-water HETATM ligand atoms were found.")

    ligand_groups, contact_summary = summarize_ligand_contacts(
        protein_atoms,
        ligand_atoms,
        thresholds,
        max_contacts=max_contacts,
    )
    notes = []
    if not protein_atoms:
        notes.append("No ATOM protein records were found; ligand contacts are empty.")
    if not contact_summary["contact_count"]:
        notes.append("No protein-ligand atom pairs were found within the contact cutoff.")

    data = {
        "input": provenance,
        "source_label": provenance["label"],
        "ligand_filter": ligand_filters,
        "parameters": {
            "thresholds_angstrom": thresholds,
            "max_contacts_per_ligand": max_contacts,
        },
        "protein_atom_count": len(protein_atoms),
        "ligand_atom_count": len(ligand_atoms),
        "ligand_group_count": len(ligand_groups),
        "ligand_groups": ligand_groups,
        "optional_tools": detect_optional_tools(),
        "notes": notes,
    }
    data.update(parsed)
    data.update(contact_summary)
    return _ok_payload(data)


def _md_escape(value: object) -> str:
    text = str(value)
    return text.replace("|", "\\|")


def _format_tool_status(tool: dict) -> str:
    executable = "yes" if tool["executable"]["ok"] else "no"
    module = "yes" if tool["module"]["ok"] else "no"
    return f"executable: {executable}, import: {module}"


def format_markdown_report(output: dict) -> str:
    if output.get("status") != "ok":
        return f"# Protein-Ligand Interaction Report\n\nERROR: {output.get('error', 'unknown error')}\n"

    data = output["data"]
    thresholds = data["parameters"]["thresholds_angstrom"]
    lines = [
        "# Protein-Ligand Interaction Report",
        "",
        f"- Source: `{_md_escape(data['source_label'])}`",
        f"- Protein atoms: {data['protein_atom_count']}",
        f"- Ligand atoms: {data['ligand_atom_count']}",
        f"- Ligand groups: {data['ligand_group_count']}",
        f"- Atom-atom contacts within {thresholds['contact_cutoff']:.2f} Angstrom: {data['contact_count']}",
    ]
    if data["ligand_filter"]:
        lines.append(f"- Ligand filter: {', '.join(data['ligand_filter'])}")
    for note in data["notes"]:
        lines.append(f"- Note: {note}")

    lines.extend([
        "",
        "## Classification Counts",
        "",
        "| Classification | Count |",
        "| --- | ---: |",
    ])
    for name in CLASSIFICATION_NAMES:
        lines.append(f"| `{name}` | {data['classification_counts'][name]} |")

    lines.extend([
        "",
        "## Optional Integrations",
        "",
        "| Tool | Detected |",
        "| --- | --- |",
    ])
    for name, tool in data["optional_tools"].items():
        lines.append(f"| {name.upper()} | {_format_tool_status(tool)} |")

    lines.extend(["", "## Ligand Groups"])
    if not data["ligand_groups"]:
        lines.extend(["", "No ligand groups were found."])
    for group in data["ligand_groups"]:
        title = f"{group['ligand']} {group['residue_label']}"
        lines.extend([
            "",
            f"### {_md_escape(title)}",
            "",
            f"- Atoms: {group['atom_count']} ({', '.join(group['atom_names']) or 'none'})",
            f"- Elements: {', '.join(group['elements']) or 'unknown'}",
            f"- Contacts: {group['contact_count']} atom pairs across {group['contacting_residue_count']} residues",
            "",
            "| Classification | Count |",
            "| --- | ---: |",
        ])
        for name in CLASSIFICATION_NAMES:
            lines.append(f"| `{name}` | {group['contact_counts'][name]} |")

        if group["closest_contacts"]:
            lines.extend([
                "",
                "| Distance | Classes | Ligand Atom | Protein Atom | Protein Residue |",
                "| ---: | --- | --- | --- | --- |",
            ])
            for contact in group["closest_contacts"]:
                ligand_atom = contact["ligand_atom"]
                protein_atom = contact["protein_atom"]
                classes = ", ".join(contact["classifications"])
                lines.append(
                    f"| {contact['distance']:.2f} | `{_md_escape(classes)}` | "
                    f"{_md_escape(ligand_atom['atom_name'])} ({_md_escape(ligand_atom['element'] or '?')}) | "
                    f"{_md_escape(protein_atom['atom_name'])} ({_md_escape(protein_atom['element'] or '?')}) | "
                    f"{_md_escape(protein_atom['resname'])} {_md_escape(protein_atom['residue_label'])} |"
                )
        else:
            lines.extend(["", "No contacts within the cutoff."])

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Summarize simple protein-ligand contacts from a PDB/mmCIF file or explicit PDB ID.",
        epilog=(
            "Examples:\n"
            "  %(prog)s complex.cif --json\n"
            "  %(prog)s complex.pdb --ligand ATP --cutoff 4.5\n"
            "  %(prog)s 1HSG --ligand MK1 --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", help="Local PDB/mmCIF path or explicit four-character PDB ID")
    parser.add_argument("--ligand", action="append",
                        help="Filter by ligand code; may be repeated or comma-separated")
    parser.add_argument("--cutoff", type=float, default=4.0,
                        help="Maximum atom-atom contact distance in Angstroms (default: 4.0)")
    parser.add_argument("--close-cutoff", type=float, default=2.0,
                        help="Distance for close contact/clash candidates (default: 2.0)")
    parser.add_argument("--polar-cutoff", type=float, default=3.5,
                        help="Distance for polar candidates between N/O/S/P atoms (default: 3.5)")
    parser.add_argument("--hydrophobic-cutoff", type=float, default=4.0,
                        help="Distance for hydrophobic candidates between C/S/halogen atoms (default: 4.0)")
    parser.add_argument("--metal-cutoff", type=float, default=3.0,
                        help="Distance for metal coordination candidates (default: 3.0)")
    parser.add_argument("--max-contacts", type=int, default=10,
                        help="Maximum closest atom contacts to list per ligand group (default: 10)")
    parser.add_argument("--outdir", default=".", help="Output directory for explicit PDB ID downloads")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    parser.add_argument("--markdown", action="store_true", help="Emit readable Markdown (default)")
    args = parser.parse_args()

    try:
        ligand_filters = _normalize_ligand_filters(args.ligand)
        output = analyze_interactions(
            args.input,
            outdir=args.outdir,
            ligand_filters=ligand_filters,
            cutoff=args.cutoff,
            close_cutoff=args.close_cutoff,
            polar_cutoff=args.polar_cutoff,
            hydrophobic_cutoff=args.hydrophobic_cutoff,
            metal_cutoff=args.metal_cutoff,
            max_contacts=args.max_contacts,
        )
    except (ValueError, OSError, InteractionReportError) as exc:
        if args.json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print(format_markdown_report(output), end="")


if __name__ == "__main__":
    main()
