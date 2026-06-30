#!/usr/bin/env python3
"""Compute an AutoDock Vina-style docking box from ligand atoms.

The helper intentionally uses only the Python standard library. It reads local
PDB/mmCIF files by default and downloads from RCSB only when the input is a
four-character PDB ID.

Usage:
    python3 docking_box.py complex.cif --json
    python3 docking_box.py 1HSG --ligand MK1 --json
    python3 docking_box.py complex.pdb --ligand ATP,FAD --padding 6 --config-out box.txt
"""

import argparse
import importlib.util
import json
import re
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path

import structure_atoms


ROOT = Path(__file__).resolve().parents[1]
RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}
OPTIONAL_EXECUTABLES = ("vina", "obabel", "smina", "gnina")


class DockingBoxError(RuntimeError):
    pass


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _looks_like_pdb_id(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()))


def _display_path(path: str | Path | None) -> str | None:
    if path is None:
        return None
    candidate = Path(path).expanduser()
    try:
        resolved = candidate.resolve(strict=False)
    except OSError:
        resolved = candidate.absolute()
    for root in (ROOT, Path.cwd()):
        try:
            relative = resolved.relative_to(root.resolve())
            return f"./{relative}" if str(relative) != "." else "."
        except ValueError:
            continue
    return f"{resolved.name} (absolute path omitted)"


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


def _download_pdb_id(pdb_id: str, outdir: str) -> tuple[Path, dict]:
    output_dir = Path(outdir)
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
                raise DockingBoxError(f"PDB ID '{pdb_id}' was not found at RCSB.") from exc
            raise DockingBoxError(f"RCSB download returned HTTP {exc.code}: {exc.reason}") from exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise DockingBoxError(f"Failed to download PDB file for {pdb_id}: {exc}") from exc
    return destination, {
        "kind": "pdb_id",
        "query": pdb_id,
        "url": url,
        "downloaded": _display_path(destination),
        "cached": cached,
    }


def _prepare_input(query: str, outdir: str) -> tuple[Path, dict]:
    path = Path(query)
    if path.exists():
        return path, {"kind": "local_file", "query": _display_path(path)}

    if not _looks_like_pdb_id(query):
        raise ValueError("Input must be a local PDB/mmCIF file or a four-character PDB ID.")

    return _download_pdb_id(query.strip().upper(), outdir)


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


def _axis_dict(x: float, y: float, z: float) -> dict:
    return {"x": x, "y": y, "z": z}


def _bounds_for_atoms(atoms: list[dict]) -> dict:
    xs = [atom["x"] for atom in atoms]
    ys = [atom["y"] for atom in atoms]
    zs = [atom["z"] for atom in atoms]
    return {
        "min": _axis_dict(min(xs), min(ys), min(zs)),
        "max": _axis_dict(max(xs), max(ys), max(zs)),
    }


def _ligand_atom_payload(parsed_atoms: dict) -> dict:
    return {
        "format": parsed_atoms["format"],
        "hetatm_records": parsed_atoms["hetatm_records"],
        "water_hetatm_records": parsed_atoms["water_hetatm_records"],
        "malformed_hetatm_records": parsed_atoms["malformed_records"],
        "ligand_atom_count": len(parsed_atoms["ligand_atoms"]),
        "ligand_atoms": parsed_atoms["ligand_atoms"],
    }


def parse_pdb_ligand_atoms(path: str | Path, ligand_filters: list[str] | None = None) -> dict:
    """Return selected non-water HETATM ligand atoms from a PDB file."""
    parsed_atoms = structure_atoms.parse_pdb_atoms(
        path,
        ligand_filters=ligand_filters,
        water_names=WATER_NAMES,
    )
    return _ligand_atom_payload(parsed_atoms)


def parse_structure_ligand_atoms(path: str | Path, ligand_filters: list[str] | None = None) -> dict:
    parsed_atoms = structure_atoms.parse_structure_atoms(
        path,
        ligand_filters=ligand_filters,
        water_names=WATER_NAMES,
    )
    return _ligand_atom_payload(parsed_atoms)


def summarize_ligands(atoms: list[dict]) -> tuple[list[dict], list[dict]]:
    groups: dict[tuple[str, str, str, str], dict] = {}
    for atom in atoms:
        key = (atom["ligand"], atom["chain"], atom["residue_id"], atom["insertion_code"])
        group = groups.setdefault(key, {
            "ligand": atom["ligand"],
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
        bounds = _bounds_for_atoms(group.pop("_atoms"))
        group["atom_names"] = sorted(set(group["atom_names"]))
        group["elements"] = sorted(set(group["elements"]))
        group["bounds"] = bounds

    by_ligand: dict[str, list[dict]] = {}
    for group in ligand_groups:
        by_ligand.setdefault(group["ligand"], []).append(group)

    components = []
    for ligand, items in sorted(by_ligand.items()):
        components.append({
            "ligand": ligand,
            "group_count": len(items),
            "atom_count": sum(item["atom_count"] for item in items),
            "chains": sorted({item["chain"] for item in items}),
            "residue_labels": [item["residue_label"] for item in sorted(items, key=_group_sort_key)],
        })

    return ligand_groups, components


def compute_docking_box(atoms: list[dict], padding: float) -> dict:
    if not atoms:
        raise DockingBoxError("No ligand atoms were selected.")
    if padding < 0:
        raise ValueError("--padding must be greater than or equal to 0.")

    bounds = _bounds_for_atoms(atoms)
    minimum = bounds["min"]
    maximum = bounds["max"]
    center = {
        axis: (minimum[axis] + maximum[axis]) / 2.0
        for axis in ("x", "y", "z")
    }
    extent = {
        axis: maximum[axis] - minimum[axis]
        for axis in ("x", "y", "z")
    }
    size = {
        axis: extent[axis] + (2.0 * padding)
        for axis in ("x", "y", "z")
    }
    return {
        "padding": padding,
        "min": minimum,
        "max": maximum,
        "extent": extent,
        "center": center,
        "size": size,
    }


def format_vina_config(box: dict, precision: int = 3) -> str:
    lines = []
    for axis in ("x", "y", "z"):
        lines.append(f"center_{axis} = {box['center'][axis]:.{precision}f}")
    for axis in ("x", "y", "z"):
        lines.append(f"size_{axis} = {box['size'][axis]:.{precision}f}")
    return "\n".join(lines) + "\n"


def detect_optional_tools() -> dict:
    tools = {}
    for name in OPTIONAL_EXECUTABLES:
        path = shutil.which(name)
        tools[name] = {"ok": bool(path), "path": _display_path(path) if path else None}

    spec = importlib.util.find_spec("rdkit")
    tools["rdkit"] = {
        "ok": spec is not None,
        "module": "rdkit",
        "origin": _display_path(getattr(spec, "origin", None)) if spec else None,
    }
    return tools


def _write_config(path: str, text: str) -> str:
    destination = Path(path)
    if destination.parent != Path("."):
        destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    return _display_path(destination)


def analyze_docking_box(query: str, outdir: str = ".", ligand_filters: list[str] | None = None,
                        padding: float = 4.0, config_out: str | None = None) -> dict:
    if padding < 0:
        raise ValueError("--padding must be greater than or equal to 0.")

    ligand_filters = ligand_filters or []
    pdb_path, provenance = _prepare_input(query, outdir)
    parsed = parse_structure_ligand_atoms(pdb_path, ligand_filters)
    atoms = parsed["ligand_atoms"]
    if not atoms:
        if ligand_filters:
            raise DockingBoxError(f"No non-water HETATM atoms matched ligand filter: {', '.join(ligand_filters)}.")
        raise DockingBoxError("No non-water HETATM ligand atoms were found.")

    ligand_groups, ligand_components = summarize_ligands(atoms)
    box = compute_docking_box(atoms, padding)
    config_text = format_vina_config(box)
    config_path = _write_config(config_out, config_text) if config_out else None

    data = {
        "input": provenance,
        "file": _display_path(pdb_path),
        "ligand_filter": ligand_filters,
        "ligand_group_count": len(ligand_groups),
        "ligand_component_count": len(ligand_components),
        "ligand_groups": ligand_groups,
        "ligand_components": ligand_components,
        "box": box,
        "vina_config": {
            "path": config_path,
            "text": config_text,
        },
        "tools": detect_optional_tools(),
    }
    data.update(parsed)
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(
        description="Compute a Vina-style docking box around ligand atoms in a PDB/mmCIF file.",
        epilog=(
            "Examples:\n"
            "  %(prog)s complex.cif --json\n"
            "  %(prog)s 1HSG --ligand MK1 --json\n"
            "  %(prog)s complex.pdb --ligand ATP,FAD --padding 6 --config-out box.txt"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", help="Local PDB/mmCIF path or four-character PDB ID")
    parser.add_argument("--ligand", action="append",
                        help="Filter by ligand code; may be repeated or comma-separated")
    parser.add_argument("--padding", type=float, default=4.0,
                        help="Angstrom padding added to each side of the ligand bounds (default: 4.0)")
    parser.add_argument("--outdir", default=".", help="Output directory for PDB ID downloads")
    parser.add_argument("--config-out", help="Write an AutoDock Vina text config to this path")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        ligand_filters = _normalize_ligand_filters(args.ligand)
        output = analyze_docking_box(
            args.input,
            outdir=args.outdir,
            ligand_filters=ligand_filters,
            padding=args.padding,
            config_out=args.config_out,
        )
    except (ValueError, OSError, DockingBoxError) as exc:
        if args.json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        data = output["data"]
        print(f"File: {data['file']}")
        print(f"Ligand groups: {data['ligand_group_count']}")
        print(f"Ligand atoms: {data['ligand_atom_count']}")
        if data["ligand_filter"]:
            print(f"Ligand filter: {', '.join(data['ligand_filter'])}")
        for group in data["ligand_groups"]:
            print(
                f"{group['ligand']} {group['residue_label']} "
                f"atoms: {group['atom_count']}"
            )
        print()
        print(data["vina_config"]["text"], end="")
        if data["vina_config"]["path"]:
            print(f"Config written: {data['vina_config']['path']}")
        found_tools = [name for name, result in data["tools"].items() if result["ok"]]
        print(f"Optional tools detected: {', '.join(found_tools) if found_tools else 'none'}")


if __name__ == "__main__":
    main()
