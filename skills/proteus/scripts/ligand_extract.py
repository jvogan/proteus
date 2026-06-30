#!/usr/bin/env python3
"""Inventory non-water HETATM ligands from local PDB/mmCIF files or PDB IDs.

This is a lightweight preflight helper: it groups HETATM records by ligand
code, chain, residue sequence, and insertion code. It does not require RDKit,
Open Babel, or any chemistry toolkit.

Usage:
    python3 ligand_extract.py complex.cif --json
    python3 ligand_extract.py 1HSG --ligand MK1 --json
    python3 ligand_extract.py complex.pdb --download-ccd both --outdir ligands --json
"""

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

import fetch_pdb
import structure_atoms


ROOT = Path(__file__).resolve().parents[1]
CCD_BASE = "https://files.rcsb.org/ligands/download"
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}


class LigandExtractError(RuntimeError):
    pass


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _looks_like_pdb_id(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()))


def _display_path(path: str | Path) -> str:
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


def _prepare_input(query: str, outdir: str) -> tuple[str, dict]:
    path = Path(query)
    if path.exists():
        return str(path), {"kind": "local_file", "query": _display_path(path)}

    if not _looks_like_pdb_id(query):
        raise ValueError("Input must be a local PDB/mmCIF file or a four-character PDB ID.")

    pdb_id = query.strip().upper()
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    url, filename = fetch_pdb.build_download_url(pdb_id, "pdb", None, False)
    destination = output_dir / filename
    cached = destination.exists()
    if not cached:
        fetch_pdb.download(url, destination)
    return str(destination), {
        "kind": "pdb_id",
        "query": pdb_id,
        "downloaded": _display_path(destination),
        "cached": cached,
    }


def _residue_label(chain: str, residue_id: str, insertion_code: str) -> str:
    suffix = insertion_code if insertion_code else ""
    return f"{chain}:{residue_id}{suffix}"


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


def _ligand_summary_from_atoms(parsed_atoms: dict) -> dict:
    groups: dict[tuple[str, str, str, str], dict] = {}
    for atom in parsed_atoms["ligand_atoms"]:
        ligand = atom["resname"]
        chain = atom["chain"]
        residue_id = atom["residue_id"]
        insertion_code = atom["insertion_code"]
        atom_name = atom["atom_name"]
        element = atom["element"]
        key = (ligand, chain, residue_id, insertion_code)
        group = groups.setdefault(key, {
            "ligand": ligand,
            "chain": chain,
            "residue_id": residue_id,
            "insertion_code": insertion_code,
            "residue_label": _residue_label(chain, residue_id, insertion_code),
            "atom_count": 0,
            "atom_names": [],
            "elements": [],
        })
        group["atom_count"] += 1
        group["atom_names"].append(atom_name)
        if element:
            group["elements"].append(element)

    ligand_groups = sorted(groups.values(), key=_group_sort_key)
    for group in ligand_groups:
        group["atom_names"] = sorted(set(group["atom_names"]))
        group["elements"] = sorted(set(group["elements"]))

    components = []
    by_ligand: dict[str, list[dict]] = {}
    for group in ligand_groups:
        by_ligand.setdefault(group["ligand"], []).append(group)

    for ligand, items in sorted(by_ligand.items()):
        components.append({
            "ligand": ligand,
            "group_count": len(items),
            "atom_count": sum(item["atom_count"] for item in items),
            "chains": sorted({item["chain"] for item in items}),
            "residue_labels": [item["residue_label"] for item in sorted(items, key=_group_sort_key)],
        })

    return {
        "format": parsed_atoms["format"],
        "hetatm_records": parsed_atoms["hetatm_records"],
        "water_hetatm_records": parsed_atoms["water_hetatm_records"],
        "malformed_hetatm_records": parsed_atoms["malformed_records"],
        "ligand_group_count": len(ligand_groups),
        "ligand_component_count": len(components),
        "ligand_atom_count": sum(group["atom_count"] for group in ligand_groups),
        "ligand_groups": ligand_groups,
        "ligand_components": components,
    }


def parse_pdb_ligands(path: str, ligand_filters: list[str] | None = None) -> dict:
    parsed_atoms = structure_atoms.parse_pdb_atoms(
        path,
        ligand_filters=ligand_filters,
        water_names=WATER_NAMES,
    )
    return _ligand_summary_from_atoms(parsed_atoms)


def parse_structure_ligands(path: str, ligand_filters: list[str] | None = None) -> dict:
    parsed_atoms = structure_atoms.parse_structure_atoms(
        path,
        ligand_filters=ligand_filters,
        water_names=WATER_NAMES,
    )
    return _ligand_summary_from_atoms(parsed_atoms)


def _ccd_targets(ligand: str, mode: str) -> list[tuple[str, str]]:
    if mode == "sdf":
        return [("sdf", f"{ligand}_model.sdf")]
    if mode == "cif":
        return [("cif", f"{ligand}.cif")]
    return [("sdf", f"{ligand}_model.sdf"), ("cif", f"{ligand}.cif")]


def download_ccd(ligands: list[str], mode: str, outdir: str) -> list[dict]:
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    downloads = []
    for ligand in ligands:
        for fmt, filename in _ccd_targets(ligand, mode):
            url = f"{CCD_BASE}/{filename}"
            destination = output_dir / filename
            cached = destination.exists()
            if not cached:
                request = urllib.request.Request(url, headers={"User-Agent": "proteus-skill/1.0"})
                try:
                    with urllib.request.urlopen(request, timeout=60) as response:
                        destination.write_bytes(response.read())
                except urllib.error.HTTPError as exc:
                    if exc.code == 404:
                        raise LigandExtractError(f"CCD file not found for ligand {ligand}: {url}") from exc
                    raise LigandExtractError(f"CCD download returned HTTP {exc.code}: {exc.reason}") from exc
                except (urllib.error.URLError, TimeoutError, OSError) as exc:
                    raise LigandExtractError(f"Failed to download CCD file {url}: {exc}") from exc
            downloads.append({
                "ligand": ligand,
                "format": fmt,
                "url": url,
                "path": _display_path(destination),
                "bytes": destination.stat().st_size,
                "cached": cached,
            })
    return downloads


def analyze_ligands(query: str, outdir: str = ".", ligand_filters: list[str] | None = None,
                    download_mode: str | None = None) -> dict:
    pdb_path, provenance = _prepare_input(query, outdir)
    ligand_filters = ligand_filters or []
    parsed = parse_structure_ligands(pdb_path, ligand_filters)
    selected_ligands = [item["ligand"] for item in parsed["ligand_components"]]
    ccd_downloads = download_ccd(selected_ligands, download_mode, outdir) if download_mode else []
    data = {
        "input": provenance,
        "file": _display_path(pdb_path),
        "ligand_filter": ligand_filters,
        "ccd_downloads": ccd_downloads,
    }
    data.update(parsed)
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(
        description="Inventory non-water HETATM ligands from a local PDB/mmCIF file or PDB ID.",
        epilog=(
            "Examples:\n"
            "  %(prog)s complex.cif --json\n"
            "  %(prog)s 1HSG --ligand MK1 --json\n"
            "  %(prog)s complex.pdb --ligand ATP,FAD --download-ccd both --outdir ligands --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", help="Local PDB/mmCIF path or four-character PDB ID")
    parser.add_argument("--ligand", action="append",
                        help="Filter by ligand code; may be repeated or comma-separated")
    parser.add_argument("--download-ccd", choices=["sdf", "cif", "both"], nargs="?", const="sdf",
                        help="Download CCD reference file(s) for selected ligands (default: sdf)")
    parser.add_argument("--outdir", default=".", help="Output directory for PDB ID and CCD downloads")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        ligand_filters = _normalize_ligand_filters(args.ligand)
        output = analyze_ligands(args.input, args.outdir, ligand_filters, args.download_ccd)
    except (ValueError, OSError, LigandExtractError, fetch_pdb.PDBFetchError) as exc:
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
        for group in data["ligand_groups"]:
            print(
                f"{group['ligand']} {group['residue_label']} "
                f"atoms: {group['atom_count']}"
            )
        for download in data["ccd_downloads"]:
            print(f"CCD {download['ligand']} {download['format']}: {download['path']}")


if __name__ == "__main__":
    main()
