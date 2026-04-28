#!/usr/bin/env python3
"""Report ligand binding-pocket residues from a PDB file or PDB ID.

The parser is intentionally PDB-format only and zero-dependency. For mmCIF or
deep chemistry, use ChimeraX after this preflight.

Usage:
    python3 pocket_report.py 1HSG --json
    python3 pocket_report.py complex.pdb --radius 4.5 --json
"""

import argparse
import json
import math
import re
import sys
from pathlib import Path

import fetch_pdb


WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _looks_like_pdb_id(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()))


def _distance(a: dict, b: dict) -> float:
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _parse_pdb_atoms(path: str) -> tuple[list[dict], list[dict]]:
    protein_atoms = []
    ligand_atoms = []
    with open(path) as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            resname = line[17:20].strip()
            chain = line[21].strip() or "?"
            resi = line[22:26].strip()
            atom = line[12:16].strip()
            try:
                parsed = {
                    "record": record,
                    "atom": atom,
                    "resname": resname,
                    "chain": chain,
                    "resi": resi,
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            except ValueError:
                continue
            if record == "ATOM":
                protein_atoms.append(parsed)
            elif resname not in WATER_NAMES:
                ligand_atoms.append(parsed)
    return protein_atoms, ligand_atoms


def _ligand_key(atom: dict) -> tuple[str, str, str]:
    return atom["chain"], atom["resi"], atom["resname"]


def _residue_key(atom: dict) -> tuple[str, str, str]:
    return atom["chain"], atom["resi"], atom["resname"]


def _prepare_input(query: str, outdir: str) -> tuple[str, dict]:
    path = Path(query)
    if path.exists():
        return str(path), {"kind": "local_file", "query": query}
    if not _looks_like_pdb_id(query):
        raise ValueError("Input must be a local .pdb file or a four-character PDB ID.")

    pdb_id = query.upper()
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    url, filename = fetch_pdb.build_download_url(pdb_id, "pdb", None, False)
    destination = output_dir / filename
    if not destination.exists():
        fetch_pdb.download(url, destination)
    return str(destination), {"kind": "pdb_id", "query": pdb_id, "downloaded": str(destination.resolve())}


def analyze_pocket(query: str, radius: float = 5.0, outdir: str = ".") -> dict:
    pdb_path, provenance = _prepare_input(query, outdir)
    protein_atoms, ligand_atoms = _parse_pdb_atoms(pdb_path)
    ligand_groups: dict[tuple[str, str, str], list[dict]] = {}
    for atom in ligand_atoms:
        ligand_groups.setdefault(_ligand_key(atom), []).append(atom)

    ligands = []
    for (chain, resi, resname), atoms in sorted(ligand_groups.items()):
        contacts = {}
        for ligand_atom in atoms:
            for protein_atom in protein_atoms:
                d = _distance(ligand_atom, protein_atom)
                if d <= radius:
                    key = _residue_key(protein_atom)
                    current = contacts.get(key)
                    if current is None or d < current["min_distance"]:
                        contacts[key] = {
                            "chain": key[0],
                            "resi": key[1],
                            "resname": key[2],
                            "min_distance": d,
                            "ligand_atom": ligand_atom["atom"],
                            "protein_atom": protein_atom["atom"],
                        }
        contact_list = sorted(contacts.values(), key=lambda item: item["min_distance"])
        for item in contact_list:
            item["min_distance"] = round(item["min_distance"], 2)
        ligands.append({
            "ligand": {"chain": chain, "resi": resi, "resname": resname, "atom_count": len(atoms)},
            "contact_residue_count": len(contact_list),
            "contact_residues": contact_list,
        })

    data = {
        "input": provenance,
        "file": str(Path(pdb_path).resolve()),
        "radius_angstrom": radius,
        "protein_atom_count": len(protein_atoms),
        "ligand_atom_count": len(ligand_atoms),
        "ligand_count": len(ligands),
        "ligands": ligands,
    }
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(description="Report ligand binding-pocket residues from a PDB file or PDB ID.")
    parser.add_argument("input", help="Local .pdb path or four-character PDB ID")
    parser.add_argument("--radius", type=float, default=5.0, help="Pocket radius in Angstroms (default: 5)")
    parser.add_argument("--outdir", default=".", help="Output directory for PDB ID downloads")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        output = analyze_pocket(args.input, args.radius, args.outdir)
    except (ValueError, OSError, fetch_pdb.PDBFetchError) as exc:
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
        print(f"Ligands: {data['ligand_count']}")
        for ligand in data["ligands"]:
            ident = ligand["ligand"]
            print(f"{ident['resname']} {ident['chain']}:{ident['resi']} contacts: {ligand['contact_residue_count']}")


if __name__ == "__main__":
    main()
