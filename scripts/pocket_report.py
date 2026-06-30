#!/usr/bin/env python3
"""Report ligand binding-pocket residues from a structure file or PDB ID.

The parser is intentionally zero-dependency and supports local PDB/mmCIF
preflights. For deep chemistry, use ChimeraX after this triage pass.

Usage:
    python3 pocket_report.py 1HSG --json
    python3 pocket_report.py complex.cif --radius 4.5 --json
"""

import argparse
import json
import math
import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import fetch_pdb
import structure_atoms


ROOT = Path(__file__).resolve().parents[1]
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}
AA_ONE_LETTER = set("ACDEFGHIKLMNPQRSTVWY")
AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
NUMBERING_WARNING = (
    "SIFTS/PDBe residue mapping was not applied; residue and variant lookups "
    "use direct structure residue numbering."
)


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
    if not candidate.is_absolute():
        return str(path)
    resolved = candidate.resolve()
    for root, prefix in ((ROOT, "."), (Path.cwd(), ".")):
        try:
            relative = resolved.relative_to(root.resolve())
        except ValueError:
            continue
        return f"{prefix}/{relative.as_posix()}"
    return f"{resolved.name} (absolute path omitted)"


def _distance(a: dict, b: dict) -> float:
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _parse_variant(value: str | None) -> dict | None:
    if not value:
        return None
    raw = value.strip()
    match = re.fullmatch(r"([A-Za-z])([0-9]+)([A-Za-z])", raw)
    if not match:
        raise ValueError("--variant must be a simple substitution like R10H.")

    from_one, index, to_one = match.groups()
    from_one = from_one.upper()
    to_one = to_one.upper()
    if from_one not in AA_ONE_LETTER or to_one not in AA_ONE_LETTER:
        raise ValueError(f"Unknown amino-acid code in variant: {value}")
    residue_index = int(index)
    return {
        "input": value,
        "short": f"{from_one}{residue_index}{to_one}",
        "from": {"one_letter": from_one},
        "to": {"one_letter": to_one},
        "residue_index": residue_index,
    }


def _parse_residue_selector(residue: str | None, variant: dict | None) -> dict | None:
    if residue:
        raw = residue.strip()
        if ":" in raw:
            chain, resi = raw.split(":", 1)
            chain = chain.strip()
            resi = resi.strip()
            if not chain or not resi:
                raise ValueError("--residue must look like A:10 or 10.")
        else:
            chain = None
            resi = raw
        if not resi:
            raise ValueError("--residue must look like A:10 or 10.")
        source = "residue"
    elif variant:
        chain = None
        resi = str(variant["residue_index"])
        raw = resi
        source = "variant"
    else:
        return None

    if variant:
        try:
            direct_index = int(resi)
        except ValueError:
            direct_index = None
        if direct_index is not None and direct_index != variant["residue_index"]:
            raise ValueError("--residue numbering does not match --variant residue index.")

    return {"input": raw, "chain": chain, "resi": resi, "source": source}


def _min_atom_distance(residue_atoms: list[dict], ligand_atoms: list[dict]) -> dict | None:
    best = None
    for protein_atom in residue_atoms:
        for ligand_atom in ligand_atoms:
            distance = _distance(protein_atom, ligand_atom)
            if best is None or distance < best["distance"]:
                best = {
                    "distance": distance,
                    "protein_atom": protein_atom["atom"],
                    "ligand_atom": ligand_atom["atom"],
                }
    return best


def _atom_from_shared(atom: dict) -> dict:
    return {
        "record": atom["record"],
        "atom": atom["atom_name"],
        "resname": atom["resname"],
        "chain": atom["chain"],
        "resi": atom["residue_id"],
        "x": atom["x"],
        "y": atom["y"],
        "z": atom["z"],
    }


def _map_parsed_atoms(parsed: dict) -> tuple[list[dict], list[dict], str]:
    protein_atoms = [_atom_from_shared(atom) for atom in parsed["protein_atoms"]]
    ligand_atoms = [_atom_from_shared(atom) for atom in parsed["ligand_atoms"]]
    return protein_atoms, ligand_atoms, parsed["format"]


def _parse_pdb_atoms(path: str) -> tuple[list[dict], list[dict]]:
    """Compatibility wrapper for callers that expect PDB-only atom lists."""
    protein_atoms, ligand_atoms, _fmt = _map_parsed_atoms(
        structure_atoms.parse_pdb_atoms(path, water_names=WATER_NAMES)
    )
    return protein_atoms, ligand_atoms


def _parse_structure_atoms(path: str) -> tuple[list[dict], list[dict], str]:
    return _map_parsed_atoms(structure_atoms.parse_structure_atoms(path, water_names=WATER_NAMES))


def _ligand_key(atom: dict) -> tuple[str, str, str]:
    return atom["chain"], atom["resi"], atom["resname"]


def _residue_key(atom: dict) -> tuple[str, str, str]:
    return atom["chain"], atom["resi"], atom["resname"]


def _prepare_input(query: str, outdir: str) -> tuple[str, dict]:
    path = Path(query)
    if path.exists():
        return str(path), {"kind": "local_file", "query": _display_path(path.resolve())}
    if not _looks_like_pdb_id(query):
        raise ValueError("Input must be a local PDB/mmCIF file or a four-character PDB ID.")

    pdb_id = query.upper()
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    url, filename = fetch_pdb.build_download_url(pdb_id, "pdb", None, False)
    destination = output_dir / filename
    if not destination.exists():
        fetch_pdb.download(url, destination)
    return str(destination), {"kind": "pdb_id", "query": pdb_id, "downloaded": _display_path(destination.resolve())}


def _build_residue_focus(protein_atoms: list[dict], ligand_groups: dict, radius: float,
                         residue: str | None, variant_text: str | None) -> dict | None:
    variant = _parse_variant(variant_text)
    selector = _parse_residue_selector(residue, variant)
    if selector is None:
        return None

    residue_groups: dict[tuple[str, str, str], list[dict]] = {}
    for atom in protein_atoms:
        residue_groups.setdefault(_residue_key(atom), []).append(atom)

    candidates = []
    for (chain, resi, resname), atoms in sorted(residue_groups.items()):
        if selector["chain"] is not None and chain != selector["chain"]:
            continue
        if resi != selector["resi"]:
            continue

        ligand_contacts = []
        for (lig_chain, lig_resi, lig_resname), ligand_atoms in sorted(ligand_groups.items()):
            closest = _min_atom_distance(atoms, ligand_atoms)
            if closest is None:
                continue
            within_cutoff = closest["distance"] <= radius
            ligand_contacts.append({
                "ligand": {
                    "chain": lig_chain,
                    "resi": lig_resi,
                    "resname": lig_resname,
                    "atom_count": len(ligand_atoms),
                },
                "min_distance": round(closest["distance"], 2),
                "ligand_atom": closest["ligand_atom"],
                "protein_atom": closest["protein_atom"],
                "within_cutoff": within_cutoff,
            })

        ligand_contacts.sort(key=lambda item: item["min_distance"])
        candidate = {
            "chain": chain,
            "resi": resi,
            "resname": resname,
            "atom_count": len(atoms),
            "within_pocket_cutoff": any(item["within_cutoff"] for item in ligand_contacts),
            "nearest_ligand_contact": ligand_contacts[0] if ligand_contacts else None,
            "nearest_ligand_contacts": ligand_contacts,
        }
        if variant:
            residue_one = AA3_TO_1.get(resname.upper())
            candidate["reference_matches_variant"] = (
                residue_one == variant["from"]["one_letter"] if residue_one else None
            )
        candidates.append(candidate)

    return {
        "query": {"residue": residue, "variant": variant_text},
        "selection": selector,
        "variant": variant,
        "pocket_cutoff_angstrom": radius,
        "residue_present": bool(candidates),
        "candidate_count": len(candidates),
        "within_pocket_cutoff": any(item["within_pocket_cutoff"] for item in candidates),
        "candidates": candidates,
    }


def analyze_pocket(query: str, radius: float = 5.0, outdir: str = ".",
                   residue: str | None = None, variant: str | None = None) -> dict:
    pdb_path, provenance = _prepare_input(query, outdir)
    protein_atoms, ligand_atoms, fmt = _parse_structure_atoms(pdb_path)
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
        "file": _display_path(Path(pdb_path).resolve()),
        "format": fmt,
        "radius_angstrom": radius,
        "protein_atom_count": len(protein_atoms),
        "ligand_atom_count": len(ligand_atoms),
        "ligand_count": len(ligands),
        "ligands": ligands,
    }
    residue_focus = _build_residue_focus(protein_atoms, ligand_groups, radius, residue, variant)
    if residue_focus is not None:
        data["warnings"] = [NUMBERING_WARNING]
        data["residue_focus"] = residue_focus
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(description="Report ligand binding-pocket residues from a structure file or PDB ID.")
    parser.add_argument("input", help="Local PDB/mmCIF path or four-character PDB ID")
    parser.add_argument("--radius", type=float, default=5.0, help="Pocket radius in Angstroms (default: 5)")
    parser.add_argument("--residue", help="Optional direct residue focus, e.g. A:10 or 10")
    parser.add_argument("--variant", help="Optional simple substitution for residue focus, e.g. R10H")
    parser.add_argument("--outdir", default=".", help="Output directory for PDB ID downloads")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        output = analyze_pocket(args.input, args.radius, args.outdir, args.residue, args.variant)
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
        for warning in data.get("warnings", []):
            print(f"WARNING: {warning}")
        focus = data.get("residue_focus")
        if focus:
            selected = focus["selection"]
            label = f"{selected['chain']}:{selected['resi']}" if selected["chain"] else selected["resi"]
            print(f"Residue focus: {label}")
            print(f"Within pocket cutoff: {'yes' if focus['within_pocket_cutoff'] else 'no'}")
            for candidate in focus["candidates"]:
                nearest = candidate.get("nearest_ligand_contact")
                if nearest:
                    ligand = nearest["ligand"]
                    print(
                        f"{candidate['chain']}:{candidate['resi']} {candidate['resname']} nearest ligand: "
                        f"{ligand['resname']} {ligand['chain']}:{ligand['resi']} "
                        f"{nearest['min_distance']} A"
                    )
                else:
                    print(f"{candidate['chain']}:{candidate['resi']} {candidate['resname']} nearest ligand: none")


if __name__ == "__main__":
    main()
