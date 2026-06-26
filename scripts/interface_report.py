#!/usr/bin/env python3
"""Report protein-protein interface residues from a structure file or PDB ID.

For each pair of polymer chains, find residues whose atoms come within a cutoff
of the other chain — the chain-chain analog of pocket_report.py (which reports
ligand pockets). The parser is zero-dependency and supports local PDB/mmCIF
preflights; for deeper chemistry, use ChimeraX after this triage pass.

The pairwise scan is O(atoms_A x atoms_B) per chain pair with a bounding-box
prefilter, which is fine for typical complexes and small assemblies.

Usage:
    python3 interface_report.py complex.cif --json
    python3 interface_report.py 1BRS --chains A,D --cutoff 4.5 --json
"""

import argparse
import json
import math
import re
import sys
from itertools import combinations
from pathlib import Path

import fetch_pdb
import structure_atoms


ROOT = Path(__file__).resolve().parents[1]
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


def _atom_from_shared(atom: dict) -> dict:
    return {
        "atom": atom["atom_name"],
        "resname": atom["resname"],
        "chain": atom["chain"],
        "resi": atom["residue_id"],
        "x": atom["x"],
        "y": atom["y"],
        "z": atom["z"],
    }


def _map_chain_atoms(parsed: dict) -> tuple[dict, str]:
    chains: dict = {}
    for atom in parsed["protein_atoms"]:
        mapped = _atom_from_shared(atom)
        chains.setdefault(mapped["chain"], []).append(mapped)
    return chains, parsed["format"]


def _parse_chain_atoms(path: str) -> dict:
    """Compatibility wrapper returning {chain: [atom, ...]} for PDB ATOM records."""
    chains, _fmt = _map_chain_atoms(structure_atoms.parse_pdb_atoms(path))
    return chains


def _parse_structure_chain_atoms(path: str) -> tuple[dict, str]:
    return _map_chain_atoms(structure_atoms.parse_structure_atoms(path))


def _bbox(atoms: list) -> tuple:
    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    zs = [a["z"] for a in atoms]
    return (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))


def _bbox_too_far(b1: tuple, b2: tuple, cutoff: float) -> bool:
    """True if two bounding boxes cannot contain atoms within cutoff."""
    return (b1[0] - b2[3] > cutoff or b2[0] - b1[3] > cutoff or
            b1[1] - b2[4] > cutoff or b2[1] - b1[4] > cutoff or
            b1[2] - b2[5] > cutoff or b2[2] - b1[5] > cutoff)


def _residue_key(atom: dict) -> tuple:
    return (atom["resi"], atom["resname"])


def analyze_pair(atoms_a: list, atoms_b: list, cutoff: float) -> tuple:
    """Return (res_a, res_b, pair_min) for one chain pair.

    res_a/res_b map residue keys to their minimum distance across the interface;
    pair_min maps (residue_a, residue_b) to the closest contact distance.
    """
    cutoff_sq = cutoff * cutoff
    res_a: dict = {}
    res_b: dict = {}
    pair_min: dict = {}
    for a in atoms_a:
        ax, ay, az = a["x"], a["y"], a["z"]
        ka = _residue_key(a)
        for b in atoms_b:
            dx = ax - b["x"]
            dy = ay - b["y"]
            dz = az - b["z"]
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= cutoff_sq:
                d = math.sqrt(d2)
                kb = _residue_key(b)
                if ka not in res_a or d < res_a[ka]:
                    res_a[ka] = d
                if kb not in res_b or d < res_b[kb]:
                    res_b[kb] = d
                pk = (ka, kb)
                if pk not in pair_min or d < pair_min[pk]:
                    pair_min[pk] = d
    return res_a, res_b, pair_min


def _residue_list(distances: dict, chain: str) -> list:
    return sorted(
        [{"chain": chain, "resi": key[0], "resname": key[1], "min_distance": round(dist, 2)}
         for key, dist in distances.items()],
        key=lambda r: r["min_distance"],
    )


def _prepare_input(query: str, outdir: str) -> tuple:
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


def _opposite_chains_for(chain: str, pairs: list[tuple[str, str]]) -> list[str]:
    opposites = []
    for ca, cb in pairs:
        if ca == chain and cb not in opposites:
            opposites.append(cb)
        elif cb == chain and ca not in opposites:
            opposites.append(ca)
    return opposites


def _nearest_opposite_contacts(residue_atoms: list[dict], opposite_atoms: list[dict],
                               cutoff: float) -> list[dict]:
    contacts: dict[tuple[str, str, str], dict] = {}
    for residue_atom in residue_atoms:
        for opposite_atom in opposite_atoms:
            distance = _distance(residue_atom, opposite_atom)
            key = (opposite_atom["chain"], opposite_atom["resi"], opposite_atom["resname"])
            current = contacts.get(key)
            if current is None or distance < current["distance"]:
                contacts[key] = {
                    "chain": key[0],
                    "resi": key[1],
                    "resname": key[2],
                    "distance": distance,
                    "residue_atom": residue_atom["atom"],
                    "opposite_atom": opposite_atom["atom"],
                    "within_cutoff": distance <= cutoff,
                }

    contact_list = sorted(contacts.values(), key=lambda item: item["distance"])
    for item in contact_list:
        item["min_distance"] = round(item.pop("distance"), 2)
    return contact_list


def _opposite_chain_summary(contacts: list[dict]) -> list[dict]:
    summary: dict[str, dict] = {}
    for contact in contacts:
        chain_summary = summary.setdefault(contact["chain"], {
            "chain": contact["chain"],
            "min_distance": contact["min_distance"],
            "within_cutoff": False,
            "contact_residue_count_within_cutoff": 0,
        })
        if contact["min_distance"] < chain_summary["min_distance"]:
            chain_summary["min_distance"] = contact["min_distance"]
        if contact["within_cutoff"]:
            chain_summary["within_cutoff"] = True
            chain_summary["contact_residue_count_within_cutoff"] += 1
    return sorted(summary.values(), key=lambda item: item["min_distance"])


def _build_residue_focus(chain_atoms: dict, pairs: list[tuple[str, str]], cutoff: float,
                         residue: str | None, variant_text: str | None) -> dict | None:
    variant = _parse_variant(variant_text)
    selector = _parse_residue_selector(residue, variant)
    if selector is None:
        return None

    residue_groups: dict[tuple[str, str, str], list[dict]] = {}
    for chain, atoms in chain_atoms.items():
        for atom in atoms:
            residue_groups.setdefault((chain, atom["resi"], atom["resname"]), []).append(atom)

    candidates = []
    for (chain, resi, resname), atoms in sorted(residue_groups.items()):
        if selector["chain"] is not None and chain != selector["chain"]:
            continue
        if resi != selector["resi"]:
            continue

        all_contacts = []
        for opposite_chain in _opposite_chains_for(chain, pairs):
            opposite_atoms = chain_atoms.get(opposite_chain, [])
            all_contacts.extend(_nearest_opposite_contacts(atoms, opposite_atoms, cutoff))
        all_contacts.sort(key=lambda item: item["min_distance"])
        participates = any(item["within_cutoff"] for item in all_contacts)

        candidate = {
            "chain": chain,
            "resi": resi,
            "resname": resname,
            "atom_count": len(atoms),
            "participates_in_interface": participates,
            "opposite_contact_count": len(all_contacts),
            "opposite_contacts_within_cutoff_count": sum(1 for item in all_contacts if item["within_cutoff"]),
            "nearest_opposite_chain_contact": all_contacts[0] if all_contacts else None,
            "nearest_opposite_chain_contacts": all_contacts[:10],
            "opposite_chain_summary": _opposite_chain_summary(all_contacts),
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
        "interface_cutoff_angstrom": cutoff,
        "residue_present": bool(candidates),
        "candidate_count": len(candidates),
        "participates_in_interface": any(item["participates_in_interface"] for item in candidates),
        "candidates": candidates,
    }


def analyze_interfaces(query: str, cutoff: float = 5.0, chains_filter: list = None,
                       outdir: str = ".", residue: str | None = None,
                       variant: str | None = None) -> dict:
    pdb_path, provenance = _prepare_input(query, outdir)
    chain_atoms, fmt = _parse_structure_chain_atoms(pdb_path)
    all_chains = sorted(chain_atoms)

    if chains_filter:
        pairs = [tuple(chains_filter)]
    else:
        pairs = list(combinations(all_chains, 2))

    bboxes = {c: _bbox(atoms) for c, atoms in chain_atoms.items() if atoms}
    interfaces = []
    for ca, cb in pairs:
        if ca not in chain_atoms or cb not in chain_atoms:
            continue
        if _bbox_too_far(bboxes[ca], bboxes[cb], cutoff):
            continue
        res_a, res_b, pair_min = analyze_pair(chain_atoms[ca], chain_atoms[cb], cutoff)
        if not res_a and not res_b:
            continue
        closest = sorted(
            [{"a": {"chain": ca, "resi": ka[0], "resname": ka[1]},
              "b": {"chain": cb, "resi": kb[0], "resname": kb[1]},
              "min_distance": round(dist, 2)}
             for (ka, kb), dist in pair_min.items()],
            key=lambda r: r["min_distance"],
        )[:10]
        interfaces.append({
            "chains": [ca, cb],
            "interface_residue_count": {ca: len(res_a), cb: len(res_b)},
            "contact_pair_count": len(pair_min),
            "residues": {ca: _residue_list(res_a, ca), cb: _residue_list(res_b, cb)},
            "closest_contacts": closest,
        })

    data = {
        "input": provenance,
        "file": _display_path(Path(pdb_path).resolve()),
        "format": fmt,
        "cutoff_angstrom": cutoff,
        "chains": all_chains,
        "interface_count": len(interfaces),
        "interfaces": interfaces,
    }
    residue_focus = _build_residue_focus(chain_atoms, pairs, cutoff, residue, variant)
    if residue_focus is not None:
        data["warnings"] = [NUMBERING_WARNING]
        data["residue_focus"] = residue_focus
    if len(all_chains) < 2:
        data["note"] = "Fewer than two polymer chains found; no interfaces to report."
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(
        description="Report protein-protein interface residues from a structure file or PDB ID.")
    parser.add_argument("input", help="Local PDB/mmCIF path or four-character PDB ID")
    parser.add_argument("--chains", help="Comma-separated chain pair, e.g. A,B (default: all pairs)")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Contact cutoff in Angstroms (default: 5)")
    parser.add_argument("--residue", help="Optional direct residue focus, e.g. A:10 or 10")
    parser.add_argument("--variant", help="Optional simple substitution for residue focus, e.g. R10H")
    parser.add_argument("--outdir", default=".", help="Output directory for PDB ID downloads")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    chains_filter = None
    if args.chains:
        parts = [c.strip() for c in args.chains.split(",") if c.strip()]
        if len(parts) != 2:
            message = "--chains must name exactly two chains, e.g. A,B"
            if args.json:
                print(json.dumps(_error_payload(message), indent=2))
            else:
                print(f"ERROR: {message}", file=sys.stderr)
            sys.exit(1)
        chains_filter = parts

    try:
        output = analyze_interfaces(args.input, args.cutoff, chains_filter, args.outdir, args.residue, args.variant)
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
        print(f"Chains: {', '.join(data['chains']) or 'none'}")
        if data.get("note"):
            print(data["note"])
        for iface in data["interfaces"]:
            ca, cb = iface["chains"]
            counts = iface["interface_residue_count"]
            print(f"{ca}-{cb}: {counts[ca]}/{counts[cb]} interface residues, "
                  f"{iface['contact_pair_count']} residue contacts")
        for warning in data.get("warnings", []):
            print(f"WARNING: {warning}")
        focus = data.get("residue_focus")
        if focus:
            selected = focus["selection"]
            label = f"{selected['chain']}:{selected['resi']}" if selected["chain"] else selected["resi"]
            print(f"Residue focus: {label}")
            print(f"Participates in interface: {'yes' if focus['participates_in_interface'] else 'no'}")
            for candidate in focus["candidates"]:
                nearest = candidate.get("nearest_opposite_chain_contact")
                if nearest:
                    print(
                        f"{candidate['chain']}:{candidate['resi']} {candidate['resname']} nearest opposite-chain: "
                        f"{nearest['chain']}:{nearest['resi']} {nearest['resname']} "
                        f"{nearest['min_distance']} A"
                    )
                else:
                    print(f"{candidate['chain']}:{candidate['resi']} {candidate['resname']} nearest opposite-chain: none")


if __name__ == "__main__":
    main()
