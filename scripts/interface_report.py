#!/usr/bin/env python3
"""Report protein-protein interface residues from a PDB file or PDB ID.

For each pair of polymer chains, find residues whose atoms come within a cutoff
of the other chain — the chain-chain analog of pocket_report.py (which reports
ligand pockets). Zero-dependency and PDB-format only; for mmCIF or deeper
chemistry, use ChimeraX after this preflight.

The pairwise scan is O(atoms_A x atoms_B) per chain pair with a bounding-box
prefilter, which is fine for typical complexes and small assemblies.

Usage:
    python3 interface_report.py complex.pdb --json
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


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _looks_like_pdb_id(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()))


def _parse_chain_atoms(path: str) -> dict:
    """Return {chain: [atom, ...]} for polymer ATOM records."""
    chains: dict = {}
    with open(path) as handle:
        for line in handle:
            if line[:6].strip() != "ATOM" or len(line) < 54:
                continue
            chain = line[21].strip() or "?"
            try:
                atom = {
                    "atom": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "resi": line[22:26].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            except ValueError:
                continue
            chains.setdefault(chain, []).append(atom)
    return chains


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


def analyze_interfaces(query: str, cutoff: float = 5.0, chains_filter: list = None,
                       outdir: str = ".") -> dict:
    pdb_path, provenance = _prepare_input(query, outdir)
    chain_atoms = _parse_chain_atoms(pdb_path)
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
        "file": str(Path(pdb_path).resolve()),
        "cutoff_angstrom": cutoff,
        "chains": all_chains,
        "interface_count": len(interfaces),
        "interfaces": interfaces,
    }
    if len(all_chains) < 2:
        data["note"] = "Fewer than two polymer chains found; no interfaces to report."
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(
        description="Report protein-protein interface residues from a PDB file or PDB ID.")
    parser.add_argument("input", help="Local .pdb path or four-character PDB ID")
    parser.add_argument("--chains", help="Comma-separated chain pair, e.g. A,B (default: all pairs)")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Contact cutoff in Angstroms (default: 5)")
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
        output = analyze_interfaces(args.input, args.cutoff, chains_filter, args.outdir)
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


if __name__ == "__main__":
    main()
