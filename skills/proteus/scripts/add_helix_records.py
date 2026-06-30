#!/usr/bin/env python3
"""Add HELIX records to a CA-only backbone so viewers draw helical cartoons.

CA-only models (RFdiffusion, Genie, and other backbone generators) render as
flat "spaghetti" in PyMOL/ChimeraX: the cartoon engine infers helices from
HELIX/SHEET records that these files don't contain. This detects helices from
backbone geometry alone — the CA(i)-CA(i+3) distance is ~5.0-5.4 A inside an
alpha-helix (vs ~10 A for a beta-strand) — and prepends standard PDB HELIX
records so the cartoon renders correctly. Pairs with PyMOL `set cartoon_trace_atoms, 1`.

Usage:
    python add_helix_records.py model.pdb                  # writes model_with_ss.pdb
    python add_helix_records.py model.pdb -o out.pdb --json
    python add_helix_records.py --help
"""

import argparse
import json
import math
import os
import sys


def read_ca_coords(path):
    """Return [(resnum, chain, x, y, z), ...] for CA atoms in a PDB file."""
    coords = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    resnum = int(line[22:26])
                except ValueError:
                    continue
                coords.append((resnum, line[21], x, y, z))
    return coords


def _dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def detect_helix_residues(coords, lo=4.7, hi=5.8):
    """Flag residues whose CA(i)-CA(i+3) distance matches one alpha-helix turn."""
    is_helix = [False] * len(coords)
    for i in range(len(coords) - 3):
        if lo <= _dist(coords[i][2:], coords[i + 3][2:]) <= hi:
            for k in range(i, i + 4):
                is_helix[k] = True
    return is_helix


def helix_segments(coords, is_helix, min_len=6):
    """Coalesce flagged residues into (chain, start, end) runs of >= min_len."""
    segments = []
    i = 0
    while i < len(coords):
        if is_helix[i]:
            j = i
            while j < len(coords) and is_helix[j]:
                j += 1
            if j - i >= min_len:
                segments.append((coords[i][1], coords[i][0], coords[j - 1][0]))
            i = j
        else:
            i += 1
    return segments


def _helix_record(n, chain, start, end):
    """Format a standard PDB HELIX record."""
    length = end - start + 1
    return (f"HELIX  {n:3d} {n:3d} ALA {chain} {start:4d}  "
            f"ALA {chain} {end:4d}  1{'':33}{length:5d}")


def add_helix_records(input_path, output_path, min_len=6):
    coords = read_ca_coords(input_path)
    if not coords:
        return {"status": "error", "error": f"No CA atoms found in {input_path}"}
    is_helix = detect_helix_residues(coords)
    segments = helix_segments(coords, is_helix, min_len)
    helix_lines = [_helix_record(n, c, s, e)
                   for n, (c, s, e) in enumerate(segments, 1)]
    with open(input_path) as fh:
        body = fh.read().splitlines()
    with open(output_path, "w") as fh:
        fh.write("\n".join(helix_lines + body) + "\n")
    helical = sum(is_helix)
    return {"status": "ok", "data": {
        "input": os.path.abspath(input_path),
        "output": os.path.abspath(output_path),
        "residues": len(coords),
        "helical_residues": helical,
        "helical_fraction": round(helical / len(coords), 3),
        "helices": [{"chain": c, "start": s, "end": e, "length": e - s + 1}
                    for c, s, e in segments],
    }}


def main():
    parser = argparse.ArgumentParser(
        description="Add HELIX records to a CA-only PDB from backbone geometry.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  %(prog)s model.pdb\n"
               "  %(prog)s model.pdb -o model_ss.pdb --json",
    )
    parser.add_argument("pdb", help="CA-only PDB file")
    parser.add_argument("-o", "--output", help="Output path (default: <input>_with_ss.pdb)")
    parser.add_argument("--min-len", type=int, default=6,
                        help="Minimum helix length in residues (default: 6)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    if not os.path.isfile(args.pdb):
        result = {"status": "error", "error": f"File not found: {args.pdb}"}
    else:
        output = args.output or (os.path.splitext(args.pdb)[0] + "_with_ss.pdb")
        result = add_helix_records(args.pdb, output, args.min_len)

    if args.json:
        print(json.dumps(result, indent=2))
    elif result["status"] == "ok":
        d = result["data"]
        print(f"{d['input']}: {d['residues']} residues, {d['helical_residues']} helical "
              f"({d['helical_fraction'] * 100:.1f}%), {len(d['helices'])} helices -> {d['output']}")
        for h in d["helices"]:
            print(f"  HELIX chain {h['chain']} {h['start']}-{h['end']} (len {h['length']})")
    else:
        print(result["error"], file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
