#!/usr/bin/env python3
"""Quick PDB file inspector — Proteus skill.

Zero external dependencies (stdlib only). Useful for pre-flight checks
before loading structures into PyMOL or ChimeraX.

Usage:
    python pdb_info.py structure.pdb
    python pdb_info.py --json structure.pdb
    python pdb_info.py --help
"""

import argparse
import json
import os
import sys
from collections import defaultdict


def parse_pdb(path: str, as_json: bool = False, force_alphafold: bool = False):
    """Parse a PDB file and report structure summary.

    Args:
        force_alphafold: Treat B-factors as pLDDT confidence scores.
            Without this flag, AlphaFold detection uses heuristics
            (filename contains 'AF-' AND B-factors in 0-100 range).
    """
    chains = set()
    residues = defaultdict(set)  # chain -> set of (resn, resi)
    atom_count = 0
    hetatm_count = 0
    b_factors = []
    title_lines = []

    try:
        fh = open(path)
    except FileNotFoundError:
        msg = f"File not found: {path}"
        if as_json:
            print(json.dumps({"error": msg}))
        else:
            print(f"ERROR: {msg}", file=sys.stderr)
        sys.exit(1)
    except PermissionError:
        msg = f"Permission denied: {path}"
        if as_json:
            print(json.dumps({"error": msg}))
        else:
            print(f"ERROR: {msg}", file=sys.stderr)
        sys.exit(1)

    with fh as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "TITLE":
                title_lines.append(line[10:].strip())
            elif rec in ("ATOM", "HETATM"):
                if len(line) < 66:
                    continue  # Malformed line
                chain = line[21]
                resn = line[17:20].strip()
                resi = line[22:26].strip()
                bfac = float(line[60:66].strip())
                chains.add(chain)
                residues[chain].add((resn, resi))
                b_factors.append(bfac)
                if rec == "ATOM":
                    atom_count += 1
                else:
                    hetatm_count += 1

    title = " ".join(title_lines) if title_lines else "(no title)"
    sorted_chains = sorted(chains)

    # B-factor / pLDDT stats
    b_stats = {}
    is_alphafold = False
    if b_factors:
        b_min = min(b_factors)
        b_max = max(b_factors)
        b_mean = sum(b_factors) / len(b_factors)
        b_stats = {"min": round(b_min, 2), "max": round(b_max, 2), "mean": round(b_mean, 2)}

        # AlphaFold detection: filename heuristic + B-factor range, or explicit flag
        filename_hint = "AF-" in os.path.basename(path).upper() or "alphafold" in path.lower()
        if force_alphafold or (filename_hint and 0 <= b_min and b_max <= 100):
            is_alphafold = True
            plddt_bins = {
                "very_high_gt90": sum(1 for b in b_factors if b > 90) / len(b_factors),
                "confident_70_90": sum(1 for b in b_factors if 70 < b <= 90) / len(b_factors),
                "low_50_70": sum(1 for b in b_factors if 50 < b <= 70) / len(b_factors),
                "very_low_lt50": sum(1 for b in b_factors if b <= 50) / len(b_factors),
            }

    # Chain details
    chain_info = {}
    for ch in sorted_chains:
        res = residues[ch]
        # Standard amino acids have 3-letter codes
        aa_count = len([r for r in res if len(r[0]) == 3])
        chain_info[ch] = {"residues": aa_count}

    if as_json:
        output = {
            "file": path,
            "title": title,
            "chains": sorted_chains,
            "chain_details": chain_info,
            "atom_records": atom_count,
            "hetatm_records": hetatm_count,
            "bfactor": b_stats,
            "likely_alphafold": is_alphafold,
        }
        if is_alphafold:
            output["plddt_distribution"] = plddt_bins
        print(json.dumps(output, indent=2))
    else:
        print(f"File: {path}")
        print(f"Title: {title}")
        print(f"Chains: {sorted_chains}")
        print(f"ATOM records: {atom_count}")
        print(f"HETATM records: {hetatm_count}")
        for ch in sorted_chains:
            print(f"  Chain {ch}: {chain_info[ch]['residues']} residues")
        if b_stats:
            print(f"B-factor range: {b_stats['min']} - {b_stats['max']} (mean {b_stats['mean']})")
        if is_alphafold:
            print(f"  Likely AlphaFold pLDDT:")
            print(f"    >90 (very high): {plddt_bins['very_high_gt90']:.1%}")
            print(f"    70-90 (confident): {plddt_bins['confident_70_90']:.1%}")
            print(f"    50-70 (low): {plddt_bins['low_50_70']:.1%}")
            print(f"    <50 (very low): {plddt_bins['very_low_lt50']:.1%}")


def main():
    parser = argparse.ArgumentParser(
        description="Quick PDB file inspector. Zero dependencies.",
        epilog="Examples:\n"
               "  %(prog)s structure.pdb\n"
               "  %(prog)s --json AF-P04637-F1-model_v6.pdb",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("pdb", help="Path to PDB file")
    parser.add_argument("--json", action="store_true", dest="as_json", help="Output as JSON")
    parser.add_argument("--alphafold", action="store_true",
                        help="Treat B-factors as pLDDT confidence (skip auto-detection)")
    args = parser.parse_args()

    parse_pdb(args.pdb, args.as_json, args.alphafold)


if __name__ == "__main__":
    main()
