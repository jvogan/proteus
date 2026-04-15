#!/usr/bin/env python3
"""Inspect PDB or mmCIF coordinate files with no external dependencies.

This supersedes pdb_info.py for new workflows because it handles both legacy
PDB and modern mmCIF files. The parser is intentionally lightweight; for deep
crystallographic work, prefer Gemmi, ChimeraX, or PyMOL.

Usage:
    python structure_info.py structure.pdb
    python structure_info.py structure.cif --json
    python structure_info.py AF-P04637-F1.cif --alphafold
"""

import argparse
import json
import os
import shlex
import sys
from collections import defaultdict
from pathlib import Path


def _bfactor_stats(values: list[float]) -> dict:
    if not values:
        return {}
    return {
        "min": round(min(values), 2),
        "max": round(max(values), 2),
        "mean": round(sum(values) / len(values), 2),
    }


def _plddt_distribution(values: list[float]) -> dict:
    return {
        "very_high_gt90": sum(1 for value in values if value > 90) / len(values),
        "confident_70_90": sum(1 for value in values if 70 < value <= 90) / len(values),
        "low_50_70": sum(1 for value in values if 50 < value <= 70) / len(values),
        "very_low_le50": sum(1 for value in values if value <= 50) / len(values),
    }


def _is_likely_alphafold(path: str, b_factors: list[float], force: bool) -> bool:
    if force:
        return True
    if not b_factors:
        return False
    filename_hint = "AF-" in os.path.basename(path).upper() or "alphafold" in path.lower()
    return filename_hint and min(b_factors) >= 0 and max(b_factors) <= 100


def parse_pdb(path: str, force_alphafold: bool = False) -> dict:
    chains = set()
    residues = defaultdict(set)
    atom_count = 0
    hetatm_count = 0
    b_factors = []
    title_lines = []

    with open(path) as handle:
        for line in handle:
            rec = line[:6].strip()
            if rec == "TITLE":
                title_lines.append(line[10:].strip())
            elif rec in {"ATOM", "HETATM"}:
                if len(line) < 66:
                    continue
                try:
                    b_factor = float(line[60:66].strip())
                except ValueError:
                    continue
                chain = line[21].strip() or "?"
                residue_name = line[17:20].strip()
                residue_id = line[22:26].strip()
                chains.add(chain)
                residues[chain].add((residue_name, residue_id))
                b_factors.append(b_factor)
                if rec == "ATOM":
                    atom_count += 1
                else:
                    hetatm_count += 1

    return _build_output(
        path=path,
        fmt="pdb",
        title=" ".join(title_lines) if title_lines else None,
        chains=chains,
        residues=residues,
        atom_count=atom_count,
        hetatm_count=hetatm_count,
        b_factors=b_factors,
        force_alphafold=force_alphafold,
    )


def _tokenize_cif_row(line: str) -> list[str]:
    return shlex.split(line, posix=True)


def _simple_cif_value(lines: list[str], key: str) -> str | None:
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith(key):
            rest = stripped[len(key):].strip()
            if rest:
                return rest.strip("'\"")
            if index + 1 < len(lines) and lines[index + 1].startswith(";"):
                values = []
                for text_line in lines[index + 2:]:
                    if text_line.startswith(";"):
                        break
                    values.append(text_line.rstrip())
                return " ".join(values).strip()
    return None


def parse_mmcif(path: str, force_alphafold: bool = False) -> dict:
    lines = Path(path).read_text(errors="replace").splitlines()
    chains = set()
    residues = defaultdict(set)
    atom_count = 0
    hetatm_count = 0
    b_factors = []

    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue

        index += 1
        headers = []
        while index < len(lines) and lines[index].strip().startswith("_"):
            headers.append(lines[index].strip())
            index += 1

        if not headers or not all(header.startswith("_atom_site.") for header in headers):
            while index < len(lines) and lines[index].strip() and not lines[index].strip().startswith(("loop_", "_", "#")):
                index += 1
            continue

        field_index = {header: pos for pos, header in enumerate(headers)}
        group_i = field_index.get("_atom_site.group_PDB")
        chain_i = field_index.get("_atom_site.auth_asym_id", field_index.get("_atom_site.label_asym_id"))
        comp_i = field_index.get("_atom_site.label_comp_id", field_index.get("_atom_site.auth_comp_id"))
        seq_i = field_index.get("_atom_site.auth_seq_id", field_index.get("_atom_site.label_seq_id"))
        b_i = field_index.get("_atom_site.B_iso_or_equiv")

        while index < len(lines):
            stripped = lines[index].strip()
            if not stripped or stripped.startswith("#"):
                index += 1
                break
            if stripped == "loop_" or stripped.startswith("_"):
                break
            try:
                row = _tokenize_cif_row(stripped)
            except ValueError:
                index += 1
                continue
            index += 1
            if len(row) < len(headers):
                continue

            group = row[group_i] if group_i is not None else "ATOM"
            chain = row[chain_i] if chain_i is not None else "?"
            comp = row[comp_i] if comp_i is not None else "?"
            seq = row[seq_i] if seq_i is not None else "?"
            chains.add(chain)
            residues[chain].add((comp, seq))
            if group == "ATOM":
                atom_count += 1
            elif group == "HETATM":
                hetatm_count += 1
            if b_i is not None:
                try:
                    b_factors.append(float(row[b_i]))
                except ValueError:
                    pass

    return _build_output(
        path=path,
        fmt="mmcif",
        title=_simple_cif_value(lines, "_struct.title"),
        chains=chains,
        residues=residues,
        atom_count=atom_count,
        hetatm_count=hetatm_count,
        b_factors=b_factors,
        force_alphafold=force_alphafold,
    )


def _build_output(path: str, fmt: str, title: str | None, chains: set[str],
                  residues: dict[str, set[tuple[str, str]]], atom_count: int,
                  hetatm_count: int, b_factors: list[float],
                  force_alphafold: bool) -> dict:
    sorted_chains = sorted(chains)
    chain_details = {
        chain: {"residues": len(residues[chain])}
        for chain in sorted_chains
    }
    likely_alphafold = _is_likely_alphafold(path, b_factors, force_alphafold)
    output = {
        "status": "ok",
        "file": str(Path(path).resolve()),
        "format": fmt,
        "title": title or "(no title)",
        "chains": sorted_chains,
        "chain_details": chain_details,
        "atom_records": atom_count,
        "hetatm_records": hetatm_count,
        "bfactor": _bfactor_stats(b_factors),
        "likely_alphafold": likely_alphafold,
    }
    if likely_alphafold and b_factors:
        output["plddt_distribution"] = _plddt_distribution(b_factors)
    return output


def inspect_structure(path: str, force_alphafold: bool = False) -> dict:
    suffix = Path(path).suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return parse_mmcif(path, force_alphafold)
    return parse_pdb(path, force_alphafold)


def main():
    parser = argparse.ArgumentParser(
        description="Inspect PDB or mmCIF coordinate files. Zero dependencies.",
        epilog=(
            "Examples:\n"
            "  %(prog)s structure.pdb\n"
            "  %(prog)s structure.cif --json\n"
            "  %(prog)s AF-P04637-F1.cif --alphafold --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("structure", help="Path to .pdb, .cif, or .mmcif file")
    parser.add_argument("--json", action="store_true", dest="as_json", help="Output as JSON")
    parser.add_argument("--alphafold", action="store_true", help="Treat B-factors as pLDDT confidence")
    args = parser.parse_args()

    if not Path(args.structure).exists():
        message = f"File not found: {args.structure}"
        if args.as_json:
            print(json.dumps({"status": "error", "error": message}, indent=2))
        else:
            print(f"ERROR: {message}", file=sys.stderr)
        sys.exit(1)

    output = inspect_structure(args.structure, args.alphafold)
    if args.as_json:
        print(json.dumps(output, indent=2))
    else:
        print(f"File: {output['file']}")
        print(f"Format: {output['format']}")
        print(f"Title: {output['title']}")
        print(f"Chains: {output['chains']}")
        print(f"ATOM records: {output['atom_records']}")
        print(f"HETATM records: {output['hetatm_records']}")
        for chain in output["chains"]:
            print(f"  Chain {chain}: {output['chain_details'][chain]['residues']} residues")
        if output["bfactor"]:
            stats = output["bfactor"]
            print(f"B-factor range: {stats['min']} - {stats['max']} (mean {stats['mean']})")
        if output["likely_alphafold"]:
            plddt = output["plddt_distribution"]
            print("Likely AlphaFold pLDDT:")
            print(f"  >90: {plddt['very_high_gt90']:.1%}")
            print(f"  70-90: {plddt['confident_70_90']:.1%}")
            print(f"  50-70: {plddt['low_50_70']:.1%}")
            print(f"  <=50: {plddt['very_low_le50']:.1%}")


if __name__ == "__main__":
    main()
