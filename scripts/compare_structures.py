#!/usr/bin/env python3
"""Compare two structures with PyMOL cealign and optional CA deviations.

Usage:
    python3 compare_structures.py reference.pdb mobile.pdb --json
    python3 compare_structures.py reference.pdb mobile.pdb --per-residue --json
"""

import argparse
import json
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import pymol_agent


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def compare(reference: str, mobile: str, per_residue: bool = False) -> dict:
    ref_path = os.path.abspath(reference)
    mob_path = os.path.abspath(mobile)
    script = f'''
import math

reference_path = {ref_path!r}
mobile_path = {mob_path!r}
cmd.load(reference_path, "reference")
cmd.load(mobile_path, "mobile")
try:
    result = cmd.cealign("reference", "mobile")
    method = "cealign"
except Exception as cealign_error:
    raw = cmd.align("mobile", "reference")
    result = {{
        "RMSD": raw[0] if raw else None,
        "alignment_length": raw[1] if len(raw) > 1 else None,
        "raw": raw,
        "cealign_error": str(cealign_error),
    }}
    method = "align"
_output["data"]["reference"] = reference_path
_output["data"]["mobile"] = mobile_path
_output["data"]["method"] = method
_output["data"]["alignment"] = result

if {bool(per_residue)!r}:
    ref_coords = {{}}
    mob_coords = {{}}
    cmd.iterate_state(1, "reference and name CA",
        "ref_coords[(chain, resi)] = (x, y, z)", space={{"ref_coords": ref_coords}})
    cmd.iterate_state(1, "mobile and name CA",
        "mob_coords[(chain, resi)] = (x, y, z)", space={{"mob_coords": mob_coords}})
    deviations = []
    for key, ref_xyz in ref_coords.items():
        mob_xyz = mob_coords.get(key)
        if mob_xyz is None:
            continue
        dx = ref_xyz[0] - mob_xyz[0]
        dy = ref_xyz[1] - mob_xyz[1]
        dz = ref_xyz[2] - mob_xyz[2]
        deviations.append({{
            "chain": key[0],
            "resi": key[1],
            "ca_distance": round(math.sqrt(dx*dx + dy*dy + dz*dz), 3),
        }})
    deviations.sort(key=lambda item: item["ca_distance"], reverse=True)
    _output["data"]["per_residue"] = {{
        "matched_ca_count": len(deviations),
        "mean_ca_distance": round(sum(item["ca_distance"] for item in deviations) / len(deviations), 3) if deviations else None,
        "max_ca_distance": deviations[0]["ca_distance"] if deviations else None,
        "largest_deviations": deviations[:20],
    }}
'''
    result = pymol_agent.run_pymol_script(script, timeout=300)
    if result.get("status") != "ok":
        return result
    return _ok_payload(result.get("data", {}))


def main():
    parser = argparse.ArgumentParser(description="Compare two structures with PyMOL cealign.")
    parser.add_argument("reference", help="Reference structure path (stays fixed)")
    parser.add_argument("mobile", help="Mobile structure path (gets aligned)")
    parser.add_argument("--per-residue", action="store_true", help="Also report matched CA deviations")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    for path in (args.reference, args.mobile):
        if not Path(path).exists():
            payload = _error_payload(f"File not found: {path}")
            if args.json:
                print(json.dumps(payload, indent=2))
            else:
                print(f"ERROR: {payload['error']}", file=sys.stderr)
            sys.exit(1)

    output = compare(args.reference, args.mobile, args.per_residue)
    if args.json:
        print(json.dumps(output, indent=2))
    elif output.get("status") != "ok":
        print(f"ERROR: {output.get('error', 'comparison failed')}", file=sys.stderr)
        sys.exit(1)
    else:
        alignment = output["data"].get("alignment")
        print(f"Alignment: {alignment}")
        if output["data"].get("per_residue"):
            per = output["data"]["per_residue"]
            print(f"Matched CA atoms: {per['matched_ca_count']}")
            print(f"Mean CA distance: {per['mean_ca_distance']}")
            print(f"Max CA distance: {per['max_ca_distance']}")


if __name__ == "__main__":
    main()
