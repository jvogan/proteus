#!/usr/bin/env python3
"""Fetch wwPDB/RCSB validation summary metrics for a PDB entry.

Usage:
    python3 validation_report.py 4HHB --json
"""

import argparse
import json
import re
import sys
import urllib.error
import urllib.request


RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


class ValidationReportError(RuntimeError):
    pass


def _pdb_id(value: str) -> str:
    pdb_id = value.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise argparse.ArgumentTypeError("PDB ID must be exactly 4 letters/digits, e.g. 4HHB")
    return pdb_id


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _fetch_entry(pdb_id: str) -> dict:
    url = RCSB_ENTRY.format(pdb_id=pdb_id)
    request = urllib.request.Request(url, headers={"User-Agent": "proteus-skill/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise ValidationReportError(f"PDB ID '{pdb_id}' was not found in RCSB PDB.") from exc
        raise ValidationReportError(f"RCSB returned HTTP {exc.code}: {exc.reason}") from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise ValidationReportError(f"Failed to fetch validation data: {exc}") from exc


def _first(items) -> dict:
    if isinstance(items, list) and items:
        return items[0] or {}
    if isinstance(items, dict):
        return items
    return {}


def fetch_validation_report(pdb_id: str) -> dict:
    entry = _fetch_entry(pdb_id)
    geometry = _first(entry.get("pdbx_vrpt_summary_geometry"))
    summary = entry.get("pdbx_vrpt_summary") or {}
    entry_info = entry.get("rcsb_entry_info") or {}
    data = {
        "pdb_id": pdb_id,
        "title": (entry.get("struct") or {}).get("title"),
        "experimental_methods": [item.get("method") for item in entry.get("exptl", []) if item.get("method")],
        "resolution": entry_info.get("resolution_combined"),
        "report": {
            "created": summary.get("report_creation_date"),
            "attempted_steps": summary.get("attempted_validation_steps"),
            "ligands_for_buster_report": summary.get("ligands_for_buster_report"),
        },
        "geometry": {
            "clashscore": geometry.get("clashscore"),
            "bond_rmsz": geometry.get("bonds_RMSZ"),
            "angle_rmsz": geometry.get("angles_RMSZ"),
            "ramachandran_outliers_percent": geometry.get("percent_ramachandran_outliers"),
            "rotamer_outliers_percent": geometry.get("percent_rotamer_outliers"),
        },
        "raw_available": {
            "pdbx_vrpt_summary": bool(summary),
            "pdbx_vrpt_summary_geometry": bool(geometry),
        },
    }
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(description="Fetch wwPDB/RCSB validation summary metrics.")
    parser.add_argument("pdb_id", type=_pdb_id, help="Four-character PDB ID, e.g. 4HHB")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        output = fetch_validation_report(args.pdb_id)
    except ValidationReportError as exc:
        if args.json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        data = output["data"]
        print(f"PDB ID: {data['pdb_id']}")
        print(f"Title: {data.get('title') or '(unknown)'}")
        print(f"Report created: {data['report'].get('created') or '(unknown)'}")
        geometry = data["geometry"]
        print(f"Clashscore: {geometry.get('clashscore')}")
        print(f"Bond RMSZ: {geometry.get('bond_rmsz')}")
        print(f"Angle RMSZ: {geometry.get('angle_rmsz')}")
        print(f"Ramachandran outliers: {geometry.get('ramachandran_outliers_percent')}%")
        print(f"Rotamer outliers: {geometry.get('rotamer_outliers_percent')}%")


if __name__ == "__main__":
    main()
