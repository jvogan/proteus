#!/usr/bin/env python3
"""Fetch experimental structures and metadata from RCSB PDB.

Uses only Python stdlib. Designed for agent workflows where a user gives a
PDB ID such as 1HSG or 4HHB and expects a local structure file plus summary
metadata.

Usage:
    python fetch_pdb.py 4HHB
    python fetch_pdb.py 1HSG --format pdb
    python fetch_pdb.py 4HHB --assembly 1
    python fetch_pdb.py 4HHB --metadata --json
"""

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_DOWNLOAD = "https://files.rcsb.org/download/{filename}"
RCSB_HEADER = "https://files.rcsb.org/header/{filename}"
RCSB_BINARY_CIF = "https://models.rcsb.org/{pdb_id}.bcif"


def _pdb_id(value: str) -> str:
    pdb_id = value.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise argparse.ArgumentTypeError("PDB ID must be exactly 4 letters/digits, e.g. 4HHB")
    return pdb_id


def _log(message: str, *, as_json: bool = False):
    print(message, file=sys.stderr if as_json else sys.stdout)


def _fetch_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=30) as response:
        return json.load(response)


def fetch_entry_metadata(pdb_id: str) -> dict:
    """Fetch selected entry metadata from RCSB Data API."""
    try:
        data = _fetch_json(RCSB_DATA_API.format(pdb_id=pdb_id))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise SystemExit(f"ERROR: PDB ID '{pdb_id}' was not found in RCSB PDB.") from exc
        raise

    info = data.get("rcsb_entry_info", {})
    title = data.get("struct", {}).get("title")
    methods = [item.get("method") for item in data.get("exptl", []) if item.get("method")]
    return {
        "pdb_id": pdb_id,
        "title": title,
        "experimental_methods": methods,
        "resolution": info.get("resolution_combined"),
        "assembly_count": info.get("assembly_count"),
        "polymer_entity_count": info.get("polymer_entity_count"),
        "nonpolymer_bound_components": info.get("nonpolymer_bound_components", []),
        "deposited_atom_count": info.get("deposited_atom_count"),
        "deposited_model_count": info.get("deposited_model_count"),
        "molecular_weight_kda": info.get("molecular_weight"),
        "selected_polymer_entity_types": info.get("selected_polymer_entity_types"),
    }


def build_download_url(pdb_id: str, fmt: str, assembly: int | None, header_only: bool) -> tuple[str, str]:
    """Return (url, filename) for an RCSB downloadable structure artifact."""
    if assembly is not None and fmt != "cif":
        raise SystemExit("ERROR: --assembly currently supports --format cif only.")
    if header_only and assembly is not None:
        raise SystemExit("ERROR: --header-only cannot be combined with --assembly.")
    if header_only and fmt not in {"cif", "xml"}:
        raise SystemExit("ERROR: --header-only supports --format cif or xml only.")

    if assembly is not None:
        filename = f"{pdb_id}-assembly{assembly}.cif"
        return RCSB_DOWNLOAD.format(filename=filename), filename

    filename = f"{pdb_id}.{fmt}"
    if fmt == "bcif":
        if header_only:
            raise SystemExit("ERROR: --header-only does not support --format bcif.")
        return RCSB_BINARY_CIF.format(pdb_id=pdb_id.lower()), filename
    if header_only:
        return RCSB_HEADER.format(filename=filename), filename
    return RCSB_DOWNLOAD.format(filename=filename), filename


def download(url: str, destination: Path):
    """Download URL to destination, converting HTTP 404s into useful errors."""
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            destination.write_bytes(response.read())
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise SystemExit(f"ERROR: RCSB file not found: {url}") from exc
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Fetch PDB structures and selected RCSB metadata.",
        epilog=(
            "Examples:\n"
            "  %(prog)s 4HHB\n"
            "  %(prog)s 1HSG --format pdb\n"
            "  %(prog)s 4HHB --assembly 1 --json\n"
            "  %(prog)s 4HHB --metadata --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("pdb_id", type=_pdb_id, help="Four-character PDB ID, e.g. 4HHB")
    parser.add_argument("--format", choices=["cif", "pdb", "xml", "bcif"], default="cif",
                        help="Download format (default: cif)")
    parser.add_argument("--assembly", type=int, help="Download biological assembly N (mmCIF only)")
    parser.add_argument("--header-only", action="store_true", help="Download header/summary coordinates only")
    parser.add_argument("--metadata", action="store_true", help="Fetch metadata only; do not download coordinates")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON to stdout")
    args = parser.parse_args()

    metadata = fetch_entry_metadata(args.pdb_id)
    output = {
        "status": "ok",
        "source": "RCSB PDB",
        "metadata": metadata,
        "download": None,
    }

    if not args.metadata:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        url, filename = build_download_url(args.pdb_id, args.format, args.assembly, args.header_only)
        path = outdir / filename
        _log(f"Downloading {args.pdb_id} from RCSB: {url}", as_json=args.json)
        download(url, path)
        output["download"] = {
            "url": url,
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "format": args.format,
            "assembly": args.assembly,
            "header_only": args.header_only,
        }

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print(f"PDB ID: {metadata['pdb_id']}")
        print(f"Title: {metadata.get('title') or '(unknown)'}")
        methods = ", ".join(metadata.get("experimental_methods") or []) or "(unknown)"
        print(f"Method: {methods}")
        if metadata.get("resolution"):
            print(f"Resolution: {metadata['resolution']}")
        if metadata.get("nonpolymer_bound_components"):
            print(f"Ligands/components: {', '.join(metadata['nonpolymer_bound_components'])}")
        if output["download"]:
            print(f"Saved: {output['download']['path']} ({output['download']['bytes']:,} bytes)")


if __name__ == "__main__":
    main()
