#!/usr/bin/env python3
"""Report biological assembly availability for an RCSB PDB entry.

Uses only Python stdlib. This helper intentionally does not reconstruct or
apply assembly transforms; it reports RCSB assembly files and can download the
selected assembly mmCIF as provided by RCSB.

Usage:
    python assembly_report.py 4HHB --json
    python assembly_report.py 4HHB --download --outdir structures --json
    python assembly_report.py 1HSG --assembly 1 --download
"""

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_DOWNLOAD = "https://files.rcsb.org/download/{filename}"
USER_AGENT = "proteus-skill/1.0"


class AssemblyReportError(RuntimeError):
    pass


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _pdb_id(value: str) -> str:
    pdb_id = value.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise argparse.ArgumentTypeError("PDB ID must be exactly 4 letters/digits, e.g. 4HHB")
    return pdb_id


def _display_path(path: str | Path) -> str:
    candidate = Path(path).expanduser()
    try:
        resolved = candidate.resolve(strict=False)
    except OSError:
        resolved = candidate.absolute()
    for root in (ROOT, Path.cwd()):
        try:
            relative = resolved.relative_to(root.resolve())
            return f"./{relative}" if str(relative) != "." else "."
        except ValueError:
            continue
    return f"{resolved.name} (absolute path omitted)"


def _fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise AssemblyReportError("RCSB entry was not found.") from exc
        raise AssemblyReportError(f"RCSB API returned HTTP {exc.code}: {exc.reason}") from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise AssemblyReportError(f"Failed to fetch RCSB metadata: {exc}") from exc


def _download(url: str, destination: Path):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            destination.write_bytes(response.read())
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise AssemblyReportError(f"RCSB assembly file was not found: {url}") from exc
        raise AssemblyReportError(f"RCSB download returned HTTP {exc.code}: {exc.reason}") from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise AssemblyReportError(f"Failed to download RCSB assembly file {url}: {exc}") from exc


def build_assembly_downloads(pdb_id: str, assembly_count: int | None) -> list[dict]:
    """Return deterministic RCSB assembly mmCIF downloads for an entry."""
    if not assembly_count:
        return []
    downloads = []
    for assembly_id in range(1, assembly_count + 1):
        filename = f"{pdb_id}-assembly{assembly_id}.cif"
        downloads.append({
            "assembly_id": str(assembly_id),
            "format": "mmcif",
            "filename": filename,
            "url": RCSB_DOWNLOAD.format(filename=filename),
        })
    return downloads


def select_assembly(assemblies: list[dict], requested: int | None = None) -> dict | None:
    """Pick requested assembly, or assembly 1 when assemblies are available."""
    if not assemblies:
        return None
    assembly_id = str(requested or 1)
    for assembly in assemblies:
        if assembly["assembly_id"] == assembly_id:
            return assembly
    raise AssemblyReportError(f"Assembly {assembly_id} is not available for this entry.")


def build_report(pdb_id: str, entry: dict, requested_assembly: int | None = None) -> dict:
    info = entry.get("rcsb_entry_info") or {}
    assembly_count = info.get("assembly_count")
    assemblies = build_assembly_downloads(pdb_id, assembly_count)
    recommended = select_assembly(assemblies)
    selected = select_assembly(assemblies, requested_assembly) if requested_assembly else recommended
    data = {
        "source": "RCSB PDB",
        "pdb_id": pdb_id,
        "title": (entry.get("struct") or {}).get("title"),
        "assembly_count": assembly_count,
        "polymer_entity_count": info.get("polymer_entity_count"),
        "deposited_atom_count": info.get("deposited_atom_count"),
        "deposited_model_count": info.get("deposited_model_count"),
        "assembly_downloads": assemblies,
        "recommended_assembly": recommended,
        "selected_assembly": selected,
        "download": None,
        "notes": [
            "Assembly files are reported from RCSB download artifacts; no local assembly transforms were applied."
        ],
    }
    return _ok_payload(data)


def fetch_assembly_report(pdb_id: str, requested_assembly: int | None = None) -> dict:
    try:
        entry = _fetch_json(RCSB_ENTRY.format(pdb_id=pdb_id))
    except AssemblyReportError as exc:
        if str(exc) == "RCSB entry was not found.":
            raise AssemblyReportError(f"PDB ID '{pdb_id}' was not found in RCSB PDB.") from exc
        raise
    return build_report(pdb_id, entry, requested_assembly)


def download_selected_assembly(report: dict, outdir: str) -> dict:
    data = report["data"]
    selected = data.get("selected_assembly")
    if not selected:
        raise AssemblyReportError("No biological assembly download is available for this entry.")

    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / selected["filename"]
    _download(selected["url"], destination)
    data["download"] = {
        "assembly_id": selected["assembly_id"],
        "format": selected["format"],
        "url": selected["url"],
        "filename": selected["filename"],
        "path": _display_path(destination),
        "bytes": destination.stat().st_size,
    }
    report["download"] = data["download"]
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Report RCSB biological assembly availability and optionally download assembly mmCIF.",
        epilog=(
            "Examples:\n"
            "  %(prog)s 4HHB --json\n"
            "  %(prog)s 4HHB --download --outdir structures --json\n"
            "  %(prog)s 1HSG --assembly 1 --download"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("pdb_id", type=_pdb_id, help="Four-character PDB ID, e.g. 4HHB")
    parser.add_argument("--assembly", type=int, help="Assembly number to recommend/download (default: 1)")
    parser.add_argument("--download", action="store_true", help="Download selected biological assembly mmCIF")
    parser.add_argument("--outdir", default=".", help="Output directory for --download (default: current)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        if args.assembly is not None and args.assembly < 1:
            raise AssemblyReportError("--assembly must be a positive integer.")
        output = fetch_assembly_report(args.pdb_id, args.assembly)
        if args.download:
            output = download_selected_assembly(output, args.outdir)
    except (AssemblyReportError, OSError) as exc:
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
        print(f"Assemblies: {data.get('assembly_count') or 0}")
        print(f"Polymer entities: {data.get('polymer_entity_count')}")
        print(f"Deposited atoms: {data.get('deposited_atom_count')}")
        print(f"Deposited models: {data.get('deposited_model_count')}")
        if data["assembly_downloads"]:
            print("Assembly mmCIF downloads:")
            for assembly in data["assembly_downloads"]:
                markers = []
                if assembly == data["recommended_assembly"]:
                    markers.append("recommended")
                if assembly == data["selected_assembly"] and assembly != data["recommended_assembly"]:
                    markers.append("selected")
                marker = f" ({', '.join(markers)})" if markers else ""
                print(f"  assembly {assembly['assembly_id']}: {assembly['filename']}{marker}")
                print(f"    {assembly['url']}")
        else:
            print("Assembly mmCIF downloads: none reported")
        if data["download"]:
            print(f"Saved: {data['download']['path']} ({data['download']['bytes']:,} bytes)")


if __name__ == "__main__":
    main()
