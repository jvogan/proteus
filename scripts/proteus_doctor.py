#!/usr/bin/env python3
"""Report local Proteus capability status.

This is the first command to run when an agent needs to know which structural
biology workflows are available on a machine.

Usage:
    python3 proteus_doctor.py --json
    python3 proteus_doctor.py --network --json
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NETWORK_ENDPOINTS = {
    "rcsb": "https://data.rcsb.org/rest/v1/core/entry/4HHB",
    "uniprot": "https://rest.uniprot.org/uniprotkb/search?query=accession:P04637&format=json&size=1",
    "alphafold": "https://alphafold.ebi.ac.uk/api/prediction/P04637",
}
NETWORK_CAPABILITIES = {
    "rcsb_fetch": "rcsb",
    "uniprot_lookup": "uniprot",
    "alphafold_fetch": "alphafold",
}


def _display_path(path: str | Path | None) -> str | None:
    if path is None:
        return None
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


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _find_pymol() -> str | None:
    found = shutil.which("pymol")
    if found:
        return found
    for path in [
        "/Applications/PyMOL.app/Contents/bin/pymol",
        os.path.expanduser("~/Applications/PyMOL.app/Contents/bin/pymol"),
        "/usr/bin/pymol",
        "/usr/local/bin/pymol",
    ]:
        if os.path.isfile(path):
            return path
    return None


def _find_chimerax() -> str | None:
    found = shutil.which("ChimeraX") or shutil.which("chimerax")
    if found:
        return found
    hits = sorted(Path("/Applications").glob("ChimeraX*.app/Contents/bin/ChimeraX"))
    if hits:
        return str(hits[-1])
    for path in [
        "/usr/bin/chimerax",
        "/usr/local/bin/chimerax",
        os.path.expanduser("~/ChimeraX/bin/ChimeraX"),
    ]:
        if os.path.isfile(path):
            return path
    return None


def _python_version() -> dict:
    return {
        "executable": _display_path(sys.executable),
        "version": platform.python_version(),
        "ok": sys.version_info >= (3, 10),
    }


def _run_quick(cmd: list[str], timeout: int = 8) -> dict:
    try:
        proc = subprocess.run(cmd, text=True, capture_output=True, timeout=timeout)
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    return {
        "ok": proc.returncode == 0,
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip()[:1000],
        "stderr": proc.stderr.strip()[:1000],
    }


def _script_smoke() -> dict:
    scripts = [
        "fetch_pdb.py",
        "pdb_search.py",
        "pdb_select.py",
        "assembly_report.py",
        "sifts_map.py",
        "ligand_extract.py",
        "dock_prep.py",
        "docking_box.py",
        "dock_vina.py",
        "interaction_report.py",
        "variant_map.py",
        "mutation_triage.py",
        "proteus_batch.py",
        "proteus_cache.py",
        "proteus_report.py",
        "target_dossier.py",
        "uniprot_lookup.py",
        "structure_info.py",
        "fetch_alphafold.py",
        "pdb_info.py",
        "pymol_agent.py",
        "chimerax_agent.py",
        "chimerax_rest.py",
        "pae_report.py",
        "resolve_structure.py",
        "validation_report.py",
        "pocket_report.py",
        "compare_structures.py",
        "add_helix_records.py",
        "interface_report.py",
        "map_info.py",
        "model_quality.py",
        "design_run.py",
        "rosetta_score.py",
        "kras_dossier.py",
    ]
    if sys.version_info < (3, 10):
        return {
            script: {"ok": False, "error": "skipped: Python 3.10+ required for helper scripts"}
            for script in scripts
        }
    results = {}
    for script in scripts:
        path = ROOT / "scripts" / script
        if not path.exists():
            results[script] = {"ok": False, "error": "missing"}
            continue
        results[script] = _run_quick([sys.executable, str(path), "--help"])
    return results


def _network_check(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "proteus-skill/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=10) as response:
            response.read(256)
        return {"ok": True}
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        return {"ok": False, "error": str(exc)}


def build_report(include_network: bool) -> dict:
    pymol = _find_pymol()
    chimerax = _find_chimerax()
    ffmpeg = shutil.which("ffmpeg")
    network = None
    if include_network:
        network = {name: _network_check(url) for name, url in NETWORK_ENDPOINTS.items()}
    capabilities = {
        "zero_dependency_inspection": True,
        "rcsb_fetch": False,
        "uniprot_lookup": False,
        "alphafold_fetch": False,
        "pymol_rendering": bool(pymol),
        "chimerax_analysis": bool(chimerax),
        "turntable_movies": bool((pymol or chimerax) and ffmpeg),
    }
    if network is not None:
        for capability, check_name in NETWORK_CAPABILITIES.items():
            capabilities[capability] = bool(network[check_name].get("ok"))
    data = {
        "root": _display_path(ROOT),
        "platform": {
            "system": platform.system(),
            "machine": platform.machine(),
        },
        "python": _python_version(),
        "tools": {
            "pymol": {"ok": bool(pymol), "path": _display_path(pymol)},
            "chimerax": {"ok": bool(chimerax), "path": _display_path(chimerax)},
            "ffmpeg": {"ok": bool(ffmpeg), "path": _display_path(ffmpeg)},
        },
        "scripts": _script_smoke(),
        "network": network,
        "capabilities": capabilities,
    }
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(
        description="Report Proteus local tool, script, and optional network readiness."
    )
    parser.add_argument("--network", action="store_true", help="Check RCSB, UniProt, and AlphaFold endpoints")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    report = build_report(args.network)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        data = report["data"]
        print(f"Proteus root: {data['root']}")
        print(f"Python: {data['python']['version']} ({'ok' if data['python']['ok'] else 'too old'})")
        print(f"PyMOL: {data['tools']['pymol']['path'] or 'not found'}")
        print(f"ChimeraX: {data['tools']['chimerax']['path'] or 'not found'}")
        print(f"ffmpeg: {data['tools']['ffmpeg']['path'] or 'not found'} (turntable movies)")
        failed = [name for name, result in data["scripts"].items() if not result["ok"]]
        print(f"Script help checks: {'ok' if not failed else 'failed: ' + ', '.join(failed)}")
        if data["network"]:
            network_failed = [name for name, result in data["network"].items() if not result["ok"]]
            print(f"Network checks: {'ok' if not network_failed else 'failed: ' + ', '.join(network_failed)}")


if __name__ == "__main__":
    main()
