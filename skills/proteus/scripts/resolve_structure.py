#!/usr/bin/env python3
"""Resolve a protein query to a local structure or structure metadata.

Accepted inputs:
  - local .pdb/.cif/.mmcif path
  - PDB ID, e.g. 4HHB
  - UniProt accession, e.g. P04637
  - gene/protein name, e.g. TP53

Usage:
    python3 resolve_structure.py TP53 --json
    python3 resolve_structure.py 4HHB --source pdb --json
    python3 resolve_structure.py P04637 --source alphafold --no-download --json
"""

import argparse
import json
import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import assembly_report
import fetch_alphafold
import fetch_pdb
import structure_info
import uniprot_lookup


ROOT = Path(__file__).resolve().parents[1]
ASYMMETRIC_UNIT_WARNING = (
    "Downloaded RCSB default mmCIF coordinates are the deposited asymmetric unit "
    "and may not represent the biological assembly."
)


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _looks_like_pdb_id(value: str) -> bool:
    # Legacy four-character PDB IDs start with a digit; this avoids treating
    # common gene symbols such as TP53 as PDB IDs in auto mode.
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()))


def _looks_like_uniprot(value: str) -> bool:
    return bool(re.fullmatch(
        r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?|[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?",
        value.strip().upper(),
    ))


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


def _safe_inspection(path: str | Path, *, force_alphafold: bool = False) -> dict:
    inspection = structure_info.inspect_structure(str(path), force_alphafold=force_alphafold)
    sanitized = dict(inspection)
    if sanitized.get("file"):
        sanitized["file"] = _display_path(sanitized["file"])
    if isinstance(sanitized.get("data"), dict):
        data = dict(sanitized["data"])
        if data.get("file"):
            data["file"] = _display_path(data["file"])
        sanitized["data"] = data
    return sanitized


def _assembly_entry_from_metadata(metadata: dict) -> dict:
    return {
        "struct": {"title": metadata.get("title")},
        "rcsb_entry_info": {
            "assembly_count": metadata.get("assembly_count"),
            "polymer_entity_count": metadata.get("polymer_entity_count"),
            "deposited_atom_count": metadata.get("deposited_atom_count"),
            "deposited_model_count": metadata.get("deposited_model_count"),
        },
    }


def _biological_assembly_metadata(pdb_id: str, metadata: dict) -> dict:
    report = assembly_report.build_report(pdb_id, _assembly_entry_from_metadata(metadata))
    data = dict(report["data"])
    data.pop("download", None)
    return data


def _assembly_warning(assembly: dict) -> str | None:
    recommended = assembly.get("recommended_assembly")
    if not recommended:
        return None
    filename = recommended.get("filename")
    if filename:
        return f"{ASYMMETRIC_UNIT_WARNING} RCSB biological assembly mmCIF is available as {filename}."
    return ASYMMETRIC_UNIT_WARNING


def _download_pdb(pdb_id: str, outdir: str, download: bool) -> dict:
    metadata = fetch_pdb.fetch_entry_metadata(pdb_id)
    assembly = _biological_assembly_metadata(pdb_id, metadata)
    result = {
        "resolved_kind": "pdb_id",
        "source": "RCSB PDB",
        "query": pdb_id,
        "metadata": metadata,
        "biological_assembly": assembly,
        "download": None,
        "structure_path": None,
        "inspection": None,
        "warnings": [],
    }
    if download:
        output_dir = Path(outdir)
        output_dir.mkdir(parents=True, exist_ok=True)
        url, filename = fetch_pdb.build_download_url(pdb_id, "cif", None, False)
        destination = output_dir / filename
        fetch_pdb.download(url, destination)
        safe_path = _display_path(destination)
        result["structure_path"] = safe_path
        result["download"] = {
            "url": url,
            "filename": filename,
            "path": safe_path,
            "format": "mmcif",
            "assembly": None,
            "coordinate_scope": "asymmetric_unit",
            "bytes": destination.stat().st_size,
        }
        result["inspection"] = _safe_inspection(destination)
        warning = _assembly_warning(assembly)
        if warning:
            result["warnings"].append(warning)
    return result


def _download_alphafold(uniprot_id: str, outdir: str, download: bool) -> dict:
    meta = fetch_alphafold.fetch_metadata(uniprot_id)
    model_id = meta.get("modelEntityId", f"AF-{uniprot_id}-F1")
    result = {
        "resolved_kind": "uniprot_accession",
        "source": "AlphaFold DB",
        "query": uniprot_id,
        "uniprot_id": uniprot_id,
        "model_id": model_id,
        "gene": meta.get("gene"),
        "global_plddt": meta.get("globalMetricValue"),
        "latest_version": meta.get("latestVersion"),
        "structure_path": None,
        "inspection": None,
    }
    if download:
        output_dir = Path(outdir)
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / f"{model_id}.pdb"
        fetch_alphafold.download(meta["pdbUrl"], destination, as_json=True)
        result["structure_path"] = _display_path(destination)
        result["inspection"] = _safe_inspection(destination, force_alphafold=True)
    return result


def _resolve_name(query: str, organism: str, outdir: str, download: bool) -> dict:
    gene_like = bool(re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]{1,15}", query.strip()))
    uniprot_query = uniprot_lookup.build_query(query, organism, True, gene_exact=gene_like)
    hits = uniprot_lookup.search_uniprot(uniprot_query, 5)
    if not hits and gene_like:
        uniprot_query = uniprot_lookup.build_query(query, organism, True, gene_exact=False)
        hits = uniprot_lookup.search_uniprot(uniprot_query, 5)
    if not hits:
        raise ValueError(f"No reviewed UniProt match found for '{query}' in organism {organism}.")
    best = hits[0]
    accession = best["accession"]
    result = _download_alphafold(accession, outdir, download)
    result["resolved_kind"] = "name_or_gene"
    result["query"] = query
    result["uniprot_match"] = best
    result["candidate_count"] = len(hits)
    return result


def resolve(query: str, source: str = "auto", outdir: str = ".", download: bool = True,
            organism: str = "9606") -> dict:
    path = Path(query)
    if path.exists():
        data = {
            "resolved_kind": "local_file",
            "source": "local",
            "query": _display_path(path),
            "structure_path": _display_path(path),
            "inspection": _safe_inspection(path),
        }
        return _ok_payload(data)

    normalized = query.strip().upper()
    if source == "pdb" or (source == "auto" and _looks_like_pdb_id(query) and not _looks_like_uniprot(query)):
        return _ok_payload(_download_pdb(normalized, outdir, download))
    if source == "alphafold" or (source == "auto" and _looks_like_uniprot(query)):
        return _ok_payload(_download_alphafold(normalized, outdir, download))
    if source != "auto":
        raise ValueError(f"Cannot resolve '{query}' as source '{source}'.")
    return _ok_payload(_resolve_name(query, organism, outdir, download))


def main():
    parser = argparse.ArgumentParser(description="Resolve a local file, PDB ID, UniProt accession, or protein name.")
    parser.add_argument("query", help="Local file, PDB ID, UniProt accession, gene, or protein name")
    parser.add_argument("--source", choices=["auto", "pdb", "alphafold"], default="auto",
                        help="Force a source instead of auto-detecting")
    parser.add_argument("--outdir", default=".", help="Output directory for downloaded structures")
    parser.add_argument("--no-download", action="store_true", help="Return metadata only; do not download coordinates")
    parser.add_argument("--organism", default="9606", help="NCBI taxonomy ID for name/gene lookup (default: human)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        output = resolve(args.query, args.source, args.outdir, not args.no_download, args.organism)
    except (ValueError, KeyError, OSError, fetch_pdb.PDBFetchError,
            fetch_alphafold.AlphaFoldFetchError, uniprot_lookup.UniProtLookupError) as exc:
        if args.json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        data = output["data"]
        print(f"Resolved: {data['resolved_kind']} via {data['source']}")
        if data.get("structure_path"):
            print(f"Structure: {data['structure_path']}")
        assembly = data.get("biological_assembly")
        if assembly:
            print(f"Biological assemblies: {assembly.get('assembly_count') or 0}")
            recommended = assembly.get("recommended_assembly")
            if recommended:
                print(f"Recommended assembly: {recommended['filename']}")
        if data.get("global_plddt") is not None:
            print(f"Global pLDDT: {data['global_plddt']}")
        for warning in data.get("warnings", []):
            print(f"WARNING: {warning}", file=sys.stderr)


if __name__ == "__main__":
    main()
