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

import fetch_alphafold
import fetch_pdb
import structure_info
import uniprot_lookup


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


def _download_pdb(pdb_id: str, outdir: str, download: bool) -> dict:
    metadata = fetch_pdb.fetch_entry_metadata(pdb_id)
    result = {
        "resolved_kind": "pdb_id",
        "source": "RCSB PDB",
        "query": pdb_id,
        "metadata": metadata,
        "structure_path": None,
        "inspection": None,
    }
    if download:
        output_dir = Path(outdir)
        output_dir.mkdir(parents=True, exist_ok=True)
        url, filename = fetch_pdb.build_download_url(pdb_id, "cif", None, False)
        destination = output_dir / filename
        fetch_pdb.download(url, destination)
        result["structure_path"] = str(destination.resolve())
        result["inspection"] = structure_info.inspect_structure(str(destination))
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
        result["structure_path"] = str(destination.resolve())
        result["inspection"] = structure_info.inspect_structure(str(destination), force_alphafold=True)
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
            "query": query,
            "structure_path": str(path.resolve()),
            "inspection": structure_info.inspect_structure(str(path)),
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
        if data.get("global_plddt") is not None:
            print(f"Global pLDDT: {data['global_plddt']}")


if __name__ == "__main__":
    main()
