#!/usr/bin/env python3
"""Search the RCSB PDB for structures and (optionally) enrich each hit.

This is the discovery counterpart to ``resolve_structure.py``: where that maps a
single query to one structure, this enumerates *many* structures matching a
free-text query and/or a UniProt accession, the tedious "click around the PDB
website" step done in one call.

Examples:
    # Free-text relevance search
    python3 pdb_search.py "KRAS G12C sotorasib" --rows 10 --details

    # Every experimental structure mapped to a UniProt accession
    python3 pdb_search.py --uniprot P01116 --rows 50 --json

    # Combine: KRAS structures whose entry text mentions GMPPNP
    python3 pdb_search.py "GMPPNP" --uniprot P01116 --details --json

Output (``--json``) is a ``{"status": "ok", "data": {...}}`` payload listing the
matching PDB IDs, with title / resolution / method / bound-ligand codes when
``--details`` is given.
"""

import argparse
import json
import sys
import urllib.error
import urllib.request

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/"


def _post(url: str, payload: dict, timeout: int = 30) -> dict:
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def _get(url: str, timeout: int = 30) -> dict:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return json.load(resp)


def _full_text_node(value: str) -> dict:
    return {"type": "terminal", "service": "full_text", "parameters": {"value": value}}


def _uniprot_nodes(accession: str) -> list:
    base = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers"
    return [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": f"{base}.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": f"{base}.database_name",
                "operator": "exact_match",
                "value": "UniProt",
            },
        },
    ]


def build_query(text: str | None, uniprot: str | None) -> dict:
    nodes = []
    if text:
        nodes.append(_full_text_node(text))
    if uniprot:
        nodes.extend(_uniprot_nodes(uniprot))
    if not nodes:
        raise ValueError("provide a free-text query and/or --uniprot")
    if len(nodes) == 1:
        return nodes[0]
    return {"type": "group", "logical_operator": "and", "nodes": nodes}


def search(text: str | None, uniprot: str | None, rows: int) -> list:
    payload = {
        "query": build_query(text, uniprot),
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_content_type": ["experimental"],
            "scoring_strategy": "combined",
        },
    }
    result = _post(SEARCH_URL, payload)
    return [hit["identifier"] for hit in result.get("result_set", [])]


def entry_details(pdb_id: str) -> dict:
    info = _get(ENTRY_URL + pdb_id)
    entry_info = info.get("rcsb_entry_info", {})
    return {
        "id": pdb_id,
        "title": info.get("struct", {}).get("title", ""),
        "resolution": entry_info.get("resolution_combined", []),
        "method": entry_info.get("experimental_method", ""),
        "ligands": entry_info.get("nonpolymer_bound_components", []) or [],
        "polymer_chains": entry_info.get("deposited_polymer_entity_instance_count"),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Search the RCSB PDB for structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("query", nargs="?", help="Free-text search terms")
    parser.add_argument("--uniprot", help="Restrict to a UniProt accession, e.g. P01116")
    parser.add_argument("--rows", type=int, default=25, help="Max hits (default: 25)")
    parser.add_argument(
        "--details",
        action="store_true",
        help="Fetch title/resolution/method/ligands for each hit (slower)",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args(argv)

    if not args.query and not args.uniprot:
        parser.error("provide a free-text query and/or --uniprot")

    try:
        ids = search(args.query, args.uniprot, args.rows)
    except (urllib.error.URLError, urllib.error.HTTPError) as exc:
        payload = {"status": "error", "error": f"RCSB search failed: {exc}"}
        print(json.dumps(payload, indent=2) if args.json else payload["error"])
        return 1

    if args.details:
        results = []
        for pdb_id in ids:
            try:
                results.append(entry_details(pdb_id))
            except (urllib.error.URLError, urllib.error.HTTPError):
                results.append({"id": pdb_id, "title": "(details unavailable)"})
    else:
        results = [{"id": pdb_id} for pdb_id in ids]

    data = {
        "query": args.query,
        "uniprot": args.uniprot,
        "count": len(results),
        "results": results,
    }

    if args.json:
        print(json.dumps({"status": "ok", "data": data}, indent=2))
    else:
        label = args.query or ""
        if args.uniprot:
            label = f"{label} [UniProt {args.uniprot}]".strip()
        print(f"{len(results)} hit(s) for {label!r}:")
        for entry in results:
            if args.details:
                res = entry.get("resolution") or ["?"]
                res_s = res[0] if res else "?"
                ligs = ",".join(entry.get("ligands", [])) or "-"
                print(f"  {entry['id']}  {res_s}A  {entry.get('method','')[:5]:5}  ligs={ligs}")
                print(f"        {entry.get('title','')[:90]}")
            else:
                print(f"  {entry['id']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
