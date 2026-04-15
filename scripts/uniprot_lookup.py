#!/usr/bin/env python3
"""Resolve protein names or gene symbols to UniProt accessions.

Uses only Python stdlib and the UniProt REST API. This is the bridge from
natural-language protein names ("human p53", "EGFR") to AlphaFold DB accessions.

Usage:
    python uniprot_lookup.py TP53
    python uniprot_lookup.py "hemoglobin alpha" --size 3
    python uniprot_lookup.py ubiquitin --all-organisms --json
"""

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_FIELDS = "accession,id,protein_name,gene_names,organism_name,reviewed,length"


def _looks_like_accession(query: str) -> bool:
    return bool(re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?|[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?", query.upper()))


def build_query(term: str, organism: str | None, reviewed_only: bool, gene_exact: bool) -> str:
    term = term.strip()
    if _looks_like_accession(term):
        query = f"accession:{term.upper()}"
    elif gene_exact:
        query = f"gene_exact:{term}"
    else:
        query = term

    parts = [f"({query})"]
    if organism:
        parts.append(f"organism_id:{organism}")
    if reviewed_only:
        parts.append("reviewed:true")
    return " AND ".join(parts)


def _protein_name(entry: dict) -> str | None:
    desc = entry.get("proteinDescription", {})
    recommended = desc.get("recommendedName", {})
    full = recommended.get("fullName", {})
    if full.get("value"):
        return full["value"]
    alternatives = desc.get("alternativeNames") or []
    for alt in alternatives:
        value = alt.get("fullName", {}).get("value")
        if value:
            return value
    return None


def _gene_names(entry: dict) -> list[str]:
    names = []
    for gene in entry.get("genes", []) or []:
        primary = gene.get("geneName", {}).get("value")
        if primary:
            names.append(primary)
        for synonym in gene.get("synonyms", []) or []:
            value = synonym.get("value")
            if value:
                names.append(value)
    return names


def normalize_entry(entry: dict) -> dict:
    return {
        "accession": entry.get("primaryAccession"),
        "id": entry.get("uniProtkbId"),
        "reviewed": entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
        "protein_name": _protein_name(entry),
        "gene_names": _gene_names(entry),
        "organism": entry.get("organism", {}).get("scientificName"),
        "taxon_id": entry.get("organism", {}).get("taxonId"),
        "length": entry.get("sequence", {}).get("length"),
    }


def search_uniprot(query: str, size: int) -> list[dict]:
    params = urllib.parse.urlencode({
        "query": query,
        "format": "json",
        "size": str(size),
        "fields": DEFAULT_FIELDS,
    })
    url = f"{UNIPROT_SEARCH}?{params}"
    with urllib.request.urlopen(url, timeout=30) as response:
        data = json.load(response)
    return [normalize_entry(entry) for entry in data.get("results", [])]


def main():
    parser = argparse.ArgumentParser(
        description="Resolve protein names or gene symbols to UniProt accessions.",
        epilog=(
            "Examples:\n"
            "  %(prog)s TP53\n"
            "  %(prog)s 'hemoglobin alpha' --size 3\n"
            "  %(prog)s ubiquitin --all-organisms --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("query", help="Gene, protein name, or UniProt accession")
    parser.add_argument("--organism", default="9606",
                        help="NCBI taxonomy ID to prefer (default: 9606, human)")
    parser.add_argument("--all-organisms", action="store_true", help="Do not filter by organism")
    parser.add_argument("--include-unreviewed", action="store_true", help="Include TrEMBL/unreviewed entries")
    parser.add_argument("--gene-exact", action="store_true", help="Treat query as an exact gene symbol")
    parser.add_argument("--size", type=int, default=5, help="Number of candidates to return (default: 5)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON to stdout")
    args = parser.parse_args()

    organism = None if args.all_organisms else args.organism
    query = build_query(args.query, organism, not args.include_unreviewed, args.gene_exact)
    results = search_uniprot(query, args.size)

    output = {
        "status": "ok",
        "query": args.query,
        "uniprot_query": query,
        "results": results,
    }

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        if not results:
            print("No UniProt matches found.", file=sys.stderr)
            sys.exit(1)
        for idx, result in enumerate(results, start=1):
            genes = ", ".join(result["gene_names"][:4]) if result["gene_names"] else "(no genes)"
            reviewed = "reviewed" if result["reviewed"] else "unreviewed"
            print(f"{idx}. {result['accession']} {result.get('id') or ''} [{reviewed}]")
            print(f"   Protein: {result.get('protein_name') or '(unknown)'}")
            print(f"   Genes: {genes}")
            print(f"   Organism: {result.get('organism') or '(unknown)'}; length={result.get('length')}")


if __name__ == "__main__":
    main()
