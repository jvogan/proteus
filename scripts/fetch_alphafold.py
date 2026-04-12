#!/usr/bin/env python3
"""Fetch AlphaFold predicted structures from the EBI database — Proteus skill.

No external dependencies — uses only Python stdlib.

Usage:
    python fetch_alphafold.py P69905              # hemoglobin alpha
    python fetch_alphafold.py P04637 --pae        # p53 + PAE matrix
    python fetch_alphafold.py P69905 --cif        # mmCIF format
    python fetch_alphafold.py P69905 --outdir ./data
    python fetch_alphafold.py --help
"""

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"


def fetch_metadata(uniprot_id: str) -> dict:
    """Fetch prediction metadata from AlphaFold DB API.

    Returns the first entry from the API response.
    NOTE: The API returns a list [{...}], not a plain dict.
    """
    url = f"{API_BASE}/{uniprot_id}"
    try:
        with urllib.request.urlopen(url) as resp:
            data = json.loads(resp.read())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"ERROR: UniProt ID '{uniprot_id}' not found in AlphaFold DB.", file=sys.stderr)
            print("  Common issue: P62988 (ubiquitin) is not in the DB. Use P0CG48 instead.", file=sys.stderr)
            sys.exit(1)
        raise
    # API returns a list — unwrap
    if isinstance(data, list):
        if not data:
            raise RuntimeError(f"No AlphaFold entries returned for '{uniprot_id}'.")
        return data[0]
    return data


def log(message: str, *, as_json: bool = False):
    """Send status output to stderr for machine-safe JSON mode."""
    stream = sys.stderr if as_json else sys.stdout
    print(message, file=stream)


def download(url: str, dest: Path, *, as_json: bool = False):
    """Download a file from URL to local path."""
    log(f"  Downloading {dest.name}...", as_json=as_json)
    urllib.request.urlretrieve(url, dest)
    log(f"  -> {dest.resolve()} ({dest.stat().st_size:,} bytes)", as_json=as_json)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch AlphaFold predicted structures from EBI database.",
        epilog="Examples:\n"
               "  %(prog)s P69905              # hemoglobin alpha\n"
               "  %(prog)s P04637 --pae        # p53 with PAE matrix\n"
               "  %(prog)s P01308 --cif        # insulin in mmCIF format\n"
               "\nKnown issues:\n"
               "  - P62988 (ubiquitin) is NOT in AlphaFold DB; use P0CG48\n"
               "  - Always queries latestVersion to avoid 404s from version drift",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("uniprot_id", help="UniProt accession (e.g., P69905, P04637)")
    parser.add_argument("--pae", action="store_true", help="Also download PAE (Predicted Aligned Error) JSON")
    parser.add_argument("--cif", action="store_true", help="Download mmCIF instead of PDB format")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current)")
    parser.add_argument("--json", action="store_true", help="Output metadata as JSON to stdout")
    args = parser.parse_args()

    uid = args.uniprot_id.upper()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Fetch metadata (always do this first to get correct version URL)
    log(f"Fetching AlphaFold metadata for {uid}...", as_json=args.json)
    meta = fetch_metadata(uid)

    gene = meta.get("gene", "N/A")
    version = meta.get("latestVersion", "N/A")
    plddt = meta.get("globalMetricValue", "N/A")
    seq_start = meta.get("sequenceStart", "?")
    seq_end = meta.get("sequenceEnd", "?")
    frac_vh = meta.get("fractionPlddtVeryHigh", 0)
    frac_c = meta.get("fractionPlddtConfident", 0)
    frac_l = meta.get("fractionPlddtLow", 0)
    frac_vl = meta.get("fractionPlddtVeryLow", 0)

    log(f"  Gene: {gene}", as_json=args.json)
    log(f"  Version: v{version}", as_json=args.json)
    log(f"  Residues: {seq_start}-{seq_end}", as_json=args.json)
    log(f"  Global pLDDT: {plddt}", as_json=args.json)
    log(
        f"  Confidence: >90={frac_vh:.1%}  70-90={frac_c:.1%}  50-70={frac_l:.1%}  <50={frac_vl:.1%}",
        as_json=args.json,
    )

    # Download structure
    model_id = meta.get("modelEntityId", f"AF-{uid}-F1")
    if args.cif:
        struct_url = meta["cifUrl"]
        struct_path = outdir / f"{model_id}.cif"
    else:
        struct_url = meta["pdbUrl"]
        struct_path = outdir / f"{model_id}.pdb"
    download(struct_url, struct_path, as_json=args.json)

    # Download PAE if requested
    pae_path = None
    if args.pae:
        pae_url = meta.get("paeDocUrl")
        if pae_url:
            pae_path = outdir / f"{model_id}_pae.json"
            download(pae_url, pae_path, as_json=args.json)
        else:
            log("  WARNING: No PAE data available for this entry.", as_json=args.json)

    # JSON output
    if args.json:
        output = {
            "uniprot_id": uid,
            "model_id": model_id,
            "gene": gene,
            "version": version,
            "global_plddt": plddt,
            "structure_path": str(struct_path.resolve()),
            "pae_path": str(pae_path.resolve()) if pae_path else None,
            "confidence": {
                "very_high": frac_vh,
                "confident": frac_c,
                "low": frac_l,
                "very_low": frac_vl,
            },
        }
        print(json.dumps(output, indent=2))
    else:
        print("Done.")


if __name__ == "__main__":
    main()
