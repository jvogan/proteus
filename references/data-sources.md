# Structural Data Sources - Proteus Skill

Use this reference when a user gives a protein name, gene symbol, PDB ID,
ligand ID, or asks for "the experimental structure" rather than a local file.

## Decision Map

| User gives... | First action | Tool/script |
|---|---|---|
| PDB ID, e.g. `4HHB` | Download experimental coordinates and metadata | `scripts/fetch_pdb.py` |
| Protein/gene name, e.g. `p53` | Resolve to UniProt accession first | `scripts/uniprot_lookup.py` |
| UniProt accession, e.g. `P04637` | Fetch AlphaFold prediction | `scripts/fetch_alphafold.py` |
| Local `.pdb` / `.cif` | Inspect locally before visualization | `scripts/structure_info.py` |
| Ligand or bound component ID | Use PDBe/RCSB metadata, then inspect pocket locally | `fetch_pdb.py` + ChimeraX/PyMOL |
| Cryo-EM map | Use ChimeraX workflows | `references/chimerax.md` |

## RCSB PDB

Best for experimental structure files and entry metadata.

Common downloads:

```text
https://files.rcsb.org/download/4HHB.cif
https://files.rcsb.org/download/4HHB.pdb
https://files.rcsb.org/download/4HHB-assembly1.cif
https://models.rcsb.org/4hhb.bcif
```

Metadata:

```text
https://data.rcsb.org/rest/v1/core/entry/4HHB
```

Use `scripts/fetch_pdb.py` rather than constructing URLs manually. It returns
selected metadata plus a local coordinate path.

## UniProt

Best for resolving natural-language proteins and gene symbols to accessions.

Common agent path:

```bash
python scripts/uniprot_lookup.py TP53 --gene-exact --json
python scripts/fetch_alphafold.py P04637 --pae --json
```

Default lookup filters to reviewed human UniProtKB entries (`organism_id:9606`).
Use `--all-organisms` or `--organism TAXON_ID` when the user specifies another
species.

## AlphaFold DB

Best for single-chain predicted structures with pLDDT and optional PAE.

Use `scripts/fetch_alphafold.py`; it queries AlphaFold DB metadata first and
uses the returned URLs instead of hardcoding model versions.

## PDBe

Best for richer entry-level annotations, ligand/cross-reference metadata,
residue mappings, and quality information. Prefer PDBe when the task asks for:

- residue-level mappings to UniProt
- ligand-centric metadata
- validation/quality summaries
- aggregated views across many PDB entries

Keep PDBe calls targeted. Do not bulk-download data unless the user asks.
