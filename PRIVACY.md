# Privacy

Proteus does not include telemetry, analytics, background services, accounts, or
credential collection.

## What Leaves Your Machine

Some helper scripts contact public scientific APIs when you ask them to resolve
or fetch data:

- `uniprot_lookup.py` sends gene/protein queries and organism filters to UniProt.
- `fetch_pdb.py` and `validation_report.py` send PDB IDs to RCSB PDB.
- `fetch_alphafold.py` and `resolve_structure.py` may send UniProt accessions or
  protein-name-derived accessions to AlphaFold DB.
- ChimeraX and PyMOL may access network resources if you run commands such as
  `fetch`, `open ... from pdb`, or `open ... from alphafold`.

Local coordinate files remain local unless you explicitly upload, commit, share,
or pass them to a tool/API that performs network access.

## Generated Files

Proteus may write downloaded structures, JSON reports, renders, or temporary
handoff files in your working directory or the output directory you choose.
Generated AlphaFold files are ignored by default through `.gitignore`.

## No Secrets Required

The default Proteus workflows do not require API keys. Do not put private tokens,
clinical data, proprietary coordinates, or unpublished sequences in repository
files unless you intentionally want to publish them.
