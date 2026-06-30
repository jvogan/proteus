# Structural Data Sources - Proteus Skill

Use this reference when a user gives a protein name, gene symbol, PDB ID,
ligand ID, or asks for "the experimental structure" rather than a local file.

## Decision Map

| User gives... | First action | Tool/script |
|---|---|---|
| Ambiguous query | Resolve to the best structure source | `scripts/resolve_structure.py` |
| Mixed target dossier | Build Markdown plus JSON provenance; add `--analyze-local` for local ligand/interface summaries | `scripts/target_dossier.py` |
| Structure landscape or target survey | Search many experimental structures | `scripts/pdb_search.py` |
| Structure selection from candidates | Rank by method, resolution, validation, assembly, and ligand fit | `scripts/pdb_select.py` |
| PDB ID, e.g. `4HHB` | Download experimental coordinates and metadata | `scripts/fetch_pdb.py` |
| Biological assembly question | Check assembly files and recommended download | `scripts/assembly_report.py` |
| Protein/gene name, e.g. `p53` | Resolve to UniProt accession first | `scripts/uniprot_lookup.py` |
| UniProt accession, e.g. `P04637` | Fetch AlphaFold prediction | `scripts/fetch_alphafold.py` |
| UniProt/PDB residue numbering | Use PDBe SIFTS mappings | `scripts/sifts_map.py` |
| Local `.pdb` / `.cif` | Inspect locally before visualization | `scripts/structure_info.py` |
| Ligand or bound component ID | Inventory HETATM ligands and fetch CCD reference files if needed | `scripts/ligand_extract.py` |
| Receptor/ligand prep question | Plan prep commands and risk checks without executing tools | `scripts/dock_prep.py` |
| Docking box from known ligand | Compute ligand-centered center/size and detect optional docking tools | `scripts/docking_box.py` |
| Vina docking handoff | Plan Vina-compatible commands or parse Vina logs without executing docking | `scripts/dock_vina.py` |
| Protein-ligand interactions | Summarize simple contacts and optional PLIP/ProLIF availability | `scripts/interaction_report.py` |
| Reused public API/download response | Cache, verify, or reuse by URL checksum | `scripts/proteus_cache.py` |
| Collected helper JSON files | Combine into Markdown and scrubbed JSON evidence pack | `scripts/proteus_report.py` |
| Design/modeling run manifest | Validate inputs, tools, stages, and dry-run commands | `scripts/design_run.py` |
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

Use `scripts/assembly_report.py` before making oligomer or interface claims from
deposited coordinates:

```bash
python3 scripts/assembly_report.py 4HHB --json
python3 scripts/assembly_report.py 4HHB --download --outdir structures --json
```

Use `scripts/pdb_search.py` when selecting among many candidate structures:

```bash
python3 scripts/pdb_search.py "KRAS G12C sotorasib" --rows 10 --details --json
python3 scripts/pdb_search.py --uniprot P01116 --rows 50 --details --json
```

Use `scripts/pdb_select.py` when candidate metadata is already available or
when a live RCSB lookup should be ranked by practical modeling criteria:

```bash
python3 scripts/pdb_select.py --input candidates.json --ligand ATP --json
python3 scripts/pdb_select.py --live --query "KRAS G12C sotorasib" --rows 20 --json
```

## UniProt

Best for resolving natural-language proteins and gene symbols to accessions.

Use `scripts/sifts_map.py` when residue numbering matters:

```bash
python3 scripts/sifts_map.py pdb 1hsg --json
python3 scripts/sifts_map.py uniprot P04637 --limit 20 --json
```

Common agent path:

```bash
python3 scripts/uniprot_lookup.py TP53 --gene-exact --json
python3 scripts/fetch_alphafold.py P04637 --pae --json
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

## Validation

RCSB exposes wwPDB validation summaries through the Data API. Use
`scripts/validation_report.py` when selecting between experimental structures or
when a user asks whether a structure is trustworthy. It reports geometry summary
fields such as clashscore, bond/angle RMSZ, Ramachandran outliers, and rotamer
outliers when available.
