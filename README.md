<p align="center">
  <img src="assets/banner.jpg" alt="Proteus — structural biology agent skill" width="100%">
</p>

# Proteus

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![CI](https://github.com/jvogan/proteus/actions/workflows/ci.yml/badge.svg)](https://github.com/jvogan/proteus/actions/workflows/ci.yml)

**Structural biology automation for AI coding agents.**

Proteus is a structural biology agent skill and stdlib-only helper toolkit for
[Codex](https://openai.com/index/codex/), [Claude Code](https://docs.anthropic.com/en/docs/claude-code),
and other AI coding agents. It teaches agents how to do AI protein structure
analysis from the terminal: resolve proteins, fetch structures, inspect
AlphaFold confidence, render PyMOL figures, automate ChimeraX analysis, compare
models, validate experimental structures, and reason about protein design
tooling.

Proteus is for computational biologists, protein engineers, educators, and
agent builders who want reliable PyMOL automation, ChimeraX automation,
AlphaFold DB / pLDDT / PAE workflows, RCSB PDB and UniProt helpers, and
Rosetta-oriented protein design guidance without building a custom plugin.

> **What is a skill?** A skill is a directory that an AI coding agent reads to gain domain-specific knowledge. Clone it into the agent's skills folder and it becomes part of the agent's working context — no code changes or plugin installs required.

## Why Proteus?

- **Turns vague protein prompts into reproducible workflows.** `TP53`, `P04637`,
  `4HHB`, and local `.pdb` / `.cif` files become local structures with provenance.
- **Keeps agents out of known tool traps.** Proteus documents hard-won PyMOL,
  ChimeraX, and AlphaFold DB gotchas that otherwise waste hours.
- **Works even before heavy tools are installed.** Zero-dependency scripts inspect
  PDB/mmCIF files, query public APIs, and emit structured JSON.
- **Scales up when PyMOL or ChimeraX are available.** Agents can render
  publication-quality PyMOL images or run ChimeraX SASA, H-bond, alignment, and
  cryo-EM workflows.
- **Produces outputs agents can chain.** Reports use machine-readable JSON for
  parallel runs, CI checks, notebooks, and downstream analysis.

## What It Provides

- **23 documented gotchas** for PyMOL, ChimeraX, AlphaFold DB, rendering, and cryo-EM maps — hard-won from real debugging
- Tool detection for PyMOL, ChimeraX, and ffmpeg across macOS and Linux installs
- Headless PyMOL rendering for publication-quality structure figures, with publication/illustration/soft presets and pLDDT coloring
- One-command annotated binding-pocket figures (ligand + pocket sticks + polar contacts + transparent context)
- Headless turntable movies (PyMOL ray-traced frames + ffmpeg), degrading gracefully when ffmpeg is absent
- Managed ChimeraX REST rendering — launches a GUI session, renders via GPU, defeats the 0-byte-PNG save race, and tears down cleanly
- ChimeraX analysis helpers for alignment, SASA, and hydrogen-bond workflows
- Zero-dependency protein-protein interface residue analysis between chains
- HELIX-record injection for CA-only backbones (RFdiffusion / Genie designs) so cartoons render correctly
- MRC/CCP4 map inspection with sigma-based contour-level suggestions
- Cryo-EM density-fit figures (model in carved density mesh; works with a real map or simulated density)
- AlphaFold DB fetch with confidence interpretation and pLDDT coloring
- RCSB PDB fetch for experimental coordinates, metadata, and biological assembly mmCIF
- RCSB PDB search by free text and/or UniProt accession
- Candidate structure ranking by method, resolution, validation hints, assemblies, and ligands
- Biological assembly availability reports and optional assembly mmCIF download
- PDBe SIFTS residue mapping between PDB chains and UniProt ranges
- Ligand inventory from local PDB/mmCIF files or PDB IDs, with optional CCD reference downloads
- Docking prep planning for receptor/ligand files without executing chemistry tools
- Ligand-centered PDB/mmCIF docking box generation and optional Vina/Open Babel/GNINA/RDKit detection
- Dry-run Vina-compatible docking command planning and log parsing
- Simple PDB/mmCIF protein-ligand interaction classes with optional PLIP/ProLIF detection
- Conservative protein substitution parsing with optional local coordinate lookup
- Local structural mutation triage against ligands, interfaces, close contacts, and PAE JSON
- Local JSON-manifest batch runner for allowlisted helper scripts
- Local URL cache and evidence-pack report combiner for reproducible outputs
- Generic target dossiers with Markdown reports, JSON provenance, and opt-in local analyses
- UniProt lookup for resolving gene/protein names before AlphaFold fetches
- PDB/mmCIF inspection via `structure_info.py`
- One-command readiness checks via `proteus_doctor.py`
- Query resolution via `resolve_structure.py` for local files, PDB IDs, UniProt accessions, and gene/protein names
- PAE, validation, ligand-pocket, and structure-comparison reports
- Optional model-quality tool detection plus USalign metric parsing/running when installed
- Dry-run design manifest planner plus Rosetta scorefile parser/scoring command planner
- KRAS G12C dossier workflow that chains fetches, analyses, and figures from public structures
- Rosetta/PyRosetta patterns plus ML alternatives (ProteinMPNN, ESM2)
- Zero-dependency PDB file inspector (`pdb_info.py` — stdlib only)
- Structured JSON output from analysis helper scripts, safe for parallel runs

## Agent Prompts That Work

Use prompts like these with `$proteus` or after installing this directory as an
agent skill:

```text
Use Proteus to resolve TP53, fetch the AlphaFold prediction, and summarize low-confidence regions.
Render the 1HSG binding pocket around indinavir in PyMOL and save a clean PNG.
Compare AF-P04637-F1 against an experimental p53 structure and report RMSD plus high-deviation residues.
Run a ChimeraX hydrogen-bond and SASA analysis for this protein-protein interface.
Check whether 4HHB has validation red flags before using it as a reference structure.
```

## Capabilities Matrix

| Capability | No local tools | PyMOL | ChimeraX | Public APIs | Rosetta/PyRosetta |
|---|---:|---:|---:|---:|---:|
| PDB/mmCIF inspection | yes | optional | optional | no | no |
| Protein/name resolution | yes | no | no | UniProt | no |
| Experimental structure fetch | yes | no | no | RCSB PDB | no |
| AlphaFold confidence, pLDDT, PAE | yes | render optional | optional | AlphaFold DB | no |
| Headless structure rendering | no | yes | limited | no | no |
| SASA, H-bonds, contacts, alignment | partial | partial | yes | no | optional |
| Protein design/scoring guidance | docs | optional | optional | optional | yes |

## Generated Outputs

Small, checked-in snapshots show the JSON shape without requiring downloads:

- [`resolve_structure.py TP53 --no-download --json`](docs/snapshots/resolve_tp53.json)
- [`pae_report.py tests/fixtures/tiny_pae.json --json`](docs/snapshots/pae_tiny.json)
- [`validation_report.py 1HSG --json`](docs/snapshots/validation_1hsg.json)
- [`pocket_report.py tests/fixtures/tiny.pdb --json`](docs/snapshots/pocket_tiny.json)

The repository also includes a curated social preview image at
[`assets/social-preview.jpg`](assets/social-preview.jpg). The larger generated
banner gallery stays ignored to keep the public repository lean.

## Try It In 60 Seconds

```bash
git clone https://github.com/jvogan/proteus.git
cd proteus
make test
python3 scripts/proteus_doctor.py --json
python3 scripts/resolve_structure.py TP53 --no-download --json
python3 scripts/pocket_report.py tests/fixtures/tiny.pdb --json
```

These commands require only Python 3.10+. Install PyMOL or ChimeraX later when
you want rendering, SASA, hydrogen bonds, structure alignment, or live molecular
visualization.

## Install

Clone into your agent's skills directory:

```bash
# Claude Code
git clone https://github.com/jvogan/proteus.git ~/.claude/skills/proteus

# Codex
git clone https://github.com/jvogan/proteus.git ~/.codex/skills/proteus
```

Or copy the directory manually into your agent's skills folder. The skill path may vary by agent version — check your agent's documentation if the above doesn't work.

## Runtime requirements

Proteus degrades gracefully — `pdb_info.py` and AlphaFold metadata fetches work with zero local tools. For full capability, install at least one:

| Tool | Role | Install |
|---|---|---|
| **PyMOL** | Headless rendering, structure inspection | [pymol.org](https://pymol.org) or `conda install -c conda-forge pymol-open-source` |
| **ChimeraX** | Analysis, GUI demos, cryo-EM visualization | [cgl.ucsf.edu/chimerax](https://www.cgl.ucsf.edu/chimerax/download.html) |
| **AlphaFold DB** | Public prediction database | No install — uses the [EBI REST API](https://alphafold.ebi.ac.uk) |
| **PyRosetta** | Scoring, energy minimization, protein design | `pip install pyrosetta-installer` (academic license required) |

Python 3.10+ is required. All helper scripts use only the standard library.

## Quick examples

With the skill installed, natural-language prompts trigger it automatically:

```
Fetch the AlphaFold prediction for p53 and show which regions look disordered.
Render the 1HSG binding pocket in PyMOL and save a clean PNG.
Compare an AlphaFold model to an experimental structure and report RMSD.
Analyze the hydrogen bonds at a protein-protein interface in ChimeraX.
```

The helper scripts also work standalone:

```bash
python3 scripts/pdb_info.py structure.pdb                          # zero-dep PDB inspection
python3 scripts/structure_info.py structure.cif --json             # PDB/mmCIF inspection
python3 scripts/fetch_pdb.py 4HHB --json                           # RCSB PDB fetch
python3 scripts/pdb_search.py "KRAS G12C" --rows 10 --details --json # RCSB search
python3 scripts/pdb_select.py --input candidates.json --json       # candidate ranking
python3 scripts/assembly_report.py 4HHB --json                     # biological assembly report
python3 scripts/sifts_map.py pdb 1hsg --json                       # PDB-to-UniProt mapping
python3 scripts/ligand_extract.py 1HSG --json                      # ligand inventory
python3 scripts/dock_prep.py --receptor receptor.pdb --ligand ligand.sdf --json # docking prep plan
python3 scripts/docking_box.py 1HSG --ligand MK1 --json            # ligand-centered docking box
python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt --json # Vina command plan
python3 scripts/interaction_report.py 1HSG --ligand MK1 --json     # ligand interactions
python3 scripts/variant_map.py "P04637 R175H" --no-download --json # substitution parser
python3 scripts/mutation_triage.py R175H --structure model.pdb --json # local variant triage
python3 scripts/proteus_batch.py manifest.json --dry-run --json    # batch manifest runner
python3 scripts/proteus_cache.py verify --json                     # local cache integrity
python3 scripts/proteus_report.py --input result.json --json       # evidence pack combiner
python3 scripts/design_run.py manifest.json --dry-run --json       # design run planner
python3 scripts/rosetta_score.py parse-scorefile score.sc --json   # Rosetta score parser
python3 scripts/target_dossier.py --uniprot P04637 --no-network --json # dossier report
python3 scripts/target_dossier.py --pdb tests/fixtures/tiny.cif --no-network --analyze-local --json # local analysis dossier
python3 scripts/uniprot_lookup.py TP53 --gene-exact --json         # UniProt lookup
python3 scripts/fetch_alphafold.py P04637 --pae --json             # AlphaFold fetch
python3 scripts/pae_report.py AF-P04637-F1_pae.json --json         # PAE/domain hints
python3 scripts/validation_report.py 4HHB --json                   # wwPDB validation metrics
python3 scripts/pocket_report.py tests/fixtures/tiny.cif --json     # ligand-pocket contacts
python3 scripts/interface_report.py 1BRS --chains A,D --json        # protein-protein interface residues
python3 scripts/model_quality.py detect --json                     # optional quality tool detection
python3 scripts/resolve_structure.py TP53 --json                   # one-command resolver
python3 scripts/pymol_agent.py render structure.pdb output.png     # headless render
python3 scripts/pymol_agent.py pocket 1HSG.pdb pocket.png --label  # annotated binding-pocket figure
python3 scripts/pymol_agent.py spin structure.pdb spin.mp4         # turntable movie (needs ffmpeg)
python3 scripts/chimerax_agent.py align reference.pdb mobile.pdb   # structure alignment
python3 scripts/chimerax_rest.py render structure.pdb out.png      # GPU render via managed REST
python3 scripts/add_helix_records.py model.pdb --json              # fix CA-only backbone cartoons
python3 scripts/map_info.py map.mrc --json                         # cryo-EM contour levels
python3 scripts/pymol_agent.py density model.pdb fit.png --simulate # model in (simulated) density
python3 scripts/kras_dossier.py --out kras_g12c_dossier --no-movie  # KRAS dossier workflow
```

See [`SHOWCASE.md`](SHOWCASE.md) for end-to-end, copy-pasteable examples on
canonical public structures (HIV protease, ACE2–RBD, hemoglobin, AlphaFold p53,
and more).

## Layout

```text
proteus/
├── SKILL.md              # Main skill — agent reads this first
├── agents/openai.yaml    # Codex discovery metadata
├── references/           # On-demand deep docs (loaded as needed)
│   ├── alphafold.md
│   ├── chimerax.md
│   ├── data-sources.md
│   ├── file-formats.md
│   ├── prediction-models.md
│   ├── pymol.md
│   └── rosetta.md
└── scripts/              # Agent helper scripts (all stdlib-only)
    ├── add_helix_records.py
    ├── chimerax_agent.py
    ├── chimerax_rest.py
    ├── compare_structures.py
    ├── design_run.py
    ├── dock_prep.py
    ├── dock_vina.py
    ├── docking_box.py
    ├── fetch_pdb.py
    ├── fetch_alphafold.py
    ├── interface_report.py
    ├── interaction_report.py
    ├── kras_dossier.py
    ├── assembly_report.py
    ├── ligand_extract.py
    ├── map_info.py
    ├── model_quality.py
    ├── mutation_triage.py
    ├── pae_report.py
    ├── pdb_info.py
    ├── pdb_select.py
    ├── pdb_search.py
    ├── pocket_report.py
    ├── proteus_batch.py
    ├── proteus_cache.py
    ├── proteus_doctor.py
    ├── proteus_report.py
    ├── resolve_structure.py
    ├── rosetta_score.py
    ├── sifts_map.py
    ├── structure_info.py
    ├── target_dossier.py
    ├── uniprot_lookup.py
    ├── variant_map.py
    ├── validation_report.py
    └── pymol_agent.py
```

## Design Principles

The tool split is deliberate:

- **PyMOL** is the default for headless image generation (software ray tracer — no display needed).
- **ChimeraX** is the default for analysis-heavy workflows and GPU-rendered GUI sessions.
- **RCSB/UniProt/AlphaFold DB** provide lightweight upstream data discovery before local visualization.
- **Rosetta/PyRosetta** are optional extensions. ML alternatives (ProteinMPNN, ESM2) are documented for when Rosetta isn't available.

Helper scripts emit machine-readable JSON, with human-readable text as a fallback. Temporary handoff files are per-process, so parallel agent runs never collide.

## Test

```bash
make test            # run the unit tests
make release-check   # pre-publish hygiene sweep (no tracked structures/maps/media/secrets, no large files)
```

## Safety And Privacy

Proteus does not include telemetry or credential collection. Some workflows call
public APIs with user-provided protein names, UniProt accessions, or PDB IDs.
Read [PRIVACY.md](PRIVACY.md), [SECURITY.md](SECURITY.md), and
[DISCLAIMER.md](DISCLAIMER.md) before using Proteus with private structures,
unpublished sequences, or regulated data.

## Contributing

Found a gotcha that isn't documented? Have a workflow that should be covered?
[Open an issue](https://github.com/jvogan/proteus/issues) or submit a PR. The
most valuable contributions are real debugging discoveries — the kind of thing
that takes hours to figure out and one sentence to explain. See
[CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT. See [LICENSE](LICENSE).
