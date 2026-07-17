---
name: proteus
description: >
  Use this skill when the user asks you to work with protein structures,
  molecular visualization, or structural biology tools. Trigger for PyMOL,
  ChimeraX, AlphaFold, Rosetta/PyRosetta, UniProt, RCSB PDB, PDBe, PDB/mmCIF
  files, molecular rendering, pLDDT, PAE, RMSD, structure alignment, binding
  pockets, interfaces, variants, docking context, cryo-EM maps, ensembles,
  restraints, or protein design. Also trigger when the user opens .pdb, .cif,
  .mmcif, .sdf, .mol2, .mrc, .map, or .ccp4 files. Do not trigger for general
  biology with no structural component, sequence-only bioinformatics, genomics,
  or transcriptomics.
---

# Proteus — Structural Biology Agent Skill

Use Proteus to turn structural-biology questions into inspectable analyses,
figures, sessions, and replayable PyMOL (`.pml`) or ChimeraX (`.cxc`) commands.
Prefer the included helper workflows over improvised command strings.

## Safety and Scientific Ground Rules

- Treat local structures, maps, annotations, sequences, and restraints as private
  unless the user identifies them as public. Do not upload private or unpublished
  inputs to a web service without explicit authorization.
- Public identifiers such as PDB, EMDB, UniProt, and AlphaFold DB accessions may
  be resolved through their public APIs when network access is needed. Say which
  source was used and preserve provenance.
- Keep secrets, credentials, arbitrary absolute paths, raw private data, and
  machine-specific details out of reports, logs, manifests, and repository files.
- Never use a protein-structure result as the sole basis for a clinical, safety,
  efficacy, or patient-specific decision. Flag medical or biosafety implications
  and recommend appropriate expert review.
- Separate observation from inference. A short distance is a geometric contact,
  not proof of a hydrogen bond, energetic hotspot, catalytic role, binding
  affinity, biological assembly, or causal variant effect.
- State uncertainty from resolution, missing atoms/residues, alternate
  conformers, occupancy, protonation, tautomerism, model confidence, chain
  mapping, alignment coverage, map processing, and assembly choice when relevant.
- Do not overwrite source coordinates. Write generated artifacts to a dedicated
  output directory and keep scripts/sessions so the result can be audited.

## Start Here

Run scripts from the skill or repository root with Python 3.10+.

```bash
python3 scripts/proteus.py --help
python3 scripts/proteus.py qc structure.cif --json
```

Do not run network checks merely because the skill was invoked. Use
`scripts/proteus_doctor.py --json` for local readiness, and add `--network` only
when the task actually needs public APIs.

Before using a helper, run its `--help`. Treat helpers as black boxes unless
debugging or changing the implementation.

## Workflow Routing

| User goal | Preferred workflow |
|---|---|
| Preflight models, conformers, occupancy, components | `proteus.py qc` |
| Reproducible publication figure or saved scene | `proteus.py figure` |
| Apo/holo, WT/mutant, or pose comparison | `proteus.py compare` |
| Residue/variant local environment | `proteus.py residue` |
| Protein-protein interface story | `proteus.py interface` |
| Per-residue annotations or score overlays | `proteus.py annotate` |
| Distance restraints/crosslinks | `proteus.py restraints` |
| Assembly and crystal-neighbor comparison | `proteus.py assembly` |
| Cryo-EM contour sweep, fit, map/model review | `proteus.py cryoem` |
| Multi-model ensemble variability | `proteus.py ensemble` |
| Qualitative electrostatic surface | `proteus.py electrostatics` |
| Optional fpocket/P2Rank detection | `proteus.py pockets` |
| Ligand, cofactor, metal, or water site | `proteus.py chemical-site` |
| Fetch/resolve public coordinates | `resolve_structure.py`, `fetch_pdb.py`, `fetch_alphafold.py` |
| PAE or model-quality metrics | `pae_report.py`, `model_quality.py` |
| Docking preparation and interaction triage | `dock_prep.py`, `docking_box.py`, `dock_vina.py`, `interaction_report.py` |
| Design planning or Rosetta score parsing | `design_run.py`, `rosetta_score.py` |

Read `references/workflows.md` for inputs, examples, output contracts, scene
manifests, and human-handoff patterns.

## Choosing PyMOL or ChimeraX

| Task | Default | Reason |
|---|---|---|
| Headless ray-traced rendering | PyMOL | Software ray tracer works without a display |
| GPU volume/surface rendering | ChimeraX managed REST | Strong map and surface tooling |
| H-bonds, contacts, clashes, SASA, buried area | ChimeraX | Built-in analysis commands |
| Scriptable publication scenes | Either | Proteus emits replayable scripts and sessions |
| Live interactive handoff | ChimeraX REST or saved session | User can continue in the GUI |
| Structural alignment | Either | Choose by metric and downstream rendering needs |

On macOS, ChimeraX `--nogui` can analyze but cannot render because it lacks an
OpenGL context. Use `scripts/chimerax_agent.py` for headless analysis and
`scripts/chimerax_rest.py` or a generated workflow with `--execute` for GUI/GPU
rendering. PyMOL headless rendering works on macOS and Linux.

## Standard Operating Pattern

1. Resolve the input identity and numbering system. Do not assume UniProt,
   author, label, and construct residue numbers are interchangeable.
2. Run `proteus.py qc` before geometry-sensitive analysis. Choose a model and
   alternate-conformer policy explicitly when multiple states are present.
3. Confirm whether the task needs the asymmetric unit, a biological assembly,
   crystal neighbors, or a predicted monomer/complex.
4. Choose the narrowest workflow that answers the question. Keep network use
   and optional external tools explicit.
5. Generate scripts and a JSON report first. Execute rendering or optional tools
   only when needed; retain `.pml`/`.cxc` and `.pse`/`.cxs` handoff artifacts.
6. Verify expected outputs are non-empty. Report measurements, provenance,
   selection/model decisions, limitations, and missing capabilities.

## Core Helpers

The unified entry point covers the main user workflows:

```bash
python3 scripts/proteus.py figure scene.json --outdir figure --execute
python3 scripts/proteus.py compare apo.pdb holo.pdb --outdir comparison
python3 scripts/proteus.py residue model.cif A:42 --outdir residue_42
python3 scripts/proteus.py interface complex.cif --chains A,B --execute
python3 scripts/proteus.py annotate model.pdb scores.csv --execute
python3 scripts/proteus.py restraints model.pdb restraints.csv --execute
python3 scripts/proteus.py assembly 4HHB --assembly 1 --execute
python3 scripts/proteus.py cryoem model.cif map.mrc --resolution 3.2 --fit --execute
python3 scripts/proteus.py ensemble models.pdb --execute
python3 scripts/proteus.py electrostatics model.pdb --execute
python3 scripts/proteus.py pockets model.pdb --detector auto --execute --render
python3 scripts/proteus.py chemical-site complex.cif --component ZN:A:501 --execute
```

Additional focused helpers include:

- `resolve_structure.py`: local file, PDB ID, UniProt accession, or name to a
  local structure with provenance.
- `pdb_search.py` and `pdb_select.py`: search and rank experimental candidates.
- `assembly_report.py`, `sifts_map.py`, and `validation_report.py`: assembly,
  residue mapping, and experimental validation context.
- `ligand_extract.py`, `pocket_report.py`, and `interaction_report.py`: ligand
  inventory and zero-dependency contact triage.
- `mutation_triage.py`: variant proximity to ligands, interfaces, contacts, and
  PAE context.
- `model_quality.py`: USalign, DockQ, and local Foldseek wrappers when installed.
- `proteus_batch.py`, `proteus_cache.py`, and `proteus_report.py`: allowlisted
  batches, cached public responses, and evidence packs.

## Interpretation Rules by Data Type

- **AlphaFold pLDDT:** confidence in local structure, not stability or function.
  Use per-residue C-alpha values when available. Low pLDDT may indicate disorder,
  context dependence, or prediction uncertainty.
- **PAE:** directional expected alignment error. Use chain-aware blocks for
  complexes and do not call high inter-chain PAE proof that chains do not bind.
- **RMSD/TM-score:** always report alignment method, aligned length/coverage,
  chain mapping, and normalization. A single global RMSD can hide local changes.
- **DockQ:** review model/native chain mapping and each interface, not only the
  aggregate score.
- **Contacts/interfaces:** distance cutoffs generate candidates. Protonation,
  donor/acceptor geometry, solvation, conformational state, and assembly choice
  determine meaning.
- **Cryo-EM:** contour levels are map-specific absolute values. Fitting and
  difference maps depend on initial placement, resolution, sharpening, masking,
  and scale; avoid overfitting.
- **Electrostatics:** ChimeraX Coulombic coloring is qualitative. For stronger
  claims, prepare protonation/charges and use a validated Poisson-Boltzmann
  workflow such as PDB2PQR/APBS.
- **Pockets/tunnels:** detector scores are hypotheses, not affinity estimates.
  Compare conformational states and add experimental/evolutionary evidence.

## Critical Execution Gotchas

- Write complex PyMOL commands to `.pml`; do not pass selections containing
  comparison operators through `pymol -d`.
- End headless PyMOL scripts with `quit`; use `iterate_state` for coordinates.
- PyMOL `cealign` argument order is target/reference first, mobile second.
- Use layered `b > 50`, `b > 70`, `b > 90` overrides for pLDDT; PyMOL has no
  `<=` selection operator.
- ChimeraX model IDs do not reset after `close session`. Generated workflows use
  a fresh managed session; otherwise discover model IDs dynamically.
- For ChimeraX REST image saves, issue `wait 1` and verify the file is non-empty.
- Avoid whole-map PyMOL isosurfaces; carve around the model or use ChimeraX.
- Deposited coordinates are usually an asymmetric unit, not automatically the
  biological assembly.

## Reading Guide

Load only the relevant reference:

| Topic | Reference |
|---|---|
| User workflows, manifests, outputs, handoff | `references/workflows.md` |
| PyMOL commands, rendering, animation | `references/pymol.md` |
| ChimeraX analysis, REST rendering, cryo-EM | `references/chimerax.md` |
| AlphaFold DB, pLDDT, PAE | `references/alphafold.md` |
| PDB/UniProt/PDBe/RCSB APIs | `references/data-sources.md` |
| PDB/mmCIF/SDF/MRC format choice | `references/file-formats.md` |
| AF3/Boltz/Chai/ColabFold and other models | `references/prediction-models.md` |
| Rosetta, PyRosetta, and ML design alternatives | `references/rosetta.md` |

## Output Contract

Helpers emit JSON with `status`, `data`, and when applicable `warnings` and
`provenance`. New workflow envelopes also include `schema_version` and
`proteus_version`. Read `data` first.

Return to the user:

- the answer and key measurements;
- generated figure/session/script paths;
- input source and model/chain/selection decisions;
- relevant uncertainty and interpretation limits;
- missing tools or incomplete outputs, without pretending success.
