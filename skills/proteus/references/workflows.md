# Reproducible PyMOL and ChimeraX Workflows

This guide describes Proteus workflows that produce a JSON report plus
replayable PyMOL (`.pml`) or ChimeraX (`.cxc`) scripts. Most also write a saved
session when executed. Run `python3 scripts/proteus.py <command> --help` for the
complete current interface.

## General output pattern

Each workflow writes into an explicit output directory:

- a JSON report with inputs, parameters, warnings, provenance, and artifacts;
- a command file that can be inspected before execution;
- PNG figures and `.pse`/`.cxs` sessions when rendering is requested;
- tool results when an optional detector or quality program is run.

Input coordinate and map files are never modified. Reports use portable path
labels and file checksums rather than publishing arbitrary absolute paths.

## Structure preflight

```bash
python3 scripts/proteus.py qc structure.cif --model first --altloc highest --json
```

Run this before distances, interfaces, comparisons, restraints, or chemical-site
analysis. It reports available/selected models, alternate sites, occupancy,
missing backbone atoms, residue-number gaps, and chemical-component roles.

Use `--model all` only for an ensemble-aware workflow. Use an explicit model
number when the deposited file contains different experimental or modeled
states. `--altloc highest` chooses the highest-occupancy conformer per atom and
records that policy.

## Declarative figures and guided views

```bash
python3 scripts/proteus.py figure scene.json --outdir figure --execute
```

Example `scene.json`:

```json
{
  "name": "active-site-story",
  "tool": "pymol",
  "width": 1400,
  "height": 1000,
  "background": "white",
  "structures": [{"id": "complex", "path": "complex.cif"}],
  "representations": [
    {"style": "cartoon", "selection": "polymer", "color": "chain"},
    {"style": "sticks", "selection": "organic or byres (polymer within 4.5 of organic)", "color": "element"}
  ],
  "labels": [{"selection": "organic and not hydro", "text": "bound ligand", "color": "black"}],
  "views": [
    {"id": "overview", "selection": "polymer", "output": "overview.png"},
    {"id": "site", "selection": "organic", "turn": {"y": 25}, "output": "site.png"}
  ]
}
```

Choose `"tool": "chimerax"` for a ChimeraX script/session. Multiple ordered
views act as deterministic scene checkpoints for a presentation or human
handoff. Keep manifest selections tool-specific and use filenames—not paths—for
view outputs.

## State comparisons

```bash
python3 scripts/proteus.py compare apo.pdb holo.pdb \
  --ligand ATP --outdir apo_holo --execute
```

This produces aligned overlay and side-by-side views. When a ligand is selected,
the report compares nearby contact residues between states. Interpret gained or
lost contacts only after checking chain mapping, ligand identity, alignment
coverage, missing residues, and alternate conformers.

## Residue and interface stories

```bash
python3 scripts/proteus.py residue model.cif A:42 --outdir residue_42 --execute
python3 scripts/proteus.py interface complex.cif --chains A,B --cutoff 5 --outdir interface_ab --execute
```

The residue workflow combines a global locator, labeled local environment, and
available ligand/interface/variant context. The interface workflow combines
local contact residues with ChimeraX buried area, H-bond, contact, and clash
commands plus PyMOL figures.

Buried area and distance contacts depend on which assembly is loaded. Do not
interpret crystal-packing contacts as a biological interface without supporting
evidence.

## Annotation overlays

```bash
python3 scripts/proteus.py annotate model.pdb annotations.csv --outdir annotated --execute
```

CSV rows can supply chain, residue, score, color, and label fields. JSON lists
use the same keys. Use this for variants, domains, PTMs, conservation,
confidence, validation, or a user-defined score. Keep the score definition and
scale in the surrounding report; colors alone are not self-explanatory.

## Experimental restraints

```bash
python3 scripts/proteus.py restraints model.pdb restraints.csv --outdir restraints --execute
```

Required CSV/JSON fields are `residue1` and `residue2`. Optional fields include
`chain1`, `atom1`, `chain2`, `atom2`, `min`, `max`, `id`, and `label`. The report
marks pairs as satisfied, too short, too long, or unresolved and colors the
rendered distance objects.

Distance bounds must match the experiment and labeling chemistry. A C-alpha
distance is not automatically comparable to an atom-to-atom crosslink, FRET, or
DEER distribution.

## Assembly and symmetry

```bash
python3 scripts/proteus.py assembly 4HHB --assembly 1 --outdir assembly --execute
```

PDB identifiers enable assembly metadata lookup. The ChimeraX workflow compares
the asymmetric unit and requested biological assembly; PyMOL expands nearby
crystal mates. For a local coordinate file, assembly availability cannot be
inferred reliably unless the necessary metadata/operators are present.

## Cryo-EM map/model review

```bash
python3 scripts/proteus.py cryoem model.cif map.mrc \
  --resolution 3.2 --fit --difference --outdir map_review --execute
```

Local MRC/CCP4 maps receive sigma-based contour suggestions. Remote EMDB maps
require explicit absolute `--levels` because Proteus cannot compute local map
statistics before download/open. The ChimeraX script saves a contour sweep,
optional rigid fit, map statistics, and a qualitative difference view.

Use several contour levels. Record resolution, map provenance, sharpening or
filtering, masking, and fit procedure. Difference density is not reliable
without appropriate scaling and independent validation.

## Ensembles

```bash
python3 scripts/proteus.py ensemble models.pdb --outdir ensemble --execute
```

Multi-model PDB/mmCIF files are compared using common C-alpha atoms. With NumPy,
Proteus performs a Kabsch superposition before RMSF; otherwise it reports the
translation-only fallback. The result includes residue-level RMSF, most-variable
residues, an all-state PyMOL view, and a ChimeraX coordinate-set session.

RMSF here measures dispersion across supplied coordinate models, not uncertainty
or dynamics unless the ensemble itself supports that interpretation.

## Electrostatics

```bash
python3 scripts/proteus.py electrostatics model.pdb \
  --selection protein --range 10 --offset 1.4 --outdir electrostatics --execute
```

This uses ChimeraX Coulombic surface coloring and reports whether PDB2PQR/APBS
are installed. It is useful for qualitative orientation and charge patches.
For quantitative comparisons, curate protonation and charges and use a validated
Poisson-Boltzmann workflow with controlled dielectric and ionic-strength choices.

## Pocket and tunnel candidates

```bash
python3 scripts/proteus.py pockets model.pdb --detector auto --execute --render
```

`auto` prefers fpocket, then P2Rank. Detection is optional and local; no
structure is uploaded. Detected pocket files are copied into the output directory
and rendered with lining residues. CAVER/HOLE availability is reported, but
tunnel execution is left to a configured tool-specific workflow because starting
points and probe settings are scientifically consequential.

## Chemical sites

```bash
python3 scripts/proteus.py chemical-site complex.cif \
  --component ZN:A:501 --radius 5 --outdir zinc_site --execute
```

Component selectors may be `RESN`, `CHAIN:RESIDUE`, or
`RESN:CHAIN:RESIDUE`. The report classifies polymer, modified polymer, ligand,
cofactor, ion, water, and common additive groups; lists nearest heavy-atom
neighbors; and highlights metal-coordination candidates. PDB `LINK`/`CONECT` or
mmCIF `struct_conn` counts are surfaced as connectivity evidence, not accepted
as proof of a specific chemistry without inspection.

## Optional structure-quality tools

```bash
python3 scripts/model_quality.py detect --json
python3 scripts/model_quality.py usalign reference.pdb mobile.pdb --json
python3 scripts/model_quality.py dockq model_complex.pdb native_complex.pdb --json
python3 scripts/model_quality.py foldseek query.pdb target_directory --json
```

DockQ writes and parses its JSON result. Foldseek uses `easy-search` with
alignment TM-scores, query/target-normalized TM-scores, lDDT, E-value, and bits.
These tools are optional and are never installed automatically.

## Human handoff checklist

Before handing a result to a user, verify:

1. the input identity, residue numbering, model, conformer, and assembly;
2. every expected image/session exists and is non-empty;
3. the JSON report names the method, selection, cutoff, and provenance;
4. replay scripts contain no credentials and were written only to the requested
   output directory;
5. interpretation distinguishes computed geometry from biological inference;
6. the user receives the `.pml`/`.cxc` and `.pse`/`.cxs` paths needed to continue.
