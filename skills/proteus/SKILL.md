---
name: proteus
description: >
  Use this skill when the user asks you to work with protein structures,
  molecular visualization, or structural biology tools. TRIGGER when:
  the user mentions PyMOL, ChimeraX, AlphaFold, Rosetta, PyRosetta,
  UniProt, RCSB PDB, PDBe, PDB files, protein structures, molecular
  rendering, pLDDT, RMSD, structure alignment, binding pockets,
  drug-target analysis, cryo-EM density maps, homology modeling, or
  protein design. Also trigger when the user opens/loads .pdb, .cif,
  .mmcif, .sdf, or .mol2 files, or .mrc density maps (see
  references/chimerax.md for cryo-EM workflows).
  DO NOT TRIGGER for: general biology questions with no structural component,
  bioinformatics sequence-only tasks (BLAST, MSA), or genomics/transcriptomics.
---

# Proteus — Structural Biology Agent Skill

You are an AI agent driving structural biology tools programmatically.
This skill teaches you how to control PyMOL, ChimeraX, AlphaFold DB,
and Rosetta/PyRosetta from the command line — including the non-obvious
gotchas that will otherwise cost hours of debugging.

## Tool Detection

Before doing anything, detect what's installed:

```python
import shutil, subprocess

PYMOL = shutil.which("pymol")
if not PYMOL:
    # macOS common locations
    import os
    for p in ["/Applications/PyMOL.app/Contents/bin/pymol",
              os.path.expanduser("~/Applications/PyMOL.app/Contents/bin/pymol")]:
        if os.path.isfile(p):
            PYMOL = p
            break

CHIMERAX = shutil.which("ChimeraX") or shutil.which("chimerax")
if not CHIMERAX:
    import glob
    hits = glob.glob("/Applications/ChimeraX*.app/Contents/bin/ChimeraX")
    if hits:
        CHIMERAX = sorted(hits)[-1]  # latest version
```

If neither is found, do not guess paths. Continue with zero-dependency workflows
(`scripts/pdb_info.py`, AlphaFold metadata fetches, file inspection) when they
fit the task; otherwise tell the user what to install and stop.

## Tool Selection — When to Use What

| Task | Best Tool | Why |
|---|---|---|
| Headless rendering (no display) | **PyMOL** | Software ray tracer works fully headless |
| Interactive demo with live GUI | **ChimeraX REST API** | HTTP control of running GUI session |
| H-bonds, SASA, clashes, contacts | **ChimeraX** | Built-in analysis commands, even in `--nogui` |
| Structure alignment + RMSD | **Either** | PyMOL `cealign` or ChimeraX `matchmaker` |
| AlphaFold confidence analysis | **PyMOL** + AlphaFold API | Fetch prediction, color by pLDDT, render headless |
| Experimental PDB download | **`fetch_pdb.py`** | RCSB metadata + coordinates |
| Structure landscape search | **`pdb_search.py`** | Enumerate PDB hits by text and/or UniProt accession |
| Structure candidate ranking | **`pdb_select.py`** | Rank hits by method, resolution, validation, assemblies, and ligands |
| Protein name -> accession | **`uniprot_lookup.py`** | Resolve names/genes before AlphaFold fetch |
| PDB/mmCIF preflight | **`structure_info.py`** | Zero-dependency file inspection |
| Target dossier/report | **`target_dossier.py`** | Markdown + JSON provenance for targets, with opt-in local analyses |
| Structural mutation triage | **`mutation_triage.py`** | Local variant proximity to ligands, interfaces, contacts, and PAE context |
| Docking box prep | **`docking_box.py`** | Ligand-centered Vina box plus optional docking-tool detection |
| Cryo-EM density map visualization | **ChimeraX REST API** | Volume rendering requires GPU/display |
| Quick legacy PDB inspection | **`pdb_info.py` script** | Backward-compatible PDB-only inspector |
| KRAS G12C dossier | **`kras_dossier.py`** | End-to-end workflow with public structures, analyses, and figures |
| Protein design / scoring | **Rosetta/PyRosetta** | Or ML alternatives (ProteinMPNN, RFdiffusion) |

**Key architectural insight:** ChimeraX `--nogui` mode has NO OpenGL context on macOS.
It can run analysis commands (H-bonds, SASA, matchmaker, info) but CANNOT render images.
For ChimeraX rendering, you must use the REST API approach with a running GUI instance.

## Reading Guide

Load reference files on demand — don't read all of them upfront:

| Working with... | Read this file |
|---|---|
| PyMOL (any task) | `references/pymol.md` |
| ChimeraX (any task) | `references/chimerax.md` |
| AlphaFold DB predictions | `references/alphafold.md` |
| PDB/UniProt/PDBe/RCSB data lookup | `references/data-sources.md` |
| File format choices (.pdb, .cif, .sdf, .mrc) | `references/file-formats.md` |
| Prediction models beyond AlphaFold DB | `references/prediction-models.md` |
| Rosetta / protein design | `references/rosetta.md` |

## Agent Helper Scripts

These scripts handle the hard parts of tool communication.

**IMPORTANT: Always run `python3 scripts/<script>.py --help` first.** Treat the
scripts as black-box utilities by default. Only read the source when you are
debugging, patching, or the help text is insufficient for the task.

| Script | Purpose | Example |
|---|---|---|
| `scripts/proteus_doctor.py` | Local readiness report for tools, scripts, and APIs | `python3 scripts/proteus_doctor.py --network --json` |
| `scripts/resolve_structure.py` | Resolve file/PDB/UniProt/name to structure + provenance | `python3 scripts/resolve_structure.py TP53 --json` |
| `scripts/fetch_pdb.py` | RCSB PDB metadata + coordinate fetcher | `python3 scripts/fetch_pdb.py 4HHB --json` |
| `scripts/assembly_report.py` | RCSB biological assembly counts, filenames, URLs, and optional assembly mmCIF download | `python3 scripts/assembly_report.py 4HHB --json` |
| `scripts/pdb_search.py` | Search RCSB PDB by free text and/or UniProt accession | `python3 scripts/pdb_search.py "KRAS G12C" --rows 10 --details --json` |
| `scripts/pdb_select.py` | Rank/select candidate structures from offline metadata or opt-in RCSB lookup | `python3 scripts/pdb_select.py --input candidates.json --json` |
| `scripts/sifts_map.py` | PDBe SIFTS residue mapping between PDB chains and UniProt ranges | `python3 scripts/sifts_map.py pdb 1hsg --json` |
| `scripts/ligand_extract.py` | Zero-dep non-water ligand inventory from PDB/mmCIF and optional CCD downloads | `python3 scripts/ligand_extract.py 1HSG --json` |
| `scripts/dock_prep.py` | Non-executing receptor/ligand prep plan with optional tool detection | `python3 scripts/dock_prep.py --receptor receptor.pdb --ligand ligand.sdf --json` |
| `scripts/docking_box.py` | Ligand-centered PDB/mmCIF Vina-style box and optional docking-tool detection | `python3 scripts/docking_box.py complex.cif --ligand ATP --json` |
| `scripts/dock_vina.py` | Dry-run Vina-compatible docking command planner and log parser | `python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt --json` |
| `scripts/interaction_report.py` | Simple PDB/mmCIF protein-ligand interaction/contact classes plus PLIP/ProLIF detection | `python3 scripts/interaction_report.py complex.cif --ligand ATP --json` |
| `scripts/variant_map.py` | Conservative protein substitution parser plus optional local coordinate lookup | `python3 scripts/variant_map.py "P04637 R175H" --no-download --json` |
| `scripts/mutation_triage.py` | Local structural triage of substitutions against ligands, chains, contacts, and PAE | `python3 scripts/mutation_triage.py R175H --structure model.pdb --json` |
| `scripts/proteus_batch.py` | Run allowlisted helper scripts from a JSON manifest | `python3 scripts/proteus_batch.py manifest.json --dry-run --json` |
| `scripts/proteus_cache.py` | Local URL cache with checksums, offline mode, verify, and garbage collection | `python3 scripts/proteus_cache.py show --json` |
| `scripts/proteus_report.py` | Combine helper JSON into Markdown and scrubbed report JSON evidence packs | `python3 scripts/proteus_report.py --input result.json --json` |
| `scripts/target_dossier.py` | Generic target dossier with Markdown, JSON provenance, and opt-in local analyses | `python3 scripts/target_dossier.py --uniprot P04637 --no-network --json` |
| `scripts/uniprot_lookup.py` | Protein/gene name to UniProt accession | `python3 scripts/uniprot_lookup.py TP53 --gene-exact --json` |
| `scripts/structure_info.py` | Zero-dep PDB/mmCIF inspector | `python3 scripts/structure_info.py structure.cif --json` |
| `scripts/fetch_alphafold.py` | AlphaFold DB fetcher | `python3 scripts/fetch_alphafold.py P04637 --pae --json` |
| `scripts/pae_report.py` | Summarize AlphaFold PAE domain/flexibility hints | `python3 scripts/pae_report.py AF-P04637-F1_pae.json --json` |
| `scripts/validation_report.py` | Fetch wwPDB/RCSB validation quality metrics | `python3 scripts/validation_report.py 4HHB --json` |
| `scripts/pocket_report.py` | Zero-dep ligand pocket contacts from local PDB/mmCIF or PDB ID | `python3 scripts/pocket_report.py complex.cif --json` |
| `scripts/interface_report.py` | Zero-dep protein-protein interface residues from local PDB/mmCIF or PDB ID | `python3 scripts/interface_report.py 1BRS --chains A,D --json` |
| `scripts/compare_structures.py` | PyMOL CE alignment + optional per-residue deviations | `python3 scripts/compare_structures.py ref.pdb mobile.pdb --json` |
| `scripts/pymol_agent.py` | Headless PyMOL driver (info, render, **pocket figure, density fit, spin movie**) | `python3 scripts/pymol_agent.py render structure.pdb out.png --color plddt` |
| `scripts/chimerax_agent.py` | Headless ChimeraX driver (analysis, `--nogui`) | `python3 scripts/chimerax_agent.py run "open 1ubq; info chains #1"` |
| `scripts/chimerax_rest.py` | Managed ChimeraX REST GUI render (GPU) + turntable | `python3 scripts/chimerax_rest.py render structure.pdb out.png --color plddt` |
| `scripts/add_helix_records.py` | Add HELIX records to CA-only backbones so cartoons render | `python3 scripts/add_helix_records.py model.pdb --json` |
| `scripts/map_info.py` | MRC/CCP4 map stats + sigma-based contour levels | `python3 scripts/map_info.py map.mrc --json` |
| `scripts/model_quality.py` | Detect optional quality tools and parse/run USalign when installed | `python3 scripts/model_quality.py detect --json` |
| `scripts/design_run.py` | Validate design manifests and write dry-run modeling/design plans | `python3 scripts/design_run.py manifest.json --dry-run --json` |
| `scripts/rosetta_score.py` | Detect Rosetta/PyRosetta, parse scorefiles, and plan dry scoring commands | `python3 scripts/rosetta_score.py detect --json` |
| `scripts/pdb_info.py` | Legacy zero-dep PDB inspector (PDB only) | `python3 scripts/pdb_info.py structure.pdb` |
| `scripts/kras_dossier.py` | KRAS G12C dossier workflow over public structures | `python3 scripts/kras_dossier.py --out kras_g12c_dossier --no-movie` |

## Critical Gotchas (Read This First)

These are hard-won discoveries. Each one represents hours of debugging that you
can skip by knowing them upfront.

### PyMOL

1. **Never use the `-d` flag for complex commands.** The shell interprets `>`, `<`,
   and `|` in PyMOL selection syntax as redirection/pipe operators. Instead, write
   a `.pml` script file and run `pymol -c -q script.pml`.

   ```python
   # WRONG — breaks on selections like "b > 90"
   subprocess.run(["pymol", "-c", "-q", "-d", 'color blue, b > 90'])

   # RIGHT — write a .pml file
   with open("/tmp/cmd.pml", "w") as f:
       f.write("color blue, b > 90\n")
   subprocess.run(["pymol", "-c", "-q", "/tmp/cmd.pml"])
   ```

   **Why:** The `-d` flag passes the string through the shell, where `>` becomes
   stdout redirection. This is never mentioned in PyMOL documentation.

2. **PyMOL stdout is unreliable.** `print()` in headless mode doesn't always
   capture to stdout. Always write results to a JSON temp file and read it back
   from the calling process.

3. **`cmd.quit()` is required** at the end of every headless script. Without it,
   the PyMOL process hangs indefinitely.

4. **`<=` doesn't exist in PyMOL selection syntax.** For pLDDT coloring, use
   layered overrides — paint the broadest range first, then override with
   narrower selections:
   ```
   color orange, all           # base: everything is low confidence
   color yellow, b > 50        # override: medium
   color cyan, b > 70          # override: high
   color blue, b > 90          # override: very high
   ```

5. **`cealign` argument order is (target, mobile).** The first argument is
   the reference that stays fixed; the second gets moved. This is opposite
   to what you might expect.

6. **`cmd.iterate` vs `cmd.iterate_state`:** Use `iterate` for molecular
   properties (chain, resname, B-factor). Use `iterate_state` for 3D
   coordinates (x, y, z). Using the wrong one silently returns nothing.

7. **GUI demos require threading.** Running a long script in PyMOL's main
   thread freezes the GUI completely. Wrap your demo in
   `threading.Thread(target=run, daemon=True).start()`.

8. **In GUI mode, `print()` goes to PyMOL's internal console**, not the
   terminal. Use `sys.stderr.write(text + "\n")` for terminal output.

9. **`Path(__file__)` is unreliable** inside PyMOL's `-r` runner. Use
   hardcoded absolute paths or resolve paths before launching PyMOL.

### ChimeraX

10. **ChimeraX `--nogui` CANNOT render images on macOS.** There is no
    OpenGL context. Use it only for analysis. For rendering, use the
    REST API with a GUI session (see `references/chimerax.md`).

11. **ChimeraX is NOT thread-safe.** Sending REST API calls from a Python
    background thread causes `EXC_BAD_ACCESS` crashes. All REST calls must
    happen from the main thread.

12. **`close session` does NOT reset model IDs.** After closing and reopening
    structures, model IDs continue incrementing (#1, #2, #3...). The only
    way to reset to #1 is a full process restart. Always use dynamic model
    ID discovery after each `open` command.

13. **Cryo-EM gotchas:** Use `lighting simple` (not `full`) for volume
    maps — `full` washes out colors on white. Never `close` a map
    mid-session — it shifts all model IDs. Use `hide`/`show` instead.
    See `references/chimerax.md` for full cryo-EM patterns.

14. **stdout lines have `INFO:`/`WARNING:`/`ERROR:` prefixes.** Parse them
    out. Also filter lines starting with `INFO: Executing:` — those are
    echoed commands, not results.

### AlphaFold DB

15. **The API returns a list, not a dict.** `json.loads(response)` gives
    `[{...}]`. You must do `data[0]["pdbUrl"]`, not `data["pdbUrl"]`.

16. **Always query `latestVersion` from the API.** Don't hardcode `v4` or
    `v6` in URLs — the version changes and old URLs return 404.

17. **P62988 (ubiquitin) is NOT in the database.** Use P0CG48
    (polyubiquitin-C, 685 residues) instead. Note: this is the full
    polyubiquitin chain, not the 76-residue monomer.

### Rendering, Animation & Maps

18. **PyMOL spin loops need `set cache_frames, 0`.** Otherwise PyMOL caches
    every ray-traced frame in RAM and the process is OOM-killed partway through
    a turntable. Set it before the render loop. (`pymol_agent.py spin` does this
    for you.)

19. **`set auto_zoom, 0` before loading when composing a manual view.** Each
    `load`/`show` otherwise re-zooms and fights your `orient`/`zoom`/`turn`
    framing. Set it first, frame last.

20. **Whole-map isosurfaces stall in headless PyMOL.** The headless build lacks
    the VTKm accelerator, so contouring a full density map can hang for minutes
    even on a small map. Carve the mesh around the model
    (`isomesh m, map, level, sele, carve=2.5`) or do the whole-map surface in
    ChimeraX. Use `scripts/map_info.py` to pick a sigma-based contour level.

21. **ChimeraX REST `save` can return before the PNG is flushed** — a 0-byte
    file, worse under heavy cartoon recompute. After `save`, issue `wait 1` and
    poll the file for non-zero size, retrying a few times.
    (`scripts/chimerax_rest.py` handles this.)

22. **ChimeraX color-name traps.** `gold`/`yellow` atom spheres often render
    visibly green — use `orange` for a true gold read. And `color #1 cartoons #...`
    silently mis-parses `cartoons` as a color: the command is `cartoon`, then
    `color #1 <color>`.

23. **Deposited coordinates are the asymmetric unit, not necessarily the
    biological assembly.** A "dimer" entry may deposit a single chain. Build the
    functional oligomer with ChimeraX `sym #1 assembly 1 copies true`, or expand
    crystal neighbors in PyMOL with `symexp mate_, obj, sele, 5`.

## Common Workflows

### Quick Structure Inspection
```bash
python3 scripts/structure_info.py structure.cif --json  # zero-dep PDB/mmCIF overview
python3 scripts/pdb_info.py structure.pdb               # legacy PDB-only overview
python3 scripts/pymol_agent.py info structure.pdb       # detailed with PyMOL
```

### Experimental Structure Fetch
```bash
python3 scripts/fetch_pdb.py 4HHB --json                # RCSB metadata + mmCIF
python3 scripts/fetch_pdb.py 1HSG --format pdb          # legacy PDB format if needed
python3 scripts/assembly_report.py 4HHB --json          # assembly counts + download URLs
python3 scripts/fetch_pdb.py 4HHB --assembly 1 --json   # biological assembly mmCIF
```

### Structure Landscape Search
```bash
python3 scripts/pdb_search.py "KRAS G12C sotorasib" --rows 10 --details --json
python3 scripts/pdb_search.py --uniprot P01116 --rows 50 --details --json
python3 scripts/pdb_select.py --input candidates.json --ligand GDP --json
```

### Residue Mapping, Variants, And Ligands
```bash
python3 scripts/sifts_map.py pdb 1hsg --json
python3 scripts/variant_map.py "P04637 R175H" --no-download --json
python3 scripts/mutation_triage.py R175H --uniprot P04637 --structure model.pdb --json
python3 scripts/ligand_extract.py 1HSG --json
python3 scripts/dock_prep.py --receptor receptor.pdb --ligand ligand.sdf --json
python3 scripts/docking_box.py complex.pdb --ligand ATP --json
python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt --json
python3 scripts/interaction_report.py complex.pdb --ligand ATP --json
```

### Dossiers And Batch Runs
```bash
python3 scripts/target_dossier.py --gene TP53 --uniprot P04637 --pdb 1TUP --json
python3 scripts/target_dossier.py --pdb structure.pdb --no-network --json
python3 scripts/target_dossier.py --pdb structure.pdb --no-network --analyze-local --json
python3 scripts/proteus_batch.py manifest.json --dry-run --json
python3 scripts/proteus_report.py --input result.json --outdir proteus_report_out --json
python3 scripts/proteus_cache.py verify --json
python3 scripts/design_run.py design_manifest.json --dry-run --json
```

### AlphaFold Confidence Analysis
```bash
python3 scripts/uniprot_lookup.py TP53 --gene-exact --json  # resolve accession if needed
python3 scripts/fetch_alphafold.py P04637 --pae --json       # fetch p53 prediction
# Output filename uses modelEntityId from API, typically AF-{UNIPROT}-F1.pdb
python3 scripts/pae_report.py AF-P04637-F1_pae.json --json   # summarize domain uncertainty
python3 scripts/pymol_agent.py render AF-P04637-F1.pdb output.png
```
Then color by pLDDT bins — see `references/alphafold.md` for the standard color scheme.

### Predicted vs Experimental Comparison
1. Fetch AlphaFold prediction for the protein's UniProt ID
2. Load both structures in PyMOL
3. Use `cealign` (not `align`) — it's purely structural, works for divergent sequences
4. Extract per-residue deviations with `iterate_state`
5. Render side-by-side or overlay

When external quality tools are installed, use `scripts/model_quality.py` to
detect them or parse/run USalign:

```bash
python3 scripts/model_quality.py detect --json
python3 scripts/model_quality.py usalign reference.pdb mobile.pdb --json
python3 scripts/rosetta_score.py detect --json
python3 scripts/rosetta_score.py parse-scorefile score.sc --json
```

### Protein-Protein Interface Analysis
```bash
# Zero-dependency: interface residues between chains, JSON for chaining
python3 scripts/interface_report.py complex.cif --json            # all chain pairs
python3 scripts/interface_report.py 1BRS --chains A,D --cutoff 4.5 --json
python3 scripts/interface_report.py complex.cif --residue A:42 --json
```
For a visual interface dissection (contacts/H-bonds/buried surface), use
ChimeraX — see `references/chimerax.md`.

### Binding Pocket Analysis
```bash
# One-command annotated figure: ligand + pocket sticks + polar contacts + context
python3 scripts/pymol_agent.py pocket 1HSG.pdb pocket.png --label
# Docking-box handoff around known ligand atoms
python3 scripts/dock_prep.py --receptor receptor.pdb --ligand ligand.sdf --json
python3 scripts/docking_box.py 1HSG --ligand MK1 --padding 6 --json
python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --config vina_box.txt --json
python3 scripts/dock_vina.py parse-log vina.log --json
# Contact classes for a ligand group
python3 scripts/interaction_report.py 1HSG --ligand MK1 --json
```
```python
# Or build it by hand. PyMOL: select residues within 5A of any ligand
cmd.select("pocket", "byres organic around 5")
# Show pocket as sticks, ligand as ball-and-stick
cmd.show("sticks", "pocket")
cmd.show("spheres", "organic")
cmd.set("sphere_scale", 0.25, "organic")
```

### ChimeraX Interactive Demo (REST API)
```bash
# 1. Launch ChimeraX with REST
/path/to/ChimeraX --cmd "remotecontrol rest start port 50888" &

# 2. Verify connection
curl -s "http://127.0.0.1:50888/run?command=version"

# 3. Send commands
curl "http://127.0.0.1:50888/run?command=open+1ubq"
curl "http://127.0.0.1:50888/run?command=cartoon"
curl "http://127.0.0.1:50888/run?command=save+/tmp/render.png+width+1200+height+900+supersample+3"
```

### Publication-Quality Rendering (PyMOL)
```
set bg_color, white
set ray_opaque_background, 1
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
ray 1200, 900
png output.png
```

Render presets are also available on the helper: `pymol_agent.py render file.pdb
out.png --preset publication|illustration|soft --color spectrum|chain|bfactor|plddt`.

### Turntable Movie (Headless)
```bash
# PyMOL ray-traces each frame (works with no display); ffmpeg encodes them.
python3 scripts/pymol_agent.py spin structure.pdb spin.mp4 --frames 60 --color plddt
# Degrades gracefully: with no ffmpeg, the frames are written and returned.
```
This is the macOS-correct path — ChimeraX needs a GPU/GUI context to render.

### ChimeraX GPU Rendering (Managed REST)
```bash
# Launches a GUI ChimeraX on an ephemeral port, renders via GPU, tears down.
python3 scripts/chimerax_rest.py render structure.pdb out.png --style surface --color bychain
python3 scripts/chimerax_rest.py spin structure.pdb out.mp4 --frames 72
```
Unlike `chimerax_agent.py` (analysis only, `--nogui`), this renders images and
defeats the 0-byte-PNG save race (gotcha 21).

### CA-only Backbone (de-novo designs)
```bash
# RFdiffusion / Genie backbones render as spaghetti without HELIX records.
python3 scripts/add_helix_records.py model.pdb -o model_ss.pdb --json
# Then render model_ss.pdb; in PyMOL also: set cartoon_trace_atoms, 1
```

### Cryo-EM Density Fit
```bash
# Sigma-based contour level (absolute level differs per map).
python3 scripts/map_info.py map.mrc --json   # -> suggested_level at 1/1.5/2/3 sigma

# Render a model in density — the mesh is carved around the model, which avoids
# the whole-map contour stall in headless PyMOL (gotcha 20).
python3 scripts/pymol_agent.py density model.pdb fit.png --map map.mrc --residue "chain A and resi 25"

# No map yet? Simulate gaussian density from the model (predicted/designed structures).
python3 scripts/pymol_agent.py density model.pdb fit.png --simulate
```

## Good Demo Proteins

| UniProt / PDB | Protein | Good for |
|---|---|---|
| P04637 | p53 | Disorder (pLDDT ~75), mixed confidence regions |
| P69905 | Hemoglobin alpha | Very high confidence (pLDDT ~98) |
| P0CG48 | Polyubiquitin-C (685 res) | Well-folded repeats; use when ubiquitin is requested |
| P01308 | Insulin | Small, well-characterized |
| 1HSG | HIV-1 protease + indinavir | Drug-target binding pocket |
| 1BRS | Barnase-barstar complex | Protein-protein interface |
| 4HHB | Hemoglobin tetramer | Multi-chain, quaternary structure |
| 1EMA | GFP | Chromophore, fluorescence |

## Output Pattern

When running analysis, always produce structured output. The agent scripts
return JSON with `{"status": "ok", "data": {...}}` or `{"status": "error", "error": "..."}` in
`--json` mode. Many scripts also mirror key fields at top level for backward
compatibility, but agents should read `data` first. Their temporary files are
per-invocation, so they are safe to call in parallel within the same workspace.

For multi-step workflows, write a summary JSON report at the end with:
- Input files and parameters
- Key measurements (RMSD, distances, counts, areas)
- Output file paths (renders, saved structures)
- Interpretation notes

## Platform Notes

- **macOS**: PyMOL and ChimeraX are in `/Applications/`. PyMOL headless
  rendering works. ChimeraX `--nogui` has no OpenGL.
- **Linux**: Both tools typically on `$PATH`. PyMOL headless works.
  ChimeraX `--nogui` may have OpenGL via virtual framebuffer (`xvfb-run`).
- **Rosetta**: Requires an academic license. Free alternatives:
  ProteinMPNN (sequence design), LocalColabFold (structure prediction),
  ESM2 (embeddings), RFdiffusion (backbone generation).

## Quick Reference

| I want to... | Do this |
|---|---|
| Check machine readiness | `python3 scripts/proteus_doctor.py --network --json` |
| Resolve any common structure query | `python3 scripts/resolve_structure.py TP53 --json` |
| Resolve protein/gene name to UniProt | `python3 scripts/uniprot_lookup.py TP53 --gene-exact --json` |
| Fetch an experimental PDB structure | `python3 scripts/fetch_pdb.py 4HHB --json` |
| Search experimental structures | `python3 scripts/pdb_search.py "KRAS G12C" --rows 10 --details --json` |
| Rank candidate structures | `python3 scripts/pdb_select.py --input candidates.json --json` |
| Check biological assemblies | `python3 scripts/assembly_report.py 4HHB --json` |
| Map PDB residues to UniProt ranges | `python3 scripts/sifts_map.py pdb 1hsg --json` |
| Inventory bound ligands | `python3 scripts/ligand_extract.py 1HSG --json` |
| Plan receptor/ligand docking prep | `python3 scripts/dock_prep.py --receptor receptor.pdb --ligand ligand.sdf --json` |
| Build a docking box around a ligand | `python3 scripts/docking_box.py complex.pdb --ligand ATP --json` |
| Plan or parse Vina-compatible docking | `python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt --json` |
| Summarize protein-ligand interactions | `python3 scripts/interaction_report.py complex.pdb --ligand ATP --json` |
| Parse/check a protein substitution | `python3 scripts/variant_map.py "P04637 R175H" --no-download --json` |
| Triage mutation structure context | `python3 scripts/mutation_triage.py R175H --structure model.pdb --json` |
| Run a local batch manifest | `python3 scripts/proteus_batch.py manifest.json --dry-run --json` |
| Cache public URL responses locally | `python3 scripts/proteus_cache.py fetch URL --json` |
| Combine JSON outputs into an evidence pack | `python3 scripts/proteus_report.py --input result.json --json` |
| Generate a target dossier | `python3 scripts/target_dossier.py --uniprot P04637 --no-network --json` |
| Generate a local analysis dossier | `python3 scripts/target_dossier.py --pdb structure.cif --no-network --analyze-local --json` |
| Inspect PDB/mmCIF file (no tools needed) | `python3 scripts/structure_info.py file.cif --json` |
| Inspect a legacy PDB file | `python3 scripts/pdb_info.py file.pdb` |
| Summarize AlphaFold PAE | `python3 scripts/pae_report.py AF-P04637-F1_pae.json --json` |
| Fetch validation metrics | `python3 scripts/validation_report.py 4HHB --json` |
| Report ligand pocket contacts | `python3 scripts/pocket_report.py complex.cif --json` |
| Report protein-protein interface residues | `python3 scripts/interface_report.py complex.cif --json` |
| Compare two structures | `python3 scripts/compare_structures.py ref.pdb mobile.pdb --per-residue --json` |
| Detect optional quality tools | `python3 scripts/model_quality.py detect --json` |
| Plan a design run | `python3 scripts/design_run.py design_manifest.json --dry-run --json` |
| Parse Rosetta scorefiles | `python3 scripts/rosetta_score.py parse-scorefile score.sc --json` |
| Get structure info via PyMOL | `python3 scripts/pymol_agent.py info file.pdb` |
| Render a structure headless | `python3 scripts/pymol_agent.py render file.pdb out.png` |
| Render an annotated binding-pocket figure | `python3 scripts/pymol_agent.py pocket file.pdb out.png --label` |
| Render a turntable movie | `python3 scripts/pymol_agent.py spin file.pdb out.mp4` |
| Render via ChimeraX GPU (REST) | `python3 scripts/chimerax_rest.py render file.pdb out.png` |
| Fix a CA-only backbone for cartoons | `python3 scripts/add_helix_records.py model.pdb` |
| Pick a cryo-EM contour level | `python3 scripts/map_info.py map.mrc --json` |
| Render a model in cryo-EM density | `python3 scripts/pymol_agent.py density model.pdb out.png --map map.mrc` |
| Fetch an AlphaFold prediction | `python3 scripts/fetch_alphafold.py UNIPROT_ID --pae --json` |
| Align two structures (ChimeraX) | `python3 scripts/chimerax_agent.py align ref.pdb mobile.pdb` |
| Measure SASA | `python3 scripts/chimerax_agent.py sasa file.pdb` |
| Find H-bonds between chains | `python3 scripts/chimerax_agent.py hbonds file.pdb --chain1 A --chain2 B` |
| Run arbitrary PyMOL commands | `python3 scripts/pymol_agent.py run "fetch 1ubq; show cartoon"` |
| Run arbitrary ChimeraX commands | `python3 scripts/chimerax_agent.py run "open 1ubq; info chains #1"` |
| Control ChimeraX GUI via REST | Read `references/chimerax.md` — REST API section |
| Color by AlphaFold confidence | Read `references/alphafold.md` — pLDDT Coloring section |
| Choose between PDB, mmCIF, SDF, MRC | Read `references/file-formats.md` |
| Use RCSB/PDBe/UniProt APIs | Read `references/data-sources.md` |
| Consider AF3/Boltz/Chai/ColabFold | Read `references/prediction-models.md` |
| Do protein design without Rosetta | Read `references/rosetta.md` — ML Alternatives section |
| Run the KRAS G12C dossier workflow | `python3 scripts/kras_dossier.py --out kras_g12c_dossier --no-movie` |
