---
name: proteus
license: MIT
description: >
  Use this skill when the user asks you to work with protein structures,
  molecular visualization, or structural biology tools. TRIGGER when:
  the user mentions PyMOL, ChimeraX, AlphaFold, Rosetta, PyRosetta,
  PDB files, protein structures, molecular rendering, pLDDT, RMSD,
  structure alignment, binding pockets, drug-target analysis, cryo-EM
  density maps, homology modeling, or protein design. Also trigger when
  the user opens/loads .pdb, .cif, .sdf, or .mol2 files, or .mrc
  density maps (see references/chimerax.md for cryo-EM workflows).
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
| Cryo-EM density map visualization | **ChimeraX REST API** | Volume rendering requires GPU/display |
| Quick PDB file inspection | **`pdb_info.py` script** | Zero dependencies, instant |
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
| Rosetta / protein design | `references/rosetta.md` |

## Agent Helper Scripts

These scripts handle the hard parts of tool communication.

**IMPORTANT: Always run `python scripts/<script>.py --help` first.** Treat the
scripts as black-box utilities by default. Only read the source when you are
debugging, patching, or the help text is insufficient for the task.

| Script | Purpose | Example |
|---|---|---|
| `scripts/pymol_agent.py` | Headless PyMOL driver | `python scripts/pymol_agent.py info structure.pdb` |
| `scripts/chimerax_agent.py` | Headless ChimeraX driver | `python scripts/chimerax_agent.py run "open 1ubq; info chains #1"` |
| `scripts/fetch_alphafold.py` | AlphaFold DB fetcher | `python scripts/fetch_alphafold.py P04637 --pae` |
| `scripts/pdb_info.py` | Zero-dep PDB inspector (PDB format only, not CIF) | `python scripts/pdb_info.py structure.pdb` |

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

## Common Workflows

### Quick Structure Inspection
```bash
python scripts/pdb_info.py structure.pdb          # zero-dep overview
python scripts/pymol_agent.py info structure.pdb   # detailed with PyMOL
```

### AlphaFold Confidence Analysis
```bash
python scripts/fetch_alphafold.py P04637 --pae     # fetch p53 prediction
# Output filename uses modelEntityId from API, typically AF-{UNIPROT}-F1.pdb
python scripts/pymol_agent.py render AF-P04637-F1.pdb output.png
```
Then color by pLDDT bins — see `references/alphafold.md` for the standard color scheme.

### Predicted vs Experimental Comparison
1. Fetch AlphaFold prediction for the protein's UniProt ID
2. Load both structures in PyMOL
3. Use `cealign` (not `align`) — it's purely structural, works for divergent sequences
4. Extract per-residue deviations with `iterate_state`
5. Render side-by-side or overlay

### Binding Pocket Analysis
```python
# PyMOL: select residues within 5A of any ligand
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
return JSON with `{"status": "ok", "data": {...}}` or `{"status": "error", "error": "..."}`.
Their temporary files are per-invocation, so they are safe to call in parallel
within the same workspace.

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
| Inspect a PDB file (no tools needed) | `python scripts/pdb_info.py file.pdb` |
| Get structure info via PyMOL | `python scripts/pymol_agent.py info file.pdb` |
| Render a structure headless | `python scripts/pymol_agent.py render file.pdb out.png` |
| Fetch an AlphaFold prediction | `python scripts/fetch_alphafold.py UNIPROT_ID --pae` |
| Align two structures (ChimeraX) | `python scripts/chimerax_agent.py align ref.pdb mobile.pdb` |
| Measure SASA | `python scripts/chimerax_agent.py sasa file.pdb` |
| Find H-bonds between chains | `python scripts/chimerax_agent.py hbonds file.pdb --chain1 A --chain2 B` |
| Run arbitrary PyMOL commands | `python scripts/pymol_agent.py run "fetch 1ubq; show cartoon"` |
| Run arbitrary ChimeraX commands | `python scripts/chimerax_agent.py run "open 1ubq; info chains #1"` |
| Control ChimeraX GUI via REST | Read `references/chimerax.md` — REST API section |
| Color by AlphaFold confidence | Read `references/alphafold.md` — pLDDT Coloring section |
| Do protein design without Rosetta | Read `references/rosetta.md` — ML Alternatives section |
