# PyMOL Reference — Proteus Skill

## Table of Contents
- [Binary Detection](#binary-detection)
- [CLI Modes](#cli-modes)
- [Output Capture Pattern](#output-capture-pattern)
- [Selection Algebra](#selection-algebra)
- [Visualization Commands](#visualization-commands)
- [Alignment](#alignment)
- [Measurement & Analysis](#measurement--analysis)
- [Rendering (Headless)](#rendering-headless)
- [Publication Settings](#publication-settings)
- [Python API (cmd module)](#python-api-cmd-module)
- [AlphaFold pLDDT Coloring](#alphafold-plddt-coloring)
- [GUI Demo Patterns](#gui-demo-patterns)
- [Headless Capabilities](#headless-capabilities)

---

## Binary Detection

```python
import shutil, os, glob

PYMOL = shutil.which("pymol")
if not PYMOL:
    for p in ["/Applications/PyMOL.app/Contents/bin/pymol",
              os.path.expanduser("~/Applications/PyMOL.app/Contents/bin/pymol"),
              "/usr/bin/pymol", "/usr/local/bin/pymol"]:
        if os.path.isfile(p):
            PYMOL = p
            break
# On Linux, pymol is usually on PATH via package manager
```

## CLI Modes

```bash
# Headless (no GUI) — for batch processing and rendering
pymol -c -q -r script.py           # Run Python script
pymol -c -q script.pml             # Run PyMOL command script

# GUI mode — for interactive demos
pymol script.py                     # Opens GUI, runs script
pymol -r script.py                  # Same, explicit -r flag
```

Flags: `-c` = command-line only (no GUI), `-q` = quiet (no splash), `-r` = run script.

**NEVER use `-d` for complex commands.** The shell interprets `>`, `<`, `|` in
selections as redirection operators. Always write a `.pml` or `.py` file instead.

## Output Capture Pattern

PyMOL's stdout is unreliable in headless mode. Use this pattern to get
structured data back from PyMOL scripts:

```python
import subprocess, json, tempfile, os

PYMOL_BIN = "/path/to/pymol"  # Detected above
OUTPUT_FILE = "/tmp/pymol_agent_output.json"

def run_pymol_script(script_content: str, timeout: int = 120) -> dict:
    """Run PyMOL Python script headlessly, return structured JSON."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        wrapper = f'''
import json, sys
_output = {{"status": "ok", "data": {{}}}}
try:
    from pymol import cmd
    # --- your script content goes here, indented ---
    _output["data"]["example"] = "value"
except Exception as e:
    _output["status"] = "error"
    _output["error"] = str(e)
finally:
    with open("{OUTPUT_FILE}", "w") as _f:
        json.dump(_output, _f, indent=2, default=str)
    try:
        cmd.quit()
    except Exception:
        pass
'''
        f.write(wrapper)
        script_path = f.name

    try:
        subprocess.run([PYMOL_BIN, "-c", "-q", "-r", script_path],
                       capture_output=True, text=True, timeout=timeout)
        if os.path.exists(OUTPUT_FILE):
            with open(OUTPUT_FILE) as f:
                return json.load(f)
        return {"status": "error", "error": "No output file produced"}
    finally:
        os.unlink(script_path)
        if os.path.exists(OUTPUT_FILE):
            os.unlink(OUTPUT_FILE)
```

**Why this pattern?** `print()` inside PyMOL headless mode gets swallowed or
mixed with PyMOL's own diagnostic output. Writing to a temp JSON file and
reading it back from the calling process is the only reliable method.

Always call `cmd.quit()` at the end — without it, the PyMOL process hangs.

## Selection Algebra

```
# Atoms
chain A                       # Chain A
resi 50                       # Residue 50
resi 50-100                   # Range
resn ALA                      # All alanines
name CA                       # Alpha carbons
name CA+CB+N+C+O              # Backbone + CB

# Logical
chain A and resi 50-100       # Intersection
chain A or chain B            # Union
not polymer                   # Complement

# Proximity
organic                       # Ligands (non-polymer, non-solvent)
solvent                       # Water molecules
polymer                       # Protein + nucleic acid
br. resi 50                   # Byres: complete residue around selection
byres organic around 5        # Residues within 5A of any ligand (binding pocket)
within 4 of organic           # Atoms within 4A of ligands

# Specifying objects
/object_name/chain/resi/atom  # Full specifier
/1ubq/A/50/CA                 # CA of residue 50, chain A, object 1ubq
```

Note: `<=` does not exist in PyMOL selection syntax. Use `b > 50` (greater than)
only. For ranges, use layered coloring (see pLDDT section).

## Visualization Commands

```
hide everything               # Clear all representations
show cartoon                  # Ribbon diagram
show sticks, resi 50-55       # Ball-and-stick for specific residues
show surface                  # Molecular surface
show spheres, name ZN         # Space-filling for metals
show lines, organic           # Wireframe for ligands

# Coloring
color red, chain A
color blue, chain B
spectrum count, rainbow       # Rainbow by sequence position
spectrum b, blue_white_red    # Color by B-factor/pLDDT
util.cbc()                    # Color by chain (Python API)
cmd.do("util.cnc organic")   # Color non-carbon by element (must use cmd.do)

# Transparency
set cartoon_transparency, 0.5
set surface_transparency, 0.3
set stick_transparency, 0.2

# Styles
cartoon putty                 # Tube width proportional to B-factor
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
```

Note: `util.cnc` (color non-carbon by element) must be called via `cmd.do("util.cnc selection")`.
It is a PyMOL built-in script function, not a Python method on the `cmd` object.

## Alignment

Three methods with different strengths:

```python
# 1. align — uses sequence + structure
rmsd = cmd.align("mobile", "target")  # Returns (RMSD, n_atoms_aligned, ...)

# 2. super — structure only superposition
rmsd = cmd.super("mobile", "target")  # Similar return format

# 3. cealign — CE algorithm, best for divergent sequences
result = cmd.cealign("target", "mobile")  # NOTE: target FIRST, mobile SECOND
# Returns dict: {"RMSD": float, "alignment_length": int, ...}
# Guard with: isinstance(result, dict) — returns -1 on failure
```

**When to use what:**
- `align` — same protein, different conformations (crystal forms, MD snapshots)
- `super` — related proteins, reasonable sequence similarity
- `cealign` — AlphaFold vs experimental, cross-species comparisons, no sequence needed

**Critical:** `cealign` argument order is **(target, mobile)** — the first argument
stays fixed, the second gets moved. This is opposite to `align` and `super`.

## Measurement & Analysis

```
# Distances
distance d1, /obj/A/50/CA, /obj/A/55/CA   # Named distance object
cmd.get_distance("/obj/A/50/CA", "/obj/A/55/CA")  # Python, returns float

# Polar contacts (H-bonds) between selection and ligand
distance hb, organic, pocket, 3.5, 2       # mode 2 = polar contacts only

# Angles
angle a1, sel1, sel2, sel3

# Binding pocket selection
select pocket, byres organic around 5       # Residues within 5A of ligand
```

### Data extraction with iterate

```python
# Properties (chain, residue, B-factor, etc.)
stored_data = []
cmd.iterate("name CA", "stored_data.append((chain, resi, resn, b))",
            space={"stored_data": stored_data})

# 3D coordinates — use iterate_state, NOT iterate
coords = []
cmd.iterate_state(1, "name CA", "coords.append((x, y, z))",
                  space={"coords": coords})
```

**`iterate` gives properties. `iterate_state` gives coordinates.** Using
the wrong one silently returns nothing.

## Rendering (Headless)

PyMOL's software ray tracer works without any display. This is the primary
advantage over ChimeraX for agent workflows.

```
# Basic render
orient                        # Auto-orient view
ray 1200, 900                 # Ray-trace at resolution
png output.png                # Save to file

# Combined (slightly faster)
png output.png, ray=1, width=1200, height=900

# View control
orient                        # Auto-orient to show all
zoom selection                # Zoom to selection
center selection              # Center on selection
turn y, 45                    # Rotate 45 degrees around Y
```

## Publication Settings

Standard block for high-quality renders:

```
set bg_color, white
set ray_opaque_background, 1
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
set ray_shadows, 1
set specular, 0.3
set ambient, 0.4
```

Ray trace modes:
- `set ray_trace_mode, 0` — normal (default, photorealistic)
- `set ray_trace_mode, 1` — quantized colors (illustration style)
- `set ray_trace_mode, 3` — quantized + black outlines (figure style)

## Python API (cmd module)

```python
from pymol import cmd

# Load
cmd.fetch("1ubq", async_=0)       # Note: async_ with underscore
cmd.load("structure.pdb", "myobj")

# Query
names = cmd.get_names()            # List of object names
count = cmd.count_atoms("chain A") # Atom count
chains = cmd.get_chains("all")     # List of chains

# Object manipulation
cmd.create("chain_a", "chain A")   # New object from selection
cmd.extract("lig", "organic")      # Extract to new object
cmd.split_chains("all")            # Split into per-chain objects

# Scenes (save/recall camera views)
cmd.scene("overview", "store")     # Save current view
cmd.scene("overview", "recall")    # Restore view

# Labels
cmd.pseudoatom("label_point", pos=[10.0, 20.0, 30.0])  # 3D label anchor
cmd.label("label_point", '"Active Site"')                 # Label text
cmd.set("label_size", 14)
cmd.set("label_color", "black")
```

## AlphaFold pLDDT Coloring

AlphaFold stores pLDDT (0-100) in the B-factor column. Two coloring approaches:

### Continuous spectrum
```python
cmd.spectrum("b", "blue_white_red", "all", minimum=50, maximum=100)
```

### Discrete bins (official AlphaFold scheme)
```
# MUST apply in this order — broadest first, narrowest last
color orange, all              # <50: very low (base layer)
color yellow, b > 50           # 50-70: low
color cyan, b > 70             # 70-90: confident
color blue, b > 90             # >90: very high (top layer)
```

**Why layered?** PyMOL's selection algebra doesn't support `<=`. You can't write
`color cyan, b > 70 and b <= 90` reliably. Instead, paint everything orange first,
then override progressively.

Standard hex colors: Blue `#0053D6`, Cyan `#65CBF3`, Yellow `#FFDB13`, Orange `#FF7D45`.

## GUI Demo Patterns

When running a live demo with PyMOL's GUI visible:

```python
import threading, sys, time
from pymol import cmd

def run_demo():
    # All demo logic here
    cmd.fetch("1ubq", async_=0)
    cmd.show("cartoon")
    cmd.refresh()        # Force GUI redraw after each change
    time.sleep(2)        # Pause for audience
    sys.stderr.write("Step 1 complete\n")  # Terminal output (not PyMOL console)

# MUST run in a daemon thread — main thread runs PyMOL's event loop
t = threading.Thread(target=run_demo, daemon=True)
t.start()
```

**Why threading?** `time.sleep()` in the main thread freezes the entire PyMOL GUI.
The daemon thread lets the GUI stay responsive while the demo runs.

**Why stderr?** In GUI mode, `print()` goes to PyMOL's internal console widget,
not the terminal. `sys.stderr.write()` always goes to the terminal.

## Headless Capabilities

| Capability | Headless? | Notes |
|---|---|---|
| Load/fetch structures | Yes | `fetch`, `load` |
| Analysis (iterate, distances) | Yes | Full Python API |
| Alignment (align, super, cealign) | Yes | Returns RMSD |
| **Image rendering** | **Yes** | `ray` + `png` (software renderer) |
| Save PDB/CIF | Yes | `save` command |
| Sessions | Yes | `save session.pse` |
| APBS electrostatics | Yes | If APBS installed |
| Morph/interpolation | Yes | `morph` command |
