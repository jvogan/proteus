# ChimeraX Reference — Proteus Skill

## Table of Contents
- [Binary Detection](#binary-detection)
- [CLI Modes](#cli-modes)
- [The No-Rendering Limitation](#the-no-rendering-limitation)
- [REST API (Recommended for Full Control)](#rest-api)
- [Output Parsing](#output-parsing)
- [Specifier Syntax](#specifier-syntax)
- [Core Commands](#core-commands)
- [Analysis Commands](#analysis-commands)
- [AlphaFold Integration](#alphafold-integration)
- [Cryo-EM Density Maps](#cryo-em-density-maps)
- [Model ID Management](#model-id-management)
- [Python API](#python-api)
- [2D Labels for Demos](#2d-labels-for-demos)
- [Headless Capabilities](#headless-capabilities)

---

## Binary Detection

```python
import shutil, os, glob

CHIMERAX = shutil.which("ChimeraX") or shutil.which("chimerax")
if not CHIMERAX:
    # macOS: search /Applications
    hits = glob.glob("/Applications/ChimeraX*.app/Contents/bin/ChimeraX")
    if hits:
        CHIMERAX = sorted(hits)[-1]  # Use latest version
    # Linux: check common paths
    for p in ["/usr/bin/chimerax", "/usr/local/bin/chimerax",
              os.path.expanduser("~/ChimeraX/bin/ChimeraX")]:
        if os.path.isfile(p):
            CHIMERAX = p
            break
```

## CLI Modes

```bash
# Headless analysis (no GUI) — analysis only, NO image rendering
ChimeraX --nogui --exit --cmd "open 1ubq; info chains #1; exit"

# Run a Python script headless
ChimeraX --nogui --exit --script analysis.py

# GUI mode with REST control (RECOMMENDED for full capability)
ChimeraX --cmd "remotecontrol rest start port 50888" &
```

Flags: `--nogui` = no display, `--exit` = quit after commands, `--cmd` = inline
commands (semicolon-separated), `--script` = run Python script.

## The No-Rendering Limitation

**ChimeraX `--nogui` CANNOT render images on macOS.** There is no OpenGL
context. The `save image.png` command will fail silently or produce an error.

**Solutions:**
1. Use **PyMOL** for headless rendering (software ray tracer, works everywhere)
2. Use the **REST API** with a GUI ChimeraX instance (renders via GPU)
3. On Linux, use `xvfb-run ChimeraX --offscreen` (virtual framebuffer)

## REST API

The REST API is the recommended approach when you need both analysis AND
rendering. It drives a running ChimeraX GUI via HTTP.

### Start the REST listener

```bash
# From command line (launches ChimeraX with GUI + REST)
ChimeraX --cmd "remotecontrol rest start port 50888" &

# Or from inside a running ChimeraX session
remotecontrol rest start port 50888
```

### Send commands via HTTP

```bash
# Simple command
curl "http://127.0.0.1:50888/run?command=open+1ubq"

# URL-encode special characters
curl "http://127.0.0.1:50888/run?command=cartoon%3B+color+bychain"

# Get JSON response
curl "http://127.0.0.1:50888/run?command=info+models&json=true"

# Save screenshot (renders via GPU)
curl "http://127.0.0.1:50888/run?command=save+/tmp/render.png+width+1200+height+900+supersample+3"
```

### Python REST client

```python
import urllib.request, urllib.parse, json

CHIMERAX_REST = "http://127.0.0.1:50888"

def cx_run(command: str, json_response: bool = True) -> dict:
    """Send a command to a running ChimeraX via REST."""
    params = urllib.parse.urlencode({
        "command": command,
        "json": "true" if json_response else "false",
    })
    url = f"{CHIMERAX_REST}/run?{params}"
    with urllib.request.urlopen(url) as resp:
        return json.loads(resp.read())
```

### Verify connection before starting work

```python
try:
    cx_run("version")
    print("ChimeraX REST is running")
except Exception:
    print("ChimeraX not running — start it with REST enabled first")
```

### Stop the REST listener

```bash
curl "http://127.0.0.1:50888/run?command=remotecontrol+rest+stop"
```

**CRITICAL: ChimeraX is NOT thread-safe.** Sending REST calls from Python
background threads causes `EXC_BAD_ACCESS` crashes on macOS. All REST calls
must happen from the main thread. If you need concurrency, use `asyncio`
with a single-threaded event loop, not `threading.Thread`.

## Output Parsing

ChimeraX prefixes stdout lines with `INFO:`, `WARNING:`, `ERROR:`, `STATUS:`.

```python
def parse_chimerax_output(stdout: str) -> list[str]:
    """Extract meaningful output lines from ChimeraX stdout."""
    results = []
    for line in stdout.splitlines():
        # Skip echoed commands
        if line.startswith("INFO: Executing:"):
            continue
        # Strip known prefixes
        for prefix in ["INFO: ", "STATUS: ", "WARNING: ", "ERROR: "]:
            if line.startswith(prefix):
                results.append(line[len(prefix):])
                break
    return results
```

### Parsing specific outputs

```python
import re

# RMSD from matchmaker
for line in output_lines:
    m = re.search(r"RMSD between (\d+) .+ is ([\d.]+)", line)
    if m:
        n_atoms, rmsd = int(m.group(1)), float(m.group(2))

# SASA value
for line in output_lines:
    m = re.search(r"area .+ = ([\d.]+)", line)
    if m:
        sasa = float(m.group(1))
```

## Specifier Syntax

ChimeraX uses a hierarchical specifier system, different from PyMOL:

```
#1                   # Model 1
#1/A                 # Chain A of model 1
#1/A:50              # Residue 50, chain A, model 1
#1/A:50-100          # Residue range
#1/A:50@CA           # Atom CA of residue 50
#1 & protein         # Protein atoms of model 1
#1 & ligand          # Ligand (non-polymer, non-water)
#1 & solvent         # Water molecules
#1 & ~protein        # Everything except protein
```

## Core Commands

### Opening structures
```
open 1ubq from pdb               # Fetch from PDB
open /path/to/file.pdb           # Local file (use absolute paths!)
open P69905 from alphafold       # AlphaFold prediction
open 21924 from emdb             # Cryo-EM density map
```

**Always use absolute paths** for local files. ChimeraX's working directory
when launched via subprocess may differ from yours.

### Visualization
```
cartoon                           # Show cartoon ribbon
hide                              # Hide all
style #1 stick                    # Stick representation
surface #1                        # Molecular surface
~surface                          # Remove surface
color #1 red                      # Color model red
color bychain #1                  # Color by chain
color bfactor #1                  # Color by B-factor
color bfactor #1 palette alphafold  # AlphaFold confidence colors
transparency #1 50                # 50% transparent
transparency #1 50 target c       # Cartoon only (target: c=cartoon, s=surface, a=atoms)
```

### Presets
```
preset interactive                # Quick visualization
preset publication 1              # Publication style 1
preset publication 2              # Publication style 2 (different lighting)
```

### Export
```
save output.pdb #1                # Save structure
save output.cif #1                # mmCIF format
save session.cxs                  # Save session
save output.png width 1200 height 900 supersample 3   # Screenshot (GUI only!)
```

### Lighting
```
lighting simple                   # Best for density maps + models
lighting soft                     # Good for proteins alone
lighting full                     # AVOID with density maps (washes out colors)
```

## Analysis Commands

```
# Distances
distance #1/A:50@CA #1/A:55@CA

# SASA
measure sasa #1                   # Total SASA
measure sasa #1/A                 # Per-chain SASA

# H-bonds between chains
hbonds #1/A restrict #1/B         # H-bonds between chains A and B
hbonds #1/A restrict #1/B reveal true log true   # Show + log to stdout

# Steric clashes
clashes #1/A restrict #1/B

# Interface contacts
interfaces #1/A contacts #1/B

# Structure alignment
matchmaker #2 to #1               # Align model 2 onto model 1
align #1 to #2                    # Superposition (pre-aligned sequences)
rmsd #1 #2                        # Compute RMSD

# Zone-based selection (binding pocket)
select zone #1/A:LIG 5 #1/A & protein residues true
# Translation: select all protein residues within 5A of ligand LIG in chain A
```

Key flags for hbonds: `reveal true` makes interacting atoms visible as sticks.
`log true` echoes results to stdout for capture.

## AlphaFold Integration

ChimeraX has built-in AlphaFold support:

```
open P69905 from alphafold        # Fetch prediction
alphafold fetch P69905            # Alternative syntax
color bfactor #1 palette alphafold  # Standard confidence coloring
alphafold pae #1                  # Show PAE heatmap
alphafold contacts #1             # Predicted contacts overlay
```

## Cryo-EM Density Maps

```
# Fetch from EMDB
open 21924 from emdb              # Downloads .map file

# Display
volume #1 style surface level 0.025 color #90CAF9 transparency 0.3 step 1

# Fit atomic model into map
fitmap #2 inMap #1 metric correlation

# Map statistics
measure mapstats #1

# Adjust contour level
volume #1 level 0.05              # Higher = less volume shown
```

**Use `lighting simple`** when combining density maps with atomic models.
`lighting full` makes volume colors look washed out on white backgrounds.

## Model ID Management

**This is the most important ChimeraX gotcha for agent workflows.**

Model IDs are global session state. They increment monotonically and
**`close session` does NOT reset them**. If you open model #1, close the
session, then open a new model, it becomes #2 — not #1.

### Dynamic model ID discovery

```python
def get_model_id(name_fragment: str, model_type: str = "AtomicStructure") -> str:
    """Find the model ID for a recently opened structure."""
    result = cx_run(f"info models")
    for line in parse_output(result):
        if name_fragment.lower() in line.lower() and model_type in line:
            # Parse model ID from the info output
            match = re.search(r"#(\d+)", line)
            if match:
                return match.group(0)  # e.g., "#3"
    return None
```

### Clean restart pattern

When you need guaranteed clean model IDs (starting from #1):

```python
import subprocess, time

def restart_chimerax_with_rest(port=50888):
    """Kill and restart ChimeraX with clean model IDs."""
    subprocess.run(["pkill", "-f", "ChimeraX"], capture_output=True)
    time.sleep(2)
    subprocess.Popen([CHIMERAX, "--cmd", f"remotecontrol rest start port {port}"])
    # Poll until REST responds
    for _ in range(30):
        try:
            cx_run("version")
            return True
        except Exception:
            time.sleep(1)
    return False
```

On Windows, replace `pkill` with `taskkill /f /im ChimeraX.exe`.

### Avoid closing models mid-session

Closing a model shifts all subsequent model IDs. Instead of `close #2`,
use `hide #2` to make it invisible. Use `show #2` to bring it back.

## Python API

For scripts run with `--script`:

```python
from chimerax.core.commands import run as rc

# Run any command
rc(session, "open 1ubq from pdb")
rc(session, "matchmaker #2 to #1")
rc(session, "hbonds #1/A restrict #1/B log true")

# Access model data
models = session.models.list()
for m in models:
    print(f"{m.name}: {m.num_atoms} atoms, {m.num_residues} residues")

# Access atoms
atoms = models[0].atoms
coords = atoms.coords        # Nx3 numpy array
bfactors = atoms.bfactors    # B-factors / pLDDT
```

## 2D Labels for Demos

Overlay text labels on the ChimeraX viewport (REST API or command):

```
# Create a label
2dlabels create name title text "Step 1: Load Structure" xpos 0.05 ypos 0.95 color black size 24 bold true

# Update it
2dlabels change name title text "Step 2: Analyze"

# Delete all labels
2dlabels delete

# Multiple labels
2dlabels create name step text "Step 1" xpos 0.05 ypos 0.95 color white size 20 bold true
2dlabels create name detail text "Loading hemoglobin" xpos 0.05 ypos 0.90 color gray size 14
```

Coordinates: `xpos`/`ypos` are 0.0-1.0 fractions of viewport.

## Spinning animation

```python
def spin(axis="y", degrees=360, frames=90):
    """Rotate the view. Must sleep to let animation complete."""
    deg_per_frame = degrees / frames
    cx_run(f"turn {axis} {deg_per_frame} {frames}")
    time.sleep(frames * 0.033)  # 30 FPS
```

ChimeraX `turn` commands execute asynchronously. Without the sleep, the
next command runs before the rotation finishes.

## `.cxc` Script Timing

In ChimeraX command scripts (`.cxc` files), `wait N` pauses for N frames.
At 30 FPS: `wait 30` = 1 second, `wait 90` = 3 seconds.

## Headless Capabilities

| Capability | `--nogui` | REST API (GUI) |
|---|---|---|
| Open/fetch structures | Yes | Yes |
| Structure alignment | Yes | Yes |
| Measurements (dist, SASA) | Yes | Yes |
| H-bonds / clashes | Yes | Yes |
| Interface analysis | Yes | Yes |
| AlphaFold fetch + PAE | Yes | Yes |
| Density map fitting | Yes | Yes |
| Save PDB/CIF | Yes | Yes |
| **Image rendering** | **No** | **Yes** |
| 2D labels / annotations | No | Yes |
| Volume visualization | No | Yes |
| Sessions (.cxs) | Yes | Yes |
