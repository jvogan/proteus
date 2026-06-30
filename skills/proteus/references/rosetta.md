# Rosetta / PyRosetta Reference — Proteus Skill

## Table of Contents
- [Installation](#installation)
- [Core Capabilities](#core-capabilities)
- [Score Functions](#score-functions)
- [Common Protocols](#common-protocols)
- [Score File Parsing](#score-file-parsing)
- [ML Alternatives](#ml-alternatives)
- [Known Pitfalls](#known-pitfalls)

---

## Installation

Rosetta requires an **academic license** (free for academic use, commercial
licenses available). PyRosetta has native Apple Silicon (ARM64) builds.

```bash
# Method 1: pip installer (simplest)
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

# Method 2: conda
conda install -y -c https://conda.rosettacommons.org -c conda-forge pyrosetta

# Verify
python -c "import pyrosetta; pyrosetta.init('-mute all'); print('OK')"
```

If Rosetta/PyRosetta is not installed, you can still:
- Parse Rosetta score files (`.sc`) with pandas
- Analyze designed PDB files with PyMOL/ChimeraX
- Use ML alternatives (ProteinMPNN, ESM2) for sequence design

## Core Capabilities

| Task | PyRosetta Function | Notes |
|---|---|---|
| Score a structure | `sfxn(pose)` | Returns Rosetta Energy Units (REU) |
| Energy minimization | `FastRelax` | Relaxes structure to local minimum |
| Point mutation ddG | `mutate_residue` + rescore | Positive = destabilizing |
| Protein-protein docking | `DockingProtocol` | Rigid-body docking |
| Sequence design | `PackRotamersMover` | Fixed backbone design |
| Loop modeling | `LoopModeler` | Ab initio loop building |
| Homology modeling | `comparative_modeling` | Template-based |

## Score Functions

```python
import pyrosetta
pyrosetta.init("-mute all")

# Load structure
pose = pyrosetta.pose_from_pdb("protein.pdb")

# Standard score function (REF2015)
sfxn = pyrosetta.create_score_function("ref2015")
total_score = sfxn(pose)  # Returns float (REU)

# Per-residue scores
for i in range(1, pose.total_residue() + 1):
    energies = pose.energies()
    res_score = energies.residue_total_energy(i)
```

**Key score terms:**
- `fa_atr` — van der Waals attractive
- `fa_rep` — van der Waals repulsive (steric clashes)
- `fa_sol` — solvation energy
- `hbond_sr_bb`, `hbond_lr_bb`, `hbond_bb_sc`, `hbond_sc` — H-bond terms
- `rama_prepro` — backbone Ramachandran
- `fa_elec` — electrostatics

## Common Protocols

### Energy minimization (FastRelax)
```python
from pyrosetta.rosetta.protocols.relax import FastRelax

relax = FastRelax()
relax.set_scorefxn(sfxn)
relax.apply(pose)
pose.dump_pdb("relaxed.pdb")
```

### Point mutation ddG
```python
from pyrosetta.toolbox.mutants import mutate_residue
from pyrosetta import Pose

# Clone the pose
mutant = Pose(pose)

# Mutate residue 50 (Rosetta numbering) to Alanine
pdb_resnum = 50
chain = "A"
rosetta_resnum = pose.pdb_info().pdb2pose(chain, pdb_resnum)
mutate_residue(mutant, rosetta_resnum, "ALA")

# Score both and compute ddG
ddg = sfxn(mutant) - sfxn(pose)
# Positive ddG = destabilizing mutation
# Negative ddG = stabilizing mutation
```

### Protein-protein docking
```python
from pyrosetta.rosetta.protocols.docking import setup_foldtree, DockingProtocol
from pyrosetta.rosetta.utility import Vector1

setup_foldtree(pose, "A_B", Vector1([1]))
docking = DockingProtocol()
docking.apply(pose)
```

## Score File Parsing

Rosetta score files (`.sc`) are whitespace-delimited text files. They have
a header row preceded by `SCORE:` labels.

```python
import pandas as pd

def parse_score_file(path: str) -> pd.DataFrame:
    """Parse a Rosetta score file into a DataFrame."""
    df = pd.read_csv(path, sep=r'\s+')
    # The first column is literally "SCORE:" — a row label, not data. Drop it.
    if 'SCORE:' in df.columns:
        df = df.drop(columns=['SCORE:'])
    return df

# Usage
df = parse_score_file("score.sc")
best = df.sort_values("total_score").head(10)
print(best[["description", "total_score", "fa_rep", "hbond_sr_bb"]])
```

## ML Alternatives

When Rosetta is not installed or you want faster results, these ML tools
cover many of the same tasks:

| Tool | Task | macOS ARM64 | License | Install |
|---|---|---|---|---|
| **ProteinMPNN** | Sequence design | Yes (CPU/MPS) | MIT | `pip install protein-mpnn-pip` |
| **LocalColabFold** | Structure prediction | Yes (CPU, slow) | Apache 2.0 | Complex setup |
| **ESM2** | Embeddings, contacts | Yes | MIT | `pip install fair-esm` |
| **RFdiffusion** | Backbone generation | No (GPU only) | BSD | Requires CUDA |
| **ESMFold** | Fast structure prediction | Yes (CPU) | MIT | Via `fair-esm` |

### ProteinMPNN (recommended for sequence design)
```bash
pip install protein-mpnn-pip

# Design 10 sequences for a backbone
protein_mpnn_run.py \
    --pdb_path scaffold.pdb \
    --out_folder output/ \
    --num_seq_per_target 10
```

### ESM2 (embeddings and zero-shot mutation effects)
```python
import torch, esm

# Load pre-trained model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.requires_grad_(False)  # Inference mode

data = [("protein", "MKTAYIAKQRQISFVKSH...")]
batch_labels, batch_strs, batch_tokens = batch_converter(data)

with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33])
    embeddings = results["representations"][33]  # Per-residue embeddings
```

## Known Pitfalls

1. **Rosetta Energy Units (REU) are NOT kcal/mol.** They are unitless scores
   from a statistical potential. Only compare REU values computed with the
   same score function. Do not mix ref2015 scores with beta_nov16 scores.

2. **Score file `SCORE:` column.** The first column is literally the string
   `SCORE:` — it's a row label, not data. Drop it before analysis.

3. **Rosetta numbering vs PDB numbering.** Rosetta renumbers residues starting
   from 1. Use `pose.pdb_info().pdb2pose(chain, pdb_resnum)` to convert
   PDB residue numbers to Rosetta numbering.

4. **FastRelax changes coordinates.** After relaxation, the structure may
   have moved significantly from the input. Always save the relaxed pose
   and compare with the original.

5. **Academic license required.** Rosetta is free for academic use but
   requires registration. Commercial use requires a separate license.
   ML alternatives (ProteinMPNN, ESM2) have permissive open-source licenses.
