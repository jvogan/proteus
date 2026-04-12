# AlphaFold DB Reference — Proteus Skill

## Table of Contents
- [API Overview](#api-overview)
- [Fetching Predictions](#fetching-predictions)
- [Understanding the Output](#understanding-the-output)
- [pLDDT Confidence Coloring](#plddt-confidence-coloring)
- [PAE (Predicted Aligned Error)](#pae-predicted-aligned-error)
- [Comparing Predicted vs Experimental](#comparing-predicted-vs-experimental)
- [Demo Proteins](#demo-proteins)
- [Known Pitfalls](#known-pitfalls)

---

## API Overview

Base URL: `https://alphafold.ebi.ac.uk`

No authentication required. No rate limits documented, but be respectful.

### Endpoints

| Endpoint | Returns |
|---|---|
| `GET /api/prediction/{UNIPROT_ID}` | Metadata JSON (URLs, confidence stats, sequence info) |
| `GET /files/AF-{ID}-F1-model_v{N}.pdb` | Structure file (PDB format) |
| `GET /files/AF-{ID}-F1-model_v{N}.cif` | Structure file (mmCIF format) |
| `GET /files/AF-{ID}-F1-predicted_aligned_error_v{N}.json` | PAE matrix (JSON) |

**Always query `latestVersion` from the metadata API first.** Don't hardcode
the version number — it changes (was v4, now v6 as of 2026) and old URLs
return 404.

## Fetching Predictions

### Quick fetch (command line)
```bash
# Metadata
curl -s "https://alphafold.ebi.ac.uk/api/prediction/P69905" | python3 -m json.tool

# Structure
curl -sL "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.pdb" -o af_hemoglobin.pdb

# PAE matrix
curl -sL "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-predicted_aligned_error_v6.json" -o af_hemoglobin_pae.json
```

### Programmatic fetch (Python, stdlib only)

```python
import urllib.request, json, os

def fetch_alphafold(uniprot_id: str, output_dir: str = ".", include_pae: bool = False):
    """Fetch AlphaFold prediction by UniProt ID. No dependencies."""

    # Step 1: Get metadata (ALWAYS do this to get correct version)
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    with urllib.request.urlopen(api_url) as resp:
        data = json.loads(resp.read())

    # GOTCHA: API returns a list, not a dict
    entry = data[0] if isinstance(data, list) else data

    version = entry.get("latestVersion", 4)
    gene = entry.get("gene", uniprot_id)

    # Step 2: Download structure
    pdb_url = entry["pdbUrl"]
    pdb_path = os.path.join(output_dir, f"AF-{uniprot_id}-{gene}.pdb")
    urllib.request.urlretrieve(pdb_url, pdb_path)

    # Step 3: Download PAE if requested
    pae_path = None
    if include_pae:
        pae_url = entry.get("paeDocUrl")
        if pae_url:
            pae_path = os.path.join(output_dir, f"AF-{uniprot_id}-{gene}_pae.json")
            urllib.request.urlretrieve(pae_url, pae_path)

    # Step 4: Report confidence stats
    stats = {
        "uniprot_id": uniprot_id,
        "gene": gene,
        "version": version,
        "pdb_path": pdb_path,
        "pae_path": pae_path,
        "global_plddt": entry.get("globalMetricValue"),
        "fraction_very_high": entry.get("fractionPlddtVeryHigh"),   # >90
        "fraction_confident": entry.get("fractionPlddtConfident"),  # 70-90
        "fraction_low": entry.get("fractionPlddtLow"),              # 50-70
        "fraction_very_low": entry.get("fractionPlddtVeryLow"),     # <50
    }
    return stats
```

### Using the bundled script

```bash
python scripts/fetch_alphafold.py P04637           # Fetch p53
python scripts/fetch_alphafold.py P69905 --pae     # Fetch hemoglobin with PAE
python scripts/fetch_alphafold.py P04637 --outdir ./data      # Custom output dir
```

## Understanding the Output

### API metadata response

The API returns a **list** containing one entry: `[{...}]`.

Key fields:
```json
{
  "uniprotAccession": "P69905",
  "modelEntityId": "AF-P69905-F1",
  "gene": "HBA1",
  "latestVersion": 6,
  "globalMetricValue": 98.06,
  "fractionPlddtVeryHigh": 0.94,
  "fractionPlddtConfident": 0.06,
  "fractionPlddtLow": 0.0,
  "fractionPlddtVeryLow": 0.0,
  "pdbUrl": "https://...",
  "cifUrl": "https://...",
  "paeDocUrl": "https://...",
  "sequenceStart": 1,
  "sequenceEnd": 142
}
```

Note: `paeDocUrl` may be absent for some entries. Always use `.get()`.

### Structure files

AlphaFold stores **pLDDT confidence** in the B-factor column of PDB files.
Values range from 0 to 100 (not the usual B-factor scale). You can detect
AlphaFold files heuristically: if all B-factors are between 0 and 100 and
the filename contains "AF-", it's likely pLDDT.

## pLDDT Confidence Coloring

### The four bins

| pLDDT | Category | Standard Color | Hex | Meaning |
|---|---|---|---|---|
| > 90 | Very high | Blue | `#0053D6` | High-confidence backbone and sidechain |
| 70-90 | Confident | Cyan | `#65CBF3` | Good backbone, sidechain less certain |
| 50-70 | Low | Yellow | `#FFDB13` | Treat with caution, may be flexible |
| < 50 | Very low | Orange | `#FF7D45` | Likely disordered, don't trust |

### PyMOL discrete coloring

```
# Layer from broadest to narrowest (PyMOL has no <= operator)
color orange, all
color yellow, b > 50
color cyan, b > 70
color blue, b > 90
```

### PyMOL continuous spectrum
```python
cmd.spectrum("b", "blue_white_red", "all", minimum=50, maximum=100)
```

### ChimeraX (built-in palette)
```
color bfactor #1 palette alphafold
```

### Interpretation guidelines

- **> 90 everywhere**: Well-folded, globular domain. Structure is reliable.
- **Mixed 70-90 and > 90**: Core is solid, surface loops less certain. Normal.
- **Regions < 50**: Intrinsically disordered regions (IDRs), flexible linkers,
  or terminal tails. These don't have a single fixed structure.
- **Global pLDDT < 70**: Either a mostly disordered protein, or AlphaFold
  struggled — check if the protein has known structure in PDB.

## PAE (Predicted Aligned Error)

The PAE is an NxN matrix where entry (i,j) represents the expected position
error of residue j when the structure is aligned on residue i.

### How to read it
- **Low PAE (dark blue) between two regions**: AlphaFold is confident about
  their relative positions — they are part of the same rigid domain.
- **High PAE (light/white) between two regions**: Uncertain relative
  orientation — likely different domains connected by a flexible linker.
- **Diagonal is always low**: Each residue is well-positioned relative to itself.
- **Block structure**: Square blocks of low PAE along the diagonal indicate
  domains. Off-diagonal low-PAE blocks indicate inter-domain contacts.

### Detecting domain boundaries from PAE

```python
# Requires: pip install numpy
import json, numpy as np

with open("pae.json") as f:
    pae_data = json.load(f)

# PAE format: list of dicts with "predicted_aligned_error" key
# or direct 2D array depending on version
if isinstance(pae_data, list) and "predicted_aligned_error" in pae_data[0]:
    pae_matrix = np.array(pae_data[0]["predicted_aligned_error"])
else:
    pae_matrix = np.array(pae_data)

# Simple domain boundary detection:
# Average PAE per residue (how uncertain is this residue's position?)
mean_pae_per_residue = pae_matrix.mean(axis=1)
# Residues with mean PAE > 15 are likely in flexible/disordered regions
```

## Comparing Predicted vs Experimental

Standard workflow for AlphaFold validation:

```python
from pymol import cmd

# Load both
cmd.load("AF-P69905-hemoglobin_alpha.pdb", "predicted")
cmd.fetch("1hba", async_=0)  # or load local PDB
cmd.create("experimental", "1hba and chain A")

# Align with cealign (structure-only, best for pred vs exp)
result = cmd.cealign("experimental", "predicted")
if isinstance(result, dict):
    rmsd = result["RMSD"]
    n_aligned = result["alignment_length"]

# Per-residue deviation analysis
exp_coords, pred_coords = {}, {}
cmd.iterate_state(1, "experimental and name CA",
    "exp_coords[int(resi)] = (x, y, z)", space={"exp_coords": exp_coords})
cmd.iterate_state(1, "predicted and name CA",
    "pred_coords[int(resi)] = (x, y, z)", space={"pred_coords": pred_coords})

import math
deviations = {}
for resi in exp_coords:
    if resi in pred_coords:
        dx = exp_coords[resi][0] - pred_coords[resi][0]
        dy = exp_coords[resi][1] - pred_coords[resi][1]
        dz = exp_coords[resi][2] - pred_coords[resi][2]
        deviations[resi] = math.sqrt(dx*dx + dy*dy + dz*dz)
```

Match residues by number (not array index) — after `cealign`, structures
are spatially aligned but residue numbering is unchanged.

## Demo Proteins

Curated set with different confidence profiles for testing:

| UniProt | Protein | pLDDT | Residues | Best for |
|---|---|---|---|---|
| P69905 | Hemoglobin alpha | 98.1 | 142 | Very high confidence baseline |
| P04637 | p53 | 75.1 | 393 | Mixed confidence, disordered regions |
| P01308 | Insulin | 52.9 | 110 | Low confidence (small, disulfide-bonded) |
| P0CG48 | Polyubiquitin-C | ~95 | 685* | Quick test (well-folded monomer) |
| P42212 | GFP | ~95 | 238 | Barrel fold with chromophore |

*P0CG48 is the full polyubiquitin chain (685 residues), not the 76-residue
monomer. The commonly cited ubiquitin ID **P62988 is NOT in the AlphaFold DB**.

## Known Pitfalls

1. **API returns a list `[{...}]`, not a dict.** Always do `data[0]` after
   `json.loads()`. Accessing `data["pdbUrl"]` directly throws `TypeError`.

2. **Version drift.** URLs contain a version (`v4`, `v5`, `v6`). Old versions
   return 404. Always query `latestVersion` from the metadata API first.

3. **P62988 (ubiquitin) is not in the database.** Use P0CG48 instead. But
   be aware P0CG48 is 685 residues (polyubiquitin), not the 76-residue monomer.

4. **`paeDocUrl` may be missing.** Not all entries have PAE data. Always use
   `.get("paeDocUrl")` with a fallback.

5. **AlphaFold 3 has no public API** (as of 2026). AF3 predictions for
   complexes are only available through the AlphaFold Server web interface
   with limited runs per user.

6. **B-factor column interpretation.** In AlphaFold PDB files, the B-factor
   column contains pLDDT (0-100), not actual thermal displacement. Don't
   interpret these as experimental B-factors.
