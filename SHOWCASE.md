# Proteus Showcase

Reproducible figures and analyses you can generate with Proteus, using only
public structures. Each entry lists the exact commands and what you should see.

Outputs (`.pdb`, `.png`, `.mp4`) are git-ignored, so running these from the repo
root never dirties the repository. Every command here is a real code path in the
helper scripts, so this doubles as a set of living examples.

> Requires Python 3.10+. Rendering entries need PyMOL (or ChimeraX for the GPU
> path); analysis entries are stdlib-only. Run `python3 scripts/proteus_doctor.py`
> to see what's available on your machine. Turntable movies additionally need
> `ffmpeg`.

## 1. Binding pocket — HIV-1 protease + indinavir (1HSG)

```bash
python3 scripts/fetch_pdb.py 1HSG --format pdb
python3 scripts/pymol_agent.py pocket 1HSG.pdb hiv_protease_pocket.png --label
```

**Shows:** the indinavir inhibitor nested in the protease active site, with the
surrounding pocket residues as labeled sticks and polar contacts dashed. The
JSON reports the ligand atom count and number of pocket residues.

## 2. Protein–protein interface — SARS-CoV-2 RBD + ACE2 (6M0J)

```bash
python3 scripts/fetch_pdb.py 6M0J --format pdb
python3 scripts/interface_report.py 6M0J.pdb --chains A,E --cutoff 4.5 --json
python3 scripts/pymol_agent.py render 6M0J.pdb ace2_rbd.png --color chain
```

**Shows:** the interface residues on each side of the ACE2 (A) / spike RBD (E)
contact — the hotspot of SARS-CoV-2 recognition — plus a chain-colored render of
the complex.

## 3. Quaternary structure — hemoglobin tetramer (4HHB)

```bash
python3 scripts/fetch_pdb.py 4HHB --format pdb
python3 scripts/pymol_agent.py render 4HHB.pdb hemoglobin.png --color chain
python3 scripts/pymol_agent.py spin 4HHB.pdb hemoglobin_spin.mp4 --frames 60 --color chain
```

**Shows:** the four globin subunits in distinct colors, and a 360° turntable
movie of the assembly.

## 4. Protein–DNA assembly — nucleosome core particle (1AOI)

```bash
python3 scripts/fetch_pdb.py 1AOI --format pdb
python3 scripts/pymol_agent.py render 1AOI.pdb nucleosome.png --color chain
```

**Shows:** the histone octamer wrapped by duplex DNA, each chain colored
separately.

## 5. AlphaFold confidence — human p53 (P04637)

```bash
python3 scripts/uniprot_lookup.py TP53 --gene-exact --json
python3 scripts/fetch_alphafold.py P04637 --pae --json
python3 scripts/pae_report.py AF-P04637-F1_pae.json --json
python3 scripts/pymol_agent.py render AF-P04637-F1.pdb p53_confidence.png --color plddt
```

**Shows:** the prediction colored by the official AlphaFold pLDDT scheme (blue
high-confidence core, orange low-confidence disordered regions), with the PAE
report flagging the inter-domain uncertainty.

## 6. Predicted vs experimental / conformational change — adenylate kinase

```bash
python3 scripts/fetch_pdb.py 4AKE --format pdb   # open conformation
python3 scripts/fetch_pdb.py 1AKE --format pdb   # closed conformation
python3 scripts/compare_structures.py 4AKE.pdb 1AKE.pdb --per-residue --json
```

**Shows:** the CE-alignment RMSD plus per-residue deviations that pinpoint the
LID and NMP domain motions between the open and closed states.

## 7. Cryo-EM density fit — sidechain in density

```bash
# Simulated density (no map download): does this sidechain sit in its density?
python3 scripts/pymol_agent.py density 1HSG.pdb asp25_density.png \
    --simulate --residue "chain A and resi 25"

# With a real map, the contour level is chosen from the map's sigma:
python3 scripts/map_info.py your_map.mrc --json
python3 scripts/pymol_agent.py density model.pdb fit.png --map your_map.mrc --residue "chain A and resi 25"
```

**Shows:** the selected residue as sticks carved inside its density mesh, with
the rest of the model as transparent cartoon for context.

## 8. Publication presets and illustration style

```bash
python3 scripts/pymol_agent.py render 4HHB.pdb hb_soft.png --preset soft --color chain
python3 scripts/pymol_agent.py render 1HSG.pdb hiv_illustration.png --preset illustration
```

**Shows:** the same structures in the neutral-background `soft` preset and the
outlined `illustration` (molecular-illustration) style.

## 9. End-to-end target dossier — KRAS G12C (the "epic" demo)

```bash
# Discover the structural landscape (200+ KRAS structures in the PDB)
python3 scripts/pdb_search.py --uniprot P01116 --rows 50 --details

# Build the full dossier: fetch ~9 structures, run every analysis,
# render the figure gallery + turntable, assemble a Markdown report
python3 scripts/kras_dossier.py --out kras_g12c_dossier
```

**Shows:** a single command that reproduces a day's worth of manual structural
review. It walks the KRAS G12C story end to end — the AlphaFold model (ordered
G-domain vs. disordered membrane-anchoring tail), the GTP/GDP conformational
switch (CE-align localizing the motion to Switch I/II), a four-structure
covalent-inhibitor timeline from chemical probe to approved drug (sotorasib,
adagrasib) in the cryptic switch-II pocket, and the SOS1/RAF protein interfaces
— then writes it all up as `DOSSIER.md` with embedded figures and a turntable
movie. A worked example of chaining the helper scripts into one analysis.

---

For the full command surface, see [`SKILL.md`](SKILL.md). For tool-specific
recipes and gotchas, see [`references/`](references/).
