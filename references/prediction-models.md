# Prediction Models - Proteus Skill

Use this reference when the task goes beyond fetching an AlphaFold DB model,
especially complexes, ligands, nucleic acids, or new sequence designs.

## Decision Map

| Need | Suggested route | Notes |
|---|---|---|
| Existing UniProt single-chain prediction | AlphaFold DB | Fastest; no local model runtime |
| Protein-ligand or protein-DNA/RNA complex prediction | AlphaFold 3 / Boltz / Chai | Heavy runtimes; check install/GPU first |
| Many sequence variants on one backbone | ProteinMPNN + structure validation | Useful without Rosetta |
| Fast rough single-chain prediction | ESMFold / ColabFold | Trade speed, accuracy, and setup complexity |
| Scoring/minimization/design protocols | Rosetta/PyRosetta | License and install may be the blocker |

## AlphaFold 3

AlphaFold 3 is relevant for proteins, nucleic acids, ligands, ions, and modified
entities. Its input is JSON-oriented and can specify ligands by CCD code, SMILES,
or user-provided CCD.

Do not assume AF3 is installed. First check for a local repo, environment, model
parameters, and databases. If the user only needs an existing single-chain model,
use AlphaFold DB instead.

## Boltz / Chai

Boltz and Chai-style models are useful open biomolecular structure prediction
options for complex inputs. Treat them as optional local runners:

1. Detect whether the CLI/package is installed.
2. Confirm GPU/runtime requirements before launching.
3. Keep generated inputs and outputs in a dedicated working directory.
4. Post-process outputs with `structure_info.py`, PyMOL, or ChimeraX.

## Output Handling

Always report:

- input sequences/entities and ligand definitions
- model/runtime used
- output coordinate path
- confidence metrics available from that model
- caveats, especially for low-confidence interfaces and flexible regions
