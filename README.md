# Proteus

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

Proteus is a **skill** for AI coding agents — a drop-in knowledge pack that teaches [Claude Code](https://docs.anthropic.com/en/docs/claude-code), [Codex](https://openai.com/index/codex/), and similar agents how to drive structural biology tools from the terminal.

It covers PyMOL, ChimeraX, AlphaFold DB, and Rosetta/PyRosetta: when to use which tool, how to call them headlessly, and the undocumented gotchas that otherwise cost hours of debugging.

Named after the shape-shifting Greek god — and after the root of the word *protein*.

> **What is a skill?** A skill is a directory that an AI coding agent reads to gain domain-specific knowledge. Clone it into the agent's skills folder and it becomes part of the agent's working context — no code changes or plugin installs required.

## What it provides

- **17 documented gotchas** for PyMOL, ChimeraX, and AlphaFold DB — hard-won from real debugging
- Tool detection for PyMOL and ChimeraX across macOS and Linux installs
- Headless PyMOL rendering for publication-quality structure figures
- ChimeraX analysis helpers for alignment, SASA, and hydrogen-bond workflows
- AlphaFold DB fetch with confidence interpretation and pLDDT coloring
- Rosetta/PyRosetta patterns plus ML alternatives (ProteinMPNN, ESM2)
- Zero-dependency PDB file inspector (`pdb_info.py` — stdlib only)
- Structured JSON output from all helper scripts, safe for parallel agent runs

## Install

Clone into your agent's skills directory:

```bash
# Claude Code
git clone https://github.com/jvogan/proteus.git ~/.claude/skills/proteus

# Codex
git clone https://github.com/jvogan/proteus.git ~/.codex/skills/proteus
```

Or copy the directory manually into your agent's skills folder. The skill path may vary by agent version — check your agent's documentation if the above doesn't work.

## Runtime requirements

Proteus degrades gracefully — `pdb_info.py` and AlphaFold metadata fetches work with zero local tools. For full capability, install at least one:

| Tool | Role | Install |
|---|---|---|
| **PyMOL** | Headless rendering, structure inspection | [pymol.org](https://pymol.org) or `conda install -c conda-forge pymol-open-source` |
| **ChimeraX** | Analysis, GUI demos, cryo-EM visualization | [cgl.ucsf.edu/chimerax](https://www.cgl.ucsf.edu/chimerax/download.html) |
| **AlphaFold DB** | Public prediction database | No install — uses the [EBI REST API](https://alphafold.ebi.ac.uk) |
| **PyRosetta** | Scoring, energy minimization, protein design | `pip install pyrosetta-installer` (academic license required) |

Python 3.8+ is required. All helper scripts use only the standard library.

## Quick examples

With the skill installed, natural-language prompts trigger it automatically:

```
Fetch the AlphaFold prediction for p53 and show which regions look disordered.
Render the 1HSG binding pocket in PyMOL and save a clean PNG.
Compare an AlphaFold model to an experimental structure and report RMSD.
Analyze the hydrogen bonds at a protein-protein interface in ChimeraX.
```

The helper scripts also work standalone:

```bash
python3 scripts/pdb_info.py structure.pdb                          # zero-dep PDB inspection
python3 scripts/fetch_alphafold.py P04637 --pae --json             # AlphaFold fetch
python3 scripts/pymol_agent.py render structure.pdb output.png     # headless render
python3 scripts/chimerax_agent.py align reference.pdb mobile.pdb   # structure alignment
```

## Layout

```text
proteus/
├── SKILL.md              # Main skill — agent reads this first
├── agents/openai.yaml    # Codex discovery metadata
├── references/           # On-demand deep docs (loaded as needed)
│   ├── alphafold.md
│   ├── chimerax.md
│   ├── pymol.md
│   └── rosetta.md
└── scripts/              # Agent helper scripts (all stdlib-only)
    ├── chimerax_agent.py
    ├── fetch_alphafold.py
    ├── pdb_info.py
    └── pymol_agent.py
```

## Design intent

The tool split is deliberate:

- **PyMOL** is the default for headless image generation (software ray tracer — no display needed).
- **ChimeraX** is the default for analysis-heavy workflows and GPU-rendered GUI sessions.
- **AlphaFold DB** is a public upstream source, not a local model runtime.
- **Rosetta/PyRosetta** are optional extensions. ML alternatives (ProteinMPNN, ESM2) are documented for when Rosetta isn't available.

Helper scripts emit machine-readable JSON, with human-readable text as a fallback. Temporary handoff files are per-process, so parallel agent runs never collide.

## Contributing

Found a gotcha that isn't documented? Have a workflow that should be covered? [Open an issue](https://github.com/jvogan/proteus/issues) or submit a PR. The most valuable contributions are real debugging discoveries — the kind of thing that takes hours to figure out and one sentence to explain.

## License

MIT. See [LICENSE](LICENSE).
