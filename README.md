# Proteus

Proteus is a structural biology skill for coding agents. It teaches the agent how to choose between PyMOL, ChimeraX, AlphaFold DB, and Rosetta-style workflows, and it bundles helper scripts that make those tools usable from terminal-first automation.

Named after the shape-shifting Greek god, and after the root of the word *protein*.

## What it provides

- Tool detection for PyMOL and ChimeraX on common macOS and Linux installs
- Headless PyMOL rendering for static structure figures
- ChimeraX analysis helpers for alignment, SASA, and hydrogen-bond workflows
- AlphaFold DB fetch and confidence interpretation
- On-demand references instead of one oversized `SKILL.md`
- Structured helper-script outputs suitable for agent use

## Install

For Codex:

```bash
git clone https://github.com/jvogan/proteus.git ~/.codex/skills/proteus
```

For Claude Code:

```bash
git clone https://github.com/jvogan/proteus.git ~/.claude/skills/proteus
```

You can also copy this directory into your local skills folder manually.

## Runtime requirements

At least one of the local structure tools is recommended:

| Tool | Role | Install |
|---|---|---|
| **PyMOL** | Headless rendering, structure inspection | [pymol.org](https://pymol.org) or `conda install -c conda-forge pymol-open-source` |
| **ChimeraX** | Analysis, GUI demos, cryo-EM work | [cgl.ucsf.edu/chimerax](https://www.cgl.ucsf.edu/chimerax/download.html) |
| **AlphaFold DB** | Public prediction source | No install required |
| **PyRosetta** | Optional scoring and design workflows | `pip install pyrosetta-installer` plus academic/commercial license |

Proteus degrades gracefully when those tools are absent. `scripts/pdb_info.py` and AlphaFold metadata fetches still work without local visualization software.

## Layout

```text
proteus/
├── SKILL.md
├── agents/openai.yaml
├── references/
│   ├── alphafold.md
│   ├── chimerax.md
│   ├── pymol.md
│   └── rosetta.md
└── scripts/
    ├── chimerax_agent.py
    ├── fetch_alphafold.py
    ├── pdb_info.py
    └── pymol_agent.py
```

## Quick examples

With the skill installed, prompts like these should trigger it:

- `Fetch the AlphaFold prediction for p53 and show which regions look disordered.`
- `Render the 1HSG binding pocket in PyMOL and save a clean PNG.`
- `Compare an AlphaFold model to an experimental structure and report RMSD.`
- `Analyze the hydrogen bonds at a protein-protein interface in ChimeraX.`

You can also run the helper scripts directly:

```bash
python3 scripts/pdb_info.py structure.pdb
python3 scripts/fetch_alphafold.py P04637 --pae --json
python3 scripts/pymol_agent.py render structure.pdb output.png
python3 scripts/chimerax_agent.py align reference.pdb mobile.pdb
```

## Design intent

The core tool split is deliberate:

- PyMOL is the default for headless image generation.
- ChimeraX is the default for analysis-heavy workflows and GUI-controlled sessions.
- AlphaFold DB is used as a public upstream source, not a local model runtime.
- Rosetta and PyRosetta are optional extensions, not hard requirements.

## Design notes

- Helper scripts emit machine-readable JSON by default, with human-readable text as a fallback.
- Temporary JSON handoff files are per-process, so parallel agent runs do not collide.
- Codex-facing metadata lives in `agents/openai.yaml`.

## License

MIT. See [LICENSE](LICENSE).
