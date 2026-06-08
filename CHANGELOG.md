# Changelog

## Unreleased

- Added a `density` subcommand to `pymol_agent.py`: render a model in cryo-EM
  density (real `--map` or `--simulate`d gaussian density), with the mesh carved
  around the model or a `--residue` selection and the contour level taken from
  the map's sigma via `map_info.py`.
- Added `SHOWCASE.md`: copy-pasteable, reproducible figures and analyses on
  canonical public structures (every command is a real, tested code path).
- Added a `make release-check` target: pre-publish hygiene sweep for tracked
  structures/maps/media/secret files and large files.
- Added `interface_report.py`: zero-dependency protein-protein interface residue
  analysis between chains (the chain-chain analog of `pocket_report.py`).
- Added a `pocket` subcommand to `pymol_agent.py` for one-command annotated
  binding-pocket figures (ligand + pocket sticks + polar contacts + transparent
  context), and a non-empty-PNG guard so renders fail loudly instead of blank.
- Improved turntable framing: orthoscopic, bounding-sphere fit, and a widened
  depth slab so the structure never clips mid-rotation.
- Added `chimerax_rest.py`: a managed ChimeraX REST renderer that launches a GUI
  session on an ephemeral port, renders via GPU, defeats the 0-byte-PNG save
  race, and guarantees teardown.
- Added headless turntable movies (`pymol_agent.py spin` and `chimerax_rest.py
  spin`): ray-traced frames encoded with ffmpeg, degrading gracefully when
  ffmpeg is absent.
- Added render presets (`--preset publication|illustration|soft`) and pLDDT
  confidence coloring (`--color plddt`) to `pymol_agent.py render`.
- Added `add_helix_records.py`: inject HELIX records into CA-only backbones so
  cartoons render for de-novo designs.
- Added `map_info.py`: MRC/CCP4 map inspection with sigma-based contour levels.
- Documented six new gotchas (rendering, animation, and cryo-EM maps) and added
  turntable, managed-REST, density-map, and exploded-comparison recipes.
- Added ffmpeg detection to `proteus_doctor.py` and hardened `.gitignore`
  against accidentally committing structures, maps, and rendered media.

## v0.1.0 - Public Launch

- Added the core Proteus structural biology agent skill.
- Added stdlib-only helpers for RCSB PDB, UniProt, AlphaFold DB, PDB/mmCIF
  inspection, PAE reports, validation summaries, pocket reports, and structure
  comparison.
- Added PyMOL and ChimeraX automation wrappers with safer path handling,
  structured JSON output, and documented gotchas.
- Added public launch documentation, privacy/security/disclaimer files, curated
  snapshots, and GitHub CI/security workflows.
