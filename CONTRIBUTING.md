# Contributing

Proteus is most useful when it captures real structural-biology workflow
knowledge: command-line patterns, gotchas, reproducible helper scripts, and
small tests that keep agents from repeating common mistakes.

## Good Contributions

- a PyMOL, ChimeraX, AlphaFold, RCSB, UniProt, or Rosetta gotcha with a tested fix
- a small stdlib-only helper workflow that returns JSON
- focused docs that help agents choose the right tool
- tiny fixtures and tests for edge cases

## Development

```bash
make test
python3 scripts/proteus_doctor.py --json
```

Keep helper scripts dependency-free unless there is a clear reason not to.
Prefer deterministic JSON output in `--json` mode:

```json
{"status": "ok", "data": {}}
```

Errors should remain parseable:

```json
{"status": "error", "error": "message"}
```

## Pull Requests

- Keep changes scoped.
- Add or update tests for script behavior.
- Do not commit downloaded structures, private data, generated AlphaFold files,
  caches, or local machine artifacts.
- Document any external API or local binary assumptions.
