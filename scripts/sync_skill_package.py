#!/usr/bin/env python3
"""Synchronize or verify the installable skills/proteus package copy."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "skills" / "proteus"
SOURCES = {
    ROOT / "SKILL.md": PACKAGE / "SKILL.md",
    ROOT / "agents" / "openai.yaml": PACKAGE / "agents" / "openai.yaml",
}
for directory in ("scripts", "references"):
    for source in sorted((ROOT / directory).iterdir()):
        if not source.is_file() or source.name == Path(__file__).name or source.suffix not in {".py", ".md"}:
            continue
        SOURCES[source] = PACKAGE / directory / source.name


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare() -> dict:
    missing = []
    changed = []
    for source, target in SOURCES.items():
        if not target.is_file():
            missing.append(str(target.relative_to(ROOT)))
        elif _digest(source) != _digest(target):
            changed.append(str(target.relative_to(ROOT)))
    expected = {target.resolve() for target in SOURCES.values()}
    extra = []
    for pattern in ("scripts/*.py", "references/*.md"):
        for target in PACKAGE.glob(pattern):
            if target.resolve() not in expected:
                extra.append(str(target.relative_to(ROOT)))
    return {"status": "ok" if not (missing or changed or extra) else "out_of_sync", "missing": missing, "changed": changed, "extra": extra}


def sync() -> dict:
    for source, target in SOURCES.items():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    expected = {target.resolve() for target in SOURCES.values()}
    removed = []
    for pattern in ("scripts/*.py", "references/*.md"):
        for target in PACKAGE.glob(pattern):
            if target.resolve() not in expected:
                removed.append(str(target.relative_to(ROOT)))
                target.unlink()
    result = compare()
    result["removed"] = removed
    result["copied"] = len(SOURCES)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify or synchronize the installable Proteus package copy.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="Fail if package files differ from canonical files")
    mode.add_argument("--sync", action="store_true", help="Copy canonical files into skills/proteus")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    report = sync() if args.sync else compare()
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
