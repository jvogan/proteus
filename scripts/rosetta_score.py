#!/usr/bin/env python3
"""Optional Rosetta/PyRosetta scoring discovery, scorefile parsing, and planning.

This adapter is intentionally stdlib-only and dry by default. Detection checks
for executables and PyRosetta import metadata without invoking Rosetta scoring.

Usage:
    python3 scripts/rosetta_score.py detect --json
    python3 scripts/rosetta_score.py parse-scorefile score.sc --json
    python3 scripts/rosetta_score.py plan model.pdb --json
"""

import argparse
import importlib.util
import json
import os
import re
import shlex
import shutil
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
NUMERIC_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$")
INTEGER_RE = re.compile(r"^[+-]?\d+$")
STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif"}
ROSETTA_CANDIDATES = {
    "rosetta_scripts": [
        "rosetta_scripts",
        "rosetta_scripts.default.linuxgccrelease",
        "rosetta_scripts.static.linuxgccrelease",
        "rosetta_scripts.default.macosclangrelease",
        "rosetta_scripts.static.macosclangrelease",
    ],
    "score_jd2": [
        "score_jd2",
        "score_jd2.default.linuxgccrelease",
        "score_jd2.static.linuxgccrelease",
        "score_jd2.default.macosclangrelease",
        "score_jd2.static.macosclangrelease",
    ],
}


class RosettaScoreError(ValueError):
    """Raised when Rosetta score planning or parsing cannot proceed."""


def _payload(status: str, data: dict | None = None, error: str | None = None) -> dict:
    output = {"status": status}
    if data is not None:
        output["data"] = data
        output.update(data)
    if error is not None:
        output["error"] = error
    return output


def _display_path(path: str | Path | None) -> str | None:
    """Return a report-safe path label without exposing absolute home paths."""

    if path is None:
        return None
    text = str(path)
    if "://" in text:
        return text
    candidate = Path(text).expanduser()
    if not candidate.is_absolute():
        return text
    try:
        resolved = candidate.resolve(strict=False)
    except OSError:
        resolved = candidate.absolute()
    for root, prefix in ((ROOT, "."), (Path.cwd(), ".")):
        try:
            relative = resolved.relative_to(root.resolve())
        except ValueError:
            continue
        return f"{prefix}/{relative.as_posix()}"
    return f"{resolved.name} (absolute path omitted)"


def _scrub_arg(value: str) -> str:
    if "=" in value:
        key, raw = value.split("=", 1)
        if raw and Path(raw).expanduser().is_absolute():
            return f"{key}={_display_path(raw)}"
        return value
    if Path(value).expanduser().is_absolute():
        return _display_path(value) or value
    return value


def _add_warning(warnings: list[str], message: str) -> None:
    if message not in warnings:
        warnings.append(message)


def _is_executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def _rosetta_env_dirs() -> list[tuple[str, Path]]:
    dirs: list[tuple[str, Path]] = []
    rosetta_bin = os.environ.get("ROSETTA_BIN")
    if rosetta_bin:
        dirs.append(("ROSETTA_BIN", Path(rosetta_bin).expanduser()))
    for name in ("ROSETTA3", "ROSETTA_HOME", "ROSETTA"):
        raw = os.environ.get(name)
        if not raw:
            continue
        base = Path(raw).expanduser()
        dirs.extend(
            [
                (name, base / "source" / "bin"),
                (name, base / "main" / "source" / "bin"),
                (name, base / "bin"),
            ]
        )
    return dirs


def _find_executable(name: str, candidates: list[str]) -> dict:
    for executable in candidates:
        path = shutil.which(executable)
        if path:
            return {
                "available": True,
                "path": _display_path(path),
                "executable": executable,
                "source": "PATH",
                "candidates": candidates,
            }

    for env_name, directory in _rosetta_env_dirs():
        for executable in candidates:
            path = directory / executable
            if _is_executable(path):
                return {
                    "available": True,
                    "path": _display_path(path),
                    "executable": executable,
                    "source": env_name,
                    "candidates": candidates,
                }

    return {
        "available": False,
        "path": None,
        "executable": None,
        "source": None,
        "candidates": candidates,
        "message": f"{name} was not found on PATH or common Rosetta environment paths.",
    }


def _detect_pyrosetta() -> dict:
    try:
        spec = importlib.util.find_spec("pyrosetta")
    except Exception as exc:  # pragma: no cover - defensive around custom import hooks
        return {
            "available": False,
            "importable": False,
            "origin": None,
            "checked_by": "importlib.util.find_spec",
            "error": f"find_spec failed with {exc.__class__.__name__}",
        }
    if spec is None:
        return {
            "available": False,
            "importable": False,
            "origin": None,
            "checked_by": "importlib.util.find_spec",
        }
    origin = spec.origin
    if origin in {None, "built-in", "namespace"}:
        origin_display = origin
    else:
        origin_display = _display_path(origin)
    return {
        "available": True,
        "importable": True,
        "origin": origin_display,
        "checked_by": "importlib.util.find_spec",
    }


def detect_tools() -> dict:
    tools = {
        name: _find_executable(name, candidates)
        for name, candidates in ROSETTA_CANDIDATES.items()
    }
    pyrosetta = _detect_pyrosetta()
    return _payload(
        "ok",
        {
            "tools": tools,
            "pyrosetta": pyrosetta,
            "capabilities": {
                "rosetta_scripts": tools["rosetta_scripts"]["available"],
                "score_jd2": tools["score_jd2"]["available"],
                "pyrosetta": pyrosetta["available"],
            },
            "notes": [
                "Detection does not run Rosetta or PyRosetta scoring.",
                "PyRosetta availability is checked with importlib.util.find_spec, not import.",
            ],
        },
    )


def _convert_token(value: str) -> int | float | str:
    if INTEGER_RE.fullmatch(value):
        return int(value)
    if NUMERIC_RE.fullmatch(value):
        return float(value)
    return value


def _is_score_header(fields: list[str]) -> bool:
    return bool(fields) and not NUMERIC_RE.fullmatch(fields[0])


def parse_scorefile(
    path: str | Path,
    score_column: str = "total_score",
    limit: int | None = None,
    descending: bool = False,
) -> dict:
    score_path = Path(path).expanduser()
    if not score_path.exists():
        raise RosettaScoreError(f"Scorefile not found: {_display_path(path)}")
    if limit is not None and limit <= 0:
        raise RosettaScoreError("--limit must be greater than zero.")

    header: list[str] | None = None
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    headers_seen = 0

    for line_number, raw_line in enumerate(score_path.read_text(errors="replace").splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#") or line.startswith("SEQUENCE:"):
            continue
        if not line.startswith("SCORE:"):
            continue

        fields = line.split()[1:]
        if not fields:
            continue
        if _is_score_header(fields):
            header = fields
            headers_seen += 1
            continue
        if header is None:
            _add_warning(warnings, f"Skipped SCORE row before a header at line {line_number}.")
            continue

        if len(fields) > len(header) and header[-1] == "description":
            fields = fields[: len(header) - 1] + [" ".join(fields[len(header) - 1 :])]
        if len(fields) != len(header):
            _add_warning(
                warnings,
                f"Skipped malformed SCORE row at line {line_number}: expected {len(header)} columns, got {len(fields)}.",
            )
            continue
        rows.append({column: _convert_token(value) for column, value in zip(header, fields)})

    if header is None:
        raise RosettaScoreError("No SCORE header was found.")
    if score_column not in header:
        raise RosettaScoreError(f"Score column '{score_column}' was not found in SCORE header.")

    scorable: list[dict[str, Any]] = []
    for row in rows:
        value = row.get(score_column)
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            scorable.append(row)
        else:
            description = row.get("description", "<unknown>")
            _add_warning(
                warnings,
                f"Skipped row {description!r}; score column '{score_column}' was not numeric.",
            )

    ranked_rows = sorted(scorable, key=lambda item: item[score_column], reverse=descending)
    if limit is not None:
        ranked_rows = ranked_rows[:limit]
    ranks = []
    for index, row in enumerate(ranked_rows, start=1):
        ranked = {"rank": index}
        ranked.update(row)
        ranks.append(ranked)

    return _payload(
        "ok",
        {
            "file": _display_path(score_path.resolve(strict=False)),
            "score_column": score_column,
            "sort": "descending" if descending else "ascending",
            "columns": header,
            "headers_seen": headers_seen,
            "row_count": len(rows),
            "ranked_count": len(ranks),
            "ranks": ranks,
            "warnings": warnings,
        },
    )


def _structure_record(value: str, warnings: list[str]) -> dict:
    path = Path(value).expanduser()
    display = _display_path(value) or value
    record = {
        "input": display,
        "exists": path.exists(),
        "suffix": path.suffix.lower(),
    }
    if path.exists():
        record["kind"] = "local_structure"
    else:
        record["kind"] = "missing_or_unverified"
        _add_warning(warnings, f"Structure input was not found locally: {display}")
    if path.suffix.lower() and path.suffix.lower() not in STRUCTURE_SUFFIXES:
        _add_warning(
            warnings,
            f"Structure input has an uncommon suffix for Rosetta scoring: {display}",
        )
    if path.is_absolute() and "absolute path omitted" in display:
        _add_warning(warnings, f"Absolute path was scrubbed for structure input: {display}")
    return record


def _planned_executable(backend: str, detection: dict) -> str:
    if backend == "pyrosetta":
        return "python3"
    info = detection["tools"][backend]
    if info.get("available") and info.get("source") == "PATH" and info.get("executable"):
        return info["executable"]
    if info.get("available") and info.get("path"):
        return info["path"]
    return backend


def _select_backend(requested: str, protocol: str | None, detection: dict) -> str:
    if requested != "auto":
        return requested
    if protocol:
        return "rosetta_scripts"
    if detection["tools"]["score_jd2"]["available"]:
        return "score_jd2"
    if detection["pyrosetta"]["available"]:
        return "pyrosetta"
    return "score_jd2"


def _build_pyrosetta_command(
    structures: list[str],
    weights: str | None,
    extra_flags: list[str],
) -> list[str]:
    init_flags = ["-mute", "all"]
    if weights:
        init_flags[0:0] = ["-score:weights", weights]
    init_flags.extend(extra_flags)
    code = (
        "import sys, pyrosetta; "
        f"pyrosetta.init({' '.join(init_flags)!r}); "
        "scorefxn=pyrosetta.get_fa_scorefxn(); "
        "[print(path, scorefxn(pyrosetta.pose_from_file(path))) for path in sys.argv[1:]]"
    )
    return ["python3", "-c", code, *structures]


def plan_scoring(
    structures: list[str],
    backend: str = "auto",
    protocol: str | None = None,
    scorefile: str = "score.sc",
    weights: str | None = "ref2015",
    extra_flags: list[str] | None = None,
) -> dict:
    if not structures:
        raise RosettaScoreError("At least one structure path is required.")

    warnings: list[str] = []
    detection = detect_tools()
    selected = _select_backend(backend, protocol, detection)
    structure_records = [_structure_record(value, warnings) for value in structures]
    structure_args = [record["input"] for record in structure_records]
    extra_display = [_scrub_arg(item) for item in (extra_flags or [])]
    scorefile_display = _display_path(scorefile) or scorefile

    protocol_record = None
    protocol_display = None
    if protocol:
        protocol_path = Path(protocol).expanduser()
        protocol_display = _display_path(protocol) or protocol
        protocol_record = {"path": protocol_display, "exists": protocol_path.exists()}
        if not protocol_path.exists():
            _add_warning(warnings, f"RosettaScripts protocol XML was not found locally: {protocol_display}")
        if protocol_path.is_absolute() and "absolute path omitted" in protocol_display:
            _add_warning(warnings, f"Absolute path was scrubbed for protocol XML: {protocol_display}")

    if selected == "rosetta_scripts" and not protocol_display:
        raise RosettaScoreError("--protocol is required when planning rosetta_scripts commands.")

    if selected in {"score_jd2", "rosetta_scripts"}:
        tool_info = detection["tools"][selected]
        if not tool_info["available"]:
            _add_warning(warnings, f"{selected} was not detected; command is a plan only.")
    elif selected == "pyrosetta":
        tool_info = detection["pyrosetta"]
        if not tool_info["available"]:
            _add_warning(warnings, "PyRosetta was not importable; command is a plan only.")
        _add_warning(
            warnings,
            "PyRosetta command is a minimal score-only one-liner; use a reviewed script for production scoring.",
        )
    else:
        raise RosettaScoreError(f"Unknown backend: {selected}")

    _add_warning(warnings, "No Rosetta or PyRosetta scoring command was executed.")
    _add_warning(
        warnings,
        "Plan only scores existing coordinates; review params, ligands, membranes, symmetry, constraints, and relax/repack needs before running.",
    )

    executable = _planned_executable(selected, detection)
    if selected == "score_jd2":
        command = [executable, "-in:file:s", *structure_args]
        if weights:
            command.extend(["-score:weights", weights])
        command.extend(["-out:file:scorefile", scorefile_display, *extra_display])
    elif selected == "rosetta_scripts":
        command = [
            executable,
            "-parser:protocol",
            protocol_display or "",
            "-in:file:s",
            *structure_args,
        ]
        if weights:
            command.extend(["-score:weights", weights])
        command.extend(["-out:file:scorefile", scorefile_display, *extra_display])
    else:
        command = _build_pyrosetta_command(structure_args, weights, extra_display)

    command_record = {
        "backend": selected,
        "available": bool(tool_info["available"]),
        "dry_run": True,
        "command": command,
        "command_line": shlex.join(command),
        "scorefile": scorefile_display,
    }
    if protocol_record is not None:
        command_record["protocol"] = protocol_record

    return _payload(
        "ok",
        {
            "dry_run": True,
            "execute": False,
            "requested_backend": backend,
            "selected_backend": selected,
            "structures": structure_records,
            "commands": [command_record],
            "provenance": {
                "script": "scripts/rosetta_score.py",
                "cwd": _display_path(Path.cwd().resolve(strict=False)),
                "detection": {
                    "tools": detection["tools"],
                    "pyrosetta": detection["pyrosetta"],
                    "capabilities": detection["capabilities"],
                },
            },
            "warnings": warnings,
        },
    )


def _print_detect_text(report: dict) -> None:
    for name, tool in report["tools"].items():
        print(f"{name}: {tool['path'] or 'unavailable'}")
    pyrosetta = report["pyrosetta"]
    print(f"pyrosetta: {pyrosetta.get('origin') or 'unavailable'}")


def _print_plan_text(report: dict) -> None:
    print(f"Dry run: {str(report['dry_run']).lower()}")
    print(f"Selected backend: {report['selected_backend']}")
    for warning in report["warnings"]:
        print(f"Warning: {warning}")
    for command in report["commands"]:
        print(f"Command: {command['command_line']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Optional stdlib-only Rosetta/PyRosetta scoring adapter.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/rosetta_score.py detect --json\n"
            "  python3 scripts/rosetta_score.py parse-scorefile score.sc --limit 10 --json\n"
            "  python3 scripts/rosetta_score.py plan model.pdb --scorefile planned.sc --json"
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    detect_parser = subparsers.add_parser("detect", help="Detect optional Rosetta/PyRosetta availability")
    detect_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    parse_parser = subparsers.add_parser("parse-scorefile", help="Parse Rosetta score.sc into JSON ranks")
    parse_parser.add_argument("scorefile", help="Rosetta score.sc file")
    parse_parser.add_argument("--score-column", default="total_score", help="Numeric score column used for ranking")
    parse_parser.add_argument("--limit", type=int, help="Maximum number of ranked rows to emit")
    parse_parser.add_argument("--descending", action="store_true", help="Rank larger scores first")
    parse_parser.add_argument("--json", action="store_true", help="Emit JSON (default for this command)")

    plan_parser = subparsers.add_parser("plan", help="Plan local Rosetta/PyRosetta scoring commands without executing them")
    plan_parser.add_argument("structures", nargs="+", help="Local PDB/mmCIF structure path(s)")
    plan_parser.add_argument(
        "--backend",
        choices=["auto", "score_jd2", "rosetta_scripts", "pyrosetta"],
        default="auto",
        help="Scoring backend to plan",
    )
    plan_parser.add_argument("--protocol", help="RosettaScripts XML protocol path")
    plan_parser.add_argument("--scorefile", default="score.sc", help="Planned Rosetta scorefile output")
    plan_parser.add_argument("--weights", default="ref2015", help="Rosetta score weights name; use empty string to omit")
    plan_parser.add_argument(
        "--extra-flag",
        action="append",
        default=[],
        help="Additional Rosetta option token to append to the planned command; repeat as needed",
    )
    plan_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    args = parser.parse_args(argv)
    try:
        if args.command == "detect":
            report = detect_tools()
            if args.json:
                print(json.dumps(report, indent=2))
            else:
                _print_detect_text(report)
            return 0

        if args.command == "parse-scorefile":
            report = parse_scorefile(
                args.scorefile,
                score_column=args.score_column,
                limit=args.limit,
                descending=args.descending,
            )
            print(json.dumps(report, indent=2))
            return 0

        if args.command == "plan":
            weights = args.weights or None
            report = plan_scoring(
                args.structures,
                backend=args.backend,
                protocol=args.protocol,
                scorefile=args.scorefile,
                weights=weights,
                extra_flags=args.extra_flag,
            )
            if args.json:
                print(json.dumps(report, indent=2))
            else:
                _print_plan_text(report)
            return 0
    except RosettaScoreError as exc:
        error = _payload("error", error=str(exc))
        print(json.dumps(error, indent=2))
        return 1

    parser.error(f"Unhandled command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
