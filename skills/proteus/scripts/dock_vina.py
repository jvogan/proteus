#!/usr/bin/env python3
"""Dry-run AutoDock Vina command planner and log parser.

This helper is intentionally stdlib-only. It detects common Vina-compatible
executables, validates local planning inputs, renders a command for review, and
parses Vina text logs. It never executes docking tools.

Usage:
    python3 scripts/dock_vina.py detect --json
    python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --config vina_box.txt --json
    python3 scripts/dock_vina.py parse-log vina.log --json
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ENGINE_ORDER = ("vina", "smina", "gnina")
PDBQT_SUFFIX = ".pdbqt"
CONFIG_FLOAT_KEYS = {
    "center_x", "center_y", "center_z",
    "size_x", "size_y", "size_z",
    "exhaustiveness", "energy_range",
    "num_modes", "seed", "cpu",
}
BOX_KEYS = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
NUMERIC_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$")
MODE_ROW_RE = re.compile(
    r"^\s*(?P<mode>\d+)\s+"
    r"(?P<affinity>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
    r"(?P<rmsd_lb>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
    r"(?P<rmsd_ub>[+-]?(?:\d+(?:\.\d*)?|\.\d+))(?:\s+.*)?$"
)
REMARK_RESULT_RE = re.compile(
    r"REMARK\s+VINA\s+RESULT:\s+"
    r"(?P<affinity>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
    r"(?P<rmsd_lb>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
    r"(?P<rmsd_ub>[+-]?(?:\d+(?:\.\d*)?|\.\d+))",
    re.IGNORECASE,
)


class DockVinaError(ValueError):
    """Raised when Vina planning or parsing cannot proceed."""


def _payload(status: str, data: dict[str, Any] | None = None,
             error: str | None = None) -> dict[str, Any]:
    output: dict[str, Any] = {"status": status}
    if data is not None:
        output["data"] = data
        output.update(data)
    if error is not None:
        output["error"] = _scrub_text(error)
    return output


def _display_path(path: str | Path | None) -> str | None:
    """Return a report-safe path label without exposing absolute local paths."""

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


def _scrub_text(text: str) -> str:
    """Remove Unix-style absolute paths from public output strings."""

    def replace(match: re.Match[str]) -> str:
        path = Path(match.group(0))
        name = path.name or "path"
        return f"{name} (absolute path omitted)"

    return re.sub(r"(?<![\w.-])/(?:[^\s'\",)]+/?)+", replace, text)


def _scrub_arg(value: str) -> str:
    if "=" in value:
        key, raw = value.split("=", 1)
        if raw and Path(raw).expanduser().is_absolute():
            return f"{key}={_display_path(raw)}"
        return _scrub_text(value)
    if Path(value).expanduser().is_absolute():
        return _display_path(value) or value
    return _scrub_text(value)


def _add_warning(warnings: list[str], message: str) -> None:
    message = _scrub_text(message)
    if message not in warnings:
        warnings.append(message)


def _convert_number(value: str) -> int | float | str:
    if not NUMERIC_RE.fullmatch(value):
        return value
    number = float(value)
    if number.is_integer() and not any(marker in value.lower() for marker in (".", "e")):
        return int(number)
    return number


def _detect_executable(name: str) -> dict[str, Any]:
    path = shutil.which(name)
    return {
        "available": bool(path),
        "executable": name,
        "path": _display_path(path) if path else None,
    }


def detect_tools() -> dict[str, Any]:
    executables = {name: _detect_executable(name) for name in ENGINE_ORDER}
    return _payload(
        "ok",
        {
            "executables": executables,
            "capabilities": {
                "vina": executables["vina"]["available"],
                "vina_compatible": any(item["available"] for item in executables.values()),
            },
            "notes": [
                "Detection uses shutil.which only.",
                "No Vina-compatible docking command was executed.",
            ],
        },
    )


def parse_vina_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).expanduser()
    if not config_path.exists():
        raise DockVinaError(f"Config file not found: {_display_path(path)}")
    if not config_path.is_file():
        raise DockVinaError(f"Config input is not a file: {_display_path(path)}")

    values: dict[str, Any] = {}
    unknown_lines: list[int] = []
    duplicate_keys: list[str] = []
    for line_number, raw_line in enumerate(config_path.read_text(errors="replace").splitlines(), start=1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            unknown_lines.append(line_number)
            continue
        key, raw_value = line.split("=", 1)
        key = key.strip()
        raw_value = raw_value.strip()
        if not key:
            unknown_lines.append(line_number)
            continue
        if key in values:
            duplicate_keys.append(key)
        values[key] = _convert_number(raw_value) if key in CONFIG_FLOAT_KEYS else raw_value

    missing_box_keys = sorted(BOX_KEYS - set(values))
    return {
        "path": _display_path(config_path.resolve(strict=False)),
        "exists": True,
        "values": values,
        "missing_box_keys": missing_box_keys,
        "duplicate_keys": sorted(set(duplicate_keys)),
        "unparsed_line_numbers": unknown_lines,
    }


def _pdbqt_summary(path: Path) -> dict[str, Any]:
    summary = {
        "atom_records": 0,
        "hetatm_records": 0,
        "model_records": 0,
        "remark_records": 0,
        "branch_records": 0,
        "torsdof_records": 0,
    }
    if not path.exists() or not path.is_file():
        return summary
    for line in path.read_text(errors="replace").splitlines():
        record = line[:6].strip().upper()
        if record == "ATOM":
            summary["atom_records"] += 1
        elif record == "HETATM":
            summary["hetatm_records"] += 1
        elif record == "MODEL":
            summary["model_records"] += 1
        elif record == "REMARK":
            summary["remark_records"] += 1
        elif record == "BRANCH":
            summary["branch_records"] += 1
        elif record == "TORSDOF":
            summary["torsdof_records"] += 1
    return summary


def _input_record(path_value: str | Path, role: str, warnings: list[str],
                  expected_suffix: str = PDBQT_SUFFIX) -> dict[str, Any]:
    path = Path(path_value).expanduser()
    display = _display_path(path_value) or str(path_value)
    record: dict[str, Any] = {
        "role": role,
        "path": display,
        "exists": path.exists(),
        "suffix": path.suffix.lower(),
    }

    if not path.exists():
        _add_warning(warnings, f"{role.title()} input was not found locally: {display}")
        record["kind"] = "missing_or_unverified"
        return record
    if not path.is_file():
        _add_warning(warnings, f"{role.title()} input is not a file: {display}")
        record["kind"] = "not_file"
        return record

    try:
        size_bytes = path.stat().st_size
    except OSError as exc:
        _add_warning(warnings, f"Could not stat {role} input: {exc}")
        size_bytes = None
    record["kind"] = role
    record["size_bytes"] = size_bytes
    if size_bytes == 0:
        _add_warning(warnings, f"{role.title()} input is empty: {display}")
    if path.suffix.lower() != expected_suffix:
        _add_warning(warnings, f"{role.title()} input has an uncommon suffix for Vina: {display}")
    record.update(_pdbqt_summary(path))
    if record["atom_records"] + record["hetatm_records"] == 0:
        _add_warning(warnings, f"{role.title()} PDBQT has no ATOM/HETATM records: {display}")
    return record


def _output_record(path_value: str | Path, role: str) -> dict[str, Any]:
    path = Path(path_value).expanduser()
    return {
        "role": role,
        "path": _display_path(path_value) or str(path_value),
        "exists": path.exists(),
        "suffix": path.suffix.lower(),
    }


def _normalize_vector(values: list[float] | tuple[float, float, float] | None,
                      name: str) -> list[float] | None:
    if values is None:
        return None
    if len(values) != 3:
        raise DockVinaError(f"--{name} requires exactly three numeric values.")
    return [float(item) for item in values]


def _select_engine(engine: str, detection: dict[str, Any]) -> str:
    if engine != "auto":
        return engine
    for name in ENGINE_ORDER:
        if detection["executables"][name]["available"]:
            return name
    return "vina"


def _planned_executable(engine: str, detection: dict[str, Any]) -> str:
    info = detection["executables"][engine]
    if info["available"] and info.get("executable"):
        return info["executable"]
    return engine


def _append_optional(command: list[str], flag: str, value: int | float | str | None) -> None:
    if value is not None:
        command.extend([flag, str(value)])


def plan_docking(
    receptor: str | Path,
    ligand: str | Path,
    config: str | Path | None = None,
    center: list[float] | tuple[float, float, float] | None = None,
    size: list[float] | tuple[float, float, float] | None = None,
    out: str | Path = "vina_out.pdbqt",
    log: str | Path = "vina.log",
    engine: str = "auto",
    exhaustiveness: int | None = None,
    num_modes: int | None = None,
    energy_range: float | None = None,
    seed: int | None = None,
    cpu: int | None = None,
    extra_args: list[str] | None = None,
    detection: dict[str, Any] | None = None,
) -> dict[str, Any]:
    warnings: list[str] = []
    center_values = _normalize_vector(center, "center")
    size_values = _normalize_vector(size, "size")
    if (center_values is None) != (size_values is None):
        raise DockVinaError("--center and --size must be supplied together.")

    detected = detection if detection is not None else detect_tools()
    selected_engine = _select_engine(engine, detected)
    if selected_engine not in ENGINE_ORDER:
        raise DockVinaError(f"Unknown Vina engine: {selected_engine}")
    if not detected["executables"][selected_engine]["available"]:
        _add_warning(warnings, f"{selected_engine} was not detected; command is a plan only.")

    receptor_record = _input_record(receptor, "receptor", warnings)
    ligand_record = _input_record(ligand, "ligand", warnings)
    out_record = _output_record(out, "output_pose")
    log_record = _output_record(log, "log")

    config_record = None
    if config:
        config_path = Path(config).expanduser()
        config_record = {
            "path": _display_path(config) or str(config),
            "exists": config_path.exists(),
        }
        if config_path.exists() and config_path.is_file():
            config_record.update(parse_vina_config(config))
            if config_record["missing_box_keys"]:
                _add_warning(
                    warnings,
                    "Config is missing Vina box keys: " + ", ".join(config_record["missing_box_keys"]),
                )
            if config_record["duplicate_keys"]:
                _add_warning(
                    warnings,
                    "Config contains duplicate keys: " + ", ".join(config_record["duplicate_keys"]),
                )
            if config_record["unparsed_line_numbers"]:
                _add_warning(
                    warnings,
                    "Config contains unparsed lines: "
                    + ", ".join(str(item) for item in config_record["unparsed_line_numbers"]),
                )
        else:
            _add_warning(warnings, f"Config file was not found locally: {config_record['path']}")
    elif center_values is None:
        _add_warning(warnings, "No config file or --center/--size box was provided.")

    if center_values is not None and size_values is not None:
        for axis, value in zip(("x", "y", "z"), size_values):
            if value <= 0:
                raise DockVinaError(f"Docking box size_{axis} must be greater than zero.")
    if exhaustiveness is not None and exhaustiveness <= 0:
        raise DockVinaError("--exhaustiveness must be greater than zero.")
    if num_modes is not None and num_modes <= 0:
        raise DockVinaError("--num-modes must be greater than zero.")
    if cpu is not None and cpu < 0:
        raise DockVinaError("--cpu must be greater than or equal to zero.")

    _add_warning(warnings, "No Vina-compatible docking command was executed.")
    _add_warning(
        warnings,
        "Review protonation, charges, grid placement, exhaustiveness, and pose validation before using planned results.",
    )

    command = [
        _planned_executable(selected_engine, detected),
        "--receptor", receptor_record["path"],
        "--ligand", ligand_record["path"],
    ]
    if config_record is not None:
        command.extend(["--config", config_record["path"]])
    if center_values is not None and size_values is not None:
        for axis, value in zip(("x", "y", "z"), center_values):
            command.extend([f"--center_{axis}", str(value)])
        for axis, value in zip(("x", "y", "z"), size_values):
            command.extend([f"--size_{axis}", str(value)])
    _append_optional(command, "--exhaustiveness", exhaustiveness)
    _append_optional(command, "--num_modes", num_modes)
    _append_optional(command, "--energy_range", energy_range)
    _append_optional(command, "--seed", seed)
    _append_optional(command, "--cpu", cpu)
    command.extend(["--out", out_record["path"], "--log", log_record["path"]])
    command.extend(_scrub_arg(item) for item in (extra_args or []))

    command_record = {
        "engine": selected_engine,
        "available": bool(detected["executables"][selected_engine]["available"]),
        "dry_run": True,
        "command": command,
        "command_line": shlex.join(command),
        "outputs": {
            "pose": out_record,
            "log": log_record,
        },
    }

    data = {
        "planner": "dock_vina",
        "dry_run": True,
        "execute": False,
        "requested_engine": engine,
        "selected_engine": selected_engine,
        "inputs": {
            "receptor": receptor_record,
            "ligand": ligand_record,
            "config": config_record,
        },
        "box": {
            "source": "inline" if center_values is not None else ("config" if config_record is not None else "missing"),
            "center": dict(zip(("x", "y", "z"), center_values)) if center_values is not None else None,
            "size": dict(zip(("x", "y", "z"), size_values)) if size_values is not None else None,
        },
        "commands": [command_record],
        "tools": detected["executables"],
        "capabilities": detected["capabilities"],
        "warnings": warnings,
    }
    return _payload("ok", data)


def parse_vina_log(path: str | Path, limit: int | None = None) -> dict[str, Any]:
    log_path = Path(path).expanduser()
    if not log_path.exists():
        raise DockVinaError(f"Vina log not found: {_display_path(path)}")
    if not log_path.is_file():
        raise DockVinaError(f"Vina log input is not a file: {_display_path(path)}")
    if limit is not None and limit <= 0:
        raise DockVinaError("--limit must be greater than zero.")

    metadata: dict[str, str] = {}
    modes: list[dict[str, Any]] = []
    remark_modes: list[dict[str, Any]] = []
    warnings: list[str] = []
    table_seen = False

    lines = log_path.read_text(errors="replace").splitlines()
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.lower().startswith("autodock vina"):
            metadata["version"] = _scrub_text(line)
            continue
        if ":" in line:
            key, value = line.split(":", 1)
            key_norm = re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")
            if key_norm in {
                "scoring_function", "rigid_receptor", "ligand", "output",
                "grid_center", "grid_size", "exhaustiveness", "cpu", "verbosity",
            }:
                metadata[key_norm] = _scrub_text(value.strip())
        if line.startswith("-----+") or ("affinity" in line.lower() and "mode" in line.lower()):
            table_seen = True
            continue

        match = MODE_ROW_RE.match(line)
        if match and table_seen:
            modes.append({
                "mode": int(match.group("mode")),
                "affinity_kcal_mol": float(match.group("affinity")),
                "rmsd_lb": float(match.group("rmsd_lb")),
                "rmsd_ub": float(match.group("rmsd_ub")),
                "line_number": line_number,
                "source": "table",
            })
            continue

        remark_match = REMARK_RESULT_RE.search(line)
        if remark_match:
            remark_modes.append({
                "mode": len(remark_modes) + 1,
                "affinity_kcal_mol": float(remark_match.group("affinity")),
                "rmsd_lb": float(remark_match.group("rmsd_lb")),
                "rmsd_ub": float(remark_match.group("rmsd_ub")),
                "line_number": line_number,
                "source": "remark",
            })

    source = "table"
    if not modes and remark_modes:
        modes = remark_modes
        source = "remark"
    if not modes:
        _add_warning(warnings, "No Vina pose table or REMARK VINA RESULT rows were parsed.")

    ranked = sorted(modes, key=lambda item: item["affinity_kcal_mol"])
    if limit is not None:
        ranked = ranked[:limit]
    ranks = []
    for rank, mode in enumerate(ranked, start=1):
        item = {"rank": rank}
        item.update(mode)
        ranks.append(item)

    data = {
        "file": _display_path(log_path.resolve(strict=False)),
        "parser": "vina_log",
        "source": source,
        "metadata": metadata,
        "mode_count": len(modes),
        "ranked_count": len(ranks),
        "best_mode": ranks[0] if ranks else None,
        "modes": modes,
        "ranks": ranks,
        "warnings": warnings,
    }
    return _payload("ok", data)


def _print_detect_text(report: dict[str, Any]) -> None:
    for name, tool in report["executables"].items():
        print(f"{name}: {tool['path'] or 'unavailable'}")


def _print_plan_text(report: dict[str, Any]) -> None:
    print(f"Dry run: {str(report['dry_run']).lower()}")
    print(f"Selected engine: {report['selected_engine']}")
    for warning in report["warnings"]:
        print(f"Warning: {warning}")
    for command in report["commands"]:
        print(f"Command: {command['command_line']}")


def _print_parse_text(report: dict[str, Any]) -> None:
    print(f"File: {report['file']}")
    print(f"Modes parsed: {report['mode_count']}")
    if report["best_mode"]:
        best = report["best_mode"]
        print(
            "Best mode: "
            f"{best['mode']} affinity {best['affinity_kcal_mol']:.3f} kcal/mol "
            f"rmsd_lb {best['rmsd_lb']:.3f} rmsd_ub {best['rmsd_ub']:.3f}"
        )
    for warning in report["warnings"]:
        print(f"Warning: {warning}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stdlib-only Vina-compatible docking planner and log parser.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/dock_vina.py detect --json\n"
            "  python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --config vina_box.txt --json\n"
            "  python3 scripts/dock_vina.py plan --receptor receptor.pdbqt --ligand ligand.pdbqt --center 0 0 0 --size 20 20 20 --json\n"
            "  python3 scripts/dock_vina.py parse-log vina.log --limit 3 --json"
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    detect_parser = subparsers.add_parser("detect", help="Detect Vina-compatible executables")
    detect_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    plan_parser = subparsers.add_parser("plan", help="Plan a Vina-compatible command without executing it")
    plan_parser.add_argument("--receptor", required=True, help="Prepared receptor PDBQT path")
    plan_parser.add_argument("--ligand", required=True, help="Prepared ligand PDBQT path")
    plan_parser.add_argument("--config", help="Vina config file with center/size fields")
    plan_parser.add_argument("--center", type=float, nargs=3, metavar=("X", "Y", "Z"), help="Inline docking box center")
    plan_parser.add_argument("--size", type=float, nargs=3, metavar=("X", "Y", "Z"), help="Inline docking box size")
    plan_parser.add_argument("--out", default="vina_out.pdbqt", help="Planned output pose path")
    plan_parser.add_argument("--log", default="vina.log", help="Planned Vina log path")
    plan_parser.add_argument("--engine", choices=["auto", *ENGINE_ORDER], default="auto", help="Vina-compatible engine to plan")
    plan_parser.add_argument("--exhaustiveness", type=int, help="Vina exhaustiveness")
    plan_parser.add_argument("--num-modes", type=int, help="Number of output poses")
    plan_parser.add_argument("--energy-range", type=float, help="Energy range in kcal/mol")
    plan_parser.add_argument("--seed", type=int, help="Random seed")
    plan_parser.add_argument("--cpu", type=int, help="CPU count; Vina uses 0 for auto")
    plan_parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Additional safe command token to append to the plan; repeat as needed",
    )
    plan_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    parse_parser = subparsers.add_parser("parse-log", help="Parse Vina log table into ranked JSON")
    parse_parser.add_argument("logfile", help="Vina text log path")
    parse_parser.add_argument("--limit", type=int, help="Maximum number of ranked modes to emit")
    parse_parser.add_argument("--json", action="store_true", help="Emit JSON (default for this command)")

    args = parser.parse_args(argv)
    try:
        if args.command == "detect":
            report = detect_tools()
            if args.json:
                print(json.dumps(report, indent=2))
            else:
                _print_detect_text(report)
            return 0

        if args.command == "plan":
            report = plan_docking(
                args.receptor,
                args.ligand,
                config=args.config,
                center=args.center,
                size=args.size,
                out=args.out,
                log=args.log,
                engine=args.engine,
                exhaustiveness=args.exhaustiveness,
                num_modes=args.num_modes,
                energy_range=args.energy_range,
                seed=args.seed,
                cpu=args.cpu,
                extra_args=args.extra_arg,
            )
            if args.json:
                print(json.dumps(report, indent=2))
            else:
                _print_plan_text(report)
            return 0

        if args.command == "parse-log":
            report = parse_vina_log(args.logfile, limit=args.limit)
            if args.json:
                print(json.dumps(report, indent=2))
            else:
                _print_parse_text(report)
            return 0
    except DockVinaError as exc:
        error = _payload("error", error=str(exc))
        print(json.dumps(error, indent=2))
        return 1

    parser.error(f"Unhandled command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
