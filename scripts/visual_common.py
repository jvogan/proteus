#!/usr/bin/env python3
"""Shared helpers for reproducible PyMOL and ChimeraX workflows."""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common


IDENTIFIER_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")
COLOR_RE = re.compile(r"^(?:[A-Za-z][A-Za-z0-9_]{0,31}|#[0-9A-Fa-f]{6})$")
STYLE_NAMES = {"cartoon", "sticks", "spheres", "surface", "lines", "ribbon"}


class VisualWorkflowError(RuntimeError):
    pass


def identifier(value: str, *, field: str = "identifier") -> str:
    if not IDENTIFIER_RE.fullmatch(value):
        raise VisualWorkflowError(f"Invalid {field}: {value!r}")
    return value


def color(value: str) -> str:
    if not COLOR_RE.fullmatch(value):
        raise VisualWorkflowError(f"Invalid color: {value!r}")
    return value


def selection(value: str, *, field: str = "selection") -> str:
    if not value or any(character in value for character in "\r\n;"):
        raise VisualWorkflowError(f"Invalid {field}: selections cannot be empty or contain newlines/semicolons.")
    return value


def style(value: str) -> str:
    if value not in STYLE_NAMES:
        raise VisualWorkflowError(f"Unsupported representation style: {value}")
    return value


def quote_pymol(value: str | Path) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def quote_chimerax(value: str | Path) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def pymol_png(path: str | Path, *, width: int | None = None, height: int | None = None) -> list[str]:
    """Write a PNG through the PyMOL Python API so quoted paths work reliably."""
    options = []
    if width is not None:
        options.append(f"width={int(width)}")
    if height is not None:
        options.append(f"height={int(height)}")
    suffix = ", " + ", ".join(options) if options else ""
    return ["python", f"cmd.png({str(Path(path).resolve())!r}{suffix})", "python end"]


def pymol_save(path: str | Path) -> list[str]:
    """Save through the PyMOL Python API so quoted/session paths stay portable."""
    return ["python", f"cmd.save({str(Path(path).resolve())!r})", "python end"]


def find_pymol() -> str | None:
    configured = os.environ.get("PYMOL_BIN", "")
    return proteus_common.find_executable(
        configured,
        "pymol",
        "/Applications/PyMOL.app/Contents/bin/pymol",
        "~/Applications/PyMOL.app/Contents/bin/pymol",
    )


def find_chimerax() -> str | None:
    configured = os.environ.get("CHIMERAX_BIN", "")
    candidates = [
        configured,
        "ChimeraX",
        "chimerax",
        "/Applications/ChimeraX.app/Contents/bin/ChimeraX",
    ]
    applications = Path("/Applications")
    if applications.is_dir():
        candidates.extend(str(item / "Contents/bin/ChimeraX") for item in sorted(applications.glob("ChimeraX*.app"), reverse=True))
    return proteus_common.find_executable(*candidates)


def load_records(path: str | Path) -> list[dict[str, Any]]:
    source = Path(path)
    if not source.is_file():
        raise VisualWorkflowError(f"File not found: {proteus_common.display_path(source)}")
    if source.suffix.lower() == ".csv":
        with source.open(newline="", encoding="utf-8-sig") as handle:
            return [dict(row) for row in csv.DictReader(handle)]
    payload = json.loads(source.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        return [dict(item) for item in payload]
    if isinstance(payload, dict):
        for key in ("records", "annotations", "restraints", "scenes"):
            value = payload.get(key)
            if isinstance(value, list):
                return [dict(item) for item in value]
    raise VisualWorkflowError("Expected a CSV or JSON list of records.")


def write_workflow(
    outdir: str | Path,
    name: str,
    *,
    report: dict[str, Any],
    pymol_lines: list[str] | None = None,
    chimerax_lines: list[str] | None = None,
) -> dict[str, str]:
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    if pymol_lines is not None:
        pml = destination / f"{name}.pml"
        pml.write_text("\n".join(pymol_lines).rstrip() + "\n", encoding="utf-8")
        outputs["pymol_script"] = proteus_common.display_path(pml)
    if chimerax_lines is not None:
        cxc = destination / f"{name}.cxc"
        cxc.write_text("\n".join(chimerax_lines).rstrip() + "\n", encoding="utf-8")
        outputs["chimerax_script"] = proteus_common.display_path(cxc)
    report_path = destination / f"{name}.json"
    outputs["report"] = proteus_common.display_path(report_path)
    data = report.get("data")
    if isinstance(data, dict):
        data["artifacts"] = dict(outputs)
    proteus_common.write_json(report_path, report)
    return outputs


def run_pymol(script_path: str | Path, *, timeout: int = 300) -> dict[str, Any]:
    executable = find_pymol()
    if not executable:
        return {"status": "unavailable", "error": "PyMOL was not found."}
    result = proteus_common.run_command([executable, "-c", "-q", str(Path(script_path).resolve())], timeout=timeout)
    return execution_summary(result)


def run_chimerax(commands: list[str], *, timeout: int = 600) -> dict[str, Any]:
    executable = find_chimerax()
    if not executable:
        return {"status": "unavailable", "error": "ChimeraX was not found."}
    env = os.environ.copy()
    env["CHIMERAX_BIN"] = executable
    command_text = ";".join(commands)
    result = proteus_common.run_command(
        [sys.executable, str(SCRIPT_DIR / "chimerax_rest.py"), "run", command_text],
        timeout=timeout,
        env=env,
    )
    if result.get("returncode") == 0:
        try:
            payload = json.loads(result.get("stdout") or "{}")
        except json.JSONDecodeError:
            payload = None
        if isinstance(payload, dict) and payload.get("status") == "error":
            result["status"] = "error"
            result["error"] = payload.get("error", "ChimeraX workflow failed.")
        elif payload is not None:
            result["payload"] = payload
    return execution_summary(result)


def execution_summary(result: dict[str, Any]) -> dict[str, Any]:
    """Return useful execution diagnostics without embedding commands or raw data."""
    cleaned = proteus_common.scrub_private(result)
    output: dict[str, Any] = {"status": cleaned.get("status", "error")}
    if "returncode" in cleaned:
        output["returncode"] = cleaned.get("returncode")
    if cleaned.get("error"):
        output["error"] = cleaned["error"]
    payload = cleaned.get("payload")
    if isinstance(payload, dict):
        history = payload.get("data", {}).get("history") if isinstance(payload.get("data"), dict) else None
        if isinstance(history, list):
            output["commands_completed"] = len(history)
            output["elapsed_seconds"] = round(sum(float(item.get("elapsed_seconds") or 0) for item in history if isinstance(item, dict)), 3)
    if output["status"] != "ok":
        for key in ("stderr", "stdout"):
            value = cleaned.get(key)
            if isinstance(value, str) and value.strip():
                output[key] = value[-4000:]
    return output


def verify_nonempty(path: str | Path, *, attempts: int = 8, delay: float = 0.25) -> bool:
    candidate = Path(path)
    for _ in range(attempts):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return True
        time.sleep(delay)
    return False


def pymol_base(*, background: str = "white", width: int = 1200, height: int = 900) -> list[str]:
    color(background)
    return [
        "reinitialize",
        "set auto_zoom, 0",
        f"bg_color {background}",
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        "set cartoon_flat_sheets, 1",
        "set ray_shadows, 1",
        "set specular, 0.25",
        "set ambient, 0.42",
        f"viewport {int(width)}, {int(height)}",
    ]


def chimerax_base(*, background: str = "white") -> list[str]:
    color(background)
    return [f"set bgColor {background}", "lighting soft"]


def finalize_pymol(output: str | Path, *, width: int = 1200, height: int = 900, session: str | Path | None = None) -> list[str]:
    lines = ["orient", f"ray {int(width)}, {int(height)}", *pymol_png(output, width=width, height=height)]
    if session is not None:
        lines.extend(pymol_save(session))
    lines.append("quit")
    return lines


def finalize_chimerax(output: str | Path, *, width: int = 1200, height: int = 900, session: str | Path | None = None) -> list[str]:
    lines = ["view", f"save {quote_chimerax(Path(output).resolve())} width {int(width)} height {int(height)} supersample 3", "wait 1"]
    if session is not None:
        lines.append(f"save {quote_chimerax(Path(session).resolve())}")
    return lines
