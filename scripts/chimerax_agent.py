#!/usr/bin/env python3
"""ChimeraX headless agent helper — Proteus skill.

Runs ChimeraX commands headlessly, capturing structured output.
Designed for AI agent workflows.

NOTE: ChimeraX --nogui cannot render images on macOS (no OpenGL).
Use PyMOL for headless rendering, or the REST API for GUI-based rendering.

Usage:
    python chimerax_agent.py run "open 1ubq from pdb; info chains #1"
    python chimerax_agent.py info structure.pdb
    python chimerax_agent.py align ref.pdb mobile.pdb
    python chimerax_agent.py sasa structure.pdb
    python chimerax_agent.py hbonds complex.pdb --chain1 A --chain2 B
    python chimerax_agent.py --help
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

def _find_chimerax() -> str:
    """Auto-detect ChimeraX binary."""
    found = shutil.which("ChimeraX") or shutil.which("chimerax")
    if found:
        return found
    # macOS: search /Applications for latest version
    hits = glob.glob("/Applications/ChimeraX*.app/Contents/bin/ChimeraX")
    if hits:
        return sorted(hits)[-1]
    # Linux common paths
    for p in ["/usr/bin/chimerax", "/usr/local/bin/chimerax",
              os.path.expanduser("~/ChimeraX/bin/ChimeraX")]:
        if os.path.isfile(p):
            return p
    return None


CHIMERAX = os.environ.get("CHIMERAX_BIN") or _find_chimerax()


def _indent(text: str, spaces: int = 4) -> str:
    prefix = " " * spaces
    return "\n".join(prefix + line for line in text.split("\n"))


def _parse_output(stdout: str) -> tuple[list[str], list[str]]:
    """Parse ChimeraX stdout into info lines and errors.

    ChimeraX prefixes output with INFO:, WARNING:, ERROR:, STATUS:.
    Lines starting with 'INFO: Executing:' are echoed commands, not results.
    """
    info_lines = []
    errors = []
    for line in stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("ERROR:"):
            errors.append(stripped[6:].strip())
        elif stripped.startswith("WARNING:"):
            # Surface warnings — ChimeraX warns about missing atoms, format issues, etc.
            info_lines.append(f"[WARNING] {stripped[8:].strip()}")
        elif stripped.startswith("INFO:"):
            content = stripped[5:].strip()
            # Skip echoed commands
            if content and not content.startswith("Executing:"):
                info_lines.append(content)
        elif stripped.startswith("STATUS:"):
            info_lines.append(stripped[7:].strip())
    return info_lines, errors


def _finalize_process_result(proc: subprocess.CompletedProcess, output_path: str) -> dict:
    """Combine wrapper JSON with subprocess diagnostics."""
    payload = None
    if os.path.exists(output_path):
        with open(output_path) as fh:
            payload = json.load(fh)

    stderr = proc.stderr.strip()

    if payload is None:
        info_lines, errors = _parse_output(proc.stdout)
        result = {"status": "error", "error": "No output file produced"}
        if info_lines:
            result["info"] = info_lines
        if errors:
            result["errors"] = errors
    else:
        result = payload

    if proc.returncode != 0:
        if result.get("status") == "ok":
            result = {
                "status": "error",
                "error": f"ChimeraX exited with code {proc.returncode}",
                "data": result.get("data", {}),
            }
        result["returncode"] = proc.returncode
        if stderr:
            result.setdefault("stderr", stderr)
    elif payload is None and stderr:
        result["stderr"] = stderr
    elif result.get("status") == "error" and stderr:
        result.setdefault("stderr", stderr)

    return result


def run_chimerax_commands(commands: str, timeout: int = 120) -> dict:
    """Run ChimeraX commands headlessly and capture output.

    Commands are semicolon-separated. Analysis commands work in --nogui mode.
    Image rendering does NOT work in --nogui mode on macOS.
    """
    if not CHIMERAX:
        return {"status": "error", "error": "ChimeraX not found. Install it or set CHIMERAX_BIN."}

    cmd_list = [c.strip() for c in commands.replace("\n", ";").split(";") if c.strip()]
    cmd_str = " ; ".join(cmd_list)

    try:
        proc = subprocess.run(
            [CHIMERAX, "--nogui", "--exit", "--cmd", cmd_str],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        info_lines, errors = _parse_output(proc.stdout)
        result = {
            "status": "error" if errors else "ok",
            "info": info_lines,
            "errors": errors if errors else None,
        }
        if proc.returncode != 0:
            result["status"] = "error"
            result["returncode"] = proc.returncode
            stderr = proc.stderr.strip()
            if stderr:
                result["stderr"] = stderr
        return result
    except subprocess.TimeoutExpired:
        return {"status": "error", "error": f"Timeout after {timeout}s"}
    except FileNotFoundError:
        return {"status": "error", "error": f"ChimeraX binary not found at {CHIMERAX}"}


def run_chimerax_python(script_content: str, timeout: int = 120) -> dict:
    """Run a Python script inside ChimeraX headlessly.

    The script has access to `session` and can use `rc(session, "command")`
    to run ChimeraX commands. Assign to _output["data"] to return values.
    """
    if not CHIMERAX:
        return {"status": "error", "error": "ChimeraX not found. Install it or set CHIMERAX_BIN."}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as out_f:
        output_path = out_f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        wrapper = f'''
import json

_output = {{"status": "ok", "data": {{}}}}
_outpath = {output_path!r}

try:
    from chimerax.core.commands import run as rc
    # `session` is a global injected by ChimeraX's --script runner

{_indent(script_content, 4)}

except Exception as e:
    _output["status"] = "error"
    _output["error"] = str(e)

with open(_outpath, "w") as _f:
    json.dump(_output, _f, indent=2, default=str)
'''
        f.write(wrapper)
        f.flush()
        script_path = f.name

    try:
        proc = subprocess.run(
            [CHIMERAX, "--nogui", "--exit", "--script", script_path],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return _finalize_process_result(proc, output_path)
    except subprocess.TimeoutExpired:
        return {"status": "error", "error": f"Timeout after {timeout}s"}
    except FileNotFoundError:
        return {"status": "error", "error": f"ChimeraX binary not found at {CHIMERAX}"}
    finally:
        os.unlink(script_path)
        if os.path.exists(output_path):
            os.unlink(output_path)


def get_structure_info(pdb_path: str) -> dict:
    """Load a structure and return basic info via ChimeraX."""
    abs_path = os.path.abspath(pdb_path)
    commands = f"open {abs_path} ; info chains #1 ; info models #1"
    return run_chimerax_commands(commands)


def align_structures(pdb1: str, pdb2: str) -> dict:
    """Align two structures and return RMSD via matchmaker."""
    abs1, abs2 = os.path.abspath(pdb1), os.path.abspath(pdb2)
    commands = f"open {abs1} ; open {abs2} ; matchmaker #2 to #1"
    return run_chimerax_commands(commands)


def measure_sasa(pdb_path: str) -> dict:
    """Measure solvent-accessible surface area."""
    abs_path = os.path.abspath(pdb_path)
    commands = f"open {abs_path} ; measure sasa #1"
    return run_chimerax_commands(commands)


def find_hbonds(pdb_path: str, chain1: str = "A", chain2: str = "B") -> dict:
    """Find hydrogen bonds between two chains."""
    abs_path = os.path.abspath(pdb_path)
    commands = f"open {abs_path} ; hbonds #1/{chain1} restrict #1/{chain2} log true"
    return run_chimerax_commands(commands)


def main():
    parser = argparse.ArgumentParser(
        description="ChimeraX headless agent helper — analysis, alignment, measurements.",
        epilog="NOTE: Image rendering requires a GUI. Use pymol_agent.py render for headless.\n\n"
               "Examples:\n"
               "  %(prog)s run 'open 1ubq from pdb; info chains #1'\n"
               "  %(prog)s info structure.pdb\n"
               "  %(prog)s align reference.pdb mobile.pdb\n"
               "  %(prog)s sasa structure.pdb\n"
               "  %(prog)s hbonds complex.pdb --chain1 A --chain2 D",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # run
    p_run = sub.add_parser("run", help="Run ChimeraX commands (headless)")
    p_run.add_argument("commands", help="Semicolon-separated ChimeraX commands")

    # info
    p_info = sub.add_parser("info", help="Inspect a structure file")
    p_info.add_argument("pdb", help="Path to PDB/CIF file")

    # align
    p_align = sub.add_parser("align", help="Align two structures (matchmaker)")
    p_align.add_argument("reference", help="Reference structure (stays fixed)")
    p_align.add_argument("mobile", help="Mobile structure (gets moved)")

    # sasa
    p_sasa = sub.add_parser("sasa", help="Measure solvent-accessible surface area")
    p_sasa.add_argument("pdb", help="Path to structure file")

    # hbonds
    p_hb = sub.add_parser("hbonds", help="Find H-bonds between two chains")
    p_hb.add_argument("pdb", help="Path to structure file")
    p_hb.add_argument("--chain1", default="A", help="First chain (default: A)")
    p_hb.add_argument("--chain2", default="B", help="Second chain (default: B)")

    args = parser.parse_args()

    if args.command == "run":
        result = run_chimerax_commands(args.commands)
    elif args.command == "info":
        result = get_structure_info(args.pdb)
    elif args.command == "align":
        result = align_structures(args.reference, args.mobile)
    elif args.command == "sasa":
        result = measure_sasa(args.pdb)
    elif args.command == "hbonds":
        result = find_hbonds(args.pdb, args.chain1, args.chain2)
    else:
        parser.print_help()
        sys.exit(1)

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
