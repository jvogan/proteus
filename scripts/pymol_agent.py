#!/usr/bin/env python3
"""PyMOL headless agent helper — Proteus skill.

Runs PyMOL commands and scripts headlessly, capturing structured JSON output.
Designed for AI agent workflows where stdout is unreliable.

Usage:
    python pymol_agent.py run "fetch 1ubq; show cartoon"
    python pymol_agent.py info structure.pdb
    python pymol_agent.py render structure.pdb output.png
    python pymol_agent.py --help
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

def _find_pymol() -> str:
    """Auto-detect PyMOL binary. Checks PATH, then common install locations."""
    found = shutil.which("pymol")
    if found:
        return found
    # macOS common locations
    for pattern in ["/Applications/PyMOL.app/Contents/bin/pymol",
                    os.path.expanduser("~/Applications/PyMOL.app/Contents/bin/pymol")]:
        for p in glob.glob(pattern):
            if os.path.isfile(p):
                return p
    # Linux
    for p in ["/usr/bin/pymol", "/usr/local/bin/pymol"]:
        if os.path.isfile(p):
            return p
    return None


PYMOL = os.environ.get("PYMOL_BIN") or _find_pymol()


def _indent(text: str, spaces: int = 4) -> str:
    prefix = " " * spaces
    return "\n".join(prefix + line for line in text.split("\n"))


def _finalize_process_result(proc: subprocess.CompletedProcess, output_path: str) -> dict:
    """Combine wrapper JSON with subprocess diagnostics."""
    payload = None
    if os.path.exists(output_path):
        with open(output_path) as fh:
            payload = json.load(fh)

    stderr = proc.stderr.strip()
    stdout = proc.stdout.strip()

    if payload is None:
        result = {"status": "error", "error": "No output file produced"}
    else:
        result = payload

    if proc.returncode != 0:
        if result.get("status") == "ok":
            result = {
                "status": "error",
                "error": f"PyMOL exited with code {proc.returncode}",
                "data": result.get("data", {}),
            }
        result["returncode"] = proc.returncode
        if stderr:
            result.setdefault("stderr", stderr)
        if stdout:
            result.setdefault("stdout", stdout)
    elif payload is None:
        if stderr:
            result["stderr"] = stderr
        if stdout:
            result["stdout"] = stdout
    elif result.get("status") == "error" and stderr:
        result.setdefault("stderr", stderr)

    return result


def run_pymol_script(script_content: str, timeout: int = 120) -> dict:
    """Run a PyMOL Python script headlessly and capture output as JSON.

    The script can assign values to _output["data"] to return structured data.
    Example: _output["data"]["rmsd"] = 1.23
    """
    if not PYMOL:
        return {"status": "error", "error": "PyMOL not found. Install it or set PYMOL_BIN."}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as out_f:
        output_path = out_f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        wrapper = f'''
import json, sys, os
_output = {{"status": "ok", "data": {{}}}}
_outpath = {output_path!r}

try:
    from pymol import cmd, util
{_indent(script_content, 4)}
except Exception as e:
    _output["status"] = "error"
    _output["error"] = str(e)
finally:
    with open(_outpath, "w") as _f:
        json.dump(_output, _f, indent=2, default=str)
    try:
        cmd.quit()
    except Exception:
        pass
'''
        f.write(wrapper)
        f.flush()
        script_path = f.name

    try:
        proc = subprocess.run(
            [PYMOL, "-c", "-q", "-r", script_path],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return _finalize_process_result(proc, output_path)
    except subprocess.TimeoutExpired:
        return {"status": "error", "error": f"Timeout after {timeout}s"}
    except FileNotFoundError:
        return {"status": "error", "error": f"PyMOL binary not found at {PYMOL}"}
    finally:
        os.unlink(script_path)
        if os.path.exists(output_path):
            os.unlink(output_path)


def run_pymol_commands(commands: str, timeout: int = 120) -> dict:
    """Run PyMOL commands (not Python — PyMOL command language).

    Wraps each line with cmd.do() for execution. For complex selections
    containing >, <, use this instead of the -d CLI flag (which breaks
    on shell metacharacters).
    """
    lines = commands.strip().split("\n")
    script = "\n".join(f"cmd.do({line.strip()!r})" for line in lines if line.strip())
    return run_pymol_script(script, timeout)


def get_structure_info(pdb_path: str) -> dict:
    """Load a structure and return basic info (chains, atoms, B-factors)."""
    abs_path = os.path.abspath(pdb_path)
    script = f'''
cmd.load("{abs_path}", "struct")
_output["data"]["names"] = cmd.get_names()
_output["data"]["atom_count"] = cmd.count_atoms("all")
_output["data"]["chains"] = cmd.get_chains("all")

# Residue count per chain
for ch in cmd.get_chains("all"):
    sel = f"chain {{ch}} and name CA"
    _output["data"][f"chain_{{ch}}_residues"] = cmd.count_atoms(sel)

# B-factor stats (pLDDT for AlphaFold structures)
stored_b = []
cmd.iterate("name CA", "stored_b.append(b)", space={{"stored_b": stored_b}})
if stored_b:
    _output["data"]["bfactor_min"] = round(min(stored_b), 2)
    _output["data"]["bfactor_max"] = round(max(stored_b), 2)
    _output["data"]["bfactor_mean"] = round(sum(stored_b) / len(stored_b), 2)
'''
    return run_pymol_script(script)


def render_structure(pdb_path: str, output_png: str, width: int = 1200, height: int = 900,
                     style: str = "cartoon", color: str = "spectrum") -> dict:
    """Load and render a structure to PNG using PyMOL's software ray tracer.

    Works fully headless — no display required.

    Args:
        style: cartoon, sticks, surface, spheres, lines
        color: spectrum (rainbow), bfactor (blue-white-red), chain, or any PyMOL color name
    """
    abs_pdb = os.path.abspath(pdb_path)
    abs_out = os.path.abspath(output_png)
    script = f'''
cmd.load("{abs_pdb}", "struct")
cmd.hide("everything")
cmd.show("{style}", "all")
if "{color}" == "spectrum":
    cmd.spectrum("count", "rainbow", "all")
elif "{color}" == "bfactor":
    cmd.spectrum("b", "blue_white_red", "all")
elif "{color}" == "chain":
    util.cbc("all")
else:
    cmd.color("{color}", "all")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.orient()
cmd.ray({width}, {height})
cmd.png("{abs_out}")
_output["data"]["rendered"] = "{abs_out}"
_output["data"]["size"] = "{width}x{height}"
'''
    return run_pymol_script(script, timeout=300)  # Rendering can take longer


def main():
    parser = argparse.ArgumentParser(
        description="PyMOL headless agent helper — run commands, inspect structures, render images.",
        epilog="Examples:\n"
               "  %(prog)s run 'fetch 1ubq; show cartoon'\n"
               "  %(prog)s info structure.pdb\n"
               "  %(prog)s render structure.pdb output.png\n"
               "  %(prog)s render structure.pdb output.png --style surface --color chain",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # run
    p_run = sub.add_parser("run", help="Run PyMOL commands")
    p_run.add_argument("commands", help="Semicolon-separated PyMOL commands")

    # info
    p_info = sub.add_parser("info", help="Inspect a structure file")
    p_info.add_argument("pdb", help="Path to PDB/CIF/SDF file")

    # render
    p_render = sub.add_parser("render", help="Render structure to PNG (headless)")
    p_render.add_argument("pdb", help="Path to structure file")
    p_render.add_argument("output", nargs="?", default="/tmp/pymol_render.png", help="Output PNG path")
    p_render.add_argument("--width", type=int, default=1200)
    p_render.add_argument("--height", type=int, default=900)
    p_render.add_argument("--style", default="cartoon", choices=["cartoon", "sticks", "surface", "spheres", "lines"])
    p_render.add_argument("--color", default="spectrum", help="spectrum, bfactor, chain, or PyMOL color name")

    args = parser.parse_args()

    if args.command == "run":
        result = run_pymol_commands(args.commands)
    elif args.command == "info":
        result = get_structure_info(args.pdb)
    elif args.command == "render":
        result = render_structure(args.pdb, args.output, args.width, args.height, args.style, args.color)
    else:
        parser.print_help()
        sys.exit(1)

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
