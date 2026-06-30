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
import re
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


def _py_literal(value: str) -> str:
    return repr(value)


def _validate_pymol_color(color: str) -> str:
    if color in {"spectrum", "bfactor", "chain", "plddt"}:
        return color
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", color):
        raise ValueError(
            "Unsafe PyMOL color name. Use spectrum, bfactor, chain, plddt, or a simple PyMOL color identifier."
        )
    return color


def _color_script(color_mode: str, selection: str = "all") -> str:
    """PyMOL command lines that apply a color mode to a selection."""
    sel = _py_literal(selection)
    if color_mode == "spectrum":
        return f'cmd.spectrum("count", "rainbow", {sel})'
    if color_mode == "bfactor":
        return f'cmd.spectrum("b", "blue_white_red", {sel})'
    if color_mode == "chain":
        return f'util.cbc({sel})'
    if color_mode == "plddt":
        # Official AlphaFold bins. Detect normalized 0-1 confidence values as
        # used by some predictors, then layer broadest-first because PyMOL
        # selection algebra has no `<=`.
        return "\n".join([
            f"_proteus_selection = {sel}",
            "_proteus_b_values = []",
            'cmd.iterate(_proteus_selection, "_proteus_b_values.append(b)", space={"_proteus_b_values": _proteus_b_values})',
            "_proteus_scale01 = bool(_proteus_b_values and max(_proteus_b_values) <= 1.5)",
            "_proteus_thresholds = (0.50, 0.70, 0.90) if _proteus_scale01 else (50.0, 70.0, 90.0)",
            'cmd.color("orange", _proteus_selection)',
            'cmd.color("yellow", f"({_proteus_selection}) and b > {_proteus_thresholds[0]}")',
            'cmd.color("cyan", f"({_proteus_selection}) and b > {_proteus_thresholds[1]}")',
            'cmd.color("blue", f"({_proteus_selection}) and b > {_proteus_thresholds[2]}")',
            "try:",
            '    _output["data"]["plddt_scale"] = "0-1" if _proteus_scale01 else "0-100"',
            "except Exception:",
            "    pass",
        ])
    return f'cmd.color({_py_literal(color_mode)}, {sel})'


def _preset_script(preset: str) -> str:
    """PyMOL command lines for a render look. Default: publication."""
    if preset == "illustration":
        return "\n".join([
            'cmd.bg_color("white")',
            'cmd.set("ray_opaque_background", 1)',
            'cmd.set("antialias", 2)',
            'cmd.set("cartoon_fancy_helices", 1)',
            'cmd.set("cartoon_smooth_loops", 1)',
            'cmd.set("cartoon_flat_sheets", 1)',
            'cmd.set("ray_trace_mode", 3)',     # quantized colors + black outlines
            'cmd.set("ray_trace_color", "black")',
        ])
    if preset == "soft":
        return "\n".join([
            'cmd.bg_color("gray90")',
            'cmd.set("ray_opaque_background", 1)',
            'cmd.set("orthoscopic", 1)',
            'cmd.set("ray_shadows", 0)',
            'cmd.set("antialias", 2)',
            'cmd.set("ambient", 0.4)',
            'cmd.set("specular", 0.15)',
            'cmd.set("cartoon_fancy_helices", 1)',
            'cmd.set("cartoon_smooth_loops", 1)',
            'cmd.set("cartoon_flat_sheets", 1)',
        ])
    return "\n".join([
        'cmd.bg_color("white")',
        'cmd.set("ray_opaque_background", 1)',
        'cmd.set("antialias", 2)',
        'cmd.set("ray_shadows", 1)',
        'cmd.set("specular", 0.25)',
        'cmd.set("ambient", 0.35)',
        'cmd.set("cartoon_fancy_helices", 1)',
        'cmd.set("cartoon_smooth_loops", 1)',
        'cmd.set("cartoon_flat_sheets", 1)',
    ])


def _encode_movie(frame_dir: str, output: str, fps: int = 30) -> str:
    """Encode frame_%04d.png frames to MP4 (or GIF if output ends .gif).

    yuv420p + even-dimension scaling keeps the MP4 web-playable; GIF uses a
    two-pass palette for clean colors. Raises CalledProcessError on failure.
    """
    out = os.path.abspath(output)
    pattern = os.path.join(frame_dir, "frame_%04d.png")
    even = "scale=trunc(iw/2)*2:trunc(ih/2)*2"
    if out.lower().endswith(".gif"):
        palette = os.path.join(frame_dir, "palette.png")
        subprocess.run(["ffmpeg", "-y", "-framerate", str(fps), "-i", pattern,
                        "-vf", f"{even},palettegen", palette],
                       check=True, capture_output=True, text=True)
        subprocess.run(["ffmpeg", "-y", "-framerate", str(fps), "-i", pattern,
                        "-i", palette, "-lavfi", f"{even} [x]; [x][1:v] paletteuse", out],
                       check=True, capture_output=True, text=True)
    else:
        subprocess.run(["ffmpeg", "-y", "-framerate", str(fps), "-i", pattern,
                        "-c:v", "libx264", "-preset", "slow", "-crf", "18",
                        "-pix_fmt", "yuv420p", "-movflags", "+faststart",
                        "-vf", even, out],
                       check=True, capture_output=True, text=True)
    return out


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
structure_path = {_py_literal(abs_path)}
cmd.load(structure_path, "struct")
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


def _verify_png(result: dict, output_path: str) -> dict:
    """Downgrade an 'ok' result to error if the PNG is missing or empty.

    PyMOL can return without writing pixels (an empty selection, a silently
    failed render). Fail loudly instead of reporting a blank success.
    """
    if result.get("status") != "ok":
        return result
    if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
        return {"status": "error",
                "error": f"Render produced no image (empty or missing file: {output_path})",
                "data": result.get("data", {})}
    return result


def render_structure(pdb_path: str, output_png: str, width: int = 1200, height: int = 900,
                     style: str = "cartoon", color: str = "spectrum",
                     preset: str = "publication") -> dict:
    """Load and render a structure to PNG using PyMOL's software ray tracer.

    Works fully headless — no display required.

    Args:
        style: cartoon, sticks, surface, spheres, lines
        color: spectrum (rainbow), bfactor (blue-white-red), chain, plddt
               (AlphaFold confidence bins), or any PyMOL color name
        preset: publication, illustration (outlined), or soft (neutral background)
    """
    try:
        color = _validate_pymol_color(color)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}

    abs_pdb = os.path.abspath(pdb_path)
    abs_out = os.path.abspath(output_png)
    script = f'''
cmd.load({_py_literal(abs_pdb)}, "struct")
cmd.hide("everything")
cmd.show({_py_literal(style)}, "all")
{_color_script(color)}
{_preset_script(preset)}
cmd.orient()
cmd.ray({width}, {height})
cmd.png({_py_literal(abs_out)})
_output["data"]["rendered"] = {_py_literal(abs_out)}
_output["data"]["size"] = "{width}x{height}"
'''
    result = run_pymol_script(script, timeout=300)  # Rendering can take longer
    return _verify_png(result, abs_out)


def render_spin(pdb_path: str, output: str, frames: int = 60, width: int = 800,
                height: int = 600, style: str = "cartoon", color: str = "spectrum",
                preset: str = "publication", fps: int = 30) -> dict:
    """Render a 360-degree y-spin as ray-traced frames, then encode to a movie.

    This is the headless-correct turntable path: PyMOL ray-traces each frame
    (works without a display) and ffmpeg encodes them. Degrades gracefully —
    if ffmpeg is missing, the frames are written and their directory returned.
    """
    if not PYMOL:
        return {"status": "error", "error": "PyMOL not found. Install it or set PYMOL_BIN."}
    if not os.path.isfile(pdb_path):
        return {"status": "error", "error": f"Structure file not found: {pdb_path}"}
    try:
        color = _validate_pymol_color(color)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}

    ffmpeg = shutil.which("ffmpeg")
    frame_dir = tempfile.mkdtemp(prefix="proteus_spin_")
    abs_pdb = os.path.abspath(pdb_path)
    script = f'''
cmd.load({_py_literal(abs_pdb)}, "struct")
cmd.hide("everything")
cmd.show({_py_literal(style)}, "all")
{_color_script(color)}
{_preset_script(preset)}
cmd.set("cache_frames", 0)   # else PyMOL caches every frame in RAM -> OOM on long spins
cmd.set("orthoscopic", 1)    # stable apparent size through the rotation
cmd.orient()
# No-clip framing: fit the bounding SPHERE (rotation-invariant) and widen the
# depth slab so nothing clips as wide axes rotate toward the camera.
cmd.zoom("all", buffer=4, complete=1)
_ext = cmd.get_extent("all")
if _ext:
    _dx = _ext[1][0] - _ext[0][0]
    _dy = _ext[1][1] - _ext[0][1]
    _dz = _ext[1][2] - _ext[0][2]
    cmd.clip("slab", (_dx * _dx + _dy * _dy + _dz * _dz) ** 0.5 * 1.6)
_n = {frames}
for _i in range(_n):
    cmd.ray({width}, {height})
    cmd.png(os.path.join({_py_literal(frame_dir)}, "frame_%04d.png" % _i))
    cmd.turn("y", 360.0 / _n)
_output["data"]["frames"] = _n
'''
    result = run_pymol_script(script, timeout=max(600, frames * 20))
    if result.get("status") != "ok":
        shutil.rmtree(frame_dir, ignore_errors=True)
        return result
    if not ffmpeg:
        return {"status": "ok", "data": {
            "movie": None, "frames_dir": frame_dir, "frame_count": frames,
            "note": "ffmpeg not found; wrote frames but did not encode a movie.",
        }}
    try:
        _encode_movie(frame_dir, output, fps)
    except subprocess.CalledProcessError as exc:
        return {"status": "error", "error": f"ffmpeg failed: {(exc.stderr or '')[:300]}",
                "data": {"frames_dir": frame_dir}}
    shutil.rmtree(frame_dir, ignore_errors=True)
    return {"status": "ok", "data": {
        "movie": os.path.abspath(output), "frame_count": frames, "fps": fps,
    }}


def render_pocket(pdb_path: str, output: str, ligand: str = "organic", radius: float = 5.0,
                  width: int = 1200, height: int = 900, label: bool = False,
                  preset: str = "publication") -> dict:
    """Render an annotated binding-pocket figure.

    Shows the ligand and the residues within `radius` as sticks with polar
    contacts drawn, the rest of the protein as a transparent cartoon for
    context, framed on the pocket. `ligand` is any PyMOL selection (default:
    `organic`, i.e. non-polymer/non-solvent ligands).
    """
    if not PYMOL:
        return {"status": "error", "error": "PyMOL not found. Install it or set PYMOL_BIN."}
    if not os.path.isfile(pdb_path):
        return {"status": "error", "error": f"Structure file not found: {pdb_path}"}

    abs_pdb = os.path.abspath(pdb_path)
    abs_out = os.path.abspath(output)
    label_cmd = 'cmd.label("pocket and name CA", "resn+resi")' if label else ""
    script = f'''
cmd.load({_py_literal(abs_pdb)}, "struct")
cmd.select("ligand", "(%s)" % {_py_literal(ligand)})
if cmd.count_atoms("ligand") == 0:
    raise RuntimeError("No ligand atoms matched selection %r; pass --ligand" % {_py_literal(ligand)})
cmd.select("pocket", "byres (polymer within {radius} of ligand)")
cmd.hide("everything")
cmd.show("cartoon", "polymer")
cmd.set("cartoon_transparency", 0.55)
cmd.show("sticks", "pocket")
cmd.show("sticks", "ligand")
cmd.set("stick_radius", 0.18, "pocket")
cmd.color("gray70", "pocket and elem C")
cmd.color("yellow", "ligand and elem C")
cmd.do("util.cnc('pocket or ligand')")
cmd.distance("contacts", "ligand", "pocket", 3.5, mode=2)
cmd.set("dash_color", "yellow")
cmd.set("dash_width", 2.5)
cmd.hide("labels", "contacts")
{label_cmd}
{_preset_script(preset)}
cmd.orient("ligand")
cmd.zoom("ligand", {radius} + 3)
cmd.ray({width}, {height})
cmd.png({_py_literal(abs_out)})
_output["data"]["rendered"] = {_py_literal(abs_out)}
_output["data"]["ligand_selection"] = {_py_literal(ligand)}
_output["data"]["ligand_atoms"] = cmd.count_atoms("ligand")
_output["data"]["pocket_residues"] = cmd.count_atoms("pocket and name CA")
'''
    result = run_pymol_script(script, timeout=300)
    return _verify_png(result, abs_out)


def render_density(pdb_path: str, output: str, map_path: str = None, simulate: bool = False,
                   level: float = None, n_sigma: float = 2.0, carve: float = 2.5,
                   residue: str = None, color: str = "chain", width: int = 1200,
                   height: int = 900, preset: str = "publication") -> dict:
    """Render a model fitted in cryo-EM density (real map or simulated).

    The density mesh is carved around the model (or around `residue`), which
    avoids the whole-map contour stall in headless PyMOL. With `residue`, shows
    that selection as sticks and zooms on it — the "is this sidechain supported
    by density?" figure. Without a map, `simulate` generates a gaussian density
    from the model so the same path works for predicted/designed structures.

    The contour level defaults to mean + n_sigma * sigma from the map (via
    map_info), or 1.0 for simulated density.
    """
    # Validate inputs before checking for the tool, so argument errors are
    # deterministic whether or not PyMOL is installed.
    if not os.path.isfile(pdb_path):
        return {"status": "error", "error": f"Structure file not found: {pdb_path}"}
    if not simulate and not map_path:
        return {"status": "error", "error": "Provide --map, or use --simulate to generate density from the model."}
    if map_path and not os.path.isfile(map_path):
        return {"status": "error", "error": f"Map file not found: {map_path}"}
    if not PYMOL:
        return {"status": "error", "error": "PyMOL not found. Install it or set PYMOL_BIN."}
    try:
        color = _validate_pymol_color(color)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}

    # Contour level: explicit > sigma-from-map > simulated default.
    auto_level = None
    if level is None:
        if simulate:
            level = 1.0
        else:
            try:
                import map_info
                _, _, mean, sigma = map_info.read_map_stats(map_path)
                level = round(mean + n_sigma * sigma, 4)
                auto_level = {"mean": round(mean, 4), "sigma": round(sigma, 4), "n_sigma": n_sigma}
            except Exception:
                level = 1.0

    abs_pdb = os.path.abspath(pdb_path)
    abs_out = os.path.abspath(output)
    sel = "(%s)" % residue if residue else None
    carve_target = sel if sel else "struct"

    lines = [f'cmd.load({_py_literal(abs_pdb)}, "struct")']
    if simulate:
        lines.append('cmd.map_new("emap", "gaussian", 1.0, "struct", 5)')
    else:
        lines.append(f'cmd.load({_py_literal(os.path.abspath(map_path))}, "emap")')
    lines += [
        f'cmd.isomesh("dens", "emap", {level}, {_py_literal(carve_target)}, carve={carve})',
        'cmd.color("gray70", "dens")',
        'cmd.hide("everything", "struct")',
        'cmd.show("cartoon", "struct")',
        _color_script(color, "struct"),
    ]
    if sel:
        lines += [
            'cmd.set("cartoon_transparency", 0.6, "struct")',
            f'cmd.show("sticks", {_py_literal(sel)})',
            f'util.cnc({_py_literal(sel)})',
        ]
    lines.append(_preset_script(preset))
    if sel:
        lines += [f'cmd.orient({_py_literal(sel)})', f'cmd.zoom({_py_literal(sel)}, 6)']
    else:
        lines.append('cmd.orient("struct")')
    lines += [
        f'cmd.ray({width}, {height})',
        f'cmd.png({_py_literal(abs_out)})',
        f'_output["data"]["rendered"] = {_py_literal(abs_out)}',
    ]
    result = run_pymol_script("\n".join(lines), timeout=300)
    result = _verify_png(result, abs_out)
    if result.get("status") == "ok":
        result["data"]["level"] = level
        result["data"]["carve"] = carve
        result["data"]["simulated"] = simulate
        if auto_level:
            result["data"]["auto_level"] = auto_level
    return result


def main():
    parser = argparse.ArgumentParser(
        description="PyMOL headless agent helper — run commands, inspect structures, render images.",
        epilog="Examples:\n"
               "  %(prog)s run 'fetch 1ubq; show cartoon'\n"
               "  %(prog)s info structure.pdb\n"
               "  %(prog)s render structure.pdb output.png --color plddt\n"
               "  %(prog)s pocket 1HSG.pdb pocket.png --label\n"
               "  %(prog)s spin structure.pdb spin.mp4 --frames 60",
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
    p_render.add_argument("--color", default="spectrum", help="spectrum, bfactor, chain, plddt, or PyMOL color name")
    p_render.add_argument("--preset", default="publication", choices=["publication", "illustration", "soft"])

    # pocket
    p_pocket = sub.add_parser("pocket", help="Render an annotated binding-pocket figure")
    p_pocket.add_argument("pdb", help="Path to structure file")
    p_pocket.add_argument("output", nargs="?", default="/tmp/pymol_pocket.png", help="Output PNG path")
    p_pocket.add_argument("--ligand", default="organic", help="Ligand selection (default: organic)")
    p_pocket.add_argument("--radius", type=float, default=5.0, help="Pocket radius in Angstroms (default: 5)")
    p_pocket.add_argument("--width", type=int, default=1200)
    p_pocket.add_argument("--height", type=int, default=900)
    p_pocket.add_argument("--label", action="store_true", help="Label pocket residues (resn+resi)")
    p_pocket.add_argument("--preset", default="publication", choices=["publication", "illustration", "soft"])

    # density
    p_density = sub.add_parser("density", help="Render a model in cryo-EM density (real map or --simulate)")
    p_density.add_argument("pdb", help="Path to model (PDB/CIF)")
    p_density.add_argument("output", nargs="?", default="/tmp/pymol_density.png", help="Output PNG path")
    p_density.add_argument("--map", dest="map_path", help="MRC/CCP4 density map (omit with --simulate)")
    p_density.add_argument("--simulate", action="store_true", help="Simulate gaussian density from the model")
    p_density.add_argument("--level", type=float, help="Contour level (default: map sigma, or 1.0 simulated)")
    p_density.add_argument("--n-sigma", type=float, default=2.0, help="Sigma multiple for auto level (default: 2)")
    p_density.add_argument("--carve", type=float, default=2.5, help="Carve radius around the model in Angstroms")
    p_density.add_argument("--residue", help="Residue selection to show as sticks and zoom on")
    p_density.add_argument("--color", default="chain", help="Model color: spectrum, chain, bfactor, plddt, or a color")
    p_density.add_argument("--width", type=int, default=1200)
    p_density.add_argument("--height", type=int, default=900)
    p_density.add_argument("--preset", default="publication", choices=["publication", "illustration", "soft"])

    # spin
    p_spin = sub.add_parser("spin", help="Render a 360-degree turntable movie (frames -> ffmpeg)")
    p_spin.add_argument("pdb", help="Path to structure file")
    p_spin.add_argument("output", nargs="?", default="/tmp/pymol_spin.mp4", help="Output movie path (.mp4 or .gif)")
    p_spin.add_argument("--frames", type=int, default=60)
    p_spin.add_argument("--width", type=int, default=800)
    p_spin.add_argument("--height", type=int, default=600)
    p_spin.add_argument("--style", default="cartoon", choices=["cartoon", "sticks", "surface", "spheres", "lines"])
    p_spin.add_argument("--color", default="spectrum", help="spectrum, bfactor, chain, plddt, or PyMOL color name")
    p_spin.add_argument("--preset", default="publication", choices=["publication", "illustration", "soft"])
    p_spin.add_argument("--fps", type=int, default=30)

    args = parser.parse_args()

    if args.command == "run":
        result = run_pymol_commands(args.commands)
    elif args.command == "info":
        result = get_structure_info(args.pdb)
    elif args.command == "render":
        result = render_structure(args.pdb, args.output, args.width, args.height,
                                  args.style, args.color, args.preset)
    elif args.command == "pocket":
        result = render_pocket(args.pdb, args.output, args.ligand, args.radius,
                               args.width, args.height, args.label, args.preset)
    elif args.command == "density":
        result = render_density(args.pdb, args.output, args.map_path, args.simulate,
                                args.level, args.n_sigma, args.carve, args.residue,
                                args.color, args.width, args.height, args.preset)
    elif args.command == "spin":
        result = render_spin(args.pdb, args.output, args.frames, args.width, args.height,
                             args.style, args.color, args.preset, args.fps)
    else:
        parser.print_help()
        sys.exit(1)

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
