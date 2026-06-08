#!/usr/bin/env python3
"""ChimeraX REST render agent — Proteus skill.

Drives a *managed* ChimeraX GUI session over its REST API so an agent can render
images (and turntable movies) from the terminal without hand-managing the
process. ChimeraX `--nogui` cannot render on macOS (no OpenGL context); this
launches a real GUI instance on an ephemeral port, talks to it over HTTP, and
guarantees teardown.

Where chimerax_agent.py does headless analysis (`--nogui`), this does GPU
rendering (REST + GUI). It also defeats the macOS "0-byte PNG" save race.

Usage:
    python chimerax_rest.py render structure.pdb out.png
    python chimerax_rest.py render structure.pdb out.png --style surface --color bychain
    python chimerax_rest.py render model.pdb out.png --color plddt
    python chimerax_rest.py spin model.pdb spin.mp4 --frames 72        # needs ffmpeg
    python chimerax_rest.py run "open 1ubq from pdb; cartoon; color bychain"
    python chimerax_rest.py --help

Environment:
    CHIMERAX_BIN   Override the ChimeraX binary path.
"""

import argparse
import glob
import http.client
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from urllib.parse import urlencode


def _find_chimerax() -> str:
    """Auto-detect ChimeraX binary (PATH, then common install locations)."""
    found = shutil.which("ChimeraX") or shutil.which("chimerax")
    if found:
        return found
    hits = glob.glob("/Applications/ChimeraX*.app/Contents/bin/ChimeraX")
    if hits:
        return sorted(hits)[-1]
    for p in ["/usr/bin/chimerax", "/usr/local/bin/chimerax",
              os.path.expanduser("~/ChimeraX/bin/ChimeraX")]:
        if os.path.isfile(p):
            return p
    return None


CHIMERAX = os.environ.get("CHIMERAX_BIN") or _find_chimerax()
DEFAULT_WIDTH = 1200
DEFAULT_HEIGHT = 900


def find_free_port() -> int:
    """Bind an ephemeral port so parallel ChimeraX sessions never collide."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


def _quote(path: str) -> str:
    return '"' + os.path.abspath(path).replace('"', '\\"') + '"'


def _validate_color(color: str) -> str:
    if color in {"rainbow", "bychain", "bfactor", "plddt"}:
        return color
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", color):
        raise ValueError(
            "Unsafe ChimeraX color. Use rainbow, bychain, bfactor, plddt, or a simple color name."
        )
    return color


class ChimeraXRest:
    """Launch a local ChimeraX GUI with its REST server and drive it over HTTP.

    Use as a context manager so the process is always torn down:

        with ChimeraXRest() as rest:
            rest.run("open 1ubq from pdb")
            rest.save_image("/tmp/out.png")
    """

    def __init__(self, chimerax: str = None, port: int = None):
        self.chimerax = chimerax or CHIMERAX
        self.port = port or find_free_port()
        self.process = None
        self.history: list = []
        self._log_path = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *exc):
        self.stop()
        return False

    def start(self, ready_timeout: int = 60) -> None:
        if not self.chimerax:
            raise RuntimeError("ChimeraX not found. Install it or set CHIMERAX_BIN.")
        log_handle = tempfile.NamedTemporaryFile(
            mode="w", suffix=".chimerax.log", delete=False)
        self._log_path = log_handle.name
        # A GUI (OpenGL) context is required to render; do NOT pass --offscreen on
        # macOS (it hangs the Qt event loop). `json true` makes /run return a JSON
        # envelope so command-level errors are detectable even on HTTP 200.
        self.process = subprocess.Popen(
            [self.chimerax, "--cmd",
             f"remotecontrol rest start port {self.port} json true log false"],
            stdout=log_handle, stderr=subprocess.STDOUT, text=True,
        )
        deadline = time.time() + ready_timeout
        last_error = ""
        while time.time() < deadline:
            if self.process.poll() is not None:
                raise RuntimeError(
                    f"ChimeraX exited before REST startup; see {self._log_path}")
            try:
                self.run("version", timeout=5)
                return
            except Exception as exc:  # startup retry loop
                last_error = str(exc)
                time.sleep(1)
        self.stop()
        raise TimeoutError(
            f"ChimeraX REST did not start on port {self.port}: {last_error}")

    def run(self, command: str, *, timeout: int = 120, soft: bool = False) -> str:
        """Run one ChimeraX command over REST.

        `soft=True` swallows command/HTTP errors so a purely decorative command
        never aborts a scene. Errors are detected from the JSON envelope's
        `error` key, not just the HTTP status.
        """
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=timeout)
        path = "/run?" + urlencode({"command": command, "json": "true"})
        started = time.time()
        try:
            conn.request("GET", path)
            response = conn.getresponse()
            body = response.read().decode("utf-8", errors="replace")
            self.history.append({
                "command": command,
                "elapsed_seconds": round(time.time() - started, 3),
                "http_status": response.status,
            })
            if response.status >= 400:
                if soft:
                    return body
                raise RuntimeError(f"ChimeraX REST {response.status}: {body[:300]}")
            try:
                payload = json.loads(body)
            except json.JSONDecodeError:
                payload = {}
            if payload.get("error"):
                if soft:
                    return body
                raise RuntimeError(
                    f"ChimeraX command failed for {command!r}: {payload['error']}")
            return body
        finally:
            conn.close()

    def run_all(self, commands, *, timeout: int = 120) -> None:
        for command in commands:
            self.run(command, timeout=timeout)

    def save_image(self, path, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT,
                   supersample=3, attempts=4) -> None:
        """Save a PNG, defeating the macOS 0-byte race.

        REST `save` can return HTTP 200 before the GL framebuffer is flushed,
        leaving a 0-byte PNG (worse under heavy cartoon recompute). Issue
        `wait 1` to force a redraw, then poll the file for non-zero size,
        retrying the whole save a few times.
        """
        path = os.path.abspath(path)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        if os.path.exists(path):
            os.unlink(path)
        last = 0
        for _ in range(attempts):
            self.run("wait 1")
            self.run(f"save {_quote(path)} width {width} height {height} "
                     f"supersample {supersample}", timeout=240)
            for _ in range(16):
                if os.path.exists(path) and os.path.getsize(path) > 0:
                    return
                time.sleep(0.5)
            last = os.path.getsize(path) if os.path.exists(path) else 0
        raise RuntimeError(
            f"ChimeraX did not write a non-empty image after {attempts} attempts "
            f"(last size={last}): {path}")

    def stop(self) -> None:
        for command in ("remotecontrol rest stop", "exit"):
            try:
                self.run(command, timeout=5)
            except Exception:
                pass
        if self.process is not None:
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.terminate()
                try:
                    self.process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    self.process.wait(timeout=5)
        if self._log_path and os.path.exists(self._log_path):
            try:
                os.unlink(self._log_path)
            except OSError:
                pass


def _scene_commands(style: str, color: str) -> list:
    """A tasteful publication default: white bg, soft light, silhouettes."""
    cmds = [
        "set bgColor white",
        "lighting soft",
        "graphics silhouettes true",
        "set subdivision 3",
    ]
    styles = {
        "cartoon": ["hide #1 atoms", "cartoon #1"],
        "surface": ["surface #1"],
        "stick": ["hide #1 cartoon", "style #1 stick", "show #1 atoms"],
        "sphere": ["hide #1 cartoon", "style #1 sphere", "show #1 atoms"],
    }
    cmds += styles.get(style, styles["cartoon"])
    color_cmd = {
        "rainbow": "rainbow #1",
        "bychain": "color bychain #1",
        "bfactor": "color bfactor #1",
        "plddt": "color bfactor #1 palette alphafold",
    }.get(color, f"color #1 {color}")
    cmds.append(color_cmd)
    return cmds


def _encode_movie(frame_dir: str, output: str, fps: int = 30) -> str:
    """Encode frame_%04d.png in frame_dir to MP4 (or GIF if output ends .gif).

    yuv420p + even-dimension scaling keeps the MP4 web-playable; GIF uses a
    two-pass palette for clean colors.
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


def rest_run(commands: str) -> dict:
    """Open a managed session, run semicolon/newline-separated commands, tear down."""
    if not CHIMERAX:
        return {"status": "error", "error": "ChimeraX not found. Install it or set CHIMERAX_BIN."}
    cmd_list = [c.strip() for c in commands.replace("\n", ";").split(";") if c.strip()]
    try:
        with ChimeraXRest() as rest:
            rest.run_all(cmd_list)
            return {"status": "ok", "data": {"history": rest.history}}
    except Exception as exc:
        return {"status": "error", "error": str(exc)}


def rest_render(structure: str, output: str, style: str = "cartoon", color: str = "rainbow",
                width: int = DEFAULT_WIDTH, height: int = DEFAULT_HEIGHT,
                supersample: int = 3) -> dict:
    """Render a structure to PNG via a managed ChimeraX GUI session (GPU)."""
    if not CHIMERAX:
        return {"status": "error", "error": "ChimeraX not found. Install it or set CHIMERAX_BIN."}
    if not os.path.isfile(structure):
        return {"status": "error", "error": f"Structure file not found: {structure}"}
    try:
        color = _validate_color(color)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}
    try:
        with ChimeraXRest() as rest:
            rest.run(f"open {_quote(structure)}")
            rest.run_all(_scene_commands(style, color))
            rest.run("view")
            rest.save_image(output, width, height, supersample)
        return {"status": "ok", "data": {
            "rendered": os.path.abspath(output),
            "size": f"{width}x{height}",
            "supersample": supersample,
        }}
    except Exception as exc:
        return {"status": "error", "error": str(exc)}


def rest_spin(structure: str, output: str, frames: int = 72, style: str = "cartoon",
              color: str = "rainbow", width: int = 800, height: int = 600,
              fps: int = 30) -> dict:
    """Render a 360-degree y-spin and encode to a movie (GPU frames + ffmpeg)."""
    if not CHIMERAX:
        return {"status": "error", "error": "ChimeraX not found. Install it or set CHIMERAX_BIN."}
    if not os.path.isfile(structure):
        return {"status": "error", "error": f"Structure file not found: {structure}"}
    try:
        color = _validate_color(color)
    except ValueError as exc:
        return {"status": "error", "error": str(exc)}
    ffmpeg = shutil.which("ffmpeg")
    frame_dir = tempfile.mkdtemp(prefix="proteus_cxspin_")
    try:
        with ChimeraXRest() as rest:
            rest.run(f"open {_quote(structure)}")
            rest.run_all(_scene_commands(style, color))
            rest.run("view")  # centers content so `turn y` spins in place
            step = 360.0 / frames
            for idx in range(frames):
                rest.save_image(os.path.join(frame_dir, f"frame_{idx:04d}.png"),
                                width, height, supersample=1)
                rest.run(f"turn y {step:.4f}")
        if not ffmpeg:
            return {"status": "ok", "data": {
                "movie": None, "frames_dir": frame_dir, "frame_count": frames,
                "note": "ffmpeg not found; wrote frames but did not encode a movie.",
            }}
        _encode_movie(frame_dir, output, fps)
        shutil.rmtree(frame_dir, ignore_errors=True)
        return {"status": "ok", "data": {
            "movie": os.path.abspath(output), "frame_count": frames, "fps": fps,
        }}
    except subprocess.CalledProcessError as exc:
        return {"status": "error", "error": f"ffmpeg failed: {(exc.stderr or '')[:300]}",
                "data": {"frames_dir": frame_dir}}
    except Exception as exc:
        shutil.rmtree(frame_dir, ignore_errors=True)
        return {"status": "error", "error": str(exc)}


def main():
    parser = argparse.ArgumentParser(
        description="ChimeraX REST render agent — GPU rendering and turntable movies via a managed GUI session.",
        epilog="Examples:\n"
               "  %(prog)s render structure.pdb out.png --color plddt\n"
               "  %(prog)s spin model.pdb spin.mp4 --frames 72\n"
               "  %(prog)s run 'open 1ubq from pdb; cartoon; color bychain'",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="Run commands in a managed REST session")
    p_run.add_argument("commands", help="Semicolon-separated ChimeraX commands")

    p_render = sub.add_parser("render", help="Render a structure to PNG (GPU via REST)")
    p_render.add_argument("structure", help="Path to structure file")
    p_render.add_argument("output", nargs="?", default="/tmp/chimerax_render.png")
    p_render.add_argument("--style", default="cartoon", choices=["cartoon", "surface", "stick", "sphere"])
    p_render.add_argument("--color", default="rainbow", help="rainbow, bychain, bfactor, plddt, or a color name")
    p_render.add_argument("--width", type=int, default=DEFAULT_WIDTH)
    p_render.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    p_render.add_argument("--supersample", type=int, default=3)

    p_spin = sub.add_parser("spin", help="Render a 360-degree turntable movie (needs ffmpeg)")
    p_spin.add_argument("structure", help="Path to structure file")
    p_spin.add_argument("output", nargs="?", default="/tmp/chimerax_spin.mp4")
    p_spin.add_argument("--frames", type=int, default=72)
    p_spin.add_argument("--style", default="cartoon", choices=["cartoon", "surface", "stick", "sphere"])
    p_spin.add_argument("--color", default="rainbow", help="rainbow, bychain, bfactor, plddt, or a color name")
    p_spin.add_argument("--width", type=int, default=800)
    p_spin.add_argument("--height", type=int, default=600)
    p_spin.add_argument("--fps", type=int, default=30)

    args = parser.parse_args()

    if args.command == "run":
        result = rest_run(args.commands)
    elif args.command == "render":
        result = rest_render(args.structure, args.output, args.style, args.color,
                             args.width, args.height, args.supersample)
    elif args.command == "spin":
        result = rest_spin(args.structure, args.output, args.frames, args.style,
                           args.color, args.width, args.height, args.fps)
    else:
        parser.print_help()
        sys.exit(1)

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
