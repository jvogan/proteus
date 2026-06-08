#!/usr/bin/env python3
"""Inspect an MRC/CCP4 density map and suggest contour levels.

Cryo-EM contour levels are in absolute map units that differ per map, so a
hard-coded `--level` rarely transfers. This reads the MRC/CCP4 header with the
standard library and derives sigma-based levels (mean + N*sigma) from the voxel
data — the protein-agnostic default ChimeraX `volume` / PyMOL `isomesh` want.

Uses numpy if available for speed; otherwise a pure-stdlib strided sampler keeps
it dependency-free.

Usage:
    python map_info.py map.mrc --json
    python map_info.py map.mrc --n-sigma 1.5
    python map_info.py --help
"""

import argparse
import json
import os
import struct
import sys


_MODE_BYTES = {0: 1, 1: 2, 2: 4, 6: 2}
_MODE_STRUCT = {0: "b", 1: "h", 2: "f", 6: "H"}
_MODE_NAME = {0: "int8", 1: "int16", 2: "float32", 6: "uint16"}


def read_map_stats(path, sample_target=200000):
    """Return (dims, mode, mean, sigma) for an MRC/CCP4 map.

    The MRC header is 1024 bytes: nx/ny/nz/mode are the first four int32s, and
    the extended-header byte count (nsymbt) is at offset 92. Voxel data starts
    at 1024 + nsymbt.
    """
    with open(path, "rb") as fh:
        header = fh.read(1024)
        if len(header) < 1024:
            raise ValueError("File too small to be an MRC/CCP4 map")
        nx, ny, nz, mode = struct.unpack("<4i", header[:16])
        nsymbt = struct.unpack("<i", header[92:96])[0]
        fh.seek(1024 + nsymbt)
        raw = fh.read()
    if mode not in _MODE_BYTES:
        raise ValueError(f"Unsupported MRC mode {mode}")
    size = _MODE_BYTES[mode]
    n = len(raw) // size
    if n == 0:
        raise ValueError("Map contains no voxel data")
    try:
        import numpy as np
        dt = {0: np.int8, 1: np.int16, 2: np.float32, 6: np.uint16}[mode]
        data = np.frombuffer(raw[: n * size], dtype=dt).astype(np.float64)
        mean, sigma = float(data.mean()), float(data.std())
    except Exception:
        fmt = "<" + _MODE_STRUCT[mode]
        step = max(1, n // sample_target)
        vals = [struct.unpack_from(fmt, raw, i * size)[0] for i in range(0, n, step)]
        mean = sum(vals) / len(vals)
        sigma = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
    return (nx, ny, nz), mode, mean, sigma


def map_info(path, n_sigma=2.0):
    dims, mode, mean, sigma = read_map_stats(path)
    return {"status": "ok", "data": {
        "map": os.path.abspath(path),
        "dimensions": {"nx": dims[0], "ny": dims[1], "nz": dims[2]},
        "mode": mode,
        "dtype": _MODE_NAME.get(mode),
        "mean": round(mean, 6),
        "sigma": round(sigma, 6),
        "n_sigma": n_sigma,
        "suggested_level": round(mean + n_sigma * sigma, 6),
        "suggested_levels": {f"{s}_sigma": round(mean + s * sigma, 6)
                             for s in (1.0, 1.5, 2.0, 3.0)},
    }}


def main():
    parser = argparse.ArgumentParser(
        description="Inspect an MRC/CCP4 density map and suggest sigma-based contour levels.",
    )
    parser.add_argument("map", help="Path to .mrc/.map/.ccp4 file")
    parser.add_argument("--n-sigma", type=float, default=2.0,
                        help="Sigma multiple for the primary suggested level (default: 2.0)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    if not os.path.isfile(args.map):
        result = {"status": "error", "error": f"File not found: {args.map}"}
    else:
        try:
            result = map_info(args.map, args.n_sigma)
        except (ValueError, struct.error, OSError) as exc:
            result = {"status": "error", "error": str(exc)}

    if args.json:
        print(json.dumps(result, indent=2))
    elif result["status"] == "ok":
        d = result["data"]
        print(d["map"])
        print(f"  dims: {d['dimensions']['nx']}x{d['dimensions']['ny']}x{d['dimensions']['nz']}"
              f"  dtype: {d['dtype']}")
        print(f"  mean: {d['mean']}  sigma: {d['sigma']}")
        print(f"  suggested level ({d['n_sigma']} sigma): {d['suggested_level']}")
    else:
        print(result["error"], file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
