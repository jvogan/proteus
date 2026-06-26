#!/usr/bin/env python3
"""Optional model-quality wrappers for local structure-analysis tools.

Usage:
    python3 scripts/model_quality.py detect --json
    python3 scripts/model_quality.py usalign reference.pdb mobile.pdb --json
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
TOOL_CANDIDATES = {
    "usalign": ["USalign", "usalign"],
    "dockq": ["DockQ", "dockq", "DockQ.py"],
    "foldseek": ["foldseek"],
}


def _payload(status: str, data: dict) -> dict:
    output = {"status": status, "data": data}
    output.update(data)
    return output


def _find_candidate(candidates: list[str]) -> tuple[str | None, str | None]:
    for executable in candidates:
        path = shutil.which(executable)
        if path:
            return executable, path
    return None, None


def detect_tools() -> dict:
    tools = {}
    for name, candidates in TOOL_CANDIDATES.items():
        executable, path = _find_candidate(candidates)
        tools[name] = {
            "available": bool(path),
            "path": path,
            "executable": executable,
            "candidates": candidates,
        }
    return _payload(
        "ok",
        {
            "tools": tools,
            "capabilities": {name: tool["available"] for name, tool in tools.items()},
        },
    )


def parse_usalign_stdout(stdout: str) -> dict:
    """Extract core metrics from USalign/TM-align style stdout."""
    metrics = {
        "aligned_length": None,
        "rmsd": None,
        "seq_id": None,
        "tm_scores": [],
        "structures": {},
    }

    aligned_pattern = re.compile(
        rf"Aligned\s+length\s*=\s*(?P<aligned_length>\d+)\s*,\s*"
        rf"RMSD\s*=\s*(?P<rmsd>{FLOAT_RE})"
        rf"(?:\s*,\s*Seq_ID=n_identical/n_aligned\s*=\s*(?P<seq_id>{FLOAT_RE}))?",
        re.IGNORECASE,
    )
    aligned_match = aligned_pattern.search(stdout)
    if aligned_match:
        metrics["aligned_length"] = int(aligned_match.group("aligned_length"))
        metrics["rmsd"] = float(aligned_match.group("rmsd"))
        seq_id = aligned_match.group("seq_id")
        if seq_id is not None:
            metrics["seq_id"] = float(seq_id)

    name_pattern = re.compile(
        r"Name of (?:Chain|Structure)_(?P<index>[12])\s*:\s*(?P<name>.+)",
        re.IGNORECASE,
    )
    length_pattern = re.compile(
        r"Length of (?:Chain|Structure)_(?P<index>[12])\s*:\s*(?P<length>\d+)\s+residues",
        re.IGNORECASE,
    )
    for match in name_pattern.finditer(stdout):
        key = f"structure_{match.group('index')}"
        metrics["structures"].setdefault(key, {})["name"] = match.group("name").strip()
    for match in length_pattern.finditer(stdout):
        key = f"structure_{match.group('index')}"
        metrics["structures"].setdefault(key, {})["length"] = int(match.group("length"))

    tm_pattern = re.compile(
        rf"TM-score\s*=\s*(?P<score>{FLOAT_RE})\s*\((?P<description>[^)]*)\)",
        re.IGNORECASE,
    )
    normalization_pattern = re.compile(
        rf"(?:if\s+)?normalized by length of (?P<label>[^,:]+)"
        rf"(?::|,\s*i\.e\.,)\s*(?:L|LN)\s*=\s*(?P<length>\d+)\s*,\s*"
        rf"d0\s*=\s*(?P<d0>{FLOAT_RE})",
        re.IGNORECASE,
    )
    for match in tm_pattern.finditer(stdout):
        item = {
            "score": float(match.group("score")),
            "description": match.group("description").strip(),
        }
        normalization_match = normalization_pattern.search(item["description"])
        if normalization_match:
            item["normalized_by"] = normalization_match.group("label").strip()
            item["length"] = int(normalization_match.group("length"))
            item["d0"] = float(normalization_match.group("d0"))
        metrics["tm_scores"].append(item)

    return metrics


def _run_usalign(reference: str, mobile: str, timeout: int) -> dict:
    reference_path = Path(reference)
    mobile_path = Path(mobile)
    for path in (reference_path, mobile_path):
        if not path.exists():
            return _payload("error", {"error": f"File not found: {path}"})

    executable, usalign_path = _find_candidate(TOOL_CANDIDATES["usalign"])
    if not usalign_path:
        return _payload(
            "unavailable",
            {
                "tool": "usalign",
                "available": False,
                "error": "USalign executable not found on PATH",
                "reference": str(reference_path),
                "mobile": str(mobile_path),
            },
        )

    command = [usalign_path, str(reference_path), str(mobile_path)]
    try:
        proc = subprocess.run(
            command,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return _payload(
            "error",
            {
                "tool": "usalign",
                "available": True,
                "path": usalign_path,
                "executable": executable,
                "error": f"USalign timed out after {timeout} seconds",
                "stdout": exc.stdout or "",
                "stderr": exc.stderr or "",
            },
        )

    if proc.returncode != 0:
        return _payload(
            "error",
            {
                "tool": "usalign",
                "available": True,
                "path": usalign_path,
                "executable": executable,
                "returncode": proc.returncode,
                "error": "USalign exited with a non-zero status",
                "stdout": proc.stdout,
                "stderr": proc.stderr,
            },
        )

    metrics = parse_usalign_stdout(proc.stdout)
    return _payload(
        "ok",
        {
            "tool": "usalign",
            "available": True,
            "path": usalign_path,
            "executable": executable,
            "reference": str(reference_path),
            "mobile": str(mobile_path),
            "metrics": metrics,
        },
    )


def _print_detect_text(report: dict) -> None:
    for name, tool in report["tools"].items():
        value = tool["path"] if tool["available"] else "unavailable"
        print(f"{name}: {value}")


def _print_usalign_text(report: dict) -> None:
    if report["status"] != "ok":
        print(f"{report['status']}: {report.get('error', 'USalign failed')}")
        return
    metrics = report["metrics"]
    print(f"Aligned length: {metrics['aligned_length']}")
    print(f"RMSD: {metrics['rmsd']}")
    for item in metrics["tm_scores"]:
        label = item.get("normalized_by", item["description"])
        print(f"TM-score ({label}): {item['score']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Optional local model-quality tool wrappers.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    detect_parser = subparsers.add_parser("detect", help="Detect optional model-quality tools")
    detect_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    usalign_parser = subparsers.add_parser("usalign", help="Run USalign if it is available")
    usalign_parser.add_argument("reference", help="Reference structure path")
    usalign_parser.add_argument("mobile", help="Mobile structure path")
    usalign_parser.add_argument("--timeout", type=int, default=300, help="USalign timeout in seconds")
    usalign_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    args = parser.parse_args(argv)
    if args.command == "detect":
        report = detect_tools()
        if args.json:
            print(json.dumps(report, indent=2))
        else:
            _print_detect_text(report)
        return 0

    if args.command == "usalign":
        report = _run_usalign(args.reference, args.mobile, args.timeout)
        if args.json:
            print(json.dumps(report, indent=2))
        else:
            _print_usalign_text(report)
        return 1 if report["status"] == "error" else 0

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
