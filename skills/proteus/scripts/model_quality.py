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
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common


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
            "path": proteus_common.display_path(path) if path else None,
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
            return _payload("error", {"error": f"File not found: {proteus_common.display_path(path)}"})

    executable, usalign_path = _find_candidate(TOOL_CANDIDATES["usalign"])
    if not usalign_path:
        return _payload(
            "unavailable",
            {
                "tool": "usalign",
                "available": False,
                "error": "USalign executable not found on PATH",
                "reference": proteus_common.display_path(reference_path),
                "mobile": proteus_common.display_path(mobile_path),
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
                "path": proteus_common.display_path(usalign_path),
                "executable": executable,
                "error": f"USalign timed out after {timeout} seconds",
                "stdout": proteus_common.scrub_text(exc.stdout or ""),
                "stderr": proteus_common.scrub_text(exc.stderr or ""),
            },
        )

    if proc.returncode != 0:
        return _payload(
            "error",
            {
                "tool": "usalign",
                "available": True,
                "path": proteus_common.display_path(usalign_path),
                "executable": executable,
                "returncode": proc.returncode,
                "error": "USalign exited with a non-zero status",
                "stdout": proteus_common.scrub_text(proc.stdout),
                "stderr": proteus_common.scrub_text(proc.stderr),
            },
        )

    metrics = parse_usalign_stdout(proc.stdout)
    return _payload(
        "ok",
        {
            "tool": "usalign",
            "available": True,
            "path": proteus_common.display_path(usalign_path),
            "executable": executable,
            "reference": proteus_common.display_path(reference_path),
            "mobile": proteus_common.display_path(mobile_path),
            "metrics": metrics,
        },
    )


def _validated_paths(*values: str) -> tuple[list[Path] | None, dict | None]:
    paths = [Path(value) for value in values]
    missing = [proteus_common.display_path(path) for path in paths if not path.exists()]
    if missing:
        return None, _payload("error", {"error": f"File or directory not found: {missing[0]}"})
    return paths, None


def _run_dockq(model: str, native: str, timeout: int) -> dict:
    paths, error = _validated_paths(model, native)
    if error:
        return error
    assert paths is not None
    executable, dockq_path = _find_candidate(TOOL_CANDIDATES["dockq"])
    if not dockq_path:
        return _payload("unavailable", {"tool": "dockq", "available": False, "error": "DockQ executable not found on PATH"})
    with tempfile.NamedTemporaryFile(suffix=".json") as output:
        try:
            proc = subprocess.run(
                [dockq_path, str(paths[0]), str(paths[1]), "--json", output.name],
                text=True, capture_output=True, timeout=timeout, check=False,
            )
        except subprocess.TimeoutExpired:
            return _payload("error", {"tool": "dockq", "available": True, "error": f"DockQ timed out after {timeout} seconds"})
        if proc.returncode != 0:
            return _payload("error", {
                "tool": "dockq", "available": True, "returncode": proc.returncode,
                "error": "DockQ exited with a non-zero status", "stderr": proteus_common.scrub_text(proc.stderr),
            })
        try:
            metrics = json.loads(Path(output.name).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            return _payload("error", {"tool": "dockq", "available": True, "error": f"DockQ JSON could not be read: {exc}"})
    return _payload("ok", {
        "tool": "dockq", "available": True, "path": proteus_common.display_path(dockq_path), "executable": executable,
        "model": proteus_common.display_path(paths[0]), "native": proteus_common.display_path(paths[1]), "metrics": metrics,
        "interpretation": "DockQ is interface- and chain-mapping-dependent; review the reported mapping before using the score.",
    })


FOLDSEEK_FIELDS = ["query", "target", "alntmscore", "qtmscore", "ttmscore", "lddt", "evalue", "bits"]


def _parse_foldseek_table(text: str) -> list[dict]:
    results = []
    numeric = set(FOLDSEEK_FIELDS[2:])
    for line in text.splitlines():
        columns = line.split("\t")
        if len(columns) != len(FOLDSEEK_FIELDS):
            continue
        item = dict(zip(FOLDSEEK_FIELDS, columns))
        for key in numeric:
            try:
                item[key] = float(item[key])
            except ValueError:
                pass
        results.append(item)
    return results


def _run_foldseek(query: str, target: str, timeout: int) -> dict:
    paths, error = _validated_paths(query, target)
    if error:
        return error
    assert paths is not None
    executable, foldseek_path = _find_candidate(TOOL_CANDIDATES["foldseek"])
    if not foldseek_path:
        return _payload("unavailable", {"tool": "foldseek", "available": False, "error": "Foldseek executable not found on PATH"})
    with tempfile.TemporaryDirectory(prefix="proteus-foldseek-") as temporary:
        root = Path(temporary)
        output = root / "results.tsv"
        work = root / "work"
        command = [
            foldseek_path, "easy-search", str(paths[0]), str(paths[1]), str(output), str(work),
            "--format-output", ",".join(FOLDSEEK_FIELDS),
        ]
        try:
            proc = subprocess.run(command, text=True, capture_output=True, timeout=timeout, check=False)
        except subprocess.TimeoutExpired:
            return _payload("error", {"tool": "foldseek", "available": True, "error": f"Foldseek timed out after {timeout} seconds"})
        if proc.returncode != 0:
            return _payload("error", {
                "tool": "foldseek", "available": True, "returncode": proc.returncode,
                "error": "Foldseek exited with a non-zero status", "stderr": proteus_common.scrub_text(proc.stderr),
            })
        if not output.is_file():
            return _payload("error", {"tool": "foldseek", "available": True, "error": "Foldseek did not create its result table"})
        hits = _parse_foldseek_table(output.read_text(encoding="utf-8", errors="replace"))
    return _payload("ok", {
        "tool": "foldseek", "available": True, "path": proteus_common.display_path(foldseek_path), "executable": executable,
        "query": proteus_common.display_path(paths[0]), "target": proteus_common.display_path(paths[1]), "hit_count": len(hits), "hits": hits,
        "interpretation": "TM-score normalization and local/global alignment choices affect comparisons; inspect coverage and both normalized scores.",
    })


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

    dockq_parser = subparsers.add_parser("dockq", help="Run DockQ on a model/native complex pair")
    dockq_parser.add_argument("model", help="Docked model structure")
    dockq_parser.add_argument("native", help="Native/reference complex")
    dockq_parser.add_argument("--timeout", type=int, default=600)
    dockq_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    foldseek_parser = subparsers.add_parser("foldseek", help="Run a local Foldseek easy-search")
    foldseek_parser.add_argument("query", help="Query structure file or directory")
    foldseek_parser.add_argument("target", help="Target structure file, directory, or database")
    foldseek_parser.add_argument("--timeout", type=int, default=1200)
    foldseek_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

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

    if args.command == "dockq":
        report = _run_dockq(args.model, args.native, args.timeout)
        print(json.dumps(report, indent=2))
        return 1 if report["status"] == "error" else 0

    if args.command == "foldseek":
        report = _run_foldseek(args.query, args.target, args.timeout)
        print(json.dumps(report, indent=2))
        return 1 if report["status"] == "error" else 0

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
