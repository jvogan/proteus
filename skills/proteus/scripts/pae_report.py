#!/usr/bin/env python3
"""Summarize AlphaFold predicted aligned error (PAE) JSON.

Usage:
    python3 pae_report.py AF-P04637-F1_pae.json --json
    python3 pae_report.py pae.json --flex-threshold 15
"""

import argparse
import json
import sys
from pathlib import Path


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _matrix_from_payload(payload) -> list[list[float]]:
    if isinstance(payload, dict):
        for key in ("predicted_aligned_error", "pae"):
            if key in payload:
                return payload[key]
        if "pae_matrix" in payload:
            return payload["pae_matrix"]
    if isinstance(payload, list):
        if payload and isinstance(payload[0], dict):
            for key in ("predicted_aligned_error", "pae"):
                if key in payload[0]:
                    return payload[0][key]
        if payload and isinstance(payload[0], list):
            return payload
    raise ValueError("Could not find a PAE matrix in the JSON payload.")


def _validate_matrix(matrix) -> list[list[float]]:
    if not isinstance(matrix, list) or not matrix:
        raise ValueError("PAE matrix is empty.")
    width = len(matrix[0])
    if width == 0:
        raise ValueError("PAE matrix rows are empty.")
    normalized = []
    for row in matrix:
        if not isinstance(row, list) or len(row) != width:
            raise ValueError("PAE matrix must be rectangular.")
        normalized.append([float(value) for value in row])
    return normalized


def _segments(values: list[float], threshold: float, min_len: int, high: bool) -> list[dict]:
    segments = []
    start = None
    for index, value in enumerate(values, start=1):
        match = value >= threshold if high else value <= threshold
        if match and start is None:
            start = index
        elif not match and start is not None:
            end = index - 1
            if end - start + 1 >= min_len:
                window = values[start - 1:end]
                segments.append({
                    "start": start,
                    "end": end,
                    "length": end - start + 1,
                    "mean_pae": round(sum(window) / len(window), 2),
                })
            start = None
    if start is not None:
        end = len(values)
        if end - start + 1 >= min_len:
            window = values[start - 1:end]
            segments.append({
                "start": start,
                "end": end,
                "length": end - start + 1,
                "mean_pae": round(sum(window) / len(window), 2),
            })
    return segments


def summarize_pae(path: str, flex_threshold: float = 15.0, rigid_threshold: float = 5.0,
                  min_segment: int = 5) -> dict:
    payload = json.loads(Path(path).read_text())
    matrix = _validate_matrix(_matrix_from_payload(payload))
    row_count = len(matrix)
    col_count = len(matrix[0])
    flat = [value for row in matrix for value in row]
    per_residue_mean = [sum(row) / len(row) for row in matrix]
    data = {
        "file": str(Path(path).resolve()),
        "size": {"rows": row_count, "columns": col_count},
        "pae": {
            "min": round(min(flat), 2),
            "max": round(max(flat), 2),
            "mean": round(sum(flat) / len(flat), 2),
        },
        "thresholds": {
            "flexible_or_uncertain_mean_pae": flex_threshold,
            "rigid_domain_mean_pae": rigid_threshold,
            "min_segment": min_segment,
        },
        "flexible_or_uncertain_regions": _segments(per_residue_mean, flex_threshold, min_segment, high=True),
        "rigid_candidate_regions": _segments(per_residue_mean, rigid_threshold, min_segment, high=False),
        "per_residue_mean_pae": [round(value, 2) for value in per_residue_mean],
    }
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(description="Summarize AlphaFold PAE JSON into domain/flexibility hints.")
    parser.add_argument("pae_json", help="Path to PAE JSON")
    parser.add_argument("--flex-threshold", type=float, default=15.0,
                        help="Mean PAE threshold for uncertain regions (default: 15)")
    parser.add_argument("--rigid-threshold", type=float, default=5.0,
                        help="Mean PAE threshold for rigid candidate regions (default: 5)")
    parser.add_argument("--min-segment", type=int, default=5,
                        help="Minimum contiguous segment length (default: 5)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        output = summarize_pae(args.pae_json, args.flex_threshold, args.rigid_threshold, args.min_segment)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        if args.json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        data = output["data"]
        print(f"File: {data['file']}")
        print(f"Matrix: {data['size']['rows']} x {data['size']['columns']}")
        print(f"PAE mean: {data['pae']['mean']} A (min {data['pae']['min']}, max {data['pae']['max']})")
        print(f"Flexible/uncertain regions: {data['flexible_or_uncertain_regions'] or 'none'}")
        print(f"Rigid candidate regions: {data['rigid_candidate_regions'] or 'none'}")


if __name__ == "__main__":
    main()
