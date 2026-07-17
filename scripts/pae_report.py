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

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common


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


def _parse_chain_lengths(value: str | None) -> list[tuple[str, int]]:
    if not value:
        return []
    output = []
    for index, item in enumerate(value.split(","), start=1):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            chain, length_text = item.rsplit(":", 1)
        else:
            chain, length_text = f"chain_{index}", item
        try:
            length = int(length_text)
        except ValueError as exc:
            raise ValueError("Chain lengths must look like A:120,B:95 or 120,95.") from exc
        if not chain or length <= 0:
            raise ValueError("Chain names must be non-empty and lengths must be positive.")
        output.append((chain, length))
    return output


def _chain_blocks(matrix: list[list[float]], chain_lengths: list[tuple[str, int]]) -> dict:
    if not chain_lengths:
        return {}
    size = len(matrix)
    if len(matrix[0]) != size:
        raise ValueError("Chain-aware PAE analysis requires a square matrix.")
    if sum(length for _chain, length in chain_lengths) != size:
        raise ValueError("Chain lengths must sum to the PAE matrix size.")
    ranges = []
    start = 0
    for chain, length in chain_lengths:
        ranges.append((chain, start, start + length))
        start += length
    blocks = []
    for row_chain, row_start, row_end in ranges:
        for column_chain, column_start, column_end in ranges:
            values = [matrix[row][column] for row in range(row_start, row_end) for column in range(column_start, column_end)]
            blocks.append({
                "row_chain": row_chain,
                "column_chain": column_chain,
                "kind": "within_chain" if row_chain == column_chain else "between_chain",
                "mean_pae": round(sum(values) / len(values), 2),
                "min_pae": round(min(values), 2),
                "max_pae": round(max(values), 2),
            })
    residue_map = []
    for chain, chain_start, chain_end in ranges:
        for index in range(chain_start, chain_end):
            residue_map.append({"index": index + 1, "chain": chain, "chain_index": index - chain_start + 1})
    return {"chains": [{"chain": chain, "length": end - start} for chain, start, end in ranges], "blocks": blocks, "residue_map": residue_map}


def summarize_pae(path: str, flex_threshold: float = 15.0, rigid_threshold: float = 5.0,
                  min_segment: int = 5, chain_lengths: str | None = None) -> dict:
    payload = json.loads(Path(path).read_text())
    matrix = _validate_matrix(_matrix_from_payload(payload))
    row_count = len(matrix)
    col_count = len(matrix[0])
    flat = [value for row in matrix for value in row]
    per_residue_mean = [sum(row) / len(row) for row in matrix]
    data = {
        "file": proteus_common.display_path(path),
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
    parsed_chains = _parse_chain_lengths(chain_lengths)
    if parsed_chains:
        data["chain_analysis"] = _chain_blocks(matrix, parsed_chains)
    return _ok_payload(data)


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description="Summarize AlphaFold PAE JSON into domain/flexibility hints.")
    parser.add_argument("pae_json", help="Path to PAE JSON")
    parser.add_argument("--flex-threshold", type=float, default=15.0,
                        help="Mean PAE threshold for uncertain regions (default: 15)")
    parser.add_argument("--rigid-threshold", type=float, default=5.0,
                        help="Mean PAE threshold for rigid candidate regions (default: 5)")
    parser.add_argument("--min-segment", type=int, default=5,
                        help="Minimum contiguous segment length (default: 5)")
    parser.add_argument("--chain-lengths", help="Chain names and lengths, e.g. A:120,B:95")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args(argv)

    try:
        output = summarize_pae(args.pae_json, args.flex_threshold, args.rigid_threshold,
                               args.min_segment, args.chain_lengths)
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
