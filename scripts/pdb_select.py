#!/usr/bin/env python3
"""Rank candidate PDB entries and pick the best structure metadata match.

The offline path reads candidate metadata from JSON. The live path is opt-in
and uses the public RCSB search and entry APIs to fetch candidate metadata
before applying the same ranking logic.

Examples:
  python3 pdb_select.py --input candidates.json --json
  python3 pdb_select.py --input candidates.json --ligand ATP --limit 5 --json
  python3 pdb_select.py --live --query "KRAS G12C sotorasib" --uniprot P01116 --rows 20 --json
"""

import argparse
import json
import math
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
USER_AGENT = "proteus-skill/1.0"

METHOD_PREFERENCES = [
    ("X-ray", ("X-RAY", "XRAY", "X RAY")),
    ("electron microscopy", ("ELECTRON MICROSCOPY", "ELECTRON CRYSTALLOGRAPHY", "CRYO-EM", "CRYO EM")),
    ("NMR", ("NMR",)),
]

VALIDATION_FIELDS = {
    "clashscore": ("clashscore", "clash_score"),
    "ramachandran_outliers_percent": (
        "percent_ramachandran_outliers",
        "ramachandran_outliers_percent",
        "rama_outliers_percent",
    ),
    "rotamer_outliers_percent": (
        "percent_rotamer_outliers",
        "rotamer_outliers_percent",
    ),
    "bond_rmsz": ("bonds_RMSZ", "bond_rmsz"),
    "angle_rmsz": ("angles_RMSZ", "angle_rmsz"),
}


class PDBSelectError(RuntimeError):
    pass


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _get_path(value: dict, *path: str) -> Any:
    current: Any = value
    for key in path:
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def _as_list(value: Any) -> list:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def _numbers(value: Any) -> list[float]:
    if value is None or isinstance(value, bool):
        return []
    if isinstance(value, (int, float)):
        if math.isfinite(float(value)):
            return [float(value)]
        return []
    if isinstance(value, str):
        found = re.findall(r"[-+]?\d+(?:\.\d+)?", value)
        return [float(item) for item in found]
    if isinstance(value, dict):
        for key in ("value", "resolution", "high", "ls_d_res_high"):
            if key in value:
                return _numbers(value[key])
        return []
    if isinstance(value, (list, tuple)):
        items: list[float] = []
        for item in value:
            items.extend(_numbers(item))
        return items
    return []


def _first_number(*values: Any) -> float | None:
    numbers: list[float] = []
    for value in values:
        numbers.extend(number for number in _numbers(value) if number > 0)
    if not numbers:
        return None
    return min(numbers)


def _first_nonnegative_number(*values: Any) -> float | None:
    numbers: list[float] = []
    for value in values:
        numbers.extend(number for number in _numbers(value) if number >= 0)
    if not numbers:
        return None
    return min(numbers)


def _first_dict(value: Any) -> dict:
    if isinstance(value, dict):
        return value
    if isinstance(value, list):
        for item in value:
            if isinstance(item, dict):
                return item
    return {}


def _unique_strings(values: list[Any]) -> list[str]:
    seen = set()
    result = []
    for value in values:
        if value is None:
            continue
        text = str(value).strip()
        if not text:
            continue
        if text not in seen:
            seen.add(text)
            result.append(text)
    return result


def normalize_ligand_filters(values: list[str] | None) -> list[str]:
    filters = []
    seen = set()
    for value in values or []:
        for item in value.split(","):
            code = item.strip().upper()
            if not code:
                continue
            if not re.fullmatch(r"[A-Z0-9]{1,10}", code):
                raise PDBSelectError(f"Invalid ligand code '{item.strip()}'.")
            if code not in seen:
                seen.add(code)
                filters.append(code)
    return filters


def _ligand_codes(value: Any) -> list[str]:
    codes = []
    for item in _as_list(value):
        if isinstance(item, dict):
            for key in ("ligand", "comp_id", "chem_comp_id", "id", "code"):
                if item.get(key):
                    codes.append(item[key])
                    break
        elif isinstance(item, str):
            codes.extend(part.strip() for part in item.split(","))
        elif item is not None:
            codes.append(item)
    normalized = []
    seen = set()
    for code in codes:
        text = str(code).strip().upper()
        if not text:
            continue
        if text not in seen:
            seen.add(text)
            normalized.append(text)
    return normalized


def _extract_ligands(candidate: dict) -> list[str]:
    values = [
        candidate.get("ligands"),
        candidate.get("ligand_codes"),
        candidate.get("nonpolymer_bound_components"),
        candidate.get("ligand_components"),
        _get_path(candidate, "rcsb_entry_info", "nonpolymer_bound_components"),
    ]
    codes = []
    for value in values:
        codes.extend(_ligand_codes(value))
    return sorted(set(codes))


def _extract_methods(candidate: dict) -> list[str]:
    values: list[Any] = [
        candidate.get("method"),
        candidate.get("experimental_method"),
        candidate.get("experimental_methods"),
        _get_path(candidate, "rcsb_entry_info", "experimental_method"),
    ]
    for item in _as_list(candidate.get("exptl")):
        if isinstance(item, dict):
            values.append(item.get("method"))
        else:
            values.append(item)
    methods: list[Any] = []
    for value in values:
        methods.extend(_as_list(value))
    return _unique_strings(methods)


def _method_preference(methods: list[str]) -> tuple[str, int]:
    if not methods:
        return ("unknown", len(METHOD_PREFERENCES) + 2)
    joined = " | ".join(methods).upper()
    for index, (label, markers) in enumerate(METHOD_PREFERENCES, start=1):
        if any(marker in joined for marker in markers):
            return (label, index)
    return ("other", len(METHOD_PREFERENCES) + 1)


def _extract_pdb_id(candidate: dict, index: int) -> str:
    values = [
        candidate.get("pdb_id"),
        candidate.get("id"),
        candidate.get("identifier"),
        candidate.get("entry_id"),
        candidate.get("rcsb_id"),
        _get_path(candidate, "rcsb_entry_container_identifiers", "entry_id"),
        _get_path(candidate, "entry", "id"),
    ]
    for value in values:
        if value:
            return str(value).strip().upper()
    return f"CANDIDATE_{index + 1}"


def _extract_title(candidate: dict) -> str | None:
    value = candidate.get("title") or _get_path(candidate, "struct", "title")
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _extract_resolution(candidate: dict) -> float | None:
    refine_resolution = None
    for item in _as_list(candidate.get("refine")):
        if isinstance(item, dict):
            refine_resolution = _first_number(refine_resolution, item.get("ls_d_res_high"))
    return _first_number(
        candidate.get("resolution"),
        candidate.get("resolution_combined"),
        _get_path(candidate, "rcsb_entry_info", "resolution_combined"),
        _get_path(candidate, "diffrn_resolution_high", "value"),
        refine_resolution,
    )


def _extract_assembly_count(candidate: dict) -> int | None:
    value = candidate.get("assembly_count")
    if value is None:
        value = _get_path(candidate, "rcsb_entry_info", "assembly_count")
    numbers = _numbers(value)
    if not numbers:
        return None
    return int(numbers[0])


def _validation_sources(candidate: dict) -> list[dict]:
    sources = []
    validation = candidate.get("validation")
    if isinstance(validation, dict):
        sources.append(validation)
        if isinstance(validation.get("geometry"), dict):
            sources.append(validation["geometry"])
    geometry = candidate.get("geometry")
    if isinstance(geometry, dict):
        sources.append(geometry)
    sources.append(_first_dict(candidate.get("pdbx_vrpt_summary_geometry")))
    return [source for source in sources if source]


def _extract_validation(candidate: dict) -> dict:
    metrics = {}
    for source in _validation_sources(candidate):
        for output_key, source_keys in VALIDATION_FIELDS.items():
            if output_key in metrics:
                continue
            value = None
            for source_key in source_keys:
                if source_key in source:
                    value = source[source_key]
                    break
            number = _first_nonnegative_number(value)
            if number is not None:
                metrics[output_key] = number
    return metrics


def normalize_candidate(candidate: dict, index: int = 0) -> dict:
    if not isinstance(candidate, dict):
        raise PDBSelectError("Candidate entries must be JSON objects.")
    methods = _extract_methods(candidate)
    method_category, method_rank = _method_preference(methods)
    normalized = {
        "pdb_id": _extract_pdb_id(candidate, index),
        "title": _extract_title(candidate),
        "resolution": _extract_resolution(candidate),
        "experimental_methods": methods,
        "method_category": method_category,
        "method_preference_rank": method_rank,
        "ligands": _extract_ligands(candidate),
        "assembly_count": _extract_assembly_count(candidate),
        "validation": _extract_validation(candidate),
    }
    if candidate.get("fetch_error"):
        normalized["fetch_error"] = str(candidate["fetch_error"])
    return normalized


def unwrap_candidates(payload: Any) -> list[dict]:
    if isinstance(payload, list):
        return payload
    if not isinstance(payload, dict):
        raise PDBSelectError("Candidate JSON must be an object or list.")

    if "data" in payload:
        try:
            return unwrap_candidates(payload["data"])
        except PDBSelectError:
            pass

    for key in ("candidates", "results", "entries", "items", "ranked", "result_set"):
        value = payload.get(key)
        if isinstance(value, list):
            return value

    if "rcsb_entry_info" in payload or "rcsb_id" in payload or "pdb_id" in payload or "id" in payload:
        return [payload]

    raise PDBSelectError("Could not find candidates in JSON; expected candidates/results/entries/items.")


def load_candidates(path: Path) -> list[dict]:
    try:
        with path.open() as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        raise PDBSelectError(f"Invalid JSON: {exc}") from exc
    except OSError as exc:
        raise PDBSelectError(f"Could not read input JSON: {exc}") from exc
    return unwrap_candidates(payload)


def _validation_penalty(metrics: dict) -> float:
    if not metrics:
        return 0.0
    score = 0.0
    score += metrics.get("clashscore", 0.0) / 100.0
    score += metrics.get("ramachandran_outliers_percent", 0.0) / 10.0
    score += metrics.get("rotamer_outliers_percent", 0.0) / 10.0
    score += metrics.get("bond_rmsz", 0.0) / 10.0
    score += metrics.get("angle_rmsz", 0.0) / 10.0
    return round(score, 6)


def _sort_key(candidate: dict) -> tuple:
    resolution = candidate["resolution"]
    return (
        1 if resolution is None else 0,
        resolution if resolution is not None else 9999.0,
        candidate["method_preference_rank"],
        1 if candidate["assembly_count"] is None else 0,
        1 if not candidate["validation"] else 0,
        _validation_penalty(candidate["validation"]),
        candidate["pdb_id"],
    )


def _selection_key(candidate: dict) -> dict:
    resolution = candidate["resolution"]
    return {
        "resolution_missing": resolution is None,
        "resolution": resolution,
        "method_preference_rank": candidate["method_preference_rank"],
        "assembly_count_missing": candidate["assembly_count"] is None,
        "validation_metrics_missing": not bool(candidate["validation"]),
        "validation_penalty": _validation_penalty(candidate["validation"]),
    }


def _candidate_summary(candidate: dict) -> dict:
    return {
        "pdb_id": candidate["pdb_id"],
        "title": candidate["title"],
        "resolution": candidate["resolution"],
        "experimental_methods": candidate["experimental_methods"],
        "method_category": candidate["method_category"],
        "ligands": candidate["ligands"],
        "assembly_count": candidate["assembly_count"],
        "validation": candidate["validation"],
        **({"fetch_error": candidate["fetch_error"]} if candidate.get("fetch_error") else {}),
    }


def _candidate_reasons(candidate: dict, ligand_filters: list[str]) -> list[str]:
    reasons = []
    resolution = candidate["resolution"]
    if resolution is None:
        reasons.append("resolution unavailable; ranked after candidates with known resolution")
    else:
        reasons.append(f"resolution {resolution:g} A; lower is better")

    methods = ", ".join(candidate["experimental_methods"]) or "unavailable"
    category = candidate["method_category"]
    if category in {"X-ray", "electron microscopy", "NMR"}:
        reasons.append(
            f"experimental method {methods} matched preference tier {candidate['method_preference_rank']} ({category})"
        )
    else:
        reasons.append(f"experimental method {methods} has preference tier {candidate['method_preference_rank']}")

    if ligand_filters:
        matched = sorted(set(candidate["ligands"]) & set(ligand_filters))
        reasons.append(f"matched ligand filter: {', '.join(matched)}")
    elif candidate["ligands"]:
        reasons.append(f"bound ligands reported: {', '.join(candidate['ligands'])}")
    else:
        reasons.append("no bound ligands reported")

    if candidate["assembly_count"] is None:
        reasons.append("assembly_count missing")
    else:
        reasons.append(f"assembly_count present: {candidate['assembly_count']}")

    if candidate["validation"]:
        fields = ", ".join(sorted(candidate["validation"]))
        reasons.append(f"optional validation metrics present: {fields}")
    else:
        reasons.append("optional validation metrics not provided")

    if candidate.get("fetch_error"):
        reasons.append(f"entry detail fetch issue: {candidate['fetch_error']}")

    return reasons


def rank_candidates(candidates: list[dict], ligand_filters: list[str] | None = None,
                    limit: int | None = None) -> dict:
    filters = ligand_filters or []
    normalized = [normalize_candidate(candidate, index) for index, candidate in enumerate(candidates)]

    eligible = []
    excluded = []
    filter_set = set(filters)
    for candidate in normalized:
        if filter_set and not filter_set.intersection(candidate["ligands"]):
            summary = _candidate_summary(candidate)
            summary["exclude_reasons"] = [
                f"missing requested ligand(s): {', '.join(filters)}"
            ]
            excluded.append(summary)
            continue
        eligible.append(candidate)

    ranked_candidates = sorted(eligible, key=_sort_key)
    limited_candidates = ranked_candidates[:limit] if limit is not None else ranked_candidates
    ranked = []
    for rank, candidate in enumerate(limited_candidates, start=1):
        item = _candidate_summary(candidate)
        item["rank"] = rank
        item["selection_key"] = _selection_key(candidate)
        item["reasons"] = _candidate_reasons(candidate, filters)
        ranked.append(item)

    return {
        "input_count": len(candidates),
        "eligible_count": len(eligible),
        "ranked_count": len(ranked),
        "excluded_count": len(excluded),
        "limit": limit,
        "ligand_filter": filters,
        "ranking_policy": {
            "sort_order": [
                "known resolution before unknown resolution",
                "lower resolution in Angstroms",
                "experimental method preference: X-ray, electron microscopy, NMR, other, unknown",
                "assembly_count present before missing",
                "optional validation metrics present before missing",
                "lower optional validation penalty when metrics are present",
                "PDB ID as deterministic tie-breaker",
            ],
        },
        "best": ranked[0] if ranked else None,
        "ranked": ranked,
        "excluded": excluded,
    }


def _full_text_node(value: str) -> dict:
    return {"type": "terminal", "service": "full_text", "parameters": {"value": value}}


def _uniprot_nodes(accession: str) -> list[dict]:
    base = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers"
    return [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": f"{base}.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": f"{base}.database_name",
                "operator": "exact_match",
                "value": "UniProt",
            },
        },
    ]


def build_search_query(text: str | None, uniprot: str | None) -> dict:
    nodes = []
    if text:
        nodes.append(_full_text_node(text))
    if uniprot:
        nodes.extend(_uniprot_nodes(uniprot.strip().upper()))
    if not nodes:
        raise PDBSelectError("Live search requires --query and/or --uniprot.")
    if len(nodes) == 1:
        return nodes[0]
    return {"type": "group", "logical_operator": "and", "nodes": nodes}


def _fetch_json(url: str, payload: dict | None = None, timeout: int = 30) -> dict:
    data = json.dumps(payload).encode() if payload is not None else None
    headers = {"User-Agent": USER_AGENT}
    if payload is not None:
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        raise PDBSelectError(f"RCSB API returned HTTP {exc.code}: {exc.reason}") from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise PDBSelectError(f"RCSB API request failed: {exc}") from exc


def search_rcsb_ids(text: str | None, uniprot: str | None, rows: int) -> list[str]:
    payload = {
        "query": build_search_query(text, uniprot),
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_content_type": ["experimental"],
            "scoring_strategy": "combined",
        },
    }
    result = _fetch_json(SEARCH_URL, payload)
    return [
        str(hit["identifier"]).upper()
        for hit in result.get("result_set", [])
        if hit.get("identifier")
    ]


def fetch_live_candidates(text: str | None, uniprot: str | None, rows: int) -> list[dict]:
    ids = search_rcsb_ids(text, uniprot, rows)
    candidates = []
    for pdb_id in ids:
        try:
            candidates.append(_fetch_json(ENTRY_URL.format(pdb_id=pdb_id)))
        except PDBSelectError as exc:
            candidates.append({"id": pdb_id, "fetch_error": str(exc)})
    return candidates


def build_selection_report(candidates: list[dict], source: dict, ligand_filters: list[str] | None = None,
                           limit: int | None = None, query: str | None = None,
                           uniprot: str | None = None) -> dict:
    data = rank_candidates(candidates, ligand_filters, limit)
    data["source"] = source
    data["query"] = query
    data["uniprot"] = uniprot
    return _ok_payload(data)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Rank candidate PDB entries by resolution, method preference, ligands, assembly metadata, and validation hints.",
        epilog=(
            "Examples:\n"
            "  %(prog)s --input candidates.json --json\n"
            "  %(prog)s --input candidates.json --ligand ATP --limit 5 --json\n"
            "  %(prog)s --live --query \"KRAS G12C sotorasib\" --uniprot P01116 --rows 20 --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, help="Offline candidate JSON file to rank")
    parser.add_argument("--live", action="store_true", help="Fetch candidate metadata from the public RCSB APIs")
    parser.add_argument("--query", help="Free-text RCSB query for --live")
    parser.add_argument("--uniprot", help="UniProt accession filter for --live, e.g. P01116")
    parser.add_argument("--rows", type=int, default=25, help="Maximum RCSB search hits for --live (default: 25)")
    parser.add_argument("--ligand", action="append", help="Require ligand code; repeat or comma-separate values")
    parser.add_argument("--limit", type=int, help="Maximum ranked candidates to return")
    parser.add_argument("--json", action="store_true", help="Accepted for consistency; output is always JSON")
    args = parser.parse_args(argv)

    try:
        if args.input and args.live:
            raise PDBSelectError("Use either --input or --live, not both.")
        if not args.input and not args.live:
            raise PDBSelectError("Provide --input for offline ranking, or --live with --query and/or --uniprot.")
        if args.rows < 1:
            raise PDBSelectError("--rows must be a positive integer.")
        if args.limit is not None and args.limit < 1:
            raise PDBSelectError("--limit must be a positive integer.")

        ligand_filters = normalize_ligand_filters(args.ligand)
        if args.live:
            candidates = fetch_live_candidates(args.query, args.uniprot, args.rows)
            source = {"kind": "rcsb_live", "rows": args.rows}
        else:
            candidates = load_candidates(args.input)
            source = {"kind": "input_file", "file": args.input.name}

        output = build_selection_report(
            candidates,
            source,
            ligand_filters=ligand_filters,
            limit=args.limit,
            query=args.query,
            uniprot=args.uniprot,
        )
    except PDBSelectError as exc:
        print(json.dumps(_error_payload(str(exc)), indent=2))
        return 1

    print(json.dumps(output, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
