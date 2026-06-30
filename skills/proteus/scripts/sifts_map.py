#!/usr/bin/env python3
"""Look up residue mappings from the PDBe SIFTS API.

The helper uses lightweight PDBe JSON endpoints and does not download SIFTS
flatfiles.

Examples:
    python sifts_map.py pdb 1hsg --json
    python sifts_map.py pdb 1hsg --all-isoforms
    python sifts_map.py uniprot P04637 --json
    python sifts_map.py uniprot P04637 --limit 20
"""

import argparse
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PDBE_BASE = "https://www.ebi.ac.uk/pdbe/api"
USER_AGENT = "proteus-sifts-map/1.0"

PDB_ID_RE = re.compile(r"[A-Za-z0-9]{4}")
UNIPROT_ACCESSION_RE = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?"
)


class SiftsLookupError(RuntimeError):
    pass


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _normalize_pdb_id(value: str) -> str:
    pdb_id = value.strip().lower()
    if not PDB_ID_RE.fullmatch(pdb_id):
        raise SiftsLookupError("PDB ID must be exactly four alphanumeric characters")
    return pdb_id


def _normalize_uniprot_accession(value: str) -> str:
    accession = value.strip().upper()
    if not UNIPROT_ACCESSION_RE.fullmatch(accession):
        raise SiftsLookupError("UniProt accession format is not recognized")
    return accession


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _get_json(url: str, timeout: int) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise SiftsLookupError("No SIFTS mapping was found for the requested identifier") from exc
        raise SiftsLookupError(f"PDBe returned HTTP {exc.code}: {exc.reason}") from exc
    except (urllib.error.URLError, TimeoutError) as exc:
        raise SiftsLookupError(f"Failed to query PDBe: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise SiftsLookupError("PDBe returned a response that was not valid JSON") from exc


def _endpoint_for_pdb(pdb_id: str, all_isoforms: bool) -> str:
    endpoint = "all_isoforms" if all_isoforms else "uniprot"
    return f"{PDBE_BASE}/mappings/{endpoint}/{urllib.parse.quote(pdb_id)}"


def _endpoint_for_uniprot(accession: str, all_isoforms: bool) -> str:
    endpoint = "all_isoforms" if all_isoforms else ""
    if endpoint:
        return f"{PDBE_BASE}/mappings/{endpoint}/{urllib.parse.quote(accession)}"
    return f"{PDBE_BASE}/mappings/{urllib.parse.quote(accession)}"


def _residue_number(value):
    if isinstance(value, dict):
        return value.get("residue_number")
    return value


def _author_residue_number(value):
    if isinstance(value, dict):
        return value.get("author_residue_number")
    return None


def _author_insertion_code(value):
    if not isinstance(value, dict):
        return None
    code = value.get("author_insertion_code")
    return code or None


def _segment_record(
    mapping: dict,
    pdb_id: str,
    uniprot_accession: str,
    uniprot_id: str | None = None,
) -> dict:
    start = mapping.get("start")
    end = mapping.get("end")
    pdb_start = mapping.get("pdb_start", _residue_number(start))
    pdb_end = mapping.get("pdb_end", _residue_number(end))

    record = {
        "pdb_id": pdb_id.lower() if pdb_id else None,
        "chain_id": mapping.get("chain_id"),
        "struct_asym_id": mapping.get("struct_asym_id"),
        "entity_id": mapping.get("entity_id"),
        "uniprot_accession": uniprot_accession,
        "uniprot_id": uniprot_id,
        "uniprot_start": mapping.get("unp_start"),
        "uniprot_end": mapping.get("unp_end"),
        "pdb_start": pdb_start,
        "pdb_end": pdb_end,
        "auth_start": _author_residue_number(start),
        "auth_end": _author_residue_number(end),
        "auth_start_ins_code": _author_insertion_code(start),
        "auth_end_ins_code": _author_insertion_code(end),
    }

    for key in ("identity", "coverage", "is_canonical"):
        if key in mapping:
            record[key] = mapping.get(key)
    return record


def normalize_pdb_mappings(raw: dict, pdb_id: str) -> list[dict]:
    pdb_key = pdb_id.lower()
    entry = raw.get(pdb_key) or raw.get(pdb_key.upper()) or {}
    uniprot_block = entry.get("UniProt") or {}
    records = []
    for accession in sorted(uniprot_block):
        accession_data = uniprot_block.get(accession) or {}
        uniprot_id = accession_data.get("identifier") or accession_data.get("name")
        for mapping in accession_data.get("mappings") or []:
            records.append(_segment_record(mapping, pdb_key, accession, uniprot_id))
    return records


def normalize_uniprot_mappings(raw: dict, accession: str) -> list[dict]:
    accession_key = accession.upper()
    entry = raw.get(accession_key) or {}
    pdb_block = entry.get("PDB") or {}
    records = []
    for pdb_id in sorted(pdb_block):
        for mapping in pdb_block.get(pdb_id) or []:
            records.append(_segment_record(mapping, pdb_id, accession_key))
    return records


def _as_int(value) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _clean_ins_code(value) -> str | None:
    if value in {None, "", ".", "?"}:
        return None
    return str(value)


def _normalize_mapping_record(record: dict) -> dict:
    item = dict(record)
    start = item.get("start")
    end = item.get("end")

    if not item.get("chain_id"):
        item["chain_id"] = item.get("chain") or item.get("auth_asym_id") or item.get("structure_chain_id")
    if not item.get("uniprot_accession"):
        item["uniprot_accession"] = (
            item.get("accession")
            or item.get("unp_accession")
            or item.get("uniprot")
        )
    if item.get("uniprot_accession"):
        item["uniprot_accession"] = str(item["uniprot_accession"]).upper()
    if item.get("pdb_id"):
        item["pdb_id"] = str(item["pdb_id"]).lower()
    if "uniprot_start" not in item:
        item["uniprot_start"] = item.get("unp_start")
    if "uniprot_end" not in item:
        item["uniprot_end"] = item.get("unp_end")
    if "pdb_start" not in item:
        item["pdb_start"] = _residue_number(start)
    if "pdb_end" not in item:
        item["pdb_end"] = _residue_number(end)
    if "auth_start" not in item:
        item["auth_start"] = _author_residue_number(start)
    if "auth_end" not in item:
        item["auth_end"] = _author_residue_number(end)
    if "auth_start_ins_code" not in item:
        item["auth_start_ins_code"] = _author_insertion_code(start)
    if "auth_end_ins_code" not in item:
        item["auth_end_ins_code"] = _author_insertion_code(end)

    for key in ("uniprot_start", "uniprot_end", "pdb_start", "pdb_end", "auth_start", "auth_end"):
        parsed = _as_int(item.get(key))
        if parsed is not None:
            item[key] = parsed
    item["auth_start_ins_code"] = _clean_ins_code(item.get("auth_start_ins_code"))
    item["auth_end_ins_code"] = _clean_ins_code(item.get("auth_end_ins_code"))
    return item


def _looks_like_mapping_record(value: dict) -> bool:
    return (
        any(key in value for key in ("uniprot_start", "unp_start"))
        and any(key in value for key in ("chain_id", "chain", "auth_asym_id", "structure_chain_id"))
    )


def normalize_sifts_json_payload(payload) -> list[dict]:
    """Normalize local SIFTS JSON into segment records without network access."""

    if isinstance(payload, dict) and isinstance(payload.get("data"), dict):
        payload = payload["data"]

    if isinstance(payload, dict) and isinstance(payload.get("mappings"), list):
        return [_normalize_mapping_record(record) for record in payload["mappings"] if isinstance(record, dict)]

    if isinstance(payload, list):
        return [_normalize_mapping_record(record) for record in payload if isinstance(record, dict)]

    if isinstance(payload, dict) and _looks_like_mapping_record(payload):
        return [_normalize_mapping_record(payload)]

    records = []
    if isinstance(payload, dict):
        for key, value in payload.items():
            if not isinstance(value, dict):
                continue
            if isinstance(value.get("UniProt"), dict):
                records.extend(normalize_pdb_mappings(payload, key))
            elif isinstance(value.get("PDB"), dict):
                records.extend(normalize_uniprot_mappings(payload, key))
    return [_normalize_mapping_record(record) for record in records]


def load_sifts_json(path: str | Path) -> list[dict]:
    """Load normalized SIFTS records from a local JSON file."""

    try:
        payload = json.loads(Path(path).read_text())
    except OSError as exc:
        reason = exc.strerror or exc.__class__.__name__
        raise SiftsLookupError(f"Failed to read SIFTS JSON file: {reason}") from exc
    except json.JSONDecodeError as exc:
        raise SiftsLookupError("SIFTS JSON file was not valid JSON") from exc

    records = normalize_sifts_json_payload(payload)
    if not records:
        raise SiftsLookupError("No SIFTS mapping records were found in JSON")
    return records


def _accession_matches(record_accession: str | None, accession: str) -> bool:
    if not record_accession:
        return False
    return str(record_accession).upper() == accession.upper()


def _record_chain_values(record: dict) -> set[str]:
    values = {
        record.get("chain_id"),
        record.get("auth_asym_id"),
        record.get("structure_chain_id"),
        record.get("struct_asym_id"),
    }
    return {str(value) for value in values if value not in {None, ""}}


def _record_matches_chain(record: dict, chain_id: str | None) -> bool:
    if chain_id is None:
        return True
    return chain_id in _record_chain_values(record)


def _interpolate_number(start, end, offset: int, segment_length: int) -> int | None:
    start_i = _as_int(start)
    end_i = _as_int(end)
    if start_i is None:
        return None
    if segment_length == 1:
        return start_i
    if end_i is None:
        return None
    delta = end_i - start_i
    if abs(delta) != segment_length - 1:
        return None
    step = 1 if delta >= 0 else -1
    return start_i + (offset * step)


def _ins_code_for_offset(record: dict, offset: int, segment_length: int) -> str | None:
    if offset == 0:
        return _clean_ins_code(record.get("auth_start_ins_code"))
    if offset == segment_length - 1:
        return _clean_ins_code(record.get("auth_end_ins_code"))
    return None


def public_mapping_candidate(candidate: dict) -> dict:
    keys = (
        "pdb_id",
        "chain_id",
        "struct_asym_id",
        "uniprot_accession",
        "uniprot_residue_index",
        "structure_chain_id",
        "structure_residue_id",
        "structure_insertion_code",
        "auth_residue_id",
        "auth_insertion_code",
        "pdb_residue_id",
    )
    return {key: candidate.get(key) for key in keys if candidate.get(key) is not None}


def map_uniprot_residue_candidates(
    records: list[dict],
    accession: str,
    residue_index: int,
    chain_id: str | None = None,
) -> list[dict]:
    """Map a UniProt residue index to structure residue IDs from SIFTS records."""

    normalized_accession = _normalize_uniprot_accession(accession)
    if residue_index < 1:
        raise SiftsLookupError("Residue index must be a positive integer")

    candidates = []
    for raw_record in records:
        record = _normalize_mapping_record(raw_record)
        if not _accession_matches(record.get("uniprot_accession"), normalized_accession):
            continue
        if not _record_matches_chain(record, chain_id):
            continue

        uniprot_start = _as_int(record.get("uniprot_start"))
        uniprot_end = _as_int(record.get("uniprot_end"))
        if uniprot_start is None or uniprot_end is None:
            continue
        lower = min(uniprot_start, uniprot_end)
        upper = max(uniprot_start, uniprot_end)
        if not (lower <= residue_index <= upper):
            continue

        offset = abs(residue_index - uniprot_start)
        segment_length = abs(uniprot_end - uniprot_start) + 1
        auth_number = _interpolate_number(record.get("auth_start"), record.get("auth_end"), offset, segment_length)
        pdb_number = _interpolate_number(record.get("pdb_start"), record.get("pdb_end"), offset, segment_length)
        structure_number = auth_number if auth_number is not None else pdb_number
        if structure_number is None:
            continue

        structure_chain_id = (
            record.get("chain_id")
            or record.get("auth_asym_id")
            or record.get("structure_chain_id")
            or record.get("struct_asym_id")
        )
        if not structure_chain_id:
            continue
        auth_ins_code = _ins_code_for_offset(record, offset, segment_length) if auth_number is not None else None
        candidate = {
            "pdb_id": record.get("pdb_id"),
            "chain_id": record.get("chain_id"),
            "struct_asym_id": record.get("struct_asym_id"),
            "uniprot_accession": normalized_accession,
            "uniprot_residue_index": residue_index,
            "structure_chain_id": str(structure_chain_id),
            "structure_residue_id": str(structure_number),
            "structure_insertion_code": auth_ins_code,
            "auth_residue_id": str(auth_number) if auth_number is not None else None,
            "auth_insertion_code": auth_ins_code,
            "pdb_residue_id": str(pdb_number) if pdb_number is not None else None,
            "mapping": record,
        }
        candidates.append(candidate)

    return candidates


def _limit_records(records: list[dict], limit: int | None) -> tuple[list[dict], bool]:
    if limit is None or len(records) <= limit:
        return records, False
    return records[:limit], True


def lookup_pdb(pdb_id: str, all_isoforms: bool = False, limit: int | None = None, timeout: int = 30) -> dict:
    normalized = _normalize_pdb_id(pdb_id)
    url = _endpoint_for_pdb(normalized, all_isoforms)
    raw = _get_json(url, timeout)
    records = normalize_pdb_mappings(raw, normalized)
    if not records:
        raise SiftsLookupError("No SIFTS mappings were found for the requested PDB ID")

    returned, truncated = _limit_records(records, limit)
    return _ok_payload({
        "mode": "pdb",
        "query": normalized,
        "source": "PDBe SIFTS API",
        "all_isoforms": all_isoforms,
        "mapping_count": len(records),
        "returned_count": len(returned),
        "truncated": truncated,
        "mappings": returned,
    })


def lookup_uniprot(accession: str, all_isoforms: bool = False, limit: int | None = None, timeout: int = 30) -> dict:
    normalized = _normalize_uniprot_accession(accession)
    url = _endpoint_for_uniprot(normalized, all_isoforms)
    raw = _get_json(url, timeout)
    records = normalize_uniprot_mappings(raw, normalized)
    if not records:
        raise SiftsLookupError("No SIFTS mappings were found for the requested UniProt accession")

    returned, truncated = _limit_records(records, limit)
    return _ok_payload({
        "mode": "uniprot",
        "query": normalized,
        "source": "PDBe SIFTS API",
        "all_isoforms": all_isoforms,
        "mapping_count": len(records),
        "returned_count": len(returned),
        "truncated": truncated,
        "mappings": returned,
    })


def _format_text(output: dict) -> str:
    mappings = output["data"]["mappings"]
    query = output["data"]["query"]
    mode = output["data"]["mode"]
    header = f"{len(mappings)} returned SIFTS mapping(s) for {mode} {query}"
    if output["data"]["truncated"]:
        header += f" ({output['data']['mapping_count']} total)"
    lines = [header + ":"]
    for record in mappings:
        unp_range = f"{record.get('uniprot_start')}-{record.get('uniprot_end')}"
        pdb_range = f"{record.get('pdb_start')}-{record.get('pdb_end')}"
        identity = record.get("identity")
        identity_text = f"; identity={identity}" if identity is not None else ""
        lines.append(
            "  "
            f"{record.get('pdb_id')} chain {record.get('chain_id')} "
            f"<-> {record.get('uniprot_accession')} {unp_range} "
            f"(PDB {pdb_range}{identity_text})"
        )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Look up residue-level PDB and UniProt mappings from PDBe SIFTS.",
        epilog=(
            "Examples:\n"
            "  %(prog)s pdb 1hsg --json\n"
            "  %(prog)s pdb 1hsg --all-isoforms\n"
            "  %(prog)s uniprot P04637 --json\n"
            "  %(prog)s uniprot P04637 --limit 20"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    common.add_argument("--all-isoforms", action="store_true", help="Include isoform mappings when PDBe provides them")
    common.add_argument("--limit", type=_positive_int, help="Limit returned mapping records")
    common.add_argument("--timeout", type=_positive_int, default=30, help="Request timeout in seconds (default: 30)")

    subparsers = parser.add_subparsers(dest="mode", required=True)
    pdb_parser = subparsers.add_parser("pdb", parents=[common], help="Map a PDB ID to UniProt segments")
    pdb_parser.add_argument("identifier", help="Four-character PDB ID, e.g. 1hsg")

    uniprot_parser = subparsers.add_parser(
        "uniprot",
        parents=[common],
        help="Map a UniProt accession to PDB segments",
    )
    uniprot_parser.add_argument("identifier", help="UniProt accession, e.g. P04637")
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.mode == "pdb":
            output = lookup_pdb(args.identifier, args.all_isoforms, args.limit, args.timeout)
        else:
            output = lookup_uniprot(args.identifier, args.all_isoforms, args.limit, args.timeout)
    except SiftsLookupError as exc:
        payload = _error_payload(str(exc))
        if args.json:
            print(json.dumps(payload, indent=2))
        else:
            print(f"ERROR: {payload['error']}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print(_format_text(output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
