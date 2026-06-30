#!/usr/bin/env python3
"""Triage protein substitutions against a local PDB/mmCIF structure.

This helper is intentionally local and conservative. It parses simple protein
substitutions through variant_map.py, finds matching structure residues by
direct structure numbering, summarizes nearby ligands and protein contacts, and
optionally annotates residues from a local AlphaFold PAE JSON file.

Usage:
    python3 mutation_triage.py R175H --structure model.pdb --json
    python3 mutation_triage.py "P04637 R175H" G245S --structure model.cif --chain A
    python3 mutation_triage.py R175H --uniprot P04637 --structure AF-P04637-F1.pdb --pae pae.json
"""

import argparse
import json
import math
import shlex
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import structure_info
import sifts_map
import variant_map


WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}
ALTLOCS_TO_KEEP = {"", "A", "1"}
NUMBERING_WARNING = (
    "Residue lookup uses direct structure residue numbering only; no SIFTS/PDBe "
    "UniProt-to-PDB residue mapping was applied."
)
PAE_WARNING = (
    "PAE rows/columns are indexed directly by residue_index; no structure-to-sequence "
    "mapping was applied."
)
ROOT = Path(__file__).resolve().parents[1]


class MutationTriageError(ValueError):
    """Raised when local mutation triage cannot be completed."""


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _coord_payload(x: float, y: float, z: float) -> dict:
    return {"x": round(x, 3), "y": round(y, 3), "z": round(z, 3)}


def _display_path(path: str | Path | None) -> str | None:
    """Return a report-safe path label without exposing absolute home paths."""

    if path is None:
        return None
    text = str(path)
    if "://" in text:
        return text
    candidate = Path(text).expanduser()
    if not candidate.is_absolute():
        return text
    resolved = candidate.resolve()
    for root, prefix in ((ROOT, "."), (Path.cwd(), ".")):
        try:
            relative = resolved.relative_to(root.resolve())
        except ValueError:
            continue
        return f"{prefix}/{relative.as_posix()}"
    return f"{resolved.name} (absolute path omitted)"


def _sanitize_report_paths(value):
    if isinstance(value, dict):
        sanitized = {}
        for key, item in value.items():
            if key in {"file", "path", "structure", "pae_json"} and isinstance(item, str):
                sanitized[key] = _display_path(item)
            else:
                sanitized[key] = _sanitize_report_paths(item)
        return sanitized
    if isinstance(value, list):
        return [_sanitize_report_paths(item) for item in value]
    return value


def _clean_cif_missing(value: str | None) -> str | None:
    if value is None or value in {".", "?"}:
        return None
    return value


def _residue_label(chain: str, residue_id: str, insertion_code: str | None) -> str:
    suffix = insertion_code or ""
    return f"{chain}:{residue_id}{suffix}"


def _sort_residue_id(value: str) -> tuple[int, int | str]:
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)


def _residue_sort_key(item: dict) -> tuple:
    return (
        item["chain"],
        _sort_residue_id(item["residue_id"]),
        item.get("insertion_code") or "",
        item["resname"],
    )


def _seq_matches(seq_id: str, residue_index: int) -> bool:
    if seq_id == str(residue_index):
        return True
    try:
        return int(seq_id) == residue_index
    except ValueError:
        return False


def _residue_id_matches(seq_id: str, target_id: str) -> bool:
    if seq_id == target_id:
        return True
    try:
        return int(seq_id) == int(target_id)
    except ValueError:
        return False


def _residue_matches_sifts_target(residue: dict, target: dict) -> bool:
    if residue["chain"] != target["structure_chain_id"]:
        return False
    target_ins = target.get("structure_insertion_code")
    if target_ins is not None and residue.get("insertion_code") != target_ins:
        return False
    return _residue_id_matches(residue["residue_id"], target["structure_residue_id"])


def _seq_int(seq_id: str) -> int | None:
    try:
        return int(seq_id)
    except ValueError:
        return None


def _distance(a: dict, b: dict) -> float:
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _min_atom_distance(left_atoms: list[dict], right_atoms: list[dict]) -> dict | None:
    best = None
    for left in left_atoms:
        for right in right_atoms:
            distance = _distance(left, right)
            if best is None or distance < best["distance"]:
                best = {
                    "distance": distance,
                    "target_atom": left["atom"],
                    "neighbor_atom": right["atom"],
                }
    if best is None:
        return None
    best["min_distance"] = round(best.pop("distance"), 2)
    return best


def _atom_payload(record: str, atom: str, resname: str, chain: str, residue_id: str,
                  insertion_code: str | None, x: float, y: float, z: float,
                  element: str = "") -> dict:
    return {
        "record": record,
        "atom": atom or "?",
        "resname": resname.upper() or "?",
        "chain": chain or "?",
        "residue_id": residue_id or "?",
        "insertion_code": insertion_code,
        "residue_label": _residue_label(chain or "?", residue_id or "?", insertion_code),
        "x": x,
        "y": y,
        "z": z,
        "element": element,
    }


def _parse_pdb_atoms(path: str) -> list[dict]:
    atoms = []
    with open(path) as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            altloc = line[16].strip()
            if altloc not in ALTLOCS_TO_KEEP:
                continue
            try:
                atoms.append(_atom_payload(
                    record=record,
                    atom=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21].strip() or "?",
                    residue_id=line[22:26].strip() or "?",
                    insertion_code=line[26].strip() or None,
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    element=line[76:78].strip() if len(line) >= 78 else "",
                ))
            except ValueError:
                continue
    return atoms


def _field_index(field_index: dict[str, int], *names: str) -> int | None:
    for name in names:
        if name in field_index:
            return field_index[name]
    return None


def _row_value(row: list[str], index: int | None, default: str | None = None) -> str | None:
    if index is None or index >= len(row):
        return default
    return row[index]


def _parse_mmcif_atoms(path: str) -> list[dict]:
    lines = Path(path).read_text(errors="replace").splitlines()
    atoms = []
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue

        index += 1
        headers = []
        while index < len(lines) and lines[index].strip().startswith("_"):
            headers.append(lines[index].strip())
            index += 1

        if not headers or not all(header.startswith("_atom_site.") for header in headers):
            while (
                index < len(lines)
                and lines[index].strip()
                and not lines[index].strip().startswith(("loop_", "_", "#"))
            ):
                index += 1
            continue

        fields = {header: pos for pos, header in enumerate(headers)}
        group_i = _field_index(fields, "_atom_site.group_PDB")
        atom_i = _field_index(fields, "_atom_site.label_atom_id", "_atom_site.auth_atom_id")
        comp_i = _field_index(fields, "_atom_site.label_comp_id", "_atom_site.auth_comp_id")
        chain_i = _field_index(fields, "_atom_site.auth_asym_id", "_atom_site.label_asym_id")
        seq_i = _field_index(fields, "_atom_site.auth_seq_id", "_atom_site.label_seq_id")
        ins_i = _field_index(fields, "_atom_site.pdbx_PDB_ins_code")
        alt_i = _field_index(fields, "_atom_site.label_alt_id")
        x_i = _field_index(fields, "_atom_site.Cartn_x")
        y_i = _field_index(fields, "_atom_site.Cartn_y")
        z_i = _field_index(fields, "_atom_site.Cartn_z")
        element_i = _field_index(fields, "_atom_site.type_symbol")

        while index < len(lines):
            stripped = lines[index].strip()
            if not stripped or stripped.startswith("#"):
                index += 1
                break
            if stripped == "loop_" or stripped.startswith("_"):
                break
            try:
                row = shlex.split(stripped, posix=True)
            except ValueError:
                index += 1
                continue
            index += 1
            if len(row) < len(headers):
                continue

            record = (_row_value(row, group_i, "ATOM") or "ATOM").upper()
            if record not in {"ATOM", "HETATM"}:
                continue
            altloc = _clean_cif_missing(_row_value(row, alt_i)) or ""
            if altloc not in ALTLOCS_TO_KEEP:
                continue
            x = _row_value(row, x_i)
            y = _row_value(row, y_i)
            z = _row_value(row, z_i)
            if x is None or y is None or z is None:
                continue
            try:
                atoms.append(_atom_payload(
                    record=record,
                    atom=_row_value(row, atom_i, "?") or "?",
                    resname=_row_value(row, comp_i, "?") or "?",
                    chain=_clean_cif_missing(_row_value(row, chain_i)) or "?",
                    residue_id=_clean_cif_missing(_row_value(row, seq_i)) or "?",
                    insertion_code=_clean_cif_missing(_row_value(row, ins_i)),
                    x=float(x),
                    y=float(y),
                    z=float(z),
                    element=_row_value(row, element_i, "") or "",
                ))
            except ValueError:
                continue
    return atoms


def _group_residues(atoms: list[dict], record: str) -> list[dict]:
    grouped: dict[tuple[str, str, str | None, str], dict] = {}
    for atom in atoms:
        if atom["record"] != record:
            continue
        key = (atom["chain"], atom["residue_id"], atom["insertion_code"], atom["resname"])
        entry = grouped.setdefault(key, {
            "chain": atom["chain"],
            "residue_id": atom["residue_id"],
            "insertion_code": atom["insertion_code"],
            "residue_label": atom["residue_label"],
            "resname": atom["resname"],
            "residue_one_letter": variant_map.AA3_TO_1.get(atom["resname"]),
            "atom_count": 0,
            "ca_coordinate": None,
            "atoms": [],
        })
        entry["atom_count"] += 1
        entry["atoms"].append(atom)
        if atom["atom"] == "CA" and entry["ca_coordinate"] is None:
            entry["ca_coordinate"] = _coord_payload(atom["x"], atom["y"], atom["z"])
    return sorted(grouped.values(), key=_residue_sort_key)


def _group_ligands(atoms: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str, str | None, str], dict] = {}
    for atom in atoms:
        if atom["record"] != "HETATM" or atom["resname"] in WATER_NAMES:
            continue
        key = (atom["chain"], atom["residue_id"], atom["insertion_code"], atom["resname"])
        entry = grouped.setdefault(key, {
            "ligand": atom["resname"],
            "chain": atom["chain"],
            "residue_id": atom["residue_id"],
            "insertion_code": atom["insertion_code"],
            "residue_label": atom["residue_label"],
            "atom_count": 0,
            "atoms": [],
        })
        entry["atom_count"] += 1
        entry["atoms"].append(atom)
    return sorted(grouped.values(), key=lambda item: (
        item["ligand"],
        item["chain"],
        _sort_residue_id(item["residue_id"]),
        item.get("insertion_code") or "",
    ))


def _public_residue(entry: dict) -> dict:
    return {
        "chain": entry["chain"],
        "residue_id": entry["residue_id"],
        "insertion_code": entry["insertion_code"],
        "residue_label": entry["residue_label"],
        "resname": entry["resname"],
        "residue_one_letter": entry["residue_one_letter"],
        "atom_count": entry["atom_count"],
        "ca_coordinate": entry["ca_coordinate"],
    }


def _public_ligand(entry: dict) -> dict:
    return {
        "ligand": entry["ligand"],
        "chain": entry["chain"],
        "residue_id": entry["residue_id"],
        "insertion_code": entry["insertion_code"],
        "residue_label": entry["residue_label"],
        "atom_count": entry["atom_count"],
    }


def _load_structure(path: str) -> dict:
    structure_path = Path(path)
    if not structure_path.exists():
        raise MutationTriageError(f"Structure file not found: {_display_path(path)}")
    suffix = structure_path.suffix.lower()
    fmt = "mmcif" if suffix in {".cif", ".mmcif"} else "pdb"
    atoms = _parse_mmcif_atoms(str(structure_path)) if fmt == "mmcif" else _parse_pdb_atoms(str(structure_path))
    protein_residues = _group_residues(atoms, "ATOM")
    ligands = _group_ligands(atoms)
    inspection = structure_info.inspect_structure(str(structure_path))
    return {
        "file": _display_path(structure_path.resolve()),
        "format": fmt,
        "inspection": _sanitize_report_paths(inspection),
        "atoms": atoms,
        "protein_residues": protein_residues,
        "ligands": ligands,
    }


def _matrix_from_payload(payload) -> list[list[float]]:
    if isinstance(payload, dict):
        for key in ("predicted_aligned_error", "pae", "pae_matrix"):
            if key in payload:
                return payload[key]
    if isinstance(payload, list):
        if payload and isinstance(payload[0], dict):
            for key in ("predicted_aligned_error", "pae", "pae_matrix"):
                if key in payload[0]:
                    return payload[0][key]
        if payload and isinstance(payload[0], list):
            return payload
    raise MutationTriageError("Could not find a PAE matrix in the JSON payload.")


def _validate_matrix(matrix) -> list[list[float]]:
    if not isinstance(matrix, list) or not matrix:
        raise MutationTriageError("PAE matrix is empty.")
    width = len(matrix[0])
    if width == 0:
        raise MutationTriageError("PAE matrix rows are empty.")
    normalized = []
    for row in matrix:
        if not isinstance(row, list) or len(row) != width:
            raise MutationTriageError("PAE matrix must be rectangular.")
        normalized.append([float(value) for value in row])
    return normalized


def _load_pae(path: str | None) -> dict | None:
    if not path:
        return None
    pae_path = Path(path)
    if not pae_path.exists():
        raise MutationTriageError(f"PAE JSON file not found: {_display_path(path)}")
    payload = json.loads(pae_path.read_text())
    matrix = _validate_matrix(_matrix_from_payload(payload))
    flat = [value for row in matrix for value in row]
    return {
        "file": _display_path(pae_path.resolve()),
        "matrix": matrix,
        "size": {"rows": len(matrix), "columns": len(matrix[0])},
        "pae": {
            "min": round(min(flat), 2),
            "max": round(max(flat), 2),
            "mean": round(sum(flat) / len(flat), 2),
        },
    }


def _pae_pair(matrix: list[list[float]], left_index: int, right_residue_id: str) -> dict | None:
    right_index = _seq_int(right_residue_id)
    if right_index is None:
        return None
    rows = len(matrix)
    columns = len(matrix[0])
    if not (1 <= left_index <= rows and 1 <= right_index <= columns):
        return None
    forward = matrix[left_index - 1][right_index - 1]
    reverse = matrix[right_index - 1][left_index - 1] if right_index <= rows and left_index <= columns else forward
    return {
        "residue_index": right_index,
        "pae_from_variant": round(forward, 2),
        "pae_to_variant": round(reverse, 2),
        "mean_pair_pae": round((forward + reverse) / 2, 2),
    }


def _pae_residue_summary(pae: dict | None, residue_index: int) -> dict | None:
    if pae is None:
        return None
    matrix = pae["matrix"]
    if not (1 <= residue_index <= len(matrix)):
        return {
            "available": False,
            "reason": f"Residue index {residue_index} is outside PAE matrix rows ({len(matrix)}).",
        }
    row = matrix[residue_index - 1]
    summary = {
        "available": True,
        "residue_index": residue_index,
        "matrix_index": residue_index - 1,
        "row_mean_pae": round(sum(row) / len(row), 2),
        "row_min_pae": round(min(row), 2),
        "row_max_pae": round(max(row), 2),
        "note": PAE_WARNING,
    }
    if residue_index <= len(row):
        summary["self_pae"] = round(row[residue_index - 1], 2)
    return summary


def _contact_payload(target: dict, neighbor: dict, pae: dict | None, residue_index: int) -> dict | None:
    contact = _min_atom_distance(target["atoms"], neighbor["atoms"])
    if contact is None:
        return None
    payload = _public_residue(neighbor)
    payload.update(contact)
    if pae is not None:
        pair = _pae_pair(pae["matrix"], residue_index, neighbor["residue_id"])
        if pair is not None:
            payload["pae_pair"] = pair
    return payload


def _chain_summary(contacts: list[dict]) -> list[dict]:
    by_chain: dict[str, dict] = {}
    for contact in contacts:
        item = by_chain.setdefault(contact["chain"], {
            "chain": contact["chain"],
            "contact_count": 0,
            "closest_distance": contact["min_distance"],
        })
        item["contact_count"] += 1
        item["closest_distance"] = min(item["closest_distance"], contact["min_distance"])
    return sorted(by_chain.values(), key=lambda item: (item["closest_distance"], item["chain"]))


def _triage_candidate(target: dict, variant: dict, residues: list[dict], ligands: list[dict],
                      cutoff: float, max_hits: int, pae: dict | None) -> dict:
    residue_index = variant["residue_index"]
    result = _public_residue(target)
    result["reference_matches_variant"] = (
        result["residue_one_letter"] == variant["from"]["one_letter"]
        if result["residue_one_letter"] else None
    )
    result["flags"] = []

    ligand_distances = []
    for ligand in ligands:
        contact = _min_atom_distance(target["atoms"], ligand["atoms"])
        if contact is None:
            continue
        payload = _public_ligand(ligand)
        payload.update(contact)
        ligand_distances.append(payload)
    ligand_distances.sort(key=lambda item: item["min_distance"])
    nearby_ligands = [item for item in ligand_distances if item["min_distance"] <= cutoff]
    if nearby_ligands:
        result["flags"].append("near_ligand")
    result["ligand_contact_count"] = len(nearby_ligands)
    result["nearby_ligands"] = nearby_ligands[:max_hits]
    result["nearest_ligands"] = ligand_distances[:max_hits]

    same_chain_contacts = []
    other_chain_contacts = []
    target_key = (target["chain"], target["residue_id"], target["insertion_code"], target["resname"])
    for residue in residues:
        residue_key = (residue["chain"], residue["residue_id"], residue["insertion_code"], residue["resname"])
        if residue_key == target_key:
            continue
        contact = _contact_payload(target, residue, pae, residue_index)
        if contact is None or contact["min_distance"] > cutoff:
            continue
        if residue["chain"] == target["chain"]:
            same_chain_contacts.append(contact)
        else:
            other_chain_contacts.append(contact)

    same_chain_contacts.sort(key=lambda item: item["min_distance"])
    other_chain_contacts.sort(key=lambda item: item["min_distance"])
    if same_chain_contacts:
        result["flags"].append("same_chain_close_contacts")
    if other_chain_contacts:
        result["flags"].append("near_other_chain")
    result["same_chain_contact_count"] = len(same_chain_contacts)
    result["same_chain_contacts"] = same_chain_contacts[:max_hits]
    result["other_chain_contact_count"] = len(other_chain_contacts)
    result["other_chain_contacts"] = other_chain_contacts[:max_hits]
    result["other_chain_summary"] = _chain_summary(other_chain_contacts)
    result["pae"] = _pae_residue_summary(pae, residue_index)
    return result


def _variant_warnings(parsed: dict, candidates: list[dict], chain: str | None) -> list[str]:
    warnings = []
    if not candidates:
        warnings.append("No structure residue matched this residue_index and chain filter.")
    if len(candidates) > 1 and chain is None:
        warnings.append("Multiple chains/residues matched this residue_index; pass --chain to disambiguate.")
    for candidate in candidates:
        if candidate.get("ca_coordinate") is None:
            warnings.append(f"No CA atom was found for {candidate['residue_label']} {candidate['resname']}.")
        if candidate.get("reference_matches_variant") is False:
            warnings.append(
                f"Structure residue {candidate['residue_label']} {candidate['resname']} does not match "
                f"variant reference {parsed['from']['three_letter']}."
            )
    return warnings


def triage_variants(structure: str, variants: list[str], *, uniprot_id: str | None = None,
                    chain: str | None = None, cutoff: float = 5.0,
                    pae_json: str | None = None, max_hits: int = 8,
                    sifts_json: str | None = None) -> dict:
    """Return a conservative local structural triage report for substitutions."""

    if cutoff <= 0:
        raise MutationTriageError("--cutoff must be positive.")
    if max_hits <= 0:
        raise MutationTriageError("--max-hits must be positive.")
    if not variants:
        raise MutationTriageError("At least one variant is required.")

    parsed_structure = _load_structure(structure)
    pae = _load_pae(pae_json)
    sifts_records = sifts_map.load_sifts_json(sifts_json) if sifts_json else None
    sifts_source = _display_path(sifts_json) if sifts_json else None
    warnings = [] if sifts_records is not None else [NUMBERING_WARNING]
    if pae is not None:
        warnings.append(PAE_WARNING)

    results = []
    any_sifts_applied = False
    any_sifts_fallback = False
    for variant_text in variants:
        parsed = variant_map.parse_variant(variant_text, uniprot_id=uniprot_id)
        variant_pre_warnings = []
        sifts_candidates = []
        variant_sifts_mapping = {
            "applied": False,
            "status": "not_supplied",
        }
        if sifts_records is not None:
            variant_sifts_mapping["source"] = sifts_source
            if not parsed.get("uniprot_id"):
                any_sifts_fallback = True
                warning = (
                    "SIFTS JSON was supplied, but no UniProt accession was available; "
                    "direct structure residue numbering was used."
                )
                variant_pre_warnings.append(warning)
                variant_sifts_mapping.update({"status": "missing_uniprot", "warning": warning})
            else:
                sifts_candidates = sifts_map.map_uniprot_residue_candidates(
                    sifts_records,
                    parsed["uniprot_id"],
                    parsed["residue_index"],
                    chain_id=chain,
                )
                if sifts_candidates:
                    any_sifts_applied = True
                    variant_sifts_mapping = {
                        "applied": True,
                        "status": "mapped",
                        "source": sifts_source,
                        "candidate_count": len(sifts_candidates),
                        "candidates": [
                            sifts_map.public_mapping_candidate(candidate)
                            for candidate in sifts_candidates
                        ],
                    }
                else:
                    any_sifts_fallback = True
                    warning = (
                        "SIFTS JSON was supplied, but no mapping matched this UniProt "
                        "residue and chain filter; direct structure residue numbering was used."
                    )
                    variant_pre_warnings.append(warning)
                    variant_sifts_mapping.update({"status": "no_match", "warning": warning})

        candidate_pairs = []
        if sifts_candidates:
            seen = set()
            for target in sifts_candidates:
                for residue in parsed_structure["protein_residues"]:
                    if not _residue_matches_sifts_target(residue, target):
                        continue
                    key = (residue["chain"], residue["residue_id"], residue.get("insertion_code"), residue["resname"])
                    if key in seen:
                        continue
                    seen.add(key)
                    candidate_pairs.append((residue, target))
        else:
            candidate_pairs = [
                (residue, None)
                for residue in parsed_structure["protein_residues"]
                if _seq_matches(residue["residue_id"], parsed["residue_index"])
                and (chain is None or residue["chain"] == chain)
            ]

        candidate_results = []
        for candidate, sifts_candidate in candidate_pairs:
            candidate_result = _triage_candidate(
                candidate,
                parsed,
                parsed_structure["protein_residues"],
                parsed_structure["ligands"],
                cutoff,
                max_hits,
                pae,
            )
            if sifts_candidate is not None:
                candidate_result["sifts_mapping"] = sifts_map.public_mapping_candidate(sifts_candidate)
            candidate_results.append(candidate_result)

        variant_warnings = variant_pre_warnings + _variant_warnings(parsed, candidate_results, chain)
        flags = sorted({flag for candidate in candidate_results for flag in candidate["flags"]})
        results.append({
            "input": variant_text,
            "variant": parsed,
            "residue_index": parsed["residue_index"],
            "sifts_mapping": variant_sifts_mapping,
            "candidate_count": len(candidate_results),
            "candidates": candidate_results,
            "flags": flags,
            "warnings": variant_warnings,
        })

    if sifts_records is not None and any_sifts_fallback:
        warnings.append(
            "One or more variants used direct structure residue numbering because no SIFTS mapping was applied."
        )

    inspection = parsed_structure["inspection"]
    structure_summary = {
        "source": "local",
        "file": parsed_structure["file"],
        "format": parsed_structure["format"],
        "chains": inspection.get("chains", []),
        "atom_records": inspection.get("atom_records"),
        "hetatm_records": inspection.get("hetatm_records"),
        "protein_residue_count": len(parsed_structure["protein_residues"]),
        "ligand_group_count": len(parsed_structure["ligands"]),
        "ligand_atom_count": sum(ligand["atom_count"] for ligand in parsed_structure["ligands"]),
        "likely_alphafold": inspection.get("likely_alphafold", False),
        "inspection": inspection,
    }
    pae_summary = None
    if pae is not None:
        pae_summary = {
            "file": pae["file"],
            "size": pae["size"],
            "pae": pae["pae"],
            "note": PAE_WARNING,
        }
    numbering_warning = None if any_sifts_applied and not any_sifts_fallback else NUMBERING_WARNING
    sifts_summary = {"applied": any_sifts_applied}
    if sifts_source is not None:
        sifts_summary["source"] = sifts_source
        sifts_summary["status"] = (
            "mapped" if any_sifts_applied and not any_sifts_fallback
            else "partial" if any_sifts_applied
            else "no_match"
        )
    data = {
        "input": {
            "variants": variants,
            "structure": _display_path(structure),
            "chain_filter": chain,
            "cutoff_angstrom": cutoff,
            "max_hits": max_hits,
            "pae_json": _display_path(pae_json),
            "sifts_json": sifts_source,
        },
        "structure": structure_summary,
        "numbering": {
            "method": (
                "sifts_uniprot_to_structure_mapping"
                if any_sifts_applied else "direct_structure_residue_numbering"
            ),
            "sifts_mapping": sifts_summary,
            "warning": numbering_warning,
        },
        "pae": pae_summary,
        "variants": results,
        "warnings": warnings,
    }
    return _ok_payload(data)


def _yes_no_unknown(value) -> str:
    if value is True:
        return "yes"
    if value is False:
        return "no"
    return "unknown"


def _format_contact_list(items: list[dict], empty_text: str) -> str:
    if not items:
        return empty_text
    chunks = []
    for item in items:
        label = item.get("residue_label") or item.get("ligand") or "?"
        name = item.get("resname") or item.get("ligand") or "?"
        chunks.append(f"{name} {label} at {item['min_distance']} A")
    return "; ".join(chunks)


def format_markdown(output: dict) -> str:
    """Render the triage output as readable Markdown."""

    data = output["data"]
    numbering = data["numbering"]
    sifts_mapping = numbering.get("sifts_mapping") or {}
    if sifts_mapping.get("applied"):
        numbering_text = "SIFTS UniProt-to-structure mapping applied"
        if sifts_mapping.get("source"):
            numbering_text += f" from {sifts_mapping['source']}"
    else:
        numbering_text = "direct structure residue numbering; SIFTS mapping was not applied"
    lines = [
        "# Mutation Triage",
        "",
        f"- Structure: {data['structure']['file']}",
        f"- Format: {data['structure']['format']}",
        f"- Cutoff: {data['input']['cutoff_angstrom']} A",
        f"- Numbering: {numbering_text}.",
    ]
    if data.get("pae"):
        pae = data["pae"]
        lines.append(f"- PAE: {pae['file']} ({pae['size']['rows']} x {pae['size']['columns']})")
    if data["warnings"]:
        lines.extend(["", "## Warnings"])
        for warning in data["warnings"]:
            lines.append(f"- {warning}")

    for item in data["variants"]:
        variant = item["variant"]
        lines.extend([
            "",
            f"## {variant['short']}",
            "",
            f"- Input: {item['input']}",
            f"- HGVS protein: {variant['hgvs_protein']}",
            f"- Residue index: {item['residue_index']}",
            f"- Candidate residues: {item['candidate_count']}",
        ])
        if item["flags"]:
            lines.append(f"- Flags: {', '.join(item['flags'])}")
        for warning in item["warnings"]:
            lines.append(f"- Warning: {warning}")
        if not item["candidates"]:
            continue

        for candidate in item["candidates"]:
            lines.extend([
                "",
                f"### {candidate['residue_label']} {candidate['resname']}",
                "",
                f"- Reference matches variant: {_yes_no_unknown(candidate['reference_matches_variant'])}",
                f"- CA coordinate: {candidate['ca_coordinate'] or 'not found'}",
                f"- Ligands within cutoff: {_format_contact_list(candidate['nearby_ligands'], 'none')}",
                f"- Other-chain contacts within cutoff: {_format_contact_list(candidate['other_chain_contacts'], 'none')}",
                f"- Same-chain close contacts: {_format_contact_list(candidate['same_chain_contacts'], 'none')}",
            ])
            if candidate.get("other_chain_summary"):
                summary = [
                    f"{entry['chain']} ({entry['contact_count']} contacts, closest {entry['closest_distance']} A)"
                    for entry in candidate["other_chain_summary"]
                ]
                lines.append(f"- Other-chain summary: {'; '.join(summary)}")
            if candidate.get("pae"):
                pae = candidate["pae"]
                if pae.get("available"):
                    lines.append(
                        f"- PAE row mean: {pae['row_mean_pae']} A "
                        f"(min {pae['row_min_pae']}, max {pae['row_max_pae']})"
                    )
                else:
                    lines.append(f"- PAE: {pae['reason']}")
    return "\n".join(lines) + "\n"


def _coalesce_variant_args(values: list[str], uniprot_id: str | None) -> list[str]:
    variants = []
    index = 0
    while index < len(values):
        if index + 1 < len(values):
            combined = f"{values[index]} {values[index + 1]}"
            try:
                variant_map.parse_variant(combined, uniprot_id=uniprot_id)
            except variant_map.VariantMapError:
                pass
            else:
                variants.append(combined)
                index += 2
                continue
        variants.append(values[index])
        index += 1
    return variants


def main():
    parser = argparse.ArgumentParser(
        description="Local stdlib-only structural triage for protein substitutions.",
        epilog=(
            "Examples:\n"
            "  %(prog)s R175H --structure model.pdb --json\n"
            "  %(prog)s 'P04637 R175H' G245S --structure model.cif --chain A\n"
            "  %(prog)s R175H --uniprot P04637 --structure AF-P04637-F1.pdb --pae pae.json\n"
            "  %(prog)s G175A --uniprot P04637 --structure model.pdb --sifts-json sifts.json --json\n"
            "\nNotes:\n"
            "  This tool never downloads structures or mappings. Without --sifts-json,\n"
            "  residues are matched by direct structure numbering; local PDBe SIFTS JSON\n"
            "  can map UniProt residue indexes to structure auth residue IDs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("variants", nargs="+", help="One or more substitutions, e.g. R175H or 'P04637 R175H'")
    parser.add_argument("--structure", required=True, help="Local .pdb, .cif, or .mmcif structure file")
    parser.add_argument("--uniprot", help="UniProt accession to apply to variants without a prefix")
    parser.add_argument("--chain", help="Optional structure chain ID filter")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Contact/proximity cutoff in Angstroms (default: 5)")
    parser.add_argument("--pae", dest="pae_json", help="Optional local AlphaFold PAE JSON file")
    parser.add_argument("--sifts-json", help="Local PDBe SIFTS-style JSON mapping file")
    parser.add_argument("--max-hits", type=int, default=8, help="Maximum contacts/ligands to show per section (default: 8)")
    parser.add_argument("--json", action="store_true", dest="as_json", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        variant_inputs = _coalesce_variant_args(args.variants, args.uniprot)
        output = triage_variants(
            args.structure,
            variant_inputs,
            uniprot_id=args.uniprot,
            chain=args.chain,
            cutoff=args.cutoff,
            pae_json=args.pae_json,
            sifts_json=args.sifts_json,
            max_hits=args.max_hits,
        )
    except (
        MutationTriageError,
        variant_map.VariantMapError,
        OSError,
        json.JSONDecodeError,
        sifts_map.SiftsLookupError,
        ValueError,
    ) as exc:
        if args.as_json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.as_json:
        print(json.dumps(output, indent=2))
    else:
        print(format_markdown(output), end="")


if __name__ == "__main__":
    main()
