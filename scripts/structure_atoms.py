#!/usr/bin/env python3
"""Small PDB/mmCIF atom-site parser shared by local Proteus helpers."""

from __future__ import annotations

import shlex
from pathlib import Path


ALTLOCS_TO_KEEP = {"", ".", "?", "A", "1"}
TWO_LETTER_ELEMENTS = {
    "AC", "AG", "AL", "AM", "AR", "AS", "AT", "AU", "BA", "BE", "BI", "BK",
    "BR", "CA", "CD", "CE", "CF", "CL", "CM", "CO", "CR", "CS", "CU", "DY",
    "ER", "ES", "EU", "FE", "FM", "FR", "GA", "GD", "GE", "HE", "HF", "HG",
    "HO", "IN", "IR", "KR", "LA", "LI", "LR", "LU", "MD", "MG", "MN", "MO",
    "NA", "ND", "NE", "NI", "NO", "NP", "OS", "PA", "PB", "PD", "PM", "PO",
    "PR", "PT", "PU", "RA", "RB", "RE", "RH", "RN", "RU", "SB", "SC", "SE",
    "SI", "SM", "SN", "SR", "TA", "TB", "TC", "TE", "TH", "TI", "TL", "TM",
    "XE", "YB", "ZN", "ZR",
}
PROTEIN_FIRST_LETTER_ELEMENTS = {"H", "C", "N", "O", "P", "S"}


def clean_cif_missing(value: str | None) -> str | None:
    if value is None or value in {"", ".", "?"}:
        return None
    return value


def residue_label(chain: str, residue_id: str, insertion_code: str | None = None) -> str:
    return f"{chain}:{residue_id}{insertion_code or ''}"


def clean_element(value: str) -> str:
    letters = "".join(ch for ch in value.strip().upper() if ch.isalpha())
    return letters[:2]


def infer_element(record: str, atom_name: str, explicit_element: str = "") -> str:
    element = clean_element(explicit_element)
    if element:
        return element
    letters = clean_element(atom_name)
    if not letters:
        return ""
    if record == "ATOM" and letters[0] in PROTEIN_FIRST_LETTER_ELEMENTS:
        return letters[0]
    if len(letters) >= 2 and letters[:2] in TWO_LETTER_ELEMENTS:
        return letters[:2]
    return letters[0]


def atom_payload(
    *,
    record: str,
    line_number: int,
    serial: str,
    atom_name: str,
    resname: str,
    chain: str,
    residue_id: str,
    insertion_code: str | None,
    x: float,
    y: float,
    z: float,
    element: str,
) -> dict:
    resname = resname.upper() or "?"
    chain = chain or "?"
    residue_id = residue_id or "?"
    insertion_code = insertion_code or ""
    return {
        "record": record,
        "line_number": line_number,
        "serial": serial,
        "atom_name": atom_name or "?",
        "resname": resname,
        "ligand": resname,
        "chain": chain,
        "residue_id": residue_id,
        "insertion_code": insertion_code,
        "residue_label": residue_label(chain, residue_id, insertion_code),
        "x": x,
        "y": y,
        "z": z,
        "element": element,
    }


def _empty_result(fmt: str) -> dict:
    return {
        "format": fmt,
        "atom_records": 0,
        "hetatm_records": 0,
        "water_hetatm_records": 0,
        "malformed_records": 0,
        "protein_atoms": [],
        "ligand_atoms": [],
    }


def _add_atom(result: dict, atom: dict, ligand_filters: set[str], water_names: set[str]) -> None:
    if atom["record"] == "ATOM":
        result["protein_atoms"].append(atom)
        return
    if atom["resname"] in water_names:
        result["water_hetatm_records"] += 1
        return
    if ligand_filters and atom["resname"] not in ligand_filters:
        return
    result["ligand_atoms"].append(atom)


def parse_pdb_atoms(path: str | Path, ligand_filters: list[str] | None = None,
                    water_names: set[str] | None = None) -> dict:
    filters = set(ligand_filters or [])
    waters = water_names or {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}
    result = _empty_result("pdb")
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            if record == "ATOM":
                result["atom_records"] += 1
            else:
                result["hetatm_records"] += 1
            if len(line) < 54:
                result["malformed_records"] += 1
                continue
            altloc = line[16].strip()
            if altloc not in ALTLOCS_TO_KEEP:
                continue
            try:
                atom = atom_payload(
                    record=record,
                    line_number=line_number,
                    serial=line[6:11].strip(),
                    atom_name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21].strip(),
                    residue_id=line[22:26].strip(),
                    insertion_code=line[26].strip(),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    element=infer_element(record, line[12:16].strip(), line[76:78] if len(line) >= 78 else ""),
                )
            except ValueError:
                result["malformed_records"] += 1
                continue
            _add_atom(result, atom, filters, waters)
    return result


def _field_index(fields: dict[str, int], *names: str) -> int | None:
    for name in names:
        if name in fields:
            return fields[name]
    return None


def _row_value(row: list[str], index: int | None, default: str | None = None) -> str | None:
    if index is None or index >= len(row):
        return default
    return row[index]


def parse_mmcif_atoms(path: str | Path, ligand_filters: list[str] | None = None,
                      water_names: set[str] | None = None) -> dict:
    filters = set(ligand_filters or [])
    waters = water_names or {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}
    result = _empty_result("mmcif")
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()

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
            continue

        fields = {header: pos for pos, header in enumerate(headers)}
        group_i = _field_index(fields, "_atom_site.group_PDB")
        serial_i = _field_index(fields, "_atom_site.id")
        atom_i = _field_index(fields, "_atom_site.label_atom_id", "_atom_site.auth_atom_id")
        comp_i = _field_index(fields, "_atom_site.label_comp_id", "_atom_site.auth_comp_id")
        chain_i = _field_index(fields, "_atom_site.auth_asym_id", "_atom_site.label_asym_id")
        seq_i = _field_index(fields, "_atom_site.auth_seq_id", "_atom_site.label_seq_id")
        ins_i = _field_index(fields, "_atom_site.pdbx_PDB_ins_code")
        alt_i = _field_index(fields, "_atom_site.label_alt_id")
        element_i = _field_index(fields, "_atom_site.type_symbol")
        x_i = _field_index(fields, "_atom_site.Cartn_x")
        y_i = _field_index(fields, "_atom_site.Cartn_y")
        z_i = _field_index(fields, "_atom_site.Cartn_z")

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
                result["malformed_records"] += 1
                index += 1
                continue
            line_number = index + 1
            index += 1
            if len(row) < len(headers):
                result["malformed_records"] += 1
                continue

            record = (_row_value(row, group_i, "ATOM") or "ATOM").upper()
            if record not in {"ATOM", "HETATM"}:
                continue
            if record == "ATOM":
                result["atom_records"] += 1
            else:
                result["hetatm_records"] += 1
            altloc = clean_cif_missing(_row_value(row, alt_i)) or ""
            if altloc not in ALTLOCS_TO_KEEP:
                continue
            try:
                atom_name = clean_cif_missing(_row_value(row, atom_i)) or "?"
                atom = atom_payload(
                    record=record,
                    line_number=line_number,
                    serial=clean_cif_missing(_row_value(row, serial_i)) or "",
                    atom_name=atom_name,
                    resname=clean_cif_missing(_row_value(row, comp_i)) or "?",
                    chain=clean_cif_missing(_row_value(row, chain_i)) or "?",
                    residue_id=clean_cif_missing(_row_value(row, seq_i)) or "?",
                    insertion_code=clean_cif_missing(_row_value(row, ins_i)) or "",
                    x=float(_row_value(row, x_i) or ""),
                    y=float(_row_value(row, y_i) or ""),
                    z=float(_row_value(row, z_i) or ""),
                    element=infer_element(record, atom_name, clean_cif_missing(_row_value(row, element_i)) or ""),
                )
            except ValueError:
                result["malformed_records"] += 1
                continue
            _add_atom(result, atom, filters, waters)
    return result


def parse_structure_atoms(path: str | Path, ligand_filters: list[str] | None = None,
                          water_names: set[str] | None = None) -> dict:
    suffix = Path(path).suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return parse_mmcif_atoms(path, ligand_filters=ligand_filters, water_names=water_names)
    return parse_pdb_atoms(path, ligand_filters=ligand_filters, water_names=water_names)
