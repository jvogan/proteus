#!/usr/bin/env python3
"""Structure preflight and quality-control report for PDB/mmCIF inputs.

The parser is intentionally dependency-free and conservative. It exposes model,
alternate-conformer, occupancy, residue, and chemical-component decisions so
downstream workflows do not silently combine incompatible coordinate states.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import shlex
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common


STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
MODIFIED_POLYMER = {
    "MSE", "SEC", "PYL", "SEP", "TPO", "PTR", "CSO", "CSD", "CME", "MLY",
    "M3L", "HYP", "KCX", "LLP", "FME", "PCA", "DAL", "DLE", "DVA", "DTH",
}
NUCLEOTIDES = {
    "A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI", "ADE", "CYT",
    "GUA", "THY", "URA",
}
WATERS = {"HOH", "WAT", "DOD", "H2O", "TIP", "T3P", "SOL"}
IONS = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE", "CO",
    "NI", "CU", "ZN", "CD", "HG", "AL", "GA", "CL", "BR", "IOD", "F", "SO4",
    "PO4", "NO3", "NH4",
}
COFACTORS = {
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "FAD", "FMN", "NAD", "NAP",
    "NDP", "SAM", "SAH", "COA", "HEM", "HEC", "PLP", "TPP", "THF", "B12",
}
ADDITIVES = {
    "GOL", "EDO", "PEG", "PG4", "PGE", "MPD", "DMS", "BME", "MES", "TRS",
    "HEP", "ACT", "FMT", "ACE", "CIT", "TAR", "EOH", "IPA", "ACN", "BOG",
}
BACKBONE = {"N", "CA", "C", "O"}


class StructureQCError(RuntimeError):
    pass


def _format_from_path(path: Path) -> str:
    name = path.name.lower()
    if name.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        return "mmcif"
    if name.endswith((".pdb", ".ent", ".pdb.gz", ".ent.gz")):
        return "pdb"
    raise StructureQCError("Input must be PDB/mmCIF, optionally gzip-compressed.")


def _read_text(path: Path) -> str:
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8", errors="replace")


def _float(value: str | None) -> float | None:
    if value is None or value in {"", ".", "?"}:
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def _clean(value: str | None, default: str = "") -> str:
    if value is None or value in {"", ".", "?"}:
        return default
    return value


def _atom(
    *, record: str, model: int, serial: str, atom_name: str, altloc: str,
    resname: str, chain: str, residue_id: str, insertion_code: str,
    x: float, y: float, z: float, occupancy: float | None,
    bfactor: float | None, element: str,
) -> dict[str, Any]:
    return {
        "record": record,
        "model": model,
        "serial": serial,
        "atom_name": atom_name,
        "altloc": altloc,
        "resname": resname.upper(),
        "chain": chain or "?",
        "residue_id": residue_id or "?",
        "insertion_code": insertion_code,
        "x": x,
        "y": y,
        "z": z,
        "occupancy": occupancy,
        "bfactor": bfactor,
        "element": element.upper(),
    }


def parse_pdb(text: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    current_model = 1
    explicit_models = False
    malformed = 0
    link_records = 0
    conect_records = 0
    for line in text.splitlines():
        record = line[:6].strip().upper()
        if record == "MODEL":
            explicit_models = True
            try:
                current_model = int(line[10:14].strip())
            except ValueError:
                current_model += 1
            continue
        if record == "LINK":
            link_records += 1
            continue
        if record == "CONECT":
            conect_records += 1
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(line) < 54:
            malformed += 1
            continue
        try:
            atoms.append(_atom(
                record=record,
                model=current_model,
                serial=line[6:11].strip(),
                atom_name=line[12:16].strip(),
                altloc=line[16].strip(),
                resname=line[17:20].strip(),
                chain=line[21].strip(),
                residue_id=line[22:26].strip(),
                insertion_code=line[26].strip(),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                occupancy=_float(line[54:60] if len(line) >= 60 else None),
                bfactor=_float(line[60:66] if len(line) >= 66 else None),
                element=line[76:78].strip() if len(line) >= 78 else "",
            ))
        except ValueError:
            malformed += 1
    return atoms, {
        "explicit_models": explicit_models,
        "malformed_atom_records": malformed,
        "link_records": link_records,
        "conect_records": conect_records,
    }


def _field(fields: dict[str, int], *names: str) -> int | None:
    for name in names:
        if name in fields:
            return fields[name]
    return None


def _value(row: list[str], index: int | None, default: str = "") -> str:
    if index is None or index >= len(row):
        return default
    return row[index]


def parse_mmcif(text: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    lines = text.splitlines()
    atoms: list[dict[str, Any]] = []
    malformed = 0
    struct_conn_rows = 0
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        headers: list[str] = []
        while index < len(lines) and lines[index].strip().startswith("_"):
            headers.append(lines[index].strip())
            index += 1
        if not headers:
            continue
        atom_loop = all(item.startswith("_atom_site.") for item in headers)
        conn_loop = all(item.startswith("_struct_conn.") for item in headers)
        fields = {name: position for position, name in enumerate(headers)}
        if atom_loop:
            group_i = _field(fields, "_atom_site.group_PDB")
            model_i = _field(fields, "_atom_site.pdbx_PDB_model_num")
            serial_i = _field(fields, "_atom_site.id")
            atom_i = _field(fields, "_atom_site.auth_atom_id", "_atom_site.label_atom_id")
            alt_i = _field(fields, "_atom_site.label_alt_id", "_atom_site.auth_alt_id")
            comp_i = _field(fields, "_atom_site.auth_comp_id", "_atom_site.label_comp_id")
            chain_i = _field(fields, "_atom_site.auth_asym_id", "_atom_site.label_asym_id")
            seq_i = _field(fields, "_atom_site.auth_seq_id", "_atom_site.label_seq_id")
            ins_i = _field(fields, "_atom_site.pdbx_PDB_ins_code")
            x_i = _field(fields, "_atom_site.Cartn_x")
            y_i = _field(fields, "_atom_site.Cartn_y")
            z_i = _field(fields, "_atom_site.Cartn_z")
            occupancy_i = _field(fields, "_atom_site.occupancy")
            bfactor_i = _field(fields, "_atom_site.B_iso_or_equiv")
            element_i = _field(fields, "_atom_site.type_symbol")

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
                malformed += 1
                index += 1
                continue
            index += 1
            if len(row) < len(headers):
                malformed += 1
                continue
            if conn_loop:
                struct_conn_rows += 1
                continue
            if not atom_loop:
                continue
            record = _clean(_value(row, group_i), "ATOM").upper()
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                atoms.append(_atom(
                    record=record,
                    model=int(_clean(_value(row, model_i), "1")),
                    serial=_clean(_value(row, serial_i)),
                    atom_name=_clean(_value(row, atom_i), "?"),
                    altloc=_clean(_value(row, alt_i)),
                    resname=_clean(_value(row, comp_i), "?"),
                    chain=_clean(_value(row, chain_i), "?"),
                    residue_id=_clean(_value(row, seq_i), "?"),
                    insertion_code=_clean(_value(row, ins_i)),
                    x=float(_value(row, x_i)),
                    y=float(_value(row, y_i)),
                    z=float(_value(row, z_i)),
                    occupancy=_float(_value(row, occupancy_i)),
                    bfactor=_float(_value(row, bfactor_i)),
                    element=_clean(_value(row, element_i)),
                ))
            except (TypeError, ValueError):
                malformed += 1
    return atoms, {
        "explicit_models": len({atom["model"] for atom in atoms}) > 1,
        "malformed_atom_records": malformed,
        "struct_conn_rows": struct_conn_rows,
    }


def _altloc_rank(atom: dict[str, Any]) -> tuple[float, int, str]:
    occupancy = atom["occupancy"] if atom["occupancy"] is not None else -1.0
    preferred = {"": 3, "A": 2, "1": 1}.get(atom["altloc"], 0)
    return occupancy, preferred, atom["altloc"]


def select_atoms(
    atoms: list[dict[str, Any]], *, model: str = "first", altloc: str = "highest",
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    models = sorted({atom["model"] for atom in atoms})
    if not models:
        return [], {"available_models": [], "selected_models": [], "alternate_sites": 0}
    if model == "all":
        selected_models = models
    elif model == "first":
        selected_models = [models[0]]
    else:
        try:
            requested = int(model)
        except ValueError as exc:
            raise StructureQCError("--model must be 'first', 'all', or an integer.") from exc
        if requested not in models:
            raise StructureQCError(f"Model {requested} is not present; available models: {models}")
        selected_models = [requested]
    candidates = [atom for atom in atoms if atom["model"] in selected_models]
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for atom in candidates:
        key = (
            atom["model"], atom["chain"], atom["residue_id"], atom["insertion_code"],
            atom["resname"], atom["atom_name"],
        )
        grouped[key].append(atom)
    alternate_sites = sum(1 for values in grouped.values() if len({item["altloc"] for item in values}) > 1)
    if altloc == "all":
        selected = candidates
    elif altloc == "highest":
        selected = [max(values, key=_altloc_rank) for values in grouped.values()]
    else:
        requested = altloc.strip()
        selected = []
        for values in grouped.values():
            exact = [item for item in values if item["altloc"] in {"", requested}]
            selected.append(max(exact or values, key=_altloc_rank))
    return selected, {
        "available_models": models,
        "selected_models": selected_models,
        "alternate_sites": alternate_sites,
        "altloc_policy": altloc,
    }


def component_role(resname: str, record: str, atom_count: int) -> str:
    if record == "ATOM" or resname in STANDARD_AMINO_ACIDS or resname in NUCLEOTIDES:
        return "polymer"
    if resname in MODIFIED_POLYMER:
        return "modified_polymer"
    if resname in WATERS:
        return "water"
    if resname in IONS:
        return "ion"
    if resname in COFACTORS:
        return "cofactor"
    if resname in ADDITIVES:
        return "additive_or_solvent"
    return "ligand_or_unknown"


def _residue_sort(value: str) -> tuple[int, int | str]:
    try:
        return 0, int(value)
    except ValueError:
        return 1, value


def build_report(path: str | Path, *, model: str = "first", altloc: str = "highest") -> dict[str, Any]:
    source = Path(path).expanduser()
    if not source.is_file():
        raise StructureQCError(f"File not found: {proteus_common.display_path(source)}")
    fmt = _format_from_path(source)
    text = _read_text(source)
    raw_atoms, parser_meta = parse_mmcif(text) if fmt == "mmcif" else parse_pdb(text)
    atoms, selection = select_atoms(raw_atoms, model=model, altloc=altloc)
    warnings: list[str] = []
    issues: list[dict[str, Any]] = []

    if len(selection["available_models"]) > 1:
        warnings.append("Multiple coordinate models are present; downstream geometry should use an explicit model.")
        issues.append({"code": "multiple_models", "count": len(selection["available_models"])})
    if selection["alternate_sites"]:
        warnings.append("Alternate conformers are present; the report records the conformer-selection policy.")
        issues.append({"code": "alternate_conformers", "count": selection["alternate_sites"]})
    if parser_meta.get("malformed_atom_records"):
        warnings.append("Some coordinate records could not be parsed by the lightweight parser.")
        issues.append({"code": "malformed_atom_records", "count": parser_meta["malformed_atom_records"]})

    residue_atoms: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    component_atoms: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for atom in atoms:
        residue_key = (
            atom["model"], atom["chain"], atom["residue_id"], atom["insertion_code"], atom["resname"],
        )
        residue_atoms[residue_key].append(atom)
        if atom["record"] == "HETATM":
            component_atoms[residue_key].append(atom)

    missing_backbone: list[dict[str, Any]] = []
    chain_residues: dict[tuple[int, str], list[tuple[str, str, str]]] = defaultdict(list)
    for (model_id, chain, residue_id, insertion_code, resname), values in residue_atoms.items():
        if values[0]["record"] == "ATOM" or resname in STANDARD_AMINO_ACIDS or resname in MODIFIED_POLYMER:
            chain_residues[(model_id, chain)].append((residue_id, insertion_code, resname))
        if resname in STANDARD_AMINO_ACIDS:
            present = {item["atom_name"] for item in values}
            missing = sorted(BACKBONE - present)
            if missing:
                missing_backbone.append({
                    "model": model_id,
                    "chain": chain,
                    "residue": f"{residue_id}{insertion_code}",
                    "resname": resname,
                    "missing": missing,
                })
    if missing_backbone:
        warnings.append("Standard amino-acid residues with missing backbone atoms were detected.")
        issues.append({"code": "missing_backbone_atoms", "count": len(missing_backbone)})

    chain_gaps: list[dict[str, Any]] = []
    for (model_id, chain), residues in sorted(chain_residues.items()):
        numeric = sorted({int(item[0]) for item in residues if item[0].lstrip("-").isdigit()})
        for left, right in zip(numeric, numeric[1:]):
            if right - left > 1:
                chain_gaps.append({"model": model_id, "chain": chain, "after": left, "before": right})
    if chain_gaps:
        warnings.append("Residue-number gaps were detected; they may represent missing density or numbering conventions.")

    components: list[dict[str, Any]] = []
    role_counts: Counter[str] = Counter()
    for (model_id, chain, residue_id, insertion_code, resname), values in sorted(
        component_atoms.items(), key=lambda item: (item[0][0], item[0][1], _residue_sort(item[0][2]), item[0][4]),
    ):
        role = component_role(resname, values[0]["record"], len(values))
        role_counts[role] += 1
        components.append({
            "model": model_id,
            "chain": chain,
            "residue": f"{residue_id}{insertion_code}",
            "component": resname,
            "role": role,
            "atom_count": len(values),
            "elements": sorted({item["element"] for item in values if item["element"]}),
        })

    occupancies = [item["occupancy"] for item in atoms if item["occupancy"] is not None]
    low_occupancy = sum(1 for value in occupancies if value < 0.5)
    invalid_occupancy = sum(1 for value in occupancies if value < 0 or value > 1.01)
    if low_occupancy:
        issues.append({"code": "low_occupancy_atoms", "count": low_occupancy})
    if invalid_occupancy:
        warnings.append("Occupancy values outside the expected 0-1 range were detected.")
        issues.append({"code": "invalid_occupancy_atoms", "count": invalid_occupancy})

    chains = sorted({item["chain"] for item in atoms})
    data = {
        "file": proteus_common.display_path(source),
        "format": fmt,
        "compressed": source.name.lower().endswith(".gz"),
        "selection": selection,
        "counts": {
            "raw_atoms": len(raw_atoms),
            "selected_atoms": len(atoms),
            "chains": len(chains),
            "residues": len(residue_atoms),
            "components": len(components),
        },
        "chains": chains,
        "component_role_counts": dict(sorted(role_counts.items())),
        "components": components,
        "occupancy": {
            "reported_atoms": len(occupancies),
            "min": min(occupancies) if occupancies else None,
            "max": max(occupancies) if occupancies else None,
            "low_below_0_5": low_occupancy,
        },
        "missing_backbone": missing_backbone,
        "numbering_gaps": chain_gaps,
        "issues": issues,
        "parser": parser_meta,
    }
    return proteus_common.ok_payload(
        data,
        warnings=warnings,
        provenance={"input": proteus_common.file_provenance(source)},
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Preflight PDB/mmCIF coordinate quality and selection decisions.")
    parser.add_argument("structure", help="PDB/mmCIF path, optionally .gz")
    parser.add_argument("--model", default="first", help="Coordinate model: first, all, or integer (default: first)")
    parser.add_argument("--altloc", default="highest", help="Alternate conformer: highest, all, or ID (default: highest)")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency; JSON is always emitted)")
    args = parser.parse_args(argv)
    try:
        report = build_report(args.structure, model=args.model, altloc=args.altloc)
    except (OSError, StructureQCError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
