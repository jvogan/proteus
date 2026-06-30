#!/usr/bin/env python3
"""Parse simple protein substitutions and map them to structure coordinates.

This helper is intentionally conservative: it supports amino-acid substitutions
such as R175H and p.Arg175His, then optionally checks a local PDB/mmCIF file for
the residue and CA coordinate. Without a local structure, a UniProt accession can
be resolved through AlphaFold DB using the existing stdlib fetch helper.

Usage:
    python3 variant_map.py "P04637 R175H" --no-download --json
    python3 variant_map.py P04637:p.Arg175His --structure AF-P04637-F1.pdb --json
    python3 variant_map.py R175H --uniprot P04637 --structure AF-P04637-F1.cif --chain A --json
"""

import argparse
import json
import re
import shlex
import sys
from pathlib import Path

import fetch_alphafold
import sifts_map
import structure_info


ROOT = Path(__file__).resolve().parents[1]
AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
AA1_TO_3 = {one: three for three, one in AA3_TO_1.items()}
AA3_TITLE = {three: three.title() for three in AA3_TO_1}

UNIPROT_RE = re.compile(
    r"^(?:[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9])(?:-[0-9]+)?$",
    re.IGNORECASE,
)
ONE_LETTER_RE = re.compile(r"^(?:p\.)?([A-Za-z])([1-9][0-9]*)([A-Za-z])$")
THREE_LETTER_RE = re.compile(r"^(?:p\.)?([A-Za-z]{3})([1-9][0-9]*)([A-Za-z]{3})$")

SIFTS_WARNING = (
    "SIFTS/PDBe UniProt-to-PDB residue mapping was not applied. For experimental "
    "PDB/mmCIF structures, residue_index is matched directly against structure "
    "residue numbering; use SIFTS mapping before treating this as authoritative."
)


class VariantMapError(ValueError):
    """Raised when a variant cannot be parsed or mapped conservatively."""


def _error_payload(message: str) -> dict:
    return {"status": "error", "error": message}


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _looks_like_uniprot(value: str) -> bool:
    return bool(UNIPROT_RE.fullmatch(value.strip()))


def _display_path(path: str | Path | None) -> str | None:
    if path is None:
        return None
    candidate = Path(path).expanduser()
    try:
        resolved = candidate.resolve(strict=False)
    except OSError:
        resolved = candidate.absolute()
    for root in (ROOT, Path.cwd()):
        try:
            relative = resolved.relative_to(root.resolve())
            return f"./{relative}" if str(relative) != "." else "."
        except ValueError:
            continue
    return f"{resolved.name} (absolute path omitted)"


def _normalize_uniprot(value: str) -> str:
    normalized = value.strip().upper()
    if not _looks_like_uniprot(normalized):
        raise VariantMapError(f"Invalid UniProt accession: {value}")
    return normalized


def _aa_payload(one_letter: str) -> dict:
    three = AA1_TO_3[one_letter]
    return {
        "one_letter": one_letter,
        "three_letter": three,
        "three_letter_title": AA3_TITLE[three],
    }


def _parse_substitution(value: str) -> dict:
    text = value.strip()
    match = THREE_LETTER_RE.fullmatch(text)
    if match:
        from3, index, to3 = match.groups()
        from3 = from3.upper()
        to3 = to3.upper()
        if from3 not in AA3_TO_1 or to3 not in AA3_TO_1:
            raise VariantMapError(f"Unknown three-letter amino-acid code in variant: {value}")
        from1 = AA3_TO_1[from3]
        to1 = AA3_TO_1[to3]
    else:
        match = ONE_LETTER_RE.fullmatch(text)
        if not match:
            raise VariantMapError(
                "Expected a simple protein substitution such as R175H or p.Arg175His."
            )
        from1, index, to1 = match.groups()
        from1 = from1.upper()
        to1 = to1.upper()
        if from1 not in AA1_TO_3 or to1 not in AA1_TO_3:
            raise VariantMapError(f"Unknown one-letter amino-acid code in variant: {value}")
        from3 = AA1_TO_3[from1]
        to3 = AA1_TO_3[to1]

    residue_index = int(index)
    return {
        "raw": value,
        "kind": "protein_substitution",
        "short": f"{from1}{residue_index}{to1}",
        "hgvs_protein": f"p.{AA3_TITLE[from3]}{residue_index}{AA3_TITLE[to3]}",
        "from": _aa_payload(from1),
        "to": _aa_payload(to1),
        "residue_index": residue_index,
    }


def parse_variant(value: str, uniprot_id: str | None = None) -> dict:
    """Parse a simple protein substitution with an optional UniProt prefix."""

    raw = " ".join(value.strip().split())
    if not raw:
        raise VariantMapError("Variant input is empty.")

    parsed_uniprot = None
    variant_text = raw

    if ":" in raw:
        left, right = raw.split(":", 1)
        if _looks_like_uniprot(left):
            parsed_uniprot = _normalize_uniprot(left)
            variant_text = right.strip()
    else:
        parts = raw.split()
        if len(parts) == 2 and _looks_like_uniprot(parts[0]):
            parsed_uniprot = _normalize_uniprot(parts[0])
            variant_text = parts[1]
        elif len(parts) > 1:
            raise VariantMapError(
                "Expected '<UniProt> <substitution>' or a single substitution plus --uniprot."
            )

    if uniprot_id:
        cli_uniprot = _normalize_uniprot(uniprot_id)
        if parsed_uniprot and parsed_uniprot != cli_uniprot:
            raise VariantMapError(
                f"UniProt accession mismatch: input has {parsed_uniprot}, --uniprot has {cli_uniprot}."
            )
        parsed_uniprot = cli_uniprot

    parsed = _parse_substitution(variant_text)
    parsed["input"] = raw
    parsed["uniprot_id"] = parsed_uniprot
    return parsed


def _coord_payload(x: float, y: float, z: float) -> dict:
    return {"x": round(x, 3), "y": round(y, 3), "z": round(z, 3)}


def _clean_cif_missing(value: str | None) -> str | None:
    if value is None or value in {".", "?"}:
        return None
    return value


def _residue_entry(chain: str, seq_id: str, insertion_code: str | None, resname: str) -> dict:
    resname = resname.upper()
    return {
        "chain": chain,
        "residue_id": seq_id,
        "insertion_code": insertion_code,
        "resname": resname,
        "residue_one_letter": AA3_TO_1.get(resname),
        "atom_count": 0,
        "ca_coordinate": None,
    }


def _parse_pdb_residues(path: str) -> list[dict]:
    residues: dict[tuple[str, str, str | None, str], dict] = {}
    with open(path) as handle:
        for line in handle:
            record = line[:6].strip()
            if record != "ATOM" or len(line) < 54:
                continue
            atom = line[12:16].strip()
            altloc = line[16].strip()
            resname = line[17:20].strip().upper()
            chain = line[21].strip() or "?"
            seq_id = line[22:26].strip()
            insertion_code = line[26].strip() or None
            key = (chain, seq_id, insertion_code, resname)
            entry = residues.setdefault(key, _residue_entry(chain, seq_id, insertion_code, resname))
            entry["atom_count"] += 1
            if atom == "CA" and entry["ca_coordinate"] is None and altloc in {"", "A", "1"}:
                try:
                    entry["ca_coordinate"] = _coord_payload(
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                except ValueError:
                    pass
    return sorted(residues.values(), key=lambda item: (item["chain"], item["residue_id"], item["resname"]))


def _field_index(field_index: dict[str, int], *names: str) -> int | None:
    for name in names:
        if name in field_index:
            return field_index[name]
    return None


def _row_value(row: list[str], index: int | None, default: str | None = None) -> str | None:
    if index is None or index >= len(row):
        return default
    return row[index]


def _parse_mmcif_residues(path: str) -> list[dict]:
    lines = Path(path).read_text(errors="replace").splitlines()
    residues: dict[tuple[str, str, str | None, str], dict] = {}

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
            while index < len(lines) and lines[index].strip() and not lines[index].strip().startswith(("loop_", "_", "#")):
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

            group = (_row_value(row, group_i, "ATOM") or "ATOM").upper()
            if group != "ATOM":
                continue
            atom = _row_value(row, atom_i)
            resname = (_row_value(row, comp_i, "?") or "?").upper()
            chain = _clean_cif_missing(_row_value(row, chain_i)) or "?"
            seq_id = _clean_cif_missing(_row_value(row, seq_i)) or "?"
            insertion_code = _clean_cif_missing(_row_value(row, ins_i))
            altloc = _clean_cif_missing(_row_value(row, alt_i)) or ""
            key = (chain, seq_id, insertion_code, resname)
            entry = residues.setdefault(key, _residue_entry(chain, seq_id, insertion_code, resname))
            entry["atom_count"] += 1
            if atom == "CA" and entry["ca_coordinate"] is None and altloc in {"", "A", "1"}:
                try:
                    x = _row_value(row, x_i)
                    y = _row_value(row, y_i)
                    z = _row_value(row, z_i)
                    if x is None or y is None or z is None:
                        raise ValueError("missing coordinate")
                    entry["ca_coordinate"] = _coord_payload(
                        float(x),
                        float(y),
                        float(z),
                    )
                except ValueError:
                    pass

    return sorted(residues.values(), key=lambda item: (item["chain"], item["residue_id"], item["resname"]))


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


def _with_variant_match(candidate: dict, variant: dict) -> dict:
    item = dict(candidate)
    residue_one = item.get("residue_one_letter")
    item["reference_matches_variant"] = (
        residue_one == variant["from"]["one_letter"] if residue_one else None
    )
    return item


def find_local_residue(
    path: str,
    variant: dict,
    chain: str | None = None,
    sifts_candidates: list[dict] | None = None,
) -> dict:
    suffix = Path(path).suffix.lower()
    fmt = "mmcif" if suffix in {".cif", ".mmcif"} else "pdb"
    residues = _parse_mmcif_residues(path) if fmt == "mmcif" else _parse_pdb_residues(path)
    residue_index = variant["residue_index"]

    if sifts_candidates:
        candidates = []
        seen = set()
        for target in sifts_candidates:
            for item in residues:
                if not _residue_matches_sifts_target(item, target):
                    continue
                key = (item["chain"], item["residue_id"], item.get("insertion_code"), item["resname"])
                if key in seen:
                    continue
                seen.add(key)
                candidate = _with_variant_match(item, variant)
                candidate["sifts_mapping"] = sifts_map.public_mapping_candidate(target)
                candidates.append(candidate)
        lookup_method = "sifts_uniprot_to_structure_mapping"
    else:
        candidates = [
            _with_variant_match(item, variant)
            for item in residues
            if _seq_matches(item["residue_id"], residue_index)
            and (chain is None or item["chain"] == chain)
        ]
        lookup_method = "direct_structure_residue_numbering"

    selected = candidates[0] if len(candidates) == 1 else None
    return {
        "checked": True,
        "format": fmt,
        "lookup_method": lookup_method,
        "residue_index": residue_index,
        "chain_filter": chain,
        "sifts_mapping_applied": bool(sifts_candidates),
        "sifts_mapping_candidate_count": len(sifts_candidates or []),
        "residue_present": bool(candidates),
        "ca_present": any(item.get("ca_coordinate") is not None for item in candidates),
        "candidate_count": len(candidates),
        "selected": selected,
        "candidates": candidates,
    }


def _inspect_structure(path: str, force_alphafold: bool = False) -> dict:
    inspection = structure_info.inspect_structure(path, force_alphafold=force_alphafold)
    if inspection.get("file"):
        inspection["file"] = _display_path(inspection["file"])
    return inspection


def _alphafold_metadata_payload(meta: dict, uniprot_id: str) -> dict:
    return {
        "source": "AlphaFold DB",
        "uniprot_id": uniprot_id,
        "model_id": meta.get("modelEntityId", f"AF-{uniprot_id}-F1"),
        "gene": meta.get("gene"),
        "global_plddt": meta.get("globalMetricValue"),
        "latest_version": meta.get("latestVersion"),
        "sequence_start": meta.get("sequenceStart"),
        "sequence_end": meta.get("sequenceEnd"),
    }


def _resolve_alphafold(uniprot_id: str, outdir: str, download: bool) -> tuple[dict, str | None]:
    meta = fetch_alphafold.fetch_metadata(uniprot_id)
    payload = _alphafold_metadata_payload(meta, uniprot_id)
    payload.update({
        "resolved_kind": "uniprot_accession",
        "structure_path": None,
        "downloaded": False,
    })
    local_path = None
    if download:
        output_dir = Path(outdir)
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / f"{payload['model_id']}.pdb"
        if not destination.exists():
            fetch_alphafold.download(meta["pdbUrl"], destination, as_json=True)
            payload["downloaded"] = True
        local_path = str(destination.resolve())
        payload["structure_path"] = _display_path(destination)
    return payload, local_path


def map_variant(value: str, *, uniprot_id: str | None = None, structure: str | None = None,
                chain: str | None = None, outdir: str = ".", download: bool = True,
                sifts_json: str | None = None) -> dict:
    variant = parse_variant(value, uniprot_id=uniprot_id)
    warnings = []
    sifts_candidates = []
    sifts_mapping = {
        "applied": False,
        "status": "not_supplied",
        "warning": None,
    }
    structure_payload = None
    local_residue = {
        "checked": False,
        "reason": "No local coordinate file was available.",
    }

    if sifts_json:
        sifts_source = _display_path(sifts_json)
        if not variant.get("uniprot_id"):
            warning = (
                "SIFTS JSON was supplied, but no UniProt accession was available; "
                "direct structure residue numbering was used."
            )
            warnings.append(warning)
            sifts_mapping = {
                "applied": False,
                "status": "missing_uniprot",
                "source": sifts_source,
                "warning": warning,
            }
        else:
            records = sifts_map.load_sifts_json(sifts_json)
            sifts_candidates = sifts_map.map_uniprot_residue_candidates(
                records,
                variant["uniprot_id"],
                variant["residue_index"],
                chain_id=chain,
            )
            if sifts_candidates:
                sifts_mapping = {
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
                warning = (
                    "SIFTS JSON was supplied, but no mapping matched this UniProt "
                    "residue and chain filter; direct structure residue numbering was used."
                )
                warnings.append(warning)
                sifts_mapping = {
                    "applied": False,
                    "status": "no_match",
                    "source": sifts_source,
                    "warning": warning,
                }

    if structure:
        path = Path(structure)
        if not path.exists():
            raise VariantMapError(f"Structure file not found: {structure}")
        structure_path = str(path.resolve())
        inspection = _inspect_structure(structure_path)
        structure_payload = {
            "source": "local",
            "resolved_kind": "local_file",
            "structure_path": _display_path(structure_path),
            "inspection": inspection,
        }
        local_residue = find_local_residue(
            structure_path,
            variant,
            chain=chain,
            sifts_candidates=sifts_candidates or None,
        )
        if variant.get("uniprot_id") and not inspection.get("likely_alphafold") and not sifts_json:
            warnings.append(SIFTS_WARNING)
            sifts_mapping = {
                "applied": False,
                "status": "deferred",
                "warning": SIFTS_WARNING,
            }
    elif variant.get("uniprot_id"):
        structure_payload, local_structure_path = _resolve_alphafold(variant["uniprot_id"], outdir, download)
        if local_structure_path:
            inspection = _inspect_structure(local_structure_path, force_alphafold=True)
            structure_payload["inspection"] = inspection
            local_residue = find_local_residue(
                local_structure_path,
                variant,
                chain=chain,
                sifts_candidates=sifts_candidates or None,
            )
        else:
            local_residue = {
                "checked": False,
                "reason": "AlphaFold metadata was resolved, but no local coordinate file was downloaded.",
            }
            warnings.append("No local coordinates checked; pass --structure or omit --no-download.")
    else:
        raise VariantMapError(
            "Provide a UniProt accession in the variant, pass --uniprot, or supply --structure for local mapping."
        )

    if local_residue.get("candidate_count", 0) > 1 and chain is None:
        warnings.append("Multiple residues matched this index; pass --chain to disambiguate.")
    if local_residue.get("residue_present") and not local_residue.get("ca_present"):
        warnings.append("Residue was present, but no CA coordinate was found.")
    selected = local_residue.get("selected")
    if selected and selected.get("reference_matches_variant") is False:
        warnings.append(
            f"Structure residue {selected['resname']} does not match variant reference "
            f"{variant['from']['three_letter']}."
        )

    data = {
        "input": value,
        "variant": variant,
        "uniprot_id": variant.get("uniprot_id"),
        "residue_index": variant["residue_index"],
        "structure": structure_payload,
        "local_residue": local_residue,
        "sifts_mapping": sifts_mapping,
        "warnings": warnings,
    }
    return _ok_payload(data)


def main():
    parser = argparse.ArgumentParser(
        description="Map a simple protein substitution to UniProt/AlphaFold and optional local coordinates.",
        epilog=(
            "Examples:\n"
            "  %(prog)s 'P04637 R175H' --no-download --json\n"
            "  %(prog)s P04637:p.Arg175His --structure AF-P04637-F1.pdb --json\n"
            "  %(prog)s R175H --uniprot P04637 --structure AF-P04637-F1.cif --chain A --json\n"
            "  %(prog)s G175A --uniprot P04637 --structure model.pdb --sifts-json sifts.json --json\n"
            "\nNotes:\n"
            "  Local experimental PDB/mmCIF lookup falls back to direct structure\n"
            "  numbering unless --sifts-json supplies offline PDBe SIFTS mappings."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("variant", nargs="+", help="Substitution, optionally prefixed by UniProt accession")
    parser.add_argument("--uniprot", help="UniProt accession when the variant does not include one")
    parser.add_argument("--structure", help="Local .pdb, .cif, or .mmcif file to inspect")
    parser.add_argument("--chain", help="Optional chain ID for local residue lookup")
    parser.add_argument("--sifts-json", help="Local PDBe SIFTS-style JSON mapping file")
    parser.add_argument("--outdir", default=".", help="Output directory for downloaded AlphaFold files")
    parser.add_argument("--no-download", action="store_true", help="Resolve AlphaFold metadata only")
    parser.add_argument("--json", action="store_true", dest="as_json", help="Emit machine-readable JSON")
    args = parser.parse_args()

    variant_text = " ".join(args.variant)
    try:
        output = map_variant(
            variant_text,
            uniprot_id=args.uniprot,
            structure=args.structure,
            chain=args.chain,
            outdir=args.outdir,
            download=not args.no_download,
            sifts_json=args.sifts_json,
        )
    except (VariantMapError, OSError, KeyError, fetch_alphafold.AlphaFoldFetchError, sifts_map.SiftsLookupError) as exc:
        if args.as_json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.as_json:
        print(json.dumps(output, indent=2))
        return

    data = output["data"]
    variant = data["variant"]
    print(f"Variant: {variant['short']} ({variant['hgvs_protein']})")
    if data.get("uniprot_id"):
        print(f"UniProt: {data['uniprot_id']}")
    print(f"Residue index: {data['residue_index']}")
    structure_payload = data.get("structure") or {}
    if structure_payload.get("source"):
        print(f"Structure source: {structure_payload['source']}")
    if structure_payload.get("structure_path"):
        print(f"Structure: {structure_payload['structure_path']}")
    local = data["local_residue"]
    if data["sifts_mapping"].get("applied"):
        print(f"SIFTS mapping: applied ({data['sifts_mapping']['candidate_count']} candidate(s))")
    if local.get("checked"):
        print(f"Residue present: {'yes' if local['residue_present'] else 'no'}")
        print(f"CA present: {'yes' if local['ca_present'] else 'no'}")
        selected = local.get("selected")
        if selected and selected.get("ca_coordinate"):
            coord = selected["ca_coordinate"]
            print(f"CA coordinate: ({coord['x']}, {coord['y']}, {coord['z']})")
    for warning in data["warnings"]:
        print(f"WARNING: {warning}", file=sys.stderr)


if __name__ == "__main__":
    main()
