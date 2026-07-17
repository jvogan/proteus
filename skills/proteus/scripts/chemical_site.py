#!/usr/bin/env python3
"""Triage ligands, metals, waters, and covalent evidence around a chemical site."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteus_common
import structure_qc
import visual_common


SPEC_RE = re.compile(r"^[A-Za-z0-9_.+?-]+$")


def _distance(left: dict[str, Any], right: dict[str, Any]) -> float:
    return math.sqrt(sum((left[axis] - right[axis]) ** 2 for axis in ("x", "y", "z")))


def _groups(atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for atom in atoms:
        key = (atom["model"], atom["chain"], atom["residue_id"], atom["insertion_code"], atom["resname"])
        grouped[key].append(atom)
    output = []
    for key, values in grouped.items():
        model, chain, residue_id, insertion_code, resname = key
        output.append({
            "model": model,
            "chain": chain,
            "residue": f"{residue_id}{insertion_code}",
            "residue_id": residue_id,
            "insertion_code": insertion_code,
            "component": resname,
            "role": structure_qc.component_role(resname, values[0]["record"], len(values)),
            "atoms": values,
        })
    return output


def _choose_site(groups: list[dict[str, Any]], component: str | None) -> dict[str, Any]:
    chemical = [item for item in groups if item["role"] not in {"polymer", "modified_polymer", "water", "additive_or_solvent"}]
    if component:
        parts = component.split(":")
        matches = []
        for item in chemical:
            if len(parts) == 1 and item["component"].upper() == parts[0].upper():
                matches.append(item)
            elif len(parts) == 2 and item["chain"] == parts[0] and item["residue"] == parts[1]:
                matches.append(item)
            elif len(parts) == 3 and item["component"].upper() == parts[0].upper() and item["chain"] == parts[1] and item["residue"] == parts[2]:
                matches.append(item)
        if len(matches) != 1:
            raise visual_common.VisualWorkflowError(f"--component matched {len(matches)} sites; use RESN:CHAIN:RESIDUE to select one.")
        return matches[0]
    preferred = [item for item in chemical if item["role"] in {"ligand_or_unknown", "cofactor"}]
    candidates = preferred or chemical
    if not candidates:
        raise visual_common.VisualWorkflowError("No ligand, cofactor, or ion site was found.")
    if len(candidates) > 1:
        labels = ", ".join(f"{item['component']}:{item['chain']}:{item['residue']}" for item in candidates[:8])
        raise visual_common.VisualWorkflowError(f"Multiple chemical sites are present ({labels}); choose one with --component.")
    return candidates[0]


def inspect_chemical_site(structure: str, outdir: str, *, component: str | None = None,
                          radius: float = 5.0, execute: bool = False,
                          width: int = 1300, height: int = 900) -> dict[str, Any]:
    source = Path(structure).expanduser()
    if not source.is_file():
        raise visual_common.VisualWorkflowError(f"Structure not found: {proteus_common.display_path(source)}")
    if radius <= 0:
        raise visual_common.VisualWorkflowError("Radius must be positive.")
    text = structure_qc._read_text(source)
    fmt = structure_qc._format_from_path(source)
    raw, parser_meta = structure_qc.parse_mmcif(text) if fmt == "mmcif" else structure_qc.parse_pdb(text)
    atoms, selection = structure_qc.select_atoms(raw, model="first", altloc="highest")
    groups = _groups(atoms)
    site = _choose_site(groups, component)
    neighbors = []
    coordination = []
    for group in groups:
        if group is site:
            continue
        distances = [
            (_distance(left, right), left, right)
            for left in site["atoms"] for right in group["atoms"]
            if left["element"] != "H" and right["element"] != "H"
        ]
        if not distances:
            continue
        closest = min(distances, key=lambda item: item[0])
        distance, site_atom, neighbor_atom = closest
        if distance <= radius:
            record = {
                "component": group["component"],
                "chain": group["chain"],
                "residue": group["residue"],
                "role": group["role"],
                "distance_angstrom": round(distance, 3),
                "site_atom": site_atom["atom_name"],
                "neighbor_atom": neighbor_atom["atom_name"],
            }
            neighbors.append(record)
            if site["role"] == "ion" and distance <= 3.2:
                coordination.append(record)
    neighbors.sort(key=lambda item: item["distance_angstrom"])
    coordination.sort(key=lambda item: item["distance_angstrom"])
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    chain, residue, resname = site["chain"], site["residue"], site["component"]
    safe_visual = all(SPEC_RE.fullmatch(str(value)) for value in (chain, residue, resname))
    warnings = [
        "Distance-based coordination and covalent-link evidence are hypotheses; verify geometry, valence, protonation, and experimental density.",
        "Waters, additives, and ions can be modeled ambiguously at limited resolution.",
    ]
    pml: list[str] | None = None
    cxc: list[str] | None = None
    if safe_visual:
        pml = visual_common.pymol_base(width=width, height=height)
        site_selection = f"resn {resname} and chain {chain} and resi {residue}"
        pml.extend([
            f"load {visual_common.quote_pymol(source.resolve())}, structure",
            "hide everything, all",
            "show cartoon, polymer",
            "color gray80, polymer",
            f"select chemical_site, {site_selection}",
            f"select site_neighbors, byres (all within {radius:g} of chemical_site)",
            "show sticks, chemical_site or site_neighbors",
            "show spheres, chemical_site and inorganic",
            "color tv_orange, chemical_site",
            "color cyan, site_neighbors and polymer",
            "distance site_contacts, chemical_site, site_neighbors, 3.6, 2",
            "label chemical_site and not hydro, resn + resi",
            "orient chemical_site or site_neighbors",
            "zoom chemical_site or site_neighbors, 4",
            *visual_common.finalize_pymol(
                destination / "chemical_site.png", width=width, height=height,
                session=destination / "chemical_site.pse",
            ),
        ])
        cx_site = f"#1/{chain}:{residue}"
        cxc = [
            *visual_common.chimerax_base(),
            f"open {visual_common.quote_chimerax(source.resolve())}",
            "cartoon #1",
            "color #1 lightgray",
            f"select {cx_site}",
            f"select zone sel {radius:g} residues true",
            "show sel atoms",
            "style sel stick",
            f"color {cx_site} orange",
            f"contacts {cx_site} restrict sel reveal true showDist true log true",
            f"hbonds {cx_site} restrict sel reveal true log true",
            "view sel",
            *visual_common.finalize_chimerax(
                destination / "chemical_site_chimerax.png", width=width, height=height,
                session=destination / "chemical_site.cxs",
            ),
        ]
    else:
        warnings.append("The selected site's identifiers were not safe to embed in replay scripts, so only the analysis report was written.")
    site_data = {key: site[key] for key in ("model", "chain", "residue", "component", "role")}
    data: dict[str, Any] = {
        "workflow": "chemical_site",
        "structure": proteus_common.display_path(source),
        "site": site_data,
        "radius_angstrom": radius,
        "neighbors": neighbors,
        "coordination_candidates": coordination,
        "connectivity_evidence": parser_meta,
        "selection": selection,
        "executed": execute,
    }
    report = proteus_common.ok_payload(
        data, warnings=warnings, provenance={"structure": proteus_common.file_provenance(source)},
    )
    data["artifacts"] = visual_common.write_workflow(
        destination, "chemical_site", report=report, pymol_lines=pml, chimerax_lines=cxc,
    )
    if execute and pml is not None and cxc is not None:
        data["execution"] = {
            "pymol": visual_common.run_pymol(destination / "chemical_site.pml"),
            "chimerax": visual_common.run_chimerax(cxc),
        }
        failures = [name for name, result in data["execution"].items() if result.get("status") not in {"ok", "unavailable"}]
        if failures:
            report["status"] = "error"
            report["error"] = f"Chemical-site rendering failed for: {', '.join(failures)}"
        proteus_common.write_json(destination / "chemical_site.json", report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Inspect a ligand, cofactor, or metal site with chemical-context warnings.")
    parser.add_argument("structure")
    parser.add_argument("--component", help="RESN, CHAIN:RESIDUE, or RESN:CHAIN:RESIDUE")
    parser.add_argument("--radius", type=float, default=5.0)
    parser.add_argument("--outdir", default="proteus_chemical_site")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--json", action="store_true", help="Emit JSON (accepted for consistency)")
    args = parser.parse_args(argv)
    try:
        report = inspect_chemical_site(
            args.structure, args.outdir, component=args.component, radius=args.radius, execute=args.execute,
        )
    except (OSError, ValueError, structure_qc.StructureQCError, visual_common.VisualWorkflowError) as exc:
        print(json.dumps(proteus_common.error_payload(str(exc)), indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0 if report.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
