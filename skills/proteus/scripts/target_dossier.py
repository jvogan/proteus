#!/usr/bin/env python3
"""Build a generic target dossier from genes, UniProt accessions, and PDB inputs.

The script is intentionally conservative: it writes a Markdown summary plus a
machine-readable provenance JSON file, and it keeps network work behind the
default online mode. Use --no-network for deterministic offline runs.

Examples:
    python3 scripts/target_dossier.py --gene KRAS --uniprot P01116 --pdb 6OIM
    python3 scripts/target_dossier.py --pdb tests/fixtures/tiny.pdb --no-network
    python3 scripts/target_dossier.py --pdb model.cif --no-network --analyze-local
    python3 scripts/target_dossier.py --gene TP53 --out tp53_dossier --json
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
PYTHON = sys.executable
REPORT_NAME = "TARGET_DOSSIER.md"
PROVENANCE_NAME = "provenance.json"
MAX_CAPTURE_CHARS = 20000
LOCAL_ANALYSIS_HELPERS = (
    ("ligands", "ligand_extract.py", []),
    ("interfaces", "interface_report.py", []),
    ("interactions", "interaction_report.py", []),
    ("docking_box", "docking_box.py", []),
)

PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")
UNIPROT_RE = re.compile(
    r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$|"
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$"
)
ABSOLUTE_PATH_TOKEN_RE = re.compile(r"(?<![:.\w~$-])(/[^\s'\"`]+)")


def payload(status: str, data: dict[str, Any] | None = None,
            error: str | None = None) -> dict[str, Any]:
    output: dict[str, Any] = {"status": status}
    if data is not None:
        output["data"] = data
        output.update(data)
    if error is not None:
        output["error"] = error
    return output


def trim(text: str) -> str:
    if len(text) <= MAX_CAPTURE_CHARS:
        return text
    return text[:MAX_CAPTURE_CHARS] + "\n... [truncated]"


def unique(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        clean = value.strip()
        if clean and clean not in seen:
            seen.add(clean)
            out.append(clean)
    return out


def looks_like_pdb_id(value: str) -> bool:
    return bool(PDB_ID_RE.fullmatch(value.strip()))


def looks_like_uniprot(value: str) -> bool:
    return bool(UNIPROT_RE.fullmatch(value.strip().upper()))


def safe_label(value: str, outdir: Path | None = None) -> str:
    """Return a human-safe label for user-visible paths and identifiers."""

    text = str(value)
    if "://" in text:
        return text

    path = Path(text).expanduser()
    if path.is_absolute():
        resolved = path.resolve()
        roots = []
        if outdir is not None:
            roots.append((outdir.resolve(), "$OUTDIR"))
        roots.append((ROOT, "."))

        for root, label in roots:
            try:
                relative = resolved.relative_to(root)
            except ValueError:
                continue
            if str(relative) == ".":
                return label
            return f"{label}/{relative.as_posix()}"
        return f"{resolved.name} (absolute path omitted)"

    return text


def sanitize_string(value: str, outdir: Path) -> str:
    """Scrub private absolute paths from provenance strings."""

    if "://" in value:
        return value

    replacements: list[tuple[str, str]] = [
        (str(outdir.resolve()), "$OUTDIR"),
        (str(ROOT), "."),
        (str(Path.home()), "~"),
    ]
    text = value
    for needle, replacement in sorted(replacements, key=lambda item: len(item[0]), reverse=True):
        if needle:
            text = text.replace(needle, replacement)

    if text.startswith("/") and "\n" not in text:
        return safe_label(text, outdir)
    text = ABSOLUTE_PATH_TOKEN_RE.sub(lambda match: safe_label(match.group(1), outdir), text)
    return text


def sanitize_value(value: Any, outdir: Path) -> Any:
    if isinstance(value, dict):
        return {str(key): sanitize_value(item, outdir) for key, item in value.items()}
    if isinstance(value, list):
        return [sanitize_value(item, outdir) for item in value]
    if isinstance(value, tuple):
        return [sanitize_value(item, outdir) for item in value]
    if isinstance(value, str):
        return sanitize_string(value, outdir)
    return value


def extract_json(text: str) -> dict[str, Any] | None:
    stripped = text.strip()
    if not stripped:
        return None
    try:
        value = json.loads(stripped)
    except json.JSONDecodeError:
        value = None
    if isinstance(value, dict):
        return value

    lines = stripped.splitlines()
    for index in range(len(lines) - 1, -1, -1):
        if lines[index].lstrip().startswith("{"):
            try:
                value = json.loads("\n".join(lines[index:]))
            except json.JSONDecodeError:
                continue
            if isinstance(value, dict):
                return value
    return None


class HelperRunner:
    def __init__(self, outdir: Path, timeout: float):
        self.outdir = outdir
        self.timeout = timeout
        self.commands: list[dict[str, Any]] = []

    def run(self, script: str, args: list[str], label: str) -> dict[str, Any]:
        display_command = ["python3", f"scripts/{script}", *[str(arg) for arg in args]]
        actual_command = [PYTHON, str(SCRIPTS / script), *[str(arg) for arg in args]]
        started = time.monotonic()
        event: dict[str, Any] = {
            "label": label,
            "command": display_command,
            "script": f"scripts/{script}",
        }
        try:
            proc = subprocess.run(
                actual_command,
                cwd=ROOT,
                text=True,
                capture_output=True,
                timeout=self.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            event.update({
                "status": "error",
                "returncode": None,
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "stdout": trim(exc.stdout or ""),
                "stderr": trim(exc.stderr or ""),
                "error": f"helper timed out after {self.timeout:g} seconds",
            })
            self.commands.append(sanitize_value(event, self.outdir))
            return payload("error", error=event["error"])

        event.update({
            "returncode": proc.returncode,
            "elapsed_seconds": round(time.monotonic() - started, 3),
        })
        stdout = proc.stdout.strip()
        stderr = proc.stderr.strip()
        parsed = extract_json(stdout)
        if stderr and proc.returncode != 0:
            event["stderr"] = trim(stderr)
        if parsed is None and stdout:
            event["stdout"] = trim(stdout)

        child_status = parsed.get("status") if isinstance(parsed, dict) else None
        if proc.returncode == 0 and child_status == "ok":
            event["status"] = "ok"
            event["result_status"] = child_status
            self.commands.append(sanitize_value(event, self.outdir))
            return payload("ok", {"result": parsed})

        error = "helper did not emit ok JSON"
        if isinstance(parsed, dict) and parsed.get("error"):
            error = str(parsed["error"])
        elif proc.returncode != 0:
            error = f"helper exited with return code {proc.returncode}"
        event["status"] = "error"
        event["result_status"] = child_status
        event["error"] = error
        self.commands.append(sanitize_value(event, self.outdir))
        return payload("error", {"result": parsed}, error)


def result_data(result: dict[str, Any]) -> dict[str, Any]:
    parsed = result.get("data", {}).get("result")
    if isinstance(parsed, dict):
        data = parsed.get("data")
        if isinstance(data, dict):
            return data
        return parsed
    return {}


def format_resolution(value: Any) -> str:
    if isinstance(value, list):
        values = [str(item) for item in value if item is not None]
        return ", ".join(values) if values else "-"
    if value is None:
        return "-"
    return str(value)


def first_gene_name(match: dict[str, Any]) -> str | None:
    genes = match.get("gene_names")
    if isinstance(genes, list) and genes:
        return str(genes[0])
    gene = match.get("gene")
    if gene:
        return str(gene)
    return None


def build_dossier(args: argparse.Namespace, outdir: Path) -> dict[str, Any]:
    runner = HelperRunner(outdir, args.timeout)
    inputs: list[dict[str, Any]] = []
    targets_by_accession: dict[str, dict[str, Any]] = {}
    structures: list[dict[str, Any]] = []
    warnings: list[str] = []
    errors: list[str] = []
    enriched_uniprots: set[str] = set()

    def add_input(kind: str, value: str, status: str, summary: str,
                  details: dict[str, Any] | None = None) -> dict[str, Any]:
        record: dict[str, Any] = {
            "kind": kind,
            "value": safe_label(value, outdir),
            "status": status,
            "summary": summary,
        }
        if details:
            record["details"] = details
        inputs.append(record)
        return record

    def upsert_target(accession: str, source_kind: str, source_value: str,
                      match: dict[str, Any] | None = None) -> dict[str, Any]:
        accession = accession.upper()
        target = targets_by_accession.setdefault(accession, {
            "accession": accession,
            "status": "pending",
            "source_inputs": [],
            "notes": [],
        })
        source = {"kind": source_kind, "value": safe_label(source_value, outdir)}
        if source not in target["source_inputs"]:
            target["source_inputs"].append(source)
        if match:
            if match.get("protein_name"):
                target["protein_name"] = match["protein_name"]
            gene = first_gene_name(match)
            if gene:
                target["gene"] = gene
            if match.get("organism"):
                target["organism"] = match["organism"]
            if match.get("length") is not None:
                target["length"] = match["length"]
            if match.get("reviewed") is not None:
                target["reviewed"] = match["reviewed"]
        return target

    def enrich_uniprot(accession: str) -> None:
        accession = accession.upper()
        target = targets_by_accession[accession]
        if accession in enriched_uniprots:
            return
        enriched_uniprots.add(accession)

        if args.no_network:
            target["status"] = "skipped"
            target["notes"].append("Network access disabled; AlphaFold/PDB metadata not queried.")
            return

        resolved = runner.run(
            "resolve_structure.py",
            [
                accession,
                "--source",
                "alphafold",
                "--no-download",
                "--outdir",
                str(outdir / "structures"),
                "--json",
            ],
            f"AlphaFold metadata for {accession}",
        )
        if resolved["status"] == "ok":
            data = result_data(resolved)
            target["status"] = "ok"
            target["alphafold"] = {
                "model_id": data.get("model_id"),
                "global_plddt": data.get("global_plddt"),
                "latest_version": data.get("latest_version"),
                "source": data.get("source"),
            }
            if data.get("gene") and "gene" not in target:
                target["gene"] = data["gene"]
        else:
            target["status"] = "error"
            target["error"] = resolved.get("error")
            errors.append(f"UniProt {accession}: {resolved.get('error')}")

        search = runner.run(
            "pdb_search.py",
            ["--uniprot", accession, "--rows", str(args.rows), "--details", "--json"],
            f"RCSB structures for {accession}",
        )
        if search["status"] == "ok":
            data = result_data(search)
            results = data.get("results") if isinstance(data.get("results"), list) else []
            target["related_pdb_count"] = data.get("count", len(results))
            for item in results:
                pdb_id = str(item.get("id", "")).upper()
                if not pdb_id:
                    continue
                structures.append({
                    "kind": "pdb_search_hit",
                    "id": pdb_id,
                    "status": "ok",
                    "source": "RCSB PDB search",
                    "linked_uniprot": accession,
                    "title": item.get("title") or "",
                    "method": item.get("method") or "",
                    "resolution": item.get("resolution"),
                    "ligands": item.get("ligands") or [],
                })
        else:
            warning = f"PDB search for {accession} failed: {search.get('error')}"
            warnings.append(warning)
            target["notes"].append(warning)

    def local_analysis_summary(name: str, result: dict[str, Any]) -> dict[str, Any]:
        if result["status"] != "ok":
            return {
                "status": "error",
                "summary": str(result.get("error") or "analysis failed"),
            }
        data = result_data(result)
        if name == "ligands":
            return {
                "status": "ok",
                "summary": f"{data.get('ligand_group_count', 0)} ligand group(s)",
                "ligand_group_count": data.get("ligand_group_count", 0),
                "ligand_atom_count": data.get("ligand_atom_count", 0),
                "ligands": [
                    item.get("ligand")
                    for item in data.get("ligand_components", [])
                    if isinstance(item, dict) and item.get("ligand")
                ],
            }
        if name == "interfaces":
            return {
                "status": "ok",
                "summary": f"{data.get('interface_count', 0)} interface(s)",
                "interface_count": data.get("interface_count", 0),
                "chains": data.get("chains") or [],
            }
        if name == "interactions":
            return {
                "status": "ok",
                "summary": f"{data.get('contact_count', 0)} ligand contact(s)",
                "contact_count": data.get("contact_count", 0),
                "ligand_group_count": data.get("ligand_group_count", 0),
                "classification_counts": data.get("classification_counts") or {},
            }
        if name == "docking_box":
            box = data.get("box") if isinstance(data.get("box"), dict) else {}
            return {
                "status": "ok",
                "summary": f"{data.get('ligand_atom_count', 0)} ligand atom(s) boxed",
                "ligand_atom_count": data.get("ligand_atom_count", 0),
                "center": box.get("center"),
                "size": box.get("size"),
            }
        return {"status": "ok", "summary": "analysis completed"}

    def analyze_local_structure(path: str, fmt: str | None) -> dict[str, Any]:
        analyses: dict[str, Any] = {}
        ligand_count = None
        for name, script, extra_args in LOCAL_ANALYSIS_HELPERS:
            if name in {"interactions", "docking_box"} and ligand_count == 0:
                analyses[name] = {
                    "status": "skipped",
                    "summary": "No ligand groups were found by ligand extraction.",
                }
                continue
            result = runner.run(
                script,
                [path, *extra_args, "--json"],
                f"{name} analysis for {safe_label(path, outdir)}",
            )
            summary = local_analysis_summary(name, result)
            analyses[name] = summary
            if name == "ligands" and summary["status"] == "ok":
                ligand_count = summary.get("ligand_group_count", 0)

        if fmt not in {"pdb", None}:
            analyses["pocket"] = {
                "status": "skipped",
                "summary": "Pocket report currently supports local PDB files; skipped for this format.",
            }
            return analyses

        if ligand_count == 0:
            analyses["pocket"] = {
                "status": "skipped",
                "summary": "No ligand groups were found by ligand extraction.",
            }
            return analyses

        pocket = runner.run(
            "pocket_report.py",
            [path, "--json"],
            f"pocket analysis for {safe_label(path, outdir)}",
        )
        analyses["pocket"] = local_analysis_summary("pocket", pocket)
        if pocket["status"] == "ok":
            data = result_data(pocket)
            analyses["pocket"].update({
                "summary": f"{data.get('ligand_count', 0)} ligand pocket(s)",
                "ligand_count": data.get("ligand_count", 0),
            })
        return analyses

    for gene in unique(args.gene):
        if args.no_network:
            add_input(
                "gene",
                gene,
                "skipped",
                "Gene lookup requires UniProt network access; skipped in --no-network mode.",
            )
            continue

        lookup = runner.run(
            "uniprot_lookup.py",
            [gene, "--gene-exact", "--size", str(args.uniprot_candidates), "--json"],
            f"UniProt lookup for gene {gene}",
        )
        if lookup["status"] != "ok":
            add_input("gene", gene, "error", str(lookup.get("error") or "lookup failed"))
            errors.append(f"Gene {gene}: {lookup.get('error')}")
            continue

        data = result_data(lookup)
        hits = data.get("results") if isinstance(data.get("results"), list) else []
        if not hits:
            add_input("gene", gene, "no_match", "No reviewed UniProt match found.")
            continue

        best = hits[0]
        accession = str(best.get("accession") or "").upper()
        add_input(
            "gene",
            gene,
            "ok",
            f"Resolved to UniProt {accession}.",
            {"candidate_count": len(hits)},
        )
        if accession:
            upsert_target(accession, "gene", gene, best)
            enrich_uniprot(accession)

    for accession in unique(args.uniprot):
        normalized = accession.upper()
        if not looks_like_uniprot(normalized):
            add_input("uniprot", accession, "error", "Invalid UniProt accession format.")
            errors.append(f"UniProt {accession}: invalid accession format")
            continue

        target = upsert_target(normalized, "uniprot", accession)
        if args.no_network:
            add_input(
                "uniprot",
                accession,
                "skipped",
                "Accepted accession; metadata queries skipped in --no-network mode.",
            )
            target["status"] = "skipped"
            target["notes"].append("Network access disabled; metadata not queried.")
        else:
            add_input("uniprot", accession, "ok", "Accepted accession for metadata lookup.")
        enrich_uniprot(normalized)

    for pdb_input in unique(args.pdb):
        path = Path(pdb_input)
        if path.exists():
            inspected = runner.run(
                "structure_info.py",
                [pdb_input, "--json"],
                f"local structure inspection for {safe_label(pdb_input, outdir)}",
            )
            if inspected["status"] == "ok":
                data = result_data(inspected)
                add_input("pdb", pdb_input, "ok", "Inspected local coordinate file.")
                structure_record = {
                    "kind": "local_structure",
                    "id": safe_label(pdb_input, outdir),
                    "status": "ok",
                    "source": "local",
                    "format": data.get("format"),
                    "title": data.get("title"),
                    "chains": data.get("chains") or [],
                    "atom_records": data.get("atom_records"),
                    "hetatm_records": data.get("hetatm_records"),
                    "bfactor": data.get("bfactor") or {},
                    "likely_alphafold": data.get("likely_alphafold"),
                    "inspection": data,
                }
                if args.analyze_local:
                    structure_record["local_analyses"] = analyze_local_structure(
                        pdb_input,
                        data.get("format"),
                    )
                structures.append(structure_record)
            else:
                add_input("pdb", pdb_input, "error", str(inspected.get("error") or "inspection failed"))
                errors.append(f"PDB input {safe_label(pdb_input, outdir)}: {inspected.get('error')}")
            continue

        normalized = pdb_input.upper()
        if not looks_like_pdb_id(normalized):
            add_input("pdb", pdb_input, "error", "Not a local file and not a four-character PDB ID.")
            errors.append(f"PDB input {safe_label(pdb_input, outdir)}: invalid or missing input")
            continue

        if args.no_network:
            add_input(
                "pdb",
                pdb_input,
                "skipped",
                "Accepted PDB ID; RCSB metadata skipped in --no-network mode.",
            )
            structures.append({
                "kind": "pdb_id",
                "id": normalized,
                "status": "skipped",
                "source": "RCSB PDB",
                "notes": ["Network access disabled; metadata not queried."],
            })
            continue

        fetched = runner.run(
            "fetch_pdb.py",
            [normalized, "--metadata", "--json"],
            f"RCSB metadata for {normalized}",
        )
        if fetched["status"] == "ok":
            data = result_data(fetched)
            metadata = data.get("metadata") if isinstance(data.get("metadata"), dict) else {}
            add_input("pdb", pdb_input, "ok", "Fetched RCSB metadata.")
            structures.append({
                "kind": "pdb_id",
                "id": normalized,
                "status": "ok",
                "source": "RCSB PDB",
                "title": metadata.get("title"),
                "method": ", ".join(metadata.get("experimental_methods") or []),
                "resolution": metadata.get("resolution"),
                "ligands": metadata.get("nonpolymer_bound_components") or [],
                "metadata": metadata,
            })
        else:
            add_input("pdb", pdb_input, "error", str(fetched.get("error") or "metadata fetch failed"))
            errors.append(f"PDB {normalized}: {fetched.get('error')}")

    status = "error" if errors else "ok"
    data = {
        "schema_version": 1,
        "generated_by": "scripts/target_dossier.py",
        "network": {
            "enabled": not args.no_network,
            "mode": "offline" if args.no_network else "online",
        },
        "inputs": inputs,
        "targets": list(targets_by_accession.values()),
        "structures": structures,
        "commands": runner.commands,
        "warnings": warnings,
        "errors": errors,
        "outputs": {
            "markdown": REPORT_NAME,
            "provenance": PROVENANCE_NAME,
        },
    }
    return payload(status, data, "; ".join(errors) if errors else None)


def markdown_table(rows: list[list[str]]) -> list[str]:
    if not rows:
        return []
    widths = [max(len(row[index]) for row in rows) for index in range(len(rows[0]))]
    lines = []
    for index, row in enumerate(rows):
        padded = [cell.ljust(widths[pos]) for pos, cell in enumerate(row)]
        lines.append("| " + " | ".join(padded) + " |")
        if index == 0:
            lines.append("| " + " | ".join("-" * width for width in widths) + " |")
    return lines


def write_markdown(report: dict[str, Any], outdir: Path) -> None:
    data = report["data"]
    lines: list[str] = [
        "# Target Dossier",
        "",
        "*Generated by `scripts/target_dossier.py`.*",
        "",
        f"Status: **{report['status']}**",
        f"Network mode: **{data['network']['mode']}**",
        "",
    ]

    if not data["network"]["enabled"]:
        lines.extend([
            "> Network access disabled; remote UniProt, AlphaFold, and RCSB metadata steps were skipped.",
            "",
        ])

    lines.extend(["## Inputs", ""])
    input_rows = [["Kind", "Value", "Status", "Summary"]]
    for item in data["inputs"]:
        input_rows.append([
            item["kind"],
            item["value"],
            item["status"],
            item["summary"],
        ])
    lines.extend(markdown_table(input_rows))
    lines.append("")

    lines.extend(["## UniProt Targets", ""])
    targets = data.get("targets") or []
    if targets:
        target_rows = [["Accession", "Gene", "Protein", "Organism", "Status"]]
        for target in targets:
            target_rows.append([
                target.get("accession", "-"),
                str(target.get("gene") or "-"),
                str(target.get("protein_name") or "-"),
                str(target.get("organism") or "-"),
                str(target.get("status") or "-"),
            ])
        lines.extend(markdown_table(target_rows))
    else:
        lines.append("No UniProt target records were resolved.")
    lines.append("")

    lines.extend(["## Structure Evidence", ""])
    structures = data.get("structures") or []
    if structures:
        structure_rows = [["Kind", "ID", "Status", "Title/summary", "Resolution", "Ligands"]]
        for structure in structures:
            ligands = structure.get("ligands")
            if isinstance(ligands, list):
                ligand_text = ", ".join(str(item) for item in ligands[:8]) or "-"
            else:
                ligand_text = "-"
            title = structure.get("title") or structure.get("source") or "-"
            if structure.get("kind") == "local_structure":
                title = (
                    f"{structure.get('format', '-')}; "
                    f"{len(structure.get('chains') or [])} chain(s); "
                    f"{structure.get('atom_records', 0)} atom records"
                )
            structure_rows.append([
                str(structure.get("kind") or "-"),
                str(structure.get("id") or "-"),
                str(structure.get("status") or "-"),
                str(title),
                format_resolution(structure.get("resolution")),
                ligand_text,
            ])
        lines.extend(markdown_table(structure_rows))
    else:
        lines.append("No structure evidence records were produced.")
    lines.append("")

    local_analysis_rows = [["Structure", "Analysis", "Status", "Summary"]]
    for structure in structures:
        analyses = structure.get("local_analyses")
        if not isinstance(analyses, dict):
            continue
        for name, analysis in analyses.items():
            if not isinstance(analysis, dict):
                continue
            local_analysis_rows.append([
                str(structure.get("id") or "-"),
                name,
                str(analysis.get("status") or "-"),
                str(analysis.get("summary") or "-"),
            ])
    if len(local_analysis_rows) > 1:
        lines.extend(["## Local Analyses", ""])
        lines.extend(markdown_table(local_analysis_rows))
        lines.append("")

    if data.get("warnings"):
        lines.extend(["## Warnings", ""])
        for warning in data["warnings"]:
            lines.append(f"- {warning}")
        lines.append("")

    if data.get("errors"):
        lines.extend(["## Errors", ""])
        for error in data["errors"]:
            lines.append(f"- {error}")
        lines.append("")

    lines.extend([
        "## Provenance",
        "",
        f"Machine-readable provenance: `{PROVENANCE_NAME}`.",
        "",
    ])
    (outdir / REPORT_NAME).write_text("\n".join(lines), encoding="utf-8")


def write_outputs(report: dict[str, Any], outdir: Path) -> dict[str, Any]:
    sanitized = sanitize_value(report, outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_markdown(sanitized, outdir)
    (outdir / PROVENANCE_NAME).write_text(
        json.dumps(sanitized, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return sanitized


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a generic target dossier from gene, UniProt, and PDB inputs.",
        epilog=(
            "Examples:\n"
            "  %(prog)s --gene KRAS --uniprot P01116 --pdb 6OIM\n"
            "  %(prog)s --pdb tests/fixtures/tiny.pdb --no-network\n"
            "  %(prog)s --pdb model.cif --no-network --analyze-local\n"
            "  %(prog)s --gene TP53 --out tp53_dossier --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--gene", action="append", default=[], help="Gene symbol or protein name; repeatable")
    parser.add_argument("--uniprot", action="append", default=[], help="UniProt accession; repeatable")
    parser.add_argument("--pdb", action="append", default=[], help="PDB ID or local coordinate file; repeatable")
    parser.add_argument("--out", "--outdir", dest="outdir", default="target_dossier",
                        help="Output directory (default: target_dossier)")
    parser.add_argument("--no-network", action="store_true",
                        help="Skip all remote UniProt, AlphaFold, and RCSB lookups")
    parser.add_argument("--analyze-local", action="store_true",
                        help="Run local ligand/interface/interaction/docking-box analyses for local coordinate files")
    parser.add_argument("--rows", type=int, default=10, help="RCSB search rows per UniProt target (default: 10)")
    parser.add_argument("--uniprot-candidates", type=int, default=5,
                        help="Candidate count for gene/name lookup (default: 5)")
    parser.add_argument("--timeout", type=float, default=120,
                        help="Per-helper timeout in seconds (default: 120)")
    parser.add_argument("--json", action="store_true", dest="as_json",
                        help="Emit machine-readable summary JSON to stdout")
    args = parser.parse_args(argv)

    if not (args.gene or args.uniprot or args.pdb):
        parser.error("provide at least one --gene, --uniprot, or --pdb input")
    if args.rows < 0:
        parser.error("--rows must be >= 0")
    if args.uniprot_candidates < 1:
        parser.error("--uniprot-candidates must be >= 1")
    if args.timeout <= 0:
        parser.error("--timeout must be > 0")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    outdir = Path(args.outdir)
    report = build_dossier(args, outdir)
    written = write_outputs(report, outdir)

    if args.as_json:
        print(json.dumps(written, indent=2, sort_keys=True))
    else:
        print(f"Wrote {outdir / REPORT_NAME}")
        print(f"Wrote {outdir / PROVENANCE_NAME}")
        if written["status"] != "ok":
            print(f"Status: {written['status']}", file=sys.stderr)
    return 0 if written["status"] == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
