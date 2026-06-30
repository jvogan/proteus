#!/usr/bin/env python3
"""Combine Proteus helper JSON reports into a portable evidence pack.

The combiner is stdlib-only and intentionally conservative: it writes a
human-readable Markdown summary plus a scrubbed machine-readable report.json
with per-input checksums and provenance. Original JSON inputs are copied into
an evidence/ directory by default; use --no-copy to write only the combined
outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT_JSON_NAME = "report.json"
REPORT_MARKDOWN_NAME = "REPORT.md"
EVIDENCE_DIR_NAME = "evidence"

MAX_EMBEDDED_STRING_CHARS = 12000
MAX_COORDINATE_STRING_CHARS = 512
MAX_COORDINATE_COLLECTION_ITEMS = 100
MAX_HUGE_COLLECTION_ITEMS = 10000

ABSOLUTE_PATH_TOKEN_RE = re.compile(r"(?<![:.\w~$-])(/[^\s'\"`]+)")
ID_CLEAN_RE = re.compile(r"[^A-Za-z0-9_.-]+")
PDB_LINE_RE = re.compile(r"^(ATOM  |HETATM|MODEL |ENDMDL|TER   )")

RAW_TEXT_KEYS = {
    "cif",
    "cif_text",
    "content",
    "contents",
    "coordinate_content",
    "coordinate_contents",
    "coordinate_text",
    "coordinates_text",
    "file_contents",
    "mmcif",
    "mmcif_text",
    "pdb",
    "pdb_contents",
    "pdb_text",
    "raw",
    "raw_text",
    "structure_contents",
    "structure_text",
}
COORDINATE_COLLECTION_KEYS = {
    "atom_site",
    "atoms",
    "coordinates",
    "frames",
    "models",
}
COMMON_SUMMARY_KEYS = (
    "source",
    "query",
    "pdb_id",
    "model_id",
    "format",
    "title",
    "resolution",
    "assembly_count",
    "ligand_count",
    "ligand_group_count",
    "interface_count",
    "chain_count",
    "chains",
    "atom_records",
    "hetatm_records",
    "warnings",
)


class CombineError(RuntimeError):
    pass


def payload(status: str, data: dict[str, Any] | None = None,
            error: str | None = None) -> dict[str, Any]:
    output: dict[str, Any] = {"status": status}
    if data is not None:
        output["data"] = data
        output.update(data)
    if error:
        output["error"] = error
    return output


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8", errors="replace"))


def stable_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def safe_label(value: str, outdir: Path | None = None) -> str:
    """Return a human-safe label for paths without exposing private prefixes."""

    text = str(value)
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9+.-]*://\S+", text):
        return text

    path = Path(text).expanduser()
    if path.is_absolute():
        resolved = path.resolve()
        roots: list[tuple[Path, str]] = []
        if outdir is not None:
            roots.append((outdir.resolve(), "$OUTDIR"))
        roots.extend([(ROOT, "."), (Path.home().resolve(), "~")])

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


def _sanitize_path_match(match: re.Match[str], outdir: Path) -> str:
    token = match.group(1)
    suffix = ""
    while token and token[-1] in ".,;:)]}":
        suffix = token[-1] + suffix
        token = token[:-1]
    if not token:
        return match.group(1)
    return safe_label(token, outdir) + suffix


def sanitize_string(value: str, outdir: Path) -> str:
    """Scrub absolute local paths from strings preserved in report outputs."""

    if re.fullmatch(r"[A-Za-z][A-Za-z0-9+.-]*://\S+", value):
        return value

    replacements: list[tuple[str, str]] = [
        (str(outdir.resolve()), "$OUTDIR"),
        (str(ROOT), "."),
        (str(Path.home().resolve()), "~"),
    ]
    text = value
    for needle, replacement in sorted(replacements, key=lambda item: len(item[0]), reverse=True):
        if needle:
            text = text.replace(needle, replacement)

    if text.startswith("/") and "\n" not in text:
        return safe_label(text, outdir)
    return ABSOLUTE_PATH_TOKEN_RE.sub(lambda match: _sanitize_path_match(match, outdir), text)


def _normalized_key(path: tuple[str, ...]) -> str:
    if not path:
        return ""
    return str(path[-1]).strip().lower().replace("-", "_")


def _dot_path(path: tuple[str, ...]) -> str:
    return ".".join(str(part) for part in path) or "$"


def looks_like_coordinate_text(text: str) -> bool:
    lines = text.splitlines()
    pdb_lines = sum(1 for line in lines[:80] if PDB_LINE_RE.match(line))
    if pdb_lines >= 3:
        return True
    return text.count("_atom_site.") >= 3


def looks_like_coordinate_item(value: Any) -> bool:
    if isinstance(value, dict):
        keys = {str(key).lower() for key in value}
        if {"x", "y", "z"}.issubset(keys):
            return True
        if {"cartn_x", "cartn_y", "cartn_z"}.issubset(keys):
            return True
        if {"atom_name", "residue_name", "chain", "x", "y", "z"}.issubset(keys):
            return True
    if isinstance(value, list) and len(value) >= 3:
        return all(isinstance(item, (int, float)) for item in value[:3])
    return False


def omitted_value(kind: str, reason: str, value: Any) -> dict[str, Any]:
    if isinstance(value, str):
        return {
            "omitted": True,
            "kind": kind,
            "reason": reason,
            "chars": len(value),
            "sha256": sha256_text(value),
        }
    encoded = stable_json_bytes(value)
    output = {
        "omitted": True,
        "kind": kind,
        "reason": reason,
        "bytes": len(encoded),
        "sha256": sha256_bytes(encoded),
    }
    if isinstance(value, list):
        output["items"] = len(value)
    if isinstance(value, dict):
        output["keys"] = len(value)
    return output


def should_omit_string(key: str, value: str) -> tuple[bool, str]:
    if len(value) > MAX_EMBEDDED_STRING_CHARS:
        return True, "large string value omitted from combined report"
    if len(value) > MAX_COORDINATE_STRING_CHARS and key in RAW_TEXT_KEYS:
        return True, "raw coordinate or file contents omitted from combined report"
    if len(value) > MAX_COORDINATE_STRING_CHARS and looks_like_coordinate_text(value):
        return True, "coordinate-like text omitted from combined report"
    return False, ""


def should_omit_collection(key: str, value: list[Any]) -> tuple[bool, str]:
    if len(value) > MAX_HUGE_COLLECTION_ITEMS:
        return True, "huge collection omitted from combined report"
    if key in COORDINATE_COLLECTION_KEYS and len(value) > MAX_COORDINATE_COLLECTION_ITEMS:
        return True, "large coordinate-like collection omitted from combined report"
    if len(value) > MAX_COORDINATE_COLLECTION_ITEMS:
        sample = value[:20]
        if sample and sum(1 for item in sample if looks_like_coordinate_item(item)) >= min(5, len(sample)):
            return True, "coordinate-like collection omitted from combined report"
    return False, ""


def scrub_value(value: Any, outdir: Path, warnings: list[str], source_label: str,
                path: tuple[str, ...] = ()) -> Any:
    key = _normalized_key(path)
    if isinstance(value, dict):
        return {
            str(item_key): scrub_value(item_value, outdir, warnings, source_label, path + (str(item_key),))
            for item_key, item_value in value.items()
        }
    if isinstance(value, list):
        omit, reason = should_omit_collection(key, value)
        if omit:
            warnings.append(f"{source_label}: omitted `{_dot_path(path)}` ({reason}).")
            return omitted_value("list", reason, value)
        return [scrub_value(item, outdir, warnings, source_label, path) for item in value]
    if isinstance(value, str):
        omit, reason = should_omit_string(key, value)
        if omit:
            warnings.append(f"{source_label}: omitted `{_dot_path(path)}` ({reason}).")
            return omitted_value("string", reason, value)
        return sanitize_string(value, outdir)
    return value


def report_id_for(path: Path, index: int, seen: set[str]) -> str:
    base = ID_CLEAN_RE.sub("-", path.stem).strip(".-") or f"input-{index}"
    candidate = base
    suffix = 2
    while candidate in seen:
        candidate = f"{base}-{suffix}"
        suffix += 1
    seen.add(candidate)
    return candidate


def copy_input(raw: bytes, source: Path, outdir: Path, report_id: str, index: int) -> str:
    evidence_dir = outdir / EVIDENCE_DIR_NAME
    evidence_dir.mkdir(parents=True, exist_ok=True)
    suffix = source.suffix if source.suffix else ".json"
    destination = evidence_dir / f"{index:02d}-{report_id}{suffix}"
    destination.write_bytes(raw)
    return f"{EVIDENCE_DIR_NAME}/{destination.name}"


def unwrap_report(raw_json: Any) -> tuple[str, str, Any, dict[str, Any]]:
    if isinstance(raw_json, dict):
        raw_status = raw_json.get("status")
        status = raw_status if isinstance(raw_status, str) and raw_status else "ok"
        if "data" in raw_json and isinstance(raw_status, str):
            data = raw_json.get("data")
            metadata = {
                key: value
                for key, value in raw_json.items()
                if key != "data" and not (isinstance(data, dict) and key in data)
            }
            return "proteus_envelope", status, data, metadata
        if isinstance(raw_status, str):
            metadata = {key: value for key, value in raw_json.items() if key in {"status", "error"}}
            data = {key: value for key, value in raw_json.items() if key not in {"status", "error"}}
            return "status_json", status, data, metadata
    return "plain_json", "ok", raw_json, {}


def summarize_report(status: str, data: Any, metadata: dict[str, Any]) -> dict[str, Any]:
    error = metadata.get("error")
    if isinstance(error, str) and error:
        return {"text": error, "fields": []}

    if isinstance(data, dict):
        parts = []
        for key in COMMON_SUMMARY_KEYS:
            if key not in data:
                continue
            value = data[key]
            if isinstance(value, (str, int, float, bool)) or value is None:
                parts.append(f"{key}={value}")
            elif isinstance(value, list):
                if key == "warnings":
                    parts.append(f"warnings={len(value)}")
                elif len(value) <= 8:
                    parts.append(f"{key}={', '.join(str(item) for item in value)}")
                else:
                    parts.append(f"{key}={len(value)} items")
            elif isinstance(value, dict):
                parts.append(f"{key}={len(value)} fields")
        if not parts:
            parts.append(f"{len(data)} top-level fields")
        return {
            "text": "; ".join(parts),
            "fields": sorted(str(key) for key in data.keys())[:40],
        }

    if isinstance(data, list):
        return {"text": f"{len(data)} list items", "fields": []}
    return {"text": f"{status} {type(data).__name__} report", "fields": []}


def load_one_report(path: Path, outdir: Path, index: int, report_id: str,
                    copy_inputs: bool, warnings: list[str], errors: list[str]) -> dict[str, Any]:
    source_label = safe_label(str(path), outdir)
    entry: dict[str, Any] = {
        "id": report_id,
        "source": source_label,
        "status": "error",
        "shape": "unreadable",
        "provenance": {},
        "summary": {"text": "input could not be read", "fields": []},
        "metadata": {},
        "data": {},
    }

    try:
        raw = path.read_bytes()
    except OSError as exc:
        message = f"{source_label}: could not read input JSON ({sanitize_string(str(exc), outdir)})"
        errors.append(message)
        entry["metadata"] = {"error": message}
        return entry

    entry["provenance"] = {
        "bytes": len(raw),
        "sha256": sha256_bytes(raw),
    }

    try:
        raw_json = json.loads(raw.decode("utf-8"))
    except UnicodeDecodeError as exc:
        message = f"{source_label}: input is not UTF-8 ({exc})"
        errors.append(message)
        entry["shape"] = "invalid_json"
        entry["summary"] = {"text": message, "fields": []}
        entry["metadata"] = {"error": message}
        if copy_inputs:
            entry["provenance"]["copied_to"] = copy_input(raw, path, outdir, report_id, index)
        return entry
    except json.JSONDecodeError as exc:
        message = f"{source_label}: input is not valid JSON ({exc})"
        errors.append(message)
        entry["shape"] = "invalid_json"
        entry["summary"] = {"text": message, "fields": []}
        entry["metadata"] = {"error": message}
        if copy_inputs:
            entry["provenance"]["copied_to"] = copy_input(raw, path, outdir, report_id, index)
        return entry

    shape, status, data, metadata = unwrap_report(raw_json)
    entry["shape"] = shape
    entry["status"] = status
    if shape == "plain_json" and not isinstance(raw_json, dict):
        warnings.append(f"{source_label}: plain JSON input is {type(raw_json).__name__}, not an object.")
    if status != "ok":
        message = f"{source_label}: source report status is {status!r}."
        if status == "error":
            errors.append(message)
        else:
            warnings.append(message)

    scrub_warnings_start = len(warnings)
    scrubbed_metadata = scrub_value(metadata, outdir, warnings, source_label, ("metadata",))
    scrubbed_data = scrub_value(data, outdir, warnings, source_label, ("data",))
    entry["metadata"] = scrubbed_metadata
    entry["data"] = scrubbed_data
    entry["summary"] = scrub_value(
        summarize_report(status, data, metadata),
        outdir,
        warnings,
        source_label,
        ("summary",),
    )
    if len(warnings) > scrub_warnings_start:
        entry["provenance"]["pruned"] = True
    if copy_inputs:
        if entry["provenance"].get("pruned"):
            entry["provenance"]["copy_skipped"] = "pruned_large_content"
            warnings.append(
                f"{source_label}: input copy skipped because large raw content was pruned; "
                "source checksum is retained in provenance."
            )
        else:
            entry["provenance"]["copied_to"] = copy_input(raw, path, outdir, report_id, index)
    return entry


def markdown_escape(value: Any) -> str:
    text = str(value).replace("\n", " ").replace("|", "\\|")
    return text if len(text) <= 180 else text[:177] + "..."


def markdown_table(rows: list[list[Any]]) -> list[str]:
    if not rows:
        return []
    widths = [max(len(markdown_escape(row[index])) for row in rows) for index in range(len(rows[0]))]
    lines = []
    for index, row in enumerate(rows):
        cells = [markdown_escape(row[pos]).ljust(widths[pos]) for pos in range(len(row))]
        lines.append("| " + " | ".join(cells) + " |")
        if index == 0:
            lines.append("| " + " | ".join("-" * width for width in widths) + " |")
    return lines


def render_markdown(report: dict[str, Any]) -> str:
    data = report["data"]
    counts = data["counts"]
    lines: list[str] = [
        f"# {data['title']}",
        "",
        "*Generated by `scripts/proteus_report.py`.*",
        "",
        f"Status: **{report['status']}**",
        f"Inputs: **{counts['inputs_total']}** total, **{counts['inputs_ok']}** ok, "
        f"**{counts['inputs_error']}** error, **{counts['inputs_warning']}** warning.",
        "",
        "## Reports",
        "",
    ]
    rows: list[list[Any]] = [["ID", "Status", "Shape", "Source", "Summary"]]
    for item in data["reports"]:
        rows.append([
            item["id"],
            item["status"],
            item["shape"],
            item["source"],
            item["summary"]["text"],
        ])
    lines.extend(markdown_table(rows))
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
        f"Machine-readable report: `{REPORT_JSON_NAME}`.",
    ])
    if data["copy_inputs"]:
        copied = counts.get("inputs_copied", 0)
        skipped = counts.get("inputs_copy_skipped", 0)
        if copied:
            lines.append(f"Copied evidence inputs: `{EVIDENCE_DIR_NAME}/` ({copied} file(s)).")
        else:
            lines.append("No input files were copied.")
        if skipped:
            lines.append(f"Skipped input copies: {skipped}.")
    else:
        lines.append("Input copying disabled by `--no-copy`.")
    lines.append("")
    return "\n".join(lines)


def combine_reports(inputs: list[Path], outdir: Path, title: str,
                    copy_inputs: bool = True) -> dict[str, Any]:
    if not inputs:
        raise CombineError("provide at least one --input JSON report")

    outdir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    errors: list[str] = []
    reports: list[dict[str, Any]] = []
    seen_ids: set[str] = set()

    for index, path in enumerate(inputs, start=1):
        report_id = report_id_for(path, index, seen_ids)
        reports.append(load_one_report(path, outdir, index, report_id, copy_inputs, warnings, errors))

    counts = {
        "inputs_total": len(reports),
        "inputs_ok": sum(1 for item in reports if item["status"] == "ok"),
        "inputs_error": sum(1 for item in reports if item["status"] == "error"),
        "inputs_warning": sum(1 for item in reports if item["status"] not in {"ok", "error"}),
        "inputs_copied": sum(1 for item in reports if "copied_to" in item.get("provenance", {})),
        "inputs_copy_skipped": sum(1 for item in reports if "copy_skipped" in item.get("provenance", {})),
        "warnings": len(warnings),
        "errors": len(errors),
    }
    data = {
        "title": sanitize_string(title, outdir),
        "created_utc": utc_now(),
        "generator": "scripts/proteus_report.py",
        "copy_inputs": copy_inputs,
        "outputs": {
            "markdown": REPORT_MARKDOWN_NAME,
            "json": REPORT_JSON_NAME,
            "evidence_dir": EVIDENCE_DIR_NAME if copy_inputs else None,
        },
        "counts": counts,
        "warnings": warnings,
        "errors": errors,
        "reports": reports,
    }
    status = "error" if errors else "ok"
    combined = payload(status, data, "; ".join(errors) if errors else None)
    (outdir / REPORT_JSON_NAME).write_text(
        json.dumps(combined, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (outdir / REPORT_MARKDOWN_NAME).write_text(render_markdown(combined), encoding="utf-8")
    return combined


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine Proteus helper JSON reports into Markdown and report.json evidence outputs.",
        epilog=(
            "Examples:\n"
            "  %(prog)s --input pocket.json --input validation.json --outdir evidence_pack\n"
            "  %(prog)s --input report.json --outdir evidence_pack --no-copy --json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", action="append", default=[], dest="inputs",
                        help="Proteus helper JSON report file; repeat for multiple inputs")
    parser.add_argument("--outdir", default="proteus_report", help="Output directory")
    parser.add_argument("--title", default="Proteus Evidence Report", help="Markdown/report title")
    parser.add_argument("--json", action="store_true", dest="as_json",
                        help="Emit the combined report JSON to stdout")
    parser.add_argument("--no-copy", action="store_true",
                        help="Do not copy input JSON reports into the evidence/ output directory")
    args = parser.parse_args(argv)
    if not args.inputs:
        parser.error("provide at least one --input JSON report")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    outdir = Path(args.outdir)
    try:
        combined = combine_reports(
            [Path(item) for item in args.inputs],
            outdir,
            args.title,
            copy_inputs=not args.no_copy,
        )
    except (CombineError, OSError) as exc:
        error = sanitize_string(str(exc), outdir)
        output = payload("error", error=error)
        if args.as_json:
            print(json.dumps(output, indent=2, sort_keys=True))
        else:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2

    if args.as_json:
        print(json.dumps(combined, indent=2, sort_keys=True))
    else:
        print(f"Wrote {safe_label(str(outdir / REPORT_MARKDOWN_NAME), outdir)}")
        print(f"Wrote {safe_label(str(outdir / REPORT_JSON_NAME), outdir)}")
        if combined["status"] != "ok":
            print(f"Status: {combined['status']}", file=sys.stderr)
    return 1 if combined["status"] == "error" else 0


if __name__ == "__main__":
    sys.exit(main())
