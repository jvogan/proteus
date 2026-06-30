#!/usr/bin/env python3
"""Plan local protein design/modeling runs from a small JSON manifest.

This helper validates a design manifest, detects optional local design tools,
and writes a run directory containing plan.json, commands.md, and
candidates.jsonl. It intentionally does not execute heavyweight modeling tools.
"""

import argparse
import importlib.util
import json
import re
import shlex
import shutil
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
ABS_PATH_TOKEN_RE = re.compile(r"(?<![\w:>])/(?:[^\s,;:'\"{}\[\]]+)")
SCALAR_TYPES = (str, int, float, bool, type(None))

TOOL_DEFINITIONS = {
    "proteinmpnn": {
        "display_name": "ProteinMPNN",
        "aliases": ["proteinmpnn", "protein_mpnn", "ProteinMPNN"],
        "executables": ["protein_mpnn_run.py", "proteinmpnn", "ProteinMPNN"],
        "imports": ["protein_mpnn", "proteinmpnn"],
    },
    "ligandmpnn": {
        "display_name": "LigandMPNN",
        "aliases": ["ligandmpnn", "ligand_mpnn", "LigandMPNN"],
        "executables": ["ligandmpnn", "ligand_mpnn", "ligand_mpnn_run.py", "LigandMPNN"],
        "imports": ["ligandmpnn", "ligand_mpnn"],
    },
    "colabfold": {
        "display_name": "ColabFold/localcolabfold",
        "aliases": ["colabfold", "localcolabfold", "colabfold_batch"],
        "executables": ["colabfold_batch", "localcolabfold", "colabfold"],
        "imports": ["colabfold", "localcolabfold"],
    },
    "boltz": {
        "display_name": "boltz",
        "aliases": ["boltz"],
        "executables": ["boltz"],
        "imports": ["boltz"],
    },
    "chai-lab": {
        "display_name": "chai-lab",
        "aliases": ["chai-lab", "chai_lab", "chai"],
        "executables": ["chai-lab", "chai_lab", "chai"],
        "imports": ["chai_lab"],
    },
}

TOOL_ALIASES = {
    alias.lower(): name
    for name, definition in TOOL_DEFINITIONS.items()
    for alias in definition["aliases"]
}
COUNT_PARAM_KEYS = (
    "candidate_count",
    "num_candidates",
    "num_sequences",
    "num_seq_per_target",
    "samples",
    "designs",
)


class ManifestError(ValueError):
    """Raised when a design manifest is syntactically valid JSON but invalid."""


def _payload(status: str, data: dict[str, Any] | None = None, error: str | None = None) -> dict[str, Any]:
    output: dict[str, Any] = {"status": status}
    if data is not None:
        output["data"] = data
        output.update(data)
    if error is not None:
        output["error"] = error
    return output


def _relative_to_root(path: Path) -> str | None:
    try:
        return str(path.resolve().relative_to(ROOT))
    except (OSError, ValueError):
        return None


def display_path(value: str | Path) -> str:
    """Return a stable, non-absolute path string for manifest/run outputs."""
    raw = str(value)
    if not raw:
        return raw
    if "://" in raw:
        return raw

    path = Path(raw).expanduser()
    if path.is_absolute():
        rel = _relative_to_root(path)
        if rel is not None:
            return "." if rel == "." else f"./{rel}"
        return f"<abs>/{path.name}"
    if raw.startswith("./") or raw.startswith("../") or raw in {".", ".."}:
        return raw
    return f"./{raw}"


def _scrub_text(text: str) -> str:
    scrubbed = text
    for path, replacement in ((ROOT, "."), (Path.home(), "~")):
        path_text = str(path)
        if path_text:
            scrubbed = scrubbed.replace(path_text, replacement)

    def replace_token(match: re.Match[str]) -> str:
        token = match.group(0)
        if "://" in token:
            return token
        return display_path(token)

    return ABS_PATH_TOKEN_RE.sub(replace_token, scrubbed)


def scrub_value(value: Any) -> Any:
    if isinstance(value, str):
        if "://" not in value and Path(value).expanduser().is_absolute():
            return display_path(value)
        return _scrub_text(value)
    if isinstance(value, list):
        return [scrub_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): scrub_value(item) for key, item in value.items()}
    return value


def _read_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ManifestError(f"Manifest is not valid JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise ManifestError("Manifest must be a JSON object.")
    return data


def _require_id(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ManifestError(f"{label} must be a non-empty string.")
    if not ID_RE.fullmatch(value):
        raise ManifestError(f"{label} may contain only letters, digits, dot, underscore, and dash.")
    return value


def _canonical_tool(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ManifestError(f"{label} must be a non-empty string.")
    canonical = TOOL_ALIASES.get(value.lower())
    if canonical is None:
        allowed = ", ".join(sorted(TOOL_DEFINITIONS))
        raise ManifestError(f"{label} must be one of: {allowed}")
    return canonical


def _ensure_string_list(value: Any, label: str, allow_empty: bool = True) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if not isinstance(value, list):
        raise ManifestError(f"{label} must be a string or list of strings.")
    if not allow_empty and not value:
        raise ManifestError(f"{label} must not be empty.")
    for index, item in enumerate(value):
        if not isinstance(item, str) or not item:
            raise ManifestError(f"{label}[{index}] must be a non-empty string.")
    return list(value)


def _validate_jsonish(value: Any, label: str) -> None:
    if isinstance(value, SCALAR_TYPES):
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _validate_jsonish(item, f"{label}[{index}]")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ManifestError(f"{label} keys must be strings.")
            _validate_jsonish(item, f"{label}.{key}")
        return
    raise ManifestError(f"{label} must contain only JSON scalar, list, or object values.")


def _infer_input_kind(path: str | None, spec: dict[str, Any]) -> str:
    if isinstance(spec.get("kind"), str) and spec["kind"]:
        return spec["kind"]
    if spec.get("sequence"):
        return "sequence"
    if not path:
        return "value"
    suffix = Path(path).suffix.lower()
    if suffix in {".pdb", ".cif", ".mmcif", ".bcif"}:
        return "structure"
    if suffix in {".fa", ".faa", ".fasta"}:
        return "sequence_file"
    if suffix in {".sdf", ".mol", ".mol2", ".pdbqt"}:
        return "ligand"
    return "file"


def _normalize_input_spec(input_id: str, raw_spec: Any, label: str) -> dict[str, Any]:
    _require_id(input_id, f"{label}.id")
    if isinstance(raw_spec, str):
        spec: dict[str, Any] = {"id": input_id, "path": raw_spec}
    elif isinstance(raw_spec, dict):
        spec = {"id": input_id, **raw_spec}
    else:
        raise ManifestError(f"{label} must be a path string or object.")

    if "id" in spec:
        spec["id"] = _require_id(spec["id"], f"{label}.id")
        if spec["id"] != input_id:
            raise ManifestError(f"{label}.id must match its inputs key.")

    path = spec.get("path")
    sequence = spec.get("sequence")
    value = spec.get("value")
    if path is not None and (not isinstance(path, str) or not path):
        raise ManifestError(f"{label}.path must be a non-empty string.")
    if sequence is not None and (not isinstance(sequence, str) or not sequence):
        raise ManifestError(f"{label}.sequence must be a non-empty string.")
    if path is None and sequence is None and value is None:
        raise ManifestError(f"{label} must contain at least one of path, sequence, or value.")
    if "description" in spec and not isinstance(spec["description"], str):
        raise ManifestError(f"{label}.description must be a string.")
    if "kind" in spec and not isinstance(spec["kind"], str):
        raise ManifestError(f"{label}.kind must be a string.")
    _validate_jsonish(spec, label)

    normalized: dict[str, Any] = {"id": input_id, "kind": _infer_input_kind(path, spec)}
    if path is not None:
        actual_path = Path(path).expanduser()
        normalized["path"] = display_path(path)
        normalized["exists"] = actual_path.exists()
    if sequence is not None:
        normalized["sequence_length"] = len(sequence)
    if value is not None:
        normalized["value"] = scrub_value(value)
    for key in ("description", "chains", "residues", "metadata"):
        if key in spec:
            normalized[key] = scrub_value(spec[key])
    return normalized


def _normalize_inputs(raw_inputs: Any) -> list[dict[str, Any]]:
    if not isinstance(raw_inputs, (dict, list)) or not raw_inputs:
        raise ManifestError("Manifest must contain a non-empty inputs object or list.")

    inputs: list[dict[str, Any]] = []
    seen: set[str] = set()
    if isinstance(raw_inputs, dict):
        iterable = raw_inputs.items()
        for input_id, raw_spec in iterable:
            if not isinstance(input_id, str):
                raise ManifestError("inputs object keys must be strings.")
            normalized = _normalize_input_spec(input_id, raw_spec, f"inputs.{input_id}")
            if normalized["id"] in seen:
                raise ManifestError(f"Duplicate input id: {normalized['id']}")
            seen.add(normalized["id"])
            inputs.append(normalized)
        return inputs

    for index, raw_spec in enumerate(raw_inputs):
        if not isinstance(raw_spec, dict):
            raise ManifestError(f"inputs[{index}] must be an object.")
        input_id = raw_spec.get("id") or raw_spec.get("name")
        if not isinstance(input_id, str):
            raise ManifestError(f"inputs[{index}].id must be a non-empty string.")
        normalized = _normalize_input_spec(input_id, raw_spec, f"inputs[{index}]")
        if normalized["id"] in seen:
            raise ManifestError(f"Duplicate input id: {normalized['id']}")
        seen.add(normalized["id"])
        inputs.append(normalized)
    return inputs


def _normalize_constraints(raw_constraints: Any) -> dict[str, Any]:
    if raw_constraints is None:
        return {}
    if not isinstance(raw_constraints, dict):
        raise ManifestError("constraints must be an object.")
    _validate_jsonish(raw_constraints, "constraints")
    constraints = scrub_value(raw_constraints)

    for key in ("max_candidates", "candidate_count"):
        if key in constraints:
            value = constraints[key]
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ManifestError(f"constraints.{key} must be a positive integer.")

    if "required_tools" in constraints:
        required_tools = _ensure_string_list(constraints["required_tools"], "constraints.required_tools")
        constraints["required_tools"] = [
            _canonical_tool(tool, "constraints.required_tools[]") for tool in required_tools
        ]
    return constraints


def _normalize_stages(raw_stages: Any, input_ids: set[str]) -> list[dict[str, Any]]:
    if not isinstance(raw_stages, list) or not raw_stages:
        raise ManifestError("Manifest must contain a non-empty stages list.")
    stages: list[dict[str, Any]] = []
    seen: set[str] = set()
    known_refs = set(input_ids)
    for index, raw_stage in enumerate(raw_stages):
        label = f"stages[{index}]"
        if not isinstance(raw_stage, dict):
            raise ManifestError(f"{label} must be an object.")
        stage_id = _require_id(raw_stage.get("id") or raw_stage.get("name"), f"{label}.id")
        if stage_id in seen:
            raise ManifestError(f"Duplicate stage id: {stage_id}")
        seen.add(stage_id)

        tool = _canonical_tool(raw_stage.get("tool"), f"{label}.tool")
        action = raw_stage.get("action", "run")
        if not isinstance(action, str) or not action:
            raise ManifestError(f"{label}.action must be a non-empty string.")

        stage_inputs = _ensure_string_list(raw_stage.get("inputs"), f"{label}.inputs")
        depends_on = _ensure_string_list(raw_stage.get("depends_on"), f"{label}.depends_on")
        for ref in [*stage_inputs, *depends_on]:
            if ref not in known_refs:
                raise ManifestError(f"{label} references unknown input or prior stage: {ref}")

        params = raw_stage.get("params", {})
        if params is None:
            params = {}
        if not isinstance(params, dict):
            raise ManifestError(f"{label}.params must be an object.")
        _validate_jsonish(params, f"{label}.params")

        stage_constraints = raw_stage.get("constraints", {})
        if stage_constraints is None:
            stage_constraints = {}
        if not isinstance(stage_constraints, dict):
            raise ManifestError(f"{label}.constraints must be an object.")
        _validate_jsonish(stage_constraints, f"{label}.constraints")

        stage = {
            "id": stage_id,
            "tool": tool,
            "tool_name": TOOL_DEFINITIONS[tool]["display_name"],
            "action": action,
            "inputs": stage_inputs,
            "depends_on": depends_on,
            "params": scrub_value(params),
            "constraints": scrub_value(stage_constraints),
        }
        stages.append(stage)
        known_refs.add(stage_id)
    return stages


def normalize_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    inputs = _normalize_inputs(manifest.get("inputs"))
    constraints = _normalize_constraints(manifest.get("constraints"))
    stages = _normalize_stages(manifest.get("stages"), {item["id"] for item in inputs})

    name = manifest.get("name", "design-run")
    if not isinstance(name, str) or not name:
        raise ManifestError("name must be a non-empty string when provided.")
    description = manifest.get("description")
    if description is not None and not isinstance(description, str):
        raise ManifestError("description must be a string when provided.")

    normalized = {
        "name": name,
        "description": description,
        "inputs": inputs,
        "constraints": constraints,
        "stages": stages,
    }
    return {key: value for key, value in normalized.items() if value is not None}


def load_manifest(path: Path) -> dict[str, Any]:
    return normalize_manifest(_read_json(path))


def _find_executable(candidates: list[str]) -> tuple[str | None, str | None]:
    for executable in candidates:
        found = shutil.which(executable)
        if found:
            return executable, found
    return None, None


def _find_import(candidates: list[str]) -> str | None:
    for import_name in candidates:
        try:
            if importlib.util.find_spec(import_name) is not None:
                return import_name
        except (ImportError, AttributeError, ValueError):
            continue
    return None


def detect_tools() -> dict[str, Any]:
    tools: dict[str, Any] = {}
    for name, definition in TOOL_DEFINITIONS.items():
        executable, executable_path = _find_executable(definition["executables"])
        import_name = _find_import(definition["imports"])
        detection = None
        if executable_path:
            detection = "executable"
        elif import_name:
            detection = "import"

        tools[name] = {
            "available": bool(executable_path or import_name),
            "detection": detection,
            "display_name": definition["display_name"],
            "executable": executable,
            "path": display_path(executable_path) if executable_path else None,
            "import_name": import_name,
            "executable_available": bool(executable_path),
            "import_available": bool(import_name),
            "executable_candidates": list(definition["executables"]),
            "import_candidates": list(definition["imports"]),
        }
    return tools


def _first_input_path(stage: dict[str, Any], inputs_by_id: dict[str, dict[str, Any]]) -> str:
    refs = stage["inputs"] or list(inputs_by_id)
    for ref in refs:
        item = inputs_by_id.get(ref)
        if item and item.get("path"):
            return item["path"]
    return "<input>"


def _first_ligand_path(inputs_by_id: dict[str, dict[str, Any]]) -> str | None:
    for item in inputs_by_id.values():
        if item.get("kind") == "ligand" and item.get("path"):
            return item["path"]
    return None


def _param_flags(params: dict[str, Any]) -> list[str]:
    flags: list[str] = []
    for key in sorted(params):
        value = params[key]
        flag = "--" + key.replace("_", "-")
        if value is True:
            flags.append(flag)
        elif value is False or value is None:
            continue
        elif isinstance(value, list):
            for item in value:
                flags.extend([flag, str(scrub_value(item))])
        elif isinstance(value, dict):
            flags.extend([flag, json.dumps(scrub_value(value), sort_keys=True, separators=(",", ":"))])
        else:
            flags.extend([flag, str(scrub_value(value))])
    return flags


def build_stage_command(stage: dict[str, Any], inputs_by_id: dict[str, dict[str, Any]], tools: dict[str, Any]) -> list[str]:
    tool = stage["tool"]
    tool_info = tools[tool]
    executable = tool_info["executable"] or TOOL_DEFINITIONS[tool]["executables"][0]
    stage_dir = f"<run>/stages/{stage['id']}"
    primary_input = _first_input_path(stage, inputs_by_id)

    if tool in {"proteinmpnn", "ligandmpnn"}:
        command = [executable, "--pdb_path", primary_input, "--out_folder", stage_dir]
        if tool == "ligandmpnn":
            ligand_path = _first_ligand_path(inputs_by_id)
            if ligand_path:
                command.extend(["--ligand_path", ligand_path])
    elif tool == "colabfold":
        sequence_input = primary_input if primary_input != "<input>" else "<run>/candidates.fasta"
        command = [executable, sequence_input, stage_dir]
    elif tool == "boltz":
        command = [executable, "predict", primary_input, "--out_dir", stage_dir]
    elif tool == "chai-lab":
        command = [executable, "fold", "--input", primary_input, "--output-dir", stage_dir]
    else:
        command = [executable, primary_input, "--out-dir", stage_dir]

    return [str(part) for part in [*command, *_param_flags(stage["params"])]]


def _candidate_count(manifest: dict[str, Any]) -> int:
    constraints = manifest.get("constraints", {})
    if "candidate_count" in constraints:
        count = constraints["candidate_count"]
    else:
        count = None
        for stage in manifest["stages"]:
            for key in COUNT_PARAM_KEYS:
                value = stage["params"].get(key)
                if isinstance(value, int) and not isinstance(value, bool) and value > 0:
                    count = value
                    break
            if count is not None:
                break
        if count is None:
            count = constraints.get("max_candidates", 1)

    max_candidates = constraints.get("max_candidates")
    if isinstance(max_candidates, int) and not isinstance(max_candidates, bool):
        count = min(count, max_candidates)
    return max(1, int(count))


def build_candidates(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    count = _candidate_count(manifest)
    pending_stages = [stage["id"] for stage in manifest["stages"]]
    return [
        {
            "candidate_id": f"candidate_{index:04d}",
            "status": "planned",
            "rank": None,
            "source": "manifest",
            "pending_stages": pending_stages,
            "constraint_status": "unchecked",
        }
        for index in range(1, count + 1)
    ]


def build_plan(manifest_path: Path, outdir: Path, dry_run: bool) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    tools = detect_tools()
    inputs_by_id = {item["id"]: item for item in manifest["inputs"]}

    stages = []
    for stage in manifest["stages"]:
        command = build_stage_command(stage, inputs_by_id, tools)
        available = tools[stage["tool"]]["available"]
        stages.append(
            {
                **stage,
                "tool_available": available,
                "command": command,
                "command_line": shlex.join(command),
                "mode": "manual",
            }
        )

    candidates = build_candidates(manifest)
    required_tools = manifest.get("constraints", {}).get("required_tools", [])
    missing_required = [tool for tool in required_tools if not tools[tool]["available"]]
    missing_stage_tools = sorted({stage["tool"] for stage in stages if not tools[stage["tool"]]["available"]})

    warnings = []
    for item in manifest["inputs"]:
        if item.get("path") and not item.get("exists"):
            warnings.append(f"Input path is not present locally: {item['id']} ({item['path']})")
    for tool in missing_stage_tools:
        warnings.append(f"Stage tool is not currently detected: {tool}")
    for tool in missing_required:
        warnings.append(f"Required tool is not currently detected: {tool}")

    readiness = "blocked" if missing_required else "ready"
    if missing_stage_tools and readiness == "ready":
        readiness = "partial"

    return {
        "status": "ok",
        "readiness": readiness,
        "mode": "plan_only",
        "dry_run": dry_run,
        "heavy_tools_executed": False,
        "manifest": display_path(manifest_path),
        "run_dir": "<run>",
        "name": manifest["name"],
        "description": manifest.get("description"),
        "inputs": manifest["inputs"],
        "constraints": manifest.get("constraints", {}),
        "tools": tools,
        "stages": stages,
        "candidates": {
            "count": len(candidates),
            "file": "<run>/candidates.jsonl",
        },
        "outputs": {
            "plan_json": "<run>/plan.json",
            "commands_markdown": "<run>/commands.md",
            "candidates_jsonl": "<run>/candidates.jsonl",
        },
        "counts": {
            "inputs": len(manifest["inputs"]),
            "stages": len(stages),
            "candidates": len(candidates),
            "tools_available": sum(1 for tool in tools.values() if tool["available"]),
            "tools_missing_for_stages": len(missing_stage_tools),
            "required_tools_missing": len(missing_required),
        },
        "warnings": warnings,
    }


def _commands_markdown(plan: dict[str, Any]) -> str:
    lines = [
        f"# Design Run: {plan['name']}",
        "",
        "This run directory is a plan. Review commands before running them manually.",
        "",
        f"Readiness: {plan['readiness']}",
        f"Dry run requested: {'yes' if plan['dry_run'] else 'no'}",
        "Heavy tools executed: no",
        "",
        "## Stages",
        "",
    ]
    for stage in plan["stages"]:
        available = "available" if stage["tool_available"] else "not detected"
        lines.extend(
            [
                f"### {stage['id']} ({stage['tool_name']}, {available})",
                "",
                "```bash",
                stage["command_line"],
                "```",
                "",
            ]
        )
    if plan["warnings"]:
        lines.extend(["## Warnings", ""])
        for warning in plan["warnings"]:
            lines.append(f"- {warning}")
        lines.append("")
    return "\n".join(lines)


def write_run(manifest_path: Path, outdir: Path, dry_run: bool) -> dict[str, Any]:
    plan = scrub_value(build_plan(manifest_path, outdir, dry_run))
    candidates = build_candidates(load_manifest(manifest_path))

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "plan.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (outdir / "commands.md").write_text(_commands_markdown(plan), encoding="utf-8")
    with (outdir / "candidates.jsonl").open("w", encoding="utf-8") as handle:
        for candidate in candidates:
            handle.write(json.dumps(scrub_value(candidate), sort_keys=True) + "\n")

    data = {
        "readiness": plan["readiness"],
        "dry_run": dry_run,
        "heavy_tools_executed": False,
        "manifest": plan["manifest"],
        "run_dir": "<run>",
        "outputs": plan["outputs"],
        "counts": plan["counts"],
        "warnings": plan["warnings"],
    }
    return _payload("ok", data)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate and plan a local protein design/modeling run without executing heavy tools.",
        epilog=(
            "Example manifest:\n"
            "{\n"
            '  "name": "kras-pocket-redesign",\n'
            '  "inputs": {\n'
            '    "backbone": {"path": "structures/target.pdb", "kind": "structure"},\n'
            '    "ligand": {"path": "ligands/gdp.sdf", "kind": "ligand"}\n'
            "  },\n"
            '  "constraints": {\n'
            '    "fixed_residues": ["A:10", "A:11"],\n'
            '    "max_candidates": 8,\n'
            '    "required_tools": ["proteinmpnn"]\n'
            "  },\n"
            '  "stages": [\n'
            '    {"id": "design", "tool": "proteinmpnn", "inputs": ["backbone"],\n'
            '     "params": {"num_seq_per_target": 8}},\n'
            '    {"id": "fold", "tool": "colabfold", "depends_on": ["design"],\n'
            '     "params": {"num_models": 1}}\n'
            "  ]\n"
            "}\n\n"
            "Supported tools: "
            + ", ".join(definition["display_name"] for definition in TOOL_DEFINITIONS.values())
            + "\n"
            "Outputs: plan.json, commands.md, candidates.jsonl."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("manifest", help="Path to a JSON design manifest")
    parser.add_argument("--outdir", default="design_run_out", help="Directory for planned run artifacts")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Record an explicit dry run; the helper always plans and never runs heavy tools",
    )
    parser.add_argument("--json", action="store_true", help="Print a machine-readable JSON summary")
    args = parser.parse_args(argv)

    try:
        summary = write_run(Path(args.manifest), Path(args.outdir), args.dry_run)
    except (ManifestError, OSError) as exc:
        output = _payload("error", error=_scrub_text(str(exc)))
        if args.json:
            print(json.dumps(output, indent=2, sort_keys=True))
        else:
            print(f"ERROR: {output['error']}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print(f"Readiness: {summary['readiness']}")
        print("Heavy tools executed: no")
        print("Plan: <run>/plan.json")
        print("Commands: <run>/commands.md")
        print("Candidates: <run>/candidates.jsonl")
    return 0


if __name__ == "__main__":
    sys.exit(main())
