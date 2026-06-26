#!/usr/bin/env python3
"""Run Proteus helper scripts from a small JSON manifest."""

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TASKS = {
    "assembly_report": "scripts/assembly_report.py",
    "docking_box": "scripts/docking_box.py",
    "structure_info": "scripts/structure_info.py",
    "ligand_extract": "scripts/ligand_extract.py",
    "interaction_report": "scripts/interaction_report.py",
    "pocket_report": "scripts/pocket_report.py",
    "interface_report": "scripts/interface_report.py",
    "mutation_triage": "scripts/mutation_triage.py",
    "validation_report": "scripts/validation_report.py",
    "resolve_structure": "scripts/resolve_structure.py",
}
ITEM_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
MAX_CAPTURE_CHARS = 20000


class ManifestError(ValueError):
    pass


def _payload(status: str, data: dict | None = None, error: str | None = None) -> dict:
    output = {"status": status}
    if data is not None:
        output["data"] = data
        output.update(data)
    if error is not None:
        output["error"] = error
    return output


def _trim(value: str) -> str:
    if len(value) <= MAX_CAPTURE_CHARS:
        return value
    return value[:MAX_CAPTURE_CHARS] + "\n... [truncated]"


def load_manifest(path: Path) -> dict:
    try:
        manifest = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ManifestError(f"Manifest is not valid JSON: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ManifestError("Manifest must be a JSON object.")
    items = manifest.get("items")
    if not isinstance(items, list) or not items:
        raise ManifestError("Manifest must contain a non-empty items list.")

    seen_ids = set()
    for index, item in enumerate(items):
        if not isinstance(item, dict):
            raise ManifestError(f"items[{index}] must be an object.")
        item_id = item.get("id")
        if not isinstance(item_id, str) or not item_id:
            raise ManifestError(f"items[{index}].id must be a non-empty string.")
        if not ITEM_ID_RE.fullmatch(item_id):
            raise ManifestError(
                f"items[{index}].id may contain only letters, digits, dot, underscore, and dash."
            )
        if item_id in seen_ids:
            raise ManifestError(f"Duplicate item id: {item_id}")
        seen_ids.add(item_id)

        item_input = item.get("input")
        if not isinstance(item_input, str) or not item_input:
            raise ManifestError(f"items[{index}].input must be a non-empty string.")

        tasks = item.get("tasks")
        if not isinstance(tasks, list) or not tasks:
            raise ManifestError(f"items[{index}].tasks must be a non-empty list.")
        for task_index, task in enumerate(tasks):
            name, extra_args = normalize_task_spec(task, f"items[{index}].tasks[{task_index}]")
            if name not in TASKS:
                allowed = ", ".join(sorted(TASKS))
                raise ManifestError(
                    f"items[{index}].tasks[{task_index}] must be one of: {allowed}"
                )
            for arg_index, arg in enumerate(extra_args):
                if not isinstance(arg, str):
                    raise ManifestError(
                        f"items[{index}].tasks[{task_index}].args[{arg_index}] must be a string."
                    )
    return manifest


def normalize_task_spec(task, label: str = "task") -> tuple[str, list[str]]:
    if isinstance(task, str):
        return task, []
    if not isinstance(task, dict):
        raise ManifestError(f"{label} must be a string or object.")
    name = task.get("name")
    if not isinstance(name, str) or not name:
        raise ManifestError(f"{label}.name must be a non-empty string.")
    args = task.get("args", [])
    if args is None:
        args = []
    if not isinstance(args, list):
        raise ManifestError(f"{label}.args must be a list of strings.")
    return name, args


def _has_option(args: list[str], option: str) -> bool:
    prefix = option + "="
    return any(arg == option or arg.startswith(prefix) for arg in args)


def build_command(task_spec, item_input: str, item_outdir: Path) -> list[str]:
    task, extra_args = normalize_task_spec(task_spec)
    script = TASKS[task]
    if task == "mutation_triage":
        if not extra_args:
            raise ManifestError("mutation_triage task requires at least one variant in args.")
        command = ["python3", script, *extra_args]
        if not _has_option(extra_args, "--structure"):
            command.extend(["--structure", item_input])
    else:
        command = ["python3", script, item_input]
    if task in {
        "docking_box",
        "ligand_extract",
        "pocket_report",
        "interface_report",
        "interaction_report",
        "resolve_structure",
        "assembly_report",
    } and not _has_option(extra_args, "--outdir"):
        command.extend(["--outdir", str(item_outdir)])
    if task != "mutation_triage":
        command.extend(extra_args)
    if not _has_option(extra_args, "--json"):
        command.append("--json")
    return command


def _actual_command(display_command: list[str]) -> list[str]:
    return [sys.executable, *display_command[1:]]


def run_task(task: str, item: dict, item_outdir: Path, dry_run: bool, timeout: float) -> dict:
    task_name, _ = normalize_task_spec(task)
    display_command = build_command(task, item["input"], item_outdir)
    base_data = {
        "item_id": item["id"],
        "input": item["input"],
        "task": task_name,
        "task_spec": task,
        "command": display_command,
        "dry_run": dry_run,
    }
    if dry_run:
        data = {
            **base_data,
            "task_status": "planned",
            "returncode": None,
            "elapsed_seconds": 0.0,
        }
        return _payload("ok", data)

    item_outdir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    try:
        proc = subprocess.run(
            _actual_command(display_command),
            cwd=ROOT,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        elapsed = round(time.monotonic() - started, 3)
    except subprocess.TimeoutExpired as exc:
        data = {
            **base_data,
            "task_status": "error",
            "returncode": None,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "stdout": _trim(exc.stdout or ""),
            "stderr": _trim(exc.stderr or ""),
        }
        return _payload("error", data, f"Task timed out after {timeout:g} seconds.")

    stdout = proc.stdout.strip()
    stderr = proc.stderr.strip()
    data = {
        **base_data,
        "returncode": proc.returncode,
        "elapsed_seconds": elapsed,
    }
    if stderr:
        data["stderr"] = _trim(stderr)

    parsed = None
    parse_error = None
    if stdout:
        try:
            parsed = json.loads(stdout)
        except json.JSONDecodeError as exc:
            parse_error = f"stdout was not valid JSON: {exc}"
            data["stdout"] = _trim(stdout)
    else:
        parse_error = "stdout was empty."

    if parsed is not None:
        data["result"] = parsed

    child_status = parsed.get("status") if isinstance(parsed, dict) else None
    if proc.returncode == 0 and child_status == "ok":
        data["task_status"] = "ok"
        return _payload("ok", data)

    data["task_status"] = "error"
    error = None
    if isinstance(parsed, dict):
        error = parsed.get("error")
    if error is None:
        if child_status and child_status != "ok":
            error = f"Child JSON status was {child_status!r}."
        elif proc.returncode != 0:
            error = f"Task exited with return code {proc.returncode}."
        else:
            error = parse_error or "Child JSON status was not ok."
    return _payload("error", data, error)


def _summarize_item(item: dict, task_results: list[dict]) -> dict:
    planned = sum(1 for result in task_results if result["data"]["task_status"] == "planned")
    errors = sum(1 for result in task_results if result["status"] == "error")
    ok = sum(1 for result in task_results if result["data"]["task_status"] == "ok")
    if errors:
        status = "error"
    elif planned == len(task_results):
        status = "planned"
    else:
        status = "ok"
    return {
        "id": item["id"],
        "status": status,
        "task_count": len(task_results),
        "tasks_ok": ok,
        "tasks_error": errors,
        "tasks_planned": planned,
        "tasks": [
            {"name": result["data"]["task"], "status": result["data"]["task_status"]}
            for result in task_results
        ],
    }


def _summary_markdown(summary: dict) -> str:
    data = summary["data"]
    counts = data["counts"]
    lines = [
        "# Proteus Batch Summary",
        "",
        f"Status: {summary['status']}",
        f"Dry run: {'yes' if data['dry_run'] else 'no'}",
        "",
        "| Item | Status | Tasks | OK | Error | Planned |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for item in data["items"]:
        lines.append(
            f"| {item['id']} | {item['status']} | {item['task_count']} | "
            f"{item['tasks_ok']} | {item['tasks_error']} | {item['tasks_planned']} |"
        )
    lines.extend([
        "",
        f"Items: {counts['items_total']} total, {counts['items_ok']} ok, "
        f"{counts['items_error']} error, {counts['items_planned']} planned.",
        f"Tasks: {counts['tasks_total']} total, {counts['tasks_ok']} ok, "
        f"{counts['tasks_error']} error, {counts['tasks_planned']} planned.",
        "",
        "Machine-readable outputs: results.jsonl and summary.json.",
        "",
    ])
    return "\n".join(lines)


def run_batch(manifest_path: Path, outdir: Path, dry_run: bool, timeout: float) -> dict:
    manifest = load_manifest(manifest_path)
    outdir.mkdir(parents=True, exist_ok=True)
    results_path = outdir / "results.jsonl"
    summary_json_path = outdir / "summary.json"
    summary_md_path = outdir / "summary.md"

    item_summaries = []
    all_results = []
    with results_path.open("w") as handle:
        for item in manifest["items"]:
            item_outdir = outdir / "items" / item["id"]
            task_results = []
            for task in item["tasks"]:
                result = run_task(task, item, item_outdir, dry_run, timeout)
                task_results.append(result)
                all_results.append(result)
                handle.write(json.dumps(result, sort_keys=True) + "\n")
                handle.flush()
            item_summaries.append(_summarize_item(item, task_results))

    counts = {
        "items_total": len(item_summaries),
        "items_ok": sum(1 for item in item_summaries if item["status"] == "ok"),
        "items_error": sum(1 for item in item_summaries if item["status"] == "error"),
        "items_planned": sum(1 for item in item_summaries if item["status"] == "planned"),
        "tasks_total": len(all_results),
        "tasks_ok": sum(1 for result in all_results if result["data"]["task_status"] == "ok"),
        "tasks_error": sum(1 for result in all_results if result["status"] == "error"),
        "tasks_planned": sum(1 for result in all_results if result["data"]["task_status"] == "planned"),
    }
    data = {
        "manifest": str(manifest_path),
        "dry_run": dry_run,
        "outputs": {
            "results_jsonl": str(results_path),
            "summary_json": str(summary_json_path),
            "summary_markdown": str(summary_md_path),
        },
        "counts": counts,
        "items": item_summaries,
    }
    status = "error" if counts["tasks_error"] else "ok"
    summary = _payload(status, data)
    summary_json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    summary_md_path.write_text(_summary_markdown(summary))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run a JSON manifest batch across selected Proteus helper scripts.",
        epilog=(
            "Manifest shape:\n"
            '  {"items":[{"id":"tiny","input":"tests/fixtures/tiny.pdb",'
            '"tasks":["structure_info"]}]}\n\n'
            "Supported tasks: " + ", ".join(sorted(TASKS))
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("manifest", help="Path to a JSON manifest")
    parser.add_argument("--outdir", default="proteus_batch_out", help="Output directory")
    parser.add_argument("--dry-run", action="store_true", help="Write planned commands without running tasks")
    parser.add_argument("--timeout", type=float, default=300.0, help="Per-task timeout in seconds")
    parser.add_argument("--json", action="store_true", help="Print the summary JSON envelope")
    args = parser.parse_args()

    try:
        summary = run_batch(Path(args.manifest), Path(args.outdir), args.dry_run, args.timeout)
    except (ManifestError, OSError) as exc:
        output = _payload("error", error=str(exc))
        if args.json:
            print(json.dumps(output, indent=2, sort_keys=True))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        counts = summary["data"]["counts"]
        print(f"Status: {summary['status']}")
        print(f"Items: {counts['items_ok']} ok, {counts['items_error']} error, {counts['items_planned']} planned")
        print(f"Tasks: {counts['tasks_ok']} ok, {counts['tasks_error']} error, {counts['tasks_planned']} planned")
        print(f"Results: {summary['data']['outputs']['results_jsonl']}")
        print(f"Summary: {summary['data']['outputs']['summary_json']}")
    return 1 if summary["status"] == "error" else 0


if __name__ == "__main__":
    sys.exit(main())
