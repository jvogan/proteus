#!/usr/bin/env python3
"""Shared runtime helpers for Proteus command-line workflows.

Keep this module dependency-free. It centralizes result envelopes, privacy-safe
path labels, checksums, executable discovery, subprocess handling, and resilient
HTTP access so workflow scripts do not each invent subtly different behavior.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1
PROTEUS_VERSION = "0.2.0"
USER_AGENT = f"proteus/{PROTEUS_VERSION}"
TRANSIENT_HTTP_CODES = {408, 425, 429, 500, 502, 503, 504}
SECRET_QUERY_NAMES = {
    "access_token", "api_key", "apikey", "auth", "authorization", "key",
    "password", "secret", "sig", "signature", "token",
}
ABSOLUTE_PATH_RE = re.compile(r"(?<![:A-Za-z0-9_./~$-])/(?!/)[^\s\"'`;]+")


class ProteusRuntimeError(RuntimeError):
    """Base error for shared Proteus runtime failures."""


@dataclass(frozen=True)
class HttpResult:
    body: bytes
    url: str
    status: int
    headers: dict[str, str]
    cached: bool = False
    attempts: int = 1


def display_path(value: str | Path, *, roots: list[Path] | None = None) -> str:
    """Return a useful path label without leaking an arbitrary absolute path."""
    candidate = Path(value).expanduser()
    try:
        resolved = candidate.resolve(strict=False)
    except OSError:
        resolved = candidate.absolute()

    search_roots = list(roots or [])
    search_roots.extend([Path.cwd(), Path.home()])
    seen: set[str] = set()
    for root in search_roots:
        try:
            resolved_root = root.expanduser().resolve(strict=False)
        except OSError:
            continue
        marker = str(resolved_root)
        if marker in seen:
            continue
        seen.add(marker)
        try:
            relative = resolved.relative_to(resolved_root)
        except ValueError:
            continue
        if resolved_root == Path.home().resolve(strict=False):
            return f"~/{relative}" if str(relative) != "." else "~"
        return f"./{relative}" if str(relative) != "." else "."
    return f"{resolved.name or 'path'} (absolute path omitted)"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def file_provenance(path: str | Path, *, kind: str = "local_file") -> dict[str, Any]:
    candidate = Path(path)
    data: dict[str, Any] = {"kind": kind, "path": display_path(candidate)}
    if candidate.is_file():
        data.update({"bytes": candidate.stat().st_size, "sha256": sha256_file(candidate)})
    return data


def scrub_text(value: str | bytes) -> str:
    """Remove arbitrary absolute local paths from captured tool output."""
    text = value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value
    replacements = [(str(Path.cwd().resolve()), "."), (str(Path.home().resolve()), "~")]
    for source, label in sorted(replacements, key=lambda item: len(item[0]), reverse=True):
        text = text.replace(source, label)

    def replace(match: re.Match[str]) -> str:
        token = match.group(0)
        suffix = ""
        while token and token[-1] in ".,:)]}":
            suffix = token[-1] + suffix
            token = token[:-1]
        return display_path(token) + suffix if token else match.group(0)

    return ABSOLUTE_PATH_RE.sub(replace, text)


def scrub_private(value: Any) -> Any:
    """Recursively scrub absolute paths from JSON-safe execution metadata."""
    if isinstance(value, dict):
        return {str(key): scrub_private(item) for key, item in value.items()}
    if isinstance(value, list):
        return [scrub_private(item) for item in value]
    if isinstance(value, tuple):
        return [scrub_private(item) for item in value]
    if isinstance(value, str):
        return scrub_text(value)
    return value


def ok_payload(
    data: dict[str, Any],
    *,
    warnings: list[str] | None = None,
    provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    output: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "proteus_version": PROTEUS_VERSION,
        "status": "ok",
        "data": data,
    }
    if warnings:
        output["warnings"] = list(dict.fromkeys(warnings))
    if provenance:
        output["provenance"] = provenance
    return output


def error_payload(message: str, *, code: str | None = None) -> dict[str, Any]:
    output: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "proteus_version": PROTEUS_VERSION,
        "status": "error",
        "error": message,
    }
    if code:
        output["code"] = code
    return output


def write_json(path: str | Path, value: Any) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return destination


def slug(value: str, default: str = "workflow") -> str:
    normalized = re.sub(r"[^a-z0-9]+", "-", value.strip().lower()).strip("-")
    return normalized or default


def find_executable(*candidates: str) -> str | None:
    for candidate in candidates:
        if not candidate:
            continue
        path = shutil.which(candidate)
        if path:
            return path
        expanded = Path(candidate).expanduser()
        if expanded.is_file() and os.access(expanded, os.X_OK):
            return str(expanded)
    return None


def runtime_metadata() -> dict[str, Any]:
    return {
        "python": platform.python_version(),
        "platform": platform.system(),
        "machine": platform.machine(),
    }


def run_command(
    command: list[str],
    *,
    timeout: int = 300,
    cwd: str | Path | None = None,
    env: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Run an argv command without a shell and return a JSON-safe record."""
    try:
        proc = subprocess.run(
            command,
            cwd=str(cwd) if cwd is not None else None,
            env=env,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "status": "error",
            "error": f"Command timed out after {timeout} seconds.",
            "returncode": None,
            "stdout": scrub_text(exc.stdout or ""),
            "stderr": scrub_text(exc.stderr or ""),
        }
    except OSError as exc:
        return {"status": "error", "error": scrub_text(str(exc)), "returncode": None, "stdout": "", "stderr": ""}
    return {
        "status": "ok" if proc.returncode == 0 else "error",
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }


def _safe_url(url: str) -> str:
    from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

    parsed = urlsplit(url)
    if parsed.scheme not in {"https", "http"}:
        raise ProteusRuntimeError("Only HTTP(S) URLs are supported.")
    if parsed.username or parsed.password:
        raise ProteusRuntimeError("Credential-bearing URLs are not allowed.")
    for name, _value in parse_qsl(parsed.query, keep_blank_values=True):
        if name.lower() in SECRET_QUERY_NAMES:
            raise ProteusRuntimeError(f"Secret-like URL query parameter is not allowed: {name}")
    return urlunsplit((parsed.scheme, parsed.netloc, parsed.path, urlencode(parse_qsl(parsed.query, keep_blank_values=True)), ""))


def _retry_delay(headers: Any, attempt: int, backoff: float) -> float:
    value = headers.get("Retry-After") if headers is not None else None
    if value:
        try:
            return max(0.0, min(float(value), 30.0))
        except ValueError:
            try:
                target = parsedate_to_datetime(value).timestamp()
                return max(0.0, min(target - time.time(), 30.0))
            except (TypeError, ValueError, OverflowError):
                pass
    return min(backoff * (2 ** max(0, attempt - 1)), 30.0)


def _cache_paths(cache_dir: str | Path, url: str) -> tuple[Path, Path]:
    root = Path(cache_dir).expanduser()
    key = sha256_bytes(url.encode("utf-8"))
    return root / f"{key}.body", root / f"{key}.json"


def request_bytes(
    url: str,
    *,
    method: str = "GET",
    data: bytes | None = None,
    headers: dict[str, str] | None = None,
    timeout: int = 30,
    retries: int = 3,
    backoff: float = 1.0,
    cache_dir: str | Path | None = None,
    offline: bool = False,
) -> HttpResult:
    safe_url = _safe_url(url)
    body_path: Path | None = None
    meta_path: Path | None = None
    if cache_dir is not None:
        body_path, meta_path = _cache_paths(cache_dir, safe_url)
        if body_path.is_file() and meta_path.is_file():
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
            body = body_path.read_bytes()
            if metadata.get("sha256") == sha256_bytes(body):
                return HttpResult(
                    body=body,
                    url=safe_url,
                    status=int(metadata.get("status", 200)),
                    headers=dict(metadata.get("headers") or {}),
                    cached=True,
                    attempts=0,
                )
        if offline:
            raise ProteusRuntimeError("Offline mode requested and no valid cached response exists.")
    elif offline:
        raise ProteusRuntimeError("Offline mode requires a cache directory.")

    request_headers = {"User-Agent": USER_AGENT, **(headers or {})}
    attempts = max(1, retries + 1)
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(safe_url, data=data, headers=request_headers, method=method)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                body = response.read()
                result = HttpResult(
                    body=body,
                    url=response.geturl(),
                    status=int(getattr(response, "status", 200)),
                    headers=dict(response.headers.items()),
                    attempts=attempt,
                )
                if body_path is not None and meta_path is not None:
                    body_path.parent.mkdir(parents=True, exist_ok=True)
                    body_path.write_bytes(body)
                    write_json(meta_path, {
                        "url": safe_url,
                        "status": result.status,
                        "headers": result.headers,
                        "sha256": sha256_bytes(body),
                        "fetched_at": int(time.time()),
                    })
                return result
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code not in TRANSIENT_HTTP_CODES or attempt >= attempts:
                raise ProteusRuntimeError(f"HTTP {exc.code} for {safe_url}: {exc.reason}") from exc
            time.sleep(_retry_delay(exc.headers, attempt, backoff))
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
            if attempt >= attempts:
                raise ProteusRuntimeError(f"Request failed for {safe_url}: {exc}") from exc
            time.sleep(_retry_delay(None, attempt, backoff))
    raise ProteusRuntimeError(f"Request failed for {safe_url}: {last_error}")


def request_json(url: str, **kwargs: Any) -> tuple[Any, HttpResult]:
    result = request_bytes(url, **kwargs)
    try:
        return json.loads(result.body), result
    except json.JSONDecodeError as exc:
        raise ProteusRuntimeError(f"Response from {result.url} was not valid JSON: {exc}") from exc
