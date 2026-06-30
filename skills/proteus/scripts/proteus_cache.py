#!/usr/bin/env python3
"""Stdlib-only local cache for Proteus public API/download workflows.

The cache is keyed by SHA256(URL) and stores two files per entry:

    .proteus-cache/<key>.body
    .proteus-cache/<key>.json

Metadata JSON records the source URL, fetch timestamp, HTTP status, content
type, body SHA256, and byte count.

Usage:
    python scripts/proteus_cache.py fetch URL --json
    python scripts/proteus_cache.py fetch URL --offline --out payload.bin
    python scripts/proteus_cache.py show URL
    python scripts/proteus_cache.py verify
    python scripts/proteus_cache.py gc --max-age-seconds 86400
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import shutil
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


DEFAULT_CACHE_DIR = ".proteus-cache"
DEFAULT_TIMEOUT = 60
USER_AGENT = "proteus-skill/1.0"
CHUNK_SIZE = 1024 * 1024

SECRET_EXACT_QUERY_NAMES = {
    "api_key",
    "apikey",
    "auth",
    "authorization",
    "bearer",
    "key",
    "password",
    "sig",
    "signature",
    "token",
}
SECRET_QUERY_MARKERS = (
    "access_token",
    "access_key",
    "credential",
    "secret",
    "session_token",
)


class CacheError(RuntimeError):
    """User-facing cache operation error."""


def _ok_payload(data: dict) -> dict:
    output = {"status": "ok", "data": data}
    output.update(data)
    return output


def _error_payload(message: str, data: dict | None = None) -> dict:
    output: dict = {"status": "error", "error": message}
    if data is not None:
        output["data"] = data
        output.update(data)
    return output


def cache_key(url: str) -> str:
    """Return the cache key for a URL."""
    return hashlib.sha256(url.encode("utf-8")).hexdigest()


def _entry_paths(cache_dir: Path, key: str) -> tuple[Path, Path]:
    return cache_dir / f"{key}.body", cache_dir / f"{key}.json"


def _body_hash_and_size(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    byte_count = 0
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(CHUNK_SIZE)
            if not chunk:
                break
            byte_count += len(chunk)
            digest.update(chunk)
    return digest.hexdigest(), byte_count


def _now_timestamp() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _parse_timestamp(value: object) -> dt.datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def _is_expired(metadata: dict, ttl_seconds: int | None) -> bool:
    if ttl_seconds is None:
        return False
    timestamp = _parse_timestamp(metadata.get("timestamp"))
    if timestamp is None:
        return True
    age = dt.datetime.now(dt.timezone.utc) - timestamp
    return age.total_seconds() > ttl_seconds


def _normalize_query_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def _looks_secret_query_name(name: str) -> bool:
    normalized = _normalize_query_name(name)
    if normalized in SECRET_EXACT_QUERY_NAMES:
        return True
    return any(marker in normalized for marker in SECRET_QUERY_MARKERS)


def validate_public_url(url: str) -> str:
    """Reject unsupported schemes and URL forms likely to contain secrets."""
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme not in {"http", "https"}:
        raise CacheError("URL must use http or https.")
    if not parsed.netloc:
        raise CacheError("URL must include a host.")
    if parsed.username or parsed.password:
        raise CacheError("Refusing URL with embedded credentials.")
    for name, _value in urllib.parse.parse_qsl(parsed.query, keep_blank_values=True):
        if _looks_secret_query_name(name):
            raise CacheError(f"Refusing URL with secret-like query parameter: {name}")
    return url


def _load_metadata(cache_dir: Path, key: str, *, strict: bool = True) -> dict | None:
    _body_path, metadata_path = _entry_paths(cache_dir, key)
    if not metadata_path.exists():
        return None
    try:
        with metadata_path.open("r", encoding="utf-8") as handle:
            metadata = json.load(handle)
    except json.JSONDecodeError as exc:
        if strict:
            raise CacheError(f"Invalid metadata JSON for cache key {key}: {exc}") from exc
        return None
    except OSError as exc:
        if strict:
            raise CacheError(f"Could not read metadata for cache key {key}: {exc}") from exc
        return None
    if not isinstance(metadata, dict):
        if strict:
            raise CacheError(f"Metadata for cache key {key} is not a JSON object.")
        return None
    return metadata


def _write_json(path: Path, data: dict):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _download_to_cache(url: str, cache_dir: Path) -> dict:
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = cache_key(url)
    body_path, metadata_path = _entry_paths(cache_dir, key)
    body_tmp = cache_dir / f".{key}.body.tmp"
    metadata_tmp = cache_dir / f".{key}.json.tmp"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    try:
        with urllib.request.urlopen(request, timeout=DEFAULT_TIMEOUT) as response:
            status = int(response.getcode())
            content_type = response.headers.get("Content-Type")
            digest = hashlib.sha256()
            byte_count = 0
            with body_tmp.open("wb") as handle:
                while True:
                    chunk = response.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    byte_count += len(chunk)
                    digest.update(chunk)
                    handle.write(chunk)
    except urllib.error.HTTPError as exc:
        raise CacheError(f"Fetch failed for URL with HTTP {exc.code}: {exc.reason}") from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise CacheError(f"Fetch failed for URL: {exc}") from exc

    metadata = {
        "url": url,
        "timestamp": _now_timestamp(),
        "status": status,
        "content_type": content_type,
        "body_sha256": digest.hexdigest(),
        "byte_count": byte_count,
    }

    try:
        _write_json(metadata_tmp, metadata)
        body_tmp.replace(body_path)
        metadata_tmp.replace(metadata_path)
    except OSError as exc:
        raise CacheError(f"Could not write cache entry: {exc}") from exc
    finally:
        for path in (body_tmp, metadata_tmp):
            try:
                if path.exists():
                    path.unlink()
            except OSError:
                pass

    return metadata


def _scrub_path_for_output(path: Path) -> str:
    """Avoid emitting machine-specific temp roots in JSON or text output."""
    text = str(path)
    try:
        resolved = path.expanduser().resolve()
        temp_root = Path(tempfile.gettempdir()).resolve()
        if resolved == temp_root:
            return "<tmp>"
        if temp_root in resolved.parents:
            return str(Path("<tmp>") / resolved.relative_to(temp_root))

        cwd = Path.cwd().resolve()
        if resolved == cwd:
            return "."
        if cwd in resolved.parents:
            return str(resolved.relative_to(cwd))
    except OSError:
        return text
    return text


def fetch_url(
    url: str,
    cache_dir: Path,
    *,
    offline: bool = False,
    ttl_seconds: int | None = None,
    out: Path | None = None,
) -> dict:
    validate_public_url(url)
    key = cache_key(url)
    body_path, _metadata_path = _entry_paths(cache_dir, key)
    metadata = _load_metadata(cache_dir, key, strict=False)
    cache_available = metadata is not None and body_path.exists()
    stale = bool(metadata and _is_expired(metadata, ttl_seconds))
    cache_hit = False

    if cache_available and (offline or not stale):
        cache_hit = True
    elif offline:
        raise CacheError(f"No cached entry for URL: {url}")
    else:
        metadata = _download_to_cache(url, cache_dir)
        stale = False

    if metadata is None or not body_path.exists():
        raise CacheError(f"No cached body for URL: {url}")

    out_display = None
    if out is not None:
        if out.parent != Path(""):
            out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(body_path, out)
        out_display = _scrub_path_for_output(out)

    data = {
        "url": url,
        "cache_key": key,
        "cache_hit": cache_hit,
        "stale": stale,
        "metadata": metadata,
        "byte_count": metadata.get("byte_count"),
        "body_sha256": metadata.get("body_sha256"),
    }
    if out_display is not None:
        data["out"] = out_display
    return data


def show_cache(cache_dir: Path, url: str | None = None) -> dict:
    if url is not None:
        validate_public_url(url)
        key = cache_key(url)
        body_path, _metadata_path = _entry_paths(cache_dir, key)
        metadata = _load_metadata(cache_dir, key)
        if metadata is None:
            raise CacheError(f"No cached entry for URL: {url}")
        return {
            "count": 1,
            "entries": [
                {
                    "cache_key": key,
                    "body_present": body_path.exists(),
                    "metadata": metadata,
                }
            ],
        }

    entries = []
    if cache_dir.exists():
        for metadata_path in sorted(cache_dir.glob("*.json")):
            key = metadata_path.stem
            body_path, _metadata_path = _entry_paths(cache_dir, key)
            metadata = _load_metadata(cache_dir, key)
            entries.append({"cache_key": key, "body_present": body_path.exists(), "metadata": metadata})
    return {"count": len(entries), "entries": entries}


def _cache_keys(cache_dir: Path) -> set[str]:
    if not cache_dir.exists():
        return set()
    keys = {path.stem for path in cache_dir.glob("*.json")}
    keys.update(path.stem for path in cache_dir.glob("*.body"))
    return keys


def verify_cache(cache_dir: Path, url: str | None = None) -> tuple[dict, bool]:
    if url is not None:
        validate_public_url(url)
        keys = [cache_key(url)]
    else:
        keys = sorted(_cache_keys(cache_dir))

    entries = []
    all_ok = True
    for key in keys:
        body_path, metadata_path = _entry_paths(cache_dir, key)
        errors = []
        metadata = None
        if not metadata_path.exists():
            errors.append("missing metadata")
        else:
            try:
                metadata = _load_metadata(cache_dir, key)
            except CacheError as exc:
                errors.append(str(exc))

        actual_sha = None
        actual_byte_count = None
        if not body_path.exists():
            errors.append("missing body")
        else:
            actual_sha, actual_byte_count = _body_hash_and_size(body_path)

        if metadata and actual_sha is not None:
            if metadata.get("body_sha256") != actual_sha:
                errors.append("body sha256 mismatch")
            if metadata.get("byte_count") != actual_byte_count:
                errors.append("byte count mismatch")

        ok = not errors
        all_ok = all_ok and ok
        entry = {
            "cache_key": key,
            "ok": ok,
            "errors": errors,
            "actual_body_sha256": actual_sha,
            "actual_byte_count": actual_byte_count,
        }
        if metadata:
            entry["url"] = metadata.get("url")
            entry["metadata"] = metadata
        entries.append(entry)

    data = {"checked": len(entries), "ok": all_ok, "entries": entries}
    return data, all_ok


def gc_cache(cache_dir: Path, max_age_seconds: int | None = None, *, dry_run: bool = False) -> dict:
    removed = []
    if not cache_dir.exists():
        return {"removed_entries": 0, "removed_files": 0, "dry_run": dry_run, "entries": removed}

    for key in sorted(_cache_keys(cache_dir)):
        body_path, metadata_path = _entry_paths(cache_dir, key)
        reason = None
        metadata = None
        if not metadata_path.exists():
            reason = "missing metadata"
        else:
            try:
                metadata = _load_metadata(cache_dir, key)
            except CacheError:
                reason = "invalid metadata"

        if reason is None and not body_path.exists():
            reason = "missing body"
        if reason is None and max_age_seconds is not None and metadata is not None:
            if _is_expired(metadata, max_age_seconds):
                reason = "expired"

        if reason is None:
            continue

        files = [path for path in (body_path, metadata_path) if path.exists()]
        removed.append({"cache_key": key, "reason": reason, "file_count": len(files)})
        if not dry_run:
            for path in files:
                path.unlink()

    removed_files = sum(entry["file_count"] for entry in removed)
    return {
        "removed_entries": len(removed),
        "removed_files": removed_files,
        "dry_run": dry_run,
        "entries": removed,
    }


def _nonnegative_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Cache public Proteus API/download responses locally using only Python stdlib.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s fetch https://files.rcsb.org/download/1HSG.cif --json\n"
            "  %(prog)s fetch https://example.org/data.json --offline --out data.json\n"
            "  %(prog)s show https://example.org/data.json --json\n"
            "  %(prog)s verify --json\n"
            "  %(prog)s gc --max-age-seconds 604800"
        ),
    )
    parser.add_argument(
        "--cache-dir",
        default=DEFAULT_CACHE_DIR,
        help=f"Cache directory (default: {DEFAULT_CACHE_DIR})",
    )
    command_parent = argparse.ArgumentParser(add_help=False)
    command_parent.add_argument(
        "--cache-dir",
        default=argparse.SUPPRESS,
        help=f"Cache directory (default: {DEFAULT_CACHE_DIR})",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch_parser = subparsers.add_parser("fetch", parents=[command_parent], help="Fetch URL into cache")
    fetch_parser.add_argument("url", help="Public http(s) URL to fetch")
    fetch_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    fetch_parser.add_argument("--offline", action="store_true", help="Use cache only; error if missing")
    fetch_parser.add_argument(
        "--ttl-seconds",
        type=_nonnegative_int,
        help="Refresh cached entry when older than this many seconds",
    )
    fetch_parser.add_argument("--out", type=Path, help="Copy cached body to this file")

    show_parser = subparsers.add_parser("show", parents=[command_parent], help="Show cache metadata")
    show_parser.add_argument("url", nargs="?", help="URL to show; omit to list all entries")
    show_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    verify_parser = subparsers.add_parser("verify", parents=[command_parent], help="Verify cached body hashes")
    verify_parser.add_argument("url", nargs="?", help="URL to verify; omit to verify all entries")
    verify_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")

    gc_parser = subparsers.add_parser("gc", parents=[command_parent], help="Remove stale or broken cache entries")
    gc_parser.add_argument(
        "--max-age-seconds",
        type=_nonnegative_int,
        help="Remove entries older than this many seconds",
    )
    gc_parser.add_argument("--dry-run", action="store_true", help="Report removals without deleting files")
    gc_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    return parser


def _print_fetch(data: dict):
    print("Cache hit." if data["cache_hit"] else "Fetched.")
    print(f"URL: {data['url']}")
    print(f"Cache key: {data['cache_key']}")
    print(f"Status: {data['metadata'].get('status')}")
    print(f"Content-Type: {data['metadata'].get('content_type') or '(unknown)'}")
    print(f"Bytes: {data['byte_count']:,}")
    print(f"Body SHA256: {data['body_sha256']}")
    if data.get("out"):
        print(f"Wrote: {data['out']}")


def _print_show(data: dict):
    print(f"Entries: {data['count']}")
    for entry in data["entries"]:
        metadata = entry["metadata"]
        print(f"{entry['cache_key']}  {metadata.get('status')}  {metadata.get('byte_count')} bytes")
        print(f"  URL: {metadata.get('url')}")
        print(f"  Timestamp: {metadata.get('timestamp')}")
        print(f"  Content-Type: {metadata.get('content_type') or '(unknown)'}")
        print(f"  Body present: {'yes' if entry['body_present'] else 'no'}")


def _print_verify(data: dict):
    if data["ok"]:
        print(f"Cache verification ok: {data['checked']} entries")
        return
    print(f"Cache verification failed: {data['checked']} entries checked")
    for entry in data["entries"]:
        if entry["ok"]:
            continue
        print(f"{entry['cache_key']}: {', '.join(entry['errors'])}")


def _print_gc(data: dict):
    action = "Would remove" if data["dry_run"] else "Removed"
    print(f"{action} {data['removed_entries']} entries ({data['removed_files']} files).")
    for entry in data["entries"]:
        print(f"{entry['cache_key']}: {entry['reason']}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    cache_dir = Path(args.cache_dir)
    as_json = bool(getattr(args, "json", False))

    try:
        if args.command == "fetch":
            data = fetch_url(
                args.url,
                cache_dir,
                offline=args.offline,
                ttl_seconds=args.ttl_seconds,
                out=args.out,
            )
            payload = _ok_payload(data)
            if as_json:
                print(json.dumps(payload, indent=2))
            else:
                _print_fetch(data)
            return 0

        if args.command == "show":
            data = show_cache(cache_dir, args.url)
            payload = _ok_payload(data)
            if as_json:
                print(json.dumps(payload, indent=2))
            else:
                _print_show(data)
            return 0

        if args.command == "verify":
            data, ok = verify_cache(cache_dir, args.url)
            payload = _ok_payload(data) if ok else _error_payload("cache verification failed", data)
            if as_json:
                print(json.dumps(payload, indent=2))
            else:
                _print_verify(data)
            return 0 if ok else 1

        if args.command == "gc":
            data = gc_cache(cache_dir, args.max_age_seconds, dry_run=args.dry_run)
            payload = _ok_payload(data)
            if as_json:
                print(json.dumps(payload, indent=2))
            else:
                _print_gc(data)
            return 0

        parser.error(f"unknown command: {args.command}")
    except CacheError as exc:
        if as_json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        if as_json:
            print(json.dumps(_error_payload(str(exc)), indent=2))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
