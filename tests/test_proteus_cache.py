#!/usr/bin/env python3
import datetime as dt
import hashlib
import http.server
import json
import subprocess
import sys
import tempfile
import threading
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import proteus_cache


def scrub_paths(text: str, *paths: Path) -> str:
    scrubbed = text
    for path in paths:
        scrubbed = scrubbed.replace(str(path), "<tmp>")
        try:
            scrubbed = scrubbed.replace(str(path.resolve()), "<tmp>")
        except OSError:
            pass
    return scrubbed


def run_cache(*args: str, check: bool = True, scrub_root: Path | None = None) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        [sys.executable, "scripts/proteus_cache.py", *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    if check and proc.returncode != 0:
        stdout = scrub_paths(proc.stdout, scrub_root) if scrub_root else proc.stdout
        stderr = scrub_paths(proc.stderr, scrub_root) if scrub_root else proc.stderr
        raise AssertionError(f"command failed: {args}\nstdout:\n{stdout}\nstderr:\n{stderr}")
    return proc


class CacheHTTPHandler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        self.server.request_count += 1
        self.send_response(self.server.status)
        self.send_header("Content-Type", self.server.content_type)
        self.send_header("Content-Length", str(len(self.server.body)))
        self.end_headers()
        self.wfile.write(self.server.body)

    def log_message(self, _format, *_args):
        pass


class LocalHTTPServer:
    def __init__(self, body: bytes, content_type: str = "text/plain", status: int = 200):
        self.body = body
        self.content_type = content_type
        self.status = status
        self.server = None
        self.thread = None
        self.url = None

    def __enter__(self):
        self.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), CacheHTTPHandler)
        self.server.body = self.body
        self.server.content_type = self.content_type
        self.server.status = self.status
        self.server.request_count = 0
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        host, port = self.server.server_address
        self.url = f"http://{host}:{port}/payload"
        return self

    def __exit__(self, _exc_type, _exc, _tb):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=5)


def write_entry(cache_dir: Path, url: str, body: bytes, timestamp: str):
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = proteus_cache.cache_key(url)
    body_path = cache_dir / f"{key}.body"
    metadata_path = cache_dir / f"{key}.json"
    body_path.write_bytes(body)
    metadata = {
        "url": url,
        "timestamp": timestamp,
        "status": 200,
        "content_type": "text/plain",
        "body_sha256": hashlib.sha256(body).hexdigest(),
        "byte_count": len(body),
    }
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    return key


class ProteusCacheTests(unittest.TestCase):
    def test_fetch_show_verify_offline_and_out(self):
        body = b'{"ok": true}\n'
        with tempfile.TemporaryDirectory() as tmp, LocalHTTPServer(body, "application/json") as server:
            tmp_path = Path(tmp)
            cache_dir = tmp_path / "cache"
            out_path = tmp_path / "copy.json"

            proc = run_cache("fetch", server.url, "--cache-dir", str(cache_dir), "--json", scrub_root=tmp_path)
            data = json.loads(proc.stdout)
            key = proteus_cache.cache_key(server.url)

            self.assertEqual(data["status"], "ok")
            self.assertFalse(data["cache_hit"])
            self.assertFalse(data["stale"])
            self.assertEqual(data["cache_key"], key)
            self.assertEqual(data["metadata"]["status"], 200)
            self.assertEqual(data["metadata"]["content_type"], "application/json")
            self.assertEqual(data["metadata"]["byte_count"], len(body))
            self.assertEqual((cache_dir / f"{key}.body").read_bytes(), body)
            self.assertEqual(server.server.request_count, 1)

            proc = run_cache("show", server.url, "--cache-dir", str(cache_dir), "--json", scrub_root=tmp_path)
            show_data = json.loads(proc.stdout)
            self.assertEqual(show_data["status"], "ok")
            self.assertEqual(show_data["count"], 1)
            self.assertTrue(show_data["entries"][0]["body_present"])
            self.assertEqual(show_data["entries"][0]["metadata"]["url"], server.url)

            proc = run_cache("verify", server.url, "--cache-dir", str(cache_dir), "--json", scrub_root=tmp_path)
            verify_data = json.loads(proc.stdout)
            self.assertEqual(verify_data["status"], "ok")
            self.assertTrue(verify_data["ok"])
            self.assertEqual(verify_data["checked"], 1)

            proc = run_cache(
                "fetch",
                server.url,
                "--cache-dir",
                str(cache_dir),
                "--offline",
                "--out",
                str(out_path),
                "--json",
                scrub_root=tmp_path,
            )
            offline_data = json.loads(proc.stdout)
            self.assertTrue(offline_data["cache_hit"])
            self.assertEqual(out_path.read_bytes(), body)
            self.assertEqual(server.server.request_count, 1)
            self.assertIn("<tmp>/", offline_data["out"])
            self.assertNotIn(str(tmp_path), proc.stdout)

    def test_ttl_refetches_expired_entries_and_keeps_fresh_entries(self):
        first = b"first"
        second = b"second"
        third = b"third"
        with tempfile.TemporaryDirectory() as tmp, LocalHTTPServer(first) as server:
            tmp_path = Path(tmp)
            cache_dir = tmp_path / "cache"
            run_cache("fetch", server.url, "--cache-dir", str(cache_dir), "--json", scrub_root=tmp_path)
            key = proteus_cache.cache_key(server.url)
            metadata_path = cache_dir / f"{key}.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["timestamp"] = "2000-01-01T00:00:00Z"
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

            server.server.body = second
            proc = run_cache(
                "fetch",
                server.url,
                "--cache-dir",
                str(cache_dir),
                "--ttl-seconds",
                "1",
                "--json",
                scrub_root=tmp_path,
            )
            refreshed = json.loads(proc.stdout)
            self.assertFalse(refreshed["cache_hit"])
            self.assertEqual((cache_dir / f"{key}.body").read_bytes(), second)
            self.assertEqual(server.server.request_count, 2)

            server.server.body = third
            proc = run_cache(
                "fetch",
                server.url,
                "--cache-dir",
                str(cache_dir),
                "--ttl-seconds",
                "86400",
                "--json",
                scrub_root=tmp_path,
            )
            cached = json.loads(proc.stdout)
            self.assertTrue(cached["cache_hit"])
            self.assertEqual((cache_dir / f"{key}.body").read_bytes(), second)
            self.assertEqual(server.server.request_count, 2)

    def test_verify_reports_corrupt_body_and_secret_urls_are_rejected(self):
        url = "https://example.org/public.txt"
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            cache_dir = tmp_path / "cache"
            key = write_entry(cache_dir, url, b"good", proteus_cache._now_timestamp())
            (cache_dir / f"{key}.body").write_bytes(b"bad")

            proc = run_cache(
                "verify",
                url,
                "--cache-dir",
                str(cache_dir),
                "--json",
                check=False,
                scrub_root=tmp_path,
            )
            data = json.loads(proc.stdout)
            self.assertEqual(proc.returncode, 1)
            self.assertEqual(data["status"], "error")
            self.assertFalse(data["ok"])
            self.assertIn("body sha256 mismatch", data["entries"][0]["errors"])

            proc = run_cache(
                "fetch",
                "https://example.org/file?api_key=secret",
                "--cache-dir",
                str(cache_dir),
                "--json",
                check=False,
                scrub_root=tmp_path,
            )
            secret_data = json.loads(proc.stdout)
            self.assertEqual(proc.returncode, 1)
            self.assertEqual(secret_data["status"], "error")
            self.assertIn("secret-like query parameter", secret_data["error"])

    def test_gc_removes_expired_entries_and_orphans(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            cache_dir = tmp_path / "cache"
            expired_key = write_entry(cache_dir, "https://example.org/old.txt", b"old", "2000-01-01T00:00:00Z")
            fresh_time = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
            fresh_key = write_entry(cache_dir, "https://example.org/fresh.txt", b"fresh", fresh_time)
            orphan_key = "f" * 64
            (cache_dir / f"{orphan_key}.body").write_bytes(b"orphan")

            proc = run_cache(
                "gc",
                "--cache-dir",
                str(cache_dir),
                "--max-age-seconds",
                "86400",
                "--json",
                scrub_root=tmp_path,
            )
            data = json.loads(proc.stdout)
            removed_keys = {entry["cache_key"] for entry in data["entries"]}

            self.assertEqual(data["status"], "ok")
            self.assertEqual(data["removed_entries"], 2)
            self.assertEqual(data["removed_files"], 3)
            self.assertEqual(removed_keys, {expired_key, orphan_key})
            self.assertFalse((cache_dir / f"{expired_key}.body").exists())
            self.assertFalse((cache_dir / f"{expired_key}.json").exists())
            self.assertFalse((cache_dir / f"{orphan_key}.body").exists())
            self.assertTrue((cache_dir / f"{fresh_key}.body").exists())
            self.assertTrue((cache_dir / f"{fresh_key}.json").exists())


if __name__ == "__main__":
    unittest.main()
