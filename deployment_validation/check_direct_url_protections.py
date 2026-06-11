"""Probe deployed Map-Kinase URL paths that should not expose files.

Usage:
    python deployment_validation/check_direct_url_protections.py https://DOMAIN_HERE

The script treats 403, 404, and common redirect-to-app/login behavior as
protected. A 2xx response with content that looks like a directly served file is
reported as FAIL.
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen


PATHS_TO_PROBE = [
    "/JSONfiles/latest_preview.json",
    "/cache/global_protein_catalog.json",
    "/output/testing_file_001/",
    "/MapKinase_WebApp/JSONfiles/latest_preview.json",
    "/MapKinase_WebApp/cache/global_protein_catalog.json",
    "/stored_pathways/",
    "/.env",
    "/requirements.txt",
    "/DEPENDENCY_SECURITY_REVIEW.md",
]

PROTECTED_STATUS_CODES = {401, 403, 404, 405}
SUSPICIOUS_2XX_CONTENT_TYPES = {
    "application/json",
    "text/plain",
    "text/markdown",
    "application/octet-stream",
}


class _TitleParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_title = False
        self.title_parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "title":
            self.in_title = True

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "title":
            self.in_title = False

    def handle_data(self, data: str) -> None:
        if self.in_title:
            self.title_parts.append(data)

    @property
    def title(self) -> str:
        return " ".join(part.strip() for part in self.title_parts if part.strip())


@dataclass
class ProbeResult:
    path: str
    url: str
    status: str
    content_type: str
    bytes_read: int
    result: str
    notes: str


def _request(url: str, timeout: int) -> tuple[str, str, bytes, str]:
    req = Request(
        url,
        headers={
            "User-Agent": "MapKinase-Deployment-Validation/1.0",
            "Accept": "*/*",
        },
        method="GET",
    )
    try:
        with urlopen(req, timeout=timeout) as response:
            status = str(response.status)
            content_type = str(response.headers.get("content-type") or "")
            body = response.read(4096)
            final_url = str(response.geturl())
            return status, content_type, body, final_url
    except HTTPError as exc:
        body = exc.read(4096)
        return str(exc.code), str(exc.headers.get("content-type") or ""), body, url
    except URLError as exc:
        return "ERROR", "", b"", f"{url} ({exc.reason})"


def _looks_like_app_html(content_type: str, body: bytes) -> tuple[bool, str]:
    lowered = content_type.lower()
    if "text/html" not in lowered:
        return False, ""
    parser = _TitleParser()
    try:
        parser.feed(body.decode("utf-8", errors="replace"))
    except Exception:
        return True, "HTML response"
    title = parser.title
    if title:
        return True, f"HTML response title={title!r}"
    return True, "HTML response"


def _classify(status: str, content_type: str, body: bytes, final_url: str, original_url: str) -> tuple[str, str]:
    try:
        status_code = int(status)
    except ValueError:
        return "REVIEW", final_url

    if status_code in PROTECTED_STATUS_CODES:
        return "PASS", f"Protected status {status_code}"

    redirected = final_url != original_url
    if 300 <= status_code < 400:
        return "PASS", f"Redirect status {status_code}"

    if 200 <= status_code < 300:
        is_html, html_note = _looks_like_app_html(content_type, body)
        if redirected and is_html:
            return "PASS", f"Redirected to app/login HTML: {html_note}"
        if is_html:
            return "REVIEW", f"2xx HTML response; verify this is the app route, not directory/file listing. {html_note}"

        lowered_type = content_type.split(";", 1)[0].strip().lower()
        if lowered_type in SUSPICIOUS_2XX_CONTENT_TYPES:
            return "FAIL", f"2xx direct-looking content type {content_type!r}"
        if body.strip().startswith((b"{", b"[", b"#", b"Uniprot", b"HMDB", b"shiny==")):
            return "FAIL", "2xx body looks like direct file content"
        return "REVIEW", f"2xx non-HTML response with content type {content_type!r}"

    return "REVIEW", f"Unexpected status {status_code}"


def probe(base_url: str, paths: Iterable[str], timeout: int = 10) -> list[ProbeResult]:
    base = base_url.rstrip("/") + "/"
    results: list[ProbeResult] = []
    for path in paths:
        url = urljoin(base, path.lstrip("/"))
        status, content_type, body, final_url = _request(url, timeout=timeout)
        result, notes = _classify(status, content_type, body, final_url, url)
        results.append(
            ProbeResult(
                path=path,
                url=url,
                status=status,
                content_type=content_type,
                bytes_read=len(body),
                result=result,
                notes=notes,
            )
        )
    return results


def main(argv: list[str]) -> int:
    if len(argv) == 2 and argv[1] in {"-h", "--help"}:
        print("Usage: python deployment_validation/check_direct_url_protections.py https://DOMAIN_HERE")
        return 0
    if len(argv) != 2:
        print("Usage: python deployment_validation/check_direct_url_protections.py https://DOMAIN_HERE")
        return 2

    base_url = argv[1].strip()
    if not base_url.startswith(("http://", "https://")):
        print("Base URL must start with http:// or https://", file=sys.stderr)
        return 2

    results = probe(base_url, PATHS_TO_PROBE)
    payload = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "base_url": base_url,
        "results": [asdict(result) for result in results],
    }

    results_dir = Path(__file__).resolve().parent / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_path = results_dir / "url_check_results.json"
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    max_label = max(len(result.path) for result in results)
    for result in results:
        print(f"{result.result:6} {result.status:>5} {result.path:<{max_label}} {result.notes}")

    print(f"\nWrote {output_path}")

    failed = [result for result in results if result.result == "FAIL"]
    review = [result for result in results if result.result == "REVIEW"]
    if failed:
        print(f"FAIL: {len(failed)} URL probe(s) looked publicly exposed.", file=sys.stderr)
        return 1
    if review:
        print(f"REVIEW: {len(review)} URL probe(s) require manual confirmation.", file=sys.stderr)
        return 3
    print("PASS: all probed paths returned protected statuses or safe redirects.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
