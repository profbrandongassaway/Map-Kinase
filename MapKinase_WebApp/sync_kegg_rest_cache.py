#!/usr/bin/env python3
"""
Preload KEGG KGML and PNG pathway files into Map-Kinase's runtime cache.

KEGG REST is available to academic users and permits at most three requests
per second. This command intentionally processes requests serially and defaults
to a 0.4-second interval.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import requests
from PIL import Image


LOGGER = logging.getLogger("sync_kegg_rest_cache")
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parent
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "stored_pathways" / "kegg"
DEFAULT_MANIFEST_PATH = DEFAULT_OUTPUT_DIR / "_rest_sync_manifest.json"
KEGG_API_BASE = "https://rest.kegg.jp"
DEFAULT_REQUEST_INTERVAL = 0.4
MIN_REQUEST_INTERVAL = 1.0 / 3.0


class KeggNotFoundError(RuntimeError):
    """Raised when a requested KEGG resource does not exist."""


class RateLimiter:
    def __init__(self, interval_seconds: float) -> None:
        self.interval = max(float(interval_seconds), 0.0)
        self._last_request_time = 0.0

    def wait(self) -> None:
        now = time.monotonic()
        elapsed = now - self._last_request_time
        if elapsed < self.interval:
            time.sleep(self.interval - elapsed)
        self._last_request_time = time.monotonic()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Preload KEGG KGML and PNG files into stored_pathways/kegg."
    )
    parser.add_argument(
        "--org",
        action="append",
        required=True,
        help="KEGG organism code to sync, such as hsa or mmu. Repeat for multiple organisms.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Destination cache directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--manifest",
        default=str(DEFAULT_MANIFEST_PATH),
        help=f"Output manifest path (default: {DEFAULT_MANIFEST_PATH}).",
    )
    parser.add_argument(
        "--max-pathways",
        type=int,
        default=None,
        help="For testing, process at most N pathways per organism.",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Redownload resources even when valid cache files already exist.",
    )
    parser.add_argument(
        "--no-images",
        action="store_true",
        help="Download KGML only and skip pathway PNG images.",
    )
    parser.add_argument(
        "--request-interval",
        type=float,
        default=DEFAULT_REQUEST_INTERVAL,
        help=(
            "Minimum seconds between KEGG requests "
            f"(default: {DEFAULT_REQUEST_INTERVAL}; minimum: {MIN_REQUEST_INTERVAL:.3f})."
        ),
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Per-request timeout in seconds (default: 30).",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=5,
        help="Maximum attempts for transient request failures (default: 5).",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print the sync manifest JSON.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def validate_args(args: argparse.Namespace) -> List[str]:
    errors: List[str] = []
    orgs = [str(org or "").strip().lower() for org in args.org]
    invalid_orgs = [org for org in orgs if not re.fullmatch(r"[a-z0-9]{2,5}", org)]
    if invalid_orgs:
        errors.append(f"Invalid KEGG organism code(s): {invalid_orgs}")
    if args.max_pathways is not None and args.max_pathways < 1:
        errors.append("--max-pathways must be >= 1")
    if args.request_interval < MIN_REQUEST_INTERVAL:
        errors.append(
            "--request-interval must be at least "
            f"{MIN_REQUEST_INTERVAL:.3f} seconds to respect KEGG's request limit"
        )
    if args.timeout <= 0:
        errors.append("--timeout must be > 0")
    if args.max_retries < 1:
        errors.append("--max-retries must be >= 1")
    return errors


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_bytes_atomic(path: Path, content: bytes) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_bytes(content)
    tmp.replace(path)


def write_json_atomic(path: Path, payload: object, pretty: bool = False) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        if pretty:
            json.dump(payload, fh, indent=2, ensure_ascii=False)
        else:
            json.dump(payload, fh, separators=(",", ":"), ensure_ascii=False)
    tmp.replace(path)


def is_valid_kgml(content: bytes) -> bool:
    if not content.strip():
        return False
    try:
        root = ET.fromstring(content)
    except ET.ParseError:
        return False
    return root.tag.rsplit("}", 1)[-1].lower() == "pathway"


def is_valid_png(content: bytes) -> bool:
    if not content:
        return False
    try:
        with Image.open(BytesIO(content)) as image:
            image.verify()
        return True
    except OSError:
        return False


def is_valid_cached_file(path: Path, resource_type: str) -> bool:
    if not path.is_file() or path.stat().st_size <= 0:
        return False
    try:
        content = path.read_bytes()
    except OSError:
        return False
    if resource_type == "kgml":
        return is_valid_kgml(content)
    if resource_type == "image":
        return is_valid_png(content)
    raise ValueError(f"Unknown resource type: {resource_type}")


def fetch_bytes(
    session: requests.Session,
    url: str,
    rate_limiter: RateLimiter,
    timeout: float,
    max_retries: int,
) -> bytes:
    last_error: Optional[Exception] = None
    for attempt in range(1, max_retries + 1):
        try:
            rate_limiter.wait()
            response = session.get(url, timeout=timeout)
            if response.status_code == 404:
                raise KeggNotFoundError(f"404 for URL: {url}")
            response.raise_for_status()
            return bytes(response.content)
        except KeggNotFoundError:
            raise
        except requests.RequestException as exc:
            last_error = exc
            if attempt >= max_retries:
                break
            delay = min(2 ** (attempt - 1), 10)
            LOGGER.warning(
                "Request failed (%s/%s) for %s: %s. Retrying in %ss",
                attempt,
                max_retries,
                url,
                exc,
                delay,
            )
            time.sleep(delay)
    raise RuntimeError(f"Failed to fetch {url}") from last_error


def parse_pathway_list(content: bytes, org: str) -> List[Dict[str, str]]:
    text = content.decode("utf-8", errors="replace")
    pathways: List[Dict[str, str]] = []
    seen = set()
    expected = re.compile(rf"^{re.escape(org)}\d{{5}}$")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            LOGGER.warning("Skipping malformed KEGG pathway-list line: %r", raw_line)
            continue
        raw_id, name = parts
        pathway_id = raw_id.strip()
        if pathway_id.startswith("path:"):
            pathway_id = pathway_id[5:]
        if not expected.fullmatch(pathway_id) or pathway_id in seen:
            continue
        seen.add(pathway_id)
        pathways.append({"pathway_id": pathway_id, "name": name.strip()})
    return pathways


def sync_resource(
    session: requests.Session,
    rate_limiter: RateLimiter,
    pathway_id: str,
    resource_type: str,
    destination: Path,
    refresh: bool,
    timeout: float,
    max_retries: int,
) -> str:
    if not refresh and is_valid_cached_file(destination, resource_type):
        return "cached"

    endpoint = "kgml" if resource_type == "kgml" else "image"
    url = f"{KEGG_API_BASE}/get/{pathway_id}/{endpoint}"
    content = fetch_bytes(
        session=session,
        url=url,
        rate_limiter=rate_limiter,
        timeout=timeout,
        max_retries=max_retries,
    )
    validator = is_valid_kgml if resource_type == "kgml" else is_valid_png
    if not validator(content):
        raise ValueError(f"KEGG returned invalid {resource_type} content for {pathway_id}")
    write_bytes_atomic(destination, content)
    return "downloaded"


def sync_organism(
    session: requests.Session,
    rate_limiter: RateLimiter,
    org: str,
    output_dir: Path,
    max_pathways: Optional[int],
    refresh: bool,
    include_images: bool,
    timeout: float,
    max_retries: int,
) -> Dict[str, Any]:
    list_url = f"{KEGG_API_BASE}/list/pathway/{org}"
    LOGGER.info("Loading KEGG pathway list for %s", org)
    list_content = fetch_bytes(
        session=session,
        url=list_url,
        rate_limiter=rate_limiter,
        timeout=timeout,
        max_retries=max_retries,
    )
    pathways = parse_pathway_list(list_content, org)
    if max_pathways is not None:
        pathways = pathways[:max_pathways]
    if not pathways:
        raise RuntimeError(f"No KEGG pathways found for organism {org}")

    org_dir = output_dir / org
    ensure_dir(org_dir)
    counts: Dict[str, int] = defaultdict(int)
    failures: List[Dict[str, str]] = []

    for index, pathway in enumerate(pathways, start=1):
        pathway_id = pathway["pathway_id"]
        LOGGER.info("[%s %s/%s] Syncing %s", org, index, len(pathways), pathway_id)
        resources = [
            ("kgml", org_dir / f"{pathway_id}.xml"),
        ]
        if include_images:
            resources.append(("image", org_dir / f"{pathway_id}.png"))

        for resource_type, destination in resources:
            try:
                result = sync_resource(
                    session=session,
                    rate_limiter=rate_limiter,
                    pathway_id=pathway_id,
                    resource_type=resource_type,
                    destination=destination,
                    refresh=refresh,
                    timeout=timeout,
                    max_retries=max_retries,
                )
                counts[f"{resource_type}_{result}"] += 1
            except Exception as exc:  # noqa: BLE001
                counts[f"{resource_type}_failed"] += 1
                failures.append(
                    {
                        "pathway_id": pathway_id,
                        "resource": resource_type,
                        "error": str(exc),
                    }
                )
                LOGGER.warning("Failed %s for %s: %s", resource_type, pathway_id, exc)

    return {
        "pathways_listed": len(pathways),
        "counts": dict(sorted(counts.items())),
        "failures": failures,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    errors = validate_args(args)
    if errors:
        for error in errors:
            LOGGER.error(error)
        return 2

    orgs = list(dict.fromkeys(str(org).strip().lower() for org in args.org))
    output_dir = Path(args.out_dir)
    manifest_path = Path(args.manifest)
    ensure_dir(output_dir)

    session = requests.Session()
    session.headers.update({"User-Agent": "Map-Kinase KEGG cache sync/1.0"})
    rate_limiter = RateLimiter(args.request_interval)
    organism_results: Dict[str, Dict[str, Any]] = {}
    top_level_failures: List[Dict[str, str]] = []

    for org in orgs:
        try:
            organism_results[org] = sync_organism(
                session=session,
                rate_limiter=rate_limiter,
                org=org,
                output_dir=output_dir,
                max_pathways=args.max_pathways,
                refresh=bool(args.refresh),
                include_images=not bool(args.no_images),
                timeout=float(args.timeout),
                max_retries=int(args.max_retries),
            )
        except Exception as exc:  # noqa: BLE001
            top_level_failures.append({"org": org, "error": str(exc)})
            LOGGER.error("Failed to sync organism %s: %s", org, exc)

    resource_failures = sum(
        len(result.get("failures", [])) for result in organism_results.values()
    )
    manifest = {
        "meta": {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "kegg_api_base": KEGG_API_BASE,
            "orgs": orgs,
            "out_dir": str(output_dir),
            "max_pathways_per_org": args.max_pathways,
            "refresh": bool(args.refresh),
            "include_images": not bool(args.no_images),
            "request_interval_seconds": float(args.request_interval),
            "academic_use_notice": (
                "KEGG REST is restricted to academic use and no more than "
                "three requests per second."
            ),
            "summary": {
                "organisms_requested": len(orgs),
                "organisms_completed": len(organism_results),
                "organism_failures": len(top_level_failures),
                "resource_failures": resource_failures,
            },
        },
        "organisms": organism_results,
        "failures": top_level_failures,
    }
    write_json_atomic(manifest_path, manifest, pretty=bool(args.pretty))
    LOGGER.info("Wrote KEGG sync manifest to %s", manifest_path)

    if top_level_failures or resource_failures:
        LOGGER.error(
            "KEGG sync completed with %s organism failure(s) and %s resource failure(s)",
            len(top_level_failures),
            resource_failures,
        )
        return 1
    LOGGER.info("KEGG sync completed successfully for: %s", ", ".join(orgs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
