#!/usr/bin/env python3
"""
Audit KEGG and WikiPathways downloads across configured organisms and write a
ban-list manifest for pathways that fail to download usable content.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import requests

CURRENT_DIR = Path(__file__).resolve().parent
PARENT_DIR = CURRENT_DIR.parent
if str(PARENT_DIR) not in sys.path:
    sys.path.insert(0, str(PARENT_DIR))

from MapKinase_WebApp.pathway_banlist import DEFAULT_BANNED_PATHWAYS_FILE, SCHEMA_VERSION


LOGGER = logging.getLogger("build_banned_pathways")
DEFAULT_RATE_LIMIT = 0.25
KEGG_API_BASE = "https://rest.kegg.jp"
WIKIPATHWAYS_JSON_LIST_URL = "https://www.wikipathways.org/json/listPathways.json"
BASE_DIR = CURRENT_DIR
PROJECT_ROOT = BASE_DIR.parent
SPECIES_REF_FILE = BASE_DIR / "species_ref_list.csv"
KEGG_PATHWAYS_FILE = BASE_DIR / "kegg_pathways.txt"
STORED_PATHWAYS_DIR = PROJECT_ROOT / "stored_pathways"
_PYWIKIPATHWAYS_CLIENT: Optional[Tuple[Any, Any]] = None
_PYWIKIPATHWAYS_IMPORT_ERROR: Optional[BaseException] = None


class RateLimiter:
    def __init__(self, interval_seconds: float) -> None:
        self.interval = max(float(interval_seconds), 0.0)
        self._last_request_time = 0.0

    def wait(self) -> None:
        if self.interval <= 0:
            return
        now = time.monotonic()
        elapsed = now - self._last_request_time
        if elapsed < self.interval:
            time.sleep(self.interval - elapsed)
        self._last_request_time = time.monotonic()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a banned-pathway manifest for KEGG and WikiPathways."
    )
    parser.add_argument(
        "--source",
        choices=["all", "kegg", "wikipathways"],
        default="all",
        help="Which pathway source to audit (default: all).",
    )
    parser.add_argument(
        "--org",
        action="append",
        default=[],
        help="Restrict to one or more organism codes such as hsa or mmu.",
    )
    parser.add_argument(
        "--out",
        default=str(DEFAULT_BANNED_PATHWAYS_FILE),
        help=f"Output manifest path (default: {DEFAULT_BANNED_PATHWAYS_FILE}).",
    )
    parser.add_argument(
        "--max-pathways",
        type=int,
        default=None,
        help="For debugging, audit at most N pathways per source/org.",
    )
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=DEFAULT_RATE_LIMIT,
        help=f"Seconds between remote requests (default: {DEFAULT_RATE_LIMIT}).",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print the output JSON.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json_atomic(path: Path, obj: object, pretty: bool = False) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        if pretty:
            json.dump(obj, fh, indent=2, ensure_ascii=False)
        else:
            json.dump(obj, fh, separators=(",", ":"), ensure_ascii=False)
    tmp.replace(path)


def normalize_species_folder(species: str) -> str:
    text = "".join(ch.lower() if ch.isalnum() else "_" for ch in str(species or "").strip())
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_") or "unknown"


def normalize_species_name(value: Any) -> str:
    return str(value or "").strip().lower().replace("_", " ")


def load_species_rows() -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not SPECIES_REF_FILE.exists():
        return rows
    with SPECIES_REF_FILE.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            cleaned = {(str(k or "").lstrip("\ufeff").strip()): str(v or "").strip() for k, v in dict(raw or {}).items()}
            org = cleaned.get("Kegg Gene ID", "").lower()
            species = cleaned.get("Species", "")
            label = cleaned.get("Common Name", "") or species or org
            if not org or not species:
                continue
            rows.append({"org": org, "species": species, "label": label})
    return rows


def load_kegg_pathway_rows() -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not KEGG_PATHWAYS_FILE.exists():
        return rows
    with KEGG_PATHWAYS_FILE.open("r", encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("Pathway_ID"):
                continue
            parts = stripped.split("\t")
            if len(parts) < 2:
                continue
            raw_id = parts[0].strip()
            name = parts[1].strip()
            digits = raw_id[3:] if raw_id.lower().startswith("map") and len(raw_id) > 3 else raw_id
            if not digits:
                continue
            rows.append({"raw_id": raw_id, "digits": digits, "name": name})
    return rows


def validate_xml_text(text: str, expected_root_fragment: str = "") -> Tuple[bool, str]:
    if not str(text or "").strip():
        return False, "empty response body"
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        return False, f"XML parse error: {exc}"
    if expected_root_fragment and expected_root_fragment not in str(root.tag or ""):
        return False, f"unexpected root tag: {root.tag}"
    return True, ""


def request_text(
    session: requests.Session,
    url: str,
    rate_limiter: RateLimiter,
    timeout: int = 30,
    max_retries: int = 3,
) -> Tuple[Optional[str], Optional[str]]:
    last_error: Optional[str] = None
    for attempt in range(1, max_retries + 1):
        try:
            rate_limiter.wait()
            response = session.get(url, timeout=timeout)
            if response.status_code == 404:
                return None, "HTTP 404"
            response.raise_for_status()
            return response.text, None
        except requests.RequestException as exc:
            last_error = str(exc)
            if attempt < max_retries:
                time.sleep(min(2 ** (attempt - 1), 5))
    return None, last_error or "request failed"


def request_json(
    session: requests.Session,
    url: str,
    rate_limiter: RateLimiter,
    timeout: int = 30,
    max_retries: int = 3,
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    text, error = request_text(
        session=session,
        url=url,
        rate_limiter=rate_limiter,
        timeout=timeout,
        max_retries=max_retries,
    )
    if error is not None:
        return None, error
    try:
        payload = json.loads(text or "")
    except json.JSONDecodeError as exc:
        return None, f"JSON decode error: {exc}"
    if not isinstance(payload, dict):
        return None, "JSON payload was not an object"
    return payload, None


def parse_wikipathways_list_for_species(payload: Dict[str, Any], species: str) -> List[Dict[str, str]]:
    organisms = payload.get("organisms")
    if not isinstance(organisms, list):
        return []
    species_norm = normalize_species_name(species)
    rows: List[Dict[str, str]] = []
    seen: set[str] = set()
    for organism_row in organisms:
        if not isinstance(organism_row, dict):
            continue
        organism_name = str(organism_row.get("latin") or organism_row.get("common") or "").strip()
        organism_norm = normalize_species_name(organism_name)
        common_norm = normalize_species_name(organism_row.get("common"))
        if species_norm and species_norm not in {organism_norm, common_norm}:
            continue
        pathways = organism_row.get("pathways")
        if not isinstance(pathways, list):
            continue
        for raw in pathways:
            if not isinstance(raw, dict):
                continue
            pathway_id = str(raw.get("id") or raw.get("pathway_id") or "").strip().upper()
            if not pathway_id or pathway_id in seen:
                continue
            seen.add(pathway_id)
            rows.append(
                {
                    "pathway_id": pathway_id,
                    "name": str(raw.get("name") or pathway_id).strip(),
                    "species": str(raw.get("species") or organism_name).strip(),
                }
            )
        break
    return rows


def find_cached_wikipathways_gpml(pathway_id: str, species: str) -> Optional[Path]:
    species_dir = STORED_PATHWAYS_DIR / "wikipathways" / normalize_species_folder(species)
    candidate = species_dir / f"{pathway_id}.gpml"
    if candidate.exists():
        return candidate
    root = STORED_PATHWAYS_DIR / "wikipathways"
    if not root.exists():
        return None
    matches = list(root.rglob(f"{pathway_id}.gpml"))
    return matches[0] if matches else None


def find_local_wikipathways_gpml(pathway_id: str) -> Optional[Path]:
    candidates = [
        Path(os.getcwd()) / f"{pathway_id}.gpml",
        Path(os.path.dirname(os.getcwd())) / f"{pathway_id}.gpml",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def extract_species_from_gpml(gpml_text: str) -> str:
    try:
        root = ET.fromstring(gpml_text)
    except ET.ParseError:
        return ""
    for key in ("Organism", "organism", "Species", "species"):
        value = root.attrib.get(key)
        if str(value or "").strip():
            return str(value).strip()
    return ""


def get_pywikipathways_client() -> Tuple[Optional[Any], Optional[Any], Optional[str]]:
    global _PYWIKIPATHWAYS_CLIENT, _PYWIKIPATHWAYS_IMPORT_ERROR
    if _PYWIKIPATHWAYS_CLIENT is not None:
        return _PYWIKIPATHWAYS_CLIENT[0], _PYWIKIPATHWAYS_CLIENT[1], None
    if _PYWIKIPATHWAYS_IMPORT_ERROR is not None:
        return None, None, f"{type(_PYWIKIPATHWAYS_IMPORT_ERROR).__name__}: {_PYWIKIPATHWAYS_IMPORT_ERROR}"
    try:
        from pywikipathways import get_pathway, get_pathway_info  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on runtime environment
        _PYWIKIPATHWAYS_IMPORT_ERROR = exc
        return None, None, f"{type(exc).__name__}: {exc}"
    _PYWIKIPATHWAYS_CLIENT = (get_pathway, get_pathway_info)
    return get_pathway, get_pathway_info, None


def fetch_wikipathways_gpml_like_mapkinase(pathway_id: str, species_hint: str) -> Tuple[Optional[Path], Optional[str]]:
    cached_path = find_cached_wikipathways_gpml(pathway_id, species_hint)
    if cached_path and cached_path.exists():
        return cached_path, None

    local_path = find_local_wikipathways_gpml(pathway_id)
    if local_path and local_path.exists():
        return local_path, None

    get_pathway, get_pathway_info, import_error = get_pywikipathways_client()
    if get_pathway is None:
        return None, import_error or "pywikipathways import failed"

    try:
        gpml_text = get_pathway(pathway_id)
    except Exception as exc:
        debug_reason = ""
        if get_pathway_info is not None:
            try:
                info = get_pathway_info(pathway_id)
                if isinstance(info, dict) and info:
                    info_keys = ", ".join(sorted(str(key) for key in info.keys()))
                    debug_reason = f"; get_pathway_info keys: {info_keys}"
            except Exception as info_exc:
                debug_reason = f"; get_pathway_info failed: {type(info_exc).__name__}: {info_exc}"
        return None, f"{type(exc).__name__}: {exc}{debug_reason}"

    species_label = extract_species_from_gpml(gpml_text or "") or str(species_hint or "").strip()
    if not species_label and get_pathway_info is not None:
        try:
            info = get_pathway_info(pathway_id)
            if isinstance(info, dict):
                species_label = str(info.get("species") or info.get("organism") or "").strip()
        except Exception:
            species_label = ""

    species_folder = normalize_species_folder(species_label)
    target_dir = STORED_PATHWAYS_DIR / "wikipathways" / species_folder
    ensure_dir(target_dir)
    file_path = target_dir / f"{pathway_id}.gpml"
    file_path.write_text(gpml_text, encoding="utf-8")
    return file_path, None


def can_audit_wikipathways(
    payload: Dict[str, Any],
    species_rows: Sequence[Dict[str, str]],
    max_pathways: Optional[int],
) -> Tuple[bool, Optional[str]]:
    get_pathway, _, import_error = get_pywikipathways_client()
    if get_pathway is not None:
        return True, None

    for species_row in species_rows:
        species = str(species_row.get("species") or "").strip()
        org = str(species_row.get("org") or "").strip().lower()
        if not species or not org:
            continue
        pathway_rows = parse_wikipathways_list_for_species(payload or {}, species)
        candidates = pathway_rows[:max_pathways] if max_pathways is not None else pathway_rows
        for row in candidates:
            pathway_id = str(row.get("pathway_id") or "").strip().upper()
            if not pathway_id:
                continue
            if find_cached_wikipathways_gpml(pathway_id, species) is not None:
                continue
            if find_local_wikipathways_gpml(pathway_id) is not None:
                continue
            return False, import_error or "pywikipathways import failed"
    return True, None


def probe_kegg_pathway(
    session: requests.Session,
    rate_limiter: RateLimiter,
    org: str,
    pathway_row: Dict[str, str],
) -> Tuple[bool, Optional[Dict[str, Any]]]:
    digits = str(pathway_row.get("digits") or "").strip()
    name = str(pathway_row.get("name") or "").strip()
    raw_id = str(pathway_row.get("raw_id") or "").strip()
    pathway_id = f"{org}{digits}"
    cached_path = STORED_PATHWAYS_DIR / "kegg" / org / f"{pathway_id}.xml"
    if cached_path.exists() and cached_path.stat().st_size > 0:
        valid, error = validate_xml_text(cached_path.read_text(encoding="utf-8", errors="replace"), "pathway")
        if valid:
            return True, None
        LOGGER.debug("Cached KEGG XML failed validation for %s: %s", pathway_id, error)

    url = f"{KEGG_API_BASE}/get/{pathway_id}/kgml"
    xml_text, request_error = request_text(session=session, url=url, rate_limiter=rate_limiter)
    if request_error is not None:
        if request_error != "HTTP 404":
            LOGGER.debug("Ignoring non-404 KEGG fetch error for %s: %s", pathway_id, request_error)
            return True, None
        return False, {
            "source": "kegg",
            "org": org,
            "pathway_id": pathway_id,
            "raw_id": raw_id,
            "name": name,
            "reason": request_error,
        }
    valid, validation_error = validate_xml_text(xml_text or "", "pathway")
    if valid:
        return True, None
    LOGGER.debug("Ignoring non-404 KEGG validation issue for %s: %s", pathway_id, validation_error)
    return True, None


def probe_wikipathways_pathway(
    session: requests.Session,
    rate_limiter: RateLimiter,
    org: str,
    species: str,
    pathway_row: Dict[str, str],
) -> Tuple[bool, Optional[Dict[str, Any]]]:
    pathway_id = str(pathway_row.get("pathway_id") or "").strip().upper()
    name = str(pathway_row.get("name") or "").strip()
    file_path, request_error = fetch_wikipathways_gpml_like_mapkinase(pathway_id=pathway_id, species_hint=species)
    if request_error is not None or file_path is None:
        return False, {
            "source": "wikipathways",
            "org": org,
            "species": species,
            "pathway_id": pathway_id,
            "name": name,
            "reason": request_error or "download failed",
        }
    if not file_path.exists() or file_path.stat().st_size <= 0:
        return False, {
            "source": "wikipathways",
            "org": org,
            "species": species,
            "pathway_id": pathway_id,
            "name": name,
            "reason": "downloaded GPML file was missing or empty",
        }
    valid, validation_error = validate_xml_text(
        file_path.read_text(encoding="utf-8", errors="replace"),
        "Pathway",
    )
    if valid:
        return True, None
    return False, {
        "source": "wikipathways",
        "org": org,
        "species": species,
        "pathway_id": pathway_id,
        "name": name,
        "reason": validation_error,
    }


def audit_kegg(
    session: requests.Session,
    rate_limiter: RateLimiter,
    species_rows: Sequence[Dict[str, str]],
    max_pathways: Optional[int],
) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    banned: List[Dict[str, Any]] = []
    summary: Dict[str, int] = defaultdict(int)
    pathway_rows = load_kegg_pathway_rows()
    for species_row in species_rows:
        org = str(species_row.get("org") or "").strip().lower()
        if not org:
            continue
        candidates = pathway_rows[:max_pathways] if max_pathways is not None else pathway_rows
        LOGGER.info("Auditing KEGG for %s (%s pathways)", org, len(candidates))
        for row in candidates:
            ok, failure = probe_kegg_pathway(session=session, rate_limiter=rate_limiter, org=org, pathway_row=row)
            summary[f"kegg:{org}:checked"] += 1
            if ok:
                continue
            summary[f"kegg:{org}:banned"] += 1
            banned.append(failure or {})
    return banned, dict(summary)


def audit_wikipathways(
    session: requests.Session,
    rate_limiter: RateLimiter,
    species_rows: Sequence[Dict[str, str]],
    max_pathways: Optional[int],
) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    banned: List[Dict[str, Any]] = []
    summary: Dict[str, int] = defaultdict(int)
    payload, request_error = request_json(
        session=session,
        url=WIKIPATHWAYS_JSON_LIST_URL,
        rate_limiter=rate_limiter,
    )
    if request_error is not None:
        LOGGER.error("Failed to load WikiPathways list: %s", request_error)
        summary["wikipathways:list_error"] += 1
        return banned, dict(summary)
    runtime_ok, runtime_error = can_audit_wikipathways(
        payload=payload or {},
        species_rows=species_rows,
        max_pathways=max_pathways,
    )
    if not runtime_ok:
        LOGGER.error(
            "Cannot audit WikiPathways in this environment: %s. "
            "pywikipathways is required for uncached pathways because the auditor now mirrors Map-Kinase's loader.",
            runtime_error,
        )
        summary["wikipathways:runtime_error"] += 1
        return banned, dict(summary)
    for species_row in species_rows:
        org = str(species_row.get("org") or "").strip().lower()
        species = str(species_row.get("species") or "").strip()
        if not org or not species:
            continue
        pathway_rows = parse_wikipathways_list_for_species(payload or {}, species)
        candidates = pathway_rows[:max_pathways] if max_pathways is not None else pathway_rows
        LOGGER.info("Auditing WikiPathways for %s (%s / %s pathways)", org, len(candidates), len(pathway_rows))
        for row in candidates:
            ok, failure = probe_wikipathways_pathway(
                session=session,
                rate_limiter=rate_limiter,
                org=org,
                species=species,
                pathway_row=row,
            )
            summary[f"wikipathways:{org}:checked"] += 1
            if ok:
                continue
            summary[f"wikipathways:{org}:banned"] += 1
            banned.append(failure or {})
    return banned, dict(summary)


def merge_summary(dest: Dict[str, int], src: Dict[str, int]) -> Dict[str, int]:
    merged = dict(dest)
    for key, value in src.items():
        merged[key] = int(merged.get(key, 0)) + int(value)
    return merged


def load_existing_manifest(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {"meta": {}, "entries": []}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {"meta": {}, "entries": []}
    if not isinstance(payload, dict):
        return {"meta": {}, "entries": []}
    if not isinstance(payload.get("meta"), dict):
        payload["meta"] = {}
    if not isinstance(payload.get("entries"), list):
        payload["entries"] = []
    return payload


def build_audited_scope(source_names: Sequence[str], species_rows: Sequence[Dict[str, str]]) -> set[tuple[str, str]]:
    scope: set[tuple[str, str]] = set()
    orgs = {
        str(row.get("org") or "").strip().lower()
        for row in species_rows
        if str(row.get("org") or "").strip()
    }
    for source_name in source_names:
        source_norm = str(source_name or "").strip().lower()
        if not source_norm:
            continue
        for org in orgs:
            scope.add((source_norm, org))
    return scope


def merge_entries_for_scopes(
    existing_entries: Sequence[Dict[str, Any]],
    new_entries: Sequence[Dict[str, Any]],
    audited_scope: set[tuple[str, str]],
) -> List[Dict[str, Any]]:
    merged: List[Dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()

    def _entry_key(entry: Dict[str, Any]) -> tuple[str, str, str]:
        return (
            str(entry.get("source") or "").strip().lower(),
            str(entry.get("org") or "").strip().lower(),
            str(entry.get("pathway_id") or "").strip(),
        )

    for raw in existing_entries:
        if not isinstance(raw, dict):
            continue
        scope_key = (
            str(raw.get("source") or "").strip().lower(),
            str(raw.get("org") or "").strip().lower(),
        )
        if scope_key in audited_scope:
            continue
        entry_key = _entry_key(raw)
        if entry_key in seen:
            continue
        seen.add(entry_key)
        merged.append(dict(raw))

    for raw in new_entries:
        if not isinstance(raw, dict):
            continue
        entry_key = _entry_key(raw)
        if entry_key in seen:
            continue
        seen.add(entry_key)
        merged.append(dict(raw))

    return merged


def merge_manifest_summary(
    existing_summary: Dict[str, Any],
    new_summary: Dict[str, int],
    source_names: Sequence[str],
    species_rows: Sequence[Dict[str, str]],
) -> Dict[str, int]:
    merged: Dict[str, int] = {}
    remove_prefixes: List[str] = []
    orgs = [
        str(row.get("org") or "").strip().lower()
        for row in species_rows
        if str(row.get("org") or "").strip()
    ]
    for source_name in source_names:
        source_norm = str(source_name or "").strip().lower()
        if not source_norm:
            continue
        remove_prefixes.append(f"{source_norm}:")
        for org in orgs:
            remove_prefixes.append(f"{source_norm}:{org}:")

    for key, value in dict(existing_summary or {}).items():
        if not isinstance(value, int):
            continue
        key_text = str(key or "")
        if any(key_text.startswith(prefix) for prefix in remove_prefixes):
            continue
        merged[key_text] = value

    for key, value in dict(new_summary or {}).items():
        merged[str(key or "")] = int(value)
    return merged


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    selected_orgs = {str(item or "").strip().lower() for item in list(args.org or []) if str(item or "").strip()}
    species_rows = load_species_rows()
    if selected_orgs:
        species_rows = [row for row in species_rows if str(row.get("org") or "").strip().lower() in selected_orgs]
    if not species_rows:
        LOGGER.error("No species rows available to audit.")
        return 2

    session = requests.Session()
    session.headers.update({"User-Agent": "build_banned_pathways.py/1.0"})
    rate_limiter = RateLimiter(args.rate_limit)

    selected_sources = ["kegg", "wikipathways"] if args.source == "all" else [str(args.source)]
    out_path = Path(args.out)
    existing_payload = load_existing_manifest(out_path)

    all_entries: List[Dict[str, Any]] = []
    summary: Dict[str, int] = {}

    if args.source in {"all", "kegg"}:
        kegg_entries, kegg_summary = audit_kegg(
            session=session,
            rate_limiter=rate_limiter,
            species_rows=species_rows,
            max_pathways=args.max_pathways,
        )
        all_entries.extend(kegg_entries)
        summary = merge_summary(summary, kegg_summary)

    if args.source in {"all", "wikipathways"}:
        wp_entries, wp_summary = audit_wikipathways(
            session=session,
            rate_limiter=rate_limiter,
            species_rows=species_rows,
            max_pathways=args.max_pathways,
        )
        all_entries.extend(wp_entries)
        summary = merge_summary(summary, wp_summary)

    if any(
        key in summary and int(summary.get(key, 0)) > 0
        for key in ("wikipathways:list_error", "wikipathways:runtime_error")
    ):
        LOGGER.error("Aborting without writing output because the WikiPathways audit did not complete successfully.")
        return 2

    audited_scope = build_audited_scope(selected_sources, species_rows)
    final_entries = merge_entries_for_scopes(
        existing_entries=list(existing_payload.get("entries") or []),
        new_entries=all_entries,
        audited_scope=audited_scope,
    )
    final_entries.sort(
        key=lambda item: (
            str(item.get("source") or ""),
            str(item.get("org") or ""),
            str(item.get("pathway_id") or ""),
        )
    )
    existing_meta = dict(existing_payload.get("meta") or {})
    existing_summary = dict(existing_meta.get("summary") or {})
    final_summary = merge_manifest_summary(
        existing_summary=existing_summary,
        new_summary=summary,
        source_names=selected_sources,
        species_rows=species_rows,
    )
    prior_sources = [str(item or "").strip().lower() for item in list(existing_meta.get("sources") or []) if str(item or "").strip()]
    final_sources = list(dict.fromkeys(prior_sources + selected_sources))
    prior_orgs = [str(item or "").strip().lower() for item in list(existing_meta.get("orgs") or []) if str(item or "").strip()]
    current_orgs = [str(row.get("org") or "").strip().lower() for row in species_rows if str(row.get("org") or "").strip()]
    final_orgs = list(dict.fromkeys(prior_orgs + current_orgs))

    output = {
        "meta": {
            "schema_version": SCHEMA_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "sources": final_sources,
            "orgs": final_orgs,
            "rate_limit": float(args.rate_limit),
            "max_pathways": args.max_pathways,
            "summary": final_summary,
        },
        "entries": final_entries,
    }

    write_json_atomic(out_path, output, pretty=bool(args.pretty))
    LOGGER.info(
        "Wrote %s banned pathway entries to %s (%s newly audited entries, %s preserved)",
        len(final_entries),
        out_path,
        len(all_entries),
        max(len(final_entries) - len(all_entries), 0),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
