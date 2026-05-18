from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set


SCHEMA_VERSION = 1
DEFAULT_BANNED_PATHWAYS_FILE = Path(__file__).resolve().parent / "index_files" / "pathway_banned_list.json"


def _empty_payload() -> Dict[str, Any]:
    return {
        "meta": {
            "schema_version": SCHEMA_VERSION,
        },
        "entries": [],
    }


def _normalize_source(value: Any) -> str:
    return str(value or "").strip().lower()


def _normalize_org(value: Any) -> str:
    return str(value or "").strip().lower()


def normalize_pathway_id_for_source(source: str, pathway_id: Any) -> str:
    text = str(pathway_id or "").strip()
    if not text:
        return ""
    if _normalize_source(source) == "wikipathways":
        return text.upper()
    return text.lower()


@lru_cache(maxsize=4)
def load_banned_pathway_payload(path: Optional[str] = None) -> Dict[str, Any]:
    payload_path = Path(path) if path else DEFAULT_BANNED_PATHWAYS_FILE
    if not payload_path.exists():
        return _empty_payload()
    try:
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
    except Exception:
        return _empty_payload()
    if not isinstance(payload, dict):
        return _empty_payload()
    entries = payload.get("entries")
    if not isinstance(entries, list):
        payload["entries"] = []
    return payload


def clear_banned_pathway_cache() -> None:
    load_banned_pathway_payload.cache_clear()


def iter_banned_entries(
    source: str,
    org: str = "",
    payload: Optional[Dict[str, Any]] = None,
) -> Iterable[Dict[str, Any]]:
    resolved = payload if isinstance(payload, dict) else load_banned_pathway_payload()
    source_norm = _normalize_source(source)
    org_norm = _normalize_org(org)
    for raw in list(resolved.get("entries") or []):
        if not isinstance(raw, dict):
            continue
        if _normalize_source(raw.get("source")) != source_norm:
            continue
        if org_norm and _normalize_org(raw.get("org")) != org_norm:
            continue
        yield raw


def get_banned_pathway_ids(
    source: str,
    org: str = "",
    payload: Optional[Dict[str, Any]] = None,
) -> Set[str]:
    blocked: Set[str] = set()
    for entry in iter_banned_entries(source=source, org=org, payload=payload):
        pathway_id = normalize_pathway_id_for_source(source, entry.get("pathway_id"))
        if pathway_id:
            blocked.add(pathway_id)
    return blocked


def is_pathway_banned(
    source: str,
    pathway_id: Any,
    org: str = "",
    payload: Optional[Dict[str, Any]] = None,
) -> bool:
    normalized = normalize_pathway_id_for_source(source, pathway_id)
    if not normalized:
        return False
    return normalized in get_banned_pathway_ids(source=source, org=org, payload=payload)


def filter_kegg_options_for_species(
    options: Iterable[Dict[str, Any]],
    species_code: str,
    payload: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    code = _normalize_org(species_code)
    blocked = get_banned_pathway_ids(source="kegg", org=code, payload=payload)
    if not blocked:
        return [dict(opt) for opt in options if isinstance(opt, dict)]
    filtered: List[Dict[str, Any]] = []
    for raw in options:
        if not isinstance(raw, dict):
            continue
        digits = str(raw.get("digits") or "").strip()
        if not digits:
            continue
        pathway_id = normalize_pathway_id_for_source("kegg", f"{code}{digits}")
        if pathway_id in blocked:
            continue
        filtered.append(dict(raw))
    return filtered


def filter_wikipathways_options(
    options: Iterable[Dict[str, Any]],
    org: str,
    payload: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    org_norm = _normalize_org(org)
    blocked = get_banned_pathway_ids(source="wikipathways", org=org_norm, payload=payload)
    if not blocked:
        return [dict(opt) for opt in options if isinstance(opt, dict)]
    filtered: List[Dict[str, Any]] = []
    for raw in options:
        if not isinstance(raw, dict):
            continue
        pathway_id = normalize_pathway_id_for_source("wikipathways", raw.get("id"))
        if pathway_id and pathway_id in blocked:
            continue
        filtered.append(dict(raw))
    return filtered
