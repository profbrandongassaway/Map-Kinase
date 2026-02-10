#!/usr/bin/env python3
"""
Build a WikiPathways index for an organism using WikiPathways web services + GPML.

This script mirrors build_kegg_index.py behavior:
- fetches pathway list and GPML with caching/retry
- parses nodes and interactions
- precomputes 1-hop / 2-hop candidate pairs
- resolves native IDs to UniProt using organism id-mapping table
- writes one dataset-independent index JSON
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple
from urllib.parse import quote
import xml.etree.ElementTree as ET

import requests


SCHEMA_VERSION = 1
PARSER_VERSION = 1
WIKIPATHWAYS_API_BASE = "https://webservice.wikipathways.org"
DEFAULT_RATE_LIMIT = 0.25

LOGGER = logging.getLogger("build_wikipathways_index")

# IDE-friendly config:
CONFIG_ORG = "sce"
CONFIG_OUT = f"MapKinase_WebApp/index_files/wikipathways_index_{CONFIG_ORG}.json"
CONFIG_CACHE = ".wikipathways_cache"
CONFIG_MAX_PATHWAYS = None
CONFIG_RATE_LIMIT = DEFAULT_RATE_LIMIT
CONFIG_PRETTY = False
CONFIG_LOG_LEVEL = "INFO"
# Optional species-name override for listPathways endpoint.
CONFIG_SPECIES_NAME = None
# If None, default is MapKinase_WebApp/annotation_files/{org}_mapping_table.txt.
CONFIG_ID_MAPPING_TABLE = None


ORG_TO_SPECIES_DEFAULTS: Dict[str, str] = {
    "hsa": "Homo sapiens",
    "mmu": "Mus musculus",
    "rno": "Rattus norvegicus",
    "dme": "Drosophila melanogaster",
    "sce": "Saccharomyces cerevisiae",
}


class WikiPathwaysFetchError(RuntimeError):
    """Raised when WikiPathways content cannot be fetched after retries."""


class WikiPathwaysNotFoundError(WikiPathwaysFetchError):
    """Raised when WikiPathways returns a missing payload for pathway."""


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


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a cached WikiPathways index for an organism."
    )
    parser.add_argument("--org", required=False, help='Organism code (e.g. "hsa", "mmu").')
    parser.add_argument("--out", required=False, help="Output JSON file path.")
    parser.add_argument(
        "--cache",
        default=".wikipathways_cache",
        help='Cache directory (default: ".wikipathways_cache").',
    )
    parser.add_argument(
        "--species-name",
        default=None,
        help="Optional species name override used by listPathways (e.g. 'Homo sapiens').",
    )
    parser.add_argument(
        "--id-mapping-table",
        default=None,
        help=(
            "Tab-delimited ID mapping table. Default: "
            "MapKinase_WebApp/annotation_files/{org}_mapping_table.txt"
        ),
    )
    parser.add_argument(
        "--max-pathways",
        type=int,
        default=None,
        help="For debugging, process at most N pathways.",
    )
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=DEFAULT_RATE_LIMIT,
        help=f"Seconds between requests (default: {DEFAULT_RATE_LIMIT}).",
    )
    parser.add_argument("--pretty", action="store_true", help="Pretty-print output JSON.")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser.parse_args(argv)


def get_runtime_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) == 0:
        if not CONFIG_ORG or not CONFIG_OUT:
            raise ValueError("CONFIG_ORG and CONFIG_OUT must be set for IDE mode.")
        return argparse.Namespace(
            org=CONFIG_ORG,
            out=CONFIG_OUT,
            cache=CONFIG_CACHE,
            species_name=CONFIG_SPECIES_NAME,
            id_mapping_table=CONFIG_ID_MAPPING_TABLE,
            max_pathways=CONFIG_MAX_PATHWAYS,
            rate_limit=CONFIG_RATE_LIMIT,
            pretty=CONFIG_PRETTY,
            log_level=CONFIG_LOG_LEVEL,
        )

    args = parse_args(argv)
    if not args.org or not args.out:
        raise ValueError("When using CLI args, both --org and --out are required.")
    return args


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text_atomic(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def write_json_atomic(path: Path, obj: object, pretty: bool = False) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        if pretty:
            json.dump(obj, fh, indent=2, ensure_ascii=False, sort_keys=False)
        else:
            json.dump(obj, fh, separators=(",", ":"), ensure_ascii=False, sort_keys=False)
    tmp.replace(path)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def fetch_json(
    session: requests.Session,
    url: str,
    rate_limiter: RateLimiter,
    timeout: int = 30,
    max_retries: int = 5,
) -> Dict[str, object]:
    last_error: Optional[Exception] = None
    for attempt in range(1, max_retries + 1):
        try:
            rate_limiter.wait()
            response = session.get(url, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, json.JSONDecodeError) as exc:
            last_error = exc
            wait_s = min(2 ** (attempt - 1), 10)
            LOGGER.warning(
                "Request failed (%s/%s) for %s: %s. Retrying in %.1fs",
                attempt,
                max_retries,
                url,
                exc,
                wait_s,
            )
            time.sleep(wait_s)
    raise WikiPathwaysFetchError(f"Failed to fetch JSON: {url}") from last_error


def load_or_fetch_json(
    session: requests.Session,
    url: str,
    cache_path: Path,
    rate_limiter: RateLimiter,
) -> Dict[str, object]:
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return json.loads(read_text(cache_path))
    payload = fetch_json(session=session, url=url, rate_limiter=rate_limiter)
    write_json_atomic(cache_path, payload, pretty=True)
    return payload


def normalize_pathway_id(raw_id: str) -> str:
    value = str(raw_id or "").strip()
    if not value:
        return ""
    if value.lower().startswith("pathway:"):
        value = value.split(":", 1)[1]
    return value.upper()


def default_mapping_table_path(org: str) -> str:
    return str(Path("MapKinase_WebApp") / "annotation_files" / f"{org}_mapping_table.txt")


def _normalize_col_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (name or "").strip().lower())


def _tokenize_mapping_cell(raw: str) -> List[str]:
    tokens: List[str] = []
    seen: Set[str] = set()
    for token in re.split(r"[,\s;|+]+", (raw or "").strip()):
        t = token.strip()
        if not t or t.upper() == "NA":
            continue
        if t not in seen:
            seen.add(t)
            tokens.append(t)
    return tokens


def load_species_lookup_table() -> Dict[str, str]:
    path = Path("MapKinase_WebApp") / "species_ref_list.csv"
    lookup: Dict[str, str] = {}
    if not path.exists():
        return lookup
    try:
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                cleaned = {(k.lstrip("\ufeff") if isinstance(k, str) else k): v for k, v in row.items()}
                org = (cleaned.get("Kegg Gene ID") or "").strip().lower()
                species = (cleaned.get("Species") or "").strip()
                if org and species:
                    lookup[org] = species
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("Could not read species_ref_list.csv: %s", exc)
    return lookup


def resolve_species_name(org: str, override: Optional[str]) -> str:
    if override and override.strip():
        return override.strip()
    species_lookup = load_species_lookup_table()
    if org in species_lookup:
        return species_lookup[org]
    return ORG_TO_SPECIES_DEFAULTS.get(org, org)


def parse_pathway_list(payload: Dict[str, object]) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    seen: Set[str] = set()
    rows = payload.get("pathways")
    if not isinstance(rows, list):
        return out
    for row in rows:
        if not isinstance(row, dict):
            continue
        pathway_id = normalize_pathway_id(str(row.get("id") or ""))
        if not pathway_id:
            continue
        if pathway_id in seen:
            continue
        seen.add(pathway_id)
        out.append(
            {
                "pathway_id": pathway_id,
                "name": str(row.get("name") or pathway_id).strip(),
                "species": str(row.get("species") or "").strip(),
                "revision": str(row.get("revision") or "").strip(),
            }
        )
    return out


def load_id_mapping_table(mapping_table_path: Path) -> Dict[str, object]:
    if not mapping_table_path.exists():
        LOGGER.warning("ID mapping table not found: %s", mapping_table_path)
        return {
            "uniprot_col": "",
            "columns_raw": {},
            "columns_norm": {},
            "token_to_uniprots": {},
            "row_count": 0,
            "mapped_cells": 0,
        }

    token_to_uniprots: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    columns_raw: Dict[str, str] = {}
    columns_norm: Dict[str, str] = {}
    row_count = 0
    mapped_cells = 0

    with mapping_table_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return {
                "uniprot_col": "",
                "columns_raw": {},
                "columns_norm": {},
                "token_to_uniprots": {},
                "row_count": 0,
                "mapped_cells": 0,
            }
        fieldnames = [(h or "").lstrip("\ufeff").strip() for h in reader.fieldnames]
        reader.fieldnames = fieldnames
        uniprot_col = fieldnames[0]
        for col in fieldnames:
            low = col.lower()
            norm = _normalize_col_key(col)
            if low not in columns_raw:
                columns_raw[low] = col
            if norm and norm not in columns_norm:
                columns_norm[norm] = col

        for row in reader:
            row_count += 1
            uniprot = str(row.get(uniprot_col) or "").strip()
            if not uniprot or uniprot.upper() == "NA":
                continue
            for col in fieldnames:
                if col == uniprot_col:
                    continue
                raw_val = str(row.get(col) or "").strip()
                if not raw_val or raw_val.upper() == "NA":
                    continue
                tokens = _tokenize_mapping_cell(raw_val)
                if not tokens:
                    continue
                mapped_cells += 1
                col_key = col.lower()
                token_map = token_to_uniprots[col_key]
                for tok in tokens:
                    token_map[tok].add(uniprot)
                    token_map[tok.upper()].add(uniprot)
                    token_map[tok.lower()].add(uniprot)
                    # Ensembl and other versioned IDs may include suffixes like ".12".
                    if "." in tok:
                        base_tok = tok.split(".", 1)[0]
                        if base_tok:
                            token_map[base_tok].add(uniprot)
                            token_map[base_tok.upper()].add(uniprot)
                            token_map[base_tok.lower()].add(uniprot)

    LOGGER.info(
        "Loaded ID mapping table %s (rows=%s, mapped_cells=%s, columns=%s)",
        mapping_table_path,
        row_count,
        mapped_cells,
        len(columns_raw),
    )
    # convert nested sets to sorted lists only when looked up (keep sets internally)
    return {
        "uniprot_col": uniprot_col,
        "columns_raw": columns_raw,
        "columns_norm": columns_norm,
        "token_to_uniprots": token_to_uniprots,
        "row_count": row_count,
        "mapped_cells": mapped_cells,
    }


def resolve_mapping_column(database: str, native_id: str, mapping_obj: Dict[str, object]) -> Optional[str]:
    db = (database or "").strip()
    xid = (native_id or "").strip()
    if not db:
        return None

    columns_raw = mapping_obj.get("columns_raw", {})
    columns_norm = mapping_obj.get("columns_norm", {})
    if not isinstance(columns_raw, dict) or not isinstance(columns_norm, dict):
        return None

    db_low = db.lower()
    db_norm = _normalize_col_key(db)

    if db_low == "ensembl":
        xid_upper = xid.upper()
        if xid_upper.startswith("ENST"):
            return columns_raw.get("ensembl_transcript") or columns_norm.get("ensembltranscript")
        if xid_upper.startswith("ENSP"):
            return columns_raw.get("ensembl_protein") or columns_norm.get("ensemblprotein")
        if xid_upper.startswith("ENSG"):
            return columns_raw.get("ensembl_gene") or columns_norm.get("ensemblgene")
        # fallback if no prefix match
        return (
            columns_raw.get("ensembl_gene")
            or columns_raw.get("ensembl_transcript")
            or columns_raw.get("ensembl_protein")
            or columns_norm.get("ensemblgene")
            or columns_norm.get("ensembltranscript")
            or columns_norm.get("ensemblprotein")
        )

    return columns_raw.get(db_low) or columns_norm.get(db_norm)


def map_native_id_to_uniprot(database: str, native_id: str, mapping_obj: Dict[str, object]) -> List[str]:
    db = (database or "").strip()
    xid = (native_id or "").strip()
    if not xid:
        return []

    if "uniprot" in db.lower():
        return [xid]

    column = resolve_mapping_column(db, xid, mapping_obj)
    if not column:
        return []

    token_to_uniprots = mapping_obj.get("token_to_uniprots", {})
    if not isinstance(token_to_uniprots, dict):
        return []
    col_lookup = token_to_uniprots.get(column.lower(), {})
    if not isinstance(col_lookup, dict):
        return []

    values: Set[str] = set()
    keys_to_try = [xid, xid.upper(), xid.lower()]
    if "." in xid:
        base_xid = xid.split(".", 1)[0]
        keys_to_try.extend([base_xid, base_xid.upper(), base_xid.lower()])
    for key in keys_to_try:
        matches = col_lookup.get(key)
        if isinstance(matches, set):
            values.update(matches)
    return sorted(v for v in values if v and v.upper() != "NA")


def parse_gpml_nodes(
    pathway_id: str,
    root: ET.Element,
    mapping_obj: Dict[str, object],
) -> Dict[str, Dict[str, object]]:
    local_nodes: Dict[str, Dict[str, object]] = {}
    for datanode in root.findall(".//{*}DataNode"):
        graph_id = str(datanode.attrib.get("GraphId") or "").strip()
        if not graph_id:
            continue
        label = str(datanode.attrib.get("TextLabel") or "").strip()
        node_type = str(datanode.attrib.get("Type") or "unknown").strip().lower()
        xref = datanode.find("{*}Xref")
        db_name = str(xref.attrib.get("Database") if xref is not None else "" or "").strip()
        native_id = str(xref.attrib.get("ID") if xref is not None else "" or "").strip()

        mapped_uniprot = map_native_id_to_uniprot(db_name, native_id, mapping_obj)
        node_id = f"{pathway_id}:{graph_id}"
        native_tokens: List[str] = []
        if db_name and native_id:
            native_tokens.append(f"{db_name}:{native_id}")
        elif native_id:
            native_tokens.append(native_id)

        local_nodes[node_id] = {
            "pathway_id": pathway_id,
            "wp_graph_id": graph_id,
            "type": node_type,
            "label": label,
            "xref": {
                "database": db_name,
                "id": native_id,
            },
            "candidates": {
                "kegg_genes": [],
                "gene_ids": [],
                "uniprot": mapped_uniprot,
                "symbols": [label] if label else [],
                "native_ids": native_tokens,
            },
            "degree": 0,
        }
    return local_nodes


def choose_interaction_endpoints(points: List[ET.Element]) -> Optional[Tuple[str, str, str, str]]:
    refs: List[Tuple[str, str]] = []
    for point in points:
        graph_ref = str(point.attrib.get("GraphRef") or "").strip()
        arrow_head = str(point.attrib.get("ArrowHead") or "").strip()
        if graph_ref:
            refs.append((graph_ref, arrow_head))
    if len(refs) < 2:
        return None
    first_ref, first_arrow = refs[0]
    last_ref, last_arrow = refs[-1]
    return first_ref, last_ref, first_arrow, last_arrow


def parse_gpml_edges(
    pathway_id: str,
    root: ET.Element,
    local_nodes: Dict[str, Dict[str, object]],
) -> Tuple[Dict[str, Dict[str, object]], Dict[str, Set[str]]]:
    local_edges: Dict[str, Dict[str, object]] = {}
    adjacency: Dict[str, Set[str]] = defaultdict(set)

    relations: List[ET.Element] = []
    relations.extend(root.findall(".//{*}Interaction"))
    relations.extend(root.findall(".//{*}GraphicalLine"))

    for idx, rel in enumerate(relations, start=1):
        rel_type = str(rel.attrib.get("Type") or rel.tag.split("}")[-1] or "interaction").strip()
        graphics = rel.find("{*}Graphics")
        if graphics is None:
            continue
        points = graphics.findall("{*}Point")
        if len(points) < 2:
            continue
        endpoints = choose_interaction_endpoints(points)
        if endpoints is None:
            continue
        first_ref, last_ref, first_arrow, last_arrow = endpoints

        directed = False
        src_ref = first_ref
        dst_ref = last_ref
        if last_arrow and not first_arrow:
            directed = True
            src_ref = first_ref
            dst_ref = last_ref
        elif first_arrow and not last_arrow:
            directed = True
            src_ref = last_ref
            dst_ref = first_ref

        src_node = f"{pathway_id}:{src_ref}"
        dst_node = f"{pathway_id}:{dst_ref}"
        if src_node not in local_nodes or dst_node not in local_nodes:
            continue

        subtypes: List[str] = []
        if first_arrow:
            subtypes.append(f"start:{first_arrow}")
        if last_arrow:
            subtypes.append(f"end:{last_arrow}")

        edge_id = f"{pathway_id}:{src_ref}->{dst_ref}:{idx}"
        local_edges[edge_id] = {
            "pathway_id": pathway_id,
            "src": src_node,
            "dst": dst_node,
            "directed": directed,
            "relation_type": rel_type,
            "subtypes": subtypes,
        }
        adjacency[src_node].add(dst_node)
        adjacency[dst_node].add(src_node)

    return local_edges, adjacency


def build_pairs(
    local_edges: Dict[str, Dict[str, object]],
    adjacency: Dict[str, Set[str]],
) -> Tuple[List[List[str]], List[List[object]]]:
    pairs1_set: Set[Tuple[str, str]] = set()
    for edge in local_edges.values():
        src = str(edge.get("src") or "")
        dst = str(edge.get("dst") or "")
        if not src or not dst or src == dst:
            continue
        pair = (src, dst) if src < dst else (dst, src)
        pairs1_set.add(pair)
    pairs1 = [[a, b] for a, b in sorted(pairs1_set)]

    pairs2_count: Dict[Tuple[str, str], int] = defaultdict(int)
    for _bridge, neighbors in adjacency.items():
        if len(neighbors) < 2:
            continue
        for u, v in combinations(sorted(neighbors), 2):
            pairs2_count[(u, v)] += 1
    pairs2 = [[u, v, w] for (u, v), w in sorted(pairs2_count.items())]
    return pairs1, pairs2


def parse_gpml_pathway(
    pathway_id: str,
    pathway_name: str,
    gpml_text: str,
    mapping_obj: Dict[str, object],
) -> Tuple[Dict[str, object], Dict[str, Dict[str, object]], Dict[str, Dict[str, object]]]:
    root = ET.fromstring(gpml_text)
    display_name = str(root.attrib.get("Name") or pathway_name or pathway_id).strip()
    organism = str(root.attrib.get("Organism") or "").strip()

    local_nodes = parse_gpml_nodes(pathway_id=pathway_id, root=root, mapping_obj=mapping_obj)
    local_edges, adjacency = parse_gpml_edges(pathway_id=pathway_id, root=root, local_nodes=local_nodes)

    for node_id, node in local_nodes.items():
        node["degree"] = len(adjacency.get(node_id, set()))

    pairs1, pairs2 = build_pairs(local_edges=local_edges, adjacency=adjacency)

    pathway_obj: Dict[str, object] = {
        "pathway_id": pathway_id,
        "name": display_name,
        "species": organism,
        "nodes": sorted(local_nodes.keys()),
        "edges": sorted(local_edges.keys()),
        "pairs1": pairs1,
        "pairs2": pairs2,
        "node_count": len(local_nodes),
        "edge_count": len(local_edges),
    }
    return pathway_obj, local_nodes, local_edges


def load_or_parse_pathway(
    pathway_id: str,
    pathway_name: str,
    gpml_text: str,
    parsed_cache_path: Path,
    mapping_obj: Dict[str, object],
) -> Tuple[Dict[str, object], Dict[str, Dict[str, object]], Dict[str, Dict[str, object]]]:
    gpml_hash = sha256_text(gpml_text)
    if parsed_cache_path.exists() and parsed_cache_path.stat().st_size > 0:
        try:
            cached = json.loads(read_text(parsed_cache_path))
            if (
                cached.get("parser_version") == PARSER_VERSION
                and cached.get("gpml_sha256") == gpml_hash
                and "pathway" in cached
                and "nodes" in cached
                and "edges" in cached
            ):
                pathway = cached["pathway"]
                if not pathway.get("name"):
                    pathway["name"] = pathway_name
                return pathway, cached["nodes"], cached["edges"]
        except Exception as exc:  # noqa: BLE001
            LOGGER.warning("Ignoring unreadable parsed cache for %s: %s", pathway_id, exc)

    pathway, nodes, edges = parse_gpml_pathway(
        pathway_id=pathway_id,
        pathway_name=pathway_name,
        gpml_text=gpml_text,
        mapping_obj=mapping_obj,
    )
    write_json_atomic(
        parsed_cache_path,
        {
            "parser_version": PARSER_VERSION,
            "gpml_sha256": gpml_hash,
            "pathway": pathway,
            "nodes": nodes,
            "edges": edges,
        },
        pretty=False,
    )
    return pathway, nodes, edges


def validate_index(
    pathways: Iterable[Dict[str, object]],
    nodes: Dict[str, Dict[str, object]],
    edges: Dict[str, Dict[str, object]],
) -> List[str]:
    errors: List[str] = []
    for edge_id, edge in edges.items():
        src = edge.get("src")
        dst = edge.get("dst")
        if src not in nodes:
            errors.append(f"Edge {edge_id} references missing src node {src}")
        if dst not in nodes:
            errors.append(f"Edge {edge_id} references missing dst node {dst}")
    for pathway in pathways:
        pathway_id = pathway.get("pathway_id", "<unknown>")
        node_ids = pathway.get("nodes", [])
        edge_ids = pathway.get("edges", [])
        if pathway.get("node_count") != len(node_ids):
            errors.append(f"Pathway {pathway_id}: node_count does not match nodes length")
        if pathway.get("edge_count") != len(edge_ids):
            errors.append(f"Pathway {pathway_id}: edge_count does not match edges length")
    return errors


def compute_stats(
    pathways: List[Dict[str, object]],
    nodes: Dict[str, Dict[str, object]],
    edges: Dict[str, Dict[str, object]],
) -> Dict[str, object]:
    uni_counts: List[int] = []
    mapped_nodes = 0
    for node in nodes.values():
        candidates = node.get("candidates", {})
        unis = candidates.get("uniprot", []) if isinstance(candidates, dict) else []
        count = len(unis) if isinstance(unis, list) else 0
        uni_counts.append(count)
        if count > 0:
            mapped_nodes += 1
    avg_uni = (sum(uni_counts) / len(uni_counts)) if uni_counts else 0.0
    return {
        "pathway_count": len(pathways),
        "node_count": len(nodes),
        "edge_count": len(edges),
        "nodes_with_uniprot": mapped_nodes,
        "avg_uniprot_candidates_per_node": round(avg_uni, 4),
    }


def extract_gpml_from_get_pathway_payload(payload: Dict[str, object], pathway_id: str) -> str:
    pathway = payload.get("pathway")
    if not isinstance(pathway, dict):
        raise WikiPathwaysNotFoundError(f"{pathway_id}: invalid getPathway payload (missing pathway object)")
    gpml = pathway.get("gpml")
    if not isinstance(gpml, str) or not gpml.strip():
        raise WikiPathwaysNotFoundError(f"{pathway_id}: missing gpml text in payload")
    return gpml


def fetch_gpml_text(
    session: requests.Session,
    pathway_id: str,
    gpml_cache_path: Path,
    rate_limiter: RateLimiter,
) -> str:
    if gpml_cache_path.exists() and gpml_cache_path.stat().st_size > 0:
        return read_text(gpml_cache_path)

    url = f"{WIKIPATHWAYS_API_BASE}/getPathway?pwId={quote(pathway_id)}&format=json"
    payload = fetch_json(session=session, url=url, rate_limiter=rate_limiter)
    gpml = extract_gpml_from_get_pathway_payload(payload=payload, pathway_id=pathway_id)
    write_text_atomic(gpml_cache_path, gpml)
    return gpml


def main(argv: Optional[List[str]] = None) -> int:
    try:
        args = get_runtime_args(argv)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    org = str(args.org or "").strip().lower()
    out_path = Path(args.out)
    cache_dir = Path(args.cache)
    species_name = resolve_species_name(org=org, override=args.species_name)
    mapping_table_path = Path(args.id_mapping_table) if args.id_mapping_table else Path(default_mapping_table_path(org))

    list_cache_dir = cache_dir / "list"
    gpml_cache_dir = cache_dir / "gpml" / org
    parsed_cache_dir = cache_dir / "parsed" / org
    ensure_dir(list_cache_dir)
    ensure_dir(gpml_cache_dir)
    ensure_dir(parsed_cache_dir)

    rate_limiter = RateLimiter(args.rate_limit)
    session = requests.Session()
    session.headers.update({"User-Agent": "build_wikipathways_index.py/1.0"})

    list_cache_path = list_cache_dir / f"pathways_{org}.json"
    list_url = f"{WIKIPATHWAYS_API_BASE}/listPathways?organism={quote(species_name)}&format=json"
    LOGGER.info("Loading WikiPathways list for org=%s species=%s", org, species_name)
    try:
        list_payload = load_or_fetch_json(
            session=session,
            url=list_url,
            cache_path=list_cache_path,
            rate_limiter=rate_limiter,
        )
    except Exception as exc:  # noqa: BLE001
        LOGGER.error("Failed to load pathway list for %s: %s", org, exc)
        return 2

    pathway_list = parse_pathway_list(list_payload)
    if not pathway_list:
        LOGGER.error("No WikiPathways found for species=%s", species_name)
        return 2

    if args.max_pathways is not None:
        if args.max_pathways < 0:
            LOGGER.error("--max-pathways must be >= 0")
            return 2
        pathway_list = pathway_list[: args.max_pathways]

    mapping_obj = load_id_mapping_table(mapping_table_path)

    total = len(pathway_list)
    LOGGER.info("Pathways to process: %s", total)

    all_pathways: List[Dict[str, object]] = []
    all_nodes: Dict[str, Dict[str, object]] = {}
    all_edges: Dict[str, Dict[str, object]] = {}
    failures: List[str] = []

    for idx, meta in enumerate(pathway_list, start=1):
        pathway_id = meta["pathway_id"]
        pathway_name = meta.get("name", pathway_id)
        LOGGER.info("[%s/%s] Processing %s", idx, total, pathway_id)
        gpml_cache_path = gpml_cache_dir / f"{pathway_id}.gpml"
        parsed_cache_path = parsed_cache_dir / f"{pathway_id}.parsed.json"

        try:
            gpml_text = fetch_gpml_text(
                session=session,
                pathway_id=pathway_id,
                gpml_cache_path=gpml_cache_path,
                rate_limiter=rate_limiter,
            )
        except WikiPathwaysNotFoundError as exc:
            msg = f"{pathway_id}: no GPML payload: {exc}"
            LOGGER.warning(msg)
            failures.append(msg)
            continue
        except Exception as exc:  # noqa: BLE001
            msg = f"{pathway_id}: failed to download GPML: {exc}"
            LOGGER.warning(msg)
            failures.append(msg)
            continue

        try:
            pathway_obj, nodes_obj, edges_obj = load_or_parse_pathway(
                pathway_id=pathway_id,
                pathway_name=pathway_name,
                gpml_text=gpml_text,
                parsed_cache_path=parsed_cache_path,
                mapping_obj=mapping_obj,
            )
        except ET.ParseError as exc:
            msg = f"{pathway_id}: XML parse error: {exc}"
            LOGGER.warning(msg)
            failures.append(msg)
            continue
        except Exception as exc:  # noqa: BLE001
            msg = f"{pathway_id}: parse failure: {exc}"
            LOGGER.warning(msg)
            failures.append(msg)
            continue

        overlap_nodes = set(all_nodes).intersection(nodes_obj)
        overlap_edges = set(all_edges).intersection(edges_obj)
        if overlap_nodes:
            msg = f"{pathway_id}: duplicate node IDs encountered ({len(overlap_nodes)})"
            LOGGER.warning(msg)
            failures.append(msg)
            continue
        if overlap_edges:
            msg = f"{pathway_id}: duplicate edge IDs encountered ({len(overlap_edges)})"
            LOGGER.warning(msg)
            failures.append(msg)
            continue

        all_pathways.append(pathway_obj)
        all_nodes.update(nodes_obj)
        all_edges.update(edges_obj)

    if not all_pathways:
        LOGGER.error("No pathways were successfully processed.")
        return 2

    validation_errors = validate_index(all_pathways, all_nodes, all_edges)
    if validation_errors:
        LOGGER.error("Validation failed with %s issue(s).", len(validation_errors))
        for err in validation_errors[:20]:
            LOGGER.error("  %s", err)
        if len(validation_errors) > 20:
            LOGGER.error("  ... and %s more", len(validation_errors) - 20)
        return 2

    stats = compute_stats(all_pathways, all_nodes, all_edges)
    notes = [
        "Index contains WikiPathways topology and node-to-UniProt candidates only.",
        "Node IDs are stable as '{pathway_id}:{GraphId}'.",
        "Each DataNode is treated as one protbox-level node.",
        "Node native IDs are mapped via Xref Database + ID using organism mapping table.",
        "For Database='Ensembl', ENST/ENSP/ENSG prefixes route to transcript/protein/gene mapping columns.",
        "pairs1 stores unique undirected 1-hop node pairs.",
        "pairs2 stores unique undirected 2-hop pairs with bridge_count.",
    ]
    if failures:
        notes.append(f"{len(failures)} pathways failed or were skipped. See failures list.")

    output = {
        "meta": {
            "schema_version": SCHEMA_VERSION,
            "pathway_source": "wikipathways",
            "org": org,
            "species": species_name,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "wikipathways_api_base": WIKIPATHWAYS_API_BASE,
            "id_mapping_table": str(mapping_table_path),
            "notes": notes,
            "stats": stats,
            "failures": failures,
        },
        "pathways": all_pathways,
        "nodes": all_nodes,
        "edges": all_edges,
    }

    LOGGER.info(
        "Writing index to %s (pathways=%s, nodes=%s, edges=%s)",
        out_path,
        len(all_pathways),
        len(all_nodes),
        len(all_edges),
    )
    write_json_atomic(out_path, output, pretty=args.pretty)
    LOGGER.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
