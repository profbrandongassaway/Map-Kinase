#!/usr/bin/env python3
"""
Build a KEGG pathway index for an organism using KEGG REST + KGML.

This script downloads pathway KGML files, parses pathway topology and node
candidate mappings, precomputes 1-hop and 2-hop node pairs, and writes a
dataset-independent index JSON for fast downstream scoring.
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
import xml.etree.ElementTree as ET

import requests


SCHEMA_VERSION = 1
PARSER_VERSION = 1
KEGG_API_BASE = "https://rest.kegg.jp"
DEFAULT_RATE_LIMIT = 0.25

LOGGER = logging.getLogger("build_kegg_index")

# IDE-friendly config:
# Edit these values and run this file directly (no terminal args required).
CONFIG_ORG = "dme"
CONFIG_OUT = f"MapKinase_WebApp/index_files/kegg_index_{CONFIG_ORG}.json"
CONFIG_CACHE = ".kegg_cache"
CONFIG_INCLUDE_CLASSES = False
CONFIG_MAX_PATHWAYS = None
CONFIG_RATE_LIMIT = DEFAULT_RATE_LIMIT
CONFIG_PRETTY = False
CONFIG_LOG_LEVEL = "INFO"
# If None, default is "MapKinase_WebApp/annotation_files/{org}_id_mapping_table.txt".
CONFIG_ID_MAPPING_TABLE = None


class KeggFetchError(RuntimeError):
    """Raised when KEGG content cannot be fetched after retries."""


class KeggNotFoundError(KeggFetchError):
    """Raised when KEGG returns 404 (resource missing)."""


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
        description="Build a cached KEGG pathway index for an organism."
    )
    parser.add_argument("--org", required=False, help='KEGG organism code (e.g. "hsa").')
    parser.add_argument("--out", required=False, help="Output JSON file path.")
    parser.add_argument(
        "--cache",
        default=".kegg_cache",
        help='Cache directory (default: ".kegg_cache").',
    )
    parser.add_argument(
        "--include-classes",
        action="store_true",
        help="Include pathway class metadata if found in KGML.",
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
        help=f"Seconds between KEGG requests (default: {DEFAULT_RATE_LIMIT}).",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print output JSON.",
    )
    parser.add_argument(
        "--id-mapping-table",
        default=None,
        help=(
            "Tab-delimited ID mapping table path (UniProt in first column, "
            "KEGG IDs in second column). Default: "
            "MapKinase_WebApp/annotation_files/{org}_mapping_table.txt"
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser.parse_args(argv)


def get_runtime_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Runtime behavior:
    - If no args are provided, use CONFIG_* values (IDE run mode).
    - If args are provided, parse CLI args and require --org/--out.
    """
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) == 0:
        if not CONFIG_ORG or not CONFIG_OUT:
            raise ValueError("CONFIG_ORG and CONFIG_OUT must be set for IDE mode.")
        return argparse.Namespace(
            org=CONFIG_ORG,
            out=CONFIG_OUT,
            cache=CONFIG_CACHE,
            include_classes=CONFIG_INCLUDE_CLASSES,
            max_pathways=CONFIG_MAX_PATHWAYS,
            rate_limit=CONFIG_RATE_LIMIT,
            pretty=CONFIG_PRETTY,
            id_mapping_table=CONFIG_ID_MAPPING_TABLE,
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


def fetch_text(
    session: requests.Session,
    url: str,
    rate_limiter: RateLimiter,
    timeout: int = 30,
    max_retries: int = 5,
) -> str:
    last_error: Optional[Exception] = None
    for attempt in range(1, max_retries + 1):
        try:
            rate_limiter.wait()
            response = session.get(url, timeout=timeout)
            if response.status_code == 404:
                raise KeggNotFoundError(f"404 for URL: {url}")
            response.raise_for_status()
            return response.text
        except KeggNotFoundError:
            raise
        except requests.RequestException as exc:
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
    raise KeggFetchError(f"Failed to fetch {url}") from last_error


def load_or_fetch_text(
    session: requests.Session,
    url: str,
    cache_path: Path,
    rate_limiter: RateLimiter,
) -> str:
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return read_text(cache_path)
    text = fetch_text(session=session, url=url, rate_limiter=rate_limiter)
    write_text_atomic(cache_path, text)
    return text


def normalize_pathway_id(raw_id: str) -> str:
    value = raw_id.strip()
    if value.startswith("path:"):
        value = value[5:]
    if ":" in value:
        value = value.split(":", 1)[1]
    return value


def default_mapping_table_path(org: str) -> str:
    return str(Path("MapKinase_WebApp") / "annotation_files" / f"{org}_mapping_table.txt")


def parse_mapping_kegg_tokens(raw_value: str, org: str) -> List[str]:
    tokens: List[str] = []
    seen: Set[str] = set()
    if not raw_value:
        return tokens
    for part in re.split(r"[,\s;|+]+", raw_value.strip()):
        part = part.strip()
        if not part or part.upper() == "NA":
            continue
        if ":" not in part and part.isdigit():
            part = f"{org}:{part}"
        if ":" in part:
            prefix, gene = part.split(":", 1)
            prefix = prefix.strip()
            gene = gene.strip()
            if not prefix or not gene:
                continue
            part = f"{prefix}:{gene}"
        if part not in seen:
            seen.add(part)
            tokens.append(part)
    return tokens


def load_kegg_to_uniprot_map(mapping_table_path: Path, org: str) -> Dict[str, List[str]]:
    if not mapping_table_path.exists():
        LOGGER.warning("Mapping table not found: %s", mapping_table_path)
        return {}

    mapping: Dict[str, List[str]] = defaultdict(list)
    with mapping_table_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            LOGGER.warning("Mapping table has no header/data: %s", mapping_table_path)
            return {}
        if len(header) < 2:
            LOGGER.warning(
                "Mapping table must have at least 2 columns (UniProt in col1, KEGG in col2): %s",
                mapping_table_path,
            )
            return {}

        row_count = 0
        mapped_rows = 0
        for row in reader:
            row_count += 1
            uniprot_id = (row[0] if len(row) > 0 else "").strip()
            kegg_raw = (row[1] if len(row) > 1 else "").strip()
            if not uniprot_id or uniprot_id.upper() == "NA":
                continue
            kegg_tokens = parse_mapping_kegg_tokens(kegg_raw, org=org)
            if not kegg_tokens:
                continue
            mapped_rows += 1
            for kegg_gene in kegg_tokens:
                existing = mapping[kegg_gene]
                if uniprot_id not in existing:
                    existing.append(uniprot_id)

    LOGGER.info(
        "Loaded KEGG->UniProt map from %s (rows=%s, mapped_rows=%s, unique_kegg=%s)",
        mapping_table_path,
        row_count,
        mapped_rows,
        len(mapping),
    )
    return dict(mapping)


def apply_uniprot_mapping_to_nodes(
    nodes: Dict[str, Dict[str, object]],
    kegg_to_uniprot: Dict[str, List[str]],
) -> Dict[str, int]:
    nodes_with_uniprot = 0
    total_uniprot_links = 0
    for node in nodes.values():
        candidates = node.get("candidates")
        if not isinstance(candidates, dict):
            continue
        kegg_genes = candidates.get("kegg_genes", [])
        if not isinstance(kegg_genes, list):
            kegg_genes = []

        seen: Set[str] = set()
        mapped_uniprots: List[str] = []
        for kegg_gene in kegg_genes:
            if not isinstance(kegg_gene, str):
                continue
            for uniprot_id in kegg_to_uniprot.get(kegg_gene, []):
                if uniprot_id not in seen:
                    seen.add(uniprot_id)
                    mapped_uniprots.append(uniprot_id)

        candidates["uniprot"] = mapped_uniprots
        if mapped_uniprots:
            nodes_with_uniprot += 1
            total_uniprot_links += len(mapped_uniprots)

    return {
        "nodes_with_uniprot": nodes_with_uniprot,
        "total_uniprot_links": total_uniprot_links,
    }


def parse_pathway_list(text: str, org: str) -> List[Dict[str, str]]:
    pathways: List[Dict[str, str]] = []
    seen: Set[str] = set()
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            LOGGER.warning("Skipping malformed pathway list line: %r", line)
            continue
        raw_id, name = parts
        pathway_id = normalize_pathway_id(raw_id)
        if not pathway_id:
            continue
        if not pathway_id.startswith(org):
            LOGGER.debug("Skipping pathway outside org %s: %s", org, pathway_id)
            continue
        if pathway_id in seen:
            continue
        seen.add(pathway_id)
        pathways.append({"pathway_id": pathway_id, "name": name.strip()})
    return pathways


def sanitize_token(token: str) -> str:
    return token.strip().strip(",;")


def parse_gene_candidates(name_attr: str) -> Tuple[List[str], List[str]]:
    """
    Parse KEGG gene identifiers from an entry name string.

    Common patterns:
    - "hsa:207 208 209"
    - "hsa:207 hsa:208"
    - "hsa:5594+hsa:5595"
    """
    kegg_genes: List[str] = []
    gene_ids: List[str] = []
    seen_full: Set[str] = set()
    seen_gene_id: Set[str] = set()

    current_prefix: Optional[str] = None
    raw_tokens = re.split(r"\s+", (name_attr or "").strip())
    for raw in raw_tokens:
        if not raw:
            continue
        subtokens = re.split(r"[+]", raw)
        for sub in subtokens:
            token = sanitize_token(sub)
            if not token:
                continue
            if ":" in token:
                prefix, gene = token.split(":", 1)
                prefix = sanitize_token(prefix)
                gene = sanitize_token(gene)
                if not prefix or not gene:
                    continue
                current_prefix = prefix
                full = f"{prefix}:{gene}"
            else:
                if not current_prefix:
                    continue
                gene = sanitize_token(token)
                if not gene:
                    continue
                full = f"{current_prefix}:{gene}"

            if full not in seen_full:
                seen_full.add(full)
                kegg_genes.append(full)
            gene_id = full.split(":", 1)[1]
            if gene_id and gene_id not in seen_gene_id:
                seen_gene_id.add(gene_id)
                gene_ids.append(gene_id)

    return kegg_genes, gene_ids


def resolve_group_candidates(
    entry_id: int,
    entries: Dict[int, Dict[str, object]],
    memo: Dict[int, Tuple[List[str], List[str]]],
    stack: Optional[Set[int]] = None,
) -> Tuple[List[str], List[str]]:
    if entry_id in memo:
        return memo[entry_id]
    if stack is None:
        stack = set()
    if entry_id in stack:
        return ([], [])
    stack.add(entry_id)

    entry = entries.get(entry_id)
    if not entry:
        stack.remove(entry_id)
        memo[entry_id] = ([], [])
        return memo[entry_id]

    if entry.get("type") != "group":
        kegg_genes = list(entry.get("kegg_genes", []))
        gene_ids = list(entry.get("gene_ids", []))
        stack.remove(entry_id)
        memo[entry_id] = (kegg_genes, gene_ids)
        return memo[entry_id]

    seen_kegg: Set[str] = set()
    seen_gene: Set[str] = set()
    merged_kegg: List[str] = []
    merged_gene: List[str] = []
    component_ids = entry.get("components", [])
    if isinstance(component_ids, list):
        for comp_id in component_ids:
            if not isinstance(comp_id, int):
                continue
            comp_kegg, comp_gene = resolve_group_candidates(comp_id, entries, memo, stack)
            for kegg_gene in comp_kegg:
                if kegg_gene not in seen_kegg:
                    seen_kegg.add(kegg_gene)
                    merged_kegg.append(kegg_gene)
            for gene_id in comp_gene:
                if gene_id not in seen_gene:
                    seen_gene.add(gene_id)
                    merged_gene.append(gene_id)

    stack.remove(entry_id)
    memo[entry_id] = (merged_kegg, merged_gene)
    return memo[entry_id]


def parse_kgml(
    pathway_id: str,
    pathway_name: str,
    kgml_text: str,
    include_classes: bool = False,
) -> Tuple[Dict[str, object], Dict[str, Dict[str, object]], Dict[str, Dict[str, object]]]:
    root = ET.fromstring(kgml_text)

    title = root.attrib.get("title") or pathway_name or pathway_id
    class_attr = root.attrib.get("class", "")
    classes = [c.strip() for c in class_attr.split(";") if c.strip()] if class_attr else []

    entry_map: Dict[int, Dict[str, object]] = {}
    entry_order: List[int] = []
    for entry in root.findall("entry"):
        raw_entry_id = entry.attrib.get("id", "").strip()
        if not raw_entry_id:
            continue
        try:
            entry_id = int(raw_entry_id)
        except ValueError:
            LOGGER.warning("Skipping non-integer entry id in %s: %r", pathway_id, raw_entry_id)
            continue
        entry_type = entry.attrib.get("type", "").strip() or "unknown"
        name_attr = entry.attrib.get("name", "").strip()
        graphics = entry.find("graphics")
        label = ""
        if graphics is not None:
            label = (graphics.attrib.get("name") or "").strip()
        if not label:
            label = name_attr

        component_ids: List[int] = []
        for comp in entry.findall("component"):
            comp_id_raw = comp.attrib.get("id", "").strip()
            if not comp_id_raw:
                continue
            try:
                component_ids.append(int(comp_id_raw))
            except ValueError:
                continue

        kegg_genes: List[str] = []
        gene_ids: List[str] = []
        if entry_type == "gene":
            kegg_genes, gene_ids = parse_gene_candidates(name_attr)

        entry_map[entry_id] = {
            "entry_id": entry_id,
            "type": entry_type,
            "name_attr": name_attr,
            "label": label,
            "components": component_ids,
            "kegg_genes": kegg_genes,
            "gene_ids": gene_ids,
        }
        entry_order.append(entry_id)

    # Resolve group entries by unioning candidate genes from component entries.
    memo: Dict[int, Tuple[List[str], List[str]]] = {}
    for entry_id, info in entry_map.items():
        if info.get("type") == "group":
            kegg_genes, gene_ids = resolve_group_candidates(entry_id, entry_map, memo)
            info["kegg_genes"] = kegg_genes
            info["gene_ids"] = gene_ids

    local_nodes: Dict[str, Dict[str, object]] = {}
    for entry_id in entry_order:
        info = entry_map[entry_id]
        node_id = f"{pathway_id}:{entry_id}"
        local_nodes[node_id] = {
            "pathway_id": pathway_id,
            "kegg_entry_id": entry_id,
            "type": info.get("type", "unknown"),
            "label": info.get("label", ""),
            "candidates": {
                "kegg_genes": list(info.get("kegg_genes", [])),
                "gene_ids": list(info.get("gene_ids", [])),
                "uniprot": [],
                "symbols": [],
            },
            "degree": 0,
        }

    local_edges: Dict[str, Dict[str, object]] = {}
    adjacency: Dict[str, Set[str]] = defaultdict(set)
    for relation_idx, relation in enumerate(root.findall("relation"), start=1):
        entry1_raw = relation.attrib.get("entry1", "").strip()
        entry2_raw = relation.attrib.get("entry2", "").strip()
        rel_type = relation.attrib.get("type", "").strip() or "unknown"
        if not entry1_raw or not entry2_raw:
            continue
        try:
            entry1 = int(entry1_raw)
            entry2 = int(entry2_raw)
        except ValueError:
            continue

        src_node = f"{pathway_id}:{entry1}"
        dst_node = f"{pathway_id}:{entry2}"
        if src_node not in local_nodes or dst_node not in local_nodes:
            LOGGER.debug(
                "Skipping relation referencing missing nodes in %s: %s -> %s",
                pathway_id,
                src_node,
                dst_node,
            )
            continue

        subtypes: List[str] = []
        for subtype in relation.findall("subtype"):
            subtype_name = (subtype.attrib.get("name") or "").strip()
            subtype_value = (subtype.attrib.get("value") or "").strip()
            if subtype_name:
                subtypes.append(subtype_name)
            elif subtype_value:
                subtypes.append(subtype_value)

        edge_id = f"{pathway_id}:{entry1}->{entry2}:{relation_idx}"
        local_edges[edge_id] = {
            "pathway_id": pathway_id,
            "src": src_node,
            "dst": dst_node,
            "directed": True,
            "relation_type": rel_type,
            "subtypes": subtypes,
        }
        adjacency[src_node].add(dst_node)
        adjacency[dst_node].add(src_node)

    for node_id, node in local_nodes.items():
        node["degree"] = len(adjacency.get(node_id, set()))

    pairs1_set: Set[Tuple[str, str]] = set()
    for edge in local_edges.values():
        src = str(edge["src"])
        dst = str(edge["dst"])
        if src == dst:
            continue
        pair = (src, dst) if src < dst else (dst, src)
        pairs1_set.add(pair)
    pairs1 = [[a, b] for a, b in sorted(pairs1_set)]

    pair2_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    for _bridge_node, neighbors in adjacency.items():
        if len(neighbors) < 2:
            continue
        sorted_neighbors = sorted(neighbors)
        for u, v in combinations(sorted_neighbors, 2):
            pair2_counts[(u, v)] += 1
    pairs2 = [[u, v, w] for (u, v), w in sorted(pair2_counts.items())]

    node_ids_sorted = sorted(
        local_nodes.keys(),
        key=lambda node_key: int(node_key.rsplit(":", 1)[1]),
    )
    edge_ids_sorted = sorted(local_edges.keys())

    pathway_obj: Dict[str, object] = {
        "pathway_id": pathway_id,
        "name": title,
        "nodes": node_ids_sorted,
        "edges": edge_ids_sorted,
        "pairs1": pairs1,
        "pairs2": pairs2,
        "node_count": len(local_nodes),
        "edge_count": len(local_edges),
    }
    if include_classes:
        pathway_obj["classes"] = classes

    return pathway_obj, local_nodes, local_edges


def load_or_parse_pathway(
    pathway_id: str,
    pathway_name: str,
    kgml_text: str,
    parsed_cache_path: Path,
    include_classes: bool,
) -> Tuple[Dict[str, object], Dict[str, Dict[str, object]], Dict[str, Dict[str, object]]]:
    kgml_hash = sha256_text(kgml_text)
    if parsed_cache_path.exists() and parsed_cache_path.stat().st_size > 0:
        try:
            cached = json.loads(read_text(parsed_cache_path))
            if (
                cached.get("parser_version") == PARSER_VERSION
                and cached.get("kgml_sha256") == kgml_hash
                and "pathway" in cached
                and "nodes" in cached
                and "edges" in cached
            ):
                pathway = cached["pathway"]
                if not pathway.get("name"):
                    pathway["name"] = pathway_name
                if include_classes and "classes" not in pathway:
                    pathway["classes"] = []
                if not include_classes and "classes" in pathway:
                    pathway.pop("classes", None)
                return pathway, cached["nodes"], cached["edges"]
        except Exception as exc:  # noqa: BLE001
            LOGGER.warning("Ignoring unreadable parsed cache for %s: %s", pathway_id, exc)

    pathway, nodes, edges = parse_kgml(
        pathway_id=pathway_id,
        pathway_name=pathway_name,
        kgml_text=kgml_text,
        include_classes=include_classes,
    )
    cached_obj = {
        "parser_version": PARSER_VERSION,
        "kgml_sha256": kgml_hash,
        "pathway": pathway,
        "nodes": nodes,
        "edges": edges,
    }
    write_json_atomic(parsed_cache_path, cached_obj, pretty=False)
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
        pathway_nodes = pathway.get("nodes", [])
        pathway_edges = pathway.get("edges", [])
        if pathway.get("node_count") != len(pathway_nodes):
            errors.append(f"Pathway {pathway_id}: node_count does not match nodes length")
        if pathway.get("edge_count") != len(pathway_edges):
            errors.append(f"Pathway {pathway_id}: edge_count does not match edges length")
        for node_id in pathway_nodes:
            if node_id not in nodes:
                errors.append(f"Pathway {pathway_id}: missing node in global map: {node_id}")
        for edge_id in pathway_edges:
            if edge_id not in edges:
                errors.append(f"Pathway {pathway_id}: missing edge in global map: {edge_id}")

    return errors


def compute_stats(
    pathways: List[Dict[str, object]],
    nodes: Dict[str, Dict[str, object]],
    edges: Dict[str, Dict[str, object]],
) -> Dict[str, object]:
    candidate_counts: List[int] = []
    for node in nodes.values():
        candidates = node.get("candidates", {})
        kegg_genes = candidates.get("kegg_genes", []) if isinstance(candidates, dict) else []
        candidate_counts.append(len(kegg_genes))

    avg_candidates = (sum(candidate_counts) / len(candidate_counts)) if candidate_counts else 0.0
    max_candidates = max(candidate_counts) if candidate_counts else 0

    return {
        "pathway_count": len(pathways),
        "node_count": len(nodes),
        "edge_count": len(edges),
        "avg_kegg_candidates_per_node": round(avg_candidates, 4),
        "max_kegg_candidates_per_node": max_candidates,
    }


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

    org = args.org.strip()
    out_path = Path(args.out)
    cache_dir = Path(args.cache)
    mapping_table_arg = args.id_mapping_table
    mapping_table_path = Path(mapping_table_arg) if mapping_table_arg else Path(
        default_mapping_table_path(org)
    )

    list_cache_dir = cache_dir / "list"
    kgml_cache_dir = cache_dir / "kgml" / org
    parsed_cache_dir = cache_dir / "parsed" / org
    ensure_dir(list_cache_dir)
    ensure_dir(kgml_cache_dir)
    ensure_dir(parsed_cache_dir)

    rate_limiter = RateLimiter(args.rate_limit)
    session = requests.Session()
    session.headers.update({"User-Agent": "build_kegg_index.py/1.0"})

    list_cache_path = list_cache_dir / f"pathway_{org}.txt"
    list_url = f"{KEGG_API_BASE}/list/pathway/{org}"
    LOGGER.info("Loading pathway list for org=%s", org)
    try:
        pathway_list_text = load_or_fetch_text(
            session=session,
            url=list_url,
            cache_path=list_cache_path,
            rate_limiter=rate_limiter,
        )
    except Exception as exc:  # noqa: BLE001
        LOGGER.error("Failed to load pathway list for %s: %s", org, exc)
        return 2

    pathways_meta = parse_pathway_list(pathway_list_text, org=org)
    if not pathways_meta:
        LOGGER.error("No pathways found for org=%s", org)
        return 2

    if args.max_pathways is not None:
        if args.max_pathways < 0:
            LOGGER.error("--max-pathways must be >= 0")
            return 2
        pathways_meta = pathways_meta[: args.max_pathways]

    total = len(pathways_meta)
    LOGGER.info("Pathways to process: %s", total)

    all_pathways: List[Dict[str, object]] = []
    all_nodes: Dict[str, Dict[str, object]] = {}
    all_edges: Dict[str, Dict[str, object]] = {}
    failures: List[str] = []

    for idx, pathway_meta in enumerate(pathways_meta, start=1):
        pathway_id = pathway_meta["pathway_id"]
        pathway_name = pathway_meta.get("name", pathway_id)
        LOGGER.info("[%s/%s] Processing %s", idx, total, pathway_id)
        kgml_url = f"{KEGG_API_BASE}/get/{pathway_id}/kgml"
        kgml_cache_path = kgml_cache_dir / f"{pathway_id}.kgml.xml"
        parsed_cache_path = parsed_cache_dir / f"{pathway_id}.parsed.json"

        try:
            kgml_text = load_or_fetch_text(
                session=session,
                url=kgml_url,
                cache_path=kgml_cache_path,
                rate_limiter=rate_limiter,
            )
        except KeggNotFoundError:
            msg = f"{pathway_id}: KGML not found (404)"
            LOGGER.warning(msg)
            failures.append(msg)
            continue
        except Exception as exc:  # noqa: BLE001
            msg = f"{pathway_id}: failed to download KGML: {exc}"
            LOGGER.warning(msg)
            failures.append(msg)
            continue

        try:
            pathway_obj, nodes_obj, edges_obj = load_or_parse_pathway(
                pathway_id=pathway_id,
                pathway_name=pathway_name,
                kgml_text=kgml_text,
                parsed_cache_path=parsed_cache_path,
                include_classes=args.include_classes,
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

    kegg_to_uniprot = load_kegg_to_uniprot_map(mapping_table_path=mapping_table_path, org=org)
    uniprot_stats = apply_uniprot_mapping_to_nodes(all_nodes, kegg_to_uniprot)

    stats = compute_stats(all_pathways, all_nodes, all_edges)
    stats["nodes_with_uniprot"] = uniprot_stats["nodes_with_uniprot"]
    stats["total_uniprot_links"] = uniprot_stats["total_uniprot_links"]
    stats["kegg_to_uniprot_keys"] = len(kegg_to_uniprot)
    notes = [
        "Index contains pathway topology and node candidate mappings only; no user data.",
        "Node IDs are stable as '{pathway_id}:{kegg_entry_id}'.",
        "pairs1 stores unique undirected 1-hop node pairs.",
        "pairs2 stores unique undirected 2-hop candidate pairs with bridge_count.",
        "Edges are stored as directed from KGML relation entry1 -> entry2.",
        f"Node candidate KEGG IDs were mapped to UniProt using: {mapping_table_path}",
    ]
    if failures:
        notes.append(f"{len(failures)} pathways failed or were skipped. See failures list.")

    output = {
        "meta": {
            "schema_version": SCHEMA_VERSION,
            "org": org,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "kegg_api_base": KEGG_API_BASE,
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
