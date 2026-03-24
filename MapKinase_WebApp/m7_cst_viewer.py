from __future__ import annotations

import base64
import html
import json
import math
import re
import zlib
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from shiny import ui

from MapKinase_WebApp.m8_pathway_label_mapper import PathwayLabelMapper
from MapKinase_WebApp.m11_cst_pathway_index import get_cst_pathway_mapping


_DEFAULT_CST_FILES = [
    {
        "id": "cst_p38_mapk_activation",
        "name": "p38 MAPK Activation",
        "filename": "p38 MAPK Activation (1).ai",
    },
    {
        "id": "cst_mtor_signaling",
        "name": "mTOR Signaling",
        "filename": "mTOR Signaling.ai",
    }
]

_CST_PAGE_WIDTH = 612.0
_CST_PAGE_HEIGHT = 699.627
_DEFAULT_NEGATIVE_COLOR = (0, 114, 178)
_DEFAULT_POSITIVE_COLOR = (0, 158, 115)
_DEFAULT_MAX_NEGATIVE = -2.0
_DEFAULT_MAX_POSITIVE = 2.0
_CST_CONTEXT_WORDS = {
    "ACTIVATION",
    "SIGNALING",
    "PATHWAY",
    "CASCADE",
    "RESPONSE",
    "RESPONSES",
    "NETWORK",
    "MODULE",
    "TRANSCRIPTION",
    "EXPRESSION",
    "REGULATION",
}
_TEXT_TOKEN_RE = re.compile(
    r"(?P<Tm>-?\d+(?:\.\d+)?\s+-?\d+(?:\.\d+)?\s+-?\d+(?:\.\d+)?\s+-?\d+(?:\.\d+)?\s+-?\d+(?:\.\d+)?\s+-?\d+(?:\.\d+)?\s+Tm)"
    r"|(?P<Td>-?\d+(?:\.\d+)?\s+-?\d+(?:\.\d+)?\s+Td)"
    r"|(?P<Tj>\((?:\\.|[^\\)])*\)\s*Tj)"
    r"|(?P<TJ>\[(?:.|\n|\r)*?\]\s*TJ)",
    re.S,
)
_ELLIPSE_BLOCK_RE = re.compile(
    r"q 1 0 0 1 ([\-0-9.]+) ([\-0-9.]+) cm\s+0 0 m\s+((?:(?:[\-0-9.]+\s+){6}c\s+){4})f",
    re.S,
)
_PDF_BOX_RE = re.compile(
    rb"/(?:CropBox|MediaBox)\s*\[\s*([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s*\]"
)


def _resolve_base_dir(base_dir: Optional[Path] = None) -> Path:
    return Path(base_dir) if base_dir is not None else Path(__file__).resolve().parent


def _coerce_rgb(value: Sequence[float] | None, fallback: Sequence[float]) -> List[int]:
    src = list(value or fallback)[:3]
    if len(src) < 3:
        src.extend(list(fallback)[: 3 - len(src)])
    output: List[int] = []
    for channel in src[:3]:
        try:
            num = int(round(float(channel)))
        except (TypeError, ValueError):
            num = 128
        output.append(max(0, min(255, num)))
    return output


def _gradient_color_from_fold(
    fold_value: Any,
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
) -> List[int]:
    try:
        if fold_value in (None, "", False):
            return [128, 128, 128]
        fold = float(fold_value)
    except (TypeError, ValueError):
        return [128, 128, 128]
    if not math.isfinite(fold):
        return [128, 128, 128]
    neg = _coerce_rgb(negative_color, _DEFAULT_NEGATIVE_COLOR)
    pos = _coerce_rgb(positive_color, _DEFAULT_POSITIVE_COLOR)
    neg_limit = abs(float(max_negative) if max_negative else 1.0) or 1.0
    pos_limit = float(max_positive) if max_positive else 1.0
    pos_limit = pos_limit if pos_limit != 0 else 1.0
    white = (255, 255, 255)
    if fold < 0:
        t = min(abs(fold) / neg_limit, 1.0)
        return [int((1 - t) * white[i] + t * neg[i]) for i in range(3)]
    t = min(fold / pos_limit, 1.0)
    return [int((1 - t) * white[i] + t * pos[i]) for i in range(3)]


def _normalize_fc_headers(raw_columns: Sequence[Any]) -> List[str]:
    headers: List[str] = []
    for item in raw_columns or []:
        if isinstance(item, (list, tuple)):
            value = item[0] if item else ""
        else:
            value = item
        text = str(value or "").strip()
        if text:
            headers.append(text)
    return headers


def _extract_pdf_literal(raw: str) -> str:
    output: List[str] = []
    idx = 0
    while idx < len(raw):
        char = raw[idx]
        if char == "\\":
            idx += 1
            if idx >= len(raw):
                break
            char = raw[idx]
            if char in {"n", "r", "t", "b", "f"}:
                output.append({"n": "\n", "r": "\r", "t": "\t", "b": "\b", "f": "\f"}[char])
            elif char in {"\\", "(", ")"}:
                output.append(char)
            elif char.isdigit():
                octal = char
                for _ in range(2):
                    if idx + 1 < len(raw) and raw[idx + 1].isdigit():
                        idx += 1
                        octal += raw[idx]
                    else:
                        break
                try:
                    output.append(chr(int(octal, 8)))
                except ValueError:
                    pass
            else:
                output.append(char)
        else:
            output.append(char)
        idx += 1
    return "".join(output)


def _clean_text_label(value: str) -> str:
    text = str(value or "")
    text = text.replace("\x1d", "-").replace("\x1e", "-").replace("\x1f", "-")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _should_merge_cst_text_labels(
    current_label: str,
    next_label: str,
    gap_x: float,
    gap_y: float,
    current_x: float,
    next_x: float,
) -> bool:
    current = _clean_text_label(current_label)
    nxt = _clean_text_label(next_label)
    if not current or not nxt:
        return False
    if gap_y <= 4.5 and -1.0 <= gap_x <= 14.0:
        pass
    elif (
        current.endswith("-")
        and gap_y <= 10.5
        and abs(next_x - current_x) <= 24.0
        and re.fullmatch(r"[A-Za-z0-9/+-]{1,12}", nxt)
    ):
        pass
    else:
        return False
    if re.fullmatch(r"[A-Za-z0-9]{1,2}", nxt):
        return True
    if re.fullmatch(r"[-/]", nxt):
        return True
    if current.endswith(("-", "/", "(", "[")):
        return True
    if re.fullmatch(r"[A-Za-z0-9]+(?:[/-][A-Za-z0-9]+)+", nxt):
        return True
    if re.fullmatch(r"[0-9]+(?:/[0-9]+)+", nxt):
        return True
    if len(current) <= 8 and re.fullmatch(r"[A-Za-z0-9/-]{1,10}", nxt):
        return True
    return False


def _cst_label_variants(label: str) -> List[str]:
    cleaned = _clean_text_label(label)
    variants: List[str] = []
    normalized = cleaned.upper()

    def _add(value: str) -> None:
        text = _clean_text_label(value)
        if text and text not in variants:
            variants.append(text)

    _add(cleaned)
    _add(cleaned.lstrip("-/ "))
    if cleaned.endswith("/"):
        return variants
    if normalized in {"MTORC1", "MTORC2"}:
        return variants
    if any(token in _CST_CONTEXT_WORDS for token in re.findall(r"[A-Z]+", normalized)):
        return variants

    for token in re.findall(r"[A-Za-z][A-Za-z0-9]*(?:[-/][A-Za-z0-9]+)*", cleaned):
        _add(token)
    for chunk in re.split(r"[\s/]+", cleaned):
        _add(chunk)
        match = re.search(r"([A-Z][A-Za-z0-9-]*)$", chunk)
        if match:
            _add(match.group(1))
    return variants


def _iter_decompressed_streams(file_path: Path) -> List[str]:
    data = file_path.read_bytes()
    streams: List[str] = []
    for match in re.finditer(rb"stream\r?\n", data):
        start = match.end()
        end = data.find(b"endstream", start)
        if end < 0:
            continue
        blob = data[start:end]
        if blob.endswith(b"\r\n"):
            blob = blob[:-2]
        elif blob.endswith(b"\n") or blob.endswith(b"\r"):
            blob = blob[:-1]
        try:
            text = zlib.decompress(blob).decode("latin1", errors="replace")
        except Exception:
            continue
        streams.append(text)
    return streams


@lru_cache(maxsize=16)
def _extract_cst_page_size(file_path_str: str) -> Dict[str, float]:
    file_path = Path(file_path_str)
    data = file_path.read_bytes()
    match = _PDF_BOX_RE.search(data)
    if not match:
        return {
            "page_width": _CST_PAGE_WIDTH,
            "page_height": _CST_PAGE_HEIGHT,
        }
    try:
        x0, y0, x1, y1 = [float(part) for part in match.groups()]
        width = max(1.0, x1 - x0)
        height = max(1.0, y1 - y0)
    except Exception:
        width = _CST_PAGE_WIDTH
        height = _CST_PAGE_HEIGHT
    return {
        "page_width": width,
        "page_height": height,
    }


@lru_cache(maxsize=16)
def _extract_cst_text_nodes(file_path_str: str) -> List[Dict[str, Any]]:
    file_path = Path(file_path_str)
    content = ""
    for stream in _iter_decompressed_streams(file_path):
        if "BT" in stream and ("Tj" in stream or "TJ" in stream):
            content = stream
            break
    if not content:
        return []

    font_size = 9.0
    matrix_a = 1.0
    matrix_b = 0.0
    matrix_c = 0.0
    matrix_d = 1.0
    matrix_e = 0.0
    matrix_f = 0.0
    nodes: List[Dict[str, Any]] = []
    order = 0
    for match in _TEXT_TOKEN_RE.finditer(content):
        token = match.group(0)
        if token.endswith("Tm"):
            nums = [float(item) for item in token[:-2].split()]
            font_size = max(abs(nums[0]), abs(nums[3]), font_size or 0.0) or 9.0
            matrix_a, matrix_b, matrix_c, matrix_d, matrix_e, matrix_f = nums
            continue
        if token.endswith("Td"):
            nums = [float(item) for item in token[:-2].split()]
            tx, ty = nums
            matrix_e = matrix_e + matrix_a * tx + matrix_c * ty
            matrix_f = matrix_f + matrix_b * tx + matrix_d * ty
            continue
        if token.endswith("Tj"):
            literal = token[:-2].strip()
            raw_text = _extract_pdf_literal(literal[1:-1]) if literal.startswith("(") and literal.endswith(")") else ""
        else:
            array_text = token[:-2].strip()
            pieces = [
                _extract_pdf_literal(piece.group(0)[1:-1])
                for piece in re.finditer(r"\((?:\\.|[^\\)])*\)", array_text, re.S)
            ]
            raw_text = "".join(pieces)
        cleaned = _clean_text_label(raw_text)
        if not cleaned:
            continue
        est_width = max(28.0, len(cleaned) * max(font_size, 6.0) * 0.56 + 12.0)
        est_height = max(18.0, max(font_size, 6.0) * 1.9)
        nodes.append(
            {
                "order": order,
                "label": cleaned,
                "pdf_x": float(matrix_e),
                "pdf_y": float(matrix_f),
                "font_size": float(max(font_size, 6.0)),
                "estimated_width": float(est_width),
                "estimated_height": float(est_height),
            }
        )
        order += 1

    merged: List[Dict[str, Any]] = []
    idx = 0
    while idx < len(nodes):
        node = dict(nodes[idx])
        while idx + 1 < len(nodes):
            nxt = nodes[idx + 1]
            gap_x = float(nxt["pdf_x"]) - (float(node["pdf_x"]) + float(node["estimated_width"]) * 0.55)
            gap_y = abs(float(nxt["pdf_y"]) - float(node["pdf_y"]))
            if not _should_merge_cst_text_labels(
                str(node["label"]),
                str(nxt["label"]),
                gap_x,
                gap_y,
                float(node["pdf_x"]),
                float(nxt["pdf_x"]),
            ):
                break
            node["label"] = f"{node['label']}{nxt['label']}"
            node["estimated_width"] = max(
                float(node["estimated_width"]),
                (float(nxt["pdf_x"]) - float(node["pdf_x"])) + float(nxt["estimated_width"]),
            )
            node["estimated_height"] = max(float(node["estimated_height"]), float(nxt["estimated_height"]))
            idx += 1
        node["normalized_label"] = _clean_text_label(node["label"]).upper()
        merged.append(node)
        idx += 1
    return merged


@lru_cache(maxsize=16)
def _extract_cst_ellipse_groups(file_path_str: str) -> List[Dict[str, float]]:
    file_path = Path(file_path_str)
    content = ""
    for stream in _iter_decompressed_streams(file_path):
        if "BT" in stream and ("Tj" in stream or "TJ" in stream):
            content = stream
            break
    if not content:
        return []

    groups: List[Dict[str, float]] = []
    for match in _ELLIPSE_BLOCK_RE.finditer(content):
        anchor_x = float(match.group(1))
        anchor_y = float(match.group(2))
        coords = [float(item) for item in re.findall(r"[\-0-9.]+", match.group(3))]
        if len(coords) < 24:
            continue
        xs = coords[0::2]
        ys = coords[1::2]
        all_x = [0.0] + xs
        all_y = [0.0] + ys
        min_x = min(all_x) if all_x else 0.0
        max_x = max(all_x) if all_x else 0.0
        min_y = min(all_y) if all_y else 0.0
        max_y = max(all_y) if all_y else 0.0
        radius_x = (max_x - min_x) * 0.5
        radius_y = (max_y - min_y) * 0.5
        center_x = anchor_x + ((min_x + max_x) * 0.5)
        center_y = anchor_y + ((min_y + max_y) * 0.5)
        if radius_x < 6.0 or radius_y < 6.0:
            continue
        merged = False
        for group in groups:
            if abs(group["center_x"] - center_x) <= 2.5 and abs(group["center_y"] - center_y) <= 2.5:
                group["center_x"] = (group["center_x"] + center_x) * 0.5
                group["center_y"] = (group["center_y"] + center_y) * 0.5
                group["radius_x"] = max(group["radius_x"], radius_x)
                group["radius_y"] = max(group["radius_y"], radius_y)
                merged = True
                break
        if not merged:
            groups.append(
                {
                    "center_x": center_x,
                    "center_y": center_y,
                    "radius_x": radius_x,
                    "radius_y": radius_y,
                }
            )
    return groups


@lru_cache(maxsize=16)
def _map_cst_text_nodes(file_path_str: str) -> List[Dict[str, Any]]:
    nodes = _extract_cst_text_nodes(file_path_str)
    file_path = Path(file_path_str)
    pathway_hint = re.sub(r"\s+\(\d+\)$", "", file_path.stem).strip()
    pathway_mapping = get_cst_pathway_mapping(pathway_hint)
    module_index: Dict[str, Dict[str, Any]] = {}
    for module in list(pathway_mapping.get("modules") or []):
        key = str(module.get("normalized_label") or "").strip().upper()
        if not key:
            continue
        existing = module_index.get(key)
        if existing is None:
            module_index[key] = {
                "suggested_uniprot_ids": list(module.get("uniprot_ids") or []),
                "suggested_gene_symbols": list(module.get("gene_symbols") or []),
                "notes": "Mapped from CST pathway protein-module index.",
            }
        else:
            for field in ("suggested_uniprot_ids", "suggested_gene_symbols"):
                current = list(existing.get(field) or [])
                for item in list(module.get("uniprot_ids" if field == "suggested_uniprot_ids" else "gene_symbols") or []):
                    if item not in current:
                        current.append(item)
                existing[field] = current

    output: List[Dict[str, Any]] = []
    mapper = PathwayLabelMapper(use_uniprot_rest=False)
    mapping_rank = {
        "cst_index": 5,
        "exact_gene": 4,
        "alias": 3,
        "family": 2,
        "ambiguous": 1,
        "unresolved": 0,
    }
    for node in nodes:
        payload = dict(node)
        indexed = None
        indexed_variant = ""
        for variant in _cst_label_variants(str(node.get("label") or "")):
            norm = _clean_text_label(variant).upper()
            candidate = module_index.get(norm)
            if candidate:
                indexed = candidate
                indexed_variant = variant
                break
        if indexed:
            payload["mapping"] = {
                "original_label": node.get("label"),
                "normalized_label": _clean_text_label(indexed_variant).upper(),
                "mapping_type": "cst_index",
                "suggested_gene_symbols": list(indexed.get("suggested_gene_symbols") or []),
                "suggested_uniprot_ids": list(indexed.get("suggested_uniprot_ids") or []),
                "notes": indexed.get("notes") or "Mapped from CST pathway protein-module index.",
                "confidence": "high",
            }
        else:
            best_result: Optional[Dict[str, Any]] = None
            best_variant = ""
            best_rank = -1
            for variant in _cst_label_variants(str(node.get("label") or "")):
                result = mapper.map_pathway_label(variant)
                rank = mapping_rank.get(str(result.get("mapping_type") or "").strip().lower(), 0)
                if rank > best_rank:
                    best_result = result
                    best_variant = variant
                    best_rank = rank
                if rank >= mapping_rank["exact_gene"]:
                    break
            fallback = dict(best_result or mapper.map_pathway_label(str(node.get("label") or "")))
            if best_variant and best_variant != str(node.get("label") or "") and str(fallback.get("mapping_type") or "").strip().lower() != "unresolved":
                fallback["original_label"] = str(node.get("label") or "")
                fallback["normalized_label"] = _clean_text_label(str(best_variant)).upper()
                notes = str(fallback.get("notes") or "").strip()
                fallback["notes"] = f"{notes} Recovered from CST text fragment '{best_variant}'.".strip()
            payload["mapping"] = fallback
        output.append(payload)
    return output


def _build_dataset_index(dataset: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not dataset:
        return {"fc_headers": [], "rows_by_uniprot": {}, "rows_by_gene": {}}
    headers = list(dataset.get("headers") or [])
    rows = list(dataset.get("rows") or [])
    if not headers:
        return {"fc_headers": [], "rows_by_uniprot": {}, "rows_by_gene": {}}
    fc_headers = _normalize_fc_headers(dataset.get("main_columns") or [h for h in headers if str(h).startswith("C:")])
    rows_by_uniprot: Dict[str, Dict[str, Any]] = {}
    rows_by_gene: Dict[str, List[Dict[str, Any]]] = {}
    uniprot_key = headers[0]
    gene_key = headers[1] if len(headers) > 1 else ""
    for row in rows:
        values = list(row)
        if len(values) < len(headers):
            values.extend([""] * (len(headers) - len(values)))
        row_map = {header: values[idx] if idx < len(values) else "" for idx, header in enumerate(headers)}
        raw_uniprot = str(row_map.get(uniprot_key, "") or "").strip().upper()
        raw_gene = str(row_map.get(gene_key, "") or "").strip().upper() if gene_key else ""
        normalized_uniprot = raw_uniprot.split("-", 1)[0] if raw_uniprot else ""
        if normalized_uniprot and normalized_uniprot not in rows_by_uniprot:
            rows_by_uniprot[normalized_uniprot] = row_map
        if raw_gene:
            rows_by_gene.setdefault(raw_gene, []).append(row_map)
    return {
        "fc_headers": fc_headers,
        "rows_by_uniprot": rows_by_uniprot,
        "rows_by_gene": rows_by_gene,
        "headers": headers,
    }


def _match_dataset_rows(mapping: Dict[str, Any], dataset_index: Dict[str, Any]) -> List[Dict[str, Any]]:
    rows_by_uniprot = dataset_index.get("rows_by_uniprot", {})
    rows_by_gene = dataset_index.get("rows_by_gene", {})
    matches: List[Dict[str, Any]] = []
    seen_ids: set[int] = set()
    for uniprot_id in mapping.get("suggested_uniprot_ids", []) or []:
        key = str(uniprot_id or "").strip().upper().split("-", 1)[0]
        row = rows_by_uniprot.get(key)
        if row is None or id(row) in seen_ids:
            continue
        matches.append(row)
        seen_ids.add(id(row))
    for gene_symbol in mapping.get("suggested_gene_symbols", []) or []:
        key = str(gene_symbol or "").strip().upper()
        for row in rows_by_gene.get(key, []):
            if id(row) in seen_ids:
                continue
            matches.append(row)
            seen_ids.add(id(row))
    return matches


def _build_cst_overlay_nodes(
    file_path: Path,
    protein_dataset: Optional[Dict[str, Any]] = None,
    negative_color: Sequence[float] = _DEFAULT_NEGATIVE_COLOR,
    positive_color: Sequence[float] = _DEFAULT_POSITIVE_COLOR,
    max_negative: float = _DEFAULT_MAX_NEGATIVE,
    max_positive: float = _DEFAULT_MAX_POSITIVE,
) -> List[Dict[str, Any]]:
    mapped_nodes = _map_cst_text_nodes(str(file_path))
    if not mapped_nodes:
        return []
    dataset_index = _build_dataset_index(protein_dataset)
    fc_headers = list(dataset_index.get("fc_headers") or [])
    overlay_nodes: List[Dict[str, Any]] = []

    for node in mapped_nodes:
        mapping = dict(node.get("mapping") or {})
        mapping_type = str(mapping.get("mapping_type") or "").strip().lower()
        if mapping_type == "unresolved":
            continue
        matches = _match_dataset_rows(mapping, dataset_index)

        est_width = float(node.get("estimated_width") or 42.0)
        est_height = float(node.get("estimated_height") or 20.0)
        radius_x = max(est_width * 0.38, 12.0)
        radius_y = max(est_height * 0.52, 8.0)

        radius_x = min(radius_x, max(est_width * 0.46, 22.0))
        radius_y = min(radius_y, max(est_height * 0.58, 11.0))

        overlay_node: Dict[str, Any] = {
            "label": node.get("label"),
            "normalized_label": mapping.get("normalized_label") or node.get("normalized_label") or _clean_text_label(str(node.get("label") or "")).upper(),
            "mapping_type": mapping.get("mapping_type"),
            "radius_x": float(radius_x),
            "radius_y": float(radius_y),
            "fallback_x": float(node.get("pdf_x") or 0.0),
            "fallback_y": float(node.get("pdf_y") or 0.0),
            "font_size": float(node.get("font_size") or 9.0),
            "estimated_width": float(est_width),
            "estimated_height": float(est_height),
            "suggested_gene_symbols": list(mapping.get("suggested_gene_symbols") or []),
            "suggested_uniprot_ids": list(mapping.get("suggested_uniprot_ids") or []),
            "has_dataset_match": bool(matches),
            "default_color": [166, 166, 166],
        }

        for idx, header in enumerate(fc_headers, 1):
            chosen_row: Optional[Dict[str, Any]] = None
            chosen_value: Optional[float] = None
            for row in matches:
                raw_value = row.get(header, "")
                try:
                    value = float(raw_value)
                except (TypeError, ValueError):
                    continue
                if chosen_value is None or abs(value) > abs(chosen_value):
                    chosen_value = value
                    chosen_row = row
            if chosen_value is None or chosen_row is None:
                continue
            overlay_node[f"fold_change_{idx}"] = chosen_value
            overlay_node[f"fc_color_{idx}"] = _gradient_color_from_fold(
                chosen_value,
                negative_color,
                positive_color,
                max_negative,
                max_positive,
            )
            if dataset_index.get("headers"):
                uniprot_key = dataset_index["headers"][0]
                gene_key = dataset_index["headers"][1] if len(dataset_index["headers"]) > 1 else ""
                overlay_node[f"matched_uniprot_{idx}"] = str(chosen_row.get(uniprot_key, "") or "")
                overlay_node[f"matched_gene_symbol_{idx}"] = str(chosen_row.get(gene_key, "") or "") if gene_key else ""

        overlay_nodes.append(overlay_node)
    return overlay_nodes


def _summarize_cst_mapping_sources(file_path: Path) -> Dict[str, int]:
    mapped_nodes = _map_cst_text_nodes(str(file_path))
    summary = {
        "total_labels": 0,
        "recognized_total": 0,
        "psp_index_count": 0,
        "backup_count": 0,
        "unresolved_count": 0,
    }
    for node in mapped_nodes:
        summary["total_labels"] += 1
        mapping = dict(node.get("mapping") or {})
        mapping_type = str(mapping.get("mapping_type") or "").strip().lower()
        if not mapping_type or mapping_type == "unresolved":
            summary["unresolved_count"] += 1
            continue
        summary["recognized_total"] += 1
        if mapping_type == "cst_index":
            summary["psp_index_count"] += 1
        else:
            summary["backup_count"] += 1
    return summary


def get_cst_pathway_catalog(base_dir: Optional[Path] = None) -> List[Dict[str, str]]:
    root = _resolve_base_dir(base_dir)
    rows: List[Dict[str, str]] = []
    for entry in _DEFAULT_CST_FILES:
        file_path = root / entry["filename"]
        if not file_path.exists():
            continue
        rows.append(
            {
                "id": entry["id"],
                "name": entry["name"],
                "filename": entry["filename"],
                "file_path": str(file_path),
                "source": "CST",
            }
        )
    return rows


def load_cst_pathway_payload(
    pathway_id: str,
    base_dir: Optional[Path] = None,
    protein_dataset: Optional[Dict[str, Any]] = None,
    negative_color: Sequence[float] = _DEFAULT_NEGATIVE_COLOR,
    positive_color: Sequence[float] = _DEFAULT_POSITIVE_COLOR,
    max_negative: float = _DEFAULT_MAX_NEGATIVE,
    max_positive: float = _DEFAULT_MAX_POSITIVE,
) -> Optional[Dict[str, Any]]:
    pathway_key = str(pathway_id or "").strip().lower()
    if not pathway_key:
        return None
    for entry in get_cst_pathway_catalog(base_dir):
        if entry["id"].strip().lower() != pathway_key:
            continue
        file_path = Path(entry["file_path"])
        pdf_bytes = file_path.read_bytes()
        page_size = _extract_cst_page_size(str(file_path))
        data_uri = "data:application/pdf;base64," + base64.b64encode(pdf_bytes).decode("ascii")
        return {
            "id": entry["id"],
            "name": entry["name"],
            "filename": entry["filename"],
            "file_path": str(file_path),
            "mime_type": "application/pdf",
            "data_uri": data_uri,
            "pdf_base64": base64.b64encode(pdf_bytes).decode("ascii"),
            "page_width": float(page_size.get("page_width") or _CST_PAGE_WIDTH),
            "page_height": float(page_size.get("page_height") or _CST_PAGE_HEIGHT),
            "ellipse_groups": _extract_cst_ellipse_groups(str(file_path)),
            "mapping_summary": _summarize_cst_mapping_sources(file_path),
            "overlay_nodes": _build_cst_overlay_nodes(
                file_path,
                protein_dataset=protein_dataset,
                negative_color=negative_color,
                positive_color=positive_color,
                max_negative=max_negative,
                max_positive=max_positive,
            ),
        }
    return None


def create_cst_pathway_viewer(payload: Optional[Dict[str, Any]]) -> Any:
    info = payload or {}
    data_uri = str(info.get("data_uri") or "").strip()
    pdf_base64 = str(info.get("pdf_base64") or "").strip()
    title = str(info.get("name") or info.get("filename") or "CST Pathway").strip()
    page_width = float(info.get("page_width") or _CST_PAGE_WIDTH)
    page_height = float(info.get("page_height") or _CST_PAGE_HEIGHT)
    active_idx = max(1, int(info.get("_active_fc_index") or 1))
    overlay_nodes = list(info.get("overlay_nodes") or [])
    ellipse_groups = list(info.get("ellipse_groups") or [])
    mapping_summary = dict(info.get("mapping_summary") or {})
    overlay_nodes_json = json.dumps(overlay_nodes)
    ellipse_groups_json = json.dumps(ellipse_groups)
    viewer_key = re.sub(r"[^A-Za-z0-9_-]+", "_", str(info.get("id") or "cst_viewer"))
    stage_id = f"cst-stage-{viewer_key}"
    canvas_id = f"cst-canvas-{viewer_key}"
    fallback_id = f"cst-fallback-{viewer_key}"
    overlay_id = f"cst-overlay-{viewer_key}"
    missing_button_id = f"cst-missing-btn-{viewer_key}"
    coord_tooltip_id = f"cst-coords-{viewer_key}"
    stats_id = f"cst-stats-{viewer_key}"
    if not data_uri or not pdf_base64:
        return ui.div({"class": "alert alert-warning"}, "No CST pathway file is available to display.")

    overlay_markup = ui.HTML(
        f'<svg id="{overlay_id}" class="cst-viewer-overlay-svg" viewBox="0 0 {page_width:.3f} {page_height:.3f}" preserveAspectRatio="xMidYMid meet" aria-hidden="true"></svg>'
    )

    return ui.TagList(
        ui.tags.style(
            """
            .cst-viewer-shell {
                width: 100%;
                position: relative;
                border-radius: 18px;
                overflow: hidden;
                background: linear-gradient(180deg, #f8fafc 0%, #eef2ff 100%);
                border: 1px solid rgba(148, 163, 184, 0.28);
            }
            .cst-viewer-stage {
                position: relative;
                width: 100%;
                aspect-ratio: 612 / 699.627;
                background: #ffffff;
            }
            .cst-viewer-canvas {
                display: block;
                width: 100%;
                height: 100%;
                background: #ffffff;
            }
            .cst-viewer-fallback {
                position: absolute;
                inset: 0;
                width: 100%;
                height: 100%;
                border: 0;
                background: #ffffff;
                display: none;
                z-index: 1;
            }
            .cst-viewer-overlay {
                position: absolute;
                inset: 0;
                pointer-events: auto;
                z-index: 2;
            }
            .cst-viewer-toolbar {
                position: absolute;
                top: 12px;
                right: 12px;
                z-index: 3;
                display: flex;
                gap: 8px;
            }
            .cst-viewer-stats {
                position: absolute;
                top: 12px;
                left: 12px;
                z-index: 3;
                display: flex;
                align-items: center;
                gap: 10px;
                padding: 8px 12px;
                border-radius: 12px;
                background: rgba(255, 255, 255, 0.92);
                border: 1px solid rgba(148, 163, 184, 0.28);
                box-shadow: 0 8px 22px rgba(15, 23, 42, 0.08);
                backdrop-filter: blur(8px);
            }
            .cst-viewer-stats-title {
                font-size: 10px;
                font-weight: 800;
                letter-spacing: 0.08em;
                text-transform: uppercase;
                color: #64748b;
                white-space: nowrap;
            }
            .cst-viewer-stats-grid {
                display: flex;
                align-items: center;
                gap: 12px;
            }
            .cst-viewer-stat {
                display: flex;
                flex-direction: column;
                gap: 1px;
                min-width: 0;
            }
            .cst-viewer-stat-label {
                font-size: 10px;
                font-weight: 700;
                color: #64748b;
                white-space: nowrap;
            }
            .cst-viewer-stat-value {
                font-size: 14px;
                font-weight: 800;
                color: #0f172a;
                line-height: 1.1;
                white-space: nowrap;
            }
            .cst-viewer-action {
                border: 1px solid rgba(100, 116, 139, 0.28);
                background: rgba(255, 255, 255, 0.88);
                color: #0f172a;
                border-radius: 10px;
                padding: 6px 10px;
                font-size: 12px;
                font-weight: 600;
                line-height: 1;
                backdrop-filter: blur(8px);
                box-shadow: 0 8px 22px rgba(15, 23, 42, 0.08);
            }
            .cst-viewer-action:hover {
                background: rgba(255, 255, 255, 0.98);
            }
            .cst-viewer-coords {
                position: absolute;
                z-index: 4;
                display: none;
                pointer-events: none;
                padding: 4px 6px;
                border-radius: 8px;
                background: rgba(15, 23, 42, 0.88);
                color: #f8fafc;
                font-size: 11px;
                line-height: 1.1;
                font-family: Consolas, "Courier New", monospace;
                white-space: nowrap;
                transform: translate(12px, 12px);
            }
            .cst-viewer-overlay-svg {
                width: 100%;
                height: 100%;
                display: block;
                overflow: visible;
                pointer-events: none;
            }
            .cst-viewer-overlay-svg ellipse {
                pointer-events: stroke;
                cursor: help;
            }
            .cst-viewer-overlay-svg .cst-missing-node {
                display: none;
            }
            .cst-viewer-stage.cst-show-missing .cst-viewer-overlay-svg .cst-missing-node {
                display: inline;
            }
            """
        ),
        ui.tags.script({"src": "https://cdnjs.cloudflare.com/ajax/libs/pdf.js/3.11.174/pdf.min.js"}),
        ui.div(
            {"class": "cst-viewer-shell"},
            ui.div(
                {"class": "cst-viewer-stage", "id": stage_id},
                ui.div(
                    {"class": "cst-viewer-stats", "id": stats_id},
                    ui.div({"class": "cst-viewer-stats-title"}, "Protein Mapping"),
                    ui.div(
                        {"class": "cst-viewer-stats-grid"},
                        ui.div(
                            {"class": "cst-viewer-stat", "title": "Protein labels recognized from the CST/PSP-backed pathway annotation file."},
                            ui.div({"class": "cst-viewer-stat-label"}, "PSP"),
                            ui.div({"class": "cst-viewer-stat-value"}, str(int(mapping_summary.get("psp_index_count") or 0))),
                        ),
                        ui.div(
                            {"class": "cst-viewer-stat", "title": "Protein labels recognized using the backup gene-label to UniProt mapper."},
                            ui.div({"class": "cst-viewer-stat-label"}, "Backup"),
                            ui.div({"class": "cst-viewer-stat-value"}, str(int(mapping_summary.get("backup_count") or 0))),
                        ),
                        ui.div(
                            {"class": "cst-viewer-stat", "title": "Total recognized protein labels in this CST pathway."},
                            ui.div({"class": "cst-viewer-stat-label"}, "Total"),
                            ui.div({"class": "cst-viewer-stat-value"}, str(int(mapping_summary.get("recognized_total") or 0))),
                        ),
                    ),
                ),
                ui.div(
                    {"class": "cst-viewer-toolbar"},
                    ui.tags.button(
                        {
                            "id": missing_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                        },
                        "Show Missing Proteins",
                    ),
                ),
                ui.div({"class": "cst-viewer-coords", "id": coord_tooltip_id}),
                ui.tags.canvas({"class": "cst-viewer-canvas", "id": canvas_id}),
                ui.tags.iframe(
                    {
                        "class": "cst-viewer-fallback",
                        "id": fallback_id,
                        "src": f"{data_uri}#toolbar=0&navpanes=0&view=FitH",
                        "title": title,
                    }
                ),
                ui.div({"class": "cst-viewer-overlay"}, overlay_markup),
            ),
        ),
        ui.tags.script(
            f"""
            (function() {{
                const stage = document.getElementById('{stage_id}');
                const canvas = document.getElementById('{canvas_id}');
                const fallback = document.getElementById('{fallback_id}');
                const overlaySvg = document.getElementById('{overlay_id}');
                const missingButton = document.getElementById('{missing_button_id}');
                const coordsTooltip = document.getElementById('{coord_tooltip_id}');
                if (!stage || !canvas) return;
                const renderStateKey = '__mapkinaseCstRenderState';
                const pdfData = '{pdf_base64}';
                const pageWidth = {page_width:.6f};
                const pageHeight = {page_height:.6f};
                const overlayData = {overlay_nodes_json};
                const ellipseGroups = {ellipse_groups_json};
                const svgNs = 'http://www.w3.org/2000/svg';
                const showFallback = () => {{
                    if (fallback) fallback.style.display = 'block';
                    canvas.style.display = 'none';
                }};
                const cleanLabel = (value) => String(value || '')
                    .replace(/[\\x1d\\x1e\\x1f]/g, '-')
                    .replace(/\\s+/g, ' ')
                    .trim();
                const labelVariants = (value) => {{
                    const cleaned = cleanLabel(value);
                    const variants = [];
                    const add = (item) => {{
                        const text = cleanLabel(item);
                        if (text && !variants.includes(text)) variants.push(text);
                    }};
                    add(cleaned);
                    add(cleaned.replace(/^[-/ ]+/, ''));
                    const tokenMatches = cleaned.match(/[A-Za-z][A-Za-z0-9]*(?:[-/][A-Za-z0-9]+)*/g) || [];
                    for (const token of tokenMatches) add(token);
                    for (const chunk of cleaned.split(/[\\s/]+/)) {{
                        add(chunk);
                        const suffix = chunk.match(/([A-Z][A-Za-z0-9-]*)$/);
                        if (suffix) add(suffix[1]);
                    }}
                    return variants.map((item) => cleanLabel(item).toUpperCase()).filter(Boolean);
                }};
                const buildMergedTextItems = (items) => {{
                    const merged = [];
                    for (const item of items || []) {{
                        const label = cleanLabel(item.str);
                        if (!label) continue;
                        const transform = Array.isArray(item.transform) ? item.transform : [1, 0, 0, 1, 0, 0];
                        const x = Number(transform[4] || 0);
                        const y = Number(transform[5] || 0);
                        const fontSize = Math.max(Math.abs(Number(transform[0] || 0)), Math.abs(Number(transform[3] || 0)), 6);
                        const width = Math.max(Number(item.width || 0), label.length * fontSize * 0.42);
                        const height = Math.max(Number(item.height || 0), fontSize);
                        const previous = merged.length ? merged[merged.length - 1] : null;
                        if (
                            previous &&
                            Math.abs(previous.y - y) <= 4 &&
                            x - (previous.x + previous.width) <= 12 &&
                            /^[A-Za-z0-9]{1,2}$/.test(label)
                        ) {{
                            previous.label += label;
                            previous.normalized = cleanLabel(previous.label).toUpperCase();
                            previous.width = Math.max(previous.width, (x - previous.x) + width);
                            previous.height = Math.max(previous.height, height);
                            continue;
                        }}
                        merged.push({{
                            label,
                            normalized: cleanLabel(label).toUpperCase(),
                            x,
                            y,
                            width,
                            height,
                            fontSize,
                        }});
                    }}
                    return merged;
                }};
                const renderOverlayNodes = (textItems) => {{
                    if (!overlaySvg) return;
                    overlaySvg.replaceChildren();
                    const available = new Map();
                    const availableEllipses = (ellipseGroups || []).map((ellipse, index) => ({{ ...ellipse, _used: false, _idx: index }}));
                    for (const item of buildMergedTextItems(textItems)) {{
                        for (const key of labelVariants(item.label || item.normalized || '')) {{
                            if (!available.has(key)) available.set(key, []);
                            available.get(key).push(item);
                        }}
                    }}
                    const getItemGeometry = (item, node) => {{
                        if (item) {{
                            const width = Number(item.width || 0);
                            const height = Number(item.height || 0);
                            const x = Number(item.x || 0);
                            const y = Number(item.y || 0);
                            const label = cleanLabel(item.label || '');
                            const fontSize = Math.max(Number(item.fontSize || 0), 6);
                            const ownerWidth = Math.min(
                                width,
                                Math.max(12, fontSize * Math.min(Math.max(label.length, 1), 3.5) * 0.55)
                            );
                            const centerX = x + (width * 0.5);
                            const centerY = y + (height * 0.45);
                            return {{
                                centerX,
                                centerY,
                                expectedWidth: width * 0.55,
                                anchorX: x,
                                anchorY: centerY,
                                left: x - 8,
                                right: x + ownerWidth,
                                top: centerY - (height * 0.72),
                                bottom: centerY + (height * 0.42),
                                hasBounds: true,
                            }};
                        }}
                        return {{
                            centerX: Number(node.fallback_x || 0),
                            centerY: Number(node.fallback_y || 0),
                            expectedWidth: Number(node.estimated_width || 0) * 0.38,
                            anchorX: Number(node.fallback_x || 0),
                            anchorY: Number(node.fallback_y || 0),
                            left: Number(node.fallback_x || 0),
                            right: Number(node.fallback_x || 0),
                            top: Number(node.fallback_y || 0),
                            bottom: Number(node.fallback_y || 0),
                            hasBounds: false,
                        }};
                    }};
                    const getFallbackEllipseGeometry = (item, node) => {{
                        const geom = getItemGeometry(item, node);
                        const rawLabel = cleanLabel((item && item.label) || node.label || '');
                        let shiftX = 0;
                        if (/^[^A-Za-z0-9]+/.test(rawLabel)) {{
                            shiftX = Math.min(
                                Math.max(Number(node.radius_x || 0), 12) * 0.45,
                                Math.max(Number(node.estimated_width || 0), Number(item && item.width || 0), 24) * 0.18
                            );
                        }} else if ((rawLabel.match(/[A-Za-z]/g) || []).length <= 4) {{
                            shiftX = 0;
                        }} else {{
                            shiftX = Math.min(4, Math.max(Number(node.radius_x || 0), 12) * 0.16);
                        }}
                        const fallbackCx = Number((item ? item.x : node.fallback_x) || 0) + shiftX;
                        const fallbackCyPdf = item
                            ? (Number(item.y || 0) + (Number(item.height || 0) * 0.42))
                            : (Number(node.fallback_y || 0) + (Number(node.estimated_height || 18) * 0.42));
                        return {{
                            cx: fallbackCx,
                            cyPdf: fallbackCyPdf,
                            rx: Math.max(Number(node.radius_x || 0), 12),
                            ry: Math.max(Number(node.radius_y || 0), 9),
                        }};
                    }};
                    const findNearestEllipse = (item, node, options = {{}}) => {{
                        const geom = getItemGeometry(item, node);
                        const centerX = geom.centerX;
                        const centerY = geom.centerY;
                        const anchorX = geom.anchorX;
                        const anchorY = geom.anchorY;
                        const expectedWidth = geom.expectedWidth;
                        if (!Number.isFinite(centerX) || !Number.isFinite(centerY)) return {{ ellipse: null, score: Infinity }};
                        let best = null;
                        let bestScore = Infinity;
                        const requireOverlap = !!options.requireOverlap;
                        for (const ellipse of availableEllipses) {{
                            if (ellipse._used) continue;
                            const dx = Number(ellipse.center_x || 0) - centerX;
                            const dy = Number(ellipse.center_y || 0) - centerY;
                            const dist = Math.hypot(dx, dy);
                            if (dist > 42) continue;
                            const rx = Number(ellipse.radius_x || 0);
                            const ry = Number(ellipse.radius_y || 0);
                            if (rx < 7 || ry < 6) continue;
                            const anchorDx = Math.abs(Number(ellipse.center_x || 0) - anchorX);
                            const anchorDy = Math.abs(Number(ellipse.center_y || 0) - anchorY);
                            let effectiveDist = dist;
                            let overlapScore = ((dx / Math.max(rx, 1)) ** 2) + ((dy / Math.max(ry, 1)) ** 2);
                            if (geom.hasBounds) {{
                                const closestX = Math.max(geom.left, Math.min(Number(ellipse.center_x || 0), geom.right));
                                const closestY = Math.max(geom.top, Math.min(Number(ellipse.center_y || 0), geom.bottom));
                                const boxDx = Number(ellipse.center_x || 0) - closestX;
                                const boxDy = Number(ellipse.center_y || 0) - closestY;
                                effectiveDist = Math.hypot(boxDx, boxDy);
                                const outsideX = Number(ellipse.center_x || 0) < geom.left
                                    ? (geom.left - Number(ellipse.center_x || 0))
                                    : Math.max(0, Number(ellipse.center_x || 0) - geom.right);
                                const outsideY = Number(ellipse.center_y || 0) < geom.top
                                    ? (geom.top - Number(ellipse.center_y || 0))
                                    : Math.max(0, Number(ellipse.center_y || 0) - geom.bottom);
                                overlapScore = ((outsideX / Math.max(rx, 1)) ** 2) + ((outsideY / Math.max(ry, 1)) ** 2);
                                if (
                                    requireOverlap &&
                                    (
                                        Number(ellipse.center_x || 0) < (geom.left - (rx * 0.6)) ||
                                        Number(ellipse.center_x || 0) > (geom.right + (rx * 0.6)) ||
                                        Number(ellipse.center_y || 0) < (geom.top - (ry * 0.7)) ||
                                        Number(ellipse.center_y || 0) > (geom.bottom + (ry * 0.7))
                                    )
                                ) {{
                                    continue;
                                }}
                            }} else if (requireOverlap && overlapScore > 1.35) {{
                                continue;
                            }}
                            const sizePenalty = Math.abs(rx - expectedWidth) * 0.18;
                            const score = (effectiveDist * 1.6) + (anchorDx * 0.55) + (anchorDy * 0.12) + sizePenalty + (overlapScore * 8) + (dist * 0.04);
                            if (score < bestScore) {{
                                bestScore = score;
                                best = ellipse;
                            }}
                        }}
                        return {{ ellipse: best, score: bestScore }};
                    }};
                    const buildNodeAssignments = () => {{
                        const candidateEdges = [];
                        const fallbackByNode = new Map();
                        const strongEdges = [];
                        for (let nodeIndex = 0; nodeIndex < overlayData.length; nodeIndex += 1) {{
                            const node = overlayData[nodeIndex];
                            const normalized = cleanLabel(node.normalized_label || node.label).toUpperCase();
                            const candidates = Array.isArray(available.get(normalized)) ? available.get(normalized) : [];
                            let bestFallback = null;
                            let bestFallbackScore = Infinity;
                            for (const item of candidates) {{
                                const strict = findNearestEllipse(item, node, {{ requireOverlap: true }});
                                const candidate = strict.ellipse ? strict : findNearestEllipse(item, node);
                                if (!candidate.ellipse) continue;
                                candidateEdges.push({{
                                    nodeIndex,
                                    item,
                                    ellipse: candidate.ellipse,
                                    score: candidate.score,
                                }});
                                if (candidate.score < bestFallbackScore) {{
                                    bestFallbackScore = candidate.score;
                                    bestFallback = {{ item, ellipse: candidate.ellipse, score: candidate.score }};
                                }}
                            }}
                            if (!bestFallback) {{
                                const fallback = findNearestEllipse(null, node);
                                if (fallback.ellipse) {{
                                    bestFallback = {{ item: null, ellipse: fallback.ellipse, score: fallback.score }};
                                }}
                            }}
                            if (bestFallback) {{
                                fallbackByNode.set(nodeIndex, bestFallback);
                                if (bestFallback.item) {{
                                    const geom = getItemGeometry(bestFallback.item, node);
                                    const strongDx = Math.abs(Number(bestFallback.ellipse.center_x || 0) - Number(geom.anchorX || 0));
                                    const strongDy = Math.abs(Number(bestFallback.ellipse.center_y || 0) - Number(geom.anchorY || 0));
                                    if (strongDx <= 8 && strongDy <= 12 && bestFallback.score <= 8) {{
                                        strongEdges.push({{
                                            nodeIndex,
                                            item: bestFallback.item,
                                            ellipse: bestFallback.ellipse,
                                            score: bestFallback.score,
                                            anchorScore: strongDx + (strongDy * 0.5),
                                        }});
                                    }}
                                }}
                            }}
                        }}
                        candidateEdges.sort((a, b) => a.score - b.score);
                        const assignedNodes = new Set();
                        const assignedEllipses = new Set();
                        const assignments = new Map();
                        strongEdges.sort((a, b) => a.anchorScore - b.anchorScore);
                        for (const edge of strongEdges) {{
                            if (assignedNodes.has(edge.nodeIndex)) continue;
                            if (assignedEllipses.has(edge.ellipse._idx)) continue;
                            assignedNodes.add(edge.nodeIndex);
                            assignedEllipses.add(edge.ellipse._idx);
                            assignments.set(edge.nodeIndex, {{ item: edge.item, ellipse: edge.ellipse }});
                        }}
                        for (const edge of candidateEdges) {{
                            if (assignedNodes.has(edge.nodeIndex)) continue;
                            if (assignedEllipses.has(edge.ellipse._idx)) continue;
                            assignedNodes.add(edge.nodeIndex);
                            assignedEllipses.add(edge.ellipse._idx);
                            assignments.set(edge.nodeIndex, {{ item: edge.item, ellipse: edge.ellipse }});
                        }}
                        for (let nodeIndex = 0; nodeIndex < overlayData.length; nodeIndex += 1) {{
                            if (assignments.has(nodeIndex)) continue;
                            const fallback = fallbackByNode.get(nodeIndex);
                            if (!fallback) continue;
                            if (assignedEllipses.has(fallback.ellipse._idx)) continue;
                            assignedEllipses.add(fallback.ellipse._idx);
                            assignments.set(nodeIndex, {{ item: fallback.item, ellipse: fallback.ellipse }});
                        }}
                        return assignments;
                    }};
                    const assignments = buildNodeAssignments();
                    for (let nodeIndex = 0; nodeIndex < overlayData.length; nodeIndex += 1) {{
                        const node = overlayData[nodeIndex];
                        const match = assignments.get(nodeIndex);
                        const item = match ? match.item : null;
                        const color = node['fc_color_{active_idx}'] || node.fc_color_1 || node.default_color || [166, 166, 166];
                        const rgb = Array.isArray(color) ? color.slice(0, 3) : [166, 166, 166];
                        const ellipseMatch = match ? match.ellipse : null;
                        const fallbackGeom = getFallbackEllipseGeometry(item, node);
                        const useFallbackEllipse = !ellipseMatch || (
                            Math.hypot(
                                Number(ellipseMatch.center_x || 0) - Number(fallbackGeom.cx || 0),
                                Number(ellipseMatch.center_y || 0) - Number(fallbackGeom.cyPdf || 0)
                            ) > Math.max(18, Math.max(Number(fallbackGeom.rx || 12), Number(fallbackGeom.ry || 9)) * 1.75)
                        );
                        const rx = useFallbackEllipse
                            ? Number(fallbackGeom.rx || 12)
                            : Number(ellipseMatch.radius_x || 0);
                        const ry = useFallbackEllipse
                            ? Number(fallbackGeom.ry || 9)
                            : Number(ellipseMatch.radius_y || 0);
                        const strokeWidth = Math.max(4, Math.min(11, Math.min(rx, ry) * 0.46));
                        const cx = useFallbackEllipse
                            ? Number(fallbackGeom.cx || 0)
                            : Number(ellipseMatch.center_x || 0);
                        const cyPdf = useFallbackEllipse
                            ? Number(fallbackGeom.cyPdf || 0)
                            : Number(ellipseMatch.center_y || 0);
                        const cy = pageHeight - cyPdf;
                        const label = String(node.label || '');
                        const matchGene = String(node['matched_gene_symbol_{active_idx}'] || node.matched_gene_symbol_1 || '');
                        const matchUniprot = String(node['matched_uniprot_{active_idx}'] || node.matched_uniprot_1 || '');
                        const foldValue = node['fold_change_{active_idx}'];
                        const foldText = Number.isFinite(Number(foldValue)) ? Number(foldValue).toFixed(3) : '';
                        const titleParts = ['Pathway label: ' + (label || matchGene || matchUniprot)];
                        if (matchGene && matchGene !== label) titleParts.push('Mapped gene: ' + matchGene);
                        if (matchUniprot) titleParts.push('UniProt: ' + matchUniprot);
                        if (foldText) titleParts.push('Fold change: ' + foldText);
                        const ellipse = document.createElementNS(svgNs, 'ellipse');
                        ellipse.setAttribute('class', node.has_dataset_match ? 'cst-overlay-ellipse' : 'cst-overlay-ellipse cst-missing-node');
                        ellipse.setAttribute('cx', cx.toFixed(3));
                        ellipse.setAttribute('cy', cy.toFixed(3));
                        ellipse.setAttribute('rx', rx.toFixed(3));
                        ellipse.setAttribute('ry', ry.toFixed(3));
                        ellipse.setAttribute('fill', 'none');
                        ellipse.setAttribute('stroke', 'rgb(' + rgb[0] + ', ' + rgb[1] + ', ' + rgb[2] + ')');
                        ellipse.setAttribute('stroke-width', strokeWidth.toFixed(3));
                        ellipse.setAttribute('stroke-linecap', 'round');
                        ellipse.setAttribute('stroke-linejoin', 'round');
                        ellipse.setAttribute('opacity', '0.98');
                        const title = document.createElementNS(svgNs, 'title');
                        title.textContent = titleParts.join(' | ');
                        ellipse.appendChild(title);
                        overlaySvg.appendChild(ellipse);
                    }}
                }};
                const setMissingState = (showMissing) => {{
                    stage.classList.toggle('cst-show-missing', !!showMissing);
                    if (missingButton) {{
                        missingButton.textContent = showMissing ? 'Hide Missing Proteins' : 'Show Missing Proteins';
                    }}
                }};
                setMissingState(false);
                if (missingButton) {{
                    missingButton.addEventListener('click', () => {{
                        setMissingState(!stage.classList.contains('cst-show-missing'));
                    }});
                }}
                const hideCoords = () => {{
                    if (coordsTooltip) coordsTooltip.style.display = 'none';
                }};
                const updateCoords = (evt) => {{
                    if (!coordsTooltip) return;
                    const isEllipse = evt.target && evt.target.tagName && String(evt.target.tagName).toLowerCase() === 'ellipse';
                    const isToolbar = evt.target && evt.target.closest && evt.target.closest('.cst-viewer-toolbar');
                    if (isEllipse || isToolbar) {{
                        hideCoords();
                        return;
                    }}
                    const rect = stage.getBoundingClientRect();
                    if (!rect.width || !rect.height) {{
                        hideCoords();
                        return;
                    }}
                    const localX = evt.clientX - rect.left;
                    const localY = evt.clientY - rect.top;
                    if (localX < 0 || localY < 0 || localX > rect.width || localY > rect.height) {{
                        hideCoords();
                        return;
                    }}
                    const pdfX = (localX / rect.width) * pageWidth;
                    const pdfY = pageHeight - ((localY / rect.height) * pageHeight);
                    coordsTooltip.textContent = 'x: ' + pdfX.toFixed(1) + '  y: ' + pdfY.toFixed(1);
                    coordsTooltip.style.left = localX + 'px';
                    coordsTooltip.style.top = localY + 'px';
                    coordsTooltip.style.display = 'block';
                }};
                stage.addEventListener('mousemove', updateCoords);
                stage.addEventListener('mouseleave', hideCoords);
                if (overlaySvg) {{
                    overlaySvg.addEventListener('mousemove', updateCoords);
                    overlaySvg.addEventListener('mouseleave', hideCoords);
                }}
                const ensurePdfJs = () => window.pdfjsLib && window.pdfjsLib.getDocument;
                const decodeBase64 = (value) => {{
                    const binary = atob(value);
                    const bytes = new Uint8Array(binary.length);
                    for (let i = 0; i < binary.length; i += 1) {{
                        bytes[i] = binary.charCodeAt(i);
                    }}
                    return bytes;
                }};
                const setupAndRender = () => {{
                    if (!ensurePdfJs()) {{
                        showFallback();
                        return;
                    }}
                    if (!window[renderStateKey]) window[renderStateKey] = {{}};
                    const state = window[renderStateKey];
                    if (!state['{viewer_key}']) {{
                        state['{viewer_key}'] = {{
                            pdfPromise: window.pdfjsLib.getDocument({{ data: decodeBase64(pdfData) }}).promise,
                            renderTask: null,
                            rafId: null,
                            resizeObserver: null,
                            renderVersion: 0,
                        }};
                    }}
                    const local = state['{viewer_key}'];
                    window.pdfjsLib.GlobalWorkerOptions.workerSrc = 'https://cdnjs.cloudflare.com/ajax/libs/pdf.js/3.11.174/pdf.worker.min.js';
                    const renderNow = () => {{
                        local.renderVersion += 1;
                        const renderVersion = local.renderVersion;
                        local.pdfPromise.then((pdf) => pdf.getPage(1)).then((page) => {{
                            const stageWidth = Math.max(stage.clientWidth || 0, 10);
                            const cssScale = stageWidth / pageWidth;
                            const cssHeight = cssScale * pageHeight;
                            stage.style.height = cssHeight + 'px';
                            const dpr = Math.max(window.devicePixelRatio || 1, 1);
                            const viewport = page.getViewport({{ scale: cssScale * dpr }});
                            const ctx = canvas.getContext('2d', {{ alpha: false }});
                            if (!ctx) {{
                                showFallback();
                                return;
                            }}
                            canvas.width = Math.max(1, Math.round(viewport.width));
                            canvas.height = Math.max(1, Math.round(viewport.height));
                            canvas.style.width = stageWidth + 'px';
                            canvas.style.height = cssHeight + 'px';
                            canvas.style.display = 'block';
                            if (fallback) fallback.style.display = 'none';
                            if (local.renderTask && local.renderTask.cancel) {{
                                try {{ local.renderTask.cancel(); }} catch (_err) {{}}
                            }}
                            local.renderTask = page.render({{ canvasContext: ctx, viewport }});
                            return Promise.all([
                                local.renderTask.promise.catch((err) => {{
                                    if (err && err.name === 'RenderingCancelledException') return null;
                                    showFallback();
                                    return null;
                                }}),
                                page.getTextContent().catch(() => ({{ items: [] }})),
                            ]).then(([, textContent]) => {{
                                if (renderVersion !== local.renderVersion) return null;
                                renderOverlayNodes((textContent && textContent.items) || []);
                                return null;
                            }}).catch((err) => {{
                                if (err && err.name === 'RenderingCancelledException') return null;
                                showFallback();
                                return null;
                            }});
                        }}).catch(() => {{
                            showFallback();
                        }});
                    }};
                    const queueRender = () => {{
                        if (local.rafId) cancelAnimationFrame(local.rafId);
                        local.rafId = requestAnimationFrame(() => {{
                            local.rafId = null;
                            renderNow();
                        }});
                    }};
                    queueRender();
                    if (!local.resizeObserver && window.ResizeObserver) {{
                        local.resizeObserver = new ResizeObserver(() => queueRender());
                        local.resizeObserver.observe(stage);
                    }}
                    window.addEventListener('resize', queueRender, {{ passive: true }});
                }};
                const waitForPdfJs = (attemptsLeft) => {{
                    if (ensurePdfJs()) {{
                        setupAndRender();
                        return;
                    }}
                    if (attemptsLeft <= 0) {{
                        showFallback();
                        return;
                    }}
                    setTimeout(() => waitForPdfJs(attemptsLeft - 1), 120);
                }};
                waitForPdfJs(40);
            }})();
            """
        ),
    )
