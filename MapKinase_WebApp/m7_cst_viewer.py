from __future__ import annotations

import base64
import contextvars
import hashlib
import html
import json
import math
import re
import zlib
from datetime import UTC, datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from shiny import ui

from MapKinase_WebApp.m4_json import (
    DEFAULT_DATA,
    _display_metabolite_name,
    _load_hmdb_metabolite_reference,
    _normalize_metabolite_token,
    _resolve_symbol_icon,
    classify_phosphosite_function,
)
from MapKinase_WebApp.m8_pathway_label_mapper import (
    PathwayLabelMapper,
    _canonicalize_label,
    load_complex_review_mapping_index,
    load_forced_review_mapping_index,
)
from MapKinase_WebApp.m11_cst_pathway_index import get_cst_pathway_mapping

_CST_PAGE_WIDTH = 612.0
_CST_PAGE_HEIGHT = 699.627
_DEFAULT_NEGATIVE_COLOR = (0, 114, 178)
_DEFAULT_POSITIVE_COLOR = (0, 158, 115)
_DEFAULT_MAX_NEGATIVE = -2.0
_DEFAULT_MAX_POSITIVE = 2.0
HUMAN_ORTHOLOG_UNIPROT_COLUMN = "Human_Ortholog_UniProt"
HUMAN_ORTHOLOG_GENE_COLUMN = "Human_Ortholog_Gene"
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
_RECT_BLOCK_RE = re.compile(
    r"(?<![\d.])([\-0-9.]+)\s+([\-0-9.]+)\s+([\-0-9.]+)\s+([\-0-9.]+)\s+re\s+(?:B\*?|b\*?|S|s|f\*?|F)\b",
    re.S,
)
_PDF_TOKEN_RE = re.compile(r"\[[^\]]*\]|-?\d+(?:\.\d+)?|[A-Za-z][A-Za-z*]*")
_PDF_BOX_RE = re.compile(
    rb"/(?:CropBox|MediaBox)\s*\[\s*([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s*\]"
)
_GRADIENT_STOPS_CTX: contextvars.ContextVar[Optional[List[Dict[str, Any]]]] = contextvars.ContextVar(
    "cst_gradient_stops", default=None
)


def _json_for_inline_script(value: Any) -> str:
    # Prevent script-tag breakouts when embedding JSON inside inline JS.
    def _sanitize_json_payload(payload: Any) -> Any:
        if payload is None:
            return None
        if isinstance(payload, float):
            return payload if math.isfinite(payload) else None
        if isinstance(payload, (str, int, bool)):
            return payload
        if isinstance(payload, dict):
            return {k: _sanitize_json_payload(v) for k, v in payload.items()}
        if isinstance(payload, (list, tuple, set)):
            return [_sanitize_json_payload(v) for v in payload]
        if isinstance(payload, Path):
            return str(payload)
        # Handle numpy scalar values without introducing a hard runtime dependency here.
        if hasattr(payload, "item"):
            try:
                return _sanitize_json_payload(payload.item())
            except Exception:
                pass
        return payload

    try:
        raw = json.dumps(value, ensure_ascii=True, allow_nan=False, default=str)
    except (TypeError, ValueError):
        cleaned = _sanitize_json_payload(value)
        raw = json.dumps(cleaned, ensure_ascii=True, allow_nan=False, default=str)
    return (
        raw
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def _resolve_base_dir(base_dir: Optional[Path] = None) -> Path:
    return Path(base_dir) if base_dir is not None else Path(__file__).resolve().parent


def _resolve_cst_dataset_keys(headers: Sequence[Any]) -> Dict[str, str]:
    header_list = [str(header or "").strip() for header in list(headers or []) if str(header or "").strip()]
    display_uniprot_key = header_list[0] if header_list else ""
    display_gene_key = header_list[1] if len(header_list) > 1 else ""
    normalized_headers = {header: header for header in header_list}
    ortholog_uniprot_key = normalized_headers.get(HUMAN_ORTHOLOG_UNIPROT_COLUMN, "")
    ortholog_gene_key = normalized_headers.get(HUMAN_ORTHOLOG_GENE_COLUMN, "")
    return {
        "display_uniprot_key": display_uniprot_key,
        "display_gene_key": display_gene_key,
        "match_uniprot_key": ortholog_uniprot_key or display_uniprot_key,
        "match_gene_key": ortholog_gene_key or (display_gene_key if not ortholog_uniprot_key else ""),
    }


def _normalize_forced_pathway_key(file_path: Path) -> str:
    return _canonicalize_label(re.sub(r"\s+\(\d+\)$", "", file_path.stem.strip()))


def _resolve_cst_pathway_dir(base_dir: Optional[Path] = None) -> Path:
    if base_dir is not None:
        candidate = Path(base_dir)
        if candidate.name.lower() == "cst" and candidate.exists():
            return candidate
    project_root = Path(__file__).resolve().parent.parent
    return project_root / "stored_pathways" / "cst"


def _cst_pathway_id_from_name(name: str) -> str:
    return "cst_" + re.sub(r"[^a-z0-9]+", "_", str(name or "").strip().lower()).strip("_")


def _cst_overlay_state_path(file_path: Path) -> Path:
    return file_path.with_suffix(".mapkinase.json")


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


def _normalize_gradient_stops(gradient_stops: Any) -> List[Dict[str, Any]]:
    if not isinstance(gradient_stops, (list, tuple)):
        return []
    parsed: List[Tuple[float, List[int]]] = []
    for item in gradient_stops:
        if not isinstance(item, dict):
            continue
        raw_value = item.get("value")
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(value):
            continue
        raw_color = item.get("color", item.get("rgb"))
        if isinstance(raw_color, str):
            raw_hex = raw_color.strip().lstrip("#")
            if len(raw_hex) != 6:
                continue
            try:
                color = [int(raw_hex[0:2], 16), int(raw_hex[2:4], 16), int(raw_hex[4:6], 16)]
            except Exception:
                continue
        elif isinstance(raw_color, (list, tuple)) and len(raw_color) >= 3:
            color = [int(max(0, min(255, round(float(ch))))) for ch in list(raw_color)[:3]]
        else:
            continue
        parsed.append((float(value), color))
    if not parsed:
        return []
    dedup: Dict[float, List[int]] = {}
    for value, color in parsed:
        dedup[value] = list(color)
    return [{"value": float(v), "color": list(c)} for v, c in sorted(dedup.items(), key=lambda kv: kv[0])]


def _default_gradient_stops(
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
) -> List[Dict[str, Any]]:
    return _normalize_gradient_stops(
        [
            {"value": float(max_negative), "color": _coerce_rgb(negative_color, _DEFAULT_NEGATIVE_COLOR)},
            {"value": 0.0, "color": [255, 255, 255]},
            {"value": float(max_positive), "color": _coerce_rgb(positive_color, _DEFAULT_POSITIVE_COLOR)},
        ]
    )


def _gradient_color_from_stops(fold_value: float, gradient_stops: List[Dict[str, Any]]) -> List[int]:
    if not gradient_stops:
        return [128, 128, 128]
    if fold_value <= float(gradient_stops[0]["value"]):
        return list(gradient_stops[0]["color"])
    if fold_value >= float(gradient_stops[-1]["value"]):
        return list(gradient_stops[-1]["color"])
    for idx in range(len(gradient_stops) - 1):
        left = gradient_stops[idx]
        right = gradient_stops[idx + 1]
        left_v = float(left["value"])
        right_v = float(right["value"])
        if fold_value < left_v or fold_value > right_v:
            continue
        if right_v == left_v:
            return list(right["color"])
        t = (fold_value - left_v) / (right_v - left_v)
        return [
            int(max(0, min(255, round((1 - t) * left["color"][c_idx] + t * right["color"][c_idx]))))
            for c_idx in range(3)
        ]
    return list(gradient_stops[-1]["color"])


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
    gradient_stops = _GRADIENT_STOPS_CTX.get()
    normalized_stops = _normalize_gradient_stops(gradient_stops)
    if len(normalized_stops) >= 2:
        return _gradient_color_from_stops(float(fold), normalized_stops)
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


def _normalize_fc_suffix(header: Any) -> str:
    value = str(header or "").strip()
    if ":" in value:
        value = value.split(":", 1)[1]
    value = re.sub(r"\s+", " ", value.strip())
    return value.lower()


def _outline_column_map(headers: Sequence[Any]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for col in list(headers or []):
        col_text = str(col or "").strip()
        if not col_text.lower().startswith("o:"):
            continue
        key = _normalize_fc_suffix(col_text)
        if key:
            mapping[key] = col_text
    return mapping


def _resolve_outline_columns(main_columns: Sequence[Any], headers: Sequence[Any]) -> List[Optional[str]]:
    outline_map = _outline_column_map(headers)
    return [outline_map.get(_normalize_fc_suffix(col)) for col in list(main_columns or [])]


def _normalize_ptm_shape_name(value: Any) -> str:
    text = str(value or "circle").strip().lower()
    if text == "square":
        return "Square"
    if text == "diamond":
        return "Diamond"
    return "Circle"


def _cst_ptm_symbol_dict(symbol_type: str, reg_site_val: str) -> Optional[Dict[str, Any]]:
    key = str(symbol_type or "").strip().lower()
    if key == "symbol_label_1":
        return {"symbol": "↑", "symbol_icon": "active.svg"}
    if key == "symbol_label_2":
        return {"symbol": "x", "symbol_icon": "inhibit.svg"}
    if key in {"symbol_label_3", "symbol_label_4"}:
        return {"symbol": "+", "symbol_icon": "regsite.svg"}
    if str(reg_site_val or "").strip():
        return {"symbol": "+", "symbol_icon": "regsite.svg"}
    return None


def _classify_cst_ptm_rows(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> List[Dict[str, Any]]:
    header_list = list(headers or [])
    if not header_list:
        return []
    normalized_rows: List[List[Any]] = []
    for row in list(rows or []):
        values = list(row)
        if len(values) < len(header_list):
            values.extend([""] * (len(header_list) - len(values)))
        normalized_rows.append(values[: len(header_list)])
    frame = pd.DataFrame(normalized_rows, columns=header_list)
    if "Phosphosite_Classification" not in frame.columns:
        frame["Phosphosite_Classification"] = "none"
    ptm_symbol_list = json.loads(json.dumps(DEFAULT_DATA["ptm"][0]["ptm_symbol_list"]))
    for item in ptm_symbol_list:
        for key, val in item.items():
            if not isinstance(val, dict):
                continue
            if key == "symbol_label_4_dict":
                val["header_to_search"] = "PSP: regulatory_site"
            else:
                val["header_to_search"] = "PSP: ON_FUNCTION"
    classified = classify_phosphosite_function(
        frame.copy(),
        ptm_symbol_list,
        reg_site_col="PSP: regulatory_site",
        reg_function_col="PSP: ON_FUNCTION",
        output_col="Phosphosite_Classification",
    )
    return classified.fillna("").to_dict(orient="records")


def _build_cst_ptm_index(
    ptm_dataset: Optional[Dict[str, Any]],
    ptm_settings: Optional[Dict[str, Any]],
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
) -> Dict[str, Any]:
    settings = dict(ptm_settings or {})
    try:
        max_display = max(0, int(settings.get("ptm_max_display", 4) or 4))
    except (TypeError, ValueError):
        max_display = 4
    if not ptm_dataset or max_display <= 0:
        return {"enabled": False, "max_display": max_display, "items_by_uniprot": {}}
    headers = list(ptm_dataset.get("headers") or [])
    rows = list(ptm_dataset.get("rows") or [])
    if len(headers) < 2:
        return {"enabled": False, "max_display": max_display, "items_by_uniprot": {}}
    dataset_keys = _resolve_cst_dataset_keys(headers)
    classified_rows = _classify_cst_ptm_rows(headers, rows)
    main_headers = _normalize_fc_headers(
        ptm_dataset.get("main_columns") or [header for header in headers if str(header).startswith("C:")]
    )
    outline_headers = list(ptm_dataset.get("outline_columns") or _resolve_outline_columns(main_headers, headers))
    tooltip_columns = _infer_ptm_tooltip_columns(ptm_dataset)
    uniprot_key = str(dataset_keys.get("match_uniprot_key") or (headers[0] if headers else "")).strip()
    display_uniprot_key = str(dataset_keys.get("display_uniprot_key") or (headers[0] if headers else "")).strip()
    site_key = headers[1]
    label_font = str(settings.get("ptm_label_font") or "Arial")
    label_color = _coerce_rgb(settings.get("ptm_label_color"), (0, 0, 0))
    try:
        label_size = max(6.0, float(settings.get("ptm_label_size", 7.0) or 7.0))
    except (TypeError, ValueError):
        label_size = 7.0
    try:
        circle_radius = max(3.0, float(settings.get("ptm_circle_radius", 4.0) or 4.0))
    except (TypeError, ValueError):
        circle_radius = 4.0
    try:
        outline_width = max(0.6, float(settings.get("ptm_outline_width", 1.0) or 1.0))
    except (TypeError, ValueError):
        outline_width = 1.0
    shape_name = _normalize_ptm_shape_name(settings.get("ptm_shape"))
    items_by_uniprot: Dict[str, List[Dict[str, Any]]] = {}
    def _normalize_ptm_uniprot(value: Any) -> str:
        token = str(value or "").strip().upper()
        if not token:
            return ""
        token = re.split(r"[|,;/]+", token, maxsplit=1)[0].strip()
        return token.split("-", 1)[0].strip() if token else ""

    for row_idx, row_map in enumerate(classified_rows):
        normalized_uniprot = _normalize_ptm_uniprot(row_map.get(uniprot_key, ""))
        if not normalized_uniprot:
            normalized_uniprot = _normalize_ptm_uniprot(row_map.get(display_uniprot_key, ""))
        site_label = str(row_map.get(site_key, "") or "").strip()
        if not normalized_uniprot or not site_label:
            continue
        reg_site_val = str(row_map.get("PSP: regulatory_site", "") or "").strip()
        symbol_class = str(row_map.get("Phosphosite_Classification", "") or "").strip()
        symbol_dict = _cst_ptm_symbol_dict(symbol_class, reg_site_val)
        symbol_text = str(symbol_dict.get("symbol") or "") if symbol_dict else ""
        symbol_icon = _resolve_symbol_icon(symbol_dict) if symbol_dict else None
        entry: Dict[str, Any] = {
            "site_label": site_label,
            "label": site_label,
            "symbol": symbol_text,
            "symbol_type": symbol_class,
            "symbol_icon": symbol_icon,
            "shape": shape_name,
            "radius": circle_radius,
            "outline_width": outline_width,
            "label_font": label_font,
            "label_color": list(label_color),
            "label_size": label_size,
            "symbol_font": "Arial",
            "symbol_color": [0, 0, 0],
            "symbol_size": max(6.0, min(11.0, (circle_radius * 1.05) + 1.0)),
            "sort_order": row_idx,
        }
        entry.update(_build_dataset_tooltip(row_map, tooltip_columns))
        for idx, header in enumerate(main_headers, 1):
            raw_value = row_map.get(header, "")
            fold_value = _coerce_float_or_none(raw_value)
            entry[f"fold_change_{idx}"] = fold_value
            entry[f"fc_color_{idx}"] = _gradient_color_from_fold(
                fold_value,
                negative_color,
                positive_color,
                max_negative,
                max_positive,
            )
            outline_header = outline_headers[idx - 1] if idx - 1 < len(outline_headers) else None
            outline_value = _coerce_float_or_none(row_map.get(outline_header, "")) if outline_header else None
            entry[f"outline_fold_change_{idx}"] = outline_value
            entry[f"outline_color_{idx}"] = (
                _gradient_color_from_fold(
                    outline_value,
                    negative_color,
                    positive_color,
                    max_negative,
                    max_positive,
                )
                if outline_value is not None
                else [0, 0, 0]
            )
        items_by_uniprot.setdefault(normalized_uniprot, []).append(entry)
    for uniprot_id, entries in list(items_by_uniprot.items()):
        entries.sort(key=lambda item: (int(item.get("sort_order") or 0), str(item.get("site_label") or "").upper()))
        if max_display > 0:
            items_by_uniprot[uniprot_id] = entries[:max_display]
    return {
        "enabled": bool(items_by_uniprot),
        "max_display": max_display,
        "shape": shape_name,
        "radius": circle_radius,
        "items_by_uniprot": items_by_uniprot,
    }


def _is_regulatory_cst_ptm_entry(entry: Optional[Dict[str, Any]]) -> bool:
    if not isinstance(entry, dict):
        return False
    symbol_type = str(entry.get("symbol_type") or "").strip().lower()
    if symbol_type in {"symbol_label_3", "symbol_label_4"}:
        return True
    symbol_icon = str(entry.get("symbol_icon") or "").strip().lower()
    if symbol_icon.endswith("regsite.svg"):
        return True
    return str(entry.get("symbol") or "").strip() == "+"


def _best_cst_regulatory_ptm_fold(
    ptm_index: Optional[Dict[str, Any]],
    uniprot_id: Any,
    fc_idx: int,
) -> Optional[float]:
    if not isinstance(ptm_index, dict) or not bool(ptm_index.get("enabled")):
        return None
    key = str(uniprot_id or "").strip().upper().split("-", 1)[0]
    if not key:
        return None
    entries = list((ptm_index.get("items_by_uniprot") or {}).get(key, []) or [])
    best_value: Optional[float] = None
    for entry in entries:
        if not _is_regulatory_cst_ptm_entry(entry):
            continue
        value = _coerce_float_or_none(entry.get(f"fold_change_{fc_idx}", ""))
        if value is None:
            continue
        if best_value is None or abs(value) > abs(best_value):
            best_value = value
    return best_value


def _choose_cst_display_row(
    matches: Sequence[Dict[str, Any]],
    headers: Sequence[Any],
    header: Any,
    ptm_index: Optional[Dict[str, Any]],
    fc_idx: int,
) -> tuple[Optional[Dict[str, Any]], Optional[float]]:
    header_text = str(header or "").strip()
    if not matches or not header_text:
        return None, None
    uniprot_key = str(headers[0] or "").strip() if headers else ""
    reg_row: Optional[Dict[str, Any]] = None
    reg_value: Optional[float] = None
    protein_row: Optional[Dict[str, Any]] = None
    protein_value: Optional[float] = None
    for row in matches:
        protein_fc = _coerce_float_or_none(row.get(header_text, ""))
        if protein_fc is not None and (protein_value is None or abs(protein_fc) > abs(protein_value)):
            protein_value = protein_fc
            protein_row = row
        regulatory_fc = _best_cst_regulatory_ptm_fold(
            ptm_index,
            row.get(uniprot_key, "") if uniprot_key else "",
            fc_idx,
        )
        if regulatory_fc is None:
            continue
        if reg_value is None or abs(regulatory_fc) > abs(reg_value):
            reg_value = regulatory_fc
            reg_row = row
    if reg_row is not None:
        chosen_value = _coerce_float_or_none(reg_row.get(header_text, ""))
        return reg_row, chosen_value if chosen_value is not None else reg_value
    return protein_row, protein_value


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
                    byte_value = int(octal, 8) & 0xFF
                    # PDF literal octal escapes are byte values (often WinAnsi encoded).
                    output.append(bytes([byte_value]).decode("cp1252", errors="replace"))
                except ValueError:
                    pass
            else:
                output.append(char)
        else:
            output.append(char)
        idx += 1
    return "".join(output)


def _normalize_pdf_text_encoding(value: Any) -> str:
    text = str(value or "")
    if text:
        normalized_chars: List[str] = []
        for ch in text:
            code = ord(ch)
            if 0x80 <= code <= 0x9F:
                normalized_chars.append(bytes([code]).decode("cp1252", errors="replace"))
            else:
                normalized_chars.append(ch)
        text = "".join(normalized_chars)
    return text.replace("\ufffd", "")


def _clean_text_label(value: str) -> str:
    text = _normalize_pdf_text_encoding(value)
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
    text_blocks: List[str] = []
    for stream in _iter_decompressed_streams(file_path):
        if "BT" not in stream or ("Tj" not in stream and "TJ" not in stream):
            continue
        for block in re.findall(r"BT(.*?)ET", stream, re.S):
            if "Tj" in block or "TJ" in block:
                text_blocks.append(block)
    if not text_blocks:
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
    for content in text_blocks:
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
    contents = [
        stream
        for stream in _iter_decompressed_streams(file_path)
        if " cm 0 0 m " in stream and " c " in stream
    ]
    if not contents:
        return []

    groups: List[Dict[str, float]] = []
    for content in contents:
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
def _extract_cst_rect_groups(file_path_str: str) -> List[Dict[str, float]]:
    file_path = Path(file_path_str)
    contents = [
        stream
        for stream in _iter_decompressed_streams(file_path)
        if " re " in stream or "\nre " in stream
    ]
    if not contents:
        return []

    groups: List[Dict[str, float]] = []
    for content in contents:
        for match in _RECT_BLOCK_RE.finditer(content):
            try:
                raw_x = float(match.group(1))
                raw_y = float(match.group(2))
                raw_w = float(match.group(3))
                raw_h = float(match.group(4))
            except Exception:
                continue
            width = abs(raw_w)
            height = abs(raw_h)
            if width < 14.0 or height < 7.0:
                continue
            if width > 120.0 or height > 50.0:
                continue
            left = min(raw_x, raw_x + raw_w)
            right = max(raw_x, raw_x + raw_w)
            bottom = min(raw_y, raw_y + raw_h)
            top = max(raw_y, raw_y + raw_h)
            center_x = (left + right) * 0.5
            center_y = (bottom + top) * 0.5
            merged = False
            for group in groups:
                if abs(group["center_x"] - center_x) <= 2.5 and abs(group["center_y"] - center_y) <= 2.5:
                    group["left"] = min(group["left"], left)
                    group["right"] = max(group["right"], right)
                    group["bottom"] = min(group["bottom"], bottom)
                    group["top"] = max(group["top"], top)
                    group["center_x"] = (group["left"] + group["right"]) * 0.5
                    group["center_y"] = (group["bottom"] + group["top"]) * 0.5
                    group["radius_x"] = (group["right"] - group["left"]) * 0.5
                    group["radius_y"] = (group["top"] - group["bottom"]) * 0.5
                    merged = True
                    break
            if not merged:
                groups.append(
                    {
                        "left": left,
                        "right": right,
                        "bottom": bottom,
                        "top": top,
                        "center_x": center_x,
                        "center_y": center_y,
                        "radius_x": width * 0.5,
                        "radius_y": height * 0.5,
                    }
                )
    return groups


def _matrix_multiply(
    left: Sequence[float], right: Sequence[float]
) -> tuple[float, float, float, float, float, float]:
    a1, b1, c1, d1, e1, f1 = [float(x) for x in left]
    a2, b2, c2, d2, e2, f2 = [float(x) for x in right]
    return (
        (a1 * a2) + (c1 * b2),
        (b1 * a2) + (d1 * b2),
        (a1 * c2) + (c1 * d2),
        (b1 * c2) + (d1 * d2),
        (a1 * e2) + (c1 * f2) + e1,
        (b1 * e2) + (d1 * f2) + f1,
    )


def _transform_point(matrix: Sequence[float], x: float, y: float) -> tuple[float, float]:
    a, b, c, d, e, f = [float(v) for v in matrix]
    return ((a * float(x)) + (c * float(y)) + e, (b * float(x)) + (d * float(y)) + f)


def _path_bbox(commands: Sequence[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    xs: List[float] = []
    ys: List[float] = []
    for command in list(commands or []):
        points = list(command.get("points") or [])
        for point in points:
            try:
                x, y = point
            except Exception:
                continue
            xs.append(float(x))
            ys.append(float(y))
    if not xs or not ys:
        return None
    return {
        "left": min(xs),
        "right": max(xs),
        "bottom": min(ys),
        "top": max(ys),
        "center_x": (min(xs) + max(xs)) * 0.5,
        "center_y": (min(ys) + max(ys)) * 0.5,
        "width": max(xs) - min(xs),
        "height": max(ys) - min(ys),
    }


def _path_end_tangent(commands: Sequence[Dict[str, Any]]) -> tuple[float, float]:
    current: Optional[tuple[float, float]] = None
    tangent = (1.0, 0.0)
    for command in list(commands or []):
        op = str(command.get("op") or "")
        points = list(command.get("points") or [])
        if op == "m":
            if points:
                current = (float(points[0][0]), float(points[0][1]))
            continue
        if op == "l" and current and points:
            end = (float(points[0][0]), float(points[0][1]))
            tangent = (end[0] - current[0], end[1] - current[1])
            current = end
            continue
        if op == "c" and current and len(points) >= 3:
            c2 = (float(points[1][0]), float(points[1][1]))
            end = (float(points[2][0]), float(points[2][1]))
            tangent = (end[0] - c2[0], end[1] - c2[1])
            current = end
            continue
    length = math.hypot(tangent[0], tangent[1]) or 1.0
    return (tangent[0] / length, tangent[1] / length)


def _sample_cubic(
    p0: tuple[float, float],
    c1: tuple[float, float],
    c2: tuple[float, float],
    p3: tuple[float, float],
    steps: int = 10,
) -> List[tuple[float, float]]:
    output: List[tuple[float, float]] = []
    for idx in range(1, max(2, steps) + 1):
        t = idx / max(2, steps)
        mt = 1.0 - t
        x = (mt**3 * p0[0]) + (3 * mt * mt * t * c1[0]) + (3 * mt * t * t * c2[0]) + (t**3 * p3[0])
        y = (mt**3 * p0[1]) + (3 * mt * mt * t * c1[1]) + (3 * mt * t * t * c2[1]) + (t**3 * p3[1])
        output.append((x, y))
    return output


def _path_polyline(commands: Sequence[Dict[str, Any]]) -> List[tuple[float, float]]:
    points: List[tuple[float, float]] = []
    current: Optional[tuple[float, float]] = None
    subpath_start: Optional[tuple[float, float]] = None
    for command in list(commands or []):
        op = str(command.get("op") or "")
        cmd_points = list(command.get("points") or [])
        if op == "m" and cmd_points:
            current = (float(cmd_points[0][0]), float(cmd_points[0][1]))
            subpath_start = current
            points.append(current)
        elif op == "l" and current and cmd_points:
            current = (float(cmd_points[0][0]), float(cmd_points[0][1]))
            points.append(current)
        elif op == "c" and current and len(cmd_points) >= 3:
            c1 = (float(cmd_points[0][0]), float(cmd_points[0][1]))
            c2 = (float(cmd_points[1][0]), float(cmd_points[1][1]))
            end = (float(cmd_points[2][0]), float(cmd_points[2][1]))
            points.extend(_sample_cubic(current, c1, c2, end))
            current = end
        elif op == "h" and current and subpath_start:
            current = subpath_start
            points.append(current)
    return points


def _polyline_length(points: Sequence[tuple[float, float]]) -> float:
    if len(points) < 2:
        return 0.0
    total = 0.0
    for idx in range(1, len(points)):
        total += math.hypot(points[idx][0] - points[idx - 1][0], points[idx][1] - points[idx - 1][1])
    return total


def _polyline_midpoint(points: Sequence[tuple[float, float]]) -> Optional[tuple[float, float]]:
    if len(points) < 2:
        return points[0] if points else None
    total = _polyline_length(points)
    if total <= 0:
        return points[len(points) // 2]
    target = total * 0.5
    traversed = 0.0
    for idx in range(1, len(points)):
        segment = math.hypot(points[idx][0] - points[idx - 1][0], points[idx][1] - points[idx - 1][1])
        if traversed + segment >= target and segment > 0:
            frac = (target - traversed) / segment
            x = points[idx - 1][0] + ((points[idx][0] - points[idx - 1][0]) * frac)
            y = points[idx - 1][1] + ((points[idx][1] - points[idx - 1][1]) * frac)
            return (x, y)
        traversed += segment
    return points[-1]


def _quadratic_control_from_midpoint(
    start: tuple[float, float], mid: tuple[float, float], end: tuple[float, float]
) -> tuple[float, float]:
    return (
        (2.0 * float(mid[0])) - (0.5 * (float(start[0]) + float(end[0]))),
        (2.0 * float(mid[1])) - (0.5 * (float(start[1]) + float(end[1]))),
    )


def _extract_cst_path_objects(file_path_str: str) -> List[Dict[str, Any]]:
    file_path = Path(file_path_str)
    shapes: List[Dict[str, Any]] = []
    for stream in _iter_decompressed_streams(file_path):
        if " m" not in stream and "\nm" not in stream:
            continue
        tokens = list(_PDF_TOKEN_RE.findall(stream))
        state_stack: List[Dict[str, Any]] = []
        state: Dict[str, Any] = {
            "matrix": (1.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            "dash": [],
            "line_width": 1.0,
        }
        path_commands: List[Dict[str, Any]] = []
        number_stack: List[float] = []
        array_token: Optional[str] = None
        for token in tokens:
            if token.startswith("[") and token.endswith("]"):
                array_token = token
                continue
            try:
                number_stack.append(float(token))
                continue
            except ValueError:
                pass
            op = token
            if op == "q":
                state_stack.append(
                    {
                        "matrix": tuple(state["matrix"]),
                        "dash": list(state["dash"]),
                        "line_width": float(state["line_width"]),
                    }
                )
                continue
            if op == "Q":
                if state_stack:
                    state = state_stack.pop()
                path_commands = []
                number_stack = []
                array_token = None
                continue
            if op == "cm" and len(number_stack) >= 6:
                vals = tuple(number_stack[-6:])
                number_stack = number_stack[:-6]
                state["matrix"] = _matrix_multiply(tuple(state["matrix"]), vals)
                continue
            if op == "w" and number_stack:
                state["line_width"] = float(number_stack[-1])
                number_stack = number_stack[:-1]
                continue
            if op == "d":
                dash: List[float] = []
                for raw_dash in re.findall(r"-?(?:\d+(?:\.\d+)?|\.\d+)", array_token or ""):
                    try:
                        dash.append(float(raw_dash))
                    except ValueError:
                        continue
                state["dash"] = [x for x in dash if x > 0]
                if number_stack:
                    number_stack = number_stack[:-1]
                array_token = None
                continue
            if op == "m" and len(number_stack) >= 2:
                x, y = number_stack[-2], number_stack[-1]
                number_stack = number_stack[:-2]
                path_commands.append({"op": "m", "points": [_transform_point(state["matrix"], x, y)]})
                continue
            if op == "l" and len(number_stack) >= 2:
                x, y = number_stack[-2], number_stack[-1]
                number_stack = number_stack[:-2]
                path_commands.append({"op": "l", "points": [_transform_point(state["matrix"], x, y)]})
                continue
            if op == "c" and len(number_stack) >= 6:
                x1, y1, x2, y2, x3, y3 = number_stack[-6:]
                number_stack = number_stack[:-6]
                path_commands.append(
                    {
                        "op": "c",
                        "points": [
                            _transform_point(state["matrix"], x1, y1),
                            _transform_point(state["matrix"], x2, y2),
                            _transform_point(state["matrix"], x3, y3),
                        ],
                    }
                )
                continue
            if op == "h":
                path_commands.append({"op": "h", "points": []})
                continue
            if op == "re" and len(number_stack) >= 4:
                x, y, w, h = number_stack[-4:]
                number_stack = number_stack[:-4]
                p1 = _transform_point(state["matrix"], x, y)
                p2 = _transform_point(state["matrix"], x + w, y)
                p3 = _transform_point(state["matrix"], x + w, y + h)
                p4 = _transform_point(state["matrix"], x, y + h)
                path_commands.extend(
                    [
                        {"op": "m", "points": [p1]},
                        {"op": "l", "points": [p2]},
                        {"op": "l", "points": [p3]},
                        {"op": "l", "points": [p4]},
                        {"op": "h", "points": []},
                    ]
                )
                continue
            if op in {"S", "s", "f", "F", "B", "b", "B*", "b*", "f*"}:
                if op in {"s", "b", "b*"}:
                    path_commands.append({"op": "h", "points": []})
                if path_commands:
                    shapes.append(
                        {
                            "paint": "stroke" if op in {"S", "s"} else ("fill" if op in {"f", "F", "f*"} else "fillstroke"),
                            "commands": json.loads(json.dumps(path_commands)),
                            "dash": list(state["dash"]),
                            "line_width": float(state["line_width"]),
                            "closed": any(str(cmd.get("op") or "") == "h" for cmd in path_commands),
                        }
                    )
                path_commands = []
                number_stack = []
                array_token = None
                continue
            if op in {"n", "W", "W*"}:
                if op == "n":
                    path_commands = []
                number_stack = []
                continue
            array_token = None
    return shapes


@lru_cache(maxsize=16)
def _extract_cst_edge_groups(file_path_str: str) -> List[Dict[str, Any]]:
    page_size = _extract_cst_page_size(file_path_str)
    page_height = float(page_size.get("page_height") or _CST_PAGE_HEIGHT)
    shapes = _extract_cst_path_objects(file_path_str)
    shafts: List[Dict[str, Any]] = []
    arrowheads: List[Dict[str, Any]] = []
    inhibitor_bars: List[Dict[str, Any]] = []
    for shape in shapes:
        commands = list(shape.get("commands") or [])
        bbox = _path_bbox(commands)
        polyline = _path_polyline(commands)
        if not bbox or not polyline:
            continue
        length = _polyline_length(polyline)
        paint = str(shape.get("paint") or "")
        closed = bool(shape.get("closed"))
        line_width = float(shape.get("line_width") or 1.0)
        dash = list(shape.get("dash") or [])
        if paint in {"fill", "fillstroke"} and closed and bbox["width"] <= 16.0 and bbox["height"] <= 16.0:
            aspect = (bbox["width"] / max(bbox["height"], 0.001)) if bbox["height"] > 0 else 999.0
            if aspect >= 2.25 or aspect <= 0.44:
                inhibitor_tangent = (1.0, 0.0) if aspect >= 1.0 else (0.0, 1.0)
                inhibitor_bars.append(
                    {
                        "start": (bbox["left"], bbox["bottom"]),
                        "end": (bbox["right"], bbox["top"]),
                        "mid_x": bbox["center_x"],
                        "mid_y": bbox["center_y"],
                        "tangent": inhibitor_tangent,
                        "length": max(bbox["width"], bbox["height"]),
                    }
                )
            else:
                arrowheads.append(
                    {
                        "center_x": bbox["center_x"],
                        "center_y": bbox["center_y"],
                        "bbox": bbox,
                    }
                )
            continue
        if paint == "stroke" and not closed and length >= 6.0 and bbox["width"] <= 18.0 and bbox["height"] <= 18.0:
            start = polyline[0]
            end = polyline[-1]
            mid = ((start[0] + end[0]) * 0.5, (start[1] + end[1]) * 0.5)
            tangent = _path_end_tangent(commands)
            inhibitor_bars.append(
                {
                    "start": start,
                    "end": end,
                    "mid_x": mid[0],
                    "mid_y": mid[1],
                    "tangent": tangent,
                    "length": length,
                }
            )
            continue
        if paint == "stroke" and not closed and length >= 12.0:
            start = polyline[0]
            end = polyline[-1]
            mid = _polyline_midpoint(polyline) or ((start[0] + end[0]) * 0.5, (start[1] + end[1]) * 0.5)
            control_x, control_y = _quadratic_control_from_midpoint(start, mid, end)
            shafts.append(
                {
                    "start_x": start[0],
                    "start_y": start[1],
                    "end_x": end[0],
                    "end_y": end[1],
                    "control_x": control_x,
                    "control_y": control_y,
                    "dash": dash,
                    "line_width": line_width,
                    "tangent": _path_end_tangent(commands),
                }
            )
    output: List[Dict[str, Any]] = []
    used_arrowheads: set[int] = set()
    used_bars: set[int] = set()
    for shaft in shafts:
        end_x = float(shaft["end_x"])
        end_y = float(shaft["end_y"])
        tangent = shaft["tangent"]
        best_arrow_idx = -1
        best_arrow_score = 999999.0
        for idx, head in enumerate(arrowheads):
            if idx in used_arrowheads:
                continue
            dist = math.hypot(float(head["center_x"]) - end_x, float(head["center_y"]) - end_y)
            if dist > 16.0:
                continue
            score = dist
            if score < best_arrow_score:
                best_arrow_score = score
                best_arrow_idx = idx
        best_bar_idx = -1
        best_bar_score = 999999.0
        for idx, bar in enumerate(inhibitor_bars):
            if idx in used_bars:
                continue
            dist = math.hypot(float(bar["mid_x"]) - end_x, float(bar["mid_y"]) - end_y)
            if dist > 14.0:
                continue
            dot = abs((float(bar["tangent"][0]) * float(tangent[0])) + (float(bar["tangent"][1]) * float(tangent[1])))
            if dot > 0.45:
                continue
            score = dist + (dot * 10.0)
            if score < best_bar_score:
                best_bar_score = score
                best_bar_idx = idx
        if best_arrow_idx < 0 and best_bar_idx < 0:
            continue
        end_type = "arrow"
        if best_bar_idx >= 0 and (best_arrow_idx < 0 or best_bar_score <= best_arrow_score):
            end_type = "inhibitor"
            used_bars.add(best_bar_idx)
        elif best_arrow_idx >= 0:
            used_arrowheads.add(best_arrow_idx)
        output.append(
            {
                "id": f"cst_ai_edge_{len(output) + 1}",
                "startX": float(shaft["start_x"]),
                "startY": page_height - float(shaft["start_y"]),
                "endX": end_x,
                "endY": page_height - end_y,
                "controlX": float(shaft["control_x"]),
                "controlY": page_height - float(shaft["control_y"]),
                "stroke": "#0f172a",
                "strokeWidth": 1.6,
                "opacity": 0.95,
                "userCreated": False,
                "dashed": any(float(x) > 0 for x in list(shaft["dash"] or [])),
                "endType": end_type,
            }
        )
    if len(output) >= 8:
        return output
    if not _should_run_cst_image_edge_fallback(file_path_str, shapes, shafts, output):
        return output
    image_edges = _extract_cst_edge_groups_image(file_path_str)
    if not image_edges:
        return output
    merged = list(output)
    for edge in image_edges:
        if len(merged) >= 14:
            break
        if _edge_matches_existing(edge, merged):
            continue
        merged.append(
            {
                **edge,
                "id": f"cst_ai_edge_{len(merged) + 1}",
            }
        )
    return merged


def _should_run_cst_image_edge_fallback(
    file_path_str: str,
    shapes: Sequence[Dict[str, Any]],
    shafts: Sequence[Dict[str, Any]],
    vector_edges: Sequence[Dict[str, Any]],
) -> bool:
    file_path = Path(file_path_str)
    try:
        file_size_mb = file_path.stat().st_size / (1024 * 1024)
    except Exception:
        file_size_mb = 0.0
    shape_count = len(list(shapes or []))
    shaft_count = len(list(shafts or []))
    edge_count = len(list(vector_edges or []))
    if file_size_mb >= 2.2:
        return False
    if shape_count >= 1800:
        return False
    if shaft_count >= 260:
        return False
    if edge_count >= 4 and shape_count >= 900:
        return False
    return True


def _edge_matches_existing(edge: Dict[str, Any], existing: Sequence[Dict[str, Any]]) -> bool:
    start = (float(edge.get("startX") or 0.0), float(edge.get("startY") or 0.0))
    end = (float(edge.get("endX") or 0.0), float(edge.get("endY") or 0.0))
    mid = (
        (float(edge.get("startX") or 0.0) + float(edge.get("endX") or 0.0)) * 0.5,
        (float(edge.get("startY") or 0.0) + float(edge.get("endY") or 0.0)) * 0.5,
    )
    for other in list(existing or []):
        other_start = (float(other.get("startX") or 0.0), float(other.get("startY") or 0.0))
        other_end = (float(other.get("endX") or 0.0), float(other.get("endY") or 0.0))
        other_mid = (
            (float(other.get("startX") or 0.0) + float(other.get("endX") or 0.0)) * 0.5,
            (float(other.get("startY") or 0.0) + float(other.get("endY") or 0.0)) * 0.5,
        )
        direct = math.hypot(start[0] - other_start[0], start[1] - other_start[1]) + math.hypot(end[0] - other_end[0], end[1] - other_end[1])
        swapped = math.hypot(start[0] - other_end[0], start[1] - other_end[1]) + math.hypot(end[0] - other_start[0], end[1] - other_start[1])
        if min(direct, swapped) <= 20.0 and math.hypot(mid[0] - other_mid[0], mid[1] - other_mid[1]) <= 10.0:
            return True
    return False


@lru_cache(maxsize=16)
def _extract_cst_edge_groups_image(file_path_str: str) -> List[Dict[str, Any]]:
    try:
        import cv2  # type: ignore
        import numpy as np  # type: ignore
        import pypdfium2 as pdfium  # type: ignore
    except Exception:
        return []
    file_path = Path(file_path_str)
    page_size = _extract_cst_page_size(file_path_str)
    page_width = float(page_size.get("page_width") or _CST_PAGE_WIDTH)
    page_height = float(page_size.get("page_height") or _CST_PAGE_HEIGHT)
    try:
        pdf = pdfium.PdfDocument(str(file_path))
        page = pdf[0]
        image = page.render(scale=2).to_numpy()
    except Exception:
        return []
    if image is None or not hasattr(image, "shape") or len(image.shape) < 2:
        return []
    img_h, img_w = int(image.shape[0]), int(image.shape[1])
    if img_h <= 0 or img_w <= 0:
        return []
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY if image.shape[2] == 3 else cv2.COLOR_BGRA2GRAY)
    _, mask = cv2.threshold(gray, 185, 255, cv2.THRESH_BINARY_INV)

    for text_node in _extract_cst_text_nodes(file_path_str):
        x = float(text_node.get("pdf_x") or 0.0)
        y = float(text_node.get("pdf_y") or 0.0)
        tw = float(text_node.get("estimated_width") or 0.0)
        th = float(text_node.get("estimated_height") or 0.0)
        x1 = max(0, int(((x - 8.0) / page_width) * img_w))
        x2 = min(img_w, int(((x + tw + 8.0) / page_width) * img_w))
        y1 = max(0, int(((page_height - (y + th) - 8.0) / page_height) * img_h))
        y2 = min(img_h, int(((page_height - (y - (th * 0.2)) + 8.0) / page_height) * img_h))
        if x2 > x1 and y2 > y1:
            cv2.rectangle(mask, (x1, y1), (x2, y2), 0, -1)

    for ellipse in _extract_cst_ellipse_groups(file_path_str):
        cx = int((float(ellipse.get("center_x") or 0.0) / page_width) * img_w)
        cy = int(((page_height - float(ellipse.get("center_y") or 0.0)) / page_height) * img_h)
        rx = max(1, int(((float(ellipse.get("radius_x") or 0.0) + 6.0) / page_width) * img_w))
        ry = max(1, int(((float(ellipse.get("radius_y") or 0.0) + 6.0) / page_height) * img_h))
        cv2.ellipse(mask, (cx, cy), (rx, ry), 0, 0, 360, 0, -1)

    for rect in _extract_cst_rect_groups(file_path_str):
        x1 = max(0, int(((float(rect.get("left") or 0.0) - 6.0) / page_width) * img_w))
        x2 = min(img_w, int(((float(rect.get("right") or 0.0) + 6.0) / page_width) * img_w))
        y1 = max(0, int(((page_height - float(rect.get("top") or 0.0) - 6.0) / page_height) * img_h))
        y2 = min(img_h, int(((page_height - float(rect.get("bottom") or 0.0) + 6.0) / page_height) * img_h))
        if x2 > x1 and y2 > y1:
            cv2.rectangle(mask, (x1, y1), (x2, y2), 0, -1)

    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((2, 2), np.uint8))
    lines = cv2.HoughLinesP(
        mask,
        1,
        np.pi / 180.0,
        threshold=42,
        minLineLength=max(160, int(img_w * 0.12)),
        maxLineGap=max(8, int(img_w * 0.005)),
    )
    if lines is None:
        return []

    def _sample_line_binary(x1: int, y1: int, x2: int, y2: int) -> List[int]:
        steps = max(12, int(math.hypot(x2 - x1, y2 - y1) / 3.0))
        values: List[int] = []
        for idx in range(steps + 1):
            t = idx / max(1, steps)
            x = int(round(x1 + ((x2 - x1) * t)))
            y = int(round(y1 + ((y2 - y1) * t)))
            if 0 <= x < img_w and 0 <= y < img_h:
                values.append(1 if mask[y, x] > 0 else 0)
        return values

    def _detect_dashed(values: Sequence[int]) -> bool:
        if not values:
            return False
        transitions = 0
        for idx in range(1, len(values)):
            if values[idx] != values[idx - 1]:
                transitions += 1
        filled_ratio = sum(values) / max(1, len(values))
        return transitions >= 4 and 0.18 <= filled_ratio <= 0.82

    def _detect_inhibitor_for_segment(x1: int, y1: int, x2: int, y2: int) -> bool:
        angle = math.atan2(y2 - y1, x2 - x1)
        length = math.hypot(x2 - x1, y2 - y1)
        if length < 12:
            return False
        for ex, ey in ((x1, y1), (x2, y2)):
            rx1 = max(0, int(ex - 14))
            rx2 = min(img_w, int(ex + 14))
            ry1 = max(0, int(ey - 14))
            ry2 = min(img_h, int(ey + 14))
            roi = mask[ry1:ry2, rx1:rx2]
            if roi.size == 0:
                continue
            roi_lines = cv2.HoughLinesP(roi, 1, np.pi / 180.0, threshold=8, minLineLength=5, maxLineGap=3)
            if roi_lines is None:
                continue
            for seg in roi_lines[:, 0, :]:
                sx1, sy1, sx2, sy2 = [int(v) for v in seg]
                seg_len = math.hypot(sx2 - sx1, sy2 - sy1)
                if seg_len < 5 or seg_len > 20:
                    continue
                seg_angle = math.atan2(sy2 - sy1, sx2 - sx1)
                angle_delta = abs(((seg_angle - angle + (math.pi * 0.5)) % math.pi) - (math.pi * 0.5))
                if angle_delta > 0.35:
                    continue
                mid_x = rx1 + ((sx1 + sx2) * 0.5)
                mid_y = ry1 + ((sy1 + sy2) * 0.5)
                if math.hypot(mid_x - ex, mid_y - ey) <= 10:
                    return True
        return False

    output: List[Dict[str, Any]] = []
    seen_bins: set[tuple[int, int, int, int]] = set()
    for line in lines[:, 0, :]:
        x1, y1, x2, y2 = [int(v) for v in line]
        length = math.hypot(x2 - x1, y2 - y1)
        if length < max(160, img_w * 0.12):
            continue
        key = tuple(sorted((int(round(x1 / 12)), int(round(y1 / 12)), int(round(x2 / 12)), int(round(y2 / 12)))))
        if key in seen_bins:
            continue
        seen_bins.add(key)
        dashed = _detect_dashed(_sample_line_binary(x1, y1, x2, y2))
        end_type = "inhibitor" if _detect_inhibitor_for_segment(x1, y1, x2, y2) else "arrow"
        start_x = (x1 / img_w) * page_width
        end_x = (x2 / img_w) * page_width
        start_y = page_height - ((y1 / img_h) * page_height)
        end_y = page_height - ((y2 / img_h) * page_height)
        control_x = (start_x + end_x) * 0.5
        control_y = (start_y + end_y) * 0.5
        output.append(
            {
                "id": f"cst_img_edge_{len(output) + 1}",
                "startX": float(start_x),
                "startY": float(start_y),
                "endX": float(end_x),
                "endY": float(end_y),
                "controlX": float(control_x),
                "controlY": float(control_y),
                "stroke": "#0f172a",
                "strokeWidth": 1.6,
                "opacity": 0.95,
                "dashed": dashed,
                "endType": end_type,
                "_length": length,
            }
        )
    output.sort(key=lambda item: float(item.get("_length") or 0.0), reverse=True)
    trimmed: List[Dict[str, Any]] = []
    for edge in output[:24]:
        clean_edge = dict(edge)
        clean_edge.pop("_length", None)
        trimmed.append(clean_edge)
    return trimmed


@lru_cache(maxsize=16)
def _get_cst_module_index(file_path_str: str) -> Dict[str, Dict[str, Any]]:
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
    return module_index


@lru_cache(maxsize=1)
def _get_cst_label_mapper() -> PathwayLabelMapper:
    return PathwayLabelMapper(use_uniprot_rest=False)


@lru_cache(maxsize=4096)
def _resolve_cst_label_mapping_cached(label: str, file_path_str: str) -> Dict[str, Any]:
    file_path = Path(file_path_str)
    module_index = _get_cst_module_index(str(file_path))
    indexed = None
    indexed_variant = ""
    for variant in _cst_label_variants(str(label or "")):
        norm = _clean_text_label(variant).upper()
        candidate = module_index.get(norm)
        if candidate:
            indexed = candidate
            indexed_variant = variant
            break
    if indexed:
        return {
            "original_label": label,
            "normalized_label": _clean_text_label(indexed_variant).upper(),
            "mapping_type": "cst_index",
            "suggested_gene_symbols": list(indexed.get("suggested_gene_symbols") or []),
            "suggested_uniprot_ids": list(indexed.get("suggested_uniprot_ids") or []),
            "notes": indexed.get("notes") or "Mapped from CST pathway protein-module index.",
            "confidence": "high",
        }

    mapper = _get_cst_label_mapper()
    mapping_rank = {
        "cst_index": 5,
        "exact_gene": 4,
        "alias": 3,
        "family": 2,
        "ambiguous": 1,
        "unresolved": 0,
    }
    best_result: Optional[Dict[str, Any]] = None
    best_variant = ""
    best_rank = -1
    for variant in _cst_label_variants(str(label or "")):
        result = mapper.map_pathway_label(variant)
        rank = mapping_rank.get(str(result.get("mapping_type") or "").strip().lower(), 0)
        if rank > best_rank:
            best_result = result
            best_variant = variant
            best_rank = rank
        if rank >= mapping_rank["exact_gene"]:
            break
    fallback = dict(best_result or mapper.map_pathway_label(str(label or "")))
    if best_variant and best_variant != str(label or "") and str(fallback.get("mapping_type") or "").strip().lower() != "unresolved":
        fallback["original_label"] = str(label or "")
        fallback["normalized_label"] = _clean_text_label(str(best_variant)).upper()
        notes = str(fallback.get("notes") or "").strip()
        fallback["notes"] = f"{notes} Recovered from CST text fragment '{best_variant}'.".strip()
    if not list(fallback.get("suggested_uniprot_ids") or []):
        forced_index = load_forced_review_mapping_index()
        pathway_key = _normalize_forced_pathway_key(file_path)
        forced_match = None
        forced_variant = ""
        for variant in _cst_label_variants(str(label or "")):
            forced_key = _canonicalize_label(_clean_text_label(variant).upper())
            candidate = forced_index.get((pathway_key, forced_key))
            if candidate:
                forced_match = dict(candidate)
                forced_variant = variant
                break
        if forced_match:
            return {
                "original_label": str(label or ""),
                "normalized_label": _clean_text_label(str(forced_variant or label or "")).upper(),
                "mapping_type": str(forced_match.get("mapping_type") or "forced_review"),
                "suggested_gene_symbols": list(forced_match.get("suggested_gene_symbols") or []),
                "suggested_uniprot_ids": list(forced_match.get("suggested_uniprot_ids") or []),
                "notes": str(forced_match.get("notes") or "Mapped from curated forced-review pathway label table."),
                "confidence": str(forced_match.get("confidence") or "medium-high"),
            }
        complex_index = load_complex_review_mapping_index()
        complex_match = None
        complex_variant = ""
        for variant in _cst_label_variants(str(label or "")):
            complex_key = _canonicalize_label(_clean_text_label(variant).upper())
            candidate = complex_index.get((pathway_key, complex_key))
            if candidate:
                complex_match = dict(candidate)
                complex_variant = variant
                break
        if complex_match:
            return {
                "original_label": str(label or ""),
                "normalized_label": _clean_text_label(str(complex_variant or label or "")).upper(),
                "mapping_type": str(complex_match.get("mapping_type") or "complex_forced"),
                "suggested_gene_symbols": list(complex_match.get("suggested_gene_symbols") or []),
                "suggested_uniprot_ids": list(complex_match.get("suggested_uniprot_ids") or []),
                "complex_components": list(complex_match.get("complex_components") or []),
                "notes": str(complex_match.get("notes") or "Mapped from curated complex/process component table."),
                "confidence": str(complex_match.get("confidence") or "medium-high"),
            }
    return fallback


def _resolve_cst_label_mapping(label: str, file_path: Path) -> Dict[str, Any]:
    return dict(_resolve_cst_label_mapping_cached(str(label or ""), str(file_path)))


@lru_cache(maxsize=16)
def _map_cst_text_nodes(file_path_str: str) -> List[Dict[str, Any]]:
    nodes = _extract_cst_text_nodes(file_path_str)
    file_path = Path(file_path_str)
    output: List[Dict[str, Any]] = []
    for node in nodes:
        payload = dict(node)
        payload["mapping"] = _resolve_cst_label_mapping(str(node.get("label") or ""), file_path)
        output.append(payload)
    return output


def _build_dataset_index(dataset: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not dataset:
        return {
            "fc_headers": [],
            "outline_headers": [],
            "rows_by_uniprot": {},
            "rows_by_gene": {},
            "tooltip_columns": [],
            "headers": [],
            "display_uniprot_key": "",
            "display_gene_key": "",
            "match_uniprot_key": "",
            "match_gene_key": "",
        }
    headers = list(dataset.get("headers") or [])
    rows = list(dataset.get("rows") or [])
    if not headers:
        return {
            "fc_headers": [],
            "outline_headers": [],
            "rows_by_uniprot": {},
            "rows_by_gene": {},
            "tooltip_columns": [],
            "headers": [],
            "display_uniprot_key": "",
            "display_gene_key": "",
            "match_uniprot_key": "",
            "match_gene_key": "",
        }
    fc_headers = _normalize_fc_headers(dataset.get("main_columns") or [h for h in headers if str(h).startswith("C:")])
    outline_headers = list(dataset.get("outline_columns") or _resolve_outline_columns(fc_headers, headers))
    tooltip_columns = list(dataset.get("tooltip_columns") or [])
    if not tooltip_columns:
        if len(headers) > 1:
            tooltip_columns = [headers[1], headers[0]]
        elif headers:
            tooltip_columns = [headers[0]]
        tooltip_columns.extend([header for header in headers if str(header).startswith("T:") and header not in tooltip_columns])
    rows_by_uniprot: Dict[str, Dict[str, Any]] = {}
    rows_by_gene: Dict[str, List[Dict[str, Any]]] = {}
    dataset_keys = _resolve_cst_dataset_keys(headers)
    uniprot_key = str(dataset_keys.get("match_uniprot_key") or "").strip()
    gene_key = str(dataset_keys.get("match_gene_key") or "").strip()
    display_uniprot_key = str(dataset_keys.get("display_uniprot_key") or "").strip()
    display_gene_key = str(dataset_keys.get("display_gene_key") or "").strip()

    def _normalize_uniprot_token(value: Any) -> str:
        token = str(value or "").strip().upper()
        if not token:
            return ""
        primary = re.split(r"[|,;/]+", token, maxsplit=1)[0].strip()
        return primary.split("-", 1)[0].strip() if primary else ""

    def _collect_gene_tokens(*values: Any) -> List[str]:
        out: List[str] = []
        seen: set[str] = set()
        for value in values:
            raw = str(value or "").strip().upper()
            if not raw:
                continue
            for token in re.split(r"[|,;/]+", raw):
                gene = str(token or "").strip().upper()
                if not gene or gene in seen:
                    continue
                seen.add(gene)
                out.append(gene)
        return out

    for row in rows:
        values = list(row)
        if len(values) < len(headers):
            values.extend([""] * (len(headers) - len(values)))
        row_map = {header: values[idx] if idx < len(values) else "" for idx, header in enumerate(headers)}
        normalized_uniprot = _normalize_uniprot_token(row_map.get(uniprot_key, ""))
        if not normalized_uniprot:
            normalized_uniprot = _normalize_uniprot_token(row_map.get(display_uniprot_key, ""))
        if normalized_uniprot and normalized_uniprot not in rows_by_uniprot:
            rows_by_uniprot[normalized_uniprot] = row_map
        for raw_gene in _collect_gene_tokens(
            row_map.get(gene_key, "") if gene_key else "",
            row_map.get(display_gene_key, "") if display_gene_key else "",
        ):
            rows_by_gene.setdefault(raw_gene, []).append(row_map)
    return {
        "fc_headers": fc_headers,
        "outline_headers": outline_headers,
        "rows_by_uniprot": rows_by_uniprot,
        "rows_by_gene": rows_by_gene,
        "tooltip_columns": tooltip_columns,
        "headers": headers,
        "display_uniprot_key": str(dataset_keys.get("display_uniprot_key") or ""),
        "display_gene_key": str(dataset_keys.get("display_gene_key") or ""),
        "match_uniprot_key": str(dataset_keys.get("match_uniprot_key") or ""),
        "match_gene_key": str(dataset_keys.get("match_gene_key") or ""),
    }


def _build_cst_search_catalog(
    protein_dataset: Optional[Dict[str, Any]],
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
    use_black_protein_outlines: bool = False,
) -> List[Dict[str, Any]]:
    dataset_index = _build_dataset_index(protein_dataset)
    headers = list(dataset_index.get("headers") or [])
    if not headers:
        return []
    uniprot_key = str(dataset_index.get("match_uniprot_key") or (headers[0] if headers else "")).strip()
    display_uniprot_key = str(dataset_index.get("display_uniprot_key") or (headers[0] if headers else "")).strip()
    gene_key = str(dataset_index.get("display_gene_key") or (headers[1] if len(headers) > 1 else "")).strip()
    match_gene_key = str(dataset_index.get("match_gene_key") or "").strip()
    fc_headers = list(dataset_index.get("fc_headers") or [])
    outline_headers = list(dataset_index.get("outline_headers") or [])
    tooltip_columns = list(dataset_index.get("tooltip_columns") or [])
    catalog: List[Dict[str, Any]] = []

    def _collect_search_aliases(row_map: Dict[str, Any], gene_symbol: str, normalized_uniprot: str) -> List[str]:
        aliases: List[str] = []
        seen: set[str] = set()

        def _add(value: Any) -> None:
            raw = str(value or "").strip()
            if not raw:
                return
            for token in re.split(r"[|,;/]+", raw):
                item = str(token or "").strip()
                if not item:
                    continue
                for candidate in (item, item.split("-", 1)[0]):
                    normalized = str(candidate or "").strip()
                    if not normalized:
                        continue
                    key = normalized.upper()
                    if key in seen:
                        continue
                    seen.add(key)
                    aliases.append(normalized)

        _add(normalized_uniprot)
        _add(row_map.get(uniprot_key, ""))
        _add(row_map.get(display_uniprot_key, ""))
        _add(gene_symbol)
        if match_gene_key:
            _add(row_map.get(match_gene_key, ""))
        return aliases

    raw_rows = list((protein_dataset or {}).get("rows") or [])
    row_maps: List[Dict[str, Any]] = []
    for row in raw_rows:
        if isinstance(row, dict):
            row_map = {str(k): ("" if pd.isna(v) else v) for k, v in row.items()}
        else:
            values = list(row or [])
            if len(values) < len(headers):
                values.extend([""] * (len(headers) - len(values)))
            row_map = {header: (values[idx] if idx < len(values) else "") for idx, header in enumerate(headers)}
        row_maps.append(row_map)
    if not row_maps:
        row_maps = [dict(row or {}) for row in list((dataset_index.get("rows_by_uniprot") or {}).values())]

    def _normalize_uniprot(value: Any) -> str:
        text = str(value or "").strip().upper()
        if not text:
            return ""
        token = re.split(r"[|,;/]+", text, maxsplit=1)[0].strip()
        token = token.split("-", 1)[0].strip()
        return token

    seen_entry_keys: set[str] = set()
    for row_map in row_maps:
        uniprot = _normalize_uniprot(row_map.get(uniprot_key, ""))
        if not uniprot:
            uniprot = _normalize_uniprot(row_map.get(display_uniprot_key, ""))
        gene_symbol = str(row_map.get(gene_key, "") or "").strip() if gene_key else ""
        if not gene_symbol and match_gene_key:
            gene_symbol = str(row_map.get(match_gene_key, "") or "").strip()
        dedupe_key = uniprot or f"GENE::{gene_symbol.upper()}"
        if not dedupe_key or dedupe_key in seen_entry_keys:
            continue
        seen_entry_keys.add(dedupe_key)
        entry: Dict[str, Any] = {
            "uniprot": str(uniprot or "").strip(),
            "geneSymbol": gene_symbol,
            "searchAliases": _collect_search_aliases(row_map, gene_symbol, str(uniprot or "").strip()),
        }
        for idx, header in enumerate(fc_headers, start=1):
            fold_value = row_map.get(header)
            entry[f"fold_change_{idx}"] = fold_value
            entry[f"fc_color_{idx}"] = _gradient_color_from_fold(
                fold_value,
                negative_color,
                positive_color,
                max_negative,
                max_positive,
            )
            outline_header = outline_headers[idx - 1] if idx - 1 < len(outline_headers) else None
            outline_value = _coerce_float_or_none(row_map.get(outline_header, "")) if outline_header else None
            entry[f"outline_fold_change_{idx}"] = outline_value
            entry[f"outline_color_{idx}"] = (
                [0, 0, 0]
                if use_black_protein_outlines
                else (
                    _gradient_color_from_fold(
                        outline_value,
                        negative_color,
                        positive_color,
                        max_negative,
                        max_positive,
                    )
                    if outline_value is not None
                    else [0, 0, 0]
                )
            )
        entry.update(_build_dataset_tooltip(row_map, tooltip_columns))
        catalog.append(entry)
    catalog.sort(key=lambda item: (str(item.get("geneSymbol") or "").upper(), str(item.get("uniprot") or "").upper()))
    return catalog


def _build_cst_metabolite_search_catalog(
    metabolite_dataset: Optional[Dict[str, Any]],
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
) -> List[Dict[str, Any]]:
    if not metabolite_dataset:
        return []
    headers = list(metabolite_dataset.get("headers") or [])
    rows = list(metabolite_dataset.get("rows") or [])
    hmdb_column = str(metabolite_dataset.get("hmdb_column") or (headers[0] if headers else "")).strip()
    if not headers or not rows or not hmdb_column:
        return []
    main_headers = list(metabolite_dataset.get("main_columns") or [header for header in headers if str(header).startswith("C:")])
    tooltip_columns = list(metabolite_dataset.get("tooltip_columns") or [header for header in headers if str(header).startswith("T:")])
    _, hmdb_reference_map = _load_hmdb_metabolite_reference()
    catalog: List[Dict[str, Any]] = []
    seen_ids: set[str] = set()
    for row in rows:
        if isinstance(row, dict):
            row_map = {str(k): ("" if pd.isna(v) else v) for k, v in row.items()}
        else:
            row_map = {
                str(headers[idx]): ("" if pd.isna(value) else value)
                for idx, value in enumerate(list(row or []))
                if idx < len(headers)
            }
        hmdb_id = str(row_map.get(hmdb_column, "") or "").strip()
        if not hmdb_id:
            continue
        hmdb_key = _normalize_metabolite_token(hmdb_id)
        if hmdb_key in seen_ids:
            continue
        seen_ids.add(hmdb_key)
        ref = hmdb_reference_map.get(hmdb_key)
        wiki_id = _display_metabolite_name((ref or {}).get("wikipedia_id") or "")
        kegg_id = str((ref or {}).get("kegg_id") or "").strip()
        display_label = _display_metabolite_name(wiki_id or (ref or {}).get("name") or hmdb_id)
        entry: Dict[str, Any] = {
            "hmdbId": hmdb_id,
            "wikipediaId": wiki_id,
            "keggId": kegg_id,
            "displayLabel": display_label,
            "name": str((ref or {}).get("name") or "").strip(),
            "synonyms": str((ref or {}).get("synonyms") or "").strip(),
            "tooltip": "",
            "tooltip_html": "",
        }
        tooltip_plain: List[str] = []
        tooltip_html: List[str] = []
        if display_label:
            tooltip_plain.append(f"Metabolite: {display_label}")
            tooltip_html.append(f"<strong>Metabolite:</strong> {html.escape(display_label)}")
        tooltip_plain.append(f"HMDB ID: {hmdb_id}")
        tooltip_html.append(f"<strong>HMDB ID:</strong> {html.escape(hmdb_id)}")
        if wiki_id:
            tooltip_plain.append(f"Wiki ID: {wiki_id}")
            tooltip_html.append(f"<strong>Wiki ID:</strong> {html.escape(wiki_id)}")
        for col in tooltip_columns:
            value = str(row_map.get(col, "") or "").strip()
            if not value:
                continue
            tooltip_plain.append(f"{col}: {value}")
            tooltip_html.append(f"<strong>{html.escape(str(col))}:</strong> {html.escape(value)}")
        entry["tooltip"] = "\n".join(tooltip_plain)
        entry["tooltip_html"] = "<br/>".join(tooltip_html)
        for idx, header in enumerate(main_headers, start=1):
            fold_value = _coerce_float_or_none(row_map.get(header))
            entry[f"fold_change_{idx}"] = fold_value
            entry[f"fc_color_{idx}"] = _gradient_color_from_fold(
                fold_value,
                negative_color,
                positive_color,
                max_negative,
                max_positive,
            )
        catalog.append(entry)
    catalog.sort(key=lambda item: (str(item.get("displayLabel") or "").upper(), str(item.get("hmdbId") or "").upper()))
    return catalog


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


def _rgb_text(rgb: Sequence[float]) -> str:
    values = [int(max(0, min(255, round(float(item))))) for item in list(rgb or [166, 166, 166])[:3]]
    while len(values) < 3:
        values.append(166)
    return f"rgb({values[0]}, {values[1]}, {values[2]})"


def _coerce_float_or_none(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _build_dataset_tooltip(row_map: Optional[Dict[str, Any]], tooltip_columns: Sequence[Any]) -> Dict[str, str]:
    if not isinstance(row_map, dict):
        return {"tooltip": "", "tooltip_html": ""}
    plain_lines: List[str] = []
    html_lines: List[str] = []
    for col in list(tooltip_columns or []):
        label = str(col or "").strip()
        if not label:
            continue
        value = str(row_map.get(label, "") or "").strip()
        if not value:
            continue
        plain_lines.append(f"{label}: {value}")
        html_lines.append(f"<strong>{html.escape(label)}:</strong> {html.escape(value)}")
    return {
        "tooltip": "\n".join(plain_lines),
        "tooltip_html": "<br/>".join(html_lines),
    }


def _select_dataset_tooltip_row(
    matches: Sequence[Dict[str, Any]],
    dataset_index: Dict[str, Any],
    matched_uniprot: str = "",
    matched_gene: str = "",
) -> Optional[Dict[str, Any]]:
    rows = [row for row in list(matches or []) if isinstance(row, dict)]
    if not rows:
        return None
    headers = list(dataset_index.get("headers") or [])
    uniprot_key = str(dataset_index.get("match_uniprot_key") or (headers[0] if headers else "")).strip()
    gene_key = str(dataset_index.get("display_gene_key") or (headers[1] if len(headers) > 1 else "")).strip()
    matched_uniprot_norm = str(matched_uniprot or "").strip().upper().split("-", 1)[0]
    matched_gene_norm = str(matched_gene or "").strip().upper()
    if matched_uniprot_norm and uniprot_key:
        for row in rows:
            row_uniprot = str(row.get(uniprot_key, "") or "").strip().upper().split("-", 1)[0]
            if row_uniprot == matched_uniprot_norm:
                return row
    if matched_gene_norm and gene_key:
        for row in rows:
            row_gene = str(row.get(gene_key, "") or "").strip().upper()
            if row_gene == matched_gene_norm:
                return row
    return rows[0]


def _build_cst_protein_options(
    mapping: Dict[str, Any],
    dataset_index: Dict[str, Any],
    tooltip_columns: Sequence[Any],
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
    use_black_protein_outlines: bool = False,
) -> List[Dict[str, Any]]:
    mapping_type = str(mapping.get("mapping_type") or "").strip().lower()
    if mapping_type in {"", "unresolved", "complex_forced"}:
        return []
    matches = _match_dataset_rows(mapping, dataset_index)
    headers = list(dataset_index.get("headers") or [])
    fc_headers = list(dataset_index.get("fc_headers") or [])
    outline_headers = list(dataset_index.get("outline_headers") or [])
    if not matches or not headers:
        return []
    match_uniprot_key = str(dataset_index.get("match_uniprot_key") or (headers[0] if headers else "")).strip()
    display_gene_key = str(dataset_index.get("display_gene_key") or (headers[1] if len(headers) > 1 else "")).strip()
    match_gene_key = str(dataset_index.get("match_gene_key") or "").strip()
    options: List[Dict[str, Any]] = []
    seen: set[str] = set()
    for row in matches:
        uniprot = str(row.get(match_uniprot_key, "") or "").strip()
        if not uniprot or uniprot in seen:
            continue
        seen.add(uniprot)
        gene_symbol = str(row.get(display_gene_key, "") or "").strip() if display_gene_key else ""
        match_gene_symbol = str(row.get(match_gene_key, "") or "").strip() if match_gene_key else gene_symbol
        option: Dict[str, Any] = {
            "uniprot": uniprot,
            "gene_symbol": gene_symbol,
            "match_gene_symbol": match_gene_symbol,
        }
        option.update(_build_dataset_tooltip(row, tooltip_columns))
        for idx, header in enumerate(fc_headers, 1):
            raw_value = row.get(header, "")
            fold_value = _coerce_float_or_none(raw_value)
            option[f"fold_change_{idx}"] = fold_value
            option[f"fc_color_{idx}"] = _gradient_color_from_fold(
                fold_value,
                negative_color,
                positive_color,
                max_negative,
                max_positive,
            )
            outline_header = outline_headers[idx - 1] if idx - 1 < len(outline_headers) else None
            outline_value = _coerce_float_or_none(row.get(outline_header, "")) if outline_header else None
            option[f"outline_fold_change_{idx}"] = outline_value
            option[f"outline_color_{idx}"] = (
                [0, 0, 0]
                if use_black_protein_outlines
                else (
                    _gradient_color_from_fold(
                        outline_value,
                        negative_color,
                        positive_color,
                        max_negative,
                        max_positive,
                    )
                    if outline_value is not None
                    else [0, 0, 0]
                )
            )
        options.append(option)
    return options


def _infer_ptm_tooltip_columns(dataset: Optional[Dict[str, Any]]) -> List[str]:
    if not dataset:
        return []
    explicit = [str(col).strip() for col in list(dataset.get("tooltip_columns") or []) if str(col).strip()]
    if explicit:
        return explicit
    headers = [str(col).strip() for col in list(dataset.get("headers") or []) if str(col).strip()]
    columns: List[str] = []
    for header in headers:
        if header.startswith("T:") and header not in columns:
            columns.append(header)
    for header in headers:
        if header.startswith("PSP:") and header != "PSP: regulatory_site" and header not in columns:
            columns.append(header)
    return columns


def _choose_complex_display_members(
    mapping: Dict[str, Any],
    dataset_index: Dict[str, Any],
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
) -> List[Dict[str, Any]]:
    headers = list(dataset_index.get("headers") or [])
    fc_headers = list(dataset_index.get("fc_headers") or [])
    rows_by_uniprot = dataset_index.get("rows_by_uniprot", {})
    rows_by_gene = dataset_index.get("rows_by_gene", {})
    match_uniprot_key = str(dataset_index.get("match_uniprot_key") or (headers[0] if headers else "")).strip()
    display_gene_key = str(dataset_index.get("display_gene_key") or (headers[1] if len(headers) > 1 else "")).strip()
    output: List[Dict[str, Any]] = []
    for component in list(mapping.get("complex_components") or []):
        component_label = str(component.get("component_label") or "").strip()
        for group in list(component.get("member_groups") or []):
            candidate_genes = [str(item or "").strip().upper() for item in list(group.get("candidate_gene_symbols") or []) if str(item or "").strip()]
            candidate_ids = [str(item or "").strip().upper().split("-", 1)[0] for item in list(group.get("candidate_uniprot_ids") or []) if str(item or "").strip()]
            chosen_gene = component_label or (candidate_genes[0] if candidate_genes else "")
            chosen_uniprot = candidate_ids[0] if candidate_ids else ""
            chosen_row: Optional[Dict[str, Any]] = None
            chosen_fc: Optional[float] = None

            for gene_symbol in candidate_genes:
                matched_rows = list(rows_by_gene.get(gene_symbol, []))
                if not matched_rows:
                    continue
                chosen_gene = gene_symbol
                chosen_row = matched_rows[0]
                break
            if chosen_row is None:
                for uniprot_id in candidate_ids:
                    matched_row = rows_by_uniprot.get(uniprot_id)
                    if matched_row is None:
                        continue
                    chosen_uniprot = uniprot_id
                    chosen_row = matched_row
                    if display_gene_key:
                        chosen_gene = str(matched_row.get(display_gene_key, "") or chosen_gene)
                    break

            fill_rgb = [229, 231, 235]
            matched_uniprot = chosen_uniprot
            matched_gene = chosen_gene
            if chosen_row is not None and fc_headers:
                best_value: Optional[float] = None
                for header in fc_headers:
                    raw_value = chosen_row.get(header, "")
                    try:
                        value = float(raw_value)
                    except (TypeError, ValueError):
                        continue
                    if best_value is None or abs(value) > abs(best_value):
                        best_value = value
                if best_value is not None:
                    chosen_fc = best_value
                    fill_rgb = _gradient_color_from_fold(
                        best_value,
                        negative_color,
                        positive_color,
                        max_negative,
                        max_positive,
                    )
                if headers:
                    matched_uniprot = str(chosen_row.get(match_uniprot_key, "") or matched_uniprot)
                    if display_gene_key:
                        matched_gene = str(chosen_row.get(display_gene_key, "") or matched_gene)

            output.append({
                "component_label": component_label,
                "label": matched_gene or chosen_gene or component_label or "Component",
                "matched_gene_symbol": matched_gene or chosen_gene,
                "matched_uniprot": matched_uniprot,
                "has_dataset_match": chosen_row is not None,
                "fold_change": chosen_fc,
                "fill_color": list(fill_rgb),
                "fill_color_css": _rgb_text(fill_rgb),
            })
    return output


def _build_cst_overlay_nodes(
    file_path: Path,
    protein_dataset: Optional[Dict[str, Any]] = None,
    ptm_index: Optional[Dict[str, Any]] = None,
    negative_color: Sequence[float] = _DEFAULT_NEGATIVE_COLOR,
    positive_color: Sequence[float] = _DEFAULT_POSITIVE_COLOR,
    max_negative: float = _DEFAULT_MAX_NEGATIVE,
    max_positive: float = _DEFAULT_MAX_POSITIVE,
    prot_outline_width: float = 1.0,
    use_black_protein_outlines: bool = False,
) -> List[Dict[str, Any]]:
    mapped_nodes = _map_cst_text_nodes(str(file_path))
    if not mapped_nodes:
        return []
    dataset_index = _build_dataset_index(protein_dataset)
    fc_headers = list(dataset_index.get("fc_headers") or [])
    outline_headers = list(dataset_index.get("outline_headers") or [])
    tooltip_columns = list(dataset_index.get("tooltip_columns") or [])
    dataset_headers = list(dataset_index.get("headers") or [])
    match_uniprot_key = str(dataset_index.get("match_uniprot_key") or (dataset_headers[0] if dataset_headers else "")).strip()
    display_uniprot_key = str(dataset_index.get("display_uniprot_key") or (dataset_headers[0] if dataset_headers else "")).strip()
    display_gene_key = str(dataset_index.get("display_gene_key") or (dataset_headers[1] if len(dataset_headers) > 1 else "")).strip()
    overlay_nodes: List[Dict[str, Any]] = []

    for node in mapped_nodes:
        mapping = dict(node.get("mapping") or {})
        mapping_type = str(mapping.get("mapping_type") or "").strip().lower()
        if mapping_type == "unresolved":
            continue
        is_complex = mapping_type == "complex_forced"
        matches = [] if is_complex else _match_dataset_rows(mapping, dataset_index)
        candidate_uniprot_ids: List[str] = []
        candidate_seen: set[str] = set()

        def _add_candidate_uniprot(value: Any) -> None:
            raw = str(value or "").strip()
            if not raw:
                return
            normalized = re.split(r"[|,;/]+", raw, maxsplit=1)[0].strip()
            if not normalized:
                return
            upper = normalized.upper()
            if upper in candidate_seen:
                return
            candidate_seen.add(upper)
            candidate_uniprot_ids.append(normalized)

        for item in list(mapping.get("suggested_uniprot_ids") or []):
            _add_candidate_uniprot(item)
        for row in matches:
            _add_candidate_uniprot(row.get(match_uniprot_key, "") if match_uniprot_key else "")
            _add_candidate_uniprot(row.get(display_uniprot_key, "") if display_uniprot_key else "")
        protein_options = _build_cst_protein_options(
            mapping,
            dataset_index,
            tooltip_columns,
            negative_color,
            positive_color,
            max_negative,
            max_positive,
            use_black_protein_outlines=use_black_protein_outlines,
        ) if not is_complex else []
        complex_members = _choose_complex_display_members(
            mapping,
            dataset_index,
            negative_color,
            positive_color,
            max_negative,
            max_positive,
        ) if is_complex else []

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
            "candidate_uniprot_ids": candidate_uniprot_ids,
            "protein_options": protein_options,
            "has_dataset_match": bool(matches),
            "default_color": [229, 231, 235] if is_complex else [166, 166, 166],
            "is_complex": is_complex,
            "complex_components": list(mapping.get("complex_components") or []),
            "complex_display_members": complex_members,
            "tooltip": "",
            "tooltip_html": "",
        }
        if matches:
            tooltip_row = _select_dataset_tooltip_row(matches, dataset_index)
            overlay_node.update(_build_dataset_tooltip(tooltip_row, tooltip_columns))

        for idx, header in enumerate(fc_headers, 1):
            chosen_row, chosen_value = _choose_cst_display_row(
                matches,
                dataset_index.get("headers") or [],
                header,
                ptm_index,
                idx,
            )
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
            outline_header = outline_headers[idx - 1] if idx - 1 < len(outline_headers) else None
            outline_value = _coerce_float_or_none(chosen_row.get(outline_header, "")) if outline_header else None
            overlay_node[f"outline_fold_change_{idx}"] = outline_value
            overlay_node[f"outline_color_{idx}"] = (
                [0, 0, 0]
                if use_black_protein_outlines
                else (
                    _gradient_color_from_fold(
                        outline_value,
                        negative_color,
                        positive_color,
                        max_negative,
                        max_positive,
                    )
                    if outline_value is not None
                    else [0, 0, 0]
                )
            )
            if dataset_headers:
                matched_uni = str(chosen_row.get(match_uniprot_key, "") or "").strip() if match_uniprot_key else ""
                if not matched_uni and display_uniprot_key:
                    matched_uni = str(chosen_row.get(display_uniprot_key, "") or "").strip()
                overlay_node[f"matched_uniprot_{idx}"] = matched_uni
                overlay_node[f"matched_gene_symbol_{idx}"] = str(chosen_row.get(display_gene_key, "") or "") if display_gene_key else ""
        overlay_node["stroke_width"] = max(0.6, float(prot_outline_width or 1.0))

        overlay_nodes.append(overlay_node)
    return overlay_nodes


def _build_editable_node_from_saved_entry(
    file_path: Path,
    entry: Dict[str, Any],
    dataset_index: Dict[str, Any],
    ptm_index: Optional[Dict[str, Any]],
    negative_color: Sequence[float],
    positive_color: Sequence[float],
    max_negative: float,
    max_positive: float,
    prot_outline_width: float = 1.0,
    use_black_protein_outlines: bool = False,
) -> Optional[Dict[str, Any]]:
    shape_type = str(entry.get("shapeType") or entry.get("shape_type") or "ellipse").strip().lower() or "ellipse"
    if bool(entry.get("isDrawingShape") or entry.get("is_drawing_shape")):
        return {
            "id": str(entry.get("id") or ""),
            "originalLabel": "",
            "displayLabel": "",
            "label": "",
            "matchedGene": "",
            "matchedUniprot": "",
            "foldText": "",
            "hasDatasetMatch": False,
            "cx": float(entry.get("cx") or 0.0),
            "cy": float(entry.get("cy") or 0.0),
            "rx": max(10.0, float(entry.get("rx") or 60.0)),
            "ry": max(10.0, float(entry.get("ry") or 40.0)),
            "shapeType": shape_type,
            "legendOrientation": str(entry.get("legendOrientation") or entry.get("legend_orientation") or "vertical"),
            "angle": float(entry.get("angle") or 0.0),
            "strokeWidth": max(1.0, float(entry.get("strokeWidth") or 1.0)),
            "stroke": str(entry.get("stroke") or "#000000"),
            "fillColor": str(entry.get("fillColor") or ("transparent" if shape_type == "bracket" else "#f5f5f5")),
            "opacity": max(0.1, min(1.0, float(entry.get("opacity") or 1.0))),
            "annotation": "",
            "annotationCommitted": True,
            "pendingAnnotation": "",
            "isCustom": True,
            "userCreated": True,
            "isDrawingShape": True,
            "className": "cst-overlay-shape",
            "title": str(entry.get("title") or "Shape"),
            "mappingType": "shape",
            "isComplex": False,
            "complexDisplayMembers": [],
        }
    if shape_type == "text":
        text_label = _normalize_pdf_text_encoding(
            entry.get("displayLabel") or entry.get("display_label") or entry.get("label") or "Text"
        )
        text_label = text_label.replace("\x1d", "-").replace("\x1e", "-").replace("\x1f", "-").strip() or "Text"
        return {
            "id": str(entry.get("id") or ""),
            "originalLabel": text_label,
            "displayLabel": text_label,
            "label": text_label,
            "matchedGene": "",
            "matchedUniprot": "",
            "foldText": "",
            "hasDatasetMatch": False,
            "cx": float(entry.get("cx") or 0.0),
            "cy": float(entry.get("cy") or 0.0),
            "rx": max(10.0, float(entry.get("rx") or 44.0)),
            "ry": max(8.0, float(entry.get("ry") or 13.0)),
            "shapeType": "text",
            "angle": float(entry.get("angle") or 0.0),
            "strokeWidth": max(1.0, float(entry.get("strokeWidth") or 1.5)),
            "stroke": str(entry.get("stroke") or "rgb(71, 85, 105)"),
            "fillColor": str(entry.get("fillColor") or "transparent"),
            "textColor": str(entry.get("textColor") or "#0f172a"),
            "opacity": max(0.1, min(1.0, float(entry.get("opacity") or 1.0))),
            "annotation": "",
            "annotationCommitted": True,
            "pendingAnnotation": "",
            "isCustom": bool(entry.get("isCustom")),
            "userCreated": bool(entry.get("userCreated")),
            "isDrawingShape": False,
            "className": "cst-overlay-text",
            "title": str(entry.get("title") or "Text box"),
            "mappingType": "text",
            "isComplex": False,
            "complexDisplayMembers": [],
            "fontSize": max(9.0, float(entry.get("fontSize") or 14.0)),
            "fontWeight": str(entry.get("fontWeight") or "600"),
            "fontFamily": str(entry.get("fontFamily") or '"Segoe UI", Arial, sans-serif'),
            "textAlign": str(entry.get("textAlign") or "center"),
        }
    original_label = _clean_text_label(str(entry.get("originalLabel") or entry.get("original_label") or entry.get("label") or ""))
    display_label = _clean_text_label(
        str(entry.get("displayLabel") or entry.get("display_label") or entry.get("annotation") or original_label)
    )
    mapping_label = display_label or original_label
    if not mapping_label and not bool(entry.get("isCustom") or entry.get("is_custom")):
        return None

    saved_gene_symbols = [str(x).strip() for x in list(entry.get("suggested_gene_symbols") or []) if str(x).strip()]
    saved_uniprot_ids = [str(x).strip() for x in list(entry.get("suggested_uniprot_ids") or []) if str(x).strip()]
    saved_complex_components = [
        dict(x)
        for x in list(entry.get("complex_components") or [])
        if isinstance(x, dict)
    ]
    saved_mapping_type = str(entry.get("mappingType") or entry.get("mapping_type") or "").strip().lower()
    if saved_mapping_type == "complex_forced" and not saved_complex_components:
        mapping = _resolve_cst_label_mapping(mapping_label, file_path) if mapping_label else {
            "mapping_type": "unresolved",
            "suggested_gene_symbols": [],
            "suggested_uniprot_ids": [],
            "complex_components": [],
        }
    elif saved_mapping_type and (saved_gene_symbols or saved_uniprot_ids or saved_complex_components):
        mapping = {
            "mapping_type": saved_mapping_type,
            "suggested_gene_symbols": saved_gene_symbols,
            "suggested_uniprot_ids": saved_uniprot_ids,
            "complex_components": saved_complex_components,
        }
    else:
        mapping = _resolve_cst_label_mapping(mapping_label, file_path) if mapping_label else {
            "mapping_type": "unresolved",
            "suggested_gene_symbols": [],
            "suggested_uniprot_ids": [],
            "complex_components": [],
        }
    mapping_type = str(mapping.get("mapping_type") or "").strip().lower() or "unresolved"
    is_complex = mapping_type == "complex_forced"
    matches = _match_dataset_rows(mapping, dataset_index) if mapping_type not in {"unresolved", "complex_forced"} else []
    fc_headers = list(dataset_index.get("fc_headers") or [])
    outline_headers = list(dataset_index.get("outline_headers") or [])
    tooltip_columns = list(dataset_index.get("tooltip_columns") or [])
    headers = list(dataset_index.get("headers") or [])
    match_uniprot_key = str(dataset_index.get("match_uniprot_key") or (headers[0] if headers else "")).strip()
    display_uniprot_key = str(dataset_index.get("display_uniprot_key") or (headers[0] if headers else "")).strip()
    display_gene_key = str(dataset_index.get("display_gene_key") or (headers[1] if len(headers) > 1 else "")).strip()
    complex_members = _choose_complex_display_members(
        mapping,
        dataset_index,
        negative_color,
        positive_color,
        max_negative,
        max_positive,
    ) if is_complex else []
    protein_options = _build_cst_protein_options(
        mapping,
        dataset_index,
        tooltip_columns,
        negative_color,
        positive_color,
        max_negative,
        max_positive,
        use_black_protein_outlines=use_black_protein_outlines,
    ) if not is_complex else []
    fill_rgb = [229, 231, 235] if is_complex else [166, 166, 166]
    stroke_rgb = [37, 99, 235] if is_complex else [0, 0, 0]
    matched_gene = ""
    matched_uniprot = ""
    fold_text = ""
    candidate_uniprot_ids: List[str] = []
    candidate_seen: set[str] = set()

    def _add_candidate_uniprot(value: Any) -> None:
        raw = str(value or "").strip()
        if not raw:
            return
        normalized = re.split(r"[|,;/]+", raw, maxsplit=1)[0].strip()
        if not normalized:
            return
        upper = normalized.upper()
        if upper in candidate_seen:
            return
        candidate_seen.add(upper)
        candidate_uniprot_ids.append(normalized)

    for item in list(mapping.get("suggested_uniprot_ids") or []):
        _add_candidate_uniprot(item)
    for row in matches:
        _add_candidate_uniprot(row.get(match_uniprot_key, "") if match_uniprot_key else "")
        _add_candidate_uniprot(row.get(display_uniprot_key, "") if display_uniprot_key else "")

    temporal_payload: Dict[str, Any] = {}
    for idx, header in enumerate(fc_headers, 1):
        chosen_row, chosen_value = _choose_cst_display_row(
            matches,
            headers,
            header,
            ptm_index,
            idx,
        )
        if chosen_value is None or chosen_row is None:
            continue
        color = _gradient_color_from_fold(
            chosen_value,
            negative_color,
            positive_color,
            max_negative,
            max_positive,
        )
        temporal_payload[f"fold_change_{idx}"] = chosen_value
        temporal_payload[f"fc_color_{idx}"] = list(color)
        outline_header = outline_headers[idx - 1] if idx - 1 < len(outline_headers) else None
        outline_value = _coerce_float_or_none(chosen_row.get(outline_header, "")) if outline_header else None
        outline_color = (
            [0, 0, 0]
            if use_black_protein_outlines
            else (
                _gradient_color_from_fold(
                    outline_value,
                    negative_color,
                    positive_color,
                    max_negative,
                    max_positive,
                )
                if outline_value is not None
                else [0, 0, 0]
            )
        )
        temporal_payload[f"outline_fold_change_{idx}"] = outline_value
        temporal_payload[f"outline_color_{idx}"] = outline_color
        if idx == 1:
            fill_rgb = list(color)
            stroke_rgb = list(outline_color)
            fold_text = f"{chosen_value:.3f}"
            if headers:
                matched_uniprot = str(chosen_row.get(match_uniprot_key, "") or "").strip() if match_uniprot_key else ""
                if not matched_uniprot and display_uniprot_key:
                    matched_uniprot = str(chosen_row.get(display_uniprot_key, "") or "").strip()
                if display_gene_key:
                    matched_gene = str(chosen_row.get(display_gene_key, "") or "")
                _add_candidate_uniprot(matched_uniprot)

    tooltip_row = _select_dataset_tooltip_row(
        matches,
        dataset_index,
        matched_uniprot=matched_uniprot,
        matched_gene=matched_gene,
    )
    tooltip_payload = _build_dataset_tooltip(tooltip_row, tooltip_columns) if tooltip_row else {"tooltip": "", "tooltip_html": ""}
    title_text = tooltip_payload.get("tooltip") or ("Mapping: complex/process module" if is_complex else "")
    if matched_gene:
        display_label = matched_gene
    is_custom = bool(entry.get("isCustom") or entry.get("is_custom"))
    payload = {
        "id": str(entry.get("id") or ""),
        "originalLabel": original_label,
        "displayLabel": display_label,
        "label": mapping_label,
        "matchedGene": matched_gene,
        "matchedUniprot": matched_uniprot,
        "candidateUniprotIds": candidate_uniprot_ids,
        "suggestedGeneSymbols": list(mapping.get("suggested_gene_symbols") or []),
        "proteinOptions": protein_options,
        "foldText": fold_text,
        "hasDatasetMatch": bool(matches),
        "cx": float(entry.get("cx") or 0.0),
        "cy": float(entry.get("cy") or 0.0),
        "rx": max(4.0, float(entry.get("rx") or 12.0)),
        "ry": max(4.0, float(entry.get("ry") or 9.0)),
        "shapeType": str(entry.get("shapeType") or entry.get("shape_type") or "ellipse").strip().lower() or "ellipse",
        "angle": float(entry.get("angle") or 0.0),
        "strokeWidth": max(0.6, float(prot_outline_width or 1.0)),
        "stroke": f"rgb({stroke_rgb[0]}, {stroke_rgb[1]}, {stroke_rgb[2]})",
        "fillColor": f"rgb({fill_rgb[0]}, {fill_rgb[1]}, {fill_rgb[2]})",
        "opacity": max(0.1, min(1.0, float(entry.get("opacity") or 1.0))),
        "annotation": "",
        "annotationCommitted": True,
        "pendingAnnotation": "",
        "isCustom": is_custom,
        "isDrawingShape": False,
        "className": "cst-overlay-ellipse cst-complex-node" if is_complex else ("cst-overlay-ellipse" if matches else "cst-overlay-ellipse cst-missing-node"),
        "title": title_text,
        "tooltip": tooltip_payload.get("tooltip", ""),
        "tooltipHtml": tooltip_payload.get("tooltip_html", ""),
        "mappingType": mapping_type,
        "isComplex": is_complex,
        "complexDisplayMembers": complex_members,
    }
    if not is_complex:
        payload[f"outline_color_1"] = list(stroke_rgb)
    payload.update(temporal_payload)
    return payload


def _summarize_saved_editable_nodes(nodes: Sequence[Dict[str, Any]]) -> Dict[str, int]:
    summary = {
        "total_labels": len(list(nodes or [])),
        "recognized_total": 0,
        "psp_index_count": 0,
        "backup_count": 0,
        "unresolved_count": 0,
    }
    for node in list(nodes or []):
        mapping_type = str(node.get("mappingType") or "").strip().lower()
        if mapping_type == "text":
            continue
        if not mapping_type or mapping_type == "unresolved":
            summary["unresolved_count"] += 1
            continue
        summary["recognized_total"] += 1
        if mapping_type == "cst_index":
            summary["psp_index_count"] += 1
        else:
            summary["backup_count"] += 1
    return summary


def load_cst_overlay_state(file_path: Path, overlay_state_path: Optional[Path] = None) -> Dict[str, Any]:
    sidecar = Path(overlay_state_path) if overlay_state_path is not None else _cst_overlay_state_path(file_path)
    if not sidecar.exists():
        return {}
    try:
        return json.loads(sidecar.read_text(encoding="utf-8"))
    except Exception:
        return {}


def save_cst_overlay_state(
    file_path: Path,
    pathway_id: str,
    pathway_name: str,
    nodes: Sequence[Dict[str, Any]],
    edges: Sequence[Dict[str, Any]] | None = None,
    groups: Sequence[Dict[str, Any]] | None = None,
    disable_pdf_reader: bool = False,
    state_file_path: Optional[Path] = None,
) -> Dict[str, Any]:
    if state_file_path is None:
        raise ValueError("state_file_path is required for CST overlay persistence.")
    sanitized_nodes: List[Dict[str, Any]] = []
    for raw in list(nodes or []):
        if not isinstance(raw, dict):
            continue
        shape_type = str(raw.get("shapeType") or "ellipse").strip().lower() or "ellipse"
        if bool(raw.get("isDrawingShape")):
            sanitized_nodes.append(
                {
                    "id": str(raw.get("id") or ""),
                    "originalLabel": "",
                    "displayLabel": "",
                    "label": "",
                    "annotation": "",
                    "mappingType": "shape",
                    "suggested_gene_symbols": [],
                    "suggested_uniprot_ids": [],
                    "isCustom": True,
                    "userCreated": True,
                    "isDrawingShape": True,
                    "cx": float(raw.get("cx") or 0.0),
                    "cy": float(raw.get("cy") or 0.0),
                    "rx": max(10.0, float(raw.get("rx") or 60.0)),
                    "ry": max(10.0, float(raw.get("ry") or 40.0)),
                    "shapeType": shape_type,
                    "legendOrientation": str(raw.get("legendOrientation") or raw.get("legend_orientation") or "vertical"),
                    "angle": float(raw.get("angle") or 0.0),
                    "opacity": max(0.1, min(1.0, float(raw.get("opacity") or 1.0))),
                    "strokeWidth": max(1.0, float(raw.get("strokeWidth") or 1.0)),
                    "stroke": str(raw.get("stroke") or "#000000"),
                    "fillColor": str(raw.get("fillColor") or ("transparent" if shape_type == "bracket" else "#f5f5f5")),
                    "title": str(raw.get("title") or "Shape"),
                }
            )
            continue
        if shape_type == "text":
            text_label = str(raw.get("displayLabel") or raw.get("label") or "Text").strip() or "Text"
            sanitized_nodes.append(
                {
                    "id": str(raw.get("id") or ""),
                    "originalLabel": text_label,
                    "displayLabel": text_label,
                    "label": text_label,
                    "annotation": "",
                    "mappingType": "text",
                    "suggested_gene_symbols": [],
                    "suggested_uniprot_ids": [],
                    "isCustom": bool(raw.get("isCustom")),
                    "userCreated": bool(raw.get("userCreated")),
                    "isDrawingShape": False,
                    "cx": float(raw.get("cx") or 0.0),
                    "cy": float(raw.get("cy") or 0.0),
                    "rx": max(10.0, float(raw.get("rx") or 44.0)),
                    "ry": max(8.0, float(raw.get("ry") or 13.0)),
                    "shapeType": "text",
                    "angle": float(raw.get("angle") or 0.0),
                    "opacity": max(0.1, min(1.0, float(raw.get("opacity") or 1.0))),
                    "strokeWidth": max(1.0, float(raw.get("strokeWidth") or 1.5)),
                    "stroke": str(raw.get("stroke") or "rgb(71, 85, 105)"),
                    "fillColor": str(raw.get("fillColor") or "transparent"),
                    "textColor": str(raw.get("textColor") or "#0f172a"),
                    "fontSize": max(9.0, float(raw.get("fontSize") or 14.0)),
                    "fontWeight": str(raw.get("fontWeight") or "600"),
                    "fontFamily": str(raw.get("fontFamily") or '"Segoe UI", Arial, sans-serif'),
                    "textAlign": str(raw.get("textAlign") or "center"),
                    "title": str(raw.get("title") or "Text box"),
                }
            )
            continue
        annotation_text = str(raw.get("annotation") or "").strip()
        display_label = str(raw.get("displayLabel") or "").strip()
        label_text = str(raw.get("label") or "").strip()
        original_label = str(raw.get("originalLabel") or "").strip()
        if annotation_text:
            committed_name = annotation_text
        elif display_label and display_label not in {"Custom Ellipse", "Custom Rectangle", "Custom Bracket"}:
            committed_name = display_label
        elif label_text and label_text not in {"Custom Ellipse", "Custom Rectangle", "Custom Bracket"}:
            committed_name = label_text
        else:
            committed_name = original_label
        mapping = _resolve_cst_label_mapping(committed_name, file_path) if committed_name else {
            "mapping_type": "unresolved",
            "suggested_gene_symbols": [],
            "suggested_uniprot_ids": [],
        }
        sanitized_nodes.append(
                {
                    "id": str(raw.get("id") or ""),
                    "originalLabel": original_label or label_text,
                    "displayLabel": committed_name,
                    "label": committed_name,
                    "annotation": "",
                    "mappingType": str(mapping.get("mapping_type") or "").strip(),
                    "suggested_gene_symbols": list(mapping.get("suggested_gene_symbols") or []),
                    "suggested_uniprot_ids": list(mapping.get("suggested_uniprot_ids") or []),
                    "complex_components": [
                        dict(component)
                        for component in list(mapping.get("complex_components") or [])
                        if isinstance(component, dict)
                    ],
                    "isCustom": bool(raw.get("isCustom")),
                    "cx": float(raw.get("cx") or 0.0),
                    "cy": float(raw.get("cy") or 0.0),
                "rx": max(4.0, float(raw.get("rx") or 12.0)),
                "ry": max(4.0, float(raw.get("ry") or 9.0)),
                "shapeType": str(raw.get("shapeType") or "ellipse").strip().lower() or "ellipse",
                "angle": float(raw.get("angle") or 0.0),
                "opacity": max(0.1, min(1.0, float(raw.get("opacity") or 1.0))),
                "strokeWidth": max(1.0, float(raw.get("strokeWidth") or 4.5)),
            }
        )
    sanitized_edges: List[Dict[str, Any]] = []
    for raw in list(edges or []):
        if not isinstance(raw, dict):
            continue
        sanitized_edges.append(
            {
                "id": str(raw.get("id") or ""),
                "startX": float(raw.get("startX") or 0.0),
                "startY": float(raw.get("startY") or 0.0),
                "endX": float(raw.get("endX") or 0.0),
                "endY": float(raw.get("endY") or 0.0),
                "controlX": float(raw.get("controlX") or 0.0),
                "controlY": float(raw.get("controlY") or 0.0),
                "stroke": str(raw.get("stroke") or "#0f172a"),
                "strokeWidth": 1.6,
                "opacity": max(0.1, min(1.0, float(raw.get("opacity") or 0.95))),
                "userCreated": bool(raw.get("userCreated")),
                "isBent": bool(raw.get("isBent")),
                "dashed": bool(raw.get("dashed")),
                "startType": "arrow" if str(raw.get("startType") or "").strip().lower() == "arrow" else ("inhibitor" if str(raw.get("startType") or "").strip().lower() == "inhibitor" else "none"),
                "endType": "arrow" if str(raw.get("endType") or "").strip().lower() == "arrow" else ("inhibitor" if str(raw.get("endType") or "").strip().lower() == "inhibitor" else "none"),
                "bendPoints": [
                    {
                        "x": float(point.get("x") or 0.0),
                        "y": float(point.get("y") or 0.0),
                    }
                    for point in list(raw.get("bendPoints") or [])
                    if isinstance(point, dict)
                ],
            }
        )
    sanitized_groups: List[Dict[str, Any]] = []
    for raw in list(groups or []):
        if not isinstance(raw, dict):
            continue
        members: List[Dict[str, str]] = []
        seen_members: set[tuple[str, str]] = set()
        for member in list(raw.get("members") or []):
            if not isinstance(member, dict):
                continue
            member_type = str(member.get("type") or "").strip().lower()
            member_id = str(member.get("id") or "").strip()
            if member_type not in {"node", "edge"} or not member_id:
                continue
            key = (member_type, member_id)
            if key in seen_members:
                continue
            seen_members.add(key)
            members.append({"type": member_type, "id": member_id})
        if len(members) < 2:
            continue
        sanitized_groups.append(
            {
                "id": str(raw.get("id") or "").strip(),
                "members": members,
            }
        )
    payload = {
        "version": 4,
        "saved_at_utc": datetime.now(UTC).isoformat(),
        "pathway_id": str(pathway_id or "").strip(),
        "pathway_name": str(pathway_name or file_path.stem).strip(),
        "pathway_file": file_path.name,
        "nodes": sanitized_nodes,
        "edges": sanitized_edges,
        "groups": sanitized_groups,
        "disable_pdf_reader": bool(disable_pdf_reader),
    }
    sidecar = Path(state_file_path)
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return payload


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
    root = _resolve_cst_pathway_dir(base_dir)
    rows: List[Dict[str, str]] = []
    if not root.exists():
        return rows
    for file_path in sorted(root.rglob("*.ai"), key=lambda p: p.name.lower()):
        if not file_path.is_file():
            continue
        pathway_name = file_path.stem.strip()
        saved_state = load_cst_overlay_state(file_path)
        has_save_file = _cst_overlay_state_path(file_path).exists()
        disable_pdf_reader = bool(saved_state.get("disable_pdf_reader")) if isinstance(saved_state, dict) else False
        has_saved_edges = (bool(list(saved_state.get("edges") or [])) if isinstance(saved_state, dict) else False) or disable_pdf_reader
        has_saved_textboxes = False
        if isinstance(saved_state, dict):
            for node in list(saved_state.get("nodes") or []):
                if not isinstance(node, dict):
                    continue
                if str(node.get("shapeType") or node.get("shape_type") or "").strip().lower() == "text":
                    has_saved_textboxes = True
                    break
        rows.append(
            {
                "id": _cst_pathway_id_from_name(pathway_name),
                "name": pathway_name,
                "filename": file_path.name,
                "file_path": str(file_path),
                "source": "CST",
                "has_save_file": has_save_file,
                "has_saved_edges": has_saved_edges,
                "has_saved_textboxes": has_saved_textboxes,
                "disable_pdf_reader": disable_pdf_reader,
            }
        )
    return rows


def load_cst_pathway_payload(
    pathway_id: str,
    base_dir: Optional[Path] = None,
    protein_dataset: Optional[Dict[str, Any]] = None,
    metabolite_dataset: Optional[Dict[str, Any]] = None,
    ptm_dataset: Optional[Dict[str, Any]] = None,
    ptm_settings: Optional[Dict[str, Any]] = None,
    negative_color: Sequence[float] = _DEFAULT_NEGATIVE_COLOR,
    positive_color: Sequence[float] = _DEFAULT_POSITIVE_COLOR,
    max_negative: float = _DEFAULT_MAX_NEGATIVE,
    max_positive: float = _DEFAULT_MAX_POSITIVE,
    gradient_stops: Optional[List[Dict[str, Any]]] = None,
    prot_outline_width: float = 1.0,
    use_black_protein_outlines: bool = False,
    simple_kegg_mode: bool = True,
    temporal_mode: bool = False,
    overlay_state_path: Optional[Path] = None,
) -> Optional[Dict[str, Any]]:
    normalized_stops = _normalize_gradient_stops(gradient_stops)
    if len(normalized_stops) < 2:
        normalized_stops = _default_gradient_stops(negative_color, positive_color, max_negative, max_positive)
    _GRADIENT_STOPS_CTX.set(normalized_stops)
    pathway_key = str(pathway_id or "").strip().lower()
    if not pathway_key:
        return None
    for entry in get_cst_pathway_catalog(base_dir):
        if entry["id"].strip().lower() != pathway_key:
            continue
        file_path = Path(entry["file_path"])
        pdf_bytes = file_path.read_bytes()
        page_size = _extract_cst_page_size(str(file_path))
        saved_state = load_cst_overlay_state(file_path, overlay_state_path=overlay_state_path)
        disable_pdf_reader = bool(saved_state.get("disable_pdf_reader")) if isinstance(saved_state, dict) else False
        auto_ptm_payload = _build_cst_ptm_index(
            ptm_dataset,
            ptm_settings,
            negative_color,
            positive_color,
            max_negative,
            max_positive,
        )
        hidden_ptm_obstacle_edges: List[Dict[str, Any]] = []
        initial_editable_nodes: List[Dict[str, Any]] = []
        initial_editable_edges: List[Dict[str, Any]] = []
        initial_groups: List[Dict[str, Any]] = []
        overlay_nodes: List[Dict[str, Any]] = []
        mapping_summary: Dict[str, int] = {}
        has_saved_edges = isinstance(saved_state, dict) and "edges" in saved_state
        simple_kegg_mode = bool(simple_kegg_mode)
        if saved_state:
            if isinstance(saved_state.get("groups"), list):
                initial_groups = [dict(group) for group in list(saved_state.get("groups") or []) if isinstance(group, dict)]
            dataset_index = _build_dataset_index(protein_dataset)
            rebuilt_nodes: List[Dict[str, Any]] = []
            for saved_node in list(saved_state.get("nodes") or []):
                if simple_kegg_mode and str(saved_node.get("shapeType") or saved_node.get("shape_type") or "").strip().lower() == "text":
                    continue
                rebuilt = _build_editable_node_from_saved_entry(
                    file_path,
                    dict(saved_node or {}),
                    dataset_index,
                    auto_ptm_payload,
                    negative_color,
                    positive_color,
                    max_negative,
                    max_positive,
                    prot_outline_width=prot_outline_width,
                    use_black_protein_outlines=use_black_protein_outlines,
                )
                if rebuilt:
                    rebuilt_nodes.append(rebuilt)
            if rebuilt_nodes:
                initial_editable_nodes = rebuilt_nodes
                mapping_summary = _summarize_saved_editable_nodes(rebuilt_nodes)
            if has_saved_edges and not simple_kegg_mode:
                for saved_edge in list(saved_state.get("edges") or []):
                    if not isinstance(saved_edge, dict):
                        continue
                    saved_edge_id = str(saved_edge.get("id") or "")
                    parsed_edge = {
                        "id": saved_edge_id,
                        "startX": float(saved_edge.get("startX") or 0.0),
                        "startY": float(saved_edge.get("startY") or 0.0),
                        "endX": float(saved_edge.get("endX") or 0.0),
                        "endY": float(saved_edge.get("endY") or 0.0),
                        "controlX": float(saved_edge.get("controlX") or ((float(saved_edge.get("startX") or 0.0) + float(saved_edge.get("endX") or 0.0)) * 0.5)),
                        "controlY": float(saved_edge.get("controlY") or ((float(saved_edge.get("startY") or 0.0) + float(saved_edge.get("endY") or 0.0)) * 0.5)),
                        "stroke": "#0f172a",
                        "strokeWidth": 1.6,
                        "opacity": max(0.1, min(1.0, float(saved_edge.get("opacity") or 0.95))),
                        "userCreated": bool(saved_edge.get("userCreated")) if "userCreated" in saved_edge else (not saved_edge_id.startswith("cst_ai_edge_")),
                        "isBent": bool(saved_edge.get("isBent")) if "isBent" in saved_edge else None,
                        "dashed": bool(saved_edge.get("dashed")),
                        "startType": "arrow" if str(saved_edge.get("startType") or "").strip().lower() == "arrow" else ("inhibitor" if str(saved_edge.get("startType") or "").strip().lower() == "inhibitor" else "none"),
                        "endType": "arrow" if str(saved_edge.get("endType") or "").strip().lower() == "arrow" else ("inhibitor" if str(saved_edge.get("endType") or "").strip().lower() == "inhibitor" else "none"),
                        "bendPoints": [
                            {
                                "x": float(point.get("x") or 0.0),
                                "y": float(point.get("y") or 0.0),
                            }
                            for point in list(saved_edge.get("bendPoints") or [])
                            if isinstance(point, dict)
                        ],
                    }
                    if simple_kegg_mode:
                        hidden_ptm_obstacle_edges.append(dict(parsed_edge))
                    else:
                        initial_editable_edges.append(parsed_edge)
        if not disable_pdf_reader:
            extracted_edges = _extract_cst_edge_groups(str(file_path))
            if simple_kegg_mode:
                hidden_ptm_obstacle_edges = extracted_edges
            elif not has_saved_edges:
                initial_editable_edges = extracted_edges
        if not initial_editable_nodes:
            overlay_nodes = _build_cst_overlay_nodes(
                file_path,
                protein_dataset=protein_dataset,
                ptm_index=auto_ptm_payload,
                negative_color=negative_color,
                positive_color=positive_color,
                max_negative=max_negative,
                max_positive=max_positive,
                prot_outline_width=prot_outline_width,
                use_black_protein_outlines=use_black_protein_outlines,
            )
            mapping_summary = _summarize_cst_mapping_sources(file_path)
        search_catalog = _build_cst_search_catalog(
            protein_dataset,
            negative_color,
            positive_color,
            max_negative,
            max_positive,
            use_black_protein_outlines=use_black_protein_outlines,
        )
        metabolite_search_catalog = _build_cst_metabolite_search_catalog(
            metabolite_dataset,
            negative_color,
            positive_color,
            max_negative,
            max_positive,
        )
        data_uri = "data:application/pdf;base64," + base64.b64encode(pdf_bytes).decode("ascii")
        return {
            "id": entry["id"],
            "name": entry["name"],
            "filename": entry["filename"],
            "mime_type": "application/pdf",
            "data_uri": data_uri,
            "pdf_base64": base64.b64encode(pdf_bytes).decode("ascii"),
            "page_width": float(page_size.get("page_width") or _CST_PAGE_WIDTH),
            "page_height": float(page_size.get("page_height") or _CST_PAGE_HEIGHT),
            "ellipse_groups": _extract_cst_ellipse_groups(str(file_path)),
            "rect_groups": _extract_cst_rect_groups(str(file_path)),
            "mapping_summary": mapping_summary,
            "overlay_nodes": overlay_nodes,
            "initial_editable_nodes": initial_editable_nodes,
            "initial_editable_edges": initial_editable_edges,
            "initial_groups": initial_groups,
            "ptm_obstacle_edges": hidden_ptm_obstacle_edges,
            "auto_ptm_payload": auto_ptm_payload,
            "search_catalog": search_catalog,
            "metabolite_search_catalog": metabolite_search_catalog,
            "legend_config": {
                "posColor": _coerce_rgb(positive_color, _DEFAULT_POSITIVE_COLOR),
                "negColor": _coerce_rgb(negative_color, _DEFAULT_NEGATIVE_COLOR),
                "maxPos": float(normalized_stops[-1]["value"]) if normalized_stops else float(max_positive),
                "maxNeg": float(normalized_stops[0]["value"]) if normalized_stops else float(max_negative),
                "stops": normalized_stops,
            },
            "prot_outline_width": max(0.6, float(prot_outline_width or 1.0)),
            "use_black_protein_outlines": bool(use_black_protein_outlines),
            "disable_pdf_reader": disable_pdf_reader,
            "simple_kegg_mode": simple_kegg_mode,
            "temporal_mode": bool(temporal_mode),
        }
    return None


def create_cst_pathway_viewer(payload: Optional[Dict[str, Any]], save_input_id: str = "", export_key: str = "") -> Any:
    info = payload or {}
    data_uri = str(info.get("data_uri") or "").strip()
    pdf_base64 = str(info.get("pdf_base64") or "").strip()
    title = str(info.get("name") or info.get("filename") or "CST Pathway").strip()
    page_width = float(info.get("page_width") or _CST_PAGE_WIDTH)
    page_height = float(info.get("page_height") or _CST_PAGE_HEIGHT)
    active_idx = max(1, int(info.get("_active_fc_index") or 1))
    overlay_nodes = list(info.get("overlay_nodes") or [])
    initial_editable_nodes = list(info.get("initial_editable_nodes") or [])
    initial_editable_edges = list(info.get("initial_editable_edges") or [])
    initial_groups = list(info.get("initial_groups") or [])
    ptm_obstacle_edges = list(info.get("ptm_obstacle_edges") or [])
    disable_pdf_reader = bool(info.get("disable_pdf_reader"))
    simple_kegg_mode = bool(info.get("simple_kegg_mode", True))
    ellipse_groups = list(info.get("ellipse_groups") or [])
    rect_groups = list(info.get("rect_groups") or [])
    auto_ptm_payload = dict(info.get("auto_ptm_payload") or {})
    search_catalog = list(info.get("search_catalog") or [])
    metabolite_search_catalog = list(info.get("metabolite_search_catalog") or [])
    legend_config = dict(info.get("legend_config") or {})
    prot_outline_width = max(0.6, float(info.get("prot_outline_width") or 1.0))
    temporal_mode = bool(info.get("temporal_mode", False))
    mapping_summary = dict(info.get("mapping_summary") or {})
    overlay_nodes_json = _json_for_inline_script(overlay_nodes)
    initial_editable_nodes_json = _json_for_inline_script(initial_editable_nodes)
    initial_editable_edges_json = _json_for_inline_script(initial_editable_edges)
    initial_groups_json = _json_for_inline_script(initial_groups)
    ptm_obstacle_edges_json = _json_for_inline_script(ptm_obstacle_edges)
    overlay_signature = hashlib.md5(
        (overlay_nodes_json + "||" + initial_editable_nodes_json + "||" + initial_editable_edges_json + "||" + initial_groups_json).encode("utf-8")
    ).hexdigest()[:12]
    ellipse_groups_json = _json_for_inline_script(ellipse_groups)
    rect_groups_json = _json_for_inline_script(rect_groups)
    auto_ptm_payload_json = _json_for_inline_script(auto_ptm_payload)
    search_catalog_json = _json_for_inline_script(search_catalog)
    metabolite_search_catalog_json = _json_for_inline_script(metabolite_search_catalog)
    legend_config_json = _json_for_inline_script(legend_config)
    viewer_key = re.sub(r"[^A-Za-z0-9_-]+", "_", str(info.get("id") or "cst_viewer"))
    export_key = str(export_key or info.get("_bookmark_key") or viewer_key).strip().lower() or viewer_key
    stage_id = f"cst-stage-{viewer_key}"
    canvas_id = f"cst-canvas-{viewer_key}"
    fallback_id = f"cst-fallback-{viewer_key}"
    overlay_id = f"cst-overlay-{viewer_key}"
    viewport_id = f"cst-viewport-{viewer_key}"
    content_id = f"cst-content-{viewer_key}"
    missing_button_id = f"cst-missing-btn-{viewer_key}"
    add_button_id = f"cst-add-btn-{viewer_key}"
    add_arrow_button_id = f"cst-add-arrow-btn-{viewer_key}"
    dashed_arrow_button_id = f"cst-dashed-arrow-btn-{viewer_key}"
    inhibitor_arrow_button_id = f"cst-inhibitor-arrow-btn-{viewer_key}"
    both_arrow_button_id = f"cst-both-arrow-btn-{viewer_key}"
    line_arrow_button_id = f"cst-line-arrow-btn-{viewer_key}"
    add_rect_button_id = f"cst-add-rect-btn-{viewer_key}"
    add_bracket_button_id = f"cst-add-bracket-btn-{viewer_key}"
    convert_rect_button_id = f"cst-convert-rect-btn-{viewer_key}"
    batch_size_button_id = f"cst-batch-size-btn-{viewer_key}"
    batch_size_panel_id = f"cst-batch-size-panel-{viewer_key}"
    batch_size_value_id = f"cst-batch-size-value-{viewer_key}"
    batch_size_mode_id = f"cst-batch-size-mode-{viewer_key}"
    batch_size_apply_id = f"cst-batch-size-apply-{viewer_key}"
    batch_size_cancel_id = f"cst-batch-size-cancel-{viewer_key}"
    delete_button_id = f"cst-delete-btn-{viewer_key}"
    bring_front_button_id = f"cst-front-btn-{viewer_key}"
    send_back_button_id = f"cst-back-btn-{viewer_key}"
    protein_oval_button_id = f"cst-protein-oval-btn-{viewer_key}"
    auto_label_button_id = f"cst-auto-label-btn-{viewer_key}"
    text_auto_label_button_id = f"cst-text-auto-label-btn-{viewer_key}"
    auto_size_button_id = f"cst-auto-size-btn-{viewer_key}"
    undo_button_id = f"cst-undo-btn-{viewer_key}"
    redo_button_id = f"cst-redo-btn-{viewer_key}"
    save_button_id = f"cst-save-btn-{viewer_key}"
    disable_pdf_reader_id = f"cst-disable-pdf-reader-{viewer_key}"
    edge_resize_button_id = f"cst-edge-resize-btn-{viewer_key}"
    add_text_button_id = f"cst-add-text-btn-{viewer_key}"
    zoom_out_button_id = f"cst-zoom-out-{viewer_key}"
    zoom_in_button_id = f"cst-zoom-in-{viewer_key}"
    zoom_reset_button_id = f"cst-zoom-reset-{viewer_key}"
    coord_tooltip_id = f"cst-coords-{viewer_key}"
    hover_tooltip_id = f"cst-hover-tooltip-{viewer_key}"
    complex_menu_id = f"cst-complex-menu-{viewer_key}"
    complex_menu_button_id = f"cst-complex-menu-button-{viewer_key}"
    edge_bend_menu_id = f"cst-edge-bend-menu-{viewer_key}"
    edge_bend_reset_button_id = f"cst-edge-bend-reset-{viewer_key}"
    text_menu_id = f"cst-text-menu-{viewer_key}"
    text_menu_font_size_button_id = f"cst-text-menu-font-size-{viewer_key}"
    text_menu_alignment_button_id = f"cst-text-menu-alignment-{viewer_key}"
    text_menu_outline_button_id = f"cst-text-menu-outline-{viewer_key}"
    text_menu_bold_button_id = f"cst-text-menu-bold-{viewer_key}"
    text_menu_delete_button_id = f"cst-text-menu-delete-{viewer_key}"
    stats_id = f"cst-stats-{viewer_key}"
    opacity_wrap_id = f"cst-opacity-wrap-{viewer_key}"
    opacity_input_id = f"cst-opacity-input-{viewer_key}"
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
                overflow: hidden;
            }
            .cst-viewer-viewport {
                position: absolute;
                inset: 0;
                overflow: auto;
                background: #ffffff;
                scrollbar-gutter: stable;
            }
            .cst-viewer-content {
                position: relative;
                width: 100%;
                height: 100%;
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
            .cst-inline-text-editor {
                position: absolute;
                z-index: 5;
                display: none;
                margin: 0;
                padding: 2px 6px;
                border: 1px solid rgba(59, 130, 246, 0.8);
                border-radius: 4px;
                background: transparent;
                color: #0f172a;
                font-family: "Segoe UI", Arial, sans-serif;
                line-height: 1.1;
                outline: none;
                box-shadow: 0 0 0 1px rgba(147, 197, 253, 0.45);
                resize: none;
                overflow: hidden;
                white-space: pre-wrap;
            }
            .cst-inline-text-editor.is-open {
                display: block;
            }
            .cst-viewer-toolbar {
                position: absolute;
                top: 12px;
                right: 12px;
                z-index: 3;
                display: none;
                flex-direction: column;
                align-items: flex-end;
                gap: 8px;
            }
            .cst-viewer-toolbar-row {
                display: flex;
                align-items: center;
                justify-content: flex-end;
                gap: 8px;
                flex-wrap: wrap;
            }
            .cst-viewer-opacity {
                display: flex;
                align-items: center;
                gap: 8px;
                padding: 6px 10px;
                border-radius: 10px;
                border: 1px solid rgba(100, 116, 139, 0.28);
                background: rgba(255, 255, 255, 0.88);
                backdrop-filter: blur(8px);
                box-shadow: 0 8px 22px rgba(15, 23, 42, 0.08);
            }
            .cst-viewer-checkbox {
                display: flex;
                align-items: center;
                gap: 8px;
                padding: 6px 10px;
                border-radius: 10px;
                border: 1px solid rgba(100, 116, 139, 0.28);
                background: rgba(255, 255, 255, 0.88);
                backdrop-filter: blur(8px);
                box-shadow: 0 8px 22px rgba(15, 23, 42, 0.08);
                font-size: 12px;
                font-weight: 700;
                color: #334155;
                white-space: nowrap;
            }
            .cst-viewer-checkbox input[type="checkbox"] {
                margin: 0;
                accent-color: #2563eb;
            }
            .cst-viewer-opacity-label {
                font-size: 12px;
                font-weight: 700;
                color: #334155;
                white-space: nowrap;
            }
            .cst-viewer-opacity input[type="range"] {
                width: 110px;
                accent-color: #0f172a;
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
            .cst-viewer-tooltip {
                position: absolute;
                z-index: 7;
                display: none;
                pointer-events: none;
                background: rgba(0, 0, 0, 0.8);
                color: white;
                padding: 5px 10px;
                border-radius: 4px;
                font-size: 12px;
                line-height: 1.35;
                max-width: 300px;
                white-space: normal;
                word-wrap: break-word;
                transform: translate(6px, -36px);
            }
            .cst-complex-menu {
                position: absolute;
                z-index: 6;
                display: none;
                min-width: 140px;
                padding: 6px;
                border-radius: 12px;
                background: rgba(255, 255, 255, 0.98);
                border: 1px solid rgba(37, 99, 235, 0.18);
                box-shadow: 0 16px 30px rgba(15, 23, 42, 0.16);
                backdrop-filter: blur(10px);
            }
            .cst-complex-menu.is-open {
                display: block;
            }
            .cst-complex-menu button {
                width: 100%;
                border: 0;
                background: transparent;
                color: #1e3a8a;
                border-radius: 9px;
                padding: 8px 10px;
                text-align: left;
                font-size: 12px;
                font-weight: 700;
                line-height: 1.1;
            }
            .cst-complex-menu button:hover {
                background: rgba(37, 99, 235, 0.10);
            }
            .cst-edge-bend-menu {
                border-color: rgba(14, 116, 144, 0.22);
            }
            .cst-edge-bend-menu button {
                color: #155e75;
            }
            .cst-edge-bend-menu button:hover {
                background: rgba(8, 145, 178, 0.10);
            }
            .cst-text-menu {
                border-color: rgba(71, 85, 105, 0.22);
                min-width: 164px;
            }
            .cst-text-menu button {
                color: #334155;
            }
            .cst-text-menu button:hover {
                background: rgba(71, 85, 105, 0.10);
            }
            .cst-text-menu .cst-text-menu-item {
                position: relative;
            }
            .cst-text-menu .cst-text-menu-submenu {
                position: absolute;
                left: calc(100% + 4px);
                top: 0;
                display: none;
                min-width: 186px;
                padding: 8px;
                border-radius: 12px;
                background: rgba(255, 255, 255, 0.98);
                border: 1px solid rgba(71, 85, 105, 0.18);
                box-shadow: 0 16px 30px rgba(15, 23, 42, 0.16);
                backdrop-filter: blur(10px);
            }
            .cst-text-menu .cst-text-menu-item:hover > .cst-text-menu-submenu,
            .cst-text-menu .cst-text-menu-item:focus-within > .cst-text-menu-submenu {
                display: block;
            }
            .cst-text-menu .cst-text-menu-label {
                font-size: 11px;
                font-weight: 700;
                color: #475569;
                margin-bottom: 6px;
            }
            .cst-text-menu .cst-text-menu-row {
                display: flex;
                align-items: center;
                gap: 6px;
                margin-bottom: 8px;
            }
            .cst-text-menu .cst-text-menu-row:last-child {
                margin-bottom: 0;
            }
            .cst-text-menu .cst-text-menu-row input[type="text"],
            .cst-text-menu .cst-text-menu-row input[type="number"],
            .cst-text-menu .cst-text-menu-row select {
                flex: 1 1 auto;
                min-width: 0;
                border: 1px solid rgba(148, 163, 184, 0.7);
                border-radius: 8px;
                padding: 6px 8px;
                font-size: 12px;
                color: #0f172a;
                background: #fff;
            }
            .cst-text-menu .cst-text-menu-row input[type="color"] {
                width: 34px;
                height: 28px;
                padding: 0;
                border: 1px solid rgba(148, 163, 184, 0.7);
                border-radius: 8px;
                background: #fff;
            }
            .cst-text-menu .cst-text-menu-option-grid {
                display: grid;
                grid-template-columns: repeat(3, minmax(0, 1fr));
                gap: 6px;
            }
            .cst-text-menu .cst-text-menu-option-grid button,
            .cst-text-menu .cst-text-menu-toggle {
                border: 1px solid rgba(148, 163, 184, 0.5);
                background: rgba(248, 250, 252, 0.98);
                color: #334155;
                border-radius: 8px;
                padding: 6px 8px;
                font-size: 12px;
                font-weight: 700;
                text-align: center;
            }
            .cst-text-menu .cst-text-menu-option-grid button.is-active,
            .cst-text-menu .cst-text-menu-toggle.is-active {
                background: rgba(71, 85, 105, 0.12);
                border-color: rgba(71, 85, 105, 0.45);
                color: #0f172a;
            }
            .cst-viewer-overlay-svg {
                width: 100%;
                height: 100%;
                display: block;
                overflow: visible;
                pointer-events: auto;
            }
            .cst-viewer-overlay-svg ellipse,
            .cst-viewer-overlay-svg rect {
                pointer-events: visiblePainted;
                cursor: help;
            }
            .cst-viewer-overlay-svg .cst-edit-node-ellipse {
                pointer-events: visiblePainted;
                cursor: move;
            }
            .cst-viewer-overlay-svg .cst-edit-node-ellipse.cst-node-selected {
                filter: drop-shadow(0 0 6px rgba(15, 23, 42, 0.18));
            }
            .cst-viewer-overlay-svg .cst-complex-node {
                fill: #e5e7eb;
                stroke: #2563eb;
            }
            .cst-viewer-overlay-svg .cst-protein-oval-label {
                pointer-events: none;
                fill: #0f172a;
                font-family: "Segoe UI", Arial, sans-serif;
                font-weight: 700;
                text-anchor: middle;
                dominant-baseline: middle;
            }
            .cst-viewer-overlay-svg .cst-protein-oval-label.cst-custom-annotation-label {
                fill: #c1121f;
            }
            .cst-viewer-overlay-svg .cst-protein-oval-label.cst-pending-annotation-label {
                fill: #2563eb;
            }
            .cst-viewer-overlay-svg .cst-auto-ptm-circle {
                pointer-events: all;
                cursor: move;
            }
            .cst-viewer-overlay-svg .cst-auto-ptm-symbol,
            .cst-viewer-overlay-svg .cst-auto-ptm-site-label {
                font-family: "Segoe UI", Arial, sans-serif;
                font-weight: 800;
                pointer-events: all;
                cursor: move;
            }
            .cst-viewer-overlay-svg .cst-auto-ptm-symbol {
                text-anchor: middle;
                dominant-baseline: middle;
            }
            .cst-viewer-overlay-svg .cst-selected-outline {
                pointer-events: none;
                stroke: rgba(15, 23, 42, 0.9);
                stroke-dasharray: 8 6;
                fill: none;
            }
            .cst-viewer-overlay-svg .cst-complex-panel {
                filter: drop-shadow(0 14px 24px rgba(15, 23, 42, 0.14));
            }
            .cst-viewer-overlay-svg .cst-complex-panel-body {
                fill: rgba(255, 255, 255, 0.97);
                stroke: rgba(37, 99, 235, 0.32);
                stroke-width: 1.2;
            }
            .cst-viewer-overlay-svg .cst-complex-panel-title {
                fill: #1e293b;
                font-family: "Segoe UI", Arial, sans-serif;
                font-size: 11px;
                font-weight: 800;
            }
            .cst-viewer-overlay-svg .cst-complex-panel-close {
                fill: rgba(239, 68, 68, 0.10);
                stroke: rgba(239, 68, 68, 0.42);
                stroke-width: 1;
                cursor: pointer;
            }
            .cst-viewer-overlay-svg .cst-complex-panel-close-text {
                fill: #b91c1c;
                font-family: "Segoe UI", Arial, sans-serif;
                font-size: 11px;
                font-weight: 800;
                text-anchor: middle;
                dominant-baseline: middle;
                pointer-events: none;
            }
            .cst-viewer-overlay-svg .cst-complex-member-rect {
                stroke: rgba(148, 163, 184, 0.65);
                stroke-width: 1;
            }
            .cst-viewer-overlay-svg .cst-complex-member-rect.cst-complex-member-matched {
                stroke: rgba(15, 23, 42, 0.24);
            }
            .cst-viewer-overlay-svg .cst-complex-member-label {
                fill: #0f172a;
                font-family: "Segoe UI", Arial, sans-serif;
                font-size: 10.5px;
                font-weight: 700;
                text-anchor: middle;
                dominant-baseline: middle;
                pointer-events: none;
            }
            .cst-viewer-overlay-svg .cst-handle-line {
                pointer-events: none;
                stroke: rgba(15, 23, 42, 0.5);
                stroke-width: 1.6;
                stroke-dasharray: 4 4;
            }
            .cst-viewer-overlay-svg .cst-edit-handle {
                pointer-events: all;
                fill: transparent;
                stroke: transparent;
                stroke-width: 0;
                opacity: 0;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-handle="x"] {
                cursor: ew-resize;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-handle="y"] {
                cursor: ns-resize;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-handle="left"] {
                cursor: ew-resize;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-handle="right"] {
                cursor: ew-resize;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-handle="top"] {
                cursor: ns-resize;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-handle="bottom"] {
                cursor: ns-resize;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-handle="rotate"] {
                cursor: grab;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-role="edge-start"],
            .cst-viewer-overlay-svg .cst-edit-handle[data-role="edge-end"] {
                fill: #dc2626;
                stroke: #ffffff;
                stroke-width: 1.1;
                opacity: 0.98;
                cursor: move;
            }
            .cst-viewer-overlay-svg .cst-edit-handle[data-role="edge-bend"] {
                fill: #2563eb;
                stroke: #ffffff;
                stroke-width: 1.1;
                opacity: 0.98;
                cursor: move;
            }
            .cst-viewer-overlay-svg .cst-edit-hitbox {
                pointer-events: all;
                fill: transparent;
                stroke: transparent;
            }
            .cst-viewer-overlay-svg .cst-missing-node {
                display: none;
            }
            .cst-viewer-stage.cst-show-missing .cst-viewer-overlay-svg .cst-missing-node {
                display: inline;
            }
            .cst-viewer-stage.cst-add-mode .cst-viewer-overlay-svg {
                cursor: crosshair;
            }
            .cst-viewer-action.is-active {
                background: rgba(37, 99, 235, 0.16);
                border-color: rgba(37, 99, 235, 0.42);
                color: #1d4ed8;
            }
            .cst-viewer-action.is-disabled {
                opacity: 0.45;
                cursor: default;
            }
            .cst-batch-size-panel {
                position: absolute;
                right: 12px;
                top: 58px;
                z-index: 5;
                display: none;
                min-width: 220px;
                padding: 12px;
                border-radius: 14px;
                background: rgba(255, 255, 255, 0.96);
                border: 1px solid rgba(148, 163, 184, 0.28);
                box-shadow: 0 16px 32px rgba(15, 23, 42, 0.14);
                backdrop-filter: blur(10px);
                gap: 10px;
            }
            .cst-batch-size-panel.is-open {
                display: flex;
                flex-direction: column;
            }
            .cst-batch-size-title {
                font-size: 11px;
                font-weight: 800;
                letter-spacing: 0.06em;
                text-transform: uppercase;
                color: #475569;
            }
            .cst-batch-size-row {
                display: flex;
                align-items: center;
                gap: 8px;
            }
            .cst-batch-size-row input,
            .cst-batch-size-row select {
                width: 100%;
                border: 1px solid rgba(148, 163, 184, 0.5);
                border-radius: 10px;
                padding: 7px 9px;
                font-size: 12px;
                color: #0f172a;
                background: #fff;
            }
            .cst-batch-size-actions {
                display: flex;
                justify-content: flex-end;
                gap: 8px;
            }
            """
        ),
        ui.tags.script({"src": "https://cdnjs.cloudflare.com/ajax/libs/pdf.js/3.11.174/pdf.min.js"}),
        ui.div(
            {"class": "cst-viewer-shell", "data-cst-viewer-key": viewer_key},
            ui.div(
                {"class": "cst-viewer-stage", "id": stage_id},
                ui.div(
                    {"class": "cst-viewer-toolbar"},
                    ui.div(
                        {"class": "cst-viewer-toolbar-row"},
                    ui.div(
                        {"class": "cst-viewer-opacity", "id": opacity_wrap_id, "title": "Change the opacity of all CST overlay ellipses."},
                        ui.div({"class": "cst-viewer-opacity-label"}, "Opacity"),
                        ui.tags.input(
                            {
                                "id": opacity_input_id,
                                "type": "range",
                                "min": "0.10",
                                "max": "1.00",
                                "step": "0.05",
                                "value": "1.00",
                            }
                        ),
                    ),
                    ui.tags.button(
                        {
                            "id": undo_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Undo the last CST ellipse change (Ctrl+Z).",
                        },
                        "Undo",
                    ),
                    ui.tags.button(
                        {
                            "id": zoom_out_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Zoom out",
                        },
                        "-",
                    ),
                    ui.tags.button(
                        {
                            "id": zoom_reset_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Reset zoom",
                        },
                        "100%",
                    ),
                    ui.tags.button(
                        {
                            "id": zoom_in_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Zoom in",
                        },
                        "+",
                    ),
                    ui.tags.button(
                        {
                            "id": protein_oval_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Switch between outline rings and filled protein ovals.",
                        },
                        "Protein Ovals",
                    ),
                    ui.tags.button(
                        {
                            "id": edge_resize_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Switch to one-sided resize handles so each side of the oval can move independently.",
                        },
                        "Edge Resize",
                    ),
                    ui.tags.button(
                        {
                            "id": add_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Add a new editable ellipse to the CST overlay.",
                        },
                        "Add Ellipse",
                    ),
                    ui.tags.button(
                        {
                            "id": add_arrow_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Add a new curved arrow to the CST overlay.",
                        },
                        "Add Arrow",
                    ),
                    ui.tags.button(
                        {
                            "id": add_text_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Add a new editable text box annotation to the CST overlay.",
                        },
                        "Add Text",
                    ),
                    ),
                    ui.div(
                        {"class": "cst-viewer-toolbar-row"},
                    ui.tags.button(
                        {
                            "id": dashed_arrow_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Toggle the selected arrow between solid and dashed.",
                        },
                        "Dashed",
                    ),
                    ui.tags.button(
                        {
                            "id": both_arrow_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Toggle the selected line between a single arrow and arrowheads on both ends.",
                        },
                        "Double Arrow",
                    ),
                    ui.tags.button(
                        {
                            "id": line_arrow_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Toggle the selected edge between a line and an arrow.",
                        },
                        "Line",
                    ),
                    ui.tags.button(
                        {
                            "id": inhibitor_arrow_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Toggle the selected arrow between an arrow head and inhibitor end.",
                        },
                        "Inhibitor",
                    ),
                    ui.tags.button(
                        {
                            "id": add_rect_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Add a new editable rectangle protein module to the CST overlay.",
                        },
                        "Add Rectangle",
                    ),
                    ui.tags.button(
                        {
                            "id": add_bracket_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Add a new editable bracket shape to the CST overlay.",
                        },
                        "Add Bracket",
                    ),
                    ui.tags.button(
                        {
                            "id": convert_rect_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Turn the selected protein modules into rectangles.",
                        },
                        "Make Rectangle",
                    ),
                    ui.tags.button(
                        {
                            "id": batch_size_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Change the width, height, or both for all selected protein modules.",
                        },
                        "Batch Size",
                    ),
                    ui.tags.button(
                        {
                            "id": delete_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Delete the selected ellipse.",
                        },
                        "Delete Ellipse",
                    ),
                    ui.tags.button(
                        {
                            "id": send_back_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Send the selected protein module behind the other protein modules.",
                        },
                        "Send To Back",
                    ),
                    ui.tags.button(
                        {
                            "id": bring_front_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                            "title": "Bring the selected protein module in front of the other protein modules.",
                        },
                        "Bring To Front",
                    ),
                    ui.tags.button(
                        {
                            "id": missing_button_id,
                            "class": "cst-viewer-action",
                            "type": "button",
                        },
                        "Show Missing Proteins",
                    ),
                    ),
                ),
                ui.div(
                    {"class": "cst-batch-size-panel", "id": batch_size_panel_id},
                    ui.div({"class": "cst-batch-size-title"}, "Batch Size"),
                    ui.div(
                        {"class": "cst-batch-size-row"},
                        ui.tags.input(
                            {
                                "id": batch_size_value_id,
                                "type": "number",
                                "min": "1",
                                "step": "1",
                                "placeholder": "Units",
                            }
                        ),
                    ),
                    ui.div(
                        {"class": "cst-batch-size-row"},
                        ui.tags.select(
                            {
                                "id": batch_size_mode_id,
                            },
                            ui.tags.option({"value": "both"}, "Width and Height"),
                            ui.tags.option({"value": "width"}, "Width Only"),
                            ui.tags.option({"value": "height"}, "Height Only"),
                        ),
                    ),
                    ui.div(
                        {"class": "cst-batch-size-actions"},
                        ui.tags.button(
                            {
                                "id": batch_size_cancel_id,
                                "class": "cst-viewer-action",
                                "type": "button",
                            },
                            "Cancel",
                        ),
                        ui.tags.button(
                            {
                                "id": batch_size_apply_id,
                                "class": "cst-viewer-action",
                                "type": "button",
                            },
                            "Apply",
                        ),
                    ),
                ),
                ui.div({"class": "cst-viewer-coords", "id": coord_tooltip_id}),
                ui.div({"class": "cst-viewer-tooltip", "id": hover_tooltip_id}),
                ui.div(
                    {"class": "cst-complex-menu", "id": complex_menu_id},
                    ui.tags.button(
                        {
                            "id": complex_menu_button_id,
                            "type": "button",
                        },
                        "Show Complex",
                    ),
                ),
                ui.div(
                    {"class": "cst-complex-menu cst-edge-bend-menu", "id": edge_bend_menu_id},
                    ui.tags.button(
                        {
                            "id": edge_bend_reset_button_id,
                            "type": "button",
                        },
                        "Reset",
                    ),
                ),
                ui.div(
                    {"class": "cst-complex-menu cst-text-menu", "id": text_menu_id},
                ),
                ui.div(
                    {"class": "cst-viewer-viewport", "id": viewport_id},
                    ui.div(
                        {"class": "cst-viewer-content", "id": content_id},
                        ui.tags.canvas({"class": "cst-viewer-canvas", "id": canvas_id}),
                        ui.tags.iframe(
                            {
                                "class": "cst-viewer-fallback",
                                "id": fallback_id,
                                "src": f"{data_uri}#toolbar=0&navpanes=0&view=FitH",
                                "title": title,
                            }
                        ),
                        ui.tags.textarea({"class": "cst-inline-text-editor", "rows": "1"}),
                        ui.div({"class": "cst-viewer-overlay"}, overlay_markup),
                    ),
                ),
            ),
        ),
        ui.tags.script(
            f"""
            (function() {{
                const stage = document.getElementById('{stage_id}');
                const canvas = document.getElementById('{canvas_id}');
                const fallback = document.getElementById('{fallback_id}');
                const overlaySvg = document.getElementById('{overlay_id}');
                const viewport = document.getElementById('{viewport_id}');
                const content = document.getElementById('{content_id}');
                const missingButton = document.getElementById('{missing_button_id}');
                const addButton = document.getElementById('{add_button_id}');
                const addArrowButton = document.getElementById('{add_arrow_button_id}');
                const dashedArrowButton = document.getElementById('{dashed_arrow_button_id}');
                const bothArrowButton = document.getElementById('{both_arrow_button_id}');
                const lineArrowButton = document.getElementById('{line_arrow_button_id}');
                const inhibitorArrowButton = document.getElementById('{inhibitor_arrow_button_id}');
                const addTextButton = document.getElementById('{add_text_button_id}');
                const addRectButton = document.getElementById('{add_rect_button_id}');
                const addBracketButton = document.getElementById('{add_bracket_button_id}');
                const convertRectButton = document.getElementById('{convert_rect_button_id}');
                const batchSizeButton = document.getElementById('{batch_size_button_id}');
                const batchSizePanel = document.getElementById('{batch_size_panel_id}');
                const batchSizeValue = document.getElementById('{batch_size_value_id}');
                const batchSizeMode = document.getElementById('{batch_size_mode_id}');
                const batchSizeApply = document.getElementById('{batch_size_apply_id}');
                const batchSizeCancel = document.getElementById('{batch_size_cancel_id}');
                const deleteButton = document.getElementById('{delete_button_id}');
                const bringFrontButton = document.getElementById('{bring_front_button_id}');
                const sendBackButton = document.getElementById('{send_back_button_id}');
                const proteinOvalButton = document.getElementById('{protein_oval_button_id}');
                const autoLabelButton = document.getElementById('{auto_label_button_id}');
                const textAutoLabelButton = document.getElementById('{text_auto_label_button_id}');
                const autoSizeButton = document.getElementById('{auto_size_button_id}');
                const edgeResizeButton = document.getElementById('{edge_resize_button_id}');
                const undoButton = document.getElementById('{undo_button_id}');
                const redoButton = document.getElementById('{redo_button_id}');
                const saveButton = document.getElementById('{save_button_id}');
                const disablePdfReaderCheckbox = document.getElementById('{disable_pdf_reader_id}');
                const zoomOutButton = document.getElementById('{zoom_out_button_id}');
                const zoomInButton = document.getElementById('{zoom_in_button_id}');
                const zoomResetButton = document.getElementById('{zoom_reset_button_id}');
                const opacityInput = document.getElementById('{opacity_input_id}');
                const coordsTooltip = document.getElementById('{coord_tooltip_id}');
                const hoverTooltip = document.getElementById('{hover_tooltip_id}');
                const complexMenu = document.getElementById('{complex_menu_id}');
                const complexMenuButton = document.getElementById('{complex_menu_button_id}');
                const edgeBendMenu = document.getElementById('{edge_bend_menu_id}');
                const edgeBendResetButton = document.getElementById('{edge_bend_reset_button_id}');
                const textMenu = document.getElementById('{text_menu_id}');
                const textMenuFontSizeButton = document.getElementById('{text_menu_font_size_button_id}');
                const textMenuAlignmentButton = document.getElementById('{text_menu_alignment_button_id}');
                const textMenuOutlineButton = document.getElementById('{text_menu_outline_button_id}');
                const textMenuBoldButton = document.getElementById('{text_menu_bold_button_id}');
                const textMenuDeleteButton = document.getElementById('{text_menu_delete_button_id}');
                const inlineTextEditor = content ? content.querySelector('.cst-inline-text-editor') : null;
                if (!stage || !canvas || !viewport || !content) return;
                const renderStateKey = '__mapkinaseCstRenderState';
                const pdfData = '{pdf_base64}';
                const pageWidth = {page_width:.6f};
                const pageHeight = {page_height:.6f};
                const overlayVersion = '{overlay_signature}_{active_idx}';
                const overlayData = {overlay_nodes_json};
                const initialEditableNodes = {initial_editable_nodes_json};
                const initialEditableEdges = {initial_editable_edges_json};
                const initialGroups = {initial_groups_json};
                const ellipseGroups = {ellipse_groups_json};
                const rectGroups = {rect_groups_json};
                const saveInputId = {_json_for_inline_script(save_input_id)};
                const svgNs = 'http://www.w3.org/2000/svg';
                const createSvg = (name) => document.createElementNS(svgNs, name);
                const pointerToViewer = (evt) => {{
                    const rect = viewport.getBoundingClientRect();
                    if (!rect.width || !rect.height) return null;
                    const localX = (evt.clientX - rect.left) + viewport.scrollLeft;
                    const localY = (evt.clientY - rect.top) + viewport.scrollTop;
                    const contentWidth = Math.max(Number(content.clientWidth || 0), 1);
                    const contentHeight = Math.max(Number(content.clientHeight || 0), 1);
                    return {{
                        overlayX: (localX / contentWidth) * pageWidth,
                        overlayY: (localY / contentHeight) * pageHeight,
                        pdfX: (localX / contentWidth) * pageWidth,
                        pdfY: pageHeight - ((localY / contentHeight) * pageHeight),
                        localX,
                        localY,
                    }};
                }};
                const getViewportCenter = () => {{
                    const contentWidth = Math.max(Number(content.clientWidth || 0), 1);
                    const contentHeight = Math.max(Number(content.clientHeight || 0), 1);
                    const localX = Number(viewport.scrollLeft || 0) + (Number(viewport.clientWidth || 0) * 0.5);
                    const localY = Number(viewport.scrollTop || 0) + (Number(viewport.clientHeight || 0) * 0.5);
                    return {{
                        overlayX: (localX / contentWidth) * pageWidth,
                        overlayY: (localY / contentHeight) * pageHeight,
                    }};
                }};
                const toViewportPosition = (evt) => {{
                    const rect = viewport.getBoundingClientRect();
                    return {{
                        x: (evt.clientX - rect.left) + viewport.scrollLeft,
                        y: (evt.clientY - rect.top) + viewport.scrollTop,
                    }};
                }};
                const rotatePoint = (x, y, angleDeg) => {{
                    const angle = (Number(angleDeg || 0) * Math.PI) / 180;
                    const cos = Math.cos(angle);
                    const sin = Math.sin(angle);
                    return {{
                        x: (x * cos) - (y * sin),
                        y: (x * sin) + (y * cos),
                    }};
                }};
                const toLocalAxes = (point, node) => {{
                    const dx = Number(point.overlayX || 0) - Number(node.cx || 0);
                    const dy = Number(point.overlayY || 0) - Number(node.cy || 0);
                    return rotatePoint(dx, dy, -Number(node.angle || 0));
                }};
                const addNodeId = (() => {{
                    let counter = 1;
                    return () => 'cst_overlay_' + Date.now().toString(36) + '_' + (counter++);
                }})();
                const cloneStateSnapshot = (snapshot) => JSON.parse(JSON.stringify(snapshot));
                const getDefaultBendPointFromQuadratic = (edge) => {{
                    const startX = Number(edge && edge.startX || 0);
                    const startY = Number(edge && edge.startY || 0);
                    const controlX = Number(edge && edge.controlX || ((startX + Number(edge && edge.endX || 0)) * 0.5));
                    const controlY = Number(edge && edge.controlY || ((startY + Number(edge && edge.endY || 0)) * 0.5));
                    const endX = Number(edge && edge.endX || 0);
                    const endY = Number(edge && edge.endY || 0);
                    return {{
                        x: (0.25 * startX) + (0.5 * controlX) + (0.25 * endX),
                        y: (0.25 * startY) + (0.5 * controlY) + (0.25 * endY),
                    }};
                }};
                const getStraightEdgeMidpoint = (edge) => {{
                    const startX = Number(edge && edge.startX || 0);
                    const startY = Number(edge && edge.startY || 0);
                    const endX = Number(edge && edge.endX || 0);
                    const endY = Number(edge && edge.endY || 0);
                    return {{
                        x: (startX + endX) * 0.5,
                        y: (startY + endY) * 0.5,
                    }};
                }};
                const isEdgeBent = (edge) => {{
                    if (!edge) return false;
                    if (typeof edge.isBent === 'boolean') return edge.isBent;
                    const straightMid = getStraightEdgeMidpoint(edge);
                    const hasStoredBends = Array.isArray(edge.bendPoints) && edge.bendPoints.length;
                    const probe = hasStoredBends
                        ? edge.bendPoints[Math.max(0, Math.min(edge.bendPoints.length - 1, Math.floor(edge.bendPoints.length * 0.5)))]
                        : getDefaultBendPointFromQuadratic(edge);
                    const dist = Math.hypot(
                        Number((probe && probe.x) || 0) - Number(straightMid.x || 0),
                        Number((probe && probe.y) || 0) - Number(straightMid.y || 0)
                    );
                    edge.isBent = dist > 1.5;
                    return edge.isBent;
                }};
                const getActiveEdgeMidpoint = (edge) => {{
                    if (!isEdgeBent(edge)) {{
                        const midpoint = getStraightEdgeMidpoint(edge);
                        edge.controlX = Number(midpoint.x || 0);
                        edge.controlY = Number(midpoint.y || 0);
                        return midpoint;
                    }}
                    if (!Array.isArray(edge.bendPoints) || !edge.bendPoints.length) {{
                        edge.bendPoints = [getDefaultBendPointFromQuadratic(edge)];
                    }}
                    edge.bendPoints = edge.bendPoints.map((point) => ({{
                        x: Number(point && point.x || 0),
                        y: Number(point && point.y || 0),
                    }}));
                    const midpointIndex = Math.floor(edge.bendPoints.length * 0.5);
                    const midpoint = edge.bendPoints[Math.max(0, Math.min(edge.bendPoints.length - 1, midpointIndex))];
                    edge.controlX = Number(midpoint && midpoint.x || edge.controlX || 0);
                    edge.controlY = Number(midpoint && midpoint.y || edge.controlY || 0);
                    return midpoint;
                }};
                const normalizeEditableNodeIds = (nodes) => {{
                    const seen = new Set();
                    for (const node of (Array.isArray(nodes) ? nodes : [])) {{
                        const rawId = String(node && node.id || '').trim();
                        if (!rawId || seen.has(rawId)) {{
                            node.id = addNodeId();
                        }} else {{
                            node.id = rawId;
                        }}
                        seen.add(node.id);
                    }}
                    return Array.isArray(nodes) ? nodes : [];
                }};
                const normalizeEditableEdgeIds = (edges) => {{
                    const seen = new Set();
                    for (const edge of (Array.isArray(edges) ? edges : [])) {{
                        const rawId = String(edge && edge.id || '').trim();
                        if (!rawId || seen.has(rawId)) {{
                            edge.id = 'cst_edge_' + Date.now() + '_' + Math.floor(Math.random() * 1000000);
                        }} else {{
                            edge.id = rawId;
                        }}
                        if (!Array.isArray(edge.bendPoints)) {{
                            edge.bendPoints = [];
                        }} else {{
                            edge.bendPoints = edge.bendPoints.map((point) => ({{
                                x: Number(point && point.x || 0),
                                y: Number(point && point.y || 0),
                            }}));
                        }}
                        isEdgeBent(edge);
                        const midpoint = getActiveEdgeMidpoint(edge);
                        edge.controlX = Number(midpoint && midpoint.x || edge.controlX || 0);
                        edge.controlY = Number(midpoint && midpoint.y || edge.controlY || 0);
                        seen.add(edge.id);
                    }}
                    return Array.isArray(edges) ? edges : [];
                }};
                const buildSnapshot = (local) => cloneStateSnapshot({{
                    editableNodes: local.editableNodes || [],
                    editableEdges: local.editableEdges || [],
                    groups: Array.isArray(local.groups) ? local.groups : [],
                    autoPtmPlacements: Array.isArray(local.autoPtmPlacements) ? local.autoPtmPlacements : null,
                    selectedNodeId: local.selectedNodeId || null,
                    selectedNodeIds: Array.isArray(local.selectedNodeIds) ? local.selectedNodeIds : [],
                    selectedEdgeId: local.selectedEdgeId || null,
                    selectedEdgeIds: Array.isArray(local.selectedEdgeIds) ? local.selectedEdgeIds : [],
                    selectedPtmId: local.selectedPtmId || null,
                    selectedPtmNodeId: local.selectedPtmNodeId || null,
                    globalOpacity: Number(local.globalOpacity || 1.0),
                    proteinOvalMode: !!local.proteinOvalMode,
                    edgeResizeMode: !!local.edgeResizeMode,
                }});
                const applySnapshot = (local, snapshot) => {{
                    const safe = cloneStateSnapshot(snapshot || {{}});
                    local.editableNodes = normalizeEditableNodeIds(Array.isArray(safe.editableNodes) ? safe.editableNodes : []);
                    local.editableEdges = normalizeEditableEdgeIds(Array.isArray(safe.editableEdges) ? safe.editableEdges : []);
                    local.groups = sanitizeGroups(local, Array.isArray(safe.groups) ? safe.groups : []);
                    local.autoPtmPlacements = Array.isArray(safe.autoPtmPlacements) ? cloneStateSnapshot(safe.autoPtmPlacements) : null;
                    local.selectedNodeId = safe.selectedNodeId || null;
                    local.selectedNodeIds = Array.isArray(safe.selectedNodeIds) ? safe.selectedNodeIds : (local.selectedNodeId ? [local.selectedNodeId] : []);
                    local.selectedEdgeId = safe.selectedEdgeId || null;
                    local.selectedEdgeIds = Array.isArray(safe.selectedEdgeIds) ? safe.selectedEdgeIds : (local.selectedEdgeId ? [local.selectedEdgeId] : []);
                    local.selectedPtmId = safe.selectedPtmId || null;
                    local.selectedPtmNodeId = safe.selectedPtmNodeId || null;
                    if (local.selectedNodeId && !local.editableNodes.some((node) => node.id === local.selectedNodeId)) {{
                        local.selectedNodeId = null;
                    }}
                    local.selectedNodeIds = Array.from(new Set((local.selectedNodeIds || []).filter((nodeId) => local.editableNodes.some((node) => node.id === nodeId))));
                    if (local.selectedNodeId && !local.selectedNodeIds.includes(local.selectedNodeId)) {{
                        local.selectedNodeIds.push(local.selectedNodeId);
                    }}
                    if (!local.selectedNodeId && local.selectedNodeIds.length) {{
                        local.selectedNodeId = local.selectedNodeIds[local.selectedNodeIds.length - 1];
                    }}
                    if (local.selectedEdgeId && !local.editableEdges.some((edge) => edge.id === local.selectedEdgeId)) {{
                        local.selectedEdgeId = null;
                    }}
                    local.selectedEdgeIds = Array.from(new Set((local.selectedEdgeIds || []).filter((edgeId) => local.editableEdges.some((edge) => edge.id === edgeId))));
                    if (local.selectedEdgeId && !local.selectedEdgeIds.includes(local.selectedEdgeId)) {{
                        local.selectedEdgeIds.push(local.selectedEdgeId);
                    }}
                    if (!local.selectedEdgeId && local.selectedEdgeIds.length) {{
                        local.selectedEdgeId = local.selectedEdgeIds[local.selectedEdgeIds.length - 1];
                    }}
                    if (
                        !local.selectedPtmId ||
                        !Array.isArray(local.autoPtmPlacements) ||
                        !local.autoPtmPlacements.some((placement) => String(placement && placement.id || '') === String(local.selectedPtmId || '')) ||
                        !local.selectedPtmNodeId ||
                        !local.editableNodes.some((node) => String(node.id || '') === String(local.selectedPtmNodeId || ''))
                    ) {{
                        local.selectedPtmId = null;
                        local.selectedPtmNodeId = null;
                    }}
                    local.globalOpacity = Math.max(0.1, Math.min(1, Number(safe.globalOpacity || 1.0)));
                    local.proteinOvalMode = !!safe.proteinOvalMode;
                    local.edgeResizeMode = !!safe.edgeResizeMode;
                    if (opacityInput) opacityInput.value = local.globalOpacity.toFixed(2);
                }};
                const cloneEditableNode = (node, overrides = {{}}) => {{
                    const copy = cloneStateSnapshot(node || {{}});
                    copy.id = addNodeId();
                    copy.pendingAnnotation = '';
                    return Object.assign(copy, overrides || {{}});
                }};
                const cloneEditableEdge = (edge, overrides = {{}}) => {{
                    const copy = cloneStateSnapshot(edge || {{}});
                    copy.id = 'cst_edge_' + Date.now() + '_' + Math.floor(Math.random() * 1000000);
                    copy.userCreated = true;
                    return Object.assign(copy, overrides || {{}});
                }};
                const updateUndoState = (local) => {{
                    if (!undoButton) return;
                    const enabled = Array.isArray(local.undoStack) && local.undoStack.length > 0;
                    undoButton.classList.toggle('is-disabled', !enabled);
                    undoButton.disabled = !enabled;
                }};
                const updateRedoState = (local) => {{
                    if (!redoButton) return;
                    const enabled = Array.isArray(local.redoStack) && local.redoStack.length > 0;
                    redoButton.classList.toggle('is-disabled', !enabled);
                    redoButton.disabled = !enabled;
                }};
                const pushUndoSnapshot = (local, snapshot) => {{
                    if (local.isRestoringHistory) return;
                    const safe = cloneStateSnapshot(snapshot);
                    const history = Array.isArray(local.undoStack) ? local.undoStack : [];
                    history.push(safe);
                    if (history.length > 100) history.shift();
                    local.undoStack = history;
                    local.redoStack = [];
                    updateUndoState(local);
                    updateRedoState(local);
                }};
                const undoLastChange = (local) => {{
                    if (!Array.isArray(local.undoStack) || !local.undoStack.length) return;
                    const snapshot = local.undoStack.pop();
                    const redoHistory = Array.isArray(local.redoStack) ? local.redoStack : [];
                    redoHistory.push(buildSnapshot(local));
                    if (redoHistory.length > 100) redoHistory.shift();
                    local.redoStack = redoHistory;
                    local.isRestoringHistory = true;
                    applySnapshot(local, snapshot);
                    local.isRestoringHistory = false;
                    updateUndoState(local);
                    updateRedoState(local);
                    renderEditableOverlay(local);
                }};
                const redoLastChange = (local) => {{
                    if (!Array.isArray(local.redoStack) || !local.redoStack.length) return;
                    const snapshot = local.redoStack.pop();
                    const undoHistory = Array.isArray(local.undoStack) ? local.undoStack : [];
                    undoHistory.push(buildSnapshot(local));
                    if (undoHistory.length > 100) undoHistory.shift();
                    local.undoStack = undoHistory;
                    local.isRestoringHistory = true;
                    applySnapshot(local, snapshot);
                    local.isRestoringHistory = false;
                    updateUndoState(local);
                    updateRedoState(local);
                    renderEditableOverlay(local);
                }};
                const showFallback = (reason, force) => {{
                    try {{
                        console.error('[MapKinase CST] showFallback', {{
                            viewerKey: '{viewer_key}',
                            reason: reason || 'unknown',
                            force: !!force,
                        }});
                    }} catch (_err) {{}}
                    if (!force && window[renderStateKey] && window[renderStateKey]['{viewer_key}'] && window[renderStateKey]['{viewer_key}'].hasSuccessfulRender) {{
                        return;
                    }}
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (local && !local.showPdfBackground) {{
                        if (fallback) fallback.style.display = 'none';
                        canvas.style.display = 'none';
                        return;
                    }}
                    if (fallback) fallback.style.display = 'block';
                    canvas.style.display = 'none';
                }};
                const cleanLabel = (value) => String(value || '')
                    .replace(/[\u0080-\u009f]/g, '')
                    .replace(/\uFFFD/g, '')
                    .replace(/[\\x1d\\x1e\\x1f]/g, '-')
                    .replace(/\\s+/g, ' ')
                    .trim();
                const blockedSuggestionWords = new Set([
                    'ACTIVATION', 'SIGNALING', 'PATHWAY', 'PATHWAYS', 'CASCADE', 'RESPONSE', 'RESPONSES',
                    'NETWORK', 'MODULE', 'TRANSCRIPTION', 'EXPRESSION', 'REGULATION', 'CYTOPLASM',
                    'NUCLEUS', 'MEMBRANE', 'EXTRACELLULAR', 'INTRACELLULAR', 'GROWTH', 'FACTOR',
                ]);
                const isSuggestibleLabel = (value) => {{
                    const label = cleanLabel(value);
                    if (!label) return false;
                    if (label.length > 28) return false;
                    if (/^[0-9.]+$/.test(label)) return false;
                    const tokens = label.toUpperCase().split(/\\s+/).filter(Boolean);
                    if (tokens.some((token) => blockedSuggestionWords.has(token))) return false;
                    return true;
                }};
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
                    for (let itemIndex = 0; itemIndex < (items || []).length; itemIndex += 1) {{
                        const item = (items || [])[itemIndex];
                        const label = cleanLabel(item.str);
                        if (!label) continue;
                        const transform = Array.isArray(item.transform) ? item.transform : [1, 0, 0, 1, 0, 0];
                        const angleDeg = (Math.atan2(Number(transform[1] || 0), Number(transform[0] || 1)) * 180) / Math.PI;
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
                            previous.itemEndIndex = itemIndex;
                            continue;
                        }}
                        merged.push({{
                            label,
                            normalized: cleanLabel(label).toUpperCase(),
                            itemStartIndex: itemIndex,
                            itemEndIndex: itemIndex,
                            x,
                            y,
                            width,
                            height,
                            fontSize,
                            angleDeg,
                            fontFamily: String(item.fontName || '"Segoe UI", Arial, sans-serif'),
                            textColor: String(item.color || item.fill || '#0f172a'),
                        }});
                    }}
                    return merged;
                }};
                const getMultilineTextLines = (value) => String(value || '').replace(/\\r\\n?/g, '\\n').split('\\n');
                const suggestAnnotationForNode = (node, mergedItems) => {{
                    const candidates = Array.isArray(mergedItems) ? mergedItems : [];
                    let bestLabel = '';
                    let bestScore = Infinity;
                    const nodeAngle = Number(node.angle || 0);
                    const shapeType = String(node.shapeType || 'ellipse');
                    const isRect = shapeType === 'rect' || shapeType === 'square';
                    const isRounded = shapeType === 'rounded';
                    const isRectLike = isRect || isRounded;
                    const isBracket = shapeType === 'bracket';
                    if (isBracket) return '';
                    for (const item of candidates) {{
                        const label = cleanLabel(item.label || item.str || '');
                        if (!isSuggestibleLabel(label)) continue;
                        const itemCenter = {{
                            x: Number(item.x || 0) + (Number(item.width || 0) * 0.5),
                            y: pageHeight - (Number(item.y || 0) + (Number(item.height || 0) * 0.45)),
                        }};
                        const localPoint = rotatePoint(
                            itemCenter.x - Number(node.cx || 0),
                            itemCenter.y - Number(node.cy || 0),
                            -nodeAngle,
                        );
                        const rx = Math.max(Number(node.rx || 12), 1);
                        const ry = Math.max(Number(node.ry || 9), 1);
                        let shapeScore = 0;
                        if (isRectLike) {{
                            shapeScore = Math.max(Math.abs(localPoint.x) / rx, Math.abs(localPoint.y) / ry);
                            if (shapeScore > 1.2) continue;
                        }} else {{
                            shapeScore = ((localPoint.x / rx) ** 2) + ((localPoint.y / ry) ** 2);
                            if (shapeScore > 1.35) continue;
                        }}
                        const dist = Math.hypot(localPoint.x, localPoint.y);
                        const widthPenalty = Math.abs((Number(item.width || 0) * 0.5) - rx) * 0.05;
                        const score = (shapeScore * 9) + dist + widthPenalty;
                        if (score < bestScore) {{
                            bestScore = score;
                            bestLabel = label;
                        }}
                    }}
                    return bestLabel;
                }};
                const autoLabelCustomNodes = (local) => {{
                    const mergedItems = Array.isArray(local.mergedTextItems) ? local.mergedTextItems : [];
                    if (!mergedItems.length) return;
                    const targetNodes = (local.editableNodes || []).filter((node) =>
                        node.isCustom &&
                        !node.isDrawingShape &&
                        String(node.shapeType || 'ellipse') !== 'bracket' &&
                        !String(node.displayLabel || '').trim()
                    );
                    if (!targetNodes.length) return;
                    const snapshot = buildSnapshot(local);
                    let changed = false;
                    for (const node of targetNodes) {{
                        const suggestion = suggestAnnotationForNode(node, mergedItems);
                        if (!suggestion) continue;
                        node.annotation = suggestion;
                        node.displayLabel = suggestion;
                        node.label = suggestion;
                        node.annotationCommitted = false;
                        node.pendingAnnotation = '';
                        changed = true;
                    }}
                    if (!changed) return;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                }};
                const textItemCenterInsideProteinNode = (local, centerX, centerY) => {{
                    for (const otherNode of (local.editableNodes || [])) {{
                        if (!otherNode) continue;
                        const shapeType = String(otherNode.shapeType || 'ellipse');
                        if (otherNode.isDrawingShape) continue;
                        if (shapeType === 'text' || shapeType === 'bracket') continue;
                        const localPoint = rotatePoint(
                            Number(centerX || 0) - Number(otherNode.cx || 0),
                            Number(centerY || 0) - Number(otherNode.cy || 0),
                            -Number(otherNode.angle || 0),
                        );
                        const rx = Math.max(Number(otherNode.rx || 12), 1);
                        const ry = Math.max(Number(otherNode.ry || 9), 1);
                        if (shapeType === 'rect' || shapeType === 'square' || shapeType === 'rounded') {{
                            if (Math.abs(Number(localPoint.x || 0)) <= rx && Math.abs(Number(localPoint.y || 0)) <= ry) return true;
                        }} else {{
                            const normX = Number(localPoint.x || 0) / rx;
                            const normY = Number(localPoint.y || 0) / ry;
                            if (((normX * normX) + (normY * normY)) <= 1) return true;
                        }}
                    }}
                    return false;
                }};
                const collectTextItemsForTextNode = (node, mergedItems) => {{
                    const hits = [];
                    const nodeAngle = Number(node.angle || 0);
                    const rx = Math.max(Number(node.rx || 12), 1);
                    const ry = Math.max(Number(node.ry || 9), 1);
                    for (const item of (mergedItems || [])) {{
                        const itemWidth = Number(item.width || 0);
                        const itemHeight = Number(item.height || 0);
                        const centerX = Number(item.x || 0) + (itemWidth * 0.5);
                        const centerY = pageHeight - (Number(item.y || 0) + (itemHeight * 0.45));
                        if (textItemCenterInsideProteinNode(window[renderStateKey] && window[renderStateKey]['{viewer_key}'], centerX, centerY)) continue;
                        const localPoint = rotatePoint(centerX - Number(node.cx || 0), centerY - Number(node.cy || 0), -nodeAngle);
                        if (Math.abs(Number(localPoint.x || 0)) > (rx + 16) || Math.abs(Number(localPoint.y || 0)) > (ry + 12)) continue;
                        hits.push({{
                            item,
                            localX: Number(localPoint.x || 0),
                            localY: Number(localPoint.y || 0),
                        }});
                    }}
                    const coreHits = hits.filter((entry) =>
                        Math.abs(Number(entry.localX || 0)) <= (rx * 0.72) &&
                        Math.abs(Number(entry.localY || 0)) <= (ry * 0.72)
                    );
                    const finalHits = coreHits.length ? coreHits : hits;
                    finalHits.sort((a, b) => (a.localY - b.localY) || (a.localX - b.localX));
                    return finalHits;
                }};
                const suggestTextBoxContent = (node, mergedItems) => {{
                    const hits = collectTextItemsForTextNode(node, mergedItems);
                    if (!hits.length) return null;
                    const groups = [];
                    const lineThreshold = Math.max(6, Math.min(...hits.map((entry) => Math.max(Number(entry.item.fontSize || 0), 6))) * 0.9);
                    for (const hit of hits) {{
                        const last = groups.length ? groups[groups.length - 1] : null;
                        if (!last || Math.abs(hit.localY - last.anchorY) > lineThreshold) {{
                            groups.push({{
                                anchorY: hit.localY,
                                items: [hit],
                            }});
                        }} else {{
                            last.items.push(hit);
                            last.anchorY = ((last.anchorY * (last.items.length - 1)) + hit.localY) / last.items.length;
                        }}
                    }}
                    const lines = groups.map((group) => {{
                        const items = group.items.slice().sort((a, b) => a.localX - b.localX);
                        return items.map((entry) => cleanLabel(entry.item.label || '')).filter(Boolean).join(' ').trim();
                    }}).filter(Boolean);
                    if (!lines.length) return null;
                    const allItems = groups.flatMap((group) => group.items);
                    const lead = allItems[0].item;
                    const avgLocalX = allItems.reduce((sum, entry) => sum + Number(entry.localX || 0), 0) / Math.max(allItems.length, 1);
                    const alignThreshold = Math.max(6, Number(node.rx || 12) * 0.18);
                    let textAlign = 'center';
                    if (avgLocalX < -alignThreshold) textAlign = 'left';
                    else if (avgLocalX > alignThreshold) textAlign = 'right';
                    return {{
                        text: lines.join('\\n'),
                        fontSize: Math.max(9, Math.min(28, Number(lead.fontSize || 11))),
                        fontFamily: String(lead.fontFamily || '"Segoe UI", Arial, sans-serif'),
                        textColor: String(lead.textColor || '#0f172a'),
                        textAlign,
                        angle: Number(lead.angleDeg || 0),
                    }};
                }};
                const autoLabelTextNodes = (local) => {{
                    const mergedItems = Array.isArray(local.mergedTextItems) ? local.mergedTextItems : [];
                    if (!mergedItems.length) return;
                    const targetNodes = (local.editableNodes || []).filter((node) => {{
                        if (String(node.shapeType || 'ellipse') !== 'text') return false;
                        const currentText = String(node.displayLabel || node.label || '').trim();
                        return !currentText || currentText === 'Text';
                    }});
                    if (!targetNodes.length) return;
                    const snapshot = buildSnapshot(local);
                    let changed = false;
                    for (const node of targetNodes) {{
                        const suggestion = suggestTextBoxContent(node, mergedItems);
                        if (!suggestion) continue;
                        node.displayLabel = suggestion.text;
                        node.label = suggestion.text;
                        node.originalLabel = suggestion.text;
                        node.fontSize = suggestion.fontSize;
                        node.fontFamily = suggestion.fontFamily;
                        node.textColor = suggestion.textColor;
                        node.textAlign = suggestion.textAlign;
                        node.angle = suggestion.angle;
                        const metrics = estimateTextNodeMetrics(suggestion.text, suggestion.fontSize);
                        node.rx = metrics.rx;
                        node.ry = metrics.ry;
                        changed = true;
                    }}
                    if (!changed) return;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                }};
                const findBestTextItemForNodeLabel = (node, mergedItems) => {{
                    const label = cleanLabel(node.displayLabel || node.annotation || '');
                    if (!label) return null;
                    const targetVariants = new Set(labelVariants(label));
                    if (!targetVariants.size) return null;
                    let bestItem = null;
                    let bestScore = Infinity;
                    for (const item of (mergedItems || [])) {{
                        const itemVariants = labelVariants(item.label || item.normalized || '');
                        if (!itemVariants.some((variant) => targetVariants.has(variant))) continue;
                        const itemCenterX = Number(item.x || 0) + (Number(item.width || 0) * 0.5);
                        const itemCenterY = pageHeight - (Number(item.y || 0) + (Number(item.height || 0) * 0.45));
                        const dist = Math.hypot(itemCenterX - Number(node.cx || 0), itemCenterY - Number(node.cy || 0));
                        const score = dist + Math.abs((Number(item.width || 0) * 0.5) - Number(node.rx || 12)) * 0.1;
                        if (score < bestScore) {{
                            bestScore = score;
                            bestItem = item;
                        }}
                    }}
                    return bestItem;
                }};
                const findBestEllipseForTextItem = (item) => {{
                    if (!item) return null;
                    const itemWidth = Number(item.width || 0);
                    const itemHeight = Number(item.height || 0);
                    const itemX = Number(item.x || 0);
                    const itemY = Number(item.y || 0);
                    const itemCenterX = itemX + (itemWidth * 0.5);
                    const itemCenterY = pageHeight - (itemY + (itemHeight * 0.45));
                    const expectedWidth = Math.max(itemWidth * 0.55, 8);
                    const ownerWidth = Math.min(
                        itemWidth,
                        Math.max(12, Math.max(Number(item.fontSize || 0), 6) * Math.min(Math.max(cleanLabel(item.label || '').length, 1), 3.5) * 0.55)
                    );
                    const left = itemX - 8;
                    const right = itemX + ownerWidth;
                    const top = itemCenterY - (itemHeight * 0.72);
                    const bottom = itemCenterY + (itemHeight * 0.42);
                    const anchorX = itemX;
                    const anchorY = itemCenterY;
                    let bestEllipse = null;
                    let bestScore = Infinity;
                    for (const ellipse of (ellipseGroups || [])) {{
                        const centerX = Number(ellipse.center_x || 0);
                        const centerY = pageHeight - Number(ellipse.center_y || 0);
                        const rx = Number(ellipse.radius_x || 0);
                        const ry = Number(ellipse.radius_y || 0);
                        if (rx < 7 || ry < 6) continue;
                        const dist = Math.hypot(centerX - itemCenterX, centerY - itemCenterY);
                        if (dist > 48) continue;
                        const anchorDx = Math.abs(centerX - anchorX);
                        const anchorDy = Math.abs(centerY - anchorY);
                        const closestX = Math.max(left, Math.min(centerX, right));
                        const closestY = Math.max(top, Math.min(centerY, bottom));
                        const boxDx = centerX - closestX;
                        const boxDy = centerY - closestY;
                        const effectiveDist = Math.hypot(boxDx, boxDy);
                        const outsideX = centerX < left ? (left - centerX) : Math.max(0, centerX - right);
                        const outsideY = centerY < top ? (top - centerY) : Math.max(0, centerY - bottom);
                        const overlapPenalty = ((outsideX / Math.max(rx, 1)) ** 2) + ((outsideY / Math.max(ry, 1)) ** 2);
                        if (
                            centerX < (left - (rx * 0.6)) ||
                            centerX > (right + (rx * 0.6)) ||
                            centerY < (top - (ry * 0.7)) ||
                            centerY > (bottom + (ry * 0.7))
                        ) {{
                            continue;
                        }}
                        const sizePenalty = Math.abs(rx - expectedWidth) * 0.18;
                        const score = (effectiveDist * 1.6) + (anchorDx * 0.55) + (anchorDy * 0.12) + sizePenalty + (overlapPenalty * 8) + (dist * 0.04);
                        if (score < bestScore) {{
                            bestScore = score;
                            bestEllipse = ellipse;
                        }}
                    }}
                    return bestEllipse;
                }};
                const findBestRectForTextItem = (item, node) => {{
                    if (!item) return null;
                    const itemWidth = Number(item.width || 0);
                    const itemHeight = Number(item.height || 0);
                    const itemX = Number(item.x || 0);
                    const itemY = Number(item.y || 0);
                    const itemCenterX = itemX + (itemWidth * 0.5);
                    const itemCenterY = pageHeight - (itemY + (itemHeight * 0.45));
                    const ownerWidth = Math.min(
                        itemWidth,
                        Math.max(12, Math.max(Number(item.fontSize || 0), 6) * Math.min(Math.max(cleanLabel(item.label || '').length, 1), 3.5) * 0.55)
                    );
                    const left = itemX - 8;
                    const right = itemX + ownerWidth;
                    const top = itemCenterY - (itemHeight * 0.72);
                    const bottom = itemCenterY + (itemHeight * 0.42);
                    const anchorX = itemX;
                    const anchorY = itemCenterY;
                    let bestRect = null;
                    let bestScore = Infinity;
                    for (const rect of (rectGroups || [])) {{
                        const centerX = Number(rect.center_x || 0);
                        const centerY = pageHeight - Number(rect.center_y || 0);
                        const rx = Number(rect.radius_x || 0);
                        const ry = Number(rect.radius_y || 0);
                        if (rx < 7 || ry < 5) continue;
                        const dist = Math.hypot(centerX - itemCenterX, centerY - itemCenterY);
                        if (dist > 56) continue;
                        const rectLeft = Number(rect.left || (centerX - rx));
                        const rectRight = Number(rect.right || (centerX + rx));
                        const rectTop = pageHeight - Number(rect.top || (Number(rect.center_y || 0) + ry));
                        const rectBottom = pageHeight - Number(rect.bottom || (Number(rect.center_y || 0) - ry));
                        const closestX = Math.max(left, Math.min(centerX, right));
                        const closestY = Math.max(top, Math.min(centerY, bottom));
                        const boxDx = centerX - closestX;
                        const boxDy = centerY - closestY;
                        const effectiveDist = Math.hypot(boxDx, boxDy);
                        const overlapX = Math.max(0, Math.min(right, rectRight) - Math.max(left, rectLeft));
                        const overlapY = Math.max(0, Math.min(bottom, rectBottom) - Math.max(top, rectTop));
                        const overlapBonus = (overlapX > 0 && overlapY > 0) ? ((overlapX * overlapY) ** 0.5) : 0;
                        const anchorDx = Math.abs(centerX - anchorX);
                        const anchorDy = Math.abs(centerY - anchorY);
                        const expectedWidth = Math.max(itemWidth * 0.6, 10);
                        const sizePenalty = Math.abs((rx * 2) - expectedWidth) * 0.1;
                        const score = (effectiveDist * 1.5) + (anchorDx * 0.45) + (anchorDy * 0.12) + sizePenalty - (overlapBonus * 0.18) + (dist * 0.03);
                        if (score < bestScore) {{
                            bestScore = score;
                            bestRect = rect;
                        }}
                    }}
                    return bestRect;
                }};
                const autoSizeCustomNodes = (local) => {{
                    const mergedItems = Array.isArray(local.mergedTextItems) ? local.mergedTextItems : [];
                    if (!mergedItems.length) return;
                    const targets = (local.editableNodes || []).filter((node) =>
                        node.isCustom &&
                        String(node.shapeType || 'ellipse') !== 'bracket' &&
                        cleanLabel(node.displayLabel || node.annotation || '')
                    );
                    if (!targets.length) return;
                    const snapshot = buildSnapshot(local);
                    let changed = false;
                    for (const node of targets) {{
                        const item = findBestTextItemForNodeLabel(node, mergedItems);
                        if (String(node.shapeType || 'ellipse') === 'rect') {{
                            const rectMatch = findBestRectForTextItem(item, node);
                            if (!rectMatch) continue;
                            node.cx = Number(rectMatch.center_x || node.cx);
                            node.cy = pageHeight - Number(rectMatch.center_y || (pageHeight - Number(node.cy || 0)));
                            node.rx = Math.max(4, Number(rectMatch.radius_x || node.rx));
                            node.ry = Math.max(4, Number(rectMatch.radius_y || node.ry));
                            node.angle = 0;
                            changed = true;
                            continue;
                        }}
                        const ellipse = findBestEllipseForTextItem(item);
                        if (!ellipse) continue;
                        node.cx = Number(ellipse.center_x || node.cx);
                        node.cy = pageHeight - Number(ellipse.center_y || (pageHeight - Number(node.cy || 0)));
                        node.rx = Math.max(4, Number(ellipse.radius_x || node.rx));
                        node.ry = Math.max(4, Number(ellipse.radius_y || node.ry));
                        node.angle = 0;
                        changed = true;
                    }}
                    if (!changed) return;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                }};
                const renderOverlayNodes = (textItems) => {{
                    if (!overlaySvg) return [];
                    const mergedItems = buildMergedTextItems(textItems);
                    const available = new Map();
                    const availableEllipses = (ellipseGroups || []).map((ellipse, index) => ({{ ...ellipse, _used: false, _idx: index }}));
                    for (const item of mergedItems) {{
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
                    const editableNodes = [];
                    for (let nodeIndex = 0; nodeIndex < overlayData.length; nodeIndex += 1) {{
                        const node = overlayData[nodeIndex];
                        const match = assignments.get(nodeIndex);
                        const item = match ? match.item : null;
                        const color = node['fc_color_{active_idx}'] || node.fc_color_1 || node.default_color || [166, 166, 166];
                        const rgb = Array.isArray(color) ? color.slice(0, 3) : [166, 166, 166];
                        const outlineColor = node['outline_color_{active_idx}'] || node.outline_color_1 || [0, 0, 0];
                        const outlineRgb = Array.isArray(outlineColor) ? outlineColor.slice(0, 3) : [0, 0, 0];
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
                        const strokeWidth = Math.max(0.6, Number(node.stroke_width || {prot_outline_width}));
                        const cx = useFallbackEllipse
                            ? Number(fallbackGeom.cx || 0)
                            : Number(ellipseMatch.center_x || 0);
                        const cyPdf = useFallbackEllipse
                            ? Number(fallbackGeom.cyPdf || 0)
                            : Number(ellipseMatch.center_y || 0);
                        const label = String(node.label || '');
                        const matchGene = String(node['matched_gene_symbol_{active_idx}'] || node.matched_gene_symbol_1 || '');
                        const matchUniprot = String(node['matched_uniprot_{active_idx}'] || node.matched_uniprot_1 || '');
                        const candidateUniprots = [];
                        const seenCandidateUniprots = new Set();
                        const pushCandidateUniprot = (value) => {{
                            const raw = String(value || '').trim();
                            if (!raw) return;
                            const normalized = raw.split(/[|,;/]/)[0].trim();
                            if (!normalized) return;
                            const key = normalized.toUpperCase();
                            if (seenCandidateUniprots.has(key)) return;
                            seenCandidateUniprots.add(key);
                            candidateUniprots.push(normalized);
                        }};
                        pushCandidateUniprot(matchUniprot);
                        (Array.isArray(node.candidate_uniprot_ids) ? node.candidate_uniprot_ids : []).forEach(pushCandidateUniprot);
                        (Array.isArray(node.suggested_uniprot_ids) ? node.suggested_uniprot_ids : []).forEach(pushCandidateUniprot);
                        const foldValue = node['fold_change_{active_idx}'];
                        const foldText = Number.isFinite(Number(foldValue)) ? Number(foldValue).toFixed(3) : '';
                        const isComplex = !!node.is_complex;
                        const tooltipText = String(node.tooltip || '');
                        editableNodes.push({{
                            id: addNodeId(),
                            originalLabel: label,
                            displayLabel: matchGene || label,
                            label,
                            matchedGene: matchGene,
                            matchedUniprot: matchUniprot,
                            candidateUniprotIds: candidateUniprots,
                            suggestedGeneSymbols: Array.isArray(node.suggested_gene_symbols) ? node.suggested_gene_symbols.slice() : [],
                            proteinOptions: Array.isArray(node.protein_options) ? node.protein_options.slice() : [],
                            foldText,
                            hasDatasetMatch: !!node.has_dataset_match,
                            cx,
                            cy: pageHeight - cyPdf,
                            rx,
                            ry,
                            shapeType: 'ellipse',
                            angle: 0,
                            strokeWidth,
                            stroke: 'rgb(' + outlineRgb[0] + ', ' + outlineRgb[1] + ', ' + outlineRgb[2] + ')',
                            fillColor: 'rgb(' + rgb[0] + ', ' + rgb[1] + ', ' + rgb[2] + ')',
                            ['outline_color_{active_idx}']: outlineRgb,
                            outline_color_1: outlineRgb,
                            opacity: 1.0,
                            annotation: '',
                            annotationCommitted: true,
                            pendingAnnotation: '',
                            isCustom: false,
                            className: isComplex ? 'cst-overlay-ellipse cst-complex-node' : (node.has_dataset_match ? 'cst-overlay-ellipse' : 'cst-overlay-ellipse cst-missing-node'),
                            title: tooltipText || (isComplex ? 'Mapping: complex/process module' : ''),
                            tooltip: tooltipText,
                            tooltipHtml: String(node.tooltip_html || ''),
                            mappingType: String(node.mapping_type || ''),
                            isComplex,
                            complexDisplayMembers: Array.isArray(node.complex_display_members) ? node.complex_display_members : [],
                        }});
                        for (const key of Object.keys(node || {{}})) {{
                            if (!/^(?:fc_color|fold_change|outline_color|outline_fold_change)_\\d+$/.test(key)) continue;
                            editableNodes[editableNodes.length - 1][key] = node[key];
                        }}
                    }}
                    return normalizeEditableNodeIds(editableNodes);
                }};
                const applyActiveFcToEditableNodes = (nodes) => {{
                    if (!Array.isArray(nodes)) return nodes;
                    for (const node of nodes) {{
                        if (!node || typeof node !== 'object') continue;
                        const shapeType = String(node.shapeType || '').toLowerCase();
                        if (shapeType === 'text' || node.isDrawingShape) continue;
                        const fcColor = node['fc_color_{active_idx}'] || node.fc_color_1 || null;
                        if (Array.isArray(fcColor) && fcColor.length >= 3) {{
                            node.fillColor = 'rgb(' + fcColor[0] + ', ' + fcColor[1] + ', ' + fcColor[2] + ')';
                        }}
                        const outlineColor = node['outline_color_{active_idx}'] || node.outline_color_1 || null;
                        if (Array.isArray(outlineColor) && outlineColor.length >= 3) {{
                            node.stroke = 'rgb(' + outlineColor[0] + ', ' + outlineColor[1] + ', ' + outlineColor[2] + ')';
                        }}
                        const foldValue = node['fold_change_{active_idx}'];
                        node.foldText = Number.isFinite(Number(foldValue)) ? Number(foldValue).toFixed(3) : '';
                        const activeGene = String(node['matched_gene_symbol_{active_idx}'] || '').trim();
                        const activeUniprot = String(node['matched_uniprot_{active_idx}'] || '').trim();
                        if (activeGene) {{
                            node.matchedGene = activeGene;
                            node.displayLabel = activeGene;
                        }}
                        if (activeUniprot) {{
                            node.matchedUniprot = activeUniprot;
                        }}
                    }}
                    return nodes;
                }};
                const setDeleteState = (local) => {{
                    const selectedCount = getSelectedNodeIds(local).length;
                    const selectedEdgeCount = getSelectedEdgeIds(local).length;
                    const enabled = selectedCount > 0 || selectedEdgeCount > 0 || !!local.selectedPtmId;
                    const multiEnabled = selectedCount > 1;
                    if (deleteButton) {{
                        deleteButton.classList.toggle('is-disabled', !enabled);
                        deleteButton.disabled = !enabled;
                    }}
                    if (convertRectButton) {{
                        convertRectButton.classList.toggle('is-disabled', !enabled);
                        convertRectButton.disabled = !enabled;
                    }}
                    if (batchSizeButton) {{
                        batchSizeButton.classList.toggle('is-disabled', !multiEnabled);
                        batchSizeButton.disabled = !multiEnabled;
                        if (!multiEnabled) batchSizeButton.classList.remove('is-active');
                    }}
                    if (!multiEnabled) local.batchSizePanelOpen = false;
                    if (batchSizePanel) {{
                        batchSizePanel.classList.toggle('is-open', !!local.batchSizePanelOpen && multiEnabled);
                    }}
                    if (bringFrontButton) {{
                        bringFrontButton.classList.toggle('is-disabled', !enabled);
                        bringFrontButton.disabled = !enabled;
                    }}
                    if (sendBackButton) {{
                        sendBackButton.classList.toggle('is-disabled', !enabled);
                        sendBackButton.disabled = !enabled;
                    }}
                    const edgeEnabled = selectedEdgeCount > 0;
                    if (dashedArrowButton) {{
                        dashedArrowButton.classList.toggle('is-disabled', !edgeEnabled);
                        dashedArrowButton.disabled = !edgeEnabled;
                    }}
                    if (bothArrowButton) {{
                        bothArrowButton.classList.toggle('is-disabled', !edgeEnabled);
                        bothArrowButton.disabled = !edgeEnabled;
                    }}
                    if (lineArrowButton) {{
                        lineArrowButton.classList.toggle('is-disabled', !edgeEnabled);
                        lineArrowButton.disabled = !edgeEnabled;
                    }}
                    if (inhibitorArrowButton) {{
                        inhibitorArrowButton.classList.toggle('is-disabled', !edgeEnabled);
                        inhibitorArrowButton.disabled = !edgeEnabled;
                    }}
                }};
                const setAddState = (local, isActive) => {{
                    local.addMode = !!isActive;
                    stage.classList.toggle('cst-add-mode', !!isActive);
                    if (addButton) addButton.classList.toggle('is-active', !!isActive && local.addShapeType === 'ellipse');
                    if (addRectButton) addRectButton.classList.toggle('is-active', !!isActive && local.addShapeType === 'rect');
                    if (addBracketButton) addBracketButton.classList.toggle('is-active', !!isActive && local.addShapeType === 'bracket');
                    if (addTextButton) addTextButton.classList.toggle('is-active', !!isActive && local.addShapeType === 'text');
                }};
                const setAddArrowState = (local, isActive) => {{
                    local.addArrowMode = !!isActive;
                    if (!isActive) {{
                        local.pendingArrowStart = null;
                        local.pendingArrowPreview = null;
                    }}
                    stage.classList.toggle('cst-add-arrow-mode', !!isActive);
                    if (addArrowButton) addArrowButton.classList.toggle('is-active', !!isActive);
                }};
                const setProteinOvalState = (local, isActive) => {{
                    local.proteinOvalMode = !!isActive;
                    if (proteinOvalButton) proteinOvalButton.classList.toggle('is-active', !!isActive);
                }};
                const setEdgeResizeState = (local, isActive) => {{
                    local.edgeResizeMode = !!isActive;
                    if (edgeResizeButton) edgeResizeButton.classList.toggle('is-active', !!isActive);
                }};
                const estimateTextNodeMetrics = (text, fontSize) => {{
                    const lines = getMultilineTextLines(text).map((line) => String(line || '').trimEnd());
                    const nonEmptyLines = lines.filter((line) => line.trim().length);
                    const label = (nonEmptyLines.length ? nonEmptyLines : ['Text']);
                    const size = Math.max(8, Number(fontSize || 11));
                    const longestLine = label.reduce((best, line) => Math.max(best, String(line || '').length), 1);
                    const lineCount = Math.max(label.length, 1);
                    const fullWidth = Math.max(26, (longestLine * size * 0.58) + 14);
                    return {{
                        rx: Math.max(13, fullWidth * 0.5),
                        ry: Math.max(7, ((lineCount * size * 0.62) * 0.5) + 5),
                    }};
                }};
                const updateInlineTextEditorLayout = (local) => {{
                    if (!inlineTextEditor || !local || !local.editingTextNodeId) return;
                    const node = findNodeById(local, local.editingTextNodeId);
                    if (!node) return;
                    const scaleX = Math.max(Number(content.clientWidth || 0), 1) / Math.max(pageWidth, 1);
                    const scaleY = Math.max(Number(content.clientHeight || 0), 1) / Math.max(pageHeight, 1);
                    inlineTextEditor.style.left = ((Number(node.cx || 0) - Number(node.rx || 12)) * scaleX) + 'px';
                    inlineTextEditor.style.top = ((Number(node.cy || 0) - Number(node.ry || 9)) * scaleY) + 'px';
                    inlineTextEditor.style.width = Math.max(12, Number(node.rx || 12) * 2 * scaleX) + 'px';
                    inlineTextEditor.style.height = Math.max(12, Number(node.ry || 9) * 2 * scaleY) + 'px';
                    inlineTextEditor.style.fontSize = Math.max(9, Number(node.fontSize || 11) * scaleY) + 'px';
                    inlineTextEditor.style.fontWeight = String(node.fontWeight || '600');
                    inlineTextEditor.style.fontFamily = String(node.fontFamily || '"Segoe UI", Arial, sans-serif');
                    inlineTextEditor.style.color = String(node.textColor || '#0f172a');
                    inlineTextEditor.style.borderColor = String(node.stroke || '#64748b');
                    inlineTextEditor.style.textAlign = String(node.textAlign || 'center');
                    inlineTextEditor.style.lineHeight = Math.max(12, Number(node.fontSize || 11) * scaleY * 1.14) + 'px';
                    inlineTextEditor.style.transform = Number(node.angle || 0)
                        ? ('rotate(' + Number(node.angle || 0).toFixed(3) + 'deg)')
                        : 'none';
                    inlineTextEditor.style.transformOrigin = 'center center';
                }};
                const finishInlineTextEdit = (local, commit) => {{
                    if (!inlineTextEditor || !local || !local.editingTextNodeId) return false;
                    const node = findNodeById(local, local.editingTextNodeId);
                    const beforeSnapshot = local.editingTextSnapshot || buildSnapshot(local);
                    if (node && commit) {{
                        const nextText = String(inlineTextEditor.value || '').replace(/\\r\\n?/g, '\\n').trim() || 'Text';
                        node.displayLabel = nextText;
                        node.label = nextText;
                        node.originalLabel = nextText;
                        const metrics = estimateTextNodeMetrics(nextText, node.fontSize || 11);
                        node.rx = metrics.rx;
                        node.ry = metrics.ry;
                    }}
                    inlineTextEditor.blur();
                    inlineTextEditor.classList.remove('is-open');
                    inlineTextEditor.value = '';
                    local.editingTextNodeId = null;
                    local.editingTextSnapshot = null;
                    if (node && commit) {{
                        const afterSnapshot = JSON.stringify(buildSnapshot(local));
                        if (JSON.stringify(beforeSnapshot) !== afterSnapshot) {{
                            pushUndoSnapshot(local, beforeSnapshot);
                        }}
                    }}
                    renderEditableOverlay(local);
                    return true;
                }};
                const startInlineTextEdit = (local, nodeId) => {{
                    if (!inlineTextEditor || !local) return;
                    const node = findNodeById(local, nodeId);
                    if (!node || String(node.shapeType || 'ellipse') !== 'text') return;
                    local.editingTextNodeId = String(nodeId || '');
                    local.editingTextSnapshot = buildSnapshot(local);
                    inlineTextEditor.value = String(node.displayLabel || node.label || 'Text');
                    inlineTextEditor.classList.add('is-open');
                    updateInlineTextEditorLayout(local);
                    window.setTimeout(() => {{
                        inlineTextEditor.focus();
                        inlineTextEditor.select();
                    }}, 0);
                    renderEditableOverlay(local);
                }};
                const findNodeById = (local, nodeId) => Array.isArray(local.editableNodes)
                    ? local.editableNodes.find((node) => node.id === nodeId)
                    : null;
                const getSelectedNodeIds = (local) => {{
                    if (!Array.isArray(local.selectedNodeIds)) local.selectedNodeIds = [];
                    const valid = Array.from(new Set(local.selectedNodeIds.filter((nodeId) => findNodeById(local, nodeId))));
                    if (local.selectedNodeId && !valid.includes(local.selectedNodeId) && findNodeById(local, local.selectedNodeId)) {{
                        valid.push(local.selectedNodeId);
                    }}
                    local.selectedNodeIds = valid;
                    return valid;
                }};
                const normalizeRectAngle = (angle) => {{
                    const normalized = ((Number(angle || 0) % 180) + 180) % 180;
                    return normalized > 90 ? normalized - 180 : normalized;
                }};
                const snapRotationAngle = (angle) => {{
                    const raw = Number(angle || 0);
                    const snapStep = 90;
                    const snapThreshold = 6;
                    const nearest = Math.round(raw / snapStep) * snapStep;
                    return Math.abs(raw - nearest) <= snapThreshold ? nearest : raw;
                }};
                const isAxisAlignedRect = (node) => {{
                    if (String(node && node.shapeType || 'ellipse') !== 'rect') return false;
                    return Math.abs(normalizeRectAngle(node && node.angle || 0)) <= 8;
                }};
                const getRectBounds = (node, cxOverride = null, cyOverride = null) => {{
                    const cx = Number(cxOverride !== null ? cxOverride : node && node.cx || 0);
                    const cy = Number(cyOverride !== null ? cyOverride : node && node.cy || 0);
                    const rx = Number(node && node.rx || 12);
                    const ry = Number(node && node.ry || 9);
                    return {{
                        left: cx - rx,
                        right: cx + rx,
                        top: cy - ry,
                        bottom: cy + ry,
                        width: rx * 2,
                        height: ry * 2,
                    }};
                }};
                const applyRectMoveSnap = (local, node, nextCx, nextCy) => {{
                    if (!isAxisAlignedRect(node)) return {{ cx: nextCx, cy: nextCy }};
                    const snapGap = 10;
                    const overlapMin = 8;
                    const moving = getRectBounds(node, nextCx, nextCy);
                    let bestX = null;
                    let bestY = null;
                    let bestXDist = Infinity;
                    let bestYDist = Infinity;
                    const selected = new Set(getSelectedNodeIds(local));
                    for (const other of (local.editableNodes || [])) {{
                        if (!other || other.id === node.id || selected.has(other.id)) continue;
                        if (!isAxisAlignedRect(other)) continue;
                        const bounds = getRectBounds(other);
                        const verticalOverlap = Math.min(moving.bottom, bounds.bottom) - Math.max(moving.top, bounds.top);
                        const horizontalOverlap = Math.min(moving.right, bounds.right) - Math.max(moving.left, bounds.left);
                        const horizontalGapLeft = Math.abs(moving.left - bounds.right);
                        const horizontalGapRight = Math.abs(moving.right - bounds.left);
                        const verticalGapTop = Math.abs(moving.top - bounds.bottom);
                        const verticalGapBottom = Math.abs(moving.bottom - bounds.top);
                        if (verticalOverlap >= overlapMin) {{
                            if (horizontalGapLeft <= snapGap && horizontalGapLeft < bestXDist) {{
                                bestXDist = horizontalGapLeft;
                                bestX = bounds.right + Number(node.rx || 12);
                            }}
                            if (horizontalGapRight <= snapGap && horizontalGapRight < bestXDist) {{
                                bestXDist = horizontalGapRight;
                                bestX = bounds.left - Number(node.rx || 12);
                            }}
                            const topFlush = Math.abs(moving.top - bounds.top);
                            const bottomFlush = Math.abs(moving.bottom - bounds.bottom);
                            if (topFlush <= snapGap && topFlush < bestYDist) {{
                                bestYDist = topFlush;
                                bestY = bounds.top + Number(node.ry || 9);
                            }}
                            if (bottomFlush <= snapGap && bottomFlush < bestYDist) {{
                                bestYDist = bottomFlush;
                                bestY = bounds.bottom - Number(node.ry || 9);
                            }}
                        }}
                        if (horizontalOverlap >= overlapMin) {{
                            if (verticalGapTop <= snapGap && verticalGapTop < bestYDist) {{
                                bestYDist = verticalGapTop;
                                bestY = bounds.bottom + Number(node.ry || 9);
                            }}
                            if (verticalGapBottom <= snapGap && verticalGapBottom < bestYDist) {{
                                bestYDist = verticalGapBottom;
                                bestY = bounds.top - Number(node.ry || 9);
                            }}
                            const leftFlush = Math.abs(moving.left - bounds.left);
                            const rightFlush = Math.abs(moving.right - bounds.right);
                            if (leftFlush <= snapGap && leftFlush < bestXDist) {{
                                bestXDist = leftFlush;
                                bestX = bounds.left + Number(node.rx || 12);
                            }}
                            if (rightFlush <= snapGap && rightFlush < bestXDist) {{
                                bestXDist = rightFlush;
                                bestX = bounds.right - Number(node.rx || 12);
                            }}
                        }}
                    }}
                    return {{
                        cx: bestX !== null ? bestX : nextCx,
                        cy: bestY !== null ? bestY : nextCy,
                    }};
                }};
                const isNodeSelected = (local, nodeId) => getSelectedNodeIds(local).includes(nodeId);
                const getSelectedEdgeIds = (local) => {{
                    if (!Array.isArray(local.selectedEdgeIds)) local.selectedEdgeIds = [];
                    const valid = Array.from(new Set(local.selectedEdgeIds.filter((edgeId) => findEdgeById(local, edgeId))));
                    if (local.selectedEdgeId && !valid.includes(local.selectedEdgeId) && findEdgeById(local, local.selectedEdgeId)) {{
                        valid.push(local.selectedEdgeId);
                    }}
                    local.selectedEdgeIds = valid;
                    return valid;
                }};
                const clearSelectedPtm = (local) => {{
                    local.selectedPtmId = null;
                    local.selectedPtmNodeId = null;
                }};
                const setSingleSelection = (local, nodeId) => {{
                    local.selectedNodeId = nodeId || null;
                    local.selectedNodeIds = nodeId ? [nodeId] : [];
                    local.selectedEdgeId = null;
                    local.selectedEdgeIds = [];
                    clearSelectedPtm(local);
                }};
                const setSinglePtmSelection = (local, nodeId, ptmId) => {{
                    local.selectedNodeId = nodeId || null;
                    local.selectedNodeIds = nodeId ? [nodeId] : [];
                    local.selectedEdgeId = null;
                    local.selectedEdgeIds = [];
                    local.selectedPtmId = ptmId || null;
                    local.selectedPtmNodeId = nodeId || null;
                }};
                const toggleNodeSelection = (local, nodeId) => {{
                    clearSelectedPtm(local);
                    const selected = getSelectedNodeIds(local).slice();
                    local.selectedEdgeId = null;
                    local.selectedEdgeIds = [];
                    const idx = selected.indexOf(nodeId);
                    if (idx >= 0) {{
                        selected.splice(idx, 1);
                        local.selectedNodeIds = selected;
                        local.selectedNodeId = selected.length ? selected[selected.length - 1] : null;
                        return;
                    }}
                    selected.push(nodeId);
                    local.selectedNodeIds = selected;
                    local.selectedNodeId = nodeId;
                }};
                const entryKey = (type, id) => `${{String(type || '').trim().toLowerCase()}}:${{String(id || '').trim()}}`;
                const normalizeGroupEntry = (entry) => {{
                    if (!entry || typeof entry !== 'object') return null;
                    const type = String(entry.type || '').trim().toLowerCase();
                    const id = String(entry.id || '').trim();
                    if ((type !== 'node' && type !== 'edge') || !id) return null;
                    return {{ type, id }};
                }};
                const sanitizeGroups = (local, groups) => {{
                    const cleaned = [];
                    const seenGroups = new Set();
                    for (const raw of (Array.isArray(groups) ? groups : [])) {{
                        if (!raw || typeof raw !== 'object') continue;
                        const groupId = String(raw.id || raw.group_id || '').trim();
                        if (!groupId || seenGroups.has(groupId)) continue;
                        seenGroups.add(groupId);
                        const members = [];
                        const seenMembers = new Set();
                        for (const member of (Array.isArray(raw.members) ? raw.members : [])) {{
                            const normalized = normalizeGroupEntry(member);
                            if (!normalized) continue;
                            const exists = normalized.type === 'node'
                                ? !!findNodeById(local, normalized.id)
                                : !!findEdgeById(local, normalized.id);
                            if (!exists) continue;
                            const key = entryKey(normalized.type, normalized.id);
                            if (seenMembers.has(key)) continue;
                            seenMembers.add(key);
                            members.push(normalized);
                        }}
                        if (members.length >= 2) {{
                            cleaned.push({{ id: groupId, members }});
                        }}
                    }}
                    return cleaned;
                }};
                const getSelectedEntries = (local) => {{
                    const entries = [];
                    getSelectedNodeIds(local).forEach((nodeId) => entries.push({{ type: 'node', id: nodeId }}));
                    getSelectedEdgeIds(local).forEach((edgeId) => entries.push({{ type: 'edge', id: edgeId }}));
                    return entries;
                }};
                const isEntrySelected = (local, type, id) => {{
                    if (String(type || '').toLowerCase() === 'edge') {{
                        return getSelectedEdgeIds(local).includes(String(id || ''));
                    }}
                    return getSelectedNodeIds(local).includes(String(id || ''));
                }};
                const getGroupIdsForEntry = (local, type, id) => {{
                    const key = entryKey(type, id);
                    return (Array.isArray(local.groups) ? local.groups : [])
                        .filter((group) => Array.isArray(group.members) && group.members.some((member) => entryKey(member.type, member.id) === key))
                        .map((group) => String(group.id || ''))
                        .filter(Boolean);
                }};
                const expandEntriesWithGroups = (local, entries) => {{
                    const expanded = [];
                    const seenEntries = new Set();
                    const queue = (Array.isArray(entries) ? entries : []).map(normalizeGroupEntry).filter(Boolean);
                    for (const entry of queue) {{
                        const key = entryKey(entry.type, entry.id);
                        if (seenEntries.has(key)) continue;
                        seenEntries.add(key);
                        expanded.push(entry);
                        const groupIds = getGroupIdsForEntry(local, entry.type, entry.id);
                        for (const groupId of groupIds) {{
                            const group = (Array.isArray(local.groups) ? local.groups : []).find((item) => String(item.id || '') === String(groupId || ''));
                            if (!group || !Array.isArray(group.members)) continue;
                            group.members.forEach((member) => {{
                                const normalized = normalizeGroupEntry(member);
                                if (!normalized) return;
                                const memberKey = entryKey(normalized.type, normalized.id);
                                if (!seenEntries.has(memberKey)) {{
                                    seenEntries.add(memberKey);
                                    expanded.push(normalized);
                                }}
                            }});
                        }}
                    }}
                    return expanded;
                }};
                const setSelectedEntries = (local, entries, primaryEntry = null) => {{
                    clearSelectedPtm(local);
                    const normalized = [];
                    const seen = new Set();
                    for (const entry of (Array.isArray(entries) ? entries : [])) {{
                        const item = normalizeGroupEntry(entry);
                        if (!item) continue;
                        const exists = item.type === 'node' ? !!findNodeById(local, item.id) : !!findEdgeById(local, item.id);
                        if (!exists) continue;
                        const key = entryKey(item.type, item.id);
                        if (seen.has(key)) continue;
                        seen.add(key);
                        normalized.push(item);
                    }}
                    const nodeIds = normalized.filter((entry) => entry.type === 'node').map((entry) => entry.id);
                    const edgeIds = normalized.filter((entry) => entry.type === 'edge').map((entry) => entry.id);
                    local.selectedNodeIds = nodeIds;
                    local.selectedEdgeIds = edgeIds;
                    const primary = normalizeGroupEntry(primaryEntry);
                    if (primary && primary.type === 'node' && nodeIds.includes(primary.id)) {{
                        local.selectedNodeId = primary.id;
                    }} else if (nodeIds.length) {{
                        local.selectedNodeId = nodeIds[nodeIds.length - 1];
                    }} else {{
                        local.selectedNodeId = null;
                    }}
                    if (primary && primary.type === 'edge' && edgeIds.includes(primary.id)) {{
                        local.selectedEdgeId = primary.id;
                    }} else if (edgeIds.length) {{
                        local.selectedEdgeId = edgeIds[edgeIds.length - 1];
                    }} else {{
                        local.selectedEdgeId = null;
                    }}
                }};
                const cleanupGroups = (local) => {{
                    local.groups = sanitizeGroups(local, Array.isArray(local.groups) ? local.groups : []);
                }};
                const createGroupFromSelection = (local) => {{
                    const selectedEntries = getSelectedEntries(local);
                    if (selectedEntries.length < 2) return false;
                    const snapshot = buildSnapshot(local);
                    const selectedKeys = new Set(selectedEntries.map((entry) => entryKey(entry.type, entry.id)));
                    const overlappingGroups = (Array.isArray(local.groups) ? local.groups : []).filter((group) =>
                        Array.isArray(group.members) && group.members.some((member) => selectedKeys.has(entryKey(member.type, member.id)))
                    );
                    const combinedEntries = expandEntriesWithGroups(local, selectedEntries.concat(
                        overlappingGroups.flatMap((group) => Array.isArray(group.members) ? group.members : [])
                    ));
                    if (combinedEntries.length < 2) return false;
                    const overlappingIds = new Set(overlappingGroups.map((group) => String(group.id || '')));
                    local.groups = (Array.isArray(local.groups) ? local.groups : []).filter((group) => !overlappingIds.has(String(group.id || '')));
                    const newGroup = {{
                        id: 'cst_group_' + Date.now() + '_' + Math.floor(Math.random() * 1000000),
                        members: combinedEntries,
                    }};
                    local.groups.push(newGroup);
                    cleanupGroups(local);
                    setSelectedEntries(local, combinedEntries, combinedEntries[combinedEntries.length - 1] || null);
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const ungroupTargetEntry = (local, type, id) => {{
                    const groupIds = new Set(getGroupIdsForEntry(local, type, id));
                    if (!groupIds.size) return false;
                    const snapshot = buildSnapshot(local);
                    local.groups = (Array.isArray(local.groups) ? local.groups : []).filter((group) => !groupIds.has(String(group.id || '')));
                    cleanupGroups(local);
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const appendGroupingMenuItems = (local, targetType, targetId, addItem) => {{
                    if (typeof addItem !== 'function') return;
                    const selectedEntries = getSelectedEntries(local);
                    if (selectedEntries.length > 1) {{
                        const groupItem = addItem('Group');
                        groupItem.addEventListener('click', () => {{
                            createGroupFromSelection(local);
                            removeNodeContextMenu();
                            closeComplexMenu();
                            closeTextMenu();
                            closeEdgeBendMenu();
                        }});
                    }}
                    if (getGroupIdsForEntry(local, targetType, targetId).length) {{
                        const ungroupItem = addItem('Ungroup');
                        ungroupItem.addEventListener('click', () => {{
                            ungroupTargetEntry(local, targetType, targetId);
                            removeNodeContextMenu();
                            closeComplexMenu();
                            closeTextMenu();
                            closeEdgeBendMenu();
                        }});
                    }}
                }};
                const populateStaticGroupingMenu = (menu, local, targetType, targetId) => {{
                    if (!menu) return;
                    Array.from(menu.querySelectorAll('[data-cst-grouping-item="true"]')).forEach((item) => item.remove());
                    const selectedEntries = getSelectedEntries(local);
                    const canGroup = selectedEntries.length > 1;
                    const canUngroup = getGroupIdsForEntry(local, targetType, targetId).length > 0;
                    if (!canGroup && !canUngroup) return;
                    const makeButton = (label, handler) => {{
                        const button = document.createElement('button');
                        button.type = 'button';
                        button.textContent = label;
                        button.setAttribute('data-cst-grouping-item', 'true');
                        button.addEventListener('click', () => {{
                            try {{ handler(); }} catch (_err) {{}}
                        }});
                        return button;
                    }};
                    const firstButton = menu.querySelector('button');
                    if (canGroup) {{
                        const groupButton = makeButton('Group', () => {{
                            createGroupFromSelection(local);
                            closeComplexMenu();
                            closeTextMenu();
                        }});
                        menu.insertBefore(groupButton, firstButton || null);
                    }}
                    if (canUngroup) {{
                        const ungroupButton = makeButton('Ungroup', () => {{
                            ungroupTargetEntry(local, targetType, targetId);
                            closeComplexMenu();
                            closeTextMenu();
                        }});
                        menu.insertBefore(ungroupButton, firstButton || null);
                    }}
                }};
                const closeComplexMenu = () => {{
                    if (complexMenu) complexMenu.classList.remove('is-open');
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (local) local.contextNodeId = null;
                }};
                const closeEdgeBendMenu = () => {{
                    if (edgeBendMenu) edgeBendMenu.classList.remove('is-open');
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (local) {{
                        local.contextEdgeId = null;
                        local.contextBendIndex = null;
                    }}
                }};
                const closeTextMenu = () => {{
                    if (textMenu) textMenu.classList.remove('is-open');
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (local) local.contextTextNodeId = null;
                }};
                const buildTextMenuContents = (local, nodeId) => {{
                    if (!textMenu) return;
                    const node = findNodeById(local, nodeId);
                    if (!node) return;
                    textMenu.innerHTML = '';
                    const makeButton = (label, onClick) => {{
                        const button = document.createElement('button');
                        button.type = 'button';
                        button.textContent = label;
                        button.addEventListener('click', (evt) => {{
                            evt.preventDefault();
                            evt.stopPropagation();
                            if (typeof onClick === 'function') onClick();
                        }});
                        return button;
                    }};
                    const makeSubmenuItem = (label, title) => {{
                        const wrapper = document.createElement('div');
                        wrapper.className = 'cst-text-menu-item';
                        const button = makeButton(label, () => {{}});
                        wrapper.appendChild(button);
                        const submenu = document.createElement('div');
                        submenu.className = 'cst-text-menu-submenu';
                        const heading = document.createElement('div');
                        heading.className = 'cst-text-menu-label';
                        heading.textContent = title;
                        submenu.appendChild(heading);
                        wrapper.appendChild(submenu);
                        textMenu.appendChild(wrapper);
                        return submenu;
                    }};
                    const applyWithUndo = (mutator) => {{
                        if (typeof mutator !== 'function') return;
                        const snapshot = buildSnapshot(local);
                        mutator();
                        pushUndoSnapshot(local, snapshot);
                        renderEditableOverlay(local);
                        buildTextMenuContents(local, nodeId);
                        textMenu.classList.add('is-open');
                    }};

                    const fontSub = makeSubmenuItem('Font Size', 'Font Size');
                    const fontRow = document.createElement('div');
                    fontRow.className = 'cst-text-menu-row';
                    const fontInput = document.createElement('input');
                    fontInput.type = 'number';
                    fontInput.min = '8';
                    fontInput.max = '48';
                    fontInput.step = '1';
                    fontInput.value = String(Math.max(8, Number(node.fontSize || 11)));
                    const applyFontSize = () => {{
                        const nextSize = Number(fontInput.value);
                        if (!Number.isFinite(nextSize)) return;
                        applyWithUndo(() => {{
                            node.fontSize = Math.max(8, Math.min(48, nextSize));
                            const metrics = estimateTextNodeMetrics(String(node.displayLabel || node.label || 'Text'), node.fontSize || 11);
                            node.rx = metrics.rx;
                            node.ry = metrics.ry;
                        }});
                    }};
                    fontInput.addEventListener('change', applyFontSize);
                    fontInput.addEventListener('keydown', (evt) => {{
                        if (evt.key === 'Enter') {{
                            evt.preventDefault();
                            applyFontSize();
                        }}
                    }});
                    fontRow.appendChild(fontInput);
                    fontSub.appendChild(fontRow);

                    const alignSub = makeSubmenuItem('Alignment', 'Alignment');
                    const alignGrid = document.createElement('div');
                    alignGrid.className = 'cst-text-menu-option-grid';
                    const currentAlign = String(node.textAlign || 'center').toLowerCase();
                    ['left', 'center', 'right'].forEach((alignValue) => {{
                        const btn = makeButton(alignValue.charAt(0).toUpperCase() + alignValue.slice(1), () => {{
                            applyWithUndo(() => {{
                                node.textAlign = alignValue;
                            }});
                        }});
                        if (currentAlign === alignValue) btn.classList.add('is-active');
                        alignGrid.appendChild(btn);
                    }});
                    alignSub.appendChild(alignGrid);

                    const outlineSub = makeSubmenuItem('Outline', 'Outline');
                    const colorRow = document.createElement('div');
                    colorRow.className = 'cst-text-menu-row';
                    const colorInput = document.createElement('input');
                    colorInput.type = 'color';
                    colorInput.value = /^#/.test(String(node.stroke || '')) ? String(node.stroke) : '#475569';
                    const colorHex = document.createElement('input');
                    colorHex.type = 'text';
                    colorHex.maxLength = 7;
                    colorHex.value = String(node.stroke || '#475569');
                    const applyOutlineColor = (rawValue) => {{
                        let nextColor = String(rawValue || '').trim();
                        if (!nextColor) return;
                        if (!nextColor.startsWith('#')) nextColor = '#' + nextColor;
                        if (!/^#([0-9a-fA-F]{{3}}|[0-9a-fA-F]{{6}})$/.test(nextColor)) return;
                        applyWithUndo(() => {{
                            node.stroke = nextColor;
                        }});
                    }};
                    colorInput.addEventListener('input', () => applyOutlineColor(colorInput.value));
                    colorHex.addEventListener('change', () => applyOutlineColor(colorHex.value));
                    colorHex.addEventListener('keydown', (evt) => {{
                        if (evt.key === 'Enter') {{
                            evt.preventDefault();
                            applyOutlineColor(colorHex.value);
                        }}
                    }});
                    colorRow.append(colorInput, colorHex);
                    outlineSub.appendChild(colorRow);
                    const thickRow = document.createElement('div');
                    thickRow.className = 'cst-text-menu-row';
                    const thickLabel = document.createElement('div');
                    thickLabel.className = 'cst-text-menu-label';
                    thickLabel.style.marginBottom = '0';
                    thickLabel.textContent = 'Thickness';
                    const thickInput = document.createElement('input');
                    thickInput.type = 'number';
                    thickInput.min = '0';
                    thickInput.max = '12';
                    thickInput.step = '0.5';
                    thickInput.value = String(Math.max(0, Number(node.strokeWidth || 1.0)));
                    const applyThickness = () => {{
                        const nextWidth = Number(thickInput.value);
                        if (!Number.isFinite(nextWidth)) return;
                        applyWithUndo(() => {{
                            node.strokeWidth = Math.max(0, Math.min(12, nextWidth));
                        }});
                    }};
                    thickInput.addEventListener('change', applyThickness);
                    thickInput.addEventListener('keydown', (evt) => {{
                        if (evt.key === 'Enter') {{
                            evt.preventDefault();
                            applyThickness();
                        }}
                    }});
                    thickRow.append(thickLabel, thickInput);
                    outlineSub.appendChild(thickRow);

                    const boldButton = makeButton('Bold', () => {{
                        applyWithUndo(() => {{
                            node.fontWeight = String(node.fontWeight || '600') === '700' ? '500' : '700';
                        }});
                    }});
                    boldButton.className = 'cst-text-menu-toggle';
                    if (String(node.fontWeight || '600') === '700') boldButton.classList.add('is-active');
                    textMenu.appendChild(boldButton);

                    const deleteButton = makeButton('Delete', () => {{
                        closeTextMenu();
                        setSingleSelection(local, nodeId);
                        deleteSelectedNode(local);
                    }});
                    textMenu.appendChild(deleteButton);

                    populateStaticGroupingMenu(textMenu, local, 'node', nodeId);
                }};
                const openShapeMenu = (local, evt, nodeId) => {{
                    const node = findNodeById(local, nodeId);
                    if (!node) return;
                    removeNodeContextMenu();
                    const menu = document.createElement('div');
                    menu.className = 'cst-node-context-menu';
                    menu.style.position = 'absolute';
                    menu.style.background = '#fff';
                    menu.style.border = '1px solid #ccc';
                    menu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                    menu.style.padding = '2px 0';
                    menu.style.minWidth = '180px';
                    menu.style.zIndex = '1200';
                    const rect = stage.getBoundingClientRect();
                    menu.style.left = `${{Math.max(8, evt.clientX - rect.left + 4)}}px`;
                    menu.style.top = `${{Math.max(8, evt.clientY - rect.top - 24)}}px`;
                    const addItem = (label) => {{
                        const item = document.createElement('div');
                        item.textContent = label;
                        item.style.padding = '4px 8px';
                        item.style.fontSize = '12px';
                        item.style.cursor = 'pointer';
                        item.style.position = 'relative';
                        item.addEventListener('mouseenter', () => {{ item.style.backgroundColor = '#eef'; }});
                        item.addEventListener('mouseleave', () => {{ item.style.backgroundColor = 'transparent'; }});
                        menu.appendChild(item);
                        return item;
                    }};
                    const buildSubmenu = (parent) => {{
                        const submenu = document.createElement('div');
                        submenu.style.position = 'absolute';
                        submenu.style.left = '100%';
                        submenu.style.top = '0';
                        submenu.style.background = '#fff';
                        submenu.style.border = '1px solid #ccc';
                        submenu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                        submenu.style.minWidth = '180px';
                        submenu.style.padding = '6px';
                        submenu.style.display = 'none';
                        submenu.style.zIndex = '1201';
                        parent.addEventListener('mouseenter', () => {{ submenu.style.display = 'block'; }});
                        parent.addEventListener('mouseleave', () => {{ submenu.style.display = 'none'; }});
                        parent.appendChild(submenu);
                        return submenu;
                    }};
                    appendGroupingMenuItems(local, 'node', nodeId, addItem);
                    const applyShapeAppearance = () => {{
                        renderEditableOverlay(local);
                    }};
                    const liOutline = addItem('Outline');
                    const outlineSub = buildSubmenu(liOutline);
                    const ocRow = document.createElement('div');
                    ocRow.style.display = 'flex';
                    ocRow.style.alignItems = 'center';
                    ocRow.style.gap = '6px';
                    const ocInput = document.createElement('input');
                    ocInput.type = 'color';
                    ocInput.value = String(node.stroke || '#000000');
                    const ocHex = document.createElement('input');
                    ocHex.type = 'text';
                    ocHex.maxLength = 7;
                    ocHex.value = String(node.stroke || '#000000');
                    ocHex.style.width = '80px';
                    const applyOutlineColor = (val) => {{
                        node.stroke = String(val || '#000000');
                        applyShapeAppearance();
                    }};
                    const tryOcHex = () => {{
                        let val = String(ocHex.value || '').trim();
                        if (!val.startsWith('#')) val = '#' + val;
                        if (/^#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})$/.test(val)) {{
                            ocInput.value = val;
                            applyOutlineColor(val);
                        }}
                    }};
                    ocInput.addEventListener('input', () => {{
                        ocHex.value = ocInput.value;
                        applyOutlineColor(ocInput.value);
                    }});
                    ocHex.addEventListener('change', tryOcHex);
                    ocHex.addEventListener('keydown', (e) => {{ if (e.key === 'Enter') tryOcHex(); }});
                    ocRow.append(ocInput, ocHex);
                    outlineSub.appendChild(ocRow);
                    const thickRow = document.createElement('div');
                    thickRow.style.display = 'flex';
                    thickRow.style.alignItems = 'center';
                    thickRow.style.gap = '6px';
                    thickRow.style.marginTop = '8px';
                    const thickLabel = document.createElement('span');
                    thickLabel.textContent = 'Thickness';
                    thickLabel.style.fontSize = '11px';
                    const thickInput = document.createElement('input');
                    thickInput.type = 'number';
                    thickInput.min = '0';
                    thickInput.max = '20';
                    thickInput.step = '0.5';
                    thickInput.value = String(Number(node.strokeWidth || 1));
                    const applyThickness = () => {{
                        const val = Number(thickInput.value);
                        if (Number.isFinite(val)) {{
                            node.strokeWidth = Math.max(0, Math.min(20, val));
                            thickInput.value = String(node.strokeWidth);
                            applyShapeAppearance();
                        }}
                    }};
                    thickInput.addEventListener('change', applyThickness);
                    thickInput.addEventListener('keydown', (e) => {{ if (e.key === 'Enter') applyThickness(); }});
                    thickRow.append(thickLabel, thickInput);
                    outlineSub.appendChild(thickRow);
                    const liBg = addItem('Color');
                    const bgSub = buildSubmenu(liBg);
                    const bgRow = document.createElement('div');
                    bgRow.style.display = 'flex';
                    bgRow.style.alignItems = 'center';
                    bgRow.style.gap = '6px';
                    bgRow.style.marginBottom = '6px';
                    const bgInput = document.createElement('input');
                    bgInput.type = 'color';
                    bgInput.value = String(node.fillColor && node.fillColor !== 'transparent' ? node.fillColor : '#f5f5f5');
                    const bgHex = document.createElement('input');
                    bgHex.type = 'text';
                    bgHex.maxLength = 7;
                    bgHex.style.width = '80px';
                    bgHex.value = String(node.fillColor && node.fillColor !== 'transparent' ? node.fillColor : '#f5f5f5');
                    const applyBg = (val) => {{
                        node.fillColor = String(val || '#f5f5f5');
                        applyShapeAppearance();
                    }};
                    const tryBgHex = () => {{
                        let val = String(bgHex.value || '').trim();
                        if (!val.startsWith('#')) val = '#' + val;
                        if (/^#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})$/.test(val)) {{
                            bgInput.value = val;
                            applyBg(val);
                        }}
                    }};
                    bgInput.addEventListener('input', () => {{
                        bgHex.value = bgInput.value;
                        applyBg(bgInput.value);
                    }});
                    bgHex.addEventListener('change', tryBgHex);
                    bgHex.addEventListener('keydown', (e) => {{ if (e.key === 'Enter') tryBgHex(); }});
                    bgRow.append(bgInput, bgHex);
                    bgSub.appendChild(bgRow);
                    if (String(node.shapeType || '') !== 'bracket') {{
                        const noneBg = document.createElement('label');
                        noneBg.style.display = 'flex';
                        noneBg.style.alignItems = 'center';
                        noneBg.style.gap = '6px';
                        const noneRadio = document.createElement('input');
                        noneRadio.type = 'radio';
                        noneRadio.checked = String(node.fillColor || '') === 'transparent';
                        noneRadio.addEventListener('click', (e) => {{
                            e.stopPropagation();
                            node.fillColor = 'transparent';
                            applyShapeAppearance();
                            menu.remove();
                        }});
                        const noneText = document.createElement('span');
                        noneText.textContent = 'Transparent';
                        noneBg.append(noneRadio, noneText);
                        bgSub.appendChild(noneBg);
                    }}
                    const liDelete = addItem('Delete');
                    liDelete.addEventListener('click', () => {{
                        setSingleSelection(local, nodeId);
                        deleteSelectedNode(local);
                        removeNodeContextMenu();
                    }});
                    stage.appendChild(menu);
                    window.setTimeout(() => {{
                        const cleanup = (clickEvt) => {{
                            if (menu.contains(clickEvt.target)) return;
                            removeNodeContextMenu();
                            document.removeEventListener('click', cleanup, true);
                        }};
                        document.addEventListener('click', cleanup, true);
                    }}, 0);
                }};
                const removeNodeContextMenu = () => {{
                    const existing = stage.querySelector('.cst-node-context-menu');
                    if (existing) existing.remove();
                }};
                const getNodeProteinOptions = (node) => {{
                    const options = Array.isArray(node && node.proteinOptions) ? node.proteinOptions.slice() : [];
                    if (options.length) return options;
                    const fallback = [];
                    const gene = String(node && node.matchedGene || '').trim();
                    const uniprot = String(node && node.matchedUniprot || '').trim();
                    if (gene || uniprot) {{
                        fallback.push({{
                            gene_symbol: gene,
                            uniprot,
                            tooltip: String(node && node.tooltip || ''),
                            tooltip_html: String(node && (node.tooltipHtml || node.tooltip_html) || ''),
                        }});
                    }}
                    const genes = Array.isArray(node && node.suggestedGeneSymbols) ? node.suggestedGeneSymbols : [];
                    const uniprots = Array.isArray(node && node.candidateUniprotIds) ? node.candidateUniprotIds : [];
                    const count = Math.max(genes.length, uniprots.length);
                    for (let index = 0; index < count; index += 1) {{
                        const geneSymbol = String(genes[index] || '').trim();
                        const uniprotId = String(uniprots[index] || '').trim();
                        if (!geneSymbol && !uniprotId) continue;
                        if (fallback.some((item) => String(item.uniprot || '') === uniprotId && String(item.gene_symbol || '') === geneSymbol)) continue;
                        fallback.push({{ gene_symbol: geneSymbol, uniprot: uniprotId, tooltip: '', tooltip_html: '' }});
                    }}
                    return fallback;
                }};
                const applyNodeProteinOption = (local, nodeId, option) => {{
                    const node = findNodeById(local, nodeId);
                    if (!node || !option) return false;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    const nextGene = String(option.gene_symbol || '').trim();
                    const nextUniprot = String(option.uniprot || '').trim();
                    const nextMatchGene = String(option.match_gene_symbol || nextGene || '').trim();
                    node.matchedGene = nextGene;
                    node.matchedUniprot = nextUniprot;
                    const currentCandidates = Array.isArray(node.candidateUniprotIds) ? node.candidateUniprotIds : [];
                    node.candidateUniprotIds = [
                        ...(nextUniprot ? [nextUniprot] : []),
                        ...currentCandidates.filter((item) => String(item || '').trim() && String(item || '').trim() !== nextUniprot),
                    ];
                    const currentGenes = Array.isArray(node.suggestedGeneSymbols) ? node.suggestedGeneSymbols : [];
                    node.suggestedGeneSymbols = [
                        ...(nextMatchGene ? [nextMatchGene] : []),
                        ...currentGenes.filter((item) => String(item || '').trim() && String(item || '').trim() !== nextMatchGene),
                    ];
                    if (nextGene) {{
                        node.displayLabel = nextGene;
                        node.label = nextGene;
                    }} else if (nextUniprot) {{
                        node.displayLabel = nextUniprot;
                        node.label = nextUniprot;
                    }}
                    node.tooltip = String(option.tooltip || '');
                    node.tooltipHtml = String(option.tooltip_html || '');
                    node.title = node.tooltip || node.title || '';
                    const nextFold = option['fold_change_{active_idx}'];
                    if (Number.isFinite(Number(nextFold))) {{
                        node.foldText = Number(nextFold).toFixed(3);
                        node.hasDatasetMatch = true;
                    }} else {{
                        node.foldText = '';
                        node.hasDatasetMatch = false;
                    }}
                    const nextColor = option['fc_color_{active_idx}'] || option.fc_color_1 || null;
                    if (Array.isArray(nextColor) && nextColor.length >= 3) {{
                        const rgb = nextColor.slice(0, 3);
                        node.fillColor = 'rgb(' + rgb[0] + ', ' + rgb[1] + ', ' + rgb[2] + ')';
                    }}
                    for (const key of Object.keys(node || {{}})) {{
                        if (!/^(?:fc_color|fold_change|outline_color|outline_fold_change)_\\d+$/.test(key)) continue;
                        delete node[key];
                    }}
                    for (const key of Object.keys(option || {{}})) {{
                        if (!/^(?:fc_color|fold_change|outline_color|outline_fold_change)_\\d+$/.test(key)) continue;
                        node[key] = option[key];
                    }}
                    const nextOutlineColor = option['outline_color_{active_idx}'] || option.outline_color_1 || null;
                    if (Array.isArray(nextOutlineColor) && nextOutlineColor.length >= 3) {{
                        const outlineRgb = nextOutlineColor.slice(0, 3);
                        node['outline_color_{active_idx}'] = outlineRgb;
                        node.outline_color_1 = outlineRgb;
                        node.stroke = 'rgb(' + outlineRgb[0] + ', ' + outlineRgb[1] + ', ' + outlineRgb[2] + ')';
                    }} else if (!node.stroke) {{
                        node.stroke = '#000000';
                    }}
                    node.strokeWidth = Math.max(0.6, Number(local.proteinOutlineWidth || {prot_outline_width}));
                    if (AUTO_PTM_RENDER_ENABLED) {{
                        local.autoPtmPlacements = buildAutoPtmPlacements(local);
                    }}
                    renderEditableOverlay(local);
                    return true;
                }};
                const isNodeAutoPtmPlaced = (local, nodeId, ptmEntry) => {{
                    const targetKey = getAutoPtmEntryKey(ptmEntry);
                    if (!targetKey) return false;
                    for (const placement of (Array.isArray(local && local.autoPtmPlacements) ? local.autoPtmPlacements : [])) {{
                        if (String(placement && placement.nodeId || '') !== String(nodeId || '')) continue;
                        if (getAutoPtmEntryKey(placement) === targetKey) return true;
                    }}
                    return false;
                }};
                const ensureNodeAutoPtmPlacement = (local, nodeId, ptmEntry) => {{
                    if (!AUTO_PTM_RENDER_ENABLED) return false;
                    const node = findNodeById(local, nodeId);
                    if (!node || !ptmEntry) return false;
                    if (isNodeAutoPtmPlaced(local, nodeId, ptmEntry)) return true;
                    const ptmKey = getAutoPtmEntryKey(ptmEntry);
                    if (!ptmKey) return false;
                    const activeUniprotKey = getNodeAutoPtmSourceKey(node);
                    if (!activeUniprotKey) return false;
                    if (!Array.isArray(local.autoPtmPlacements)) {{
                        local.autoPtmPlacements = buildAutoPtmPlacements(local);
                    }}
                    const snapshot = buildSnapshot(local);
                    removeNodeHiddenAutoPtmKey(node, activeUniprotKey, ptmKey);
                    addNodeForcedAutoPtmKey(node, activeUniprotKey, ptmKey);
                    const existingPlacements = Array.isArray(local.autoPtmPlacements) ? local.autoPtmPlacements : [];
                    const placedPoints = [];
                    for (const placement of existingPlacements) {{
                        const ownerNode = findNodeById(local, placement && placement.nodeId);
                        if (!ownerNode) continue;
                        const resolved = resolveAutoPtmPlacementPoint(ownerNode, placement);
                        placedPoints.push({{
                            x: Number(resolved.x || 0),
                            y: Number(resolved.y || 0),
                            radius: Number(placement && placement.radius || AUTO_PTM_RADIUS),
                        }});
                    }}
                    const nextOrder = existingPlacements
                        .filter((placement) => String(placement && placement.nodeId || '') === String(nodeId || ''))
                        .reduce((maxOrder, placement) => Math.max(maxOrder, Number(placement && placement.order || 0)), 0) + 1;
                    const built = buildAutoPtmPlacementForEntry(local, node, ptmEntry, nextOrder, placedPoints, null, null, true);
                    if (!built || !built.placement) return false;
                    existingPlacements.push(built.placement);
                    local.autoPtmPlacements = existingPlacements;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const removeNodeAutoPtmPlacement = (local, nodeId, ptmId) => {{
                    const node = findNodeById(local, nodeId);
                    const placements = Array.isArray(local.autoPtmPlacements) ? local.autoPtmPlacements : [];
                    const placement = placements.find((item) => String(item && item.id || '') === String(ptmId || ''));
                    if (!node || !placement) return false;
                    const activeUniprotKey = getNodeAutoPtmSourceKey(node);
                    const ptmKey = getAutoPtmEntryKey(placement);
                    const snapshot = buildSnapshot(local);
                    local.autoPtmPlacements = placements.filter((item) => String(item && item.id || '') !== String(ptmId || ''));
                    if (activeUniprotKey && ptmKey) {{
                        removeNodeForcedAutoPtmKey(node, activeUniprotKey, ptmKey);
                        addNodeHiddenAutoPtmKey(node, activeUniprotKey, ptmKey);
                    }}
                    if (String(local.selectedPtmId || '') === String(ptmId || '')) {{
                        clearSelectedPtm(local);
                    }}
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const resetNodeAutoPtmLabelPosition = (local, nodeId, ptmId) => {{
                    const placement = Array.isArray(local.autoPtmPlacements)
                        ? local.autoPtmPlacements.find((item) => String(item && item.id || '') === String(ptmId || ''))
                        : null;
                    if (!placement) return false;
                    if (Math.abs(Number(placement.labelOffsetX || 0)) < 0.001 && Math.abs(Number(placement.labelOffsetY || 0)) < 0.001) return false;
                    const snapshot = buildSnapshot(local);
                    placement.labelOffsetX = 0;
                    placement.labelOffsetY = 0;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const showAutoPtmContextMenu = (local, evt, nodeId, ptmId) => {{
                    const node = findNodeById(local, nodeId);
                    const placement = Array.isArray(local.autoPtmPlacements)
                        ? local.autoPtmPlacements.find((item) => String(item && item.id || '') === String(ptmId || ''))
                        : null;
                    if (!node || !placement) return;
                    removeNodeContextMenu();
                    const menu = document.createElement('div');
                    menu.className = 'cst-node-context-menu';
                    menu.style.position = 'absolute';
                    menu.style.background = '#fff';
                    menu.style.border = '1px solid #ccc';
                    menu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                    menu.style.padding = '2px 0';
                    menu.style.minWidth = '170px';
                    menu.style.zIndex = '1200';
                    const rect = stage.getBoundingClientRect();
                    menu.style.left = `${{Math.max(8, evt.clientX - rect.left + 4)}}px`;
                    menu.style.top = `${{Math.max(8, evt.clientY - rect.top - 24)}}px`;
                    const addItem = (label, onClick) => {{
                        const item = document.createElement('div');
                        item.textContent = label;
                        item.style.padding = '4px 8px';
                        item.style.fontSize = '12px';
                        item.style.cursor = 'pointer';
                        item.addEventListener('mouseenter', () => {{ item.style.backgroundColor = '#eef'; }});
                        item.addEventListener('mouseleave', () => {{ item.style.backgroundColor = 'transparent'; }});
                        item.addEventListener('click', () => {{
                            onClick();
                            removeNodeContextMenu();
                        }});
                        menu.appendChild(item);
                        return item;
                    }};
                    const isKsActive = !!document.querySelector('#bookmark_selector a.nav-link.active[data-value=\"ks\"]');
                    if (isKsActive && typeof Shiny !== 'undefined' && Shiny && typeof Shiny.setInputValue === 'function') {{
                        addItem(\"Show Substrate's Kinases\", () => {{
                            try {{
                                Shiny.setInputValue('ks_spawn_ptm_kinases', {{
                                    protbox_id: String(nodeId || ''),
                                    uniprot: String(node.matchedUniprot || ''),
                                    ptm_key: String(placement.siteLabel || placement.label || ''),
                                    ts: Date.now(),
                                }}, {{ priority: 'event' }});
                            }} catch (_err) {{}}
                        }});
                    }}
                    addItem('Reset Label', () => {{
                        resetNodeAutoPtmLabelPosition(local, nodeId, ptmId);
                    }});
                    addItem('Remove', () => {{
                        removeNodeAutoPtmPlacement(local, nodeId, ptmId);
                    }});
                    stage.appendChild(menu);
                    window.setTimeout(() => {{
                        const cleanup = (clickEvt) => {{
                            if (menu.contains(clickEvt.target)) return;
                            removeNodeContextMenu();
                            document.removeEventListener('click', cleanup, true);
                        }};
                        document.addEventListener('click', cleanup, true);
                    }}, 0);
                }};
                const showNodeContextMenu = (local, evt, nodeId) => {{
                    const node = findNodeById(local, nodeId);
                    if (!node) return;
                    removeNodeContextMenu();
                    const menu = document.createElement('div');
                    menu.className = 'cst-node-context-menu';
                    menu.style.position = 'absolute';
                    menu.style.background = '#fff';
                    menu.style.border = '1px solid #ccc';
                    menu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                    menu.style.padding = '2px 0';
                    menu.style.minWidth = '150px';
                    menu.style.zIndex = '1200';
                    const rect = stage.getBoundingClientRect();
                    menu.style.left = `${{Math.max(8, evt.clientX - rect.left + 4)}}px`;
                    menu.style.top = `${{Math.max(8, evt.clientY - rect.top - 24)}}px`;
                    const addItem = (label) => {{
                        const item = document.createElement('div');
                        item.textContent = label;
                        item.style.padding = '4px 8px';
                        item.style.fontSize = '12px';
                        item.style.cursor = 'pointer';
                        item.style.position = 'relative';
                        item.addEventListener('mouseenter', () => {{ item.style.backgroundColor = '#eef'; }});
                        item.addEventListener('mouseleave', () => {{ item.style.backgroundColor = 'transparent'; }});
                        menu.appendChild(item);
                        return item;
                    }};
                    const buildSubmenu = (parent) => {{
                        const submenu = document.createElement('div');
                        submenu.style.position = 'absolute';
                        submenu.style.left = '100%';
                        submenu.style.top = '0';
                        submenu.style.background = '#fff';
                        submenu.style.border = '1px solid #ccc';
                        submenu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                        submenu.style.minWidth = '170px';
                        submenu.style.display = 'none';
                        submenu.style.zIndex = '1201';
                        parent.addEventListener('mouseenter', () => {{ submenu.style.display = 'block'; }});
                        parent.addEventListener('mouseleave', () => {{ submenu.style.display = 'none'; }});
                        parent.appendChild(submenu);
                        return submenu;
                    }};
                    appendGroupingMenuItems(local, 'node', nodeId, addItem);
                    const proteinItem = addItem('Protein');
                    const proteinSub = buildSubmenu(proteinItem);
                    const proteinOptions = getNodeProteinOptions(node);
                    if (!proteinOptions.length) {{
                        const empty = document.createElement('div');
                        empty.textContent = 'No proteins available';
                        empty.style.padding = '4px 8px';
                        empty.style.fontSize = '12px';
                        empty.style.color = '#666';
                        proteinSub.appendChild(empty);
                    }} else {{
                        proteinOptions.forEach((option) => {{
                            const sub = document.createElement('div');
                            sub.style.padding = '4px 8px';
                            sub.style.fontSize = '12px';
                            sub.style.cursor = 'pointer';
                            sub.style.display = 'flex';
                            sub.style.alignItems = 'center';
                            sub.style.gap = '8px';
                            const swatch = document.createElement('span');
                            const swatchColor = option['fc_color_{active_idx}'] || option.fc_color_1 || [166, 166, 166];
                            swatch.style.display = 'inline-block';
                            swatch.style.width = '12px';
                            swatch.style.height = '12px';
                            swatch.style.border = '1px solid #ccc';
                            swatch.style.flexShrink = '0';
                            if (Array.isArray(swatchColor) && swatchColor.length >= 3) {{
                                swatch.style.backgroundColor = `rgb(${{swatchColor[0]}}, ${{swatchColor[1]}}, ${{swatchColor[2]}})`;
                            }}
                            const labelSpan = document.createElement('span');
                            labelSpan.textContent = String(option.gene_symbol || option.uniprot || 'Protein');
                            sub.appendChild(swatch);
                            sub.appendChild(labelSpan);
                            if (String(option.uniprot || '') === String(node.matchedUniprot || '') && String(option.gene_symbol || '') === String(node.matchedGene || '')) {{
                                sub.style.fontWeight = 'bold';
                            }} else {{
                                sub.addEventListener('mouseenter', () => {{ sub.style.backgroundColor = '#eef'; }});
                                sub.addEventListener('mouseleave', () => {{ sub.style.backgroundColor = 'transparent'; }});
                            }}
                            sub.addEventListener('click', () => {{
                                applyNodeProteinOption(local, nodeId, option);
                                removeNodeContextMenu();
                            }});
                            proteinSub.appendChild(sub);
                        }});
                    }}
                    const ptmItem = addItem('PTMs');
                    const ptmSub = buildSubmenu(ptmItem);
                    const ptmEntries = getNodeAutoPtmEntries(node, {{ includeHidden: true }});
                    if (!ptmEntries.length) {{
                        const empty = document.createElement('div');
                        empty.textContent = 'No PTMs available';
                        empty.style.padding = '4px 8px';
                        empty.style.fontSize = '12px';
                        empty.style.color = '#666';
                        ptmSub.appendChild(empty);
                    }} else {{
                        ptmEntries.forEach((ptm) => {{
                            const sub = document.createElement('div');
                            sub.style.padding = '4px 8px';
                            sub.style.fontSize = '12px';
                            sub.style.cursor = 'pointer';
                            sub.style.display = 'flex';
                            sub.style.alignItems = 'center';
                            sub.style.gap = '8px';
                            const swatch = document.createElement('span');
                            const swatchColor = ptm['fc_color_{active_idx}'] || ptm.fc_color_1 || [209, 213, 219];
                            swatch.style.display = 'inline-flex';
                            swatch.style.alignItems = 'center';
                            swatch.style.justifyContent = 'center';
                            swatch.style.width = '14px';
                            swatch.style.height = '14px';
                            swatch.style.border = '1px solid #222';
                            swatch.style.borderRadius = String(ptm.shape || 'Circle').toLowerCase() === 'circle' ? '50%' : '3px';
                            swatch.style.fontSize = '9px';
                            swatch.style.fontWeight = 'bold';
                            if (Array.isArray(swatchColor) && swatchColor.length >= 3) {{
                                swatch.style.backgroundColor = `rgb(${{swatchColor[0]}}, ${{swatchColor[1]}}, ${{swatchColor[2]}})`;
                            }}
                            swatch.textContent = String(ptm.symbol || '');
                            const labelSpan = document.createElement('span');
                            labelSpan.textContent = String(ptm.site_label || ptm.label || 'PTM');
                            sub.appendChild(swatch);
                            sub.appendChild(labelSpan);
                            const isSpawned = isNodeAutoPtmPlaced(local, nodeId, ptm);
                            if (isSpawned) {{
                                sub.style.backgroundColor = '#d1d5db';
                                sub.style.cursor = 'not-allowed';
                            }} else {{
                                sub.addEventListener('mouseenter', () => {{ sub.style.backgroundColor = '#eef'; }});
                                sub.addEventListener('mouseleave', () => {{ sub.style.backgroundColor = 'transparent'; }});
                                sub.addEventListener('click', () => {{
                                    ensureNodeAutoPtmPlacement(local, nodeId, ptm);
                                    removeNodeContextMenu();
                                }});
                            }}
                            ptmSub.appendChild(sub);
                        }});
                    }}
                    const linksItem = addItem('Links');
                    const linksSub = buildSubmenu(linksItem);
                    const createLinkItem = (label, urlBuilder) => {{
                        const sub = document.createElement('div');
                        sub.textContent = label;
                        sub.style.padding = '4px 8px';
                        sub.style.fontSize = '12px';
                        sub.style.cursor = 'pointer';
                        const activeUniprot = String(node.matchedUniprot || '').trim();
                        if (!activeUniprot) {{
                            sub.style.color = '#666';
                            sub.style.cursor = 'default';
                        }} else {{
                            sub.addEventListener('mouseenter', () => {{ sub.style.backgroundColor = '#eef'; }});
                            sub.addEventListener('mouseleave', () => {{ sub.style.backgroundColor = 'transparent'; }});
                            sub.addEventListener('click', () => {{
                                const safeId = encodeURIComponent(activeUniprot);
                                window.open(urlBuilder(safeId), '_blank');
                                removeNodeContextMenu();
                            }});
                        }}
                        return sub;
                    }};
                    linksSub.appendChild(createLinkItem('UniProt', (safeId) => `https://www.uniprot.org/uniprotkb/${{safeId}}/entry`));
                    linksSub.appendChild(createLinkItem('PSP', (safeId) => `https://www.phosphosite.org/uniprotAccAction?id=${{safeId}}`));
                    stage.appendChild(menu);
                    window.setTimeout(() => {{
                        const cleanup = (clickEvt) => {{
                            if (menu.contains(clickEvt.target)) return;
                            removeNodeContextMenu();
                            document.removeEventListener('click', cleanup, true);
                        }};
                        document.addEventListener('click', cleanup, true);
                    }}, 0);
                }};
                const showEdgeContextMenu = (local, evt, edgeId) => {{
                    const edge = findEdgeById(local, edgeId);
                    if (!edge) return;
                    removeNodeContextMenu();
                    const menu = document.createElement('div');
                    menu.className = 'cst-node-context-menu';
                    menu.style.position = 'absolute';
                    menu.style.background = '#fff';
                    menu.style.border = '1px solid #ccc';
                    menu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                    menu.style.padding = '2px 0';
                    menu.style.minWidth = '150px';
                    menu.style.zIndex = '1200';
                    const rect = stage.getBoundingClientRect();
                    menu.style.left = `${{Math.max(8, evt.clientX - rect.left + 4)}}px`;
                    menu.style.top = `${{Math.max(8, evt.clientY - rect.top - 24)}}px`;
                    const addItem = (label) => {{
                        const item = document.createElement('div');
                        item.textContent = label;
                        item.style.padding = '4px 8px';
                        item.style.fontSize = '12px';
                        item.style.cursor = 'pointer';
                        item.addEventListener('mouseenter', () => {{ item.style.backgroundColor = '#eef'; }});
                        item.addEventListener('mouseleave', () => {{ item.style.backgroundColor = 'transparent'; }});
                        menu.appendChild(item);
                        return item;
                    }};
                    appendGroupingMenuItems(local, 'edge', edgeId, addItem);
                    if (!menu.childNodes.length) {{
                        return;
                    }}
                    stage.appendChild(menu);
                    window.setTimeout(() => {{
                        const cleanup = (clickEvt) => {{
                            if (menu.contains(clickEvt.target)) return;
                            removeNodeContextMenu();
                            document.removeEventListener('click', cleanup, true);
                        }};
                        document.addEventListener('click', cleanup, true);
                    }}, 0);
                }};
                const openComplexMenu = (local, evt, nodeId) => {{
                    if (!complexMenu) return;
                    closeTextMenu();
                    closeEdgeBendMenu();
                    removeNodeContextMenu();
                    populateStaticGroupingMenu(complexMenu, local, 'node', nodeId);
                    const rect = stage.getBoundingClientRect();
                    complexMenu.style.left = Math.max(8, evt.clientX - rect.left) + 'px';
                    complexMenu.style.top = Math.max(8, evt.clientY - rect.top) + 'px';
                    complexMenu.classList.add('is-open');
                    local.contextNodeId = nodeId;
                }};
                const openEdgeBendMenu = (local, evt, edgeId, bendIndex) => {{
                    if (!edgeBendMenu) return;
                    closeComplexMenu();
                    closeTextMenu();
                    removeNodeContextMenu();
                    const rect = stage.getBoundingClientRect();
                    edgeBendMenu.style.left = Math.max(8, evt.clientX - rect.left) + 'px';
                    edgeBendMenu.style.top = Math.max(8, evt.clientY - rect.top) + 'px';
                    edgeBendMenu.classList.add('is-open');
                    local.contextEdgeId = String(edgeId || '');
                    local.contextBendIndex = Number.isFinite(Number(bendIndex)) ? Number(bendIndex) : 0;
                }};
                const openTextMenu = (local, evt, nodeId) => {{
                    if (!textMenu) return;
                    closeComplexMenu();
                    closeEdgeBendMenu();
                    removeNodeContextMenu();
                    buildTextMenuContents(local, nodeId);
                    const rect = stage.getBoundingClientRect();
                    textMenu.style.left = Math.max(8, evt.clientX - rect.left) + 'px';
                    textMenu.style.top = Math.max(8, evt.clientY - rect.top) + 'px';
                    textMenu.classList.add('is-open');
                    local.contextTextNodeId = String(nodeId || '');
                }};
                const estimateComplexMemberWidth = (text) => {{
                    const label = String(text || '').trim();
                    return Math.max(38, Math.min(140, (label.length * 7.1) + 18));
                }};
                const buildComplexPanelLayout = (node) => {{
                    const members = Array.isArray(node.complexDisplayMembers) ? node.complexDisplayMembers : [];
                    if (!members.length) return null;
                    const headerHeight = 24;
                    const rowHeight = 15;
                    const rowGap = 6;
                    const colGap = 8;
                    const padding = 10;
                    const closeSize = 14;
                    const memberWidths = members.map((item) => estimateComplexMemberWidth(item.label || item.matched_gene_symbol || item.component_label || 'Component'));
                    const maxMemberWidth = Math.max(86, ...memberWidths);
                    const desiredCols = Math.max(1, Math.min(4, Math.ceil(Math.sqrt((members.length * 0.9)))));
                    const cols = Math.max(1, Math.min(desiredCols, members.length));
                    const rows = Math.max(1, Math.ceil(members.length / cols));
                    const contentWidth = (cols * maxMemberWidth) + (Math.max(0, cols - 1) * colGap);
                    const contentHeight = (rows * rowHeight) + (Math.max(0, rows - 1) * rowGap);
                    const panelWidth = contentWidth + (padding * 2);
                    const panelHeight = headerHeight + padding + contentHeight + padding;
                    return {{
                        panelWidth,
                        panelHeight,
                        contentWidth,
                        contentHeight,
                        padding,
                        headerHeight,
                        rowHeight,
                        rowGap,
                        colGap,
                        closeSize,
                        cols,
                        rows,
                        memberWidths,
                        maxMemberWidth,
                    }};
                }};
                const openComplexPanelForNode = (local, node) => {{
                    const layout = buildComplexPanelLayout(node);
                    if (!layout) return false;
                    local.complexPanels = Array.isArray(local.complexPanels) ? local.complexPanels.filter((panel) => panel.nodeId !== node.id) : [];
                    local.complexPanels.push({{
                        id: 'complex_panel_' + node.id,
                        nodeId: node.id,
                        x: Number(node.cx || 0) - (layout.panelWidth * 0.5),
                        y: Number(node.cy || 0) - (layout.panelHeight * 0.5),
                        width: layout.panelWidth,
                        height: layout.panelHeight,
                    }});
                    renderEditableOverlay(local);
                    return true;
                }};
                const closeComplexPanel = (local, panelId) => {{
                    if (!Array.isArray(local.complexPanels)) return false;
                    const before = local.complexPanels.length;
                    local.complexPanels = local.complexPanels.filter((panel) => panel.id !== panelId);
                    if (local.complexPanels.length === before) return false;
                    renderEditableOverlay(local);
                    return true;
                }};
                const isShapeUnderlayNode = (node) => {{
                    if (!node) return false;
                    const shapeType = String(node.shapeType || 'ellipse');
                    const mappingType = String(node.mappingType || '').toLowerCase();
                    if (shapeType === 'text') return false;
                    if (node.isCompoundNode || mappingType === 'metabolite_dataset') return false;
                    return !!node.isDrawingShape || mappingType === 'shape' || shapeType === 'legend' || shapeType === 'bracket';
                }};
                const renderEditableOverlay = (local) => {{
                    if (!overlaySvg) return;
                    overlaySvg.replaceChildren();
                    const defs = createSvg('defs');
                    const addLegendGradient = (id, isHorizontal) => {{
                        const gradient = createSvg('linearGradient');
                        gradient.setAttribute('id', id);
                        gradient.setAttribute('x1', '0%');
                        gradient.setAttribute('y1', '0%');
                        gradient.setAttribute('x2', isHorizontal ? '100%' : '0%');
                        gradient.setAttribute('y2', isHorizontal ? '0%' : '100%');
                        const rawStops = Array.isArray(legendConfig.stops) ? legendConfig.stops : [];
                        let stops = rawStops
                            .map((entry) => {{
                                const value = Number(entry && entry.value);
                                const color = Array.isArray(entry && entry.color) ? entry.color.slice(0, 3).map((c) => Number(c)) : null;
                                if (!Number.isFinite(value) || !color || color.some((c) => !Number.isFinite(c))) return null;
                                return {{ value, color }};
                            }})
                            .filter((entry) => !!entry)
                            .sort((a, b) => a.value - b.value);
                        if (stops.length < 2) {{
                            const neg = Array.isArray(legendConfig.negColor) ? legendConfig.negColor : [0, 114, 178];
                            const pos = Array.isArray(legendConfig.posColor) ? legendConfig.posColor : [0, 158, 115];
                            stops = [
                                {{ value: Number(legendConfig.maxNeg ?? -2), color: neg }},
                                {{ value: 0, color: [255, 255, 255] }},
                                {{ value: Number(legendConfig.maxPos ?? 2), color: pos }},
                            ].sort((a, b) => a.value - b.value);
                        }}
                        const minValue = Number(stops[0].value);
                        const maxValue = Number(stops[stops.length - 1].value);
                        const span = (maxValue - minValue) || 1;
                        stops.forEach((entry) => {{
                            const rawOffset = (Number(entry.value) - minValue) / span;
                            const normalized = Math.max(0, Math.min(1, isHorizontal ? rawOffset : (1 - rawOffset)));
                            const offset = `${{(normalized * 100).toFixed(3)}}%`;
                            const color = `rgb(${{Math.max(0, Math.min(255, Math.round(entry.color[0]))) }}, ${{Math.max(0, Math.min(255, Math.round(entry.color[1]))) }}, ${{Math.max(0, Math.min(255, Math.round(entry.color[2]))) }})`;
                            const stop = createSvg('stop');
                            stop.setAttribute('offset', offset);
                            stop.setAttribute('stop-color', color);
                            gradient.appendChild(stop);
                        }});
                        defs.appendChild(gradient);
                    }};
                    addLegendGradient('cst-legend-grad-v-{viewer_key}', false);
                    addLegendGradient('cst-legend-grad-h-{viewer_key}', true);
                    overlaySvg.appendChild(defs);
                    const temporalModeEnabled = !!local.temporalMode;
                    const rgbCss = (value, fallback = 'rgb(0, 0, 0)') => {{
                        return Array.isArray(value) && value.length >= 3
                            ? ('rgb(' + Number(value[0] || 0) + ', ' + Number(value[1] || 0) + ', ' + Number(value[2] || 0) + ')')
                            : fallback;
                    }};
                    const getNodeTemporalSegmentCount = (node) => {{
                        let maxIndex = 0;
                        for (const key of Object.keys(node || {{}})) {{
                            const match = String(key || '').match(/^(?:fc_color|fold_change)_([0-9]+)$/);
                            if (!match) continue;
                            const idx = Number(match[1]);
                            if (Number.isFinite(idx) && idx > maxIndex) maxIndex = idx;
                        }}
                        return Math.max(1, maxIndex);
                    }};
                    const hasSegmentedOutlineData = (node, segmentCount) => {{
                        for (let idx = 1; idx <= segmentCount; idx += 1) {{
                            const value = node ? node['outline_fold_change_' + idx] : null;
                            if (value === 0 || Number.isFinite(Number(value))) return true;
                        }}
                        return false;
                    }};
                    const applyNodeGeometry = (shape, node) => {{
                        const nodeShapeType = String(node.shapeType || 'ellipse');
                        const nodeIsRect = nodeShapeType === 'rect' || nodeShapeType === 'square';
                        const nodeIsRounded = nodeShapeType === 'rounded';
                        const nodeIsLegend = nodeShapeType === 'legend';
                        const nodeIsRectLike = nodeIsRect || nodeIsRounded || nodeIsLegend;
                        const nodeIsBracket = nodeShapeType === 'bracket';
                        const nodeIsText = nodeShapeType === 'text';
                        if (nodeIsBracket) {{
                            shape.setAttribute('d', buildBracketPath(node));
                        }} else if (nodeIsRectLike || nodeIsText) {{
                            shape.setAttribute('x', (Number(node.cx || 0) - Number(node.rx || 12)).toFixed(3));
                            shape.setAttribute('y', (Number(node.cy || 0) - Number(node.ry || 9)).toFixed(3));
                            shape.setAttribute('width', (Number(node.rx || 12) * 2).toFixed(3));
                            shape.setAttribute('height', (Number(node.ry || 9) * 2).toFixed(3));
                            if (nodeIsText) {{
                                shape.setAttribute('rx', '8');
                                shape.setAttribute('ry', '8');
                            }} else if (nodeIsLegend) {{
                                shape.setAttribute('rx', '0');
                                shape.setAttribute('ry', '0');
                            }} else if (nodeIsRounded) {{
                                shape.setAttribute('rx', '10');
                                shape.setAttribute('ry', '10');
                            }}
                        }} else {{
                            shape.setAttribute('cx', Number(node.cx || 0).toFixed(3));
                            shape.setAttribute('cy', Number(node.cy || 0).toFixed(3));
                            shape.setAttribute('rx', Number(node.rx || 12).toFixed(3));
                            shape.setAttribute('ry', Number(node.ry || 9).toFixed(3));
                        }}
                        if (Number(node.angle || 0)) {{
                            shape.setAttribute('transform', 'rotate(' + Number(node.angle || 0).toFixed(3) + ' ' + Number(node.cx || 0).toFixed(3) + ' ' + Number(node.cy || 0).toFixed(3) + ')');
                        }}
                    }};
                    const appendTemporalNodeSegments = (node, group, bodyShape) => {{
                        if (!temporalModeEnabled || !local.proteinOvalMode || local.tabPreviewMode) return false;
                        if (!node || node.isComplex || node.isDrawingShape) return false;
                        const nodeShapeType = String(node.shapeType || 'ellipse');
                        if (nodeShapeType === 'text' || nodeShapeType === 'legend' || nodeShapeType === 'bracket') return false;
                        const segmentCount = getNodeTemporalSegmentCount(node);
                        if (segmentCount <= 1) return false;
                        const strokeWidth = Math.max(0.6, Number(node.strokeWidth || local.proteinOutlineWidth || {prot_outline_width}));
                        const bboxX = Number(node.cx || 0) - Number(node.rx || 12);
                        const bboxY = Number(node.cy || 0) - Number(node.ry || 9);
                        const bboxWidth = Number(node.rx || 12) * 2;
                        const bboxHeight = Number(node.ry || 9) * 2;
                        const hasSegmentedOutline = hasSegmentedOutlineData(node, segmentCount);
                        for (let idx = 1; idx <= segmentCount; idx += 1) {{
                            const segmentStart = bboxX + (((idx - 1) * bboxWidth) / segmentCount);
                            const segmentEnd = idx === segmentCount ? (bboxX + bboxWidth) : (bboxX + ((idx * bboxWidth) / segmentCount));
                            const fillClipId = 'cst_temporal_clip_{viewer_key}_' + String(node.id || '') + '_' + String(idx);
                            const fillClipPath = createSvg('clipPath');
                            fillClipPath.setAttribute('id', fillClipId);
                            fillClipPath.setAttribute('clipPathUnits', 'userSpaceOnUse');
                            const fillClipRect = createSvg('rect');
                            fillClipRect.setAttribute('x', segmentStart.toFixed(3));
                            fillClipRect.setAttribute('y', bboxY.toFixed(3));
                            fillClipRect.setAttribute('width', Math.max(0.5, (segmentEnd - segmentStart)).toFixed(3));
                            fillClipRect.setAttribute('height', Math.max(0.5, bboxHeight).toFixed(3));
                            fillClipPath.appendChild(fillClipRect);
                            defs.appendChild(fillClipPath);
                            const fillShape = createSvg(bodyShape.tagName.toLowerCase());
                            applyNodeGeometry(fillShape, node);
                            fillShape.setAttribute('fill', rgbCss(node['fc_color_' + idx], String(node.fillColor || 'rgb(166, 166, 166)')));
                            fillShape.setAttribute('stroke', 'none');
                            fillShape.setAttribute('pointer-events', 'none');
                            fillShape.setAttribute('clip-path', 'url(#' + fillClipId + ')');
                            fillShape.setAttribute('opacity', Number(node.opacity || 1.0).toFixed(2));
                            group.appendChild(fillShape);
                            if (hasSegmentedOutline) {{
                                const outlineClipId = 'cst_temporal_outline_clip_{viewer_key}_' + String(node.id || '') + '_' + String(idx);
                                const outlineClipPath = createSvg('clipPath');
                                outlineClipPath.setAttribute('id', outlineClipId);
                                outlineClipPath.setAttribute('clipPathUnits', 'userSpaceOnUse');
                                const outlineClipRect = createSvg('rect');
                                outlineClipRect.setAttribute('x', segmentStart.toFixed(3));
                                outlineClipRect.setAttribute('y', bboxY.toFixed(3));
                                outlineClipRect.setAttribute('width', Math.max(0.5, (segmentEnd - segmentStart)).toFixed(3));
                                outlineClipRect.setAttribute('height', Math.max(0.5, bboxHeight).toFixed(3));
                                outlineClipPath.appendChild(outlineClipRect);
                                defs.appendChild(outlineClipPath);
                                const outlineShape = createSvg(bodyShape.tagName.toLowerCase());
                                applyNodeGeometry(outlineShape, node);
                                outlineShape.setAttribute('fill', 'none');
                                outlineShape.setAttribute('stroke', rgbCss(node['outline_color_' + idx], String(node.stroke || '#000000')));
                                outlineShape.setAttribute('stroke-width', strokeWidth.toFixed(3));
                                outlineShape.setAttribute('stroke-linecap', 'round');
                                outlineShape.setAttribute('stroke-linejoin', 'round');
                                outlineShape.setAttribute('pointer-events', 'none');
                                outlineShape.setAttribute('clip-path', 'url(#' + outlineClipId + ')');
                                outlineShape.setAttribute('opacity', Number(node.opacity || 1.0).toFixed(2));
                                group.appendChild(outlineShape);
                            }}
                        }}
                        bodyShape.setAttribute('fill', 'transparent');
                        bodyShape.setAttribute('stroke', hasSegmentedOutline ? 'transparent' : String(node.stroke || '#000000'));
                        bodyShape.setAttribute('stroke-width', hasSegmentedOutline ? '0' : strokeWidth.toFixed(3));
                        bodyShape.setAttribute('pointer-events', 'all');
                        bodyShape.style.pointerEvents = 'all';
                        return true;
                    }};
                    const getPtmTemporalSegmentCount = (ptm) => {{
                        let maxIndex = 0;
                        for (const key of Object.keys(ptm || {{}})) {{
                            const match = String(key || '').match(/^(?:fc_color|fold_change)_([0-9]+)$/);
                            if (!match) continue;
                            const idx = Number(match[1]);
                            if (Number.isFinite(idx) && idx > maxIndex) maxIndex = idx;
                        }}
                        return Math.max(1, maxIndex);
                    }};
                    const hasPtmSegmentedOutlineData = (ptm, segmentCount) => {{
                        for (let idx = 1; idx <= segmentCount; idx += 1) {{
                            const foldValue = ptm ? ptm['outline_fold_change_' + idx] : null;
                            if (foldValue === 0 || Number.isFinite(Number(foldValue))) return true;
                            if (idx <= 1) continue;
                            const colorKey = 'outline_color_' + idx;
                            const colorValue = ptm ? ptm[colorKey] : null;
                            if (ptm && Object.prototype.hasOwnProperty.call(ptm, colorKey) && Array.isArray(colorValue) && colorValue.length >= 3) return true;
                        }}
                        return false;
                    }};
                    const getPtmTemporalRenderSpec = (ptm, ptmPoint, ptmRadius) => {{
                        const visualRadius = Math.max(2.2, Number(ptmRadius || 0) * 0.78);
                        const nativeShape = String(ptm && ptm.shape || 'Circle').toLowerCase();
                        const segmentCount = temporalModeEnabled ? getPtmTemporalSegmentCount(ptm) : 1;
                        const renderShape = segmentCount >= 3 ? 'pill' : nativeShape;
                        const width = renderShape === 'pill'
                            ? Math.max(visualRadius * 3.35, segmentCount * visualRadius * 1.3)
                            : (visualRadius * 2);
                        const height = visualRadius * 2;
                        return {{
                            nativeShape,
                            renderShape,
                            segmentCount,
                            visualRadius,
                            width,
                            height,
                            left: Number(ptmPoint.x || 0) - (width * 0.5),
                            top: Number(ptmPoint.y || 0) - visualRadius,
                            centerX: Number(ptmPoint.x || 0),
                            centerY: Number(ptmPoint.y || 0),
                        }};
                    }};
                    const applyPtmShapeGeometry = (shape, spec) => {{
                        if (spec.renderShape === 'circle') {{
                            shape.setAttribute('cx', spec.centerX.toFixed(3));
                            shape.setAttribute('cy', spec.centerY.toFixed(3));
                            shape.setAttribute('r', (spec.height * 0.5).toFixed(3));
                            shape.removeAttribute('transform');
                            return;
                        }}
                        shape.setAttribute('x', spec.left.toFixed(3));
                        shape.setAttribute('y', spec.top.toFixed(3));
                        shape.setAttribute('width', spec.width.toFixed(3));
                        shape.setAttribute('height', spec.height.toFixed(3));
                        if (spec.renderShape === 'pill') {{
                            shape.setAttribute('rx', (spec.height * 0.5).toFixed(3));
                            shape.setAttribute('ry', (spec.height * 0.5).toFixed(3));
                            shape.removeAttribute('transform');
                        }} else {{
                            shape.removeAttribute('rx');
                            shape.removeAttribute('ry');
                            if (spec.renderShape === 'diamond') {{
                                shape.setAttribute('transform', 'rotate(45 ' + spec.centerX.toFixed(3) + ' ' + spec.centerY.toFixed(3) + ')');
                            }} else {{
                                shape.removeAttribute('transform');
                            }}
                        }}
                    }};
                    const makeCstCirclePtmSegmentPath = (centerX, centerY, radius, side) => {{
                        const topY = centerY - radius;
                        const bottomY = centerY + radius;
                        const sweep = side === 'left' ? 0 : 1;
                        return 'M ' + centerX.toFixed(3) + ' ' + topY.toFixed(3)
                            + ' A ' + radius.toFixed(3) + ' ' + radius.toFixed(3) + ' 0 0 ' + sweep + ' ' + centerX.toFixed(3) + ' ' + bottomY.toFixed(3)
                            + ' L ' + centerX.toFixed(3) + ' ' + topY.toFixed(3) + ' Z';
                    }};
                    const makeCstPillEndSegmentPath = (left, top, width, height, side) => {{
                        const radius = height * 0.5;
                        const right = left + width;
                        const bottom = top + height;
                        const centerLeft = left + radius;
                        const centerRight = left + width - radius;
                        if (side === 'left') {{
                            return 'M ' + right.toFixed(3) + ' ' + top.toFixed(3)
                                + ' L ' + centerLeft.toFixed(3) + ' ' + top.toFixed(3)
                                + ' A ' + radius.toFixed(3) + ' ' + radius.toFixed(3) + ' 0 0 0 ' + centerLeft.toFixed(3) + ' ' + bottom.toFixed(3)
                                + ' L ' + right.toFixed(3) + ' ' + bottom.toFixed(3) + ' Z';
                        }}
                        return 'M ' + left.toFixed(3) + ' ' + top.toFixed(3)
                            + ' L ' + centerRight.toFixed(3) + ' ' + top.toFixed(3)
                            + ' A ' + radius.toFixed(3) + ' ' + radius.toFixed(3) + ' 0 0 1 ' + centerRight.toFixed(3) + ' ' + bottom.toFixed(3)
                            + ' L ' + left.toFixed(3) + ' ' + bottom.toFixed(3) + ' Z';
                    }};
                    const createPtmSegmentedBody = (local, group, ptm, ptmPoint, ptmRadius, bodyShape) => {{
                        if (!temporalModeEnabled || !ptm || !bodyShape) return null;
                        const spec = getPtmTemporalRenderSpec(ptm, ptmPoint, ptmRadius);
                        const segmentCount = spec.segmentCount;
                        const strokeWidth = Math.max(0.6, Number(ptm.outlineWidth || 1.0));
                        const opacity = Number(ptm.opacity || 1.0).toFixed(2);
                        const fallbackFill = rgbCss(ptm['fc_color_{active_idx}'] || ptm.fc_color_1, 'rgb(209, 213, 219)');
                        const fallbackOutline = rgbCss(ptm['outline_color_{active_idx}'] || ptm.outline_color_1, 'rgb(0, 0, 0)');
                        const hasSegmentedOutline = segmentCount > 1 && hasPtmSegmentedOutlineData(ptm, segmentCount);
                        const makeSegmentClipId = (idx) => 'cst_temporal_ptm_clip_{viewer_key}_' + String(ptm.id || '') + '_' + String(idx);
                        const segmentBounds = [];
                        if (segmentCount > 1) {{
                            for (let idx = 1; idx <= segmentCount; idx += 1) {{
                                const segmentStart = spec.left + (((idx - 1) * spec.width) / segmentCount);
                                const segmentEnd = idx === segmentCount ? (spec.left + spec.width) : (spec.left + ((idx * spec.width) / segmentCount));
                                segmentBounds.push({{ idx, x1: segmentStart, x2: segmentEnd }});
                                let fillShape = null;
                                if (spec.renderShape === 'circle' && segmentCount === 2) {{
                                    fillShape = createSvg('path');
                                    fillShape.setAttribute('d', makeCstCirclePtmSegmentPath(spec.centerX, spec.centerY, spec.height * 0.5, idx === 1 ? 'left' : 'right'));
                                }} else if (spec.renderShape === 'pill') {{
                                    const actualWidth = Math.max(0.5, segmentEnd - segmentStart);
                                    if (idx === 1) {{
                                        fillShape = createSvg('path');
                                        fillShape.setAttribute('d', makeCstPillEndSegmentPath(segmentStart, spec.top, actualWidth, spec.height, 'left'));
                                    }} else if (idx === segmentCount) {{
                                        fillShape = createSvg('path');
                                        fillShape.setAttribute('d', makeCstPillEndSegmentPath(segmentStart, spec.top, actualWidth, spec.height, 'right'));
                                    }} else {{
                                        fillShape = createSvg('rect');
                                        fillShape.setAttribute('x', segmentStart.toFixed(3));
                                        fillShape.setAttribute('y', spec.top.toFixed(3));
                                        fillShape.setAttribute('width', actualWidth.toFixed(3));
                                        fillShape.setAttribute('height', spec.height.toFixed(3));
                                    }}
                                }} else {{
                                    const clipPath = createSvg('clipPath');
                                    clipPath.setAttribute('id', makeSegmentClipId(idx));
                                    const clipRect = createSvg('rect');
                                    clipRect.setAttribute('x', (segmentStart - strokeWidth - 1).toFixed(3));
                                    clipRect.setAttribute('y', (spec.top - strokeWidth - 1).toFixed(3));
                                    clipRect.setAttribute('width', Math.max(0.5, (segmentEnd - segmentStart) + (strokeWidth * 2) + 2).toFixed(3));
                                    clipRect.setAttribute('height', Math.max(0.5, spec.height + (strokeWidth * 2) + 2).toFixed(3));
                                    clipPath.appendChild(clipRect);
                                    defs.appendChild(clipPath);
                                    fillShape = createSvg(bodyShape.tagName.toLowerCase());
                                    applyPtmShapeGeometry(fillShape, spec);
                                    fillShape.setAttribute('clip-path', 'url(#' + makeSegmentClipId(idx) + ')');
                                }}
                                if (!fillShape) continue;
                                fillShape.setAttribute('fill', rgbCss(ptm['fc_color_' + idx], fallbackFill));
                                fillShape.setAttribute('stroke', 'none');
                                fillShape.setAttribute('pointer-events', 'none');
                                fillShape.setAttribute('opacity', opacity);
                                group.appendChild(fillShape);
                            }}
                        }}
                        if (hasSegmentedOutline) {{
                            if (spec.renderShape === 'pill') {{
                                const radius = spec.height * 0.5;
                                const straightLeft = spec.left + radius;
                                const straightRight = spec.left + spec.width - radius;
                                segmentBounds.forEach((segment) => {{
                                    const color = rgbCss(ptm['outline_color_' + segment.idx], fallbackOutline);
                                    const topStart = Math.max(straightLeft, segment.x1);
                                    const topEnd = Math.min(straightRight, segment.x2);
                                    if (topEnd > topStart) {{
                                        const topLine = createSvg('line');
                                        topLine.setAttribute('x1', topStart.toFixed(3));
                                        topLine.setAttribute('y1', spec.top.toFixed(3));
                                        topLine.setAttribute('x2', topEnd.toFixed(3));
                                        topLine.setAttribute('y2', spec.top.toFixed(3));
                                        topLine.setAttribute('stroke', color);
                                        topLine.setAttribute('stroke-width', strokeWidth.toFixed(3));
                                        topLine.setAttribute('stroke-linecap', 'butt');
                                        topLine.setAttribute('stroke-linejoin', 'round');
                                        topLine.setAttribute('fill', 'none');
                                        topLine.setAttribute('pointer-events', 'none');
                                        topLine.setAttribute('opacity', opacity);
                                        group.appendChild(topLine);
                                        const bottomLine = createSvg('line');
                                        bottomLine.setAttribute('x1', topStart.toFixed(3));
                                        bottomLine.setAttribute('y1', (spec.top + spec.height).toFixed(3));
                                        bottomLine.setAttribute('x2', topEnd.toFixed(3));
                                        bottomLine.setAttribute('y2', (spec.top + spec.height).toFixed(3));
                                        bottomLine.setAttribute('stroke', color);
                                        bottomLine.setAttribute('stroke-width', strokeWidth.toFixed(3));
                                        bottomLine.setAttribute('stroke-linecap', 'butt');
                                        bottomLine.setAttribute('stroke-linejoin', 'round');
                                        bottomLine.setAttribute('fill', 'none');
                                        bottomLine.setAttribute('pointer-events', 'none');
                                        bottomLine.setAttribute('opacity', opacity);
                                        group.appendChild(bottomLine);
                                    }}
                                }});
                                const leftArc = createSvg('path');
                                leftArc.setAttribute('d', 'M ' + (spec.left + radius).toFixed(3) + ' ' + spec.top.toFixed(3) + ' A ' + radius.toFixed(3) + ' ' + radius.toFixed(3) + ' 0 0 0 ' + (spec.left + radius).toFixed(3) + ' ' + (spec.top + spec.height).toFixed(3));
                                leftArc.setAttribute('fill', 'none');
                                leftArc.setAttribute('stroke', rgbCss(ptm['outline_color_1'], fallbackOutline));
                                leftArc.setAttribute('stroke-width', strokeWidth.toFixed(3));
                                leftArc.setAttribute('stroke-linecap', 'butt');
                                leftArc.setAttribute('stroke-linejoin', 'round');
                                leftArc.setAttribute('pointer-events', 'none');
                                leftArc.setAttribute('opacity', opacity);
                                group.appendChild(leftArc);
                                const rightArc = createSvg('path');
                                rightArc.setAttribute('d', 'M ' + (spec.left + spec.width - radius).toFixed(3) + ' ' + spec.top.toFixed(3) + ' A ' + radius.toFixed(3) + ' ' + radius.toFixed(3) + ' 0 0 1 ' + (spec.left + spec.width - radius).toFixed(3) + ' ' + (spec.top + spec.height).toFixed(3));
                                rightArc.setAttribute('fill', 'none');
                                rightArc.setAttribute('stroke', rgbCss(ptm['outline_color_' + segmentCount], fallbackOutline));
                                rightArc.setAttribute('stroke-width', strokeWidth.toFixed(3));
                                rightArc.setAttribute('stroke-linecap', 'butt');
                                rightArc.setAttribute('stroke-linejoin', 'round');
                                rightArc.setAttribute('pointer-events', 'none');
                                rightArc.setAttribute('opacity', opacity);
                                group.appendChild(rightArc);
                            }} else {{
                                for (let idx = 1; idx <= segmentCount; idx += 1) {{
                                    const outlineShape = createSvg(bodyShape.tagName.toLowerCase());
                                    applyPtmShapeGeometry(outlineShape, spec);
                                    outlineShape.setAttribute('fill', 'none');
                                    outlineShape.setAttribute('stroke', rgbCss(ptm['outline_color_' + idx], fallbackOutline));
                                    outlineShape.setAttribute('stroke-width', strokeWidth.toFixed(3));
                                    outlineShape.setAttribute('stroke-linecap', 'round');
                                    outlineShape.setAttribute('stroke-linejoin', 'round');
                                    outlineShape.setAttribute('pointer-events', 'none');
                                    outlineShape.setAttribute('clip-path', 'url(#' + makeSegmentClipId(idx) + ')');
                                    outlineShape.setAttribute('opacity', opacity);
                                    group.appendChild(outlineShape);
                                }}
                            }}
                        }}
                        bodyShape.setAttribute('fill', segmentCount > 1 ? 'transparent' : fallbackFill);
                        bodyShape.setAttribute('stroke', hasSegmentedOutline ? 'transparent' : fallbackOutline);
                        bodyShape.setAttribute('stroke-width', hasSegmentedOutline ? '0' : strokeWidth.toFixed(3));
                        bodyShape.setAttribute('pointer-events', 'all');
                        bodyShape.style.pointerEvents = 'all';
                        bodyShape.setAttribute('opacity', opacity);
                        return spec;
                    }};
                    const underlayFragment = document.createDocumentFragment();
                    const fragment = document.createDocumentFragment();
                    const autoPtmPlacements = (!AUTO_PTM_RENDER_ENABLED || local.tabPreviewMode) ? [] : (Array.isArray(local.autoPtmPlacements) ? local.autoPtmPlacements : []);
                    const autoPtmsByNode = new Map();
                    for (const placement of autoPtmPlacements) {{
                        const key = String(placement.nodeId || '');
                        if (!autoPtmsByNode.has(key)) autoPtmsByNode.set(key, []);
                        autoPtmsByNode.get(key).push(placement);
                    }}
                    if (!local.tabPreviewMode) for (const edge of (local.editableEdges || [])) {{
                        try {{
                            if (local.showArrows === false) {{
                                continue;
                            }}
                            const edgeDisplayStroke = String(edge.stroke || '#0f172a');
                            const edgeGroup = createSvg('g');
                            edgeGroup.setAttribute('class', 'cst-edit-edge-group');
                            edgeGroup.setAttribute('data-edge-id', edge.id);
                            const edgePath = createSvg('path');
                            edgePath.setAttribute('d', buildEdgePath(edge));
                            edgePath.setAttribute('fill', 'none');
                            edgePath.setAttribute('stroke', edgeDisplayStroke);
                            edgePath.setAttribute('stroke-width', Number(edge.strokeWidth || 3).toFixed(3));
                            edgePath.setAttribute('stroke-linecap', 'round');
                            edgePath.setAttribute('stroke-linejoin', 'round');
                            if (edge.dashed) edgePath.setAttribute('stroke-dasharray', '7 5');
                            edgePath.setAttribute('opacity', Number(edge.opacity || 0.95).toFixed(2));
                            edgePath.setAttribute('pointer-events', 'none');
                            edgeGroup.appendChild(edgePath);
                        if (String(edge.startType || 'none') === 'inhibitor') {{
                            const startBar = getStartInhibitorBarPoints(edge);
                            const startInhibitorShape = createSvg('line');
                            startInhibitorShape.setAttribute('x1', Number(startBar.left.x || 0).toFixed(3));
                            startInhibitorShape.setAttribute('y1', Number(startBar.left.y || 0).toFixed(3));
                            startInhibitorShape.setAttribute('x2', Number(startBar.right.x || 0).toFixed(3));
                            startInhibitorShape.setAttribute('y2', Number(startBar.right.y || 0).toFixed(3));
                            startInhibitorShape.setAttribute('stroke', edgeDisplayStroke);
                            startInhibitorShape.setAttribute('stroke-width', Number(edge.strokeWidth || 1.6).toFixed(3));
                            startInhibitorShape.setAttribute('stroke-linecap', 'round');
                            startInhibitorShape.setAttribute('opacity', Number(edge.opacity || 0.95).toFixed(2));
                            startInhibitorShape.setAttribute('pointer-events', 'none');
                            edgeGroup.appendChild(startInhibitorShape);
                        }} else if (String(edge.startType || 'none') === 'arrow') {{
                            const startHead = getStartArrowHeadPoints(edge);
                            const startHeadShape = createSvg('path');
                            startHeadShape.setAttribute('d', 'M ' + startHead.left.x.toFixed(3) + ' ' + startHead.left.y.toFixed(3) + ' L ' + startHead.tip.x.toFixed(3) + ' ' + startHead.tip.y.toFixed(3) + ' L ' + startHead.right.x.toFixed(3) + ' ' + startHead.right.y.toFixed(3) + ' Z');
                            startHeadShape.setAttribute('fill', edgeDisplayStroke);
                            startHeadShape.setAttribute('stroke', 'none');
                            startHeadShape.setAttribute('stroke-linejoin', 'round');
                            startHeadShape.setAttribute('opacity', Number(edge.opacity || 0.95).toFixed(2));
                            startHeadShape.setAttribute('pointer-events', 'none');
                            edgeGroup.appendChild(startHeadShape);
                        }}
                        if (String(edge.endType || 'arrow') === 'inhibitor') {{
                            const bar = getInhibitorBarPoints(edge);
                            const inhibitorShape = createSvg('line');
                            inhibitorShape.setAttribute('x1', Number(bar.left.x || 0).toFixed(3));
                            inhibitorShape.setAttribute('y1', Number(bar.left.y || 0).toFixed(3));
                            inhibitorShape.setAttribute('x2', Number(bar.right.x || 0).toFixed(3));
                            inhibitorShape.setAttribute('y2', Number(bar.right.y || 0).toFixed(3));
                            inhibitorShape.setAttribute('stroke', edgeDisplayStroke);
                            inhibitorShape.setAttribute('stroke-width', Number(edge.strokeWidth || 1.6).toFixed(3));
                            inhibitorShape.setAttribute('stroke-linecap', 'round');
                            inhibitorShape.setAttribute('opacity', Number(edge.opacity || 0.95).toFixed(2));
                            inhibitorShape.setAttribute('pointer-events', 'none');
                            edgeGroup.appendChild(inhibitorShape);
                        }} else if (String(edge.endType || 'arrow') === 'arrow') {{
                            const head = getArrowHeadPoints(edge);
                            const headShape = createSvg('path');
                            headShape.setAttribute('d', 'M ' + head.left.x.toFixed(3) + ' ' + head.left.y.toFixed(3) + ' L ' + head.tip.x.toFixed(3) + ' ' + head.tip.y.toFixed(3) + ' L ' + head.right.x.toFixed(3) + ' ' + head.right.y.toFixed(3) + ' Z');
                            headShape.setAttribute('fill', edgeDisplayStroke);
                            headShape.setAttribute('stroke', 'none');
                            headShape.setAttribute('stroke-linejoin', 'round');
                            headShape.setAttribute('opacity', Number(edge.opacity || 0.95).toFixed(2));
                            headShape.setAttribute('pointer-events', 'none');
                            edgeGroup.appendChild(headShape);
                        }}
                        const edgeHitbox = createSvg('path');
                        edgeHitbox.setAttribute('d', buildFullEdgePath(edge));
                        edgeHitbox.setAttribute('fill', 'none');
                        edgeHitbox.setAttribute('stroke', 'transparent');
                        edgeHitbox.setAttribute('stroke-width', Math.max(8, Number(edge.strokeWidth || 3) + 5).toFixed(3));
                        edgeHitbox.setAttribute('stroke-linecap', 'round');
                        edgeHitbox.setAttribute('data-role', 'edge-body');
                        edgeHitbox.setAttribute('data-edge-id', edge.id);
                        edgeHitbox.setAttribute('pointer-events', 'stroke');
                        edgeGroup.appendChild(edgeHitbox);
                        if (getSelectedEdgeIds(local).includes(edge.id)) {{
                            const bends = getEdgeBendPoints(edge);
                            const getArrowHandleRadius = (baseRadius) => {{
                                const zoom = Math.max(0.4, Math.min(3, Number(local.zoom || 1)));
                                return Math.max(1.1, Math.min(14, Number(baseRadius || 4.5) / Math.pow(zoom, 1.25)));
                            }};
                            const edgeOutline = createSvg('path');
                            edgeOutline.setAttribute('d', buildFullEdgePath(edge));
                            edgeOutline.setAttribute('fill', 'none');
                            edgeOutline.setAttribute('stroke', 'rgba(15, 23, 42, 0.85)');
                            edgeOutline.setAttribute('stroke-width', (Number(edge.strokeWidth || 3) + 1.4).toFixed(3));
                            edgeOutline.setAttribute('stroke-dasharray', '8 6');
                            edgeOutline.setAttribute('pointer-events', 'none');
                            edgeGroup.appendChild(edgeOutline);
                            const makeEdgeHandle = (kind, x, y, radius, bendIndex = null) => {{
                                const handle = createSvg('circle');
                                handle.setAttribute('class', 'cst-edit-handle');
                                handle.setAttribute('data-role', kind);
                                handle.setAttribute('data-edge-id', edge.id);
                                if (bendIndex !== null && bendIndex !== undefined) {{
                                    handle.setAttribute('data-bend-index', String(bendIndex));
                                }}
                                handle.setAttribute('cx', Number(x || 0).toFixed(3));
                                handle.setAttribute('cy', Number(y || 0).toFixed(3));
                                handle.setAttribute('r', String(radius));
                                if (kind === 'edge-start' || kind === 'edge-end') {{
                                    handle.setAttribute('fill', '#dc2626');
                                    handle.setAttribute('stroke', '#ffffff');
                                    handle.setAttribute('stroke-width', '1.1');
                                    handle.setAttribute('opacity', '0.98');
                                }} else if (kind === 'edge-bend') {{
                                    handle.setAttribute('fill', isEdgeBent(edge) ? '#2563eb' : '#94a3b8');
                                    handle.setAttribute('stroke', '#ffffff');
                                    handle.setAttribute('stroke-width', '1.1');
                                    handle.setAttribute('opacity', '0.98');
                                }}
                                return handle;
                            }};
                            edgeGroup.appendChild(makeEdgeHandle('edge-start', edge.startX, edge.startY, getArrowHandleRadius(4.5)));
                            edgeGroup.appendChild(makeEdgeHandle('edge-end', edge.endX, edge.endY, getArrowHandleRadius(4.5)));
                            for (let bendIndex = 0; bendIndex < bends.length; bendIndex += 1) {{
                                const bendPoint = bends[bendIndex];
                                edgeGroup.appendChild(makeEdgeHandle('edge-bend', bendPoint.x, bendPoint.y, getArrowHandleRadius(5.0), bendIndex));
                            }}
                        }}
                            fragment.appendChild(edgeGroup);
                        }} catch (err) {{
                            try {{
                                console.error('[MapKinase CST] edge render error', edge && edge.id, err);
                            }} catch (_err) {{}}
                        }}
                    }}
                    if (local.addArrowMode && local.pendingArrowStart && local.pendingArrowPreview) {{
                        const previewGroup = createSvg('g');
                        previewGroup.setAttribute('class', 'cst-edit-edge-preview');
                        const previewPath = createSvg('path');
                        previewPath.setAttribute('d', buildEdgePath(local.pendingArrowPreview));
                        previewPath.setAttribute('fill', 'none');
                        previewPath.setAttribute('stroke', String(local.pendingArrowPreview.stroke || '#64748b'));
                        previewPath.setAttribute('stroke-width', Number(local.pendingArrowPreview.strokeWidth || 2.5).toFixed(3));
                        previewPath.setAttribute('stroke-linecap', 'round');
                        previewPath.setAttribute('stroke-linejoin', 'round');
                        previewPath.setAttribute('stroke-dasharray', '6 5');
                        previewPath.setAttribute('opacity', Number(local.pendingArrowPreview.opacity || 0.75).toFixed(2));
                        previewPath.setAttribute('pointer-events', 'none');
                        previewGroup.appendChild(previewPath);
                        const previewHead = getArrowHeadPoints(local.pendingArrowPreview);
                        const previewHeadShape = createSvg('path');
                        previewHeadShape.setAttribute('d', 'M ' + previewHead.left.x.toFixed(3) + ' ' + previewHead.left.y.toFixed(3) + ' L ' + previewHead.tip.x.toFixed(3) + ' ' + previewHead.tip.y.toFixed(3) + ' L ' + previewHead.right.x.toFixed(3) + ' ' + previewHead.right.y.toFixed(3) + ' Z');
                        previewHeadShape.setAttribute('fill', String(local.pendingArrowPreview.stroke || '#64748b'));
                        previewHeadShape.setAttribute('stroke', 'none');
                        previewHeadShape.setAttribute('stroke-linejoin', 'round');
                        previewHeadShape.setAttribute('opacity', Number(local.pendingArrowPreview.opacity || 0.75).toFixed(2));
                        previewHeadShape.setAttribute('pointer-events', 'none');
                        previewGroup.appendChild(previewHeadShape);
                        fragment.appendChild(previewGroup);
                    }}
                    for (const node of (local.editableNodes || [])) {{
                        const group = createSvg('g');
                        group.setAttribute('data-node-id', node.id);
                        group.setAttribute('class', 'cst-edit-node-group');
                        const shapeType = String(node.shapeType || 'ellipse');
                        const isRect = shapeType === 'rect' || shapeType === 'square';
                        const isRounded = shapeType === 'rounded';
                        const isLegend = shapeType === 'legend';
                        const isRectLike = isRect || isRounded || isLegend;
                        const isBracket = shapeType === 'bracket';
                        const isText = shapeType === 'text';
                        const isDrawingLike = !!node.isDrawingShape || isLegend || isBracket;
                        if (isText && local.showGroups === false) {{
                            continue;
                        }}
                        if (!isText && isDrawingLike && local.showCompounds === false) {{
                            continue;
                        }}
                        if (!isText && !isDrawingLike && local.showProtboxes === false) {{
                            continue;
                        }}
                        const bodyShape = createSvg((isBracket ? 'path' : ((isRectLike || isText) ? 'rect' : 'ellipse')));
                        bodyShape.setAttribute('class', 'cst-edit-node-ellipse ' + String(node.className || ''));
                        if (isBracket) {{
                            bodyShape.setAttribute('d', buildBracketPath(node));
                        }} else if (isRectLike || isText) {{
                            bodyShape.setAttribute('x', (Number(node.cx || 0) - Number(node.rx || 12)).toFixed(3));
                            bodyShape.setAttribute('y', (Number(node.cy || 0) - Number(node.ry || 9)).toFixed(3));
                            bodyShape.setAttribute('width', (Number(node.rx || 12) * 2).toFixed(3));
                            bodyShape.setAttribute('height', (Number(node.ry || 9) * 2).toFixed(3));
                            if (isText) {{
                                bodyShape.setAttribute('rx', '8');
                                bodyShape.setAttribute('ry', '8');
                            }} else if (isLegend) {{
                                bodyShape.setAttribute('rx', '0');
                                bodyShape.setAttribute('ry', '0');
                            }} else if (isRounded) {{
                                bodyShape.setAttribute('rx', '10');
                                bodyShape.setAttribute('ry', '10');
                            }}
                        }} else {{
                            bodyShape.setAttribute('cx', Number(node.cx || 0).toFixed(3));
                            bodyShape.setAttribute('cy', Number(node.cy || 0).toFixed(3));
                            bodyShape.setAttribute('rx', Number(node.rx || 12).toFixed(3));
                            bodyShape.setAttribute('ry', Number(node.ry || 9).toFixed(3));
                        }}
                        if (local.tabPreviewMode) {{
                            bodyShape.setAttribute('fill', 'none');
                            bodyShape.setAttribute('stroke', '#dc2626');
                            bodyShape.setAttribute('stroke-width', '1.8');
                            bodyShape.setAttribute('pointer-events', isBracket ? 'none' : 'all');
                            bodyShape.style.pointerEvents = isBracket ? 'none' : 'all';
                        }} else if (isText) {{
                            bodyShape.setAttribute('fill', String(node.fillColor || 'transparent'));
                            if (node.userCreated) {{
                                bodyShape.setAttribute('stroke', String(node.stroke || '#475569'));
                                bodyShape.setAttribute('stroke-width', Number(node.strokeWidth || 1.5).toFixed(3));
                            }} else {{
                                bodyShape.setAttribute('stroke', 'none');
                                bodyShape.setAttribute('stroke-width', '0');
                            }}
                            bodyShape.setAttribute('pointer-events', 'all');
                            bodyShape.style.pointerEvents = 'all';
                        }} else if (isLegend) {{
                            const gradId = String(node.legendOrientation || 'vertical').toLowerCase() === 'horizontal'
                                ? 'cst-legend-grad-h-{viewer_key}'
                                : 'cst-legend-grad-v-{viewer_key}';
                            bodyShape.setAttribute('fill', `url(#${{gradId}})`);
                            bodyShape.setAttribute('stroke', String(node.stroke || '#000000'));
                            bodyShape.setAttribute('stroke-width', Number(node.strokeWidth || 1.0).toFixed(3));
                            bodyShape.setAttribute('pointer-events', 'all');
                            bodyShape.style.pointerEvents = 'all';
                        }} else if (node.isComplex) {{
                            bodyShape.setAttribute('fill', '#e5e7eb');
                            bodyShape.setAttribute('stroke', '#2563eb');
                            bodyShape.setAttribute('stroke-width', local.proteinOvalMode ? '1.6' : Number(node.strokeWidth || 3.2).toFixed(3));
                            bodyShape.setAttribute('pointer-events', isBracket ? 'none' : 'all');
                            bodyShape.style.pointerEvents = isBracket ? 'none' : 'all';
                        }} else if (local.proteinOvalMode) {{
                            const activeOutline = node['outline_color_{active_idx}'] || node.outline_color_1 || null;
                            const outlineStroke = Array.isArray(activeOutline) && activeOutline.length >= 3
                                ? ('rgb(' + activeOutline[0] + ', ' + activeOutline[1] + ', ' + activeOutline[2] + ')')
                                : String(node.stroke || '#000000');
                            bodyShape.setAttribute('fill', isBracket ? 'none' : String(node.fillColor || node.stroke || 'rgb(166, 166, 166)'));
                            bodyShape.setAttribute('stroke', outlineStroke);
                            bodyShape.setAttribute('stroke-width', Math.max(0.6, Number(node.strokeWidth || local.proteinOutlineWidth || {prot_outline_width})).toFixed(3));
                            bodyShape.setAttribute('pointer-events', isBracket ? 'none' : 'all');
                            bodyShape.style.pointerEvents = isBracket ? 'none' : 'all';
                        }} else {{
                            bodyShape.setAttribute('fill', 'none');
                            bodyShape.setAttribute('stroke', String(node.stroke || 'rgb(166, 166, 166)'));
                            bodyShape.setAttribute('stroke-width', Number(node.strokeWidth || 4).toFixed(3));
                            bodyShape.setAttribute('pointer-events', isBracket ? 'none' : (isRectLike ? 'all' : 'stroke'));
                            bodyShape.style.pointerEvents = isBracket ? 'none' : (isRectLike ? 'all' : 'stroke');
                        }}
                        bodyShape.setAttribute('stroke-linecap', 'round');
                        bodyShape.setAttribute('stroke-linejoin', 'round');
                        bodyShape.setAttribute('opacity', Number(node.opacity || 1.0).toFixed(2));
                        bodyShape.setAttribute('data-node-id', node.id);
                        bodyShape.setAttribute('data-role', 'body');
                        if (Number(node.angle || 0)) {{
                            bodyShape.setAttribute('transform', 'rotate(' + Number(node.angle || 0).toFixed(3) + ' ' + Number(node.cx || 0).toFixed(3) + ' ' + Number(node.cy || 0).toFixed(3) + ')');
                        }}
                        if (!isBracket && !isText && !isLegend && !node.isComplex && !isDrawingLike && local.proteinOvalMode) {{
                            appendTemporalNodeSegments(node, group, bodyShape);
                        }}
                        applyTooltipAttrs(bodyShape, String(node.tooltip || node.title || ''), String(node.tooltipHtml || node.tooltip_html || ''));
                        group.appendChild(bodyShape);
                        if (isLegend && !local.tabPreviewMode) {{
                            const fontSize = 11;
                            const makeLegendText = (x, y, textValue, anchor = 'start') => {{
                                const text = createSvg('text');
                                text.setAttribute('x', Number(x).toFixed(3));
                                text.setAttribute('y', Number(y).toFixed(3));
                                text.setAttribute('font-size', String(fontSize));
                                text.setAttribute('font-family', 'Arial');
                                text.setAttribute('text-anchor', anchor);
                                text.setAttribute('fill', '#000000');
                                text.textContent = String(textValue);
                                group.appendChild(text);
                            }};
                            if (String(node.legendOrientation || 'vertical').toLowerCase() === 'horizontal') {{
                                const labelY = Number(node.cy || 0) + Number(node.ry || 9) + fontSize + 2;
                                makeLegendText(Number(node.cx || 0) - Number(node.rx || 12), labelY, String(legendConfig.maxNeg ?? ''), 'start');
                                makeLegendText(Number(node.cx || 0), labelY, '0', 'middle');
                                makeLegendText(Number(node.cx || 0) + Number(node.rx || 12) - 2, labelY, String(legendConfig.maxPos ?? ''), 'end');
                            }} else {{
                                const labelX = Number(node.cx || 0) + Number(node.rx || 12) + 6;
                                makeLegendText(labelX, Number(node.cy || 0) - Number(node.ry || 9) + fontSize - 2, String(legendConfig.maxPos ?? ''), 'start');
                                makeLegendText(labelX, Number(node.cy || 0) + (fontSize * 0.35), '0', 'start');
                                makeLegendText(labelX, Number(node.cy || 0) + Number(node.ry || 9) - 2, String(legendConfig.maxNeg ?? ''), 'start');
                            }}
                        }}
                        if (isBracket) {{
                            const bracketHitbox = createSvg('path');
                            bracketHitbox.setAttribute('d', buildBracketPath(node));
                            bracketHitbox.setAttribute('fill', 'none');
                            bracketHitbox.setAttribute('stroke', 'transparent');
                            bracketHitbox.setAttribute('stroke-width', Math.max(10, Number(node.strokeWidth || 3) + 8).toFixed(3));
                            bracketHitbox.setAttribute('stroke-linecap', 'round');
                            bracketHitbox.setAttribute('stroke-linejoin', 'round');
                            bracketHitbox.setAttribute('opacity', '0.01');
                            bracketHitbox.setAttribute('pointer-events', 'stroke');
                            bracketHitbox.setAttribute('data-node-id', node.id);
                            bracketHitbox.setAttribute('data-role', 'body');
                            if (Number(node.angle || 0)) {{
                                bracketHitbox.setAttribute('transform', 'rotate(' + Number(node.angle || 0).toFixed(3) + ' ' + Number(node.cx || 0).toFixed(3) + ' ' + Number(node.cy || 0).toFixed(3) + ')');
                            }}
                            group.appendChild(bracketHitbox);
                        }}
                        const nodePtms = autoPtmsByNode.get(String(node.id || '')) || [];
                        for (const ptm of (local.showPtms === false ? [] : nodePtms)) {{
                            const ptmPoint = resolveAutoPtmPlacementPoint(node, ptm);
                            const ptmRadius = Math.max(3, Number(ptm.radius || AUTO_PTM_RADIUS));
                            const ptmShape = String(ptm.shape || 'Circle').toLowerCase();
                            const temporalPtmSpec = temporalModeEnabled ? getPtmTemporalRenderSpec(ptm, ptmPoint, ptmRadius) : null;
                            const useInternalTemporalPtmText = !!(temporalPtmSpec && Number(temporalPtmSpec.segmentCount || 1) >= 3);
                            const effectivePtmShape = temporalPtmSpec ? temporalPtmSpec.renderShape : ptmShape;
                            const shapeNode = createSvg(effectivePtmShape === 'circle' ? 'circle' : 'rect');
                            shapeNode.setAttribute('class', 'cst-auto-ptm-circle' + (ptm.isMissing ? ' cst-missing-node' : ''));
                            if (temporalPtmSpec) {{
                                applyPtmShapeGeometry(shapeNode, temporalPtmSpec);
                            }} else if (ptmShape === 'circle') {{
                                shapeNode.setAttribute('cx', Number(ptmPoint.x || 0).toFixed(3));
                                shapeNode.setAttribute('cy', Number(ptmPoint.y || 0).toFixed(3));
                                shapeNode.setAttribute('r', ptmRadius.toFixed(3));
                            }} else {{
                                shapeNode.setAttribute('x', (Number(ptmPoint.x || 0) - ptmRadius).toFixed(3));
                                shapeNode.setAttribute('y', (Number(ptmPoint.y || 0) - ptmRadius).toFixed(3));
                                shapeNode.setAttribute('width', (ptmRadius * 2).toFixed(3));
                                shapeNode.setAttribute('height', (ptmRadius * 2).toFixed(3));
                                if (ptmShape === 'diamond') {{
                                    shapeNode.setAttribute('transform', 'rotate(45 ' + Number(ptmPoint.x || 0).toFixed(3) + ' ' + Number(ptmPoint.y || 0).toFixed(3) + ')');
                                }}
                            }}
                            if (!temporalPtmSpec) {{
                                const fillColor = ptm['fc_color_{active_idx}'] || ptm.fc_color_1 || [209, 213, 219];
                                const outlineColor = ptm['outline_color_{active_idx}'] || ptm.outline_color_1 || [0, 0, 0];
                                const fillRgb = Array.isArray(fillColor) ? fillColor.slice(0, 3) : [209, 213, 219];
                                const outlineRgb = Array.isArray(outlineColor) ? outlineColor.slice(0, 3) : [0, 0, 0];
                                shapeNode.setAttribute('fill', 'rgb(' + fillRgb[0] + ', ' + fillRgb[1] + ', ' + fillRgb[2] + ')');
                                shapeNode.setAttribute('stroke', 'rgb(' + outlineRgb[0] + ', ' + outlineRgb[1] + ', ' + outlineRgb[2] + ')');
                                shapeNode.setAttribute('stroke-width', Math.max(0.6, Number(ptm.outlineWidth || 1.0)).toFixed(2));
                                shapeNode.setAttribute('opacity', Number(ptm.opacity || 1.0).toFixed(2));
                            }} else {{
                                createPtmSegmentedBody(local, group, ptm, ptmPoint, ptmRadius, shapeNode);
                            }}
                            shapeNode.setAttribute('data-role', 'ptm-body');
                            shapeNode.setAttribute('data-ptm-id', String(ptm.id || ''));
                            shapeNode.setAttribute('data-node-id', String(node.id || ''));
                            applyTooltipAttrs(shapeNode, String(ptm.tooltip || ''), String(ptm.tooltip_html || ''));
                            group.appendChild(shapeNode);
                            const symbolIcon = String(ptm.symbolIcon || '').trim();
                            const symbolText = String(ptm.symbol || '').trim();
                            const siteLabel = String(ptm.siteLabel || ptm.label || '').trim();
                            if (useInternalTemporalPtmText) {{
                                const hasSymbol = local.showPtmSymbols !== false && (!!symbolIcon || !!symbolText);
                                const hasLabel = local.showPtmLabels !== false && !!siteLabel;
                                const slotLeftX = Number(temporalPtmSpec.left || 0) + (Number(temporalPtmSpec.width || 0) * 0.25);
                                const slotRightX = Number(temporalPtmSpec.left || 0) + (Number(temporalPtmSpec.width || 0) * 0.75);
                                const centerX = Number(temporalPtmSpec.centerX || 0);
                                const centerY = Number(temporalPtmSpec.centerY || 0);
                                if (hasSymbol && symbolIcon) {{
                                    const iconSize = Math.max(5.6, Math.min(Number(temporalPtmSpec.height || (ptmRadius * 2)) * 0.78, Number(temporalPtmSpec.width || (ptmRadius * 2)) * (hasLabel ? 0.36 : 0.58)));
                                    const imageX = hasLabel ? slotLeftX : centerX;
                                    const ptmSymbolImage = createSvg('image');
                                    ptmSymbolImage.setAttribute('class', 'cst-auto-ptm-symbol' + (ptm.isMissing ? ' cst-missing-node' : ''));
                                    ptmSymbolImage.setAttribute('href', symbolIcon);
                                    ptmSymbolImage.setAttribute('x', (imageX - (iconSize * 0.5)).toFixed(3));
                                    ptmSymbolImage.setAttribute('y', (centerY - (iconSize * 0.5)).toFixed(3));
                                    ptmSymbolImage.setAttribute('width', iconSize.toFixed(2));
                                    ptmSymbolImage.setAttribute('height', iconSize.toFixed(2));
                                    ptmSymbolImage.setAttribute('opacity', Number(ptm.opacity || 1.0).toFixed(2));
                                    ptmSymbolImage.setAttribute('pointer-events', 'none');
                                    group.appendChild(ptmSymbolImage);
                                }} else if (hasSymbol && symbolText) {{
                                    const symbolColor = Array.isArray(ptm.symbolColor) ? ptm.symbolColor.slice(0, 3) : [0, 0, 0];
                                    const ptmSymbol = createSvg('text');
                                    ptmSymbol.setAttribute('class', 'cst-auto-ptm-symbol' + (ptm.isMissing ? ' cst-missing-node' : ''));
                                    ptmSymbol.setAttribute('x', (hasLabel ? slotLeftX : centerX).toFixed(3));
                                    ptmSymbol.setAttribute('y', centerY.toFixed(3));
                                    ptmSymbol.setAttribute('font-size', Math.max(5.4, Math.min(Number(temporalPtmSpec.height || (ptmRadius * 2)) * 0.7, hasLabel ? (Number(temporalPtmSpec.width || (ptmRadius * 2)) * 0.2) : (Number(temporalPtmSpec.width || (ptmRadius * 2)) * 0.38), Number(ptm.symbolSize || (ptmRadius * 1.55)) * 0.82)).toFixed(2));
                                    ptmSymbol.setAttribute('font-family', String(ptm.symbolFont || 'Arial'));
                                    ptmSymbol.setAttribute('fill', 'rgb(' + symbolColor[0] + ', ' + symbolColor[1] + ', ' + symbolColor[2] + ')');
                                    ptmSymbol.setAttribute('text-anchor', 'middle');
                                    ptmSymbol.setAttribute('dominant-baseline', 'middle');
                                    ptmSymbol.setAttribute('opacity', Number(ptm.opacity || 1.0).toFixed(2));
                                    ptmSymbol.setAttribute('pointer-events', 'none');
                                    ptmSymbol.textContent = symbolText;
                                    group.appendChild(ptmSymbol);
                                }}
                                if (hasLabel) {{
                                    const labelColor = Array.isArray(ptm.labelColor) ? ptm.labelColor.slice(0, 3) : [0, 0, 0];
                                    const siteText = createSvg('text');
                                    siteText.setAttribute('class', 'cst-auto-ptm-site-label' + (ptm.isMissing ? ' cst-missing-node' : ''));
                                    siteText.setAttribute('x', ((hasSymbol ? slotRightX : centerX)).toFixed(3));
                                    siteText.setAttribute('y', centerY.toFixed(3));
                                    siteText.setAttribute('font-size', Math.max(4.9, Math.min(Number(temporalPtmSpec.height || (ptmRadius * 2)) * 0.64, hasSymbol ? (Number(temporalPtmSpec.width || (ptmRadius * 2)) * 0.18) : (Number(temporalPtmSpec.width || (ptmRadius * 2)) * 0.34), Number(ptm.labelSize || 7) * 0.82)).toFixed(2));
                                    siteText.setAttribute('font-family', String(ptm.labelFont || 'Arial'));
                                    siteText.setAttribute('fill', 'rgb(' + labelColor[0] + ', ' + labelColor[1] + ', ' + labelColor[2] + ')');
                                    siteText.setAttribute('text-anchor', 'middle');
                                    siteText.setAttribute('dominant-baseline', 'middle');
                                    siteText.setAttribute('opacity', Number(ptm.opacity || 1.0).toFixed(2));
                                    siteText.setAttribute('pointer-events', 'none');
                                    siteText.textContent = siteLabel;
                                    group.appendChild(siteText);
                                }}
                            }} else if (local.showPtmSymbols !== false && symbolIcon) {{
                                const ptmSymbolImage = createSvg('image');
                                const iconSize = Math.max(8.8, Number(ptm.symbolSize || (ptmRadius * 1.55)));
                                ptmSymbolImage.setAttribute('href', symbolIcon);
                                ptmSymbolImage.setAttribute('x', (Number(ptmPoint.x || 0) - (iconSize * 0.5)).toFixed(3));
                                ptmSymbolImage.setAttribute('y', (Number(ptmPoint.y || 0) - (iconSize * 0.5)).toFixed(3));
                                ptmSymbolImage.setAttribute('width', iconSize.toFixed(2));
                                ptmSymbolImage.setAttribute('height', iconSize.toFixed(2));
                                ptmSymbolImage.setAttribute('opacity', Number(ptm.opacity || 1.0).toFixed(2));
                                ptmSymbolImage.setAttribute('data-role', 'ptm-body');
                                ptmSymbolImage.setAttribute('data-ptm-id', String(ptm.id || ''));
                                ptmSymbolImage.setAttribute('data-node-id', String(node.id || ''));
                                applyTooltipAttrs(ptmSymbolImage, String(ptm.tooltip || ''), String(ptm.tooltip_html || ''));
                                group.appendChild(ptmSymbolImage);
                            }} else if (local.showPtmSymbols !== false && symbolText) {{
                                const ptmSymbol = createSvg('text');
                                const symbolColor = Array.isArray(ptm.symbolColor) ? ptm.symbolColor.slice(0, 3) : [0, 0, 0];
                                ptmSymbol.setAttribute('class', 'cst-auto-ptm-symbol' + (ptm.isMissing ? ' cst-missing-node' : ''));
                                ptmSymbol.setAttribute('x', Number(ptmPoint.x || 0).toFixed(3));
                                ptmSymbol.setAttribute('y', Number(ptmPoint.y || 0).toFixed(3));
                                ptmSymbol.setAttribute('font-size', Math.max(8.4, Number(ptm.symbolSize || (ptmRadius * 1.55))).toFixed(2));
                                ptmSymbol.setAttribute('font-family', String(ptm.symbolFont || 'Arial'));
                                ptmSymbol.setAttribute('fill', 'rgb(' + symbolColor[0] + ', ' + symbolColor[1] + ', ' + symbolColor[2] + ')');
                                ptmSymbol.setAttribute('opacity', Number(ptm.opacity || 1.0).toFixed(2));
                                ptmSymbol.setAttribute('data-role', 'ptm-body');
                                ptmSymbol.setAttribute('data-ptm-id', String(ptm.id || ''));
                                ptmSymbol.setAttribute('data-node-id', String(node.id || ''));
                                applyTooltipAttrs(ptmSymbol, String(ptm.tooltip || ''), String(ptm.tooltip_html || ''));
                                ptmSymbol.textContent = symbolText;
                                group.appendChild(ptmSymbol);
                            }}
                            if (!useInternalTemporalPtmText && local.showPtmLabels !== false && siteLabel) {{
                                const dx = Number(ptmPoint.x || 0) - Number(node.cx || 0);
                                const dy = Number(ptmPoint.y || 0) - Number(node.cy || 0);
                                const distance = Math.max(0.001, Math.hypot(dx, dy));
                                const ux = dx / distance;
                                const uy = dy / distance;
                                const labelDistance = ptmRadius + Math.max(1.2, Number(ptm.labelSize || 7) * 0.18);
                                const labelBaseX = Number(ptmPoint.x || 0) + (ux * labelDistance);
                                const labelBaseY = Number(ptmPoint.y || 0) + (uy * labelDistance);
                                const labelX = labelBaseX + Number(ptm.labelOffsetX || 0);
                                const labelY = labelBaseY + Number(ptm.labelOffsetY || 0);
                                const labelColor = Array.isArray(ptm.labelColor) ? ptm.labelColor.slice(0, 3) : [0, 0, 0];
                                const siteText = createSvg('text');
                                siteText.setAttribute('class', 'cst-auto-ptm-site-label' + (ptm.isMissing ? ' cst-missing-node' : ''));
                                siteText.setAttribute('x', labelX.toFixed(3));
                                siteText.setAttribute('y', labelY.toFixed(3));
                                siteText.setAttribute('font-size', Math.max(6.0, Number(ptm.labelSize || 7)).toFixed(2));
                                siteText.setAttribute('font-family', String(ptm.labelFont || 'Arial'));
                                siteText.setAttribute('fill', 'rgb(' + labelColor[0] + ', ' + labelColor[1] + ', ' + labelColor[2] + ')');
                                siteText.setAttribute('text-anchor', Math.abs(ux) < 0.24 ? 'middle' : (ux > 0 ? 'start' : 'end'));
                                siteText.setAttribute('dominant-baseline', Math.abs(uy) > 0.58 ? (uy > 0 ? 'hanging' : 'ideographic') : 'middle');
                                siteText.setAttribute('opacity', Number(ptm.opacity || 1.0).toFixed(2));
                                siteText.setAttribute('data-role', 'ptm-label');
                                siteText.setAttribute('data-ptm-id', String(ptm.id || ''));
                                siteText.setAttribute('data-node-id', String(node.id || ''));
                                applyTooltipAttrs(siteText, String(ptm.tooltip || ''), String(ptm.tooltip_html || ''));
                                siteText.textContent = siteLabel;
                                group.appendChild(siteText);
                            }}
                        }}
                        const isCompoundNode = !!node.isCompoundNode;
                        const centerLabel = node.isDrawingShape ? '' : String(node.displayLabel || (node.isCustom ? node.annotation : node.label) || '').trim();
                        const compoundLabel = isCompoundNode ? String(node.displayLabel || node.label || '').trim() : '';
                        const isEditingTextNode = isText && String(local.editingTextNodeId || '') === String(node.id || '');
                        const shouldShowCenterLabel = !node.isDrawingShape && !isBracket && !local.tabPreviewMode && !!centerLabel && (isText || local.proteinOvalMode || node.isCustom || String(node.annotation || '').trim());
                        if (shouldShowCenterLabel && !isEditingTextNode) {{
                            const text = createSvg('text');
                            text.setAttribute('class', 'cst-protein-oval-label');
                            text.setAttribute('data-node-id', String(node.id || ''));
                            text.setAttribute('data-role', isText ? 'label' : 'body');
                            text.setAttribute('pointer-events', 'all');
                            if (!isText && !node.annotationCommitted && (node.isCustom || String(node.annotation || '').trim())) {{
                                text.classList.add('cst-custom-annotation-label');
                            }}
                            text.setAttribute('x', Number(node.cx || 0).toFixed(3));
                            text.setAttribute('y', Number(node.cy || 0).toFixed(3));
                            const labelText = centerLabel;
                            const labelLength = Math.max(labelText.length, 1);
                            const fontSize = isText
                                ? Math.max(9, Math.min(28, Number(node.fontSize || Math.max(12, Number(node.ry || 12)))))
                                : Math.max(
                                    6.5,
                                    Math.min(
                                        15,
                                        Number(node.ry || 9) * 0.65,
                                        ((Number(node.rx || 12) * 1.55) / labelLength) * 2.2
                                    )
                                );
                            text.setAttribute('font-size', fontSize.toFixed(2));
                            if (isText) {{
                                text.setAttribute('font-weight', String(node.fontWeight || '600'));
                                text.setAttribute('fill', String(node.textColor || '#0f172a'));
                                text.setAttribute('font-family', String(node.fontFamily || '"Segoe UI", Arial, sans-serif'));
                                const textAlign = String(node.textAlign || 'center').toLowerCase();
                                const labelLines = getMultilineTextLines(labelText);
                                const lineHeight = fontSize * 1.14;
                                const baseY = Number(node.cy || 0) - (((Math.max(labelLines.length, 1) - 1) * lineHeight) * 0.5);
                                if (textAlign === 'right') {{
                                    text.setAttribute('text-anchor', 'end');
                                    text.setAttribute('x', (Number(node.cx || 0) + Number(node.rx || 12) - 10).toFixed(3));
                                }} else if (textAlign === 'center') {{
                                    text.setAttribute('text-anchor', 'middle');
                                    text.setAttribute('x', Number(node.cx || 0).toFixed(3));
                                }} else {{
                                    text.setAttribute('text-anchor', 'start');
                                    text.setAttribute('x', (Number(node.cx || 0) - Number(node.rx || 12) + 10).toFixed(3));
                                }}
                                text.setAttribute('y', baseY.toFixed(3));
                                while (text.firstChild) text.removeChild(text.firstChild);
                                labelLines.forEach((line, index) => {{
                                    const tspan = createSvg('tspan');
                                    tspan.setAttribute('x', text.getAttribute('x'));
                                    tspan.setAttribute('dy', index === 0 ? '0' : lineHeight.toFixed(3));
                                    tspan.setAttribute('data-node-id', String(node.id || ''));
                                    tspan.setAttribute('data-role', 'label');
                                    tspan.textContent = line || ' ';
                                    text.appendChild(tspan);
                                }});
                            }}
                            text.setAttribute('opacity', Number(node.opacity || 1.0).toFixed(2));
                            const labelAngle = isText ? Number(node.angle || 0) : (Number(node.angle || 0) + (Number(node.ry || 0) > Number(node.rx || 0) ? 90 : 0));
                            if (labelAngle) {{
                                text.setAttribute('transform', 'rotate(' + labelAngle.toFixed(3) + ' ' + Number(node.cx || 0).toFixed(3) + ' ' + Number(node.cy || 0).toFixed(3) + ')');
                            }}
                            if (!isText) text.textContent = labelText;
                            group.appendChild(text);
                        }}
                        if (isCompoundNode && !local.tabPreviewMode && compoundLabel) {{
                            const labelText = createSvg('text');
                            labelText.setAttribute('class', 'cst-protein-oval-label');
                            labelText.setAttribute('data-node-id', String(node.id || ''));
                            labelText.setAttribute('data-role', 'label');
                            labelText.setAttribute('pointer-events', 'none');
                            labelText.setAttribute('x', Number(node.cx || 0).toFixed(3));
                            labelText.setAttribute('y', (Number(node.cy || 0) + Number(node.ry || 3) + 7).toFixed(3));
                            labelText.setAttribute('font-size', '8.00');
                            labelText.setAttribute('fill', '#0f172a');
                            labelText.setAttribute('text-anchor', 'middle');
                            labelText.setAttribute('opacity', Number(node.opacity || 1.0).toFixed(2));
                            labelText.textContent = compoundLabel;
                            group.appendChild(labelText);
                        }}
                        if (isNodeSelected(local, node.id)) {{
                            bodyShape.classList.add('cst-node-selected');
                        }}
                        if (isNodeSelected(local, node.id)) {{
                            const outline = createSvg((isRectLike || isBracket || isText) ? 'rect' : 'ellipse');
                            outline.setAttribute('class', 'cst-selected-outline');
                            if (isRectLike || isBracket || isText) {{
                                const outlinePad = isText ? 1.5 : 3;
                                outline.setAttribute('x', (Number(node.cx || 0) - Number(node.rx || 12) - outlinePad).toFixed(3));
                                outline.setAttribute('y', (Number(node.cy || 0) - Number(node.ry || 9) - outlinePad).toFixed(3));
                                outline.setAttribute('width', ((Number(node.rx || 12) + outlinePad) * 2).toFixed(3));
                                outline.setAttribute('height', ((Number(node.ry || 9) + outlinePad) * 2).toFixed(3));
                                if (isRounded) {{
                                    outline.setAttribute('rx', '12');
                                    outline.setAttribute('ry', '12');
                                }}
                            }} else {{
                                outline.setAttribute('cx', Number(node.cx || 0).toFixed(3));
                                outline.setAttribute('cy', Number(node.cy || 0).toFixed(3));
                                outline.setAttribute('rx', (Number(node.rx || 12) + 3).toFixed(3));
                                outline.setAttribute('ry', (Number(node.ry || 9) + 3).toFixed(3));
                            }}
                            if (Number(node.angle || 0)) {{
                                outline.setAttribute('transform', 'rotate(' + Number(node.angle || 0).toFixed(3) + ' ' + Number(node.cx || 0).toFixed(3) + ' ' + Number(node.cy || 0).toFixed(3) + ')');
                            }}
                            group.appendChild(outline);
                        }}
                        if (local.selectedNodeId === node.id) {{
                            const xHandlePos = rotatePoint(Number(node.rx || 12), 0, Number(node.angle || 0));
                            const leftHandlePos = rotatePoint(-Number(node.rx || 12), 0, Number(node.angle || 0));
                            const topHandlePos = rotatePoint(0, -Number(node.ry || 9), Number(node.angle || 0));
                            const bottomHandlePos = rotatePoint(0, Number(node.ry || 9), Number(node.angle || 0));
                            const rotateHandlePos = rotatePoint(0, -(Number(node.ry || 9) + 18), Number(node.angle || 0));
                            const line = createSvg('line');
                            line.setAttribute('class', 'cst-handle-line');
                            line.setAttribute('x1', Number(node.cx || 0).toFixed(3));
                            line.setAttribute('y1', Number(node.cy || 0).toFixed(3));
                            line.setAttribute('x2', (Number(node.cx || 0) + rotateHandlePos.x).toFixed(3));
                            line.setAttribute('y2', (Number(node.cy || 0) + rotateHandlePos.y).toFixed(3));
                            group.appendChild(line);
                            const makeHandle = (kind, pos, radius) => {{
                                const handle = createSvg('circle');
                                handle.setAttribute('class', 'cst-edit-handle');
                                handle.setAttribute('data-node-id', node.id);
                                handle.setAttribute('data-role', kind);
                                handle.setAttribute('data-handle', kind);
                                handle.setAttribute('cx', (Number(node.cx || 0) + pos.x).toFixed(3));
                                handle.setAttribute('cy', (Number(node.cy || 0) + pos.y).toFixed(3));
                                handle.setAttribute('r', String(radius));
                                return handle;
                            }};
                            if (local.edgeResizeMode) {{
                                group.appendChild(makeHandle('left', leftHandlePos, 4.2));
                                group.appendChild(makeHandle('right', xHandlePos, 4.2));
                                group.appendChild(makeHandle('top', topHandlePos, 4.2));
                                group.appendChild(makeHandle('bottom', bottomHandlePos, 4.2));
                            }} else {{
                                group.appendChild(makeHandle('x', xHandlePos, 4.2));
                                group.appendChild(makeHandle('y', topHandlePos, 4.2));
                            }}
                            group.appendChild(makeHandle('rotate', rotateHandlePos, 5.2));
                        }}
                        if (isShapeUnderlayNode(node)) {{
                            underlayFragment.appendChild(group);
                        }} else {{
                            fragment.appendChild(group);
                        }}
                    }}
                    for (const panel of (Array.isArray(local.complexPanels) ? local.complexPanels : [])) {{
                        const node = findNodeById(local, panel.nodeId);
                        if (!node) continue;
                        const members = Array.isArray(node.complexDisplayMembers) ? node.complexDisplayMembers : [];
                        if (!members.length) continue;
                        const layout = buildComplexPanelLayout(node);
                        if (!layout) continue;
                        const panelGroup = createSvg('g');
                        panelGroup.setAttribute('class', 'cst-complex-panel');
                        const panelX = Number(panel.x || 0);
                        const panelY = Number(panel.y || 0);
                        const panelBody = createSvg('rect');
                        panelBody.setAttribute('class', 'cst-complex-panel-body');
                        panelBody.setAttribute('x', panelX.toFixed(3));
                        panelBody.setAttribute('y', panelY.toFixed(3));
                        panelBody.setAttribute('rx', '10');
                        panelBody.setAttribute('ry', '10');
                        panelBody.setAttribute('width', Number(layout.panelWidth || 140).toFixed(3));
                        panelBody.setAttribute('height', Number(layout.panelHeight || 80).toFixed(3));
                        panelGroup.appendChild(panelBody);

                        const headerHitbox = createSvg('rect');
                        headerHitbox.setAttribute('x', panelX.toFixed(3));
                        headerHitbox.setAttribute('y', panelY.toFixed(3));
                        headerHitbox.setAttribute('width', Number(layout.panelWidth || 140).toFixed(3));
                        headerHitbox.setAttribute('height', Number(layout.headerHeight || 24).toFixed(3));
                        headerHitbox.setAttribute('rx', '10');
                        headerHitbox.setAttribute('ry', '10');
                        headerHitbox.setAttribute('fill', 'transparent');
                        headerHitbox.setAttribute('stroke', 'transparent');
                        headerHitbox.setAttribute('data-role', 'complex-panel-drag');
                        headerHitbox.setAttribute('data-panel-id', String(panel.id || ''));
                        panelGroup.appendChild(headerHitbox);

                        const title = createSvg('text');
                        title.setAttribute('class', 'cst-complex-panel-title');
                        title.setAttribute('x', (panelX + layout.padding).toFixed(3));
                        title.setAttribute('y', (panelY + 16).toFixed(3));
                        title.textContent = String(node.displayLabel || node.label || 'Complex');
                        panelGroup.appendChild(title);

                        const closeRect = createSvg('rect');
                        closeRect.setAttribute('class', 'cst-complex-panel-close');
                        closeRect.setAttribute('data-role', 'complex-close');
                        closeRect.setAttribute('data-panel-id', String(panel.id || ''));
                        closeRect.setAttribute('x', (panelX + layout.panelWidth - layout.padding - layout.closeSize).toFixed(3));
                        closeRect.setAttribute('y', (panelY + 7).toFixed(3));
                        closeRect.setAttribute('rx', '4');
                        closeRect.setAttribute('ry', '4');
                        closeRect.setAttribute('width', Number(layout.closeSize).toFixed(3));
                        closeRect.setAttribute('height', Number(layout.closeSize).toFixed(3));
                        panelGroup.appendChild(closeRect);

                        const closeText = createSvg('text');
                        closeText.setAttribute('class', 'cst-complex-panel-close-text');
                        closeText.setAttribute('x', (panelX + layout.panelWidth - layout.padding - (layout.closeSize * 0.5)).toFixed(3));
                        closeText.setAttribute('y', (panelY + 14).toFixed(3));
                        closeText.textContent = '×';
                        panelGroup.appendChild(closeText);

                        const gridStartX = panelX + layout.padding;
                        const gridStartY = panelY + layout.headerHeight + layout.padding;
                        for (let memberIndex = 0; memberIndex < members.length; memberIndex += 1) {{
                            const member = members[memberIndex];
                            const memberLabel = String(member.label || member.matched_gene_symbol || member.component_label || 'Component');
                            const memberWidth = Number(layout.maxMemberWidth || estimateComplexMemberWidth(memberLabel));
                            const colIndex = memberIndex % Number(layout.cols || 1);
                            const rowIndex = Math.floor(memberIndex / Math.max(1, Number(layout.cols || 1)));
                            const memberX = gridStartX + (colIndex * (memberWidth + Number(layout.colGap || 8)));
                            const memberY = gridStartY + (rowIndex * (Number(layout.rowHeight || 15) + Number(layout.rowGap || 6)));
                            const memberRect = createSvg('rect');
                            memberRect.setAttribute('class', 'cst-complex-member-rect' + (member.has_dataset_match ? ' cst-complex-member-matched' : ''));
                            memberRect.setAttribute('x', memberX.toFixed(3));
                            memberRect.setAttribute('y', memberY.toFixed(3));
                            memberRect.setAttribute('rx', '5');
                            memberRect.setAttribute('ry', '5');
                            memberRect.setAttribute('width', Number(memberWidth).toFixed(3));
                            memberRect.setAttribute('height', Number(layout.rowHeight).toFixed(3));
                            const memberFill = Array.isArray(member.fill_color) && member.fill_color.length >= 3
                                ? member.fill_color.slice(0, 3)
                                : null;
                            memberRect.setAttribute(
                                'fill',
                                memberFill
                                    ? ('rgb(' + Number(memberFill[0] || 229) + ', ' + Number(memberFill[1] || 231) + ', ' + Number(memberFill[2] || 235) + ')')
                                    : String(member.fill_color_css || 'rgb(229, 231, 235)')
                            );
                            panelGroup.appendChild(memberRect);

                            const memberText = createSvg('text');
                            memberText.setAttribute('class', 'cst-complex-member-label');
                            memberText.setAttribute('x', (memberX + (memberWidth * 0.5)).toFixed(3));
                            memberText.setAttribute('y', (memberY + (layout.rowHeight * 0.56)).toFixed(3));
                            memberText.textContent = memberLabel;
                            panelGroup.appendChild(memberText);
                        }}
                        fragment.appendChild(panelGroup);
                    }}
                    if (local.dragState && local.dragState.mode === 'marquee') {{
                        const marquee = createSvg('rect');
                        const x = Math.min(Number(local.dragState.startPoint.overlayX || 0), Number(local.dragState.currentPoint && local.dragState.currentPoint.overlayX || 0));
                        const y = Math.min(Number(local.dragState.startPoint.overlayY || 0), Number(local.dragState.currentPoint && local.dragState.currentPoint.overlayY || 0));
                        const width = Math.abs(Number(local.dragState.startPoint.overlayX || 0) - Number(local.dragState.currentPoint && local.dragState.currentPoint.overlayX || 0));
                        const height = Math.abs(Number(local.dragState.startPoint.overlayY || 0) - Number(local.dragState.currentPoint && local.dragState.currentPoint.overlayY || 0));
                        marquee.setAttribute('x', x.toFixed(3));
                        marquee.setAttribute('y', y.toFixed(3));
                        marquee.setAttribute('width', width.toFixed(3));
                        marquee.setAttribute('height', height.toFixed(3));
                        marquee.setAttribute('fill', 'rgba(77, 163, 255, 0.12)');
                        marquee.setAttribute('stroke', '#4da3ff');
                        marquee.setAttribute('stroke-width', '1');
                        marquee.setAttribute('stroke-dasharray', '4 2');
                        marquee.setAttribute('pointer-events', 'none');
                        fragment.appendChild(marquee);
                    }}
                    overlaySvg.appendChild(underlayFragment);
                    overlaySvg.appendChild(fragment);
                    setDeleteState(local);
                    updateUndoState(local);
                    if (opacityInput) opacityInput.value = Number(local.globalOpacity || 1.0).toFixed(2);
                    if (proteinOvalButton) proteinOvalButton.classList.toggle('is-active', !!local.proteinOvalMode);
                    if (edgeResizeButton) edgeResizeButton.classList.toggle('is-active', !!local.edgeResizeMode);
                    updateInlineTextEditorLayout(local);
                    emitExternalControlState(local);
                }};
                const deleteSelectedNode = (local) => {{
                    if (local.selectedPtmId && local.selectedPtmNodeId) {{
                        if (removeNodeAutoPtmPlacement(local, local.selectedPtmNodeId, local.selectedPtmId)) {{
                            return;
                        }}
                    }}
                    const selectedEdgeIds = getSelectedEdgeIds(local);
                    if (selectedEdgeIds.length && Array.isArray(local.editableEdges)) {{
                        const beforeCount = local.editableEdges.length;
                        const snapshot = buildSnapshot(local);
                        local.editableEdges = local.editableEdges.filter((edge) => !selectedEdgeIds.includes(edge.id));
                        if (local.editableEdges.length !== beforeCount) {{
                            cleanupGroups(local);
                            pushUndoSnapshot(local, snapshot);
                            local.selectedEdgeId = null;
                            local.selectedEdgeIds = [];
                            renderEditableOverlay(local);
                        }}
                    }}
                    const selectedIds = getSelectedNodeIds(local);
                    if (!selectedIds.length) return;
                    const beforeCount = Array.isArray(local.editableNodes) ? local.editableNodes.length : 0;
                    const snapshot = buildSnapshot(local);
                    local.editableNodes = (local.editableNodes || []).filter((node) => !selectedIds.includes(node.id));
                    if (local.editableNodes.length === beforeCount) return;
                    cleanupGroups(local);
                    pushUndoSnapshot(local, snapshot);
                    setSingleSelection(local, null);
                    renderEditableOverlay(local);
                }};
                const copySelectedNode = (local) => {{
                    const selectedIds = getSelectedNodeIds(local);
                    if (!selectedIds.length) return false;
                    const nodes = (local.editableNodes || []).filter((node) => selectedIds.includes(node.id));
                    if (!nodes.length) return false;
                    local.copiedNodes = cloneStateSnapshot(nodes);
                    local.copiedEdges = [];
                    return true;
                }};
                const copySelectedEdge = (local) => {{
                    const edge = findEdgeById(local, local.selectedEdgeId);
                    if (!edge) return false;
                    local.copiedEdges = [cloneStateSnapshot(edge)];
                    local.copiedNodes = [];
                    return true;
                }};
                const pasteCopiedNode = (local) => {{
                    if (!Array.isArray(local.copiedNodes) || !local.copiedNodes.length) return false;
                    const snapshot = buildSnapshot(local);
                    const pastedNodes = cloneStateSnapshot(local.copiedNodes).map((source) => cloneEditableNode(source, {{
                        cx: Number(source.cx || 0) + 14,
                        cy: Number(source.cy || 0) + 14,
                    }}));
                    local.editableNodes = Array.isArray(local.editableNodes)
                        ? local.editableNodes.concat(pastedNodes)
                        : pastedNodes;
                    local.selectedNodeIds = pastedNodes.map((node) => node.id);
                    local.selectedNodeId = local.selectedNodeIds.length ? local.selectedNodeIds[local.selectedNodeIds.length - 1] : null;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const pasteCopiedEdge = (local) => {{
                    if (!Array.isArray(local.copiedEdges) || !local.copiedEdges.length) return false;
                    const snapshot = buildSnapshot(local);
                    const pastedEdges = cloneStateSnapshot(local.copiedEdges).map((source) => {{
                        const bends = Array.isArray(source.bendPoints) ? source.bendPoints : [];
                        return cloneEditableEdge(source, {{
                            startX: Number(source.startX || 0) + 14,
                            startY: Number(source.startY || 0) + 14,
                            endX: Number(source.endX || 0) + 14,
                            endY: Number(source.endY || 0) + 14,
                            controlX: Number(source.controlX || 0) + 14,
                            controlY: Number(source.controlY || 0) + 14,
                            bendPoints: bends.map((point) => ({{
                                x: Number(point.x || 0) + 14,
                                y: Number(point.y || 0) + 14,
                            }})),
                        }});
                    }});
                    local.editableEdges = Array.isArray(local.editableEdges)
                        ? local.editableEdges.concat(pastedEdges)
                        : pastedEdges;
                    local.selectedEdgeId = pastedEdges.length ? pastedEdges[pastedEdges.length - 1].id : null;
                    local.selectedEdgeIds = pastedEdges.map((edge) => edge.id);
                    local.selectedNodeId = null;
                    local.selectedNodeIds = [];
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const moveSelectedNodeToFront = (local) => {{
                    if (local.selectedEdgeId || !local.selectedNodeId) return false;
                    if (!local.selectedNodeId || !Array.isArray(local.editableNodes)) return false;
                    const index = local.editableNodes.findIndex((node) => node.id === local.selectedNodeId);
                    if (index < 0 || index === local.editableNodes.length - 1) return false;
                    const snapshot = buildSnapshot(local);
                    const [node] = local.editableNodes.splice(index, 1);
                    local.editableNodes.push(node);
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const moveSelectedNodeToBack = (local) => {{
                    if (local.selectedEdgeId || !local.selectedNodeId) return false;
                    if (!local.selectedNodeId || !Array.isArray(local.editableNodes)) return false;
                    const index = local.editableNodes.findIndex((node) => node.id === local.selectedNodeId);
                    if (index <= 0) return false;
                    const snapshot = buildSnapshot(local);
                    const [node] = local.editableNodes.splice(index, 1);
                    local.editableNodes.unshift(node);
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const createDefaultEdge = (startX, startY, endX, endY) => {{
                    const sx = Number(startX || 0);
                    const sy = Number(startY || 0);
                    const ex = Number(endX || 0);
                    const ey = Number(endY || 0);
                    const midX = (sx + ex) * 0.5;
                    const midY = (sy + ey) * 0.5;
                    return {{
                        id: 'cst_edge_' + Date.now() + '_' + Math.floor(Math.random() * 1000000),
                        startX: sx,
                        startY: sy,
                        endX: ex,
                        endY: ey,
                        controlX: midX,
                        controlY: midY,
                        bendPoints: [],
                        isBent: false,
                        stroke: '#0f172a',
                        strokeWidth: 1.6,
                        opacity: 0.95,
                        userCreated: true,
                        dashed: false,
                        startType: 'none',
                        endType: 'arrow',
                    }};
                }};
                const createLegendNode = (orientation = 'vertical', x = 0, y = 0) => {{
                    const isHorizontal = String(orientation || 'vertical').toLowerCase() === 'horizontal';
                    const width = isHorizontal ? 160 : 26;
                    const height = isHorizontal ? 26 : 160;
                    const node = {{
                        id: addNodeId(),
                        originalLabel: '',
                        displayLabel: '',
                        label: '',
                        matchedGene: '',
                        matchedUniprot: '',
                        foldText: '',
                        hasDatasetMatch: false,
                        cx: Number(x || 0),
                        cy: Number(y || 0),
                        rx: width * 0.5,
                        ry: height * 0.5,
                        shapeType: 'legend',
                        legendOrientation: isHorizontal ? 'horizontal' : 'vertical',
                        angle: 0,
                        strokeWidth: 1.0,
                        stroke: '#000000',
                        fillColor: 'transparent',
                        opacity: 1.0,
                        annotation: '',
                        annotationCommitted: true,
                        pendingAnnotation: '',
                        isCustom: true,
                        userCreated: true,
                        isDrawingShape: true,
                        className: 'cst-overlay-shape',
                        title: 'Legend',
                        mappingType: 'shape',
                        isComplex: false,
                        complexDisplayMembers: [],
                    }};
                }};
                const createProteinModuleNode = (entry, x = 0, y = 0, localCtx = null) => {{
                    const uniprot = String(entry && entry.uniprot || '').trim();
                    const geneSymbol = String(entry && entry.geneSymbol || '').trim() || uniprot || 'Protein';
                    const fcColor = (entry && (entry['fc_color_{active_idx}'] || entry.fc_color_1)) || [166, 166, 166];
                    const rgb = Array.isArray(fcColor) && fcColor.length >= 3 ? fcColor.slice(0, 3) : [166, 166, 166];
                    const outlineColor = (entry && (entry['outline_color_{active_idx}'] || entry.outline_color_1)) || [0, 0, 0];
                    const outlineRgb = Array.isArray(outlineColor) && outlineColor.length >= 3 ? outlineColor.slice(0, 3) : [0, 0, 0];
                    const outlineWidth = Math.max(
                        0.6,
                        Number((localCtx && localCtx.proteinOutlineWidth) || {prot_outline_width})
                    );
                    const node = {{
                        id: addNodeId(),
                        originalLabel: geneSymbol,
                        displayLabel: geneSymbol,
                        label: geneSymbol,
                        matchedGene: geneSymbol,
                        matchedUniprot: uniprot,
                        candidateUniprotIds: uniprot ? [uniprot] : [],
                        suggestedGeneSymbols: geneSymbol ? [geneSymbol] : [],
                        proteinOptions: [Object.assign({{}}, entry, {{ gene_symbol: geneSymbol, uniprot }})],
                        foldText: Number.isFinite(Number(entry && entry['fold_change_{active_idx}'])) ? Number(entry['fold_change_{active_idx}']).toFixed(3) : '',
                        hasDatasetMatch: !!uniprot,
                        cx: Number(x || 0),
                        cy: Number(y || 0),
                        rx: 18,
                        ry: 11,
                        shapeType: 'ellipse',
                        angle: 0,
                        strokeWidth: outlineWidth,
                        stroke: `rgb(${{outlineRgb[0]}}, ${{outlineRgb[1]}}, ${{outlineRgb[2]}})`,
                        fillColor: `rgb(${{rgb[0]}}, ${{rgb[1]}}, ${{rgb[2]}})`,
                        ['outline_color_{active_idx}']: outlineRgb,
                        outline_color_1: outlineRgb,
                        opacity: 1.0,
                        annotation: '',
                        annotationCommitted: true,
                        pendingAnnotation: '',
                        isCustom: true,
                        userCreated: true,
                        isDrawingShape: false,
                        className: 'cst-overlay-ellipse',
                        title: String(entry && entry.tooltip || ''),
                        tooltip: String(entry && entry.tooltip || ''),
                        tooltipHtml: String(entry && (entry.tooltip_html || entry.tooltipHtml) || ''),
                        mappingType: 'manual_dataset',
                        isComplex: false,
                        complexDisplayMembers: [],
                    }};
                    for (const key of Object.keys(entry || {{}})) {{
                        if (!/^(?:fc_color|fold_change|outline_color|outline_fold_change)_\\d+$/.test(key)) continue;
                        node[key] = entry[key];
                    }}
                    return node;
                }};
                const createMetaboliteNode = (entry, x = 0, y = 0) => {{
                    const hmdbId = String(entry && entry.hmdbId || '').trim();
                    const displayLabel = String(entry && (entry.displayLabel || entry.wikipediaId || entry.hmdbId) || hmdbId || 'Metabolite').trim();
                    const fcColor = (entry && (entry['fc_color_{active_idx}'] || entry.fc_color_1)) || [166, 166, 166];
                    const rgb = Array.isArray(fcColor) && fcColor.length >= 3 ? fcColor.slice(0, 3) : [166, 166, 166];
                    return {{
                        id: addNodeId(),
                        originalLabel: displayLabel,
                        displayLabel: displayLabel,
                        label: displayLabel,
                        matchedGene: '',
                        matchedUniprot: '',
                        foldText: Number.isFinite(Number(entry && entry['fold_change_{active_idx}'])) ? Number(entry['fold_change_{active_idx}']).toFixed(3) : '',
                        hasDatasetMatch: !!hmdbId,
                        cx: Number(x || 0),
                        cy: Number(y || 0),
                        rx: 3,
                        ry: 3,
                        shapeType: 'ellipse',
                        angle: 0,
                        strokeWidth: 1.0,
                        stroke: '#000000',
                        fillColor: `rgb(${{rgb[0]}}, ${{rgb[1]}}, ${{rgb[2]}})`,
                        opacity: 1.0,
                        annotation: '',
                        annotationCommitted: true,
                        pendingAnnotation: '',
                        isCustom: true,
                        userCreated: true,
                        isDrawingShape: true,
                        isCompoundNode: true,
                        compoundHmdbId: hmdbId,
                        compoundWikiId: String(entry && entry.wikipediaId || ''),
                        compoundKeggId: String(entry && entry.keggId || ''),
                        className: 'cst-overlay-shape',
                        title: String(entry && entry.tooltip || ''),
                        tooltip: String(entry && entry.tooltip || ''),
                        tooltipHtml: String(entry && (entry.tooltip_html || entry.tooltipHtml) || ''),
                        mappingType: 'metabolite_dataset',
                        isComplex: false,
                        complexDisplayMembers: [],
                    }};
                }};
                const getConnectableSelectedNodes = (local) => getSelectedNodeIds(local)
                    .map((nodeId) => findNodeById(local, nodeId))
                    .filter((node) => node && !node.isDrawingShape && String(node.shapeType || 'ellipse') !== 'text' && String(node.shapeType || 'ellipse') !== 'bracket');
                const getNodeAnchorPairs = (fromNode, toNode) => {{
                    const fromBounds = getNodeBounds(fromNode);
                    const toBounds = getNodeBounds(toNode);
                    if (!fromBounds || !toBounds) return [];
                    const fromCx = Number(fromNode.cx || 0);
                    const fromCy = Number(fromNode.cy || 0);
                    const toCx = Number(toNode.cx || 0);
                    const toCy = Number(toNode.cy || 0);
                    const fromAnchors = [
                        {{ x: fromBounds.left, y: fromCy }},
                        {{ x: fromBounds.right, y: fromCy }},
                        {{ x: fromCx, y: fromBounds.top }},
                        {{ x: fromCx, y: fromBounds.bottom }},
                    ];
                    const toAnchors = [
                        {{ x: toBounds.left, y: toCy }},
                        {{ x: toBounds.right, y: toCy }},
                        {{ x: toCx, y: toBounds.top }},
                        {{ x: toCx, y: toBounds.bottom }},
                    ];
                    const pairs = [];
                    fromAnchors.forEach((start) => {{
                        toAnchors.forEach((end) => {{
                            const dx = Number(end.x || 0) - Number(start.x || 0);
                            const dy = Number(end.y || 0) - Number(start.y || 0);
                            pairs.push({{
                                start,
                                end,
                                distance: Math.hypot(dx, dy),
                            }});
                        }});
                    }});
                    return pairs;
                }};
                const hasEdgeBetweenNodes = (local, fromNodeId, toNodeId) => {{
                    const a = String(fromNodeId || '');
                    const b = String(toNodeId || '');
                    return Array.isArray(local.editableEdges) && local.editableEdges.some((edge) => {{
                        const startNodeId = String(edge.startNodeId || '');
                        const endNodeId = String(edge.endNodeId || '');
                        return (startNodeId === a && endNodeId === b) || (startNodeId === b && endNodeId === a);
                    }});
                }};
                const createEdgeBetweenNodes = (local, fromNode, toNode, shortestOnly = false) => {{
                    if (!local || !fromNode || !toNode) return null;
                    if (String(fromNode.id || '') === String(toNode.id || '')) return null;
                    if (hasEdgeBetweenNodes(local, fromNode.id, toNode.id)) return null;
                    const anchorPairs = getNodeAnchorPairs(fromNode, toNode);
                    if (!anchorPairs.length) return null;
                    let bestPair = anchorPairs[0];
                    if (shortestOnly) {{
                        anchorPairs.forEach((pair) => {{
                            if (pair.distance < bestPair.distance) bestPair = pair;
                        }});
                    }} else {{
                        const preferred = anchorPairs.slice().sort((a, b) => {{
                            const aDx = Math.abs(Number(a.end.x || 0) - Number(a.start.x || 0));
                            const bDx = Math.abs(Number(b.end.x || 0) - Number(b.start.x || 0));
                            const aDy = Math.abs(Number(a.end.y || 0) - Number(a.start.y || 0));
                            const bDy = Math.abs(Number(b.end.y || 0) - Number(b.start.y || 0));
                            const aScore = a.distance + Math.min(aDx, aDy) * 0.2;
                            const bScore = b.distance + Math.min(bDx, bDy) * 0.2;
                            return aScore - bScore;
                        }});
                        bestPair = preferred[0];
                    }}
                    const edge = createDefaultEdge(bestPair.start.x, bestPair.start.y, bestPair.end.x, bestPair.end.y);
                    edge.startNodeId = String(fromNode.id || '');
                    edge.endNodeId = String(toNode.id || '');
                    return edge;
                }};
                const orderConnectNodes = (a, b) => {{
                    if (Number(a.cx || 0) < Number(b.cx || 0)) return {{ from: a, to: b }};
                    if (Number(a.cx || 0) > Number(b.cx || 0)) return {{ from: b, to: a }};
                    return Number(a.cy || 0) <= Number(b.cy || 0) ? {{ from: a, to: b }} : {{ from: b, to: a }};
                }};
                const autoConnectSelectedNodes = (local, shortestOnly = false) => {{
                    const nodes = getConnectableSelectedNodes(local);
                    if (!Array.isArray(nodes) || nodes.length < 2) return 0;
                    const snapshot = buildSnapshot(local);
                    const createdEdges = [];
                    const tryConnect = (firstNode, secondNode) => {{
                        const ordered = orderConnectNodes(firstNode, secondNode);
                        const edge = createEdgeBetweenNodes(local, ordered.from, ordered.to, shortestOnly);
                        if (!edge) return false;
                        local.editableEdges = Array.isArray(local.editableEdges) ? local.editableEdges.concat([edge]) : [edge];
                        createdEdges.push(edge.id);
                        return true;
                    }};
                    if (nodes.length === 2) {{
                        tryConnect(nodes[0], nodes[1]);
                    }} else {{
                        const rowThreshold = 30;
                        const rows = [];
                        nodes.slice().sort((a, b) => Number(a.cy || 0) - Number(b.cy || 0)).forEach((node) => {{
                            const row = rows.find((item) => Math.abs(Number(item.centerY || 0) - Number(node.cy || 0)) <= rowThreshold);
                            if (row) {{
                                row.nodes.push(node);
                                row.centerY = row.nodes.reduce((sum, item) => sum + Number(item.cy || 0), 0) / row.nodes.length;
                            }} else {{
                                rows.push({{ nodes: [node], centerY: Number(node.cy || 0) }});
                            }}
                        }});
                        rows.forEach((row) => {{
                            row.nodes.sort((a, b) => Number(a.cx || 0) - Number(b.cx || 0) || Number(a.cy || 0) - Number(b.cy || 0));
                            for (let index = 0; index < row.nodes.length - 1; index += 1) {{
                                tryConnect(row.nodes[index], row.nodes[index + 1]);
                            }}
                        }});
                    }}
                    if (!createdEdges.length) return 0;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return createdEdges.length;
                }};
                const applyEdgePreset = (edge, presetType = 'arrow') => {{
                    if (!edge) return edge;
                    const preset = String(presetType || 'arrow').toLowerCase();
                    edge.dashed = preset.startsWith('dashed_');
                    if (preset.endsWith('inhibition')) {{
                        edge.startType = 'none';
                        edge.endType = 'inhibitor';
                    }} else if (preset.endsWith('line')) {{
                        edge.startType = 'none';
                        edge.endType = 'none';
                    }} else {{
                        edge.startType = 'none';
                        edge.endType = 'arrow';
                    }}
                    return edge;
                }};
                const findEdgeById = (local, edgeId) => Array.isArray(local.editableEdges)
                    ? local.editableEdges.find((edge) => edge.id === edgeId)
                    : null;
                const getEdgeBendPoints = (edge) => {{
                    if (!edge) return [];
                    const midpoint = getActiveEdgeMidpoint(edge);
                    return midpoint ? [midpoint] : [];
                }};
                const getEdgeBounds = (edge) => {{
                    const points = [
                        {{ x: Number(edge.startX || 0), y: Number(edge.startY || 0) }},
                        {{ x: Number(edge.endX || 0), y: Number(edge.endY || 0) }},
                        ...getEdgeBendPoints(edge).map((point) => ({{
                            x: Number(point.x || 0),
                            y: Number(point.y || 0),
                        }})),
                    ];
                    const xs = points.map((point) => point.x);
                    const ys = points.map((point) => point.y);
                    const pad = Math.max(8, Number(edge.strokeWidth || 1.6) * 5);
                    return {{
                        left: Math.min(...xs) - pad,
                        right: Math.max(...xs) + pad,
                        top: Math.min(...ys) - pad,
                        bottom: Math.max(...ys) + pad,
                    }};
                }};
                const getNodeBounds = (node) => {{
                    if (!node) return null;
                    const shapeType = String(node.shapeType || 'ellipse');
                    const cx = Number(node.cx || 0);
                    const cy = Number(node.cy || 0);
                    const rx = Math.max(0, Number(node.rx || 12));
                    const ry = Math.max(0, Number(node.ry || 9));
                    if (shapeType === 'bracket') {{
                        const left = cx - rx;
                        const right = cx + rx;
                        const top = cy - ry;
                        const bottom = cy + ry;
                        const pad = Math.max(6, Number(node.strokeWidth || 3) + 4);
                        return {{
                            left: left - pad,
                            right: right + pad,
                            top: top - pad,
                            bottom: bottom + pad,
                        }};
                    }}
                    const angle = (Number(node.angle || 0) * Math.PI) / 180;
                    const cos = Math.cos(angle);
                    const sin = Math.sin(angle);
                    const xExtent = Math.abs(rx * cos) + Math.abs(ry * sin);
                    const yExtent = Math.abs(rx * sin) + Math.abs(ry * cos);
                    return {{
                        left: cx - xExtent,
                        right: cx + xExtent,
                        top: cy - yExtent,
                        bottom: cy + yExtent,
                    }};
                }};
                const setSingleEdgeSelection = (local, edgeId) => {{
                    local.selectedEdgeId = edgeId || null;
                    local.selectedEdgeIds = edgeId ? [edgeId] : [];
                    local.selectedNodeId = null;
                    local.selectedNodeIds = [];
                    clearSelectedPtm(local);
                }};
                const getQuadraticPoint = (edge, t) => {{
                    if (!isEdgeBent(edge)) {{
                        const startX = Number(edge.startX || 0);
                        const startY = Number(edge.startY || 0);
                        const endX = Number(edge.endX || 0);
                        const endY = Number(edge.endY || 0);
                        const tt = Number(t || 0);
                        return {{
                            x: startX + ((endX - startX) * tt),
                            y: startY + ((endY - startY) * tt),
                        }};
                    }}
                    const mt = 1 - Number(t || 0);
                    const bends = getEdgeBendPoints(edge);
                    const midpoint = bends.length ? bends[0] : getDefaultBendPointFromQuadratic(edge);
                    const startX = Number(edge.startX || 0);
                    const startY = Number(edge.startY || 0);
                    const endX = Number(edge.endX || 0);
                    const endY = Number(edge.endY || 0);
                    const controlX = (2 * Number(midpoint.x || 0)) - (0.5 * (startX + endX));
                    const controlY = (2 * Number(midpoint.y || 0)) - (0.5 * (startY + endY));
                    const tt = Number(t || 0);
                    return {{
                        x: (mt * mt * startX) + (2 * mt * tt * controlX) + (tt * tt * endX),
                        y: (mt * mt * startY) + (2 * mt * tt * controlY) + (tt * tt * endY),
                    }};
                }};
                const getQuadraticTangentAtEnd = (edge) => {{
                    if (!isEdgeBent(edge)) {{
                        const dx = Number(edge.endX || 0) - Number(edge.startX || 0);
                        const dy = Number(edge.endY || 0) - Number(edge.startY || 0);
                        const len = Math.hypot(dx, dy) || 1;
                        return {{ x: dx / len, y: dy / len }};
                    }}
                    const midpoint = (getEdgeBendPoints(edge)[0]) || getDefaultBendPointFromQuadratic(edge);
                    const controlX = (2 * Number(midpoint.x || 0)) - (0.5 * (Number(edge.startX || 0) + Number(edge.endX || 0)));
                    const controlY = (2 * Number(midpoint.y || 0)) - (0.5 * (Number(edge.startY || 0) + Number(edge.endY || 0)));
                    const dx = Number(edge.endX || 0) - controlX;
                    const dy = Number(edge.endY || 0) - controlY;
                    const len = Math.hypot(dx, dy) || 1;
                    return {{ x: dx / len, y: dy / len }};
                }};
                const getQuadraticTangentAtStart = (edge) => {{
                    if (!isEdgeBent(edge)) {{
                        const dx = Number(edge.endX || 0) - Number(edge.startX || 0);
                        const dy = Number(edge.endY || 0) - Number(edge.startY || 0);
                        const len = Math.hypot(dx, dy) || 1;
                        return {{ x: dx / len, y: dy / len }};
                    }}
                    const midpoint = (getEdgeBendPoints(edge)[0]) || getDefaultBendPointFromQuadratic(edge);
                    const controlX = (2 * Number(midpoint.x || 0)) - (0.5 * (Number(edge.startX || 0) + Number(edge.endX || 0)));
                    const controlY = (2 * Number(midpoint.y || 0)) - (0.5 * (Number(edge.startY || 0) + Number(edge.endY || 0)));
                    const dx = controlX - Number(edge.startX || 0);
                    const dy = controlY - Number(edge.startY || 0);
                    const len = Math.hypot(dx, dy) || 1;
                    return {{ x: dx / len, y: dy / len }};
                }};
                const getArrowHeadPoints = (edge) => {{
                    const tangent = getQuadraticTangentAtEnd(edge);
                    const headLength = Math.max(7.5, Number(edge.strokeWidth || 2.2) * 3.25);
                    const headWidth = Math.max(4.5, Number(edge.strokeWidth || 2.2) * 1.9);
                    const endX = Number(edge.endX || 0);
                    const endY = Number(edge.endY || 0);
                    const baseX = endX - (tangent.x * headLength);
                    const baseY = endY - (tangent.y * headLength);
                    const normalX = -tangent.y;
                    const normalY = tangent.x;
                    return {{
                        tip: {{ x: endX, y: endY }},
                        left: {{ x: baseX + (normalX * (headWidth * 0.5)), y: baseY + (normalY * (headWidth * 0.5)) }},
                        right: {{ x: baseX - (normalX * (headWidth * 0.5)), y: baseY - (normalY * (headWidth * 0.5)) }},
                        base: {{ x: baseX, y: baseY }},
                    }};
                }};
                const getStartArrowHeadPoints = (edge) => {{
                    const tangent = getQuadraticTangentAtStart(edge);
                    const headLength = Math.max(7.5, Number(edge.strokeWidth || 2.2) * 3.25);
                    const headWidth = Math.max(4.5, Number(edge.strokeWidth || 2.2) * 1.9);
                    const startX = Number(edge.startX || 0);
                    const startY = Number(edge.startY || 0);
                    const baseX = startX + (tangent.x * headLength);
                    const baseY = startY + (tangent.y * headLength);
                    const normalX = -tangent.y;
                    const normalY = tangent.x;
                    return {{
                        tip: {{ x: startX, y: startY }},
                        left: {{ x: baseX + (normalX * (headWidth * 0.5)), y: baseY + (normalY * (headWidth * 0.5)) }},
                        right: {{ x: baseX - (normalX * (headWidth * 0.5)), y: baseY - (normalY * (headWidth * 0.5)) }},
                        base: {{ x: baseX, y: baseY }},
                    }};
                }};
                const getInhibitorBarPoints = (edge) => {{
                    const tangent = getQuadraticTangentAtEnd(edge);
                    const barHalf = Math.max(4, Number(edge.strokeWidth || 1.6) * 2.35);
                    const normalX = -tangent.y;
                    const normalY = tangent.x;
                    const endX = Number(edge.endX || 0);
                    const endY = Number(edge.endY || 0);
                    return {{
                        left: {{ x: endX + (normalX * barHalf), y: endY + (normalY * barHalf) }},
                        right: {{ x: endX - (normalX * barHalf), y: endY - (normalY * barHalf) }},
                    }};
                }};
                const getStartInhibitorBarPoints = (edge) => {{
                    const tangent = getQuadraticTangentAtStart(edge);
                    const barHalf = Math.max(4, Number(edge.strokeWidth || 1.6) * 2.35);
                    const normalX = -tangent.y;
                    const normalY = tangent.x;
                    const startX = Number(edge.startX || 0);
                    const startY = Number(edge.startY || 0);
                    return {{
                        left: {{ x: startX + (normalX * barHalf), y: startY + (normalY * barHalf) }},
                        right: {{ x: startX - (normalX * barHalf), y: startY - (normalY * barHalf) }},
                    }};
                }};
                const buildEdgePath = (edge) => {{
                    if (!isEdgeBent(edge)) {{
                        let startX = Number(edge.startX || 0);
                        let startY = Number(edge.startY || 0);
                        let endX = Number(edge.endX || 0);
                        let endY = Number(edge.endY || 0);
                        if (String(edge.startType || 'none') === 'arrow') {{
                            const startHead = getStartArrowHeadPoints(edge);
                            startX = Number(startHead.base.x || startX);
                            startY = Number(startHead.base.y || startY);
                        }}
                        if (String(edge.endType || 'arrow') === 'arrow') {{
                            const head = getArrowHeadPoints(edge);
                            endX = Number(head.base.x || endX);
                            endY = Number(head.base.y || endY);
                        }}
                        return 'M ' + startX.toFixed(3)
                            + ' ' + startY.toFixed(3)
                            + ' L ' + endX.toFixed(3)
                            + ' ' + endY.toFixed(3);
                    }}
                    const midpoint = (getEdgeBendPoints(edge)[0]) || getDefaultBendPointFromQuadratic(edge);
                    const rawStartX = Number(edge.startX || 0);
                    const rawStartY = Number(edge.startY || 0);
                    const rawEndX = Number(edge.endX || 0);
                    const rawEndY = Number(edge.endY || 0);
                    const controlX = (2 * Number(midpoint.x || 0)) - (0.5 * (rawStartX + rawEndX));
                    const controlY = (2 * Number(midpoint.y || 0)) - (0.5 * (rawStartY + rawEndY));
                    let startX = rawStartX;
                    let startY = rawStartY;
                    let endX = rawEndX;
                    let endY = rawEndY;
                    if (String(edge.startType || 'none') === 'arrow') {{
                        const startHead = getStartArrowHeadPoints(edge);
                        startX = Number(startHead.base.x || startX);
                        startY = Number(startHead.base.y || startY);
                    }}
                    if (String(edge.endType || 'arrow') === 'arrow') {{
                        const head = getArrowHeadPoints(edge);
                        endX = Number(head.base.x || endX);
                        endY = Number(head.base.y || endY);
                    }}
                    return 'M ' + startX.toFixed(3)
                        + ' ' + startY.toFixed(3)
                        + ' Q ' + controlX.toFixed(3)
                        + ' ' + controlY.toFixed(3)
                        + ' ' + endX.toFixed(3)
                        + ' ' + endY.toFixed(3);
                }};
                const buildFullEdgePath = (edge) => {{
                    if (!isEdgeBent(edge)) {{
                        const startX = Number(edge.startX || 0);
                        const startY = Number(edge.startY || 0);
                        const endX = Number(edge.endX || 0);
                        const endY = Number(edge.endY || 0);
                        return 'M ' + startX.toFixed(3)
                            + ' ' + startY.toFixed(3)
                            + ' L ' + endX.toFixed(3)
                            + ' ' + endY.toFixed(3);
                    }}
                    const midpoint = (getEdgeBendPoints(edge)[0]) || getDefaultBendPointFromQuadratic(edge);
                    const startX = Number(edge.startX || 0);
                    const startY = Number(edge.startY || 0);
                    const endX = Number(edge.endX || 0);
                    const endY = Number(edge.endY || 0);
                    const controlX = (2 * Number(midpoint.x || 0)) - (0.5 * (startX + endX));
                    const controlY = (2 * Number(midpoint.y || 0)) - (0.5 * (startY + endY));
                    return 'M ' + startX.toFixed(3)
                        + ' ' + startY.toFixed(3)
                        + ' Q ' + controlX.toFixed(3)
                        + ' ' + controlY.toFixed(3)
                        + ' ' + endX.toFixed(3)
                        + ' ' + endY.toFixed(3);
                }};
                const AUTO_PTM_PAYLOAD = {auto_ptm_payload_json};
                const proteinSearchCatalog = {search_catalog_json};
                const metaboliteSearchCatalog = {metabolite_search_catalog_json};
                const legendConfig = {legend_config_json};
                const AUTO_PTM_RENDER_ENABLED = !!(AUTO_PTM_PAYLOAD && AUTO_PTM_PAYLOAD.enabled);
                const AUTO_PTM_RADIUS = Math.max(3, Number((AUTO_PTM_PAYLOAD && AUTO_PTM_PAYLOAD.default_radius) || (AUTO_PTM_PAYLOAD && AUTO_PTM_PAYLOAD.radius) || 4.0));
                const AUTO_PTM_ANGLE_PRIORITY = [270, 315, 225, 0, 180, 45, 135, 90, 292.5, 247.5, 337.5, 22.5, 202.5, 157.5, 67.5, 112.5, 258, 282, 348, 12];
                const AUTO_PTM_OUTWARD_FACTORS = [0.14, 0.34, 0.56, 0.78, 0.94];
                const AUTO_PTM_MIN_OUTWARD_FACTOR = 0.08;
                const AUTO_PTM_PULL_IN_STEP = 0.02;
                const normalizeAutoPtmUniprotKey = (value) => String(value || '').trim().toUpperCase().split('-')[0];
                const getAutoPtmEntryKey = (entry) => String((entry && (entry.site_label || entry.siteLabel || entry.label)) || '').trim().toUpperCase();
                const getNodeAutoPtmSourceKey = (node) => {{
                    const matchedKey = normalizeAutoPtmUniprotKey(node && node.matchedUniprot);
                    if (matchedKey) return matchedKey;
                    for (const candidate of (Array.isArray(node && node.candidateUniprotIds) ? node.candidateUniprotIds : [])) {{
                        const normalized = normalizeAutoPtmUniprotKey(candidate);
                        if (normalized) return normalized;
                    }}
                    return '';
                }};
                const getNodeForcedAutoPtmKeys = (node, uniprotKey) => {{
                    if (!node || !uniprotKey) return [];
                    const keyMap = (node && typeof node.forcedAutoPtmKeysByUniprot === 'object' && node.forcedAutoPtmKeysByUniprot) || null;
                    const forced = keyMap ? keyMap[uniprotKey] : null;
                    return Array.isArray(forced)
                        ? forced.map((item) => String(item || '').trim().toUpperCase()).filter(Boolean)
                        : [];
                }};
                const getNodeHiddenAutoPtmKeys = (node, uniprotKey) => {{
                    if (!node || !uniprotKey) return [];
                    const keyMap = (node && typeof node.hiddenAutoPtmKeysByUniprot === 'object' && node.hiddenAutoPtmKeysByUniprot) || null;
                    const hidden = keyMap ? keyMap[uniprotKey] : null;
                    return Array.isArray(hidden)
                        ? hidden.map((item) => String(item || '').trim().toUpperCase()).filter(Boolean)
                        : [];
                }};
                const addNodeForcedAutoPtmKey = (node, uniprotKey, ptmKey) => {{
                    if (!node || !uniprotKey || !ptmKey) return false;
                    if (!node.forcedAutoPtmKeysByUniprot || typeof node.forcedAutoPtmKeysByUniprot !== 'object') {{
                        node.forcedAutoPtmKeysByUniprot = {{}};
                    }}
                    const normalizedKey = String(ptmKey || '').trim().toUpperCase();
                    if (!normalizedKey) return false;
                    const current = Array.isArray(node.forcedAutoPtmKeysByUniprot[uniprotKey])
                        ? node.forcedAutoPtmKeysByUniprot[uniprotKey].map((item) => String(item || '').trim().toUpperCase()).filter(Boolean)
                        : [];
                    if (current.includes(normalizedKey)) return false;
                    current.unshift(normalizedKey);
                    node.forcedAutoPtmKeysByUniprot[uniprotKey] = current;
                    return true;
                }};
                const removeNodeForcedAutoPtmKey = (node, uniprotKey, ptmKey) => {{
                    if (!node || !uniprotKey || !ptmKey || !node.forcedAutoPtmKeysByUniprot || typeof node.forcedAutoPtmKeysByUniprot !== 'object') return false;
                    const normalizedKey = String(ptmKey || '').trim().toUpperCase();
                    const current = Array.isArray(node.forcedAutoPtmKeysByUniprot[uniprotKey])
                        ? node.forcedAutoPtmKeysByUniprot[uniprotKey].map((item) => String(item || '').trim().toUpperCase()).filter(Boolean)
                        : [];
                    const next = current.filter((item) => item !== normalizedKey);
                    if (next.length === current.length) return false;
                    if (next.length) {{
                        node.forcedAutoPtmKeysByUniprot[uniprotKey] = next;
                    }} else {{
                        delete node.forcedAutoPtmKeysByUniprot[uniprotKey];
                    }}
                    return true;
                }};
                const addNodeHiddenAutoPtmKey = (node, uniprotKey, ptmKey) => {{
                    if (!node || !uniprotKey || !ptmKey) return false;
                    if (!node.hiddenAutoPtmKeysByUniprot || typeof node.hiddenAutoPtmKeysByUniprot !== 'object') {{
                        node.hiddenAutoPtmKeysByUniprot = {{}};
                    }}
                    const normalizedKey = String(ptmKey || '').trim().toUpperCase();
                    if (!normalizedKey) return false;
                    const current = Array.isArray(node.hiddenAutoPtmKeysByUniprot[uniprotKey])
                        ? node.hiddenAutoPtmKeysByUniprot[uniprotKey].map((item) => String(item || '').trim().toUpperCase()).filter(Boolean)
                        : [];
                    if (current.includes(normalizedKey)) return false;
                    current.unshift(normalizedKey);
                    node.hiddenAutoPtmKeysByUniprot[uniprotKey] = current;
                    return true;
                }};
                const removeNodeHiddenAutoPtmKey = (node, uniprotKey, ptmKey) => {{
                    if (!node || !uniprotKey || !ptmKey || !node.hiddenAutoPtmKeysByUniprot || typeof node.hiddenAutoPtmKeysByUniprot !== 'object') return false;
                    const normalizedKey = String(ptmKey || '').trim().toUpperCase();
                    const current = Array.isArray(node.hiddenAutoPtmKeysByUniprot[uniprotKey])
                        ? node.hiddenAutoPtmKeysByUniprot[uniprotKey].map((item) => String(item || '').trim().toUpperCase()).filter(Boolean)
                        : [];
                    const next = current.filter((item) => item !== normalizedKey);
                    if (next.length === current.length) return false;
                    if (next.length) {{
                        node.hiddenAutoPtmKeysByUniprot[uniprotKey] = next;
                    }} else {{
                        delete node.hiddenAutoPtmKeysByUniprot[uniprotKey];
                    }}
                    return true;
                }};
                const autoPtmItemsByUniprot = (AUTO_PTM_PAYLOAD && AUTO_PTM_PAYLOAD.items_by_uniprot) || {{}};
                const searchProteinCatalog = (query, limit = 40) => {{
                    const rawQuery = String(query || '').trim();
                    if (!rawQuery) return {{ results: [] }};
                    let regex = null;
                    try {{
                        regex = new RegExp(rawQuery, 'i');
                    }} catch (_err) {{
                        return {{ error: 'Invalid regex', results: [] }};
                    }}
                    const results = [];
                    for (const entry of (Array.isArray(proteinSearchCatalog) ? proteinSearchCatalog : [])) {{
                        const aliasText = Array.isArray(entry.searchAliases) ? entry.searchAliases.join(' ') : '';
                        const searchText = `${{String(entry.uniprot || '')}} ${{String(entry.geneSymbol || '')}} ${{aliasText}}`;
                        if (!regex.test(searchText)) continue;
                        results.push(entry);
                        if (results.length >= Math.max(1, Number(limit || 40))) break;
                    }}
                    return {{ results }};
                }};
                const searchObjectCatalog = (query, limit = 40) => {{
                    const rawQuery = String(query || '').trim();
                    if (!rawQuery) return {{ results: [] }};
                    let regex = null;
                    try {{
                        regex = new RegExp(rawQuery, 'i');
                    }} catch (_err) {{
                        return {{ error: 'Invalid regex', results: [] }};
                    }}
                    const maxResults = Math.max(1, Number(limit || 40));
                    const proteinResults = [];
                    for (const entry of (Array.isArray(proteinSearchCatalog) ? proteinSearchCatalog : [])) {{
                        const aliasText = Array.isArray(entry.searchAliases) ? entry.searchAliases.join(' ') : '';
                        const searchText = `${{String(entry.uniprot || '')}} ${{String(entry.geneSymbol || '')}} ${{aliasText}}`;
                        if (!regex.test(searchText)) continue;
                        proteinResults.push(Object.assign({{ kind: 'protein', label: `${{String(entry.uniprot || '')}} - ${{String(entry.geneSymbol || '')}}` }}, entry));
                    }}
                    const metaboliteResults = [];
                    for (const entry of (Array.isArray(metaboliteSearchCatalog) ? metaboliteSearchCatalog : [])) {{
                        const searchText = [
                            String(entry.hmdbId || ''),
                            String(entry.wikipediaId || ''),
                            String(entry.keggId || ''),
                            String(entry.displayLabel || ''),
                            String(entry.name || ''),
                            String(entry.synonyms || '')
                        ].join(' ');
                        if (!regex.test(searchText)) continue;
                        const labelParts = [String(entry.hmdbId || '').trim(), String(entry.wikipediaId || entry.displayLabel || '').trim()].filter(Boolean);
                        metaboliteResults.push(Object.assign({{
                            kind: 'metabolite',
                            label: labelParts.join(' - ')
                        }}, entry));
                    }}
                    const results = [];
                    const primary = /hmdb|^c\\d{{5}}$|wikipedia|wiki/i.test(rawQuery) ? metaboliteResults : proteinResults;
                    const secondary = primary === metaboliteResults ? proteinResults : metaboliteResults;
                    let idx = 0;
                    while (results.length < maxResults && (idx < primary.length || idx < secondary.length)) {{
                        if (idx < primary.length && results.length < maxResults) results.push(primary[idx]);
                        if (idx < secondary.length && results.length < maxResults) results.push(secondary[idx]);
                        idx += 1;
                    }}
                    return {{ results }};
                }};
                const isAutoPtmEligibleNode = (node) => {{
                    if (!node) return false;
                    const shapeType = String(node.shapeType || 'ellipse');
                    if (shapeType === 'text' || shapeType === 'bracket') return false;
                    if (node.isDrawingShape) return false;
                    if (node.isComplex) return false;
                    return !!String(node.displayLabel || node.label || node.annotation || '').trim();
                }};
                const getNodeAutoPtmEntries = (node, options = {{}}) => {{
                    const includeHidden = !!(options && options.includeHidden);
                    const candidateKeys = [];
                    const pushKey = (value) => {{
                        const key = normalizeAutoPtmUniprotKey(value);
                        if (!key || candidateKeys.includes(key)) return;
                        candidateKeys.push(key);
                    }};
                    const matchedKey = normalizeAutoPtmUniprotKey(node && node.matchedUniprot);
                    if (matchedKey) pushKey(matchedKey);
                    for (const item of (Array.isArray(node && node.candidateUniprotIds) ? node.candidateUniprotIds : [])) {{
                        pushKey(item);
                    }}
                    for (const key of candidateKeys) {{
                        const matched = autoPtmItemsByUniprot[key];
                        if (Array.isArray(matched) && matched.length) {{
                            const forcedKeys = getNodeForcedAutoPtmKeys(node, key);
                            const hiddenKeys = includeHidden ? [] : getNodeHiddenAutoPtmKeys(node, key);
                            const visibleEntries = hiddenKeys.length
                                ? matched.filter((entry) => !hiddenKeys.includes(getAutoPtmEntryKey(entry)))
                                : matched.slice();
                            if (!forcedKeys.length) return visibleEntries;
                            const forcedEntries = [];
                            const regularEntries = [];
                            for (const entry of visibleEntries) {{
                                if (forcedKeys.includes(getAutoPtmEntryKey(entry))) {{
                                    forcedEntries.push(entry);
                                }} else {{
                                    regularEntries.push(entry);
                                }}
                            }}
                            forcedEntries.sort((a, b) => forcedKeys.indexOf(getAutoPtmEntryKey(a)) - forcedKeys.indexOf(getAutoPtmEntryKey(b)));
                            return forcedEntries.concat(regularEntries);
                        }}
                    }}
                    return [];
                }};
                const isNodeObstacleForAutoPtm = (node) => {{
                    if (!node) return false;
                    const shapeType = String(node.shapeType || 'ellipse');
                    return shapeType !== 'text' && shapeType !== 'bracket';
                }};
                const pointSegmentDistance = (px, py, ax, ay, bx, by) => {{
                    const dx = Number(bx || 0) - Number(ax || 0);
                    const dy = Number(by || 0) - Number(ay || 0);
                    const lenSq = (dx * dx) + (dy * dy);
                    if (lenSq <= 0.0001) return Math.hypot(px - Number(ax || 0), py - Number(ay || 0));
                    const t = Math.max(0, Math.min(1, (((px - Number(ax || 0)) * dx) + ((py - Number(ay || 0)) * dy)) / lenSq));
                    const projX = Number(ax || 0) + (dx * t);
                    const projY = Number(ay || 0) + (dy * t);
                    return Math.hypot(px - projX, py - projY);
                }};
                const buildAutoPtmEdgeSegments = (local) => {{
                    const segments = [];
                    const obstacleEdges = [];
                    for (const edge of (Array.isArray(local.editableEdges) ? local.editableEdges : [])) {{
                        obstacleEdges.push(edge);
                    }}
                    for (const edge of (Array.isArray(local.autoPtmObstacleEdges) ? local.autoPtmObstacleEdges : [])) {{
                        obstacleEdges.push(edge);
                    }}
                    for (const edge of obstacleEdges) {{
                        const strokePad = Math.max(3.5, Number(edge.strokeWidth || 1.6) * 0.9);
                        let prev = getQuadraticPoint(edge, 0);
                        for (let step = 1; step <= 14; step += 1) {{
                            const next = getQuadraticPoint(edge, step / 14);
                            segments.push({{
                                ax: Number(prev.x || 0),
                                ay: Number(prev.y || 0),
                                bx: Number(next.x || 0),
                                by: Number(next.y || 0),
                                pad: strokePad,
                            }});
                            prev = next;
                        }}
                        if (String(edge.startType || 'none') === 'arrow') {{
                            const head = getStartArrowHeadPoints(edge);
                            segments.push({{ ax: head.tip.x, ay: head.tip.y, bx: head.left.x, by: head.left.y, pad: 2.8 }});
                            segments.push({{ ax: head.tip.x, ay: head.tip.y, bx: head.right.x, by: head.right.y, pad: 2.8 }});
                            segments.push({{ ax: head.left.x, ay: head.left.y, bx: head.right.x, by: head.right.y, pad: 2.4 }});
                        }} else if (String(edge.startType || 'none') === 'inhibitor') {{
                            const bar = getStartInhibitorBarPoints(edge);
                            segments.push({{ ax: bar.left.x, ay: bar.left.y, bx: bar.right.x, by: bar.right.y, pad: 2.6 }});
                        }}
                        if (String(edge.endType || 'arrow') === 'arrow') {{
                            const head = getArrowHeadPoints(edge);
                            segments.push({{ ax: head.tip.x, ay: head.tip.y, bx: head.left.x, by: head.left.y, pad: 2.8 }});
                            segments.push({{ ax: head.tip.x, ay: head.tip.y, bx: head.right.x, by: head.right.y, pad: 2.8 }});
                            segments.push({{ ax: head.left.x, ay: head.left.y, bx: head.right.x, by: head.right.y, pad: 2.4 }});
                        }} else if (String(edge.endType || 'arrow') === 'inhibitor') {{
                            const bar = getInhibitorBarPoints(edge);
                            segments.push({{ ax: bar.left.x, ay: bar.left.y, bx: bar.right.x, by: bar.right.y, pad: 2.6 }});
                        }}
                    }}
                    return segments;
                }};
                const circleOverlapsAutoPtmEdge = (x, y, radius, segments) => {{
                    for (const segment of (segments || [])) {{
                        if (pointSegmentDistance(x, y, segment.ax, segment.ay, segment.bx, segment.by) <= (radius + Number(segment.pad || 0))) {{
                            return true;
                        }}
                    }}
                    return false;
                }};
                const circleOverlapsAutoPtmNode = (x, y, radius, node, ownerNodeId = null) => {{
                    if (!isNodeObstacleForAutoPtm(node)) return false;
                    if (ownerNodeId && String(node.id || '') === String(ownerNodeId || '')) return false;
                    const localPoint = rotatePoint(
                        Number(x || 0) - Number(node.cx || 0),
                        Number(y || 0) - Number(node.cy || 0),
                        -Number(node.angle || 0)
                    );
                    const rx = Math.max(2, Number(node.rx || 12));
                    const ry = Math.max(2, Number(node.ry || 9));
                    const pad = radius + 2.2;
                    if (String(node.shapeType || 'ellipse') === 'rect') {{
                        return Math.abs(Number(localPoint.x || 0)) <= (rx + pad) && Math.abs(Number(localPoint.y || 0)) <= (ry + pad);
                    }}
                    const normX = Number(localPoint.x || 0) / (rx + pad);
                    const normY = Number(localPoint.y || 0) / (ry + pad);
                    return ((normX * normX) + (normY * normY)) <= 1;
                }};
                const isValidAutoPtmCandidate = (candidate, radius, node, nodeObstacles, edgeSegments, placed) => {{
                    if (
                        Number(candidate.x || 0) < (radius + 3) ||
                        Number(candidate.y || 0) < (radius + 3) ||
                        Number(candidate.x || 0) > (pageWidth - radius - 3) ||
                        Number(candidate.y || 0) > (pageHeight - radius - 3)
                    ) {{
                        return false;
                    }}
                    if (circleOverlapsAutoPtmEdge(candidate.x, candidate.y, radius, edgeSegments)) return false;
                    for (const otherNode of (nodeObstacles || [])) {{
                        if (circleOverlapsAutoPtmNode(candidate.x, candidate.y, radius, otherNode, node.id)) {{
                            return false;
                        }}
                    }}
                    for (const otherPtm of (placed || [])) {{
                        if (Math.hypot(Number(candidate.x || 0) - Number(otherPtm.x || 0), Number(candidate.y || 0) - Number(otherPtm.y || 0)) < (radius + Number(otherPtm.radius || radius) + 1.8)) {{
                            return false;
                        }}
                    }}
                    return true;
                }};
                const getAutoPtmCandidate = (node, angleDeg, outwardFactor, radius) => {{
                    const radians = (Number(angleDeg || 0) * Math.PI) / 180;
                    const dirWorld = {{ x: Math.cos(radians), y: Math.sin(radians) }};
                    const dirLocal = rotatePoint(dirWorld.x, dirWorld.y, -Number(node.angle || 0));
                    const rx = Math.max(2, Number(node.rx || 12));
                    const ry = Math.max(2, Number(node.ry || 9));
                    let boundaryDistance = 0;
                    if (String(node.shapeType || 'ellipse') === 'rect') {{
                        const dx = Math.abs(Number(dirLocal.x || 0));
                        const dy = Math.abs(Number(dirLocal.y || 0));
                        const distX = dx > 0.0001 ? (rx / dx) : Infinity;
                        const distY = dy > 0.0001 ? (ry / dy) : Infinity;
                        boundaryDistance = Math.min(distX, distY);
                    }} else {{
                        const denom = Math.sqrt(
                            ((Number(dirLocal.x || 0) * Number(dirLocal.x || 0)) / Math.max(rx * rx, 1)) +
                            ((Number(dirLocal.y || 0) * Number(dirLocal.y || 0)) / Math.max(ry * ry, 1))
                        );
                        boundaryDistance = denom > 0.0001 ? (1 / denom) : Math.max(rx, ry);
                    }}
                    const boundaryLocal = {{
                        x: Number(dirLocal.x || 0) * boundaryDistance,
                        y: Number(dirLocal.y || 0) * boundaryDistance,
                    }};
                    const boundaryWorld = rotatePoint(boundaryLocal.x, boundaryLocal.y, Number(node.angle || 0));
                    return {{
                        x: Number(node.cx || 0) + Number(boundaryWorld.x || 0) + (Number(dirWorld.x || 0) * radius * Number(outwardFactor || 0)),
                        y: Number(node.cy || 0) + Number(boundaryWorld.y || 0) + (Number(dirWorld.y || 0) * radius * Number(outwardFactor || 0)),
                    }};
                }};
                const resolveAutoPtmPlacementPoint = (node, placement) => {{
                    return getAutoPtmCandidate(
                        node,
                        Number(placement && placement.angleDeg || 0),
                        Number(placement && placement.outwardFactor || 0.34),
                        Number(placement && placement.radius || AUTO_PTM_RADIUS),
                    );
                }};
                const buildAutoPtmPlacementForEntry = (local, node, ptmEntry, order, placed, nodeObstacles = null, edgeSegments = null, allowFallback = false) => {{
                    const ptmRadius = Math.max(3, Number(ptmEntry && ptmEntry.radius || AUTO_PTM_RADIUS));
                    const resolvedNodeObstacles = Array.isArray(nodeObstacles)
                        ? nodeObstacles
                        : (Array.isArray(local.editableNodes) ? local.editableNodes : []).filter((item) => isNodeObstacleForAutoPtm(item));
                    const resolvedEdgeSegments = Array.isArray(edgeSegments)
                        ? edgeSegments
                        : buildAutoPtmEdgeSegments(local);
                    let fallbackAccepted = null;
                    for (const angleDeg of AUTO_PTM_ANGLE_PRIORITY) {{
                        let accepted = null;
                        for (const outwardFactor of AUTO_PTM_OUTWARD_FACTORS) {{
                            const candidate = getAutoPtmCandidate(node, angleDeg, outwardFactor, ptmRadius);
                            if (!fallbackAccepted) {{
                                fallbackAccepted = {{
                                    x: Number(candidate.x || 0),
                                    y: Number(candidate.y || 0),
                                    angleDeg: Number(angleDeg || 0),
                                    outwardFactor: Number(outwardFactor || 0),
                                }};
                            }}
                            if (!isValidAutoPtmCandidate(candidate, ptmRadius, node, resolvedNodeObstacles, resolvedEdgeSegments, placed)) continue;
                            accepted = {{
                                x: Number(candidate.x || 0),
                                y: Number(candidate.y || 0),
                                outwardFactor: Number(outwardFactor || 0),
                            }};
                            break;
                        }}
                        if (!accepted) continue;
                        for (let inwardFactor = Number(accepted.outwardFactor || 0) - AUTO_PTM_PULL_IN_STEP; inwardFactor >= AUTO_PTM_MIN_OUTWARD_FACTOR; inwardFactor -= AUTO_PTM_PULL_IN_STEP) {{
                            const inwardCandidate = getAutoPtmCandidate(node, angleDeg, inwardFactor, ptmRadius);
                            if (!isValidAutoPtmCandidate(inwardCandidate, ptmRadius, node, resolvedNodeObstacles, resolvedEdgeSegments, placed)) break;
                            accepted = {{
                                x: Number(inwardCandidate.x || 0),
                                y: Number(inwardCandidate.y || 0),
                                outwardFactor: Number(inwardFactor || 0),
                            }};
                        }}
                        const placement = {{
                            id: 'cst_auto_ptm_' + String(node.id || '') + '_' + String(order),
                            nodeId: String(node.id || ''),
                            radius: ptmRadius,
                            order,
                            angleDeg: Number(angleDeg || 0),
                            outwardFactor: Number(accepted.outwardFactor || 0),
                            isMissing: String(node.className || '').includes('cst-missing-node'),
                            opacity: Math.max(0.18, Math.min(1.0, Number(node.opacity || local.globalOpacity || 1.0))),
                            siteLabel: String((ptmEntry && (ptmEntry.site_label || ptmEntry.label)) || order),
                            symbol: String((ptmEntry && ptmEntry.symbol) || ''),
                            shape: String((ptmEntry && ptmEntry.shape) || 'Circle'),
                            labelFont: String((ptmEntry && ptmEntry.label_font) || 'Arial'),
                            labelColor: Array.isArray(ptmEntry && ptmEntry.label_color) ? ptmEntry.label_color.slice(0, 3) : [0, 0, 0],
                            labelSize: Math.max(6, Number((ptmEntry && ptmEntry.label_size) || 7)),
                            symbolFont: String((ptmEntry && ptmEntry.symbol_font) || 'Arial'),
                            symbolColor: Array.isArray(ptmEntry && ptmEntry.symbol_color) ? ptmEntry.symbol_color.slice(0, 3) : [0, 0, 0],
                            symbolIcon: String((ptmEntry && ptmEntry.symbol_icon) || ''),
                            symbolSize: Math.max(8.4, Number((ptmEntry && ptmEntry.symbol_size) || (ptmRadius * 1.55))),
                            outlineWidth: Math.max(0.6, Number((ptmEntry && ptmEntry.outline_width) || 1.0)),
                            tooltip: String((ptmEntry && ptmEntry.tooltip) || ''),
                            tooltip_html: String((ptmEntry && ptmEntry.tooltip_html) || ''),
                            fc_color_1: (ptmEntry && ptmEntry.fc_color_1) || [166, 166, 166],
                            outline_color_1: (ptmEntry && ptmEntry.outline_color_1) || [0, 0, 0],
                            fold_change_1: (ptmEntry && ptmEntry.fold_change_1),
                        }};
                        for (const key of Object.keys(ptmEntry || {{}})) {{
                            if (!/^fc_color_\\d+$/.test(key) && !/^outline_color_\\d+$/.test(key) && !/^fold_change_\\d+$/.test(key) && !/^outline_fold_change_\\d+$/.test(key)) continue;
                            placement[key] = ptmEntry[key];
                        }}
                        return {{
                            placement,
                            placedPoint: {{
                                x: Number(accepted.x || 0),
                                y: Number(accepted.y || 0),
                                radius: ptmRadius,
                            }},
                        }};
                    }}
                    if (allowFallback && fallbackAccepted) {{
                        const placement = {{
                            id: 'cst_auto_ptm_' + String(node.id || '') + '_' + String(order),
                            nodeId: String(node.id || ''),
                            radius: ptmRadius,
                            order,
                            angleDeg: Number(fallbackAccepted.angleDeg || 0),
                            outwardFactor: Number(fallbackAccepted.outwardFactor || 0),
                            isMissing: String(node.className || '').includes('cst-missing-node'),
                            opacity: Math.max(0.18, Math.min(1.0, Number(node.opacity || local.globalOpacity || 1.0))),
                            siteLabel: String((ptmEntry && (ptmEntry.site_label || ptmEntry.label)) || order),
                            symbol: String((ptmEntry && ptmEntry.symbol) || ''),
                            shape: String((ptmEntry && ptmEntry.shape) || 'Circle'),
                            labelFont: String((ptmEntry && ptmEntry.label_font) || 'Arial'),
                            labelColor: Array.isArray(ptmEntry && ptmEntry.label_color) ? ptmEntry.label_color.slice(0, 3) : [0, 0, 0],
                            labelSize: Math.max(6, Number((ptmEntry && ptmEntry.label_size) || 7)),
                            symbolFont: String((ptmEntry && ptmEntry.symbol_font) || 'Arial'),
                            symbolColor: Array.isArray(ptmEntry && ptmEntry.symbol_color) ? ptmEntry.symbol_color.slice(0, 3) : [0, 0, 0],
                            symbolIcon: String((ptmEntry && ptmEntry.symbol_icon) || ''),
                            symbolSize: Math.max(8.4, Number((ptmEntry && ptmEntry.symbol_size) || (ptmRadius * 1.55))),
                            outlineWidth: Math.max(0.6, Number((ptmEntry && ptmEntry.outline_width) || 1.0)),
                            tooltip: String((ptmEntry && ptmEntry.tooltip) || ''),
                            tooltip_html: String((ptmEntry && ptmEntry.tooltip_html) || ''),
                            fc_color_1: (ptmEntry && ptmEntry.fc_color_1) || [166, 166, 166],
                            outline_color_1: (ptmEntry && ptmEntry.outline_color_1) || [0, 0, 0],
                            fold_change_1: (ptmEntry && ptmEntry.fold_change_1),
                        }};
                        for (const key of Object.keys(ptmEntry || {{}})) {{
                            if (!/^fc_color_\\d+$/.test(key) && !/^outline_color_\\d+$/.test(key) && !/^fold_change_\\d+$/.test(key) && !/^outline_fold_change_\\d+$/.test(key)) continue;
                            placement[key] = ptmEntry[key];
                        }}
                        return {{
                            placement,
                            placedPoint: {{
                                x: Number(fallbackAccepted.x || 0),
                                y: Number(fallbackAccepted.y || 0),
                                radius: ptmRadius,
                            }},
                        }};
                    }}
                    return null;
                }};
                const buildAutoPtmPlacements = (local) => {{
                    const placements = [];
                    const placed = [];
                    const nodeObstacles = (Array.isArray(local.editableNodes) ? local.editableNodes : []).filter((node) => isNodeObstacleForAutoPtm(node));
                    const edgeSegments = buildAutoPtmEdgeSegments(local);
                    for (const node of (Array.isArray(local.editableNodes) ? local.editableNodes : [])) {{
                        if (!isAutoPtmEligibleNode(node)) continue;
                        const nodePtms = getNodeAutoPtmEntries(node);
                        if (!Array.isArray(nodePtms) || !nodePtms.length) continue;
                        let order = 1;
                        for (const ptmEntry of nodePtms) {{
                            const built = buildAutoPtmPlacementForEntry(local, node, ptmEntry, order, placed, nodeObstacles, edgeSegments, false);
                            if (built && built.placement) {{
                                placements.push(built.placement);
                                placed.push(built.placedPoint);
                            }}
                            order += 1;
                        }}
                    }}
                    return placements;
                }};
                const setBatchSizePanelOpen = (local, isOpen) => {{
                    local.batchSizePanelOpen = !!isOpen;
                    if (batchSizeButton) batchSizeButton.classList.toggle('is-active', !!isOpen && getSelectedNodeIds(local).length > 1);
                    if (batchSizePanel) batchSizePanel.classList.toggle('is-open', !!isOpen && getSelectedNodeIds(local).length > 1);
                }};
                const applyBatchSize = (local) => {{
                    const selectedIds = getSelectedNodeIds(local);
                    if (selectedIds.length < 2 || !batchSizeValue || !batchSizeMode) return false;
                    const units = Number(batchSizeValue.value || 0);
                    if (!Number.isFinite(units) || units <= 0) return false;
                    const mode = String(batchSizeMode.value || 'both');
                    const snapshot = buildSnapshot(local);
                    let changed = false;
                    for (const node of (local.editableNodes || [])) {{
                        if (!selectedIds.includes(node.id)) continue;
                        if (mode === 'both' || mode === 'width') {{
                            node.rx = Math.max(4, units * 0.5);
                            changed = true;
                        }}
                        if (mode === 'both' || mode === 'height') {{
                            node.ry = Math.max(4, units * 0.5);
                            changed = true;
                        }}
                    }}
                    if (!changed) return false;
                    pushUndoSnapshot(local, snapshot);
                    setBatchSizePanelOpen(local, false);
                    renderEditableOverlay(local);
                    return true;
                }};
                const convertSelectedNodesToRect = (local) => {{
                    const selectedIds = getSelectedNodeIds(local);
                    if (!selectedIds.length) return false;
                    const snapshot = buildSnapshot(local);
                    let changed = false;
                    for (const node of (local.editableNodes || [])) {{
                        if (!selectedIds.includes(node.id)) continue;
                        if (String(node.shapeType || 'ellipse') === 'rect' || String(node.shapeType || 'ellipse') === 'bracket') continue;
                        node.shapeType = 'rect';
                        changed = true;
                    }}
                    if (!changed) return false;
                    pushUndoSnapshot(local, snapshot);
                    renderEditableOverlay(local);
                    return true;
                }};
                const buildBracketPath = (node) => {{
                    const cx = Number(node.cx || 0);
                    const cy = Number(node.cy || 0);
                    const rx = Math.max(4, Number(node.rx || 12));
                    const ry = Math.max(4, Number(node.ry || 9));
                    const leftX = cx - rx;
                    const rightX = cx + rx;
                    const topY = cy - ry;
                    const bottomY = cy + ry;
                    return 'M ' + leftX.toFixed(3)
                        + ' ' + topY.toFixed(3)
                        + ' L ' + rightX.toFixed(3)
                        + ' ' + topY.toFixed(3)
                        + ' L ' + rightX.toFixed(3)
                        + ' ' + bottomY.toFixed(3)
                        + ' L ' + leftX.toFixed(3)
                        + ' ' + bottomY.toFixed(3);
                }};
                const applyGlobalOpacity = (local, value) => {{
                    const opacity = Math.max(0.1, Math.min(1, Number(value || 1.0)));
                    if (Math.abs(opacity - Number(local.globalOpacity || 1.0)) < 0.001) return;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    local.globalOpacity = opacity;
                    for (const node of (local.editableNodes || [])) {{
                        node.opacity = opacity;
                    }}
                    for (const edge of (local.editableEdges || [])) {{
                        edge.opacity = opacity;
                    }}
                    renderEditableOverlay(local);
                }};
                const toggleSelectedEdgeDashed = (local) => {{
                    const edge = findEdgeById(local, local.selectedEdgeId);
                    if (!edge) return false;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    edge.dashed = !edge.dashed;
                    renderEditableOverlay(local);
                    return true;
                }};
                const toggleSelectedEdgeBothArrows = (local) => {{
                    const edge = findEdgeById(local, local.selectedEdgeId);
                    if (!edge) return false;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    const startType = String(edge.startType || 'none');
                    const endType = String(edge.endType || 'arrow');
                    if (endType === 'inhibitor' || startType === 'inhibitor') {{
                        if (startType === 'inhibitor' && endType === 'inhibitor') {{
                            edge.startType = 'none';
                            edge.endType = 'inhibitor';
                        }} else {{
                            edge.startType = 'inhibitor';
                            edge.endType = 'inhibitor';
                        }}
                    }} else if (startType === 'arrow' && endType === 'arrow') {{
                        edge.startType = 'none';
                        edge.endType = 'arrow';
                    }} else {{
                        edge.startType = 'arrow';
                        edge.endType = 'arrow';
                    }}
                    renderEditableOverlay(local);
                    return true;
                }};
                const toggleSelectedEdgeLine = (local) => {{
                    const edge = findEdgeById(local, local.selectedEdgeId);
                    if (!edge) return false;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    if (String(edge.startType || 'none') === 'none' && String(edge.endType || 'arrow') === 'none') {{
                        edge.startType = 'none';
                        edge.endType = 'arrow';
                    }} else {{
                        edge.startType = 'none';
                        edge.endType = 'none';
                    }}
                    renderEditableOverlay(local);
                    return true;
                }};
                const resetSelectedEdgeBend = (local, edgeId) => {{
                    const edge = findEdgeById(local, edgeId);
                    if (!edge) return false;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    edge.isBent = false;
                    edge.bendPoints = [];
                    const midpoint = getStraightEdgeMidpoint(edge);
                    edge.controlX = Number(midpoint.x || 0);
                    edge.controlY = Number(midpoint.y || 0);
                    renderEditableOverlay(local);
                    return true;
                }};
                const toggleSelectedEdgeEndType = (local) => {{
                    const edge = findEdgeById(local, local.selectedEdgeId);
                    if (!edge) return false;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    edge.endType = String(edge.endType || 'arrow') === 'inhibitor' ? 'arrow' : 'inhibitor';
                    renderEditableOverlay(local);
                    return true;
                }};
                const cycleSelectedEdgeInteractorType = (local) => {{
                    const edge = findEdgeById(local, local.selectedEdgeId);
                    if (!edge) return false;
                    const startType = String(edge.startType || 'none');
                    const endType = String(edge.endType || 'arrow');
                    const hasStart = startType !== 'none';
                    const startIsArrow = startType === 'arrow';
                    const startIsInhibitor = startType === 'inhibitor';
                    const endIsArrow = endType === 'arrow';
                    const endIsInhibitor = endType === 'inhibitor';
                    pushUndoSnapshot(local, buildSnapshot(local));
                    if ((hasStart ? startIsArrow : true) && endIsArrow) {{
                        edge.startType = hasStart ? 'inhibitor' : 'none';
                        edge.endType = 'inhibitor';
                    }} else if ((hasStart ? startIsInhibitor : true) && endIsInhibitor) {{
                        edge.startType = 'none';
                        edge.endType = 'none';
                    }} else {{
                        edge.startType = 'none';
                        edge.endType = 'arrow';
                    }}
                    renderEditableOverlay(local);
                    return true;
                }};
                const flipSelectedEdgeDirection = (local) => {{
                    const edge = findEdgeById(local, local.selectedEdgeId);
                    if (!edge) return false;
                    pushUndoSnapshot(local, buildSnapshot(local));
                    const startX = Number(edge.startX || 0);
                    const startY = Number(edge.startY || 0);
                    edge.startX = Number(edge.endX || 0);
                    edge.startY = Number(edge.endY || 0);
                    edge.endX = startX;
                    edge.endY = startY;
                    const startType = String(edge.startType || 'none');
                    edge.startType = String(edge.endType || 'arrow');
                    edge.endType = startType;
                    if (Array.isArray(edge.bendPoints) && edge.bendPoints.length) {{
                        edge.bendPoints = edge.bendPoints.slice().reverse().map((point) => ({{
                            x: Number(point.x || 0),
                            y: Number(point.y || 0),
                        }}));
                    }}
                    renderEditableOverlay(local);
                    return true;
                }};
                const setMissingState = (showMissing) => {{
                    stage.classList.toggle('cst-show-missing', !!showMissing);
                    if (missingButton) {{
                        missingButton.textContent = showMissing ? 'Hide Missing Proteins' : 'Show Missing Proteins';
                    }}
                }};
                const buildOverlayVisibilityState = (local) => {{
                    return {{
                        background: !!local.showPdfBackground,
                        protboxes: local.showProtboxes !== false,
                        ptms: local.showPtms !== false,
                        ptmLabels: local.showPtmLabels !== false,
                        ptmSymbols: local.showPtmSymbols !== false,
                        arrows: local.showArrows !== false,
                        compounds: local.showCompounds !== false,
                        groups: local.showGroups !== false,
                        tooltips: local.showTooltips !== false,
                    }};
                }};
                const applyOverlayVisibilityState = (local, key, value) => {{
                    const enabled = value !== false;
                    if (key === 'background') {{
                        local.showPdfBackground = enabled;
                        if (canvas) {{
                            canvas.style.opacity = enabled ? '1' : '0';
                        }}
                        if (fallback) {{
                            fallback.style.opacity = enabled ? '1' : '0';
                        }}
                        queueRender();
                    }} else if (key === 'protboxes') {{
                        local.showProtboxes = enabled;
                    }} else if (key === 'ptms') {{
                        local.showPtms = enabled;
                    }} else if (key === 'ptmLabels') {{
                        local.showPtmLabels = enabled;
                    }} else if (key === 'ptmSymbols') {{
                        local.showPtmSymbols = enabled;
                    }} else if (key === 'arrows') {{
                        local.showArrows = enabled;
                    }} else if (key === 'compounds') {{
                        local.showCompounds = enabled;
                    }} else if (key === 'groups') {{
                        local.showGroups = enabled;
                    }} else if (key === 'tooltips') {{
                        local.showTooltips = enabled;
                        if (!enabled) hideHoverTooltip();
                    }}
                    renderEditableOverlay(local);
                    return buildOverlayVisibilityState(local);
                }};
                const escapeExportHtml = (value) => {{
                    return String(value == null ? '' : value)
                        .replace(/&/g, '&amp;')
                        .replace(/</g, '&lt;')
                        .replace(/>/g, '&gt;')
                        .replace(/"/g, '&quot;')
                        .replace(/'/g, '&#39;');
                }};
                const isExportTextNode = (node) => String(node && node.shapeType || 'ellipse') === 'text';
                const isExportDrawingNode = (node) => !!(node && node.isDrawingShape);
                const isExportProteinNode = (node) => {{
                    if (!node) return false;
                    const shapeType = String(node.shapeType || 'ellipse');
                    if (shapeType === 'text' || shapeType === 'legend' || shapeType === 'bracket') return false;
                    return !node.isDrawingShape;
                }};
                const nodeRectForExport = (node) => {{
                    const cx = Number(node && node.cx || 0);
                    const cy = Number(node && node.cy || 0);
                    const rx = Math.max(1, Number(node && node.rx || 1));
                    const ry = Math.max(1, Number(node && node.ry || 1));
                    return {{
                        x: cx - rx,
                        y: cy - ry,
                        width: rx * 2,
                        height: ry * 2,
                        cx,
                        cy,
                    }};
                }};
                const exportArrowSideForPoint = (node, x, y) => {{
                    const rect = nodeRectForExport(node);
                    const dx = Number(x || 0) - rect.cx;
                    const dy = Number(y || 0) - rect.cy;
                    if (Math.abs(dx) > Math.abs(dy)) {{
                        return dx >= 0 ? 'East' : 'West';
                    }}
                    return dy >= 0 ? 'South' : 'North';
                }};
                const closestExportNodeForPoint = (nodes, x, y) => {{
                    let best = null;
                    let bestDist = Infinity;
                    for (const node of (nodes || [])) {{
                        if (!node || !String(node.id || '').trim()) continue;
                        const rect = nodeRectForExport(node);
                        const dx = rect.cx - Number(x || 0);
                        const dy = rect.cy - Number(y || 0);
                        const dist = Math.hypot(dx, dy);
                        if (dist < bestDist) {{
                            bestDist = dist;
                            best = node;
                        }}
                    }}
                    return best;
                }};
                const buildCustomExportSnapshot = (local) => {{
                    const settings = {{
                        pathway_id: {_json_for_inline_script(str(info.get("id") or ""))},
                        pathway_source: 'cst',
                        pathway_name: {_json_for_inline_script(str(info.get("name") or ""))},
                    }};
                    const snapshot = {{
                        general_data: {{ settings }},
                        protbox_data: [],
                        compound_data: [],
                        text_data: [],
                        arrows: [],
                        groups: [],
                    }};
                    const exportableNodes = [];
                    for (const node of (local.editableNodes || [])) {{
                        if (!node || !String(node.id || '').trim()) continue;
                        const rect = nodeRectForExport(node);
                        if (isExportTextNode(node)) {{
                            snapshot.text_data.push({{
                                text_id: String(node.id || ''),
                                x: rect.x,
                                y: rect.y,
                                width: rect.width,
                                height: rect.height,
                                label: String(node.displayLabel || node.label || 'Text'),
                                html: escapeExportHtml(String(node.displayLabel || node.label || 'Text')).replace(/\\n/g, '<br/>'),
                                text_style: {{
                                    fontFamily: String(node.fontFamily || '"Segoe UI", Arial, sans-serif'),
                                    fontSize: Number(node.fontSize || 14),
                                    fontWeight: String(node.fontWeight || '600'),
                                    color: String(node.textColor || '#0f172a'),
                                    align: String(node.textAlign || 'center'),
                                }},
                                bgcolor: String(node.fillColor || 'transparent'),
                                fgcolor: String(node.textColor || '#0f172a'),
                                border_color: String(node.stroke || '#475569'),
                                angle: Number(node.angle || 0),
                                opacity: Number(node.opacity || 1.0),
                            }});
                            exportableNodes.push(node);
                            continue;
                        }}
                        if (isExportProteinNode(node)) {{
                            snapshot.protbox_data.push({{
                                protbox_id: String(node.id || ''),
                                x: rect.x,
                                y: rect.y,
                                width: rect.width,
                                height: rect.height,
                                label: String(node.displayLabel || node.label || ''),
                                proteins: Array.isArray(node.candidateUniprotIds) && node.candidateUniprotIds.length
                                    ? cloneStateSnapshot(node.candidateUniprotIds)
                                    : (String(node.matchedUniprot || '').trim() ? [String(node.matchedUniprot || '').trim()] : []),
                                uniprot_ids: Array.isArray(node.candidateUniprotIds) && node.candidateUniprotIds.length
                                    ? cloneStateSnapshot(node.candidateUniprotIds)
                                    : (String(node.matchedUniprot || '').trim() ? [String(node.matchedUniprot || '').trim()] : []),
                                matched_uniprot: String(node.matchedUniprot || ''),
                                matched_gene: String(node.matchedGene || ''),
                                shape_type: String(node.shapeType || 'ellipse'),
                                angle: Number(node.angle || 0),
                                stroke: String(node.stroke || ''),
                                fillColor: String(node.fillColor || ''),
                                opacity: Number(node.opacity || 1.0),
                                mapping_type: String(node.mappingType || ''),
                                is_complex: !!node.isComplex,
                            }});
                            exportableNodes.push(node);
                            continue;
                        }}
                        snapshot.compound_data.push({{
                            compound_id: String(node.id || ''),
                            x: rect.x,
                            y: rect.y,
                            width: rect.width,
                            height: rect.height,
                            label: String(node.displayLabel || node.label || ''),
                            shape_type: String(node.shapeType || 'rect'),
                            angle: Number(node.angle || 0),
                            stroke: String(node.stroke || ''),
                            fillColor: String(node.fillColor || ''),
                            strokeWidth: Number(node.strokeWidth || 1),
                            opacity: Number(node.opacity || 1.0),
                        }});
                        exportableNodes.push(node);
                    }}
                    for (const edge of (local.editableEdges || [])) {{
                        if (!edge) continue;
                        const startNode = closestExportNodeForPoint(exportableNodes, edge.startX, edge.startY);
                        const endNode = closestExportNodeForPoint(exportableNodes, edge.endX, edge.endY);
                        if (!startNode || !endNode) continue;
                        snapshot.arrows.push({{
                            protbox_id_1: String(startNode.id || ''),
                            protbox_id_2: String(endNode.id || ''),
                            protbox_id_1_side: exportArrowSideForPoint(startNode, edge.startX, edge.startY),
                            protbox_id_2_side: exportArrowSideForPoint(endNode, edge.endX, edge.endY),
                            line: String(edge.dashed ? 'dashed' : 'solid'),
                            type: (String(edge.endType || 'arrow') === 'inhibitor')
                                ? 'inhibition'
                                : (String(edge.endType || 'arrow') === 'none' ? 'line' : 'arrow'),
                            startType: String(edge.startType || 'none'),
                            endType: String(edge.endType || 'arrow'),
                            dashed: !!edge.dashed,
                            strokeWidth: Number(edge.strokeWidth || 1.6),
                            stroke: String(edge.stroke || '#0f172a'),
                            opacity: Number(edge.opacity || 0.95),
                            x1: Number(edge.startX || 0),
                            y1: Number(edge.startY || 0),
                            x2: Number(edge.endX || 0),
                            y2: Number(edge.endY || 0),
                            controlX: Number(edge.controlX || 0),
                            controlY: Number(edge.controlY || 0),
                            bendPoints: cloneStateSnapshot(edge.bendPoints || []),
                        }});
                    }}
                    try {{
                        return JSON.parse(JSON.stringify(snapshot));
                    }} catch (_err) {{
                        return snapshot;
                    }}
                }};
                const buildExternalControlState = (local) => {{
                    const selectedNodeIds = getSelectedNodeIds(local);
                    const selectedEdgeIds = getSelectedEdgeIds(local);
                    return {{
                        canUndo: Array.isArray(local.undoStack) && local.undoStack.length > 0,
                        canRedo: Array.isArray(local.redoStack) && local.redoStack.length > 0,
                        canDelete: selectedNodeIds.length > 0 || selectedEdgeIds.length > 0,
                        canAutoConnect: getConnectableSelectedNodes(local).length >= 2,
                        canEdgeEdit: selectedEdgeIds.length > 0,
                        canArrange: selectedNodeIds.length > 0,
                        proteinOvalMode: !!local.proteinOvalMode,
                        edgeResizeMode: !!local.edgeResizeMode,
                        showMissing: stage.classList.contains('cst-show-missing'),
                        mouseMode: String(local.mouseMode || 'drag'),
                        addArrowMode: !!local.addArrowMode,
                        addShapeType: String(local.addShapeType || 'ellipse'),
                        selectedNodeCount: selectedNodeIds.length,
                        selectedEdgeCount: selectedEdgeIds.length,
                    }};
                }};
                const emitExternalControlState = (local) => {{
                    try {{
                        window.dispatchEvent(new CustomEvent('mk-cst-viewer-state', {{
                            detail: {{
                                viewerKey: '{viewer_key}',
                                state: buildExternalControlState(local),
                            }},
                        }}));
                    }} catch (_err) {{}}
                    try {{
                        window.dispatchEvent(new CustomEvent('mk-viewer-controls-ready', {{
                            detail: {{
                                key: '{viewer_key}',
                                state: buildOverlayVisibilityState(local),
                            }},
                        }}));
                    }} catch (_err) {{}}
                }};
                setMissingState(true);
                if (missingButton) {{
                    missingButton.addEventListener('click', () => {{
                        setMissingState(!stage.classList.contains('cst-show-missing'));
                    }});
                }}
                const tooltipShowDelayMs = 1000;
                let hoverTooltipTimer = null;
                let hoverTooltipPending = null;
                const clearHoverTooltipTimer = () => {{
                    if (hoverTooltipTimer) {{
                        window.clearTimeout(hoverTooltipTimer);
                        hoverTooltipTimer = null;
                    }}
                }};
                const hideHoverTooltip = () => {{
                    clearHoverTooltipTimer();
                    hoverTooltipPending = null;
                    if (!hoverTooltip) return;
                    hoverTooltip.style.display = 'none';
                    hoverTooltip.innerHTML = '';
                }};
                const sanitizeTooltipHtml = (rawHtml) => {{
                    const raw = String(rawHtml || '');
                    if (!raw) return '';
                    const template = document.createElement('template');
                    template.innerHTML = raw;
                    const allowed = new Set(['STRONG', 'EM', 'I', 'B', 'U', 'BR', 'SPAN', 'DIV']);
                    const walk = (node) => {{
                        const children = Array.from(node.childNodes || []);
                        for (const child of children) {{
                            if (child.nodeType === Node.ELEMENT_NODE) {{
                                const tag = String(child.tagName || '').toUpperCase();
                                if (!allowed.has(tag)) {{
                                    const replacement = document.createTextNode(child.textContent || '');
                                    child.replaceWith(replacement);
                                    continue;
                                }}
                                const attrs = Array.from(child.attributes || []);
                                for (const attr of attrs) {{
                                    child.removeAttribute(attr.name);
                                }}
                                walk(child);
                            }} else if (child.nodeType === Node.COMMENT_NODE) {{
                                child.remove();
                            }}
                        }}
                    }};
                    walk(template.content);
                    return template.innerHTML;
                }};
                const tooltipInfoFromElement = (element) => {{
                    if (!element || !element.getAttribute) return null;
                    const htmlTip = element.getAttribute('data-tooltip-html');
                    const tip = htmlTip || element.getAttribute('data-tooltip');
                    if (!tip) return null;
                    return {{ tip, useHtml: !!htmlTip, target: element }};
                }};
                const applyTooltipAttrs = (element, textTip, htmlTip) => {{
                    if (!element || !element.setAttribute) return;
                    const plain = String(textTip || '').trim();
                    const rich = String(htmlTip || '').trim();
                    if (plain) {{
                        element.setAttribute('data-tooltip', plain);
                    }} else {{
                        element.removeAttribute('data-tooltip');
                    }}
                    if (rich) {{
                        element.setAttribute('data-tooltip-html', rich);
                    }} else {{
                        element.removeAttribute('data-tooltip-html');
                    }}
                }};
                const resolveCstTooltipTarget = (target) => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (local && local.showTooltips === false) return null;
                    if (!target || !target.closest) return null;
                    return target.closest('[data-tooltip],[data-tooltip-html]');
                }};
                const showHoverTooltipNow = (pending) => {{
                    if (!hoverTooltip || !pending || !pending.tip) return;
                    if (pending.useHtml) {{
                        hoverTooltip.innerHTML = sanitizeTooltipHtml(pending.tip);
                    }} else {{
                        hoverTooltip.textContent = pending.tip;
                    }}
                    hoverTooltip.style.left = pending.left + 'px';
                    hoverTooltip.style.top = pending.top + 'px';
                    hoverTooltip.style.display = 'block';
                }};
                const scheduleHoverTooltip = (info, evt) => {{
                    if (!hoverTooltip || !info || !evt) {{
                        hideHoverTooltip();
                        return false;
                    }}
                    const pos = toViewportPosition(evt);
                    hoverTooltipPending = {{
                        tip: info.tip,
                        useHtml: !!info.useHtml,
                        target: info.target || null,
                        clientX: evt.clientX,
                        clientY: evt.clientY,
                        left: pos.x,
                        top: pos.y,
                    }};
                    clearHoverTooltipTimer();
                    hoverTooltipTimer = window.setTimeout(() => {{
                        hoverTooltipTimer = null;
                        const pending = hoverTooltipPending;
                        if (!pending || !pending.tip) return;
                        const liveTarget = (typeof document.elementFromPoint === 'function')
                            ? document.elementFromPoint(pending.clientX, pending.clientY)
                            : null;
                        const resolvedTarget = resolveCstTooltipTarget(liveTarget);
                        if (!resolvedTarget || (pending.target && resolvedTarget !== pending.target && !pending.target.contains(resolvedTarget))) {{
                            hoverTooltipPending = null;
                            return;
                        }}
                        showHoverTooltipNow(pending);
                    }}, tooltipShowDelayMs);
                    return true;
                }};
                const moveHoverTooltip = (evt) => {{
                    if (!hoverTooltip || !evt) return false;
                    if (hoverTooltipPending) {{
                        const pos = toViewportPosition(evt);
                        hoverTooltipPending.left = pos.x;
                        hoverTooltipPending.top = pos.y;
                        hoverTooltipPending.clientX = evt.clientX;
                        hoverTooltipPending.clientY = evt.clientY;
                        hoverTooltipPending.target = resolveCstTooltipTarget(evt.target) || hoverTooltipPending.target;
                    }}
                    if (hoverTooltip.style.display !== 'block') return false;
                    const pos = toViewportPosition(evt);
                    hoverTooltip.style.left = pos.x + 'px';
                    hoverTooltip.style.top = pos.y + 'px';
                    return true;
                }};
                const updateHoverTooltipFromEvent = (evt) => {{
                    if (!evt || (typeof evt.buttons === 'number' && evt.buttons > 0)) {{
                        hideHoverTooltip();
                        return false;
                    }}
                    const target = resolveCstTooltipTarget(evt.target);
                    const info = tooltipInfoFromElement(target);
                    if (!info || !info.tip) {{
                        hideHoverTooltip();
                        return false;
                    }}
                    if (hoverTooltip.style.display === 'block') {{
                        return moveHoverTooltip(evt);
                    }}
                    return scheduleHoverTooltip(info, evt);
                }};
                if (overlaySvg) {{
                    overlaySvg.addEventListener('mousemove', updateHoverTooltipFromEvent);
                    overlaySvg.addEventListener('mouseleave', hideHoverTooltip);
                    overlaySvg.addEventListener('mousedown', hideHoverTooltip);
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
                            hasSuccessfulRender: false,
                            lastStageWidth: 0,
                            initialRenderAttempts: 0,
                            overlayVersion: null,
                            editableNodes: null,
                            editableEdges: null,
                            groups: [],
                            selectedNodeId: null,
                            selectedNodeIds: [],
                            selectedEdgeId: null,
                            selectedEdgeIds: [],
                            selectedPtmId: null,
                            selectedPtmNodeId: null,
                            disablePdfReader: {str(disable_pdf_reader).lower()},
                            simpleKeggMode: {str(simple_kegg_mode).lower()},
                            showPdfBackground: {str(simple_kegg_mode).lower()},
                            addMode: false,
                            addShapeType: 'ellipse',
                            addArrowMode: false,
                            addArrowPreset: 'arrow',
                            pendingArrowStart: null,
                            pendingArrowPreview: null,
                            batchSizePanelOpen: false,
                            proteinOvalMode: true,
                            proteinOutlineWidth: {prot_outline_width},
                            temporalMode: {str(temporal_mode).lower()},
                            edgeResizeMode: false,
                            dragState: null,
                            globalOpacity: 1.0,
                            zoom: 1,
                            undoStack: [],
                            redoStack: [],
                            isRestoringHistory: false,
                            copiedNodes: [],
                            copiedEdges: [],
                            tabPreviewMode: false,
                            mouseMode: 'drag',
                            mergedTextItems: [],
                            editingTextNodeId: null,
                            editingTextSnapshot: null,
                            contextNodeId: null,
                            contextTextNodeId: null,
                            contextEdgeId: null,
                            contextBendIndex: null,
                            complexPanels: [],
                            autoPtmPlacements: null,
                            autoPtmObstacleEdges: cloneStateSnapshot({ptm_obstacle_edges_json}),
                            showProtboxes: true,
                            showPtms: true,
                            showPtmLabels: true,
                            showPtmSymbols: true,
                            showArrows: true,
                            showCompounds: true,
                            showGroups: true,
                            showTooltips: true,
                        }};
                    }}
                    const local = state['{viewer_key}'];
                    local.proteinOutlineWidth = {prot_outline_width};
                    local.temporalMode = {str(temporal_mode).lower()};
                    local.disablePdfReader = {str(disable_pdf_reader).lower()};
                    local.simpleKeggMode = {str(simple_kegg_mode).lower()};
                    local.showPdfBackground = {str(simple_kegg_mode).lower()};
                    local.autoPtmObstacleEdges = cloneStateSnapshot({ptm_obstacle_edges_json});
                    if (canvas) {{
                        canvas.style.opacity = local.showPdfBackground ? '1' : '0';
                    }}
                    if (fallback) {{
                        fallback.style.opacity = local.showPdfBackground ? '1' : '0';
                    }}
                    if (disablePdfReaderCheckbox) {{
                        disablePdfReaderCheckbox.checked = !!local.disablePdfReader;
                        disablePdfReaderCheckbox.addEventListener('change', () => {{
                            local.disablePdfReader = !!disablePdfReaderCheckbox.checked;
                        }});
                    }}
                    const updateZoomUi = () => {{
                        if (zoomResetButton) {{
                            zoomResetButton.textContent = Math.round((local.zoom || 1) * 100) + '%';
                        }}
                    }};
                    const updateContentLayout = () => {{
                        const stageWidth = Math.max(stage.clientWidth || 0, 10);
                        const stageHeight = Math.max(stage.clientHeight || 0, 10);
                        const baseCssScale = stageWidth / pageWidth;
                        const baseCssHeight = baseCssScale * pageHeight;
                        const zoom = Math.max(0.4, Math.min(3, Number(local.zoom || 1)));
                        const contentWidth = stageWidth * zoom;
                        const contentHeight = baseCssHeight * zoom;
                        stage.style.height = baseCssHeight + 'px';
                        content.style.width = contentWidth + 'px';
                        content.style.height = contentHeight + 'px';
                        canvas.style.width = contentWidth + 'px';
                        canvas.style.height = contentHeight + 'px';
                        if (fallback) {{
                            fallback.style.width = contentWidth + 'px';
                            fallback.style.height = contentHeight + 'px';
                        }}
                        updateZoomUi();
                        return {{ stageWidth, stageHeight, baseCssScale, baseCssHeight, contentWidth, contentHeight, zoom }};
                    }};
                    const applyZoom = (newZoom) => {{
                        const nextZoom = Math.max(0.4, Math.min(3, Number(newZoom || 1)));
                        const prevZoom = Math.max(0.4, Math.min(3, Number(local.zoom || 1)));
                        if (Math.abs(nextZoom - prevZoom) < 0.001) return;
                        const rect = viewport.getBoundingClientRect();
                        const focusX = (viewport.scrollLeft + (rect.width * 0.5)) / Math.max(prevZoom, 0.001);
                        const focusY = (viewport.scrollTop + (rect.height * 0.5)) / Math.max(prevZoom, 0.001);
                        local.zoom = nextZoom;
                        updateContentLayout();
                        viewport.scrollLeft = Math.max(0, (focusX * nextZoom) - (rect.width * 0.5));
                        viewport.scrollTop = Math.max(0, (focusY * nextZoom) - (rect.height * 0.5));
                    }};
                    window.pdfjsLib.GlobalWorkerOptions.workerSrc = 'https://cdnjs.cloudflare.com/ajax/libs/pdf.js/3.11.174/pdf.worker.min.js';
                    const renderNow = () => {{
                        try {{
                            console.info('[MapKinase CST] renderNow start', {{
                                viewerKey: '{viewer_key}',
                                renderVersion: Number(local.renderVersion || 0) + 1,
                                stageWidth: Math.max(stage.clientWidth || 0, 0),
                                hasSuccessfulRender: !!local.hasSuccessfulRender,
                            }});
                        }} catch (_err) {{}}
                        local.renderVersion += 1;
                        const renderVersion = local.renderVersion;
                        local.pdfPromise.then((pdf) => pdf.getPage(1)).then((page) => {{
                            const layout = updateContentLayout();
                            const stageWidth = layout.stageWidth;
                            const baseCssScale = layout.baseCssScale;
                            const dpr = Math.max(window.devicePixelRatio || 1, 1);
                            const renderScaleBoost = 4;
                            const pdfViewport = page.getViewport({{ scale: baseCssScale * dpr * renderScaleBoost }});
                            const ctx = canvas.getContext('2d', {{ alpha: false }});
                            if (!ctx) {{
                                showFallback('canvas_context_missing', true);
                                return;
                            }}
                            local.lastStageWidth = stageWidth;
                            canvas.width = Math.max(1, Math.round(pdfViewport.width));
                            canvas.height = Math.max(1, Math.round(pdfViewport.height));
                            canvas.style.display = 'block';
                            canvas.style.opacity = local.showPdfBackground ? '1' : '0';
                            if (fallback) fallback.style.display = 'none';
                            if (local.renderTask && local.renderTask.cancel) {{
                                try {{ local.renderTask.cancel(); }} catch (_err) {{}}
                            }}
                            local.renderTask = page.render({{ canvasContext: ctx, viewport: pdfViewport }});
                            const renderPromise = local.renderTask.promise.catch((err) => {{
                                    if (err && err.name === 'RenderingCancelledException') return null;
                                    try {{
                                        console.error('[MapKinase CST] renderTask error', err);
                                    }} catch (_err) {{}}
                                    showFallback('render_task_error', false);
                                    return null;
                                }});
                            const textPromise = page.getTextContent().catch(() => ({{ items: [] }}));
                            return Promise.all([renderPromise, textPromise]).then(([, textContent]) => {{
                                if (renderVersion !== local.renderVersion) return null;
                                if (local.overlayVersion !== overlayVersion) {{
                                    local.overlayVersion = overlayVersion;
                                    local.editableNodes = Array.isArray(initialEditableNodes) && initialEditableNodes.length
                                        ? cloneStateSnapshot(initialEditableNodes)
                                        : null;
                                    local.editableEdges = Array.isArray(initialEditableEdges)
                                        ? normalizeEditableEdgeIds(cloneStateSnapshot(initialEditableEdges))
                                        : [];
                                    local.groups = sanitizeGroups(local, Array.isArray(initialGroups) ? cloneStateSnapshot(initialGroups) : []);
                                    local.selectedNodeId = null;
                                    local.selectedNodeIds = [];
                                    local.selectedEdgeId = null;
                                    local.selectedEdgeIds = [];
                                    local.addMode = false;
                                    local.addArrowMode = false;
                                    local.pendingArrowStart = null;
                                    local.pendingArrowPreview = null;
                                    local.autoPtmPlacements = null;
                                    local.undoStack = [];
                                    if (Math.abs(Number(local.globalOpacity || 1.0) - 0.98) < 0.001) {{
                                        local.globalOpacity = 1.0;
                                    }}
                                    if (Array.isArray(local.editableNodes)) {{
                                        for (const node of local.editableNodes) {{
                                            const nodeOpacity = Number(node && node.opacity);
                                            if (Math.abs(nodeOpacity - 0.98) < 0.001) {{
                                                node.opacity = Number(local.globalOpacity || 1.0);
                                            }}
                                        }}
                                    }}
                                }}
                                if (textContent) {{
                                    local.mergedTextItems = buildMergedTextItems((textContent && textContent.items) || []);
                                }} else if (!Array.isArray(local.mergedTextItems)) {{
                                    local.mergedTextItems = [];
                                }}
                                if (!Array.isArray(local.editableNodes) || !local.editableNodes.length) {{
                                    local.editableNodes = renderOverlayNodes((textContent && textContent.items) || []);
                                    for (const node of (local.editableNodes || [])) {{
                                        node.opacity = local.globalOpacity;
                                    }}
                                }}
                                local.editableNodes = applyActiveFcToEditableNodes(local.editableNodes || []);
                                if (AUTO_PTM_RENDER_ENABLED && !Array.isArray(local.autoPtmPlacements)) {{
                                    local.autoPtmPlacements = buildAutoPtmPlacements(local);
                                }}
                                local.hasSuccessfulRender = true;
                                try {{
                                    console.info('[MapKinase CST] renderNow success', {{
                                        viewerKey: '{viewer_key}',
                                        renderVersion,
                                        textItems: Array.isArray((textContent && textContent.items) || []) ? ((textContent && textContent.items) || []).length : 0,
                                        editableNodes: Array.isArray(local.editableNodes) ? local.editableNodes.length : 0,
                                    }});
                                }} catch (_err) {{}}
                                renderEditableOverlay(local);
                                return null;
                            }}).catch((err) => {{
                                if (err && err.name === 'RenderingCancelledException') return null;
                                try {{
                                    console.error('[MapKinase CST] renderNow catch', err);
                                }} catch (_err) {{}}
                                showFallback('render_now_catch', false);
                                return null;
                            }});
                        }}).catch(() => {{
                            try {{
                                console.error('[MapKinase CST] pdfPromise/page load failed');
                            }} catch (_err) {{}}
                            showFallback('pdf_page_load_failed', false);
                        }});
                    }};
                    const queueRender = () => {{
                        if (local.rafId) cancelAnimationFrame(local.rafId);
                        local.rafId = requestAnimationFrame(() => {{
                            local.rafId = null;
                            renderNow();
                        }});
                    }};
                    const startWhenReady = () => {{
                        const stageWidth = Math.max(stage.clientWidth || 0, 0);
                        local.initialRenderAttempts = Number(local.initialRenderAttempts || 0) + 1;
                        if (stageWidth > 20 || local.initialRenderAttempts >= 10) {{
                            queueRender();
                            return;
                        }}
                        window.setTimeout(startWhenReady, 120);
                    }};
                    startWhenReady();
                    if (!local.resizeObserver && window.ResizeObserver) {{
                        local.resizeObserver = new ResizeObserver(() => {{
                            const stageWidth = Math.max(stage.clientWidth || 0, 0);
                            if (Math.abs(stageWidth - Number(local.lastStageWidth || 0)) > 1) {{
                                queueRender();
                            }}
                        }});
                        local.resizeObserver.observe(stage);
                    }}
                    window.addEventListener('resize', queueRender, {{ passive: true }});
                    if (zoomInButton) zoomInButton.addEventListener('click', () => applyZoom((local.zoom || 1) + 0.2));
                    if (zoomOutButton) zoomOutButton.addEventListener('click', () => applyZoom((local.zoom || 1) - 0.2));
                    if (zoomResetButton) zoomResetButton.addEventListener('click', () => applyZoom(1));
                    viewport.addEventListener('wheel', (evt) => {{
                        if (!evt.ctrlKey) return;
                        evt.preventDefault();
                        const step = evt.deltaY < 0 ? 0.12 : -0.12;
                        applyZoom((local.zoom || 1) + step);
                    }}, {{ passive: false }});
                }};
                const waitForPdfJs = (attemptsLeft) => {{
                    if (ensurePdfJs()) {{
                        setupAndRender();
                        return;
                    }}
                    if (attemptsLeft <= 0) {{
                        showFallback('pdfjs_not_available', true);
                        return;
                    }}
                    setTimeout(() => waitForPdfJs(attemptsLeft - 1), 120);
                }};
                const beginDrag = (evt, local, nodeId, kind) => {{
                    const node = findNodeById(local, nodeId);
                    const point = pointerToViewer(evt);
                    if (!node || !point) return;
                    const selectedNodeIds = getSelectedNodeIds(local);
                    const selectedEdgeIds = getSelectedEdgeIds(local);
                    const keepMultiSelection = kind === 'move' && isNodeSelected(local, nodeId) && (selectedNodeIds.length + selectedEdgeIds.length) > 1;
                    const startLocalPoint = rotatePoint(
                        Number(point.overlayX || 0) - Number(node.cx || 0),
                        Number(point.overlayY || 0) - Number(node.cy || 0),
                        -Number(node.angle || 0),
                    );
                    if (!keepMultiSelection) {{
                        setSingleSelection(local, nodeId);
                    }}
                    local.dragState = {{
                        mode: kind,
                        nodeId,
                        startPoint: point,
                        startCx: Number(node.cx || 0),
                        startCy: Number(node.cy || 0),
                        startRx: Number(node.rx || 12),
                        startRy: Number(node.ry || 9),
                        startAngle: Number(node.angle || 0),
                        startLocalPoint,
                        moveNodeStarts: keepMultiSelection
                            ? selectedNodeIds
                                .map((selectedId) => {{
                                    const selectedNode = findNodeById(local, selectedId);
                                    if (!selectedNode) return null;
                                    return {{
                                        id: String(selectedNode.id || ''),
                                        cx: Number(selectedNode.cx || 0),
                                        cy: Number(selectedNode.cy || 0),
                                    }};
                                }})
                                .filter(Boolean)
                            : [],
                        moveEdgeStarts: keepMultiSelection
                            ? selectedEdgeIds
                                .map((selectedId) => {{
                                    const selectedEdge = findEdgeById(local, selectedId);
                                    if (!selectedEdge) return null;
                                    return {{
                                        id: String(selectedEdge.id || ''),
                                        startX: Number(selectedEdge.startX || 0),
                                        startY: Number(selectedEdge.startY || 0),
                                        endX: Number(selectedEdge.endX || 0),
                                        endY: Number(selectedEdge.endY || 0),
                                        controlX: Number(selectedEdge.controlX || 0),
                                        controlY: Number(selectedEdge.controlY || 0),
                                        bendPoints: getEdgeBendPoints(selectedEdge).map((item) => ({{
                                            x: Number(item.x || 0),
                                            y: Number(item.y || 0),
                                        }})),
                                    }};
                                }})
                                .filter(Boolean)
                            : [],
                        beforeSnapshot: buildSnapshot(local),
                    }};
                    setAddState(local, false);
                    renderEditableOverlay(local);
                    evt.preventDefault();
                    evt.stopPropagation();
                }};
                const beginEdgeDrag = (evt, local, edgeId, kind, bendIndex = null) => {{
                    const edge = findEdgeById(local, edgeId);
                    const point = pointerToViewer(evt);
                    if (!edge || !point) return;
                    const bends = getEdgeBendPoints(edge);
                    const selectedNodeIds = getSelectedNodeIds(local);
                    const selectedEdgeIds = getSelectedEdgeIds(local);
                    const keepMultiSelection = kind === 'edge-move' && isEntrySelected(local, 'edge', edgeId) && (selectedNodeIds.length + selectedEdgeIds.length) > 1;
                    if (!keepMultiSelection) {{
                        setSingleEdgeSelection(local, edgeId);
                    }}
                    local.dragState = {{
                        mode: kind,
                        edgeId,
                        bendIndex: bendIndex !== null && bendIndex !== undefined ? Number(bendIndex) : null,
                        startPoint: point,
                        startStartX: Number(edge.startX || 0),
                        startStartY: Number(edge.startY || 0),
                        startEndX: Number(edge.endX || 0),
                        startEndY: Number(edge.endY || 0),
                        startControlX: Number(edge.controlX || 0),
                        startControlY: Number(edge.controlY || 0),
                        startControlOffsetFromStartX: Number(edge.controlX || 0) - Number(edge.startX || 0),
                        startControlOffsetFromStartY: Number(edge.controlY || 0) - Number(edge.startY || 0),
                        startControlOffsetFromEndX: Number(edge.controlX || 0) - Number(edge.endX || 0),
                        startControlOffsetFromEndY: Number(edge.controlY || 0) - Number(edge.endY || 0),
                        startBendPoints: bends.map((item) => ({{
                            x: Number(item.x || 0),
                            y: Number(item.y || 0),
                        }})),
                        moveNodeStarts: keepMultiSelection
                            ? selectedNodeIds
                                .map((selectedId) => {{
                                    const selectedNode = findNodeById(local, selectedId);
                                    if (!selectedNode) return null;
                                    return {{
                                        id: String(selectedNode.id || ''),
                                        cx: Number(selectedNode.cx || 0),
                                        cy: Number(selectedNode.cy || 0),
                                    }};
                                }})
                                .filter(Boolean)
                            : [],
                        moveEdgeStarts: keepMultiSelection
                            ? selectedEdgeIds
                                .map((selectedId) => {{
                                    const selectedEdge = findEdgeById(local, selectedId);
                                    if (!selectedEdge) return null;
                                    return {{
                                        id: String(selectedEdge.id || ''),
                                        startX: Number(selectedEdge.startX || 0),
                                        startY: Number(selectedEdge.startY || 0),
                                        endX: Number(selectedEdge.endX || 0),
                                        endY: Number(selectedEdge.endY || 0),
                                        controlX: Number(selectedEdge.controlX || 0),
                                        controlY: Number(selectedEdge.controlY || 0),
                                        bendPoints: getEdgeBendPoints(selectedEdge).map((item) => ({{
                                            x: Number(item.x || 0),
                                            y: Number(item.y || 0),
                                        }})),
                                    }};
                                }})
                                .filter(Boolean)
                            : [],
                        beforeSnapshot: buildSnapshot(local),
                    }};
                    setAddState(local, false);
                    setAddArrowState(local, false);
                    renderEditableOverlay(local);
                    evt.preventDefault();
                    evt.stopPropagation();
                }};
                const beginComplexPanelDrag = (evt, local, panelId) => {{
                    const point = pointerToViewer(evt);
                    if (!point) return;
                    const panel = Array.isArray(local.complexPanels)
                        ? local.complexPanels.find((item) => String(item.id || '') === String(panelId || ''))
                        : null;
                    if (!panel) return;
                    local.dragState = {{
                        mode: 'complex-panel-move',
                        panelId: String(panel.id || ''),
                        startPoint: point,
                        startPanelX: Number(panel.x || 0),
                        startPanelY: Number(panel.y || 0),
                        beforeSnapshot: null,
                    }};
                    evt.preventDefault();
                    evt.stopPropagation();
                }};
                const beginAutoPtmDrag = (evt, local, ptmId, nodeId) => {{
                    const point = pointerToViewer(evt);
                    const node = findNodeById(local, nodeId);
                    const placement = Array.isArray(local.autoPtmPlacements)
                        ? local.autoPtmPlacements.find((item) => String(item.id || '') === String(ptmId || ''))
                        : null;
                    if (!point || !node || !placement) return;
                    local.dragState = {{
                        mode: 'ptm-move',
                        ptmId: String(ptmId || ''),
                        nodeId: String(nodeId || ''),
                        startPoint: point,
                        startPtmAngleDeg: Number(placement.angleDeg || 0),
                        startPtmOutwardFactor: Number(placement.outwardFactor || 0.34),
                        beforeSnapshot: null,
                    }};
                    evt.preventDefault();
                    evt.stopPropagation();
                }};
                const beginAutoPtmLabelDrag = (evt, local, ptmId, nodeId) => {{
                    const point = pointerToViewer(evt);
                    const placement = Array.isArray(local.autoPtmPlacements)
                        ? local.autoPtmPlacements.find((item) => String(item && item.id || '') === String(ptmId || ''))
                        : null;
                    if (!point || !placement) return;
                    local.dragState = {{
                        mode: 'ptm-label-move',
                        ptmId: String(ptmId || ''),
                        nodeId: String(nodeId || ''),
                        startPoint: point,
                        startLabelOffsetX: Number(placement.labelOffsetX || 0),
                        startLabelOffsetY: Number(placement.labelOffsetY || 0),
                        beforeSnapshot: buildSnapshot(local),
                    }};
                    evt.preventDefault();
                    evt.stopPropagation();
                }};
                const updatePendingArrowPreview = (local, point) => {{
                    if (!local.addArrowMode || !local.pendingArrowStart || !point) return;
                    const preview = createDefaultEdge(
                        Number(local.pendingArrowStart.x || 0),
                        Number(local.pendingArrowStart.y || 0),
                        Number(point.overlayX || 0),
                        Number(point.overlayY || 0),
                    );
                    applyEdgePreset(preview, String(local.addArrowPreset || 'arrow'));
                    preview.id = 'cst_edge_preview';
                    preview.stroke = '#0f172a';
                    preview.strokeWidth = 1.0;
                    preview.opacity = 0.75;
                    local.pendingArrowPreview = preview;
                }};
                const beginMarqueeSelection = (evt, local) => {{
                    const point = pointerToViewer(evt);
                    if (!point) return;
                    local.dragState = {{
                        mode: 'marquee',
                        nodeId: null,
                        startPoint: point,
                        currentPoint: point,
                        beforeSnapshot: null,
                    }};
                    setAddState(local, false);
                    renderEditableOverlay(local);
                    evt.preventDefault();
                    evt.stopPropagation();
                }};
                const updateMarqueeSelection = (local) => {{
                    const drag = local.dragState;
                    if (!drag || drag.mode !== 'marquee') return;
                    const minX = Math.min(Number(drag.startPoint.overlayX || 0), Number(drag.currentPoint && drag.currentPoint.overlayX || 0));
                    const maxX = Math.max(Number(drag.startPoint.overlayX || 0), Number(drag.currentPoint && drag.currentPoint.overlayX || 0));
                    const minY = Math.min(Number(drag.startPoint.overlayY || 0), Number(drag.currentPoint && drag.currentPoint.overlayY || 0));
                    const maxY = Math.max(Number(drag.startPoint.overlayY || 0), Number(drag.currentPoint && drag.currentPoint.overlayY || 0));
                    const selectedIds = (local.editableNodes || [])
                        .filter((node) => {{
                            const bounds = getNodeBounds(node);
                            if (!bounds) return false;
                            return !(bounds.right < minX || bounds.left > maxX || bounds.bottom < minY || bounds.top > maxY);
                        }})
                        .map((node) => node.id);
                    const selectedEdgeIds = (local.editableEdges || [])
                        .filter((edge) => {{
                            const bounds = getEdgeBounds(edge);
                            return !(bounds.right < minX || bounds.left > maxX || bounds.bottom < minY || bounds.top > maxY);
                        }})
                        .map((edge) => edge.id);
                    local.selectedNodeIds = selectedIds;
                    local.selectedNodeId = selectedIds.length ? selectedIds[selectedIds.length - 1] : null;
                    local.selectedEdgeIds = selectedEdgeIds;
                    local.selectedEdgeId = selectedEdgeIds.length ? selectedEdgeIds[selectedEdgeIds.length - 1] : null;
                }};
                const handlePointerMove = (evt) => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    const point = pointerToViewer(evt);
                    if (!local) return;
                    if (local.addArrowMode && local.pendingArrowStart && point) {{
                        updatePendingArrowPreview(local, point);
                        renderEditableOverlay(local);
                    }}
                    if (!local.dragState) return;
                    const drag = local.dragState;
                    if (!point) return;
                    if (drag.mode === 'marquee') {{
                        drag.currentPoint = point;
                        updateMarqueeSelection(local);
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if (drag.mode === 'complex-panel-move') {{
                        const panel = Array.isArray(local.complexPanels)
                            ? local.complexPanels.find((item) => String(item.id || '') === String(drag.panelId || ''))
                            : null;
                        if (!panel) return;
                        panel.x = Number(drag.startPanelX || 0) + (Number(point.overlayX || 0) - Number(drag.startPoint.overlayX || 0));
                        panel.y = Number(drag.startPanelY || 0) + (Number(point.overlayY || 0) - Number(drag.startPoint.overlayY || 0));
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if (drag.mode === 'ptm-move') {{
                        const node = findNodeById(local, drag.nodeId);
                        const placement = Array.isArray(local.autoPtmPlacements)
                            ? local.autoPtmPlacements.find((item) => String(item.id || '') === String(drag.ptmId || ''))
                            : null;
                        if (!node || !placement) return;
                        const dx = Number(point.overlayX || 0) - Number(node.cx || 0);
                        const dy = Number(point.overlayY || 0) - Number(node.cy || 0);
                        if (Math.hypot(dx, dy) > 0.001) {{
                            placement.angleDeg = (Math.atan2(dy, dx) * 180) / Math.PI;
                        }}
                        placement.outwardFactor = Math.max(0.08, Math.min(1.05, Number(placement.outwardFactor || drag.startPtmOutwardFactor || 0.34)));
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if (drag.mode === 'ptm-label-move') {{
                        const placement = Array.isArray(local.autoPtmPlacements)
                            ? local.autoPtmPlacements.find((item) => String(item.id || '') === String(drag.ptmId || ''))
                            : null;
                        if (!placement) return;
                        placement.labelOffsetX = Number(drag.startLabelOffsetX || 0) + (Number(point.overlayX || 0) - Number(drag.startPoint.overlayX || 0));
                        placement.labelOffsetY = Number(drag.startLabelOffsetY || 0) + (Number(point.overlayY || 0) - Number(drag.startPoint.overlayY || 0));
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if (drag.mode === 'edge-move') {{
                        const dx = Number(point.overlayX || 0) - Number(drag.startPoint.overlayX || 0);
                        const dy = Number(point.overlayY || 0) - Number(drag.startPoint.overlayY || 0);
                        if ((Array.isArray(drag.moveNodeStarts) && drag.moveNodeStarts.length) || (Array.isArray(drag.moveEdgeStarts) && drag.moveEdgeStarts.length > 1)) {{
                            (Array.isArray(drag.moveNodeStarts) ? drag.moveNodeStarts : []).forEach((item) => {{
                                const moveNode = findNodeById(local, item.id);
                                if (!moveNode) return;
                                moveNode.cx = Number(item.cx || 0) + dx;
                                moveNode.cy = Number(item.cy || 0) + dy;
                            }});
                            (Array.isArray(drag.moveEdgeStarts) ? drag.moveEdgeStarts : []).forEach((item) => {{
                                const moveEdge = findEdgeById(local, item.id);
                                if (!moveEdge) return;
                                moveEdge.startX = Number(item.startX || 0) + dx;
                                moveEdge.startY = Number(item.startY || 0) + dy;
                                moveEdge.endX = Number(item.endX || 0) + dx;
                                moveEdge.endY = Number(item.endY || 0) + dy;
                                moveEdge.controlX = Number(item.controlX || 0) + dx;
                                moveEdge.controlY = Number(item.controlY || 0) + dy;
                                moveEdge.bendPoints = (Array.isArray(item.bendPoints) ? item.bendPoints : []).map((pointItem) => ({{
                                    x: Number(pointItem.x || 0) + dx,
                                    y: Number(pointItem.y || 0) + dy,
                                }}));
                                if (!isEdgeBent(moveEdge)) {{
                                    const midpoint = getStraightEdgeMidpoint(moveEdge);
                                    moveEdge.controlX = Number(midpoint.x || 0);
                                    moveEdge.controlY = Number(midpoint.y || 0);
                                    moveEdge.bendPoints = [];
                                }}
                            }});
                        }} else {{
                            const edge = findEdgeById(local, drag.edgeId);
                            if (!edge) return;
                            edge.startX = Number(drag.startStartX || 0) + dx;
                            edge.startY = Number(drag.startStartY || 0) + dy;
                            edge.endX = Number(drag.startEndX || 0) + dx;
                            edge.endY = Number(drag.startEndY || 0) + dy;
                            edge.controlX = Number(drag.startControlX || 0) + dx;
                            edge.controlY = Number(drag.startControlY || 0) + dy;
                            edge.bendPoints = (Array.isArray(drag.startBendPoints) ? drag.startBendPoints : []).map((item) => ({{
                                x: Number(item.x || 0) + dx,
                                y: Number(item.y || 0) + dy,
                            }}));
                            if (!isEdgeBent(edge)) {{
                                const midpoint = getStraightEdgeMidpoint(edge);
                                edge.controlX = Number(midpoint.x || 0);
                                edge.controlY = Number(midpoint.y || 0);
                                edge.bendPoints = [];
                            }}
                        }}
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if (drag.mode === 'edge-start') {{
                        const edge = findEdgeById(local, drag.edgeId);
                        if (!edge) return;
                        edge.startX = Number(point.overlayX || 0);
                        edge.startY = Number(point.overlayY || 0);
                        edge.controlX = edge.startX + Number(drag.startControlOffsetFromStartX || 0);
                        edge.controlY = edge.startY + Number(drag.startControlOffsetFromStartY || 0);
                        edge.bendPoints = (Array.isArray(drag.startBendPoints) ? drag.startBendPoints : []).map((item) => ({{
                            x: edge.startX + (Number(item.x || 0) - Number(drag.startStartX || 0)),
                            y: edge.startY + (Number(item.y || 0) - Number(drag.startStartY || 0)),
                        }}));
                        if (!isEdgeBent(edge)) {{
                            const midpoint = getStraightEdgeMidpoint(edge);
                            edge.controlX = Number(midpoint.x || 0);
                            edge.controlY = Number(midpoint.y || 0);
                            edge.bendPoints = [];
                        }}
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if (drag.mode === 'edge-end') {{
                        const edge = findEdgeById(local, drag.edgeId);
                        if (!edge) return;
                        edge.endX = Number(point.overlayX || 0);
                        edge.endY = Number(point.overlayY || 0);
                        edge.controlX = edge.endX + Number(drag.startControlOffsetFromEndX || 0);
                        edge.controlY = edge.endY + Number(drag.startControlOffsetFromEndY || 0);
                        edge.bendPoints = (Array.isArray(drag.startBendPoints) ? drag.startBendPoints : []).map((item) => ({{
                            x: edge.endX + (Number(item.x || 0) - Number(drag.startEndX || 0)),
                            y: edge.endY + (Number(item.y || 0) - Number(drag.startEndY || 0)),
                        }}));
                        if (!isEdgeBent(edge)) {{
                            const midpoint = getStraightEdgeMidpoint(edge);
                            edge.controlX = Number(midpoint.x || 0);
                            edge.controlY = Number(midpoint.y || 0);
                            edge.bendPoints = [];
                        }}
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if (drag.mode === 'edge-bend') {{
                        const edge = findEdgeById(local, drag.edgeId);
                        if (!edge) return;
                        const bends = getEdgeBendPoints(edge);
                        const bendIndex = Math.max(0, Math.min(bends.length - 1, Number(drag.bendIndex || 0)));
                        edge.isBent = true;
                        bends[bendIndex] = {{
                            x: Number(point.overlayX || 0),
                            y: Number(point.overlayY || 0),
                        }};
                        edge.bendPoints = bends;
                        edge.controlX = Number(bends[bendIndex].x || 0);
                        edge.controlY = Number(bends[bendIndex].y || 0);
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    const node = findNodeById(local, local.dragState.nodeId);
                    if (!node) return;
                    if (drag.mode === 'move') {{
                        const dx = Number(point.overlayX || 0) - Number(drag.startPoint.overlayX || 0);
                        const dy = Number(point.overlayY || 0) - Number(drag.startPoint.overlayY || 0);
                        if ((Array.isArray(drag.moveNodeStarts) && drag.moveNodeStarts.length > 1) || (Array.isArray(drag.moveEdgeStarts) && drag.moveEdgeStarts.length)) {{
                            (Array.isArray(drag.moveNodeStarts) ? drag.moveNodeStarts : []).forEach((item) => {{
                                const moveNode = findNodeById(local, item.id);
                                if (!moveNode) return;
                                moveNode.cx = Number(item.cx || 0) + dx;
                                moveNode.cy = Number(item.cy || 0) + dy;
                            }});
                            (Array.isArray(drag.moveEdgeStarts) ? drag.moveEdgeStarts : []).forEach((item) => {{
                                const moveEdge = findEdgeById(local, item.id);
                                if (!moveEdge) return;
                                moveEdge.startX = Number(item.startX || 0) + dx;
                                moveEdge.startY = Number(item.startY || 0) + dy;
                                moveEdge.endX = Number(item.endX || 0) + dx;
                                moveEdge.endY = Number(item.endY || 0) + dy;
                                moveEdge.controlX = Number(item.controlX || 0) + dx;
                                moveEdge.controlY = Number(item.controlY || 0) + dy;
                                moveEdge.bendPoints = (Array.isArray(item.bendPoints) ? item.bendPoints : []).map((pointItem) => ({{
                                    x: Number(pointItem.x || 0) + dx,
                                    y: Number(pointItem.y || 0) + dy,
                                }}));
                                if (!isEdgeBent(moveEdge)) {{
                                    const midpoint = getStraightEdgeMidpoint(moveEdge);
                                    moveEdge.controlX = Number(midpoint.x || 0);
                                    moveEdge.controlY = Number(midpoint.y || 0);
                                    moveEdge.bendPoints = [];
                                }}
                            }});
                        }} else {{
                            const nextCx = drag.startCx + dx;
                            const nextCy = drag.startCy + dy;
                            if (String(node.shapeType || 'ellipse') === 'text') {{
                                node.cx = nextCx;
                                node.cy = nextCy;
                            }} else if (evt.shiftKey) {{
                                node.cx = nextCx;
                                node.cy = nextCy;
                            }} else {{
                                const snapped = applyRectMoveSnap(local, node, nextCx, nextCy);
                                node.cx = snapped.cx;
                                node.cy = snapped.cy;
                            }}
                        }}
                    }} else if (drag.mode === 'x') {{
                        const localPoint = toLocalAxes(point, node);
                        node.rx = Math.max(4, Math.abs(localPoint.x));
                    }} else if (drag.mode === 'y') {{
                        const localPoint = toLocalAxes(point, node);
                        node.ry = Math.max(4, Math.abs(localPoint.y));
                    }} else if (drag.mode === 'right') {{
                        const localPoint = rotatePoint(
                            Number(point.overlayX || 0) - drag.startCx,
                            Number(point.overlayY || 0) - drag.startCy,
                            -drag.startAngle,
                        );
                        const leftLocal = -drag.startRx;
                        const rightLocal = Math.max(4, localPoint.x);
                        const localCenter = (leftLocal + rightLocal) * 0.5;
                        const centerShift = rotatePoint(localCenter, 0, drag.startAngle);
                        node.cx = drag.startCx + centerShift.x;
                        node.cy = drag.startCy + centerShift.y;
                        node.rx = Math.max(4, (rightLocal - leftLocal) * 0.5);
                    }} else if (drag.mode === 'left') {{
                        const localPoint = rotatePoint(
                            Number(point.overlayX || 0) - drag.startCx,
                            Number(point.overlayY || 0) - drag.startCy,
                            -drag.startAngle,
                        );
                        const rightLocal = drag.startRx;
                        const leftLocal = Math.min(-4, localPoint.x);
                        const localCenter = (leftLocal + rightLocal) * 0.5;
                        const centerShift = rotatePoint(localCenter, 0, drag.startAngle);
                        node.cx = drag.startCx + centerShift.x;
                        node.cy = drag.startCy + centerShift.y;
                        node.rx = Math.max(4, (rightLocal - leftLocal) * 0.5);
                    }} else if (drag.mode === 'top') {{
                        const localPoint = rotatePoint(
                            Number(point.overlayX || 0) - drag.startCx,
                            Number(point.overlayY || 0) - drag.startCy,
                            -drag.startAngle,
                        );
                        const bottomLocal = drag.startRy;
                        const topLocal = Math.min(-4, localPoint.y);
                        const localCenter = (topLocal + bottomLocal) * 0.5;
                        const centerShift = rotatePoint(0, localCenter, drag.startAngle);
                        node.cx = drag.startCx + centerShift.x;
                        node.cy = drag.startCy + centerShift.y;
                        node.ry = Math.max(4, (bottomLocal - topLocal) * 0.5);
                    }} else if (drag.mode === 'bottom') {{
                        const localPoint = rotatePoint(
                            Number(point.overlayX || 0) - drag.startCx,
                            Number(point.overlayY || 0) - drag.startCy,
                            -drag.startAngle,
                        );
                        const topLocal = -drag.startRy;
                        const bottomLocal = Math.max(4, localPoint.y);
                        const localCenter = (topLocal + bottomLocal) * 0.5;
                        const centerShift = rotatePoint(0, localCenter, drag.startAngle);
                        node.cx = drag.startCx + centerShift.x;
                        node.cy = drag.startCy + centerShift.y;
                        node.ry = Math.max(4, (bottomLocal - topLocal) * 0.5);
                    }} else if (drag.mode === 'rotate') {{
                        const angle = (Math.atan2(point.overlayY - Number(node.cy || 0), point.overlayX - Number(node.cx || 0)) * 180) / Math.PI;
                        node.angle = evt.shiftKey ? (angle + 90) : snapRotationAngle(angle + 90);
                    }}
                    renderEditableOverlay(local);
                    evt.preventDefault();
                }};
                const endPointerInteraction = () => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (!local) return;
                    const drag = local.dragState;
                    if (drag && drag.beforeSnapshot) {{
                        const before = JSON.stringify(drag.beforeSnapshot);
                        const after = JSON.stringify(buildSnapshot(local));
                        if (before !== after) {{
                            pushUndoSnapshot(local, drag.beforeSnapshot);
                        }}
                    }}
                    local.dragState = null;
                    if (drag && drag.mode === 'marquee') {{
                        renderEditableOverlay(local);
                    }}
                }};
                overlaySvg.addEventListener('pointerdown', (evt) => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (!local) return;
                    if (local.editingTextNodeId) {{
                        finishInlineTextEdit(local, true);
                    }}
                    if (complexMenu && !complexMenu.contains(evt.target)) {{
                        closeComplexMenu();
                    }}
                    if (edgeBendMenu && !edgeBendMenu.contains(evt.target)) {{
                        closeEdgeBendMenu();
                    }}
                    if (textMenu && !textMenu.contains(evt.target)) {{
                        closeTextMenu();
                    }}
                    if (!(evt.target && evt.target.closest && evt.target.closest('.cst-node-context-menu'))) {{
                        removeNodeContextMenu();
                    }}
                    if (evt.button === 2) {{
                        return;
                    }}
                    const target = evt.target;
                    const role = target && target.dataset ? String(target.dataset.role || '') : '';
                    const panelId = target && target.dataset ? String(target.dataset.panelId || '') : '';
                    const nodeId = target && target.dataset ? String(target.dataset.nodeId || '') : '';
                    const edgeId = target && target.dataset ? String(target.dataset.edgeId || '') : '';
                    const ptmId = target && target.dataset ? String(target.dataset.ptmId || '') : '';
                    const bendIndex = target && target.dataset ? Number(target.dataset.bendIndex || 0) : 0;
                    const isMultiToggle = !!(evt.ctrlKey || evt.metaKey);
                    const hadSelectedEdge = !!local.selectedEdgeId;
                    if (local.addArrowMode && local.pendingArrowStart) {{
                        setBatchSizePanelOpen(local, false);
                        const point = pointerToViewer(evt);
                        if (!point) return;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        const newEdge = createDefaultEdge(
                            Number(local.pendingArrowStart.x || 0),
                            Number(local.pendingArrowStart.y || 0),
                            Number(point.overlayX || 0),
                            Number(point.overlayY || 0),
                        );
                        applyEdgePreset(newEdge, String(local.addArrowPreset || 'arrow'));
                        newEdge.opacity = Number(local.globalOpacity || 1.0);
                        local.editableEdges = Array.isArray(local.editableEdges) ? local.editableEdges.concat([newEdge]) : [newEdge];
                        setSingleEdgeSelection(local, newEdge.id);
                        local.pendingArrowStart = null;
                        local.pendingArrowPreview = null;
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        evt.stopPropagation();
                        return;
                    }}
                    if (role === 'complex-close' && panelId) {{
                        closeComplexPanel(local, panelId);
                        evt.preventDefault();
                        evt.stopPropagation();
                        return;
                    }}
                    if (role === 'complex-panel-drag' && panelId) {{
                        beginComplexPanelDrag(evt, local, panelId);
                        return;
                    }}
                    if (role === 'ptm-body' && ptmId && nodeId) {{
                        setSinglePtmSelection(local, nodeId, ptmId);
                        beginAutoPtmDrag(evt, local, ptmId, nodeId);
                        return;
                    }}
                    if (role === 'ptm-label' && ptmId && nodeId) {{
                        setSinglePtmSelection(local, nodeId, ptmId);
                        beginAutoPtmLabelDrag(evt, local, ptmId, nodeId);
                        return;
                    }}
                    if ((role === 'edge-body' || role === 'edge-start' || role === 'edge-end' || role === 'edge-bend') && edgeId) {{
                        if (evt.button === 2 && role === 'edge-bend') {{
                            return;
                        }}
                        if (role === 'edge-body' && !isEntrySelected(local, 'edge', edgeId)) {{
                            const expanded = expandEntriesWithGroups(local, [{{ type: 'edge', id: edgeId }}]);
                            setSelectedEntries(local, expanded, {{ type: 'edge', id: edgeId }});
                            renderEditableOverlay(local);
                        }}
                        setBatchSizePanelOpen(local, false);
                        if (role === 'edge-body') beginEdgeDrag(evt, local, edgeId, 'edge-move');
                        if (role === 'edge-start') beginEdgeDrag(evt, local, edgeId, 'edge-start');
                        if (role === 'edge-end') beginEdgeDrag(evt, local, edgeId, 'edge-end');
                        if (role === 'edge-bend') beginEdgeDrag(evt, local, edgeId, 'edge-bend', bendIndex);
                        return;
                    }}
                    if (!role && (isMultiToggle || String(local.mouseMode || 'drag') === 'selection') && !local.addMode) {{
                        setBatchSizePanelOpen(local, false);
                        beginMarqueeSelection(evt, local);
                        return;
                    }}
                    if (role === 'body' && nodeId && isMultiToggle) {{
                        toggleNodeSelection(local, nodeId);
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        evt.stopPropagation();
                        return;
                    }}
                    if (role === 'body' && nodeId) {{
                        if (!isNodeSelected(local, nodeId)) {{
                            const expanded = expandEntriesWithGroups(local, [{{ type: 'node', id: nodeId }}]);
                            setSelectedEntries(local, expanded, {{ type: 'node', id: nodeId }});
                            renderEditableOverlay(local);
                        }}
                        beginDrag(evt, local, nodeId, 'move');
                        return;
                    }}
                    if ((role === 'x' || role === 'y' || role === 'rotate') && nodeId) {{
                        beginDrag(evt, local, nodeId, role);
                        return;
                    }}
                    if ((role === 'left' || role === 'right' || role === 'top' || role === 'bottom') && nodeId) {{
                        beginDrag(evt, local, nodeId, role);
                        return;
                    }}
                    if (local.addArrowMode) {{
                        setBatchSizePanelOpen(local, false);
                        const point = pointerToViewer(evt);
                        if (!point) return;
                        local.pendingArrowStart = {{ x: Number(point.overlayX || 0), y: Number(point.overlayY || 0) }};
                        updatePendingArrowPreview(local, point);
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        evt.stopPropagation();
                        return;
                    }}
                    if (local.addMode) {{
                        setBatchSizePanelOpen(local, false);
                        const point = pointerToViewer(evt);
                        if (!point) return;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        const shapeType = String(local.addShapeType || 'ellipse');
                        const isRect = shapeType === 'rect' || shapeType === 'square';
                        const isRounded = shapeType === 'rounded';
                        const isRectLike = isRect || isRounded;
                        const isBracket = shapeType === 'bracket';
                        const isText = shapeType === 'text';
                        const newNode = {{
                            id: addNodeId(),
                            originalLabel: '',
                            displayLabel: isText ? 'Text' : '',
                            label: isText ? 'Text' : '',
                            matchedGene: '',
                            matchedUniprot: '',
                            foldText: '',
                            hasDatasetMatch: false,
                            cx: point.overlayX,
                            cy: point.overlayY,
                            rx: isText ? 24 : 60,
                            ry: isText ? 8 : 40,
                            shapeType: isText ? 'text' : (isBracket ? 'bracket' : (isRounded ? 'rounded' : (isRect ? 'square' : 'ellipse'))),
                            angle: 0,
                            strokeWidth: isText ? 1.0 : 1.0,
                            stroke: isText ? 'rgb(71, 85, 105)' : '#000000',
                            fillColor: isText ? 'transparent' : (isBracket ? 'transparent' : '#f5f5f5'),
                            opacity: local.globalOpacity || 1.0,
                            annotation: '',
                            annotationCommitted: false,
                            isCustom: true,
                            userCreated: true,
                            isDrawingShape: !isText,
                            className: 'cst-overlay-ellipse',
                            title: isText ? 'Custom text box' : 'Shape',
                            fontSize: isText ? 11 : undefined,
                            fontWeight: isText ? '600' : undefined,
                            fontFamily: isText ? '"Segoe UI", Arial, sans-serif' : undefined,
                            textColor: isText ? '#0f172a' : undefined,
                            textAlign: isText ? 'center' : undefined,
                        }};
                        local.editableNodes = local.editableNodes || [];
                        if (newNode.isDrawingShape) {{
                            const firstNonShapeIndex = local.editableNodes.findIndex((item) => !(item && item.isDrawingShape));
                            if (firstNonShapeIndex >= 0) {{
                                local.editableNodes.splice(firstNonShapeIndex, 0, newNode);
                            }} else {{
                                local.editableNodes.push(newNode);
                            }}
                        }} else {{
                            local.editableNodes.push(newNode);
                        }}
                        setSingleSelection(local, newNode.id);
                        if (isText) {{
                            local.edgeResizeMode = false;
                        }}
                        renderEditableOverlay(local);
                        if (isText) {{
                            window.setTimeout(() => {{
                                const latestLocal = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                                if (!latestLocal) return;
                                startInlineTextEdit(latestLocal, newNode.id);
                            }}, 0);
                        }}
                        evt.preventDefault();
                        evt.stopPropagation();
                        return;
                    }}
                    if (!isMultiToggle && (getSelectedNodeIds(local).length || hadSelectedEdge)) {{
                        setBatchSizePanelOpen(local, false);
                        setSingleSelection(local, null);
                        local.selectedEdgeId = null;
                        renderEditableOverlay(local);
                    }}
                }});
                overlaySvg.addEventListener('dblclick', (evt) => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (!local) return;
                    let target = evt.target;
                    let role = '';
                    let nodeId = '';
                    while (target && target !== overlaySvg) {{
                        role = target && target.dataset ? String(target.dataset.role || '') : '';
                        nodeId = target && target.dataset ? String(target.dataset.nodeId || '') : '';
                        if (nodeId && (role === 'body' || role === 'label')) break;
                        target = target.parentNode;
                        role = '';
                        nodeId = '';
                    }}
                    if (!nodeId || (role !== 'body' && role !== 'label')) return;
                    const node = findNodeById(local, nodeId);
                    if (!node || String(node.shapeType || 'ellipse') !== 'text') return;
                    evt.preventDefault();
                    evt.stopPropagation();
                    setSingleSelection(local, nodeId);
                    startInlineTextEdit(local, nodeId);
                }});
                window.addEventListener('pointermove', handlePointerMove, {{ passive: false }});
                window.addEventListener('pointerup', endPointerInteraction, {{ passive: true }});
                if (addButton) {{
                    addButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        const wasActive = !!local.addMode && String(local.addShapeType || 'ellipse') === 'ellipse';
                        local.addShapeType = 'ellipse';
                        setAddArrowState(local, false);
                        setAddState(local, !wasActive);
                    }});
                }}
                if (addRectButton) {{
                    addRectButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        const wasActive = !!local.addMode && String(local.addShapeType || 'ellipse') === 'rect';
                        local.addShapeType = 'rect';
                        setAddArrowState(local, false);
                        setAddState(local, !wasActive);
                    }});
                }}
                if (addBracketButton) {{
                    addBracketButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        const wasActive = !!local.addMode && String(local.addShapeType || 'ellipse') === 'bracket';
                        local.addShapeType = 'bracket';
                        setAddArrowState(local, false);
                        setAddState(local, !wasActive);
                    }});
                }}
                if (addTextButton) {{
                    addTextButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        const wasActive = !!local.addMode && String(local.addShapeType || 'ellipse') === 'text';
                        local.addShapeType = 'text';
                        setAddArrowState(local, false);
                        setAddState(local, !wasActive);
                    }});
                }}
                    if (addArrowButton) {{
                        addArrowButton.addEventListener('click', () => {{
                            const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                            if (!local) return;
                            const wasActive = !!local.addArrowMode;
                            local.addArrowPreset = 'arrow';
                            setAddState(local, false);
                            setAddArrowState(local, !wasActive);
                            renderEditableOverlay(local);
                        }});
                    }}
                if (textAutoLabelButton) {{
                    textAutoLabelButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        autoLabelTextNodes(local);
                    }});
                }}
                if (dashedArrowButton) {{
                    dashedArrowButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        toggleSelectedEdgeDashed(local);
                    }});
                }}
                if (bothArrowButton) {{
                    bothArrowButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        toggleSelectedEdgeBothArrows(local);
                    }});
                }}
                if (lineArrowButton) {{
                    lineArrowButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        toggleSelectedEdgeLine(local);
                    }});
                }}
                if (inhibitorArrowButton) {{
                    inhibitorArrowButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        toggleSelectedEdgeEndType(local);
                    }});
                }}
                if (deleteButton) {{
                    deleteButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        deleteSelectedNode(local);
                    }});
                }}
                if (convertRectButton) {{
                    convertRectButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        convertSelectedNodesToRect(local);
                    }});
                }}
                if (batchSizeButton) {{
                    batchSizeButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || getSelectedNodeIds(local).length < 2) return;
                        setBatchSizePanelOpen(local, !local.batchSizePanelOpen);
                        if (local.batchSizePanelOpen && batchSizeValue) batchSizeValue.focus();
                    }});
                }}
                if (batchSizeCancel) {{
                    batchSizeCancel.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        setBatchSizePanelOpen(local, false);
                    }});
                }}
                if (batchSizeApply) {{
                    batchSizeApply.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        applyBatchSize(local);
                    }});
                }}
                if (sendBackButton) {{
                    sendBackButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        moveSelectedNodeToBack(local);
                    }});
                }}
                if (bringFrontButton) {{
                    bringFrontButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        moveSelectedNodeToFront(local);
                    }});
                }}
                if (undoButton) {{
                    undoButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        undoLastChange(local);
                    }});
                }}
                if (redoButton) {{
                    redoButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        redoLastChange(local);
                    }});
                }}
                if (saveButton) {{
                    saveButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !saveInputId) return;
                        if (!(window.Shiny && typeof window.Shiny.setInputValue === 'function')) return;
                        window.Shiny.setInputValue(saveInputId, {{
                            pathway_id: {_json_for_inline_script(str(info.get("id") or ""))},
                            pathway_name: {_json_for_inline_script(str(info.get("name") or ""))},
                            nodes: cloneStateSnapshot(local.editableNodes || []),
                            edges: cloneStateSnapshot(local.editableEdges || []),
                            groups: cloneStateSnapshot(local.groups || []),
                            disable_pdf_reader: !!local.disablePdfReader,
                            ts: Date.now(),
                        }}, {{ priority: 'event' }});
                    }});
                }}
                if (proteinOvalButton) {{
                    proteinOvalButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        setProteinOvalState(local, !local.proteinOvalMode);
                        renderEditableOverlay(local);
                    }});
                }}
                if (autoLabelButton) {{
                    autoLabelButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        autoLabelCustomNodes(local);
                    }});
                }}
                if (autoSizeButton) {{
                    autoSizeButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        autoSizeCustomNodes(local);
                    }});
                }}
                if (edgeResizeButton) {{
                    edgeResizeButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        setEdgeResizeState(local, !local.edgeResizeMode);
                        renderEditableOverlay(local);
                    }});
                }}
                if (opacityInput) {{
                    opacityInput.addEventListener('input', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        applyGlobalOpacity(local, opacityInput.value);
                    }});
                }}
                overlaySvg.addEventListener('contextmenu', (evt) => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (!local) return;
                    const target = evt.target;
                    const nodeId = target && target.dataset ? String(target.dataset.nodeId || '') : '';
                    const role = target && target.dataset ? String(target.dataset.role || '') : '';
                    const edgeId = target && target.dataset ? String(target.dataset.edgeId || '') : '';
                    const ptmId = target && target.dataset ? String(target.dataset.ptmId || '') : '';
                    const bendIndex = target && target.dataset ? Number(target.dataset.bendIndex || 0) : 0;
                    if (role === 'edge-bend' && edgeId) {{
                        evt.preventDefault();
                        evt.stopPropagation();
                        setSingleEdgeSelection(local, edgeId);
                        openEdgeBendMenu(local, evt, edgeId, bendIndex);
                        renderEditableOverlay(local);
                        return;
                    }}
                    if (role === 'edge-body' && edgeId) {{
                        evt.preventDefault();
                        evt.stopPropagation();
                        if (!isEntrySelected(local, 'edge', edgeId)) {{
                            const expanded = expandEntriesWithGroups(local, [{{ type: 'edge', id: edgeId }}]);
                            setSelectedEntries(local, expanded, {{ type: 'edge', id: edgeId }});
                        }}
                        showEdgeContextMenu(local, evt, edgeId);
                        renderEditableOverlay(local);
                        return;
                    }}
                    if ((role === 'ptm-body' || role === 'ptm-label') && ptmId && nodeId) {{
                        evt.preventDefault();
                        evt.stopPropagation();
                        setSinglePtmSelection(local, nodeId, ptmId);
                        showAutoPtmContextMenu(local, evt, nodeId, ptmId);
                        renderEditableOverlay(local);
                        return;
                    }}
                    if (role !== 'body' || !nodeId) return;
                    const node = findNodeById(local, nodeId);
                    if (!node) return;
                    evt.preventDefault();
                    evt.stopPropagation();
                    if (!isEntrySelected(local, 'node', nodeId)) {{
                        const expanded = expandEntriesWithGroups(local, [{{ type: 'node', id: nodeId }}]);
                        setSelectedEntries(local, expanded, {{ type: 'node', id: nodeId }});
                    }}
                    if (String(node.shapeType || 'ellipse') === 'text') {{
                        openTextMenu(local, evt, nodeId);
                        renderEditableOverlay(local);
                        return;
                    }}
                    if (node.isDrawingShape) {{
                        openShapeMenu(local, evt, nodeId);
                        renderEditableOverlay(local);
                        return;
                    }}
                    if (node.isComplex) {{
                        openComplexMenu(local, evt, nodeId);
                        renderEditableOverlay(local);
                        return;
                    }}
                    if (String(node.shapeType || 'ellipse') === 'bracket') {{
                        renderEditableOverlay(local);
                        return;
                    }}
                    showNodeContextMenu(local, evt, nodeId);
                    renderEditableOverlay(local);
                }});
                if (complexMenuButton) {{
                    complexMenuButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.contextNodeId) return;
                        const node = findNodeById(local, local.contextNodeId);
                        closeComplexMenu();
                        if (!node) return;
                        openComplexPanelForNode(local, node);
                    }});
                }}
                if (edgeBendResetButton) {{
                    edgeBendResetButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.contextEdgeId) return;
                        const edgeId = String(local.contextEdgeId || '');
                        closeEdgeBendMenu();
                        resetSelectedEdgeBend(local, edgeId);
                    }});
                }}
                if (textMenuFontSizeButton) {{
                    textMenuFontSizeButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.contextTextNodeId) return;
                        const node = findNodeById(local, local.contextTextNodeId);
                        closeTextMenu();
                        if (!node) return;
                        const currentSize = Math.max(8, Number(node.fontSize || 11));
                        const response = window.prompt('Text font size:', String(currentSize));
                        if (response === null) return;
                        const nextSize = Number(response);
                        if (!Number.isFinite(nextSize)) return;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        node.fontSize = Math.max(8, Math.min(48, nextSize));
                        const metrics = estimateTextNodeMetrics(String(node.displayLabel || node.label || 'Text'), node.fontSize || 11);
                        node.rx = metrics.rx;
                        node.ry = metrics.ry;
                        renderEditableOverlay(local);
                    }});
                }}
                if (textMenuAlignmentButton) {{
                    textMenuAlignmentButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.contextTextNodeId) return;
                        const node = findNodeById(local, local.contextTextNodeId);
                        closeTextMenu();
                        if (!node) return;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        const currentAlign = String(node.textAlign || 'center').toLowerCase();
                        node.textAlign = currentAlign === 'left' ? 'center' : (currentAlign === 'center' ? 'right' : 'left');
                        renderEditableOverlay(local);
                    }});
                }}
                if (textMenuOutlineButton) {{
                    textMenuOutlineButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.contextTextNodeId) return;
                        const node = findNodeById(local, local.contextTextNodeId);
                        closeTextMenu();
                        if (!node) return;
                        const currentColor = String(node.stroke || '#475569');
                        const response = window.prompt('Text box outline color:', currentColor);
                        if (response === null) return;
                        const nextColor = String(response || '').trim();
                        if (!nextColor) return;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        node.stroke = nextColor;
                        renderEditableOverlay(local);
                    }});
                }}
                if (textMenuBoldButton) {{
                    textMenuBoldButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.contextTextNodeId) return;
                        const node = findNodeById(local, local.contextTextNodeId);
                        closeTextMenu();
                        if (!node) return;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        node.fontWeight = String(node.fontWeight || '600') === '700' ? '500' : '700';
                        renderEditableOverlay(local);
                    }});
                }}
                if (textMenuDeleteButton) {{
                    textMenuDeleteButton.addEventListener('click', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.contextTextNodeId) return;
                        const nodeId = String(local.contextTextNodeId || '');
                        closeTextMenu();
                        setSingleSelection(local, nodeId);
                        deleteSelectedNode(local);
                    }});
                }}
                const buildCurrentExportSnapshot = () => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    return local ? buildCustomExportSnapshot(local) : {{
                        general_data: {{ settings: {{ pathway_id: {_json_for_inline_script(str(info.get("id") or ""))}, pathway_source: 'cst' }} }},
                        protbox_data: [],
                        compound_data: [],
                        text_data: [],
                        arrows: [],
                        groups: [],
                    }};
                }};
                const buildCurrentExportSvg = () => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    const sourceOverlay = document.getElementById('{overlay_id}');
                    if (!sourceOverlay) return null;
                    const exportSvg = sourceOverlay.cloneNode(true);
                    const inlineStyleProps = [
                        'display',
                        'visibility',
                        'opacity',
                        'fill',
                        'fill-opacity',
                        'fill-rule',
                        'stroke',
                        'stroke-opacity',
                        'stroke-width',
                        'stroke-dasharray',
                        'stroke-dashoffset',
                        'stroke-linecap',
                        'stroke-linejoin',
                        'stroke-miterlimit',
                        'paint-order',
                        'vector-effect',
                        'shape-rendering',
                        'text-rendering',
                        'font-family',
                        'font-size',
                        'font-weight',
                        'font-style',
                        'letter-spacing',
                        'word-spacing',
                        'text-anchor',
                        'dominant-baseline',
                        'alignment-baseline',
                        'white-space',
                        'overflow',
                        'filter',
                    ];
                    const applyInlineStyles = (sourceEl, targetEl) => {{
                        if (!sourceEl || !targetEl || sourceEl.nodeType !== 1 || targetEl.nodeType !== 1) return;
                        try {{
                            const computed = window.getComputedStyle(sourceEl);
                            inlineStyleProps.forEach((prop) => {{
                                const value = computed.getPropertyValue(prop);
                                if (value) {{
                                    targetEl.style.setProperty(prop, value);
                                }}
                            }});
                        }} catch (_err) {{}}
                        const sourceChildren = sourceEl.children || [];
                        const targetChildren = targetEl.children || [];
                        const count = Math.min(sourceChildren.length, targetChildren.length);
                        for (let idx = 0; idx < count; idx += 1) {{
                            applyInlineStyles(sourceChildren[idx], targetChildren[idx]);
                        }}
                    }};
                    applyInlineStyles(sourceOverlay, exportSvg);
                    exportSvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
                    exportSvg.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
                    exportSvg.setAttribute('xml:space', 'preserve');
                    const canvasEl = document.getElementById('{canvas_id}');
                    const contentEl = document.getElementById('{content_id}');
                    const contentRect = contentEl && typeof contentEl.getBoundingClientRect === 'function'
                        ? contentEl.getBoundingClientRect()
                        : null;
                    const exportWidth = Math.max(
                        1,
                        Number((canvasEl && canvasEl.width) || 0) ||
                        Math.round(Number((contentRect && contentRect.width) || 0)) ||
                        Math.round({page_width} * 2)
                    );
                    const exportHeight = Math.max(
                        1,
                        Number((canvasEl && canvasEl.height) || 0) ||
                        Math.round(Number((contentRect && contentRect.height) || 0)) ||
                        Math.round({page_height} * 2)
                    );
                    exportSvg.setAttribute('width', String(exportWidth));
                    exportSvg.setAttribute('height', String(exportHeight));
                    exportSvg.setAttribute('viewBox', '0 0 {page_width} {page_height}');
                    exportSvg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
                    const backgroundVisible = !local || !!local.showPdfBackground;
                    if (backgroundVisible) {{
                        if (canvasEl && typeof canvasEl.toDataURL === 'function') {{
                            try {{
                                const href = canvasEl.toDataURL('image/png');
                                if (href && href.indexOf('data:image/png') === 0) {{
                                    const imageEl = document.createElementNS('http://www.w3.org/2000/svg', 'image');
                                    imageEl.setAttribute('x', '0');
                                    imageEl.setAttribute('y', '0');
                                    imageEl.setAttribute('width', String({page_width}));
                                    imageEl.setAttribute('height', String({page_height}));
                                    imageEl.setAttributeNS('http://www.w3.org/1999/xlink', 'href', href);
                                    imageEl.setAttribute('href', href);
                                    exportSvg.insertBefore(imageEl, exportSvg.firstChild);
                                }}
                            }} catch (err) {{
                                console.warn('[MapKinase CST] export background capture failed', err);
                            }}
                        }}
                    }}
                    try {{
                        const debugSelector = 'text.cst-protein-oval-label, text.cst-auto-ptm-site-label, text.cst-auto-ptm-symbol';
                        const sourceLabels = Array.from(sourceOverlay.querySelectorAll(debugSelector));
                        const exportLabels = Array.from(exportSvg.querySelectorAll(debugSelector));
                        const n = Math.min(sourceLabels.length, exportLabels.length);
                        for (let i = 0; i < n; i += 1) {{
                            const src = sourceLabels[i];
                            const dst = exportLabels[i];
                            if (!src || !dst || typeof src.getBBox !== 'function') continue;
                            try {{
                                const b = src.getBBox();
                                if (!b) continue;
                                const bx = Number(b.x || 0);
                                const by = Number(b.y || 0);
                                const bw = Number(b.width || 0);
                                const bh = Number(b.height || 0);
                                if (![bx, by, bw, bh].every((v) => Number.isFinite(v))) continue;
                                dst.setAttribute('data-mk-bbox', `${{bx.toFixed(3)}},${{by.toFixed(3)}},${{bw.toFixed(3)}},${{bh.toFixed(3)}}`);
                            }} catch (_bboxErr) {{}}
                        }}
                    }} catch (_bboxAnnotateErr) {{}}
                    try {{
                        const serializer = new XMLSerializer();
                        return {{
                            text: serializer.serializeToString(exportSvg),
                            width: exportWidth,
                            height: exportHeight,
                        }};
                    }} catch (err) {{
                        console.error('[MapKinase CST] export svg serialize failed', err);
                        return null;
                    }}
                }};
                window.__mkExportSnapshotMap = window.__mkExportSnapshotMap || {{}};
                window.__mkExportSnapshotMap[{_json_for_inline_script(export_key)}] = buildCurrentExportSnapshot;
                window.__mkExportSvgMap = window.__mkExportSvgMap || {{}};
                window.__mkExportSvgMap[{_json_for_inline_script(export_key)}] = buildCurrentExportSvg;
                window.__mkExportSvgMap[{_json_for_inline_script(viewer_key.lower())}] = buildCurrentExportSvg;
                window.__mkExportSvgMap[{_json_for_inline_script(str(info.get("id") or "").strip().lower())}] = buildCurrentExportSvg;
                if (window.Shiny && Shiny.addCustomMessageHandler && !window.__mkExportHandlerInstalled) {{
                    window.__mkExportHandlerInstalled = true;
                    Shiny.addCustomMessageHandler('request_export_snapshot', function(msg) {{
                        const key = (msg && msg.prefix ? String(msg.prefix).toLowerCase() : '');
                        const snapshotFn = (window.__mkExportSnapshotMap && window.__mkExportSnapshotMap[key]) || buildCurrentExportSnapshot;
                        const payload = snapshotFn();
                        const prefix = msg && msg.prefix ? msg.prefix : '';
                        Shiny.setInputValue('export_snapshot', {{ prefix, payload, ts: Date.now() }}, {{ priority: 'event' }});
                    }});
                }}
                window.__mkCstViewerControls = window.__mkCstViewerControls || {{}};
                window.__mkCstViewerControls['{viewer_key}'] = {{
                    searchObjects: (query, limit = 40) => searchObjectCatalog(query, limit),
                    searchProteins: (query, limit = 40) => searchProteinCatalog(query, limit),
                    addProtbox: (uniprot) => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        const token = String(uniprot || '').trim();
                        if (!token) return false;
                        const catalog = Array.isArray(proteinSearchCatalog) ? proteinSearchCatalog : [];
                        const normalizeUni = (value) => String(value || '').trim().toUpperCase().split('-')[0];
                        const tokenUpper = token.toUpperCase();
                        const tokenUni = normalizeUni(token);
                        let entry = catalog.find((item) => normalizeUni(item && item.uniprot) === tokenUni);
                        if (!entry) {{
                            entry = catalog.find((item) => String(item && item.geneSymbol || '').trim().toUpperCase() === tokenUpper);
                        }}
                        if (!entry) {{
                            const payload = searchProteinCatalog(token, 1) || {{}};
                            const results = Array.isArray(payload.results) ? payload.results : [];
                            entry = results.length ? results[0] : null;
                        }}
                        if (!entry) return false;
                        const center = getViewportCenter();
                        const node = createProteinModuleNode(entry, center.overlayX, center.overlayY, local);
                        pushUndoSnapshot(local, buildSnapshot(local));
                        local.editableNodes = Array.isArray(local.editableNodes) ? local.editableNodes.concat([node]) : [node];
                        setSingleSelection(local, node.id);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    addCompound: (hmdbId) => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        const entry = (Array.isArray(metaboliteSearchCatalog) ? metaboliteSearchCatalog : []).find((item) => String(item.hmdbId || '') === String(hmdbId || ''));
                        if (!entry) return false;
                        const center = getViewportCenter();
                        const node = createMetaboliteNode(entry, center.overlayX, center.overlayY);
                        pushUndoSnapshot(local, buildSnapshot(local));
                        local.editableNodes = Array.isArray(local.editableNodes) ? local.editableNodes.concat([node]) : [node];
                        setSingleSelection(local, node.id);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    addArrow: (type = 'arrow') => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        const preset = String(type || 'arrow').toLowerCase();
                        const wasActive = !!local.addArrowMode && String(local.addArrowPreset || 'arrow') === preset;
                        local.addArrowPreset = preset;
                        setAddState(local, false);
                        setAddArrowState(local, !wasActive);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    addShape: (shapeType = 'ellipse') => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        const normalized = String(shapeType || 'ellipse').toLowerCase();
                        const nextType = normalized === 'square'
                            ? 'square'
                            : (normalized === 'rounded'
                                ? 'rounded'
                                : (normalized === 'bracket'
                                    ? 'bracket'
                                    : 'ellipse'));
                        const wasActive = !!local.addMode && String(local.addShapeType || 'ellipse') === nextType;
                        local.addShapeType = nextType;
                        setAddArrowState(local, false);
                        setAddState(local, !wasActive);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    addText: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        const wasActive = !!local.addMode && String(local.addShapeType || 'ellipse') === 'text';
                        local.addShapeType = 'text';
                        setAddArrowState(local, false);
                        setAddState(local, !wasActive);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    addLegend: (orientation = 'vertical') => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        const center = getViewportCenter();
                        const node = createLegendNode(orientation, center.overlayX, center.overlayY);
                        pushUndoSnapshot(local, buildSnapshot(local));
                        const firstNonShapeIndex = Array.isArray(local.editableNodes)
                            ? local.editableNodes.findIndex((item) => !(item && item.isDrawingShape))
                            : -1;
                        local.editableNodes = local.editableNodes || [];
                        if (firstNonShapeIndex >= 0) {{
                            local.editableNodes.splice(firstNonShapeIndex, 0, node);
                        }} else {{
                            local.editableNodes.push(node);
                        }}
                        setSingleSelection(local, node.id);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    undo: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        undoLastChange(local);
                        return true;
                    }},
                    redo: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        redoLastChange(local);
                        return true;
                    }},
                    deleteSelected: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        deleteSelectedNode(local);
                        return true;
                    }},
                    autoConnectEdges: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? autoConnectSelectedNodes(local, false) : 0;
                    }},
                    autoConnectShortestEdges: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? autoConnectSelectedNodes(local, true) : 0;
                    }},
                    cycleSelectedEdgeType: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? cycleSelectedEdgeInteractorType(local) : false;
                    }},
                    toggleSelectedEdgeDash: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? toggleSelectedEdgeDashed(local) : false;
                    }},
                    flipSelectedEdgeDirection: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? flipSelectedEdgeDirection(local) : false;
                    }},
                    toggleProteinOvals: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        pushUndoSnapshot(local, buildSnapshot(local));
                        setProteinOvalState(local, !local.proteinOvalMode);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    toggleShowMissing: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        setMissingState(!stage.classList.contains('cst-show-missing'));
                        emitExternalControlState(local);
                        return true;
                    }},
                    toggleEdgeResize: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return false;
                        setEdgeResizeState(local, !local.edgeResizeMode);
                        renderEditableOverlay(local);
                        return true;
                    }},
                    bringSelectedToFront: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? moveSelectedNodeToFront(local) : false;
                    }},
                    sendSelectedToBack: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? moveSelectedNodeToBack(local) : false;
                    }},
                    setMouseMode: (mode = 'drag') => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return 'drag';
                        local.mouseMode = String(mode || 'drag') === 'selection' ? 'selection' : 'drag';
                        emitExternalControlState(local);
                        return local.mouseMode;
                    }},
                    getMouseMode: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? String(local.mouseMode || 'drag') : 'drag';
                    }},
                    getState: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? buildExternalControlState(local) : {{}};
                    }},
                }};
                window.__mkViewerControls = window.__mkViewerControls || {{}};
                window.__mkViewerControls['{viewer_key}'] = {{
                    getState: () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        return local ? buildOverlayVisibilityState(local) : {{}};
                    }},
                    toggle: (key) => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return {{}};
                        const state = buildOverlayVisibilityState(local);
                        const nextValue = state[String(key || '')] === false;
                        return applyOverlayVisibilityState(local, String(key || ''), nextValue);
                    }},
                }};
                window.setTimeout(() => {{
                    try {{
                        window.dispatchEvent(new CustomEvent('mk-cst-controls-ready', {{
                            detail: {{
                                viewerKey: '{viewer_key}',
                            }},
                        }}));
                    }} catch (_err) {{}}
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (local) emitExternalControlState(local);
                }}, 0);
                if (inlineTextEditor) {{
                    inlineTextEditor.addEventListener('pointerdown', (evt) => {{
                        evt.stopPropagation();
                    }});
                    inlineTextEditor.addEventListener('keydown', (evt) => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local) return;
                        if (evt.key === 'Enter' && (evt.ctrlKey || evt.metaKey)) {{
                            finishInlineTextEdit(local, true);
                            evt.preventDefault();
                            return;
                        }}
                        if (evt.key === 'Escape') {{
                            finishInlineTextEdit(local, false);
                            evt.preventDefault();
                        }}
                    }});
                    inlineTextEditor.addEventListener('blur', () => {{
                        const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                        if (!local || !local.editingTextNodeId) return;
                        finishInlineTextEdit(local, true);
                    }});
                }}
                window.addEventListener('keydown', (evt) => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (!local) return;
                    const activeEl = document.activeElement;
                    const isTyping = !!(activeEl && /input|textarea|select/i.test(String(activeEl.tagName || '')));
                    if (evt.key === 'Tab') {{
                        if (isTyping) return;
                        if (!local.tabPreviewMode) {{
                            local.tabPreviewMode = true;
                            renderEditableOverlay(local);
                        }}
                        evt.preventDefault();
                        return;
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() === 'a') {{
                        if (isTyping) return;
                        const wasActive = !!local.addArrowMode;
                        setAddState(local, false);
                        setAddArrowState(local, !wasActive);
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() === 'q') {{
                        if (isTyping) return;
                        const wasActive = !!local.addMode && String(local.addShapeType || 'ellipse') === 'text';
                        local.addShapeType = 'text';
                        setAddArrowState(local, false);
                        setAddState(local, !wasActive);
                        evt.preventDefault();
                        return;
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() == 'z') {{
                        if (isTyping) return;
                        undoLastChange(local);
                        evt.preventDefault();
                        return;
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() === 'y') {{
                        if (isTyping) return;
                        redoLastChange(local);
                        evt.preventDefault();
                        return;
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() === 'd') {{
                        if (isTyping) return;
                        if (local.selectedEdgeId) {{
                            toggleSelectedEdgeDashed(local);
                            evt.preventDefault();
                            return;
                        }}
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() === 's') {{
                        if (isTyping) return;
                        if (local.selectedEdgeId) {{
                            cycleSelectedEdgeInteractorType(local);
                            evt.preventDefault();
                            return;
                        }}
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() === 'c') {{
                        if (isTyping) return;
                        if (local.selectedEdgeId) {{
                            if (copySelectedEdge(local)) evt.preventDefault();
                        }} else if (copySelectedNode(local)) {{
                            evt.preventDefault();
                        }}
                        return;
                    }}
                    if ((evt.ctrlKey || evt.metaKey) && !evt.shiftKey && String(evt.key || '').toLowerCase() === 'v') {{
                        if (isTyping) return;
                        if (Array.isArray(local.copiedEdges) && local.copiedEdges.length) {{
                            if (pasteCopiedEdge(local)) evt.preventDefault();
                        }} else if (pasteCopiedNode(local)) {{
                            evt.preventDefault();
                        }}
                        return;
                    }}
                    if (evt.key === 'Escape' && local.batchSizePanelOpen) {{
                        setBatchSizePanelOpen(local, false);
                        evt.preventDefault();
                        return;
                    }}
                    if (evt.key === 'Escape') {{
                        closeComplexMenu();
                        closeEdgeBendMenu();
                        closeTextMenu();
                    }}
                    const isBatchSizeField = activeEl && (activeEl === batchSizeValue || activeEl === batchSizeMode);
                    if (evt.key === 'Enter' && local.batchSizePanelOpen && (!isTyping || isBatchSizeField)) {{
                        if (applyBatchSize(local)) evt.preventDefault();
                        return;
                    }}
                    if (!local.selectedNodeId && !local.selectedEdgeId && !local.selectedPtmId) return;
                    if ((evt.key === 'Delete' || evt.key === 'Backspace') && !evt.metaKey && !evt.ctrlKey && !evt.altKey) {{
                        if (isTyping) return;
                        deleteSelectedNode(local);
                        evt.preventDefault();
                    }}
                }});
                window.addEventListener('keyup', (evt) => {{
                    const local = window[renderStateKey] && window[renderStateKey]['{viewer_key}'];
                    if (!local) return;
                    if (evt.key === 'Tab' && local.tabPreviewMode) {{
                        local.tabPreviewMode = false;
                        renderEditableOverlay(local);
                        evt.preventDefault();
                        return;
                    }}
                }});
                waitForPdfJs(40);
            }})();
            """
        ),
    )
