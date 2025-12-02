import atexit
import base64
import io
import json
import os
import re
import threading
import time
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple


from m4_json import DEFAULT_SETTINGS, get_default_json
from factory import get_pathway_api
from m2_protein_catalog import ensure_global_protein_catalog

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

from MapKinase_WebApp.m3_shiny_shortened_prevday import create_pathway_svg

try:
    import uvicorn  # type: ignore
except ImportError:  # pragma: no cover - runtime guard
    uvicorn = None

try:
    import webview  # type: ignore
except ImportError:  # pragma: no cover - runtime guard
    webview = None

HOST = os.environ.get("M5_HOST", "127.0.0.1")
PORT = int(os.environ.get("M5_PORT", os.environ.get("PORT", 8004)))
STATUS_READY = "Ready"
DISPLAY_TYPE_CHOICES = sorted({
    "prot_box",
    "gene",
    "geneproduct",
    "compound",
    "group",
    "map",
    *DEFAULT_SETTINGS.get("display_types", [])
})

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
KEGG_PATHWAYS_FILE = os.path.join(BASE_DIR, "kegg_pathways.txt")
KEGG_PATHWAY_MAX_MATCHES = 12
DEFAULT_BG_OPACITY = 0.9
DEFAULT_BG_SCALE = 1.0
DEFAULT_BG_OFFSET_X = 0.0
DEFAULT_BG_OFFSET_Y = 0.0
DEFAULT_BOX_Y_STRETCH = 1.0
JSON_PREVIEW_DIR = os.path.join(BASE_DIR, "JSONfiles")
JSON_PREVIEW_FILE = os.path.join(JSON_PREVIEW_DIR, "latest_preview.json")
CUSTOM_LAYOUT_EXPORT_FILE = os.path.join(JSON_PREVIEW_DIR, "custom_pathway_export.json")
GLOBAL_CATALOG_INFO = ensure_global_protein_catalog()
# Flip to True to mirror terminal stdout/stderr into TERMINAL_LOG_FILE by default.
TERMINAL_LOG_DEFAULT = False
TERMINAL_LOG_FILE = os.environ.get(
    "M5_TERMINAL_LOG_FILE", os.path.join(BASE_DIR, "m5_terminal_output.txt")
)

SPECIES_CHOICES: Dict[str, Dict[str, str]] = {
    "human": {"label": "Human", "code": "hsa"},
}
DEFAULT_SPECIES = "human"


def _env_var_truthy(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _load_kegg_pathways(path: str) -> List[Dict[str, str]]:
    options: List[Dict[str, str]] = []
    if not os.path.exists(path):
        return options
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped.startswith("Pathway_ID"):
                    continue
                parts = stripped.split("\t")
                if len(parts) < 2:
                    continue
                raw_id, name = parts[0].strip(), parts[1].strip()
                if not raw_id:
                    continue
                digits = raw_id[3:] if raw_id.lower().startswith("map") and len(raw_id) > 3 else raw_id
                options.append(
                    {
                        "raw_id": raw_id,
                        "digits": digits,
                        "name": name,
                    }
                )
    except Exception:
        return []
    return options


KEGG_PATHWAY_OPTIONS = _load_kegg_pathways(KEGG_PATHWAYS_FILE)

TERMINAL_LOG_ENABLED = _env_var_truthy("M5_TERMINAL_LOG", TERMINAL_LOG_DEFAULT)
CUSTOM_LAYOUT_SCHEMA_VERSION = 1


def _coerce_float(value: Any) -> Optional[float]:
    if value in (None, "", False):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _build_custom_layout_export(payload: Dict[str, Any]) -> Dict[str, Any]:
    settings = payload.get("general_data", {}).get("settings", {})
    export_data: Dict[str, Any] = {
        "schema_version": CUSTOM_LAYOUT_SCHEMA_VERSION,
        "pathway_id": settings.get("pathway_id"),
        "pathway_source": settings.get("pathway_source"),
        "protbox_data": [],
        "compound_data": [],
        "text_data": [],
        "arrows": [],
    }

    def _append_shape_section(items: Sequence[Dict[str, Any]], target: List[Dict[str, Any]], id_key: str) -> None:
        for item in items or []:
            entry_id = str(item.get(id_key) or "").strip()
            if not entry_id:
                continue
            x_val = _coerce_float(item.get("x"))
            y_val = _coerce_float(item.get("y"))
            if x_val is None or y_val is None:
                continue
            entry = {id_key: entry_id, "x": x_val, "y": y_val}
            width = _coerce_float(item.get("width"))
            height = _coerce_float(item.get("height"))
            if width is not None:
                entry["width"] = width
            if height is not None:
                entry["height"] = height
            target.append(entry)

    _append_shape_section(payload.get("protbox_data", []), export_data["protbox_data"], "protbox_id")
    _append_shape_section(payload.get("compound_data", []), export_data["compound_data"], "compound_id")
    _append_shape_section(payload.get("text_data", []), export_data["text_data"], "text_id")

    for arrow in payload.get("arrows", []) or []:
        first = str(arrow.get("protbox_id_1") or "").strip()
        second = str(arrow.get("protbox_id_2") or "").strip()
        if not first or not second:
            continue
        arrow_entry: Dict[str, Any] = {
            "protbox_id_1": first,
            "protbox_id_2": second,
        }
        for key in ("protbox_id_1_side", "protbox_id_2_side", "line", "type"):
            if key in arrow and arrow[key] not in (None, ""):
                arrow_entry[key] = arrow[key]
        export_data["arrows"].append(arrow_entry)
    return export_data


def _sanitize_custom_layout(raw_data: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(raw_data, dict):
        raise ValueError("Custom pathway file must be a JSON object.")
    sanitized: Dict[str, Any] = {
        "schema_version": int(raw_data.get("schema_version") or CUSTOM_LAYOUT_SCHEMA_VERSION),
        "pathway_id": str(raw_data.get("pathway_id") or ""),
        "pathway_source": str(raw_data.get("pathway_source") or ""),
        "protbox_data": [],
        "compound_data": [],
        "text_data": [],
        "arrows": [],
    }

    def _ingest_shape_section(section_name: str, id_key: str) -> None:
        for item in raw_data.get(section_name, []) or []:
            if not isinstance(item, dict):
                continue
            entry_id = str(item.get(id_key) or item.get("id") or item.get("entry_id") or "").strip()
            if not entry_id:
                continue
            x_val = _coerce_float(item.get("x"))
            y_val = _coerce_float(item.get("y"))
            if x_val is None or y_val is None:
                continue
            entry = {id_key: entry_id, "x": x_val, "y": y_val}
            width = _coerce_float(item.get("width"))
            height = _coerce_float(item.get("height"))
            if width is not None:
                entry["width"] = width
            if height is not None:
                entry["height"] = height
            sanitized[section_name].append(entry)

    _ingest_shape_section("protbox_data", "protbox_id")
    _ingest_shape_section("compound_data", "compound_id")
    _ingest_shape_section("text_data", "text_id")

    for arrow in raw_data.get("arrows", []) or []:
        if not isinstance(arrow, dict):
            continue
        first = str(arrow.get("protbox_id_1") or arrow.get("entry1") or "").strip()
        second = str(arrow.get("protbox_id_2") or arrow.get("entry2") or "").strip()
        if not first or not second:
            continue
        arrow_entry: Dict[str, Any] = {
            "protbox_id_1": first,
            "protbox_id_2": second,
        }
        for key in ("protbox_id_1_side", "protbox_id_2_side", "line", "type"):
            value = arrow.get(key)
            if value is None:
                continue
            arrow_entry[key] = str(value)
        sanitized["arrows"].append(arrow_entry)
    return sanitized


def _apply_custom_layout(payload: Dict[str, Any], layout: Optional[Dict[str, Any]]) -> None:
    if not payload or not layout:
        return

    def _apply_section(section_name: str, id_key: str) -> None:
        overrides = {entry[id_key]: entry for entry in layout.get(section_name, []) if entry.get(id_key)}
        if not overrides:
            return
        for item in payload.get(section_name, []):
            entry_id = str(item.get(id_key) or "").strip()
            if not entry_id or entry_id not in overrides:
                continue
            override = overrides[entry_id]
            for key in ("x", "y", "width", "height"):
                if key in override:
                    item[key] = override[key]

    _apply_section("protbox_data", "protbox_id")
    _apply_section("compound_data", "compound_id")
    _apply_section("text_data", "text_id")

    if layout.get("arrows") is not None:
        payload["arrows"] = [dict(arrow) for arrow in layout["arrows"]]


def _enable_terminal_logging(path: str) -> None:
    if not path:
        return
    log_dir = os.path.dirname(path)
    if log_dir:
        try:
            os.makedirs(log_dir, exist_ok=True)
        except OSError:
            # Directory might be read-only; skip enabling logging in that case.
            return
    try:
        log_handle = open(path, "w", encoding="utf-8")
    except OSError:
        return

    class _TeeStream:
        def __init__(self, original, log_stream):
            self._original = original
            self._log_stream = log_stream

        def write(self, data: str):
            self._log_stream.write(data)
            self._log_stream.flush()
            return self._original.write(data)

        def flush(self):
            self._log_stream.flush()
            return self._original.flush()

        def isatty(self):
            # Uvicorn logging checks this flag to decide whether to use colors.
            original_isatty = getattr(self._original, "isatty", None)
            if callable(original_isatty):
                return bool(original_isatty())
            return False

        def __getattr__(self, name):
            return getattr(self._original, name)

    atexit.register(log_handle.close)
    sys.stdout = _TeeStream(sys.stdout, log_handle)
    sys.stderr = _TeeStream(sys.stderr, log_handle)


if TERMINAL_LOG_ENABLED:
    _enable_terminal_logging(TERMINAL_LOG_FILE)


def _attach_kegg_background_image(data: Any, force: bool = False) -> Tuple[Any, bool]:
    if not isinstance(data, dict):
        return data, False

    if not force and data.get("kegg_bg_image"):
        return data, False

    settings = data.get("general_data", {}).get("settings", {})
    pathway_source = str(settings.get("pathway_source", "")).lower()
    pathway_id = settings.get("pathway_id") or data.get("pathway_id")
    if pathway_source != "kegg" or not pathway_id:
        return data, False

    try:
        api = get_pathway_api("kegg")
        img = api.download_pathway_image(pathway_id)
        if img is None:
            return data, False
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("ascii")
        updated = dict(data)
        updated["kegg_bg_image"] = f"data:image/png;base64,{b64}"
        updated["kegg_bg_size"] = {"width": img.width, "height": img.height}
        return updated, True
    except Exception:
        return data, False


def _rgb_tuple_to_hex(value: Sequence[int]) -> str:
    if isinstance(value, (list, tuple)) and len(value) == 3:
        try:
            r, g, b = (int(max(0, min(255, c))) for c in value)
            return f"#{r:02X}{g:02X}{b:02X}"
        except Exception:
            pass
    return "#000000"


def _rgb_tuple_to_list(value: Sequence[int]) -> List[int]:
    if isinstance(value, (list, tuple)) and len(value) == 3:
        return [int(max(0, min(255, c))) for c in value]
    return [0, 0, 0]


def _coerce_float(value: Any) -> Optional[float]:
    if value in (None, "", False):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _compute_gradient_color(value: Optional[float], neg: Sequence[int], pos: Sequence[int], max_neg: float, max_pos: float) -> Optional[List[int]]:
    if value is None:
        return None
    white = (255, 255, 255)
    neg = _rgb_tuple_to_list(neg)
    pos = _rgb_tuple_to_list(pos)
    if value < 0:
        denom = abs(max_neg) if max_neg else 1.0
        t = min(abs(value) / denom, 1.0)
        r = int((1 - t) * white[0] + t * neg[0])
        g = int((1 - t) * white[1] + t * neg[1])
        b = int((1 - t) * white[2] + t * neg[2])
    else:
        denom = max_pos if max_pos else 1.0
        if denom <= 0:
            denom = 1.0
        t = min(value / denom, 1.0)
        r = int((1 - t) * white[0] + t * pos[0])
        g = int((1 - t) * white[1] + t * pos[1])
        b = int((1 - t) * white[2] + t * pos[2])
    return [r, g, b]


def _apply_color_overrides(payload: Dict[str, Any], color_override: Dict[str, Any]) -> None:
    protein_data = payload.get("protein_data") or {}
    if not protein_data:
        return
    neg = color_override.get("negative_color")
    pos = color_override.get("positive_color")
    max_neg = _coerce_float(color_override.get("max_negative")) or DEFAULT_SETTINGS["max_negative"]
    max_pos = _coerce_float(color_override.get("max_positive")) or DEFAULT_SETTINGS["max_positive"]
    for protein in protein_data.values():
        _recolor_entry(protein, neg, pos, max_neg, max_pos)
        ptms = protein.get("PTMs") or {}
        for ptm in ptms.values():
            _recolor_entry(ptm, neg, pos, max_neg, max_pos)


def _recolor_entry(entry: Dict[str, Any], neg: Sequence[int], pos: Sequence[int], max_neg: float, max_pos: float) -> None:
    if not isinstance(entry, dict):
        return
    for key, value in list(entry.items()):
        if not key.startswith("fold_change_"):
            continue
        try:
            idx = key.rsplit("_", 1)[1]
        except IndexError:
            continue
        fc_value = _coerce_float(value)
        new_color = _compute_gradient_color(fc_value, neg, pos, max_neg, max_pos)
        if new_color is not None:
            entry[f"fc_color_{idx}"] = new_color


def _hex_to_rgb(value: str, fallback: Tuple[int, int, int]) -> Tuple[int, int, int]:
    if not value:
        return fallback
    value = value.strip().lstrip("#")
    if len(value) == 6:
        try:
            r = int(value[0:2], 16)
            g = int(value[2:4], 16)
            b = int(value[4:6], 16)
            return (r, g, b)
        except Exception:
            return fallback
    return fallback


def _color_picker_input(input_id: str, label: str, default_hex: str):
    picker_id = f"{input_id}_picker"
    hidden_input = ui.input_text(input_id, label, value=default_hex, width="100%")
    hidden_wrapper = ui.div({"style": "display:none;"}, hidden_input)
    picker = ui.tags.input(
        {
            "id": picker_id,
            "type": "color",
            "value": default_hex,
            "class": "form-control gradient-color-input",
            "style": "width:100%;height:42px;padding:0;",
        }
    )
    script = ui.tags.script(
        f"""
        (function(){{
            const picker = document.getElementById('{picker_id}');
            const hidden = document.getElementById('{input_id}');
            if (!picker || !hidden) return;
            const sync = (val) => {{
                hidden.value = val;
                hidden.dispatchEvent(new Event('input', {{ bubbles: true }}));
            }};
            picker.addEventListener('input', () => sync(picker.value || ''));
            // initialize hidden value on load
            sync(picker.value || hidden.value || '{default_hex}');
        }})();
        """
    )
    return ui.div(
        {"class": "gradient-color-field"},
        ui.tags.label({"for": picker_id}, label),
        picker,
        hidden_wrapper,
        script,
    )


def _to_int(value: Any, fallback: int) -> int:
    try:
        if value in (None, ""):
            return fallback
        return int(value)
    except (TypeError, ValueError):
        return fallback


def _to_float(value: Any, fallback: float) -> float:
    try:
        if value in (None, ""):
            return fallback
        return float(value)
    except (TypeError, ValueError):
        return fallback


def _to_bool(value: Any, fallback: bool) -> bool:
    if value in (None, ""):
        return fallback
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes", "on", "t"}:
            return True
        if lowered in {"false", "0", "no", "off", "f"}:
            return False
    return fallback


def _get_input_value(input, name: str):  # type: ignore[override]
    try:
        return getattr(input, name)()
    except Exception:
        return None


def _resolve_species(selection: Any) -> Tuple[str, Dict[str, str]]:
    key = str(selection).strip().lower() if selection else DEFAULT_SPECIES
    info = SPECIES_CHOICES.get(key, SPECIES_CHOICES[DEFAULT_SPECIES])
    return key, info


def _normalize_pathway_id(raw_value: str, species_code: str) -> str:
    if not raw_value:
        return ""
    cleaned = raw_value.strip()
    if "|" in cleaned:
        cleaned = cleaned.split("|", 1)[0].strip()
    if " " in cleaned:
        cleaned = cleaned.split(" ", 1)[0].strip()
    digit_match = re.search(r"(\d{5})", cleaned)
    if digit_match:
        return f"{species_code}{digit_match.group(1)}"
    try:
        pattern = re.compile(cleaned, re.IGNORECASE)
    except re.error:
        return cleaned
    for opt in KEGG_PATHWAY_OPTIONS:
        species_id = f"{species_code}{opt['digits']}"
        search_target = f"{species_id} | {opt['raw_id']} | {opt['name']}"
        if pattern.search(search_target):
            return species_id
    return cleaned

def collect_settings(input) -> Dict[str, Any]:  # type: ignore[override]
    def _get(name: str, fallback: Any) -> Any:
        value = _get_input_value(input, name)
        if value in (None, ""):
            return fallback
        return value

    species_key, species_info = _resolve_species(DEFAULT_SPECIES)
    species_code = species_info["code"]

    raw_pathway_input = str(_get("pathway_id", DEFAULT_SETTINGS["pathway_id"]))
    pathway_id = _normalize_pathway_id(raw_pathway_input, species_code) or DEFAULT_SETTINGS["pathway_id"]

    overrides: Dict[str, Any] = {}
    overrides["pathway_id"] = pathway_id
    overrides["pathway_source"] = DEFAULT_SETTINGS["pathway_source"]
    overrides["protein_selection_option"] = DEFAULT_SETTINGS["protein_selection_option"]
    overrides["ptm_selection_option"] = DEFAULT_SETTINGS["ptm_selection_option"]
    overrides["ptm_max_display"] = _to_int(_get("ptm_max_display", DEFAULT_SETTINGS["ptm_max_display"]), DEFAULT_SETTINGS["ptm_max_display"])
    overrides["show_background_image"] = _to_bool(_get("show_background_image", DEFAULT_SETTINGS["show_background_image"]), DEFAULT_SETTINGS["show_background_image"])
    overrides["display_types"] = list(DEFAULT_SETTINGS.get("display_types", []))
    overrides["show_groups"] = _to_bool(_get("show_groups", DEFAULT_SETTINGS["show_groups"]), DEFAULT_SETTINGS["show_groups"])
    overrides["show_multi_protein_indicator"] = _to_bool(_get("show_multi_protein_indicator", DEFAULT_SETTINGS["show_multi_protein_indicator"]), DEFAULT_SETTINGS["show_multi_protein_indicator"])
    overrides["show_arrows"] = _to_bool(_get("show_arrows", DEFAULT_SETTINGS.get("show_arrows", True)), DEFAULT_SETTINGS.get("show_arrows", True))
    overrides["show_text_boxes"] = _to_bool(_get("show_text_boxes", DEFAULT_SETTINGS.get("show_text_boxes", True)), DEFAULT_SETTINGS.get("show_text_boxes", True))
    overrides["negative_color"] = _hex_to_rgb(_get("negative_color", _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"])), DEFAULT_SETTINGS["negative_color"])
    overrides["positive_color"] = _hex_to_rgb(_get("positive_color", _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"])), DEFAULT_SETTINGS["positive_color"])
    overrides["max_negative"] = _to_float(_get("max_negative", DEFAULT_SETTINGS["max_negative"]), DEFAULT_SETTINGS["max_negative"])
    overrides["max_positive"] = _to_float(_get("max_positive", DEFAULT_SETTINGS["max_positive"]), DEFAULT_SETTINGS["max_positive"])
    overrides["prot_label_font"] = DEFAULT_SETTINGS["prot_label_font"]
    overrides["prot_label_size"] = _to_int(_get("prot_label_size", DEFAULT_SETTINGS["prot_label_size"]), DEFAULT_SETTINGS["prot_label_size"])
    overrides["ptm_label_font"] = DEFAULT_SETTINGS["ptm_label_font"]
    overrides["ptm_label_color"] = DEFAULT_SETTINGS["ptm_label_color"]
    overrides["ptm_label_size"] = DEFAULT_SETTINGS["ptm_label_size"]
    overrides["ptm_circle_radius"] = DEFAULT_SETTINGS["ptm_circle_radius"]
    overrides["ptm_circle_spacing"] = DEFAULT_SETTINGS["ptm_circle_spacing"]
    overrides["protein_tooltip_columns"] = list(DEFAULT_SETTINGS["protein_tooltip_columns"])
    overrides["species"] = species_key
    overrides["species_code"] = species_code
    return overrides


def collect_data_override(input) -> Dict[str, Any]:  # type: ignore[override]
    return {}


CUSTOM_STYLES = ui.tags.style(
    """
    .pathway-search-wrapper { position: relative; margin-bottom: 0.25rem; }
    .pathway-search-wrapper input.form-control { padding-right: 2.5rem; border-radius: 999px; }
    .pathway-results { position: absolute; top: calc(100% + 4px); left: 0; right: 0; max-height: 220px; overflow-y: auto;
        background: #ffffff; border: 1px solid #cbd3dd; border-radius: 12px; list-style: none; margin: 0;
        padding: 0; box-shadow: 0 18px 32px rgba(17, 24, 39, 0.18); z-index: 1030; }
    .pathway-result-item { padding: 0.55rem 0.9rem; cursor: pointer; border-top: 1px solid #eef1f5; font-size: 0.92rem; }
    .pathway-result-item:first-child { border-top: none; border-radius: 12px 12px 0 0; }
    .pathway-result-item:last-child { border-radius: 0 0 12px 12px; }
    .pathway-result-item:hover { background-color: #eef5ff; }
    .pathway-result-item.active { background-color: #dbe8ff; font-weight: 600; }
    .pathway-search-empty, .pathway-search-error { margin-top: 0.35rem; font-size: 0.85rem; color: #6b7280; padding-left: 0.2rem; }
    .pathway-search-error { color: #b91c1c; }
    .gradient-form label { font-weight: 600; font-size: 0.85rem; margin-bottom: 0.3rem; display: block; }
    .gradient-form input[type="color"] { width: 100%; height: 42px; padding: 0; border: 1px solid #cbd3dd;
        border-radius: 12px; background: none; box-shadow: inset 0 0 6px rgba(17, 24, 39, 0.08); }
    .gradient-preview { margin-top: 1rem; border-radius: 12px; border: 1px solid #cbd3dd; height: 56px;
        display: flex; align-items: center; justify-content: space-between; padding: 0 0.85rem; color: #111827;
        font-weight: 600; letter-spacing: 0.02em; }
    .gradient-preview span { background-color: rgba(255, 255, 255, 0.8); padding: 0.2rem 0.45rem;
        border-radius: 6px; font-size: 0.85rem; }
    .mode-controls { margin-top: 0.75rem; padding-top: 0.75rem; padding-bottom: 0.75rem;
        border-top: 1px solid #e5e7eb; border-bottom: 1px solid #e5e7eb; margin-bottom: 0.75rem; }
    .mode-controls .shiny-input-container { margin-bottom: 0; }
    .mode-controls label.control-label { font-weight: 700; font-size: 0.92rem; color: #1f2937; }
    .mode-controls .form-check { display: flex; align-items: center; gap: 0.5rem; margin-bottom: 0.35rem;
        padding: 0.35rem 0.55rem; border-radius: 10px; border: 1px solid transparent; transition: background 0.15s ease, border 0.15s ease; }
    .mode-controls .form-check:hover { background-color: #f3f4f6; border-color: #d1d5db; }
    .mode-controls .form-check-input { width: 1.1rem; height: 1.1rem; border-radius: 999px; }
    .mode-controls .form-check-input:checked { background-color: #2563eb; border-color: #2563eb; }
    .mode-control-disabled { opacity: 0.5; }
    """
)

MODE_BEHAVIOR_SCRIPT = ui.tags.script(
    """
    (function(){
        const MODE_INPUT_ID = "mode_selector";
        const targetIds = [
            "show_background_image",
            "show_groups",
            "show_multi_protein_indicator",
            "show_arrows",
            "show_text_boxes"
        ];
        const CUSTOM_DEFAULTS = {
            show_background_image: false,
            show_groups: false,
            show_multi_protein_indicator: false,
            show_arrows: true,
            show_text_boxes: true
        };
        const MODE_PRESETS = {
            analysis: {
                show_background_image: true,
                show_groups: true,
                show_multi_protein_indicator: true,
                show_arrows: false,
                show_text_boxes: false
            },
            figure: {
                show_background_image: false,
                show_groups: false,
                show_multi_protein_indicator: false,
                show_arrows: true,
                show_text_boxes: true
            }
        };
        let batchDepth = 0;
        const pendingValues = [];
        function pushValue(id, value){
            if (window.Shiny && typeof window.Shiny.setInputValue === "function") {
                window.Shiny.setInputValue(id, value, {priority: "event"});
            } else {
                pendingValues.push([id, value]);
            }
        }
        function beginBatch(){
            batchDepth += 1;
            if (batchDepth === 1){
                pushValue("mode_batch_active", true);
            }
        }
        function endBatch(){
            if (batchDepth === 0){
                return;
            }
            batchDepth -= 1;
            if (batchDepth === 0){
                pushValue("mode_batch_active", false);
                pushValue("mode_batch_flush", Date.now());
            }
        }
        function flushPending(){
            if (!(window.Shiny && typeof window.Shiny.setInputValue === "function")) {
                return;
            }
            while (pendingValues.length){
                const [id, value] = pendingValues.shift();
                window.Shiny.setInputValue(id, value, {priority: "event"});
            }
        }
        document.addEventListener("shiny:connected", flushPending);

        function ready(fn){
            if (document.readyState === "complete" || document.readyState === "interactive") {
                fn();
            } else {
                document.addEventListener("DOMContentLoaded", fn);
            }
        }

        ready(() => {
            const radioWrapper = document.getElementById(MODE_INPUT_ID);
            if (!radioWrapper) {
                return;
            }
            const radioInputs = Array.from(
                radioWrapper.querySelectorAll('input[type="radio"]')
            );
            if (!radioInputs.length) {
                return;
            }
            const targetEls = targetIds
                .map((id) => document.getElementById(id))
                .filter((el) => el);
            if (targetEls.length === 0) {
                return;
            }

            function setTargets(values, lock){
                targetEls.forEach((el) => {
                    const desired = Object.prototype.hasOwnProperty.call(values, el.id) ? values[el.id] : undefined;
                    if (typeof desired === "boolean" && el.checked !== desired){
                        el.checked = desired;
                        pushValue(el.id, desired);
                    }
                    el.disabled = !!lock;
                    const wrapper = el.closest(".form-check");
                    if (wrapper){
                        wrapper.classList.toggle("mode-control-disabled", !!lock);
                    }
                });
            }

            function collectCurrentValues(){
                const snapshot = {};
                targetEls.forEach((el) => {
                    snapshot[el.id] = !!el.checked;
                });
                return snapshot;
            }

            function applyValues(values, lock){
                beginBatch();
                setTargets(values, lock);
                endBatch();
            }

            function getSelectedMode(){
                const active = radioInputs.find((input) => input.checked);
                return active ? active.value : null;
            }

            function sanitizeMode(candidate){
                if (!candidate) {
                    return "custom";
                }
                if (candidate === "custom") {
                    return "custom";
                }
                return Object.prototype.hasOwnProperty.call(MODE_PRESETS, candidate) ? candidate : "custom";
            }

            function syncModeInput(mode){
                radioInputs.forEach((input) => {
                    const shouldCheck = input.value === mode;
                    if (input.checked !== shouldCheck){
                        input.checked = shouldCheck;
                    }
                });
                pushValue(MODE_INPUT_ID, mode);
            }

            let currentMode = sanitizeMode(getSelectedMode());
            let customValues = Object.assign({}, CUSTOM_DEFAULTS);

            function enterCustom(useDefaults){
                currentMode = "custom";
                if (useDefaults) {
                    customValues = Object.assign({}, CUSTOM_DEFAULTS);
                }
                applyValues(customValues, false);
                syncModeInput("custom");
            }

            function applyPreset(mode){
                if (!Object.prototype.hasOwnProperty.call(MODE_PRESETS, mode)){
                    return;
                }
                currentMode = mode;
                applyValues(MODE_PRESETS[mode], true);
                syncModeInput(mode);
            }

            function handleModeSelection(mode, fromUser){
                if (mode === "custom"){
                    if (currentMode !== "custom"){
                        customValues = collectCurrentValues();
                    }
                    enterCustom(!fromUser);
                    return;
                }
                if (!Object.prototype.hasOwnProperty.call(MODE_PRESETS, mode)){
                    return;
                }
                if (currentMode === "custom"){
                    customValues = collectCurrentValues();
                }
                applyPreset(mode);
            }

            radioInputs.forEach((input) => {
                input.addEventListener("change", () => {
                    if (!input.checked){
                        return;
                    }
                    handleModeSelection(input.value, true);
                });
            });

            targetEls.forEach((el) => {
                el.addEventListener("change", () => {
                    if (currentMode === "custom"){
                        customValues[el.id] = !!el.checked;
                    }
                });
            });

            customValues = collectCurrentValues();
            if (currentMode === "custom"){
                enterCustom(true);
            } else {
                applyPreset(currentMode);
            }
            pushValue("mode_batch_active", false);
            flushPending();
        });
    })();
    """
)

pathway_search_ui = ui.div(
    {"class": "pathway-search-wrapper"},
    ui.input_text(
        "pathway_id",
        "Pathway ID",
        value=DEFAULT_SETTINGS["pathway_id"],
        placeholder="Type regex to search KEGG pathways..."
    ),
    ui.output_ui("pathway_search_results"),
)


gradient_controls = ui.div(
    {"class": "gradient-form"},
    ui.layout_columns(
        ui.column(
            6,
            ui.input_numeric(
                "max_negative",
                "Max negative",
                value=DEFAULT_SETTINGS["max_negative"],
                step=0.1
            ),
            _color_picker_input(
                "negative_color",
                "Negative color",
                _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]),
            ),
        ),
        ui.column(
            6,
            ui.input_numeric(
                "max_positive",
                "Max positive",
                value=DEFAULT_SETTINGS["max_positive"],
                step=0.1
            ),
            _color_picker_input(
                "positive_color",
                "Positive color",
                _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]),
            ),
        ),
    ),
    ui.output_ui("gradient_preview"),
)

mode_controls = ui.div(
    {"class": "mode-controls"},
    ui.input_radio_buttons(
        "mode_selector",
        "Mode",
        choices={
            "analysis": "Data Analysis",
            "figure": "Figure Creation",
            "custom": "Custom",
        },
        selected="analysis",
    ),
)


settings_panel = ui.card(
    ui.card_header("Pathway Settings"),
    pathway_search_ui,
    ui.input_action_button("generate", "Generate Pathway", width="100%"),
    gradient_controls,
    ui.download_button("download_json", "Download JSON", width="100%"),
    ui.download_button("download_custom_pathway", "Export Custom Pathway", width="100%"),
    ui.input_action_button("save_custom_pathway_disk", "Save Custom Pathway to Disk", width="100%"),
    ui.input_file("upload_custom_pathway", "Import Custom Pathway", accept=[".json"], multiple=False),
    mode_controls,
    ui.input_checkbox("show_arrows", "Show arrows", value=True),
    ui.input_checkbox("show_text_boxes", "Show text boxes", value=True),
    ui.input_checkbox("show_background_image", "Background Image", value=False),
    ui.input_checkbox("show_multi_protein_indicator", "Show multi-protein indicator", value=False),
    ui.input_checkbox("show_groups", "Show groups", value=False),
    ui.input_numeric("prot_label_size", "Protein label size", value=DEFAULT_SETTINGS["prot_label_size"], min=1, max=72, step=1),
    ui.input_numeric("ptm_max_display", "PTM max display", value=DEFAULT_SETTINGS["ptm_max_display"], min=0),
)


preview_panel = ui.card(
    ui.card_header("Preview"),
    ui.output_ui("pathway_preview"),
    ui.hr(),
    ui.output_text("status_message"),
    ui.output_text("json_summary"),
)


app_ui = ui.page_fluid(
    CUSTOM_STYLES,
    MODE_BEHAVIOR_SCRIPT,
    ui.h2("MapKinase Settings (Shiny)"),
    ui.row(
        ui.column(4, settings_panel),
        ui.column(8, preview_panel),
    ),
)


def server(input, output, session):  # type: ignore[override]
    json_data = reactive.Value(None)
    status_msg = reactive.Value(STATUS_READY)
    custom_layout = reactive.Value(None)
    initial_build_done = reactive.Value(False)
    _event_initialized: Dict[str, bool] = {}

    def build_json():
        with reactive.isolate():
            settings_override = collect_settings(input)
            data_override = collect_data_override(input)
        try:
            payload = get_default_json(
                data_override=data_override if data_override else None,
                settings_override=settings_override,
                skip_disk_write=True,
            )
            show_bg_flag = bool(settings_override.get("show_background_image", False))
            payload, _ = _attach_kegg_background_image(payload, force=show_bg_flag)
            payload = dict(payload)
            payload["_global_protein_catalog"] = dict(GLOBAL_CATALOG_INFO)
            color_override = {
                "negative_color": _rgb_tuple_to_list(settings_override.get("negative_color", DEFAULT_SETTINGS["negative_color"])),
                "positive_color": _rgb_tuple_to_list(settings_override.get("positive_color", DEFAULT_SETTINGS["positive_color"])),
                "max_negative": settings_override.get("max_negative", DEFAULT_SETTINGS["max_negative"]),
                "max_positive": settings_override.get("max_positive", DEFAULT_SETTINGS["max_positive"]),
            }
            general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
            general_settings["negative_color"] = color_override["negative_color"]
            general_settings["positive_color"] = color_override["positive_color"]
            general_settings["max_negative"] = color_override["max_negative"]
            general_settings["max_positive"] = color_override["max_positive"]
            payload["_color_preview_override"] = color_override
            _apply_color_overrides(payload, color_override)
            try:
                os.makedirs(JSON_PREVIEW_DIR, exist_ok=True)
                with open(JSON_PREVIEW_FILE, "w", encoding="utf-8") as preview_fh:
                    json.dump(payload, preview_fh, indent=2)
            except OSError as write_err:
                print(f"Warning: could not write preview JSON to {JSON_PREVIEW_FILE}: {write_err}")
            payload["_kegg_preview"] = {
                "show": show_bg_flag,
                "opacity": DEFAULT_BG_OPACITY,
                "offset_x": DEFAULT_BG_OFFSET_X,
                "offset_y": DEFAULT_BG_OFFSET_Y,
                "scale": DEFAULT_BG_SCALE,
            }
            payload["_box_preview"] = {"y_stretch": DEFAULT_BOX_Y_STRETCH}
            layout_override = custom_layout.get()
            if layout_override:
                _apply_custom_layout(payload, layout_override)
                payload["_custom_layout_applied"] = True
            else:
                payload["_custom_layout_applied"] = False
            json_data.set(payload)
            status_msg.set("Pathway JSON generated successfully.")
        except Exception as exc:  # pragma: no cover - runtime error path
            status_msg.set(f"Error generating JSON: {exc}")
            json_data.set(None)

    def _ready_for_build(allow_batch: bool = False) -> bool:
        if not initial_build_done.get():
            return False
        if not allow_batch and _mode_batch_active():
            return False
        return True

    def _skip_initial_event(name: str) -> bool:
        if not _event_initialized.get(name):
            _event_initialized[name] = True
            return True
        return False

    def _mode_batch_active() -> bool:
        flag = _get_input_value(input, "mode_batch_active")
        return bool(flag)

    def _build_if_not_batch():
        if not _ready_for_build():
            return
        build_json()

    @reactive.Effect
    def _initial_build():
        if initial_build_done.get():
            return
        build_json()
        initial_build_done.set(True)

    @reactive.Effect
    @reactive.event(input.generate)
    def _build_on_click():
        if _skip_initial_event("generate"):
            return
        if not _ready_for_build():
            return
        build_json()

    @reactive.Effect
    @reactive.event(input.show_background_image)
    def _show_background_changed():
        if _skip_initial_event("show_background_image"):
            return
        _build_if_not_batch()

    @reactive.Effect
    @reactive.event(input.show_arrows)
    def _show_arrows_changed():
        if _skip_initial_event("show_arrows"):
            return
        _build_if_not_batch()

    @reactive.Effect
    @reactive.event(input.show_text_boxes)
    def _show_text_boxes_changed():
        if _skip_initial_event("show_text_boxes"):
            return
        _build_if_not_batch()

    @reactive.Effect
    @reactive.event(input.prot_label_size)
    def _prot_label_size_changed():
        if _skip_initial_event("prot_label_size"):
            return
        _build_if_not_batch()

    @reactive.Effect
    @reactive.event(input.mode_batch_flush)
    def _mode_batch_flushed():
        if _skip_initial_event("mode_batch_flush"):
            return
        if not _ready_for_build(allow_batch=True):
            return
        build_json()

    @reactive.Effect
    @reactive.event(input.pathway_id_choice)
    def _apply_pathway_choice():
        if _skip_initial_event("pathway_id_choice"):
            return
        if not _ready_for_build(allow_batch=True):
            return
        value = _get_input_value(input, "pathway_id_choice")
        if value:
            session.send_input_message("pathway_id", {"value": value})
            build_json()

    @reactive.Effect
    @reactive.event(input.upload_custom_pathway)
    def _import_custom_pathway():
        upload = input.upload_custom_pathway()
        if not upload:
            return
        file_info = upload[0]
        datapath = getattr(file_info, "datapath", None)
        if datapath is None and isinstance(file_info, dict):
            datapath = file_info.get("datapath")
        if not datapath:
            status_msg.set("Could not locate the uploaded custom pathway file.")
            return
        try:
            with open(datapath, "r", encoding="utf-8") as fh:
                raw_data = json.load(fh)
            layout = _sanitize_custom_layout(raw_data)
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            status_msg.set(f"Failed to import custom pathway: {exc}")
            return
        custom_layout.set(layout)
        status_msg.set("Custom pathway imported. Regenerating preview.")
        build_json()

    @reactive.Effect
    @reactive.event(input.save_custom_pathway_disk)
    def _save_custom_pathway_disk():
        if _skip_initial_event("save_custom_pathway_disk"):
            return
        data = json_data.get()
        if not data:
            status_msg.set("No pathway data available to export.")
            return
        export_payload = _build_custom_layout_export(data)
        try:
            os.makedirs(JSON_PREVIEW_DIR, exist_ok=True)
            with open(CUSTOM_LAYOUT_EXPORT_FILE, "w", encoding="utf-8") as fh:
                json.dump(export_payload, fh, indent=2)
        except OSError as exc:
            status_msg.set(f"Failed to write custom pathway: {exc}")
            return
        status_msg.set(f"Custom pathway saved to {CUSTOM_LAYOUT_EXPORT_FILE}")

    @output
    @render.ui
    def pathway_search_results():
        if not KEGG_PATHWAY_OPTIONS:
            return ui.div({"class": "pathway-search-empty"}, "No KEGG catalogue found.")
        _, species_info = _resolve_species(DEFAULT_SPECIES)
        species_code = species_info["code"]
        query = (str(_get_input_value(input, "pathway_id")) or "").strip()
        if not query:
            return ui.div({"class": "pathway-search-empty"}, "Type a regex to search the KEGG catalogue.")
        try:
            pattern = re.compile(query, re.IGNORECASE)
        except re.error:
            return ui.div({"class": "pathway-search-error"}, "Invalid regex pattern")
        matches: List[Dict[str, str]] = []
        for opt in KEGG_PATHWAY_OPTIONS:
            species_id = f"{species_code}{opt['digits']}"
            search_target = f"{species_id} | {opt['raw_id']} | {opt['name']}"
            if pattern.search(search_target):
                matches.append({"value": species_id, "label": f"{species_id} | {opt['name']}"})
            if len(matches) >= KEGG_PATHWAY_MAX_MATCHES:
                break
        if not matches:
            return ui.div({"class": "pathway-search-empty"}, "No pathways matched your search.")
        active_value = (str(_get_input_value(input, "pathway_id")) or "").strip()
        items = []
        for opt in matches:
            classes = ["pathway-result-item"]
            if active_value and active_value == opt["value"]:
                classes.append("active")
            items.append(
                ui.tags.li(
                    {
                        "class": " ".join(classes),
                        "data-value": opt["value"],
                        "onclick": "Shiny.setInputValue('pathway_id_choice', this.dataset.value, {priority: 'event'})",
                    },
                    opt["label"],
                )
            )
        return ui.TagList(
            ui.tags.ul({"class": "pathway-results"}, *items),
            ui.tags.script(
                """
                (function(){
                    const wrapper = document.querySelector('.pathway-search-wrapper');
                    if (!wrapper) return;
                    const results = wrapper.querySelector('.pathway-results');
                    if (!results) return;
                    results.style.display = '';
                    if (wrapper.dataset.bound === '1') {
                        return;
                    }
                    const input = wrapper.querySelector('input');
                    function hideResults(){
                        if (results){
                            results.style.display = 'none';
                        }
                    }
                    function handleDocClick(ev){
                        if (!wrapper.contains(ev.target) && results){
                            hideResults();
                        }
                    }
                    function showResults(){
                        if (results){
                            results.style.display = '';
                        }
                    }
                    document.addEventListener('click', handleDocClick, true);
                    if (input){
                        input.addEventListener('focus', showResults);
                        input.addEventListener('input', showResults);
                        input.addEventListener('blur', () => {
                            setTimeout(hideResults, 120);
                        });
                    }
                    wrapper.dataset.bound = '1';
                })();
                """
            ),
        )

    @output
    @render.ui
    def gradient_preview():
        neg_hex = str(_get_input_value(input, "negative_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]))
        pos_hex = str(_get_input_value(input, "positive_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]))
        neg_val = _to_float(_get_input_value(input, "max_negative"), DEFAULT_SETTINGS["max_negative"])
        pos_val = _to_float(_get_input_value(input, "max_positive"), DEFAULT_SETTINGS["max_positive"])
        style = (
            "background: linear-gradient(90deg, "
            f"{neg_hex}, #ffffff 50%, {pos_hex});"
        )
        return ui.div(
            {"class": "gradient-preview", "style": style},
            ui.tags.span(f"{neg_val}"),
            ui.tags.span("0"),
            ui.tags.span(f"{pos_val}"),
        )

    @output
    @render.ui
    def pathway_preview():
        data = json_data.get()
        if not data:
            return ui.div({"class": "alert alert-warning"}, "No pathway data to display.")
        settings = data.get("general_data", {}).get("settings", {})
        show_bg = bool(settings.get("show_background_image", False))
        return create_pathway_svg(data, show_kegg_bg=show_bg)

    @output
    @render.text
    def status_message() -> str:
        return status_msg.get()

    @output
    @render.text
    def json_summary() -> str:
        data = json_data.get()
        if not data:
            return "Summary: not available"
        protbox_count = len(data.get("protbox_data", []))
        group_count = len(data.get("groups", []))
        arrow_count = len(data.get("arrows", []))
        return f"Protboxes: {protbox_count} | Groups: {group_count} | Arrows: {arrow_count}"

    @output
    @render.download(filename="mapkinase_pathway.json")
    def download_json():
        data = json_data.get()
        if not data:
            yield b""
            return
        yield json.dumps(data, indent=2).encode("utf-8")

    @output
    @render.download(filename="mapkinase_custom_pathway.json")
    def download_custom_pathway():
        data = json_data.get()
        if not data:
            yield b""
            return
        export_payload = _build_custom_layout_export(data)
        yield json.dumps(export_payload, indent=2).encode("utf-8")


app = App(app_ui, server)


def _run_uvicorn_app():
    if uvicorn is None:  # pragma: no cover - runtime guard
        raise RuntimeError("uvicorn is required to run the Shiny app server")
    uvicorn.run(app, host=HOST, port=PORT, log_level="info")


def _launch_desktop_window(url: str):
    if webview is None:  # pragma: no cover - runtime guard
        print(f"pywebview is not installed. Open the Shiny app manually at {url}")
        return
    webview.create_window("MapKinase Settings", url)
    webview.start()


if __name__ == "__main__":
    server_thread = threading.Thread(target=_run_uvicorn_app, daemon=True)
    server_thread.start()
    time.sleep(1.2)
    app_url = f"http://{HOST}:{PORT}"
    try:
        _launch_desktop_window(app_url)
    except KeyboardInterrupt:
        pass
    finally:
        if webview is None:
            try:
                server_thread.join()
            except KeyboardInterrupt:
                pass
