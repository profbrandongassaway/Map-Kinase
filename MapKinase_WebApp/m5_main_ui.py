import atexit
import base64
import copy
import io
import json
import os
import re
import threading
import time
import sys
import csv
import asyncio
import html
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

from shiny import App, reactive, render, ui

from MapKinase_WebApp.m4_json import DEFAULT_DATA, DEFAULT_SETTINGS, get_default_json
from MapKinase_WebApp.a1_factory import get_pathway_api
from MapKinase_WebApp.m2_protein_catalog import ensure_global_protein_catalog
from MapKinase_WebApp.m1_file_processor import validate_protein_file, validate_ptm_file
from MapKinase_WebApp.d2_psp_regulatorysites import load_regulatory_sites, annotate_ptm_dataset
from MapKinase_WebApp.d2_psp_kinasesubstrates import load_kinase_substrate_map, annotate_ptm_dataset_with_kinases
from MapKinase_WebApp.d1_transfer_kegg_annotations import load_kegg_map, annotate_protein_with_kegg

from MapKinase_WebApp.m3_svg_viewer import create_pathway_svg, _build_blank_canvas

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


def _resolve_kegg_pathways_file() -> str:
    candidates = [os.path.join(BASE_DIR, "kegg_pathways.txt")]
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        candidates.extend([
            os.path.join(sys._MEIPASS, "MapKinase_WebApp", "kegg_pathways.txt"),
            os.path.join(sys._MEIPASS, "kegg_pathways.txt"),
        ])
    for path in candidates:
        if os.path.exists(path):
            return path
    print(f"Warning: kegg_pathways.txt not found (searched: {candidates})")
    return candidates[0]


KEGG_PATHWAYS_FILE = _resolve_kegg_pathways_file()
KEGG_PATHWAY_MAX_MATCHES = 12
WIKIPATHWAYS_MAX_MATCHES = 25
WEB_PATHWAY_SOURCES = ["WikiPathways", "KEGG"]
DEFAULT_BG_OPACITY = 0.9
DEFAULT_BG_SCALE = 1.0
DEFAULT_BG_OFFSET_X = 0.0
DEFAULT_BG_OFFSET_Y = 0.0
DEFAULT_BOX_Y_STRETCH = 1.0
JSON_PREVIEW_DIR = os.path.join(BASE_DIR, "JSONfiles")
JSON_PREVIEW_FILE = os.path.join(JSON_PREVIEW_DIR, "latest_preview.json")
CUSTOM_LAYOUT_EXPORT_FILE = os.path.join(JSON_PREVIEW_DIR, "custom_pathway_export.json")
RESOURCE_ROOT = getattr(sys, "_MEIPASS", PARENT_DIR)
SAMPLE_DATA_DIR = os.path.join(RESOURCE_ROOT, "sample_input_files")
SAMPLE_PROTEIN_FILE = os.path.join(SAMPLE_DATA_DIR, "853paz_Prot_mapk.csv")
SAMPLE_PTM_FILE = os.path.join(SAMPLE_DATA_DIR, "853paz_Phos_mapk.csv")
SPECIES_REF_PATH = os.path.join(BASE_DIR, "species_ref_list.csv")
UPLOAD_ACCEPT_TYPES = [
    ".txt",
    ".tsv",
    ".csv",
    "text/plain",
    "text/csv",
    "application/csv",
    "application/vnd.ms-excel",
]
# Flip to True to mirror terminal stdout/stderr into TERMINAL_LOG_FILE by default.
TERMINAL_LOG_DEFAULT = False
TERMINAL_LOG_FILE = os.environ.get(
    "M5_TERMINAL_LOG_FILE", os.path.join(BASE_DIR, "m5_terminal_output.txt")
)
MANUAL_BUILD_ONLY = True
debug_var = False

def _load_species_choices() -> Dict[str, Dict[str, str]]:
    choices: Dict[str, Dict[str, str]] = {}
    try:
        with open(SPECIES_REF_PATH, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                cleaned = { (k.lstrip("\ufeff") if isinstance(k, str) else k): v for k, v in row.items() }
                common = (cleaned.get("Common Name") or "").strip()
                code = (cleaned.get("Kegg Gene ID") or "").strip()
                species = (cleaned.get("Species") or "").strip()
                if not common or not code:
                    continue
                key = common.strip().lower().replace(" ", "_")
                choices[key] = {"label": common, "code": code, "species": species or common}
    except Exception as exc:
        print(f"Warning: could not load species_ref_list.csv: {exc}")
    if not choices:
        choices["human"] = {"label": "Human", "code": "hsa", "species": "Homo sapiens"}
    return choices

SPECIES_CHOICES: Dict[str, Dict[str, str]] = _load_species_choices()
DEFAULT_SPECIES = "human" if "human" in SPECIES_CHOICES else next(iter(SPECIES_CHOICES.keys()))
SPECIES_OPTIONS = {key: cfg["label"] for key, cfg in SPECIES_CHOICES.items() if cfg.get("label")}

MODE_PRESET_DEFAULTS: Dict[str, Dict[str, bool]] = {
    "analysis": {
        "show_background_image": True,
        "show_groups": True,
        "show_multi_protein_indicator": True,
        "show_arrows": True,
        "show_text_boxes": False,
    },
    "figure": {
        "show_background_image": False,
        "show_groups": False,
        "show_multi_protein_indicator": False,
        "show_arrows": True,
        "show_text_boxes": True,
    },
}

BOOKMARK_CONFIGS: List[Dict[str, Any]] = [
    {
        "key": "simple",
        "label": "Simple (KEGG)",
        "mode": "analysis",
        "show_search": True,
        "start_blank": False,
        "show_downloads": False,
        "show_custom_io": False,
        "show_toggles": False,
    },
    {
        "key": "web",
        "label": "Web Pathways",
        "mode": "figure",
        "show_search": True,
        "start_blank": False,
        "show_downloads": True,
        "show_custom_io": True,
        "show_toggles": True,
    },
    {
        "key": "figure",
        "label": "Figure Creation",
        "mode": "figure",
        "show_search": False,
        "start_blank": True,
        "show_downloads": True,
        "show_custom_io": True,
        "show_toggles": True,
    },
    {
        "key": "ks",
        "label": "Kinase Substrates",
        "mode": "figure",
        "show_search": False,
        "start_blank": True,
        "show_downloads": False,
        "show_custom_io": False,
        "show_toggles": False,
    },
]


def _prefixed_id(prefix: str, name: str) -> str:
    return f"{prefix}_{name}"


def _bool_default(cfg: Dict[str, Any], key: str) -> bool:
    return MODE_PRESET_DEFAULTS.get(cfg.get("mode"), {}).get(key, DEFAULT_SETTINGS.get(key, False))


def _env_var_truthy(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}

LAUNCH_DESKTOP_GUI = _env_var_truthy("M5_DESKTOP_GUI", default=False)

BUILD_GLOBAL_CATALOG_ON_STARTUP = _env_var_truthy(
    "M5_BUILD_GLOBAL_CATALOG_ON_STARTUP",
    default=not getattr(sys, "frozen", False),
)


def _init_global_catalog() -> Dict[str, Any]:
    default_path = os.environ.get(
        "GLOBAL_PROTEIN_CATALOG_PATH",
        os.path.join(BASE_DIR, "cache", "global_protein_catalog.json"),
    )
    if BUILD_GLOBAL_CATALOG_ON_STARTUP:
        return ensure_global_protein_catalog()
    metadata: Dict[str, Any] = {}
    if os.path.exists(default_path):
        try:
            with open(default_path, "r", encoding="utf-8") as fh:
                payload = json.load(fh)
            metadata = payload.get("metadata", {}) if isinstance(payload, dict) else {}
        except (OSError, json.JSONDecodeError):
            metadata = {}
    os.environ.setdefault("GLOBAL_PROTEIN_CATALOG_PATH", default_path)
    return {"path": default_path, "metadata": metadata}


GLOBAL_CATALOG_INFO = _init_global_catalog()


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
WIKIPATHWAYS_CACHE: Dict[str, List[Dict[str, str]]] = {}
WIKIPATHWAYS_ORG_CACHE: Optional[set] = None


def _wp_supported_orgs() -> set:
    global WIKIPATHWAYS_ORG_CACHE
    if WIKIPATHWAYS_ORG_CACHE is not None:
        return WIKIPATHWAYS_ORG_CACHE
    try:
        import pywikipathways as pwp  # type: ignore
        orgs = pwp.list_organisms()
        if isinstance(orgs, (list, tuple, set)):
            WIKIPATHWAYS_ORG_CACHE = {str(o).strip().lower() for o in orgs if o}
        else:
            WIKIPATHWAYS_ORG_CACHE = set()
    except Exception:
        WIKIPATHWAYS_ORG_CACHE = set()
    return WIKIPATHWAYS_ORG_CACHE


def _load_wikipathways_catalog(organism: str, fallback: Optional[str] = None) -> List[Dict[str, str]]:
    cache_key = organism.strip().lower() if organism else "__all__"
    if cache_key in WIKIPATHWAYS_CACHE:
        return WIKIPATHWAYS_CACHE[cache_key]
    options: List[Dict[str, str]] = []

    def _fetch_df(name: str, quiet: bool = False):
        try:
            import pywikipathways as pwp  # type: ignore
            df = pwp.list_pathways(name or "")
            # Convert list responses to DataFrame for consistent handling
            if isinstance(df, list):
                try:
                    import pandas as pd  # type: ignore
                    df = pd.DataFrame(df)
                except Exception:
                    return None
            if df is None:
                return None
            if hasattr(df, "empty") and df.empty:
                return None
            return df
        except Exception as exc:  # pragma: no cover - network/service issues
            if not quiet:
                print(f"Warning: list_pathways failed for '{name}': {exc}")
            return None

    try:
        names_to_try: List[str] = []
        if organism:
            names_to_try.append(organism)
        if fallback and (not organism or fallback.lower() != organism.lower()):
            names_to_try.append(fallback)
        # Prefer the first name that is explicitly supported
        supported = _wp_supported_orgs()
        if supported:
            for idx, nm in enumerate(names_to_try):
                if nm and nm.strip().lower() in supported:
                    # Move supported name to front
                    names_to_try.insert(0, names_to_try.pop(idx))
                    break

        df = None
        for idx, nm in enumerate(names_to_try):
            quiet = idx < len(names_to_try) - 1
            df = _fetch_df(nm, quiet=quiet)
            if df is not None:
                break
        if df is None:
            WIKIPATHWAYS_CACHE[cache_key] = []
            return []
        records = df.to_dict("records") if hasattr(df, "to_dict") else []
        for row in records:
            path_id = str(row.get("id") or row.get("wpid") or "").strip()
            name = str(row.get("name") or "").strip()
            org = str(row.get("species") or organism or fallback or "").strip()
            if not path_id or not name:
                continue
            path_id = path_id.upper()
            label = f"{path_id} | {name}"
            if org:
                label = f"{label} ({org})"
                options.append({"id": path_id, "name": name, "species": org, "label": label})
            options.sort(key=lambda opt: opt.get("name", "").lower())
    except Exception as exc:
        print(f"Warning: failed to load WikiPathways catalogue for '{organism or fallback or 'all'}': {exc}")
    WIKIPATHWAYS_CACHE[cache_key] = options
    return options

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
        "pathway_source": "custom",
        "protbox_data": [],
        "compound_data": [],
        "text_data": [],
        "arrows": [],
        "groups": [],
    }

    def _append_shape_section(items: Sequence[Dict[str, Any]], target: List[Dict[str, Any]], id_key: str, extra_keys: Optional[Sequence[str]] = None) -> None:
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
            for extra in extra_keys or []:
                if extra in item and item[extra] not in (None, ""):
                    entry[extra] = item[extra]
            # Preserve simple metadata for export
            for key, val in item.items():
                if key in entry or (extra_keys and key in extra_keys):
                    continue
                if isinstance(val, (str, int, float, list, dict)):
                    entry[key] = val
            target.append(entry)

    _append_shape_section(
        payload.get("protbox_data", []),
        export_data["protbox_data"],
        "protbox_id",
        extra_keys=("proteins", "uniprot_ids", "uniprot", "protein_ids", "label"),
    )
    _append_shape_section(payload.get("compound_data", []), export_data["compound_data"], "compound_id")
    _append_shape_section(
        payload.get("text_data", []),
        export_data["text_data"],
        "text_id",
        extra_keys=("html", "text_style", "bgcolor", "fgcolor", "border_color", "label"),
    )

    for arrow in payload.get("arrows", []) or []:
        first = str(arrow.get("protbox_id_1") or "").strip()
        second = str(arrow.get("protbox_id_2") or "").strip()
        if not first or not second:
            continue
        arrow_entry: Dict[str, Any] = {
            "protbox_id_1": first,
            "protbox_id_2": second,
        }
        for key, val in arrow.items():
            if key in arrow_entry:
                continue
            if val in (None, ""):
                continue
            if isinstance(val, (str, int, float, bool, list, dict)):
                arrow_entry[key] = val
        export_data["arrows"].append(arrow_entry)
    for group in payload.get("groups", []) or []:
        if not isinstance(group, dict):
            continue
        entry: Dict[str, Any] = {}
        for key, val in group.items():
            if val in (None, ""):
                continue
            if isinstance(val, (str, int, float, bool, list, dict)):
                entry[key] = val
        if entry:
            export_data["groups"].append(entry)
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
        "groups": [],
    }

    def _ingest_shape_section(section_name: str, id_key: str, extra_keys: Optional[Sequence[str]] = None) -> None:
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
            for extra in extra_keys or []:
                if extra in item and item[extra] not in (None, ""):
                    entry[extra] = item[extra]
            sanitized[section_name].append(entry)

    _ingest_shape_section("protbox_data", "protbox_id", extra_keys=("proteins", "uniprot_ids", "uniprot", "protein_ids", "label"))
    _ingest_shape_section("compound_data", "compound_id")
    _ingest_shape_section("text_data", "text_id", extra_keys=("html", "text_style", "bgcolor", "fgcolor", "border_color", "label"))

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
    for group in raw_data.get("groups", []) or []:
        if not isinstance(group, dict):
            continue
        entry: Dict[str, Any] = {}
        for key, val in group.items():
            if val in (None, ""):
                continue
            if isinstance(val, (str, int, float, bool, list, dict)):
                entry[key] = val
        if entry:
            sanitized["groups"].append(entry)
    return sanitized


def _apply_custom_layout(payload: Dict[str, Any], layout: Optional[Dict[str, Any]]) -> None:
    if not payload or not layout:
        return

    def _apply_section(section_name: str, id_key: str, extra_keys: Optional[Sequence[str]] = None) -> None:
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
            for extra in extra_keys or []:
                if extra in override:
                    item[extra] = override[extra]

    _apply_section("protbox_data", "protbox_id", extra_keys=("proteins", "uniprot_ids", "uniprot", "protein_ids", "label"))
    _apply_section("compound_data", "compound_id")
    _apply_section("text_data", "text_id", extra_keys=("html", "text_style", "bgcolor", "fgcolor", "border_color", "label"))

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
            const copyHiddenToPicker = () => {{
                if ((hidden.value || '') !== (picker.value || '')) {{
                    picker.value = hidden.value || '{default_hex}';
                }}
            }};
            hidden.addEventListener('input', copyHiddenToPicker);
            hidden.addEventListener('change', copyHiddenToPicker);
            // Fallback poll in case Shiny updates the value without firing input/change
            let lastVal = hidden.value;
            setInterval(() => {{
                if (hidden.value !== lastVal) {{
                    lastVal = hidden.value;
                    copyHiddenToPicker();
                }}
            }}, 300);
            // initialize both directions on load
            const startVal = hidden.value || picker.value || '{default_hex}';
            picker.value = startVal;
            sync(startVal);
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


def _send_custom_message(session, name: str, payload: Dict[str, Any]) -> None:
    try:
        result = session.send_custom_message(name, payload)
        if asyncio.iscoroutine(result):
            asyncio.create_task(result)
    except Exception:
        return


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

def collect_settings(input, cfg: Dict[str, Any]) -> Dict[str, Any]:  # type: ignore[override]
    prefix = cfg["key"]

    def _get(name: str, fallback: Any) -> Any:
        value = _get_input_value(input, _prefixed_id(prefix, name))
        if value in (None, ""):
            return fallback
        return value

    species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
    species_code = species_info["code"]
    species_label = species_info.get("label") or species_key
    species_full = species_info.get("species") or species_label

    raw_pathway_value = _get_input_value(input, _prefixed_id(prefix, "pathway_id"))
    if raw_pathway_value in (None, "") and prefix != "web":
        raw_pathway_value = DEFAULT_SETTINGS["pathway_id"]
    raw_pathway_input = str(raw_pathway_value or "")
    pathway_source = DEFAULT_SETTINGS["pathway_source"]
    if prefix == "web":
        pathway_source = str(_get_input_value(input, _prefixed_id(prefix, "pathway_source_choice")) or "wikipathways").strip().lower()
        cleaned_input = raw_pathway_input.strip()
        if pathway_source == "kegg":
            normalized = _normalize_pathway_id(cleaned_input, species_code)
            pathway_id = (normalized or cleaned_input or DEFAULT_SETTINGS["pathway_id"]).lower()
        else:
            if cleaned_input.lower().startswith("wp"):
                cleaned_input = f"WP{cleaned_input[2:]}"
            pathway_id = cleaned_input.upper()
    else:
        pathway_source = "kegg"
        pathway_id = _normalize_pathway_id(raw_pathway_input, species_code) or DEFAULT_SETTINGS["pathway_id"]

    overrides: Dict[str, Any] = {}
    overrides["pathway_id"] = pathway_id
    overrides["pathway_source"] = pathway_source
    overrides["protein_selection_option"] = DEFAULT_SETTINGS["protein_selection_option"]
    overrides["ptm_selection_option"] = DEFAULT_SETTINGS["ptm_selection_option"]
    global_ptm_max = _to_int(_get_input_value(input, "settings_ptm_max_display"), DEFAULT_SETTINGS["ptm_max_display"])
    overrides["ptm_max_display"] = global_ptm_max
    overrides["use_original_protbox_size"] = _to_bool(
        _get_input_value(input, "settings_use_original_protbox_size"),
        DEFAULT_SETTINGS.get("use_original_protbox_size", False),
    )
    overrides["show_background_image"] = _to_bool(_get("show_background_image", _bool_default(cfg, "show_background_image")), _bool_default(cfg, "show_background_image"))
    overrides["display_types"] = list(DEFAULT_SETTINGS.get("display_types", []))
    overrides["show_groups"] = _to_bool(
        _get_input_value(input, "settings_show_groups"),
        _bool_default(cfg, "show_groups"),
    )
    overrides["debug_mode"] = _to_bool(
        _get_input_value(input, "settings_debug_mode"),
        False,
    )
    overrides["show_multi_protein_indicator"] = _to_bool(
        _get_input_value(input, "settings_show_multi_protein_indicator"),
        _bool_default(cfg, "show_multi_protein_indicator"),
    )
    overrides["show_arrows"] = _to_bool(_get("show_arrows", _bool_default(cfg, "show_arrows")), _bool_default(cfg, "show_arrows"))
    overrides["show_text_boxes"] = _to_bool(_get("show_text_boxes", _bool_default(cfg, "show_text_boxes")), _bool_default(cfg, "show_text_boxes"))
    overrides["mode"] = str(cfg.get("mode", "analysis")).strip().lower()
    overrides["negative_color"] = _hex_to_rgb(
        _get_input_value(input, "settings_negative_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]),
        DEFAULT_SETTINGS["negative_color"],
    )
    overrides["positive_color"] = _hex_to_rgb(
        _get_input_value(input, "settings_positive_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]),
        DEFAULT_SETTINGS["positive_color"],
    )
    overrides["max_negative"] = _to_float(
        _get_input_value(input, "settings_max_negative"), DEFAULT_SETTINGS["max_negative"]
    )
    overrides["max_positive"] = _to_float(
        _get_input_value(input, "settings_max_positive"), DEFAULT_SETTINGS["max_positive"]
    )
    overrides["prot_label_font"] = DEFAULT_SETTINGS["prot_label_font"]
    global_prot_label = _to_int(_get_input_value(input, "settings_prot_label_size"), DEFAULT_SETTINGS["prot_label_size"])
    overrides["prot_label_size"] = global_prot_label
    overrides["ptm_label_font"] = DEFAULT_SETTINGS["ptm_label_font"]
    overrides["ptm_label_color"] = DEFAULT_SETTINGS["ptm_label_color"]
    overrides["ptm_label_size"] = DEFAULT_SETTINGS["ptm_label_size"]
    overrides["ptm_circle_radius"] = DEFAULT_SETTINGS["ptm_circle_radius"]
    overrides["ptm_circle_spacing"] = DEFAULT_SETTINGS["ptm_circle_spacing"]
    overrides["prot_outline_width"] = _to_float(
        _get_input_value(input, "settings_prot_outline_width"), DEFAULT_SETTINGS.get("prot_outline_width", 1)
    )
    overrides["ptm_outline_width"] = _to_float(
        _get_input_value(input, "settings_ptm_outline_width"), DEFAULT_SETTINGS.get("ptm_outline_width", 1)
    )
    overrides["protein_tooltip_columns"] = list(DEFAULT_SETTINGS["protein_tooltip_columns"])
    overrides["include_psp_tooltips"] = bool(_get_input_value(input, "settings_include_psp_tooltips"))
    overrides["species"] = species_label
    overrides["species_code"] = species_code
    overrides["_species_full_name"] = species_full
    # Default to KEGG gene IDs; will override to UniProt for non-KEGG later when dataset headers are known.
    overrides["hsa_id_column"] = "KEGG_Gene_ID"
    # Use KEGG column added during upload-time annotation
    overrides["hsa_id_column"] = "KEGG_Gene_ID"
    return overrides


def _empty_ks_index() -> Dict[str, Any]:
    return {
        "kinases": {},
        "substrates": {},
        "ptms_by_uniprot": {},
        "ptm_headers": [],
        "prot_gene_map": {},
    }


PTM_POSITION_PRIORITY = ["W1", "W2", "E1", "E2", "N1", "S1", "N2", "S2", "N3", "S3", ]
PTM_LABEL_DEFAULTS: Dict[str, Tuple[float, float, str]] = {
    "N1": (-5, -5, "right"), # -5, -5,
    "N2": (0, -11, "center"),
    "N3": (5, -5, "left"), # 5, -5,
    "S1": (-3, 5, "right"),
    "S2": (0, 12, "center"),
    "S3": (3, 5, "left"),
    "W1": (-5, -2, "right"),
    "W2": (-5, 2, "right"),
    "E1": (5, -2, "left"), # 3, -2
    "E2": (5, 2, "left"),   # 3, 2
}

KS_VERTICAL_SPACING = 25  


def _compute_ptm_position(
    pos_key: str, x: float, y: float, width: float, height: float, spacing: float
) -> Tuple[float, float]:
    half = spacing * 0.5
    positions = {
        "N1": (x + width * 0.2, y - spacing),
        "N2": (x + width * 0.5, y - spacing),
        "N3": (x + width * 0.8, y - spacing),
        "S1": (x + width * 0.2, y + height + spacing),
        "S2": (x + width * 0.5, y + height + spacing),
        "S3": (x + width * 0.8, y + height + spacing),
        "W1": (x - spacing, y + height * 0.33 - half),
        "W2": (x - spacing, y + height * 0.66 + half),
        "E1": (x + width + spacing, y + height * 0.33 - half),
        "E2": (x + width + spacing, y + height * 0.66 + half),
    }
    return positions.get(pos_key, (x, y))


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
    neg = [int(v) for v in negative_color][:3]
    pos = [int(v) for v in positive_color][:3]
    neg_limit = abs(max_negative) if max_negative else 1.0
    pos_limit = max_positive if max_positive else 1.0
    pos_limit = pos_limit if pos_limit != 0 else 1.0
    white = (255, 255, 255)
    if fold < 0:
        t = min(abs(fold) / neg_limit, 1.0)
        return [int((1 - t) * white[i] + t * neg[i]) for i in range(3)]
    t = min(fold / pos_limit, 1.0)
    return [int((1 - t) * white[i] + t * pos[i]) for i in range(3)]

def _normalize_fc_suffix(header: str) -> str:
    value = str(header or "").strip()
    if ":" in value:
        value = value.split(":", 1)[1]
    return re.sub(r"\s+", " ", value.strip()).lower()

def _outline_column_map(headers: Sequence[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for col in headers or []:
        col_str = str(col)
        if col_str.strip().lower().startswith("o:"):
            key = _normalize_fc_suffix(col_str)
            if key:
                mapping[key] = col_str
    return mapping

def _resolve_outline_columns(main_columns: Sequence[str], headers: Sequence[str]) -> List[Optional[str]]:
    outline_map = _outline_column_map(headers)
    return [outline_map.get(_normalize_fc_suffix(col)) for col in main_columns or []]


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
    .gradient-preview { margin-top: 1rem; border-radius: 12px; border: 1px solid #cbd3dd; height: 42px;
        display: flex; align-items: center; justify-content: space-between; padding: 0 0.85rem; color: #111827;
        font-weight: 600; letter-spacing: 0.02em; }
    .gradient-preview span { background-color: rgba(255, 255, 255, 0.8); padding: 0.2rem 0.45rem;
        border-radius: 6px; font-size: 0.85rem; }
    #settings_max_negative, #settings_max_positive { height: 48px; font-size: 1rem; }
    #settings_gradient_preview .gradient-preview { margin-top: 0.5rem; }
    .mk-export-spinner { display: inline-flex; align-items: center; gap: 6px; font-size: 12px; color: #1f2937; }
    .mk-export-spinner::before { content: ""; width: 14px; height: 14px; border: 2px solid #cbd5e1;
        border-top-color: #1f2937; border-radius: 50%; animation: mk-spin 0.8s linear infinite; }
    @keyframes mk-spin { to { transform: rotate(360deg); } }
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
    .ks-search-row { display: flex; gap: 8px; align-items: flex-end; margin-bottom: 0.5rem; }
    .ks-search-row .selectize-control { flex: 1; }
    .ks-search-row .form-control { width: 100%; }
    .ks-mode-toggle { display: flex; flex-direction: column; gap: 4px; }
    .ks-mode-toggle .form-check { display: flex; align-items: center; padding: 10px 14px; margin: 0;
        border: 1px solid #cbd3dd; border-radius: 10px; background: #f8fafc; cursor: pointer; }
    .ks-mode-toggle .form-check-input { width: 1.2rem; height: 1.2rem; margin-right: 8px; }
    .ks-mode-toggle .form-check-input:checked { background-color: #2563eb; border-color: #2563eb; }
    .ks-mode-toggle .form-check:hover { border-color: #2563eb; }
    .ks-filter-button { display: inline-flex; align-items: center; gap: 6px; padding: 8px 12px; border-radius: 10px;
        border: 1px solid #cbd3dd; background: #eef2ff; color: #1e3a8a; font-weight: 600; cursor: pointer; }
    .ks-filter-button:hover { background: #e0e7ff; }
    .ks-filter-panel { position: relative; margin-top: 8px; }
    .ks-filter-popup { position: absolute; z-index: 20; min-width: 320px; background: #fff; border: 1px solid #d4d7dd;
        box-shadow: 0 14px 32px rgba(0,0,0,0.18); border-radius: 12px; padding: 12px; display: none; }
    .ks-filter-popup.active { display: block; }
    .ks-filter-row { display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 8px; }
    .ks-filter-row > * { flex: 1; min-width: 140px; }
    .ks-filter-disabled { opacity: 0.45; pointer-events: none; }
    .pathway-table { width: 100%; }
    .pathway-viewer-card { position: relative; }
    .svg-download-row { margin-top: 10px; display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
    .svg-download-row .svg-download-label { font-weight: 700; color: #1f2937; margin-right: 4px; }
    .svg-download-row .btn { font-weight: 600; padding: 6px 12px; }
    """
)

NAV_LOCK_SCRIPT = ui.tags.script(
    """
    (function(){
        if (!(window.Shiny && Shiny.addCustomMessageHandler)) return;
        Shiny.addCustomMessageHandler("toggle_nav_lock", function(msg){
            var locked = !!(msg && msg.locked);
            var nav = document.getElementById("bookmark_selector");
            if (!nav) return;
            var links = nav.querySelectorAll("a.nav-link");
            links.forEach(function(link){
                var val = link.getAttribute("data-value");
                if (val === "input") return;
                if (locked){
                    if (!link.dataset.prevHref){
                        link.dataset.prevHref = link.getAttribute("href") || "";
                    }
                    link.removeAttribute("href");
                    link.style.pointerEvents = "none";
                    link.style.opacity = "0.55";
                    link.style.cursor = "not-allowed";
                } else {
                    if (link.dataset.prevHref){
                        link.setAttribute("href", link.dataset.prevHref);
                    }
                    link.style.pointerEvents = "";
                    link.style.opacity = "";
                    link.style.cursor = "";
                }
            });
        });
    })();
    """
)

EXPORT_DOWNLOAD_SCRIPT = ui.tags.script(
    """
    (function(){
        if (!(window.Shiny && Shiny.addCustomMessageHandler)) return;
        Shiny.addCustomMessageHandler("download_payload", function(msg){
            var filename = msg && msg.filename ? msg.filename : "custom_pathway.json";
            var content = msg && typeof msg.content === "string" ? msg.content : "";
            var blob = new Blob([content], { type: "application/json;charset=utf-8" });
            var url = URL.createObjectURL(blob);
            var link = document.createElement("a");
            link.href = url;
            link.download = filename;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            URL.revokeObjectURL(url);
            var btn = msg && msg.button_id ? document.getElementById(msg.button_id) : null;
            if (btn) btn.disabled = false;
            var spinner = msg && msg.spinner_id ? document.getElementById(msg.spinner_id) : null;
            if (spinner) spinner.style.display = "none";
        });
        Shiny.addCustomMessageHandler("export_failed", function(msg){
            var btn = msg && msg.button_id ? document.getElementById(msg.button_id) : null;
            if (btn) btn.disabled = false;
            var spinner = msg && msg.spinner_id ? document.getElementById(msg.spinner_id) : null;
            if (spinner) spinner.style.display = "none";
        });
    })();
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
                show_arrows: true,
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

SVG_DOWNLOAD_SCRIPT = ui.tags.script(
    """
    (function(){
        function downloadHref(href, filename){
            var link = document.createElement("a");
            link.href = href;
            link.download = filename;
            document.body.appendChild(link);
            link.click();
            link.remove();
        }
        function serializeSvg(svg){
            if (!svg) return null;
            var clone = svg.cloneNode(true);
            if (!clone.getAttribute("xmlns")){
                clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
            }
            var vb = clone.viewBox && clone.viewBox.baseVal;
            var rect = svg.getBoundingClientRect();
            var width = (vb && vb.width) ? vb.width : (rect.width || 1200);
            var height = (vb && vb.height) ? vb.height : (rect.height || 900);
            clone.setAttribute("width", width);
            clone.setAttribute("height", height);
            var serializer = new XMLSerializer();
            var text = serializer.serializeToString(clone);
            return { text: text, width: width, height: height };
        }
        function ensureJsPdf(){
            return new Promise(function(resolve, reject){
                if (window.jspdf && window.jspdf.jsPDF){
                    resolve(window.jspdf.jsPDF);
                    return;
                }
                var existing = document.querySelector("script[data-mk-jspdf]");
                if (existing && existing.dataset.loading === "1"){
                    existing.addEventListener("load", function(){ resolve(window.jspdf && window.jspdf.jsPDF); });
                    existing.addEventListener("error", function(){ reject(new Error("Failed to load jsPDF")); });
                    return;
                }
                var script = document.createElement("script");
                script.src = "https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js";
                script.async = true;
                script.dataset.mkJspdf = "1";
                script.dataset.loading = "1";
                script.onload = function(){ resolve(window.jspdf && window.jspdf.jsPDF); };
                script.onerror = function(){ reject(new Error("Failed to load jsPDF")); };
                document.head.appendChild(script);
            });
        }
        function findSvg(btn){
            if (!btn) return null;
            var scoped = btn.closest(".pathway-viewer-card");
            if (scoped){
                var found = scoped.querySelector(".svg-container svg");
                if (found) return found;
            }
            return document.querySelector(".svg-container svg");
        }
        function exportSvg(svgEl, format, name){
            var payload = serializeSvg(svgEl);
            if (!payload) return;
            var text = payload.text;
            var width = payload.width;
            var height = payload.height;
            if (format === "svg"){
                var blob = new Blob([text], { type: "image/svg+xml" });
                var href = URL.createObjectURL(blob);
                downloadHref(href, name + ".svg");
                setTimeout(function(){ URL.revokeObjectURL(href); }, 1500);
                return;
            }
            var canvas = document.createElement("canvas");
            canvas.width = width;
            canvas.height = height;
            var ctx = canvas.getContext("2d");
            var img = new Image();
            var svgBlob = new Blob([text], { type: "image/svg+xml" });
            var url = URL.createObjectURL(svgBlob);
            img.onload = function(){
                try {
                    ctx.fillStyle = "white";
                    ctx.fillRect(0, 0, width, height);
                    ctx.drawImage(img, 0, 0, width, height);
                    if (format === "pdf"){
                        ensureJsPdf().then(function(jsPDF){
                            if (!jsPDF){ console.error("jsPDF unavailable"); return; }
                            var orientation = width >= height ? "l" : "p";
                            var pdf = new jsPDF({ orientation: orientation, unit: "pt", format: [width, height] });
                            pdf.addImage(canvas.toDataURL("image/png"), "PNG", 0, 0, width, height);
                            pdf.save(name + ".pdf");
                        }).catch(function(err){ console.error("PDF export failed", err); });
                    } else {
                        var mime = format === "jpeg" ? "image/jpeg" : "image/png";
                        var ext = format === "jpeg" ? ".jpeg" : ".png";
                        var data = canvas.toDataURL(mime, 1.0);
                        downloadHref(data, name + ext);
                    }
                } finally {
                    URL.revokeObjectURL(url);
                }
            };
            img.onerror = function(err){
                console.error("SVG export failed", err);
                URL.revokeObjectURL(url);
            };
            img.src = url;
        }
        document.addEventListener("click", function(ev){
            var btn = ev.target.closest("[data-mk-download]");
            if (!btn) return;
            ev.preventDefault();
            var format = btn.getAttribute("data-mk-download");
            var svg = findSvg(btn);
            if (!svg){
                console.warn("No SVG element found for download");
                return;
            }
            var name = btn.getAttribute("data-mk-name") || "pathway";
            exportSvg(svg, format, name);
        });
        window.mkExportPathwaySvg = exportSvg;
    })();
    """
)

def _pathway_search_ui(prefix: str) -> ui.div:
    placeholder = "Type regex to search KEGG pathways..." if prefix != "web" else "Search pathways (regex or text)..."
    default_value = "" if prefix == "web" else DEFAULT_SETTINGS["pathway_id"]
    filter_btn_id = _prefixed_id(prefix, "pathway_filter_btn")
    filter_popup_id = _prefixed_id(prefix, "pathway_filter_popup")
    filter_choices = {opt: opt for opt in WEB_PATHWAY_SOURCES}
    children = []
    search_input = ui.input_text(
        _prefixed_id(prefix, "pathway_id"),
        "Pathway ID",
        value=default_value,
        placeholder=placeholder,
        width="100%",
    )
    if prefix == "web":
        children.append(
            ui.div(
                {"style": "display:flex; flex-direction:column; gap:8px; position:relative;", "class": "pathway-search-row"},
                ui.div({"style": "flex:1;"}, search_input),
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.tags.button(
                        ui.tags.i({"class": "fa fa-filter"}), "Filters",
                        {"id": filter_btn_id, "type": "button", "class": "ks-filter-button pathway-filter-button"}
                    ),
                ),
                ui.div(
                    {"id": filter_popup_id, "class": "ks-filter-popup"},
                    ui.div({"class": "ks-filter-row"},
                        ui.input_checkbox_group(
                            _prefixed_id(prefix, "pathway_filter_sources"),
                            "Sources",
                            choices=filter_choices,
                            selected=list(filter_choices.keys()),
                            inline=False,
                        ),
                    ),
                ),
                ui.tags.style(
                    """
                    .ks-filter-popup { display:none; }
                    .ks-filter-popup.active { display:block; }
                    .pathway-filter-button { gap: 0; }
                    .pathway-filter-button { gap: 0; }
                    .pathway-filter-button .fa { margin-right: 0; }
                    """
                ),
                ui.tags.script(
                    f"""
                    (function(){{
                        const btn = document.getElementById('{filter_btn_id}');
                        const popup = document.getElementById('{filter_popup_id}');
                        if (!btn || !popup) return;
                        btn.addEventListener('click', () => {{
                            popup.classList.toggle('active');
                        }});
                        document.addEventListener('click', (ev) => {{
                            if (!popup.contains(ev.target) && ev.target !== btn && !btn.contains(ev.target)) {{
                                popup.classList.remove('active');
                            }}
                        }}, true);
                    }})();
                    """
                ),
            )
        )
    else:
        children.append(search_input)
    if prefix == "web":
        # Hidden input to carry the selected pathway source (wikipathways/kegg)
        children.append(
            ui.div(
                {"style": "display:none;"},
                ui.input_text(
                    _prefixed_id(prefix, "pathway_source_choice"),
                    None,
                    value="wikipathways",
                    placeholder="",
                    width="100%",
                ),
            )
        )
    # Only non-web bookmarks show inline search results; web shows a full table instead.
    if prefix != "web":
        children.append(ui.output_ui(_prefixed_id(prefix, "pathway_search_results")))
    return ui.div(
        {"class": "pathway-search-wrapper", "id": _prefixed_id(prefix, "pathway_search_wrapper")},
        *children,
    )


def _simple_pathway_search_ui(prefix: str) -> ui.div:
    # Legacy KEGG-only search UI used for the Simple bookmark
    return ui.div(
        {"class": "pathway-search-wrapper", "id": _prefixed_id(prefix, "pathway_search_wrapper")},
        ui.input_text(
            _prefixed_id(prefix, "pathway_id"),
            "Pathway ID",
            value=DEFAULT_SETTINGS["pathway_id"],
            placeholder="Type regex to search KEGG pathways...",
        ),
        ui.output_ui(_prefixed_id(prefix, "pathway_search_results")),
    )


def _web_pathway_table_card(prefix: str) -> ui.card:
    return ui.card(
        ui.card_header("Pathway Catalogue"),
        _pathway_search_ui(prefix),
        ui.div(
            {"style": "margin-bottom:6px; color:#4b5563; font-size:0.9rem;"},
            "Click a row to select. Regex search above filters this table.",
        ),
        ui.output_ui(_prefixed_id(prefix, "pathway_table")),
    )


def _gradient_controls(prefix: str) -> ui.div:
    # Bookmark-level view only shows the shared gradient preview (values come from global Settings tab).
    return ui.div(
        {"class": "gradient-form"},
        ui.output_ui(_prefixed_id(prefix, "gradient_preview")),
    )


def _bookmark_gradient_preview(prefix: str) -> ui.div:
    return ui.div(
        {"class": "gradient-form"},
        ui.output_ui(_prefixed_id(prefix, "gradient_preview")),
    )


def _add_objects_panel(prefix: str) -> ui.card:
    settings_id = _prefixed_id(prefix, "spawn_settings_panel")
    mode_id = _prefixed_id(prefix, "spawn_layout_mode")
    grid_x_id = _prefixed_id(prefix, "spawn_grid_x")
    grid_y_id = _prefixed_id(prefix, "spawn_grid_y")
    grid_row_id = _prefixed_id(prefix, "spawn_grid_row")
    conc_arrows_id = _prefixed_id(prefix, "spawn_conc_arrows")
    conc_radius_mode_id = _prefixed_id(prefix, "spawn_conc_radius_mode")
    conc_radius_fixed_id = _prefixed_id(prefix, "spawn_conc_radius_fixed")
    conc_space_id = _prefixed_id(prefix, "spawn_conc_space")
    conc_arrow_stop_id = _prefixed_id(prefix, "spawn_conc_arrow_stop")
    conc_use_tooltip_id = _prefixed_id(prefix, "spawn_conc_use_tooltip")
    conc_tooltip_col_id = _prefixed_id(prefix, "spawn_conc_tooltip_col")

    return ui.card(
        ui.card_header("Add Objects"),
        ui.div(
            {"style": "display:flex; flex-direction:column; gap:8px;"},
            ui.div(
                {"style": "display:flex; gap:10px; align-items:center;"},
                ui.input_action_button(_prefixed_id(prefix, "spawn_protboxes"), "Add Protboxes", width="100%"),
            ),
            ui.input_text_area(
                _prefixed_id(prefix, "spawn_protboxes_ids"),
                "Uniprot IDs (comma or newline separated)",
                rows=8,
                placeholder="e.g. P04637, Q9Y243\nP31749\nQ02750",
                width="100%",
            ),
            ui.div(
                {"id": settings_id, "style": "border:1px solid #ddd; padding:10px; border-radius:8px; background:#fafafa;"},
                ui.input_radio_buttons(
                    mode_id,
                    "Layout Mode",
                    choices={"grid": "Grid layout", "concentric": "Concentric layout"},
                    selected="grid",
                    inline=True,
                ),
                ui.panel_conditional(
                    f"input.{mode_id} == 'grid'",
                    ui.div(
                        {"style": "display:flex; gap:10px; flex-wrap:wrap;"},
                        ui.input_numeric(grid_x_id, "Horizontal spacing (px)", value=75, min=10, max=400, width="180px"),
                        ui.input_numeric(grid_y_id, "Vertical spacing (px)", value=45, min=10, max=400, width="180px"),
                        ui.input_numeric(grid_row_id, "Boxes per row", value=5, min=1, max=20, width="140px"),
                    ),
                ),
                ui.panel_conditional(
                    f"input.{mode_id} == 'concentric'",
                    ui.div(
                        {"style": "display:flex; flex-direction:column; gap:8px;"},
                        ui.input_checkbox(conc_arrows_id, "Add arrows to center", value=True),
                        ui.input_radio_buttons(
                            conc_radius_mode_id,
                            "Ring radius",
                            choices={"auto": "Auto", "fixed": "Fixed"},
                            selected="auto",
                            inline=True,
                        ),
                        ui.panel_conditional(
                            f"input.{conc_radius_mode_id} == 'fixed'",
                            ui.input_numeric(conc_radius_fixed_id, "Fixed radius (px)", value=220, min=20, max=2000, width="200px"),
                        ),
                        ui.panel_conditional(
                            f"input.{conc_radius_mode_id} == 'auto'",
                            ui.input_numeric(conc_space_id, "Protbox space (auto)", value=70, min=10, max=400, width="200px"),
                        ),
                        ui.input_numeric(conc_arrow_stop_id, "Arrow stop radius (px)", value=50, min=0, max=500, width="220px"),
                        ui.input_checkbox(conc_use_tooltip_id, "Use tooltip column for encircling protboxes", value=False),
                        ui.input_select(
                            conc_tooltip_col_id,
                            "Tooltip column",
                            choices={c: c for c in (DEFAULT_SETTINGS.get("protein_tooltip_columns") or [])},
                            selected=None,
                        ),
                    ),
                ),
            ),
        ),
    )


def _settings_panel(cfg: Dict[str, Any]) -> ui.card:
    prefix = cfg["key"]
    controls: List[Any] = []
    top_controls: List[Any] = []
    is_ks = prefix == "ks"
    is_simple_web = prefix in {"simple", "web"}
    is_simple = prefix == "simple"
    ks_mode_toggle = None
    if is_ks:
        ks_sub_label = ui.span(
            "Substrates ",
            ui.span({"style": "color:#777; font-style:italic;"}, "(& Upstream Kinases)"),
        )
        ks_kin_label = ui.span(
            "Kinases ",
            ui.span({"style": "color:#777; font-style:italic;"}, "(& Downstream Substrates)"),
        )
        ks_mode_toggle = ui.div(
            {"class": "ks-mode-toggle", "style": "margin-bottom:8px;"},
            ui.input_radio_buttons(
                _prefixed_id(prefix, "ks_entity_mode"),
                None,
                choices={"substrate": ks_sub_label, "kinase": ks_kin_label},
                selected="substrate",
                inline=False,
            ),
        )
        controls.append(
            ui.div(
                {"class": "ks-filter-panel"},
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.tags.button(
                        ui.tags.i({"class": "fa fa-filter"}), "Filters",
                        {"id": _prefixed_id(prefix, "ks_filter_btn"), "type": "button", "class": "ks-filter-button"}
                    ),
                    ui.tags.button(
                        ui.tags.span("Refresh"),
                        {
                            "id": _prefixed_id(prefix, "ks_filter_refresh"),
                            "type": "button",
                            "class": "ks-filter-button",
                            "style": "background:#e5f4ff; color:#0f3d64;",
                            "title": "Refresh filters",
                            "onclick": f"Shiny.setInputValue('{_prefixed_id(prefix, 'ks_filter_refresh_evt')}', Date.now(), {{priority:'event'}});"
                        }
                    ),
                ),
                ui.div(
                    {"id": _prefixed_id(prefix, "ks_filter_popup"), "class": "ks-filter-popup"},
                    ui.div({"class": "ks-filter-row"},
                        ui.input_text(
                            _prefixed_id(prefix, "ks_filter_regex"),
                            "Regex search",
                            placeholder="Type regex...",
                        ),
                    ),
                    ui.tags.hr({"style": "margin:6px 0;"}),
                    ui.div(
                        {"style": "margin-top:6px;"},
                        ui.input_checkbox(
                            _prefixed_id(prefix, "ks_filter_reg_only"),
                            "Known regulatory sites only (substrates)",
                            value=False,
                        )
                    ),
                    ui.tags.hr({"style": "margin:6px 0;"}),
                    ui.div({"class": "ks-filter-row"},
                        ui.input_select(
                            _prefixed_id(prefix, "ks_filter_evidence"),
                            "Evidence",
                            choices={"both": "Both", "in_vivo": "in vivo", "in_vitro": "in vitro"},
                            selected="both",
                        ),
                    ),
                    ui.tags.hr({"style": "margin:6px 0;"}),
                    ui.div({"class": "ks-filter-row"},
                        ui.input_select(
                            _prefixed_id(prefix, "ks_filter_fc_op"),
                            "FC filter",
                            choices={ "": "Any", "gt": ">", "lt": "<", "ge": "â‰¥", "le": "â‰¤", "eq": "=", "ne": "â‰ "},
                            selected="",
                        ),
                        ui.input_text(
                            _prefixed_id(prefix, "ks_filter_fc_val"),
                            "Value",
                            placeholder="Enter number",
                        ),
                    ),
                    ui.tags.hr({"style": "margin:6px 0;"}),
                ),
                ui.tags.script(
                    f"""
                    (function(){{
                        const btn = document.getElementById('{_prefixed_id(prefix, "ks_filter_btn")}');
                        const popup = document.getElementById('{_prefixed_id(prefix, "ks_filter_popup")}');
                        const regexInput = document.getElementById('{_prefixed_id(prefix, "ks_filter_regex")}');
                        if (!btn || !popup) return;
                        const openPopup = () => {{
                            popup.classList.add('active');
                            if (regexInput) {{
                                regexInput.focus();
                                regexInput.select();
                            }}
                        }};
                        btn.addEventListener('click', () => {{
                            popup.classList.toggle('active');
                            if (popup.classList.contains('active')) openPopup();
                        }});
                        document.addEventListener('click', (ev) => {{
                            if (!popup.contains(ev.target) && ev.target !== btn && !btn.contains(ev.target)) {{
                                popup.classList.remove('active');
                            }}
                        }}, true);
                        document.addEventListener('keydown', (ev) => {{
                            if (ev.ctrlKey && ev.key.toLowerCase() === 'f') {{
                                const activeKs = document.querySelector('#bookmark_selector a.nav-link.active[data-value=\"ks\"]');
                                if (activeKs) {{
                                    ev.preventDefault();
                                    openPopup();
                                }}
                            }}
                        }});
                    }})();
                    """
                ),
            )
        )
        controls.append(ui.output_ui(_prefixed_id(prefix, "ks_table")))
    if not is_ks:
        if cfg.get("show_search", False):
            if is_simple_web:
                load_btn = ui.input_action_button(
                    _prefixed_id(prefix, "generate"),
                    "Load Pathway",
                    width="auto",
                    style="white-space: nowrap; min-width: 140px;"
                )
                if is_simple:
                    search_ui = _simple_pathway_search_ui(prefix)
                    top_controls.append(
                        ui.div(
                            {"style": "display:flex; align-items:center; gap:10px; flex-wrap:wrap;"},
                            search_ui,
                            load_btn,
                        )
                    )
                else:
                    selected_label = ui.output_text(_prefixed_id(prefix, "selected_pathway_label"))
                    top_controls.append(
                        ui.div(
                            {"style": "display:flex; align-items:flex-end; gap:10px; flex-wrap:wrap;"},
                            load_btn,
                            ui.div({"style": "font-weight:600; color:#1f2937;"}, selected_label),
                        )
                    )
            else:
                controls.append(_pathway_search_ui(prefix))
        elif cfg.get("key") != "figure":
            controls.append(ui.div({"class": "pathway-search-empty"}, "Pathway search is not available for this bookmark."))
    # Buttons and I/O ordering/titles
    clear_btn: Optional[Any] = None
    if not is_ks:
        if cfg.get("start_blank"):
            if cfg.get("key") == "figure":
                clear_btn = ui.input_action_button(_prefixed_id(prefix, "generate"), "Clear Canvas", width="160px")
            else:
                controls.append(ui.input_action_button(_prefixed_id(prefix, "generate"), "Create Blank Canvas", width="100%"))
        else:
            if not is_simple_web:
                controls.append(ui.input_action_button(_prefixed_id(prefix, "generate"), "Generate Pathway", width="100%"))
    export_btn: Optional[Any] = None
    import_btn: Optional[Any] = None
    export_script: Optional[Any] = None
    if cfg.get("show_downloads", True):
        if cfg.get("key") not in {"figure", "web"}:
            controls.append(ui.download_button(_prefixed_id(prefix, "download_json"), "Download JSON", width="100%"))
        export_label = "Export Pathway" if cfg.get("key") == "figure" else "Export Custom Pathway"
        export_btn_id = _prefixed_id(prefix, "export_custom_pathway")
        export_spinner_id = _prefixed_id(prefix, "export_custom_pathway_spinner")
        export_action_btn = ui.input_action_button(export_btn_id, export_label, width="auto")
        export_btn = ui.div(
            {"style": "display:flex; align-items:center; gap:8px;"},
            export_action_btn,
            ui.tags.span({"id": export_spinner_id, "class": "mk-export-spinner", "style": "display:none;"}, "Preparing..."),
        )
        export_script = ui.tags.script(
            f"""
            (function(){{
                const btn = document.getElementById('{export_btn_id}');
                const spinner = document.getElementById('{export_spinner_id}');
                if (!btn || !spinner) return;
                btn.addEventListener('click', () => {{
                    spinner.style.display = 'inline-flex';
                    btn.disabled = true;
                }});
            }})();
            """
        )
    if cfg.get("show_custom_io", True) and cfg.get("key") != "web":
        import_btn = ui.input_file(_prefixed_id(prefix, "upload_custom_pathway"), "Import Custom Pathway", accept=[".json"], multiple=False)
    if import_btn or export_btn:
        io_children = [btn for btn in (import_btn, export_btn) if btn is not None]
        controls.append(
            ui.div(
                {"style": "display:flex; gap:10px; flex-wrap:wrap; justify-content:flex-start; align-items:center;"},
                *io_children,
            )
        )
        if export_script is not None:
            controls.append(export_script)
    if cfg.get("show_toggles", True):
        toggle_controls = []
        if cfg.get("key") not in {"web", "figure"}:
            toggle_controls.extend([
                ui.input_checkbox(_prefixed_id(prefix, "show_arrows"), "Show arrows", value=_bool_default(cfg, "show_arrows")),
                ui.input_checkbox(_prefixed_id(prefix, "show_text_boxes"), "Show text boxes", value=_bool_default(cfg, "show_text_boxes")),
                ui.input_checkbox(_prefixed_id(prefix, "show_background_image"), "Background Image", value=_bool_default(cfg, "show_background_image")),
            ])
        controls.extend(toggle_controls)
    header_rows: List[Any] = []
    if is_ks:
        if ks_mode_toggle:
            header_rows.append(ks_mode_toggle)
        header_rows.append(
            ui.div(
                {
                    "style": (
                        "display:flex; gap:12px; align-items:flex-start; flex-wrap:wrap; "
                        "margin-bottom:8px;"
                    )
                },
                ui.div({"style": "flex:1; min-width:260px;"}, _bookmark_gradient_preview(prefix)),
                ui.div({"style": "width:220px; min-width:200px;"}, ui.output_ui(_prefixed_id(prefix, "fc_selector"))),
            )
        )
    else:
        if cfg.get("key") == "figure" and cfg.get("start_blank"):
            header_rows.append(
                ui.div(
                    {
                        "style": (
                            "display:flex; gap:12px; align-items:center; flex-wrap:wrap; "
                            "margin-bottom:8px;"
                        )
                    },
                    ui.div({"style": "flex:1; min-width:260px;"}, _bookmark_gradient_preview(prefix)),
                    clear_btn if clear_btn else ui.input_action_button(_prefixed_id(prefix, "generate"), "Clear Canvas", width="160px"),
                )
            )
            header_rows.append(
                ui.div(
                    {"style": "margin-bottom:8px;"},
                    ui.output_ui(_prefixed_id(prefix, "fc_selector")),
                )
            )
        else:
            header_rows.append(_bookmark_gradient_preview(prefix))
            header_rows.append(ui.output_ui(_prefixed_id(prefix, "fc_selector")))

    return ui.card(
        ui.card_header(f"{cfg['label']} Settings"),
        *top_controls,
        *header_rows,
        *controls,
    )


def _preview_panel(cfg: Dict[str, Any]) -> ui.card:
    prefix = cfg["key"]
    download_name = f"{prefix}_pathway"
    download_controls = ui.div(
        {"class": "svg-download-row"},
        ui.span({"class": "svg-download-label"}, "Download:"),
        ui.tags.button(
            {
                "class": "btn btn-primary btn-sm",
                "data-mk-download": "svg",
                "data-mk-name": download_name,
                "type": "button",
            },
            "SVG",
        ),
        ui.tags.button(
            {
                "class": "btn btn-outline-primary btn-sm",
                "data-mk-download": "pdf",
                "data-mk-name": download_name,
                "type": "button",
            },
            "PDF",
        ),
        ui.tags.button(
            {
                "class": "btn btn-outline-primary btn-sm",
                "data-mk-download": "png",
                "data-mk-name": download_name,
                "type": "button",
            },
            "PNG",
        ),
        ui.tags.button(
            {
                "class": "btn btn-outline-primary btn-sm",
                "data-mk-download": "jpeg",
                "data-mk-name": download_name,
                "type": "button",
            },
            "JPEG",
        ),
    )
    return ui.card(
        ui.card_header(f"{cfg['label']} Viewer"),
        ui.div(
            {
                "class": "pathway-viewer-card",
                "data-prefix": prefix,
                "data-download-name": download_name,
            },
            ui.output_ui(_prefixed_id(prefix, "pathway_preview")),
            download_controls,
        ),
        ui.hr(),
        ui.output_text(_prefixed_id(prefix, "status_message")),
        ui.output_text(_prefixed_id(prefix, "json_summary")),
    )


def _bookmark_tab(cfg: Dict[str, Any]) -> ui.nav_panel:
    if cfg.get("key") == "figure":
        settings_col = 4
        preview_col = 12 - settings_col
        return ui.nav_panel(
            cfg["label"],
            ui.row(
                ui.column(settings_col, ui.TagList(_settings_panel(cfg), _add_objects_panel(cfg["key"]))),
                ui.column(preview_col, _preview_panel(cfg)),
            ),
            value=cfg["key"],
        )
    settings_col = 6 if cfg.get("key") == "ks" else 4
    preview_col = 12 - settings_col
    settings_children: Any = _settings_panel(cfg)
    if cfg.get("key") == "web":
        settings_children = ui.TagList(settings_children, _web_pathway_table_card(cfg["key"]))
    return ui.nav_panel(
        cfg["label"],
        ui.row(
            ui.column(settings_col, settings_children),
            ui.column(preview_col, _preview_panel(cfg)),
        ),
        value=cfg["key"],
    )


input_tab = ui.nav_panel(
    "Input",
    ui.div(
        {"class": "input-bg"},
        ui.div(
            {"class": "input-wrapper"},
            ui.output_ui("input_data_inputs_panel"),
        ),
    ),
    value="input",
)


app_ui = ui.page_fluid(
    CUSTOM_STYLES,
    ui.tags.style(
        """
        .input-bg { background: #0b4f9c; min-height: 100vh; padding: 32px 24px; }
        .input-wrapper { min-height: 100vh; display: flex; align-items: flex-start; justify-content: flex-start; padding-top: 16px; }
        .data-input-card { max-width: 520px; width: 100%; box-shadow: 0 18px 40px rgba(0,0,0,0.25); border-radius: 16px; }
        """
    ),
    NAV_LOCK_SCRIPT,
    EXPORT_DOWNLOAD_SCRIPT,
    SVG_DOWNLOAD_SCRIPT,
    ui.h2("MapKinase Demo"),
    ui.navset_tab(
        input_tab,
        *[_bookmark_tab(cfg) for cfg in BOOKMARK_CONFIGS],
        ui.nav_panel(
            "Settings",
            ui.div(
                {"style": "display:flex; gap:16px; align-items:flex-start; flex-wrap:wrap;"},
                ui.div(
                    {"style": "flex:1; min-width:320px;"},
                    ui.card(
                    ui.card_header("Global Settings"),
                    ui.div(
                        {"style": "margin-bottom: 10px; color:#555;"},
                        "These settings apply to all bookmarks. Changes take effect when you load or generate a new pathway."
                    ),
                    ui.input_numeric(
                        "settings_prot_label_size",
                        "Protein Label Size",
                        value=DEFAULT_SETTINGS["prot_label_size"],
                        min=1,
                        max=72,
                        step=1,
                    ),
                    ui.input_numeric(
                        "settings_ptm_max_display",
                        "PTM Max Display",
                        value=DEFAULT_SETTINGS["ptm_max_display"],
                        min=0,
                    ),
                    ui.input_numeric(
                        "settings_prot_outline_width",
                        "Protbox Outline Width",
                        value=DEFAULT_SETTINGS.get("prot_outline_width", 1),
                        min=0,
                        step=0.25,
                    ),
                    ui.input_numeric(
                        "settings_ptm_outline_width",
                        "PTM Outline Width",
                        value=DEFAULT_SETTINGS.get("ptm_outline_width", 1),
                        min=0,
                        step=0.25,
                    ),
                    ui.input_checkbox(
                        "settings_use_original_protbox_size",
                        "Use original protbox size (for PTM placement)",
                        value=DEFAULT_SETTINGS.get("use_original_protbox_size", False),
                    ),
                    ui.input_checkbox(
                        "settings_include_psp_tooltips",
                        "Include PSP annotations as PTM tooltips",
                        value=True,
                    ),
                    ui.input_checkbox(
                        "settings_show_multi_protein_indicator",
                        "Show multi-protein indicator (all bookmarks)",
                        value=_bool_default({"mode": "analysis"}, "show_multi_protein_indicator"),
                    ),
                    ui.input_checkbox(
                        "settings_show_groups",
                        "Show groups (all bookmarks)",
                        value=_bool_default({"mode": "analysis"}, "show_groups"),
                    ),
                    ui.input_checkbox(
                        "settings_debug_mode",
                        "Debug Mode",
                        value=DEFAULT_SETTINGS.get("debug_mode", False),
                    ),
                    ),
                ),
                ui.div(
                    {"style": "flex:1; min-width:320px;"},
                    ui.card(
                    ui.card_header("Gradient Settings"),
                    ui.div({"style": "font-weight:700; margin-bottom:6px;"}, "Gradient (applies to all bookmarks)"),
                    ui.input_select(
                        "settings_gradient_preset",
                        "Color preset",
                        choices={
                            "main_default": "Main Default (Yellow / Blue)",
                            "tan_turquoise": "Tan / Turquoise",
                            "orange_purple": "Orange / Purple",
                            "green_purple": "Green / Purple",
                            "red_blue": "Red / Blue",
                            "pink_blue": "Pink / Blue",
                            "yellow_pink": "Yellow / Pink",
                            "brown_blue": "Brown / Blue",
                            "custom": "Custom",
                        },
                        selected="main_default",
                    ),
                    ui.layout_columns(
                        ui.column(
                            3,
                            ui.input_numeric(
                                "settings_max_negative",
                                "Max negative",
                                value=DEFAULT_SETTINGS["max_negative"],
                                step=0.1,
                                width="140px",
                            ),
                        ),
                        ui.column(
                            3,
                            _color_picker_input(
                                "settings_negative_color",
                                "Negative color",
                                _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]),
                            ),
                        ),
                        ui.column(
                            3,
                            ui.input_numeric(
                                "settings_max_positive",
                                "Max positive",
                                value=DEFAULT_SETTINGS["max_positive"],
                                step=0.1,
                                width="140px",
                            ),
                        ),
                        ui.column(
                            3,
                            _color_picker_input(
                                "settings_positive_color",
                                "Positive color",
                                _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]),
                            ),
                        ),
                    ),
                    ui.output_ui("settings_gradient_preview"),
                    ),
                ),
            ),
            value="settings",
        ),
        id="bookmark_selector",
        selected="input",
    ),
)


def server(input, output, session):  # type: ignore[override]
    bookmark_state: Dict[str, Dict[str, reactive.Value]] = {}
    for cfg in BOOKMARK_CONFIGS:
        bookmark_state[cfg["key"]] = {
            "json": reactive.Value(None),
            "status": reactive.Value(STATUS_READY),
            "custom_layout": reactive.Value(None),
            "fc_index": reactive.Value(1),
            "export_snapshot": reactive.Value(None),
            "export_pending": reactive.Value(False),
        }

    protein_validation = reactive.Value(
        {"status": "Upload a protein file to begin.", "errors": [], "valid": False, "comparisons": []}  # type: ignore[var-annotated]
    )
    ptm_validation = reactive.Value(
        {"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False}  # type: ignore[var-annotated]
    )
    nav_lock_status = reactive.Value("User mode: upload valid Protein file to unlock other tabs. PTM optional.")
    extra_datasets = reactive.Value([])  # type: ignore[var-annotated]
    protein_dataset = reactive.Value(None)
    ptm_dataset = reactive.Value(None)
    protein_dataset_path = reactive.Value(None)
    ptm_dataset_path = reactive.Value(None)
    psp_cache: Dict[str, Any] = {}
    kegg_cache: Dict[str, Dict[str, str]] = {}
    protein_kegg_warning = reactive.Value("")
    ks_index = reactive.Value(_empty_ks_index())

    def _active_bookmark() -> str:
        active = _get_input_value(input, "bookmark_selector")
        if active in {cfg["key"] for cfg in BOOKMARK_CONFIGS}:
            return active
        return BOOKMARK_CONFIGS[0]["key"]

    def _current_mode() -> str:
        mode = _get_input_value(input, "input_mode")
        return str(mode or "user").lower()

    def _finalize_custom_export(prefix: str, snapshot: Dict[str, Any]) -> None:
        if prefix not in bookmark_state:
            return
        state = bookmark_state[prefix]
        try:
            export_payload = _build_custom_layout_export(snapshot)
            export_json = json.dumps(export_payload, indent=2, ensure_ascii=True)
            _send_custom_message(
                session,
                "download_payload",
                {
                    "filename": f"{prefix}_custom_pathway.json",
                    "content": export_json,
                    "button_id": _prefixed_id(prefix, "export_custom_pathway"),
                    "spinner_id": _prefixed_id(prefix, "export_custom_pathway_spinner"),
                },
            )
        except Exception as exc:
            state["status"].set(f"Failed to export custom pathway: {exc}")
            print(f"_finalize_custom_export error: {exc}")
            _send_custom_message(
                session,
                "export_failed",
                {
                    "button_id": _prefixed_id(prefix, "export_custom_pathway"),
                    "spinner_id": _prefixed_id(prefix, "export_custom_pathway_spinner"),
                },
            )
        finally:
            state["export_pending"].set(False)

    def _detect_delimiter(path: str) -> str:
        ext = os.path.splitext(path)[1].lower()
        return "," if ext == ".csv" else "\t"

    def _load_dataset(path: str) -> Dict[str, Any]:
        delimiter = _detect_delimiter(path)
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                raw_lines = fh.read().splitlines()
            reader = csv.reader(raw_lines, delimiter=delimiter)
            parsed = list(reader)
            if not parsed:
                return {"error": "No rows parsed from file."}
            headers = []
            for idx, h in enumerate(parsed[0]):
                cleaned = h or ""
                if idx == 0:
                    cleaned = cleaned.lstrip("\ufeff")
                headers.append(cleaned.strip())
            rows = parsed[1:]
            return {"headers": headers, "rows": rows}
        except Exception as exc:
            return {"error": f"Debug load failed: {exc}"}

    def _write_debug_dump(filename: str, payload: Dict[str, Any]) -> None:
        if not debug_var:
            return
        try:
            debug_path = os.path.join(BASE_DIR, filename)
            with open(debug_path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
        except Exception as exc:
            print(f"Failed to write debug dump {filename}: {exc}")

    @reactive.Effect
    @reactive.event(input.export_snapshot)
    def _handle_export_snapshot():
        payload = _get_input_value(input, "export_snapshot") or {}
        prefix = str(payload.get("prefix") or "")
        snapshot = payload.get("payload")
        if not prefix or prefix not in bookmark_state:
            return
        state = bookmark_state[prefix]
        if not state.get("export_pending") or not state["export_pending"].get():
            return
        if not isinstance(snapshot, dict):
            state["status"].set("Export snapshot was invalid.")
            _send_custom_message(
                session,
                "export_failed",
                {
                    "button_id": _prefixed_id(prefix, "export_custom_pathway"),
                    "spinner_id": _prefixed_id(prefix, "export_custom_pathway_spinner"),
                },
            )
            state["export_pending"].set(False)
            return
        state["export_snapshot"].set(snapshot)
        _finalize_custom_export(prefix, snapshot)

    def _split_multi_ids(raw: Any) -> List[str]:
        if raw in (None, ""):
            return []
        if isinstance(raw, (list, tuple)):
            values = []
            for item in raw:
                values.extend(_split_multi_ids(item))
            return values
        text = str(raw)
        parts = re.split(r"[;,]", text)
        return [p.strip() for p in parts if p and p.strip()]

    def _clean_site_label(site_raw: str) -> str:
        if not site_raw:
            return ""
        site_str = str(site_raw).strip()
        # If it looks like an Excel-style month conversion (e.g., 2-Sep), fall back to digit+optional AA letter.
        if re.match(r"^\d{1,2}[-/][A-Za-z]{3}$", site_str):
            match = re.search(r"([A-Za-z])?(\d+)", site_str)
            if match:
                letter, digits = match.group(1), match.group(2)
                return f"{letter or ''}{digits}"
            return re.sub(r"[-/]", "", site_str)
        match = re.search(r"([A-Za-z])?(\d+)", site_str)
        if match:
            letter, digits = match.group(1), match.group(2)
            return f"{letter or ''}{digits}"
        return site_str

    def _clean_protein_label(raw: Any, fallback: str) -> str:
        text = str(raw or "").strip()
        if not text:
            return fallback
        if re.match(r"^\d{1,2}[-/][A-Za-z]{3}$", text):
            return fallback
        return text

    def _build_ks_index(ptm_data: Optional[Dict[str, Any]], prot_data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        if not ptm_data or not ptm_data.get("headers"):
            return _empty_ks_index()
        headers = list(ptm_data.get("headers") or [])
        if len(headers) < 2:
            return _empty_ks_index()
        idx_map = {h: i for i, h in enumerate(headers)}
        vivo_col = "PSP: uniprot_in_vivo_kinases"
        vitro_col = "PSP: uniprot_in_vitro_kinases"
        if vivo_col not in idx_map and vitro_col not in idx_map:
            return _empty_ks_index()
        prot_gene_map: Dict[str, str] = {}
        if prot_data:
            prot_headers = prot_data.get("headers") or []
            if len(prot_headers) >= 2:
                for row in prot_data.get("rows") or []:
                    uid = (row[0] if len(row) > 0 else "").strip()
                    gene = (row[1] if len(row) > 1 else "").strip()
                    if uid:
                        prot_gene_map[uid] = gene or uid
        kinases: Dict[str, Dict[str, Any]] = {}
        substrates: Dict[str, Dict[str, Any]] = {}
        ptms_by_uniprot: Dict[str, List[Dict[str, Any]]] = {}
        ptm_rows = ptm_data.get("rows") or []
        for row in ptm_rows:
            row_vals = list(row)
            row_map = {h: (row_vals[idx] if idx < len(row_vals) else "") for h, idx in idx_map.items()}
            uniprot = str(row_map.get(headers[0], "") or "").strip()
            site = str(row_map.get(headers[1], "") or "").strip()
            site_label = _clean_site_label(site)
            if not uniprot:
                continue
            sub_key = f"{uniprot}:{site}" if site else uniprot
            ptms_by_uniprot.setdefault(uniprot, []).append(row_map)
            vivo_ids = _split_multi_ids(row_map.get(vivo_col, ""))
            vitro_ids = _split_multi_ids(row_map.get(vitro_col, ""))
            if not (vivo_ids or vitro_ids):
                continue
            sub_entry = substrates.setdefault(
                sub_key,
                {
                    "uniprot": uniprot,
                    "site": site,
                    "site_label": site_label,
                    "types": set(),
                    "gene": prot_gene_map.get(uniprot, uniprot),
                    "kinases": {},
                    "row": row_map,
                },
            )
            if vivo_ids:
                sub_entry["types"].add("in_vivo")
            if vitro_ids:
                sub_entry["types"].add("in_vitro")
            for kin_id in vivo_ids:
                kin_entry = kinases.setdefault(
                    kin_id,
                    {"types": set(), "gene": prot_gene_map.get(kin_id, kin_id), "substrates": set()},
                )
                kin_entry["types"].add("in_vivo")
                kin_entry["substrates"].add(sub_key)
                rel = sub_entry["kinases"].setdefault(kin_id, {"types": set()})
                rel["types"].add("in_vivo")
            for kin_id in vitro_ids:
                kin_entry = kinases.setdefault(
                    kin_id,
                    {"types": set(), "gene": prot_gene_map.get(kin_id, kin_id), "substrates": set()},
                )
                kin_entry["types"].add("in_vitro")
                kin_entry["substrates"].add(sub_key)
                rel = sub_entry["kinases"].setdefault(kin_id, {"types": set()})
                rel["types"].add("in_vitro")
        return {
            "kinases": kinases,
            "substrates": substrates,
            "ptms_by_uniprot": ptms_by_uniprot,
            "ptm_headers": headers,
            "prot_gene_map": prot_gene_map,
        }

    def _update_ks_index(reset: bool = False):
        if reset:
            ks_index.set(_empty_ks_index())
            return
        with reactive.isolate():
            ptm_data = ptm_dataset.get()
            prot_data = protein_dataset.get()
        ks_index.set(_build_ks_index(ptm_data, prot_data))

    def _get_psp_map():
        cached = psp_cache.get("data")
        if not cached:
            psp_cache["data"] = load_regulatory_sites(BASE_DIR)
        return psp_cache.get("data") or {}

    def _get_psp_ks_map():
        cached = psp_cache.get("ks")
        if not cached:
            psp_cache["ks"] = load_kinase_substrate_map(BASE_DIR)
        return psp_cache.get("ks") or {}

    def _get_kegg_map(species_code: str):
        key = (species_code or "").strip().lower()
        if key not in kegg_cache:
            kegg_cache[key] = load_kegg_map(BASE_DIR, key)
        return kegg_cache.get(key) or {}

    def _send_custom_message_safe(msg_type: str, payload: Dict[str, Any]) -> None:
        try:
            result = session.send_custom_message(msg_type, payload)
            if asyncio.iscoroutine(result):
                asyncio.create_task(result)
        except Exception as exc:
            print(f"Warning: send_custom_message failed: {exc}")

    def _load_demo_datasets():
        """Load bundled demo protein/PTM files so Demo mode behaves like preloaded uploads."""
        missing = []
        if not os.path.exists(SAMPLE_PROTEIN_FILE):
            missing.append(f"Protein file not found: {SAMPLE_PROTEIN_FILE}")
        if not os.path.exists(SAMPLE_PTM_FILE):
            missing.append(f"PTM file not found: {SAMPLE_PTM_FILE}")
        if missing:
            protein_validation.set({"status": "Demo mode files missing.", "errors": missing, "valid": False, "comparisons": []})
            ptm_validation.set({"status": "Demo mode files missing.", "errors": missing, "valid": False})
            protein_dataset.set(None)
            ptm_dataset.set(None)
            nav_lock_status.set("Demo mode: sample files missing. Navigation locked.")
            return

        try:
            prot_result = validate_protein_file(SAMPLE_PROTEIN_FILE)
            if not prot_result.valid:
                protein_validation.set({"status": "Demo protein file failed validation.", "errors": prot_result.errors, "valid": False, "comparisons": []})
                ptm_validation.set({"status": "PTM upload disabled until demo protein is valid.", "errors": [], "valid": False})
                protein_dataset.set(None)
                ptm_dataset.set(None)
                nav_lock_status.set("Demo mode: sample protein validation failed.")
                return

            species_choice, species_info = _resolve_species(_get_input_value(input, "input_species"))
            if not species_info:
                species_choice = DEFAULT_SPECIES
                species_info = SPECIES_CHOICES.get(DEFAULT_SPECIES, {})
            species_code = species_info.get("code", "") or SPECIES_CHOICES.get(DEFAULT_SPECIES, {}).get("code", "hsa")
            print(f"Demo mode: loading sample datasets for species '{species_choice}' code '{species_code}'")

            demo_prot_payload = _load_dataset(SAMPLE_PROTEIN_FILE)
            try:
                demo_prot_payload = annotate_protein_with_kegg(demo_prot_payload, species_code, _get_kegg_map(species_code))
            except Exception as exc:
                print(f"Warning: demo protein KEGG annotation failed: {exc}")
            protein_dataset.set(demo_prot_payload)
            protein_dataset_path.set(SAMPLE_PROTEIN_FILE)
            protein_validation.set({
                "status": (
                    f"Demo protein loaded. Rows: {prot_result.summary.get('rows', 0)}, "
                    f"Comparison columns: {prot_result.summary.get('comparisons', 0)}, "
                    f"Tooltip columns: {prot_result.summary.get('tooltips', 0)}."
                ),
                "errors": [],
                "valid": True,
                "comparisons": prot_result.comparisons,
            })
            _write_debug_dump("user_protein_dataset_debug.txt", demo_prot_payload)

            ptm_result = validate_ptm_file(SAMPLE_PTM_FILE, prot_result.comparisons)
            if not ptm_result.valid:
                ptm_validation.set({"status": "Demo PTM file failed validation.", "errors": ptm_result.errors, "valid": False})
                ptm_dataset.set(None)
                nav_lock_status.set("Demo mode: sample PTM validation failed.")
                return

            demo_ptm_payload = _load_dataset(SAMPLE_PTM_FILE)
            try:
                demo_ptm_payload = annotate_ptm_dataset(demo_ptm_payload, species_choice, _get_psp_map())
                demo_ptm_payload = annotate_ptm_dataset_with_kinases(demo_ptm_payload, species_choice, _get_psp_ks_map())
            except Exception as exc:
                print(f"Warning: demo PTM annotation failed: {exc}")
            ptm_dataset.set(demo_ptm_payload)
            ptm_dataset_path.set(SAMPLE_PTM_FILE)
            ptm_validation.set({
                "status": (
                    f"Demo PTM loaded. Rows: {ptm_result.summary.get('rows', 0)}, "
                    f"Comparison columns: {ptm_result.summary.get('comparisons', 0)}, "
                    f"Tooltip columns: {ptm_result.summary.get('tooltips', 0)}."
                ),
                "errors": [],
                "valid": True,
            })
            _write_debug_dump("user_ptm_dataset_debug.txt", demo_ptm_payload)
            _update_ks_index()
            nav_lock_status.set("Demo mode: sample datasets loaded. Navigation unlocked.")
            print("Demo mode: sample protein/PTM loaded successfully.")
        except Exception as exc:
            print(f"Warning: demo mode encountered an unexpected error: {exc}")
        except Exception as exc:
            print(f"Warning: demo mode encountered an unexpected error: {exc}")
            protein_validation.set({"status": "Demo mode failed to load sample protein.", "errors": [str(exc)], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "Demo mode failed to load sample PTM.", "errors": [str(exc)], "valid": False})
            protein_dataset.set(None)
            ptm_dataset.set(None)
            nav_lock_status.set("Demo mode: sample datasets failed to load.")

    def collect_data_override() -> Dict[str, Any]:  # type: ignore[override]
        mode = _get_input_value(input, "input_mode") or "user"
        if str(mode).lower() not in {"user", "demo"}:
            return {}
        with reactive.isolate():
            protein_valid = protein_validation.get().get("valid")
            ptm_valid = ptm_validation.get().get("valid")
            protein_data = protein_dataset.get()
            ptm_data = ptm_dataset.get()
        if not protein_valid:
            return {}
        if not protein_data:
            return {}
        species_choice, species_info = _resolve_species(_get_input_value(input, "input_species"))
        species_code = species_info.get("code", "")
        shape_choice = (_get_input_value(input, "input_ptm_shape") or "circle").strip().title()

        # Protein columns: first = uniprot, second = gene, KEGG column added by annotator
        prot_headers = list(protein_data.get("headers") or [])
        if len(prot_headers) < 2:
            return {}
        prot_uniprot = prot_headers[0]
        prot_gene = prot_headers[1]
        prot_kegg = "KEGG_Gene_ID"
        if prot_kegg not in prot_headers:
            prot_headers.append(prot_kegg)
            padded_rows: List[List[str]] = []
            for row in protein_data.get("rows") or []:
                new_row = list(row)
                new_row.append("")
                padded_rows.append(new_row)
            protein_data["rows"] = padded_rows
            protein_data["headers"] = prot_headers
        prot_main = [h for h in prot_headers if h.startswith("C:")]
        prot_tooltips = [h for h in prot_headers if h.startswith("T:")]
        prot_outline = _resolve_outline_columns(prot_main, prot_headers)

        # PTM columns
        ptm_headers = list(ptm_data.get("headers") or []) if ptm_data else []
        ptm_entries: List[Dict[str, Any]] = []
        ptm_tooltips: List[str] = []
        modulation_col = "PSP: regulatory_site"
        ptm_symbol_list = copy.deepcopy(DEFAULT_DATA["ptm"][0]["ptm_symbol_list"])
        # Update symbol headers: most use PSP: ON_FUNCTION, but label_4 should use PSP: regulatory_site
        for item in ptm_symbol_list:
            for key, val in item.items():
                if not isinstance(val, dict) or not val.get("header_to_search"):
                    continue
                if key == "symbol_label_4_dict":
                    val["header_to_search"] = "PSP: regulatory_site"
                else:
                    val["header_to_search"] = "PSP: ON_FUNCTION"
        if ptm_headers and len(ptm_headers) >= 2 and ptm_data:
            ptm_uniprot = ptm_headers[0]
            ptm_site = ptm_headers[1]
            added_ptm_cols: List[str] = []
            for needed in ["PSP: regulatory_site", "PSP: ON_FUNCTION"]:
                if needed not in ptm_headers:
                    ptm_headers.append(needed)
                    added_ptm_cols.append(needed)
            if added_ptm_cols:
                padded_rows = []
                for row in ptm_data.get("rows") or []:
                    new_row = list(row)
                    new_row.extend([""] * len(added_ptm_cols))
                    padded_rows.append(new_row)
                ptm_data["rows"] = padded_rows
                ptm_data["headers"] = ptm_headers
            ptm_main = [(h, h) for h in ptm_headers if h.startswith("C:")]
            ptm_main_cols = [col for col, _ in ptm_main]
            ptm_outline = _resolve_outline_columns(ptm_main_cols, ptm_headers)
            ptm_tooltips = [h for h in ptm_headers if h.startswith("T:")]
            ptm_entries = [
                {
                    "type": "Phosphorylation",
                    "data_headers": ptm_headers,
                    "data_rows": ptm_data.get("rows") or [],
                    "uniprot_column": ptm_uniprot,
                    "site_column": ptm_site,
                    "shape": shape_choice,
                    "main_columns": ptm_main,
                    "outline_columns": ptm_outline,
                    "modulation_column": modulation_col,
                    "tooltip_columns": ptm_tooltips,
                    "ptm_symbol_list": ptm_symbol_list,
                }
            ]
        else:
            ptm_entries = []

        # Optionally add PSP annotation columns to tooltip list (exclude modulation_col)
        include_psp = bool(_get_input_value(input, "settings_include_psp_tooltips"))
        if include_psp:
            psp_cols = [
                "PSP: DOMAIN",
                "PSP: ON_FUNCTION",
                "PSP: ON_PROCESS",
                "PSP: ON_PROT_INTERACT",
                "PSP: ON_OTHER_INTERACT",
                "PSP: NOTES",
                "PSP: in_vivo_kinases",
                "PSP: in_vitro_kinases",
            ]
            for col in psp_cols:
                if col in ptm_headers and col not in ptm_tooltips and col != modulation_col:
                    ptm_tooltips.append(col)

        data_override = {
            "protein": {
                "data_headers": prot_headers,
                "data_rows": protein_data.get("rows") or [],
                "file_path": protein_dataset_path.get() or "",
                "uniprot_column": prot_uniprot,
                "kegg_column": prot_kegg,
                "gene_column": prot_gene,
                "main_columns": prot_main,
                "outline_columns": prot_outline,
                "tooltip_columns": prot_tooltips,
            },
            "ptm": ptm_entries,
        }
        # Attach file_path to each PTM entry so downstream processors/catalog use the correct dataset
        for entry in data_override["ptm"]:
            entry["file_path"] = ptm_dataset_path.get() or ""
        return data_override

    @reactive.Effect
    @reactive.event(input.input_protein_upload)
    def _process_protein_upload():
        if str((_get_input_value(input, "input_mode") or "user")).lower() == "demo":
            protein_validation.set({"status": "Demo mode uses bundled sample files; switch to User mode to upload.", "errors": [], "valid": False, "comparisons": []})
            protein_kegg_warning.set("")
            return
        upload = input.input_protein_upload()
        if not upload:
            protein_validation.set({"status": "Upload a protein file to begin.", "errors": [], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "PTM upload disabled until a valid protein file is uploaded.", "errors": [], "valid": False})
            protein_dataset.set(None)
            _write_debug_dump("user_protein_dataset_debug.txt", {"info": "No dataset loaded"})
            protein_kegg_warning.set("")
            _update_ks_index(reset=True)
            return
        file_info = upload[0]
        protein_validation.set({"status": "Processing protein upload...", "errors": [], "valid": False, "comparisons": []})
        datapath = getattr(file_info, "datapath", None)
        if datapath is None and isinstance(file_info, dict):
            datapath = file_info.get("datapath")
        if not datapath:
            protein_validation.set({"status": "Could not locate the uploaded protein file.", "errors": [], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
            protein_dataset.set(None)
            _write_debug_dump("user_protein_dataset_debug.txt", {"info": "No dataset loaded"})
            protein_kegg_warning.set("")
            _update_ks_index(reset=True)
            return
        result = validate_protein_file(datapath)
        if not result.valid:
            protein_validation.set({"status": "Protein file failed validation. See errors below.", "errors": result.errors, "valid": False, "comparisons": []})
            ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
            protein_dataset.set(None)
            protein_dataset_path.set(None)
            _write_debug_dump("user_protein_dataset_debug.txt", {"errors": result.errors})
            protein_kegg_warning.set("")
            _update_ks_index(reset=True)
            return

        status = (
            f"Protein file valid. Rows: {result.summary.get('rows', 0)}, "
            f"Comparison columns: {result.summary.get('comparisons', 0)}, "
            f"Tooltip columns: {result.summary.get('tooltips', 0)}."
        )
        protein_validation.set({"status": status, "errors": [], "valid": True, "comparisons": result.comparisons})
        dataset_payload = _load_dataset(datapath)
        species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
        species_code = species_info.get("code", "")
        try:
            kegg_map = _get_kegg_map(species_code)
            dataset_payload = annotate_protein_with_kegg(dataset_payload, species_code, kegg_map)
            headers = dataset_payload.get("headers") or []
            kegg_idx = headers.index("KEGG_Gene_ID") if "KEGG_Gene_ID" in headers else -1
            kegg_matches = 0
            if kegg_idx >= 0:
                for row in dataset_payload.get("rows") or []:
                    if len(row) > kegg_idx and (row[kegg_idx] or "").strip():
                        kegg_matches += 1
            if kegg_matches <= 5:
                protein_kegg_warning.set("Few KEGG matches detected. Check species selection or mapping file.")
            else:
                protein_kegg_warning.set("")
        except Exception as exc:
            print(f"Warning: KEGG annotation failed: {exc}")
            protein_kegg_warning.set("")
        protein_dataset.set(dataset_payload)
        protein_dataset_path.set(datapath)
        _write_debug_dump("user_protein_dataset_debug.txt", dataset_payload)
        ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})

    @output
    @render.text
    def input_protein_status():
        data = protein_validation.get()
        return data.get("status", "")

    @output
    @render.ui
    def input_protein_errors():
        data = protein_validation.get()
        errs = data.get("errors") or []
        if not errs:
            return ui.div()
        items = [ui.tags.li(err) for err in errs[:50]]
        if len(errs) > 50:
            items.append(ui.tags.li(f"...and {len(errs) - 50} more errors."))
        return ui.tags.ul(*items, style="color:#b00020;")

    @output
    @render.text
    def input_protein_kegg_warning():
        return protein_kegg_warning.get()

    @output
    @render.text
    def input_ptm_status():
        data = ptm_validation.get()
        return data.get("status", "")

    @output
    @render.ui
    def input_ptm_errors():
        data = ptm_validation.get()
        errs = data.get("errors") or []
        if not errs:
            return ui.div()
        items = [ui.tags.li(err) for err in errs[:50]]
        if len(errs) > 50:
            items.append(ui.tags.li(f"...and {len(errs) - 50} more errors."))
        return ui.tags.ul(*items, style="color:#b00020;")

    @output
    @render.text
    def input_nav_lock_status():
        return nav_lock_status.get()

    @reactive.Effect
    @reactive.event(input.input_ptm_upload)
    def _process_ptm_upload():
        if not protein_validation.get().get("valid"):
            ptm_validation.set({"status": "Upload a valid protein file first.", "errors": ["Protein file not validated yet."], "valid": False})
            _update_ks_index(reset=True)
            return
        upload = input.input_ptm_upload()
        if not upload:
            ptm_validation.set({"status": "No PTM file uploaded (optional).", "errors": [], "valid": False})
            ptm_dataset.set(None)
            ptm_dataset_path.set(None)
            _write_debug_dump("user_ptm_dataset_debug.txt", {"info": "No dataset loaded"})
            _update_ks_index(reset=True)
            return
        file_info = upload[0]
        datapath = getattr(file_info, "datapath", None)
        if datapath is None and isinstance(file_info, dict):
            datapath = file_info.get("datapath")
        if not datapath:
            ptm_validation.set({"status": "Could not locate the uploaded PTM file.", "errors": [], "valid": False})
            ptm_dataset.set(None)
            ptm_dataset_path.set(None)
            _write_debug_dump("user_ptm_dataset_debug.txt", {"info": "No dataset loaded"})
            _update_ks_index(reset=True)
            return
        protein_comparisons = protein_validation.get().get("comparisons") or []
        result = validate_ptm_file(datapath, protein_comparisons)
        if result.valid:
            status = (
                f"PTM file valid. Rows: {result.summary.get('rows', 0)}, "
                f"Comparison columns: {result.summary.get('comparisons', 0)}, "
                f"Tooltip columns: {result.summary.get('tooltips', 0)}."
            )
            ptm_validation.set({"status": status, "errors": [], "valid": True})
            dataset_payload = _load_dataset(datapath)
            try:
                species_choice, _ = _resolve_species(_get_input_value(input, "input_species"))
                dataset_payload = annotate_ptm_dataset(dataset_payload, species_choice, _get_psp_map())
                dataset_payload = annotate_ptm_dataset_with_kinases(dataset_payload, species_choice, _get_psp_ks_map())
            except Exception as exc:
                print(f"Warning: PSP annotation failed: {exc}")
            ptm_dataset.set(dataset_payload)
            ptm_dataset_path.set(datapath)
            _write_debug_dump("user_ptm_dataset_debug.txt", dataset_payload)
            _update_ks_index()
        else:
            ptm_validation.set({"status": "PTM file failed validation. See errors below.", "errors": result.errors, "valid": False})
            ptm_dataset.set(None)
            ptm_dataset_path.set(None)
            _write_debug_dump("user_ptm_dataset_debug.txt", {"errors": result.errors})
            _update_ks_index(reset=True)

    @reactive.Effect
    @reactive.event(input.bookmark_selector, input.input_mode, input.input_protein_upload, input.input_ptm_upload)
    def _enforce_nav_lock():
        mode = _current_mode()
        selected = _get_input_value(input, "bookmark_selector") or "input"
        protein_ok = bool(protein_validation.get().get("valid"))
        ptm_ok = bool(ptm_validation.get().get("valid"))
        if mode == "demo":
            status_msg = "Demo mode: sample datasets loaded. Navigation unlocked." if protein_ok else "Demo mode: loading sample datasets..."
            nav_lock_status.set(status_msg)
            locked = False
            _send_custom_message_safe("toggle_nav_lock", {"locked": locked})
            return
        ready = protein_ok
        locked = not ready
        if locked and selected != "input":
            session.send_input_message("bookmark_selector", {"value": "input"})
        if ready:
            nav_lock_status.set("User mode: protein validated. Navigation unlocked. PTM optional.")
        else:
            nav_lock_status.set("User mode: upload valid Protein file to unlock other tabs. PTM optional.")
        _send_custom_message_safe("toggle_nav_lock", {"locked": locked})

    @reactive.Effect
    @reactive.event(input.input_mode)
    def _sync_mode_status():
        mode = str((_get_input_value(input, "input_mode") or "user")).lower()
        if mode == "demo":
            try:
                _load_demo_datasets()
            except Exception as exc:
                print(f"Warning: demo mode failed to load sample datasets: {exc}")
                protein_dataset.set(None)
                ptm_dataset.set(None)
                protein_validation.set({"status": "Demo mode failed to load sample protein.", "errors": [str(exc)], "valid": False, "comparisons": []})
                ptm_validation.set({"status": "Demo mode failed to load sample PTM.", "errors": [str(exc)], "valid": False})
                nav_lock_status.set("Demo mode: sample datasets failed to load.")
        else:
            protein_dataset.set(None)
            ptm_dataset.set(None)
            protein_kegg_warning.set("")
            protein_validation.set({"status": "Upload a protein file to begin.", "errors": [], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
            nav_lock_status.set("User mode: upload valid Protein file to unlock other tabs. PTM optional.")
            _update_ks_index(reset=True)

    @reactive.Effect
    @reactive.event(input.settings_gradient_preset)
    def _apply_gradient_preset():
        presets = {
            "main_default": {"neg": (2, 81, 150), "pos": (253, 179, 56)},       # swapped
            "tan_turquoise": {"neg": (64, 176, 166), "pos": (255, 190, 106)},   # swapped
            "orange_purple": {"neg": (81, 40, 136), "pos": (235, 97, 35)},      # swapped
            "green_purple": {"neg": (88, 9, 79), "pos": (41, 94, 17)},          # swapped
            "red_blue": {"neg": (47, 103, 177), "pos": (191, 44, 35)},          # swapped
            "pink_blue": {"neg": (16, 85, 154), "pos": (219, 76, 119)},         # swapped
            "yellow_pink": {"neg": (219, 16, 72), "pos": (244, 179, 1)},        # swapped
            "brown_blue": {"neg": (15, 101, 161), "pos": (106, 74, 60)},        # swapped
        }
        choice = str(_get_input_value(input, "settings_gradient_preset") or "").strip().lower()
        preset = presets.get(choice)
        if not preset:
            return
        neg_hex = _rgb_tuple_to_hex(preset["neg"])
        pos_hex = _rgb_tuple_to_hex(preset["pos"])
        try:
            session.send_input_message("settings_negative_color", {"value": neg_hex})
            session.send_input_message("settings_positive_color", {"value": pos_hex})
        except Exception:
            pass

    @output
    @render.ui
    def settings_gradient_preview():
        neg_hex = str(_get_input_value(input, "settings_negative_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]))
        pos_hex = str(_get_input_value(input, "settings_positive_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]))
        neg_val = _to_float(_get_input_value(input, "settings_max_negative"), DEFAULT_SETTINGS["max_negative"])
        pos_val = _to_float(_get_input_value(input, "settings_max_positive"), DEFAULT_SETTINGS["max_positive"])
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

    @reactive.Effect
    @reactive.event(input.input_dataset_add)
    def _add_dynamic_dataset_card():
        with reactive.isolate():
            current = list(extra_datasets.get() or [])
            next_id = max([d.get("id", 0) for d in current], default=0) + 1
            current.append({"id": next_id, "label": f"Dataset {next_id}"})
            extra_datasets.set(current)

    @reactive.Effect
    @reactive.event(input.input_dataset_remove)
    def _remove_dynamic_dataset_card():
        raw = _get_input_value(input, "input_dataset_remove")
        try:
            dataset_id = int(raw)
        except (TypeError, ValueError):
            return
        with reactive.isolate():
            remaining = [d for d in (extra_datasets.get() or []) if d.get("id") != dataset_id]
            extra_datasets.set(remaining)
        try:
            session.send_input_message("input_dataset_remove", {"value": None})
        except Exception:
            pass

    @output
    @render.ui
    def input_data_inputs_panel():
        mode = _current_mode()
        protein_ok = bool(protein_validation.get().get("valid"))
        disabled = mode == "demo"
        ptm_disabled = disabled or not protein_ok
        species_wrap_style = "opacity:0.55; pointer-events:none;" if disabled else ""
        protein_wrap_style = "opacity:0.55; pointer-events:none;" if disabled else ""
        ptm_wrap_style = "opacity:0.55; pointer-events:none;" if ptm_disabled else ""
        add_wrap_style = "opacity:0.55; pointer-events:none;" if disabled else ""
        current_species = _get_input_value(input, "input_species") or DEFAULT_SPECIES
        if current_species not in SPECIES_OPTIONS:
            current_species = DEFAULT_SPECIES

        pv = protein_validation.get()
        ptv = ptm_validation.get()
        protein_valid = bool(pv.get("valid"))
        ptm_valid = bool(ptv.get("valid"))
        error_flag = False
        for entry in (pv.get("errors") or []):
            if entry:
                error_flag = True
                break
        if not error_flag:
            for entry in (ptv.get("errors") or []):
                if entry:
                    error_flag = True
                    break
        # Also treat statuses containing "failed" or "missing" as errors
        statuses_to_check = [str(pv.get("status", "")), str(ptv.get("status", ""))]
        if not error_flag and any("failed" in s.lower() or "missing" in s.lower() for s in statuses_to_check):
            error_flag = True

        data_card = ui.card(
            {"class": "data-input-card", "style": "max-width: 540px;"},
            ui.card_header(
                ui.div(
                    {"style": "display:flex; align-items:center; justify-content:space-between; gap:12px;"},
                    ui.tags.span("Data Inputs", style="font-weight:700; color:#0b1f33; font-size:16px;"),
                    ui.div(
                        {"style": "display:flex; align-items:center; gap:8px;"},
                        ui.tags.span("Mode:", style="font-weight:600; color:#0b1f33;"),
                        ui.input_radio_buttons(
                            "input_mode",
                            None,
                            choices={"user": "User", "demo": "Demo"},
                            selected=mode or "user",
                            inline=True,
                        ),
                    ),
                )
            ),
            ui.card_body(
                ui.div(
                    {"style": species_wrap_style},
                    ui.input_select(
                        "input_species",
                        "Select Species",
                        choices=SPECIES_OPTIONS,
                        selected=current_species,
                    ),
                ),
            ),
        )

        add_card = ui.div(
            {
                "style": (
                    f"{add_wrap_style}display:flex; align-items:center; justify-content:center; "
                    "min-width: 140px; min-height: 120px; border:1px dashed #2563eb; "
                    "border-radius:12px; background:rgba(37,99,235,0.08); padding:12px;"
                )
            },
            ui.input_action_button("input_dataset_add", "+ Add dataset", width="100%"),
        )

        ptm_card = ui.card(
            {"class": "data-input-card", "style": "max-width: 400px; min-width: 350px; width: 380px;"},
            ui.card_header(
                # Always render a span for the checkmark to avoid None in the UI list; toggle visibility via style.
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.tags.span("PTM Dataset", style="font-weight:700; color:#0b1f33;"),
                    ui.tags.span(
                        "OK",
                        style="color:#16a34a; font-weight:700;" if ptm_valid else "color:transparent; font-weight:700;",
                    ),
                )
            ),
            ui.card_body(
                    ui.div(
                        {"style": ptm_wrap_style},
                        ui.input_file(
                            "input_ptm_upload",
                            "PTM Dataset Upload",
                        multiple=False,
                        accept=UPLOAD_ACCEPT_TYPES,
                    ),
                    ui.div(
                        {"style": "display:flex; align-items:center; gap:8px; margin-top:4px;"},
                        ui.input_select(
                            "input_ptm_shape",
                            None,
                            choices={"circle": "Circle", "diamond": "Diamond", "triangle": "Triangle"},
                            selected="circle",
                            width="120px",
                        ),
                    ),
                    ui.div(
                        "- The first column must contain UniProt IDs. The header name does not matter, but every row must be populated.\n"
                        "- The second column must contain PTM site positions. The header name does not matter, and all values must be positive integers.\n"
                        "- The third column is the first (and required) comparison column. Its header must start with C: and may include any descriptive text (e.g., C: TNBC vs ER+). All values must be numeric (integer or float, positive or negative).\n"
                        "- Additional comparison columns may appear after the third column. These headers must start with C:, and every row must contain a numeric value.\n"
                        "- Optional outline comparison columns may be included after the required columns. These headers must start with O: and match a C: header label (e.g., O: TNBC vs ER+). Values must be numeric or NA (NA keeps the outline black).\n"
                        "- Optional tooltip columns may be included after the required columns. These headers must start with T: (e.g., T: GOCC). Tooltip values are optional and may be blank.",
                        style="margin-top:10px; font-size:11px; color:#444; white-space:pre-line;",
                    ),
                )
            ),
        )

        protein_card = ui.card(
            {"class": "data-input-card", "style": "max-width: 400px; min-width: 350px; width: 380px;"},
            ui.card_header(
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.tags.span("Protein Dataset", style="font-weight:700; color:#0b1f33;"),
                    ui.tags.span(
                        "OK",
                        style="color:#16a34a; font-weight:700;" if protein_valid else "color:transparent; font-weight:700;",
                    ),
                )
            ),
            ui.card_body(
                ui.div(
                    {"style": protein_wrap_style},
                    ui.div(
                        {"style": "display: flex; align-items: center; gap: 8px;"},
                        ui.input_file(
                            "input_protein_upload",
                            "Protein Dataset Upload",
                            multiple=False,
                            accept=UPLOAD_ACCEPT_TYPES,
                        ),
                        ui.tags.span(
                            ui.output_text("input_protein_kegg_warning"),
                            style="color: #b00020; font-weight:600; font-size: 12px;",
                        ),
                    ),
                    ui.div(
                        "- The first column must contain UniProt IDs. The header name does not matter, but every row must have a valid UniProt ID.\n"
                        "- The second column must contain gene symbols. The header name does not matter, and all rows must be populated.\n"
                        "- The third column is the first (and required) comparison column. Its header must start with C: and may include any descriptive text (e.g., C: TNBC vs ER+). All values must be numeric (integer or float, positive or negative).\n"
                        "- Any additional comparison columns may appear after the third column. These columns must also have headers starting with C: and contain only numeric values.\n"
                        "- Optional outline comparison columns may be included after the required columns. These headers must start with O: and match a C: header label (e.g., O: TNBC vs ER+). Values must be numeric or NA (NA keeps the outline black).\n"
                        "- Optional tooltip columns may be included after the required columns. These headers must start with T: (e.g., T: GOCC). Tooltip values are optional and may be blank.",
                        style="margin-top:10px; font-size:11px; color:#444; white-space:pre-line;",
                    ),
                    ui.div(
                        "Demo mode uses bundled sample files. Switch to User mode to upload your own.",
                        style="color: #b00020; font-size: 12px;" if disabled else "display:none;",
                    ),
                )
            ),
        )

        error_block = ui.div(
            ui.output_ui("input_protein_errors"),
            ui.output_ui("input_ptm_errors"),
            style="margin-top:10px; color:#b00020; font-size:13px;",
        )

        return ui.TagList(
            ui.div(
                {"style": "display:flex; align-items:flex-start; gap:16px; width:100%;"},
                ui.div(
                    error_block,
                    style="min-width: 260px;" if error_flag else "display:none; min-width:0; width:0;",
                ),
                ui.div(
                    {"style": "display:flex; flex-direction:column; gap:16px; flex:1;"},
                    data_card,
                    ui.div(
                        {"style": "display:flex; flex-wrap:wrap; gap:16px; align-items:flex-start;"},
                        protein_card,
                        ptm_card,
                        ui.output_ui("input_extra_datasets_panel"),
                        add_card,
                    ),
                ),
            )
        )

    @output
    @render.ui
    def input_extra_datasets_panel():
        cards: List[Any] = []
        for entry in extra_datasets.get() or []:
            dataset_id = entry.get("id")
            label = entry.get("label") or f"Dataset {dataset_id}"
            cards.append(
                ui.card(
                    {"class": "data-input-card", "style": "max-width: 400px; min-width: 350px; width: 380px;"},
                    ui.card_header(
                        ui.div(
                            {"style": "display:flex; align-items:center; justify-content:space-between; gap:8px;"},
                            ui.tags.span(label, style="font-weight:700; color:#0b1f33;"),
                            ui.tags.button(
                                "Remove",
                                type="button",
                                style=(
                                    "padding:4px 10px; border:1px solid #d1d5db; background:#f9fafb; "
                                    "border-radius:8px; font-size:12px; cursor:pointer;"
                                ),
                                onclick=f"Shiny.setInputValue('input_dataset_remove', {dataset_id}, {{priority:'event'}});",
                            ),
                        )
                    ),
                    ui.card_body(
                        ui.input_text(f"input_dyn_{dataset_id}_name", "Dataset name", value=label),
                        ui.input_select(
                            f"input_dyn_{dataset_id}_type",
                            "Dataset type",
                            choices={
                                "proteomics": "Proteomics",
                                "phospho": "Phospho/PTM",
                                "transcript": "Transcriptomics",
                                "metabolomics": "Metabolomics",
                                "other": "Other",
                            },
                            selected="other",
                        ),
                        ui.input_file(
                            f"input_dyn_{dataset_id}_file",
                            "Upload dataset",
                            multiple=False,
                            accept=UPLOAD_ACCEPT_TYPES,
                        ),
                    ),
                )
            )
        return ui.TagList(*cards)

    def _color_override_from_settings(settings_override: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "negative_color": _rgb_tuple_to_list(settings_override.get("negative_color", DEFAULT_SETTINGS["negative_color"])),
            "positive_color": _rgb_tuple_to_list(settings_override.get("positive_color", DEFAULT_SETTINGS["positive_color"])),
            "max_negative": settings_override.get("max_negative", DEFAULT_SETTINGS["max_negative"]),
            "max_positive": settings_override.get("max_positive", DEFAULT_SETTINGS["max_positive"]),
        }

    def _apply_color_metadata(payload: Dict[str, Any], color_override: Dict[str, Any]) -> None:
        general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
        general_settings["negative_color"] = color_override["negative_color"]
        general_settings["positive_color"] = color_override["positive_color"]
        general_settings["max_negative"] = color_override["max_negative"]
        general_settings["max_positive"] = color_override["max_positive"]
        payload["_color_preview_override"] = color_override
        _apply_color_overrides(payload, color_override)

    def _register_bookmark(cfg: Dict[str, Any]) -> None:
        prefix = cfg["key"]
        state = bookmark_state[prefix]
        web_sort_col = reactive.Value("pathway") if prefix == "web" else None
        web_sort_dir = reactive.Value("asc") if prefix == "web" else None
        web_selected_label = reactive.Value("")
        web_selected_id = reactive.Value("") if prefix == "web" else None
        web_filter_options = reactive.Value(list(WEB_PATHWAY_SOURCES)) if prefix == "web" else None
        web_filter_selected = reactive.Value(list(WEB_PATHWAY_SOURCES)) if prefix == "web" else None
        web_filter_tick = reactive.Value(0) if prefix == "web" else None
        web_filter_refresh_evt = reactive.Value(0) if prefix == "web" else None

        def _get_custom_layout():
            # Access reactive value safely even outside a reactive context
            with reactive.isolate():
                return state["custom_layout"].get()

        def _fc_choices() -> List[str]:
            if _current_mode() in {"user", "demo"}:
                with reactive.isolate():
                    prot_data = protein_dataset.get()
                if prot_data:
                    return [h for h in (prot_data.get("headers") or []) if h.startswith("C:")]
                if _current_mode() == "user":
                    return [c for c in (protein_validation.get().get("comparisons") or []) if c]
            return [c for c in (DEFAULT_DATA.get("protein", {}).get("main_columns") or []) if c]

        is_ks_bookmark = prefix == "ks"
        ks_sort_col = reactive.Value("gene")
        ks_sort_dir = reactive.Value("asc")

        def _ks_mode_value() -> str:
            return "both"

        def _ks_allowed_types(mode: str) -> set:
            return {"in_vivo", "in_vitro"}

        def _ks_context() -> Dict[str, Any]:
            data_override = collect_data_override()
            prot_cfg = data_override.get("protein") or {}
            with reactive.isolate():
                prot_data_val = protein_dataset.get()
                ptm_data_val = ptm_dataset.get()
            # Prefer validated protein/PTM data (matches other bookmarks) when available
            if prot_data_val and prot_data_val.get("headers") and prot_data_val.get("rows"):
                prot_cfg = {
                    **prot_cfg,
                    "data_headers": prot_data_val.get("headers") or prot_cfg.get("data_headers") or [],
                    "data_rows": prot_data_val.get("rows") or prot_cfg.get("data_rows") or [],
                    "main_columns": prot_data_val.get("main_columns") or prot_cfg.get("main_columns") or [],
                    "tooltip_columns": prot_data_val.get("tooltip_columns") or prot_cfg.get("tooltip_columns") or [],
                }
            ptm_cfg_list = data_override.get("ptm") or []
            ptm_cfg = ptm_cfg_list[0] if ptm_cfg_list else {}
            if ptm_data_val and ptm_data_val.get("headers") and ptm_data_val.get("rows"):
                ptm_cfg = {
                    **ptm_cfg,
                    "data_headers": ptm_data_val.get("headers") or ptm_cfg.get("data_headers") or [],
                    "data_rows": ptm_data_val.get("rows") or ptm_cfg.get("data_rows") or [],
                    "main_columns": ptm_data_val.get("main_columns") or ptm_cfg.get("main_columns") or [],
                    "tooltip_columns": ptm_data_val.get("tooltip_columns") or ptm_cfg.get("tooltip_columns") or [],
                }
            prot_headers = prot_cfg.get("data_headers") or []
            prot_rows = prot_cfg.get("data_rows") or []
            prot_row_map: Dict[str, Dict[str, Any]] = {}
            for row in prot_rows:
                if not prot_headers:
                    continue
                uid = (row[0] if len(row) > 0 else "") or ""
                if not str(uid).strip():
                    continue
                row_dict = {h: (row[i] if i < len(row) else "") for i, h in enumerate(prot_headers)}
                prot_row_map[str(uid).strip()] = row_dict
            ptm_headers = ptm_cfg.get("data_headers") or []
            ptm_rows = ptm_cfg.get("data_rows") or []
            ptm_main = [col for col, _ in (ptm_cfg.get("main_columns") or [])]
            prot_main = prot_cfg.get("main_columns") or [h for h in prot_headers if h.startswith("C:")]
            prot_outline = prot_cfg.get("outline_columns") or _resolve_outline_columns(prot_main, prot_headers)
            ptm_outline = ptm_cfg.get("outline_columns") or _resolve_outline_columns(ptm_main, ptm_headers)
            prot_label_map: Dict[str, str] = {}
            if prot_headers:
                for uid, row_map in prot_row_map.items():
                    gene_raw = row_map.get(prot_headers[1], "") if len(prot_headers) > 1 else ""
                    prot_label_map[uid] = _clean_protein_label(gene_raw, uid)
            ctx = {
                "prot_headers": prot_headers,
                "prot_rows_by_uniprot": prot_row_map,
                "protein_main_columns": prot_main,
                "protein_outline_columns": prot_outline,
                "protein_tooltips": prot_cfg.get("tooltip_columns") or [],
                "ptm_headers": ptm_headers,
                "ptm_rows": ptm_rows,
                "ptm_main_columns": ptm_main,
                "ptm_outline_columns": ptm_outline,
                "ptm_tooltips": ptm_cfg.get("tooltip_columns") or [],
                "ptm_shape": ptm_cfg.get("shape", "Circle"),
                "ptm_type": ptm_cfg.get("type", "ptm_0"),
                "prot_label_map": prot_label_map,
            }
            with reactive.isolate():
                ks_data = ks_index.get() or _empty_ks_index()
            ctx["ptms_by_uniprot"] = ks_data.get("ptms_by_uniprot", {})
            ctx["gene_map"] = ks_data.get("prot_gene_map", {})
            ctx["ks_data"] = ks_data
            return ctx

        def _prepare_base_payload(settings_override: Dict[str, Any], main_columns: List[str]) -> Dict[str, Any]:
            settings_override = dict(settings_override)
            settings_override["ptm_max_display"] = max(10, int(settings_override.get("ptm_max_display", 10) or 10))
            settings_override["main_columns"] = main_columns or settings_override.get("main_columns", [])
            settings_override["show_arrows"] = True
            color_override = _color_override_from_settings(settings_override)
            payload = _build_blank_canvas(GLOBAL_CATALOG_INFO)
            _apply_color_metadata(payload, color_override)
            general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
            general_settings["show_background_image"] = False
            general_settings["show_groups"] = bool(settings_override.get("show_groups", False))
            general_settings["show_multi_protein_indicator"] = bool(settings_override.get("show_multi_protein_indicator", False))
            general_settings["show_arrows"] = True
            general_settings["show_text_boxes"] = True
            general_settings["debug_mode"] = bool(settings_override.get("debug_mode", False))
            general_settings["pathway_id"] = settings_override.get("pathway_id")
            general_settings["pathway_source"] = settings_override.get("pathway_source")
            general_settings["ptm_max_display"] = settings_override["ptm_max_display"]
            general_settings["prot_label_size"] = settings_override.get("prot_label_size", DEFAULT_SETTINGS["prot_label_size"])
            general_settings["ptm_label_size"] = settings_override.get("ptm_label_size", DEFAULT_SETTINGS["ptm_label_size"])
            general_settings["ptm_label_font"] = settings_override.get("ptm_label_font", DEFAULT_SETTINGS["ptm_label_font"])
            base_radius = settings_override.get("ptm_circle_radius", DEFAULT_SETTINGS["ptm_circle_radius"])
            if cfg.get("key") == "ks":
                base_radius = max(1, base_radius - 1)
            general_settings["ptm_circle_radius"] = base_radius
            general_settings["ptm_circle_spacing"] = settings_override.get("ptm_circle_spacing", DEFAULT_SETTINGS["ptm_circle_spacing"])
            general_settings["prot_outline_width"] = settings_override.get("prot_outline_width", DEFAULT_SETTINGS.get("prot_outline_width", 1))
            general_settings["ptm_outline_width"] = settings_override.get("ptm_outline_width", DEFAULT_SETTINGS.get("ptm_outline_width", 1))
            general_settings["negative_color"] = color_override["negative_color"]
            general_settings["positive_color"] = color_override["positive_color"]
            general_settings["max_negative"] = color_override["max_negative"]
            general_settings["max_positive"] = color_override["max_positive"]
            general_settings["main_columns"] = settings_override["main_columns"]
            general_settings["protein_tooltip_columns"] = settings_override.get("protein_tooltip_columns", [])
            general_settings["species"] = settings_override.get("species")
            general_settings["species_code"] = settings_override.get("species_code")
            payload["_kegg_preview"] = {
                "show": False,
                "opacity": DEFAULT_BG_OPACITY,
                "offset_x": DEFAULT_BG_OFFSET_X,
                "offset_y": DEFAULT_BG_OFFSET_Y,
                "scale": DEFAULT_BG_SCALE,
            }
            payload["_box_preview"] = {"y_stretch": DEFAULT_BOX_Y_STRETCH}
            payload["_active_mode"] = settings_override.get("mode", cfg["mode"])
            payload["_full_width_canvas"] = True
            payload["_custom_layout_applied"] = False
            payload["_global_protein_catalog"] = dict(GLOBAL_CATALOG_INFO)
            return payload

        def _make_protein_entry(
            uniprot_id: str,
            protein_row: Dict[str, Any],
            ctx: Dict[str, Any],
            settings_override: Dict[str, Any],
            ptms: Dict[str, Any],
        ) -> Dict[str, Any]:
            label = (
                ctx.get("gene_map", {}).get(uniprot_id)
                or ctx.get("prot_label_map", {}).get(uniprot_id)
                or protein_row.get(ctx.get("prot_headers", [None, None])[1] if ctx.get("prot_headers") else "", "")
                or uniprot_id
            )
            label = _clean_protein_label(label, uniprot_id)
            entry: Dict[str, Any] = {
                "label": label,
                "gene_symbol": label,
                "label_color": [0, 0, 0],
                "transcriptomic_color": [],
                "annotations": "",
                "tooltip": "",
                "tooltip_html": "",
                "PTMs": ptms,
            }
            tooltip_plain: List[str] = []
            tooltip_html: List[str] = []
            for col in ctx.get("protein_tooltips", []):
                val = str(protein_row.get(col, "") or "").strip()
                if not val:
                    continue
                tooltip_plain.append(f"{col}: {val}")
                tooltip_html.append(f"<strong>{html.escape(col)}:</strong> {html.escape(val)}")
            entry["tooltip"] = "\n".join(tooltip_plain)
            entry["tooltip_html"] = "<br/>".join(tooltip_html)
            ann_values = []
            for col in ctx.get("protein_tooltips", []):
                val = str(protein_row.get(col, "") or "").strip()
                if val:
                    ann_values.append(f'"{val}"')
            entry["annotations"] = ",".join(ann_values)
            neg = settings_override.get("negative_color", DEFAULT_SETTINGS["negative_color"])
            pos = settings_override.get("positive_color", DEFAULT_SETTINGS["positive_color"])
            max_neg = settings_override.get("max_negative", DEFAULT_SETTINGS["max_negative"])
            max_pos = settings_override.get("max_positive", DEFAULT_SETTINGS["max_positive"])
            outline_cols = ctx.get("protein_outline_columns", [])
            for idx, col in enumerate(ctx.get("protein_main_columns", []), 1):
                raw = protein_row.get(col, "")
                try:
                    fc_val = float(raw)
                except (TypeError, ValueError):
                    fc_val = None
                entry[f"fold_change_{idx}"] = fc_val
                entry[f"fc_color_{idx}"] = _gradient_color_from_fold(fc_val, neg, pos, max_neg, max_pos)
                outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                outline_raw = protein_row.get(outline_col, "") if outline_col else ""
                try:
                    outline_val = float(outline_raw)
                except (TypeError, ValueError):
                    outline_val = None
                entry[f"outline_fold_change_{idx}"] = outline_val
                entry[f"outline_color_{idx}"] = _gradient_color_from_fold(outline_val, neg, pos, max_neg, max_pos) if outline_val is not None else [0, 0, 0]
            if not ctx.get("protein_main_columns"):
                entry["fold_change_1"] = None
                entry["fc_color_1"] = [128, 128, 128]
                entry["outline_fold_change_1"] = None
                entry["outline_color_1"] = [0, 0, 0]
            return entry

        def _make_ptm_entry(
            row_map: Dict[str, Any],
            ptm_key: str,
            pos_key: str,
            box_x: float,
            box_y: float,
            box_width: float,
            box_height: float,
            ctx: Dict[str, Any],
            settings_override: Dict[str, Any],
            display_label: Optional[str] = None,
        ) -> Dict[str, Any]:
            spacing = settings_override.get("ptm_circle_spacing", DEFAULT_SETTINGS["ptm_circle_spacing"])
            shape_x, shape_y = _compute_ptm_position(pos_key, box_x, box_y, box_width, box_height, spacing)
            label_dx, label_dy, label_center = PTM_LABEL_DEFAULTS.get(pos_key, (0, -6, "center"))
            tooltip_plain: List[str] = []
            tooltip_html: List[str] = []
            for col in ctx.get("ptm_tooltips", []):
                val = str(row_map.get(col, "") or "").strip()
                if not val:
                    continue
                tooltip_plain.append(f"{col}: {val}")
                tooltip_html.append(f"<strong>{html.escape(col)}:</strong> {html.escape(val)}")
            reg_site_val = str(row_map.get("PSP: regulatory_site", "") or "").strip()
            symbol_class = str(row_map.get("Phosphosite_Classification", "") or "").strip()
            symbol_text = ""
            annotated_flag = ""
            if reg_site_val:
                symbol_text = "+"
                annotated_flag = "+"
            elif symbol_class and symbol_class.lower() != "none":
                symbol_text = symbol_class[:3]
                annotated_flag = "+"
            neg = settings_override.get("negative_color", DEFAULT_SETTINGS["negative_color"])
            pos = settings_override.get("positive_color", DEFAULT_SETTINGS["positive_color"])
            max_neg = settings_override.get("max_negative", DEFAULT_SETTINGS["max_negative"])
            max_pos = settings_override.get("max_positive", DEFAULT_SETTINGS["max_positive"])
            payload: Dict[str, Any] = {
                "ptm_type": ctx.get("ptm_type", "ptm_0"),
                "uniprot_id": row_map.get(ctx.get("ptm_headers", [None, None])[0] if ctx.get("ptm_headers") else "", ""),
                "fc_color_1": [128, 128, 128],
                "outline_color_1": [0, 0, 0],
                "outline_fold_change_1": None,
                "label": display_label or str(row_map.get(ctx.get("ptm_headers", [None, None])[1] if ctx.get("ptm_headers") else "", "") or ""),
                "label_color": [0, 0, 0],
                "symbol_type": symbol_class,
                "annotated": annotated_flag,
                "shape": ctx.get("ptm_shape", "Circle"),
                "symbol": symbol_text,
                "symbol_color": [0, 0, 0],
                "symbol_font": "Arial",
                "symbol_size": 6,
                "tooltip": "\n".join(tooltip_plain),
                "tooltip_html": "<br/>".join(tooltip_html),
                "ptm_position": pos_key,
                "shape_x": float(shape_x),
                "shape_y": float(shape_y),
                "label_x": float(shape_x + label_dx),
                "label_y": float(shape_y + label_dy),
                "label_centering": label_center,
                "symbol_x": float(shape_x),
                "symbol_y": float(shape_y),
            }
            outline_cols = ctx.get("ptm_outline_columns", [])
            for idx, col in enumerate(ctx.get("ptm_main_columns", []), 1):
                raw = row_map.get(col, "")
                try:
                    fc_val = float(raw)
                except (TypeError, ValueError):
                    fc_val = None
                payload[f"fold_change_{idx}"] = fc_val
                payload[f"fc_color_{idx}"] = _gradient_color_from_fold(fc_val, neg, pos, max_neg, max_pos)
                outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                outline_raw = row_map.get(outline_col, "") if outline_col else ""
                try:
                    outline_val = float(outline_raw)
                except (TypeError, ValueError):
                    outline_val = None
                payload[f"outline_fold_change_{idx}"] = outline_val
                payload[f"outline_color_{idx}"] = _gradient_color_from_fold(outline_val, neg, pos, max_neg, max_pos) if outline_val is not None else [0, 0, 0]
            return payload

        def _first_fc_value(row_map: Dict[str, Any], columns: List[str]) -> Optional[float]:
            for col in columns:
                raw = row_map.get(col, "")
                try:
                    return float(raw)
                except (TypeError, ValueError):
                    continue
            return None

        def _resolve_gene_label(uniprot_id: str, info: Dict[str, Any], ctx: Dict[str, Any]) -> str:
            gene = info.get("gene")
            if not gene:
                gene = ctx.get("prot_label_map", {}).get(uniprot_id)
            if not gene:
                row = ctx.get("prot_rows_by_uniprot", {}).get(uniprot_id, {})
                headers = ctx.get("prot_headers", [])
                if len(headers) > 1:
                    gene = row.get(headers[1], "")
            return _clean_protein_label(gene, uniprot_id)

        def _ks_thresholds() -> Tuple[float, float]:
            max_pos = _to_float(_get_input_value(input, "settings_max_positive"), DEFAULT_SETTINGS["max_positive"])
            max_neg = _to_float(_get_input_value(input, "settings_max_negative"), DEFAULT_SETTINGS["max_negative"])
            return max_pos, max_neg

        def _build_view_for_kinase(
            kinase_id: str, mode: str, ctx: Dict[str, Any], settings_override: Dict[str, Any]
        ) -> Dict[str, Any]:
            ks_data = ctx.get("ks_data") or _empty_ks_index()
            kin_entry = ks_data.get("kinases", {}).get(kinase_id)
            allowed_types = _ks_allowed_types(mode)
            if not kin_entry:
                state["status"].set(f"No kinase data found for {kinase_id}.")
                return _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            sub_keys: List[str] = []
            for sub_key in kin_entry.get("substrates", []):
                rel_types = (
                    ks_data.get("substrates", {})
                    .get(sub_key, {})
                    .get("kinases", {})
                    .get(kinase_id, {})
                    .get("types", set())
                )
                if rel_types and rel_types.intersection(allowed_types):
                    sub_keys.append(sub_key)
            if not sub_keys:
                state["status"].set("No substrates matched this kinase and evidence filter.")
                return _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            sub_keys = sorted(sub_keys)
            grouped: Dict[str, List[str]] = {}
            for sk in sub_keys:
                if ":" in sk:
                    uid, site_part = sk.split(":", 1)
                else:
                    uid, site_part = sk, ""
                grouped.setdefault(uid, []).append(sk)
            grouped_items = sorted(grouped.items(), key=lambda kv: kv[0])
            protbox_width = 46
            protbox_height = 17
            spacing_y = KS_VERTICAL_SPACING
            base_x_left = 140
            base_x_right = 380
            total_span = spacing_y * (len(grouped_items) - 1) if grouped_items else 0
            start_y = 200 - (total_span / 2.0)
            kinase_y = start_y + (total_span / 2.0)
            payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            protein_data: Dict[str, Any] = {}
            protboxes: List[Dict[str, Any]] = []
            arrows: List[Dict[str, Any]] = []
            kinase_ptms: Dict[str, Any] = {}
            kinase_ptm_rows = ctx.get("ptms_by_uniprot", {}).get(kinase_id, [])
            for idx, row_map in enumerate(kinase_ptm_rows[: len(PTM_POSITION_PRIORITY)]):
                pos_key = PTM_POSITION_PRIORITY[idx]
                ptm_key = f"{kinase_id}_{str(row_map.get(ctx.get('ptm_headers', [None, None])[1] if ctx.get('ptm_headers') else '', '')).strip()}"
                kinase_ptms[ptm_key] = _make_ptm_entry(
                    row_map,
                    ptm_key,
                    pos_key,
                    base_x_left,
                    kinase_y,
                    protbox_width,
                    protbox_height,
                    ctx,
                    settings_override,
                )
            kinase_row = ctx.get("prot_rows_by_uniprot", {}).get(kinase_id, {})
            protein_data[kinase_id] = _make_protein_entry(kinase_id, kinase_row, ctx, settings_override, kinase_ptms)
            kin_label = protein_data[kinase_id].get("label") or kin_entry.get("gene") or kinase_id
            kin_pb_id = f"{prefix}_kinase"
            kin_tooltip = protein_data[kinase_id].get("tooltip") or ""
            protboxes.append(
                {
                    "protbox_id": kin_pb_id,
                    "proteins": [kinase_id],
                    "backup_label": kin_label,
                    "x": base_x_left,
                    "y": kinase_y,
                    "width": protbox_width,
                    "height": protbox_height,
                    "tooltip": kin_tooltip,
                }
            )
            for idx, (sub_uid, sub_list) in enumerate(grouped_items):
                y_pos = start_y + idx * spacing_y
                ptm_entries: Dict[str, Any] = {}
                for ptm_idx, sub_key in enumerate(sub_list[: len(PTM_POSITION_PRIORITY)]):
                    row_map = ks_data.get("substrates", {}).get(sub_key, {}).get("row", {}) or {}
                    if not row_map:
                        continue
                    sub_info = ks_data.get("substrates", {}).get(sub_key, {}) or {}
                    site_label = sub_info.get("site_label") or sub_info.get("site") or str(
                        row_map.get(ctx.get("ptm_headers", [None, None])[1] if ctx.get("ptm_headers") else "", "")
                    )
                    pos_key = PTM_POSITION_PRIORITY[ptm_idx]
                    ptm_key = f"{sub_uid}_{site_label}"
                    ptm_entries[ptm_key] = _make_ptm_entry(
                        row_map,
                        ptm_key,
                        pos_key,
                        base_x_right,
                        y_pos,
                        protbox_width,
                        protbox_height,
                        ctx,
                        settings_override,
                        display_label=site_label,
                    )
                sub_row = ctx.get("prot_rows_by_uniprot", {}).get(sub_uid, {})
                protein_data[sub_uid] = _make_protein_entry(sub_uid, sub_row, ctx, settings_override, ptm_entries)
                sub_label = protein_data[sub_uid].get("label") or ks_data.get("substrates", {}).get(sub_list[0], {}).get("gene") or sub_uid
                pb_id = f"{prefix}_sub_{idx + 1}"
                protboxes.append(
                    {
                        "protbox_id": pb_id,
                        "proteins": [sub_uid],
                        "backup_label": sub_label,
                        "x": base_x_right,
                        "y": y_pos,
                        "width": protbox_width,
                        "height": protbox_height,
                    }
                )
                arrows.append(
                    {
                        "protbox_id_1": kin_pb_id,
                        "protbox_id_2": pb_id,
                        "protbox_id_1_side": "East",
                        "protbox_id_2_side": "West",
                        "line": "arrow",
                        "type": "",
                    }
                )
            payload["protein_data"] = protein_data
            payload["protbox_data"] = protboxes
            payload["arrows"] = arrows
            state["status"].set(
                f"Loaded kinase {kin_label} with {len(grouped_items)} substrate protein(s) and {len(sub_keys)} substrate PTM(s)."
            )
            return payload

        def _build_view_for_ptm(
            sub_key: str, mode: str, ctx: Dict[str, Any], settings_override: Dict[str, Any]
        ) -> Dict[str, Any]:
            ks_data = ctx.get("ks_data") or _empty_ks_index()
            sub_entry = ks_data.get("substrates", {}).get(sub_key)
            allowed_types = _ks_allowed_types(mode)
            if not sub_entry or not (sub_entry.get("types", set()) & allowed_types):
                state["status"].set("No PTM entry matched the evidence filter.")
                return _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            kin_ids = []
            for kin_id, rel in (sub_entry.get("kinases") or {}).items():
                rel_types = rel.get("types", set())
                if rel_types and rel_types.intersection(allowed_types):
                    kin_ids.append(kin_id)
            if not kin_ids:
                state["status"].set("No known kinases found for this PTM with the selected evidence filter.")
                return _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            kin_ids = sorted(set(kin_ids))
            protbox_width = 46
            protbox_height = 17
            spacing_y = KS_VERTICAL_SPACING
            base_x_left = 140
            base_x_right = 380
            total_span = spacing_y * (len(kin_ids) - 1) if kin_ids else 0
            start_y = 200 - (total_span / 2.0)
            ptm_y = start_y + (total_span / 2.0)
            payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            protein_data: Dict[str, Any] = {}
            protboxes: List[Dict[str, Any]] = []
            arrows: List[Dict[str, Any]] = []
            ptm_row = sub_entry.get("row", {}) or {}
            ptm_label = sub_entry.get("site_label") or sub_entry.get("site") or str(
                ptm_row.get(ctx.get("ptm_headers", [None, None])[1] if ctx.get("ptm_headers") else "", "")
            )
            ptm_key = f"{sub_entry.get('uniprot', '')}_{ptm_label}"
            ptm_entries = {
                ptm_key: _make_ptm_entry(
                    ptm_row,
                    ptm_key,
                    "W1",
                    base_x_right,
                    ptm_y,
                    protbox_width,
                    protbox_height,
                    ctx,
                    settings_override,
                    display_label=ptm_label,
                )
            }
            sub_uid = sub_entry.get("uniprot", "")
            sub_row = ctx.get("prot_rows_by_uniprot", {}).get(sub_uid, {})
            protein_data[sub_uid] = _make_protein_entry(sub_uid, sub_row, ctx, settings_override, ptm_entries)
            sub_label = protein_data[sub_uid].get("label") or sub_entry.get("gene") or sub_uid
            ptm_pb_id = f"{prefix}_ptm"
            protboxes.append(
                {
                    "protbox_id": ptm_pb_id,
                    "proteins": [sub_uid],
                    "backup_label": sub_label,
                    "x": base_x_right,
                    "y": ptm_y,
                    "width": protbox_width,
                    "height": protbox_height,
                }
            )
            for idx, kin_id in enumerate(kin_ids):
                kin_y = start_y + idx * spacing_y
                kin_ptms: Dict[str, Any] = {}
                kin_ptm_rows = ctx.get("ptms_by_uniprot", {}).get(kin_id, [])
                for ptm_idx, row_map in enumerate(kin_ptm_rows[: len(PTM_POSITION_PRIORITY)]):
                    pos_key = PTM_POSITION_PRIORITY[ptm_idx]
                    kin_ptm_key = f"{kin_id}_{str(row_map.get(ctx.get('ptm_headers', [None, None])[1] if ctx.get('ptm_headers') else '', '')).strip()}"
                    kin_ptms[kin_ptm_key] = _make_ptm_entry(
                        row_map,
                        kin_ptm_key,
                        pos_key,
                        base_x_left,
                        kin_y,
                        protbox_width,
                        protbox_height,
                        ctx,
                        settings_override,
                    )
                kin_row = ctx.get("prot_rows_by_uniprot", {}).get(kin_id, {})
                protein_data[kin_id] = _make_protein_entry(kin_id, kin_row, ctx, settings_override, kin_ptms)
                kin_label = protein_data[kin_id].get("label") or ks_data.get("kinases", {}).get(kin_id, {}).get("gene") or kin_id
                kin_pb_id = f"{prefix}_kin_{idx + 1}"
                protboxes.append(
                    {
                        "protbox_id": kin_pb_id,
                        "proteins": [kin_id],
                        "backup_label": kin_label,
                        "x": base_x_left,
                        "y": kin_y,
                        "width": protbox_width,
                        "height": protbox_height,
                    }
                )
                arrows.append(
                    {
                        "protbox_id_1": kin_pb_id,
                        "protbox_id_2": ptm_pb_id,
                        "protbox_id_1_side": "East",
                        "protbox_id_2_side": "West",
                        "line": "arrow",
                        "type": "",
                    }
                )
            payload["protein_data"] = protein_data
            payload["protbox_data"] = protboxes
            payload["arrows"] = arrows
            state["status"].set(
                f"Loaded PTM {sub_label} with {len(kin_ids)} known kinase(s) ({mode})."
            )
            return payload

        def build_json():
            settings_override = collect_settings(input, cfg)
            if settings_override.get("pathway_source") == "wikipathways" and not settings_override.get("pathway_id"):
                state["json"].set(None)
                state["status"].set("Select a WikiPathways entry, then click Load Pathway.")
                return
            if cfg.get("key") == "web":
                selected_id = ""
                if web_selected_id is not None:
                    try:
                        selected_id = web_selected_id.get() or ""
                    except Exception:
                        selected_id = ""
                current_id = str(settings_override.get("pathway_id") or "").strip()
                if not selected_id or not current_id or selected_id.strip().lower() != current_id.lower():
                    state["json"].set(None)
                    state["status"].set("Select a pathway from the table, then click Load Pathway.")
                    return
            if is_ks_bookmark:
                ctx = _ks_context()
                if ctx.get("protein_main_columns"):
                    settings_override["main_columns"] = list(ctx.get("protein_main_columns", []))
                if ctx.get("protein_tooltips"):
                    settings_override["protein_tooltip_columns"] = list(ctx.get("protein_tooltips", []))
                mode = _ks_mode_value()
                selection = _get_input_value(input, _prefixed_id(prefix, "ks_choice"))
                if not ctx.get("prot_headers") or not ctx.get("ptm_headers"):
                    payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
                    state["json"].set(payload)
                    state["status"].set("Upload protein and PTM datasets to use Kinase Substrates.")
                    return
                if not selection:
                    payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
                    state["json"].set(payload)
                    state["status"].set("Select a kinase or substrate, then click Load.")
                    return
                selection_str = str(selection)
                if selection_str.startswith("kinase|"):
                    payload = _build_view_for_kinase(selection_str.split("|", 1)[1], mode, ctx, settings_override)
                elif selection_str.startswith("ptm|"):
                    payload = _build_view_for_ptm(selection_str.split("|", 1)[1], mode, ctx, settings_override)
                else:
                    payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
                    state["status"].set("Choose a kinase or substrate option, then click Load.")
                state["json"].set(payload)
                return
            data_override = collect_data_override()
            if data_override:
                protein_cfg = data_override.get("protein", {})
                if protein_cfg.get("main_columns"):
                    settings_override["main_columns"] = list(protein_cfg["main_columns"])
                with reactive.isolate():
                    prot_data = protein_dataset.get()
                prot_headers = list(prot_data.get("headers") or []) if prot_data else []
                if len(prot_headers) >= 2:
                    settings_override["prot_uniprot_column"] = prot_headers[0]
                    settings_override["gene_name_column"] = prot_headers[1]
                    tooltip_cols = [prot_headers[1], prot_headers[0]]
                    tooltip_cols.extend([h for h in prot_headers if h.startswith("T:")])
                    settings_override["protein_tooltip_columns"] = tooltip_cols
                    # Set ID column for matching based on source: KEGG uses KEGG_Gene_ID; others use UniProt (first column)
                    if settings_override.get("pathway_source", "").lower() == "kegg":
                        settings_override["hsa_id_column"] = "KEGG_Gene_ID"
                    else:
                        settings_override["hsa_id_column"] = prot_headers[0]
            # Align column names in settings with the uploaded datasets when in User mode
            if data_override and (_get_input_value(input, "input_mode") or "user") == "user":
                with reactive.isolate():
                    prot_data = protein_dataset.get()
                if prot_data:
                    headers = prot_data.get("headers") or []
                    if len(headers) >= 2:
                        settings_override["prot_uniprot_column"] = headers[0]
                        settings_override["gene_name_column"] = headers[1]
                        if settings_override.get("pathway_source", "").lower() == "kegg":
                            settings_override["hsa_id_column"] = "KEGG_Gene_ID"
                        else:
                            settings_override["hsa_id_column"] = headers[0]
                # KEGG column already forced to KEGG_Gene_ID in collect_settings; keep consistent
            catalog_info = dict(GLOBAL_CATALOG_INFO)
            color_override = _color_override_from_settings(settings_override)
            if cfg.get("start_blank"):
                payload = _build_blank_canvas(GLOBAL_CATALOG_INFO)
                payload = dict(payload)
                seed_payload = None
                # Only seed from data override for non-figure modes; figure clear should start empty
                if data_override and cfg.get("key") != "figure":
                    try:
                        seed_payload = get_default_json(
                            data_override=data_override,
                            settings_override=settings_override,
                            skip_disk_write=True,
                            debug_write=debug_var,
                        )
                    except Exception as exc:
                        print(f"Warning: failed to build seed payload for blank canvas: {exc}")
                merged_settings = dict(DEFAULT_SETTINGS)
                merged_settings.update(settings_override)
                merged_settings["pathway_id"] = settings_override.get("pathway_id", DEFAULT_SETTINGS["pathway_id"])
                merged_settings["pathway_source"] = settings_override.get("pathway_source", DEFAULT_SETTINGS["pathway_source"])
                merged_settings["mode"] = settings_override.get("mode", cfg["mode"])
                merged_settings["show_background_image"] = False
                merged_settings["show_groups"] = bool(settings_override.get("show_groups", False))
                merged_settings["show_multi_protein_indicator"] = bool(settings_override.get("show_multi_protein_indicator", False))
                merged_settings["show_arrows"] = bool(settings_override.get("show_arrows", True))
                merged_settings["show_text_boxes"] = bool(settings_override.get("show_text_boxes", True))
                merged_settings["use_original_protbox_size"] = bool(settings_override.get("use_original_protbox_size", DEFAULT_SETTINGS.get("use_original_protbox_size", False)))
                merged_settings["ptm_max_display"] = max(
                    1, int(settings_override.get("ptm_max_display", DEFAULT_SETTINGS["ptm_max_display"]) or DEFAULT_SETTINGS["ptm_max_display"])
                )
                merged_settings["prot_label_size"] = settings_override.get("prot_label_size", DEFAULT_SETTINGS["prot_label_size"])
                merged_settings["ptm_label_size"] = settings_override.get("ptm_label_size", DEFAULT_SETTINGS["ptm_label_size"])
                merged_settings["ptm_label_font"] = settings_override.get("ptm_label_font", DEFAULT_SETTINGS["ptm_label_font"])
                base_radius = settings_override.get("ptm_circle_radius", DEFAULT_SETTINGS["ptm_circle_radius"])
                merged_settings["ptm_circle_radius"] = max(1, base_radius)
                merged_settings["ptm_circle_spacing"] = settings_override.get("ptm_circle_spacing", DEFAULT_SETTINGS["ptm_circle_spacing"])
                merged_settings["prot_outline_width"] = settings_override.get("prot_outline_width", DEFAULT_SETTINGS.get("prot_outline_width", 1))
                merged_settings["ptm_outline_width"] = settings_override.get("ptm_outline_width", DEFAULT_SETTINGS.get("ptm_outline_width", 1))
                merged_settings["main_columns"] = settings_override.get("main_columns", [])
                merged_settings["protein_tooltip_columns"] = settings_override.get("protein_tooltip_columns", DEFAULT_SETTINGS.get("protein_tooltip_columns", []))
                merged_settings["species"] = settings_override.get("species")
                merged_settings["species_code"] = settings_override.get("species_code")
                if merged_settings.get("pathway_source", "").lower() == "kegg":
                    merged_settings["hsa_id_column"] = "KEGG_Gene_ID"
                else:
                    # Default to first protein column when available; otherwise fall back to UniProt_ID
                    if prot_data and prot_data.get("headers"):
                        merged_settings["hsa_id_column"] = prot_data["headers"][0]
                    else:
                        merged_settings["hsa_id_column"] = "Uniprot_ID"
                general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
                general_settings.clear()
                general_settings.update(merged_settings)
                if seed_payload:
                    payload["protein_data"] = seed_payload.get("protein_data", {}) or {}
                    payload["compound_data"] = seed_payload.get("compound_data", []) or []
                    payload["text_data"] = seed_payload.get("text_data", []) or []
                _apply_color_metadata(payload, color_override)
                payload["_bookmark_key"] = cfg["key"]
                payload["_persist_token"] = time.time()
                payload["_global_protein_catalog"] = dict(GLOBAL_CATALOG_INFO)
                payload["_kegg_preview"] = {
                    "show": False,
                    "opacity": DEFAULT_BG_OPACITY,
                    "offset_x": DEFAULT_BG_OFFSET_X,
                    "offset_y": DEFAULT_BG_OFFSET_Y,
                    "scale": DEFAULT_BG_SCALE,
                }
                payload["_box_preview"] = {"y_stretch": DEFAULT_BOX_Y_STRETCH}
                payload["_active_mode"] = settings_override.get("mode", cfg["mode"])
                if cfg["key"] in {"figure", "web"}:
                    payload["_full_width_canvas"] = True
                layout_override = _get_custom_layout()
                if layout_override:
                    _apply_custom_layout(payload, layout_override)
                    payload["_custom_layout_applied"] = True
                else:
                    payload["_custom_layout_applied"] = False
                state["json"].set(payload)
                state["status"].set("Blank canvas ready. Use the viewer to add protboxes and elements.")
                return
            try:
                payload = get_default_json(
                    data_override=data_override if data_override else None,
                    settings_override=settings_override,
                    skip_disk_write=True,
                    debug_write=debug_var,
                )
                show_bg_flag = bool(settings_override.get("show_background_image", False))
                payload, _ = _attach_kegg_background_image(payload, force=show_bg_flag)
                payload = dict(payload)
                payload["_bookmark_key"] = cfg["key"]
                payload["_persist_token"] = time.time()
                payload["_global_protein_catalog"] = dict(GLOBAL_CATALOG_INFO)
                _apply_color_metadata(payload, color_override)
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
                active_mode = settings_override.get("mode", cfg["mode"])
                payload["_active_mode"] = active_mode
                if cfg["key"] in {"figure", "web"}:
                    payload["_full_width_canvas"] = True
                if str(active_mode).lower() == "analysis":
                    payload["arrows"] = []
                if cfg["key"] == "simple":
                    for group in payload.get("groups", []) or []:
                        if isinstance(group, dict):
                            group["show_box"] = False
                layout_override = _get_custom_layout()
                if layout_override:
                    _apply_custom_layout(payload, layout_override)
                    payload["_custom_layout_applied"] = True
                else:
                    payload["_custom_layout_applied"] = False
                state["json"].set(payload)
                state["status"].set(f"{cfg['label']} pathway generated.")
            except Exception as exc:
                state["status"].set(f"Error generating JSON: {exc}")
                state["json"].set(None)

        if is_ks_bookmark:
            def _ks_table_rows() -> List[Dict[str, Any]]:
                ctx = _ks_context()
                ks_data = ctx.get("ks_data") or _empty_ks_index()
                entity_mode = str(_get_input_value(input, _prefixed_id(prefix, "ks_entity_mode")) or "substrate").lower()
                allowed = {"in_vivo", "in_vitro"}
                fc_choices_local = _fc_choices()
                fc_idx = state["fc_index"].get() or 1
                selected_fc = fc_choices_local[fc_idx - 1] if fc_idx >= 1 and fc_idx <= len(fc_choices_local) else None
                search_term = str(_get_input_value(input, _prefixed_id(prefix, "ks_filter_regex")) or "").strip()
                search_re = None
                if search_term:
                    try:
                        search_re = re.compile(search_term, re.IGNORECASE)
                    except re.error:
                        search_re = None
                max_pos, max_neg = _ks_thresholds()
                evidence_filter = str(_get_input_value(input, _prefixed_id(prefix, "ks_filter_evidence")) or "both").lower()
                fc_op = str(_get_input_value(input, _prefixed_id(prefix, "ks_filter_fc_op")) or "")
                fc_val_raw = _get_input_value(input, _prefixed_id(prefix, "ks_filter_fc_val"))
                try:
                    fc_val = float(fc_val_raw) if fc_val_raw not in (None, "", False) else None
                except (TypeError, ValueError):
                    fc_val = None
                reg_only = bool(_get_input_value(input, _prefixed_id(prefix, "ks_filter_reg_only")))
                def _is_sig(fc: Optional[float]) -> bool:
                    return fc is not None and (fc >= max_pos or fc <= max_neg)

                def _evidence_matches(label: str) -> bool:
                    if evidence_filter == "both":
                        return True
                    if evidence_filter == "in_vivo":
                        return label in {"in vivo", "both"}
                    if evidence_filter == "in_vitro":
                        return label in {"in vitro", "both"}
                    return True

                def _fc_matches(value: Optional[float]) -> bool:
                    if fc_val is None or not fc_op:
                        return True
                    if value is None:
                        return False
                    if fc_op == "gt":
                        return value > fc_val
                    if fc_op == "lt":
                        return value < fc_val
                    if fc_op == "ge":
                        return value >= fc_val
                    if fc_op == "le":
                        return value <= fc_val
                    if fc_op == "eq":
                        return value == fc_val
                    if fc_op == "ne":
                        return value != fc_val
                    return True
                rows: List[Dict[str, Any]] = []
                if entity_mode == "kinase":
                    for kin_id, info in (ks_data.get("kinases") or {}).items():
                        types = info.get("types", set())
                        if not types:
                            continue
                        gene = _resolve_gene_label(kin_id, info, ctx)
                        prot_row = ctx.get("prot_rows_by_uniprot", {}).get(kin_id, {})
                        prot_fc_cols = [selected_fc] if selected_fc else ctx.get("protein_main_columns", [])
                        fc_val = _first_fc_value(prot_row, prot_fc_cols)
                        sig_count = 0
                        sub_count = 0
                        for sub_key in info.get("substrates", []):
                            sub_entry = (ks_data.get("substrates") or {}).get(sub_key, {})
                            if sub_entry.get("types", set()) & allowed:
                                sub_count += 1
                                sub_fc = _first_fc_value(sub_entry.get("row", {}), ctx.get("ptm_main_columns", []))
                                if _is_sig(sub_fc):
                                    sig_count += 1
                        evid_label = "both" if types >= {"in_vivo", "in_vitro"} else ("in vivo" if "in_vivo" in types else ("in vitro" if "in_vitro" in types else ""))
                        if not _evidence_matches(evid_label):
                            continue
                        if not _fc_matches(fc_val):
                            continue
                        has_reg = ""
                        for row_map in ctx.get("ptms_by_uniprot", {}).get(kin_id, []):
                            if str(row_map.get("PSP: regulatory_site", "")).strip():
                                has_reg = "+"
                                break
                        reg_display = has_reg
                        display = gene
                        if search_re:
                            if not (search_re.search(display) or search_re.search(str(gene)) or search_re.search(kin_id)):
                                continue
                        rows.append(
                            {
                                "id": f"kinase|{kin_id}",
                                "gene": gene,
                                "label": display,
                                "fc": fc_val,
                                "reg": reg_display,
                                "sig_pct": (sig_count / sub_count * 100.0) if sub_count else None,
                                "evidence": evid_label,
                                "count": sub_count,
                            }
                        )
                else:
                    for sub_key, info in (ks_data.get("substrates") or {}).items():
                        types = info.get("types", set())
                        if not types:
                            continue
                        sub_uid = info.get("uniprot", "")
                        gene = _resolve_gene_label(sub_uid, info, ctx)
                        site = info.get("site_label") or info.get("site") or ""
                        ptm_fc_cols = []
                        if selected_fc and selected_fc in (ctx.get("ptm_main_columns") or []):
                            ptm_fc_cols = [selected_fc]
                        else:
                            ptm_fc_cols = ctx.get("ptm_main_columns", [])
                        ptm_fc_val = _first_fc_value(info.get("row", {}), ptm_fc_cols)
                        kin_count = len(info.get("kinases") or {})
                        reg_val = str(info.get("row", {}).get("PSP: regulatory_site", "")).strip()
                        evid_label = "both" if types >= {"in_vivo", "in_vitro"} else ("in vivo" if "in_vivo" in types else ("in vitro" if "in_vitro" in types else ""))
                        if not _evidence_matches(evid_label):
                            continue
                        if not _fc_matches(ptm_fc_val):
                            continue
                        if reg_only:
                            if not reg_val:
                                continue
                        display = f"{gene} {site}".strip()
                        if search_re:
                            if not (search_re.search(display) or search_re.search(str(gene)) or search_re.search(sub_uid) or search_re.search(site)):
                                continue
                        rows.append(
                            {
                                "id": f"ptm|{sub_key}",
                                "gene": gene,
                                "site": site,
                                "fc": ptm_fc_val,
                                "reg": reg_val,
                                "evidence": evid_label,
                                "count": kin_count,
                            }
                        )
                sort_col = ks_sort_col.get() or "gene"
                sort_dir = ks_sort_dir.get() or "asc"
                reverse = sort_dir == "desc"
                def sort_key(entry):
                    val = entry.get(sort_col)
                    if isinstance(val, str):
                        return val.lower()
                    return val if val is not None else -float("inf")
                rows.sort(key=sort_key, reverse=reverse)
                return rows

            @output(id=_prefixed_id(prefix, "ks_table"))
            @render.ui
            def ks_table():
                rows = _ks_table_rows()
                entity_mode = str(_get_input_value(input, _prefixed_id(prefix, "ks_entity_mode")) or "substrate").lower()
                sort_col = ks_sort_col.get() or "gene"
                sort_dir = ks_sort_dir.get() or "asc"
                if not rows:
                    return ui.div({"style": "font-size:12px; color:#555; margin-top:6px;"}, "No matches. Upload data and adjust search or filters.")
                def header(label: str, col_key: str) -> Any:
                    arrow = "▲" if sort_col == col_key and sort_dir == "asc" else ("▼" if sort_col == col_key else "")
                    return ui.tags.th(
                        {
                            "style": "cursor:pointer; user-select:none;",
                            "onclick": f"Shiny.setInputValue('{_prefixed_id(prefix, 'ks_sort')}', {{col:'{col_key}', ts:Date.now()}}, {{priority:'event'}});"
                        },
                        ui.tags.span(label),
                        ui.tags.span(f" {arrow}" if arrow else ""),
                    )
                headers = []
                if entity_mode == "kinase":
                    headers = [
                        header("Protein", "gene"),
                        header("Regulatory site", "reg"),
                        header("FC", "fc"),
                        header("Evidence", "evidence"),
                        header("Significant PTMs (%)", "sig_pct"),
                        header("Substrates", "count"),
                    ]
                else:
                    headers = [
                        header("Protein", "gene"),
                        header("Site", "site"),
                        header("FC", "fc"),
                        header("Regulatory site", "reg"),
                        header("Evidence", "evidence"),
                        header("Kinases", "count"),
                    ]
                body_rows = []
                selected_id = str(_get_input_value(input, _prefixed_id(prefix, "ks_choice")) or "")
                for entry in rows:
                    row_id = entry.get("id", "")
                    row_classes = []
                    if selected_id and selected_id == row_id:
                        row_classes.append("table-active")
                    row_attrs = {
                        "style": "cursor:pointer;",
                        "onclick": (
                            f"Shiny.setInputValue('{_prefixed_id(prefix, 'ks_choice')}', '{row_id}', {{priority:'event'}});"
                        ),
                    }
                    if row_classes:
                        row_attrs["class"] = " ".join(row_classes)
                    if entity_mode == "kinase":
                        cells = [
                            ui.tags.td(entry.get("gene", "")),
                            ui.tags.td(entry.get("reg", "")),
                            ui.tags.td("N/A" if entry.get("fc") is None else f"{entry['fc']:.3g}"),
                            ui.tags.td(entry.get("evidence", "")),
                            ui.tags.td("N/A" if entry.get("sig_pct") is None else f"{entry['sig_pct']:.1f}%"),
                            ui.tags.td(str(entry.get("count", 0))),
                        ]
                    else:
                        cells = [
                            ui.tags.td(entry.get("gene", "")),
                            ui.tags.td(entry.get("site", "")),
                            ui.tags.td("" if entry.get("fc") is None else f"{entry['fc']:.3g}"),
                            ui.tags.td(entry.get("reg", "")),
                            ui.tags.td(entry.get("evidence", "")),
                            ui.tags.td(str(entry.get("count", 0))),
                        ]
                    body_rows.append(ui.tags.tr(*cells, **row_attrs))
                return ui.tags.table(
                    {"class": "table table-sm table-hover", "style": "margin-top:6px;"},
                    ui.tags.thead(ui.tags.tr(*headers)),
                    ui.tags.tbody(*body_rows),
                )

            @reactive.Effect
            @reactive.event(getattr(input, _prefixed_id(prefix, "ks_sort")))
            def _update_ks_sort():
                payload = _get_input_value(input, _prefixed_id(prefix, "ks_sort")) or {}
                col = str(payload.get("col") or "")
                if col not in {"gene", "site", "fc", "count", "reg", "sig_pct", "evidence"}:
                    return
                current_col = ks_sort_col.get()
                current_dir = ks_sort_dir.get()
                if col == current_col:
                    ks_sort_dir.set("desc" if current_dir == "asc" else "asc")
                else:
                    ks_sort_col.set(col)
                    ks_sort_dir.set("asc")

            @reactive.Effect
            @reactive.event(getattr(input, _prefixed_id(prefix, "ks_filter_refresh_evt")))
            def _reset_ks_filters():
                session.send_input_message(_prefixed_id(prefix, "ks_filter_regex"), {"value": ""})
                session.send_input_message(_prefixed_id(prefix, "ks_filter_reg_only"), {"value": False})
                session.send_input_message(_prefixed_id(prefix, "ks_filter_evidence"), {"value": "both"})
                session.send_input_message(_prefixed_id(prefix, "ks_filter_fc_op"), {"value": ""})
                session.send_input_message(_prefixed_id(prefix, "ks_filter_fc_val"), {"value": ""})

            @reactive.Effect
            @reactive.event(getattr(input, _prefixed_id(prefix, "ks_choice")))
            def _build_kinase_substrate_view():
                build_json()

            @reactive.Effect
            @reactive.event(getattr(input, "ks_spawn_ptm_kinases"))
            def _append_ks_ptm_kinases():
                payload = _get_input_value(input, "ks_spawn_ptm_kinases") or {}
                target_pb_id = str(payload.get("protbox_id") or "").strip()
                ptm_key_raw = str(payload.get("ptm_key") or "").strip()
                sub_uniprot = str(payload.get("uniprot") or "").strip()
                if not (target_pb_id and ptm_key_raw and sub_uniprot):
                    return
                ctx = _ks_context()
                ks_data = ctx.get("ks_data") or {}
                mode = _ks_mode_value()
                allowed_types = _ks_allowed_types(mode)
                def _match_sub_key() -> Optional[str]:
                    suffix = ptm_key_raw.split("_", 1)[-1] if "_" in ptm_key_raw else ptm_key_raw
                    for key, entry in (ks_data.get("substrates") or {}).items():
                        if entry.get("uniprot") != sub_uniprot:
                            continue
                        if suffix and (entry.get("site_label") == suffix or entry.get("site") == suffix):
                            return key
                    for key, entry in (ks_data.get("substrates") or {}).items():
                        if entry.get("uniprot") == sub_uniprot:
                            return key
                    return None
                sub_key = None
                if ptm_key_raw in (ks_data.get("substrates") or {}):
                    sub_key = ptm_key_raw
                else:
                    sub_key = _match_sub_key()
                sub_entry = (ks_data.get("substrates") or {}).get(sub_key) if sub_key else None
                if not sub_entry:
                    state["status"].set("No matching substrate found for this PTM.")
                    return
                kin_ids: List[str] = []
                for kin_id, rel in (sub_entry.get("kinases") or {}).items():
                    rel_types = rel.get("types", set())
                    if rel_types and rel_types.intersection(allowed_types):
                        kin_ids.append(kin_id)
                kin_ids = sorted(set(kin_ids))
                if not kin_ids:
                    state["status"].set("No known kinases for this substrate.")
                    return
                settings_override = collect_settings(input, cfg)
                if ctx.get("protein_main_columns"):
                    settings_override["main_columns"] = list(ctx.get("protein_main_columns", []))
                if ctx.get("protein_tooltips"):
                    settings_override["protein_tooltip_columns"] = list(ctx.get("protein_tooltips", []))
                if ctx.get("ptm_main_columns"):
                    settings_override["ptm_main_columns"] = list(ctx.get("ptm_main_columns", []))
                current_payload = state["json"].get()
                if not current_payload:
                    current_payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
                protein_data = current_payload.setdefault("protein_data", {})
                protboxes = current_payload.setdefault("protbox_data", [])
                arrows = current_payload.setdefault("arrows", [])
                target_pb = None
                for pb in protboxes:
                    if str(pb.get("protbox_id")) == target_pb_id:
                        target_pb = pb
                        break
                if not target_pb:
                    target_pb = {
                        "protbox_id": target_pb_id,
                        "proteins": [sub_uniprot],
                        "backup_label": sub_entry.get("gene") or sub_uniprot,
                        "x": 380,
                        "y": 200,
                        "width": 46,
                        "height": 17,
                    }
                    protboxes.append(target_pb)
                pb_width = float(target_pb.get("width") or 46)
                pb_height = float(target_pb.get("height") or 17)
                pb_x = float(target_pb.get("x") or 380)
                pb_y = float(target_pb.get("y") or 200)
                protbox_by_uniprot: Dict[str, Dict[str, Any]] = {}
                for pb in protboxes:
                    for uid in pb.get("proteins") or []:
                        protbox_by_uniprot[str(uid)] = pb
                next_idx = len(protboxes) + 1
                spacing = KS_VERTICAL_SPACING
                total_span = spacing * (len(kin_ids) - 1) if kin_ids else 0
                start_y = pb_y - (total_span / 2.0)
                base_x_left = pb_x - 200
                ptm_headers = ctx.get("ptm_headers", [])
                for idx, kin_id in enumerate(kin_ids):
                    y_pos = start_y + idx * spacing
                    kin_pb = protbox_by_uniprot.get(kin_id)
                    if not kin_pb:
                        kin_ptms: Dict[str, Any] = {}
                        kin_ptm_rows = ctx.get("ptms_by_uniprot", {}).get(kin_id, [])
                        for ptm_idx, row_map in enumerate(kin_ptm_rows[: len(PTM_POSITION_PRIORITY)]):
                            pos_key = PTM_POSITION_PRIORITY[ptm_idx]
                            site_val = str(row_map.get(ptm_headers[1] if len(ptm_headers) > 1 else "", "")).strip()
                            kin_ptm_key = f"{kin_id}_{site_val}"
                            kin_ptms[kin_ptm_key] = _make_ptm_entry(
                                row_map,
                                kin_ptm_key,
                                pos_key,
                                base_x_left,
                                y_pos,
                                pb_width,
                                pb_height,
                                ctx,
                                settings_override,
                            )
                        kin_row = ctx.get("prot_rows_by_uniprot", {}).get(kin_id, {})
                        protein_data.setdefault(kin_id, _make_protein_entry(kin_id, kin_row, ctx, settings_override, kin_ptms))
                        kin_label = protein_data[kin_id].get("label") or kin_id
                        kin_pb = {
                            "protbox_id": f"{prefix}_ks_add_{next_idx}",
                            "proteins": [kin_id],
                            "backup_label": kin_label,
                            "x": base_x_left,
                            "y": y_pos,
                            "width": pb_width,
                            "height": pb_height,
                        }
                        next_idx += 1
                        protboxes.append(kin_pb)
                        protbox_by_uniprot[kin_id] = kin_pb
                    target_id = target_pb.get("protbox_id")
                    kin_pb_id = kin_pb.get("protbox_id")
                    if target_id and kin_pb_id:
                        exists = False
                        for ar in arrows:
                            if (
                                str(ar.get("protbox_id_1")) == str(kin_pb_id) and str(ar.get("protbox_id_2")) == str(target_id)
                            ) or (
                                str(ar.get("protbox_id_1")) == str(target_id) and str(ar.get("protbox_id_2")) == str(kin_pb_id)
                            ):
                                exists = True
                                break
                        if not exists:
                            arrows.append(
                                {
                                    "protbox_id_1": kin_pb_id,
                                    "protbox_id_2": target_id,
                                    "protbox_id_1_side": "East",
                                    "protbox_id_2_side": "West",
                                    "line": "arrow",
                                    "type": "",
                                }
                            )
                current_payload["protbox_data"] = protboxes
                current_payload["protein_data"] = protein_data
                current_payload["arrows"] = arrows
                state["json"].set(current_payload)
                state["status"].set(f"Added {len(kin_ids)} kinase(s) for this substrate.")

        if cfg.get("show_search"):
            @output(id=_prefixed_id(prefix, "pathway_search_results"))
            @render.ui
            def pathway_search_results():
                if cfg.get("key") == "web":
                    return ui.div({"class": "pathway-search-empty"}, "Use the table below to pick a pathway.")
                query = (str(_get_input_value(input, _prefixed_id(prefix, "pathway_id"))) or "").strip()
                matches: List[Dict[str, str]] = []
                if cfg.get("key") == "web":
                    species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
                    species_label = species_info.get("label") or species_key
                    species_full = species_info.get("species") or species_label
                    wp_options = _load_wikipathways_catalog(species_label, fallback=species_full)
                    if not wp_options:
                        msg = f"No WikiPathways catalogue found for {species_label or 'selected species'}."
                        return ui.div({"class": "pathway-search-empty"}, msg)
                    def _wp_match_entry(opt: Dict[str, Any]) -> Optional[Dict[str, str]]:
                        path_id = str(opt.get("id") or "").strip().upper()
                        if not path_id:
                            return None
                        label = opt.get("label") or f"{path_id} | {opt.get('name', '')}".strip()
                        return {"value": path_id, "label": label}

                    if not query:
                        for opt in wp_options[:WIKIPATHWAYS_MAX_MATCHES]:
                            entry = _wp_match_entry(opt)
                            if entry:
                                matches.append(entry)
                    else:
                        try:
                            pattern = re.compile(query, re.IGNORECASE)
                        except re.error:
                            return ui.div({"class": "pathway-search-error"}, "Invalid regex pattern")
                        for opt in wp_options:
                            search_target = f"{opt.get('id', '')} | {opt.get('name', '')} | {opt.get('species', '')}"
                            if pattern.search(search_target):
                                entry = _wp_match_entry(opt)
                                if entry:
                                    matches.append(entry)
                            if len(matches) >= WIKIPATHWAYS_MAX_MATCHES:
                                break
                    if not matches:
                        msg = f"No pathways matched your search for {species_label or 'selection'}."
                        return ui.div({"class": "pathway-search-empty"}, msg)
                else:
                    if not KEGG_PATHWAY_OPTIONS:
                        return ui.div({"class": "pathway-search-empty"}, "No KEGG catalogue found.")
                    species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
                    species_code = species_info["code"]
                    if not query:
                        # Show all options when empty; list is scrollable in UI.
                        for opt in KEGG_PATHWAY_OPTIONS:
                            species_id = f"{species_code}{opt['digits']}"
                            matches.append({"value": species_id, "label": f"{species_id} | {opt['name']}"})
                    else:
                        try:
                            pattern = re.compile(query, re.IGNORECASE)
                        except re.error:
                            return ui.div({"class": "pathway-search-error"}, "Invalid regex pattern")
                        for opt in KEGG_PATHWAY_OPTIONS:
                            species_id = f"{species_code}{opt['digits']}"
                            search_target = f"{species_id} | {opt['raw_id']} | {opt['name']}"
                            if pattern.search(search_target):
                                matches.append({"value": species_id, "label": f"{species_id} | {opt['name']}"})
                            if len(matches) >= KEGG_PATHWAY_MAX_MATCHES:
                                break
                    if not matches:
                        return ui.div({"class": "pathway-search-empty"}, "No pathways matched your search.")
                active_value = (str(_get_input_value(input, _prefixed_id(prefix, "pathway_id"))) or "").strip()
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
                                "onclick": f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_id_choice')}', this.dataset.value, {{priority: 'event'}})",
                            },
                            opt["label"],
                        )
                    )
                return ui.TagList(
                    ui.tags.ul({"class": "pathway-results"}, *items),
                    ui.tags.script(
                        f"""
                        (function(){{
                            const wrapper = document.getElementById('{_prefixed_id(prefix, "pathway_search_wrapper")}');
                            if (!wrapper) return;
                            const results = wrapper.querySelector('.pathway-results');
                            if (!results) return;
                            results.style.display = '';
                            if (wrapper.dataset.bound === '1') {{
                                return;
                            }}
                            const input = wrapper.querySelector('input');
                            function hideResults(){{
                                if (results){{
                                    results.style.display = 'none';
                                }}
                            }}
                            function handleDocClick(ev){{
                                if (!wrapper.contains(ev.target) && results){{
                                    hideResults();
                                }}
                            }}
                            function showResults(){{
                                if (results){{
                                    results.style.display = '';
                                }}
                            }}
                            document.addEventListener('click', handleDocClick, true);
                            if (input){{
                                input.addEventListener('focus', showResults);
                                input.addEventListener('input', showResults);
                                input.addEventListener('blur', () => {{
                                    setTimeout(hideResults, 120);
                                }});
                            }}
                            if (results){{
                                results.addEventListener('click', hideResults);
                            }}
                            wrapper.dataset.bound = '1';
                        }})();
                        """
                    ),
                )

            @reactive.Effect
            @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_id_choice")))
            def _apply_pathway_choice():
                selection = _get_input_value(input, _prefixed_id(prefix, "pathway_id_choice"))
                value = selection
                source_val = _get_input_value(input, _prefixed_id(prefix, "pathway_source_choice"))
                name_val = ""
                if isinstance(selection, dict):
                    value = selection.get("value") or selection.get("id")
                    if selection.get("source"):
                        source_val = selection.get("source")
                    name_val = selection.get("name") or ""
                if value:
                    session.send_input_message(_prefixed_id(prefix, "pathway_id"), {"value": value})
                    if web_selected_id is not None and cfg.get("key") == "web":
                        web_selected_id.set(str(value))
                if source_val and cfg.get("key") == "web":
                    session.send_input_message(_prefixed_id(prefix, "pathway_source_choice"), {"value": source_val})
                if cfg.get("key") == "web":
                    label_src = str(source_val or "").upper()
                    label_val = str(value or "").upper()
                    label = f"{label_src}: {label_val}" if value else ""
                    if name_val:
                        label = f"{label} | {name_val}" if label else name_val
                    if web_selected_label is not None:
                        web_selected_label.set(label)
                state["status"].set("Pathway selected. Click Load Pathway to load with current settings.")

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_id")))
                def _sync_web_selection_input():
                    current = (str(_get_input_value(input, _prefixed_id(prefix, "pathway_id"))) or "").strip()
                    selected = web_selected_id.get() if web_selected_id is not None else ""
                    if selected and current and selected.strip().lower() == current.lower():
                        return
                    if web_selected_id is not None:
                        web_selected_id.set("")
                    if web_selected_label is not None:
                        web_selected_label.set("")

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_table_sort")))
                def _handle_web_sort():
                    payload = _get_input_value(input, _prefixed_id(prefix, "pathway_table_sort")) or {}
                    col = str(payload.get("col") or "").lower()
                    if col not in {"source", "pathway"}:
                        return
                    current_col = web_sort_col.get()
                    current_dir = web_sort_dir.get()
                    if col == current_col:
                        web_sort_dir.set("desc" if current_dir == "asc" else "asc")
                    else:
                        web_sort_col.set(col)
                        web_sort_dir.set("asc")

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_filter_sources")))
                def _sync_web_filter_sources():
                    raw = _get_input_value(input, _prefixed_id(prefix, "pathway_filter_sources"))
                    selected: List[str] = []
                    canonical = {s.lower(): s for s in WEB_PATHWAY_SOURCES}
                    if isinstance(raw, (list, tuple, set)):
                        selected = [canonical.get(str(r).lower(), str(r)) for r in raw if r]
                    elif raw:
                        selected = [canonical.get(str(raw).lower(), str(raw))]
                    if web_filter_options is not None:
                        with reactive.isolate():
                            options = web_filter_options.get() or []
                        selected = [s for s in selected if s in options]
                    if not selected:
                        selected = list(WEB_PATHWAY_SOURCES)
                    if web_filter_selected is not None:
                        web_filter_selected.set(selected)
                    if web_filter_tick is not None:
                        try:
                            current = web_filter_tick.get() or 0
                        except Exception:
                            current = 0
                        web_filter_tick.set(current + 1)

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_filter_refresh_evt")))
                def _sync_web_filter_refresh():
                    # Treat refresh as a manual bump to the tick used by the table render
                    if web_filter_tick is not None:
                        try:
                            current = web_filter_tick.get() or 0
                        except Exception:
                            current = 0
                        web_filter_tick.set(current + 1)
                    if web_filter_refresh_evt is not None:
                        try:
                            current = web_filter_refresh_evt.get() or 0
                        except Exception:
                            current = 0
                        web_filter_refresh_evt.set(current + 1)

            if cfg.get("key") == "web":
                @output(id=_prefixed_id(prefix, "pathway_table"))
                @render.ui
                def pathway_table():
                    if web_filter_tick is not None:
                        # Depend on tick so table refreshes when checkboxes change
                        _ = web_filter_tick.get()
                    species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
                    species_label = species_info.get("label") or species_key
                    species_full = species_info.get("species") or species_label
                    options = _load_wikipathways_catalog(species_label, fallback=species_full)
                    species_code = species_info.get("code", "")
                    kegg_rows: List[Dict[str, str]] = []
                    for opt in KEGG_PATHWAY_OPTIONS:
                        species_id = f"{species_code}{opt['digits']}"
                        kegg_rows.append({"id": species_id, "pathway": opt["name"], "source": "KEGG"})
                    wp_rows: List[Dict[str, str]] = []
                    if not options:
                        if not kegg_rows:
                            return ui.div({"class": "pathway-search-empty"}, "No pathway catalogue found for this species.")
                    else:
                        wp_rows = [
                            {
                                "id": str(opt.get("id") or "").strip(),
                                "pathway": str(opt.get("name") or "").strip(),
                                "source": "WikiPathways",
                                "species": str(opt.get("species") or ""),
                            }
                            for opt in options
                            if opt.get("id") and opt.get("name")
                        ]
                    query = (str(_get_input_value(input, _prefixed_id(prefix, "pathway_id"))) or "").strip()
                    combined_rows: List[Dict[str, str]] = []
                    pattern = None
                    if query:
                        try:
                            pattern = re.compile(query, re.IGNORECASE)
                        except re.error:
                            return ui.div({"class": "pathway-search-error"}, "Invalid regex pattern")
                    for opt in wp_rows:
                        search_target = f"{opt['id']} | {opt['pathway']} | {opt.get('species', '')} | {opt['source']}"
                        if pattern and not pattern.search(search_target):
                            continue
                        combined_rows.append({"id": opt["id"], "pathway": opt["pathway"], "source": opt["source"]})
                    for opt in kegg_rows:
                        search_target = f"{opt['id']} | {opt['pathway']} | {opt['source']}"
                        if pattern and not pattern.search(search_target):
                            continue
                        combined_rows.append(opt)
                    available_sources = list(WEB_PATHWAY_SOURCES)
                    if web_filter_options is not None:
                        web_filter_options.set(available_sources)
                    raw_selected = _get_input_value(input, _prefixed_id(prefix, "pathway_filter_sources"))
                    selected_sources: Optional[List[str]] = None
                    if isinstance(raw_selected, (list, tuple, set)):
                        selected_sources = [str(s) for s in raw_selected if s]
                    elif isinstance(raw_selected, str) and raw_selected:
                        selected_sources = [str(raw_selected)]
                    if selected_sources is None and web_filter_selected is not None:
                        try:
                            selected_sources = web_filter_selected.get()
                        except Exception:
                            selected_sources = None
                    if selected_sources is None:
                        selected_sources = available_sources
                    selected_sources = [s for s in (selected_sources or []) if s in available_sources]
                    if not selected_sources:
                        selected_sources = available_sources
                        if web_filter_selected is not None:
                            web_filter_selected.set(selected_sources)
                    if selected_sources:
                        selected_lower = {s.lower() for s in selected_sources}
                        combined_rows = [row for row in combined_rows if str(row.get("source", "")).lower() in selected_lower]
                    if not combined_rows:
                        return ui.div({"class": "pathway-search-empty"}, "No pathways matched your search.")
                    sort_col = (web_sort_col.get() or "pathway") if web_sort_col else "pathway"
                    sort_dir = (web_sort_dir.get() or "asc") if web_sort_dir else "asc"
                    reverse = sort_dir == "desc"
                    if sort_col == "source":
                        combined_rows.sort(
                            key=lambda r: (
                                str(r.get("source", "")).lower(),
                                str(r.get("pathway", "")).lower(),
                            ),
                            reverse=reverse,
                        )
                    else:
                        combined_rows.sort(
                            key=lambda r: (
                                str(r.get("pathway", "")).lower(),
                                str(r.get("source", "")).lower(),
                            ),
                            reverse=reverse,
                        )
                    active_value = (str(_get_input_value(input, _prefixed_id(prefix, "pathway_id"))) or "").strip().upper()
                    # Update selection label to include name if still selected
                    for entry in combined_rows:
                        if active_value and active_value == entry["id"]:
                            label_src = entry.get("source", "").upper()
                            label_val = entry["id"].upper()
                            name_val = entry.get("pathway") or ""
                            label = f"{label_src}: {label_val}"
                            if name_val:
                                label = f"{label} | {name_val}"
                            web_selected_label.set(label)
                            break
                    table_rows: List[Any] = []
                    for entry in combined_rows:
                        classes = []
                        if active_value and active_value == entry["id"]:
                            classes.append("table-active")
                        table_rows.append(
                            ui.tags.tr(
                                {
                                    "data-value": entry["id"],
                                    "data-source": entry.get("source", ""),
                                    "data-name": entry.get("pathway", ""),
                                    "class": " ".join(classes),
                                    "onclick": (
                                        "(() => {"
                                        "const data={value:this.dataset.value,source:this.dataset.source||'wikipathways',name:this.dataset.name||''};"
                                        f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_id_choice')}', data, {{priority:'event'}});"
                                        f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_source_choice')}', data.source, {{priority:'event'}});"
                                        "})();"
                                    ),
                                },
                                ui.tags.td(entry["source"]),
                                ui.tags.td(
                                    ui.tags.div(entry["pathway"]),
                                    ui.tags.small({"style": "color:#6b7280;"}, entry["id"]),
                                ),
                            )
                        )
                    source_arrow = "▲" if sort_col == "source" and sort_dir == "asc" else ("▼" if sort_col == "source" else "")
                    pathway_arrow = "▲" if sort_col == "pathway" and sort_dir == "asc" else ("▼" if sort_col == "pathway" else "")
                    table = ui.tags.table(
                        {"class": "table table-sm table-hover pathway-table"},
                        ui.tags.thead(
                            ui.tags.tr(
                                ui.tags.th(
                                    {"style": "cursor:pointer;", "onclick": f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_table_sort')}', {{col:'source', ts:Date.now()}}, {{priority:'event'}});"},
                                    ui.tags.span("Source"),
                                    ui.tags.span(f" {source_arrow}" if source_arrow else ""),
                                ),
                                ui.tags.th(
                                    {"style": "cursor:pointer;", "onclick": f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_table_sort')}', {{col:'pathway', ts:Date.now()}}, {{priority:'event'}});"},
                                    ui.tags.span("Pathway"),
                                    ui.tags.span(f" {pathway_arrow}" if pathway_arrow else ""),
                                ),
                            )
                        ),
                        ui.tags.tbody(*table_rows),
                    )
                    return ui.TagList(
                        table,
                        ui.tags.script(
                            """
                            (function(){
                                const table = document.querySelector('.pathway-table');
                                if(!table) return;
                                table.addEventListener('click', (ev) => {
                                    const row = ev.target.closest('tr');
                                    if(!row || !row.dataset.value) return;
                                    const rows = table.querySelectorAll('tr');
                                    rows.forEach(r => r.classList.remove('table-active'));
                                    row.classList.add('table-active');
                                });
                            })();
                            """
                        ),
                    )

                @output(id=_prefixed_id(prefix, "selected_pathway_label"))
                @render.text
                def selected_pathway_label():
                    label = web_selected_label.get() if web_selected_label else ""
                    if not label:
                        raw_value = (str(_get_input_value(input, _prefixed_id(prefix, "pathway_id"))) or "").strip()
                        if raw_value:
                            if cfg.get("key") == "web":
                                source_hint = str(_get_input_value(input, _prefixed_id(prefix, "pathway_source_choice")) or "wikipathways").strip().upper() or "WIKIPATHWAYS"
                            else:
                                source_hint = "KEGG"
                            label = f"{source_hint}: {raw_value}"
                    return label or "Selected: none"

        @output(id=_prefixed_id(prefix, "gradient_preview"))
        @render.ui
        def gradient_preview():
            neg_hex = str(_get_input_value(input, "settings_negative_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]))
            pos_hex = str(_get_input_value(input, "settings_positive_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]))
            neg_val = _to_float(_get_input_value(input, "settings_max_negative"), DEFAULT_SETTINGS["max_negative"])
            pos_val = _to_float(_get_input_value(input, "settings_max_positive"), DEFAULT_SETTINGS["max_positive"])
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

        @output(id=_prefixed_id(prefix, "fc_selector"))
        @render.ui
        def fc_selector():
            choices = _fc_choices()
            if not choices:
                return ui.div({"style": "font-size:12px; color:#666;"}, "No comparison columns available yet.")
            current_idx = state["fc_index"].get() or 1
            if current_idx > len(choices):
                current_idx = 1
                state["fc_index"].set(current_idx)
            select_choices = {str(i + 1): label for i, label in enumerate(choices)}
            return ui.input_select(
                _prefixed_id(prefix, "fc_select"),
                "Fold-change column",
                choices=select_choices,
                selected=str(current_idx),
            )

        @reactive.Effect
        @reactive.event(getattr(input, _prefixed_id(prefix, "fc_select")))
        def _sync_fc_choice():
            val = _get_input_value(input, _prefixed_id(prefix, "fc_select"))
            choices = _fc_choices()
            if not choices:
                return
            try:
                idx = int(val)
            except Exception:
                idx = 1
            idx = max(1, min(idx, len(choices)))
            state["fc_index"].set(idx)
            # Update persist token so client can reuse layout but refresh colors for new FC
            current = state["json"].get()
            if current:
                updated = dict(current)
                updated["_active_fc_index"] = idx
                updated["_persist_token"] = time.time()
                state["json"].set(updated)

        @reactive.Effect
        @reactive.event(
            input.settings_negative_color,
            input.settings_positive_color,
            input.settings_max_negative,
            input.settings_max_positive,
        )
        def _sync_colors_from_settings():
            current = state["json"].get()
            if not current:
                return
            try:
                settings_override = collect_settings(input, cfg)
                color_override = _color_override_from_settings(settings_override)
                payload = copy.deepcopy(current)
                _apply_color_metadata(payload, color_override)
                payload["_persist_token"] = time.time()
                state["json"].set(payload)
            except Exception as exc:
                print(f"Warning: failed to sync colors from settings for bookmark {prefix}: {exc}")

        @reactive.Effect
        @reactive.event(input.settings_debug_mode)
        def _sync_debug_mode():
            current = state["json"].get()
            if not current:
                return
            try:
                debug_mode = _to_bool(_get_input_value(input, "settings_debug_mode"), False)
                payload = copy.deepcopy(current)
                general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
                general_settings["debug_mode"] = debug_mode
                payload["_persist_token"] = time.time()
                state["json"].set(payload)
            except Exception as exc:
                print(f"Warning: failed to sync debug mode for bookmark {prefix}: {exc}")

        @reactive.Effect
        @reactive.event(getattr(input, _prefixed_id(prefix, "spawn_protboxes")))
        def _spawn_protboxes():
            if cfg.get("key") != "figure":
                return
            raw = _get_input_value(input, _prefixed_id(prefix, "spawn_protboxes_ids")) or ""
            ids: List[str] = []
            for part in re.split(r"[,\n]+", str(raw)):
                cleaned = part.strip()
                if cleaned:
                    ids.append(cleaned)
            if not ids:
                state["status"].set("No Uniprot IDs provided.")
                return
            current = state["json"].get()
            if not current:
                state["status"].set("Canvas not ready. Click Clear Canvas first.")
                return
            payload = copy.deepcopy(current)
            prot_data = payload.setdefault("protein_data", {})
            protboxes = payload.setdefault("protbox_data", [])
            existing_ids = {str(pb.get("protbox_id")) for pb in protboxes}
            catalog_info = payload.get("_global_protein_catalog") or GLOBAL_CATALOG_INFO

            def _catalog_lookup(uniprot: str) -> Dict[str, Any]:
                if not isinstance(catalog_info, dict):
                    return {}
                cat_map = catalog_info.get("protein_catalog")
                if not isinstance(cat_map, dict):
                    path = catalog_info.get("path")
                    if path and os.path.exists(path):
                        try:
                            with open(path, "r", encoding="utf-8") as fh:
                                cat_payload = json.load(fh)
                            cat_map = cat_payload.get("protein_catalog", {})
                            if isinstance(cat_map, dict):
                                catalog_info["protein_catalog"] = cat_map
                        except Exception:
                            cat_map = {}
                if isinstance(cat_map, dict) and uniprot in cat_map:
                    return cat_map.get(uniprot) or {}
                return {}

            def _resolve_protein_entry(uniprot: str) -> Dict[str, Any]:
                existing_entry = payload.get("protein_data", {}).get(uniprot) if isinstance(payload, dict) else None
                cat_entry = _catalog_lookup(uniprot)
                base_entry: Dict[str, Any] = {}
                if isinstance(existing_entry, dict):
                    base_entry = copy.deepcopy(existing_entry)
                elif isinstance(cat_entry, dict) and cat_entry:
                    base_entry = copy.deepcopy(cat_entry)
                if not base_entry:
                    base_entry = {"label": uniprot, "gene_symbol": uniprot, "PTMs": {}}
                else:
                    # If existing PTMs are empty but catalog has them, merge them in
                    base_ptms = base_entry.get("PTMs") if isinstance(base_entry, dict) else {}
                    cat_ptms = cat_entry.get("PTMs") if isinstance(cat_entry, dict) else {}
                    if (not base_ptms or len(base_ptms) == 0) and isinstance(cat_ptms, dict) and cat_ptms:
                        base_entry["PTMs"] = copy.deepcopy(cat_ptms)
                return base_entry

            start_x, start_y = 80, 120
            layout_mode = str(_get_input_value(input, _prefixed_id(prefix, "spawn_layout_mode")) or "grid").lower()
            grid_step_x = _to_float(_get_input_value(input, _prefixed_id(prefix, "spawn_grid_x")), 75) or 75
            grid_step_y = _to_float(_get_input_value(input, _prefixed_id(prefix, "spawn_grid_y")), 45) or 45
            grid_max_per_row = max(1, int(_to_float(_get_input_value(input, _prefixed_id(prefix, "spawn_grid_row")), 5) or 5))
            added = 0
            general_settings = payload.get("general_data", {}).get("settings", {}) if isinstance(payload, dict) else {}
            ptm_max_display = int(general_settings.get("ptm_max_display", DEFAULT_SETTINGS["ptm_max_display"]) or DEFAULT_SETTINGS["ptm_max_display"])
            # Layout helpers
            def _append_protbox(uni: str, idx: int, x: float, y: float) -> None:
                nonlocal added
                pb_id = f"user_pb_{len(protboxes) + added + 1}"
                while pb_id in existing_ids:
                    added += 1
                    pb_id = f"user_pb_{len(protboxes) + added + 1}"
                prot_entry = _resolve_protein_entry(uni)
                gene = prot_entry.get("gene_symbol") or prot_entry.get("label") or uni
                ptms_map = prot_entry.get("PTMs") if isinstance(prot_entry, dict) else None
                if isinstance(ptms_map, dict):
                    ordered_ptms = list(ptms_map.items())
                    if ptm_max_display > 0:
                        ordered_ptms = ordered_ptms[:ptm_max_display]
                    prot_entry["PTMs"] = {k: v for k, v in ordered_ptms}
                else:
                    prot_entry["PTMs"] = {}
                prot_data[uni] = prot_entry
                protboxes.append(
                    {
                        "protbox_id": pb_id,
                        "proteins": [uni],
                        "backup_label": gene,
                        "x": x,
                        "y": y,
                        "width": 46,
                        "height": 17,
                    }
                )
                existing_ids.add(pb_id)
                added += 1

            if layout_mode == "concentric" and ids:
                center_id = ids[0]
                ring_ids = ids[1:]
                center_x, center_y = start_x, start_y
                pb_w, pb_h = 46, 17
                _append_protbox(center_id, 0, center_x - pb_w * 0.5, center_y - pb_h * 0.5)
                # Determine ring ids from tooltip if requested
                use_tooltip = bool(_get_input_value(input, _prefixed_id(prefix, "spawn_conc_use_tooltip")))
                tooltip_col = _get_input_value(input, _prefixed_id(prefix, "spawn_conc_tooltip_col")) or ""
                if use_tooltip and tooltip_col:
                    center_entry = prot_data.get(center_id) or _resolve_protein_entry(center_id)
                    tooltip_val = ""
                    if isinstance(center_entry, dict):
                        tooltip_val = center_entry.get(tooltip_col) or center_entry.get("annotations", "")
                    extracted: List[str] = []
                    if tooltip_val:
                        for part in re.split(r"[,;\\s]+", str(tooltip_val)):
                            part = part.strip()
                            if part:
                                extracted.append(part)
                    if extracted:
                        ring_ids = extracted
                count = len(ring_ids)
                if count:
                    add_arrows = bool(_get_input_value(input, _prefixed_id(prefix, "spawn_conc_arrows")))
                    radius_mode = str(_get_input_value(input, _prefixed_id(prefix, "spawn_conc_radius_mode")) or "auto").lower()
                    if radius_mode == "fixed":
                        radius = _to_float(_get_input_value(input, _prefixed_id(prefix, "spawn_conc_radius_fixed")), 220) or 220.0
                    else:
                        space = _to_float(_get_input_value(input, _prefixed_id(prefix, "spawn_conc_space")), 70) or 70.0
                        radius = max(60.0, (count * space) / (2 * math.pi) * 2.0)
                    arrow_stop = _to_float(_get_input_value(input, _prefixed_id(prefix, "spawn_conc_arrow_stop")), 24.0) or 0.0
                    for i, uni in enumerate(ring_ids):
                        angle = (2 * math.pi * i) / max(1, count)
                        cx = center_x + radius * math.cos(angle)
                        cy = center_y + radius * math.sin(angle)
                        x = cx - pb_w * 0.5
                        y = cy - pb_h * 0.5
                        _append_protbox(uni, i + 1, x, y)
                        if add_arrows:
                            start_cx = x + pb_w * 0.5
                            start_cy = y + pb_h * 0.5
                            vx = center_x - start_cx
                            vy = center_y - start_cy
                            dist = math.hypot(vx, vy) or 1.0
                            stop_r = max(0.0, arrow_stop)
                            end_x = center_x - (vx / dist) * stop_r if stop_r else center_x
                            end_y = center_y - (vy / dist) * stop_r if stop_r else center_y
                            arrows = payload.setdefault("arrows", [])
                            arrows.append(
                                {
                                    "x1": start_cx,
                                    "y1": start_cy,
                                    "x2": end_x,
                                    "y2": end_y,
                                    "line": "arrow",
                                    "type": "",
                                }
                            )
            else:
                for idx, uni in enumerate(ids):
                    row = idx // grid_max_per_row
                    col = idx % grid_max_per_row
                    x = start_x + col * grid_step_x
                    y = start_y + row * grid_step_y
                    _append_protbox(uni, idx, x, y)
            if not added:
                state["status"].set("No new protboxes added (possibly duplicates).")
                return
            payload["_persist_token"] = time.time()
            state["json"].set(payload)
            state["status"].set(f"Added {added} protbox(es).")

        @output(id=_prefixed_id(prefix, "pathway_preview"))
        @render.ui
        def pathway_preview():
            if _active_bookmark() != prefix:
                return ui.div({"class": "alert alert-info"}, "Select this tab to view its pathway.")
            fc_idx = state["fc_index"].get()
            data = state["json"].get()
            if not data:
                if cfg.get("start_blank"):
                    prompt = "Click Reset Canvas to start building a pathway." if cfg.get("key") == "figure" else "Click Create Blank Canvas to start building a pathway."
                    return ui.div({"class": "alert alert-warning"}, prompt)
                return ui.div({"class": "alert alert-warning"}, "No pathway data to display.")
            settings = data.get("general_data", {}).get("settings", {})
            show_bg = bool(settings.get("show_background_image", False))
            data_with_fc = dict(data)
            data_with_fc["_active_fc_index"] = fc_idx
            return create_pathway_svg(data_with_fc, show_kegg_bg=show_bg)

        @output(id=_prefixed_id(prefix, "status_message"))
        @render.text
        def status_message() -> str:
            return state["status"].get()

        @output(id=_prefixed_id(prefix, "json_summary"))
        @render.text
        def json_summary() -> str:
            data = state["json"].get()
            if not data:
                return "Summary: not available"
            protbox_count = len(data.get("protbox_data", []))
            group_count = len(data.get("groups", []))
            arrow_count = len(data.get("arrows", []))
            return f"Protboxes: {protbox_count} | Groups: {group_count} | Arrows: {arrow_count}"

        @output(id=_prefixed_id(prefix, "download_json"))
        @render.download(filename=f"{prefix}_pathway.json")
        def download_json():
            data = state["json"].get()
            if not data:
                yield b""
                return
            yield json.dumps(data, indent=2).encode("utf-8")

        @output(id=_prefixed_id(prefix, "download_custom_pathway"))
        @render.download(filename=f"{prefix}_custom_pathway.json")
        def download_custom_pathway():
            export_file = os.path.join(JSON_PREVIEW_DIR, f"{prefix}_custom_pathway_export.json")
            try:
                if os.path.exists(export_file):
                    with open(export_file, "rb") as export_fh:
                        yield export_fh.read()
                        return
                data = state["json"].get()
                if not data:
                    yield json.dumps(
                        {"schema_version": CUSTOM_LAYOUT_SCHEMA_VERSION, "pathway_source": "custom"},
                        ensure_ascii=True,
                    ).encode("utf-8")
                    return
                export_payload = _build_custom_layout_export(data)
                yield json.dumps(export_payload, indent=2, ensure_ascii=True).encode("utf-8")
            except Exception as exc:
                print(f"download_custom_pathway error: {exc}")
                yield json.dumps({"error": str(exc), "schema_version": CUSTOM_LAYOUT_SCHEMA_VERSION}, ensure_ascii=True).encode("utf-8")

        @reactive.Effect
        @reactive.event(getattr(input, _prefixed_id(prefix, "export_custom_pathway")))
        def _prepare_custom_export():
            try:
                state["export_pending"].set(True)
                _send_custom_message(
                    session,
                    "request_export_snapshot",
                    {
                        "prefix": prefix,
                    },
                )
            except Exception as exc:
                state["status"].set(f"Failed to export custom pathway: {exc}")
                print(f"_prepare_custom_export error: {exc}")
                state["export_pending"].set(False)
                _send_custom_message(
                    session,
                    "export_failed",
                    {
                        "button_id": _prefixed_id(prefix, "export_custom_pathway"),
                        "spinner_id": _prefixed_id(prefix, "export_custom_pathway_spinner"),
                    },
                )

        @reactive.Effect
        @reactive.event(getattr(input, _prefixed_id(prefix, "generate")))
        def _build_on_click():
            state["status"].set("Generating pathway...")
            build_json()

        @reactive.Effect
        @reactive.event(getattr(input, _prefixed_id(prefix, "upload_custom_pathway")))
        def _import_custom_pathway():
            upload = _get_input_value(input, _prefixed_id(prefix, "upload_custom_pathway"))
            if not upload:
                return
            file_info = upload[0]
            datapath = getattr(file_info, "datapath", None)
            if datapath is None and isinstance(file_info, dict):
                datapath = file_info.get("datapath")
            if not datapath:
                state["status"].set("Could not locate the uploaded custom pathway file.")
                return
            try:
                with open(datapath, "r", encoding="utf-8") as fh:
                    raw_data = json.load(fh)
                layout = _sanitize_custom_layout(raw_data)
            except (OSError, json.JSONDecodeError, ValueError) as exc:
                state["status"].set(f"Failed to import custom pathway: {exc}")
                return
            state["custom_layout"].set(layout)
            # Apply immediately
            if prefix == "web":
                try:
                    session.send_input_message(_prefixed_id(prefix, "pathway_id"), {"value": "Custom"})
                except Exception:
                    pass
            state["status"].set("Custom pathway imported and applied.")
            build_json()

        if cfg.get("start_blank"):
            build_json()

    for cfg in BOOKMARK_CONFIGS:
        _register_bookmark(cfg)


app = App(app_ui, server)


def _run_uvicorn_app():
    if uvicorn is None:  # pragma: no cover - runtime guard
        raise RuntimeError("uvicorn is required to run the Shiny app server")
    uvicorn.run(app, host=HOST, port=PORT, log_level="info")


def _launch_desktop_window(url: str):
    if webview is None:  # pragma: no cover - runtime guard
        print(f"pywebview is not installed. Open the Shiny app manually at {url}")
        return
    # Prefer the Edge backend on Windows; fall back to auto-selection.
    try:
        webview.config.gui = "edgechromium"  # type: ignore[attr-defined]
    except Exception:
        pass
    webview.create_window("MapKinase Settings", url)
    webview.start()


if __name__ == "__main__":
    if LAUNCH_DESKTOP_GUI:
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
    else:
        _run_uvicorn_app()










