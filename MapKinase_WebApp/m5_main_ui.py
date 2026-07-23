import atexit
import ast
import base64
import copy
import io
import json
import os
import re
import threading
import tempfile
import time
import sys
import csv
import asyncio
import html
import hashlib
import math
import xml.etree.ElementTree as ET
from argparse import Namespace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from PIL import Image

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

from shiny import App, reactive, render, ui

from MapKinase_WebApp.m4_json import DEFAULT_DATA, DEFAULT_SETTINGS, get_default_json
from MapKinase_WebApp.a1_factory import get_pathway_api
from MapKinase_WebApp.m2_protein_catalog import build_global_protein_catalog, ensure_global_protein_catalog
from MapKinase_WebApp.m1_file_processor import validate_protein_file, validate_ptm_file, validate_uploaded_file
from MapKinase_WebApp.mk_security_limits import (
    MAX_RUNS_PER_MINUTE,
    MAX_TOTAL_UPLOAD_SIZE_MB,
    MAX_UPLOAD_SIZE_MB,
    MIN_SECONDS_BETWEEN_RUNS,
    evaluate_run_throttle,
    log_run_guard_event,
    validate_total_upload_size_from_sizes,
)
from MapKinase_WebApp.mk_session_workspace import (
    SESSION_TMP_TTL_HOURS,
    cleanup_old_session_workspaces,
    cleanup_session_workspace,
    get_session_workspace,
    safe_session_path,
    safely_delete_temp_file,
)
from MapKinase_WebApp.mk_runtime_env import env_truthy, is_production_mode, runtime_environment
from MapKinase_WebApp.d2_psp_regulatorysites import load_regulatory_sites, annotate_ptm_dataset
from MapKinase_WebApp.d2_psp_kinasesubstrates import load_kinase_substrate_map, annotate_ptm_dataset_with_kinases
from MapKinase_WebApp.d1_transfer_kegg_annotations import load_kegg_map, annotate_protein_with_kegg
from MapKinase_WebApp.m6_rank_pathways import (
    build_gene_to_uniprot_map,
    candidate_uniprots_for_node,
    build_protein_lookup,
    compute_single_protein_scores,
    load_kegg_index,
    normalize_uniprot,
    parse_weights,
    rank_all_pathways,
    resolve_node_scores,
)

from MapKinase_WebApp.m3_svg_viewer import create_pathway_svg, _build_blank_canvas
from MapKinase_WebApp.m7_cst_viewer import (
    create_cst_pathway_viewer,
    get_cst_pathway_catalog,
    load_cst_pathway_payload,
    save_cst_overlay_state,
)
from MapKinase_WebApp.m11_cst_pathway_index import iter_effective_cst_modules
from MapKinase_WebApp.pathway_banlist import (
    filter_kegg_options_for_species,
    filter_wikipathways_options,
)

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
INDEX_FILES_DIR = os.path.join(BASE_DIR, "index_files")
ANNOTATION_FILES_DIR = os.path.join(BASE_DIR, "annotation_files")
HUMAN_ORTHOLOG_UNIPROT_COLUMN = "Human_Ortholog_UniProt"
HUMAN_ORTHOLOG_GENE_COLUMN = "Human_Ortholog_Gene"
_HUMAN_ORTHOLOG_LOOKUP_CACHE: Optional[Dict[str, Dict[str, Dict[str, str]]]] = None


def _resolve_human_ortholog_map_file() -> Optional[Path]:
    candidates = [
        Path(ANNOTATION_FILES_DIR) / "hsa_uniprot_map_ortholog_map_long.tsv",
    ]
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        candidates.extend([
            Path(sys._MEIPASS) / "MapKinase_WebApp" / "annotation_files" / "hsa_uniprot_map_ortholog_map_long.tsv",
            Path(sys._MEIPASS) / "annotation_files" / "hsa_uniprot_map_ortholog_map_long.tsv",
        ])
    for path in candidates:
        if path.exists():
            return path
    return None


def _ortholog_priority(entry: Dict[str, str]) -> Tuple[int, float]:
    orthology_type = str(entry.get("orthology_type") or "").strip().lower()
    try:
        identity = float(entry.get("identity") or 0.0)
    except (TypeError, ValueError):
        identity = 0.0
    if orthology_type == "ortholog_one2one":
        return (3, identity)
    if orthology_type == "ortholog_one2many":
        return (2, identity)
    if orthology_type == "ortholog_many2many":
        return (1, identity)
    return (0, identity)


def _load_human_ortholog_lookup() -> Dict[str, Dict[str, Dict[str, str]]]:
    global _HUMAN_ORTHOLOG_LOOKUP_CACHE
    if _HUMAN_ORTHOLOG_LOOKUP_CACHE is not None:
        return _HUMAN_ORTHOLOG_LOOKUP_CACHE
    lookup: Dict[str, Dict[str, Dict[str, str]]] = {"mmu": {}, "rno": {}}
    ortholog_path = _resolve_human_ortholog_map_file()
    if ortholog_path is None:
        print("Warning: human ortholog mapping file not found.")
        _HUMAN_ORTHOLOG_LOOKUP_CACHE = lookup
        return lookup
    try:
        with ortholog_path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for raw_row in reader:
                row = {str(k or "").strip(): str(v or "").strip() for k, v in dict(raw_row or {}).items()}
                human_uniprot = normalize_uniprot(
                    row.get("Human UniProt (resolved)") or row.get("Human UniProt (input)") or ""
                )
                human_gene = str(row.get("Human gene") or "").strip()
                if not human_uniprot:
                    continue
                for species_code, prefix in (("mmu", "Mouse"), ("rno", "Rat")):
                    source_uniprot = normalize_uniprot(row.get(f"{prefix} UniProt") or "")
                    if not source_uniprot:
                        continue
                    candidate = {
                        "human_uniprot": human_uniprot,
                        "human_gene": human_gene,
                        "source_gene": str(row.get(f"{prefix} gene") or "").strip(),
                        "orthology_type": str(row.get(f"{prefix} orthology type") or "").strip(),
                        "identity": str(row.get(f"{prefix} % identity") or "").strip(),
                    }
                    current = lookup[species_code].get(source_uniprot)
                    if current is None or _ortholog_priority(candidate) > _ortholog_priority(current):
                        lookup[species_code][source_uniprot] = candidate
    except Exception as exc:
        print(f"Warning: failed to load human ortholog mappings: {exc}")
    _HUMAN_ORTHOLOG_LOOKUP_CACHE = lookup
    return lookup


def _annotate_dataset_with_human_orthologs(dataset: Dict[str, Any], species_code: str) -> Dict[str, Any]:
    code = str(species_code or "").strip().lower()
    if code not in {"mmu", "rno"}:
        return dataset
    headers = list(dataset.get("headers") or [])
    rows = list(dataset.get("rows") or [])
    if not headers:
        return dataset
    ortholog_lookup = _load_human_ortholog_lookup().get(code, {})
    out_headers = list(headers)
    if HUMAN_ORTHOLOG_UNIPROT_COLUMN in out_headers:
        human_uniprot_idx = out_headers.index(HUMAN_ORTHOLOG_UNIPROT_COLUMN)
    else:
        human_uniprot_idx = len(out_headers)
        out_headers.append(HUMAN_ORTHOLOG_UNIPROT_COLUMN)
    if HUMAN_ORTHOLOG_GENE_COLUMN in out_headers:
        human_gene_idx = out_headers.index(HUMAN_ORTHOLOG_GENE_COLUMN)
    else:
        human_gene_idx = len(out_headers)
        out_headers.append(HUMAN_ORTHOLOG_GENE_COLUMN)
    out_rows: List[List[str]] = []
    for row in rows:
        values = list(row or [])
        if len(values) < len(headers):
            values.extend([""] * (len(headers) - len(values)))
        if len(values) < len(out_headers):
            values.extend([""] * (len(out_headers) - len(values)))
        source_uniprot = normalize_uniprot(values[0] if values else "")
        ortholog_entry = ortholog_lookup.get(source_uniprot, {})
        values[human_uniprot_idx] = str(ortholog_entry.get("human_uniprot") or "")
        values[human_gene_idx] = str(ortholog_entry.get("human_gene") or "")
        out_rows.append(values)
    output = dict(dataset)
    output["headers"] = out_headers
    output["rows"] = out_rows
    return output


def _resolve_cst_ortholog_columns(columns: Sequence[Any]) -> Tuple[str, str]:
    normalized = {str(col or "").strip(): str(col or "").strip() for col in list(columns or []) if str(col or "").strip()}
    return (
        normalized.get(HUMAN_ORTHOLOG_UNIPROT_COLUMN, ""),
        normalized.get(HUMAN_ORTHOLOG_GENE_COLUMN, ""),
    )


def _resolve_icons_dir() -> str:
    candidates = [os.path.join(BASE_DIR, "icons")]
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        candidates.extend([
            os.path.join(sys._MEIPASS, "MapKinase_WebApp", "icons"),
            os.path.join(sys._MEIPASS, "icons"),
        ])
    for path in candidates:
        if os.path.isdir(path):
            return path
    print(f"Warning: icons directory not found (searched: {candidates})")
    return candidates[0]


ICONS_DIR = _resolve_icons_dir()


def _resolve_outline_font_file(filename: str) -> Optional[Path]:
    candidates: List[Path] = [
        Path(BASE_DIR) / "fonts" / filename,
        Path(BASE_DIR) / "icons" / "fonts" / filename,
    ]
    windir = Path(os.environ.get("WINDIR", r"C:\Windows"))
    candidates.append(windir / "Fonts" / filename)
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        candidates.extend([
            Path(sys._MEIPASS) / "MapKinase_WebApp" / "fonts" / filename,
            Path(sys._MEIPASS) / "fonts" / filename,
            Path(sys._MEIPASS) / "MapKinase_WebApp" / "icons" / "fonts" / filename,
            Path(sys._MEIPASS) / "icons" / "fonts" / filename,
        ])
    for path in candidates:
        try:
            if path.exists() and path.is_file():
                return path
        except Exception:
            continue
    return None


def _load_outline_font_base64(filename: str) -> str:
    font_path = _resolve_outline_font_file(filename)
    if font_path is None:
        print(f"Warning: outline font '{filename}' not found.")
        return ""
    try:
        raw = font_path.read_bytes()
    except Exception as exc:
        print(f"Warning: failed to read outline font '{font_path}': {exc}")
        return ""
    return base64.b64encode(raw).decode("ascii")


OUTLINE_FONT_REGULAR_B64 = _load_outline_font_base64("segoeui.ttf")
OUTLINE_FONT_BOLD_B64 = _load_outline_font_base64("segoeuib.ttf")


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


def _resolve_species_ref_file() -> str:
    candidates = [os.path.join(BASE_DIR, "species_ref_list.csv")]
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        candidates.extend([
            os.path.join(sys._MEIPASS, "MapKinase_WebApp", "species_ref_list.csv"),
            os.path.join(sys._MEIPASS, "species_ref_list.csv"),
        ])
    for path in candidates:
        if os.path.exists(path):
            return path
    print(f"Warning: species_ref_list.csv not found (searched: {candidates})")
    return candidates[0]


def _icon_markup(icon_name: str, class_name: str = "mk-inline-icon") -> Any:
    icon_path = Path(ICONS_DIR) / f"{icon_name}.svg"
    fallback = ui.tags.span(icon_name[:1].upper(), {"aria-hidden": "true"})
    try:
        svg_text = icon_path.read_text(encoding="utf-8")
    except OSError:
        return fallback
    svg_text = re.sub(r"<\?xml[^>]*\?>", "", svg_text, count=1).strip()
    svg_text = svg_text.replace("#000000", "currentColor").replace("#000", "currentColor")
    escaped_class = html.escape(class_name, quote=True)
    if "<svg" in svg_text:
        svg_text = re.sub(
            r"<svg\b",
            f'<svg class="{escaped_class}" aria-hidden="true" focusable="false"',
            svg_text,
            count=1,
        )
    return ui.HTML(svg_text)


def _asset_data_uri(filename: str) -> str:
    asset_path = Path(ICONS_DIR) / filename
    suffix = asset_path.suffix.lower()
    mime_type = {
        ".gif": "image/gif",
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".webp": "image/webp",
        ".svg": "image/svg+xml",
    }.get(suffix, "application/octet-stream")
    try:
        raw = asset_path.read_bytes()
    except OSError:
        return ""
    return f"data:{mime_type};base64,{base64.b64encode(raw).decode('ascii')}"


def _asset_data_uri_from_path(path: Path) -> str:
    suffix = path.suffix.lower()
    mime_type = {
        ".gif": "image/gif",
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".webp": "image/webp",
        ".svg": "image/svg+xml",
    }.get(suffix, "application/octet-stream")
    try:
        raw = path.read_bytes()
    except OSError:
        return ""
    return f"data:{mime_type};base64,{base64.b64encode(raw).decode('ascii')}"


KEGG_PATHWAYS_FILE = _resolve_kegg_pathways_file()
KEGG_PATHWAY_MAX_MATCHES = 12
WIKIPATHWAYS_MAX_MATCHES = 25
WEB_PATHWAY_SOURCES = ["WikiPathways", "KEGG", "CST"]
DEFAULT_BG_OPACITY = 0.9
DEFAULT_BG_SCALE = 1.0
DEFAULT_BG_OFFSET_X = 0.0
DEFAULT_BG_OFFSET_Y = 0.0
DEFAULT_WIKIPATHWAYS_BG_SCALE = 1.0
DEFAULT_WIKIPATHWAYS_BG_OFFSET_X = 0.0
DEFAULT_WIKIPATHWAYS_BG_OFFSET_Y = 0.0
DEFAULT_BOX_Y_STRETCH = 1.0
SESSION_PREVIEW_FILENAME = "latest_preview.json"
RESOURCE_ROOT = getattr(sys, "_MEIPASS", PARENT_DIR)


def _resolve_sample_data_dir() -> str:
    candidates = [
        os.path.join(RESOURCE_ROOT, "Sample_input_files"),
        os.path.join(RESOURCE_ROOT, "sample_input_files"),
    ]
    for path in candidates:
        if os.path.isdir(path):
            return path
    return candidates[0]


SAMPLE_DATA_DIR = _resolve_sample_data_dir()
SAMPLE_PROTEIN_FILE = os.path.join(SAMPLE_DATA_DIR, "TCA_protein_for_mapk.txt")
SAMPLE_PTM_FILE = os.path.join(SAMPLE_DATA_DIR, "TCA_phospho_for_mapk_singlesite.csv")


SAMPLE_METABOLITE_FILE = os.path.join(SAMPLE_DATA_DIR, "example_metabolite_file.txt")
SPECIES_REF_PATH = _resolve_species_ref_file()
UPLOAD_ACCEPT_TYPES = [
    ".csv",
    ".tsv",
    ".txt",
    "text/plain",
    "text/csv",
]
DATA_USE_WARNING_TEXT = (
    "Map-Kinase is intended for non-sensitive research data only. Do not upload regulated, confidential, "
    "identifiable, clinical, FERPA, HIPAA, export-controlled, or otherwise sensitive data."
)
# Flip to True to mirror terminal stdout/stderr into TERMINAL_LOG_FILE by default.
TERMINAL_LOG_DEFAULT = False
TERMINAL_LOG_DIR_DEFAULT = os.path.join(tempfile.gettempdir(), "mapkinase_logs")
TERMINAL_LOG_FILE = os.environ.get(
    "M5_TERMINAL_LOG_FILE",
    os.path.join(TERMINAL_LOG_DIR_DEFAULT, "m5_terminal_output.txt"),
)
TERMINAL_LOG_MAX_BYTES = max(1, int(os.environ.get("M5_TERMINAL_LOG_MAX_BYTES", str(5 * 1024 * 1024))))
TERMINAL_LOG_BACKUPS = max(1, int(os.environ.get("M5_TERMINAL_LOG_BACKUPS", "5")))
MANUAL_BUILD_ONLY = True
GUI_POPUP = False

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


def _resolve_demo_species_key() -> str:
    for key, cfg in SPECIES_CHOICES.items():
        if str(cfg.get("label") or "").strip().lower() == "mouse":
            return key
    return "mouse" if "mouse" in SPECIES_CHOICES else DEFAULT_SPECIES


DEMO_SPECIES = _resolve_demo_species_key()

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
        "label": "Kinase-Substrate",
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


RUNTIME_ENVIRONMENT = runtime_environment()
IS_PRODUCTION_MODE = is_production_mode()
DEBUG_UI_ENABLED = env_truthy("MAPKINASE_ENABLE_DEBUG_UI", default=not IS_PRODUCTION_MODE)
DEBUG_EXPORT_ENABLED = env_truthy(
    "MAPKINASE_ENABLE_DEBUG_EXPORT",
    default=(not IS_PRODUCTION_MODE and DEBUG_UI_ENABLED),
)

BUILD_GLOBAL_CATALOG_ON_STARTUP = _env_var_truthy(
    "M5_BUILD_GLOBAL_CATALOG_ON_STARTUP",
    default=False,
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


RAW_KEGG_PATHWAY_OPTIONS = _load_kegg_pathways(KEGG_PATHWAYS_FILE)
WIKIPATHWAYS_CACHE: Dict[str, List[Dict[str, str]]] = {}
WIKIPATHWAYS_ORG_CACHE: Optional[set] = None


def _get_kegg_pathway_options_for_species(species_code: str) -> List[Dict[str, str]]:
    return filter_kegg_options_for_species(RAW_KEGG_PATHWAY_OPTIONS, species_code)


def _resolve_species_code_for_catalog(organism: str, fallback: Optional[str] = None) -> str:
    names = [str(v).strip().lower() for v in (organism, fallback) if v]
    for key, info in SPECIES_CHOICES.items():
        candidates = {
            key.strip().lower(),
            str(info.get("label") or "").strip().lower(),
            str(info.get("species") or "").strip().lower(),
            str(info.get("code") or "").strip().lower(),
        }
        if any(name in candidates for name in names):
            return str(info.get("code") or "").strip().lower()
    return ""


def _normalize_catalog_species_name(value: Any) -> str:
    return str(value or "").strip().lower().replace("_", " ")


def _load_wikipathways_catalog_from_local_index(organism: str, fallback: Optional[str] = None) -> List[Dict[str, str]]:
    species_code = _resolve_species_code_for_catalog(organism, fallback)
    if not species_code:
        return []
    index_path = Path(INDEX_FILES_DIR) / f"wikipathways_index_{species_code}.json"
    if not index_path.exists():
        return []
    try:
        payload = json.loads(index_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"Warning: failed to read local WikiPathways index '{index_path.name}': {exc}")
        return []
    pathways = payload.get("pathways") if isinstance(payload, dict) else None
    if not isinstance(pathways, list):
        return []
    meta = payload.get("meta") if isinstance(payload, dict) else {}
    species_name = ""
    if isinstance(meta, dict):
        species_name = str(meta.get("species") or organism or fallback or "").strip()
    options: List[Dict[str, str]] = []
    for row in pathways:
        if not isinstance(row, dict):
            continue
        path_id = str(row.get("pathway_id") or row.get("id") or row.get("wpid") or "").strip().upper()
        name = str(row.get("name") or "").strip()
        if not path_id or not name:
            continue
        org = str(row.get("species") or species_name).strip()
        label = f"{path_id} | {name}"
        if org:
            label = f"{label} ({org})"
        options.append({"id": path_id, "name": name, "species": org, "label": label})
    options = filter_wikipathways_options(options, species_code)
    options.sort(key=lambda opt: opt.get("name", "").lower())
    return options


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
    cache_key = "||".join(
        str(v).strip().lower()
        for v in (organism or "", fallback or "")
    ) or "__all__"
    if cache_key in WIKIPATHWAYS_CACHE:
        return WIKIPATHWAYS_CACHE[cache_key]
    local_options = _load_wikipathways_catalog_from_local_index(organism, fallback=fallback)
    if local_options:
        WIKIPATHWAYS_CACHE[cache_key] = local_options
        return local_options
    options: List[Dict[str, str]] = []

    def _extract_rows(payload: Any) -> List[Dict[str, Any]]:
        if isinstance(payload, list):
            return [row for row in payload if isinstance(row, dict)]
        if not isinstance(payload, dict):
            return []
        candidates = [
            payload.get("pathways"),
            payload.get("result"),
            payload.get("listPathwaysResult"),
            payload.get("pathway"),
        ]
        for candidate in candidates:
            if isinstance(candidate, list):
                return [row for row in candidate if isinstance(row, dict)]
            if isinstance(candidate, dict):
                nested = candidate.get("pathways") or candidate.get("pathway")
                if isinstance(nested, list):
                    return [row for row in nested if isinstance(row, dict)]
                return [candidate]
        return [payload]

    def _fetch_rows(name: str, quiet: bool = False) -> List[Dict[str, Any]]:
        try:
            import requests
            response = requests.get("https://www.wikipathways.org/json/listPathways.json", timeout=8)
            response.raise_for_status()
            payload = response.json()
            organism_norm = _normalize_catalog_species_name(name)
            rows = []
            organisms = payload.get("organisms") if isinstance(payload, dict) else None
            if isinstance(organisms, list):
                for organism_row in organisms:
                    if not isinstance(organism_row, dict):
                        continue
                    candidates = {
                        _normalize_catalog_species_name(organism_row.get("latin")),
                        _normalize_catalog_species_name(organism_row.get("common")),
                    }
                    if organism_norm and organism_norm not in candidates:
                        continue
                    pathways = organism_row.get("pathways")
                    if isinstance(pathways, list):
                        rows = [row for row in pathways if isinstance(row, dict)]
                    break
            if not rows:
                rows = _extract_rows(payload)
            return rows
        except Exception as exc:  # pragma: no cover - network/service issues
            if not quiet:
                print(f"Warning: list_pathways failed for '{name}': {exc}")
            return []

    try:
        names_to_try: List[str] = []
        if organism:
            names_to_try.append(organism)
        if fallback and (not organism or fallback.lower() != organism.lower()):
            names_to_try.append(fallback)

        records: List[Dict[str, Any]] = []
        for idx, nm in enumerate(names_to_try):
            quiet = idx < len(names_to_try) - 1
            records = _fetch_rows(nm, quiet=quiet)
            if records:
                break
        if not records:
            WIKIPATHWAYS_CACHE[cache_key] = []
            return []

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
        species_code = _resolve_species_code_for_catalog(organism, fallback)
        if species_code:
            options = filter_wikipathways_options(options, species_code)
        options.sort(key=lambda opt: opt.get("name", "").lower())
    except Exception as exc:
        print(f"Warning: failed to load WikiPathways catalogue for '{organism or fallback or 'all'}': {exc}")
    WIKIPATHWAYS_CACHE[cache_key] = options
    return options

TERMINAL_LOG_ENABLED = _env_var_truthy("M5_TERMINAL_LOG", TERMINAL_LOG_DEFAULT)
DEBUG_FILE_OUTPUT_ENABLED = _env_var_truthy("MAPKINASE_DEBUG_FILE_OUTPUT", False)
PERSISTENT_CST_SAVE_ENABLED = _env_var_truthy("MAPKINASE_ENABLE_PERSISTENT_CST_SAVE", False)
CUSTOM_LAYOUT_SCHEMA_VERSION = 1

try:
    _removed_workspace_count = cleanup_old_session_workspaces(SESSION_TMP_TTL_HOURS)
    if _removed_workspace_count:
        print(
            f"Removed {_removed_workspace_count} stale session workspace(s) "
            f"(TTL={SESSION_TMP_TTL_HOURS}h)."
        )
except Exception as _workspace_cleanup_exc:
    print(f"Warning: failed to run session workspace TTL cleanup: {_workspace_cleanup_exc}")


def _coerce_float(value: Any) -> Optional[float]:
    if value in (None, "", False):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _background_preview_from_settings(
    settings: Optional[Dict[str, Any]] = None,
    show: bool = False,
) -> Dict[str, float | bool]:
    cfg = dict(settings or {})
    pathway_source = str(cfg.get("pathway_source") or "").strip().lower()
    if pathway_source == "wikipathways":
        scale = float(cfg.get("bg_scale", DEFAULT_WIKIPATHWAYS_BG_SCALE))
        offset_x = float(cfg.get("bg_offset_x", DEFAULT_WIKIPATHWAYS_BG_OFFSET_X))
        offset_y = float(cfg.get("bg_offset_y", DEFAULT_WIKIPATHWAYS_BG_OFFSET_Y))
    else:
        scale = DEFAULT_BG_SCALE
        offset_x = DEFAULT_BG_OFFSET_X
        offset_y = DEFAULT_BG_OFFSET_Y
    return {
        "show": bool(show),
        "opacity": float(cfg.get("bg_opacity", DEFAULT_BG_OPACITY)),
        "offset_x": offset_x,
        "offset_y": offset_y,
        "scale": scale,
    }


_SPREADSHEET_FORMULA_PREFIXES = ("=", "+", "-", "@", "\t", "\r")
_HTML_BR_RE = re.compile(r"(?is)<br\s*/?>")
_HTML_TAG_RE = re.compile(r"(?is)<[^>]*>")
_UNIPROT_ACCESSION_RE = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])(?:-\d+)?\b", re.IGNORECASE)


def _sanitize_layout_text_plain(value: Any, max_chars: int = 4000) -> str:
    raw = str(value or "")
    raw = raw.replace("\r\n", "\n").replace("\r", "\n").replace("\x00", "")
    if len(raw) > max_chars:
        raw = raw[:max_chars]
    return raw


def _sanitize_layout_text_html(value: Any, max_chars: int = 4000) -> str:
    # Treat imported/exported text-box content as inert text, not executable markup.
    raw = _sanitize_layout_text_plain(value, max_chars=max_chars)
    raw = _HTML_BR_RE.sub("\n", raw)
    raw = _HTML_TAG_RE.sub("", raw)
    raw = html.unescape(raw)
    raw = _sanitize_layout_text_plain(raw, max_chars=max_chars)
    return "<br/>".join(html.escape(line) for line in raw.split("\n"))


def _sanitize_layout_scalar(value: Any, *, max_chars: int = 256) -> Any:
    if isinstance(value, (int, float, bool)) or value is None:
        return value
    return _sanitize_layout_text_plain(value, max_chars=max_chars)


def _sanitize_layout_metadata(value: Any, *, max_chars: int = 800, depth: int = 0) -> Any:
    if isinstance(value, (int, float, bool)) or value is None:
        return value
    if isinstance(value, str):
        return _sanitize_layout_text_plain(value, max_chars=max_chars)
    if depth >= 3:
        return _sanitize_layout_text_plain(value, max_chars=max_chars)
    if isinstance(value, list):
        return [_sanitize_layout_metadata(item, max_chars=max_chars, depth=depth + 1) for item in value[:200]]
    if isinstance(value, dict):
        cleaned: Dict[str, Any] = {}
        for idx, (key, val) in enumerate(value.items()):
            if idx >= 200:
                break
            clean_key = _sanitize_layout_text_plain(key, max_chars=120)
            if clean_key:
                cleaned[clean_key] = _sanitize_layout_metadata(val, max_chars=max_chars, depth=depth + 1)
        return cleaned
    return _sanitize_layout_text_plain(value, max_chars=max_chars)


def _custom_layout_protbox_uniprots(protbox: Dict[str, Any]) -> List[str]:
    if not isinstance(protbox, dict):
        return []
    raw_values = [
        protbox.get("selected_uniprot"),
        protbox.get("proteins"),
        protbox.get("uniprot_ids"),
        protbox.get("uniprot"),
        protbox.get("protein_ids"),
    ]
    found: List[str] = []
    seen = set()

    def _add_candidate(value: Any) -> None:
        normalized = normalize_uniprot(value)
        if normalized and normalized not in seen:
            found.append(normalized)
            seen.add(normalized)

    def _ingest(value: Any, depth: int = 0) -> None:
        if value in (None, "") or depth > 4:
            return
        if isinstance(value, (list, tuple, set)):
            for item in value:
                _ingest(item, depth + 1)
            return
        if isinstance(value, dict):
            for item in value.values():
                _ingest(item, depth + 1)
            return
        text = str(value).strip()
        if not text:
            return
        try:
            parsed = ast.literal_eval(text)
        except (ValueError, SyntaxError):
            parsed = None
        if parsed is not None and parsed != text and isinstance(parsed, (list, tuple, set, dict, str)):
            _ingest(parsed, depth + 1)
            return
        regex_matches = _UNIPROT_ACCESSION_RE.findall(text)
        if regex_matches:
            for match in regex_matches:
                _add_candidate(match)
            return
        for part in re.split(r"[;,\s]+", text):
            cleaned = part.strip().strip("()[]{}'\"")
            if cleaned:
                _add_candidate(cleaned)

    for raw_value in raw_values:
        _ingest(raw_value)
    return found


def _neutralize_spreadsheet_formula_cell(value: Any) -> Any:
    if not isinstance(value, str):
        return value
    if not value:
        return value
    check = value.lstrip("\ufeff")
    if not check:
        return value
    if check.startswith("'"):
        return value
    if check[:1] in _SPREADSHEET_FORMULA_PREFIXES:
        return f"'{value}"
    return value


def _build_custom_layout_export(payload: Dict[str, Any]) -> Dict[str, Any]:
    settings = payload.get("general_data", {}).get("settings", {})
    export_data: Dict[str, Any] = {
        "schema_version": CUSTOM_LAYOUT_SCHEMA_VERSION,
        "pathway_id": _sanitize_layout_text_plain(settings.get("pathway_id"), max_chars=128),
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
                    if id_key == "text_id" and extra == "html":
                        entry[extra] = _sanitize_layout_text_html(item[extra], max_chars=12000)
                    elif extra == "label":
                        entry[extra] = _sanitize_layout_text_plain(item[extra], max_chars=800)
                    elif extra in {"text_style", "bgcolor", "fgcolor", "border_color"}:
                        entry[extra] = _sanitize_layout_scalar(item[extra], max_chars=400)
                    else:
                        entry[extra] = _sanitize_layout_scalar(item[extra], max_chars=800)
            # Preserve simple metadata for export
            for key, val in item.items():
                if key in entry or (extra_keys and key in extra_keys):
                    continue
                if isinstance(val, (str, int, float, bool, list, dict)):
                    if isinstance(val, str):
                        entry[key] = _sanitize_layout_text_plain(val, max_chars=1200)
                    else:
                        entry[key] = _sanitize_layout_metadata(val, max_chars=1200)
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
                arrow_entry[key] = _sanitize_layout_scalar(val, max_chars=400) if isinstance(val, str) else val
        export_data["arrows"].append(arrow_entry)
    for group in payload.get("groups", []) or []:
        if not isinstance(group, dict):
            continue
        entry: Dict[str, Any] = {}
        for key, val in group.items():
            if val in (None, ""):
                continue
            if isinstance(val, (str, int, float, bool, list, dict)):
                entry[key] = _sanitize_layout_scalar(val, max_chars=1200) if isinstance(val, str) else val
        if entry:
            export_data["groups"].append(entry)
    return export_data


def _sanitize_custom_layout(raw_data: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(raw_data, dict):
        raise ValueError("Custom pathway file must be a JSON object.")
    try:
        schema_version = int(raw_data.get("schema_version") or CUSTOM_LAYOUT_SCHEMA_VERSION)
    except (TypeError, ValueError):
        schema_version = CUSTOM_LAYOUT_SCHEMA_VERSION
    sanitized: Dict[str, Any] = {
        "schema_version": schema_version,
        "pathway_id": _sanitize_layout_text_plain(raw_data.get("pathway_id"), max_chars=128),
        "pathway_source": _sanitize_layout_text_plain(raw_data.get("pathway_source"), max_chars=64),
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
                    if section_name == "text_data" and extra == "html":
                        entry[extra] = _sanitize_layout_text_html(item[extra], max_chars=12000)
                    elif extra == "label":
                        entry[extra] = _sanitize_layout_text_plain(item[extra], max_chars=800)
                    elif extra in {"text_style", "bgcolor", "fgcolor", "border_color"}:
                        entry[extra] = _sanitize_layout_scalar(item[extra], max_chars=400)
                    else:
                        entry[extra] = _sanitize_layout_metadata(item[extra], max_chars=800)
            for key, val in item.items():
                if key in entry or key in {id_key, "id", "entry_id", "x", "y", "width", "height"}:
                    continue
                if val in (None, ""):
                    continue
                if not isinstance(val, (str, int, float, bool, list, dict)):
                    continue
                if section_name == "text_data" and key == "html":
                    entry[key] = _sanitize_layout_text_html(val, max_chars=12000)
                elif key == "label":
                    entry[key] = _sanitize_layout_text_plain(val, max_chars=800)
                elif key in {"text_style", "bgcolor", "fgcolor", "border_color"}:
                    entry[key] = _sanitize_layout_scalar(val, max_chars=400)
                else:
                    entry[key] = _sanitize_layout_metadata(val, max_chars=1200)
            if section_name == "text_data":
                if "html" not in entry:
                    entry["html"] = _sanitize_layout_text_html(item.get("label", ""), max_chars=12000)
                if "label" not in entry:
                    entry["label"] = _sanitize_layout_text_plain(item.get("label", ""), max_chars=800)
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
            arrow_entry[key] = _sanitize_layout_text_plain(value, max_chars=120)
        sanitized["arrows"].append(arrow_entry)
    for group in raw_data.get("groups", []) or []:
        if not isinstance(group, dict):
            continue
        entry: Dict[str, Any] = {}
        for key, val in group.items():
            if val in (None, ""):
                continue
            if isinstance(val, (str, int, float, bool, list, dict)):
                entry[key] = _sanitize_layout_text_plain(val, max_chars=1200) if isinstance(val, str) else val
        if entry:
            sanitized["groups"].append(entry)
    return sanitized


def _apply_custom_layout(payload: Dict[str, Any], layout: Optional[Dict[str, Any]], *, include_missing: bool = False) -> None:
    if not payload or not layout:
        return

    def _apply_section(section_name: str, id_key: str, extra_keys: Optional[Sequence[str]] = None) -> None:
        layout_entries = [entry for entry in layout.get(section_name, []) if isinstance(entry, dict) and entry.get(id_key)]
        overrides = {str(entry[id_key]).strip(): entry for entry in layout_entries if str(entry[id_key]).strip()}
        existing_ids = set()
        for item in payload.get(section_name, []) or []:
            entry_id = str(item.get(id_key) or "").strip()
            if not entry_id:
                continue
            existing_ids.add(entry_id)
            if entry_id not in overrides:
                continue
            override = overrides[entry_id]
            for key in ("x", "y", "width", "height"):
                if key in override:
                    item[key] = override[key]
            for extra in extra_keys or []:
                if extra in override:
                    if section_name == "text_data" and extra == "html":
                        item[extra] = _sanitize_layout_text_html(override[extra], max_chars=12000)
                    elif extra == "label":
                        item[extra] = _sanitize_layout_text_plain(override[extra], max_chars=800)
                    elif extra in {"text_style", "bgcolor", "fgcolor", "border_color"}:
                        item[extra] = _sanitize_layout_scalar(override[extra], max_chars=400)
                    else:
                        item[extra] = _sanitize_layout_metadata(override[extra], max_chars=800)
        if include_missing:
            target = payload.setdefault(section_name, [])
            for entry in layout_entries:
                entry_id = str(entry.get(id_key) or "").strip()
                if entry_id and entry_id not in existing_ids:
                    target.append(dict(entry))
                    existing_ids.add(entry_id)

    _apply_section("protbox_data", "protbox_id", extra_keys=("proteins", "uniprot_ids", "uniprot", "protein_ids", "label"))
    _apply_section("compound_data", "compound_id")
    _apply_section("text_data", "text_id", extra_keys=("html", "text_style", "bgcolor", "fgcolor", "border_color", "label"))

    if layout.get("arrows") is not None:
        payload["arrows"] = [dict(arrow) for arrow in layout["arrows"]]
    if include_missing and layout.get("groups") is not None:
        payload["groups"] = [dict(group) for group in layout["groups"] if isinstance(group, dict)]


def _enable_terminal_logging(path: str) -> None:
    if not path:
        return
    target = Path(str(path)).expanduser()
    forbidden_tokens = {"www", "static", "public", "jsonfiles", "cache", "stored_pathways", "output"}
    if {part.lower() for part in target.parts} & forbidden_tokens:
        print("Warning: refusing to enable terminal file logging inside public/static/cache paths.")
        return
    log_dir = os.path.dirname(str(target))
    if log_dir:
        try:
            os.makedirs(log_dir, exist_ok=True)
        except OSError:
            # Directory might be read-only; skip enabling logging in that case.
            return
    if target.exists():
        try:
            if target.stat().st_size >= TERMINAL_LOG_MAX_BYTES:
                for idx in range(TERMINAL_LOG_BACKUPS, 0, -1):
                    src = target.with_name(f"{target.name}.{idx}")
                    dst = target.with_name(f"{target.name}.{idx + 1}")
                    if idx == TERMINAL_LOG_BACKUPS and src.exists():
                        src.unlink(missing_ok=True)
                    elif src.exists():
                        src.replace(dst)
                target.replace(target.with_name(f"{target.name}.1"))
        except OSError:
            pass
    try:
        log_handle = open(target, "a", encoding="utf-8")
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


def _attach_pathway_background_image(data: Any, force: bool = False) -> Tuple[Any, bool]:
    if not isinstance(data, dict):
        return data, False

    if not force and data.get("kegg_bg_image"):
        return data, False

    settings = data.get("general_data", {}).get("settings", {})
    pathway_source = str(settings.get("pathway_source", "")).lower()
    pathway_id = settings.get("pathway_id") or data.get("pathway_id")
    if pathway_source not in {"kegg", "wikipathways"} or not pathway_id:
        return data, False

    try:
        api = get_pathway_api(pathway_source)
        img = api.download_pathway_image(pathway_id)
        if img is None:
            return data, False
        output_width = img.width
        output_height = img.height
        if pathway_source == "wikipathways":
            gpml_path = api.download_pathway_data(pathway_id, species_hint=settings.get("species"))
            if gpml_path and os.path.exists(gpml_path):
                try:
                    tree = ET.parse(gpml_path)
                    root = tree.getroot()
                    board_graphics = root.find("{http://pathvisio.org/GPML/2013a}Graphics")
                    if board_graphics is not None:
                        board_width = float(board_graphics.get("BoardWidth", "0") or 0)
                        board_height = float(board_graphics.get("BoardHeight", "0") or 0)
                        if board_width > 0 and board_height > 0 and img.width > 0 and img.height > 0:
                            render_scale = min(img.width / board_width, img.height / board_height)
                            render_width = board_width * render_scale
                            render_height = board_height * render_scale
                            margin_x = max(0.0, (img.width - render_width) / 2.0)
                            margin_y = max(0.0, (img.height - render_height) / 2.0)
                            left = max(0, int(round(margin_x)))
                            top = max(0, int(round(margin_y)))
                            right = min(img.width, int(round(img.width - margin_x)))
                            bottom = min(img.height, int(round(img.height - margin_y)))
                            if right > left and bottom > top:
                                img = img.crop((left, top, right, bottom))
                            output_width = max(1, int(round(board_width)))
                            output_height = max(1, int(round(board_height)))
                            if img.width != output_width or img.height != output_height:
                                try:
                                    resampling = Image.Resampling.LANCZOS
                                except AttributeError:
                                    resampling = Image.LANCZOS
                                img = img.resize((output_width, output_height), resampling)
                except Exception:
                    output_width = img.width
                    output_height = img.height
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("ascii")
        updated = dict(data)
        updated["kegg_bg_image"] = f"data:image/png;base64,{b64}"
        updated["kegg_bg_size"] = {"width": output_width, "height": output_height}
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
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _normalize_gradient_stop_entries(
    gradient_stops: Any,
) -> List[Dict[str, Any]]:
    if not isinstance(gradient_stops, (list, tuple)):
        return []
    parsed: List[Tuple[float, List[int]]] = []
    fallback_color = [128, 128, 128]
    for entry in gradient_stops:
        if not isinstance(entry, dict):
            continue
        raw_value = entry.get("value")
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(value):
            continue
        raw_color = entry.get("color", entry.get("rgb"))
        color: List[int]
        if isinstance(raw_color, str):
            color = list(_hex_to_rgb(raw_color, tuple(fallback_color)))
        elif isinstance(raw_color, (list, tuple)) and len(raw_color) >= 3:
            color = [int(max(0, min(255, round(float(c))))) for c in list(raw_color)[:3]]
        else:
            continue
        parsed.append((value, color))
    if not parsed:
        return []
    dedup: Dict[float, List[int]] = {}
    for value, color in parsed:
        dedup[float(value)] = list(color)
    ordered = sorted(dedup.items(), key=lambda kv: kv[0])
    return [{"value": float(v), "color": list(c)} for v, c in ordered]


def _default_gradient_stops() -> List[Dict[str, Any]]:
    from_defaults = _normalize_gradient_stop_entries(DEFAULT_SETTINGS.get("gradient_stops"))
    if len(from_defaults) >= 2:
        return from_defaults
    return _normalize_gradient_stop_entries(
        [
            {"value": float(DEFAULT_SETTINGS.get("max_negative", -2)), "color": list(DEFAULT_SETTINGS.get("negative_color", (179, 21, 41)))},
            {"value": 0.0, "color": [255, 255, 255]},
            {"value": float(DEFAULT_SETTINGS.get("max_positive", 2)), "color": list(DEFAULT_SETTINGS.get("positive_color", (16, 101, 171)))},
        ]
    )


def _normalize_gradient_stops_with_fallback(gradient_stops: Any) -> List[Dict[str, Any]]:
    normalized = _normalize_gradient_stop_entries(gradient_stops)
    if len(normalized) >= 2:
        return normalized
    return _default_gradient_stops()


def _gradient_color_from_stops(value: float, gradient_stops: List[Dict[str, Any]]) -> List[int]:
    if not gradient_stops:
        return [128, 128, 128]
    if value <= float(gradient_stops[0]["value"]):
        return list(gradient_stops[0]["color"])
    if value >= float(gradient_stops[-1]["value"]):
        return list(gradient_stops[-1]["color"])
    for idx in range(len(gradient_stops) - 1):
        left = gradient_stops[idx]
        right = gradient_stops[idx + 1]
        left_v = float(left["value"])
        right_v = float(right["value"])
        if value < left_v or value > right_v:
            continue
        if right_v == left_v:
            return list(right["color"])
        t = (value - left_v) / (right_v - left_v)
        left_c = list(left["color"])
        right_c = list(right["color"])
        return [
            int(max(0, min(255, round((1 - t) * left_c[c_idx] + t * right_c[c_idx]))))
            for c_idx in range(3)
        ]
    return list(gradient_stops[-1]["color"])


def _compute_gradient_color(
    value: Optional[float],
    neg: Sequence[int],
    pos: Sequence[int],
    max_neg: float,
    max_pos: float,
    gradient_stops: Optional[List[Dict[str, Any]]] = None,
) -> Optional[List[int]]:
    if value is None or not math.isfinite(value):
        return None
    stops = _normalize_gradient_stop_entries(gradient_stops or [])
    if len(stops) >= 2:
        return _gradient_color_from_stops(float(value), stops)
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
    gradient_stops = _normalize_gradient_stop_entries(color_override.get("gradient_stops"))
    for protein in protein_data.values():
        _recolor_entry(protein, neg, pos, max_neg, max_pos, gradient_stops=gradient_stops)
        ptms = protein.get("PTMs") or {}
        for ptm in ptms.values():
            _recolor_entry(ptm, neg, pos, max_neg, max_pos, gradient_stops=gradient_stops)


def _recolor_entry(
    entry: Dict[str, Any],
    neg: Sequence[int],
    pos: Sequence[int],
    max_neg: float,
    max_pos: float,
    gradient_stops: Optional[List[Dict[str, Any]]] = None,
) -> None:
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
        new_color = _compute_gradient_color(fc_value, neg, pos, max_neg, max_pos, gradient_stops=gradient_stops)
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
    for opt in _get_kegg_pathway_options_for_species(species_code):
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
    protein_match_mode = str(
        _get_input_value(input, "settings_protein_match_mode") or "uniprot_id"
    ).strip().lower()
    overrides["protein_match_mode"] = (
        protein_match_mode
        if (
            prefix == "web"
            and pathway_source in {"kegg", "wikipathways"}
            and protein_match_mode in {"uniprot_id", "gene_symbol"}
        )
        else "uniprot_id"
    )
    overrides["protein_selection_option"] = DEFAULT_SETTINGS["protein_selection_option"]
    overrides["ptm_selection_option"] = DEFAULT_SETTINGS["ptm_selection_option"]
    global_ptm_max = _to_int(_get_input_value(input, "settings_ptm_max_display"), DEFAULT_SETTINGS["ptm_max_display"])
    overrides["ptm_max_display"] = global_ptm_max
    overrides["use_original_protbox_size"] = _to_bool(
        _get_input_value(input, "settings_use_original_protbox_size"),
        DEFAULT_SETTINGS.get("use_original_protbox_size", True),
    )
    overrides["temporal_mode"] = _to_bool(
        _get_input_value(input, "settings_temporal_mode"),
        DEFAULT_SETTINGS.get("temporal_mode", False),
    )
    overrides["show_background_image"] = _to_bool(_get("show_background_image", _bool_default(cfg, "show_background_image")), _bool_default(cfg, "show_background_image"))
    overrides["display_types"] = list(DEFAULT_SETTINGS.get("display_types", []))
    overrides["show_groups"] = _to_bool(
        _get_input_value(input, "settings_show_groups"),
        _bool_default(cfg, "show_groups"),
    )
    if DEBUG_UI_ENABLED:
        overrides["debug_mode"] = _to_bool(
            _get_input_value(input, "settings_debug_mode"),
            False,
        )
    else:
        overrides["debug_mode"] = False
    overrides["show_multi_protein_indicator"] = _to_bool(
        _get_input_value(input, "settings_show_multi_protein_indicator"),
        _bool_default(cfg, "show_multi_protein_indicator"),
    )
    overrides["show_arrows"] = _to_bool(_get("show_arrows", _bool_default(cfg, "show_arrows")), _bool_default(cfg, "show_arrows"))
    overrides["show_text_boxes"] = _to_bool(_get("show_text_boxes", _bool_default(cfg, "show_text_boxes")), _bool_default(cfg, "show_text_boxes"))
    overrides["mode"] = str(cfg.get("mode", "analysis")).strip().lower()
    overrides["simple_kegg_mode"] = _to_bool(_get_input_value(input, _prefixed_id(prefix, "simple_kegg_mode")), prefix == "web")
    overrides["bg_scale"] = DEFAULT_WIKIPATHWAYS_BG_SCALE
    overrides["bg_offset_x"] = DEFAULT_WIKIPATHWAYS_BG_OFFSET_X
    overrides["bg_offset_y"] = DEFAULT_WIKIPATHWAYS_BG_OFFSET_Y
    legacy_neg = _hex_to_rgb(
        _get_input_value(input, "settings_negative_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]),
        DEFAULT_SETTINGS["negative_color"],
    )
    legacy_pos = _hex_to_rgb(
        _get_input_value(input, "settings_positive_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]),
        DEFAULT_SETTINGS["positive_color"],
    )
    legacy_max_neg = _to_float(
        _get_input_value(input, "settings_max_negative"), DEFAULT_SETTINGS["max_negative"]
    )
    legacy_max_pos = _to_float(
        _get_input_value(input, "settings_max_positive"), DEFAULT_SETTINGS["max_positive"]
    )
    stops_payload = _get_input_value(input, "settings_gradient_stops_json")
    parsed_stops: List[Dict[str, Any]] = []
    if isinstance(stops_payload, dict):
        parsed_stops = _normalize_gradient_stop_entries(stops_payload.get("rows"))
    elif isinstance(stops_payload, str) and stops_payload.strip():
        try:
            payload_obj = json.loads(stops_payload)
            if isinstance(payload_obj, dict):
                parsed_stops = _normalize_gradient_stop_entries(payload_obj.get("rows"))
            elif isinstance(payload_obj, list):
                parsed_stops = _normalize_gradient_stop_entries(payload_obj)
        except Exception:
            parsed_stops = []
    if len(parsed_stops) < 2:
        parsed_stops = _normalize_gradient_stop_entries(
            [
                {"value": legacy_max_neg, "color": list(legacy_neg)},
                {"value": 0.0, "color": [255, 255, 255]},
                {"value": legacy_max_pos, "color": list(legacy_pos)},
            ]
        )
    if len(parsed_stops) < 2:
        parsed_stops = _default_gradient_stops()
    min_stop = parsed_stops[0]
    max_stop = parsed_stops[-1]
    overrides["gradient_stops"] = parsed_stops
    overrides["negative_color"] = tuple(int(v) for v in list(min_stop.get("color", legacy_neg))[:3])
    overrides["positive_color"] = tuple(int(v) for v in list(max_stop.get("color", legacy_pos))[:3])
    overrides["max_negative"] = float(min_stop.get("value", legacy_max_neg))
    overrides["max_positive"] = float(max_stop.get("value", legacy_max_pos))
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
    current_mode = str(_get_input_value(input, "input_mode") or "user").strip().lower()
    overrides["use_black_protein_outlines"] = (
        True
        if current_mode == "demo"
        else _to_bool(
            _get_input_value(input, "settings_use_black_protein_outlines"),
            DEFAULT_SETTINGS.get("use_black_protein_outlines", False),
        )
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
KS_LAYOUT_GAP_X_DEFAULT = 240.0
KS_LAYOUT_CENTER_X_DEFAULT = 260.0
KS_LAYOUT_CENTER_Y_DEFAULT = 200.0
KS_CONCENTRIC_SPACE_DEFAULT = 70.0
KS_CONCENTRIC_RADIUS_DEFAULT = 220.0
KS_CONCENTRIC_ARROW_STOP_DEFAULT = 50.0


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
    gradient_stops: Optional[List[Dict[str, Any]]] = None,
) -> List[int]:
    try:
        if fold_value in (None, "", False):
            return [128, 128, 128]
        fold = float(fold_value)
    except (TypeError, ValueError):
        return [128, 128, 128]
    if not math.isfinite(fold):
        return [128, 128, 128]
    stops = _normalize_gradient_stop_entries(gradient_stops or [])
    if len(stops) >= 2:
        return _gradient_color_from_stops(float(fold), stops)
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


def _has_outline_columns(headers: Sequence[str]) -> bool:
    return any(str(header or "").strip().lower().startswith("o:") for header in list(headers or []))


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
    .gradient-table-wrap { margin-top: 8px; border: 1px solid #d5dde7; border-radius: 12px; overflow: hidden; background: #fff; }
    .gradient-table { width: 100%; border-collapse: collapse; table-layout: fixed; }
    .gradient-table thead th { background: #f8fafc; color: #334155; font-size: 12px; font-weight: 700; padding: 8px 10px; border-bottom: 1px solid #e5e7eb; }
    .gradient-table tbody td { padding: 8px 10px; border-bottom: 1px solid #eef2f7; vertical-align: middle; }
    .gradient-table tbody tr:last-child td { border-bottom: none; }
    .gradient-table .mk-grad-value { width: 100%; height: 34px; border: 1px solid #cbd5e1; border-radius: 8px; padding: 0 8px; font-size: 13px; }
    .gradient-table .mk-grad-color { width: 100%; height: 34px; border: 1px solid #cbd5e1; border-radius: 8px; padding: 0; background: transparent; }
    .mk-export-spinner { display: inline-flex; align-items: center; gap: 6px; font-size: 12px; color: #1f2937; }
    .mk-export-spinner::before { content: ""; width: 14px; height: 14px; border: 2px solid #cbd5e1;
        border-top-color: #1f2937; border-radius: 50%; animation: mk-spin 0.8s linear infinite; }
    .web-load-spinner { display: none; align-items: center; gap: 6px; font-size: 12px; color: #1f2937; white-space: nowrap; }
    .web-load-spinner::before { content: ""; width: 14px; height: 14px; border: 2px solid #cbd5e1;
        border-top-color: #1f2937; border-radius: 50%; animation: mk-spin 0.8s linear infinite; }
    .web-load-spinner.active { display: inline-flex; }
    .web-load-error { display: inline-flex; align-items: center; font-size: 12px; font-weight: 700; color: #b91c1c; max-width: 520px;
        line-height: 1.35; }
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
    .ks-evidence-toggle { margin-bottom: 8px; }
    .ks-evidence-toggle .shiny-input-checkboxgroup { margin-bottom: 0; }
    .ks-evidence-toggle .shiny-options-group { display: flex; flex-wrap: wrap; gap: 8px; }
    .ks-evidence-toggle .checkbox { margin: 0; }
    .ks-evidence-toggle .checkbox label {
        display: inline-flex; align-items: center; gap: 6px; padding: 8px 12px;
        border: 1px solid #d4d7dd; border-radius: 10px; background: #f8fafc; cursor: pointer;
        font-weight: 600; color: #1f2937; transition: background 0.18s ease, border-color 0.18s ease;
    }
    .ks-evidence-toggle .checkbox label:hover { background: #eef2ff; border-color: #b7c4f6; }
    .ks-evidence-toggle input[type="checkbox"] { margin: 0; accent-color: #2563eb; }
    .ks-filter-row { display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 8px; }
    .ks-filter-row > * { flex: 1; min-width: 140px; }
    .ks-filter-disabled { opacity: 0.45; pointer-events: none; }
    .pathway-table { width: 100%; table-layout: fixed; }
    .pathway-table th:nth-child(1),
    .pathway-table td:nth-child(1) { width: 110px; }
    .pathway-table th:nth-child(2),
    .pathway-table td:nth-child(2) { width: 280px; }
    .pathway-table th:nth-child(3),
    .pathway-table td:nth-child(3),
    .pathway-table th:nth-child(4),
    .pathway-table td:nth-child(4) { width: 120px; }
    .pathway-table td:nth-child(2) > div {
        display: block;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }
    .pathway-table td:nth-child(2) > small {
        display: block;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }
    .pathway-table-hover-tooltip { display: none; position: fixed; left: 0; top: 0;
        max-width: 420px; padding: 8px 10px; border-radius: 10px; background: rgba(15, 23, 42, 0.96);
        color: #ffffff; font-size: 12px; line-height: 1.35; box-shadow: 0 14px 30px rgba(15, 23, 42, 0.28); z-index: 1200;
        pointer-events: none; white-space: normal; word-break: break-word; }
    .input-page-stack { display: flex; flex-direction: column; row-gap: 16px; width: 100%; }
    .input-page-row { width: 100%; }
    .mk-input-busy-overlay { position: fixed; inset: 0; z-index: 2200; display: none; align-items: center; justify-content: center;
        background: rgba(11, 31, 51, 0.28); backdrop-filter: blur(2px); }
    .mk-input-busy-overlay.is-visible { display: flex; }
    .mk-input-busy-card { min-width: 280px; max-width: 360px; padding: 22px 24px; border-radius: 18px; background: rgba(255,255,255,0.98);
        border: 1px solid #dbe2ea; box-shadow: 0 28px 60px rgba(15, 23, 42, 0.20); display: flex; align-items: center; gap: 14px; }
    .mk-input-busy-spinner { width: 28px; height: 28px; border: 3px solid #dbe5f0; border-top-color: #2563eb; border-radius: 50%;
        animation: mk-spin 0.8s linear infinite; flex: 0 0 auto; }
    .mk-input-busy-copy { min-width: 0; }
    .mk-input-busy-title { font-size: 15px; font-weight: 800; color: #0f172a; line-height: 1.2; }
    .mk-input-busy-text { margin-top: 4px; font-size: 12px; color: #475569; line-height: 1.45; }
    .input-run-btn {
        width: 100%;
        min-height: 48px;
        border: 0;
        border-radius: 14px;
        background: linear-gradient(180deg, #22c55e 0%, #16a34a 100%);
        color: #ffffff;
        font-size: 16px;
        font-weight: 800;
        letter-spacing: 0.01em;
        box-shadow: 0 12px 28px rgba(22, 163, 74, 0.24);
    }
    .input-run-btn:hover,
    .input-run-btn:focus {
        background: linear-gradient(180deg, #1fb956 0%, #15803d 100%);
        color: #ffffff;
    }
    .input-run-btn:disabled {
        opacity: 0.55;
        cursor: not-allowed;
        box-shadow: none;
    }
    .input-run-status {
        margin-top: 10px;
        font-size: 12px;
        line-height: 1.45;
        color: #475569;
        min-height: 34px;
        white-space: normal;
    }
    .input-run-status.is-error {
        color: #b91c1c;
        font-weight: 700;
    }
    .input-sample-download-row { display: inline-flex; align-items: center; gap: 8px; margin-top: 2px; }
    .input-sample-download-btn {
        width: 34px; height: 34px; min-width: 34px; padding: 0; border-radius: 10px;
        display: inline-flex; align-items: center; justify-content: center;
    }
    .input-sample-download-label { font-size: 12px; font-weight: 600; color: #334155; line-height: 1.2; }
    .input-preview-labeled-row { display: flex; align-items: stretch; gap: 14px; width: 100%; }
    .input-preview-side-label {
        flex: 0 0 160px;
        min-width: 160px;
        display: flex;
        align-items: center;
        justify-content: center;
        border-radius: 16px;
        border: 1px solid #d7e3f4;
        background: linear-gradient(180deg, #f8fbff 0%, #e8f0ff 100%);
        box-shadow: 0 12px 28px rgba(15, 23, 42, 0.08);
        color: #0f172a;
        font-size: 16px;
        font-weight: 800;
        letter-spacing: 0.12em;
        text-transform: uppercase;
    }
    .input-preview-side-label.is-protein { color: #1d4ed8; }
    .input-preview-side-label.is-ptm { color: #0f766e; }
    .input-preview-side-label.is-metabolite { color: #7c3aed; }
    .input-preview-table-wrap { flex: 1 1 auto; min-width: 0; }
    .input-preview-section { display: flex; flex-wrap: wrap; gap: 16px; align-items: flex-start; width: 100%; }
    .input-preview-panel { flex: 0 0 auto; min-width: 0; max-width: 100%; }
    .input-preview-wrap { display: inline-block; border: 1px solid #dbe2ea; border-radius: 14px; background: #ffffff; overflow-x: auto;
        box-shadow: 0 12px 28px rgba(15, 23, 42, 0.08); }
    .input-preview-table { margin: 0; min-width: max-content; width: auto; table-layout: fixed; }
    .input-preview-table thead th { background: #f8fafc; color: #0f172a; font-weight: 700; white-space: nowrap; }
    .input-preview-table thead th.input-preview-header-error { color: #b91c1c; }
    .input-preview-table tbody td { color: #334155; white-space: nowrap; }
    .input-preview-table th, .input-preview-table td { padding: 10px 12px; border-color: #e5e7eb;
        overflow: hidden; text-overflow: ellipsis; }
    .input-preview-guide-card { width: 1380px; min-width: 1080px; height: 252px; border: 1px solid #dbe2ea;
        border-radius: 16px; background: linear-gradient(180deg, #f8fbff 0%, #eef5ff 100%);
        box-shadow: 0 18px 36px rgba(15, 23, 42, 0.10); overflow: hidden; }
    .input-preview-guide-layout { display: block; height: 100%; }
    .input-preview-guide-body { display: grid; grid-template-columns: repeat(8, minmax(0, 1fr)); gap: 12px; padding: 16px 18px 16px 18px; height: 100%; }
    .input-preview-guide-pill { background: rgba(255,255,255,0.94); border: 1px solid #d7e3f4; border-radius: 12px; padding: 12px 14px; }
    .input-preview-guide-pill-title { font-size: 13px; font-weight: 700; color: #0f172a; margin-bottom: 6px; }
    .input-preview-guide-required { color: #dc2626; margin-left: 3px; font-weight: 800; }
    .input-preview-guide-pill-text { font-size: 12px; line-height: 1.45; color: #334155; }
    .mk-mode-help-row { position: relative; overflow: visible; }
    .mk-mode-help-row .shiny-input-container { margin-bottom: 0; overflow: visible; }
    .mk-mode-help-row .form-check { overflow: visible; }
    .mk-mode-help-row label { overflow: visible; }
    .mk-inline-help-wrap { position: relative; display: inline-flex; align-items: center; vertical-align: middle; margin-left: 6px; overflow: visible; }
    .mk-inline-help { display: inline-flex; align-items: center; justify-content: center; width: 18px; height: 18px;
        border-radius: 999px; background: #2563eb; color: #ffffff; font-size: 12px; font-weight: 800; line-height: 1;
        cursor: help; user-select: none; box-shadow: 0 2px 6px rgba(37, 99, 235, 0.28); }
    .mk-floating-help-tooltip { display: none; position: fixed; left: 0; top: 0;
        min-width: 310px; max-width: 360px; padding: 10px 12px; border-radius: 12px; background: rgba(15, 23, 42, 0.96);
        color: #ffffff; font-size: 12px; line-height: 1.45; box-shadow: 0 18px 38px rgba(15, 23, 42, 0.28); z-index: 1200; }
    .mk-floating-help-tooltip::before { content: ""; position: absolute; left: -6px; top: 50%; width: 12px; height: 12px;
        background: rgba(15, 23, 42, 0.96); transform: translateY(-50%) rotate(45deg); }
    .web-load-ready { background-color: #16a34a !important; border-color: #16a34a !important; color: #ffffff !important; }
    .web-load-ready:hover, .web-load-ready:focus { background-color: #15803d !important; border-color: #15803d !important; color: #ffffff !important; }
    .pathway-viewer-card { position: relative; }
    .viewer-fullscreen-btn, .viewer-overlay-btn { border: 1px solid rgba(107, 114, 128, 0.45);
        background: rgba(107, 114, 128, 0.22); color: #111827; border-radius: 10px; padding: 6px 10px; font-size: 12px;
        font-weight: 600; backdrop-filter: blur(4px); transition: background 0.18s ease, border-color 0.18s ease, color 0.18s ease; }
    .viewer-fullscreen-btn { position: absolute; top: 10px; right: 10px; z-index: 45; }
    .viewer-fullscreen-btn:hover, .viewer-fullscreen-btn:focus, .viewer-overlay-btn:hover, .viewer-overlay-btn:focus {
        background: rgba(75, 85, 99, 0.82); border-color: rgba(75, 85, 99, 0.95); color: #ffffff; }
    .viewer-overlay-panel { position: absolute; top: 52px; right: 10px; z-index: 15; display: flex; flex-direction: column; align-items: flex-end; gap: 5px; }
    .viewer-overlay-settings { min-width: 92px; padding: 6px 12px; font-size: 12px !important; line-height: 1.2; letter-spacing: 0; }
    .viewer-overlay-panel.is-open .viewer-overlay-settings { background: rgba(75, 85, 99, 0.82); border-color: rgba(75, 85, 99, 0.95); color: #ffffff; }
    .viewer-overlay-controls { display: none; flex-direction: column; gap: 8px; align-items: flex-end; }
    .viewer-overlay-panel.is-open .viewer-overlay-controls { display: flex; }
    .viewer-overlay-btn { min-width: 36px; min-height: 32px; display: inline-flex; align-items: center; justify-content: center; padding: 4px 8px; font-size: 16px; line-height: 1; }
    .viewer-overlay-btn.is-off { background: rgba(31, 41, 55, 0.72); border-color: rgba(31, 41, 55, 0.9); color: #ffffff; }
    .viewer-overlay-btn.viewer-overlay-text { font-size: 13px; font-weight: 700; }
    .viewer-overlay-icon, .viewer-create-icon, .viewer-create-select-icon { width: 18px; height: 18px; display: block; color: currentColor; }
    .viewer-create-icon { width: 20px; height: 20px; }
    .viewer-create-panel { position: relative; z-index: 24; width: 100%; max-width: 100%; min-width: 0;
        margin-bottom: 12px; overflow: visible; }
    .viewer-create-scroll { width: 100%; overflow-x: auto; overflow-y: visible; scrollbar-width: thin; }
    .viewer-create-card { display: flex; align-items: flex-end; gap: 12px; flex-wrap: nowrap; width: max-content; min-width: 100%;
        background: rgba(255,255,255,0.94); border: 1px solid rgba(203, 213, 225, 0.9); border-radius: 16px;
        box-shadow: none; backdrop-filter: blur(10px); padding: 10px 12px; }
    .viewer-create-group { display: flex; flex-direction: column; gap: 5px; min-width: 0; flex: 0 0 auto; }
    .viewer-create-group-label { display: none; }
    .viewer-create-group-label.viewer-create-main-label { display: block; font-size: 11px; font-weight: 800; letter-spacing: 0.08em;
        text-transform: uppercase; color: #64748b; margin-bottom: 2px; }
    .viewer-create-section { display: flex; flex-direction: column; gap: 5px; min-width: 0; flex: 0 0 auto; margin-left: 4px; padding-left: 12px;
        border-left: 1px solid #e2e8f0; }
    .viewer-create-section.is-disabled { opacity: 0.5; }
    .viewer-create-section-label { display: block; font-size: 11px; font-weight: 800; letter-spacing: 0.08em; text-transform: uppercase; color: #64748b; margin-bottom: 2px; }
    .viewer-create-inline { display: flex; align-items: center; gap: 8px; min-width: 0; }
    .viewer-create-select-wrap { position: relative; display: inline-flex; align-items: center; }
    .viewer-create-select-wrap .viewer-create-select-icon { position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%); pointer-events: none; color: #0f172a; }
    .viewer-create-select-wrap .viewer-create-select-icon.viewer-create-shape-icon { width: 24px; height: 24px; }
    .viewer-create-field { height: 34px; border: 1px solid #cbd5e1; border-radius: 12px; background: #ffffff; color: #0f172a;
        font-size: 12px; padding: 0 10px; box-sizing: border-box; }
    .viewer-create-field:focus { outline: none; border-color: #94a3b8; box-shadow: 0 0 0 3px rgba(148, 163, 184, 0.18); }
    .viewer-create-field.viewer-create-search { width: 190px; }
    .viewer-create-field.viewer-create-select { min-width: 51px; width: 51px; padding-left: 22px; padding-right: 10px; color: transparent; text-shadow: none; }
    .viewer-create-select-wrap.is-active .viewer-create-field.viewer-create-select { background: linear-gradient(180deg, #dbeafe, #bfdbfe);
        border-color: #60a5fa; box-shadow: 0 0 0 3px rgba(96, 165, 250, 0.18); }
    .viewer-create-select-wrap.is-active .viewer-create-select-icon { color: #1d4ed8; }
    .viewer-create-field.viewer-create-select-interaction { width: 51px; min-width: 51px; }
    .viewer-create-field.viewer-create-select-shape { width: 51px; min-width: 51px; }
    .viewer-create-field.viewer-create-select option { color: #0f172a; }
    .viewer-create-field.viewer-create-select-protein { min-width: 220px; max-width: 260px; }
    .viewer-create-action { height: 34px; border: 1px solid #cbd5e1; border-radius: 12px; background: linear-gradient(180deg, #f8fafc, #eef2f7);
        color: #0f172a; padding: 0 12px; font-size: 12px; font-weight: 700; white-space: nowrap; cursor: pointer;
        transition: background 0.16s ease, border-color 0.16s ease, transform 0.16s ease; }
    .viewer-create-action:hover, .viewer-create-action:focus { background: linear-gradient(180deg, #eef2f7, #e2e8f0); border-color: #94a3b8; transform: translateY(-1px); }
    .viewer-create-action:disabled { opacity: 0.55; cursor: default; transform: none; }
    .viewer-create-action.viewer-create-icon-action { width: 44px; min-width: 44px; padding: 0; display: inline-flex; align-items: center; justify-content: center; }
    .viewer-create-action.viewer-create-mode-btn.is-active { background: linear-gradient(180deg, #dbeafe, #bfdbfe); border-color: #60a5fa; color: #1d4ed8; }
    .viewer-create-icon-pair { position: relative; width: 20px; height: 20px; display: inline-flex; align-items: center; justify-content: center; }
    .viewer-create-icon-pair .viewer-create-icon { position: absolute; inset: 0; margin: auto; }
    .viewer-create-icon-secondary { display: none; }
    .viewer-create-action:disabled .viewer-create-icon-primary { display: none; }
    .viewer-create-action:disabled .viewer-create-icon-secondary { display: block; }
    .viewer-create-action.viewer-create-snap-btn .viewer-create-icon-secondary { display: none; }
    .viewer-create-action.viewer-create-snap-btn.is-active .viewer-create-icon-primary { display: block; }
    .viewer-create-action.viewer-create-snap-btn.is-active .viewer-create-icon-secondary { display: none; }
    .viewer-create-action.viewer-create-snap-btn:not(.is-active) .viewer-create-icon-primary { display: none; }
    .viewer-create-action.viewer-create-snap-btn:not(.is-active) .viewer-create-icon-secondary { display: block; }
    .viewer-create-item { position: relative; }
    .viewer-create-popover { display: none; position: absolute; top: calc(100% + 8px); left: 0; z-index: 40; min-width: 150px; background: rgba(255,255,255,0.98);
        border: 1px solid #d7dee8; border-radius: 14px; box-shadow: 0 18px 44px rgba(15, 23, 42, 0.18); backdrop-filter: blur(12px); padding: 8px; }
    .viewer-create-item.is-open .viewer-create-popover { display: block; }
    .viewer-create-option-stack { display: flex; flex-direction: column; gap: 6px; }
    .viewer-create-option { border: 1px solid #d7dee8; border-radius: 12px; background: #f8fafc; color: #0f172a; padding: 9px 10px;
        font-size: 12px; font-weight: 700; text-align: left; cursor: pointer; transition: background 0.16s ease, border-color 0.16s ease; }
    .viewer-create-option:hover, .viewer-create-option:focus { background: #eef2f7; border-color: #94a3b8; }
    .viewer-create-protein-popover-overlay { display: none; position: absolute; z-index: 40; min-width: 240px; max-width: 300px; max-height: 220px; overflow-y: auto; }
    .viewer-create-legend-popover-overlay { display: none; position: absolute; z-index: 50; min-width: 150px; }
    .pathway-viewer-body { min-width: 0; position: relative; }
    .viewer-create-panel.is-hidden-for-cst,
    .pathway-viewer-card:fullscreen .viewer-create-panel.is-hidden-for-cst { display: none !important; }
    .viewer-create-option.viewer-create-protein-option { width: 100%; display: flex; align-items: center; justify-content: space-between; gap: 10px; }
    .viewer-create-protein-label { min-width: 0; flex: 1 1 auto; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
    .viewer-create-protein-swatch { width: 12px; height: 12px; border-radius: 3px; border: 1px solid rgba(15, 23, 42, 0.16); flex: 0 0 auto; }
    .viewer-create-metabolite-swatch { width: 12px; height: 12px; border-radius: 999px; border: 1px solid rgba(15, 23, 42, 0.16); flex: 0 0 auto; }
    .pathway-viewer-card:fullscreen { background: #ffffff; padding: 16px; overflow: auto; box-sizing: border-box;
        width: 100vw; height: 100vh; max-width: none; display: flex; flex-direction: column; align-items: stretch; }
    .pathway-viewer-card:fullscreen > *:not(.viewer-fullscreen-btn):not(.viewer-overlay-panel) { width: 100%; max-width: none; }
    .pathway-viewer-card:fullscreen .viewer-create-panel { display: block !important; visibility: visible !important; flex: 0 0 auto;
        align-self: stretch; width: 100%; max-width: 100%; min-width: 0; }
    .pathway-viewer-card:fullscreen .pathway-viewer-body { flex: 1 1 auto; min-height: 0; display: flex; flex-direction: column; overflow: auto; }
    .pathway-viewer-card .svg-container #svgCanvas > svg { width: 100%; height: 100%; display: block; position: absolute; inset: 0; }
    .pathway-viewer-card:fullscreen .pathway-viewer-body > .shiny-html-output { width: 100%; min-width: 0; flex: 1 1 auto; }
    .pathway-viewer-card:fullscreen .svg-container { width: 100% !important; min-width: 100% !important; max-width: none !important; flex: 1 1 auto; min-height: 0; overflow: auto; }
    .pathway-viewer-card:fullscreen #svgCanvas { width: 100% !important; max-width: none !important; }
    .pathway-viewer-card:fullscreen .viewer-fullscreen-btn { top: 16px; right: 16px; z-index: 60; }
    .pathway-viewer-card:fullscreen .viewer-overlay-panel { top: 58px; right: 16px; }
    .pathway-viewer-card:fullscreen .cst-viewer-shell { flex: 1 1 auto; min-height: 0; display: flex; flex-direction: column; overflow: hidden; }
    .pathway-viewer-card:fullscreen .cst-viewer-stage { flex: 1 1 auto; min-height: 0; height: 100% !important; }
    .pathway-viewer-card:fullscreen .cst-viewer-viewport { height: 100%; max-height: 100%; overflow: auto; }
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
                if (val === "home" || val === "input") return;
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
            var mimeType = msg && msg.mime_type ? msg.mime_type : "application/json;charset=utf-8";
            var blob = new Blob([content], { type: mimeType });
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

INPUT_BUSY_SCRIPT = ui.tags.script(
    """
    (function(){
        let datasetLoadPending = false;
        let forcedBusy = false;
        function overlay(){
            return document.getElementById("mk-input-busy-overlay");
        }
        function titleEl(){
            return document.querySelector('#mk-input-busy-overlay .mk-input-busy-title');
        }
        function textEl(){
            return document.querySelector('#mk-input-busy-overlay .mk-input-busy-text');
        }
        function inputTabActive(){
            return !!document.querySelector('#bookmark_selector a.nav-link.active[data-value="input"]');
        }
        function showOverlay(){
            const el = overlay();
            if (!el) return;
            el.classList.add("is-visible");
        }
        function hideOverlay(){
            const el = overlay();
            if (!el) return;
            el.classList.remove("is-visible");
        }
        if (window.Shiny && Shiny.addCustomMessageHandler){
            Shiny.addCustomMessageHandler("set_input_busy_state", function(msg){
                forcedBusy = !!(msg && msg.active);
                const nextTitle = msg && msg.title ? String(msg.title) : '';
                const nextText = msg && msg.text ? String(msg.text) : '';
                const title = titleEl();
                const text = textEl();
                if (nextTitle && title) title.textContent = nextTitle;
                if (nextText && text) text.textContent = nextText;
                if (forcedBusy && inputTabActive()){
                    showOverlay();
                } else {
                    datasetLoadPending = false;
                    hideOverlay();
                }
            });
        }
        document.addEventListener("change", function(ev){
            const target = ev && ev.target;
            if (!target) return;
            if (target.id !== "input_protein_upload" && target.id !== "input_ptm_upload" && target.id !== "input_metabolite_upload") return;
            datasetLoadPending = !!(target.files && target.files.length);
            if (datasetLoadPending && inputTabActive()){
                showOverlay();
            } else if (!datasetLoadPending) {
                hideOverlay();
            }
        }, true);
        document.addEventListener("click", function(ev){
            const target = ev && ev.target && ev.target.closest ? ev.target.closest("#input_run_pipeline") : null;
            if (!target) return;
            if (target.disabled) return;
            forcedBusy = true;
            if (inputTabActive()){
                showOverlay();
            }
        }, true);
        document.addEventListener("shiny:busy", function(){
            if ((datasetLoadPending || forcedBusy) && inputTabActive()){
                showOverlay();
            }
        });
        document.addEventListener("shiny:idle", function(){
            datasetLoadPending = false;
            if (!forcedBusy){
                hideOverlay();
            }
        });
        document.addEventListener("shown.bs.tab", function(){
            if (!((datasetLoadPending || forcedBusy) && inputTabActive())){
                hideOverlay();
            }
        });
    })();
    """
)

INLINE_HELP_TOOLTIP_SCRIPT = ui.tags.script(
    """
    (function(){
        let tooltipEl = null;
        let activeTarget = null;
        function ensureTooltip(){
            if (tooltipEl) return tooltipEl;
            tooltipEl = document.createElement('div');
            tooltipEl.className = 'mk-floating-help-tooltip';
            tooltipEl.setAttribute('role', 'tooltip');
            document.body.appendChild(tooltipEl);
            return tooltipEl;
        }
        function hideTooltip(){
            const el = ensureTooltip();
            el.style.display = 'none';
            el.innerHTML = '';
            activeTarget = null;
        }
        function positionTooltip(target){
            const el = ensureTooltip();
            const rect = target.getBoundingClientRect();
            const margin = 12;
            const desiredLeft = rect.right + margin;
            const maxLeft = window.innerWidth - el.offsetWidth - 12;
            const left = Math.max(12, Math.min(desiredLeft, maxLeft));
            const centerTop = rect.top + (rect.height * 0.5) - (el.offsetHeight * 0.5);
            const top = Math.max(12, Math.min(centerTop, window.innerHeight - el.offsetHeight - 12));
            el.style.left = left + 'px';
            el.style.top = top + 'px';
        }
        function showTooltip(target){
            const html = target.getAttribute('data-help-tooltip-html');
            if (!html) return;
            const el = ensureTooltip();
            activeTarget = target;
            el.innerHTML = html;
            el.style.display = 'block';
            positionTooltip(target);
        }
        document.addEventListener('mouseenter', function(ev){
            const target = ev.target && ev.target.closest ? ev.target.closest('.mk-inline-help[data-help-tooltip-html]') : null;
            if (!target) return;
            showTooltip(target);
        }, true);
        document.addEventListener('mouseleave', function(ev){
            const target = ev.target && ev.target.closest ? ev.target.closest('.mk-inline-help[data-help-tooltip-html]') : null;
            if (!target) return;
            hideTooltip();
        }, true);
        document.addEventListener('focusin', function(ev){
            const target = ev.target && ev.target.closest ? ev.target.closest('.mk-inline-help[data-help-tooltip-html]') : null;
            if (!target) return;
            showTooltip(target);
        });
        document.addEventListener('focusout', function(ev){
            const target = ev.target && ev.target.closest ? ev.target.closest('.mk-inline-help[data-help-tooltip-html]') : null;
            if (!target) return;
            hideTooltip();
        });
        window.addEventListener('scroll', function(){
            if (activeTarget) positionTooltip(activeTarget);
        }, true);
        window.addEventListener('resize', function(){
            if (activeTarget) positionTooltip(activeTarget);
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

GRADIENT_TABLE_SCRIPT = ui.tags.script(
    """
    (function () {
        function setShinyInputValue(id, value, opts) {
            if (!window.Shiny) return false;
            if (typeof window.Shiny.setInputValue === 'function') {
                window.Shiny.setInputValue(id, value, opts || { priority: 'event' });
                return true;
            }
            if (typeof window.Shiny.onInputChange === 'function') {
                window.Shiny.onInputChange(id, value);
                return true;
            }
            return false;
        }
        function collectGradientRows() {
            const table = document.getElementById('settings_gradient_table_editor');
            if (!table) return [];
            const rows = [];
            const trs = table.querySelectorAll('tbody tr.mk-gradient-row');
            trs.forEach((tr) => {
                const valueEl = tr.querySelector('input.mk-grad-value');
                const colorEl = tr.querySelector('input.mk-grad-color');
                if (!valueEl || !colorEl) return;
                const rawValueText = String(valueEl.value || '').trim();
                if (!rawValueText) return;
                const rawValue = Number(rawValueText);
                const rawColor = String(colorEl.value || '').trim();
                if (!Number.isFinite(rawValue)) return;
                if (!/^#[0-9a-fA-F]{6}$/.test(rawColor)) return;
                rows.push({ value: rawValue, color: rawColor });
            });
            return rows;
        }
        function syncLegacyFromRows(rows) {
            if (!Array.isArray(rows) || rows.length < 2) return;
            const sorted = [...rows].sort((a, b) => Number(a.value) - Number(b.value));
            const minRow = sorted[0];
            const maxRow = sorted[sorted.length - 1];
            setShinyInputValue('settings_max_negative', Number(minRow.value), { priority: 'event' });
            setShinyInputValue('settings_max_positive', Number(maxRow.value), { priority: 'event' });
            setShinyInputValue('settings_negative_color', String(minRow.color || '#000000'), { priority: 'event' });
            setShinyInputValue('settings_positive_color', String(maxRow.color || '#000000'), { priority: 'event' });
        }
        function pushGradientRows() {
            const rows = collectGradientRows();
            syncLegacyFromRows(rows);
            setShinyInputValue(
                'settings_gradient_stops_json',
                { rows, ts: Date.now() },
                { priority: 'event' }
            );
        }
        window.mkGradientTableUpdate = pushGradientRows;
        document.addEventListener('change', function (evt) {
            const target = evt && evt.target;
            if (!target || !target.closest) return;
            if (target.closest('#settings_gradient_table_editor')) {
                pushGradientRows();
            }
        });
        document.addEventListener('DOMContentLoaded', function () {
            pushGradientRows();
        });
        setTimeout(pushGradientRows, 300);
    })();
    """
)

SVG_OUTLINE_FONT_BOOTSTRAP_SCRIPT = ui.tags.script(
    f"""
    window.__mkOutlineFontData = {{
        regularB64: {json.dumps(OUTLINE_FONT_REGULAR_B64)},
        boldB64: {json.dumps(OUTLINE_FONT_BOLD_B64)},
        source: "local-segoe"
    }};
    """
)

SVG_DOWNLOAD_SCRIPT = ui.tags.script(
    """
    (function(){
        function ensureCitationPopup(){
            var existing = document.getElementById("mk-citation-popup-overlay");
            if (existing) return existing;
            var overlay = document.createElement("div");
            overlay.id = "mk-citation-popup-overlay";
            overlay.setAttribute("aria-hidden", "true");
            overlay.style.cssText = [
                "position:fixed",
                "inset:0",
                "background:rgba(15,23,42,0.52)",
                "display:none",
                "align-items:center",
                "justify-content:center",
                "z-index:9999",
                "padding:16px"
            ].join(";");
            var panel = document.createElement("div");
            panel.style.cssText = [
                "position:relative",
                "width:min(680px,calc(100vw - 24px))",
                "background:linear-gradient(180deg,#ffffff 0%,#f8fbff 100%)",
                "border:1px solid #cbd5e1",
                "border-radius:16px",
                "box-shadow:0 30px 70px rgba(15,23,42,0.28)",
                "padding:18px 20px 18px 20px",
                "font-family:Segoe UI,Tahoma,sans-serif",
                "color:#0f172a"
            ].join(";");
            var closeBtn = document.createElement("button");
            closeBtn.type = "button";
            closeBtn.setAttribute("aria-label", "Close citation popup");
            closeBtn.textContent = "×";
            closeBtn.style.cssText = [
                "position:absolute",
                "top:10px",
                "right:10px",
                "width:30px",
                "height:30px",
                "border:1px solid #cbd5e1",
                "background:#ffffff",
                "border-radius:999px",
                "font-size:20px",
                "line-height:1",
                "cursor:pointer",
                "color:#334155"
            ].join(";");
            closeBtn.addEventListener("click", function(){
                overlay.style.display = "none";
                overlay.setAttribute("aria-hidden", "true");
            });
            var heading = document.createElement("div");
            heading.textContent = "Please Cite MapKinase";
            heading.style.cssText = "font-size:22px;font-weight:800;margin:0 34px 8px 0;color:#0b1f33;";
            var sub = document.createElement("div");
            sub.textContent = "Your file download is in progress. Please include this citation in manuscripts and presentations.";
            sub.style.cssText = "font-size:14px;line-height:1.45;margin:0 0 12px 0;color:#334155;";
            var citationWrap = document.createElement("div");
            citationWrap.style.cssText = [
                "background:#f1f5f9",
                "border:1px solid #dbeafe",
                "border-left:4px solid #2563eb",
                "border-radius:12px",
                "padding:12px 13px"
            ].join(";");
            var citationTitle = document.createElement("div");
            citationTitle.textContent = "bioRxiv";
            citationTitle.style.cssText = "font-size:12px;font-weight:700;letter-spacing:0.02em;text-transform:uppercase;color:#1d4ed8;margin:0 0 6px 0;";
            var citationText = document.createElement("div");
            citationText.textContent = "(placeholder)";
            citationText.style.cssText = "font-size:14px;line-height:1.5;color:#0f172a;word-break:break-word;";
            citationWrap.appendChild(citationTitle);
            citationWrap.appendChild(citationText);
            panel.appendChild(closeBtn);
            panel.appendChild(heading);
            panel.appendChild(sub);
            panel.appendChild(citationWrap);
            overlay.appendChild(panel);
            overlay.addEventListener("click", function(ev){
                if (ev.target === overlay){
                    overlay.style.display = "none";
                    overlay.setAttribute("aria-hidden", "true");
                }
            });
            document.addEventListener("keydown", function(ev){
                if (ev.key === "Escape" && overlay.style.display !== "none"){
                    overlay.style.display = "none";
                    overlay.setAttribute("aria-hidden", "true");
                }
            });
            document.body.appendChild(overlay);
            return overlay;
        }
        function showCitationPopup(){
            var overlay = ensureCitationPopup();
            overlay.style.display = "flex";
            overlay.setAttribute("aria-hidden", "false");
        }
        function downloadHref(href, filename){
            var link = document.createElement("a");
            link.href = href;
            link.download = filename;
            document.body.appendChild(link);
            link.click();
            link.remove();
        }
        function blobToDataUrl(blob){
            return new Promise(function(resolve, reject){
                var reader = new FileReader();
                reader.onload = function(){ resolve(reader.result); };
                reader.onerror = function(err){ reject(err || new Error("Failed to read blob")); };
                reader.readAsDataURL(blob);
            });
        }
        function coerceNumber(value, fallback){
            var num = Number(value);
            return Number.isFinite(num) ? num : fallback;
        }
        function htmlToExportLines(html){
            var source = String(html || "");
            if (!source) return [];
            var normalized = source
                .replace(/<br\\s*\\/?>/gi, "\\n")
                .replace(/<\\/div\\s*>/gi, "\\n")
                .replace(/<div[^>]*>/gi, "")
                .replace(/<\\/p\\s*>/gi, "\\n")
                .replace(/<p[^>]*>/gi, "")
                .replace(/<\\/li\\s*>/gi, "\\n")
                .replace(/<li[^>]*>/gi, "• ");
            var tmp = document.createElement("div");
            tmp.innerHTML = normalized;
            var text = (tmp.textContent || tmp.innerText || "").replace(/\\r\\n?/g, "\\n");
            return text
                .split("\\n")
                .map(function(line){ return line.replace(/\\s+/g, " ").trim(); })
                .filter(function(line){ return line.length > 0; });
        }
        function replaceForeignObjectsForExport(svg){
            if (!svg || !svg.querySelectorAll) return;
            var SVG_NS = "http://www.w3.org/2000/svg";
            var groups = Array.prototype.slice.call(svg.querySelectorAll('g[data-type="text-box"]'));
            groups.forEach(function(group){
                var foreign = group.querySelector('[data-role="text-fo"]');
                if (!foreign) return;
                var editor = foreign.querySelector('.mk-text-editor');
                if (!editor) {
                    foreign.remove();
                    return;
                }
                var lines = htmlToExportLines(editor.innerHTML);
                if (!lines.length) {
                    lines = htmlToExportLines(editor.textContent || "");
                }
                if (!lines.length) {
                    foreign.remove();
                    return;
                }
                var x = coerceNumber(foreign.getAttribute('x'), 0);
                var y = coerceNumber(foreign.getAttribute('y'), 0);
                var width = coerceNumber(foreign.getAttribute('width'), 160);
                var height = coerceNumber(foreign.getAttribute('height'), 80);
                var fontSize = coerceNumber(editor.style.fontSize, 14);
                var lineHeight = fontSize * 1.2;
                var align = String(editor.style.textAlign || 'left').toLowerCase();
                var justify = String(editor.style.justifyContent || 'flex-start').toLowerCase();
                var xPad = 4;
                var xPos = x + xPad;
                var anchor = 'start';
                if (align === 'center') {
                    xPos = x + width / 2;
                    anchor = 'middle';
                } else if (align === 'right' || align === 'end') {
                    xPos = x + width - xPad;
                    anchor = 'end';
                }
                var textHeight = Math.max(lineHeight, lines.length * lineHeight);
                var startY = y + fontSize;
                if (justify === 'center') {
                    startY = y + (height - textHeight) / 2 + fontSize;
                } else if (justify === 'flex-end' || justify === 'end') {
                    startY = y + height - textHeight + fontSize;
                }
                var textNode = document.createElementNS(SVG_NS, 'text');
                textNode.setAttribute('x', String(xPos));
                textNode.setAttribute('y', String(startY));
                textNode.setAttribute('fill', editor.style.color || '#000000');
                textNode.setAttribute('font-family', editor.style.fontFamily || 'Arial');
                textNode.setAttribute('font-size', String(fontSize));
                textNode.setAttribute('text-anchor', anchor);
                textNode.setAttribute('dominant-baseline', 'alphabetic');
                textNode.setAttribute('pointer-events', 'none');
                if (editor.style.fontWeight) textNode.setAttribute('font-weight', editor.style.fontWeight);
                if (editor.style.fontStyle) textNode.setAttribute('font-style', editor.style.fontStyle);
                if (editor.style.textDecoration && editor.style.textDecoration !== 'none') {
                    textNode.setAttribute('text-decoration', editor.style.textDecoration);
                }
                lines.forEach(function(line, idx){
                    var tspan = document.createElementNS(SVG_NS, 'tspan');
                    tspan.setAttribute('x', String(xPos));
                    if (idx === 0) {
                        tspan.setAttribute('dy', '0');
                    } else {
                        tspan.setAttribute('dy', String(lineHeight));
                    }
                    tspan.textContent = line;
                    textNode.appendChild(tspan);
                });
                group.appendChild(textNode);
                foreign.remove();
            });
        }
        function annotatePublicationLabelMetrics(sourceSvg, cloneSvg){
            if (!sourceSvg || !cloneSvg || !sourceSvg.querySelectorAll || !cloneSvg.querySelectorAll) return;
            var selector = "text.cst-protein-oval-label, text.cst-auto-ptm-site-label, text.cst-auto-ptm-symbol";
            var sourceNodes = Array.prototype.slice.call(sourceSvg.querySelectorAll(selector));
            var cloneNodes = Array.prototype.slice.call(cloneSvg.querySelectorAll(selector));
            var count = Math.min(sourceNodes.length, cloneNodes.length);
            for (var i = 0; i < count; i += 1){
                var src = sourceNodes[i];
                var dst = cloneNodes[i];
                if (!src || !dst || typeof src.getBBox !== "function") continue;
                try {
                    var box = src.getBBox();
                    if (!box) continue;
                    var x = Number(box.x || 0);
                    var y = Number(box.y || 0);
                    var w = Number(box.width || 0);
                    var h = Number(box.height || 0);
                    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(w) || !Number.isFinite(h)) continue;
                    dst.setAttribute("data-mk-bbox", x.toFixed(3) + "," + y.toFixed(3) + "," + w.toFixed(3) + "," + h.toFixed(3));
                } catch (_bboxErr) {}
            }
        }
        function serializeSvg(svg){
            if (!svg) return null;
            var clone = svg.cloneNode(true);
            if (!clone.getAttribute("xmlns")){
                clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
            }
            replaceForeignObjectsForExport(clone);
            annotatePublicationLabelMetrics(svg, clone);
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
        function inlineSvgImageResources(text){
            return new Promise(function(resolve){
                try {
                    var parser = new DOMParser();
                    var doc = parser.parseFromString(text, "image/svg+xml");
                    var svgEl = doc.documentElement;
                    if (!svgEl || svgEl.nodeName.toLowerCase() !== "svg"){
                        resolve(text);
                        return;
                    }
                    var images = Array.prototype.slice.call(svgEl.querySelectorAll("image"));
                    if (!images.length){
                        resolve(text);
                        return;
                    }
                    var tasks = images.map(function(node){
                        var rawHref = node.getAttribute("href") || node.getAttribute("xlink:href") || "";
                        rawHref = String(rawHref || "").trim();
                        if (!rawHref || rawHref.indexOf("data:") === 0 || rawHref.indexOf("blob:") === 0){
                            return Promise.resolve();
                        }
                        var resolvedHref = rawHref;
                        try {
                            resolvedHref = String(new URL(rawHref, window.location.href));
                        } catch (_urlErr) {}
                        return fetch(resolvedHref)
                            .then(function(resp){
                                if (!resp.ok) throw new Error("Failed to fetch SVG image resource: " + resolvedHref);
                                return resp.blob();
                            })
                            .then(function(blob){ return blobToDataUrl(blob); })
                            .then(function(dataUrl){
                                node.setAttribute("href", dataUrl);
                                node.setAttribute("xlink:href", dataUrl);
                            })
                            .catch(function(err){
                                console.warn("Failed to inline SVG image resource", resolvedHref, err);
                            });
                    });
                    Promise.all(tasks).then(function(){
                        try {
                            var serializer = new XMLSerializer();
                            resolve(serializer.serializeToString(svgEl));
                        } catch (_serializeErr){
                            resolve(text);
                        }
                    });
                } catch (_err){
                    resolve(text);
                }
            });
        }
        function decodeDataSvgMarkup(uri){
            var raw = String(uri || "").trim();
            if (!raw) return null;
            var prefix = "data:image/svg+xml";
            if (raw.slice(0, prefix.length).toLowerCase() !== prefix){
                return null;
            }
            var commaIndex = raw.indexOf(",");
            if (commaIndex < 0){
                return null;
            }
            var meta = String(raw.slice(prefix.length, commaIndex) || "");
            var payload = String(raw.slice(commaIndex + 1) || "");
            var isBase64 = /;base64/i.test(meta);
            try {
                if (isBase64){
                    var binary = atob(payload.replace(/\\s+/g, ""));
                    try {
                        if (window.TextDecoder){
                            var bytes = new Uint8Array(binary.length);
                            for (var i = 0; i < binary.length; i += 1){
                                bytes[i] = binary.charCodeAt(i);
                            }
                            return new TextDecoder("utf-8").decode(bytes);
                        }
                    } catch (_decodeErr) {}
                    return decodeURIComponent(escape(binary));
                }
                return decodeURIComponent(payload);
            } catch (_err){
                return null;
            }
        }
        function coerceSvgLength(value, fallback){
            if (value == null) return fallback;
            var num = Number(String(value).replace(/px$/i, "").trim());
            return Number.isFinite(num) ? num : fallback;
        }
        function parseSvgViewBox(value){
            var raw = String(value || "").trim();
            if (!raw){
                return { minX: 0, minY: 0, width: 24, height: 24 };
            }
            var parts = raw.split(/[\\s,]+/).filter(function(token){ return token.length > 0; });
            if (parts.length !== 4){
                return { minX: 0, minY: 0, width: 24, height: 24 };
            }
            var nums = parts.map(function(token){ return Number(token); });
            if (!nums.every(function(n){ return Number.isFinite(n); })){
                return { minX: 0, minY: 0, width: 24, height: 24 };
            }
            return { minX: nums[0], minY: nums[1], width: nums[2], height: nums[3] };
        }
        function vectorizeEmbeddedSvgImages(text){
            try {
                var parser = new DOMParser();
                var doc = parser.parseFromString(text, "image/svg+xml");
                var svgEl = doc.documentElement;
                if (!svgEl || String(svgEl.nodeName || "").toLowerCase() !== "svg"){
                    return text;
                }
                var SVG_NS = "http://www.w3.org/2000/svg";
                var XLINK_NS = "http://www.w3.org/1999/xlink";
                var images = Array.prototype.slice.call(svgEl.querySelectorAll("image"));
                images.forEach(function(imageNode){
                    var href = imageNode.getAttribute("href") || imageNode.getAttributeNS(XLINK_NS, "href") || imageNode.getAttribute("xlink:href") || "";
                    href = String(href || "").trim();
                    if (href.slice(0, 18).toLowerCase() !== "data:image/svg+xml"){
                        return;
                    }
                    var embeddedMarkup = decodeDataSvgMarkup(href);
                    if (!embeddedMarkup) return;
                    var embeddedDoc = parser.parseFromString(embeddedMarkup, "image/svg+xml");
                    var embeddedSvg = embeddedDoc && embeddedDoc.documentElement;
                    if (!embeddedSvg || String(embeddedSvg.nodeName || "").toLowerCase() !== "svg"){
                        return;
                    }
                    var imageX = coerceSvgLength(imageNode.getAttribute("x"), 0);
                    var imageY = coerceSvgLength(imageNode.getAttribute("y"), 0);
                    var imageWidth = coerceSvgLength(imageNode.getAttribute("width"), 0);
                    var imageHeight = coerceSvgLength(imageNode.getAttribute("height"), 0);
                    if (!(imageWidth > 0) || !(imageHeight > 0)){
                        return;
                    }
                    var vb = parseSvgViewBox(embeddedSvg.getAttribute("viewBox"));
                    if (!(vb.width > 0) || !(vb.height > 0)){
                        vb.width = 24;
                        vb.height = 24;
                    }
                    var scaleX = imageWidth / vb.width;
                    var scaleY = imageHeight / vb.height;
                    var groupNode = doc.createElementNS(SVG_NS, "g");
                    Array.prototype.slice.call(imageNode.attributes || []).forEach(function(attr){
                        var name = String(attr && attr.name || "");
                        if (!name) return;
                        if (name === "href" || name === "xlink:href" || name === "x" || name === "y" || name === "width" || name === "height"){
                            return;
                        }
                        groupNode.setAttribute(name, attr.value);
                    });
                    groupNode.setAttribute(
                        "transform",
                        "translate(" + imageX + " " + imageY + ") scale(" + scaleX + " " + scaleY + ") translate(" + (-vb.minX) + " " + (-vb.minY) + ")"
                    );
                    Array.prototype.slice.call(embeddedSvg.childNodes || []).forEach(function(child){
                        try {
                            groupNode.appendChild(doc.importNode(child, true));
                        } catch (_childErr) {}
                    });
                    if (groupNode.childNodes.length < 1){
                        return;
                    }
                    if (imageNode.parentNode){
                        imageNode.parentNode.replaceChild(groupNode, imageNode);
                    }
                });
                var serializer = new XMLSerializer();
                return serializer.serializeToString(svgEl);
            } catch (_err){
                return text;
            }
        }
        function makeOfficeSafeSvgText(text){
            return inlineSvgImageResources(text).then(function(inlinedText){
                return vectorizeEmbeddedSvgImages(inlinedText);
            }).catch(function(){
                return text;
            });
        }
        function parseInlineStyleValue(styleText, prop){
            var source = String(styleText || "");
            if (!source) return "";
            var rx = new RegExp("(?:^|;)\\s*" + prop.replace(/[.*+?^${}()|[\\]\\\\]/g, "\\$&") + "\\s*:\\s*([^;]+)", "i");
            var match = source.match(rx);
            return match && match[1] ? String(match[1]).trim() : "";
        }
        function ensureOpenTypeJs(){
            return new Promise(function(resolve, reject){
                if (window.opentype && typeof window.opentype.parse === "function"){
                    resolve(window.opentype);
                    return;
                }
                var existing = document.querySelector("script[data-mk-opentype]");
                if (existing && existing.dataset.loading === "1"){
                    existing.addEventListener("load", function(){
                        if (window.opentype && typeof window.opentype.parse === "function"){
                            resolve(window.opentype);
                        } else {
                            reject(new Error("opentype.js unavailable after load"));
                        }
                    });
                    existing.addEventListener("error", function(){ reject(new Error("Failed to load opentype.js")); });
                    return;
                }
                var script = document.createElement("script");
                script.src = "https://cdnjs.cloudflare.com/ajax/libs/opentype.js/1.3.4/opentype.min.js";
                script.async = true;
                script.dataset.mkOpentype = "1";
                script.dataset.loading = "1";
                script.onload = function(){
                    if (window.opentype && typeof window.opentype.parse === "function"){
                        resolve(window.opentype);
                    } else {
                        reject(new Error("opentype.js loaded but parse API missing"));
                    }
                };
                script.onerror = function(){ reject(new Error("Failed to load opentype.js")); };
                document.head.appendChild(script);
            });
        }
        function loadMkOutlineFonts(){
            if (window.__mkOutlineFontsPromise) return window.__mkOutlineFontsPromise;
            window.__mkOutlineFontsPromise = ensureOpenTypeJs().then(function(opentypeLib){
                var sourceData = window.__mkOutlineFontData || {};
                var regularB64 = String(sourceData.regularB64 || "").trim();
                var boldB64 = String(sourceData.boldB64 || "").trim();
                if (!regularB64 || !boldB64){
                    throw new Error("Local Segoe outline fonts were not found. Expected segoeui.ttf and segoeuib.ttf.");
                }
                var parseOne = function(b64){
                    var clean = String(b64 || "").replace(/\\s+/g, "");
                    if (!clean) throw new Error("Outline font payload was empty.");
                    var binary = atob(clean);
                    var bytes = new Uint8Array(binary.length);
                    for (var i = 0; i < binary.length; i += 1){
                        bytes[i] = binary.charCodeAt(i);
                    }
                    return opentypeLib.parse(bytes.buffer);
                };
                return { regular: parseOne(regularB64), bold: parseOne(boldB64) };
            });
            return window.__mkOutlineFontsPromise;
        }
        function fontWeightToNumber(value){
            var token = String(value || "").trim().toLowerCase();
            var numeric = Number(token);
            if (Number.isFinite(numeric)) return numeric;
            if (token === "bold") return 700;
            if (token === "bolder") return 800;
            if (token === "normal") return 400;
            if (token === "lighter") return 300;
            return 400;
        }
        function svgPathCommandsToD(commands){
            if (!Array.isArray(commands) || !commands.length) return "";
            var parts = [];
            commands.forEach(function(cmd){
                if (!cmd || !cmd.type) return;
                var t = String(cmd.type || "").toUpperCase();
                if (t === "M" || t === "L"){
                    parts.push(t + " " + Number(cmd.x || 0).toFixed(3) + " " + Number(cmd.y || 0).toFixed(3));
                    return;
                }
                if (t === "C"){
                    parts.push(
                        t + " "
                        + Number(cmd.x1 || 0).toFixed(3) + " " + Number(cmd.y1 || 0).toFixed(3) + " "
                        + Number(cmd.x2 || 0).toFixed(3) + " " + Number(cmd.y2 || 0).toFixed(3) + " "
                        + Number(cmd.x || 0).toFixed(3) + " " + Number(cmd.y || 0).toFixed(3)
                    );
                    return;
                }
                if (t === "Q"){
                    parts.push(
                        t + " "
                        + Number(cmd.x1 || 0).toFixed(3) + " " + Number(cmd.y1 || 0).toFixed(3) + " "
                        + Number(cmd.x || 0).toFixed(3) + " " + Number(cmd.y || 0).toFixed(3)
                    );
                    return;
                }
                if (t === "Z"){
                    parts.push("Z");
                }
            });
            return parts.join(" ");
        }
        function parseExportBBoxAttr(node){
            if (!node || !node.getAttribute) return null;
            var raw = String(node.getAttribute("data-mk-bbox") || "").trim();
            if (!raw) return null;
            var parts = raw.split(",");
            if (parts.length !== 4) return null;
            var nums = parts.map(function(v){ return Number(v); });
            if (!nums.every(function(v){ return Number.isFinite(v); })) return null;
            return { x: nums[0], y: nums[1], width: nums[2], height: nums[3] };
        }
        function outlinePublicationLabelsToPaths(text, fontPair, options){
            try {
                if (!fontPair || !fontPair.regular){
                    throw new Error("Publication outline fonts were not available.");
                }
                var debugOverlay = !!(options && options.debugOverlay);
                var parser = new DOMParser();
                var doc = parser.parseFromString(text, "image/svg+xml");
                var svgEl = doc.documentElement;
                if (!svgEl || String(svgEl.nodeName || "").toLowerCase() !== "svg"){
                    throw new Error("Invalid SVG payload for publication outlining.");
                }
                var SVG_NS = "http://www.w3.org/2000/svg";
                var selector = "text.cst-protein-oval-label, text.cst-auto-ptm-site-label, text.cst-auto-ptm-symbol";
                var nodes = Array.prototype.slice.call(svgEl.querySelectorAll(selector));
                if (!nodes.length) return text;
                var unsupportedFamilies = [];
                var convertedCount = 0;
                nodes.forEach(function(node){
                    var value = String(node.textContent || "").replace(/\\s+/g, " ").trim();
                    if (!value) return;
                    var className = String(node.getAttribute("class") || "");
                    var isProteinLabel = className.indexOf("cst-protein-oval-label") >= 0;
                    var isPtmSiteLabel = className.indexOf("cst-auto-ptm-site-label") >= 0;
                    var styleText = String(node.getAttribute("style") || "");
                    var fontSize = coerceSvgLength(node.getAttribute("font-size"), NaN);
                    if (!Number.isFinite(fontSize) || fontSize <= 0){
                        fontSize = coerceSvgLength(parseInlineStyleValue(styleText, "font-size"), 0);
                    }
                    if (!(fontSize > 0)) fontSize = 12;
                    var sourceFamily = String(
                        node.getAttribute("font-family")
                        || parseInlineStyleValue(styleText, "font-family")
                        || ""
                    );
                    var sourceFamilyNormalized = sourceFamily.toLowerCase();
                    if (
                        sourceFamilyNormalized
                        && sourceFamilyNormalized.indexOf("noto sans") < 0
                        && sourceFamilyNormalized.indexOf("arial") < 0
                        && sourceFamilyNormalized.indexOf("sans-serif") < 0
                    ){
                        unsupportedFamilies.push(sourceFamily || "(unknown)");
                        return;
                    }
                    var fontWeightValue = String(
                        node.getAttribute("font-weight")
                        || parseInlineStyleValue(styleText, "font-weight")
                        || "400"
                    );
                    var weight = fontWeightToNumber(fontWeightValue);
                    // CST protein/PTM labels are intended to render bold in publication exports.
                    var font = fontPair.bold || ((weight >= 600 && fontPair.bold) ? fontPair.bold : fontPair.regular);
                    var fill = String(
                        node.getAttribute("fill")
                        || parseInlineStyleValue(styleText, "fill")
                        || "#0f172a"
                    );
                    var x = coerceSvgLength(node.getAttribute("x"), 0);
                    var y = coerceSvgLength(node.getAttribute("y"), 0);
                    var anchor = String(
                        node.getAttribute("text-anchor")
                        || parseInlineStyleValue(styleText, "text-anchor")
                        || "start"
                    ).toLowerCase();
                    var dominant = String(
                        node.getAttribute("dominant-baseline")
                        || parseInlineStyleValue(styleText, "dominant-baseline")
                        || "alphabetic"
                    ).toLowerCase();
                    var bboxHint = parseExportBBoxAttr(node);
                    var startX = x;
                    var baselineY = y;
                    var transformPrefix = "";
                    var glyphPath;
                    if (bboxHint){
                        startX = 0;
                        glyphPath = font.getPath(value, startX, 0, fontSize, { kerning: true });
                        var hintedBox = glyphPath && typeof glyphPath.getBoundingBox === "function" ? glyphPath.getBoundingBox() : null;
                        if (hintedBox){
                            var hbX1 = Number(hintedBox.x1 || 0);
                            var hbY1 = Number(hintedBox.y1 || 0);
                            var hbW = Number(hintedBox.x2 || 0) - hbX1;
                            var hbH = Number(hintedBox.y2 || 0) - hbY1;
                            if (isPtmSiteLabel && hbW > 0.001 && hbH > 0.001 && Number(bboxHint.width || 0) > 0.001 && Number(bboxHint.height || 0) > 0.001){
                                var fitW = Number(bboxHint.width || 0);
                                var fitH = Number(bboxHint.height || 0);
                                var uniformScale = Math.min(fitW / hbW, fitH / hbH);
                                if (!Number.isFinite(uniformScale) || uniformScale <= 0){
                                    uniformScale = fitW / hbW;
                                }
                                // PTM labels look visually larger after bold outlining; trim slightly.
                                uniformScale = uniformScale * 0.94;
                                var scaledW = hbW * uniformScale;
                                var scaledH = hbH * uniformScale;
                                var txBox = Number(bboxHint.x || 0) + ((fitW - scaledW) * 0.5);
                                var tyBox = Number(bboxHint.y || 0) + ((fitH - scaledH) * 0.5);
                                transformPrefix =
                                    "translate(" + txBox.toFixed(3) + " " + tyBox.toFixed(3) + ") "
                                    + "scale(" + uniformScale.toFixed(6) + " " + uniformScale.toFixed(6) + ") "
                                    + "translate(" + (-hbX1).toFixed(3) + " " + (-hbY1).toFixed(3) + ")";
                            } else {
                                var dx = Number(bboxHint.x || 0) - hbX1;
                                var dy = Number(bboxHint.y || 0) - hbY1;
                                if (Number.isFinite(dx) && Number.isFinite(dy)){
                                    transformPrefix = "translate(" + dx.toFixed(3) + " " + dy.toFixed(3) + ")";
                                }
                            }
                        }
                    } else {
                        if (anchor === "middle"){
                            startX = x;
                        } else if (anchor === "end"){
                            startX = x;
                        }
                        var probePath = font.getPath(value, 0, 0, fontSize, { kerning: true });
                        var probeBbox = (probePath && typeof probePath.getBoundingBox === "function") ? probePath.getBoundingBox() : null;
                        if (probeBbox){
                            var centerX = (Number(probeBbox.x1 || 0) + Number(probeBbox.x2 || 0)) * 0.5;
                            if (anchor === "middle"){
                                startX = x - centerX;
                            } else if (anchor === "end"){
                                startX = x - Number(probeBbox.x2 || 0);
                            } else {
                                startX = x - Number(probeBbox.x1 || 0);
                            }
                            var topY = Number(probeBbox.y1 || 0);
                            var bottomY = Number(probeBbox.y2 || 0);
                            var centerY = (topY + bottomY) * 0.5;
                            if (dominant === "middle" || dominant === "central"){
                                baselineY = y - centerY;
                            } else if (dominant === "hanging" || dominant === "text-before-edge"){
                                baselineY = y - topY;
                            } else if (dominant === "ideographic" || dominant === "text-after-edge"){
                                baselineY = y - bottomY;
                            }
                        }
                        glyphPath = font.getPath(value, startX, baselineY, fontSize, { kerning: true });
                    }
                    if (!glyphPath || !Array.isArray(glyphPath.commands) || !glyphPath.commands.length){
                        return;
                    }
                    var d = svgPathCommandsToD(glyphPath.commands);
                    if (!d) return;
                    var proteinShiftY = isProteinLabel ? (fontSize * 0.36) : 0;
                    var opacity = Number(node.getAttribute("opacity"));
                    if (!Number.isFinite(opacity)){
                        opacity = Number(parseInlineStyleValue(styleText, "opacity"));
                    }
                    if (!Number.isFinite(opacity)) opacity = 1;
                    var pathNode = doc.createElementNS(SVG_NS, "path");
                    pathNode.setAttribute("d", d);
                    if (debugOverlay){
                        pathNode.setAttribute("fill", "#ef4444");
                        pathNode.setAttribute("opacity", "0.70");
                    } else {
                        pathNode.setAttribute("fill", fill || "#0f172a");
                        pathNode.setAttribute("opacity", Math.max(0, Math.min(1, opacity)).toFixed(3));
                    }
                    pathNode.setAttribute("stroke", "none");
                    pathNode.setAttribute("pointer-events", "none");
                    pathNode.setAttribute("data-mk-label-outlined", "1");
                    pathNode.setAttribute("data-role", String(node.getAttribute("data-role") || ""));
                    pathNode.setAttribute("data-node-id", String(node.getAttribute("data-node-id") || ""));
                    var ptmId = String(node.getAttribute("data-ptm-id") || "");
                    if (ptmId) pathNode.setAttribute("data-ptm-id", ptmId);
                    var nodeTransform = String(node.getAttribute("transform") || "").trim();
                    var shiftTransform = proteinShiftY ? ("translate(0 " + proteinShiftY.toFixed(3) + ")") : "";
                    var joinedTransform = [transformPrefix, shiftTransform, nodeTransform].filter(function(v){ return !!String(v || "").trim(); }).join(" ").trim();
                    if (joinedTransform){
                        pathNode.setAttribute("transform", joinedTransform);
                    }
                    if (node.parentNode){
                        if (debugOverlay){
                            node.setAttribute("fill", "#2563eb");
                            node.setAttribute("opacity", "1.000");
                            node.parentNode.insertBefore(pathNode, node.nextSibling);
                        } else {
                            node.parentNode.replaceChild(pathNode, node);
                        }
                        convertedCount += 1;
                    }
                });
                if (unsupportedFamilies.length){
                    var familyList = Array.from(new Set(unsupportedFamilies)).slice(0, 6).join(", ");
                    throw new Error("Unsupported label font family for outline conversion: " + familyList + ". Use a supported family or keep standard SVG.");
                }
                if (convertedCount < 1){
                    throw new Error("No publication labels were converted to outlines.");
                }
                var serializer = new XMLSerializer();
                return serializer.serializeToString(svgEl);
            } catch (_err){
                throw _err;
            }
        }
        function makePublicationSafeSvgText(text){
            return makeOfficeSafeSvgText(text).then(function(nextText){
                return loadMkOutlineFonts().then(function(fontPair){
                    return outlinePublicationLabelsToPaths(nextText, fontPair, { debugOverlay: false });
                });
            });
        }
        function makePublicationDebugSvgText(text){
            return makeOfficeSafeSvgText(text).then(function(nextText){
                return loadMkOutlineFonts().then(function(fontPair){
                    return outlinePublicationLabelsToPaths(nextText, fontPair, { debugOverlay: true });
                });
            }).catch(function(){
                return text;
            });
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
        function ensureCanvg(){
            return new Promise(function(resolve, reject){
                if (window.canvg && window.canvg.Canvg){
                    resolve(window.canvg);
                    return;
                }
                var existing = document.querySelector("script[data-mk-canvg]");
                if (existing && existing.dataset.loading === "1"){
                    existing.addEventListener("load", function(){ resolve(window.canvg); });
                    existing.addEventListener("error", function(){ reject(new Error("Failed to load canvg")); });
                    return;
                }
                var script = document.createElement("script");
                script.src = "https://cdnjs.cloudflare.com/ajax/libs/canvg/3.0.11/umd.min.js";
                script.async = true;
                script.dataset.mkCanvg = "1";
                script.dataset.loading = "1";
                script.onload = function(){ resolve(window.canvg); };
                script.onerror = function(){ reject(new Error("Failed to load canvg")); };
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
        function exportSerializedSvg(payload, format, name){
            if (!payload) return;
            var text = payload.text;
            var width = payload.width;
            var height = payload.height;
            var rasterScale = 4;
            var maxRasterPixels = format === "pdf" ? 24000000 : 48000000;
            if (format === "svg"){
                makeOfficeSafeSvgText(text).then(function(svgText){
                    var blob = new Blob([svgText], { type: "image/svg+xml" });
                    var href = URL.createObjectURL(blob);
                    downloadHref(href, name + ".svg");
                    setTimeout(function(){ URL.revokeObjectURL(href); }, 1500);
                }).catch(function(err){
                    console.warn("Office-safe SVG export fallback", err);
                    var blob = new Blob([text], { type: "image/svg+xml" });
                    var href = URL.createObjectURL(blob);
                    downloadHref(href, name + ".svg");
                    setTimeout(function(){ URL.revokeObjectURL(href); }, 1500);
                });
                return;
            }
            if (format === "svgprint"){
                makePublicationSafeSvgText(text).then(function(svgText){
                    var blob = new Blob([svgText], { type: "image/svg+xml" });
                    var href = URL.createObjectURL(blob);
                    downloadHref(href, name + "_print.svg");
                    setTimeout(function(){ URL.revokeObjectURL(href); }, 1500);
                }).catch(function(err){
                    console.error("Publication-safe SVG export failed", err);
                    try {
                        window.alert("SVG Print export failed: " + (err && err.message ? err.message : "outline conversion failed"));
                    } catch (_alertErr) {}
                });
                return;
            }
            if (format === "svgdebug"){
                if (window.MAPKINASE_DEBUG_EXPORTS !== true){
                    console.warn("SVG debug export is disabled in this deployment.");
                    return;
                }
                makePublicationDebugSvgText(text).then(function(svgText){
                    var blob = new Blob([svgText], { type: "image/svg+xml" });
                    var href = URL.createObjectURL(blob);
                    downloadHref(href, name + "_debug.svg");
                    setTimeout(function(){ URL.revokeObjectURL(href); }, 1500);
                }).catch(function(err){
                    console.error("Publication-debug SVG export failed", err);
                    try {
                        window.alert("SVG Debug export failed: " + (err && err.message ? err.message : "debug outline conversion failed"));
                    } catch (_alertErr) {}
                });
                return;
            }
            inlineSvgImageResources(text).then(function(inlinedText){
                var exportScale = rasterScale;
                var scaledWidth = Math.max(1, Math.round(width * exportScale));
                var scaledHeight = Math.max(1, Math.round(height * exportScale));
                var pixelCount = scaledWidth * scaledHeight;
                if (pixelCount > maxRasterPixels){
                    exportScale = Math.sqrt(maxRasterPixels / Math.max(width * height, 1));
                    scaledWidth = Math.max(1, Math.round(width * exportScale));
                    scaledHeight = Math.max(1, Math.round(height * exportScale));
                }
                var canvas = document.createElement("canvas");
                canvas.width = scaledWidth;
                canvas.height = scaledHeight;
                var ctx = canvas.getContext("2d");
                ensureCanvg().then(function(canvgLib){
                    if (!ctx || !canvgLib || !canvgLib.Canvg){
                        throw new Error("canvg unavailable");
                    }
                    ctx.save();
                    ctx.setTransform(1, 0, 0, 1, 0, 0);
                    ctx.fillStyle = "white";
                    ctx.fillRect(0, 0, canvas.width, canvas.height);
                    ctx.scale(exportScale, exportScale);
                    return canvgLib.Canvg.fromString(ctx, inlinedText, {
                        ignoreDimensions: true,
                        ignoreClear: true,
                    }).render().then(function(){
                        ctx.restore();
                        if (format === "pdf"){
                            ensureJsPdf().then(function(jsPDF){
                                if (!jsPDF){ console.error("jsPDF unavailable"); return; }
                                var orientation = width >= height ? "l" : "p";
                                var pdf = new jsPDF({ orientation: orientation, unit: "pt", format: [width, height] });
                                pdf.addImage(canvas.toDataURL("image/jpeg", 0.92), "JPEG", 0, 0, width, height);
                                pdf.save(name + ".pdf");
                            }).catch(function(err){ console.error("PDF export failed", err); });
                        } else {
                            var mime = format === "jpeg" ? "image/jpeg" : "image/png";
                            var ext = format === "jpeg" ? ".jpeg" : ".png";
                            var data = canvas.toDataURL(mime, 1.0);
                            downloadHref(data, name + ext);
                        }
                    });
                }).catch(function(err){
                    console.error("SVG raster export failed", err);
                });
            });
        }
        function exportSvg(svgEl, format, name){
            var payload = serializeSvg(svgEl);
            exportSerializedSvg(payload, format, name);
        }
        function findExportPayload(btn){
            if (!btn) return null;
            var scoped = btn.closest(".pathway-viewer-card");
            if (!scoped) return null;
            var svg = scoped.querySelector(".svg-container svg");
            if (svg){
                return { kind: "svg", value: svg };
            }
            var key = String(btn.getAttribute("data-mk-prefix") || ((scoped.dataset && scoped.dataset.prefix) || "")).toLowerCase();
            if (key && window.__mkExportSvgMap && typeof window.__mkExportSvgMap[key] === "function"){
                try {
                    var payload = window.__mkExportSvgMap[key]();
                    if (payload && payload.text){
                        return { kind: "payload", value: payload };
                    }
                } catch (err){
                    console.error("Custom export payload failed", err);
                }
            }
            svg = document.querySelector(".svg-container svg");
            if (svg){
                return { kind: "svg", value: svg };
            }
            return null;
        }
        document.addEventListener("click", function(ev){
            var btn = ev.target.closest("[data-mk-download]");
            if (!btn) return;
            ev.preventDefault();
            showCitationPopup();
            var format = btn.getAttribute("data-mk-download");
            var target = findExportPayload(btn);
            if (!target){
                console.warn("No export target found for download");
                return;
            }
            var name = btn.getAttribute("data-mk-name") || "pathway";
            if (target.kind === "payload"){
                exportSerializedSvg(target.value, format, name);
            } else {
                exportSvg(target.value, format, name);
            }
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
    fisher_btn_id = _prefixed_id(prefix, "run_fisher_test")
    fisher_settings_btn_id = _prefixed_id(prefix, "fisher_settings_btn")
    fisher_settings_popup_id = _prefixed_id(prefix, "fisher_settings_popup")
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
                    {"style": "display:flex; align-items:center; justify-content:space-between; gap:12px;"},
                    ui.tags.button(
                        ui.tags.i({"class": "fa fa-filter"}), "Filters",
                        {"id": filter_btn_id, "type": "button", "class": "ks-filter-button pathway-filter-button"}
                    ),
                    ui.div(
                        {"class": "pathway-action-buttons"},
                        ui.input_action_button(
                            _prefixed_id(prefix, "download_pathway_table"),
                            _icon_markup("download", "mk-inline-icon pathway-action-icon"),
                            class_="btn btn-outline-secondary btn-sm pathway-action-btn pathway-icon-btn",
                        ),
                        ui.input_action_button(
                            fisher_btn_id,
                            ui.span({"class": "pathway-run-label"}, "Fisher's Exact Test"),
                            class_="btn btn-primary btn-sm pathway-action-btn pathway-run-btn",
                        ),
                        ui.tags.button(
                            _icon_markup("settings", "mk-inline-icon pathway-action-icon"),
                            {
                                "id": fisher_settings_btn_id,
                                "type": "button",
                                "class": "btn btn-outline-secondary btn-sm pathway-action-btn pathway-icon-btn fisher-settings-trigger",
                                "title": "Fisher exact test settings",
                                "aria-label": "Fisher exact test settings",
                            },
                        ),
                    ),
                ),
                ui.div(
                    {"id": _prefixed_id(prefix, "fisher_run_progress"), "class": "fisher-run-progress"},
                    ui.div({"class": "fisher-run-progress-bar"}),
                    ui.div({"class": "fisher-run-progress-text"}, "Running Fisher's Exact Test..."),
                ),
                ui.output_ui(_prefixed_id(prefix, "fisher_run_state")),
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
                ui.div(
                    {"id": fisher_settings_popup_id, "class": "ks-filter-popup fisher-settings-popup"},
                    ui.div({"class": "fisher-settings-title"}, "Fisher's Exact Test Settings"),
                    ui.div(
                        {"class": "fisher-settings-copy"},
                        "Choose which direction of significance to test and set the fold-change cutoffs used for proteins and PTMs.",
                    ),
                    ui.div(
                        {"class": "ks-filter-row"},
                        ui.input_radio_buttons(
                            _prefixed_id(prefix, "fisher_sig_mode"),
                            "Significance Mode",
                            choices={
                                "both": "Positive + Negative",
                                "positive": "Positive only",
                                "negative": "Negative only",
                            },
                            selected="both",
                            inline=False,
                        ),
                    ),
                    ui.div(
                        {"class": "fisher-settings-thresholds"},
                        ui.input_numeric(
                            _prefixed_id(prefix, "fisher_sig_pos"),
                            "Positive cutoff",
                            value=1.5,
                            step=0.1,
                            width="140px",
                        ),
                        ui.input_numeric(
                            _prefixed_id(prefix, "fisher_sig_neg"),
                            "Negative cutoff",
                            value=-1.5,
                            step=0.1,
                            width="140px",
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
                    .pathway-action-buttons { display:flex; align-items:center; justify-content:flex-end; gap:8px; margin-left:auto; }
                    .pathway-action-btn {
                        height: 38px;
                        min-height: 38px;
                        display:flex;
                        align-items:center;
                        justify-content:center;
                        flex-wrap: nowrap;
                        gap: 8px;
                        font-weight: 600;
                    }
                    .pathway-icon-btn { width: 38px; min-width: 38px; padding: 0; }
                    .pathway-run-btn { min-width: 92px; padding: 0 12px; white-space: nowrap; }
                    .pathway-run-btn .btn { white-space: nowrap; }
                    .pathway-run-btn .shiny-bound-input,
                    .pathway-run-btn span,
                    .pathway-run-btn svg { vertical-align: middle; }
                    .pathway-run-label { white-space: nowrap; line-height: 1; display:inline-block; }
                    .pathway-action-icon { width: 17px; height: 17px; display:block; }
                    .fisher-run-progress {
                        display: none;
                        align-items: center;
                        gap: 10px;
                        margin-top: -2px;
                    }
                    .fisher-run-progress.active { display:flex; }
                    .fisher-run-progress-bar {
                        position: relative;
                        flex: 1 1 auto;
                        height: 6px;
                        border-radius: 999px;
                        overflow: hidden;
                        background: rgba(37, 99, 235, 0.14);
                    }
                    .fisher-run-progress-bar::before {
                        content: "";
                        position: absolute;
                        inset: 0;
                        width: 42%;
                        border-radius: inherit;
                        background: linear-gradient(90deg, #2563eb 0%, #60a5fa 100%);
                        animation: fisher-run-slide 1.1s ease-in-out infinite;
                    }
                    .fisher-run-progress-text {
                        flex: 0 0 auto;
                        color: #1f4fa3;
                        font-size: 0.84rem;
                        font-weight: 600;
                        white-space: nowrap;
                    }
                    @keyframes fisher-run-slide {
                        0% { transform: translateX(-120%); }
                        100% { transform: translateX(260%); }
                    }
                    .fisher-settings-popup {
                        min-width: 320px;
                        max-width: 420px;
                        padding: 14px 16px 12px;
                        border-radius: 14px;
                        box-shadow: 0 16px 34px rgba(15, 23, 42, 0.16);
                    }
                    .fisher-settings-title { font-size: 0.96rem; font-weight: 700; color: #111827; margin-bottom: 4px; }
                    .fisher-settings-copy { font-size: 0.84rem; color: #4b5563; line-height: 1.4; margin-bottom: 12px; }
                    .fisher-settings-popup .shiny-input-radiogroup > label { font-weight: 700; color: #1f2937; margin-bottom: 8px; }
                    .fisher-settings-popup .radio { margin-bottom: 6px; }
                    .fisher-settings-popup .radio label { color: #374151; font-size: 0.9rem; }
                    .fisher-settings-thresholds {
                        display:flex;
                        gap: 12px;
                        align-items:flex-end;
                        flex-wrap:wrap;
                        margin-top: 10px;
                        padding-top: 10px;
                        border-top: 1px solid #e5e7eb;
                    }
                    """
                ),
                ui.tags.script(
                    f"""
                    (function(){{
                        const btn = document.getElementById('{filter_btn_id}');
                        const popup = document.getElementById('{filter_popup_id}');
                        const runBtn = document.getElementById('{fisher_btn_id}');
                        const runProgress = document.getElementById('{_prefixed_id(prefix, "fisher_run_progress")}');
                        const fisherBtn = document.getElementById('{fisher_settings_btn_id}');
                        const fisherPopup = document.getElementById('{fisher_settings_popup_id}');
                        if (!btn || !popup) return;
                        if (runBtn && runProgress && !runBtn.dataset.bound) {{
                            runBtn.addEventListener('click', () => {{
                                runBtn.disabled = true;
                                runProgress.classList.add('active');
                            }});
                            runBtn.dataset.bound = '1';
                        }}
                        btn.addEventListener('click', () => {{
                            popup.classList.toggle('active');
                            if (fisherPopup) fisherPopup.classList.remove('active');
                        }});
                        if (fisherBtn && fisherPopup) {{
                            fisherBtn.addEventListener('click', () => {{
                                fisherPopup.classList.toggle('active');
                                popup.classList.remove('active');
                            }});
                        }}
                        document.addEventListener('click', (ev) => {{
                            if (!popup.contains(ev.target) && ev.target !== btn && !btn.contains(ev.target)) {{
                                popup.classList.remove('active');
                            }}
                            if (fisherPopup && !fisherPopup.contains(ev.target) && ev.target !== fisherBtn && !(fisherBtn && fisherBtn.contains(ev.target))) {{
                                fisherPopup.classList.remove('active');
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
                {"class": "ks-evidence-toggle"},
                ui.input_checkbox_group(
                    _prefixed_id(prefix, "ks_evidence_types"),
                    None,
                    choices={"in_vivo": "in vivo", "in_vitro": "in vitro"},
                    selected=["in_vivo"],
                    inline=True,
                ),
            )
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
        ks_layout_mode_id = _prefixed_id(prefix, "ks_layout_mode")
        ks_layout_radius_mode_id = _prefixed_id(prefix, "ks_conc_radius_mode")
        ks_two_col_side_id = _prefixed_id(prefix, "ks_two_col_focus_side")
        ks_two_col_columns_id = _prefixed_id(prefix, "ks_two_col_columns")
        ks_dataset_only_id = _prefixed_id(prefix, "ks_dataset_kinases_only")
        controls.append(
            ui.div(
                {"class": "ks-filter-panel"},
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.tags.button(
                        ui.tags.i({"class": "fa fa-sliders"}), "Settings",
                        {"id": _prefixed_id(prefix, "ks_layout_btn"), "type": "button", "class": "ks-filter-button"}
                    ),
                ),
                ui.div(
                    {"id": _prefixed_id(prefix, "ks_layout_popup"), "class": "ks-filter-popup"},
                    ui.div(
                        {"style": "font-weight:600; margin-bottom:6px;"},
                        "Data filters",
                    ),
                    ui.input_checkbox(
                        ks_dataset_only_id,
                        "Only include kinases present in uploaded protein dataset (table + viewer)",
                        value=False,
                    ),
                    ui.tags.hr({"style": "margin:8px 0;"}),
                    ui.div(
                        {"style": "font-weight:600; margin-bottom:6px;"},
                        "Layout",
                    ),
                    ui.input_radio_buttons(
                        ks_layout_mode_id,
                        "Protbox layout",
                        choices={"two_column": "Two-column layout", "concentric": "Concentric layout"},
                        selected="two_column",
                        inline=True,
                    ),
                    ui.panel_conditional(
                        f"input.{ks_layout_mode_id} == 'two_column'",
                        ui.div(
                            {"class": "ks-filter-row"},
                            ui.input_numeric(
                                _prefixed_id(prefix, "ks_layout_gap_x"),
                                "Horizontal gap (px)",
                                value=int(KS_LAYOUT_GAP_X_DEFAULT),
                                min=80,
                                max=1000,
                                width="180px",
                            ),
                            ui.input_numeric(
                                _prefixed_id(prefix, "ks_layout_spacing_y"),
                                "Vertical spacing (px)",
                                value=int(KS_VERTICAL_SPACING),
                                min=10,
                                max=300,
                                width="180px",
                            ),
                        ),
                        ui.div(
                            {"class": "ks-filter-row"},
                            ui.input_radio_buttons(
                                ks_two_col_side_id,
                                "Single-node side",
                                choices={"left": "Left", "right": "Right", "top": "Top", "bottom": "Bottom"},
                                selected="left",
                                inline=True,
                            ),
                            ui.input_numeric(
                                ks_two_col_columns_id,
                                "Columns on multi-node side",
                                value=1,
                                min=1,
                                max=6,
                                width="180px",
                            ),
                        ),
                        ui.div(
                            {"style": "font-size:11px; color:#666; margin-top:2px;"},
                            "Only nodes in the first column nearest the single node receive arrows.",
                        ),
                    ),
                    ui.panel_conditional(
                        f"input.{ks_layout_mode_id} == 'concentric'",
                        ui.div(
                            {"style": "display:flex; flex-direction:column; gap:8px;"},
                            ui.input_checkbox(_prefixed_id(prefix, "ks_conc_arrows"), "Add arrows to center", value=True),
                            ui.input_radio_buttons(
                                ks_layout_radius_mode_id,
                                "Ring radius",
                                choices={"auto": "Auto", "fixed": "Fixed"},
                                selected="auto",
                                inline=True,
                            ),
                            ui.panel_conditional(
                                f"input.{ks_layout_radius_mode_id} == 'fixed'",
                                ui.input_numeric(
                                    _prefixed_id(prefix, "ks_conc_radius_fixed"),
                                    "Fixed radius (px)",
                                    value=int(KS_CONCENTRIC_RADIUS_DEFAULT),
                                    min=20,
                                    max=2000,
                                    width="220px",
                                ),
                            ),
                            ui.panel_conditional(
                                f"input.{ks_layout_radius_mode_id} == 'auto'",
                                ui.input_numeric(
                                    _prefixed_id(prefix, "ks_conc_space"),
                                    "Protbox space (auto)",
                                    value=int(KS_CONCENTRIC_SPACE_DEFAULT),
                                    min=10,
                                    max=400,
                                    width="220px",
                                ),
                            ),
                            ui.input_numeric(
                                _prefixed_id(prefix, "ks_conc_arrow_stop"),
                                "Arrow stop radius (px)",
                                value=int(KS_CONCENTRIC_ARROW_STOP_DEFAULT),
                                min=0,
                                max=500,
                                width="220px",
                            ),
                        ),
                    ),
                ),
                ui.tags.script(
                    f"""
                    (function(){{
                        const btn = document.getElementById('{_prefixed_id(prefix, "ks_layout_btn")}');
                        const popup = document.getElementById('{_prefixed_id(prefix, "ks_layout_popup")}');
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
        controls.append(ui.output_ui(_prefixed_id(prefix, "ks_table")))
    if not is_ks:
        if cfg.get("show_search", False):
            if is_simple_web:
                load_btn_control = ui.input_action_button(
                    _prefixed_id(prefix, "generate"),
                    "Load Pathway",
                    width="auto",
                    style="white-space: nowrap; min-width: 140px;"
                )
                if cfg.get("key") == "web":
                    load_btn = ui.div(
                        {"style": "display:flex; align-items:center; gap:10px;"},
                        load_btn_control,
                        ui.tags.span(
                            {
                                "id": _prefixed_id(prefix, "generate_spinner"),
                                "class": "web-load-spinner",
                                "aria-live": "polite",
                            },
                            "Loading pathway..."
                        ),
                        ui.output_ui(_prefixed_id(prefix, "generate_feedback")),
                        ui.output_ui(_prefixed_id(prefix, "generate_button_state")),
                    )
                else:
                    load_btn = load_btn_control
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
        if cfg.get("key") == "web":
            toggle_controls.append(
                ui.div(
                    {"class": "mk-mode-help-row"},
                    ui.input_checkbox(
                        _prefixed_id(prefix, "simple_kegg_mode"),
                        ui.TagList(
                            "Data Analysis Mode",
                            ui.tags.span(
                                {"class": "mk-inline-help-wrap"},
                                ui.tags.span(
                                    {
                                        "class": "mk-inline-help",
                                        "aria-label": "Data Analysis Mode help",
                                        "tabindex": "0",
                                        "data-help-tooltip-html": (
                                            "<strong>When on:</strong> uses the simpler data-analysis view for KEGG, WikiPathways, and CST. "
                                            "KEGG and WikiPathways show the pathway image background, while CST shows the PDF background and keeps imported CST arrows/text hidden.<br/><br/>"
                                            "<strong>When off:</strong> uses the full layout with arrows and text boxes, "
                                            "while the pathway background starts hidden but can be shown again with the Background Image button."
                                            "<br/><br/>Applies to KEGG, WikiPathways, and CST pathways."
                                        ),
                                    },
                                    "i",
                                ),
                            ),
                        ),
                        value=True,
                    ),
                )
            )
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
    fullscreen_btn = None
    viewer_controls = None
    create_panel = None
    fullscreen_script = None
    toolbar_script = None
    if cfg.get("key") in {"web", "figure", "ks"}:
        btn_id = _prefixed_id(prefix, "fullscreen_btn")
        card_id = _prefixed_id(prefix, "viewer_card")
        panel_id = _prefixed_id(prefix, "viewer_overlay_panel")
        settings_btn_id = _prefixed_id(prefix, "viewer_settings_btn")
        create_panel_id = _prefixed_id(prefix, "viewer_create_panel")
        protein_search_id = _prefixed_id(prefix, "viewer_create_protein_search")
        protein_item_id = _prefixed_id(prefix, "viewer_create_protein_item")
        protein_popover_id = _prefixed_id(prefix, "viewer_create_protein_popover")
        arrow_select_id = _prefixed_id(prefix, "viewer_create_arrow_select")
        shape_select_id = _prefixed_id(prefix, "viewer_create_shape_select")
        text_add_id = _prefixed_id(prefix, "viewer_create_text_add")
        legend_add_id = _prefixed_id(prefix, "viewer_create_legend_add")
        legend_item_id = _prefixed_id(prefix, "viewer_create_legend_item")
        legend_popover_id = _prefixed_id(prefix, "viewer_create_legend_popover")
        undo_btn_id = _prefixed_id(prefix, "viewer_create_undo")
        redo_btn_id = _prefixed_id(prefix, "viewer_create_redo")
        delete_btn_id = _prefixed_id(prefix, "viewer_create_delete")
        interaction_section_id = _prefixed_id(prefix, "viewer_create_interaction_section")
        edge_auto_btn_id = _prefixed_id(prefix, "viewer_create_edge_auto")
        edge_shortest_btn_id = _prefixed_id(prefix, "viewer_create_edge_shortest")
        edge_type_btn_id = _prefixed_id(prefix, "viewer_create_edge_type")
        edge_dash_btn_id = _prefixed_id(prefix, "viewer_create_edge_dash")
        edge_flip_btn_id = _prefixed_id(prefix, "viewer_create_edge_flip")
        drag_mode_btn_id = _prefixed_id(prefix, "viewer_create_drag_mode")
        select_mode_btn_id = _prefixed_id(prefix, "viewer_create_select_mode")
        protbox_snap_btn_id = _prefixed_id(prefix, "viewer_create_protbox_snap")
        protbox_align_btn_id = _prefixed_id(prefix, "viewer_create_protbox_align")
        arrow_snap_btn_id = _prefixed_id(prefix, "viewer_create_arrow_snap")
        control_specs = [
            ("background", "media-image", "Hide or show the pathway background image."),
            ("protboxes", "xmark-square", "Hide or show protein boxes."),
            ("ptms", "xmark-circle", "Hide or show PTM circles and labels."),
            ("ptmLabels", "phoscircle_label", "Hide or show PTM site labels."),
            ("ptmSymbols", "phoscircle_symbol", "Hide or show PTM symbols."),
            ("arrows", "arrow-right", "Hide or show pathway interactions."),
            ("compounds", "select-point3d", "Hide or show compounds."),
            ("groups", "combine", "Hide or show group shapes."),
            ("tooltips", "edit-pencil", "Turn pathway tooltips on or off."),
        ]
        fullscreen_btn = ui.tags.button(
            {
                "id": btn_id,
                "class": "viewer-fullscreen-btn",
                "type": "button",
            },
            "Full Screen",
        )
        viewer_controls = ui.div(
            {"class": "viewer-overlay-panel", "id": panel_id},
            ui.tags.button(
                {
                    "id": settings_btn_id,
                    "class": "viewer-overlay-btn viewer-overlay-settings",
                    "type": "button",
                    "title": "Show or hide viewer display controls.",
                },
                "Hide Objects ^",
            ),
            ui.div(
                {"class": "viewer-overlay-controls"},
                *[
                    ui.tags.button(
                        {
                            "id": _prefixed_id(prefix, f"viewer_toggle_{key}"),
                            "class": "viewer-overlay-btn",
                            "type": "button",
                            "title": title,
                            "aria-label": title,
                        },
                        _icon_markup(label, "viewer-overlay-icon"),
                    )
                    for key, label, title in control_specs
                ],
            ),
        )
        create_panel = ui.div(
            {"class": "viewer-create-panel", "id": create_panel_id},
            ui.div(
                {"class": "viewer-create-scroll"},
                ui.div(
                    {"class": "viewer-create-card"},
                    ui.div(
                        {"class": "viewer-create-group viewer-create-item viewer-create-protein-item", "id": protein_item_id},
                        ui.div({"class": "viewer-create-group-label viewer-create-main-label"}, "Add Objects"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.input({"id": protein_search_id, "class": "viewer-create-field viewer-create-search", "type": "text", "placeholder": "Regex UniProt, gene symbol, HMDB ID, or Wiki ID"}),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group"},
                        ui.div({"class": "viewer-create-group-label"}, "Interaction"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.div(
                                {"class": "viewer-create-select-wrap"},
                                _icon_markup("arrow-right", "viewer-create-select-icon"),
                                ui.tags.select(
                                    {"id": arrow_select_id, "class": "viewer-create-field viewer-create-select viewer-create-select-interaction", "aria-label": "Add Edge"},
                                    ui.tags.option(" ", value="", selected=True, disabled=True),
                                    ui.tags.option("→ Arrow", value="arrow"),
                                    ui.tags.option("⊣ Inhibitor", value="inhibition"),
                                    ui.tags.option("─ Line", value="line"),
                                    ui.tags.option("⇢ Dashed Arrow", value="dashed_arrow"),
                                    ui.tags.option("╌⊣ Dashed Inhibitor", value="dashed_inhibition"),
                                    ui.tags.option("╌ Line", value="dashed_line"),
                                ),
                            ),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group"},
                        ui.div({"class": "viewer-create-group-label"}, "Shape"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.div(
                                {"class": "viewer-create-select-wrap"},
                                _icon_markup("shapes", "viewer-create-select-icon viewer-create-shape-icon"),
                                ui.tags.select(
                                    {"id": shape_select_id, "class": "viewer-create-field viewer-create-select viewer-create-select-shape", "aria-label": "Add Shape"},
                                    ui.tags.option(" ", value="", selected=True, disabled=True),
                                    ui.tags.option("□ Square", value="square"),
                                    ui.tags.option("▢ Rounded", value="rounded"),
                                    ui.tags.option("○ Circle", value="circle"),
                                    ui.tags.option("] Bracket", value="bracket"),
                                ),
                            ),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group"},
                        ui.div({"class": "viewer-create-group-label"}, "Text"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button({"id": text_add_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Add Textbox", "aria-label": "Add Textbox"}, _icon_markup("text-box", "viewer-create-icon")),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group viewer-create-item", "id": legend_item_id},
                        ui.div({"class": "viewer-create-group-label"}, "Legend"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button({"id": legend_add_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Add Legend", "aria-label": "Add Legend"}, _icon_markup("key-plus", "viewer-create-icon")),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-section", "id": interaction_section_id},
                        ui.div({"class": "viewer-create-section-label"}, "Interactions"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button({"id": edge_shortest_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Auto-add arrows using the shortest side-to-side path. Shortcut: Ctrl+Shift+A", "aria-label": "Shortest Auto Add Arrows", "disabled": "disabled"}, _icon_markup("divide-three-solid", "viewer-create-icon")),
                            ui.tags.button({"id": edge_type_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Cycle the selected interactor type. Shortcut: Ctrl+S", "aria-label": "Cycle Interactor Type", "disabled": "disabled"}, _icon_markup("arrow_head_switch", "viewer-create-icon")),
                            ui.tags.button({"id": edge_dash_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Toggle dashed or solid style for the selected interactor. Shortcut: Ctrl+D", "aria-label": "Toggle Interactor Dash", "disabled": "disabled"}, _icon_markup("data-transfer-up", "viewer-create-icon")),
                            ui.tags.button({"id": edge_flip_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Flip the selected interactor direction. Shortcut: Ctrl+F", "aria-label": "Flip Interactor Direction", "disabled": "disabled"}, _icon_markup("ruler-arrows", "viewer-create-icon")),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-section"},
                        ui.div({"class": "viewer-create-section-label"}, "Edit Controls"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button(
                                {"id": undo_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Undo", "aria-label": "Undo"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("undo", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("undo (1)", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button(
                                {"id": redo_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Redo", "aria-label": "Redo"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("redo", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("redo (1)", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button(
                                {"id": delete_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Select an object to delete.", "aria-label": "Delete Selected Object", "disabled": "disabled"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("xmark (1)", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("xmark (2)", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button({"id": drag_mode_btn_id, "class": "viewer-create-action viewer-create-icon-action viewer-create-mode-btn", "type": "button", "title": "Drag Mode", "aria-label": "Drag Mode"}, _icon_markup("drag-hand-gesture", "viewer-create-icon")),
                            ui.tags.button({"id": select_mode_btn_id, "class": "viewer-create-action viewer-create-icon-action viewer-create-mode-btn", "type": "button", "title": "Selection Mode", "aria-label": "Selection Mode"}, _icon_markup("frame-select", "viewer-create-icon")),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-section"},
                        ui.div({"class": "viewer-create-section-label"}, "Snapping"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button(
                                {"id": protbox_snap_btn_id, "class": "viewer-create-action viewer-create-icon-action viewer-create-snap-btn", "type": "button", "aria-label": "Protbox Snapping"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("protboxes_align", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("protboxes_not_align", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button(
                                {"id": protbox_align_btn_id, "class": "viewer-create-action viewer-create-icon-action viewer-create-snap-btn", "type": "button", "aria-label": "Protbox Alignment"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("align_on", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("align_off", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button(
                                {"id": arrow_snap_btn_id, "class": "viewer-create-action viewer-create-icon-action viewer-create-snap-btn", "type": "button", "aria-label": "Arrow Snapping"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("Arrow_snap", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("Arrow_not_snap", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                        ),
                    ),
                ),
            ),
            ui.div(
                {"id": protein_popover_id, "class": "viewer-create-popover viewer-create-protein-popover-overlay"},
            ),
            ui.div(
                {"id": legend_popover_id, "class": "viewer-create-popover viewer-create-legend-popover-overlay"},
                ui.div(
                    {"class": "viewer-create-option-stack"},
                    ui.tags.button({"class": "viewer-create-option", "type": "button", "data-legend-type": "horizontal"}, "Horizontal"),
                    ui.tags.button({"class": "viewer-create-option", "type": "button", "data-legend-type": "vertical"}, "Vertical"),
                ),
            ),
        )
        fullscreen_script = ui.tags.script(
            f"""
            (function(){{
                const btn = document.getElementById('{btn_id}');
                const card = document.getElementById('{card_id}');
                const panel = document.getElementById('{panel_id}');
                const settingsBtn = document.getElementById('{settings_btn_id}');
                const toggleButtons = {{
                    background: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_background")}'),
                    protboxes: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_protboxes")}'),
                    ptms: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_ptms")}'),
                    ptmLabels: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_ptmLabels")}'),
                    ptmSymbols: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_ptmSymbols")}'),
                    arrows: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_arrows")}'),
                    compounds: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_compounds")}'),
                    groups: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_groups")}'),
                    tooltips: document.getElementById('{_prefixed_id(prefix, "viewer_toggle_tooltips")}'),
                }};
                if (!btn || !card || !panel || !settingsBtn || btn.dataset.bound === '1') return;
                const getViewerKey = () => {{
                    const cstShell = card.querySelector('.cst-viewer-shell[data-cst-viewer-key]');
                    if (cstShell) {{
                        return String(cstShell.dataset.cstViewerKey || '');
                    }}
                    return '{prefix}';
                }};
                const getApi = () => {{
                    const viewerKey = getViewerKey();
                    return viewerKey && window.__mkViewerControls ? window.__mkViewerControls[viewerKey] : null;
                }};
                const updateSettingsLabel = () => {{
                    settingsBtn.textContent = panel.classList.contains('is-open') ? 'Hide Objects ⌄' : 'Hide Objects ^';
                }};
                const setPanelOpen = (open) => {{
                    panel.classList.toggle('is-open', !!open);
                    updateSettingsLabel();
                }};
                const syncToggleButtons = (state) => {{
                    if (!state) return;
                    Object.entries(toggleButtons).forEach(([key, el]) => {{
                        if (!el) return;
                        const enabled = state[key] !== false;
                        el.classList.toggle('is-off', !enabled);
                    }});
                }};
                function syncLabel(){{
                    const active = document.fullscreenElement === card;
                    btn.textContent = active ? 'Exit Full Screen' : 'Full Screen';
                    window.requestAnimationFrame(() => window.requestAnimationFrame(() => {{
                        const api = getApi();
                        if (api && typeof api.syncCanvasSize === 'function') api.syncCanvasSize();
                        window.dispatchEvent(new Event('resize'));
                    }}));
                }}
                btn.addEventListener('click', async () => {{
                    try {{
                        if (document.fullscreenElement === card) {{
                            await document.exitFullscreen();
                        }} else {{
                            await card.requestFullscreen();
                        }}
                    }} catch (err) {{
                        console.log('fullscreen toggle failed', err);
                    }}
                    syncLabel();
                }});
                document.addEventListener('fullscreenchange', syncLabel);
                settingsBtn.addEventListener('click', () => {{
                    setPanelOpen(!panel.classList.contains('is-open'));
                }});
                Object.entries(toggleButtons).forEach(([key, el]) => {{
                    if (!el) return;
                    el.addEventListener('click', () => {{
                        const api = getApi();
                        if (!api || typeof api.toggle !== 'function') return;
                        syncToggleButtons(api.toggle(key));
                    }});
                }});
                document.addEventListener('mk-viewer-controls-ready', (evt) => {{
                    if (evt?.detail?.key === getViewerKey()) {{
                        syncToggleButtons(evt.detail.state || {{}});
                    }}
                }});
                btn.dataset.bound = '1';
                setPanelOpen(false);
                syncLabel();
                window.setTimeout(() => {{
                    const api = getApi();
                    if (api && typeof api.getState === 'function') {{
                        syncToggleButtons(api.getState());
                    }}
                }}, 0);
            }})();
            """
        )
        toolbar_script = ui.tags.script(
            f"""
            (function(){{
                const panel = document.getElementById('{create_panel_id}');
                if (!panel || panel.dataset.bound === '1') return;
                const viewerKey = '{prefix}';
                const getApi = () => window.__mkViewerControls && window.__mkViewerControls[viewerKey];
                const createScroll = panel.querySelector('.viewer-create-scroll');
                const proteinSearch = document.getElementById('{protein_search_id}');
                const proteinItem = document.getElementById('{protein_item_id}');
                const proteinPopover = document.getElementById('{protein_popover_id}');
                const arrowSelect = document.getElementById('{arrow_select_id}');
                const arrowSelectWrap = arrowSelect?.closest('.viewer-create-select-wrap') || null;
                const shapeSelect = document.getElementById('{shape_select_id}');
                const textAdd = document.getElementById('{text_add_id}');
                const legendAdd = document.getElementById('{legend_add_id}');
                const legendItem = document.getElementById('{legend_item_id}');
                const legendPopover = document.getElementById('{legend_popover_id}');
                const undoBtn = document.getElementById('{undo_btn_id}');
                const redoBtn = document.getElementById('{redo_btn_id}');
                const deleteBtn = document.getElementById('{delete_btn_id}');
                const interactionSection = document.getElementById('{interaction_section_id}');
                const edgeAutoBtn = document.getElementById('{edge_auto_btn_id}');
                const edgeShortestBtn = document.getElementById('{edge_shortest_btn_id}');
                const edgeTypeBtn = document.getElementById('{edge_type_btn_id}');
                const edgeDashBtn = document.getElementById('{edge_dash_btn_id}');
                const edgeFlipBtn = document.getElementById('{edge_flip_btn_id}');
                const dragModeBtn = document.getElementById('{drag_mode_btn_id}');
                const selectModeBtn = document.getElementById('{select_mode_btn_id}');
                const protboxSnapBtn = document.getElementById('{protbox_snap_btn_id}');
                const protboxAlignBtn = document.getElementById('{protbox_align_btn_id}');
                const arrowSnapBtn = document.getElementById('{arrow_snap_btn_id}');
                const proteinMatches = new Map();
                const proteinSearchDelayMs = 1500;
                let proteinSearchTimer = null;
                let activeArrowPlacementType = '';
                const clearProteinSearchTimer = () => {{
                    if (proteinSearchTimer) {{
                        window.clearTimeout(proteinSearchTimer);
                        proteinSearchTimer = null;
                    }}
                }};
                const closeProteinResults = () => {{
                    clearProteinSearchTimer();
                    proteinMatches.clear();
                    if (proteinPopover) {{
                        proteinPopover.innerHTML = '';
                        proteinPopover.style.display = 'none';
                    }}
                    proteinItem?.classList.remove('is-open');
                }};
                const positionProteinResults = () => {{
                    if (!panel || !proteinSearch || !proteinPopover) return;
                    const panelRect = panel.getBoundingClientRect();
                    const searchRect = proteinSearch.getBoundingClientRect();
                    const left = Math.max(0, searchRect.left - panelRect.left);
                    const top = Math.max(0, searchRect.bottom - panelRect.top + 8);
                    proteinPopover.style.left = `${{left}}px`;
                    proteinPopover.style.top = `${{top}}px`;
                    proteinPopover.style.minWidth = `${{Math.max(240, Math.round(searchRect.width))}}px`;
                }};
                const syncHistoryButtons = (state = {{}}) => {{
                    const canUndo = !!state.canUndo;
                    const canRedo = !!state.canRedo;
                    if (undoBtn) undoBtn.disabled = !canUndo;
                    if (redoBtn) redoBtn.disabled = !canRedo;
                }};
                const syncSelectionButtons = (state = {{}}) => {{
                    const canDelete = !!state.canDelete;
                    const canAutoConnect = !!state.canAutoConnect;
                    const canModifyInteractor = !!state.canModifyInteractor;
                    if (deleteBtn) {{
                        deleteBtn.disabled = !canDelete;
                        deleteBtn.title = canDelete
                            ? 'Delete selected object.'
                            : 'Select an object to delete.';
                    }}
                    interactionSection?.classList.toggle('is-disabled', !canAutoConnect && !canModifyInteractor);
                    if (edgeAutoBtn) edgeAutoBtn.disabled = !canAutoConnect;
                    if (edgeShortestBtn) edgeShortestBtn.disabled = !canAutoConnect;
                    if (edgeTypeBtn) edgeTypeBtn.disabled = !canModifyInteractor;
                    if (edgeDashBtn) edgeDashBtn.disabled = !canModifyInteractor;
                    if (edgeFlipBtn) edgeFlipBtn.disabled = !canModifyInteractor;
                }};
                const syncMouseModeButtons = (mode = 'drag') => {{
                    dragModeBtn?.classList.toggle('is-active', mode === 'drag');
                    selectModeBtn?.classList.toggle('is-active', mode === 'selection');
                }};
                const syncArrowPlacementState = (state = {{}}) => {{
                    activeArrowPlacementType = state && state.active ? String(state.type || 'arrow') : '';
                    if (arrowSelect) {{
                        arrowSelect.classList.toggle('is-active', !!activeArrowPlacementType);
                        arrowSelect.title = activeArrowPlacementType
                            ? 'Arrow placement is active. Click again to cancel.'
                            : 'Add Arrow';
                    }}
                    arrowSelectWrap?.classList.toggle('is-active', !!activeArrowPlacementType);
                }};
                const syncSnapButtons = (state = {{}}) => {{
                    const protboxSnap = state.protboxSnap !== false;
                    const protboxAlign = state.protboxAlign !== false;
                    const arrowSnap = state.arrowSnap !== false;
                    protboxSnapBtn?.classList.toggle('is-active', protboxSnap);
                    protboxAlignBtn?.classList.toggle('is-active', protboxAlign);
                    arrowSnapBtn?.classList.toggle('is-active', arrowSnap);
                    if (protboxSnapBtn) {{
                        protboxSnapBtn.title = protboxSnap
                            ? 'Protbox snapping is on. Nearby protboxes snap together; hold Shift while dragging to temporarily disable it.'
                            : 'Protbox snapping is off. Nearby protboxes will not snap together; hold Shift while dragging to temporarily enable it.';
                    }}
                    if (protboxAlignBtn) {{
                        protboxAlignBtn.title = protboxAlign
                            ? 'Protbox alignment is on. Protboxes align to matching horizontal or vertical positions; hold Shift while dragging to temporarily disable it.'
                            : 'Protbox alignment is off. Protboxes will not align to shared horizontal or vertical positions; hold Shift while dragging to temporarily enable it.';
                    }}
                    if (arrowSnapBtn) {{
                        arrowSnapBtn.title = arrowSnap
                            ? 'Arrow snapping is on. Interactor ends snap to protbox sides and snap guides are shown.'
                            : 'Arrow snapping is off. Interactor ends do not snap to protbox sides and snap guides are hidden.';
                    }}
                }};
                const positionLegendPopover = () => {{
                    if (!panel || !legendAdd || !legendItem || !legendPopover) return;
                    const panelRect = panel.getBoundingClientRect();
                    const btnRect = legendAdd.getBoundingClientRect();
                    const left = Math.max(0, btnRect.left - panelRect.left);
                    const top = Math.max(0, btnRect.bottom - panelRect.top + 8);
                    legendPopover.style.left = `${{left}}px`;
                    legendPopover.style.top = `${{top}}px`;
                }};
                const closeLegendPopover = () => {{
                    if (legendPopover) {{
                        legendPopover.style.display = 'none';
                    }}
                    legendItem?.classList.remove('is-open');
                }};
                const proteinColorToCss = (value) => {{
                    if (Array.isArray(value) && value.length >= 3) {{
                        const parts = value.slice(0, 3).map((part) => {{
                            const num = Number(part);
                            if (!Number.isFinite(num)) return 0;
                            return Math.max(0, Math.min(255, Math.round(num)));
                        }});
                        return `rgb(${{parts[0]}}, ${{parts[1]}}, ${{parts[2]}})`;
                    }}
                    if (typeof value === 'string' && value.trim()) {{
                        return value.trim();
                    }}
                    return '#d1d5db';
                }};
                const renderProteinResults = () => {{
                    if (!proteinSearch || !proteinPopover) return;
                    proteinMatches.clear();
                    proteinPopover.innerHTML = '';
                    proteinItem?.classList.remove('is-open');
                    const api = getApi();
                    if (!api || (typeof api.searchObjects !== 'function' && typeof api.searchProteins !== 'function')) {{
                        return;
                    }}
                    const query = proteinSearch.value.trim();
                    if (!query) {{
                        return;
                    }}
                    const payload = (typeof api.searchObjects === 'function'
                        ? api.searchObjects(query, 40)
                        : api.searchProteins(query, 40)) || {{}};
                    if (payload.error) {{
                        return;
                    }}
                    const results = Array.isArray(payload.results) ? payload.results : [];
                    if (!results.length) {{
                        return;
                    }}
                    const stack = document.createElement('div');
                    stack.className = 'viewer-create-option-stack';
                    results.forEach((entry) => {{
                        const isMetabolite = entry.kind === 'metabolite';
                        const label = isMetabolite
                            ? (entry.label || `${{entry.hmdbId || ''}} - ${{entry.wikipediaId || entry.displayLabel || ''}}`.replace(/^\\s*-\\s*|\\s*-\\s*$/g, ''))
                            : `${{entry.uniprot}} - ${{entry.geneSymbol}}`;
                        proteinMatches.set(label, {{
                            kind: isMetabolite ? 'metabolite' : 'protein',
                            value: isMetabolite ? entry.hmdbId : (entry.uniprot || entry.geneSymbol || '')
                        }});
                        const option = document.createElement('button');
                        option.className = 'viewer-create-option viewer-create-protein-option';
                        option.type = 'button';
                        const text = document.createElement('span');
                        text.className = 'viewer-create-protein-label';
                        text.textContent = label;
                        const swatch = document.createElement('span');
                        swatch.className = isMetabolite
                            ? 'viewer-create-metabolite-swatch'
                            : 'viewer-create-protein-swatch';
                        swatch.style.background = proteinColorToCss(entry.color);
                        option.appendChild(text);
                        option.appendChild(swatch);
                        option.addEventListener('click', () => {{
                            const api = getApi();
                            if (isMetabolite) {{
                                if (api && typeof api.addCompound === 'function') {{
                                    api.addCompound(entry.hmdbId);
                                }}
                            }} else if (api && typeof api.addProtbox === 'function') {{
                                api.addProtbox(entry.uniprot || entry.geneSymbol || '');
                            }}
                            if (proteinSearch) {{
                                proteinSearch.value = '';
                            }}
                            closeProteinResults();
                        }});
                        stack.appendChild(option);
                    }});
                    proteinPopover.appendChild(stack);
                    positionProteinResults();
                    proteinPopover.style.display = 'block';
                    proteinItem?.classList.add('is-open');
                }};
                const scheduleProteinResults = () => {{
                    if (!proteinSearch) return;
                    const query = proteinSearch.value.trim();
                    if (!query) {{
                        closeProteinResults();
                        return;
                    }}
                    clearProteinSearchTimer();
                    proteinSearchTimer = window.setTimeout(() => {{
                        proteinSearchTimer = null;
                        renderProteinResults();
                    }}, proteinSearchDelayMs);
                }};
                proteinSearch?.addEventListener('input', scheduleProteinResults);
                proteinSearch?.addEventListener('keydown', (evt) => {{
                    if (evt.key === 'Escape') {{
                        evt.preventDefault();
                        closeProteinResults();
                        return;
                    }}
                    if (evt.key === 'Enter') {{
                        evt.preventDefault();
                        const raw = proteinSearch?.value?.trim() || '';
                        const match = proteinMatches.get(raw) || null;
                        const api = getApi();
                        if (api && match && match.value) {{
                            if (match.kind === 'metabolite' && typeof api.addCompound === 'function') {{
                                api.addCompound(match.value);
                            }} else if (match.kind === 'protein' && typeof api.addProtbox === 'function') {{
                                api.addProtbox(match.value);
                            }} else {{
                                renderProteinResults();
                                return;
                            }}
                            proteinSearch.value = '';
                            closeProteinResults();
                        }} else if (raw) {{
                            renderProteinResults();
                        }}
                    }}
                }});
                proteinSearch?.addEventListener('focus', () => {{
                    if ((proteinSearch?.value || '').trim()) {{
                        scheduleProteinResults();
                    }}
                }});
                createScroll?.addEventListener('scroll', () => {{
                    if (proteinItem?.classList.contains('is-open')) {{
                        positionProteinResults();
                    }}
                    if (legendItem?.classList.contains('is-open')) {{
                        positionLegendPopover();
                    }}
                }});
                window.addEventListener('resize', () => {{
                    if (proteinItem?.classList.contains('is-open')) {{
                        positionProteinResults();
                    }}
                    if (legendItem?.classList.contains('is-open')) {{
                        positionLegendPopover();
                    }}
                }});
                arrowSelect?.addEventListener('change', () => {{
                    const value = arrowSelect?.value || '';
                    if (!value) return;
                    const api = getApi();
                    if (api && typeof api.addArrow === 'function') {{
                        api.addArrow(value);
                    }}
                    arrowSelect.value = '';
                }});
                arrowSelect?.addEventListener('mousedown', (evt) => {{
                    if (!activeArrowPlacementType) return;
                    evt.preventDefault();
                    evt.stopPropagation();
                    const api = getApi();
                    if (api && typeof api.addArrow === 'function') {{
                        api.addArrow(activeArrowPlacementType);
                    }}
                }});
                shapeSelect?.addEventListener('change', () => {{
                    const value = shapeSelect?.value || '';
                    if (!value) return;
                    const api = getApi();
                    if (api && typeof api.addShape === 'function') {{
                        api.addShape(value);
                    }}
                    shapeSelect.value = '';
                }});
                undoBtn?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.undo === 'function') {{
                        syncHistoryButtons(api.undo() || {{}});
                    }}
                }});
                redoBtn?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.redo === 'function') {{
                        syncHistoryButtons(api.redo() || {{}});
                    }}
                }});
                edgeAutoBtn?.addEventListener('click', () => {{
                    if (edgeAutoBtn?.disabled) return;
                    const api = getApi();
                    if (api && typeof api.autoConnectEdges === 'function') {{
                        api.autoConnectEdges();
                    }}
                }});
                edgeShortestBtn?.addEventListener('click', () => {{
                    if (edgeShortestBtn?.disabled) return;
                    const api = getApi();
                    if (api && typeof api.autoConnectShortestEdges === 'function') {{
                        api.autoConnectShortestEdges();
                    }}
                }});
                edgeTypeBtn?.addEventListener('click', () => {{
                    if (edgeTypeBtn?.disabled) return;
                    const api = getApi();
                    if (api && typeof api.cycleSelectedEdgeType === 'function') {{
                        api.cycleSelectedEdgeType();
                    }}
                }});
                edgeDashBtn?.addEventListener('click', () => {{
                    if (edgeDashBtn?.disabled) return;
                    const api = getApi();
                    if (api && typeof api.toggleSelectedEdgeDash === 'function') {{
                        api.toggleSelectedEdgeDash();
                    }}
                }});
                edgeFlipBtn?.addEventListener('click', () => {{
                    if (edgeFlipBtn?.disabled) return;
                    const api = getApi();
                    if (api && typeof api.flipSelectedEdgeDirection === 'function') {{
                        api.flipSelectedEdgeDirection();
                    }}
                }});
                deleteBtn?.addEventListener('click', () => {{
                    if (deleteBtn?.disabled) return;
                    const api = getApi();
                    if (api && typeof api.deleteSelected === 'function') {{
                        syncSelectionButtons(api.deleteSelected() || {{}});
                    }}
                }});
                dragModeBtn?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.setMouseMode === 'function') {{
                        syncMouseModeButtons(api.setMouseMode('drag') || 'drag');
                    }} else {{
                        syncMouseModeButtons('drag');
                    }}
                }});
                selectModeBtn?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.setMouseMode === 'function') {{
                        syncMouseModeButtons(api.setMouseMode('selection') || 'selection');
                    }} else {{
                        syncMouseModeButtons('selection');
                    }}
                }});
                protboxSnapBtn?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.setProtboxSnapEnabled === 'function' && typeof api.getProtboxSnapEnabled === 'function') {{
                        syncSnapButtons(api.setProtboxSnapEnabled(!api.getProtboxSnapEnabled()) || {{}});
                    }}
                }});
                protboxAlignBtn?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.setProtboxAlignEnabled === 'function' && typeof api.getProtboxAlignEnabled === 'function') {{
                        syncSnapButtons(api.setProtboxAlignEnabled(!api.getProtboxAlignEnabled()) || {{}});
                    }}
                }});
                arrowSnapBtn?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.setArrowSnapEnabled === 'function' && typeof api.getArrowSnapEnabled === 'function') {{
                        syncSnapButtons(api.setArrowSnapEnabled(!api.getArrowSnapEnabled()) || {{}});
                    }}
                }});
                textAdd?.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.addText === 'function') {{
                        api.addText();
                    }}
                }});
                legendAdd?.addEventListener('click', () => {{
                    const opening = !legendItem?.classList.contains('is-open');
                    if (!opening) {{
                        closeLegendPopover();
                        return;
                    }}
                    positionLegendPopover();
                    if (legendPopover) {{
                        legendPopover.style.display = 'block';
                    }}
                    legendItem?.classList.add('is-open');
                }});
                proteinItem?.addEventListener('click', (evt) => evt.stopPropagation());
                legendItem?.addEventListener('click', (evt) => evt.stopPropagation());
                panel.querySelectorAll('[data-legend-type]').forEach((el) => {{
                    el.addEventListener('click', () => {{
                        const api = getApi();
                        if (api && typeof api.addLegend === 'function') {{
                            api.addLegend(el.getAttribute('data-legend-type') || 'vertical');
                        }}
                        closeLegendPopover();
                    }});
                }});
                document.addEventListener('click', () => {{
                    closeProteinResults();
                    closeLegendPopover();
                }});
                document.addEventListener('mk-viewer-history-state', (evt) => {{
                    if (evt?.detail?.key === viewerKey) {{
                        syncHistoryButtons(evt.detail || {{}});
                    }}
                }});
                document.addEventListener('mk-viewer-selection-state', (evt) => {{
                    if (evt?.detail?.key === viewerKey) {{
                        syncSelectionButtons(evt.detail || {{}});
                    }}
                }});
                document.addEventListener('mk-viewer-mouse-mode', (evt) => {{
                    if (evt?.detail?.key === viewerKey) {{
                        syncMouseModeButtons(evt.detail.mode || 'drag');
                    }}
                }});
                document.addEventListener('mk-viewer-arrow-placement-state', (evt) => {{
                    if (evt?.detail?.key === viewerKey) {{
                        syncArrowPlacementState(evt.detail || {{}});
                    }}
                }});
                document.addEventListener('mk-viewer-snap-state', (evt) => {{
                    if (evt?.detail?.key === viewerKey) {{
                        syncSnapButtons(evt.detail || {{}});
                    }}
                }});
                const api = getApi();
                if (api && typeof api.getHistoryState === 'function') {{
                    syncHistoryButtons(api.getHistoryState() || {{}});
                }} else {{
                    syncHistoryButtons({{}});
                }}
                if (api && typeof api.getSelectionState === 'function') {{
                    syncSelectionButtons(api.getSelectionState() || {{}});
                }} else {{
                    syncSelectionButtons({{}});
                }}
                if (api && typeof api.getMouseMode === 'function') {{
                    syncMouseModeButtons(api.getMouseMode() || 'drag');
                }} else {{
                    syncMouseModeButtons('drag');
                }}
                if (api && typeof api.getArrowPlacementMode === 'function') {{
                    const arrowPlacementType = String(api.getArrowPlacementMode() || '');
                    syncArrowPlacementState({{ active: !!arrowPlacementType, type: arrowPlacementType }});
                }} else {{
                    syncArrowPlacementState({{}});
                }}
                if (api && typeof api.getSnapState === 'function') {{
                    syncSnapButtons(api.getSnapState() || {{}});
                }} else {{
                    syncSnapButtons({{}});
                }}
                panel.dataset.bound = '1';
            }})();
            """
        )
    download_controls = ui.div(
        {"class": "svg-download-row"},
        ui.span({"class": "svg-download-label"}, "Download:"),
        ui.tags.button(
            {
                "class": "btn btn-primary btn-sm",
                "data-mk-download": "svg",
                "data-mk-prefix": prefix,
                "data-mk-name": download_name,
                "type": "button",
            },
            "SVG",
        ),
        ui.tags.button(
            {
                "class": "btn btn-outline-primary btn-sm",
                "data-mk-download": "svgprint",
                "data-mk-prefix": prefix,
                "data-mk-name": download_name,
                "type": "button",
                "title": "Print-ready SVG with outlined label paths for consistent publisher rendering",
            },
            "SVG Print",
        ),
        ui.tags.button(
            {
                "class": "btn btn-outline-primary btn-sm",
                "data-mk-download": "pdf",
                "data-mk-prefix": prefix,
                "data-mk-name": download_name,
                "type": "button",
            },
            "PDF",
        ),
        ui.tags.button(
            {
                "class": "btn btn-outline-primary btn-sm",
                "data-mk-download": "png",
                "data-mk-prefix": prefix,
                "data-mk-name": download_name,
                "type": "button",
            },
            "PNG",
        ),
        ui.tags.button(
            {
                "class": "btn btn-outline-primary btn-sm",
                "data-mk-download": "jpeg",
                "data-mk-prefix": prefix,
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
                "id": _prefixed_id(prefix, "viewer_card"),
                "data-prefix": prefix,
                "data-download-name": download_name,
            },
            *( [create_panel] if create_panel is not None else [] ),
            ui.div(
                {"class": "pathway-viewer-body"},
                *( [fullscreen_btn] if fullscreen_btn is not None else [] ),
                *( [viewer_controls] if viewer_controls is not None else [] ),
                ui.output_ui(_prefixed_id(prefix, "pathway_preview")),
            ),
            download_controls,
            *( [fullscreen_script] if fullscreen_script is not None else [] ),
            *( [toolbar_script] if toolbar_script is not None else [] ),
        ),
        ui.hr(),
        ui.output_text(_prefixed_id(prefix, "status_message")),
        ui.output_text(_prefixed_id(prefix, "json_summary")),
    )


def _cst_create_panel(prefix: str) -> ui.TagList:
    card_id = _prefixed_id(prefix, "viewer_card")
    panel_id = _prefixed_id(prefix, "cst_create_panel")
    protein_item_id = _prefixed_id(prefix, "cst_create_protein_item")
    protein_search_id = _prefixed_id(prefix, "cst_create_protein_search")
    protein_popover_id = _prefixed_id(prefix, "cst_create_protein_popover")
    legend_item_id = _prefixed_id(prefix, "cst_create_legend_item")
    legend_add_id = _prefixed_id(prefix, "cst_create_legend_add")
    legend_popover_id = _prefixed_id(prefix, "cst_create_legend_popover")
    arrow_select_id = _prefixed_id(prefix, "cst_create_arrow_select")
    shape_select_id = _prefixed_id(prefix, "cst_create_shape_select")
    text_add_id = _prefixed_id(prefix, "cst_create_text_add")
    undo_btn_id = _prefixed_id(prefix, "cst_create_undo")
    redo_btn_id = _prefixed_id(prefix, "cst_create_redo")
    delete_btn_id = _prefixed_id(prefix, "cst_create_delete")
    edge_shortest_btn_id = _prefixed_id(prefix, "cst_create_edge_shortest")
    edge_type_btn_id = _prefixed_id(prefix, "cst_create_edge_type")
    edge_dash_btn_id = _prefixed_id(prefix, "cst_create_edge_dash")
    edge_flip_btn_id = _prefixed_id(prefix, "cst_create_edge_flip")
    select_mode_btn_id = _prefixed_id(prefix, "cst_create_select_mode")
    front_btn_id = _prefixed_id(prefix, "cst_create_front")
    back_btn_id = _prefixed_id(prefix, "cst_create_back")
    return ui.TagList(
        ui.div(
            {"class": "viewer-create-panel", "id": panel_id},
            ui.div(
                {"class": "viewer-create-scroll"},
                ui.div(
                    {"class": "viewer-create-card"},
                    ui.div(
                        {"class": "viewer-create-group viewer-create-item", "id": protein_item_id},
                        ui.div({"class": "viewer-create-group-label viewer-create-main-label"}, "Add Objects"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.input({"id": protein_search_id, "class": "viewer-create-field viewer-create-search", "type": "text", "placeholder": "Regex UniProt, gene symbol, HMDB ID, or Wiki ID"}),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group"},
                        ui.div({"class": "viewer-create-group-label"}, "Interaction"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.div(
                                {"class": "viewer-create-select-wrap"},
                                _icon_markup("arrow-right", "viewer-create-select-icon"),
                                ui.tags.select(
                                    {"id": arrow_select_id, "class": "viewer-create-field viewer-create-select viewer-create-select-interaction", "aria-label": "Add Edge"},
                                    ui.tags.option(" ", value="", selected=True, disabled=True),
                                    ui.tags.option("→ Arrow", value="arrow"),
                                    ui.tags.option("⊣ Inhibitor", value="inhibition"),
                                    ui.tags.option("─ Line", value="line"),
                                    ui.tags.option("⇢ Dashed Arrow", value="dashed_arrow"),
                                    ui.tags.option("╌⊣ Dashed Inhibitor", value="dashed_inhibition"),
                                    ui.tags.option("╌ Line", value="dashed_line"),
                                ),
                            ),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group"},
                        ui.div({"class": "viewer-create-group-label"}, "Shape"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.div(
                                {"class": "viewer-create-select-wrap"},
                                _icon_markup("shapes", "viewer-create-select-icon viewer-create-shape-icon"),
                                ui.tags.select(
                                    {"id": shape_select_id, "class": "viewer-create-field viewer-create-select viewer-create-select-shape", "aria-label": "Add Shape"},
                                    ui.tags.option(" ", value="", selected=True, disabled=True),
                                    ui.tags.option("□ Square", value="square"),
                                    ui.tags.option("▢ Rounded", value="rounded"),
                                    ui.tags.option("○ Circle", value="circle"),
                                    ui.tags.option("] Bracket", value="bracket"),
                                ),
                            ),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group"},
                        ui.div({"class": "viewer-create-group-label"}, "Text"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button({"id": text_add_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Add Textbox", "aria-label": "Add Textbox"}, _icon_markup("text-box", "viewer-create-icon")),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-group viewer-create-item", "id": legend_item_id},
                        ui.div({"class": "viewer-create-group-label"}, "Legend"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button({"id": legend_add_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Add Legend", "aria-label": "Add Legend"}, _icon_markup("key-plus", "viewer-create-icon")),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-section"},
                        ui.div({"class": "viewer-create-section-label"}, "Interactions"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button({"id": edge_shortest_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Auto-add arrows using the shortest side-to-side path.", "aria-label": "Shortest Auto Add Arrows", "disabled": "disabled"}, _icon_markup("divide-three-solid", "viewer-create-icon")),
                            ui.tags.button({"id": edge_type_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Cycle the selected interactor type. Shortcut: Ctrl+S", "aria-label": "Cycle Interactor Type", "disabled": "disabled"}, _icon_markup("arrow_head_switch", "viewer-create-icon")),
                            ui.tags.button({"id": edge_dash_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Toggle dashed or solid style for the selected interactor. Shortcut: Ctrl+D", "aria-label": "Toggle Interactor Dash", "disabled": "disabled"}, _icon_markup("data-transfer-up", "viewer-create-icon")),
                            ui.tags.button({"id": edge_flip_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Flip the selected interactor direction.", "aria-label": "Flip Interactor Direction", "disabled": "disabled"}, _icon_markup("ruler-arrows", "viewer-create-icon")),
                        ),
                    ),
                    ui.div(
                        {"class": "viewer-create-section"},
                        ui.div({"class": "viewer-create-section-label"}, "Edit Controls"),
                        ui.div(
                            {"class": "viewer-create-inline"},
                            ui.tags.button(
                                {"id": undo_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Undo", "aria-label": "Undo"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("undo", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("undo (1)", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button(
                                {"id": redo_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Redo", "aria-label": "Redo"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("redo", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("redo (1)", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button(
                                {"id": delete_btn_id, "class": "viewer-create-action viewer-create-icon-action", "type": "button", "title": "Delete selected object", "aria-label": "Delete Selected Object", "disabled": "disabled"},
                                ui.span({"class": "viewer-create-icon-pair"}, _icon_markup("xmark (1)", "viewer-create-icon viewer-create-icon-primary"), _icon_markup("xmark (2)", "viewer-create-icon viewer-create-icon-secondary")),
                            ),
                            ui.tags.button({"id": select_mode_btn_id, "class": "viewer-create-action viewer-create-icon-action viewer-create-mode-btn", "type": "button", "title": "Selection Mode", "aria-label": "Selection Mode"}, _icon_markup("frame-select", "viewer-create-icon")),
                            ui.tags.button({"id": back_btn_id, "class": "viewer-create-action", "type": "button", "title": "Send selected module behind the others", "disabled": "disabled"}, "Back"),
                            ui.tags.button({"id": front_btn_id, "class": "viewer-create-action", "type": "button", "title": "Bring selected module in front of the others", "disabled": "disabled"}, "Front"),
                        ),
                    ),
                ),
            ),
            ui.div(
                {"id": protein_popover_id, "class": "viewer-create-popover viewer-create-protein-popover-overlay"},
            ),
            ui.div(
                {"id": legend_popover_id, "class": "viewer-create-popover viewer-create-legend-popover-overlay"},
                ui.div(
                    {"class": "viewer-create-option-stack"},
                    ui.tags.button({"class": "viewer-create-option", "type": "button", "data-legend-type": "horizontal"}, "Horizontal"),
                    ui.tags.button({"class": "viewer-create-option", "type": "button", "data-legend-type": "vertical"}, "Vertical"),
                ),
            ),
        ),
        ui.tags.script(
            f"""
            (function(){{
                const panel = document.getElementById('{panel_id}');
                const card = document.getElementById('{card_id}');
                const createScroll = panel?.querySelector('.viewer-create-scroll');
                const proteinSearch = document.getElementById('{protein_search_id}');
                const proteinItem = document.getElementById('{protein_item_id}') || panel?.querySelector('.viewer-create-item');
                const proteinPopover = document.getElementById('{protein_popover_id}');
                const legendAdd = document.getElementById('{legend_add_id}');
                const legendItem = document.getElementById('{legend_item_id}');
                const legendPopover = document.getElementById('{legend_popover_id}');
                const arrowSelect = document.getElementById('{arrow_select_id}');
                const shapeSelect = document.getElementById('{shape_select_id}');
                const textAdd = document.getElementById('{text_add_id}');
                const undoBtn = document.getElementById('{undo_btn_id}');
                const redoBtn = document.getElementById('{redo_btn_id}');
                const deleteBtn = document.getElementById('{delete_btn_id}');
                const edgeShortestBtn = document.getElementById('{edge_shortest_btn_id}');
                const edgeTypeBtn = document.getElementById('{edge_type_btn_id}');
                const edgeDashBtn = document.getElementById('{edge_dash_btn_id}');
                const edgeFlipBtn = document.getElementById('{edge_flip_btn_id}');
                const selectModeBtn = document.getElementById('{select_mode_btn_id}');
                const frontBtn = document.getElementById('{front_btn_id}');
                const backBtn = document.getElementById('{back_btn_id}');
                if (!panel || !card || panel.dataset.bound === '1') return;
                const proteinMatches = new Map();
                const proteinSearchDelayMs = 1500;
                let proteinSearchTimer = null;
                const getViewerKey = () => {{
                    const shell = card.querySelector('.cst-viewer-shell[data-cst-viewer-key]');
                    return shell ? String(shell.dataset.cstViewerKey || '') : '';
                }};
                const getApi = () => {{
                    const viewerKey = getViewerKey();
                    return viewerKey && window.__mkCstViewerControls ? window.__mkCstViewerControls[viewerKey] : null;
                }};
                const setEnabled = (el, enabled) => {{
                    if (!el) return;
                    el.disabled = !enabled;
                    el.classList.toggle('is-disabled', !enabled);
                }};
                const clearProteinSearchTimer = () => {{
                    if (proteinSearchTimer) {{
                        window.clearTimeout(proteinSearchTimer);
                        proteinSearchTimer = null;
                    }}
                }};
                const closeProteinResults = () => {{
                    clearProteinSearchTimer();
                    proteinMatches.clear();
                    if (proteinPopover) {{
                        proteinPopover.innerHTML = '';
                        proteinPopover.style.display = 'none';
                    }}
                    proteinItem?.classList.remove('is-open');
                }};
                const positionProteinResults = () => {{
                    if (!panel || !proteinSearch || !proteinPopover) return;
                    const panelRect = panel.getBoundingClientRect();
                    const searchRect = proteinSearch.getBoundingClientRect();
                    const left = Math.max(0, searchRect.left - panelRect.left);
                    const top = Math.max(0, searchRect.bottom - panelRect.top + 8);
                    proteinPopover.style.left = `${{left}}px`;
                    proteinPopover.style.top = `${{top}}px`;
                    proteinPopover.style.minWidth = `${{Math.max(240, Math.round(searchRect.width))}}px`;
                }};
                const positionLegendPopover = () => {{
                    if (!panel || !legendAdd || !legendItem || !legendPopover) return;
                    const panelRect = panel.getBoundingClientRect();
                    const btnRect = legendAdd.getBoundingClientRect();
                    const left = Math.max(0, btnRect.left - panelRect.left);
                    const top = Math.max(0, btnRect.bottom - panelRect.top + 8);
                    legendPopover.style.left = `${{left}}px`;
                    legendPopover.style.top = `${{top}}px`;
                }};
                const closeLegendPopover = () => {{
                    if (legendPopover) {{
                        legendPopover.style.display = 'none';
                    }}
                    legendItem?.classList.remove('is-open');
                }};
                const proteinColorToCss = (value) => {{
                    if (Array.isArray(value) && value.length >= 3) {{
                        const parts = value.slice(0, 3).map((part) => {{
                            const num = Number(part);
                            if (!Number.isFinite(num)) return 0;
                            return Math.max(0, Math.min(255, Math.round(num)));
                        }});
                        return `rgb(${{parts[0]}}, ${{parts[1]}}, ${{parts[2]}})`;
                    }}
                    if (typeof value === 'string' && value.trim()) {{
                        return value.trim();
                    }}
                    return '#d1d5db';
                }};
                const renderProteinResults = () => {{
                    if (!proteinSearch || !proteinPopover) return;
                    proteinMatches.clear();
                    proteinPopover.innerHTML = '';
                    proteinItem?.classList.remove('is-open');
                    const api = getApi();
                    if (!api || (typeof api.searchObjects !== 'function' && typeof api.searchProteins !== 'function')) return;
                    const query = proteinSearch.value.trim();
                    if (!query) return;
                    const payload = (typeof api.searchObjects === 'function'
                        ? api.searchObjects(query, 40)
                        : api.searchProteins(query, 40)) || {{}};
                    if (payload.error) return;
                    const results = Array.isArray(payload.results) ? payload.results : [];
                    if (!results.length) return;
                    const stack = document.createElement('div');
                    stack.className = 'viewer-create-option-stack';
                    results.forEach((entry) => {{
                        const isMetabolite = entry.kind === 'metabolite';
                        const label = isMetabolite
                            ? (entry.label || `${{entry.hmdbId || ''}} - ${{entry.wikipediaId || entry.displayLabel || ''}}`.replace(/^\\s*-\\s*|\\s*-\\s*$/g, ''))
                            : `${{entry.uniprot}} - ${{entry.geneSymbol}}`;
                        proteinMatches.set(label, {{
                            kind: isMetabolite ? 'metabolite' : 'protein',
                            value: isMetabolite ? entry.hmdbId : entry.uniprot
                        }});
                        const option = document.createElement('button');
                        option.className = 'viewer-create-option viewer-create-protein-option';
                        option.type = 'button';
                        const text = document.createElement('span');
                        text.className = 'viewer-create-protein-label';
                        text.textContent = label;
                        const swatch = document.createElement('span');
                        swatch.className = isMetabolite
                            ? 'viewer-create-metabolite-swatch'
                            : 'viewer-create-protein-swatch';
                        swatch.style.background = proteinColorToCss(entry.color || entry.fc_color_1 || entry['fc_color_1']);
                        option.appendChild(text);
                        option.appendChild(swatch);
                        option.addEventListener('click', () => {{
                            const api = getApi();
                            if (isMetabolite) {{
                                if (api && typeof api.addCompound === 'function') {{
                                    api.addCompound(entry.hmdbId);
                                }}
                            }} else if (api && typeof api.addProtbox === 'function') {{
                                api.addProtbox(entry.uniprot || entry.geneSymbol || '');
                            }}
                            if (proteinSearch) proteinSearch.value = '';
                            closeProteinResults();
                            syncState();
                        }});
                        stack.appendChild(option);
                    }});
                    proteinPopover.appendChild(stack);
                    positionProteinResults();
                    proteinPopover.style.display = 'block';
                    proteinItem?.classList.add('is-open');
                }};
                const scheduleProteinResults = () => {{
                    if (!proteinSearch) return;
                    const query = proteinSearch.value.trim();
                    if (!query) {{
                        closeProteinResults();
                        return;
                    }}
                    clearProteinSearchTimer();
                    proteinSearchTimer = window.setTimeout(() => {{
                        proteinSearchTimer = null;
                        renderProteinResults();
                    }}, proteinSearchDelayMs);
                }};
                const syncState = () => {{
                    const api = getApi();
                    const state = api && typeof api.getState === 'function' ? (api.getState() || {{}}) : {{}};
                    setEnabled(undoBtn, !!state.canUndo);
                    setEnabled(redoBtn, !!state.canRedo);
                    setEnabled(deleteBtn, !!state.canDelete);
                    setEnabled(edgeShortestBtn, !!state.canAutoConnect);
                    setEnabled(edgeTypeBtn, !!state.canEdgeEdit);
                    setEnabled(edgeDashBtn, !!state.canEdgeEdit);
                    setEnabled(edgeFlipBtn, !!state.canEdgeEdit);
                    setEnabled(frontBtn, !!state.canArrange);
                    setEnabled(backBtn, !!state.canArrange);
                    if (selectModeBtn) selectModeBtn.classList.toggle('is-active', String(state.mouseMode || 'drag') === 'selection');
                }};
                proteinSearch?.addEventListener('input', scheduleProteinResults);
                proteinSearch?.addEventListener('keydown', (evt) => {{
                    if (evt.key === 'Escape') {{
                        evt.preventDefault();
                        closeProteinResults();
                        return;
                    }}
                    if (evt.key === 'Enter') {{
                        evt.preventDefault();
                        const raw = proteinSearch?.value?.trim() || '';
                        const match = proteinMatches.get(raw) || null;
                        const api = getApi();
                        if (api && match && match.value) {{
                            if (match.kind === 'metabolite' && typeof api.addCompound === 'function') {{
                                api.addCompound(match.value);
                            }} else if (match.kind === 'protein' && typeof api.addProtbox === 'function') {{
                                api.addProtbox(match.value);
                            }} else {{
                                renderProteinResults();
                                return;
                            }}
                            proteinSearch.value = '';
                            closeProteinResults();
                            syncState();
                        }} else if (raw) {{
                            renderProteinResults();
                        }}
                    }}
                }});
                proteinSearch?.addEventListener('focus', () => {{
                    if ((proteinSearch?.value || '').trim()) {{
                        scheduleProteinResults();
                    }}
                }});
                createScroll?.addEventListener('scroll', () => {{
                    if (proteinItem?.classList.contains('is-open')) positionProteinResults();
                    if (legendItem?.classList.contains('is-open')) positionLegendPopover();
                }});
                window.addEventListener('resize', () => {{
                    if (proteinItem?.classList.contains('is-open')) positionProteinResults();
                    if (legendItem?.classList.contains('is-open')) positionLegendPopover();
                }});
                if (arrowSelect) {{
                    arrowSelect.addEventListener('change', () => {{
                        const api = getApi();
                        if (api && typeof api.addArrow === 'function' && arrowSelect.value) api.addArrow(arrowSelect.value);
                        arrowSelect.value = '';
                        syncState();
                    }});
                }}
                if (shapeSelect) {{
                    shapeSelect.addEventListener('change', () => {{
                        const api = getApi();
                        if (api && typeof api.addShape === 'function' && shapeSelect.value) api.addShape(shapeSelect.value);
                        shapeSelect.value = '';
                        syncState();
                    }});
                }}
                if (textAdd) textAdd.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.addText === 'function') api.addText(); syncState(); }});
                if (legendAdd) {{
                    legendAdd.addEventListener('click', () => {{
                        const opening = !legendItem?.classList.contains('is-open');
                        if (!opening) {{
                            closeLegendPopover();
                            return;
                        }}
                        positionLegendPopover();
                        if (legendPopover) legendPopover.style.display = 'block';
                        legendItem?.classList.add('is-open');
                    }});
                }}
                if (undoBtn) undoBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.undo === 'function') api.undo(); syncState(); }});
                if (redoBtn) redoBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.redo === 'function') api.redo(); syncState(); }});
                if (deleteBtn) deleteBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.deleteSelected === 'function') api.deleteSelected(); syncState(); }});
                if (edgeShortestBtn) edgeShortestBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.autoConnectShortestEdges === 'function') api.autoConnectShortestEdges(); syncState(); }});
                if (edgeTypeBtn) edgeTypeBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.cycleSelectedEdgeType === 'function') api.cycleSelectedEdgeType(); syncState(); }});
                if (edgeDashBtn) edgeDashBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.toggleSelectedEdgeDash === 'function') api.toggleSelectedEdgeDash(); syncState(); }});
                if (edgeFlipBtn) edgeFlipBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.flipSelectedEdgeDirection === 'function') api.flipSelectedEdgeDirection(); syncState(); }});
                if (selectModeBtn) selectModeBtn.addEventListener('click', () => {{
                    const api = getApi();
                    if (api && typeof api.setMouseMode === 'function' && typeof api.getMouseMode === 'function') {{
                        const nextMode = api.getMouseMode() === 'selection' ? 'drag' : 'selection';
                        api.setMouseMode(nextMode);
                    }}
                    syncState();
                }});
                if (frontBtn) frontBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.bringSelectedToFront === 'function') api.bringSelectedToFront(); syncState(); }});
                if (backBtn) backBtn.addEventListener('click', () => {{ const api = getApi(); if (api && typeof api.sendSelectedToBack === 'function') api.sendSelectedToBack(); syncState(); }});
                proteinItem?.addEventListener('click', (evt) => evt.stopPropagation());
                legendItem?.addEventListener('click', (evt) => evt.stopPropagation());
                panel.querySelectorAll('[data-legend-type]').forEach((el) => {{
                    el.addEventListener('click', () => {{
                        const api = getApi();
                        if (api && typeof api.addLegend === 'function') {{
                            api.addLegend(el.getAttribute('data-legend-type') || 'vertical');
                        }}
                        closeLegendPopover();
                        syncState();
                    }});
                }});
                document.addEventListener('click', () => {{
                    closeProteinResults();
                    closeLegendPopover();
                }});
                window.addEventListener('mk-cst-controls-ready', () => syncState());
                window.addEventListener('mk-cst-viewer-state', (evt) => {{
                    const viewerKey = getViewerKey();
                    if (viewerKey && evt && evt.detail && String(evt.detail.viewerKey || '') === viewerKey) syncState();
                }});
                panel.dataset.bound = '1';
                window.setTimeout(syncState, 50);
            }})();
            """
        ),
    )


def _apply_simple_pathway_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    payload = dict(payload)
    payload["_active_mode"] = "analysis"
    payload.pop("_full_width_canvas", None)
    payload["arrows"] = []
    for group in payload.get("groups", []) or []:
        if isinstance(group, dict):
            group["show_box"] = False
    general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
    general_settings["mode"] = "analysis"
    general_settings["show_background_image"] = True
    general_settings["show_text_boxes"] = False
    return payload


def _home_fake_button(label: Any, icon_name: Optional[str] = None, title: str = "") -> Any:
    contents: List[Any] = []
    if icon_name:
        contents.append(_icon_markup(icon_name, "home-inline-icon"))
    if label not in (None, ""):
        contents.append(ui.tags.span(label))
    attrs = {
        "class": "home-static-button",
        "type": "button",
        "disabled": "disabled",
        "aria-disabled": "true",
    }
    if title:
        attrs["title"] = title
    return ui.tags.button(*contents, attrs)


def _home_shortcut_row(keys: str, action: str, description: str, button: Optional[Any] = None) -> Any:
    button_cell = button if button is not None else ui.div({"class": "home-shortcut-empty"}, "Keyboard only")
    return ui.div(
        {"class": "home-shortcut-row"},
        ui.div({"class": "home-shortcut-keys"}, keys),
        ui.div({"class": "home-shortcut-button-cell"}, button_cell),
        ui.div(
            {"class": "home-shortcut-copy"},
            ui.div({"class": "home-shortcut-action"}, action),
            ui.div({"class": "home-shortcut-description"}, description),
        ),
    )


def _home_bookmark_card(
    title: str,
    subtitle: str,
    details: Sequence[str],
    image_label: Optional[str] = None,
    image_file: Optional[str] = None,
) -> Any:
    image_block: Any = ""
    image_src = _asset_data_uri(image_file) if image_file else ""
    if image_src:
        image_block = ui.div(
            {"class": "home-shot-image-wrap"},
            ui.tags.img(
                {
                    "src": image_src,
                    "alt": f"{title} screenshot",
                    "class": "home-shot-image",
                }
            ),
        )
    elif image_label:
        image_block = ui.div(
            {"class": "home-shot-placeholder"},
            ui.div({"class": "home-shot-placeholder-title"}, "Screenshot Placeholder"),
            ui.div({"class": "home-shot-placeholder-text"}, image_label),
        )
    return ui.div(
        {"class": "home-detail-card"},
        ui.div({"class": "home-detail-card-title"}, title),
        ui.div({"class": "home-detail-card-subtitle"}, subtitle),
        ui.tags.ul(
            {"class": "home-detail-list"},
            *[ui.tags.li(item) for item in details],
        ),
        image_block,
    )


def _home_tab() -> ui.nav_panel:
    hero_gif_src = _asset_data_uri("gif_arrows_add.gif")
    mapkinase_logo_src = _asset_data_uri("MapKinase_logo.png")
    logos_dir = Path(ICONS_DIR) / "Logos"
    combined_logos_src = _asset_data_uri_from_path(logos_dir / "Logos.png")
    return ui.nav_panel(
        "Home",
        ui.div(
            {"class": "home-page-shell"},
            ui.div(
                {"class": "home-hero-card"},
                ui.div(
                    {"class": "home-hero-main"},
                    ui.div(
                        {"class": "home-hero-brand-row"},
                        ui.tags.img(
                            {
                                "src": mapkinase_logo_src,
                                "alt": "MapKinase logo",
                                "class": "home-brand-logo",
                            }
                        ) if mapkinase_logo_src else "",
                    ),
                    ui.div({"class": "home-hero-title"}, "Pathway scoring, editing, and figure design in one workspace."),
                    ui.div(
                        {"class": "home-hero-text"},
                        "MapKinase is an interactive pathway analysis and figure-building environment for protein and PTM datasets. "
                        "It lets you score pathways from uploaded comparisons, open KEGG, WikiPathways, and CST maps, add or edit objects directly in the viewer, "
                        "and turn pathway views into publication-ready figures without leaving the program.",
                    ),
                    ui.div(
                        {"class": "home-hero-text"},
                        "The main workflow is: upload protein and PTM data, choose a pathway source or start from a blank figure, inspect pathway scores, "
                        "and then refine the view with protboxes, PTMs, edges, legends, groups, and text. The viewer supports both mouse-driven editing and keyboard shortcuts for fast work.",
                    ),
                    ui.div(
                        {"class": "home-credit-strip"},
                        ui.div(
                            {"class": "home-credit-copy"},
                            ui.div({"class": "home-credit-label"}, "Credit"),
                            ui.div(
                                {"class": "home-credit-text"},
                                "PTM annotations and kinase-substrate information in MapKinase draw on PhosphoSitePlus resources.",
                            ),
                        ),
                        ui.div(
                            {"class": "home-credit-logo-row"},
                            ui.tags.img(
                                {
                                    "src": combined_logos_src,
                                    "alt": "Project and data source logos",
                                    "class": "home-credit-logo home-credit-logo-combined",
                                }
                            ) if combined_logos_src else "",
                        ),
                    ),
                ),
                ui.div(
                    {"class": "home-hero-side"},
                    ui.div(
                        {"class": "home-hero-media-shell"},
                        ui.div({"class": "home-hero-media-label"}, "Viewer Demo"),
                        ui.div(
                            {"class": "home-hero-media-frame"},
                            ui.tags.img(
                                {
                                    "src": hero_gif_src,
                                    "alt": "MapKinase pathway editing demo",
                                    "class": "home-hero-media",
                                }
                            ) if hero_gif_src else ui.div({"class": "home-shot-placeholder"}, "Add gif_arrows_add.gif here"),
                        ),
                    ),
                    ui.div({"class": "home-side-card-title"}, "Quick Start"),
                    ui.tags.ol(
                        {"class": "home-quickstart-list"},
                        ui.tags.li("Upload a protein dataset and, if available, a PTM dataset in the Input tab."),
                        ui.tags.li("Use Data Analysis Mode (Web Pathways) to load scored KEGG or WikiPathways maps, or Figure Creation to start from a blank canvas."),
                        ui.tags.li("Edit the pathway with Add Objects, Edit Controls, Snapping, and Hide Objects."),
                        ui.tags.li("Export figures or custom pathway layouts once the view looks the way you want."),
                    ),
                ),
            ),
            ui.div(
                {"class": "home-section-grid"},
                _home_bookmark_card(
                    "Data Analysis Mode: Web Pathways",
                    "Search, score, and load KEGG or WikiPathways pathway maps.",
                    [
                        "Pathways are ranked using Over Representation Analysis (ORA) results from the uploaded dataset.",
                        "The Web Pathways viewer includes Full Screen, Hide Objects, Add Objects, Edit Controls, Snapping, and interaction editing tools.",
                        "Data Analysis Mode is the pathway-analysis workflow for loading and editing KEGG, WikiPathways, or CST pathways.",
                    ],
                    image_file="web_pathway (1).png",
                ),
                _home_bookmark_card(
                    "Figure Creation",
                    "Start from a blank viewer and build a custom figure manually.",
                    [
                        "Use this tab when you want to place protboxes, text, shapes, legends, and interactions without starting from an existing pathway.",
                        "The same viewer controls used in Web Pathways are available here, so editing behavior stays consistent.",
                        "Custom layouts can be exported once the figure is complete.",
                    ],
                    image_file="fig_create_pic.png",
                ),
                _home_bookmark_card(
                    "Kinase-Substrate",
                    "Inspect kinase-substrate relationships from the uploaded PTM dataset.",
                    [
                        "The in vivo and in vitro checkboxes determine which evidence types are considered in both the table and the viewer.",
                        "This tab focuses on upstream kinase relationships rather than pathway maps.",
                        "The viewer shares the same editing and visibility controls as the other figure-style tabs.",
                    ],
                    image_file="Substrates_pic.png",
                ),
                _home_bookmark_card(
                    "Input Format",
                    "Protein and PTM uploads follow a header-based format.",
                    [
                        "Protein files use the first column as the protein ID and the second column as the gene label, followed by one or more comparison columns beginning with C:.",
                        "PTM files use the first column as the protein ID and the second column as the PTM site, followed by comparison columns that match the protein comparison headers exactly.",
                        "Optional tool-tip columns begin with T:, and optional Dual-Comparison columns begin with O: and should match the related comparison names.",
                        "Accepted upload types: .txt, .tsv, and .csv.",
                    ],
                ),
            ),
            ui.div(
                {"class": "home-shortcuts-card"},
                ui.div({"class": "home-section-title"}, "Viewer Shortcuts"),
                ui.div(
                    {"class": "home-section-subtitle"},
                    "Keyboard shortcuts below use the same actions as the viewer controls. The button column shows the matching tool when a visible button exists.",
                ),
                ui.div(
                    {"class": "home-shortcut-list"},
                    _home_shortcut_row(
                        "Ctrl + A",
                        "Auto-connect selected protboxes",
                        "Creates standard interactor connections between the selected protboxes using the default auto-connect behavior.",
                    ),
                    _home_shortcut_row(
                        "Ctrl + Shift + A",
                        "Shortest auto-connect",
                        "Adds interactions between selected protboxes using the shortest available path.",
                        _home_fake_button("", "divide-three-solid", "Shortest auto-connect"),
                    ),
                    _home_shortcut_row(
                        "Ctrl + S",
                        "Cycle interactor type",
                        "Cycles the selected interactor, or the interactor between two selected protboxes, through the available interaction types.",
                        _home_fake_button("", "arrow_head_switch", "Cycle interactor type"),
                    ),
                    _home_shortcut_row(
                        "Ctrl + D",
                        "Toggle dashed style",
                        "Switches the selected interactor between solid and dashed styling.",
                        _home_fake_button("", "data-transfer-up", "Toggle interactor dash"),
                    ),
                    _home_shortcut_row(
                        "Ctrl + F",
                        "Flip interactor direction",
                        "Reverses the direction of the selected interactor.",
                        _home_fake_button("", "ruler-arrows", "Flip interactor direction"),
                    ),
                    _home_shortcut_row(
                        "Ctrl + Z",
                        "Undo",
                        "Undoes the most recent viewer action.",
                        _home_fake_button("", "undo", "Undo"),
                    ),
                    _home_shortcut_row(
                        "Ctrl + Y",
                        "Redo",
                        "Redoes the most recently undone viewer action.",
                        _home_fake_button("", "redo", "Redo"),
                    ),
                    _home_shortcut_row(
                        "Delete / Backspace",
                        "Delete selected object",
                        "Deletes the currently selected viewer object when that object can be removed.",
                        _home_fake_button("", "xmark (1)", "Delete selected object"),
                    ),
                    _home_shortcut_row(
                        "Escape",
                        "Clear selection",
                        "Exits group editing if active and clears the current selection.",
                    ),
                ),
            ),
        ),
        value="home",
    )


def _ptm_shape_picker_for(input_id: str = "input_ptm_shape", width: str = "120px") -> Any:
    return ui.input_select(
        input_id,
        None,
        choices={
            "circle": "Circle",
            "square": "Square",
            "diamond": "Diamond",
        },
        selected="circle",
        width=width,
    )


def _ptm_shape_picker() -> Any:
    return _ptm_shape_picker_for("input_ptm_shape", "120px")


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
    settings_col = 4
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
    ui.tags.script(
        f"window.MAPKINASE_DEBUG_EXPORTS = {'true' if DEBUG_EXPORT_ENABLED else 'false'};"
    ),
    ui.tags.style(
        """
        .input-bg { background: #0b4f9c; min-height: 100vh; padding: 32px 24px; }
        .input-wrapper { min-height: 100vh; display: flex; align-items: flex-start; justify-content: flex-start; padding-top: 16px; }
        .data-input-card { max-width: 520px; width: 100%; box-shadow: 0 18px 40px rgba(0,0,0,0.25); border-radius: 16px; }
        .home-page-shell { min-height: 100vh; padding: 28px 18px 40px; background:
            radial-gradient(circle at top left, rgba(14,116,144,0.18), transparent 30%),
            radial-gradient(circle at top right, rgba(30,64,175,0.16), transparent 28%),
            linear-gradient(180deg, #f8fbff 0%, #eef4fb 100%); }
        .home-hero-card, .home-shortcuts-card, .home-detail-card { background: rgba(255,255,255,0.94); border: 1px solid rgba(148,163,184,0.22);
            box-shadow: 0 24px 54px rgba(15,23,42,0.08); }
        .home-hero-card { border-radius: 26px; padding: 28px; display: grid; grid-template-columns: minmax(0, 1.8fr) minmax(280px, 0.9fr); gap: 22px; margin-bottom: 22px; }
        .home-hero-main { display: flex; flex-direction: column; min-height: 100%; }
        .home-hero-brand-row { display: flex; align-items: center; justify-content: space-between; gap: 16px; margin-bottom: 8px; }
        .home-eyebrow { display: inline-flex; align-items: center; gap: 8px; font-size: 12px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase;
            color: #0f766e; background: rgba(204,251,241,0.9); border-radius: 999px; padding: 8px 12px; margin-bottom: 14px; }
        .home-brand-logo { display: block; max-width: 190px; width: 34%; min-width: 120px; height: auto; object-fit: contain; }
        .home-hero-title { font-size: clamp(28px, 3.1vw, 44px); line-height: 1.02; font-weight: 900; color: #0f172a; max-width: 860px; margin-bottom: 18px; }
        .home-hero-text { color: #334155; font-size: 15px; line-height: 1.72; max-width: 900px; margin-bottom: 12px; }
        .home-credit-strip { margin-top: auto; border-radius: 18px; border: 1px solid rgba(148,163,184,0.26); background: #ffffff;
            padding: 12px 16px; display: flex; align-items: stretch; justify-content: space-between; gap: 18px; min-height: 150px; }
        .home-credit-label { font-size: 11px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #0f766e; }
        .home-credit-copy { min-width: 0; max-width: 360px; align-self: flex-start; }
        .home-credit-text { color: #475569; font-size: 13px; line-height: 1.55; margin-top: 4px; }
        .home-credit-logo-row { display: flex; align-items: stretch; justify-content: flex-end; gap: 16px; flex-shrink: 0; }
        .home-credit-logo { display: block; height: auto; object-fit: contain; }
        .home-credit-logo-combined { max-width: 460px; width: 100%; min-width: 220px; height: calc(100% - 8px); margin: 4px 0; object-fit: contain; }
        .home-hero-side { border-radius: 22px; padding: 22px; background: linear-gradient(180deg, #eff6ff 0%, #e0f2fe 100%); border: 1px solid rgba(96,165,250,0.25); }
        .home-hero-media-shell { margin-bottom: 18px; }
        .home-hero-media-label { font-size: 11px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #0369a1; margin-bottom: 8px; }
        .home-hero-media-frame { border-radius: 18px; overflow: hidden; border: 1px solid rgba(148,163,184,0.35); background:
            linear-gradient(135deg, rgba(15,23,42,0.10), rgba(255,255,255,0.72)); box-shadow: inset 0 1px 0 rgba(255,255,255,0.55); }
        .home-hero-media { display: block; width: 100%; height: auto; object-fit: cover; background: #ffffff; }
        .home-side-card-title, .home-section-title { font-size: 16px; font-weight: 800; color: #0f172a; margin-bottom: 8px; }
        .home-quickstart-list { margin: 0; padding-left: 18px; color: #334155; line-height: 1.6; font-size: 13px; }
        .home-section-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 18px; margin-bottom: 22px; }
        .home-detail-card { border-radius: 22px; padding: 20px; }
        .home-detail-card-title { font-size: 18px; font-weight: 800; color: #0f172a; margin-bottom: 8px; }
        .home-detail-card-subtitle { color: #0f4c81; font-size: 14px; font-weight: 700; margin-bottom: 10px; }
        .home-detail-list { margin: 0; padding-left: 18px; color: #334155; line-height: 1.65; }
        .home-detail-list li + li { margin-top: 6px; }
        .home-shot-image-wrap { margin-top: 16px; border-radius: 18px; overflow: hidden; border: 1px solid rgba(148,163,184,0.32);
            background: linear-gradient(180deg, #ffffff, #f8fafc); box-shadow: inset 0 1px 0 rgba(255,255,255,0.65); }
        .home-shot-image { display: block; width: 100%; height: auto; object-fit: cover; }
        .home-shot-placeholder { margin-top: 16px; border-radius: 18px; border: 1px dashed #94a3b8; min-height: 160px; background:
            linear-gradient(135deg, rgba(226,232,240,0.55), rgba(248,250,252,0.95)); display: flex; flex-direction: column; align-items: center; justify-content: center;
            text-align: center; padding: 18px; color: #475569; }
        .home-shot-placeholder-title { font-weight: 800; color: #0f172a; margin-bottom: 6px; }
        .home-shot-placeholder-text { max-width: 220px; line-height: 1.55; }
        .home-shortcuts-card { border-radius: 26px; padding: 24px; }
        .home-section-subtitle { color: #475569; line-height: 1.65; margin-bottom: 16px; }
        .home-shortcut-list { display: flex; flex-direction: column; gap: 10px; }
        .home-shortcut-row { display: grid; grid-template-columns: 180px 92px minmax(0, 1fr); align-items: center; gap: 12px; border-radius: 18px;
            background: linear-gradient(180deg, #f8fafc 0%, #f1f5f9 100%); border: 1px solid #dbe5f0; padding: 12px 14px; }
        .home-shortcut-keys { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 13px; font-weight: 800; color: #0f172a;
            background: #ffffff; border: 1px solid #dbe5f0; border-radius: 12px; padding: 10px 12px; text-align: center; }
        .home-shortcut-button-cell { display: flex; justify-content: center; }
        .home-shortcut-copy { min-width: 0; }
        .home-shortcut-action { font-weight: 800; color: #0f172a; margin-bottom: 4px; }
        .home-shortcut-description { color: #475569; line-height: 1.55; }
        .home-static-button { min-width: 44px; height: 38px; border: 1px solid #cbd5e1; border-radius: 12px; background: linear-gradient(180deg, #ffffff, #eef2f7);
            color: #0f172a; display: inline-flex; align-items: center; justify-content: center; gap: 6px; padding: 0 10px; font-size: 12px; font-weight: 800; opacity: 1; }
        .home-shortcut-empty { color: #64748b; font-size: 12px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em; }
        .home-inline-icon { width: 18px; height: 18px; display: inline-block; color: currentColor; }
        .ptm-shape-picker { position: relative; width: 42px; min-width: 42px; }
        .ptm-shape-picker-trigger { width: 42px; height: 38px; border: 1px solid #cbd5e1; border-radius: 12px;
            background: linear-gradient(180deg, #ffffff, #eef2f7); color: #0f172a; display: inline-flex; align-items: center; justify-content: center;
            padding: 0; cursor: pointer; box-shadow: inset 0 1px 0 rgba(255,255,255,0.8); }
        .ptm-shape-picker-trigger:hover, .ptm-shape-picker.is-open .ptm-shape-picker-trigger { border-color: #94a3b8; background: linear-gradient(180deg, #f8fafc, #e2e8f0); }
        .ptm-shape-trigger-icon { display: none; align-items: center; justify-content: center; }
        .ptm-shape-trigger-icon.is-active { display: inline-flex; }
        .ptm-shape-picker-menu { position: absolute; top: calc(100% + 6px); right: 0; display: none; flex-direction: column; gap: 6px;
            padding: 8px; border-radius: 14px; border: 1px solid #d7dee8; background: rgba(255,255,255,0.98);
            box-shadow: 0 18px 44px rgba(15,23,42,0.18); z-index: 20; }
        .ptm-shape-picker.is-open .ptm-shape-picker-menu { display: flex; }
        .ptm-shape-picker-option { width: 42px; height: 38px; border: 1px solid #d7dee8; border-radius: 12px;
            background: #f8fafc; color: #0f172a; display: inline-flex; align-items: center; justify-content: center; padding: 0; cursor: pointer; }
        .ptm-shape-picker-option:hover { background: #e2e8f0; border-color: #94a3b8; }
        .ptm-shape-picker-icon { width: 18px; height: 18px; display: inline-block; color: currentColor; }
        @media (max-width: 900px) {
            .home-hero-card { grid-template-columns: 1fr; }
            .home-hero-brand-row { flex-direction: column; align-items: flex-start; }
            .home-brand-logo { width: auto; max-width: 180px; }
            .home-shortcut-row { grid-template-columns: 1fr; align-items: flex-start; }
            .home-shortcut-button-cell { justify-content: flex-start; }
            .home-credit-strip { min-height: 0; flex-direction: column; align-items: flex-start; }
            .home-credit-logo-row { flex-wrap: wrap; justify-content: flex-start; }
            .home-credit-logo-combined { width: auto; max-width: 320px; }
        }
        #bookmark_selector a.nav-link[data-value="simple"] { display: none; }
        """
    ),
    NAV_LOCK_SCRIPT,
    EXPORT_DOWNLOAD_SCRIPT,
    INPUT_BUSY_SCRIPT,
    INLINE_HELP_TOOLTIP_SCRIPT,
    GRADIENT_TABLE_SCRIPT,
    SVG_OUTLINE_FONT_BOOTSTRAP_SCRIPT,
    SVG_DOWNLOAD_SCRIPT,
    ui.navset_tab(
        _home_tab(),
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
                        {"class": "mk-mode-help-row"},
                        ui.input_checkbox(
                            "settings_temporal_mode",
                            ui.TagList(
                                "Multi-Comparison Mode",
                                ui.tags.span(
                                    {"class": "mk-inline-help-wrap"},
                                    ui.tags.span(
                                        {
                                            "class": "mk-inline-help",
                                            "aria-label": "Multi-Comparison Mode help",
                                            "tabindex": "0",
                                            "data-help-tooltip-html": (
                                                "<strong>When on:</strong> nodes and PTMs display all comparison values as segmented shapes for multi-comparison datasets. "
                                                "Comparisons are shown left-to-right in the same order as the Protein file comparison columns (C: columns).<br/><br/>"
                                                "<strong>When off:</strong> nodes and PTMs use the single-comparison style."
                                            ),
                                        },
                                        "i",
                                    ),
                                ),
                            ),
                            value=DEFAULT_SETTINGS.get("temporal_mode", False),
                        ),
                    ),
                    ui.div(
                        {"class": "mk-mode-help-row"},
                        ui.input_radio_buttons(
                            "settings_protein_match_mode",
                            ui.TagList(
                                "Web pathway protein matching",
                                ui.tags.span(
                                    {"class": "mk-inline-help-wrap"},
                                    ui.tags.span(
                                        {
                                            "class": "mk-inline-help",
                                            "aria-label": "Web pathway protein matching help",
                                            "tabindex": "0",
                                            "data-help-tooltip-html": (
                                                "<strong>UniProt ID:</strong> uses the existing pathway identifier matching. "
                                                "WikiPathways nodes are resolved to UniProt IDs, while KEGG node identifiers use the uploaded KEGG annotation.<br/><br/>"
                                                "<strong>Gene symbol:</strong> first compares each visible pathway node label with the uploaded Gene Symbol column, ignoring letter case. "
                                                "When no label in a protein box matches, it falls back to the same identifier/UniProt matching.<br/><br/>"
                                                "This setting only affects KEGG and WikiPathways loaded in Web Pathways."
                                            ),
                                        },
                                        "i",
                                    ),
                                ),
                            ),
                            choices={
                                "uniprot_id": "UniProt ID",
                                "gene_symbol": "Gene symbol",
                            },
                            selected=DEFAULT_SETTINGS.get("protein_match_mode", "uniprot_id"),
                            inline=True,
                        ),
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
                        value=DEFAULT_SETTINGS.get("use_original_protbox_size", True),
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
                    *(
                        [
                            ui.input_checkbox(
                                "settings_debug_mode",
                                "Debug Mode",
                                value=DEFAULT_SETTINGS.get("debug_mode", False),
                            ),
                        ]
                        if DEBUG_UI_ENABLED
                        else []
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
                    ui.output_ui("settings_gradient_table"),
                    ui.div(
                        {"style": "display:none;"},
                        ui.input_text("settings_gradient_stops_json", "Gradient stops payload", value="", width="100%"),
                        ui.input_numeric(
                            "settings_max_negative",
                            "Max negative",
                            value=DEFAULT_SETTINGS["max_negative"],
                            step=0.1,
                            width="100%",
                        ),
                        _color_picker_input(
                            "settings_negative_color",
                            "Negative color",
                            _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]),
                        ),
                        ui.input_numeric(
                            "settings_max_positive",
                            "Max positive",
                            value=DEFAULT_SETTINGS["max_positive"],
                            step=0.1,
                            width="100%",
                        ),
                        _color_picker_input(
                            "settings_positive_color",
                            "Positive color",
                            _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]),
                        ),
                    ),
                    ui.output_ui("settings_gradient_preview"),
                    ),
                ),
            ),
            value="settings",
        ),
        id="bookmark_selector",
        selected="home",
    ),
)


def server(input, output, session):  # type: ignore[override]
    session_workspace = get_session_workspace(session)
    print(f"Session workspace initialized: {session_workspace}")
    try:
        session.on_ended(lambda: cleanup_session_workspace(session))
    except Exception as exc:
        print(f"Warning: failed to register session cleanup hook: {exc}")

    bookmark_state: Dict[str, Dict[str, reactive.Value]] = {}
    for cfg in BOOKMARK_CONFIGS:
        bookmark_state[cfg["key"]] = {
            "json": reactive.Value(None),
            "status": reactive.Value(STATUS_READY),
            "loading": reactive.Value(False),
            "load_feedback": reactive.Value(""),
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
    metabolite_validation = reactive.Value(
        {"status": "Metabolite upload optional. Provide after protein if available.", "errors": [], "valid": False, "comparisons": []}  # type: ignore[var-annotated]
    )
    nav_lock_status = reactive.Value("User mode: upload and validate files, then press Run to unlock other tabs.")
    validated_protein_dataset = reactive.Value(None)
    validated_ptm_dataset = reactive.Value(None)
    validated_metabolite_dataset = reactive.Value(None)
    protein_dataset = reactive.Value(None)
    ptm_dataset = reactive.Value(None)
    metabolite_dataset = reactive.Value(None)
    protein_preview_dataset = reactive.Value(None)
    ptm_preview_dataset = reactive.Value(None)
    metabolite_preview_dataset = reactive.Value(None)
    global_catalog_info = reactive.Value(dict(GLOBAL_CATALOG_INFO))
    session_workspace_path = str(session_workspace)
    protein_dataset_path = reactive.Value(None)
    ptm_dataset_path = reactive.Value(None)
    metabolite_dataset_path = reactive.Value(None)
    protein_upload_size_bytes = reactive.Value(0)
    ptm_upload_size_bytes = reactive.Value(0)
    metabolite_upload_size_bytes = reactive.Value(0)
    psp_cache: Dict[str, Any] = {}
    kegg_cache: Dict[str, Dict[str, str]] = {}
    protein_kegg_warning = reactive.Value("")
    ks_index = reactive.Value(_empty_ks_index())
    pathway_score_cache = reactive.Value(
        {
            "status": "Pathway scoring pending.",
            "species_code": "",
            "index_file": "",
            "index_files": {},
            "fc_columns": [],
            "selected_fc": "",
            "results_by_fc": {},
            "download_rows_by_fc": {},
            "significance_mode": "both",
            "positive_cutoff": 1.5,
            "negative_cutoff": -1.5,
            "updated_at": 0.0,
        }
    )
    pipeline_ready = reactive.Value(False)
    pipeline_running = reactive.Value(False)
    mode_sync_in_progress = reactive.Value(False)
    last_applied_input_mode = reactive.Value(None)
    demo_species_sync_target = reactive.Value(None)
    run_guard_message = reactive.Value("")
    last_accepted_run_monotonic = reactive.Value(None)
    run_attempt_times_monotonic = reactive.Value([])
    cst_session_state_paths: Dict[str, str] = {}

    def _safe_session_id() -> str:
        raw_id = str(getattr(session, "id", "") or getattr(session, "_session_id", "") or "").strip()
        if not raw_id:
            return "unknown"
        return hashlib.sha256(raw_id.encode("utf-8")).hexdigest()[:12]

    def _current_session_workspace() -> Path:
        nonlocal session_workspace_path
        raw = str(session_workspace_path or "").strip()
        if raw:
            return Path(raw)
        workspace = get_session_workspace(session)
        session_workspace_path = str(workspace)
        return workspace

    def _get_cst_session_state_path(pathway_id: Any) -> Optional[Path]:
        key = str(pathway_id or "").strip().lower()
        if not key:
            return None
        raw = str(cst_session_state_paths.get(key) or "").strip()
        return Path(raw) if raw else None

    def _set_cst_session_state_path(pathway_id: Any, path: Path) -> None:
        key = str(pathway_id or "").strip().lower()
        if not key:
            return
        cst_session_state_paths[key] = str(path)

    def _resolve_cst_catalog_entry(pathway_id: Any) -> Optional[Dict[str, Any]]:
        key = str(pathway_id or "").strip().lower()
        if not key:
            return None
        try:
            for row in get_cst_pathway_catalog(Path(BASE_DIR)):
                row_id = str(row.get("id") or "").strip().lower()
                if row_id == key:
                    return dict(row)
        except Exception as exc:
            print(f"Warning: failed to resolve CST pathway catalog entry: {exc}")
        return None

    def _safe_source_label(path_or_name: Any) -> str:
        raw = str(path_or_name or "").strip()
        if not raw:
            return ""
        return Path(raw).name

    def _upload_file_size_bytes(path: str) -> int:
        try:
            return int(Path(path).stat().st_size)
        except OSError:
            return 0

    def _extract_upload_datapath(file_info: Any) -> str:
        if isinstance(file_info, dict):
            return str(file_info.get("datapath") or "")
        return str(getattr(file_info, "datapath", "") or "")

    def _cleanup_upload_temp_file(path: str, reason: str) -> None:
        safely_delete_temp_file(path, reason=reason)

    def _write_session_preview_json(payload: Dict[str, Any]) -> None:
        try:
            preview_path = safe_session_path(session, SESSION_PREVIEW_FILENAME)
            with open(preview_path, "w", encoding="utf-8") as preview_fh:
                json.dump(payload, preview_fh, indent=2)
        except OSError as write_err:
            print(f"Warning: could not write session preview JSON: {write_err}")

    def _clear_pathway_scores(status: str = "Pathway scoring pending.") -> None:
        pathway_score_cache.set(
            {
                "status": status,
                "species_code": "",
                "index_file": "",
                "index_files": {},
                "fc_columns": [],
                "selected_fc": "",
                "results_by_fc": {},
                "download_rows_by_fc": {},
                "significance_mode": "both",
                "positive_cutoff": 1.5,
                "negative_cutoff": -1.5,
                "updated_at": time.time(),
            }
        )

    def _mark_pathway_scores_pending(status: str = "Run Fisher's Exact Test to score pathways.") -> None:
        _clear_pathway_scores(status)

    def _status_has_blocking_error(status: Any) -> bool:
        text = str(status or "").strip().lower()
        if not text:
            return False
        return any(token in text for token in ("failed", "missing", "could not", "error"))

    def _validation_has_explicit_error(state: Optional[Dict[str, Any]]) -> bool:
        payload = dict(state or {})
        if payload.get("errors"):
            return True
        return _status_has_blocking_error(payload.get("status"))

    def _validation_has_blocking_error(state: Optional[Dict[str, Any]], optional: bool = False) -> bool:
        payload = dict(state or {})
        if _validation_has_explicit_error(payload):
            return True
        if optional and not payload.get("valid"):
            return False
        return not bool(payload.get("valid"))

    def _current_run_upload_sizes() -> List[Tuple[str, int]]:
        sizes: List[Tuple[str, int]] = []
        protein_size = int(protein_upload_size_bytes.get() or 0)
        ptm_size = int(ptm_upload_size_bytes.get() or 0)
        metabolite_size = int(metabolite_upload_size_bytes.get() or 0)
        if protein_size > 0:
            sizes.append(("protein", protein_size))
        if ptm_size > 0:
            sizes.append(("ptm", ptm_size))
        if metabolite_size > 0:
            sizes.append(("metabolite", metabolite_size))
        return sizes

    def _clear_live_user_datasets(score_status: str = "Validation complete. Click Run to process datasets.") -> None:
        protein_dataset.set(None)
        ptm_dataset.set(None)
        metabolite_dataset.set(None)
        protein_kegg_warning.set("")
        _reset_global_catalog_from_default()
        _write_debug_dump("user_protein_dataset_debug.txt", {"info": "Processed dataset not generated yet"})
        _write_debug_dump("user_ptm_dataset_debug.txt", {"info": "Processed dataset not generated yet"})
        _update_ks_index(reset=True)
        _clear_pathway_scores(score_status)

    def _invalidate_user_run(score_status: str = "Validation complete. Click Run to process datasets.") -> None:
        if _current_mode() == "demo":
            return
        pipeline_ready.set(False)
        pipeline_running.set(False)
        run_guard_message.set("")
        _clear_live_user_datasets(score_status)

    def _dataset_to_df(dataset: Optional[Dict[str, Any]]) -> pd.DataFrame:
        if not dataset:
            return pd.DataFrame()
        headers = list(dataset.get("headers") or [])
        rows = list(dataset.get("rows") or [])
        if not headers:
            return pd.DataFrame()
        if not rows:
            return pd.DataFrame(columns=headers)
        normalized_rows: List[List[str]] = []
        width = len(headers)
        for row in rows:
            vals = list(row)
            if len(vals) < width:
                vals.extend([""] * (width - len(vals)))
            elif len(vals) > width:
                vals = vals[:width]
            normalized_rows.append(vals)
        return pd.DataFrame(normalized_rows, columns=headers, dtype=str)


    def _fisher_right_tail(total_n: int, pathway_n: int, sig_n: int, sig_in_pathway: int) -> float:
        if total_n <= 0 or pathway_n < 0 or sig_n < 0 or sig_in_pathway < 0:
            return 1.0
        if pathway_n > total_n or sig_n > total_n:
            return 1.0
        max_hits = min(pathway_n, sig_n)
        if sig_in_pathway > max_hits:
            return 1.0
        min_hits = max(0, pathway_n - (total_n - sig_n))
        if sig_in_pathway < min_hits:
            return 1.0

        def _log_choose(n: int, k: int) -> float:
            if k < 0 or k > n:
                return float("-inf")
            return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

        denom = _log_choose(total_n, pathway_n)
        if not math.isfinite(denom):
            return 1.0
        p_val = 0.0
        for hit_count in range(sig_in_pathway, max_hits + 1):
            log_prob = _log_choose(sig_n, hit_count) + _log_choose(total_n - sig_n, pathway_n - hit_count) - denom
            if math.isfinite(log_prob):
                p_val += math.exp(log_prob)
        return min(max(p_val, 0.0), 1.0)

    def _benjamini_hochberg(values: List[Optional[float]]) -> List[Optional[float]]:
        indexed = [(idx, float(val)) for idx, val in enumerate(values) if val is not None and math.isfinite(float(val))]
        if not indexed:
            return [None] * len(values)
        indexed.sort(key=lambda item: item[1])
        adjusted: List[Optional[float]] = [None] * len(values)
        total = len(indexed)
        running = 1.0
        for rev_rank, (original_idx, p_val) in enumerate(reversed(indexed), start=1):
            rank = total - rev_rank + 1
            bh_val = min(running, (p_val * total) / rank)
            running = bh_val
            adjusted[original_idx] = min(max(bh_val, 0.0), 1.0)
        return adjusted

    def _coerce_numeric_series(series: pd.Series) -> pd.Series:
        return pd.to_numeric(series.astype(str).str.strip(), errors="coerce")

    def _fisher_significance_mask(
        series: pd.Series,
        positive_cutoff: float,
        negative_cutoff: float,
        mode: str = "both",
    ) -> pd.Series:
        active_mode = str(mode or "both").strip().lower()
        if active_mode == "positive":
            return series >= positive_cutoff
        if active_mode == "negative":
            return series <= negative_cutoff
        return (series >= positive_cutoff) | (series <= negative_cutoff)

    def _pathway_uniprots(pathway: Dict[str, Any], index_nodes: Dict[str, Dict[str, Any]], gene_to_uniprot: Dict[str, List[str]]) -> set[str]:
        pathway_uniprots: set[str] = set()
        modules = list(pathway.get("modules", []))
        if modules:
            prebuilt_uniprots = list(pathway.get("resolved_uniprots") or [])
            prebuilt_gene_symbols = list(pathway.get("gene_symbols") or [])
            if prebuilt_uniprots or prebuilt_gene_symbols:
                for uni in prebuilt_uniprots:
                    normalized = normalize_uniprot(uni)
                    if normalized:
                        pathway_uniprots.add(normalized)
                for gene_symbol in prebuilt_gene_symbols:
                    for uni in gene_to_uniprot.get(str(gene_symbol or "").strip().upper(), []):
                        normalized = normalize_uniprot(uni)
                        if normalized:
                            pathway_uniprots.add(normalized)
                return pathway_uniprots
            for module in modules:
                if not isinstance(module, dict):
                    continue
                raw_ids = list(module.get("uniprot_ids") or [])
                if raw_ids:
                    for uni in raw_ids:
                        normalized = normalize_uniprot(uni)
                        if normalized:
                            pathway_uniprots.add(normalized)
                for gene_symbol in list(module.get("gene_symbols") or []):
                    for uni in gene_to_uniprot.get(str(gene_symbol or "").strip().upper(), []):
                        normalized = normalize_uniprot(uni)
                        if normalized:
                            pathway_uniprots.add(normalized)
            return pathway_uniprots
        for node_id in list(pathway.get("nodes", [])):
            node_obj = index_nodes.get(str(node_id), {})
            for uni in candidate_uniprots_for_node(node_obj, gene_to_uniprot):
                normalized = normalize_uniprot(uni)
                if normalized:
                    pathway_uniprots.add(normalized)
        return pathway_uniprots

    def _cst_pathway_id_from_name(name: str) -> str:
        return "cst_" + re.sub(r"[^a-z0-9]+", "_", str(name or "").strip().lower()).strip("_")

    def _build_cst_fisher_index(base_dir: Path) -> Dict[str, Any]:
        grouped: Dict[str, Dict[str, Any]] = {}
        for row in iter_effective_cst_modules(base_dir):
            if not isinstance(row, dict):
                continue
            pathway_name = str(row.get("pathway") or "").strip()
            protein_name = str(row.get("protein_module_name") or "").strip()
            if not pathway_name or not protein_name:
                continue
            pathway_id = _cst_pathway_id_from_name(pathway_name)
            pathway_entry = grouped.setdefault(
                pathway_id,
                {
                    "pathway_id": pathway_id,
                    "name": pathway_name,
                    "modules": [],
                },
            )
            uniprot_ids: List[str] = []
            for field_name in ("psp_uniprot_ids", "backup_uniprot_ids"):
                raw_val = str(row.get(field_name) or "")
                for item in raw_val.split(";"):
                    normalized = normalize_uniprot(item)
                    if normalized and normalized not in uniprot_ids:
                        uniprot_ids.append(normalized)
            gene_symbols: List[str] = []
            raw_genes = str(row.get("gene_symbols") or "")
            for item in raw_genes.split(";"):
                symbol = str(item or "").strip().upper()
                if symbol and symbol not in gene_symbols:
                    gene_symbols.append(symbol)
            pathway_entry["modules"].append(
                {
                    "label": protein_name,
                    "uniprot_ids": uniprot_ids,
                    "gene_symbols": gene_symbols,
                }
            )
        return {
            "pathways": list(grouped.values()),
            "nodes": {},
        }

    def _compute_fisher_pathway_rows(
        prot_df: pd.DataFrame,
        site_df: Optional[pd.DataFrame],
        fc_col: str,
        positive_cutoff: float,
        negative_cutoff: float,
        significance_mode: str,
        index_sources: List[Tuple[str, Dict[str, Any], str]],
        gene_map: Dict[str, List[str]],
        site_fc_col: Optional[str] = None,
    ) -> Tuple[Dict[str, Dict[str, Dict[str, Any]]], List[Dict[str, Any]]]:
        protein_rows = prot_df.copy()
        protein_id_col = protein_rows.columns[0]
        protein_rows["_uniprot"] = protein_rows[protein_id_col].map(normalize_uniprot)
        protein_rows = protein_rows[protein_rows["_uniprot"] != ""].copy()
        protein_gene_col = protein_rows.columns[1] if len(protein_rows.columns) > 1 else None
        protein_cst_uniprot_col, protein_cst_gene_col = _resolve_cst_ortholog_columns(list(protein_rows.columns))
        if protein_gene_col:
            protein_rows["_gene_symbol"] = (
                protein_rows[protein_gene_col]
                .astype(str)
                .str.strip()
                .str.upper()
            )
        else:
            protein_rows["_gene_symbol"] = ""
        if protein_cst_uniprot_col:
            protein_rows["_cst_uniprot"] = protein_rows[protein_cst_uniprot_col].map(normalize_uniprot)
        else:
            protein_rows["_cst_uniprot"] = protein_rows["_uniprot"]
        if protein_cst_gene_col:
            protein_rows["_cst_gene_symbol"] = (
                protein_rows[protein_cst_gene_col]
                .astype(str)
                .str.strip()
                .str.upper()
            )
        elif protein_cst_uniprot_col:
            protein_rows["_cst_gene_symbol"] = ""
        else:
            protein_rows["_cst_gene_symbol"] = protein_rows["_gene_symbol"]
        protein_rows["_fc"] = _coerce_numeric_series(protein_rows[fc_col])
        protein_rows["_significant"] = _fisher_significance_mask(
            protein_rows["_fc"],
            positive_cutoff,
            negative_cutoff,
            significance_mode,
        )
        protein_rows = protein_rows.drop_duplicates(subset=["_uniprot"], keep="first")

        protein_uniprots = set(protein_rows["_uniprot"].tolist())
        significant_protein_uniprots = set(protein_rows.loc[protein_rows["_significant"], "_uniprot"].tolist())
        total_proteins = int(len(protein_uniprots))
        significant_proteins = int(len(significant_protein_uniprots)) if total_proteins else 0
        protein_cst_uniprots = set(
            str(item or "").strip()
            for item in protein_rows["_cst_uniprot"].tolist()
            if str(item or "").strip()
        )
        significant_protein_cst_uniprots = set(
            str(item or "").strip()
            for item in protein_rows.loc[protein_rows["_significant"], "_cst_uniprot"].tolist()
            if str(item or "").strip()
        )
        total_cst_proteins = int(len(protein_cst_uniprots))
        significant_cst_proteins = int(len(significant_protein_cst_uniprots)) if total_cst_proteins else 0

        effective_gene_map: Dict[str, List[str]] = {
            str(key or "").strip().upper(): list(dict.fromkeys(str(item or "").strip().upper() for item in values if str(item or "").strip()))
            for key, values in (gene_map or {}).items()
            if str(key or "").strip()
        }
        if protein_gene_col:
            for row in protein_rows.itertuples(index=False):
                row_map = dict(zip(protein_rows.columns, row))
                symbol = str(row_map.get("_gene_symbol") or "").strip().upper()
                uni = normalize_uniprot(row_map.get("_uniprot"))
                if not symbol or not uni:
                    continue
                bucket = effective_gene_map.setdefault(symbol, [])
                if uni not in bucket:
                    bucket.append(uni)
        effective_cst_gene_map: Dict[str, List[str]] = {}
        if not protein_cst_uniprot_col:
            effective_cst_gene_map = {key: list(values) for key, values in effective_gene_map.items()}
        for row in protein_rows.itertuples(index=False):
            row_map = dict(zip(protein_rows.columns, row))
            symbol = str(row_map.get("_cst_gene_symbol") or "").strip().upper()
            uni = normalize_uniprot(row_map.get("_cst_uniprot"))
            if not symbol or not uni:
                continue
            bucket = effective_cst_gene_map.setdefault(symbol, [])
            if uni not in bucket:
                bucket.append(uni)

        site_rows = pd.DataFrame()
        total_sites = 0
        significant_sites = 0
        total_cst_sites = 0
        significant_cst_sites = 0
        if site_df is not None and not site_df.empty:
            site_rows = site_df.copy()
            site_uniprot_col = site_rows.columns[0]
            site_rows["_parent_uniprot"] = site_rows[site_uniprot_col].map(normalize_uniprot)
            site_rows = site_rows[site_rows["_parent_uniprot"] != ""].copy()
            site_cst_uniprot_col, _site_cst_gene_col = _resolve_cst_ortholog_columns(list(site_rows.columns))
            if site_cst_uniprot_col:
                site_rows["_parent_cst_uniprot"] = site_rows[site_cst_uniprot_col].map(normalize_uniprot)
            else:
                site_rows["_parent_cst_uniprot"] = site_rows["_parent_uniprot"]
            active_site_fc_col = site_fc_col if site_fc_col and site_fc_col in site_rows.columns else None
            if active_site_fc_col:
                site_rows["_fc"] = _coerce_numeric_series(site_rows[active_site_fc_col])
                site_rows["_significant"] = _fisher_significance_mask(
                    site_rows["_fc"],
                    positive_cutoff,
                    negative_cutoff,
                    significance_mode,
                )
            else:
                site_rows["_significant"] = False
            total_sites = int(len(site_rows))
            significant_sites = int(site_rows["_significant"].sum()) if total_sites else 0
            cst_site_mask = site_rows["_parent_cst_uniprot"].astype(str).str.strip() != ""
            total_cst_sites = int(cst_site_mask.sum())
            significant_cst_sites = int((site_rows["_significant"] & cst_site_mask).sum()) if total_cst_sites else 0

        combined_rows: List[Dict[str, Any]] = []
        nested_rows: Dict[str, Dict[str, Dict[str, Any]]] = {}
        for source_key, source_index, _source_file in index_sources:
            source_rows: Dict[str, Dict[str, Any]] = {}
            is_cst_source = str(source_key or "").strip().lower() == "cst"
            active_gene_map = effective_cst_gene_map if is_cst_source else effective_gene_map
            active_protein_uniprots = protein_cst_uniprots if is_cst_source else protein_uniprots
            active_significant_protein_uniprots = significant_protein_cst_uniprots if is_cst_source else significant_protein_uniprots
            active_total_proteins = total_cst_proteins if is_cst_source else total_proteins
            active_significant_proteins = significant_cst_proteins if is_cst_source else significant_proteins
            active_site_key = "_parent_cst_uniprot" if is_cst_source else "_parent_uniprot"
            active_total_sites = total_cst_sites if is_cst_source else total_sites
            active_significant_sites = significant_cst_sites if is_cst_source else significant_sites
            for pathway in list(source_index.get("pathways", [])):
                pathway_id = str(pathway.get("pathway_id", "")).strip().lower()
                if not pathway_id:
                    continue
                pathway_name = str(pathway.get("name") or pathway_id).strip()
                pathway_uniprots = _pathway_uniprots(pathway, source_index.get("nodes", {}), active_gene_map)
                proteins_in_dataset = active_protein_uniprots & pathway_uniprots
                sig_proteins_in_path = active_significant_protein_uniprots & pathway_uniprots
                pathway_protein_total = len(proteins_in_dataset)
                sig_protein_in_pathway = len(sig_proteins_in_path)
                protein_p = _fisher_right_tail(active_total_proteins, pathway_protein_total, active_significant_proteins, sig_protein_in_pathway) if active_total_proteins else None

                pathway_site_total = 0
                sig_site_in_pathway = 0
                site_p: Optional[float] = None
                if active_total_sites and not site_rows.empty:
                    valid_site_mask = site_rows[active_site_key].astype(str).str.strip() != ""
                    in_path_mask = valid_site_mask & site_rows[active_site_key].isin(pathway_uniprots)
                    pathway_site_total = int(in_path_mask.sum())
                    if pathway_site_total:
                        sig_site_in_pathway = int((site_rows["_significant"] & in_path_mask).sum())
                    site_p = _fisher_right_tail(active_total_sites, pathway_site_total, active_significant_sites, sig_site_in_pathway)

                row = {
                    "pathway_source": source_key,
                    "pathway_id": pathway_id,
                    "name": pathway_name,
                    "prot_dataset_total": active_total_proteins,
                    "prot_pathway_total": pathway_protein_total,
                    "prot_significant_total": active_significant_proteins,
                    "prot_significant_in_pathway": sig_protein_in_pathway,
                    "prot_fisher_p": protein_p,
                    "phos_dataset_total": active_total_sites,
                    "phos_pathway_total": pathway_site_total,
                    "phos_significant_total": active_significant_sites,
                    "phos_significant_in_pathway": sig_site_in_pathway,
                    "phos_fisher_p": site_p,
                    "positive_cutoff": positive_cutoff,
                    "negative_cutoff": negative_cutoff,
                    "significance_mode": significance_mode,
                    "comparison_column": fc_col,
                }
                combined_rows.append(row)
                source_rows[pathway_id] = row
            nested_rows[source_key] = source_rows

        return nested_rows, combined_rows

    def _resolve_kegg_index_file_for_species(species_code: str) -> str:
        code = (species_code or "").strip().lower()
        if not code:
            return ""
        candidates = [
            os.path.join(INDEX_FILES_DIR, f"kegg_index_{code}.json"),
            os.path.join(INDEX_FILES_DIR, f"{code}_kegg_index.json"),
            os.path.join(INDEX_FILES_DIR, f"kegg_index_{code}_v1.json"),
        ]
        for path in candidates:
            if os.path.exists(path):
                return path
        # Final fallback: if there is only one KEGG index file, use it.
        try:
            all_index_files = [
                os.path.join(INDEX_FILES_DIR, name)
                for name in os.listdir(INDEX_FILES_DIR)
                if name.lower().startswith("kegg_index_") and name.lower().endswith(".json")
            ]
            if len(all_index_files) == 1:
                return all_index_files[0]
        except Exception:
            pass
        return ""

    def _resolve_wikipathways_index_file_for_species(species_code: str) -> str:
        code = (species_code or "").strip().lower()
        if not code:
            return ""
        candidates = [
            os.path.join(INDEX_FILES_DIR, f"wikipathways_index_{code}.json"),
            os.path.join(INDEX_FILES_DIR, f"{code}_wikipathways_index.json"),
            os.path.join(INDEX_FILES_DIR, f"wikipathways_index_{code}_v1.json"),
        ]
        for path in candidates:
            if os.path.exists(path):
                return path
        # Final fallback: if there is only one WikiPathways index file, use it.
        try:
            all_index_files = [
                os.path.join(INDEX_FILES_DIR, name)
                for name in os.listdir(INDEX_FILES_DIR)
                if name.lower().startswith("wikipathways_index_") and name.lower().endswith(".json")
            ]
            if len(all_index_files) == 1:
                return all_index_files[0]
        except Exception:
            pass
        return ""

    def _resolve_cst_index_file_for_species(species_code: str) -> str:
        code = (species_code or "").strip().lower()
        candidates: List[str] = []
        if code:
            candidates.extend([
                os.path.join(INDEX_FILES_DIR, f"cst_index_{code}.json"),
                os.path.join(INDEX_FILES_DIR, f"{code}_cst_index.json"),
            ])
        candidates.append(os.path.join(INDEX_FILES_DIR, "cst_index.json"))
        for path in candidates:
            if os.path.exists(path):
                return path
        try:
            all_index_files = [
                os.path.join(INDEX_FILES_DIR, name)
                for name in os.listdir(INDEX_FILES_DIR)
                if name.lower().startswith("cst_index_") and name.lower().endswith(".json")
            ]
            if len(all_index_files) == 1:
                return all_index_files[0]
        except Exception:
            pass
        return ""

    def _load_cst_fisher_index(path: Path) -> Dict[str, Any]:
        with open(path, "r", encoding="utf-8") as fh:
            payload = json.load(fh)
        if not isinstance(payload, dict):
            raise ValueError(f"Invalid CST Fisher index format in {path}")
        pathways = payload.get("pathways")
        if not isinstance(pathways, list):
            raise KeyError("Pathway index missing key: pathways")
        return payload

    def _resolve_gene_to_uniprot_file_for_species(species_code: str) -> str:
        code = (species_code or "").strip().lower()
        if not code:
            return ""
        candidates = [
            os.path.join(ANNOTATION_FILES_DIR, f"{code}_id_mapping_table.txt"),
            os.path.join(ANNOTATION_FILES_DIR, f"{code}_mapping_table.txt"),
        ]
        for path in candidates:
            if os.path.exists(path):
                return path
        return ""

    def _refresh_pathway_scores(
        selected_fc: Optional[str] = None,
        *,
        species_selection: Optional[str] = None,
    ) -> None:
        with reactive.isolate():
            protein_ok = bool(protein_validation.get().get("valid"))
            protein_data = protein_dataset.get()
            ptm_ok = bool(ptm_validation.get().get("valid"))
            ptm_data = ptm_dataset.get() if ptm_ok else None
            resolved_species_selection = (
                species_selection
                if species_selection is not None
                else _get_input_value(input, "input_species")
            )
            _species_choice, species_info = _resolve_species(resolved_species_selection)
            current_mode = _current_mode()
            significance_mode = str(_get_input_value(input, "web_fisher_sig_mode") or "both").strip().lower()
            if significance_mode not in {"both", "positive", "negative"}:
                significance_mode = "both"
            positive_cutoff = abs(_to_float(_get_input_value(input, "web_fisher_sig_pos"), 1.5))
            negative_cutoff = -abs(_to_float(_get_input_value(input, "web_fisher_sig_neg"), -1.5))
            existing_score_cache = dict(pathway_score_cache.get() or {})

        if not protein_ok or not protein_data:
            _clear_pathway_scores("Pathway scoring waiting for validated inputs and Run.")
            return

        species_code = (species_info.get("code") or "").strip().lower()
        kegg_index_file = _resolve_kegg_index_file_for_species(species_code)
        wikipathways_index_file = _resolve_wikipathways_index_file_for_species(species_code)
        cst_index_file = _resolve_cst_index_file_for_species(species_code)
        cst_ai_dir = Path(BASE_DIR).resolve().parent / "stored_pathways" / "cst"
        cst_available = bool(cst_index_file) or cst_ai_dir.exists()
        if not kegg_index_file and not wikipathways_index_file and not cst_available:
            _clear_pathway_scores(f"No pathway index files found for species '{species_code}' in index_files.")
            return

        prot_df = _dataset_to_df(protein_data)
        site_df = _dataset_to_df(ptm_data) if ptm_data else None
        if prot_df.empty:
            _clear_pathway_scores("Protein dataset is empty after parsing; pathway scoring skipped.")
            return

        fc_columns = [c for c in list(prot_df.columns) if str(c).startswith("C:")]
        if not fc_columns:
            _clear_pathway_scores("No protein main columns (C:) were found for pathway scoring.")
            return

        if selected_fc and selected_fc in fc_columns:
            ordered_fc_cols = [selected_fc]
            active_selected_fc = selected_fc
        else:
            ordered_fc_cols = [fc_columns[0]]
            active_selected_fc = ordered_fc_cols[0]

        try:
            gene_map_file = _resolve_gene_to_uniprot_file_for_species(species_code)
            gene_map = build_gene_to_uniprot_map(Path(gene_map_file)) if gene_map_file else {}
            index_sources: List[Tuple[str, Dict[str, Any], str]] = []
            if kegg_index_file:
                index_sources.append(("kegg", load_kegg_index(Path(kegg_index_file)), kegg_index_file))
            if wikipathways_index_file:
                index_sources.append(("wikipathways", load_kegg_index(Path(wikipathways_index_file)), wikipathways_index_file))
            if cst_index_file:
                index_sources.append(("cst", _load_cst_fisher_index(Path(cst_index_file)), cst_index_file))
            elif cst_ai_dir.exists():
                index_sources.append(("cst", _build_cst_fisher_index(cst_ai_dir), str(cst_ai_dir)))

            same_context = (
                str(existing_score_cache.get("species_code") or "").strip().lower() == species_code
                and str(existing_score_cache.get("significance_mode") or "both").strip().lower() == significance_mode
                and abs(abs(_to_float(existing_score_cache.get("positive_cutoff"), 1.5)) - positive_cutoff) < 1e-9
                and abs((-abs(_to_float(existing_score_cache.get("negative_cutoff"), -1.5))) - negative_cutoff) < 1e-9
            )
            results_by_fc: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = (
                dict(existing_score_cache.get("results_by_fc") or {}) if same_context else {}
            )
            download_rows_by_fc: Dict[str, List[Dict[str, Any]]] = (
                dict(existing_score_cache.get("download_rows_by_fc") or {}) if same_context else {}
            )
            site_headers = list(site_df.columns) if site_df is not None and not site_df.empty else []

            for fc_col in ordered_fc_cols:
                if fc_col in results_by_fc and fc_col in download_rows_by_fc:
                    continue
                site_fc_col = fc_col if site_df is not None and fc_col in site_headers else None
                source_maps, flat_rows = _compute_fisher_pathway_rows(
                    prot_df=prot_df,
                    site_df=site_df,
                    fc_col=fc_col,
                    positive_cutoff=positive_cutoff,
                    negative_cutoff=negative_cutoff,
                    significance_mode=significance_mode,
                    index_sources=index_sources,
                    gene_map=gene_map,
                    site_fc_col=site_fc_col,
                )
                results_by_fc[fc_col] = source_maps
                download_rows_by_fc[fc_col] = flat_rows

            index_files_obj: Dict[str, str] = {}
            if kegg_index_file:
                index_files_obj["kegg"] = kegg_index_file
            if wikipathways_index_file:
                index_files_obj["wikipathways"] = wikipathways_index_file
            if cst_index_file:
                index_files_obj["cst"] = cst_index_file
            elif cst_ai_dir.exists():
                index_files_obj["cst"] = str(cst_ai_dir)
            source_label = ", ".join(sorted(index_files_obj.keys())) if index_files_obj else "none"

            pathway_score_cache.set(
                {
                    "status": (
                        f"Pathway scoring complete ({len(results_by_fc)} of {len(fc_columns)} main columns cached, "
                        f"sources={source_label}, mode={current_mode}, "
                        f"fisher={significance_mode}, cutoffs={positive_cutoff:g}/{negative_cutoff:g})."
                    ),
                    "species_code": species_code,
                    "index_file": kegg_index_file or wikipathways_index_file or "",
                    "index_files": index_files_obj,
                    "fc_columns": list(fc_columns),
                    "selected_fc": active_selected_fc,
                    "results_by_fc": results_by_fc,
                    "download_rows_by_fc": download_rows_by_fc,
                    "significance_mode": significance_mode,
                    "positive_cutoff": positive_cutoff,
                    "negative_cutoff": negative_cutoff,
                    "updated_at": time.time(),
                }
            )
        except Exception as exc:
            print(f"Warning: pathway scoring failed: {exc}")
            _clear_pathway_scores("Pathway scoring failed.")

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
            state["status"].set("Failed to export custom pathway.")
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
            print(f"Warning: dataset preview load failed: {exc}")
            return {"error": "File could not be parsed as tabular text."}

    def _fill_blank_gene_with_uniprot(dataset: Dict[str, Any]) -> Dict[str, Any]:
        """For protein datasets, fill blank gene symbol cells with UniProt IDs."""
        headers = list(dataset.get("headers") or [])
        rows = list(dataset.get("rows") or [])
        if len(headers) < 2 or not rows:
            return dataset
        for row in rows:
            if len(row) < 2:
                row.extend([""] * (2 - len(row)))
            uniprot = str(row[0] or "").strip()
            gene = str(row[1] or "").strip()
            if uniprot and not gene:
                row[1] = uniprot
        dataset["rows"] = rows
        return dataset


    def _validate_metabolite_file(path: str, protein_comparisons: Sequence[str]) -> Namespace:
        errors: List[str] = []
        dataset = _load_dataset(path)
        if dataset.get("error"):
            errors.append(str(dataset.get("error")))
            return Namespace(valid=False, errors=errors, summary={"rows": 0, "comparisons": 0, "tooltips": 0}, comparisons=[])
        headers = [str(h or "").strip() for h in list(dataset.get("headers") or [])]
        rows = list(dataset.get("rows") or [])
        if not headers:
            errors.append("No headers found in metabolite file.")
            return Namespace(valid=False, errors=errors, summary={"rows": 0, "comparisons": 0, "tooltips": 0}, comparisons=[])
        fc_columns = [h for h in headers if str(h).startswith("C:")]
        tooltip_columns = [h for h in headers if str(h).startswith("T:")]
        if not fc_columns:
            errors.append("Metabolite file requires at least one comparison column beginning with 'C:'.")
        expected = [str(c or "").strip() for c in list(protein_comparisons or []) if str(c or "").strip()]
        if expected:
            missing = [col for col in expected if col not in fc_columns]
            extra = [col for col in fc_columns if col not in expected]
            if missing:
                errors.append("Missing Comparison columns present in protein file: " + ", ".join(missing) + ".")
            if extra:
                errors.append("Metabolite comparison columns must match the Protein dataset exactly. Unexpected columns: " + ", ".join(extra) + ".")
        summary = {
            "rows": len(rows),
            "comparisons": len(fc_columns),
            "tooltips": len(tooltip_columns),
        }
        return Namespace(valid=not errors, errors=errors, summary=summary, comparisons=fc_columns)

    def _get_input_preview_content(
        dataset: Optional[Dict[str, Any]],
        fallback_headers: List[str],
    ) -> Tuple[List[str], List[str]]:
        preview_char_limit = 20
        headers: List[str] = []
        row_values: List[str] = []
        if dataset and not dataset.get("error"):
            raw_headers = dataset.get("headers") or []
            raw_rows = dataset.get("rows") or []
            first_data_row = raw_rows[0] if raw_rows else []
            width = max(len(raw_headers), len(first_data_row))
            for idx in range(width):
                header = str(raw_headers[idx]).strip() if idx < len(raw_headers) else ""
                value = str(first_data_row[idx]).strip() if idx < len(first_data_row) else ""
                if header or value:
                    headers.append(header or f"Column {idx + 1}")
                    row_values.append(value[:preview_char_limit])
        if not headers:
            headers = list(fallback_headers)
            row_values = [""] * len(headers)
        return headers, row_values

    def _get_input_preview_widths(*table_contents: Tuple[List[str], List[str]]) -> List[str]:
        max_columns = max((len(headers) for headers, _ in table_contents), default=0)
        widths: List[str] = []
        for idx in range(max_columns):
            placeholder_only = True
            max_chars = 10
            for headers, row_values in table_contents:
                if idx < len(headers):
                    header_text = str(headers[idx]).strip()
                    max_chars = max(max_chars, len(header_text))
                    if header_text != "...":
                        placeholder_only = False
                if idx < len(row_values):
                    cell_text = str(row_values[idx]).strip()
                    max_chars = max(max_chars, len(cell_text))
                    if cell_text:
                        placeholder_only = False
            if placeholder_only:
                widths.append("5ch")
            else:
                widths.append(f"{min(max(max_chars + 2, 10), 24)}ch")
        return widths

    def _extract_problem_column_indexes(headers: List[str], errors: List[str]) -> List[int]:
        if not headers or not errors:
            return []
        problem_indexes: set[int] = set()
        header_map = {str(h).strip(): idx for idx, h in enumerate(headers)}
        for err in errors:
            text = str(err or "").strip()
            if not text:
                continue
            col_match = re.search(r"\bcolumn\s+(\d+)\b", text, flags=re.IGNORECASE)
            if col_match:
                idx = int(col_match.group(1)) - 1
                if 0 <= idx < len(headers):
                    problem_indexes.add(idx)
            invalid_header_match = re.search(r"Invalid header in column\s+(\d+)", text, flags=re.IGNORECASE)
            if invalid_header_match:
                idx = int(invalid_header_match.group(1)) - 1
                if 0 <= idx < len(headers):
                    problem_indexes.add(idx)
            outline_match = re.search(r"Outline header '([^']+)'", text)
            if outline_match:
                header_name = outline_match.group(1).strip()
                idx = header_map.get(header_name)
                if idx is not None:
                    problem_indexes.add(idx)
            missing_match = re.search(r"Missing Comparison columns present in protein file:\s*(.+)\.", text)
            if missing_match:
                for part in missing_match.group(1).split(","):
                    idx = header_map.get(part.strip())
                    if idx is not None:
                        problem_indexes.add(idx)
        return sorted(problem_indexes)

    def _build_input_preview_table(
        headers: List[str],
        row_values: List[str],
        column_widths: List[str],
        error_columns: Optional[List[int]] = None,
        fixed_width: Optional[str] = None,
    ) -> Any:
        error_set = set(error_columns or [])
        wrap_style = ""
        table_style = ""
        if fixed_width:
            wrap_style = f"display:block; width:{fixed_width}; max-width:none;"
            table_style = "width:100%;"
        colgroup = ui.tags.colgroup(
            *[
                ui.tags.col(style=f"width:{column_widths[idx]}; min-width:{column_widths[idx]};")
                for idx in range(len(headers))
            ]
        )
        return ui.div(
            {"class": "input-preview-wrap", "style": wrap_style} if wrap_style else {"class": "input-preview-wrap"},
            ui.tags.table(
                {"class": "table table-sm input-preview-table", "style": table_style} if table_style else {"class": "table table-sm input-preview-table"},
                colgroup,
                ui.tags.thead(
                    ui.tags.tr(
                        *[
                            ui.tags.th({"class": "input-preview-header-error"} if idx in error_set else {}, h)
                            for idx, h in enumerate(headers)
                        ]
                    )
                ),
                ui.tags.tbody(ui.tags.tr(*[ui.tags.td(v) for v in row_values])),
            ),
        )

    def _build_input_preview_guide() -> Any:
        return ui.div(
            {"class": "input-preview-guide-card"},
            ui.div(
                {"class": "input-preview-guide-layout"},
                ui.div(
                    {"class": "input-preview-guide-body"},
                    ui.div(
                        {"class": "input-preview-guide-pill"},
                        ui.div(
                            {"class": "input-preview-guide-pill-title"},
                            "UniProt ID",
                            ui.tags.span({"class": "input-preview-guide-required"}, "*"),
                        ),
                        ui.div(
                            {"class": "input-preview-guide-pill-text"},
                            "Put the protein identifier for each Protein or PTM row in the first column. The header name can be anything.",
                        ),
                    ),
                    ui.div(
                        {"class": "input-preview-guide-pill"},
                        ui.div(
                            {"class": "input-preview-guide-pill-title"},
                            "Gene Symbol",
                            ui.tags.span({"class": "input-preview-guide-required"}, "*"),
                        ),
                        ui.div(
                            {"class": "input-preview-guide-pill-text"},
                            "Put the matching gene symbol or display label in the second column of the Protein file. The header name can be anything.",
                        ),
                    ),
                    ui.div(
                        {"class": "input-preview-guide-pill"},
                        ui.div(
                            {"class": "input-preview-guide-pill-title"},
                            "PTM Site",
                            ui.tags.span({"class": "input-preview-guide-required"}, "*"),
                        ),
                        ui.div(
                            {"class": "input-preview-guide-pill-text"},
                            "Put the modified residue position for each PTM entry in the second column of the PTM file. The header name can be anything.",
                        ),
                    ),
                    ui.div(
                        {"class": "input-preview-guide-pill"},
                        ui.div(
                            {"class": "input-preview-guide-pill-title"},
                            "Comparisons",
                            ui.tags.span({"class": "input-preview-guide-required"}, "*"),
                        ),
                        ui.div(
                            {"class": "input-preview-guide-pill-text"},
                            "Use log2 fold-change values comparing two conditions. Each comparison header must start with 'C:' and comparison headers must match between datasets.",
                        ),
                    ),
                    ui.div(
                        {"class": "input-preview-guide-pill"},
                        ui.div(
                            {"class": "input-preview-guide-pill-title"},
                            "HMDB ID",
                            ui.tags.span({"class": "input-preview-guide-required"}, "*"),
                        ),
                        ui.div(
                            {"class": "input-preview-guide-pill-text"},
                            "Put the HMDB identifier for each metabolite in the first column of the metabolite file. The header name can be anything.",
                        ),
                    ),
                    ui.div(
                        {"class": "input-preview-guide-pill"},
                        ui.div({"class": "input-preview-guide-pill-title"}, "Tooltip"),
                        ui.div(
                            {"class": "input-preview-guide-pill-text"},
                            "Optional tooltip or descriptive columns. These headers should start with 'T:'.",
                        ),
                    ),
                    ui.div(
                        {"class": "input-preview-guide-pill"},
                        ui.div({"class": "input-preview-guide-pill-title"}, "Dual-Comparison"),
                        ui.div(
                            {"class": "input-preview-guide-pill-text"},
                            "Optional outline columns. These headers should start with 'O:' and correspond to a comparison header. Quantitative values are displayed on the outline of nodes.",
                        ),
                    ),
                ),
            ),
        )

    def _write_debug_dump(filename: str, payload: Dict[str, Any]) -> None:
        if not DEBUG_FILE_OUTPUT_ENABLED:
            return
        try:
            debug_path = safe_session_path(session, f"debug/{Path(filename).name}")
            with open(debug_path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
        except Exception as exc:
            print(f"Failed to write debug dump {filename}: {exc}")

    def _sync_catalog_into_open_payloads(catalog_info: Dict[str, Any]) -> None:
        for state in bookmark_state.values():
            current = state["json"].get()
            if not current:
                continue
            updated = dict(current)
            updated["_global_protein_catalog"] = dict(catalog_info)
            updated["_persist_token"] = time.time()
            state["json"].set(updated)

    def _reset_global_catalog_from_default() -> None:
        info = dict(GLOBAL_CATALOG_INFO)
        global_catalog_info.set(info)
        _sync_catalog_into_open_payloads(info)

    def _global_catalog_build_inputs(
        species_selection: Optional[str] = None,
    ) -> Optional[Tuple[Dict[str, Any], Dict[str, Any]]]:
        data_override = collect_data_override()
        if not data_override:
            return None
        resolved_species_selection = (
            species_selection
            if species_selection is not None
            else _get_input_value(input, "input_species")
        )
        species_choice, species_info = _resolve_species(resolved_species_selection)
        protein_cfg = data_override.get("protein", {}) if isinstance(data_override, dict) else {}
        settings_override = {
            "species": species_choice,
            "species_code": species_info.get("code", ""),
            "prot_uniprot_column": protein_cfg.get("uniprot_column", "Uniprot_ID"),
            "gene_name_column": protein_cfg.get("gene_column", "Gene Symbol"),
            "main_columns": protein_cfg.get("main_columns", []),
            "outline_columns": protein_cfg.get("outline_columns", []),
            "protein_tooltip_columns": protein_cfg.get("tooltip_columns", []),
            "hsa_id_column": protein_cfg.get("kegg_column", protein_cfg.get("uniprot_column", "Uniprot_ID")),
        }
        return data_override, settings_override

    def _apply_global_catalog_payload(payload: Dict[str, Any]) -> None:
        info = {
            "path": "",
            "metadata": payload.get("metadata", {}),
            "protein_catalog": payload.get("protein_catalog", {}),
            "scope": "session",
        }
        global_catalog_info.set(dict(info))
        _sync_catalog_into_open_payloads(dict(info))

    async def _refresh_global_catalog_from_current_async(
        species_selection: Optional[str] = None,
    ) -> None:
        build_inputs = _global_catalog_build_inputs(species_selection=species_selection)
        if build_inputs is None:
            _reset_global_catalog_from_default()
            return
        data_override, settings_override = build_inputs
        try:
            payload = await asyncio.to_thread(
                build_global_protein_catalog,
                data_override=data_override,
                settings_override=settings_override,
            )
        except Exception as exc:
            print(f"Warning: failed to refresh global protein catalog: {exc}")
            return
        _apply_global_catalog_payload(payload)

    def _current_global_catalog_info() -> Dict[str, Any]:
        try:
            with reactive.isolate():
                info = global_catalog_info.get()
        except RuntimeError:
            info = None
        info_dict = dict(info or GLOBAL_CATALOG_INFO)
        if info_dict.get("protein_catalog"):
            return info_dict
        path = info_dict.get("path")
        if path and os.path.exists(path):
            try:
                with open(path, "r", encoding="utf-8") as fh:
                    payload = json.load(fh)
                protein_catalog = payload.get("protein_catalog")
                if isinstance(protein_catalog, dict):
                    info_dict["protein_catalog"] = protein_catalog
            except Exception as exc:
                print(f"Warning: failed to load global protein catalog from {path}: {exc}")
        return info_dict

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

    def _clear_input_busy_state(message: str = "") -> None:
        _send_custom_message_safe(
            "set_input_busy_state",
            {
                "active": False,
                "title": "Processing input",
                "text": message or "Validation finished. Review the input status before running again.",
            },
        )

    def _apply_outline_width_defaults(protein_headers: Optional[Sequence[str]] = None, ptm_headers: Optional[Sequence[str]] = None) -> None:
        try:
            if protein_headers is not None:
                session.send_input_message(
                    "settings_prot_outline_width",
                    {"value": 3 if _has_outline_columns(protein_headers) else 1},
                )
            if ptm_headers is not None:
                session.send_input_message(
                    "settings_ptm_outline_width",
                    {"value": 3 if _has_outline_columns(ptm_headers) else 1},
                )
        except Exception as exc:
            print(f"Warning: failed to apply outline width defaults: {exc}")

    async def _load_demo_datasets():
        """Load bundled demo protein/PTM files so Demo mode behaves like preloaded uploads."""
        missing = []
        if not os.path.exists(SAMPLE_PROTEIN_FILE):
            missing.append(f"Protein file not found: {SAMPLE_PROTEIN_FILE}")
        if not os.path.exists(SAMPLE_PTM_FILE):
            missing.append(f"PTM file not found: {SAMPLE_PTM_FILE}")
        if missing:
            protein_validation.set({"status": "Demo mode files missing.", "errors": missing, "valid": False, "comparisons": []})
            ptm_validation.set({"status": "Demo mode files missing.", "errors": missing, "valid": False})
            validated_protein_dataset.set(None)
            validated_ptm_dataset.set(None)
            validated_metabolite_dataset.set(None)
            protein_dataset.set(None)
            ptm_dataset.set(None)
            metabolite_dataset.set(None)
            protein_preview_dataset.set(None)
            ptm_preview_dataset.set(None)
            metabolite_preview_dataset.set(None)
            _reset_global_catalog_from_default()
            pipeline_running.set(False)
            pipeline_ready.set(False)
            nav_lock_status.set("Demo mode: sample files missing. Navigation locked.")
            _clear_pathway_scores("Demo datasets missing; pathway scoring unavailable.")
            return

        try:
            prot_result = validate_protein_file(SAMPLE_PROTEIN_FILE)
            if not prot_result.valid:
                protein_validation.set({"status": "Demo protein file failed validation.", "errors": prot_result.errors, "valid": False, "comparisons": []})
                ptm_validation.set({"status": "PTM upload disabled until demo protein is valid.", "errors": [], "valid": False})
                metabolite_validation.set({"status": "Metabolite upload optional. Sample metabolite file available below.", "errors": [], "valid": False, "comparisons": []})
                validated_protein_dataset.set(None)
                validated_ptm_dataset.set(None)
                validated_metabolite_dataset.set(None)
                protein_dataset.set(None)
                ptm_dataset.set(None)
                metabolite_dataset.set(None)
                protein_preview_dataset.set(None)
                ptm_preview_dataset.set(None)
                metabolite_preview_dataset.set(None)
                _reset_global_catalog_from_default()
                pipeline_running.set(False)
                pipeline_ready.set(False)
                nav_lock_status.set("Demo mode: sample protein validation failed.")
                _clear_pathway_scores("Demo protein validation failed; pathway scoring unavailable.")
                return

            species_choice = DEMO_SPECIES
            species_info = SPECIES_CHOICES.get(DEMO_SPECIES, SPECIES_CHOICES.get(DEFAULT_SPECIES, {}))
            current_species = str(_get_input_value(input, "input_species") or "").strip().lower()
            if species_choice and current_species != species_choice:
                demo_species_sync_target.set(species_choice)
                try:
                    session.send_input_message("input_species", {"value": species_choice})
                except Exception:
                    demo_species_sync_target.set(None)
                    raise
            species_code = species_info.get("code", "") or SPECIES_CHOICES.get(DEFAULT_SPECIES, {}).get("code", "hsa")
            print(f"Demo mode: loading sample datasets for species '{species_choice}' code '{species_code}'")

            demo_prot_payload = _fill_blank_gene_with_uniprot(_load_dataset(SAMPLE_PROTEIN_FILE))
            protein_preview_dataset.set(demo_prot_payload)
            _apply_outline_width_defaults(protein_headers=demo_prot_payload.get("headers") or [])
            try:
                demo_prot_payload = annotate_protein_with_kegg(demo_prot_payload, species_code, _get_kegg_map(species_code))
            except Exception as exc:
                print(f"Warning: demo protein KEGG annotation failed: {exc}")
            try:
                demo_prot_payload = _annotate_dataset_with_human_orthologs(demo_prot_payload, species_code)
            except Exception as exc:
                print(f"Warning: demo protein human ortholog annotation failed: {exc}")
            protein_dataset.set(demo_prot_payload)
            protein_dataset_path.set(os.path.basename(SAMPLE_PROTEIN_FILE))
            protein_upload_size_bytes.set(_upload_file_size_bytes(SAMPLE_PROTEIN_FILE))
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
                validated_protein_dataset.set(None)
                validated_ptm_dataset.set(None)
                validated_metabolite_dataset.set(None)
                ptm_dataset.set(None)
                ptm_preview_dataset.set(None)
                metabolite_dataset.set(None)
                metabolite_preview_dataset.set(None)
                metabolite_dataset_path.set(None)
                metabolite_upload_size_bytes.set(0)
                metabolite_validation.set({"status": "Metabolite upload optional. Sample metabolite file available below.", "errors": [], "valid": False, "comparisons": []})
                pipeline_running.set(False)
                pipeline_ready.set(False)
                nav_lock_status.set("Demo mode: sample PTM validation failed.")
                _refresh_pathway_scores()
                return

            demo_ptm_payload = _load_dataset(SAMPLE_PTM_FILE)
            ptm_preview_dataset.set(demo_ptm_payload)
            _apply_outline_width_defaults(ptm_headers=demo_ptm_payload.get("headers") or [])
            try:
                demo_ptm_payload = annotate_ptm_dataset(demo_ptm_payload, species_choice, _get_psp_map())
                demo_ptm_payload = annotate_ptm_dataset_with_kinases(demo_ptm_payload, species_choice, _get_psp_ks_map())
            except Exception as exc:
                print(f"Warning: demo PTM annotation failed: {exc}")
            try:
                demo_ptm_payload = _annotate_dataset_with_human_orthologs(demo_ptm_payload, species_code)
            except Exception as exc:
                print(f"Warning: demo PTM human ortholog annotation failed: {exc}")
            ptm_dataset.set(demo_ptm_payload)
            ptm_dataset_path.set(os.path.basename(SAMPLE_PTM_FILE))
            ptm_upload_size_bytes.set(_upload_file_size_bytes(SAMPLE_PTM_FILE))
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
            metabolite_dataset.set(None)
            metabolite_preview_dataset.set(None)
            metabolite_dataset_path.set(None)
            metabolite_upload_size_bytes.set(0)
            metabolite_validation.set({
                "status": "Metabolite upload optional. Sample metabolite file available below.",
                "errors": [],
                "valid": False,
                "comparisons": [],
            })
            await _refresh_global_catalog_from_current_async(species_selection=species_choice)
            _update_ks_index()
            _refresh_pathway_scores(species_selection=species_choice)
            pipeline_running.set(False)
            pipeline_ready.set(True)
            nav_lock_status.set("Demo mode: sample datasets loaded. Navigation unlocked.")
            print("Demo mode: sample protein/PTM loaded successfully.")
        except Exception as exc:
            print(f"Warning: demo mode encountered an unexpected error: {exc}")
            protein_validation.set({"status": "Demo mode failed to load sample protein.", "errors": [str(exc)], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "Demo mode failed to load sample PTM.", "errors": [str(exc)], "valid": False})
            metabolite_validation.set({"status": "Metabolite upload optional. Sample metabolite file available below.", "errors": [], "valid": False, "comparisons": []})
            protein_dataset.set(None)
            ptm_dataset.set(None)
            metabolite_dataset.set(None)
            protein_preview_dataset.set(None)
            ptm_preview_dataset.set(None)
            metabolite_preview_dataset.set(None)
            _reset_global_catalog_from_default()
            pipeline_running.set(False)
            pipeline_ready.set(False)
            nav_lock_status.set("Demo mode: sample datasets failed to load.")
            _clear_pathway_scores("Demo datasets failed to load; pathway scoring unavailable.")

    def collect_data_override() -> Dict[str, Any]:  # type: ignore[override]
        mode = _get_input_value(input, "input_mode") or "user"
        if str(mode).lower() not in {"user", "demo"}:
            return {}
        with reactive.isolate():
            protein_valid = protein_validation.get().get("valid")
            ptm_valid = ptm_validation.get().get("valid")
            metabolite_valid = metabolite_validation.get().get("valid")
            protein_data = protein_dataset.get()
            ptm_data = ptm_dataset.get()
            metabolite_data = metabolite_dataset.get()
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
                "file_path": _safe_source_label(protein_dataset_path.get() or ""),
                "uniprot_column": prot_uniprot,
                "kegg_column": prot_kegg,
                "gene_column": prot_gene,
                "main_columns": prot_main,
                "outline_columns": prot_outline,
                "tooltip_columns": prot_tooltips,
            },
            "ptm": ptm_entries,
        }
        if metabolite_valid and metabolite_data:
            metabolite_headers = list(metabolite_data.get("headers") or [])
            if metabolite_headers:
                metabolite_hmdb = metabolite_headers[0]
                data_override["metabolite"] = {
                    "data_headers": metabolite_headers,
                    "data_rows": metabolite_data.get("rows") or [],
                    "file_path": _safe_source_label(metabolite_dataset_path.get() or ""),
                    "hmdb_column": metabolite_hmdb,
                    "main_columns": [h for h in metabolite_headers if h.startswith("C:")],
                    "tooltip_columns": [h for h in metabolite_headers if h.startswith("T:")],
                }
        # Attach file_path to each PTM entry so downstream processors/catalog use the correct dataset
        for entry in data_override["ptm"]:
            entry["file_path"] = _safe_source_label(ptm_dataset_path.get() or "")
        return data_override

    @reactive.Effect
    @reactive.event(input.input_protein_upload)
    def _process_protein_upload():
        if str((_get_input_value(input, "input_mode") or "user")).lower() == "demo":
            protein_validation.set({"status": "Demo mode uses bundled sample files; switch to User mode to upload.", "errors": [], "valid": False, "comparisons": []})
            validated_protein_dataset.set(None)
            protein_kegg_warning.set("")
            _clear_pathway_scores("Demo mode active in User upload handler; pathway scoring deferred.")
            return
        upload = input.input_protein_upload()
        if not upload:
            protein_validation.set({"status": "Upload a protein file to begin.", "errors": [], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "PTM upload disabled until a valid protein file is uploaded.", "errors": [], "valid": False})
            metabolite_validation.set({"status": "Metabolite upload optional. Provide after protein if available.", "errors": [], "valid": False, "comparisons": []})
            validated_protein_dataset.set(None)
            validated_ptm_dataset.set(None)
            validated_metabolite_dataset.set(None)
            protein_preview_dataset.set(None)
            ptm_preview_dataset.set(None)
            metabolite_preview_dataset.set(None)
            protein_dataset_path.set(None)
            ptm_dataset_path.set(None)
            metabolite_dataset_path.set(None)
            protein_upload_size_bytes.set(0)
            ptm_upload_size_bytes.set(0)
            metabolite_upload_size_bytes.set(0)
            _reset_global_catalog_from_default()
            _write_debug_dump("user_protein_dataset_debug.txt", {"info": "No dataset loaded"})
            protein_kegg_warning.set("")
            _update_ks_index(reset=True)
            pipeline_ready.set(False)
            pipeline_running.set(False)
            _clear_pathway_scores("Pathway scoring waiting for validated inputs and Run.")
            return
        file_info = upload[0]
        upload_temp_path = _extract_upload_datapath(file_info)
        protein_validation.set({"status": "Validating protein upload...", "errors": [], "valid": False, "comparisons": []})
        try:
            upload_validation = validate_uploaded_file(file_info, expected_role="protein", session_id=_safe_session_id())
        except Exception:
            _cleanup_upload_temp_file(upload_temp_path, "protein upload validation exception")
            msg = "Invalid analysis dataset upload: file could not be validated."
            protein_validation.set({"status": msg, "errors": [msg], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
            metabolite_validation.set({"status": "Metabolite upload optional. Provide after protein if available.", "errors": [], "valid": False, "comparisons": []})
            validated_protein_dataset.set(None)
            validated_ptm_dataset.set(None)
            validated_metabolite_dataset.set(None)
            protein_preview_dataset.set(None)
            ptm_preview_dataset.set(None)
            metabolite_preview_dataset.set(None)
            protein_dataset_path.set(None)
            ptm_dataset_path.set(None)
            metabolite_dataset_path.set(None)
            protein_upload_size_bytes.set(0)
            ptm_upload_size_bytes.set(0)
            metabolite_upload_size_bytes.set(0)
            _reset_global_catalog_from_default()
            _write_debug_dump("user_protein_dataset_debug.txt", {"info": "No dataset loaded"})
            protein_kegg_warning.set("")
            _update_ks_index(reset=True)
            pipeline_ready.set(False)
            pipeline_running.set(False)
            _clear_pathway_scores("Pathway scoring waiting for validated inputs and Run.")
            return
        if not upload_validation.valid:
            _cleanup_upload_temp_file(upload_validation.datapath, "invalid protein upload validation")
            protein_validation.set({"status": upload_validation.user_message, "errors": [upload_validation.user_message], "valid": False, "comparisons": []})
            ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
            metabolite_validation.set({"status": "Metabolite upload optional. Provide after protein if available.", "errors": [], "valid": False, "comparisons": []})
            validated_protein_dataset.set(None)
            validated_ptm_dataset.set(None)
            validated_metabolite_dataset.set(None)
            protein_preview_dataset.set(None)
            ptm_preview_dataset.set(None)
            metabolite_preview_dataset.set(None)
            protein_dataset_path.set(None)
            ptm_dataset_path.set(None)
            metabolite_dataset_path.set(None)
            protein_upload_size_bytes.set(0)
            ptm_upload_size_bytes.set(0)
            metabolite_upload_size_bytes.set(0)
            _reset_global_catalog_from_default()
            _write_debug_dump("user_protein_dataset_debug.txt", {"info": "No dataset loaded"})
            protein_kegg_warning.set("")
            _update_ks_index(reset=True)
            pipeline_ready.set(False)
            pipeline_running.set(False)
            _clear_pathway_scores("Pathway scoring waiting for validated inputs and Run.")
            return
        datapath = upload_validation.datapath
        upload_size_bytes = _upload_file_size_bytes(datapath)
        try:
            result = validate_protein_file(datapath)
            if not result.valid:
                protein_preview_dataset.set(_fill_blank_gene_with_uniprot(_load_dataset(datapath)))
                protein_validation.set({"status": "Protein file failed validation. See errors below.", "errors": result.errors, "valid": False, "comparisons": []})
                ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
                metabolite_validation.set({"status": "Metabolite upload optional. Provide after protein if available.", "errors": [], "valid": False, "comparisons": []})
                validated_protein_dataset.set(None)
                validated_ptm_dataset.set(None)
                validated_metabolite_dataset.set(None)
                protein_dataset_path.set(None)
                ptm_dataset_path.set(None)
                metabolite_dataset_path.set(None)
                protein_upload_size_bytes.set(0)
                ptm_upload_size_bytes.set(0)
                metabolite_upload_size_bytes.set(0)
                ptm_preview_dataset.set(None)
                metabolite_preview_dataset.set(None)
                _reset_global_catalog_from_default()
                _write_debug_dump("user_protein_dataset_debug.txt", {"errors": result.errors})
                protein_kegg_warning.set("")
                _update_ks_index(reset=True)
                pipeline_ready.set(False)
                pipeline_running.set(False)
                _clear_pathway_scores("Pathway scoring waiting for validated inputs and Run.")
                return

            status = (
                f"Protein file valid. Rows: {result.summary.get('rows', 0)}, "
                f"Comparison columns: {result.summary.get('comparisons', 0)}, "
                f"Tooltip columns: {result.summary.get('tooltips', 0)}."
            )
            protein_validation.set({"status": status, "errors": [], "valid": True, "comparisons": result.comparisons})
            dataset_payload = _fill_blank_gene_with_uniprot(_load_dataset(datapath))
            protein_preview_dataset.set(dataset_payload)
            _apply_outline_width_defaults(protein_headers=dataset_payload.get("headers") or [])
            validated_protein_dataset.set(dataset_payload)
            protein_kegg_warning.set("")
            protein_dataset_path.set(upload_validation.sanitized_filename)
            protein_upload_size_bytes.set(upload_size_bytes)
            ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
            metabolite_validation.set({"status": "Metabolite upload optional. Provide after protein if available.", "errors": [], "valid": False, "comparisons": []})
            validated_ptm_dataset.set(None)
            validated_metabolite_dataset.set(None)
            ptm_dataset_path.set(None)
            metabolite_dataset_path.set(None)
            ptm_upload_size_bytes.set(0)
            metabolite_upload_size_bytes.set(0)
            ptm_preview_dataset.set(None)
            metabolite_preview_dataset.set(None)
            _invalidate_user_run("Protein validated. Upload optional PTM/Metabolite files, then press Run.")
        finally:
            _cleanup_upload_temp_file(datapath, "protein upload parsed")

    @output
    @render.text
    def input_protein_status():
        data = protein_validation.get()
        return data.get("status", "")

    @output
    @render.ui
    def input_protein_preview_guide():
        return _build_input_preview_guide()

    @output
    @render.ui
    def input_protein_preview():
        protein_preview_data = protein_preview_dataset.get()
        protein_headers, protein_row = _get_input_preview_content(
            protein_preview_data,
            ["Uniprot", "GeneSymbol", "C: Comparison (1)", "...", "T: Tooltip(1)", "...", "O: Outline(1)", "..."],
        )
        ptm_headers, ptm_row = _get_input_preview_content(
            ptm_preview_dataset.get(),
            ["Uniprot", "PTM Site", "C: Comparison (1)", "...", "T: Tooltip(1)", "...", "O: Outline(1)", "..."],
        )
        metabolite_headers, metabolite_row = _get_input_preview_content(
            metabolite_preview_dataset.get(),
            ["HMDB ID", "C: Comparison (1)", "...", "T: Tooltip(1)", "..."],
        )
        return ui.div(
            {"class": "input-preview-labeled-row"},
            ui.div({"class": "input-preview-side-label is-protein"}, "Protein"),
            ui.div(
                {"class": "input-preview-table-wrap"},
                _build_input_preview_table(
                    protein_headers,
                    protein_row,
                    _get_input_preview_widths((protein_headers, protein_row), (ptm_headers, ptm_row), (metabolite_headers, metabolite_row)),
                    _extract_problem_column_indexes(protein_headers, protein_validation.get().get("errors") or []),
                    None if protein_preview_data else "1104px",
                ),
            ),
        )

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
    def input_ptm_preview():
        ptm_preview_data = ptm_preview_dataset.get()
        protein_headers, protein_row = _get_input_preview_content(
            protein_preview_dataset.get(),
            ["Uniprot", "GeneSymbol", "C: Comparison (1)", "...", "T: Tooltip(1)", "...", "O: Outline(1)", "..."],
        )
        ptm_headers, ptm_row = _get_input_preview_content(
            ptm_preview_data,
            ["Uniprot", "PTM Site", "C: Comparison (1)", "...", "T: Tooltip(1)", "...", "O: Outline(1)", "..."],
        )
        metabolite_headers, metabolite_row = _get_input_preview_content(
            metabolite_preview_dataset.get(),
            ["HMDB ID", "C: Comparison (1)", "...", "T: Tooltip(1)", "..."],
        )
        return ui.div(
            {"class": "input-preview-labeled-row"},
            ui.div({"class": "input-preview-side-label is-ptm"}, "PTM"),
            ui.div(
                {"class": "input-preview-table-wrap"},
                _build_input_preview_table(
                    ptm_headers,
                    ptm_row,
                    _get_input_preview_widths((protein_headers, protein_row), (ptm_headers, ptm_row), (metabolite_headers, metabolite_row)),
                    _extract_problem_column_indexes(ptm_headers, ptm_validation.get().get("errors") or []),
                    None if ptm_preview_data else "1104px",
                ),
            ),
        )

    @output
    @render.ui
    def input_metabolite_preview():
        protein_headers, protein_row = _get_input_preview_content(
            protein_preview_dataset.get(),
            ["Uniprot", "GeneSymbol", "C: Comparison (1)", "...", "T: Tooltip(1)", "...", "O: Outline(1)", "..."],
        )
        ptm_headers, ptm_row = _get_input_preview_content(
            ptm_preview_dataset.get(),
            ["Uniprot", "PTM Site", "C: Comparison (1)", "...", "T: Tooltip(1)", "...", "O: Outline(1)", "..."],
        )
        metabolite_preview_data = metabolite_preview_dataset.get()
        metabolite_headers, metabolite_row = _get_input_preview_content(
            metabolite_preview_data,
            ["HMDB ID", "C: Comparison (1)", "...", "T: Tooltip(1)", "..."],
        )
        metabolite_headers = [
            metabolite_headers[0] if len(metabolite_headers) > 0 else "HMDB ID",
            "",
            metabolite_headers[1] if len(metabolite_headers) > 1 else "C: Comparison (1)",
            metabolite_headers[2] if len(metabolite_headers) > 2 else "...",
            metabolite_headers[3] if len(metabolite_headers) > 3 else "T: Tooltip(1)",
            metabolite_headers[4] if len(metabolite_headers) > 4 else "...",
            "",
            "",
        ]
        metabolite_row = [
            metabolite_row[0] if len(metabolite_row) > 0 else "",
            "",
            metabolite_row[1] if len(metabolite_row) > 1 else "",
            metabolite_row[2] if len(metabolite_row) > 2 else "",
            metabolite_row[3] if len(metabolite_row) > 3 else "",
            metabolite_row[4] if len(metabolite_row) > 4 else "",
            "",
            "",
        ]
        return ui.div(
            {"class": "input-preview-labeled-row"},
            ui.div({"class": "input-preview-side-label is-metabolite"}, "Metabolite"),
            ui.div(
                {"class": "input-preview-table-wrap"},
                _build_input_preview_table(
                    metabolite_headers,
                    metabolite_row,
                    _get_input_preview_widths((protein_headers, protein_row), (ptm_headers, ptm_row), (metabolite_headers, metabolite_row)),
                    None,
                    None if metabolite_preview_data else "1104px",
                ),
            ),
        )

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
    @render.ui
    def input_metabolite_errors():
        data = metabolite_validation.get()
        errs = data.get("errors") or []
        if not errs:
            return ui.div()
        items = [ui.tags.li(err) for err in errs[:50]]
        if len(errs) > 50:
            items.append(ui.tags.li(f"...and {len(errs) - 50} more errors."))
        return ui.tags.ul(*items, style="color:#b00020;")

    @reactive.Effect
    @reactive.event(input.input_metabolite_upload)
    def _process_metabolite_upload():
        if not protein_validation.get().get("valid"):
            metabolite_validation.set({"status": "Upload a valid protein file first.", "errors": ["Protein file not validated yet."], "valid": False, "comparisons": []})
            validated_metabolite_dataset.set(None)
            metabolite_preview_dataset.set(None)
            metabolite_dataset_path.set(None)
            metabolite_upload_size_bytes.set(0)
            _invalidate_user_run("Upload and validate a Protein dataset before running.")
            return
        if str((_get_input_value(input, "input_mode") or "user")).lower() == "demo":
            metabolite_validation.set({"status": "Demo mode uses bundled sample files; switch to User mode to upload.", "errors": [], "valid": False, "comparisons": []})
            validated_metabolite_dataset.set(None)
            return
        upload = input.input_metabolite_upload()
        if not upload:
            metabolite_validation.set({"status": "No metabolite file uploaded (optional).", "errors": [], "valid": False, "comparisons": []})
            validated_metabolite_dataset.set(None)
            metabolite_preview_dataset.set(None)
            metabolite_dataset_path.set(None)
            metabolite_upload_size_bytes.set(0)
            _invalidate_user_run("Validation updated. Press Run to process datasets.")
            return
        file_info = upload[0]
        upload_temp_path = _extract_upload_datapath(file_info)
        try:
            upload_validation = validate_uploaded_file(file_info, expected_role="metabolite", session_id=_safe_session_id())
        except Exception:
            _cleanup_upload_temp_file(upload_temp_path, "metabolite upload validation exception")
            msg = "Invalid analysis dataset upload: file could not be validated."
            metabolite_validation.set({"status": msg, "errors": [msg], "valid": False, "comparisons": []})
            validated_metabolite_dataset.set(None)
            metabolite_preview_dataset.set(None)
            metabolite_dataset_path.set(None)
            metabolite_upload_size_bytes.set(0)
            _invalidate_user_run("Fix the Metabolite dataset issue, then press Run.")
            return
        if not upload_validation.valid:
            _cleanup_upload_temp_file(upload_validation.datapath, "invalid metabolite upload validation")
            metabolite_validation.set({"status": upload_validation.user_message, "errors": [upload_validation.user_message], "valid": False, "comparisons": []})
            validated_metabolite_dataset.set(None)
            metabolite_preview_dataset.set(None)
            metabolite_dataset_path.set(None)
            metabolite_upload_size_bytes.set(0)
            _invalidate_user_run("Fix the Metabolite dataset issue, then press Run.")
            return
        datapath = upload_validation.datapath
        upload_size_bytes = _upload_file_size_bytes(datapath)
        try:
            protein_comparisons = protein_validation.get().get("comparisons") or []
            result = _validate_metabolite_file(datapath, protein_comparisons)
            preview_payload = _load_dataset(datapath)
            metabolite_preview_dataset.set(preview_payload)
            if result.valid:
                validated_metabolite_dataset.set(preview_payload)
                metabolite_dataset_path.set(upload_validation.sanitized_filename)
                metabolite_upload_size_bytes.set(upload_size_bytes)
                metabolite_validation.set({
                    "status": (
                        f"Metabolite file valid. Rows: {result.summary.get('rows', 0)}, "
                        f"Comparison columns: {result.summary.get('comparisons', 0)}, "
                        f"Tooltip columns: {result.summary.get('tooltips', 0)}."
                    ),
                    "errors": [],
                    "valid": True,
                    "comparisons": result.comparisons,
                })
                _invalidate_user_run("Metabolite validated. Press Run to process datasets.")
            else:
                validated_metabolite_dataset.set(None)
                metabolite_dataset_path.set(None)
                metabolite_upload_size_bytes.set(0)
                metabolite_validation.set({
                    "status": "Metabolite file failed validation. See errors below.",
                    "errors": result.errors,
                    "valid": False,
                    "comparisons": [],
                })
                _invalidate_user_run("Fix the Metabolite dataset errors before running.")
        finally:
            _cleanup_upload_temp_file(datapath, "metabolite upload parsed")

    @output
    @render.text
    def input_nav_lock_status():
        return nav_lock_status.get()

    @output
    @render.ui
    def input_run_button_state():
        mode = _current_mode()
        running = bool(pipeline_running.get())
        protein_valid = bool(protein_validation.get().get("valid"))
        protein_error = _validation_has_explicit_error(protein_validation.get())
        ptm_error = _validation_has_explicit_error(ptm_validation.get())
        metabolite_error = _validation_has_explicit_error(metabolite_validation.get())
        disabled = mode == "demo" or running or (not protein_valid) or protein_error or ptm_error or metabolite_error
        run_guard_active = bool(str(run_guard_message.get() or "").strip())
        return ui.tags.script(
            f"""
            (function(){{
                const btn = document.getElementById('input_run_pipeline');
                if (btn) btn.disabled = {str(disabled).lower()};
                const status = document.querySelector('.input-run-status');
                if (status) status.classList.toggle('is-error', {str(run_guard_active).lower()});
            }})();
            """
        )

    @reactive.Effect
    @reactive.event(input.input_run_pipeline)
    async def _run_validated_input_pipeline():
        if _current_mode() == "demo":
            _clear_input_busy_state("Demo mode uses bundled sample data.")
            return
        if pipeline_running.get():
            return
        run_guard_message.set("")
        if not bool(protein_validation.get().get("valid")):
            msg = "Upload and validate a Protein dataset before running."
            _invalidate_user_run(msg)
            _clear_input_busy_state(msg)
            return
        if _validation_has_explicit_error(protein_validation.get()):
            msg = "Fix the Protein dataset errors before running."
            _invalidate_user_run(msg)
            _clear_input_busy_state(msg)
            return
        if _validation_has_explicit_error(ptm_validation.get()):
            msg = "Fix the PTM dataset errors before running."
            _invalidate_user_run(msg)
            _clear_input_busy_state(msg)
            return
        if _validation_has_explicit_error(metabolite_validation.get()):
            msg = "Fix the Metabolite dataset errors before running."
            _invalidate_user_run(msg)
            _clear_input_busy_state(msg)
            return
        raw_protein = validated_protein_dataset.get() or {}
        raw_ptm = validated_ptm_dataset.get() or {}
        raw_metabolite = validated_metabolite_dataset.get() or {}
        if not raw_protein or not raw_protein.get("headers"):
            msg = "Upload and validate a Protein dataset before running."
            _invalidate_user_run(msg)
            _clear_input_busy_state(msg)
            return

        # H) Apply per-session analysis-run throttling before expensive processing.
        # This guard tracks attempted runs (not only accepted runs) so repeated rapid clicks are limited.
        now_monotonic = time.monotonic()
        allow_run, throttle_message, updated_attempts = evaluate_run_throttle(
            now_monotonic=now_monotonic,
            last_accepted_run_monotonic=last_accepted_run_monotonic.get(),
            prior_attempt_times_monotonic=run_attempt_times_monotonic.get() or [],
            min_seconds_between_runs=MIN_SECONDS_BETWEEN_RUNS,
            max_runs_per_minute=MAX_RUNS_PER_MINUTE,
        )
        run_attempt_times_monotonic.set(updated_attempts)
        if not allow_run:
            run_guard_message.set(throttle_message)
            nav_lock_status.set(throttle_message)
            _clear_input_busy_state(throttle_message)
            log_run_guard_event(
                session_id=_safe_session_id(),
                reason=throttle_message,
                run_attempt_count=len(updated_attempts),
            )
            return
        run_guard_message.set("")

        # G) Validate combined upload size before expensive processing.
        upload_sizes = _current_run_upload_sizes()
        _size_ok, total_upload_size = validate_total_upload_size_from_sizes([size for _role, size in upload_sizes])
        if total_upload_size is None:
            msg = "Upload rejected: one or more uploaded files could not be read."
            run_guard_message.set(msg)
            nav_lock_status.set(msg)
            _clear_input_busy_state(msg)
            log_run_guard_event(session_id=_safe_session_id(), reason="uploaded file size missing during run-size validation")
            return
        if not _size_ok:
            msg = f"Upload rejected: combined upload size cannot exceed {MAX_TOTAL_UPLOAD_SIZE_MB} MB per analysis run."
            run_guard_message.set(msg)
            nav_lock_status.set(msg)
            _clear_input_busy_state(msg)
            log_run_guard_event(
                session_id=_safe_session_id(),
                reason="combined upload size exceeded",
                total_size_bytes=total_upload_size,
                run_attempt_count=len(updated_attempts),
            )
            return
        last_accepted_run_monotonic.set(now_monotonic)

        pipeline_running.set(True)
        pipeline_ready.set(False)
        _send_custom_message_safe(
            "set_input_busy_state",
            {
                "active": True,
                "title": "Running pipeline",
                "text": "Annotating datasets, rebuilding indexes, and scoring pathways.",
            },
        )
        try:
            _species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
            species_choice = str(species_info.get("label") or _species_key or "").strip().lower()
            species_code = str(species_info.get("code") or "").strip().lower()

            processed_protein = raw_protein
            try:
                kegg_map = _get_kegg_map(species_code)
                processed_protein = annotate_protein_with_kegg(processed_protein, species_code, kegg_map)
                headers = processed_protein.get("headers") or []
                kegg_idx = headers.index("KEGG_Gene_ID") if "KEGG_Gene_ID" in headers else -1
                kegg_matches = 0
                if kegg_idx >= 0:
                    for row in processed_protein.get("rows") or []:
                        if len(row) > kegg_idx and str(row[kegg_idx] or "").strip():
                            kegg_matches += 1
                if kegg_matches <= 5:
                    protein_kegg_warning.set("Few KEGG matches detected. Check species selection or mapping file.")
                else:
                    protein_kegg_warning.set("")
            except Exception as exc:
                print(f"Warning: KEGG annotation failed: {exc}")
                protein_kegg_warning.set("")
            try:
                processed_protein = _annotate_dataset_with_human_orthologs(processed_protein, species_code)
            except Exception as exc:
                print(f"Warning: human ortholog annotation failed: {exc}")

            processed_ptm: Optional[Dict[str, Any]] = None
            if raw_ptm and raw_ptm.get("headers"):
                processed_ptm = raw_ptm
                try:
                    processed_ptm = annotate_ptm_dataset(processed_ptm, species_choice, _get_psp_map())
                    processed_ptm = annotate_ptm_dataset_with_kinases(processed_ptm, species_choice, _get_psp_ks_map())
                except Exception as exc:
                    print(f"Warning: PSP annotation failed: {exc}")
                try:
                    processed_ptm = _annotate_dataset_with_human_orthologs(processed_ptm, species_code)
                except Exception as exc:
                    print(f"Warning: PTM human ortholog annotation failed: {exc}")

            processed_metabolite: Optional[Dict[str, Any]] = raw_metabolite if raw_metabolite and raw_metabolite.get("headers") else None

            protein_dataset.set(processed_protein)
            ptm_dataset.set(processed_ptm)
            metabolite_dataset.set(processed_metabolite)
            _write_debug_dump("user_protein_dataset_debug.txt", processed_protein)
            _write_debug_dump("user_ptm_dataset_debug.txt", processed_ptm or {"info": "No PTM dataset loaded"})
            await _refresh_global_catalog_from_current_async()
            _update_ks_index()
            _refresh_pathway_scores()
            pipeline_ready.set(True)
        except Exception as exc:
            print(f"Warning: input processing run failed: {exc}")
            _invalidate_user_run("Run failed during dataset processing.")
        finally:
            pipeline_running.set(False)
            _send_custom_message_safe(
                "set_input_busy_state",
                {
                    "active": False,
                    "title": "Processing input",
                    "text": "Validating files or preparing annotations, pathway scoring, and lookup indexes. This may take up to one minute.",
                },
            )

    @reactive.Effect
    @reactive.event(input.input_ptm_upload)
    def _process_ptm_upload():
        if not protein_validation.get().get("valid"):
            ptm_validation.set({"status": "Upload a valid protein file first.", "errors": ["Protein file not validated yet."], "valid": False})
            validated_ptm_dataset.set(None)
            ptm_upload_size_bytes.set(0)
            _invalidate_user_run("Upload and validate a Protein dataset before running.")
            return
        upload = input.input_ptm_upload()
        if not upload:
            ptm_validation.set({"status": "No PTM file uploaded (optional).", "errors": [], "valid": False})
            validated_ptm_dataset.set(None)
            ptm_preview_dataset.set(None)
            ptm_dataset_path.set(None)
            ptm_upload_size_bytes.set(0)
            _write_debug_dump("user_ptm_dataset_debug.txt", {"info": "No dataset loaded"})
            _invalidate_user_run("Validation updated. Press Run to process datasets.")
            return
        file_info = upload[0]
        upload_temp_path = _extract_upload_datapath(file_info)
        try:
            upload_validation = validate_uploaded_file(file_info, expected_role="ptm", session_id=_safe_session_id())
        except Exception:
            _cleanup_upload_temp_file(upload_temp_path, "PTM upload validation exception")
            msg = "Invalid analysis dataset upload: file could not be validated."
            ptm_validation.set({"status": msg, "errors": [msg], "valid": False})
            validated_ptm_dataset.set(None)
            ptm_preview_dataset.set(None)
            ptm_dataset_path.set(None)
            ptm_upload_size_bytes.set(0)
            _write_debug_dump("user_ptm_dataset_debug.txt", {"info": "No dataset loaded"})
            _invalidate_user_run("Fix the PTM dataset issue, then press Run.")
            return
        if not upload_validation.valid:
            _cleanup_upload_temp_file(upload_validation.datapath, "invalid PTM upload validation")
            ptm_validation.set({"status": upload_validation.user_message, "errors": [upload_validation.user_message], "valid": False})
            validated_ptm_dataset.set(None)
            ptm_preview_dataset.set(None)
            ptm_dataset_path.set(None)
            ptm_upload_size_bytes.set(0)
            _write_debug_dump("user_ptm_dataset_debug.txt", {"info": "No dataset loaded"})
            _invalidate_user_run("Fix the PTM dataset issue, then press Run.")
            return
        datapath = upload_validation.datapath
        upload_size_bytes = _upload_file_size_bytes(datapath)
        try:
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
                ptm_preview_dataset.set(dataset_payload)
                _apply_outline_width_defaults(ptm_headers=dataset_payload.get("headers") or [])
                validated_ptm_dataset.set(dataset_payload)
                ptm_dataset_path.set(upload_validation.sanitized_filename)
                ptm_upload_size_bytes.set(upload_size_bytes)
                _invalidate_user_run("PTM validated. Press Run to process datasets.")
            else:
                ptm_preview_dataset.set(_load_dataset(datapath))
                ptm_validation.set({"status": "PTM file failed validation. See errors below.", "errors": result.errors, "valid": False})
                validated_ptm_dataset.set(None)
                ptm_dataset_path.set(None)
                ptm_upload_size_bytes.set(0)
                _write_debug_dump("user_ptm_dataset_debug.txt", {"errors": result.errors})
                _invalidate_user_run("Fix the PTM dataset errors before running.")
        finally:
            _cleanup_upload_temp_file(datapath, "PTM upload parsed")

    @reactive.Effect
    def _enforce_nav_lock():
        mode = _current_mode()
        selected = _get_input_value(input, "bookmark_selector") or "home"
        protein_ok = bool(protein_validation.get().get("valid"))
        protein_blocked = _validation_has_explicit_error(protein_validation.get())
        ptm_blocked = _validation_has_explicit_error(ptm_validation.get())
        metabolite_blocked = _validation_has_explicit_error(metabolite_validation.get())
        ready = bool(pipeline_ready.get())
        if mode == "demo":
            status_msg = "Demo mode: sample datasets loaded. Navigation unlocked." if protein_ok else "Demo mode: loading sample datasets..."
            nav_lock_status.set(status_msg)
            locked = False
            _send_custom_message_safe("toggle_nav_lock", {"locked": locked})
            return
        locked = not ready
        if locked and selected not in {"input", "home"}:
            session.send_input_message("bookmark_selector", {"value": "home"})
        guard_text = str(run_guard_message.get() or "").strip()
        if guard_text:
            nav_lock_status.set(guard_text)
        elif ready:
            nav_lock_status.set("User mode: run complete. Navigation unlocked.")
        elif pipeline_running.get():
            nav_lock_status.set("User mode: running dataset processing. Tabs will unlock when processing finishes.")
        elif protein_blocked or ptm_blocked or metabolite_blocked:
            nav_lock_status.set("User mode: fix input errors, then press Run.")
        elif protein_ok:
            nav_lock_status.set("User mode: validation complete. Press Run to process datasets and unlock other tabs.")
        else:
            nav_lock_status.set("User mode: upload and validate a Protein dataset, then press Run. PTM and Metabolite are optional.")
        _send_custom_message_safe("toggle_nav_lock", {"locked": locked})

    @reactive.Effect
    @reactive.event(input.input_mode)
    async def _sync_mode_status():
        mode = str((_get_input_value(input, "input_mode") or "user")).lower()
        if mode_sync_in_progress.get():
            return
        previous_mode = str(last_applied_input_mode.get() or "").lower()
        if mode == previous_mode:
            return
        mode_sync_in_progress.set(True)
        if mode == "demo":
            try:
                try:
                    session.send_input_message("settings_use_black_protein_outlines", {"value": True})
                except Exception:
                    pass
                await _load_demo_datasets()
            except Exception as exc:
                print(f"Warning: demo mode failed to load sample datasets: {exc}")
                protein_dataset.set(None)
                ptm_dataset.set(None)
                protein_preview_dataset.set(None)
                ptm_preview_dataset.set(None)
                _reset_global_catalog_from_default()
                protein_validation.set({"status": "Demo mode failed to load sample protein.", "errors": [str(exc)], "valid": False, "comparisons": []})
                ptm_validation.set({"status": "Demo mode failed to load sample PTM.", "errors": [str(exc)], "valid": False})
                validated_protein_dataset.set(None)
                validated_ptm_dataset.set(None)
                validated_metabolite_dataset.set(None)
                run_guard_message.set("")
                last_accepted_run_monotonic.set(None)
                run_attempt_times_monotonic.set([])
                pipeline_running.set(False)
                pipeline_ready.set(False)
                nav_lock_status.set("Demo mode: sample datasets failed to load.")
                _clear_pathway_scores("Demo datasets failed to load; pathway scoring unavailable.")
            finally:
                last_applied_input_mode.set(mode)
                mode_sync_in_progress.set(False)
        else:
            try:
                demo_species_sync_target.set(None)
                try:
                    session.send_input_message("settings_use_black_protein_outlines", {"value": False})
                except Exception:
                    pass
                protein_dataset.set(None)
                ptm_dataset.set(None)
                metabolite_dataset.set(None)
                validated_protein_dataset.set(None)
                validated_ptm_dataset.set(None)
                validated_metabolite_dataset.set(None)
                protein_preview_dataset.set(None)
                ptm_preview_dataset.set(None)
                metabolite_preview_dataset.set(None)
                _reset_global_catalog_from_default()
                protein_kegg_warning.set("")
                protein_validation.set({"status": "Upload a protein file to begin.", "errors": [], "valid": False, "comparisons": []})
                ptm_validation.set({"status": "PTM upload optional. Provide after protein if available.", "errors": [], "valid": False})
                metabolite_validation.set({"status": "Metabolite upload optional. Provide after protein if available.", "errors": [], "valid": False, "comparisons": []})
                run_guard_message.set("")
                last_accepted_run_monotonic.set(None)
                run_attempt_times_monotonic.set([])
                pipeline_running.set(False)
                pipeline_ready.set(False)
                nav_lock_status.set("User mode: upload and validate files, then press Run to unlock other tabs.")
                _update_ks_index(reset=True)
                _clear_pathway_scores("Pathway scoring waiting for validated inputs and Run.")
            finally:
                last_applied_input_mode.set(mode)
                mode_sync_in_progress.set(False)

    @reactive.Effect
    @reactive.event(input.input_species)
    def _refresh_scores_for_species():
        with reactive.isolate():
            protein_ok = bool(protein_validation.get().get("valid"))
            pending_demo_species = str(demo_species_sync_target.get() or "").strip().lower()
        if pending_demo_species:
            demo_species_sync_target.set(None)
            current_species = str(_get_input_value(input, "input_species") or "").strip().lower()
            if _current_mode() == "demo" and current_species == pending_demo_species:
                return
        if protein_ok:
            if _current_mode() == "demo":
                _refresh_pathway_scores()
            else:
                _invalidate_user_run("Species changed. Press Run to reprocess datasets for the selected organism.")

    def _gradient_stops_from_settings_inputs() -> List[Dict[str, Any]]:
        payload = _get_input_value(input, "settings_gradient_stops_json")
        parsed_stops: List[Dict[str, Any]] = []
        if isinstance(payload, dict):
            parsed_stops = _normalize_gradient_stop_entries(payload.get("rows"))
        elif isinstance(payload, str) and payload.strip():
            try:
                parsed_obj = json.loads(payload)
                if isinstance(parsed_obj, dict):
                    parsed_stops = _normalize_gradient_stop_entries(parsed_obj.get("rows"))
                elif isinstance(parsed_obj, list):
                    parsed_stops = _normalize_gradient_stop_entries(parsed_obj)
            except Exception:
                parsed_stops = []
        if len(parsed_stops) >= 2:
            return parsed_stops
        neg_hex = str(_get_input_value(input, "settings_negative_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["negative_color"]))
        pos_hex = str(_get_input_value(input, "settings_positive_color") or _rgb_tuple_to_hex(DEFAULT_SETTINGS["positive_color"]))
        neg_val = _to_float(_get_input_value(input, "settings_max_negative"), DEFAULT_SETTINGS["max_negative"])
        pos_val = _to_float(_get_input_value(input, "settings_max_positive"), DEFAULT_SETTINGS["max_positive"])
        return _normalize_gradient_stop_entries(
            [
                {"value": neg_val, "color": neg_hex},
                {"value": 0.0, "color": "#ffffff"},
                {"value": pos_val, "color": pos_hex},
            ]
        ) or _default_gradient_stops()

    def _gradient_editor_rows() -> List[Dict[str, Any]]:
        stops = _gradient_stops_from_settings_inputs()
        rows = sorted(stops, key=lambda item: float(item["value"]), reverse=True)
        next_idx = len(rows) + 1
        rows.append({"value": None, "color": "#ffffff", "_row_id": f"add_{next_idx}"})
        for idx, row in enumerate(rows, start=1):
            if "_row_id" not in row:
                row["_row_id"] = f"stop_{idx}"
        return rows

    def _gradient_preview_style_and_labels() -> Tuple[str, List[str]]:
        stops = _gradient_stops_from_settings_inputs()
        if len(stops) < 2:
            stops = _default_gradient_stops()
        min_val = float(stops[0]["value"])
        max_val = float(stops[-1]["value"])
        span = (max_val - min_val) if max_val != min_val else 1.0
        gradient_parts: List[str] = []
        for stop in stops:
            value = float(stop["value"])
            pct = max(0.0, min(100.0, ((value - min_val) / span) * 100.0))
            color_hex = _rgb_tuple_to_hex(tuple(int(c) for c in list(stop["color"])[:3]))
            gradient_parts.append(f"{color_hex} {pct:.3f}%")
        style = "background: linear-gradient(90deg, " + ", ".join(gradient_parts) + ");"
        labels = [f"{float(stop['value']):g}" for stop in sorted(stops, key=lambda item: float(item["value"]))]
        return style, labels

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
        preset_rows = [
            {"value": float(DEFAULT_SETTINGS.get("max_negative", -2)), "color": neg_hex},
            {"value": 0.0, "color": "#ffffff"},
            {"value": float(DEFAULT_SETTINGS.get("max_positive", 2)), "color": pos_hex},
        ]
        try:
            session.send_input_message("settings_negative_color", {"value": neg_hex})
            session.send_input_message("settings_positive_color", {"value": pos_hex})
            session.send_input_message(
                "settings_gradient_stops_json",
                {"value": json.dumps({"rows": preset_rows})},
            )
        except Exception:
            pass

    @output
    @render.ui
    def settings_gradient_table():
        rows = _gradient_editor_rows()
        body_rows: List[Any] = []
        for row in rows:
            row_id = str(row.get("_row_id") or "")
            raw_value = row.get("value")
            value_str = "" if raw_value in (None, "") else f"{float(raw_value):g}"
            raw_color = row.get("color")
            color_val = "#ffffff"
            if isinstance(raw_color, str) and re.match(r"^#[0-9a-fA-F]{6}$", raw_color):
                color_val = raw_color
            elif isinstance(raw_color, (list, tuple)) and len(raw_color) >= 3:
                try:
                    color_val = _rgb_tuple_to_hex(tuple(int(float(c)) for c in list(raw_color)[:3]))
                except Exception:
                    color_val = "#ffffff"
            value_input = ui.tags.input(
                {
                    "type": "number",
                    "step": "any",
                    "class": "mk-grad-value",
                    "value": value_str,
                    "placeholder": "Add value",
                    "onchange": "window.mkGradientTableUpdate && window.mkGradientTableUpdate();",
                    "data-row-id": row_id,
                }
            )
            color_input = ui.tags.input(
                {
                    "type": "color",
                    "class": "mk-grad-color",
                    "value": color_val,
                    "onchange": "window.mkGradientTableUpdate && window.mkGradientTableUpdate();",
                    "data-row-id": row_id,
                }
            )
            body_rows.append(ui.tags.tr({"class": "mk-gradient-row"}, ui.tags.td(value_input), ui.tags.td(color_input)))
        return ui.div(
            {"class": "gradient-table-wrap"},
            ui.tags.table(
                {"id": "settings_gradient_table_editor", "class": "gradient-table"},
                ui.tags.thead(ui.tags.tr(ui.tags.th("Value"), ui.tags.th("Color"))),
                ui.tags.tbody(*body_rows),
            ),
        )

    @output
    @render.ui
    def settings_gradient_preview():
        style, labels = _gradient_preview_style_and_labels()
        if not labels:
            labels = ["-2", "0", "2"]
        return ui.div(
            {"class": "gradient-preview", "style": style},
            *[ui.tags.span(label) for label in labels],
        )

    @output
    @render.ui
    def input_data_inputs_panel():
        mode = _current_mode()
        protein_ok = bool(protein_validation.get().get("valid"))
        disabled = mode == "demo"
        ptm_disabled = disabled or not protein_ok
        metabolite_disabled = disabled or not protein_ok
        species_wrap_style = "opacity:0.55; pointer-events:none;" if disabled else ""
        protein_wrap_style = "opacity:0.55; pointer-events:none;" if disabled else ""
        ptm_wrap_style = "opacity:0.55; pointer-events:none;" if ptm_disabled else ""
        metabolite_wrap_style = "opacity:0.55; pointer-events:none;" if metabolite_disabled else ""
        current_species = _get_input_value(input, "input_species") or DEFAULT_SPECIES
        if current_species not in SPECIES_OPTIONS:
            current_species = DEFAULT_SPECIES

        pv = protein_validation.get()
        ptv = ptm_validation.get()
        protein_valid = bool(pv.get("valid"))
        ptm_valid = bool(ptv.get("valid"))
        metabolite_valid = bool(metabolite_validation.get().get("valid"))
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
        if not error_flag:
            for entry in (metabolite_validation.get().get("errors") or []):
                if entry:
                    error_flag = True
                    break
        # Also treat statuses containing "failed" or "missing" as errors
        statuses_to_check = [str(pv.get("status", "")), str(ptv.get("status", ""))]
        if not error_flag and any("failed" in s.lower() or "missing" in s.lower() for s in statuses_to_check):
            error_flag = True

        data_card = ui.card(
            {"class": "data-input-card", "style": "max-width: 400px; min-width: 350px; width: 380px; height: 224px;"},
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
                    ui.div(
                        ui.input_select(
                            "input_species",
                            "Select Species",
                            choices=SPECIES_OPTIONS,
                            selected=current_species,
                        ),
                        ui.input_action_button(
                            "input_run_pipeline",
                            "Run",
                            class_="btn input-run-btn",
                        ),
                        ui.div(
                            {"class": "input-run-status"},
                            ui.output_text("input_nav_lock_status"),
                        ),
                        ui.output_ui("input_run_button_state"),
                    ),
                ),
            ),
        )

        add_card = ui.div(style="display:none;")

        ptm_card = ui.card(
            {"class": "data-input-card", "style": "max-width: 400px; min-width: 350px; width: 380px; height: 224px;"},
            ui.card_header(
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.div(
                        {"style": "display:flex; align-items:center; gap:8px;"},
                        ui.tags.span("PTM Dataset", style="font-weight:700; color:#0b1f33;"),
                        ui.tags.span(
                            "OK",
                            style="color:#16a34a; font-weight:700;" if ptm_valid else "color:transparent; font-weight:700;",
                        ),
                    ),
                )
            ),
            ui.card_body(
                ui.div(
                    {"style": ptm_wrap_style},
                    ui.div(
                        {"style": "display:flex; align-items:flex-end; gap:8px;"},
                        ui.input_file(
                            "input_ptm_upload",
                            ui.TagList("PTM Upload ", ui.tags.em("(.txt, .tsv, .csv)")),
                            multiple=False,
                            accept=UPLOAD_ACCEPT_TYPES,
                        ),
                        _ptm_shape_picker(),
                    ),
                ),
                ui.div(
                    {"class": "input-sample-download-row", "style": "margin-top:-18px;"},
                    ui.download_button(
                        "download_sample_ptm_dataset",
                        _icon_markup("download", "viewer-create-icon"),
                        class_="btn btn-outline-secondary btn-sm input-sample-download-btn",
                    ),
                    ui.tags.span({"class": "input-sample-download-label"}, "Sample PTM Dataset"),
                )
            ),
        )

        protein_card = ui.card(
            {"class": "data-input-card", "style": "max-width: 400px; min-width: 350px; width: 380px; height: 224px;"},
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
                            ui.TagList("Protein Dataset Upload ", ui.tags.em("(.txt, .tsv, .csv)")),
                            multiple=False,
                            accept=UPLOAD_ACCEPT_TYPES,
                        ),
                        ui.tags.span(
                            ui.output_text("input_protein_kegg_warning"),
                            style="color: #b00020; font-weight:600; font-size: 12px;",
                        ),
                    ),
                    ui.div(
                        "Demo mode uses bundled sample files. Switch to User mode to upload your own.",
                        style="color: #b00020; font-size: 12px;" if disabled else "display:none;",
                    ),
                ),
                ui.div(
                    {"class": "input-sample-download-row", "style": "margin-top:-18px;"},
                    ui.download_button(
                        "download_sample_protein_dataset",
                        _icon_markup("download", "viewer-create-icon"),
                        class_="btn btn-outline-secondary btn-sm input-sample-download-btn",
                    ),
                    ui.tags.span({"class": "input-sample-download-label"}, "Sample Protein Dataset"),
                ),
            ),
        )

        metabolite_card = ui.card(
            {"class": "data-input-card", "style": "max-width: 400px; min-width: 350px; width: 380px; height: 224px;"},
            ui.card_header(
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.tags.span("Metabolite Dataset", style="font-weight:700; color:#0b1f33;"),
                    ui.tags.span(
                        "OK",
                        style="color:#16a34a; font-weight:700;" if metabolite_valid else "color:transparent; font-weight:700;",
                    ),
                )
            ),
            ui.card_body(
                ui.div(
                    {"style": metabolite_wrap_style},
                    ui.input_file(
                        "input_metabolite_upload",
                        ui.TagList("Metabolite Upload ", ui.tags.em("(.txt, .tsv, .csv)")),
                        multiple=False,
                        accept=UPLOAD_ACCEPT_TYPES,
                    ),
                    ui.div(
                        "Demo mode uses bundled sample files. Switch to User mode to upload your own.",
                        style="margin-top:6px; color:#b00020; font-size:12px;" if disabled else "display:none;",
                    ),
                ),
                ui.div(
                    {"class": "input-sample-download-row", "style": "margin-top:-18px;"},
                    ui.download_button(
                        "download_sample_metabolite_dataset",
                        _icon_markup("download", "viewer-create-icon"),
                        class_="btn btn-outline-secondary btn-sm input-sample-download-btn",
                    ),
                    ui.tags.span({"class": "input-sample-download-label"}, "Sample Metabolite Dataset"),
                ),
            ),
        )

        error_block = ui.card(
            {"class": "data-input-card", "style": "max-width: 420px; min-width: 350px; width: 380px; border-color:#fecaca;"},
            ui.card_header(
                ui.div(
                    {"style": "display:flex; align-items:center; gap:8px;"},
                    ui.tags.span("Input Errors", style="font-weight:700; color:#7f1d1d;"),
                )
            ),
            ui.card_body(
                ui.div(
                    ui.output_ui("input_protein_errors"),
                    ui.output_ui("input_ptm_errors"),
                    ui.output_ui("input_metabolite_errors"),
                    style="color:#b00020; font-size:13px;",
                )
            ),
        )

        return ui.TagList(
            ui.div(
                {"class": "input-page-stack"},
                ui.div(
                    {"id": "mk-input-busy-overlay", "class": "mk-input-busy-overlay", "aria-live": "polite", "aria-hidden": "true"},
                    ui.div(
                        {"class": "mk-input-busy-card", "role": "status"},
                        ui.div({"class": "mk-input-busy-spinner"}),
                        ui.div(
                            {"class": "mk-input-busy-copy"},
                            ui.div({"class": "mk-input-busy-title"}, "Processing input"),
                            ui.div(
                                {"class": "mk-input-busy-text"},
                                "Validating files or preparing annotations, pathway scoring, and lookup indexes. This may take up to one minute.",
                            ),
                        ),
                    ),
                ),
                ui.div(
                    {
                        "class": "input-page-row",
                        "style": (
                            "border-left:4px solid #f59e0b; background:#fffbeb; color:#78350f; "
                            "padding:5px 8px; font-size:11px; font-weight:600; line-height:1.25; "
                            "max-width:760px; box-sizing:border-box;"
                        ),
                    },
                    DATA_USE_WARNING_TEXT,
                ),
                ui.div(
                    {"class": "input-page-row", "style": "display:flex; flex-wrap:wrap; gap:16px; align-items:flex-start;"},
                    data_card,
                    protein_card,
                    ptm_card,
                    metabolite_card,
                    add_card,
                ),
                ui.div({"class": "input-page-row"}, ui.output_ui("input_protein_preview_guide")),
                ui.div({"class": "input-page-row input-preview-section"}, ui.output_ui("input_protein_preview")),
                ui.div(
                    {"class": "input-page-row"},
                    error_block,
                    style="display:block;" if error_flag else "display:none;",
                ),
                ui.div({"class": "input-page-row input-preview-section"}, ui.output_ui("input_ptm_preview")),
                ui.div({"class": "input-page-row input-preview-section"}, ui.output_ui("input_metabolite_preview")),
            )
        )

    def _color_override_from_settings(settings_override: Dict[str, Any]) -> Dict[str, Any]:
        gradient_stops = _normalize_gradient_stops_with_fallback(settings_override.get("gradient_stops"))
        min_stop = gradient_stops[0]
        max_stop = gradient_stops[-1]
        return {
            "negative_color": _rgb_tuple_to_list(settings_override.get("negative_color", tuple(min_stop["color"]))),
            "positive_color": _rgb_tuple_to_list(settings_override.get("positive_color", tuple(max_stop["color"]))),
            "max_negative": settings_override.get("max_negative", float(min_stop["value"])),
            "max_positive": settings_override.get("max_positive", float(max_stop["value"])),
            "gradient_stops": gradient_stops,
        }

    def _apply_color_metadata(payload: Dict[str, Any], color_override: Dict[str, Any]) -> None:
        general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
        general_settings["negative_color"] = color_override["negative_color"]
        general_settings["positive_color"] = color_override["positive_color"]
        general_settings["max_negative"] = color_override["max_negative"]
        general_settings["max_positive"] = color_override["max_positive"]
        general_settings["gradient_stops"] = color_override.get("gradient_stops", _default_gradient_stops())
        payload["_color_preview_override"] = color_override
        _apply_color_overrides(payload, color_override)

    def _register_bookmark(cfg: Dict[str, Any]) -> None:
        prefix = cfg["key"]
        state = bookmark_state[prefix]
        web_sort_col = reactive.Value("prot_fisher_p") if prefix == "web" else None
        web_sort_dir = reactive.Value("asc") if prefix == "web" else None
        web_selected_label = reactive.Value("")
        web_selected_id = reactive.Value("") if prefix == "web" else None
        web_selected_source = reactive.Value("wikipathways") if prefix == "web" else None
        web_filter_options = reactive.Value(list(WEB_PATHWAY_SOURCES)) if prefix == "web" else None
        web_filter_selected = reactive.Value(list(WEB_PATHWAY_SOURCES)) if prefix == "web" else None
        web_filter_tick = reactive.Value(0) if prefix == "web" else None
        web_filter_refresh_evt = reactive.Value(0) if prefix == "web" else None
        web_page = reactive.Value(1) if prefix == "web" else None

        def _set_load_feedback(message: str = "") -> None:
            state["load_feedback"].set(str(message or "").strip())

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
        ks_last_choice = reactive.Value("")
        ks_last_choice_cache: Dict[str, str] = {"value": ""}

        def _remember_ks_choice(value: Any) -> None:
            choice_val = str(value or "")
            ks_last_choice_cache["value"] = choice_val
            try:
                ks_last_choice.set(choice_val)
            except Exception:
                pass

        def _last_ks_choice() -> str:
            try:
                return str(ks_last_choice.get() or "")
            except Exception:
                return str(ks_last_choice_cache.get("value") or "")

        def _ks_selected_types() -> set:
            raw = _get_input_value(input, _prefixed_id(prefix, "ks_evidence_types"))
            if raw in (None, "", False):
                return {"in_vivo"}
            values = raw if isinstance(raw, (list, tuple, set)) else [raw]
            selected = {str(value).strip().lower() for value in values if value not in (None, "", False)}
            return selected.intersection({"in_vivo", "in_vitro"})

        def _ks_mode_value() -> str:
            selected = _ks_selected_types()
            if selected == {"in_vivo", "in_vitro"}:
                return "both"
            if selected == {"in_vitro"}:
                return "in_vitro"
            if selected == {"in_vivo"}:
                return "in_vivo"
            return "none"

        def _ks_allowed_types(mode: str) -> set:
            if mode == "both":
                return {"in_vivo", "in_vitro"}
            if mode == "in_vivo":
                return {"in_vivo"}
            if mode == "in_vitro":
                return {"in_vitro"}
            return set()

        def _ks_evidence_label(types: set) -> str:
            clean_types = set(types or set())
            if {"in_vivo", "in_vitro"}.issubset(clean_types):
                return "both"
            if "in_vivo" in clean_types:
                return "in vivo"
            if "in_vitro" in clean_types:
                return "in vitro"
            return ""

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
            catalog_info = _current_global_catalog_info()
            payload = _build_blank_canvas(catalog_info)
            _apply_color_metadata(payload, color_override)
            general_settings = payload.setdefault("general_data", {}).setdefault("settings", {})
            general_settings["show_background_image"] = False
            general_settings["show_groups"] = bool(settings_override.get("show_groups", False))
            general_settings["show_multi_protein_indicator"] = bool(settings_override.get("show_multi_protein_indicator", False))
            general_settings["show_arrows"] = True
            general_settings["show_text_boxes"] = True
            general_settings["debug_mode"] = bool(settings_override.get("debug_mode", False))
            general_settings["temporal_mode"] = bool(settings_override.get("temporal_mode", DEFAULT_SETTINGS.get("temporal_mode", False)))
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
            general_settings["use_black_protein_outlines"] = bool(settings_override.get("use_black_protein_outlines", DEFAULT_SETTINGS.get("use_black_protein_outlines", False)))
            general_settings["negative_color"] = color_override["negative_color"]
            general_settings["positive_color"] = color_override["positive_color"]
            general_settings["max_negative"] = color_override["max_negative"]
            general_settings["max_positive"] = color_override["max_positive"]
            general_settings["main_columns"] = settings_override["main_columns"]
            general_settings["protein_tooltip_columns"] = settings_override.get("protein_tooltip_columns", [])
            general_settings["species"] = settings_override.get("species")
            general_settings["species_code"] = settings_override.get("species_code")
            general_settings["bg_scale"] = settings_override.get("bg_scale", DEFAULT_WIKIPATHWAYS_BG_SCALE)
            general_settings["bg_offset_x"] = settings_override.get("bg_offset_x", DEFAULT_WIKIPATHWAYS_BG_OFFSET_X)
            general_settings["bg_offset_y"] = settings_override.get("bg_offset_y", DEFAULT_WIKIPATHWAYS_BG_OFFSET_Y)
            payload["_kegg_preview"] = _background_preview_from_settings(settings_override, show=False)
            payload["_box_preview"] = {"y_stretch": DEFAULT_BOX_Y_STRETCH}
            payload["_active_mode"] = settings_override.get("mode", cfg["mode"])
            payload["_full_width_canvas"] = True
            payload["_custom_layout_applied"] = False
            payload["_bookmark_key"] = cfg.get("key", prefix)
            payload["_persist_token"] = time.time()
            payload["_global_protein_catalog"] = dict(catalog_info)
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
            gradient_stops = _normalize_gradient_stops_with_fallback(settings_override.get("gradient_stops"))
            outline_cols = ctx.get("protein_outline_columns", [])
            use_black_protein_outlines = bool(settings_override.get("use_black_protein_outlines", DEFAULT_SETTINGS.get("use_black_protein_outlines", False)))
            for idx, col in enumerate(ctx.get("protein_main_columns", []), 1):
                raw = protein_row.get(col, "")
                try:
                    fc_val = float(raw)
                except (TypeError, ValueError):
                    fc_val = None
                entry[f"fold_change_{idx}"] = fc_val
                entry[f"fc_color_{idx}"] = _gradient_color_from_fold(
                    fc_val, neg, pos, max_neg, max_pos, gradient_stops=gradient_stops
                )
                outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                outline_raw = protein_row.get(outline_col, "") if outline_col else ""
                try:
                    outline_val = float(outline_raw)
                except (TypeError, ValueError):
                    outline_val = None
                entry[f"outline_fold_change_{idx}"] = outline_val
                entry[f"outline_color_{idx}"] = (
                    [0, 0, 0]
                    if use_black_protein_outlines
                    else (
                        _gradient_color_from_fold(
                            outline_val,
                            neg,
                            pos,
                            max_neg,
                            max_pos,
                            gradient_stops=gradient_stops,
                        )
                        if outline_val is not None
                        else [0, 0, 0]
                    )
                )
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
            gradient_stops = _normalize_gradient_stops_with_fallback(settings_override.get("gradient_stops"))
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
                payload[f"fc_color_{idx}"] = _gradient_color_from_fold(
                    fc_val, neg, pos, max_neg, max_pos, gradient_stops=gradient_stops
                )
                outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                outline_raw = row_map.get(outline_col, "") if outline_col else ""
                try:
                    outline_val = float(outline_raw)
                except (TypeError, ValueError):
                    outline_val = None
                payload[f"outline_fold_change_{idx}"] = outline_val
                payload[f"outline_color_{idx}"] = (
                    _gradient_color_from_fold(
                        outline_val,
                        neg,
                        pos,
                        max_neg,
                        max_pos,
                        gradient_stops=gradient_stops,
                    )
                    if outline_val is not None
                    else [0, 0, 0]
                )
            return payload

        def _hydrate_custom_layout_from_current_dataset(payload: Dict[str, Any], ctx: Dict[str, Any], settings_override: Dict[str, Any]) -> Dict[str, int]:
            if not isinstance(payload, dict):
                return {"boxes_with_ids": 0, "matched_boxes": 0, "matched_proteins": 0}
            protboxes = payload.get("protbox_data") or []
            if not isinstance(protboxes, list):
                return {"boxes_with_ids": 0, "matched_boxes": 0, "matched_proteins": 0}
            protein_rows_raw = ctx.get("prot_rows_by_uniprot", {}) or {}
            ptm_rows_raw = ctx.get("ptms_by_uniprot", {}) or {}
            protein_rows_by_norm: Dict[str, Tuple[str, Dict[str, Any]]] = {}
            ptm_rows_by_norm: Dict[str, List[Dict[str, Any]]] = {}
            for uid, row_map in protein_rows_raw.items():
                norm_uid = normalize_uniprot(uid)
                if norm_uid and norm_uid not in protein_rows_by_norm:
                    protein_rows_by_norm[norm_uid] = (str(uid), row_map)
            for uid, rows in ptm_rows_raw.items():
                norm_uid = normalize_uniprot(uid)
                if norm_uid and norm_uid not in ptm_rows_by_norm:
                    ptm_rows_by_norm[norm_uid] = list(rows or [])

            protein_data: Dict[str, Any] = {}
            boxes_with_ids = 0
            matched_boxes = 0
            matched_proteins = set()
            ptm_limit = int(settings_override.get("ptm_max_display", DEFAULT_SETTINGS["ptm_max_display"]) or DEFAULT_SETTINGS["ptm_max_display"])
            if ptm_limit <= 0:
                ptm_limit = len(PTM_POSITION_PRIORITY)
            ptm_limit = min(ptm_limit, len(PTM_POSITION_PRIORITY))

            for pb in protboxes:
                if not isinstance(pb, dict):
                    continue
                candidates = _custom_layout_protbox_uniprots(pb)
                if candidates:
                    boxes_with_ids += 1
                matched_ids: List[str] = []
                for candidate in candidates:
                    match = protein_rows_by_norm.get(candidate)
                    if not match:
                        continue
                    dataset_uid, protein_row = match
                    if dataset_uid not in matched_ids:
                        matched_ids.append(dataset_uid)
                    if dataset_uid in protein_data:
                        continue
                    ptm_entries: Dict[str, Any] = {}
                    ptm_rows = ptm_rows_by_norm.get(candidate, [])
                    box_x = float(_coerce_float(pb.get("x")) or 0.0)
                    box_y = float(_coerce_float(pb.get("y")) or 0.0)
                    box_width = float(_coerce_float(pb.get("width")) or 46.0)
                    box_height = float(_coerce_float(pb.get("height")) or 17.0)
                    for ptm_idx, row_map in enumerate(ptm_rows[:ptm_limit]):
                        site_label = _clean_site_label(
                            str(row_map.get(ctx.get("ptm_headers", [None, None])[1] if ctx.get("ptm_headers") else "", "") or "")
                        )
                        ptm_key = f"{dataset_uid}_{site_label or (ptm_idx + 1)}"
                        ptm_entries[ptm_key] = _make_ptm_entry(
                            row_map,
                            ptm_key,
                            PTM_POSITION_PRIORITY[ptm_idx],
                            box_x,
                            box_y,
                            box_width,
                            box_height,
                            ctx,
                            settings_override,
                            display_label=site_label or None,
                        )
                    protein_data[dataset_uid] = _make_protein_entry(dataset_uid, protein_row, ctx, settings_override, ptm_entries)
                    matched_proteins.add(dataset_uid)
                if matched_ids:
                    matched_boxes += 1
                    pb["proteins"] = matched_ids
                    pb["selected_uniprot"] = matched_ids[0]
                    primary_entry = protein_data.get(matched_ids[0]) or {}
                    pb["backup_label"] = primary_entry.get("label") or pb.get("backup_label") or matched_ids[0]
                    pb["tooltip"] = primary_entry.get("tooltip") or pb.get("tooltip") or ""
                    pb["tooltip_html"] = primary_entry.get("tooltip_html") or pb.get("tooltip_html") or ""
                elif candidates:
                    pb["proteins"] = candidates
                    pb["selected_uniprot"] = candidates[0]
            payload["protein_data"] = protein_data
            return {
                "boxes_with_ids": boxes_with_ids,
                "matched_boxes": matched_boxes,
                "matched_proteins": len(matched_proteins),
            }

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

        def _ks_dataset_kinases_only() -> bool:
            return bool(_get_input_value(input, _prefixed_id(prefix, "ks_dataset_kinases_only")))

        def _ks_protein_in_dataset(uniprot_id: str, ctx: Dict[str, Any]) -> bool:
            prot_rows = ctx.get("prot_rows_by_uniprot", {}) or {}
            return str(uniprot_id or "") in prot_rows

        def _ks_kinase_in_dataset(kinase_id: str, ctx: Dict[str, Any]) -> bool:
            return _ks_protein_in_dataset(kinase_id, ctx)

        def _ks_layout_config() -> Dict[str, Any]:
            mode = str(_get_input_value(input, _prefixed_id(prefix, "ks_layout_mode")) or "two_column").strip().lower()
            gap_x = _to_float(_get_input_value(input, _prefixed_id(prefix, "ks_layout_gap_x")), KS_LAYOUT_GAP_X_DEFAULT)
            spacing_y = _to_float(_get_input_value(input, _prefixed_id(prefix, "ks_layout_spacing_y")), KS_VERTICAL_SPACING)
            radius_mode = str(_get_input_value(input, _prefixed_id(prefix, "ks_conc_radius_mode")) or "auto").strip().lower()
            fixed_radius = _to_float(_get_input_value(input, _prefixed_id(prefix, "ks_conc_radius_fixed")), KS_CONCENTRIC_RADIUS_DEFAULT)
            auto_space = _to_float(_get_input_value(input, _prefixed_id(prefix, "ks_conc_space")), KS_CONCENTRIC_SPACE_DEFAULT)
            arrow_stop = _to_float(_get_input_value(input, _prefixed_id(prefix, "ks_conc_arrow_stop")), KS_CONCENTRIC_ARROW_STOP_DEFAULT)
            two_col_focus_side = str(_get_input_value(input, _prefixed_id(prefix, "ks_two_col_focus_side")) or "left").strip().lower()
            two_col_columns = int(
                max(1, _to_float(_get_input_value(input, _prefixed_id(prefix, "ks_two_col_columns")), 1) or 1)
            )
            return {
                "mode": mode if mode in {"two_column", "concentric"} else "two_column",
                "gap_x": max(80.0, float(gap_x if gap_x is not None else KS_LAYOUT_GAP_X_DEFAULT)),
                "spacing_y": max(10.0, float(spacing_y if spacing_y is not None else KS_VERTICAL_SPACING)),
                "center_x": KS_LAYOUT_CENTER_X_DEFAULT,
                "center_y": KS_LAYOUT_CENTER_Y_DEFAULT,
                "two_col_focus_side": two_col_focus_side if two_col_focus_side in {"left", "right", "top", "bottom"} else "left",
                "two_col_columns": min(6, max(1, two_col_columns)),
                "conc_add_arrows": bool(_get_input_value(input, _prefixed_id(prefix, "ks_conc_arrows"))),
                "conc_radius_mode": radius_mode if radius_mode in {"auto", "fixed"} else "auto",
                "conc_radius_fixed": max(20.0, float(fixed_radius if fixed_radius is not None else KS_CONCENTRIC_RADIUS_DEFAULT)),
                "conc_space": max(10.0, float(auto_space if auto_space is not None else KS_CONCENTRIC_SPACE_DEFAULT)),
                "conc_arrow_stop": max(0.0, float(arrow_stop if arrow_stop is not None else KS_CONCENTRIC_ARROW_STOP_DEFAULT)),
            }

        def _ks_two_col_arrow_sides(
            focus_side: str, direction: str = "single_to_multi"
        ) -> Tuple[str, str]:
            side = str(focus_side or "left").strip().lower()
            single_to_multi = {
                "left": ("East", "West"),
                "right": ("West", "East"),
                "top": ("South", "North"),
                "bottom": ("North", "South"),
            }
            source_side, target_side = single_to_multi.get(side, ("East", "West"))
            if str(direction or "").strip().lower() == "multi_to_single":
                return target_side, source_side
            return source_side, target_side

        def _ks_stable_persist_token(
            selection_value: str,
            mode: str,
            ctx: Dict[str, Any],
            settings_override: Dict[str, Any],
        ) -> str:
            token_payload = {
                "selection": str(selection_value or ""),
                "mode": str(mode or ""),
                "layout": _ks_layout_config(),
                "dataset_kinases_only": _ks_dataset_kinases_only(),
                "ptm_max_display": int(settings_override.get("ptm_max_display", DEFAULT_SETTINGS["ptm_max_display"]) or DEFAULT_SETTINGS["ptm_max_display"]),
                "protein_main_columns": list(ctx.get("protein_main_columns") or []),
                "ptm_main_columns": list(ctx.get("ptm_main_columns") or []),
                "protein_outline_columns": list(ctx.get("protein_outline_columns") or []),
                "ptm_outline_columns": list(ctx.get("ptm_outline_columns") or []),
                "use_black_protein_outlines": bool(
                    settings_override.get(
                        "use_black_protein_outlines",
                        DEFAULT_SETTINGS.get("use_black_protein_outlines", False),
                    )
                ),
            }
            digest = hashlib.sha1(
                json.dumps(token_payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
            ).hexdigest()
            return f"ks:{digest}"

        def _ks_concentric_radius(node_count: int, layout_cfg: Dict[str, Any]) -> float:
            if node_count <= 0:
                return KS_CONCENTRIC_RADIUS_DEFAULT
            if layout_cfg.get("conc_radius_mode") == "fixed":
                return float(layout_cfg.get("conc_radius_fixed", KS_CONCENTRIC_RADIUS_DEFAULT))
            space = float(layout_cfg.get("conc_space", KS_CONCENTRIC_SPACE_DEFAULT))
            return max(60.0, (node_count * space) / (2 * math.pi) * 2.0)

        def _ks_two_column_positions(
            node_ids: List[str],
            layout_cfg: Dict[str, Any],
            center_x: float,
            center_y: float,
            half_w: float,
            half_h: float,
            protbox_width: float,
        ) -> Tuple[float, float, Dict[str, Tuple[float, float]], set]:
            spacing_y = float(layout_cfg.get("spacing_y", KS_VERTICAL_SPACING))
            gap_x = float(layout_cfg.get("gap_x", KS_LAYOUT_GAP_X_DEFAULT))
            focus_side = str(layout_cfg.get("two_col_focus_side", "left"))
            columns = int(layout_cfg.get("two_col_columns", 1) or 1)
            columns = max(1, columns)
            # Keep additional columns much tighter so multi-column KS layouts remain readable.
            if columns > 1:
                col_spacing = max(protbox_width + 10.0, min(120.0, gap_x * 0.32))
            else:
                col_spacing = max(protbox_width + 16.0, min(220.0, gap_x * 0.75))

            if focus_side == "right":
                rows = max(1, math.ceil(len(node_ids) / columns)) if node_ids else 1
                total_span = spacing_y * (rows - 1)
                start_y = center_y - (total_span / 2.0) - half_h
                focus_x = center_x + (gap_x * 0.5) - half_w
                near_x = center_x - (gap_x * 0.5) - half_w
                col_x = lambda col_idx: near_x - (col_idx * col_spacing)
                row_y = lambda row_idx: start_y + row_idx * spacing_y
            else:
                if focus_side == "top":
                    depth_rows = columns
                    cols_across = max(1, math.ceil(len(node_ids) / depth_rows)) if node_ids else 1
                    across_spacing = max(protbox_width + 8.0, min(160.0, spacing_y))
                    depth_spacing = max(protbox_width + 10.0, min(120.0, gap_x * 0.32))
                    focus_y = center_y - (gap_x * 0.5) - half_h
                    near_y = center_y + (gap_x * 0.5) - half_h
                    start_x = center_x - (across_spacing * (cols_across - 1) * 0.5) - half_w
                    col_x = lambda col_idx: start_x + col_idx * across_spacing
                    row_y = lambda row_idx: near_y + (row_idx * depth_spacing)
                    focus_x = center_x - half_w
                elif focus_side == "bottom":
                    depth_rows = columns
                    cols_across = max(1, math.ceil(len(node_ids) / depth_rows)) if node_ids else 1
                    across_spacing = max(protbox_width + 8.0, min(160.0, spacing_y))
                    depth_spacing = max(protbox_width + 10.0, min(120.0, gap_x * 0.32))
                    focus_y = center_y + (gap_x * 0.5) - half_h
                    near_y = center_y - (gap_x * 0.5) - half_h
                    start_x = center_x - (across_spacing * (cols_across - 1) * 0.5) - half_w
                    col_x = lambda col_idx: start_x + col_idx * across_spacing
                    row_y = lambda row_idx: near_y - (row_idx * depth_spacing)
                    focus_x = center_x - half_w
                else:
                    rows = max(1, math.ceil(len(node_ids) / columns)) if node_ids else 1
                    total_span = spacing_y * (rows - 1)
                    start_y = center_y - (total_span / 2.0) - half_h
                    focus_x = center_x - (gap_x * 0.5) - half_w
                    near_x = center_x + (gap_x * 0.5) - half_w
                    col_x = lambda col_idx: near_x + (col_idx * col_spacing)
                    row_y = lambda row_idx: start_y + row_idx * spacing_y

            positions: Dict[str, Tuple[float, float]] = {}
            arrow_targets: set = set()
            for idx, node_id in enumerate(node_ids):
                if focus_side in {"top", "bottom"}:
                    depth_rows = columns
                    cols_across = max(1, math.ceil(len(node_ids) / depth_rows)) if node_ids else 1
                    row_idx = idx // cols_across
                    col_idx = idx % cols_across
                else:
                    rows = max(1, math.ceil(len(node_ids) / columns)) if node_ids else 1
                    col_idx = idx // rows
                    row_idx = idx % rows
                positions[str(node_id)] = (col_x(col_idx), row_y(row_idx))
                nearest_band = row_idx == 0 if focus_side in {"top", "bottom"} else col_idx == 0
                if nearest_band:
                    arrow_targets.add(str(node_id))
            focus_y_out = center_y - half_h
            if focus_side == "top":
                focus_y_out = focus_y
            elif focus_side == "bottom":
                focus_y_out = focus_y
            return focus_x, focus_y_out, positions, arrow_targets

        def _build_view_for_kinase(
            kinase_id: str, mode: str, ctx: Dict[str, Any], settings_override: Dict[str, Any]
        ) -> Dict[str, Any]:
            ks_data = ctx.get("ks_data") or _empty_ks_index()
            kin_entry = ks_data.get("kinases", {}).get(kinase_id)
            allowed_types = _ks_allowed_types(mode)
            if _ks_dataset_kinases_only() and not _ks_kinase_in_dataset(kinase_id, ctx):
                state["status"].set("Selected kinase is not present in the uploaded protein dataset.")
                return _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            if not kin_entry:
                state["status"].set(f"No kinase data found for {kinase_id}.")
                return _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            sub_keys: List[str] = []
            for sub_key in kin_entry.get("substrates", []):
                sub_entry = ks_data.get("substrates", {}).get(sub_key, {}) or {}
                sub_uid = str(sub_entry.get("uniprot") or "").strip()
                if _ks_dataset_kinases_only() and (not sub_uid or not _ks_protein_in_dataset(sub_uid, ctx)):
                    continue
                rel_types = (
                    sub_entry
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
            layout_cfg = _ks_layout_config()
            mode_layout = layout_cfg.get("mode", "two_column")
            center_x = float(layout_cfg.get("center_x", KS_LAYOUT_CENTER_X_DEFAULT))
            center_y = float(layout_cfg.get("center_y", KS_LAYOUT_CENTER_Y_DEFAULT))
            half_w = protbox_width * 0.5
            half_h = protbox_height * 0.5
            payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            protein_data: Dict[str, Any] = {}
            protboxes: List[Dict[str, Any]] = []
            arrows: List[Dict[str, Any]] = []
            substrate_positions: Dict[str, Tuple[float, float]] = {}
            first_column_targets: set = set()
            kin_pb_id = f"{prefix}_kinase"

            if mode_layout == "concentric":
                kinase_x = center_x - half_w
                kinase_y = center_y - half_h
                radius = _ks_concentric_radius(len(grouped_items), layout_cfg)
                for idx, (sub_uid, _) in enumerate(grouped_items):
                    angle = (2.0 * math.pi * idx) / max(1, len(grouped_items))
                    sub_cx = center_x + radius * math.cos(angle)
                    sub_cy = center_y + radius * math.sin(angle)
                    substrate_positions[sub_uid] = (sub_cx - half_w, sub_cy - half_h)
            else:
                ordered_sub_uids = [sub_uid for (sub_uid, _) in grouped_items]
                kinase_x, kinase_y, substrate_positions, first_column_targets = _ks_two_column_positions(
                    ordered_sub_uids,
                    layout_cfg,
                    center_x,
                    center_y,
                    half_w,
                    half_h,
                    protbox_width,
                )

            kinase_ptms: Dict[str, Any] = {}
            kinase_ptm_rows = ctx.get("ptms_by_uniprot", {}).get(kinase_id, [])
            for idx, row_map in enumerate(kinase_ptm_rows[: len(PTM_POSITION_PRIORITY)]):
                pos_key = PTM_POSITION_PRIORITY[idx]
                ptm_key = f"{kinase_id}_{str(row_map.get(ctx.get('ptm_headers', [None, None])[1] if ctx.get('ptm_headers') else '', '')).strip()}"
                kinase_ptms[ptm_key] = _make_ptm_entry(
                    row_map,
                    ptm_key,
                    pos_key,
                    kinase_x,
                    kinase_y,
                    protbox_width,
                    protbox_height,
                    ctx,
                    settings_override,
                )
            kinase_row = ctx.get("prot_rows_by_uniprot", {}).get(kinase_id, {})
            protein_data[kinase_id] = _make_protein_entry(kinase_id, kinase_row, ctx, settings_override, kinase_ptms)
            kin_label = protein_data[kinase_id].get("label") or kin_entry.get("gene") or kinase_id
            kin_tooltip = protein_data[kinase_id].get("tooltip") or ""
            protboxes.append(
                {
                    "protbox_id": kin_pb_id,
                    "proteins": [kinase_id],
                    "backup_label": kin_label,
                    "x": kinase_x,
                    "y": kinase_y,
                    "width": protbox_width,
                    "height": protbox_height,
                    "tooltip": kin_tooltip,
                }
            )
            for idx, (sub_uid, sub_list) in enumerate(grouped_items):
                sub_x, sub_y = substrate_positions.get(sub_uid, (center_x - half_w, center_y - half_h))
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
                        sub_x,
                        sub_y,
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
                        "x": sub_x,
                        "y": sub_y,
                        "width": protbox_width,
                        "height": protbox_height,
                    }
                )
                if mode_layout == "concentric":
                    if layout_cfg.get("conc_add_arrows"):
                        target_x = sub_x + half_w
                        target_y = sub_y + half_h
                        center_x_node = kinase_x + half_w
                        center_y_node = kinase_y + half_h
                        vx = target_x - center_x_node
                        vy = target_y - center_y_node
                        dist = math.hypot(vx, vy) or 1.0
                        stop_r = float(layout_cfg.get("conc_arrow_stop", KS_CONCENTRIC_ARROW_STOP_DEFAULT))
                        stop = min(stop_r, dist * 0.45) if stop_r else 0.0
                        source_x = center_x_node + (vx / dist) * stop if stop else center_x_node
                        source_y = center_y_node + (vy / dist) * stop if stop else center_y_node
                        end_x = target_x - (vx / dist) * stop if stop else target_x
                        end_y = target_y - (vy / dist) * stop if stop else target_y
                        arrows.append(
                            {
                                "x1": source_x,
                                "y1": source_y,
                                "x2": end_x,
                                "y2": end_y,
                                "line": "arrow",
                                "type": "",
                            }
                        )
                else:
                    if str(sub_uid) not in first_column_targets:
                        continue
                    focus_side = str(layout_cfg.get("two_col_focus_side", "left"))
                    kin_side, sub_side = _ks_two_col_arrow_sides(focus_side, "single_to_multi")
                    arrows.append(
                        {
                            "protbox_id_1": kin_pb_id,
                            "protbox_id_2": pb_id,
                            "protbox_id_1_side": kin_side,
                            "protbox_id_2_side": sub_side,
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
            if _ks_dataset_kinases_only():
                kin_ids = [kin_id for kin_id in kin_ids if _ks_kinase_in_dataset(kin_id, ctx)]
            if not kin_ids:
                state["status"].set("No known kinases found for this PTM with the selected evidence filter.")
                return _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            kin_ids = sorted(set(kin_ids))
            protbox_width = 46
            protbox_height = 17
            layout_cfg = _ks_layout_config()
            mode_layout = layout_cfg.get("mode", "two_column")
            center_x = float(layout_cfg.get("center_x", KS_LAYOUT_CENTER_X_DEFAULT))
            center_y = float(layout_cfg.get("center_y", KS_LAYOUT_CENTER_Y_DEFAULT))
            half_w = protbox_width * 0.5
            half_h = protbox_height * 0.5
            payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
            protein_data: Dict[str, Any] = {}
            protboxes: List[Dict[str, Any]] = []
            arrows: List[Dict[str, Any]] = []
            kinase_positions: Dict[str, Tuple[float, float]] = {}
            first_column_targets: set = set()
            ptm_x = center_x - half_w
            ptm_y = center_y - half_h

            if mode_layout == "concentric":
                radius = _ks_concentric_radius(len(kin_ids), layout_cfg)
                for idx, kin_id in enumerate(kin_ids):
                    angle = (2.0 * math.pi * idx) / max(1, len(kin_ids))
                    kin_cx = center_x + radius * math.cos(angle)
                    kin_cy = center_y + radius * math.sin(angle)
                    kinase_positions[kin_id] = (kin_cx - half_w, kin_cy - half_h)
            else:
                ptm_x, ptm_y, kinase_positions, first_column_targets = _ks_two_column_positions(
                    kin_ids,
                    layout_cfg,
                    center_x,
                    center_y,
                    half_w,
                    half_h,
                    protbox_width,
                )

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
                    ptm_x,
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
                    "x": ptm_x,
                    "y": ptm_y,
                    "width": protbox_width,
                    "height": protbox_height,
                }
            )
            for idx, kin_id in enumerate(kin_ids):
                kin_x, kin_y = kinase_positions.get(kin_id, (center_x - half_w, center_y - half_h))
                kin_ptms: Dict[str, Any] = {}
                kin_ptm_rows = ctx.get("ptms_by_uniprot", {}).get(kin_id, [])
                for ptm_idx, row_map in enumerate(kin_ptm_rows[: len(PTM_POSITION_PRIORITY)]):
                    pos_key = PTM_POSITION_PRIORITY[ptm_idx]
                    kin_ptm_key = f"{kin_id}_{str(row_map.get(ctx.get('ptm_headers', [None, None])[1] if ctx.get('ptm_headers') else '', '')).strip()}"
                    kin_ptms[kin_ptm_key] = _make_ptm_entry(
                        row_map,
                        kin_ptm_key,
                        pos_key,
                        kin_x,
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
                        "x": kin_x,
                        "y": kin_y,
                        "width": protbox_width,
                        "height": protbox_height,
                    }
                )
                if mode_layout == "concentric":
                    if layout_cfg.get("conc_add_arrows"):
                        source_x = kin_x + half_w
                        source_y = kin_y + half_h
                        target_x = ptm_x + half_w
                        target_y = ptm_y + half_h
                        vx = target_x - source_x
                        vy = target_y - source_y
                        dist = math.hypot(vx, vy) or 1.0
                        stop_r = float(layout_cfg.get("conc_arrow_stop", KS_CONCENTRIC_ARROW_STOP_DEFAULT))
                        stop = min(stop_r, dist * 0.45) if stop_r else 0.0
                        source_x = source_x + (vx / dist) * stop if stop else source_x
                        source_y = source_y + (vy / dist) * stop if stop else source_y
                        end_x = target_x - (vx / dist) * stop if stop else target_x
                        end_y = target_y - (vy / dist) * stop if stop else target_y
                        arrows.append(
                            {
                                "x1": source_x,
                                "y1": source_y,
                                "x2": end_x,
                                "y2": end_y,
                                "line": "arrow",
                                "type": "",
                            }
                        )
                else:
                    if str(kin_id) not in first_column_targets:
                        continue
                    focus_side = str(layout_cfg.get("two_col_focus_side", "left"))
                    kin_side, ptm_side = _ks_two_col_arrow_sides(focus_side, "multi_to_single")
                    arrows.append(
                        {
                            "protbox_id_1": kin_pb_id,
                            "protbox_id_2": ptm_pb_id,
                            "protbox_id_1_side": kin_side,
                            "protbox_id_2_side": ptm_side,
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
            settings_override["session_workspace_dir"] = str(_current_session_workspace())
            if DEBUG_FILE_OUTPUT_ENABLED:
                try:
                    settings_override["debug_output_dir"] = str(safe_session_path(session, "debug"))
                except Exception:
                    settings_override["debug_output_dir"] = ""
            web_simple_pathway = (
                cfg.get("key") == "web"
                and str(settings_override.get("pathway_source", "")).lower() in {"kegg", "wikipathways"}
                and bool(settings_override.get("simple_kegg_mode"))
            )
            if web_simple_pathway:
                settings_override["mode"] = "analysis"
                settings_override["show_background_image"] = True
                settings_override["show_text_boxes"] = False
            if cfg.get("key") == "web":
                selected_id = ""
                selected_source = ""
                if web_selected_id is not None:
                    try:
                        selected_id = web_selected_id.get() or ""
                    except Exception:
                        selected_id = ""
                if web_selected_source is not None:
                    try:
                        selected_source = str(web_selected_source.get() or "").strip().lower()
                    except Exception:
                        selected_source = ""
                if selected_id:
                    settings_override["pathway_id"] = selected_id
                if selected_source:
                    settings_override["pathway_source"] = selected_source
                current_id = str(settings_override.get("pathway_id") or "").strip()
                if not selected_id or not current_id or selected_id.strip().lower() != current_id.lower():
                    state["json"].set(None)
                    state["status"].set("Select a pathway from the table, then click Load Pathway.")
                    return
                if str(settings_override.get("pathway_source") or "").strip().lower() == "cst":
                    with reactive.isolate():
                        prot_data = protein_dataset.get()
                        metabolite_data = metabolite_dataset.get()
                        ptm_data = ptm_dataset.get()
                    cst_state_path = _get_cst_session_state_path(selected_id)
                    cst_payload = load_cst_pathway_payload(
                        selected_id,
                        Path(BASE_DIR),
                        protein_dataset=prot_data,
                        metabolite_dataset=metabolite_data,
                        ptm_dataset=ptm_data,
                        ptm_settings={
                            "ptm_max_display": settings_override.get("ptm_max_display", DEFAULT_SETTINGS["ptm_max_display"]),
                            "ptm_label_font": settings_override.get("ptm_label_font", DEFAULT_SETTINGS["ptm_label_font"]),
                            "ptm_label_color": settings_override.get("ptm_label_color", DEFAULT_SETTINGS["ptm_label_color"]),
                            "ptm_label_size": settings_override.get("ptm_label_size", DEFAULT_SETTINGS["ptm_label_size"]),
                            "ptm_circle_radius": settings_override.get("ptm_circle_radius", DEFAULT_SETTINGS["ptm_circle_radius"]),
                            "ptm_circle_spacing": settings_override.get("ptm_circle_spacing", DEFAULT_SETTINGS["ptm_circle_spacing"]),
                            "ptm_outline_width": settings_override.get("ptm_outline_width", DEFAULT_SETTINGS["ptm_outline_width"]),
                            "ptm_shape": (_get_input_value(input, "input_ptm_shape") or "circle"),
                        },
                        negative_color=settings_override.get("negative_color", DEFAULT_SETTINGS["negative_color"]),
                        positive_color=settings_override.get("positive_color", DEFAULT_SETTINGS["positive_color"]),
                        max_negative=settings_override.get("max_negative", DEFAULT_SETTINGS["max_negative"]),
                        max_positive=settings_override.get("max_positive", DEFAULT_SETTINGS["max_positive"]),
                        gradient_stops=settings_override.get("gradient_stops"),
                        prot_outline_width=settings_override.get("prot_outline_width", DEFAULT_SETTINGS["prot_outline_width"]),
                        use_black_protein_outlines=bool(settings_override.get("use_black_protein_outlines", DEFAULT_SETTINGS.get("use_black_protein_outlines", False))),
                        simple_kegg_mode=bool(settings_override.get("simple_kegg_mode", True)),
                        temporal_mode=bool(settings_override.get("temporal_mode", DEFAULT_SETTINGS.get("temporal_mode", False))),
                        overlay_state_path=cst_state_path,
                    )
                    if not cst_payload:
                        state["json"].set(None)
                        state["status"].set("The selected CST pathway file could not be loaded.")
                        return
                    colored_count = len(cst_payload.get("overlay_nodes") or [])
                    payload = {
                        "_viewer_kind": "cst_pdf",
                        "_cst_payload": cst_payload,
                        "_bookmark_key": cfg["key"],
                        "_persist_token": time.time(),
                    }
                    state["json"].set(payload)
                    if colored_count:
                        state["status"].set(
                            f"Loaded CST pathway: {cst_payload.get('name', selected_id)}. Applied protein fold-change colors to {colored_count} CST nodes."
                        )
                    else:
                        state["status"].set(
                            f"Loaded CST pathway: {cst_payload.get('name', selected_id)}. No mapped protein fold-change colors were found for the current dataset."
                        )
                    return
            if settings_override.get("pathway_source") == "wikipathways" and not settings_override.get("pathway_id"):
                state["json"].set(None)
                state["status"].set("Select a WikiPathways entry, then click Load Pathway.")
                return
            if is_ks_bookmark:
                ctx = _ks_context()
                if ctx.get("protein_main_columns"):
                    settings_override["main_columns"] = list(ctx.get("protein_main_columns", []))
                if ctx.get("protein_tooltips"):
                    settings_override["protein_tooltip_columns"] = list(ctx.get("protein_tooltips", []))
                mode = _ks_mode_value()
                selection_input = _get_input_value(input, _prefixed_id(prefix, "ks_choice"))
                if selection_input:
                    _remember_ks_choice(selection_input)
                selection = selection_input or _last_ks_choice()
                selection_str = str(selection or "")

                def _set_ks_payload(next_payload: Dict[str, Any]) -> None:
                    payload_with_token = dict(next_payload or {})
                    payload_with_token["_bookmark_key"] = cfg.get("key", prefix)
                    payload_with_token["_persist_token"] = _ks_stable_persist_token(
                        selection_str,
                        mode,
                        ctx,
                        settings_override,
                    )
                    payload_with_token["_ks_selection"] = selection_str
                    state["json"].set(payload_with_token)

                if not ctx.get("prot_headers") or not ctx.get("ptm_headers"):
                    payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
                    _set_ks_payload(payload)
                    state["status"].set("Upload protein and PTM datasets to use Kinase-Substrate.")
                    return
                if not selection:
                    payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
                    _set_ks_payload(payload)
                    state["status"].set("Select a kinase or substrate, then click Load.")
                    return
                if selection_str.startswith("kinase|"):
                    payload = _build_view_for_kinase(selection_str.split("|", 1)[1], mode, ctx, settings_override)
                elif selection_str.startswith("ptm|"):
                    payload = _build_view_for_ptm(selection_str.split("|", 1)[1], mode, ctx, settings_override)
                else:
                    payload = _prepare_base_payload(settings_override, ctx.get("protein_main_columns", []))
                    state["status"].set("Choose a kinase or substrate option, then click Load.")
                _set_ks_payload(payload)
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
            catalog_info = _current_global_catalog_info()
            color_override = _color_override_from_settings(settings_override)
            if cfg.get("start_blank"):
                payload = _build_blank_canvas(catalog_info)
                payload = dict(payload)
                seed_payload = None
                # Only seed from data override for non-figure modes; figure clear should start empty
                if data_override and cfg.get("key") != "figure":
                    try:
                        seed_payload = get_default_json(
                            data_override=data_override,
                            settings_override=settings_override,
                            skip_disk_write=True,
                            debug_write=DEBUG_FILE_OUTPUT_ENABLED,
                        )
                    except Exception as exc:
                        print(f"Warning: failed to build seed payload for blank canvas: {exc}")
                merged_settings = dict(DEFAULT_SETTINGS)
                merged_settings.update(settings_override)
                merged_settings["pathway_id"] = settings_override.get("pathway_id", DEFAULT_SETTINGS["pathway_id"])
                merged_settings["pathway_source"] = settings_override.get("pathway_source", DEFAULT_SETTINGS["pathway_source"])
                merged_settings["mode"] = settings_override.get("mode", cfg["mode"])
                merged_settings["show_background_image"] = False
                merged_settings["temporal_mode"] = bool(settings_override.get("temporal_mode", DEFAULT_SETTINGS.get("temporal_mode", False)))
                merged_settings["show_groups"] = bool(settings_override.get("show_groups", False))
                merged_settings["show_multi_protein_indicator"] = bool(settings_override.get("show_multi_protein_indicator", False))
                merged_settings["show_arrows"] = bool(settings_override.get("show_arrows", True))
                merged_settings["show_text_boxes"] = bool(settings_override.get("show_text_boxes", True))
                merged_settings["use_original_protbox_size"] = bool(settings_override.get("use_original_protbox_size", DEFAULT_SETTINGS.get("use_original_protbox_size", True)))
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
                merged_settings["use_black_protein_outlines"] = bool(settings_override.get("use_black_protein_outlines", DEFAULT_SETTINGS.get("use_black_protein_outlines", False)))
                merged_settings["main_columns"] = settings_override.get("main_columns", [])
                merged_settings["protein_tooltip_columns"] = settings_override.get("protein_tooltip_columns", DEFAULT_SETTINGS.get("protein_tooltip_columns", []))
                merged_settings["species"] = settings_override.get("species")
                merged_settings["species_code"] = settings_override.get("species_code")
                merged_settings["bg_scale"] = settings_override.get("bg_scale", DEFAULT_WIKIPATHWAYS_BG_SCALE)
                merged_settings["bg_offset_x"] = settings_override.get("bg_offset_x", DEFAULT_WIKIPATHWAYS_BG_OFFSET_X)
                merged_settings["bg_offset_y"] = settings_override.get("bg_offset_y", DEFAULT_WIKIPATHWAYS_BG_OFFSET_Y)
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
                payload["_global_protein_catalog"] = dict(catalog_info)
                payload["_kegg_preview"] = _background_preview_from_settings(settings_override, show=False)
                payload["_box_preview"] = {"y_stretch": DEFAULT_BOX_Y_STRETCH}
                payload["_active_mode"] = settings_override.get("mode", cfg["mode"])
                if cfg["key"] in {"figure", "web"}:
                    payload["_full_width_canvas"] = True
                layout_override = _get_custom_layout()
                if layout_override:
                    _apply_custom_layout(payload, layout_override, include_missing=(cfg["key"] == "figure"))
                    if cfg["key"] == "figure":
                        match_stats = _hydrate_custom_layout_from_current_dataset(payload, _ks_context(), settings_override)
                        payload["_custom_layout_dataset_match_stats"] = match_stats
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
                    debug_write=DEBUG_FILE_OUTPUT_ENABLED,
                )
                generation_error = str(payload.get("_generation_error") or "").strip() if isinstance(payload, dict) else ""
                if generation_error:
                    if prefix == "web":
                        _set_load_feedback("Pathway generation failed.")
                    state["status"].set("Error generating pathway data.")
                    return
                show_bg_flag = bool(settings_override.get("show_background_image", False))
                payload, _ = _attach_pathway_background_image(payload, force=show_bg_flag)
                payload = dict(payload)
                payload["_bookmark_key"] = cfg["key"]
                payload["_persist_token"] = time.time()
                payload["_global_protein_catalog"] = dict(catalog_info)
                _apply_color_metadata(payload, color_override)
                payload["_kegg_preview"] = _background_preview_from_settings(settings_override, show=show_bg_flag)
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
                if web_simple_pathway:
                    payload = _apply_simple_pathway_payload(payload)
                    payload["_bookmark_key"] = cfg["key"]
                    payload["_persist_token"] = time.time()
                    payload["_global_protein_catalog"] = dict(catalog_info)
                    payload["_kegg_preview"] = _background_preview_from_settings(settings_override, show=True)
                layout_override = _get_custom_layout()
                if layout_override:
                    _apply_custom_layout(payload, layout_override, include_missing=(cfg["key"] == "figure"))
                    if cfg["key"] == "figure":
                        match_stats = _hydrate_custom_layout_from_current_dataset(payload, _ks_context(), settings_override)
                        payload["_custom_layout_dataset_match_stats"] = match_stats
                    payload["_custom_layout_applied"] = True
                else:
                    payload["_custom_layout_applied"] = False
                _write_session_preview_json(payload)
                state["json"].set(payload)
                if prefix == "web":
                    _set_load_feedback("")
                state["status"].set(f"{cfg['label']} pathway generated.")
            except Exception as exc:
                if prefix == "web":
                    _set_load_feedback("Pathway generation failed.")
                print(f"Warning: pathway JSON generation failed: {exc}")
                state["status"].set("Error generating pathway data.")

        if is_ks_bookmark:
            def _ks_table_rows() -> List[Dict[str, Any]]:
                ctx = _ks_context()
                ks_data = ctx.get("ks_data") or _empty_ks_index()
                entity_mode = str(_get_input_value(input, _prefixed_id(prefix, "ks_entity_mode")) or "substrate").lower()
                allowed_types = _ks_allowed_types(_ks_mode_value())
                dataset_kinase_only = _ks_dataset_kinases_only()
                if not allowed_types:
                    return []
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
                fc_op = str(_get_input_value(input, _prefixed_id(prefix, "ks_filter_fc_op")) or "")
                fc_val_raw = _get_input_value(input, _prefixed_id(prefix, "ks_filter_fc_val"))
                try:
                    fc_val = float(fc_val_raw) if fc_val_raw not in (None, "", False) else None
                except (TypeError, ValueError):
                    fc_val = None
                reg_only = bool(_get_input_value(input, _prefixed_id(prefix, "ks_filter_reg_only")))
                def _is_sig(fc: Optional[float]) -> bool:
                    return fc is not None and (fc >= max_pos or fc <= max_neg)

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
                        if dataset_kinase_only and not _ks_kinase_in_dataset(kin_id, ctx):
                            continue
                        gene = _resolve_gene_label(kin_id, info, ctx)
                        prot_row = ctx.get("prot_rows_by_uniprot", {}).get(kin_id, {})
                        prot_fc_cols = [selected_fc] if selected_fc else ctx.get("protein_main_columns", [])
                        fc_val = _first_fc_value(prot_row, prot_fc_cols)
                        sig_count = 0
                        sub_count = 0
                        visible_types: set = set()
                        for sub_key in info.get("substrates", []):
                            sub_entry = (ks_data.get("substrates") or {}).get(sub_key, {})
                            sub_uid = str(sub_entry.get("uniprot") or "").strip()
                            if dataset_kinase_only and (not sub_uid or not _ks_protein_in_dataset(sub_uid, ctx)):
                                continue
                            rel_types = (sub_entry.get("kinases", {}).get(kin_id, {}) or {}).get("types", set())
                            matched_types = set(rel_types or set()).intersection(allowed_types)
                            if not matched_types:
                                continue
                            visible_types.update(matched_types)
                            sub_count += 1
                            sub_fc = _first_fc_value(sub_entry.get("row", {}), ctx.get("ptm_main_columns", []))
                            if _is_sig(sub_fc):
                                sig_count += 1
                        if not sub_count:
                            continue
                        evid_label = _ks_evidence_label(visible_types)
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
                        sub_uid = info.get("uniprot", "")
                        gene = _resolve_gene_label(sub_uid, info, ctx)
                        site = info.get("site_label") or info.get("site") or ""
                        ptm_fc_cols = []
                        if selected_fc and selected_fc in (ctx.get("ptm_main_columns") or []):
                            ptm_fc_cols = [selected_fc]
                        else:
                            ptm_fc_cols = ctx.get("ptm_main_columns", [])
                        ptm_fc_val = _first_fc_value(info.get("row", {}), ptm_fc_cols)
                        kin_count = 0
                        visible_types: set = set()
                        for kin_id, rel in (info.get("kinases") or {}).items():
                            if dataset_kinase_only and not _ks_kinase_in_dataset(kin_id, ctx):
                                continue
                            matched_types = set((rel or {}).get("types", set()) or set()).intersection(allowed_types)
                            if not matched_types:
                                continue
                            kin_count += 1
                            visible_types.update(matched_types)
                        if not kin_count:
                            continue
                        reg_val = str(info.get("row", {}).get("PSP: regulatory_site", "")).strip()
                        evid_label = _ks_evidence_label(visible_types)
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
                        header("Significant PTMs (%)", "sig_pct"),
                        header("Substrates", "count"),
                    ]
                else:
                    headers = [
                        header("Protein", "gene"),
                        header("Site", "site"),
                        header("FC", "fc"),
                        header("Regulatory site", "reg"),
                        header("Kinases", "count"),
                    ]
                body_rows = []
                selected_id = str(_get_input_value(input, _prefixed_id(prefix, "ks_choice")) or "")
                if not selected_id:
                    selected_id = _last_ks_choice()
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
                            ui.tags.td("N/A" if entry.get("sig_pct") is None else f"{entry['sig_pct']:.1f}%"),
                            ui.tags.td(str(entry.get("count", 0))),
                        ]
                    else:
                        cells = [
                            ui.tags.td(entry.get("gene", "")),
                            ui.tags.td(entry.get("site", "")),
                            ui.tags.td("" if entry.get("fc") is None else f"{entry['fc']:.3g}"),
                            ui.tags.td(entry.get("reg", "")),
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
                session.send_input_message(_prefixed_id(prefix, "ks_filter_fc_op"), {"value": ""})
                session.send_input_message(_prefixed_id(prefix, "ks_filter_fc_val"), {"value": ""})

            @reactive.Effect
            @reactive.event(getattr(input, _prefixed_id(prefix, "ks_choice")))
            def _build_kinase_substrate_view():
                choice = _get_input_value(input, _prefixed_id(prefix, "ks_choice"))
                if choice:
                    _remember_ks_choice(choice)
                    build_json()

            @reactive.Effect
            @reactive.event(getattr(input, _prefixed_id(prefix, "ks_evidence_types")))
            def _sync_ks_evidence_types():
                current_choice = _get_input_value(input, _prefixed_id(prefix, "ks_choice")) or _last_ks_choice()
                if current_choice:
                    build_json()

            @reactive.Effect
            @reactive.event(
                getattr(input, _prefixed_id(prefix, "ks_layout_mode")),
                getattr(input, _prefixed_id(prefix, "ks_layout_gap_x")),
                getattr(input, _prefixed_id(prefix, "ks_layout_spacing_y")),
                getattr(input, _prefixed_id(prefix, "ks_two_col_focus_side")),
                getattr(input, _prefixed_id(prefix, "ks_two_col_columns")),
                getattr(input, _prefixed_id(prefix, "ks_conc_arrows")),
                getattr(input, _prefixed_id(prefix, "ks_conc_radius_mode")),
                getattr(input, _prefixed_id(prefix, "ks_conc_radius_fixed")),
                getattr(input, _prefixed_id(prefix, "ks_conc_space")),
                getattr(input, _prefixed_id(prefix, "ks_conc_arrow_stop")),
                getattr(input, _prefixed_id(prefix, "ks_dataset_kinases_only")),
            )
            def _sync_ks_layout_and_filters():
                current_choice = _get_input_value(input, _prefixed_id(prefix, "ks_choice")) or _last_ks_choice()
                if current_choice:
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
                if _ks_dataset_kinases_only():
                    kin_ids = [kin_id for kin_id in kin_ids if _ks_kinase_in_dataset(kin_id, ctx)]
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
                    species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
                    species_code = species_info["code"]
                    kegg_options = _get_kegg_pathway_options_for_species(species_code)
                    if not kegg_options:
                        return ui.div({"class": "pathway-search-empty"}, "No KEGG catalogue found.")
                    if not query:
                        # Show all options when empty; list is scrollable in UI.
                        for opt in kegg_options:
                            species_id = f"{species_code}{opt['digits']}"
                            matches.append({"value": species_id, "label": f"{species_id} | {opt['name']}"})
                    else:
                        try:
                            pattern = re.compile(query, re.IGNORECASE)
                        except re.error:
                            return ui.div({"class": "pathway-search-error"}, "Invalid regex pattern")
                        for opt in kegg_options:
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
                    if cfg.get("key") == "web":
                        if web_selected_id is not None:
                            web_selected_id.set(str(value))
                        if web_selected_source is not None:
                            web_selected_source.set(str(source_val or "wikipathways").strip().lower())
                    else:
                        session.send_input_message(_prefixed_id(prefix, "pathway_id"), {"value": value})
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
                    _set_load_feedback("")
                state["status"].set("Pathway selected. Click Load Pathway to load with current settings.")

            if cfg.get("key") == "web":
                @output(id=_prefixed_id(prefix, "generate_feedback"))
                @render.ui
                def generate_feedback():
                    feedback = str(state.get("load_feedback").get() or "").strip()
                    loading = bool(state.get("loading").get()) if state.get("loading") is not None else False
                    if loading or not feedback:
                        return None
                    return ui.tags.span(
                        {
                            "class": "web-load-error",
                            "aria-live": "polite",
                        },
                        feedback,
                    )

                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_id")))
                def _sync_web_selection_input():
                    return

                @output(id=_prefixed_id(prefix, "generate_button_state"))
                @render.ui
                def generate_button_state():
                    selected = (web_selected_id.get() if web_selected_id is not None else "") or ""
                    loading = bool(state.get("loading").get()) if state.get("loading") is not None else False
                    class_name = "web-load-ready" if str(selected).strip() else ""
                    return ui.tags.script(
                        f"""
                        (function(){{
                            const btn = document.getElementById('{_prefixed_id(prefix, "generate")}');
                            const spinner = document.getElementById('{_prefixed_id(prefix, "generate_spinner")}');
                            if (!btn || !spinner) return;
                            if (btn.dataset.mkLoadBound !== '1') {{
                                btn.addEventListener('click', function(){{
                                    if (btn.disabled) return;
                                    spinner.classList.add('active');
                                    btn.disabled = true;
                                }});
                                btn.dataset.mkLoadBound = '1';
                            }}
                            btn.classList.toggle('web-load-ready', {str(bool(class_name)).lower()});
                            btn.disabled = {str(loading).lower()};
                            spinner.classList.toggle('active', {str(loading).lower()});
                        }})();
                        """
                    )

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "run_fisher_test")))
                def _handle_run_fisher_test():
                    fc_choices = _fc_choices()
                    fc_idx = state["fc_index"].get() or 1
                    if fc_idx < 1 or fc_idx > len(fc_choices):
                        fc_idx = 1
                    selected_fc = fc_choices[fc_idx - 1] if fc_choices and fc_idx >= 1 else None
                    _refresh_pathway_scores(selected_fc)

            if cfg.get("key") == "web":
                @output(id=_prefixed_id(prefix, "fisher_run_state"))
                @render.ui
                def fisher_run_state():
                    score_cache = pathway_score_cache.get() or {}
                    _ = score_cache.get("updated_at")
                    return ui.tags.script(
                        f"""
                        (function(){{
                            const runBtn = document.getElementById('{_prefixed_id(prefix, "run_fisher_test")}');
                            const progress = document.getElementById('{_prefixed_id(prefix, "fisher_run_progress")}');
                            if (runBtn) runBtn.disabled = false;
                            if (progress) progress.classList.remove('active');
                        }})();
                        """
                    )

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(
                    getattr(input, _prefixed_id(prefix, "fisher_sig_mode")),
                    getattr(input, _prefixed_id(prefix, "fisher_sig_pos")),
                    getattr(input, _prefixed_id(prefix, "fisher_sig_neg")),
                )
                def _mark_fisher_scores_stale():
                    with reactive.isolate():
                        protein_ok = bool(protein_validation.get().get("valid"))
                        score_cache = pathway_score_cache.get() or {}
                    if protein_ok:
                        cached_results = score_cache.get("results_by_fc", {}) if isinstance(score_cache, dict) else {}
                        if not cached_results:
                            return
                        cached_mode = str(score_cache.get("significance_mode") or "both").strip().lower()
                        cached_pos = abs(_to_float(score_cache.get("positive_cutoff"), 1.5))
                        cached_neg = -abs(_to_float(score_cache.get("negative_cutoff"), -1.5))
                        current_mode = str(_get_input_value(input, _prefixed_id(prefix, "fisher_sig_mode")) or "both").strip().lower()
                        current_pos = abs(_to_float(_get_input_value(input, _prefixed_id(prefix, "fisher_sig_pos")), 1.5))
                        current_neg = -abs(_to_float(_get_input_value(input, _prefixed_id(prefix, "fisher_sig_neg")), -1.5))
                        if (
                            current_mode == cached_mode
                            and abs(current_pos - cached_pos) < 1e-9
                            and abs(current_neg - cached_neg) < 1e-9
                        ):
                            return
                        _mark_pathway_scores_pending("Fisher settings changed. Run Fisher's Exact Test to refresh pathway scores.")

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_table_sort")))
                def _handle_web_sort():
                    payload = _get_input_value(input, _prefixed_id(prefix, "pathway_table_sort")) or {}
                    col = str(payload.get("col") or "").lower()
                    if col not in {
                        "source",
                        "pathway",
                        "prot_fisher_p",
                        "phos_fisher_p",
                    }:
                        return
                    numeric_cols = {
                        "prot_fisher_p",
                        "phos_fisher_p",
                    }
                    current_col = web_sort_col.get()
                    current_dir = web_sort_dir.get()
                    if col == current_col:
                        web_sort_dir.set("desc" if current_dir == "asc" else "asc")
                    else:
                        web_sort_col.set(col)
                        web_sort_dir.set("desc" if col in numeric_cols else "asc")
                    if web_page is not None:
                        web_page.set(1)

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
                    if web_page is not None:
                        web_page.set(1)

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_id")))
                def _reset_web_page_for_search():
                    if web_page is not None:
                        web_page.set(1)

            if cfg.get("key") == "web":
                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "pathway_table_page_nav")))
                def _handle_web_page_nav():
                    if web_page is None:
                        return
                    payload = _get_input_value(input, _prefixed_id(prefix, "pathway_table_page_nav")) or {}
                    action = str(payload.get("action") or "").strip().lower()
                    current_page = max(1, int(web_page.get() or 1))
                    if action == "prev":
                        web_page.set(max(1, current_page - 1))
                    elif action == "next":
                        web_page.set(current_page + 1)

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
                    if web_page is not None:
                        web_page.set(1)

            if cfg.get("key") == "web":
                @output(id=_prefixed_id(prefix, "pathway_table"))
                @render.ui
                def pathway_table():
                    if web_filter_tick is not None:
                        # Depend on tick so table refreshes when checkboxes change
                        _ = web_filter_tick.get()
                    fc_choices = _fc_choices()
                    fc_idx = state["fc_index"].get() or 1
                    if fc_idx < 1 or fc_idx > len(fc_choices):
                        fc_idx = 1
                    selected_fc = fc_choices[fc_idx - 1] if fc_choices and fc_idx >= 1 else ""
                    score_cache = pathway_score_cache.get() or {}
                    results_by_fc = score_cache.get("results_by_fc", {}) if isinstance(score_cache, dict) else {}
                    selected_score_bundle = results_by_fc.get(selected_fc, {}) if isinstance(results_by_fc, dict) else {}
                    if not selected_score_bundle and isinstance(results_by_fc, dict):
                        # Fallback to first available scoring result
                        for _, value in results_by_fc.items():
                            if isinstance(value, dict):
                                selected_score_bundle = value
                                break
                    source_score_maps: Dict[str, Dict[str, Dict[str, Any]]] = {"kegg": {}, "wikipathways": {}, "cst": {}}
                    if isinstance(selected_score_bundle, dict):
                        for source_key in ("kegg", "wikipathways", "cst"):
                            source_obj = selected_score_bundle.get(source_key)
                            if isinstance(source_obj, dict):
                                source_score_maps[source_key] = source_obj
                        # Backward compatibility: legacy flat KEGG-only map.
                        if not source_score_maps["kegg"] and not source_score_maps["wikipathways"] and not source_score_maps["cst"]:
                            maybe_pathway_id = next(iter(selected_score_bundle.keys()), "")
                            maybe_row = selected_score_bundle.get(maybe_pathway_id) if maybe_pathway_id else None
                            if isinstance(maybe_row, dict) and "final_score" in maybe_row:
                                source_score_maps["kegg"] = selected_score_bundle
                    species_key, species_info = _resolve_species(_get_input_value(input, "input_species"))
                    species_label = species_info.get("label") or species_key
                    species_full = species_info.get("species") or species_label
                    options = _load_wikipathways_catalog(species_label, fallback=species_full)
                    cst_options = get_cst_pathway_catalog(Path(BASE_DIR))
                    species_code = species_info.get("code", "")
                    kegg_options = _get_kegg_pathway_options_for_species(species_code)
                    kegg_rows: List[Dict[str, Any]] = []
                    for opt in kegg_options:
                        species_id = f"{species_code}{opt['digits']}"
                        kegg_rows.append(
                            {
                                "id": species_id,
                                "pathway": opt["name"],
                                "source": "KEGG",
                                "prot_fisher_p": None,
                                "phos_fisher_p": None,
                            }
                        )
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
                                "prot_fisher_p": None,
                                "phos_fisher_p": None,
                            }
                            for opt in options
                            if opt.get("id") and opt.get("name")
                        ]
                    cst_rows: List[Dict[str, Any]] = [
                        {
                            "id": str(opt.get("id") or "").strip(),
                            "pathway": str(opt.get("name") or "").strip(),
                            "source": "CST",
                            "has_save_file": bool(opt.get("has_save_file")),
                            "has_saved_edges": bool(opt.get("has_saved_edges")),
                            "has_saved_textboxes": bool(opt.get("has_saved_textboxes")),
                            "prot_fisher_p": None,
                            "phos_fisher_p": None,
                        }
                        for opt in cst_options
                        if opt.get("id") and opt.get("name")
                    ]
                    query = (str(_get_input_value(input, _prefixed_id(prefix, "pathway_id"))) or "").strip()
                    combined_rows: List[Dict[str, Any]] = []
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
                        combined_rows.append(dict(opt))
                    for opt in kegg_rows:
                        search_target = f"{opt['id']} | {opt['pathway']} | {opt['source']}"
                        if pattern and not pattern.search(search_target):
                            continue
                        combined_rows.append(opt)
                    for opt in cst_rows:
                        search_target = f"{opt['id']} | {opt['pathway']} | {opt['source']}"
                        if pattern and not pattern.search(search_target):
                            continue
                        combined_rows.append(dict(opt))
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

                    if isinstance(source_score_maps, dict):
                        for row in combined_rows:
                            source_key = str(row.get("source", "")).strip().lower()
                            pid = str(row.get("id", "")).strip().lower()
                            source_map = source_score_maps.get(source_key, {})
                            score_entry = source_map.get(pid) if isinstance(source_map, dict) else None
                            if not isinstance(score_entry, dict):
                                continue
                            row["prot_fisher_p"] = score_entry.get("prot_fisher_p")
                            row["phos_fisher_p"] = score_entry.get("phos_fisher_p")
                    if not combined_rows:
                        return ui.div({"class": "pathway-search-empty"}, "No pathways matched your search.")
                    sort_col = (web_sort_col.get() or "prot_fisher_p") if web_sort_col else "prot_fisher_p"
                    sort_dir = (web_sort_dir.get() or "asc") if web_sort_dir else "asc"
                    sort_field_map = {
                        "source": "source",
                        "pathway": "pathway",
                        "prot_fisher_p": "prot_fisher_p",
                        "phos_fisher_p": "phos_fisher_p",
                    }
                    field = sort_field_map.get(sort_col, "pathway")
                    numeric_fields = {
                        "prot_fisher_p",
                        "phos_fisher_p",
                    }
                    if field in numeric_fields:
                        ascending = sort_dir == "asc"

                        def _numeric_sort_key(row: Dict[str, Any]) -> Tuple[int, float, str]:
                            raw_val = row.get(field)
                            parsed_val: Optional[float] = None
                            if raw_val not in (None, ""):
                                try:
                                    parsed_val = float(raw_val)
                                except (TypeError, ValueError):
                                    parsed_val = None
                            if parsed_val is None:
                                # Missing values always sort after real pathway statistics.
                                missing_rank = 1
                                sort_val = 0.0
                            else:
                                missing_rank = 0
                                sort_val = parsed_val if ascending else -parsed_val
                            return (missing_rank, sort_val, str(row.get("pathway", "")).lower())

                        combined_rows.sort(key=_numeric_sort_key)
                    else:
                        reverse = sort_dir == "desc"
                        combined_rows.sort(
                            key=lambda r: (
                                str(r.get(field, "")).lower(),
                                str(r.get("pathway", "")).lower(),
                            ),
                            reverse=reverse,
                        )
                    rows_per_page = 20
                    total_rows = len(combined_rows)
                    total_pages = max(1, math.ceil(total_rows / rows_per_page))
                    current_page = max(1, int(web_page.get() or 1)) if web_page is not None else 1
                    current_page = min(current_page, total_pages)
                    page_start = (current_page - 1) * rows_per_page
                    page_end = page_start + rows_per_page
                    page_rows = combined_rows[page_start:page_end]
                    active_value = ""
                    def _fmt_metric(raw: Any, digits: int = 4) -> str:
                        if raw is None or raw == "":
                            return ""
                        try:
                            return f"{float(raw):.{digits}g}"
                        except (TypeError, ValueError):
                            return str(raw)

                    table_rows: List[Any] = []
                    def _metric_title(entry: Dict[str, Any], prefix_name: str) -> str:
                        label = "Protein" if prefix_name == "prot" else "Phospho"
                        fisher_val = entry.get(f"{prefix_name}_fisher_p")
                        dataset_total = entry.get(f"{prefix_name}_dataset_total")
                        pathway_total = entry.get(f"{prefix_name}_pathway_total")
                        sig_total = entry.get(f"{prefix_name}_significant_total")
                        sig_path = entry.get(f"{prefix_name}_significant_in_pathway")
                        return (
                            f"{label} Fisher p: {_fmt_metric(fisher_val, 6) or 'NA'}\n"
                            f"Dataset total: {dataset_total if dataset_total is not None else 'NA'}\n"
                            f"In pathway: {pathway_total if pathway_total is not None else 'NA'}\n"
                            f"Significant total: {sig_total if sig_total is not None else 'NA'}\n"
                            f"Significant in pathway: {sig_path if sig_path is not None else 'NA'}"
                        )

                    for entry in page_rows:
                        classes = []
                        if active_value and active_value == entry["id"]:
                            classes.append("table-active")
                        table_rows.append(
                            ui.tags.tr(
                                {
                                    "data-value": entry["id"],
                                    "data-source": entry.get("source", ""),
                                    "data-name": entry.get("pathway", ""),
                                    "data-pathway-full-name": entry.get("pathway", ""),
                                    "class": " ".join(classes),
                                    "onclick": (
                                        "(() => {"
                                        "const row=this;"
                                        "const table=row.closest('table');"
                                        "if(table){table.querySelectorAll('tbody tr.table-active').forEach((el)=>el.classList.remove('table-active'));}"
                                        "row.classList.add('table-active');"
                                        "const data={value:this.dataset.value,source:this.dataset.source||'wikipathways',name:this.dataset.name||''};"
                                        f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_id_choice')}', data, {{priority:'event'}});"
                                        f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_source_choice')}', data.source, {{priority:'event'}});"
                                        "})();"
                                    ),
                                },
                                ui.tags.td(entry["source"]),
                                ui.tags.td(
                                    ui.tags.div(
                                        {},
                                        entry["pathway"],
                                    ),
                                    ui.tags.small({"style": "color:#6b7280;"}, entry["id"]),
                                ),
                                ui.tags.td({"title": _metric_title(entry, "prot")}, _fmt_metric(entry.get("prot_fisher_p"), 5)),
                                ui.tags.td({"title": _metric_title(entry, "phos")}, _fmt_metric(entry.get("phos_fisher_p"), 5)),
                            )
                        )
                    source_arrow = "▲" if sort_col == "source" and sort_dir == "asc" else ("▼" if sort_col == "source" else "")
                    pathway_arrow = "▲" if sort_col == "pathway" and sort_dir == "asc" else ("▼" if sort_col == "pathway" else "")
                    prot_arrow = "▲" if sort_col == "prot_fisher_p" and sort_dir == "asc" else ("▼" if sort_col == "prot_fisher_p" else "")
                    phos_arrow = "▲" if sort_col == "phos_fisher_p" and sort_dir == "asc" else ("▼" if sort_col == "phos_fisher_p" else "")
                    page_nav_input_id = _prefixed_id(prefix, "pathway_table_page_nav")
                    table_id = _prefixed_id(prefix, "pathway_table_element")

                    def _pagination_controls(location: str) -> Any:
                        prev_attrs: Dict[str, Any] = {
                            "type": "button",
                            "class": "btn btn-outline-secondary btn-sm",
                            "onclick": (
                                f"Shiny.setInputValue('{page_nav_input_id}', "
                                "{action:'prev', ts:Date.now()}, {priority:'event'});"
                            ),
                        }
                        next_attrs: Dict[str, Any] = {
                            "type": "button",
                            "class": "btn btn-outline-secondary btn-sm",
                            "onclick": (
                                f"Shiny.setInputValue('{page_nav_input_id}', "
                                "{action:'next', ts:Date.now()}, {priority:'event'});"
                            ),
                        }
                        if current_page <= 1:
                            prev_attrs["disabled"] = True
                        if current_page >= total_pages:
                            next_attrs["disabled"] = True
                        return ui.div(
                            {
                                "class": f"pathway-table-pagination pathway-table-pagination-{location}",
                                "style": "display:flex; align-items:center; justify-content:space-between; gap:12px; margin:8px 0;",
                            },
                            ui.tags.div(
                                {"style": "color:#4b5563; font-size:0.9rem;"},
                                f"Showing {page_start + 1}-{min(page_end, total_rows)} of {total_rows}",
                            ),
                            ui.tags.div(
                                {"style": "display:flex; align-items:center; gap:8px;"},
                                ui.tags.button("Previous", prev_attrs),
                                ui.tags.span(
                                    {"style": "min-width:92px; text-align:center; color:#374151; font-size:0.9rem;"},
                                    f"Page {current_page} / {total_pages}",
                                ),
                                ui.tags.button("Next", next_attrs),
                            ),
                        )

                    table = ui.tags.table(
                        {"id": table_id, "class": "table table-sm table-hover pathway-table"},
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
                                ui.tags.th(
                                    {
                                        "style": "cursor:pointer;",
                                        "title": "Fisher exact test p-value for protein-level pathway enrichment.",
                                        "onclick": f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_table_sort')}', {{col:'prot_fisher_p', ts:Date.now()}}, {{priority:'event'}});",
                                    },
                                    ui.tags.span("Prot p-value"),
                                    ui.tags.span(f" {prot_arrow}" if prot_arrow else ""),
                                ),
                                ui.tags.th(
                                    {
                                        "style": "cursor:pointer;",
                                        "title": "Fisher exact test p-value for phosphosite-level pathway enrichment.",
                                        "onclick": f"Shiny.setInputValue('{_prefixed_id(prefix, 'pathway_table_sort')}', {{col:'phos_fisher_p', ts:Date.now()}}, {{priority:'event'}});",
                                    },
                                    ui.tags.span("PTM p-value"),
                                    ui.tags.span(f" {phos_arrow}" if phos_arrow else ""),
                                ),
                            )
                        ),
                        ui.tags.tbody(*table_rows),
                    )
                    return ui.TagList(
                        _pagination_controls("top"),
                        table,
                        _pagination_controls("bottom"),
                        ui.tags.script(
                            f"""
                            (function(){{
                                const table = document.getElementById('{table_id}');
                                if(!table) return;
                                const tooltipId = '{table_id}_hover_tooltip';
                                let tooltip = document.getElementById(tooltipId);
                                if(!tooltip){{
                                    tooltip = document.createElement('div');
                                    tooltip.id = tooltipId;
                                    tooltip.className = 'pathway-table-hover-tooltip';
                                    document.body.appendChild(tooltip);
                                }}
                                let hoverTimer = null;
                                let activeRow = null;
                                let lastX = 0;
                                let lastY = 0;

                                const hideTooltip = () => {{
                                    if (hoverTimer) {{
                                        clearTimeout(hoverTimer);
                                        hoverTimer = null;
                                    }}
                                    activeRow = null;
                                    tooltip.style.display = 'none';
                                }};

                                const positionTooltip = (x, y) => {{
                                    const margin = 14;
                                    const maxLeft = window.innerWidth - tooltip.offsetWidth - margin;
                                    const maxTop = window.innerHeight - tooltip.offsetHeight - margin;
                                    const left = Math.max(margin, Math.min(x + 16, maxLeft));
                                    const top = Math.max(margin, Math.min(y + 18, maxTop));
                                    tooltip.style.left = left + 'px';
                                    tooltip.style.top = top + 'px';
                                }};

                                const scheduleTooltip = (row) => {{
                                    hideTooltip();
                                    if (!row) return;
                                    const fullName = row.getAttribute('data-pathway-full-name') || '';
                                    if (!fullName) return;
                                    activeRow = row;
                                    hoverTimer = window.setTimeout(() => {{
                                        if (activeRow !== row) return;
                                        tooltip.textContent = fullName;
                                        tooltip.style.display = 'block';
                                        positionTooltip(lastX, lastY);
                                    }}, 500);
                                }};

                                table.addEventListener('click', (ev) => {{
                                    const row = ev.target.closest('tr');
                                    if(!row || !row.dataset.value) return;
                                    const rows = table.querySelectorAll('tr');
                                    rows.forEach(r => r.classList.remove('table-active'));
                                    row.classList.add('table-active');
                                    hideTooltip();
                                }});
                                table.addEventListener('mousemove', (ev) => {{
                                    lastX = ev.clientX;
                                    lastY = ev.clientY;
                                    const row = ev.target.closest('tbody tr[data-pathway-full-name]');
                                    if (row !== activeRow) {{
                                        scheduleTooltip(row);
                                        return;
                                    }}
                                    if (tooltip.style.display === 'block') {{
                                        positionTooltip(lastX, lastY);
                                    }}
                                }});
                                table.addEventListener('mouseleave', hideTooltip);
                                table.addEventListener('mousedown', hideTooltip);
                                window.addEventListener('scroll', hideTooltip, true);
                                window.addEventListener('resize', hideTooltip);
                                document.addEventListener('keydown', hideTooltip);
                                document.addEventListener('mousedown', (ev) => {{
                                    if (!table.contains(ev.target)) hideTooltip();
                                }});
                            }})();
                            """
                        ),
                    )

                @reactive.Effect
                @reactive.event(getattr(input, _prefixed_id(prefix, "download_pathway_table")))
                def _download_pathway_table():
                    fc_choices = _fc_choices()
                    fc_idx = state["fc_index"].get() or 1
                    if fc_idx < 1 or fc_idx > len(fc_choices):
                        fc_idx = 1
                    selected_fc = fc_choices[fc_idx - 1] if fc_choices and fc_idx >= 1 else ""
                    score_cache = pathway_score_cache.get() or {}
                    download_rows_by_fc = score_cache.get("download_rows_by_fc", {}) if isinstance(score_cache, dict) else {}
                    rows = download_rows_by_fc.get(selected_fc, []) if isinstance(download_rows_by_fc, dict) else []
                    if not rows and isinstance(download_rows_by_fc, dict):
                        for value in download_rows_by_fc.values():
                            if isinstance(value, list) and value:
                                rows = value
                                break
                    columns = [
                        "pathway_source",
                        "pathway_id",
                        "name",
                        "comparison_column",
                        "significance_mode",
                        "positive_cutoff",
                        "negative_cutoff",
                        "prot_dataset_total",
                        "prot_pathway_total",
                        "prot_significant_total",
                        "prot_significant_in_pathway",
                        "prot_fisher_p",
                        "phos_dataset_total",
                        "phos_pathway_total",
                        "phos_significant_total",
                        "phos_significant_in_pathway",
                        "phos_fisher_p",
                    ]
                    buffer = io.StringIO()
                    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore")
                    writer.writeheader()
                    for row in rows:
                        safe_row = {
                            key: _neutralize_spreadsheet_formula_cell(row.get(key))
                            for key in columns
                        }
                        writer.writerow(safe_row)
                    _send_custom_message_safe(
                        "download_payload",
                        {
                            "filename": "web_pathway_table.csv",
                            "content": buffer.getvalue(),
                            "mime_type": "text/csv;charset=utf-8",
                        },
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
            style, labels = _gradient_preview_style_and_labels()
            if not labels:
                labels = ["-2", "0", "2"]
            return ui.div(
                {"class": "gradient-preview", "style": style},
                *[ui.tags.span(label) for label in labels],
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
            # Keep the existing persist token so client-side edits remain cached.
            # Only update the active fold-change index when it actually changes.
            current = state["json"].get()
            if current:
                current_active_idx = current.get("_active_fc_index")
                try:
                    current_active_idx = int(current_active_idx)
                except Exception:
                    current_active_idx = None
                if current_active_idx != idx:
                    updated = dict(current)
                    updated["_active_fc_index"] = idx
                    state["json"].set(updated)
            if (
                cfg.get("key") == "web"
                and bool(pipeline_ready.get())
                and choices
            ):
                selected_fc = choices[idx - 1]
                score_cache = pathway_score_cache.get() or {}
                cached_results = score_cache.get("results_by_fc", {}) if isinstance(score_cache, dict) else {}
                if selected_fc not in cached_results:
                    _refresh_pathway_scores(selected_fc)

        @reactive.Effect
        @reactive.event(
            input.settings_negative_color,
            input.settings_positive_color,
            input.settings_max_negative,
            input.settings_max_positive,
            input.settings_gradient_stops_json,
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
                if is_ks_bookmark:
                    ctx = _ks_context()
                    mode = _ks_mode_value()
                    selection = str(_get_input_value(input, _prefixed_id(prefix, "ks_choice")) or "")
                    payload["_bookmark_key"] = cfg.get("key", prefix)
                    payload["_persist_token"] = _ks_stable_persist_token(selection, mode, ctx, settings_override)
                else:
                    payload["_persist_token"] = time.time()
                state["json"].set(payload)
            except Exception as exc:
                print(f"Warning: failed to sync colors from settings for bookmark {prefix}: {exc}")

        if DEBUG_UI_ENABLED:
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
                    if is_ks_bookmark:
                        settings_override = collect_settings(input, cfg)
                        ctx = _ks_context()
                        mode = _ks_mode_value()
                        selection = str(_get_input_value(input, _prefixed_id(prefix, "ks_choice")) or "")
                        payload["_bookmark_key"] = cfg.get("key", prefix)
                        payload["_persist_token"] = _ks_stable_persist_token(selection, mode, ctx, settings_override)
                    else:
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
            catalog_info = payload.get("_global_protein_catalog") or _current_global_catalog_info()

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
            if data.get("_viewer_kind") == "cst_pdf":
                cst_payload = dict(data.get("_cst_payload") or {})
                cst_payload["_active_fc_index"] = fc_idx or 1
                return ui.TagList(
                    ui.tags.script(
                        f"""
                        (function(){{
                            const overlay = document.getElementById('{_prefixed_id(prefix, "viewer_overlay_panel")}');
                            const editor = document.getElementById('{_prefixed_id(prefix, "viewer_create_panel")}');
                            if (overlay) {{
                                overlay.style.display = '';
                            }}
                            if (editor) editor.classList.add('is-hidden-for-cst');
                        }})();
                        """
                    ),
                    _cst_create_panel(prefix),
                    create_cst_pathway_viewer(
                        cst_payload,
                        save_input_id=_prefixed_id(prefix, "cst_save_state"),
                        export_key=prefix,
                    ),
                )
            hide_cst_script = ui.tags.script(
                f"""
                (function(){{
                    const overlay = document.getElementById('{_prefixed_id(prefix, "viewer_overlay_panel")}');
                    const editor = document.getElementById('{_prefixed_id(prefix, "viewer_create_panel")}');
                    if (overlay) {{
                        overlay.style.display = '';
                    }}
                    if (editor) editor.classList.remove('is-hidden-for-cst');
                }})();
                """
            )
            settings = data.get("general_data", {}).get("settings", {})
            show_bg = bool(settings.get("show_background_image", False))
            data_with_fc = dict(data)
            data_with_fc["_active_fc_index"] = fc_idx
            return ui.TagList(
                hide_cst_script,
                create_pathway_svg(data_with_fc, show_kegg_bg=show_bg),
            )

        @output(id=_prefixed_id(prefix, "status_message"))
        @render.text
        def status_message() -> str:
            return state["status"].get()

        @reactive.effect
        @reactive.event(getattr(input, _prefixed_id(prefix, "cst_save_state")))
        def _handle_cst_save_state() -> None:
            payload = _get_input_value(input, _prefixed_id(prefix, "cst_save_state"))
            if not payload or not isinstance(payload, dict):
                return
            if not PERSISTENT_CST_SAVE_ENABLED:
                state["status"].set("CST disk persistence is disabled for this deployment.")
                return
            pathway_id = str(payload.get("pathway_id") or "").strip()
            pathway_name = str(payload.get("pathway_name") or "").strip()
            nodes = list(payload.get("nodes") or [])
            edges = list(payload.get("edges") or [])
            groups = list(payload.get("groups") or [])
            disable_pdf_reader = bool(payload.get("disable_pdf_reader"))
            if not pathway_id:
                return

            catalog_entry = _resolve_cst_catalog_entry(pathway_id)
            if not catalog_entry:
                state["status"].set("Unable to save CST pathway edits.")
                return
            file_path_raw = str(catalog_entry.get("file_path") or "").strip()
            if not file_path_raw:
                state["status"].set("Unable to save CST pathway edits.")
                return
            file_path = Path(file_path_raw)
            if not file_path.exists():
                state["status"].set("Unable to save CST pathway edits.")
                return
            resolved_pathway_name = str(pathway_name or catalog_entry.get("name") or file_path.stem).strip()

            session_state_path: Optional[Path] = None
            try:
                safe_pathway_id = re.sub(r"[^a-z0-9._-]+", "_", str(pathway_id or "").strip().lower())
                safe_pathway_id = safe_pathway_id.strip("._-") or "cst_pathway"
                session_state_path = safe_session_path(
                    session,
                    f"cst_overlay_state/{safe_pathway_id}.mapkinase.json",
                )
                save_cst_overlay_state(
                    file_path,
                    pathway_id=pathway_id,
                    pathway_name=resolved_pathway_name,
                    nodes=nodes,
                    edges=edges,
                    groups=groups,
                    disable_pdf_reader=disable_pdf_reader,
                    state_file_path=session_state_path,
                )
                _set_cst_session_state_path(pathway_id, session_state_path)
            except Exception as exc:
                print(f"Warning: failed to save CST pathway edits for {pathway_id}: {exc}")
                state["status"].set("Unable to save CST pathway edits.")
                return

            with reactive.isolate():
                prot_data = protein_dataset.get()
                metabolite_data = metabolite_dataset.get()
                ptm_data = ptm_dataset.get()
            settings_override = collect_settings(input, cfg)
            refreshed = load_cst_pathway_payload(
                pathway_id,
                Path(BASE_DIR),
                protein_dataset=prot_data,
                metabolite_dataset=metabolite_data,
                ptm_dataset=ptm_data,
                ptm_settings={
                    "ptm_max_display": settings_override.get("ptm_max_display", DEFAULT_SETTINGS["ptm_max_display"]),
                    "ptm_label_font": settings_override.get("ptm_label_font", DEFAULT_SETTINGS["ptm_label_font"]),
                    "ptm_label_color": settings_override.get("ptm_label_color", DEFAULT_SETTINGS["ptm_label_color"]),
                    "ptm_label_size": settings_override.get("ptm_label_size", DEFAULT_SETTINGS["ptm_label_size"]),
                    "ptm_circle_radius": settings_override.get("ptm_circle_radius", DEFAULT_SETTINGS["ptm_circle_radius"]),
                    "ptm_circle_spacing": settings_override.get("ptm_circle_spacing", DEFAULT_SETTINGS["ptm_circle_spacing"]),
                    "ptm_outline_width": settings_override.get("ptm_outline_width", DEFAULT_SETTINGS["ptm_outline_width"]),
                    "ptm_shape": (_get_input_value(input, "input_ptm_shape") or "circle"),
                },
                negative_color=settings_override.get("negative_color", DEFAULT_SETTINGS["negative_color"]),
                positive_color=settings_override.get("positive_color", DEFAULT_SETTINGS["positive_color"]),
                max_negative=settings_override.get("max_negative", DEFAULT_SETTINGS["max_negative"]),
                max_positive=settings_override.get("max_positive", DEFAULT_SETTINGS["max_positive"]),
                gradient_stops=settings_override.get("gradient_stops"),
                prot_outline_width=settings_override.get("prot_outline_width", DEFAULT_SETTINGS["prot_outline_width"]),
                use_black_protein_outlines=bool(settings_override.get("use_black_protein_outlines", DEFAULT_SETTINGS.get("use_black_protein_outlines", False))),
                simple_kegg_mode=bool(settings_override.get("simple_kegg_mode", True)),
                temporal_mode=bool(settings_override.get("temporal_mode", DEFAULT_SETTINGS.get("temporal_mode", False))),
                overlay_state_path=session_state_path,
            )
            current_json = dict(state["json"].get() or {})
            if (
                refreshed
                and current_json.get("_viewer_kind") == "cst_pdf"
                and not bool(settings_override.get("simple_kegg_mode", True))
            ):
                current_json["_cst_payload"] = refreshed
                current_json["_persist_token"] = time.time()
                state["json"].set(current_json)
            state["status"].set(f"Saved CST pathway edits for {resolved_pathway_name}.")

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
            try:
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
                yield json.dumps(
                    {"error": "custom pathway export failed", "schema_version": CUSTOM_LAYOUT_SCHEMA_VERSION},
                    ensure_ascii=True,
                ).encode("utf-8")

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
                state["status"].set("Failed to export custom pathway.")
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
            if prefix == "web":
                _set_load_feedback("")
            state["loading"].set(True)
            try:
                state["status"].set("Generating pathway...")
                build_json()
            finally:
                state["loading"].set(False)

        @reactive.Effect
        @reactive.event(getattr(input, _prefixed_id(prefix, "upload_custom_pathway")))
        def _import_custom_pathway():
            upload = _get_input_value(input, _prefixed_id(prefix, "upload_custom_pathway"))
            if not upload:
                return
            file_info = upload[0]
            upload_temp_path = _extract_upload_datapath(file_info)
            try:
                upload_validation = validate_uploaded_file(file_info, expected_role="custom_pathway", session_id=_safe_session_id())
            except Exception:
                _cleanup_upload_temp_file(upload_temp_path, "custom pathway import validation exception")
                state["status"].set("Invalid custom pathway import: file could not be validated.")
                return
            if not upload_validation.valid:
                _cleanup_upload_temp_file(upload_validation.datapath, "invalid custom pathway import validation")
                state["status"].set(upload_validation.user_message)
                return
            datapath = upload_validation.datapath
            try:
                with open(datapath, "r", encoding="utf-8") as fh:
                    raw_data = json.load(fh)
                # Custom pathway imports are parsed strictly as JSON data (no code execution).
                # The sanitized layout stays in reactive memory and is not written to public/static paths.
                layout = _sanitize_custom_layout(raw_data)
            except (OSError, json.JSONDecodeError, ValueError):
                state["status"].set("Invalid custom pathway import: file must be a valid pathway JSON object.")
                return
            finally:
                _cleanup_upload_temp_file(datapath, "custom pathway import parsed")
            state["custom_layout"].set(layout)
            # Apply immediately
            if prefix == "web":
                try:
                    session.send_input_message(_prefixed_id(prefix, "pathway_id"), {"value": "Custom"})
                except Exception:
                    pass
            build_json()
            with reactive.isolate():
                built_payload = state["json"].get()
            if isinstance(built_payload, dict) and built_payload.get("_custom_layout_applied"):
                match_stats = built_payload.get("_custom_layout_dataset_match_stats") or {}
                if isinstance(match_stats, dict) and "matched_boxes" in match_stats:
                    state["status"].set(
                        "Custom pathway imported and applied. "
                        f"Matched {int(match_stats.get('matched_boxes') or 0)} of "
                        f"{int(match_stats.get('boxes_with_ids') or 0)} protein box(es) to the current dataset."
                    )
                else:
                    state["status"].set("Custom pathway imported and applied.")

        if cfg.get("start_blank"):
            build_json()

    @output(id="download_sample_protein_dataset")
    @render.download(filename=os.path.basename(SAMPLE_PROTEIN_FILE))
    def download_sample_protein_dataset():
        try:
            with open(SAMPLE_PROTEIN_FILE, "rb") as sample_fh:
                yield sample_fh.read()
        except Exception as exc:
            print(f"download_sample_protein_dataset error: {exc}")
            yield b""

    @output(id="download_sample_ptm_dataset")
    @render.download(filename=os.path.basename(SAMPLE_PTM_FILE))
    def download_sample_ptm_dataset():
        try:
            with open(SAMPLE_PTM_FILE, "rb") as sample_fh:
                yield sample_fh.read()
        except Exception as exc:
            print(f"download_sample_ptm_dataset error: {exc}")
            yield b""

    @output(id="download_sample_metabolite_dataset")
    @render.download(filename=os.path.basename(SAMPLE_METABOLITE_FILE))
    def download_sample_metabolite_dataset():
        try:
            with open(SAMPLE_METABOLITE_FILE, "rb") as sample_fh:
                yield sample_fh.read()
        except Exception as exc:
            print(f"download_sample_metabolite_dataset error: {exc}")
            yield b""

    for cfg in BOOKMARK_CONFIGS:
        _register_bookmark(cfg)


app = App(
    app_ui,
    server,
    debug=(not IS_PRODUCTION_MODE and env_truthy("MAPKINASE_APP_DEBUG", False)),
)

# --- Health check wrapper ---
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import JSONResponse
from starlette.routing import Mount, Route

async def health(request: Request):
    return JSONResponse({"status": "ok"})

_shiny_app = app
app = Starlette(
    routes=[
        Route("/health", health),
        Route("/health/", health),
        Mount("/", app=_shiny_app),
    ]
)

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
    if GUI_POPUP:
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
