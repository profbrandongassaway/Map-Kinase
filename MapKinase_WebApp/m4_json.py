import sys
import base64
import pandas as pd
import numpy as np
import os
import re
import math
from datetime import datetime
from collections import defaultdict
import json
import html
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from MapKinase_WebApp.a1_factory import get_pathway_api


def _safe_debug_print(*parts):
    text = " ".join(str(part) for part in parts)
    try:
        print(text)
    except UnicodeEncodeError:
        encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
        safe_text = text.encode(encoding, errors="replace").decode(encoding, errors="replace")
        print(safe_text)


_SYMBOL_ICON_CACHE = {}
_DEBUG_TRUE_VALUES = {"1", "true", "yes", "on"}


def _debug_file_output_enabled() -> bool:
    return str(os.getenv("MAPKINASE_DEBUG_FILE_OUTPUT", "false") or "").strip().lower() in _DEBUG_TRUE_VALUES


def _resolve_debug_output_dir(settings: dict | None = None) -> Path | None:
    cfg = settings or {}
    candidates = [
        str(cfg.get("session_workspace_dir") or "").strip(),
        str(cfg.get("debug_output_dir") or "").strip(),
        str(os.getenv("MAPKINASE_DEBUG_OUTPUT_DIR", "") or "").strip(),
    ]
    for raw in candidates:
        if not raw:
            continue
        path = Path(raw).expanduser()
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        return path
    return None


def _resolve_icons_dir() -> Path:
    candidates = [Path(__file__).resolve().parent / "icons"]
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        candidates.extend([
            Path(sys._MEIPASS) / "MapKinase_WebApp" / "icons",
            Path(sys._MEIPASS) / "icons",
        ])
    for path in candidates:
        if path.is_dir():
            return path
    return candidates[0]


_ICONS_DIR = _resolve_icons_dir()


def normalize_uniprot(value):
    token = str(value or "").strip().upper()
    if not token:
        return ""
    token = token.split(";")[0].strip()
    token = token.split(",")[0].strip()
    if "-" in token:
        base, suffix = token.split("-", 1)
        if suffix.isdigit():
            token = base
    return token


def _load_symbol_icon_data_uri(icon_name):
    if not icon_name:
        return None
    cache_key = str(icon_name).strip()
    if not cache_key:
        return None
    if cache_key in _SYMBOL_ICON_CACHE:
        return _SYMBOL_ICON_CACHE[cache_key]
    icon_path = _ICONS_DIR / cache_key
    if not icon_path.exists():
        _SYMBOL_ICON_CACHE[cache_key] = None
        return None
    try:
        svg_bytes = icon_path.read_bytes()
        encoded = base64.b64encode(svg_bytes).decode("ascii")
        data_uri = f"data:image/svg+xml;base64,{encoded}"
        _SYMBOL_ICON_CACHE[cache_key] = data_uri
        return data_uri
    except Exception:
        _SYMBOL_ICON_CACHE[cache_key] = None
        return None


def _resolve_symbol_icon(symbol_dict):
    if not isinstance(symbol_dict, dict):
        return None
    explicit = str(symbol_dict.get("symbol_icon") or "").strip()
    if explicit:
        return _load_symbol_icon_data_uri(explicit)
    symbol_text = str(symbol_dict.get("symbol") or "").strip().lower()
    icon_name = {
        "↑": "active.svg",
        "â†‘": "active.svg",
        "x": "inhibit.svg",
        "+": "regsite.svg",
    }.get(symbol_text)
    return _load_symbol_icon_data_uri(icon_name) if icon_name else None


def calculate_angle(x1, y1, x2, y2):
    return np.arctan2(y2 - y1, x2 - x1)

def determine_arrow_side(angle):
    ANGLE_THRESHOLDS = {
        'right': (-np.pi/4, np.pi/4),
        'bottom': (np.pi/4, 3*np.pi/4),
        'top': (-3*np.pi/4, -np.pi/4),
        'left': (3*np.pi/4, np.pi)
    }
    for side, (lower, upper) in ANGLE_THRESHOLDS.items():
        if lower <= angle < upper:
            return side
    return 'left'

def map_side(side):
    if side == 'right':
        return 'East'
    elif side == 'left':
        return 'West'
    elif side == 'top':
        return 'North'
    elif side == 'bottom':
        return 'South'
    return 'East'

def get_possible_outgoing_sides(dx, dy):
    if dx > 0:
        if dy > 0:
            return ['East', 'South']
        elif dy < 0:
            return ['East', 'North']
        else:
            return ['East']
    elif dx < 0:
        if dy > 0:
            return ['West', 'South']
        elif dy < 0:
            return ['West', 'North']
        else:
            return ['West']
    else:
        if dy > 0:
            return ['South']
        elif dy < 0:
            return ['North']
        else:
            return []

def adjust_arrow_endpoint(x, y, width, height, side):
    if side == 'right':
        return x - 10, y + height / 2
    elif side == 'left':
        return x + 10, y + height / 2
    elif side == 'top':
        return x + width / 2, y + 16
    elif side == 'bottom':
        return x + width / 2, y - 16
    return x + width / 2, y + height / 2

def _normalize_fc_suffix(header: str) -> str:
    value = str(header or "").strip()
    if ":" in value:
        value = value.split(":", 1)[1]
    value = re.sub(r"\s+", " ", value.strip())
    return value.lower()

def _outline_column_lookup(headers):
    outline_map = {}
    header_list = headers if headers is not None else []
    for col in header_list:
        if not isinstance(col, str):
            continue
        if col.strip().lower().startswith("o:"):
            key = _normalize_fc_suffix(col)
            if key:
                outline_map[key] = col
    return outline_map

def _resolve_outline_columns(fold_change_columns, headers):
    outline_map = _outline_column_lookup(headers)
    resolved = []
    fc_list = fold_change_columns if fold_change_columns is not None else []
    for col in fc_list:
        key = _normalize_fc_suffix(col)
        resolved.append(outline_map.get(key))
    return resolved


HMDB_METABOLITES_PARSED_FILE = Path(__file__).resolve().parent.parent / "stored_pathways" / "hmdb_metabolites_parsed.csv"


def _normalize_metabolite_token(value):
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    if not text:
        return ""
    text = text.replace("cpd:", "").replace("CPD:", "")
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text)
    return text.casefold()


def _display_metabolite_name(value):
    raw = str(value or "").strip()
    if not raw:
        return ""
    return re.sub(r"\s+", " ", raw.replace("_", " ")).strip()


def _split_metabolite_values(value):
    if value is None:
        return []
    try:
        if pd.isna(value):
            return []
    except Exception:
        pass
    text = str(value).strip()
    if not text:
        return []
    if "|" in text:
        return [part.strip() for part in text.split("|") if part.strip()]
    if ";" in text:
        return [part.strip() for part in text.split(";") if part.strip()]
    return [text]


@lru_cache(maxsize=1)
def _load_hmdb_metabolite_reference():
    if not HMDB_METABOLITES_PARSED_FILE.exists():
        return pd.DataFrame(), {}
    try:
        df = pd.read_csv(HMDB_METABOLITES_PARSED_FILE, dtype=str).fillna("")
    except Exception as exc:
        print(f"Warning: could not load HMDB metabolite reference {HMDB_METABOLITES_PARSED_FILE}: {exc}")
        return pd.DataFrame(), {}
    hmdb_map = {}
    for _, row in df.iterrows():
        accession = str(row.get("accession", "") or "").strip()
        for candidate in [accession, * _split_metabolite_values(row.get("secondary_accessions", ""))]:
            key = _normalize_metabolite_token(candidate)
            if key and key not in hmdb_map:
                hmdb_map[key] = row.to_dict()
    return df, hmdb_map


@lru_cache(maxsize=4)
def _build_cached_reference_metabolite_index(source_name):
    source = str(source_name or 'kegg').lower()
    ref_df, _ = _load_hmdb_metabolite_reference()
    if ref_df is None or ref_df.empty:
        return {}
    index = {}
    for _, row in ref_df.iterrows():
        ref = row.to_dict()
        accession = str(ref.get('accession', '') or '').strip()
        display_label = _display_metabolite_name(ref.get('wikipedia_id') or ref.get('name') or accession)
        native_keys = set()
        if source == 'kegg':
            for raw_id in _split_metabolite_values(ref.get('kegg_id', '')):
                key = _normalize_metabolite_token(raw_id)
                if key:
                    native_keys.add(key)
        if source == 'wikipathways':
            for raw_name in _split_metabolite_values(ref.get('wikipedia_id', '')):
                key = _normalize_metabolite_token(raw_name)
                if key:
                    native_keys.add(key)
        for raw_name in [ref.get('name', ''), ref.get('wikipedia_id', ''), *_split_metabolite_values(ref.get('synonyms', ''))]:
            key = _normalize_metabolite_token(raw_name)
            if key:
                native_keys.add(key)
        for raw_id in [
            ref.get('chebi_id', ''),
            ref.get('pubchem_compound_id', ''),
            ref.get('drugbank_id', ''),
            ref.get('foodb_id', ''),
        ]:
            key = _normalize_metabolite_token(raw_id)
            if key:
                native_keys.add(key)
        if not native_keys:
            continue
        match_payload = {
            'hmdb_id': accession,
            'display_label': display_label,
            'pathway_ids': sorted(native_keys),
            'reference': ref,
            'row_map': {},
        }
        for native_key in native_keys:
            existing = index.get(native_key)
            if not existing or (not existing.get('display_label') and display_label):
                index[native_key] = match_payload
    return index

def parse_exceptions_file(file_path='exceptions_file.txt'):
    exceptions = {
        'global': {'proximity': [], 'specific': {}},
        'hsa_specific': {}
    }
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        candidate_dirs = [
            os.path.join(sys._MEIPASS, 'MapKinase_WebApp'),
            os.path.join(sys._MEIPASS, 'Scripts'),
            sys._MEIPASS,
        ]
    else:
        candidate_dirs = [os.path.dirname(os.path.abspath(__file__))]
    full_path = None
    for base_dir in candidate_dirs:
        candidate = os.path.join(base_dir, file_path)
        if os.path.exists(candidate):
            full_path = candidate
            break
    try:
        if not full_path:
            print(f"Exceptions file not found (searched: {candidate_dirs})")
            return exceptions
        with open(full_path, 'r') as f:
            current_hsa = None
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('HSA:'):
                    current_hsa = line.split('HSA:')[1].strip()
                    exceptions['hsa_specific'][current_hsa] = {'proximity': [], 'specific': {}}
                    continue
                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 2:
                    continue
                rule_type = parts[0].lower()
                if rule_type == 'proximity':
                    conditions = parts[1].split('&&')
                    conditions = [cond.strip() for cond in conditions]
                    action = parts[2].lower()
                    num_conditions = len(conditions)
                    threshold_start_idx = 3
                    threshold_end_idx = threshold_start_idx + num_conditions
                    if len(parts) < threshold_end_idx:
                        continue
                    thresholds = parts[threshold_start_idx:threshold_end_idx]
                    try:
                        thresholds = [float(t) for t in thresholds]
                    except ValueError:
                        continue
                    remaining_values = parts[threshold_end_idx:]
                    if action == 'move_circle' and len(remaining_values) != 3:
                        continue
                    elif action == 'move_label' and len(remaining_values) != 4:
                        continue
                    rule = {
                        'conditions': list(zip(conditions, thresholds)),
                        'action': action,
                        'values': remaining_values,
                        'priority': num_conditions
                    }
                    if current_hsa:
                        exceptions['hsa_specific'][current_hsa]['proximity'].append(rule)
                    else:
                        exceptions['global']['proximity'].append(rule)
                elif rule_type == 'specific':
                    prot_id = parts[1]
                    rule = parts[2:]
                    if current_hsa:
                        exceptions['hsa_specific'][current_hsa]['specific'].setdefault(prot_id, []).append(rule)
                    else:
                        exceptions['global']['specific'].setdefault(prot_id, []).append(rule)
        return exceptions
    except (IOError, PermissionError) as e:
        print(f"Error reading exceptions file {full_path}: {str(e)}")
        return exceptions

def calculate_proximities(genes, proteomic_data):
    protein_coords = {}
    proximities = {}
    protein_to_entry_ids = defaultdict(list)
    for gene in genes:
        entry_id = gene["id"]
        proteins = gene["name"].split()
        protein_coords[entry_id] = {
            'x': gene["x"],
            'y': gene["y"],
            'width': gene["width"],
            'height': gene["height"],
            'proteins': proteins
        }
        for protein in proteins:
            protein_to_entry_ids[protein].append(entry_id)
    for entry_id1, coords1 in protein_coords.items():
        proximities[entry_id1] = {}
        for entry_id2, coords2 in protein_coords.items():
            if entry_id1 == entry_id2:
                continue
            left1, right1 = coords1['x'], coords1['x'] + coords1['width']
            top1, bottom1 = coords1['y'], coords1['y'] + coords1['height']
            left2, right2 = coords2['x'], coords2['x'] + coords2['width']
            top2, bottom2 = coords2['y'], coords2['y'] + coords2['height']
            dx = min(abs(left1 - right2), abs(right1 - left2)) if left1 < right2 < right1 or left2 < right1 < left2 else max(0, min(abs(left1 - left2), abs(right1 - right2)))
            dy = min(abs(top1 - bottom2), abs(bottom1 - top2)) if top1 < bottom2 < bottom1 or top2 < bottom1 < top2 else max(0, min(abs(top1 - top2), abs(bottom1 - bottom2)))
            dx_signed = coords2['x'] - coords1['x']
            dy_signed = coords2['y'] - coords1['y']
            dist = np.sqrt(dx ** 2 + dy ** 2)
            proximities[entry_id1][entry_id2] = {'dx': dx_signed, 'dy': dy_signed, 'dist': dist, 'edge_dx': dx, 'edge_dy': dy}
    return proximities, protein_coords, protein_to_entry_ids

def savefile(input_file, directory, suffix, extension=".json", include_timestamp=True):
    try:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        script_dir = os.path.dirname(os.path.abspath(__file__))
        phosmap_dir = os.path.dirname(script_dir)
        output_dir = directory if os.path.isabs(str(directory or "")) else os.path.join(phosmap_dir, directory)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if include_timestamp:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            file_name = f"{base_name}{suffix}_{timestamp}{extension}"
        else:
            file_name = f"{base_name}{suffix}{extension}"
        return os.path.join(output_dir, file_name)
    except (OSError, PermissionError) as e:
        print(f"Error creating output file path for {input_file}: {e}")
        return None

def find_entry_or_group(entry_id, entries, groups):
    for entry in entries:
        if entry["id"] == entry_id:
            return entry
    for group in groups:
        if group["id"] == entry_id:
            group["type"] = "group"
            return group
    return None

def classify_phosphosite_function(phospho_data, ptm_symbol_list, reg_site_col='C: Regulatory site',
                                  reg_function_col='C: Regulatory site function',
                                  output_col='Phosphosite_Classification'):
    phospho_data[output_col] = 'none'
    _safe_debug_print(f"Classifying {len(phospho_data)} PTM sites with symbol list: {ptm_symbol_list}")
    _safe_debug_print(f"PTM data columns: {phospho_data.columns.tolist()}")

    valid_symbol_list = []
    for symbol_dict in ptm_symbol_list:
        if not symbol_dict or len(symbol_dict) != 1:
            _safe_debug_print(f"Warning: Skipping invalid rule with incorrect structure: {symbol_dict}")
            continue
        inner_dict_key = list(symbol_dict.keys())[0]
        inner_dict = symbol_dict[inner_dict_key]
        header = inner_dict.get('header_to_search', '')
        if not header:
            _safe_debug_print(f"Warning: Skipping rule with empty header_to_search: {symbol_dict}")
            continue
        if header not in phospho_data.columns:
            _safe_debug_print(f"Warning: Header {header} not in PTM data columns: {phospho_data.columns}")
            continue
        valid_symbol_list.append(symbol_dict)

    if not valid_symbol_list:
        _safe_debug_print("Error: No valid rules in ptm_symbol_list. All PTMs will be classified as 'none'.")
        return phospho_data

    for idx in phospho_data.index:
        ptm_label_type = 'none'
        for symbol_dict in [d for d in valid_symbol_list if
                            d[list(d.keys())[0]].get('statement_type') != 'if_and_notlabeled_statement']:
            inner_dict = symbol_dict[list(symbol_dict.keys())[0]]
            statement_type = inner_dict.get('statement_type', 'NA')
            header = inner_dict.get('header_to_search', '')
            search_text_1 = inner_dict.get('search_text_1', '')
            search_text_2 = inner_dict.get('search_text_2', '')
            search_text_3 = inner_dict.get('search_text_3', '')
            search_text_4 = inner_dict.get('search_text_4', '')
            label_key = list(symbol_dict.keys())[0]
            symbol_label = label_key.replace('_dict', '')
            function_text = str(phospho_data.loc[idx, header]).lower()
            assign_label = False
            if statement_type == 'NA':
                continue
            elif statement_type == 'if_statment' and search_text_1:
                assign_label = search_text_1.lower() in function_text
            elif statement_type == 'if_not_statement' and search_text_1:
                assign_label = search_text_1.lower() not in function_text
            elif statement_type == 'if_and_statement' and search_text_1 and search_text_2:
                assign_label = (search_text_1.lower() in function_text and
                                search_text_2.lower() in function_text)
            elif statement_type == 'if_or_statement' and search_text_1 and search_text_2:
                assign_label = (search_text_1.lower() in function_text or
                                search_text_2.lower() in function_text)
            elif statement_type == 'if_and_not_statement' and search_text_1 and search_text_2:
                assign_label = (search_text_1.lower() in function_text and
                                search_text_2.lower() not in function_text)
            elif statement_type == 'if_and_and_not_statement' and search_text_1 and search_text_2 and search_text_3:
                assign_label = (search_text_1.lower() in function_text and
                                search_text_2.lower() in function_text and
                                search_text_3.lower() not in function_text)
            elif statement_type == 'if_or_and_not_statement' and search_text_1 and search_text_2 and search_text_3:
                assign_label = ((search_text_1.lower() in function_text or
                                 search_text_2.lower() in function_text) and
                                search_text_3.lower() not in function_text)
            elif statement_type == 'if_or_and_or_statement' and search_text_1 and search_text_2 and search_text_3 and search_text_4:
                assign_label = ((search_text_1.lower() in function_text or
                                 search_text_2.lower() in function_text) and
                                (search_text_3.lower() in function_text or
                                 search_text_4.lower() in function_text))
            if assign_label:
                ptm_label_type = symbol_label
        for symbol_dict in [d for d in valid_symbol_list if
                            d[list(d.keys())[0]].get('statement_type') == 'if_and_notlabeled_statement']:
            inner_dict = symbol_dict[list(symbol_dict.keys())[0]]
            statement_type = inner_dict.get('statement_type', 'NA')
            header = inner_dict.get('header_to_search', '')
            search_text_1 = inner_dict.get('search_text_1', '')
            label_key = list(symbol_dict.keys())[0]
            symbol_label = label_key.replace('_dict', '')
            function_text = str(phospho_data.loc[idx, header]).lower()
            if statement_type == 'if_and_notlabeled_statement' and search_text_1:
                if ptm_label_type == 'none' and search_text_1.lower() in function_text:
                    ptm_label_type = symbol_label
        phospho_data.loc[idx, output_col] = ptm_label_type
    _safe_debug_print(f"Classification complete. Value counts:\n{phospho_data[output_col].value_counts()}")
    return phospho_data

class PathwayProcessor:
    def __init__(self, entries, proteomic_data, ptm_datasets, settings, metabolite_data=None, metabolite_config=None):
        self.proteomic_data = proteomic_data
        self.ptm_datasets = {f"ptm_{i}": dataset for i, dataset in enumerate(ptm_datasets)}
        self.protein_data_map = {}  # Initialize protein_data_map
        self.current_proteins = {}  # Initialize current_proteins
        for dataset_id, dataset in self.ptm_datasets.items():
            if 'dataframe' in dataset and isinstance(dataset['dataframe'], pd.DataFrame):
                dataset['data'] = dataset['dataframe'].copy()
            else:
                dataset['data'] = pd.read_csv(dataset['file_path'], sep="\t")
            _safe_debug_print(f"PTM dataset {dataset_id}: {len(dataset['data'])} rows, columns: {dataset['data'].columns.tolist()}")
            main_cols = []
            for col in dataset.get('main_columns', []):
                if isinstance(col, (list, tuple)) and len(col) > 1:
                    main_cols.append(str(col[1]))
                else:
                    main_cols.append(str(col))
            dataset['outline_columns'] = _resolve_outline_columns(main_cols, dataset['data'].columns)
            symbol_list = dataset.get('ptm_symbol_list', [])
            reg_site_col = dataset.get('modulation_column', 'C: Regulatory site')
            reg_function_col = dataset.get('tooltip_columns', ['C: Regulatory site function'])[0]
            _safe_debug_print("Calling classify_phosphosite_function with:")
            _safe_debug_print(f"  ptm_symbol_list: {symbol_list}")
            _safe_debug_print(f"  reg_site_col: {reg_site_col}")
            _safe_debug_print(f"  reg_function_col: {reg_function_col}")
            if not symbol_list:
                _safe_debug_print(f"Warning: ptm_symbol_list is empty for dataset {dataset_id}. PTMs will be classified as 'none'.")
            if reg_site_col not in dataset['data'].columns:
                _safe_debug_print(f"Warning: reg_site_col {reg_site_col} not in PTM data columns.")
            if reg_function_col not in dataset['data'].columns:
                _safe_debug_print(f"Warning: reg_function_col {reg_function_col} not in PTM data columns.")
            dataset['data'] = classify_phosphosite_function(
                dataset['data'],
                ptm_symbol_list=symbol_list,
                reg_site_col=reg_site_col,
                reg_function_col=reg_function_col,
                output_col='Phosphosite_Classification'
            )
        self.settings = settings
        self.protein_selection_option = settings.get('protein_selection_option', 2)
        self.ptm_selection_option = settings.get('ptm_selection_option', 2)
        self.fold_change_columns = settings.get('main_columns', ['ER+_Est-_x-y_TNBC'])
        self.outline_columns = _resolve_outline_columns(self.fold_change_columns, self.proteomic_data.columns)
        self.hsa_id_column = settings.get('hsa_id_column', 'KEGG_hsa')
        self.prot_uniprot_column = settings.get('prot_uniprot_column', 'Uniprot_ID')
        self.gene_name_column = settings.get('gene_name_column', 'Gene Symbol')
        self.protein_tooltip_columns = settings.get('protein_tooltip_columns', ['Gene Symbol', 'Uniprot_ID'])
        self.negative_color = settings.get('negative_color', (255, 0, 0))
        self.positive_color = settings.get('positive_color', (0, 0, 255))
        self.max_negative = settings.get('max_negative', -2)
        self.max_positive = settings.get('max_positive', 2)
        self.gradient_stops = self._normalize_gradient_stops(settings.get('gradient_stops'))
        if len(self.gradient_stops) < 2:
            self.gradient_stops = self._normalize_gradient_stops(
                [
                    {'value': float(self.max_negative), 'color': list(self.negative_color)},
                    {'value': 0.0, 'color': [255, 255, 255]},
                    {'value': float(self.max_positive), 'color': list(self.positive_color)},
                ]
            )
        if len(self.gradient_stops) >= 2:
            self.max_negative = float(self.gradient_stops[0]['value'])
            self.max_positive = float(self.gradient_stops[-1]['value'])
            self.negative_color = tuple(int(v) for v in self.gradient_stops[0]['color'])
            self.positive_color = tuple(int(v) for v in self.gradient_stops[-1]['color'])
        self.ptm_label_color = settings.get('ptm_label_color', (0, 0, 0))
        self.ptm_circle_radius = settings.get('ptm_circle_radius', 5)
        self.ptm_circle_spacing = settings.get('ptm_circle_spacing', 4)
        self.display_types = settings.get('display_types', ['prot_box'])
        default_workers = os.cpu_count() or 4
        self._ptm_worker_count = max(1, min(settings.get('ptm_worker_count', default_workers), 16))
        self._ptm_summary_cache = {}
        self.metabolite_data = metabolite_data.copy() if isinstance(metabolite_data, pd.DataFrame) else pd.DataFrame()
        self.metabolite_config = dict(metabolite_config or {})
        self.metabolite_main_columns = list(self.metabolite_config.get('main_columns') or [])
        self.metabolite_tooltip_columns = list(self.metabolite_config.get('tooltip_columns') or [])
        self.metabolite_hmdb_column = self.metabolite_config.get('hmdb_column') or (self.metabolite_data.columns[0] if not self.metabolite_data.empty else 'HMDB_ID')
        self.pathway_source = str(self.settings.get('pathway_source', 'kegg') or 'kegg').lower()
        self.metabolite_reference_df, self.hmdb_reference_map = _load_hmdb_metabolite_reference()
        self.metabolite_index = self._build_metabolite_index()
        self.reference_metabolite_index = _build_cached_reference_metabolite_index(self.pathway_source)
        self._protein_id_index = self._build_protein_match_index(self.hsa_id_column)
        self._protein_gene_index = self._build_protein_match_index(self.gene_name_column)
        self._protein_uniprot_index = self._build_protein_match_index(self.prot_uniprot_column, normalize_uniprot_values=True)

    def _normalize_gradient_stops(self, gradient_stops):
        if not isinstance(gradient_stops, (list, tuple)):
            return []
        parsed = []
        for item in gradient_stops:
            if not isinstance(item, dict):
                continue
            raw_value = item.get('value')
            try:
                value = float(raw_value)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(value):
                continue
            raw_color = item.get('color', item.get('rgb'))
            if isinstance(raw_color, str):
                color_text = raw_color.strip().lstrip('#')
                if len(color_text) != 6:
                    continue
                try:
                    color = [int(color_text[0:2], 16), int(color_text[2:4], 16), int(color_text[4:6], 16)]
                except Exception:
                    continue
            elif isinstance(raw_color, (list, tuple)) and len(raw_color) >= 3:
                try:
                    color = [int(max(0, min(255, round(float(c))))) for c in list(raw_color)[:3]]
                except Exception:
                    continue
            else:
                continue
            parsed.append((value, color))
        if not parsed:
            return []
        dedup = {}
        for value, color in parsed:
            dedup[float(value)] = list(color)
        return [{'value': float(v), 'color': list(c)} for v, c in sorted(dedup.items(), key=lambda kv: kv[0])]

    def _color_from_gradient_stops(self, fold_change):
        if pd.isna(fold_change):
            return [128, 128, 128]
        if not self.gradient_stops or len(self.gradient_stops) < 2:
            return None
        try:
            fold = float(fold_change)
        except (TypeError, ValueError):
            return [128, 128, 128]
        if not math.isfinite(fold):
            return [128, 128, 128]
        if fold <= float(self.gradient_stops[0]['value']):
            return list(self.gradient_stops[0]['color'])
        if fold >= float(self.gradient_stops[-1]['value']):
            return list(self.gradient_stops[-1]['color'])
        for idx in range(len(self.gradient_stops) - 1):
            left = self.gradient_stops[idx]
            right = self.gradient_stops[idx + 1]
            left_v = float(left['value'])
            right_v = float(right['value'])
            if fold < left_v or fold > right_v:
                continue
            if right_v == left_v:
                return list(right['color'])
            t = (fold - left_v) / (right_v - left_v)
            return [
                int(max(0, min(255, round((1 - t) * left['color'][c_idx] + t * right['color'][c_idx]))))
                for c_idx in range(3)
            ]
        return list(self.gradient_stops[-1]['color'])

    def _build_protein_match_index(self, column_name, normalize_uniprot_values=False):
        index = {}
        if (
            self.proteomic_data is None
            or self.proteomic_data.empty
            or not column_name
            or column_name not in self.proteomic_data.columns
        ):
            return index
        for row_idx, raw_value in self.proteomic_data[column_name].items():
            if raw_value is None or pd.isna(raw_value):
                continue
            raw_text = str(raw_value).strip()
            if not raw_text:
                continue
            if normalize_uniprot_values:
                keys = [normalize_uniprot(raw_text)]
            else:
                keys = [raw_text.upper()]
            for key in keys:
                if key:
                    index.setdefault(key, []).append(row_idx)
        return index

    def _lookup_protein_row_indexes(self, protein_token):
        token = str(protein_token or '').strip()
        if not token:
            return []
        upper_token = token.upper()
        indexes = list(self._protein_id_index.get(upper_token, []))
        if indexes:
            return indexes
        indexes = list(self._protein_gene_index.get(upper_token, []))
        if indexes:
            return indexes
        uni_token = normalize_uniprot(token)
        if uni_token:
            indexes = list(self._protein_uniprot_index.get(uni_token, []))
            if indexes:
                return indexes
        return []

    def _pick_best_row_index(self, row_indexes):
        if not row_indexes:
            return None
        best_idx = None
        best_score = -1.0
        for row_idx in row_indexes:
            try:
                row = self.proteomic_data.loc[row_idx]
            except Exception:
                continue
            candidate_score = -1.0
            for fc_col in self.fold_change_columns:
                if fc_col not in self.proteomic_data.columns:
                    continue
                value = self._safe_float(row.get(fc_col))
                if value is not None and math.isfinite(value):
                    candidate_score = max(candidate_score, abs(value))
            if best_idx is None or candidate_score > best_score:
                best_idx = row_idx
                best_score = candidate_score
        return best_idx if best_idx is not None else row_indexes[0]

    def _build_protein_match_payload(self, protein_token):
        row_indexes = self._lookup_protein_row_indexes(protein_token)
        if not row_indexes:
            return None
        row_idx = self._pick_best_row_index(row_indexes)
        if row_idx is None:
            return None
        try:
            row = self.proteomic_data.loc[row_idx]
        except Exception:
            return None
        uniprot_id = row.get(self.prot_uniprot_column) if self.prot_uniprot_column in self.proteomic_data.columns else None
        annotations = []
        for col in self.protein_tooltip_columns:
            value = row.get(col, '') if col in self.proteomic_data.columns else ''
            if value is None or pd.isna(value):
                value = ''
            annotations.append(str(value))
        return {
            'protein': protein_token,
            'uniprot_id': uniprot_id,
            'annotations': annotations,
            'row_index': row_idx,
            'row': row,
        }

    def _choose_metabolite_match(self, matches):
        if not matches:
            return None
        if len(matches) == 1:
            return matches[0]
        best = None
        best_score = -1
        for match in matches:
            row_map = match.get('row_map') or {}
            score = -1
            for col in self.metabolite_main_columns:
                value = self._safe_float(row_map.get(col))
                if value is not None:
                    score = max(score, abs(value))
            if score > best_score:
                best = match
                best_score = score
        return best or matches[0]

    def _build_metabolite_index(self):
        source = self.pathway_source
        if self.metabolite_data.empty or self.metabolite_hmdb_column not in self.metabolite_data.columns:
            return {}
        index = {}
        for _, row in self.metabolite_data.iterrows():
            hmdb_raw = row.get(self.metabolite_hmdb_column, '')
            hmdb_key = _normalize_metabolite_token(hmdb_raw)
            if not hmdb_key:
                continue
            ref = self.hmdb_reference_map.get(hmdb_key)
            if not ref:
                continue
            display_label = _display_metabolite_name(ref.get('wikipedia_id') or ref.get('name') or hmdb_raw)
            native_keys = set()
            if source == 'kegg':
                for raw_id in _split_metabolite_values(ref.get('kegg_id', '')):
                    key = _normalize_metabolite_token(raw_id)
                    if key:
                        native_keys.add(key)
            elif source == 'wikipathways':
                for raw_name in _split_metabolite_values(ref.get('wikipedia_id', '')):
                    key = _normalize_metabolite_token(raw_name)
                    if key:
                        native_keys.add(key)
            if source == 'wikipathways':
                for raw_name in [ref.get('name', ''), * _split_metabolite_values(ref.get('synonyms', ''))]:
                    key = _normalize_metabolite_token(raw_name)
                    if key:
                        native_keys.add(key)
            if not native_keys:
                continue
            row_map = row.to_dict()
            row_map = {str(k): ('' if pd.isna(v) else v) for k, v in row_map.items()}
            match_payload = {
                'hmdb_id': str(ref.get('accession') or hmdb_raw or '').strip(),
                'display_label': display_label,
                'pathway_ids': sorted(native_keys),
                'reference': ref,
                'row_map': row_map,
            }
            for native_key in native_keys:
                index.setdefault(native_key, []).append(match_payload)
        return index

    def _resolve_metabolite_match(self, entry):
        if not isinstance(entry, dict):
            return None
        candidates = []
        seen = set()
        for raw_value in [
            entry.get('name', ''),
            entry.get('first_name', ''),
            entry.get('label', ''),
            ((entry.get('xref') or {}).get('ID', '') if isinstance(entry.get('xref'), dict) else ''),
        ]:
            key = _normalize_metabolite_token(raw_value)
            if key and key not in seen:
                seen.add(key)
                candidates.append(key)
        for key in candidates:
            chosen = self._choose_metabolite_match(self.metabolite_index.get(key) or [])
            if chosen:
                return dict(chosen, pathway_match_id=key)
        for key in candidates:
            chosen = self.reference_metabolite_index.get(key)
            if chosen:
                return dict(chosen, pathway_match_id=key)
        return None

    def _build_metabolite_tooltip(self, match):
        if not match:
            return "", ""
        tooltip_plain = []
        tooltip_html = []
        display_label = str(match.get('display_label') or '').strip()
        hmdb_id = str(match.get('hmdb_id') or '').strip()
        pathway_match_id = str(match.get('pathway_match_id') or '').strip()
        if display_label:
            tooltip_plain.append(f"Metabolite: {display_label}")
            tooltip_html.append(f"<strong>Metabolite:</strong> {html.escape(display_label)}")
        if hmdb_id:
            tooltip_plain.append(f"HMDB ID: {hmdb_id}")
            tooltip_html.append(f"<strong>HMDB ID:</strong> {html.escape(hmdb_id)}")
        if pathway_match_id:
            tooltip_plain.append(f"Pathway Match: {pathway_match_id}")
            tooltip_html.append(f"<strong>Pathway Match:</strong> {html.escape(pathway_match_id)}")
        row_map = match.get('row_map') or {}
        for col in self.metabolite_tooltip_columns:
            value = str(row_map.get(col, '') or '').strip()
            if not value:
                continue
            tooltip_plain.append(f"{col}: {value}")
            tooltip_html.append(f"<strong>{html.escape(col)}:</strong> {html.escape(value)}")
        return "\n".join(tooltip_plain), "<br/>".join(tooltip_html)

    def get_color(self, fold_change):
        stop_color = self._color_from_gradient_stops(fold_change)
        if isinstance(stop_color, list):
            return stop_color
        if pd.isna(fold_change):
            return [128, 128, 128]
        white = (255, 255, 255)
        r0, g0, b0 = white
        if fold_change < 0:
            r1, g1, b1 = self.negative_color
            t = min(abs(fold_change) / abs(self.max_negative), 1)
        else:
            r1, g1, b1 = self.positive_color
            t = min(fold_change / self.max_positive, 1)
        r = int((1 - t) * r0 + t * r1)
        g = int((1 - t) * g0 + t * g1)
        b = int((1 - t) * b0 + t * b1)
        return [r, g, b]

    def choose_protein(self, proteins):
        valid_proteins = []
        for protein in proteins:
            payload = self._build_protein_match_payload(protein)
            if payload is not None:
                valid_proteins.append(payload)
        if not valid_proteins:
            return None
        try:
            if self.protein_selection_option == 1:
                selected_protein = valid_proteins[0]
            elif self.protein_selection_option == 2:
                max_fold_change = -float('inf')
                selected_protein = None
                for protein in valid_proteins:
                    row = protein.get('row')
                    if row is None:
                        continue
                    for fc_col in self.fold_change_columns:
                        try:
                            fc_val = float(row.get(fc_col))
                        except (TypeError, ValueError):
                            fc_val = None
                        if fc_val is not None and abs(fc_val) > max_fold_change:
                            max_fold_change = abs(fc_val)
                            selected_protein = protein
                if not selected_protein:
                    selected_protein = valid_proteins[0]
            else:
                selected_protein = valid_proteins[0]
            return selected_protein
        except (KeyError, IndexError) as e:
            return None

    def prioritize_ptm_sites(self, uniprot_id):
        all_ptms = []
        for dataset_id, dataset in self.ptm_datasets.items():
            uniprot_col = dataset['uniprot_column']
            fold_change_cols = [col[1] for col in dataset['main_columns']]
            ptm_data = dataset['data']
            ptm_subset = ptm_data[ptm_data[uniprot_col] == uniprot_id].copy()
            if ptm_subset.empty:
                continue
            ptm_subset['dataset_id'] = dataset_id
            ptm_subset['is_phospho'] = dataset['type'] == 'Phosphorylation'
            ptm_subset['is_modulating'] = ptm_subset.get(dataset.get('modulation_column', ''), '') == '+' if dataset['type'] == 'Phosphorylation' else True
            all_ptms.append(ptm_subset)
        if not all_ptms:
            return pd.DataFrame()
        combined_ptms = pd.concat(all_ptms, ignore_index=True)
        if combined_ptms.empty:
            return combined_ptms
        if 'dataset_id' not in combined_ptms.columns:
            return pd.DataFrame()
        if combined_ptms['dataset_id'].isna().any():
            combined_ptms = combined_ptms.dropna(subset=['dataset_id'])
        if combined_ptms.empty:
            return combined_ptms
        # Normalize fold-change columns to numeric
        fc_cols = [col[1] for dataset in self.ptm_datasets.values() for col in dataset['main_columns']]
        for fc in fc_cols:
            if fc in combined_ptms.columns:
                combined_ptms[fc] = pd.to_numeric(combined_ptms[fc], errors='coerce')
        try:
            if self.ptm_selection_option == 1:
                print("Applying PTM selection option 1: Highest absolute fold change")
                combined_ptms = combined_ptms.sort_values(
                    by=fc_cols,
                    key=lambda x: pd.to_numeric(x, errors="coerce").abs(),
                    ascending=False
                )
            elif self.ptm_selection_option == 2:
                print("Applying PTM selection option 2: Prefer modulating phospho PTMs")
                def sort_key(row):
                    try:
                        dataset_id = row['dataset_id']
                        dataset = self.ptm_datasets[dataset_id]
                        fold_change = row.get(dataset['main_columns'][0][1], 0)
                        if pd.isna(fold_change):
                            fold_change = 0.0
                        try:
                            fold_change = float(fold_change)
                        except (TypeError, ValueError):
                            fold_change = 0.0
                        is_modulating_val = row.get('is_modulating')
                        is_modulating = bool(is_modulating_val is True or is_modulating_val == '+' or str(is_modulating_val).lower() == 'true')
                        is_phospho = bool(row.get('is_phospho'))
                        return (-is_phospho, -is_modulating, -abs(fold_change))
                    except (KeyError, TypeError) as e:
                        print(f"Error in sort_key for row {row}: {e}")
                        return (0, 0, 0)
                combined_ptms['sort_key'] = combined_ptms.apply(sort_key, axis=1)
                combined_ptms = combined_ptms.sort_values(by='sort_key')
                combined_ptms = combined_ptms.drop(columns=['sort_key'])
            elif self.ptm_selection_option == 3:
                print("Applying PTM selection option 3: Modulating phospho and all non-phospho PTMs")
                combined_ptms = combined_ptms[
                    (combined_ptms['is_phospho'] & combined_ptms['is_modulating']) |
                    (~combined_ptms['is_phospho'])
                ]
                combined_ptms = combined_ptms.sort_values(
                    by=fc_cols,
                    key=lambda x: pd.to_numeric(x, errors="coerce").abs(),
                    ascending=False
                )
            print(f"Found {len(combined_ptms)} PTM sites")
            return combined_ptms
        except Exception as e:
            print(f"Error prioritizing PTM sites for UniProt ID {uniprot_id}: {str(e)}")
            return pd.DataFrame()

    def process_pathway(self, entries, groups, arrows, proteomic_data, ptm_datasets, skip_disk_write=False):
        try:
            exceptions = parse_exceptions_file()
            proximities, protein_coords, protein_to_entry_ids = calculate_proximities(entries, self.proteomic_data)
            hsa_id = self.settings.get('pathway_id', 'hsa04010')
            hsa_exceptions = exceptions['hsa_specific'].get(hsa_id,
                                                            {'proximity': [], 'specific': {}})
            global_exceptions = exceptions['global']
            chosen_proteins = set()
            defaults = {
                'N1': (-5, -5, 'right'), 'N2': (0, -11, 'center'), 'N3': (5, -5, 'left'),
                'S1': (-3, 5, 'right'), 'S2': (0, 12, 'center'), 'S3': (3, 5, 'left'),
                'W1': (-3, -2, 'right'), 'W2': (-3, 2, 'right'),
                'E1': (3, -2, 'left'), 'E2': (3, 2, 'left')
            }
            json_data = {
                'general_data': {
                    'settings': {
                        'pathway_id': self.settings.get('pathway_id', 'hsa04010'),
                        'pathway_source': self.settings.get('pathway_source', 'kegg'),
                        'protein_selection_option': self.settings.get('protein_selection_option', 2),
                        'ptm_selection_option': self.settings.get('ptm_selection_option', 2),
                        'ptm_max_display': self.settings.get('ptm_max_display', 4),
                        'show_background_image': self.settings.get('show_background_image', True),
                        'display_types': self.settings.get('display_types', ['prot_box']),
                        'show_groups': self.settings.get('show_groups', False),
                        'show_multi_protein_indicator': self.settings.get('show_multi_protein_indicator', True),
                        'show_arrows': self.settings.get('show_arrows', True),
                        'show_text_boxes': self.settings.get('show_text_boxes', True),
                        'debug_mode': self.settings.get('debug_mode', False),
                        'temporal_mode': self.settings.get('temporal_mode', False),
                        'negative_color': self.settings.get('negative_color', (255, 0, 0)),
                        'positive_color': self.settings.get('positive_color', (0, 0, 255)),
                        'max_negative': self.settings.get('max_negative', -2),
                        'max_positive': self.settings.get('max_positive', 2),
                        'gradient_stops': self.settings.get('gradient_stops', self.gradient_stops),
                        'prot_label_font': self.settings.get('prot_label_font', 'Arial'),
                        'prot_label_size': self.settings.get('prot_label_size', 12),
                        'ptm_label_font': self.settings.get('ptm_label_font', 'Arial'),
                        'ptm_label_color': self.settings.get('ptm_label_color', (0, 0, 0)),
                        'ptm_label_size': self.settings.get('ptm_label_size', 10),
                        'ptm_circle_radius': self.settings.get('ptm_circle_radius', 5),
                        'ptm_circle_spacing': self.settings.get('ptm_circle_spacing', 4),
                        'prot_outline_width': self.settings.get('prot_outline_width', 1),
                        'ptm_outline_width': self.settings.get('ptm_outline_width', 1),
                        'protein_tooltip_columns': self.settings.get('protein_tooltip_columns', ['Gene Symbol', 'Uniprot_ID']),
                        'protein_file_path': self.settings.get('protein_file_path', ''),
                        'metabolite_file_path': self.settings.get('metabolite_file_path', '')
                    },
                    'data': {
                        'protein': {
                            'file_path': self.settings.get('protein_file_path', ''),
                            'uniprot_column': self.settings.get('prot_uniprot_column', 'Uniprot_ID'),
                            'kegg_column': self.settings.get('hsa_id_column', 'KEGG_hsa'),
                            'gene_column': self.settings.get('gene_name_column', 'Gene Symbol'),
                            'main_columns': self.settings.get('main_columns', ['ER+_Est-_x-y_TNBC']),
                            'tooltip_columns': self.settings.get('protein_tooltip_columns', ['Gene Symbol', 'Uniprot_ID'])
                        },
                        'ptm': [
                            {
                                'type': dataset['type'],
                                'file_path': dataset.get('file_path', ''),
                                'uniprot_column': dataset['uniprot_column'],
                                'site_column': dataset['site_column'],
                                'shape': dataset['shape'],
                                'main_columns': dataset['main_columns'],
                                'modulation_column': dataset.get('modulation_column', ''),
                                'tooltip_columns': dataset.get('tooltip_columns', []),
                                'ptm_symbol_list': dataset.get('ptm_symbol_list', [])
                            } for dataset in ptm_datasets
                        ],
                        'metabolite': {
                            'file_path': self.metabolite_config.get('file_path', ''),
                            'hmdb_column': self.metabolite_hmdb_column,
                            'main_columns': self.metabolite_main_columns,
                            'tooltip_columns': self.metabolite_tooltip_columns,
                        }
                    }
                },
                'protein_data': {},
                'protbox_data': [],
                'groups': [],
                'arrows': [],
                'compound_data' : [],
                'text_data' : []
            }
            def _find_entry_by_id(entry_id):
                for entry in entries:
                    if str(entry.get('id')) == str(entry_id):
                        return entry
                return None

            for group in groups:
                members = []
                protbox_ids = set()
                component_ids = group.get('components') or group.get('protbox_ids') or []
                for comp_id in component_ids:
                    entry = _find_entry_by_id(comp_id)
                    if not entry:
                        continue
                    entry_type = entry.get('type')
                    if entry_type == 'prot_box':
                        members.append({'type': 'prot-box', 'id': comp_id})
                        protbox_ids.add(str(comp_id))
                    elif entry_type == 'compound':
                        members.append({'type': 'compound', 'id': comp_id})
                    elif entry_type in {'map', 'label', 'text'} or entry.get('graphics_type') in {'roundrectangle', 'label'} or entry.get('shape_type'):
                        members.append({'type': 'text-box', 'id': comp_id})
                group_entry = {
                    'group_id': group['id'],
                    'members': members,
                    'protbox_ids': list(protbox_ids)
                }
                if members and str(self.settings.get('pathway_source', 'kegg')).lower() == 'kegg':
                    group_entry['show_box'] = True
                    group_entry['box_padding'] = 10
                    group_entry['box_radius'] = 8
                json_data['groups'].append(group_entry)

            binding_edges = []
            if str(self.settings.get('pathway_source', 'kegg')).lower() == 'kegg':
                filtered_arrows = []
                for arrow in arrows:
                    if arrow.get('binding') or str(arrow.get('type', '')).lower() == 'binding/association':
                        entry1_id = arrow.get('entry1')
                        entry2_id = arrow.get('entry2')
                        if entry1_id and entry2_id:
                            binding_edges.append((str(entry1_id), str(entry2_id)))
                        continue
                    filtered_arrows.append(arrow)
                arrows = filtered_arrows
            if binding_edges:
                neighbors = defaultdict(set)
                for a_id, b_id in binding_edges:
                    neighbors[a_id].add(b_id)
                    neighbors[b_id].add(a_id)
                existing_ids = {str(g.get('group_id')) for g in json_data.get('groups', []) if isinstance(g, dict)}
                for hub_id, linked in neighbors.items():
                    if len(linked) < 2:
                        continue
                    hub_entry = _find_entry_by_id(hub_id)
                    if not hub_entry or hub_entry.get('type') != 'prot_box':
                        continue
                    members = []
                    protbox_ids = set()
                    for pid in [hub_id] + sorted(linked):
                        entry = _find_entry_by_id(pid)
                        if not entry or entry.get('type') != 'prot_box':
                            continue
                        members.append({'type': 'prot-box', 'id': pid})
                        protbox_ids.add(str(pid))
                    if len(members) < 2:
                        continue
                    group_id = f"bind_assoc_{hub_id}"
                    if group_id in existing_ids:
                        continue
                    json_data['groups'].append(
                        {
                            'group_id': group_id,
                            'members': members,
                            'protbox_ids': list(protbox_ids),
                            'show_box': True,
                            'box_padding': 10,
                            'box_radius': 8,
                            'anchor_member': hub_id,
                        }
                    )
                    existing_ids.add(group_id)
            prot_entries = [e for e in entries if e.get('type') == 'prot_box']
            prot_ids = {e.get("id") for e in prot_entries}
            opposite_side = {'East': 'West', 'West': 'East', 'North': 'South', 'South': 'North'}
            is_kegg = str(self.settings.get('pathway_source', 'kegg')).lower() == 'kegg'
            def _entry_center(entry):
                if not entry:
                    return None, None
                x_val = entry.get("x", 0)
                y_val = entry.get("y", 0)
                if is_kegg:
                    return x_val, y_val
                return x_val + entry.get("width", 0) / 2, y_val + entry.get("height", 0) / 2

            def _shorten_end(start_x, start_y, end_x, end_y, amount):
                dx = end_x - start_x
                dy = end_y - start_y
                dist = (dx ** 2 + dy ** 2) ** 0.5
                if dist <= 0:
                    return end_x, end_y
                return end_x - dx / dist * amount, end_y - dy / dist * amount

            def _shorten_start(start_x, start_y, end_x, end_y, amount):
                dx = end_x - start_x
                dy = end_y - start_y
                dist = (dx ** 2 + dy ** 2) ** 0.5
                if dist <= 0:
                    return start_x, start_y
                return start_x + dx / dist * amount, start_y + dy / dist * amount
            for arrow in arrows:
                entry1_id = arrow.get("entry1")
                entry2_id = arrow.get("entry2")
                entry1 = find_entry_or_group(entry1_id, entries, groups) if entry1_id else None
                entry2 = find_entry_or_group(entry2_id, entries, groups) if entry2_id else None
                entry1_is_protbox = bool(entry1 and entry1.get("id") in prot_ids)
                entry2_is_protbox = bool(entry2 and entry2.get("id") in prot_ids)
                entry1_is_compound = bool(entry1 and entry1.get("type") == "compound")
                entry2_is_compound = bool(entry2 and entry2.get("type") == "compound")
                center1_x, center1_y = _entry_center(entry1)
                center2_x, center2_y = _entry_center(entry2)
                # If exactly one endpoint is a protbox, attach to the protbox handle and use coords for the other end.
                if entry1 and entry2 and entry1_is_protbox != entry2_is_protbox:
                    dx = (center2_x or 0) - (center1_x or 0)
                    dy = (center2_y or 0) - (center1_y or 0)
                    angle = calculate_angle(center1_x or 0, center1_y or 0, center2_x or 0, center2_y or 0)
                    primary_out_old = determine_arrow_side(angle)
                    primary_out = map_side(primary_out_old)
                    side_out = primary_out
                    angle_rev = np.arctan2(-dy, -dx)
                    primary_in_old = determine_arrow_side(angle_rev)
                    primary_in = map_side(primary_in_old)
                    side_in = primary_in
                    arrow_entry = {
                        'line': arrow.get('line', 'arrow'),
                        'type': arrow.get('type', ''),
                        'control_points': arrow.get('control_points', []),
                        'connector_type': arrow.get('connector_type', ''),
                    }
                    if entry1_is_protbox:
                        arrow_entry['protbox_id_1'] = entry1_id
                        arrow_entry['protbox_id_1_side'] = side_out
                        if center2_x is not None and center2_y is not None:
                            end_x, end_y = center2_x, center2_y
                            if entry2_is_compound:
                                end_x, end_y = _shorten_end(center1_x or 0, center1_y or 0, end_x, end_y, 10)
                            arrow_entry['x2'] = end_x
                            arrow_entry['y2'] = end_y
                    else:
                        arrow_entry['protbox_id_2'] = entry2_id
                        arrow_entry['protbox_id_2_side'] = side_in
                        if center1_x is not None and center1_y is not None:
                            start_x, start_y = center1_x, center1_y
                            if entry1_is_compound:
                                start_x, start_y = _shorten_start(start_x, start_y, center2_x or 0, center2_y or 0, 10)
                            arrow_entry['x1'] = start_x
                            arrow_entry['y1'] = start_y
                    json_data['arrows'].append(arrow_entry)
                    continue
                # Fallback: if either endpoint is not a protbox, use provided coords or derive from entry centers
                if not (entry1 and entry2 and entry1_is_protbox and entry2_is_protbox):
                    x1 = arrow.get("x1")
                    y1 = arrow.get("y1")
                    x2 = arrow.get("x2")
                    y2 = arrow.get("y2")
                    if (x1 is None or y1 is None) and center1_x is not None and center1_y is not None:
                        x1 = center1_x
                        y1 = center1_y
                    if (x2 is None or y2 is None) and center2_x is not None and center2_y is not None:
                        x2 = center2_x
                        y2 = center2_y
                    if x1 is not None and y1 is not None and x2 is not None and y2 is not None:
                        if entry1_is_compound:
                            x1, y1 = _shorten_start(x1, y1, x2, y2, 10)
                        if entry2_is_compound:
                            x2, y2 = _shorten_end(x1, y1, x2, y2, 10)
                        json_data['arrows'].append(
                            {
                                'x1': x1,
                                'y1': y1,
                                'x2': x2,
                                'y2': y2,
                                'line': arrow.get('line', 'arrow'),
                                'type': arrow.get('type', ''),
                                'control_points': arrow.get('control_points', []),
                                'connector_type': arrow.get('connector_type', ''),
                            }
                        )
                    continue
                if center1_x is None:
                    center1_x = entry1['x'] + entry1['width'] / 2
                if center1_y is None:
                    center1_y = entry1['y'] + entry1['height'] / 2
                if center2_x is None:
                    center2_x = entry2['x'] + entry2['width'] / 2
                if center2_y is None:
                    center2_y = entry2['y'] + entry2['height'] / 2
                dx = center2_x - center1_x
                dy = center2_y - center1_y
                angle = calculate_angle(center1_x, center1_y, center2_x, center2_y)
                possible_out = get_possible_outgoing_sides(dx, dy)
                primary_out_old = determine_arrow_side(angle)
                primary_out = map_side(primary_out_old)
                side_out = primary_out  # Use primary without choose_side

                angle_rev = np.arctan2(-dy, -dx)
                possible_in = [opposite_side[s] for s in possible_out if s in opposite_side]
                primary_in_old = determine_arrow_side(angle_rev)
                primary_in = map_side(primary_in_old)
                side_in = primary_in  # Use primary without choose_side

                arrow_entry = {
                    'protbox_id_1': arrow['entry1'],
                    'protbox_id_2': arrow['entry2'],
                    'protbox_id_1_side': side_out,
                    'protbox_id_2_side': side_in,
                    'line': arrow['line'],
                    'type': arrow['type'],
                    'control_points': arrow.get('control_points', []),
                    'connector_type': arrow.get('connector_type', ''),
                }
                json_data['arrows'].append(arrow_entry)
            if json_data.get('groups') and json_data.get('arrows'):
                entry_map = {str(entry.get('id')): entry for entry in entries if entry.get('type') == 'prot_box'}
                group_members = {}
                group_anchor = {}
                for group in json_data.get('groups', []):
                    if not isinstance(group, dict):
                        continue
                    if not group.get('show_box'):
                        continue
                    gid = str(group.get('group_id'))
                    prot_ids = []
                    if isinstance(group.get('protbox_ids'), list):
                        prot_ids.extend([str(pid) for pid in group.get('protbox_ids') if pid is not None])
                    if (not prot_ids) and isinstance(group.get('members'), list):
                        for member in group.get('members'):
                            if isinstance(member, dict) and member.get('type') == 'prot-box' and member.get('id') is not None:
                                prot_ids.append(str(member.get('id')))
                    prot_ids = [pid for pid in prot_ids if pid in entry_map]
                    if len(prot_ids) < 2:
                        continue
                    group_members[gid] = sorted(set(prot_ids))
                    if group.get('anchor_member') is not None:
                        group_anchor[gid] = str(group.get('anchor_member'))
                if group_members:
                    membership = defaultdict(set)
                    for gid, ids in group_members.items():
                        for pid in ids:
                            membership[pid].add(gid)

                    def _entry_center(entry):
                        if not entry:
                            return None, None
                        return (
                            entry.get('x', 0) + entry.get('width', 0) / 2,
                            entry.get('y', 0) + entry.get('height', 0) / 2,
                        )

                    def _pick_representative(gid, target_id):
                        anchor = group_anchor.get(gid)
                        if anchor and anchor in group_members.get(gid, []):
                            return anchor
                        target_entry = entry_map.get(str(target_id))
                        tx, ty = _entry_center(target_entry)
                        best_id = None
                        best_dist = None
                        for pid in group_members.get(gid, []):
                            entry = entry_map.get(pid)
                            cx, cy = _entry_center(entry)
                            if cx is None or cy is None or tx is None or ty is None:
                                continue
                            dist = (cx - tx) ** 2 + (cy - ty) ** 2
                            if best_dist is None or dist < best_dist:
                                best_dist = dist
                                best_id = pid
                        return best_id or (group_members.get(gid, []) or [None])[0]

                    def _make_arrow_entry(pid1, pid2, template):
                        entry1 = entry_map.get(str(pid1))
                        entry2 = entry_map.get(str(pid2))
                        if not entry1 or not entry2:
                            return None
                        center1_x = entry1['x'] + entry1['width'] / 2
                        center1_y = entry1['y'] + entry1['height'] / 2
                        center2_x = entry2['x'] + entry2['width'] / 2
                        center2_y = entry2['y'] + entry2['height'] / 2
                        dx = center2_x - center1_x
                        dy = center2_y - center1_y
                        angle = calculate_angle(center1_x, center1_y, center2_x, center2_y)
                        primary_out = map_side(determine_arrow_side(angle))
                        angle_rev = np.arctan2(-dy, -dx)
                        primary_in = map_side(determine_arrow_side(angle_rev))
                        return {
                            'protbox_id_1': pid1,
                            'protbox_id_2': pid2,
                            'protbox_id_1_side': primary_out,
                            'protbox_id_2_side': primary_in,
                            'line': template.get('line', 'arrow'),
                            'type': template.get('type', ''),
                            'control_points': [],
                            'connector_type': template.get('connector_type', ''),
                        }

                    group_links = {}
                    for idx, arrow in enumerate(list(json_data.get('arrows', []))):
                        pid1 = arrow.get('protbox_id_1')
                        pid2 = arrow.get('protbox_id_2')
                        if pid1 is None or pid2 is None:
                            continue
                        pid1 = str(pid1)
                        pid2 = str(pid2)
                        line = arrow.get('line', 'arrow')
                        rel_type = arrow.get('type', '')
                        if pid1 in membership:
                            for gid in membership[pid1]:
                                if pid2 in group_members.get(gid, []):
                                    continue
                                key = (gid, pid2, 'out', line, rel_type)
                                entry = group_links.setdefault(key, {'members': set(), 'indices': [], 'template': arrow})
                                entry['members'].add(pid1)
                                entry['indices'].append(idx)
                        if pid2 in membership:
                            for gid in membership[pid2]:
                                if pid1 in group_members.get(gid, []):
                                    continue
                                key = (gid, pid1, 'in', line, rel_type)
                                entry = group_links.setdefault(key, {'members': set(), 'indices': [], 'template': arrow})
                                entry['members'].add(pid2)
                                entry['indices'].append(idx)

                    remove_indices = set()
                    collapsed_arrows = []
                    for key, info in group_links.items():
                        gid, external_id, direction, line, rel_type = key
                        group_ids = set(group_members.get(gid, []))
                        if not group_ids or info['members'] != group_ids:
                            continue
                        indices = info['indices']
                        if not indices:
                            continue
                        if any(idx in remove_indices for idx in indices):
                            continue
                        representative = _pick_representative(gid, external_id)
                        if not representative:
                            continue
                        if direction == 'out':
                            new_arrow = _make_arrow_entry(representative, external_id, info['template'])
                        else:
                            new_arrow = _make_arrow_entry(external_id, representative, info['template'])
                        if new_arrow:
                            collapsed_arrows.append(new_arrow)
                        remove_indices.update(indices)
                    if remove_indices:
                        kept = [a for idx, a in enumerate(json_data['arrows']) if idx not in remove_indices]
                        kept.extend(collapsed_arrows)
                        json_data['arrows'] = kept
            # Populate protein_data_map for all entries first
            for entry in entries:
                entry_type = entry.get("type", "")
                if entry_type != 'prot_box':
                    continue
                proteins = entry["name"].split()
                self.protein_data_map[entry["id"]] = {
                    'proteins': proteins,
                    'x': entry["x"],
                    'y': entry["y"],
                    'width': entry["width"],
                    'height': entry["height"],
                    'first_name': entry.get("first_name", entry["name"].split(",")[0].strip()),
                    'label': entry.get("label") or entry.get("backup_label") or entry.get("first_name"),
                }
            protbox_protein_cache = {}
            all_uniprot_ids = set()
            for entry in entries:
                if entry.get("type", "") != 'prot_box':
                    continue
                proteins = entry["name"].split()
                valid_protein_dicts = []
                uniprot_ids = []
                for protein in proteins:
                    protein_payload = self._build_protein_match_payload(protein)
                    if not protein_payload:
                        continue
                    uniprot_id = protein_payload.get('uniprot_id')
                    if not uniprot_id:
                        continue
                    valid_protein_dicts.append(protein_payload)
                    uniprot_ids.append(uniprot_id)
                    all_uniprot_ids.add(uniprot_id)
                protbox_protein_cache[entry["id"]] = {
                    'proteins': proteins,
                    'valid_protein_dicts': valid_protein_dicts,
                    'uniprot_ids': uniprot_ids
                }
            # Prefetch PTM summaries once after collecting all UniProt IDs
            self._prefetch_ptm_summaries(all_uniprot_ids)
            # Process protein boxes
            for entry in entries:
                entry_type = entry.get("type", "")
                if entry_type != 'prot_box':
                    continue
                cache_payload = protbox_protein_cache.get(entry["id"], {})
                proteins = cache_payload.get('proteins', entry["name"].split())
                protbox_ptm_overrides = {}
                valid_protein_dicts = list(cache_payload.get('valid_protein_dicts', []))
                uniprot_ids = list(cache_payload.get('uniprot_ids', []))
                use_original_size = bool(self.settings.get("use_original_protbox_size", True))
                if valid_protein_dicts:
                    # Determine the chosen (highest priority) protein
                    chosen_protein = self.choose_protein(proteins)
                    if chosen_protein and chosen_protein['uniprot_id']:
                        chosen_uniprot = chosen_protein['uniprot_id']
                        # Reorder uniprot_ids with chosen first
                        uniprot_ids = [chosen_uniprot] + [u for u in uniprot_ids if u != chosen_uniprot]
                        # Process all valid proteins
                        for protein_dict in valid_protein_dicts:
                            protein_entry = self.process_protein_box(entry["id"], protein_dict, proximities, protein_to_entry_ids, exceptions)
                            json_data['protein_data'][protein_dict['uniprot_id']] = protein_entry
                            snapshot = self._extract_ptm_override_block(protein_entry.get('PTMs', {}))
                            if snapshot and protein_dict.get('uniprot_id'):
                                protbox_ptm_overrides[protein_dict['uniprot_id']] = snapshot
                            chosen_proteins.add(protein_dict['protein'])
                        self.current_proteins[entry["id"]] = chosen_protein
                    else:
                        # Fallback if no chosen
                        for protein_dict in valid_protein_dicts:
                            protein_entry = self.process_protein_box(entry["id"], protein_dict, proximities, protein_to_entry_ids, exceptions)
                            json_data['protein_data'][protein_dict['uniprot_id']] = protein_entry
                            snapshot = self._extract_ptm_override_block(protein_entry.get('PTMs', {}))
                            if snapshot and protein_dict.get('uniprot_id'):
                                protbox_ptm_overrides[protein_dict['uniprot_id']] = snapshot
                            chosen_proteins.add(protein_dict['protein'])
                        self.current_proteins[entry["id"]] = valid_protein_dicts[0]
                else:
                    # No valid proteins
                    fallback_protein = proteins[0] if proteins else "Unknown"
                    self.current_proteins[entry["id"]] = {'protein': fallback_protein, 'uniprot_id': None, 'annotations': []}
                    protein_entry = self.process_protein_box(entry["id"], {'protein': fallback_protein, 'uniprot_id': None, 'annotations': []}, proximities, protein_to_entry_ids, exceptions)
                    json_data['protein_data'][f"unknown_{entry['id']}"] = protein_entry
                is_kegg = str(self.settings.get("pathway_source", "kegg")).lower() == "kegg"
                center_x = entry["x"] if is_kegg else entry["x"] + entry["width"] / 2
                center_y = entry["y"] if is_kegg else entry["y"] + entry["height"] / 2
                protbox_entry = {
                    'protbox_id': entry["id"],
                    'proteins': uniprot_ids,
                    'backup_label': self.protein_data_map[entry["id"]]['first_name']
                }
                if use_original_size:
                    width_val = entry["width"]
                    height_val = entry["height"]
                else:
                    width_val = 46
                    height_val = 17
                protbox_entry.update({
                    'x': center_x - width_val / 2,
                    'y': center_y - height_val / 2,
                    'width': width_val,
                    'height': height_val,
                })
                # Attach tooltip from the primary protein if available
                if uniprot_ids:
                    primary = uniprot_ids[0]
                    primary_payload = json_data['protein_data'].get(primary, {})
                    protbox_entry['tooltip'] = primary_payload.get('tooltip', '')
                    protbox_entry['tooltip_html'] = primary_payload.get('tooltip_html', '')
                if protbox_ptm_overrides:
                    protbox_entry['ptm_overrides'] = protbox_ptm_overrides
                json_data['protbox_data'].append(protbox_entry)
            for entry in entries:
                if entry.get('type') == 'compound':
                    # Normalize compound size to a small fixed circle-like footprint, but keep the original center.
                    pathway_source = str(self.settings.get('pathway_source', 'kegg')).lower()
                    # KEGG graphics x/y are already centers; WikiPathways entries store top-left, so adjust accordingly.
                    if pathway_source == 'kegg':
                        center_x = entry['x']
                        center_y = entry['y']
                    else:
                        center_x = entry['x'] + entry['width'] / 2
                        center_y = entry['y'] + entry['height'] / 2
                    # Match KEGG default compound size (8x8 circle); viewer treats x/y as center
                    default_size = 8
                    metabolite_match = self._resolve_metabolite_match(entry)
                    tooltip, tooltip_html = self._build_metabolite_tooltip(metabolite_match)
                    display_label = (
                        metabolite_match.get('display_label')
                        if metabolite_match else
                        (entry.get('first_name', '') or entry.get('label') or '')
                    )
                    compound = {
                        'compound_id': entry['id'],  # KGML numeric id
                        'kegg_compound': entry['name'],  # e.g., 'cpd:C00165'
                        'label': display_label,
                        'x': center_x,
                        'y': center_y,
                        'width': default_size,
                        'height': default_size,
                        'fgcolor': entry.get('fgcolor', '#000000'),
                        'bgcolor': entry.get('bgcolor', '#FFFFFF'),
                        'graphics_type': entry.get('graphics_type', ''),  # likely 'circle'
                        'link': entry.get('link', ''),
                        'tooltip': tooltip,
                        'tooltip_html': tooltip_html,
                    }
                    if metabolite_match:
                        compound['hmdb_id'] = metabolite_match.get('hmdb_id', '')
                        compound['pathway_match_id'] = metabolite_match.get('pathway_match_id', '')
                        row_map = metabolite_match.get('row_map') or {}
                        for idx, fc_col in enumerate(self.metabolite_main_columns, 1):
                            fc_value = self._safe_float(row_map.get(fc_col))
                            compound[f'fold_change_{idx}'] = fc_value
                            compound[f'fc_color_{idx}'] = self.get_color(fc_value) if fc_value is not None else [128, 128, 128]
                    json_data['compound_data'].append(compound)
            for entry in entries:
                entry_type = entry.get('type')
                gfx_type = entry.get('graphics_type')
                has_shape = bool(entry.get('shape_type'))
                if entry_type in {'map', 'label', 'text'} or gfx_type in {'roundrectangle', 'label'} or has_shape:
                    label = entry.get('first_name', '')
                    # Optionally strip 'TITLE:' prefix if present
                    if label.upper().startswith('TITLE:'):
                        label = label.split(':', 1)[1].strip()

                    bgcolor = entry.get('bgcolor', '#FFFFFF')
                    # Make default white backgrounds transparent so labels appear see-through
                    if isinstance(bgcolor, str) and bgcolor.lower().lstrip('#') in {'ffffff', 'fff'}:
                        bgcolor = 'transparent'
                    border_color = entry.get('border_color', 'transparent') or 'transparent'
                    border_width = entry.get('border_width', 0) or 0
                    is_kegg = str(self.settings.get('pathway_source', 'kegg')).lower() == 'kegg'
                    width = entry.get('width', 0)
                    height = entry.get('height', 0)
                    x_val = entry.get('x', 0)
                    y_val = entry.get('y', 0)
                    if is_kegg:
                        x_val = x_val - (width / 2)
                        y_val = y_val - (height / 2)
                    text_item = {
                        'text_id': entry['id'],
                        'label': label,
                        'x': x_val,
                        'y': y_val,
                        'width': width,
                        'height': height,
                        'fgcolor': entry.get('fgcolor', '#000000'),
                        'bgcolor': bgcolor,
                        'graphics_type': entry.get('graphics_type', ''),  # often 'roundrectangle'
                        'text_style': entry.get('text_style', {}),
                        'border_color': border_color,
                        'border_width': border_width,
                        'is_background': entry.get('is_background', False),
                        'shape_type': entry.get('shape_type', ''),
                        'link': entry.get('link', '')
                    }
                    json_data['text_data'].append(text_item)

            if not skip_disk_write:
                if not _debug_file_output_enabled():
                    print("Skipping pathway data save to disk (MAPKINASE_DEBUG_FILE_OUTPUT is disabled).")
                else:
                    debug_output_dir = _resolve_debug_output_dir(self.settings)
                    if debug_output_dir is None:
                        print("Skipping pathway data save to disk (no debug output directory configured).")
                    else:
                        temp_file = savefile(
                            self.settings.get('pathway_id', 'hsa04010'),
                            str(debug_output_dir),
                            '_visprots',
                            extension='.txt',
                            include_timestamp=False,
                        )
                        try:
                            with open(temp_file, 'w', encoding='utf-8') as f:
                                for protein in chosen_proteins:
                                    f.write(f"{protein}\n")
                        except (IOError, PermissionError) as e:
                            print(f"Error writing to temporary file {temp_file}: {e}")
                        output_file = savefile(
                            self.settings.get('pathway_id', 'hsa04010'),
                            str(debug_output_dir),
                            '_pathway_data',
                            extension='.json',
                        )
                        with open(output_file, 'w', encoding='utf-8') as f:
                            json.dump(json_data, f, indent=4)
                        print(f"Pathway data saved to {output_file}")
            else:
                print("Skipping pathway data save to disk (skip_disk_write=True)")
            return json_data
        except Exception as e:
            print(f"Error in process_pathway: {e}")
            return None

    def process_protein_box(self, entry_id, protein_dict, proximities, protein_to_entry_ids, exceptions):
        try:
            entry_data = self.protein_data_map[entry_id]
            use_original_size = bool(self.settings.get("use_original_protbox_size", True))
            is_kegg = str(self.settings.get("pathway_source", "kegg")).lower() == "kegg"
            center_x = entry_data['x'] if is_kegg else entry_data['x'] + entry_data['width'] / 2
            center_y = entry_data['y'] if is_kegg else entry_data['y'] + entry_data['height'] / 2
            if use_original_size:
                width = entry_data['width']
                height = entry_data['height']
            else:
                width = 46
                height = 17
            x = center_x - width / 2
            y = center_y - height / 2
            proteins = entry_data['proteins']
            protein = protein_dict['protein']
            uniprot_id = protein_dict['uniprot_id']
            annotations = protein_dict.get('annotations', [])
            protein_row = protein_dict.get('row')
            protein_in_data = protein_row is not None
            valid_proteins = [p for p in proteins if self._build_protein_match_payload(p) is not None]
            protein_entry = {
                'label': '',
                'label_color': [0, 0, 0],
                'transcriptomic_color': [],
                'annotations': ','.join(f'"{ann}"' for ann in annotations if ann),
                'tooltip': '',
                'tooltip_html': '',
                'PTMs': {}
            }
            for idx, fc_col in enumerate(self.fold_change_columns, 1):
                protein_entry[f'fc_color_{idx}'] = [128, 128, 128]
                protein_entry[f'fold_change_{idx}'] = None
                protein_entry[f'outline_color_{idx}'] = [0, 0, 0]
                protein_entry[f'outline_fold_change_{idx}'] = None
            if protein_in_data:
                gene_name = protein_row.get(self.gene_name_column, protein)
                protein_entry['label'] = gene_name
                outline_cols = self.outline_columns or []
                use_black_outlines = bool(self.settings.get('use_black_protein_outlines', False))
                for idx, fc_col in enumerate(self.fold_change_columns, 1):
                    raw_value = protein_row.get(fc_col) if fc_col in self.proteomic_data.columns else None
                    fc_value = None
                    if raw_value is not None and not pd.isna(raw_value):
                        try:
                            fc_value = float(raw_value)
                        except (TypeError, ValueError):
                            fc_value = None
                    protein_entry[f'fold_change_{idx}'] = fc_value
                    if fc_value is not None:
                        protein_entry[f'fc_color_{idx}'] = self.get_color(fc_value)
                    else:
                        protein_entry[f'fc_color_{idx}'] = [128, 128, 128]
                    outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                    outline_value = None
                    if outline_col and outline_col in self.proteomic_data.columns:
                        outline_value = self._safe_float(protein_row.get(outline_col))
                    protein_entry[f'outline_fold_change_{idx}'] = outline_value
                    if use_black_outlines:
                        protein_entry[f'outline_color_{idx}'] = [0, 0, 0]
                    elif outline_value is not None:
                        protein_entry[f'outline_color_{idx}'] = self.get_color(outline_value)
                    else:
                        protein_entry[f'outline_color_{idx}'] = [0, 0, 0]
            else:
                gene_name = entry_data['first_name'].split(",")[0].strip()
                protein_entry['label'] = gene_name
            # Build protein tooltip from configured columns (skip blanks)
            tooltip_plain_lines = []
            tooltip_html_lines = []
            if protein_in_data:
                for col in self.protein_tooltip_columns:
                    label = str(col)
                    value = ''
                    if col in self.proteomic_data.columns:
                        cell = protein_row.get(col, '')
                        if cell is not None and not pd.isna(cell):
                            value = str(cell)
                    if str(value).strip() == '':
                        continue
                    tooltip_plain_lines.append(f"{label}: {value}")
                    safe_label = html.escape(label)
                    safe_value = html.escape(value)
                    tooltip_html_lines.append(f"<strong>{safe_label}:</strong> {safe_value}")
            protein_entry['tooltip'] = '\n'.join(tooltip_plain_lines)
            protein_entry['tooltip_html'] = '<br/>'.join(tooltip_html_lines)
            if protein == 'MAPK1' or 'MAPK1' in proteins:
                print(
                    f"Entry {entry_id}: Proteins = {proteins}, Chosen = {protein}, UniProt = {uniprot_id}, Valid = {valid_proteins}")
            if protein_in_data and uniprot_id:
                ptm_sites = self.prioritize_ptm_sites(uniprot_id)
                if not ptm_sites.empty:
                    if protein == 'MAPK1':
                        print(ptm_sites[[dataset['site_column'] for dataset in self.ptm_datasets.values() if
                                         'site_column' in dataset] + ['dataset_id', 'Phosphosite_Classification']])
                    positions = {
                        'N1': (x + width * 0.2, y - self.ptm_circle_spacing),
                        'N2': (x + width * 0.5, y - self.ptm_circle_spacing),
                        'N3': (x + width * 0.8, y - self.ptm_circle_spacing),
                        'S1': (x + width * 0.2, y + height + self.ptm_circle_spacing),
                        'S2': (x + width * 0.5, y + height + self.ptm_circle_spacing),
                        'S3': (x + width * 0.8, y + height + self.ptm_circle_spacing),
                        'W1': (x - self.ptm_circle_spacing, y + height * 0.33),
                        'W2': (x - self.ptm_circle_spacing, y + height * 0.66),
                        'E1': (x + width + self.ptm_circle_spacing, y + height * 0.33),
                        'E2': (x + width + self.ptm_circle_spacing, y + height * 0.66)
                    }
                    position_priority = ['N1', 'N2', 'N3', 'S1', 'S2', 'S3', 'W1', 'W2', 'E1', 'E2']
                    pathway_id = self.settings.get('pathway_id', 'hsa04010')
                    hsa_exceptions = exceptions.get('hsa_specific', {}).get(pathway_id,
                                                                            {'proximity': [], 'specific': {}})
                    global_exceptions = exceptions.get('global', {'proximity': [], 'specific': {}})
                    entry_proximities = proximities.get(entry_id, {})
                    if not entry_proximities:
                        print(f"Warning: No proximity data for entry {entry_id} (protein: {protein}).")
                    specific_rules = hsa_exceptions['specific'].get(protein, []) + global_exceptions['specific'].get(
                        protein, [])
                    blocked_positions = set()
                    custom_priority = None
                    move_circle = {}
                    move_label = {}
                    for rule in specific_rules:
                        action = rule[0].lower()
                        if action == 'block' and len(rule) > 1:
                            blocked_positions.add(rule[1])
                        elif action == 'priority':
                            custom_priority = rule[1:]
                        elif action == 'move_circle' and len(rule) == 4:
                            pos, dx, dy = rule[1], float(rule[2]), float(rule[3])
                            move_circle[pos] = (dx, dy)
                        elif action == 'move_label' and len(rule) == 5:
                            pos, dx, dy, align = rule[1], float(rule[2]), float(rule[3]), rule[4]
                            move_label[pos] = (dx, dy, align)
                    all_proximity_rules = global_exceptions['proximity'] + hsa_exceptions['proximity']
                    proximity_rules = sorted(all_proximity_rules,
                                             key=lambda r: (r['priority'], all_proximity_rules.index(r)))
                    passed_rules_by_action = {'block': [], 'priority': [], 'move_circle': [], 'move_label': []}
                    for rule in proximity_rules:
                        conditions_satisfied = {cond: False for cond, _ in rule['conditions']}
                        for other_entry_id, distances in entry_proximities.items():
                            dx = distances['dx']
                            dy = distances['dy']
                            for condition, threshold in rule['conditions']:
                                if condition == 'adjacent_y_south' and dy > 0 and abs(dy) < threshold:
                                    if abs(dx) <= 23:
                                        conditions_satisfied[condition] = True
                                elif condition == 'adjacent_y_north' and dy < 0 and abs(dy) < threshold:
                                    if abs(dx) <= 23:
                                        conditions_satisfied[condition] = True
                                elif condition == 'adjacent_x_west' and dx < 0 and abs(dx) < threshold:
                                    if abs(dy) <= 8:
                                        conditions_satisfied[condition] = True
                                elif condition == 'adjacent_x_east' and dx > 0 and abs(dx) < threshold:
                                    if abs(dy) <= 8:
                                        conditions_satisfied[condition] = True
                        if all(conditions_satisfied.values()):
                            passed_rules_by_action[rule['action']].append(rule)
                    for action in ['block', 'priority']:
                        if passed_rules_by_action[action]:
                            highest_rule = max(passed_rules_by_action[action], key=lambda r: r['priority'])
                            if action == 'block':
                                for pos in highest_rule['values']:
                                    blocked_positions.add(pos)
                            elif action == 'priority':
                                custom_priority = highest_rule['values']
                    for action in ['move_circle', 'move_label']:
                        for rule in sorted(passed_rules_by_action[action],
                                           key=lambda r: (r['priority'], proximity_rules.index(r))):
                            if action == 'move_circle' and len(rule['values']) == 3:
                                pos, dx, dy = rule['values'][0], float(rule['values'][1]), float(rule['values'][2])
                                move_circle[pos] = (dx, dy)
                            elif action == 'move_label' and len(rule['values']) == 4:
                                pos, dx, dy, align = rule['values'][0], float(rule['values'][1]), float(
                                    rule['values'][2]), rule['values'][3]
                                move_label[pos] = (dx, dy, align)
                    available_positions = [p for p in (custom_priority or position_priority) if
                                           p not in blocked_positions]
                    defaults = {
                        'N1': (-5, -5, 'right'), 'N2': (0, -8, 'center'), 'N3': (5, -5, 'left'),
                        # yaxis N2 = -11 for Program, -8 for website
                        # xaxis N1/3 = -3 (or 3) for Program, and -5 (or 5) for Website
                        'S1': (-5, 5, 'right'), 'S2': (0, 10, 'center'), 'S3': (5, 5, 'left'),
                        # yaxis S2 = 12 for Program, 10 for website
                        # xaxis S1/3 = -3 (or 3) for Program, and -5 (or 5) for Website
                        'W1': (-5, -2, 'right'), 'W2': (-5, 2, 'right'),  # x = -3 for Program, -5 for website
                        'E1': (5, -2, 'left'), 'E2': (5, 2, 'left')     # x = 3 for Program, 5 for website
                    }
                    max_display = self.settings.get('ptm_max_display', 4)
                    num_to_visualize = min(len(available_positions), max_display)
                    for i, (_, row) in enumerate(ptm_sites.iterrows()):
                        dataset_id = row['dataset_id']
                        dataset = self.ptm_datasets[dataset_id]
                        ptm_shape = dataset.get('shape', 'Circle')
                        classification = row.get('Phosphosite_Classification', 'none')
                        site_position = str(row[dataset['site_column']]) if dataset['site_column'] in row else ''
                        ptm_label = site_position
                        ptm_key = f"{uniprot_id}_{ptm_label}" if uniprot_id and ptm_label else f"ptm_{i}"

                        visualize = i < num_to_visualize

                        # Common fields for all PTMs
                        ptm_entry = {
                            'ptm_type': dataset_id,
                            'uniprot_id': uniprot_id,
                            'fc_color_1': [128, 128, 128],
                            'label': ptm_label,
                            'label_color': list(self.ptm_label_color),
                            'symbol_type': classification,
                            'annotated': '+' if classification != 'none' else '',
                            'shape': ptm_shape,
                            'symbol': '',
                            'symbol_color': [0, 0, 0],
                            'symbol_font': 'Arial',
                            'symbol_size': 6,
                        }
                        # Calculate tooltip values
                        tooltip_plain_lines = []
                        tooltip_html_lines = []
                        for col in dataset.get('tooltip_columns', []):
                            label = str(col)
                            if col in row:
                                raw_value = row[col]
                                value = '' if pd.isna(raw_value) else str(raw_value)
                            else:
                                print(f"Warning: Tooltip column {col} not found in PTM dataset {dataset_id}")
                                value = ''
                            if value.strip() == '':
                                continue  # skip empty values entirely
                            tooltip_plain_lines.append(f"{label}: {value}")
                            safe_label = html.escape(label)
                            safe_value = html.escape(value)
                            tooltip_html_lines.append(f"<strong>{safe_label}:</strong> {safe_value}")
                        ptm_entry['tooltip'] = '\n'.join(tooltip_plain_lines)
                        ptm_entry['tooltip_html'] = '<br/>'.join(tooltip_html_lines)
                        outline_cols = dataset.get('outline_columns') or []
                        for idx, (_, fc_col) in enumerate(dataset['main_columns'], 1):
                            ptm_entry.setdefault(f'fc_color_{idx}', [128, 128, 128])
                            ptm_entry[f'fold_change_{idx}'] = None
                            ptm_entry.setdefault(f'outline_color_{idx}', [0, 0, 0])
                            ptm_entry[f'outline_fold_change_{idx}'] = None
                            ptm_fold_change = row.get(fc_col)
                            fc_value = None
                            if ptm_fold_change is not None and not pd.isna(ptm_fold_change):
                                try:
                                    fc_value = float(ptm_fold_change)
                                except (TypeError, ValueError):
                                    fc_value = None
                            ptm_entry[f'fold_change_{idx}'] = fc_value
                            if fc_value is not None:
                                ptm_entry[f'fc_color_{idx}'] = self.get_color(fc_value)
                            outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                            outline_value = self._safe_float(row.get(outline_col)) if outline_col else None
                            ptm_entry[f'outline_fold_change_{idx}'] = outline_value
                            if outline_value is not None:
                                ptm_entry[f'outline_color_{idx}'] = self.get_color(outline_value)

                        if visualize:
                            pos_key = available_positions[i]
                            ptm_entry['ptm_position'] = pos_key

                            # Calculate base position for shape
                            pos_x, pos_y = positions.get(pos_key, (x, y))
                            shape_x, shape_y = pos_x, pos_y

                            # Apply move_circle adjustments to shape only
                            if pos_key in move_circle:
                                circle_dx, circle_dy = move_circle[pos_key]
                                shape_x += circle_dx
                                shape_y += circle_dy

                            # Apply label positioning
                            if pos_key in move_label:
                                dx, dy, align = move_label[pos_key]
                                label_x = shape_x + dx
                                label_y = shape_y + dy
                                label_centering = align.lower()
                            else:
                                x_offset, y_offset, label_centering = defaults[pos_key]
                                label_x = shape_x + x_offset
                                label_y = shape_y + y_offset

                            ptm_entry['shape_x'] = float(shape_x)
                            ptm_entry['shape_y'] = float(shape_y)
                            ptm_entry['label_x'] = float(label_x)
                            ptm_entry['label_y'] = float(label_y)
                            ptm_entry['label_centering'] = label_centering

                            # Calculate symbol positioning
                            symbol_type = classification if classification != 'none' else ''
                            symbol = ''
                            symbol_x = shape_x
                            symbol_y = shape_y
                            symbol_color = [0, 0, 0]
                            symbol_font = 'Arial'
                            symbol_size = 6
                            symbol_icon = None
                            if symbol_type:
                                try:
                                    ptm_type_idx = int(dataset_id.split('_')[1])
                                    ptm_datasets = self.ptm_datasets
                                    if ptm_type_idx < len(ptm_datasets):
                                        symbol_dict = next(
                                            (item[f"{symbol_type}_dict"] for item in ptm_datasets[dataset_id].get('ptm_symbol_list', []) if f"{symbol_type}_dict" in item),
                                            None
                                        )
                                        if symbol_dict:
                                            symbol = symbol_dict.get('symbol', '')
                                            symbol_color = symbol_dict.get('symbol_color', [0, 0, 0])
                                            # Anchor the symbol to the PTM shape center rather than per-side offsets.
                                            symbol_x = shape_x
                                            symbol_y = shape_y
                                            symbol_font = symbol_dict.get('symbol_font', 'Arial')
                                            symbol_size = symbol_dict.get('symbol_size', 6)
                                            symbol_icon = _resolve_symbol_icon(symbol_dict)
                                            _safe_debug_print(f"PTM {ptm_key}: symbol={symbol}, symbol_x={symbol_x}, symbol_y={symbol_y}, pos_key={pos_key}")
                                        else:
                                            _safe_debug_print(f"Warning: No symbol dict for {symbol_type} in dataset {dataset_id}")
                                    else:
                                        _safe_debug_print(f"Warning: ptm_type_idx {ptm_type_idx} out of range for ptm_datasets")
                                except (IndexError, ValueError) as e:
                                    _safe_debug_print(f"Error processing symbol for PTM {ptm_key}: {e}")
                            ptm_entry['symbol'] = symbol
                            ptm_entry['symbol_icon'] = symbol_icon
                            ptm_entry['symbol_x'] = float(symbol_x)
                            ptm_entry['symbol_y'] = float(symbol_y)
                            ptm_entry['symbol_color'] = symbol_color
                            ptm_entry['symbol_font'] = symbol_font
                            ptm_entry['symbol_size'] = symbol_size
                            _safe_debug_print(f"PTM {ptm_key}: shape=({shape_x}, {shape_y}), label=({label_x}, {label_y}), centering={label_centering}")
                        else:
                            ptm_entry['ptm_position'] = ""

                        protein_entry['PTMs'][ptm_key] = ptm_entry
            return protein_entry
        except Exception as e:
            print(f"Error processing protein box for {entry_id}: {e}")
            return {
                'label': entry_data['first_name'].split(",")[0].strip(),
                'label_color': [0, 0, 0],
                'transcriptomic_color': [],
                'annotations': '',
                'PTMs': {}
            }

    def _safe_float(self, value):
        if value is None or pd.isna(value):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _extract_ptm_override_block(self, ptm_dict):
        overrides = {}
        if not isinstance(ptm_dict, dict):
            return overrides
        for key, payload in ptm_dict.items():
            if not isinstance(payload, dict):
                continue
            override_entry = {
                'shape_x': payload.get('shape_x'),
                'shape_y': payload.get('shape_y'),
                'label_x': payload.get('label_x'),
                'label_y': payload.get('label_y'),
                'label_centering': payload.get('label_centering'),
                'symbol_x': payload.get('symbol_x'),
                'symbol_y': payload.get('symbol_y'),
                'ptm_position': payload.get('ptm_position'),
            }
            overrides[key] = override_entry
        # Filter out entries that provided no positional data
        return {
            key: value
            for key, value in overrides.items()
            if any(coord is not None for coord in value.values())
        }

    def _prefetch_ptm_summaries(self, uniprot_ids):
        valid_ids = sorted({uid for uid in uniprot_ids if uid})
        if not valid_ids:
            return
        cache = getattr(self, "_ptm_summary_cache", {})
        self._ptm_summary_cache = cache
        missing = [uid for uid in valid_ids if uid not in cache]
        if not missing:
            return
        worker_count = min(max(1, self._ptm_worker_count), len(missing))
        if worker_count <= 1:
            for uid in missing:
                cache[uid] = self._build_ptm_summary(uid)
            return
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            future_to_uid = {executor.submit(self._build_ptm_summary, uid): uid for uid in missing}
            for future in as_completed(future_to_uid):
                uid = future_to_uid[future]
                try:
                    cache[uid] = future.result()
                except Exception as exc:
                    print(f"Warning: PTM worker failed for {uid}: {exc}")
                    cache[uid] = {}

    def _build_ptm_summary(self, uniprot_id):
        if not uniprot_id:
            return {}
        cache = getattr(self, "_ptm_summary_cache", None)
        if cache is not None and uniprot_id in cache:
            return cache[uniprot_id]
        ptm_entries = {}
        try:
            ptm_sites = self.prioritize_ptm_sites(uniprot_id)
        except Exception as exc:
            print(f"Warning: unable to build PTM summary for {uniprot_id}: {exc}")
            return ptm_entries
        if ptm_sites.empty:
            if cache is not None:
                cache[uniprot_id] = ptm_entries
            return ptm_entries
        max_display = self.settings.get('ptm_max_display', 4)
        for i, (_, row) in enumerate(ptm_sites.iterrows()):
            if i >= max_display:
                break
            dataset_id = row.get('dataset_id')
            dataset = self.ptm_datasets.get(dataset_id)
            if not dataset:
                continue
            site_column = dataset.get('site_column')
            ptm_label = ''
            if site_column and site_column in row and not pd.isna(row[site_column]):
                ptm_label = str(row[site_column])
            ptm_key = f"{uniprot_id}_{ptm_label}" if ptm_label else f"{uniprot_id}_ptm_{i}"
            ptm_entry = {
                'ptm_type': dataset_id,
                'uniprot_id': uniprot_id,
                'shape': dataset.get('shape', 'Circle'),
                'label': ptm_label,
                'label_color': list(self.ptm_label_color),
                'symbol_type': row.get('Phosphosite_Classification', 'none'),
                'annotated': '+' if str(row.get('Phosphosite_Classification', 'none')).lower() != 'none' else '',
                'symbol': '',
                'symbol_icon': None,
                'symbol_color': [0, 0, 0],
                'symbol_font': 'Arial',
                'symbol_size': 6,
                'tooltip': '',
                'tooltip_html': ''
            }
            tooltip_plain_lines = []
            tooltip_html_lines = []
            for col in dataset.get('tooltip_columns', []):
                label = str(col)
                if col in row and not pd.isna(row[col]):
                    value = str(row[col])
                else:
                    value = ''
                if value.strip() == '':
                    continue
                tooltip_plain_lines.append(f"{label}: {value}")
                safe_label = html.escape(label)
                safe_value = html.escape(value)
                tooltip_html_lines.append(f"<strong>{safe_label}:</strong> {safe_value}")
            ptm_entry['tooltip'] = '\n'.join(tooltip_plain_lines)
            ptm_entry['tooltip_html'] = '<br/>'.join(tooltip_html_lines)
            outline_cols = dataset.get('outline_columns') or []
            for idx, (_, fc_col) in enumerate(dataset.get('main_columns', []), 1):
                ptm_entry.setdefault(f'fc_color_{idx}', [128, 128, 128])
                ptm_entry[f'fold_change_{idx}'] = None
                ptm_entry.setdefault(f'outline_color_{idx}', [0, 0, 0])
                ptm_entry[f'outline_fold_change_{idx}'] = None
                fc_value = self._safe_float(row.get(fc_col))
                ptm_entry[f'fold_change_{idx}'] = fc_value
                if fc_value is not None:
                    ptm_entry[f'fc_color_{idx}'] = self.get_color(fc_value)
                outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                outline_value = self._safe_float(row.get(outline_col)) if outline_col else None
                ptm_entry[f'outline_fold_change_{idx}'] = outline_value
                if outline_value is not None:
                    ptm_entry[f'outline_color_{idx}'] = self.get_color(outline_value)
            symbol_type = ptm_entry['symbol_type']
            if symbol_type and str(symbol_type).lower() != 'none':
                symbol_dict = next(
                    (item.get(f"{symbol_type}_dict") for item in dataset.get('ptm_symbol_list', []) if f"{symbol_type}_dict" in item),
                    None
                )
                if symbol_dict:
                    ptm_entry['symbol'] = symbol_dict.get('symbol', '')
                    ptm_entry['symbol_icon'] = _resolve_symbol_icon(symbol_dict)
                    ptm_entry['symbol_color'] = symbol_dict.get('symbol_color', [0, 0, 0])
                    ptm_entry['symbol_font'] = symbol_dict.get('symbol_font', 'Arial')
                    ptm_entry['symbol_size'] = symbol_dict.get('symbol_size', 6)
            ptm_entries[ptm_key] = ptm_entry
        if cache is not None:
            cache[uniprot_id] = ptm_entries
        return ptm_entries

    def build_full_protein_catalog(self):
        catalog = {}
        seen = set()
        for _, row in self.proteomic_data.iterrows():
            if self.prot_uniprot_column not in row.index or pd.isna(row[self.prot_uniprot_column]):
                continue
            uniprot_id = str(row[self.prot_uniprot_column]).strip()
            if not uniprot_id or uniprot_id in seen:
                continue
            seen.add(uniprot_id)
            label_value = row.get(self.gene_name_column, '')
            label = str(label_value).strip() if label_value is not None and not pd.isna(label_value) else uniprot_id
            annotations = []
            for col in self.protein_tooltip_columns:
                value = row.get(col, '')
                if value is None or pd.isna(value):
                    annotations.append('')
                else:
                    annotations.append(str(value))
            protein_entry = {
                'label': label,
                'gene_symbol': label,
                'label_color': [0, 0, 0],
                'transcriptomic_color': [],
                'annotations': ','.join(f'\"{ann}\"' for ann in annotations if ann),
                'PTMs': self._build_ptm_summary(uniprot_id)
            }
            outline_cols = self.outline_columns or []
            use_black_outlines = bool(self.settings.get('use_black_protein_outlines', False))
            for idx, fc_col in enumerate(self.fold_change_columns, 1):
                protein_entry[f'fc_color_{idx}'] = [128, 128, 128]
                protein_entry[f'fold_change_{idx}'] = None
                protein_entry[f'outline_color_{idx}'] = [0, 0, 0]
                protein_entry[f'outline_fold_change_{idx}'] = None
                fc_value = self._safe_float(row.get(fc_col))
                protein_entry[f'fold_change_{idx}'] = fc_value
                if fc_value is not None:
                    protein_entry[f'fc_color_{idx}'] = self.get_color(fc_value)
                outline_col = outline_cols[idx - 1] if idx - 1 < len(outline_cols) else None
                outline_value = self._safe_float(row.get(outline_col)) if outline_col else None
                protein_entry[f'outline_fold_change_{idx}'] = outline_value
                if use_black_outlines:
                    protein_entry[f'outline_color_{idx}'] = [0, 0, 0]
                elif outline_value is not None:
                    protein_entry[f'outline_color_{idx}'] = self.get_color(outline_value)
            catalog[uniprot_id] = protein_entry
        return catalog

    def build_full_metabolite_catalog(self):
        catalog = {}
        if self.metabolite_data.empty or self.metabolite_hmdb_column not in self.metabolite_data.columns:
            return catalog
        seen = set()
        for _, row in self.metabolite_data.iterrows():
            hmdb_raw = row.get(self.metabolite_hmdb_column, '')
            hmdb_key = _normalize_metabolite_token(hmdb_raw)
            if not hmdb_key or hmdb_key in seen:
                continue
            seen.add(hmdb_key)
            ref = self.hmdb_reference_map.get(hmdb_key)
            if not ref:
                continue
            display_label = _display_metabolite_name(ref.get('wikipedia_id') or ref.get('name') or hmdb_raw)
            row_map = row.to_dict()
            row_map = {str(k): ('' if pd.isna(v) else v) for k, v in row_map.items()}
            match_payload = {
                'hmdb_id': str(ref.get('accession') or hmdb_raw or '').strip(),
                'display_label': display_label,
                'pathway_ids': [],
                'reference': ref,
                'row_map': row_map,
            }
            tooltip, tooltip_html = self._build_metabolite_tooltip(match_payload)
            entry = {
                'hmdb_id': str(ref.get('accession') or hmdb_raw or '').strip(),
                'display_label': display_label,
                'wikipedia_id': _display_metabolite_name(ref.get('wikipedia_id') or ''),
                'kegg_id': str(ref.get('kegg_id') or '').strip(),
                'name': str(ref.get('name') or '').strip(),
                'synonyms': str(ref.get('synonyms') or '').strip(),
                'tooltip': tooltip,
                'tooltip_html': tooltip_html,
            }
            for idx, fc_col in enumerate(self.metabolite_main_columns, 1):
                fc_value = self._safe_float(row_map.get(fc_col))
                entry[f'fold_change_{idx}'] = fc_value
                entry[f'fc_color_{idx}'] = self.get_color(fc_value) if fc_value is not None else [128, 128, 128]
            catalog[entry['hmdb_id'] or hmdb_key] = entry
        return catalog

def generate_pathway_json(pathway_id, data, settings, skip_disk_write=False, debug_write=False):
    try:
        debug_output_dir = _resolve_debug_output_dir(settings) if (debug_write and _debug_file_output_enabled()) else None
        # Debug logging: capture what enters m4
        if debug_output_dir is not None:
            try:
                debug_in_path = debug_output_dir / "loaded_pathway_m4.txt"
                with debug_in_path.open("w", encoding="utf-8") as fh:
                    json.dump(
                        {
                            "pathway_id": pathway_id,
                            "settings": settings,
                            "data_keys": list(data.keys()) if isinstance(data, dict) else None,
                        },
                        fh,
                        indent=2,
                        default=str,
                    )
            except Exception as log_exc:  # pragma: no cover - debug helper
                print(f"Debug log (m4 input) failed: {log_exc}")

        def _load_df(entry):
            if entry is None:
                return None
            if isinstance(entry, pd.DataFrame):
                return entry.copy()
            headers = entry.get('data_headers')
            rows = entry.get('data_rows')
            if headers is not None and rows is not None:
                try:
                    return pd.DataFrame(rows, columns=headers)
                except Exception:
                    pass
            file_path = entry.get('file_path')
            if file_path:
                return pd.read_csv(file_path, sep="\t")
            return None

        proteomic_data = _load_df(data['protein'])
        if proteomic_data is None:
            raise RuntimeError("No protein data provided.")

        ptm_datasets = data['ptm']
        metabolite_dataset = data.get('metabolite') or {}
        settings['protein_file_path'] = data['protein'].get('file_path', '')
        settings['main_columns'] = data['protein']['main_columns']
        settings['protein_tooltip_columns'] = data['protein'].get('tooltip_columns', ['Gene Symbol', 'Uniprot_ID'])
        settings['metabolite_file_path'] = metabolite_dataset.get('file_path', '')
        settings['metabolite_main_columns'] = metabolite_dataset.get('main_columns', [])
        settings['metabolite_tooltip_columns'] = metabolite_dataset.get('tooltip_columns', [])
        pathway_api = get_pathway_api(settings.get('pathway_source', 'kegg'))
        species_code = settings.get("species_code") or ""
        if species_code:
            setattr(pathway_api, "species_code", species_code)
        species_hint = settings.get("_species_full_name") or settings.get("species")
        pathway_file = pathway_api.download_pathway_data(pathway_id, species_hint=species_hint)
        print(f"Pathway file downloaded: {pathway_file}")
        entries, groups, arrows = pathway_api.parse_pathway(pathway_file)
        print(f"Parsed {len(entries)} entries, {len(groups)} groups, and {len(arrows)} arrows")
        # Debug: capture what we got from the pathway API (a1_factory output)
        if debug_output_dir is not None:
            try:
                debug_in_path = debug_output_dir / "loaded_pathway_m4.txt"
                with debug_in_path.open("w", encoding="utf-8") as fh:
                    json.dump(
                        {
                            "pathway_id": pathway_id,
                            "pathway_source": settings.get("pathway_source"),
                            "api": type(pathway_api).__name__,
                            "pathway_file": pathway_file,
                            "counts": {"entries": len(entries), "groups": len(groups), "arrows": len(arrows)},
                            "entries": entries,
                            "groups": groups,
                            "arrows": arrows,
                        },
                        fh,
                        indent=2,
                        default=str,
                    )
            except Exception as log_exc:  # pragma: no cover - debug helper
                print(f"Debug log (m4 input after parse) failed: {log_exc}")

        loaded_ptm = []
        for dataset in ptm_datasets:
            df = _load_df(dataset)
            if df is None:
                raise RuntimeError(f"No PTM data provided for dataset {dataset.get('type', 'PTM')}")
            ds = dict(dataset)
            ds['dataframe'] = df
            if 'file_path' not in ds:
                ds['file_path'] = dataset.get('file_path', '')
            loaded_ptm.append(ds)

        metabolite_df = _load_df(metabolite_dataset) if metabolite_dataset else None
        processor = PathwayProcessor(entries, proteomic_data, loaded_ptm, settings, metabolite_data=metabolite_df, metabolite_config=metabolite_dataset)
        print("PathwayProcessor created")
        json_data = processor.process_pathway(entries, groups, arrows, proteomic_data, loaded_ptm, skip_disk_write=skip_disk_write)
        metabolite_catalog = processor.build_full_metabolite_catalog()
        if metabolite_catalog:
            json_data["_global_metabolite_catalog"] = metabolite_catalog
        print("Pathway JSON generated")

        # Debug logging: capture what m4 returns (consumed by m5/m3)
        if debug_output_dir is not None:
            try:
                debug_out_path = debug_output_dir / "loaded_pathway_m3.txt"
                with debug_out_path.open("w", encoding="utf-8") as fh:
                    json.dump(json_data, fh, indent=2, default=str)
            except Exception as log_exc:  # pragma: no cover - debug helper
                print(f"Debug log (m4 output) failed: {log_exc}")

        return json_data
    except Exception as e:
        print(f"Error: {e}")
        return None
# Expose default settings and data so other scripts can import them and call generate_pathway_json
DEFAULT_SETTINGS = {
    'pathway_id': 'hsa04010',
    'pathway_source': 'kegg',
    'protein_selection_option': 2,
    'ptm_selection_option': 2,
    'ptm_max_display': 4,
    'show_background_image': True,
    'display_types': ['prot_box'],
    'show_groups': False,
    'show_multi_protein_indicator': True,
    'show_arrows': True,
    'show_text_boxes': True,
    'debug_mode': False,
    'negative_color': (179, 21, 41),
    'positive_color': (16, 101, 171),
    'max_negative': -2,
    'max_positive': 2,
    'gradient_stops': [
        {'value': -2.0, 'color': [179, 21, 41]},
        {'value': 0.0, 'color': [255, 255, 255]},
        {'value': 2.0, 'color': [16, 101, 171]},
    ],
    'prot_label_font': 'Arial',
    'prot_label_size': 9,
    'ptm_label_font': 'Arial',
    'ptm_label_color': (0, 0, 0),
    'ptm_label_size': 6,
    'ptm_circle_radius': 5,
    'ptm_circle_spacing': 4,
    'prot_outline_width': 1,
    'ptm_outline_width': 1,
    'use_black_protein_outlines': False,
    'temporal_mode': False,
    'protein_tooltip_columns': ['Gene Symbol', 'Uniprot_ID', 'ER+_Est-_x-y_TNBC'],
    'use_original_protbox_size': True
}

DEFAULT_DATA = {
    "protein": {
        "file_path": r"C:\\Users\\clayt\\OneDrive - Brigham Young University\\Desktop\\Graduate_Documents\\data_analysis\\2025.03.27 - BCp2_data\\Phosmap\\BCp2_ProtMaphsaanno.txt",
        "uniprot_column": "Uniprot_ID",
        "kegg_column": "KEGG_hsa",
        "gene_column": "Gene Symbol",
        "main_columns": ["ER+_Est-_x-y_ER+", "ER+_Est-_x-y_TNBC", "ER+_x-y_TNBC"],
        "tooltip_columns": ["Gene Symbol", "Uniprot_ID", "ER+_Est-_x-y_TNBC"]
    },
    "ptm": [
        {
            "type": "Phosphorylation",
            "file_path": r"C:\\Users\\clayt\\OneDrive - Brigham Young University\\Desktop\\Graduate_Documents\\data_analysis\\2025.03.27 - BCp2_data\\Phosmap\\BCp2_PhosMaphsanno.txt",
            "uniprot_column": "T: Uniprot_ID",
            "site_column": "T: Site Position",
            "shape": "Circle",
            "main_columns": [("ER+_Est-_x-y_ER+", "ER+_Est-_x-y_ER+"),
                             ("ER+_Est-_x-y_TNBC", "ER+_Est-_x-y_TNBC"),
                             ("ER+_x-y_TNBC", "ER+_x-y_TNBC")],
            "modulation_column": "C: Regulatory site",
            "tooltip_columns": ["C: Regulatory site function", "C: Regulatory site process"],
            "ptm_symbol_list": [
                {
                    "symbol_label_1_dict": {
                        "symbol": "↑",
                        "symbol_font": "Segoe UI Symbol",
                        "symbol_icon": "active.svg",
                        "symbol_size": 8,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -0.1,
                        "symbol_y_offset": -0.5,   # was -0.3
                        "statement_type": "if_or_statement",
                        "header_to_search": "C: Regulatory site function",
                        "search_text_1": "activity, induced",
                        "search_text_2": "enzymatic activity, induced",
                        "search_text_3": "",
                        "search_text_4": ""
                    }
                },
                {
                    "symbol_label_2_dict": {
                        "symbol": "x",
                        "symbol_font": "Arial",
                        "symbol_icon": "inhibit.svg",
                        "symbol_size": 8,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -0.3,    #-1.2
                        "symbol_y_offset": 0,   # was -1
                        "statement_type": "if_or_statement",
                        "header_to_search": "C: Regulatory site function",
                        "search_text_1": "activity, inhibited",
                        "search_text_2": "enzymatic activity, inhibited",
                        "search_text_3": "",
                        "search_text_4": ""
                    }
                },
                {
                    "symbol_label_3_dict": {
                        "symbol": "+",
                        "symbol_font": "Arial",
                        "symbol_icon": "regsite.svg",
                        "symbol_size": 9,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -0.1,
                        "symbol_y_offset": 0,
                        "statement_type": "if_or_and_or_statement",
                        "header_to_search": "C: Regulatory site function",
                        "search_text_1": "activity, inhibited",
                        "search_text_2": "enzymatic activity, inhibited",
                        "search_text_3": "activity, induced",
                        "search_text_4": "enzymatic activity, induced"
                    }
                },
                {
                    "symbol_label_4_dict": {
                        "symbol": "+",
                        "symbol_font": "Arial",
                        "symbol_icon": "regsite.svg",
                        "symbol_size": 9,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -0.1,
                        "symbol_y_offset": 0,
                        "statement_type": "if_and_notlabeled_statement",
                        "header_to_search": "PSP: regulatory_site",
                        "search_text_1": "+",
                        "search_text_2": "",
                        "search_text_3": "",
                        "search_text_4": ""
                    }
                }
            ]
        }
    ]
}


def get_default_json(data_override=None, settings_override=None, skip_disk_write=False, debug_write=False):
    """Return the JSON pathway data using the DEFAULT_DATA and DEFAULT_SETTINGS.

    Optional overrides may be supplied to change the DEFAULT_DATA or DEFAULT_SETTINGS
    (for example to point at different input file paths) without editing this file.

    Args:
        data_override (dict|None): partial dict to merge on top of DEFAULT_DATA
        settings_override (dict|None): partial dict to merge on top of DEFAULT_SETTINGS
        skip_disk_write (bool): if True, avoid writing intermediate/output files to disk
        debug_write (bool): if True, write debug snapshot files for m4/m3 handoff

    Returns:
        dict: generated pathway JSON (or a minimal fallback structure on error)
    """
    # Merge settings
    settings = dict(DEFAULT_SETTINGS)
    if settings_override:
        settings.update(settings_override)

    # Merge data (shallow merge at top-level is sufficient for file path overrides)
    import copy
    data = copy.deepcopy(DEFAULT_DATA)
    if data_override:
        for k, v in data_override.items():
            data[k] = v

    generation_error = ""
    try:
        json_data = generate_pathway_json(
            settings.get('pathway_id', 'hsa04010'),
            data,
            settings,
            skip_disk_write=skip_disk_write,
            debug_write=debug_write,
        )
        if json_data is None:
            raise RuntimeError('generate_pathway_json returned None')
        return json_data
    except FileNotFoundError as e:
        generation_error = str(e)
        print(f"Warning: data file not found while generating default JSON: {e}")
    except Exception as e:
        generation_error = str(e)
        print(f"Warning: could not generate default JSON ({e}). Falling back to minimal structure.")

    # Fallback minimal structure so consumers can still render/debug the UI
    return {
        'general_data': {'settings': settings},
        'protein_data': {},
        'protbox_data': [],
        'groups': [],
        'arrows': [],
        'compound_data': [],
        'text_data': [],
        '_generation_error': generation_error,
    }


if __name__ == "__main__":
    # preserve original behavior when running m4 directly
    json_data = get_default_json()
    if json_data:
        print("Generated JSON data successfully (not printed).")
