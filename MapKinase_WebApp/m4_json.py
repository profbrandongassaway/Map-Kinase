import sys
import pandas as pd
import numpy as np
import os
from datetime import datetime
from collections import defaultdict
import json
import html
from concurrent.futures import ThreadPoolExecutor, as_completed
from factory import get_pathway_api

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

def parse_exceptions_file(file_path='exceptions_file.txt'):
    exceptions = {
        'global': {'proximity': [], 'specific': {}},
        'hsa_specific': {}
    }
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        base_dir = os.path.join(sys._MEIPASS, 'Scripts')
    else:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(base_dir, file_path)
    try:
        if not os.path.exists(full_path):
            print(f"Exceptions file not found at {full_path}")
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
        output_dir = os.path.join(phosmap_dir, directory)
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
    print(f"Classifying {len(phospho_data)} PTM sites with symbol list: {ptm_symbol_list}")
    print(f"PTM data columns: {phospho_data.columns.tolist()}")

    valid_symbol_list = []
    for symbol_dict in ptm_symbol_list:
        if not symbol_dict or len(symbol_dict) != 1:
            print(f"Warning: Skipping invalid rule with incorrect structure: {symbol_dict}")
            continue
        inner_dict_key = list(symbol_dict.keys())[0]
        inner_dict = symbol_dict[inner_dict_key]
        header = inner_dict.get('header_to_search', '')
        if not header:
            print(f"Warning: Skipping rule with empty header_to_search: {symbol_dict}")
            continue
        if header not in phospho_data.columns:
            print(f"Warning: Header {header} not in PTM data columns: {phospho_data.columns}")
            continue
        valid_symbol_list.append(symbol_dict)

    if not valid_symbol_list:
        print("Error: No valid rules in ptm_symbol_list. All PTMs will be classified as 'none'.")
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
    print(f"Classification complete. Value counts:\n{phospho_data[output_col].value_counts()}")
    return phospho_data

class PathwayProcessor:
    def __init__(self, entries, proteomic_data, ptm_datasets, settings):
        self.proteomic_data = proteomic_data
        self.ptm_datasets = {f"ptm_{i}": dataset for i, dataset in enumerate(ptm_datasets)}
        self.protein_data_map = {}  # Initialize protein_data_map
        self.current_proteins = {}  # Initialize current_proteins
        for dataset_id, dataset in self.ptm_datasets.items():
            dataset['data'] = pd.read_csv(dataset['file_path'], sep="\t")
            print(f"PTM dataset {dataset_id}: {len(dataset['data'])} rows, columns: {dataset['data'].columns.tolist()}")
            symbol_list = dataset.get('ptm_symbol_list', [])
            reg_site_col = dataset.get('modulation_column', 'C: Regulatory site')
            reg_function_col = dataset.get('tooltip_columns', ['C: Regulatory site function'])[0]
            print(f"Calling classify_phosphosite_function with:")
            print(f"  ptm_symbol_list: {symbol_list}")
            print(f"  reg_site_col: {reg_site_col}")
            print(f"  reg_function_col: {reg_function_col}")
            if not symbol_list:
                print(f"Warning: ptm_symbol_list is empty for dataset {dataset_id}. PTMs will be classified as 'none'.")
            if reg_site_col not in dataset['data'].columns:
                print(f"Warning: reg_site_col {reg_site_col} not in PTM data columns.")
            if reg_function_col not in dataset['data'].columns:
                print(f"Warning: reg_function_col {reg_function_col} not in PTM data columns.")
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
        self.hsa_id_column = settings.get('hsa_id_column', 'KEGG_hsa')
        self.prot_uniprot_column = settings.get('prot_uniprot_column', 'Uniprot_ID')
        self.gene_name_column = settings.get('gene_name_column', 'Gene Symbol')
        self.protein_tooltip_columns = settings.get('protein_tooltip_columns', ['Gene Symbol', 'Uniprot_ID'])
        self.negative_color = settings.get('negative_color', (255, 0, 0))
        self.positive_color = settings.get('positive_color', (0, 0, 255))
        self.max_negative = settings.get('max_negative', -2)
        self.max_positive = settings.get('max_positive', 2)
        self.ptm_label_color = settings.get('ptm_label_color', (0, 0, 0))
        self.ptm_circle_radius = settings.get('ptm_circle_radius', 5)
        self.ptm_circle_spacing = settings.get('ptm_circle_spacing', 4)
        self.display_types = settings.get('display_types', ['prot_box'])
        default_workers = os.cpu_count() or 4
        self._ptm_worker_count = max(1, min(settings.get('ptm_worker_count', default_workers), 16))
        self._ptm_summary_cache = {}

    def get_color(self, fold_change):
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
        valid_proteins = [p for p in proteins if p in self.proteomic_data[self.hsa_id_column].values]
        if not valid_proteins:
            return None
        try:
            if self.protein_selection_option == 1:
                selected_protein = valid_proteins[0]
            elif self.protein_selection_option == 2:
                max_fold_change = -float('inf')
                selected_protein = None
                for protein in valid_proteins:
                    for fc_col in self.fold_change_columns:
                        fold_change = self.proteomic_data.loc[
                            self.proteomic_data[self.hsa_id_column] == protein, fc_col].values
                        if len(fold_change) > 0 and abs(fold_change[0]) > max_fold_change:
                            max_fold_change = abs(fold_change[0])
                            selected_protein = protein
                if not selected_protein:
                    selected_protein = valid_proteins[0]
            else:
                selected_protein = valid_proteins[0]
            uniprot_id = self.proteomic_data.loc[
                self.proteomic_data[self.hsa_id_column] == selected_protein, self.prot_uniprot_column].values
            uniprot_id = uniprot_id[0] if len(uniprot_id) > 0 else None
            annotations = []
            for col in self.protein_tooltip_columns:
                value = self.proteomic_data.loc[
                    self.proteomic_data[self.hsa_id_column] == selected_protein, col].values
                annotations.append(str(value[0]) if len(value) > 0 else '')
            return {'protein': selected_protein, 'uniprot_id': uniprot_id, 'annotations': annotations}
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
        try:
            if self.ptm_selection_option == 1:
                print("Applying PTM selection option 1: Highest absolute fold change")
                combined_ptms = combined_ptms.sort_values(
                    by=[col[1] for dataset in self.ptm_datasets.values() for col in dataset['main_columns']],
                    key=lambda x: x.abs(),
                    ascending=False
                )
            elif self.ptm_selection_option == 2:
                print("Applying PTM selection option 2: Prefer modulating phospho PTMs")
                def sort_key(row):
                    try:
                        dataset_id = row['dataset_id']
                        dataset = self.ptm_datasets[dataset_id]
                        fold_change = row[dataset['main_columns'][0][1]]
                        is_modulating = row['is_modulating']
                        is_phospho = row['is_phospho']
                        return (-is_phospho, -is_modulating, -abs(fold_change) if not pd.isna(fold_change) else 0)
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
                    by=[col[1] for dataset in self.ptm_datasets.values() for col in dataset['main_columns']],
                    key=lambda x: x.abs(),
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
                        'negative_color': self.settings.get('negative_color', (255, 0, 0)),
                        'positive_color': self.settings.get('positive_color', (0, 0, 255)),
                        'max_negative': self.settings.get('max_negative', -2),
                        'max_positive': self.settings.get('max_positive', 2),
                        'prot_label_font': self.settings.get('prot_label_font', 'Arial'),
                        'prot_label_size': self.settings.get('prot_label_size', 12),
                        'ptm_label_font': self.settings.get('ptm_label_font', 'Arial'),
                        'ptm_label_color': self.settings.get('ptm_label_color', (0, 0, 0)),
                        'ptm_label_size': self.settings.get('ptm_label_size', 10),
                        'ptm_circle_radius': self.settings.get('ptm_circle_radius', 5),
                        'ptm_circle_spacing': self.settings.get('ptm_circle_spacing', 4),
                        'protein_tooltip_columns': self.settings.get('protein_tooltip_columns', ['Gene Symbol', 'Uniprot_ID'])
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
                                'file_path': dataset['file_path'],
                                'uniprot_column': dataset['uniprot_column'],
                                'site_column': dataset['site_column'],
                                'shape': dataset['shape'],
                                'main_columns': dataset['main_columns'],
                                'modulation_column': dataset.get('modulation_column', ''),
                                'tooltip_columns': dataset.get('tooltip_columns', []),
                                'ptm_symbol_list': dataset.get('ptm_symbol_list', [])
                            } for dataset in ptm_datasets
                        ]
                    }
                },
                'protein_data': {},
                'protbox_data': [],
                'groups': [],
                'arrows': [],
                'compound_data' : [],
                'text_data' : []
            }
            for group in groups:
                group_entry = {
                    'group_id': group['id'],
                    'protbox_ids': group.get('protbox_ids', [])
                }
                json_data['groups'].append(group_entry)
            prot_entries = [e for e in entries if e.get('type') == 'prot_box']
            opposite_side = {'East': 'West', 'West': 'East', 'North': 'South', 'South': 'North'}
            for arrow in arrows:
                entry1 = find_entry_or_group(arrow['entry1'], entries, groups)
                entry2 = find_entry_or_group(arrow['entry2'], entries, groups)
                if not (entry1 and entry2) or entry1.get('type') != 'prot_box' or entry2.get('type') != 'prot_box':
                    continue
                center1_x = entry1['x'] + entry1['width'] / 2
                center1_y = entry1['y'] + entry1['height'] / 2
                center2_x = entry2['x'] + entry2['width'] / 2
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
                    'type': arrow['type']
                }
                json_data['arrows'].append(arrow_entry)
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
                    'first_name': entry.get("first_name", entry["name"].split(",")[0].strip())
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
                    if protein not in self.proteomic_data[self.hsa_id_column].values:
                        continue
                    mask = self.proteomic_data[self.hsa_id_column] == protein
                    uniprot_vals = self.proteomic_data.loc[mask, self.prot_uniprot_column].values
                    uniprot_id = uniprot_vals[0] if len(uniprot_vals) > 0 else None
                    if not uniprot_id:
                        continue
                    annotations = []
                    for col in self.protein_tooltip_columns:
                        value = self.proteomic_data.loc[mask, col].values
                        annotations.append(str(value[0]) if len(value) > 0 else '')
                    valid_protein_dicts.append({'protein': protein, 'uniprot_id': uniprot_id, 'annotations': annotations})
                    uniprot_ids.append(uniprot_id)
                    all_uniprot_ids.add(uniprot_id)
                protbox_protein_cache[entry["id"]] = {
                    'proteins': proteins,
                    'valid_protein_dicts': valid_protein_dicts,
                    'uniprot_ids': uniprot_ids
                }
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
                protbox_entry = {
                    'protbox_id': entry["id"],
                    'proteins': uniprot_ids,
                    'backup_label': self.protein_data_map[entry["id"]]['first_name'],
                    'x': entry["x"],
                    'y': entry["y"],
                    'width': 46,
                    'height': 17
                }
                if protbox_ptm_overrides:
                    protbox_entry['ptm_overrides'] = protbox_ptm_overrides
                json_data['protbox_data'].append(protbox_entry)
            for entry in entries:
                if entry.get('type') == 'compound':
                    compound = {
                        'compound_id': entry['id'],  # KGML numeric id
                        'kegg_compound': entry['name'],  # e.g., 'cpd:C00165'
                        'label': entry.get('first_name', ''),  # usually the C-number in MAPK
                        'x': entry['x'],
                        'y': entry['y'],
                        'width': entry['width'],
                        'height': entry['height'],
                        'fgcolor': entry.get('fgcolor', '#000000'),
                        'bgcolor': entry.get('bgcolor', '#FFFFFF'),
                        'graphics_type': entry.get('graphics_type', ''),  # likely 'circle'
                        'link': entry.get('link', '')
                    }
                    json_data['compound_data'].append(compound)
            for entry in entries:
                if entry.get('type') == 'map' or entry.get('graphics_type') == 'roundrectangle':
                    label = entry.get('first_name', '')
                    # Optionally strip 'TITLE:' prefix if present
                    if label.upper().startswith('TITLE:'):
                        label = label.split(':', 1)[1].strip()

                    text_item = {
                        'text_id': entry['id'],
                        'label': label,
                        'x': entry['x'],
                        'y': entry['y'],
                        'width': entry['width'],
                        'height': entry['height'],
                        'fgcolor': entry.get('fgcolor', '#000000'),
                        'bgcolor': entry.get('bgcolor', '#FFFFFF'),
                        'graphics_type': entry.get('graphics_type', ''),  # often 'roundrectangle'
                        'link': entry.get('link', '')
                    }
                    json_data['text_data'].append(text_item)

            if not skip_disk_write:
                temp_file = savefile(self.settings.get('pathway_id', 'hsa04010'),
                                     self.settings.get('output_subdir', 'output/testing_file_001'),
                                     '_visprots', extension='.txt', include_timestamp=False)
                try:
                    with open(temp_file, 'w') as f:
                        for protein in chosen_proteins:
                            f.write(f"{protein}\n")
                except (IOError, PermissionError) as e:
                    print(f"Error writing to temporary file {temp_file}: {e}")
                output_file = savefile(self.settings.get('pathway_id', 'hsa04010'),
                                       self.settings.get('output_subdir', 'output/testing_file_001'),
                                       '_pathway_data', extension='.json')
                with open(output_file, 'w') as f:
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
            x = entry_data['x']
            y = entry_data['y']
            width = entry_data['width']
            height = entry_data['height']
            proteins = entry_data['proteins']
            protein = protein_dict['protein']
            uniprot_id = protein_dict['uniprot_id']
            annotations = protein_dict.get('annotations', [])
            protein_in_data = protein in self.proteomic_data[self.hsa_id_column].values
            valid_proteins = [p for p in proteins if p in self.proteomic_data[self.hsa_id_column].values]
            protein_entry = {
                'label': '',
                'label_color': [0, 0, 0],
                'transcriptomic_color': [],
                'annotations': ','.join(f'"{ann}"' for ann in annotations if ann),
                'PTMs': {}
            }
            for idx, fc_col in enumerate(self.fold_change_columns, 1):
                protein_entry[f'fc_color_{idx}'] = [128, 128, 128]
                protein_entry[f'fold_change_{idx}'] = None
            if protein_in_data:
                gene_name = self.proteomic_data.loc[
                    self.proteomic_data[self.hsa_id_column] == protein, self.gene_name_column].values[0]
                protein_entry['label'] = gene_name
                for idx, fc_col in enumerate(self.fold_change_columns, 1):
                    fold_change = self.proteomic_data.loc[
                        self.proteomic_data[self.hsa_id_column] == protein, fc_col].values
                    if len(fold_change) > 0:
                        raw_value = fold_change[0]
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
            else:
                gene_name = entry_data['first_name'].split(",")[0].strip()
                protein_entry['label'] = gene_name
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
                            tooltip_plain_lines.append(f"{label}: {value}")
                            safe_label = html.escape(label)
                            safe_value = html.escape(value)
                            tooltip_html_lines.append(f"<strong>{safe_label}:</strong> {safe_value}")
                        ptm_entry['tooltip'] = '\n'.join(tooltip_plain_lines)
                        ptm_entry['tooltip_html'] = '<br/>'.join(tooltip_html_lines)
                        for idx, (_, fc_col) in enumerate(dataset['main_columns'], 1):
                            ptm_entry.setdefault(f'fc_color_{idx}', [128, 128, 128])
                            ptm_entry[f'fold_change_{idx}'] = None
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
                                            symbol_x = shape_x + symbol_dict.get('symbol_x_offset', 0)
                                            symbol_y = shape_y + symbol_dict.get('symbol_y_offset', 0)
                                            symbol_font = symbol_dict.get('symbol_font', 'Arial')
                                            symbol_size = symbol_dict.get('symbol_size', 6)
                                            print(f"PTM {ptm_key}: symbol={symbol}, symbol_x={symbol_x}, symbol_y={symbol_y}, pos_key={pos_key}")
                                        else:
                                            print(f"Warning: No symbol dict for {symbol_type} in dataset {dataset_id}")
                                    else:
                                        print(f"Warning: ptm_type_idx {ptm_type_idx} out of range for ptm_datasets")
                                except (IndexError, ValueError) as e:
                                    print(f"Error processing symbol for PTM {ptm_key}: {e}")
                            ptm_entry['symbol'] = symbol
                            ptm_entry['symbol_x'] = float(symbol_x)
                            ptm_entry['symbol_y'] = float(symbol_y)
                            ptm_entry['symbol_color'] = symbol_color
                            ptm_entry['symbol_font'] = symbol_font
                            ptm_entry['symbol_size'] = symbol_size
                            print(f"PTM {ptm_key}: shape=({shape_x}, {shape_y}), label=({label_x}, {label_y}), centering={label_centering}")
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
                tooltip_plain_lines.append(f"{label}: {value}")
                safe_label = html.escape(label)
                safe_value = html.escape(value)
                tooltip_html_lines.append(f"<strong>{safe_label}:</strong> {safe_value}")
            ptm_entry['tooltip'] = '\n'.join(tooltip_plain_lines)
            ptm_entry['tooltip_html'] = '<br/>'.join(tooltip_html_lines)
            for idx, (_, fc_col) in enumerate(dataset.get('main_columns', []), 1):
                ptm_entry.setdefault(f'fc_color_{idx}', [128, 128, 128])
                ptm_entry[f'fold_change_{idx}'] = None
                fc_value = self._safe_float(row.get(fc_col))
                ptm_entry[f'fold_change_{idx}'] = fc_value
                if fc_value is not None:
                    ptm_entry[f'fc_color_{idx}'] = self.get_color(fc_value)
            symbol_type = ptm_entry['symbol_type']
            if symbol_type and str(symbol_type).lower() != 'none':
                symbol_dict = next(
                    (item.get(f"{symbol_type}_dict") for item in dataset.get('ptm_symbol_list', []) if f"{symbol_type}_dict" in item),
                    None
                )
                if symbol_dict:
                    ptm_entry['symbol'] = symbol_dict.get('symbol', '')
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
            for idx, fc_col in enumerate(self.fold_change_columns, 1):
                protein_entry[f'fc_color_{idx}'] = [128, 128, 128]
                protein_entry[f'fold_change_{idx}'] = None
                fc_value = self._safe_float(row.get(fc_col))
                protein_entry[f'fold_change_{idx}'] = fc_value
                if fc_value is not None:
                    protein_entry[f'fc_color_{idx}'] = self.get_color(fc_value)
            catalog[uniprot_id] = protein_entry
        return catalog

def generate_pathway_json(pathway_id, data, settings, skip_disk_write=False):
    try:
        proteomic_data = pd.read_csv(data['protein']['file_path'], sep="\t")
        ptm_datasets = data['ptm']
        settings['protein_file_path'] = data['protein']['file_path']
        settings['main_columns'] = data['protein']['main_columns']
        settings['protein_tooltip_columns'] = data['protein'].get('tooltip_columns', ['Gene Symbol', 'Uniprot_ID'])
        pathway_api = get_pathway_api(settings.get('pathway_source', 'kegg'))
        pathway_file = pathway_api.download_pathway_data(pathway_id)
        print(f"Pathway file downloaded: {pathway_file}")
        entries, groups, arrows = pathway_api.parse_pathway(pathway_file)
        print(f"Parsed {len(entries)} entries, {len(groups)} groups, and {len(arrows)} arrows")
        processor = PathwayProcessor(entries, proteomic_data, ptm_datasets, settings)
        print("PathwayProcessor created")
        json_data = processor.process_pathway(entries, groups, arrows, proteomic_data, ptm_datasets, skip_disk_write=skip_disk_write)
        print("Pathway JSON generated")
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
    'negative_color': (179, 21, 41),
    'positive_color': (16, 101, 171),
    'max_negative': -2,
    'max_positive': 2,
    'prot_label_font': 'Arial',
    'prot_label_size': 9,
    'ptm_label_font': 'Arial',
    'ptm_label_color': (0, 0, 0),
    'ptm_label_size': 6,
    'ptm_circle_radius': 5,
    'ptm_circle_spacing': 4,
    'protein_tooltip_columns': ['Gene Symbol', 'Uniprot_ID', 'ER+_Est-_x-y_TNBC']
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
                        "symbol_size": 8,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -0.1,
                        "symbol_y_offset": -0.3,   # was -1
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
                        "symbol_size": 8,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -1.2,
                        "symbol_y_offset": -0.6,   # was -1
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
                        "symbol_size": 9,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -0.1,
                        "symbol_y_offset": -0.8,
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
                        "symbol_size": 9,
                        "symbol_color": (0, 0, 0),
                        "symbol_x_offset": -0.1,
                        "symbol_y_offset": -0.8,
                        "statement_type": "if_and_notlabeled_statement",
                        "header_to_search": "C: Regulatory site",
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


def get_default_json(data_override=None, settings_override=None, skip_disk_write=False):
    """Return the JSON pathway data using the DEFAULT_DATA and DEFAULT_SETTINGS.

    Optional overrides may be supplied to change the DEFAULT_DATA or DEFAULT_SETTINGS
    (for example to point at different input file paths) without editing this file.

    Args:
        data_override (dict|None): partial dict to merge on top of DEFAULT_DATA
        settings_override (dict|None): partial dict to merge on top of DEFAULT_SETTINGS
        skip_disk_write (bool): if True, avoid writing intermediate/output files to disk

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

    try:
        json_data = generate_pathway_json(settings.get('pathway_id', 'hsa04010'), data, settings, skip_disk_write=skip_disk_write)
        if json_data is None:
            raise RuntimeError('generate_pathway_json returned None')
        return json_data
    except FileNotFoundError as e:
        print(f"Warning: data file not found while generating default JSON: {e}")
    except Exception as e:
        print(f"Warning: could not generate default JSON ({e}). Falling back to minimal structure.")

    # Fallback minimal structure so consumers can still render/debug the UI
    return {
        'general_data': {'settings': settings},
        'protein_data': {},
        'protbox_data': [],
        'groups': [],
        'arrows': [],
        'compound_data': [],
        'text_data': []
    }


if __name__ == "__main__":
    # preserve original behavior when running m4 directly
    json_data = get_default_json()
    if json_data:
        print("Generated JSON data successfully (not printed).")
