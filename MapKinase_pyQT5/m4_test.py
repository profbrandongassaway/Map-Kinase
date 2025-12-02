import sys
import requests
import xml.etree.ElementTree as ET
from PIL import Image
import io
import pandas as pd
import numpy as np
import os
from datetime import datetime
from collections import defaultdict
from PyQt5.QtWidgets import (
    QApplication, QGraphicsView, QGraphicsScene, QGraphicsRectItem,
    QGraphicsEllipseItem, QGraphicsTextItem, QGraphicsItemGroup,
    QGraphicsPixmapItem, QGraphicsPolygonItem, QMenu, QAction, QGraphicsItem,
    QPushButton, QFileDialog, QMessageBox
)
from PyQt5.QtGui import (
    QPixmap, QColor, QFont, QPainter, QImage, QPolygonF, QPen,
    QPdfWriter, QPagedPaintDevice, QPageSize, QBrush
)
from PyQt5.QtCore import Qt, QRectF, QPointF, QSize, QSizeF, QMarginsF, QTimer

from e1_popup_errorwindow import show_error_popup

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
            show_error_popup(f"Exceptions file not found at {full_path}", "File Error")
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
        show_error_popup(f"Error reading exceptions file {full_path}: {str(e)}", "File Error")
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

def savefile(input_file, directory, suffix, extension=".xml", include_timestamp=True):
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
        show_error_popup(f"Error creating output file path for {input_file}: {e}", "File Error")
        return None

def download_kegg_kgml(pathway_id):
    try:
        kgml_url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
        response = requests.get(kgml_url, timeout=10)
        response.raise_for_status()
        output_file = savefile(pathway_id, 'kegg_pathway', '', extension=".xml", include_timestamp=False)
        with open(output_file, "w") as f:
            f.write(response.text)
        return output_file
    except requests.exceptions.RequestException as e:
        show_error_popup(f"Failed to download KGML file for {pathway_id}: {e}", "Error")
        return None
    except IOError as e:
        show_error_popup(f"Failed to save KGML file: {e}", "Error")
        return None

def download_kegg_image(pathway_id):
    try:
        img_url = f"https://rest.kegg.jp/get/{pathway_id}/image"
        response = requests.get(img_url, timeout=10)
        response.raise_for_status()
        return Image.open(io.BytesIO(response.content))
    except requests.exceptions.RequestException as e:
        show_error_popup(f"Failed to download pathway image for {pathway_id}: {e}", "File Error")
        return None
    except Image.UnidentifiedImageError as e:
        show_error_popup(f"Invalid image data for {pathway_id}: {e}", "File Error")
        return None

def parse_kgml(kgml_file):
    entries = []
    groups = []
    try:
        tree = ET.parse(kgml_file)
        root = tree.getroot()
        for entry in root.findall("entry"):
            try:
                entry_id = entry.get("id")
                entry_type = entry.get("type")
                name = entry.get("name", "")
                graphics = entry.find("graphics")
                if graphics is None:
                    continue
                x = int(graphics.get("x", "0"))
                y = int(graphics.get("y", "0"))
                width = int(graphics.get("width", "0"))
                height = int(graphics.get("height", "0"))
                graphics_name = graphics.get("name", "")
                first_name = graphics_name.split(",")[0].strip() if graphics_name else name
                if entry_type == "group":
                    fgcolor = graphics.get("fgcolor", "#000000")
                    bgcolor = graphics.get("bgcolor", "#FFFFFF")
                    components = [comp.get("id") for comp in entry.findall("component")]
                    groups.append({
                        "id": entry_id,
                        "x": x,
                        "y": y,
                        "width": width,
                        "height": height,
                        "fgcolor": fgcolor,
                        "bgcolor": bgcolor,
                        "components": components
                    })
                else:
                    entries.append({
                        "id": entry_id,
                        "name": name,
                        "x": x,
                        "y": y,
                        "width": width,
                        "height": height,
                        "first_name": first_name,
                        "type": entry_type
                    })
            except (AttributeError, ValueError) as e:
                print(f"Error processing entry in KGML: {e}")
                continue
    except ET.ParseError as e:
        show_error_popup(f"Error parsing KGML file {kgml_file}: {e}")
        return [], []
    return entries, groups

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
    """
    Classifies PTM sites based on the ptm_symbol_list dictionary structure.
    Assigns a ptm_label_type to each PTM based on conditional statements.
    """
    phospho_data[output_col] = 'none'  # Default to 'none' for no label
    print(f"Classifying {len(phospho_data)} PTM sites with symbol list: {ptm_symbol_list}")
    print(f"PTM data columns: {phospho_data.columns.tolist()}")
    print(f"Sample data (first 5 rows):\n{phospho_data[[reg_site_col, reg_function_col]].head().to_string()}")

    # Validate ptm_symbol_list
    valid_symbol_list = []
    for symbol_dict in ptm_symbol_list:
        # Access the inner dictionary (e.g., symbol_label_1_dict)
        if not symbol_dict or len(symbol_dict) != 1:
            print(f"Warning: Skipping invalid rule with incorrect structure: {symbol_dict}")
            continue
        inner_dict_key = list(symbol_dict.keys())[0]  # e.g., 'symbol_label_1_dict'
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
        # First pass: Apply all rules except 'if_and_notlabeled_statement'
        for symbol_dict in [d for d in valid_symbol_list if
                            d[list(d.keys())[0]].get('statement_type') != 'if_and_notlabeled_statement']:
            inner_dict = symbol_dict[list(symbol_dict.keys())[0]]
            statement_type = inner_dict.get('statement_type', 'NA')
            header = inner_dict.get('header_to_search', '')
            search_text_1 = inner_dict.get('search_text_1', '')
            search_text_2 = inner_dict.get('search_text_2', '')
            search_text_3 = inner_dict.get('search_text_3', '')
            search_text_4 = inner_dict.get('search_text_4', '')
            label_key = list(symbol_dict.keys())[0]  # e.g., 'symbol_label_1_dict'
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

        # Second pass: Apply 'if_and_notlabeled_statement' rules
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

class ArrowItem(QGraphicsPolygonItem):
    def __init__(self, start_point, end_point):
        super().__init__()
        self.start_point = start_point
        self.end_point = end_point
        self.setPen(QPen(Qt.black, 1))
        self.setBrush(Qt.black)
        self.update_arrow()

    def update_arrow(self):
        line = QLineF(self.start_point, self.end_point)
        arrow_size = 10
        angle = np.arctan2(-line.dy(), line.dx())
        arrow_p1 = line.p2() - QPointF(np.sin(angle + np.pi / 3) * arrow_size,
                                       np.cos(angle + np.pi / 3) * arrow_size)
        arrow_p2 = line.p2() - QPointF(np.sin(angle + 2 * np.pi / 3) * arrow_size,
                                       np.cos(angle + 2 * np.pi / 3) * arrow_size)
        polygon = QPolygonF([line.p2(), arrow_p1, arrow_p2])
        self.setPolygon(polygon)

class PTMItem(QGraphicsItem):
    def __init__(self, x, y, width, height, ptm_row, dataset_id, viewer, shape='Circle', parent=None):
        super().__init__(parent)
        self.ptm_row = ptm_row
        self.dataset_id = dataset_id
        self.viewer = viewer
        self.shape = shape
        self.rect = QRectF(x, y, width, height)
        self.setAcceptHoverEvents(True)
        self.tooltip = None
        self.bg_rect = None
        self.hover_timer = QTimer()
        self.hover_timer.setSingleShot(True)
        self.hover_timer.timeout.connect(self.show_tooltip)

    def boundingRect(self):
        return self.rect

    def paint(self, painter, option, widget=None):
        painter.setRenderHint(QPainter.Antialiasing)
        color = self.viewer.get_color(self.ptm_row[self.viewer.ptm_datasets[self.dataset_id]['fold_change_column']])
        painter.setBrush(color)
        painter.setPen(QColor(0, 0, 0))
        if self.shape == 'Circle':
            painter.drawEllipse(self.rect)
        elif self.shape == 'Diamond':
            polygon = QPolygonF([
                QPointF(self.rect.center().x(), self.rect.top()),
                QPointF(self.rect.right(), self.rect.center().y()),
                QPointF(self.rect.center().x(), self.rect.bottom()),
                QPointF(self.rect.left(), self.rect.center().y())
            ])
            painter.drawPolygon(polygon)
        elif self.shape == 'Triangle':
            polygon = QPolygonF([
                QPointF(self.rect.center().x(), self.rect.top()),
                QPointF(self.rect.right(), self.rect.bottom()),
                QPointF(self.rect.left(), self.rect.bottom())
            ])
            painter.drawPolygon(polygon)

    def show_tooltip(self):
        try:
            dataset = self.viewer.ptm_datasets[self.dataset_id]
            tooltip_columns = dataset.get('tooltip_columns', [])
            if not tooltip_columns:
                print(f"No tooltip columns defined for dataset {self.dataset_id}")
                return
            annotations = []
            all_empty = True
            for col in tooltip_columns:
                if col in self.ptm_row.index:
                    value = self.ptm_row[col]
                    value_str = str(value) if not pd.isna(value) else ''
                    value_str = value_str.strip()
                    if value_str:
                        all_empty = False
                        annotations.append(f"<b>{col}</b>: {value_str}")
                    else:
                        annotations.append(f"<b>{col}</b>: N/A")
                else:
                    #print(f"Warning: Column {col} not found in ptm_row")
                    annotations.append(f"<b>{col}</b>: N/A")
            if all_empty:
                #print(f"All tooltip columns empty for PTM in dataset {self.dataset_id}, skipping tooltip")
                return
            if self.tooltip is not None:
                #print(f"Tooltip already exists for dataset {self.dataset_id}, skipping")
                return
            tooltip_text = "<br>".join(annotations)
            self.tooltip = QGraphicsTextItem()
            self.tooltip.setHtml(tooltip_text)
            self.tooltip.setFont(QFont("Arial", 8))
            self.tooltip.setDefaultTextColor(QColor(0, 0, 0))
            max_width = 300
            self.tooltip.setTextWidth(max_width)
            tooltip_pos = self.mapToScene(self.rect.center() + QPointF(self.rect.width() + 5, 0))
            self.tooltip.setPos(tooltip_pos)
            bounding_rect = self.tooltip.boundingRect()
            self.bg_rect = QGraphicsRectItem(bounding_rect.adjusted(-2, -2, 2, 2))
            self.bg_rect.setBrush(QColor(255, 255, 200, 200))
            self.bg_rect.setPen(QColor(0, 0, 0))
            self.bg_rect.setPos(tooltip_pos)
            self.bg_rect.setZValue(10)
            if self.scene():
                self.scene().addItem(self.bg_rect)
                self.scene().addItem(self.tooltip)
                self.tooltip.setZValue(11)
        except Exception as e:
            print(f"Error in PTM show_tooltip: {e}")
            self.cleanup_tooltip()

    def hoverEnterEvent(self, event):
        self.hover_timer.start(100)

    def hoverLeaveEvent(self, event):
        try:
            self.hover_timer.stop()
            self.cleanup_tooltip()
        except Exception as e:
            print(f"Error in PTM hoverLeaveEvent: {e}")

    def cleanup_tooltip(self):
        try:
            if self.tooltip is not None and self.scene():
                if self.bg_rect and self.bg_rect.scene():
                    self.scene().removeItem(self.bg_rect)
                if self.tooltip.scene():
                    self.scene().removeItem(self.tooltip)
            self.tooltip = None
            self.bg_rect = None
        except Exception as e:
            print(f"Error in PTM cleanup_tooltip: {e}")

class HoverRectItem(QGraphicsRectItem):
    def __init__(self, x, y, width, height, entry_id, viewer, parent=None):
        super().__init__(x, y, width, height, parent)
        self.entry_id = entry_id
        self.viewer = viewer
        self.setAcceptHoverEvents(True)
        self.setZValue(5)  # Higher than protein box
        self.setBrush(QBrush(Qt.NoBrush))  # Invisible
        self.setPen(QPen(Qt.NoPen))  # No border
        self.tooltip = None
        self.bg_rect = None
        self.hover_timer = QTimer()
        self.hover_timer.setSingleShot(True)
        self.hover_timer.timeout.connect(self.show_tooltip)

    def show_tooltip(self):
        try:
            if self.tooltip is not None:
                #print(f"Tooltip already exists for protein box {self.entry_id}")
                return
            protein = self.viewer.current_proteins[self.entry_id]['protein']
            if protein not in self.viewer.proteomic_data[self.viewer.hsa_id_column].values:
                #print(f"Protein {protein} not in proteomic data, skipping tooltip")
                return
            tooltip_columns = self.viewer.settings.get('protein_tooltip_columns', [])
            if not tooltip_columns:
                #print(f"No tooltip columns defined for protein {protein}")
                return
            protein_row = self.viewer.proteomic_data[
                self.viewer.proteomic_data[self.viewer.hsa_id_column] == protein
            ].iloc[0]
            annotations = []
            all_empty = True
            for col in tooltip_columns:
                if col in protein_row.index:
                    value = protein_row[col]
                    value_str = str(value) if not pd.isna(value) else ''
                    value_str = value_str.strip()
                    if value_str:
                        all_empty = False
                        annotations.append(f"<b>{col}</b>: {value_str}")
                    else:
                        annotations.append(f"<b>{col}</b>: N/A")
                else:
                    #print(f"Warning: Column {col} not found in protein_row for {protein}")
                    annotations.append(f"<b>{col}</b>: N/A")
            if all_empty:
                #print(f"All tooltip columns empty for protein {protein}, skipping tooltip")
                return
            tooltip_text = "<br>".join(annotations)
            self.tooltip = QGraphicsTextItem()
            self.tooltip.setHtml(tooltip_text)
            self.tooltip.setFont(QFont("Arial", 8))
            self.tooltip.setDefaultTextColor(QColor(0, 0, 0))
            max_width = 300
            self.tooltip.setTextWidth(max_width)
            tooltip_pos = self.mapToScene(QPointF(self.rect().right() + 5, self.rect().center().y()))
            self.tooltip.setPos(tooltip_pos)
            bounding_rect = self.tooltip.boundingRect()
            self.bg_rect = QGraphicsRectItem(bounding_rect.adjusted(-2, -2, 2, 2))
            self.bg_rect.setBrush(QColor(255, 255, 200, 200))
            self.bg_rect.setPen(QColor(0, 0, 0))
            self.bg_rect.setPos(tooltip_pos)
            self.bg_rect.setZValue(10)
            if self.scene():
                self.scene().addItem(self.bg_rect)
                self.scene().addItem(self.tooltip)
                self.tooltip.setZValue(11)
        except Exception as e:
            print(f"Error in HoverRectItem show_tooltip: {e}")
            self.cleanup_tooltip()

    def hoverEnterEvent(self, event):
        self.hover_timer.start(100)

    def hoverLeaveEvent(self, event):
        try:
            self.hover_timer.stop()
            self.cleanup_tooltip()
        except Exception as e:
            print(f"Error in HoverRectItem hoverLeaveEvent: {e}")

    def cleanup_tooltip(self):
        try:
            if self.tooltip is not None and self.scene():
                if self.bg_rect and self.bg_rect.scene():
                    self.scene().removeItem(self.bg_rect)
                if self.tooltip.scene():
                    self.scene().removeItem(self.tooltip)
            self.tooltip = None
            self.bg_rect = None
        except Exception as e:
            print(f"Error in HoverRectItem cleanup_tooltip: {e}")

class ProteinBoxItem(QGraphicsRectItem):
    def __init__(self, x, y, width, height, entry_id, viewer):
        super().__init__(x, y, width, height)
        self.entry_id = entry_id
        self.viewer = viewer
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setAcceptHoverEvents(False)  # Hover handled by hover_rect

        # Create an invisible hover rectangle (2/3 size, centered)
        hover_width = width * 2 / 3
        hover_height = height * 2 / 3
        hover_x = x + (width - hover_width) / 2
        hover_y = y + (height - hover_height) / 2
        self.hover_rect = HoverRectItem(hover_x, hover_y, hover_width, hover_height, entry_id, viewer, parent=self)

    def contextMenuEvent(self, event):
        menu = QMenu()
        entry_data = self.viewer.protein_data_map[self.entry_id]
        current_prot = self.viewer.current_proteins.get(self.entry_id)
        for protein in entry_data['proteins']:
            if protein in self.viewer.proteomic_data[self.viewer.hsa_id_column].values:
                action = menu.addAction(protein)
                action.setCheckable(True)
                action.setChecked(protein == current_prot['protein'])
                action.triggered.connect(lambda checked, p=protein:
                                         self.viewer.switch_protein(self.entry_id, p))
        menu.exec_(event.screenPos())

class PathwayViewer(QGraphicsView):
    def __init__(self, entries, image, proteomic_data, ptm_datasets, settings):
        super().__init__()
        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.setRenderHint(QPainter.Antialiasing)
        self.setRenderHint(QPainter.SmoothPixmapTransform)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorUnderMouse)
        self.ptm_items = {}
        self.multi_protein_indicators = {}
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.setStyleSheet("""
            QScrollBar:horizontal {
                border: none;
                background: #f0f0f0;
                height: 12px;
                margin: 0px 20px 0 20px;
            }
            QScrollBar::handle:horizontal {
                background: #666666;
                min-width: 20px;
                border-radius: 6px;
            }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
                background: none;
                width: 0px;
            }
            QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {
                background: none;
            }
            QScrollBar:vertical {
                border: none;
                background: #f0f0f0;
                width: 12px;
                margin: 20px 0 20px 0;
            }
            QScrollBar::handle:vertical {
                background: #666666;
                min-height: 20px;
                border-radius: 6px;
            }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                background: none;
                height: 0px;
            }
            QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
                background: none;
            }
        """)
        self.background_item = None
        self.protein_group = None
        self.proteomic_data = proteomic_data
        self.ptm_datasets = {f"ptm_{i}": dataset for i, dataset in enumerate(ptm_datasets)}
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
        self.protein_data_map = {}
        self.current_proteins = {}
        self.ptm_items = {}
        self.settings = settings
        self.protein_selection_option = settings.get('protein_selection_option', 2)
        self.ptm_selection_option = settings.get('ptm_selection_option', 2)
        self.fold_change_column = settings.get('fold_change_column', 'ER+_Est-_x-y_TNBC')
        self.hsa_id_column = settings.get('hsa_id_column', 'KEGG_hsa')
        self.prot_uniprot_column = settings.get('prot_uniprot_column', 'Uniprot_ID')
        self.gene_name_column = settings.get('gene_name_column', 'Gene Symbol')
        self.negative_color = settings.get('negative_color', (255, 0, 0))
        self.positive_color = settings.get('positive_color', (0, 0, 255))
        self.max_negative = settings.get('max_negative', -2)
        self.max_positive = settings.get('max_positive', 2)
        self.prot_label_font = settings.get('prot_label_font', 'Arial')
        self.prot_label_size = settings.get('prot_label_size', 6)
        self.ptm_label_font = settings.get('ptm_label_font', 'Arial')
        self.ptm_label_color = settings.get('ptm_label_color', QColor(0, 0, 0))
        self.ptm_label_size = settings.get('ptm_label_size', 4)
        self.ptm_circle_radius = settings.get('ptm_circle_radius', 5)
        self.ptm_circle_spacing = settings.get('ptm_circle_spacing', 2)
        self.show_background_image = settings.get('show_background_image', True)
        self.display_types = settings.get('display_types', ['gene'])
        self.show_groups = settings.get('show_groups', False)
        self.show_multi_protein_indicator = settings.get('show_multi_protein_indicator', True)
        if self.show_background_image and image:
            image = image.convert("RGBA")
            data = image.tobytes("raw", "RGBA")
            qimage = QImage(data, image.size[0], image.size[1], QImage.Format_RGBA8888)
            qpixmap = QPixmap.fromImage(qimage)
            self.background_item = QGraphicsPixmapItem(qpixmap)
            self.background_item.setPos(22, 9)
            self.scene.addItem(self.background_item)
            self.scene.setSceneRect(0, 0, image.size[0], image.size[1])
            self.image_width = image.size[0]
            self.image_height = image.size[1]
        else:
            self.scene.setSceneRect(0, 0, 0, 0)
            self.image_width = 0
            self.image_height = 0
        self.protein_group = QGraphicsItemGroup()
        self.scene.addItem(self.protein_group)
        self.save_button = QPushButton("Save as PDF", self)
        self.save_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                padding: 5px 10px;
                border: none;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        try:
            self.save_button.clicked.connect(self.save_to_pdf)
        except Exception as e:
            print(f"Error connecting save button: {e}")
        self.update_button_position()
        if self.show_background_image and image:
            self.fitInView(self.background_item, Qt.KeepAspectRatio)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.update_button_position()
        if self.show_background_image and self.background_item and \
                (self.isMaximized() or self.isFullScreen()):
            self.fitInView(self.background_item, Qt.KeepAspectRatio)

    def update_button_position(self):
        button_width = self.save_button.width()
        button_height = self.save_button.height()
        self.save_button.move(
            self.width() - button_width - 20,
            self.height() - button_height - 20
        )

    def save_to_pdf(self):
        try:
            pathway_id = self.settings.get('pathway_id', 'hsa04010')
            default_filename = f"{pathway_id}_pathway_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
            output_file, _ = QFileDialog.getSaveFileName(
                self, "Save Pathway as PDF", default_filename, "PDF Files (*.pdf)")
            if not output_file:
                return
            pdf_writer = QPdfWriter(output_file)
            pdf_writer.setResolution(300)
            content_rect = None
            for item in self.scene.items():
                if not isinstance(item, QGraphicsItem):
                    continue
                item_rect = item.sceneBoundingRect()
                content_rect = item_rect if content_rect is None else content_rect.united(item_rect)
            if content_rect is None:
                print("No content to save.")
                return
            if self.show_background_image and self.background_item:
                bg_pos = self.background_item.pos()
                content_rect.adjust(-bg_pos.x(), -bg_pos.y(), 0, 0)
                bg_rect = self.background_item.sceneBoundingRect()
                content_rect = content_rect.united(bg_rect)
            content_width = content_rect.width()
            content_height = content_rect.height()
            margin_mm = 5
            mm_per_pixel = 25.4 / 300
            content_width_mm = content_width * mm_per_pixel
            content_height_mm = content_height * mm_per_pixel
            pdf_width_mm = content_width_mm + 2 * margin_mm
            pdf_height_mm = content_height_mm + 2 * margin_mm
            custom_size = QPageSize(QSizeF(pdf_width_mm, pdf_height_mm), QPageSize.Millimeter, "Custom")
            pdf_writer.setPageSize(custom_size)
            pdf_writer.setPageMargins(QMarginsF(0, 0, 0, 0))
            painter = QPainter(pdf_writer)
            painter.setRenderHint(QPainter.Antialiasing)
            target_width = (pdf_width_mm - 2 * margin_mm) / mm_per_pixel
            target_height = (pdf_height_mm - 2 * margin_mm) / mm_per_pixel
            target_rect = QRectF(margin_mm / mm_per_pixel, margin_mm / mm_per_pixel, target_width, target_height)
            adjusted_content_rect = QRectF(0, 0, content_width, content_height)
            self.scene.render(painter, target_rect, adjusted_content_rect)
            painter.end()
            print(f"Pathway saved as PDF: {output_file}")
        except (IOError, PermissionError) as e:
            show_error_popup(f"Failed to save PDF: {str(e)}", "Save Error", parent=self)
            print(f"Error saving PDF (file issue): {e}")
        except Exception as e:
            show_error_popup(f"Unexpected error: {str(e)}", "Save Error", parent=self)
            print(f"Unexpected error saving PDF: {e}")

    def get_color(self, fold_change):
        if pd.isna(fold_change):
            return QColor(128, 128, 128)
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
        return QColor(r, g, b)

    def choose_protein(self, proteins):
        #print(f"Choosing protein from: {proteins}")
        #print(f"Using hsa_id_column: {self.hsa_id_column}, prot_uniprot_column: {self.prot_uniprot_column}")
        valid_proteins = [p for p in proteins if p in self.proteomic_data[self.hsa_id_column].values]
        #print(f"Valid proteins found: {valid_proteins}")
        if not valid_proteins:
            #print("No valid proteins found in proteomic data.")
            return None
        try:
            if self.protein_selection_option == 1:
                selected_protein = valid_proteins[0]
            elif self.protein_selection_option == 2:
                max_fold_change = -float('inf')
                selected_protein = None
                for protein in valid_proteins:
                    fold_change = self.proteomic_data.loc[
                        self.proteomic_data[self.hsa_id_column] == protein, self.fold_change_column].values
                    if len(fold_change) > 0 and abs(fold_change[0]) > max_fold_change:
                        max_fold_change = abs(fold_change[0])
                        selected_protein = protein
                if not selected_protein:
                    selected_protein = valid_proteins[0]
            else:
                selected_protein = valid_proteins[0]
            #print(f"Selected protein: {selected_protein}")
            uniprot_id = self.proteomic_data.loc[
                self.proteomic_data[self.hsa_id_column] == selected_protein, self.prot_uniprot_column].values
            uniprot_id = uniprot_id[0] if len(uniprot_id) > 0 else None
            #print(f"Retrieved UniProt ID: {uniprot_id}")
            return {'protein': selected_protein, 'uniprot_id': uniprot_id}
        except (KeyError, IndexError) as e:
            #print(f"Error choosing protein: {str(e)}")
            return None

    def prioritize_ptm_sites(self, uniprot_id):
        #print(f"Prioritizing PTM sites for UniProt ID: {uniprot_id}")
        all_ptms = []
        for dataset_id, dataset in self.ptm_datasets.items():
            uniprot_col = dataset['uniprot_column']
            fold_change_col = dataset['fold_change_column']
            ptm_data = dataset['data']
            #print(f"Checking dataset {dataset_id} for UniProt ID {uniprot_id} in column {uniprot_col}")
            ptm_subset = ptm_data[ptm_data[uniprot_col] == uniprot_id].copy()
            if ptm_subset.empty:
                #print(f"No PTM sites found for {uniprot_id} in dataset {dataset_id}")
                continue
            #print(f"Found {len(ptm_subset)} PTM sites for {uniprot_id} in dataset {dataset_id}")
            #print(f"PTM subset for {uniprot_id}:\n{ptm_subset[['T: Uniprot_ID', 'T: Site Position', fold_change_col]]}")
            ptm_subset['dataset_id'] = dataset_id
            ptm_subset['is_phospho'] = dataset['type'] == 'Phosphorylation'
            ptm_subset['is_modulating'] = ptm_subset.get(dataset.get('modulation_column', ''), '') == '+' if dataset['type'] == 'Phosphorylation' else True
            all_ptms.append(ptm_subset)
        if not all_ptms:
            #print("No PTM sites found for this UniProt ID across all datasets.")
            return pd.DataFrame()
        combined_ptms = pd.concat(all_ptms, ignore_index=True)
        #print(f"Combined PTM dataframe has {len(combined_ptms)} rows")
        #print(f"Combined PTM columns: {combined_ptms.columns.tolist()}")
        if combined_ptms.empty:
            #print("Combined PTM dataframe is empty after concatenation.")
            return combined_ptms
        if 'dataset_id' not in combined_ptms.columns:
            #print("Error: 'dataset_id' column missing in combined_ptms")
            return pd.DataFrame()
        if combined_ptms['dataset_id'].isna().any():
            #print(f"Warning: {combined_ptms['dataset_id'].isna().sum()} rows have NaN in dataset_id")
            combined_ptms = combined_ptms.dropna(subset=['dataset_id'])
        if combined_ptms.empty:
            #print("Combined PTM dataframe is empty after dropping NaN dataset_id rows.")
            return combined_ptms
        try:
            if self.ptm_selection_option == 1:
                print("Applying PTM selection option 1: Highest absolute fold change")
                combined_ptms = combined_ptms.sort_values(
                    by=[dataset['fold_change_column'] for dataset in self.ptm_datasets.values()],
                    key=lambda x: x.abs(),
                    ascending=False
                )
            elif self.ptm_selection_option == 2:
                print("Applying PTM selection option 2: Prefer modulating phospho PTMs")
                def sort_key(row):
                    try:
                        dataset_id = row['dataset_id']
                        dataset = self.ptm_datasets[dataset_id]
                        fold_change = row[dataset['fold_change_column']]
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
                    by=[dataset['fold_change_column'] for dataset in self.ptm_datasets.values()],
                    key=lambda x: x.abs(),
                    ascending=False
                )
            print(f"Selected {len(combined_ptms.head(self.settings.get('ptm_max_display', 4)))} PTM sites")
            return combined_ptms.head(self.settings.get('ptm_max_display', 4))
        except Exception as e:
            print(f"Error prioritizing PTM sites for UniProt ID {uniprot_id}: {str(e)}")
            return pd.DataFrame()
        print(f"Combined PTM dataframe has {len(combined_ptms)} rows")
        print(f"PTM selection option: {self.ptm_selection_option}")
        print(f"Final selected PTM sites: {len(combined_ptms.head(self.settings.get('ptm_max_display', 4)))}")
        print(combined_ptms.head(self.settings.get('ptm_max_display', 4)))

    def display_pathway(self, entries, proteomic_data, ptm_datasets):
        COMPOUND_X_OFFSET = 0
        COMPOUND_Y_OFFSET = 0
        MAP_X_OFFSET = 0
        MAP_Y_OFFSET = 0
        try:
            exceptions = parse_exceptions_file()
            proximities, protein_coords, protein_to_entry_ids = calculate_proximities(entries, self.proteomic_data)
            self.settings['parsed_exceptions'] = exceptions
            self.settings['proximities'] = proximities
            self.settings['protein_to_entry_ids'] = protein_to_entry_ids
            hsa_exceptions = exceptions['hsa_specific'].get(self.settings.get('pathway_id', 'hsa04010'),
                                                            {'proximity': [], 'specific': {}})
            global_exceptions = exceptions['global']
            self.save_button.raise_()
            self.save_button.show()
            chosen_proteins = set()
            kgml_file = download_kegg_kgml(self.settings.get('pathway_id', 'hsa04010'))
            tree = ET.parse(kgml_file)
            root = tree.getroot()
            entries = []
            groups = []
            for entry in root.findall("entry"):
                entry_id = entry.get("id")
                entry_type = entry.get("type")
                name = entry.get("name")
                graphics = entry.find("graphics")
                if graphics is not None:
                    x = int(graphics.get("x"))
                    y = int(graphics.get("y"))
                    width = int(graphics.get("width"))
                    height = int(graphics.get("height"))
                    graphics_name = graphics.get("name", "")
                    first_name = graphics_name.split(",")[0].strip() if graphics_name else name
                    if entry_type == "group":
                        fgcolor = graphics.get("fgcolor", "#000000")
                        bgcolor = graphics.get("bgcolor", "#FFFFFF")
                        components = [comp.get("id") for comp in entry.findall("component")]
                        groups.append({
                            "id": entry_id,
                            "x": x,
                            "y": y,
                            "width": width,
                            "height": height,
                            "fgcolor": fgcolor,
                            "bgcolor": bgcolor,
                            "components": components
                        })
                    else:
                        entries.append({
                            "id": entry_id,
                            "name": name,
                            "x": x,
                            "y": y,
                            "width": width,
                            "height": height,
                            "first_name": first_name,
                            "type": entry_type
                        })
            if self.show_groups:
                for group in groups:
                    rect = QGraphicsRectItem(group["x"], group["y"], group["width"], group["height"])
                    rect.setBrush(QColor(group["bgcolor"]))
                    rect.setPen(QColor(group["fgcolor"]))
                    rect.setZValue(0)
                    self.scene.addItem(rect)
            relations = root.findall("relation")
            if isinstance(self.display_types, str):
                display_types_set = {self.display_types}
            else:
                display_types_set = set(self.display_types)
            if 'compound' in display_types_set or 'map' in display_types_set:
                for relation in relations:
                    entry1 = relation.get("entry1")
                    entry2 = relation.get("entry2")
                    entry1_data = find_entry_or_group(entry1, entries, groups)
                    entry2_data = find_entry_or_group(entry2, entries, groups)
                    if entry1_data and entry2_data:
                        if (entry1_data["type"] not in display_types_set and entry1_data["type"] != "group") or \
                                (entry2_data["type"] not in display_types_set and entry2_data["type"] != "group"):
                            continue
                        if entry1_data["type"] == "compound":
                            x1 = entry1_data["x"] + entry1_data["width"] / 2 + COMPOUND_X_OFFSET
                            y1 = entry1_data["y"] + entry1_data["height"] / 2 + COMPOUND_Y_OFFSET
                        elif entry1_data["type"] == "map":
                            x1 = entry1_data["x"] + entry1_data["width"] / 2 + MAP_X_OFFSET
                            y1 = entry1_data["y"] + entry1_data["height"] / 2 + MAP_Y_OFFSET
                        else:
                            x1 = entry1_data["x"] + entry1_data["width"] / 2
                            y1 = entry1_data["y"] + entry1_data["height"] / 2
                        if entry2_data["type"] == "compound":
                            x2 = entry2_data["x"] + entry2_data["width"] / 2 + COMPOUND_X_OFFSET
                            y2 = entry2_data["y"] + entry2_data["height"] / 2 + COMPOUND_Y_OFFSET
                        elif entry2_data["type"] == "map":
                            x2 = entry2_data["x"] + entry2_data["width"] / 2 + MAP_X_OFFSET
                            y2 = entry2_data["y"] + entry2_data["height"] / 2 + MAP_Y_OFFSET
                        else:
                            x2 = entry2_data["x"] + entry2_data["width"] / 2
                            y2 = entry2_data["y"] + entry2_data["height"] / 2
                        angle = calculate_angle(x1, y1, x2, y2)
                        side = determine_arrow_side(angle)
                        x2_adj, y2_adj = adjust_arrow_endpoint(entry2_data["x"], entry2_data["y"],
                                                               entry2_data["width"], entry2_data["height"], side)
                        start_point = QPointF(x1, y1)
                        end_point = QPointF(x2_adj, y2_adj)
                        arrow = ArrowItem(start_point, end_point)
                        self.scene.addItem(arrow)
            for entry in entries:
                entry_type = entry["type"]
                if entry_type not in display_types_set:
                    continue
                if entry_type == 'gene':
                    proteins = entry["name"].split()
                    self.protein_data_map[entry["id"]] = {
                        'proteins': proteins,
                        'x': entry["x"],
                        'y': entry["y"],
                        'width': entry["width"],
                        'height': entry["height"],
                        'first_name': entry["first_name"]
                    }
                    chosen_protein = self.choose_protein(proteins)
                    if chosen_protein:
                        chosen_proteins.add(chosen_protein['protein'])
                        self.current_proteins[entry["id"]] = chosen_protein
                        self.draw_protein_box(entry["id"], chosen_protein)
                    else:
                        fallback_protein = proteins[0] if proteins else "Unknown"
                        self.current_proteins[entry["id"]] = {'protein': fallback_protein, 'uniprot_id': None}
                        self.draw_protein_box(entry["id"], {'protein': fallback_protein, 'uniprot_id': None})
            if self.show_background_image:
                max_x = max(entry["x"] + entry["width"] for entry in entries if entry["type"] in display_types_set)
                max_y = max(entry["y"] + entry["height"] for entry in entries if entry["type"] in display_types_set)
                if max_x > self.image_width or max_y > self.image_height:
                    print(f"Warning: Some elements (max_x={max_x}, max_y={max_y}) exceed image bounds "
                          f"({self.image_width}x{self.image_height}). They may be clipped.")
            temp_file = savefile(self.settings.get('pathway_id', 'hsa04010'),
                                 self.settings.get('output_subdir', 'output/testing_file_001'),
                                 '_visprots', extension='.txt', include_timestamp=False)
            try:
                with open(temp_file, 'w') as f:
                    for protein in chosen_proteins:
                        f.write(f"{protein}\n")
            except (IOError, PermissionError) as e:
                print(f"Error writing to temporary file {temp_file}: {e}")
        except Exception as e:
            show_error_popup(f"Error in display_pathway: {e}", "Error", parent=self)
            print(f"Error in display_pathway: {e}")

    def draw_protein_box(self, entry_id, protein_dict):
        try:
            # Clean up existing PTM items and indicators
            if entry_id in self.ptm_items:
                for item in self.ptm_items[entry_id]:
                    self.scene.removeItem(item)
                del self.ptm_items[entry_id]
            if entry_id in self.multi_protein_indicators:
                for item in self.multi_protein_indicators[entry_id]:
                    self.scene.removeItem(item)
                del self.multi_protein_indicators[entry_id]

            entry_data = self.protein_data_map[entry_id]
            x = entry_data['x']
            y = entry_data['y']
            width = entry_data['width']
            height = entry_data['height']
            proteins = entry_data['proteins']
            protein = protein_dict['protein']
            uniprot_id = protein_dict['uniprot_id']

            # Create the protein box
            box = ProteinBoxItem(x, y, width, height, entry_id, self)
            protein_in_data = protein in self.proteomic_data[self.hsa_id_column].values
            valid_proteins = [p for p in proteins if p in self.proteomic_data[self.hsa_id_column].values]

            if protein_in_data:
                fold_change = self.proteomic_data.loc[
                    self.proteomic_data[self.hsa_id_column] == protein, self.fold_change_column].values[0]
                gene_name = self.proteomic_data.loc[
                    self.proteomic_data[self.hsa_id_column] == protein, self.gene_name_column].values[0]
                color = self.get_color(fold_change)
            else:
                gene_name = entry_data['first_name'].split(",")[0].strip()
                color = QColor(128, 128, 128)

            box.setBrush(color)
            box.setPen(QColor(0, 0, 0))
            self.scene.addItem(box)

            # Add protein label
            label = QGraphicsTextItem(gene_name)
            label.setFont(QFont(self.prot_label_font, self.prot_label_size))
            label.setDefaultTextColor(QColor(0, 0, 0))
            bounding_rect = label.boundingRect()
            adjusted_x = x + (width / 2) - (bounding_rect.width() / 2)
            adjusted_y = y + (height / 2) - (bounding_rect.height() / 2)
            label.setPos(adjusted_x, adjusted_y)
            label.setParentItem(box)

            if protein == 'MAPK1' or 'MAPK1' in proteins:
                print(
                    f"Entry {entry_id}: Proteins = {proteins}, Chosen = {protein}, UniProt = {uniprot_id}, Valid = {valid_proteins}")

            # Initialize PTM items list
            self.ptm_items[entry_id] = []

            # Render PTMs if available
            if protein_in_data and uniprot_id:
                ptm_sites = self.prioritize_ptm_sites(uniprot_id)
                if not ptm_sites.empty:
                    #print(f"Rendering {len(ptm_sites)} PTM sites for {protein} (UniProt: {uniprot_id})")
                    if protein == 'MAPK1':
                        #print(f"MAPK1 PTM sites to render ({len(ptm_sites)}):")
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
                    exceptions = self.settings.get('parsed_exceptions', {})
                    proximities = self.settings.get('proximities', {})
                    protein_to_entry_ids = self.settings.get('protein_to_entry_ids', {})
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
                            #print(f"Specific rule: Blocking {rule[1]} for {protein}")
                        elif action == 'priority':
                            custom_priority = rule[1:]
                            #print(f"Specific rule: Setting priority {custom_priority} for {protein}")
                        elif action == 'move_circle' and len(rule) == 4:
                            pos, dx, dy = rule[1], float(rule[2]), float(rule[3])
                            move_circle[pos] = (dx, dy)
                            #print(f"Specific rule: Moving circle {pos} by ({dx}, {dy}) for {protein}")
                        elif action == 'move_label' and len(rule) == 5:
                            pos, dx, dy, align = rule[1], float(rule[2]), float(rule[3]), rule[4]
                            move_label[pos] = (dx, dy, align)
                            #print(f"Specific rule: Moving label {pos} by ({dx}, {dy}, {align}) for {protein}")
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
                            #print(f"Rule passed for {entry_id} ({protein}): {rule}")
                            passed_rules_by_action[rule['action']].append(rule)
                    for action in ['block', 'priority']:
                        if passed_rules_by_action[action]:
                            highest_rule = max(passed_rules_by_action[action], key=lambda r: r['priority'])
                            if action == 'block':
                                for pos in highest_rule['values']:
                                    blocked_positions.add(pos)
                                #print(f"Proximity rule: Blocking {highest_rule['values']} for {protein}")
                            elif action == 'priority':
                                custom_priority = highest_rule['values']
                                #print(f"Proximity rule: Setting priority {custom_priority} for {protein}")
                    for action in ['move_circle', 'move_label']:
                        for rule in sorted(passed_rules_by_action[action],
                                           key=lambda r: (r['priority'], proximity_rules.index(r))):
                            if action == 'move_circle' and len(rule['values']) == 3:
                                pos, dx, dy = rule['values'][0], float(rule['values'][1]), float(rule['values'][2])
                                move_circle[pos] = (dx, dy)
                                #print(f"Proximity rule: Moving circle {pos} by ({dx}, {dy}) for {protein}")
                            elif action == 'move_label' and len(rule['values']) == 4:
                                pos, dx, dy, align = rule['values'][0], float(rule['values'][1]), float(
                                    rule['values'][2]), rule['values'][3]
                                move_label[pos] = (dx, dy, align)
                                #print(f"Proximity rule: Moving label {pos} by ({dx}, {dy}, {align}) for {protein}")
                    available_positions = [p for p in (custom_priority or position_priority) if
                                           p not in blocked_positions]
                    #print(f"Available positions for {entry_id} ({protein}): {available_positions}")
                    for i, (_, row) in enumerate(ptm_sites.iterrows()):
                        if i >= len(available_positions):
                            break
                        pos_key = available_positions[i]
                        circle_x, circle_y = positions[pos_key]
                        dataset_id = row['dataset_id']
                        dataset = self.ptm_datasets[dataset_id]
                        if pos_key in move_circle:
                            circle_x += move_circle[pos_key][0]
                            circle_y += move_circle[pos_key][1]
                            #print(f"Applying move_circle to {pos_key}: ({move_circle[pos_key][0]}, {move_circle[pos_key][1]})")
                        ptm_fold_change = row[dataset['fold_change_column']]
                        ptm_shape = dataset.get('shape', 'Circle')
                        ptm = PTMItem(
                            circle_x - self.ptm_circle_radius,
                            circle_y - self.ptm_circle_radius,
                            self.ptm_circle_radius * 2,
                            self.ptm_circle_radius * 2,
                            row,
                            dataset_id,
                            self,
                            shape=ptm_shape,
                            parent=box
                        )
                        self.ptm_items[entry_id].append(ptm)
                        self.scene.addItem(ptm)
                        classification = row.get('Phosphosite_Classification', 'none')
                        #print(f"PTM site {row[dataset['site_column']]} classification: {classification}")
                        if classification != 'none':
                            symbol_dict_key = f"{classification}_dict"
                            symbol_dict = next((d[symbol_dict_key] for d in dataset.get('ptm_symbol_list', []) if
                                                symbol_dict_key in d), None)
                            if symbol_dict:
                                symbol = symbol_dict.get('symbol', '+')
                                font = symbol_dict.get('symbol_font', self.ptm_label_font)
                                size = symbol_dict.get('symbol_size', self.ptm_label_size)
                                color = QColor(*symbol_dict.get('symbol_color', (0, 0, 0)))
                                x_offset = symbol_dict.get('symbol_x_offset', 0)
                                y_offset = symbol_dict.get('symbol_y_offset', 0)
                                print(
                                    f"Rendering symbol '{symbol}' for classification {classification} at position {pos_key}")
                                try:
                                    sign = QGraphicsTextItem(symbol, parent=box)
                                    sign.setFont(QFont(font, size, QFont.Bold))
                                    sign.setDefaultTextColor(color)
                                    sign.setAcceptHoverEvents(False)
                                    sign.setZValue(12)
                                    text_rect = sign.boundingRect()
                                    if text_rect.width() == 0 or text_rect.height() == 0:
                                        raise ValueError(f"Font '{font}' does not support symbol '{symbol}'")
                                    adjusted_x = circle_x - (text_rect.width() / 2) + x_offset
                                    adjusted_y = circle_y - (text_rect.height() / 2) + y_offset
                                    sign.setPos(adjusted_x, adjusted_y)
                                    self.ptm_items[entry_id].append(sign)
                                    self.scene.addItem(sign)
                                except Exception as e:
                                    error_msg = f"Font Error: {str(e)}. Please choose a font that supports '{symbol}' (e.g., Arial)."
                                    QMessageBox.critical(
                                        self,
                                        "Font Unsupported",
                                        error_msg,
                                        QMessageBox.Ok
                                    )
                                    sign = QGraphicsTextItem(symbol, parent=box)
                                    sign.setFont(QFont("Arial", size, QFont.Bold))
                                    sign.setDefaultTextColor(color)
                                    sign.setAcceptHoverEvents(False)
                                    sign.setZValue(12)
                                    text_rect = sign.boundingRect()
                                    adjusted_x = circle_x - (text_rect.width() / 2) + x_offset
                                    adjusted_y = circle_y - (text_rect.height() / 2) + y_offset
                                    sign.setPos(adjusted_x, adjusted_y)
                                    self.ptm_items[entry_id].append(sign)
                                    self.scene.addItem(sign)
                            else:
                                print(
                                    f"Warning: No symbol dict found for classification {classification} in dataset {dataset_id}")
                        site_number = row[dataset['site_column']]
                        if pos_key in move_label:
                            x_offset, y_offset, ha = move_label[pos_key]
                            #print(f"Applying move_label to {pos_key}: ({x_offset}, {y_offset}, {ha})")
                        else:
                            defaults = {
                                'N1': (-3, -5, 'right'), 'N2': (0, -11, 'center'), 'N3': (3, -5, 'left'),
                                'S1': (-3, 5, 'right'), 'S2': (0, 12, 'center'), 'S3': (3, 5, 'left'),
                                'W1': (-3, -2, 'right'), 'W2': (-3, 2, 'right'),
                                'E1': (3, -2, 'left'), 'E2': (3, 2, 'left')
                            }
                            x_offset, y_offset, ha = defaults[pos_key]
                        text = QGraphicsTextItem(str(site_number), parent=box)
                        text.setFont(QFont(self.ptm_label_font, self.ptm_label_size))
                        label_color_tuple = self.settings.get('ptm_label_color', (0, 0, 0))
                        if not isinstance(label_color_tuple, tuple) or len(label_color_tuple) != 3:
                            print(f"Warning: Invalid ptm_label_color {label_color_tuple}, using default black")
                            label_color_tuple = (0, 0, 0)
                        try:
                            label_color = QColor(*label_color_tuple)
                            if not label_color.isValid():
                                raise ValueError("Invalid color values")
                        except (ValueError, TypeError) as e:
                            print(
                                f"Error: Failed to create QColor from {label_color_tuple}: {str(e)}, using default black")
                            label_color = (0, 0, 0)
                        text.setDefaultTextColor(label_color)
                        bounding_rect = text.boundingRect()
                        if ha == 'center':
                            adjusted_x = circle_x + x_offset - (bounding_rect.width() / 2)
                            adjusted_y = circle_y + y_offset - (bounding_rect.height() / 2)
                        elif ha == 'right':
                            adjusted_x = circle_x + x_offset - bounding_rect.width()
                            adjusted_y = circle_y + y_offset - (bounding_rect.height() / 2)
                        else:
                            adjusted_x = circle_x + x_offset
                            adjusted_y = circle_y + y_offset - (bounding_rect.height() / 2)
                        text.setPos(adjusted_x, adjusted_y)
                        self.ptm_items[entry_id].append(text)
                        self.scene.addItem(text)

            # Add multi-protein indicator regardless of PTM presence
            if self.show_multi_protein_indicator and len(valid_proteins) > 1:
                self.multi_protein_indicators[entry_id] = []
                multi_indicator = QGraphicsTextItem("*")
                multi_indicator.setFont(QFont(self.prot_label_font, self.prot_label_size + 2, QFont.Bold))
                multi_indicator.setDefaultTextColor(QColor(0, 0, 0))
                multi_indicator.setZValue(0)
                bounding_rect = multi_indicator.boundingRect()
                adjusted_x = x + width - bounding_rect.width() + 3
                adjusted_y = y - 6
                multi_indicator.setPos(adjusted_x, adjusted_y)
                multi_indicator.setParentItem(box)
                self.multi_protein_indicators[entry_id].append(multi_indicator)
                self.scene.addItem(multi_indicator)
                if protein == 'MAPK1' or 'MAPK1' in proteins:
                    print(f"Added multi-protein indicator (*) for entry {entry_id} to multi_protein_indicators")

        except Exception as e:
            print(f"Error drawing protein box for {entry_id}: {e}")
            return

    def switch_protein(self, entry_id, new_protein):
        uniprot_id = self.proteomic_data.loc[
            self.proteomic_data[self.hsa_id_column] == new_protein, self.prot_uniprot_column].values
        uniprot_id = uniprot_id[0] if len(uniprot_id) > 0 else None
        self.current_proteins[entry_id] = {'protein': new_protein, 'uniprot_id': uniprot_id}
        self.draw_protein_box(entry_id, self.current_proteins[entry_id])

    def wheelEvent(self, event):
        current_scale = self.transform().m11()
        content_rect = self.scene.itemsBoundingRect()
        viewport_width = self.viewport().width()
        viewport_height = self.viewport().height()
        content_width = content_rect.width() * current_scale
        content_height = content_rect.height() * current_scale
        if event.modifiers() & Qt.ControlModifier:
            zoom_in_factor = 1.1
            zoom_out_factor = 1 / zoom_in_factor
            if event.angleDelta().y() > 0:
                self.scale(zoom_in_factor, zoom_in_factor)
            elif event.angleDelta().y() < 0:
                new_scale = current_scale * zoom_out_factor
                if self.show_background_image and new_scale < 1.0:
                    self.fitInView(self.background_item, Qt.KeepAspectRatio)
                else:
                    self.scale(zoom_out_factor, zoom_out_factor)
            event.accept()
            return
        scrolled = False
        if event.angleDelta().x() != 0 and content_width > viewport_width:
            scroll_step = -event.angleDelta().x() * 0.5
            horizontal_bar = self.horizontalScrollBar()
            current_value = horizontal_bar.value()
            new_value = current_value + int(scroll_step)
            new_value = max(horizontal_bar.minimum(), min(new_value, horizontal_bar.maximum()))
            horizontal_bar.setValue(new_value)
            scrolled = True
        if event.angleDelta().y() != 0 and content_height > viewport_height:
            scroll_step = -event.angleDelta().y() * 0.5
            vertical_bar = self.verticalScrollBar()
            current_value = vertical_bar.value()
            new_value = current_value + int(scroll_step)
            new_value = max(vertical_bar.minimum(), min(new_value, vertical_bar.maximum()))
            vertical_bar.setValue(new_value)
            scrolled = True
        if scrolled:
            event.accept()

def visualize_kegg_pathway(pathway_id, data, settings):
    try:
        proteomic_data = pd.read_csv(data['protein']['file_path'], sep="\t")
        ptm_datasets = data['ptm']
        for dataset in ptm_datasets:
            dataset['fold_change_column'] = dataset['main_columns'][0][1] if dataset['main_columns'] else 'ER+_Est-_x-y_TNBC'
        kgml_file = download_kegg_kgml(pathway_id)
        print(f"KGML file downloaded: {kgml_file}")
        pathway_image = download_kegg_image(pathway_id) if settings.get('show_background_image', True) else None
        print("Pathway image downloaded" if pathway_image else "No background image")
        entries, groups = parse_kgml(kgml_file)
        print(f"Parsed {len(entries)} entries and {len(groups)} groups")
        viewer = PathwayViewer(entries, pathway_image, proteomic_data, ptm_datasets, settings)
        print("PathwayViewer created")
        viewer.display_pathway(entries, proteomic_data, ptm_datasets)
        print("Pathway displayed")
        return viewer
    except Exception as e:
        show_error_popup(f"Error: {e}", "File Error")
        print(f"Error: {e}")
        return None

if __name__ == "__main__":
    settings = {
        'pathway_id': 'hsa04150',
        'protein_selection_option': 2,
        'ptm_selection_option': 2,
        'ptm_max_display': 4,
        'show_background_image': True,
        'display_types': ['gene'],
        'show_groups': False,
        'show_multi_protein_indicator': True,
        'negative_color': (255, 0, 0),
        'positive_color': (0, 0, 255),
        'max_negative': -2,
        'max_positive': 2,
        'prot_label_font': 'Arial',
        'prot_label_size': 6,
        'ptm_label_font': 'Arial',
        'ptm_label_color': (0, 0, 0),  # Changed from QColor to tuple
        'ptm_label_size': 4,
        'ptm_circle_radius': 5,
        'ptm_circle_spacing': 2,
        'protein_tooltip_columns': ['Gene Symbol', 'Uniprot_ID', 'ER+_Est-_x-y_TNBC']
    }
    data = {
        "protein": {
            "file_path": r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\Phosmap\BCp2_ProtMaphsaanno.txt",
            "uniprot_column": "Uniprot_ID",
            "kegg_column": "KEGG_hsa",
            "gene_column": "Gene Symbol",
            "main_columns": ["ER+_Est-_x-y_TNBC"]
        },
        "ptm": [
            {
                "type": "Phosphorylation",
                "file_path": r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\Phosmap\BCp2_PhosMaphsanno.txt",
                "uniprot_column": "T: Uniprot_ID",
                "site_column": "T: Site Position",
                "shape": "Circle",
                "main_columns": [("ER+_Est-_x-y_TNBC", "ER+_Est-_x-y_TNBC")],
                "modulation_column": "C: Regulatory site",
                "tooltip_columns": ["C: Regulatory site function", "C: Regulatory site process"],
                "ptm_symbol_list": [
                    {
                        "symbol_label_1_dict": {
                            "symbol": "↑",
                            "symbol_font": "Segoe UI Symbol",
                            "symbol_size": 6,
                            "symbol_color": (0, 0, 0),
                            "symbol_x_offset": 0,
                            "symbol_y_offset": -1,
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
                            "symbol_size": 6,
                            "symbol_color": (0, 0, 0),
                            "symbol_x_offset": 0,
                            "symbol_y_offset": -1,
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
                            "symbol_size": 6,
                            "symbol_color": (0, 0, 0),
                            "symbol_x_offset": 0,
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
                            "symbol_size": 6,
                            "symbol_color": (0, 0, 0),
                            "symbol_x_offset": 0,
                            "symbol_y_offset": 0,
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
    app = QApplication(sys.argv)
    viewer = visualize_kegg_pathway(settings['pathway_id'], data, settings)
    if viewer:
        viewer.show()
        sys.exit(app.exec_())