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
    QPushButton, QFileDialog, QGraphicsLineItem  # Added QGraphicsLineItem
)
from PyQt5.QtGui import (
    QPixmap, QColor, QFont, QPainter, QImage, QPolygonF, QPen,
    QPdfWriter, QPagedPaintDevice, QPageSize, QBrush
)
from PyQt5.QtCore import Qt, QRectF, QPointF, QSize, QSizeF, QMarginsF, QLineF  # Added QMarginsF here
import pywikipathways as pwp
import requests
from io import BytesIO



import requests
from bs4 import BeautifulSoup  # Add this import at the top of your file


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

def parse_exceptions_file(file_path):
    exceptions = {
        'global': {'proximity': [], 'specific': {}},
        'hsa_specific': {}
    }
    if not os.path.exists(file_path):
        return exceptions
    with open(file_path, 'r') as f:
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
    print("Parsed exceptions:", exceptions)
    return exceptions

def calculate_proximities(genes, proteomic_data):
    protein_coords = {}  # Keyed by entry_id
    proximities = {}     # Keyed by entry_id
    protein_to_entry_ids = defaultdict(list)  # Map protein names to their entry_ids

    # Populate coordinates and mapping
    for gene in genes:
        entry_id = gene["id"]
        proteins = gene["name"].split()
        protein_coords[entry_id] = {
            'x': gene["x"],
            'y': gene["y"],
            'width': gene["width"],
            'height': gene["height"],
            'proteins': proteins  # Store all proteins for this entry
        }
        for protein in proteins:
            protein_to_entry_ids[protein].append(entry_id)

    # Calculate proximities between all entry_ids
    for entry_id1, coords1 in protein_coords.items():
        proximities[entry_id1] = {}
        for entry_id2, coords2 in protein_coords.items():
            if entry_id1 == entry_id2:
                continue
            # Calculate edge-to-edge distances (unchanged logic)
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

def download_kegg_kgml(pathway_id):
    kgml_url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
    response = requests.get(kgml_url)
    if response.status_code == 200:
        output_file = savefile(pathway_id, 'kegg_pathway', '', extension=".xml", include_timestamp=False)
        with open(output_file, "w") as f:
            f.write(response.text)
        return output_file
    else:
        raise Exception(f"Failed to download KGML file: {response.status_code}")

def download_kegg_image(pathway_id):
    img_url = f"https://rest.kegg.jp/get/{pathway_id}/image"
    response = requests.get(img_url)
    if response.status_code == 200:
        return Image.open(io.BytesIO(response.content))
    else:
        raise Exception(f"Failed to download pathway image: {response.status_code}")

def parse_kgml(kgml_file):
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
    return entries, groups


def find_entry_or_group(entry_id, entries, groups):
    # Check entries first
    for entry in entries:
        if entry["id"] == entry_id:
            return entry

    # Check groups
    for group in groups:
        if group["id"] == entry_id:
            # If the ID matches a group, find a representative DataNode from its components
            if group["components"]:
                # Find the first DataNode in the group
                for entry in entries:
                    if entry["id"] == group["components"][0]:
                        print(f"Resolved group {entry_id} to DataNode {entry['id']}")
                        return entry
            print(f"Group {entry_id} has no components, cannot resolve to a DataNode")
            return None

    print(f"Could not find entry or group with ID {entry_id}")
    return None


def download_wikipathways_image(wp_id):
    # WikiPathways REST API
    image_url = f"https://webservice.wikipathways.org/getPathwayAs?fileType=png&pwId={wp_id}"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'image/png'
    }
    response = requests.get(image_url, headers=headers, allow_redirects=True)
    print(f"Image URL: {image_url}")
    print(f"Response status code: {response.status_code}")
    print(f"Response content type: {response.headers.get('Content-Type')}")
    print(f"Response content length: {len(response.content)} bytes")
    print(f"Response history: {len(response.history)} redirects")
    if response.history:
        for r in response.history:
            print(f"Redirect from: {r.status_code} - {r.url}")

    if response.status_code == 200:
        if "image/png" not in response.headers.get('Content-Type', ''):
            print("Response content preview:", response.text[:200])
            with open(f"error_{wp_id}.html", "w", encoding="utf-8") as f:
                f.write(response.text)
            print(f"Error HTML saved to: error_{wp_id}.html")
            raise Exception(f"Received non-image content: {response.headers.get('Content-Type')}")
        image_data = BytesIO(response.content)
        try:
            img = Image.open(image_data)
            img.verify()
            image_data.seek(0)
            return Image.open(image_data)
        except Exception as e:
            raise Exception(f"Failed to open image from BytesIO: {e}")
    else:
        raise Exception(f"Failed to download WikiPathways image: {response.status_code}")


def parse_gpml_to_entries(gpml_content):
    root = ET.fromstring(gpml_content)
    namespace = "{http://pathvisio.org/GPML/2013a}"

    entries = []
    groups = []
    relations = []

    # Debug: Print all top-level tags to verify structure
    print("Top-level tags in GPML:", [elem.tag for elem in root])

    # Parse DataNodes (similar to KEGG entries)
    datanodes = root.findall(f".//{namespace}DataNode")
    print(f"Found {len(datanodes)} DataNode elements")
    for datanode in datanodes:
        graph_id = datanode.get("GraphId", "Unknown")
        print(f"Processing DataNode with GraphId: {graph_id}")

        graphics = datanode.find(f"{namespace}Graphics")
        if graphics is None:
            print(f"Warning: DataNode {graph_id} is missing Graphics element, skipping")
            continue

        center_x = graphics.get("CenterX", 0)
        center_y = graphics.get("CenterY", 0)
        width = graphics.get("Width", 0)
        height = graphics.get("Height", 0)

        if not all([center_x, center_y, width, height]):
            print(
                f"Warning: DataNode {graph_id} is missing required Graphics attributes (CenterX={center_x}, CenterY={center_y}, Width={width}, Height={height}), skipping")
            continue

        entry = {
            "id": graph_id,
            "name": datanode.get("TextLabel", ""),
            "type": datanode.get("Type", "gene").lower(),
            "x": float(center_x) - float(width) / 2,
            "y": float(center_y) - float(height) / 2,
            "width": float(width),
            "height": float(height),
            "first_name": datanode.get("TextLabel", "").split(",")[0].strip(),
            "group_ref": datanode.get("GroupRef")  # Store GroupRef in the entry
        }
        print(f"Parsed DataNode {graph_id}: {entry}")

        # Add Xref information if needed
        xref = datanode.find(f"{namespace}Xref")
        if xref is not None:
            entry["xref"] = {
                "Database": xref.get("Database", ""),
                "ID": xref.get("ID", "")
            }

        entries.append(entry)

    # Parse Groups
    group_elements = root.findall(f".//{namespace}Group")
    print(f"Found {len(group_elements)} Group elements")
    for group in group_elements:
        group_id = group.get("GraphId", "")
        if not group_id:
            print("Warning: Group element missing GraphId, skipping")
            continue

        graphics = group.find(f"{namespace}Graphics")
        if graphics is None:
            print(f"Warning: Group {group_id} is missing Graphics element, skipping")
            continue

        group_data = {
            "id": group_id,
            "x": float(graphics.get("CenterX", 0)) - float(graphics.get("Width", 0)) / 2,
            "y": float(graphics.get("CenterY", 0)) - float(graphics.get("Height", 0)) / 2,
            "width": float(graphics.get("Width", 0)),
            "height": float(graphics.get("Height", 0)),
            "fgcolor": graphics.get("Color", "#000000"),
            "bgcolor": graphics.get("FillColor", "#FFFFFF"),
            "components": []  # In GPML, group membership is via GroupRef in DataNodes
        }
        groups.append(group_data)
        print(f"Added group {group_id}")

    # Link DataNodes to Groups via GroupRef
    for entry in entries:
        group_ref = entry.get("group_ref")
        if group_ref:
            for group in groups:
                if group["id"] == group_ref:
                    group["components"].append(entry["id"])
                    print(f"Linked DataNode {entry['id']} to group {group_ref}")

    # Parse Interactions (similar to KEGG relations)
    for interaction in root.findall(f".//{namespace}Interaction"):
        graphics = interaction.find(f"{namespace}Graphics")
        if graphics is not None:
            points = graphics.findall(f"{namespace}Point")
            if len(points) >= 2:
                start_point = points[0]
                end_point = points[-1]
                relation = {
                    "entry1": start_point.get("GraphRef"),
                    "entry2": end_point.get("GraphRef"),
                    "type": "interaction"
                }
                if relation["entry1"] and relation["entry2"]:
                    relations.append(relation)
                    print(f"Added relation: {relation}")

    print(f"Parsed {len(entries)} entries, {len(groups)} groups, and {len(relations)} relations")
    return entries, groups, relations


def visualize_wikipathways_pathway(wp_id, proteomic_data, phospho_data, settings):
    try:
        # Download GPML
        gpml_content = pwp.get_pathway(wp_id)
        output_file = savefile(wp_id, 'wikipathways', '', extension=".gpml", include_timestamp=False)
        with open(output_file, "w") as f:
            f.write(gpml_content)
        print(f"GPML file downloaded: {output_file}")

        # Skip image download
        pathway_image = None  # No background image
        print("No background image (skipped)")

        # Parse GPML
        entries, groups, relations = parse_gpml_to_entries(gpml_content)
        print(f"Parsed {len(entries)} entries, {len(groups)} groups, and {len(relations)} relations")

        # Update settings with WikiPathways-specific info
        settings['pathway_id'] = wp_id
        settings['show_background_image'] = False  # Ensure no image is expected

        # Create and display viewer
        viewer = PathwayViewer(entries, pathway_image, proteomic_data, phospho_data, settings)
        print("PathwayViewer created")

        # Store relations and groups for display_pathway to use
        viewer.relations = relations
        viewer.groups = groups
        viewer.display_pathway(entries, proteomic_data, phospho_data)
        print("Pathway displayed")

        return viewer
    except Exception as e:
        print(f"Error in visualize_wikipathways_pathway: {e}")
        return None


class ArrowItem(QGraphicsPolygonItem):
    def __init__(self, start_point, end_point, parent=None):
        super().__init__(parent)
        self.start_point = start_point
        self.end_point = end_point

        # Draw the line
        self.line = QGraphicsLineItem(start_point.x(), start_point.y(), end_point.x(), end_point.y())
        self.line.setPen(QPen(Qt.black, 1))

        # Draw the arrowhead
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

    def add_to_scene(self, scene):
        scene.addItem(self.line)
        scene.addItem(self)

class PhosphoCircleItem(QGraphicsEllipseItem):
    def __init__(self, x, y, width, height, phospho_row, viewer, parent=None):
        super().__init__(x, y, width, height, parent)
        self.phospho_row = phospho_row
        self.viewer = viewer
        self.setAcceptHoverEvents(True)
        self.tooltip = None
        self.bg_rect = None  # Explicitly track bg_rect separately

    def hoverEnterEvent(self, event):
        try:
            if self.phospho_row[self.viewer.modulation_column] != '+':
                return

            # If a tooltip already exists, don’t recreate it
            if self.tooltip is not None:
                return

            annotation_columns = self.viewer.phos_annotation_columns
            annotations = [str(self.phospho_row[col]) for col in annotation_columns
                           if col in self.phospho_row.index]
            tooltip_text = "\n\n".join(annotations)

            # Create tooltip
            self.tooltip = QGraphicsTextItem(tooltip_text)
            self.tooltip.setFont(QFont("Arial", 8))
            self.tooltip.setDefaultTextColor(QColor(0, 0, 0))
            max_width = 300
            self.tooltip.setTextWidth(max_width)

            tooltip_pos = self.mapToScene(self.rect().center() + QPointF(self.rect().width() + 5, 0))
            self.tooltip.setPos(tooltip_pos)

            # Create background rectangle
            bounding_rect = self.tooltip.boundingRect()
            self.bg_rect = QGraphicsRectItem(bounding_rect.adjusted(-2, -2, 2, 2))
            self.bg_rect.setBrush(QColor(255, 255, 200, 200))
            self.bg_rect.setPen(QColor(0, 0, 0))
            self.bg_rect.setPos(tooltip_pos)
            self.bg_rect.setZValue(10)

            # Add items to scene safely
            if self.scene():
                self.scene().addItem(self.bg_rect)
                self.scene().addItem(self.tooltip)
                self.tooltip.setZValue(11)
        except Exception as e:
            print(f"Error in hoverEnterEvent: {e}")
            self.tooltip = None  # Reset on error
            self.bg_rect = None

    def hoverLeaveEvent(self, event):
        try:
            if self.tooltip is not None and self.scene():
                # Only remove items if they’re still in the scene
                if self.bg_rect in self.scene().items():
                    self.scene().removeItem(self.bg_rect)
                if self.tooltip in self.scene().items():
                    self.scene().removeItem(self.tooltip)
            self.tooltip = None
            self.bg_rect = None
        except Exception as e:
            print(f"Error in hoverLeaveEvent: {e}")
            self.tooltip = None  # Reset on error
            self.bg_rect = None

class ProteinBoxItem(QGraphicsRectItem):
    def __init__(self, x, y, width, height, entry_id, viewer):
        super().__init__(x, y, width, height)
        self.entry_id = entry_id
        self.viewer = viewer
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setAcceptHoverEvents(True)

    def contextMenuEvent(self, event):
        menu = QMenu()
        entry_data = self.viewer.protein_data_map[self.entry_id]
        current_prot = self.viewer.current_proteins.get(self.entry_id)
        for protein in entry_data['proteins']:
            if protein in self.viewer.proteomic_data[self.viewer.hsa_id_column].values:
                action = menu.addAction(protein)
                action.setCheckable(True)
                action.setChecked(protein == current_prot)
                action.triggered.connect(lambda checked, p=protein:
                                         self.viewer.switch_protein(self.entry_id, p))
        menu.exec_(event.screenPos())

class PathwayViewer(QGraphicsView):
    def __init__(self, entries, image, proteomic_data, phospho_data, settings):
        super().__init__()
        self.entries = entries  # Store entries for scene sizing
        self.scene = QGraphicsScene(self)
        self.setScene(self.scene)

        # Estimate scene size based on entries
        if entries:
            max_x = max(entry["x"] + entry["width"] for entry in entries)
            max_y = max(entry["y"] + entry["height"] for entry in entries)
            min_x = min(entry["x"] for entry in entries)
            min_y = min(entry["y"] for entry in entries)
            self.scene.setSceneRect(min_x - 50, min_y - 50, (max_x - min_x) + 100, (max_y - min_y) + 100)
        else:
            self.scene.setSceneRect(0, 0, 1000, 1000)

        # Ensure the view fits the scene
        self.setRenderHint(QPainter.Antialiasing)
        self.setRenderHint(QPainter.SmoothPixmapTransform)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorUnderMouse)
        self.scale(0.5, 0.5)  # Start zoomed out to see the whole pathway

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
        self.phospho_data = phospho_data
        self.protein_data_map = {}
        self.current_proteins = {}
        self.phospho_items = {}
        self.settings = settings

        self.protein_selection_option = settings.get('protein_selection_option', 2)
        self.phos_selection_option = settings.get('phos_selection_option', 2)
        self.fold_change_column = settings.get('fold_change_column', 'ER+_Est-_x-y_TNBC')
        self.phos_fold_change_column = settings.get('phos_fold_change_column', 'ER+_Est-_x-y_TNBC')
        self.hsa_id_column = settings.get('hsa_id_column', 'KEGG_hsa')
        self.gene_name_column = settings.get('gene_name_column', 'Gene Symbol')
        self.modulation_column = settings.get('modulation_column', 'C: Regulatory site')
        self.phos_site_column = settings.get('phos_site_column', 'T: Site Position')
        self.phos_max_display = settings.get('phos_max_display', 4)
        self.negative_color = settings.get('negative_color', (255, 0, 0))
        self.positive_color = settings.get('positive_color', (0, 0, 255))
        self.max_negative = settings.get('max_negative', -2)
        self.max_positive = settings.get('max_positive', 2)
        self.prot_label_font = settings.get('prot_label_font', 'Arial')
        self.prot_label_size = settings.get('prot_label_size', 6)
        self.phos_label_font = settings.get('phos_label_font', 'Arial')
        self.phos_label_color = settings.get('phos_label_color', QColor(0, 0, 0))
        self.phos_label_size = settings.get('phos_label_size', 4)
        self.phos_circle_radius = settings.get('phos_circle_radius', 5)
        self.phos_circle_spacing = settings.get('phos_circle_spacing', 2)
        self.show_background_image = settings.get('show_background_image', True)
        self.display_types = settings.get('display_types', ['gene'])
        self.show_groups = settings.get('show_groups', False)
        self.phos_annotation_columns = settings['phos_annotation_columns']
        self.settings = settings  # Store settings for use in methods

        # Background image setup
        if self.show_background_image and image:
            image = image.convert("RGBA")
            data = image.tobytes("raw", "RGBA")
            qimage = QImage(data, image.size[0], image.size[1], QImage.Format_RGBA8888)
            qpixmap = QPixmap.fromImage(qimage)
            self.background_item = QGraphicsPixmapItem(qpixmap)
            self.background_item.setPos(22, 9)
            self.scene.addItem(self.background_item)
            # Set scene rect to exactly match the image size
            self.scene.setSceneRect(0, 0, image.size[0], image.size[1])
            # Store image dimensions for later use
            self.image_width = image.size[0]
            self.image_height = image.size[1]
        else:
            self.image_width = 0
            self.image_height = 0

        self.protein_group = QGraphicsItemGroup()
        self.scene.addItem(self.protein_group)

        # Add Save to PDF button
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
        self.save_button.clicked.connect(self.save_to_pdf)
        self.update_button_position()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.update_button_position()
        # Automatically zoom to fit the pathway image when maximized or fullscreen
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
        """Export the scene as a PDF with minimal whitespace, fitting the pathway content."""
        try:
            # Open save file dialog
            pathway_id = self.settings.get('pathway_id', 'hsa04010')
            default_filename = f"{pathway_id}_pathway_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
            output_file, _ = QFileDialog.getSaveFileName(
                self,
                "Save Pathway as PDF",
                default_filename,
                "PDF Files (*.pdf)"
            )
            if not output_file:
                print("Save operation canceled by user.")
                return

            # Create PDF writer with high resolution
            pdf_writer = QPdfWriter(output_file)
            pdf_writer.setResolution(300)  # 300 DPI for high quality

            # Get all scene items and compute content bounds manually, excluding non-graphics items
            all_items = self.scene.items()
            if not all_items:
                print("No content to save.")
                return

            # Initialize content rect with the first graphics item's bounding rect
            content_rect = None
            for item in all_items:
                # Skip non-QGraphicsItem objects (e.g., QPushButton like save_button)
                if not isinstance(item, QGraphicsItem):
                    continue
                item_rect = item.sceneBoundingRect()
                if content_rect is None:
                    content_rect = item_rect
                else:
                    content_rect = content_rect.united(item_rect)

            if content_rect is None:
                print("No valid graphics items found in scene.")
                return

            # Adjust content_rect to account for background image offset (22, 9)
            if self.show_background_image and self.background_item:
                bg_pos = self.background_item.pos()  # Should be (22, 9)
                content_rect.adjust(-bg_pos.x(), -bg_pos.y(), 0, 0)  # Shift left/top to include offset
                # Ensure right/bottom extend to full image size if background defines bounds
                bg_rect = self.background_item.sceneBoundingRect()
                content_rect = content_rect.united(bg_rect)

            content_width = content_rect.width()
            content_height = content_rect.height()
            print(f"Content dimensions: {content_width} x {content_height} pixels")

            # Define minimal margins (e.g., 5 mm on each side)
            margin_mm = 5
            mm_per_pixel = 25.4 / 300  # mm per pixel at 300 DPI
            content_width_mm = content_width * mm_per_pixel
            content_height_mm = content_height * mm_per_pixel
            pdf_width_mm = content_width_mm + 2 * margin_mm
            pdf_height_mm = content_height_mm + 2 * margin_mm
            print(f"PDF page size: {pdf_width_mm:.2f} x {pdf_height_mm:.2f} mm")

            # Set custom page size
            custom_size = QPageSize(QSizeF(pdf_width_mm, pdf_height_mm), QPageSize.Millimeter, "Custom")
            pdf_writer.setPageSize(custom_size)
            pdf_writer.setPageMargins(QMarginsF(0, 0, 0, 0))  # No additional margins beyond our calculation

            # Create painter and set up rendering
            painter = QPainter(pdf_writer)
            painter.setRenderHint(QPainter.Antialiasing)

            # Define target rect to fit content within the page, accounting for margins
            target_width = (pdf_width_mm - 2 * margin_mm) / mm_per_pixel  # Back to pixels
            target_height = (pdf_height_mm - 2 * margin_mm) / mm_per_pixel
            target_rect = QRectF(margin_mm / mm_per_pixel, margin_mm / mm_per_pixel, target_width, target_height)
            print(
                f"Target rect: {target_width} x {target_height} at ({margin_mm / mm_per_pixel}, {margin_mm / mm_per_pixel})")

            # Adjust content_rect to start at (0, 0) for rendering
            adjusted_content_rect = QRectF(0, 0, content_width, content_height)
            self.scene.render(painter, target_rect, adjusted_content_rect)
            painter.end()

            print(f"Pathway saved as PDF: {output_file} (Size: {pdf_width_mm:.2f} x {pdf_height_mm:.2f} mm)")
        except Exception as e:
            print(f"Error saving PDF: {e}")

    def get_color(self, fold_change):
        """Get color based on fold change, interpolating from white (0) to negative/positive colors."""
        if pd.isna(fold_change):
            return QColor(128, 128, 128)  # Grey for missing data

        # White as the base color (fold change = 0)
        white = (255, 255, 255)
        r0, g0, b0 = white  # Starting point

        if fold_change < 0:
            # Interpolate from white to negative_color (e.g., red)
            r1, g1, b1 = self.negative_color
            t = min(abs(fold_change) / abs(self.max_negative), 1)  # Normalize between 0 and 1
        else:
            # Interpolate from white to positive_color (e.g., blue)
            r1, g1, b1 = self.positive_color
            t = min(fold_change / self.max_positive, 1)  # Normalize between 0 and 1

        # Linear interpolation: color = (1 - t) * white + t * target_color
        r = int((1 - t) * r0 + t * r1)
        g = int((1 - t) * g0 + t * g1)
        b = int((1 - t) * b0 + t * b1)

        return QColor(r, g, b)

    def choose_protein(self, proteins):
        # Debug: Inspect proteomic data
        print(f"Proteomic data columns: {list(self.proteomic_data.columns)}")
        print(f"Sample KEGG_hsa values: {self.proteomic_data['KEGG_hsa'].head().tolist()}")

        # Get Xref ID (Entrez Gene) from the current entry
        xref_id = self.protein_data_map.get(self.current_entry_id, {}).get('xref_id', None)

        # Match using KEGG_hsa (strip 'hsa:' to match Entrez Gene ID)
        if xref_id:
            hsa_ids = self.proteomic_data.get('KEGG_hsa', []).str.replace('hsa:', '')
            print(f"Looking for Xref ID {xref_id} in KEGG_hsa (after stripping 'hsa:')")
            if xref_id in hsa_ids.values:
                # Find the corresponding protein (e.g., Gene Symbol for display)
                idx = hsa_ids[hsa_ids == xref_id].index[0]
                protein = self.proteomic_data.loc[idx, 'Gene Symbol']
                print(f"Matched Xref {xref_id} to KEGG_hsa, protein: {protein}")
                return protein
            else:
                print(f"Xref ID {xref_id} not found in KEGG_hsa")
        else:
            print(f"No Xref ID for entry {self.current_entry_id}")

        print(f"No match for entry {self.current_entry_id}, Xref {xref_id}")
        return None

    def prioritize_phospho_sites(self, protein):
        protein_phospho = self.phospho_data[self.phospho_data['Gene Symbol'] == protein]
        if self.phos_selection_option == 1:
            protein_phospho = protein_phospho.sort_values(by=self.phos_fold_change_column, key=abs, ascending=False)
        elif self.phos_selection_option == 2:
            protein_phospho = protein_phospho.sort_values(
                by=[self.modulation_column, self.phos_fold_change_column],
                key=lambda x: x == '+', ascending=False)
        elif self.phos_selection_option == 3:
            protein_phospho = protein_phospho[protein_phospho[self.modulation_column] == '+']
            protein_phospho = protein_phospho.sort_values(by=self.phos_fold_change_column, key=abs, ascending=False)
        protein_phospho['is_modulating'] = protein_phospho[self.modulation_column] == '+'
        return protein_phospho.head(self.phos_max_display)

    def display_pathway(self, entries, proteomic_data, phospho_data):
        try:
            print(f"Entry types in GPML: {set(entry['type'] for entry in entries)}")

            # Relations (arrows)
            if hasattr(self, 'relations'):
                display_types_set = set(self.display_types) if not isinstance(self.display_types, str) else {
                    self.display_types}
                for relation in self.relations:
                    entry1_data = find_entry_or_group(relation["entry1"], entries, self.groups)
                    entry2_data = find_entry_or_group(relation["entry2"], entries, self.groups)
                    if entry1_data and entry2_data and entry1_data["type"] in display_types_set and entry2_data[
                        "type"] in display_types_set:
                        x1 = entry1_data["x"] + entry1_data["width"] / 2
                        y1 = entry1_data["y"] + entry1_data["height"] / 2
                        x2 = entry2_data["x"] + entry2_data["width"] / 2
                        y2 = entry2_data["y"] + entry2_data["height"] / 2
                        angle = calculate_angle(x1, y1, x2, y2)
                        side = determine_arrow_side(angle)
                        x2_adj, y2_adj = adjust_arrow_endpoint(entry2_data["x"], entry2_data["y"],
                                                               entry2_data["width"], entry2_data["height"], side)
                        arrow = ArrowItem(QPointF(x1, y1), QPointF(x2_adj, y2_adj))
                        arrow.add_to_scene(self.scene)
                        print(f"Added arrow from {entry1_data['id']} to {entry2_data['id']}")

            # Entries with data matching
            chosen_proteins = set()
            for entry in entries:
                if entry["type"].lower() in ["geneproduct", "gene", "protein"]:
                    x = entry["x"]
                    y = entry["y"]
                    width = entry["width"]
                    height = entry["height"]
                    rect = QGraphicsRectItem(x, y, width, height)
                    rect.setBrush(QBrush(QColor("lightgray")))
                    rect.setPen(QPen(QColor("black")))
                    self.scene.addItem(rect)
                    label = QGraphicsTextItem(entry["name"])
                    label.setPos(x, y - 20)  # Move label above the box for clarity
                    self.scene.addItem(label)
                    print(f"Added box for {entry['id']} - {entry['name']} at ({x}, {y})")

                    # Store for protein matching
                    self.protein_data_map[entry["id"]] = {
                        'proteins': entry["name"].split(),
                        'xref_id': entry.get("xref", {}).get("ID", None),
                        'x': x, 'y': y, 'width': width, 'height': height,
                        'first_name': entry["name"]  # Ensure first_name is set
                    }
                    self.current_entry_id = entry["id"]
                    chosen_protein = self.choose_protein(entry["name"].split())
                    if chosen_protein:
                        chosen_proteins.add(chosen_protein)
                        self.current_proteins[entry["id"]] = chosen_protein
                        self.draw_protein_box(entry["id"], chosen_protein)

            # Save visualized proteins
            temp_file = savefile(self.settings.get('pathway_id', 'WP382'),
                                 self.settings.get('output_subdir', 'output/testing_file_001'),
                                 '_visprots', extension='.txt', include_timestamp=False)
            with open(temp_file, 'w') as f:
                for protein in chosen_proteins:
                    f.write(f"{protein}\n")

            self.save_button.raise_()
            self.save_button.show()
            self.fitInView(self.scene.sceneRect(), Qt.KeepAspectRatio)

        except Exception as e:
            print(f"Error in display_pathway: {e}")

    def draw_protein_box(self, entry_id, protein):
        if entry_id in self.phospho_items:
            for item in self.phospho_items[entry_id]:
                self.scene.removeItem(item)
            del self.phospho_items[entry_id]

        entry_data = self.protein_data_map[entry_id]
        x = entry_data['x']
        y = entry_data['y']
        width = entry_data['width']
        height = entry_data['height']

        box = ProteinBoxItem(x, y, width, height, entry_id, self)
        protein_in_data = protein in self.proteomic_data[self.hsa_id_column].values

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

        label = QGraphicsTextItem(gene_name)
        label.setFont(QFont(self.prot_label_font, self.prot_label_size))
        label.setDefaultTextColor(QColor(0, 0, 0))
        bounding_rect = label.boundingRect()
        adjusted_x = x + (width / 2) - (bounding_rect.width() / 2)
        adjusted_y = y + (height / 2) - (bounding_rect.height() / 2)
        label.setPos(adjusted_x, adjusted_y)
        label.setParentItem(box)

        if protein_in_data and not self.phospho_data[self.phospho_data[self.hsa_id_column] == protein].empty:
            phospho_sites = self.prioritize_phospho_sites(protein)
            if not phospho_sites.empty:
                self.phospho_items[entry_id] = []
                positions = {
                    'N1': (x + width * 0.2, y - self.phos_circle_spacing),
                    'N2': (x + width * 0.5, y - self.phos_circle_spacing),
                    'N3': (x + width * 0.8, y - self.phos_circle_spacing),
                    'S1': (x + width * 0.2, y + height + self.phos_circle_spacing),
                    'S2': (x + width * 0.5, y + height + self.phos_circle_spacing),
                    'S3': (x + width * 0.8, y + height + self.phos_circle_spacing),
                    'W1': (x - self.phos_circle_spacing, y + height * 0.33),
                    'W2': (x - self.phos_circle_spacing, y + height * 0.66),
                    'E1': (x + width + self.phos_circle_spacing, y + height * 0.33),
                    'E2': (x + width + self.phos_circle_spacing, y + height * 0.66)
                }
                position_priority = ['N1', 'N2', 'N3', 'S1', 'S2', 'S3', 'W1', 'W2', 'E1', 'E2']

                exceptions = self.settings.get('parsed_exceptions', {})
                proximities = self.settings.get('proximities', {})
                protein_to_entry_ids = self.settings.get('protein_to_entry_ids', {})
                pathway_id = self.settings.get('pathway_id', 'hsa04010')
                hsa_exceptions = exceptions.get('hsa_specific', {}).get(pathway_id, {'proximity': [], 'specific': {}})
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

                # Apply specific rules
                for rule in specific_rules:
                    action = rule[0].lower()
                    if action == 'block' and len(rule) > 1:
                        blocked_positions.add(rule[1])
                        print(f"Specific rule: Blocking {rule[1]} for {protein}")
                    elif action == 'priority':
                        custom_priority = rule[1:]
                        print(f"Specific rule: Setting priority {custom_priority} for {protein}")
                    elif action == 'move_circle' and len(rule) == 4:
                        pos, dx, dy = rule[1], float(rule[2]), float(rule[3])
                        move_circle[pos] = (dx, dy)
                        print(f"Specific rule: Moving circle {pos} by ({dx}, {dy}) for {protein}")
                    elif action == 'move_label' and len(rule) == 5:
                        pos, dx, dy, align = rule[1], float(rule[2]), float(rule[3]), rule[4]
                        move_label[pos] = (dx, dy, align)
                        print(f"Specific rule: Moving label {pos} by ({dx}, {dy}, {align}) for {protein}")

                # Apply proximity rules
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
                        print(f"Rule passed for {entry_id} ({protein}): {rule}")
                        passed_rules_by_action[rule['action']].append(rule)

                # Apply block and priority (one per protein box)
                for action in ['block', 'priority']:
                    if passed_rules_by_action[action]:
                        highest_rule = max(passed_rules_by_action[action], key=lambda r: r['priority'])
                        if action == 'block':
                            for pos in highest_rule['values']:
                                blocked_positions.add(pos)
                            print(f"Proximity rule: Blocking {highest_rule['values']} for {protein}")
                        elif action == 'priority':
                            custom_priority = highest_rule['values']
                            print(f"Proximity rule: Setting priority {custom_priority} for {protein}")

                # Apply move_circle and move_label (one per phos circle position)
                for action in ['move_circle', 'move_label']:
                    for rule in sorted(passed_rules_by_action[action],
                                       key=lambda r: (r['priority'], proximity_rules.index(r))):
                        if action == 'move_circle' and len(rule['values']) == 3:
                            pos, dx, dy = rule['values'][0], float(rule['values'][1]), float(rule['values'][2])
                            move_circle[pos] = (dx, dy)
                            print(f"Proximity rule: Moving circle {pos} by ({dx}, {dy}) for {protein}")
                        elif action == 'move_label' and len(rule['values']) == 4:
                            pos, dx, dy, align = rule['values'][0], float(rule['values'][1]), float(rule['values'][2]), \
                            rule['values'][3]
                            move_label[pos] = (dx, dy, align)
                            print(f"Proximity rule: Moving label {pos} by ({dx}, {dy}, {align}) for {protein}")

                available_positions = [p for p in (custom_priority or position_priority) if p not in blocked_positions]
                print(f"Available positions for {entry_id} ({protein}): {available_positions}")

                for i, (_, row) in enumerate(phospho_sites.iterrows()):
                    if i >= len(available_positions):
                        break
                    pos_key = available_positions[i]
                    circle_x, circle_y = positions[pos_key]
                    if pos_key in move_circle:
                        circle_x += move_circle[pos_key][0]
                        circle_y += move_circle[pos_key][1]
                        print(
                            f"Applying move_circle to {pos_key}: ({move_circle[pos_key][0]}, {move_circle[pos_key][1]})")

                    phos_fold_change = row[self.phos_fold_change_column]
                    phos_color = self.get_color(phos_fold_change)
                    circle = PhosphoCircleItem(
                        circle_x - self.phos_circle_radius,
                        circle_y - self.phos_circle_radius,
                        self.phos_circle_radius * 2,
                        self.phos_circle_radius * 2,
                        row,
                        self,
                        parent=box
                    )
                    circle.setBrush(phos_color)
                    circle.setPen(QColor(0, 0, 0))
                    self.phospho_items[entry_id].append(circle)
                    self.scene.addItem(circle)

                    if row[self.modulation_column] == '+':
                        plus_sign = QGraphicsTextItem("+", parent=box)
                        plus_sign.setFont(QFont(self.phos_label_font, self.phos_label_size + 2, QFont.Bold))
                        plus_sign.setDefaultTextColor(QColor(0, 0, 0))
                        plus_sign.setAcceptHoverEvents(False)
                        plus_sign.setZValue(12)
                        text_rect = plus_sign.boundingRect()
                        adjusted_x = circle_x - text_rect.width() / 2
                        adjusted_y = circle_y - text_rect.height() / 2
                        plus_sign.setPos(adjusted_x, adjusted_y)
                        self.phospho_items[entry_id].append(plus_sign)
                        self.scene.addItem(plus_sign)

                    site_number = row[self.phos_site_column]
                    if pos_key in move_label:
                        x_offset, y_offset, ha = move_label[pos_key]
                        print(f"Applying move_label to {pos_key}: ({x_offset}, {y_offset}, {ha})")
                    else:
                        defaults = {
                            'N1': (-3, -5, 'right'), 'N2': (0, -11, 'center'), 'N3': (3, -5, 'left'),
                            'S1': (-3, 5, 'right'), 'S2': (0, 12, 'center'), 'S3': (3, 5, 'left'),
                            'W1': (-3, -2, 'right'), 'W2': (-3, 2, 'right'),
                            'E1': (3, -2, 'left'), 'E2': (3, 2, 'left')
                        }
                        x_offset, y_offset, ha = defaults[pos_key]

                    text = QGraphicsTextItem(str(site_number), parent=box)
                    text.setFont(QFont(self.phos_label_font, self.phos_label_size))
                    text.setDefaultTextColor(self.phos_label_color)
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
                    self.phospho_items[entry_id].append(text)
                    self.scene.addItem(text)

    def switch_protein(self, entry_id, new_protein):
        self.current_proteins[entry_id] = new_protein
        self.draw_protein_box(entry_id, new_protein)

    def wheelEvent(self, event):
        # Get current scale and content dimensions
        current_scale = self.transform().m11()  # Horizontal scale factor (assumes uniform scaling)
        content_rect = self.scene.itemsBoundingRect()
        viewport_width = self.viewport().width()
        viewport_height = self.viewport().height()
        content_width = content_rect.width() * current_scale
        content_height = content_rect.height() * current_scale

        # Zoom handling (only with Ctrl modifier)
        if event.modifiers() & Qt.ControlModifier:  # Check if Ctrl is held
            zoom_in_factor = 1.1
            zoom_out_factor = 1 / zoom_in_factor
            if event.angleDelta().y() > 0:
                self.scale(zoom_in_factor, zoom_in_factor)
            elif event.angleDelta().y() < 0:
                new_scale = current_scale * zoom_out_factor
                if self.show_background_image and new_scale < 1.0:  # Prevent zooming out beyond fit
                    self.fitInView(self.background_item, Qt.KeepAspectRatio)
                else:
                    self.scale(zoom_out_factor, zoom_out_factor)
            event.accept()
            return

        # Scroll handling when zoomed in
        scrolled = False

        # Horizontal scroll
        if event.angleDelta().x() != 0 and content_width > viewport_width:
            scroll_step = -event.angleDelta().x() * 0.5  # Adjust sensitivity
            horizontal_bar = self.horizontalScrollBar()
            current_value = horizontal_bar.value()
            new_value = current_value + int(scroll_step)
            new_value = max(horizontal_bar.minimum(), min(new_value, horizontal_bar.maximum()))
            horizontal_bar.setValue(new_value)
            scrolled = True

        # Vertical scroll
        if event.angleDelta().y() != 0 and content_height > viewport_height:
            scroll_step = -event.angleDelta().y() * 0.5  # Adjust sensitivity (negative for natural scrolling)
            vertical_bar = self.verticalScrollBar()
            current_value = vertical_bar.value()
            new_value = current_value + int(scroll_step)
            new_value = max(vertical_bar.minimum(), min(new_value, vertical_bar.maximum()))
            vertical_bar.setValue(new_value)
            scrolled = True

        if scrolled:
            event.accept()  # Mark event as handled if scrolling occurred

def visualize_kegg_pathway(pathway_id, proteomic_data, phospho_data, settings):
    try:
        kgml_file = download_kegg_kgml(pathway_id)
        print(f"KGML file downloaded: {kgml_file}")
        pathway_image = download_kegg_image(pathway_id) if settings.get('show_background_image', True) else None
        print("Pathway image downloaded" if pathway_image else "No background image")
        entries, groups = parse_kgml(kgml_file)
        print(f"Parsed {len(entries)} entries and {len(groups)} groups")
        viewer = PathwayViewer(entries, pathway_image, proteomic_data, phospho_data, settings)
        print("PathwayViewer created")
        viewer.display_pathway(entries, proteomic_data, phospho_data)
        print("Pathway displayed")
        return viewer  # Return instead of showing
    except Exception as e:
        print(f"Error: {e}")
        return None


if __name__ == "__main__":
    settings = {
        'pathway_id': 'WP382',  # Example WikiPathways ID
        'prot_inputfile': r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\Phosmap\BCp2_ProtMaphsaanno.txt",
        'phos_inputfile': r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\Phosmap\BCp2_PhosMaphsanno.txt",
        'phos_annotation_columns': ['C: Regulatory site function', 'C: Regulatory site process'],
        'display_types': ['geneproduct'],  # Match WikiPathways type
        'fold_change_column': 'ER+_Est-_x-y_ER+',  # Match your data column
        'phos_fold_change_column': 'ER+_Est-_x-y_ER+',  # Match your data column
        'hsa_id_column': 'KEGG_hsa',  # Match your data column
        'gene_name_column': 'Gene Symbol',  # Match your data column
    }
    proteomic_data = pd.read_csv(settings['prot_inputfile'], sep="\t")
    phospho_data = pd.read_csv(settings['phos_inputfile'], sep="\t")
    app = QApplication(sys.argv)

    # Use WikiPathways visualization instead of KEGG
    viewer = visualize_wikipathways_pathway(settings['pathway_id'], proteomic_data, phospho_data, settings)
    if viewer:
        viewer.show()
        sys.exit(app.exec_())