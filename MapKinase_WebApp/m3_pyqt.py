import sys
import json
import os
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QToolTip
from PyQt5.QtGui import QPainter, QPen, QBrush, QColor, QFont, QFontMetrics
from PyQt5.QtCore import Qt, QRectF, QPoint
from PyQt5.QtCore import QRect as QRectImported
QRect = QRectImported

# Path to the JSON file
json_file_path = r"C:\Users\clayt\OneDrive - Brigham Young University\pycharm\ttk_projects\Phosmap\scripts\output\testing_file_001\hsa04010_pathway_data_20250723_180525.json"

class PathwayWidget(QWidget):
    def __init__(self, json_data, parent=None):
        super().__init__(parent)
        self.json_data = json_data
        self.settings = json_data.get('general_data', {}).get('settings', {})
        self.protein_data = json_data.get('protein_data', {})
        self.protbox_data = json_data.get('protbox_data', [])
        self.groups = json_data.get('groups', [])
        self.arrows = json_data.get('arrows', [])
        self.setMouseTracking(True)
        self.tooltip = None
        print("Initializing PathwayWidget")
        print(f"Settings: {self.settings}")
        print(f"Protein data keys: {list(self.protein_data.keys())}")
        print(f"Protbox count: {len(self.protbox_data)}")
        print(f"Group count: {len(self.groups)}")
        print(f"Arrow count: {len(self.arrows)}")
        self.init_ui()

    def init_ui(self):
        try:
            if not self.protbox_data:
                print("Warning: No protbox_data found. Setting default window size.")
                self.setMinimumSize(800, 600)
                return
            max_x = max(pb.get('x', 0) + pb.get('width', 0) for pb in self.protbox_data if pb.get('x') is not None and pb.get('width') is not None) + 50
            max_y = max(pb.get('y', 0) + pb.get('height', 0) for pb in self.protbox_data if pb.get('y') is not None and pb.get('height') is not None) + 50
            max_x = max(max_x, 800)
            max_y = max(max_y, 600)
            self.setFixedSize(int(max_x), int(max_y))
            print(f"Set window size to {int(max_x)}x{int(max_y)}")
        except Exception as e:
            print(f"Error in init_ui: {e}")
            self.setMinimumSize(800, 600)

    def paintEvent(self, event):
        try:
            painter = QPainter(self)
            painter.setRenderHint(QPainter.Antialiasing)
            self.draw_protein_boxes(painter)
            self.draw_ptms(painter)
            if self.settings.get('show_groups', False):
                self.draw_groups(painter)
            if 'arrows' in self.settings.get('display_types', ['prot_box']):
                self.draw_arrows(painter)
        except Exception as e:
            print(f"Error in paintEvent: {e}")

    def draw_protein_boxes(self, painter):
        try:
            for protbox in self.protbox_data:
                protbox_id = protbox.get('protbox_id', '')
                x = protbox.get('x', 0)
                y = protbox.get('y', 0)
                width = protbox.get('width', 46)
                height = protbox.get('height', 17)
                if not all(isinstance(v, (int, float)) for v in [x, y, width, height]):
                    print(f"Warning: Invalid coordinates for protbox {protbox_id}: x={x}, y={y}, width={width}, height={height}")
                    continue
                uniprot_ids = protbox.get('proteins', [])
                label = protbox.get('backup_label', 'Unknown')
                fc_color = [128, 128, 128]
                annotations = ''
                selected_uniprot_id = None
                for uniprot_id in uniprot_ids:
                    if uniprot_id in self.protein_data:
                        protein = self.protein_data[uniprot_id]
                        label = protein.get('label', label)
                        fc_color = protein.get('fc_color_1', fc_color)
                        annotations = protein.get('annotations', '')
                        selected_uniprot_id = uniprot_id
                        break
                if not all(isinstance(c, int) and 0 <= c <= 255 for c in fc_color):
                    print(f"Warning: Invalid fc_color for protbox {protbox_id}: {fc_color}")
                    fc_color = [128, 128, 128]
                painter.setPen(QPen(QColor(0, 0, 0), 1))
                painter.setBrush(QBrush(QColor(*fc_color)))
                painter.drawRect(QRect(int(x), int(y), int(width), int(height)))
                font = QFont(self.settings.get('prot_label_font', 'Arial'), self.settings.get('prot_label_size', 6))
                painter.setFont(font)
                label_color = self.protein_data.get(selected_uniprot_id, {}).get('label_color', [0, 0, 0]) if selected_uniprot_id else [0, 0, 0]
                if not all(isinstance(c, int) and 0 <= c <= 255 for c in label_color):
                    print(f"Warning: Invalid label_color for protbox {protbox_id}: {label_color}")
                    label_color = [0, 0, 0]
                painter.setPen(QPen(QColor(*label_color)))
                fm = QFontMetrics(font)
                text_width = fm.width(label)
                text_height = fm.height()
                text_x = x + (width - text_width) / 2
                text_y = y + (height + text_height) / 2 - fm.descent()
                painter.drawText(int(text_x), int(text_y), label)
                protbox['tooltip'] = annotations
                protbox['rect'] = QRect(int(x), int(y), int(width), int(height))
                print(f"Drew protein box {protbox_id} at ({x}, {y}) with label {label}")
        except Exception as e:
            print(f"Error in draw_protein_boxes for protbox {protbox_id}: {e}")

    def draw_ptms(self, painter):
        try:
            for uniprot_id, protein in self.protein_data.items():
                ptms = protein.get('PTMs', {})
                print(f"Processing PTMs for protein {uniprot_id}: {list(ptms.keys())}")
                for ptm_key, ptm in ptms.items():
                    x = ptm.get('label_x', 0)
                    y = ptm.get('label_y', 0)
                    if not all(isinstance(v, (int, float)) for v in [x, y]):
                        print(f"Warning: Invalid coordinates for PTM {ptm_key}: x={x}, y={y}")
                        continue
                    shape = ptm.get('shape', 'circle').lower()
                    fc_color = ptm.get('fc_color_1', [128, 128, 128])
                    if not all(isinstance(c, int) and 0 <= c <= 255 for c in fc_color):
                        print(f"Warning: Invalid fc_color for PTM {ptm_key}: {fc_color}")
                        fc_color = [128, 128, 128]
                    label = ptm.get('label', '')
                    label_color = ptm.get('label_color', [0, 0, 0])
                    if not all(isinstance(c, int) and 0 <= c <= 255 for c in label_color):
                        print(f"Warning: Invalid label_color for PTM {ptm_key}: {label_color}")
                        label_color = [0, 0, 0]
                    label_centering = ptm.get('label_centering', 'center').lower()
                    symbol_type = ptm.get('symbol_type', '')
                    ptm_type = ptm.get('ptm_type', 'ptm_0')
                    tooltip = ptm.get('tooltip', '')
                    if not ptm_type or 'ptm_' not in ptm_type:
                        print(f"Warning: Invalid or missing ptm_type for PTM {ptm_key}: {ptm_type}")
                        continue
                    radius = self.settings.get('ptm_circle_radius', 5)
                    painter.setPen(QPen(QColor(0, 0, 0), 1))
                    painter.setBrush(QBrush(QColor(*fc_color)))
                    if shape == 'circle':
                        painter.drawEllipse(QRectF(x - radius, y - radius, radius * 2, radius * 2))
                    elif shape == 'square':
                        painter.drawRect(QRectF(x - radius, y - radius, radius * 2, radius * 2))
                    font = QFont(self.settings.get('ptm_label_font', 'Arial'), self.settings.get('ptm_label_size', 4))
                    painter.setFont(font)
                    painter.setPen(QPen(QColor(*label_color)))
                    fm = QFontMetrics(font)
                    text_width = fm.width(label)
                    text_height = fm.height()
                    if label_centering == 'left':
                        text_x = x + radius
                    elif label_centering == 'right':
                        text_x = x - radius - text_width
                    else:  # center
                        text_x = x - text_width / 2
                    text_y = y + text_height / 2 - fm.descent()
                    painter.drawText(int(text_x), int(text_y), label)
                    if symbol_type:
                        try:
                            ptm_type_idx = int(ptm_type.split('_')[1])
                            ptm_datasets = self.json_data.get('general_data', {}).get('data', {}).get('ptm', [])
                            if ptm_type_idx < len(ptm_datasets):
                                symbol_dict = next(
                                    (item[f"{symbol_type}_dict"] for item in ptm_datasets[ptm_type_idx].get('ptm_symbol_list', []) if f"{symbol_type}_dict" in item),
                                    None
                                )
                                if symbol_dict:
                                    symbol = symbol_dict.get('symbol', '')
                                    symbol_font = QFont(symbol_dict.get('symbol_font', 'Arial'), symbol_dict.get('symbol_size', 6))
                                    painter.setFont(symbol_font)
                                    symbol_color = symbol_dict.get('symbol_color', [0, 0, 0])
                                    if not all(isinstance(c, int) and 0 <= c <= 255 for c in symbol_color):
                                        print(f"Warning: Invalid symbol_color for PTM {ptm_key}: {symbol_color}")
                                        symbol_color = [0, 0, 0]
                                    painter.setPen(QPen(QColor(*symbol_color)))
                                    # Calculate text size for centering
                                    fm = QFontMetrics(symbol_font)
                                    text_rect = fm.boundingRect(symbol)
                                    text_width = text_rect.width()
                                    text_height = text_rect.height()
                                    # Center the symbol and apply offsets with hardcoded downward adjustment
                                    symbol_x = x - (text_width / 2) + symbol_dict.get('symbol_x_offset', 0)
                                    symbol_y = y - (text_height / 2) + symbol_dict.get('symbol_y_offset', 0) + fm.descent() + 8
                                    painter.drawText(int(symbol_x), int(symbol_y), symbol)
                                    print(f"Drew PTM symbol {symbol} for {ptm_key} at ({symbol_x}, {symbol_y})")
                            else:
                                print(f"Warning: ptm_type_idx {ptm_type_idx} out of range for ptm_datasets (length {len(ptm_datasets)})")
                        except (IndexError, ValueError) as e:
                            print(f"Error processing ptm_type {ptm_type} for PTM {ptm_key}: {e}")
                            continue
                    ptm['rect'] = QRectF(x - radius, y - radius, radius * 2, radius * 2)
                    print(f"Drew PTM {ptm_key} at ({x}, {y}) with label {label}, tooltip: {tooltip}")
        except Exception as e:
            print(f"Error in draw_ptms for protein {uniprot_id}: {e}")

    def draw_groups(self, painter):
        try:
            for group in self.groups:
                protbox_ids = group.get('protbox_ids', [])
                if not protbox_ids:
                    continue
                min_x = float('inf')
                max_x = float('-inf')
                min_y = float('inf')
                max_y = float('-inf')
                for protbox_id in protbox_ids:
                    protbox = next((pb for pb in self.protbox_data if pb.get('protbox_id') == protbox_id), None)
                    if protbox:
                        x = protbox.get('x', 0)
                        y = protbox.get('y', 0)
                        width = protbox.get('width', 0)
                        height = protbox.get('height', 0)
                        if all(isinstance(v, (int, float)) for v in [x, y, width, height]):
                            min_x = min(min_x, x)
                            max_x = max(max_x, x + width)
                            min_y = min(min_y, y)
                            max_y = max(max_y, y + height)
                        else:
                            print(f"Warning: Invalid coordinates for protbox {protbox_id} in group {group.get('group_id')}: x={x}, y={y}, width={width}, height={height}")
                if min_x != float('inf'):
                    painter.setPen(QPen(QColor(0, 0, 0), 1, Qt.DashLine))
                    painter.setBrush(Qt.NoBrush)
                    padding = 5
                    painter.drawRect(QRect(int(min_x - padding), int(min_y - padding),
                                           int(max_x - min_x + 2 * padding), int(max_y - min_y + 2 * padding)))
                    print(f"Drew group {group.get('group_id')} with protbox_ids {protbox_ids}")
        except Exception as e:
            print(f"Error in draw_groups: {e}")

    def draw_arrows(self, painter):
        try:
            painter.setPen(QPen(QColor(0, 0, 0), 1))
            for arrow in self.arrows:
                start_x = arrow.get('start_x', 0)
                start_y = arrow.get('start_y', 0)
                point_x = arrow.get('point_x', 0)
                point_y = arrow.get('point_y', 0)
                if not all(isinstance(v, (int, float)) for v in [start_x, start_y, point_x, point_y]):
                    print(f"Warning: Invalid coordinates for arrow from ({start_x}, {start_y}) to ({point_x}, {point_y})")
                    continue
                painter.drawLine(int(start_x), int(start_y), int(point_x), int(point_y))
                dx = point_x - start_x
                dy = point_y - start_y
                angle = np.arctan2(dy, dx)
                arrow_size = 10
                arrow_p1 = QPoint(int(point_x - arrow_size * np.cos(angle + np.pi / 6)),
                                  int(point_y - arrow_size * np.sin(angle + np.pi / 6)))
                arrow_p2 = QPoint(int(point_x - arrow_size * np.cos(angle - np.pi / 6)),
                                  int(point_y - arrow_size * np.sin(angle - np.pi / 6)))
                painter.drawLine(int(point_x), int(point_y), arrow_p1.x(), arrow_p1.y())
                painter.drawLine(int(point_x), int(point_y), arrow_p2.x(), arrow_p2.y())
                print(f"Drew arrow from ({start_x}, {start_y}) to ({point_x}, {point_y})")
        except Exception as e:
            print(f"Error in draw_arrows: {e}")

    def mouseMoveEvent(self, event):
        try:
            pos = event.pos()
            tooltip = None
            for protbox in self.protbox_data:
                if 'rect' in protbox and protbox['rect'].contains(pos):
                    tooltip = protbox.get('tooltip', '')
                    break
            if not tooltip:
                for uniprot_id, protein in self.protein_data.items():
                    for ptm_key, ptm in protein.get('PTMs', {}).items():
                        if 'rect' in ptm and ptm['rect'].contains(pos):
                            tooltip = ptm.get('tooltip', '')
                            print(f"PTM {ptm_key} tooltip displayed: {tooltip}")
                            break
                    if tooltip:
                        break
            if tooltip:
                QToolTip.showText(event.globalPos(), tooltip, self)
            else:
                QToolTip.hideText()
        except Exception as e:
            print(f"Error in mouseMoveEvent: {e}")

class PathwayWindow(QMainWindow):
    def __init__(self, json_data):
        super().__init__()
        self.setWindowTitle("Pathway Visualization")
        try:
            self.widget = PathwayWidget(json_data, self)
            self.setCentralWidget(self.widget)
            print("PathwayWindow initialized")
        except Exception as e:
            print(f"Error initializing PathwayWindow: {e}")

def visualize_pathway():
    try:
        if not os.path.exists(json_file_path):
            print(f"Error: JSON file not found at {json_file_path}")
            return None
        print(f"Loading JSON from {json_file_path}")
        with open(json_file_path, 'r') as f:
            json_data = json.load(f)
        print("JSON data loaded successfully")
        app = QApplication(sys.argv)
        window = PathwayWindow(json_data)
        window.show()
        print("Window shown")
        sys.exit(app.exec_())
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        return None
    except Exception as e:
        print(f"Error visualizing pathway: {e}")
        return None

if __name__ == "__main__":
    visualize_pathway()