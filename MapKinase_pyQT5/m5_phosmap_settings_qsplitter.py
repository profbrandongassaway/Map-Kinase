import sys
import os
import re
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFormLayout, QLabel, QLineEdit,
    QComboBox, QSpinBox, QPushButton, QColorDialog, QSplitter,
    QHBoxLayout, QVBoxLayout, QFrame, QScrollArea, QCheckBox, QSizePolicy,
    QListWidget, QListWidgetItem, QFileDialog, QMessageBox, QDoubleSpinBox
)
from PyQt5.QtGui import QColor, QFocusEvent, QIcon, QPixmap, QPainter, QBrush, QPolygon, QLinearGradient, QFont
from PyQt5.QtCore import Qt, QRect, QEvent, QPoint, QTimer
import pandas as pd
import m4_test as pathway_visualizer
from e1_popup_errorwindow import show_error_popup
from e5_select_organism_popup import select_organism
import logging


def resource_path(relative_path):
    if hasattr(sys, '_MEIPASS'):
        base_path = sys._MEIPASS
    else:
        base_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_path, relative_path)


def get_user_config_dir():
    return os.path.join(os.path.expanduser("~"), ".mapkinase")


class NoScrollSpinBox(QSpinBox):
    def wheelEvent(self, event):
        event.ignore()


class NoScrollComboBox(QComboBox):
    def wheelEvent(self, event):
        event.ignore()


class CustomLineEdit(QLineEdit):
    def __init__(self, parent, results_list, variables, favorite_button):
        super().__init__(parent.centralWidget())
        self.parent_window = parent
        self.results_list = results_list
        self.variables = variables
        self.favorite_button = favorite_button
        self.favorites = parent.favorites
        self.setPlaceholderText("Click or type to search pathways...")

    def focusInEvent(self, event: QFocusEvent):
        super().focusInEvent(event)
        self.show_results_list()
        self.update_results(self.text())

    def focusOutEvent(self, event: QFocusEvent):
        super().focusOutEvent(event)
        if not self.results_list.hasFocus():
            self.results_list.hide()

    def show_results_list(self):
        # Use the central widget of VariableSetterUI as the reference parent
        parent_widget = self.parent_window.centralWidget()
        # Map the bottom-left of the line edit to the parent widget
        pos = self.mapTo(parent_widget, self.rect().bottomLeft())

        list_width = self.width()
        list_height = min(self.results_list.count() * 30 + 10, 200)
        if list_height < 50:
            list_height = 50

        # Get the parent widget's geometry
        parent_rect = parent_widget.rect()

        # Adjust position if the dropdown would extend beyond the parent widget
        if pos.y() + list_height > parent_rect.height():
            pos.setY(pos.y() - list_height - self.height())
        if pos.x() + list_width > parent_rect.width():
            pos.setX(parent_rect.width() - list_width)

        # Set the geometry of the results list
        self.results_list.setGeometry(QRect(pos.x(), pos.y(), list_width, list_height))

        # Show and raise the results list
        self.results_list.show()
        self.results_list.raise_()
        self.results_list.setFocusPolicy(Qt.StrongFocus)

    def update_results(self, text):
        self.results_list.clear()
        if not text:
            favorited_vars = [var for var in self.variables if var['id'] in self.favorites]
            non_favorited_vars = [var for var in self.variables if var['id'] not in self.favorites]
            for var in favorited_vars:
                item = QListWidgetItem(var['display'])
                item.setData(Qt.UserRole, var['id'])
                pixmap = QPixmap(16, 16)
                pixmap.fill(Qt.transparent)
                painter = QPainter(pixmap)
                painter.setBrush(QBrush(QColor("yellow")))
                star_points = QPolygon([
                    QPoint(8, 0), QPoint(10, 5), QPoint(16, 5), QPoint(11, 8), QPoint(13, 13),
                    QPoint(8, 10), QPoint(3, 13), QPoint(5, 8), QPoint(0, 5), QPoint(6, 5)
                ])
                painter.drawPolygon(star_points)
                painter.end()
                item.setIcon(QIcon(pixmap))
                self.results_list.addItem(item)
            for var in non_favorited_vars:
                item = QListWidgetItem(var['display'])
                item.setData(Qt.UserRole, var['id'])
                self.results_list.addItem(item)
            self.show_results_list()
            return
        try:
            pattern = re.compile(text, re.IGNORECASE)
        except re.error:
            self.results_list.addItem("Invalid regex pattern")
            self.show_results_list()
            return
        matches = [var for var in self.variables if pattern.search(var['id']) or pattern.search(var['display'])]
        if matches:
            favorited_matches = [var for var in matches if var['id'] in self.favorites]
            non_favorited_matches = [var for var in matches if var['id'] not in self.favorites]
            for var in favorited_matches:
                item = QListWidgetItem(var['display'])
                item.setData(Qt.UserRole, var['id'])
                pixmap = QPixmap(16, 16)
                pixmap.fill(Qt.transparent)
                painter = QPainter(pixmap)
                painter.setBrush(QBrush(QColor("yellow")))
                star_points = QPolygon([
                    QPoint(8, 0), QPoint(10, 5), QPoint(16, 5), QPoint(11, 8), QPoint(13, 13),
                    QPoint(8, 10), QPoint(3, 13), QPoint(5, 8), QPoint(0, 5), QPoint(6, 5)
                ])
                painter.drawPolygon(star_points)
                painter.end()
                item.setIcon(QIcon(pixmap))
                self.results_list.addItem(item)
            for var in non_favorited_matches:
                item = QListWidgetItem(var['display'])
                item.setData(Qt.UserRole, var['id'])
                self.results_list.addItem(item)
        else:
            self.results_list.addItem("No matches found")
        self.show_results_list()


class GradientWidget(QWidget):
    def __init__(self, positive_color, negative_color, max_positive, max_negative, parent=None):
        super().__init__(parent)
        self.positive_color = positive_color
        self.negative_color = negative_color
        self.max_positive = max_positive
        self.max_negative = max_negative
        self.setMinimumHeight(35)
        self.setMaximumHeight(35)

    def update_colors(self, positive_color, negative_color):
        self.positive_color = positive_color
        self.negative_color = negative_color
        self.update()

    def update_values(self, max_positive, max_negative):
        self.max_positive = max_positive
        self.max_negative = max_negative
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        gradient_rect = QRect(0, 18, self.width(), 17)
        gradient = QLinearGradient(gradient_rect.left(), 0, gradient_rect.right(), 0)
        gradient.setColorAt(0.0, self.positive_color)
        gradient.setColorAt(0.5, Qt.white)
        gradient.setColorAt(1.0, self.negative_color)
        painter.fillRect(gradient_rect, gradient)

        painter.setPen(Qt.black)
        painter.drawRect(gradient_rect.adjusted(0, 0, -1, -1))

        font = QFont("Arial", 8)
        painter.setFont(font)
        painter.setPen(Qt.black)

        pos_text = f"{self.max_positive}"
        pos_rect = painter.boundingRect(QRect(0, 0, self.width(), 18), Qt.AlignLeft | Qt.AlignTop, pos_text)
        pos_rect.moveLeft(5)
        pos_rect.moveTop(2)
        painter.drawText(pos_rect, Qt.AlignLeft | Qt.AlignTop, pos_text)

        zero_text = "0"
        zero_rect = painter.boundingRect(QRect(0, 0, self.width(), 18), Qt.AlignHCenter | Qt.AlignTop, zero_text)
        zero_rect.moveTop(2)
        painter.drawText(zero_rect, Qt.AlignHCenter | Qt.AlignTop, zero_text)

        neg_text = f"{self.max_negative}"
        neg_rect = painter.boundingRect(QRect(0, 0, self.width(), 18), Qt.AlignRight | Qt.AlignTop, neg_text)
        neg_rect.moveRight(self.width() - 5)
        neg_rect.moveTop(2)
        painter.drawText(neg_rect, Qt.AlignRight | Qt.AlignTop, neg_text)


class VariableSetterUI(QMainWindow):
    def __init__(self, use_splitter=True, parent=None, combined_ui_parent=None):
        super().__init__(parent)
        self.setWindowTitle("KEGG Pathway Visualizer Settings")
        self.setGeometry(100, 100, 1200, 800)
        self.use_splitter = use_splitter
        self.combined_ui_parent = combined_ui_parent
        self.custom_colors_file = os.path.join(get_user_config_dir(), "custom_colors.txt")
        self.script_dir = resource_path("Scripts")
        self.favorites = set()
        self.variables = []
        self.selected_main_column = None
        self.input_data = None
        self.instance_id = id(self)
        print(f"VariableSetterUI initialized, instance_id: {self.instance_id}")
        logging.debug(f"VariableSetterUI initialized, instance_id: {self.instance_id}")

        self.organism_mapping = {
            "hsa": "Homo sapiens (human)",
            "mmu": "Mus musculus (house mouse)",
            "rno": "Rattus norvegicus (rat)",
            "sce": "Saccharomyces cerevisiae (yeast)",
            "sma": "Schistosoma mansoni",
            "dme": "Drosophila melanogaster"
        }
        self.selected_organism_code = "hsa"
        pathways_path = os.path.join(self.script_dir, "kegg_pathways.txt")
        try:
            if os.path.exists(pathways_path):
                with open(pathways_path, "r", encoding="utf-8") as f:
                    next(f)  # Skip header
                    for line in f:
                        if line.strip():
                            parts = line.strip().split("\t")
                            if len(parts) >= 2:
                                display = f"{parts[0]} - {parts[1]}"
                                self.variables.append({
                                    "name": parts[1],
                                    "id": parts[0],
                                    "display": display
                                })
            else:
                show_error_popup(f"KEGG pathways file not found: {pathways_path}", "File Error", parent=self)
                self.variables = [{"name": "No pathways available", "id": "", "display": "No pathways available"}]
        except Exception as e:
            show_error_popup(f"Error reading kegg_pathways.txt: {str(e)}", "File Error", parent=self)
            self.variables = [
                {"name": f"Error reading file: {str(e)}", "id": "", "display": f"Error reading file: {str(e)}"}]

        if self.use_splitter:
            central_widget = QWidget()
            self.setCentralWidget(central_widget)
            layout = QHBoxLayout(central_widget)
            layout.setContentsMargins(0, 0, 0, 0)
            self.splitter = QSplitter(Qt.Horizontal)
            layout.addWidget(self.splitter)
            self.settings_container = QWidget()
            settings_layout = QVBoxLayout(self.settings_container)
            settings_layout.setContentsMargins(0, 0, 0, 0)
            self.kegg_accession_widget = QWidget()
            self.setup_kegg_accession_ui(self.kegg_accession_widget)
            settings_layout.addWidget(self.kegg_accession_widget)
            scroll_area = QScrollArea()
            scroll_area.setWidgetResizable(True)
            scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
            scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
            self.settings_widget = QWidget()
            self.setup_settings_ui(self.settings_widget)
            scroll_area.setWidget(self.settings_widget)
            settings_layout.addWidget(scroll_area)
            self.splitter.addWidget(self.settings_container)
            scroll_area.setStyleSheet("""
                QScrollBar:vertical {
                    border: none;
                    background: #f0f0f0;
                    width: 12px;
                    margin: 0px 0 0px 0;
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
            self.viewer_container = QWidget()
            viewer_layout = QVBoxLayout(self.viewer_container)
            viewer_layout.setContentsMargins(0, 0, 0, 0)
            self.viewer = None
            self.splitter.addWidget(self.viewer_container)
            screen_width = QApplication.primaryScreen().availableGeometry().width()
            min_width = self.settings_container.minimumSizeHint().width()
            self.splitter.setSizes([min_width, screen_width - min_width])
        else:
            central_widget = QWidget()
            self.setCentralWidget(central_widget)
            layout = QVBoxLayout(central_widget)
            layout.setContentsMargins(0, 0, 0, 0)
            self.kegg_accession_widget = QWidget()
            self.setup_kegg_accession_ui(self.kegg_accession_widget)
            layout.addWidget(self.kegg_accession_widget)
            scroll_area = QScrollArea()
            scroll_area.setWidgetResizable(True)
            scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
            scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
            content_widget = QWidget()
            self.setup_settings_ui(content_widget)
            scroll_area.setWidget(content_widget)
            layout.addWidget(scroll_area)
            scroll_area.setStyleSheet("""
                QScrollBar:vertical {
                    border: none;
                    background: #f0f0f0;
                    width: 12px;
                    margin: 0px 0 0px 0;
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
        QTimer.singleShot(0, lambda: self.run_button.setFocus())

        # Initialize with parent data if available
        initial_data = None
        if self.combined_ui_parent and hasattr(self.combined_ui_parent, 'selected_data'):
            initial_data = self.combined_ui_parent.selected_data
        elif self.parent() and hasattr(self.parent(), 'selected_data'):
            initial_data = self.parent().selected_data
        self.update_data(initial_data)

    def update_data(self, data):
        print(f"VariableSetterUI.update_data called, instance_id: {self.instance_id}, data: {data}")
        logging.debug(f"VariableSetterUI.update_data called, instance_id: {self.instance_id}, data: {data}")

        self.input_data = data

        if not hasattr(self, 'main_column_selector') or self.main_column_selector is None:
            logging.error("main_column_selector is not initialized")
            print("Error: main_column_selector is not initialized, instance_id:", self.instance_id)
            return

        self.main_column_selector.blockSignals(True)
        logging.debug("Blocked signals for main_column_selector")

        try:
            self.main_column_selector.clear()
            logging.debug("Cleared main_column_selector")

            if not data or not isinstance(data, dict) or not data.get('protein', {}).get('main_columns', []):
                self.main_column_selector.addItem("No columns available")
                self.selected_main_column = None
                print(f"No valid protein main columns, set to 'No columns available', instance_id: {self.instance_id}")
                logging.debug(
                    f"No valid protein main columns, set to 'No columns available', instance_id: {self.instance_id}")
            else:
                main_columns = [str(col) for col in data['protein']['main_columns']]
                self.main_column_selector.addItems(main_columns)
                if not self.selected_main_column or self.selected_main_column not in main_columns:
                    self.selected_main_column = main_columns[0]
                self.main_column_selector.setCurrentText(self.selected_main_column)
                print(f"Updated main_column_selector with columns: {main_columns}, instance_id: {self.instance_id}")
                logging.debug(
                    f"Updated main_column_selector with columns: {main_columns}, instance_id: {self.instance_id}")

            # Process PTM data, skip ptm_symbol_list validation
            if data and 'ptm' in data:
                for ptm in data['ptm']:
                    # Log ptm_symbol_list for debugging
                    if 'ptm_symbol_list' in ptm and ptm['ptm_symbol_list']:
                        logging.debug(
                            f"Preserving ptm_symbol_list for PTM type {ptm.get('type')}: {ptm['ptm_symbol_list']}")
                    # Remove deprecated fields
                    ptm.pop('phos_annotation_column1', None)
                    ptm.pop('phos_annotation_column2', None)

        except Exception as e:
            logging.error(f"Error updating main_column_selector: {str(e)}")
            print(f"Error updating main_column_selector: {str(e)}, instance_id: {self.instance_id}")
            self.main_column_selector.clear()
            self.main_column_selector.addItem("No columns available")
            self.selected_main_column = None

        finally:
            self.main_column_selector.blockSignals(False)
            logging.debug("Unblocked signals for main_column_selector")

    def setup_kegg_accession_ui(self, widget):
        layout = QFormLayout(widget)
        layout.addRow(QLabel("<b>KEGG Accession</b>"))
        pathway_widget = QWidget()
        pathway_layout = QHBoxLayout(pathway_widget)
        pathway_layout.setSpacing(10)
        self.results_list = QListWidget(self.centralWidget())
        self.results_list.setStyleSheet("""
            QListWidget {
                border: 1px solid gray;
                background: white;
                border-radius: 4px;
                box-shadow: 2px 2px 5px rgba(0,0,0,0.2);
            }
            QListWidget::item {
                padding: 5px;
            }
            QListWidget::item:hover {
                background: #e0e0e0;
            }
        """)
        self.results_list.hide()
        self.results_list.focusOutEvent = lambda \
            event: self.results_list.hide() if not self.pathway_id.hasFocus() else None
        self.results_list.itemClicked.connect(self.on_item_clicked)
        self.favorite_button = QPushButton(self.centralWidget())
        self.favorite_button.setFixedSize(24, 24)
        self.update_favorite_button("")
        self.favorite_button.clicked.connect(self.toggle_favorite)
        self.pathway_id = CustomLineEdit(self, self.results_list, self.variables, self.favorite_button)
        self.pathway_id.setMinimumWidth(350)
        self.pathway_id.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.pathway_id.textChanged.connect(self.on_search_text_changed)
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_visualization)
        self.run_button.setMaximumWidth(70)
        pathway_layout.addWidget(self.favorite_button)
        pathway_layout.addWidget(self.pathway_id, stretch=1)
        pathway_layout.addWidget(self.run_button)
        pathway_layout.addStretch()
        layout.addRow(QLabel("Pathway ID:"), pathway_widget)
        button_layout = QHBoxLayout()
        button_layout.setSpacing(10)
        self.save_button = QPushButton("Save Settings")
        self.save_button.setMaximumWidth(120)
        self.save_button.clicked.connect(self.save_settings)
        button_layout.addWidget(self.save_button)
        self.load_button = QPushButton("Load Settings")
        self.load_button.setMaximumWidth(120)
        self.load_button.clicked.connect(self.load_settings)
        button_layout.addWidget(self.load_button)
        button_layout.addStretch()
        layout.addRow(QLabel("Settings:"), button_layout)

    def on_item_clicked(self, item):
        pathway_id = item.data(Qt.UserRole)
        if pathway_id:
            self.pathway_id.setText(pathway_id)
            self.results_list.hide()
            self.pathway_id.setFocus()
            self.update_favorite_button(pathway_id)

    def on_search_text_changed(self, text):
        self.pathway_id.update_results(text)
        self.update_favorite_button(text)

    def update_favorite_button(self, pathway_id):
        is_favorited = pathway_id in self.favorites
        pixmap = QPixmap(16, 16)
        pixmap.fill(Qt.transparent)
        painter = QPainter(pixmap)
        painter.setBrush(QBrush(QColor("yellow" if is_favorited else "gray")))
        star_points = QPolygon([
            QPoint(8, 0), QPoint(10, 5), QPoint(16, 5), QPoint(11, 8), QPoint(13, 13),
            QPoint(8, 10), QPoint(3, 13), QPoint(5, 8), QPoint(0, 5), QPoint(6, 5)
        ])
        painter.drawPolygon(star_points)
        painter.end()
        self.favorite_button.setIcon(QIcon(pixmap))
        self.favorite_button.setEnabled(bool(pathway_id))

    def toggle_favorite(self):
        pathway_id = self.pathway_id.text()
        if not pathway_id:
            return
        if pathway_id in self.favorites:
            self.favorites.remove(pathway_id)
        else:
            self.favorites.add(pathway_id)
        self.update_favorite_button(pathway_id)
        self.pathway_id.update_results(self.pathway_id.text())

    def select_organism_and_update(self):
        selected_code = select_organism()
        if selected_code:
            self.selected_organism_code = selected_code
            organism_name = self.organism_mapping.get(selected_code, "Custom: " + selected_code)
            self.organism_label.setText(organism_name)

    def setup_settings_ui(self, widget):
        self.layout = QFormLayout(widget)
        self.layout.setSpacing(3)
        self.layout.addRow(QLabel("<b>Input</b>"))
        organism_widget = QWidget()
        organism_layout = QHBoxLayout(organism_widget)
        organism_layout.setContentsMargins(0, 0, 0, 0)
        organism_layout.setSpacing(10)
        self.organism_button = QPushButton("Select Organism")
        self.organism_button.setMaximumWidth(120)
        self.organism_button.clicked.connect(self.select_organism_and_update)
        self.organism_label = QLabel(self.organism_mapping.get(self.selected_organism_code, ""))
        organism_layout.addWidget(self.organism_button, alignment=Qt.AlignLeft)
        organism_layout.addWidget(self.organism_label, alignment=Qt.AlignLeft)
        organism_layout.addStretch()
        self.layout.addRow(QLabel("Organism:"), organism_widget)
        self.layout.addRow(QLabel(""))

        self.layout.addRow(QLabel("<b>Settings</b>"))
        self.main_column_selector = NoScrollComboBox()
        self.main_column_selector.addItem("No columns available")
        self.main_column_selector.setMaximumWidth(300)
        self.main_column_selector.currentTextChanged.connect(self.update_selected_main_column)
        self.layout.addRow(QLabel("Main Column to Visualize:"), self.main_column_selector)
        self.protein_selection_option = NoScrollComboBox()
        self.protein_selection_option.addItems([
            "First in List (1)",
            "Protein Fold Change (2)",
            "Regulatory Phosphosite Fold Change (3)"
        ])
        self.protein_selection_option.setCurrentIndex(2)
        self.layout.addRow(QLabel("Protein Selection Option:"), self.protein_selection_option)
        self.protein_selection_option.setMaximumWidth(300)
        self.ptm_selection_option = NoScrollComboBox()
        self.ptm_selection_option.addItems([
            "Largest Fold Change (1)",
            "Modulating Sites Prioritized (2)",
            "Only Modulating Sites (3)"
        ])
        self.ptm_selection_option.setCurrentIndex(1)
        self.layout.addRow(QLabel("PTM Selection Option:"), self.ptm_selection_option)
        self.ptm_selection_option.setMaximumWidth(300)
        self.ptm_max_display = NoScrollSpinBox()
        self.ptm_max_display.setRange(1, 10)
        self.ptm_max_display.setValue(4)
        self.layout.addRow(QLabel("Max PTM Sites to Display:"), self.ptm_max_display)
        self.ptm_max_display.setMaximumWidth(50)

        colors_widget = QWidget()
        colors_layout = QVBoxLayout(colors_widget)
        colors_layout.setSpacing(5)

        self.max_positive = NoScrollSpinBox()
        self.max_positive.setRange(0, 10)
        self.max_positive.setValue(2)
        self.max_positive.setMaximumWidth(60)
        self.max_positive.valueChanged.connect(self.update_gradient_values)
        self.max_negative = NoScrollSpinBox()
        self.max_negative.setRange(-10, 0)
        self.max_negative.setValue(-2)
        self.max_negative.setMaximumWidth(60)
        self.max_negative.valueChanged.connect(self.update_gradient_values)

        gradient_row = QWidget()
        gradient_layout = QHBoxLayout(gradient_row)
        gradient_layout.setSpacing(10)
        self.positive_color = QColor(205, 41, 41)
        self.negative_color = QColor(0, 143, 144)
        self.gradient_widget = GradientWidget(
            self.positive_color,
            self.negative_color,
            self.max_positive.value(),
            self.max_negative.value()
        )
        self.positive_color_button = QPushButton("")
        self.positive_color_button.setFixedSize(30, 30)
        self.positive_color_button.setStyleSheet(
            f"QPushButton {{ background-color: {self.positive_color.name()}; "
            "border: 2px solid #555555; "
            "border-radius: 5px; }}"
            "QPushButton:hover {{ border: 2px solid #333333; }}"
            "QPushButton:pressed {{ border: 2px solid #111111; }}"
        )
        self.positive_color_button.clicked.connect(self.pick_positive_color)
        self.negative_color_button = QPushButton("")
        self.negative_color_button.setFixedSize(30, 30)
        self.negative_color_button.setStyleSheet(
            f"QPushButton {{ background-color: {self.negative_color.name()}; "
            "border: 2px solid #555555; "
            "border-radius: 5px; }}"
            "QPushButton:hover {{ border: 2px solid #333333; }}"
            "QPushButton:pressed {{ border: 2px solid #111111; }}"
        )
        self.negative_color_button.clicked.connect(self.pick_negative_color)
        gradient_layout.addWidget(self.positive_color_button)
        gradient_layout.addWidget(self.gradient_widget, stretch=1)
        gradient_layout.addWidget(self.negative_color_button)
        colors_layout.addWidget(gradient_row)

        values_row = QWidget()
        values_layout = QHBoxLayout(values_row)
        values_layout.setSpacing(10)
        values_layout.addStretch()
        values_layout.addWidget(QLabel("Max Positive:"))
        values_layout.addWidget(self.max_positive)
        values_layout.addWidget(QLabel("Max Negative:"))
        values_layout.addWidget(self.max_negative)
        values_layout.addStretch()
        colors_layout.addWidget(values_row)

        self.layout.addRow(colors_widget)
        self.show_multi_protein_indicator = QCheckBox("Show Multi-Protein Indicator (*)")
        self.show_multi_protein_indicator.setChecked(True)
        self.layout.addRow(self.show_multi_protein_indicator)

    def update_gradient_values(self):
        self.gradient_widget.update_values(
            self.max_positive.value(),
            self.max_negative.value()
        )

    def update_selected_main_column(self, text):
        if text and text != "No columns available":
            self.selected_main_column = text
        else:
            self.selected_main_column = None

    def load_custom_colors(self, color_dialog):
        if not os.path.exists(self.custom_colors_file):
            bundled_colors_file = resource_path("Scripts/custom_colors.txt")
            if os.path.exists(bundled_colors_file):
                os.makedirs(os.path.dirname(self.custom_colors_file), exist_ok=True)
                import shutil
                shutil.copy(bundled_colors_file, self.custom_colors_file)
            else:
                return
        try:
            with open(self.custom_colors_file, 'r') as f:
                lines = f.readlines()
        except (IOError, PermissionError) as e:
            show_error_popup(f"Failed to read custom colors file: {str(e)}", "File Error", parent=self)
            return
        for i, line in enumerate(lines[:16]):
            try:
                r, g, b = map(int, line.strip().split(','))
                if not (0 <= r <= 255 and 0 <= g <= 255 and 0 <= b <= 255):
                    raise ValueError("Color values out of range")
                color = QColor(r, g, b)
                color_dialog.setCustomColor(i, color)
            except (ValueError, IndexError) as e:
                show_error_popup(f"Invalid color format in line {i + 1}: {str(e)}", "Color Parse Error", parent=self)
                continue

    def save_custom_colors(self, color_dialog):
        try:
            os.makedirs(os.path.dirname(self.custom_colors_file), exist_ok=True)
            with open(self.custom_colors_file, 'w', encoding='utf-8') as f:
                for i in range(16):
                    color = color_dialog.customColor(i)
                    f.write(f"{color.red()},{color.green()},{color.blue()}\n")
        except (IOError, PermissionError) as e:
            show_error_popup(f"Failed to save custom colors: {str(e)}", "File Error", parent=self)

    def pick_positive_color(self):
        try:
            color_dialog = QColorDialog(self.positive_color, self)
            self.load_custom_colors(color_dialog)
            if color_dialog.exec_():
                self.positive_color = color_dialog.selectedColor()
                self.positive_color_button.setStyleSheet(
                    f"QPushButton {{ background-color: {self.positive_color.name()}; "
                    "border: 2px solid #555555; "
                    "border-radius: 5px; }}"
                    "QPushButton:hover {{ border: 2px solid #333333; }}"
                    "QPushButton:pressed {{ border: 2px solid #111111; }}"
                )
                self.gradient_widget.update_colors(self.positive_color, self.negative_color)
                self.save_custom_colors(color_dialog)
        except RuntimeError as e:
            show_error_popup(f"Error in color picker: {str(e)}", "GUI Error", parent=self)

    def pick_negative_color(self):
        try:
            color_dialog = QColorDialog(self.negative_color, self)
            self.load_custom_colors(color_dialog)
            if color_dialog.exec_():
                self.negative_color = color_dialog.selectedColor()
                self.negative_color_button.setStyleSheet(
                    f"QPushButton {{ background-color: {self.negative_color.name()}; "
                    "border: 2px solid #555555; "
                    "border-radius: 5px; }}"
                    "QPushButton:hover {{ border: 2px solid #333333; }}"
                    "QPushButton:pressed {{ border: 2px solid #111111; }}"
                )
                self.gradient_widget.update_colors(self.positive_color, self.negative_color)
                self.save_custom_colors(color_dialog)
        except RuntimeError as e:
            show_error_popup(f"Error in color picker: {str(e)}", "GUI Error", parent=self)

    def get_current_settings(self):
        vis_settings = {
            'prot_label_font': 'Arial',
            'prot_label_size': 6,
            'ptm_label_font': 'Arial',
            'ptm_label_size': 4,
            'ptm_label_color': (0, 0, 0),
            'ptm_circle_radius': 5,
            'ptm_circle_spacing': 2,
            'protein_tooltip_columns': []
        }
        parent = self.combined_ui_parent if self.combined_ui_parent else self.parent()
        input_data = parent.selected_data if parent and hasattr(parent, 'selected_data') else None
        if input_data and input_data.get('visualization_settings'):
            vis_settings.update(input_data['visualization_settings'])
        settings = {
            'selected_organism_code': self.selected_organism_code,
            'pathway_id': self.selected_pathway_id if hasattr(self, 'selected_pathway_id') else self.construct_pathway_id(self.pathway_id.text()),
            'favorites': list(self.favorites),
            'selected_main_column': self.selected_main_column if self.selected_main_column else '',
            'protein_selection_option': self.protein_selection_option.currentIndex() + 1,
            'ptm_selection_option': self.ptm_selection_option.currentIndex() + 1,
            'ptm_max_display': self.ptm_max_display.value(),
            'negative_color': (self.negative_color.red(), self.negative_color.green(), self.negative_color.blue()),
            'positive_color': (self.positive_color.red(), self.positive_color.green(), self.positive_color.blue()),
            'max_negative': self.max_negative.value(),
            'max_positive': self.max_positive.value(),
            'prot_label_font': vis_settings['prot_label_font'],
            'prot_label_size': vis_settings['prot_label_size'],
            'ptm_label_font': vis_settings['ptm_label_font'],
            'ptm_label_color': vis_settings['ptm_label_color'],
            'ptm_circle_radius': vis_settings['ptm_circle_radius'],
            'show_multi_protein_indicator': self.show_multi_protein_indicator.isChecked(),
            'protein_tooltip_columns': vis_settings['protein_tooltip_columns']
        }
        if input_data:
            protein_main_columns = input_data.get('protein', {}).get('main_columns', [])
            settings.update({
                'protein_file_path': input_data.get('protein', {}).get('file_path', ''),
                'protein_main_columns': [str(col) for col in protein_main_columns],
                'protein_kegg_column': input_data.get('protein', {}).get('kegg_column', ''),
                'protein_uniprot_column': input_data.get('protein', {}).get('uniprot_column', ''),
                'protein_gene_column': input_data.get('protein', {}).get('gene_column', '')
            })
            ptm_datasets = []
            for ptm in input_data.get('ptm', []):
                ptm_main_columns = ptm.get('main_columns', [])
                ptm_main_columns_str = [f"{p}:{t}" for p, t in ptm_main_columns] if ptm_main_columns else []
                ptm_symbol_list = ptm.get('ptm_symbol_list', [])
                ptm_entry = [
                    ptm.get('file_path', ''),
                    ptm.get('type', ''),
                    ','.join(ptm_main_columns_str),
                    ptm.get('site_column', ''),
                    ptm.get('modulation_column', ''),
                    str(ptm_symbol_list) if ptm_symbol_list else ''
                ]
                ptm_datasets.append('|'.join([str(field) for field in ptm_entry]))
            settings['ptm_datasets'] = ptm_datasets
        else:
            settings.update({
                'protein_file_path': '',
                'protein_main_columns': [],
                'ptm_datasets': [],
                'protein_tooltip_columns': []
            })
        return settings

    def save_settings(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Settings", "", "Text Files (*.txt);;All Files (*)"
        )
        if file_path:
            try:
                settings = self.get_current_settings()
                with open(file_path, 'w', encoding='utf-8') as f:
                    for key, value in settings.items():
                        try:
                            if key in ['negative_color', 'positive_color', 'ptm_label_color']:
                                value = f"{value[0]},{value[1]},{value[2]}"
                            elif key == 'favorites':
                                value = ";".join(str(item) for item in value)
                            elif key == 'protein_main_columns':
                                value = ",".join(str(item) for item in value)
                            elif key == 'ptm_datasets':
                                value = ";".join(str(item) for item in value)
                            else:
                                value = str(value)
                            f.write(f"{key}={value}\n")
                        except TypeError as e:
                            show_error_popup(
                                f"Serialization error for key '{key}': {str(e)}. Value: {value}",
                                "Settings Error",
                                parent=self
                            )
                            raise
            except (IOError, PermissionError) as e:
                show_error_popup(f"Failed to save settings: {str(e)}", "File Error", parent=self)
            except UnicodeEncodeError as e:
                show_error_popup(f"Cannot save settings due to invalid characters: {str(e)}", "Encoding Error",
                                 parent=self)
            except Exception as e:
                show_error_popup(f"Unexpected error saving settings: {str(e)}", "Settings Error", parent=self)

    def load_settings(self):
        print("DEBUG: load_settings method called")
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Settings", "", "Text Files (*.txt);;All Files (*)"
        )
        if not file_path:
            return
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                settings = {}
                for line in f:
                    if '=' in line:
                        try:
                            key, value = line.strip().split('=', 1)
                            settings[key] = value
                        except ValueError:
                            continue
        except (IOError, PermissionError) as e:
            show_error_popup(f"Failed to load settings file: {str(e)}", "File Error", parent=self)
            return
        try:
            input_data = {'protein': None, 'ptm': [], 'visualization_settings': {}}

            if 'protein_file_path' in settings and settings['protein_file_path']:
                protein_main_columns = settings.get('protein_main_columns', '').split(',') if settings.get(
                    'protein_main_columns') else []
                protein_main_columns = [col for col in protein_main_columns if col]
                input_data['protein'] = {
                    'file_path': settings['protein_file_path'],
                    'main_columns': protein_main_columns,
                    'kegg_column': settings.get('protein_kegg_column', ''),
                    'uniprot_column': settings.get('protein_uniprot_column', ''),
                    'gene_column': settings.get('protein_gene_column', '')
                }
                if not os.path.exists(settings['protein_file_path']):
                    show_error_popup(f"Protein file not found: {settings['protein_file_path']}", "File Warning",
                                     parent=self)

            if 'ptm_datasets' in settings and settings['ptm_datasets']:
                ptm_datasets = settings['ptm_datasets'].split(';')
                for ptm_entry in ptm_datasets:
                    fields = ptm_entry.split('|')
                    if len(fields) >= 5:
                        ptm_main_columns_str = fields[2].split(',') if fields[2] else []
                        ptm_main_columns = []
                        for col_str in ptm_main_columns_str:
                            if ':' in col_str:
                                prot_col, ptm_col = col_str.split(':', 1)
                                ptm_main_columns.append([prot_col, ptm_col])
                            else:
                                ptm_main_columns.append([col_str, col_str])
                        modulation_column = fields[4] if fields[4] else 'Regulatory site'
                        print("MODUCOLUMNHERE", modulation_column)
                        ptm_dataset = {
                            'file_path': fields[0],
                            'type': fields[1],
                            'main_columns': ptm_main_columns,
                            'site_column': fields[3],
                            'modulation_column': modulation_column,
                            'ptm_symbol_list': eval(fields[5]) if len(fields) > 5 and fields[5] else []
                        }
                        logging.debug(f"Parsed PTM dataset: {ptm_dataset}")
                        if not os.path.exists(fields[0]):
                            show_error_popup(f"PTM file not found: {fields[0]}", "File Warning", parent=self)
                        input_data['ptm'].append(ptm_dataset)

            vis_settings = {}
            for key in ['prot_label_font', 'ptm_label_font']:
                if key in settings:
                    vis_settings[key] = settings[key]
            for key in ['prot_label_size', 'ptm_label_size', 'ptm_circle_radius']:
                if key in settings and settings[key].isdigit():
                    vis_settings[key] = int(settings[key])
            if 'ptm_label_color' in settings:
                try:
                    r, g, b = map(int, settings['ptm_label_color'].split(','))
                    vis_settings['ptm_label_color'] = (r, g, b)
                except (ValueError, IndexError):
                    vis_settings['ptm_label_color'] = (0, 0, 0)
            vis_settings['ptm_circle_spacing'] = 2
            input_data['visualization_settings'] = vis_settings

            self.main_column_selector.clear()
            if input_data and input_data.get('protein', {}).get('main_columns', []):
                main_columns = input_data['protein']['main_columns']
                self.main_column_selector.addItems([str(col) for col in main_columns])
                self.selected_main_column = settings.get('selected_main_column', str(main_columns[0]))
                if self.selected_main_column in main_columns:
                    self.main_column_selector.setCurrentText(self.selected_main_column)
                else:
                    self.selected_main_column = str(main_columns[0])
                    self.main_column_selector.setCurrentText(self.selected_main_column)
            else:
                self.main_column_selector.addItem("No columns available")
                self.selected_main_column = None

            if 'selected_organism_code' in settings:
                self.selected_organism_code = settings['selected_organism_code']
                organism_name = self.organism_mapping.get(self.selected_organism_code,
                                                          "Custom: " + self.selected_organism_code)
                self.organism_label.setText(organism_name)

            if 'pathway_id' in settings:
                self.pathway_id.setText(settings['pathway_id'])
                self.update_favorite_button(settings['pathway_id'])

            if 'favorites' in settings:
                self.favorites = set(settings['favorites'].split(';')) if settings['favorites'] else set()
                self.pathway_id.favorites = self.favorites
                self.pathway_id.update_results(self.pathway_id.text())

            if 'protein_selection_option' in settings and settings['protein_selection_option'].isdigit():
                index = int(settings['protein_selection_option']) - 1
                if 0 <= index < self.protein_selection_option.count():
                    self.protein_selection_option.setCurrentIndex(index)

            if 'ptm_selection_option' in settings and settings['ptm_selection_option'].isdigit():
                index = int(settings['ptm_selection_option']) - 1
                if 0 <= index < self.ptm_selection_option.count():
                    self.ptm_selection_option.setCurrentIndex(index)

            if 'ptm_max_display' in settings and settings['ptm_max_display'].isdigit():
                self.ptm_max_display.setValue(int(settings['ptm_max_display']))

            if 'negative_color' in settings:
                try:
                    r, g, b = map(int, settings['negative_color'].split(','))
                    if 0 <= r <= 255 and 0 <= g <= 255 and 0 <= b <= 255:
                        self.negative_color = QColor(r, g, b)
                        self.negative_color_button.setStyleSheet(
                            f"QPushButton {{ background-color: {self.negative_color.name()}; "
                            "border: 2px solid #555555; "
                            "border-radius: 5px; }}"
                            "QPushButton:hover {{ border: 2px solid #333333; }}"
                            "QPushButton:pressed {{ border: 2px solid #111111; }}"
                        )
                except (ValueError, IndexError):
                    pass

            if 'positive_color' in settings:
                try:
                    r, g, b = map(int, settings['positive_color'].split(','))
                    if 0 <= r <= 255 and 0 <= g <= 255 and 0 <= b <= 255:
                        self.positive_color = QColor(r, g, b)
                        self.positive_color_button.setStyleSheet(
                            f"QPushButton {{ background-color: {self.positive_color.name()}; "
                            "border: 2px solid #555555; "
                            "border-radius: 5px; }}"
                            "QPushButton:hover {{ border: 2px solid #333333; }}"
                            "QPushButton:pressed {{ border: 2px solid #111111; }}"
                        )
                except (ValueError, IndexError):
                    pass

            self.gradient_widget.update_colors(self.positive_color, self.negative_color)

            if 'max_negative' in settings and settings['max_negative'].lstrip('-').isdigit():
                self.max_negative.setValue(int(settings['max_negative']))
            if 'max_positive' in settings and settings['max_positive'].isdigit():
                self.max_positive.setValue(int(settings['max_positive']))

            if 'show_multi_protein_indicator' in settings:
                self.show_multi_protein_indicator.setChecked(
                    settings['show_multi_protein_indicator'].lower() == 'true')

        except Exception as e:
            show_error_popup(f"Error applying settings: {str(e)}", "Settings Error", parent=self)

    def on_pathway_changed(self, text):
        selected = next((v for v in self.variables if v["display"] == text), None)
        if selected and selected["id"]:
            self.selected_pathway_id = selected["id"]
            print(f"Selected pathway: {self.selected_pathway_id}")
            logging.debug(f"Selected pathway: {self.selected_pathway_id}")

    def on_main_column_changed(self, text):
        self.selected_main_column = text if text != "No columns available" else None
        print(f"Selected main column: {self.selected_main_column}")
        logging.debug(f"Selected main column: {self.selected_main_column}")

    def construct_pathway_id(self, input_id):
        """
        Constructs a full KEGG pathway ID by combining the selected organism code
        with the input pathway ID if necessary. Does not validate against kegg_pathways.txt
        during development.

        Args:
            input_id (str): The user-entered pathway ID (e.g., '04150', 'map04150', 'hsa04150')

        Returns:
            str: The constructed full pathway ID (e.g., 'hsa04150') or empty string if invalid
        """
        if not input_id:
            print("Input pathway ID is empty")
            logging.debug("Input pathway ID is empty")
            return ""

        # Remove leading/trailing whitespace
        input_id = input_id.strip()
        print(f"Constructing pathway ID from input: {input_id}, organism: {self.selected_organism_code}")
        logging.debug(f"Constructing pathway ID from input: {input_id}, organism: {self.selected_organism_code}")

        # List of valid organism codes (from self.organism_mapping or extended list)
        valid_organisms = list(self.organism_mapping.keys())  # e.g., ['hsa', 'mmu', 'rno', ...]

        # Check if input starts with a valid organism code followed by 5 digits
        if re.match(r'^[a-z]{3,4}\d{5}$', input_id):
            organism_code = input_id[:3] if len(input_id) == 8 else input_id[:4]
            if organism_code in valid_organisms:
                full_id = input_id
                print(f"Input is a valid full KEGG ID: {full_id}")
                logging.debug(f"Input is a valid full KEGG ID: {full_id}")
            else:
                # Treat as partial ID (e.g., 'map04150') and use selected organism
                numeric_part = re.sub(r'^[a-z]+', '', input_id, flags=re.IGNORECASE).strip()
                if not numeric_part.isdigit():
                    print(f"Invalid pathway ID: {input_id} (no numeric part)")
                    logging.debug(f"Invalid pathway ID: {input_id} (no numeric part)")
                    return ""
                full_id = f"{self.selected_organism_code}{numeric_part.zfill(5)}"
                print(f"Constructed from invalid organism prefix: {full_id}")
                logging.debug(f"Constructed from invalid organism prefix: {full_id}")
        else:
            # Handle partial IDs (e.g., '04150', 'map04150')
            numeric_part = re.sub(r'^map', '', input_id, flags=re.IGNORECASE).strip()
            if not numeric_part.isdigit():
                print(f"Invalid pathway ID: {input_id} (no numeric part)")
                logging.debug(f"Invalid pathway ID: {input_id} (no numeric part)")
                return ""
            full_id = f"{self.selected_organism_code}{numeric_part.zfill(5)}"
            print(f"Constructed from partial ID: {full_id}")
            logging.debug(f"Constructed from partial ID: {full_id}")

        # Validate format of constructed ID
        if not re.match(r'^[a-z]{3,4}\d{5}$', full_id) or full_id[:3] not in valid_organisms and full_id[
                                                                                                 :4] not in valid_organisms:
            print(f"Invalid constructed pathway ID format or organism: {full_id}")
            logging.debug(f"Invalid constructed pathway ID format or organism: {full_id}")
            return ""

        print(f"Final constructed pathway ID: {full_id}")
        logging.debug(f"Final constructed pathway ID: {full_id}")
        return full_id

    def run_visualization(self):
        self.results_list.hide()
        print(f"Running visualization, instance_id: {self.instance_id}")
        logging.debug(f"Running visualization, instance_id: {self.instance_id}")

        input_data = self.input_data
        if input_data is None and self.combined_ui_parent and hasattr(self.combined_ui_parent, 'selected_data'):
            input_data = self.combined_ui_parent.selected_data
            print(f"Fallback: Retrieved data from combined_ui_parent: {input_data}")
            logging.debug(f"Fallback: Retrieved data from combined_ui_parent: {input_data}")
            self.input_data = input_data

        print(f"Input data in run_visualization: {input_data}")
        logging.debug(f"Input data in run_visualization: {input_data}")

        parent = self.parent()
        parent_has_data = parent and hasattr(parent, 'selected_data') and parent.selected_data is not None
        print(f"Parent context: {parent}, has selected_data: {parent_has_data}")
        logging.debug(f"Parent context: {parent}, has selected_data: {parent_has_data}")

        if not input_data:
            show_error_popup("No data available. Please configure data in the Data Selection tab.", "Data Error",
                             parent=self)
            print("Error: No input data available")
            return

        protein_data = input_data.get('protein', {})
        ptm_data = input_data.get('ptm', [])
        if not protein_data.get('file_path'):
            show_error_popup("No protein file selected.", "Data Error", parent=self)
            print("Error: No protein file selected")
            return

        try:
            protein_df = pd.read_csv(protein_data['file_path'], sep='\t')
            required_columns = [protein_data.get('uniprot_column'), protein_data.get('kegg_column'),
                                protein_data.get('gene_column')] + protein_data.get('main_columns', [])
            required_columns = [c for c in required_columns if c]
            missing_columns = [c for c in required_columns if c not in protein_df.columns]
            if missing_columns:
                show_error_popup(f"Missing columns in protein file: {missing_columns}", "Data Error", parent=self)
                print(f"Error: Missing columns in protein file: {missing_columns}")
                return
            print("Protein file validated successfully")
            logging.debug(f"Protein file columns: {list(protein_df.columns)}")
        except Exception as e:
            show_error_popup(f"Error reading protein file: {str(e)}", "File Error", parent=self)
            print(f"Error reading protein file: {str(e)}")
            return

        if ptm_data:
            for ptm in ptm_data:
                try:
                    ptm_df = pd.read_csv(ptm['file_path'], sep='\t')
                    required_columns = [ptm.get('uniprot_column'), ptm.get('site_column')]
                    required_columns = [c for c in required_columns if c]
                    missing_columns = [c for c in required_columns if c not in ptm_df.columns]
                    if missing_columns:
                        show_error_popup(f"Missing columns in PTM file: {missing_columns}", "Data Error", parent=self)
                        print(f"Error: Missing columns in PTM file: {missing_columns}")
                        return
                    print("PTM file validated successfully")
                    logging.debug(f"PTM file columns: {list(ptm_df.columns)}")
                except Exception as e:
                    show_error_popup(f"Error reading PTM file: {str(e)}", "File Error", parent=self)
                    print(f"Error reading PTM file: {str(e)}")
                    return

        input_pathway_id = self.pathway_id.text()
        print(f"Raw pathway input: {input_pathway_id}")
        logging.debug(f"Raw pathway input: {input_pathway_id}")
        pathway_id = self.construct_pathway_id(input_pathway_id)
        if not pathway_id:
            show_error_popup(
                f"Invalid pathway ID: {input_pathway_id}. Please enter a valid KEGG pathway ID (e.g., 'hsa04150' or '04150').",
                "Input Error", parent=self)
            print(f"Error: Invalid pathway ID: {input_pathway_id}")
            return
        self.selected_pathway_id = pathway_id
        print(f"Final pathway ID for visualization: {self.selected_pathway_id}")
        logging.debug(f"Final pathway ID for visualization: {self.selected_pathway_id}")

        self.update_favorite_button(self.selected_pathway_id)

        settings = self.get_current_settings()
        settings['fold_change_column'] = self.selected_main_column if self.selected_main_column else ''
        settings['hsa_id_column'] = input_data.get('protein', {}).get('kegg_column', 'KEGG_hsa')
        settings['prot_uniprot_column'] = input_data.get('protein', {}).get('uniprot_column', 'Uniprot_ID')
        settings['gene_name_column'] = input_data.get('protein', {}).get('gene_column', 'Gene Symbol')
        settings['show_background_image'] = True
        settings['display_types'] = ['gene']
        settings['show_groups'] = False
        print(f"Visualization settings prepared: {settings}")
        logging.debug(f"Visualization settings prepared: {settings}")

        data = {
            'protein': {
                'file_path': protein_data.get('file_path', ''),
                'uniprot_column': protein_data.get('uniprot_column', 'Uniprot_ID'),
                'kegg_column': protein_data.get('kegg_column', 'KEGG_hsa'),
                'gene_column': protein_data.get('gene_column', 'Gene Symbol'),
                'main_columns': protein_data.get('main_columns', []),
                'tooltip_columns': protein_data.get('tooltip_columns', [])
            },
            'ptm': [
                {
                    'type': ptm.get('type', ''),
                    'file_path': ptm.get('file_path', ''),
                    'uniprot_column': ptm.get('uniprot_column', ''),
                    'site_column': ptm.get('site_column', ''),
                    'shape': ptm.get('shape', 'Circle'),
                    'main_columns': ptm.get('main_columns', []),
                    'tooltip_columns': ptm.get('tooltip_columns', []),
                    'modulation_column': ptm.get('modulation_column', ''),
                    'ptm_symbol_list': ptm.get('ptm_symbol_list', [])
                } for ptm in ptm_data
            ]
        }
        print(f"Visualization data prepared: {data}")
        logging.debug(f"Visualization data prepared: {data}")

        try:
            print(f"Calling visualize_kegg_pathway with pathway_id: {self.selected_pathway_id}")
            logging.debug(f"Calling visualize_kegg_pathway with pathway_id: {self.selected_pathway_id}")

            viewer = pathway_visualizer.visualize_kegg_pathway(
                pathway_id=self.selected_pathway_id,
                data=data,
                settings=settings
            )

            if viewer:
                print(f"Visualization completed successfully, viewer created")
                logging.debug(f"Visualization completed successfully, viewer created")
                if self.use_splitter:
                    if self.viewer:
                        self.viewer_container.layout().removeWidget(self.viewer)
                        self.viewer.deleteLater()
                    self.viewer = viewer
                    self.viewer_container.layout().addWidget(self.viewer)
                    self.viewer.show()
                else:
                    viewer.show()
                self.run_button.setFocus()
            else:
                show_error_popup("Visualization failed to generate output.", "Visualization Error", parent=self)
                print("Error: Visualization returned no output")
                return

        except Exception as e:
            show_error_popup(f"Error during visualization: {str(e)}", "Visualization Error", parent=self)
            print(f"Error during visualization: {str(e)}")
            self.run_button.setFocus()
            return


    def get_current_data(self):
        """Return the current state of VariableSetterUI for saving."""
        logging.debug(f"VariableSetterUI.get_current_data called, instance_id: {self.instance_id}")
        print(f"VariableSetterUI.get_current_data called, instance_id: {self.instance_id}")

        try:
            settings = self.get_current_settings()
            data = {
                'settings': settings,
                'input_data': self.input_data if self.input_data else {},
            }
            logging.debug(f"Collected data: {data}")
            return data
        except Exception as e:
            logging.error(f"Error in get_current_data: {str(e)}")
            print(f"Error in get_current_data: {str(e)}")
            return {'error': str(e)}

    def set_current_data(self, data):
        logging.debug(f"VariableSetterUI.set_current_data called with data: {data}")
        print(f"VariableSetterUI.set_current_data called with data: {data}")

        try:
            self.input_data = data.get('input_data', {})
            if self.input_data:
                self.update_data(self.input_data)
                logging.debug("Applied input_data to VariableSetterUI")
                for ptm in self.input_data.get('ptm', []):
                    logging.debug(f"PTM dataset in set_current_data: {ptm}")
                    if not ptm.get('modulation_column'):
                        logging.warning(f"Missing modulation_column in PTM dataset: {ptm.get('file_path')}")
            else:
                logging.warning("No input_data found in loaded data")

            settings = data.get('settings', {})
            if not settings:
                logging.warning("No settings found in loaded data")
                return

            for key, value in settings.items():
                if hasattr(self, key):
                    widget = getattr(self, key)
                    try:
                        if isinstance(widget, QLineEdit):
                            widget.setText(str(value))
                        elif isinstance(widget, QSpinBox):
                            widget.setValue(int(float(value)))
                        elif isinstance(widget, QDoubleSpinBox):
                            widget.setValue(float(value))
                        elif isinstance(widget, QComboBox):
                            widget.setCurrentText(str(value))
                        elif isinstance(widget, QCheckBox):
                            widget.setChecked(value.lower() == 'true' if isinstance(value, str) else bool(value))
                        elif key.endswith('_color') and isinstance(value, (list, tuple)) and len(value) == 3:
                            color = QColor(*value)
                            setattr(self, key, color)
                            if key in ['positive_color', 'negative_color']:
                                color_button = getattr(self, f"{key}_button", None)
                                if color_button:
                                    color_button.setStyleSheet(
                                        f"QPushButton {{ background-color: {color.name()}; "
                                        "border: 2px solid #555555; "
                                        "border-radius: 5px; }}"
                                        "QPushButton:hover {{ border: 2px solid #333333; }}"
                                        "QPushButton:pressed {{ border: 2px solid #111111; }}"
                                    )
                            if key in ['positive_color', 'negative_color']:
                                self.gradient_widget.update_colors(
                                    self.positive_color,
                                    self.negative_color
                                )
                        logging.debug(f"Applied setting {key}: {value}")
                    except (ValueError, TypeError) as e:
                        logging.warning(f"Failed to apply setting {key} with value {value}: {str(e)}")
                else:
                    logging.warning(f"UI element {key} not found in VariableSetterUI")

            if hasattr(self, 'results_list'):
                self.results_list.hide()
            if hasattr(self, 'run_button'):
                self.run_button.setFocus()
                logging.debug("Set focus to run_button and hid results_list")

        except Exception as e:
            logging.error(f"Error in set_current_data: {str(e)}")
            print(f"Error in set_current_data: {str(e)}")
            raise


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = VariableSetterUI(use_splitter=True)
    window.show()
    sys.exit(app.exec_())