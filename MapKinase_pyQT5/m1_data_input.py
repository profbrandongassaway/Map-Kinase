import sys
import os
import logging
import pandas as pd
import re
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QFrame, QComboBox, QPushButton,
    QLabel, QScrollArea, QFileDialog, QDialog, QListWidget, QGridLayout, QMessageBox, QLineEdit,
    QRadioButton, QButtonGroup, QGroupBox, QSpinBox, QColorDialog, QTabWidget
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor

# Helper function to get resource path
def resource_path(relative_path):
    try:
        if hasattr(sys, '_MEIPASS'):
            return os.path.join(sys._MEIPASS, relative_path)
        return os.path.join(os.path.abspath("."), relative_path)
    except Exception as e:
        logging.error(f"Error resolving resource path for {relative_path}: {str(e)}")
        return relative_path

# Get user config directory
def get_user_config_dir():
    return os.path.join(os.path.expanduser("~"), ".mapkinase")

# Find the best header match
def find_best_header_match(file_path: str, search_strings: list, parent=None) -> str:
    try:
        df = pd.read_csv(file_path, sep="\t", nrows=1)
        headers = df.columns.tolist()
    except (FileNotFoundError, PermissionError) as e:
        QMessageBox.critical(parent, "File Error", f"Cannot access file {file_path}: {str(e)}")
        return ""
    except pd.errors.ParserError as e:
        QMessageBox.critical(parent, "File Error", f"Invalid file format in {file_path}: {str(e)}")
        return ""
    except Exception as e:
        QMessageBox.critical(parent, "File Error", f"Error reading headers from {file_path}: {str(e)}")
        return ""

    normalized_headers = [str(h).lower().replace(" ", "_") for h in headers]
    for search_str in search_strings:
        if not search_str:
            continue
        norm_search = str(search_str).lower().replace(" ", "_")
        for norm_header, orig_header in zip(normalized_headers, headers):
            if norm_search == norm_header:
                return orig_header
        for norm_header, orig_header in zip(normalized_headers, headers):
            if norm_search in norm_header:
                return orig_header
    return ""

# NoScrollSpinBox
class NoScrollSpinBox(QSpinBox):
    def wheelEvent(self, event):
        event.ignore()

# NoScrollComboBox
class NoScrollComboBox(QComboBox):
    def wheelEvent(self, event):
        event.ignore()

class PTMTypeDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select PTM Dataset Type")
        self.setGeometry(200, 200, 300, 200)
        self.selected_type = None

        layout = QVBoxLayout(self)
        layout.setAlignment(Qt.AlignTop)

        self.button_group = QButtonGroup(self)
        ptm_types = ['Phosphorylation', 'Cysteine Oxidation', 'Metabolite', 'Acetylation', 'Ubiquitination',
                     'Methylation']
        for i, ptm_type in enumerate(ptm_types):
            radio = QRadioButton(ptm_type)
            self.button_group.addButton(radio, i)
            layout.addWidget(radio)
            if i == 0:
                radio.setChecked(True)
        confirm_button = QPushButton("Confirm")
        confirm_button.clicked.connect(self.confirm_selection)
        layout.addWidget(confirm_button)

    def confirm_selection(self):
        selected_button = self.button_group.checkedButton()
        if selected_button:
            self.selected_type = selected_button.text()
            self.accept()
        else:
            QMessageBox.critical(self, "Selection Error", "Please select a PTM type.")

class MainColumnsDialog(QDialog):
    def __init__(self, file_path, dataset_type, protein_columns=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Select Main Columns ({dataset_type})")
        self.setGeometry(200, 200, 600, 400)
        self.file_path = file_path
        self.dataset_type = dataset_type
        self.protein_columns = protein_columns or []
        self.selected_columns = []
        self.headers = []
        logging.debug(
            f"MainColumnsDialog initialized for {dataset_type} with file_path: {file_path}, protein_columns: {self.protein_columns}")

        self.layout = QVBoxLayout(self)
        self.layout.setAlignment(Qt.AlignTop)

        try:
            if dataset_type == "Protein":
                self.setup_protein_ui()
            else:
                self.setup_ptm_ui()
        except Exception as e:
            logging.error(f"Error setting up MainColumnsDialog UI for {dataset_type}: {str(e)}")
            QMessageBox.critical(self, "UI Error", f"Failed to initialize dialog: {str(e)}")
            self.reject()

    def setup_protein_ui(self):
        logging.debug(f"Setting up Protein UI with file_path: {self.file_path}")
        try:
            df = pd.read_csv(self.file_path, sep="\t", nrows=1)
            self.headers = df.columns.tolist()
        except Exception as e:
            logging.error(f"Error reading headers from {self.file_path}: {str(e)}")
            QMessageBox.critical(self, "File Error", f"Error reading headers: {str(e)}")
            self.headers = []
            self.reject()
            return

        main_frame = QFrame()
        main_layout = QHBoxLayout(main_frame)

        headers_frame = QFrame()
        headers_layout = QVBoxLayout(headers_frame)
        headers_layout.addWidget(QLabel("Available Headers"))
        self.available_list = QListWidget()
        self.available_list.setSelectionMode(QListWidget.ExtendedSelection)
        headers_layout.addWidget(self.available_list)
        main_layout.addWidget(headers_frame)

        buttons_frame = QFrame()
        buttons_layout = QVBoxLayout(buttons_frame)
        buttons_layout.addStretch()
        add_button = QPushButton("Add →")
        add_button.clicked.connect(self.move_to_selected)
        buttons_layout.addWidget(add_button)
        remove_button = QPushButton("← Remove")
        remove_button.clicked.connect(self.remove_from_selected)
        buttons_layout.addWidget(remove_button)
        buttons_layout.addStretch()
        main_layout.addWidget(buttons_frame)

        selected_frame = QFrame()
        selected_layout = QVBoxLayout(selected_frame)
        selected_layout.addWidget(QLabel("Selected Columns"))
        self.selected_list = QListWidget()
        self.selected_list.setSelectionMode(QListWidget.ExtendedSelection)
        selected_layout.addWidget(self.selected_list)
        main_layout.addWidget(selected_frame)

        self.layout.addWidget(main_frame)

        confirm_button = QPushButton("Confirm Selection")
        confirm_button.clicked.connect(self.confirm_selection)
        self.layout.addWidget(confirm_button)

        self.populate_available_list()

    def setup_ptm_ui(self):
        logging.debug(f"Setting up PTM UI with protein_columns: {self.protein_columns}")
        if not self.protein_columns:
            logging.warning("No Protein main columns available for PTM dialog")
            QMessageBox.critical(self, "Selection Error",
                                 "No Protein main columns available. Please select Protein main columns first.")
            self.reject()
            return

        try:
            df = pd.read_csv(self.file_path, sep="\t", nrows=1)
            self.headers = df.columns.tolist()
        except Exception as e:
            logging.error(f"Error reading headers from {self.file_path}: {str(e)}")
            QMessageBox.critical(self, "File Error", f"Error reading headers: {str(e)}")
            self.headers = []
            self.reject()
            return

        if not self.headers:
            logging.warning(f"No headers found in {self.file_path}")
            QMessageBox.critical(self, "File Error", "No valid headers found in the selected file.")
            self.reject()
            return

        grid_frame = QFrame()
        grid_layout = QGridLayout(grid_frame)
        grid_layout.addWidget(QLabel("Protein Main Column"), 0, 0)
        grid_layout.addWidget(QLabel("PTM Matching Column"), 0, 1)

        self.ptm_combos = {}
        for i, prot_col in enumerate(self.protein_columns, 1):
            grid_layout.addWidget(QLabel(str(prot_col)), i, 0)
            combo = NoScrollComboBox()
            combo.addItem("None")
            combo.addItems(self.headers)
            prot_col_str = str(prot_col).lower()
            for header in self.headers:
                header_lower = header.lower()
                if prot_col_str == header_lower or re.search(re.escape(prot_col_str), header_lower):
                    combo.setCurrentText(header)
                    break
            grid_layout.addWidget(combo, i, 1)
            self.ptm_combos[prot_col] = combo

        self.layout.addWidget(grid_frame)

        confirm_button = QPushButton("Confirm Selection")
        confirm_button.clicked.connect(self.confirm_ptm_selection)
        self.layout.addWidget(confirm_button)

    def populate_available_list(self):
        self.available_list.clear()
        for header in self.headers:
            self.available_list.addItem(header)

    def move_to_selected(self):
        selected_items = self.available_list.selectedItems()
        for item in selected_items:
            header = item.text()
            if header not in [self.selected_list.item(i).text() for i in range(self.selected_list.count())]:
                self.selected_list.addItem(header)
                self.available_list.takeItem(self.available_list.row(item))

    def remove_from_selected(self):
        selected_items = self.selected_list.selectedItems()
        for item in selected_items:
            header = item.text()
            self.available_list.addItem(header)
            self.selected_list.takeItem(self.selected_list.row(item))

    def confirm_selection(self):
        self.selected_columns = [self.selected_list.item(i).text() for i in range(self.selected_list.count())]
        if not self.selected_columns:
            QMessageBox.critical(self, "Selection Error", "No columns selected.")
            return
        self.accept()

    def confirm_ptm_selection(self):
        self.selected_columns = []
        for prot_col, combo in self.ptm_combos.items():
            ptm_col = combo.currentText()
            if ptm_col != "None":
                self.selected_columns.append((str(prot_col), ptm_col))
        if not self.selected_columns:
            QMessageBox.critical(self, "Selection Error", "At least one PTM column must be matched.")
            return
        unmatched = [col for col in self.protein_columns if not any(p == col for p, _ in self.selected_columns)]
        if unmatched:
            reply = QMessageBox.warning(
                self, "Warning",
                f"Not all Protein columns have matching PTM columns: {', '.join(unmatched)}. Proceed anyway?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        self.accept()

class TooltipColumnsDialog(QDialog):
    def __init__(self, file_path, dataset_type, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Select Tooltip Columns ({dataset_type})")
        self.setGeometry(200, 200, 600, 400)
        self.file_path = file_path
        self.dataset_type = dataset_type
        self.selected_columns = []
        self.headers = []
        logging.debug(
            f"TooltipColumnsDialog initialized for {dataset_type} with file_path: {file_path}")

        self.layout = QVBoxLayout(self)
        self.layout.setAlignment(Qt.AlignTop)

        try:
            self.setup_ui()
        except Exception as e:
            logging.error(f"Error setting up TooltipColumnsDialog UI for {dataset_type}: {str(e)}")
            QMessageBox.critical(self, "UI Error", f"Failed to initialize dialog: {str(e)}")
            self.reject()

    def setup_ui(self):
        logging.debug(f"Setting up TooltipColumns UI with file_path: {self.file_path}")
        try:
            df = pd.read_csv(self.file_path, sep="\t", nrows=1)
            self.headers = df.columns.tolist()
        except Exception as e:
            logging.error(f"Error reading headers from {self.file_path}: {str(e)}")
            QMessageBox.critical(self, "File Error", f"Error reading headers: {str(e)}")
            self.headers = []
            self.reject()
            return

        main_frame = QFrame()
        main_layout = QHBoxLayout(main_frame)

        headers_frame = QFrame()
        headers_layout = QVBoxLayout(headers_frame)
        headers_layout.addWidget(QLabel("Available Headers"))
        self.available_list = QListWidget()
        self.available_list.setSelectionMode(QListWidget.ExtendedSelection)
        headers_layout.addWidget(self.available_list)
        main_layout.addWidget(headers_frame)

        buttons_frame = QFrame()
        buttons_layout = QVBoxLayout(buttons_frame)
        buttons_layout.addStretch()
        add_button = QPushButton("Add →")
        add_button.clicked.connect(self.move_to_selected)
        buttons_layout.addWidget(add_button)
        remove_button = QPushButton("← Remove")
        remove_button.clicked.connect(self.remove_from_selected)
        buttons_layout.addWidget(remove_button)
        buttons_layout.addStretch()
        main_layout.addWidget(buttons_frame)

        selected_frame = QFrame()
        selected_layout = QVBoxLayout(selected_frame)
        selected_layout.addWidget(QLabel("Selected Tooltip Columns"))
        self.selected_list = QListWidget()
        self.selected_list.setSelectionMode(QListWidget.ExtendedSelection)
        selected_layout.addWidget(self.selected_list)
        main_layout.addWidget(selected_frame)

        self.layout.addWidget(main_frame)

        confirm_button = QPushButton("Confirm Selection")
        confirm_button.clicked.connect(self.confirm_selection)
        self.layout.addWidget(confirm_button)

        self.populate_available_list()

    def populate_available_list(self):
        self.available_list.clear()
        for header in self.headers:
            self.available_list.addItem(header)

    def move_to_selected(self):
        selected_items = self.available_list.selectedItems()
        for item in selected_items:
            header = item.text()
            if header not in [self.selected_list.item(i).text() for i in range(self.selected_list.count())]:
                self.selected_list.addItem(header)
                self.available_list.takeItem(self.available_list.row(item))

    def remove_from_selected(self):
        selected_items = self.selected_list.selectedItems()
        for item in selected_items:
            header = item.text()
            self.available_list.addItem(header)
            self.selected_list.takeItem(self.selected_list.row(item))

    def confirm_selection(self):
        self.selected_columns = [self.selected_list.item(i).text() for i in range(self.selected_list.count())]
        self.accept()


# Helper function to get resource path
def resource_path(relative_path):
    try:
        if hasattr(sys, '_MEIPASS'):
            return os.path.join(sys._MEIPASS, relative_path)
        return os.path.join(os.path.abspath("."), relative_path)
    except Exception as e:
        logging.error(f"Error resolving resource path for {relative_path}: {str(e)}")
        return relative_path

# Get user config directory
def get_user_config_dir():
    return os.path.join(os.path.expanduser("~"), ".mapkinase")

# Find the best header match
def find_best_header_match(file_path: str, search_strings: list, parent=None) -> str:
    try:
        df = pd.read_csv(file_path, sep="\t", nrows=1)
        headers = df.columns.tolist()
    except (FileNotFoundError, PermissionError) as e:
        QMessageBox.critical(parent, "File Error", f"Cannot access file {file_path}: {str(e)}")
        return ""
    except pd.errors.ParserError as e:
        QMessageBox.critical(parent, "File Error", f"Invalid file format in {file_path}: {str(e)}")
        return ""
    except Exception as e:
        QMessageBox.critical(parent, "File Error", f"Error reading headers from {file_path}: {str(e)}")
        return ""

    normalized_headers = [str(h).lower().replace(" ", "_") for h in headers]
    for search_str in search_strings:
        if not search_str:
            continue
        norm_search = str(search_str).lower().replace(" ", "_")
        for norm_header, orig_header in zip(normalized_headers, headers):
            if norm_search == norm_header:
                return orig_header
        for norm_header, orig_header in zip(normalized_headers, headers):
            if norm_search in norm_header:
                return orig_header
    return ""

# NoScrollSpinBox
class NoScrollSpinBox(QSpinBox):
    def wheelEvent(self, event):
        event.ignore()

# NoScrollComboBox
class NoScrollComboBox(QComboBox):
    def wheelEvent(self, event):
        event.ignore()

class SymbolListDialog(QDialog):
    MAX_SYMBOLS = 8  # Maximum number of symbol tabs

    def __init__(self, file_path, dataset_type, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Configure PTM Symbol Labels ({dataset_type})")
        self.setGeometry(200, 200, 800, 500)
        self.file_path = file_path
        self.dataset_type = dataset_type
        self.headers = []
        self.ptm_dict = {"ptm_symbol_list": []}
        logging.debug(f"SymbolListDialog initialized for {dataset_type} with file_path: {file_path}")

        self.layout = QVBoxLayout(self)
        self.setup_ui()

    def setup_ui(self):
        try:
            df = pd.read_csv(self.file_path, sep="\t", nrows=1)
            self.headers = df.columns.tolist()
            logging.debug(f"Headers loaded: {self.headers}")
        except Exception as e:
            logging.error(f"Error reading headers from {self.file_path}: {str(e)}")
            QMessageBox.critical(self, "File Error", f"Error reading headers: {str(e)}")
            self.reject()
            return

        # Layout selector
        layout_selector_frame = QFrame()
        layout_selector_layout = QHBoxLayout(layout_selector_frame)
        layout_selector_layout.addWidget(QLabel("Default Layout:"))
        self.layout_selector = NoScrollComboBox()
        self.layout_selector.addItems(["PhosphoSitePlus", "Custom"])
        self.layout_selector.setCurrentText("PhosphoSitePlus")
        layout_selector_layout.addWidget(self.layout_selector)
        layout_selector_layout.addStretch()
        self.layout.addWidget(layout_selector_frame)

        # Tab widget for symbol configurations
        self.tab_widget = QTabWidget()
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self.remove_tab)
        self.layout.addWidget(self.tab_widget)

        # Add tab button (for Custom layout)
        self.add_tab_button = QPushButton("+ Add Symbol Label")
        self.add_tab_button.setStyleSheet("""
            QPushButton {
                font-size: 12pt;
                font-weight: bold;
                border: 2px solid #555555;
                border-radius: 8px;
                background-color: #e0e0e0;
                padding: 5px;
            }
            QPushButton:hover {
                background-color: #d0d0d0;
                border: 2px solid #333333;
            }
            QPushButton:pressed {
                background-color: #c0c0c0;
                border: 2px solid #111111;
            }
        """)
        self.add_tab_button.clicked.connect(self.add_tab)
        self.layout.addWidget(self.add_tab_button)

        # Confirm button
        confirm_button = QPushButton("Confirm")
        confirm_button.clicked.connect(self.confirm_selection)
        self.layout.addWidget(confirm_button)

        # Initialize tabs based on layout
        self.layout_selector.currentTextChanged.connect(self.update_symbol_layout)
        self.update_symbol_layout("PhosphoSitePlus")

    def update_symbol_layout(self, layout_type):
        # Clear existing tabs
        while self.tab_widget.count():
            self.tab_widget.removeTab(0)
        logging.debug(f"Cleared {self.tab_widget.count()} existing tabs for layout: {layout_type}")

        if layout_type == "PhosphoSitePlus":
            reg_site_header = find_best_header_match(self.file_path, ["Regulatory site"], self)
            reg_function_header = find_best_header_match(self.file_path, ["Regulatory site function", "ON_FUNCTION"], self)
            logging.debug(f"Headers found - Regulatory site: {reg_site_header}, Regulatory site function: {reg_function_header}")
            if not reg_site_header or not reg_function_header:
                QMessageBox.warning(self, "Header Warning", "Some expected headers not found in the file. Please verify the file format.")
            # Define PhosphoSitePlus symbol settings
            symbol_settings = [
                {
                    "symbol": "↑",
                    "symbol_font": "Segoe UI Symbol",
                    "symbol_size": 6,
                    "symbol_color": (0, 0, 0),
                    "symbol_x_offset": 0,
                    "symbol_y_offset": -1,
                    "statement_type": "if_or_statement",
                    "header_to_search": reg_function_header if reg_function_header in self.headers else "None",
                    "search_text_1": "activity, induced",
                    "search_text_2": "enzymatic activity, induced",
                    "search_text_3": "",
                    "search_text_4": ""
                },
                {
                    "symbol": "x",
                    "symbol_font": "Arial",
                    "symbol_size": 6,
                    "symbol_color": (0, 0, 0),
                    "symbol_x_offset": 0,
                    "symbol_y_offset": -1,
                    "statement_type": "if_or_statement",
                    "header_to_search": reg_function_header if reg_function_header in self.headers else "None",
                    "search_text_1": "activity, inhibited",
                    "search_text_2": "enzymatic activity, inhibited",
                    "search_text_3": "",
                    "search_text_4": ""
                },
                {
                    "symbol": "+",
                    "symbol_font": "Arial",
                    "symbol_size": 6,
                    "symbol_color": (0, 0, 0),
                    "symbol_x_offset": 0,
                    "symbol_y_offset": -1,
                    "statement_type": "if_or_and_or_statement",
                    "header_to_search": reg_function_header if reg_function_header in self.headers else "None",
                    "search_text_1": "activity, induced",
                    "search_text_2": "enzymatic activity, induced",
                    "search_text_3": "activity, inhibited",
                    "search_text_4": "enzymatic activity, inhibited"
                },
                {
                    "symbol": "+",
                    "symbol_font": "Arial",
                    "symbol_size": 6,
                    "symbol_color": (0, 0, 0),
                    "symbol_x_offset": 0,
                    "symbol_y_offset": 0,
                    "statement_type": "if_and_notlabeled_statement",
                    "header_to_search": reg_site_header if reg_site_header in self.headers else "None",
                    "search_text_1": "+",
                    "search_text_2": "",
                    "search_text_3": "",
                    "search_text_4": ""
                }
            ]
            for i, setting in enumerate(symbol_settings):
                self.add_tab(symbol_dict=setting, removable=False)
                logging.debug(f"Added PhosphoSitePlus tab {i+1}: {setting}")
            self.add_tab_button.setVisible(False)
            self.tab_widget.setTabsClosable(False)
        else:
            # Custom layout: Start with one empty tab
            self.add_tab()
            self.add_tab_button.setVisible(True)
            self.tab_widget.setTabsClosable(True)
            logging.debug("Custom layout selected, added one initial tab")

    def add_tab(self, symbol_dict=None, removable=True):
        if self.tab_widget.count() >= self.MAX_SYMBOLS:
            QMessageBox.warning(self, "Limit Reached", f"Cannot add more than {self.MAX_SYMBOLS} symbols.")
            return

        # Create tab widget
        tab = QWidget()
        tab_layout = QVBoxLayout(tab)
        form_frame = QFrame()
        form_layout = QGridLayout(form_frame)
        row = 0

        # Symbol settings
        form_layout.addWidget(QLabel("Symbol:"), row, 0)
        symbol_input = QLineEdit(symbol_dict.get("symbol", "+") if symbol_dict else "+")
        symbol_input.setMaximumWidth(50)
        form_layout.addWidget(symbol_input, row, 1)
        row += 1

        form_layout.addWidget(QLabel("Font:"), row, 0)
        font_input = QLineEdit(symbol_dict.get("symbol_font", "Arial") if symbol_dict else "Arial")
        form_layout.addWidget(font_input, row, 1)
        row += 1

        form_layout.addWidget(QLabel("Size:"), row, 0)
        size_input = NoScrollSpinBox()
        size_input.setRange(1, 20)
        size_input.setValue(symbol_dict.get("symbol_size", 6) if symbol_dict else 6)
        size_input.setMaximumWidth(50)
        form_layout.addWidget(size_input, row, 1)
        row += 1

        form_layout.addWidget(QLabel("Color:"), row, 0)
        color = QColor(*symbol_dict.get("symbol_color", (0, 0, 0))) if symbol_dict else QColor(0, 0, 0)
        color_button = QPushButton("Pick")
        color_button.setMaximumWidth(60)
        color_button.setStyleSheet(f"background-color: {color.name()};")
        color_button.color = color
        color_button.clicked.connect(lambda: self.pick_color(color_button))
        form_layout.addWidget(color_button, row, 1)
        row += 1

        form_layout.addWidget(QLabel("X Offset:"), row, 0)
        x_offset_input = NoScrollSpinBox()
        x_offset_input.setRange(-20, 20)
        x_offset_input.setValue(symbol_dict.get("symbol_x_offset", 0) if symbol_dict else 0)
        x_offset_input.setMaximumWidth(50)
        form_layout.addWidget(x_offset_input, row, 1)
        row += 1

        form_layout.addWidget(QLabel("Y Offset:"), row, 0)
        y_offset_input = NoScrollSpinBox()
        y_offset_input.setRange(-20, 20)
        y_offset_input.setValue(symbol_dict.get("symbol_y_offset", 0) if symbol_dict else 0)
        y_offset_input.setMaximumWidth(50)
        form_layout.addWidget(y_offset_input, row, 1)
        row += 1

        # Statement settings
        form_layout.addWidget(QLabel("Statement Type:"), row, 0)
        statement_combo = NoScrollComboBox()
        statement_types = [
            "NA", "if_statement", "if_not_statement", "if_and_statement", "if_or_statement",
            "if_and_not_statement", "if_and_and_not_statement", "if_or_and_not_statement",
            "if_or_and_or_statement", "if_and_notlabeled_statement"
        ]
        statement_combo.addItems(statement_types)
        statement_combo.setCurrentText(symbol_dict.get("statement_type", "NA") if symbol_dict else "NA")
        form_layout.addWidget(statement_combo, row, 1)
        row += 1

        # Header to search
        header_label = QLabel("Header to Search:")
        form_layout.addWidget(header_label, row, 0)
        header_combo = NoScrollComboBox()
        header_combo.addItems(["None"] + self.headers)
        header_combo.setCurrentText(symbol_dict.get("header_to_search", "None") if symbol_dict else "None")
        form_layout.addWidget(header_combo, row, 1)
        row += 1

        # Search text fields (up to 4, hide as needed)
        self.search_text_1 = QLineEdit(symbol_dict.get("search_text_1", "") if symbol_dict else "")
        self.search_text_2 = QLineEdit(symbol_dict.get("search_text_2", "") if symbol_dict else "")
        self.search_text_3 = QLineEdit(symbol_dict.get("search_text_3", "") if symbol_dict else "")
        self.search_text_4 = QLineEdit(symbol_dict.get("search_text_4", "") if symbol_dict else "")
        search_labels = [QLabel("Search Text 1:"), QLabel("Search Text 2:"), QLabel("Search Text 3:"), QLabel("Search Text 4:")]
        self.search_widgets = [(lbl, txt) for lbl, txt in zip(search_labels, [self.search_text_1, self.search_text_2, self.search_text_3, self.search_text_4])]

        for i, (lbl, txt) in enumerate(self.search_widgets):
            form_layout.addWidget(lbl, row, 0)
            form_layout.addWidget(txt, row, 1)
            row += 1

        # NA message
        self.na_label = QLabel("No search statement. Symbol will not be applied.")
        form_layout.addWidget(self.na_label, row, 0, 1, 2)
        row += 1

        tab_layout.addWidget(form_frame)
        tab_layout.addStretch()

        # Store widgets for later access
        tab.symbol_input = symbol_input
        tab.font_input = font_input
        tab.size_input = size_input
        tab.color_button = color_button
        tab.x_offset_input = x_offset_input
        tab.y_offset_input = y_offset_input
        tab.statement_combo = statement_combo
        tab.header_combo = header_combo
        tab.search_text_1 = self.search_text_1
        tab.search_text_2 = self.search_text_2
        tab.search_text_3 = self.search_text_3
        tab.search_text_4 = self.search_text_4
        tab.na_label = self.na_label

        # Update visibility based on statement type
        def update_visibility():
            statement_type = statement_combo.currentText()
            logging.debug(f"Updating visibility for statement_type: {statement_type}")

            # Hide all search fields and header by default
            header_label.setVisible(statement_type != "NA")
            header_combo.setVisible(statement_type != "NA")
            for lbl, txt in self.search_widgets:
                lbl.setVisible(False)
                txt.setVisible(False)
            self.na_label.setVisible(statement_type == "NA")

            # Show relevant fields based on statement type
            if statement_type in ["if_statement", "if_not_statement", "if_and_notlabeled_statement"]:
                self.search_widgets[0][0].setVisible(True)
                self.search_widgets[0][1].setVisible(True)
            elif statement_type in ["if_and_statement", "if_or_statement", "if_and_not_statement"]:
                self.search_widgets[0][0].setVisible(True)
                self.search_widgets[0][1].setVisible(True)
                self.search_widgets[1][0].setVisible(True)
                self.search_widgets[1][1].setVisible(True)
            elif statement_type in ["if_and_and_not_statement", "if_or_and_not_statement"]:
                self.search_widgets[0][0].setVisible(True)
                self.search_widgets[0][1].setVisible(True)
                self.search_widgets[1][0].setVisible(True)
                self.search_widgets[1][1].setVisible(True)
                self.search_widgets[2][0].setVisible(True)
                self.search_widgets[2][1].setVisible(True)
            elif statement_type == "if_or_and_or_statement":
                for i in range(4):
                    self.search_widgets[i][0].setVisible(True)
                    self.search_widgets[i][1].setVisible(True)

        # Block signals during initialization
        statement_combo.blockSignals(True)
        update_visibility()
        statement_combo.blockSignals(False)

        # Connect statement type change
        statement_combo.currentTextChanged.connect(update_visibility)

        # Add tab to widget
        tab_index = self.tab_widget.addTab(tab, f"Symbol {self.tab_widget.count() + 1}")
        self.tab_widget.setCurrentIndex(tab_index)
        tab.removable = removable
        logging.debug(f"Added tab {tab_index + 1}, removable: {removable}, settings: {symbol_dict}")

        # Update add button state
        self.add_tab_button.setEnabled(self.tab_widget.count() < self.MAX_SYMBOLS)

    def remove_tab(self, index):
        if not self.tab_widget.widget(index).removable:
            QMessageBox.warning(self, "Cannot Remove", "PhosphoSitePlus symbols cannot be removed.")
            return
        self.tab_widget.removeTab(index)
        self.add_tab_button.setEnabled(self.tab_widget.count() < self.MAX_SYMBOLS)
        logging.debug(f"Removed tab at index {index}, remaining tabs: {self.tab_widget.count()}")

    def pick_color(self, button):
        color_dialog = QColorDialog(button.color, self)
        if color_dialog.exec_():
            new_color = color_dialog.selectedColor()
            button.setStyleSheet(f"background-color: {new_color.name()};")
            button.color = new_color

    def confirm_selection(self):
        if not self.tab_widget.count():
            logging.error("No symbol tabs defined")
            QMessageBox.critical(self, "Error", "At least one symbol must be defined.")
            return

        symbol_list = []
        for i in range(self.tab_widget.count()):
            tab = self.tab_widget.widget(i)
            symbol_dict = {}
            try:
                # Symbol settings
                symbol_dict.update({
                    "symbol": tab.symbol_input.text(),
                    "symbol_font": tab.font_input.text(),
                    "symbol_size": tab.size_input.value(),
                    "symbol_color": (tab.color_button.color.red(), tab.color_button.color.green(), tab.color_button.color.blue()),
                    "symbol_x_offset": tab.x_offset_input.value(),
                    "symbol_y_offset": tab.y_offset_input.value(),
                    "statement_type": tab.statement_combo.currentText(),
                    "header_to_search": tab.header_combo.currentText() if tab.header_combo.currentText() != "None" else "",
                    "search_text_1": tab.search_text_1.text(),
                    "search_text_2": tab.search_text_2.text(),
                    "search_text_3": tab.search_text_3.text(),
                    "search_text_4": tab.search_text_4.text()
                })
                symbol_list.append(symbol_dict)
                logging.debug(f"Tab {i+1}: Processed symbol_dict: {symbol_dict}")
            except Exception as e:
                logging.error(f"Tab {i+1}: Error processing symbol: {str(e)}")
                continue

        if not symbol_list:
            logging.error("No valid symbols defined after processing tabs")
            QMessageBox.critical(self, "Error", "At least one symbol must be defined.")
            return

        self.ptm_dict = {"ptm_symbol_list": [{f"symbol_label_{i+1}_dict": d} for i, d in enumerate(symbol_list)]}
        logging.debug(f"Final PTM dict: {self.ptm_dict}")
        self.accept()


class DatasetSelector(QFrame):
    custom_colors_file = os.path.join(get_user_config_dir(), "custom_colors.txt")

    def __init__(self, dataset_type, removable=False, parent=None, data_selection_ui=None):
        super().__init__(parent)
        self.dataset_type = dataset_type
        self.data_selection_ui = data_selection_ui
        self.setFrameStyle(QFrame.Box | QFrame.Raised)
        self.setLineWidth(1)
        self.setFixedWidth(450)

        self.file_path = ""
        self.headers = []
        self.main_columns = []
        self.tooltip_columns = []
        self.ptm_symbol_list = []

        layout = QVBoxLayout()

        if removable:
            remove_btn = QPushButton("X")
            remove_btn.setFixedSize(20, 20)
            remove_btn.setStyleSheet("background-color: red; color: white;")
            remove_btn.clicked.connect(self.remove_selector)
            btn_layout = QHBoxLayout()
            btn_layout.addStretch()
            btn_layout.addWidget(remove_btn)
            layout.addLayout(btn_layout)

        label = QLabel(f"{dataset_type} Dataset")
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)

        input_layout = QHBoxLayout()
        self.input_text = QLineEdit()
        self.input_text.setPlaceholderText("Select a file...")
        self.input_text.setReadOnly(True)
        input_layout.addWidget(self.input_text)
        browse_btn = QPushButton("Browse")
        browse_btn.setMaximumWidth(80)
        browse_btn.clicked.connect(self.browse_file)
        input_layout.addWidget(browse_btn)
        layout.addLayout(input_layout)

        self.main_columns_btn = QPushButton("Select Main Columns")
        self.main_columns_btn.clicked.connect(self.select_main_columns)
        if dataset_type != "Protein":
            parent_columns = self.data_selection_ui.protein_main_columns if self.data_selection_ui else []
            self.main_columns_btn.setEnabled(bool(parent_columns))
        layout.addWidget(self.main_columns_btn)

        self.main_columns_label = QLabel("No main columns selected")
        self.main_columns_label.setWordWrap(True)
        layout.addWidget(self.main_columns_label)

        self.tooltip_columns_btn = QPushButton("Select Tooltip Columns")
        self.tooltip_columns_btn.clicked.connect(self.select_tooltip_columns)
        layout.addWidget(self.tooltip_columns_btn)

        self.tooltip_columns_label = QLabel("No tooltip columns selected")
        self.tooltip_columns_label.setWordWrap(True)
        layout.addWidget(self.tooltip_columns_label)

        self.uniprot_combo = NoScrollComboBox()
        self.uniprot_combo.addItem("Select UniProt ID column")
        layout.addWidget(QLabel("UniProt ID Column"))
        layout.addWidget(self.uniprot_combo)

        if dataset_type == "Protein":
            self.kegg_combo = NoScrollComboBox()
            self.kegg_combo.addItem("Select KEGG Gene ID column")
            layout.addWidget(QLabel("KEGG Gene ID Column"))
            layout.addWidget(self.kegg_combo)

            self.gene_combo = NoScrollComboBox()
            self.gene_combo.addItem("Select Gene column")
            layout.addWidget(QLabel("Gene Column"))
            layout.addWidget(self.gene_combo)
        else:
            self.site_combo = NoScrollComboBox()
            self.site_combo.addItem("Select Site ID column")
            layout.addWidget(QLabel("Site ID Column"))
            layout.addWidget(self.site_combo)

            self.shape_combo = NoScrollComboBox()
            self.shape_combo.addItems(["Circle", "Diamond", "Triangle"])
            layout.addWidget(QLabel("PTM Shape"))
            layout.addWidget(self.shape_combo)

            self.modulation_combo = NoScrollComboBox()
            self.modulation_combo.addItem("Select Modulation column")
            layout.addWidget(QLabel("Modulation Column"))
            layout.addWidget(self.modulation_combo)

            self.symbol_list_btn = QPushButton("Configure Symbol Labels")
            self.symbol_list_btn.clicked.connect(self.configure_symbol_list)
            layout.addWidget(self.symbol_list_btn)

            self.symbol_list_label = QLabel("No symbols configured")
            self.symbol_list_label.setWordWrap(True)
            layout.addWidget(self.symbol_list_label)

        layout.addStretch()
        self.setLayout(layout)

        # Connect signals to update data on change
        self.uniprot_combo.currentTextChanged.connect(self.data_selection_ui.update_data)
        if dataset_type == "Protein":
            self.kegg_combo.currentTextChanged.connect(self.data_selection_ui.update_data)
            self.gene_combo.currentTextChanged.connect(self.data_selection_ui.update_data)
        else:
            self.site_combo.currentTextChanged.connect(self.data_selection_ui.update_data)
            self.shape_combo.currentTextChanged.connect(self.data_selection_ui.update_data)
            self.modulation_combo.currentTextChanged.connect(self.data_selection_ui.update_data)

    def configure_symbol_list(self):
        if not self.file_path or not os.path.exists(self.file_path):
            QMessageBox.critical(self, "File Error", "Please select a valid file first.")
            return
        dialog = SymbolListDialog(self.file_path, self.dataset_type, self)
        if dialog.exec_():
            self.ptm_symbol_list = dialog.ptm_dict.get('ptm_symbol_list', [])
            if self.ptm_symbol_list:
                symbol_texts = []
                for i, symbol_item in enumerate(self.ptm_symbol_list):
                    symbol_dict = symbol_item.get(f"symbol_label_{i+1}_dict", {})
                    symbol = symbol_dict.get("symbol", "N/A")
                    statement_type = symbol_dict.get("statement_type", "N/A")
                    header = symbol_dict.get("header_to_search", "None")
                    symbol_texts.append(f"{symbol}\t{statement_type}\t{header}")
                self.symbol_list_label.setText("\n".join(symbol_texts))
            else:
                self.symbol_list_label.setText("No symbols configured")
            self.data_selection_ui.update_data()

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
            QMessageBox.critical(self, "File Error", f"Failed to read custom colors file: {str(e)}")
            return
        for i, line in enumerate(lines[:16]):
            try:
                r, g, b = map(int, line.strip().split(','))
                if not (0 <= r <= 255 and 0 <= g <= 255 and 0 <= b <= 255):
                    raise ValueError("Color values out of range")
                color = QColor(r, g, b)
                color_dialog.setCustomColor(i, color)
            except ValueError as e:
                QMessageBox.critical(self, "Color Parse Error", f"Invalid color format in line {i + 1}: {str(e)}")
                continue

    def save_custom_colors(self, color_dialog):
        try:
            os.makedirs(os.path.dirname(self.custom_colors_file), exist_ok=True)
            with open(self.custom_colors_file, 'w', encoding='utf-8') as f:
                for i in range(16):
                    color = color_dialog.customColor(i)
                    f.write(f"{color.red()},{color.green()},{color.blue()}\n")
        except (IOError, PermissionError) as e:
            QMessageBox.critical(self, "File Error", f"Failed to save custom colors: {str(e)}")

    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select File", "", "Text Files (*.txt);;All Files (*)"
        )
        if file_path:
            self.file_path = file_path
            self.input_text.setText(file_path)
            self.update_columns()
            self.data_selection_ui.update_data()

    def update_columns(self):
        if not os.path.exists(self.file_path):
            QMessageBox.critical(self, "File Error", "Selected file does not exist.")
            return
        try:
            df = pd.read_csv(self.file_path, sep="\t", nrows=1)
            self.headers = df.columns.tolist()
        except Exception as e:
            logging.error(f"Error loading file {self.file_path}: {str(e)}")
            QMessageBox.critical(self, "File Error", f"Error loading file: {str(e)}")
            return

        self.uniprot_combo.clear()
        self.uniprot_combo.addItems(self.headers)
        uniprot_match = find_best_header_match(self.file_path, ["UniProt_ID", "UniProt"], self)
        if uniprot_match:
            self.uniprot_combo.setCurrentText(uniprot_match)

        if self.dataset_type == "Protein":
            self.kegg_combo.clear()
            self.kegg_combo.addItems(self.headers)
            kegg_match = find_best_header_match(self.file_path, ["KEGG_Gene_ID", "Kegg_hsa"], self)
            if kegg_match:
                self.kegg_combo.setCurrentText(kegg_match)

            self.gene_combo.clear()
            self.gene_combo.addItems(self.headers)
            gene_match = find_best_header_match(
                self.file_path, ["gene_symbol", "gene_name", "gene_id", "gene"], self
            )
            if gene_match:
                self.gene_combo.setCurrentText(gene_match)
        else:
            self.site_combo.clear()
            self.site_combo.addItems(self.headers)
            site_match = find_best_header_match(
                self.file_path, ["site_position", "site_id", "phosphosite"], self
            )
            if site_match:
                self.site_combo.setCurrentText(site_match)

            self.modulation_combo.clear()
            self.modulation_combo.addItems(self.headers)
            modulation_match = find_best_header_match(
                self.file_path, ["Regulatory site", "modulation"], self
            )
            if modulation_match:
                self.modulation_combo.setCurrentText(modulation_match)

    def update_ptm_buttons(self):
        if not self.data_selection_ui:
            logging.warning("No DataSelectionUI reference for update_ptm_buttons")
            return
        try:
            protein_main_columns = self.data_selection_ui.protein_main_columns
            logging.debug(f"Updating PTM buttons with protein_main_columns: {protein_main_columns}")
            for selector in self.data_selection_ui.selector_widget.findChildren(DatasetSelector):
                if selector.dataset_type != "Protein":
                    logging.debug(f"Enabling main_columns_btn for {selector.dataset_type}")
                    selector.main_columns_btn.setEnabled(bool(protein_main_columns))
        except AttributeError as e:
            logging.error(f"Error accessing protein_main_columns in update_ptm_buttons: {str(e)}")

    def select_main_columns(self):
        logging.debug(f"select_main_columns called for {self.dataset_type}")
        if not self.file_path or not os.path.exists(self.file_path):
            logging.warning(f"No valid file selected for {self.dataset_type}: {self.file_path}")
            QMessageBox.critical(self, "File Error", "Please select a valid file first.")
            return

        protein_main_columns = []
        if self.data_selection_ui:
            try:
                protein_main_columns = self.data_selection_ui.protein_main_columns
            except AttributeError as e:
                logging.error(f"Error accessing data_selection_ui.protein_main_columns: {str(e)}")
                QMessageBox.critical(self, "Internal Error", "Failed to access Protein main columns.")
                return
        else:
            logging.warning("No DataSelectionUI reference for select_main_columns")

        logging.debug(
            f"Opening MainColumnsDialog for {self.dataset_type} with protein_main_columns: {protein_main_columns}")
        if self.dataset_type != "Protein" and not protein_main_columns:
            logging.warning("No Protein main columns selected for PTM dataset")
            QMessageBox.critical(self, "Selection Error", "Please select Protein main columns first.")
            return

        try:
            dialog = MainColumnsDialog(
                self.file_path,
                self.dataset_type,
                protein_main_columns if self.dataset_type != "Protein" else None,
                self
            )
        except Exception as e:
            logging.error(f"Error creating MainColumnsDialog for {self.dataset_type}: {str(e)}")
            QMessageBox.critical(self, "Dialog Error", f"Failed to open column selection dialog: {str(e)}")
            return

        if dialog.exec_():
            self.main_columns = dialog.selected_columns
            if self.dataset_type == "Protein":
                try:
                    if self.data_selection_ui:
                        self.data_selection_ui.protein_main_columns = self.main_columns
                    display_text = ", ".join(str(col) for col in self.main_columns)
                    self.main_columns_label.setText(f"Main Columns: {display_text}")
                    logging.debug(f"Protein main columns selected: {self.main_columns}")
                    self.update_ptm_buttons()
                except AttributeError as e:
                    logging.error(f"Error updating protein_main_columns for Protein dataset: {str(e)}")
                    QMessageBox.critical(self, "Internal Error", "Failed to save Protein main columns.")
            else:
                display_text = ", ".join(f"{p} → {t}" for p, t in self.main_columns)
                self.main_columns_label.setText(f"Main Columns: {display_text}")
                logging.debug(f"PTM main columns selected: {self.main_columns}")
            self.data_selection_ui.update_data()

    def select_tooltip_columns(self):
        logging.debug(f"select_tooltip_columns called for {self.dataset_type}")
        if not self.file_path or not os.path.exists(self.file_path):
            logging.warning(f"No valid file selected for {self.dataset_type}: {self.file_path}")
            QMessageBox.critical(self, "File Error", "Please select a valid file first.")
            return

        try:
            dialog = TooltipColumnsDialog(self.file_path, self.dataset_type, self)
        except Exception as e:
            logging.error(f"Error creating TooltipColumnsDialog for {self.dataset_type}: {str(e)}")
            QMessageBox.critical(self, "Dialog Error", f"Failed to open tooltip column selection dialog: {str(e)}")
            return

        if dialog.exec_():
            self.tooltip_columns = dialog.selected_columns
            display_text = ", ".join(str(col) for col in self.tooltip_columns) if self.tooltip_columns else "No tooltip columns selected"
            self.tooltip_columns_label.setText(f"Tooltip Columns: {display_text}")
            logging.debug(f"Tooltip columns selected for {self.dataset_type}: {self.tooltip_columns}")
            self.data_selection_ui.update_data()

    def remove_selector(self):
        logging.debug(f"Removing {self.dataset_type} dataset selector")
        self.setParent(None)
        self.deleteLater()
        self.data_selection_ui.update_data()

class AddPTMSelector(QWidget):
    def __init__(self, add_callback, parent=None):
        super().__init__(parent)
        self.add_callback = add_callback
        layout = QHBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(5)

        layout.addStretch()
        self.add_button = QPushButton("+")
        self.add_button.setMinimumSize(40, 40)
        self.add_button.setStyleSheet("""
            QPushButton {
                font-size: 16pt;
                font-weight: bold;
                border: 2px solid #555555;
                border-radius: 8px;
                background-color: #e0e0e0;
                padding: 5px;
            }
            QPushButton:hover {
                background-color: #d0d0d0;
                border: 2px solid #333333;
            }
            QPushButton:pressed {
                background-color: #c0c0c0;
                border: 2px solid #111111;
            }
        """)
        self.add_button.clicked.connect(self.show_ptm_type_dialog)
        layout.addWidget(self.add_button, alignment=Qt.AlignCenter)
        layout.addStretch()

        self.setMinimumWidth(400)
        logging.debug(f"AddPTMSelector initialized, button size: {self.add_button.size()}, widget size: {self.size()}")

    def show_ptm_type_dialog(self):
        logging.debug("Entering AddPTMSelector.show_ptm_type_dialog")
        try:
            dialog = PTMTypeDialog(self)
            logging.debug("PTMTypeDialog created")
            if dialog.exec_() and dialog.selected_type:
                logging.debug(f"PTM type selected: {dialog.selected_type}")
                self.add_callback(dialog.selected_type)
                logging.debug("add_callback executed")
            else:
                logging.debug("PTMTypeDialog cancelled or no type selected")
        except Exception as e:
            logging.error(f"Error in show_ptm_type_dialog: {str(e)}")
        logging.debug("Exiting AddPTMSelector.show_ptm_type_dialog")

class DataSelectionUI(QWidget):
    def __init__(self, input_data=None, parent=None):
        super().__init__(parent)
        self.input_data = input_data or {'protein': None, 'ptm': [], 'visualization_settings': {}}
        self.protein_main_columns = []
        self.visualization_settings = {
            'prot_label_font': 'Arial',
            'prot_label_size': 6,
            'ptm_label_font': 'Arial',
            'ptm_label_size': 4,
            'ptm_label_color': (0, 0, 0),
            'ptm_circle_radius': 5
        }
        self.combined_ui_parent = parent
        logging.info("Initializing DataSelectionUI")

        main_layout = QVBoxLayout()

        vis_group = QGroupBox("Visualization Settings")
        vis_layout = QHBoxLayout()
        vis_layout.setSpacing(10)

        self.prot_label_font = QLineEdit("Arial")
        self.prot_label_font.setMaximumWidth(80)
        self.prot_label_size = NoScrollSpinBox()
        self.prot_label_size.setRange(1, 20)
        self.prot_label_size.setValue(6)
        self.prot_label_size.setMaximumWidth(40)
        vis_layout.addWidget(QLabel("Prot Label Font:"))
        vis_layout.addWidget(self.prot_label_font)
        vis_layout.addWidget(QLabel("Size:"))
        vis_layout.addWidget(self.prot_label_size)

        self.ptm_label_font = QLineEdit("Arial")
        self.ptm_label_font.setMaximumWidth(80)
        self.ptm_label_color = QColor(0, 0, 0)
        self.ptm_label_color_button = QPushButton("Pick")
        self.ptm_label_color_button.setMaximumWidth(60)
        self.ptm_label_color_button.clicked.connect(self.pick_ptm_label_color)
        self.ptm_label_color_square = QLabel()
        self.ptm_label_color_square.setFixedSize(20, 20)
        self.ptm_label_color_square.setStyleSheet(
            f"background-color: {self.ptm_label_color.name()}; border: 1px solid gray;")
        self.ptm_label_size = NoScrollSpinBox()
        self.ptm_label_size.setRange(1, 20)
        self.ptm_label_size.setValue(4)
        self.ptm_label_size.setMaximumWidth(40)
        vis_layout.addWidget(QLabel("PTM Label Font:"))
        vis_layout.addWidget(self.ptm_label_font)
        vis_layout.addWidget(QLabel("Color:"))
        vis_layout.addWidget(self.ptm_label_color_button)
        vis_layout.addWidget(self.ptm_label_color_square)
        vis_layout.addWidget(QLabel("Size:"))
        vis_layout.addWidget(self.ptm_label_size)

        self.ptm_circle_radius = NoScrollSpinBox()
        self.ptm_circle_radius.setRange(1, 20)
        self.ptm_circle_radius.setValue(5)
        self.ptm_circle_radius.setMaximumWidth(40)
        vis_layout.addWidget(QLabel("PTM Circle Radius:"))
        vis_layout.addWidget(self.ptm_circle_radius)

        vis_layout.addStretch()
        vis_group.setLayout(vis_layout)
        main_layout.addWidget(vis_group)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.selector_widget = QWidget()
        self.selector_layout = QHBoxLayout(self.selector_widget)
        self.selector_layout.setAlignment(Qt.AlignLeft)

        self.protein_selector = DatasetSelector("Protein", data_selection_ui=self)
        self.selector_layout.addWidget(self.protein_selector)

        if self.input_data.get('ptm'):
            for ptm_data in self.input_data['ptm']:
                ptm_selector = DatasetSelector(ptm_data.get('type', 'Phosphorylation'), removable=True,
                                               data_selection_ui=self)
                self.selector_layout.addWidget(ptm_selector)
        else:
            self.ptm_selector = DatasetSelector("Phosphorylation", removable=True, data_selection_ui=self)
            self.selector_layout.addWidget(self.ptm_selector)

        self.add_ptm_selector = AddPTMSelector(self.add_ptm_selector)
        self.selector_layout.addWidget(self.add_ptm_selector)
        self.selector_layout.addStretch()

        scroll_area.setWidget(self.selector_widget)
        main_layout.addWidget(scroll_area)

        self.setLayout(main_layout)
        logging.debug("DataSelectionUI layout initialized")

        # Connect visualization settings signals
        self.prot_label_font.textChanged.connect(self.update_data)
        self.prot_label_size.valueChanged.connect(self.update_data)
        self.ptm_label_font.textChanged.connect(self.update_data)
        self.ptm_label_size.valueChanged.connect(self.update_data)
        self.ptm_circle_radius.valueChanged.connect(self.update_data)
        self.ptm_label_color_button.clicked.connect(self.update_data_after_color)

        self.populate_existing_data()

    def get_current_data(self):
        logging.debug("DataSelectionUI.get_current_data called")
        try:
            data = self.collect_data()
            if not data:
                data = {
                    "protein": None,
                    "ptm": [],
                    "visualization_settings": {
                        'prot_label_font': self.prot_label_font.text(),
                        'prot_label_size': self.prot_label_size.value(),
                        'ptm_label_font': self.ptm_label_font.text(),
                        'ptm_label_size': self.ptm_label_size.value(),
                        'ptm_label_color': (
                            self.ptm_label_color.red(),
                            self.ptm_label_color.green(),
                            self.ptm_label_color.blue()
                        ),
                        'ptm_circle_radius': self.ptm_circle_radius.value()
                    }
                }
            logging.debug(f"Collected data: {data}")
            return data
        except Exception as e:
            logging.error(f"Error in get_current_data: {str(e)}")
            return {'error': str(e)}

    def set_current_data(self, data):
        logging.debug(f"DataSelectionUI.set_current_data called with data: {data}")
        try:
            for selector in self.selector_widget.findChildren(DatasetSelector):
                if selector.dataset_type != "Protein":
                    selector.remove_selector()

            vis_settings = data.get('visualization_settings', {})
            self.prot_label_font.setText(vis_settings.get('prot_label_font', 'Arial'))
            self.ptm_label_font.setText(vis_settings.get('ptm_label_font', 'Arial'))
            self.ptm_label_size.setValue(vis_settings.get('ptm_label_size', 4))
            ptm_label_color = vis_settings.get('ptm_label_color', (0, 0, 0))
            self.ptm_label_color = QColor(*ptm_label_color)
            self.ptm_label_color_square.setStyleSheet(
                f"background-color: {self.ptm_label_color.name()}; border: 1px solid gray;"
            )
            self.ptm_circle_radius.setValue(vis_settings.get('ptm_circle_radius', 5))
            self.visualization_settings.update({
                'prot_label_font': self.prot_label_font.text(),
                'prot_label_size': self.prot_label_size.value(),
                'ptm_label_font': self.ptm_label_font.text(),
                'ptm_label_size': self.ptm_label_size.value(),
                'ptm_label_color': (
                    self.ptm_label_color.red(),
                    self.ptm_label_color.green(),
                    self.ptm_label_color.blue()
                ),
                'ptm_circle_radius': self.ptm_circle_radius.value(),
                'protein_tooltip_columns': vis_settings.get('protein_tooltip_columns', [])
            })

            protein_data = data.get('protein')
            if protein_data and protein_data.get('file_path'):
                selector = self.protein_selector
                selector.file_path = protein_data['file_path']
                selector.input_text.setText(protein_data['file_path'])
                selector.update_columns()
                selector.main_columns = protein_data.get('main_columns', [])
                selector.tooltip_columns = vis_settings.get('protein_tooltip_columns', [])
                self.protein_main_columns = selector.main_columns
                display_text = ", ".join(str(col) for col in selector.main_columns)
                selector.main_columns_label.setText(f"Main Columns: {display_text}")
                tooltip_display_text = ", ".join(str(col) for col in
                                                 selector.tooltip_columns) if selector.tooltip_columns else "No tooltip columns selected"
                selector.tooltip_columns_label.setText(f"Tooltip Columns: {tooltip_display_text}")
                selector.uniprot_combo.setCurrentText(protein_data.get('uniprot_column', ''))
                selector.kegg_combo.setCurrentText(protein_data.get('kegg_column', ''))
                selector.gene_combo.setCurrentText(protein_data.get('gene_column', ''))
                selector.update_ptm_buttons()

            for i, ptm_data in enumerate(data.get('ptm', [])):
                selector = DatasetSelector(
                    ptm_data.get('type', 'Phosphorylation'),
                    removable=True,
                    data_selection_ui=self
                )
                self.selector_layout.insertWidget(
                    self.selector_layout.indexOf(self.add_ptm_selector), selector
                )
                selector.file_path = ptm_data.get('file_path', '')
                selector.input_text.setText(ptm_data['file_path'])
                selector.update_columns()
                selector.main_columns = [(str(p), str(t)) for p, t in ptm_data.get('main_columns', [])]
                selector.tooltip_columns = ptm_data.get('tooltip_columns', [])
                display_text = ", ".join(f"{p} → {t}" for p, t in selector.main_columns)
                selector.main_columns_label.setText(f"Main Columns: {display_text}")
                tooltip_display_text = ", ".join(str(col) for col in
                                                 selector.tooltip_columns) if selector.tooltip_columns else "No tooltip columns selected"
                selector.tooltip_columns_label.setText(f"Tooltip Columns: {tooltip_display_text}")
                selector.uniprot_combo.setCurrentText(ptm_data.get('uniprot_column', ''))
                selector.site_combo.setCurrentText(ptm_data.get('site_column', ''))
                selector.shape_combo.setCurrentText(ptm_data.get('shape', 'Circle'))
                selector.modulation_combo.setCurrentText(ptm_data.get('modulation_column', ''))
                selector.ptm_symbol_list = ptm_data.get('ptm_symbol_list', [])
                logging.debug(f"Set PTM symbol list for {selector.dataset_type}: {selector.ptm_symbol_list}")

            self.update_data()
        except Exception as e:
            logging.error(f"Error in set_current_data: {str(e)}")
            raise

    def update_data_after_color(self):
        self.pick_ptm_label_color()
        self.update_data()

    def pick_ptm_label_color(self):
        try:
            color_dialog = QColorDialog(self.ptm_label_color, self)
            self.load_custom_colors(color_dialog)
            if color_dialog.exec_():
                self.ptm_label_color = color_dialog.selectedColor()
                self.ptm_label_color_square.setStyleSheet(
                    f"background-color: {self.ptm_label_color.name()}; border: 1px solid gray;")
                self.save_custom_colors(color_dialog)
        except RuntimeError as e:
            QMessageBox.critical(self, "GUI Error", f"Error in color picker: {str(e)}")

    def load_custom_colors(self, color_dialog):
        if not os.path.exists(DatasetSelector.custom_colors_file):
            bundled_colors_file = resource_path("Scripts/custom_colors.txt")
            if os.path.exists(bundled_colors_file):
                os.makedirs(os.path.dirname(DatasetSelector.custom_colors_file), exist_ok=True)
                import shutil
                shutil.copy(bundled_colors_file, DatasetSelector.custom_colors_file)
            else:
                return
        try:
            with open(DatasetSelector.custom_colors_file, 'r') as f:
                lines = f.readlines()
        except (IOError, PermissionError) as e:
            QMessageBox.critical(self, "File Error", f"Failed to read custom colors file: {str(e)}")
            return
        for i, line in enumerate(lines[:16]):
            try:
                r, g, b = map(int, line.strip().split(','))
                if not (0 <= r <= 255 and 0 <= g <= 255 and 0 <= b <= 255):
                    raise ValueError("Color values out of range")
                color = QColor(r, g, b)
                color_dialog.setCustomColor(i, color)
            except (ValueError, IndexError) as e:
                QMessageBox.critical(self, "Color Parse Error", f"Invalid color format in line {i + 1}: {str(e)}")
                continue

    def save_custom_colors(self, color_dialog):
        try:
            os.makedirs(os.path.dirname(DatasetSelector.custom_colors_file), exist_ok=True)
            with open(DatasetSelector.custom_colors_file, 'w', encoding='utf-8') as f:
                for i in range(16):
                    color = color_dialog.customColor(i)
                    f.write(f"{color.red()},{color.green()},{color.blue()}\n")
        except (IOError, PermissionError) as e:
            QMessageBox.critical(self, "File Error", f"Failed to save custom colors: {str(e)}")

    def populate_existing_data(self):
        logging.debug(f"Populating DataSelectionUI with input_data: {self.input_data}")

        vis_settings = self.input_data.get('visualization_settings', {})
        if vis_settings:
            self.prot_label_font.setText(vis_settings.get('prot_label_font', 'Arial'))
            self.prot_label_size.setValue(vis_settings.get('prot_label_size', 6))
            self.ptm_label_font.setText(vis_settings.get('ptm_label_font', 'Arial'))
            self.ptm_label_size.setValue(vis_settings.get('ptm_label_size', 4))
            ptm_label_color = vis_settings.get('ptm_label_color', (0, 0, 0))
            self.ptm_label_color = QColor(*ptm_label_color)
            self.ptm_label_color_square.setStyleSheet(
                f"background-color: {self.ptm_label_color.name()}; border: 1px solid gray;"
            )
            self.ptm_circle_radius.setValue(vis_settings.get('ptm_circle_radius', 5))
            self.visualization_settings.update({
                'prot_label_font': self.prot_label_font.text(),
                'prot_label_size': self.prot_label_size.value(),
                'ptm_label_font': self.ptm_label_font.text(),
                'ptm_label_size': self.ptm_label_size.value(),
                'ptm_label_color': (
                    self.ptm_label_color.red(),
                    self.ptm_label_color.green(),
                    self.ptm_label_color.blue()
                ),
                'ptm_circle_radius': self.ptm_circle_radius.value(),
                'protein_tooltip_columns': vis_settings.get('protein_tooltip_columns', [])
            })

        protein_data = self.input_data.get('protein')
        if protein_data and protein_data.get('file_path'):
            selector = self.protein_selector
            selector.file_path = protein_data['file_path']
            selector.input_text.setText(protein_data['file_path'])
            selector.update_columns()
            selector.main_columns = protein_data.get('main_columns', [])
            selector.tooltip_columns = vis_settings.get('protein_tooltip_columns', [])
            self.protein_main_columns = selector.main_columns
            display_text = ", ".join(str(col) for col in selector.main_columns)
            selector.main_columns_label.setText(f"Main Columns: {display_text}")
            tooltip_display_text = ", ".join(str(col) for col in
                                             selector.tooltip_columns) if selector.tooltip_columns else "No tooltip columns selected"
            selector.tooltip_columns_label.setText(f"Tooltip Columns: {tooltip_display_text}")
            selector.uniprot_combo.setCurrentText(protein_data.get('uniprot_column', ''))
            selector.kegg_combo.setCurrentText(protein_data.get('kegg_column', ''))
            selector.gene_combo.setCurrentText(protein_data.get('gene_column', ''))
            selector.update_ptm_buttons()

        ptm_selectors = [s for s in self.selector_widget.findChildren(DatasetSelector) if s.dataset_type != "Protein"]
        for i, ptm_data in enumerate(self.input_data.get('ptm', [])):
            if i < len(ptm_selectors):
                selector = ptm_selectors[i]
            else:
                selector = DatasetSelector(ptm_data.get('type', 'Phosphorylation'), removable=True,
                                           data_selection_ui=self)
                self.selector_layout.insertWidget(self.selector_layout.indexOf(self.add_ptm_selector), selector)

            selector.file_path = ptm_data.get('file_path', '')
            selector.input_text.setText(ptm_data['file_path'])
            selector.update_columns()

            main_columns = ptm_data.get('main_columns', [])
            if main_columns and isinstance(main_columns[0], str):
                protein_cols = self.protein_main_columns or []
                selector.main_columns = []
                for j, ptm_col in enumerate(main_columns):
                    prot_col = protein_cols[j] if j < len(protein_cols) else protein_cols[
                        -1] if protein_cols else ptm_col
                    selector.main_columns.append((prot_col, ptm_col))
            else:
                selector.main_columns = [(str(p), str(t)) for p, t in main_columns]

            selector.tooltip_columns = ptm_data.get('tooltip_columns', [])
            display_text = ", ".join(f"{p} → {t}" for p, t in selector.main_columns)
            selector.main_columns_label.setText(f"Main Columns: {display_text}")
            tooltip_display_text = ", ".join(str(col) for col in
                                             selector.tooltip_columns) if selector.tooltip_columns else "No tooltip columns selected"
            selector.tooltip_columns_label.setText(f"Tooltip Columns: {tooltip_display_text}")
            selector.uniprot_combo.setCurrentText(ptm_data.get('uniprot_column', ''))
            selector.site_combo.setCurrentText(ptm_data.get('site_column', ''))
            selector.shape_combo.setCurrentText(ptm_data.get('shape', 'Circle'))
            selector.modulation_combo.setCurrentText(ptm_data.get('modulation_column', ''))
            selector.ptm_symbol_list = ptm_data.get('ptm_symbol_list', [])
            logging.debug(f"Set PTM symbol list for {selector.dataset_type}: {selector.ptm_symbol_list}")

        self.update_data()

    def add_ptm_selector(self, ptm_type):
        logging.debug(f"Entering DataSelectionUI.add_ptm_selector for type: {ptm_type}")
        try:
            if not hasattr(self, 'protein_main_columns') or self.protein_main_columns is None:
                logging.error("protein_main_columns is None or not set")
                QMessageBox.critical(self, "Error", "Protein main columns not set. Please select Protein data first.")
                return
            new_ptm_selector = DatasetSelector(ptm_type, removable=True, data_selection_ui=self)
            logging.debug(f"New DatasetSelector created for {ptm_type}")
            new_ptm_selector.main_columns_btn.setEnabled(bool(self.protein_main_columns))
            logging.debug(f"main_columns_btn enabled: {bool(self.protein_main_columns)}")
            if self.selector_layout is None:
                logging.error("selector_layout is None")
                return
            self.selector_layout.insertWidget(self.selector_layout.indexOf(self.add_ptm_selector), new_ptm_selector)
            logging.debug("New selector added to layout")
        except Exception as e:
            logging.error(f"Error in add_ptm_selector: {str(e)}")
        logging.debug("Exiting DataSelectionUI.add_ptm_selector")

    def collect_data(self):
        data = {
            "protein": None,
            "ptm": [],
            "visualization_settings": {}
        }

        data["visualization_settings"] = {
            'prot_label_font': self.prot_label_font.text(),
            'prot_label_size': self.prot_label_size.value(),
            'ptm_label_font': self.ptm_label_font.text(),
            'ptm_label_size': self.ptm_label_size.value(),
            'ptm_label_color': (
                self.ptm_label_color.red(),
                self.ptm_label_color.green(),
                self.ptm_label_color.blue()
            ),
            'ptm_circle_radius': self.ptm_circle_radius.value(),
            'protein_tooltip_columns': []
        }

        valid_data = False
        for selector in self.selector_widget.findChildren(DatasetSelector):
            logging.debug(f"Processing {selector.dataset_type} selector with file_path: {selector.file_path}")

            if not selector.file_path:
                logging.warning(f"Skipping {selector.dataset_type} selector: No file selected")
                continue
            if selector.uniprot_combo.currentText() == "Select UniProt ID column":
                logging.warning(f"Skipping {selector.dataset_type} selector: No UniProt ID column selected")
                continue
            if not selector.main_columns:
                logging.warning(f"Skipping {selector.dataset_type} selector: No main columns selected")
                continue

            selector_data = {
                "file_path": selector.file_path,
                "uniprot_column": selector.uniprot_combo.currentText(),
                "main_columns": selector.main_columns,
                "tooltip_columns": selector.tooltip_columns if selector.dataset_type != "Protein" else []
            }

            if selector.dataset_type == "Protein":
                if (selector.kegg_combo.currentText() == "Select KEGG Gene ID column" or
                        selector.gene_combo.currentText() == "Select Gene column"):
                    logging.warning("Skipping Protein selector: Kegg or Gene column not selected")
                    continue
                selector_data.update({
                    "kegg_column": selector.kegg_combo.currentText(),
                    "gene_column": selector.gene_combo.currentText()
                })
                data["protein"] = selector_data
                data["visualization_settings"]["protein_tooltip_columns"] = selector.tooltip_columns
                valid_data = True
            else:
                if selector.site_combo.currentText() == "Select Site ID column":
                    logging.warning(f"Skipping {selector.dataset_type} selector: No Site ID column selected")
                    continue
                selector_data.update({
                    "type": selector.dataset_type,
                    "site_column": selector.site_combo.currentText(),
                    "shape": selector.shape_combo.currentText(),
                    "modulation_column": selector.modulation_combo.currentText() if selector.modulation_combo.currentText() != "Select Modulation column" else "",
                    "ptm_symbol_list": selector.ptm_symbol_list
                })
                data["ptm"].append(selector_data)
                print("m1:", data["ptm"])
                valid_data = True

        if not valid_data:
            logging.warning("No valid data collected from any selectors")
            return None

        logging.debug(f"Collected PTM symbol list for {selector.dataset_type}: {selector.ptm_symbol_list}")
        return data

    def update_data(self):
        logging.debug(f"update_data called with Qt parent: {self.parent()}, Stored parent: {self.combined_ui_parent}")
        data = self.collect_data()
        logging.debug(f"Collected data: {data}")
        if data:
            logging.debug(f"Updating data from DataSelectionUI: {data}")
            if self.combined_ui_parent and hasattr(self.combined_ui_parent, 'update_selected_data'):
                self.combined_ui_parent.update_selected_data(data)
            else:
                logging.error("No valid stored parent or update_selected_data method")
                QMessageBox.critical(self, "Internal Error", "Cannot update data: Parent reference lost.")
        else:
            logging.debug("No valid data to update")
            if self.combined_ui_parent and hasattr(self.combined_ui_parent, 'update_selected_data'):
                self.combined_ui_parent.update_selected_data(None)
            else:
                logging.error("No valid stored parent or update_selected_data method when sending None")
                QMessageBox.critical(self, "Internal Error", "Cannot update data: Parent reference lost.")
        logging.debug("update_data completed")

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("debug.log")
        ]
    )
    app = QApplication(sys.argv)
    window = QMainWindow()
    window.setWindowTitle("Data Selection UI - Debug")
    window.setGeometry(100, 100, 1200, 800)
    ui = DataSelectionUI()
    window.setCentralWidget(ui)
    window.show()
    sys.exit(app.exec_())