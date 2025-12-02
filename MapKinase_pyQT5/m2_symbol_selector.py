import sys
import os
import json
import logging
import pandas as pd
import re
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QFrame, QComboBox, QPushButton,
    QLabel, QScrollArea, QFileDialog, QDialog, QListWidget, QGridLayout, QMessageBox, QLineEdit,
    QRadioButton, QButtonGroup, QGroupBox, QSpinBox, QColorDialog, QTabWidget, QInputDialog
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor

def resource_path(relative_path):
    try:
        if hasattr(sys, '_MEIPASS'):
            return os.path.join(sys._MEIPASS, relative_path)
        return os.path.join(os.path.abspath("."), relative_path)
    except Exception as e:
        logging.error(f"Error resolving resource path for {relative_path}: {str(e)}")
        return relative_path

def get_user_config_dir():
    return os.path.join(os.path.expanduser("~"), ".mapkinase")

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

class NoScrollSpinBox(QSpinBox):
    def wheelEvent(self, event):
        event.ignore()

class NoScrollComboBox(QComboBox):
    def wheelEvent(self, event):
        event.ignore()

class SymbolListDialog(QDialog):
    MAX_SYMBOLS = 8  # Maximum number of symbol tabs
    PRESET_DIR = "presets"  # Directory for preset files

    def __init__(self, file_path, dataset_type, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Configure PTM Symbol Labels ({dataset_type})")
        self.setGeometry(200, 200, 800, 500)
        self.file_path = file_path
        self.dataset_type = dataset_type
        self.headers = []
        self.ptm_dict = {"ptm_symbol_list": []}
        self.preset_files = []
        logging.debug(f"SymbolListDialog initialized for {dataset_type} with file_path: {file_path}")

        os.makedirs(self.PRESET_DIR, exist_ok=True)
        self.load_preset_files()

        self.layout = QVBoxLayout(self)
        self.setup_ui()

    def load_preset_files(self):
        self.preset_files = [f for f in os.listdir(self.PRESET_DIR) if f.endswith('.json')]
        logging.debug(f"Found preset files: {self.preset_files}")

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

        layout_selector_frame = QFrame()
        layout_selector_layout = QHBoxLayout(layout_selector_frame)
        layout_selector_layout.addWidget(QLabel("Preset:"))
        self.layout_selector = NoScrollComboBox()
        self.layout_selector.addItems(["Custom"] + [os.path.splitext(f)[0] for f in self.preset_files])
        self.layout_selector.setCurrentText("Custom")
        layout_selector_layout.addWidget(self.layout_selector)

        load_button = QPushButton("Load")
        load_button.clicked.connect(self.load_preset)
        layout_selector_layout.addWidget(load_button)

        save_button = QPushButton("Save")
        save_button.clicked.connect(self.save_preset)
        layout_selector_layout.addWidget(save_button)

        self.add_tab_button = QPushButton("+ Add Symbol Label")
        self.add_tab_button.clicked.connect(self.add_tab)
        layout_selector_layout.addWidget(self.add_tab_button)

        layout_selector_layout.addStretch()
        self.layout.addWidget(layout_selector_frame)

        self.tab_widget = QTabWidget()
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self.remove_tab)
        self.layout.addWidget(self.tab_widget)

        confirm_button = QPushButton("Confirm")
        confirm_button.clicked.connect(self.confirm_selection)
        self.layout.addWidget(confirm_button)

        self.update_symbol_layout("Custom")

    def load_preset(self):
        preset_name = self.layout_selector.currentText()
        if preset_name == "Custom":
            self.update_symbol_layout("Custom")
            return

        preset_path = os.path.join(self.PRESET_DIR, f"{preset_name}.json")
        try:
            with open(preset_path, 'r', encoding='utf-8') as f:
                symbol_settings = json.load(f)
            logging.debug(f"Loaded preset {preset_name}: {symbol_settings}")
        except Exception as e:
            logging.error(f"Error loading preset {preset_name}: {str(e)}")
            QMessageBox.critical(self, "Preset Error", f"Failed to load preset {preset_name}: {str(e)}")
            return

        while self.tab_widget.count():
            self.tab_widget.removeTab(0)

        for setting in symbol_settings:
            reg_site_header = find_best_header_match(self.file_path, ["Regulatory site"], self)
            reg_function_header = find_best_header_match(self.file_path, ["Regulatory site function", "ON_FUNCTION"], self)
            if setting.get("header_to_search") in ["Regulatory site", "Regulatory site function", "ON_FUNCTION"]:
                setting["header_to_search"] = reg_function_header if "function" in setting.get("header_to_search", "").lower() else reg_site_header
                if not setting["header_to_search"]:
                    setting["header_to_search"] = "None"
            self.add_tab(symbol_dict=setting, removable=True)
        self.add_tab_button.setVisible(True)
        self.tab_widget.setTabsClosable(True)
        logging.debug(f"Loaded {len(symbol_settings)} tabs from preset {preset_name}")

    def save_preset(self):
        preset_name, ok = QInputDialog.getText(self, "Save Preset", "Enter preset name:")
        if not ok or not preset_name:
            return

        preset_name = re.sub(r'[^\w\-_\. ]', '', preset_name)
        preset_path = os.path.join(self.PRESET_DIR, f"{preset_name}.json")

        symbol_list = []
        for i in range(self.tab_widget.count()):
            tab = self.tab_widget.widget(i)
            symbol_dict = {
                "symbol": tab.symbol_input.text(),
                "symbol_font": tab.font_input.text(),
                "symbol_size": tab.size_input.value(),
                "symbol_color": [tab.color_button.color.red(), tab.color_button.color.green(), tab.color_button.color.blue()],
                "symbol_x_offset": tab.x_offset_input.value(),
                "symbol_y_offset": tab.y_offset_input.value(),
                "statement_type": tab.statement_combo.currentText(),
                "header_to_search": tab.header_combo.currentText() if tab.header_combo.currentText() != "None" else "",
                "search_text_1": tab.search_text_1.text(),
                "search_text_2": tab.search_text_2.text(),
                "search_text_3": tab.search_text_3.text(),
                "search_text_4": tab.search_text_4.text()
            }
            symbol_list.append(symbol_dict)

        try:
            with open(preset_path, 'w', encoding='utf-8') as f:
                json.dump(symbol_list, f, indent=4, ensure_ascii=False)
            logging.debug(f"Saved preset to {preset_path}")
            QMessageBox.information(self, "Success", f"Preset '{preset_name}' saved successfully.")
            self.load_preset_files()
            self.layout_selector.clear()
            self.layout_selector.addItems(["Custom"] + [os.path.splitext(f)[0] for f in self.preset_files])
            self.layout_selector.setCurrentText(preset_name)
        except Exception as e:
            logging.error(f"Error saving preset {preset_name}: {str(e)}")
            QMessageBox.critical(self, "Save Error", f"Failed to save preset: {str(e)}")

    def update_symbol_layout(self, layout_type):
        while self.tab_widget.count():
            self.tab_widget.removeTab(0)
        logging.debug(f"Cleared {self.tab_widget.count()} existing tabs for layout: {layout_type}")
        self.add_tab()
        self.add_tab_button.setVisible(True)
        self.tab_widget.setTabsClosable(True)
        logging.debug("Custom layout selected, added one initial tab")

    def add_tab(self, symbol_dict=None, removable=True):
        if self.tab_widget.count() >= self.MAX_SYMBOLS:
            QMessageBox.warning(self, "Limit Reached", f"Cannot add more than {self.MAX_SYMBOLS} symbols.")
            return

        tab = QWidget()
        tab_layout = QVBoxLayout(tab)
        form_frame = QFrame()
        form_layout = QGridLayout(form_frame)
        row = 0

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

        header_label = QLabel("Header to Search:")
        form_layout.addWidget(header_label, row, 0)
        header_combo = NoScrollComboBox()
        header_combo.addItems(["None"] + self.headers)
        header_combo.setCurrentText(symbol_dict.get("header_to_search", "None") if symbol_dict else "None")
        form_layout.addWidget(header_combo, row, 1)
        row += 1

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

        self.na_label = QLabel("No search statement. Symbol will not be applied.")
        form_layout.addWidget(self.na_label, row, 0, 1, 2)
        row += 1

        tab_layout.addWidget(form_frame)
        tab_layout.addStretch()

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

        def update_visibility():
            statement_type = statement_combo.currentText()
            logging.debug(f"Updating visibility for statement_type: {statement_type}")
            header_label.setVisible(statement_type != "NA")
            header_combo.setVisible(statement_type != "NA")
            for lbl, txt in self.search_widgets:
                lbl.setVisible(False)
                txt.setVisible(False)
            self.na_label.setVisible(statement_type == "NA")
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

        statement_combo.blockSignals(True)
        update_visibility()
        statement_combo.blockSignals(False)
        statement_combo.currentTextChanged.connect(update_visibility)

        tab_index = self.tab_widget.addTab(tab, f"Symbol {self.tab_widget.count() + 1}")
        self.tab_widget.setCurrentIndex(tab_index)
        tab.removable = removable
        logging.debug(f"Added tab {tab_index + 1}, removable: {removable}, settings: {symbol_dict}")
        self.add_tab_button.setEnabled(self.tab_widget.count() < self.MAX_SYMBOLS)

    def remove_tab(self, index):
        if not self.tab_widget.widget(index).removable:
            QMessageBox.warning(self, "Cannot Remove", "This symbol cannot be removed.")
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
                symbol_dict.update({
                    "symbol": tab.symbol_input.text(),
                    "symbol_font": tab.font_input.text(),
                    "symbol_size": tab.size_input.value(),
                    "symbol_color": [tab.color_button.color.red(), tab.color_button.color.green(), tab.color_button.color.blue()],
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

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    app = QApplication(sys.argv)
    test_file_path = r"C:\Users\clayt\OneDrive - Brigham Young University\pycharm\ttk_projects\Phosmap\output\testing_file_001\pY Data2.txt"  # Replace with a valid TSV file path
    dataset_type = "Test Dataset"
    dialog = SymbolListDialog(file_path=test_file_path, dataset_type=dataset_type)
    dialog.show()
    sys.exit(app.exec_())