import sys
import logging
import os
import json
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QLabel, QMessageBox, QMenuBar, QAction, QFileDialog
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from e1_popup_errorwindow import show_error_popup

try:
    logging.basicConfig(
        filename="mapkinase_log.txt",
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        filemode="w"
    )
except (IOError, PermissionError) as e:
    show_error_popup(f"Failed to initialize logging: {str(e)}. Using console logging.", "Logging Error")
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
logging.info("Program started")

from m5_phosmap_settings_qsplitter import VariableSetterUI
from m6_annotations_ui import KEGGConverterUI
from m1_data_input import DataSelectionUI

def resource_path(relative_path):
    try:
        if hasattr(sys, '_MEIPASS'):
            return os.path.join(sys._MEIPASS, relative_path)
        return os.path.join(os.path.abspath("."), relative_path)
    except Exception as e:
        logging.error(f"Error resolving resource path for {relative_path}: {str(e)}")
        return relative_path

class CombinedUI(QMainWindow):
    def __init__(self):
        super().__init__()
        logging.info("Initializing CombinedUI")
        print("CombinedUI initialized")
        self.setWindowTitle("MapKinase 0.1.4-beta")
        self.setGeometry(100, 100, 1200, 800)

        self.selected_data = None
        self.variable_setter_ui = None
        self.is_initialized = False  # Flag to track initialization
        self.pending_data_update = None  # Store pending data update

        icon_path = resource_path("App.ico")
        try:
            self.setWindowIcon(QIcon(icon_path))
            logging.debug(f"Set window icon to {icon_path}")
        except Exception as e:
            show_error_popup(f"Failed to load window icon: {str(e)}", "Resource Error", parent=self)
            logging.error(f"Failed to load icon {icon_path}: {str(e)}")

        # Create menu bar
        self.menu_bar = self.menuBar()
        self.menu_bar.setStyleSheet("""...""")  # Unchanged stylesheet
        self.file_menu = self.menu_bar.addMenu("&File")

        # Add Save action
        self.save_action = QAction("&Save", self)
        self.save_action.setShortcut("Ctrl+S")
        self.save_action.triggered.connect(self.save_all_data)
        self.file_menu.addAction(self.save_action)
        logging.debug("Added File menu with Save action")

        # Add Load action
        self.load_action = QAction("&Load", self)
        self.load_action.setShortcut("Ctrl+O")
        self.load_action.triggered.connect(self.load_all_data)
        self.file_menu.addAction(self.load_action)
        logging.debug("Added File menu with Load action")

        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)

        try:
            self.data_selection_ui = DataSelectionUI(parent=self)
            self.tab_widget.addTab(self.data_selection_ui, "Data Selection")
            logging.debug("Added DataSelectionUI tab")
            print("DataSelectionUI tab added")
        except Exception as e:
            show_error_popup(f"Failed to load Data Selection tab: {str(e)}", "Initialization Error", parent=self)
            logging.error(f"DataSelectionUI initialization failed: {str(e)}")
            placeholder = QLabel("Data Selection failed to load.")
            placeholder.setAlignment(Qt.AlignCenter)
            self.tab_widget.addTab(placeholder, "Data Selection")

        try:
            self.variable_setter_ui = VariableSetterUI(use_splitter=True, parent=self, combined_ui_parent=self)
            self.tab_widget.addTab(self.variable_setter_ui, "Pathway Visualizer")
            logging.debug("Added VariableSetterUI tab")
            print("VariableSetterUI tab added")
        except Exception as e:
            show_error_popup(f"Failed to load Pathway Visualizer tab: {str(e)}", "Initialization Error", parent=self)
            logging.error(f"VariableSetterUI initialization failed: {str(e)}")
            self.variable_setter_ui = None
            placeholder = QLabel("Pathway Visualizer failed to load.")
            placeholder.setAlignment(Qt.AlignCenter)
            self.tab_widget.addTab(placeholder, "Pathway Visualizer")

        try:
            self.kegg_converter_ui = KEGGConverterUI()
            self.tab_widget.addTab(self.kegg_converter_ui, "Add Annotations")
            logging.debug("Added KEGGConverterUI tab")
            print("KEGGConverterUI tab added")
        except Exception as e:
            show_error_popup(f"Failed to load Add Annotations tab: {str(e)}", "Initialization Error", parent=self)
            logging.error(f"KEGGConverterUI initialization failed: {str(e)}")
            placeholder = QLabel("Add Annotations failed to load.")
            placeholder.setAlignment(Qt.AlignCenter)
            self.tab_widget.addTab(placeholder, "Add Annotations")

        self.is_initialized = True  # Mark initialization complete
        self.process_pending_updates()  # Process any queued updates
        logging.debug("CombinedUI initialization completed")

    def update_selected_data(self, data):
        """Store data from DataSelectionUI and notify VariableSetterUI."""
        self.selected_data = data
        logging.debug(f"Updated selected_data in CombinedUI: {data}")
        print("CombinedUI.update_selected_data called with data:", data)

        if not self.is_initialized:
            logging.debug("Deferring VariableSetterUI update until initialization complete")
            self.pending_data_update = data
            return

        # Notify VariableSetterUI
        try:
            if self.variable_setter_ui is None:
                logging.warning("VariableSetterUI is not initialized; skipping update")
                print("Warning: VariableSetterUI is not initialized; skipping update")
                show_error_popup("Pathway Visualizer is not initialized yet.", "Initialization Error", parent=self)
                return
            if hasattr(self.variable_setter_ui, 'update_data'):
                print("Notifying VariableSetterUI with data:", data)
                self.variable_setter_ui.update_data(data)
                logging.debug("Successfully notified VariableSetterUI")
            else:
                logging.warning("VariableSetterUI has no update_data method")
                print("Warning: VariableSetterUI has no update_data method")
        except Exception as e:
            logging.error(f"Error notifying VariableSetterUI: {str(e)}")
            print(f"Error notifying VariableSetterUI: {str(e)}")
            show_error_popup(f"Failed to update Pathway Visualizer: {str(e)}", "Update Error", parent=self)

    def process_pending_updates(self):
        """Process any queued data updates after initialization."""
        if self.pending_data_update is not None:
            logging.debug(f"Processing pending data update: {self.pending_data_update}")
            print(f"Processing pending data update: {self.pending_data_update}")
            self.update_selected_data(self.pending_data_update)
            self.pending_data_update = None


    def save_all_data(self):
        """Save the current state of all tabs to a file."""
        logging.info("Attempting to save all data")
        print("Saving all data...")

        # Collect data from all tabs
        save_data = {}

        # DataSelectionUI
        try:
            if hasattr(self.data_selection_ui, 'get_current_data'):
                save_data['data_selection'] = self.data_selection_ui.get_current_data()
            else:
                save_data['data_selection'] = self.selected_data if self.selected_data else {}
            logging.debug("Collected data from DataSelectionUI")
        except Exception as e:
            logging.error(f"Failed to collect data from DataSelectionUI: {str(e)}")
            show_error_popup(f"Failed to collect data from Data Selection: {str(e)}", "Save Error", parent=self)
            save_data['data_selection'] = {"error": str(e)}

        # VariableSetterUI
        try:
            if hasattr(self.variable_setter_ui, 'get_current_data'):
                save_data['variable_setter'] = self.variable_setter_ui.get_current_data()
            else:
                save_data['variable_setter'] = {}
            logging.debug("Collected data from VariableSetterUI")
        except Exception as e:
            logging.error(f"Failed to collect data from VariableSetterUI: {str(e)}")
            show_error_popup(f"Failed to collect data from Pathway Visualizer: {str(e)}", "Save Error", parent=self)
            save_data['variable_setter'] = {"error": str(e)}

        # KEGGConverterUI
        try:
            if hasattr(self.kegg_converter_ui, 'get_current_data'):
                save_data['kegg_converter'] = self.kegg_converter_ui.get_current_data()
            else:
                save_data['kegg_converter'] = {}
            logging.debug("Collected data from KEGGConverterUI")
        except Exception as e:
            logging.error(f"Failed to collect data from KEGGConverterUI: {str(e)}")
            show_error_popup(f"Failed to collect data from Add Annotations: {str(e)}", "Save Error", parent=self)
            save_data['kegg_converter'] = {"error": str(e)}

        # Prompt user for save location
        try:
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save MapKinase Data", "", "JSON Files (*.json);;All Files (*)"
            )
            if file_path:
                if not file_path.endswith('.json'):
                    file_path += '.json'
                with open(file_path, 'w') as f:
                    json.dump(save_data, f, indent=4)
                logging.info(f"Data saved to {file_path}")
                QMessageBox.information(self, "Save Successful", f"Data saved to {file_path}")
            else:
                logging.info("Save operation cancelled by user")
        except Exception as e:
            logging.error(f"Failed to save data: {str(e)}")
            show_error_popup(f"Failed to save data: {str(e)}", "Save Error", parent=self)

    def load_all_data(self):
        """Load the state of all tabs from a file."""
        logging.info("Attempting to load all data")
        print("Loading all data...")

        # Prompt user for load location
        try:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Load MapKinase Data", "", "JSON Files (*.json);;All Files (*)"
            )
            if not file_path:
                logging.info("Load operation cancelled by user")
                return

            with open(file_path, 'r') as f:
                load_data = json.load(f)
            logging.debug(f"Loaded data from {file_path}: {load_data}")
            print(f"Loaded data from {file_path}")
            print("DEBUG: Loaded JSON variable_setter.input_data.ptm:",
                  load_data.get('variable_setter', {}).get('input_data', {}).get('ptm', []))

            # Apply data to DataSelectionUI
            try:
                if hasattr(self.data_selection_ui, 'set_current_data'):
                    self.data_selection_ui.set_current_data(load_data.get('data_selection', {}))
                    logging.debug("Applied data to DataSelectionUI")
                else:
                    logging.warning("DataSelectionUI has no set_current_data method")
                    print("Warning: DataSelectionUI has no set_current_data method")
            except Exception as e:
                logging.error(f"Failed to apply data to DataSelectionUI: {str(e)}")
                show_error_popup(f"Failed to load data for Data Selection: {str(e)}", "Load Error", parent=self)

            # Apply data to VariableSetterUI
            try:
                if hasattr(self.variable_setter_ui, 'set_current_data'):
                    self.variable_setter_ui.set_current_data(load_data.get('variable_setter', {}))
                    logging.debug("Applied data to VariableSetterUI")
                else:
                    logging.warning("VariableSetterUI has no set_current_data method")
                    print("Warning: VariableSetterUI has no set_current_data method")
            except Exception as e:
                logging.error(f"Failed to apply data to VariableSetterUI: {str(e)}")
                show_error_popup(f"Failed to load data for Pathway Visualizer: {str(e)}", "Load Error", parent=self)

            # Apply data to KEGGConverterUI
            try:
                if hasattr(self.kegg_converter_ui, 'set_current_data'):
                    self.kegg_converter_ui.set_current_data(load_data.get('kegg_converter', {}))
                    logging.debug("Applied data to KEGGConverterUI")
                else:
                    logging.warning("KEGGConverterUI has no set_current_data method")
                    print("Warning: KEGGConverterUI has no set_current_data method")
            except Exception as e:
                logging.error(f"Failed to apply data to KEGGConverterUI: {str(e)}")
                show_error_popup(f"Failed to load data for Add Annotations: {str(e)}", "Load Error", parent=self)

            # Sync m1 UI state to m5
            try:
                if hasattr(self.data_selection_ui, 'get_current_data'):
                    selected_data = self.data_selection_ui.get_current_data()
                    self.update_selected_data(selected_data)
                    logging.debug("Synced DataSelectionUI state to VariableSetterUI after JSON load")
                    print("DEBUG: Synced data_selection.ptm:", selected_data.get('ptm', []))
            except Exception as e:
                logging.error(f"Failed to sync DataSelectionUI state to VariableSetterUI: {str(e)}")
                show_error_popup(f"Failed to sync Data Selection state: {str(e)}", "Load Error", parent=self)

            QMessageBox.information(self, "Load Successful", f"Data loaded from {file_path}")
        except Exception as e:
            logging.error(f"Failed to load data: {str(e)}")
            show_error_popup(f"Failed to load data: {str(e)}", "Load Error", parent=self)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = CombinedUI()
    window.show()
    logging.info("Main window displayed")
    sys.exit(app.exec_())