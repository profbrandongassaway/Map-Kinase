import sys
import os
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFormLayout, QLabel, QLineEdit,
    QPushButton, QFileDialog, QScrollArea, QProgressDialog, QMessageBox,
    QComboBox, QFrame, QSplitter, QRadioButton, QButtonGroup
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
import pandas as pd
import logging

# Import functions from d2_kegg_annotations.py and d3_transfer_kegg_annotations.py
from d2_kegg_annotations import read_uniprot_ids, uniprot_to_kegg_batch, save_results
from d3_transfer_kegg_annotations import transfer_kegg_annotations

# Import PSP functions from d1_psp_regulatorysites.py
from d1_psp_regulatorysites import read_regulatory_data, annotate_dataset



def resource_path(relative_path):
    """Get the absolute path to a resource, works for dev and PyInstaller."""
    if hasattr(sys, '_MEIPASS'):
        # PyInstaller creates a temp folder and stores files in sys._MEIPASS
        base_path = sys._MEIPASS
    else:
        # In development, use the script's directory
        base_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_path, relative_path)

# Worker thread for KEGG conversion (modified for single dataset)
class ConversionWorker(QThread):
    progress = pyqtSignal(int, int)  # Signal for progress updates (current, total)
    finished = pyqtSignal(str)  # Signal for completion with output file path
    error = pyqtSignal(str)  # Signal for errors

    def __init__(self, input_file, column_name, output_file, species, kegg_files_dir, kegg_annotation_col):
        super().__init__()
        self.input_file = input_file
        self.column_name = column_name
        self.output_file = output_file
        self.species = species
        self.kegg_files_dir = kegg_files_dir
        self.kegg_annotation_col = kegg_annotation_col
        self._is_canceled = False

    def run(self):
        try:
            if self.species == "Other":
                # Step 1: Convert UniProt IDs to KEGG IDs via API
                uniprot_ids, headers = read_uniprot_ids(self.input_file, self.column_name)
                total_batches = (len(uniprot_ids) + 100 - 1) // 100  # Ceiling division for batch_size=100
                hsa_mappings, ko_mappings = uniprot_to_kegg_batch(
                    uniprot_ids,
                    progress_callback=self.update_progress if not self._is_canceled else None
                )
                if not self._is_canceled:
                    save_results(
                        self.input_file, uniprot_ids, hsa_mappings, ko_mappings, headers,
                        self.output_file, self.column_name
                    )
                    self.progress.emit(total_batches, total_batches)
                    self.finished.emit(self.output_file)
            else:
                # Step 1: Use pre-annotated KEGG file for the species
                kegg_file = os.path.join(self.kegg_files_dir, f"{self.species}_KEGG_Conversion.txt")
                if not os.path.exists(kegg_file):
                    raise FileNotFoundError(f"KEGG file for {self.species} not found: {kegg_file}")

                transfer_kegg_annotations(
                    kegg_file,  # Pre-annotated KEGG file
                    self.input_file,
                    self.output_file,
                    "Uniprot_ID",  # Assuming UniProt ID column in KEGG file
                    self.column_name,
                    self.kegg_annotation_col
                )
                self.progress.emit(1, 1)
                self.finished.emit(self.output_file)
        except Exception as e:
            self.error.emit(str(e))

    def update_progress(self, current, total):
        if not self._is_canceled:
            self.progress.emit(current, total)

    def cancel(self):
        self._is_canceled = True


# Worker thread for PSP annotation (unchanged)
class PSPAnnotationWorker(QThread):
    progress = pyqtSignal(int, int)  # Signal for progress updates (current, total)
    finished = pyqtSignal(str)  # Signal for completion with output file path
    error = pyqtSignal(str)  # Signal for errors

    def __init__(self, psp_folder, input_file, output_file, uniprot_col, site_col):
        super().__init__()
        self.psp_folder = psp_folder
        self.input_file = input_file
        self.output_file = output_file
        self.uniprot_col = uniprot_col
        self.site_col = site_col
        self._is_canceled = False

    def run(self):
        try:
            total_steps = 2  # Loading PSP data and annotating dataset
            self.progress.emit(0, total_steps)

            regulatory_data = read_regulatory_data(self.psp_folder)
            if self._is_canceled:
                return

            self.progress.emit(1, total_steps)
            annotate_dataset(
                self.input_file,
                self.output_file,
                self.uniprot_col,
                self.site_col,
                regulatory_data
            )
            if self._is_canceled:
                return

            self.progress.emit(total_steps, total_steps)
            self.finished.emit(self.output_file)
        except Exception as e:
            self.error.emit(str(e))

    def cancel(self):
        self._is_canceled = True


class KEGGConverterUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("UniProt to KEGG and PSP Annotation Converter")
        self.setGeometry(100, 100, 900, 600)

        # Directory for KEGG files (same as script)
        self.kegg_files_dir = resource_path("Scripts")

        # Metadata file path
        self.metadata_file = "metadata.txt"

        # Create splitter as central widget
        splitter = QSplitter(Qt.Horizontal)
        self.setCentralWidget(splitter)

        # Left side: KEGG Annotation Section
        left_widget = QWidget()
        left_layout = QFormLayout(left_widget)

        kegg_scroll_area = QScrollArea()
        kegg_scroll_area.setWidgetResizable(True)
        kegg_scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        kegg_scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        kegg_content_widget = QWidget()
        self.kegg_layout = QFormLayout(kegg_content_widget)
        kegg_scroll_area.setWidget(kegg_content_widget)

        # KEGG Input Section
        self.kegg_layout.addRow(QLabel("<b>KEGG Annotation Input</b>"))

        # Species Selection (Radio Buttons)
        self.species_group = QButtonGroup(self)
        species_layout = QFormLayout()
        self.species_radios = {}
        for species in ["Human", "Mouse", "Rat", "Yeast", "DROME", "Other"]:
            radio = QRadioButton(species)
            self.species_radios[species] = radio
            self.species_group.addButton(radio)
            species_layout.addRow(radio)
        self.species_radios["Human"].setChecked(True)  # Default to Human
        self.kegg_layout.addRow(QLabel("Select Species:"), species_layout)

        # Input File (Proteomic or PTM)
        self.input_file = QLineEdit("")
        self.input_file_button = QPushButton("Browse")
        self.input_file_button.clicked.connect(lambda: self.browse_file(self.input_file, self.update_columns))
        self.kegg_layout.addRow(QLabel("Proteomic or PTM File (.txt tab-delimited):"), self.input_file)
        self.kegg_layout.addRow("", self.input_file_button)

        # UniProt Column Name (QComboBox)
        self.column_name = QComboBox()
        self.column_name.addItem("Select column after loading file")
        self.kegg_layout.addRow(QLabel("UniProt Column:"), self.column_name)

        self.kegg_layout.addRow(QLabel(""))  # Spacing
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet("QFrame { border: 1px solid gray; }")
        self.kegg_layout.addRow(separator)
        self.kegg_layout.addRow(QLabel(""))  # Spacing

        # KEGG Output Section
        self.kegg_layout.addRow(QLabel("<b>KEGG Annotation Output</b>"))

        # Output File
        self.output_file = QLineEdit("")
        self.output_file_button = QPushButton("Browse")
        self.output_file_button.clicked.connect(lambda: self.browse_output_file(self.output_file))
        self.kegg_layout.addRow(QLabel("Output File:"), self.output_file)
        self.kegg_layout.addRow("", self.output_file_button)

        # KEGG Annotation Column Name
        self.kegg_annotation_col = QLineEdit("KEGG_Gene_ID")
        self.kegg_layout.addRow(QLabel("Added KEGG Annotation Column:"), self.kegg_annotation_col)

        self.kegg_layout.addRow(QLabel(""))  # Spacing
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet("QFrame { border: 1px solid gray; }")
        self.kegg_layout.addRow(separator)
        self.kegg_layout.addRow(QLabel(""))  # Spacing

        # KEGG Run Button
        self.run_button = QPushButton("Run KEGG Annotation")
        self.run_button.clicked.connect(self.run_conversion)
        self.kegg_layout.addRow(self.run_button)

        left_layout.addRow(kegg_scroll_area)

        # Middle: PSP Annotation Section
        middle_widget = QWidget()
        middle_layout = QFormLayout(middle_widget)

        psp_scroll_area = QScrollArea()
        psp_scroll_area.setWidgetResizable(True)
        psp_scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        psp_scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        psp_content_widget = QWidget()
        self.psp_layout = QFormLayout(psp_content_widget)
        psp_scroll_area.setWidget(psp_content_widget)

        # PSP Input Section
        self.psp_layout.addRow(QLabel("<b>PSP Annotation Input</b>"))

        # PSP Folder Path
        self.psp_folder = QLineEdit(self.load_psp_folder_default())
        self.psp_folder_button = QPushButton("Browse")
        self.psp_folder_button.clicked.connect(lambda: self.browse_directory(self.psp_folder))
        self.psp_layout.addRow(QLabel("PSP Folder (containing Regulatory_sites.gz):"), self.psp_folder)
        self.psp_layout.addRow("", self.psp_folder_button)

        # Phosphoproteomic Input File
        self.psp_input_file = QLineEdit("")
        self.psp_input_file_button = QPushButton("Browse")
        self.psp_input_file_button.clicked.connect(
            lambda: self.browse_file(self.psp_input_file, self.update_psp_columns))
        self.psp_layout.addRow(QLabel("Phospho Input File (.txt tab-delimited):"), self.psp_input_file)
        self.psp_layout.addRow("", self.psp_input_file_button)

        # Phosphoproteomic UniProt Column Name (QComboBox)
        self.psp_uniprot_col = QComboBox()
        self.psp_uniprot_col.addItem("Select column after loading phospho file")
        self.psp_layout.addRow(QLabel("(Phos File) UniProt Column:"), self.psp_uniprot_col)

        # Phosphoproteomic Site Column Name (QComboBox)
        self.psp_site_col = QComboBox()
        self.psp_site_col.addItem("Select column after loading phospho file")
        self.psp_layout.addRow(QLabel("(Phos File) Phosphosite Column:"), self.psp_site_col)

        self.psp_layout.addRow(QLabel(""))  # Spacing
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet("QFrame { border: 1px solid gray; }")
        self.psp_layout.addRow(separator)
        self.psp_layout.addRow(QLabel(""))  # Spacing

        # PSP Output Section
        self.psp_layout.addRow(QLabel("<b>PSP Annotation Output</b>"))

        # PSP Output File
        self.psp_output_file = QLineEdit("")
        self.psp_output_file_button = QPushButton("Browse")
        self.psp_output_file_button.clicked.connect(lambda: self.browse_output_file(self.psp_output_file))
        self.psp_layout.addRow(QLabel("PSP Annotated Output File:"), self.psp_output_file)
        self.psp_layout.addRow("", self.psp_output_file_button)

        self.psp_layout.addRow(QLabel(""))  # Spacing
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet("QFrame { border: 1px solid gray; }")
        self.psp_layout.addRow(separator)
        self.psp_layout.addRow(QLabel(""))  # Spacing

        # PSP Run Button
        self.psp_run_button = QPushButton("Run PSP Annotation")
        self.psp_run_button.clicked.connect(self.run_psp_annotation)
        self.psp_layout.addRow(self.psp_run_button)

        middle_layout.addRow(psp_scroll_area)

        # Right side: Placeholder Section
        right_widget = QWidget()
        right_widget.setStyleSheet("background-color: #f0f0f0;")

        # Add widgets to splitter
        splitter.addWidget(left_widget)
        splitter.addWidget(middle_widget)
        splitter.addWidget(right_widget)

        # Set initial sizes (1/3 each)
        total_width = 900
        splitter.setSizes([total_width // 3, total_width // 3, total_width // 3])

        # Apply scrollbar styling
        for scroll_area in [kegg_scroll_area, psp_scroll_area]:
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

    def get_current_data(self):
        """Return the current state of KEGGConverterUI for saving."""
        print("KEGGConverterUI.get_current_data called")

        try:
            # Collect KEGG annotation settings
            species = next(
                species for species, radio in self.species_radios.items() if radio.isChecked()
            )
            kegg_data = {
                'input_file': self.input_file.text(),
                'column_name': self.column_name.currentText(),
                'output_file': self.output_file.text(),
                'kegg_annotation_col': self.kegg_annotation_col.text(),
                'species': species
            }

            # Collect PSP annotation settings
            psp_data = {
                'psp_folder': self.psp_folder.text(),
                'input_file': self.psp_input_file.text(),
                'output_file': self.psp_output_file.text(),
                'uniprot_col': self.psp_uniprot_col.currentText(),
                'site_col': self.psp_site_col.currentText()
            }

            # Combine all data
            data = {
                'kegg_annotation': kegg_data,
                'psp_annotation': psp_data
            }

            print(f"Collected data: {data}")
            return data
        except Exception as e:
            print(f"Error in get_current_data: {str(e)}")
            return {'error': str(e)}

    def set_current_data(self, data):
        """Restore the state of KEGGConverterUI from saved data."""
        logging.debug(f"KEGGConverterUI.set_current_data called with data: {data}")
        print(f"KEGGConverterUI.set_current_data called with data: {data}")

        try:
            # Apply KEGG annotation settings
            kegg_data = data.get('kegg_annotation', {})
            self.input_file.setText(kegg_data.get('input_file', ''))
            self.update_columns()  # Update column combo box based on input file
            self.column_name.setCurrentText(kegg_data.get('column_name', ''))
            self.output_file.setText(kegg_data.get('output_file', ''))
            self.kegg_annotation_col.setText(kegg_data.get('kegg_annotation_col', 'KEGG_Gene_ID'))
            species = kegg_data.get('species', 'Human')
            if species in self.species_radios:
                self.species_radios[species].setChecked(True)

            # Apply PSP annotation settings
            psp_data = data.get('psp_annotation', {})
            self.psp_folder.setText(psp_data.get('psp_folder', self.load_psp_folder_default()))
            self.psp_input_file.setText(psp_data.get('input_file', ''))
            self.update_psp_columns()  # Update column combo boxes based on input file
            self.psp_output_file.setText(psp_data.get('output_file', ''))
            self.psp_uniprot_col.setCurrentText(psp_data.get('uniprot_col', ''))
            self.psp_site_col.setCurrentText(psp_data.get('site_col', ''))

        except Exception as e:
            logging.error(f"Error in set_current_data: {str(e)}")
            print(f"Error in set_current_data: {str(e)}")
            raise

    def browse_file(self, line_edit, update_function=None):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Input File", "", "Text Files (*.txt);;All Files (*)")
        if file_path:
            line_edit.setText(file_path)
            if update_function:
                update_function()

    def browse_output_file(self, line_edit):
        file_path, _ = QFileDialog.getSaveFileName(self, "Select Output File", "", "Text Files (*.txt);;All Files (*)")
        if file_path:
            line_edit.setText(file_path)

    def browse_directory(self, line_edit):
        directory = QFileDialog.getExistingDirectory(self, "Select PSP Folder")
        if directory:
            line_edit.setText(directory)

    def load_psp_folder_default(self):
        default_folder = "PSP"
        if os.path.exists(self.metadata_file):
            try:
                with open(self.metadata_file, 'r') as f:
                    for line in f:
                        if line.startswith("psp_folder="):
                            default_folder = line.strip().split("=", 1)[1]
                            break
            except Exception as e:
                print(f"Error reading metadata.txt: {e}")
        return default_folder

    def save_psp_folder(self, folder_path):
        metadata = {}
        if os.path.exists(self.metadata_file):
            try:
                with open(self.metadata_file, 'r') as f:
                    for line in f:
                        if "=" in line:
                            key, value = line.strip().split("=", 1)
                            metadata[key] = value
            except Exception as e:
                print(f"Error reading metadata.txt: {e}")

        metadata["psp_folder"] = folder_path
        try:
            with open(self.metadata_file, 'w') as f:
                for key, value in metadata.items():
                    f.write(f"{key}={value}\n")
        except Exception as e:
            print(f"Error writing to metadata.txt: {e}")

    def update_columns(self):
        input_file = self.input_file.text()
        if os.path.exists(input_file):
            try:
                data = pd.read_csv(input_file, sep="\t")
                columns = list(data.columns)
                self.column_name.clear()
                self.column_name.addItems(columns)
                if "Uniprot_ID" in columns:
                    self.column_name.setCurrentText("Uniprot_ID")
            except Exception as e:
                QMessageBox.warning(self, "Warning", f"Error loading file columns: {e}")

    def update_psp_columns(self):
        psp_input_file = self.psp_input_file.text()
        if os.path.exists(psp_input_file):
            try:
                psp_data = pd.read_csv(psp_input_file, sep="\t")
                columns = list(psp_data.columns)
                self.psp_uniprot_col.clear()
                self.psp_uniprot_col.addItems(columns)
                self.psp_site_col.clear()
                self.psp_site_col.addItems(columns)
                if "Uniprot_ID" in columns:
                    self.psp_uniprot_col.setCurrentText("Uniprot_ID")
                if "Phosphosite" in columns:
                    self.psp_site_col.setCurrentText("Phosphosite")
            except Exception as e:
                QMessageBox.warning(self, "Warning", f"Error loading PSP phospho file columns: {e}")

    def run_conversion(self):
        input_file = self.input_file.text()
        column_name = self.column_name.currentText()
        output_file = self.output_file.text()
        kegg_annotation_col = self.kegg_annotation_col.text()
        species = next(species for species, radio in self.species_radios.items() if radio.isChecked())

        if not input_file or not output_file:
            QMessageBox.warning(self, "Input Error", "Please fill in all input and output file fields.")
            return
        if column_name == "Select column after loading file":
            QMessageBox.warning(self, "Input Error", "Please select a valid UniProt column.")
            return

        self.run_button.setEnabled(False)

        self.worker = ConversionWorker(
            input_file, column_name, output_file, species, self.kegg_files_dir, kegg_annotation_col
        )

        self.progress_dialog = QProgressDialog("Processing KEGG Annotation...", "Cancel", 0, 100, self)
        self.progress_dialog.setWindowModality(Qt.WindowModal)
        self.progress_dialog.setAutoClose(True)
        self.progress_dialog.canceled.connect(self.worker.cancel)

        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.on_conversion_finished)
        self.worker.error.connect(self.on_conversion_error)

        self.worker.start()

    def run_psp_annotation(self):
        psp_folder = self.psp_folder.text()
        input_file = self.psp_input_file.text()
        output_file = self.psp_output_file.text()
        uniprot_col = self.psp_uniprot_col.currentText()
        site_col = self.psp_site_col.currentText()

        if not psp_folder or not input_file or not output_file:
            QMessageBox.warning(self, "Input Error", "Please fill in all input and output file fields.")
            return
        if uniprot_col == "Select column after loading phospho file" or \
                site_col == "Select column after loading phospho file":
            QMessageBox.warning(self, "Input Error", "Please select valid UniProt and Phosphosite columns.")
            return

        self.psp_run_button.setEnabled(False)

        self.psp_worker = PSPAnnotationWorker(
            psp_folder, input_file, output_file, uniprot_col, site_col
        )

        self.psp_progress_dialog = QProgressDialog("Processing PSP Annotation...", "Cancel", 0, 100, self)
        self.psp_progress_dialog.setWindowModality(Qt.WindowModal)
        self.psp_progress_dialog.setAutoClose(True)
        self.psp_progress_dialog.canceled.connect(self.psp_worker.cancel)

        self.psp_worker.progress.connect(self.update_psp_progress)
        self.psp_worker.finished.connect(lambda output_file: self.on_psp_annotation_finished(output_file, psp_folder))
        self.psp_worker.error.connect(self.on_psp_annotation_error)

        self.psp_worker.start()

    def update_progress(self, current, total):
        progress = int((current / total) * 100)
        self.progress_dialog.setValue(progress)

    def update_psp_progress(self, current, total):
        progress = int((current / total) * 100)
        self.psp_progress_dialog.setValue(progress)

    def on_conversion_finished(self, output_file):
        self.run_button.setEnabled(True)
        QMessageBox.information(self, "Success", f"KEGG annotated results saved to {output_file}")

    def on_psp_annotation_finished(self, output_file, psp_folder):
        self.psp_run_button.setEnabled(True)
        self.save_psp_folder(psp_folder)
        QMessageBox.information(self, "Success", f"PSP annotated results saved to {output_file}")

    def on_conversion_error(self, error_message):
        self.run_button.setEnabled(True)
        QMessageBox.critical(self, "Error", f"Error during KEGG annotation: {error_message}")

    def on_psp_annotation_error(self, error_message):
        self.psp_run_button.setEnabled(True)
        QMessageBox.critical(self, "Error", f"Error during PSP annotation: {error_message}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = KEGGConverterUI()
    window.show()
    sys.exit(app.exec_())