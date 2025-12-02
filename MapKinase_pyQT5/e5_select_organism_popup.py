from PyQt5.QtWidgets import QDialog, QVBoxLayout, QRadioButton, QPushButton, QLabel, QLineEdit
from PyQt5.QtCore import Qt


class OrganismSelector(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_DeleteOnClose)  # Ensure dialog is deleted after closing
        self.organism_code = None
        self.organisms = {
            "hsa": "Homo sapiens (human)",
            "mmu": "Mus musculus (house mouse)",
            "rno": "Rattus norvegicus (rat)",
            "sce": "Saccharomyces cerevisiae (yeast)",
            "sma": "Schistosoma mansoni",
            "dme": "Drosophila melanogaster"
        }
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Select Organism")
        self.setFixedSize(400, 350)  # Increased height for better spacing

        # Create layout
        layout = QVBoxLayout()

        # Create label and place it at the top
        label = QLabel("Select Organism for Proteomic Dataset:")
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)
        layout.addSpacing(10)  # Add spacing after the label

        # Create radio buttons for predefined organisms
        self.radio_buttons = {}
        for code, name in self.organisms.items():
            rb = QRadioButton(f"{code} - {name}")
            rb.setObjectName(code)
            self.radio_buttons[code] = rb
            layout.addWidget(rb)

        # Add "Other" radio button and text input
        self.other_rb = QRadioButton("Other")
        self.other_input = QLineEdit()
        self.other_input.setPlaceholderText("Enter KEGG organism code")
        self.other_input.setEnabled(False)
        self.other_input.setFixedWidth(350)

        # Connect radio buttons to enable/disable text input
        self.other_rb.toggled.connect(self.toggle_other_input)
        for rb in self.radio_buttons.values():
            rb.toggled.connect(self.clear_other_input)

        # Create OK button with proper size
        self.ok_button = QPushButton("OK")
        self.ok_button.setFixedWidth(100)  # Set a reasonable width for the button
        self.ok_button.clicked.connect(self.accept_selection)

        # Add "Other" radio button, text input, and OK button to layout
        layout.addSpacing(10)  # Add spacing before "Other"
        layout.addWidget(self.other_rb)
        layout.addWidget(self.other_input)
        layout.addSpacing(20)  # Add spacing before the OK button
        layout.addWidget(self.ok_button, alignment=Qt.AlignCenter)  # Center the OK button

        self.setLayout(layout)

    def toggle_other_input(self, checked):
        """Enable/disable the text input when 'Other' is selected"""
        self.other_input.setEnabled(checked)
        if not checked:
            self.other_input.clear()

    def clear_other_input(self):
        """Clear the text input when a predefined organism is selected"""
        if not self.other_rb.isChecked():
            self.other_input.clear()
            self.other_input.setEnabled(False)
            self.other_rb.setChecked(False)

    def accept_selection(self):
        """Store selected organism code and close dialog"""
        for code, rb in self.radio_buttons.items():
            if rb.isChecked():
                self.organism_code = code
                self.accept()
                return
        if self.other_rb.isChecked():
            text = self.other_input.text().strip()
            if text:
                self.organism_code = text
                self.accept()
            else:
                self.other_input.setPlaceholderText("Please enter a valid KEGG code")
        else:
            self.ok_button.setText("Please select an option")

    def get_organism_code(self):
        """Return the selected organism code"""
        return self.organism_code


def select_organism(parent=None):
    """Function to show the dialog and return the selected organism code"""
    dialog = OrganismSelector(parent)
    if dialog.exec_():
        return dialog.get_organism_code()
    return None