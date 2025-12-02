# File: error_popup.py
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtCore import Qt


def show_error_popup(error_message: str, title: str = "Error", parent=None) -> None:
    """
    Displays a pop-up error window using PyQt5's QMessageBox.

    Args:
        error_message (str): The error message to display.
        title (str, optional): The title of the error window. Defaults to "Error".
        parent (QWidget, optional): The parent widget for the dialog. Defaults to None.
    """
    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)
    msg_box.setText(error_message)
    msg_box.setIcon(QMessageBox.Critical)
    msg_box.setStandardButtons(QMessageBox.Ok)
    msg_box.setTextInteractionFlags(Qt.TextSelectableByMouse)  # Allows copying text
    msg_box.exec_()