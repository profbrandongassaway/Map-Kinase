import sys
import re
import logging
import os
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QLineEdit, QListWidget, QListWidgetItem, \
    QPushButton, QHBoxLayout
from PyQt5.QtCore import Qt, QRect, QEvent, QPoint
from PyQt5.QtGui import QFocusEvent, QIcon, QPixmap, QColor, QPainter, QBrush, QPolygon

# Set up logging to a file for debugging
logging.basicConfig(
    filename="../scripts/search_bar_log.txt",
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


class CustomLineEdit(QLineEdit):
    def __init__(self, parent, results_list, variables, favorite_button):
        super().__init__(parent.centralWidget())  # Set parent to the central widget for widget hierarchy
        self.parent_search_bar = parent  # Store reference to RegexSearchBar instance
        self.results_list = results_list
        self.variables = variables
        self.favorite_button = favorite_button
        self.setPlaceholderText("Click or type to search pathways...")
        self.favorites = parent.favorites  # Reference to favorites set from RegexSearchBar

    def focusInEvent(self, event: QFocusEvent):
        super().focusInEvent(event)
        logging.debug("Search bar gained focus")
        self.show_results_list()
        self.update_results(self.text())

    def focusOutEvent(self, event: QFocusEvent):
        super().focusOutEvent(event)
        if not self.results_list.hasFocus():
            logging.debug("Search bar lost focus, hiding results")
            self.results_list.hide()

    def show_results_list(self):
        # Map the position relative to the top-level window
        top_level_window = self.window()  # Get the top-level window (VariableSetterUI or RegexSearchBar)
        pos = self.mapTo(top_level_window, self.rect().bottomLeft())
        # Ensure the results list has a reasonable width and height
        list_width = max(self.width(), 300)  # Increased minimum width for better visibility
        list_height = min(self.results_list.count() * 30 + 10, 200)  # Dynamic height, max 200, min 50
        if list_height < 50:
            list_height = 50  # Ensure minimum height to avoid clipping
        # Adjust position to ensure the list stays within the top-level window
        window_rect = top_level_window.rect()
        if pos.y() + list_height > window_rect.height():
            pos.setY(pos.y() - list_height - self.height())  # Move above the search bar if it would go off-screen
        if pos.x() + list_width > window_rect.width():
            pos.setX(window_rect.width() - list_width)  # Adjust x to stay within window
        self.results_list.setGeometry(QRect(pos.x(), pos.y(), list_width, list_height))
        self.results_list.show()
        self.results_list.raise_()
        self.results_list.setFocusPolicy(Qt.StrongFocus)  # Ensure it can receive focus
        logging.debug(f"Showing results list at position: {pos.x()}, {pos.y()} with size: {list_width}x{list_height}")
        logging.debug(f"Results list visible: {self.results_list.isVisible()}, items: {self.results_list.count()}")
        logging.debug(f"Top-level window size: {window_rect.width()}x{window_rect.height()}")

    def update_results(self, text):
        try:
            logging.debug(f"Updating results for text: {text}")
            self.results_list.clear()
            if not text:
                # Separate favorited and non-favorited pathways
                favorited_vars = [var for var in self.variables if var['id'] in self.favorites]
                non_favorited_vars = [var for var in self.variables if var['id'] not in self.favorites]

                # Add favorited pathways first
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

                # Add non-favorited pathways
                for var in non_favorited_vars:
                    item = QListWidgetItem(var['display'])
                    item.setData(Qt.UserRole, var['id'])
                    self.results_list.addItem(item)

                logging.debug(f"Showing all {len(self.variables)} pathways ({len(favorited_vars)} favorited)")
                self.show_results_list()
                return

            try:
                pattern = re.compile(text, re.IGNORECASE)
            except re.error as e:
                logging.warning(f"Invalid regex pattern: {text}, error: {str(e)}")
                self.results_list.addItem("Invalid regex pattern")
                self.show_results_list()
                return

            # Filter variables that match the regex in either the ID or the display text
            matches = [var for var in self.variables if pattern.search(var['id']) or pattern.search(var['display'])]
            if matches:
                # Separate favorited and non-favorited matches
                favorited_matches = [var for var in matches if var['id'] in self.favorites]
                non_favorited_matches = [var for var in matches if var['id'] not in self.favorites]

                # Add favorited matches first
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

                # Add non-favorited matches
                for var in non_favorited_matches:
                    item = QListWidgetItem(var['display'])
                    item.setData(Qt.UserRole, var['id'])
                    self.results_list.addItem(item)

                logging.debug(f"Found {len(matches)} matches for text: {text} ({len(favorited_matches)} favorited)")
            else:
                logging.debug("No matches found")
                self.results_list.addItem("No matches found")
            self.show_results_list()
        except Exception as e:
            logging.error(f"Error in update_results: {str(e)}")


class RegexSearchBar(QMainWindow):
    def __init__(self):
        super().__init__()
        try:
            logging.info("Initializing RegexSearchBar")
            self.setWindowTitle("Regex Search Bar - PyQt5")
            self.setGeometry(100, 100, 400, 400)

            # Determine the directory of search_bar.py
            self.script_dir = os.path.dirname(os.path.abspath(__file__))
            logging.debug(f"Script directory: {self.script_dir}")

            # Load favorites from file
            self.favorites = set()
            favorites_path = os.path.join(self.script_dir, "favorites.txt")
            logging.debug(f"Looking for favorites.txt at: {favorites_path}")
            try:
                if os.path.exists(favorites_path):
                    with open(favorites_path, "r", encoding="utf-8") as f:
                        for line in f:
                            if line.strip():
                                self.favorites.add(line.strip())
                    logging.info(f"Loaded {len(self.favorites)} favorites from {favorites_path}")
                else:
                    logging.warning(f"favorites.txt not found at {favorites_path}")
            except Exception as e:
                logging.error(f"Error reading favorites.txt: {str(e)}")

            # Read pathway names and IDs from kegg_pathways.txt
            self.variables = []
            pathways_path = os.path.join(self.script_dir, "kegg_pathways.txt")
            logging.debug(f"Looking for kegg_pathways.txt at: {pathways_path}")
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
                        logging.info(f"Loaded {len(self.variables)} pathway entries from {pathways_path}")
                else:
                    logging.error(f"kegg_pathways.txt not found at {pathways_path}")
                    self.variables = [{"name": f"Error: kegg_pathways.txt not found at {pathways_path}", "id": "",
                                       "display": f"Error: kegg_pathways.txt not found at {pathways_path}"}]
            except Exception as e:
                logging.error(f"Error reading kegg_pathways.txt: {str(e)}")
                self.variables = [
                    {"name": f"Error reading file: {str(e)}", "id": "", "display": f"Error reading file: {str(e)}"}]

            # Log the number of variables loaded for debugging
            logging.debug(f"Total variables loaded: {len(self.variables)}")
            if self.variables and "Error" in self.variables[0]["display"]:
                logging.debug("Variables contain error message, likely due to file not found")

            # Main widget and layout
            main_widget = QWidget()
            self.setCentralWidget(main_widget)
            main_widget.setFocusPolicy(Qt.StrongFocus)
            layout = QVBoxLayout()
            layout.setSpacing(5)  # Reduce spacing between elements
            main_widget.setLayout(layout)

            # Results list, now a child of main_widget
            self.results_list = QListWidget(main_widget)
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
            self.results_list.focusOutEvent = lambda event: self.hide_results_if_no_focus()
            self.results_list.itemClicked.connect(self.on_item_clicked)

            # Search bar and favorite button layout
            search_layout = QHBoxLayout()

            # Favorite button
            self.favorite_button = QPushButton(main_widget)
            self.favorite_button.setFixedSize(24, 24)
            self.update_favorite_button("")  # Initialize with empty pathway
            self.favorite_button.clicked.connect(self.toggle_favorite)
            search_layout.addWidget(self.favorite_button)

            # Search bar
            self.search_bar = CustomLineEdit(self, self.results_list, self.variables, self.favorite_button)
            self.search_bar.textChanged.connect(self.on_search_text_changed)
            search_layout.addWidget(self.search_bar)

            layout.addLayout(search_layout)
            logging.info("UI setup complete")
        except Exception as e:
            logging.error(f"Error in __init__: {str(e)}")
            raise

    def eventFilter(self, obj, event):
        if obj == self.centralWidget() and event.type() == QEvent.MouseButtonPress:
            logging.debug("Main widget clicked, clearing search bar focus")
            self.search_bar.clearFocus()
            return True
        return super().eventFilter(obj, event)

    def hide_results_if_no_focus(self):
        if not self.search_bar.hasFocus() and not self.results_list.hasFocus():
            logging.debug("Results list lost focus, hiding")
            self.results_list.hide()

    def on_item_clicked(self, item):
        pathway_id = item.data(Qt.UserRole)
        if pathway_id:
            logging.debug(f"Item clicked: {item.text()}, setting search bar to: {pathway_id}")
            self.search_bar.setText(pathway_id)
            self.results_list.hide()
            self.search_bar.setFocus()

    def on_search_text_changed(self, text):
        self.search_bar.update_results(text)
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
        self.favorite_button.setEnabled(bool(pathway_id))  # Disable if no pathway ID

    def toggle_favorite(self):
        pathway_id = self.search_bar.text()
        if not pathway_id:
            return
        if pathway_id in self.favorites:
            self.favorites.remove(pathway_id)
            logging.debug(f"Removed favorite: {pathway_id}")
        else:
            self.favorites.add(pathway_id)
            logging.debug(f"Added favorite: {pathway_id}")
        # Save favorites to file
        favorites_path = os.path.join(self.script_dir, "favorites.txt")
        try:
            with open(favorites_path, "w", encoding="utf-8") as f:
                for fav in self.favorites:
                    f.write(f"{fav}\n")
            logging.info(f"Favorites saved to {favorites_path}")
        except Exception as e:
            logging.error(f"Error saving favorites.txt: {str(e)}")
        self.update_favorite_button(pathway_id)
        self.search_bar.update_results(self.search_bar.text())  # Refresh results list


if __name__ == '__main__':
    try:
        logging.info("Starting application")
        app = QApplication(sys.argv)
        window = RegexSearchBar()
        window.show()
        logging.info("Window shown")
        sys.exit(app.exec_())
    except ImportError as e:
        logging.error(f"ImportError: {str(e)} - Ensure PyQt5 is installed")
        print("Error: PyQt5 not found. Install it using 'pip install PyQt5'")
    except Exception as e:
        logging.error(f"Application error: {str(e)}")
        print(f"An error occurred: {str(e)}. Check search_bar_log.txt for details.")