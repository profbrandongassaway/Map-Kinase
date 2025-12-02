import os
import json
import requests
import io
from PIL import Image
from lxml import etree
import re

class PathwayCombiner:
    def __init__(self, svg_file):
        self.svg_file = svg_file
        self.pathbank_base_url = "https://pathbank.org"

    def download_pathway_image(self, pathway_id):
        """Download SVG or PNG image from PathBank for visualization."""
        try:
            # Try SVG first for coordinate extraction
            svg_url = f"{self.pathbank_base_url}/downloads/pathways/{pathway_id}.svg"
            response = requests.get(svg_url)
            response.raise_for_status()
            svg_path = f"{pathway_id}.svg"
            os.makedirs(os.path.dirname(svg_path) or ".", exist_ok=True)
            with open(svg_path, "wb") as f:
                f.write(response.content)
            print(f"Successfully downloaded SVG image: {svg_url}")
            # Convert SVG to PNG
            from svglib.svglib import svg2rlg
            from reportlab.graphics import renderPM
            svg_io = io.BytesIO(response.content)
            drawing = svg2rlg(svg_io)
            if drawing is None:
                print(f"Failed to convert SVG to drawing: {svg_path}")
                return None, None, None
            img_io = io.BytesIO()
            renderPM.drawToFile(drawing, img_io, fmt="PNG")
            img_io.seek(0)
            png_path = f"{pathway_id}.png"
            with open(png_path, "wb") as f:
                f.write(img_io.getvalue())
            print(f"Converted SVG to PNG: {png_path}")
            return Image.open(img_io), svg_path, png_path
        except requests.exceptions.RequestException as e:
            print(f"Error downloading SVG image for {pathway_id}: {e}")
            # Fallback to PNG
            try:
                png_url = f"{self.pathbank_base_url}/downloads/pathways/{pathway_id}.png"
                response = requests.get(png_url)
                response.raise_for_status()
                png_path = f"{pathway_id}.png"
                os.makedirs(os.path.dirname(png_path) or ".", exist_ok=True)
                with open(png_path, "wb") as f:
                    f.write(response.content)
                print(f"Successfully downloaded PNG image: {png_url}")
                return Image.open(io.BytesIO(response.content)), None, png_path
            except requests.exceptions.RequestException as e:
                print(f"Error downloading PNG image for {pathway_id}: {e}")
                return None, None, None
        except Exception as e:
            print(f"Unexpected error downloading or converting image for {pathway_id}: {e}")
            return None, None, None

    def parse_svg(self):
        """Parse SVG file to extract entries with coordinates and labels."""
        if not os.path.exists(self.svg_file):
            print(f"SVG file not found: {self.svg_file}")
            return []

        try:
            tree = etree.parse(self.svg_file)
            root = tree.getroot()
            entries = []
            ns = {"svg": "http://www.w3.org/2000/svg"}
            element_count = 0

            # Find graphical elements (<text>, <rect>, <g>) with coordinates
            for elem in root.xpath("//svg:text | //svg:rect | //svg:g", namespaces=ns):
                x, y = None, None
                label = None
                # Extract coordinates from x, y attributes
                if "x" in elem.attrib and "y" in elem.attrib:
                    try:
                        x = float(elem.attrib.get("x", 0.0))
                        y = float(elem.attrib.get("y", 0.0))
                    except ValueError:
                        print(f"Invalid x, y attributes in element: {elem.attrib}")
                        continue
                # Extract coordinates from transform attribute
                elif "transform" in elem.attrib:
                    transform = elem.attrib.get("transform", "")
                    # Match translate(x, y) or translate(x y)
                    match = re.match(r"translate\(([-]?\d*\.?\d*)\s*,?\s*([-]?\d*\.?\d*)\)", transform)
                    if match:
                        try:
                            x = float(match.group(1))
                            y = float(match.group(2))
                        except ValueError:
                            print(f"Invalid transform coordinates: {transform}")
                            continue
                    else:
                        print(f"Skipping invalid transform: {transform}")
                        continue
                # Extract label
                if elem.tag.endswith("text"):
                    label = elem.text.strip() if elem.text else None
                elif "id" in elem.attrib:
                    label = elem.attrib.get("id", None)
                if x is not None and y is not None and label:
                    element_count += 1
                    entry = {
                        "id": f"node_{element_count}",
                        "name": label,
                        "type": "node",
                        "x": x,
                        "y": y,
                        "width": 50.0,
                        "height": 20.0,
                        "first_name": label.split(",")[0].strip() if label else f"node_{element_count}",
                        "fgcolor": "#000000",
                        "bgcolor": "#FFFFFF",
                        "xref": []  # No BioPAX, so empty xref
                    }
                    entries.append(entry)
                    print(f"Parsed entry: {entry['id']} ({entry['name']}): ({x}, {y})")

            print(f"Parsed {len(entries)} entries from SVG")
            return entries
        except Exception as e:
            print(f"Error parsing SVG file {self.svg_file}: {e}")
            return []

    def save_combined_data(self, entries, png_path, output_file):
        """Save combined data to a JSON file."""
        try:
            data = {
                "entries": entries,
                "groups": [],  # No BioPAX, so empty
                "arrows": [],  # No BioPAX, so empty
                "image_path": png_path if png_path else None
            }
            output_dir = os.path.dirname(output_file)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
            with open(output_file, "w") as f:
                json.dump(data, f, indent=4)
            print(f"Combined data saved to {output_file}")
        except Exception as e:
            print(f"Error saving combined data to {output_file}: {e}")

    def combine_pathway_data(self, pathway_id, output_file=None):
        """Combine SVG data, returning entries and image, and optionally save to file."""
        # Parse SVG
        entries = self.parse_svg()

        # Download image and SVG (if not already provided)
        image, svg_path, png_path = None, self.svg_file, None
        if not svg_path or not os.path.exists(svg_path):
            image, svg_path, png_path = self.download_pathway_image(pathway_id)

        # Save to output file if specified
        if output_file:
            self.save_combined_data(entries, png_path, output_file)

        return entries, image

def main():
    # Hardcoded file paths
    svg_file = "PW000146.svg"
    pathway_id = "PW000146"
    output_file = "output.json"

    # Initialize combiner
    combiner = PathwayCombiner(svg_file=svg_file)

    # Combine data and save to file
    entries, image = combiner.combine_pathway_data(pathway_id, output_file)
    print(f"Processing complete. Output saved to {output_file}")

if __name__ == "__main__":
    main()