import requests
from PIL import Image
import io
import os
import pybiopax
from base_api import BasePathwayAPI
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM


class PathBankAPI(BasePathwayAPI):
    def __init__(self):
        self.pathway_commons_base_url = "http://www.pathwaycommons.org/pc2"
        self.pathbank_base_url = "https://pathbank.org"

    def download_pathway_data(self, pathway_id):
        """Download BioPAX file for the given PathBank pathway ID from Pathway Commons."""
        try:
            # Construct Pathway Commons BioPAX URL for PathBank pathway
            uri = f"http://bioregistry.io/pathbank:{pathway_id}"
            url = f"{self.pathway_commons_base_url}/get?uri={uri}&format=BIO_PAX"
            response = requests.get(url)
            response.raise_for_status()
            file_path = f"{pathway_id}.biopax"
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            with open(file_path, "wb") as f:
                f.write(response.content)
            print(f"Successfully downloaded PathBank BioPAX file from Pathway Commons: {file_path}")
            return file_path
        except requests.exceptions.RequestException as e:
            print(f"Error downloading PathBank data for {pathway_id} from Pathway Commons: {e}")
            return None
        except (IOError, PermissionError) as e:
            print(f"Error saving BioPAX file for {pathway_id}: {e}")
            return None

    def download_pathway_image(self, pathway_id):
        """Attempt to download PNG or SVG image from PathBank, as Pathway Commons does not provide images."""
        try:
            # Try PathBank PNG first (Pathway Commons does not provide images)
            url = f"{self.pathbank_base_url}/downloads/pathways/{pathway_id}.png"
            response = requests.get(url)
            response.raise_for_status()
            print(f"Successfully downloaded PNG image: {url}")
            return Image.open(io.BytesIO(response.content))
        except requests.exceptions.RequestException as e:
            print(f"Error downloading PNG image for {pathway_id}: {e}")
            # Fallback to SVG
            try:
                url = f"{self.pathbank_base_url}/downloads/pathways/{pathway_id}.svg"
                response = requests.get(url)
                response.raise_for_status()
                print(f"Successfully downloaded SVG image: {url}")
                # Convert SVG to PIL Image (requires svglib)
                svg_io = io.BytesIO(response.content)
                drawing = svg2rlg(svg_io)
                img_io = io.BytesIO()
                renderPM.drawToFile(drawing, img_io, fmt="PNG")
                img_io.seek(0)
                return Image.open(img_io)
            except Exception as svg_e:
                print(f"Error downloading or converting SVG image for {pathway_id}: {svg_e}")
                return None
        except Exception as e:
            print(f"Unexpected error downloading image for {pathway_id}: {e}")
            return None

    def parse_pathway(self, file_path):
        """Parse BioPAX file to extract entries, groups, and arrows."""
        if not file_path or not os.path.exists(file_path):
            print(f"No valid file to parse for {file_path}")
            return [], [], []

        try:
            # Parse BioPAX file using pybiopax
            biopax_model = pybiopax.biopax.model_from_owl_file(file_path)
            entries = []
            groups = []
            arrows = []

            # Extract PhysicalEntities (e.g., proteins, small molecules)
            for obj in biopax_model.objects.values():
                if isinstance(obj, pybiopax.biopax.PhysicalEntity):
                    entity_type = "prot_box" if isinstance(obj, pybiopax.biopax.Protein) else "compound"
                    entry = {
                        "id": obj.uid,
                        "name": obj.display_name or obj.uid,
                        "type": entity_type,
                        "x": 0.0,  # No coordinates in BioPAX; assign defaults
                        "y": 0.0,
                        "width": 50.0,  # Default size; adjust in PathwayViewer
                        "height": 20.0,
                        "first_name": obj.display_name.split(",")[0].strip() if obj.display_name else obj.uid,
                        "fgcolor": "#000000",
                        "bgcolor": "#FFFFFF",
                        "xref": {
                            "Database": getattr(obj, "xref_database", ""),
                            "ID": getattr(obj, "xref_id", "")
                        }
                    }
                    entries.append(entry)
                    print(f"Parsed entry: {entry['id']} ({entry['type']})")

                # Extract Complexes as groups
                elif isinstance(obj, pybiopax.biopax.Complex):
                    group = {
                        "id": obj.uid,
                        "type": "group",
                        "name": obj.display_name or obj.uid,
                        "x": 0.0,
                        "y": 0.0,
                        "width": 100.0,  # Larger default size for complexes
                        "height": 50.0,
                        "fgcolor": "#000000",
                        "bgcolor": "#FFFFFF",
                        "components": [comp.uid for comp in obj.components]
                    }
                    groups.append(group)
                    print(f"Parsed group: {group['id']} with {len(group['components'])} components")

                # Extract BiochemicalReactions as arrows
                elif isinstance(obj, pybiopax.biopax.BiochemicalReaction):
                    source = obj.left[0].uid if obj.left else None
                    target = obj.right[0].uid if obj.right else None
                    if source and target:
                        arrow = {
                            "entry1": source,
                            "entry2": target,
                            "type": "reaction"
                        }
                        arrows.append(arrow)
                        print(f"Parsed arrow: {source} -> {target}")

            print(f"Parsed {len(entries)} entries, {len(groups)} groups, and {len(arrows)} arrows")
            return entries, groups, arrows
        except Exception as e:
            print(f"Error parsing BioPAX file {file_path}: {e}")
            return [], [], []