import os
import requests
import xml.etree.ElementTree as ET
from PIL import Image
import io
from base_api import BasePathwayAPI
from pywikipathways import get_pathway, get_pathway_info

class WikiPathwaysAPI(BasePathwayAPI):
    def download_pathway_data(self, pathway_id):
        local_file = f"{pathway_id}.gpml"
        if os.path.exists(local_file):
            print(f"Using local file: {local_file}")
            return local_file
        try:
            print(f"Fetching pathway {pathway_id} using pywikipathways")
            gpml_content = get_pathway(pathway_id)
            with open(local_file, "w", encoding="utf-8") as f:
                f.write(gpml_content)
            print(f"Successfully downloaded {local_file}")
            return local_file
        except Exception as e:
            print(f"Error downloading {pathway_id} with pywikipathways: {e}")
            try:
                print(f"Fetching pathway info for debugging: {pathway_id}")
                info = get_pathway_info(pathway_id)
                print(f"Pathway info: {info}")
            except Exception as info_e:
                print(f"Failed to get pathway info: {info_e}")
            if os.path.exists(local_file):
                print(f"Falling back to local file: {local_file}")
                return local_file
            raise Exception(f"Failed to download pathway {pathway_id} from WikiPathways")

    def download_pathway_image(self, pathway_id):
        try:
            url = f"https://www.wikipathways.org/wpi/wpi.php?action=downloadFile&type=png&pwTitle=Pathway:{pathway_id}"
            response = requests.get(url)
            response.raise_for_status()
            return Image.open(io.BytesIO(response.content))
        except requests.exceptions.RequestException as e:
            print(f"Error downloading image for {pathway_id}: {e}")
            return None

    def parse_pathway(self, file_path):
        try:
            tree = ET.parse(file_path)
            root = tree.getroot()
            ns = {'gpml': 'http://pathvisio.org/GPML/2013a'}
            print("Top-level tags in GPML:", [elem.tag for elem in root])

            entries = []
            groups = []
            arrows = []

            # Parse DataNode elements
            datanodes = root.findall(f".//gpml:DataNode", namespaces=ns)
            print(f"Found {len(datanodes)} DataNode elements")
            for datanode in datanodes:
                graph_id = datanode.get("GraphId", "Unknown")
                print(f"Processing DataNode with GraphId: {graph_id}")

                graphics = datanode.find(f"gpml:Graphics", namespaces=ns)
                if graphics is None:
                    print(f"Warning: DataNode {graph_id} is missing Graphics element, skipping")
                    continue

                center_x = graphics.get("CenterX", 0)
                center_y = graphics.get("CenterY", 0)
                width = graphics.get("Width", 0)
                height = graphics.get("Height", 0)

                if not all([center_x, center_y, width, height]):
                    print(f"Warning: DataNode {graph_id} is missing required Graphics attributes, skipping")
                    continue

                node_type = datanode.get("Type", "GeneProduct").lower()
                entry_type = 'prot_box' if node_type in ["geneproduct", "protein", "molecule"] else node_type

                text_label = datanode.get("TextLabel", "").strip()
                entry = {
                    "id": graph_id,
                    "name": text_label,
                    "type": entry_type,
                    "x": float(center_x) - float(width) / 2,
                    "y": float(center_y) - float(height) / 2,
                    "width": float(width),
                    "height": float(height),
                    "first_name": text_label.split(",")[0].strip() if text_label else "",
                    "fgcolor": graphics.get("Color", "#000000"),
                    "bgcolor": graphics.get("FillColor", "#FFFFFF"),
                    "group_ref": datanode.get("GroupRef"),
                    "xref": {}
                }

                # Add Xref information
                xref = datanode.find(f"gpml:Xref", namespaces=ns)
                if xref is not None:
                    entry["xref"] = {
                        "Database": xref.get("Database", ""),
                        "ID": xref.get("ID", "")
                    }

                entries.append(entry)
                print(f"Parsed DataNode {graph_id}: {entry}")

            # Parse Group elements
            group_elements = root.findall(f".//gpml:Group", namespaces=ns)
            print(f"Found {len(group_elements)} Group elements")
            for group in group_elements:
                group_id = group.get("GraphId", "")
                if not group_id:
                    print("Warning: Group element missing GraphId, skipping")
                    continue

                graphics = group.find(f"gpml:Graphics", namespaces=ns)
                if graphics is None:
                    print(f"Warning: Group {group_id} is missing Graphics element, skipping")
                    continue

                group_data = {
                    "id": group_id,
                    "type": "group",
                    "x": float(graphics.get("CenterX", 0)) - float(graphics.get("Width", 0)) / 2,
                    "y": float(graphics.get("CenterY", 0)) - float(graphics.get("Height", 0)) / 2,
                    "width": float(graphics.get("Width", 0)),
                    "height": float(graphics.get("Height", 0)),
                    "fgcolor": graphics.get("Color", "#000000"),
                    "bgcolor": graphics.get("FillColor", "#FFFFFF"),
                    "components": []
                }
                groups.append(group_data)
                print(f"Added group {group_id}")

            # Link DataNodes to Groups via GroupRef
            for entry in entries:
                group_ref = entry.get("group_ref")
                if group_ref:
                    for group in groups:
                        if group["id"] == group_ref:
                            group["components"].append(entry["id"])
                            print(f"Linked DataNode {entry['id']} to group {group_ref}")

            # Parse Interaction elements for arrows
            interactions = root.findall(f".//gpml:Interaction", namespaces=ns)
            print(f"Found {len(interactions)} Interaction elements")
            for interaction in interactions:
                graphics = interaction.find(f"gpml:Graphics", namespaces=ns)
                if graphics is not None:
                    points = graphics.findall(f"gpml:Point", namespaces=ns)
                    if len(points) >= 2:
                        start_point = points[0]
                        end_point = points[-1]
                        arrow = {
                            "entry1": start_point.get("GraphRef"),
                            "entry2": end_point.get("GraphRef"),
                            "type": interaction.get("Type", "interaction")
                        }
                        if arrow["entry1"] and arrow["entry2"]:
                            arrows.append(arrow)
                            print(f"Added arrow: {arrow}")

            print(f"Parsed {len(entries)} entries, {len(groups)} groups, and {len(arrows)} arrows")
            return entries, groups, arrows

        except ET.ParseError as e:
            print(f"XML Parse Error for pathway file {file_path}: {e}")
            return [], [], []
        except Exception as e:
            print(f"Error parsing pathway file {file_path}: {e}")
            return [], [], []