import pywikipathways as pwp
import xml.etree.ElementTree as ET
from collections import defaultdict
import requests

def download_wikipathway(wp_id):
    """
    Download GPML file and image for a given WikiPathways ID
    """
    # Get pathway information
    pathway_info = pwp.get_pathway_info(wp_id)

    # Download GPML
    gpml_content = pwp.get_pathway(wp_id)

    # Download image (PNG)
    image_url = f"http://www.wikipathways.org/wpi/wpi.php?action=downloadFile&type=png&pwTitle=Pathway:{wp_id}"
    image_response = requests.get(image_url)

    # Save files locally
    with open(f"{wp_id}.gpml", "w") as gpml_file:
        gpml_file.write(gpml_content)

    with open(f"{wp_id}.png", "wb") as image_file:
        image_file.write(image_response.content)

    return gpml_content

def gpml_to_dict(gpml_content):
    """
    Convert GPML content to a nested dictionary
    """
    # Parse XML content
    root = ET.fromstring(gpml_content)

    # Initialize main dictionary
    pathway_dict = defaultdict(dict)

    # Extract pathway attributes
    pathway_dict["Pathway"]["Attributes"] = {
        key: value for key, value in root.attrib.items()
    }

    # Extract Comments
    comments = []
    for comment in root.findall(".//{http://pathvisio.org/GPML/2013a}Comment"):
        comments.append({
            "Source": comment.get("Source", ""),
            "Text": comment.text.strip() if comment.text else ""
        })
    pathway_dict["Pathway"]["Comments"] = comments

    # Extract DataNodes
    data_nodes = []
    for datanode in root.findall(".//{http://pathvisio.org/GPML/2013a}DataNode"):
        node_dict = {
            "TextLabel": datanode.get("TextLabel", ""),
            "GraphId": datanode.get("GraphId", ""),
            "Type": datanode.get("Type", ""),
            "GroupRef": datanode.get("GroupRef", "")
        }

        # Get Attributes
        attributes = {}
        for attr in datanode.findall(".//{http://pathvisio.org/GPML/2013a}Attribute"):
            attributes[attr.get("Key")] = attr.get("Value")
        node_dict["Attributes"] = attributes

        # Get Graphics
        graphics = datanode.find(".//{http://pathvisio.org/GPML/2013a}Graphics")
        if graphics is not None:
            node_dict["Graphics"] = {key: value for key, value in graphics.attrib.items()}

        # Get Xref
        xref = datanode.find(".//{http://pathvisio.org/GPML/2013a}Xref")
        if xref is not None:
            node_dict["Xref"] = {
                "Database": xref.get("Database", ""),
                "ID": xref.get("ID", "")
            }

        data_nodes.append(node_dict)
    pathway_dict["DataNodes"] = data_nodes

    # Extract Interactions
    interactions = []
    for interaction in root.findall(".//{http://pathvisio.org/GPML/2013a}Interaction"):
        int_dict = {
            "GraphId": interaction.get("GraphId", "")
        }

        # Get Graphics
        graphics = interaction.find(".//{http://pathvisio.org/GPML/2013a}Graphics")
        if graphics is not None:
            int_dict["Graphics"] = {
                "Attributes": {key: value for key, value in graphics.attrib.items()},
                "Points": [
                    {key: value for key, value in point.attrib.items()}
                    for point in graphics.findall(".//{http://pathvisio.org/GPML/2013a}Point")
                ]
            }

        interactions.append(int_dict)
    pathway_dict["Interactions"] = interactions

    # Extract Groups
    groups = []
    for group in root.findall(".//{http://pathvisio.org/GPML/2013a}Group"):
        groups.append({
            "GroupId": group.get("GroupId", ""),
            "Style": group.get("Style", ""),
            "GraphId": group.get("GraphId", "")
        })
    pathway_dict["Groups"] = groups

    # Convert defaultdict to regular dict
    return dict(pathway_dict)


if __name__ == "__main__":
    wp_id = "WP382"

    # Download the pathway
    gpml_content = download_wikipathway(wp_id)

    # Convert to dictionary
    pathway_dict = gpml_to_dict(gpml_content)

    # Print some example contents
    print("Pathway Attributes:", pathway_dict["Pathway"]["Attributes"])
    print("\nFirst Comment:", pathway_dict["Pathway"]["Comments"][0])
    print("\nFirst DataNode:", pathway_dict["DataNodes"][0])
    print("\nNumber of DataNodes:", len(pathway_dict["DataNodes"]))
    print("\nNumber of Interactions:", len(pathway_dict["Interactions"]))