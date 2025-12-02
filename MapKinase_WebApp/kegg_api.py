import requests
import xml.etree.ElementTree as ET
from PIL import Image
import io
from base_api import BasePathwayAPI

class KeggAPI(BasePathwayAPI):
    def download_pathway_data(self, pathway_id):
        url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
        response = requests.get(url)
        response.raise_for_status()
        file_path = f"{pathway_id}.xml"
        with open(file_path, "w") as f:
            f.write(response.text)
        return file_path

    def download_pathway_image(self, pathway_id):
        url = f"https://rest.kegg.jp/get/{pathway_id}/image"
        response = requests.get(url)
        response.raise_for_status()
        return Image.open(io.BytesIO(response.content))

    def parse_pathway(self, file_path):
        try:
            tree = ET.parse(file_path)
            root = tree.getroot()
            entries = []
            groups = []
            arrows = []

            for entry in root.findall('entry'):
                entry_data = {
                    'id': entry.get('id'),
                    'name': entry.get('name', ''),
                    'type': 'prot_box' if entry.get('type') == 'gene' else entry.get('type', ''),
                    'x': float(entry.find('graphics').get('x', '0')),
                    'y': float(entry.find('graphics').get('y', '0')),
                    'width': float(entry.find('graphics').get('width', '0')),
                    'height': float(entry.find('graphics').get('height', '0')),
                    'first_name': entry.find('graphics').get('name', '').split(',')[0].strip(),
                    'fgcolor': entry.find('graphics').get('fgcolor', '#000000'),
                    'bgcolor': entry.find('graphics').get('bgcolor', '#FFFFFF'),
                    'graphics_type': entry.find('graphics').get('type', ''),
                    'link': entry.get('link', '')
                }
                if entry.get('type') == 'group':
                    groups.append(entry_data)
                else:
                    entries.append(entry_data)

            for relation in root.findall('relation'):
                line = 'arrow'  # default
                rel_type = ''
                for subtype in relation.findall('subtype'):
                    value = subtype.get('value')
                    name = subtype.get('name')
                    if value == '-->':
                        line = 'arrow'
                    elif value == '--|':
                        line = 'inhibition'
                    elif value == '.>' or value == '...>':
                        line = 'dashed_arrow'
                    elif value == '---':
                        line = 'dashed_arrow'
                    elif value == '-o':
                        line = 'line'
                    elif value == '-/-':
                        line = 'dashed_line'
                    if name == 'phosphorylation':
                        rel_type = 'phosphorylation' if value == '+p' else rel_type
                    elif name == 'dephosphorylation':
                        rel_type = 'dephosphorylation' if value == '-p' else rel_type
                    elif name == 'glycosylation':
                        rel_type = 'glycosylation' if value == '+g' else rel_type
                    elif name == 'ubiquitination':
                        rel_type = 'ubiquitination' if value == '+u' else rel_type
                    elif name == 'methylation':
                        rel_type = 'methylation' if value == '+m' else rel_type
                    # Add more modification types as needed
                arrows.append({
                    'entry1': relation.get('entry1'),
                    'entry2': relation.get('entry2'),
                    'line': line,
                    'type': rel_type
                })

            return entries, groups, arrows
        except Exception as e:
            print(f"Error parsing pathway file {file_path}: {e}")
            return [], [], []