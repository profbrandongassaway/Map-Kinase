import requests
import xml.etree.ElementTree as ET
from PIL import Image
import io
from pathlib import Path
from MapKinase_WebApp.a1_base_api import BasePathwayAPI


def _derive_species_folder(pathway_id: str) -> str:
    # Pathway IDs are typically like hsa04010, mmu04010, etc.
    for idx, ch in enumerate(pathway_id):
        if ch.isdigit():
            return pathway_id[:idx].lower() or "unknown"
    return pathway_id.lower() or "unknown"

class KeggAPI(BasePathwayAPI):
    def download_pathway_data(self, pathway_id, species_hint=None):
        species_folder = _derive_species_folder(pathway_id)
        base_dir = Path(__file__).resolve().parent.parent / "stored_pathways" / "kegg" / species_folder
        base_dir.mkdir(parents=True, exist_ok=True)
        file_path = base_dir / f"{pathway_id}.xml"
        if file_path.exists():
            return str(file_path)
        url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
        response = requests.get(url)
        response.raise_for_status()
        with file_path.open("w", encoding="utf-8") as f:
            f.write(response.text)
        return str(file_path)

    def download_pathway_image(self, pathway_id):
        url = f"https://rest.kegg.jp/get/{pathway_id}/image"
        response = requests.get(url)
        response.raise_for_status()
        species_folder = _derive_species_folder(pathway_id)
        base_dir = Path(__file__).resolve().parent.parent / "stored_pathways" / "kegg" / species_folder
        base_dir.mkdir(parents=True, exist_ok=True)
        file_path = base_dir / f"{pathway_id}.png"
        with file_path.open("wb") as fh:
            fh.write(response.content)
        return Image.open(file_path)

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
                    components = [comp.get('id') for comp in entry.findall('component') if comp.get('id')]
                    if components:
                        entry_data['components'] = components
                    groups.append(entry_data)
                else:
                    entries.append(entry_data)

            for relation in root.findall('relation'):
                line = 'arrow'  # default
                rel_type = ''
                compound_id = None
                is_binding = False
                for subtype in relation.findall('subtype'):
                    value = subtype.get('value')
                    name_raw = subtype.get('name') or ''
                    name = name_raw.lower()
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
                    if name == 'compound' and value:
                        compound_id = value
                    elif name == 'binding/association':
                        rel_type = 'binding/association'
                        is_binding = True
                    elif name == 'phosphorylation':
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
                entry1 = relation.get('entry1')
                entry2 = relation.get('entry2')
                if is_binding and entry1 and entry2:
                    arrows.append({
                        'entry1': entry1,
                        'entry2': entry2,
                        'line': 'line',
                        'type': rel_type,
                        'binding': True
                    })
                    continue
                if compound_id:
                    if entry1 and compound_id:
                        arrows.append({
                            'entry1': entry1,
                            'entry2': compound_id,
                            'line': line,
                            'type': rel_type
                        })
                    if compound_id and entry2:
                        arrows.append({
                            'entry1': compound_id,
                            'entry2': entry2,
                            'line': line,
                            'type': rel_type
                        })
                else:
                    arrows.append({
                        'entry1': entry1,
                        'entry2': entry2,
                        'line': line,
                        'type': rel_type
                    })

            return entries, groups, arrows
        except Exception as e:
            print(f"Error parsing pathway file {file_path}: {e}")
            return [], [], []
