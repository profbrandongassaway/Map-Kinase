from m4_json import get_default_json, DEFAULT_DATA
import json, os

print('Running get_default_json() with no overrides...')
d = get_default_json()
print('Type:', type(d))
if isinstance(d, dict):
    print('keys:', list(d.keys()))
    print('protbox count:', len(d.get('protbox_data', [])))
else:
    print('Returned non-dict:', d)

print('\nNow running with explicit file overrides (Graduate_Documents paths)...')
# Adjust these paths if your file names differ
prot_path = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_ProtMaphsaanno.txt"
ptm_path = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_PhosMaphsanno.txt"
if not os.path.exists(prot_path):
    print('Protein file not found at prot_path:', prot_path)
if not os.path.exists(ptm_path):
    print('PTM file not found at ptm_path:', ptm_path)

override = dict(DEFAULT_DATA)
override['protein'] = dict(override.get('protein', {}))
override['protein']['file_path'] = prot_path
override['ptm'] = list(override.get('ptm', []))
if override['ptm']:
    override['ptm'][0] = dict(override['ptm'][0])
    override['ptm'][0]['file_path'] = ptm_path

d2 = get_default_json(data_override=override)
print('Type:', type(d2))
if isinstance(d2, dict):
    print('keys:', list(d2.keys()))
    print('protbox count:', len(d2.get('protbox_data', [])))
    if d2.get('protbox_data'):
        print('Sample protbox[0]:', json.dumps(d2['protbox_data'][0], indent=2))
else:
    print('Returned non-dict:', d2)
