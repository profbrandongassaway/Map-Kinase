from MapKinase_WebApp.m4_json import get_default_json, DEFAULT_DATA
import json, traceback

try:
    # Use explicit overrides pointing to the actual files under Graduate_Documents
    data_override = dict(DEFAULT_DATA)
    data_override['protein'] = dict(data_override.get('protein', {}))
    data_override['protein']['file_path'] = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_ProtMaphsaanno.txt"
    ptm_list = list(data_override.get('ptm', []))
    if ptm_list:
        ptm_list[0] = dict(ptm_list[0])
        ptm_list[0]['file_path'] = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\data_analysis\BCp2_data_03272025\Phosmap\BCp2_PhosMaphsanno.txt"
        data_override['ptm'] = ptm_list

    jd = get_default_json(data_override=data_override)
    print('TYPE:', type(jd))
    if isinstance(jd, dict):
        print('keys:', list(jd.keys()))
        pbd = jd.get('protbox_data', None)
        print('protbox_data type:', type(pbd), 'count:', len(pbd) if isinstance(pbd, list) else 'N/A')
        if isinstance(pbd, list) and len(pbd) > 0:
            print('first protbox sample:')
            print(json.dumps(pbd[0], indent=2)[:2000])
    else:
        print('Returned non-dict:', repr(jd)[:2000])
except Exception as e:
    print('ERROR calling get_default_json:')
    traceback.print_exc()
