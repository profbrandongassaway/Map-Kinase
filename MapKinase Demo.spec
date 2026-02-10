# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all

datas = [('MapKinase_WebApp\\annotation_files', 'MapKinase_WebApp\\annotation_files'), ('MapKinase_WebApp\\sample_files', 'MapKinase_WebApp\\sample_files'), ('MapKinase_WebApp\\cache', 'MapKinase_WebApp\\cache'), ('MapKinase_WebApp\\exceptions_file.txt', 'MapKinase_WebApp'), ('MapKinase_WebApp\\kegg_pathways.txt', 'MapKinase_WebApp'), ('MapKinase_WebApp\\species_ref_list.csv', 'MapKinase_WebApp'), ('sample_input_files', 'sample_input_files')]
binaries = []
hiddenimports = []
tmp_ret = collect_all('shiny')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]
tmp_ret = collect_all('webview')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]


a = Analysis(
    ['MapKinase_WebApp\\m5_main_ui.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='MapKinase Demo',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
