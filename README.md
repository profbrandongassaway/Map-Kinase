Created by Clayton Tacker

# MapKinase
MapKinase is a Python Shiny web app for visualizing proteomics, PTM modifications, and (in the future) other omics-related data on pathway diagrams. It can pull pathways from KEGG and WikiPathways, and export pathway JSON for downstream use.

## Features
- Fetch and cache pathway layouts/images from KEGG and WikiPathways, and overlay user fold-change data.
- Incorporate PTMs into pathway maps using overlaid shapes (color-coded by fold-change) plus PhosphoSitePlus symbols for PTM functionality.
- PowerPoint-like tools in the viewer: marquee grouping, auto protbox alignment, and keyboard shortcuts for faster figure creation.
- Export pathway JSON and custom layouts.

## Quick start
1. `python -m venv .venv`
2. `.\.venv\Scripts\Activate.ps1`
3. `pip install -r requirements.txt`
4. `python MapKinase_WebApp\m5_main_ui.py`
5. Open `http://127.0.0.1:8004`

### Optional desktop window
If you want the app to open in a native window (pywebview):
```
$env:M5_DESKTOP_GUI = "1"
python MapKinase_WebApp\m5_main_ui.py
```

## Input file formats
Protein data (CSV/TSV):
- Col 1: Uniprot ID (required)
- Col 2: Gene Symbol (required)
- Col 3+: Comparison columns, headers must start with `C:` (at least one)
- Optional tooltip columns, headers start with `T:`

Example:
```
Uniprot ID,Gene Symbol,C:Control_vs_Treated,T:Notes
P27361,MAPK3,1.2,example
```

PTM data (CSV/TSV):
- Col 1: Uniprot ID (required)
- Col 2: Site Position (positive integer, required)
- Col 3+: Comparison columns, headers must start with `C:` (must match protein file comparisons)
- Optional tooltip columns, headers start with `T:`

Example:
```
Uniprot ID,Site Position,C:Control_vs_Treated
P27361,202,0.9
```

Sample inputs live in `sample_input_files` and `MapKinase_WebApp\sample_files`.

## Configuration
Environment variables:
- `M5_HOST`: bind address (default `127.0.0.1`)
- `M5_PORT`: port (default `8004`)
- `M5_DESKTOP_GUI`: set to `1` to open a desktop window
- `M5_BUILD_GLOBAL_CATALOG_ON_STARTUP`: set to `1` to rebuild the protein catalog at startup
- `M5_TERMINAL_LOG_FILE`: path to write terminal logs


## Outputs and caches
- Downloaded pathways are cached under `stored_pathways\`.
- Exported JSON previews land in `MapKinase_WebApp\JSONfiles\`.
