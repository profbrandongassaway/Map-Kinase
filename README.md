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

## Temporary password protection
The main web app now requires a username and password before the site loads. The current built-in login is:

- Username: `gassawaylab`
- Password: `Qs1xULa5XCrB`

Run the app normally:

```powershell
python MapKinase_WebApp\m5_main_ui.py
```

The separate protected launcher also still works:

```powershell
python MapKinase_WebApp\m5_secure_ui.py
```

If you want to override the credentials later without editing code, set environment variables first:

```powershell
$env:MAPKINASE_LOGIN_USERNAME = "gassawaylab"
$env:MAPKINASE_LOGIN_PASSWORD = "Qs1xULa5XCrB"
python MapKinase_WebApp\m5_secure_ui.py
```

Optional environment variables:
- `MAPKINASE_AUTH_SECRET`: overrides the cookie-signing secret. If omitted, a random secret is generated each time the app starts.
- `MAPKINASE_AUTH_MAX_AGE_SECONDS`: login duration in seconds before re-authentication is required. Default is `43200` (12 hours).
- `MAPKINASE_AUTH_COOKIE_SECURE`: set to `1` when serving over HTTPS so the auth cookie is marked `Secure`.
- `MAPKINASE_ENABLE_LOGIN`: set to `0` only if you intentionally want to disable the login gate.

The protected launcher blocks both normal page requests and the Shiny websocket until the login succeeds.

## Required PSP annotations
To run the PTM annotation features properly, download the PhosphoSitePlus datasets from `https://www.phosphosite.org/staticDownloads` and place the compressed files (do not unzip) in `MapKinase_WebApp\annotation_files\`:
- `Phosphorylation_site_dataset.gz`
- `Regulatory_sites.gz`
- `Kinase_Substrate_Dataset.gz`

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

## Banned pathway audit
To build a manifest of pathways that appear in KEGG or WikiPathways catalogs but fail to download usable content, run:

```powershell
python MapKinase_WebApp\build_banned_pathways.py --source all --pretty
```

This writes `MapKinase_WebApp\index_files\pathway_banned_list.json`. The app uses that manifest to hide banned KEGG and WikiPathways entries from the pathway selectors.

## WikiPathways GitHub sync
To preload WikiPathways from the upstream GitHub repositories instead of downloading pathways one at a time through the API, run:

```powershell
python MapKinase_WebApp\sync_wikipathways_github_cache.py --pretty
```

This clones or updates:
- `wikipathways/wikipathways-database`
- `wikipathways/wikipathways-assets`

and syncs the pathway files into `stored_pathways\wikipathways\<species>\...`.

Useful options:
- `--org hsa --org mmu` to limit the sync to selected organisms from `species_ref_list.csv`
- `--max-pathways 25` for a small test run
- `--skip-fetch` to reuse an existing local clone without pulling updates first

After the sync, Map-Kinase will use local cached WikiPathways GPML files automatically, and it will use cached PNG files when available.
