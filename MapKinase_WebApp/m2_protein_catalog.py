import copy
import json
import os
from datetime import datetime
from typing import Any, Dict, Optional

import pandas as pd

from MapKinase_WebApp.m4_json import DEFAULT_DATA, DEFAULT_SETTINGS, PathwayProcessor

CATALOG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
GLOBAL_PROTEIN_CATALOG_PATH = os.path.join(CATALOG_DIR, "global_protein_catalog.json")


def _detect_sep(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    return "," if ext == ".csv" else "\t"


def _load_dataframe(entry: Any) -> Optional[pd.DataFrame]:
    if entry is None:
        return None
    if isinstance(entry, pd.DataFrame):
        return entry.copy()
    if isinstance(entry, dict):
        headers = entry.get("data_headers")
        rows = entry.get("data_rows")
        if headers is not None and rows is not None:
            try:
                return pd.DataFrame(rows, columns=headers)
            except Exception:
                pass
        path = entry.get("file_path")
        if path:
            try:
                return pd.read_csv(path, sep=_detect_sep(path))
            except Exception:
                return None
    return None


def _deep_merge(base: Dict[str, Any], override: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not override:
        return base
    for key, value in override.items():
        if (
            key in base
            and isinstance(base[key], dict)
            and isinstance(value, dict)
        ):
            base[key] = _deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def build_global_protein_catalog(
    data_override: Optional[Dict[str, Any]] = None,
    settings_override: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    data_cfg = copy.deepcopy(DEFAULT_DATA)
    _deep_merge(data_cfg, data_override)
    settings = dict(DEFAULT_SETTINGS)
    if settings_override:
        settings.update(settings_override)

    protein_cfg = data_cfg["protein"]
    proteomic_data = _load_dataframe(protein_cfg)
    ptm_datasets = []
    for dataset in data_cfg["ptm"]:
        df = _load_dataframe(dataset)
        ds_copy = copy.deepcopy(dataset)
        if df is not None:
            ds_copy["dataframe"] = df
        ptm_datasets.append(ds_copy)

    processor = PathwayProcessor([], proteomic_data, ptm_datasets, settings)
    catalog = processor.build_full_protein_catalog()
    metadata = {
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "protein_count": len(catalog),
        "protein_file": protein_cfg["file_path"],
        "ptm_files": [dataset["file_path"] for dataset in data_cfg["ptm"]],
    }
    return {"metadata": metadata, "protein_catalog": catalog}


def ensure_global_protein_catalog(
    force: bool = False,
    data_override: Optional[Dict[str, Any]] = None,
    settings_override: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    os.makedirs(CATALOG_DIR, exist_ok=True)
    if not force and os.path.exists(GLOBAL_PROTEIN_CATALOG_PATH):
        try:
            with open(GLOBAL_PROTEIN_CATALOG_PATH, "r", encoding="utf-8") as fh:
                payload = json.load(fh)
            metadata = payload.get("metadata", {})
        except (OSError, json.JSONDecodeError):
            metadata = {}
        info = {"path": GLOBAL_PROTEIN_CATALOG_PATH, "metadata": metadata}
        os.environ.setdefault("GLOBAL_PROTEIN_CATALOG_PATH", GLOBAL_PROTEIN_CATALOG_PATH)
        return info

    payload = build_global_protein_catalog(data_override=data_override, settings_override=settings_override)
    with open(GLOBAL_PROTEIN_CATALOG_PATH, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    info = {"path": GLOBAL_PROTEIN_CATALOG_PATH, "metadata": payload.get("metadata", {})}
    os.environ["GLOBAL_PROTEIN_CATALOG_PATH"] = GLOBAL_PROTEIN_CATALOG_PATH
    return info


if __name__ == "__main__":
    info = ensure_global_protein_catalog(force=True)
    print(f"Global protein catalog written to {info['path']} ({info['metadata'].get('protein_count', 0)} proteins).")
