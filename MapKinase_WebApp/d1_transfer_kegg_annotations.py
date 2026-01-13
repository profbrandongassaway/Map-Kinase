import csv
import os
import re
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

KEGG_COL = "KEGG_Gene_ID"


def _normalize_species_name(name: str) -> str:
    return re.sub(r"[\s_]+", "", (name or "").lower())


def load_kegg_map(base_dir: str, species: str) -> Dict[str, str]:
    """
    Load KEGG conversion map by scanning annotation_files for files named like
    <KEGG_CODE>_KEGG_Conversion*. No cross-species fallback is applied.
    """
    target = _normalize_species_name(species)
    ann_dirs: List[Path] = [Path(base_dir) / "annotation_files"]
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        ann_dirs.extend([
            Path(sys._MEIPASS) / "MapKinase_WebApp" / "annotation_files",
            Path(sys._MEIPASS) / "annotation_files",
        ])
    candidates: List[tuple[str, Path]] = []
    for ann_dir in ann_dirs:
        if not ann_dir.exists():
            continue
        for path in ann_dir.iterdir():
            if not path.is_file():
                continue
            name_lower = path.name.lower()
            if "kegg_conversion" not in name_lower:
                continue
            prefix = name_lower.split("kegg_conversion", 1)[0]
            norm_prefix = _normalize_species_name(prefix)
            candidates.append((norm_prefix, path))

    def _pick(norm: str) -> Path | None:
        for cand_norm, cand_path in candidates:
            if cand_norm == norm:
                return cand_path
        return None

    path = _pick(target)
    if path is None:
        searched = [str(p) for p in ann_dirs]
        print(f"KEGG map: no conversion file found for species '{species}' (normalized='{target}'). Searched: {searched}")
        return {}

    mapping: Dict[str, str] = {}
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            fieldnames = [f.strip() if f else "" for f in (reader.fieldnames or [])]
            kegg_col = next((c for c in fieldnames if c.lower().startswith("kegg")), "KEGG_Gene_ID")
            uni_col = next((c for c in fieldnames if "uniprot" in c.lower()), "Uniprot_ID")
            print(f"KEGG map: loading '{path.name}' using columns uniprot='{uni_col}', kegg='{kegg_col}'")
            for row in reader:
                uni = (row.get(uni_col) or "").strip()
                kegg = (row.get(kegg_col) or "").strip()
                if not uni:
                    continue
                mapping[uni] = kegg
        print(f"KEGG map: loaded {len(mapping)} entries for species '{species}' from '{path.name}'")
    except Exception as exc:
        print(f"Warning: failed to load KEGG map {path}: {exc}")
    return mapping


def annotate_protein_with_kegg(dataset: Dict[str, Any], species: str, kegg_map: Dict[str, str]) -> Dict[str, Any]:
    """
    Given a dataset dict {'headers': [...], 'rows': [...]}, append KEGG_Gene_ID column using provided map.
    """
    headers: List[str] = list(dataset.get("headers") or [])
    rows: List[List[str]] = list(dataset.get("rows") or [])
    if KEGG_COL not in headers:
        headers.append(KEGG_COL)
    kegg_idx = headers.index(KEGG_COL)
    out_rows: List[List[str]] = []
    for row in rows:
        row_vals = list(row)
        if len(row_vals) < len(headers):
            row_vals.extend([""] * (len(headers) - len(row_vals)))
        uni = (row_vals[0] or "").strip()
        row_vals[kegg_idx] = kegg_map.get(uni, "")
        out_rows.append(row_vals)
    return {"headers": headers, "rows": out_rows}
