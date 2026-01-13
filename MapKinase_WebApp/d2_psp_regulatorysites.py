import csv
import gzip
import os
import re
import sys
from pathlib import Path
from typing import Dict, Tuple, Any, List, Optional


PSP_FIELDS = ["DOMAIN", "ON_FUNCTION", "ON_PROCESS", "ON_PROT_INTERACT", "ON_OTHER_INTERACT", "NOTES"]


def _parse_mod_rsd(mod_rsd: str) -> Optional[str]:
    if not mod_rsd:
        return None
    match = re.match(r"^\D+(\d+)", mod_rsd)
    if not match:
        return None
    return match.group(1)


def load_regulatory_sites(base_dir: str, compressed_file: str = "Regulatory_sites.gz") -> Dict[Tuple[str, str, str], Dict[str, str]]:
    """
    Load PSP regulatory site data from the compressed file in annotation_files.
    Returns mapping: (ACC_ID, site_number, organism_lower) -> row dict.
    """
    base_path = Path(base_dir)
    candidates = [base_path / "annotation_files"]
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        candidates.extend([
            Path(sys._MEIPASS) / "MapKinase_WebApp" / "annotation_files",
            Path(sys._MEIPASS) / "annotation_files",
        ])
    full_path = None
    for folder in candidates:
        candidate = folder / compressed_file
        if candidate.exists():
            full_path = candidate
            break
    data: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    if not full_path:
        print(f"PSP regulatory sites: file not found (searched: {[str(p) for p in candidates]})")
        return data
    print(f"PSP regulatory sites: loading {full_path}")
    try:
        with gzip.open(full_path, "rt", encoding="utf-8", errors="replace") as fh:
            header: List[str] = []
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if not header:
                    # Find the real header row (contains GENE and PROTEIN)
                    if "GENE" in parts and "PROTEIN" in parts:
                        header = [col.strip() for col in parts]
                    continue
                if len(parts) < len(header):
                    parts.extend([""] * (len(header) - len(parts)))
                entry = dict(zip(header, parts))
                acc_id = (entry.get("ACC_ID") or "").strip()
                mod_rsd_num = _parse_mod_rsd(entry.get("MOD_RSD", ""))
                organism = (entry.get("ORGANISM") or "").strip().lower()
                if not acc_id or not mod_rsd_num or not organism:
                    continue
                key = (acc_id, mod_rsd_num, organism)
                data[key] = entry
        print(f"PSP regulatory sites: loaded {len(data)} entries")
    except Exception as exc:
        print(f"Warning: failed to load PSP regulatory sites: {exc}")
    return data


def annotate_ptm_dataset(dataset: Dict[str, Any], species: str, regulatory_map: Dict[Tuple[str, str, str], Dict[str, str]]) -> Dict[str, Any]:
    """
    Annotate PTM dataset (dict with 'headers' and 'rows') with PSP info.
    Adds PSP: fields and 'regulatory_site' column.
    Only applies for species in {'human','mouse','rat'}.
    """
    species_key = (species or "").strip().lower()
    if species_key not in {"human", "mouse", "rat"}:
        return dataset
    headers: List[str] = list(dataset.get("headers") or [])
    rows: List[List[str]] = list(dataset.get("rows") or [])
    psp_headers = [f"PSP: {field}" for field in PSP_FIELDS]
    new_headers = headers + ["PSP: regulatory_site"] + psp_headers
    out_rows: List[List[str]] = []
    for row in rows:
        row_vals = list(row)
        # Ensure row has base columns
        if len(row_vals) < len(headers):
            row_vals.extend([""] * (len(headers) - len(row_vals)))
        uniprot = (row_vals[0] or "").strip()
        site_raw = (row_vals[1] or "").strip()
        site_num = None
        try:
            if site_raw:
                site_num = str(int(float(site_raw)))
        except ValueError:
            site_num = None
        # Defaults
        reg_mark = ""
        psp_values = {field: "" for field in PSP_FIELDS}
        if uniprot and site_num:
            hit = regulatory_map.get((uniprot, site_num, species_key))
            if hit:
                reg_mark = "+"
                for field in PSP_FIELDS:
                    psp_values[field] = hit.get(field, "")
        row_vals.append(reg_mark)
        row_vals.extend([psp_values[f] for f in PSP_FIELDS])
        out_rows.append(row_vals)
    return {"headers": new_headers, "rows": out_rows}
