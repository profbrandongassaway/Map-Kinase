import csv
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, Tuple, Any, List, Optional

PSP_KINASE_COLUMNS = [
    "PSP: in_vivo_kinases",
    "PSP: in_vitro_kinases",
    "PSP: uniprot_in_vivo_kinases",
    "PSP: uniprot_in_vitro_kinases",
]


def _parse_sub_mod_rsd(mod_rsd: str) -> Optional[str]:
    if not mod_rsd:
        return None
    match = re.match(r"^\D*(\d+)", mod_rsd)
    if not match:
        return None
    return match.group(1)


def load_kinase_substrate_map(base_dir: str, compressed_file: str = "Kinase_Substrate_Dataset.gz") -> Dict[Tuple[str, str, str], List[Dict[str, str]]]:
    """
    Load PSP kinase-substrate data from compressed file in annotation_files.
    Returns mapping: (SUB_ACC_ID, site_number, sub_organism_lower) -> list of entries.
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
    data: Dict[Tuple[str, str, str], List[Dict[str, str]]] = {}
    if not full_path:
        print(f"PSP kinase-substrate: file not found (searched: {[str(p) for p in candidates]})")
        return data
    print(f"PSP kinase-substrate: loading {full_path}")
    try:
        with gzip.open(full_path, "rt", encoding="utf-8", errors="replace") as fh:
            header: List[str] = []
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if not header:
                    if "SUB_ACC_ID" in parts and "SUB_MOD_RSD" in parts:
                        header = [col.strip() for col in parts]
                    continue
                if len(parts) < len(header):
                    parts.extend([""] * (len(header) - len(parts)))
                entry = dict(zip(header, parts))
                acc = (entry.get("SUB_ACC_ID") or "").strip()
                mod_num = _parse_sub_mod_rsd(entry.get("SUB_MOD_RSD", ""))
                organism = (entry.get("SUB_ORGANISM") or "").strip().lower()
                if not acc or not mod_num or not organism:
                    continue
                key = (acc, mod_num, organism)
                data.setdefault(key, []).append(entry)
        print(f"PSP kinase-substrate: loaded {len(data)} entries")
    except Exception as exc:
        print(f"Warning: failed to load PSP kinase-substrate data: {exc}")
    return data


def annotate_ptm_dataset_with_kinases(dataset: Dict[str, Any], species: str, ks_map: Dict[Tuple[str, str, str], List[Dict[str, str]]]) -> Dict[str, Any]:
    """
    Annotate PTM dataset with PSP kinase-substrate associations.
    Adds PSP: in_vivo_kinases, PSP: in_vitro_kinases, PSP: uniprot_in_vivo_kinases, PSP: uniprot_in_vitro_kinases.
    """
    species_key = (species or "").strip().lower()
    if species_key not in {"human", "mouse", "rat"}:
        return dataset
    headers: List[str] = list(dataset.get("headers") or [])
    rows: List[List[str]] = list(dataset.get("rows") or [])
    new_headers = headers + PSP_KINASE_COLUMNS
    out_rows: List[List[str]] = []
    for row in rows:
        row_vals = list(row)
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
        vivo_kin = []
        vitro_kin = []
        vivo_uni = []
        vitro_uni = []
        if uniprot and site_num:
            matches = ks_map.get((uniprot, site_num, species_key), [])
            for hit in matches:
                kinase = (hit.get("KINASE") or "").strip()
                kin_acc = (hit.get("KIN_ACC_ID") or "").strip()
                if (hit.get("IN_VIVO_RXN") or "").strip().upper() == "X":
                    if kinase:
                        vivo_kin.append(kinase)
                    if kin_acc:
                        vivo_uni.append(kin_acc)
                if (hit.get("IN_VITRO_RXN") or "").strip().upper() == "X":
                    if kinase:
                        vitro_kin.append(kinase)
                    if kin_acc:
                        vitro_uni.append(kin_acc)
        row_vals.extend([
            "; ".join(vivo_kin),
            "; ".join(vitro_kin),
            "; ".join(vivo_uni),
            "; ".join(vitro_uni),
        ])
        out_rows.append(row_vals)
    return {"headers": new_headers, "rows": out_rows}
