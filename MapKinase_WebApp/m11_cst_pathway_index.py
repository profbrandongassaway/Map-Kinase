from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import UTC, datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_CST_PATHWAYS_DIR = BASE_DIR / "cache" / "CST_Pathways"
DEFAULT_CST_INDEX_FILE = BASE_DIR / "cache" / "CST_pathway_module_index.json"
DEFAULT_CST_LOADED_MODULES_CSV = BASE_DIR / "cache" / "cst_loaded_protein_modules.csv"
DEFAULT_CST_UNIQUE_NAMES_CSV = BASE_DIR / "cache" / "cst_unique_protein_module_names.csv"
DEFAULT_INDEX_FILES_DIR = BASE_DIR / "index_files"
_PATHWAY_STOPWORDS = {"pathway", "diagram", "interactive", "overview"}


def _normalize_text(value: str) -> str:
    text = str(value or "")
    text = text.replace("\xa0", " ")
    text = text.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _normalize_key(value: str) -> str:
    text = _normalize_text(value)
    text = re.sub(r"\s+\(\d+\)$", "", text)
    text = text.lower()
    text = text.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    return re.sub(r"[^a-z0-9]+", " ", text).strip()


def _tokenize_key(value: str) -> List[str]:
    key = _normalize_key(value)
    return [tok for tok in key.split() if tok and tok not in _PATHWAY_STOPWORDS]


def build_cst_pathway_index(
    cst_pathways_dir: Path = DEFAULT_CST_PATHWAYS_DIR,
    output_file: Path = DEFAULT_CST_INDEX_FILE,
) -> Dict[str, Any]:
    pathway_rows: List[Dict[str, Any]] = []
    for path in sorted(cst_pathways_dir.glob("*.json"), key=lambda p: p.name.lower()):
        if not path.is_file() or path.name.lower() == "_manifest.json":
            continue
        obj = json.loads(path.read_text(encoding="utf-8"))
        pathway_name = _normalize_text(str(obj.get("pathway_name") or path.stem))
        modules_out: List[Dict[str, Any]] = []
        for module in list(obj.get("protein_modules") or []):
            uniprot_ids = [str(item).strip().upper() for item in list(module.get("uniprot_ids") or []) if str(item).strip()]
            if not uniprot_ids:
                continue
            gene_symbols = [str(item).strip().upper() for item in list(module.get("psp_gene_symbols") or []) if str(item).strip()]
            if not gene_symbols:
                gene_symbols = [str(item).strip().upper() for item in list(module.get("suggested_gene_symbols") or []) if str(item).strip()]
            label = _normalize_text(str(module.get("label") or ""))
            if not label:
                continue
            modules_out.append(
                {
                    "label": label,
                    "normalized_label": _normalize_key(label).upper(),
                    "uniprot_ids": uniprot_ids,
                    "gene_symbols": gene_symbols,
                    "uniprot_urls": list(module.get("uniprot_urls") or []),
                    "manual_uniprot_override": bool(module.get("manual_uniprot_override")),
                }
            )
        pathway_rows.append(
            {
                "pathway_name": pathway_name,
                "normalized_pathway_name": _normalize_key(pathway_name).upper(),
                "pathway_file": path.name,
                "source_sections": list(obj.get("source_sections") or []),
                "protein_module_count": len(list(obj.get("protein_modules") or [])),
                "mapped_module_count": len(modules_out),
                "modules": modules_out,
            }
        )

    payload = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "source_dir": str(cst_pathways_dir.resolve()),
        "pathway_count": len(pathway_rows),
        "pathways": pathway_rows,
    }
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return payload


def _display_name_from_saved_node(node: Dict[str, Any]) -> str:
    for key in ("displayLabel", "label", "originalLabel", "annotation"):
        text = _normalize_text(str(node.get(key) or ""))
        if text:
            return text
    return ""


def _mapping_for_saved_node(node: Dict[str, Any], file_path: Path) -> Dict[str, Any]:
    from MapKinase_WebApp.m7_cst_viewer import _resolve_cst_label_mapping

    saved_ids = [str(item).strip().upper() for item in list(node.get("suggested_uniprot_ids") or []) if str(item).strip()]
    saved_genes = [str(item).strip().upper() for item in list(node.get("suggested_gene_symbols") or []) if str(item).strip()]
    saved_type = str(node.get("mappingType") or node.get("mapping_type") or "").strip().lower()
    if saved_type or saved_ids or saved_genes:
        return {
            "mapping_type": saved_type or ("cst_index" if saved_ids else "unresolved"),
            "suggested_uniprot_ids": saved_ids,
            "suggested_gene_symbols": saved_genes,
        }
    display_name = _display_name_from_saved_node(node)
    return _resolve_cst_label_mapping(display_name, file_path) if display_name else {
        "mapping_type": "unresolved",
        "suggested_uniprot_ids": [],
        "suggested_gene_symbols": [],
    }


def _module_row(
    pathway_name: str,
    protein_name: str,
    mapping: Dict[str, Any],
) -> Dict[str, str]:
    mapping_type = str(mapping.get("mapping_type") or "").strip().lower()
    uniprot_ids = [str(item).strip().upper() for item in list(mapping.get("suggested_uniprot_ids") or []) if str(item).strip()]
    gene_symbols = [str(item).strip().upper() for item in list(mapping.get("suggested_gene_symbols") or []) if str(item).strip()]
    if mapping_type == "cst_index":
        psp_ids = uniprot_ids
        backup_ids: List[str] = []
    elif mapping_type == "unresolved":
        psp_ids = []
        backup_ids = []
    else:
        psp_ids = []
        backup_ids = uniprot_ids
    return {
        "pathway": pathway_name,
        "protein_module_name": protein_name,
        "psp_uniprot_ids": "; ".join(psp_ids),
        "backup_uniprot_ids": "; ".join(backup_ids),
        "gene_symbols": "; ".join(gene_symbols),
        "mapping_type": mapping_type,
    }


def iter_effective_cst_modules(base_dir: Path) -> Iterable[Dict[str, str]]:
    from MapKinase_WebApp.m7_cst_viewer import (
        get_cst_pathway_catalog,
        load_cst_overlay_state,
    )

    for row in get_cst_pathway_catalog(base_dir):
        file_path = Path(str(row.get("file_path") or ""))
        pathway_name = _normalize_text(str(row.get("name") or file_path.stem))
        sidecar = load_cst_overlay_state(file_path)
        saved_nodes = list(sidecar.get("nodes") or [])
        if saved_nodes:
            for node in saved_nodes:
                protein_name = _display_name_from_saved_node(node)
                mapping = _mapping_for_saved_node(node, file_path)
                yield _module_row(pathway_name, protein_name, mapping)
            continue
        indexed = get_cst_pathway_mapping(pathway_name)
        for module in list(indexed.get("modules") or []):
            protein_name = _normalize_text(str(module.get("label") or ""))
            mapping = {
                "mapping_type": "cst_index",
                "suggested_uniprot_ids": [str(item).strip().upper() for item in list(module.get("uniprot_ids") or []) if str(item).strip()],
                "suggested_gene_symbols": [str(item).strip().upper() for item in list(module.get("gene_symbols") or []) if str(item).strip()],
            }
            yield _module_row(pathway_name, protein_name, mapping)


def build_cst_loaded_module_reports(
    cst_ai_dir: Path,
    loaded_modules_csv: Path = DEFAULT_CST_LOADED_MODULES_CSV,
    unique_names_csv: Path = DEFAULT_CST_UNIQUE_NAMES_CSV,
) -> Dict[str, Any]:
    rows = list(iter_effective_cst_modules(cst_ai_dir))
    rows.sort(key=lambda item: (str(item.get("pathway") or "").lower(), str(item.get("protein_module_name") or "").lower()))

    loaded_modules_csv.parent.mkdir(parents=True, exist_ok=True)
    with loaded_modules_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["pathway", "protein_module_name", "psp_uniprot_ids", "backup_uniprot_ids"],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    unique_map: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        protein_name = _normalize_text(str(row.get("protein_module_name") or ""))
        key = protein_name.upper()
        entry = unique_map.setdefault(
            key,
            {
                "protein_module_name": protein_name,
                "psp_ids": set(),
                "backup_ids": set(),
            },
        )
        for item in str(row.get("psp_uniprot_ids") or "").split(";"):
            token = item.strip().upper()
            if token:
                entry["psp_ids"].add(token)
        for item in str(row.get("backup_uniprot_ids") or "").split(";"):
            token = item.strip().upper()
            if token:
                entry["backup_ids"].add(token)

    unique_rows = [
        {
            "protein_module_name": entry["protein_module_name"],
            "psp_uniprot_count": len(entry["psp_ids"]),
            "backup_uniprot_count": len(entry["backup_ids"]),
        }
        for entry in unique_map.values()
    ]
    unique_rows.sort(key=lambda item: str(item.get("protein_module_name") or "").lower())

    with unique_names_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["protein_module_name", "psp_uniprot_count", "backup_uniprot_count"],
        )
        writer.writeheader()
        writer.writerows(unique_rows)

    return {
        "pathway_rows": len(rows),
        "unique_protein_names": len(unique_rows),
        "loaded_modules_csv": str(loaded_modules_csv.resolve()),
        "unique_names_csv": str(unique_names_csv.resolve()),
    }


def build_cst_fisher_index(
    cst_ai_dir: Path,
    output_file: Optional[Path] = None,
    species_code: str = "",
) -> Dict[str, Any]:
    grouped: Dict[str, Dict[str, Any]] = {}
    for row in iter_effective_cst_modules(cst_ai_dir):
        if not isinstance(row, dict):
            continue
        pathway_name = _normalize_text(str(row.get("pathway") or ""))
        protein_name = _normalize_text(str(row.get("protein_module_name") or ""))
        if not pathway_name or not protein_name:
            continue
        pathway_id = "cst_" + re.sub(r"[^a-z0-9]+", "_", pathway_name.strip().lower()).strip("_")
        pathway_entry = grouped.setdefault(
            pathway_id,
            {
                "pathway_id": pathway_id,
                "name": pathway_name,
                "modules": [],
                "resolved_uniprots": [],
                "gene_symbols": [],
            },
        )
        module_uniprots: List[str] = []
        for field_name in ("psp_uniprot_ids", "backup_uniprot_ids"):
            raw_val = str(row.get(field_name) or "")
            for item in raw_val.split(";"):
                token = str(item or "").strip().upper()
                if token and token not in module_uniprots:
                    module_uniprots.append(token)
                if token and token not in pathway_entry["resolved_uniprots"]:
                    pathway_entry["resolved_uniprots"].append(token)
        module_gene_symbols: List[str] = []
        raw_genes = str(row.get("gene_symbols") or "")
        for item in raw_genes.split(";"):
            symbol = str(item or "").strip().upper()
            if symbol and symbol not in module_gene_symbols:
                module_gene_symbols.append(symbol)
            if symbol and symbol not in pathway_entry["gene_symbols"]:
                pathway_entry["gene_symbols"].append(symbol)
        pathway_entry["modules"].append(
            {
                "label": protein_name,
                "uniprot_ids": module_uniprots,
                "gene_symbols": module_gene_symbols,
            }
        )

    payload = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "species_code": str(species_code or "").strip().lower(),
        "source_dir": str(cst_ai_dir.resolve()),
        "pathway_count": len(grouped),
        "pathways": list(grouped.values()),
        "nodes": {},
    }
    if output_file is not None:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return payload


@lru_cache(maxsize=4)
def load_cst_pathway_index(index_file: str = str(DEFAULT_CST_INDEX_FILE)) -> Dict[str, Any]:
    path = Path(index_file)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def get_cst_pathway_mapping(
    pathway_name_or_file: str,
    index_file: str = str(DEFAULT_CST_INDEX_FILE),
) -> Dict[str, Any]:
    index_obj = load_cst_pathway_index(index_file)
    pathways = list(index_obj.get("pathways") or [])
    if not pathways:
        return {}
    needle_norm = _normalize_key(pathway_name_or_file).upper()
    needle_tokens = set(_tokenize_key(pathway_name_or_file))
    if not needle_norm and not needle_tokens:
        return {}

    for row in pathways:
        if str(row.get("normalized_pathway_name") or "") == needle_norm:
            return row

    best_row: Optional[Dict[str, Any]] = None
    best_score = 0.0
    for row in pathways:
        hay_tokens = set(_tokenize_key(str(row.get("pathway_name") or "")))
        if not hay_tokens or not needle_tokens:
            continue
        overlap = len(needle_tokens & hay_tokens)
        if overlap <= 0:
            continue
        score = overlap / max(1.0, float(len(needle_tokens | hay_tokens)))
        if overlap >= 2 and score > best_score:
            best_score = score
            best_row = row
    return best_row or {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a consolidated CST pathway protein-module index.")
    parser.add_argument("--input-dir", default=str(DEFAULT_CST_PATHWAYS_DIR), help="Directory containing per-pathway CST JSON files.")
    parser.add_argument("--output", default=str(DEFAULT_CST_INDEX_FILE), help="Output JSON index file.")
    parser.add_argument("--report-cst-ai-dir", default="", help="If set, build the CST loaded-module reports from the .ai pathways in this directory.")
    parser.add_argument("--loaded-modules-csv", default=str(DEFAULT_CST_LOADED_MODULES_CSV), help="Output CSV for the pathway-by-pathway loaded protein modules.")
    parser.add_argument("--unique-names-csv", default=str(DEFAULT_CST_UNIQUE_NAMES_CSV), help="Output CSV for the unique loaded protein module names.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if str(args.report_cst_ai_dir or "").strip():
        result = build_cst_loaded_module_reports(
            Path(args.report_cst_ai_dir),
            Path(args.loaded_modules_csv),
            Path(args.unique_names_csv),
        )
        print(f"Pathway rows: {result.get('pathway_rows')}")
        print(f"Unique protein names: {result.get('unique_protein_names')}")
        print(f"Loaded modules CSV: {result.get('loaded_modules_csv')}")
        print(f"Unique names CSV: {result.get('unique_names_csv')}")
        return
    result = build_cst_pathway_index(Path(args.input_dir), Path(args.output))
    print(f"Pathways indexed: {result.get('pathway_count')}")
    print(f"Output: {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
