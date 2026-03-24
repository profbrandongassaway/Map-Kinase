from __future__ import annotations

import argparse
import json
import re
from datetime import UTC, datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_CST_PATHWAYS_DIR = BASE_DIR / "cache" / "CST_Pathways"
DEFAULT_CST_INDEX_FILE = BASE_DIR / "cache" / "CST_pathway_module_index.json"
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build_cst_pathway_index(Path(args.input_dir), Path(args.output))
    print(f"Pathways indexed: {result.get('pathway_count')}")
    print(f"Output: {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
