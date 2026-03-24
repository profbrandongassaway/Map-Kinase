from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

from MapKinase_WebApp.m8_pathway_label_mapper import PathwayLabelMapper
from MapKinase_WebApp.m11_cst_pathway_index import DEFAULT_CST_INDEX_FILE


def compare_cst_mapping_sources(index_file: Path) -> Dict[str, Any]:
    index_obj = json.loads(index_file.read_text(encoding="utf-8"))
    pathways = list(index_obj.get("pathways") or [])
    mapper = PathwayLabelMapper(use_uniprot_rest=False)

    summary = Counter()
    per_pathway: List[Dict[str, Any]] = []
    examples = {
        "psp_only": [],
        "mapper_only": [],
        "neither": [],
    }

    for pathway in pathways:
        pathway_name = str(pathway.get("pathway_name") or "")
        modules = list(pathway.get("modules") or [])
        path_counter = Counter()

        for module in modules:
            label = str(module.get("label") or "").strip()
            psp_uniprots = [str(x).strip().upper() for x in list(module.get("uniprot_ids") or []) if str(x).strip()]
            mapper_result = mapper.map_pathway_label(label)
            mapper_uniprots = [
                str(x).strip().upper()
                for x in list(mapper_result.get("suggested_uniprot_ids") or [])
                if str(x).strip()
            ]

            summary["total_modules"] += 1
            path_counter["total_modules"] += 1

            psp_hit = bool(psp_uniprots)
            mapper_hit = bool(mapper_uniprots)

            if psp_hit:
                summary["psp_recognized"] += 1
                path_counter["psp_recognized"] += 1
            if mapper_hit:
                summary["mapper_recognized"] += 1
                path_counter["mapper_recognized"] += 1
            if psp_hit and mapper_hit:
                summary["both"] += 1
                path_counter["both"] += 1
            elif psp_hit:
                summary["psp_only"] += 1
                path_counter["psp_only"] += 1
                if len(examples["psp_only"]) < 20:
                    examples["psp_only"].append(
                        {
                            "pathway": pathway_name,
                            "label": label,
                            "psp_uniprot_ids": psp_uniprots,
                            "mapper_uniprot_ids": mapper_uniprots,
                        }
                    )
            elif mapper_hit:
                summary["mapper_only"] += 1
                path_counter["mapper_only"] += 1
                if len(examples["mapper_only"]) < 20:
                    examples["mapper_only"].append(
                        {
                            "pathway": pathway_name,
                            "label": label,
                            "psp_uniprot_ids": psp_uniprots,
                            "mapper_uniprot_ids": mapper_uniprots,
                            "mapper_type": mapper_result.get("mapping_type"),
                        }
                    )
            else:
                summary["neither"] += 1
                path_counter["neither"] += 1
                if len(examples["neither"]) < 20:
                    examples["neither"].append(
                        {
                            "pathway": pathway_name,
                            "label": label,
                        }
                    )

        per_pathway.append(
            {
                "pathway_name": pathway_name,
                "total_modules": path_counter["total_modules"],
                "psp_recognized": path_counter["psp_recognized"],
                "mapper_recognized": path_counter["mapper_recognized"],
                "both": path_counter["both"],
                "psp_only": path_counter["psp_only"],
                "mapper_only": path_counter["mapper_only"],
                "neither": path_counter["neither"],
            }
        )

    per_pathway.sort(key=lambda row: (-int(row["psp_only"]), row["pathway_name"]))
    return {
        "index_file": str(index_file.resolve()),
        "summary": dict(summary),
        "per_pathway": per_pathway,
        "examples": examples,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare PSP-backed CST mappings against the gene/label->UniProt mapper.")
    parser.add_argument("--index-file", default=str(DEFAULT_CST_INDEX_FILE), help="Path to CST_pathway_module_index.json")
    parser.add_argument("--output", default="", help="Optional JSON output path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = compare_cst_mapping_sources(Path(args.index_file))
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
