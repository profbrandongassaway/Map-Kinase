#!/usr/bin/env python3
"""
m6_rank_pathways.py

Fast pathway ranking for MapKinase using:
1) one or more prebuilt pathway index JSON files (KEGG and/or WikiPathways), and
2) user protein/site evidence tables.

This module is intentionally standalone so m5 can call into it later without
changing scoring logic.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd


LOGGER = logging.getLogger("m6_rank_pathways")
BASE_DIR = Path(__file__).resolve().parent
INDEX_FILES_DIR = BASE_DIR / "index_files"

DEFAULT_WEIGHTS: Dict[str, float] = {
    "sig_scale": 5.0,
    "eff_scale": 2.0,
    "w_ann": 1.0,
    "ptm_site_scale": 0.3,
    "ptm_weight": 1.0,
    "epsilon": 0.2,
    "reg_gate": 0.15,
    "two_hop_base": 0.7,
    "conn2_weight": 1.0,
    "alpha": 0.5,
    "node_mass_weight": 0.2,
    "node_mass_top_k": 10.0,
    "site_top_k": 2.0,
    "top_edges_n": 10.0,
}

UNIPROT_SIMPLE_RE = re.compile(r"^[A-Z0-9]{6,10}(?:-\d+)?$")


@dataclass
class ProteinLookup:
    exact: Dict[str, Dict[str, Any]]
    by_base: Dict[str, Dict[str, Any]]

    def get(self, uniprot: str) -> Optional[Dict[str, Any]]:
        key = normalize_uniprot(uniprot)
        if not key:
            return None
        direct = self.exact.get(key)
        if direct is not None:
            return direct
        base = key.split("-", 1)[0]
        return self.by_base.get(base)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rank pathway indices (KEGG/WikiPathways) for a user dataset.")
    parser.add_argument(
        "--kegg_index",
        default=None,
        help=(
            "KEGG index JSON path or filename. If a filename is provided, it is "
            "resolved under MapKinase_WebApp/index_files."
        ),
    )
    parser.add_argument(
        "--wikipathways_index",
        default=None,
        help=(
            "Optional WikiPathways index JSON path or filename. If a filename is "
            "provided, it is resolved under MapKinase_WebApp/index_files."
        ),
    )
    parser.add_argument("--protein_table", required=True, help="Protein-level CSV/TSV table.")
    parser.add_argument("--site_table", default=None, help="Optional site-level CSV/TSV table.")
    parser.add_argument("--out", required=True, help="Output path for ranked pathways.")
    parser.add_argument("--format", default="json", choices=["json", "csv"], help="Output format.")

    parser.add_argument("--protein_id_col", default="Uniprot_ID", help="Protein UniProt column in protein table.")
    parser.add_argument("--p_col_prot", default="p_value", help="Protein-level p-value column.")
    parser.add_argument("--fc_col_prot", default="fold_change", help="Protein-level fold-change column.")
    parser.add_argument("--p_col_phospho", default=None, help="Phospho aggregate p-value column (typically in protein table).")
    parser.add_argument("--fc_col_phospho", default=None, help="Phospho aggregate fold-change column.")
    parser.add_argument("--p_col_site", default=None, help="Site-level p-value column (optional; fallback is --p_col_phospho).")
    parser.add_argument("--fc_col_site", default=None, help="Site-level fold-change column (optional; fallback is --fc_col_phospho).")

    parser.add_argument("--site_uniprot_col", default=None, help="UniProt column in site table (default: auto/first column).")
    parser.add_argument("--site_key_col", default="site_key", help='Existing site key column (default: "site_key").')
    parser.add_argument(
        "--site_key_cols",
        default=None,
        help='Build site_key from columns. Preferred format: "uniprot_col,residue_col,position_col".',
    )
    parser.add_argument("--reg_annot_col", default="PSP: regulatory_site", help="Regulatory annotation column in site table.")
    parser.add_argument("--locprob_col", default=None, help="Localization probability column in site table.")
    parser.add_argument("--locprob_min", type=float, default=0.75, help="Minimum localization probability.")

    parser.add_argument(
        "--gene_to_uniprot",
        default=None,
        help=(
            "Optional gene->UniProt mapping table path or filename. Supports columns "
            "[gene_id,uniprot] or [KEGG_Gene_ID,Uniprot_ID]. If a filename is "
            "provided, it is resolved under MapKinase_WebApp/index_files."
        ),
    )
    parser.add_argument("--weights", default=None, help="JSON object to override scoring weights.")
    parser.add_argument("--max_pathways", type=int, default=None, help="Optional debug limit.")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args(argv)


def infer_sep(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".txt"}:
        return "\t"
    return ","


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing table: {path}")
    sep = infer_sep(path)
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, na_values=[""])
    df.columns = [str(c).strip() for c in df.columns]
    return df


def normalize_uniprot(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return ""
    token = re.split(r"[;,|\s]+", text)[0].strip()
    return token.upper()


def looks_like_uniprot(value: str) -> bool:
    token = normalize_uniprot(value)
    if not token:
        return False
    return bool(UNIPROT_SIMPLE_RE.match(token))


def parse_bool(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"", "0", "false", "f", "no", "n", "none", "nan"}:
        return False
    return True


def safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def clamp01(x: float) -> float:
    if not math.isfinite(x):
        return 0.0
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    return x


def sig_transform(p_value: Any, sig_scale: float) -> float:
    p = safe_float(p_value)
    if not math.isfinite(p) or p <= 0:
        if p == 0:
            p = 1e-300
        else:
            return 0.0
    return clamp01((-math.log10(p)) / sig_scale)


def eff_transform(fc_value: Any, eff_scale: float) -> float:
    fc = safe_float(fc_value)
    if not math.isfinite(fc):
        return 0.0
    if fc == 0:
        return 0.0
    if fc < 0:
        fc = abs(fc)
    if fc <= 0:
        return 0.0
    return clamp01(abs(math.log2(fc)) / eff_scale)


def resolve_column_name(df: pd.DataFrame, requested: Optional[str], fallbacks: Sequence[str] = ()) -> Optional[str]:
    if df is None or df.empty:
        return None
    lower_map = {c.lower(): c for c in df.columns}
    if requested:
        req = requested.strip()
        if req in df.columns:
            return req
        maybe = lower_map.get(req.lower())
        if maybe:
            return maybe
    for fb in fallbacks:
        maybe = lower_map.get(fb.lower())
        if maybe:
            return maybe
    return None


def parse_site_key_cols(raw: Optional[str]) -> List[str]:
    if not raw:
        return []
    return [part.strip() for part in raw.split(",") if part.strip()]


def build_site_key(row_map: Dict[str, Any], site_key_col: Optional[str], site_key_cols: Sequence[str], site_uniprot_col: Optional[str]) -> str:
    if site_key_col and site_key_col in row_map:
        key = str(row_map.get(site_key_col, "") or "").strip()
        if key:
            return key
    if len(site_key_cols) >= 3:
        uni_col, res_col, pos_col = site_key_cols[0], site_key_cols[1], site_key_cols[2]
        uni = normalize_uniprot(row_map.get(uni_col, ""))
        residue = str(row_map.get(res_col, "") or "").strip()
        position = str(row_map.get(pos_col, "") or "").strip()
        if uni and residue and position:
            return f"{uni}_{residue}{position}"
    if len(site_key_cols) == 2:
        uni_col, site_col = site_key_cols[0], site_key_cols[1]
        uni = normalize_uniprot(row_map.get(uni_col, ""))
        site = str(row_map.get(site_col, "") or "").strip()
        if uni and site:
            return f"{uni}_{site}"
    if site_uniprot_col and site_uniprot_col in row_map:
        uni = normalize_uniprot(row_map.get(site_uniprot_col, ""))
        if uni:
            return uni
    return ""


def load_kegg_index(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        index = json.load(fh)
    if not isinstance(index, dict):
        raise ValueError("Pathway index JSON must be an object.")
    for key in ("pathways", "nodes", "edges"):
        if key not in index:
            raise ValueError(f"Pathway index missing key: {key}")
    return index


def load_user_tables(protein_table: Path, site_table: Optional[Path]) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    prot_df = load_table(protein_table)
    site_df = load_table(site_table) if site_table else None
    LOGGER.info("Loaded protein table: %s rows, %s cols", len(prot_df), len(prot_df.columns))
    if site_df is not None:
        LOGGER.info("Loaded site table: %s rows, %s cols", len(site_df), len(site_df.columns))
    return prot_df, site_df


def resolve_index_file_path(raw_path: Optional[str]) -> Optional[Path]:
    if raw_path is None:
        return None
    raw = str(raw_path).strip()
    if not raw:
        return None
    path = Path(raw)
    # For bare filenames, prefer MapKinase_WebApp/index_files.
    if path.parent == Path("."):
        fallback = INDEX_FILES_DIR / raw
        if fallback.exists():
            return fallback
        if path.exists():
            return path
        return fallback
    if path.exists():
        return path
    fallback = INDEX_FILES_DIR / path.name
    if fallback.exists():
        return fallback
    return path


def parse_weights(raw_weights: Optional[str]) -> Dict[str, float]:
    merged = dict(DEFAULT_WEIGHTS)
    if not raw_weights:
        return merged
    try:
        override = json.loads(raw_weights)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid --weights JSON: {exc}") from exc
    if not isinstance(override, dict):
        raise ValueError("--weights must decode to a JSON object.")
    for key, val in override.items():
        if key not in merged:
            LOGGER.warning("Ignoring unknown weight key: %s", key)
            continue
        merged[key] = float(val)
    return merged


def compute_single_protein_scores(
    protein_df: pd.DataFrame,
    site_df: Optional[pd.DataFrame],
    args: argparse.Namespace,
    weights: Dict[str, float],
) -> Dict[str, Dict[str, Any]]:
    protein_id_col = resolve_column_name(protein_df, args.protein_id_col, ("uniprot", "uniprot_id"))
    if not protein_id_col:
        protein_id_col = protein_df.columns[0]
        LOGGER.warning("Protein UniProt column not found; using first column: %s", protein_id_col)

    prot_p_col = resolve_column_name(protein_df, args.p_col_prot)
    prot_fc_col = resolve_column_name(protein_df, args.fc_col_prot)
    phos_p_col = resolve_column_name(protein_df, args.p_col_phospho) if args.p_col_phospho else None
    phos_fc_col = resolve_column_name(protein_df, args.fc_col_phospho) if args.fc_col_phospho else None

    if args.p_col_prot and not prot_p_col:
        LOGGER.warning("Protein p-value column missing: %s (sig(prot_p) set to 0)", args.p_col_prot)
    if args.p_col_phospho and not phos_p_col:
        LOGGER.warning("Phospho aggregate p-value column missing: %s (sig(phospho_p) set to 0)", args.p_col_phospho)

    accum: Dict[str, Dict[str, Any]] = {}

    for row in protein_df.itertuples(index=False, name=None):
        row_map = dict(zip(protein_df.columns, row))
        uniprot = normalize_uniprot(row_map.get(protein_id_col))
        if not uniprot:
            continue
        rec = accum.setdefault(
            uniprot,
            {
                "uniprot": uniprot,
                "reg_site_scores": [],
                "ptm_site_scores": [],
                "prot_sig": 0.0,
                "phospho_sig": 0.0,
                "prot_eff": 0.0,
                "phospho_eff": 0.0,
            },
        )

        prot_sig = sig_transform(row_map.get(prot_p_col), weights["sig_scale"]) if prot_p_col else 0.0
        phos_sig = sig_transform(row_map.get(phos_p_col), weights["sig_scale"]) if phos_p_col else 0.0
        prot_eff = eff_transform(row_map.get(prot_fc_col), weights["eff_scale"]) if prot_fc_col else 0.0
        phos_eff = eff_transform(row_map.get(phos_fc_col), weights["eff_scale"]) if phos_fc_col else 0.0

        rec["prot_sig"] = max(rec["prot_sig"], prot_sig)
        rec["phospho_sig"] = max(rec["phospho_sig"], phos_sig)
        rec["prot_eff"] = max(rec["prot_eff"], prot_eff)
        rec["phospho_eff"] = max(rec["phospho_eff"], phos_eff)

    if site_df is not None and not site_df.empty:
        site_uniprot_col = resolve_column_name(site_df, args.site_uniprot_col, ("uniprot", "uniprot_id", "uniprot id"))
        if not site_uniprot_col:
            site_uniprot_col = site_df.columns[0]
            LOGGER.warning("Site UniProt column not found; using first column: %s", site_uniprot_col)

        site_key_col = resolve_column_name(site_df, args.site_key_col)
        site_key_cols = parse_site_key_cols(args.site_key_cols)
        site_key_cols = [resolve_column_name(site_df, col) or col for col in site_key_cols]
        if site_key_cols and any(col not in site_df.columns for col in site_key_cols):
            LOGGER.warning("One or more --site_key_cols were not found in site table; will fallback when possible.")

        reg_col = resolve_column_name(site_df, args.reg_annot_col) if args.reg_annot_col else None
        locprob_col = resolve_column_name(site_df, args.locprob_col) if args.locprob_col else None
        site_p_col = resolve_column_name(site_df, args.p_col_site or args.p_col_phospho)
        site_fc_col = resolve_column_name(site_df, args.fc_col_site or args.fc_col_phospho or args.fc_col_prot)

        if (args.p_col_site or args.p_col_phospho) and not site_p_col:
            LOGGER.warning("Site p-value column missing; site sig(p) defaults to 0.")
        if not site_fc_col:
            LOGGER.warning("Site fold-change column missing; site eff(fc) defaults to 0.")
        if args.locprob_col and not locprob_col:
            LOGGER.warning("Localization column missing: %s (no localization filter applied).", args.locprob_col)

        for row in site_df.itertuples(index=False, name=None):
            row_map = dict(zip(site_df.columns, row))
            uniprot = normalize_uniprot(row_map.get(site_uniprot_col))
            if not uniprot:
                continue

            if locprob_col:
                locprob = safe_float(row_map.get(locprob_col))
                if not math.isfinite(locprob) or locprob < args.locprob_min:
                    continue
            else:
                locprob = math.nan

            sig_site = sig_transform(row_map.get(site_p_col), weights["sig_scale"]) if site_p_col else 0.0
            eff_site = eff_transform(row_map.get(site_fc_col), weights["eff_scale"]) if site_fc_col else 0.0
            combined = 0.8 * sig_site + 0.2 * eff_site
            if combined <= 0:
                continue

            rec = accum.setdefault(
                uniprot,
                {
                    "uniprot": uniprot,
                    "reg_site_scores": [],
                    "ptm_site_scores": [],
                    "prot_sig": 0.0,
                    "phospho_sig": 0.0,
                    "prot_eff": 0.0,
                    "phospho_eff": 0.0,
                },
            )
            site_key = build_site_key(
                row_map=row_map,
                site_key_col=site_key_col,
                site_key_cols=site_key_cols,
                site_uniprot_col=site_uniprot_col,
            )
            is_reg = parse_bool(row_map.get(reg_col)) if reg_col else False
            site_payload = {
                "site_key": site_key or uniprot,
                "sig": sig_site,
                "eff": eff_site,
                "combined": combined,
                "p_value": row_map.get(site_p_col) if site_p_col else None,
                "fold_change": row_map.get(site_fc_col) if site_fc_col else None,
                "locprob": None if not math.isfinite(locprob) else locprob,
            }
            if is_reg:
                site_payload["contribution"] = weights["w_ann"] * combined
                rec["reg_site_scores"].append(site_payload)
            else:
                site_payload["contribution"] = weights["ptm_site_scale"] * combined
                rec["ptm_site_scores"].append(site_payload)

    site_top_k = int(weights["site_top_k"])
    protein_scores: Dict[str, Dict[str, Any]] = {}
    for uniprot, rec in accum.items():
        reg_sites = sorted(
            rec["reg_site_scores"],
            key=lambda x: (-float(x.get("contribution", 0.0)), str(x.get("site_key", ""))),
        )
        ptm_sites = sorted(
            rec["ptm_site_scores"],
            key=lambda x: (-float(x.get("contribution", 0.0)), str(x.get("site_key", ""))),
        )

        reg_evidence = sum(float(x["contribution"]) for x in reg_sites[:site_top_k])
        ptm_evidence = sum(float(x["contribution"]) for x in ptm_sites[:site_top_k])
        ab_evidence = 0.5 * float(rec["prot_sig"]) + 0.5 * float(rec["phospho_sig"])

        single_score = reg_evidence + weights["ptm_weight"] * ptm_evidence + weights["epsilon"] * ab_evidence
        has_reg = reg_evidence >= weights["reg_gate"]

        protein_scores[uniprot] = {
            "uniprot": uniprot,
            "single_score": float(single_score),
            "reg_evidence": float(reg_evidence),
            "ptm_evidence": float(ptm_evidence),
            "ab_evidence": float(ab_evidence),
            "has_reg": bool(has_reg),
            "prot_sig": float(rec["prot_sig"]),
            "phospho_sig": float(rec["phospho_sig"]),
            "prot_eff": float(rec["prot_eff"]),
            "phospho_eff": float(rec["phospho_eff"]),
            "top_reg_sites": reg_sites[: max(3, site_top_k)],
        }

    LOGGER.info("Computed SingleProteinScore for %s proteins", len(protein_scores))
    return protein_scores


def extract_numeric_gene_id(value: str) -> str:
    token = str(value).strip()
    if not token:
        return ""
    if ":" in token:
        token = token.split(":", 1)[1]
    token = token.strip()
    m = re.match(r"^(\d+)$", token)
    return m.group(1) if m else ""


def build_gene_to_uniprot_map(mapping_path: Optional[Path]) -> Dict[str, List[str]]:
    if mapping_path is None:
        return {}
    if not mapping_path.exists():
        LOGGER.warning("gene_to_uniprot mapping file not found: %s", mapping_path)
        return {}

    df = load_table(mapping_path)
    gene_col = resolve_column_name(df, "gene_id", ("KEGG_Gene_ID",))
    uni_col = resolve_column_name(df, "uniprot", ("Uniprot_ID",))

    if gene_col is None:
        LOGGER.warning("No gene_id/KEGG_Gene_ID column in %s; node mapping may be sparse.", mapping_path)
        return {}
    if uni_col is None:
        LOGGER.warning("No uniprot/Uniprot_ID column in %s; node mapping may be sparse.", mapping_path)
        return {}

    mapping: Dict[str, set[str]] = {}
    for row in df.itertuples(index=False, name=None):
        row_map = dict(zip(df.columns, row))
        uni = normalize_uniprot(row_map.get(uni_col))
        if not uni:
            continue
        raw_gene = str(row_map.get(gene_col, "") or "").strip()
        for part in re.split(r"[,\s;+|]+", raw_gene):
            gid = extract_numeric_gene_id(part)
            if not gid:
                continue
            mapping.setdefault(gid, set()).add(uni)

    fixed = {gid: sorted(unis) for gid, unis in mapping.items()}
    LOGGER.info("Loaded gene->UniProt map: %s gene IDs", len(fixed))
    return fixed


def build_protein_lookup(protein_scores: Dict[str, Dict[str, Any]]) -> ProteinLookup:
    by_base: Dict[str, Dict[str, Any]] = {}
    for uid, rec in protein_scores.items():
        base = uid.split("-", 1)[0]
        prev = by_base.get(base)
        if prev is None:
            by_base[base] = rec
            continue
        prev_score = float(prev.get("single_score", 0.0))
        cur_score = float(rec.get("single_score", 0.0))
        if cur_score > prev_score or (cur_score == prev_score and uid < str(prev.get("uniprot", ""))):
            by_base[base] = rec
    return ProteinLookup(exact=protein_scores, by_base=by_base)


def candidate_uniprots_for_node(node_obj: Dict[str, Any], gene_to_uniprot: Dict[str, List[str]]) -> List[str]:
    candidates = node_obj.get("candidates", {}) if isinstance(node_obj, dict) else {}
    kegg_genes = candidates.get("kegg_genes", []) if isinstance(candidates, dict) else []
    gene_ids = candidates.get("gene_ids", []) if isinstance(candidates, dict) else []
    pre_mapped_unis = candidates.get("uniprot", []) if isinstance(candidates, dict) else []

    out: set[str] = set()

    for uni in pre_mapped_unis if isinstance(pre_mapped_unis, list) else []:
        normalized = normalize_uniprot(uni)
        if normalized:
            out.add(normalized)

    all_gene_ids: set[str] = set()
    if isinstance(gene_ids, list):
        for g in gene_ids:
            gid = extract_numeric_gene_id(g)
            if gid:
                all_gene_ids.add(gid)
            if looks_like_uniprot(str(g)):
                out.add(normalize_uniprot(g))

    if isinstance(kegg_genes, list):
        for kg in kegg_genes:
            gid = extract_numeric_gene_id(kg)
            if gid:
                all_gene_ids.add(gid)
            if looks_like_uniprot(str(kg)):
                out.add(normalize_uniprot(kg))

    for gid in sorted(all_gene_ids):
        for uni in gene_to_uniprot.get(gid, []):
            normalized = normalize_uniprot(uni)
            if normalized:
                out.add(normalized)

    return sorted(out)


def resolve_node_scores(
    index_nodes: Dict[str, Dict[str, Any]],
    protein_lookup: ProteinLookup,
    gene_to_uniprot: Dict[str, List[str]],
) -> Dict[str, Dict[str, Any]]:
    node_state: Dict[str, Dict[str, Any]] = {}

    for node_id, node_obj in index_nodes.items():
        candidate_unis = candidate_uniprots_for_node(node_obj, gene_to_uniprot)
        candidate_records: List[Tuple[str, Dict[str, Any]]] = []
        for uni in candidate_unis:
            rec = protein_lookup.get(uni)
            if rec is not None:
                candidate_records.append((uni, rec))

        rep_uniprot = None
        rep_record: Optional[Dict[str, Any]] = None
        best_score = 0.0
        for uni, rec in candidate_records:
            score = float(rec.get("single_score", 0.0))
            if rep_record is None or score > best_score or (score == best_score and str(rec.get("uniprot", uni)) < str(rep_uniprot)):
                rep_uniprot = str(rec.get("uniprot", uni))
                rep_record = rec
                best_score = score

        node_has_reg = any(bool(rec.get("has_reg", False)) for _, rec in candidate_records)

        node_state[node_id] = {
            "node_id": node_id,
            "node_score": float(best_score if rep_record is not None else 0.0),
            "node_has_reg": bool(node_has_reg),
            "rep_uniprot": rep_uniprot,
            "rep_has_reg": bool(rep_record.get("has_reg", False)) if rep_record else False,
            "rep_top_reg_sites": (rep_record.get("top_reg_sites", []) if rep_record else []),
            "rep_score": float(rep_record.get("single_score", 0.0)) if rep_record else 0.0,
            "candidate_uniprot_count": len(candidate_unis),
            "present_candidate_count": len(candidate_records),
            "rep_reason": "max_single_protein_score" if rep_record else "no_mapped_candidates",
        }

    return node_state


def _node_edge_payload(state: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "representative_uniprot": state.get("rep_uniprot"),
        "node_score": float(state.get("node_score", 0.0)),
        "node_has_reg": bool(state.get("node_has_reg", False)),
        "representative_has_reg": bool(state.get("rep_has_reg", False)),
        "top_reg_sites": state.get("rep_top_reg_sites", []),
    }


def score_pathway(
    pathway: Dict[str, Any],
    node_state: Dict[str, Dict[str, Any]],
    weights: Dict[str, float],
    pathway_source: str = "kegg",
) -> Dict[str, Any]:
    pathway_id = str(pathway.get("pathway_id", ""))
    node_ids: List[str] = list(pathway.get("nodes", []))
    pairs1: List[List[Any]] = list(pathway.get("pairs1", []))
    pairs2: List[List[Any]] = list(pathway.get("pairs2", []))
    node_count = int(pathway.get("node_count", len(node_ids)))
    edge_count = int(pathway.get("edge_count", len(pathway.get("edges", []))))

    conn1 = 0.0
    conn2 = 0.0
    top1: List[Dict[str, Any]] = []
    top2: List[Dict[str, Any]] = []
    top_n = int(weights["top_edges_n"])

    for pair in pairs1:
        if len(pair) < 2:
            continue
        a, b = str(pair[0]), str(pair[1])
        sa = node_state.get(a)
        sb = node_state.get(b)
        if sa is None or sb is None:
            continue
        if not (sa["node_has_reg"] and sb["node_has_reg"]):
            continue
        contrib = float(sa["node_score"]) * float(sb["node_score"])
        if contrib <= 0:
            continue
        conn1 += contrib
        top1.append(
            {
                "node_a": a,
                "node_b": b,
                "contribution": contrib,
                "node_a_details": _node_edge_payload(sa),
                "node_b_details": _node_edge_payload(sb),
            }
        )

    for pair in pairs2:
        if len(pair) < 3:
            continue
        a, b = str(pair[0]), str(pair[1])
        bridge_count = safe_float(pair[2])
        if not math.isfinite(bridge_count) or bridge_count <= 0:
            continue
        sa = node_state.get(a)
        sb = node_state.get(b)
        if sa is None or sb is None:
            continue
        if not (sa["node_has_reg"] and sb["node_has_reg"]):
            continue

        bridge_weight = weights["two_hop_base"] * math.log1p(bridge_count)
        contrib = float(sa["node_score"]) * float(sb["node_score"]) * bridge_weight
        if contrib <= 0:
            continue
        conn2 += contrib
        top2.append(
            {
                "node_a": a,
                "node_b": b,
                "bridge_count": int(bridge_count),
                "bridge_weight": bridge_weight,
                "contribution": contrib,
                "node_a_details": _node_edge_payload(sa),
                "node_b_details": _node_edge_payload(sb),
            }
        )

    top1 = sorted(top1, key=lambda x: (-x["contribution"], x["node_a"], x["node_b"]))[:top_n]
    top2 = sorted(top2, key=lambda x: (-x["contribution"], x["node_a"], x["node_b"]))[:top_n]

    node_scores = [float(node_state.get(nid, {}).get("node_score", 0.0)) for nid in node_ids]
    non_zero_scores = [s for s in node_scores if s > 0]
    if non_zero_scores:
        top_k = int(weights["node_mass_top_k"])
        node_mass = sum(sorted(non_zero_scores, reverse=True)[:top_k]) / min(top_k, len(non_zero_scores))
    else:
        node_mass = 0.0

    alpha = float(weights["alpha"])
    norm = float(node_count ** alpha) if node_count > 0 else 1.0
    connection_score = (conn1 + weights["conn2_weight"] * conn2) / norm if norm > 0 else 0.0
    final_score = connection_score + weights["node_mass_weight"] * node_mass

    scored_node_count = sum(1 for nid in node_ids if float(node_state.get(nid, {}).get("node_score", 0.0)) > 0)
    reg_node_count = sum(1 for nid in node_ids if bool(node_state.get(nid, {}).get("node_has_reg", False)))
    mapped_node_count = sum(1 for nid in node_ids if int(node_state.get(nid, {}).get("present_candidate_count", 0)) > 0)

    return {
        "pathway_source": str(pathway_source or "kegg").strip().lower(),
        "pathway_id": pathway_id,
        "name": pathway.get("name", pathway_id),
        "final_score": float(final_score),
        "connection_score": float(connection_score),
        "conn1": float(conn1),
        "conn2": float(conn2),
        "node_mass": float(node_mass),
        "node_count": node_count,
        "edge_count": edge_count,
        "scored_node_count": int(scored_node_count),
        "reg_node_count": int(reg_node_count),
        "mapped_node_count": int(mapped_node_count),
        "top_edges_1hop": top1,
        "top_edges_2hop": top2,
    }


def rank_all_pathways(
    kegg_index: Dict[str, Any],
    node_state: Dict[str, Dict[str, Any]],
    weights: Dict[str, float],
    max_pathways: Optional[int] = None,
    pathway_source: Optional[str] = None,
) -> List[Dict[str, Any]]:
    pathways: List[Dict[str, Any]] = list(kegg_index.get("pathways", []))
    if max_pathways is not None:
        pathways = pathways[:max_pathways]
    resolved_source = str(pathway_source or kegg_index.get("meta", {}).get("pathway_source", "kegg")).strip().lower() or "kegg"

    results: List[Dict[str, Any]] = []
    total = len(pathways)
    for idx, pathway in enumerate(pathways, start=1):
        if idx % 25 == 0 or idx == total:
            LOGGER.info("Scoring %s pathway %s/%s", resolved_source, idx, total)
        results.append(score_pathway(pathway, node_state=node_state, weights=weights, pathway_source=resolved_source))

    results.sort(key=lambda x: (-x["final_score"], x["pathway_source"], x["pathway_id"]))
    return results


def write_output(results: List[Dict[str, Any]], out_path: Path, fmt: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "json":
        with out_path.open("w", encoding="utf-8") as fh:
            json.dump(results, fh, ensure_ascii=False, indent=2)
        return

    flat_rows: List[Dict[str, Any]] = []
    for row in results:
        flat_rows.append(
            {
                "pathway_source": row.get("pathway_source", "kegg"),
                "pathway_id": row["pathway_id"],
                "name": row["name"],
                "final_score": row["final_score"],
                "connection_score": row["connection_score"],
                "conn1": row["conn1"],
                "conn2": row["conn2"],
                "node_mass": row["node_mass"],
                "node_count": row["node_count"],
                "edge_count": row["edge_count"],
                "scored_node_count": row["scored_node_count"],
                "reg_node_count": row["reg_node_count"],
                "mapped_node_count": row["mapped_node_count"],
                "top_edges_1hop": json.dumps(row.get("top_edges_1hop", []), ensure_ascii=False),
                "top_edges_2hop": json.dumps(row.get("top_edges_2hop", []), ensure_ascii=False),
            }
        )
    df = pd.DataFrame(flat_rows)
    sep = infer_sep(out_path)
    df.to_csv(out_path, sep=sep, index=False)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s | %(levelname)s | %(message)s")

    weights = parse_weights(args.weights)
    if args.max_pathways is not None and args.max_pathways < 0:
        raise ValueError("--max_pathways must be >= 0")

    kegg_index_path = resolve_index_file_path(args.kegg_index) if args.kegg_index else None
    wp_index_path = resolve_index_file_path(args.wikipathways_index) if args.wikipathways_index else None
    if kegg_index_path is None and wp_index_path is None:
        raise ValueError("At least one index must be provided: --kegg_index and/or --wikipathways_index.")
    gene_map_path = resolve_index_file_path(args.gene_to_uniprot) if args.gene_to_uniprot else None
    LOGGER.info("Using index directory fallback: %s", INDEX_FILES_DIR)
    if kegg_index_path is not None:
        LOGGER.info("Resolved KEGG index path: %s", kegg_index_path)
    if wp_index_path is not None:
        LOGGER.info("Resolved WikiPathways index path: %s", wp_index_path)
    if gene_map_path is not None:
        LOGGER.info("Resolved gene->UniProt map path: %s", gene_map_path)

    indices_to_score: List[Tuple[str, Dict[str, Any]]] = []
    if kegg_index_path is not None:
        indices_to_score.append(("kegg", load_kegg_index(kegg_index_path)))
    if wp_index_path is not None:
        indices_to_score.append(("wikipathways", load_kegg_index(wp_index_path)))

    protein_df, site_df = load_user_tables(Path(args.protein_table), Path(args.site_table) if args.site_table else None)

    protein_scores = compute_single_protein_scores(protein_df=protein_df, site_df=site_df, args=args, weights=weights)
    protein_lookup = build_protein_lookup(protein_scores)

    gene_map = build_gene_to_uniprot_map(gene_map_path)
    ranked: List[Dict[str, Any]] = []
    for source_key, pathway_index in indices_to_score:
        node_state = resolve_node_scores(
            index_nodes=pathway_index.get("nodes", {}),
            protein_lookup=protein_lookup,
            gene_to_uniprot=gene_map,
        )
        ranked.extend(
            rank_all_pathways(
                kegg_index=pathway_index,
                node_state=node_state,
                weights=weights,
                max_pathways=args.max_pathways,
                pathway_source=source_key,
            )
        )
    ranked.sort(key=lambda x: (-x.get("final_score", 0.0), str(x.get("pathway_source", "")), str(x.get("pathway_id", ""))))
    write_output(ranked, Path(args.out), args.format)

    if args.max_pathways is not None:
        LOGGER.info("Debug subset mode: scored first %s pathways.", args.max_pathways)
    LOGGER.info("Done. Wrote %s ranked pathways to %s", len(ranked), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
