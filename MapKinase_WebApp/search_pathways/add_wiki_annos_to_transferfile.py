#!/usr/bin/env python3
"""
Add UniProt cross-reference annotations to a tabular file of UniProt IDs.
"""

from __future__ import annotations

import argparse
import csv
import time
from typing import Dict, List, Optional, Tuple

import requests


INPUT_FILE = ""
OUTPUT_FILE = ""
INPUT_DELIMITER = "\t"


DEFAULT_COLUMNS = [
    "BRENDA",
    "Complex Portal",
    "EMBL",
    "EcoGene",
    "Ensembl",
    "Entrez Gene",
    "Enzyme Nomenclature",
    "GenBank",
    "HGNC",
    "HGNC Accession number",
    "HomoloGene",
    "InterPro",
    "KEGG Genes",
    "NCBI Protein",
    "Pfam",
    "Reactome",
    "RefSeq",
    "Rfam",
    "Uniprot-TrEMBL",
    "WikiPathways",
    "Wikidata",
    "miRBase Sequence",
    "miRBase mature sequence",
    "pato",
]


DB_MAP = {
    "BRENDA": "BRENDA",
    "Complex Portal": "ComplexPortal",
    "EMBL": "EMBL",
    "EcoGene": "EcoGene",
    "Ensembl": "Ensembl",
    "Entrez Gene": "GeneID",
    "Enzyme Nomenclature": "EnzymeNomenclature",
    "GenBank": "GenBank",
    "HGNC": "HGNC",
    "HGNC Accession number": "HGNC",
    "HomoloGene": "HomoloGene",
    "InterPro": "InterPro",
    "KEGG Genes": "KEGG",
    "NCBI Protein": "NCBI Protein",
    "Pfam": "Pfam",
    "Reactome": "Reactome",
    "RefSeq": "RefSeq",
    "Rfam": "Rfam",
    "Uniprot-TrEMBL": "UniProtKB-TrEMBL",
    "WikiPathways": "WikiPathways",
    "Wikidata": "Wikidata",
    "miRBase Sequence": "miRBase",
    "miRBase mature sequence": "miRBase",
    "pato": "PATO",
}


def _pick_uniprot_column(headers: List[str], hint: Optional[str]) -> int:
    if hint:
        if hint.isdigit():
            idx = int(hint)
            if 0 <= idx < len(headers):
                return idx
        else:
            for i, h in enumerate(headers):
                if h.strip().lower() == hint.strip().lower():
                    return i
    for i, h in enumerate(headers):
        if "uniprot" in h.strip().lower():
            return i
    return 0


def _fetch_uniprot_json(session: requests.Session, accession: str) -> Optional[Dict]:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = session.get(url, timeout=30)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        return resp.json()
    except requests.RequestException:
        return None


def _fetch_uniprot_batch(session: requests.Session, accessions: List[str]) -> Dict[str, Dict]:
    if not accessions:
        return {}
    query = " OR ".join(f"accession:{acc}" for acc in accessions if acc)
    if not query:
        return {}
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": query, "format": "json"}
    try:
        resp = session.get(url, params=params, timeout=60)
        resp.raise_for_status()
        payload = resp.json()
    except requests.RequestException:
        return {}
    results = payload.get("results") or []
    out: Dict[str, Dict] = {}
    for item in results:
        acc = item.get("primaryAccession")
        if acc:
            out[acc] = item
    return out


def _collect_xrefs(payload: Dict) -> Dict[str, List[str]]:
    xrefs: Dict[str, List[str]] = {}
    items = payload.get("uniProtKBCrossReferences") or []
    for item in items:
        db = item.get("database")
        if not db:
            continue
        xrefs.setdefault(db, []).append(item.get("id") or "")
    return xrefs


def _collect_mirbase(payload: Dict) -> Tuple[List[str], List[str]]:
    seqs: List[str] = []
    mats: List[str] = []
    items = payload.get("uniProtKBCrossReferences") or []
    for item in items:
        if item.get("database") != "miRBase":
            continue
        props = item.get("properties") or {}
        for key, value in props.items():
            key_l = str(key).strip().lower()
            if "mature" in key_l:
                mats.append(str(value))
            elif "sequence" in key_l:
                seqs.append(str(value))
        if not props:
            seqs.append(item.get("id") or "")
    return seqs, mats


def _sanitize_hgnc_accession(values: List[str]) -> List[str]:
    cleaned: List[str] = []
    for val in values:
        if val.startswith("HGNC:"):
            cleaned.append(val.split(":", 1)[-1])
        else:
            cleaned.append(val)
    return cleaned


def _join(values: List[str]) -> str:
    return ",".join(sorted({v for v in values if v}))


def build_annotations(payload: Dict) -> Dict[str, str]:
    annotations: Dict[str, str] = {col: "" for col in DEFAULT_COLUMNS}
    if not payload:
        return annotations

    xrefs = _collect_xrefs(payload)
    for col, db in DB_MAP.items():
        if col.startswith("miRBase"):
            continue
        vals = xrefs.get(db, [])
        if col == "HGNC Accession number":
            vals = _sanitize_hgnc_accession(vals)
        annotations[col] = _join(vals)

    mir_seq, mir_mat = _collect_mirbase(payload)
    annotations["miRBase Sequence"] = _join(mir_seq)
    annotations["miRBase mature sequence"] = _join(mir_mat)
    return annotations


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate a file of UniProt IDs via UniProt REST.")
    parser.add_argument("--input", default=INPUT_FILE, help="Input file path.")
    parser.add_argument("--output", default=OUTPUT_FILE, help="Output file path.")
    parser.add_argument("--delimiter", default=INPUT_DELIMITER, help="Input/output delimiter (default: tab).")
    parser.add_argument("--uniprot-column", default="", help="Column name or 0-based index.")
    parser.add_argument("--sleep", type=float, default=0.0, help="Seconds to sleep between requests.")
    parser.add_argument("--batch-size", type=int, default=50, help="Accessions per batch request.")
    args = parser.parse_args()

    if not args.input:
        raise SystemExit("Input file path is required (set INPUT_FILE or use --input).")
    if not args.output:
        raise SystemExit("Output file path is required (set OUTPUT_FILE or use --output).")

    with open(args.input, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=args.delimiter)
        rows = list(reader)

    if not rows:
        raise SystemExit("Input file is empty.")

    headers = rows[0]
    data_rows = rows[1:]
    col_idx = _pick_uniprot_column(headers, args.uniprot_column)

    out_headers = list(headers)
    for col in DEFAULT_COLUMNS:
        if col not in out_headers:
            out_headers.append(col)

    session = requests.Session()
    cache: Dict[str, Dict[str, str]] = {}

    accessions = []
    for row in data_rows:
        if len(row) > col_idx:
            acc = row[col_idx].strip()
            if acc and acc not in accessions:
                accessions.append(acc)

    if args.batch_size > 0:
        for i in range(0, len(accessions), args.batch_size):
            batch = accessions[i : i + args.batch_size]
            batch_payloads = _fetch_uniprot_batch(session, batch)
            for acc in batch:
                payload = batch_payloads.get(acc)
                if payload:
                    cache[acc] = build_annotations(payload)
            if args.sleep:
                time.sleep(args.sleep)

    output_rows: List[List[str]] = [out_headers]
    for row in data_rows:
        while len(row) < len(headers):
            row.append("")
        accession = row[col_idx].strip()
        if accession:
            if accession not in cache:
                payload = _fetch_uniprot_json(session, accession)
                cache[accession] = build_annotations(payload or {})
                if args.sleep:
                    time.sleep(args.sleep)
            annos = cache[accession]
        else:
            annos = {col: "" for col in DEFAULT_COLUMNS}

        out_row = list(row)
        for col in DEFAULT_COLUMNS:
            out_row.append(annos.get(col, ""))
        output_rows.append(out_row)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=args.delimiter)
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
