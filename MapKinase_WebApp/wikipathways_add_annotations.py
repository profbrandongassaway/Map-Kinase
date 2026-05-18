#!/usr/bin/env python3
"""
Add UniProt cross-reference annotations to a tabular file of UniProt IDs.
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from typing import Dict, List, Optional, Tuple

import requests

ORG = "sce"
INPUT_FILE = f"proteomes/{ORG}_KEGG_Conversion.txt"
OUTPUT_FILE = f"id_mapping_tables/{ORG}_mapping_table.txt"
INPUT_DELIMITER = "\t"


DEFAULT_COLUMNS = [
    "BRENDA",
    "Complex Portal",
    "EMBL",
    "EcoGene",
    "Ensembl_Transcript",
    "Ensembl_Protein",
    "Ensembl_Gene",
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
    "Ensembl_Transcript": "Ensembl",
    "Ensembl_Protein": "Ensembl",
    "Ensembl_Gene": "Ensembl",
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

def _props_to_dict(props):
    """
    UniProt JSON 'properties' can be either:
      - a list of {"key": ..., "value": ...} dicts (most common), OR
      - a dict of key->value (less common)

    Normalize to dict[str, str].
    """
    if not props:
        return {}
    if isinstance(props, dict):
        return {str(k): str(v) for k, v in props.items()}
    if isinstance(props, list):
        out = {}
        for p in props:
            if isinstance(p, dict) and "key" in p and "value" in p:
                out[str(p["key"])] = str(p["value"])
        return out
    return {}




def _collect_ensembl_ids(payload: Dict) -> Tuple[List[str], List[str], List[str]]:
    """
    Extract Ensembl Gene/Transcript/Protein IDs from UniProt 'Ensembl' cross-references.

    Returns:
      (genes_ENSG, transcripts_ENST, proteins_ENSP)

    UniProt is protein-centric, so the main Ensembl xref 'id' is often ENST or ENSP.
    ENSG frequently appears in 'properties' (GeneId / gene ID).
    """
    genes: List[str] = []
    transcripts: List[str] = []
    proteins: List[str] = []

    items = payload.get("uniProtKBCrossReferences") or []
    for item in items:
        if item.get("database") != "Ensembl":
            continue

        # Sometimes the main ID is an ENST or ENSP (occasionally ENSG)
        main_id = item.get("id") or ""
        if main_id.startswith("ENSG"):
            genes.append(main_id)
        elif main_id.startswith("ENST"):
            transcripts.append(main_id)
        elif main_id.startswith("ENSP"):
            proteins.append(main_id)

        props = _props_to_dict(item.get("properties"))
        for k, v in props.items():
            kl = k.lower().replace(" ", "").replace("_", "")

            # Common-ish variants UniProt uses in Ensembl xrefs
            if ("geneid" in kl or "gene" == kl) and v.startswith("ENSG"):
                genes.append(v)
            elif ("transcriptid" in kl or "transcript" == kl) and v.startswith("ENST"):
                transcripts.append(v)
            elif ("proteinid" in kl or "protein" == kl) and v.startswith("ENSP"):
                proteins.append(v)

    return genes, transcripts, proteins



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
    """
    Batch query using UniProtKB search. Returns mapping primaryAccession -> payload item.

    Change vs your original:
      - Adds 'size' so UniProt returns up to the number of accessions requested
        (otherwise the API often defaults to a small page size).
    """
    if not accessions:
        return {}
    query = " OR ".join(f"accession:{acc}" for acc in accessions if acc)
    if not query:
        return {}

    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": query, "format": "json", "size": len(accessions)}

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
    """
    Extract miRBase 'sequence' and 'mature' IDs from UniProt xrefs.
    Handles UniProt 'properties' being list-of-dicts or dict.
    """
    seqs: List[str] = []
    mats: List[str] = []

    items = payload.get("uniProtKBCrossReferences") or []
    for item in items:
        if item.get("database") != "miRBase":
            continue

        props = _props_to_dict(item.get("properties"))
        if props:
            for key, value in props.items():
                key_l = str(key).strip().lower()
                if "mature" in key_l:
                    mats.append(str(value))
                elif "sequence" in key_l:
                    seqs.append(str(value))
        else:
            # Fallback: store the main xref ID if no properties are present
            maybe = item.get("id") or ""
            if maybe:
                seqs.append(maybe)

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
    """
    Build output-column -> comma-joined string annotations for a UniProtKB JSON payload.

    Changes vs your original:
      - Skips Ensembl_* in the DB_MAP loop so we don't overwrite those columns.
      - Adds Ensembl special parsing to populate ENSG/ENST/ENSP correctly.
      - Uses the updated _collect_mirbase() which handles properties list/dict.
    """
    annotations: Dict[str, str] = {col: "" for col in DEFAULT_COLUMNS}
    if not payload:
        return annotations

    xrefs = _collect_xrefs(payload)

    # Fill "simple" xref columns (by database name)
    for col, db in DB_MAP.items():
        if col.startswith("miRBase") or col.startswith("Ensembl_"):
            continue

        vals = xrefs.get(db, [])
        if col == "HGNC Accession number":
            vals = _sanitize_hgnc_accession(vals)

        annotations[col] = _join(vals)

    # Ensembl special handling: split ENSG/ENST/ENSP out of the "Ensembl" xref
    ens_g, ens_t, ens_p = _collect_ensembl_ids(payload)
    annotations["Ensembl_Gene"] = _join(ens_g)
    annotations["Ensembl_Transcript"] = _join(ens_t)
    annotations["Ensembl_Protein"] = _join(ens_p)

    # miRBase special handling (properties structure varies)
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
    parser.add_argument("--batch-size", type=int, default=100, help="Accessions per batch request (default: 100).")
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

    # Build unique accession list (preserve first-seen order)
    accessions: List[str] = []
    seen = set()
    for row in data_rows:
        if len(row) > col_idx:
            acc = row[col_idx].strip()
            if acc and acc not in seen:
                accessions.append(acc)
                seen.add(acc)

    total = len(accessions)
    if total == 0:
        raise SystemExit("No UniProt accessions found in the selected column.")

    # -------------------------
    # Batch fetch with progress
    # -------------------------
    if args.batch_size > 0:
        total_batches = math.ceil(total / args.batch_size)
        for batch_idx, start in enumerate(range(0, total, args.batch_size), start=1):
            batch = accessions[start : start + args.batch_size]

            batch_payloads = _fetch_uniprot_batch(session, batch)
            for acc in batch:
                payload = batch_payloads.get(acc)
                if payload:
                    cache[acc] = build_annotations(payload)
                # Optional: leave uncached here; it can fall back to per-ID later if you want.
                # else:
                #     print(f"[batch {batch_idx}] Missing in batch response: {acc}")

            processed = min(start + len(batch), total)
            remaining_batches = total_batches - batch_idx
            print(
                f"Batch {batch_idx}/{total_batches} done "
                f"({processed}/{total} IDs processed). "
                f"{remaining_batches} batch queries left."
            )

            if args.sleep:
                time.sleep(args.sleep)

    # -------------------------
    # Build output rows
    # -------------------------
    output_rows: List[List[str]] = [out_headers]
    for row in data_rows:
        while len(row) < len(headers):
            row.append("")

        accession = row[col_idx].strip()
        if accession:
            if accession not in cache:
                # Fallback: per-accession query if it wasn't returned in batch results
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
