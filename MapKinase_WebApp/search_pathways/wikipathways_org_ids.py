#!/usr/bin/env python3
"""
Download organism-wide WikiPathways GPML zip to disk, iterate GPML files inside the zip,
parse GeneProduct Xrefs (Database/ID), and write a dynamic-column TSV.

Requires: pip install requests
"""

from __future__ import annotations

import csv
import re
import sys
import zipfile
from pathlib import Path
from typing import Dict, Set, Optional, List, Tuple

import requests
import xml.etree.ElementTree as ET


# ==========================
# USER SETTINGS
# ==========================
ORGANISM = "Homo sapiens"  # e.g., "Mus musculus"
OUT_TSV = "wikipathways_geneproduct_ids.tsv"
OUT_IDS = "wikipathways_matching_ids.txt"

# Where to save the downloaded zip (folder or full path):
DOWNLOAD_DIR = "."  # set to e.g. r"C:\\data\\wikipathways"
# If you set ZIP_PATH explicitly, it will use that (and skip auto-naming):
ZIP_PATH: Optional[str] = None

INCLUDE_PATHWAY_NAME = False  # adds PathwayName column after PathwayID
TIMEOUT_S = 60

# Search behavior (empty list means no filtering)
SEARCH_TERMS = []  # e.g., ["MAPK", "ERK1"]
SEARCH_IN = "all"  # "pathway_name", "geneproduct_xref", "text_label", "all"
CASE_INSENSITIVE = True


# ==========================
# CONSTANTS
# ==========================
DATA_BASE = "https://data.wikipathways.org"
CURRENT_GPML_DIR = f"{DATA_BASE}/current/gpml/"


# ==========================
# HELPERS
# ==========================
def organism_to_slug(org: str) -> str:
    return re.sub(r"\s+", "_", org.strip())


def fetch_current_gpml_zip_url_for_species(species_slug: str) -> Tuple[str, str]:
    """
    Parse the current GPML directory listing and locate:
      wikipathways-YYYYMMDD-gpml-<species_slug>.zip
    """
    r = requests.get(CURRENT_GPML_DIR, timeout=TIMEOUT_S)
    r.raise_for_status()

    html = r.text
    pattern = rf'wikipathways-(\d{{8}})-gpml-{re.escape(species_slug)}\.zip'
    m = re.search(pattern, html)
    if not m:
        raise RuntimeError(
            f"Could not find a GPML zip for organism '{species_slug}' at:\n{CURRENT_GPML_DIR}\n"
            f"Double-check the organism name (spelling/case)."
        )

    date = m.group(1)
    filename = f"wikipathways-{date}-gpml-{species_slug}.zip"
    return date, f"{CURRENT_GPML_DIR}{filename}"


def stream_download(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=TIMEOUT_S) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)


def extract_wpid_from_zip_member(member: str) -> Optional[str]:
    base = member.split("/")[-1]
    m = re.search(r"(WP\d+)", base, flags=re.IGNORECASE)
    return m.group(1).upper() if m else None


def safe_strip(s: Optional[str]) -> Optional[str]:
    if s is None:
        return None
    s2 = s.strip()
    return s2 if s2 else None


def local_name(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def parse_geneproducts_from_gpml(
    gpml_bytes: bytes,
) -> Tuple[Optional[str], Dict[str, Set[str]], Set[str], Set[str]]:
    """
    Returns (pathway_name, db_to_ids, text_labels, xref_ids) where:
      db_to_ids[Database] = set of IDs (from Xref) for all DataNode types
    """
    root = ET.fromstring(gpml_bytes)

    pathway_name = None
    if local_name(root.tag) == "Pathway":
        pathway_name = safe_strip(root.attrib.get("Name"))

    db_to_ids: Dict[str, Set[str]] = {}
    text_labels: Set[str] = set()
    xref_ids: Set[str] = set()

    for elem in root.iter():
        if local_name(elem.tag) != "DataNode":
            continue
        if elem.attrib.get("Type") != "GeneProduct":
            continue

        label = safe_strip(elem.attrib.get("TextLabel"))
        if label:
            text_labels.add(label)

        found_xref = False
        for child in list(elem):
            if local_name(child.tag) != "Xref":
                continue
            xref_db = safe_strip(child.attrib.get("Database"))
            xref_id = safe_strip(child.attrib.get("ID"))
            if not xref_db or not xref_id:
                continue
            db_to_ids.setdefault(xref_db, set()).add(xref_id)
            xref_ids.add(xref_id)
            found_xref = True
        if not found_xref:
            continue

    return pathway_name, db_to_ids, text_labels, xref_ids


def normalize_text(value: str) -> str:
    return value.lower() if CASE_INSENSITIVE else value


def pathway_matches(
    pathway_name: Optional[str],
    text_labels: Set[str],
    xref_ids: Set[str],
    terms: List[str],
    mode: str,
) -> bool:
    if not terms:
        return True
    haystacks: List[str] = []
    if mode in {"pathway_name", "all"} and pathway_name:
        haystacks.append(pathway_name)
    if mode in {"text_label", "all"} and text_labels:
        haystacks.extend(sorted(text_labels))
    if mode in {"geneproduct_xref", "all"} and xref_ids:
        haystacks.extend(sorted(xref_ids))
    if not haystacks:
        return False
    norm_hay = [normalize_text(h) for h in haystacks if h]
    for term in terms:
        if not term:
            continue
        t = normalize_text(term)
        if any(t in h for h in norm_hay):
            return True
    return False


# ==========================
# MAIN
# ==========================
def main() -> None:
    slug = organism_to_slug(ORGANISM)
    release_date, zip_url = fetch_current_gpml_zip_url_for_species(slug)

    if ZIP_PATH:
        zip_path = Path(ZIP_PATH)
    else:
        zip_filename = f"wikipathways-{release_date}-gpml-{slug}.zip"
        zip_path = Path(DOWNLOAD_DIR) / zip_filename

    print(f"[1/3] Organism: {ORGANISM} ({slug})")
    print(f"      Release date: {release_date}")
    print(f"      Download URL: {zip_url}")
    print(f"      Saving to: {zip_path}")

    if not zip_path.exists() or zip_path.stat().st_size == 0:
        print("[2/3] Downloading zip to disk...")
        stream_download(zip_url, zip_path)
        print(f"      Done. Size: {zip_path.stat().st_size:,} bytes")
    else:
        print("[2/3] Zip already exists on disk — skipping download.")

    print("[3/3] Reading GPML files inside zip and building TSV rows...")
    rows: List[Dict[str, str]] = []
    matched_ids: List[str] = []
    all_dbs: Set[str] = set()
    failures = 0
    parsed_ok = 0
    skipped = 0

    with zipfile.ZipFile(zip_path, "r") as zf:
        members = [i for i in zf.infolist() if (not i.is_dir()) and i.filename.lower().endswith(".gpml")]
        total = len(members)
        print(f"      Found {total:,} GPML files in zip.")

        for idx, info in enumerate(members, start=1):
            wpid = extract_wpid_from_zip_member(info.filename)
            if not wpid:
                continue

            try:
                gpml_bytes = zf.read(info)
                pname, db_to_ids, text_labels, xref_ids = parse_geneproducts_from_gpml(gpml_bytes)
                parsed_ok += 1
                if not pathway_matches(pname, text_labels, xref_ids, SEARCH_TERMS, SEARCH_IN):
                    skipped += 1
                    continue

                all_dbs.update(db_to_ids.keys())

                row: Dict[str, str] = {"PathwayID": wpid}
                if INCLUDE_PATHWAY_NAME:
                    row["PathwayName"] = pname or ""

                for db, ids in db_to_ids.items():
                    row[db] = ",".join(sorted(ids))

                rows.append(row)
                matched_ids.append(wpid)

                if idx % 250 == 0 or idx == total:
                    print(f"      Parsed {idx:,}/{total:,} GPML entries...")

            except Exception as e:
                failures += 1
                print(f"      WARNING: failed to parse {info.filename}: {e}", file=sys.stderr)

    # Write TSV with dynamic columns
    headers = ["PathwayID"]
    if INCLUDE_PATHWAY_NAME:
        headers.append("PathwayName")
    headers.extend(sorted(all_dbs))

    with open(OUT_TSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(headers)
        for r in rows:
            w.writerow([r.get(h, "") for h in headers])

    print(f"\nDone. Wrote TSV: {OUT_TSV}")
    print(f"Summary: parsed={parsed_ok:,}, skipped={skipped:,}, failures={failures:,}.")
    if failures:
        print(f"Note: {failures} GPML files failed to parse (see warnings).")
    if matched_ids:
        with open(OUT_IDS, "w", encoding="utf-8") as f:
            f.write("\n".join(sorted(set(matched_ids))))
        print(f"Done. Wrote matching IDs: {OUT_IDS}")
    else:
        print("No pathway IDs matched the search terms.")


if __name__ == "__main__":
    main()
