#!/usr/bin/env python3
"""
Download UniProt function-related annotations for the human reference proteome.

This script uses the UniProt REST API "stream" endpoint to download a TSV table
for all *reviewed* human proteome entries (proteome:UP000005640 AND reviewed:true)
with the following annotation columns:

- Function                -> cc_function
- Absorption             -> absorption
- Active site            -> ft_act_site
- Binding site           -> ft_binding
- Catalytic activity     -> cc_catalytic_activity
- Cofactor               -> cc_cofactor
- DNA binding            -> ft_dna_bind
- EC number              -> ec
- Activity regulation    -> cc_activity_regulation
- Function [CC]          -> cc_function (same as Function)
- Kinetics               -> kinetics
- Pathway                -> cc_pathway
- pH dependence          -> ph_dependence
- Redox potential        -> redox_potential
- Rhea ID                -> rhea
- Site                   -> ft_site
- Temperature dependence -> temp_dependence

You get one row per UniProtKB entry, with at least:
    accession, id, gene_names, and all fields above.

Usage:
    python download_uniprot_human_function_annotations.py  [output.tsv]

If output.tsv is omitted, the default is:
    uniprot_human_function_annotations.tsv
"""

import sys
import time
from typing import Optional

try:
    import requests  # type: ignore
except ImportError:  # pragma: no cover - fallback handles this at runtime.
    requests = None


UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"

# UniProt query for the *reviewed human reference proteome*
# See: proteome:UP000005640 AND reviewed:true
UNIPROT_QUERY = "proteome:UP000005640"

# Columns/fields to retrieve from UniProtKB.
# These are the current UniProt REST "field IDs" for the annotations of interest.
UNIPROT_FIELDS = [
    "accession",            # UniProt accession
    "id",                   # Entry name (e.g. P53_HUMAN)
    "gene_names",           # Gene names (space-separated)
    # Function-related annotations:
    "cc_function",          # Function / Function [CC]
    "absorption",           # Absorption
    "ft_act_site",          # Active site
    "ft_binding",           # Binding site
    "cc_catalytic_activity",# Catalytic activity
    "cc_cofactor",          # Cofactor
    "ft_dna_bind",          # DNA binding
    "ec",                   # EC number(s)
    "cc_activity_regulation",  # Activity regulation
    "kinetics",             # Kinetics
    "cc_pathway",           # Pathway
    "ph_dependence",        # pH dependence
    "redox_potential",      # Redox potential
    "rhea",                 # Rhea ID(s)
    "ft_site",              # Site
    "temp_dependence",      # Temperature dependence
]


def download_uniprot_human_annotations(
    output_path: str = r"C:\Users\clayt\Downloads\uniprot_human_function_annotations.tsv",
    query: str = UNIPROT_QUERY,
    fields: Optional[list[str]] = None,
) -> None:
    """
    Download UniProtKB annotations for the given query into a TSV file.

    Parameters
    ----------
    output_path : str
        Path to the output TSV file.
    query : str
        UniProtKB query string (default is reviewed human proteome).
    fields : list of str or None
        List of UniProt "field IDs" to retrieve. If None, uses UNIPROT_FIELDS.
    """
    if fields is None:
        fields = UNIPROT_FIELDS

    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
        # You *can* add "compressed": "true" if you want gzip,
        # but then you'd need to decompress when reading.
    }

    print("UniProt REST request:")
    print(f"  URL    : {UNIPROT_STREAM_URL}")
    print(f"  query  : {query}")
    print(f"  fields : {','.join(fields)}")
    print(f"  output : {output_path}")
    print("Starting download from UniProt (this may take a bit)...", flush=True)

    t0 = time.time()
    try:
        if requests is not None:
            line_count = _stream_download_with_requests(params, output_path)
        else:
            print("The 'requests' package is not installed; falling back to urllib.", flush=True)
            line_count = _stream_download_with_urllib(params, output_path)

        elapsed = time.time() - t0
        print(f"Download complete. ~{max(line_count - 1, 0)} entries written.")
        print(f"Saved to: {output_path}")
        print(f"Elapsed: {elapsed:0.1f} seconds")

    except Exception as e:  # noqa: BLE001 - want to show message to the user.
        print("ERROR: Problem while querying UniProt REST API.", file=sys.stderr)
        print(f"       {e}", file=sys.stderr)
        sys.exit(1)


def _stream_download_with_requests(params: dict[str, str], output_path: str) -> int:
    assert requests is not None, "requests must be available"
    with requests.get(UNIPROT_STREAM_URL, params=params, stream=True, timeout=60) as r:
        r.raise_for_status()
        line_count = 0

        with open(output_path, "wb") as out_f:
            for chunk in r.iter_content(chunk_size=8192):
                if not chunk:
                    continue
                out_f.write(chunk)
                line_count += chunk.count(b"\n")
    return line_count


def _stream_download_with_urllib(params: dict[str, str], output_path: str) -> int:
    from urllib import error as urllib_error
    from urllib import parse as urllib_parse
    from urllib import request as urllib_request

    query_string = urllib_parse.urlencode(params)
    url = f"{UNIPROT_STREAM_URL}?{query_string}"
    try:
        with urllib_request.urlopen(url, timeout=60) as resp, open(output_path, "wb") as out_f:
            line_count = 0
            while True:
                chunk = resp.read(8192)
                if not chunk:
                    break
                out_f.write(chunk)
                line_count += chunk.count(b"\n")
    except urllib_error.URLError as exc:
        raise RuntimeError(f"UniProt request failed: {exc}") from exc

    return line_count


def main(argv: list[str]) -> None:
    if len(argv) > 2:
        print("Usage: python download_uniprot_human_function_annotations.py [output.tsv]")
        sys.exit(1)

    if len(argv) == 2:
        out = argv[1]
    else:
        out = "uniprot_human_function_annotations.tsv"

    download_uniprot_human_annotations(output_path=out)


if __name__ == "__main__":
    main(sys.argv)
