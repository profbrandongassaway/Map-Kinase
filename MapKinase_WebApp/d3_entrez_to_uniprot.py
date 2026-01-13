import requests

# ✅ Put your Entrez Gene ID here
ENTREZ_GENE_ID = "2475"  # example: MTOR is 2475

def entrez_to_uniprot(entrez_id: str, organism_id: int = 9606) -> str | None:
    """
    Convert an Entrez Gene ID -> UniProt primary accession (best hit),
    by querying UniProt's REST Search endpoint.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        # Use the UniProt cross-reference query for GeneID (Entrez)
        "query": f'(xref:GeneID-{entrez_id}) AND (organism_id:{organism_id})',
        "format": "tsv",
        "fields": "accession",
        "size": 1,  # just the top hit
    }
    headers = {
        "User-Agent": "EntrezToUniProt/1.0 (requests)"
    }

    r = requests.get(url, params=params, headers=headers, timeout=30)
    r.raise_for_status()

    lines = [ln.strip() for ln in r.text.splitlines() if ln.strip()]
    # TSV format: header line + 1 data line (if found)
    if len(lines) < 2:
        return None

    # data line should be just the accession (since fields=accession)
    return lines[1].split("\t")[0].strip()

def ensembl_to_uniprot(ensembl_id: str, organism_id: int = 9606) -> str | None:
    """
    Convert an Ensembl Gene ID -> UniProt primary accession (best hit),
    by querying UniProt's REST Search endpoint.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        # Use the UniProt cross-reference query for Ensembl
        "query": f'(xref:Ensembl-{ensembl_id}) AND (organism_id:{organism_id})',
        "format": "tsv",
        "fields": "accession",
        "size": 1,
    }
    headers = {"User-Agent": "EnsemblToUniProt/1.0 (requests)"}
    r = requests.get(url, params=params, headers=headers, timeout=30)
    r.raise_for_status()
    lines = [ln.strip() for ln in r.text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return None
    return lines[1].split("\t")[0].strip()

if __name__ == "__main__":
    uniprot_id = entrez_to_uniprot(ENTREZ_GENE_ID)

    if uniprot_id:
        print(uniprot_id)
    else:
        print(f"No UniProt ID found for Entrez Gene ID {ENTREZ_GENE_ID} (human).")
