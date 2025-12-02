import requests

def fetch_kegg_pathways():
    """Fetch all KEGG pathway IDs and names from the KEGG API."""
    url = "http://rest.kegg.jp/list/pathway"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        pathways = []
        for line in response.text.splitlines():
            pathway_id, pathway_name = line.split("\t")
            # Clean up pathway_id (remove 'path:') but keep full pathway_name
            pathway_id = pathway_id.replace("path:", "")
            pathways.append((pathway_id, pathway_name))
        return pathways
    except requests.RequestException as e:
        print(f"Error fetching KEGG pathways: {e}")
        return []

def save_pathways_to_file(pathways, filename="kegg_pathways.txt"):
    """Save pathway IDs and names to a text file."""
    with open(filename, "w", encoding="utf-8") as f:
        f.write("Pathway_ID\tPathway_Name\n")
        for pid, pname in pathways:
            f.write(f"{pid}\t{pname}\n")
    print(f"Saved {len(pathways)} pathways to {filename}")

if __name__ == "__main__":
    pathways = fetch_kegg_pathways()
    if pathways:
        save_pathways_to_file(pathways)