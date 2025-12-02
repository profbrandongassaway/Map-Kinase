import requests

def fetch_kegg_organisms():
    """Fetch all KEGG organism codes and names from the KEGG API."""
    url = "http://rest.kegg.jp/list/organism"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        lines = response.text.splitlines()
        print(f"Received {len(lines)} lines from API")
        organisms = []
        for line in lines:
            try:
                # Format: T number, code, name, taxonomy
                parts = line.split("\t")
                if len(parts) < 3:
                    print(f"Skipping malformed line: {line}")
                    continue
                org_code = parts[1]  # e.g., 'hsa'
                org_name = parts[2]  # e.g., 'Homo sapiens (human)'
                organisms.append((org_code, org_name))
            except Exception as e:
                print(f"Error processing line '{line}': {e}")
                continue
        print(f"Processed {len(organisms)} organism codes")
        return organisms
    except requests.RequestException as e:
        print(f"Error fetching KEGG organisms: {e}")
        return []

def save_organisms_to_file(organisms, filename="kegg_organism_codes.txt"):
    """Save organism codes and names to a text file."""
    try:
        with open(filename, "w", encoding="utf-8") as f:
            f.write("Organism_Code\tOrganism_Name\n")
            for code, name in organisms:
                f.write(f"{code}\t{name}\n")
        print(f"Saved {len(organisms)} organism codes to {filename}")
    except Exception as e:
        print(f"Error saving file: {e}")

if __name__ == "__main__":
    organisms = fetch_kegg_organisms()
    if organisms:
        save_organisms_to_file(organisms)