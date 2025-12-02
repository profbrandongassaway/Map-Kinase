import requests
import os
from datetime import datetime

def read_uniprot_ids(file_path, column_name):
    """Read UniProt IDs from a specific column in a tab-delimited file."""
    uniprot_ids = []
    with open(file_path, "r") as file:
        headers = file.readline().strip().split("\t")
        try:
            col_index = headers.index(column_name)
        except ValueError:
            raise ValueError(f"Column '{column_name}' not found in the file.")

        for line in file:
            parts = line.strip().split("\t")
            if len(parts) < len(headers):
                parts.extend([""] * (len(headers) - len(parts)))
            if len(parts) > col_index:
                uniprot_ids.append(parts[col_index])
    return uniprot_ids, headers

def uniprot_to_kegg_batch(uniprot_ids, batch_size=100, progress_callback=None):
    """Convert a list of UniProt IDs to KEGG IDs (hsa and ko) using the KEGG API."""
    hsa_mappings = {}
    ko_mappings = {}
    total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size  # Ceiling division
    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i:i + batch_size]
        hsa_query = "+".join([f"uniprot:{id}" for id in batch])
        hsa_url = f"https://rest.kegg.jp/conv/genes/{hsa_query}"
        hsa_response = requests.get(hsa_url)
        if hsa_response.status_code == 200:
            for line in hsa_response.text.strip().split("\n"):
                parts = line.split("\t")
                if len(parts) == 2:
                    uniprot_id = parts[0].split(":")[1]
                    hsa_id = parts[1]
                    hsa_mappings[uniprot_id] = hsa_id

        if hsa_mappings:
            hsa_ids = list(hsa_mappings.values())
            ko_query = "+".join(hsa_ids)
            ko_url = f"https://rest.kegg.jp/link/ko/{ko_query}"
            ko_response = requests.get(ko_url)
            if ko_response.status_code == 200:
                for line in ko_response.text.strip().split("\n"):
                    parts = line.split("\t")
                    if len(parts) == 2:
                        hsa_id = parts[0]
                        ko_id = parts[1]
                        for uniprot_id, mapped_hsa_id in hsa_mappings.items():
                            if mapped_hsa_id == hsa_id:
                                ko_mappings[uniprot_id] = ko_id

        # Update progress if callback is provided
        if progress_callback:
            current_batch = i // batch_size + 1
            progress_callback(current_batch, total_batches)

    return hsa_mappings, ko_mappings

def save_results(file_path, uniprot_ids, hsa_mappings, ko_mappings, headers, output_file, column_name):
    """Save the results to a new file with new columns for hsa and ko IDs."""
    with open(file_path, "r") as infile, open(output_file, "w") as outfile:
        outfile.write("\t".join(headers) + "\tKEGG_Gene_ID\tKEGG_ko\n")
        next(infile)
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) < len(headers):
                parts.extend([""] * (len(headers) - len(parts)))
            uniprot_id = parts[headers.index(column_name)]
            hsa_id = hsa_mappings.get(uniprot_id, "NA")
            ko_id = ko_mappings.get(uniprot_id, "NA")
            parts.extend([hsa_id, ko_id])
            outfile.write("\t".join(parts) + "\n")

def savefile(input_file, directory, suffix):
    base_name = os.path.splitext(input_file)[0]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    phosmap_dir = os.path.dirname(script_dir)
    temp_folder = os.path.join(phosmap_dir, directory)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_name = f"{base_name}_{suffix}_{timestamp}.txt"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    output_path = os.path.join(temp_folder, file_name)
    return output_path

if __name__ == "__main__":
    input_file = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\ProtMap.txt"
    column_name = "Uniprot_ID"
    output_file = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\ProtMaphsaanno.txt"

    uniprot_ids, headers = read_uniprot_ids(input_file, column_name)
    hsa_mappings, ko_mappings = uniprot_to_kegg_batch(uniprot_ids)
    save_results(input_file, uniprot_ids, hsa_mappings, ko_mappings, headers, output_file, column_name)

    print(f"Results saved to {output_file}")