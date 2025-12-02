import requests
import os
import time
from datetime import datetime


def read_uniprot_ids(file_path):
    """Read UniProt IDs from a file containing only UniProt IDs (one per line)."""
    with open(file_path, "r") as file:
        uniprot_ids = [line.strip() for line in file if line.strip()]
    return uniprot_ids


def initialize_output_file(output_file):
    """Initialize the output file with headers if it doesn't exist."""
    if not os.path.exists(output_file):
        with open(output_file, "w") as f:
            f.write("Uniprot_ID\tKEGG_hsa\n")


def get_processed_ids(output_file):
    """Get the set of already processed UniProt IDs and their count from the output file."""
    processed = set()
    count = 0
    if os.path.exists(output_file):
        with open(output_file, "r") as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 1:
                    processed.add(parts[0])
                    count += 1
    return processed, count


def uniprot_to_kegg_batch(uniprot_ids, output_file, input_file, batch_size=100, max_retries=2):
    """
    Convert UniProt IDs to KEGG hsa IDs with:
    - Incremental saving
    - Automatic reconnection
    - Resume capability from last processed point
    - Proper handling of empty responses
    """
    initialize_output_file(output_file)
    processed_ids, processed_count = get_processed_ids(output_file)

    # Validate input file consistency by checking if processed IDs are at the start of input
    if processed_ids:
        for i, uniprot_id in enumerate(uniprot_ids[:len(processed_ids)]):
            if uniprot_id not in processed_ids:
                raise ValueError(f"Mismatch in input file at ID {uniprot_id}. Ensure the input file ({input_file}) is the same as in previous runs.")

    # Filter out already processed IDs
    uniprot_ids = [id for id in uniprot_ids if id not in processed_ids]

    total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size
    batch_count = 0
    total_hsa_mapped = 0

    # Set headers with user-agent
    headers = {
        'User-Agent': 'UniProt-to-KEGG-Converter/1.0 (contact: your-email@example.com)'
    }

    while batch_count < total_batches:
        start_idx = batch_count * batch_size
        end_idx = start_idx + batch_size
        current_batch = uniprot_ids[start_idx:end_idx]

        retry_count = 0
        success = False

        while not success and retry_count < max_retries:
            try:
                # Initialize mappings for this batch
                hsa_mappings = {}

                # Get hsa mappings
                hsa_query = "+".join([f"uniprot:{id}" for id in current_batch])
                hsa_url = f"https://rest.kegg.jp/conv/genes/{hsa_query}"
                hsa_response = requests.get(hsa_url, headers=headers, timeout=60)

                if hsa_response.status_code == 403:
                    raise requests.exceptions.RequestException("HTTP 403: Access forbidden. Check KEGG API restrictions or contact KEGG support.")
                elif hsa_response.status_code != 200:
                    raise requests.exceptions.RequestException(f"HTTP {hsa_response.status_code}")

                # Process hsa response only if not empty
                if hsa_response.text.strip():
                    for line in hsa_response.text.strip().split("\n"):
                        if line:
                            parts = line.split("\t")
                            if len(parts) == 2:
                                uniprot_id = parts[0].split(":")[1]
                                hsa_id = parts[1]
                                hsa_mappings[uniprot_id] = hsa_id

                # Save results for this batch - including IDs with no mappings
                with open(output_file, "a") as f:
                    for uniprot_id in current_batch:
                        hsa_id = hsa_mappings.get(uniprot_id, "NA")
                        f.write(f"{uniprot_id}\t{hsa_id}\n")

                # Update total mappings
                batch_hsa_mapped = len(hsa_mappings)
                total_hsa_mapped += batch_hsa_mapped

                success = True
                batch_count += 1
                time.sleep(3)  # Add delay to respect API rate limits

                # Print progress
                print(f"Processed batch {batch_count}/{total_batches} ({len(current_batch)} IDs)")
                print(f"  KEGG hsa IDs identified in this batch: {batch_hsa_mapped}")
                print(f"  Total KEGG hsa IDs identified so far: {total_hsa_mapped}")

            except (requests.exceptions.RequestException, requests.exceptions.Timeout) as e:
                retry_count += 1
                wait_time = 5  # Fixed 5-second wait time
                if retry_count < max_retries:
                    print(f"Error: {e}. Retrying ({retry_count}/{max_retries}) in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    print(f"Failed to process batch {batch_count + 1} after {max_retries} attempts. Saving unmapped IDs...")
                    with open("failed_batches.log", "a") as log:
                        log.write(f"Batch {batch_count + 1}: {','.join(current_batch)}\n")
                    with open(output_file, "a") as f:
                        for uniprot_id in current_batch:
                            f.write(f"{uniprot_id}\tNA\n")
                    batch_count += 1


def savefile(input_file, directory, suffix):
    base_name = os.path.splitext(os.path.basename(input_file))[0]
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
    # For a file containing only UniProt IDs (one per line)
    input_file = r"human_uniprot_ids.txt"
    output_file = r"human_uni_to_kegg.txt"

    uniprot_ids = read_uniprot_ids(input_file)
    print(f"Total UniProt IDs to process: {len(uniprot_ids)}")

    uniprot_to_kegg_batch(uniprot_ids, output_file, input_file, batch_size=100)

    print(f"Results saved to {output_file}")
