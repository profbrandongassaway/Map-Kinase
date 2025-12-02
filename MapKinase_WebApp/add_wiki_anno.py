import requests
import pandas as pd
import os
import time
from datetime import datetime


def initialize_output_file(output_file):
    """Initialize the output file with headers if it doesn't exist."""
    if not os.path.exists(output_file):
        with open(output_file, "w") as f:
            f.write("Uniprot_ID\tKEGG_Gene_ID\tEntrez_Gene_ID\tEnsembl_ID\n")


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


def get_annotations(uniprot_id):
    """Query UniProt API for Entrez Gene and Ensembl IDs."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    headers = {'User-Agent': 'UniProt-Annotation-Converter/1.0 (contact: your-email@example.com)'}
    try:
        response = requests.get(url, headers=headers, timeout=60)
        response.raise_for_status()
        data = response.json()

        entrez_id = "NA"
        ensembl_id = "NA"

        if "uniProtKBCrossReferences" in data:
            for ref in data["uniProtKBCrossReferences"]:
                if ref["database"] == "GeneID":
                    entrez_id = ref["id"]
                elif ref["database"] == "Ensembl":
                    ensembl_id = ref["id"]

        return entrez_id, ensembl_id
    except requests.exceptions.RequestException as e:
        print(f"Error querying UniProt for {uniprot_id}: {e}")
        return "NA", "NA"


def annotate_uniprot_ids(input_file, output_file, batch_size=500, max_retries=2):
    """Process UniProt IDs in batches to add Entrez Gene and Ensembl annotations."""
    # Read input file
    try:
        df = pd.read_csv(input_file, sep="\t")
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        exit()

    # Verify column names
    expected_columns = ["Uniprot_ID", "KEGG_Gene_ID"]
    if not all(col in df.columns for col in expected_columns):
        print(f"Error: Input file must contain columns: {expected_columns}. Found: {df.columns.tolist()}")
        exit()

    # Initialize output file
    initialize_output_file(output_file)
    processed_ids, processed_count = get_processed_ids(output_file)

    # Validate input file consistency
    if processed_ids:
        for i, uniprot_id in enumerate(df["Uniprot_ID"].iloc[:len(processed_ids)]):
            if uniprot_id not in processed_ids:
                raise ValueError(
                    f"Mismatch in input file at ID {uniprot_id}. Ensure the input file ({input_file}) is the same as in previous runs.")

    # Filter out already processed IDs
    df_to_process = df[~df["Uniprot_ID"].isin(processed_ids)]
    uniprot_ids = df_to_process["Uniprot_ID"].tolist()

    total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size
    batch_count = 0
    total_annotated = 0

    while batch_count < total_batches:
        start_idx = batch_count * batch_size
        end_idx = start_idx + batch_size
        current_batch = uniprot_ids[start_idx:end_idx]
        batch_df = df_to_process.iloc[start_idx:end_idx]

        retry_count = 0
        success = False

        while not success and retry_count < max_retries:
            try:
                # Initialize mappings for this batch
                annotations = {uniprot_id: {"entrez": "NA", "ensembl": "NA"} for uniprot_id in current_batch}

                # Query annotations for the batch
                for uniprot_id in current_batch:
                    entrez_id, ensembl_id = get_annotations(uniprot_id)
                    annotations[uniprot_id] = {"entrez": entrez_id, "ensembl": ensembl_id}
                    if entrez_id != "NA" or ensembl_id != "NA":
                        total_annotated += 1

                # Save results for this batch
                with open(output_file, "a") as f:
                    for index, row in batch_df.iterrows():
                        uniprot_id = row["Uniprot_ID"]
                        kegg_id = row["KEGG_Gene_ID"]
                        entrez_id = annotations[uniprot_id]["entrez"]
                        ensembl_id = annotations[uniprot_id]["ensembl"]
                        f.write(f"{uniprot_id}\t{kegg_id}\t{entrez_id}\t{ensembl_id}\n")

                success = True
                batch_count += 1
                time.sleep(3)  # Delay to respect API rate limits

                # Print progress
                print(f"Processed batch {batch_count}/{total_batches} ({len(current_batch)} IDs)")
                print(
                    f"  Annotations added in this batch: {sum(1 for uniprot_id in current_batch if annotations[uniprot_id]['entrez'] != 'NA' or annotations[uniprot_id]['ensembl'] != 'NA')}")
                print(f"  Total annotations added so far: {total_annotated}")

            except (requests.exceptions.RequestException, requests.exceptions.Timeout) as e:
                retry_count += 1
                wait_time = 5
                if retry_count < max_retries:
                    print(f"Error: {e}. Retrying ({retry_count}/{max_retries}) in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    print(
                        f"Failed to process batch {batch_count + 1} after {max_retries} attempts. Saving unmapped IDs...")
                    with open("failed_batches.log", "a") as log:
                        log.write(f"Batch {batch_count + 1}: {','.join(current_batch)}\n")
                    with open(output_file, "a") as f:
                        for index, row in batch_df.iterrows():
                            uniprot_id = row["Uniprot_ID"]
                            kegg_id = row["KEGG_Gene_ID"]
                            f.write(f"{uniprot_id}\t{kegg_id}\tNA\tNA\n")
                    batch_count += 1


def savefile(input_file, directory, suffix):
    """Generate output file path with timestamp."""
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
    input_file = r"/Phosmap/scripts/Human_KEGG_Conversion.txt"
    output_file = savefile(input_file, "output", "annotated")
    annotate_uniprot_ids(input_file, output_file, batch_size=100)
    print(f"Results saved to {output_file}")