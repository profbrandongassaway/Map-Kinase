import gzip
import os
import csv
import re
from pathlib import Path


def read_regulatory_data(psp_folder_path, compressed_file="Regulatory_sites.gz"):
    """
    Reads tab-separated data from a Regulatory_sites.gz file inside a PSP folder and stores it in a dictionary.

    Args:
        psp_folder_path (str): Path to the PSP folder (e.g., 'PSP')
        compressed_file (str): Name of the compressed file (default: 'Regulatory_sites.gz')

    Returns:
        dict: Dictionary with keys as 'ACC_ID:MOD_RSD_NUM' (e.g., 'P41181:270') and values as dictionaries of all fields
    """
    # Construct the full path to the compressed file
    full_path = os.path.join(psp_folder_path, compressed_file)

    if not os.path.exists(full_path):
        raise FileNotFoundError(f"The compressed file {full_path} does not exist")

    # Dictionary to store the data
    data_dict = {}

    try:
        with gzip.open(full_path, 'rt', encoding='utf-8') as gz_ref:
            lines = gz_ref.readlines()

            if len(lines) < 3:
                raise ValueError("The file is too short; expected at least 3 lines (metadata, license, and header)")

            # Find the header by looking for a line that contains 'GENE' and 'PROTEIN'
            header_line = None
            header_index = None
            for i, line in enumerate(lines):
                if 'GENE' in line and 'PROTEIN' in line:
                    header_line = line
                    header_index = i
                    break

            if header_line is None:
                raise ValueError("Could not find the header line (expected a line containing 'GENE' and 'PROTEIN')")

            # Use the identified header line
            header = header_line.strip().split('\t')

            # Process each data row (skip lines up to and including the header)
            for i, line in enumerate(lines[header_index + 1:], start=header_index + 2):
                fields = line.strip().split('\t')

                if len(fields) != len(header):
                    if len(fields) < len(header):
                        fields.extend([''] * (len(header) - len(fields)))
                    else:
                        fields = fields[:len(header)]
                    print(
                        f"Warning: Adjusted line {i} with {len(fields)} fields to match header ({len(header)} fields): {line.strip()}")

                row_dict = dict(zip(header, fields))

                # Extract the numerical part of MOD_RSD (e.g., 'K270-ub' -> '270')
                mod_rsd = row_dict.get('MOD_RSD', '')
                mod_rsd_num = ''
                if mod_rsd:
                    # Use regex to extract the number between the letter and the modification type
                    match = re.match(r'^\D+(\d+)-.*$', mod_rsd)
                    if match:
                        mod_rsd_num = match.group(1)
                    else:
                        print(f"Warning: Could not parse MOD_RSD '{mod_rsd}' in line {i}; skipping")
                        continue

                if not row_dict.get('ACC_ID') or not mod_rsd_num:
                    print(f"Warning: Skipping line {i} due to missing ACC_ID or MOD_RSD number: {line.strip()}")
                    continue

                # Create the key using ACC_ID and the numerical part of MOD_RSD
                key = f"{row_dict['ACC_ID']}:{mod_rsd_num}"
                data_dict[key] = row_dict

        return data_dict

    except gzip.BadGzipFile:
        raise ValueError(f"The file {full_path} is not a valid GZIP file.")
    except Exception as e:
        raise Exception(f"An error occurred while processing the file: {str(e)}")


def annotate_dataset(input_file, output_file, uniprot_col, site_col, regulatory_data):
    """
    Annotates a tab-delimited dataset with regulatory site information.

    Args:
        input_file (str): Path to the input tab-delimited file
        output_file (str): Path to the output tab-delimited file
        uniprot_col (str): Column header containing UniProt IDs
        site_col (str): Column header containing site IDs (phosphosites)
        regulatory_data (dict): Dictionary of regulatory data with keys as 'ACC_ID:MOD_RSD_NUM'
    """
    # Read the input dataset
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames
        if not headers:
            raise ValueError("Input file is empty or has no headers")

        if uniprot_col not in headers:
            raise ValueError(f"UniProt column '{uniprot_col}' not found in input file headers: {headers}")
        if site_col not in headers:
            raise ValueError(f"Site column '{site_col}' not found in input file headers: {headers}")

        # Add new columns to headers
        new_headers = headers + ['Regulatory site', 'Regulatory site function', 'Regulatory site process']
        rows = list(reader)

    # Annotate each row
    for row in rows:
        uniprot_id = row.get(uniprot_col, '')
        site_id = row.get(site_col, '')

        # Default values for new columns
        row['Regulatory site'] = ''
        row['Regulatory site function'] = ''
        row['Regulatory site process'] = ''

        # Skip if either UniProt ID or site ID is missing
        if not uniprot_id or not site_id:
            continue

        # Create the key to match against regulatory_data
        key = f"{uniprot_id}:{site_id}"

        # Check for a match
        if key in regulatory_data:
            regulatory_entry = regulatory_data[key]
            row['Regulatory site'] = '+'
            row['Regulatory site function'] = regulatory_entry.get('ON_FUNCTION', 'na') or 'na'
            row['Regulatory site process'] = regulatory_entry.get('ON_PROCESS', 'na') or 'na'

    # Write the updated dataset to the output file
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=new_headers, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(f"Annotated dataset written to {output_file}")


# Example usage
if __name__ == "__main__":
    # Paths and parameters
    psp_path = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Graduate_Documents\Perseus\AA Perseus Backend\PSP"  # Path to the PSP folder
    input_file = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\Phosmap\BCp2_PhosMaphsanno.txt"  # Path to your tab-delimited dataset
    output_file = "your_dataset_annotated.txt"  # Path to the output annotated dataset
    uniprot_col = "T: Uniprot_ID"  # Column header for UniProt IDs in your dataset
    site_col = "T: Site Position"  # Column header for site IDs in your dataset

    try:
        # Step 1: Load the regulatory data
        regulatory_data = read_regulatory_data(psp_path)
        print(f"Loaded {len(regulatory_data)} regulatory records")

        # Step 2: Annotate your dataset
        annotate_dataset(input_file, output_file, uniprot_col, site_col, regulatory_data)

    except Exception as e:
        print(f"Error: {e}")