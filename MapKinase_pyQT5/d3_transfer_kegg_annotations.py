import pandas as pd
import os
from datetime import datetime

def transfer_kegg_annotations(prot_file, phos_file, output_file, prot_uniprot_col, phos_uniprot_col, kegg_annotation_col):
    """
    Transfer KEGG annotations from a proteomic dataset to a phosphoproteomic dataset based on UniProt IDs.

    Args:
        prot_file (str): Path to the proteomic dataset file.
        phos_file (str): Path to the phosphoproteomic dataset file.
        output_file (str): Path to save the updated phosphoproteomic dataset.
        prot_uniprot_col (str): Name of the UniProt ID column in the proteomic dataset.
        phos_uniprot_col (str): Name of the UniProt ID column in the phosphoproteomic dataset.
        kegg_annotation_col (str): Name of the KEGG annotation column in the proteomic dataset.
    """
    # Load the datasets
    prot_data = pd.read_csv(prot_file, delimiter='\t')
    phos_data = pd.read_csv(phos_file, delimiter='\t')

    # Ensure the UniProt ID columns are of the same type (string)
    prot_data[prot_uniprot_col] = prot_data[prot_uniprot_col].astype(str)
    phos_data[phos_uniprot_col] = phos_data[phos_uniprot_col].astype(str)

    # Create a dictionary to map UniProt IDs to KEGG annotations
    kegg_map = prot_data.set_index(prot_uniprot_col)[kegg_annotation_col].to_dict()

    # Add a new column to the phosphoproteomic dataset for KEGG annotations
    phos_data[kegg_annotation_col] = phos_data[phos_uniprot_col].map(kegg_map)

    # Save the updated phosphoproteomic dataset
    phos_data.to_csv(output_file, sep='\t', index=False)
    print(f"Updated phosphoproteomic dataset saved to: {output_file}")

if __name__ == "__main__":
    # Define variables only when running directly
    prot_file = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\BCp2_data_03272025\Phosmap\BCp2_ProtMaphsaanno.txt"
    phos_file = r"C:\Users\clayt\OneDrive - Brigham Young University\pycharm\ttk_projects\Phosmap\output\testing_file_001\pY Data.txt"
    output_file = r"C:\Users\clayt\OneDrive - Brigham Young University\pycharm\ttk_projects\Phosmap\output\testing_file_001\pY Data2.txt"
    prot_uniprot_col = "Uniprot_ID"
    phos_uniprot_col = "T: Uniprot_ID"
    kegg_annotation_col = "KEGG_Gene_ID"

    transfer_kegg_annotations(prot_file, phos_file, output_file, prot_uniprot_col, phos_uniprot_col, kegg_annotation_col)