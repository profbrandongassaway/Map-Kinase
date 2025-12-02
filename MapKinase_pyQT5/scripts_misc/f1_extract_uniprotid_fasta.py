# extract_uniprot_ids.py

def extract_uniprot_ids(fasta_path, output_path):
    with open(fasta_path, 'r') as infile, open(output_path, 'w') as outfile:
        outfile.write("UniProtID\n")  # header
        for line in infile:
            if line.startswith('>'):
                # Format: >tr|A0A087WZ06|A0A087WZ06_HUMAN ...
                parts = line.split('|')
                if len(parts) > 1:
                    uniprot_id = parts[1]
                    outfile.write(f"{uniprot_id}\n")

# Example usage
fasta_file = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Uniprot_files\uniprotkb_DROME_proteome_UP000000803_2025_04_26.fasta.txt"
output_file = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\Uniprot_files\DROME_uniprot_ids.txt"
extract_uniprot_ids(fasta_file, output_file)
