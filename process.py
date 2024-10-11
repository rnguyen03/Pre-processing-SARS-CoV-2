import os
from Bio import SeqIO
import re

# Function to clean and preprocess sequences
def clean_sequence(sequence):
    sequence = sequence.upper()  # Convert to uppercase
    sequence = re.sub(r'[^ACGT]', '', sequence)  # Keep only A, C, G, T
    return sequence

# Function to read all fasta files in a given directory and clean sequences
def read_and_clean_fasta(directory):
    print(f"Checking directory: {directory}")
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return []
    
    cleaned_sequences = []
    # Iterate over each fasta file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            file_path = os.path.join(directory, filename)
            for record in SeqIO.parse(file_path, "fasta"):
                cleaned_seq = clean_sequence(str(record.seq))
                cleaned_sequences.append(cleaned_seq)
    return cleaned_sequences

# Correct the base path to point to your current working directory structure
base_path = "./Asn1DataSet 2"

# Provide correct paths to the directories
alpha_path = os.path.join(base_path, "Alpha")
beta_path = os.path.join(base_path, "Beta")
delta_path = os.path.join(base_path, "Delta")
gamma_path = os.path.join(base_path, "Gamma")

# Print the paths to check if they are correct
print(f"Alpha path: {os.path.abspath(alpha_path)}")
print(f"Beta path: {os.path.abspath(beta_path)}")
print(f"Delta path: {os.path.abspath(delta_path)}")
print(f"Gamma path: {os.path.abspath(gamma_path)}")

# Try reading the fasta files
alpha_sequences = read_and_clean_fasta(alpha_path)
beta_sequences = read_and_clean_fasta(beta_path)
delta_sequences = read_and_clean_fasta(delta_path)
gamma_sequences = read_and_clean_fasta(gamma_path)

# Now you have all sequences cleaned and stored in their respective lists
print(f"Alpha sequences: {len(alpha_sequences)}")
print(f"Beta sequences: {len(beta_sequences)}")
print(f"Delta sequences: {len(delta_sequences)}")
print(f"Gamma sequences: {len(gamma_sequences)}")
