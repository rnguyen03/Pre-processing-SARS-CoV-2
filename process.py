import os
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
import matplotlib.patches as mpatches

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

# Chaos Game Representation (CGR) function to calculate CGR matrix for k-mers
def cgr(seq, order, k):
    ln = len(seq)
    pw = 2**k
    out = [[0 for i in range(pw)] for j in range(pw)]
    
    x = 2**(k-1)
    y = 2**(k-1)

    for i in range(ln):
        x = x // 2
        y = y // 2
        if seq[i] == order[2] or seq[i] == order[3]:
            x = x + (2**(k-1))
        if seq[i] == order[0] or seq[i] == order[3]:
            y = y + (2**(k-1))
        if i >= k-1:
            out[y][x] = out[y][x] + 1
    
    return out

# Function to flatten CGR matrix into a 1D vector
def flatten_cgr_matrix(cgr_matrix):
    return np.array(cgr_matrix).flatten()

# Function to compute pairwise distance matrix for a set of sequences' CGRs
def compute_distance_matrix(cgr_matrices):
    flattened_matrices = [flatten_cgr_matrix(cgr_matrix) for cgr_matrix in cgr_matrices]
    pairwise_distances = pdist(flattened_matrices, metric='euclidean')  # Using Euclidean distance
    distance_matrix = squareform(pairwise_distances)  # Convert to symmetric matrix
    return distance_matrix

# Function to plot CGR
def plot_cgr(cgr_matrix, sequence_name, k, save_path):
    plt.figure(figsize=(6, 6))
    plt.imshow(cgr_matrix, cmap='gray', interpolation='nearest')  # Use grayscale
    plt.title(f'CGR Plot for {sequence_name} (k={k})')
    plt.colorbar(label='Frequency')
    plt.savefig(save_path)
    plt.close()

# Set base paths and directories
base_path = "./Asn1DataSet 2"

alpha_path = os.path.join(base_path, "Alpha")
beta_path = os.path.join(base_path, "Beta")
delta_path = os.path.join(base_path, "Delta")
gamma_path = os.path.join(base_path, "Gamma")

# Read the sequences from the directories
alpha_sequences = read_and_clean_fasta(alpha_path)
beta_sequences = read_and_clean_fasta(beta_path)
delta_sequences = read_and_clean_fasta(delta_path)
gamma_sequences = read_and_clean_fasta(gamma_path)

# Count the number of sequences in each class
num_alpha = len(alpha_sequences)
num_beta = len(beta_sequences)
num_delta = len(delta_sequences)
num_gamma = len(gamma_sequences)

# Generate a list of labels for color-coding
labels = (['Alpha'] * num_alpha +
          ['Beta'] * num_beta +
          ['Delta'] * num_delta +
          ['Gamma'] * num_gamma)

# Define colors for each class
color_map = {'Alpha': 'r', 'Beta': 'g', 'Delta': 'b', 'Gamma': 'c'}
colors = [color_map[label] for label in labels]

# Parameters
k_values = [3, 7]
alphabet = "ACGT"
output_dir = "./cgr_sequence_plots"
os.makedirs(output_dir, exist_ok=True)

# Function to generate CGR plots for each sequence
def generate_and_save_cgr_plots_for_sequences(sequences, genome_name):
    for idx, seq in enumerate(sequences):
        for k in k_values:
            cgr_matrix = cgr(seq, alphabet, k)
            sequence_name = f"{genome_name}_seq_{idx + 1}"
            save_path = os.path.join(output_dir, f"{sequence_name}_k{k}.png")
            plot_cgr(cgr_matrix, sequence_name, k, save_path)

# Combine all sequences into one list for distance matrix computation
all_sequences = alpha_sequences + beta_sequences + delta_sequences + gamma_sequences

# Function to compute and save pairwise distance matrices for given k values
def compute_and_save_distance_matrices(sequences):
    for k in k_values:
        # Compute CGRs for all sequences for the current k value
        cgr_matrices = [cgr(seq, alphabet, k) for seq in sequences]
        
        # Compute the pairwise distance matrix
        distance_matrix = compute_distance_matrix(cgr_matrices)
        
        # Perform MDS with 3 components
        mds = MDS(n_components=3, dissimilarity='precomputed', random_state=42)
        mds_coords = mds.fit_transform(distance_matrix)
        
        # Plot the MDS result in 3D
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        for i in range(len(mds_coords)):
            ax.scatter(mds_coords[i, 0], mds_coords[i, 1], mds_coords[i, 2],
                       color=colors[i], label=labels[i])
        
        # Create custom legend handles to avoid duplicate labels
        handles = []
        for label in color_map:
            handles.append(mpatches.Patch(color=color_map[label], label=label))
        ax.legend(handles=handles)
        ax.set_title(f'3D MDS Plot for k={k}')
        ax.set_xlabel('Dimension 1')
        ax.set_ylabel('Dimension 2')
        ax.set_zlabel('Dimension 3')
        plt.savefig(f'mds_3d_k{k}.png')
        plt.close()
        
        # Optionally, print or save the distance matrix
        print(f"Distance Matrix for k={k}:")
        print(distance_matrix)

# Generate CGR plots for each class (Alpha, Beta, Delta, Gamma)
generate_and_save_cgr_plots_for_sequences(alpha_sequences, "Alpha")
generate_and_save_cgr_plots_for_sequences(beta_sequences, "Beta")
generate_and_save_cgr_plots_for_sequences(delta_sequences, "Delta")
generate_and_save_cgr_plots_for_sequences(gamma_sequences, "Gamma")

# Compute and display distance matrices
compute_and_save_distance_matrices(all_sequences)

print("CGR sequence plots and 3D distance matrices generated successfully.")
