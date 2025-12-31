#!/usr/bin/env python3
"""
Script to build a phylogenetic tree of PDE proteins
using the Neighbor-Joining algorithm.
"""

import os
from pathlib import Path
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt


def parse_fasta_files(fasta_dir):
    """Parse all FASTA files and return a dict {name: sequence}."""
    sequences = {}
    fasta_dir = Path(fasta_dir)

    for fasta_file in fasta_dir.glob('*.fasta'):
        with open(fasta_file, 'r') as f:
            current_name = None
            current_seq = []

            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_name and current_seq:
                        sequences[current_name] = ''.join(current_seq)
                    # Extract short name (e.g., PDE4A_HUMAN)
                    parts = line.split('|')
                    if len(parts) >= 3:
                        name = parts[2].split()[0]
                    else:
                        name = line[1:].split()[0]
                    current_name = name
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_name and current_seq:
                sequences[current_name] = ''.join(current_seq)

    return sequences


def calculate_kmer_similarity(seq1, seq2, k=3):
    """Calculate similarity between two sequences based on k-mers."""
    def get_kmers(seq, k):
        return set(seq[i:i+k] for i in range(len(seq) - k + 1))

    kmers1 = get_kmers(seq1, k)
    kmers2 = get_kmers(seq2, k)

    if not kmers1 or not kmers2:
        return 0.0

    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)

    return intersection / union if union > 0 else 0.0


def calculate_distance_matrix(sequences):
    """Calculate the distance matrix between all sequences."""
    names = list(sequences.keys())
    n = len(names)
    distances = np.zeros((n, n))

    print(f"Computing distance matrix ({n}x{n})...")

    for i in range(n):
        for j in range(i + 1, n):
            similarity = calculate_kmer_similarity(
                sequences[names[i]],
                sequences[names[j]],
                k=3
            )
            # Convert similarity to distance
            distance = 1.0 - similarity
            distances[i, j] = distance
            distances[j, i] = distance

    return names, distances


def neighbor_joining(names, distance_matrix):
    """
    Implementation of the Neighbor-Joining algorithm.
    Returns a tree as a nested structure.
    """
    n = len(names)
    if n < 2:
        return names[0] if names else None

    # Copy to avoid modifying the original
    D = distance_matrix.copy()
    nodes = list(names)

    while len(nodes) > 2:
        n = len(nodes)

        # Calculate Q-matrix
        Q = np.zeros((n, n))
        row_sums = D.sum(axis=1)

        for i in range(n):
            for j in range(i + 1, n):
                Q[i, j] = (n - 2) * D[i, j] - row_sums[i] - row_sums[j]
                Q[j, i] = Q[i, j]

        # Find pair with minimum Q
        np.fill_diagonal(Q, np.inf)
        min_idx = np.unravel_index(Q.argmin(), Q.shape)
        i, j = min(min_idx), max(min_idx)

        # Calculate distances to new nodes
        if n > 2:
            delta = (row_sums[i] - row_sums[j]) / (n - 2)
        else:
            delta = 0

        dist_i = (D[i, j] + delta) / 2
        dist_j = (D[i, j] - delta) / 2

        # Create new node
        new_node = (nodes[i], nodes[j], dist_i, dist_j)

        # Calculate distances from new node to others
        new_distances = []
        for k in range(n):
            if k != i and k != j:
                new_dist = (D[i, k] + D[j, k] - D[i, j]) / 2
                new_distances.append(new_dist)

        # Update matrix and nodes
        remaining = [k for k in range(n) if k != i and k != j]
        new_D = np.zeros((len(remaining) + 1, len(remaining) + 1))

        # Copy existing distances
        for new_i, old_i in enumerate(remaining):
            for new_j, old_j in enumerate(remaining):
                new_D[new_i, new_j] = D[old_i, old_j]

        # Add distances for new node
        for new_i, dist in enumerate(new_distances):
            new_D[new_i, -1] = dist
            new_D[-1, new_i] = dist

        D = new_D
        nodes = [nodes[k] for k in remaining] + [new_node]

    # Connect the last two nodes
    if len(nodes) == 2:
        return (nodes[0], nodes[1], D[0, 1] / 2, D[0, 1] / 2)
    return nodes[0]


def tree_to_newick(tree, branch_length=None):
    """Convert tree to Newick format."""
    if isinstance(tree, str):
        if branch_length is not None:
            return f"{tree}:{branch_length:.4f}"
        return tree

    if isinstance(tree, tuple) and len(tree) == 4:
        left, right, dist_left, dist_right = tree
        left_str = tree_to_newick(left, dist_left)
        right_str = tree_to_newick(right, dist_right)
        result = f"({left_str},{right_str})"
        if branch_length is not None:
            result += f":{branch_length:.4f}"
        return result

    return str(tree)


def plot_dendrogram(names, distance_matrix, output_file):
    """Create a dendrogram with scipy and save it."""
    from scipy.cluster.hierarchy import set_link_color_palette, fcluster

    # Convert to condensed format for scipy
    condensed = squareform(distance_matrix)

    # Create linkage (average method = UPGMA-like)
    Z = linkage(condensed, method='average')

    # Calculate optimal number of colors based on PDE families
    # Extract family names (e.g., PDE4 from PDE4A_HUMAN)
    families = set()
    for name in names:
        # Extract PDE family number (e.g., "PDE4" from "PDE4A_HUMAN")
        import re
        match = re.match(r'(PDE\d+)', name)
        if match:
            families.add(match.group(1))

    n_families = len(families)
    n_sequences = len(names)

    # Optimal colors: number of families (11 for PDE1-11)
    n_colors = max(n_families, 5)  # At least 5 colors
    print(f"  Detected {n_families} PDE families, using {n_colors} colors")

    # Define a vibrant color palette with enough distinct colors
    color_palette = [
        '#e6194B',  # Red
        '#3cb44b',  # Green
        '#ffe119',  # Yellow
        '#4363d8',  # Blue
        '#f58231',  # Orange
        '#911eb4',  # Purple
        '#42d4f4',  # Cyan
        '#f032e6',  # Magenta
        '#bfef45',  # Lime
        '#fabed4',  # Pink
        '#469990',  # Teal
        '#dcbeff',  # Lavender
        '#9A6324',  # Brown
        '#fffac8',  # Beige
        '#800000',  # Maroon
        '#aaffc3',  # Mint
    ]

    # Set custom color palette
    set_link_color_palette(color_palette[:n_colors])

    # Calculate color threshold to get approximately n_families clusters
    # Sort merge distances and find threshold that gives desired clusters
    heights = Z[:, 2]
    sorted_heights = np.sort(heights)[::-1]

    # Find threshold that creates approximately n_families clusters
    if n_families > 1 and len(sorted_heights) >= n_families:
        color_threshold = sorted_heights[n_families - 2] * 0.99
    else:
        color_threshold = 0.5 * max(heights)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 9))
    plt.title('Phylogenetic Tree of PDE Proteins\n(Neighbor-Joining based on k-mer similarity)',
              fontsize=16, fontweight='bold', pad=20)

    dendrogram(
        Z,
        labels=names,
        leaf_rotation=45,
        leaf_font_size=11,
        color_threshold=color_threshold,
        above_threshold_color='#808080',  # Gray for links above threshold
        ax=ax
    )

    plt.xlabel('Proteins', fontsize=13, fontweight='bold')
    plt.ylabel('Distance (1 - k-mer similarity)', fontsize=13, fontweight='bold')

    # Add grid for better readability
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()

    # Reset color palette to default
    set_link_color_palette(None)

    print(f"Dendrogram saved: {output_file}")


def print_distance_matrix(names, distances):
    """Display the distance matrix in a readable format."""
    print("\nDistance matrix:")
    print("-" * 60)

    # Header
    header = "          " + "  ".join(f"{n[:8]:>8}" for n in names)
    print(header)

    for i, name in enumerate(names):
        row = f"{name[:8]:>8}  " + "  ".join(f"{distances[i,j]:8.3f}" for j in range(len(names)))
        print(row)


def main():
    script_dir = Path(__file__).parent
    fasta_dir = script_dir / "Fasta"
    output_file = script_dir / "phylogenetic_tree.png"
    newick_file = script_dir / "phylogenetic_tree.nwk"

    print("=" * 60)
    print("Phylogenetic Tree Construction")
    print("=" * 60)

    # Load sequences
    print(f"\nLoading sequences from: {fasta_dir}")
    sequences = parse_fasta_files(fasta_dir)
    print(f"Sequences loaded: {len(sequences)}")

    for name, seq in sequences.items():
        print(f"  - {name}: {len(seq)} amino acids")

    if len(sequences) < 2:
        print("Error: At least 2 sequences are required to build a tree")
        return

    # Calculate distance matrix
    names, distances = calculate_distance_matrix(sequences)
    print_distance_matrix(names, distances)

    # Build tree with Neighbor-Joining
    print("\nBuilding tree (Neighbor-Joining)...")
    tree = neighbor_joining(names.copy(), distances.copy())

    # Save in Newick format
    newick = tree_to_newick(tree) + ";"
    with open(newick_file, 'w') as f:
        f.write(newick)
    print(f"Newick tree saved: {newick_file}")
    print(f"Newick: {newick[:100]}...")

    # Create dendrogram
    print("\nCreating dendrogram...")
    plot_dendrogram(names, distances, output_file)

    print("\n" + "=" * 60)
    print("Done!")
    print(f"- Image: {output_file}")
    print(f"- Newick: {newick_file}")


if __name__ == "__main__":
    main()
