#!/usr/bin/env python3
"""
Script to clean FASTA files:
- Keep only human proteins (HUMAN)
- Remove duplicates (same sequence)
"""

import os
from pathlib import Path


def parse_fasta(file_path):
    """Parse a FASTA file and return a list of tuples (header, sequence)."""
    proteins = []
    current_header = None
    current_sequence = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    proteins.append((current_header, ''.join(current_sequence)))
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header is not None:
            proteins.append((current_header, ''.join(current_sequence)))

    return proteins


def is_human_protein(header):
    """Check if the protein is human."""
    return '_HUMAN' in header or 'Homo sapiens' in header


def clean_proteins(proteins):
    """
    Clean proteins:
    - Keep only human proteins
    - Remove duplicates (same sequence)
    """
    seen_sequences = set()
    clean = []
    removed_non_human = 0
    removed_duplicates = 0

    for header, sequence in proteins:
        # Check if human
        if not is_human_protein(header):
            removed_non_human += 1
            print(f"  [NON-HUMAN] Removed: {header[:70]}...")
            continue

        # Check if duplicate
        if sequence in seen_sequences:
            removed_duplicates += 1
            print(f"  [DUPLICATE] Removed: {header[:70]}...")
            continue

        seen_sequences.add(sequence)
        clean.append((header, sequence))

    return clean, removed_non_human, removed_duplicates


def write_fasta(proteins, output_path):
    """Write proteins to a FASTA file."""
    with open(output_path, 'w') as f:
        for header, sequence in proteins:
            f.write(f"{header}\n")
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")


def process_fasta_files(fasta_dir):
    """Process all FASTA files in the directory."""
    fasta_dir = Path(fasta_dir)
    total_removed_non_human = 0
    total_removed_duplicates = 0
    total_kept = 0

    # Global set to detect duplicates across files
    global_seen_sequences = set()

    for fasta_file in sorted(fasta_dir.glob('*.fasta')):
        print(f"\nProcessing: {fasta_file.name}")

        proteins = parse_fasta(fasta_file)
        original_count = len(proteins)

        # Filter non-human first
        human_proteins = []
        removed_nh = 0
        for header, sequence in proteins:
            if not is_human_protein(header):
                removed_nh += 1
                print(f"  [NON-HUMAN] Removed: {header[:70]}...")
            else:
                human_proteins.append((header, sequence))

        # Filter global duplicates
        clean = []
        removed_dup = 0
        for header, sequence in human_proteins:
            if sequence in global_seen_sequences:
                removed_dup += 1
                print(f"  [GLOBAL DUPLICATE] Removed: {header[:60]}...")
            else:
                global_seen_sequences.add(sequence)
                clean.append((header, sequence))

        total_removed_non_human += removed_nh
        total_removed_duplicates += removed_dup
        total_kept += len(clean)

        if removed_nh > 0 or removed_dup > 0:
            if len(clean) > 0:
                write_fasta(clean, fasta_file)
                print(f"  Result: {original_count} -> {len(clean)} proteins")
            else:
                # Delete file if no proteins remain
                os.remove(fasta_file)
                print(f"  File deleted (no proteins remaining)")
        else:
            print(f"  OK: {len(proteins)} human protein(s)")

    return total_kept, total_removed_non_human, total_removed_duplicates


if __name__ == "__main__":
    fasta_directory = Path(__file__).parent / "Fasta"

    if not fasta_directory.exists():
        print(f"Error: Directory {fasta_directory} does not exist")
        exit(1)

    print("=" * 70)
    print("FASTA File Cleanup")
    print("- Keep only HUMAN proteins")
    print("- Remove duplicates")
    print("=" * 70)

    kept, non_human, duplicates = process_fasta_files(fasta_directory)

    print("\n" + "=" * 70)
    print("SUMMARY:")
    print(f"  Proteins kept:         {kept}")
    print(f"  Non-human removed:     {non_human}")
    print(f"  Duplicates removed:    {duplicates}")
    print("=" * 70)
