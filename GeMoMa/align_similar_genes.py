#!/usr/bin/env python3

# ======================
# Import required modules
# ======================
import argparse  # For parsing command-line arguments
import os        # For file path operations
import subprocess  # For running external alignment programs
from Bio import SeqIO  # For reading and writing FASTA files

# ======================
# Define the main function
# ======================
def main():
    # ----------------------
    # Set up argument parser
    # ----------------------
    parser = argparse.ArgumentParser(description="Align query-target sequence pairs from two FASTA files using MAFFT or MUSCLE.")
    parser.add_argument("--table", required=True, help="Path to TSV/CSV table with at least two columns: query_id and target_id.")
    parser.add_argument("--query_fasta", required=True, help="FASTA file containing query sequences.")
    parser.add_argument("--target_fasta", required=True, help="FASTA file containing target sequences.")
    parser.add_argument("--outdir", required=True, help="Output directory where alignments will be saved.")
    parser.add_argument("--aligner", choices=["mafft", "muscle"], default="mafft", help="Alignment program to use (default: mafft).")
    args = parser.parse_args()

    # ----------------------
    # Create output directory if it does not exist
    # ----------------------
    os.makedirs(args.outdir, exist_ok=True)

    # ----------------------
    # Load sequences from query FASTA into a dictionary
    # ----------------------
    query_dict = {record.id: record for record in SeqIO.parse(args.query_fasta, "fasta")}

    # ----------------------
    # Load sequences from target FASTA into a dictionary
    # ----------------------
    target_dict = {record.id: record for record in SeqIO.parse(args.target_fasta, "fasta")}

    # ----------------------
    # Read the table line by line
    # ----------------------
    with open(args.table, "r") as infile:
        for line in infile:
            # Skip empty lines
            if not line.strip():
                continue

            # Split line into columns
            cols = line.strip().split()
            query_id = cols[0]   # First column: query sequence ID
            target_id = cols[1]  # Second column: target sequence ID

            # Check if both sequences exist in dictionaries
            if query_id not in query_dict:
                print(f"Warning: query ID {query_id} not found in query FASTA. Skipping.")
                continue
            if target_id not in target_dict:
                print(f"Warning: target ID {target_id} not found in target FASTA. Skipping.")
                continue

            # Create temporary FASTA file with the two sequences
            temp_fasta = os.path.join(args.outdir, f"{query_id}_{target_id}_temp.fasta")
            SeqIO.write([query_dict[query_id], target_dict[target_id]], temp_fasta, "fasta")

            # Define output alignment file path
            out_alignment = os.path.join(args.outdir, f"{query_id}_{target_id}.fasta")

            # ----------------------
            # Load HPC modules for aligner
            # ----------------------
            if args.aligner == "mafft":
                subprocess.run(["module", "load", "MAFFT/7.520-GCC-12.3.0-with-extensions"], shell=True, executable="/bin/bash")
            elif args.aligner == "muscle":
                subprocess.run(f"module load MUSCLE/3.8.1551-GCC-9.3.0 && muscle -in {temp_fasta} -out {out_alignment}",
               shell=True, executable="/bin/bash", check=True)

            # Run the chosen aligner
            if args.aligner == "mafft":
                subprocess.run(["mafft", "--auto", temp_fasta], stdout=open(out_alignment, "w"), check=True)
            elif args.aligner == "muscle":
                subprocess.run(["muscle", "-in", temp_fasta, "-out", out_alignment], check=True)

            # Remove the temporary FASTA file
            os.remove(temp_fasta)

# ======================
# Run main function
# ======================
if __name__ == "__main__":
    main()

