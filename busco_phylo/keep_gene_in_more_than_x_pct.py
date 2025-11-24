#!/usr/bin/env python3  # Interpreter directive

import os  # Import os for directory operations
import argparse  # Import argparse for external arguments

def main():  # Define main function
    parser = argparse.ArgumentParser()  # Create argument parser
    parser.add_argument("--input_folder", required=True, help="Folder containing all species directories")  # Input folder argument
    parser.add_argument("--output_file", required=True, help="Output text file for the kept genes")  # Output file argument
    parser.add_argument("--threshold", type=float, default=0.8, help="Proportion cutoff")  # Threshold argument
    args = parser.parse_args()  # Parse arguments

    base = args.input_folder  # Store the base directory path
    species_dirs = [d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d))]  # List species directories

    presence = {}  # Dictionary for gene presence counts

    for species in species_dirs:  # Loop over species
        path = os.path.join(base, species)  # Path to species folder
        files = [f for f in os.listdir(path) if f.endswith(".fna")]  # List .fna files
        genes = set(os.path.splitext(f)[0] for f in files)  # Extract gene names without extension

        for g in genes:  # Loop over each gene
            if g not in presence:  # Initialize count if needed
                presence[g] = 0  # Set count to zero
            presence[g] += 1  # Increase count for this gene

    n_species = len(species_dirs)  # Count total species
    cutoff = int(args.threshold * n_species)  # Compute minimum number of species needed

    kept_genes = [g for g, count in presence.items() if count >= cutoff]  # List genes that meet threshold

    with open(args.output_file, "w") as out:  # Open output file
        for gene in kept_genes:  # Loop over kept genes
            out.write(gene + "\n")  # Write gene name

if __name__ == "__main__":  # Main guard
    main()  # Run main function

