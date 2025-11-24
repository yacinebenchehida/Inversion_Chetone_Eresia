#!/usr/bin/env python3  # Interpreter directive

import os  # Import os for directory and file operations
import argparse  # Import argparse for external argument handling

def main():  # Define main function
    parser = argparse.ArgumentParser()  # Create argument parser
    parser.add_argument("input_folder", help="Folder containing all species directories")  # Add input folder argument
    parser.add_argument("--threshold", type=float, default=0.8, help="Proportion cutoff")  # Add threshold argument
    args = parser.parse_args()  # Parse arguments

    base = args.input_folder  # Store path to base folder
    species_dirs = [d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d))]  # List species directories

    presence = {}  # Create dictionary for gene presence
    species_gene_lists = {}  # Create dictionary to store gene sets per species

    for species in species_dirs:  # Loop over each species
        path = os.path.join(base, species)  # Construct species directory path
        files = [f for f in os.listdir(path) if f.endswith(".fna")]  # List files that end with .fna
        genes = set([os.path.splitext(f)[0] for f in files])  # Extract gene names without extension
        species_gene_lists[species] = genes  # Store gene set for species

        for g in genes:  # Loop over each gene
            if g not in presence:  # Check if gene not in presence dictionary
                presence[g] = 0  # Initialize count for gene
            presence[g] += 1  # Increment count for gene

    n_species = len(species_dirs)  # Count number of species
    cutoff = int(args.threshold * n_species)  # Compute cutoff number of species

    kept_genes = set([g for g, count in presence.items() if count >= cutoff])  # Create set of kept genes

    for species in species_dirs:  # Loop again over species
        path = os.path.join(base, species)  # Path to species directory
        for f in os.listdir(path):  # Loop over files in species directory
            if f.endswith(".fna"):  # Check file extension
                gene = os.path.splitext(f)[0]  # Extract gene name
                if gene not in kept_genes:  # Check if gene is not kept
                    os.remove(os.path.join(path, f))  # Delete file

if __name__ == "__main__":  # Main guard
    main()  # Run main function

