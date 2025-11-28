

#!/usr/bin/env python3  # Use Python 3 interpreter

import os               # Import for path operations
from Bio import SeqIO    # Import for FASTA parsing

BASE = "/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/busco_phylo/Results/alignement"   # Base folder containing all gene subdirectories

genes = sorted(os.listdir(BASE))  # Get sorted list of all gene folders

species = set()  # Create an empty set to collect all species names

for g in genes:   # Loop over each gene folder
    aln = f"{BASE}/{g}/{g}_trimmed_protein.aln"  # Path to trimmed alignment file
    if os.path.exists(aln):  # Check if this alignment file exists
        for rec in SeqIO.parse(aln, "fasta"):  # Parse each sequence record
            species.add(rec.id)  # Add species header to the species set

species = sorted(species)  # Sort species names alphabetically

supermatrix = {sp: "" for sp in species}  # Initialize a dictionary: species → empty alignment string

for g in genes:  # Loop again over genes
    aln_path = f"{BASE}/{g}/{g}_trimmed_protein.aln"  # Path to trimmed alignment
    if not os.path.exists(aln_path):  # Skip if file not found
        continue

    records = list(SeqIO.parse(aln_path, "fasta"))  # Load all sequences in alignment

    lengths = {rec.id: len(rec.seq) for rec in records}  # Record length per species present
    L = max(lengths.values())  # Determine alignment length for this gene

    seqs = {rec.id: str(rec.seq) for rec in records}  # Store sequences as strings

    for sp in species:  # Loop over all species
        if sp in seqs:  # If species present in this alignment
            supermatrix[sp] += seqs[sp]  # Append its sequence
        else:  # If species absent
            supermatrix[sp] += "-" * L   # Append gap block of length L

with open("supermatrix.fasta", "w") as out:  # Write final combined matrix
    for sp in species:  # Loop over species
        out.write(f">{sp}\n{supermatrix[sp]}\n")  # Write header and sequence

