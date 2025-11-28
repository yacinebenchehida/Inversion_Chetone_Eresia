#!/usr/bin/env python3                           # use python3 interpreter

import os                                        # file system operations
from Bio import SeqIO                            # read fasta files
from tqdm import tqdm                             # progress bars

aligned_dir = "/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/busco_phylo/Results/alignement" # directory containing gene folders
genes = sorted(os.listdir(aligned_dir))           # list all gene folders in sorted order

species = set()                                   # set to collect all species names

for g in tqdm(genes, desc="Scanning species", unit="genes"):    # loop through genes with progress bar
    aln = f"{aligned_dir}/{g}/{g}_trimmed_protein.aln"              # path to trimmed alignment for gene g
    if os.path.exists(aln):                                     # skip if trimmed alignment does not exist
        for rec in SeqIO.parse(aln, "fasta"):                   # read all sequences from this alignment
            species.add(rec.id)                                 # add species ID to set

species = sorted(species)                                       # sort species list for stable ordering
supermatrix = {sp: "" for sp in species}                        # map species to their growing concatenated sequence

for g in tqdm(genes, desc="Concatenating", unit="genes"):       # second loop: build supermatrix
    aln_path = f"{aligned_dir}/{g}/{g}_trimmed_protein.aln"         # path to trimmed alignment
    if not os.path.exists(aln_path):                            # skip missing alignments
        continue

    records = list(SeqIO.parse(aln_path, "fasta"))              # read alignment as list of SeqRecord
    lengths = {rec.id: len(rec.seq) for rec in records}         # map species to sequence length
    L = max(lengths.values())                                   # alignment length (max)

    seqs = {rec.id: str(rec.seq) for rec in records}            # map species to sequence string

    for sp in species:                                          # iterate in fixed species order
        if sp in seqs:                                          # species present in this alignment
            supermatrix[sp] += seqs[sp]                         # append actual sequence
        else:                                                   # species missing for this gene
            supermatrix[sp] += "-" * L                          # append gap padding

with open("supermatrix.fasta", "w") as out:                     # write final concatenated supermatrix
    for sp in species:                                          # loop species in order
        out.write(f">{sp}\n{supermatrix[sp]}\n")                # fasta header and full concatenated sequence
