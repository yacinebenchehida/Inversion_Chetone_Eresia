#!/usr/bin/env python3

import sys
from Bio import SeqIO

path = "/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/busco_phylo/tree/supermatrix.fasta"

print("Checking:", path)

lengths = set()
names = set()
duplicates = []
illegal = {}

valid = set("ACGTNacgtn-.?")

count = 0

for rec in SeqIO.parse(path, "fasta"):
    count += 1
    seq = str(rec.seq)

    # record lengths
    lengths.add(len(seq))

    # check duplicate names
    if rec.id in names:
        duplicates.append(rec.id)
    else:
        names.add(rec.id)

    # scan illegal chars (fast)
    bad = set(c for c in seq if c not in valid)
    if bad:
        illegal[rec.id] = bad

print("\nNumber of sequences:", count)
print("Unique sequence lengths:", lengths)

if len(lengths) != 1:
    print("WARNING: Not all sequences have equal length!")
else:
    print("All sequences same length.")

if duplicates:
    print("\nDUPLICATE NAMES FOUND:")
    for d in duplicates[:20]:
        print(" ", d)
    print(f"(total duplicates: {len(duplicates)})")
else:
    print("\nNo duplicate names.")

if illegal:
    print("\nILLEGAL CHARACTERS FOUND:")
    for k,v in list(illegal.items())[:20]:
        print(" ", k, ":", v)
    print(f"(total sequences with illegal chars: {len(illegal)})")
else:
    print("\nNo illegal characters.")

print("\nDone.")
