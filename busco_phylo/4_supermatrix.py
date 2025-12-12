#!/usr/bin/env python3                                # use python3

import os                                            # filesystem ops
from Bio import SeqIO                                # fasta parsing
from tqdm import tqdm                                 # progress bars
import re                                             # regex species extraction

aligned_dir = "/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/busco_phylo/Results/alignement"  # alignment dir

genes = sorted(os.listdir(aligned_dir))              # list of gene folders

species = []                                         # final species list
seen = set()                                         # track unique species

def extract_species_id(header):                      # extract species/accession from FASTA header
    m = re.match(r"^(GC[AF]_\d+\.\d+)", header)      # match GCA/GCF accession at start
    if m:
        return m.group(1)
    m = re.match(r"^(Chetone|Kallima|Eresia)_", header)  # match genus-prefixed scaffolds
    if m:
        return m.group(1)
    return None

# -------------------- PASS 1: extract species FROM FASTA HEADERS --------------------
for g in tqdm(genes, desc="Scanning species"):        # loop genes
    aln_path = f"{aligned_dir}/{g}/{g}_trimmed.aln"   # alignment file path
    if not os.path.exists(aln_path):                  # skip missing
        continue
    for rec in SeqIO.parse(aln_path, "fasta"):        # loop sequences
        sid = extract_species_id(rec.id)              # extract species id
        if sid:                                       # if match
            if sid not in seen:                       # new species
                seen.add(sid)                         # add to set
                species.append(sid)                   # preserve order

# -------------------- CREATE TEMP FILES --------------------
tmp_dir = "tmp_supermatrix_parts"                    # folder for partial sequences
os.makedirs(tmp_dir, exist_ok=True)                  # create folder

tmp_files = {}                                       # species -> file handle
for sp in species:                                   # loop species
    tmp_path = f"{tmp_dir}/{sp}.part"                # partial file path
    tmp_files[sp] = open(tmp_path, "w")              # open file

# -------------------- PASS 2: build supermatrix --------------------
for g in tqdm(genes, desc="Concatenating genes"):     # loop genes
    aln_path = f"{aligned_dir}/{g}/{g}_trimmed.aln"   # alignment file
    if not os.path.exists(aln_path):                  # skip missing
        continue

    records = list(SeqIO.parse(aln_path, "fasta"))    # load alignment
    if not records:                                   # skip empty
        continue

    L = len(records[0].seq)                           # alignment length
    seqs = {}                                         # species -> sequence

    for rec in records:                               # loop sequences
        sid = extract_species_id(rec.id)              # extract species id
        if sid:                                       # if match
            seqs[sid] = str(rec.seq)                  # store sequence

    for sp in species:                                # loop species
        fragment = seqs.get(sp, "-" * L)              # real seq or gaps
        tmp_files[sp].write(fragment)                 # write fragment

# -------------------- CLOSE TEMP FILES --------------------
for f in tmp_files.values():                          # loop file handles
    f.close()                                         # close file

# -------------------- MERGE INTO FINAL FASTA --------------------
with open("supermatrix.fasta", "w") as out:           # output FASTA
    for sp in species:                                # loop species
        out.write(f">{sp}\n")                         # write header
        tmp_path = f"{tmp_dir}/{sp}.part"             # partial path
        with open(tmp_path) as tf:                    # open partial
            out.write(tf.read() + "\n")               # write sequence

# -------------------- CLEANUP --------------------
for sp in species:                                    # loop species
    os.remove(f"{tmp_dir}/{sp}.part")                 # delete partial

os.rmdir(tmp_dir)                                    # delete tmp folder
