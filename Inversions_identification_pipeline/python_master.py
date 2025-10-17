#!/usr/bin/env python3
import subprocess, json
from itertools import islice

# 1. Run datasets command
subprocess.run([
    "/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/datasets",
    "summary", "genome", "taxon", "Lepidoptera[Organism]",
    "--assembly-level", "chromosome",
    "--reference"
], stdout=open("jason_list_lepi_ref_genome.txt", "w"))

# 2. Parse JSON and extract accessions
with open("jason_list_lepi_ref_genome.txt") as f:
    data = json.load(f)
accessions = [r["accession"] for r in data["reports"]]
with open("chromosome_assembled_lepi.txt", "w") as f:
    f.write("\n".join(accessions))

# 3. Batch accessions (10 per sbatch job), only first 2 batches
script = "inversion_script.sh"
batch_size = 3
max_batches = 2
batch_count = 0

with open("chromosome_assembled_lepi.txt") as f:
    while batch_count < max_batches:
        batch = list(islice(f, batch_size))
        if not batch:
            break
        batch = [x.strip() for x in batch if x.strip()]
        subprocess.run(["sbatch", script] + batch)
        batch_count += 1
