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

# 2. Parse JSON and extract accessions and metadata
with open("jason_list_lepi_ref_genome.txt") as f:
    data = json.load(f)

accessions = []
with open("metadata_ref_genome.txt", "w") as meta_out:
    for report in data["reports"]:
        accession = report.get("accession", "")
        organism = report.get("organism", {}).get("organism_name", "")
        assembly_type = report.get("assembly_info", {}).get("assembly_type", "")
        meta_out.write(f"{accession}\t{organism}\t{assembly_type}\n")
        accessions.append(accession)

# 3. Batch accessions (10 per sbatch job), only first 2 batches
script = "inversion_script.sh"
batch_size = 2
max_batches = 5
batch_count = 0

with open("chromosome_assembled_lepi.txt", "w") as f:
    f.write("\n".join(accessions))

with open("chromosome_assembled_lepi.txt") as f:
    while batch_count < max_batches:
        batch = list(islice(f, batch_size))
        if not batch:
            break
        batch = [x.strip() for x in batch if x.strip()]
        subprocess.run(["sbatch", script] + batch)
        batch_count += 1
