#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=mast_inv

###########################
# 0 - Load useful modules #
###########################
module load MMseqs2/17-b804f-gompi-2023b
module load R/4.2.1-foss-2022a

DATASETS="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/datasets"
MEL_GENES="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Syntheny/Inputs/Hmel_genes.txt"
BASE_RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Results"

#################################
# Function to process one genome
#################################
process_species() {
    local ACCESSION="$1"

    # 1 - Retrieve species and assembly type
    SPECIES=$("$DATASETS" summary genome accession "$ACCESSION" | jq -r '.reports[0].organism.organism_name')
    SPECIES_SAFE=$(echo "$SPECIES" | perl -pe 's/ /_/g')
    ASSEMBLY_TYPE_RAW=$("$DATASETS" summary genome accession "$ACCESSION" | jq -r '.reports[0].assembly_info.assembly_type')
    if [[ "$ASSEMBLY_TYPE_RAW" == "haploid" ]]; then
        ASSEMBLY="primary_assembly"
    elif [[ "$ASSEMBLY_TYPE_RAW" == "alternate-pseudohaplotype" ]]; then
        ASSEMBLY="alternate_assembly"
    else
        ASSEMBLY="$ASSEMBLY_TYPE_RAW"
    fi

    WORKDIR="$BASE_RESULTS/${ACCESSION}_${SPECIES}_${ASSEMBLY}"
    mkdir -p "$WORKDIR"

    echo "=== Processing $ACCESSION | $SPECIES | $ASSEMBLY ==="

    # 2 - Download and extract genome
    "$DATASETS" download genome accession "$ACCESSION" --include genome --filename "$WORKDIR/${ACCESSION}.zip"
    unzip -q "$WORKDIR/${ACCESSION}.zip" -d "$WORKDIR"
    REF=$(find "$WORKDIR" -type f -name "*.fna" | head -n 1)

    # 3 - Create MMseqs2 databases
    mmseqs createdb "$REF" "$WORKDIR/genome_db"
    mmseqs createdb "$MEL_GENES" "$WORKDIR/mel_genes_db"

    # 4 - Run search
    mmseqs search "$WORKDIR/mel_genes_db" "$WORKDIR/genome_db" "$WORKDIR/tmp_result" "$WORKDIR/tmp_dir" \
        --search-type 2 --threads 4

    # 5 - Convert results to TSV
    mmseqs convertalis "$WORKDIR/mel_genes_db" "$WORKDIR/genome_db" "$WORKDIR/tmp_result" \
        "$WORKDIR/${ACCESSION}_raw_blast_results.tsv" \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    # 6 - Keep best hits
    awk '!seen[$1]++' "$WORKDIR/${ACCESSION}_raw_blast_results.tsv" > "$WORKDIR/${ACCESSION}_best_hits.tsv"

    # 7 - Cleanup genome files to save space
    rm -rf "$WORKDIR/ncbi_dataset" "$WORKDIR/${ACCESSION}.zip" "$REF"

    echo "=== Finished $ACCESSION | $SPECIES | $ASSEMBLY ==="
}

##########################################
# Loop over all accessions passed as args
##########################################
for ACCESSION in "$@"; do
    process_species "$ACCESSION"
done
