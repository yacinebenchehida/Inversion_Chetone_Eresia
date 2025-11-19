#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=busco_chunk

module load SciPy-bundle/2025.07-gfbf-2025b
module load MMseqs2/17-b804f-gompi-2023b


GENE_DB="$1"
CHEMIN="$2"
shift 2


#################################
# Function to process one genome
#################################
process_one() {
    local FILE="$1"
    local BASE=$(basename "$FILE" | perl -pe 's/_genomic\.fna//')
    local RESULTS_PATH="$CHEMIN/$BASE"

    mkdir -p "$RESULTS_PATH"
    rm -rf "$RESULTS_PATH"/*

    mmseqs createdb "$FILE" "$RESULTS_PATH/genome_db"
    mmseqs createindex "$RESULTS_PATH/genome_db" "$RESULTS_PATH/tmp_dir"

    mmseqs search "$GENE_DB" "$RESULTS_PATH/genome_db" "$RESULTS_PATH/tmp_result" "$RESULTS_PATH/tmp_dir" \
        --search-type 2 --threads 4

    mmseqs convertalis "$GENE_DB" "$RESULTS_PATH/genome_db" "$RESULTS_PATH/tmp_result" \
        "$RESULTS_PATH/busco_raw_blast_results.tsv" \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    awk '!seen[$1]++' "$RESULTS_PATH/busco_raw_blast_results.tsv" > "$RESULTS_PATH/busco_best_hits.tsv"

    python3 /mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/reorder_busco_in_numerical_order.py \
        "$RESULTS_PATH/busco_best_hits.tsv" "$RESULTS_PATH/${BASE}_busco.tsv"

    rm -rf "$RESULTS_PATH"/tmp* \
       "$RESULTS_PATH"/genome_db* \
       "$RESULTS_PATH"/mel_genes_db*
    
}


#################################
# Run all genomes passed as args
#################################
for GENOME in "$@"; do
    process_one "$GENOME"
done
