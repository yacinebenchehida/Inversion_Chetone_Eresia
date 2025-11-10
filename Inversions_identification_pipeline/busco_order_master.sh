#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=26:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=busco_order

module load  SciPy-bundle/2025.07-gfbf-2025b
module load MMseqs2/17-b804f-gompi-2023b

# Create the BUSCO genes database once
mmseqs createdb ../../Inversions_identification_pipeline/Inputs/busco_genes_on_Hmel215003o.faa mel_genes_db
mmseqs createindex mel_genes_db tmp_dir_global

readlink -f /mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Inputs/ref_genome/GC*/*fna | while read line
do
    nom=$(basename "$line")
    chemin=$(dirname "$line")

    # Remove any previous temporary data
    rm -rf genome_db tmp_dir tmp_result* busco_raw_blast_results.tsv busco_best_hits.tsv genome* 

    # Create genome database fresh
    mmseqs createdb "${chemin}/${nom}" genome_db
    mmseqs createindex genome_db tmp_dir

    # Protein vs nucleotide search
     mmseqs search mel_genes_db genome_db tmp_result tmp_dir --search-type 2 --threads 4

     mmseqs convertalis mel_genes_db genome_db tmp_result busco_raw_blast_results.tsv  --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

     awk '!seen[$1]++' busco_raw_blast_results.tsv > busco_best_hits.tsv

     python3 ../../Inversions_identification_pipeline/Scripts/reorder_busco_in_numerical_order.py busco_best_hits.tsv "${nom}_busco.tsv"

done
