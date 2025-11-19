#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=busco_launch

module load MMseqs2/17-b804f-gompi-2023b

GENE_DB="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Inputs/busco_genes_on_Hmel215003o.faa"
CHEMIN="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Results/busco"
CHUNK=20
FILE=list_all_genomes.txt

# Create BUSCO DB once here
mmseqs createdb "$GENE_DB" "$CHEMIN/mel_genes_db"
mmseqs createindex "$CHEMIN/mel_genes_db" "$CHEMIN/tmp_dir_global"

# Split genomes into chunks and submit jobs
awk -v chunk="$CHUNK" '{
    i = int((NR-1)/chunk)
    groups[i] = groups[i] " " $0
}
END {
    for (i in groups) print groups[i]
}' "$FILE" | while read GROUP; do
    sbatch busco_batch.sh $CHEMIN/mel_genes_db $CHEMIN $GROUP
done
