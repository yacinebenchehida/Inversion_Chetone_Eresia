#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=Exo


########################################
# 0 - Load needed modules from the hpc #
########################################
module load parallel/20230722-GCCcore-12.3.0
module load Exonerate/2.4.0-GCC-11.3.0
module load Biopython/1.83-foss-2023a
module load AUGUSTUS/3.5.0-foss-2022a
module load MAFFT/7.505-GCC-11.3.0-with-extensions

MELPO_PROT="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/SnpEff/Inputs/melpomene_proteins_sequences.fasta"
REF_GENOME="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/Chetone_histrio/chetone_histrio_mtDNA_05_02_23.fasta"

#################################################
# 2 - Split each fasta entry in a separate file #
#################################################
csplit -z -f seq_ "$MELPO_PROT" '/^>/' '{*}'
for f in seq_*; do
    name=$(grep '^>' "$f" | sed 's/^>//;s/ .*//')
    mv "$f" "${name}.fa"
done

#########################################
# 3 - Function for parallel execution   #
#########################################
run_exonerate() {
    local i="$1"
    echo "Processing FASTA: $i"
    # running exonerate
    exonerate --model protein2genome "$i" "$REF_GENOME" \
        --querytype protein --bestn 1 --showvulgar no --showtargetgff yes \
        --showalignment no > "${i}_exonerate_output.txt"
    
    # extract only genes (for augustus later)
    python3 gene_exonerate2fasta.py \
       -f "${i}_exonerate_output.txt" \
        -g "$REF_GENOME" \
        -o "${i}_tmp_genes"

    # Get CDS directly from exonerate
    python3 exonerate2fasta_cds.py \
        -f "${i}_exonerate_output.txt" \
        -g "$REF_GENOME" \
        -o "${i}_cds"
    
    # Infer cds using prediction (from exonerate gene prediction)
    augustus --progress=true --strand=both \
             --species=heliconius_melpomene1 \
             "${i}_tmp_genes" > "${i}_augustus.gff"
    getAnnoFasta.pl "${i}_augustus.gff"

    # Compare (align) augustus infered by augustus with the initial melpomene protein
    mafft --maxiterate 1000 --localpair <(cat "${i}_augustus.aa" "$i")> "${i}_comparison.fa"

}
export -f run_exonerate
export REF_GENOME

#########################################
# 4 - Launch all jobs with parallel     #
#########################################
parallel -j 12 run_exonerate ::: *.fa
