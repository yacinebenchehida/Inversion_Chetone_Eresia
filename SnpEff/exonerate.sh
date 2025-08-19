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

REF_GENOME="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/Chetone_histrio/chetone_histrio_mtDNA_05_02_23.fasta"

#########################################
# 3 - Function for parallel execution   #
#########################################
run_exonerate() {
    local i="$1"
    echo "Processing FASTA: $i"
    # running exonerate
    exonerate --model protein2genome "$i" "$REF_GENOME" \
        --querytype protein --bestn 1 --showvulgar no --showtargetgff yes \
        --showalignment no > "${i}.pseudogff"
}
export -f run_exonerate
export REF_GENOME

#########################################
# 4 - Launch all jobs with parallel     #
#########################################
parallel -j 12 run_exonerate ::: *.fa
