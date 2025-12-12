#!/bin/bash
#SBATCH --job-name=busco       
#SBATCH --partition=nodes
#SBATCH --ntasks=1                       
#SBATCH --cpus-per-task=4                
#SBATCH --mem=15GB                      
#SBATCH --time=6:00:00                  
#SBATCH --account=BIOL-SPECGEN-2018      

############################################
# 0 - DEFINE WORKING DIRECTORIES AND PATHS #
############################################
INPUT="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Inputs/ref_genome/${1}"
REF=$(ls $INPUT|grep -E "*.fa$|*.fasta$|*.fna$"|grep -E $1)
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/busco_phylo/Results/${1}"

mkdir -p $RESULTS

#######################################
# 1 - LOAD BUSCO AND ASSOCIATED TOOLS #
#######################################
module load MetaEuk/4-GCC-10.2.0
module load AUGUSTUS/3.5.0-foss-2022a
module load BLAST+/2.14.0-gompi-2022b
module load BBMap/38.90-GCC-10.2.0
module load BUSCO/5.4.3-foss-2020b

#################       
# 2 - RUN BUSCO #
#################
cd $RESULTS
busco -m genome -i $INPUT/$REF -o $1 --out_path $RESULTS -l lepidoptera_odb10 -c 4

#############################################    
# 3 - KEEP GENES AND PROTEIN SEQUENCES ONLY #
#############################################
find "$RESULTS/$1/run_lepidoptera_odb10/busco_sequences/single_copy_busco_sequences/" -name "*.fna" -exec mv {} "$RESULTS" \;
find "$RESULTS/$1/run_lepidoptera_odb10/busco_sequences/single_copy_busco_sequences/" -name "*.faa" -exec mv {} "$RESULTS" \;

rm -rf $RESULTS/busco_downloads
rm -rf $RESULTS/$1
