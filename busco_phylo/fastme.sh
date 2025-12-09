#!/bin/bash
#SBATCH --time=2-00:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=fastme

############################
# load tools and set paths #
############################
module load FastME/2.1.6.1-GCC-10.2.0
module load Armadillo/12.6.2-foss-2023a
PXRR="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/phyx/src/pxrr"


#############################
# Get NJ tree using fastme  #
#############################
fastme -i ../tree/supermatrix.phy -o njtree_10_bootstraps.nwk -n -d -T 32 -b 10



##############################
# Root tree with caddisflies #
##############################
$PXRR -t njtree_10_bootstraps.nwk -g GCA_947579605.1,GCA_965644405.1,GCA_964276685.1 -o rooted_tree_NJ_10_boot.nwk
