#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=Plot

module load R/4.2.1-foss-2022a

Rscript ./plot.R ${1}/${2} $2 $3 $4 $5
