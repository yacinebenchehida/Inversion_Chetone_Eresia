#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=supmatrix

module load Biopython/1.84-foss-2024a
module load tqdm/4.66.5-GCCcore-13.3.0

python3 ./supermatrix.py
