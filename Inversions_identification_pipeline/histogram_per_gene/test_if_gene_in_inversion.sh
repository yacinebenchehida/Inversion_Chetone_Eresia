#!/bin/bash
#SBATCH --job-name=geninv              
#SBATCH --ntasks=1                          
#SBATCH --cpus-per-task=1                   
#SBATCH --mem=5gb                          
#SBATCH --time=0-00:20:00                   
#SBATCH --account=BIOL-SPECGEN-2018        

module load R/4.4.1-gfbf-2023b

Rscript test_if_gene_in_inversion.R $@
