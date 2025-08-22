#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=05:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=snpeff

# Input variables
SNPEFF="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/old_snpeff/snpEff.jar"
WD="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/SnpEff/Inputs/Chetone_histrio"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/SnpEff/Results"
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/Chetone_histrio/Chetone_histrio.vcf.gz"
SCAF=$1

# Build database
#java -jar $SNPEFF build -gff3 -v Chetone_histrio

# Create VCF for the inversion region
#module load BCFtools/1.19-GCC-13.2.0
# Inversion
#bcftools view --regions ctg001860_1_np1212:9835233-10873574 -e 'F_MISSING > 0.25' $VCF -o $WD/ctg001860_1_np1212_9835233_10873574_filtered_output.vcf
#bgzip $WD/ctg001860_1_np1212_9835233_10873574_filtered_output.vcf
#tabix -f $WD/ctg001860_1_np1212_9835233_10873574_filtered_output.vcf.gz
#VCF="$WD/ctg001860_1_np1212_9835233_10873574_filtered_output.vcf.gz"

# unique scaffold
#bcftools view --regions $SCAF -e 'F_MISSING > 0.25' $VCF -o $WD/${SCAF}_filtered_output.vcf
#bgzip $WD/${SCAF}_filtered_output.vcf
#tabix $WD/${SCAF}_filtered_output.vcf.gz

# whole genome
# Annotate the VCF with snpeff
module load Java/21.0.2
mkdir -p $RESULTS/$SCAF

java -jar $SNPEFF -lof -stats $RESULTS/snpeff_summary_${SCAF}.html Chetone_histrio $VCF > $RESULTS/$SCAF/annotated_$SCAF.vcf
