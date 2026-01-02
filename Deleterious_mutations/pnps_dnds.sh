#!/bin/bash
#SBATCH --job-name=snpeff
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018 

# Define snpeff path and load modules
SNPEFF="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/old_snpeff/snpEff.jar"
module load Java/21.0.2
module load BCFtools/1.19-GCC-13.2.0

# Parse input parameters
while getopts "v:o:w:i:p:" opt; do
  case $opt in
    v) VCF="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    w) WINDOW_SIZE="$OPTARG" ;;
    i) SAMPLES="$OPTARG" ;;
    p) PATTERN="$OPTARG" ;;
  esac
done

# Determine scaffold to analyse + scaffold size
SCAFFOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $OUTDIR/contigs_size.txt| awk '{print $1}')
SIZE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $OUTDIR/contigs_size.txt| awk '{print $2}')

# Create working directories
mkdir -p $OUTDIR/$SCAFFOLD
OUTDIR="$OUTDIR/$SCAFFOLD"
echo working directories created

# Build database
java -jar $SNPEFF build -gff3 -v Arctia_plantaginis # must be adjusted to the species of intersest
echo snpeff data base created

# Create a vcf for the region to analyse
cat $SAMPLES| grep $PATTERN| awk '{print $1}'  > "$OUTDIR/samples.txt"

#cat $SAMPLES|grep -E $GROUP|awk '{print $1}'| head -n 1 > $OUTDIR/samples.txt
SAMPLES="$OUTDIR/samples.txt"
bcftools view --regions $SCAFFOLD -S $SAMPLES $VCF -Oz -o $OUTDIR/${SCAFFOLD}_raw_filtered_output.vcf.gz
tabix $OUTDIR/${SCAFFOLD}_raw_filtered_output.vcf.gz
VCF="$OUTDIR/${SCAFFOLD}_raw_filtered_output.vcf.gz"
echo raw vcf ready

# Create a snpeff vcf (with mutations annotation) for the region of interest s
java -jar $SNPEFF -lof -stats  $OUTDIR/snpeff_summary_${SCAFFOLD}.html Heliconius_numata $VCF > $OUTDIR/annotated_$SCAFFOLD.vcf #  must be changed to the species of interest
#bcftools view -i 'FMT/GQ>=10 && FMT/DP>6' $OUTDIR/annotated_$SCAFFOLD.vcf -o $OUTDIR/GQ10_annotated_$SCAFFOLD.vcf -O v
bgzip $OUTDIR/annotated_$SCAFFOLD.vcf
tabix $OUTDIR/annotated_$SCAFFOLD.vcf.gz
VCF="$OUTDIR/annotated_$SCAFFOLD.vcf.gz"

# Keep genes larger than 5Kb
Genes_2_keep="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/SnpEff/Input/gene_regions_more5kb_numata.txt"
bcftools view -R $Genes_2_keep $VCF > $OUTDIR/5kb_plus_genes_$SCAFFOLD.vcf
bgzip $OUTDIR/5kb_plus_genes_$SCAFFOLD.vcf
tabix $OUTDIR/5kb_plus_genes_$SCAFFOLD.vcf.gz
VCF="$OUTDIR/5kb_plus_genes_$SCAFFOLD.vcf.gz"

rm $OUTDIR/*html $OUTDIR/*raw_filtered_output*  $OUTDIR/annotated_$SCAFFOLD.vcf*
echo snpeff mutation vcf ready

# Transform to a genotype file
./VCF2GenoNumb.pl -i "$VCF" -o "$OUTDIR/$SCAFFOLD.geno"
awk '{count=gsub(/\./,"."); if(count<4) print $0}' "$OUTDIR/$SCAFFOLD.geno" > "$OUTDIR/filtered_$SCAFFOLD.geno" # Keep only sites with less than 40% of missing data
# rm "$OUTDIR/$SCAFFOLD.geno"

#  Functio that processes the region of interest
process_region() {
    local VCF="$1"
    local OUTDIR="$2"
    local SCAFFOLD="$3"
    local START="$4"
    local END="$5"

    # Get the positions of genes of interest (exons)
    echo -e "${SCAFFOLD}\t${START}\t${END}" > "$OUTDIR/gene.pos"

    # Count the number of mutations and their type in the genes
    ./CountMutSubSlidingGeneV2.pl -i "$OUTDIR/$SCAFFOLD.geno" -o "$OUTDIR/${SCAFFOLD}_${START}_${END}.CountDnDs.tmp" -g "$OUTDIR/gene.pos"
    sed -n '2p' "$OUTDIR/${SCAFFOLD}_${START}_${END}.CountDnDs.tmp" > "$OUTDIR/${SCAFFOLD}_${START}_${END}.CountDnDs"

    # Cleaning folder
    #rm "$OUTDIR/${SCAFFOLD}_${START}_${END}.CountDnDs.tmp" 

}

# Loop over scaffold in non-overlapping windows
for (( START_POS=1; START_POS<=SIZE; START_POS+=WINDOW_SIZE )); do
    END_POS=$(( START_POS + WINDOW_SIZE - 1 ))
    if (( END_POS > SIZE )); then
        END_POS=$SIZE
    fi

    # Call the function for this window
    process_region "$VCF" "$OUTDIR" "$SCAFFOLD" "$START_POS" "$END_POS" 
done

# Handle last window explicitly if SIZE is not a multiple of WINDOW_SIZE
REMAINDER=$(( SIZE % WINDOW_SIZE ))
if (( REMAINDER != 0 )); then
    LAST_START=$(( SIZE - REMAINDER + 1 ))
    LAST_END=$SIZE
    process_region "$VCF" "$OUTDIR" "$SCAFFOLD" "$LAST_START" "$LAST_END" 
fi

#  Combine results for the scaffold
awk -v p="$PATTERN" '{print $0, p}' "$OUTDIR"/*CountDnDs > "$OUTDIR"/"$SCAFFOLD"_combined_CountDnDs.txt

# Remove tmp large vcf from the working folder
rm "$OUTDIR"/*gz "$OUTDIR"/*tbi "$OUTDIR"/snpeff_summary* "$OUTDIR"/*CountDnDs "$OUTDIR"/*geno

#  Plot the data with the R script
# Usage: ./pnps_dnds.sh -v /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/Arctia_plantaginis/Arctia_plantaginis.DP3_GQ3_QUAL3.indels_snps.vcf.gz -o ../Results/CADEBD010000171.1 -s CADEBD010000171.1 -w 100000 -i pop_info.txt -p histrio-Peru
