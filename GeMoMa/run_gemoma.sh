#!/bin/bash
#SBATCH --job-name=gemoma
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00
#SBATCH --mem=50G

module load MMseqs2/17-b804f-gompi-2023b
module load snakemake/8.4.2-foss-2023a
module load BLAST+/2.14.1-gompi-2023a
module load parallel/20230722-GCCcore-12.3.0

GEMOMA="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/gemoma/bin/GeMoMa"
export GEMOMA

############
# Run STAR #
############
#snakemake --snakefile star.smk --cores 8 --printshellcmds

#############################
# Run GeMoMA on RNAseq data #
#############################
BAM="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/GeMoMa/Results/star/aln.sorted.bam"
WORKING_DIR="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/GeMoMa/Results/GeMoMa"
GENOME="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/Chetone_histrio/chetone_histrio_mtDNA_05_02_23.fasta"

echo STARTING
#Extract introns (extract splice junctions/intron evidence from your RNA-seq BAM (creates introns.gff)
#$GEMOMA ERE s=FR_UNSTRANDED m=${BAM} v=LENIENT c=true mmq=40 mc=1 e=0 u=false

#Check introns (sanity-checks intron coordinates against the target genome and reports distributions/metrics.)
#$GEMOMA CheckIntrons t=${GENOME} i=introns.gff > intron_dist.txt

#Denoise  (filters low-confidence/noisy introns (produces denoised_introns.gff) to keep only reliable splice evidence.)
#$GEMOMA DenoiseIntrons i=introns.gff c=UNSTRANDED coverage_unstranded=coverage.bedgraph m=500000

#Search positions of reference transcripts in target genome
#make mmseqs2 database for target genome
#mkdir -p $WORKING_DIR/Chetone_histrio
#mmseqs createdb ${GENOME} $WORKING_DIR/Chetone_histrio/GenomeDB -v 2

############################################################################
# Run GeMoMA on proteomes (species: bombyx, melpomene, arctia plantagiris) #
############################################################################
# Define the function that runs GeMoMa analysis for one species
Gemoma_proteom() {
    echo SOMETHING HAPPENING
    # Read three arguments: species name, fasta file path, gff file path
    local species=$1  # species name
    local fasta=$2    # path to reference fasta file
    local gff=$3      # path to reference gff file

    # Define the output directory by combining WORKING_DIR and species name
    local OUTDIR="${WORKING_DIR}/${species}"

    # Create the output directory if it does not already exist
    mkdir -p "${OUTDIR}"

    # Run GeMoMa Extractor with annotation GFF and genome FASTA
    echo starting GEMOMA Extractor
     #$GEMOMA Extractor a="${gff}" g="${fasta}" outdir="${OUTDIR}"
    echo GEMOMA Extractor RAN
    # Create a QUERY subdirectory inside output directory for mmseqs2 database
    mkdir -p "${OUTDIR}/QUERY"

    # Create an mmseqs2 database from cds-parts.fasta located in output directory
     #mmseqs createdb "${OUTDIR}/cds-parts.fasta" "${OUTDIR}/QUERY/cdsDB" -v 2
    echo mmseqs create db $species

    # Run mmseqs search: query is cdsDB, target is GenomeDB
     #mmseqs search "${OUTDIR}/QUERY/cdsDB" $WORKING_DIR/Chetone_histrio/GenomeDB \
       #"${OUTDIR}/${species}_mmseqs.out" "${OUTDIR}/mmseqs_tmp" \
       #-e 100.0 --threads 4 -s 8.5 -a --comp-bias-corr 0 \
       #--max-seqs 500 --mask 0 --orf-start-mode 1 -v 2
    echo mmseqs search $species

    # Convert mmseqs output to a formatted text file with specified columns
     #mmseqs convertalis "${OUTDIR}/QUERY/cdsDB" $WORKING_DIR/Chetone_histrio/GenomeDB \
      # "${OUTDIR}/${species}_mmseqs.out" "${OUTDIR}/${species}_search.txt" \
       #--threads 4 \
       #--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen" -v 2
    echo mmseqs convertalis $species

    # Run the main GeMoMa gene prediction command with parameters and inputs
    java_gemoma="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/gemoma/share/gemoma-1.9-0/GeMoMa-1.9.jar"
    java -Xmx50G -jar $java_gemoma CLI GeMoMa \
      s="${OUTDIR}/${species}_search.txt" \
      c="${OUTDIR}/cds-parts.fasta" a="${OUTDIR}/assignment.tabular" \
      t="${GENOME}" sort=true Score=ReAlign i=denoised_introns.gff \
      coverage=UNSTRANDED coverage_unstranded=coverage.bedgraph \
      outdir="${OUTDIR}"
      echo "gemoma on proteom done"
}

# Export the function so GNU parallel can access it in subshells
export -f Gemoma_proteom

# Export variables used inside the function to environment
export WORKING_DIR GENOME

# Run the function in parallel using GNU parallel
cat /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/GeMoMa/Inputs/references.tsv |awk 'NR > 1'|parallel -j 4 --line-buffer --colsep '\t' Gemoma_proteom {1} {2} {3} 

####################################################
# Combine GeMoMa proteom evidence in a single file #
####################################################
echo STARTING MERGING EVIDENCE FROM PROTEOM FOR $@
GAF_ARGS=()
for sp in "$@"; do
    GAF_ARGS+=("p=${sp}" "g=${WORKING_DIR}/${sp}/predicted_annotation.gff")
done
echo "${GAF_ARGS[@]}"
$GEMOMA GAF "${GAF_ARGS[@]}"

#######################
# Add UTR information #
#######################
echo Adding UTR
$GEMOMA AnnotationFinalizer \
    g="${GENOME}" \
    a=filtered_predictions.gff \
    u=YES \
    i=denoised_introns.gff \
    c=UNSTRANDED \
    coverage_unstranded=coverage.bedgraph \
    rename=NO

cat final_annotation.gff |grep -P  "gene\t"|grep ctg001860|awk '$4 > 9850688'| awk '$4 < 10873586'
cat final_annotation.gff |grep -P  "gene\t"|grep ctg001860|awk '$4 > 9850688'| awk '$4 < 10873586'|awk '{print $5-$4}
