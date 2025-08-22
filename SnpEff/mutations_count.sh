#!/bin/bash
#SBATCH --job-name=worker_${SCAFFOLD}_${WIN_START}_${WIN_END}
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018 

while getopts "v:s:t:e:o:m:c:g:" opt; do
  case $opt in
    v) VCF="$OPTARG" ;;
    s) SCAFFOLD="$OPTARG" ;;
    t) WIN_START="$OPTARG" ;;
    e) WIN_END="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    m) SAMPLES="$OPTARG" ;;
    g) GFF="$OPTARG" ;;
  esac
done

module load BCFtools/1.19-GCC-13.2.0
module load rnaQUAST/2.3.0-foss-2023a
# Extract window and subset samples
bcftools view \
  -r ${SCAFFOLD}:${WIN_START}-${WIN_END} \
  -S $SAMPLES \
  $VCF \
  -Oz \
  -o $OUTDIR/vcf/${SCAFFOLD}_${WIN_START}_${WIN_END}.vcf.gz 

python3 ./mut_counter.py --vcf $VCF \
  --impact LOW \
  --scaffold $SCAFFOLD \
  --start $WIN_START  \
  --end $WIN_END \
  --gff $GFF \
  --outdir $OUTDIR
