#!/bin/bash

while getopts "v:o:s:a:e:w:i:" opt; do
  case $opt in
    v) VCF="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    s) SCAFFOLD="$OPTARG" ;;
    a) START="$OPTARG" ;;
    e) END="$OPTARG" ;;
    w) WIN_SIZE="$OPTARG" ;;
    i) SAMPLES="$OPTARG" ;;
    g) GFF="$optarg" ;;
  esac
done

mkdir -p $OUTDIR/logs $OUTDIR/vcf $OUTDIR/results

POS=$START
while [ $POS -lt $END ]; do
    WIN_START=$POS
    WIN_END=$(( POS + WIN_SIZE - 1 ))
    [ $WIN_END -gt $END ] && WIN_END=$END

    sbatch mutations_count.sh \
           -v $VCF \
           -s $SCAFFOLD \
           -t $WIN_START \
           -e $WIN_END \
           -o $OUTDIR \
           -m $SAMPLES \
           -g $GFF \

    POS=$(( WIN_START + WIN_SIZE ))
done

# Submit R plotting job after all workers finish
#sbatch --dependency=singleton plot_results.sh \
#       -o $OUTDIR \
#       -s $SCAFFOLD


  #./master_snpeff.sh -v ../Results/ctg001430_1_np1212/annotated_ctg001430_1_np1212.vcf.gz -o ../Results/montest -s ctg001430_1_np1212 -a 1000000 -e 5000000 -w 500000 -i all_samples.txt -G /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/SnpEff/Inputs/Chetone_histrio/genes.gff