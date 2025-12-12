#!/usr/bin/env bash

GENE_LIST=genes_2_keep  
BATCH_SIZE=10

mapfile -t GENES < "$GENE_LIST"
TOTAL=${#GENES[@]}

for ((i=0; i<$TOTAL; i+=BATCH_SIZE)); do
    chunk=("${GENES[@]:i:BATCH_SIZE}")
    sbatch 3_align.sh "${chunk[@]}"
done
