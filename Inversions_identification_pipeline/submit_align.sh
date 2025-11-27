#!/usr/bin/env bash

GENE_LIST=busco_genes_99.txt  
BATCH_SIZE=20

mapfile -t GENES < "$GENE_LIST"
TOTAL=${#GENES[@]}

for ((i=0; i<$TOTAL; i+=BATCH_SIZE)); do
    chunk=("${GENES[@]:i:BATCH_SIZE}")
    sbatch align.sh "${chunk[@]}"
done
