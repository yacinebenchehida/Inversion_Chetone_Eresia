#!/bin/bash

module load MAFFT/7.520-GCC-12.3.0-with-extensions

fasta1=$1
fasta2=$2
outdir=$3
mkdir -p "$outdir"

# Split input1
awk '/^>/{if(f){close(f)}; f=sprintf("a_%03d.fa", ++n)} {print > f}' "$fasta1"
n1=$(ls a_*.fa | wc -l)

# Split input2
awk '/^>/{if(f){close(f)}; f=sprintf("b_%03d.fa", ++m)} {print > f}' "$fasta2"
n2=$(ls b_*.fa | wc -l)

# Pairwise align all combinations: a_i vs b_j
for i in $(seq -f "%03g" 1 $n1); do
  for j in $(seq -f "%03g" 1 $n2); do
    (cat a_"$i".fa; cat b_"$j".fa) > tmp.fa
    mafft --maxiterate 1000 --localpair tmp.fa > "$outdir"/pair_"$i"_"$j".aln.fa
  done
done

# Cleanup
rm a_*.fa b_*.fa tmp.fa
