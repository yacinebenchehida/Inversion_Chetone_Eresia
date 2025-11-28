#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=align

module load MAFFT/7.520-GCC-12.3.0-with-extensions
module load trimAl/1.4.1-GCC-9.3.0
PAL2NAL="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/pal2nal.v14/pal2nal.pl"
GENES=("$@")
WD="/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/busco_phylo/Results"
count_genes=0

cd $WD

for gene in "${GENES[@]}"; do

    outdir="$WD/alignement/$gene"
    mkdir -p "$outdir"
    count_genes=$((count_genes + 1))
    echo "$count_genes $gene"

    # collect protein sequences
    for sp in GC*; do
        if [[ -f $sp/$gene.faa ]]; then
            sed "s/>/>${sp}_/" "$sp/$gene.faa" | sed 's/:.*//' >> "$outdir/protein.faa"
        fi
    done
    echo $gene protein sequences extracted

    # collect dna sequences
    for sp in GC*; do
        if [[ -f $sp/$gene.fna ]]; then
            sed "s/>/>${sp}_/" "$sp/$gene.fna" | sed 's/:.*//' >> "$outdir/dna.fna"
        fi
    done
    echo $gene dna sequences extracted

    # skip empty genes
    if [[ ! -s "$outdir/protein.faa" ]]; then
        continue
    fi

    # align protein
    echo $gene starting alignment
    mafft --thread 1 --auto "$outdir/protein.faa" > "$outdir/protein.aln"
    echo $gene Alignement performed

    # dna alignment via pal2nal
    echo $gene starting pal2nal
    $PAL2NAL "$outdir/protein.aln" "$outdir/dna.fna" -output fasta -codontable 1 > "$outdir/codon.aln"
    echo $gene pal2nal finished

    # trimal alignment cleaning
    echo $gene Starting trimming
    trimal -in "$outdir/codon.aln" -out "$outdir/${gene}_trimmed.aln" -gt 0.8 -st 0.001 -resoverlap 0.75 -seqoverlap 80
    echo $gene Trimming done

    # Clean folder
    rm "$outdir/codon.aln" "$outdir/protein.aln" "$outdir/dna.fna" "$outdir/protein.faa" 
done
