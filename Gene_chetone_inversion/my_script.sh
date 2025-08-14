# EXTRACT GENES FASTA IN INVERSION FROM THE GFF OF MELPOMENE
INPUTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Inputs/Heliconius_melpomene"
GFF="Heliconius_melpomene.gff3"
PROTEOME="Heliconius_melpomene.fa"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/Gene_inversion_Chetone/Results"

cat $INPUTS/$GFF|grep -E "Hmel215003o"|grep -P "\tgene\t"| awk '$4 > 1312092 && $4 < 1796339'|awk '{print $9}'|perl -pe 's/ID=//g'|while read line; do grep -A 1 $line $INPUTS/$PROTEOME; done  > $RESULTS/Melpomene_gene_inversion.fasta

# CONVERT NUCLEOTIDE TO PROTEIN FOR BETTER BLAST
module load EMBOSS/6.6.0-GCC-10.2.0-Java-13
transeq $RESULTS/Melpomene_gene_inversion.fasta $RESULTS/Protein_Melpomene_gene_inversion.fasta

# MAKE CHETONE BLAST DB
CHETONE_GENOME="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/Chetone_histrio/chetone_histrio_mtDNA_05_02_23.fasta"
module load BLAST+/2.14.1-gompi-2023a
makeblastdb -in $CHETONE_GENOME -dbtype nucl -input_type fasta -out $RESULTS/chetone -title chetone

# GET INVERSION GENE SIZE IN MELPOMENE
cat $INPUTS/$GFF|grep -E "Hmel215003o"|grep -P "\tgene\t"| awk '$4 > 1312092 && $4 < 1796339 {print $9"\t"$5-$4}'|perl -pe 's/ID=//g' > $RESULTS/gene_size.txt

# BLAST EACH INVERSION PROTEIN TO CHETONE + USE gene SIZE IN melpomene
module load Biopython/1.83-foss-2023a


