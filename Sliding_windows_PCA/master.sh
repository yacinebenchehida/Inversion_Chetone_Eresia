module load Biopython/1.79-foss-2022a

REF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/$1"
FASTA_REF=$(ls $REF|grep -E "*.fa$|*.fasta$")

mkdir -p ../Inputs/$1
python -W ignore scaffold_size.py $REF/$FASTA_REF > ../Inputs/$1/$1_scaffold_size.txt


cat ../Inputs/$1/$1_scaffold_size.txt|grep -v -E "MT|loc"|while read line
	do 
		SCAFFOLD=$(echo $line|awk '{print $1}')
		SIZE=$(echo $line|awk '{print $2}')
		sbatch --job-name=PCA_"$SCAFFOLD" ./PCA_sliding.sh $1 $SCAFFOLD $SIZE $2 $3
	done

# $1 = SPECIES ; $2 = window size ; $3 = sliding size
