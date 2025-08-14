#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=10:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=PCA

###########################
# 1 - DEFINE USEFUL PATHS #
###########################
VCF="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/$1/*.vcf.gz"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/Sliding_window_PCA/Results/$1/$6"
PHENO="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/PCA/Inputs/Chetone_histrio/INFORMATION_chetone.txt"
PLINK="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/plink_linux_x86_64_20231018/plink"
SCAFFOLD="$2"
SIZE="$3"
WINDOWS="$4"
SLIDE="$5"

########################   
# 2 - LOAD PLINK and R #
######################## 
module load R/4.2.1-foss-2022a
module load BCFtools/1.15-GCC-11.2.0

mkdir -p $RESULTS

#######################
# 3 - PERFORM THE PCA #
#######################

cat "$PHENO" | grep -E "hydra-Peru|histrio-Columbia|histrio-Peru|hydra-hom|histrio-Ecuador" | cut -f1 > sample_list_"$6".txt
bcftools view -S sample_list_"$6".txt  --regions $SCAFFOLD $VCF -o filtered_output.vcf
bgzip filtered_output.vcf
tabix filtered_output.vcf.gz

VCF="filtered_output.vcf.gz"

echo -e IND"\t"SCAFFOLD"\t"START"\t"END"\t"PC1"\t"PC2"\t"PC3 > $RESULTS/"$SCAFFOLD"_PCA_results.txt
for ((i = 9700000, j = i + $WINDOWS; i < $SIZE && j < $SIZE; i = i + $SLIDE, j=j+$SLIDE)); do 
	#bcftools view -s ^Sample_46-NR14_20_C_h_histrio,Sample_48-NR15-614,Sample_50-NR15-616,Sample_51-NR15-617,Sample_53-NR15-620,Sample_55-NR15-622,Sample_57-NR16-572,Sample_58-NR16-573,Sample_66-M3984,Sample_67-M3985,Sample_70-M3988,Sample_78-NR15_504_C_h_hydra,Sample_79-NR16_16_C_histrio_histrio --regions  "$SCAFFOLD":"$i"-"$j" $VCF > $RESULTS/"$SCAFFOLD"_"$i"_"$j".vcf
	bcftools view --regions "$SCAFFOLD":"$i"-"$j" $VCF > $RESULTS/"$SCAFFOLD"_"$i"_"$j".vcf
	$PLINK --vcf $RESULTS/"$SCAFFOLD"_"$i"_"$j".vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out $RESULTS/"$SCAFFOLD"_"$i"_"$j" --geno 0.08
	cat  $RESULTS/"$SCAFFOLD"_"$i"_"$j".eigenvec  |awk '{print $1"\t"$3"\t"$4"\t"$5}' >> $RESULTS/"$SCAFFOLD"_"$i"_"$j"_TMP_PCA_results.txt
	cat $RESULTS/"$SCAFFOLD"_"$i"_"$j"_TMP_PCA_results.txt| while read line 
		do IND=$(echo $line| awk '{print $1}')
		PC1=$(echo $line| awk '{print $2}')
		PC2=$(echo $line| awk '{print $3}')
		PC3=$(echo $line| awk '{print $4}')
		echo -e $IND"\t"$SCAFFOLD"\t"$i"\t"$j"\t"$PC1"\t"$PC2"\t"$PC3 >> $RESULTS/"$SCAFFOLD"_PCA_results.txt 
		done
	rm $RESULTS/"$SCAFFOLD"_"$i"_"$j"* 
	done
# To handle last bit
j=$SIZE # Do the same for the last window
#bcftools view -s ^Sample_46-NR14_20_C_h_histrio,Sample_48-NR15-614,Sample_50-NR15-616,Sample_51-NR15-617,Sample_53-NR15-620,Sample_55-NR15-622,Sample_57-NR16-572,Sample_58-NR16-573,Sample_66-M3984,Sample_67-M3985,Sample_70-M3988,Sample_78-NR15_504_C_h_hydra,Sample_79-NR16_16_C_histrio_histrio --regions  "$SCAFFOLD":"$i"-"$j" $VCF > $RESULTS/$SCAFFOLD"_"$i"_"$j".vcf
bcftools view --regions "$SCAFFOLD":"$i"-"$j" $VCF > $RESULTS/"$SCAFFOLD"_"$i"_"$j".vcf
$PLINK --vcf $RESULTS/"$SCAFFOLD"_"$i"_"$j".vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --pca --out $RESULTS/"$SCAFFOLD"_"$i"_"$j" --geno 0.1
cat $RESULTS/"$SCAFFOLD"_"$i"_"$j".eigenvec  |awk '{print $1"\t"$3"\t"$4"\t"$5}' >> $RESULTS/"$SCAFFOLD"_"$i"_"$j"_TMP_PCA_results.txt
cat $RESULTS"$SCAFFOLD"_"$i"_"$j"_TMP_PCA_results.txt| while read line 
	do IND=$(echo $line| awk '{print $1}')
	PC1=$(echo $line| awk '{print $2}')
	PC2=$(echo $line| awk '{print $3}')
	PC3=$(echo $line| awk '{print $4}')
	echo -e $IND"\t"$SCAFFOLD"\t"$i"\t"$j"\t"$PC1"\t"$PC2"\t"$PC3 >> $RESULTS/"$SCAFFOLD"_PCA_results.txt 
	 done
rm $RESULTS/"$SCAFFOLD"_"$i"_"$j".*vcf $RESULTS/"$SCAFFOLD"_"$i"_"$j".*
rm filtered_output.vcf.gz

# FOR SLIDING WINDOWS: for ((i = 1, j = $WINDOWS; i < $SIZE && j < $SIZE; i = i + $SLIDE, j=j+$SLIDE)) 
# FOR NON OVERLAPPING WINDOWS: for ((i = 1, j = $WINDOWS; i < $SIZE && j < $SIZE; i = j + 1, j=j+$SLIDE))
# FOR NON OVERLAPPING WINDOWS to start at a non 1 start: for ((i = value, j = value + $WINDOWS; i < $SIZE && j < $SIZE; i = j + 1, j=j+$SLIDE)) 
# sbatch ./PCA_sliding.sh Chetone_histrio ctg001860_1_np1212 12062479 25000 10000 All_25000_10000_geno_0.1

