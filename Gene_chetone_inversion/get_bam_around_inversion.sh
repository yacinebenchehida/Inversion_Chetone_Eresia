module load SAMtools/1.17-GCC-12.2.0

#cat ref_genome.txt |grep menophilus|awk '{print $1}'|while read line
#do
#        echo $line
#       samtools view -b -o /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam/$line/CAKKCZ020000024_"$line".bam /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam/$line/$line*.bam CAKKCZ020000024.1
#        mv /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam/$line/CAKKCZ020000024_"$line".bam  /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/GWAS/ANGSD/Inputs/Melinaea_menophilus/BAM_files
#done
#SAMPLE="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam/Sample_56-NR15-823/Sample_56-NR15-823_sorted_dedup.bam"
#SAMPLE="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam/Sample_NR15-519_merged_C_histrio/Sample_NR15-519_merged_C_histrio_sorted_dedup.bam"
SAMPLE="/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/1_mapping/Results/1_sorted_dedup_bam/Sample_NR16-17_merged_C_histrio/Sample_NR16-17_merged_C_histrio_sorted_dedup.bam"
samtools view -b $SAMPLE -o Sample_NR16-17_ctg1860_around_inversion.bam ctg001860_1_np1212:9845688-10883586
samtools index Sample_NR16-17_ctg1860_around_inversion.bam