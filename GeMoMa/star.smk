GENOME = "/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/0_Data/reference_genomes/Chetone_histrio/chetone_histrio_mtDNA_05_02_23.fasta"
READS_R1 = "/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/Data/Chetone_histrio/RNA_seq_Chetone_R1.fastq"
READS_R2 = "/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/Data/Chetone_histrio/RNA_seq_Chetone_R2.fastq"
THREADS = 4
STAR_INDEX = "/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/GeMoMa/Results/star/index"
ALIGN_PREFIX = "/mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/7_genome_annotation/GeMoMa/Results/star/"

rule all:
    input:
        f"{ALIGN_PREFIX}aln.sorted.bam"

rule star_index:
    input:
        genome=GENOME
    output:
        directory(STAR_INDEX)
    threads: THREADS
    shell:
        """
        module load STAR/2.7.10b-GCC-11.3.0
        STAR --runThreadN {threads} --runMode genomeGenerate \
             --genomeDir {output} --genomeFastaFiles {input.genome}
        """

rule star_align:
    input:
        index=STAR_INDEX,
        r1=READS_R1,
        r2=READS_R2
    output:
        bam=f"{ALIGN_PREFIX}aln.sorted.bam"
    threads: THREADS
    shell:
        """
        module load STAR/2.7.10b-GCC-11.3.0
        module load SAMtools/1.20-GCC-13.2.0

        STAR --runThreadN {threads} --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {ALIGN_PREFIX}

        mv {ALIGN_PREFIX}Aligned.sortedByCoord.out.bam {output.bam}
        samtools index {output.bam}
        """