#!/bin/bash
# Author: Yacine Ben Chehida

###################################
# Set usage and script parameters #
###################################
# set usage information
usage()
{
cat << EOF
./master.sh --pattern <subset_pattern> -v <vcf> -f <pop_file> -c <chromosome> -s <start_position> -e <end_position> --w <window> --sl <slide> -p <prefix> --is <inversion_start> --ie <inversion_end> -o <path_output>  -h

OPTIONS:
  --pattern        pattern to look for to determine the species to analyse
  -v		       vcf 
  -f		       population group file
  -c               chromosome name
  -s		       start position
  -e	           end position
  -w               window size to analyse
  --sl		       sliding window size
  -p		       output prefix
  --is		       start position for the inversion
  --ie		       end position for the inversion
  --st             statistics to compute (Fst,Dxy,r2,r2_decay,pi)
  -o               output path
  -h               usage information and help (this message)
EOF
exit 0
}

# Parse command-line options
options=$(getopt -o v:f:c:s:e:w:p:o:h --long pattern:,sl:,ie:,is:,st: -- "$@")

# Check for getopt errors
if [ $? -ne 0 ]; then
  usage
fi

# Set the parsed options back to positional parameters
eval set -- "$options"

# Extract options and their arguments
while [[ $# -gt 0 ]]
do
    case $1 in
		-v)
		  VCF=$2
		  shift 2
		  ;;
		-c)
		  CHR=$2
		  shift 2
		  ;;
		-s)
		  START=$2
		  shift 2
		  ;;
		-e)
		  END=$2
		  shift 2
		  ;;
		-f)
		  POP=$2
		  shift 2
		  ;;
		-w)
		  WINDOWS=$2
		  shift 2
		  ;;
        --sl)
		  SLIDE=$2
		  shift 2
		  ;;
		--pattern)
          PATTERN=$2
          shift 2
           ;;
		-o)
		  OUTPUT_PATH=$2
		  shift 2
		  ;;
		-p)
		  PREFIX=$2
		  shift 2
		  ;;
        --ie)
          INV_START=$2
          shift 2
           ;;
        --is)
          INV_END=$2
          shift 2
           ;;
		--st)
          STATS=$2
          shift 2
           ;;
		-h)
		  usage
		  ;;
		--)
		  shift
		  break
		  ;;
		*)
		  usage
		  ;;
	esac
done

# Check if all mandatory arguments are provided and print which one are missing
if [[ -z $VCF ]] || [[ -z $CHR ]] || [[ -z $START ]] || [[ -z $END ]] || [[ -z $WINDOWS ]] || [[ -z $POP ]] || [[ -z $SLIDE ]] || [[ -z $OUTPUT_PATH ]] || [[ -z $INV_START ]] || [[ -z $INV_END ]] || [[ -z $PATTERN ]] || [[ -z $PREFIX ]] ||  [[ -z $STATS ]]; then
    missing=()

    for var in VCF CHR START END WINDOWS POP SLIDE OUTPUT_PATH INV_END INV_START PATTERN PREFIX STATS; do
        if [[ -z "${!var}" ]]; then
            missing+=("$var")
        fi
    done

    if (( ${#missing[@]} > 0 )); then
        echo "Error: Missing required arguments: ${missing[*]}"
        usage
        exit 1
    fi
fi
# Display arguments value
echo -e VCF file:"${VCF}"
echo -e Population file:"${POP}"
echo -e Population pattern to know which population to analyse:"${PATTERN}"
echo -e Chromosome:"${CHR}"
echo -e Start:"${START}"
echo -e End:"${END}"
echo -e Winwos size:"${WINDOWS}"
echo -e Sliding windows size:"${SLIDING}"
echo -e Inversion start:"${INV_START}"
echo -e Inversion end:"${INV_END}"
echo -e Computed statistics:"${STATS}"
echo -e Output path:"${OUTPUT_PATH}"
echo -e Prefix:"${PREFIX}"

# Load modules 
module load R/4.2.1-foss-2022a
module load BCFtools/1.19-GCC-13.2.0
module load VCFtools/0.1.16-GCC-11.2.0

# Set samples to analyse
cat "$POP" | grep -E $PATTERN > sample_list.txt
IND=$(cat sample_list.txt|wc -l)

# Create working directory
mkdir -p $OUTPUT_PATH/$PREFIX

# Run LD decay and R2
for ((i = $START, j = i + $WINDOWS; i < $END && j < $END; i = i + $SLIDE, j=j+$SLIDE)); do
    sbatch ./LD_script.sh $PATTERN sample_list.txt $IND $i $j $OUTPUT_PATH/$PREFIX $VCF $WINDOWS $CHR $STATS
done
# RUN LAST WINDOWS
j=$END
sbatch ./LD_script.sh $PATTERN sample_list.txt $IND $i $j $OUTPUT_PATH/$PREFIX $VCF $WINDOWS $CHR $STATS

# Plot
running_jobs_alignments=$(squeue|grep $(whoami)| grep -P "PrLD"| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
sbatch --job-name=plot --dependency=aftercorr:$running_jobs_alignments ./plot.sh $OUTPUT_PATH $PREFIX $INV_START $INV_END $STATS


#./master_option.sh  --pattern "Col|histrio-Peru" -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/PCA/Inputs/Chetone_histrio/INFORMATION_chetone.txt -c ctg001860_1_np1212 -s 7000000 -e 12000000 -w 100000 --sl 100000 --is 9855688 --ie 10873586 -o ./ -p sixiemesens -v /mnt/scratch/projects/biol-specgen-2018/yacine/Bioinformatics/6_Combine_intervals/Results/Chetone_histrio/Chetone_histrio.vcf.gz --st Fst,Dxy,r2,r2_decay,pi
#./master_option.sh  --pattern "hydra|Hydra" -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/grouping_phased.txt -c ctg001860_1_np1212 -s 9835233 -e 10873574 -w 50000 --sl 10000 --ie 9855688 --is 10873586 -o ./ -p sixiemesens -v /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/Inversion_Chetone_histrio_diploid_haplotypes.vcf.gz --st Fst,Dxy,r2,r2_decay,pi
#./master_option.sh  --pattern "Hydra_haplotype" -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/grouping_phased.txt -c ctg001860_1_np1212 -s 9835233 -e 10873574 -w 50000 --sl 5000 --ie 9855688 --is 10873586 -o . -p inverted_hydra_haplotype -v /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/Inversion_Chetone_histrio_diploid_haplotypes.vcf.gz --st r2,r2_decay
#./master_option.sh  --pattern "Histrio_haplotype_Peru_het_hydra" -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/grouping_phased.txt -c ctg001860_1_np1212 -s 9835233 -e 10873574 -w 50000 --sl 5000 --ie 9855688 --is 10873586 -o . -p histrio_haplotype_Peru_from_hetero_hydra -v /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/Inversion_Chetone_histrio_diploid_haplotypes.vcf.gz --st r2,r2_decay
#./master_option.sh  --pattern "Histrio_haplotype_Columbia|Histrio_haplotype_Peru" -f /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/grouping_phased.txt -c ctg001860_1_np1212 -s 9835233 -e 10873574 -w 50000 --sl 5000 --ie 9855688 --is 10873586 -o . -p histrio_peru_colombia -v /mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Inputs/Inversion_Chetone_histrio_diploid_haplotypes.vcf.gz --st r2,r2_decay

