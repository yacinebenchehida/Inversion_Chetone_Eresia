#!/bin/bash
# Author: Yacine Ben Chehida

###################################
# Set usage and script parameters #
###################################
# set usage information
usage()
{
cat << EOF
./master.sh --pattern <subset_pattern> -v <vcf> -f <pop_file> -c <chromosome> -s <start_position> -e <end_position> -p <prefix> --is <inversion_start> --ie <inversion_end> -o <path_output> --window_mode <bp|snp> [window options] --st <statistics> -h

OPTIONS:
  --pattern        pattern to look for to determine the species to analyse
  -v               vcf
  -f               population group file
  -c               chromosome name
  -s               start position
  -e               end position
  -p               output prefix
  --is             start position for the inversion
  --ie             end position for the inversion
  --st             statistics to compute (Fst,Dxy,Da,r2,r2_decay,pi,Tajima)
  --window_mode    bp or snp
  -o               output path
  -h               usage information and help (this message)

WINDOW OPTIONS (bp mode):
  -w               window size to analyse (base pairs)
  --sl             sliding window size (base pairs)

WINDOW OPTIONS (snp mode):
  --snp_window     number of SNPs per window
  --sl             sliding window size (number of SNPs)
EOF
exit 0
}

####################################
# Parse command-line options        #
####################################
options=$(getopt -o v:f:c:s:e:w:p:o:h --long pattern:,sl:,ie:,is:,st:,window_mode:,snp_window: -- "$@")

# Check for getopt errors
if [ $? -ne 0 ]; then
  usage
fi

# Set the parsed options back to positional parameters
eval set -- "$options"

# Extract options and their arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v)
          VCF=$2
          shift 2
          ;;
        -f)
          POP=$2
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
        -w)
          WINDOWS=$2
          shift 2
          ;;
        --sl)
          SLIDE=$2
          shift 2
          ;;
        -p)
          PREFIX=$2
          shift 2
          ;;
        -o)
          OUTPUT_PATH=$2
          shift 2
          ;;
        --pattern)
          PATTERN=$2
          shift 2
          ;;
        --is)
          INV_START=$2
          shift 2
          ;;
        --ie)
          INV_END=$2
          shift 2
          ;;
        --st)
          STATS=$2
          shift 2
          ;;
        --window_mode)
          WINDOW_MODE=$2
          shift 2
          ;;
        --snp_window)
          SNP_WINDOW=$2
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

####################################
# Set defaults and sanity checks    #
####################################
# Set default window mode
if [[ -z "$WINDOW_MODE" ]]; then
  WINDOW_MODE="bp"
fi

# Check mandatory arguments (mode-independent)
missing=()
for var in VCF CHR START END POP OUTPUT_PATH INV_START INV_END PATTERN PREFIX STATS; do
    if [[ -z "${!var}" ]]; then
        missing+=("$var")
    fi
done

# Check mode-specific required arguments
if [[ "$WINDOW_MODE" == "bp" ]]; then
    if [[ -z "$WINDOWS" ]]; then
        missing+=("WINDOWS (-w)")
    fi
    if [[ -z "$SLIDE" ]]; then
        missing+=("SLIDE (--sl)")
    fi
fi

if [[ "$WINDOW_MODE" == "snp" ]]; then
    if [[ -z "$SNP_WINDOW" ]]; then
        missing+=("SNP_WINDOW (--snp_window)")
    fi
    if [[ -z "$SLIDE" ]]; then
        missing+=("SLIDE (--sl)")
    fi
fi

if (( ${#missing[@]} > 0 )); then
    echo "Error: Missing required arguments: ${missing[*]}"
    usage
    exit 1
fi

####################################
# Display run configuration         #
####################################
echo -e "VCF file: ${VCF}"
echo -e "Population file: ${POP}"
echo -e "Population pattern to know which population to analyse: ${PATTERN}"
echo -e "Chromosome: ${CHR}"
echo -e "Start: ${START}"
echo -e "End: ${END}"
echo -e "Window mode: ${WINDOW_MODE}"

if [[ "$WINDOW_MODE" == "bp" ]]; then
  echo -e "Window size (bp): ${WINDOWS}"
  echo -e "Sliding window size (bp): ${SLIDE}"
fi

if [[ "$WINDOW_MODE" == "snp" ]]; then
  echo -e "Window size (SNPs): ${SNP_WINDOW}"
  echo -e "Sliding window size (SNPs): ${SLIDE}"
fi

echo -e "Inversion start: ${INV_START}"
echo -e "Inversion end: ${INV_END}"
echo -e "Computed statistics: ${STATS}"
echo -e "Output path: ${OUTPUT_PATH}"
echo -e "Prefix: ${PREFIX}"

####################################
# Load required modules             #
####################################
module load R/4.2.1-foss-2022a
module load BCFtools/1.19-GCC-13.2.0
module load VCFtools/0.1.16-GCC-11.2.0

####################################
# Prepare sample list               #
####################################
cat "$POP" | grep -P "$PATTERN$" > sample_list.txt
IND=$(wc -l < sample_list.txt)

####################################
# Create working directory          #
####################################
mkdir -p "$OUTPUT_PATH/$PREFIX"

####################################
# Submit LD and statistics jobs     #
####################################
if [[ "$WINDOW_MODE" == "bp" ]]; then
  for ((i = START, j = i + WINDOWS; i < END && j < END; i = i + SLIDE, j = j + SLIDE)); do
      sbatch ./LD_script.sh "$PATTERN" sample_list.txt "$IND" "$i" "$j" "$OUTPUT_PATH/$PREFIX" "$VCF" "$WINDOWS" "$CHR" "$STATS"
  done
  j=$END
  sbatch ./LD_script.sh "$PATTERN" sample_list.txt "$IND" "$i" "$j" "$OUTPUT_PATH/$PREFIX" "$VCF" "$WINDOWS" "$CHR" "$STATS"
fi

if [[ "$WINDOW_MODE" == "snp" ]]; then
  bcftools view --regions "$CHR:$START-$END" -S <(cut -f 1 sample_list.txt) "$VCF" \
  | bcftools query -f '%POS\n' \
  | awk -v n="$SNP_WINDOW" '
      NR==1 {s=$1}
      NR%n==0 {print s"\t"$1; s=""}
      NR%n==1 {s=$1}
      END {if(s!="") print s"\t"$1}
    ' | while read i j; do
        sbatch ./LD_script.sh "$PATTERN" sample_list.txt "$IND" "$i" "$j" "$OUTPUT_PATH/$PREFIX" "$VCF" "$WINDOWS" "$CHR" "$STATS"
    done
fi

####################################
# Submit plotting job               #
####################################
running_jobs_alignments=$(squeue | grep "$(whoami)" | grep -P "PrLD" | awk '{print $1}' | perl -pe 's/\n/,/g' | sed 's/,$//g')
sbatch --job-name=plot --dependency=aftercorr:$running_jobs_alignments ./plot.sh "$OUTPUT_PATH" "$PREFIX" "$INV_START" "$INV_END" "$STATS"
