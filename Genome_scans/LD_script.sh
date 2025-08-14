#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G 
#SBATCH --account=BIOL-SPECGEN-2018 
#SBATCH --job-name=PrLD

# Input variables
PLINK="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/plink_linux_x86_64_20231018/plink"
FIT_DECAY="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Scripts/Old/ngsLD-1.1.1/scripts/fit_LDdecay.R"
GENOMICS_GENERAL="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/genomics_general"

SUBSET=$1
SAMPLE_LIST=$2
IND=$3
START=$4
END=$5
WD=$6
VCF=$7
WINDOW=$8
CHR=$9
STATS=${10}
IFS=',' read -ra STATS_ARR <<< "$STATS"

mkdir -p $WD
echo "Working on: $CHR:$START-$END"

# Subset and index
bcftools view --regions $CHR:$START-$END -S <(cut -f 1 $SAMPLE_LIST) -e 'F_MISSING > 0.25' $VCF -o $WD/${START}_${END}_filtered_output.vcf
bgzip $WD/${START}_${END}_filtered_output.vcf
tabix $WD/${START}_${END}_filtered_output.vcf.gz

POP1=$(cut -f2 $SAMPLE_LIST | sort -u | sed -n 1p)
POP2=$(cut -f2 $SAMPLE_LIST | sort -u | sed -n 2p)

# Run selected statistics
for stat in "${STATS_ARR[@]}"; do
  echo $stat
  case "$stat" in
    r2)
      echo "Running r2..."
      $PLINK --vcf $WD/${START}_${END}_filtered_output.vcf.gz --r2 --ld-window-kb 20 \
        --out $WD/${START}_${END}_ld_results --double-id --allow-extra-chr \
        --set-missing-var-ids @:# --ld-window 999999 --ld-window-r2 0
      if [ "$(wc -l < $WD/${START}_${END}_ld_results.ld)" -ge 100 ]; then
        awk '$5 - $2 > 200 && $5 - $2 < 5000 {sum += $7} END {print sum/NR}' $WD/${START}_${END}_ld_results.ld |
        paste <(echo -e "$START\t$END") - >> $WD/r2
      fi
      ;;
    r2_decay)
      echo "Running r2_decay..."
      if [ -f $WD/${START}_${END}_ld_results.ld ]; then
        awk -v chr="$CHR" '$5 - $2 > 0 && $5 - $2 < 20000 {print chr ":" $2 "\t" chr ":" $5 "\t" $5 - $2 "\t" $7 "\t" $7 "\t" $7 "\t" $7}' \
        $WD/${START}_${END}_ld_results.ld > $WD/${START}_${END}_ld_decay_input
        Rscript --vanilla --slave $FIT_DECAY --ld_file <(ls $WD/${START}_${END}_ld_decay_input) \
          --fit_level 2 --fit_bin_size 500 --max_kb_dist 20 --n_ind $IND | \
          grep -A 1 LDmax | awk 'NR > 1 {print $4}' |
          paste <(echo -e "$START\t$END") - >> $WD/r2_decay
      fi
      ;;
    Fst|Dxy|pi)
      echo "Running Dxy, Fst, pi..."
      # Make windows file
      echo -e $CHR"\t"$START"\t"$END > $WD/${START}_${END}_windows.txt
      $GENOMICS_GENERAL/VCF_processing/parseVCF.py -i $WD/${START}_${END}_filtered_output.vcf.gz \
        -o $WD/${START}_${END}_filtered_output.geno.gz
      python $GENOMICS_GENERAL/popgenWindows.py \
        --windType predefined --windCoords $WD/${START}_${END}_windows.txt \
        -g $WD/${START}_${END}_filtered_output.geno.gz \
        -o $WD/${START}_${END}_sliding_windows.txt \
        -p $POP1 -p $POP2 -f phased --popsFile $SAMPLE_LIST --addWindowID
      ;;
  esac
done

# Parse sliding window outputs
for stat in "${STATS_ARR[@]}"; do
  case "$stat" in
    Dxy)
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$9}' $WD/${START}_${END}_sliding_windows.txt >> $WD/Dxy
      ;;
    Fst)
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$10}' $WD/${START}_${END}_sliding_windows.txt >> $WD/Fst
      ;;
    pi)
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$7}' $WD/${START}_${END}_sliding_windows.txt >> $WD/pi_${POP1}
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$8}' $WD/${START}_${END}_sliding_windows.txt >> $WD/pi_${POP2}
      ;;
  esac
done

# Clean intermediate files
rm -f $WD/${START}_${END}_*
