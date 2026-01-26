#!/bin/bash
# Author: Yacine Ben Chehida                 # Script author

#SBATCH --time=00:50:00                      # Set job time limit
#SBATCH --nodes=1                            # Request one node
#SBATCH --ntasks=1                           # Request one task
#SBATCH --cpus-per-task=1                    # Request one CPU per task
#SBATCH --mem=2G                             # Set memory requirement
#SBATCH --account=BIOL-SPECGEN-2018          # Set SLURM account
#SBATCH --job-name=PrLD                      # Set SLURM job name

PLINK="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/plink_linux_x86_64_20231018/plink"     # Path to PLINK binary
FIT_DECAY="/mnt/scratch/projects/biol-specgen-2018/yacine/Conv_Evol/LD_Decay/Scripts/Old/ngsLD-1.1.1/scripts/fit_LDdecay.R"   # Path to LD decay fitting script
GENOMICS_GENERAL="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/genomics_general"           # Path to genomics general scripts

SUBSET=$1                                   # Store argument for population subset label
SAMPLE_LIST=$2                              # Store sample list file
IND=$3                                      # Store number of individuals
START=$4                                    # Store window start coordinate
END=$5                                      # Store window end coordinate
WD=$6                                       # Store working directory path
VCF=$7                                      # Store input VCF path
WINDOW=$8                                   # Store window size
CHR=$9                                      # Store chromosome name
STATS=${10}                                 # Store requested statistics as comma separated list

if echo "$STATS" | grep -q "Da"; then       # Check if Da is requested
    if [[ $STATS != *"Dxy"* ]]; then        # If Dxy is not already in the list
        STATS="$STATS,Dxy"                  # Append Dxy to the stats list
    fi                                      # End Dxy check
    if [[ $STATS != *"pi"* ]]; then         # If pi is not already in the list
        STATS="$STATS,pi"                   # Append pi to the stats list
    fi                                      # End pi check
fi                                          # End Da dependency block

IFS=',' read -ra STATS_ARR <<< "$STATS"     # Split statistics list into array entries

mkdir -p "$WD"                              # Create working directory if missing
echo "Working on: $CHR:$START-$END"         # Print progress message

bcftools view --regions "$CHR:$START-$END" -S <(cut -f 1 "$SAMPLE_LIST") -e 'F_MISSING > 0.25' "$VCF" -o "$WD/${START}_${END}_filtered_output.vcf"   # Subset VCF by region and sample list
bgzip "$WD/${START}_${END}_filtered_output.vcf"   # Compress subset VCF
tabix "$WD/${START}_${END}_filtered_output.vcf.gz" # Index compressed VCF

POP1=$(cut -f2 "$SAMPLE_LIST" | sort -u | sed -n 1p)    # Extract first population label
POP2=$(cut -f2 "$SAMPLE_LIST" | sort -u | sed -n 2p)    # Extract second population label

for stat in "${STATS_ARR[@]}"; do              # Loop over requested statistics
  case "$stat" in                               # Switch on statistic name
    r2)                                         # If r2 is requested
      "$PLINK" --vcf "$WD/${START}_${END}_filtered_output.vcf.gz" --r2 --ld-window-kb 20 \
        --out "$WD/${START}_${END}_ld_results" --double-id --allow-extra-chr \
        --set-missing-var-ids @:# --ld-window 999999 --ld-window-r2 0              # Run PLINK r2 computation
      if [ "$(wc -l < "$WD/${START}_${END}_ld_results.ld")" -ge 100 ]; then        # Check minimum number of LD pairs
        awk '$5 - $2 > 200 && $5 - $2 < 5000 {sum += $7} END {print sum/NR}' \
          "$WD/${START}_${END}_ld_results.ld" | paste <(echo -e "$START\t$END") - >> "$WD/r2"   # Compute mean r2 in distance range
      fi
      ;;
    r2_decay)                                    # If r2_decay is requested
      if [ -f "$WD/${START}_${END}_ld_results.ld" ]; then   # Ensure LD file exists
        awk -v chr="$CHR" '$5 - $2 > 0 && $5 - $2 < 20000 {print chr ":" $2 "\t" chr ":" $5 "\t" $5 - $2 "\t" $7 "\t" $7 "\t" $7 "\t" $7}' \
        "$WD/${START}_${END}_ld_results.ld" > "$WD/${START}_${END}_ld_decay_input"     # Format file for LD decay script
        Rscript --vanilla --slave "$FIT_DECAY" --ld_file <(ls "$WD/${START}_${END}_ld_decay_input") \
          --fit_level 2 --fit_bin_size 500 --max_kb_dist 20 --n_ind "$IND" | grep -A 1 LDmax | awk 'NR > 1 {print $4}' \
          | paste <(echo -e "$START\t$END") - >> "$WD/r2_decay"                     # Run LD decay and extract LDmax
      fi
      ;;
    Fst|Dxy|pi)                                  # If Fst, Dxy, or pi is requested
      echo -e "$CHR\t$START\t$END" > "$WD/${START}_${END}_windows.txt"             # Create predefined window file
      "$GENOMICS_GENERAL/VCF_processing/parseVCF.py" -i "$WD/${START}_${END}_filtered_output.vcf.gz" \
        -o "$WD/${START}_${END}_filtered_output.geno.gz"                           # Convert VCF to geno format
      python "$GENOMICS_GENERAL/popgenWindows.py" \
        --windType predefined --windCoords "$WD/${START}_${END}_windows.txt" \
        -g "$WD/${START}_${END}_filtered_output.geno.gz" \
        -o "$WD/${START}_${END}_sliding_windows.txt" \
        -p "$POP1" -p "$POP2" -f phased --popsFile "$SAMPLE_LIST" --addWindowID    # Compute Dxy Fst pi per window
      ;;
    Tajima)                                     # If Tajima is requested
      module load VCFtools/0.1.16-GCC-11.2.0    # Load VCFtools module
      for POP in $(cut -f2 "$SAMPLE_LIST" | sort -u); do     # Loop over populations
        bcftools view --regions "$CHR:$START-$END" -S <(grep "$POP" "$SAMPLE_LIST" | cut -f1) \
          "$VCF" -o "$WD/${START}_${END}_filtered_output_${POP}.vcf"               # Subset VCF for each population
        bgzip "$WD/${START}_${END}_filtered_output_${POP}.vcf"                     # Compress population VCF
        vcftools --gzvcf "$WD/${START}_${END}_filtered_output_${POP}.vcf.gz" --out "$WD/${START}_${END}_${POP}" --TajimaD 10000000000000   # Compute Tajima D
        paste <(echo -e "$START\t$END") <(awk 'NR==2 {print $4}' "$WD/${START}_${END}_${POP}.Tajima.D") >> "$WD/Tajima_${POP}"             # Store Tajima output
      done
      ;;
    Da)                                          # If Da is requested, no direct computation here
      :                                          # Da is computed later after Dxy and pi are parsed
      ;;
  esac                                           # End case on statistic name
done                                             # End loop over statistics

for stat in "${STATS_ARR[@]}"; do                # Parse sliding window outputs
  case "$stat" in                                 # Switch on statistic name
    Dxy)                                          # Parse Dxy
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$9*$6/($4-$3)}' \
        "$WD/${START}_${END}_sliding_windows.txt" >> "$WD/Dxy"                     # Write Dxy file
      ;;
    Fst)                                          # Parse Fst
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$10}' \
        "$WD/${START}_${END}_sliding_windows.txt" >> "$WD/Fst"                    # Write Fst file
      ;;
    pi)                                           # Parse pi for both populations
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$7*$6/($4-$3)}' \
        "$WD/${START}_${END}_sliding_windows.txt" >> "$WD/pi_${POP1}"             # Write pi for population one
      awk -F"," 'NR > 1 {print $3"\t"$4"\t"$8*$6/($4-$3)}' \
        "$WD/${START}_${END}_sliding_windows.txt" >> "$WD/pi_${POP2}"             # Write pi for population two
      ;;
  esac                                           # End case for parsing
done                                             # End loop over parsed statistics

if echo "$STATS" | grep -q "Da"; then            # Check if Da should be generated
  if [ -f "$WD/Dxy" ] && [ -f "$WD/pi_${POP1}" ] && [ -f "$WD/pi_${POP2}" ]; then  # Ensure required inputs exist
    paste "$WD/Dxy" "$WD/pi_${POP1}" "$WD/pi_${POP2}" |                           # Combine Dxy and both pi files
    awk '{print $1"\t"$2"\t"$3 - (($6 + $9)/2)}' > "$WD/Da"                       # Compute Da and write output
  fi                                           # End input existence check
fi                                             # End Da generation block

rm -f "$WD/${START}_${END}_*gz*"                  # Remove temporary files
