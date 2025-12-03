GENES_LIST="genes_order_busco.txt"             # file containing busco genes
CHUNK_SIZE=5                               # number of genes per job
job_ids=()                                  # array storing job IDs

mapfile -t GENES < "$GENES_LIST"        # read genes into array
TOTAL=${#GENES[@]}                        # total genes count

for ((i=0; i<$TOTAL; i+=CHUNK_SIZE)); do    # loop in chunks
    chunk=("${GENES[@]:i:CHUNK_SIZE}")    # extract chunk of genes
    sbatch_output=$(sbatch test_if_gene_in_inversion.sh "${chunk[@]}")  # submit job
    jid=$(echo "$sbatch_output" | awk '{print $4}')               # capture job ID
    job_ids+=("$jid")                                              # store job ID
done

dep_str=""                                  # dependency string
for j in "${job_ids[@]}"; do
    dep_str="$dep_str:$j"                   # append each job ID
done

sbatch --dependency=afterok${dep_str} collect_and_plot_histo.R   # submit collector
