#!/bin/bash
#SBATCH --job-name=lau_bus       
#SBATCH --ntasks=1                       
#SBATCH --cpus-per-task=1                
#SBATCH --mem=5gb                       
#SBATCH --time=1-00:00:00            
#SBATCH --account=BIOL-SPECGEN-2018 

MAXJOBS=100   # maximum allowed concurrent jobs
SLEEP=30      # how often to re-check
count=0

while read line; do
    count=$((count + 1))
    echo "$count $line"

    # throttle submission until you have fewer than $MAXJOBS running/pending
    while true; do
        # count your running or pending jobs
        njobs=$(squeue -u "$USER" | wc -l)
        # subtract header line
        njobs=$((njobs - 1))

        if [ "$njobs" -lt "$MAXJOBS" ]; then
            break
        fi

        # wait before checking again
        sleep "$SLEEP"
    done

    # submit the job
    sbatch ./1_busco.sh "$line"

    # your existing delay
    sleep 75s

done < 2bedone.txt
