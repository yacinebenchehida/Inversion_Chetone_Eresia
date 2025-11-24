!/bin/bash
#SBATCH --job-name=lau_bus       
#SBATCH --partition=nodes
#SBATCH --ntasks=1                       
#SBATCH --cpus-per-task=1                
#SBATCH --mem=1gb                       
#SBATCH --time=01:00:00                  
#SBATCH --account=BIOL-SPECGEN-2018 

count=0

cat list_genome.txt | while read line; do
    count=$((count + 1))
    echo "$count $line"
    sbatch ./busco.sh "$line"
    sleep 30s
done
