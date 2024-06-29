#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=admix_test
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load angsd

reference_genome="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
fai_file="/scratch/negishi/allen715/Gray_whales/reference/ref.fa.fai"
beagle="/scratch/negishi/allen715/Gray_whales/final_admix/whale.beagle.gz" # Path to the Beagle file

cd ADX

# Loop over K values
for K in 1 2 3 4 5; do
  # Loop to run NGSadmix 20 times for each K
  for i in {1..20}; do
    # Construct the output filename
    output_file="admix_K${K}_run${i}"
    # Run NGSadmix
    NGSadmix -likes "$beagle" -K $K -o $output_file -P 64
  done
done

