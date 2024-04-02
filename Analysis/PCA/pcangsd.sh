#!/bin/bash
#SBATCH --job-name=cov_matrix
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH --error=cov.err
#SBATCH --output=cov.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load biocontainers
module load angsd
module load pcangsd
 
cd /scratch/bell/allen715/Chicken_turtles/chrysemys_BEAGLE2/cov/

reference_genome="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
fai_file="/scratch/negishi/allen715/Gray_whales/reference/ref.fa.fai"
beagle_file="/scratch/negishi/allen715/Gray_whales/PCA/whole_genome.beagle.gz" # Path to the Beagle file

pcangsd -b $beagle_file -o whales -t 40 
