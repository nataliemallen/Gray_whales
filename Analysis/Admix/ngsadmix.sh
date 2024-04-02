#!/bin/bash
#SBATCH --job-name=admix
#SBATCH -A johnwayne
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH --error=admix_ins.err
#SBATCH --output=admix_ins.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load biocontainers
module load angsd

reference_genome="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
fai_file="/scratch/negishi/allen715/Gray_whales/reference/ref.fa.fai"
##used beagle from PCA
beagle="/scratch/negishi/allen715/Gray_whales/PCA/whole_genome.beagle.gz" # Path to the Beagle file

#with beagle file, run ngsadmix for k=2
#output_base_ngsadmix.qopt will contain 1 row for each individual and one column for each k 
mkdir k2
cd k2

NGSadmix -likes "$beagle" -K 2 -o output_base_ngsadmix -P 64

cd ..
mkdir k3
cd k3
NGSadmix -likes "$beagle" -K 3 -o output_base_ngsadmix -P 64

cd ..
mkdir k4
cd k4
NGSadmix -likes "$beagle" -K 4 -o output_base_ngsadmix -P 64
