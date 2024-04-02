#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=beagle
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

### generating Beagle

module load biocontainers
module load angsd
 
cd /scratch/negishi/allen715/Gray_whales/PCA/

# Set your input BAM file list
bam_list="/scratch/negishi/allen715/Gray_whales/final_bams/merged/merged_bams.list"

# Set the reference genome FASTA file path
ref_fasta="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"

# Set output directory
output_dir="/scratch/negishi/allen715/Gray_whales/PCA/"

#generate beagle file
angsd -GL 1 -out "$output_dir/whole_genome" -minQ 20 -P 10 \
-doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
-bam "$bam_list" -ref "$ref_fasta"
