#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=pca
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# generate the genome-wide beagle genotype likelihoods

module --force purge
module load biocontainers
module load angsd
 
cd /scratch/gautschi/allen715/2026_whales/

# input BAM file list
bam_list="/scratch/gautschi/allen715/2026_whales/bams_list.txt"

# reference genome
ref_fasta="/scratch/gautschi/allen715/2026_whales/reference/ref.fa"

# output directory
output_dir="/scratch/gautschi/allen715/2026_whales/PCA"

# generate beagle file
angsd -GL 1 -out "$output_dir/whale" -minQ 20 -P 10 \
-doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
-bam "$bam_list" -ref "$ref_fasta"