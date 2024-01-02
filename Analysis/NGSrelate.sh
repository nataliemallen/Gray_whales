#!/bin/bash
#SBATCH -A highmem
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 5:00:00
#SBATCH --job-name=NgsRelate
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load angsd/0.935

# Set your input BAM file list
#bam_list="/scratch/negishi/allen715/Gray_whales/final_bams/merged/merged_bams.list"

# Set the reference genome FASTA file path
#ref_fasta="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"

# Set output directory
#output_dir="/scratch/negishi/allen715/Gray_whales/PCA/"

### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
#angsd -bam "$bam_list" -GL 2 -doMaf 1 -SNP_pval 1e-6 -doMajorMinor 1 -doGlf 3 -minmaf 0.05 -nThreads 64 \
# -out /scratch/negishi/allen715/Gray_whales/kinship/nobo

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
#zcat nobo.mafs.gz | cut -f5 |sed 1d >freq

### run NgsRelate -n is number of individuals
/home/allen715/ngsrelate/ngsRelate/ngsRelate -g nobo.glf.gz -n 71 -f freq -O nobo_relate
