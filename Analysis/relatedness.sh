#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=ngsrelate
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 5-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# pairwise relatedness with ngsrelate v2 (Fig 3, S1)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/relatedness"
mkdir -p "$OUT"; cd "$OUT"

bam_list="$PROJ/bams_list.txt"
N_IND=71

# Individual ID labels (basename, in bam-list order) for NgsRelate output
sed 's#.*/##; s#\.bam$##' "$bam_list" > ids.txt

# 1) genotype likelihoods (glf3) + allele frequencies
angsd -bam "$bam_list" -ref "$REF" -GL 2 -doMajorMinor 1 -doMaf 1 \
  -SNP_pval 1e-6 -minMaf 0.05 -minMapQ 20 -minQ 20 -doGlf 3 \
  -nThreads 64 -out whale

# 2) frequency column (no header) for NgsRelate
zcat whale.mafs.gz | cut -f6 | sed 1d > freq

module load anaconda
conda activate ngsrelate_env

# 3) NgsRelate (path to your local install)
ngsRelate -g whale.glf.gz -n "$N_IND" -f freq -z ids.txt -O whale_relate.tsv

