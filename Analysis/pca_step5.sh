#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=pca_persite
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 6-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# per-sample-site pca for east and west (Fig S3)

module --force purge
module load biocontainers
module load angsd
module load pcangsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/PCA/persite"
mkdir -p "$OUT"; cd "$OUT"

for site in east west; do
  angsd -GL 1 -out "$site" -minQ 20 -minMapQ 20 -P 64 \
    -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
    -bam "$PROJ/lists/${site}_bams.txt" -ref "$REF"
  pcangsd -b "${site}.beagle.gz" -o "${site}" -t 64
done
