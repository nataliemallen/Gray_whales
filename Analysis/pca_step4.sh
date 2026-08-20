#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=pca_norel
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# relatives-excluded pca (Fig S2)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/PCA/norel"
mkdir -p "$OUT"; cd "$OUT"

# Build a bam list excluding flagged relatives
grep -vFf "$PROJ/lists/relatives_remove.txt" "$PROJ/bams_list.txt" > norel_bams.txt
echo "Individuals after removing relatives: $(wc -l < norel_bams.txt)"

angsd -GL 1 -out whale_norel -minQ 20 -minMapQ 20 -P 64 \
  -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
  -bam norel_bams.txt -ref "$REF"

module load pcangsd
pcangsd -b whale_norel.beagle.gz -o whale_norel -t 64
