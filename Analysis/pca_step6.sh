#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=pca_perchr
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 6-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# per-chromosome pca (Fig S4)

module --force purge
module load biocontainers
module load angsd
module load pcangsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/PCA/perchr"
mkdir -p "$OUT"; cd "$OUT"

while read -r chr; do
  [ -z "$chr" ] && continue
  angsd -GL 1 -out "chr_${chr}" -r "${chr}:" -minQ 20 -minMapQ 20 -P 64 \
    -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
    -bam "$PROJ/bams_list.txt" -ref "$REF"
  # only run PCAngsd if SNPs were found
  if [ "$(zcat chr_${chr}.beagle.gz | wc -l)" -gt 1 ]; then
    pcangsd -b "chr_${chr}.beagle.gz" -o "chr_${chr}" -t 64
  fi
done < "$PROJ/reference/chrs.txt"
