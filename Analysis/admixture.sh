#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=ngsadmix
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# ngsadmix admixture, k=1-5, 20 runs each (Fig 4A, S6; Table S1)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
beagle="$PROJ/PCA/whale.beagle.gz"
OUT="$PROJ/admixture/ADX"
mkdir -p "$OUT"; cd "$OUT"

for K in 1 2 3 4 5; do
  for i in $(seq 1 20); do
    NGSadmix -likes "$beagle" -K "$K" -o "admix_K${K}_run${i}" -P 32
  done
done
