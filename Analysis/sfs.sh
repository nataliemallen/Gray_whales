#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=sfs_tajima
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# folded sfs and tajima's d per sample site, downsampled dataset (Fig 5C, S10)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/sfs"
mkdir -p "$OUT"; cd "$OUT"

declare -A LIST=( [east]="$PROJ/lists/east_bams.txt" [west]="$PROJ/lists/west_ds_bams.txt" )

for pop in east west; do
 # 1) site allele frequency likelihoods (reference as ancestral)
  angsd -bam "${LIST[$pop]}" -ref "$REF" -anc "$REF" -P 64 \
    -minMapQ 20 -minQ 20 -GL 1 -doSaf 1 -doCounts 1 -out "$pop"

  #2) folded SFS
  realSFS "$pop.saf.idx" -fold 1 -P 64 > "$pop.sfs"

  # 3) thetas
  realSFS saf2theta "$pop.saf.idx" -sfs "$pop.sfs" -fold 1 -outname "$pop" -P 64

  # 4) genome-wide + sliding-window (50 kb / 10 kb) Tajima's D and thetas
  thetaStat do_stat "$pop.thetas.idx" > "${pop}_global_theta.txt"
  thetaStat do_stat "$pop.thetas.idx" -win 50000 -step 10000 -outnames "${pop}_50kb"
done
