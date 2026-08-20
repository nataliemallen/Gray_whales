#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=fst
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# pairwise fst and dxy between ENP and WNP, downsampled dataset

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/fst"
mkdir -p "$OUT"; cd "$OUT"

declare -A LIST=( [east]="$PROJ/lists/east_bams.txt" [west]="$PROJ/lists/west_ds_bams.txt" )

#1) per-population folded SAF
# for pop in east west; do
#   angsd -bam "${LIST[$pop]}" -ref "$REF" -anc "$REF" -P 64 \
#     -minMapQ 20 -minQ 20 -GL 1 -doSaf 1 -doCounts 1 -out "$pop"
# done
# 
# # 2) 2D SFS prior
# realSFS east.saf.idx west.saf.idx -fold 1 -P 64 > east_west.sfs

# 3) FST index + global estimate
realSFS fst index east.saf.idx west.saf.idx -sfs east_west.sfs -fstout east_west -fold 1 -P 96
realSFS fst stats  east_west.fst.idx -P 96 > east_west_global_fst.txt

# 4) 500 kb sliding-window FST
realSFS fst stats2 east_west.fst.idx -win 500000 -step 500000 -P 96 > east_west_fst_500kb.txt

# 5) per-site alpha/beta (for jackknife + subsampling significance in fst_step2.sh)
realSFS fst print east_west.fst.idx > east_west_fst_persite.txt

# 6) per-population allele frequencies (same SNP set) for DXY
SITES_SNP="$PROJ/PCA/whale.beagle.gz"   # reference SNP set (info only)
for pop in east west; do
  angsd -bam "${LIST[$pop]}" -ref "$REF" -P 96 -minMapQ 20 -minQ 20 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -out "${pop}_maf"
done

cat east_west_global_fst.txt

