#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=fst_extreme
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# genome-wide fst between low-admixture ENP and WNP individuals (Fig 3)

module --force purge
module load biocontainers
module load angsd
 
PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/fst/extreme"
mkdir -p "$OUT"; cd "$OUT"
 
declare -A LIST=( [enp6]="$PROJ/lists/enp6_bams.txt" [wnp6]="$PROJ/lists/wnp6_extreme_bams.txt" )
 
# for pop in enp6 wnp6; do
#   angsd -bam "${LIST[$pop]}" -ref "$REF" -anc "$REF" -P 64 \
#     -minMapQ 20 -minQ 20 -GL 1 -doSaf 1 -doCounts 1 -out "$pop"
# done
 
realSFS enp6.saf.idx wnp6.saf.idx -fold 1 -P 64 > enp6_wnp6.sfs
realSFS fst index enp6.saf.idx wnp6.saf.idx -sfs enp6_wnp6.sfs -fstout enp6_wnp6 -fold 1 -P 64
realSFS fst stats enp6_wnp6.fst.idx -P 64 > enp6_wnp6_global_fst.txt
cat enp6_wnp6_global_fst.txt
 
