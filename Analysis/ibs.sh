#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=ibs
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 6-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# ibs distance matrix and neighbour-joining tree (Fig S8, 5D)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/ibs"
mkdir -p "$OUT"; cd "$OUT"

bam_list="$PROJ/bams_list.txt"
sed 's#.*/##; s#\.bam$##' "$bam_list" > ids.txt

# single-read sampling IBS -> pairwise distance matrix
angsd -bam "$bam_list" -ref "$REF" -P 64 -minMapQ 20 -minQ 20 \
  -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
  -doIBS 1 -makeMatrix 1 -doCounts 1 -out whale_ibs

# build a PHYLIP distance matrix (taxa count + labels) for FastME
N=71
{ echo "  $N"; paste -d' ' ids.txt whale_ibs.ibsMat; } > whale_ibs.phy

module load anaconda
conda activate fastme 

fastme -i whale_ibs.phy -o whale_ibs_nj.nwk -m N   # N = NJ
echo "NJ tree: $OUT/whale_ibs_nj.nwk"
