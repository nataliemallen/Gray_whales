#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=phylo_geno
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# angsd genotype calls (min depth 10x) for the nuclear phylogeny bcf (Fig 4A)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/phylogeny"
mkdir -p "$OUT"; cd "$OUT"

# cohort + outgroup bam list
cat "$PROJ/bams_list.txt" > tree_bams.list
echo "$PROJ/reference/outgroup/finwhale.bam" >> tree_bams.list
N=$(wc -l < tree_bams.list)

angsd -GL 2 -P 64 -bam tree_bams.list -ref "$REF" \
  -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 \
  -doBcf 1 -doPost 1 -doGeno 5 -doCounts 1 -geno_minDepth 10 \
  -doHWE 1 -minHWEpval 0.05 -remove_bads 1 -only_proper_pairs 1 \
  -uniqueOnly 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
  -out tree_bams
