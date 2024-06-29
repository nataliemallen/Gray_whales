#!/bin/bash
#SBATCH --job-name=admix
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH --error=admix_ins.err
#SBATCH --output=admix_ins.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load biocontainers
module load angsd

#beagle for whale admixture 

# Set your input BAM file list
BAM_LIST="/scratch/negishi/allen715/Gray_whales/final_bams/merged/new_merged_bams.list"
# Set the reference genome FASTA file path
REF="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
# Set output directory
OUT1="/scratch/negishi/allen715/Gray_whales/admix/whale"

#changed min depth, need to double check 
#mindepth = 1*37 = 37 maxdepth =25*74 ~= 2000
angsd -bam "$BAM_LIST" -P 64 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 4 -doMaf 2 -SNP_pval 1e-6 \
-doIBS 2 -doCov 1 -doHWE 1 -minHWEpval 0.05 -remove_bads 1 -only_proper_pairs 1 -ref "$REF" \
-minInd 37 -uniqueOnly 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepth 37 -setMaxDepth 2000 -out "$OUT"

#Get list of sites:
# gzip -d whale.beagle.gz
# awk '{print $1}' whale.beagle | sed 's/NW_/NW./g' | sed 's/_/\t/g' | sed 's/NW./NW_/g' | tail -n +2 > FINAL.SITES
# wc -l whale.beagle | awk '{print $1 - 1}' #total sites 3026989
# gzip whale.beagle
