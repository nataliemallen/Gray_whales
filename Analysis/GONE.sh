#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=whale_PLINK
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# *incomplete* GONE script

module load biocontainers
module load angsd
module load bcftools

####STEP1 first prep region and positions files
cut -f1 /scratch/negishi/allen715/Gray_whales/reference/ref.fa.fai > chrs-plink.txt
 
cd /scratch/negishi/allen715/Gray_whales/trees_new/trial2
bcftools query -f '%CHROM\t%POS\n' whale_tree.vcf.gz > /scratch/negishi/allen715/Gray_whales/parentage/angsd.file
# 
#index positions file
angsd sites index /scratch/negishi/allen715/Gray_whales/parentage/angsd.file

###STEP2 run angsd 
# bam list
BAM_LIST="/scratch/negishi/allen715/Gray_whales/final_bams/merged/new_merged_bams.list"
# ref genome path
REF="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
OUT="/scratch/negishi/allen715/Gray_whales/parentage/for_plink"
ANGSDFILE="/scratch/negishi/allen715/Gray_whales/parentage/angsd.file"
 
# plink file for GONE
angsd -bam "$BAM_LIST" -ref "$REF" -rf chrs-plink.txt -sites "$ANGSDFILE" -out "$OUT" \
-doPlink 2 -doGeno -4 -doPost 1 -doMajorMinor 1 -GL 2 -doCounts 1 -doMaf 1 -postCutoff 0.99 -SNP_pval 1e-6 -geno_minDepth 5 \
-minMapQ 30 -minQ 30 -minInd 46 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -baq 2 -P 64

###STEP3 Use PLINK to convert file types
# define the input and output files
module load biocontainers
module load plink/1.90b6.21
# 
TPED_FILE="for_plink.tped"  
TFAM_FILE="for_plink.tfam"  
OUT_PREFIX="plink_output"     #prefix for the .ped and .map files
GONE_EXEC="/scratch/negishi/allen715/chicken_turtles/GONE/Linux/PROGRAMMES/GONE"
OUTPUT_DIR="/scratch/negishi/allen715/chicken_turtles/GONE/"
# 
plink --tped $TPED_FILE --tfam $TFAM_FILE --recode --allow-extra-chr --out $OUT_PREFIX
# 

# convert .ped and .map files to binary format (.bed, .bim, .fam)
plink --file plink_output --make-bed --allow-extra-chr --out plink_output
