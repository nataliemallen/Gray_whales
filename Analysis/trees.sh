#!/bin/bash
#SBATCH -A highmem
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 1-00:00:00
#SBATCH --job-name=tree_trial
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module purge
#module load biocontainers
#module load angsd

######
#step 1 beagle
######

# Set your input BAM file list
# BAM_LIST="/scratch/negishi/allen715/Gray_whales/final_bams/merged/new_merged_bams.list"
# # Set the reference genome FASTA file path
# REF="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
# # Set output directory
# OUT="/scratch/negishi/allen715/Gray_whales/admixture/"
# 
# # angsd -GL 2 \
# #  -P 64 \
# #  -bam "$BAM_LIST" \
# #  -doGlf 2 \
# #  -doMajorMinor 1 \
# #  -doMaf 1 \
# #  -ref "$REF" \
# #  -doBcf 1 \
# #  -doPost 1 \
# #  -docounts 1 \
# #  -dogeno 5
# 
# module load bcftools
# module load samtools
# module load vcftools
# module load htslib
# 
# # Convert BCF to VCF
# bcftools view -O v -o angsdput.vcf angsdput.bcf
# 
# # Compress the VCF file
# bgzip -c angsdput.vcf > angsdput.vcf.gz
# 
# # Index the compressed VCF file (optional)
# tabix -p vcf angsdput.vcf.gz

######
#step 2 java
######
module load openjdk/11.0.17_8
java  -Xss4m -Xmx60g -jar /home/allen715/beagle.27Jan18.7e1.jar gl=angsdput.vcf.gz out=whales_use

#used this
#python3 vcf2phylip.py --input angsdput.vcf.gz  --fasta --nexus --output-prefix  angsdput.not.iupac.resolved

# module purge
# module load biocontainers
# module load iqtree/2.2.2.2

#step 1 TREE with IUPAC not resolved ##use this
#iqtree2 -s  angsdput.not.iupac.resolved.min4.fasta -m GTR+ASC -B 1000 -T 64 -pre angsdput.not.iupac.resolved -st DNA

#trying fast bootstrapping for memory
#iqtree2 -s angsdput.not.iupac.resolved.min4.fasta -m GTR+ASC -bb 1000 -T 64 -pre angsdput.not.iupac.resolved -st DNA

#step 2
#iqtree2 -s angsdput.not.iupac.resolved.varsites.phy -m GTR+ASC -B 1000 -T 64 -pre angsdput.not.iupac.resolved.FINAL -st DNA
