#!/bin/bash
#SBATCH --job-name=whale_tree
#SBATCH -A johnwayne
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

# module load biocontainers
# module load angsd
# 
# # Set your input BAM file list
# BAM_LIST="/scratch/negishi/allen715/Gray_whales/final_bams/merged/final_merged_use/two_outgroup_bams.list"
# # Set the reference genome FASTA file path
# REF="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
# # Set output directory
# OUT="/scratch/negishi/allen715/Gray_whales/final_tree/final_bams"
# 
# angsd -GL 2 -P 64 -bam "$BAM_LIST" -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -doBcf 1 -SNP_pval 1e-6 \
# -ref "$REF" -doPost 1 -docounts 1 -dogeno 5 -doHWE 1 -minHWEpval 0.05 -remove_bads 1 -only_proper_pairs 1 \
# -minInd 37 -uniqueOnly 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepth 37 -setMaxDepth 2000 -out "$OUT"

#Get list of sites:
# gzip -d final_bams.beagle.gz
# awk '{print $1}' final_bams.beagle | sed 's/NW_/NW./g' | sed 's/_/\t/g' | sed 's/NW./NW_/g' | tail -n +2 > FINAL.SITES
# wc -l final_bams.beagle | awk '{print $1 - 1}' #total sites 2192991
# gzip final_bams.beagle

# module load bcftools
# module load samtools
# module load vcftools
# module load htslib

# # Convert BCF to VCF
#bcftools view -O v -o final_bams.vcf final_bams.bcf
# 
# # Compress the VCF file
#bgzip -c final_bams.vcf > whales.vcf.gz
# 
# # Index the compressed VCF file 
#tabix -p vcf whales.vcf.gz

#module load openjdk/11.0.17_8
#java  -Xss4m -Xmx256g -jar /home/allen715/beagle.27Jan18.7e1.jar gl=whales.vcf.gz out=whales_beagled

# module load anaconda
# # 
# # #vcf2phylip
# python3 vcf2phylip.py --input whales_beagled.vcf.gz  --fasta --nexus --output-prefix  whales_use.not.iupac.resolved

module purge
module load biocontainers
module load iqtree/2.2.2.2

#step 1 TREE with IUPAC not resolved ##use this
#iqtree2 -s  whales_use.not.iupac.resolved.min4.fasta -m GTR+ASC -B 1000 -T 64 -pre step1_ASC.not.iupac.resolved -st DNA
iqtree2 -s whales_use.not.iupac.resolved.min4.fasta -m GTR -B 1000 -T 64 -pre step1.not.iupac.resolved -st DNA

#trying fast bootstrapping for memory
#iqtree2 -s angsdput.not.iupac.resolved.min4.fasta -m GTR+ASC -bb 1000 -T 64 -pre angsdput.not.iupac.resolved -st DNA

#step 2
#iqtree2 -s angsdput.not.iupac.resolved.varsites.phy -m GTR+ASC -B 1000 -T 64 -pre angsdput.not.iupac.resolved.FINAL -st DNA

