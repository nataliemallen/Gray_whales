#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=whale_rohs_bcf
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

### identifying ROHs with plink

module load biocontainers
module load plink/1.90b6.21
# 
cd /scratch/negishi/allen715/Gray_whales/parentage/
# 
# run PLINK to calculate ROHs
plink --bfile plink_output --homozyg --homozyg-snp 50 --homozyg-kb 100 --homozyg-density 50 --allow-extra-chr --out ./ROHs/roh_output

### identifying ROHs with bcftools

module load biocontainers
module load bcftools
module load htslib

cp /scratch/negishi/allen715/Gray_whales/trees_new/trial2/whale_tree.vcf /scratch/negishi/allen715/Gray_whales/parentage/ROHs/
# 
bcftools sort -o whale_tree.vcf sorted_whale_tree.vcf
bgzip -@ 20 -c sorted_whale_tree.vcf > sorted_whale_tree.vcf.gz
tabix -p vcf --csi sorted_whale_tree.vcf.gz
# 
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' sorted_whale_tree.vcf.gz | bgzip -c > AF.tab.gz
tabix -s1 -b2 -e2 --csi AF.tab.gz

bcftools roh -G 30 --AF-file AF.tab.gz -o bcf_roh_out.txt sorted_whale_tree.vcf.gz

#####

cp /scratch/negishi/allen715/Gray_whales/trees_new/trial2/whale_tree.bcf /scratch/negishi/allen715/Gray_whales/parentage/ROHs/
# 
bcftools index whale_tree.bcf

# bcftools +fill-tags whale_tree.bcf -- -t AF
bcftools query --allow-extra-chr -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' whale_tree.bcf > allele_frequencies.txt
# 
awk '{OFS="\t"; print $1, $2, $3","$4, $5}' allele_frequencies.txt > allele_frequencies_fixed.txt
# 
gzip allele_frequencies_fixed.txt
tabix -s 1 -b 2 -e 2 allele_frequencies_fixed.txt.gz

bcftools roh --AF-file allele_frequencies_fixed.txt.gz -o bcf_roh_out.txt whale_tree.bcf

###

bcftools roh -G 30 --AF-file allele_frequencies_fixed.txt.gz -o bcf_roh_out.txt whale_tree.bcf

# https://github.com/avril-m-harder/roh_inference_testing/blob/2957f9600f66b4e2bd59acca528a41439f014ee3/empirical/bash/04a_bcftoolsROH.sh#L3
# http://www.htslib.org/doc/bcftools.html#roh
