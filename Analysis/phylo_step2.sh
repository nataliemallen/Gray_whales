#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=iqtree
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 10-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# bcf -> vcf -> fasta -> ml tree with 1000 bootstraps (iqtree2)

PROJ="/scratch/gautschi/allen715/2026_whales"
cd "$PROJ/phylogeny"

# 1) BCF -> compressed, indexed VCF
module --force purge
module load biocontainers
module load bcftools
module load htslib
# bcftools view -O v -o tree_bams.vcf tree_bams.bcf
# bgzip -c tree_bams.vcf > tree_bams.vcf.gz
# tabix -p vcf tree_bams.vcf.gz
# 
# # 2) VCF -> FASTA alignment (IUPAC not resolved), keep sites with >=4 samples
# module load anaconda 2>/dev/null || true
# python3 vcf2phylip.py --input tree_bams.vcf.gz --fasta --nexus \
#   --output-prefix tree_bams.not.iupac.resolved

# 3) ML tree, GTR, 1000 standard bootstraps, outgroup = finwhale
module --force purge
module load biocontainers
module load iqtree/2.2.2.2
iqtree2 -s tree_bams.not.iupac.resolved.min4.fasta -m GTR -B 1000 -T 64 \
  -o "/scratch/gautschi/allen715/2026_whales/reference/outgroup/finwhale.bam" -pre tree_final -st DNA
