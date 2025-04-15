#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=whale_tree
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

NGSNGS=/home/allen715/NGSNGS

BLUE=/scratch/negishi/allen715/Gray_whales/final_tree/blue_whale_GCF_009873245.2_mBalMus1.pri.v3_genomic.fna
MINKE=/scratch/negishi/allen715/Gray_whales/final_tree/minke_whale_GCF_949987535.1_mBalAcu1.1_genomic.fna
HUMPBACK=/scratch/negishi/allen715/Gray_whales/final_tree/humpback_GCA_041834305.1_ASM4183430v1_genomic.fna

${NGSNGS}/ngsngs -i ${BLUE} -c 5 -ld Norm,500,50 -seq PE -f fq.gz -o blue_sim -qs 40 -cl 150 -s 3 -t 64 -t2 64
${NGSNGS}/ngsngs -i ${MINKE} -c 5 -ld Norm,500,50 -seq PE -f fq.gz -o minke_sim -qs 40 -cl 150 -s 3 -t 64 -t2 64
${NGSNGS}/ngsngs -i ${HUMPBACK} -c 5 -ld Norm,500,50 -seq PE -f fq.gz -o humpback_sim -qs 40 -cl 150 -s 3 -t 64 -t2 64
