#!/bin/bash
#SBATCH --job-name=whale_dictionary
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 12-00:00:00
#SBATCH --error=dict.err
#SBATCH --output=dict.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load biocontainers
module load bwa
module load picard
module load samtools

#REF=/scratch/bell/allen715/Gray_whales/Reference/original.fa

#bwa index $REF
#samtools faidx $REF
#PicardCommandLine CreateSequenceDictionary reference=$REF output=original.fa.dict

REF=/scratch/negishi/allen715/Gray_whales/reference/ref_100kb.fa

bwa index $REF
samtools faidx $REF
PicardCommandLine CreateSequenceDictionary reference=$REF output=ref.fa.dict
